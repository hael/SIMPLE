!@descr: kernel PCG reconstruct3D routes: shared-memory solve, distributed worker accumulation and master reduction and solve
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,             only: builder
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_reconstructor_pcg,   only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_LAMBDA, PCG_WEIGHT_THRESHOLD, &
    &pcg_raw_accum_compatible, measure_closed_form_agreement, handle_cold_restart_outcome, report_pcg_solve, &
    &report_closed_form_agreement, write_closed_form_diagnostics, validate_solved_map, read_pcg_raw_accum_header
use simple_matcher_ptcl_io,     only: prepimgbatch, ptcl_read_session, killimgbatch, &
    &ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, rec_state_half, rec_obs_oris, &
    &rec_obs_inv_sig2, read_rec_obs_batch
use simple_matcher_3Drec,       only: write_rec3D_bench
use simple_sigma2_files,        only: load_sigma2
use simple_image,               only: image
use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, evaluate_halfmap_pair, write_halfmap_diagnostics, &
    &write_support_provenance
use simple_image_msk,           only: image_msk
use simple_pcg_solvent_sidecar, only: solvent_prior_cross_half_objective, PCG_SOLVENT_LAMBDA_GRID, &
    &write_pcg_solvent_pair, prepare_solvent_prior_on_pair, solvent_prior_provenance, build_solvent_check_support
use simple_nu_state_filter,     only: nonuniform_filter_state, nu_aux_member, nu_state_filter_timings
use simple_refine3D_fnames,     only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname, &
    &refine3D_pcg_trail_accum_fname, refine3D_pcg_trail_manifest_fname
use simple_trail_chain_manifest, only: trail_chain_manifest, TRAIL_MANIFEST_OK
use simple_oris,                only: population_blend_weights
use simple_state_weight_set,    only: state_weight_set
use simple_pcg_atomic_prior,    only: pcg_atomic_prior, pcg_atomic_prior_stats
use simple_atoms,               only: atoms
!$ use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_get_max_active_levels, &
!$     &omp_set_max_active_levels, omp_set_num_threads
implicit none

public :: execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
public :: rec3D_master_nthr, validate_pcg_common
public :: set_pcg_solve_support, calculate_pcg_state_diagnostics, accumulate_pcg_state_half
public :: pcg_assemble_bench, write_pcg_assemble_bench
private
#include "simple_local_flags.inc"

integer, parameter :: PCG_MASTER_NTHR_CAP    = 32   !< master-phase thread-boost ceiling
! Solve support follows automsk (reconstruct3D_pcg_policy.md, section 3):
! the sphere for automsk=no, the lag-one density envelope for yes and nu
! (never the NU-evidence envelope), an explicit pcg_mskfile over both.
logical, parameter :: DEBUG = .false.
! atomicity prior: fitted B/occupancy/posterior intensity with a species table; otherwise free amplitudes
integer, parameter :: ATM_AMP_MAXITS = 30
real,    parameter :: ATM_AMP_RTOL   = 1.0e-3

!> wall-clock phases of one distributed PCG master assembly (BENCH/PCG_ASSEMBLE_BENCH.txt); index
!! 1 of the pairs is the base solve kind, 2 the regularized (ml); the caller times the handoff and total
type :: pcg_assemble_bench
    real(dp) :: rt_init = 0._dp, rt_support = 0._dp, rt_op_setup = 0._dp, rt_chain = 0._dp
    real(dp) :: rt_read = 0._dp, rt_sum = 0._dp, rt_finalize(2) = 0._dp, rt_solve(2) = 0._dp
    real(dp) :: rt_prior_resolve = 0._dp, rt_finish = 0._dp
    real(dp) :: rt_fsc = 0._dp, rt_map_output = 0._dp, rt_nu_filter = 0._dp
    real(dp) :: rt_proj_write = 0._dp, rt_cleanup = 0._dp, rt_nu_handoff = 0._dp, rt_tot = 0._dp
    real(dp) :: rt_half_solve(2,2) = 0._dp !< (even/odd, kind): each half's own solve time
    type(nu_state_filter_timings) :: nu
    integer(int64) :: raw_bytes = 0_int64
    integer  :: nraw = 0, niters(2,2) = 0, nstates_solved = 0, nthr_master = 0, nthr_half = 0
end type pcg_assemble_bench

contains

    ! Solve contract (reconstruct3D_pcg_policy.md, 'ML two-map contract and
    ! starts'): the base solve starts from zero every iteration
    ! (reconstructor_pcg%solve_with_cold_restart); the regularized half is the
    ! closed-form shrink of the base solution followed by maxits_ml coupled
    ! iterations from it (reconstructor_pcg%solve_regularized).

    !> Half-map diagnostics of a state through the backend-neutral evaluator,
    !! with the PCG mask policy (spherical msk_crop on the base pair) and the
    !! automask write on the envfsc path; both routes call it.
    subroutine calculate_pcg_state_diagnostics( params, state_here, context, even, odd, avg, diagnostics, &
        &l_pair_support_constrained, support_kind )
        class(parameters),                intent(in)  :: params
        integer,                          intent(in)  :: state_here
        character(len=*),                 intent(in)  :: context
        class(image),                     intent(in)  :: even, odd, avg
        type(halfmap_diagnostics_result), intent(out) :: diagnostics
        logical,                          intent(in)  :: l_pair_support_constrained
        character(len=*),                 intent(in)  :: support_kind
        type(image) :: envmask
        if( params%l_envfsc .and. trim(params%automsk) /= 'nu' )then
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, diagnostics, 'pcg', envmask=envmask, &
                &l_pair_support_constrained=l_pair_support_constrained, support_kind=support_kind, mask_kind='density')
            call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
            call envmask%kill
        else
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, diagnostics, 'pcg', &
                &l_pair_support_constrained=l_pair_support_constrained, support_kind=support_kind, mask_kind='none')
        endif
        ! the resolution document names the solvent prior beside the FSC mode
        if( params%l_pcg_solvent ) diagnostics%fsc_mode = &
            &trim(diagnostics%fsc_mode)//' solvent_prior=soft(per_half,base_pair)'
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.500 RESOLUTION = ', diagnostics%res_fsc05
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.143 RESOLUTION = ', diagnostics%res_fsc0143
    end subroutine calculate_pcg_state_diagnostics

    !> The resolution document's name, as gridding's volassemble names it: an
    !! explicit outfile, else the state name tagged with which_iter.
    function resolve_pcg_fsc_txt_fname( params, cline, state ) result( fname )
        class(parameters), intent(in) :: params
        class(cmdline),    intent(in) :: cline
        integer,           intent(in) :: state
        type(string) :: fname, ext
        if( cline%defined('outfile') )then
            fname = params%outfile
            ext   = fname2ext(fname)
            select case(ext%to_char())
                case('txt','simple')
                    fname = get_fbody(fname, ext)
            end select
            fname = fname//'_STATE'//int2str_pad(state,2)
            call ext%kill
        else if( cline%defined('which_iter') )then
            fname = refine3D_resolution_txt_fbody(state, params%which_iter)
        else
            fname = refine3D_resolution_txt_fbody(state)
        endif
    end function resolve_pcg_fsc_txt_fname

    !> Install the solve support P of (P H P) u = P b: an explicit pcg_mskfile,
    !! else the state's density envelope when built, else the mskdiam sphere.
    subroutine set_pcg_solve_support( pcgop, params, state_support, l_state_support )
        class(reconstructor_pcg),  intent(inout) :: pcgop
        class(parameters),         intent(in)    :: params
        class(image),    optional, intent(in)    :: state_support
        logical,         optional, intent(in)    :: l_state_support
        type(image) :: mskvol
        logical     :: l_state
        if( params%pcg_mskfile%is_allocated() )then
            if( len_trim(params%pcg_mskfile%to_char()) > 0 )then
                call mskvol%read_and_crop(params%pcg_mskfile, params%smpd, params%box_crop, params%smpd_crop)
                call pcgop%set_window_volume(mskvol)
                call mskvol%kill
                return
            endif
        endif
        l_state = .false.
        if( present(l_state_support) ) l_state = l_state_support
        if( l_state .and. present(state_support) )then
            call pcgop%set_window_volume(state_support)
            return
        endif
        call pcgop%set_window_sphere(params%msk_crop)
    end subroutine set_pcg_solve_support

    !> The state's solve support: an explicit pcg_mskfile, else under automsk
    !! the density envelope of the lag-one reference; without a reference the
    !! base bootstraps on the sphere (reconstruct3D_pcg_policy.md, section 3).
    subroutine build_pcg_state_support( params, state_here, support, l_have, support_kind )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        type(image_msk),   intent(inout) :: support
        logical,           intent(out)   :: l_have
        character(len=*),  intent(out)   :: support_kind
        type(image)  :: vol_prev
        l_have = .false.
        support_kind = 'sphere'
        call support%kill_bimg
        ! an explicit, non-empty pcg_mskfile constrains every solve whatever
        ! automsk, and is reported as the state support
        if( params%pcg_mskfile%is_allocated() )then
            if( len_trim(params%pcg_mskfile%to_char()) > 0 )then
                call support%read_and_crop(params%pcg_mskfile, params%smpd, params%box_crop, params%smpd_crop)
                write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                    &', explicit pcg_mskfile '//params%pcg_mskfile%to_char()//' constrains base and replay'
                l_have = .true.
                support_kind = 'explicit'
                return
            endif
        endif
        ! without automsk every solve runs on the sphere
        if( .not. pcg_density_support_enabled(params) )then
            write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                &', spherical support for base and replay (automsk=no; density envelope disabled)'
            return
        endif
        if( .not. params%l_envfsc ) &
            &THROW_HARD('active automsk requires envfsc=yes (derived in parameters); the coupling was bypassed')
        ! the NU-evidence envelope is never the solve support, under automsk=nu
        ! either (reconstruct3D_pcg_policy.md, section 3)
        if( state_here < 1 .or. state_here > size(params%vols) ) then
            call handle_missing_reference('no reference volume slot')
            return
        endif
        if( len_trim(params%vols(state_here)%to_char()) == 0 )then
            call handle_missing_reference('no reference volume recorded')
            return
        endif
        if( .not. file_exists(params%vols(state_here)) )then
            call handle_missing_reference('missing reference volume')
            return
        endif
        call vol_prev%read_and_crop(params%vols(state_here), params%smpd, params%box_crop, params%smpd_crop)
        call build_pcg_density_support(params, state_here, vol_prev, support, 'lag-one reference')
        call vol_prev%kill
        l_have = .true.
        support_kind = 'density'

    contains

        subroutine handle_missing_reference( why )
            character(len=*), intent(in) :: why
            write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                &' has '//trim(why)//'; bootstrap base uses the sphere and replay support derives from its density fallback'
        end subroutine handle_missing_reference

    end subroutine build_pcg_state_support

    !> Automatic solve support is enabled for yes and nu. An explicit
    !! pcg_mskfile is the development override and is installed regardless.
    logical function pcg_density_support_enabled( params ) result( l_enabled )
        class(parameters), intent(in) :: params
        l_enabled = trim(params%automsk) .ne. 'no'
    end function pcg_density_support_enabled

    !> The conservative density support from a volume: the lag-one reference,
    !! or the completed base pair when no reference exists yet.
    subroutine build_pcg_density_support( params, state_here, volume, support, source )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        class(image),      intent(in)    :: volume
        type(image_msk),   intent(inout) :: support
        character(len=*),  intent(in)    :: source
        call support%automask3D(params, volume, .false., lp_override=params%envmsklp, l_report=.false.)
        write(logfhandle,'(A,I0,A,F6.1,A,A,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
            &', conservative density envelope at ', params%envmsklp, ' A from ', trim(source), &
            &' (replaces the spherical support)'
    end subroutine build_pcg_density_support

    !> Master-phase thread budget of both backends: on local execution the
    !! workers are idle, so the full allocation up to PCG_MASTER_NTHR_CAP;
    !! nthr_master on a cluster.
    integer function rec3D_master_nthr( params, nthr_master ) result( nthr )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nthr_master
        nthr = nthr_master
        if( trim(params%qsys_name) == 'local' )then
            nthr = max(params%nthr, min(PCG_MASTER_NTHR_CAP, max(1, params%nparts) * params%nthr))
            !$ nthr = min(omp_get_num_procs(), nthr)
        endif
    end function rec3D_master_nthr

    !> Distributed worker: accumulate and atomically publish raw full-range B
    !! and real D for every local (state,half). Workers never call end_accum;
    !! folding and every nonlinear finalization step belong to the master.
    !! With a state weight set, a particle is a member of every state it weighs
    !! into above PCG_WEIGHT_THRESHOLD and enters with its noise spectrum over w.
    subroutine execute_rec3D_pcg_worker( params, build, cline, selected_pinds, wset, wtab )
        type(parameters),                  intent(inout) :: params
        type(builder),                     intent(inout) :: build
        class(cmdline),                    intent(inout) :: cline
        integer,                           intent(in)    :: selected_pinds(:)
        class(state_weight_set), optional, intent(inout) :: wset
        real,                    optional, intent(in)    :: wtab(:,:)   !< caller's weight table, rows of selected_pinds
        integer, allocatable :: half_pinds(:)
        real,    allocatable :: half_w(:)
        logical :: l_weighted
        character(len=256) :: provenance
        integer(int64)     :: wset_id(2)
        integer :: state, eo, n_half, nimgs
        logical :: l_sigma_loaded
        type(ptcl_read_counts)  :: reads0, reads
        type(ptcl_read_session) :: rsess
        integer(timer_int_kind) :: t0, t_half
        real(timer_int_kind)    :: rt_init, rt_half, rt_prep

        if( L_BENCH_GLOB )then
            t0      = tic()
            reads0  = ptcl_read_counts_now()
            rt_half = 0.
            rt_prep = 0.
        endif
        nimgs = 0
        call validate_pcg_common(params, check_solver=.false.)
        provenance = pcg_raw_provenance(params)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            l_sigma_loaded = allocated(build%esig%sigma2_noise)
            if( .not. l_sigma_loaded )then
                call load_sigma2(params, build%pftc, build%esig, build%spproj, &
                    &build%spproj_field, l_sigma_loaded)
            endif
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif
        if( size(selected_pinds) > 0 )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        call rsess%new(params)
        wset_id = 0_int64
        if( present(wset) )then
            if( wset%get_nstates() /= params%nstates ) THROW_HARD('nstates differs from the state weight set; PCG worker')
            wset_id = [wset%get_generation(), wset%get_layout_digest()]
        endif
        if( present(wtab) )then
            if( size(wtab,1) /= size(selected_pinds) .or. size(wtab,2) /= params%nstates ) THROW_HARD('weight table shape; PCG worker')
        endif
        l_weighted = present(wset) .or. present(wtab)
        if( L_BENCH_GLOB ) rt_init = toc(t0)
        do state = 1, params%nstates
            do eo = 0, 1
                if( l_weighted )then
                    call collect_weighted_state_half(state, eo, selected_pinds, half_pinds, half_w)
                else
                    call rec_state_half(build, state, eo, selected_pinds, half_pinds)
                endif
                n_half = size(half_pinds)
                nimgs  = nimgs + n_half
                if( L_BENCH_GLOB ) t_half = tic()
                call accumulate_worker_state_half(state, eo, half_pinds, provenance)
                if( L_BENCH_GLOB ) rt_half = rt_half + toc(t_half)
                deallocate(half_pinds)
                if( allocated(half_w) ) deallocate(half_w)
                if( DEBUG )then
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> PCG RAW WORKER: PART ', params%part, &
                        &' STATE ', state, ' HALF ', eo, ' PARTICLES ', n_half
                endif
            enddo
        enddo
        call rsess%kill
        if( size(selected_pinds) > 0 ) call killimgbatch(build)
        if( L_BENCH_GLOB )then
            reads = ptcl_read_counts_since(reads0)
            call write_rec3D_bench(params, 'pcg', size(selected_pinds), nimgs, &
                &reads, rt_init, rt_prep, rt_half - reads%rt_read - rt_prep, toc(t0))
        endif

    contains

        !> the state's weighted members of one half (only params%state when state= was given)
        subroutine collect_weighted_state_half( state_here, eo_here, pinds, selected, weights )
            integer,              intent(in)  :: state_here, eo_here, pinds(:)
            integer, allocatable, intent(out) :: selected(:)
            real,    allocatable, intent(out) :: weights(:)
            integer, allocatable :: members(:)
            real,    allocatable :: wmem(:)
            logical, allocatable :: l_half(:)
            integer :: i
            if( params%l_state_defined .and. state_here /= params%state )then
                allocate(selected(0), weights(0))
                return
            endif
            if( present(wset) )then
                call wset%get_members(state_here, PCG_WEIGHT_THRESHOLD, pinds, members, wmem)
            else
                members = pack(pinds,              wtab(:,state_here) > PCG_WEIGHT_THRESHOLD)
                wmem    = pack(wtab(:,state_here), wtab(:,state_here) > PCG_WEIGHT_THRESHOLD)
            endif
            allocate(l_half(size(members)))
            do i = 1, size(members)
                l_half(i) = build%spproj_field%get_eo(members(i)) == eo_here
            enddo
            selected = pack(members, l_half)
            weights  = pack(wmem,    l_half)
            deallocate(members, wmem, l_half)
        end subroutine collect_weighted_state_half

        subroutine accumulate_worker_state_half( state_here, eo_here, pinds, provenance_here )
            integer,          intent(in) :: state_here, eo_here, pinds(:)
            character(len=*), intent(in) :: provenance_here
            type(reconstructor_pcg) :: pcgop
            type(string)    :: fname

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            fname = refine3D_pcg_raw_accum_fname(state_here, params%part, params%numlen, &
                &merge('odd ', 'even', eo_here == 1))
            if( size(pinds) == 0 )then
                call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                    &params%nparts, 0, provenance_here, mass=0._dp, wset_id=wset_id)
                call pcgop%kill
                call fname%kill
                return
            endif
            if( l_weighted )then
                call accumulate_pcg_state_half(params, build, rsess, eo_here, pinds, pcgop, rt_prep, half_w)
            else
                call accumulate_pcg_state_half(params, build, rsess, eo_here, pinds, pcgop, rt_prep)
            endif
            if( l_weighted )then
                call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                    &params%nparts, size(pinds), provenance_here, mass=sum(real(half_w,dp)), wset_id=wset_id)
            else
                call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                    &params%nparts, size(pinds), provenance_here)
            endif
            call pcgop%kill
            call fname%kill
        end subroutine accumulate_worker_state_half

    end subroutine execute_rec3D_pcg_worker

    !> The accumulation of the particles pinds, all in half eo_here, into pcgop (made by the caller): the
    !! point group, the half's noise spectrum at the native shells (objfun=euclid), the particles' orientations
    !! and CTF at the cropped box, and one pass over their observation batches through rsess. The frame the
    !! PCG worker publishes and reconstruct3D_wave's start finalizes. weights: the state weights, an inverse
    !! noise scale that multiplies the B term and |T|^2 alike
    subroutine accumulate_pcg_state_half( params, build, rsess, eo_here, pinds, pcgop, rt_prep, weights )
        class(parameters),        intent(in)    :: params
        class(builder),           intent(inout) :: build
        class(ptcl_read_session), intent(inout) :: rsess
        integer,                  intent(in)    :: eo_here, pinds(:)
        class(reconstructor_pcg), intent(inout) :: pcgop
        real(timer_int_kind),     intent(inout) :: rt_prep
        real, optional,           intent(in)    :: weights(:)
        type(oris)  :: selection
        type(image) :: obs
        complex, allocatable :: y_batch(:,:,:)
        real,    allocatable :: inv_sig2(:)
        integer :: lims2(2,2), R, batchlims(2), batchsz, ibatch
        call pcgop%set_sym(build%pgrpsyms)
        lims2 = pcgop%get_lims2()
        R     = lims2(1,2)
        ! every particle of the call is in half eo_here: one inverse spectrum at the native shells 0:R
        if( params%cc_objfun == OBJFUN_EUCLID ) call rec_obs_inv_sig2(params, build, eo_here, R, inv_sig2)
        call rec_obs_oris(params, build, pinds, selection)
        if( allocated(inv_sig2) )then
            if( present(weights) )then
                call pcgop%prep_particles(selection, use_ctf=.true., inv_sig2=inv_sig2, weights=weights)
            else
                call pcgop%prep_particles(selection, use_ctf=.true., inv_sig2=inv_sig2)
            endif
        else
            if( present(weights) )then
                call pcgop%prep_particles(selection, use_ctf=.true., weights=weights)
            else
                call pcgop%prep_particles(selection, use_ctf=.true.)
            endif
        endif
        allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
        call pcgop%begin_accum
        call obs%new([params%box_crop,params%box_crop,1], params%smpd_crop)
        do ibatch = 1, size(pinds), MAXIMGBATCHSZ
            batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            ! the search's records of a cropped pass are the batch's observations, else the batch is read
            call read_rec_obs_batch(params, build, rsess, pinds, batchlims, obs, lims2, R*R, y_batch, rt_prep)
            call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
        enddo
        call obs%kill
        call selection%kill
        deallocate(y_batch)
        if( allocated(inv_sig2) ) deallocate(inv_sig2)
    end subroutine accumulate_pcg_state_half

    !> Distributed master: reduce raw worker B,D artifacts in ascending part
    !! order, then perform all folding, finalization and PCG locally. For each
    !! even/odd pair, construction and teardown stay serial while the two fully
    !! prepared PCG solves execute concurrently with disjoint thread budgets.
    subroutine execute_rec3D_pcg_distributed_master( params, build, cline, nu_align_lps, bench )
        type :: distributed_half_job
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            real, allocatable :: x(:,:,:), x_cf(:,:,:), rel_res_hist(:)
            ! the base operator kept for the regularized pair, with its raw density and the count and mass it represents
            real, allocatable :: d_raw(:,:,:)
            integer  :: nptcls_raw = 0
            real(dp) :: mass_raw = 0.0_dp
            logical  :: l_keep_op = .false.
            integer :: state = 0, eo = 0, nptcls = 0, niters = 0
            real(dp) :: mass = 0.0_dp  !< applied weight mass of the reduction (nptcls under hard labels)
            integer :: band_shell = 0 !< the pair's FSC=0.143 shell (regularized solve; agreement diagnostic)
            integer :: prior_npositive = 0
            character(len=8) :: half = '', solve_kind = ''
            real(dp) :: time_reduce = 0.0_dp, time_finalize = 0.0_dp, time_solve = 0.0_dp
            real :: prior_positive_min = 0.0, prior_positive_max = 0.0
            real :: prior_to_khat_l1 = 0.0, prior_to_khat_rms = 0.0
            logical :: l_ml_solve = .false., ready = .false., l_concurrent = .false.
            logical :: l_nonzero = .false. !< nonzero initial guess: eligible for the cold restart
            logical :: l_atm_solve = .false. !< the atomicity-prior re-solve: maxits_atm/rtol_atm budget
        end type distributed_half_job
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        real,    optional, intent(out)  :: nu_align_lps(:)
        type(pcg_assemble_bench), optional, intent(inout) :: bench
        type(image)  :: half_even, half_odd, ml_even, ml_odd, merged, solvent_even, solvent_odd, atm_even, atm_odd
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt, eonames(2)
        type(halfmap_diagnostics_result) :: hm_diag
        real, allocatable :: fsc(:), res0143s(:), res05s(:), cfars(:), align_lps(:)
        real, allocatable :: realized_fractions(:), update_weights(:), chain_weights(:), current_scales(:)
        integer, allocatable :: nrep(:), nsmp(:), rep_rows(:), smp_rows(:), rec_rows(:)
        real,    allocatable :: rnrep(:), rnsmp(:) !< N(s), n(s) of the population rule: applied masses
        logical, allocatable :: state_written(:), state_carried(:), state_empty(:)
        type(state_weight_set) :: wset
        integer(int64) :: wset_id(2)
        real(dp)       :: chain_mass(2)
        integer        :: chain_contrib(2), chain_gen_prev
        logical        :: chain_written(2), l_weighted
        character(len=:), allocatable :: errmsg
        character(len=256) :: provenance, chain_provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        real    :: n_active_state, n_sampled_state
        integer :: pcg_master_nthreads, pcg_half_nthreads
        type(distributed_half_job) :: even_job, odd_job
        type(image_msk) :: state_support_msk
        type(image)     :: solvent_weight(2)
        logical :: l_state_support, l_base_support_constrained, l_nu_base_constrained
        logical :: l_has_updates, l_sample_seed, l_chain_carry, l_even_chain, l_odd_chain
        logical :: l_shipped_support_constrained, l_solvent_weight
        real    :: res0143_prior_free
        real    :: solvent_lambda_eff !< the solvent ridge in force for the current state
        ! atm_prior=yes: the prior shared by both half operators, its strength and the products of the current state
        type(pcg_atomic_prior), target :: atm_prior
        type(image)  :: atm_res
        real(dp), allocatable :: atm_c(:,:)       !< (natoms,2) amplitudes fitted per half
        real    :: atm_lambda_rel, atm_lambda_abs(2), atm_amp_rres(2), atm_resid(2), atm_mres(2)
        integer :: atm_amp_nits(2), atm_nits(2)
        logical :: l_atm_resolved
        character(len=16) :: state_support_kind, base_support_kind, shipped_support_kind
        integer(timer_int_kind) :: t_state_phase, t_bench
        real(dp) :: time_map_output, time_fsc_output, time_nu_filter
        type(pcg_assemble_bench) :: bt !< this assembly's phases, handed to bench

        if( L_BENCH_GLOB ) t_bench = tic()
        call validate_pcg_common(params)
        solvent_lambda_eff = 0.
        ! add-on chains reconstruct by gridding only: the frozen-set handshakes have no PCG side
        if( cline%defined('frozen_rec') .or. cline%defined('frozen_seed') )then
            THROW_HARD('frozen sets are gridding only; rec_backend=pcg takes no frozen_rec or frozen_seed')
        endif
        ! m_estimator=flex: masses come from the state weight set (membership above PCG_WEIGHT_THRESHOLD)
        l_weighted = params%l_m_estimator_flex
        wset_id    = 0_int64
        if( l_weighted )then
            call wset%new(build%spproj, build%spproj_field)
            if( wset%get_nstates() /= params%nstates ) THROW_HARD('nstates differs from the state weight set; PCG master')
            wset_id = [wset%get_generation(), wset%get_layout_digest()]
        endif
        ! master phase: the full local allocation while the workers are idle
        ! (rec3D_master_nthr), restored to nthr before returning
        !$ if( trim(params%qsys_name) == 'local' ) &
        !$ &call omp_set_num_threads(rec3D_master_nthr(params, params%nthr))
        pcg_master_nthreads = 1
        !$ pcg_master_nthreads = omp_get_max_threads()
        pcg_half_nthreads = max(1, pcg_master_nthreads / 2)
        if( pcg_master_nthreads >= 2 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> PCG DISTRIBUTED: EVEN/ODD SOLVES RUN CONCURRENTLY (', &
                &pcg_half_nthreads, ' THREADS PER HALF; ', pcg_master_nthreads, ' MASTER THREADS AVAILABLE)'
        endif
        if( present(nu_align_lps) )then
            if( size(nu_align_lps) /= params%nstates ) &
                &THROW_HARD('PCG NU matching low-pass output has invalid size')
            nu_align_lps = 0.0
        endif
        allocate(align_lps(params%nstates), source=0.0)
        provenance = pcg_raw_provenance(params)
        chain_provenance = pcg_chain_provenance(params)
        l_has_updates = .false.
        do iptcl = params%fromp, params%top
            if( build%spproj_field%get_updatecnt(iptcl) > 0 )then
                l_has_updates = .true.
                exit
            endif
        enddo
        allocate(realized_fractions(params%nstates), source=1.0)
        allocate(update_weights(params%nstates), source=1.0)
        allocate(chain_weights(params%nstates),  source=0.0)
        allocate(current_scales(params%nstates), source=1.0)
        allocate(rnrep(params%nstates), rnsmp(params%nstates), source=0.)
        if( params%l_trail_rec )then
            ! N(s) and n(s) of the population rule; f = n/N
            if( l_weighted )then
                call build%spproj%os_ptcl3D%get_update_rows(rep_rows, smp_rows, rec_rows)
                do istate = 1, params%nstates
                    rnrep(istate) = real(wset%applied_mass(istate, PCG_WEIGHT_THRESHOLD, rep_rows))
                    rnsmp(istate) = real(wset%applied_mass(istate, PCG_WEIGHT_THRESHOLD, smp_rows))
                    realized_fractions(istate) = 0.
                    if( rnrep(istate) > 0. ) realized_fractions(istate) = rnsmp(istate) / rnrep(istate)
                enddo
            else
                call build%spproj%os_ptcl3D%get_group_update_counts('state', params%nstates, nrep, nsmp)
                call build%spproj%os_ptcl3D%get_state_update_fracs(params%nstates, realized_fractions)
                rnrep = real(nrep)
                rnsmp = real(nsmp)
            endif
            update_weights = realized_fractions
            if( params%l_ufrac_trec_defined )then
                if( params%nstates == 1 )then
                    update_weights(1) = params%ufrac_trec
                else
                    THROW_WARN('ufrac_trec ignored for multi-state PCG; using realized state update fractions')
                endif
            endif
        endif
        allocate(res0143s(params%nstates), source=0.0)
        allocate(res05s(params%nstates),   source=0.0)
        allocate(cfars(params%nstates),    source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        allocate(state_carried(params%nstates), source=.false.)
        call find_empty_states(state_empty)
        if( L_BENCH_GLOB )then
            bt%rt_init     = toc(t_bench)
            bt%nthr_master = pcg_master_nthreads
            bt%nthr_half   = pcg_half_nthreads
        endif
        do state = 1, params%nstates
            ! state= reconstructs that state alone: the other states' parts are empty by construction
            if( params%l_state_defined .and. state /= params%state ) cycle
            ! the gridding rule (volassemble's dropped states): a state without particles, while
            ! another state has some, is dropped for this iteration and its previous maps are
            ! carried forward; a trailing chain it still has is carried by the chain path below
            if( state_empty(state) )then
                l_even_chain = .false.
                if( params%l_trail_rec )then
                    call discard_stale_trail_chain_pair(state)
                    raw_fname    = refine3D_pcg_trail_accum_fname(state, 'even')
                    l_even_chain = file_exists(raw_fname)
                endif
                if( .not. l_even_chain )then
                    if( .not. previous_state_maps_exist(state) )then
                        errmsg = 'state '//int2str(state)//' has no particles and no previous maps in this '//&
                            &'directory to carry forward; distributed PCG'
                        THROW_HARD(errmsg)
                    endif
                    write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                        &' HAS NO PARTICLES; DROPPING AND CARRYING PREVIOUS VOLUME FORWARD'
                    state_carried(state) = .true.
                    cycle
                endif
            endif
            ! l_sample_seed: no chain yet, so the chain is seeded from the current
            ! sample at full mass and this iteration ships the current sample's map;
            ! l_chain_carry: a sample below 0.001 of the state, so the chain carries
            ! unchanged and this iteration is solved from it (gridding volassemble)
            l_sample_seed = .false.
            l_chain_carry = .false.
            chain_gen_prev = 0
            if( params%l_trail_rec )then
                ! a chain pair of another identity is discarded and re-seeded
                call discard_stale_trail_chain_pair(state)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'even')
                l_even_chain = file_exists(raw_fname)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'odd')
                l_odd_chain = file_exists(raw_fname)
                if( l_even_chain .neqv. l_odd_chain ) THROW_HARD('PCG trailing chain pair is incomplete')
                l_sample_seed = .not. l_even_chain
                n_active_state  = rnrep(state)
                n_sampled_state = rnsmp(state)
                if( n_sampled_state > 0 .and. realized_fractions(state) < 0.001 )then
                    if( .not. l_sample_seed )then
                        l_chain_carry = .true.
                        write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                            &' HAS NO SAMPLE; TRAILING CHAIN CARRIED UNCHANGED'
                    else if( previous_state_maps_exist(state) )then
                        ! nothing to restore: the previous maps stay as they are on disk
                        write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                            &' HAS NO SAMPLE AND NO TRAILING CHAIN; CARRYING PREVIOUS VOLUME FORWARD'
                        state_carried(state) = .true.
                        cycle
                    else
                        ! as gridding: a map of the sample alone would be degenerate
                        errmsg = 'state '//int2str(state)//': realized update fraction below 0.001, no trailing '// &
                            &'chain and no previous volume in this directory; nothing to restore (raise update_frac)'
                        THROW_HARD(errmsg)
                    endif
                else if( .not. l_sample_seed )then
                    call set_chain_blend_weights(state)
                endif
                if( l_sample_seed .and. n_sampled_state > 0 )then
                    write(logfhandle,'(A,I0,A,F12.2,A)') &
                        &'>>> PCG DISTRIBUTED: SEEDED FULL-MASS TRAILING CHAIN, STATE ', state, &
                        &', REPRESENTED POPULATION ', rnrep(state), '; THIS ITERATION USES THE CURRENT SAMPLE ONLY'
                endif
            endif
            ! one solve support per state (build_pcg_state_support)
            if( L_BENCH_GLOB ) t_bench = tic()
            call build_pcg_state_support(params, state, state_support_msk, l_state_support, state_support_kind)
            if( L_BENCH_GLOB ) bt%rt_support = bt%rt_support + toc(t_bench)
            l_solvent_weight = .false.
            l_base_support_constrained = l_state_support
            base_support_kind = 'sphere'
            if( l_base_support_constrained ) base_support_kind = state_support_kind
            chain_written = .false.
            chain_mass    = 0.0_dp
            chain_contrib = 0
            l_atm_resolved = .false.
            call reduce_solve_state_pair(state, half_even, half_odd, n_even, n_odd, 'base', &
                &solvent_even=solvent_even, solvent_odd=solvent_odd, atm_even=atm_even, atm_odd=atm_odd)
            ! the chain pair written by the base reductions is published by its manifest, last
            if( all(chain_written) ) call write_chain_manifest(state)
            if( params%l_trail_rec .and. .not. l_chain_carry )then
                if( l_weighted )then
                    if( abs(real(even_job%mass+odd_job%mass) - n_sampled_state) > 1.0e-3*max(1.,n_sampled_state) ) &
                        &THROW_HARD('PCG raw applied mass does not match the latest sampled cohort')
                else
                    if( real(n_even+n_odd) /= n_sampled_state ) THROW_HARD('PCG raw particles do not match the latest sampled cohort')
                    if( n_active_state > 0. )then
                        if( abs(n_sampled_state/n_active_state-realized_fractions(state)) > 1.0e-6 )then
                            THROW_HARD('PCG realized fraction disagrees with gridding sampling bookkeeping')
                        endif
                    endif
                endif
            endif
            if( n_even == 0 .and. n_odd == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                    &' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even < 1 .or. n_odd < 1 )then
                ! a state left with members in one halfset only (the halves of class averages are
                ! assigned independently): its previous maps are carried forward, as for an empty
                ! state; gridding assembles a zero half instead and continues as well
                if( params%l_trail_rec .or. .not. previous_state_maps_exist(state) )&
                    &THROW_HARD('distributed PCG requires both halfsets')
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                    &' HAS PARTICLES IN ONE HALFSET ONLY; CARRYING PREVIOUS VOLUME FORWARD'
                state_carried(state) = .true.
                cycle
            endif
            fname_even = refine3D_state_halfvol_fname(state, 'even')
            fname_odd  = refine3D_state_halfvol_fname(state, 'odd')
            fname_vol  = refine3D_state_vol_fname(state)
            fname_fsc  = refine3D_fsc_fname(state)
            call merged%copy(half_even)
            call merged%add(half_odd)
            call merged%mul(0.5)
            if( params%l_ml_reg .and. .not. l_state_support .and. pcg_density_support_enabled(params) )then
                if( L_BENCH_GLOB ) t_bench = tic()
                call build_pcg_density_support(params, state, merged, state_support_msk, 'current base pair')
                if( L_BENCH_GLOB ) bt%rt_support = bt%rt_support + toc(t_bench)
                l_state_support = .true.
                state_support_kind = 'density'
            endif
            time_map_output = 0.0_dp
            time_nu_filter  = 0.0_dp
            if( params%l_ml_reg )then
                fname_even_unfil = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
                fname_odd_unfil  = refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.)
                t_state_phase = tic()
                call half_even%write(fname_even_unfil, del_if_exists=.true.)
                call half_odd%write(fname_odd_unfil, del_if_exists=.true.)
                time_map_output = real(toc(t_state_phase),dp)
            endif
            if( l_solvent_weight )then
                if( L_BENCH_GLOB ) t_bench = tic()
                call write_pcg_solvent_pair(params, state, solvent_even, solvent_odd)
                if( L_BENCH_GLOB ) bt%rt_map_output = bt%rt_map_output + toc(t_bench)
            endif
            t_state_phase = tic()
            ! the FSC pair is the current base pair, which also seeds the NU bank
            call calculate_pcg_state_diagnostics(params, state, 'DISTRIBUTED', half_even, &
                &half_odd, merged, hm_diag, l_base_support_constrained, base_support_kind)
            fsc             = hm_diag%fsc
            res0143s(state) = hm_diag%res_fsc0143
            res05s(state)   = hm_diag%res_fsc05
            cfars(state)    = hm_diag%cfar
            call arr2file(fsc, fname_fsc)
            fname_restxt = resolve_pcg_fsc_txt_fname(params, cline, state)
            call write_halfmap_diagnostics(hm_diag, params%box_crop, params%smpd_crop, fname_restxt)
            call fname_restxt%kill
            call hm_diag%kill
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                ! the ordinary global-ML replay (P_tau from the FSC pair), of
                ! the solvent-prior'd pair when the prior is on
                if( l_solvent_weight )then
                    call reduce_solve_state_pair(state, ml_even, ml_odd, n_even, n_odd, 'ml', fsc, &
                        &solvent_even, solvent_odd)
                else
                    call reduce_solve_state_pair(state, ml_even, ml_odd, n_even, n_odd, 'ml', fsc, &
                        &half_even, half_odd)
                endif
                call merged%kill
                call merged%copy(ml_even)
                call merged%add(ml_odd)
                call merged%mul(0.5)
            endif
            t_state_phase = tic()
            l_shipped_support_constrained = merge(l_state_support, l_base_support_constrained, params%l_ml_reg)
            shipped_support_kind = base_support_kind
            if( params%l_ml_reg ) shipped_support_kind = state_support_kind
            if( params%l_ml_reg )then
                call ml_even%write(fname_even, del_if_exists=.true.)
                call ml_odd%write(fname_odd, del_if_exists=.true.)
            else if( l_atm_resolved )then
                ! the prior-free pair keeps its products under _base, the prior'd pair ships
                call write_atomic_prior_products(state)
                call atm_even%write(fname_even, del_if_exists=.true.)
                call atm_odd%write(fname_odd, del_if_exists=.true.)
                call merged%kill
                call merged%copy(atm_even)
                call merged%add(atm_odd)
                call merged%mul(0.5)
            else
                call half_even%write(fname_even, del_if_exists=.true.)
                call half_odd%write(fname_odd, del_if_exists=.true.)
            endif
            call merged%write(fname_vol, del_if_exists=.true.)
            ! the sidecar follows the published map, never precedes it
            if( params%l_ml_reg .and. l_solvent_weight )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'regularized', shipped_support_kind, &
                    &solvent_prior=solvent_prior_provenance(params, solvent_lambda_eff))
            else if( params%l_ml_reg )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'regularized', shipped_support_kind)
            else if( l_solvent_weight )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'base', shipped_support_kind, &
                    &solvent_prior=solvent_prior_provenance(params, solvent_lambda_eff))
            else if( l_atm_resolved )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'base', shipped_support_kind, &
                    &atomic_prior=atomic_prior_provenance())
            else
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'base', shipped_support_kind)
            endif
            time_map_output = time_map_output + real(toc(t_state_phase),dp)
            if( params%l_nonuniform )then
                ! NU competition as gridding's volassemble: the bank from the current
                ! base pair, the ML pair as the auxiliary member
                t_state_phase = tic()
                eonames(1) = fname_even
                eonames(2) = fname_odd
                ! a constrained base pair hands its support over (the evidence null
                ! on the dilation ring)
                l_nu_base_constrained = l_base_support_constrained
                if( l_solvent_weight )then
                    ! label field from the base pair, applied to the prior'd pair
                    call nonuniform_filter_state(params, state, half_even, half_odd, &
                        &ml_even, ml_odd, nu_aux_member(params), &
                        &res0143s(state), fname_vol, eonames, align_lps(state), timings=bt%nu, &
                        &base_support=state_support_msk, l_base_constrained=l_nu_base_constrained, &
                        &vol_apply_even=solvent_even, vol_apply_odd=solvent_odd)
                else
                    call nonuniform_filter_state(params, state, half_even, half_odd, &
                        &ml_even, ml_odd, nu_aux_member(params), &
                        &res0143s(state), fname_vol, eonames, align_lps(state), timings=bt%nu, &
                        &base_support=state_support_msk, l_base_constrained=l_nu_base_constrained)
                endif
                time_nu_filter = real(toc(t_state_phase),dp)
            endif
            call write_output_diagnostics(state, 'distributed', time_map_output, time_fsc_output, time_nu_filter)
            if( L_BENCH_GLOB )then
                bt%rt_map_output  = bt%rt_map_output + time_map_output
                bt%rt_fsc         = bt%rt_fsc        + time_fsc_output
                bt%rt_nu_filter   = bt%rt_nu_filter  + time_nu_filter
                bt%nstates_solved = bt%nstates_solved + 1
            endif
            params%vols(state)      = fname_vol
            params%vols_even(state) = fname_even
            params%vols_odd(state)  = fname_odd
            call cline%set('vol'//int2str(state), fname_vol)
            state_written(state) = .true.
            call half_even%kill
            call half_odd%kill
            call solvent_even%kill
            call solvent_odd%kill
            call atm_even%kill
            call atm_odd%kill
            call atm_res%kill
            call atm_prior%kill
            if( allocated(atm_c) ) deallocate(atm_c)
            if( params%l_ml_reg )then
                call ml_even%kill
                call ml_odd%kill
                call fname_even_unfil%kill
                call fname_odd_unfil%kill
            endif
            call merged%kill
            call solvent_weight(1)%kill
            call solvent_weight(2)%kill
            call fname_even%kill
            call fname_odd%kill
            call fname_vol%kill
            call fname_fsc%kill
            if( allocated(fsc) ) deallocate(fsc)
        enddo
        if( present(nu_align_lps) ) nu_align_lps = align_lps
        if( .not. any(state_written .or. state_carried) ) THROW_HARD('distributed PCG produced no populated states')
        if( L_BENCH_GLOB ) t_bench = tic()
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res',   res0143s(1))
            call build%spproj_field%set_all2single('res05', res05s(1))
            call build%spproj_field%set_all2single('cfar',  cfars(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) )then
                        call build%spproj_field%set(iptcl, 'res',   res0143s(istate))
                        call build%spproj_field%set(iptcl, 'res05', res05s(istate))
                        call build%spproj_field%set(iptcl, 'cfar',  cfars(istate))
                    endif
                endif
            enddo
        endif
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        if( L_BENCH_GLOB )then
            bt%rt_proj_write = toc(t_bench)
            t_bench = tic()
        endif
        ! Delete raw artifacts only after every state has completed. Until this
        ! point they remain a restart/debug boundary for a failed master solve.
        do state = 1, params%nstates
            do eo = 0, 1
                do part = 1, params%nparts
                    raw_fname = refine3D_pcg_raw_accum_fname(state, part, params%numlen, &
                        &merge('odd ', 'even', eo == 1))
                    call del_file(raw_fname)
                enddo
            enddo
        enddo
        if( .not. params%l_trail_rec .and. .not. pcg_trail_seed_requested(cline) )then
            do state = 1, params%nstates
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'even')
                call del_file(raw_fname)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'odd')
                call del_file(raw_fname)
            enddo
        endif
        call raw_fname%kill
        call state_support_msk%kill_bimg
        call wset%kill
        ! operators a skipped state kept for its regularized pair
        call even_job%pcgop%kill
        call odd_job%pcgop%kill
        deallocate(res0143s, res05s, cfars, state_written, state_carried, state_empty, realized_fractions, update_weights,&
            &align_lps)
        deallocate(chain_weights, current_scales)
        if( allocated(nrep) ) deallocate(nrep, nsmp)
        !$ call omp_set_num_threads(params%nthr)
        if( L_BENCH_GLOB ) bt%rt_cleanup = toc(t_bench)
        if( present(bench) ) bench = bt

    contains

        !> Population-rule weights of one state's chain pair (class-average and
        !! reconstruct3D partials note, Section 4.1). The chain's represented
        !! population M(s) is read from the pair's manifest (a mass under a state
        !! weight set); with N(s) active
        !! updated rows and n(s) in the current sample, current *= s = u/f and
        !! chain *= w = (1-u)*N/M, so the blended mass is N whatever joined or left
        !! the state, and the current-map coefficient stays u. With N = M this is
        !! the former chain weight 1 - u.
        subroutine set_chain_blend_weights( state_here )
            integer, intent(in) :: state_here
            type(trail_chain_manifest) :: man
            type(string) :: manifest
            real    :: mnew, mrep
            integer :: status, nnew
            ! M: the mass the chain pair represents, from its manifest (validated by discard_stale_trail_chain_pair)
            manifest = refine3D_pcg_trail_manifest_fname(state_here)
            call man%read(manifest, status)
            if( status /= TRAIL_MANIFEST_OK ) THROW_HARD('unreadable PCG trailing chain manifest')
            mrep = man%get_mrep()
            call man%kill
            call manifest%kill
            call population_blend_weights(rnrep(state_here), rnsmp(state_here), mrep, &
                &current_scales(state_here), chain_weights(state_here), mnew, ufrac=update_weights(state_here))
            nnew = count_first_time(state_here)
            write(logfhandle,'(A,I0,A,F8.4,A,F8.4,A,F8.4,A,4F12.2,A,I9)') '>>> PCG TRAILING BLEND, STATE ', state_here, &
                &', PREVIOUS-CHAIN WEIGHT ', chain_weights(state_here), ', CURRENT SCALE ', current_scales(state_here), &
                &', FORMER WEIGHT 1-U ', 1.0 - update_weights(state_here), ', N n M MNEW', rnrep(state_here), &
                &rnsmp(state_here), mrep, mnew, ', FIRST-TIME', nnew
        end subroutine set_chain_blend_weights

        !> first-time rows (updatecnt = 1) of the current sample of a state
        integer function count_first_time( state_here ) result( nnew )
            integer, intent(in) :: state_here
            integer :: p, sample_ind
            sample_ind = 0
            do p = 1, build%spproj_field%get_noris()
                sample_ind = max(sample_ind, build%spproj_field%get_sampled(p))
            enddo
            nnew = 0
            do p = 1, build%spproj_field%get_noris()
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_sampled(p) /= sample_ind ) cycle
                if( build%spproj_field%get_updatecnt(p) == 1 ) nnew = nnew + 1
            enddo
        end function count_first_time

        !> Discard a chain pair of another identity: provenance, field of view, a
        !! larger crop than the current, or an unreadable or old format; constant-
        !! FOV crop growth is not stale (zero-extension on read). Both halves go
        !! together; the caller re-seeds the chain from the current sample.
        subroutine discard_stale_trail_chain_pair( state_here )
            integer, intent(in) :: state_here
            type(string) :: even_fname, odd_fname
            logical      :: l_even_here, l_odd_here, l_stale
            even_fname  = refine3D_pcg_trail_accum_fname(state_here, 'even')
            odd_fname   = refine3D_pcg_trail_accum_fname(state_here, 'odd')
            l_even_here = file_exists(even_fname)
            l_odd_here  = file_exists(odd_fname)
            l_stale     = .false.
            if( l_even_here ) l_stale = .not. pcg_raw_accum_compatible(even_fname, &
                &params%box_crop, params%smpd_crop, chain_provenance)
            if( .not. l_stale .and. l_odd_here ) l_stale = .not. pcg_raw_accum_compatible(odd_fname, &
                &params%box_crop, params%smpd_crop, chain_provenance)
            ! the manifest publishes the pair: missing, unreadable, of another weight set or of other sizes = stale
            if( .not. l_stale .and. (l_even_here .or. l_odd_here) ) l_stale = .not. chain_manifest_valid(state_here)
            if( l_stale )then
                if( l_even_here ) call del_file(even_fname)
                if( l_odd_here  ) call del_file(odd_fname)
                even_fname = refine3D_pcg_trail_manifest_fname(state_here)
                call del_file(even_fname)
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: DISCARDING STALE TRAILING CHAIN, STATE ', &
                    &state_here, ' (GEOMETRY/IDENTITY CHANGE); RE-SEEDING FROM THE CURRENT SAMPLE'
            endif
            call even_fname%kill
            call odd_fname%kill
        end subroutine discard_stale_trail_chain_pair

        !> The states without particles in this iteration (no particle in the raw accumulator headers
        !! of any part, either half) while another state has some. None in a single-state run or a
        !! state= run.
        subroutine find_empty_states( l_empty )
            logical, allocatable, intent(out) :: l_empty(:)
            integer :: nptcls_state(params%nstates), istate_here, part_here, eo_here, status
            integer :: hstate, heo, hpart, hnparts, hnptcls, hbox
            real    :: hsmpd
            character(len=256) :: hprov
            type(string) :: fname
            allocate(l_empty(params%nstates), source=.false.)
            if( params%nstates <= 1 .or. params%l_state_defined ) return
            nptcls_state = 0
            do istate_here = 1, params%nstates
                do eo_here = 0, 1
                    do part_here = 1, params%nparts
                        fname = refine3D_pcg_raw_accum_fname(istate_here, part_here, params%numlen, &
                            &merge('odd ', 'even', eo_here == 1))
                        call read_pcg_raw_accum_header(fname, hstate, heo, hpart, hnparts, hnptcls, hbox, hsmpd, &
                            &hprov, status)
                        if( status == 0 ) nptcls_state(istate_here) = nptcls_state(istate_here) + max(0, hnptcls)
                    enddo
                enddo
            enddo
            call fname%kill
            if( .not. any(nptcls_state > 0) ) return
            l_empty = nptcls_state == 0
        end subroutine find_empty_states

        !> the shipped volume and half maps a carried-forward state keeps on disk
        logical function previous_state_maps_exist( state_here ) result( l_exists )
            integer, intent(in) :: state_here
            l_exists = file_exists(refine3D_state_vol_fname(state_here)) .and. &
                &file_exists(refine3D_state_halfvol_fname(state_here, 'even')) .and. &
                &file_exists(refine3D_state_halfvol_fname(state_here, 'odd')) .and. &
                &file_exists(refine3D_fsc_fname(state_here))
        end function previous_state_maps_exist

        !> contributors of the full population of a state half (the rows a chain represents)
        integer function count_full_state_half( state_here, eo_here ) result(n)
            integer, intent(in) :: state_here, eo_here
            real(dp) :: mass
            call full_state_half(state_here, eo_here, n, mass)
        end function count_full_state_half

        !> contributors and applied mass of the full population of a state half: the rows of the half with
        !! the state (label, or weight above PCG_WEIGHT_THRESHOLD), updated once any row was updated
        subroutine full_state_half( state_here, eo_here, n, mass )
            integer,  intent(in)  :: state_here, eo_here
            integer,  intent(out) :: n
            real(dp), intent(out) :: mass
            integer, allocatable :: rows(:), members(:)
            real,    allocatable :: wmem(:)
            logical, allocatable :: l_row(:)
            integer :: p
            n    = 0
            mass = 0.0_dp
            if( l_weighted )then
                allocate(l_row(params%fromp:params%top), source=.false.)
                do p = params%fromp, params%top
                    if( build%spproj_field%get_state(p) <= 0 ) cycle
                    if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                    if( l_has_updates .and. build%spproj_field%get_updatecnt(p) < 1 ) cycle
                    l_row(p) = .true.
                enddo
                rows = pack([(p, p = params%fromp, params%top)], l_row)
                call wset%get_members(state_here, PCG_WEIGHT_THRESHOLD, rows, members, wmem)
                n    = size(members)
                mass = sum(real(wmem,dp))
                deallocate(rows, members, wmem, l_row)
                return
            endif
            do p = params%fromp, params%top
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                if( l_has_updates .and. build%spproj_field%get_updatecnt(p) < 1 ) cycle
                n = n + 1
            enddo
            mass = real(n, dp)
        end subroutine full_state_half

        !> write one half of the chain pair at the full population's contributors and mass; the first half
        !! deletes the pair's manifest so an interrupted write never validates
        subroutine write_chain_half( pcgop, fname, state_here, eo_here )
            class(reconstructor_pcg), intent(inout) :: pcgop
            class(string),            intent(in)    :: fname
            integer,                  intent(in)    :: state_here, eo_here
            type(string) :: manifest
            integer  :: n
            real(dp) :: mass
            if( .not. any(chain_written) )then
                manifest = refine3D_pcg_trail_manifest_fname(state_here)
                call del_file(manifest)
                call manifest%kill
            endif
            call full_state_half(state_here, eo_here, n, mass)
            call pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, n, chain_provenance, mass=mass, wset_id=wset_id)
            chain_written(eo_here+1) = .true.
            chain_contrib(eo_here+1) = n
            chain_mass(eo_here+1)    = mass
        end subroutine write_chain_half

        !> publish the chain pair of a state: represented mass, contributors, identity and component sizes
        subroutine write_chain_manifest( state_here )
            integer, intent(in) :: state_here
            type(trail_chain_manifest) :: man
            type(string)    :: manifest, half_fname
            integer(int64)  :: sizes(4)
            integer         :: io_stat, gen
            integer(timer_int_kind) :: t_manifest
            if( L_BENCH_GLOB ) t_manifest = tic()
            sizes = 0_int64
            half_fname = refine3D_pcg_trail_accum_fname(state_here, 'even')
            inquire(file=half_fname%to_char(), size=sizes(1))
            half_fname = refine3D_pcg_trail_accum_fname(state_here, 'odd')
            inquire(file=half_fname%to_char(), size=sizes(2))
            call half_fname%kill
            manifest = refine3D_pcg_trail_manifest_fname(state_here)
            gen = chain_gen_prev + 1
            call man%new(params%box_crop, params%smpd_crop, build%spproj_field%get_noris(), params%nstates, state_here, &
                &gen, sizes, real(sum(chain_mass)), sum(chain_contrib), wset_id)
            call man%write(manifest, io_stat)
            if( io_stat /= 0 ) THROW_WARN('failed to write the PCG trailing chain manifest; the chain will be re-seeded')
            call man%kill
            call manifest%kill
            if( L_BENCH_GLOB ) bt%rt_chain = bt%rt_chain + toc(t_manifest)
        end subroutine write_chain_manifest

        !> the pair's manifest exists, parses, belongs to this weight set and state layout, and matches both sizes
        logical function chain_manifest_valid( state_here ) result( l_valid )
            integer, intent(in) :: state_here
            type(trail_chain_manifest) :: man
            type(string)   :: manifest, half_fname
            integer(int64) :: sz_even, sz_odd
            integer        :: status
            l_valid  = .false.
            manifest = refine3D_pcg_trail_manifest_fname(state_here)
            call man%read(manifest, status)
            call manifest%kill
            if( status /= TRAIL_MANIFEST_OK ) return
            sz_even = -1_int64
            sz_odd  = -1_int64
            half_fname = refine3D_pcg_trail_accum_fname(state_here, 'even')
            if( file_exists(half_fname) ) inquire(file=half_fname%to_char(), size=sz_even)
            half_fname = refine3D_pcg_trail_accum_fname(state_here, 'odd')
            if( file_exists(half_fname) ) inquire(file=half_fname%to_char(), size=sz_odd)
            call half_fname%kill
            l_valid = all(man%get_wset_id() == wset_id) .and. man%get_state() == state_here .and. &
                &man%get_nstates() == params%nstates .and. man%get_size(1) == sz_even .and. man%get_size(2) == sz_odd
            ! a blended chain continues the generation count of the chain it replaces; a re-seed restarts at one
            if( l_valid ) chain_gen_prev = man%get_gen()
            call man%kill
        end function chain_manifest_valid

        !> check_solvent_lambda_by_resolve through the prepared half jobs, which
        !! are left zeroed with the production strength installed
        subroutine check_solvent_lambda_by_resolve_distributed( state_here, even_job, odd_job, base_even, base_odd )
            integer,                    intent(in)    :: state_here
            type(distributed_half_job), intent(inout) :: even_job, odd_job
            type(image),                intent(in)    :: base_even, base_odd
            type(image) :: cand_even, cand_odd, support
            real    :: j, jref
            integer :: ig
            call build_solvent_check_support(params, state_support_msk, l_state_support, support)
            call cand_even%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call cand_odd%new( [params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            jref = solvent_prior_cross_half_objective(base_even, base_odd, base_even, base_odd, support)
            write(logfhandle,'(A,I0,A)') '>>> PCG SOLVENT PRIOR LAMBDA CHECK: STATE ', state_here, &
                &', the same grid by real re-solves (pcg_solvent_check=yes)'
            write(logfhandle,'(A)') '    lambda_rel     J/J(0) (re-solve)   RESID even/odd   MRES even/odd'
            do ig = 1, size(PCG_SOLVENT_LAMBDA_GRID)
                call even_job%pcgop%set_solvent_prior(solvent_weight(1), PCG_SOLVENT_LAMBDA_GRID(ig))
                call odd_job%pcgop%set_solvent_prior( solvent_weight(2), PCG_SOLVENT_LAMBDA_GRID(ig))
                even_job%x = 0.0
                odd_job%x  = 0.0
                even_job%l_nonzero = .false.
                odd_job%l_nonzero  = .false.
                call solve_distributed_half_pair(even_job, odd_job)
                call cand_even%set_rmat(even_job%x, .false.)
                call cand_odd%set_rmat( odd_job%x,  .false.)
                j = solvent_prior_cross_half_objective(base_even, base_odd, cand_even, cand_odd, support)
                write(logfhandle,'(A,F8.3,A,F10.4,A,2(1X,ES9.3),A,2(1X,ES9.3))') '     ', PCG_SOLVENT_LAMBDA_GRID(ig), &
                    &'  ', j / max(TINY, jref), '        ', even_job%result%final_rel_residual, &
                    &odd_job%result%final_rel_residual, '  ', even_job%result%final_rel_residual_m, &
                    &odd_job%result%final_rel_residual_m
            enddo
            ! restore the production strength, leave the jobs cold
            call even_job%pcgop%set_solvent_prior(solvent_weight(1), solvent_lambda_eff)
            call odd_job%pcgop%set_solvent_prior( solvent_weight(2), solvent_lambda_eff)
            even_job%x = 0.0
            odd_job%x  = 0.0
            even_job%l_nonzero = .false.
            odd_job%l_nonzero  = .false.
            call cand_even%kill
            call cand_odd%kill
            call support%kill
        end subroutine check_solvent_lambda_by_resolve_distributed

        !> atm_prior=yes, between the prior-free base solve and its re-solve: the model of pdbfile on the solve support,
        !! its kernel width fitted on the merged base map (or bfac_data), the strength from the base pair's moments (or
        !! lambda_atm), amplitudes fitted per half through that half's operator, then the prior on both operators
        subroutine prepare_atomic_prior_on_pair( state_here )
            integer, intent(in) :: state_here
            type(atoms)       :: model
            type(pcg_atomic_prior_stats) :: astats
            real, allocatable :: pos(:,:), supp(:,:,:), xb(:,:,:), model_map(:,:,:)
            real, allocatable :: atom_bfac(:), atom_occ(:), atom_prior_intensity(:)
            real    :: bdata, lambda_rel_est, shift_a
            integer :: i, nsh
            logical :: l_fitted_model
            call model%new(params%pdbfile)
            if( model%get_n() < 1 ) THROW_HARD('the atom model pdbfile holds no atoms; atm_prior')
            ! PDB A = 0-based native voxel * smpd; the crop keeps the box centre and the field of view
            shift_a = real(params%box/2) * params%smpd
            allocate(pos(3,model%get_n()))
            do i = 1, model%get_n()
                pos(:,i) = (model%get_coord(i) - shift_a) / params%smpd_crop + real(params%box_crop/2 + 1)
            end do
            call read_species_prior_fields(model%get_n(), atom_bfac, atom_occ, atom_prior_intensity, l_fitted_model)
            call model%kill
            call even_job%pcgop%get_solve_support(supp)
            if( l_fitted_model )then
                call atm_prior%new(params%box_crop, params%smpd_crop, pos, supp, bfacs=atom_bfac,&
                    &occupancies=atom_occ, intensities=atom_prior_intensity)
            else
                call atm_prior%new(params%box_crop, params%smpd_crop, pos, supp)
            endif
            deallocate(pos, supp)
            xb  = 0.5 * (even_job%x + odd_job%x)
            nsh = 0
            if( l_fitted_model )then
                bdata = 0.
            elseif( params%l_bfac_data_fit )then
                call atm_prior%fit_width(xb, bdata, nsh)
            else
                bdata = params%bfac_data
            endif
            if( .not. l_fitted_model )then
                if( bdata > 0. ) call atm_prior%build(atm_prior%get_bfac_min() + bdata)
            endif
            lambda_rel_est = atm_prior%estimate_lambda(xb, even_job%x, odd_job%x)
            atm_lambda_rel = lambda_rel_est
            if( params%lambda_atm > 0. ) atm_lambda_rel = params%lambda_atm
            astats = atm_prior%get_stats()
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') '>>> PCG ATOMIC PRIOR: STATE ', state_here, ', ', &
                &astats%natoms, ' of ', astats%natoms_in, ' atoms (', astats%ndup, ' duplicates, ', astats%noutside, &
                &' off the support)'
            if( astats%fitted_model )then
                write(logfhandle,'(A,3F8.3,A,F7.4,A,ES11.3)') '>>> PCG ATOMIC PRIOR: fitted B min/mean/max ', &
                    &astats%bfac_model_min, astats%bfac_model_mean, astats%bfac_model_max, &
                    &' A^2, mean occupancy ', astats%occupancy_mean, ', mean posterior intensity ', astats%intensity_mean
            else
                write(logfhandle,'(A,F8.3,A,F8.3,A,I0,A,F6.3,A,I0)') '>>> PCG ATOMIC PRIOR: B_min ', astats%bfac_min, &
                    &' A^2 + B_data ', astats%bfac_data, ' A^2 (fit over ', nsh, ' shells), sigma ', astats%sigma_vox, &
                    &' voxels, footprint radius ', astats%foot_rad
            endif
            write(logfhandle,'(A,4ES11.3,A,F7.4,A,F9.5,A)') '>>> PCG ATOMIC PRIOR: N S E C^2 ', astats%noise, &
                &astats%signal, astats%explained, astats%c2, ' -> lambda_rel estimate ', astats%lambda_rel_est, &
                &' (merged-map ratio ', astats%lambda_rel_merged, ')'
            ! amplitudes per half through its own operator, before the prior is attached (it vanishes on span(G))
            allocate(atm_c(atm_prior%get_natoms(),2), source=0._dp)
            call fit_amplitudes_pair
            if( atm_prior%uses_fitted_model() )then
                model_map = atm_prior%project(xb)
            else
                model_map = 0.5 * (even_job%pcgop%atomic_model_map(atm_prior, atm_c(:,1)) + &
                    &odd_job%pcgop%atomic_model_map(atm_prior, atm_c(:,2)))
            endif
            call atm_res%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call atm_res%set_rmat(xb - model_map, .false.)
            deallocate(xb, model_map)
            write(logfhandle,'(A,2(I0,A,ES10.3,A))') '>>> PCG ATOMIC PRIOR: free amplitudes fitted for diagnostics, even ', &
                &atm_amp_nits(1), ' iterations (residual ', atm_amp_rres(1), '), odd ', atm_amp_nits(2), &
                &' iterations (residual ', atm_amp_rres(2), ')'
            call even_job%pcgop%set_atomic_prior(atm_prior, atm_lambda_rel)
            call odd_job%pcgop%set_atomic_prior(atm_prior, atm_lambda_rel)
            atm_lambda_abs = [even_job%pcgop%get_atomic_lambda(), odd_job%pcgop%get_atomic_lambda()]
            write(logfhandle,'(A,F9.5,A,2ES11.3,A,2ES11.3)') '>>> PCG ATOMIC PRIOR: lambda_rel ', atm_lambda_rel, &
                &' | lambda_A even/odd ', atm_lambda_abs, ' | realspace diagonal even/odd ', &
                &even_job%pcgop%get_realspace_diagonal(), odd_job%pcgop%get_realspace_diagonal()
        end subroutine prepare_atomic_prior_on_pair

        !> The species CSV beside its PDB carries fitted B, occupancy and posterior-weighted prior intensity.
        subroutine read_species_prior_fields( natoms, bfacs, occupancies, prior_intensities, found )
            integer,           intent(in)  :: natoms
            real, allocatable, intent(out) :: bfacs(:), occupancies(:), prior_intensities(:)
            logical,           intent(out) :: found
            type(string) :: csvfile
            character(len=STDLEN) :: line
            real, allocatable :: row(:)
            integer :: funit, io_stat, i, j, ncol
            found   = .false.
            csvfile = get_fbody(params%pdbfile, 'pdb')//'.csv'
            if( .not. file_exists(csvfile) ) return
            call fopen(funit, file=csvfile, status='old', action='read', iostat=io_stat)
            call fileiochk('atomic prior opening fitted model table '//csvfile%to_char(), io_stat)
            read(funit,'(A)',iostat=io_stat) line
            call fileiochk('atomic prior reading fitted model table header', io_stat)
            if( index(trim(line), 'MODEL_INTENSITY,POSTERIOR_INTENSITY,OCCUPANCY') == 0 )then
                call fclose(funit)
                return
            endif
            ncol = count([(line(j:j) == ',', j=1,len_trim(line))]) + 1
            allocate(row(ncol), bfacs(natoms), occupancies(natoms), prior_intensities(natoms))
            do i = 1, natoms
                read(funit,'(A)',iostat=io_stat) line
                call fileiochk('atomic prior reading fitted model table row', io_stat)
                do j = 1, len_trim(line)
                    if( line(j:j) == ',' ) line(j:j) = ' '
                enddo
                read(line,*,iostat=io_stat) row
                call fileiochk('atomic prior parsing fitted model table row', io_stat)
                bfacs(i)      = row(13)
                prior_intensities(i) = row(ncol-1)
                occupancies(i) = row(ncol)
            enddo
            read(funit,'(A)',iostat=io_stat) line
            if( io_stat == 0 ) THROW_HARD('fitted model table has more rows than its PDB; atomic prior')
            call fclose(funit)
            if( any(.not. ieee_is_finite(bfacs)) .or. any(.not. ieee_is_finite(occupancies)) .or. &
                &any(.not. ieee_is_finite(prior_intensities)) ) THROW_HARD('fitted atomic model fields must be finite')
            found = .true.
        end subroutine read_species_prior_fields

        !> the per-half amplitude fits, concurrently as the half solves
        subroutine fit_amplitudes_pair
            real(dp), allocatable :: c_even(:), c_odd(:)
            integer :: previous_max_active_levels
            if( pcg_master_nthreads >= 2 )then
                previous_max_active_levels = 1
                !$ previous_max_active_levels = omp_get_max_active_levels()
                !$ call omp_set_max_active_levels(max(2, previous_max_active_levels))
                !$omp parallel sections num_threads(2) default(shared)
                !$omp section
                !$ call omp_set_num_threads(pcg_half_nthreads)
                call even_job%pcgop%fit_atom_amplitudes(atm_prior, c_even, ATM_AMP_MAXITS, ATM_AMP_RTOL, &
                    &atm_amp_nits(1), atm_amp_rres(1))
                !$omp section
                !$ call omp_set_num_threads(pcg_half_nthreads)
                call odd_job%pcgop%fit_atom_amplitudes(atm_prior, c_odd, ATM_AMP_MAXITS, ATM_AMP_RTOL, &
                    &atm_amp_nits(2), atm_amp_rres(2))
                !$omp end parallel sections
                !$ call omp_set_max_active_levels(previous_max_active_levels)
                !$ call omp_set_num_threads(pcg_master_nthreads)
            else
                call even_job%pcgop%fit_atom_amplitudes(atm_prior, c_even, ATM_AMP_MAXITS, ATM_AMP_RTOL, &
                    &atm_amp_nits(1), atm_amp_rres(1))
                call odd_job%pcgop%fit_atom_amplitudes(atm_prior, c_odd, ATM_AMP_MAXITS, ATM_AMP_RTOL, &
                    &atm_amp_nits(2), atm_amp_rres(2))
            endif
            atm_c(:,1) = c_even
            atm_c(:,2) = c_odd
        end subroutine fit_amplitudes_pair

        !> atm_prior=yes: the prior-free pair and its average under _base, the residual to the active atom model under _res, the
        !! sidecar _atmprior.txt; half_even/half_odd and merged still hold the base pair here
        subroutine write_atomic_prior_products( state_here )
            integer, intent(in) :: state_here
            type(pcg_atomic_prior_stats) :: astats
            type(string) :: fname
            integer :: funit
            fname = add2fbody(refine3D_state_halfvol_fname(state_here, 'even'), MRC_EXT, '_base')
            call half_even%write(fname, del_if_exists=.true.)
            fname = add2fbody(refine3D_state_halfvol_fname(state_here, 'odd'), MRC_EXT, '_base')
            call half_odd%write(fname, del_if_exists=.true.)
            fname = add2fbody(refine3D_state_vol_fname(state_here), MRC_EXT, '_base')
            call merged%write(fname, del_if_exists=.true.)
            fname = add2fbody(refine3D_state_vol_fname(state_here), MRC_EXT, '_res')
            call atm_res%write(fname, del_if_exists=.true.)
            astats = atm_prior%get_stats()
            fname  = get_fbody(refine3D_state_vol_fname(state_here), 'mrc')//'_atmprior.txt'
            call fopen(funit, file=fname, status='replace', action='write')
            write(funit,'(A,A)')      'pdbfile=',            params%pdbfile%to_char()
            write(funit,'(A,I0)')     'natoms_in=',          astats%natoms_in
            write(funit,'(A,I0)')     'natoms=',             astats%natoms
            write(funit,'(A,I0)')     'natoms_duplicate=',   astats%ndup
            write(funit,'(A,I0)')     'natoms_off_support=', astats%noutside
            write(funit,'(A,I0)')     'support_voxels=',     astats%nsupport
            write(funit,'(A,F10.4)')  'bfac_min=',           astats%bfac_min
            write(funit,'(A,F10.4)')  'bfac_data=',          astats%bfac_data
            write(funit,'(A,L1)')     'bfac_data_fitted=',   params%l_bfac_data_fit .and. .not. astats%fitted_model
            write(funit,'(A,L1)')     'fitted_model=',       astats%fitted_model
            write(funit,'(A,F10.4)')  'bfac_model_min=',     astats%bfac_model_min
            write(funit,'(A,F10.4)')  'bfac_model_mean=',    astats%bfac_model_mean
            write(funit,'(A,F10.4)')  'bfac_model_max=',     astats%bfac_model_max
            write(funit,'(A,F10.6)')  'occupancy_mean=',     astats%occupancy_mean
            write(funit,'(A,ES14.6)') 'posterior_intensity_mean=', astats%intensity_mean
            write(funit,'(A,I0)')     'width_fit_shells=',   astats%nshells_fit
            write(funit,'(A,F10.4)')  'sigma_voxels=',       astats%sigma_vox
            write(funit,'(A,I0)')     'footprint_radius=',   astats%foot_rad
            write(funit,'(A,I0)')     'gram_nnz=',           astats%gram_nnz
            write(funit,'(A,ES14.6)') 'N=',                  astats%noise
            write(funit,'(A,ES14.6)') 'S=',                  astats%signal
            write(funit,'(A,ES14.6)') 'E=',                  astats%explained
            write(funit,'(A,F10.6)')  'C2=',                 astats%c2
            write(funit,'(A,F12.6)')  'lambda_rel_merged=',  astats%lambda_rel_merged
            write(funit,'(A,F12.6)')  'lambda_rel_estimate=', astats%lambda_rel_est
            write(funit,'(A,F12.6)')  'lambda_rel=',         atm_lambda_rel
            write(funit,'(A,L1)')     'lambda_rel_given=',   params%lambda_atm > 0.
            write(funit,'(A,2ES14.6)') 'lambda_A_even_odd=', atm_lambda_abs
            write(funit,'(A,2I6)')    'amplitude_its_even_odd=',   atm_amp_nits
            write(funit,'(A,2ES14.6)') 'amplitude_resid_even_odd=', atm_amp_rres
            write(funit,'(A,2I6)')    'solve_its_even_odd=',       atm_nits
            write(funit,'(A,2ES14.6)') 'solve_resid_l2_even_odd=',  atm_resid
            write(funit,'(A,2ES14.6)') 'solve_resid_m_even_odd=',   atm_mres
            write(funit,'(A,I0)')     'maxits_atm=',         params%maxits_atm
            write(funit,'(A,ES14.6)') 'rtol_atm=',           params%rtol_atm
            call fclose(funit)
            call fname%kill
            write(logfhandle,'(A,I0,A)') '>>> PCG ATOMIC PRIOR: STATE ', state_here, ', prior-free pair written '//&
                &'under _base (the FSC is theirs), residual to the active atom model under _res, prior''d pair ships'
        end subroutine write_atomic_prior_products

        function atomic_prior_provenance() result( prov )
            character(len=:), allocatable :: prov
            prov = 'lambda_rel='//trim(real2str(atm_lambda_rel))//' natoms='//int2str(atm_prior%get_natoms())//&
                &' bfac='//trim(real2str(atm_prior%get_bfac()))//' fitted_model='//&
                &merge('yes','no ',atm_prior%uses_fitted_model())
        end function atomic_prior_provenance

        !> pcg_solvent=yes: solvent_even/odd receive the ridge re-solve, even/odd
        !! the prior-free base pair (FSC, NU, evidence, _unfil)
        subroutine reduce_solve_state_pair( state_here, even, odd, n_even_here, n_odd_here, solve_kind, &
                &fsc_prior, warm_even, warm_odd, solvent_even, solvent_odd, atm_even, atm_odd )
            integer,          intent(in)    :: state_here
            character(len=*), intent(in)    :: solve_kind
            type(image),      intent(inout) :: even, odd
            integer,          intent(out)   :: n_even_here, n_odd_here
            real, optional,   intent(in)    :: fsc_prior(:)
            type(image), optional, intent(in) :: warm_even, warm_odd
            type(image), optional, intent(inout) :: solvent_even, solvent_odd
            type(image), optional, intent(inout) :: atm_even, atm_odd !< atm_prior=yes: the prior'd pair
            integer(timer_int_kind) :: t_pair
            integer :: ikind
            logical :: l_resolved
            l_resolved = .false.
            ikind = merge(2, 1, present(fsc_prior))

            if( present(fsc_prior) )then
                ! regularized pair: closed form of the base pair (warm_even/odd
                ! are the base solutions it is derived from)
                if( .not. present(warm_even) .or. .not. present(warm_odd) ) &
                    &THROW_HARD('distributed PCG regularized pair requires both base half maps')
                call prepare_distributed_half_job(state_here, 0, 'even', solve_kind, even_job, &
                    &fsc_prior, warm_even)
                call prepare_distributed_half_job(state_here, 1, 'odd', solve_kind, odd_job, &
                    &fsc_prior, warm_odd)
            else
                if( present(warm_even) .or. present(warm_odd) ) &
                    &THROW_HARD('distributed PCG base solve cannot take replay warm starts')
                call prepare_distributed_half_job(state_here, 0, 'even', solve_kind, even_job)
                call prepare_distributed_half_job(state_here, 1, 'odd', solve_kind, odd_job)
            endif
            n_even_here = even_job%nptcls
            n_odd_here  = odd_job%nptcls

            if( L_BENCH_GLOB ) t_pair = tic()
            call solve_distributed_half_pair(even_job, odd_job)
            if( L_BENCH_GLOB )then
                bt%rt_solve(ikind) = bt%rt_solve(ikind) + toc(t_pair)
                if( even_job%ready ) bt%rt_half_solve(1,ikind) = bt%rt_half_solve(1,ikind) + even_job%time_solve
                if( odd_job%ready  ) bt%rt_half_solve(2,ikind) = bt%rt_half_solve(2,ikind) + odd_job%time_solve
                bt%niters(:,ikind) = bt%niters(:,ikind) + [even_job%niters, odd_job%niters]
                t_pair = tic()
            endif

            if( .not. present(fsc_prior) .and. params%l_pcg_solvent )then
                ! pcg_solvent=yes: weights from the prior-free pair, the ridge on
                ! both operators, the same cold solve again into the solvent pair
                if( even_job%ready .and. odd_job%ready )then
                    call report_pcg_solve('DISTRIBUTED', state_here, 'even', 'pre', even_job%nptcls, &
                        &even_job%niters, even_job%time_solve, even_job%result)
                    call report_pcg_solve('DISTRIBUTED', state_here, 'odd', 'pre', odd_job%nptcls, &
                        &odd_job%niters, odd_job%time_solve, odd_job%result)
                    call prepare_solvent_prior_on_pair(params, state_here, even_job%pcgop, odd_job%pcgop, &
                        &even_job%x, odd_job%x, state_support_msk, l_state_support, &
                        &solvent_weight, l_solvent_weight, res0143_prior_free, solvent_lambda_eff)
                    if( l_solvent_weight )then
                        if( .not.(present(solvent_even) .and. present(solvent_odd)) ) &
                            &THROW_HARD('the solvent-prior re-solve needs its output pair; reduce_solve_state_pair')
                        ! the prior-free pair is the base pair; the re-solve
                        ! with the ridge goes to the solvent pair
                        call validate_solved_map(even_job%x, 'distributed', state_here, 'even', 'pre')
                        call validate_solved_map(odd_job%x,  'distributed', state_here, 'odd',  'pre')
                        call even%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                        call odd%new( [params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                        call even%set_rmat(even_job%x, .false.)
                        call odd%set_rmat( odd_job%x,  .false.)
                        even_job%x = 0.0
                        odd_job%x  = 0.0
                        even_job%l_nonzero = .false.
                        odd_job%l_nonzero  = .false.
                        if( params%l_pcg_solvent_check ) &
                            &call check_solvent_lambda_by_resolve_distributed(state_here, even_job, odd_job, even, odd)
                        call solve_distributed_half_pair(even_job, odd_job)
                        l_resolved = .true.
                    endif
                endif
            endif

            if( .not. present(fsc_prior) .and. params%l_atm_prior )then
                ! atm_prior=yes: the model on the prior-free pair, the prior on both operators, the cold solve
                ! again into the atm pair; the prior-free pair is the base pair (FSC, products under _base)
                if( even_job%ready .and. odd_job%ready )then
                    if( .not.(present(atm_even) .and. present(atm_odd)) ) &
                        &THROW_HARD('the atomicity-prior re-solve needs its output pair; reduce_solve_state_pair')
                    call validate_solved_map(even_job%x, 'distributed', state_here, 'even', 'base')
                    call validate_solved_map(odd_job%x,  'distributed', state_here, 'odd',  'base')
                    call report_pcg_solve('DISTRIBUTED', state_here, 'even', 'base', even_job%nptcls, &
                        &even_job%niters, even_job%time_solve, even_job%result)
                    call report_pcg_solve('DISTRIBUTED', state_here, 'odd', 'base', odd_job%nptcls, &
                        &odd_job%niters, odd_job%time_solve, odd_job%result)
                    call write_distributed_diagnostics(state_here, 'even', 'base', even_job%l_concurrent, &
                        &even_job%nptcls, even_job%result, even_job%rel_res_hist, even_job%time_reduce, &
                        &even_job%time_finalize, even_job%time_solve, even_job%pcgop%get_data_scale(), &
                        &even_job%pcgop%get_effective_lambda(), pcgop=even_job%pcgop)
                    call write_distributed_diagnostics(state_here, 'odd', 'base', odd_job%l_concurrent, &
                        &odd_job%nptcls, odd_job%result, odd_job%rel_res_hist, odd_job%time_reduce, &
                        &odd_job%time_finalize, odd_job%time_solve, odd_job%pcgop%get_data_scale(), &
                        &odd_job%pcgop%get_effective_lambda(), pcgop=odd_job%pcgop)
                    call even%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                    call odd%new( [params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
                    call even%set_rmat(even_job%x, .false.)
                    call odd%set_rmat( odd_job%x,  .false.)
                    call prepare_atomic_prior_on_pair(state_here)
                    even_job%x = 0.0
                    odd_job%x  = 0.0
                    even_job%l_nonzero  = .false.
                    odd_job%l_nonzero   = .false.
                    even_job%l_atm_solve = .true.
                    odd_job%l_atm_solve  = .true.
                    even_job%solve_kind = 'atm'
                    odd_job%solve_kind  = 'atm'
                    call solve_distributed_half_pair(even_job, odd_job)
                    atm_nits  = [even_job%niters, odd_job%niters]
                    atm_resid = [even_job%result%final_rel_residual, odd_job%result%final_rel_residual]
                    atm_mres  = [even_job%result%final_rel_residual_m, odd_job%result%final_rel_residual_m]
                    call finish_distributed_half_job(even_job, atm_even)
                    call finish_distributed_half_job(odd_job, atm_odd)
                    l_atm_resolved = .true.
                    if( L_BENCH_GLOB ) bt%rt_prior_resolve = bt%rt_prior_resolve + toc(t_pair)
                    return
                endif
            endif

            if( L_BENCH_GLOB )then
                bt%rt_prior_resolve = bt%rt_prior_resolve + toc(t_pair)
                t_pair = tic()
            endif
            if( present(fsc_prior) )then
                call finish_distributed_half_job(even_job, even, warm_even)
                call finish_distributed_half_job(odd_job, odd, warm_odd)
            else if( l_resolved )then
                call finish_distributed_half_job(even_job, solvent_even)
                call finish_distributed_half_job(odd_job, solvent_odd)
            else
                call finish_distributed_half_job(even_job, even)
                call finish_distributed_half_job(odd_job, odd)
            endif
            if( L_BENCH_GLOB ) bt%rt_finish = bt%rt_finish + toc(t_pair)
        end subroutine reduce_solve_state_pair

        ! Preparation is deliberately serial. It owns mask memoization, FFTW
        ! planning, raw-accumulator I/O, trail writes, prior attachment and
        ! initial-guess construction. Only fully prepared, half-owned operators
        ! cross the OpenMP sections boundary below.
        subroutine prepare_distributed_half_job( state_here, eo_here, half, solve_kind, job, &
                &fsc_prior, warm_start )
            integer,          intent(in)    :: state_here, eo_here
            character(len=*), intent(in)    :: half, solve_kind
            type(distributed_half_job), intent(inout) :: job
            real, optional,   intent(in)    :: fsc_prior(:)
            type(image), optional, intent(in) :: warm_start
            type(string) :: fname
            integer :: part_here, n_part, n_full_half
            real(dp) :: part_mass
            integer(int64) :: part_wset_id(2)
            integer(timer_int_kind) :: t_phase, t_bench_half
            real(dp) :: rt_sum_part, rt_sum_half
            real :: realized_fraction, update_weight, current_scale
            logical :: l_chain_exists, l_seed_chain

            job%state = state_here
            job%eo = eo_here
            job%half = half
            job%solve_kind = solve_kind
            job%nptcls = 0
            job%mass = 0.0_dp
            job%niters = 0
            job%ready = .false.
            job%l_atm_solve = .false.
            job%l_ml_solve = present(fsc_prior)
            if( job%l_ml_solve .neqv. present(warm_start) ) &
                &THROW_HARD('distributed PCG regularized pair requires both FSC and base map')
            if( allocated(job%x) ) deallocate(job%x)
            if( allocated(job%rel_res_hist) ) deallocate(job%rel_res_hist)
            job%prior_npositive = 0
            job%prior_positive_min = 0.0
            job%prior_positive_max = 0.0
            job%prior_to_khat_l1 = 0.0
            job%prior_to_khat_rms = 0.0
            if( job%l_ml_solve )then
                call prepare_regularized_half_job(job, fsc_prior, warm_start)
                return
            endif
            ! an operator kept by a state without a regularized pair is stale (new kills it)
            if( allocated(job%d_raw) ) deallocate(job%d_raw)
            job%l_keep_op = .false.

            if( L_BENCH_GLOB ) t_bench_half = tic()
            call job%pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA, &
                &fft_nthreads=pcg_half_nthreads)
            ! serial: the spherical mask memoizes coordinates at module scope
            call set_pcg_solve_support(job%pcgop, params, state_support_msk, l_state_support)
            call job%pcgop%begin_reduction
            if( L_BENCH_GLOB ) bt%rt_op_setup = bt%rt_op_setup + toc(t_bench_half)
            ! the reduction's summation and scaling are timed apart from its read I/O
            rt_sum_half = 0._dp
            t_phase = tic()
            if( params%l_trail_rec .and. l_chain_carry )then
                ! a carried chain is the whole of this iteration's data
                fname = refine3D_pcg_trail_accum_fname(state_here, half)
                if( L_BENCH_GLOB ) call add_raw_read_size(fname)
                call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                    &chain_provenance, 1.0, job%nptcls, mass=job%mass, rt_sum=rt_sum_half)
                call fname%kill
            else
                do part_here = 1, params%nparts
                    fname = refine3D_pcg_raw_accum_fname(state_here, part_here, params%numlen, half)
                    if( L_BENCH_GLOB ) call add_raw_read_size(fname)
                    call job%pcgop%add_raw_accum(fname, state_here, eo_here, part_here, params%nparts, &
                        &provenance, n_part, mass=part_mass, wset_id=part_wset_id, rt_sum=rt_sum_part)
                    if( any(part_wset_id /= wset_id) ) THROW_HARD('a raw PCG part belongs to another state weight set')
                    job%nptcls = job%nptcls + n_part
                    job%mass   = job%mass   + part_mass
                    rt_sum_half = rt_sum_half + rt_sum_part
                    call fname%kill
                enddo
            endif
            job%time_reduce = real(toc(t_phase),dp)
            if( L_BENCH_GLOB )then
                bt%rt_sum  = bt%rt_sum  + rt_sum_half
                bt%rt_read = bt%rt_read + max(0._dp, job%time_reduce - rt_sum_half)
                t_bench_half = tic()
            endif
            if( job%nptcls == 0 )then
                call job%pcgop%kill
                return
            endif
            if( .not. l_chain_carry )then
                realized_fraction = realized_fractions(state_here)
                update_weight = update_weights(state_here)
                l_seed_chain = pcg_trail_seed_requested(cline)
                ! the full population is what a chain represents: checked where a chain is written (a
                ! caller-owned weight table reconstructs without one, and the project cannot count it)
                if( params%l_trail_rec .or. l_seed_chain )then
                    n_full_half = count_full_state_half(state_here, eo_here)
                    if( n_full_half < job%nptcls ) &
                        &THROW_HARD('PCG current half population exceeds its full population')
                endif
                fname = refine3D_pcg_trail_accum_fname(state_here, half)
                l_chain_exists = file_exists(fname)
                if( params%l_trail_rec )then
                    if( realized_fraction <= 0.0 ) &
                        &THROW_HARD('PCG trailing update has zero realized state fraction')
                    if( l_chain_exists .and. 1.0-update_weight > 0.01 )then
                        ! population rule: set_chain_blend_weights
                        current_scale = current_scales(state_here)
                        call job%pcgop%scale_raw_accum(current_scale)
                        call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                            &chain_provenance, chain_weights(state_here), n_part)
                    else
                        current_scale = 1.0 / realized_fraction
                        call job%pcgop%scale_raw_accum(current_scale)
                    endif
                    call write_chain_half(job%pcgop, fname, state_here, eo_here)
                    if( .not. l_chain_exists .or. 1.0-update_weight <= 0.01 ) &
                        &call job%pcgop%scale_raw_accum(realized_fraction)
                    write(logfhandle,'(A,I0,A,A,A,F8.4,A,F8.4)') '>>> PCG TRAIL | STATE=', state_here, &
                        &' | HALF=', trim(half), ' | F=', realized_fraction, ' | U=', update_weight
                else if( l_seed_chain )then
                    call write_chain_half(job%pcgop, fname, state_here, eo_here)
                endif
                call fname%kill
            endif
            if( L_BENCH_GLOB ) bt%rt_chain = bt%rt_chain + toc(t_bench_half)

            t_phase = tic()
            if( params%l_ml_reg )then
                ! the regularized pair that follows is this operator, regularized (prepare_regularized_half_job)
                call job%pcgop%end_accum(.true., job%d_raw)
                job%l_keep_op = .true.
                ! what the reduction represents: the chain half written above, else the base reduction
                job%nptcls_raw = job%nptcls
                job%mass_raw   = job%mass
                if( params%l_trail_rec .and. chain_written(eo_here+1) )then
                    job%nptcls_raw = chain_contrib(eo_here+1)
                    job%mass_raw   = chain_mass(eo_here+1)
                endif
            else
                call job%pcgop%end_accum(.true.)
            endif
            call job%pcgop%set_op_mode(PCG_OP_KERNEL)
            job%time_finalize = real(toc(t_phase),dp)
            if( L_BENCH_GLOB ) bt%rt_finalize(1) = bt%rt_finalize(1) + job%time_finalize
            ! the base solve starts from zero
            allocate(job%x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            job%l_nonzero = .false.
            job%ready = .true.
        end subroutine prepare_distributed_half_job

        !> The regularized pair is its base pair's operator, kept by finish_distributed_half_job and
        !! regularized: the support in force now, P_tau from the base pair's FSC, the preconditioner
        !! rebuilt from the base reduction's raw density. Neither the raw artifacts nor the trailing
        !! chain are read again, and the RHS fold and Khat are not rebuilt
        subroutine prepare_regularized_half_job( job, fsc_prior, warm_start )
            type(distributed_half_job), intent(inout) :: job
            real,                       intent(in)    :: fsc_prior(:)
            type(image),                intent(in)    :: warm_start
            integer(timer_int_kind) :: t_phase
            if( .not. job%l_keep_op .or. .not. allocated(job%d_raw) ) &
                &THROW_HARD('the regularized pair needs the finalized operator of its base pair')
            job%l_keep_op   = .false.
            job%nptcls      = job%nptcls_raw
            job%mass        = job%mass_raw
            job%time_reduce = 0.0_dp
            if( L_BENCH_GLOB ) t_phase = tic()
            call set_pcg_solve_support(job%pcgop, params, state_support_msk, l_state_support)
            if( L_BENCH_GLOB ) bt%rt_op_setup = bt%rt_op_setup + toc(t_phase)
            t_phase = tic()
            call job%pcgop%set_ml_prior(fsc_prior, params%tau, params%hp)
            job%band_shell = get_find_at_crit(fsc_prior, 0.143)
            call job%pcgop%regularize_finalized(job%d_raw)
            call job%pcgop%set_op_mode(PCG_OP_KERNEL)
            job%time_finalize = real(toc(t_phase),dp)
            if( L_BENCH_GLOB ) bt%rt_finalize(2) = bt%rt_finalize(2) + job%time_finalize
            call job%pcgop%get_ml_prior_stats(job%prior_npositive, job%prior_positive_min, &
                &job%prior_positive_max, job%prior_to_khat_l1, job%prior_to_khat_rms)
            ! the regularized map is the closed-form P_tau optimum of the
            ! current base solution (shrink_by_ml_prior); nothing is solved
            job%x = warm_start%get_rmat()
            job%l_nonzero = .false.
            job%ready = .true.
        end subroutine prepare_regularized_half_job

        !> bench: one raw accumulator about to be read and its size
        subroutine add_raw_read_size( fname )
            class(string), intent(in) :: fname
            integer(int64) :: nbytes
            nbytes = -1_int64
            inquire(file=fname%to_char(), size=nbytes)
            bt%raw_bytes = bt%raw_bytes + max(0_int64, nbytes)
            bt%nraw      = bt%nraw + 1
        end subroutine add_raw_read_size

        subroutine solve_distributed_half_pair( even, odd )
            type(distributed_half_job), intent(inout) :: even, odd
            integer :: previous_max_active_levels

            even%l_concurrent = .false.
            odd%l_concurrent  = .false.
            if( even%ready .and. odd%ready .and. pcg_master_nthreads >= 2 )then
                even%l_concurrent = .true.
                odd%l_concurrent  = .true.
                previous_max_active_levels = 1
                !$ previous_max_active_levels = omp_get_max_active_levels()
                !$ call omp_set_max_active_levels(max(2, previous_max_active_levels))
                !$omp parallel sections num_threads(2) default(shared)
                !$omp section
                !$ call omp_set_num_threads(pcg_half_nthreads)
                call solve_prepared_half_job(even)
                !$omp section
                !$ call omp_set_num_threads(pcg_half_nthreads)
                call solve_prepared_half_job(odd)
                !$omp end parallel sections
                !$ call omp_set_max_active_levels(previous_max_active_levels)
                !$ call omp_set_num_threads(pcg_master_nthreads)
            else
                call solve_prepared_half_job(even)
                call solve_prepared_half_job(odd)
            endif
        end subroutine solve_distributed_half_pair

        subroutine solve_prepared_half_job( job )
            type(distributed_half_job), intent(inout) :: job
            integer(timer_int_kind) :: t_phase, t_end, t_rate
            if( .not. job%ready ) return
            call system_clock(count=t_phase)
            if( job%l_ml_solve .and. l_solvent_weight )then
                ! closed form, then coupled iterations with the soft solvent prior
                call job%pcgop%solve_regularized(job%x, params%maxits_ml, job%rel_res_hist, &
                    &job%niters, job%result, job%x_cf, solvent_weight=solvent_weight(job%eo+1), &
                    &solvent_lambda_rel=solvent_lambda_eff)
            else if( job%l_ml_solve )then
                ! closed form, then maxits_ml coupled iterations from it
                call job%pcgop%solve_regularized(job%x, params%maxits_ml, job%rel_res_hist, &
                    &job%niters, job%result, job%x_cf)
            else if( job%l_atm_solve )then
                call job%pcgop%solve_with_cold_restart(job%x, job%l_nonzero, params%maxits_atm, params%rtol_atm, &
                    &job%rel_res_hist, job%niters, job%result)
            else
                call job%pcgop%solve_with_cold_restart(job%x, job%l_nonzero, params%maxits_pcg, params%rtol, &
                    &job%rel_res_hist, job%niters, job%result)
            endif
            call system_clock(count=t_end, count_rate=t_rate)
            job%time_solve = real(t_end-t_phase,dp) / real(t_rate,dp)
        end subroutine solve_prepared_half_job

        ! Finalization is serial for deterministic logging, diagnostics and
        ! image/FFTW lifecycle management.
        subroutine finish_distributed_half_job( job, volume, warm_start )
            type(distributed_half_job), intent(inout) :: job
            type(image), intent(inout) :: volume
            type(image), optional, intent(in) :: warm_start
            if( .not. job%ready ) return
            call handle_cold_restart_outcome(job%result, 'distributed', job%half, job%solve_kind)
            call validate_solved_map(job%x, 'distributed', job%state, job%half, job%solve_kind)
            if( allocated(job%x_cf) )then
                call measure_closed_form_agreement(job%x_cf, job%x, params%box_crop, params%smpd_crop, &
                    &job%band_shell, job%result)
                deallocate(job%x_cf)
            endif
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(job%x, .false.)
            call report_beyond_band_excess(volume, params, job%state, job%half, job%solve_kind)
            if( job%l_ml_solve )then
                call write_distributed_diagnostics(job%state, job%half, job%solve_kind, job%l_concurrent, &
                    &job%nptcls, job%result, job%rel_res_hist, job%time_reduce, job%time_finalize, &
                    &job%time_solve, job%pcgop%get_data_scale(), job%pcgop%get_effective_lambda(), &
                    &job%prior_npositive, job%prior_positive_min, job%prior_positive_max, &
                    &job%prior_to_khat_l1, job%prior_to_khat_rms, job%pcgop)
            else
                call write_distributed_diagnostics(job%state, job%half, job%solve_kind, job%l_concurrent, &
                    &job%nptcls, job%result, job%rel_res_hist, job%time_reduce, job%time_finalize, &
                    &job%time_solve, job%pcgop%get_data_scale(), job%pcgop%get_effective_lambda(), &
                    &pcgop=job%pcgop)
            endif
            call report_pcg_solve('DISTRIBUTED', job%state, job%half, job%solve_kind, job%nptcls, &
                &job%niters, job%time_solve, job%result)
            call report_closed_form_agreement('DISTRIBUTED', job%state, job%half, job%result)
            if( .not. job%l_keep_op ) call job%pcgop%kill
            if( allocated(job%x) ) deallocate(job%x)
            if( allocated(job%x_cf) ) deallocate(job%x_cf)
            if( allocated(job%rel_res_hist) ) deallocate(job%rel_res_hist)
            job%ready = .false.
        end subroutine finish_distributed_half_job

        subroutine write_distributed_diagnostics( state_here, half, solve_kind, l_concurrent, nptcls, result, &
                &history, reduce_time, finalize_time, solve_time, data_scale, lambda_eff, prior_npositive, &
                &prior_positive_min, prior_positive_max, prior_to_khat_l1, prior_to_khat_rms, pcgop )
            integer,                  intent(in) :: state_here, nptcls
            character(len=*),         intent(in) :: half, solve_kind
            logical,                  intent(in) :: l_concurrent
            type(pcg_solver_outcome), intent(in) :: result
            real,                     intent(in) :: history(:)
            real(dp),                 intent(in) :: reduce_time, finalize_time, solve_time
            real,                     intent(in) :: data_scale, lambda_eff
            integer, optional,        intent(in) :: prior_npositive
            real, optional,           intent(in) :: prior_positive_min, prior_positive_max
            real, optional,           intent(in) :: prior_to_khat_l1, prior_to_khat_rms
            type(reconstructor_pcg),   intent(in) :: pcgop
            type(string) :: fname
            integer :: funit, i
            fname = 'reconstruct3D_pcg_state'//int2str_pad(state_here,2)//'_'//trim(half)//'_'// &
                &trim(solve_kind)//'_diagnostics.txt'
            call fopen(funit, file=fname, status='replace', action='write')
            write(funit,'(A,A)')      'execution_mode=',        'distributed'
            write(funit,'(A,A)')      'solve_kind=',             trim(solve_kind)
            write(funit,'(A,L1)')     'half_pair_parallel=',     l_concurrent
            write(funit,'(A,I0)')     'threads_per_half=',       merge(pcg_half_nthreads, pcg_master_nthreads, l_concurrent)
            write(funit,'(A,I0)')     'nparts=',                params%nparts
            write(funit,'(A,I0)')     'nptcls=',                nptcls
            write(funit,'(A,I0)')     'requested_maxits=',      result%requested_maxits
            write(funit,'(A,ES14.6)') 'requested_rtol=',        merge(params%rtol_atm, params%rtol, trim(solve_kind) == 'atm')
            write(funit,'(A,I0)')     'iteration_count=',       result%iteration_count
            write(funit,'(A,A)')      'stop_reason=',           trim(result%stop_reason)
            write(funit,'(A,L1)')     'converged=',             result%converged
            write(funit,'(A,L1)')     'cold_restart_used=',     result%cold_restart_used
            if( result%cold_restart_used )then
                write(funit,'(A,I0)')     'restart_trigger_iteration=', result%restart_trigger_iteration
                write(funit,'(A,ES14.6)') 'restart_trigger_curvature=', result%restart_trigger_curvature
            endif
            write(funit,'(A,ES14.6)') 'initial_rel_resid_l2=',  result%initial_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_resid_l2=',    result%final_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_resid_m=',     result%final_rel_residual_m
            write(funit,'(A,ES14.6)') 'final_rel_update=',      result%final_rel_update
            call write_closed_form_diagnostics(funit, result)
            write(funit,'(A,ES14.6)') 'pcg_data_scale=',        data_scale
            write(funit,'(A,ES14.6)') 'pcg_lambda_effective=',  lambda_eff
            if( present(prior_npositive) )then
                if( .not. present(prior_positive_min) .or. .not. present(prior_positive_max) .or. &
                    &.not. present(prior_to_khat_l1) .or. .not. present(prior_to_khat_rms) )then
                    THROW_HARD('incomplete distributed PCG ML prior diagnostics')
                endif
                write(funit,'(A,I0)')     'ml_prior_nonzero_bins=',       prior_npositive
                write(funit,'(A,ES14.6)') 'ml_prior_positive_min=',       prior_positive_min
                write(funit,'(A,ES14.6)') 'ml_prior_positive_max=',       prior_positive_max
                write(funit,'(A,ES14.6)') 'ml_prior_to_data_khat_l1=',    prior_to_khat_l1
                write(funit,'(A,ES14.6)') 'ml_prior_to_data_khat_rms=',   prior_to_khat_rms
            endif
            write(funit,'(A,F12.6)')  'raw_reduce_seconds=',    reduce_time
            write(funit,'(A,F12.6)')  'master_finalize_seconds=', finalize_time
            write(funit,'(A,F12.6)')  'solve_seconds=',         solve_time
            do i = 1, size(history)
                write(funit,'(A,I0,A,ES14.6)') 'iter', i, '_rel_resid_l2=', history(i)
                write(funit,'(A,I0,A,ES14.6)') 'iter', i, '_rel_update=', result%rel_update_history(i)
                if( result%preconditioned_residual_history(i) >= 0.0 )then
                    write(funit,'(A,I0,A,ES14.6)') 'iter', i, '_preconditioned_resid=', &
                        &result%preconditioned_residual_history(i)
                endif
                write(funit,'(A,I0,A,F12.6)') 'iter', i, '_seconds=', result%iteration_seconds(i)
            enddo
            call pcgop%report_finalize_profile(funit)
            call pcgop%report_profile(result%iteration_count, funit)
            call fclose(funit)
            call fname%kill
        end subroutine write_distributed_diagnostics

    end subroutine execute_rec3D_pcg_distributed_master

    !> The current matching band in native crop-box shells (params%kfromto(2)),
    !! or 0 when it does not describe a usable band of this volume.
    integer function matched_band_kstop( params, volume ) result( kstop )
        type(parameters), intent(in) :: params
        type(image),      intent(in) :: volume
        kstop = params%kfromto(2)
        if( kstop < 1 .or. kstop >= volume%get_filtsz() ) kstop = 0
    end function matched_band_kstop

    !> Diagnostic only: report when the shells beyond the matching band carry
    !! more RMS amplitude than the band-edge shell, the regression signal of
    !! the fixed-iteration solve's beyond-band excess (reconstruct3D_pcg_policy.md,
    !! 'Beyond-band diagnostic'). Silent without a matching band.
    subroutine report_beyond_band_excess( volume, params, state, half, solve_kind )
        type(image),      intent(in) :: volume
        type(parameters), intent(in) :: params
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half, solve_kind
        real, parameter :: EXCESS_REPORT_RATIO = 10.
        type(image) :: tmpvol
        complex  :: comp
        real(dp) :: sumsq_edge, sumsq_beyond
        real     :: ratio
        integer  :: lims(3,2), phys(3), h, k, l, sh, n_edge, n_beyond, kstop
        kstop = matched_band_kstop(params, volume)
        if( kstop < 1 ) return
        call tmpvol%copy(volume)
        call tmpvol%fft()
        lims         = tmpvol%loop_lims(2)
        sumsq_edge   = 0.0_dp
        sumsq_beyond = 0.0_dp
        n_edge       = 0
        n_beyond     = 0
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,comp) &
        !$omp reduction(+:sumsq_edge,sumsq_beyond,n_edge,n_beyond) schedule(static) proc_bind(close)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do l = lims(3,1),lims(3,2)
                    sh = nint(sqrt(real(h*h + k*k + l*l)))
                    if( sh < kstop ) cycle
                    phys = tmpvol%comp_addr_phys(h,k,l)
                    comp = tmpvol%get_cmat_at(phys(1),phys(2),phys(3))
                    if( sh == kstop )then
                        sumsq_edge = sumsq_edge + real(comp,dp)**2 + real(aimag(comp),dp)**2
                        n_edge     = n_edge + 1
                    else
                        sumsq_beyond = sumsq_beyond + real(comp,dp)**2 + real(aimag(comp),dp)**2
                        n_beyond     = n_beyond + 1
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        call tmpvol%kill
        if( n_edge < 1 .or. n_beyond < 1 .or. sumsq_edge <= 0.0_dp ) return
        ratio = real(sqrt( (sumsq_beyond / real(n_beyond,dp)) / (sumsq_edge / real(n_edge,dp)) ))
        if( ratio >= EXCESS_REPORT_RATIO )then
            write(logfhandle,'(A,I0,A,A,A,A,A,I0,A,ES9.2)') '>>> PCG BEYOND-BAND EXCESS: STATE ', state, &
                &' | HALF=', trim(half), ' | KIND=', trim(solve_kind), ' | BAND EDGE k=', kstop, &
                &' | BEYOND/EDGE RMS RATIO=', ratio
        endif
    end subroutine report_beyond_band_excess

    subroutine write_output_diagnostics( state, execution_mode, map_time, fsc_time, nu_filter_time )
        integer,            intent(in) :: state
        character(len=*),   intent(in) :: execution_mode
        real(dp),           intent(in) :: map_time, fsc_time
        real(dp), optional, intent(in) :: nu_filter_time
        type(string) :: fname
        integer :: funit
        fname = 'reconstruct3D_pcg_state'//int2str_pad(state,2)//'_output_diagnostics.txt'
        call fopen(funit, file=fname, status='replace', action='write')
        write(funit,'(A,A)')     'execution_mode=', trim(execution_mode)
        write(funit,'(A,F12.6)') 'halfmap_merged_output_seconds=', map_time
        write(funit,'(A,F12.6)') 'fsc_cfar_summary_seconds=', fsc_time
        if( present(nu_filter_time) ) write(funit,'(A,F12.6)') 'nu_filter_seconds=', nu_filter_time
        call fclose(funit)
        call fname%kill
    end subroutine write_output_diagnostics

    !> one record per distributed PCG master assembly in BENCH/PCG_ASSEMBLE_BENCH.txt: the phases are
    !! wall-clock and disjoint, the split timings are parts of them
    subroutine write_pcg_assemble_bench( params, bench )
        use simple_bench_io, only: open_bench_report
        type(parameters),         intent(in) :: params
        type(pcg_assemble_bench), intent(in) :: bench
        real(dp) :: rt_accounted
        integer  :: fnr
        if( .not. L_BENCH_GLOB ) return
        rt_accounted = 0._dp
        if( bench%rt_tot > TINY ) rt_accounted = 100._dp * (bench%rt_init + bench%rt_support + bench%rt_op_setup + &
            &bench%rt_read + bench%rt_sum + bench%rt_chain + sum(bench%rt_finalize) + sum(bench%rt_solve) + &
            &bench%rt_prior_resolve + bench%rt_finish + bench%rt_fsc + bench%rt_map_output + bench%rt_nu_filter + &
            &bench%rt_proj_write + bench%rt_cleanup + bench%rt_nu_handoff) / bench%rt_tot
        call open_bench_report(fnr, 'PCG_ASSEMBLE', iter=params%which_iter)
        write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
        write(fnr,'(a)')    'pcg_assemble execution mode             : distributed'
        write(fnr,'(a,i0)') 'pcg_assemble box_crop                   : ', params%box_crop
        write(fnr,'(a,i0)') 'pcg_assemble nstates                    : ', params%nstates
        write(fnr,'(a,i0)') 'pcg_assemble states solved              : ', bench%nstates_solved
        write(fnr,'(a,i0)') 'pcg_assemble nparts                     : ', params%nparts
        write(fnr,'(a,i0)') 'pcg_assemble master threads             : ', bench%nthr_master
        write(fnr,'(a,i0)') 'pcg_assemble threads per half           : ', bench%nthr_half
        write(fnr,'(a,a)')  'pcg_assemble ml_reg                     : ', trim(merge('yes','no ',params%l_ml_reg))
        write(fnr,'(a,a)')  'pcg_assemble trail_rec                  : ', trim(merge('yes','no ',params%l_trail_rec))
        write(fnr,'(a,a)')  'pcg_assemble nonuniform                 : ', trim(merge('yes','no ',params%l_nonuniform))
        write(fnr,'(a,i0)') 'pcg_assemble maxits_pcg                 : ', params%maxits_pcg
        write(fnr,'(a,i0)') 'pcg_assemble maxits_ml                  : ', params%maxits_ml
        write(fnr,'(a,i0)') 'pcg_assemble raw files read             : ', bench%nraw
        write(fnr,'(a,i0)') 'pcg_assemble raw bytes read             : ', bench%raw_bytes
        write(fnr,'(a,i0)') 'pcg_assemble base solve iters even      : ', bench%niters(1,1)
        write(fnr,'(a,i0)') 'pcg_assemble base solve iters odd       : ', bench%niters(2,1)
        write(fnr,'(a,i0)') 'pcg_assemble ml solve iters even        : ', bench%niters(1,2)
        write(fnr,'(a,i0)') 'pcg_assemble ml solve iters odd         : ', bench%niters(2,2)
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** TIMINGS (s) ***'
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble init                       :', bench%rt_init
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble solve support              :', bench%rt_support
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble operator setup             :', bench%rt_op_setup
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble raw read                   :', bench%rt_read
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble raw summation              :', bench%rt_sum
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble trailing chain             :', bench%rt_chain
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble operator finalize          :', sum(bench%rt_finalize)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble base solve                 :', bench%rt_solve(1)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble prior re-solve             :', bench%rt_prior_resolve
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble ml solve                   :', bench%rt_solve(2)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble solve finish               :', bench%rt_finish
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble fsc and cfar               :', bench%rt_fsc
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble map output                 :', bench%rt_map_output
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nonuniform filter          :', bench%rt_nu_filter
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble project write              :', bench%rt_proj_write
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble cleanup                    :', bench%rt_cleanup
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu lp handoff              :', bench%rt_nu_handoff
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble total time                 :', bench%rt_tot
        write(fnr,'(a,1x,f0.2)') 'pcg_assemble % accounted for            :', rt_accounted
        write(fnr,'(a)') ''
        write(fnr,'(a)') '*** SPLIT TIMINGS (S) ***'
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble base operator finalize     :', bench%rt_finalize(1)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble ml operator finalize       :', bench%rt_finalize(2)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble base solve even half       :', bench%rt_half_solve(1,1)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble base solve odd half        :', bench%rt_half_solve(2,1)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble ml solve even half         :', bench%rt_half_solve(1,2)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble ml solve odd half          :', bench%rt_half_solve(2,2)
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu filter envelope         :', bench%nu%envmask
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu filter candidate filter :', bench%nu%filtering
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu filter objective        :', bench%nu%objective
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu filter optimisation     :', bench%nu%optimise
        write(fnr,'(a,1x,f0.3)') 'pcg_assemble nu filter composition      :', bench%nu%compose
        call fclose(fnr)
    end subroutine write_pcg_assemble_bench

    logical function pcg_trail_seed_requested( cline ) result(l_seed)
        class(cmdline), intent(in) :: cline
        type(string) :: value
        l_seed = .false.
        if( .not. cline%defined('trail_seed') ) return
        value = cline%get_carg('trail_seed')
        l_seed = trim(value%to_char()) == 'yes'
        call value%kill
    end function pcg_trail_seed_requested

    subroutine validate_pcg_common( params, check_solver )
        type(parameters), intent(in) :: params
        logical, optional, intent(in) :: check_solver
        logical :: l_check_solver
        l_check_solver = .true.
        if( present(check_solver) ) l_check_solver = check_solver
        if( l_check_solver )then
            if( params%maxits_pcg < 1 .or. params%maxits_pcg > 100 ) &
                &THROW_HARD('PCG requires 1<=maxits_pcg<=100')
            if( params%maxits_ml < 0 .or. params%maxits_ml > 100 ) &
                &THROW_HARD('PCG requires 0<=maxits_ml<=100 (0: closed-form regularized pair only)')
            if( params%maxits_pcg > 8 )then
                THROW_WARN('maxits_pcg exceeds the production refinement budget (8); appropriate for offline converged solves only')
            endif
        endif
        if( abs(real(params%box)*params%smpd - real(params%box_crop)*params%smpd_crop) > &
            &1.0e-5*real(params%box)*params%smpd )then
            THROW_HARD('PCG crop must preserve the native physical box extent')
        endif
        if( params%msk <= 0.5 .or. params%msk_crop <= 0.5 ) THROW_HARD('rec_backend=pcg requires mskdiam')
        if( .not. ieee_is_finite(params%rtol) ) THROW_HARD('PCG rtol must be finite')
    end subroutine validate_pcg_common

    function pcg_raw_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgraw-v2|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|iter='//trim(int2str(params%which_iter))// &
            &'|box='//trim(int2str(params%box))//'|smpd='//trim(real2str(params%smpd))// &
            &'|box_crop='//trim(int2str(params%box_crop))// &
            &'|smpd_crop='//trim(real2str(params%smpd_crop))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_raw_provenance

    !> Chain identity: native geometry and objective; neither
    !! the iteration nor the crop, so the chain survives stage transitions.
    !! The represented population M(s) of the population rule is in the pair's
    !! manifest; an older chain is discarded and re-seeded
    function pcg_chain_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgtrail-v3|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_chain_provenance

end module simple_rec3D_pcg_strategy
