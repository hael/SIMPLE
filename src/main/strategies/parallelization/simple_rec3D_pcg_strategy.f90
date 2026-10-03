!@descr: kernel PCG reconstruct3D routes: shared-memory solve, distributed worker accumulation and master reduction and solve
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,             only: builder
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_reconstructor_pcg,   only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_LAMBDA, &
    &pcg_raw_accum_compatible, measure_closed_form_agreement, handle_cold_restart_outcome, report_pcg_solve, &
    &report_closed_form_agreement, write_closed_form_diagnostics, validate_solved_map, read_pcg_raw_accum_header
use simple_matcher_ptcl_io,     only: prepimgbatch, discrete_read_imgbatch, killimgbatch, prep_rec_observation
use simple_sigma2_files,        only: load_sigma2_groups
use simple_math_ft,             only: resample_sigma2
use simple_image,               only: image
use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, evaluate_halfmap_pair, write_halfmap_diagnostics, &
    &write_support_provenance, read_support_provenance
use simple_image_msk,           only: image_msk
use simple_pcg_solvent_sidecar, only: solvent_prior_cross_half_objective, PCG_SOLVENT_LAMBDA_GRID, &
    &write_pcg_solvent_pair, prepare_solvent_prior_on_pair, solvent_prior_provenance, build_solvent_check_support
use simple_nu_state_filter,     only: nonuniform_filter_state, nu_aux_member
use simple_refine3D_fnames,     only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname, &
    &refine3D_pcg_trail_accum_fname
use simple_frozen_accum,        only: frozen_accum
use simple_oris,                only: population_blend_weights
!$ use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_get_max_active_levels, &
!$     &omp_set_max_active_levels, omp_set_num_threads
implicit none

public :: execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
public :: rec3D_master_nthr, validate_pcg_common
private
#include "simple_local_flags.inc"

integer, parameter :: PCG_MASTER_NTHR_CAP    = 32   !< master-phase thread-boost ceiling
! Solve support follows automsk (reconstruct3D_pcg_policy.md, section 3):
! the sphere for automsk=no, the lag-one density envelope for yes and nu
! (never the NU-evidence envelope), an explicit pcg_mskfile over both.
logical, parameter :: DEBUG = .false.

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
                call pcgop%set_mask_volume(mskvol)
                call mskvol%kill
                return
            endif
        endif
        l_state = .false.
        if( present(l_state_support) ) l_state = l_state_support
        if( l_state .and. present(state_support) )then
            call pcgop%set_mask_volume(state_support)
            return
        endif
        call pcgop%set_mask(params%msk_crop)
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
    subroutine execute_rec3D_pcg_worker( params, build, cline, selected_pinds )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: selected_pinds(:)
        integer, allocatable :: half_pinds(:)
        character(len=256) :: provenance
        integer :: state, eo, n_half
        logical :: l_sigma_loaded

        call validate_pcg_common(params, check_solver=.false.)
        provenance = pcg_raw_provenance(params)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            l_sigma_loaded = allocated(build%esig%sigma2_noise)
            if( .not. l_sigma_loaded )then
                call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, &
                    &build%spproj_field, l_sigma_loaded)
            endif
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif
        if( size(selected_pinds) > 0 )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
        endif
        do state = 1, params%nstates
            do eo = 0, 1
                call collect_worker_state_half(state, eo, selected_pinds, half_pinds)
                n_half = size(half_pinds)
                call accumulate_worker_state_half(state, eo, half_pinds, provenance)
                deallocate(half_pinds)
                if( DEBUG )then
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> PCG RAW WORKER: PART ', params%part, &
                        &' STATE ', state, ' HALF ', eo, ' PARTICLES ', n_half
                endif
            enddo
        enddo
        if( size(selected_pinds) > 0 ) call killimgbatch(build)

    contains

        subroutine collect_worker_state_half( state_here, eo_here, pinds, selected )
            integer,              intent(in)  :: state_here, eo_here, pinds(:)
            integer, allocatable, intent(out) :: selected(:)
            integer :: i, n, p
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                n = n + 1
            enddo
            allocate(selected(n))
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                n = n + 1
                selected(n) = p
            enddo
        end subroutine collect_worker_state_half

        subroutine accumulate_worker_state_half( state_here, eo_here, pinds, provenance_here )
            integer,          intent(in) :: state_here, eo_here, pinds(:)
            character(len=*), intent(in) :: provenance_here
            type(reconstructor_pcg) :: pcgop
            type(oris)      :: selection
            type(ori)       :: orientation
            type(ctfparams) :: ctfparms
            type(string)    :: fname
            type(image)     :: obs
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch
            real    :: shift(2), crop_factor

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            fname = refine3D_pcg_raw_accum_fname(state_here, params%part, params%numlen, &
                &merge('odd ', 'even', eo_here == 1))
            if( size(pinds) == 0 )then
                call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                    &params%nparts, 0, provenance_here)
                call pcgop%kill
                call fname%kill
                return
            endif
            call pcgop%set_sym(build%pgrpsyms)
            lims2 = pcgop%get_lims2()
            R     = lims2(1,2)
            allocate(sig2(0:R,size(pinds)), source=1.0)
            if( params%cc_objfun == OBJFUN_EUCLID )then
                kfromto = build%esig%get_kfromto()
                do i = 1, size(pinds)
                    call resample_sigma2(kfromto(1), kfromto(2), &
                        &build%esig%sigma2_noise(kfromto(1):kfromto(2),pinds(i)), R, 1.0, sig2(0:R,i))
                enddo
            endif
            call selection%new(size(pinds), .true.)
            call orientation%new(.false.)
            crop_factor = real(params%box_crop) / real(params%box)
            do i = 1, size(pinds)
                iptcl = pinds(i)
                call build%spproj_field%get_ori(iptcl, orientation)
                ctfparms      = build%spproj%get_ctfparams(params%oritype, iptcl)
                ctfparms%smpd = params%smpd_crop
                shift         = build%spproj_field%get_2Dshift(iptcl) * crop_factor
                call orientation%set_ctfvars(ctfparms)
                call orientation%set_shift(shift)
                call selection%set_ori(i, orientation)
            enddo
            call pcgop%prep_particles(selection, use_ctf=.true., sig2=sig2)
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call pcgop%begin_accum
            call obs%new([params%box_crop,params%box_crop,1], params%smpd_crop)
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                do ii = 1, batchsz
                    ! the backend-neutral observation (normalize, crop, taper), see prep_rec_observation
                    call prep_rec_observation(build%imgbatch(ii), build%lmsk, obs, .true.)
                    call obs%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(obs)
                enddo
                call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
            enddo
            call obs%kill
            call pcgop%write_raw_accum(fname, state_here, eo_here, params%part, &
                &params%nparts, size(pinds), provenance_here)
            call pcgop%kill
            call selection%kill
            call orientation%kill
            call fname%kill
            deallocate(y_batch, sig2)
        end subroutine accumulate_worker_state_half

    end subroutine execute_rec3D_pcg_worker

    !> Distributed master: reduce raw worker B,D artifacts in ascending part
    !! order, then perform all folding, finalization and PCG locally. For each
    !! even/odd pair, construction and teardown stay serial while the two fully
    !! prepared PCG solves execute concurrently with disjoint thread budgets.
    subroutine execute_rec3D_pcg_distributed_master( params, build, cline, trail_bootstrap_states, &
            &nu_align_lps )
        type :: distributed_half_job
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            real, allocatable :: x(:,:,:), x_cf(:,:,:), rel_res_hist(:)
            integer :: state = 0, eo = 0, nptcls = 0, niters = 0
            integer :: nfrozen = 0     !< frozen particles added to the reduction (solve3D_addon)
            integer :: band_shell = 0 !< the pair's FSC=0.143 shell (regularized solve; agreement diagnostic)
            integer :: prior_npositive = 0
            character(len=8) :: half = '', solve_kind = ''
            real(dp) :: time_reduce = 0.0_dp, time_finalize = 0.0_dp, time_solve = 0.0_dp
            real :: prior_positive_min = 0.0, prior_positive_max = 0.0
            real :: prior_to_khat_l1 = 0.0, prior_to_khat_rms = 0.0
            logical :: l_ml_solve = .false., ready = .false., l_concurrent = .false.
            logical :: l_nonzero = .false. !< nonzero initial guess: eligible for the cold restart
        end type distributed_half_job
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        logical, optional, intent(out)  :: trail_bootstrap_states(:)
        real,    optional, intent(out)  :: nu_align_lps(:)
        type(image), target  :: half_even, half_odd, ml_even, ml_odd, merged, solvent_even, solvent_odd
        type(image), target  :: previous_even, previous_odd, previous_merged
        type(image), pointer :: fsc_pair_even, fsc_pair_odd, fsc_pair_merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt, eonames(2)
        type(halfmap_diagnostics_result) :: hm_diag
        real, allocatable :: fsc(:), res0143s(:), res05s(:), cfars(:), align_lps(:)
        real, allocatable :: realized_fractions(:), update_weights(:), chain_weights(:), current_scales(:)
        integer, allocatable :: nrep(:), nsmp(:)
        logical, allocatable :: state_written(:)
        character(len=256) :: provenance, chain_provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        integer :: n_active_state, n_sampled_state
        integer :: pcg_master_nthreads, pcg_half_nthreads
        type(distributed_half_job) :: even_job, odd_job
        type(image_msk) :: state_support_msk
        type(image)     :: solvent_weight(2)
        logical :: l_state_support, l_base_support_constrained, l_nu_base_constrained
        logical :: l_has_updates, l_bootstrap, l_even_chain, l_odd_chain
        logical :: l_fsc_pair_support_constrained, l_prev_support_constrained, l_prev_provenance_found
        logical :: l_shipped_support_constrained, l_solvent_weight
        real    :: res0143_prior_free
        real    :: solvent_lambda_eff !< the solvent ridge in force for the current state
        character(len=16) :: state_support_kind, base_support_kind, fsc_support_kind
        character(len=16) :: previous_support_kind, shipped_support_kind
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output, time_nu_filter
        type(frozen_accum) :: frozen_ctx
        logical :: l_frozen_rec, l_frozen_seed
        integer :: nfz_even, nfz_odd

        call validate_pcg_common(params)
        solvent_lambda_eff = 0.
        ! solve3D_addon handshakes (in-process only): the master owns every
        ! frozen read and write; workers only ever see their own particles
        l_frozen_rec  = cline%defined('frozen_rec')
        l_frozen_seed = cline%defined('frozen_seed')
        if( l_frozen_rec .and. l_frozen_seed ) THROW_HARD('a reconstruction cannot both produce and consume a frozen set')
        if( l_frozen_rec ) call frozen_ctx%load(cline%get_carg('frozen_rec'), 'pcg', &
            &params%nstates, build%spproj_field%get_noris(), params%cc_objfun, producer=.false.)
        if( l_frozen_seed ) call frozen_ctx%load(cline%get_carg('frozen_seed'), 'pcg', &
            &params%nstates, build%spproj_field%get_noris(), params%cc_objfun, producer=.true.)
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
        if( present(trail_bootstrap_states) )then
            if( size(trail_bootstrap_states) /= params%nstates ) &
                &THROW_HARD('PCG trailing-bootstrap state output has invalid size')
            trail_bootstrap_states = .false.
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
        if( params%l_trail_rec )then
            ! N(s) and n(s) of the population rule; f = n/N
            call build%spproj%os_ptcl3D%get_group_update_counts('state', params%nstates, nrep, nsmp)
            call build%spproj%os_ptcl3D%get_state_update_fracs(params%nstates, realized_fractions)
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
        do state = 1, params%nstates
            l_bootstrap = .false.
            if( params%l_trail_rec )then
                ! a chain pair of another identity is discarded and re-seeded
                call discard_stale_trail_chain_pair(state)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'even')
                l_even_chain = file_exists(raw_fname)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'odd')
                l_odd_chain = file_exists(raw_fname)
                if( l_even_chain .neqv. l_odd_chain ) THROW_HARD('PCG trailing chain pair is incomplete')
                l_bootstrap = .not. l_even_chain
                ! add-on mode never enters the legacy union-volume bootstrap
                if( l_bootstrap .and. l_frozen_rec ) &
                    &THROW_HARD('solve3D_addon trailing assembly requires a seeded cohort chain')
                if( .not. l_bootstrap ) call set_chain_blend_weights(state)
            endif
            if( present(trail_bootstrap_states) ) trail_bootstrap_states(state) = l_bootstrap
            ! one solve support per state (build_pcg_state_support)
            call build_pcg_state_support(params, state, state_support_msk, l_state_support, state_support_kind)
            l_solvent_weight = .false.
            l_base_support_constrained = l_state_support
            base_support_kind = 'sphere'
            if( l_base_support_constrained ) base_support_kind = state_support_kind
            call reduce_solve_state_pair(state, half_even, half_odd, n_even, n_odd, 'base', &
                &solvent_even=solvent_even, solvent_odd=solvent_odd)
            nfz_even = even_job%nfrozen
            nfz_odd  = odd_job%nfrozen
            if( params%l_trail_rec )then
                call count_state_sampling(state, n_active_state, n_sampled_state)
                if( n_even+n_odd /= n_sampled_state ) THROW_HARD('PCG raw particles do not match the latest sampled cohort')
                if( n_active_state > 0 )then
                    if( abs(real(n_sampled_state)/real(n_active_state)-realized_fractions(state)) > 1.0e-6 )then
                        THROW_HARD('PCG realized fraction disagrees with gridding sampling bookkeeping')
                    endif
                endif
            endif
            ! n_even/n_odd count this run's own particles; an add-on half also
            ! carries its frozen particles
            if( n_even + nfz_even == 0 .and. n_odd + nfz_odd == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                    &' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even + nfz_even < 1 .or. n_odd + nfz_odd < 1 ) THROW_HARD('distributed PCG requires both halfsets')
            fname_even = refine3D_state_halfvol_fname(state, 'even')
            fname_odd  = refine3D_state_halfvol_fname(state, 'odd')
            fname_vol  = refine3D_state_vol_fname(state)
            fname_fsc  = refine3D_fsc_fname(state)
            call merged%copy(half_even)
            call merged%add(half_odd)
            call merged%mul(0.5)
            if( params%l_ml_reg .and. .not. l_state_support .and. pcg_density_support_enabled(params) )then
                call build_pcg_density_support(params, state, merged, state_support_msk, 'current base pair')
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
            if( l_solvent_weight ) call write_pcg_solvent_pair(params, state, solvent_even, solvent_odd)
            t_state_phase = tic()
            ! the FSC pair: the current base pair, or the previous shipped pair in
            ! the trailing bootstrap (lag-one, as gridding); the NU bank is always
            ! seeded from the current base pair
            if( l_bootstrap )then
                call load_previous_state_halves(state, previous_even, previous_odd, previous_merged)
                fsc_pair_even        => previous_even
                fsc_pair_odd         => previous_odd
                fsc_pair_merged      => previous_merged
                ! the previous pair's own solve support, from its provenance
                ! sidecar (none: unconstrained)
                previous_support_kind = 'sphere'
                call read_support_provenance(params%vols(state), l_prev_support_constrained, &
                    &l_prev_provenance_found, support_kind=previous_support_kind)
                l_fsc_pair_support_constrained = l_prev_provenance_found .and. l_prev_support_constrained
                fsc_support_kind = 'sphere'
                if( l_fsc_pair_support_constrained ) fsc_support_kind = previous_support_kind
                if( .not. l_prev_provenance_found ) write(logfhandle,'(A,I0,A)') &
                    &'>>> PCG DISTRIBUTED: STATE ', state, &
                    &' previous pair has no solve-support provenance; treated as unconstrained'
            else
                fsc_pair_even        => half_even
                fsc_pair_odd         => half_odd
                fsc_pair_merged      => merged
                l_fsc_pair_support_constrained = l_base_support_constrained
                fsc_support_kind = base_support_kind
            endif
            call calculate_pcg_state_diagnostics(params, state, 'DISTRIBUTED', fsc_pair_even, &
                &fsc_pair_odd, fsc_pair_merged, hm_diag, l_fsc_pair_support_constrained, fsc_support_kind)
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
            if( l_bootstrap .and. update_weights(state) < 0.99 )then
                ! bootstrap blend, as gridding's trail_restored_halves_if_needed: the
                ! previous pair scaled once and added to the base and the ML pair
                call previous_even%mul(1.0-update_weights(state))
                call previous_odd%mul( 1.0-update_weights(state))
                call blend_bootstrap_half(half_even, previous_even, update_weights(state))
                call blend_bootstrap_half(half_odd,  previous_odd,  update_weights(state))
                if( params%l_ml_reg )then
                    call blend_bootstrap_half(ml_even, previous_even, update_weights(state))
                    call blend_bootstrap_half(ml_odd,  previous_odd,  update_weights(state))
                endif
                if( l_solvent_weight )then
                    call blend_bootstrap_half(solvent_even, previous_even, update_weights(state))
                    call blend_bootstrap_half(solvent_odd,  previous_odd,  update_weights(state))
                endif
                if( params%l_lpset )then
                    call merged%kill
                    if( params%l_ml_reg )then
                        call merged%copy(ml_even)
                        call merged%add(ml_odd)
                    else
                        call merged%copy(half_even)
                        call merged%add(half_odd)
                    endif
                    call merged%mul(0.5)
                endif
            endif
            t_state_phase = tic()
            l_shipped_support_constrained = merge(l_state_support, l_base_support_constrained, params%l_ml_reg)
            shipped_support_kind = base_support_kind
            if( params%l_ml_reg ) shipped_support_kind = state_support_kind
            ! a bootstrap blend carries the previous pair's support into the
            ! shipped pair: constrained only if both contributions were
            if( l_bootstrap .and. update_weights(state) < 0.99 )then
                l_shipped_support_constrained = l_shipped_support_constrained .and. l_fsc_pair_support_constrained
                if( l_shipped_support_constrained )then
                    if( trim(shipped_support_kind) /= trim(fsc_support_kind) ) shipped_support_kind = 'mixed'
                else
                    shipped_support_kind = 'sphere'
                endif
            endif
            if( params%l_ml_reg )then
                call ml_even%write(fname_even, del_if_exists=.true.)
                call ml_odd%write(fname_odd, del_if_exists=.true.)
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
            else if( l_bootstrap .and. update_weights(state) < 0.99 .and. l_solvent_weight )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'mixed', shipped_support_kind, &
                    &solvent_prior=solvent_prior_provenance(params, solvent_lambda_eff))
            else if( l_bootstrap .and. update_weights(state) < 0.99 )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'mixed', shipped_support_kind)
            else if( l_solvent_weight )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'base', shipped_support_kind, &
                    &solvent_prior=solvent_prior_provenance(params, solvent_lambda_eff))
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
                ! on the dilation ring); a bootstrap blend only if both parts were
                l_nu_base_constrained = l_base_support_constrained
                if( l_bootstrap .and. update_weights(state) < 0.99 ) &
                    &l_nu_base_constrained = l_nu_base_constrained .and. l_fsc_pair_support_constrained
                if( l_solvent_weight )then
                    ! label field from the base pair, applied to the prior'd pair
                    call nonuniform_filter_state(params, state, half_even, half_odd, &
                        &ml_even, ml_odd, nu_aux_member(params), &
                        &res0143s(state), fname_vol, eonames, align_lps(state), &
                        &base_support=state_support_msk, l_base_constrained=l_nu_base_constrained, &
                        &vol_apply_even=solvent_even, vol_apply_odd=solvent_odd)
                else
                    call nonuniform_filter_state(params, state, half_even, half_odd, &
                        &ml_even, ml_odd, nu_aux_member(params), &
                        &res0143s(state), fname_vol, eonames, align_lps(state), &
                        &base_support=state_support_msk, l_base_constrained=l_nu_base_constrained)
                endif
                time_nu_filter = real(toc(t_state_phase),dp)
            endif
            call write_output_diagnostics(state, 'distributed', time_map_output, time_fsc_output, time_nu_filter)
            params%vols(state)      = fname_vol
            params%vols_even(state) = fname_even
            params%vols_odd(state)  = fname_odd
            call cline%set('vol'//int2str(state), fname_vol)
            state_written(state) = .true.
            call half_even%kill
            call half_odd%kill
            call solvent_even%kill
            call solvent_odd%kill
            if( params%l_ml_reg )then
                call ml_even%kill
                call ml_odd%kill
                call fname_even_unfil%kill
                call fname_odd_unfil%kill
            endif
            if( l_bootstrap )then
                call previous_even%kill
                call previous_odd%kill
                call previous_merged%kill
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
        if( .not. any(state_written) ) THROW_HARD('distributed PCG produced no populated states')
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
        call frozen_ctx%kill
        deallocate(res0143s, res05s, cfars, state_written, realized_fractions, update_weights, align_lps)
        deallocate(chain_weights, current_scales)
        if( allocated(nrep) ) deallocate(nrep, nsmp)
        !$ call omp_set_num_threads(params%nthr)

    contains

        !> Population-rule weights of one state's chain pair (class-average and
        !! reconstruct3D partials note, Section 4.1). The chain's represented
        !! population M(s) is the sum of the particle counts in its two headers,
        !! written as the represented population of each half; with N(s) active
        !! updated rows and n(s) in the current sample, current *= s = u/f and
        !! chain *= w = (1-u)*N/M, so the blended mass is N whatever joined or left
        !! the state, and the current-map coefficient stays u. With N = M this is
        !! the former chain weight 1 - u.
        subroutine set_chain_blend_weights( state_here )
            integer, intent(in) :: state_here
            type(string) :: chain_fname
            character(len=256) :: prov_here
            real    :: smpd_here, mnew
            integer :: ieo, st_here, eo_here, part_here, nparts_here, npop, box_here, status, mrep, nnew
            character(len=4), parameter :: HALVES(2) = ['even', 'odd ']
            mrep = 0
            do ieo = 1, 2
                chain_fname = refine3D_pcg_trail_accum_fname(state_here, trim(HALVES(ieo)))
                call read_pcg_raw_accum_header(chain_fname, st_here, eo_here, part_here, nparts_here, npop, &
                    &box_here, smpd_here, prov_here, status)
                if( status /= 0 ) THROW_HARD('unreadable PCG trailing chain header')
                mrep = mrep + npop
                call chain_fname%kill
            enddo
            call population_blend_weights(nrep(state_here), nsmp(state_here), real(mrep), &
                &current_scales(state_here), chain_weights(state_here), mnew, ufrac=update_weights(state_here))
            nnew = count_first_time(state_here)
            write(logfhandle,'(A,I0,A,F8.4,A,F8.4,A,F8.4,A,4I9,A,I9)') '>>> PCG TRAILING BLEND, STATE ', state_here, &
                &', PREVIOUS-CHAIN WEIGHT ', chain_weights(state_here), ', CURRENT SCALE ', current_scales(state_here), &
                &', FORMER WEIGHT 1-U ', 1.0 - update_weights(state_here), ', N n M MNEW', nrep(state_here), &
                &nsmp(state_here), mrep, nint(mnew), ', FIRST-TIME', nnew
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

        subroutine count_state_sampling( state_here, n_active, n_sampled )
            integer, intent(in)  :: state_here
            integer, intent(out) :: n_active, n_sampled
            integer :: p, sample_ind
            n_active  = 0
            n_sampled = 0
            ! Match get_state_update_fracs without reaching into the private
            ! sampling API: the current cohort is the largest sampled index.
            sample_ind = 0
            do p = 1, build%spproj_field%get_noris()
                sample_ind = max(sample_ind, build%spproj_field%get_sampled(p))
            enddo
            do p = 1, build%spproj_field%get_noris()
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_updatecnt(p) < 1 ) cycle
                n_active = n_active + 1
                if( build%spproj_field%get_sampled(p) == sample_ind ) n_sampled = n_sampled + 1
            enddo
        end subroutine count_state_sampling

        !> Discard a chain pair of another identity: provenance, field of view, a
        !! larger crop than the current, or an unreadable or old format; constant-
        !! FOV crop growth is not stale (zero-extension on read). Both halves go
        !! together; the caller re-enters the bootstrap blend.
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
            if( l_stale )then
                if( l_even_here ) call del_file(even_fname)
                if( l_odd_here  ) call del_file(odd_fname)
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: DISCARDING STALE TRAILING CHAIN, STATE ', &
                    &state_here, ' (GEOMETRY/IDENTITY CHANGE); RE-SEEDING VIA BOOTSTRAP'
            endif
            call even_fname%kill
            call odd_fname%kill
        end subroutine discard_stale_trail_chain_pair

        subroutine load_previous_state_halves( state_here, even, odd, avg )
            integer,     intent(in)    :: state_here
            type(image), intent(inout) :: even, odd, avg
            type(string) :: previous_volume, previous_even_fname, previous_odd_fname
            previous_volume     = params%vols(state_here)
            previous_even_fname = add2fbody(previous_volume, MRC_EXT, '_even')
            previous_odd_fname  = add2fbody(previous_volume, MRC_EXT, '_odd')
            if( .not. file_exists(previous_even_fname) ) THROW_HARD('PCG trailing bootstrap requires the previous even halfmap')
            if( .not. file_exists(previous_odd_fname) ) THROW_HARD('PCG trailing bootstrap requires the previous odd halfmap')
            call even%read_and_crop(previous_even_fname, params%smpd, params%box_crop, params%smpd_crop)
            call odd%read_and_crop(previous_odd_fname, params%smpd, params%box_crop, params%smpd_crop)
            call avg%copy(even)
            call avg%add(odd)
            call avg%mul(0.5)
            call previous_volume%kill
            call previous_even_fname%kill
            call previous_odd_fname%kill
        end subroutine load_previous_state_halves

        !> current := weight_current * current + previous, where previous has
        !! already been scaled by (1 - weight_current) exactly once by the
        !! caller so the same previous pair can anchor several current pairs
        subroutine blend_bootstrap_half( current, previous, weight_current )
            type(image), intent(inout) :: current
            type(image), intent(in)    :: previous
            real,        intent(in)    :: weight_current
            call current%mul(weight_current)
            call current%add(previous)
        end subroutine blend_bootstrap_half

        integer function count_full_state_half( state_here, eo_here ) result(n)
            integer, intent(in) :: state_here, eo_here
            integer :: p
            n = 0
            do p = params%fromp, params%top
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                if( l_has_updates .and. build%spproj_field%get_updatecnt(p) < 1 ) cycle
                n = n + 1
            enddo
        end function count_full_state_half

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

        !> pcg_solvent=yes: solvent_even/odd receive the ridge re-solve, even/odd
        !! the prior-free base pair (FSC, NU, evidence, _unfil)
        subroutine reduce_solve_state_pair( state_here, even, odd, n_even_here, n_odd_here, solve_kind, &
                &fsc_prior, warm_even, warm_odd, solvent_even, solvent_odd )
            integer,          intent(in)    :: state_here
            character(len=*), intent(in)    :: solve_kind
            type(image),      intent(inout) :: even, odd
            integer,          intent(out)   :: n_even_here, n_odd_here
            real, optional,   intent(in)    :: fsc_prior(:)
            type(image), optional, intent(in) :: warm_even, warm_odd
            type(image), optional, intent(inout) :: solvent_even, solvent_odd
            logical :: l_resolved
            l_resolved = .false.

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

            call solve_distributed_half_pair(even_job, odd_job)

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
            integer(timer_int_kind) :: t_phase
            real :: realized_fraction, update_weight, current_scale
            logical :: l_chain_exists, l_seed_chain

            job%state = state_here
            job%eo = eo_here
            job%half = half
            job%solve_kind = solve_kind
            job%nptcls = 0
            job%nfrozen = 0
            job%niters = 0
            job%ready = .false.
            job%l_ml_solve = present(fsc_prior)
            if( job%l_ml_solve .neqv. present(warm_start) ) &
                &THROW_HARD('distributed PCG regularized pair requires both FSC and base map')
            if( allocated(job%x) ) deallocate(job%x)
            if( allocated(job%rel_res_hist) ) deallocate(job%rel_res_hist)

            call job%pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA, &
                &fft_nthreads=pcg_half_nthreads)
            ! serial: the spherical mask memoizes coordinates at module scope
            call set_pcg_solve_support(job%pcgop, params, state_support_msk, l_state_support)
            call job%pcgop%begin_reduction
            t_phase = tic()
            if( job%l_ml_solve .and. params%l_trail_rec )then
                fname = refine3D_pcg_trail_accum_fname(state_here, half)
                call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                    &chain_provenance, 1.0, job%nptcls)
                call fname%kill
                if( l_bootstrap .or. 1.0-update_weights(state_here) <= 0.01 )then
                    call job%pcgop%scale_raw_accum(realized_fractions(state_here))
                endif
            else
                do part_here = 1, params%nparts
                    fname = refine3D_pcg_raw_accum_fname(state_here, part_here, params%numlen, half)
                    call job%pcgop%add_raw_accum(fname, state_here, eo_here, part_here, params%nparts, &
                        &provenance, n_part)
                    job%nptcls = job%nptcls + n_part
                    call fname%kill
                enddo
            endif
            job%time_reduce = real(toc(t_phase),dp)
            if( job%nptcls == 0 .and. .not. l_frozen_rec )then
                call job%pcgop%kill
                return
            endif
            if( job%nptcls == 0 )then
                ! an add-on half without a cohort particle follows the trailing
                ! recurrence (solve3D_addon_policy.md, section 6)
                if( .not. job%l_ml_solve )then
                    fname = refine3D_pcg_trail_accum_fname(state_here, half)
                    if( params%l_trail_rec )then
                        if( .not. file_exists(fname) ) &
                            &THROW_HARD('solve3D_addon trailing assembly requires a seeded cohort chain')
                        if( realized_fractions(state_here) < 0.001 )then
                            call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                                &chain_provenance, 1.0, n_part)
                        else
                            if( 1.0-update_weights(state_here) > 0.01 ) &
                                &call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                                &chain_provenance, chain_weights(state_here), n_part)
                            call job%pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, &
                                &count_full_state_half(state_here, eo_here), chain_provenance)
                        endif
                    else if( pcg_trail_seed_requested(cline) )then
                        call job%pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, &
                            &count_full_state_half(state_here, eo_here), chain_provenance)
                    endif
                    call fname%kill
                endif
            else if( .not. job%l_ml_solve )then
                n_full_half = count_full_state_half(state_here, eo_here)
                if( n_full_half < job%nptcls ) &
                    &THROW_HARD('PCG current half population exceeds its full population')
                realized_fraction = realized_fractions(state_here)
                update_weight = update_weights(state_here)
                l_seed_chain = pcg_trail_seed_requested(cline)
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
                    call job%pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, n_full_half, &
                        &chain_provenance)
                    if( .not. l_chain_exists .or. 1.0-update_weight <= 0.01 ) &
                        &call job%pcgop%scale_raw_accum(realized_fraction)
                    write(logfhandle,'(A,I0,A,A,A,F8.4,A,F8.4)') '>>> PCG TRAIL | STATE=', state_here, &
                        &' | HALF=', trim(half), ' | F=', realized_fraction, ' | U=', update_weight
                else if( l_seed_chain )then
                    call job%pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, n_full_half, &
                        &chain_provenance)
                endif
                call fname%kill
                ! solve3D_addon producer: the frozen particles' raw pair at this box
                if( l_frozen_seed ) call frozen_ctx%write_pcg_half(state_here, eo_here, params%box_crop, &
                    &job%pcgop, job%nptcls)
            endif
            ! solve3D_addon: the frozen raw pair joins both reductions after the
            ! chain write, before end_accum (solve3D_addon_policy.md, section 7)
            if( l_frozen_rec )then
                call frozen_ctx%add_pcg_half(state_here, eo_here, params%box_crop, params%smpd_crop, &
                    &job%pcgop, job%nfrozen)
                if( job%nptcls + job%nfrozen == 0 )then
                    call job%pcgop%kill
                    return
                endif
            endif

            t_phase = tic()
            if( job%l_ml_solve )then
                call job%pcgop%set_ml_prior(fsc_prior, params%tau, params%hp)
                job%band_shell = get_find_at_crit(fsc_prior, 0.143)
            endif
            call job%pcgop%end_accum(.true.)
            call job%pcgop%set_op_mode(PCG_OP_KERNEL)
            job%time_finalize = real(toc(t_phase),dp)
            job%prior_npositive = 0
            job%prior_positive_min = 0.0
            job%prior_positive_max = 0.0
            job%prior_to_khat_l1 = 0.0
            job%prior_to_khat_rms = 0.0
            if( job%l_ml_solve )then
                call job%pcgop%get_ml_prior_stats(job%prior_npositive, job%prior_positive_min, &
                    &job%prior_positive_max, job%prior_to_khat_l1, job%prior_to_khat_rms)
                ! the regularized map is the closed-form P_tau optimum of the
                ! current base solution (shrink_by_ml_prior); nothing is solved
                job%x = warm_start%get_rmat()
                job%l_nonzero = .false.
            else
                ! the base solve starts from zero
                allocate(job%x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
                job%l_nonzero = .false.
            endif
            job%ready = .true.
        end subroutine prepare_distributed_half_job

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
            call job%pcgop%kill
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
            write(funit,'(A,ES14.6)') 'requested_rtol=',        params%rtol
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
        if( trim(params%pcgop) /= 'kernel' ) THROW_HARD('production rec_backend=pcg requires pcgop=kernel')
        if( l_check_solver )then
            if( params%maxits_pcg < 1 .or. params%maxits_pcg > 100 ) &
                &THROW_HARD('PCG requires 1<=maxits_pcg<=100')
            if( params%maxits_ml < 0 .or. params%maxits_ml > 100 ) &
                &THROW_HARD('PCG requires 0<=maxits_ml<=100 (0: closed-form regularized pair only)')
            if( params%maxits_pcg > 8 )then
                THROW_WARN('maxits_pcg exceeds the production refinement budget (8); appropriate for offline converged solves only')
            endif
        endif
        if( trim(params%projrec) /= 'no' ) THROW_HARD('rec_backend=pcg does not yet support projrec=yes')
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
    !! v3: the header particle count of each half is its represented population,
    !! the M(s) of the population rule; an older chain is discarded and re-seeded
    function pcg_chain_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgtrail-v3|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_chain_provenance

end module simple_rec3D_pcg_strategy
