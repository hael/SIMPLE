!@descr: shared-memory production strategy body for kernel PCG reconstruct3D
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_STOP_INDEFINITE, &
    &pcg_raw_accum_compatible
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, &
    &discrete_read_imgbatch_source, killimgbatch, prep_rec_observation
use simple_sigma2_files,      only: load_sigma2_groups
use simple_math_ft,           only: resample_sigma2
use simple_estimate_ssnr,     only: fsc2shrink_filter
use simple_image,             only: image
use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, evaluate_halfmap_pair, &
    &write_halfmap_diagnostics, &
    &write_support_provenance, read_support_provenance
use simple_image_msk,         only: image_msk
use simple_nu_filter,         only: NU_DEV_OUTPUT
use simple_nu_state_filter,   only: nonuniform_filter_state, nu_static_aux_replacement
use simple_refine3D_fnames,   only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname, &
    &refine3D_pcg_trail_accum_fname
!$ use omp_lib, only: omp_get_max_threads, omp_get_num_procs, omp_get_max_active_levels, &
!$     &omp_set_max_active_levels, omp_set_num_threads
implicit none

public :: execute_rec3D_pcg_shared, execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
public :: validate_rec3D_pcg_fractional_updates, rec3D_master_nthr
private
#include "simple_local_flags.inc"

real,    parameter :: PCG_LAMBDA = 1.0e-3
!> a nonzero start whose initial relative residual exceeds this is discarded
!! by the solver for a zero start (which has exactly 1.0), before iterating
real,    parameter :: PCG_START_MAX_REL_RESID = 1.0
integer, parameter :: PCG_MASTER_NTHR_CAP    = 32   !< master-phase thread-boost ceiling
! Solve-support envelope (pcg_priors_history.md dev item 5): the conservative density
! envelope (automask3D at envmsklp) replaces the spherical support in the PCG
! solves, so the mask constrains the ESTIMATOR rather than post-processing.
! The envelope is an automsk=yes feature (policy 2026-09-06) and, once a
! prior reconstruction exists, constrains BOTH the base/unfil and the
! regularized pass (policy 2026-09-09); envfsc=yes is derived from
! automsk=yes in parameters (both backends), so here the FSC pair is
! reported as support-constrained and the phase-randomized correction is
! skipped, while gridding applies the same envelope post hoc with it.
! A first reconstruction has no density source and necessarily bootstraps
! the base on the sphere; its current pair then supplies the replay
! support. With automsk=no every solve runs on the sphere. The envelope
! is generous by construction --
! protein plus micelle, dilated, soft-edged, and about half the spherical
! support on the datasets measured so far -- so it removes the far solvent,
! where deapodization amplifies hardest, while leaving the NU evidence a
! populated support. The evidence readiness contract (null_fraction bounds)
! is the guard if a specimen's envelope ever proves too tight for the null
! calibration.
logical, parameter :: DEBUG = .false.


contains


    ! No cross-iteration warm starts (policy 2026-09-10, PfCRT regression
    ! gridding_vs_pcg): the base solve starts from zero and the ML replay from
    ! the shell-shrunk current base solution, every iteration. Warm-starting
    ! from the previous iteration's half maps carried an unconverged transient
    ! of a slightly different system (new FSC prior, new poses) into a
    ! 2-iteration budget that cannot pull it back; the relative residual of the
    ! replay then grew iteration over iteration to 10^3-10^4 in the failed
    ! runs and never fell below 1 in any run. A fixed small budget from a
    ! state-free start is what the simulated-data calibration validated
    ! (2 iterations beat gridding; beyond 5 the residual moves but nothing
    ! interpretable in the map does). The support-provenance solve kind is
    ! still written for the trailing bootstrap's lag-one FSC pair.

    !> Solve with one cold restart: a solve from a NONZERO start (l_nonzero,
    !! the ML replay's shell-shrunk base) that loses positive-definiteness is
    !! retried once from zero. A nonzero start that is WORSE THAN ZERO
    !! (initial relative residual above PCG_START_MAX_REL_RESID; zero has
    !! exactly 1) is discarded by the solver itself before the first
    !! iteration, at no cost (2026-09-10, streptavidin canonical set: in the
    !! low-resolution stages the shrinkage start has INIT 2-4 and two
    !! iterations left the replay at RESID 0.4-2, a half-converged transient
    !! that was the matching reference and the symmetry-search input). This
    !! procedure is compute-only so it is safe inside the concurrent half
    !! sections: it touches only the job's own operator, iterate, and
    !! outcomes. Reporting and fatal handling are deferred to
    !! handle_cold_restart_outcome on the serial side.
    subroutine solve_with_cold_restart( pcgop, x, l_nonzero, maxits, rtol, rel_res_hist, niters, outcome )
        type(reconstructor_pcg),  intent(inout) :: pcgop
        real,                     intent(inout) :: x(:,:,:)
        logical,                  intent(in)    :: l_nonzero
        integer,                  intent(in)    :: maxits
        real,                     intent(in)    :: rtol
        real, allocatable,        intent(out)   :: rel_res_hist(:)
        integer,                  intent(out)   :: niters
        type(pcg_solver_outcome), intent(out)   :: outcome
        type(pcg_solver_outcome) :: first_failure
        if( l_nonzero )then
            call pcgop%solve_accum(x, maxits=maxits, rtol=rtol, rel_res_hist=rel_res_hist, niters=niters, &
                &outcome=outcome, start_max_rel_resid=PCG_START_MAX_REL_RESID)
        else
            call pcgop%solve_accum(x, maxits=maxits, rtol=rtol, rel_res_hist=rel_res_hist, niters=niters, &
                &outcome=outcome)
        endif
        if( trim(outcome%stop_reason) /= PCG_STOP_INDEFINITE ) return
        if( .not. l_nonzero ) return
        if( outcome%start_rejected ) return  ! already a zero start: fatal below
        first_failure = outcome
        x = 0.0
        if( allocated(rel_res_hist) ) deallocate(rel_res_hist)
        call pcgop%solve_accum(x, maxits=maxits, rtol=rtol, rel_res_hist=rel_res_hist, niters=niters, &
            &outcome=outcome)
        outcome%cold_restart_used = .true.
        outcome%restart_trigger_curvature = first_failure%failure_curvature
        outcome%restart_trigger_iteration = first_failure%failure_iteration
    end subroutine solve_with_cold_restart

    !> Serial reporting/failure boundary for solve_with_cold_restart. Keeping
    !! log I/O and THROW_HARD outside the distributed OpenMP sections preserves
    !! the prepared-job concurrency contract.
    subroutine handle_cold_restart_outcome( outcome, context, half, solve_kind )
        type(pcg_solver_outcome), intent(in) :: outcome
        character(len=*), intent(in) :: context, half, solve_kind
        character(len=256) :: error_message
        logical :: l_restarted
        if( outcome%start_rejected )then
            write(logfhandle,'(A,ES10.3,A)') '>>> PCG '//trim(context)//' ('//trim(half)//'/'//&
                &trim(solve_kind)//'): start worse than zero (INIT=', outcome%rejected_start_initial, &
                &'); discarded, solved from zero'
        endif
        l_restarted = outcome%cold_restart_used
        if( l_restarted )then
            write(logfhandle,'(A,I0,A,ES12.4,A)') '>>> PCG '//trim(context)//' ('//trim(half)//'/'//&
                &trim(solve_kind)//'): CG from the nonzero start lost positive-definiteness at iteration ', &
                &outcome%restart_trigger_iteration, ' (dot(p,Hp)=', outcome%restart_trigger_curvature, &
                &'); restarted from zero'
        endif
        if( trim(outcome%stop_reason) /= PCG_STOP_INDEFINITE ) return
        if( l_restarted )then
            error_message = 'PCG lost positive-definiteness from the cold restart as well; '//&
                &trim(context)//'/'//trim(half)//'/'//trim(solve_kind)
        else
            write(logfhandle,'(A,I0,A,ES12.4,A)') '>>> PCG '//trim(context)//' ('//trim(half)//'/'//&
                &trim(solve_kind)//'): cold CG lost positive-definiteness at iteration ', &
                &outcome%failure_iteration, ' (dot(p,Hp)=', outcome%failure_curvature, ')'
            error_message = 'PCG lost positive-definiteness from a cold start; '//trim(context)//'/'//&
                &trim(half)//'/'//trim(solve_kind)
        endif
        THROW_HARD(error_message)
    end subroutine handle_cold_restart_outcome

    !> Regularized initial guess of every ML replay (no cross-iteration warm
    !! starts, 2026-09-10). The documented cold-start gap
    !! is that the regularized optimum differs from ANY unregularized map in
    !! slowly-converging directions -- beyond-band/high-shell noise the ML
    !! prior shrinks -- so a
    !! small iteration budget leaves the solve transient-dominated. Apply
    !! the regularizers' expected effect to the base solution in closed form
    !! instead: shrink each shell's amplitude by the FSC (the Wiener
    !! shrinkage the ML prior's optimum implies; no shrinkage below the hp
    !! no-prior limit, zero beyond the measured band). CG then corrects an
    !! approximate optimum rather than constructing it. Pure initialization:
    !! the quadratic objective has a unique optimum, so this changes the
    !! convergence path, never the converged solution.
    subroutine regularized_ml_initial_guess( params, fsc, x, context, half )
        class(parameters), intent(in)    :: params
        real,              intent(in)    :: fsc(:)
        real,              intent(inout) :: x(:,:,:)
        character(len=*),  intent(in)    :: context, half
        type(image) :: img
        real, allocatable :: filt(:)
        call img%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call img%set_rmat(x, .false.)
        call ml_shrinkage_filter(params, fsc, img%get_filtsz(), filt)
        call img%fft()
        call img%apply_filter(filt)
        call img%ifft()
        x = img%get_rmat()
        call img%kill
        deallocate(filt)
        write(logfhandle,'(A)') '>>> PCG ML REGULARIZED INIT ('//trim(context)//'/'//trim(half)//&
            &'): shell-shrunk base'
    end subroutine regularized_ml_initial_guess

    !> The strategy-side wrapper for the ML shrinkage filter (fsc2shrink_filter
    !! in simple_estimate_ssnr, alongside the other FSC-derived filters):
    !! resolves the no-shrinkage high-pass index from the workflow parameters.
    subroutine ml_shrinkage_filter( params, fsc, nyq, filt )
        class(parameters), intent(in)  :: params
        real,              intent(in)  :: fsc(:)
        integer,           intent(in)  :: nyq
        real, allocatable, intent(out) :: filt(:)
        integer :: k_hp
        k_hp = max(1, calc_fourier_index(params%hp, params%box_crop, params%smpd_crop))
        call fsc2shrink_filter(fsc, k_hp, nyq, filt)
    end subroutine ml_shrinkage_filter


    !> PCG half-map diagnostics through the backend-neutral common evaluator:
    !! this owns the PCG mask policy (params%msk_crop spherical radius on the
    !! restored real-space base pair), the automask artifact write on the
    !! envfsc path, and the execution-context resolution log lines. The shared
    !! and distributed paths pass identical scientific policy and differ only
    !! in the context label.
    subroutine calculate_pcg_state_diagnostics( params, state_here, context, even, odd, avg, diagnostics, &
        &l_pair_support_constrained )
        class(parameters),                intent(in)  :: params
        integer,                          intent(in)  :: state_here
        character(len=*),                 intent(in)  :: context
        class(image),                     intent(in)  :: even, odd, avg
        type(halfmap_diagnostics_result), intent(out) :: diagnostics
        logical,                          intent(in)  :: l_pair_support_constrained
        type(image) :: envmask
        if( params%l_envfsc )then
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, diagnostics, envmask=envmask, &
                &l_pair_support_constrained=l_pair_support_constrained)
            call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
            call envmask%kill
        else
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, diagnostics, &
                &l_pair_support_constrained=l_pair_support_constrained)
        endif
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.500 RESOLUTION = ', diagnostics%res_fsc05
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.143 RESOLUTION = ', diagnostics%res_fsc0143
    end subroutine calculate_pcg_state_diagnostics




    ! The support-provenance sidecar (write/read_support_provenance) lives in
    ! simple_halfmap_diagnostics since both backends write it (2026-09-09).

    !> Resolution-text naming, mirroring the gridding volassemble contract
    !! (resolve_fsc_txt_fname in simple_commanders_rec_distr): an explicit
    !! outfile names the document (the final-reconstruction RESOLUTION_FINAL
    !! case), which_iter tags per-iteration documents so refinement
    !! iterations do not overwrite one another, and the plain state name is
    !! the fallback.
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

    !> Install the solve support constraint, in precedence order: an explicit
    !! pcg_mskfile envelope, the per-state density-envelope support when one
    !! was built for this state, else the spherical mskdiam
    !! support. The projected system (P H P) u = P b is what set_mask already
    !! implements; this only chooses P (pcg_priors_history.md dev item 5).
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

    !> Build the per-state solve-support envelope from the reference volume
    !! this iteration matched against (lag-one, the same lag the matching
    !! references carry). A start volume is a valid density source. With no
    !! usable reconstruction yet, the base solve necessarily bootstraps on the
    !! sphere; its current pair then supplies density support for the replay.
    subroutine build_pcg_state_support( params, state_here, support, l_have )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        type(image_msk),   intent(inout) :: support
        logical,           intent(out)   :: l_have
        type(image)  :: vol_prev
        l_have = .false.
        call support%kill_bimg
        ! An explicit pcg_mskfile (development/focused mode) constrains every
        ! solve regardless of automsk, and is therefore reported as the state
        ! support so the FSC mode, the support provenance and the NU evidence
        ! null regime all see a constrained pair (review 2026-09-09). An
        ! allocated empty string is not an override.
        if( params%pcg_mskfile%is_allocated() )then
            if( len_trim(params%pcg_mskfile%to_char()) > 0 )then
                call support%read_and_crop(params%pcg_mskfile, params%smpd, params%box_crop, params%smpd_crop)
                write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                    &', explicit pcg_mskfile '//params%pcg_mskfile%to_char()//' constrains base and replay'
                l_have = .true.
                return
            endif
        endif
        ! Density solve support is an automsk feature (policy 2026-09-06,
        ! reversing the 2026-09-02 "independent of automsk" review item): no
        ! envelope constrains any PCG solve unless automsk=yes. With automsk
        ! enabled BOTH the base/unfil and the regularized solve take the
        ! envelope (policy 2026-09-09); envfsc=yes is derived from that in
        ! parameters, never requested separately. Without a lagged reference
        ! the replay support is built from the completed current base pair
        ! (the base then bootstraps on the sphere). With automsk=no every
        ! solve, base and replay, runs on the sphere.
        if( .not. pcg_density_support_enabled(params) )then
            write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                &', spherical support for base and replay (automsk=no; density envelope disabled)'
            return
        endif
        if( .not. params%l_envfsc ) &
            &THROW_HARD('automsk=yes requires envfsc=yes (derived in parameters); the coupling was bypassed')
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

    contains

        subroutine handle_missing_reference( why )
            character(len=*), intent(in) :: why
            write(logfhandle,'(A,I0,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
                &' has '//trim(why)//'; bootstrap base uses the sphere and replay support derives from that base pair'
        end subroutine handle_missing_reference

    end subroutine build_pcg_state_support

    !> Density solve support is permitted only under automsk=yes (policy
    !! 2026-09-06). An explicit pcg_mskfile is the development override and is
    !! installed by set_pcg_solve_support regardless of this gate.
    logical function pcg_density_support_enabled( params ) result( l_enabled )
        class(parameters), intent(in) :: params
        l_enabled = trim(params%automsk) .ne. 'no'
    end function pcg_density_support_enabled

    !> Construct the conservative density support from an explicit volume.
    !! Used for the normal lag-one path and, when no prior reference exists,
    !! for the regularized replay after the spherical base pair is available.
    !! Both callers are gated on automsk=yes. The NU-evidence envelope never
    !! enters this routine.
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





    !> Master-phase thread budget shared by both reconstruction backends
    !! (2026-09-09): on local execution the partition workers are idle while
    !! the master assembles (PCG solves, or the gridding restoration and the
    !! NU competition), so the full allocation is used, capped where OpenMP
    !! scaling saturates at these box sizes; on a cluster the master's slot is
    !! fixed at nthr_master. One rule for both backends keeps their
    !! master-phase timings comparable.
    integer function rec3D_master_nthr( params, nthr_master ) result( nthr )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nthr_master
        nthr = nthr_master
        if( trim(params%qsys_name) == 'local' )then
            nthr = max(params%nthr, min(PCG_MASTER_NTHR_CAP, max(1, params%nparts) * params%nthr))
            !$ nthr = min(omp_get_num_procs(), nthr)
        endif
    end function rec3D_master_nthr

    subroutine execute_rec3D_pcg_shared( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(image) :: half_even, half_odd, ml_even, ml_odd, merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt
        type(halfmap_diagnostics_result) :: hm_diag
        integer, allocatable :: selected_pinds(:), half_pinds(:)
        real, allocatable :: fsc(:), res0143s(:), res05s(:), nu_align_lps(:)
        logical, allocatable :: state_written(:)
        type(string) :: eonames(2)
        integer :: nselected, state, n_state, n_even, n_odd, iptcl, istate
        type(image_msk) :: state_support_msk
        logical :: l_state_support, l_base_support_constrained
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output, time_nu_filter
        real :: align_lp
        logical :: l_sigma_loaded

        call validate_supported_mode()
        nselected = 0
        call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
        if( nselected < 1 ) THROW_HARD('no active particles selected for PCG reconstruct3D')

        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, &
                &build%spproj_field, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif

        allocate(res0143s(params%nstates), source=0.0)
        allocate(res05s(params%nstates),   source=0.0)
        allocate(nu_align_lps(params%nstates), source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)

        do state = 1, params%nstates
            n_state = count_state(state)
            n_even = count_state_half(state, 0)
            n_odd  = count_state_half(state, 1)
            if( n_state == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG RECONSTRUCT3D: STATE ', state, ' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even + n_odd /= n_state ) THROW_HARD('PCG reconstruct3D found invalid halfset labels')
            if( n_even < 1 .or. n_odd < 1 ) THROW_HARD('PCG reconstruct3D requires particles in both halfsets')

            ! One density-envelope support per state, shared with the
            ! distributed owner, and only under automsk=yes: both the base
            ! and the regularized replay take it once a prior reconstruction
            ! exists (envfsc=yes follows). automsk=no: sphere throughout.
            call build_pcg_state_support(params, state, state_support_msk, l_state_support)
            l_base_support_constrained = l_state_support
            call collect_state_half(state, 0, n_even, half_pinds)
            call solve_state_half(state, 0, 'even', half_pinds, half_even)
            deallocate(half_pinds)
            call collect_state_half(state, 1, n_odd, half_pinds)
            call solve_state_half(state, 1, 'odd', half_pinds, half_odd)
            deallocate(half_pinds)

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

            t_state_phase = tic()
            call calculate_pcg_state_diagnostics(params, state, 'RECONSTRUCT3D', half_even, half_odd, &
                &merged, hm_diag, l_base_support_constrained)
            fsc             = hm_diag%fsc
            res0143s(state) = hm_diag%res_fsc0143
            res05s(state)   = hm_diag%res_fsc05
            call arr2file(fsc, fname_fsc)
            fname_restxt = resolve_pcg_fsc_txt_fname(params, cline, state)
            call write_halfmap_diagnostics(hm_diag, params%box_crop, params%smpd_crop, fname_restxt)
            call fname_restxt%kill
            call hm_diag%kill
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                ! the ordinary global-ML replay (P_tau from the current base
                ! pair); nonuniform filtering is assembly-owned and runs
                ! post hoc below, exactly as on the gridding backend
                call regularize_state_half(state, 0, 'even', fsc, half_even, ml_even)
                call regularize_state_half(state, 1, 'odd',  fsc, half_odd,  ml_odd)
                call merged%kill
                call merged%copy(ml_even)
                call merged%add(ml_odd)
                call merged%mul(0.5)
            endif
            t_state_phase = tic()
            if( params%l_ml_reg )then
                call ml_even%write(fname_even, del_if_exists=.true.)
                call ml_odd%write(fname_odd, del_if_exists=.true.)
            else
                call half_even%write(fname_even, del_if_exists=.true.)
                call half_odd%write(fname_odd, del_if_exists=.true.)
            endif
            call merged%write(fname_vol, del_if_exists=.true.)
            ! the sidecar follows the published map, never precedes it
            if( params%l_ml_reg )then
                call write_support_provenance(fname_vol, l_state_support, 'regularized')
            else
                call write_support_provenance(fname_vol, l_base_support_constrained, 'base')
            endif
            time_map_output = time_map_output + real(toc(t_state_phase),dp)
            if( params%l_nonuniform )then
                ! NU competition, mirroring the gridding volassemble: the
                ! candidate bank from the base (unfil) pair, the ML pair as
                ! the auxiliary member; consumes both pairs
                t_state_phase = tic()
                eonames(1) = fname_even
                eonames(2) = fname_odd
                ! an envelope-constrained base pair hands its support over so
                ! the evidence null is designated on the dilation ring
                call nonuniform_filter_state(params, state, half_even, half_odd, &
                    &ml_even, ml_odd, params%l_ml_reg .and. nu_static_aux_replacement(params), &
                    &res0143s(state), fname_vol, eonames, nu_align_lps(state), &
                    &base_support=state_support_msk, l_base_constrained=l_base_support_constrained)
                time_nu_filter = real(toc(t_state_phase),dp)
            endif
            call write_output_diagnostics(state, 'shared', time_map_output, time_fsc_output, time_nu_filter)

            params%vols(state)      = fname_vol
            params%vols_even(state) = fname_even
            params%vols_odd(state)  = fname_odd
            call cline%set('vol'//int2str(state), fname_vol)
            state_written(state) = .true.
            call half_even%kill
            call half_odd%kill
            if( params%l_ml_reg )then
                call ml_even%kill
                call ml_odd%kill
                raw_fname = refine3D_pcg_raw_accum_fname(state, 1, params%numlen, 'even')
                call del_file(raw_fname)
                raw_fname = refine3D_pcg_raw_accum_fname(state, 1, params%numlen, 'odd')
                call del_file(raw_fname)
                call fname_even_unfil%kill
                call fname_odd_unfil%kill
                call raw_fname%kill
            endif
            call merged%kill
            call fname_even%kill
            call fname_odd%kill
            call fname_vol%kill
            call fname_fsc%kill
            if( allocated(fsc) ) deallocate(fsc)
        enddo

        call state_support_msk%kill_bimg
        call killimgbatch(build)
        if( .not. any(state_written) ) THROW_HARD('PCG reconstruct3D produced no populated states')
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res',   res0143s(1))
            call build%spproj_field%set_all2single('res05', res05s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) )then
                        call build%spproj_field%set(iptcl, 'res',   res0143s(istate))
                        call build%spproj_field%set(iptcl, 'res05', res05s(istate))
                    endif
                endif
            enddo
        endif
        if( params%l_nonuniform )then
            ! the best resolved populated state sets the single matching band
            if( any(state_written .and. nu_align_lps > TINY) )then
                align_lp = minval(nu_align_lps, mask=state_written .and. nu_align_lps > TINY)
                call build%spproj_field%set_all2single('lp',     align_lp)
                call build%spproj_field%set_all2single('lp_est', align_lp)
                write(logfhandle,'(A,F8.3,A)') '>>> NU filter matching low-pass limit for next iteration: ', &
                    &align_lp, ' A'
            endif
        endif
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call register_project_outputs()

        deallocate(selected_pinds, res0143s, res05s, nu_align_lps, state_written)

    contains

        subroutine validate_supported_mode()
            if( params%nparts > 1 .or. cline%defined('part') )then
                THROW_HARD('shared PCG entry received distributed parameters')
            endif
            if( trim(params%pcgop) /= 'kernel' ) THROW_HARD('production rec_backend=pcg requires pcgop=kernel')
            ! refinement stays within the small fixed-iteration budget the
            ! simulated-data calibration validated; offline harness/benchmark runs legitimately ask
            ! for rtol-terminated converged solves (pcg_priors_history.md R3), so above
            ! the production budget we warn rather than refuse
            if( params%maxits_pcg > 100 ) THROW_HARD('maxits_pcg > 100 is not supported')
            if( params%maxits_pcg > 8 )then
                THROW_WARN('maxits_pcg exceeds the production refinement budget (8); appropriate for offline converged solves only')
            endif
            if( trim(params%projrec) /= 'no' ) THROW_HARD('rec_backend=pcg does not yet support projrec=yes')
            if( abs(real(params%box)*params%smpd - real(params%box_crop)*params%smpd_crop) > &
                &1.0e-5*real(params%box)*params%smpd )then
                THROW_HARD('PCG crop must preserve the native physical box extent')
            endif
            if( params%l_update_frac .or. params%l_trail_rec )then
                THROW_HARD('PCG fractional and trailing reconstruction require accumulator-domain integration')
            endif
            if( trim(params%conical_fsc) == 'yes' ) THROW_HARD('PCG conical FSC integration is not implemented')
            if( params%msk <= 0.5 .or. params%msk_crop <= 0.5 )then
                THROW_HARD('rec_backend=pcg requires mskdiam for normalization and solve support')
            endif
            if( params%maxits_pcg < 1 ) THROW_HARD('PCG maxits_pcg must be at least 1')
            if( .not. ieee_is_finite(params%rtol) ) THROW_HARD('PCG rtol must be finite')
        end subroutine validate_supported_mode

        integer function count_state( state_here ) result(n)
            integer, intent(in) :: state_here
            integer :: i
            n = 0
            do i = 1, nselected
                if( build%spproj_field%get_state(selected_pinds(i)) == state_here ) n = n + 1
            enddo
        end function count_state

        integer function count_state_half( state_here, eo_here ) result(n)
            integer, intent(in) :: state_here, eo_here
            integer :: i, p
            n = 0
            do i = 1, nselected
                p = selected_pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) == eo_here ) n = n + 1
            enddo
        end function count_state_half

        subroutine collect_state_half( state_here, eo_here, n, pinds )
            integer,              intent(in)  :: state_here, eo_here, n
            integer, allocatable, intent(out) :: pinds(:)
            integer :: i, p, cnt
            allocate(pinds(n))
            cnt = 0
            do i = 1, nselected
                p = selected_pinds(i)
                if( build%spproj_field%get_state(p) /= state_here ) cycle
                if( build%spproj_field%get_eo(p) /= eo_here ) cycle
                cnt = cnt + 1
                pinds(cnt) = p
            enddo
            if( cnt /= n ) THROW_HARD('inconsistent PCG state-half particle count')
        end subroutine collect_state_half

        subroutine solve_state_half( state_here, eo_here, half, pinds, volume, outcome )
            integer,          intent(in)    :: state_here, eo_here
            character(len=*), intent(in)    :: half
            integer,          intent(in)    :: pinds(:)
            type(image),      intent(inout) :: volume
            type(pcg_solver_outcome), optional, intent(out) :: outcome
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            type(oris)      :: selection
            type(ori)       :: orientation
            type(ctfparams) :: ctfparms
            type(string)    :: raw_fname_here
            type(image)     :: obs
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:), x(:,:,:), rel_res_hist(:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch, niters
            real    :: shift(2), crop_factor
            integer(timer_int_kind) :: t_half, t_phase
            real(dp) :: time_metadata, time_particles, time_accum_init, time_accum
            real(dp) :: time_finalize, time_solve, time_total

            t_half = tic()
            time_metadata  = 0.0_dp
            time_particles = 0.0_dp
            time_accum_init = 0.0_dp
            time_accum     = 0.0_dp
            time_finalize  = 0.0_dp
            time_solve     = 0.0_dp
            t_phase = tic()
            crop_factor = real(params%box_crop) / real(params%box)
            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call pcgop%set_sym(build%pgrpsyms)
            ! base/unregularized solve: the density support whenever one was
            ! built for this state (automsk=yes with a prior reconstruction)
            call set_pcg_solve_support(pcgop, params, state_support_msk, l_state_support)
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
            time_metadata = real(toc(t_phase),dp)
            t_phase = tic()
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call pcgop%begin_accum
            time_accum_init = real(toc(t_phase),dp)
            call obs%new([params%box_crop,params%box_crop,1], params%smpd_crop)
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                t_phase = tic()
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
                do ii = 1, batchsz
                    ! the backend-neutral observation (normalize, crop, taper), see prep_rec_observation
                    call prep_rec_observation(build%imgbatch(ii), build%lmsk, obs, .true.)
                    call obs%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(obs)
                enddo
                time_particles = time_particles + real(toc(t_phase),dp)
                t_phase = tic()
                call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
                time_accum = time_accum + real(toc(t_phase),dp)
            enddo
            if( params%l_ml_reg )then
                raw_fname_here = refine3D_pcg_raw_accum_fname(state_here, 1, params%numlen, half)
                call pcgop%write_raw_accum(raw_fname_here, state_here, eo_here, 1, 1, &
                    &size(pinds), pcg_raw_provenance(params))
                call raw_fname_here%kill
            endif
            if( pcg_trail_seed_requested(cline) )then
                raw_fname_here = refine3D_pcg_trail_accum_fname(state_here, half)
                call pcgop%write_raw_accum(raw_fname_here, state_here, eo_here, 1, 1, &
                    &size(pinds), pcg_chain_provenance(params))
                call raw_fname_here%kill
            endif
            call obs%kill
            t_phase = tic()
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            time_finalize = real(toc(t_phase),dp)
            ! the base solve starts from zero (no cross-iteration warm start)
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            t_phase = tic()
            call solve_with_cold_restart(pcgop, x, .false., params%maxits_pcg, params%rtol, &
                &rel_res_hist, niters, result)
            time_solve = real(toc(t_phase),dp)
            call handle_cold_restart_outcome(result, 'shared', half, 'base')
            call validate_solved_map(x, 'shared', state_here, half, 'base')
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            call report_beyond_band_excess(volume, params, state_here, half, 'base')
            time_total = real(toc(t_half),dp)
            call write_half_diagnostics(state_here, half, 'base', size(pinds), result, rel_res_hist, &
                &time_metadata, time_particles, time_accum_init, time_accum, time_finalize, time_solve, time_total, &
                &pcgop%get_data_scale(), pcgop%get_effective_lambda(), pcgop=pcgop)
            call report_solve_summary('SHARED', state_here, half, 'base', size(pinds), niters, &
                &result%final_rel_residual, time_solve, result%stop_reason, result%initial_rel_residual)
            if( present(outcome) ) outcome = result

            call pcgop%kill
            call selection%kill
            call orientation%kill
            deallocate(y_batch, sig2, x, rel_res_hist)
        end subroutine solve_state_half

        !> Reopen the exact raw statistics used for the base half-map, add the
        !! FSC/SSNR prior only on the master/shared owner, and start from the
        !! corresponding unregularized solution (shell-shrunk). No particle data are read
        !! a second time and the base solve remains the FSC oracle.
        subroutine regularize_state_half( state_here, eo_here, half, fsc_here, base_volume, volume )
            integer,          intent(in)    :: state_here, eo_here
            character(len=*), intent(in)    :: half
            real,             intent(in)    :: fsc_here(:)
            type(image),      intent(in)    :: base_volume
            type(image),      intent(inout) :: volume
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            type(string) :: fname
            real, allocatable :: x(:,:,:), rel_res_hist(:)
            integer :: nptcls, niters, prior_npositive
            integer(timer_int_kind) :: t_phase
            real(dp) :: time_reduce, time_finalize, time_solve, time_total
            real :: prior_positive_min, prior_positive_max, prior_to_khat_l1, prior_to_khat_rms

            t_phase = tic()
            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            ! regularized replay: always takes the density support when built
            call set_pcg_solve_support(pcgop, params, state_support_msk, l_state_support)
            call pcgop%begin_reduction
            fname = refine3D_pcg_raw_accum_fname(state_here, 1, params%numlen, half)
            call pcgop%add_raw_accum(fname, state_here, eo_here, 1, 1, &
                &pcg_raw_provenance(params), nptcls)
            time_reduce = real(toc(t_phase),dp)
            if( nptcls < 1 ) THROW_HARD('PCG ML replay requires a populated raw half accumulator')
            t_phase = tic()
            ! the global FSC/SSNR replay precision P_tau; its effective
            ! strength is derived from the data scale in end_accum alongside
            ! the relative ridge lambda
            call pcgop%set_ml_prior(fsc_here, params%tau, params%hp)
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            time_finalize = real(toc(t_phase),dp)
            prior_npositive    = 0
            prior_positive_min = 0.0
            prior_positive_max = 0.0
            prior_to_khat_l1   = 0.0
            prior_to_khat_rms  = 0.0
            call pcgop%get_ml_prior_stats(prior_npositive, prior_positive_min, prior_positive_max, &
                &prior_to_khat_l1, prior_to_khat_rms)
            ! the replay starts from the current base solution with the
            ! closed-form shrinkage initial guess (encodes the P_tau optimum);
            ! no cross-iteration warm start
            x = base_volume%get_rmat()
            call regularized_ml_initial_guess(params, fsc_here, x, 'shared', half)
            t_phase = tic()
            ! the replay iterate is never zero: restart-eligible
            call solve_with_cold_restart(pcgop, x, .true., params%maxits_pcg, params%rtol, &
                &rel_res_hist, niters, result)
            time_solve = real(toc(t_phase),dp)
            call handle_cold_restart_outcome(result, 'shared', half, 'ml')
            call validate_solved_map(x, 'shared', state_here, half, 'ml')
            time_total = time_reduce + time_finalize + time_solve
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            call report_beyond_band_excess(volume, params, state_here, half, 'ml')
            call write_half_diagnostics(state_here, half, 'ml', nptcls, result, rel_res_hist, &
                &0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, time_finalize, time_solve, time_total, &
                &pcgop%get_data_scale(), pcgop%get_effective_lambda(), reduce_time=time_reduce, &
                &prior_npositive=prior_npositive, prior_positive_min=prior_positive_min, &
                &prior_positive_max=prior_positive_max, prior_to_khat_l1=prior_to_khat_l1, &
                &prior_to_khat_rms=prior_to_khat_rms, pcgop=pcgop)
            call report_solve_summary('SHARED', state_here, half, 'ml', nptcls, niters, &
                &result%final_rel_residual, time_solve, result%stop_reason, result%initial_rel_residual)
            call pcgop%kill
            call fname%kill
            deallocate(x, rel_res_hist)
        end subroutine regularize_state_half

        subroutine write_half_diagnostics( state_here, half, solve_kind, nptcls, result, history, &
                &metadata, particles, accum_init, accum, finalize, solve_time, total, data_scale, lambda_eff, &
                &reduce_time, prior_npositive, prior_positive_min, prior_positive_max, prior_to_khat_l1, &
                &prior_to_khat_rms, pcgop )
            integer,                  intent(in) :: state_here, nptcls
            character(len=*),         intent(in) :: half, solve_kind
            type(pcg_solver_outcome), intent(in) :: result
            real,                     intent(in) :: history(:)
            real(dp),                 intent(in) :: metadata, particles, accum_init, accum
            real(dp),                 intent(in) :: finalize, solve_time, total
            real,                     intent(in) :: data_scale, lambda_eff
            real(dp), optional,       intent(in) :: reduce_time
            integer, optional,        intent(in) :: prior_npositive
            real, optional,           intent(in) :: prior_positive_min, prior_positive_max
            real, optional,           intent(in) :: prior_to_khat_l1, prior_to_khat_rms
            type(reconstructor_pcg),   intent(in) :: pcgop
            type(string) :: fname
            integer :: funit, i
            real(dp) :: other_time
            fname = 'reconstruct3D_pcg_state'//int2str_pad(state_here,2)//'_'//trim(half)//'_'// &
                &trim(solve_kind)//'_diagnostics.txt'
            call fopen(funit, file=fname, status='replace', action='write')
            write(funit,'(A)')        'execution_mode=shared'
            write(funit,'(A,I0)')     'nptcls=',               nptcls
            write(funit,'(A,A)')      'solve_kind=',            trim(solve_kind)
            write(funit,'(A,I0)')     'requested_maxits=',     result%requested_maxits
            write(funit,'(A,ES14.6)') 'requested_rtol=',       params%rtol
            write(funit,'(A,I0)')     'iteration_count=',      result%iteration_count
            write(funit,'(A,A)')      'stop_reason=',          trim(result%stop_reason)
            write(funit,'(A,L1)')     'converged=',            result%converged
            write(funit,'(A,L1)')     'cold_restart_used=',    result%cold_restart_used
            if( result%cold_restart_used )then
                write(funit,'(A,I0)')     'restart_trigger_iteration=', result%restart_trigger_iteration
                write(funit,'(A,ES14.6)') 'restart_trigger_curvature=', result%restart_trigger_curvature
            endif
            write(funit,'(A,ES14.6)') 'initial_rel_resid_l2=', result%initial_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_resid_l2=',   result%final_rel_residual
            write(funit,'(A,ES14.6)') 'final_rel_update=',     result%final_rel_update
            write(funit,'(A,ES14.6)') 'pcg_data_scale=',       data_scale
            write(funit,'(A,ES14.6)') 'pcg_lambda_effective=', lambda_eff
            if( present(prior_npositive) )then
                if( .not. present(prior_positive_min) .or. .not. present(prior_positive_max) .or. &
                    &.not. present(prior_to_khat_l1) .or. .not. present(prior_to_khat_rms) )then
                    THROW_HARD('incomplete shared PCG ML prior diagnostics')
                endif
                write(funit,'(A,I0)')     'ml_prior_nonzero_bins=',       prior_npositive
                write(funit,'(A,ES14.6)') 'ml_prior_positive_min=',       prior_positive_min
                write(funit,'(A,ES14.6)') 'ml_prior_positive_max=',       prior_positive_max
                write(funit,'(A,ES14.6)') 'ml_prior_to_data_khat_l1=',    prior_to_khat_l1
                write(funit,'(A,ES14.6)') 'ml_prior_to_data_khat_rms=',   prior_to_khat_rms
            endif
            write(funit,'(A,F12.6)')  'metadata_seconds=',     metadata
            write(funit,'(A,F12.6)')  'particle_io_prep_seconds=', particles
            write(funit,'(A,F12.6)')  'accum_init_seconds=',   accum_init
            write(funit,'(A,F12.6)')  'fused_B_D_accum_seconds=', accum
            write(funit,'(A,F12.6)')  'accum_finalize_seconds=', finalize
            write(funit,'(A,F12.6)')  'solve_seconds=',        solve_time
            write(funit,'(A,F12.6)')  'total_half_seconds=',   total
            if( present(reduce_time) ) write(funit,'(A,F12.6)') 'raw_replay_seconds=', reduce_time
            other_time = total - metadata - particles - accum_init - accum - finalize - solve_time
            if( present(reduce_time) ) other_time = other_time - reduce_time
            write(funit,'(A,F12.6)')  'other_seconds=', max(0.0_dp, other_time)
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
        end subroutine write_half_diagnostics

        subroutine register_project_outputs()
            character(len=16) :: imgkind
            integer :: s
            if( trim(params%mkdir) /= 'yes' ) return
            if( trim(params%oritype) == 'cls3D' )then
                imgkind = 'vol_cavg'
            else
                imgkind = 'vol'
            endif
            do s = 1, params%nstates
                if( .not. state_written(s) ) cycle
                fname_vol = refine3D_state_vol_fname(s)
                fname_fsc = refine3D_fsc_fname(s)
                call build%spproj%add_vol2os_out(fname_vol, params%smpd_crop, s, trim(imgkind), box=params%box_crop)
                call build%spproj%add_fsc2os_out(fname_fsc, s, params%box_crop)
                call fname_vol%kill
                call fname_fsc%kill
            enddo
            call build%spproj%write_segment_inside('out', params%projfile)
        end subroutine register_project_outputs

    end subroutine execute_rec3D_pcg_shared

    !> Project-backed pre-integration gate for accumulator-domain fractional
    !! updates. The caller supplies exactly the reconstruction project and
    !! geometry; complementary subsets, realized fractions and continuation
    !! artifacts are generated deterministically here.
    subroutine validate_rec3D_pcg_fractional_updates( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        real, parameter :: RAW_TOL = 2.0e-5, REPLAY_TOL = 2.0e-6
        integer, allocatable :: selected_pinds(:), state_pinds(:), state_subset1(:), state_subset2(:)
        integer, allocatable :: half_pinds(:), subset1(:), subset2(:)
        logical :: l_sigma_loaded
        integer :: nselected, state, eo, n_half, nhalves_tested
        real :: state_f1, state_f2
        character(len=4) :: half
        integer :: funit
        type(string) :: diag_fname

        call validate_pcg_common(params)
        if( cline%defined('part') ) THROW_HARD('pcg_frac_update cannot run as a worker part')
        nselected = 0
        call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
        if( nselected < 1 ) THROW_HARD('no active particles selected for PCG fractional-update validation')
        if( params%cc_objfun == OBJFUN_EUCLID )then
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, &
            &build%spproj_field, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG fractional-update validation requires sigma2 for objfun=euclid')
        endif
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        diag_fname = 'pcg_fractional_update_validation.txt'
        call fopen(funit, file=diag_fname, status='replace', action='write')
        write(funit,'(A)') 'test=pcg_frac_update'
        write(funit,'(A,I0)') 'selected_particles=', nselected
        write(funit,'(A,I0)') 'box_crop=', params%box_crop
        write(funit,'(A,F12.6)') 'smpd_crop=', params%smpd_crop
        write(funit,'(A,I0)') 'maxits_pcg=', params%maxits_pcg
        write(funit,'(A,ES14.6)') 'rtol=', params%rtol
        nhalves_tested = 0

        do state = 1, params%nstates
            call collect_state(state, selected_pinds, state_pinds)
            if( size(state_pinds) == 0 )then
                deallocate(state_pinds)
                cycle
            endif
            call split_complementary(state_pinds, state_subset1, state_subset2)
            state_f1 = real(size(state_subset1)) / real(size(state_pinds))
            state_f2 = real(size(state_subset2)) / real(size(state_pinds))
            do eo = 0, 1
                call collect_state_half(state, eo, state_pinds, half_pinds)
                call collect_state_half(state, eo, state_subset1, subset1)
                call collect_state_half(state, eo, state_subset2, subset2)
                n_half = size(half_pinds)
                if( n_half == 0 )then
                    deallocate(half_pinds, subset1, subset2)
                    cycle
                endif
                if( n_half < 2 ) THROW_HARD('PCG fractional test needs two particles per state/half')
                if( size(subset1) < 1 .or. size(subset2) < 1 )then
                    THROW_HARD('PCG fractional test needs both state subsets represented in each half')
                endif
                half = merge('odd ', 'even', eo == 1)
                call validate_half(state, eo, trim(half), half_pinds, subset1, subset2, &
                    &state_f1, state_f2, funit)
                nhalves_tested = nhalves_tested + 1
                deallocate(half_pinds, subset1, subset2)
            enddo
            deallocate(state_pinds, state_subset1, state_subset2)
        enddo
        if( nhalves_tested < 1 ) THROW_HARD('PCG fractional-update validation found no populated state/half')
        call fclose(funit)
        call killimgbatch(build)
        deallocate(selected_pinds)
        write(logfhandle,'(A,1X,A)') '>>> PCG FRACTIONAL-UPDATE VALIDATION: PASS; DIAGNOSTICS:', diag_fname%to_char()
        call diag_fname%kill

    contains

        subroutine collect_state( state_here, pinds, selected )
            integer,              intent(in)  :: state_here, pinds(:)
            integer, allocatable, intent(out) :: selected(:)
            integer :: i, n
            n = 0
            do i = 1, size(pinds)
                if( build%spproj_field%get_state(pinds(i)) == state_here ) n = n + 1
            enddo
            allocate(selected(n))
            n = 0
            do i = 1, size(pinds)
                if( build%spproj_field%get_state(pinds(i)) /= state_here ) cycle
                n = n + 1
                selected(n) = pinds(i)
            enddo
        end subroutine collect_state

        subroutine collect_state_half( state_here, eo_here, pinds, selected )
            integer,              intent(in)  :: state_here, eo_here, pinds(:)
            integer, allocatable, intent(out) :: selected(:)
            integer :: i, n, p
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) == state_here .and. build%spproj_field%get_eo(p) == eo_here ) n = n + 1
            enddo
            allocate(selected(n))
            n = 0
            do i = 1, size(pinds)
                p = pinds(i)
                if( build%spproj_field%get_state(p) /= state_here .or. build%spproj_field%get_eo(p) /= eo_here ) cycle
                n = n + 1
                selected(n) = p
            enddo
        end subroutine collect_state_half

        subroutine split_complementary( pinds, first, second )
            integer,              intent(in)  :: pinds(:)
            integer, allocatable, intent(out) :: first(:), second(:)
            integer :: i, i1, i2
            ! Deliberately use an unequal one-third/two-thirds split. A 50/50
            ! split with the default u=0.5 would make u/f=1 and fail to exercise
            ! the mass-renormalization step at all.
            allocate(first((size(pinds)+2)/3), second(size(pinds)-(size(pinds)+2)/3))
            i1 = 0
            i2 = 0
            do i = 1, size(pinds)
                if( mod(i-1,3) == 0 )then
                    i1 = i1 + 1
                    first(i1) = pinds(i)
                else
                    i2 = i2 + 1
                    second(i2) = pinds(i)
                endif
            enddo
        end subroutine split_complementary

        subroutine validate_half( state_here, eo_here, half_here, full_pinds, pinds1, pinds2, f1, f2, unit )
            integer,          intent(in) :: state_here, eo_here, full_pinds(:), pinds1(:), pinds2(:), unit
            character(len=*), intent(in) :: half_here
            real,             intent(in) :: f1, f2
            type(reconstructor_pcg) :: op_full, op_subset, op_sum, op_ref
            type(reconstructor_pcg) :: op_blend, op_oracle, op_replay, op_ensemble
            type(string) :: f_full, f_sub1, f_sub2, f_chain1, f_chain2
            real, allocatable :: x_direct(:,:,:), x_replay(:,:,:), x_full(:,:,:), hist(:)
            character(len=256) :: provenance
            real :: u, berr, derr, solve_err, sample_map_err
            integer :: nraw, nraw_total, niters

            provenance = 'pcg-frac-update-v1'
            u  = 0.5
            if( params%l_ufrac_trec_defined ) u = params%ufrac_trec
            if( u <= 0.0 .or. u > 1.0 ) THROW_HARD('pcg_frac_update requires 0 < ufrac_trec <= 1')
            f_full   = frac_fname(state_here, half_here, 'full')
            f_sub1   = frac_fname(state_here, half_here, 'subset1')
            f_sub2   = frac_fname(state_here, half_here, 'subset2')
            f_chain1 = frac_fname(state_here, half_here, 'chain1')
            f_chain2 = frac_fname(state_here, half_here, 'chain2')
            write(unit,'(A,I0,A,A)') 'state=', state_here, ' half=', trim(half_here)
            write(unit,'(A,I0,A,I0,A,I0)') 'n_full=', size(full_pinds), &
                &' n_subset1=', size(pinds1), ' n_subset2=', size(pinds2)
            write(unit,'(A,F10.6,A,F10.6,A,F10.6)') 'state_f1=', f1, ' state_f2=', f2, ' u=', u

            call accumulate_raw(full_pinds, op_full)
            call op_full%write_raw_accum(f_full, state_here, eo_here, 1, 1, size(full_pinds), provenance)
            call op_full%kill
            call accumulate_raw(pinds1, op_subset)
            call op_subset%write_raw_accum(f_sub1, state_here, eo_here, 1, 2, size(pinds1), provenance)
            call op_subset%kill
            call accumulate_raw(pinds2, op_subset)
            call op_subset%write_raw_accum(f_sub2, state_here, eo_here, 2, 2, size(pinds2), provenance)
            call op_subset%kill

            call load_weighted(op_blend, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
            call op_blend%scale_raw_accum(1.0/f1)
            call op_blend%write_raw_accum(f_chain1, state_here, eo_here, 1, 1, size(full_pinds), provenance)
            call op_blend%scale_raw_accum(f1)
            call load_weighted(op_ref, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
            call op_blend%compare_raw_accum(op_ref, berr, derr)
            call require_raw('bootstrap working-mass restore', berr, derr, REPLAY_TOL)
            call op_blend%kill
            call op_ref%kill
            call load_weighted(op_replay, f_chain1, state_here, eo_here, 1, 1, provenance, 1.0)
            call load_weighted(op_oracle, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
            call op_oracle%scale_raw_accum(1.0/f1)
            call op_replay%compare_raw_accum(op_oracle, berr, derr)
            call require_raw('bootstrap full-mass chain seed', berr, derr, REPLAY_TOL)
            call op_replay%kill
            call op_oracle%kill

            call new_reduction(op_sum)
            call op_sum%add_raw_accum(f_sub1, state_here, eo_here, 1, 2, provenance, nraw)
            nraw_total = nraw
            call op_sum%add_raw_accum(f_sub2, state_here, eo_here, 2, 2, provenance, nraw)
            nraw_total = nraw_total + nraw
            if( nraw_total /= size(full_pinds) ) THROW_HARD('complementary raw PCG particle counts do not close')
            call load_weighted(op_ref, f_full, state_here, eo_here, 1, 1, provenance, 1.0)
            call op_sum%compare_raw_accum(op_ref, berr, derr)
            call require_raw('complementary additivity', berr, derr, RAW_TOL)
            call op_sum%kill
            call op_ref%kill

            call build_blend(f_sub1, 1, f1, f_full, state_here, eo_here, provenance, u, op_blend)
            call load_weighted(op_oracle, f_sub1, state_here, eo_here, 1, 2, provenance, 1.0)
            call op_oracle%scale_raw_accum(u/f1)
            call add_weighted(op_oracle, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
            call op_blend%compare_raw_accum(op_oracle, berr, derr)
            call require_raw('subset1 u/f blend', berr, derr, REPLAY_TOL)
            call op_blend%write_raw_accum(f_chain1, state_here, eo_here, 1, 1, size(full_pinds), provenance)
            call op_oracle%kill

            call load_weighted(op_replay, f_chain1, state_here, eo_here, 1, 1, provenance, 1.0)
            call op_blend%compare_raw_accum(op_replay, berr, derr)
            call require_raw('continuation write/read', berr, derr, REPLAY_TOL)
            call finalize_and_solve(op_blend, x_direct, hist, niters)
            deallocate(hist)
            call finalize_and_solve(op_replay, x_replay, hist, niters)
            deallocate(hist)
            solve_err = volume_rel_error(x_direct, x_replay)
            if( solve_err > REPLAY_TOL ) THROW_HARD('PCG continuation replay changed the reconstructed solution')
            call op_blend%kill
            call op_replay%kill
            deallocate(x_replay)

            call build_blend(f_sub2, 2, f2, f_full, state_here, eo_here, provenance, u, op_blend)
            call load_weighted(op_oracle, f_sub2, state_here, eo_here, 2, 2, provenance, 1.0)
            call op_oracle%scale_raw_accum(u/f2)
            call add_weighted(op_oracle, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
            call op_blend%compare_raw_accum(op_oracle, berr, derr)
            call require_raw('subset2 u/f blend', berr, derr, REPLAY_TOL)
            call op_blend%write_raw_accum(f_chain2, state_here, eo_here, 1, 1, size(full_pinds), provenance)
            call op_blend%kill
            call op_oracle%kill

            call load_weighted(op_ensemble, f_chain1, state_here, eo_here, 1, 1, provenance, f1)
            call add_weighted(op_ensemble, f_chain2, state_here, eo_here, 1, 1, provenance, f2)
            call load_weighted(op_ref, f_full, state_here, eo_here, 1, 1, provenance, 1.0)
            call op_ensemble%compare_raw_accum(op_ref, berr, derr)
            call require_raw('full-mass weighted ensemble', berr, derr, RAW_TOL)
            call op_ensemble%kill
            call finalize_and_solve(op_ref, x_full, hist, niters)
            deallocate(hist)
            sample_map_err = volume_rel_error(x_direct, x_full)
            call op_ref%kill

            write(unit,'(A,ES14.6)') 'continuation_replay_solution_relerr=', solve_err
            write(unit,'(A,ES14.6)') 'single_subset_vs_full_map_relerr_diagnostic_only=', sample_map_err
            write(unit,'(A)') 'status=PASS'
            write(logfhandle,'(A,I0,A,A,A,F7.4,A,ES10.3)') '>>> PCG FRAC | STATE=', state_here, &
                &' | HALF=', trim(half_here), ' | U=', u, ' | REPLAY=', solve_err

            deallocate(x_direct, x_full)
            call del_file(f_full)
            call del_file(f_sub1)
            call del_file(f_sub2)
            call del_file(f_chain1)
            call del_file(f_chain2)
            call f_full%kill
            call f_sub1%kill
            call f_sub2%kill
            call f_chain1%kill
            call f_chain2%kill
        end subroutine validate_half

        subroutine accumulate_raw( pinds, op )
            integer,                 intent(in)    :: pinds(:)
            type(reconstructor_pcg), intent(inout) :: op
            type(oris) :: selection
            type(ori) :: orientation
            type(ctfparams) :: ctfparms
            type(image) :: obs
            complex, allocatable :: y_batch(:,:,:)
            real, allocatable :: sig2(:,:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch
            real :: shift(2), crop_factor
            call op%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call op%set_sym(build%pgrpsyms)
            call op%set_mask(params%msk_crop)
            lims2 = op%get_lims2()
            R = lims2(1,2)
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
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                ctfparms%smpd = params%smpd_crop
                shift = build%spproj_field%get_2Dshift(iptcl) * crop_factor
                call orientation%set_ctfvars(ctfparms)
                call orientation%set_shift(shift)
                call selection%set_ori(i, orientation)
            enddo
            call op%prep_particles(selection, use_ctf=.true., sig2=sig2)
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call op%begin_accum
            call obs%new([params%box_crop,params%box_crop,1], params%smpd_crop)
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz = batchlims(2)-batchlims(1)+1
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
                do ii = 1, batchsz
                    ! the backend-neutral observation (normalize, crop, taper), see prep_rec_observation
                    call prep_rec_observation(build%imgbatch(ii), build%lmsk, obs, .true.)
                    call obs%fft
                    y_batch(:,:,ii) = op%extract_native_plane(obs)
                enddo
                call op%accumulate_batch(y_batch, batchsz, batchlims(1))
            enddo
            call obs%kill
            call selection%kill
            call orientation%kill
            deallocate(y_batch, sig2)
        end subroutine accumulate_raw

        subroutine new_reduction( op )
            type(reconstructor_pcg), intent(inout) :: op
            call op%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call op%set_mask(params%msk_crop)
            call op%begin_reduction
        end subroutine new_reduction

        subroutine load_weighted( op, fname, state_here, eo_here, part, nparts, provenance, weight )
            type(reconstructor_pcg), intent(inout) :: op
            class(string),           intent(in)    :: fname
            integer,                 intent(in)    :: state_here, eo_here, part, nparts
            character(len=*),        intent(in)    :: provenance
            real,                    intent(in)    :: weight
            call new_reduction(op)
            call add_weighted(op, fname, state_here, eo_here, part, nparts, provenance, weight)
        end subroutine load_weighted

        subroutine add_weighted( op, fname, state_here, eo_here, part, nparts, provenance, weight )
            type(reconstructor_pcg), intent(inout) :: op
            class(string),           intent(in)    :: fname
            integer,                 intent(in)    :: state_here, eo_here, part, nparts
            character(len=*),        intent(in)    :: provenance
            real,                    intent(in)    :: weight
            integer :: nraw
            call op%add_raw_accum_weighted(fname, state_here, eo_here, part, nparts, provenance, weight, nraw)
        end subroutine add_weighted

        subroutine build_blend( f_subset, subset_part, frac, f_full, state_here, eo_here, provenance, u, op )
            class(string),           intent(in)    :: f_subset, f_full
            integer,                 intent(in)    :: subset_part, state_here, eo_here
            real,                    intent(in)    :: frac, u
            character(len=*),        intent(in)    :: provenance
            type(reconstructor_pcg), intent(inout) :: op
            call load_weighted(op, f_subset, state_here, eo_here, subset_part, 2, provenance, u/frac)
            call add_weighted(op, f_full, state_here, eo_here, 1, 1, provenance, 1.0-u)
        end subroutine build_blend

        subroutine finalize_and_solve( op, x, history, niters )
            type(reconstructor_pcg), intent(inout) :: op
            real, allocatable,       intent(out)   :: x(:,:,:), history(:)
            integer,                 intent(out)   :: niters
            type(pcg_solver_outcome) :: outcome
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            call op%end_accum(.true.)
            call op%set_op_mode(PCG_OP_KERNEL)
            call op%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, rel_res_hist=history, &
                &niters=niters, outcome=outcome)
            if( trim(outcome%stop_reason) == PCG_STOP_INDEFINITE ) &
                &THROW_HARD('PCG lost positive-definiteness in the harness solve')
        end subroutine finalize_and_solve

        subroutine require_raw( label, b_err, d_err, tolerance )
            character(len=*), intent(in) :: label
            real,             intent(in) :: b_err, d_err, tolerance
            write(funit,'(A,A)') 'gate=', trim(label)
            write(funit,'(A,ES14.6)') 'B_relerr=', b_err
            write(funit,'(A,ES14.6)') 'D_relerr=', d_err
            if( b_err > tolerance .or. d_err > tolerance )then
                write(logfhandle,'(A,A,A,ES12.4,A,ES12.4)') '>>> PCG FRAC FAIL | ', trim(label), &
                    &' | B=', b_err, ' | D=', d_err
                THROW_HARD('PCG fractional-update raw accumulator validation failed')
            endif
        end subroutine require_raw

        pure real function volume_rel_error( lhs, rhs ) result(err)
            real, intent(in) :: lhs(:,:,:), rhs(:,:,:)
            real(dp) :: numerator, denominator
            integer :: i, j, k
            numerator = 0.0_dp
            denominator = 0.0_dp
            do k = 1, size(lhs,3)
                do j = 1, size(lhs,2)
                    do i = 1, size(lhs,1)
                        numerator = numerator + real(lhs(i,j,k)-rhs(i,j,k),dp)**2
                        denominator = denominator + real(rhs(i,j,k),dp)**2
                    enddo
                enddo
            enddo
            err = real(sqrt(numerator) / max(1.0_dp, sqrt(denominator)))
        end function volume_rel_error

        function frac_fname( state_here, half_here, tag ) result(fname)
            integer,          intent(in) :: state_here
            character(len=*), intent(in) :: half_here, tag
            type(string) :: fname
            fname = 'pcg_frac_state'//int2str_pad(state_here,2)//'_'//trim(half_here)//'_'//trim(tag)//'.raw'
        end function frac_fname

    end subroutine validate_rec3D_pcg_fractional_updates

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
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
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
            real, allocatable :: x(:,:,:), rel_res_hist(:)
            integer :: state = 0, eo = 0, nptcls = 0, niters = 0
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
        type(image), target  :: half_even, half_odd, ml_even, ml_odd, merged
        type(image), target  :: previous_even, previous_odd, previous_merged
        type(image), pointer :: fsc_pair_even, fsc_pair_odd, fsc_pair_merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt, eonames(2)
        type(halfmap_diagnostics_result) :: hm_diag
        real, allocatable :: fsc(:), res0143s(:), res05s(:), align_lps(:)
        real, allocatable :: realized_fractions(:), update_weights(:)
        logical, allocatable :: state_written(:)
        character(len=256) :: provenance, chain_provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        integer :: n_active_state, n_sampled_state
        integer :: pcg_master_nthreads, pcg_half_nthreads
        type(distributed_half_job) :: even_job, odd_job
        type(image_msk) :: state_support_msk
        logical :: l_state_support, l_base_support_constrained, l_nu_base_constrained
        logical :: l_has_updates, l_bootstrap, l_even_chain, l_odd_chain
        logical :: l_fsc_pair_support_constrained, l_prev_support_constrained, l_prev_provenance_found
        logical :: l_shipped_support_constrained
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output, time_nu_filter

        call validate_pcg_common(params)
        ! the partition workers are idle during the master-side solve and NU
        ! filtering phases: on local execution use the full
        ! allocation, restored to nthr before returning to the matching
        ! phase (never on a cluster, where the master's slot is fixed).
        ! Capped: OpenMP scaling saturates well before large core counts
        ! at these box sizes
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
        if( params%l_trail_rec )then
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
        allocate(state_written(params%nstates), source=.false.)
        do state = 1, params%nstates
            l_bootstrap = .false.
            if( params%l_trail_rec )then
                ! constant-FOV crop growth keeps the chain continuous (the
                ! reader embeds a smaller previous grid by zero-extension);
                ! only a true identity change -- provenance, field of view,
                ! a shrinking box, or a pre-v2 chain format -- discards the
                ! pair and re-seeds through the bootstrap blend, mirroring
                ! the gridding backend's manifest validation policy
                call discard_stale_trail_chain_pair(state)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'even')
                l_even_chain = file_exists(raw_fname)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'odd')
                l_odd_chain = file_exists(raw_fname)
                if( l_even_chain .neqv. l_odd_chain ) THROW_HARD('PCG trailing chain pair is incomplete')
                l_bootstrap = .not. l_even_chain
            endif
            if( present(trail_bootstrap_states) ) trail_bootstrap_states(state) = l_bootstrap
            ! Build one conservative support per state, only under automsk=yes.
            ! It is installed in both solves once a prior reconstruction
            ! exists (envfsc=yes follows); otherwise the base bootstraps on
            ! the sphere and the replay uses density. automsk=no: sphere throughout.
            call build_pcg_state_support(params, state, state_support_msk, l_state_support)
            l_base_support_constrained = l_state_support
            call reduce_solve_state_pair(state, half_even, half_odd, n_even, n_odd, 'base')
            if( params%l_trail_rec )then
                call count_state_sampling(state, n_active_state, n_sampled_state)
                if( n_even+n_odd /= n_sampled_state ) THROW_HARD('PCG raw particles do not match the latest sampled cohort')
                if( n_active_state > 0 )then
                    if( abs(real(n_sampled_state)/real(n_active_state)-realized_fractions(state)) > 1.0e-6 )then
                        THROW_HARD('PCG realized fraction disagrees with gridding sampling bookkeeping')
                    endif
                endif
            endif
            if( n_even == 0 .and. n_odd == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG DISTRIBUTED: STATE ', state, &
                    &' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even < 1 .or. n_odd < 1 ) THROW_HARD('distributed PCG requires both halfsets')
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
            t_state_phase = tic()
            ! the FSC pair is selected once here for the FSC and its summary
            ! only: the current base pair ordinarily (in trailing mode the
            ! blended chain solution the replay re-reads), the previous shipped
            ! pair in the trailing bootstrap (lag-one FSC prior, exactly as the
            ! gridding bootstrap). The NU candidate bank is ALWAYS seeded from
            ! the current base pair (volume-blended with the previous pair in
            ! the bootstrap), never from the lag-one FSC pair -- same recipe
            ! as gridding, review 2026-09-06.
            if( l_bootstrap )then
                call load_previous_state_halves(state, previous_even, previous_odd, previous_merged)
                fsc_pair_even        => previous_even
                fsc_pair_odd         => previous_odd
                fsc_pair_merged      => previous_merged
                ! Current support availability says nothing about how the
                ! lagged previous pair was reconstructed: read the solve-support
                ! provenance persisted beside it. An imported pair without a
                ! sidecar is treated as unconstrained so an envfsc request
                ! receives the phase-randomized correction.
                call read_support_provenance(params%vols(state), l_prev_support_constrained, &
                    &l_prev_provenance_found)
                l_fsc_pair_support_constrained = l_prev_provenance_found .and. l_prev_support_constrained
                if( .not. l_prev_provenance_found ) write(logfhandle,'(A,I0,A)') &
                    &'>>> PCG DISTRIBUTED: STATE ', state, &
                    &' previous pair has no solve-support provenance; treated as unconstrained'
            else
                fsc_pair_even        => half_even
                fsc_pair_odd         => half_odd
                fsc_pair_merged      => merged
                l_fsc_pair_support_constrained = l_base_support_constrained
            endif
            call calculate_pcg_state_diagnostics(params, state, 'DISTRIBUTED', fsc_pair_even, &
                &fsc_pair_odd, fsc_pair_merged, hm_diag, l_fsc_pair_support_constrained)
            fsc             = hm_diag%fsc
            res0143s(state) = hm_diag%res_fsc0143
            res05s(state)   = hm_diag%res_fsc05
            call arr2file(fsc, fname_fsc)
            fname_restxt = resolve_pcg_fsc_txt_fname(params, cline, state)
            call write_halfmap_diagnostics(hm_diag, params%box_crop, params%smpd_crop, fname_restxt)
            call fname_restxt%kill
            call hm_diag%kill
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                ! the ordinary global-ML replay (P_tau from the FSC pair)
                call reduce_solve_state_pair(state, ml_even, ml_odd, n_even, n_odd, 'ml', fsc, &
                    &half_even, half_odd)
                call merged%kill
                call merged%copy(ml_even)
                call merged%add(ml_odd)
                call merged%mul(0.5)
            endif
            if( l_bootstrap .and. update_weights(state) < 0.99 )then
                ! legacy volume-domain blend of the bootstrap iteration, the
                ! same recipe as the gridding trail_restored_halves_if_needed
                ! applied to this solver's own maps: the previous shipped pair
                ! is scaled once and added to BOTH the base (NU bank source)
                ! and the regularized (shipped, NU aux) pair. Same recipe and
                ! the same kinds of inputs on both backends; the maps differ
                ! because the estimators differ
                call previous_even%mul(1.0-update_weights(state))
                call previous_odd%mul( 1.0-update_weights(state))
                call blend_bootstrap_half(half_even, previous_even, update_weights(state))
                call blend_bootstrap_half(half_odd,  previous_odd,  update_weights(state))
                if( params%l_ml_reg )then
                    call blend_bootstrap_half(ml_even, previous_even, update_weights(state))
                    call blend_bootstrap_half(ml_odd,  previous_odd,  update_weights(state))
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
            ! a bootstrap blend carries the previous pair's support into the
            ! shipped pair: constrained only if both contributions were
            if( l_bootstrap .and. update_weights(state) < 0.99 ) &
                &l_shipped_support_constrained = l_shipped_support_constrained .and. l_fsc_pair_support_constrained
            if( params%l_ml_reg )then
                call ml_even%write(fname_even, del_if_exists=.true.)
                call ml_odd%write(fname_odd, del_if_exists=.true.)
            else
                call half_even%write(fname_even, del_if_exists=.true.)
                call half_odd%write(fname_odd, del_if_exists=.true.)
            endif
            call merged%write(fname_vol, del_if_exists=.true.)
            ! the sidecar follows the published map, never precedes it
            if( params%l_ml_reg )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'regularized')
            else if( l_bootstrap .and. update_weights(state) < 0.99 )then
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'mixed')
            else
                call write_support_provenance(fname_vol, l_shipped_support_constrained, 'base')
            endif
            time_map_output = time_map_output + real(toc(t_state_phase),dp)
            if( params%l_nonuniform )then
                ! NU competition, mirroring the gridding volassemble: the
                ! candidate bank from the current base pair (bootstrap: blended
                ! with the previous pair above, never the lag-one FSC pair),
                ! the shipped ML pair as the auxiliary member; consumes both
                t_state_phase = tic()
                eonames(1) = fname_even
                eonames(2) = fname_odd
                ! an envelope-constrained base pair hands its support over so
                ! the evidence null is designated on the dilation ring; a
                ! bootstrap blend counts as constrained only if both
                ! contributions were
                l_nu_base_constrained = l_base_support_constrained
                if( l_bootstrap .and. update_weights(state) < 0.99 ) &
                    &l_nu_base_constrained = l_nu_base_constrained .and. l_fsc_pair_support_constrained
                call nonuniform_filter_state(params, state, half_even, half_odd, &
                    &ml_even, ml_odd, params%l_ml_reg .and. nu_static_aux_replacement(params), &
                    &res0143s(state), fname_vol, eonames, align_lps(state), &
                    &base_support=state_support_msk, l_base_constrained=l_nu_base_constrained)
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
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) )then
                        call build%spproj_field%set(iptcl, 'res',   res0143s(istate))
                        call build%spproj_field%set(iptcl, 'res05', res05s(istate))
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
        deallocate(res0143s, res05s, state_written, realized_fractions, update_weights, align_lps)
        !$ call omp_set_num_threads(params%nthr)

    contains

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

        !> Discard a trailing chain pair that is genuinely unusable: a
        !! continuation-identity (provenance) change, a field-of-view change,
        !! a chain persisted at a LARGER crop than the current one, or an
        !! unreadable/old-format artifact. Constant-FOV crop growth (the
        !! abinitio3D stage ramp) is NOT stale -- the weighted reader embeds
        !! the smaller previous grid by index-aligned zero-extension. Both
        !! halves are removed together so a stale half can never combine with
        !! a fresh one; the caller then re-enters the bootstrap blend.
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

        subroutine reduce_solve_state_pair( state_here, even, odd, n_even_here, n_odd_here, solve_kind, &
                &fsc_prior, warm_even, warm_odd )
            integer,          intent(in)    :: state_here
            character(len=*), intent(in)    :: solve_kind
            type(image),      intent(inout) :: even, odd
            integer,          intent(out)   :: n_even_here, n_odd_here
            real, optional,   intent(in)    :: fsc_prior(:)
            type(image), optional, intent(in) :: warm_even, warm_odd

            if( present(fsc_prior) )then
                if( .not. present(warm_even) .or. .not. present(warm_odd) ) &
                    &THROW_HARD('distributed PCG ML replay requires both half-map warm starts')
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

            if( present(fsc_prior) )then
                call finish_distributed_half_job(even_job, even, warm_even)
                call finish_distributed_half_job(odd_job, odd, warm_odd)
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
            job%niters = 0
            job%ready = .false.
            job%l_ml_solve = present(fsc_prior)
            if( job%l_ml_solve .neqv. present(warm_start) ) &
                &THROW_HARD('distributed PCG ML replay requires both FSC and warm start')
            if( allocated(job%x) ) deallocate(job%x)
            if( allocated(job%rel_res_hist) ) deallocate(job%rel_res_hist)

            call job%pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA, &
                &fft_nthreads=pcg_half_nthreads)
            ! Under automsk=yes both the base and the regularized pass take
            ! the density envelope once it exists; with automsk=no
            ! l_state_support is false and both run on the sphere. This call
            ! stays outside the parallel region because spherical-mask
            ! construction memoizes coordinates at module scope.
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
            if( job%nptcls == 0 )then
                call job%pcgop%kill
                return
            endif
            if( .not. job%l_ml_solve )then
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
                        current_scale = update_weight / realized_fraction
                        call job%pcgop%scale_raw_accum(current_scale)
                        call job%pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                            &chain_provenance, 1.0-update_weight, n_part)
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
            endif

            t_phase = tic()
            if( job%l_ml_solve ) call job%pcgop%set_ml_prior(fsc_prior, params%tau, params%hp)
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
                ! the replay starts from the current base solution with the
                ! shrinkage initial guess; no cross-iteration warm start
                job%x = warm_start%get_rmat()
                call regularized_ml_initial_guess(params, fsc_prior, job%x, 'distributed', half)
                ! the replay iterate is never zero: restart-eligible
                job%l_nonzero = .true.
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
            call solve_with_cold_restart(job%pcgop, job%x, job%l_nonzero, params%maxits_pcg, params%rtol, &
                &job%rel_res_hist, job%niters, job%result)
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
            call report_solve_summary('DISTRIBUTED', job%state, job%half, job%solve_kind, job%nptcls, &
                &job%niters, job%result%final_rel_residual, job%time_solve, job%result%stop_reason, &
                &job%result%initial_rel_residual)
            call job%pcgop%kill
            if( allocated(job%x) ) deallocate(job%x)
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
            write(funit,'(A,ES14.6)') 'final_rel_update=',      result%final_rel_update
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

    subroutine validate_solved_map( x, execution_mode, state, half, solve_kind )
        real,             intent(in) :: x(:,:,:)
        character(len=*), intent(in) :: execution_mode, half, solve_kind
        integer,          intent(in) :: state
        character(len=256) :: error_message
        real :: peak
        if( any(.not. ieee_is_finite(x)) )then
            error_message = 'PCG '//trim(execution_mode)//' solve produced a non-finite map; state='// &
                &int2str(state)//' half='//trim(half)//' kind='//trim(solve_kind)
            THROW_HARD(error_message)
        endif
        peak = maxval(abs(x))
        if( peak <= 0.0 )then
            error_message = 'PCG '//trim(execution_mode)//' solve produced an empty map; state='// &
                &int2str(state)//' half='//trim(half)//' kind='//trim(solve_kind)
            THROW_HARD(error_message)
        endif
    end subroutine validate_solved_map

    !> Amplitude convention: maps are stored at the particle coefficient scale
    !! (the plain data quotient, the solver's native convention), identical to
    !! the gridding backend since the matched box-size pair (reconstructor
    !! division, projector multiplication) was retired; see
    !! doc/implementation_notes/drop_legacy_box_division.md S0/S5.3a.


    !> The current matching band in native crop-box shells (params%kfromto(2)),
    !! or 0 when it does not describe a usable band of this volume.
    integer function matched_band_kstop( params, volume ) result( kstop )
        type(parameters), intent(in) :: params
        type(image),      intent(in) :: volume
        kstop = params%kfromto(2)
        if( kstop < 1 .or. kstop >= volume%get_filtsz() ) kstop = 0
    end function matched_band_kstop

    !> Diagnostic only (the map is not modified): compare the RMS amplitude of
    !! the shells beyond the matching band with that of the band-edge shell and
    !! report when the former dominates. Reconstructing beyond the matching
    !! band is essential (resolution extending past the matching limit is the
    !! validation that the map is right), but the fixed-iteration PCG solve has
    !! been observed to leave beyond-band content orders of magnitude above the
    !! band edge under bootstrap-scale sigma2 (see
    !! doc/implementation_notes/pcg_priors_history.md §2, beyond-band excess). Such
    !! content is invisible to the matcher until a stage transition extends the
    !! band into it, so this line is the regression signal for that solver
    !! defect. Silent when no matching band is known.
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

    !> One summary line per solve; INIT is the relative residual of the start
    !! (1.0 from zero) so a start that is worse than nothing is visible in the
    !! log (the PfCRT regression hid a 10^4 residual in the sidecar files)
    subroutine report_solve_summary( execution_mode, state, half, solve_kind, nptcls, niters, &
            &residual, solve_time, stop_reason, initial_residual )
        character(len=*), intent(in) :: execution_mode, half, solve_kind, stop_reason
        integer,          intent(in) :: state, nptcls, niters
        real,             intent(in) :: residual
        real(dp),         intent(in) :: solve_time
        real,             intent(in) :: initial_residual
        character(len=4) :: half_label, kind_label
        half_label = adjustl(half)
        kind_label = adjustl(solve_kind)
        write(logfhandle,'(4A,I2,A,A4,A,A4,A,I6,A,I2,A,ES10.3,A,ES10.3,A,F7.2,2A)') &
            &'>>> PCG ', trim(execution_mode), ' | ', 'STATE=', state, ' | HALF=', half_label, &
            &' | KIND=', kind_label, ' | N=', nptcls, ' | ITS=', niters, ' | INIT=', initial_residual, &
            &' | RESID=', residual, ' | TIME=', solve_time, ' s | STOP=', trim(stop_reason)
        call flush(logfhandle)
    end subroutine report_solve_summary

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
            if( params%maxits_pcg > 8 )then
                THROW_WARN('maxits_pcg exceeds the production refinement budget (8); appropriate for offline converged solves only')
            endif
        endif
        if( trim(params%projrec) /= 'no' ) THROW_HARD('rec_backend=pcg does not yet support projrec=yes')
        if( abs(real(params%box)*params%smpd - real(params%box_crop)*params%smpd_crop) > &
            &1.0e-5*real(params%box)*params%smpd )then
            THROW_HARD('PCG crop must preserve the native physical box extent')
        endif
        if( trim(params%conical_fsc) == 'yes' ) THROW_HARD('PCG conical FSC integration is not implemented')
        if( params%msk <= 0.5 .or. params%msk_crop <= 0.5 ) THROW_HARD('rec_backend=pcg requires mskdiam')
        if( .not. ieee_is_finite(params%rtol) ) THROW_HARD('PCG rtol must be finite')
    end subroutine validate_pcg_common

    function pcg_raw_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgraw-v2|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|ptcl_src='//trim(params%ptcl_src)//'|iter='//trim(int2str(params%which_iter))// &
            &'|box='//trim(int2str(params%box))//'|smpd='//trim(real2str(params%smpd))// &
            &'|box_crop='//trim(int2str(params%box_crop))// &
            &'|smpd_crop='//trim(real2str(params%smpd_crop))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_raw_provenance

    ! Continuation identity deliberately excludes which_iter AND the crop
    ! geometry: the chain must survive iteration boundaries and constant-FOV
    ! crop changes (stage transitions; the reader embeds a smaller previous
    ! grid by zero-extension and validates the field of view structurally),
    ! while native geometry, objective and particle source changes must
    ! invalidate it.
    function pcg_chain_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgtrail-v2|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|ptcl_src='//trim(params%ptcl_src)//'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_chain_provenance

end module simple_rec3D_pcg_strategy
