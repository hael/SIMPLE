!@descr: shared-memory production strategy body for kernel PCG reconstruct3D
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, &
    &discrete_read_imgbatch_source, killimgbatch
use simple_ptcl_cache,        only: ptcl_cache_read_batch
use simple_sigma2_files,      only: load_sigma2_groups
use simple_math_ft,           only: resample_sigma2
use simple_fsc,               only: phase_rand_fsc, fsc_area_score_result
use simple_estimate_ssnr,     only: fsc2shrink_filter
use simple_image,             only: image
use simple_image_msk,         only: image_msk
use simple_refine3D_fnames,   only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname, &
    &refine3D_pcg_trail_accum_fname
implicit none

public :: execute_rec3D_pcg_shared, execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
public :: validate_rec3D_pcg_fractional_updates
private
#include "simple_local_flags.inc"

real,    parameter :: PCG_LAMBDA = 1.0e-3
logical, parameter :: DEBUG = .false.

contains

    !> Cross-iteration ML warm start (drop_legacy_box_division.md S7/S11.2).
    !! The ML replay used to start from the unregularized base solution, which
    !! carries full-amplitude beyond-band noise the ML prior drives toward
    !! zero; on large boxes a 2-4 iteration budget cannot close that gap
    !! (bgal, box 256: relative residual 6.5 after 2 iterations). The previous
    !! refinement iteration's ML half map is close to the current ML solution
    !! (small per-iteration orientation change), and under the data-quotient
    !! convention a constant-FOV box change is a factor-free Fourier pad/clip,
    !! so read_and_crop makes the warm start valid across crop changes. Rules:
    !! each half warm-starts strictly from its own previous half (gold-standard
    !! FSC independence), the soft support mask is re-applied after resampling
    !! (Fourier padding rings slightly outside it), and the noise starting
    !! volume of the first iteration is excluded by its workflow-contract name
    !! (warm-starting the ML system from noise is worse than the base
    !! solution). When no usable previous half exists, x keeps the base
    !! solution it already holds.
    subroutine override_ml_warm_start_from_previous( params, state_here, half, x, context, l_found )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        character(len=*),  intent(in)    :: half, context
        real,              intent(inout) :: x(:,:,:)
        logical,           intent(out)   :: l_found
        type(string) :: prev_fname
        type(image)  :: prev
        l_found = .false.
        if( state_here < 1 .or. state_here > size(params%vols) ) return
        if( len_trim(params%vols(state_here)%to_char()) == 0 ) return
        prev_fname = add2fbody(params%vols(state_here), MRC_EXT, '_'//trim(half))
        if( index(prev_fname%to_char(), 'startvol') > 0 )then
            call prev_fname%kill
            return
        endif
        if( .not. file_exists(prev_fname) )then
            call prev_fname%kill
            return
        endif
        call prev%read_and_crop(prev_fname, params%smpd, params%box_crop, params%smpd_crop)
        call prev%mask3D_soft(params%msk_crop, backgr=0.)
        x = prev%get_rmat()
        call prev%kill
        l_found = .true.
        write(logfhandle,'(A)') '>>> PCG ML WARM START ('//trim(context)//'/'//trim(half)//&
            &'): previous-iteration ML half map '//prev_fname%to_char()
        call prev_fname%kill
    end subroutine override_ml_warm_start_from_previous

    !> Regularized initial guess for a COLD ML replay (no previous ML half
    !! map to warm-start from: the standalone harness, the first refinement
    !! iteration, abinitio3D stage handoffs). The documented cold-start gap
    !! is that the regularized optimum differs from ANY unregularized map in
    !! slowly-converging directions -- beyond-band/high-shell noise the ML
    !! prior shrinks, solvent texture the solvent prior flattens -- so a
    !! small iteration budget leaves the solve transient-dominated. Apply
    !! the regularizers' expected effect to the base solution in closed form
    !! instead: shrink each shell's amplitude by the FSC (the Wiener
    !! shrinkage the ML prior's optimum implies; no shrinkage below the hp
    !! no-prior limit, zero beyond the measured band), then, when a
    !! validated solvent envelope is attached, blend the solvent domain
    !! toward its weighted mean -- Q_s's fixed point. CG then corrects an
    !! approximate optimum rather than constructing it. Pure initialization:
    !! the quadratic objective has a unique optimum, so this changes the
    !! convergence path, never the converged solution.
    subroutine regularized_ml_initial_guess( params, fsc, m_env, l_env, x, context, half )
        class(parameters), intent(in)    :: params
        real,              intent(in)    :: fsc(:)
        real, allocatable, intent(in)    :: m_env(:,:,:)
        logical,           intent(in)    :: l_env
        real,              intent(inout) :: x(:,:,:)
        character(len=*),  intent(in)    :: context, half
        type(image) :: img
        real, allocatable :: filt(:)
        real(dp) :: wsum, xsum
        real     :: mu
        call img%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
        call img%set_rmat(x, .false.)
        call ml_shrinkage_filter(params, fsc, img%get_filtsz(), filt)
        call img%fft()
        call img%apply_filter(filt)
        call img%ifft()
        x = img%get_rmat()
        call img%kill
        deallocate(filt)
        if( l_env )then
            wsum = sum(1.0_dp - real(m_env,dp))
            if( wsum > 0.0_dp )then
                xsum = sum((1.0_dp - real(m_env,dp)) * real(x,dp))
                mu   = real(xsum / wsum)
                x    = m_env * x + (1.0 - m_env) * mu
            endif
            write(logfhandle,'(A)') '>>> PCG ML REGULARIZED INIT ('//trim(context)//'/'//trim(half)//&
                &'): shell-shrunk base + envelope-flattened solvent'
        else
            write(logfhandle,'(A)') '>>> PCG ML REGULARIZED INIT ('//trim(context)//'/'//trim(half)//&
                &'): shell-shrunk base'
        endif
    end subroutine regularized_ml_initial_guess

    !> The strategy-side wrapper for the ML shrinkage filter (fsc2shrink_filter
    !! in simple_estimate_ssnr, alongside the other FSC-derived filters):
    !! resolves the no-shrinkage high-pass index from the workflow parameters.
    !! Shared by the regularized initial guess and the solvent suppression
    !! reference.
    subroutine ml_shrinkage_filter( params, fsc, nyq, filt )
        class(parameters), intent(in)  :: params
        real,              intent(in)  :: fsc(:)
        integer,           intent(in)  :: nyq
        real, allocatable, intent(out) :: filt(:)
        integer :: k_hp
        k_hp = max(1, calc_fourier_index(params%hp, params%box_crop, params%smpd_crop))
        call fsc2shrink_filter(fsc, k_hp, nyq, filt)
    end subroutine ml_shrinkage_filter

    !> Reference level for the solvent suppression readout: the weighted
    !! solvent RMS of the shell-shrunk base solution. The FSC/Wiener shrinkage
    !! is the part of the ML replay that reduces solvent variation WITHOUT the
    !! solvent prior, so measuring the final ML map against this reference
    !! isolates the prior's own contribution -- an inert prior reads ~0%
    !! suppression regardless of tau.
    real function solvent_suppression_reference( pcgop, params, fsc, base_volume ) result( rms_ref )
        type(reconstructor_pcg), intent(in) :: pcgop
        class(parameters),       intent(in) :: params
        real,                    intent(in) :: fsc(:)
        type(image),             intent(in) :: base_volume
        type(image) :: img
        real, allocatable :: filt(:), xr(:,:,:)
        real :: mean_s, penalty
        call img%copy(base_volume)
        call ml_shrinkage_filter(params, fsc, img%get_filtsz(), filt)
        call img%fft()
        call img%apply_filter(filt)
        call img%ifft()
        xr = img%get_rmat()
        call pcgop%get_solvent_stats(xr, mean_s, rms_ref, penalty)
        call img%kill
        deallocate(filt, xr)
    end function solvent_suppression_reference

    !> Persist the per-state solvent suppression readout for the convergence
    !! reporter (simple_convergence prints it with the other iteration stats
    !! and advises on pcg_solvent_lambda_rel). The file is rewritten on every
    !! ML volassemble and deleted first, so an iteration in which the prior is
    !! skipped never leaves stale values behind.
    subroutine write_solvent_convergence_stats( params, supps, rmss, cnts )
        class(parameters), intent(in) :: params
        real,              intent(in) :: supps(:), rmss(:)
        integer,           intent(in) :: cnts(:)
        type(oris)   :: os
        type(string) :: key
        integer :: state
        call del_file(PCG_SOLVENT_STATS_FILE)
        if( .not. any(cnts > 0) ) return
        call os%new(1, is_ptcl=.false.)
        call os%set(1, 'PCG_SOLVENT_LAMBDA_REL', params%pcg_solvent_lambda_rel)
        do state = 1, size(supps)
            if( cnts(state) < 1 ) cycle
            key = 'PCG_SOLVENT_SUPP_STATE'//int2str_pad(state,2)
            call os%set(1, key%to_char(), supps(state) / real(cnts(state)))
            key = 'PCG_SOLVENT_RMS_STATE'//int2str_pad(state,2)
            call os%set(1, key%to_char(), rmss(state) / real(cnts(state)))
        end do
        call os%write(string(PCG_SOLVENT_STATS_FILE))
        call os%kill
        call key%kill
    end subroutine write_solvent_convergence_stats

    !> Solvent-envelope loader for the graded solvent-flatness prior
    !! (pcg_priors.md S4, S6.2). Resolves the state-specific NU evidence
    !! envelope, runs the full S6.2 validation contract -- presence, cubic
    !! lattice, physical-extent identity, constant-FOV resample when the
    !! lattice differs, finiteness, [0,1] range with re-clip after resampling,
    !! nonzero solvent evidence -- and emits the pcg_solvent_* diagnostic
    !! block. When the envelope validates AND pcg_solvent_lambda_rel > 0 the
    !! clipped envelope is returned with l_ready true, and the caller attaches
    !! it to BOTH halves' ML replays (one fixed mask identity per state and
    !! iteration, S6.2). In every other case l_ready is false with a clear
    !! skip reason and the solve path is untouched (R5: strength zero is
    !! bit-identical). Never substitutes a density or spherical mask.
    subroutine resolve_solvent_envelope( params, state_here, context, m_env, l_ready )
        class(parameters), intent(in)  :: params
        integer,           intent(in)  :: state_here
        character(len=*),  intent(in)  :: context
        real, allocatable, intent(out) :: m_env(:,:,:)
        logical,           intent(out) :: l_ready
        real, parameter :: EXTENT_RELTOL = 1.0e-3 !< constant-FOV identity tolerance
        real, parameter :: VALUE_SLACK   = 1.0e-2 !< clip-vs-reject boundary around [0,1]
        type(string) :: env_fname
        type(image)  :: env
        real,    allocatable :: m(:,:,:)
        character(len=:), allocatable :: skip_reason
        real     :: smpd_env, extent_env, extent_rec, mask_min, mask_max
        real(dp) :: weight_sum, weight_fraction
        integer  :: ldim_env(3), nimgs
        logical  :: l_resampled, l_env_valid
        l_ready     = .false.
        l_env_valid = .false.
        l_resampled = .false.
        mask_min    = 0.0
        mask_max    = 0.0
        weight_sum  = 0.0_dp
        env_fname   = string(NU_ENVMASK_FBODY)//int2str_pad(state_here,2)//string(MRC_EXT)
        if( .not. file_exists(env_fname) )then
            skip_reason = 'envelope_absent'
        else
            call find_ldim_nptcls(env_fname, ldim_env, nimgs)
            smpd_env = find_img_smpd(env_fname)
            if( ldim_env(1) /= ldim_env(2) .or. ldim_env(1) /= ldim_env(3) )then
                skip_reason = 'envelope_not_cubic'
            else
                extent_env = real(ldim_env(1))     * smpd_env
                extent_rec = real(params%box_crop) * params%smpd_crop
                if( abs(extent_env - extent_rec) > EXTENT_RELTOL * extent_rec )then
                    skip_reason = 'physical_extent_mismatch'
                    ! the numbers make a wrong-smpd envelope immediately
                    ! visible: a header smpd typo shifts extent_env off the
                    ! reconstruction FOV by exactly the smpd error
                    write(logfhandle,'(A,I0,A,F8.4,A,F10.3,A,I0,A,F8.4,A,F10.3,A)') &
                        &'>>> PCG SOLVENT ENVELOPE EXTENT MISMATCH: envelope box ', ldim_env(1), &
                        &' x smpd ', smpd_env, ' = ', extent_env, ' A vs reconstruction box ', &
                        &params%box_crop, ' x smpd ', params%smpd_crop, ' = ', extent_rec, ' A'
                else
                    ! constant-FOV lattice change: factor-free Fourier pad/clip
                    ! under the data-quotient convention (S3 item 8)
                    l_resampled = ldim_env(1) /= params%box_crop
                    call env%read_and_crop(env_fname, smpd_env, params%box_crop, params%smpd_crop)
                    m = env%get_rmat()
                    call env%kill
                    if( .not. all(ieee_is_finite(m)) )then
                        skip_reason = 'envelope_not_finite'
                    else
                        mask_min = minval(m)
                        mask_max = maxval(m)
                        if( mask_min < -VALUE_SLACK .or. mask_max > 1.0 + VALUE_SLACK )then
                            skip_reason = 'values_outside_unit_range'
                        else
                            ! S3: re-clip to [0,1]; resampling rings slightly
                            m = min(1.0, max(0.0, m))
                            ! solvent confidence is the molecular-envelope complement
                            weight_sum = sum(1.0_dp - real(m,dp))
                            if( weight_sum <= 0.0_dp )then
                                skip_reason = 'no_solvent_evidence'
                            else
                                l_env_valid = .true.
                                if( params%pcg_solvent_lambda_rel > 0.0 )then
                                    skip_reason = 'none'
                                else
                                    skip_reason = 'strength_zero'
                                endif
                            endif
                        endif
                    endif
                endif
            endif
        endif
        l_ready = l_env_valid .and. params%pcg_solvent_lambda_rel > 0.0
        if( l_ready )then
            call move_alloc(m, m_env)
        else if( params%pcg_solvent_lambda_rel > 0.0 .and. skip_reason /= 'envelope_absent' )then
            ! the prior is on by default, so an ABSENT envelope is the normal
            ! lag-one state (early iterations, gridding-era projects) and only
            ! rates the diagnostic block below; a PRESENT-but-unusable envelope
            ! is a broken input and warns loudly. Never substitute a mask type.
            THROW_WARN('NU evidence envelope is unusable ('//skip_reason//'); solvent prior skipped for this iteration')
        endif
        weight_fraction = weight_sum / real(params%box_crop,dp)**3
        write(logfhandle,'(A,I0,A)') '>>> PCG SOLVENT ENVELOPE ('//trim(context)//'/state ', &
            &state_here, ')'
        if( l_ready )then
            write(logfhandle,'(A,ES11.4)') '    pcg_solvent_prior_enabled=T pcg_solvent_lambda_rel=', &
                &params%pcg_solvent_lambda_rel
        else
            write(logfhandle,'(A,ES11.4)') '    pcg_solvent_prior_enabled=F pcg_solvent_lambda_rel=', &
                &params%pcg_solvent_lambda_rel
        endif
        write(logfhandle,'(A)')      '    pcg_solvent_mask_file='//env_fname%to_char()
        if( l_resampled )then
            write(logfhandle,'(A,I0,A,I0,A)') '    pcg_solvent_mask_resampled=T (box ', ldim_env(1), &
                &' -> ', params%box_crop, ', constant FOV)'
        else
            write(logfhandle,'(A)') '    pcg_solvent_mask_resampled=F'
        endif
        write(logfhandle,'(A,ES11.4,A,ES11.4)') '    pcg_solvent_mask_min=', mask_min, &
            &' pcg_solvent_mask_max=', mask_max
        write(logfhandle,'(A,ES11.4,A,ES11.4)') '    pcg_solvent_weight_sum=', real(weight_sum), &
            &' pcg_solvent_weight_fraction=', real(weight_fraction)
        ! lambda_eff and the mean/rms/penalty_final values are per-half solve
        ! outputs; written by the PCG SOLVENT PRIOR line after each ML replay
        write(logfhandle,'(A)')      '    pcg_solvent_skip_reason='//skip_reason
        call env_fname%kill
    end subroutine resolve_solvent_envelope

    !> The solve-output half of the pcg_solvent_* diagnostic block (S6.2):
    !! effective absolute strength and the final map's weighted solvent mean,
    !! RMS deviation and penalty, written after each priored ML replay.
    subroutine report_solvent_solve_stats( pcgop, x, rms_ref, context, half, supp_pct, rms_out )
        type(reconstructor_pcg), intent(in)  :: pcgop
        real,                    intent(in)  :: x(:,:,:)
        real,                    intent(in)  :: rms_ref
        character(len=*),        intent(in)  :: context, half
        real,                    intent(out) :: supp_pct, rms_out
        real :: mean_s, rms_s, penalty
        call pcgop%get_solvent_stats(x, mean_s, rms_s, penalty)
        ! suppression is measured against the shell-shrunk base reference
        ! (solvent_suppression_reference), which carries the Wiener component
        ! of the replay but no solvent prior: ~0% means the prior is inert
        supp_pct = 0.0
        if( rms_ref > 1.0e-12 ) supp_pct = 100.0 * (1.0 - rms_s / rms_ref)
        rms_out = rms_s
        write(logfhandle,'(A,ES11.4,A,ES11.4,A,ES11.4,A,ES11.4)') &
            &'>>> PCG SOLVENT PRIOR ('//trim(context)//'/'//trim(half)//'): pcg_solvent_lambda_eff=', &
            &pcgop%get_effective_solvent_lambda(), ' pcg_solvent_mean_final=', mean_s, &
            &' pcg_solvent_rms_final=', rms_s, ' pcg_solvent_penalty_final=', penalty
        write(logfhandle,'(A,ES11.4,A,F8.2)') '    pcg_solvent_rms_ref=', rms_ref, &
            &' pcg_solvent_suppression_pct=', supp_pct
    end subroutine report_solvent_solve_stats

    subroutine execute_rec3D_pcg_shared( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(image) :: half_even, half_odd, ml_even, ml_odd, merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        integer, allocatable :: selected_pinds(:), half_pinds(:)
        real, allocatable :: fsc(:), res0143s(:), solvent_env(:,:,:), solvent_supps(:), solvent_rmss(:)
        integer, allocatable :: solvent_supp_cnts(:)
        logical, allocatable :: state_written(:)
        integer :: nselected, state, n_state, n_even, n_odd, iptcl, istate
        real :: res05, cfar
        logical :: l_solvent_ready
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output
        logical :: l_sigma_loaded

        call validate_supported_mode()
        nselected = 0
        call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
        if( nselected < 1 ) THROW_HARD('no active particles selected for PCG reconstruct3D')

        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif

        allocate(res0143s(params%nstates), source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        allocate(solvent_supps(params%nstates), solvent_rmss(params%nstates), source=0.0)
        allocate(solvent_supp_cnts(params%nstates), source=0)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)

        do state = 1, params%nstates
            l_solvent_ready = .false.
            n_state = count_state(state)
            n_even = count_state_half(state, 0)
            n_odd  = count_state_half(state, 1)
            if( n_state == 0 )then
                write(logfhandle,'(A,I0,A)') '>>> PCG RECONSTRUCT3D: STATE ', state, ' HAS NO SELECTED PARTICLES; SKIPPING'
                cycle
            endif
            if( n_even + n_odd /= n_state ) THROW_HARD('PCG reconstruct3D found invalid halfset labels')
            if( n_even < 1 .or. n_odd < 1 ) THROW_HARD('PCG reconstruct3D requires particles in both halfsets')

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
            time_map_output = 0.0_dp
            if( params%l_ml_reg )then
                fname_even_unfil = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
                fname_odd_unfil  = refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.)
                t_state_phase = tic()
                call half_even%write(fname_even_unfil, del_if_exists=.true.)
                call half_odd%write(fname_odd_unfil, del_if_exists=.true.)
                time_map_output = real(toc(t_state_phase),dp)
            endif

            t_state_phase = tic()
            call calculate_state_fsc(state, half_even, half_odd, merged, fsc, res05, res0143s(state), cfar)
            call arr2file(fsc, fname_fsc)
            call write_fsc_summary(state, merged, fsc, res05, res0143s(state), cfar)
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                call resolve_solvent_envelope(params, state, 'shared', solvent_env, l_solvent_ready)
                call regularize_state_half(state, 0, 'even', fsc, half_even, ml_even)
                call regularize_state_half(state, 1, 'odd',  fsc, half_odd,  ml_odd)
                if( allocated(solvent_env) ) deallocate(solvent_env)
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
            time_map_output = time_map_output + real(toc(t_state_phase),dp)
            call write_output_diagnostics(state, 'shared', time_map_output, time_fsc_output)

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

        call killimgbatch(build)
        if( params%l_ml_reg ) call write_solvent_convergence_stats(params, solvent_supps, solvent_rmss, solvent_supp_cnts)
        if( .not. any(state_written) ) THROW_HARD('PCG reconstruct3D produced no populated states')
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res', res0143s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                endif
            enddo
        endif
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call register_project_outputs()

        deallocate(selected_pinds, res0143s, state_written, solvent_supps, solvent_rmss, solvent_supp_cnts)

    contains

        subroutine validate_supported_mode()
            if( params%nparts > 1 .or. cline%defined('part') )then
                THROW_HARD('shared PCG entry received distributed parameters')
            endif
            if( trim(params%pcgop) /= 'kernel' ) THROW_HARD('production rec_backend=pcg requires pcgop=kernel')
            ! refinement stays within the small-iteration budget the warm start
            ! makes sufficient; offline harness/benchmark runs legitimately ask
            ! for rtol-terminated converged solves (pcg_priors.md R3), so above
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
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:), x(:,:,:), rel_res_hist(:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch, niters
            real    :: shift(2), crop_factor, sdev_noise, edge_mean
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
            call pcgop%set_mask(params%msk_crop)
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
            sdev_noise = 0.0
            edge_mean  = 0.0
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
                    call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                    call build%imgbatch(ii)%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(build%imgbatch(ii))
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
            t_phase = tic()
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            time_finalize = real(toc(t_phase),dp)
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call validate_solved_map(x, 'shared', state_here, half, 'base')
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            call report_beyond_band_excess(volume, params, state_here, half, 'base')
            time_total = real(toc(t_half),dp)
            call write_half_diagnostics(state_here, half, 'base', size(pinds), result, rel_res_hist, &
                &time_metadata, time_particles, time_accum_init, time_accum, time_finalize, time_solve, time_total, &
                &pcgop%get_data_scale(), pcgop%get_effective_lambda(), pcgop=pcgop)
            call report_solve_summary('SHARED', state_here, half, 'base', size(pinds), niters, &
                &result%final_rel_residual, time_solve, result%stop_reason)
            if( present(outcome) ) outcome = result

            call pcgop%kill
            call selection%kill
            call orientation%kill
            deallocate(y_batch, sig2, x, rel_res_hist)
        end subroutine solve_state_half

        !> Reopen the exact raw statistics used for the base half-map, add the
        !! FSC/SSNR prior only on the master/shared owner, and warm-start from
        !! the corresponding unregularized solution. No particle data are read
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
            real :: rms_supp_ref, supp_pct, rms_half
            logical :: l_warm

            t_phase = tic()
            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call pcgop%set_mask(params%msk_crop)
            call pcgop%begin_reduction
            fname = refine3D_pcg_raw_accum_fname(state_here, 1, params%numlen, half)
            call pcgop%add_raw_accum(fname, state_here, eo_here, 1, 1, &
                &pcg_raw_provenance(params), nptcls)
            time_reduce = real(toc(t_phase),dp)
            if( nptcls < 1 ) THROW_HARD('PCG ML replay requires a populated raw half accumulator')
            t_phase = tic()
            call pcgop%set_ml_prior(fsc_here, params%tau, params%hp)
            ! before end_accum: the effective solvent strength is derived from
            ! the data scale alongside the relative ridge lambda
            if( l_solvent_ready ) call pcgop%set_solvent_prior(solvent_env, params%pcg_solvent_lambda_rel)
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            call pcgop%assert_prior_attachment_mode
            time_finalize = real(toc(t_phase),dp)
            call pcgop%get_ml_prior_stats(prior_npositive, prior_positive_min, prior_positive_max, &
                &prior_to_khat_l1, prior_to_khat_rms)
            x = base_volume%get_rmat()
            call override_ml_warm_start_from_previous(params, state_here, half, x, 'shared', l_warm)
            if( .not. l_warm ) call regularized_ml_initial_guess(params, fsc_here, solvent_env, &
                &l_solvent_ready, x, 'shared', half)
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call validate_solved_map(x, 'shared', state_here, half, 'ml')
            if( l_solvent_ready )then
                rms_supp_ref = solvent_suppression_reference(pcgop, params, fsc_here, base_volume)
                call report_solvent_solve_stats(pcgop, x, rms_supp_ref, 'shared', half, supp_pct, rms_half)
                solvent_supps(state_here)     = solvent_supps(state_here) + supp_pct
                solvent_rmss(state_here)      = solvent_rmss(state_here)  + rms_half
                solvent_supp_cnts(state_here) = solvent_supp_cnts(state_here) + 1
            endif
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
                &result%final_rel_residual, time_solve, result%stop_reason)
            call pcgop%kill
            call fname%kill
            deallocate(x, rel_res_hist)
        end subroutine regularize_state_half

        subroutine calculate_state_fsc( state_here, even, odd, avg, fsc_out, res05, res0143, cfar )
            integer,     intent(in)    :: state_here
            type(image), intent(in)    :: even, odd, avg
            real, allocatable, intent(out) :: fsc_out(:)
            real,                    intent(out) :: res05, res0143, cfar
            type(image)     :: work_even, work_odd
            type(image_msk) :: envmask
            type(fsc_area_score_result) :: cones
            real, allocatable :: fsc_t(:), fsc_n(:), res(:)
            integer :: nyq
            nyq = even%get_filtsz()
            if( params%l_envfsc )then
                call envmask%automask3D(params, avg, .false., lp_override=params%envmsklp)
                call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
                call phase_rand_fsc(even, odd, envmask, state_here, nyq, fsc_out, fsc_t, fsc_n)
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%zero_env_background(envmask)
                call work_odd%zero_env_background(envmask)
                call work_even%mul(envmask)
                call work_odd%mul(envmask)
                deallocate(fsc_t, fsc_n)
                call envmask%kill_bimg
            else
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%mask3D_soft(params%msk_crop, backgr=0.0)
                call work_odd%mask3D_soft(params%msk_crop, backgr=0.0)
                allocate(fsc_out(nyq), source=0.0)
            endif
            call cones%new(work_even, 256, 20.0, 0.143, 1)
            call cones%calc_fsc_area_score(work_even, work_odd, state=state_here)
            cfar = cones%cfar
            if( .not. params%l_envfsc ) call work_even%fsc(work_odd, fsc_out)
            res = avg%get_res()
            call get_resolution(fsc_out, res, res05, res0143)
            res05   = max(res05,   2.0 * params%smpd_crop)
            res0143 = max(res0143, 2.0 * params%smpd_crop)
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, &
                &' FSC=0.500 RESOLUTION = ', res05
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG RECONSTRUCT3D: STATE ', state_here, &
                &' FSC=0.143 RESOLUTION = ', res0143
            call cones%kill
            call work_even%kill
            call work_odd%kill
            deallocate(res)
        end subroutine calculate_state_fsc

        subroutine write_fsc_summary( state_here, avg, fsc_in, res05, res0143, cfar )
            integer,     intent(in) :: state_here
            type(image), intent(in) :: avg
            real,        intent(in) :: fsc_in(:), res05, res0143, cfar
            type(string) :: fname
            real, allocatable :: res(:)
            integer :: funit, k
            fname = refine3D_resolution_txt_fbody(state_here)
            res   = avg%get_res()
            call fopen(funit, file=fname, status='replace', action='write')
            do k = 1, min(size(res),size(fsc_in))
                write(funit,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                    &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', fsc_in(k)
            enddo
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(funit,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', cfar
            call fclose(funit)
            call fname%kill
            deallocate(res)
        end subroutine write_fsc_summary

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
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
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
            complex, allocatable :: y_batch(:,:,:)
            real, allocatable :: sig2(:,:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch
            real :: shift(2), crop_factor, sdev_noise, edge_mean
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
            sdev_noise = 0.0
            edge_mean = 0.0
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
                    call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                    call build%imgbatch(ii)%fft
                    y_batch(:,:,ii) = op%extract_native_plane(build%imgbatch(ii))
                enddo
                call op%accumulate_batch(y_batch, batchsz, batchlims(1))
            enddo
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
            allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            call op%end_accum(.true.)
            call op%set_op_mode(PCG_OP_KERNEL)
            call op%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, rel_res_hist=history, niters=niters)
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
    subroutine execute_rec3D_pcg_worker( params, build, cline, selected_pinds, cached )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: selected_pinds(:)
        logical, optional, intent(in)   :: cached
        integer, allocatable :: half_pinds(:)
        character(len=256) :: provenance
        integer :: state, eo, n_half
        logical :: l_sigma_loaded, l_cached

        call validate_pcg_common(params, check_solver=.false.)
        l_cached = .false.
        if( present(cached) ) l_cached = cached
        provenance = pcg_raw_provenance(params)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            l_sigma_loaded = allocated(build%esig%sigma2_noise)
            if( .not. l_sigma_loaded )then
                call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
            endif
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif
        if( size(selected_pinds) > 0 )then
            if( l_cached )then
                call prepimgbatch(params, build, MAXIMGBATCHSZ, box=params%box_crop, smpd=params%smpd_crop)
            else
                call prepimgbatch(params, build, MAXIMGBATCHSZ)
            endif
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
            complex, allocatable :: y_batch(:,:,:)
            real,    allocatable :: sig2(:,:)
            integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz
            integer :: i, ii, iptcl, ibatch
            real    :: shift(2), crop_factor, sdev_noise, edge_mean

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
            sdev_noise = 0.0
            edge_mean  = 0.0
            do ibatch = 1, size(pinds), MAXIMGBATCHSZ
                batchlims = [ibatch, min(size(pinds),ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                if( l_cached )then
                    call ptcl_cache_read_batch(params, build, size(pinds), pinds, batchlims)
                else if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', size(pinds), pinds, &
                        &batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, size(pinds), pinds, batchlims)
                endif
                do ii = 1, batchsz
                    if( .not. l_cached ) call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                    call build%imgbatch(ii)%fft
                    y_batch(:,:,ii) = pcgop%extract_native_plane(build%imgbatch(ii))
                enddo
                call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
            enddo
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
    !! order, then perform all folding, finalization and PCG locally. Independent
    !! state/half reductions are completed and released one at a time.
    subroutine execute_rec3D_pcg_distributed_master( params, build, cline, trail_bootstrap_states )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        logical, optional, intent(out)  :: trail_bootstrap_states(:)
        type(image) :: half_even, half_odd, ml_even, ml_odd, merged
        type(image) :: previous_even, previous_odd, previous_merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        real, allocatable :: fsc(:), res0143s(:), solvent_env(:,:,:), solvent_supps(:), solvent_rmss(:)
        real, allocatable :: realized_fractions(:), update_weights(:)
        integer, allocatable :: solvent_supp_cnts(:)
        logical, allocatable :: state_written(:)
        character(len=256) :: provenance, chain_provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        integer :: n_active_state, n_sampled_state
        real :: res05, cfar
        logical :: l_has_updates, l_bootstrap, l_even_chain, l_odd_chain, l_solvent_ready
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output

        call validate_pcg_common(params)
        if( present(trail_bootstrap_states) )then
            if( size(trail_bootstrap_states) /= params%nstates ) &
                &THROW_HARD('PCG trailing-bootstrap state output has invalid size')
            trail_bootstrap_states = .false.
        endif
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
        allocate(state_written(params%nstates), source=.false.)
        allocate(solvent_supps(params%nstates), solvent_rmss(params%nstates), source=0.0)
        allocate(solvent_supp_cnts(params%nstates), source=0)
        do state = 1, params%nstates
            l_solvent_ready = .false.
            l_bootstrap = .false.
            if( params%l_trail_rec )then
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'even')
                l_even_chain = file_exists(raw_fname)
                raw_fname = refine3D_pcg_trail_accum_fname(state, 'odd')
                l_odd_chain = file_exists(raw_fname)
                if( l_even_chain .neqv. l_odd_chain ) THROW_HARD('PCG trailing chain pair is incomplete')
                l_bootstrap = .not. l_even_chain
            endif
            if( present(trail_bootstrap_states) ) trail_bootstrap_states(state) = l_bootstrap
            call reduce_solve_state_half(state, 0, 'even', half_even, n_even, 'base')
            call reduce_solve_state_half(state, 1, 'odd',  half_odd,  n_odd,  'base')
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
            time_map_output = 0.0_dp
            if( params%l_ml_reg )then
                fname_even_unfil = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
                fname_odd_unfil  = refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.)
                t_state_phase = tic()
                call half_even%write(fname_even_unfil, del_if_exists=.true.)
                call half_odd%write(fname_odd_unfil, del_if_exists=.true.)
                time_map_output = real(toc(t_state_phase),dp)
            endif
            t_state_phase = tic()
            if( l_bootstrap )then
                call load_previous_state_halves(state, previous_even, previous_odd, previous_merged)
                call calculate_distributed_fsc(state, previous_even, previous_odd, previous_merged, &
                    &fsc, res05, res0143s(state), cfar)
            else
                call calculate_distributed_fsc(state, half_even, half_odd, merged, &
                    &fsc, res05, res0143s(state), cfar)
            endif
            call arr2file(fsc, fname_fsc)
            if( l_bootstrap )then
                call write_distributed_fsc_summary(state, previous_merged, fsc, res05, res0143s(state), cfar)
            else
                call write_distributed_fsc_summary(state, merged, fsc, res05, res0143s(state), cfar)
            endif
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                call resolve_solvent_envelope(params, state, 'distributed', solvent_env, l_solvent_ready)
                call reduce_solve_state_half(state, 0, 'even', ml_even, n_even, 'ml', fsc, half_even)
                call reduce_solve_state_half(state, 1, 'odd',  ml_odd,  n_odd,  'ml', fsc, half_odd)
                if( allocated(solvent_env) ) deallocate(solvent_env)
                call merged%kill
                call merged%copy(ml_even)
                call merged%add(ml_odd)
                call merged%mul(0.5)
            endif
            if( l_bootstrap .and. update_weights(state) < 0.99 )then
                if( params%l_ml_reg )then
                    call blend_bootstrap_half(ml_even, previous_even, update_weights(state))
                    call blend_bootstrap_half(ml_odd,  previous_odd,  update_weights(state))
                else
                    call blend_bootstrap_half(half_even, previous_even, update_weights(state))
                    call blend_bootstrap_half(half_odd,  previous_odd,  update_weights(state))
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
            if( params%l_ml_reg )then
                call ml_even%write(fname_even, del_if_exists=.true.)
                call ml_odd%write(fname_odd, del_if_exists=.true.)
            else
                call half_even%write(fname_even, del_if_exists=.true.)
                call half_odd%write(fname_odd, del_if_exists=.true.)
            endif
            call merged%write(fname_vol, del_if_exists=.true.)
            time_map_output = time_map_output + real(toc(t_state_phase),dp)
            call write_output_diagnostics(state, 'distributed', time_map_output, time_fsc_output)
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
        if( params%l_ml_reg ) call write_solvent_convergence_stats(params, solvent_supps, solvent_rmss, solvent_supp_cnts)
        if( .not. any(state_written) ) THROW_HARD('distributed PCG produced no populated states')
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res', res0143s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    if( state_written(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
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
        deallocate(res0143s, state_written, realized_fractions, update_weights, solvent_supps, solvent_rmss, &
            &solvent_supp_cnts)

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

        subroutine blend_bootstrap_half( current, previous, weight_current )
            type(image), intent(inout) :: current, previous
            real,        intent(in)    :: weight_current
            call current%mul(weight_current)
            call previous%mul(1.0-weight_current)
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

        subroutine reduce_solve_state_half( state_here, eo_here, half, volume, nptcls, solve_kind, &
                &fsc_prior, warm_start )
            integer,          intent(in)    :: state_here, eo_here
            character(len=*), intent(in)    :: half, solve_kind
            type(image),      intent(inout) :: volume
            integer,          intent(out)   :: nptcls
            real, optional,   intent(in)    :: fsc_prior(:)
            type(image), optional, intent(in) :: warm_start
            type(reconstructor_pcg) :: pcgop
            type(pcg_solver_outcome) :: result
            type(string) :: fname
            real, allocatable :: x(:,:,:), rel_res_hist(:)
            integer :: part_here, n_part, niters, prior_npositive, n_full_half
            integer(timer_int_kind) :: t_phase
            real(dp) :: time_reduce, time_finalize, time_solve
            real :: prior_positive_min, prior_positive_max, prior_to_khat_l1, prior_to_khat_rms
            real :: realized_fraction, update_weight, current_scale, rms_supp_ref, supp_pct, rms_half
            logical :: l_ml_solve, l_chain_exists, l_seed_chain, l_warm

            l_ml_solve = present(fsc_prior)
            if( l_ml_solve .neqv. present(warm_start) )then
                THROW_HARD('distributed PCG ML replay requires both FSC and warm start')
            endif

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call pcgop%set_mask(params%msk_crop)
            call pcgop%begin_reduction
            nptcls = 0
            t_phase = tic()
            if( l_ml_solve .and. params%l_trail_rec )then
                fname = refine3D_pcg_trail_accum_fname(state_here, half)
                call pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                    &chain_provenance, 1.0, nptcls)
                call fname%kill
                if( l_bootstrap .or. 1.0-update_weights(state_here) <= 0.01 )then
                    call pcgop%scale_raw_accum(realized_fractions(state_here))
                endif
            else
                ! Without trailing, replay the current worker statistics
                ! directly for ML regularization. The persistent trail file is
                ! reserved for full-mass continuation data and must not double
                ! as an ML scratch artifact carrying fractional mass.
                do part_here = 1, params%nparts
                    fname = refine3D_pcg_raw_accum_fname(state_here, part_here, params%numlen, half)
                    call pcgop%add_raw_accum(fname, state_here, eo_here, part_here, params%nparts, &
                        &provenance, n_part)
                    nptcls = nptcls + n_part
                    call fname%kill
                enddo
            endif
            time_reduce = real(toc(t_phase),dp)
            if( nptcls == 0 )then
                call pcgop%kill
                return
            endif
            if( .not. l_ml_solve )then
                n_full_half = count_full_state_half(state_here, eo_here)
                if( n_full_half < nptcls ) THROW_HARD('PCG current half population exceeds its full population')
                realized_fraction = realized_fractions(state_here)
                update_weight = update_weights(state_here)
                l_seed_chain = pcg_trail_seed_requested(cline)
                fname = refine3D_pcg_trail_accum_fname(state_here, half)
                l_chain_exists = file_exists(fname)
                if( params%l_trail_rec )then
                    if( realized_fraction <= 0.0 ) THROW_HARD('PCG trailing update has zero realized state fraction')
                    if( l_chain_exists .and. 1.0-update_weight > 0.01 )then
                        current_scale = update_weight / realized_fraction
                        call pcgop%scale_raw_accum(current_scale)
                        call pcgop%add_raw_accum_weighted(fname, state_here, eo_here, 1, 1, &
                            &chain_provenance, 1.0-update_weight, n_part)
                    else
                        current_scale = 1.0 / realized_fraction
                        call pcgop%scale_raw_accum(current_scale)
                    endif
                    call pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, n_full_half, chain_provenance)
                    if( .not. l_chain_exists .or. 1.0-update_weight <= 0.01 )then
                        call pcgop%scale_raw_accum(realized_fraction)
                    endif
                    write(logfhandle,'(A,I0,A,A,A,F8.4,A,F8.4)') '>>> PCG TRAIL | STATE=', state_here, &
                        &' | HALF=', trim(half), ' | F=', realized_fraction, ' | U=', update_weight
                else if( l_seed_chain )then
                    call pcgop%write_raw_accum(fname, state_here, eo_here, 1, 1, n_full_half, chain_provenance)
                endif
                call fname%kill
            endif
            t_phase = tic()
            if( l_ml_solve )then
                call pcgop%set_ml_prior(fsc_prior, params%tau, params%hp)
                ! before end_accum: effective strength derives from the data scale
                if( l_solvent_ready ) call pcgop%set_solvent_prior(solvent_env, params%pcg_solvent_lambda_rel)
            endif
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            if( l_ml_solve ) call pcgop%assert_prior_attachment_mode
            time_finalize = real(toc(t_phase),dp)
            if( l_ml_solve )then
                call pcgop%get_ml_prior_stats(prior_npositive, prior_positive_min, prior_positive_max, &
                    &prior_to_khat_l1, prior_to_khat_rms)
            endif
            if( l_ml_solve )then
                x = warm_start%get_rmat()
                call override_ml_warm_start_from_previous(params, state_here, half, x, 'distributed', l_warm)
                if( .not. l_warm ) call regularized_ml_initial_guess(params, fsc_prior, solvent_env, &
                    &l_solvent_ready, x, 'distributed', half)
            else
                allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            endif
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call validate_solved_map(x, 'distributed', state_here, half, solve_kind)
            if( l_ml_solve .and. l_solvent_ready )then
                ! warm_start is the base half solution at both ML call sites
                rms_supp_ref = solvent_suppression_reference(pcgop, params, fsc_prior, warm_start)
                call report_solvent_solve_stats(pcgop, x, rms_supp_ref, 'distributed', half, supp_pct, rms_half)
                solvent_supps(state_here)     = solvent_supps(state_here) + supp_pct
                solvent_rmss(state_here)      = solvent_rmss(state_here)  + rms_half
                solvent_supp_cnts(state_here) = solvent_supp_cnts(state_here) + 1
            endif
            call volume%new([params%box_crop,params%box_crop,params%box_crop], params%smpd_crop)
            call volume%set_rmat(x, .false.)
            call report_beyond_band_excess(volume, params, state_here, half, solve_kind)
            if( l_ml_solve )then
                call write_distributed_diagnostics(state_here, half, solve_kind, nptcls, result, rel_res_hist, &
                    &time_reduce, time_finalize, time_solve, pcgop%get_data_scale(), &
                    &pcgop%get_effective_lambda(), prior_npositive, prior_positive_min, prior_positive_max, &
                    &prior_to_khat_l1, prior_to_khat_rms, pcgop)
            else
                call write_distributed_diagnostics(state_here, half, solve_kind, nptcls, result, rel_res_hist, &
                    &time_reduce, time_finalize, time_solve, pcgop%get_data_scale(), &
                    &pcgop%get_effective_lambda(), pcgop=pcgop)
            endif
            call report_solve_summary('DISTRIBUTED', state_here, half, solve_kind, nptcls, niters, &
                &result%final_rel_residual, time_solve, result%stop_reason)
            call pcgop%kill
            deallocate(x, rel_res_hist)
        end subroutine reduce_solve_state_half

        subroutine calculate_distributed_fsc( state_here, even, odd, avg, fsc_out, res05, res0143, cfar )
            integer,     intent(in)    :: state_here
            type(image), intent(in)    :: even, odd, avg
            real, allocatable, intent(out) :: fsc_out(:)
            real,                    intent(out) :: res05, res0143, cfar
            type(image)     :: work_even, work_odd
            type(image_msk) :: envmask
            type(fsc_area_score_result) :: cones
            real, allocatable :: fsc_t(:), fsc_n(:), res(:)
            integer :: nyq
            nyq = even%get_filtsz()
            if( params%l_envfsc )then
                call envmask%automask3D(params, avg, .false., lp_override=params%envmsklp)
                call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
                call phase_rand_fsc(even, odd, envmask, state_here, nyq, fsc_out, fsc_t, fsc_n)
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%zero_env_background(envmask)
                call work_odd%zero_env_background(envmask)
                call work_even%mul(envmask)
                call work_odd%mul(envmask)
                deallocate(fsc_t, fsc_n)
                call envmask%kill_bimg
            else
                call work_even%copy(even)
                call work_odd%copy(odd)
                call work_even%mask3D_soft(params%msk_crop, backgr=0.0)
                call work_odd%mask3D_soft(params%msk_crop, backgr=0.0)
                allocate(fsc_out(nyq), source=0.0)
            endif
            call cones%new(work_even, 256, 20.0, 0.143, 1)
            call cones%calc_fsc_area_score(work_even, work_odd, state=state_here)
            cfar = cones%cfar
            if( .not. params%l_envfsc ) call work_even%fsc(work_odd, fsc_out)
            res = avg%get_res()
            call get_resolution(fsc_out, res, res05, res0143)
            res05   = max(res05,   2.0 * params%smpd_crop)
            res0143 = max(res0143, 2.0 * params%smpd_crop)
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG DISTRIBUTED: STATE ', state_here, &
                &' FSC=0.500 RESOLUTION = ', res05
            write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG DISTRIBUTED: STATE ', state_here, &
                &' FSC=0.143 RESOLUTION = ', res0143
            call cones%kill
            call work_even%kill
            call work_odd%kill
            deallocate(res)
        end subroutine calculate_distributed_fsc

        subroutine write_distributed_fsc_summary( state_here, avg, fsc_in, res05, res0143, cfar )
            integer,     intent(in) :: state_here
            type(image), intent(in) :: avg
            real,        intent(in) :: fsc_in(:), res05, res0143, cfar
            type(string) :: fname
            real, allocatable :: res(:)
            integer :: funit, k
            fname = refine3D_resolution_txt_fbody(state_here)
            res   = avg%get_res()
            call fopen(funit, file=fname, status='replace', action='write')
            do k = 1, min(size(res),size(fsc_in))
                write(funit,'(A,1X,F6.2,1X,A,1X,F7.3)') &
                    &'>>> RESOLUTION:', res(k), '>>> CORRELATION:', fsc_in(k)
            enddo
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res05
            write(funit,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res0143
            write(funit,'(A,1X,F6.2)') '>>> CONICAL FSC AREA RATIO (cFAR) SCORE  :', cfar
            call fclose(funit)
            call fname%kill
            deallocate(res)
        end subroutine write_distributed_fsc_summary

        subroutine write_distributed_diagnostics( state_here, half, solve_kind, nptcls, result, history, &
                &reduce_time, finalize_time, solve_time, data_scale, lambda_eff, prior_npositive, &
                &prior_positive_min, prior_positive_max, prior_to_khat_l1, prior_to_khat_rms, pcgop )
            integer,                  intent(in) :: state_here, nptcls
            character(len=*),         intent(in) :: half, solve_kind
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
            write(funit,'(A,I0)')     'nparts=',                params%nparts
            write(funit,'(A,I0)')     'nptcls=',                nptcls
            write(funit,'(A,I0)')     'requested_maxits=',      result%requested_maxits
            write(funit,'(A,ES14.6)') 'requested_rtol=',        params%rtol
            write(funit,'(A,I0)')     'iteration_count=',       result%iteration_count
            write(funit,'(A,A)')      'stop_reason=',           trim(result%stop_reason)
            write(funit,'(A,L1)')     'converged=',             result%converged
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
    !! doc/implementation_notes/pcg_priors.md §2, beyond-band excess). Such
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

    subroutine report_solve_summary( execution_mode, state, half, solve_kind, nptcls, niters, &
            &residual, solve_time, stop_reason )
        character(len=*), intent(in) :: execution_mode, half, solve_kind, stop_reason
        integer,          intent(in) :: state, nptcls, niters
        real,             intent(in) :: residual
        real(dp),         intent(in) :: solve_time
        character(len=4) :: half_label, kind_label
        half_label = adjustl(half)
        kind_label = adjustl(solve_kind)
        write(logfhandle,'(4A,I2,A,A4,A,A4,A,I6,A,I2,A,ES10.3,A,F7.2,2A)') &
            &'>>> PCG ', trim(execution_mode), ' | ', 'STATE=', state, ' | HALF=', half_label, &
            &' | KIND=', kind_label, ' | N=', nptcls, ' | ITS=', niters, ' | RESID=', residual, &
            &' | TIME=', solve_time, ' s | STOP=', trim(stop_reason)
        call flush(logfhandle)
    end subroutine report_solve_summary

    subroutine write_output_diagnostics( state, execution_mode, map_time, fsc_time )
        integer,          intent(in) :: state
        character(len=*), intent(in) :: execution_mode
        real(dp),         intent(in) :: map_time, fsc_time
        type(string) :: fname
        integer :: funit
        fname = 'reconstruct3D_pcg_state'//int2str_pad(state,2)//'_output_diagnostics.txt'
        call fopen(funit, file=fname, status='replace', action='write')
        write(funit,'(A,A)')     'execution_mode=', trim(execution_mode)
        write(funit,'(A,F12.6)') 'halfmap_merged_output_seconds=', map_time
        write(funit,'(A,F12.6)') 'fsc_cfar_summary_seconds=', fsc_time
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

    ! Continuation identity deliberately excludes which_iter: the chain must
    ! survive iteration boundaries, while geometry, objective and particle
    ! source changes must invalidate it.
    function pcg_chain_provenance( params ) result(provenance)
        type(parameters), intent(in) :: params
        character(len=256) :: provenance
        provenance = 'pcgtrail-v1|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|ptcl_src='//trim(params%ptcl_src)//'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))//'|box_crop='//trim(int2str(params%box_crop))// &
            &'|smpd_crop='//trim(real2str(params%smpd_crop))// &
            &'|msk='//trim(real2str(params%msk))//'|ctf='//trim(params%ctf)
    end function pcg_chain_provenance

end module simple_rec3D_pcg_strategy
