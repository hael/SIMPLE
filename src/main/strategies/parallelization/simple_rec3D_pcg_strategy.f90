!@descr: shared-memory production strategy body for kernel PCG reconstruct3D
module simple_rec3D_pcg_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, &
    &pcg_raw_accum_compatible
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, &
    &discrete_read_imgbatch_source, killimgbatch
use simple_sigma2_files,      only: load_sigma2_groups
use simple_math_ft,           only: resample_sigma2
use simple_estimate_ssnr,     only: fsc2shrink_filter, get_resolution
use simple_image,             only: image
use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, evaluate_halfmap_pair, &
    &write_halfmap_diagnostics
use simple_image_msk,         only: image_msk
use simple_nu_filter,         only: setup_nu_dmats, optimize_nu_cutoff_finds, cleanup_nu_filter, &
    &set_nu_solvent_envelope, retain_nu_filter_setup, NU_DEV_OUTPUT, &
    &extend_nu_filter_highres_shells, get_nu_filter_bank_finest_lp, &
    &build_nu_evidence_state, nu_evidence_state, nu_evidence_summary, get_nu_evidence_summary, &
    &expand_nu_evidence_band_weights, &
    &print_nu_evidence_summary, print_nu_evidence_lowpass_histogram, assert_nu_evidence_replay_ready, &
    &nu_evidence_finest_supported_lp, &
    &write_nu_evidence_envmask, &
    &NU_EVIDENCE_SOURCE_BASE, NU_EVIDENCE_SOURCE_PREV
use simple_vol_pproc_policy,  only: vol_pproc_plan, plan_state_postprocess, &
    &NU_ENVMASK_ACTION_REGENERATE
use simple_refine3D_fnames,   only: refine3D_state_halfvol_fname, refine3D_state_vol_fname, &
    &refine3D_fsc_fname, refine3D_resolution_txt_fbody, refine3D_pcg_raw_accum_fname, &
    &refine3D_pcg_trail_accum_fname
implicit none

public :: execute_rec3D_pcg_shared, execute_rec3D_pcg_worker, execute_rec3D_pcg_distributed_master
public :: validate_rec3D_pcg_fractional_updates, get_pcg_nu_evidence_bench_seconds
private
#include "simple_local_flags.inc"

real,    parameter :: PCG_LAMBDA = 1.0e-3
! Suppression-targeted auto-lambda controller (pcg_priors.md dev item 2).
! Plant model validated on the PfCRT/1WCM scan: amplitude suppression
! follows the one-pole law s = g*lambda/(1+g*lambda), so one secant step
! lands near the setpoint. The controller is memoryless: it reads the
! PREVIOUS iteration's stats file (lambda used + measured suppression),
! identifies g, and solves for the setpoint. Active only when
! pcg_nu_lambda_rel was left to its dynamic default (l_pcg_nu_autolambda);
! an explicit strength pins lambda and keeps every recorded control
! reproducible. The setpoint is owned by the auto-target outer loop below
! unless pinned via pcg_nu_supp_target.
real,    parameter :: NU_AUTOLAMBDA_DEADBAND_PCT = 5.0  !< hold inside target +/- deadband
real,    parameter :: NU_AUTOLAMBDA_STEP_CLAMP   = 5.0  !< max multiplicative step per iteration
real,    parameter :: NU_AUTOLAMBDA_LAMBDA_MIN   = 0.01 !< hard strength bounds
real,    parameter :: NU_AUTOLAMBDA_LAMBDA_MAX   = 30.0
real,    parameter :: NU_AUTOLAMBDA_SUPP_FLOOR   = 0.1  !< % floor so a near-zero readout cannot divide the model
! Auto-target outer loop (pcg_priors.md dev item 2): no fixed suppression
! setpoint transfers across datasets, so the setpoint itself is controlled.
! AIMD on the shipped-pair FSC=0.143 trajectory (the persisted
! over-regularization diagnostic): the target ratchets up additively only
! while the shipped-pair crossing improves, backs off multiplicatively when
! it degrades, and holds when it stalls (lp-limited stage plateaus therefore
! hold). An explicit pcg_nu_supp_target pins the setpoint.
real,    parameter :: NU_AUTOTARGET_MIN      = 5.0  !< setpoint bounds; low bound = the banner's inert threshold
real,    parameter :: NU_AUTOTARGET_MAX      = 75.0
real,    parameter :: NU_AUTOTARGET_STEP_ADD = 5.0  !< additive setpoint increase per improving iteration (% points)
real,    parameter :: NU_AUTOTARGET_BACKOFF  = 0.6  !< multiplicative setpoint decrease on shipped-pair degradation
real,    parameter :: NU_AUTOTARGET_SHIP_TOL = 0.02 !< relative deadband on the shipped-pair crossing (stall band)
integer, parameter :: PCG_MASTER_NTHR_CAP    = 32   !< master-phase thread-boost ceiling; OpenMP scaling saturates beyond this at typical box sizes
! Solve-support envelope (pcg_priors.md dev item 5): the conservative density
! envelope (automask3D at envmsklp) replaces the spherical support in the PCG
! solves, so the mask constrains the ESTIMATOR rather than post-processing.
! automsk=yes enables it; the REGULARIZED pass always takes it, the UNFIL
! pass only under envfsc=yes (envfsc=no keeps that pass spherical, so the
! FSC pair stays unconstrained). That envelope is generous by construction --
! protein plus micelle, dilated, soft-edged, and about half the spherical
! support on the datasets measured so far -- so it removes the far solvent,
! where deapodization amplifies hardest, while leaving the NU evidence a
! populated support. The evidence readiness contract (null_fraction bounds)
! is the guard if a specimen's envelope ever proves too tight for the null
! calibration.
logical, parameter :: DEBUG = .false.

! Rebuild-on-advance cadence for the NU replay evidence (streptavidin
! benchmark 2026-08-30: the per-iteration evidence build is a leading term of
! the pcg assembly overhead). The frozen compact evidence state rides across
! refinement iterations and is rebuilt when the resolution state advances
! (the FSC=0.143 crossing reaches a finer shell than the ridden evidence was
! built at), when the alignment band is the binding constraint (the crossing
! sits within one shell of the matching low-pass the evidence handed off --
! riding then would seal a search-bandwidth stall, PfCRT regression record
! 2026-08-31), when the evidence geometry or source changes, when the
! envmask cadence requires regeneration from live margins, or when the age
! cap is hit. The cache is in-memory and master-side only: a restart or a
! fresh process rebuilds on first use. The cache may only ever trade speed
! for staleness in the Q_NU band weights; it must never govern the search
! bandwidth from a frozen statistic.
integer, parameter :: NU_EVIDENCE_REBUILD_MAX_LAG = 5 !< iterations a frozen evidence state may ride before a forced rebuild
type :: nu_evidence_cache_entry
    logical :: valid = .false.
    integer :: box = 0
    integer :: fsc_find = 0
    integer :: handoff_find = 0 !< Fourier index of the matching low-pass this evidence handed off
    integer :: iter_built = -1
    character(len=32) :: source = ''
    type(nu_evidence_state) :: evstate
end type nu_evidence_cache_entry
type(nu_evidence_cache_entry), allocatable :: nu_evidence_cache(:)
! per-invocation NU-evidence wall-clock accumulator, surfaced to the
! standardized refine3D benchmark report (plot_refine3d_bench.py)
real(dp) :: nu_evidence_bench_seconds = 0.0_dp

contains

    !> Total NU replay evidence-phase seconds of the most recent PCG
    !! execution (all states), for the strategy-level benchmark report.
    real function get_pcg_nu_evidence_bench_seconds() result( seconds )
        seconds = real(nu_evidence_bench_seconds)
    end function get_pcg_nu_evidence_bench_seconds

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
        if( NU_DEV_OUTPUT ) &
        write(logfhandle,'(A)') '>>> PCG ML WARM START ('//trim(context)//'/'//trim(half)//&
            &'): previous-iteration ML half map '//prev_fname%to_char()
        call prev_fname%kill
    end subroutine override_ml_warm_start_from_previous

    !> Regularized initial guess for a COLD ML replay (no previous ML half
    !! map to warm-start from: the standalone harness, the first refinement
    !! iteration, abinitio3D stage handoffs). The documented cold-start gap
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

    !> Resolution readout of the shipped (regularized) half pair, measured the
    !! same way as the harness diagnostic (soft spherical mask, standard FSC
    !! crossings). The shipped pair shares its regularizers between halves, so
    !! this is NEVER a resolution claim — its crossing pulling materially
    !! finer than the base pair's is the portable over-regularization signal
    !! (calibrated on the retired solvent prior's strength ladders; retained
    !! as the NU replay's over-regularization diagnostic).
    subroutine shipped_pair_res( params, even_in, odd_in, res05, res0143 )
        class(parameters), intent(in)  :: params
        type(image),       intent(in)  :: even_in, odd_in
        real,              intent(out) :: res05, res0143
        type(image) :: he, ho
        real, allocatable :: corrs(:), res_arr(:)
        integer :: nyq
        call he%copy(even_in)
        call ho%copy(odd_in)
        call he%mask3D_soft(params%msk_crop, backgr=0.0)
        call ho%mask3D_soft(params%msk_crop, backgr=0.0)
        call he%fft()
        call ho%fft()
        nyq = he%get_filtsz()
        allocate(corrs(nyq), source=0.0)
        call he%fsc(ho, corrs)
        res_arr = he%get_res()
        call get_resolution(corrs, res_arr, res05, res0143)
        res05   = max(res05,   2.0*params%smpd_crop)
        res0143 = max(res0143, 2.0*params%smpd_crop)
        call he%kill
        call ho%kill
        deallocate(corrs, res_arr)
    end subroutine shipped_pair_res

    !> PCG half-map diagnostics through the backend-neutral common evaluator:
    !! this owns the PCG mask policy (params%msk_crop spherical radius on the
    !! restored real-space base pair), the automask artifact write on the
    !! envfsc path, and the execution-context resolution log lines. The shared
    !! and distributed paths pass identical scientific policy and differ only
    !! in the context label.
    subroutine calculate_pcg_state_diagnostics( params, state_here, context, even, odd, avg, diagnostics )
        class(parameters),                intent(in)  :: params
        integer,                          intent(in)  :: state_here
        character(len=*),                 intent(in)  :: context
        class(image),                     intent(in)  :: even, odd, avg
        type(halfmap_diagnostics_result), intent(out) :: diagnostics
        type(image) :: envmask
        if( params%l_envfsc )then
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, params%msk_crop, &
                &diagnostics, envmask=envmask)
            call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state_here,2)//MRC_EXT))
            call envmask%kill
        else
            call evaluate_halfmap_pair(params, state_here, even, odd, avg, params%msk_crop, diagnostics)
        endif
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.500 RESOLUTION = ', diagnostics%res_fsc05
        write(logfhandle,'(A,I0,A,F8.3)') '>>> PCG '//trim(context)//': STATE ', state_here, &
            &' FSC=0.143 RESOLUTION = ', diagnostics%res_fsc0143
    end subroutine calculate_pcg_state_diagnostics

    !> Stage-6 direct NU-evidence replay (pcg_priors.md S5-S6): construct the
    !! frozen compact evidence state and expand it into the graded band
    !! lack-of-evidence weights the solver consumes. The evidence pair is
    !! always the FSC pair of the volassemble: the CURRENT unregularized base
    !! half pair (source=base_unfil; in trailing mode that pair is the
    !! full-mass blended base solution, i.e. exactly the statistics the ML
    !! replay re-reads), or -- in the trailing bootstrap only -- the previous
    !! iteration's shipped half pair (source=previous_shipped), the same
    !! lag-one pair the bootstrap FSC is computed from. Built once per state
    !! after both base solves and before either replay; the two half replays
    !! share the one immutable evidence identity. No envelope artifact is
    !! read, and no silent fallback exists: evidence-construction failure is
    !! a hard error. With automsk enabled the NU evidence envelope artifact
    !! is regenerated here (policy 2026-08-29): the raw per-voxel evidence
    !! margins are live at this point, so the matching-reference envelope
    !! derives from the same frozen evidence as the Q_NU precision, without
    !! a second NU analysis; cadence and artifact naming follow the same
    !! plan_state_postprocess contract as the post-hoc NU paths.
    subroutine build_nu_replay_evidence( params, state_here, context, vol_even, vol_odd, &
            &res0143, band_w, band_limits, finest_lp, evidence_source, evidence_seconds, &
            &l_retain_setup )
        class(parameters),          intent(in)  :: params
        integer,                    intent(in)  :: state_here
        character(len=*),           intent(in)  :: context
        type(image),                intent(in)  :: vol_even, vol_odd
        real,                       intent(in)  :: res0143 !< FSC=0.143 crossing of the evidence pair; steers the rebuild cadence
        real, allocatable,          intent(out) :: band_w(:,:,:,:)
        real, allocatable,          intent(out) :: band_limits(:)
        real, optional,             intent(out) :: finest_lp
        character(len=*), optional, intent(in)  :: evidence_source
        real(dp), optional,         intent(out) :: evidence_seconds
        logical, optional,          intent(in)  :: l_retain_setup !< leave the optimized setup armed for the matching-reference pass
        type(nu_evidence_state)   :: evstate
        type(nu_evidence_summary) :: evsumm
        type(vol_pproc_plan)      :: pp_plan
        type(image)               :: vol_env_avg
        type(image_msk)           :: solvent_env
        character(len=32) :: source_here
        integer :: nsteps_ext, ldim_here(3), fsc_find
        real    :: handoff_lp
        logical :: l_rebuild, l_envmask_regen, l_retain
        integer(timer_int_kind) :: t_evidence
        real(dp) :: seconds_here
        t_evidence  = tic()
        source_here = NU_EVIDENCE_SOURCE_BASE
        if( present(evidence_source) ) source_here = trim(evidence_source)
        ldim_here = vol_even%get_ldim()
        fsc_find  = calc_fourier_index(res0143, ldim_here(1), vol_even%get_smpd())
        ! envmask cadence decides first: a required regeneration forces a
        ! rebuild, because the envelope derives from the live per-voxel
        ! margins that only exist during evidence construction
        call plan_state_postprocess(params, state_here, params%which_iter, pp_plan)
        l_envmask_regen = pp_plan%nu_envmask_action == NU_ENVMASK_ACTION_REGENERATE
        l_rebuild       = nu_evidence_needs_rebuild(params, state_here, ldim_here(1), source_here, &
            &fsc_find, l_envmask_regen)
        if( l_rebuild )then
            write(logfhandle,'(A,I0,A)') '>>> PCG NU REPLAY ('//trim(context)//'): BUILDING EVIDENCE FROM THE '//&
                &'FSC HALF PAIR OF STATE ', state_here, ' (source='//trim(source_here)//')'
            call setup_nu_dmats(vol_even, vol_odd, params%mskdiam, [real ::], evidence_source=trim(source_here))
            ! solvent constraint: conservative density envelope from the
            ! evidence pair average, computed on the fly and frozen into the
            ! evidence state (pcg_priors.md dev item 4)
            call vol_env_avg%copy(vol_even)
            call vol_env_avg%add(vol_odd)
            call vol_env_avg%mul(0.5)
            call solvent_env%automask3D(params, vol_env_avg, .false., lp_override=params%envmsklp, l_report=.false.)
            call set_nu_solvent_envelope(solvent_env)
            call solvent_env%kill_bimg
            call vol_env_avg%kill
            call optimize_nu_cutoff_finds()
            ! Stage 6.6 (pcg_priors.md): with nu_refine=yes the evidence candidate
            ! bank is extended by the proven high-resolution shell walk -- one
            ! Fourier shell at a time from the populated frontier, accepted only
            ! on strict unary win-fraction, exactly as on the gridding path
            ! (refine3D_auto mirrors its gridding bootstrap; abinitio3D keeps the
            ! discrete static ladder by default, with an explicit pcg-only
            ! nu_refine=yes opt-in enabling this extension on its NU stages).
            ! The band ladder then grows over ACCEPTED candidates only, inside
            ! build_nu_evidence_state, and rides the frozen state.
            if( params%l_nu_refine )then
                call extend_nu_filter_highres_shells(vol_even, vol_odd, nsteps=nsteps_ext)
                if( nsteps_ext > 0 )then
                    write(logfhandle,'(A,I0,A,F8.3,A)') '>>> PCG NU REPLAY ('//trim(context)//&
                        &'): EVIDENCE BANK EXTENDED BY ', nsteps_ext, ' ACCEPTED SHELL STEP(S) TO ', &
                        &get_nu_filter_bank_finest_lp(), ' A'
                endif
            endif
            if( pp_plan%l_nu_envmask_incompatible )then
                write(logfhandle,'(A,1X,A)') &
                    &'>>> Existing NU evidence envelope incompatible with current box/sampling, regenerating:', &
                    &pp_plan%nu_envmask_file%to_char()
            endif
            if( l_envmask_regen )then
                call write_nu_evidence_envmask(params%nu_msk_sig, params%amsklp, vol_even%get_smpd(), &
                    &state_here, pp_plan%nu_envmask_file)
            endif
            call build_nu_evidence_state(vol_even, vol_odd, evstate)
            l_retain = .false.
            if( present(l_retain_setup) ) l_retain = l_retain_setup
            if( l_retain )then
                ! single setup serves both per-iteration NU consumers: the
                ! optimized, extended, solvent-clamped setup stays armed for
                ! the matching-reference pass, which consumes and cleans up
                call retain_nu_filter_setup(state_here)
            else
                call cleanup_nu_filter()
            endif
            ! store the frozen state for the rebuild-on-advance cadence,
            ! together with the matching low-pass it hands off: when the
            ! FSC crossing later reaches this band, the alignment search is
            ! band-limited by the ridden evidence and the cache must not be
            ! allowed to seal the stall (PfCRT regression record 2026-08-31)
            handoff_lp = nu_evidence_finest_supported_lp(evstate, 0.)
            nu_evidence_cache(state_here)%valid        = .true.
            nu_evidence_cache(state_here)%box          = ldim_here(1)
            nu_evidence_cache(state_here)%fsc_find     = fsc_find
            nu_evidence_cache(state_here)%handoff_find = 0
            if( handoff_lp > TINY ) nu_evidence_cache(state_here)%handoff_find = &
                &calc_fourier_index(handoff_lp, ldim_here(1), vol_even%get_smpd())
            nu_evidence_cache(state_here)%iter_built   = params%which_iter
            nu_evidence_cache(state_here)%source       = source_here
            nu_evidence_cache(state_here)%evstate      = evstate
        else
            evstate = nu_evidence_cache(state_here)%evstate
            write(logfhandle,'(A,I0,A,I0,A)') '>>> PCG NU REPLAY ('//trim(context)//'): STATE ', state_here, &
                &' RIDING FROZEN EVIDENCE (age ', params%which_iter - nu_evidence_cache(state_here)%iter_built, &
                &' iterations; no resolution advance)'
        endif
        call pp_plan%nu_envmask_file%kill
        ! readiness contract: a valid state with an inadequate null population
        ! (the observed zero-null calibration failure) must hard-error before
        ! either replay, never attach silently
        call assert_nu_evidence_replay_ready(evstate)
        if( l_rebuild ) call print_nu_evidence_summary(evstate)
        ! local resolution assignment table: the merged-reference topology
        ! generates no filter bank, so this is the only place the per-voxel
        ! low-pass distribution is reported there; printed every iteration
        ! (cached states included) so the trajectory stays continuous
        call print_nu_evidence_lowpass_histogram(evstate)
        write(logfhandle,'(A)') '    pcg_replay_prior_mode=nu_evidence'
        call expand_nu_evidence_band_weights(evstate, band_w)
        ! the active band ladder rides the frozen state; hand it to the caller
        ! for set_nu_prior so operator and evidence can never disagree
        call get_nu_evidence_summary(evstate, evsumm)
        band_limits = evsumm%band_limits
        ! finest evidenced local cutoff, for the LP-set matching handoff --
        ! the compact state replaces the retired second NU analysis of the
        ! postprocess (pcg_priors.md S6.3 evidence-state sharing). RAW finest
        ! selected cutoff (min_pct=0), NOT the assignment-support percentile:
        ! on small membrane proteins only a small core carries fine-resolution
        ! evidence, and the support gate collapsed the matching bandwidth to
        ! the FSC crossing, sealing the alignment/evidence stall (PfCRT
        ! regression record 2026-08-31). Sparse-but-real fine evidence must be
        ! allowed to pull the search band forward; over-extension merely
        ! widens the search, over-restriction deadlocks it.
        if( present(finest_lp) )then
            finest_lp = nu_evidence_finest_supported_lp(evstate, 0.)
        endif
        seconds_here = real(toc(t_evidence),dp)
        nu_evidence_bench_seconds = nu_evidence_bench_seconds + seconds_here
        write(logfhandle,'(A,I0,A,A,A,F9.3,A)') '>>> PCG NU REPLAY ('//trim(context)//'): STATE ', state_here, &
            &' EVIDENCE PHASE (', merge('rebuilt','cached ', l_rebuild), ') took ', real(seconds_here), ' s'
        if( present(evidence_seconds) ) evidence_seconds = seconds_here
    end subroutine build_nu_replay_evidence

    !> Rebuild-on-advance decision for the frozen NU replay evidence. Rebuild
    !! when there is no valid cached state, the evidence geometry or source
    !! changed, the envmask cadence requires regeneration from live margins,
    !! the FSC=0.143 crossing reached a FINER shell than the ridden evidence
    !! was built at (an advance; the oscillating tail around a plateau does
    !! not retrigger), the alignment band is the BINDING constraint (the FSC
    !! crossing has come within one shell of the matching low-pass this
    !! evidence handed off, so riding further would seal a search-bandwidth
    !! stall -- PfCRT regression record 2026-08-31), the iteration counter is
    !! non-monotone (restart or stage reset), or the age cap
    !! NU_EVIDENCE_REBUILD_MAX_LAG is hit.
    logical function nu_evidence_needs_rebuild( params, state_here, box_here, source_here, &
            &fsc_find, l_envmask_regen ) result( l_rebuild )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: state_here, box_here, fsc_find
        character(len=*),  intent(in) :: source_here
        logical,           intent(in) :: l_envmask_regen
        l_rebuild = .true.
        if( .not. allocated(nu_evidence_cache) )then
            allocate(nu_evidence_cache(params%nstates))
            return
        endif
        if( size(nu_evidence_cache) /= params%nstates )then
            deallocate(nu_evidence_cache)
            allocate(nu_evidence_cache(params%nstates))
            return
        endif
        if( state_here < 1 .or. state_here > params%nstates ) THROW_HARD('invalid state for NU evidence cache')
        if( .not. nu_evidence_cache(state_here)%valid                          ) return
        if( nu_evidence_cache(state_here)%box /= box_here                      ) return
        if( trim(nu_evidence_cache(state_here)%source) /= trim(source_here)    ) return
        if( l_envmask_regen                                                    ) return
        if( fsc_find > nu_evidence_cache(state_here)%fsc_find                  ) return
        ! binding-band test: the matching low-pass handed off by the ridden
        ! evidence caps the alignment search; once the FSC crossing sits
        ! within one shell of that cap, no advance can be measured under the
        ! frozen evidence and riding it would deadlock the stall
        if( nu_evidence_cache(state_here)%handoff_find > 0 .and. &
            &fsc_find >= nu_evidence_cache(state_here)%handoff_find - 1        ) return
        if( params%which_iter <= nu_evidence_cache(state_here)%iter_built      ) return
        if( params%which_iter - nu_evidence_cache(state_here)%iter_built >= NU_EVIDENCE_REBUILD_MAX_LAG ) return
        l_rebuild = .false.
    end function nu_evidence_needs_rebuild

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
    !! was built for this state (automsk=yes), else the spherical mskdiam
    !! support. The projected system (P H P) u = P b is what set_mask already
    !! implements; this only chooses P (pcg_priors.md dev item 5).
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
    !! references carry). Absent/startvol reference => no support, and the
    !! caller falls back to the spherical one.
    subroutine build_pcg_state_support( params, state_here, support, l_have )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        type(image_msk),   intent(inout) :: support
        logical,           intent(out)   :: l_have
        type(image)  :: vol_prev
        real         :: lp_support
        l_have = .false.
        call support%kill_bimg
        ! automsk gates using the density envelope in the solves at all;
        ! which passes consume it is decided per solve (see the l_use_support
        ! rule in reduce_solve_state_half): the regularized pass always, the
        ! unfil pass only under envfsc=yes
        if( trim(params%automsk) == 'no'          ) return
        if( params%pcg_mskfile%is_allocated()     ) return ! explicit mask wins
        if( state_here < 1 .or. state_here > size(params%vols) ) return
        if( len_trim(params%vols(state_here)%to_char()) == 0    ) return
        if( index(params%vols(state_here)%to_char(), 'startvol') > 0 ) return
        if( .not. file_exists(params%vols(state_here)) ) return
        lp_support = params%envmsklp
        call vol_prev%read_and_crop(params%vols(state_here), params%smpd, params%box_crop, params%smpd_crop)
        call support%automask3D(params, vol_prev, .false., lp_override=lp_support, l_report=.false.)
        call vol_prev%kill
        l_have = .true.
        write(logfhandle,'(A,I0,A,F6.1,A)') '>>> PCG SOLVE SUPPORT: STATE ', state_here, &
            &', conservative density envelope at ', lp_support, ' A (replaces the spherical support)'
    end subroutine build_pcg_state_support

    !> Hard activation contract for the direct NU replay (no silent fallback):
    !! a defined pcg_nu_lambda_rel must be finite and non-negative, and a
    !! POSITIVE strength is only meaningful when the euclid ML replay actually
    !! runs -- otherwise the request would be silently ignored, which R10-style
    !! explicitness forbids. Called from every PCG execution entry.
    subroutine validate_nu_replay_request( params )
        class(parameters), intent(in) :: params
        if( .not. ieee_is_finite(params%pcg_nu_lambda_rel) .or. params%pcg_nu_lambda_rel < 0.0 )then
            THROW_HARD('pcg_nu_lambda_rel must be finite and non-negative')
        endif
        if( params%pcg_nu_lambda_rel > 0.0 .and. .not. params%l_ml_reg )then
            THROW_HARD('pcg_nu_lambda_rel > 0 requires the regularized replay: objfun=euclid ml_reg=yes')
        endif
        ! setpoint contract: an explicit pcg_nu_supp_target is only meaningful
        ! when the auto-lambda controller tracks it -- an explicit lambda pins
        ! the strength outright and would silently ignore the setpoint
        if( params%pcg_nu_supp_target > 0.0 .and. .not. params%l_pcg_nu_autotarget .and. &
            &.not. params%l_pcg_nu_autolambda )then
            THROW_HARD('pcg_nu_supp_target requires the auto-lambda controller; unset pcg_nu_lambda_rel')
        endif
    end subroutine validate_nu_replay_request

    !> NU-replay firing readout of a final map: the prior energy of the ML
    !! solution against the prior energy of the unregularized base solution of
    !! the SAME half (the replay's own reference; with P_tau absent in NU mode
    !! a vanishing pcg_nu_lambda_rel reproduces the base solution, so an inert
    !! prior reads ~0% suppression). The amplitude-domain ratio
    !! sqrt(E_ML/E_base) mirrors the rms ratio of the retired solvent
    !! suppression readout; the lambda_nu factor inside the penalty cancels.
    !! Costs two full Q_NU applications (~26 padded FFTs) -- timed as
    !! diagnostic overhead, material at small iteration budgets.
    subroutine report_nu_solve_stats( pcgop, x, base_volume, context, half, supp_pct, overhead_s )
        type(reconstructor_pcg), intent(inout) :: pcgop
        real,                    intent(in)    :: x(:,:,:)
        type(image),             intent(in)    :: base_volume
        character(len=*),        intent(in)    :: context, half
        real,                    intent(out)   :: supp_pct
        real,                    intent(out)   :: overhead_s
        real, allocatable :: xbase(:,:,:)
        real :: nu_penalty, nu_penalty_base, nu_penalty_evid, nu_penalty_base_evid
        integer(timer_int_kind) :: t_stats
        t_stats = tic()
        call pcgop%get_nu_prior_stats(x, nu_penalty, penalty_evidenced=nu_penalty_evid)
        xbase = base_volume%get_rmat()
        call pcgop%get_nu_prior_stats(xbase, nu_penalty_base, penalty_evidenced=nu_penalty_base_evid)
        deallocate(xbase)
        ! the suppression readout the auto-lambda/auto-target controllers
        ! consume is restricted to the evidenced (molecular) region: the
        ! solvent clamp holds solvent at maximum weight, and total-energy
        ! suppression would let the solvent term flood the readout and drive
        ! lambda to the floor (weakening the prior where evidence grades it)
        supp_pct = 0.0
        if( nu_penalty_base_evid > 1.0e-12 ) &
            &supp_pct = 100.0 * (1.0 - sqrt(max(nu_penalty_evid,0.0) / nu_penalty_base_evid))
        overhead_s = real(real(toc(t_stats),dp))
        write(logfhandle,'(A,ES10.3,A,F6.1,A)') '>>> PCG NU REPLAY ('//trim(half)//'): lambda_eff=', &
            &pcgop%get_effective_nu_lambda(), '  evidenced suppression ', supp_pct, ' %'
        if( NU_DEV_OUTPUT )then
            write(logfhandle,'(A,ES12.4,A,ES12.4,A,F9.3)') '    pcg_nu_prior_energy_final=', &
                &nu_penalty, '  pcg_nu_prior_energy_base=', nu_penalty_base, '  stats_overhead_s=', overhead_s
            write(logfhandle,'(A,ES12.4,A,ES12.4)') '    evidenced energies final/base=', &
                &nu_penalty_evid, ' /', nu_penalty_base_evid
        endif
    end subroutine report_nu_solve_stats

    !> Persist the per-state NU-replay firing readout for the convergence
    !! reporter (simple_convergence prints it with the other iteration stats
    !! and advises on pcg_nu_lambda_rel), alongside the shipped-pair FSC=0.143
    !! crossing (the over-regularization diagnostic, never a resolution claim).
    !! The file is rewritten on every NU-replay volassemble and deleted first,
    !! so an iteration in which the replay is skipped never leaves stale
    !! values behind.
    subroutine write_nu_convergence_stats( params, supps, cnts, ship0143s )
        class(parameters), intent(in) :: params
        real,              intent(in) :: supps(:)
        integer,           intent(in) :: cnts(:)
        real,              intent(in) :: ship0143s(:)
        type(oris)   :: os_old, os
        type(string) :: key
        integer :: state, nship
        real    :: ship_avg, ship_prev
        logical :: l_ship_prev
        ! the auto-target outer loop needs one step of shipped-pair history:
        ! lift the previous iteration's average crossing before the rewrite
        ship_prev   = 0.0
        l_ship_prev = .false.
        if( file_exists(PCG_NU_STATS_FILE) )then
            call os_old%new(1, is_ptcl=.false.)
            call os_old%read(string(PCG_NU_STATS_FILE))
            if( os_old%isthere('PCG_NU_SHIP0143_AVG') )then
                ship_prev   = os_old%get(1, 'PCG_NU_SHIP0143_AVG')
                l_ship_prev = .true.
            endif
            call os_old%kill
        endif
        call del_file(PCG_NU_STATS_FILE)
        if( .not. any(cnts > 0) ) return
        ship_avg = 0.0
        nship    = 0
        do state = 1, size(supps)
            if( cnts(state) < 1 ) cycle
            ship_avg = ship_avg + ship0143s(state)
            nship    = nship + 1
        end do
        ship_avg = ship_avg / real(nship)
        call os%new(1, is_ptcl=.false.)
        call os%set(1, 'PCG_NU_LAMBDA_REL',   params%pcg_nu_lambda_rel)
        call os%set(1, 'PCG_NU_AUTOLAMBDA',   merge(1.0, 0.0, params%l_pcg_nu_autolambda))
        call os%set(1, 'PCG_NU_SUPP_TARGET',  params%pcg_nu_supp_target)
        call os%set(1, 'PCG_NU_AUTOTARGET',   merge(1.0, 0.0, params%l_pcg_nu_autotarget))
        call os%set(1, 'PCG_NU_SHIP0143_AVG', ship_avg)
        if( l_ship_prev ) call os%set(1, 'PCG_NU_SHIP0143_PREV', ship_prev)
        do state = 1, size(supps)
            if( cnts(state) < 1 ) cycle
            key = 'PCG_NU_SUPP_STATE'//int2str_pad(state,2)
            call os%set(1, key%to_char(), supps(state) / real(cnts(state)))
            key = 'PCG_NU_SHIP0143_STATE'//int2str_pad(state,2)
            call os%set(1, key%to_char(), ship0143s(state))
        end do
        call os%write(string(PCG_NU_STATS_FILE))
        call os%kill
        call key%kill
    end subroutine write_nu_convergence_stats

    !> Suppression-targeted auto-lambda: one secant step on the one-pole
    !! plant model, driven by the previous iteration's persisted readout.
    !! Memoryless by construction -- the stats file carries the lambda that
    !! produced the measurement, so continuity needs no extra state. Cold
    !! start (no stats file yet) keeps the dynamic default. Single-shot
    !! reconstructions in a fresh directory are therefore deterministic.
    !! The setpoint itself is resolved first: the AIMD outer loop when the
    !! target was left to its dynamic default, the pinned value otherwise.
    subroutine resolve_nu_autolambda( params )
        class(parameters), intent(inout) :: params
        type(oris)   :: os
        type(string) :: key
        real    :: lam_prev, lam_new, s_frac, s_t, supp_sum, ratio, target
        real    :: target_prev, target_new, ship_curr, ship_prev
        integer :: state, nsupp
        logical :: l_ship_curr, l_ship_prev, l_ship_shells_apart
        if( .not. params%l_pcg_nu_autolambda   ) return
        if( params%pcg_nu_lambda_rel <= 0.0    ) return
        if( .not. file_exists(PCG_NU_STATS_FILE) ) return
        call os%new(1, is_ptcl=.false.)
        call os%read(string(PCG_NU_STATS_FILE))
        if( .not. os%isthere('PCG_NU_LAMBDA_REL') )then
            call os%kill
            return
        endif
        lam_prev = os%get(1, 'PCG_NU_LAMBDA_REL')
        supp_sum = 0.0
        nsupp    = 0
        do state = 1, params%nstates
            key = 'PCG_NU_SUPP_STATE'//int2str_pad(state,2)
            if( os%isthere(key%to_char()) )then
                supp_sum = supp_sum + os%get(1, key%to_char())
                nsupp    = nsupp + 1
            endif
        enddo
        ! setpoint continuity and the shipped-pair trajectory for the outer loop
        target_prev = params%pcg_nu_supp_target
        if( os%isthere('PCG_NU_SUPP_TARGET') ) target_prev = os%get(1, 'PCG_NU_SUPP_TARGET')
        ship_curr   = 0.0
        ship_prev   = 0.0
        l_ship_curr = os%isthere('PCG_NU_SHIP0143_AVG')
        l_ship_prev = os%isthere('PCG_NU_SHIP0143_PREV')
        if( l_ship_curr ) ship_curr = os%get(1, 'PCG_NU_SHIP0143_AVG')
        if( l_ship_prev ) ship_prev = os%get(1, 'PCG_NU_SHIP0143_PREV')
        call os%kill
        call key%kill
        if( nsupp < 1 .or. lam_prev <= 0.0 ) return
        if( params%l_pcg_nu_autotarget )then
            ! AIMD outer loop on the setpoint, gated by the shipped-pair
            ! FSC=0.143 trajectory (crossings in Angstrom, finer = smaller):
            ! ratchet up only while the crossing improves, back off fast when
            ! it degrades, hold when it stalls (lp-limited plateaus hold).
            ! Crossings are Fourier-shell-quantized, so an adjacent-shell
            ! flip is plateau noise, not a signal: any change of less than
            ! two shells stalls regardless of the relative tolerance
            target_new = min(NU_AUTOTARGET_MAX, max(NU_AUTOTARGET_MIN, target_prev))
            if( l_ship_curr .and. l_ship_prev .and. ship_prev > 0.0 .and. ship_curr > 0.0 )then
                l_ship_shells_apart = abs(calc_fourier_index(ship_curr, params%box_crop, params%smpd_crop) - &
                    &calc_fourier_index(ship_prev, params%box_crop, params%smpd_crop)) >= 2
                if( .not. l_ship_shells_apart )then
                    write(logfhandle,'(A,F5.1)') &
                        &'>>> PCG NU AUTO-TARGET: SHIPPED PAIR STALLED; HOLDING SUPP TARGET ', target_new
                else if( ship_curr > ship_prev * (1.0 + NU_AUTOTARGET_SHIP_TOL) )then
                    target_new = max(NU_AUTOTARGET_MIN, target_prev * NU_AUTOTARGET_BACKOFF)
                    write(logfhandle,'(A,F6.2,A,F6.2,A,F5.1,A,F5.1)') &
                        &'>>> PCG NU AUTO-TARGET: SHIPPED PAIR DEGRADED (', ship_prev, ' ->', ship_curr, &
                        &' A); SUPP TARGET ', target_prev, ' ->', target_new
                else if( ship_curr < ship_prev * (1.0 - NU_AUTOTARGET_SHIP_TOL) )then
                    target_new = min(NU_AUTOTARGET_MAX, target_prev + NU_AUTOTARGET_STEP_ADD)
                    write(logfhandle,'(A,F6.2,A,F6.2,A,F5.1,A,F5.1)') &
                        &'>>> PCG NU AUTO-TARGET: SHIPPED PAIR IMPROVED (', ship_prev, ' ->', ship_curr, &
                        &' A); SUPP TARGET ', target_prev, ' ->', target_new
                else
                    write(logfhandle,'(A,F5.1)') &
                        &'>>> PCG NU AUTO-TARGET: SHIPPED PAIR STALLED; HOLDING SUPP TARGET ', target_new
                endif
            endif
            params%pcg_nu_supp_target = target_new
        endif
        target = params%pcg_nu_supp_target
        ! measured suppression as an amplitude fraction, floored and capped
        ! so the one-pole identification stays finite
        s_frac = max(NU_AUTOLAMBDA_SUPP_FLOOR, min(99.0, supp_sum / real(nsupp))) / 100.0
        s_t    = target / 100.0
        if( abs(s_frac*100.0 - target) <= NU_AUTOLAMBDA_DEADBAND_PCT )then
            ! on target: hold the previous strength (continuity, not the default)
            params%pcg_nu_lambda_rel = lam_prev
            write(logfhandle,'(A,F7.1,A,F5.1,A,F8.4)') '>>> PCG NU AUTO-LAMBDA: ON TARGET (supp ', &
                &s_frac*100.0, ' %, target ', target, ' %); HOLDING LAMBDA_REL', lam_prev
            return
        endif
        ! secant on s = g*lambda/(1+g*lambda): lambda_new = lambda_prev *
        ! [s_t(1-s)] / [s(1-s_t)], step-clamped and bounded
        ratio   = (s_t * (1.0 - s_frac)) / (s_frac * (1.0 - s_t))
        ratio   = min(NU_AUTOLAMBDA_STEP_CLAMP, max(1.0/NU_AUTOLAMBDA_STEP_CLAMP, ratio))
        lam_new = min(NU_AUTOLAMBDA_LAMBDA_MAX, max(NU_AUTOLAMBDA_LAMBDA_MIN, lam_prev * ratio))
        write(logfhandle,'(A,F7.1,A,F5.1,A,F8.4,A,F8.4)') '>>> PCG NU AUTO-LAMBDA: MEASURED supp ', &
            &s_frac*100.0, ' % (target ', target, ' %); LAMBDA_REL ', lam_prev, ' ->', lam_new
        params%pcg_nu_lambda_rel = lam_new
    end subroutine resolve_nu_autolambda

    subroutine execute_rec3D_pcg_shared( params, build, cline )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(image) :: half_even, half_odd, ml_even, ml_odd, merged
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt
        type(halfmap_diagnostics_result) :: hm_diag
        integer, allocatable :: selected_pinds(:), half_pinds(:)
        real, allocatable :: fsc(:), res0143s(:)
        real, allocatable :: ship05s(:), ship0143s(:), nu_band_w(:,:,:,:), nu_band_limits(:), nu_supps(:)
        integer, allocatable :: nu_supp_cnts(:)
        logical, allocatable :: state_written(:)
        integer :: nselected, state, n_state, n_even, n_odd, iptcl, istate
        logical :: l_nu_replay
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output, time_evidence
        logical :: l_sigma_loaded

        call validate_supported_mode()
        call validate_nu_replay_request(params)
        call resolve_nu_autolambda(params)
        nu_evidence_bench_seconds = 0.0_dp
        ! replay precision mode: pcg_nu_lambda_rel > 0 selects the direct
        ! NU-evidence replay (Q_NU), replacing P_tau (mode-exclusive,
        ! pcg_priors.md R10); validated above, so a positive strength here
        ! implies the euclid ML replay is active
        l_nu_replay = params%pcg_nu_lambda_rel > 0.0
        nselected = 0
        call build%spproj_field%sample4rec([params%fromp,params%top], nselected, selected_pinds)
        if( nselected < 1 ) THROW_HARD('no active particles selected for PCG reconstruct3D')

        if( params%cc_objfun == OBJFUN_EUCLID )then
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
            if( .not. l_sigma_loaded ) THROW_HARD('PCG objfun=euclid requires sigma2 files')
        endif

        allocate(res0143s(params%nstates), source=0.0)
        allocate(state_written(params%nstates), source=.false.)
        allocate(ship05s(params%nstates), ship0143s(params%nstates), source=0.0)
        allocate(nu_supps(params%nstates), source=0.0)
        allocate(nu_supp_cnts(params%nstates), source=0)
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
            time_evidence   = 0.0_dp
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
                &merged, hm_diag)
            fsc             = hm_diag%fsc
            res0143s(state) = hm_diag%res_fsc0143
            call arr2file(fsc, fname_fsc)
            fname_restxt = resolve_pcg_fsc_txt_fname(params, cline, state)
            call write_halfmap_diagnostics(hm_diag, params%box_crop, params%smpd_crop, fname_restxt)
            call fname_restxt%kill
            call hm_diag%kill
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                ! NU mode: evidence from the current base pair, frozen before
                ! either replay; with nu_refine=yes the evidence bank is
                ! extended by the shell walk (Stage 6.6), and with automsk
                ! enabled the evidence envelope is regenerated inside from
                ! the same frozen evidence (no envelope artifact is ever read)
                if( l_nu_replay ) call build_nu_replay_evidence(params, state, 'shared', &
                    &half_even, half_odd, res0143s(state), nu_band_w, nu_band_limits, &
                    &evidence_seconds=time_evidence)
                call regularize_state_half(state, 0, 'even', fsc, half_even, ml_even)
                call regularize_state_half(state, 1, 'odd',  fsc, half_odd,  ml_odd)
                if( allocated(nu_band_w)      ) deallocate(nu_band_w)
                if( allocated(nu_band_limits) ) deallocate(nu_band_limits)
                ! shipped-pair crossing: the over-regularization diagnostic
                ! (never a resolution claim)
                if( l_nu_replay )then
                    call shipped_pair_res(params, ml_even, ml_odd, ship05s(state), ship0143s(state))
                endif
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
            call write_output_diagnostics(state, 'shared', time_map_output, time_fsc_output, time_evidence)

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
        ! rewritten (or deleted) every volassemble so the convergence reporter
        ! never reads a stale NU firing readout
        call write_nu_convergence_stats(params, nu_supps, nu_supp_cnts, ship0143s)
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

        deallocate(selected_pinds, res0143s, state_written, ship05s, ship0143s, nu_supps, nu_supp_cnts)

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
            call set_pcg_solve_support(pcgop, params)
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
            real(dp) :: time_reduce, time_finalize, time_solve, time_total, time_nu_stats
            real :: prior_positive_min, prior_positive_max, prior_to_khat_l1, prior_to_khat_rms
            real :: supp_pct, nu_stats_overhead
            logical :: l_warm

            t_phase = tic()
            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            call set_pcg_solve_support(pcgop, params)
            call pcgop%begin_reduction
            fname = refine3D_pcg_raw_accum_fname(state_here, 1, params%numlen, half)
            call pcgop%add_raw_accum(fname, state_here, eo_here, 1, 1, &
                &pcg_raw_provenance(params), nptcls)
            time_reduce = real(toc(t_phase),dp)
            if( nptcls < 1 ) THROW_HARD('PCG ML replay requires a populated raw half accumulator')
            t_phase = tic()
            ! mode-exclusive replay precision (R10): Q_NU replaces P_tau; the
            ! reconstructor hard-errors if both are requested. Effective
            ! strengths are derived from the data scale in end_accum alongside
            ! the relative ridge lambda.
            if( l_nu_replay )then
                if( .not. allocated(nu_band_w) ) THROW_HARD('NU replay evidence was not constructed before the replay')
                if( .not. allocated(nu_band_limits) ) THROW_HARD('NU replay band ladder was not constructed before the replay')
                call pcgop%set_nu_prior(nu_band_w, nu_band_limits, params%pcg_nu_lambda_rel)
            else
                call pcgop%set_ml_prior(fsc_here, params%tau, params%hp)
            endif
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            call pcgop%assert_prior_attachment_mode
            time_finalize = real(toc(t_phase),dp)
            prior_npositive    = 0
            prior_positive_min = 0.0
            prior_positive_max = 0.0
            prior_to_khat_l1   = 0.0
            prior_to_khat_rms  = 0.0
            if( .not. l_nu_replay )then
                call pcgop%get_ml_prior_stats(prior_npositive, prior_positive_min, prior_positive_max, &
                    &prior_to_khat_l1, prior_to_khat_rms)
            endif
            x = base_volume%get_rmat()
            call override_ml_warm_start_from_previous(params, state_here, half, x, 'shared', l_warm)
            ! the closed-form shrinkage initial guess encodes the P_tau
            ! optimum; the NU replay has no such closed form yet and cold-starts
            ! from the base solution
            if( .not. l_warm .and. .not. l_nu_replay )then
                call regularized_ml_initial_guess(params, fsc_here, x, 'shared', half)
            endif
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call validate_solved_map(x, 'shared', state_here, half, 'ml')
            time_nu_stats = 0.0_dp
            if( l_nu_replay )then
                call report_nu_solve_stats(pcgop, x, base_volume, 'shared', half, supp_pct, nu_stats_overhead)
                time_nu_stats = real(nu_stats_overhead,dp)
                nu_supps(state_here)     = nu_supps(state_here) + supp_pct
                nu_supp_cnts(state_here) = nu_supp_cnts(state_here) + 1
            endif
            time_total = time_reduce + time_finalize + time_solve + time_nu_stats
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
                call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sigma_loaded)
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
    subroutine execute_rec3D_pcg_distributed_master( params, build, cline, trail_bootstrap_states, &
            &nu_replay_finest_lps )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        logical, optional, intent(out)  :: trail_bootstrap_states(:)
        real,    optional, intent(out)  :: nu_replay_finest_lps(:)
        type(image), target  :: half_even, half_odd, ml_even, ml_odd, merged
        type(image), target  :: previous_even, previous_odd, previous_merged
        type(image), pointer :: fsc_pair_even, fsc_pair_odd, fsc_pair_merged
        character(len=32)    :: evidence_source_here
        type(string) :: fname_even, fname_odd, fname_even_unfil, fname_odd_unfil, fname_vol, fname_fsc, raw_fname
        type(string) :: fname_restxt
        type(halfmap_diagnostics_result) :: hm_diag
        real, allocatable :: fsc(:), res0143s(:)
        real, allocatable :: realized_fractions(:), update_weights(:), ship05s(:), ship0143s(:)
        real, allocatable :: nu_band_w(:,:,:,:), nu_band_limits(:), nu_supps(:)
        integer, allocatable :: nu_supp_cnts(:)
        logical, allocatable :: state_written(:)
        character(len=256) :: provenance, chain_provenance
        integer :: state, part, eo, n_even, n_odd, iptcl, istate
        integer :: n_active_state, n_sampled_state
        type(image_msk) :: state_support_msk
        logical :: l_state_support
        logical :: l_has_updates, l_bootstrap, l_even_chain, l_odd_chain, l_nu_replay
        integer(timer_int_kind) :: t_state_phase
        real(dp) :: time_map_output, time_fsc_output, time_evidence

        call validate_pcg_common(params)
        call validate_nu_replay_request(params)
        call resolve_nu_autolambda(params)
        ! the partition workers are idle during the master-side solve, NU
        ! filtering and evidence phases: on local execution use the full
        ! allocation, restored to nthr before returning to the matching
        ! phase (never on a cluster, where the master's slot is fixed).
        ! Capped: OpenMP scaling saturates well before large core counts
        ! at these box sizes
        !$ if( trim(params%qsys_name) == 'local' ) &
        !$ &call omp_set_num_threads(min(omp_get_num_procs(), &
        !$ &max(params%nthr, min(PCG_MASTER_NTHR_CAP, max(1, params%nparts) * params%nthr))))
        nu_evidence_bench_seconds = 0.0_dp
        ! replay precision mode: Q_NU replaces P_tau (mode-exclusive,
        ! pcg_priors.md R10); same rule as the shared path, validated above
        ! so a positive strength implies the ML replay is active
        l_nu_replay = params%pcg_nu_lambda_rel > 0.0
        ! NU replay + trailing (policy 2026-08-28): the evidence pair is
        ! always the FSC pair. In trailing mode the base solves are the
        ! full-mass blended chain solutions -- the very statistics the ML
        ! replay re-reads -- so current-iteration evidence from that pair
        ! satisfies the evidence contract; the accumulator arithmetic is the
        ! test=pcg_frac_update gated path and is untouched by the prior. The
        ! bootstrap iteration (no chain yet) uses lag-one evidence from the
        ! previous shipped pair, exactly as its FSC does.
        if( present(trail_bootstrap_states) )then
            if( size(trail_bootstrap_states) /= params%nstates ) &
                &THROW_HARD('PCG trailing-bootstrap state output has invalid size')
            trail_bootstrap_states = .false.
        endif
        if( present(nu_replay_finest_lps) )then
            if( size(nu_replay_finest_lps) /= params%nstates ) &
                &THROW_HARD('PCG NU-replay finest-lp output has invalid size')
            nu_replay_finest_lps = 0.0
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
        allocate(ship05s(params%nstates), ship0143s(params%nstates), source=0.0)
        allocate(nu_supps(params%nstates), source=0.0)
        allocate(nu_supp_cnts(params%nstates), source=0)
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
            ! one support envelope per state, used by BOTH the base and the
            ! regularized solves: the mask constrains the estimator instead of
            ! post-processing it, and a single support keeps the base pair
            ! (FSC, NU evidence, B-factor reference) and the shipped ML pair
            ! on the same footing
            call build_pcg_state_support(params, state, state_support_msk, l_state_support)
            ! NOTE: the halves were briefly solved as two concurrent OpenMP
            ! sections; that regressed (both halves threw at the box-300
            ! final reconstruction) and is reverted pending a diagnosis of
            ! what in the solve path is not master-thread-safe or not
            ! affordable at twice the peak workspace
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
            time_evidence   = 0.0_dp
            if( params%l_ml_reg )then
                fname_even_unfil = refine3D_state_halfvol_fname(state, 'even', unfil=.true.)
                fname_odd_unfil  = refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.)
                t_state_phase = tic()
                call half_even%write(fname_even_unfil, del_if_exists=.true.)
                call half_odd%write(fname_odd_unfil, del_if_exists=.true.)
                time_map_output = real(toc(t_state_phase),dp)
            endif
            t_state_phase = tic()
            ! the FSC pair is selected once here and reused for the FSC, its
            ! summary, and the NU-replay evidence, keeping the
            ! evidence-pair-is-the-FSC-pair contract structural: the current
            ! base pair ordinarily (in trailing mode the blended chain solution
            ! the replay re-reads), the previous shipped pair in the trailing
            ! bootstrap (lag-one, exactly as its FSC)
            if( l_bootstrap )then
                call load_previous_state_halves(state, previous_even, previous_odd, previous_merged)
                fsc_pair_even        => previous_even
                fsc_pair_odd         => previous_odd
                fsc_pair_merged      => previous_merged
                evidence_source_here =  NU_EVIDENCE_SOURCE_PREV
            else
                fsc_pair_even        => half_even
                fsc_pair_odd         => half_odd
                fsc_pair_merged      => merged
                evidence_source_here =  NU_EVIDENCE_SOURCE_BASE
            endif
            call calculate_pcg_state_diagnostics(params, state, 'DISTRIBUTED', fsc_pair_even, &
                &fsc_pair_odd, fsc_pair_merged, hm_diag)
            fsc             = hm_diag%fsc
            res0143s(state) = hm_diag%res_fsc0143
            call arr2file(fsc, fname_fsc)
            fname_restxt = resolve_pcg_fsc_txt_fname(params, cline, state)
            call write_halfmap_diagnostics(hm_diag, params%box_crop, params%smpd_crop, fname_restxt)
            call fname_restxt%kill
            call hm_diag%kill
            time_fsc_output = real(toc(t_state_phase),dp)

            if( params%l_ml_reg )then
                if( l_nu_replay )then
                    ! the evidence pair is the FSC pair, selected once above;
                    ! its FSC=0.143 crossing steers the adaptive band ladder
                    if( present(nu_replay_finest_lps) )then
                        ! the finest-lp output signals a matching-reference
                        ! volassemble follows: retain the setup for it when
                        ! it would be rebuilt identically (same base pair, no
                        ! aux under nu_refine; single-state, as the module
                        ! holds one setup)
                        call build_nu_replay_evidence(params, state, 'distributed', fsc_pair_even, &
                            &fsc_pair_odd, res0143s(state), nu_band_w, nu_band_limits, &
                            &finest_lp=nu_replay_finest_lps(state), &
                            &evidence_source=trim(evidence_source_here), &
                            &evidence_seconds=time_evidence, &
                            &l_retain_setup=(params%nstates == 1 .and. params%l_nu_refine .and. &
                            &trim(evidence_source_here) == NU_EVIDENCE_SOURCE_BASE))
                    else
                        call build_nu_replay_evidence(params, state, 'distributed', fsc_pair_even, &
                            &fsc_pair_odd, res0143s(state), nu_band_w, nu_band_limits, &
                            &evidence_source=trim(evidence_source_here), &
                            &evidence_seconds=time_evidence)
                    endif
                endif
                ! sequential; see the note on the base pair above
                call reduce_solve_state_half(state, 0, 'even', ml_even, n_even, 'ml', fsc, half_even)
                call reduce_solve_state_half(state, 1, 'odd',  ml_odd,  n_odd,  'ml', fsc, half_odd)
                if( allocated(nu_band_w)      ) deallocate(nu_band_w)
                if( allocated(nu_band_limits) ) deallocate(nu_band_limits)
                ! shipped-pair crossing: the over-regularization diagnostic
                ! (never a resolution claim)
                if( l_nu_replay )then
                    call shipped_pair_res(params, ml_even, ml_odd, ship05s(state), ship0143s(state))
                endif
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
            call write_output_diagnostics(state, 'distributed', time_map_output, time_fsc_output, time_evidence)
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
        ! rewritten (or deleted) every volassemble so the convergence reporter
        ! never reads a stale NU firing readout
        call write_nu_convergence_stats(params, nu_supps, nu_supp_cnts, ship0143s)
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
        call state_support_msk%kill_bimg
        deallocate(res0143s, state_written, realized_fractions, update_weights, ship05s, ship0143s, &
            &nu_supps, nu_supp_cnts)
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
            real :: realized_fraction, update_weight, current_scale
            real :: supp_pct, nu_stats_overhead
            logical :: l_ml_solve, l_chain_exists, l_seed_chain, l_warm

            l_ml_solve = present(fsc_prior)
            if( l_ml_solve .neqv. present(warm_start) )then
                THROW_HARD('distributed PCG ML replay requires both FSC and warm start')
            endif

            call pcgop%new(params%box_crop, params%smpd_crop, PCG_LAMBDA)
            ! the regularized pass always takes the density-envelope support;
            ! the unfil pass takes it only under envfsc=yes, so with
            ! envfsc=no the FSC pair is estimated on the spherical support
            call set_pcg_solve_support(pcgop, params, state_support_msk, &
                &l_state_support .and. (l_ml_solve .or. params%l_envfsc))
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
                ! mode-exclusive replay precision (R10): Q_NU replaces P_tau;
                ! the reconstructor hard-errors if both are requested. The
                ! effective strengths derive from the data scale in end_accum.
                if( l_nu_replay )then
                    if( .not. allocated(nu_band_w) ) &
                        &THROW_HARD('NU replay evidence was not constructed before the replay')
                    if( .not. allocated(nu_band_limits) ) &
                        &THROW_HARD('NU replay band ladder was not constructed before the replay')
                    call pcgop%set_nu_prior(nu_band_w, nu_band_limits, params%pcg_nu_lambda_rel)
                else
                    call pcgop%set_ml_prior(fsc_prior, params%tau, params%hp)
                endif
            endif
            call pcgop%end_accum(.true.)
            call pcgop%set_op_mode(PCG_OP_KERNEL)
            if( l_ml_solve ) call pcgop%assert_prior_attachment_mode
            time_finalize = real(toc(t_phase),dp)
            prior_npositive    = 0
            prior_positive_min = 0.0
            prior_positive_max = 0.0
            prior_to_khat_l1   = 0.0
            prior_to_khat_rms  = 0.0
            if( l_ml_solve .and. .not. l_nu_replay )then
                call pcgop%get_ml_prior_stats(prior_npositive, prior_positive_min, prior_positive_max, &
                    &prior_to_khat_l1, prior_to_khat_rms)
            endif
            if( l_ml_solve )then
                x = warm_start%get_rmat()
                call override_ml_warm_start_from_previous(params, state_here, half, x, 'distributed', l_warm)
                ! the closed-form shrinkage initial guess encodes the P_tau
                ! optimum; the NU replay cold-starts from the base solution
                if( .not. l_warm .and. .not. l_nu_replay )then
                    call regularized_ml_initial_guess(params, fsc_prior, x, 'distributed', half)
                endif
            else
                allocate(x(params%box_crop,params%box_crop,params%box_crop), source=0.0)
            endif
            t_phase = tic()
            call pcgop%solve_accum(x, maxits=params%maxits_pcg, rtol=params%rtol, &
                &rel_res_hist=rel_res_hist, niters=niters, outcome=result)
            time_solve = real(toc(t_phase),dp)
            call validate_solved_map(x, 'distributed', state_here, half, solve_kind)
            if( l_ml_solve .and. l_nu_replay )then
                ! warm_start is this iteration's unregularized base half solve,
                ! the reference the firing readout is measured against
                call report_nu_solve_stats(pcgop, x, warm_start, 'distributed', half, supp_pct, nu_stats_overhead)
                ! the two halves may execute as concurrent sections
                !$omp critical(pcg_nu_supp_accum)
                nu_supps(state_here)     = nu_supps(state_here) + supp_pct
                nu_supp_cnts(state_here) = nu_supp_cnts(state_here) + 1
                !$omp end critical(pcg_nu_supp_accum)
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

    subroutine write_output_diagnostics( state, execution_mode, map_time, fsc_time, evidence_time )
        integer,            intent(in) :: state
        character(len=*),   intent(in) :: execution_mode
        real(dp),           intent(in) :: map_time, fsc_time
        real(dp), optional, intent(in) :: evidence_time
        type(string) :: fname
        integer :: funit
        fname = 'reconstruct3D_pcg_state'//int2str_pad(state,2)//'_output_diagnostics.txt'
        call fopen(funit, file=fname, status='replace', action='write')
        write(funit,'(A,A)')     'execution_mode=', trim(execution_mode)
        write(funit,'(A,F12.6)') 'halfmap_merged_output_seconds=', map_time
        write(funit,'(A,F12.6)') 'fsc_cfar_summary_seconds=', fsc_time
        if( present(evidence_time) ) write(funit,'(A,F12.6)') 'nu_evidence_phase_seconds=', evidence_time
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
