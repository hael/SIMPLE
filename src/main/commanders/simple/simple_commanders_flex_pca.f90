!@descr: projection-aware covariance heterogeneity commander
module simple_commanders_flex_pca
use simple_commanders_api
implicit none

#include "simple_local_flags.inc"

! 2*smpd_crop is the working Nyquist; the 1.25 safety factor keeps the covariance band clear of it
real,    parameter :: COV_LP_OVER_NYQUIST   = 2.5
!> default target sampling (flex_pca_envelope_support.md 2.2): working Nyquist 4.4 A, so helical pitch
!! (5.4 A) and strand separation (4.8 A) sit inside the lattice; open question 6 of the note
real,    parameter :: COV_LP_DEFAULT = 16.0   !< default variance resolution (A): box 64 on every project of the 2026-09 campaign; the helix regime (4.5 A) waits for the memory work
!> smallest working box (the former fixed default)
integer, parameter :: COV_MINBOX = 64

!> rec_backend=pcg defaults: iteration cap and true relative-residual tolerance of the flex solves.
!! These MATCH the flex_pca UI declaration, so a UI-driven launch and a command line without the keys
!! run the same estimator (they did not before: the UI said 20/1e-3 while this said 2/0).
!!
!! Production's refine3D/abinitio3D budget of 2 does NOT transfer. Both flex solve kinds converge far
!! more slowly than a single-volume refinement solve, measured on PfCRT (box 100 covariance, box 300
!! states) with the cold start: the coupled M-step leaves a relative residual of 0.78 after 2
!! iterations and 0.028 after 20, and the state solves leave 0.07-0.09 after 5 and 0.033-0.036 after
!! 20. The FSC prior in the operator does not change that (0.81 vs 0.78 at 2 iterations). So the cap
!! is 20 and the two stop rules do the work: rtol on the true residual and the FLEX_PCG_XTOL=1.5e-2
!! diminishing-returns stop on dx/x, which is armed only when rtol > 0. On PfCRT neither fires before
!! the cap (dx/x is still 3-4% at 20); on an easier problem they end the solve early.
!! Cold state solves still take max(maxits_pcg, FINAL_PCG_MAXITS_FLOOR=5), which now only binds when
!! maxits_pcg=0 is passed explicitly to select the gridding estimator for the basis.
integer, parameter :: FLEX_PCG_MAXITS_DEFAULT = 4      ! 2026-09-16: warm-started; ladder {2,4,8} pending (plan section 6)
real,    parameter :: FLEX_PCG_RTOL_DEFAULT   = 0.0

type, extends(commander_base) :: commander_flex_pca
  contains
    procedure :: execute => exec_flex_pca
end type commander_flex_pca

contains

    subroutine exec_flex_pca( self, cline )
        use simple_flex_pca_strategy, only: flex_pca_strategy, create_flex_pca_strategy
        use omp_lib, only: omp_set_num_threads
        class(commander_flex_pca), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(flex_pca_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        ! defaults are the commander's; role logic is the strategy's (part= -> worker, nparts>1 ->
        ! master, else shared memory), exactly as rec3D/refine3D
        call apply_flex_pca_defaults(cline)
        strategy = create_flex_pca_strategy(cline)
        call strategy%initialize(params, build, cline)
        ! canonical sigma state (origin/master, 2026-09): validated or rebuilt from particle power
        ! before any sigma read; the master process only, workers consume the state it wrote
        if( .not. cline%defined('part') )then
            call ensure_canonical_sigma_state(params, build, cline)
            ! the fallback's decision travels to the part scripts (job_descr was captured at initialize)
            if( params%l_sigma_glob ) call strategy%set_worker_key('sigma_est', 'global')
            ! the nested calc_pspec runs its own params%new, which resets the OpenMP budget and
            ! nthr_glob to the command line's per-worker nthr; restore the strategy's decision
            ! (the master thread boost) so the master tails keep their combined budget
            call omp_set_num_threads(params%nthr)
            nthr_glob = params%nthr
        endif
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        call build%kill_general_tbox
        if( cline%defined('part') )then
            call simple_end('**** SIMPLE_FLEX_PCA WORKER NORMAL STOP ****',print_simple=.false.)
        else
            call simple_end('**** SIMPLE_FLEX_PCA NORMAL STOP ****')
        endif
    end subroutine exec_flex_pca

    !> Command-line defaults (a worker keeps mkdir=no: it runs in the master's directory and must
    !! not descend into one of its own). NOT merge('no ','yes',..): merge pads the shorter branch.
    subroutine apply_flex_pca_defaults( cline )
        class(cmdline), intent(inout) :: cline
        ! the pickup and derivations below read the project before params%new discovers one in cwd
        if( .not. cline%defined('projfile') ) THROW_HARD('projfile is required for flex_pca')
        if( .not.cline%defined('mkdir') )then
            if( cline%defined('part') )then
                call cline%set('mkdir','no')
            else
                call cline%set('mkdir','yes')
            endif
        endif
        if( .not.cline%defined('oritype') )     call cline%set('oritype','ptcl3D')
        if( .not.cline%defined('nstates') )     call cline%set('nstates',1)
        ! npreimages is a CEILING, not a target: the two-gate merge collapses indistinct states below
        ! it. Under preimage_auto the key is deliberately LEFT UNDEFINED when the user did not pin
        ! one, because run_flex_pca reads cline%defined('npreimages') to decide whether it may raise
        ! the ceiling to the auto value -- injecting a default here would make that test always true.
        if( .not.cline%defined('npreimages') .and. .not.flex_pca_auto_states(cline) ) &
            &call cline%set('npreimages',4)
        if( .not.cline%defined('neigs') )       call cline%set('neigs',10)
        ! an UPPER BOUND, not an iteration count: the probe stops itself on COV_PROBE_CONV
        if( .not.cline%defined('n_probe_iters') ) call cline%set('n_probe_iters',4)
        call pickup_project_consensus_volume(cline)
        call derive_flex_pca_sampling(cline)
        call derive_flex_pca_band(cline)
        if( .not.cline%defined('ptcl_src') )    call cline%set('ptcl_src','raw')
        if( .not.cline%defined('objfun') )      call cline%set('objfun','euclid')
        if( .not.cline%defined('outvol') )      call cline%set('outvol','flex_pca_state_001.mrc')
        call apply_flex_pca_pcg_defaults(cline)
    end subroutine apply_flex_pca_defaults

    !> rec_backend=pcg (doc/implementation_notes/flex_pca_envelope_support.md, step 1): the state maps
    !! and the coupled M-step run on reconstructor_pcg, every positive budget from a ZERO start, capped
    !! at maxits_pcg with rtol and the dx/x stop doing the work; maxits_pcg=0 is the budget-0 gate that
    !! ships the gridding solution unchanged.
    subroutine apply_flex_pca_pcg_defaults( cline )
        class(cmdline), intent(inout) :: cline
        type(string) :: backend
        ! the masked PCG M-step is the default basis estimator (2026-09-16: 10 vs 3-5 reproducible
        ! components on the 10180 selection, sharper leading modes on 10028); gridding stays an option
        if( .not.cline%defined('rec_backend') ) call cline%set('rec_backend', 'pcg')
        backend = cline%get_carg('rec_backend')
        if( trim(backend%to_char()) /= 'pcg' )then
            call backend%kill
            return
        endif
        call backend%kill
        ! the solvent mask is a user input (the support the basis is solved on), never derived here;
        ! without one the PCG solve runs on the spherical mskdiam support
        if( .not.cline%defined('pcg_mskfile') ) write(logfhandle,'(A)') '>>> FLEX_PCA rec_backend=pcg without pcg_mskfile: &
            &the basis is solved on the spherical mskdiam support (pass pcg_mskfile=<solvent mask> for the envelope)'
        if( .not.cline%defined('pcgop') )      call cline%set('pcgop','kernel')
        if( .not.cline%defined('maxits_pcg') ) call cline%set('maxits_pcg', FLEX_PCG_MAXITS_DEFAULT)
        if( .not.cline%defined('rtol') )       call cline%set('rtol', FLEX_PCG_RTOL_DEFAULT)
        if( .not.cline%defined('projrec') )    call cline%set('projrec','no')
    end subroutine apply_flex_pca_pcg_defaults


    subroutine ensure_canonical_sigma_state( params, build, cline )
        use, intrinsic :: iso_fortran_env, only: int64
        use simple_commanders_euclid, only: commander_calc_pspec
        use simple_sigma2_state, only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
        use simple_sigma2_state_file, only: sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, &
            &SIGMA2_GROUP_STACK, SIGMA2_STATE_COMMITTED
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline) :: cline_pspec
        type(string) :: state_path
        integer(int64) :: layout_digest
        integer, allocatable :: cnt(:,:)
        integer :: iptcl, ngroups, status, g, e, nempty
        logical :: found, rebuild
        character(len=STDLEN) :: message
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        ! Per-stack (group) sigma2 needs particles in BOTH halves of every stack. A 2D selection deselects
        ! whole stacks, and the canonical reduce then refuses ("empty even/odd half"). Decide that here
        ! and fall back to the pooled spectrum, so sigma_est=global never has to be typed for a subset.
        if( .not. params%l_sigma_glob )then
            ngroups = 0
            do iptcl = 1, params%nptcls
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                ngroups = max(ngroups, build%spproj_field%get_int(iptcl, 'stkind'))
            enddo
            if( ngroups > 0 )then
                allocate(cnt(0:1,ngroups), source=0)
                do iptcl = 1, params%nptcls
                    if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                    g = build%spproj_field%get_int(iptcl, 'stkind')
                    e = build%spproj_field%get_eo(iptcl)
                    if( g < 1 .or. g > ngroups .or. e < 0 .or. e > 1 ) cycle
                    cnt(e,g) = cnt(e,g) + 1
                enddo
                nempty = 0
                do g = 1, ngroups
                    if( cnt(0,g) == 0 .or. cnt(1,g) == 0 ) nempty = nempty + 1
                enddo
                if( nempty > 0 )then
                    if( sum(cnt(0,:)) == 0 .or. sum(cnt(1,:)) == 0 )then
                        write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: the project carries no even/odd assignment, so &
                            &per-stack sigma2 halves cannot be formed; using the pooled (global) noise spectrum &
                            &(sigma_est=global)'
                    else
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA SIGMA: ', nempty, ' of ', ngroups, &
                            &' sigma2 groups (stacks) have no particles in one half after the selection; &
                            &using the pooled (global) noise spectrum (sigma_est=global)'
                    endif
                    params%sigma_est   = 'global'
                    params%l_sigma_glob = .true.
                    call cline%set('sigma_est', 'global')
                endif
                deallocate(cnt)
            endif
        endif
        rebuild = .true.
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( found )then
            call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
            if( status == 0 )then
                layout_digest = sigma2_state_project_layout_digest(build%spproj, build%spproj_field)
                if( params%l_sigma_glob )then
                    call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, 1, &
                        &fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                        &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_GLOBAL, &
                        &expected_ngroups=1)
                else
                    ngroups = 0
                    do iptcl = 1, params%nptcls
                        if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                        ngroups = max(ngroups, build%spproj_field%get_int(iptcl, 'stkind'))
                    enddo
                    call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, 1, &
                        &fdim(params%box)-1, params%nptcls, layout_digest, status, message, &
                        &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=SIGMA2_GROUP_STACK, &
                        &expected_ngroups=ngroups)
                endif
                rebuild = status /= 0
            endif
        endif
        if( rebuild )then
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: initializing missing or stale canonical state from particle power'
            cline_pspec = cline
            call cline_pspec%set('prg', 'calc_pspec')
            call cline_pspec%set('mkdir', 'no')
            call cline_pspec%delete('part')
            call xcalc_pspec%execute(cline_pspec)
            call build%spproj%read_segment('projinfo', params%projfile)
            call cline_pspec%kill
        endif
        call state_path%kill
    end subroutine ensure_canonical_sigma_state

    !> Whether this run lets the data set the state count. Read straight off the cmdline because it
    !! is needed before params%new.
    logical function flex_pca_auto_states( cline ) result( auto )
        class(cmdline), intent(inout) :: cline
        type(string) :: val
        auto = .false.
        if( .not. cline%defined('preimage_auto') ) return
        val  = cline%get_carg('preimage_auto')
        auto = trim(val%to_char()) == 'yes'
        call val%kill
    end function flex_pca_auto_states

    ! Pick up the project consensus map (out segment, imgkind=vol, state 1) when vol1 is not on the
    ! command line. Called on the master before params%new so the workers inherit the resolved path.
    ! The mean is read at the native particle sampling, so a consensus map registered at a stage
    ! crop is rejected here and must be passed explicitly instead.
    subroutine pickup_project_consensus_volume( cline )
        class(cmdline), intent(inout) :: cline
        type(sp_project) :: spproj
        type(string)     :: projfile, vol1
        real             :: smpd, vol_smpd
        integer          :: box, vol_box
        if( cline%defined('vol1') ) return
        projfile = cline%get_carg('projfile')
        call spproj%read_segment('out', projfile)
        if( .not. spproj%isthere_in_osout('vol', 1) )then
            call spproj%kill
            THROW_HARD('flex_pca requires a consensus mean map: pass vol1 or register one in the project out segment')
        endif
        call spproj%get_vol('vol', 1, vol1, vol_smpd, vol_box)
        call spproj%kill
        if( .not. file_exists(vol1) ) THROW_HARD('flex_pca project consensus map does not exist: '//vol1%to_char())
        call spproj%read_segment('stk', projfile)
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        call spproj%kill
        if( vol_box /= box .or. vol_smpd <= 0. .or. abs(vol_smpd - smpd) > 1.e-6 )then
            THROW_HARD('flex_pca project consensus map must match the native particle sampling; pass vol1 explicitly')
        endif
        vol1 = simple_abspath(vol1)
        call cline%set('vol1', vol1)
        write(logfhandle,'(A)') '>>> FLEX_PCA consensus map from project: '//vol1%to_char()
        call projfile%kill
        call vol1%kill
    end subroutine pickup_project_consensus_volume

    ! Resolve lp and box_rec from project geometry; neither overrides an explicit command-line value.
    !> Step 0 of flex_pca_envelope_support.md: the working box follows a target sampling distance
    !! (refine3D's autoscale on magic boxes, floor COV_MINBOX, never finer than the data) instead of a
    !! fixed box_crop=64. Precedence: an explicit box_crop is an override (tests); else smpd_target;
    !! else an explicit lp sets smpd_target = lp / COV_LP_OVER_NYQUIST so a requested band gets the
    !! lattice that carries it; else COV_SMPD_TARGET_DEFAULT. Called on the master before params%new
    !! so the workers inherit box_crop through the command line.
    subroutine derive_flex_pca_sampling( cline )
        use simple_magic_boxes, only: autoscale
        class(cmdline), intent(inout) :: cline
        type(sp_project) :: spproj
        type(string)     :: projfile
        real    :: smpd, smpd_target, smpd_crop, scale
        integer :: box, box_crop
        character(len=:), allocatable :: src
        if( cline%defined('box_crop') )then
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA working box: explicit box_crop=', &
                &cline%get_iarg('box_crop'), ' (override; smpd_target ignored)'
            return
        endif
        projfile = cline%get_carg('projfile')
        call spproj%read_segment('stk', projfile)
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        call spproj%kill
        call projfile%kill
        if( box < 1 .or. smpd <= 0. ) return
        ! lp is the one user-facing key (the resolution to which the variance is resolved); smpd_target
        ! stays as an expert alias and box_crop as the override (2026-09-16)
        if( cline%defined('lp') )then
            smpd_target = cline%get_rarg('lp') / COV_LP_OVER_NYQUIST; src = 'lp / 2.5'
        else if( cline%defined('smpd_target') )then
            smpd_target = cline%get_rarg('smpd_target'); src = 'smpd_target'
        else
            smpd_target = COV_LP_DEFAULT / COV_LP_OVER_NYQUIST; src = 'default lp'
            call cline%set('lp', COV_LP_DEFAULT)
        endif
        if( smpd_target <= smpd )then
            ! the target is at or finer than the data: no down-sampling, native lattice
            box_crop  = box
            smpd_crop = smpd
            scale     = 1.0
        else
            call autoscale(box, smpd, smpd_target, box_crop, smpd_crop, scale, minbox=COV_MINBOX)
            box_crop = min(box_crop, box)
        endif
        call cline%set('box_crop', box_crop)
        call cline%set('smpd_target', smpd_target)
        write(logfhandle,'(A,F6.3,A,A,A,I0,A,F6.3,A,I0,A,F6.3,A)') '>>> FLEX_PCA working box from smpd_target=', &
            &smpd_target, ' A (', src, '): box ', box, ' @ ', smpd, ' A -> box_crop ', box_crop, ' @ ', &
            &real(box)/real(box_crop)*smpd, ' A'
        if( box_crop > 64 ) write(logfhandle,'(A,F5.1,A)') '>>> FLEX_PCA NOTE: M-step accumulators scale with &
            &box_crop^3; this box costs ', (real(box_crop)/64.)**3, 'x the box-64 footprint'
        call flush(logfhandle)
    end subroutine derive_flex_pca_sampling

    ! Called on the master before params%new so the workers inherit the resolved numbers.
    subroutine derive_flex_pca_band( cline )
        class(cmdline), intent(inout) :: cline
        type(sp_project) :: spproj
        type(string)     :: projfile
        real             :: smpd, smpd_crop, lp_here
        integer          :: box, box_crop
        if( cline%defined('lp') .and. cline%defined('box_rec') ) return
        projfile = cline%get_carg('projfile')
        call spproj%read_segment('stk', projfile)
        box  = spproj%get_box()
        smpd = spproj%get_smpd()
        call spproj%kill
        if( box < 1 .or. smpd <= 0. ) return
        box_crop = box
        if( cline%defined('box_crop') ) box_crop = cline%get_iarg('box_crop')
        if( box_crop < 1 ) box_crop = box
        if( .not. cline%defined('lp') )then
            smpd_crop = real(box) / real(box_crop) * smpd
            lp_here   = COV_LP_OVER_NYQUIST * smpd_crop
            call cline%set('lp', lp_here)
            write(logfhandle,'(A,F7.3,A,F6.3,A,I0,A)') '>>> FLEX_PCA derived lp=', lp_here, &
                &' A from smpd_crop=', smpd_crop, ' A (box ', box, ' -> box_crop)'
        endif
        ! the covariance is fitted at box_crop, but the state maps need not inherit that limit
        if( .not. cline%defined('box_rec') )then
            call cline%set('box_rec', box)
            write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA derived box_rec=', box, &
                &' (native box; state maps not capped at the covariance Nyquist)'
        endif
        call projfile%kill
    end subroutine derive_flex_pca_band

end module simple_commanders_flex_pca
