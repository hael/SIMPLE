!@descr: projection-aware covariance heterogeneity commander
module simple_commanders_flex_pca
use simple_commanders_api, only: autoscale, builder, cmdline, commander_base, file_exists, logfhandle, &
    &parameters, simple_end, simple_exception, sp_project, simple_abspath, string
implicit none

#include "simple_local_flags.inc"

! 2*smpd_crop is the working Nyquist; the 1.25 safety factor keeps the covariance band clear of it
real,    parameter :: COV_LP_OVER_NYQUIST   = 2.5
real,    parameter :: COV_LP_DEFAULT = 16.0   !< default lp (A); smpd_target = lp / COV_LP_OVER_NYQUIST; the flex_pca UI declares the same value
!> smallest working box
integer, parameter :: COV_MINBOX = 64

!> rec_backend=pcg defaults: fixed 4-iteration budget of the warm-started basis M-step (rtol=0 disarms the rtol
!! and FLEX_PCG_XTOL stops). Matches the flex_pca UI declaration. State maps are reconstructed by the
!! reconstruction service (reconstruct3D) with the run's maxits_pcg.
integer, parameter :: FLEX_PCG_MAXITS_DEFAULT = 4
real,    parameter :: FLEX_PCG_RTOL_DEFAULT   = 0.0

type, extends(commander_base) :: commander_flex_pca
  contains
    procedure :: execute => exec_flex_pca
end type commander_flex_pca

!> one part of a distributed run (simple_private_exec prg=flex_pca part=N), launched by the master
type, extends(commander_base) :: commander_flex_pca_worker
  contains
    procedure :: execute => exec_flex_pca_worker
end type commander_flex_pca_worker

contains

    subroutine exec_flex_pca( self, cline )
        use simple_flex_pca_strategy, only: flex_pca_strategy, create_flex_pca_strategy
        class(commander_flex_pca), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        class(flex_pca_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        ! defaults are the commander's; role logic is the strategy's (nparts>1 -> master, else shared
        ! memory), exactly as rec3D/refine3D; the canonical sigma2 state is settled in initialize
        call apply_flex_pca_defaults(cline)
        strategy = create_flex_pca_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_PCA NORMAL STOP ****')
    end subroutine exec_flex_pca

    !> One part: the stage the master's job description names, over the part's particle list
    subroutine exec_flex_pca_worker( self, cline )
        use simple_flex_pca_strategy, only: flex_pca_worker_strategy
        class(commander_flex_pca_worker), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(flex_pca_worker_strategy) :: strategy
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('part') ) THROW_HARD('a flex_pca worker needs part=')
        call apply_flex_pca_defaults(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_FLEX_PCA WORKER NORMAL STOP ****',print_simple=.false.)
    end subroutine exec_flex_pca_worker

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
        ! half-set EM iterations; the merge adds one joint iteration (cnga1 optimum, 2026-09-14)
        if( .not.cline%defined('n_probe_iters') ) call cline%set('n_probe_iters',4)
        call pickup_project_consensus_volume(cline)
        call derive_flex_pca_sampling(cline)
        call derive_flex_pca_band(cline)
        if( .not.cline%defined('objfun') )      call cline%set('objfun','euclid')
        call apply_flex_pca_pcg_defaults(cline)
    end subroutine apply_flex_pca_defaults

    !> rec_backend=pcg defaults: the coupled M-step on the PCG operator, warm-started from the gridding solution
    !! and capped at maxits_pcg; maxits_pcg=0 ships the gridding solution unchanged. The state maps follow
    !! rec_states_backend (default gridding), not rec_backend.
    subroutine apply_flex_pca_pcg_defaults( cline )
        class(cmdline), intent(inout) :: cline
        type(string) :: backend
        ! the masked PCG M-step is the default basis estimator; gridding stays an option
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
        if( .not.cline%defined('maxits_pcg') ) call cline%set('maxits_pcg', FLEX_PCG_MAXITS_DEFAULT)
        if( .not.cline%defined('rtol') )       call cline%set('rtol', FLEX_PCG_RTOL_DEFAULT)
    end subroutine apply_flex_pca_pcg_defaults


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

    !> Working box from a target sampling distance (refine3D's autoscale on magic boxes, floor COV_MINBOX, never
    !! finer than the data). Precedence: explicit box_crop, else lp / COV_LP_OVER_NYQUIST, else smpd_target, else
    !! COV_LP_DEFAULT / COV_LP_OVER_NYQUIST (also sets lp). Master-only, before params%new; workers inherit box_crop.
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
        ! is an expert alias and box_crop the override
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

    ! Resolve lp from project geometry; it never overrides an explicit command-line value.
    ! Called on the master before params%new so the workers inherit the resolved numbers.
    subroutine derive_flex_pca_band( cline )
        class(cmdline), intent(inout) :: cline
        type(sp_project) :: spproj
        type(string)     :: projfile
        real             :: smpd, smpd_crop, lp_here
        integer          :: box, box_crop
        if( cline%defined('lp') ) return
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
        call projfile%kill
    end subroutine derive_flex_pca_band

end module simple_commanders_flex_pca
