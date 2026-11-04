!@descr: utility routines for solve2D refine2D staging and limits
module simple_solve2D_controller
use simple_core_module_api
use simple_parameters,      only: parameters
use simple_refine2D_fnames, only: refine2D_cavgs_fname
implicit none

public :: stage_params, determine_solve2D_stages, mskdiam2lplimits_refine2D, set_cline_refine2D_stage
public :: set_solve2D_sampling_policy, solve2D_seed_pass_iter, set_cline_refine2D_seed_pass
public :: SMPD_TARGET, MINBOXSZ, NSTAGES_CLS, ITS_INCR, PHASES, EXTR_LIM_LOCAL, EO_STAGE, NSAMPLE_DEFAULT_2D
public :: PROBREFINE_STAGE, STOCH_SAMPL_STAGE, STICKY_SAMPL_STAGE, FRAC_UPDATE_STAGE, FIRST_SHIFT_STAGE
public :: solve2D_first_shift_stage, solve2D_mskdiam_estimate
public :: frc05_ptcl_quantile, solve2D_class_frc05_quantile, solve2D_stage_behind_frc05
public :: SKIP_MIN_STAGE, FRC05_SKIP_PTCL_FRAC
private
#include "simple_local_flags.inc"

! Dimensions
real,             parameter :: SMPD_TARGET        = 2.67
integer,          parameter :: MINBOXSZ           = 88
! Stages
integer,          parameter :: NSTAGES_CLS        = 6
integer,          parameter :: ITS_INCR           = 5
integer,          parameter :: PHASES(2)          = [4, 6]
integer,          parameter :: EXTR_LIM_LOCAL     = 20
integer,          parameter :: PROBREFINE_STAGE   = 3
integer,          parameter :: STOCH_SAMPL_STAGE  = PROBREFINE_STAGE ! switch from sticky to stochastic sampling when prob starts
integer,          parameter :: STICKY_SAMPL_STAGE = 1                ! sticky random subset stage
integer,          parameter :: FRAC_UPDATE_STAGE  = 2                ! fractional class-average carry-over starts here
integer,          parameter :: FIRST_SHIFT_STAGE  = 2                ! the first stage with shifts: the dense table under the key, the box's mask for an automatic run's mask estimate
integer,          parameter :: NSAMPLE_DEFAULT_2D = 200000
integer,          parameter :: COVERAGE_NITS_CAP  = 10
character(len=3), parameter :: EO_STAGE           = 'yes'
! Stage omission: a stage whose planned limit the class averages already resolve past is omitted
integer,          parameter :: SKIP_MIN_STAGE       = 2   ! never the first two stages (the second is the centring stage)
real,             parameter :: FRC05_SKIP_PTCL_FRAC = 0.5 ! the fraction of the particles whose classes must resolve past the limit at FRC=0.5

! convenience type
type stage_params
    real    :: lp=0., smpd_crop=0., trslim=0.
    real    :: update_frac=1.
    real    :: mskdiam=0.       !< the stage's mask diameter (A; 0: the command line's), set under the polar shift map
    integer :: box_crop = 0, nptcls=0
    logical :: l_lpset=.false., l_update_frac=.false.
    logical :: l_sticky_sampling=.false., l_frac_restore=.false.
end type stage_params

type refine2D_stage_cfg
    type(string) :: refine, center, objfun, refs, gauref, ml_reg
    integer      :: iphase=0, iter=0, imaxits=0, minits=0, extr_iter=0
    real         :: trs=0., gaufreq=0., mskdiam=0., shift_sigma=0.
    logical      :: l_prob_beta_end_user=.false. !< the user set prob_beta_end: the stages keep it
    logical      :: l_anneal_tail=.false.        !< the last probabilistic stage anneals its draws to the physical posterior
end type refine2D_stage_cfg

contains

    subroutine determine_solve2D_stages( params, nstages )
        class(parameters), intent(in)  :: params
        integer,           intent(out) :: nstages
        select case(trim(params%refine))
            case('snhc','snhc_smpl','snhc_smpl_many','prob','prob_neigh')
                nstages = NSTAGES_CLS
            case DEFAULT
                THROW_HARD('Unsupported REFINE argument: '//trim(params%refine))
        end select
        if( trim(params%eo_stage).ne.'yes' ) nstages = nstages-1
    end subroutine determine_solve2D_stages

    subroutine mskdiam2lplimits_refine2D( mskdiam, lpstart, lpstop, lpcen )
        real, intent(in)    :: mskdiam
        real, intent(inout) :: lpstart, lpstop, lpcen
        lpstart = min(max(mskdiam/12., 15.), 20.)
        lpstop  = min(max(mskdiam/22.,  6.),  8.)
        lpcen   = min(max(mskdiam/6.,  20.), 30.)
    end subroutine mskdiam2lplimits_refine2D

    subroutine set_cline_refine2D_stage( cline_refine2D, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        class(cmdline),              intent(inout) :: cline_refine2D
        class(cmdline),              intent(in)    :: cline
        class(parameters),           intent(in)    :: params
        type(stage_params),          intent(in)    :: stage_parms(:)
        integer,                     intent(in)    :: maxits, istage
        type(refine2D_stage_cfg) :: cfg
        call build_refine2D_stage_cfg( cfg, cline_refine2D, cline, params, stage_parms, maxits, istage )
        call emit_refine2D_stage_cfg( cline_refine2D, cfg, stage_parms, istage )
    end subroutine set_cline_refine2D_stage

    subroutine set_solve2D_sampling_policy( params, stage_parms, nstages, nptcls_eff, nsample_target_2D )
        class(parameters),  intent(in)    :: params
        type(stage_params), intent(inout) :: stage_parms(:)
        integer,            intent(in)    :: nstages, nptcls_eff
        integer,            intent(out)   :: nsample_target_2D
        if( params%nsample > 0 )then
            nsample_target_2D = params%nsample
        else
            nsample_target_2D = NSAMPLE_DEFAULT_2D
        endif
        if( nsample_target_2D < 1 ) THROW_HARD('nsample must be >= 1 for solve2D sampled update')
        stage_parms(:)%nptcls      = min(nptcls_eff, nsample_target_2D)
        if( nptcls_eff > 0 )then
            stage_parms(:)%update_frac   = min(1.0, real(stage_parms(1)%nptcls) / real(nptcls_eff))
            stage_parms(:)%l_update_frac = stage_parms(:)%update_frac < 0.99
        else
            stage_parms(:)%update_frac   = 1.0
            stage_parms(:)%l_update_frac = .false.
        endif
        stage_parms(:)%l_sticky_sampling = .false.
        stage_parms(:)%l_frac_restore    = stage_parms(:)%l_update_frac
        if( nstages >= STICKY_SAMPL_STAGE )then
            stage_parms(STICKY_SAMPL_STAGE:min(STOCH_SAMPL_STAGE-1,nstages))%l_sticky_sampling = &
                stage_parms(STICKY_SAMPL_STAGE:min(STOCH_SAMPL_STAGE-1,nstages))%l_update_frac
        endif
        if( FRAC_UPDATE_STAGE > 1 )then
            stage_parms(1:min(FRAC_UPDATE_STAGE-1,nstages))%l_frac_restore = .false.
        endif
    end subroutine set_solve2D_sampling_policy

    !> Iteration index of the seed pass of a cls_init=prev run: the iteration at
    !! which the last pre-probabilistic stage would have ended, so that the first
    !! probabilistic stage (PROBREFINE_STAGE) runs exactly the iterations it runs
    !! on an unseeded run.
    integer function solve2D_seed_pass_iter( maxits, update_frac ) result( it_pass )
        integer, intent(in) :: maxits
        real,    intent(in) :: update_frac
        it_pass = phase1_stage_endit(maxits, PROBREFINE_STAGE - 1, update_frac)
        it_pass = max(2, it_pass) ! never a fresh start: the seed partition must survive
    end function solve2D_seed_pass_iter

    !> Command line of the seed pass: one dense probabilistic iteration of every
    !! active particle against the seed references, at the low-pass limit of
    !! PROBREFINE_STAGE, without extremal neighbourhood, sampling or fill-in.
    !! Mirrors the terminal all-particle pass of the commander.
    subroutine set_cline_refine2D_seed_pass( cline_pass, params, stage_parms, it_pass, refs )
        use simple_cmdline, only: cmdline
        class(cmdline),     intent(inout) :: cline_pass
        class(parameters),  intent(in)    :: params
        type(stage_params), intent(in)    :: stage_parms(:)
        integer,            intent(in)    :: it_pass
        character(len=*),   intent(in)    :: refs
        integer :: istage
        istage = min(PROBREFINE_STAGE, size(stage_parms))
        call cline_pass%delete('which_iter')
        call cline_pass%delete('endit')
        call cline_pass%delete('update_frac')
        call cline_pass%delete('nsample')
        call cline_pass%delete('fillin')
        call cline_pass%delete('assign_all')
        call cline_pass%delete('gaufreq')
        ! Seeded restarts skip the centring stage and keep the ordinary shift window (the key off explicitly: the
        ! parameter's default is the key)
        call cline_pass%set('shift_sigma', 0.)
        if( stage_parms(istage)%l_lpset )then
            call cline_pass%set('lp', stage_parms(istage)%lp)
        else
            call cline_pass%delete('lp')
        endif
        call cline_pass%set('refs',          trim(refs))
        call cline_pass%set('startit',       it_pass)
        call cline_pass%set('minits',        1)
        call cline_pass%set('maxits',        1)
        call cline_pass%set('extr_iter',     params%extr_lim + 1)
        call cline_pass%set('refine',        'prob')
        if( params%cc_objfun == OBJFUN_CC )then
            call cline_pass%set('objfun', 'cc')
        else
            call cline_pass%set('objfun', 'euclid')
        endif
        call cline_pass%set('trs',           stage_parms(istage)%trslim)
        call cline_pass%set('center',        trim(params%center))
        call cline_pass%set('ml_reg',        params%ml_reg)
        call cline_pass%set('filt_mode',     'none')
        call cline_pass%set('gauref',        'no')
        call cline_pass%set('restore_cavgs', 'yes')
        call cline_pass%set('box_crop',      stage_parms(istage)%box_crop)
        call cline_pass%set('smpd_crop',     stage_parms(istage)%smpd_crop)
    end subroutine set_cline_refine2D_seed_pass

    subroutine build_refine2D_stage_cfg( cfg, cline_refine2D, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline_refine2D
        class(cmdline),            intent(in)    :: cline
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: maxits, istage
        call set_refine2D_stage_objfun_policy( cfg, params )
        call set_refine2D_stage_iteration_policy( cfg, cline_refine2D )
        call set_refine2D_stage_phase_policy( cfg, istage )
        call set_refine2D_stage_reference_policy( cfg, cline, params, stage_parms, maxits, istage )
        call set_refine2D_stage_search_policy( cfg, params, istage, size(stage_parms) )
        call set_refine2D_stage_mask_policy( cfg, params, stage_parms, istage )
        ! the last stage anneals its draws from prob_beta to the physical posterior (prob_beta_end=1) once its class
        ! overlap converges (the tail of PROB_BETA_TAIL iterations, then the stage ends), as solve3D's last stage;
        ! the other stages run flat, unless the user set the end
        cfg%l_prob_beta_end_user = cline%defined('prob_beta_end')
        cfg%l_anneal_tail        = istage == size(stage_parms) .and. istage >= PROBREFINE_STAGE
    end subroutine build_refine2D_stage_cfg

    subroutine set_refine2D_stage_objfun_policy( cfg, params )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        if( params%cc_objfun == OBJFUN_CC )then
            cfg%objfun = 'cc'
        else
            cfg%objfun = 'euclid'
        endif
    end subroutine set_refine2D_stage_objfun_policy

    subroutine set_refine2D_stage_iteration_policy( cfg, cline_refine2D )
        use simple_cmdline, only: cmdline
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline_refine2D
        cfg%iter = 0
        if( cline_refine2D%defined('endit') ) cfg%iter = cline_refine2D%get_iarg('endit')
        cfg%iter = cfg%iter + 1
    end subroutine set_refine2D_stage_iteration_policy

    subroutine set_refine2D_stage_phase_policy( cfg, istage )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        integer,                   intent(in)    :: istage
        if(      istage <= PHASES(1) )then
            cfg%iphase = 1
        else if( istage <= PHASES(2) )then
            cfg%iphase = 2
        else
            THROW_HARD('Invalid istage index')
        endif
    end subroutine set_refine2D_stage_phase_policy

    subroutine set_refine2D_stage_reference_policy( cfg, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: maxits, istage
        logical      :: l_gaufreq_input
        l_gaufreq_input = cline%defined('gaufreq')
        select case(cfg%iphase)
            case(1)
                cfg%extr_iter = 0
                cfg%imaxits = phase1_stage_endit(maxits, istage, stage_parms(istage)%update_frac)
                cfg%minits  = cfg%imaxits
                select case(istage)
                    case(1)
                        cfg%trs      = 0.
                        cfg%center   = 'no'
                        if( cline%defined('refs') )then
                            cfg%refs = params%refs
                        else
                            cfg%refs = NIL
                        endif
                        cfg%ml_reg   = 'no'
                        cfg%gauref   = 'yes'
                        if( l_gaufreq_input )then
                            cfg%gaufreq = params%gaufreq
                        else
                            cfg%gaufreq = stage_parms(istage)%lp
                        endif
                    case(2,3,4)
                        call set_refine2D_stage_regular_refs( cfg, params, stage_parms, istage )
                end select
            case(2)
                cfg%imaxits   = cfg%iter + sampled_stage_nits(stage_parms(istage)%update_frac, params%nits_per_stage) - 1
                cfg%trs       = stage_parms(istage)%trslim
                cfg%center    = trim(params%center)
                cfg%extr_iter = params%extr_lim+1
                cfg%refs      = refine2D_cavgs_fname(cfg%iter-1, params%ext%to_char())
                cfg%ml_reg    = params%ml_reg
                cfg%gauref    = 'no'
                cfg%minits    = cfg%iter + 1
                if( cfg%imaxits >= cfg%iter + params%nits_per_stage ) cfg%minits = cfg%imaxits
        end select
    end subroutine set_refine2D_stage_reference_policy

    integer function phase1_stage_endit( maxits, istage, update_frac ) result( endit )
        integer, intent(in) :: maxits, istage
        real,    intent(in) :: update_frac
        integer :: default_endit, first_endit, post_default_nits, guarded_nits
        default_endit = nint(real(istage) * real(maxits) / real(PHASES(1)))
        if( istage <= STICKY_SAMPL_STAGE )then
            endit = default_endit
            return
        endif
        first_endit       = nint(real(STICKY_SAMPL_STAGE) * real(maxits) / real(PHASES(1)))
        post_default_nits = max(1, nint(real(maxits - first_endit) / real(max(1, PHASES(1) - STICKY_SAMPL_STAGE))))
        guarded_nits      = sample_coverage_nits(update_frac)
        if( guarded_nits <= post_default_nits )then
            endit = default_endit
        else
            endit = first_endit + (istage - STICKY_SAMPL_STAGE) * guarded_nits
        endif
    end function phase1_stage_endit

    integer function sampled_stage_nits( update_frac, default_nits ) result( nits )
        real,    intent(in) :: update_frac
        integer, intent(in) :: default_nits
        nits = max(default_nits, sample_coverage_nits(update_frac))
    end function sampled_stage_nits

    integer function sample_coverage_nits( update_frac ) result( nits )
        real, intent(in) :: update_frac
        integer :: coverage_stages
        real    :: nupdates_needed
        coverage_stages = max(1, PHASES(1) - STICKY_SAMPL_STAGE)
        if( update_frac <= 0. )then
            nits = COVERAGE_NITS_CAP
            return
        endif
        nupdates_needed = max(0., (1. / update_frac) - 1.)
        nits = ceiling(nupdates_needed / real(coverage_stages))
        nits = max(1, min(COVERAGE_NITS_CAP, nits))
    end function sample_coverage_nits

    subroutine set_refine2D_stage_regular_refs( cfg, params, stage_parms, istage )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        cfg%trs      = stage_parms(istage)%trslim
        cfg%center   = trim(params%center)
        cfg%refs     = refine2D_cavgs_fname(cfg%iter-1, params%ext%to_char())
        cfg%ml_reg   = params%ml_reg
        cfg%gauref   = 'no'
    end subroutine set_refine2D_stage_regular_refs

    subroutine set_refine2D_stage_search_policy( cfg, params, istage, nstages )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        integer,                   intent(in)    :: istage, nstages
        ! the first stage with shifts runs the dense table under the key: every particle not yet centred takes the disk
        if( solve2D_first_shift_stage(params, istage) )then
            cfg%refine = 'prob'
        else if( istage < PROBREFINE_STAGE )then
            cfg%refine = 'snhc_smpl'
        else if( trim(params%refine) == 'prob' )then
            cfg%refine = 'prob'
        else if( trim(params%refine) == 'prob_neigh' .and. istage == nstages )then
            cfg%refine = 'prob'
        else
            cfg%refine = 'prob_neigh'
        endif
    end subroutine set_refine2D_stage_search_policy

    !> the stage's mask and the polar shift map's key: the stage's mask (the user's, the box's for the mask estimate
    !! of an automatic run's first shift stage, or that estimate, set by the commander); the key goes to every stage
    !! from the first with shifts on (stage 1 is the SNHC window stage, which has no table)
    subroutine set_refine2D_stage_mask_policy( cfg, params, stage_parms, istage )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        cfg%mskdiam     = stage_parms(istage)%mskdiam
        cfg%shift_sigma = 0.
        if( params%l_shift_map .and. istage >= FIRST_SHIFT_STAGE .and. stage_parms(istage)%trslim > 0. )&
            &cfg%shift_sigma = params%shift_sigma
    end subroutine set_refine2D_stage_mask_policy

    logical function refine2D_stage_requires_full_assignment( cfg, stage_parms, istage ) result( l_required )
        type(refine2D_stage_cfg), intent(in) :: cfg
        type(stage_params),        intent(in) :: stage_parms(:)
        integer,                   intent(in) :: istage
        l_required = stage_parms(istage)%l_update_frac .and.&
            &(cfg%iphase == size(PHASES) .or. istage == size(stage_parms))
    end function refine2D_stage_requires_full_assignment

    subroutine emit_refine2D_stage_cfg( cline_refine2D, cfg, stage_parms, istage )
        use simple_cmdline, only: cmdline
        class(cmdline),            intent(inout) :: cline_refine2D
        type(refine2D_stage_cfg), intent(in)     :: cfg
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        logical :: l_require_full_assignment
        l_require_full_assignment = refine2D_stage_requires_full_assignment(cfg, stage_parms, istage)
        call cline_refine2D%delete('which_iter')
        if( stage_parms(istage)%l_lpset )then
            call cline_refine2D%set('lp', stage_parms(istage)%lp)
        else
            call cline_refine2D%delete('lp')
        endif
        if( cfg%refs .ne. NIL ) call cline_refine2D%set('refs', cfg%refs)
        if( cfg%extr_iter > 0 )then
            call cline_refine2D%set('extr_iter', cfg%extr_iter)
        else
            call cline_refine2D%delete('extr_iter')
        endif
        call cline_refine2D%set('ml_reg', cfg%ml_reg)
        call cline_refine2D%set('filt_mode', 'none')
        call cline_refine2D%set('gauref', cfg%gauref)
        if( cfg%gauref.eq.'yes' )then
            call cline_refine2D%set('gaufreq', cfg%gaufreq)
        else
            call cline_refine2D%delete('gaufreq')
        endif
        ! the stage's absolute iterations as refine2D's run lengths from startit (minits and maxits count from it)
        call cline_refine2D%set('minits',     max(0, cfg%minits  - cfg%iter + 1))
        call cline_refine2D%set('maxits',     max(1, cfg%imaxits - cfg%iter + 1))
        call cline_refine2D%set('startit',    cfg%iter)
        call cline_refine2D%set('refine',     cfg%refine)
        call cline_refine2D%set('objfun',     cfg%objfun)
        if( .not. cfg%l_prob_beta_end_user )then
            if( cfg%l_anneal_tail )then
                call cline_refine2D%set('prob_beta_end', 1.)
            else
                call cline_refine2D%delete('prob_beta_end')
            endif
        endif
        call cline_refine2D%set('trs',        cfg%trs)
        ! zero switches the key off: the parameter's default is the key
        call cline_refine2D%set('shift_sigma', cfg%shift_sigma)
        if( cfg%mskdiam > 0. ) call cline_refine2D%set('mskdiam', cfg%mskdiam)
        call cline_refine2D%set('center',     cfg%center)
        call cline_refine2D%set('restore_cavgs', 'yes')
        call cline_refine2D%set('box_crop',   stage_parms(istage)%box_crop)
        call cline_refine2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
        if( stage_parms(istage)%l_update_frac )then
            call cline_refine2D%set('update_frac', stage_parms(istage)%update_frac)
        else
            call cline_refine2D%delete('update_frac')
        endif
        ! the stages that need every particle assigned take the coverage guard assign_all=yes (2D's fillin=yes before
        ! 2026-11-03); fillin is 3D's periodic sampling of the particles never updated, which 2D does not have
        if( l_require_full_assignment )then
            call cline_refine2D%set('assign_all', 'yes')
        else
            call cline_refine2D%delete('assign_all')
        endif
        call cline_refine2D%delete('fillin')
        call cline_refine2D%delete('endit')
    end subroutine emit_refine2D_stage_cfg

    !> the first stage with shifts under the polar shift map's key (shift_sigma > 0): the dense table, where every
    !! particle not yet centred takes the disk; an automatic run estimates its mask from this stage's class averages
    logical function solve2D_first_shift_stage( params, istage )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        solve2D_first_shift_stage = params%l_shift_map .and. istage == FIRST_SHIFT_STAGE
    end function solve2D_first_shift_stage

    !> The FRC=0.5 resolution (A) the classes of a fraction frac of the particles reach: the classes best first by
    !! their FRC=0.5 crossing, the populations cumulated to frac of their total (the lower weighted quantile; 0.5 is
    !! the median particle's class resolution, which a few overfit or empty classes cannot move). frac_past is the
    !! population fraction in classes past lp. A class with pop < 1 is left out; 0 without a populated class.
    subroutine frc05_ptcl_quantile( res05s, pops, frac, res05_q, lp, frac_past )
        real,           intent(in)  :: res05s(:) !< FRC=0.5 resolution per class (A)
        integer,        intent(in)  :: pops(:)   !< population per class
        real,           intent(in)  :: frac
        real,           intent(out) :: res05_q
        real, optional, intent(in)  :: lp
        real, optional, intent(out) :: frac_past
        real,    allocatable :: keys(:)
        integer, allocatable :: order(:)
        integer :: ncls, i, ntot, ncum
        res05_q = 0.
        if( present(frac_past) ) frac_past = 0.
        ncls = size(res05s)
        if( ncls < 1 .or. size(pops) /= ncls ) return
        ntot = sum(pops, mask=pops > 0)
        if( ntot < 1 ) return
        if( present(lp) .and. present(frac_past) )then
            frac_past = real(sum(pops, mask=pops > 0 .and. res05s < lp)) / real(ntot)
        endif
        keys  = res05s
        order = [(i, i=1,ncls)]
        call hpsort(keys, order)
        ncum = 0
        do i = 1, ncls
            if( pops(order(i)) < 1 ) cycle
            ncum    = ncum + pops(order(i))
            res05_q = keys(i)
            if( real(ncum) >= frac * real(ntot) ) exit
        end do
        deallocate(keys, order)
    end subroutine frc05_ptcl_quantile

    !> frc05_ptcl_quantile over the project's selected populated classes (cls2D state > 0, pop > 0) and their
    !! FRCs (the out segment's frc2D, read from the file); 0 without them
    subroutine solve2D_class_frc05_quantile( spproj, frac, res05_q, lp, frac_past )
        use simple_sp_project, only: sp_project
        use simple_class_frcs, only: class_frcs
        class(sp_project), intent(inout) :: spproj
        real,              intent(in)    :: frac
        real,              intent(out)   :: res05_q
        real, optional,    intent(in)    :: lp
        real, optional,    intent(out)   :: frac_past
        type(class_frcs)     :: clsfrcs
        type(string)         :: frcs_fname
        real,    allocatable :: res05s(:)
        integer, allocatable :: pops(:)
        integer :: ncls, icls
        real    :: res05, res0143
        res05_q = 0.
        if( present(frac_past) ) frac_past = 0.
        ncls = spproj%os_cls2D%get_noris()
        if( ncls < 1 ) return
        call spproj%get_frcs(frcs_fname, 'frc2D', fail=.false.)
        if( frcs_fname%to_char() == NIL )then
            call frcs_fname%kill
            return
        endif
        if( .not. file_exists(frcs_fname) )then
            call frcs_fname%kill
            return
        endif
        call clsfrcs%read(frcs_fname)
        if( clsfrcs%get_ncls() /= ncls ) THROW_HARD('class FRCs and cls2D disagree on the class count; solve2D_class_frc05_quantile')
        allocate(res05s(ncls), source=0.)
        allocate(pops(ncls),   source=0)
        do icls = 1, ncls
            if( spproj%os_cls2D%get_state(icls) < 1 ) cycle
            pops(icls) = nint(spproj%os_cls2D%get(icls, 'pop'))
            if( pops(icls) < 1 ) cycle
            call clsfrcs%estimate_res(icls, res05, res0143)
            res05s(icls) = res05
        end do
        call frc05_ptcl_quantile(res05s, pops, frac, res05_q, lp, frac_past)
        call clsfrcs%kill
        call frcs_fname%kill
        deallocate(res05s, pops)
    end subroutine solve2D_class_frc05_quantile

    !> Whether stage istage is omitted: its planned limit is behind the FRC=0.5 resolution the classes of
    !! FRC05_SKIP_PTCL_FRAC of the particles reached at the ends of the last two executed stages (the estimate
    !! held across a stage). Never the entry stage, the first SKIP_MIN_STAGE stages, the last stage of the
    !! invocation or a stage without a planned limit (the gold-standard stage). Phase 1 ends at fixed iterations
    !! (the extremal schedule's): an omitted phase-1 stage's iterations go to the phase-1 stage that follows; the
    !! last phase-1 stage's and a phase-2 stage's are dropped.
    logical function solve2D_stage_behind_frc05( stage_parms, istage, start_stage, stop_stage, res05_hist, nhist ) result( l_behind )
        type(stage_params), intent(in) :: stage_parms(:)
        integer,            intent(in) :: istage, start_stage, stop_stage
        real,               intent(in) :: res05_hist(:) !< the quantile at the ends of the executed stages, in order
        integer,            intent(in) :: nhist         !< # stages executed so far
        l_behind = .false.
        if( istage <= start_stage )                 return
        if( istage <= SKIP_MIN_STAGE )              return
        if( istage >= stop_stage )                  return
        if( istage > size(stage_parms) )            return
        if( .not. stage_parms(istage)%l_lpset )     return
        if( nhist < 2 .or. nhist > size(res05_hist) ) return
        if( any(res05_hist(nhist-1:nhist) < TINY) ) return
        l_behind = all(res05_hist(nhist-1:nhist) < stage_parms(istage)%lp)
    end function solve2D_stage_behind_frc05

    !> the mask diameter (A) the stages after the centring stage take from its class averages: automask2D's
    !! diameters of the populated classes, their population-weighted 90th percentile, as the diameter itself (a
    !! user's mskdiam is the molecule's diameter; the picking box's padding and factor overshot by 40 %)
    subroutine solve2D_mskdiam_estimate( cavgs, cls_field, box, smpd, mskdiam )
        use simple_image,          only: image
        use simple_image_msk,      only: automask2D_settings, automask2D
        use simple_default_clines, only: AUTOMASK2D_NGROW, AUTOMASK2D_WINSZ, AUTOMASK2D_AMSKLP, AUTOMASK2D_EDGE
        class(string), intent(in)  :: cavgs
        class(oris),   intent(in)  :: cls_field
        integer,       intent(in)  :: box
        real,          intent(in)  :: smpd
        real,          intent(out) :: mskdiam
        type(automask2D_settings) :: settings
        type(image), allocatable  :: imgs(:)
        real,        allocatable  :: diams(:), shifts(:,:), pops(:)
        integer,     allocatable  :: order(:)
        integer :: ldim(3), ncls, icls
        real    :: wtot, wcum, diam_q
        call find_ldim_nptcls(cavgs, ldim, ncls)
        if( ncls /= cls_field%get_noris() ) THROW_HARD('class averages and class field disagree on the class count; solve2D_mskdiam_estimate')
        settings%amsklp   = AUTOMASK2D_AMSKLP
        settings%ngrow    = AUTOMASK2D_NGROW
        settings%winsz    = nint(AUTOMASK2D_WINSZ)
        settings%edge     = AUTOMASK2D_EDGE
        settings%l_writer = .false.
        allocate(imgs(ncls), pops(ncls), order(ncls))
        do icls = 1, ncls
            call imgs(icls)%new([ldim(1),ldim(2),1], smpd)
            call imgs(icls)%read(cavgs, icls)
            pops(icls) = 0.
            if( cls_field%get_state(icls) > 0 ) pops(icls) = max(0., cls_field%get(icls, 'pop'))
        enddo
        call automask2D(settings, imgs, diams, shifts, write2disk=.false., verbose=.false.)
        where( diams < TINY ) pops = 0. ! no component found
        wtot = sum(pops)
        if( .not. (wtot > 0.) ) THROW_HARD('no populated class to estimate the mask diameter from; solve2D_mskdiam_estimate')
        ! the population-weighted 90th percentile of the diameters
        order = [(icls, icls = 1, ncls)]
        call hpsort(diams, order)
        wcum   = 0.
        diam_q = diams(ncls)
        do icls = 1, ncls
            wcum = wcum + pops(order(icls))
            if( wcum >= 0.9 * wtot )then
                diam_q = diams(icls)
                exit
            endif
        enddo
        mskdiam = min(diam_q, (real(box) - COSMSKHALFWIDTH) * smpd)
        do icls = 1, ncls
            call imgs(icls)%kill
        enddo
        deallocate(imgs, pops, order, diams, shifts)
    end subroutine solve2D_mskdiam_estimate

end module simple_solve2D_controller
