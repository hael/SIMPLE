!@descr: utility routines for solve2D refine2D staging and limits
module simple_solve2D_controller
use simple_core_module_api
use simple_parameters,      only: parameters
use simple_refine2D_fnames, only: refine2D_cavgs_fname
implicit none

public :: stage_params, determine_solve2D_stages, mskdiam2lplimits_refine2D, set_cline_refine2D_stage
public :: set_solve2D_sampling_policy, solve2D_seed_pass_iter, set_cline_refine2D_seed_pass
public :: SMPD_TARGET, MINBOXSZ, NSTAGES_CLS, ITS_INCR, PHASES, EXTR_LIM_LOCAL, EO_STAGE, NSAMPLE_DEFAULT_2D
public :: PROBREFINE_STAGE, STOCH_SAMPL_STAGE, STICKY_SAMPL_STAGE, FRAC_UPDATE_STAGE, CENTRING_STAGE, SHIFT_SIGMA_DEFAULT
public :: solve2D_centring_stage, solve2D_mskdiam_estimate
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
integer,          parameter :: CENTRING_STAGE     = 2                ! the polar shift map's centring stage (shift_sigma > 0)
real,             parameter :: SHIFT_SIGMA_DEFAULT = 10.             ! solve2D's shift_sigma when none is given (A); 0 = the window search
integer,          parameter :: NSAMPLE_DEFAULT_2D = 200000
integer,          parameter :: COVERAGE_NITS_CAP  = 10
character(len=3), parameter :: EO_STAGE           = 'yes'

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
    real         :: trs=0., gaufreq=0., mskdiam=0., shift_sigma=0., trs_window=0.
    logical      :: l_prob_beta_end_user=.false. !< the user set prob_beta_end: the stages keep it
    logical      :: l_anneal_tail=.false.        !< the last probabilistic stage anneals its draws to the physical posterior
end type refine2D_stage_cfg

contains

    subroutine determine_solve2D_stages( params, nstages )
        class(parameters), intent(in)  :: params
        integer,           intent(out) :: nstages
        select case(trim(params%refine))
            case('snhc','snhc_smpl','snhc_smpl_many','prob','prob_snhc')
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
        call cline_pass%delete('gaufreq')
        ! Seeded restarts skip the centring stage and keep the ordinary shift window.
        call cline_pass%delete('shift_sigma')
        call cline_pass%delete('trs_window')
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
                ! the centring stage ends when its shifts stabilise (check_conv2D's shift exit), within its budget
                if( solve2D_centring_stage(params, istage) ) cfg%minits = cfg%iter + 1
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
        ! the centring stage runs the dense table, the only fill with the map
        if( solve2D_centring_stage(params, istage) )then
            cfg%refine = 'prob'
        else if( istage < PROBREFINE_STAGE )then
            cfg%refine = 'snhc_smpl'
        else if( trim(params%refine) == 'prob' )then
            cfg%refine = 'prob'
        else if( trim(params%refine) == 'prob_snhc' .and. istage == nstages )then
            cfg%refine = 'prob'
        else
            cfg%refine = 'prob_snhc'
        endif
    end subroutine set_refine2D_stage_search_policy

    !> the stage's mask and the polar shift map's key: the centring stage's mask is the box's, the others' the
    !! stage's (the user's or the centring stage's estimate, set by the commander); the key goes to the centring
    !! stage alone (the dense table, where the map lives)
    subroutine set_refine2D_stage_mask_policy( cfg, params, stage_parms, istage )
        type(refine2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        cfg%mskdiam     = stage_parms(istage)%mskdiam
        cfg%shift_sigma = 0.
        cfg%trs_window  = 0.
        if( solve2D_centring_stage(params, istage) ) cfg%shift_sigma = params%shift_sigma
        ! the centring stage measures its posterior mass beyond the window the next stage searches
        if( solve2D_centring_stage(params, istage) ) cfg%trs_window = stage_parms(min(istage + 1, size(stage_parms)))%trslim
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
        call cline_refine2D%set('minits',     cfg%minits)
        call cline_refine2D%set('maxits',     cfg%imaxits)
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
        if( cfg%shift_sigma > 0. )then
            call cline_refine2D%set('shift_sigma', cfg%shift_sigma)
        else
            call cline_refine2D%delete('shift_sigma')
        endif
        if( cfg%mskdiam > 0. ) call cline_refine2D%set('mskdiam', cfg%mskdiam)
        if( cfg%trs_window > 0. )then
            call cline_refine2D%set('trs_window', cfg%trs_window)
        else
            call cline_refine2D%delete('trs_window')
        endif
        call cline_refine2D%set('center',     cfg%center)
        call cline_refine2D%set('restore_cavgs', 'yes')
        call cline_refine2D%set('box_crop',   stage_parms(istage)%box_crop)
        call cline_refine2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
        if( stage_parms(istage)%l_update_frac )then
            call cline_refine2D%set('update_frac', stage_parms(istage)%update_frac)
        else
            call cline_refine2D%delete('update_frac')
        endif
        ! refine2D fillin=yes requests full-assignment coverage; particle
        ! selection still follows the normal sampled-update path.
        if( l_require_full_assignment )then
            call cline_refine2D%set('fillin', 'yes')
        else
            call cline_refine2D%delete('fillin')
        endif
        call cline_refine2D%delete('endit')
    end subroutine emit_refine2D_stage_cfg

    !> the centring stage of the polar shift map (shift_sigma > 0): the first stage with references, where the table
    !! maps the disk of three sigma under the box's mask and the stage exits once the shifts stabilise
    logical function solve2D_centring_stage( params, istage )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: istage
        solve2D_centring_stage = params%l_shift_map .and. istage == CENTRING_STAGE
    end function solve2D_centring_stage

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
