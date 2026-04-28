!@descr: utility routines for ab initio 2D cluster2D staging and limits
module simple_abinitio2D_controller
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

public :: stage_params, determine_abinitio2D_stages, mskdiam2lplimits_cluster2D, set_cline_cluster2D_stage
public :: set_abinitio2D_sampling_policy
public :: SMPD_TARGET, MINBOXSZ, NSTAGES_CLS, ITS_INCR, PHASES, EXTR_LIM_LOCAL, EO_STAGE, NPTCLS2SAMPLE_2D
public :: PROBREFINE_STAGE, STOCH_SAMPL_STAGE, STICKY_SAMPL_STAGE, FRAC_UPDATE_STAGE
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
integer,          parameter :: PROBREFINE_STAGE   = 5
integer,          parameter :: STOCH_SAMPL_STAGE  = PROBREFINE_STAGE ! switch from sticky to stochastic sampling when prob starts
integer,          parameter :: STICKY_SAMPL_STAGE = 1               ! sticky random subset stage
integer,          parameter :: FRAC_UPDATE_STAGE  = 2                ! fractional class-average carry-over starts here
integer,          parameter :: NPTCLS2SAMPLE_2D   = 200000
character(len=3), parameter :: EO_STAGE           = 'yes'

! convenience type
type stage_params
    real    :: lp=0., smpd_crop=0., trslim=0.
    real    :: update_frac=1.
    integer :: box_crop = 0, max_cls_pop=0, nptcls=0
    logical :: l_lpset=.false., l_update_frac=.false.
    logical :: l_sticky_sampling=.false., l_frac_restore=.false.
end type stage_params

type cluster2D_stage_cfg
    type(string) :: refine, center, objfun, refs, gauref, ml_reg
    integer      :: iphase=0, iter=0, imaxits=0, minits=0, extr_iter=0
    real         :: trs=0., gaufreq=0.
end type cluster2D_stage_cfg

contains

    subroutine determine_abinitio2D_stages( params, nstages )
        class(parameters), intent(in)  :: params
        integer,           intent(out) :: nstages
        select case(trim(params%refine))
            case('snhc','snhc_smpl','prob')
                nstages = NSTAGES_CLS
            case DEFAULT
                THROW_HARD('Unsupported REFINE argument: '//trim(params%refine))
        end select
        if( trim(params%eo_stage).ne.'yes' ) nstages = nstages-1
    end subroutine determine_abinitio2D_stages

    subroutine mskdiam2lplimits_cluster2D( mskdiam, lpstart, lpstop, lpcen )
        real, intent(in)    :: mskdiam
        real, intent(inout) :: lpstart, lpstop, lpcen
        lpstart = min(max(mskdiam/12., 15.), 20.)
        lpstop  = min(max(mskdiam/22.,  6.),  8.)
        lpcen   = min(max(mskdiam/6.,  20.), 30.)
    end subroutine mskdiam2lplimits_cluster2D

    subroutine set_cline_cluster2D_stage( cline_cluster2D, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        class(cmdline),              intent(inout) :: cline_cluster2D
        class(cmdline),              intent(in)    :: cline
        class(parameters),           intent(in)    :: params
        type(stage_params),          intent(in)    :: stage_parms(:)
        integer,                     intent(in)    :: maxits, istage
        type(cluster2D_stage_cfg) :: cfg
        call build_cluster2D_stage_cfg( cfg, cline_cluster2D, cline, params, stage_parms, maxits, istage )
        call emit_cluster2D_stage_cfg( cline_cluster2D, cfg, stage_parms, istage )
    end subroutine set_cline_cluster2D_stage

    subroutine set_abinitio2D_sampling_policy( params, stage_parms, nstages, nptcls_eff, nsample_target_2D )
        class(parameters),  intent(in)    :: params
        type(stage_params), intent(inout) :: stage_parms(:)
        integer,            intent(in)    :: nstages, nptcls_eff
        integer,            intent(out)   :: nsample_target_2D
        if( params%nsample > 0 )then
            nsample_target_2D = params%nsample
        else
            nsample_target_2D = NPTCLS2SAMPLE_2D
        endif
        if( nsample_target_2D < 1 ) THROW_HARD('nsample must be >= 1 for abinitio2D sampled update')
        stage_parms(:)%max_cls_pop = 0
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
    end subroutine set_abinitio2D_sampling_policy

    subroutine build_cluster2D_stage_cfg( cfg, cline_cluster2D, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline_cluster2D
        class(cmdline),            intent(in)    :: cline
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: maxits, istage
        call set_cluster2D_stage_objfun_policy( cfg, params )
        call set_cluster2D_stage_iteration_policy( cfg, cline_cluster2D )
        call set_cluster2D_stage_phase_policy( cfg, istage )
        call set_cluster2D_stage_reference_policy( cfg, cline, params, stage_parms, maxits, istage )
        call set_cluster2D_stage_search_policy( cfg, params, istage )
    end subroutine build_cluster2D_stage_cfg

    subroutine set_cluster2D_stage_objfun_policy( cfg, params )
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        if( params%cc_objfun == OBJFUN_CC )then
            cfg%objfun = 'cc'
        else
            cfg%objfun = 'euclid'
        endif
    end subroutine set_cluster2D_stage_objfun_policy

    subroutine set_cluster2D_stage_iteration_policy( cfg, cline_cluster2D )
        use simple_cmdline, only: cmdline
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline_cluster2D
        cfg%iter = 0
        if( cline_cluster2D%defined('endit') ) cfg%iter = cline_cluster2D%get_iarg('endit')
        cfg%iter = cfg%iter + 1
    end subroutine set_cluster2D_stage_iteration_policy

    subroutine set_cluster2D_stage_phase_policy( cfg, istage )
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        integer,                   intent(in)    :: istage
        if(      istage <= PHASES(1) )then
            cfg%iphase = 1
        else if( istage <= PHASES(2) )then
            cfg%iphase = 2
        else
            THROW_HARD('Invalid istage index')
        endif
    end subroutine set_cluster2D_stage_phase_policy

    subroutine set_cluster2D_stage_reference_policy( cfg, cline, params, stage_parms, maxits, istage )
        use simple_cmdline, only: cmdline
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(cmdline),            intent(in)    :: cline
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: maxits, istage
        logical      :: l_gaufreq_input
        l_gaufreq_input = cline%defined('gaufreq')
        select case(cfg%iphase)
            case(1)
                cfg%extr_iter = 0
                cfg%imaxits = nint(real(istage)*real(maxits)/real(PHASES(1)))
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
                        call set_cluster2D_stage_regular_refs( cfg, params, stage_parms, istage )
                end select
            case(2)
                cfg%imaxits   = cfg%iter + params%nits_per_stage - 1
                cfg%trs       = stage_parms(istage)%trslim
                cfg%center    = trim(params%center)
                cfg%extr_iter = params%extr_lim+1
                cfg%refs      = CAVGS_ITER_FBODY//int2str_pad(cfg%iter-1,3)//params%ext%to_char()
                cfg%ml_reg    = params%ml_reg
                cfg%gauref    = 'no'
                cfg%minits    = cfg%iter + 1
        end select
    end subroutine set_cluster2D_stage_reference_policy

    subroutine set_cluster2D_stage_regular_refs( cfg, params, stage_parms, istage )
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        cfg%trs      = stage_parms(istage)%trslim
        cfg%center   = trim(params%center)
        cfg%refs     = CAVGS_ITER_FBODY//int2str_pad(cfg%iter-1,3)//params%ext%to_char()
        cfg%ml_reg   = params%ml_reg
        cfg%gauref   = 'no'
    end subroutine set_cluster2D_stage_regular_refs

    subroutine set_cluster2D_stage_search_policy( cfg, params, istage )
        type(cluster2D_stage_cfg), intent(inout) :: cfg
        class(parameters),         intent(in)    :: params
        integer,                   intent(in)    :: istage
        cfg%refine = trim(params%refine)
        if( istage < STOCH_SAMPL_STAGE )then
            cfg%refine = 'snhc_smpl'
        else
            cfg%refine = 'prob'
        endif
    end subroutine set_cluster2D_stage_search_policy

    subroutine emit_cluster2D_stage_cfg( cline_cluster2D, cfg, stage_parms, istage )
        use simple_cmdline, only: cmdline
        class(cmdline),            intent(inout) :: cline_cluster2D
        type(cluster2D_stage_cfg), intent(in)    :: cfg
        type(stage_params),        intent(in)    :: stage_parms(:)
        integer,                   intent(in)    :: istage
        call cline_cluster2D%delete('which_iter')
        if( stage_parms(istage)%l_lpset )then
            call cline_cluster2D%set('lp', stage_parms(istage)%lp)
        else
            call cline_cluster2D%delete('lp')
        endif
        if( cfg%refs .ne. NIL ) call cline_cluster2D%set('refs', cfg%refs)
        if( cfg%extr_iter > 0 )then
            call cline_cluster2D%set('extr_iter', cfg%extr_iter)
        else
            call cline_cluster2D%delete('extr_iter')
        endif
        call cline_cluster2D%set('ml_reg', cfg%ml_reg)
        call cline_cluster2D%set('gauref', cfg%gauref)
        if( cfg%gauref.eq.'yes' )then
            call cline_cluster2D%set('gaufreq', cfg%gaufreq)
        else
            call cline_cluster2D%delete('gaufreq')
        endif
        call cline_cluster2D%set('minits',    cfg%minits)
        call cline_cluster2D%set('maxits',    cfg%imaxits)
        call cline_cluster2D%set('startit',   cfg%iter)
        call cline_cluster2D%set('refine',    cfg%refine)
        call cline_cluster2D%set('objfun',    cfg%objfun)
        call cline_cluster2D%set('trs',       cfg%trs)
        call cline_cluster2D%set('center',    cfg%center)
        call cline_cluster2D%set('box_crop',  stage_parms(istage)%box_crop)
        call cline_cluster2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
        if( stage_parms(istage)%l_update_frac )then
            call cline_cluster2D%set('update_frac', stage_parms(istage)%update_frac)
        else
            call cline_cluster2D%delete('update_frac')
        endif
        call cline_cluster2D%delete('endit')
    end subroutine emit_cluster2D_stage_cfg

end module simple_abinitio2D_controller
