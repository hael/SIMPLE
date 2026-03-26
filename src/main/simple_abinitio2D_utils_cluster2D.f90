!@descr: utility routines for ab initio 2D cluster2D staging and limits
module simple_abinitio2D_utils_cluster2D
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

public :: stage_params, determine_abinitio2D_stages, mskdiam2lplimits_cluster2D, set_cline_cluster2D_stage
public :: SMPD_TARGET, MINBOXSZ, NSTAGES_CLS, NSTAGES_SINGLE, ITS_INCR_SINGLE, ITS_INCR, PHASES, EXTR_LIM_LOCAL
private
#include "simple_local_flags.inc"

! Dimensions
real,    parameter :: SMPD_TARGET      = 2.67
integer, parameter :: MINBOXSZ         = 88
! Stages
integer, parameter :: NSTAGES_CLS      = 6
integer, parameter :: NSTAGES_SINGLE   = 1
integer, parameter :: ITS_INCR_SINGLE  = 4
integer, parameter :: ITS_INCR         = 5
integer, parameter :: PHASES(2)        = [4, 6]
integer, parameter :: EXTR_LIM_LOCAL   = 20
integer, parameter :: PROBREFINE_STAGE = 5

! convenience type
type stage_params
    real    :: lp=0., smpd_crop=0., trslim=0.
    integer :: box_crop = 0, max_cls_pop=0, nptcls=0
    logical :: l_lpset=.false.
end type stage_params

contains

    subroutine determine_abinitio2D_stages( refine, nstages, l_inpl )
        character(len=*), intent(in)  :: refine
        integer,          intent(out) :: nstages
        logical,          intent(out) :: l_inpl
        l_inpl = .false.
        select case(trim(refine))
            case('inpl','inpl_smpl')
                nstages = NSTAGES_SINGLE
                l_inpl  = .true.
            case('prob')
                nstages = NSTAGES_CLS
            case('snhc','snhc_smpl')
                nstages = NSTAGES_CLS
            case DEFAULT
                THROW_HARD('Unsupported REFINE argument: '//trim(refine))
        end select
    end subroutine determine_abinitio2D_stages

    subroutine mskdiam2lplimits_cluster2D( mskdiam, lpstart, lpstop, lpcen )
        real, intent(in)    :: mskdiam
        real, intent(inout) :: lpstart, lpstop, lpcen
        lpstart = min(max(mskdiam/12., 15.), 20.)
        lpstop  = min(max(mskdiam/22.,  6.),  8.)
        lpcen   = min(max(mskdiam/6.,  20.), 30.)
    end subroutine mskdiam2lplimits_cluster2D

    subroutine set_cline_cluster2D_stage( cline_cluster2D, cline, params, stage_parms, nstages, maxits, istage )
        use simple_cmdline, only: cmdline
        class(cmdline),              intent(inout) :: cline_cluster2D
        class(cmdline),              intent(in)    :: cline
        class(parameters),           intent(in)    :: params
        type(stage_params),          intent(in)    :: stage_parms(:)
        integer,                     intent(in)    :: nstages, maxits, istage
        type(string) :: refine, center, objfun, refs, gauref, ml_reg
        integer      :: iphase, iter, imaxits, minits, extr_iter
        real         :: trs, gaufreq
        logical      :: l_gaufreq_input
        refine          = trim(params%refine)
        l_gaufreq_input = cline%defined('gaufreq')
        if( params%cc_objfun == OBJFUN_CC )then
            objfun = 'cc'
        else
            objfun = 'euclid'
        endif
        iter = 0
        if( cline_cluster2D%defined('endit') ) iter = cline_cluster2D%get_iarg('endit')
        iter = iter + 1
        call cline_cluster2D%delete('which_iter')
        if( nstages == 1 )then
            minits    = ITS_INCR_SINGLE
            imaxits   = ITS_INCR_SINGLE+2
            extr_iter = 999
            trs       = stage_parms(1)%trslim
            center    = trim(params%center)
            objfun    = 'euclid'
            gauref    = 'no'
            ml_reg    = params%ml_reg
            refs      = params%refs
            call cline_cluster2D%set('lpstart', params%lpstart)
            call cline_cluster2D%set('lpstop',  params%lpstop)
            if( stage_parms(1)%l_lpset ) call cline_cluster2D%set('lp', params%lp)
        else
            if(      istage <= PHASES(1) )then
                iphase = 1
            else if( istage <= PHASES(2) )then
                iphase = 2
            else
                THROW_HARD('Invalid istage index')
            endif
            select case(iphase)
                case(1)
                    extr_iter = 0
                    imaxits = nint(real(istage)*real(maxits)/real(PHASES(1)))
                    minits  = imaxits
                    select case(istage)
                        case(1)
                            trs      = 0.
                            center   = 'no'
                            if( cline%defined('refs') )then
                                refs = params%refs
                            else
                                refs = NIL
                            endif
                            ml_reg   = 'no'
                            gauref   = 'yes'
                            if( l_gaufreq_input )then
                                gaufreq = params%gaufreq
                            else
                                gaufreq = stage_parms(istage)%lp
                            endif
                        case(2)
                            trs      = stage_parms(istage)%trslim
                            center   = trim(params%center)
                            refs     = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                            ml_reg   = params%ml_reg
                            gauref   = 'no'
                        case(3)
                            trs      = stage_parms(istage)%trslim
                            center   = trim(params%center)
                            refs     = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                            ml_reg   = params%ml_reg
                            gauref   = 'no'
                        case(4)
                            trs      = stage_parms(istage)%trslim
                            center   = trim(params%center)
                            refs     = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                            ml_reg   = params%ml_reg
                            gauref   = 'no'
                    end select
                case(2)
                    imaxits   = iter+ITS_INCR-1
                    trs       = stage_parms(istage)%trslim
                    center    = trim(params%center)
                    extr_iter = params%extr_lim+1
                    refs      = CAVGS_ITER_FBODY//int2str_pad(iter-1,3)//params%ext%to_char()
                    ml_reg    = params%ml_reg
                    gauref    = 'no'
                    minits    = iter+1
            end select
        endif
        if( stage_parms(istage)%l_lpset )then
            call cline_cluster2D%set('lp', stage_parms(istage)%lp)
        else
            call cline_cluster2D%delete('lp')
        endif
        if( refs .ne. NIL ) call cline_cluster2D%set('refs', refs)
        if( extr_iter > 0 )then
            call cline_cluster2D%set('extr_iter', extr_iter)
        else
            call cline_cluster2D%delete('extr_iter')
        endif
        call cline_cluster2D%set('ml_reg', ml_reg)
        call cline_cluster2D%set('gauref', gauref)
        if( gauref.eq.'yes' )then
            call cline_cluster2D%set('gaufreq', gaufreq)
        else
            call cline_cluster2D%delete('gaufreq')
        endif
        if( trim(params%refine) == 'prob' )then
            if( istage < PROBREFINE_STAGE )then
                refine = 'shc_smpl'
            else
                refine = 'prob'
            endif
        else
            refine = trim(params%refine)
        endif
        call cline_cluster2D%set('minits',    minits)
        call cline_cluster2D%set('maxits',    imaxits)
        call cline_cluster2D%set('startit',   iter)
        call cline_cluster2D%set('refine',    refine)
        call cline_cluster2D%set('objfun',    objfun)
        call cline_cluster2D%set('trs',       trs)
        call cline_cluster2D%set('center',    center)
        call cline_cluster2D%set('box_crop',  stage_parms(istage)%box_crop)
        call cline_cluster2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
        call cline_cluster2D%delete('update_frac')
        call cline_cluster2D%delete('endit')
    end subroutine set_cline_cluster2D_stage

end module simple_abinitio2D_utils_cluster2D
