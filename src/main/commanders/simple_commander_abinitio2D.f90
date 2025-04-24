! concrete commander: high-level workflows
module simple_commander_abinitio2D
include 'simple_lib.f08'
use simple_commander_base,      only: commander_base
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_exec_helpers,        only: set_shmem_flag
use simple_commander_cluster2D
use simple_commander_euclid
use simple_qsys_funs
implicit none

public :: abinitio2D_commander, autosample2D, abinitio_cleanup2D_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: abinitio2D_commander
    contains
    procedure :: execute => exec_abinitio2D
end type abinitio2D_commander

type, extends(commander_base) :: abinitio_cleanup2D_commander
    contains
    procedure :: execute => exec_abinitio_cleanup2D
end type abinitio_cleanup2D_commander

! Dimensions
real,    parameter :: SMPD_TARGET    = 2.67
integer, parameter :: MINBOXSZ       = 88
! Stages
integer, parameter :: NSTAGES        = 6
integer, parameter :: ITS_INCR       = 5
integer, parameter :: PHASES(2)      = [4, 6]
! Stages variables
real,    parameter :: ICM_LAMBDA     = 1.0
integer, parameter :: EXTR_LIM_LOCAL = 20

! convenience type
type stage_params
    real    :: lp=0., smpd_crop=0., scale=1., trslim=0.
    integer :: box_crop = 0, max_cls_pop=0, nptcls=0
    logical :: l_lpset=.false.
end type stage_params

! class variables
type(stage_params) :: stage_parms(NSTAGES)

contains

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio2D( self, cline )
        class(abinitio2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)  :: xcluster2D_distr
        type(cluster2D_commander)        :: xcluster2D
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        ! command lines
        type(cmdline)                    :: cline_cluster2D, cline_calc_pspec
        ! other
        type(parameters)                 :: params
        type(sp_project)                 :: spproj
        class(oris),             pointer :: spproj_field
        integer :: maxits(2), istage, last_iter, nptcls_eff
        logical :: l_shmem
        call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('masscen')   ) call cline%set('masscen',   'yes')
        if( .not. cline%defined('center')    ) call cline%set('center',    'yes')
        if( .not. cline%defined('sh_first')  ) call cline%set('sh_first',  'no')
        if( .not. cline%defined('cls_init')  ) call cline%set('cls_init',  'rand')
        if( .not. cline%defined('icm')       ) call cline%set('icm',       'yes')
        if( .not. cline%defined('lambda')    ) call cline%set('lambda',    ICM_LAMBDA)
        if( .not. cline%defined('extr_lim')  ) call cline%set('extr_lim',  EXTR_LIM_LOCAL)
        if( .not. cline%defined('rank_cavgs')) call cline%set('rank_cavgs','yes')
        if( .not. cline%defined('sigma_est') ) call cline%set('sigma_est', 'global')
        if( .not. cline%defined('stats')     ) call cline%set('stats',     'no')
        if( .not. cline%defined('refine')    ) call cline%set('refine',    'snhc_smpl')
        ! shared memory execution
        l_shmem = set_shmem_flag(cline)
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        maxits = [params%extr_lim, params%extr_lim+2*ITS_INCR]
        call cline%delete('stats')
        ! set downscaling
        call set_dims
        ! set resolutions limits
        call set_lplims
        ! prepare class command lines
        call prep_command_lines(cline)
        ! read project
        call spproj%read(params%projfile)
        ! sampling
        call set_sampling
        ! summary
        do istage = 1,NSTAGES
            write(logfhandle,'(A,I2,A,L1,F6.1,2I8)')'>>> STAGE ', istage,' LPSET LP MAXCLSPOP NPTCLS: ',&
            &stage_parms(istage)%l_lpset,stage_parms(istage)%lp, stage_parms(istage)%max_cls_pop,&
            &stage_parms(istage)%nptcls
        end do
        ! prep particles field
        call spproj_field%set_all2single('w',1.)
        call spproj_field%delete_2Dclustering
        if( spproj_field%get_nevenodd() == 0 ) call spproj_field%partition_eo
        call spproj%write_segment_inside(params%oritype, params%projfile)
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            if( stage_parms(istage)%l_lpset )then
                write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', stage_parms(istage)%lp
            else
                write(logfhandle,'(A,I3,A)')'>>> STAGE ', istage,' WITH GOLD STANDARD E/O'
            endif
            ! parameters update
            call set_cline_cluster2D(istage)
            ! classify
            call execute_cluster2D
        enddo
        ! transfer 2D shifts to 3D field
        call spproj%read_segment(params%oritype,params%projfile)
        call spproj%read_segment('ptcl3D',params%projfile)
        call spproj%os_ptcl3D%transfer_2Dshifts(spproj_field)
        call spproj%write_segment_inside('ptcl3D', params%projfile)
        ! weights & final mapping of particles
        if( trim(params%autosample).eq.'yes' )then
            call spproj_field%set_all2single('w',1.)
            call spproj%write_segment_inside(params%oritype, params%projfile)
            ! mapping of all particles
            if( stage_parms(NSTAGES)%nptcls < nptcls_eff )then
                call set_final_mapping
                call execute_cluster2D
                if( trim(params%stats).eq.'yes' )then
                    call spproj%read_segment(params%oritype,params%projfile)
                    call output_stats('final')
                endif
            else
                if( trim(params%stats).eq.'yes' ) call output_stats('final')
            endif
        else
            if( trim(params%stats).eq.'yes' ) call output_stats('final')
        endif
        ! final class generation & ranking
        if ( trim(params%rank_cavgs).eq.'yes' )then
            last_iter = cline_cluster2D%get_iarg('endit')
            call gen_final_cavgs(last_iter)
        endif
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        call del_file('start2Drefs_even'//params%ext)
        call del_file('start2Drefs_odd'//params%ext)
        call del_files(DIST_FBODY,      params%nparts,ext='.dat')
        call del_files(ASSIGNMENT_FBODY,params%nparts,ext='.dat')
        call del_file(DIST_FBODY      //'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
        call spproj%kill
        nullify(spproj_field)
        call qsys_cleanup
        call simple_touch(ABINITIO2D_FINISHED)
        call simple_end('**** SIMPLE_ABINITIO2D NORMAL STOP ****')
      contains

        ! Downscaling/cropping dimensions used throughout
        subroutine set_dims
            real :: smpd_target_eff, scale_factor
            if( cline%defined('box_crop') )then
                scale_factor = real(params%box_crop) / real(params%box)
                if( .not.cline%defined('smpd_crop') ) params%smpd_crop = params%smpd / scale_factor
                if( .not.cline%defined('msk_crop')  ) params%msk_crop  = round2even(params%msk * scale_factor)
                params%l_autoscale = params%box_crop < params%box
            else
                smpd_target_eff  = max(SMPD_TARGET, params%smpd)
                scale_factor     = 1.0
                params%smpd_crop = params%smpd
                params%box_crop  = params%box
                params%msk_crop  = params%msk
                if( params%l_autoscale .and. params%box >= MINBOXSZ )then
                    call autoscale(params%box, params%smpd, smpd_target_eff, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOXSZ)
                    params%l_autoscale = params%box_crop < params%box
                endif
                if( params%l_autoscale ) params%msk_crop = round2even(params%msk * scale_factor)
            endif
            if( params%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
            endif
            stage_parms(:)%box_crop  = params%box_crop
            stage_parms(:)%smpd_crop = params%smpd_crop
            stage_parms(:)%scale     = scale_factor      ! unused
            stage_parms(:)%trslim    = min(5.,max(2.0, AHELIX_WIDTH/params%smpd_crop))
        end subroutine set_dims

        ! Set resolution limits
        subroutine set_lplims
            real    :: lpstart, lpstop, cenlp
            integer :: istage
            ! Resolution limits
            call mskdiam2lplimits_here(params%mskdiam, lpstart, lpstop, cenlp)
            lpstart = max(lpstart, 2.*params%smpd_crop)
            lpstop  = max(lpstop,  2.*params%smpd_crop)
            cenlp   = max(cenlp,   2.*params%smpd_crop)
            if( .not. cline%defined('lpstart') ) params%lpstart = lpstart
            if( .not. cline%defined('lpstop')  ) params%lpstop  = lpstop
            if( .not. cline%defined('cenlp') )   params%cenlp   = cenlp
            write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', params%lpstart
            write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', params%lpstop
            write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params%cenlp
            ! Stages resolution limits
            stage_parms(1)%lp      = params%lpstart
            stage_parms(1)%l_lpset = .true.
            do istage = 2, NSTAGES-1
                stage_parms(istage)%lp      = stage_parms(istage-1)%lp - (stage_parms(istage-1)%lp - params%lpstop)/2.0
                stage_parms(istage)%l_lpset = .true.
            end do
            stage_parms(NSTAGES-1)%lp    = params%lpstop
            stage_parms(NSTAGES)%l_lpset = .false.
            stage_parms(NSTAGES)%lp      = params%lpstop
        end subroutine set_lplims

        subroutine set_sampling
            integer :: i
            nptcls_eff = spproj%count_state_gt_zero()
            stage_parms(:)%max_cls_pop = 0
            stage_parms(:)%nptcls      = nptcls_eff
            if( trim(params%autosample).eq.'yes' )then
                do i = 1,NSTAGES
                    call autosample2D(cline, nptcls_eff, params%ncls,&
                    &stage_parms(i)%max_cls_pop, stage_parms(i)%nptcls)
                    if( (params%stage > 0) .and. (params%stage < NSTAGES))then
                        if( i<=params%stage )then
                            call autosample2D(cline, nptcls_eff, params%ncls,&
                            &stage_parms(i)%max_cls_pop, stage_parms(i)%nptcls, popfac=0.5)
                        endif
                    endif
                enddo
            endif
        end subroutine set_sampling

        subroutine prep_command_lines( cline )
            use simple_default_clines, only: set_automask2D_defaults
            class(cmdline), intent(in) :: cline
            cline_cluster2D  = cline
            cline_calc_pspec = cline
            ! initial sigma2
            call cline_calc_pspec%set('prg',      'calc_pspec')
            ! cluster2D
            call cline_cluster2D%set('prg',       'cluster2D')
            call cline_cluster2D%set('wiener',    'full')
            call cline_cluster2D%set('ptclw',     'no')
            call cline_cluster2D%set('ml_reg',    'no')
            call cline_cluster2D%set('kweight',   'default')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call set_automask2D_defaults( cline_cluster2D )
        end subroutine prep_command_lines

        subroutine set_cline_cluster2D( istage )
            integer,          intent(in)  :: istage
            character(len=:), allocatable :: sh_first, refine, center, objfun, refs, icm
            integer :: iphase, iter, imaxits, cc_iters, minits, extr_iter
            real    :: trs, lambda
            select case(trim(params%refine))
            case('snhc','snhc_smpl')
                ! usual suspects
                refine = trim(params%refine)
            case('prob','prob_smpl','prob_smpl_shc','prob_greedy')
                ! prob family of algorithm
                refine = trim(params%refine)
            case DEFAULT
                THROW_HARD('Unsupported REFINE argument: '//trim(params%refine))
            end select
            ! iteration number book-keeping
            iter = 0
            if( cline_cluster2D%defined('endit') ) iter = cline_cluster2D%get_iarg('endit')
            iter = iter + 1
            call cline_cluster2D%delete('which_iter')
            ! phase logics
            if(      istage <= PHASES(1) )then
                iphase = 1
            else if( istage <= PHASES(2) )then
                iphase = 2
            else
                iphase = 0
                THROW_HARD('Invalid istage index')
            endif
            ! phase control parameters
            select case(iphase)
            case(1)
                ! phase constants
                extr_iter = 0
                ! phase variables
                imaxits   = nint(real(istage)*real(maxits(1))/real(PHASES(1)))
                minits    = imaxits
                select case(istage)
                case(1)
                    trs      = 0.
                    sh_first = 'no'
                    center   = 'no'
                    if( cline%defined('refs') )then
                        refs     = trim(params%refs)
                    else
                        refs     = NIL
                    endif
                    cc_iters     = 0
                    if( params%cc_objfun == OBJFUN_CC )then
                        objfun   = 'cc'
                    else
                        objfun   = 'euclid'
                    endif
                    if( params%l_icm )then
                        icm      = 'yes'
                        lambda   = params%lambda
                    else
                        icm      = 'no'
                    endif
                case(2)
                    trs          = stage_parms(istage)%trslim
                    sh_first     = trim(params%sh_first)
                    center       = trim(params%center)
                    refs         = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    if( params%cc_objfun == OBJFUN_CC )then
                        objfun   = 'cc'
                        cc_iters = imaxits
                    else
                        objfun   = 'euclid'
                        cc_iters = 0
                    endif
                    if( params%l_icm )then
                        icm      = 'yes'
                        lambda   = params%lambda/2.
                    else
                        icm      = 'no'
                    endif
                case(3)
                    trs          = stage_parms(istage)%trslim
                    sh_first     = trim(params%sh_first)
                    center       = trim(params%center)
                    refs         = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    cc_iters     = 0
                    objfun       = 'euclid'
                    if( params%l_icm )then
                        icm      = 'yes'
                        lambda   = params%lambda/4.
                    else
                        icm      = 'no'
                    endif
                case(4)
                    trs          = stage_parms(istage)%trslim
                    sh_first     = trim(params%sh_first)
                    center       = trim(params%center)
                    refs         = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    cc_iters     = 0
                    objfun       = 'euclid'
                    icm          = 'no'
                end select
            case(2)
                ! phase constants
                imaxits   = iter+ITS_INCR-1
                sh_first  = trim(params%sh_first)
                trs       = stage_parms(istage)%trslim
                center    = trim(params%center)
                cc_iters  = 0
                objfun    = 'euclid'
                extr_iter = params%extr_lim+1
                refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                icm       = 'no'
                minits    = iter+1
            end select
            ! command line update
            call cline_cluster2D%set('startit',   iter)
            call cline_cluster2D%set('minits',    minits)
            call cline_cluster2D%set('maxits',    imaxits)
            call cline_cluster2D%set('refine',    refine)
            call cline_cluster2D%set('objfun',    objfun)
            call cline_cluster2D%set('cc_iters',  cc_iters)
            call cline_cluster2D%set('trs',       trs)
            call cline_cluster2D%set('sh_first',  sh_first)
            call cline_cluster2D%set('center',    center)
            call cline_cluster2D%set('box_crop',  stage_parms(istage)%box_crop)
            call cline_cluster2D%set('smpd_crop', stage_parms(istage)%smpd_crop)
            if( stage_parms(istage)%l_lpset )then
                call cline_cluster2D%set('lp',    stage_parms(istage)%lp)
            else
                call cline_cluster2D%delete('lp')
            endif
            if( trim(refs) /= NIL ) call cline_cluster2D%set('refs', refs)
            if( extr_iter > 0 )then
                call cline_cluster2D%set('extr_iter', extr_iter)
            else
                call cline_cluster2D%delete('extr_iter')
            endif
            call cline_cluster2D%set('icm',    icm)
            if( trim(icm).eq.'yes' )then
                call cline_cluster2D%set('lambda', lambda)
            else
                call cline_cluster2D%delete('lambda')
            endif
            call cline_cluster2D%delete('maxpop')
            call cline_cluster2D%delete('nsample_max')
            call cline_cluster2D%delete('nsample')
            call cline_cluster2D%delete('autosample')
            call cline_cluster2D%delete('update_frac')
            if( stage_parms(istage)%max_cls_pop > 0 )then
                call cline_cluster2D%set('maxpop', stage_parms(istage)%max_cls_pop)
                if( cline%defined('update_frac') )then
                    call cline_cluster2D%set('update_frac', params%update_frac)
                endif
            endif
            call cline_cluster2D%delete('endit')
        end subroutine set_cline_cluster2D

        subroutine execute_cluster2D
            call del_file(CLUSTER2D_FINISHED)
            ! Initial sigma2
            if( istage == 1 )then
                call xcalc_pspec_distr%execute_safe(cline_calc_pspec)
            endif
            ! Classification
            if( l_shmem )then
                call xcluster2D%execute_safe(cline_cluster2D)
            else
                call xcluster2D_distr%execute_safe(cline_cluster2D)
            endif
        end subroutine execute_cluster2D

        subroutine output_stats( prefix )
            character(len=*) :: prefix
            real, allocatable :: M(:,:)
            allocate(M(nptcls_eff,2))
            M(:,1) = spproj_field%get_all('class', nonzero=.true.)
            M(:,2) = spproj_field%get_all('corr',  nonzero=.true.)
            call rmat2file(M, trim(prefix)//'_class_scores.mat')
        end subroutine output_stats

        subroutine set_final_mapping
            character(len=:), allocatable :: refs
            integer :: iter, minits
            iter      = cline_cluster2D%get_iarg('endit') + 1
            refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
            minits    = iter
            if( params%l_autoscale ) call cline_cluster2D%set('restore_cavgs', 'no')
            call cline_cluster2D%set('startit',   iter)
            call cline_cluster2D%set('minits',    minits)
            call cline_cluster2D%set('maxits',    minits)
            call cline_cluster2D%set('refine',    'greedy_smpl')
            call cline_cluster2D%set('objfun',    'euclid')
            call cline_cluster2D%set('cc_iters',  0)
            call cline_cluster2D%set('trs',       stage_parms(NSTAGES)%trslim)
            call cline_cluster2D%set('sh_first',  params%sh_first)
            call cline_cluster2D%set('center',    'no')
            call cline_cluster2D%set('box_crop',  stage_parms(NSTAGES)%box_crop)
            call cline_cluster2D%set('smpd_crop', stage_parms(NSTAGES)%smpd_crop)
            call cline_cluster2D%set('refs',      refs)
            call cline_cluster2D%set('extr_iter', params%extr_lim+1)
            call cline_cluster2D%set('icm',       'no')
            call cline_cluster2D%delete('lp')
            call cline_cluster2D%delete('lambda')
            call cline_cluster2D%delete('nsample_max')
            call cline_cluster2D%delete('nsample')
            call cline_cluster2D%delete('autosample')
            call cline_cluster2D%delete('update_frac')
            call cline_cluster2D%delete('maxpop')
            call cline_cluster2D%delete('endit')
            write(logfhandle,'(A,/,A)')'>>>','>>> FINAL MAPPING'
        end subroutine set_final_mapping

        subroutine gen_final_cavgs( iter )
            type(make_cavgs_commander_distr) :: xmake_cavgs_distr
            type(make_cavgs_commander)       :: xmake_cavgs
            type(rank_cavgs_commander)       :: xrank_cavgs
            type(cmdline)                    :: cline_make_cavgs, cline_rank_cavgs
            character(len=:),    allocatable :: finalcavgs, finalcavgs_ranked
            integer :: iter
            finalcavgs = trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//params%ext
            ! classes generation
            if( params%l_autoscale )then
                cline_make_cavgs = cline ! ncls is transferred here
                call cline_make_cavgs%delete('autoscale')
                call cline_make_cavgs%delete('balance')
                call cline_make_cavgs%delete('smpd_crop')
                call cline_make_cavgs%delete('box_crop')
                call cline_make_cavgs%set('prg',        'make_cavgs')
                call cline_make_cavgs%set('refs',       finalcavgs)
                call cline_make_cavgs%set('which_iter', iter)
                if( l_shmem )then
                    call xmake_cavgs%execute_safe(cline_make_cavgs)
                else
                    call xmake_cavgs_distr%execute_safe(cline_make_cavgs)
                endif
            endif
            ! adding cavgs & FRCs to project
            call spproj%read_segment('out', params%projfile)
            call spproj%add_frcs2os_out( trim(FRCS_FILE), 'frc2D')
            call spproj%add_cavgs2os_out(trim(finalcavgs), params%smpd, imgkind='cavg')
            call spproj%write_segment_inside('out', params%projfile)
            ! rank based on gold-standard resolution estimates
            finalcavgs_ranked = trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//'_ranked'//params%ext
            call cline_rank_cavgs%set('projfile', params%projfile)
            call cline_rank_cavgs%set('flag',     'res')
            call cline_rank_cavgs%set('stk',      finalcavgs)
            call cline_rank_cavgs%set('outstk',   finalcavgs_ranked)
            call xrank_cavgs%execute( cline_rank_cavgs )
        end subroutine gen_final_cavgs

        ! attempt at defining resolution limits
        subroutine mskdiam2lplimits_here( mskdiam, lpstart,lpstop, lpcen )
            real, intent(in)    :: mskdiam
            real, intent(inout) :: lpstart,lpstop, lpcen
            lpstart = min(max(mskdiam/12., 15.), 20.)
            lpstop  = min(max(mskdiam/22.,  6.),  8.)
            lpcen   = min(max(mskdiam/6.,  20.), 30.)
            ! was:
            ! lpstart = max(min(mskdiam/12., 15.),  8.)
            ! lpstop  = min(max(mskdiam/22.,  5.),  8.)
            ! lpcen   = min(max(mskdiam/6.,  20.), 30.)
        end subroutine mskdiam2lplimits_here

    end subroutine exec_abinitio2D

    subroutine autosample2D( cline, nptcls, ncls, max_pop, max_nptcls, popfac )
        use simple_parameters, only: params_glob
        class(cmdline), intent(in)  :: cline
        integer,        intent(in)  :: nptcls, ncls
        integer,        intent(out) :: max_pop, max_nptcls
        real, optional, intent(in)  :: popfac ! testing purpose
        max_pop    = 0
        max_nptcls = nptcls
        if( trim(params_glob%autosample).ne.'yes' ) return
        ! Class population limit
        max_pop = MAXPOP_CLS
        if( cline%defined('nsample') ) max_pop = params_glob%nsample
        ! Total population limit
        max_nptcls = min(nptcls, 2*max_pop*ncls)
        if( cline%defined('nsample_max') )then
            max_nptcls = min(params_glob%nsample_max, max_nptcls)
        else
            if( present(popfac) )then
                max_nptcls = min(nint(popfac*real(MAXPOP_PTCLS)), nptcls)
            else
                max_nptcls = min(MAXPOP_PTCLS, nptcls)
            endif
        endif
    end subroutine autosample2D

    subroutine exec_abinitio_cleanup2D( self, cline )
        class(abinitio_cleanup2D_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        ! commanders
        type(abinitio2D_commander)       :: xabinitio2D
        type(autoselect_cavgs_commander) :: xautoselect_cavgs
        ! others
        type(parameters)                 :: params
        type(sp_project)                 :: spproj, work_proj
        integer, allocatable             :: states_autosel(:), states_autosel_inv(:), tmpinds(:), states_cavgs(:), clsinds(:)
        integer, allocatable             :: pinds(:), states_prank(:), pinds_bad_ptcls(:), pinds_good_ptcls(:)
        type(class_sample),  allocatable :: clssmp(:)
        character(len=*),      parameter :: work_projfile = 'abinitio_cleanup2D_tmpproj.simple'
        type(cmdline)                    :: cline_autosel_cavgs
        integer :: nptcls
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! run first 2D
        call xabinitio2D%execute_safe(cline)
        ! read project
        call spproj%read(params%projfile)
        ! make work project
        call simple_copy_file(params%projfile, work_projfile)
        ! make automatic class selection
        call cline_autosel_cavgs%set('mskdiam',  params%mskdiam)
        call cline_autosel_cavgs%set('projfile', work_projfile)
        call cline_autosel_cavgs%set('prune',    'no')
        call xautoselect_cavgs%execute_safe(cline_autosel_cavgs)
        call work_proj%read(work_projfile)
        clsinds = work_proj%get_selected_clsinds()
        call work_proj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
        ! divide into two parts by splitting good classes into top/bottom ranked particles
        nptcls  = work_proj%get_nptcls()
        allocate(states_prank(nptcls), states_autosel(nptcls), states_autosel_inv(nptcls), source=0)
        call work_proj%os_ptcl2D%sample_ranked_parts(clssmp, 2, states_prank)
        where( states_prank.eq.1 )
            states_autosel     = 1
            states_autosel_inv = 0
        else where
            states_autosel     = 0
            states_autosel_inv = 0
        endwhere
        deallocate(clsinds)
        call simple_end('**** SIMPLE_ABINITIO_CLEANUP2D NORMAL STOP ****')
    end subroutine exec_abinitio_cleanup2D

end module simple_commander_abinitio2D
