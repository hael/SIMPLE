! concrete commander: high-level workflows
module simple_commander_abinitio2D
include 'simple_lib.f08'
use simple_commander_base,      only: commander_base
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_commander_cluster2D
use simple_commander_euclid

! use simple_euclid_sigma2
use simple_qsys_funs
implicit none

public :: abinitio2D_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: abinitio2D_commander
    contains
    procedure :: execute => exec_abinitio2D
end type abinitio2D_commander

! class constants
real,    parameter :: SMPD_TARGET    = 2.67
real,    parameter :: ICM_LAMBDA     = 0.5
integer, parameter :: NSTAGES        = 6
integer, parameter :: ITS_INCR       = 5
integer, parameter :: PHASES(2)      = [4, 6]
integer, parameter :: MINBOXSZ       = 88
integer, parameter :: EXTR_LIM_LOCAL = 20

! class variables
type(lp_crop_inf)  :: lpinfo(NSTAGES)

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
        integer :: maxits(2), istage, last_iter
        logical :: l_shmem
        call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('autoscale')) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('masscen')  ) call cline%set('masscen',   'yes')
        if( .not. cline%defined('center')   ) call cline%set('center',    'yes')
        if( .not. cline%defined('sh_first') ) call cline%set('sh_first',  'no')
        if( .not. cline%defined('cls_init') ) call cline%set('cls_init',  'rand')
        if( .not. cline%defined('icm')      ) call cline%set('icm',       'yes')
        if( .not. cline%defined('lambda')   ) call cline%set('lambda',    ICM_LAMBDA)
        if( .not. cline%defined('extr_lim') ) call cline%set('extr_lim',  EXTR_LIM_LOCAL)
        if( cline%defined('nparts') )then
            l_shmem = nint(cline%get_rarg('nparts')) == 1
        else
            l_shmem = .true.
        endif
        if( l_shmem ) call cline%delete('nparts')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        maxits = [params%extr_lim, params%extr_lim+2*ITS_INCR]
        ! set downscaling
        call set_dims
        ! set resolutions limits
        call set_lplims
        ! prepare class command lines
        call prep_command_lines(cline)
        ! read project
        call spproj%read(params%projfile)
        call spproj%update_projinfo(cline)
        call spproj%write_segment_inside('projinfo', params%projfile)
        ! prep particles field
        call spproj_field%set_all2single('w',1.)
        call spproj_field%delete_2Dclustering
        if( spproj_field%get_nevenodd() == 0 ) call spproj_field%partition_eo
        call spproj%write_segment_inside(params%oritype, params%projfile)
        call spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            if( lpinfo(istage)%l_lpset )then
                write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
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
        ! final class generation & ranking
        last_iter = nint(cline_cluster2D%get_rarg('endit'))
        call gen_final_cavgs(last_iter)
        ! cleanup
        call del_file('start2Drefs'//params%ext)
        call del_file('start2Drefs_even'//params%ext)
        call del_file('start2Drefs_odd'//params%ext)
        call spproj%kill
        nullify(spproj_field)
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO2D NORMAL STOP ****')
      contains

        ! Downscaling/cropping dimensions used throughout
        subroutine set_dims
            real :: smpd_target_eff, scale_factor
            smpd_target_eff  = max(SMPD_TARGET, params%smpd)
            scale_factor     = 1.0
            params%smpd_crop = params%smpd
            params%box_crop  = params%box
            params%msk_crop  = params%msk
            if( params%l_autoscale .and. params%box >= MINBOXSZ )then
                call autoscale(params%box, params%smpd, smpd_target_eff, params%box_crop, params%smpd_crop, scale_factor, minbox=MINBOXSZ)
                params%l_autoscale = params%box_crop < params%box
                if( params%l_autoscale )then
                    params%msk_crop = round2even(params%msk * scale_factor)
                    write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',params%box_crop
                endif
            endif
            lpinfo(:)%box_crop    = params%box_crop
            lpinfo(:)%smpd_crop   = params%smpd_crop
            lpinfo(:)%scale       = scale_factor      ! unused
            lpinfo(:)%l_autoscale = .false.           ! unused
            lpinfo(:)%trslim      = min(5.,max(2.0, AHELIX_WIDTH/params%smpd_crop))
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
            if( .not. cline%defined('lpstart') ) params%cenlp   = cenlp
            write(logfhandle,'(A,F5.1)') '>>> DID SET STARTING  LOW-PASS LIMIT (IN A) TO: ', params%lpstart
            write(logfhandle,'(A,F5.1)') '>>> DID SET HARD      LOW-PASS LIMIT (IN A) TO: ', params%lpstop
            write(logfhandle,'(A,F5.1)') '>>> DID SET CENTERING LOW-PASS LIMIT (IN A) TO: ', params%cenlp
            ! Stages resolution limits
            lpinfo(1)%lp      = params%lpstart
            lpinfo(1)%l_lpset = .true.
            do istage = 2, NSTAGES-1
                lpinfo(istage)%lp      = lpinfo(istage-1)%lp - (lpinfo(istage-1)%lp - params%lpstop)/2.0
                lpinfo(istage)%l_lpset = .true.
            end do
            lpinfo(NSTAGES-1)%lp    = params%lpstop
            lpinfo(NSTAGES)%l_lpset = .false.
            lpinfo(NSTAGES)%lp      = params%lpstop
            do istage = 1,NSTAGES
                print *, 'stage lpset lp: ', istage,lpinfo(istage)%l_lpset,lpinfo(istage)%lp
            end do
        end subroutine set_lplims

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
            call cline_cluster2D%set('sigma_est', 'group')
            call cline_cluster2D%set('ptclw',     'no')
            call cline_cluster2D%set('ml_reg',    'no')
            call cline_cluster2D%set('kweight',   'default')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call set_automask2D_defaults( cline_cluster2D )
        end subroutine prep_command_lines

        subroutine set_cline_cluster2D( istage )
            integer,          intent(in)  :: istage
            character(len=:), allocatable :: sh_first, refine, center, objfun, refs, icm
            integer :: iphase, iter, imaxits, maxits_glob, cc_iters, minits, extr_iter
            real    :: trs, snr_noise_reg, lambda
            refine = 'snhc_smpl' ! not optional
            ! iteration number bookkeeping
            iter = 0
            if( cline_cluster2D%defined('endit') ) iter = nint(cline_cluster2D%get_rarg('endit'))
            iter = iter + 1
            call cline_cluster2D%delete('which_iter')
            ! phase logics
            if(      istage <= PHASES(1) )then
                iphase = 1
            else if( istage <= PHASES(2) )then
                iphase = 2
            else 
                THROW_HARD('Invalid istage index')
            endif
            ! phase control parameters
            select case(iphase)
            case(1)
                ! phase constants
                maxits_glob   = params%extr_lim+ITS_INCR
                snr_noise_reg = params%snr_noise_reg
                extr_iter     = 0
                ! phase variables
                imaxits       = nint(real(istage)*real(maxits(1))/real(PHASES(1)))
                minits        = imaxits
                select case(istage)
                case(1)
                    trs          = 0.
                    sh_first     = 'no'
                    center       = 'no'
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
                    trs          = lpinfo(istage)%trslim
                    sh_first     = trim(params%sh_first)
                    center       = trim(params%center)
                    refs         = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    cc_iters     = imaxits
                    if( params%cc_objfun == OBJFUN_CC )then
                        objfun   = 'cc'
                    else
                        objfun   = 'euclid'
                    endif
                    if( params%l_icm )then
                        icm      = 'yes'
                        lambda   = params%lambda/2.
                    else
                        icm      = 'no'
                    endif
                case(3)
                    trs          = lpinfo(istage)%trslim
                    sh_first     = trim(params%sh_first)
                    center       = trim(params%center)
                    refs         = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    cc_iters      = 0
                    objfun        = 'euclid'
                    if( params%l_icm )then
                        icm      = 'yes'
                        lambda   = params%lambda/4.
                    else
                        icm      = 'no'
                    endif
                case(4)
                    trs           = lpinfo(istage)%trslim
                    sh_first      = trim(params%sh_first)
                    center        = trim(params%center)
                    refs          = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                    cc_iters      = 0
                    objfun        = 'euclid'
                    icm           = 'no'
                end select
            case(2)
                ! phase constants
                imaxits           = iter+ITS_INCR-1
                sh_first          = trim(params%sh_first)
                trs               = lpinfo(istage)%trslim
                center            = trim(params%center)
                cc_iters          = 0
                objfun            = 'euclid'
                extr_iter         = params%extr_lim+1
                refs              = trim(CAVGS_ITER_FBODY)//int2str_pad(iter-1,3)//params%ext
                icm               = 'no'
                ! phase variables
                select case(istage)
                case(5)
                    minits        = iter + 1
                    maxits_glob   = params%extr_lim+ITS_INCR
                    snr_noise_reg = params%snr_noise_reg
                case(6)
                    minits        = iter
                    maxits_glob   = 0
                    snr_noise_reg = 0.
                end select
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
            call cline_cluster2D%set('box_crop',  lpinfo(istage)%box_crop)
            call cline_cluster2D%set('smpd_crop', lpinfo(istage)%smpd_crop)
            if( lpinfo(istage)%l_lpset )then
                call cline_cluster2D%set('lp',    lpinfo(istage)%lp)
            else
                call cline_cluster2D%delete('lp')
            endif
            if( trim(refs) /= NIL ) call cline_cluster2D%set('refs', refs)
            if( params%l_noise_reg .and. snr_noise_reg > 0.001)then
                call cline_cluster2D%set('snr_noise_reg', snr_noise_reg)
                call cline_cluster2D%set('maxits_glob',   maxits_glob)
            else
                call cline_cluster2D%delete('snr_noise_reg')
                call cline_cluster2D%delete('maxits_glob')
            endif
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
            call cline_cluster2D%delete('endit')
        end subroutine set_cline_cluster2D

        subroutine execute_cluster2D
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

end module simple_commander_abinitio2D
