! concrete commander: high-level workflows
module simple_commander_abinitio2D
include 'simple_lib.f08'
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_qsys_env,            only: qsys_env
use simple_commander_base,      only: commander_base
use simple_image,               only: image
use simple_class_frcs,          only: class_frcs
use simple_convergence,         only: convergence
use simple_commander_cluster2D, only: cluster2D_commander_distr
use simple_commander_euclid
use simple_euclid_sigma2
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
real,             parameter :: LPSTART_DEFAULT      = 15.
real,             parameter :: CENLP_DEFAULT        = 25.
real,             parameter :: SMPD_TARGET          = 2.67
integer,          parameter :: NSTAGES              = 6
integer,          parameter :: PHASES(2)            = [5,  6]
integer,          parameter :: MAXITS(2)            = [15,20]
integer,          parameter :: MINBOXSZ             = 88
! class variables
type(lp_crop_inf) :: lpinfo(NSTAGES)
type(cmdline)     :: cline_cluster2D, cline_calc_pspec_distr

contains

    !> for generation of an initial 3d model from particles
    subroutine exec_abinitio2D( self, cline )
        class(abinitio2D_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        ! commanders
        type(cluster2D_commander_distr)  :: xcluster2D_distr
        type(calc_pspec_commander_distr) :: xcalc_pspec_distr
        ! other
        type(parameters)     :: params
        type(sp_project)     :: spproj
        class(oris), pointer :: spproj_field 
        integer              :: istage
        call cline%set('oritype',   'ptcl2D')
        if( .not. cline%defined('autoscale')) call cline%set('autoscale', 'yes')
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('masscen')  ) call cline%set('masscen',   'yes')
        if( .not. cline%defined('center')   ) call cline%set('center',    'yes')
        if( .not. cline%defined('cls_init') ) call cline%set('cls_init', 'rand')
        ! make master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%ptr2oritype(params%oritype, spproj_field)
        ! prepare downscaling, resolution limits
        call prep_dims_lplims
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
        ! Starting references
        call prep_refs
        ! Initial sigma2
        call xcalc_pspec_distr%execute_safe(cline_calc_pspec_distr)
        ! Frequency marching
        do istage = 1, NSTAGES
            write(logfhandle,'(A)')'>>>'
            write(logfhandle,'(A,I3,A9,F5.1)')'>>> STAGE ', istage,' WITH LP =', lpinfo(istage)%lp
            ! parameters update
            call set_cline_cluster2D(istage)
            if( lpinfo(istage)%l_autoscale )then
                write(logfhandle,'(A,I3,A1,I3)')'>>> ORIGINAL/CROPPED IMAGE SIZE (pixels): ',params%box,'/',lpinfo(istage)%box_crop
            endif
            ! Search
            call cline_cluster2D%delete('endit')
            call xcluster2D_distr%execute_safe(cline_cluster2D)
        enddo
        ! TODO: final class generation, ranking, clustering
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_ABINITIO2D NORMAL STOP ****')
      contains

        ! calculates:
        ! 1. Downscaling/cropping dimensions used throughout
        ! 2. resolution limits
        subroutine prep_dims_lplims
            use simple_estimate_ssnr, only: mskdiam2lplimits
            real    :: lpstart, lpstop, cenlp, smpd_target_eff, scale_factor
            integer :: istage
            ! Downscaling, same for ALL STAGES
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
            ! Resolution limits
            call mskdiam2lplimits(params%mskdiam, lpstart, lpstop, cenlp)
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
            lpinfo(:)%box_crop    = params%box_crop
            lpinfo(:)%smpd_crop   = params%smpd_crop
            lpinfo(:)%scale       = scale_factor      ! unused
            lpinfo(:)%l_autoscale = .false.           ! unused
            lpinfo(:)%trslim      = min(8.,max(2.0, AHELIX_WIDTH/params%smpd))
            lpinfo(1)%lp      = params%lpstart
            lpinfo(1)%l_lpset = .true.
            do istage = 2, NSTAGES-2
                lpinfo(istage)%lp      = lpinfo(istage-1)%lp - (lpinfo(istage-1)%lp - params%lpstop)/2.0
                lpinfo(istage)%l_lpset = .true.
            end do
            lpinfo(NSTAGES-1:)%l_lpset = .false.
            lpinfo(NSTAGES-1:)%lp      = lpinfo(NSTAGES-2)%lp
            do istage = 1,NSTAGES
                print *, 'stage lpset lp: ', istage,lpinfo(istage)%l_lpset,lpinfo(istage)%lp
            end do
        end subroutine prep_dims_lplims

        subroutine prep_command_lines( cline )
            use simple_default_clines, only: set_automask2D_defaults
            class(cmdline),   intent(in) :: cline
            cline_cluster2D        = cline
            cline_calc_pspec_distr = cline
            ! initial sigma2
            call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )
            ! cluster2D
            call cline_cluster2D%set('prg',       'cluster2D_distr')
            call cline_cluster2D%set('wiener',    'full')
            call cline_cluster2D%set('objfun',    'euclid')
            call cline_cluster2D%set('cc_iters',  0) ! euclid from the top
            call cline_cluster2D%set('sigma_est', 'group')
            call cline_cluster2D%set('ptclw',     'no')
            call cline_cluster2D%set('ml_reg',    'no')
            call cline_cluster2D%set('kweight',   'default')
            call cline_cluster2D%set('cenlp',     params%cenlp)
            call set_automask2D_defaults( cline_cluster2D )
            ! Single downscaling stage
            call cline_cluster2D%set('box_crop',  params%box_crop)
            call cline_cluster2D%set('smpd_crop', params%smpd_crop)
        end subroutine prep_command_lines

        ! generates starting references
        subroutine prep_refs
            use simple_commander_imgproc,   only: scale_commander
            use simple_commander_cluster2D, only: make_cavgs_commander_distr
            use simple_procimgstk,          only: random_selection_from_imgfile
            type(make_cavgs_commander_distr) :: xmake_cavgs
            type(scale_commander)            :: xscale
            type(cmdline)                    :: cline_make_cavgs, cline_scalerefs
            type(image)                      :: noiseimg
            character(len=:), allocatable    :: refs_sc
            integer :: ldim(3), iptcl, n, icls
            logical :: l_scale
            l_scale = .false.
            if( cline%defined('refs') )then
                ! existing references
                call find_ldim_nptcls(params%refs, ldim, n)
                if( n /= params%ncls ) THROW_HARD('Inconsistent number of classes in '//trim(params%refs))
                l_scale = ldim(1) == params%box_crop
            else
                params%refs = 'start2Drefs.mrc'
                select case(trim(params%cls_init))
                case('ptcl')
                    ! random selection of raw particles
                    call random_selection_from_imgfile(spproj, params%refs, params%box, params%ncls)
                    l_scale = params%l_autoscale
                case('randcls')
                    ! particles summed in random orientations
                    do iptcl=1,params%nptcls
                        if( spproj_field%get_state(iptcl) == 0 ) cycle
                        call spproj_field%set(iptcl, 'class', real(irnd_uni(params%ncls)))
                        call spproj_field%e3set(iptcl,ran3()*360.0)
                    end do
                    call spproj%write_segment_inside(params%oritype, params%projfile)
                    cline_make_cavgs = cline
                    call cline_make_cavgs%set('refs',      params%refs)
                    call cline_make_cavgs%set('box_crop',  params%box_crop)
                    call cline_make_cavgs%set('smpd_crop', params%smpd_crop)
                    call xmake_cavgs%execute_safe(cline_make_cavgs)
                    call cline_make_cavgs%kill
                case('rand')
                    ! from noise
                    call noiseimg%new([params%box_crop,params%box_crop,1], params%smpd_crop)
                    do icls = 1,params%ncls
                        call noiseimg%ran
                        call noiseimg%write(params%refs,icls)
                    enddo
                    call noiseimg%kill
                case DEFAULT
                    THROW_HARD('Unsupported mode of starting refernces generation: '//trim(params%cls_init))
                end select
            endif
            if( l_scale )then
                refs_sc = 'tmp'//params%ext
                call cline_scalerefs%set('stk',    params%refs)
                call cline_scalerefs%set('outstk', refs_sc)
                call cline_scalerefs%set('smpd',   params%smpd)
                call cline_scalerefs%set('newbox', params%box_crop)
                call xscale%execute_safe(cline_scalerefs)
                call cline_scalerefs%kill
                params%refs = 'start2Drefs.mrc'
                call simple_rename(refs_sc, params%refs)
                call del_file(refs_sc)
            endif
        end subroutine prep_refs

        subroutine set_cline_cluster2D( istage )
            integer, intent(in) :: istage
            character(len=:), allocatable :: sh_first, refine, center
            integer :: iphase, iter, imaxits
            real    :: trs, snr_noise_reg    
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
            refine = 'snhc_smpl'
            select case(iphase)
            case(1)
                snr_noise_reg = real(istage)
                select case(istage)
                case(1)
                    trs           = 0.
                    sh_first      = 'no'
                    imaxits       = 3
                    center        = 'no'
                    call cline_cluster2D%set('refs', params%refs)
                case(2,3,4,5)
                    imaxits       = 3 + (istage-1)*params%extr_lim/(NSTAGES-1)
                    trs           = lpinfo(istage)%trslim
                    sh_first      = 'yes'
                    center        = trim(params%center)
                end select
            case(2)
                imaxits       = MAXITS(2)
                sh_first      = 'yes'
                snr_noise_reg = 6.0
                center        = trim(params%center)
                ! deactivates stochastic withdrawal
                call cline_cluster2D%set('extr_iter', params%extr_lim+1)
            end select
            ! command line update
            call cline_cluster2D%set('startit',       iter)
            call cline_cluster2D%set('refine',        refine)
            if( lpinfo(istage)%l_lpset )then
                ! call cline_cluster2D%set('lp',        lpinfo(istage)%lp)
                call cline_cluster2D%set('lp',        lpinfo(1)%lp)
            else
                call cline_cluster2D%delete('lp')
            endif
            call cline_cluster2D%set('maxits',        imaxits)
            call cline_cluster2D%set('trs',           trs)
            call cline_cluster2D%set('sh_first',      sh_first)
            call cline_cluster2D%set('snr_noise_reg', snr_noise_reg)
            call cline_cluster2D%set('center',        center)
        end subroutine set_cline_cluster2D

    end subroutine exec_abinitio2D

end module simple_commander_abinitio2D
