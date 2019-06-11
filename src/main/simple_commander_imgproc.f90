! concrete commander: general image processing routines
module simple_commander_imgproc
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
implicit none

public :: binarise_commander
public :: edge_detector_commander
public :: convert_commander
public :: ctfops_commander
public :: filter_commander
public :: normalize_commander
public :: pspec_stats_commander
public :: scale_commander
public :: stack_commander
public :: stackops_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: binarise_commander
  contains
    procedure :: execute      => exec_binarise
end type binarise_commander
type, extends(commander_base) :: edge_detector_commander
  contains
    procedure :: execute      => exec_edge_detector
end type edge_detector_commander
type, extends(commander_base) :: convert_commander
  contains
    procedure :: execute      => exec_convert
end type convert_commander
type, extends(commander_base) :: ctfops_commander
  contains
    procedure :: execute      => exec_ctfops
end type ctfops_commander
type, extends(commander_base) :: filter_commander
  contains
    procedure :: execute      => exec_filter
end type filter_commander
type, extends(commander_base) :: normalize_commander
  contains
    procedure :: execute      => exec_normalize
end type normalize_commander
type, extends(commander_base) :: pspec_stats_commander
  contains
    procedure :: execute      => exec_pspec_stats
end type pspec_stats_commander
type, extends(commander_base) :: scale_commander
  contains
    procedure :: execute      => exec_scale
end type scale_commander
type, extends(commander_base) :: stack_commander
  contains
    procedure :: execute      => exec_stack
end type stack_commander
type, extends(commander_base) :: stackops_commander
  contains
    procedure :: execute      => exec_stackops
end type stackops_commander

contains

    !> for binarisation of stacks and volumes
    subroutine exec_binarise( self, cline )
        use simple_image, only: image
        class(binarise_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer :: igrow, iptcl
        ! error check
        if( .not. cline%defined('stk') .and. .not. cline%defined('vol1') )then
            THROW_HARD('ERROR! stk or vol1 needs to be present; simple_binarise')
        endif
        if( cline%defined('stk') .and. cline%defined('vol1') )then
            THROW_HARD('either stk or vol1 key can be present, not both; simple_binarise')
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            THROW_HARD('either thres-based or npix-based binarisation; both keys cannot be present; simple_binarise')
        endif
        if( cline%defined('stk') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            do iptcl=1,params%nptcls
                call build%img%read(params%stk, iptcl)
                call doit(build%img)
                call build%img%write(params%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            call build%vol%read(params%vols(1))
            call doit(build%vol)
            call build%vol%write(params%outvol)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_BINARISE NORMAL STOP ****')

        contains

            subroutine doit( img_or_vol )
                class(image), intent(inout) :: img_or_vol
                real :: ave, sdev, maxv, minv
                if( cline%defined('thres') )then
                    call img_or_vol%bin(params%thres)
                else if( cline%defined('npix') )then
                    call img_or_vol%bin(params%npix)
                else if( cline%defined('ndev') )then
                    call img_or_vol%stats( ave, sdev, maxv, minv)
                    call img_or_vol%bin(ave + params%ndev * sdev)
                else
                    call img_or_vol%bin_kmeans
                endif
                write(logfhandle,'(a,1x,i9)') '# FOREGROUND PIXELS:', img_or_vol%nforeground()
                write(logfhandle,'(a,1x,i9)') '# BACKGROUND PIXELS:', img_or_vol%nbackground()
                if( cline%defined('grow') )then
                    do igrow=1,params%grow
                        call img_or_vol%grow_bin
                    end do
                endif
                if( cline%defined('edge') ) call img_or_vol%cos_edge(params%edge)
                if( cline%defined('neg')  ) call img_or_vol%bin_inv
            end subroutine doit

    end subroutine exec_binarise

    !> for edge detection of stacks
    subroutine exec_edge_detector( self, cline )
        use simple_image, only: image
        class(edge_detector_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: igrow, iptcl
        ! error check
        if( .not. cline%defined('stk') )then
            THROW_HARD('ERROR! stk needs to be present; simple_segmentation')
        endif
        if( .not. cline%defined('detector') )then
            THROW_HARD('ERROR! detector needs to be present; simple_segmentation')
        endif
        if( .not. cline%defined('automatic') )then
            THROW_HARD('ERROR! automatic needs to be present; simple_segmentation')
        endif
        if(.not. cline%defined('outstk') )then
            params%outstk = 'Outstk.mrc'
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            THROW_HARD('either thres-based or npix-based edge detection; both keys cannot be present; simple_segmentation')
        endif
        if( cline%defined('thres') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres in automatic mode; simple_segmentation')
        endif
        if( cline%defined('npix') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose npix in automatic mode; simple_segmentation')
        endif
        if( cline%defined('thres_low') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_low in automatic mode; simple_segmentation')
        endif
        if( cline%defined('thres_up') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_up in automatic mode; simple_segmentation')
        endif
        if( cline%defined('thres') .and. params%detector .eq. 'canny') then
            THROW_HARD('canny needs double thresholding; simple_segmentation')
        endif

        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        do iptcl=1,params%nptcls
            call build%img%new([params%box,params%box,1],params%smpd,wthreads=.false.)
            call build%img%read(params%stk, iptcl)
            call doit(build%img)
            call build%img%write(params%outstk, iptcl)
            call build%img%kill
        end do
        ! end gracefully
        call simple_end('**** simple_segmentation NORMAL STOP ****')

    contains

            subroutine doit( img )
                use simple_segmentation
                class(image), intent(inout) :: img
                real, allocatable :: grad(:,:,:)
                real    :: thresh(1), ave, sdev, maxv, minv, lp(1)
                integer :: ldim(3)
                thresh = 0. !initialise
                ldim = img%get_ldim()
                allocate(grad(ldim(1), ldim(2), ldim(3)), source = 0.)
                if(cline%defined('lp')) lp(1) = params%lp
                select case ( params%detector )
                case ('sobel')
                    if( cline%defined('thres') )then
                        thresh(1) = params%thres
                        call sobel(img,thresh)
                    else if( cline%defined('npix') )then
                        call automatic_thresh_sobel(img,real(params%npix)/(real(ldim(1)*ldim(2))))
                    elseif( params%automatic .eq. 'yes') then
                        call img%scale_pixels([0.,255.])
                        call img%calc_gradient(grad)
                        call img%stats( ave, sdev, maxv, minv )
                        thresh(1) = ave + 0.7*sdev
                        write(logfhandle,*) 'Selected threshold: ', thresh
                        call sobel(img,thresh)
                        deallocate(grad)
                    else
                        THROW_HARD('If not automatic threshold needed; simple_segmentation')
                    endif
                case('canny')
                    if(ldim(3) .ne. 1) THROW_HARD('Canny for vol is not implemented; simple_segmentation')
                    if( params%automatic .eq. 'no' ) then
                        if(.not. cline%defined('thres_low') .or. .not. cline%defined('thres_up') )then
                            THROW_HARD('both upper and lower threshold needed; simple_segmentation')
                        else
                            if( cline%defined('lp')) then
                                call canny(img,[params%thres_low, params%thres_up],lp(1))
                            else
                                call canny(img,[params%thres_low, params%thres_up])
                            endif
                        endif
                    elseif( params%automatic .eq. 'yes') then
                        if(cline%defined('thres_low') .or. cline%defined('thres_up')) then
                            THROW_HARD('cannot define thresholds in automatic mode; simple_segmentation')
                        else
                            if( cline%defined('lp')) then
                                call canny(img,lp = lp(1))
                            else
                                call canny(img)
                            endif
                        endif
                    else
                        call canny(img)
                    endif
                ! case ('bin')
                !   call img%stats( ave, sdev, maxv, minv )
                !   call img%bin(ave+.7*sdev)
                case DEFAULT
                    THROW_HARD('Unknown detector argument; simple_segmentation')
               end select
            end subroutine

    end subroutine exec_edge_detector

    !> for converting between SPIDER and MRC formats
    subroutine exec_convert( self, cline )
        class(convert_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: iptcl
        if( cline%defined('stk') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            do iptcl=1,params%nptcls
                call progress(iptcl, params%nptcls)
                call build%img%read(params%stk, iptcl)
                call build%img%write(params%outstk, iptcl)
            end do
        else if( cline%defined('vol1') )then
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            call build%vol%read(params%vols(1))
            call build%img%write(params%outvol)
        else
            THROW_HARD('either vol1 or stk argument required to execute simple_convert')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CONVERT NORMAL STOP ****')
    end subroutine exec_convert

    !> for applying CTF to stacked images
    subroutine exec_ctfops( self, cline )
        use simple_procimgfile, only: apply_ctf_imgfile
        class(ctfops_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( cline%defined('oritab') .or. cline%defined('deftab') )then
        else
            THROW_HARD('oritab/deftab with CTF info needed for phase flipping/multiplication/CTF image generation')
        endif
        if( params%ctf .ne. 'no' )then
            select case( params%ctf )
                case( 'flip' )
                    if( params%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'flipneg')
                    else
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'flip')
                    endif
                case( 'yes' )
                    if( params%neg .eq. 'yes' )then
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'neg')
                    else
                        call apply_ctf_imgfile(params%stk, params%outstk, build%spproj_field, params%smpd, 'ctf')
                    endif
                case DEFAULT
                    THROW_HARD('Unknown ctf argument')
            end select
        else
            THROW_HARD('Nothing to do!')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CTFOPS NORMAL STOP ****')
    end subroutine exec_ctfops

    subroutine exec_filter( self, cline )
        use simple_procimgfile
        use simple_estimate_ssnr, only: fsc2optlp
        class(filter_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        real, allocatable :: fsc(:), optlp(:), res(:)
        real              :: width, fsc05, fsc0143
        integer           :: find
        width = 10.
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            if( cline%defined('width') ) width = params%width
            if( .not.file_exists(params%stk) ) THROW_HARD('cannot find input stack (stk)')
            if( cline%defined('filter') )then
                select case(trim(params%filter))
                case('tv')
                    call tvfilter_imgfile(params%stk, params%outstk, params%smpd, params%lam)
                case('nlmean')
                    call nlmean_imgfile(params%stk, params%outstk, params%smpd)
                case DEFAULT
                    THROW_HARD('Unknown filter!')
                end select
            else if( params%phrand .eq. 'no')then
                ! projection_frcs filtering
                if( cline%defined('frcs') )then
                    call matchfilt_imgfile(params%stk, params%outstk, params%frcs, params%smpd)
                ! Band pass
                else if( cline%defined('lp') .and. cline%defined('hp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, params%hp, params%lp, width=width)
                else if( cline%defined('lp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, 0., params%lp, width=width)
                else if( cline%defined('hp') )then
                    call bp_imgfile(params%stk, params%outstk, params%smpd, params%hp, 0., width=width)
                ! B-factor
            else if( cline%defined('bfac') )then
                    call apply_bfac_imgfile(params%stk, params%outstk, params%smpd, params%bfac)
                ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') .and. trim(params%real_filter) .ne. 'NLmean') THROW_HARD('need winsz input for real-space filtering')
                    call real_filter_imgfile(params%stk, params%outstk, params%smpd, trim(params%real_filter), nint(params%winsz))
                else
                    THROW_HARD('Nothing to do!')
                endif
            else if ( params%phrand.eq.'yes' )then
                ! Phase randomization
                if( .not. cline%defined('lp') ) THROW_HARD('low-pass limit needed 4 phase randomization')
                call phase_rand_imgfile(params%stk, params%outstk, params%smpd, params%lp)
            endif
        else
            ! 3D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
            if( cline%defined('width') ) width = params%width
            if( .not.file_exists(params%vols(1)) ) THROW_HARD('Cannot find input volume (vol1)')
            call build%vol%read(params%vols(1))
            if( params%phrand.eq.'no')then
                if( cline%defined('bfac') )then
                ! bfactor
                    call build%vol%apply_bfac(params%bfac)
                ! Band pass
                else if( cline%defined('hp') .and. cline%defined('lp') )then
                    call build%vol%bp(params%hp, params%lp, width=width)
                else if( cline%defined('hp') )then
                    call build%vol%bp(params%hp, 0., width=width)
                else if( params%tophat .eq. 'yes' .and. cline%defined('find') )then
                    call build%vol%tophat(params%find)
                else if( cline%defined('lp') )then
                    call build%vol%bp(0., params%lp, width=width)
                ! real-space
                else if( cline%defined('real_filter') )then
                    if( .not. cline%defined('winsz') ) THROW_HARD('need winsz input for real-space filtering')
                    call build%vol%real_space_filter(nint(params%winsz), params%real_filter)
                else if( cline%defined('vol_filt') )then
                    if( .not.file_exists(params%vol_filt)) THROW_HARD('cannot find volume filter (vol_filt)')
                    call build%vol2%read(params%vol_filt)
                    call build%vol%fft
                    call build%vol%apply_filter(build%vol2)
                    call build%vol2%kill
                    call build%vol%ifft
                else if( cline%defined('fsc') )then
                    if( .not.file_exists(params%fsc)) THROW_HARD('Cannot find FSC filter (vol_filt)')
                    ! resolution & optimal low-pass filter from FSC
                    fsc   = file2rarr(params%fsc)
                    optlp = fsc2optlp(fsc)
                    res   = build%vol%get_res()
                    call get_resolution( fsc, res, fsc05, fsc0143 )
                    where(res < TINY) optlp = 0.
                    call build%vol%apply_filter(optlp)
                else
                    THROW_HARD('Nothing to do!')
                endif
            else
                if( .not. cline%defined('lp') )THROW_HARD('low-pass limit needed 4 phase randomization')
                call build%vol%phase_rand(params%lp)
            endif
            if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FILTER NORMAL STOP ****')
    end subroutine exec_filter

    !> normalize is a program for normalization of MRC or SPIDER stacks and volumes.
    !! If you want to normalize your images inputted with stk, set norm=yes.
    !! hfun (e.g. hfun=sigm) controls the normalization function. If you want to
    !! perform noise normalization of the images set noise_norm=yes given a mask
    !! radius msk (pixels). If you want to normalize your images or volume
    !! (vol1) with respect to their power spectrum set shell_norm=yes
    subroutine exec_normalize( self, cline )
        use simple_procimgfile, only: norm_imgfile, noise_norm_imgfile, shellnorm_imgfile
        class(normalize_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        real, allocatable :: spec(:)
        if( cline%defined('stk')  .and. cline%defined('vol1') )THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            if( params%norm.eq.'yes' )then
                ! Normalization
                call norm_imgfile(params%stk, params%outstk, params%smpd)
            else if( params%noise_norm.eq.'yes' )then
                ! Noise normalization
                if( cline%defined('msk') )then
                    call noise_norm_imgfile(params%stk, params%msk, params%outstk, params%smpd)
                else
                    THROW_HARD('need msk parameter for noise normalization')
                endif
            else if( params%shellnorm.eq.'yes' )then
                ! shell normalization
                call shellnorm_imgfile( params%stk, params%outstk, params%smpd)
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true., boxmatch_off=.true.)
            if( .not.file_exists(params%vols(1)) )THROW_HARD('Cannot find input volume')
            call build%vol%read(params%vols(1))
            if( params%norm.eq.'yes' )then
                call build%vol%norm()
                call build%vol%write(params%outvol, del_if_exists=.true.)
            else if( params%shellnorm.eq.'yes' )then
                ! shell normalization
                call build%vol%shellnorm()
                call build%vol%spectrum('power', spec)
                call build%vol%write(params%outvol, del_if_exists=.true.)
            else
                THROW_HARD('Normalization type not implemented yet')
            endif
        else
            THROW_HARD('No input images(s) or volume provided')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_NORMALIZE NORMAL STOP ****')
    end subroutine exec_normalize

    !> for online power spectra analysis through statistic calculation
    !! to perform prior ctf estimation.
    subroutine exec_pspec_stats( self, cline )
        use simple_genpspec_and_statistics, only: pspec_statistics
        class(pspec_stats_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)       :: params
        type(pspec_statistics) :: pspec
        integer :: ldim(3), nptcls
        real    :: smpd
        call params%new(cline)
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_pspec_stats')
        endif
        if( .not. cline%defined('fname') )then
            THROW_HARD('ERROR! fname needs to be present; exec_pspec_stats')
        endif
        call pspec%new(params%fname,8) !it also create a new instance of the class powerspectrum
        call find_ldim_nptcls (params%fname,ldim, nptcls, smpd)
        call pspec%run()
        ! end gracefully
        call simple_end('**** SIMPLE_PSPEC_STATS NORMAL STOP ****')
    end subroutine exec_pspec_stats

    !> provides re-scaling and clipping routines for MRC or SPIDER stacks and volumes
    recursive subroutine exec_scale( self, cline )
        use simple_procimgfile, only: resize_imgfile_double, resize_and_clip_imgfile, resize_imgfile, pad_imgfile, clip_imgfile
        use simple_image,       only: image
        use simple_qsys_funs,   only: qsys_job_finished
        class(scale_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: vol2, img, img2
        real             :: ave, sdev, var, med, smpd_new, scale
        integer          :: ldim(3), ldim_scaled(3), nfiles, nframes, iframe, ifile
        !integer          :: ldims_scaled(2,3)
        character(len=:), allocatable :: fname
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        if( cline%defined('stk') .and. cline%defined('vol1') ) THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            call img%new([params%box,params%box,1],params%smpd)
            call img2%new([params%box,params%box,1],params%smpd)
            if( cline%defined('scale') .or. cline%defined('newbox') )then
                call del_file(params%outstk)
                ! Rescaling
                ldim_scaled = [params%newbox,params%newbox,1] ! dimension of scaled
                if( cline%defined('clip') )then
                    call resize_and_clip_imgfile(params%stk,params%outstk,params%smpd,ldim_scaled,&
                    [params%clip,params%clip,1],smpd_new)
                else
                    call resize_imgfile(params%stk,params%outstk,params%smpd,ldim_scaled,smpd_new)
                    write(logfhandle,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', ldim_scaled(1)
                endif
                write(logfhandle,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
            else if( cline%defined('clip') )then
                call del_file(params%outstk)
                ! Clipping
                if( params%clip > params%box )then
                    call pad_imgfile(params%stk,params%outstk,[params%clip,params%clip,1],params%smpd)
                else
                    call clip_imgfile(params%stk,params%outstk,[params%clip,params%clip,1],params%smpd)
                endif
            endif
        else if( cline%defined('vol1') )then
            ! 3D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true., boxmatch_off=.true.)
            if( .not.file_exists(params%vols(1)) ) THROW_HARD('Cannot find input volume')
            call build%vol%read(params%vols(1))
            if( cline%defined('scale') .or. cline%defined('newbox') )then
                ! Rescaling
                call vol2%new([params%newbox,params%newbox,params%newbox],params%smpd)
                call build%vol%fft()
                call vol2%set_ft(.true.)
                if( params%newbox < params%box )then
                    call build%vol%clip(vol2)
                else if( params%newbox > params%box )then
                    call build%vol%pad(vol2)
                else
                    call vol2%copy(build%vol)
                endif
                call build%vol%copy(vol2)
                call build%vol%ifft()
                scale = real(params%newbox)/real(params%box)
                params%box = params%newbox
                smpd_new = params%smpd/scale
                write(logfhandle,'(a,1x,f9.4)') 'SAMPLING DISTANCE AFTER SCALING:', smpd_new
            endif
            if( cline%defined('clip') )then
                ! Clipping
                call vol2%new([params%clip,params%clip,params%clip],params%smpd)
                if( params%clip < params%box )then
                    call build%vol%clip(vol2)
                else
                    call build%vol%stats( 'background', ave, sdev, var, med )
                    call build%vol%pad(vol2, backgr=med)
                endif
                call build%vol%copy(vol2)
            else
                 write(logfhandle,'(a,1x,i5)') 'BOX SIZE AFTER SCALING:', params%newbox
            endif
            if( params%outvol .ne. '' )call build%vol%write(params%outvol, del_if_exists=.true.)
        else if( cline%defined('filetab') )then
            call params%new(cline)
            call read_filetable(params%filetab, filenames)
            nfiles = size(filenames)
            call find_ldim_nptcls(filenames(1),ldim,nframes)
            ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
            if( cline%defined('newbox') )then
                ldim_scaled = [params%newbox,params%newbox,1] ! dimension of scaled
                params%scale = real(params%newbox) / real(ldim(1))
            else if( cline%defined('scale') )then
                ldim_scaled(1) = round2even(real(ldim(1))*params%scale)
                ldim_scaled(2) = round2even(real(ldim(2))*params%scale)
                ldim_scaled(3)   = 1
            else
                THROW_HARD('filetab key only in combination with scale or newbox!')
            endif
            call img%new(ldim,params%smpd)
            call img2%new(ldim_scaled,params%smpd/params%scale)
            do ifile=1,nfiles
                call progress(ifile, nfiles)
                if( cline%defined('dir_target') )then
                    fname = filepath(trim(params%dir_target),add2fbody(basename(filenames(ifile)), params%ext, SCALE_SUFFIX))
                else
                    fname = add2fbody(trim(filenames(ifile)), params%ext, '_sc')
                endif
                call find_ldim_nptcls(filenames(ifile),ldim,nframes)
                ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
                do iframe= 1, nframes
                    call img%read(filenames(ifile), iframe)
                    call img%fft()
                    if( ldim_scaled(1) <= ldim(1) .and. ldim_scaled(2) <= ldim(2) .and. ldim_scaled(3) <= ldim(3) )then
                        call img%clip(img2)
                    else
                        call img%pad(img2)
                    endif
                    call img2%ifft()
                    call img2%write(fname, iframe)
                end do
                deallocate(fname)
            end do
        else
            THROW_HARD('SIMPLE_SCALE needs input image(s) or volume or filetable!')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****', print_simple=.false.)
        call qsys_job_finished(  'simple_commander_imgproc :: exec_scale' )
    end subroutine exec_scale

    !>  for stacking individual images or multiple stacks into one
    subroutine exec_stack( self, cline )
        use simple_image,  only: image
        class(stack_commander), intent(inout)  :: self
        class(cmdline),         intent(inout)  :: cline
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        type(parameters) :: params
        type(builder)    :: build
        integer          :: nfiles, ldim(3), ifile, ifoo, cnt
        integer          :: lfoo(3), nimgs, iimg
        type(image)      :: tmp
        real             :: mm(2)
        if( cline%defined('lp') )then
            if( .not. cline%defined('smpd') ) THROW_HARD('smpd (sampling distance) needs to be defined if lp is')
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        call find_ldim_nptcls(filenames(1),ldim,ifoo)
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        params%box = ldim(1)
        ! prepare build%img and tmp for reading
        call build%img%new([params%box,params%box,1], params%smpd)
        if( cline%defined('clip') ) call tmp%new([params%clip,params%clip,1], params%smpd)
        ! loop over files
        cnt = 0
        do ifile=1,nfiles
            if( .not. file_exists(filenames(ifile)) )then
                write(logfhandle,*) 'inputted spec file does not exist: ', trim(adjustl(filenames(ifile)))
            endif
            call find_ldim_nptcls(filenames(ifile),lfoo,nimgs)
            do iimg=1,nimgs
                cnt = cnt+1
                call build%img%read(filenames(ifile), iimg, readhead=.false.)
                if( cline%defined('clip') )then
                    call build%img%clip(tmp)
                    mm = tmp%minmax()
                    DebugPrint 'min/max: ', mm(1), mm(2)
                    call tmp%write(params%outstk, cnt)
                else
                    call build%img%write(params%outstk, cnt)
                endif
            end do
            call progress(ifile, nfiles)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_STACK NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_stack

    !> provides standard single-particle image
    !> processing routines that are applied to MRC or SPIDER stacks.
    subroutine exec_stackops( self, cline )
        use simple_oris, only: oris
        use simple_ori,  only: ori
        use simple_stackops
        use simple_procimgfile
        class(stackops_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(ran_tabu)                :: rt
        type(oris)                    :: o_here, os_ran
        type(ori)                     :: o, o2
        integer,          allocatable :: pinds(:)
        character(len=:), allocatable :: fname
        integer :: i, s, cnt, nincl
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.,boxmatch_off=.true.)
        ! random selection
        if( cline%defined('nran') )then
            write(logfhandle,'(a)') '>>> RANDOMLY SELECTING IMAGES'
            allocate( pinds(params%nran), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('In: simple_commander; stackops',alloc_stat)
            rt = ran_tabu(params%nptcls)
            call rt%ne_ran_iarr(pinds)
            if( cline%defined('oritab') .or. cline%defined('deftab') )then
                call del_file(params%outfile)
            endif
            call os_ran%new(params%nran)
            do i=1,params%nran
                call progress(i, params%nran)
                call build%img%read(params%stk, pinds(i))
                call build%img%write(params%outstk, i)
                if( cline%defined('oritab') .or. cline%defined('deftab') )then
                    call build%spproj_field%get_ori(pinds(i), o)
                    call os_ran%set_ori(i, o)
                    call o%kill
                endif
            end do
            if( cline%defined('oritab') .or. cline%defined('deftab') )then
                call build%spproj_field%write(params%outfile, [1,params%nran])
            endif
            goto 999
        endif
        ! fishing expeditions
        ! frac only
        if( cline%defined('frac') )then
            if( params%oritab == '' ) THROW_HARD('need input orientation doc for fishing expedition; simple_stackops')
            ! determine how many particles to include
            if( params%frac < 0.99 )then
                nincl = nint(real(params%nptcls)*params%frac)
            else
                nincl = params%nptcls
            endif
            ! order the particles
            pinds = build%spproj_field%order()
            ! fish the best ones out
            if( cline%defined('state') )then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'state'))
                    if( s == params%state )then
                        cnt = cnt+1
                        call build%img%read(params%stk, pinds(i))
                        call build%img%write(params%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'state'))
                    if( s == params%state )then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(pinds(i), o2)
                        call o_here%set_ori(cnt, o2)
                    endif
                end do
                allocate(fname, source='extracted_oris_state'//int2str_pad(params%state,2)//trim(TXT_EXT))
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'class'))
                    if( s == params%class )then
                        cnt = cnt+1
                        call build%img%read(params%stk, pinds(i))
                        call build%img%write(params%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,nincl
                    call progress(i, nincl)
                    s = nint(build%spproj_field%get(pinds(i), 'class'))
                    if( s == params%class )then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(pinds(i), o2)
                        call o_here%set_ori(cnt, o2)
                    endif
                end do
                allocate(fname, source='extracted_oris_class'//int2str_pad(params%class,5)//trim(TXT_EXT))
            else
                o_here = oris(nincl)
                do i=1,nincl
                    call progress(i, nincl)
                    call build%img%read(params%stk, pinds(i))
                    call build%img%write(params%outstk, i)
                    call build%spproj_field%get_ori(pinds(i), o2)
                    call o_here%set_ori(i, o2)
                end do
                allocate(fname, source='extracted_oris'//trim(TXT_EXT))
            endif
            if( cline%defined('outfile') )then
                call o_here%write(params%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            endif
            goto 999
        endif
        ! state/class + frac
        if( (cline%defined('state').or.cline%defined('class')) .and. .not.cline%defined('frac') )then
            if( params%oritab == '' ) THROW_HARD('need input orientation doc for fishing expedition; simple_stackops')
            if( cline%defined('state') )then
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'state'))
                    if( s == params%state )then
                        cnt = cnt+1
                        call build%img%read(params%stk, i)
                        call build%img%write(params%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the extracted ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'state'))
                    if( s == params%state )then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(i, o2)
                        call o_here%set_ori(cnt, o2)
                    endif
                end do
                allocate(fname, source='extracted_oris_state'//int2str_pad(params%state,2)//trim(TXT_EXT))
            else if( cline%defined('class') )then
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'class'))
                    if( s == params%class )then
                        cnt = cnt+1
                        call build%img%read(params%stk, i)
                        call build%img%write(params%outstk, cnt)
                    endif
                end do
                ! make orientation structure for the best ones
                o_here = oris(cnt)
                cnt = 0
                do i=1,params%nptcls
                    call progress(i, params%nptcls)
                    s = nint(build%spproj_field%get(i, 'class'))
                    if( s == params%class )then
                        cnt = cnt+1
                        call build%spproj_field%get_ori(i, o2)
                        call o_here%set_ori(cnt, o2)
                    endif
                end do
                allocate(fname, source='extracted_oris_class'//int2str_pad(params%class,5)//trim(TXT_EXT))
            endif
            if( cline%defined('outfile') )then
                call o_here%write(params%outfile, [1,o_here%get_noris()])
            else
                call o_here%write(fname, [1,o_here%get_noris()])
            endif
            goto 999
        endif
        ! invert contrast
        if( params%neg .eq. 'yes' )then
            call neg_imgfile(params%stk, params%outstk, params%smpd)
            goto 999
        endif
        ! auto correlation function
        if( params%acf .eq. 'yes' )then
            call acf_stack(params%stk, params%outstk)
            goto 999
        endif
        ! produce image statistics
        if( params%stats .eq. 'yes' )then
            call stats_imgfile(params%stk, build%spproj_field)
            call build%spproj_field%write('image_statistics.txt')
            goto 999
        endif
        ! create frame averages
        if( params%nframesgrp > 0 )then
            call frameavg_stack(params%stk, params%outstk, params%nframesgrp, params%smpd)
            goto 999
        endif
        ! visualize
        if( params%vis .eq. 'yes' )then
            do i=1,params%nptcls
                call build%img%read(params%stk, i)
                call build%img%vis()
            end do
            goto 999
        endif
        ! average
        if( params%avg .eq. 'yes' )then
            call make_avg_stack(params%stk, params%outstk, params%smpd)
            goto 999
        endif
        ! add noise
        if( cline%defined('snr') )then
            call add_noise_imgfile(params%stk, params%outstk, params%snr, params%smpd)
            goto 999
        endif
        ! copy
        if( cline%defined('top') .and. .not. cline%defined('part') )then
            call copy_imgfile(params%stk, params%outstk, params%smpd, fromto=[params%fromp,params%top])
            goto 999
        endif
        ! mirror
        if( params%mirr .ne. 'no' )then
            call mirror_imgfile(params%stk, params%outstk, params%mirr, params%smpd)
            goto 999
        endif
        if( params%subtr_backgr .eq. 'yes' )then
            if( .not. cline%defined('lp') ) THROW_HARD('need lp arg for background subtraction; exec_stackops')
            call subtr_backgr_imgfile(params%stk, params%outstk, params%smpd, params%lp)
            goto 999
        endif
        ! default
        write(logfhandle,*)'Nothing to do!'
        ! end gracefully
999     call o_here%kill
        call os_ran%kill
        call o%kill
        call o2%kill
        call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
    end subroutine exec_stackops

end module simple_commander_imgproc
