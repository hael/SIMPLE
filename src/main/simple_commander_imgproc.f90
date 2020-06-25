! concrete commander: general image processing routines
module simple_commander_imgproc
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_binimage,       only: binimage
implicit none

public :: binarise_commander
public :: edge_detect_commander
public :: convert_commander
public :: ctfops_commander
public :: filter_commander
public :: normalize_commander
public :: pspec_stats_commander
public :: scale_commander
public :: stack_commander
public :: stackops_commander
public :: pspec_int_rank_commander
public :: estimate_diam_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: binarise_commander
  contains
    procedure :: execute      => exec_binarise
end type binarise_commander
type, extends(commander_base) :: edge_detect_commander
  contains
    procedure :: execute      => exec_edge_detect
end type edge_detect_commander
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
type, extends(commander_base) :: pspec_int_rank_commander
  contains
    procedure :: execute      => exec_pspec_int_rank
end type pspec_int_rank_commander
type, extends(commander_base) :: estimate_diam_commander
  contains
    procedure :: execute      => exec_estimate_diam
end type estimate_diam_commander

contains

    !> for binarisation of stacks and volumes
    subroutine exec_binarise( self, cline )
        use simple_segmentation, only: otsu_robust_fast
        class(binarise_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters) :: params
        type(binimage)   :: img_or_vol
        integer :: igrow, iptcl
        logical :: is2D, otsu, fill_holes
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
        call params%new(cline)
        if( cline%defined('fill_holes') )then
          if(params%fill_holes .eq. 'yes') then
            fill_holes = .true.
          else
            fill_holes = .false.
          endif
        endif
        if( cline%defined('stk') )then
            is2D = .true.
            call img_or_vol%new_bimg([params%box,params%box,1], params%smpd)
            do iptcl=1,params%nptcls
                call img_or_vol%read(params%stk, iptcl)
                call doit(otsu)
                if(fill_holes) then
                  call img_or_vol%fill_holes
                endif
                call img_or_vol%write(params%outstk, iptcl)
            end do
            if(otsu) write(logfhandle,*) 'Method applied for binarisation: OTSU algorithm'
        else if( cline%defined('vol1') )then
            is2D = .false.
            call img_or_vol%new_bimg([params%box,params%box,params%box], params%smpd)
            call img_or_vol%read(params%vols(1))
            call doit(otsu)
            call img_or_vol%write(params%outvol)
            if(otsu) write(logfhandle,*) 'Method applied for binarisation: OTSU algorithm'
        endif
        call img_or_vol%kill
        ! end gracefully
        call simple_end('**** SIMPLE_BINARISE NORMAL STOP ****')

        contains

            subroutine doit(otsu)
                type(binimage) :: cos_img
                real :: ave, sdev, maxv, minv, thresh(3)
                logical, optional, intent(inout) :: otsu
                otsu =  .false.
                if( cline%defined('thres') )then
                    call img_or_vol%binarize(params%thres)
                else if( cline%defined('npix') )then
                    call img_or_vol%binarize(params%npix)
                else if( cline%defined('ndev') )then
                    call img_or_vol%stats( ave, sdev, maxv, minv)
                    call img_or_vol%binarize(ave + params%ndev * sdev)
                else
                    call otsu_robust_fast(img_or_vol, is2D, noneg=.false., thresh=thresh)
                    otsu = .true.
                endif
                write(logfhandle,'(a,1x,i9)') '# FOREGROUND PIXELS:', img_or_vol%nforeground()
                write(logfhandle,'(a,1x,i9)') '# BACKGROUND PIXELS:', img_or_vol%nbackground()
                if( cline%defined('grow') ) call img_or_vol%grow_bins(params%grow)
                if( cline%defined('edge') ) then
                    call img_or_vol%cos_edge(params%edge,cos_img)
                    call img_or_vol%copy_bimg(cos_img)
                    call cos_img%kill
                endif
                if( cline%defined('neg')  ) call img_or_vol%bin_inv
            end subroutine doit

    end subroutine exec_binarise

    !> for edge detection of stacks
    subroutine exec_edge_detect( self, cline )
        class(edge_detect_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: igrow, iptcl
        ! error check
        if( .not. cline%defined('stk') )then
            THROW_HARD('ERROR! stk needs to be present; exec_edge_detect')
        endif
        if( .not. cline%defined('detector') )then
            THROW_HARD('ERROR! detector needs to be present; exec_edge_detect')
        endif
        if( .not. cline%defined('automatic') )then
            THROW_HARD('ERROR! automatic needs to be present; exec_edge_detect')
        endif
        if(.not. cline%defined('outstk') )then
            params%outstk = 'Outstk.mrc'
        endif
        if( cline%defined('thres') .and. cline%defined('npix') )then
            THROW_HARD('either thres-based or npix-based edge detection; both keys cannot be present; exec_edge_detect')
        endif
        if( cline%defined('thres') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('npix') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose npix in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres_low') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_low in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres_up') .and. params%automatic .eq. 'yes') then
            THROW_HARD('cannot chose thres_up in automatic mode; exec_edge_detect')
        endif
        if( cline%defined('thres') .and. params%detector .eq. 'canny') then
            THROW_HARD('canny needs double thresholding; exec_edge_detect')
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
        call simple_end('**** SIMPLE_EDGE_DETECT NORMAL STOP ****')

    contains

            subroutine doit( img )
                use simple_segmentation
                class(image), intent(inout) :: img
                type (image)      :: img_grad
                real, allocatable :: grad(:,:,:)
                real    :: thresh(1), ave, sdev, maxv, minv, lp(1)
                real    :: smpd
                integer :: ldim(3)
                thresh = 0. !initialise
                ldim = img%get_ldim()
                smpd = img%get_smpd()
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
                        call img_grad%new(ldim, smpd)
                        call img_grad%set_rmat(grad)
                        call img_grad%stats( ave, sdev, maxv, minv )
                        call img_grad%kill
                        thresh(1) = ave + 0.7*sdev
                        write(logfhandle,*) 'Selected threshold: ', thresh
                        call sobel(img,thresh)
                        deallocate(grad)
                    else
                        THROW_HARD('If not automatic threshold needed; exec_edge_detect')
                    endif
                case('canny')
                    if(ldim(3) .ne. 1) THROW_HARD('Canny for vol is not implemented; exec_edge_detect')
                    if( params%automatic .eq. 'no' ) then
                        if(.not. cline%defined('thres_low') .or. .not. cline%defined('thres_up') )then
                            THROW_HARD('both upper and lower threshold needed; exec_edge_detect')
                        else
                            if( cline%defined('lp')) then
                                call canny(img,img,thresh=[params%thres_low, params%thres_up],lp=lp(1)) !inout/output image coincide
                            else
                                call canny(img,img,thresh=[params%thres_low, params%thres_up])
                            endif
                        endif
                    elseif( params%automatic .eq. 'yes') then
                        if(cline%defined('thres_low') .or. cline%defined('thres_up')) then
                            THROW_HARD('cannot define thresholds in automatic mode; exec_edge_detect')
                        else
                            if( .not. cline%defined('lp')) then
                              THROW_HARD('Canny automatic requires lp in input; exec_edge_detect')
                            else
                                call canny(img,lp = lp(1))
                            endif
                        endif
                    else
                      THROW_HARD('Wrong input for automatic parameter!; exec_edge_detect')
                    endif
                case DEFAULT
                    THROW_HARD('Unknown detector argument; exec_edge_detect')
               end select
            end subroutine

    end subroutine exec_edge_detect

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
        ! cleanup
        call build%kill_general_tbox
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
        if( cline%defined('oritab') .or. cline%defined('deftab') )then
        else
            THROW_HARD('oritab/deftab with CTF info needed for phase flipping/multiplication/CTF image generation')
        endif
        if( .not. cline%defined('stk') ) call cline%set('box', 256.)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
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
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_CTFOPS NORMAL STOP ****')
    end subroutine exec_ctfops

    subroutine exec_filter( self, cline )
        use simple_procimgfile
        use simple_estimate_ssnr, only: fsc2optlp
        use simple_tvfilter,      only: tvfilter
        class(filter_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(tvfilter)    :: tvfilt
        real, allocatable :: fsc(:), optlp(:), res(:)
        real, parameter   :: SIGMA_DEFAULT=1.0
        real              :: width, fsc05, fsc0143, sigma
        integer           :: find
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'filtered.mrcs')
        if( .not. cline%defined('outvol') ) call cline%set('outvol', 'filtered.mrc')
        width = 10.
        if( cline%defined('stk') )then
            ! 2D
            call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
            if( cline%defined('width') ) width = params%width
            if( .not.file_exists(params%stk) ) THROW_HARD('cannot find input stack (stk)')
            if( cline%defined('filter') )then
                select case(trim(params%filter))
                    case('tv')
                        if( trim(params%tseries) .eq. 'yes' )then
                            call tvfilt_tseries_imgfile(params%stk, params%outstk, params%smpd, params%lambda)
                        else
                            call tvfilter_imgfile(params%stk, params%outstk, params%smpd, params%lambda)
                        endif
                    case('nlmean')
                        if( cline%defined('sigma') )then
                            call nlmean_imgfile(params%stk, params%outstk, params%smpd, params%sigma)
                        else
                            call nlmean_imgfile(params%stk, params%outstk, params%smpd)
                        endif
                    case('corr')
                        if( .not. cline%defined('lp') ) THROW_HARD('correlation filter requires lp (low-pass limit) input; exec_filter')
                        sigma = SIGMA_DEFAULT
                        if( cline%defined('sigma') ) sigma = params%sigma
                        if( trim(params%tseries) .eq. 'yes' )then
                            call corrfilt_tseries_imgfile(params%stk, params%outstk, params%smpd, sigma, params%lp)
                        else
                            call corrfilt_imgfile(params%stk, params%outstk, params%smpd, sigma, params%lp, isvol=.false.)
                        endif
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
                    if(cline%defined('filter') .and. (trim(params%filter) .eq. 'corr')) then
                        if(cline%defined('sigma')) then
                            call corrfilt_imgfile(params%vols(1), params%outvol, params%smpd, params%sigma, lp=params%lp, isvol=.true.)
                            call build%vol%read(params%outvol)
                        elseif(cline%defined('element')) then
                            call corrfilt_imgfile(params%vols(1), params%outvol, params%element(1:2), params%smpd, lp=params%lp)
                            call build%vol%read(params%outvol)
                        else
                            call corrfilt_imgfile(params%vols(1),params%outvol, params%smpd, SIGMA_DEFAULT, lp=params%lp, isvol=.true.)
                            call build%vol%read(params%outvol)
                        endif
                    else
                        call build%vol%bp(0., params%lp, width=width)
                    endif
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
                else if( cline%defined('filter') )then
                    select case(trim(params%filter))
                        case('tv')
                            call tvfilt%new
                            call tvfilt%apply_filter_3d(build%vol, params%lambda)
                            call tvfilt%kill
                        case('nlmean')
                            call build%vol%NLmean3D()
                        case('corr')
                            if(.not. cline%defined('lp')) then
                                THROW_HARD('Correlation filter needs lp as input!; exec_filter')
                            endif
                        case DEFAULT
                            THROW_HARD('Unknown filter!')
                    end select
                else
                    THROW_HARD('Nothing to do!')
                endif
            else
                if( .not. cline%defined('lp') )THROW_HARD('low-pass limit needed 4 phase randomization')
                call build%vol%phase_rand(params%lp)
            endif
            if( params%outvol .ne. '' ) call build%vol%write(params%outvol, del_if_exists=.true.)
        endif
        ! cleanup
        call build%kill_general_tbox
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
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
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
        ! cleanup
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_NORMALIZE NORMAL STOP ****')
    end subroutine exec_normalize

    !> for online power spectra analysis through statistic calculation
    !! to perform prior ctf estimation.
    subroutine exec_pspec_stats( self, cline )
        use simple_pspec_stats, only: pspec_stats
        class(pspec_stats_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline !< command line input
        type(pspec_stats) :: pspec
        type(parameters)  :: params
        integer :: nfiles,n_mic
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        real, allocatable    :: score(:)
        character(len = 100) :: iom
        integer :: status
        logical :: keep
        ! sanity checks
        if( .not. cline%defined('filetab') )then
            THROW_HARD('ERROR! filetab needs to be present; exec_pspec_stats')
        endif
        if( .not. cline%defined('smpd') )then
            THROW_HARD('ERROR! smpd needs to be present; exec_pspec_stats')
        endif
        call params%new(cline)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        allocate(score(nfiles), source = 0.)
        do n_mic = 1, nfiles
            call pspec%new(filenames(n_mic),params%smpd)
            call pspec%process_ps()
            call pspec%get_output(keep)
            write(logfhandle,*) trim(filenames(n_mic)), ' keep: ', keep , '  CTF_ccscore: ', pspec%get_score(), '  CTF_curvature: ', pspec%get_curvature(), 'CTF_value: ', pspec%get_value()
            call pspec%kill()
        enddo
        ! end gracefully
        call simple_end('**** SIMPLE_PSPEC_STATS NORMAL STOP ****')
    end subroutine exec_pspec_stats

    !> provides re-scaling and clipping routines for MRC or SPIDER stacks and volumes
    subroutine exec_scale( self, cline )
        use simple_procimgfile, only: resize_imgfile_double, resize_and_clip_imgfile, resize_imgfile, pad_imgfile, clip_imgfile
        use simple_qsys_funs, only: qsys_job_finished
        class(scale_commander), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(image)      :: vol2, img, img2
        integer, allocatable :: states(:)
        logical, allocatable :: ptcl_msk(:)
        real             :: ave, sdev, var, med, smpd_new, scale
        integer          :: ldim(3), ldim_scaled(3), nfiles, nframes, iframe, ifile
        integer          :: istk, nstks, ptcl_fromp, ptcl_top
        character(len=:), allocatable :: fname, stkin, stkout
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        character(len=STDLEN)         :: ext
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( cline%defined('stk') .and. cline%defined('vol1') ) THROW_HARD('Cannot operate on images AND volume at once')
        if( cline%defined('projfile') .and. cline%defined('fromp') .and. cline%defined('top')&
            &.and..not.cline%defined('stk') .and. .not.cline%defined('vol1') .and. .not.cline%defined('filetab'))then
            ! input comes from distributed execution scale_project
            call cline%set('oritype', 'ptcl2D')
            call build%init_params_and_build_spproj(cline, params)
            call build%spproj%read_segment('stk',params%projfile)
            nstks = build%spproj%os_stk%get_noris()
            if( nstks == 0 ) THROW_HARD('NO EXISTING STACK IN PROJECT')
            if( cline%defined('newbox') )then
                ldim_scaled  = [params%newbox,params%newbox,1]
                params%scale = real(params%newbox) / real(params%box)
            else if( cline%defined('scale') )then
                ldim_scaled(1) = round2even(real(params%box)*params%scale)
                ldim_scaled(2) = round2even(real(params%box)*params%scale)
                ldim_scaled(3) = 1
            else
                THROW_HARD('MISSING NEW DIMENSION!')
            endif
            ! state handling
            states = nint(build%spproj_field%get_all('state'))
            if( cline%defined('state') )then
                if( .not.any(states == params%state) )then
                    THROW_HARD('NO PARTICLE MATCHING STATE FOUND. STATE: '//int2str(params%state))
                endif
                where( states /= params%state ) states = 0
            else
                ! all inclusive by default
                params%state = 1
                states = 1
            endif
            do istk = params%fromp, params%top
                call progress(istk, nstks)
                stkin      = build%spproj%get_stkname(istk)
                ptcl_fromp = nint(build%spproj%os_stk%get(istk,'fromp'))
                ptcl_top   = nint(build%spproj%os_stk%get(istk,'top'))
                ptcl_msk   = (states(ptcl_fromp:ptcl_top) == 1)
                ext        = '.'//fname2ext(stkin)
                if( cline%defined('dir_target') )then
                    stkout = filepath(trim(params%dir_target),add2fbody(basename(stkin),ext,SCALE_SUFFIX))
                else
                    stkout = add2fbody(trim(stkin),ext,SCALE_SUFFIX)
                endif
                call resize_imgfile(stkin, stkout, params%smpd, ldim_scaled, smpd_new, mask=ptcl_msk)
            enddo
            deallocate(states,ptcl_msk)
        else
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
                        fname = add2fbody(trim(filenames(ifile)), params%ext, SCALE_SUFFIX)
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
        endif
        ! cleanup
        call vol2%kill
        call img%kill
        call img2%kill
        call build%kill_general_tbox
        ! end gracefully
        call simple_end('**** SIMPLE_SCALE NORMAL STOP ****', print_simple=.false.)
        call qsys_job_finished(  'simple_commander_imgproc :: exec_scale' )
    end subroutine exec_scale

    !>  for stacking individual images or multiple stacks into one
    subroutine exec_stack( self, cline )
        class(stack_commander), intent(inout)  :: self
        class(cmdline),         intent(inout)  :: cline
        character(len=LONGSTRLEN), allocatable :: filenames(:)
        type(parameters) :: params
        integer          :: nfiles, ldim(3), ifile, ifoo, cnt
        integer          :: lfoo(3), nimgs, iimg
        type(image)      :: img, tmp
        real             :: mm(2)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call read_filetable(params%filetab, filenames)
        nfiles = size(filenames)
        call find_ldim_nptcls(filenames(1),ldim,ifoo)
        ldim(3) = 1 ! to correct for the stupid 3:d dim of mrc stacks
        params%box = ldim(1)
        ! prepare img and tmp for reading
        call img%new(ldim, params%smpd)
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
                call img%read(filenames(ifile), iimg, readhead=.false.)
                if( cline%defined('clip') )then
                    call img%clip(tmp)
                    mm = tmp%minmax()
                    call tmp%write(params%outstk, cnt)
                else
                    call img%write(params%outstk, cnt)
                endif
            end do
            call progress(ifile, nfiles)
        end do
        ! cleanup
        call img%kill
        call tmp%kill
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
        if( .not. cline%defined('outfile') ) call cline%set('outfile', 'outfile.txt')
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
        ! rotationally average
        if( params%roavg .eq. 'yes' )then
            call roavg_imgfile(params%stk, params%outstk, params%angstep, params%smpd)
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
        if( params%makemovie .eq. 'yes' )then
            call prep_imgfile4movie(params%stk, params%smpd)
            goto 999
        endif
        ! default
        write(logfhandle,*)'Nothing to do!'
        ! end gracefully
999     call o_here%kill
        call os_ran%kill
        call o%kill
        call o2%kill
        call rt%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_STACKOPS NORMAL STOP ****')
    end subroutine exec_stackops

    subroutine exec_pspec_int_rank( self, cline )
        use simple_segmentation, only: otsu_robust_fast
        class(pspec_int_rank_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! varables
        type(parameters)     :: params
        type(image)          :: img, img_pspec, graphene_mask
        character(len=:), allocatable :: ranked_fname, good_fname, bad_fname
        real,             allocatable :: peakvals(:), peakvals_copy(:)
        integer,          allocatable :: order(:)
        logical,          allocatable :: good_bad_msk(:)
        integer :: i, cnt_good, cnt_bad, funit
        real    :: thresh(3), ave, sdev, minv, mskrad
        if( .not. cline%defined('lp_backgr') ) call cline%set('lp_backgr', 7.)
        call params%new(cline)
        mskrad = params%moldiam / params%smpd / 2.
        ranked_fname = add2fbody(trim(params%stk), '.'//fname2ext(trim(params%stk)), '_ranked')
        good_fname   = add2fbody(trim(params%stk), '.'//fname2ext(trim(params%stk)), '_high_quality')
        bad_fname    = add2fbody(trim(params%stk), '.'//fname2ext(trim(params%stk)), '_low_quality')
        call graphene_mask%pspec_graphene_mask([params%box,params%box,1], params%smpd)
        call graphene_mask%write('graphene_mask.mrc')
        call img%new([params%box,params%box,1], params%smpd)
        call img_pspec%new([params%box,params%box,1], params%smpd)
        allocate(peakvals(params%nptcls), order(params%nptcls))
        call fopen(funit, file='power_spectrum_stats.txt', status='replace')
        do i=1,params%nptcls
            call progress(i,params%nptcls)
            call img%read(params%stk, i)
            call img%norm
            call img%mask(mskrad, 'soft')
            call img%img2spec('sqrt', params%lp_backgr, img_pspec)
            call img_pspec%write('pspecs.mrc', i)
            call img_pspec%mul(graphene_mask)
            call img_pspec%write('pspecs_graphene_msk.mrc', i)
            call img_pspec%nlmean
            call img_pspec%write('pspecs_graphene_msk_nlmean.mrc', i)
            call img_pspec%stats(ave, sdev, peakvals(i), minv)
            write(funit, '(A,1X,F7.3,1X,F7.3,1X,F7.3,1X,F7.3)') 'AVE/SDEV/MAXV/MINV: ', ave, sdev, peakvals(i), minv
        end do
        call fclose(funit)
        ! write ranked cavgs
        allocate(peakvals_copy(params%nptcls), source=peakvals)
        order = (/(i,i=1,params%nptcls)/)
        call hpsort(peakvals, order)
        call reverse(order) ! largest first
        do i=1,params%nptcls
            call img%read(params%stk, order(i))
            call img%write(ranked_fname, i)
        end do
        ! make an automatic good/bad classification
        call otsu(peakvals_copy, good_bad_msk)
        cnt_good = 0
        cnt_bad  = 0
        do i=1,params%nptcls
            call img%read(params%stk, i)
            if( good_bad_msk(i) )then
                cnt_good = cnt_good + 1
                call img%write(good_fname, cnt_good)
            else
                cnt_bad = cnt_bad + 1
                call img%write(bad_fname, cnt_bad)
            endif
        end do
        if( allocated(good_bad_msk) ) deallocate(good_bad_msk)
        call img%kill
        call img_pspec%kill
        call graphene_mask%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PSPEC_INT_RANK NORMAL STOP ****')
    end subroutine exec_pspec_int_rank

    subroutine exec_estimate_diam( self, cline )
        use simple_segmentation, only: otsu_robust_fast
        class(estimate_diam_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! constants
        character(len=*), parameter :: FILT   = 'nlmean_filtered.mrc'
        character(len=*), parameter :: MASKED = 'masked.mrc'
        ! varables
        type(parameters)            :: params
        type(binimage), allocatable :: imgs_mask(:) ! images mask
        type(image),    allocatable :: imgs(:)      ! images
        type(image)                 :: roavg        ! rotational average
        type(stats_struct)          :: diamstats    ! stats struct
        real,           allocatable :: diams(:)     ! diameters
        integer :: funit, i, loc(1)
        real    :: med_diam, thresh(3), msk_rad
        if( .not. cline%defined('lp')    ) call cline%set('lp',     7.0)
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir','yes')
        call params%new(cline)
        ! set radius for hard mask of binary image
        if(cline%defined('msk') )  then
            msk_rad = params%msk
        else
            msk_rad = 0.45*params%box
        endif
        ! allocate & read cavgs
        allocate(imgs_mask(params%nptcls),diams(params%nptcls),imgs(params%nptcls))
        diams = 0.
        do i=1,params%nptcls
            call imgs_mask(i)%new_bimg([params%box,params%box,1],params%smpd)
            call imgs_mask(i)%read(params%stk, i)
            call imgs(i)%copy(imgs_mask(i))
        end do
        call roavg%new([params%box,params%box,1],params%smpd)
        ! prepare thread safe images in image class
        call imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        write(logfhandle,'(A)') '>>> ESTIMATING DIAMETERS THROUGH BINARY IMAGE PROCESSING'
        call fopen(funit, file='diameters_in_Angstroms.txt', status='replace')
        do i=1,params%nptcls
            call progress(i,params%nptcls)
            call imgs_mask(i)%zero_edgeavg
            call imgs_mask(i)%bp(0.,params%lp)
            ! non-local mneans filter for denoising
            call imgs_mask(i)%NLmean
            call imgs_mask(i)%write(FILT, i)
            ! rotational averaging
            if( params%roavg .eq. 'yes' )then
                call imgs_mask(i)%roavg(params%angstep, roavg)
                call imgs_mask(i)%copy(roavg)
            endif
            ! binarise with Otsu
            call otsu_robust_fast(imgs_mask(i), is2D=.true., noneg=.false., thresh=thresh)
            ! hard mask for removing noise, default diameter 90% of the box size
            call imgs_mask(i)%mask(msk_rad,'hard')
            call imgs_mask(i)%set_imat     ! integer matrix set in binimage instance
            call imgs_mask(i)%grow_bins(9) ! to compensate low-pass erosion/binarization
            call imgs_mask(i)%update_img_rmat
            call imgs_mask(i)%diameter_bin(diams(i))
            call imgs(i)%mul(imgs_mask(i))
            call imgs(i)%write(MASKED, i)
            ! estimate diameter
            diams(i) = diams(i) * params%smpd
            write(funit,'(F6.1)') diams(i)
        end do
        call fclose(funit)
        call calc_stats(diams, diamstats)
        ! output
        med_diam = median(diams)
        write(logfhandle,'(A,2F6.1)') '>>> AVG    DIAMETER (IN A & pix): ', diamstats%avg,  diamstats%avg/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> SDEV   DIAMETER (IN A & pix): ', diamstats%sdev, diamstats%sdev/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MEDIAN DIAMETER (IN A & pix): ', med_diam,       med_diam/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MAX    DIAMETER (IN A & pix): ', diamstats%maxv, diamstats%maxv/params%smpd
        write(logfhandle,'(A,2F6.1)') '>>> MIN    DIAMETER (IN A & pix): ', diamstats%minv, diamstats%minv/params%smpd
        call fopen(funit, file='diameter_stats.txt', status='replace')
        write(funit,     '(A,2F6.1)') '>>> AVG    DIAMETER (IN A & pix): ', diamstats%avg,  diamstats%avg/params%smpd
        write(funit,     '(A,2F6.1)') '>>> SDEV   DIAMETER (IN A & pix): ', diamstats%sdev, diamstats%sdev/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MEDIAN DIAMETER (IN A & pix): ', med_diam,       med_diam/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MAX    DIAMETER (IN A & pix): ', diamstats%maxv, diamstats%maxv/params%smpd
        write(funit,     '(A,2F6.1)') '>>> MIN    DIAMETER (IN A & pix): ', diamstats%minv, diamstats%minv/params%smpd
        call fclose(funit)
        ! destruct
        do i=1,size(imgs)
            call imgs_mask(i)%kill
            call imgs_mask(i)%kill_bimg
            call imgs(i)%kill
        end do
        call roavg%kill
        deallocate(imgs, diams, imgs_mask)
        ! end gracefully
        call simple_end('**** SIMPLE_ESTIMATE_DIAM NORMAL STOP ****')
    end subroutine exec_estimate_diam

end module simple_commander_imgproc
