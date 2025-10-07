! operations on micrographs
module simple_micproc
include 'simple_lib.f08'
use simple_image,    only: image
use simple_binimage, only: binimage
implicit none
#include "simple_local_flags.inc"

contains

    subroutine read_mic_subtr_backgr_shrink( micname, smpd, scale, pcontrast, mic_raw, mic_shrink )
        character(len=*), intent(in)    :: micname !< micrograph file name
        real,             intent(in)    :: smpd    !< sampling distance in A
        real,             intent(in)    :: scale   !< scale factor
        character(len=*), intent(in)    :: pcontrast
        class(image),     intent(inout) :: mic_raw, mic_shrink
        integer :: nframes, ldim(3), ldim_shrink(3)
        ! set micrograph info
        call find_ldim_nptcls(micname, ldim, nframes)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        ! set shrunked dims
        ldim_shrink(1) = round2even(real(ldim(1)) * scale)
        ldim_shrink(2) = round2even(real(ldim(2)) * scale)
        ldim_shrink(3) = 1
        ! make shrunken micrograph
        call mic_shrink%new(ldim_shrink, smpd/scale)
        ! read micrograph
        call mic_raw%new(ldim, smpd)
        call mic_raw%read(micname)
        call mic_raw%subtract_background(HP_BACKGR_SUBTR)
        call mic_raw%fft
        select case(trim(pcontrast))
            case('black')
                ! flip contrast (assuming black particle contrast on input)
                call mic_raw%mul(-1.)
            case('white')
                ! nothing to do
            case DEFAULT
                THROW_HARD('uknown pcontrast parameter, use (black|white)')
        end select
        call mic_raw%mul(real(product(ldim))) ! to prevent numerical underflow when performing FFT
        ! shrink micrograph
        call mic_shrink%set_ft(.true.)
        call mic_raw%clip(mic_shrink)
        call mic_raw%ifft
        call mic_raw%div(real(product(ldim)))
        call mic_shrink%ifft
    end subroutine read_mic_subtr_backgr_shrink

    subroutine read_mic( micname, mic_out )
        character(len=*), intent(in)    :: micname !< micrograph file name
        class(image),     intent(inout) :: mic_out
        integer :: nframes, ldim(3)
        real    :: smpd
        call find_ldim_nptcls(micname, ldim, nframes, smpd=smpd)
        if( ldim(3) /= 1 .or. nframes /= 1 ) THROW_HARD('Only for 2D images')
        call mic_out%new(ldim, smpd)
        call mic_out%read(micname)
    end subroutine read_mic

    subroutine cascade_filter_biomol( mic_shrink )
        class(image), intent(inout) :: mic_shrink
        real,    parameter :: SMPD_SHRINK1  = 4.0, LP_UB = 15., LAM_ICM = 100., DAMP = 10., FRAC_FG = 0.17, LAM_TV = 7.
        integer, parameter :: WINSZ_MED = 3
        call cascade_filter(mic_shrink, DAMP, LP_UB, LAM_TV, LAM_ICM, WINSZ_MED)
    end subroutine cascade_filter_biomol

    subroutine cascade_filter( mic_shrink, damp_below_zero, lp, lam_tv, lam_icm, winsz_med )
        use simple_tvfilter, only: tvfilter
        class(image), intent(inout) :: mic_shrink
        real,         intent(in)    :: damp_below_zero, lp, lam_tv, lam_icm
        integer,      intent(in)    :: winsz_med
        type(tvfilter) :: tvf
        call mic_shrink%zero_edgeavg
        ! dampens below zero
        call mic_shrink%div_below(0.,damp_below_zero)
        ! low-pass filter
        call mic_shrink%bp(0.,lp)
        ! TV denoising
        call tvf%new()
        call tvf%apply_filter(mic_shrink, lam_tv)
        call tvf%kill
        ! Non-local-means denoising
        call mic_shrink%NLmean2D
        ! Iterated conditional modes denoising
        call mic_shrink%ICM2D(lam_icm)
        ! Median filter
        call mic_shrink%real_space_filter(winsz_med, 'median')
    end subroutine cascade_filter

    subroutine binarize_mic_den( mic_den, frac_fg, mic_bin )
        class(image),    intent(in)  :: mic_den
        real,            intent(in)  :: frac_fg
        class(binimage), intent(out) :: mic_bin
        real :: bin_t
        call mic_bin%transfer2bimg(mic_den)
        call mic_bin%calc_bin_thres(frac_fg, bin_t)
        call mic_bin%binarize(bin_t)
        call mic_bin%set_imat            
        call mic_bin%erode() ! -4 A
        call mic_bin%erode() ! -4 A
        call mic_bin%set_largestcc2background
        call mic_bin%inv_bimg()
    end subroutine binarize_mic_den

    subroutine identify_masscens( mic_bin, masscens )
        use simple_segmentation
        class(binimage),   intent(inout) :: mic_bin
        real, allocatable, intent(inout) :: masscens(:,:)
        type(binimage) :: img_cc
        integer :: i, nccs
        real    :: px(3)
        ! identify connected components
        call mic_bin%find_ccs(img_cc)
        call img_cc%get_nccs(nccs)
        if( allocated(masscens) ) deallocate(masscens)
        allocate(masscens(nccs,2), source=0.)
        do i = 1, nccs
            px = center_mass_cc(i)
            masscens(i,:2) = px(:2)
        enddo


        contains

            function center_mass_cc( i_cc ) result( px )
                integer, intent(in) :: i_cc
                real :: px(3)
                integer, allocatable :: pos(:,:)
                integer, allocatable :: imat_cc(:,:,:)
                call img_cc%get_imat(imat_cc)
                where( imat_cc .ne. i_cc ) imat_cc = 0
                call get_pixel_pos(imat_cc,pos)
                px(1) = sum(pos(1,:))/real(size(pos,dim = 2))
                px(2) = sum(pos(2,:))/real(size(pos,dim = 2))
                px(3) = 1.
                if( allocated(imat_cc) ) deallocate(imat_cc)
            end function center_mass_cc

    end subroutine identify_masscens

end module simple_micproc
