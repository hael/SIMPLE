!@descr: common routines used by the high-level strategy 2D and 3D matchers
module simple_matcher_2Dprep
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_pftc_srch_api
use simple_timer
use simple_builder,                 only: builder
use simple_butterworth,             only: butterworth_filter
use simple_matcher_ptcl_io,         only: prepimgbatch, killimgbatch
use simple_matcher_smpl_and_lplims, only: set_bp_range3D, set_bp_range2D
use simple_projector,               only: projector
use simple_strategy2D_utils,        only: calc_cavg_offset
use simple_memoize_ft_maps,         only: ft_map_get_farray_shape
implicit none

! Particle image processing for alignment
public :: prepimg4align, prepimg4align_cached, prepimg4align_cart
! Reference processing for alignment
public :: prep2Dref, calc_2Dref_offset
private
#include "simple_local_flags.inc"

contains

    !>  \ brief  prepares one particle image for alignment
    !          serial routine
    subroutine prepimg4align( params, build, iptcl, img, img_out, img_out_pd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        class(image),      intent(inout) :: img_out
        class(image),      intent(inout) :: img_out_pd
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        crop_factor = real(params%box_crop) / real(params%box)
        shvec(1)    = -build%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused noise normalization, FFT, clip & shift
                call img%norm_noise_fft_clip_shift(build%lmsk, img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused noise normalization, FFT, clip, shift & CTF flip
                call img%norm_noise_fft_clip_shift_ctf_flip(build%lmsk, img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align')
        end select
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params%msk_crop, img_out_pd)
    end subroutine prepimg4align

    !>  \ brief  prepares one cached particle image for alignment
    !          serial routine
    !>  img is a particle from the downscaled cache: already noise-normalized against
    !>  the full-box lmsk and already Fourier-cropped, so only the iteration-dependent
    !>  remainder is left to do. Every stage below is the same operation, in the same
    !>  order, as its counterpart in prepimg4align.
    subroutine prepimg4align_cached( params, build, iptcl, img, img_out, img_out_pd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        class(image),      intent(inout) :: img_out
        class(image),      intent(inout) :: img_out_pd
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        crop_factor = real(params%box_crop) / real(params%box)
        shvec(1)    = -build%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        select case(ctfparms%ctfflag)
            case(CTFFLAG_NO, CTFFLAG_FLIP)
                ! fused FFT & shift
                call img%fft_clip_shift(img_out, shvec)
            case(CTFFLAG_YES)
                ctfparms%smpd = ctfparms%smpd / crop_factor != smpd_crop
                tfun          = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                ! fused FFT, shift & CTF flip
                call img%fft_clip_shift_ctf_flip(img_out, shvec, tfun, ctfparms)
            case DEFAULT
                THROW_HARD('unsupported CTF flag: '//int2str(ctfparms%ctfflag)//' prepimg4align_cached')
        end select
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params%msk_crop, img_out_pd)
    end subroutine prepimg4align_cached

    !> The Cartesian sibling of prepimg4align (plan section 6.4): the full-disk observation of one
    !! particle as prepimg4align prepares the polar particle
    !! up to its mask (O4): noise normalization, Fourier clip and the shift by minus the stored
    !! shift (shift_crop, cropped pixels), the phase flip for CTFFLAG_YES at the cropped
    !! sampling, the soft mask at mskrad; without the padding. The masked image is then
    !! multiplied by the separable taper of the gather stencil (O5 (a); cartft_calc%get_ptcl_taper)
    !! before its transform. ctfparms_out carries the cropped sampling. For CTFFLAG_YES the
    !! caller memoizes the Fourier maps of the cropped box (memoize_ft_maps), as for
    !! prepimg4align; THROW_HARD when they are absent.
    subroutine prepimg4align_cart(raw_img, noise_mask, work_img, mskrad, smpd_crop, &
        &shift_crop, taper, ctfparms_in, observed, ctfparms_out)
        class(image), intent(inout) :: raw_img, work_img
        logical, intent(in) :: noise_mask(:, :, :)
        real, intent(in) :: mskrad, smpd_crop, shift_crop(2), taper(:)
        type(ctfparams), intent(in) :: ctfparms_in
        complex, allocatable, intent(out) :: observed(:, :)
        type(ctfparams), intent(out) :: ctfparms_out
        type(ctf) :: tfun
        real, pointer :: rmat(:, :, :)
        integer :: ldim(3), i, j

        ldim = work_img%get_ldim()
        if (ldim(3) /= 1 .or. ldim(1) < 2 .or. mod(ldim(1), 2) /= 0 .or. ldim(2) /= ldim(1)) &
            &THROW_HARD('pose_cont observation workspace must be an even square image')
        if (any(shape(noise_mask) /= raw_img%get_ldim())) &
            &THROW_HARD('pose_cont observation noise mask has incompatible dimensions')
        if (smpd_crop <= TINY .or. .not. ieee_is_finite(smpd_crop)) &
            &THROW_HARD('pose_cont observation requires positive finite sampling')
        if (size(taper) /= ldim(1)) THROW_HARD('pose_cont observation taper does not match the box')
        ctfparms_out = ctfparms_in
        ctfparms_out%smpd = smpd_crop
        select case (ctfparms_in%ctfflag)
        case (CTFFLAG_NO, CTFFLAG_FLIP)
            call raw_img%norm_noise_fft_clip_shift(noise_mask, work_img, -shift_crop)
        case (CTFFLAG_YES)
            if (any(ft_map_get_farray_shape() /= [ldim(1)/2 + 1, ldim(2), 1])) &
                &THROW_HARD('pose_cont observation requires the Fourier maps of the cropped box')
            tfun = ctf(smpd_crop, ctfparms_in%kv, ctfparms_in%cs, ctfparms_in%fraca)
            call raw_img%norm_noise_fft_clip_shift_ctf_flip(noise_mask, work_img, -shift_crop, tfun, ctfparms_out)
        case default
            THROW_HARD('unsupported CTF flag for pose_cont observation')
        end select
        call work_img%ifft_mask_fft(mskrad)
        ! the stencil's taper on the masked image, centred at ldim/2+1 as the image is after ifft
        call work_img%ifft()
        call work_img%get_rmat_ptr(rmat)
        do j = 1, ldim(2)
            do i = 1, ldim(1)
                rmat(i, j, 1) = rmat(i, j, 1)*taper(i)*taper(j)
            end do
        end do
        nullify (rmat)
        call work_img%fft()
        allocate (observed(-ldim(1)/2:ldim(1)/2, -ldim(2)/2:ldim(2)/2))
        observed = work_img%expand_ft()
    end subroutine prepimg4align_cart

    !>  \brief  Calculates the centering offset of the input cavg
    !>          cavg & particles centering is not performed
    subroutine calc_2Dref_offset( params, build, img, icls, centype, xyz )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        class(image),      intent(inout) :: img
        integer,           intent(in)    :: icls
        integer,           intent(in)    :: centype
        real,              intent(out)   :: xyz(3)
        real :: xy_cavg(2), crop_factor
        crop_factor = real(params%box_crop) / real(params%box)
        select case(centype)
            case(PARAMS_CEN)
                call build%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                if( arg(xy_cavg) < CENTHRESH )then
                    xyz = 0.
                else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                    xyz(1:2) = xy_cavg * crop_factor
                    xyz(3)   = 0.
                else
                    xyz = img%calc_shiftcen_serial(params%cenlp, params%msk_crop)
                    if( arg(xyz(1:2)/crop_factor - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                endif
            case(SEG_CEN)
                call calc_cavg_offset(img, params%cenlp, params%msk_crop, xy_cavg)
                xyz = [xy_cavg(1), xy_cavg(2), 0.]
            case(MASS_CEN)
                xyz = img%calc_shiftcen_serial(params%cenlp, params%msk_crop)
        end select
        if( arg(xyz) < CENTHRESH ) xyz = 0.0
    end subroutine calc_2Dref_offset

    !>  \brief  Prepares one cluster centre image for alignment
    subroutine prep2Dref( params, build, img, icls, xyz, img_pd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        class(image),      intent(inout) :: img, img_pd
        integer,           intent(in)    :: icls
        real,              intent(in)    :: xyz(3)
        integer :: filtsz
        real    :: frc(img%get_filtsz()), filter(img%get_filtsz()), crop_factor
        ! Cavg & particles centering
        if( arg(xyz) > CENTHRESH )then
            call img%fft()
            call img%shift2Dserial(xyz(1:2))
            crop_factor = real(params%box_crop) / real(params%box)
            call build%spproj_field%add_shift2class(icls, -xyz(1:2) / crop_factor)
        endif
        ! Filtering
        if( params%l_ml_reg )then
            ! no filtering
        else
            if( params%l_lpset.and.params%l_gauref )then
                ! Gaussian filter only applied when lp is set, FRC filtering turned off
                call img%fft
                call img%bpgau2d(0., params%gaufreq)
            else
                ! FRC-based filtering
                call build%clsfrcs%frc_getter(icls, frc)
                if( any(frc > 0.143) )then
                    filtsz = img%get_filtsz()
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params%l_lpset)
                    call img%fft    ! in case the shift was not applied above
                    call img%apply_filter_serial(filter)
                endif
            endif
        endif
        ! ensure we are in real-space
        call img%ifft()
        ! mask, pad, FFT
        call img%mask_pad_fft(params%msk_crop, img_pd)
    end subroutine prep2Dref

end module simple_matcher_2Dprep
