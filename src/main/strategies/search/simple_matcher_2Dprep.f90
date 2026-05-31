!@descr: common routines used by the high-level strategy 2D and 3D matchers
module simple_matcher_2Dprep
use simple_pftc_srch_api
use simple_timer
use simple_builder,                 only: builder
use simple_butterworth,             only: butterworth_filter
use simple_matcher_ptcl_io,         only: prepimgbatch, killimgbatch
use simple_matcher_smpl_and_lplims, only: set_bp_range3D, set_bp_range2D
use simple_projector,               only: projector
use simple_strategy2D_utils,        only: calc_cavg_offset
implicit none

! Particle image processing for alignment
public :: prepimg4align, prepimg4align_bench
public :: prepimg4align_reset_ctf_model_audit, prepimg4align_get_nctf_models_seen
public :: prepimg4align_disable_ctf_model_audit
! Reference processing for alignment
public :: prep2Dref, calc_2Dref_offset
private
#include "simple_local_flags.inc"

logical                       :: prepimg4align_ctf_model_audit = .false.
integer                       :: prepimg4align_nctf_models_seen = 0
type(ctfparams), allocatable  :: prepimg4align_ctf_models_seen(:)

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
        if( prepimg4align_ctf_model_audit ) call audit_ctf_model(ctfparms)
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

    subroutine prepimg4align_bench( params, build, iptcl, img, img_out, img_out_pd, rt_prep1, rt_prep2, rt_prep )
        class(parameters),    intent(in)     :: params
        class(builder),       intent(inout) :: build
        integer,              intent(in)    :: iptcl
        class(image),         intent(inout) :: img
        class(image),         intent(inout) :: img_out
        class(image),         intent(inout) :: img_out_pd
        real(timer_int_kind), intent(inout) :: rt_prep1, rt_prep2, rt_prep
        integer(timer_int_kind) :: t_prep1, t_prep2, t_prep
        type(ctf)       :: tfun
        type(ctfparams) :: ctfparms
        real            :: shvec(2), crop_factor
        t_prep1 = tic()
        t_prep   = t_prep1
        crop_factor = real(params%box_crop) / real(params%box)
        shvec(1)    = -build%spproj_field%get(iptcl, 'x') * crop_factor
        shvec(2)    = -build%spproj_field%get(iptcl, 'y') * crop_factor
        ! Phase-flipping
        ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
        if( prepimg4align_ctf_model_audit ) call audit_ctf_model(ctfparms)
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
        rt_prep1 = rt_prep1 + toc(t_prep1)
        t_prep2 = tic()
        ! fused IFFT, mask, FFT
        call img_out%ifft_mask_pad_fft(params%msk_crop, img_out_pd)
        rt_prep2 = rt_prep2 + toc(t_prep2)
        rt_prep   = rt_prep   + toc(t_prep)
    end subroutine prepimg4align_bench

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

    subroutine prepimg4align_reset_ctf_model_audit(capacity)
        integer, optional, intent(in) :: capacity
        integer :: capacity_here
        capacity_here = 16
        if( present(capacity) ) capacity_here = max(1, capacity)
        prepimg4align_ctf_model_audit = .true.
        prepimg4align_nctf_models_seen = 0
        if( allocated(prepimg4align_ctf_models_seen) ) deallocate(prepimg4align_ctf_models_seen)
        allocate(prepimg4align_ctf_models_seen(capacity_here))
    end subroutine prepimg4align_reset_ctf_model_audit

    integer function prepimg4align_get_nctf_models_seen()
        prepimg4align_get_nctf_models_seen = prepimg4align_nctf_models_seen
    end function prepimg4align_get_nctf_models_seen

    subroutine prepimg4align_disable_ctf_model_audit()
        prepimg4align_ctf_model_audit = .false.
        prepimg4align_nctf_models_seen = 0
        if( allocated(prepimg4align_ctf_models_seen) ) deallocate(prepimg4align_ctf_models_seen)
    end subroutine prepimg4align_disable_ctf_model_audit

    subroutine audit_ctf_model(ctfparms)
        type(ctfparams), intent(in) :: ctfparms
        integer :: imodel
        logical :: model_seen
        if( .not. prepimg4align_ctf_model_audit ) return
        !$omp critical(prepimg4align_ctf_audit)
        if( .not. allocated(prepimg4align_ctf_models_seen) ) allocate(prepimg4align_ctf_models_seen(16))
        model_seen = .false.
        do imodel = 1,prepimg4align_nctf_models_seen
            if( same_ctf_model(ctfparms, prepimg4align_ctf_models_seen(imodel)) )then
                model_seen = .true.
                exit
            endif
        enddo
        if( .not. model_seen )then
            if( prepimg4align_nctf_models_seen == size(prepimg4align_ctf_models_seen) ) call grow_ctf_model_audit
            prepimg4align_nctf_models_seen = prepimg4align_nctf_models_seen + 1
            prepimg4align_ctf_models_seen(prepimg4align_nctf_models_seen) = ctfparms
        endif
        !$omp end critical(prepimg4align_ctf_audit)
    end subroutine audit_ctf_model

    subroutine grow_ctf_model_audit()
        type(ctfparams), allocatable :: tmp(:)
        integer :: oldsz, newsz
        oldsz = size(prepimg4align_ctf_models_seen)
        newsz = max(1, 2 * oldsz)
        allocate(tmp(newsz))
        tmp(1:oldsz) = prepimg4align_ctf_models_seen
        call move_alloc(tmp, prepimg4align_ctf_models_seen)
    end subroutine grow_ctf_model_audit

    logical function same_ctf_model(lhs, rhs)
        type(ctfparams), intent(in) :: lhs, rhs
        real, parameter :: CTFTOL = 1.0e-5
        same_ctf_model = (lhs%ctfflag == rhs%ctfflag) .and. &
            (lhs%l_phaseplate .eqv. rhs%l_phaseplate) .and. &
            (abs(lhs%kv    - rhs%kv)    <= CTFTOL) .and. &
            (abs(lhs%cs    - rhs%cs)    <= CTFTOL) .and. &
            (abs(lhs%fraca - rhs%fraca) <= CTFTOL)
    end function same_ctf_model

end module simple_matcher_2Dprep
