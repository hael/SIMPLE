!@descr: simulation of single-particle images
module simple_simulator
use simple_core_module_api
implicit none
#include "simple_local_flags.inc"

contains

    subroutine simimg( img, orientation, tfun, ctfflag, snr, bfac, apply_ctf )
        use simple_image, only: image
        use simple_ctf,   only: ctf
        class(image),      intent(inout) :: img
        class(ori),        intent(inout) :: orientation
        class(ctf),        intent(inout) :: tfun
        character(len=*),  intent(in)    :: ctfflag
        real,              intent(in)    :: snr
        real,    optional, intent(in)    :: bfac
        logical, optional, intent(in)    :: apply_ctf
        type(ctfparams) :: ctfparms
        logical         :: aapply_ctf
        aapply_ctf = .true.
        if( present(apply_ctf) ) aapply_ctf = apply_ctf
        ! back FT (to make sure)
        call img%ifft()
        call img%fft()
        if( aapply_ctf )then
            ctfparms = orientation%get_ctfvars()
            call img%apply_ctf(tfun, 'ctf', ctfparms )
        endif
        if( present(bfac) ) call img%apply_bfac(bfac)
        ! add detector noise
        call img%ifft()
        if( snr < 5. )then
            call img%add_gauran(snr)
        endif
        if( .not. aapply_ctf ) return
        if( ctfflag .eq. 'flip' )then
            ! simulate phase-flipped images
            call img%fft()
            call img%apply_ctf(tfun, 'flip', ctfparms )
            call img%ifft()
        else if( ctfflag .eq. 'mul' )then
            ! simulate CTF multiplied images
            call img%fft()
            call img%apply_ctf(tfun, 'ctf', ctfparms )
            call img%ifft()
        endif
    end subroutine simimg

    !> Wave mode's detector on img, the relative intensity I of the native box: shot noise at the dose on the
    !! counts D a^2 I (unless shot_noise is false), Gaussian detector noise of standard deviation sqrt(var / snr)
    !! for snr < 5 (add_gauran's rule), then SIMPLE's particles: x = -(counts - D a^2) / sqrt(D a^2)
    subroutine simimg_wave( img, total_dose, snr, shot_noise )
        use simple_image,        only: image
        use simple_wave_imaging, only: detect_intensity
        class(image),      intent(inout) :: img
        real,              intent(in)    :: total_dose, snr
        logical, optional, intent(in)    :: shot_noise
        real, allocatable :: rmat(:,:,:)
        real :: dose_px, sdev, var, dose
        if( total_dose <= 0. ) THROW_HARD('total_dose must be positive; simimg_wave')
        dose_px = total_dose * img%get_smpd()**2
        ! the noise of a relative intensity: the counts' standard deviation over D a^2
        sdev = 0.
        if( snr < 5. )then
            var = img%variance()
            if( var <= TINY ) THROW_HARD('variance of image is zero; simimg_wave')
            sdev = sqrt(var / snr)
        endif
        dose = total_dose
        if( present(shot_noise) )then
            if( .not. shot_noise ) dose = 0.
        endif
        call detect_intensity(img, dose, sdev)
        rmat = img%get_rmat()
        rmat = sqrt(dose_px) * (1. - rmat)
        call img%set_rmat(rmat, .false.)
    end subroutine simimg_wave

end module simple_simulator
