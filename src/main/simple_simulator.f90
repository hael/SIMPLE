!@descr: simulation of single-particle images
module simple_simulator
use simple_core_module_api
implicit none

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
        real            :: dfx, dfy, angast
        aapply_ctf = .true.
        if( present(apply_ctf) ) aapply_ctf = apply_ctf
        ! back FT (to make sure)
        call img%ifft()
        call img%fft()
        ctfparms = orientation%get_ctfvars()
        call img%apply_ctf(tfun, 'ctf', ctfparms )
        if( present(bfac) ) call img%apply_bfac(bfac)
        ! add detector noise
        call img%ifft()
        if( snr < 5 )then
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

end module simple_simulator
