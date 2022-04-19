! simulation of single-particle images
module simple_simulator
implicit none

contains

    subroutine simimg( img, orientation, tfun, ctfflag, snr, bfac, apply_ctf )
        use simple_image, only: image
        use simple_ori,   only: ori
        use simple_ctf,   only: ctf
        class(image),      intent(inout) :: img
        class(ori),        intent(inout) :: orientation
        class(ctf),        intent(inout) :: tfun
        character(len=*),  intent(in)    :: ctfflag
        real,              intent(in)    :: snr
        real,    optional, intent(in)    :: bfac
        logical, optional, intent(in)    :: apply_ctf
        logical :: aapply_ctf
        real :: dfx, dfy, angast
        aapply_ctf = .true.
        if( present(apply_ctf) ) aapply_ctf = apply_ctf
        ! back FT (to make sure)
        call img%ifft()
        call img%fft()
        ! apply ctf/bfactor
        if( orientation%isthere('dfx') .and. orientation%isthere('dfy') .and. aapply_ctf )then
            dfx = orientation%get_dfx()
            dfy = orientation%get_dfy()
            angast = orientation%get('angast')
            call tfun%apply(img, dfx, 'ctf', dfy, angast, bfac=bfac)
        else if( orientation%isthere('dfx') .and. aapply_ctf )then
            dfx = orientation%get_dfx()
            dfy = dfx
            angast = 0.
            call tfun%apply(img, orientation%get_dfx(), 'ctf', bfac=bfac)
        else
            if( present(bfac) ) call img%apply_bfac(bfac)
        endif
        ! add detector noise
        call img%ifft()
        call img%add_gauran(snr)
        if( .not. aapply_ctf ) return
        if( ctfflag .eq. 'flip' )then
            ! simulate phase-flipped images
            call img%fft()
            call tfun%apply(img, dfx, 'flip', dfy, angast)
            call img%ifft()
        else if( ctfflag .eq. 'mul' )then
            ! simulate CTF multiplied images
            call img%fft()
            call tfun%apply(img, dfx, 'ctf', dfy, angast)
            call img%ifft()
        endif
    end subroutine simimg

end module simple_simulator
