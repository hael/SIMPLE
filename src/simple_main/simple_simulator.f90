module simple_simulator
implicit none

contains

    subroutine simimg( img, orientation, tfun, ctfflag, snr, snr_pink, snr_detector, bfac )
        use simple_image, only: image
        use simple_ori,   only: ori
        use simple_ctf,   only: ctf
        class(image),     intent(inout) :: img
        class(ori),       intent(inout) :: orientation
        class(ctf),       intent(inout) :: tfun
        character(len=*), intent(in)    :: ctfflag
        real,             intent(in)    :: snr, snr_pink, snr_detector
        real, optional,   intent(in)    :: bfac
        real :: dfx, dfy, angast
        ! back FT (to make sure)
        call img%bwd_ft
        ! add pink noise
        if( snr < 3. ) call img%add_gauran(snr_pink)
        call img%fwd_ft
        ! apply ctf/bfactor
        if( orientation%isthere('dfx') .and. orientation%isthere('dfy') )then
            dfx = orientation%get('dfx')
            dfy = orientation%get('dfy')
            angast = orientation%get('angast')
            if( present(bfac) )then
                call tfun%apply(img, dfx, 'ctf', dfy, angast, bfac=bfac)
            else
                call tfun%apply(img, dfx, 'ctf', dfy, angast)
            endif
        else if( orientation%isthere('dfx') )then
            dfx = orientation%get('dfx')
            dfy = dfx
            angast = 0.
            if( present(bfac) )then
                call tfun%apply(img, orientation%get('dfx'), 'ctf', bfac=bfac)
            else
                call tfun%apply(img, orientation%get('dfx'), 'ctf')
            endif
        else
            if( present(bfac) ) call img%apply_bfac(bfac)
        endif
        ! add detector noise
        call img%bwd_ft
        if( snr < 3. ) call img%add_gauran(snr_detector)
        if( ctfflag .eq. 'flip' )then
            ! simulate phase-flipped images
            call img%fwd_ft
            call tfun%apply(img, dfx, 'flip', dfy, angast)
            call img%bwd_ft
        else if( ctfflag .eq. 'mul' )then
            ! simulate CTF multiplied images
            call img%fwd_ft
            call tfun%apply(img, dfx, 'ctf', dfy, angast)
            call img%bwd_ft
        endif
    end subroutine simimg

end module simple_simulator
