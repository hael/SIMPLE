module simple_estimate_ssnr
use simple_defs
use simple_image, only: image
implicit none

contains
    
    ! ESTIMATION ROUTINES FOR THE 2D ANALYSIS
    
    !>  \brief make a fake ssnr estimate, useful for testing purposes
    function fake_ssnr( img, lplim, width ) result( ssnr )
        class(image), intent(in)   :: img
        real, intent(in)           :: lplim
        real, intent(in), optional :: width
        real, allocatable :: ssnr(:), fsc(:), res(:)
        integer           :: nyq, k
        real              :: freq, wwidth, lplim_freq
        wwidth = 7.
        if( present(width) ) wwidth = width
        nyq = img%get_nyq()
        allocate(fsc(nyq), res(nyq))
        lplim_freq = real(img%get_find(lplim)) ! assuming square 4 now
        fsc = 1.0
        do k=1,nyq
            freq = real(k)
            res(k) = img%get_lp(k)
            if(freq .gt. lplim_freq)then
                fsc(k) = 0.
            else if(freq .ge. lplim_freq-wwidth)then
                fsc(k) = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
            endif
        end do
        ssnr = fsc2ssnr(fsc)
        deallocate(fsc,res)
    end function fake_ssnr
    
    !>  \brief estimate the SSNR using Unser's method
    !!         this routine is too slow in practice but it was useful for prototyping
    function estimate_ssnr( imgs, os, msk, tfun, tfplan ) result( ssnr )
        use simple_oris, only: oris
        use simple_ctf, only: ctf
        class(image),     intent(inout) :: imgs(:)
        class(oris),      intent(inout) :: os
        real,             intent(in)    :: msk
        class(ctf),       intent(inout) :: tfun
        type(ctfplan),    intent(in)    :: tfplan
        type(image)       :: favg, fdiff
        integer           :: ldim(3), iptcl
        real, allocatable :: spec(:), specsum(:), ssnr(:), specdiff(:), specdiffsum(:), s(:)
        real              :: dfx, dfy, angast, smpd, rnimgs, x, y, sc
        ! make Fourier average
        ldim   = imgs(1)%get_ldim()
        smpd   = imgs(1)%get_smpd()
        call favg%new(ldim, smpd)
        favg   = cmplx(0.,0.)
        rnimgs = real(size(imgs))
        ! FT & shellnorm & apply CTF
        do iptcl=1,size(imgs)
            call imgs(iptcl)%fwd_ft
            call imgs(iptcl)%shellnorm
            ! set CTF parameters
            select case(tfplan%mode)
                case('astig') ! astigmatic CTF
                    dfx    = os%get(iptcl,'dfx')
                    dfy    = os%get(iptcl,'dfy')
                    angast = os%get(iptcl,'angast')
                case('noastig') ! non-astigmatic CTF
                    dfx    = os%get(iptcl,'dfx')
                    dfy    = dfx
                    angast = 0.
            end select
            select case(tfplan%flag)
                case('mul','no') 
                    ! do nothing
                case('yes')  ! multiply with CTF
                    call tfun%apply(imgs(iptcl), dfx, 'ctf', dfy, angast)
                case('flip') ! multiply with abs(CTF)
                    call tfun%apply(imgs(iptcl), dfx, 'abs', dfy, angast)
                case DEFAULT
                    write(*,*) 'Unsupported ctfflag: ', tfplan%flag
                    stop 'simple_filterer :: estimate_ssnr'
            end select
            x = os%get(iptcl,'x')
            y = os%get(iptcl,'y')
            call imgs(iptcl)%shift(-x, -y)
            call imgs(iptcl)%bwd_ft
            call imgs(iptcl)%rtsq(-os%e3get(iptcl), 0., 0.)
            call imgs(iptcl)%mask(msk, 'soft')
            call imgs(iptcl)%fwd_ft
            call favg%add(imgs(iptcl))
        end do
        call favg%div(rnimgs)
        ! assemble spectra
        do iptcl=1,size(imgs)
            spec = imgs(iptcl)%spectrum('power')
            if( allocated(specsum) )then
                specsum = specsum+spec
            else
                allocate(specsum(size(spec)), source=spec)
            endif
            fdiff = imgs(iptcl)
            call fdiff%subtr(favg)
            specdiff = fdiff%spectrum('power')
            if( allocated(specdiffsum) )then
                specdiffsum = specdiffsum+specdiff
            else
                allocate(specdiffsum(size(specdiff)), source=specdiff)
            endif
            deallocate(spec,specdiff)
        end do
        ! calculate ssnr
        allocate(s(size(specsum)), ssnr(size(specsum)))
        s    = 0.
        ssnr = 0.
        sc   = rnimgs/(rnimgs-1.)
        where( abs(specdiffsum) > 1e-6 )
            s = specsum/(specdiffsum*sc)
        end where
        where( s > 1.)
            ssnr = s-1.
        end where
        deallocate(s,specsum,specdiffsum)
        call fdiff%kill
        call favg%kill
    end function estimate_ssnr
    
    !>  \brief estimate the noise power spectrum in an highly optimised online fashion
    !!         it is assumed that img is rotated, shifted and multiplied with the CTF
    !!         Hence, the reference should be multiplied with the square of the CTF
    subroutine estimate_specnoise_online( avg, img, msk, o, tfun, tfplan, specnoisesum, inner_width )
        use simple_ori, only: ori
        use simple_ctf, only: ctf
        class(image),     intent(inout) :: avg
        class(image),     intent(inout) :: img
        real,             intent(in)    :: msk
        class(ori),       intent(inout) :: o
        class(ctf),       intent(inout) :: tfun
        type(ctfplan),    intent(in)    :: tfplan
        real,             intent(inout) :: specnoisesum(:)
        real, optional,   intent(in)    :: inner_width(2)
        type(image)       :: fdiff, favg, fimg
        real, allocatable :: specnoise(:)
        real              :: dfx, dfy, angast
        if( img%is_ft() )&
        &stop 'input particle image (img) should not be FTed; simple_filterer :: estimate_specnoise_online'
        if( avg%is_ft() )&
        &stop 'input average (avg) should not be FTed; simple_filterer :: estimate_specnoise_online'
        ! copy inputs
        favg = avg
        fimg = img
        ! mask
        if( present(inner_width) )then
            call favg%mask(msk, 'soft', inner=inner_width(1), width=inner_width(2))
            call fimg%mask(msk, 'soft', inner=inner_width(1), width=inner_width(2))
        else
            call favg%mask(msk, 'soft')
            call fimg%mask(msk, 'soft')
        endif
        ! FT
        call favg%fwd_ft
        call fimg%fwd_ft
        ! deal with CTF
        select case(tfplan%mode)
            case('astig') ! astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = o%get('dfy')
                angast = o%get('angast')
            case('noastig') ! non-astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = dfx
                angast = 0.
        end select
        if( tfplan%flag .ne. 'no' ) call tfun%apply(img, dfx, 'square', dfy, angast)
        ! calculate difference image (fimg-favg)
        fdiff = fimg
        call fdiff%subtr(favg)
        ! calculate noise spectrum
        specnoise    = fdiff%spectrum('power')
        specnoisesum = specnoisesum+specnoise
        ! kill
        deallocate(specnoise)
        call favg%kill
        call fimg%kill
        call fdiff%kill        
    end subroutine estimate_specnoise_online
    
    ! ESTIMATION ROUTINES FOR THE 3D ANALYSIS
    
    !> \brief  converts the FSC to SSNR (the 2.* is because of the division of the data)
    function fsc2ssnr( corrs ) result( ssnr )
        real, intent(in)  :: corrs(:)
        real, allocatable :: ssnr(:)
        integer :: nyq, k
        real    :: fsc
        nyq = size(corrs)
        allocate( ssnr(nyq) )
        do k=1,nyq
            fsc = min(abs(corrs(k)),0.999)
            ssnr(k) = (2.*fsc)/(1.-fsc)
        end do
    end function fsc2ssnr
    
    !> \brief  converts the FSC to the optimal low-pass filter
    function fsc2optlp( corrs ) result( filt )
        real, intent(in)           :: corrs(:) !< fsc plot (correlations)
        real, allocatable          :: filt(:)  !< output filter coefficients
        integer :: nyq, k
        nyq = size(corrs)
        allocate( filt(nyq) )
        filt = 0.
        where( corrs > 0. ) filt = sqrt( 2. * corrs / (corrs + 1.) )
        where( filt  > 0.9999 ) filt = 0.9999
    end function fsc2optlp

    !> \brief  converts the SSNR to FSC
    function ssnr2fsc( ssnr ) result( corrs )
        real, intent(in)  :: ssnr(:)
        real, allocatable :: corrs(:)
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( corrs(nyq) )
        do k=1,nyq
            corrs(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2fsc
    
    !> \brief  converts the SSNR 2 the optimal low-pass filter
    function ssnr2optlp( ssnr ) result( w )
        real, intent(in)  :: ssnr(:)
        real, allocatable :: w(:)
        integer :: nyq, k
        nyq = size(ssnr)
        allocate( w(nyq) )
        do k=1,nyq
            w(k) = ssnr(k)/(ssnr(k)+1.)
        end do
    end function ssnr2optlp
    
    !> \brief  calculates the particle SSNR in 2D (Grigorieff)
    function estimate_pssnr2D( avr, fsc ) result( pssnr )
        use simple_AVratios, only: AVratios
        class(AVratios), intent(in) :: avr
        real, intent(in)            :: fsc(:)
        real, allocatable :: pssnr(:)
        pssnr = fsc2ssnr(fsc)
        pssnr = pssnr*avr%Abox_o_Aptcl()*avr%Vmsk_o_Vbox()
    end function estimate_pssnr2D
    
    !> \brief  calculates the particle SSNR in 3D (Grigorieff)
    function estimate_pssnr3D( avr, fsc ) result( pssnr )
        use simple_AVratios, only: AVratios
        class(AVratios), intent(in) :: avr
        real, intent(in)            :: fsc(:)
        real, allocatable :: pssnr(:)
        pssnr = estimate_pssnr2D(avr, fsc)
        pssnr = pssnr*avr%Vbox_o_Vptcl()*avr%Aptcl_o_Abox()
    end function estimate_pssnr3D
    
    !> \brief  calculates the particle SSNR in 2D (Grigorieff)
    function from_ssnr_estimate_pssnr2D( avr, ssnr ) result( pssnr )
        use simple_AVratios, only: AVratios
        class(AVratios), intent(in) :: avr
        real, intent(in)            :: ssnr(:)
        real, allocatable :: pssnr(:)
        allocate(pssnr(size(ssnr)), source=ssnr)
        pssnr = pssnr*avr%Abox_o_Aptcl()*avr%Vmsk_o_Vbox()
    end function from_ssnr_estimate_pssnr2D
    
    !> \brief  calculates the particle SSNR in 3D (Grigorieff)
    function from_ssnr_estimate_pssnr3D( avr, ssnr ) result( pssnr )
        use simple_AVratios, only: AVratios
        class(AVratios), intent(in) :: avr
        real, intent(in)            :: ssnr(:)
        real, allocatable :: pssnr(:)
        pssnr = from_ssnr_estimate_pssnr2D(avr, ssnr)
        pssnr = pssnr*avr%Vbox_o_Vptcl()*avr%Aptcl_o_Abox()
    end function from_ssnr_estimate_pssnr3D

end module simple_estimate_ssnr          
        