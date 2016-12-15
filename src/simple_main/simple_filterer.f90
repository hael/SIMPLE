module simple_filterer
use simple_defs       ! singleton
use simple_image,     only: image
use simple_projector, only: projector
implicit none

real, parameter :: SHTHRESH=0.0001

contains

    ! DOSE FILTERING (Grant, Grigorieff eLife 2015)

    !>  \brief input is template image, accumulative dose (in e/A2) and acceleration voltage
    !!         output is filter coefficients
    function acc_dose2filter( img, acc_dose, kV ) result( filter )
        type(image), intent(in) :: img
        real,        intent(in) :: acc_dose, kV
        real, allocatable       :: filter(:)
        integer :: find, sz
        sz = img%get_filtsz()
        allocate(filter(sz))
        do find=1,sz
            filter(find) = dose_weight(acc_dose, img%get_spat_freq(find), kV)
        end do
    end function acc_dose2filter

    !>  \brief input is accumulative dose (in e/A2) and spatial frequency (in 1/A)
    !!         output is resolution dependent weight applied to individual frames
    !!         before correlation search and averaging
    real function dose_weight( acc_dose, spat_freq, kV )
        real, intent(in) :: acc_dose, spat_freq, kV
        real, parameter  :: A=0.245, B=-1.665, C=2.81, kV_factor=0.75
        real :: critical_exp ! critical exposure (only depends on spatial frequency)
        critical_exp = A*(spat_freq**B)+C
        if( abs(kV-300.) < 0.001 )then
            ! critical exposure does not need modification
        else if( abs(kV-200.) < 0.001 )then
            ! critical exposure at 200 kV expected to be ~25% lower
            critical_exp = critical_exp*kV_factor
        else
            stop 'unsupported kV (acceleration voltage); simple_filterer :: dose_weight'
        endif
        dose_weight = exp(-acc_dose/(2.0*critical_exp))
    end function dose_weight

    ! SHELL-WEIGHT ROUTINES

    !> \brief  read in shellweights from file
    subroutine read_shellweights( img, img_pad, nptcls, doshellweight, res, res_pad, wmat )
        use simple_jiffys, only: file_exists, get_fileunit
        class(image),      intent(in)  :: img, img_pad
        integer,           intent(in)  :: nptcls
        logical,           intent(out) :: doshellweight
        real, allocatable, intent(out) :: res(:), res_pad(:), wmat(:,:)
        integer :: filtsz, filnum, io_stat, alloc_stat
        if( file_exists('cont3D_shellweights.bin') )then
            if( allocated(res)     ) deallocate(res)
            if( allocated(res_pad) ) deallocate(res_pad)
            if( allocated(wmat)    ) deallocate(wmat)
            filtsz  = img%get_filtsz()
            res     = img%get_res()
            res_pad = img_pad%get_res()
            allocate( wmat(nptcls,filtsz), stat=alloc_stat)
            filnum = get_fileunit()
            open(unit=filnum, status='OLD', action='READ', file='cont3D_shellweights.bin', access='STREAM')
            read(unit=filnum,pos=1,iostat=io_stat) wmat
            ! check if the read was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(eorec): I/O error ',&
                io_stat, ' when reading cont3D_shellweights.bin'
                stop 'I/O error; read_shellweights; simple_filterer'
            endif
            close(filnum)
            doshellweight = .true.
        else
            doshellweight = .false.
        endif
    end subroutine read_shellweights

    !> \brief  re-samples a filter array
    function resample_filter( filt_orig, res_orig, res_new ) result( filt_resamp )
        use simple_math, only: find
        real, intent(in)  :: filt_orig(:), res_orig(:), res_new(:)
        real, allocatable :: filt_resamp(:)
        integer :: filtsz_orig, filtsz_resamp, k, ind
        real    :: dist
        filtsz_orig   = size(filt_orig)
        filtsz_resamp = size(res_new)
        allocate(filt_resamp(filtsz_resamp))
        do k=1,filtsz_resamp
            call find(res_orig, filtsz_resamp, res_new(k), ind, dist)
            filt_resamp(k) = filt_orig(ind)
        end do
    end function resample_filter
    
    ! WIENER RESTORATION ROUTINES
    
    !>  \brief does the Wiener restoration of aligned images in 2D 
    subroutine wiener_restore2D( img_set, o_set, tfun, tfplan, img_rec, msk, shellw )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        use simple_ctf,   only: ctf
        class(image),     intent(inout) :: img_set(:)
        class(oris),      intent(inout) :: o_set
        class(ctf),       intent(inout) :: tfun
        type(ctfplan),    intent(in)    :: tfplan
        class(image),     intent(out)   :: img_rec
        real,             intent(in)    :: msk
        real, optional,   intent(in)    :: shellw(:,:)
        integer           :: ldim(3), nimgs, iptcl, k
        type(ori)         :: o
        type(image)       :: ctfsqsum
        real              :: smpd
        logical           :: doshellw
        if( o_set%get_noris() /= size(img_set) )&
        stop 'nr of imgs and oris not consistent; simple_filterer :: wiener_restore2D_1'
        doshellw = present(shellw)
        ! set constants
        ldim  = img_set(1)%get_ldim()
        smpd  = img_set(1)%get_smpd()
        nimgs = size(img_set)
        ! create & init objs
        call img_rec%new(ldim, smpd)
        call ctfsqsum%new(ldim, smpd)
        ctfsqsum = cmplx(0.,0.)
        ! average in the assumption of infinite signal
        do iptcl=1,nimgs
            o = o_set%get_ori(iptcl)
            if( doshellw )then
                call wiener_restore2D_online(img_set(iptcl), o, tfun, tfplan, img_rec, msk, shellw(iptcl,:))
                call assemble_ctfsqsum_online(img_set(iptcl), o, tfun, tfplan, ctfsqsum, shellw(iptcl,:))
            else
                call wiener_restore2D_online(img_set(iptcl), o, tfun, tfplan, img_rec, msk)
                call assemble_ctfsqsum_online(img_set(iptcl), o, tfun, tfplan, ctfsqsum)
            endif
        end do
        ! do the density correction
        call img_rec%fwd_ft
        call img_rec%ctf_dens_correct(ctfsqsum)
        call img_rec%bwd_ft
        ! destroy objects
        call ctfsqsum%kill
    end subroutine wiener_restore2D

    !>  \brief does the online Wiener restoration of 2D images, including shift+rotations
    !!         the image is left shifted and Fourier transformed on output
    subroutine wiener_restore2D_online( img, o, tfun, tfplan, img_rec, msk, shellw, add )
        use simple_ori,       only: ori
        use simple_ctf,       only: ctf
        use simple_projector, only: projector
        class(image),      intent(inout) :: img
        class(ori),        intent(inout) :: o
        class(ctf),        intent(inout) :: tfun
        type(ctfplan),     intent(in)    :: tfplan
        class(image),      intent(inout) :: img_rec
        real,              intent(in)    :: msk
        real,    optional, intent(in)    :: shellw(:)
        logical, optional, intent(in)    :: add
        type(projector) :: proj
        type(image)     :: roimg
        real            :: angast, dfx, dfy, x, y
        logical         :: aadd
        aadd = .true.
        if( present(add) ) aadd = add
        proj = projector('kb')
        ! set CTF parameters
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
        call img%fwd_ft
        ! take care of the nominator
        select case(tfplan%flag)
            case('yes')  ! multiply with CTF
                call tfun%apply(img, dfx, 'ctf', dfy, angast)
            case('flip') ! multiply with abs(CTF)
                call tfun%apply(img, dfx, 'abs', dfy, angast)
            case('mul','no')  ! do nothing
                ! doing nothing
            case DEFAULT
                write(*,*) 'Unsupported ctfflag: ', tfplan%flag
                stop 'simple_filterer :: wiener_restore2D_online'
        end select
        ! shift image
        x = o%get('x')
        y = o%get('y')
        if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH )call img%shift(-x, -y)
        ! griding-based image rotation
        call proj%rotimg(img, -o%e3get(), msk, roimg)
        if( present(shellw) ) call roimg%apply_filter(shellw)
        ! assemble img_rec sum
        if( aadd )then
            call img_rec%add(roimg)
        else
            call img_rec%subtr(roimg)
        endif
        call roimg%kill
    end subroutine wiener_restore2D_online
    
    !>  \brief assembles the ctfsqsums in an online fashion
    subroutine assemble_ctfsqsum_online( img, o, tfun, tfplan, ctfsqsum, shellw, add )
        use simple_ori, only: ori
        use simple_ctf, only: ctf
        class(image),      intent(inout) :: img
        class(ori),        intent(inout) :: o
        class(ctf),        intent(inout) :: tfun
        type(ctfplan),     intent(in)    :: tfplan
        class(image),      intent(inout) :: ctfsqsum
        real,    optional, intent(in)    :: shellw(:)
        logical, optional, intent(in)    :: add
        integer     :: ldim(3)
        type(image) :: ctfsq
        real        :: angast, dfx, dfy, smpd
        logical     :: aadd
        aadd = .true.
        if( present(add) ) aadd = add
        ! set constants
        ldim = img%get_ldim()
        smpd = img%get_smpd()
        ! create & init objs
        call ctfsq%new(ldim, smpd)
        ! set CTF parameters
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
        if( tfplan%flag .ne. 'no' )then
            call tfun%ctf2img(ctfsq, dfx, 'square', dfy, angast)
        else
            ctfsq = cmplx(1.,0.)
        endif
        if( present(shellw) ) call ctfsq%apply_filter(shellw)
        if( aadd )then
            call ctfsqsum%add(ctfsq)
        else
            call ctfsqsum%subtr(ctfsq)
        endif
        call ctfsq%kill
    end subroutine assemble_ctfsqsum_online

end module simple_filterer
