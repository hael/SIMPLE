!@descr: Fourier transform-related math routines
module simple_math_ft
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs
use simple_error, only: simple_exception
use simple_is_check_assert
use simple_srch_sort_loc
implicit none

interface csq_fast
    module procedure csq_fast_1, csq_fast_2
end interface

interface csq
    module procedure csq_1, csq_2
end interface

contains

    !>   returns the Fourier index of resolution
    integer pure function calc_fourier_index( res, box, smpd )
        real, intent(in)    :: res, smpd
        integer, intent(in) :: box
        real :: res_here
        res_here = max(res,2.*smpd)
        calc_fourier_index = min(box/2,nint((real(box)*smpd)/res_here))
    end function calc_fourier_index

    !> \brief logical mask over the resolution shells of get_resarr(box, smpd): .false. at the three
    !! shells nearest each band (a resolution in Angstrom, e.g. the graphene bands) and .true.
    !! elsewhere. A band finer than Nyquist (2*smpd) is not in the spectrum and masks nothing.
    function calc_graphene_mask( box, smpd, bands ) result( mask )
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd, bands(:)
        integer, parameter   :: NSHELLS = 3
        real,    allocatable :: res(:), sqdiff(:)
        logical, allocatable :: mask(:)
        integer              :: loc(NSHELLS), n, i, iband
        res = get_resarr( box, smpd )
        n   = size(res)
        allocate(mask(n), source=.true.)
        allocate(sqdiff(n))
        do iband = 1,size(bands)
            if( bands(iband) < 2.*smpd ) cycle ! beyond Nyquist
            sqdiff = (res - bands(iband))**2
            loc    = minnloc(sqdiff, NSHELLS)
            do i = 1,NSHELLS
                mask(loc(i)) = .false.
            end do
        end do
    end function calc_graphene_mask
    
    !>   returns the Fourier index of res
    real pure function calc_lowpass_lim( find, box, smpd )
        integer, intent(in) :: find, box !< box size
        real, intent(in)    :: smpd      !< smpd pixel size \f$ (\si{\angstrom}) \f$
        calc_lowpass_lim = real(box)*smpd/real(find)
    end function calc_lowpass_lim

    ! faster unsafe complex magnitude, single precision
    real(sp) elemental function csq_fast_1( comp )
        complex(sp), intent(in) :: comp
        csq_fast_1 = real(comp*conjg(comp))
    end function csq_fast_1

    ! faster unsafe complex magnitude, double precision
    real(dp) elemental function csq_fast_2( comp )
        complex(dp), intent(in) :: comp
        csq_fast_2 = real(comp*conjg(comp),dp)
    end function csq_fast_2

    !>   is for complex squaring
    elemental function csq_1( a ) result( sq )
        complex(sp), intent(in) :: a !< complx component
        real(sp) :: sq, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x)) then
            sq = y*y
        else if( is_zero(y) ) then
           sq = x*x
        else if( x > y ) then
            frac = y/x
            sq = x*x*(1.+frac*frac)
        else
            frac = x/y
            sq = y*y*(1.+frac*frac)
        endif
    end function csq_1

    !>  is for double complex squaring
    elemental function csq_2( a ) result( sq )
        complex(dp), intent(in) :: a !< complx component
        real(dp) :: sq, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x)) then
            sq = y*y
        else if( is_zero(y) ) then
            sq = x*x
        else if( x > y ) then
            frac = y/x
            sq = x*x*(1.+frac*frac)
        else
            frac = x/y
            sq = y*y*(1.+frac*frac)
        endif
    end function csq_2

    !>   one-dimensional cyclic index generation
    !! \param lims limits
    pure function cyci_1d( lims, i ) result( ind )
        integer, intent(in) :: lims(2), i !< input var
        integer :: ind, del
        ind = i
        if( ind > lims(2) )then
            del = ind-lims(2)
            ind = lims(1)+del-1
        else if( ind < lims(1) )then
            del = lims(1)-ind
            ind = lims(2)-del+1
        endif
    end function cyci_1d

    !>  as per cyci_1d, assumes lower limit equals 1
    pure function cyci_1d_static( ulim, i ) result( ind )
        integer, intent(in) :: ulim, i !< input vars upper limit and index
        integer :: ind ! return index
        ind = merge(ulim+i, i , i<1)
        ind = merge(i-ulim,ind ,i>ulim)
    end function cyci_1d_static

    !>   is for working out the fourier dimension
    pure function fdim( d ) result( fd )
        integer, intent(in) :: d !< dimension
        integer :: fd
        if(is_even(d))then
            fd = d/2+1
        else
            fd = (d-1)/2+1
        endif
    end function fdim

    !=============================================================
    ! Accessors for half-plane stored fplane (native lattice, k<=0 stored)
    !=============================================================
    pure function fplane_get_cmplx( fpl, h, k ) result(c)
        use, intrinsic :: iso_c_binding, only: c_float_complex
        use simple_type_defs, only: fplane_type
        implicit none
        type(fplane_type), intent(in) :: fpl
        integer,           intent(in) :: h, k
        complex(c_float_complex) :: c
        if (k <= 0) then
            c = fpl%cmplx_plane(h,k)
        else
            c = conjg( fpl%cmplx_plane(-h,-k) )
        end if
    end function fplane_get_cmplx

    !> The plane of a particle at another shift: the numerator re-shifted by dshift (the new shift minus the one
    !! the plane was generated at, in the pixels of the plane's padded box, the cropped ones when cropped) with
    !! the phase gen_fplane4rec applies, exp(i (h p1 + k p2)) at padded sample (h,k) with p = -dshift shconst,
    !! the stored index (hn,kn) holding the padded sample (hn stride, kn stride); the CTF^2 plane is shift-free and
    !! copied, the forward transfer (when stored) carries the same phase. The weighted insertion (topk > 1)
    !! inserts a particle at its extra poses from this plane, without reading or transforming the particle again
    subroutine fplane_shift( fpl, dshift, stride, fpl_out )
        use, intrinsic :: iso_c_binding, only: c_float_complex, c_float
        use simple_type_defs, only: fplane_type
        type(fplane_type), intent(in)    :: fpl
        real,              intent(in)    :: dshift(2)
        integer,           intent(in)    :: stride
        type(fplane_type), intent(inout) :: fpl_out
        complex(c_float_complex) :: w1, w2, ph0, ph_h, ph_k
        real(dp) :: pshift(2)
        integer  :: hlo, hhi, klo, h, k
        if( .not. allocated(fpl%cmplx_plane) ) return
        hlo = lbound(fpl%cmplx_plane,1); hhi = ubound(fpl%cmplx_plane,1)
        klo = lbound(fpl%cmplx_plane,2)
        fpl_out%frlims  = fpl%frlims
        fpl_out%shconst = fpl%shconst
        fpl_out%nyq     = fpl%nyq
        if( allocated(fpl_out%cmplx_plane) )then
            if( any(lbound(fpl_out%cmplx_plane) /= lbound(fpl%cmplx_plane)) .or. &
               &any(ubound(fpl_out%cmplx_plane) /= ubound(fpl%cmplx_plane)) ) deallocate(fpl_out%cmplx_plane)
        endif
        if( .not. allocated(fpl_out%cmplx_plane) ) allocate(fpl_out%cmplx_plane(hlo:hhi,klo:0))
        if( allocated(fpl_out%ctfsq_plane) )then
            if( any(lbound(fpl_out%ctfsq_plane) /= lbound(fpl%ctfsq_plane)) .or. &
               &any(ubound(fpl_out%ctfsq_plane) /= ubound(fpl%ctfsq_plane)) ) deallocate(fpl_out%ctfsq_plane)
        endif
        if( .not. allocated(fpl_out%ctfsq_plane) ) allocate(fpl_out%ctfsq_plane(hlo:hhi,klo:0))
        fpl_out%ctfsq_plane = fpl%ctfsq_plane
        if( allocated(fpl%transfer_plane) )then
            if( allocated(fpl_out%transfer_plane) )then
                if( any(lbound(fpl_out%transfer_plane) /= lbound(fpl%transfer_plane)) .or. &
                   &any(ubound(fpl_out%transfer_plane) /= ubound(fpl%transfer_plane)) ) deallocate(fpl_out%transfer_plane)
            endif
            if( .not. allocated(fpl_out%transfer_plane) ) allocate(fpl_out%transfer_plane(hlo:hhi,klo:0))
        else if( allocated(fpl_out%transfer_plane) )then
            deallocate(fpl_out%transfer_plane)
        endif
        ! the phase recurrence of gen_fplane4rec on the stored lattice
        pshift = real(-dshift * fpl%shconst(1:2), dp)
        w1   = cmplx(real(cos(real(stride,dp)*pshift(1)), c_float), real(sin(real(stride,dp)*pshift(1)), c_float), kind=c_float_complex)
        w2   = cmplx(real(cos(real(stride,dp)*pshift(2)), c_float), real(sin(real(stride,dp)*pshift(2)), c_float), kind=c_float_complex)
        ph0  = cmplx(real(cos(real(hlo*stride,dp)*pshift(1)), c_float), real(sin(real(hlo*stride,dp)*pshift(1)), c_float), kind=c_float_complex)
        ph_k = cmplx(real(cos(real(klo*stride,dp)*pshift(2)), c_float), real(sin(real(klo*stride,dp)*pshift(2)), c_float), kind=c_float_complex)
        do k = klo, 0
            ph_h = ph0
            do h = hlo, hhi
                fpl_out%cmplx_plane(h,k) = fpl%cmplx_plane(h,k) * (ph_k * ph_h)
                if( allocated(fpl%transfer_plane) ) fpl_out%transfer_plane(h,k) = fpl%transfer_plane(h,k) * (ph_k * ph_h)
                ph_h = ph_h * w1
            enddo
            ph_k = ph_k * w2
        enddo
    end subroutine fplane_shift

    pure function fplane_get_ctfsq( fpl, h, k ) result(v)
        use, intrinsic :: iso_c_binding, only: c_float
         use simple_type_defs, only: fplane_type
        implicit none
        type(fplane_type), intent(in) :: fpl
        integer,           intent(in) :: h, k
        real(c_float) :: v
        if (k <= 0) then
            v = fpl%ctfsq_plane(h,k)
        else
            v = fpl%ctfsq_plane(-h,-k)
        end if
    end function fplane_get_ctfsq

    function get_find_at_res( res, crit_res ) result( find )
        real,    intent(in) :: res(:), crit_res
        integer :: n, h, find
        n = size(res)
        find = n - 1
        do h=3,n-1
            if( res(h) >= crit_res )then
                cycle
            else
                find = h - 1
                exit
            endif
        end do
    end function get_find_at_res

    function get_find_at_crit( corrs, crit_corr, incrreslim ) result( find )
        real,    intent(in) :: corrs(:), crit_corr
        logical, intent(in), optional :: incrreslim
        integer :: n, h, find
        logical :: iincrreslim
        iincrreslim = .false.
        if( present(incrreslim) ) iincrreslim = incrreslim
        n = size(corrs)
        find = n - 1
        do h=3,n-1
            if( corrs(h) >= crit_corr )then
                cycle
            else
                find = h - 1
                exit
            endif
        end do
        if( iincrreslim )then
            find = min(find+10,n-1)
        else
            find = min(find,n-1)
        endif
    end function get_find_at_crit

    !>  \brief get array of resolution steps
    function get_resarr( box, smpd ) result( res )
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        real, allocatable   :: res(:)
        integer :: n, k
        n = fdim(box) - 1
        allocate( res(n) )
        do k=1,n
            res(k) = calc_lowpass_lim(k, box, smpd)
        end do
    end function get_resarr

    !>   is for calculating complex arg/abs/modulus, from numerical recipes
    pure elemental function mycabs( a ) result( myabs )
        complex, intent(in) :: a      !< complx component
        real                :: myabs, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x) ) then
            myabs = y
        else if( is_zero(y)  ) then
           myabs = x
        else if( x > y ) then
            frac = y/x
            myabs = x*sqrt(1.+frac*frac)
        else
            frac = x/y
            myabs = y*sqrt(1.+frac*frac)
        endif
    end function mycabs

    !>   is for calculating the phase angle of a Fourier component
    !! \return phase phase angle only meaningful when cabs(comp) is well above the noise level
    elemental function phase_angle( comp ) result( phase )
        complex, intent(in) :: comp !< Fourier component
        real :: nom, denom, phase
        nom = aimag(comp)
        denom = real(comp)
        if( is_zero(denom) )then
           if( is_zero(nom) )then
                phase = 0.
            else if( nom > 0. )then
                phase = pi/2.
            else
                phase = -pi/2.
            endif
            return
        endif
        phase = atan(nom/denom)
    end function phase_angle

    !>   get the resolution in angstrom, given angle and diameter
    !! \param ang,diam angular resolution (degrees) and diameter (\f$\si{\angstrom}\f$)
    pure function resang( ang, diam ) result( res )
        real, intent(in)  :: ang, diam
        real :: res                      !< spatial resolution (\f$\si{\per\angstrom}\f$)
        res = (ang/360.)*(pi*diam)
    end function resang

    !> Map a positive shell-wise noise spectrum onto another Fourier lattice.
    !! input_shell_step gives the input-shell coordinate increment represented
    !! by one output shell. Use 1 for a physical-extent-preserving Fourier crop;
    !! use nyq/nyq_out when padding refines the Fourier-shell spacing. Values
    !! below kstart and above nyq are extended from the nearest valid endpoint.
    subroutine resample_sigma2( kstart, nyq, sigma2, nyq_out, input_shell_step, sigma2_out )
        integer, intent(in)  :: kstart, nyq, nyq_out
        real,    intent(in)  :: sigma2(kstart:nyq), input_shell_step
        real,    intent(out) :: sigma2_out(0:nyq_out)
        real    :: ksc, d
        integer :: k, f, c
        if( kstart < 0 .or. nyq < kstart .or. nyq_out < 0 )then
            call simple_exception('invalid sigma2 shell bounds', __FILENAME__, __LINE__)
        endif
        if( .not. ieee_is_finite(input_shell_step) .or. input_shell_step <= 0.0 )then
            call simple_exception('invalid sigma2 shell mapping', __FILENAME__, __LINE__)
        endif
        do k = kstart, nyq
            if( .not. ieee_is_finite(sigma2(k)) .or. sigma2(k) <= 0.0 )then
                call simple_exception('sigma2 values must be finite and positive', __FILENAME__, __LINE__)
            endif
        enddo
        do k = 0, nyq_out
            ksc = input_shell_step * real(k)
            if( ksc <= real(kstart) )then
                sigma2_out(k) = sigma2(kstart)
            else if( ksc >= real(nyq) )then
                sigma2_out(k) = sigma2(nyq)
            else
                f = floor(ksc)
                c = f + 1
                d = ksc - real(f)
                sigma2_out(k) = (1.0-d)*sigma2(f) + d*sigma2(c)
            endif
        enddo
    end subroutine resample_sigma2

end module simple_math_ft
