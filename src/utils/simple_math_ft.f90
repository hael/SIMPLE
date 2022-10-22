module simple_math_ft
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
        calc_fourier_index = nint((real(box)*smpd)/res)
    end function calc_fourier_index

    !> \brief calculate logical mask filtering out the Graphene bands
    function calc_graphene_mask( box, smpd ) result( mask )
        integer, intent(in)  :: box
        real,    intent(in)  :: smpd
        real,    allocatable :: res(:), sqdiff_band1(:), sqdiff_band2(:)
        logical, allocatable :: mask(:)
        integer, parameter   :: NBANDS = 3
        integer              :: loc(NBANDS), n, i
        res = get_resarr( box, smpd )
        n   = size(res)
        allocate(sqdiff_band1(n), source=(res - GRAPHENE_BAND1)**2.0)
        allocate(sqdiff_band2(n), source=(res - GRAPHENE_BAND2)**2.0)
        allocate(mask(n), source=.true.)
        loc = minnloc(sqdiff_band1, NBANDS)
        do i=1,NBANDS
            mask(loc(i)) = .false.
        end do
        loc = minnloc(sqdiff_band2, NBANDS)
        do i=1,NBANDS
            mask(loc(i)) = .false.
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

    pure subroutine get_find_at_crit( n, corrs, crit, find )
        integer, intent(in)  :: n
        real,    intent(in)  :: corrs(n), crit
        integer, intent(out) :: find
        find = 1
        do while( find <= n )
            if( corrs(find) >= crit )then
                find = find + 1
                cycle
            else
                exit
            endif
        end do
    end subroutine get_find_at_crit

    !>   returns the Fourier index of the resolution limit at corr
    integer function get_lplim_at_corr( fsc, corr, incrreslim )
        real,    intent(in)           :: fsc(:), corr
        logical, intent(in), optional :: incrreslim
        integer :: n, h
        n = size(fsc)
        if( n < 3 )then
            call simple_exception('nonconforming size of fsc array; get_lplim_at_corr', __FILENAME__ , __LINE__)
        endif
        get_lplim_at_corr = n-1
        if( present(incrreslim) )then
            if( incrreslim )then
                do h = n-1,3,-1
                    if( fsc(h) > corr )then
                        get_lplim_at_corr = h
                        exit
                    endif
                enddo
                get_lplim_at_corr = min(get_lplim_at_corr+10,n-1)
            else
                do h=3,n-1
                    if( fsc(h) >= corr )then
                        cycle
                    else
                        get_lplim_at_corr = h - 1
                        exit
                    endif
                end do
            endif
        else
            do h=3,n-1
                if( fsc(h) >= corr )then
                    cycle
                else
                    get_lplim_at_corr = h - 1
                    exit
                endif
            end do
        endif
    end function get_lplim_at_corr

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

end module simple_math_ft
