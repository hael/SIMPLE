!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
module simple_kbinterpol_elemental
use simple_defs
implicit none

public :: init_kbiterpol, kb_apod, kb_instr, get_kb_winsz, get_kb_alpha
private

! parameters for the bessel function evaluaation
double precision, parameter :: ps(7) = [1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2]
double precision, parameter :: qs(9) = [0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,0.916281d-2,&
                                      &-0.2057706d-1,0.2635537d-1,-0.1647633d-1,0.392377d-2]
double precision, parameter :: thresh = 3.75d0

! module variables
real :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf

contains

    ! INITIALIZER

    !>  \brief  initializer (should be in params or builder)
    subroutine init_kbiterpol( Whalf_in, alpha_in )
        real, intent(in) :: Whalf_in, alpha_in
        Whalf  = Whalf_in
        alpha  = alpha_in
        W      = 2.0 * Whalf
        piW    = pi * W
        if( Whalf <= 1.5 )then
            beta = 7.4
        else
            beta=pi * sqrt((W**2.0 / alpha**2.0) * (alpha - 0.5)**2.0 - 0.8)
        endif
        betasq = beta * beta
        twooW  = 2.0 / W
        oneoW  = 1.0 / W 
    end subroutine init_kbiterpol

    ! GETTERS

    real function get_kb_winsz()
        get_kb_winsz = Whalf
    end function get_kb_winsz

    real function get_kb_alpha()
        get_kb_alpha = alpha
    end function get_kb_alpha

    ! ELEMENTAL APODIZATION AND INSTRUMENT FUNCTIONS

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    elemental function kb_apod( x ) result( r )
        real, intent(in) :: x
        real :: r, arg
        if( abs(x) > Whalf )then
            r = 0.
            return
        endif
        arg = twooW * x
        arg = 1.-arg * arg
        r   = oneoW * bessi0(beta * sqrt(arg))
    end function kb_apod

    !>  \brief  is the Kaiser-Bessel instrument function
    elemental function kb_instr( x ) result( r )
        real, intent(in) :: x
        real :: r, arg1, arg2
        arg1 = piW * x
        arg1 = betasq - arg1 * arg1
        if( arg1 > 0. )then
            arg2 = sqrt(arg1)
            if( abs(arg2) <= TINY ) then
                r = 1.
            else
                r = sinh(arg2) / (arg2)
            endif
        else
            r = 1
        endif
    end function kb_instr

    !>  \brief returns the modified Bessel function I0(x) for any real x
    elemental function bessi0(x) result(bess)
        real, intent(in) :: x
        real             :: bess
        double precision :: y, ax ! accumulate polynomials in double precision
        if( abs(x) .lt. thresh)then
            y = x/thresh
            y = y*y
            bess = real(ps(1)+y*(ps(3)+y*(ps(4)+y*(ps(5)+y*(ps(6)+y*ps(7))))))
        else
            ax = dble(abs(x))
            y = thresh/ax
            bess = real(exp(ax)/sqrt(ax))*(qs(1)+y*(qs(2)+y*(qs(3)+&
                &y*(qs(4)+y*(qs(5)+y*(qs(6)+y*(qs(7)+y*(qs(8)+y*qs(9)))))))))
        endif
    end function bessi0
    
end module simple_kbinterpol_elemental
