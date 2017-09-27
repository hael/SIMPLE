! Kaiser-Bessel interpolation kernel
module simple_kbinterpol
use simple_defs
implicit none

public :: kbinterpol
private

type :: kbinterpol
    private
#ifndef USETINY
    double precision   :: ps(7) = 0.d0
    double precision   :: qs(9) = 0.d0
    double precision   :: thresh = 0.d0
#endif
    real :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf
 contains
    procedure          :: new
    procedure          :: get_winsz
    procedure          :: get_alpha
    procedure          :: apod
    procedure          :: instr
#if 0
    procedure, private :: bessi0
#endif
end type kbinterpol

interface kbinterpol
   module procedure constructor
end interface kbinterpol

contains

    function constructor( Whalf_in, alpha_in ) result( self )
        real, intent(in) :: Whalf_in, alpha_in
        type(kbinterpol) :: self
        call self%new(Whalf_in, alpha_in)
    end function constructor

    subroutine new( self, Whalf_in, alpha_in )
        class(kbinterpol), intent(inout) :: self
        real,              intent(in)    :: Whalf_in, alpha_in
#if 0
        self%ps = [1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2]
        self%qs = [0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,0.916281d-2,&
            &-0.2057706d-1,0.2635537d-1,-0.1647633d-1,0.392377d-2]
        self%thresh = 3.75d0
#endif
        self%Whalf  = Whalf_in
        self%alpha  = alpha_in
        self%W      = 2.0 * self%Whalf
        self%piW    = pi * self%W
        if( self%Whalf <= 1.5 )then
            self%beta = 7.4
        else
            self%beta = pi * sqrt((self%W**2.0 / self%alpha**2.0) * &
                &(self%alpha - 0.5)**2.0 - 0.8)
        endif
        self%betasq = self%beta * self%beta
        self%twooW  = 2.0 / self%W
        self%oneoW  = 1.0 / self%W
    end subroutine new

    real function get_winsz( self )
        class(kbinterpol), intent(in) :: self
        get_winsz = self%Whalf
    end function get_winsz

    real function get_alpha( self )
        class(kbinterpol), intent(in) :: self
        get_alpha = self%alpha
    end function get_alpha

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    function apod( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r, arg
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        arg = self%twooW * x
        arg = 1. - arg * arg
#if 0
        r = self%oneoW * self%bessi0(self%beta * sqrt(arg))
#else
        r = self%oneoW * bessi0(self%beta * sqrt(arg))
#endif
    end function apod

    !>  \brief  is the Kaiser-Bessel instrument function
    function instr( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r, arg1, arg2
        arg1 = self%piW * x
        arg1 = self%betasq - arg1 * arg1
        if( arg1 > 0. )then
            arg2 = sqrt(arg1)
            if( abs(arg2) <= TINY ) then
                r = 1.0
            else
                r = sinh(arg2) / (arg2)
            endif
        else
            r = 1.0
        endif
    end function instr

#if 0

    !>  \brief returns the modified Bessel function I0(x) for any real x
    function bessi0( self, x ) result( bess )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        double precision :: y, ax !< accumulate polynomials in double precision
        real             :: bess
        if( abs(x) .lt. self%thresh)then
            y = x/self%thresh
            y = y*y
            bess = real(self%ps(1)+y*(self%ps(3)+y*(self%ps(4)+y*(self%ps(5)+y*(self%ps(6)+y*self%ps(7))))))
        else
            ax = dble(abs(x))
            y = self%thresh/ax
            bess = real(exp(ax)/sqrt(ax),sp)*real(self%qs(1)+y*(self%qs(2)+y*(self%qs(3)+&
                &y*(self%qs(4)+y*(self%qs(5)+y*(self%qs(6)+y*(self%qs(7)+y*(self%qs(8)+y*self%qs(9)))))))),sp)
        endif
    end function bessi0

#else

    !>  \brief returns the modified Bessel function I0(x) for any real x
    !! p.378 Handbook of Mathematical Functions, Abramowitz and Stegun
    elemental pure real(dp) function bessi0( x )
        real(sp), intent(in) :: x
        real(dp) :: y, ax
        ax = abs(x)
        if (ax < 3.75d0) then
            y= real(x,dp) / 3.75d0
            y=y*y
            bessi0=1.0d0+&
                y*(3.5156229d0 + y*(3.0899424d0 + y*(1.2067492d0 +&
                y*(0.2659732d0 + y*(0.0360768d0 + y* 0.0045813d0)))))        
        else
            y=3.75d0/ax
            bessi0=( 0.39894228d0 + y*(  0.01328592d0 +&
                y*( 0.00225319d0 + y*( -0.00157565d0 + y*( 0.00916281d0 +&
                y*(-0.02057706d0 + y*(  0.02635537d0 + y*(-0.01647633d0 +&
                y*  0.00392377d0)))))))) * exp( ax ) / sqrt( ax )
        end if
    end function bessi0

#endif

end module simple_kbinterpol
