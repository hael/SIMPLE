! Kaiser-Bessel interpolation kernel
module simple_kbinterpol
use simple_defs
implicit none

public :: kbinterpol
private

type :: kbinterpol
   private
   real :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf
 contains
    procedure :: new
    procedure :: get_winsz
    procedure :: get_alpha
    procedure :: apod
    procedure :: instr
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
        r = self%oneoW * bessi0(self%beta * sqrt(arg))
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

end module simple_kbinterpol
