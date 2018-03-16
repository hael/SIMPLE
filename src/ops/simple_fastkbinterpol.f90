! Kaiser-Bessel interpolation kernel
module simple_fastkbinterpol
use simple_defs
implicit none

public :: kbinterpol
private

type :: kbinterpol
   private
   real :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf, threshInstr
 contains
    procedure :: new
    procedure :: get_winsz
    procedure :: get_alpha
    procedure :: get_wdim
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
        self%threshInstr = self%betasq/(self%piW ** 2) - TINY**2
    end subroutine new

    pure real function get_winsz( self )
        class(kbinterpol), intent(in) :: self
        get_winsz = self%Whalf
    end function get_winsz

    pure real function get_alpha( self )
        class(kbinterpol), intent(in) :: self
        get_alpha = self%alpha
    end function get_alpha

    pure integer function get_wdim( self )
        use simple_math, only: sqwin_1d
        class(kbinterpol), intent(in) :: self
        integer :: win(2)
        call sqwin_1d(0., self%Whalf, win)
        get_wdim = win(2) - win(1) + 1
    end function get_wdim

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    pure function apod( self, x ) result( r )
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
    elemental function instr( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r, arg1, arg2
        ! arg1 = self%piW * x
        ! arg1 = self%betasq - arg1 * arg1
        ! if( arg1 > 0. )then
        !     arg2 = sqrt(arg1)
        !     ! if( abs(arg2) <= TINY ) then !! arg2 is already positive
        !     if(arg2 < TINY)then
        !         r = 1.0
        !     else
        !         r = sinhfme(arg2) / (arg2)
        !     endif
        ! else
        !     r = 1.0
        ! endif
        if ( abs(x) < self%threshInstr)then
            arg2 = sqrt(self%betasq - (self%piW * x)**2)
            ! if( abs(arg2) <= TINY ) then !! arg2 is already positive
            if(arg2 < TINY) then
                r = 1.0
            else
                r = sinhc( arg2 )
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
        ax = x ! abs(x)  !! Assumption 1:  beta * sqrt(arg) is always positive
        if ( ax < 3.75d0 ) then
            y= x / 3.75d0
            y=y*y
#ifdef USE_FMA
            bessi0= fma(y,fma(y,fma(y,fma(y,fma(y,fma(y,0.0045813d0,&
                 0.0360768d0),0.2659732d0),1.2067492d0),3.0899424d0),&
                 3.5156229d0),1.0d0)
#else
            bessi0=1.0d0+&
                y*(3.5156229d0 + y*(3.0899424d0 + y*(1.2067492d0 +&
                y*(0.2659732d0 + y*(0.0360768d0 + y* 0.0045813d0)))))
#endif
                   else
            y=3.75d0/ax
            bessi0=( 0.39894228d0 + y*(  0.01328592d0 +&
                y*( 0.00225319d0 + y*( -0.00157565d0 + y*( 0.00916281d0 +&
                y*(-0.02057706d0 + y*(  0.02635537d0 + y*(-0.01647633d0 +&
                y*  0.00392377d0)))))))) * exp( ax ) / sqrt( ax )
        end if
    end function bessi0

    elemental pure real(dp) function bessi0f( x )
        real(sp), intent(in) :: x
        real(dp) :: y, ax
        ax = abs(x)
        if (ax < 3.75) then
            y= x / 3.75
            y=y*y
            bessi0f=1.0+&
                y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +&
                y*(0.2659732 + y*(0.0360768 + y* 0.0045813)))))
        else
            y=3.75/ax
           ! bessi0f=( 0.39894228 + y*(  0.01328592 + y*( 0.00225319 + y*( -0.00157565 + y*( 0.00916281 +&
           !     y*(-0.02057706 + y*(  0.02635537 + y*(-0.01647633 + y*  0.00392377)))))))) * exp( ax ) / sqrt( ax )

            bessi0f=(exp(ax) *(y *(y *(y *(y *(y *(y *((0.00202623 *y - 0.00850834)* y + 0.0136099) &
                - 0.0106259) + 0.00473165) - 0.000813662) + 0.00116354) + 0.00686082) + 0.206013))/sqrt(1/y)
        end if
    end function bessi0f

    elemental pure function sinhc(xin) result (y)
        !! The coefficients are #2029 from Hart & Cheney. (20.36D)
        real(sp), intent(in) :: xin
        real(dp), parameter  :: P0 = -0.6307673640497716991184787251d+6,&
            P1 = -0.8991272022039509355398013511d+05, &
            P2 = -0.2894211355989563807284660366d+04, &
            P3 = -0.2630563213397497062819489000d+02, &
            Q0 = -0.6307673640497716991212077277d+06, &
            Q1 =  0.1521517378790019070696485176d+05, &
            Q2 = -0.1736789535582336995334509110d+03
       ! logical ::sign
        real(dp) :: y,x,xsq
        x=xin
       ! sign = .false. !! Assumption 1:  input is always positive
        ! if (x < 0 ) then
        !     x = -x
        !     sign = .true.
        ! end if
        ! if (x > 21) then  !! Assumption 2:  input range is less than 12
        !    y = exp(x) / 2
        ! else
        if (x > 0.5) then
            y = (exp(x) - exp(-x)) / 2
        else
            xsq = x * x
            y = (((P3*xsq+P2)*xsq+P1)*xsq + P0)
            y = y / (((xsq+Q2)*xsq+Q1)*xsq + Q0)
        end if

        ! if (sign) then
        !     y = -y
        ! end if

    end function sinhc
end module simple_fastkbinterpol
