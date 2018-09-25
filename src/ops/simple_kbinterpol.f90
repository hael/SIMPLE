! Kaiser-Bessel interpolation kernel
module simple_kbinterpol
use simple_defs
use iso_c_binding
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
    procedure :: apod_dp
    procedure :: dapod
    procedure :: apod_memo
    procedure :: apod_memo_dp
    procedure :: dapod_memo
    procedure :: instr
    procedure :: memoize
end type kbinterpol

interface kbinterpol
   module procedure constructor
end interface kbinterpol

interface
    subroutine kbinterp_memo_set(Whalf_in, alpha_in, Nx_in) bind(c)
        use iso_c_binding
        real(c_double), value :: Whalf_in, alpha_in
        integer(c_int), value :: Nx_in
    end subroutine kbinterp_memo_set
    subroutine kbinterp_memo_memoize() bind(c)
    end subroutine kbinterp_memo_memoize
    subroutine kbinterp_memo_kill() bind(c)
    end subroutine kbinterp_memo_kill
    function apod_nointerp(x) result(r) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        real(c_double)        :: r
    end function apod_nointerp
    function apod_lininterp(x) result(r) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        real(c_double)        :: r
    end function apod_lininterp
    function dapod_nointerp(x) result(r) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        real(c_double)        :: r
    end function dapod_nointerp
    function dapod_lininterp(x) result(r) bind(c)
        use iso_c_binding
        real(c_double), value :: x
        real(c_double)        :: r
    end function dapod_lininterp
    function apod_nointerp_sp(x) result(r) bind(c)
        use iso_c_binding
        real(c_float), value :: x
        real(c_float)        :: r
    end function apod_nointerp_sp
    function apod_lininterp_sp(x) result(r) bind(c)
        use iso_c_binding
        real(c_float), value :: x
        real(c_float)        :: r
    end function apod_lininterp_sp
    function dapod_nointerp_sp(x) result(r) bind(c)
        use iso_c_binding
        real(c_float), value :: x
        real(c_float)        :: r
    end function dapod_nointerp_sp
    function dapod_lininterp_sp(x) result(r) bind(c)
        use iso_c_binding
        real(c_float), value :: x
        real(c_float)        :: r
    end function dapod_lininterp_sp
end interface

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
        self%threshInstr = self%beta/(self%piW) - TINY**2
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
            return ! for insignificant values return as soon as possible
        endif
        arg = self%twooW * x
        arg = 1. - arg * arg
        r = self%oneoW * bessi0(self%beta * sqrt(arg))
    end function apod

    function apod_memo( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real                          :: r
        !r = apod_nointerp_sp(x)
        r = apod_lininterp_sp(x)
    end function apod_memo

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    pure function apod_dp( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real(dp),          intent(in) :: x
        real(dp) :: r, arg
        if( abs(x) > self%Whalf )then
            r = 0.
            return ! for insignificant values return as soon as possible
        endif
        arg = self%twooW * x
        arg = 1._dp - arg * arg
        r = self%oneoW * bessi0_dp(self%beta * sqrt(arg))
    end function apod_dp

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    function apod_memo_dp( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real(dp),          intent(in) :: x
        real(dp)                      :: r
        !r = apod_nointerp(x)
        r = apod_lininterp(x)
    end function apod_memo_dp

    !>  \brief  is the derivative of the Kaiser-Bessel apodization function, abs(x) <= Whalf
    pure function dapod( self, x ) result(r)
        class(kbinterpol), intent(in) :: self
        real(dp),          intent(in) :: x
        real(dp) :: r, arg, sqrtarg
        arg  = self%twooW * x
        arg  = 1._dp - arg * arg
        if (arg <= 0._dp) then
            r = 0._dp
            return
        end if
        sqrtarg = sqrt(arg)
        r    = - 4._dp * self%beta * x * bessi1(self%beta * sqrtarg) / &
                sqrtarg / self%W**3
    end function dapod

    !>  \brief  is the derivative of the Kaiser-Bessel apodization function, abs(x) <= Whalf
    function dapod_memo( self, x ) result(r)
        class(kbinterpol), intent(in) :: self
        real(dp),          intent(in) :: x
        real(dp)                      :: r
        !r = dapod_nointerp(x)
        r = dapod_lininterp(x)
    end function dapod_memo

    subroutine memoize( self, N_in )
        class(kbinterpol), intent(in) :: self
        integer,           intent(in) :: N_in
        call kbinterp_memo_set( real(self%Whalf, kind=dp), real(self%alpha, kind=dp), N_in )
        call kbinterp_memo_memoize()
    end subroutine memoize

    !>  \brief  is the Kaiser-Bessel instrument function
    elemental function instr( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r, arg2
        if ( abs(x) < self%threshInstr)then
            arg2 = sqrt(self%betasq - (self%piW * x)**2)
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

    !>  \brief returns the modified Bessel function I0(x) for any real x
    !! p.378 Handbook of Mathematical Functions, Abramowitz and Stegun
    elemental pure real(dp) function bessi0_dp( x )
        real(dp), intent(in) :: x
        real(dp) :: y, ax
        ax = x ! abs(x)  !! Assumption 1:  beta * sqrt(arg) is always positive
        if ( ax < 3.75d0 ) then
            y= x / 3.75d0
            y=y*y
#ifdef USE_FMA
            bessi0_dp= fma(y,fma(y,fma(y,fma(y,fma(y,fma(y,0.0045813d0,&
                 0.0360768d0),0.2659732d0),1.2067492d0),3.0899424d0),&
                 3.5156229d0),1.0d0)
#else
            bessi0_dp=1.0d0+&
                y*(3.5156229d0 + y*(3.0899424d0 + y*(1.2067492d0 +&
                y*(0.2659732d0 + y*(0.0360768d0 + y* 0.0045813d0)))))
#endif
                   else
            y=3.75d0/ax
            bessi0_dp=( 0.39894228d0 + y*(  0.01328592d0 +&
                y*( 0.00225319d0 + y*( -0.00157565d0 + y*( 0.00916281d0 +&
                y*(-0.02057706d0 + y*(  0.02635537d0 + y*(-0.01647633d0 +&
                y*  0.00392377d0)))))))) * exp( ax ) / sqrt( ax )
        end if
    end function bessi0_dp


    elemental pure real(dp) function bessi1( x )
        real(dp), intent(in) :: x
        real(dp) :: y, ax, bx
        ax = x
        if ( ax < 3.75d0) then
            y = x / 3.75d0
            y = y*y
#ifdef USE_FMA
            bessi1 = x*( &
                fma(y,fma(y,fma(y,fma(y,fma(y,fma(y,0.32411d-3, &
                0.301532d-2),0.2658733d-1),0.15084934d0), &
                0.51498869d0),0.87890594d0),0.5d0) )
#else
            bessi1 = x*(                                             &
                0.5d0 + y*(0.87890594d0 + y*(0.51498869d0 +          &
                y*(0.15084934d0 + y*(0.2658733d-1 + y*(0.301532d-2 + &
                y* 0.32411d-3))))))
#endif
        else
            y   = 3.75d0 / ax
            bx  = exp(ax) / sqrt(ax)
            ax  = 0.39894228d0   + y*(-0.3988024d-1 + y*(-0.362018d-2 + &
                y*(0.163801d-2   + y*(-0.1031555d-1 + y*(0.2282967d-1 + &
                y*(-0.2895312d-1 + y*( 0.1787654d-1 + y*(-0.420059d-2))))))))
            bessi1 = ax*bx
        end if
    end function bessi1


    elemental pure real(dp) function bessi0f( x )
        real(dp), intent(in) :: x
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
        real(dp) :: y,x,xsq
        x=xin
        if (x > 0.5) then
            y = (exp(x) - exp(-x)) / (2. * xin)
        else
            xsq = x * x
            y = (((P3*xsq+P2)*xsq+P1)*xsq + P0)
            y = y / (((xsq+Q2)*xsq+Q1)*xsq + Q0)
        end if
    end function sinhc

end module simple_kbinterpol
