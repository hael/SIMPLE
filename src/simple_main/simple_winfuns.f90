!==Class simple_winfuns
!
! the simple_winfuns provides interfaces and definitions for the window/instrument functions used in the SIMPLE library. 
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution 
! or modification is regulated by the GNU General Public License. *Author:* Hans Elmlund, 2013-08-12.
! 
!==Changes are documented below
!
! RELEVANT INFO:
! http://mathworld.wolfram.com/ApodizationFunction.html
!
! Jackson 1991, Least Relative Aliased Energy on a 1X grid
!           Gaussian Kaiser-Bessel
!             var       beta
! W=3       0.4241    1.9980     
! W=4       0.4927    2.3934
! W=5       0.4839    3.3800

! Beatty 2005, 2X grid
!           Jackson   Beatty
!             beta     beta
! W=3        6.6875   6.4861     
! W=4        9.1375   8.9962
! W=5       11.5250  11.4410
!
module simple_winfuns
use simple_defs
implicit none

public :: winfuns
private

type :: winfuns
    private
    character(STDLEN) :: wfun_str=''                      !< wfun string descriptor
    procedure(wfun), pointer, nopass :: apod_fun=>null()  !< apodization function
    procedure(ifun), pointer, nopass :: instr_fun=>null() !< instrument function
    real :: Whalf=0.                                      !< window halfwidth
    real :: alpha=0.                                      !< oversampling (padding) fatcor
    real :: W=0.                                          !< window fullwidth      
    real :: beta=0.                                       !< KB shape factor
    real :: betasq, twooW, twoW, oneoW, piW, pioW
    real :: bmanc1, bmanc2, bmanc3, bmanc4, bmanc5, bmanc6, bmanc7
  contains
    procedure :: which
    procedure :: get_Whalf
    procedure :: eval_apod
    procedure :: eval_instr
    procedure, private :: bman_apod
    procedure, private :: bman_instr
    procedure, private :: hann_apod
    procedure, private :: hann_instr
    procedure, private :: kb_apod
    procedure, private :: kb_instr
    procedure, private :: lin_apod
    procedure, private :: lin_instr
    procedure, private :: nn_apod
    procedure, private :: nn_instr
    procedure, private :: sinc_apod
    procedure, private :: sinc_instr
end type

!>  \brief  defines the interface for the window
abstract interface
    function wfun( self, x ) result( w )
        import winfuns
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: w
    end function
end interface

!>  \brief  defines the interface for the instrument function
abstract interface
    function ifun( self, x ) result( w )
        import winfuns
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: w
    end function
end interface

interface winfuns
    module procedure constructor
end interface

! parameters for the bessel function evaluaation
double precision, parameter :: ps(7) = [1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2]
double precision, parameter :: qs(9) = [0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,0.916281d-2,&
                                      &-0.2057706d-1,0.2635537d-1,-0.1647633d-1,0.392377d-2]
double precision, parameter :: thresh = 3.75d0

contains

    !> \brief  is a constructor
    function constructor( wfun_str, Whalf, alpha ) result( self )
        character(len=*), intent(in) :: wfun_str
        real, intent(in)             :: Whalf, alpha
        type(winfuns)                :: self
        ! set window function constants
        self%wfun_str = wfun_str
        self%Whalf = Whalf
        self%alpha = alpha
        self%W=2.*self%Whalf
        if( self%Whalf <= 1.5 )then
            self%beta = 7.4
        else
            self%beta=pi*sqrt((self%W**2./self%alpha**2.)*(self%alpha-0.5)**2.-0.8)
        endif
        self%betasq = self%beta*self%beta
        self%twooW = 2./self%W 
        self%twoW = 2.*self%W 
        self%oneoW = 1./self%W 
        self%piW = pi*self%W
        self%pioW = pi/self%W
        self%bmanc1=21./50.            
        self%bmanc2=2./25.
        self%bmanc3=pi/self%Whalf
        self%bmanc4=twopi/self%Whalf
        self%bmanc5=21./25.
        self%bmanc6=9./25.
        self%bmanc7=twopi*self%Whalf
        ! set window function pointers
        select case(self%wfun_str)
            case('')
                self%apod_fun  => kb_apod
                self%instr_fun => kb_instr
            case('bman')
                self%apod_fun  => bman_apod
                self%instr_fun => bman_instr
            case('hann')
                self%apod_fun  => hann_apod
                self%instr_fun => hann_instr
            case('kb')
                self%apod_fun  => kb_apod
                self%instr_fun => kb_instr
            case('lin')
                self%apod_fun  => lin_apod
                self%instr_fun => lin_instr
            case('nn')
                self%apod_fun  => nn_apod
                self%instr_fun => nn_instr
            case('sinc')
                self%apod_fun  => sinc_apod
                self%instr_fun => sinc_instr
            case DEFAULT
                write(*,*) 'window function:', trim(self%wfun_str), 'Unsupprted constructor; simple_winfuns'
                stop
        end select
    end function
    
    !> \brief  2 check which window function
    function which( self ) result( this )
        class(winfuns), intent(in) :: self
        character(len=STDLEN) :: this
        this = self%wfun_str
    end function
    
    !> \brief  2 get window half-width
    function get_Whalf( self ) result( Whalf )
        class(winfuns), intent(in) :: self
        real :: Whalf
        Whalf = self%Whalf
    end function
    
    !> \brief  2 evaluate the apodization function
    function eval_apod( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: r
        r = self%apod_fun(self,x)
    end function
    
    !> \brief  2 evaluate the instrument function
    function eval_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: r
        r = self%instr_fun(self,x)
    end function
    
    ! APODIZATION/INSTRUMENT FUNCTION PAIRS
    
    ! BLACKMANN
    
    !>  \brief  is the Blackmann apodization function
    function bman_apod( self, x ) result( r ) ! OK
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        r = self%bmanc1+0.5*cos(self%bmanc3*x)+self%bmanc2*cos(self%bmanc4*x)
    end function
    
    !>  \brief  is the Blackmann instrument function
    function bman_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r, sincarg, arg, div
        sincarg = self%bmanc7*x
        if( abs(sincarg) < 0.00000001 ) then
            r = 1.
        else
            r = sin(sincarg)/(sincarg)
        endif
        arg = self%Whalf*x
        arg = arg*arg
        div = (1.-arg)*(1.-4.*arg)
        if( abs(div) < 0.00000001 ) then
            r = 1.
        else
            r = (self%Whalf*(self%bmanc5-self%bmanc6*arg)*r)/div
        endif
    end function
    
    ! HANNING
    
    !>  \brief  is the Hanning apodization function
    function hann_apod( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        r = cos(self%pioW*x)
        r = r*r
    end function
    
    !>  \brief  is the Hanning instrument function
    function hann_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r, arg1, arg2, div
        arg1 = self%piW*x
        if( abs(arg1) < 0.00000001 ) then
            r = 1.
        else
            r = sin(arg1)/(arg1)
        endif
        arg2 = 2.*self%Whalf*x
        arg2 = arg2*arg2
        div  = 1.-arg2
        if( abs(div) < 0.00000001 ) then
            r = 1.
        else
            r = (self%Whalf*r)/div
        endif
    end function
    
    ! KAISER-BESSEL
    
    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    function kb_apod( self, x ) result( r ) ! OK
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r, arg
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        arg = self%twooW*x
        arg = 1.-arg*arg
        r = self%oneoW*bessi0(self%beta*sqrt(arg))
    end function
    
    !>  \brief  is the Kaiser-Bessel instrument function
    function kb_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r, arg1, arg2
        arg1 = self%piW*x
        arg1 = self%betasq-arg1*arg1
        if( arg1 > 0. )then
            arg2 = sqrt(arg1)
            if( abs(arg2) < 0.00000001 ) then
                r = 1.
            else
                r = sinh(arg2)/(arg2)
            endif
        else
            r = 1
        endif
    end function
    
    ! LINEAR
    
    !>  \brief  returns the normalized linear apodization function of x
    !!          if x = 0., linear = 1. if x .ge. 1.0, linear = 0.
    !!          also called Bartlett window
    function lin_apod( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) ::  x
        real :: ax, r
        ax = abs(x)
        if( ax .gt. self%Whalf )then
            r = 0.
        else
            r = 1.-ax/self%Whalf
        endif
    end function 

    !>  \brief  the instrument function to linear interpolation kernel
    function lin_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) ::  x
        real :: r, arg
        arg = self%sinc_apod(x*self%Whalf)
        r = self%Whalf*arg*arg
    end function
    
    ! NEAREST NEIGHBOR
    
    !>  \brief  nearest neighbor apodization function, returns 0.0 if abs(x) .gt. 0.5, else 1.0
    !!          warning: make sure that x is positive!
    function nn_apod( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: r
        if( abs(x) .gt. self%Whalf )then
            r = 0
        else
            r = 1.
        endif
    end function

    !>  \brief  instrument fun of nearest neighbor kernel
    function nn_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real :: r
        r = self%sinc_apod(x*self%Whalf)
    end function
    
    ! SINC

    !>  \brief  is a normalized sinc apodization function
    function sinc_apod( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real, intent(in) :: x
        real             :: r, arg
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        if( abs(x) < 0.00000001 ) then
            r = 1.
        else
            arg = pi*x
            r = sin(arg)/(arg)
        endif
    end function
    
    !>  \brief  returns the instrument function to sinc, which is a rectangle
    function sinc_instr( self, x ) result( r )
        class(winfuns), intent(in) :: self
        real,           intent(in) :: x
        real :: r
        if( abs(x) .gt. self%Whalf )then
            r = 0
        else
            r = 1.
        endif
    end function
    
    ! BESSEL FUNCTION FOR THE KAISER-BESSEL WINDOW

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
    
    !>  \brief returns the modified Bessel function I0(x) for any real x
    ! function bessi0(x) result(bess)
    !     real, intent(in) :: x
    !     real :: bess, ax
    !     double precision :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y ! accumulate polynomials in double precision
    !     save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    !     data p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
    !     data q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,0.225319d-2,-0.157565d-2,0.916281d-2,&
    !     -0.2057706d-1,0.2635537d-1,-0.1647633d-1,0.392377d-2/
    !     if( abs(x) .lt. 3.75 )then
    !         y = x/3.75
    !         y = y*y
    !         bess = p1+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))
    !     else
    !         ax = abs(x)
    !         y = 3.75/ax
    !         bess = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
    !     endif
    ! end function
    
    !>  \brief returns the modified Bessel function I0(x) for positive real x
    ! function bessk0(x) result(bess)
    !     real, intent(in) :: x
    !     real :: bess
    !     double precision :: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,y ! accumulate polynomials in double precision
    !     save p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
    !     data p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
    !     data q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,-0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
    !     if( x .le. 2.0 )then ! polynomial fit
    !         y = x*x/4.0
    !         bess = (-log(x/2.0)*bessi0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    !     else
    !         y = (2.0/x)
    !         bess = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7)))))))
    !     endif
    ! end function
  
end module simple_winfuns
