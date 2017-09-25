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
#ifndef USETINY
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
#ifndef USETINY
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
            (self%alpha - 0.5)**2.0 - 0.8)
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
#ifndef USETINY
    r   = self%oneoW * self%bessi0(self%beta * sqrt(arg))
#else
    r   = self%oneoW * bessi0(self%beta * sqrt(arg))
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
#ifndef USETINY
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
  elemental pure real(sp) function bessi0( x )
    ! class(kbinterpol), intent(in) :: self
    real(sp),              intent(in) :: x
    real(dp) :: y, ax
    ax = abs(x)
    if (ax < 3.75) then
       y=x/3.75
       y=y*y
       bessi0=10+y*(3.5156229 + y*(3.0899424 + y*( 1.2067492 +&
                 y*(0.2659732 + y*(0.0360768 + y*0.0045813)))))
    else
       y=3.75/ax
       bessi0=(exp(ax)/sqrt(ax))*(0.39894228 + y*( 0.01328592 +&
            y*( 0.00225319 + y*( -0.00157565 + y*( 0.00916281 +&
            y*(-0.02057706 + y*(  0.02635537 + y*(-0.01647633 +&
            y*  0.00392377))))))))
    end if
  end function bessi0
#endif
  ! subroutine  KBDWindow( window, wsize, alpha)
  !   real(dp), intent(inout) :: window(:)
  !   real(dp), intent(in)    :: wsize, alpha
  !   real(dp)         :: sumvalue
  !   integer          :: i
  !   sumvalue = 0.0
  !   do i=1, INT(wsize/2)
  !      sumvalue = sumvalue +  DblBesselI0(PI * alpha * sqrt(1.0 - 2**(4.0*i/wsize - 1.0)))
  !      window(i) = sumvalue
  !   end do

  !   !! need to add one more value to the nomalization factor at size/2
  !   sumvalue = sumvalue +  DblBesselI0(PI * alpha * sqrt(1.0 - (4.0*(wsize/2)/wsize-1.0)**2))

  !   !! normalize the window and fill in the righthand side of the window:
  !   do i=1, INT(wsize/2)
  !      window(i) = sqrt(window(i)/sumvalue)
  !      window(INT(wsize)-1-i) = window(i)
  !   end do
  ! end subroutine KBDWindow

  !! BesselI0 -- Regular Modified Cylindrical Bessel Function (Bessel I).
  elemental pure real(dp) function DblBesselI0(x)
    real(dp),intent(in) :: x
    real(dp) :: den,num,z
     DblBesselI0=0.
    if (abs(x) < TINY) return
    z = x * x
    num = (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* (z* &
         (z* 0.210580722890567e-22  + 0.380715242345326e-19 ) +&
         0.479440257548300e-16) + 0.435125971262668e-13 ) +&
         0.300931127112960e-10) + 0.160224679395361e-7  ) +&
         0.654858370096785e-5)  + 0.202591084143397e-2  ) +&
         0.463076284721000e0)   + 0.754337328948189e2   ) +&
         0.830792541809429e4)   + 0.571661130563785e6   ) +&
         0.216415572361227e8)   + 0.356644482244025e9   ) +&
         0.144048298227235e10)
    den = (z*(z*(z-0.307646912682801e4) + 0.347626332405882e7) -&
         0.144048298227235e10)
     DblBesselI0= -num/den
  end function DblBesselI0

end module simple_kbinterpol
