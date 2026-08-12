!@descr: Kaiser-Bessel interpolation kernel
!
! meaning of α: the oversampling factor
! | α  | Padding | Meaning                                       |
! | -- | ------- | --------------------------------------------- |
! | 1  | none    | KB tries to suppress aliasing it cannot avoid |
! | 2  | 2×      | KB is tuned for a wide guard band             |
! | >2 | heavy   | diminishing returns                           |
!
! meaning of β: controls how aggressively the KB kernel suppresses high-frequency content in its Fourier transform
! | Parameter | Controls               | What breaks if wrong     |
! | --------- | ---------------------- | ------------------------ |
! | α         | Guard band size        | Unavoidable aliasing     |
! | W         | Cost vs smoothness     | Runtime / memory         |
! | β         | Spectral concentration | Leakage *or* instability |
!
module simple_kbinterpol
use simple_defs
use iso_c_binding
use simple_edges_sqwins, only: sqwin_1d
implicit none

public :: kbinterpol, kb_windim
public :: apod_device, apod_fast_device, apod_kb15_a2
private

type :: kbinterpol
    private
    real :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf, threshInstr
  contains
    procedure          :: new
    procedure, private :: copy
    generic            :: assignment(=) => copy
    ! Getters
    procedure          :: get_winsz
    procedure          :: get_alpha
    procedure          :: get_wdim
    ! Calculators
    procedure          :: apod
    procedure          :: apod_fast
    procedure          :: apod_fast_value_deriv
    procedure          :: apod_mat_2d
    procedure          :: apod_mat_2d_fast
    procedure          :: apod_mat_3d
    procedure          :: apod_mat_3d_fast
    procedure          :: apod_mat_3d_fast_grad
    procedure          :: instr
end type kbinterpol

interface kbinterpol
   module procedure constructor
end interface kbinterpol

#ifdef USE_OPENMP_OFFLOAD
!$omp declare target(apod_device, apod_fast_device, apod_kb15_a2, bessi0)
#endif

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
        if( abs(self%Whalf - KBWINSZ) <= 1.e-6 .and. abs(self%alpha - KBALPHA) <= 1.e-6 )then
            self%beta = KB_BETA_KB15_A2
        else
            self%beta = pi * sqrt((self%W**2.0 / self%alpha**2.0) * &
                &(self%alpha - 0.5)**2.0 - 0.8)
        endif
        self%betasq = self%beta * self%beta
        self%twooW  = 2.0 / self%W
        self%oneoW  = 1.0 / self%W
        self%threshInstr = self%beta/(self%piW) - TINY**2
    end subroutine new

    subroutine copy( self, self2copy )
        class(kbinterpol), intent(inout) :: self
        class(kbinterpol), intent(in)    :: self2copy
        self%alpha       = self2copy%alpha
        self%beta        = self2copy%beta
        self%betasq      = self2copy%betasq
        self%oneoW       = self2copy%oneoW
        self%piW         = self2copy%piW
        self%twooW       = self2copy%twooW
        self%W           = self2copy%W
        self%Whalf       = self2copy%Whalf
        self%threshInstr = self2copy%threshInstr
    end subroutine copy

    ! Getters

    pure real function get_winsz( self )
        class(kbinterpol), intent(in) :: self
        get_winsz = self%Whalf
    end function get_winsz

    pure real function get_alpha( self )
        class(kbinterpol), intent(in) :: self
        get_alpha = self%alpha
    end function get_alpha

    pure integer function get_wdim( self )
        class(kbinterpol), intent(in) :: self
        integer :: win(2)
        call sqwin_1d(0., self%Whalf, win(1), win(2))
        get_wdim = win(2) - win(1) + 1
    end function get_wdim

    ! Calculators

    !>  \brief  is the Kaiser-Bessel apodization function, abs(x) <= Whalf
    pure elemental function apod( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r, arg
        if( abs(x) > self%Whalf )then
            r = 0.
            return ! for insignificant values return as soon as possible
        endif
        arg = self%twooW * x
        arg = 1. - arg * arg
        r = real(self%oneoW * bessi0(self%beta * sqrt(arg)), kind=kind(r))
    end function apod

    !>  \brief  fast KB apodization for the common KBWINSZ=1.5, KBALPHA=2 kernel
    !!
    !! Method used:
    !!   apod(x) = (1/W) * I0(beta * sqrt(u)), u = 1 - (2x/W)^2.
    !! For the configured KBWINSZ=1.5, KBALPHA=2 beta this is evaluated as a
    !! degree-14 polynomial in u using
    !! the modified Bessel I0 power series,
    !!   I0(beta*sqrt(u)) = sum_k ((beta^2/4)^k / (k!)^2) * u^k.
    !! The coefficients below are those terms for k=0..14, evaluated with
    !! Horner's rule in single precision because apod returns real(sp).
    !! This is not a lookup table and performs no linear interpolation.  For
    !! any other KB parameters the routine falls back to the original apod.
    pure elemental function apod_fast( self, x ) result( r )
        class(kbinterpol), intent(in) :: self
        real,              intent(in) :: x
        real :: r
        real(sp), parameter :: coeffs(0:14) = [&
            1.0000000_sp, 1.0517297e1_sp, 2.7653385e1_sp, 3.2315430e1_sp, 2.1241936e1_sp,&
            8.9363103_sp, 2.6107175_sp, 5.6036106e-1_sp, 9.2085685e-2_sp, 1.1956698e-2_sp,&
            1.2575214e-3_sp, 1.0930353e-4_sp, 7.9831782e-6_sp, 4.9681336e-7_sp, 2.6658846e-8_sp]
        real(sp) :: p, q, u
        integer  :: i
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        if( abs(self%Whalf - 1.5) > 1.e-6 .or. abs(self%alpha - 2.0) > 1.e-6 )then
            r = self%apod(x)
            return
        endif
        q = self%twooW * x
        u = max(0._sp, 1._sp - q*q)
        p = coeffs(14)
        do i = 13,0,-1
            p = p * u + coeffs(i)
        enddo
        r = self%oneoW * p
    end function apod_fast

    !>  \brief  Joint value and derivative of the standard fast KB polynomial.
    !!
    !! The derivative is with respect to the signed kernel argument x. Value
    !! and derivative share one Horner traversal, so they use the same support
    !! decision, clamping, coefficients and single-precision arithmetic. This
    !! routine intentionally supports only the KBWINSZ=1.5, KBALPHA=2 branch:
    !! apod_fast falls back to the separately approximated Bessel path for any
    !! other parameters, whose code-consistent derivative is not implemented.
    pure elemental subroutine apod_fast_value_deriv( self, x, value, dvalue_dx )
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: x
        real(sp),          intent(out) :: value, dvalue_dx
        real(sp), parameter :: coeffs(0:14) = [&
            1.0000000_sp, 1.0517297e1_sp, 2.7653385e1_sp, 3.2315430e1_sp, 2.1241936e1_sp,&
            8.9363103_sp, 2.6107175_sp, 5.6036106e-1_sp, 9.2085685e-2_sp, 1.1956698e-2_sp,&
            1.2575214e-3_sp, 1.0930353e-4_sp, 7.9831782e-6_sp, 4.9681336e-7_sp, 2.6658846e-8_sp]
        real(sp) :: dpdu, p, q, u
        integer  :: i
        if( abs(self%Whalf - 1.5_sp) > 1.e-6_sp .or. abs(self%alpha - 2.0_sp) > 1.e-6_sp )then
            error stop 'apod_fast_value_deriv supports only KBWINSZ=1.5, KBALPHA=2'
        endif
        if( abs(x) > self%Whalf )then
            value     = 0._sp
            dvalue_dx = 0._sp
            return
        endif
        q    = self%twooW * x
        u    = max(0._sp, 1._sp - q*q)
        p    = coeffs(14)
        dpdu = 0._sp
        do i = 13,0,-1
            ! Accumulate dP/du by Horner's rule.
            dpdu = dpdu * u + p
            ! Accumulate P(u) by Horner's rule.
            p    = p * u + coeffs(i)
        enddo
        ! Scale P(u) to the KB value.
        value     = self%oneoW * p
        ! Apply the chain rule dK/dx = (dK/du)(du/dx).
        dvalue_dx = self%oneoW * dpdu * (-2._sp * self%twooW * q)
    end subroutine apod_fast_value_deriv

    ! generate apodization fuction matrix in window
    pure subroutine apod_mat_2d(self, loc, iwinsz, wdim, kbw)
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: loc(2)
        integer,           intent(in)  :: iwinsz, wdim
        real(sp),          intent(out) :: kbw(wdim,wdim)
        integer  :: win_lo(2)
        real(sp) :: base(2), ww(2)
        real(sp) :: wrow(wdim), wcol(wdim)
        real(sp) :: recip
        integer  :: i
        ! Window lower corner
        win_lo = nint(loc) - iwinsz
        base   = real(win_lo, sp) - loc
        ! 1D weights along each axis (self%apod is pure elemental)
        do i = 1, wdim
            ww      = self%apod(base + real(i-1)) ! length-2 vector -> length-2 weights
            wrow(i) = ww(1)
            wcol(i) = ww(2)
        end do
        ! weights normalization
        wrow = wrow * (1.0_sp / sum(wrow))
        wcol = wcol * (1.0_sp / sum(wcol))
        ! Separable 2D outer product: kbw(i,j) = wrow(i) * wcol(j)
        do i = 1, wdim
            kbw(:,i) = wrow * wcol(i)
        end do
        ! Always normalize
        recip = 1.0_sp / sum(kbw)
        kbw   = kbw * recip
    end subroutine apod_mat_2d

    pure subroutine apod_mat_3d(self, loc, iwinsz, wdim, kbw)
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: loc(3)
        integer,           intent(in)  :: iwinsz, wdim
        real(sp),          intent(out) :: kbw(wdim,wdim,wdim)
        integer  :: win_lo(3)
        real(sp) :: base(3), ww(3)
        real(sp) :: wx(wdim), wy(wdim), wz(wdim)
        real(sp) :: recip
        integer  :: i, j
        ! Window lower corner
        win_lo = nint(loc) - iwinsz
        base   = real(win_lo, sp) - loc
        ! 1D weights along each axis (self%apod is pure elemental)
        do i = 1, wdim
            ww    = self%apod(base + real(i-1)) ! length-3 vector -> length-3 weights
            wx(i) = ww(1)
            wy(i) = ww(2)
            wz(i) = ww(3)
        end do
        ! Separable 3D outer product: kbw(i,j,k) = wx(i) * wy(j) * wz(k)
        wx = wx * (1.0_sp / sum(wx))
        wy = wy * (1.0_sp / sum(wy))
        wz = wz * (1.0_sp / sum(wz))
        do j = 1, wdim
            do i = 1, wdim
                kbw(:,i,j) = wx(:) * (wy(i) * wz(j))
            end do
        end do
        ! Always normalize
        recip = 1.0_sp / sum(kbw)
        kbw   = kbw * recip
    end subroutine apod_mat_3d

    !> Fast variant of apod_mat_2d: uses apod_fast (polynomial Bessel) and
    !! drops the redundant final renormalization (after 1-D normalizations
    !! sum(wrow)=sum(wcol)=1, so sum(kbw)=1 by construction).
    pure subroutine apod_mat_2d_fast(self, loc, iwinsz, wdim, kbw)
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: loc(2)
        integer,           intent(in)  :: iwinsz, wdim
        real(sp),          intent(out) :: kbw(wdim,wdim)
        integer  :: win_lo(2)
        real(sp) :: base(2), ww(2)
        real(sp) :: wrow(wdim), wcol(wdim)
        integer  :: i
        win_lo = nint(loc) - iwinsz
        base   = real(win_lo, sp) - loc
        do i = 1, wdim
            ww      = self%apod_fast(base + real(i-1))
            wrow(i) = ww(1)
            wcol(i) = ww(2)
        end do
        wrow = wrow * (1.0_sp / sum(wrow))
        wcol = wcol * (1.0_sp / sum(wcol))
        do i = 1, wdim
            kbw(:,i) = wrow * wcol(i)
        end do
    end subroutine apod_mat_2d_fast

    !> Fast variant of apod_mat_3d: uses apod_fast and drops the redundant
    !! final renormalization (separable 1-D normalizations already give
    !! sum(kbw) = sum(wx)*sum(wy)*sum(wz) = 1).
    pure subroutine apod_mat_3d_fast(self, loc, iwinsz, wdim, kbw)
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: loc(3)
        integer,           intent(in)  :: iwinsz, wdim
        real(sp),          intent(out) :: kbw(wdim,wdim,wdim)
        integer  :: win_lo(3)
        real(sp) :: base(3), ww(3)
        real(sp) :: wx(wdim), wy(wdim), wz(wdim)
        integer  :: i, j
        win_lo = nint(loc) - iwinsz
        base   = real(win_lo, sp) - loc
        do i = 1, wdim
            ww    = self%apod_fast(base + real(i-1))
            wx(i) = ww(1)
            wy(i) = ww(2)
            wz(i) = ww(3)
        end do
        wx = wx * (1.0_sp / sum(wx))
        wy = wy * (1.0_sp / sum(wy))
        wz = wz * (1.0_sp / sum(wz))
        do j = 1, wdim
            do i = 1, wdim
                kbw(:,i,j) = wx(:) * (wy(i) * wz(j))
            end do
        end do
    end subroutine apod_mat_3d_fast

    !>  \brief  Standard fast normalized 3-D KB stencil and fixed-cell gradient.
    !!
    !! dw_dloc(:,:,:,a) differentiates the normalized stencil with respect to
    !! loc(a), holding i0 fixed. switch_margin reports the distance to the next
    !! half-integer nint switch on each axis; it is zero at a switch. No
    !! derivative across a switch is implied.
    pure subroutine apod_mat_3d_fast_grad(self, loc, iwinsz, wdim, i0, switch_margin, kbw, dw_dloc)
        class(kbinterpol), intent(in)  :: self
        real(sp),          intent(in)  :: loc(3)
        integer,           intent(in)  :: iwinsz, wdim
        integer,           intent(out) :: i0(3)
        real(sp),          intent(out) :: switch_margin(3)
        real(sp),          intent(out) :: kbw(wdim,wdim,wdim)
        real(sp),          intent(out) :: dw_dloc(wdim,wdim,wdim,3)
        real(sp) :: base(3), draw(3), raw(3)
        real(sp) :: dx(wdim), dy(wdim), dz(wdim)
        real(sp) :: wx(wdim), wy(wdim), wz(wdim)
        real(sp) :: sx, sy, sz, dsx, dsy, dsz
        integer  :: i, j
        ! Select the fixed stencil lower corner.
        i0            = nint(loc) - iwinsz
        ! Measure the distance to the next stencil switch.
        switch_margin = 0.5_sp - abs(loc - real(nint(loc), sp))
        ! Form signed kernel arguments grid - loc.
        base           = real(i0, sp) - loc
        do i = 1, wdim
            ! Evaluate raw axis weights and dK/dx.
            call self%apod_fast_value_deriv(base + real(i-1, sp), raw, draw)
            wx(i) = raw(1)
            wy(i) = raw(2)
            wz(i) = raw(3)
            ! Convert dK/dx to dK/dloc.
            dx(i) = -draw(1)
            dy(i) = -draw(2)
            dz(i) = -draw(3)
        enddo
        ! Sum raw weights and derivatives on each axis.
        sx  = sum(wx); sy  = sum(wy); sz  = sum(wz)
        dsx = sum(dx); dsy = sum(dy); dsz = sum(dz)
        ! Normalize the three one-dimensional weight vectors.
        wx  = wx / sx; wy = wy / sy; wz = wz / sz
        ! Differentiate each normalized weight vector.
        dx  = (dx - wx * dsx) / sx
        dy  = (dy - wy * dsy) / sy
        dz  = (dz - wz * dsz) / sz
        do j = 1, wdim
            do i = 1, wdim
                ! Form the separable normalized 3-D stencil.
                kbw(:,i,j)          = wx(:) * (wy(i) * wz(j))
                ! Differentiate the x-axis factor.
                dw_dloc(:,i,j,1)    = dx(:) * (wy(i) * wz(j))
                ! Differentiate the y-axis factor.
                dw_dloc(:,i,j,2)    = wx(:) * (dy(i) * wz(j))
                ! Differentiate the z-axis factor.
                dw_dloc(:,i,j,3)    = wx(:) * (wy(i) * dz(j))
            enddo
        enddo
    end subroutine apod_mat_3d_fast_grad

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
                r = real(sinhc(arg2), kind=kind(r))
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

    !>  \brief returns the modified Bessel function I0(x) for any real x
    !! p.378 Handbook of Mathematical Functions, Abramowitz and Stegun
    elemental pure real(dp) function bessi0_dp( x )
        real(dp), intent(in) :: x
        real(dp) :: y, ax
        ax = x ! abs(x)  !! Assumption 1:  beta * sqrt(arg) is always positive
        if ( ax < 3.75d0 ) then
            y= x / 3.75d0
            y=y*y
            bessi0_dp=1.0d0+&
                y*(3.5156229d0 + y*(3.0899424d0 + y*(1.2067492d0 +&
                y*(0.2659732d0 + y*(0.0360768d0 + y* 0.0045813d0)))))
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
            bessi1 = x*(                                             &
                0.5d0 + y*(0.87890594d0 + y*(0.51498869d0 +          &
                y*(0.15084934d0 + y*(0.2658733d-1 + y*(0.301532d-2 + &
                y* 0.32411d-3))))))
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

    ! PUBLIC FUNCTIONS

    !>  \brief  is the Kaiser-Bessel apodization function
    !>  It must be kept identical to apod! Because it uses a non polymorphic
    !>  input it is compatible with OpenMP offload.
    !>  It is assumed the kbinterpol object has been constructed prior to call.
    pure real elemental function apod_device( self, x )
        type(kbinterpol), intent(in) :: self
        real,             intent(in) :: x
        real :: arg
        if( abs(x) > self%Whalf )then
            apod_device = 0.
            return ! for insignificant values return as soon as possible
        endif
        arg = self%twooW * x
        arg = 1. - arg * arg
        apod_device = real(self%oneoW * bessi0(self%beta * sqrt(arg)), kind=kind(apod_device))
    end function apod_device

    !>  \brief  is the Kaiser-Bessel apodization function
    !>  It must be kept identical to apod_device! Because it uses a non polymorphic
    !>  input it is compatible with OpenMP offload.
    !>  It is assumed the kbinterpol object has been constructed prior to call.
    pure elemental function apod_fast_device( self, x ) result( r )
        type(kbinterpol), intent(in) :: self
        real,             intent(in) :: x
        real :: r
        real(sp), parameter :: coeffs(0:14) = [&
            1.0000000_sp, 1.0517297e1_sp, 2.7653385e1_sp, 3.2315430e1_sp, 2.1241936e1_sp,&
            8.9363103_sp, 2.6107175_sp, 5.6036106e-1_sp, 9.2085685e-2_sp, 1.1956698e-2_sp,&
            1.2575214e-3_sp, 1.0930353e-4_sp, 7.9831782e-6_sp, 4.9681336e-7_sp, 2.6658846e-8_sp]
        real(sp) :: p, q, u
        integer  :: i
        if( abs(x) > self%Whalf )then
            r = 0.
            return
        endif
        if( abs(self%Whalf - 1.5) > 1.e-6 .or. abs(self%alpha - 2.0) > 1.e-6 )then
            r = apod_device(self, x)
            return
        endif
        q = self%twooW * x
        u = max(0._sp, 1._sp - q*q)
        p = coeffs(14)
        do i = 13,0,-1
            p = p * u + coeffs(i)
        enddo
        r = self%oneoW * p
    end function apod_fast_device

    !>  \brief  fast KB apodization for the common KBWINSZ=1.5, KBALPHA=2 kernel
    pure elemental function apod_kb15_a2( x ) result( r )
        real,    intent(in) :: x
        real(sp), parameter :: twooW = 2.0_sp / 3.0_sp
        real(sp), parameter :: oneoW = 1.0_sp / 3.0_sp
        real(sp), parameter :: coeffs(0:14) = [&
            1.0000000_sp, 1.0517297e1_sp, 2.7653385e1_sp, 3.2315430e1_sp, 2.1241936e1_sp,&
            8.9363103_sp, 2.6107175_sp, 5.6036106e-1_sp, 9.2085685e-2_sp, 1.1956698e-2_sp,&
            1.2575214e-3_sp, 1.0930353e-4_sp, 7.9831782e-6_sp, 4.9681336e-7_sp, 2.6658846e-8_sp]
        real     :: r
        real(sp) :: q, u
        if( x*x > 2.25_sp )then
            r = 0.
            return
        endif
        q = twooW * x
        u = max(0._sp, 1._sp - q*q)
        r = coeffs(14)
        r = r * u + coeffs(13)
        r = r * u + coeffs(12)
        r = r * u + coeffs(11)
        r = r * u + coeffs(10)
        r = r * u + coeffs(9)
        r = r * u + coeffs(8)
        r = r * u + coeffs(7)
        r = r * u + coeffs(6)
        r = r * u + coeffs(5)
        r = r * u + coeffs(4)
        r = r * u + coeffs(3)
        r = r * u + coeffs(2)
        r = r * u + coeffs(1)
        r = r * u + coeffs(0)
        r = r * oneoW
    end function apod_kb15_a2

    ! Public utility to get the windows pixel size without instantiating the object
    pure integer function kb_windim( winsz )
        real, intent(in) :: winsz
        integer          :: lb, ub
        call sqwin_1d(0., winsz, lb, ub)
        kb_windim = ub - lb  + 1
    end function kb_windim

end module simple_kbinterpol
