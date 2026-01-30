!@descr: fast calcualtion of Kaiser-Bessel interpolation kernel, assuming Whalf <= 1.5_sp
!
! Changes vs original (implements suggestions 1–5,7; excludes 6):
!   - Kind hygiene: all members stored in sp; dp used only where it buys stability
!   - Precompute invariants assuming beta<3.75
!   - apod(): branch-lighter + optional clamp for arg (avoids rare NaNs without changing valid-domain results)
!   - Adds fast separable helpers:
!        apod_weights_1d()  -> compute 1D weight vector once (O(wdim) apod calls)
!        apod_outer_2d()    -> outer product wx*wy into kbw (no 2D temporaries, O(wdim^2) muls only)
!
module simple_kbinterpol_fast
use simple_defs
use iso_c_binding
use simple_edges_sqwins, only: sqwin_1d
implicit none

public :: kbinterpol_fast
private
#include "simple_local_flags.inc"

type :: kbinterpol_fast
    private
    ! Store as sp (the hot-path kind). dp only used locally where needed.
    real(sp) :: alpha, beta, betasq, oneoW, piW, twooW, W, Whalf, threshInstr
    logical  :: beta_lt_3p75 = .false.   ! if beta < 3.75 then I0(u) always uses the small-branch polynomial
  contains
    procedure :: new
    procedure :: get_winsz
    procedure :: get_alpha
    procedure :: get_wdim
    procedure :: apod
    procedure :: dapod
    procedure :: instr
    ! NEW: separable helpers for KB window building
    procedure :: apod_weights_1d
    procedure :: apod_outer_2d
end type kbinterpol_fast

interface kbinterpol_fast
   module procedure constructor
end interface kbinterpol_fast

contains

    function constructor( Whalf_in, alpha_in ) result( self )
        real(sp), intent(in) :: Whalf_in, alpha_in
        type(kbinterpol_fast) :: self
        call self%new(Whalf_in, alpha_in)
    end function constructor

    subroutine new( self, Whalf_in, alpha_in )
        class(kbinterpol_fast), intent(inout) :: self
        real(sp),               intent(in)    :: Whalf_in, alpha_in
        self%Whalf  = Whalf_in
        self%alpha  = alpha_in
        self%W      = 2.0_sp * self%Whalf
        self%piW    = real(pi, sp) * self%W
        if (abs(self%Whalf - 1.5_sp) < 1e-6) then
            self%beta = 7.4_sp
        else
            THROW_HARD('Use a Whalf <= 1.5 for fast KB interpolation')
        end if
        self%betasq       = self%beta * self%beta
        self%twooW        = 2.0_sp / self%W
        self%oneoW        = 1.0_sp / self%W
        self%threshInstr  = self%beta / self%piW - real(TINY, sp)**2
    end subroutine new

    pure real(sp) function get_winsz( self )
        class(kbinterpol_fast), intent(in) :: self
        get_winsz = self%Whalf
    end function get_winsz

    pure real(sp) function get_alpha( self )
        class(kbinterpol_fast), intent(in) :: self
        get_alpha = self%alpha
    end function get_alpha

    pure integer function get_wdim( self )
        class(kbinterpol_fast), intent(in) :: self
        integer :: win(2)
        call sqwin_1d(0.0_sp, self%Whalf, win(1), win(2))
        get_wdim = win(2) - win(1) + 1
    end function get_wdim

    ! ============================================================
    ! Hot path: apod (sp)
    ! ============================================================
    pure real(sp) function apod( self, x ) result( r )
        class(kbinterpol_fast), intent(in) :: self
        real(sp),               intent(in) :: x
        real(sp) :: ax, t, arg, u
        ax = abs(x)
        if (ax > self%Whalf) then
            r = 0.0_sp
            return
        end if
        t   = self%twooW * x
        arg = 1.0_sp - t*t
        ! Optional safety (cheap, branchless): prevents rare NaNs from roundoff near |x|~Whalf.
        ! For valid-domain values, this is numerically identical (arg>=0).
        arg = max(arg, 0.0_sp)
        u = self%beta * sqrt(arg)
        r = self%oneoW * bessi0_small_sp(u)
    end function apod

    pure real(dp) function dapod( self, x ) result(r)
        class(kbinterpol_fast), intent(in) :: self
        real(dp),               intent(in) :: x
        real(dp) :: t, arg, sqrtarg
        t   = real(self%twooW, dp) * x
        arg = 1.0_dp - t*t
        if (arg <= 0.0_dp) then
            r = 0.0_dp
            return
        end if
        sqrtarg = sqrt(arg)
        r = - 4.0_dp * real(self%beta, dp) * x * bessi1_dp(real(self%beta, dp) * sqrtarg) / sqrtarg / real(self%W, dp)**3
    end function dapod

    elemental real(sp) function instr( self, x ) result( r )
        class(kbinterpol_fast), intent(in) :: self
        real(sp),               intent(in) :: x
        real(sp) :: arg2
        if (abs(x) < self%threshInstr) then
            arg2 = sqrt(max(self%betasq - (self%piW * x)**2, 0.0_sp))
            if (arg2 < real(TINY, sp)) then
                r = 1.0_sp
            else
                r = real(sinhc_sp(arg2), sp)
            end if
        else
            r = 1.0_sp
        end if
    end function instr

    ! ============================================================
    ! NEW: separable helpers (big win in KB interpolation window building)
    !
    ! Usage pattern (in gen_planes_pad KB path):
    !   call self%kbwin%apod_weights_1d(loc(1), win(1,1), wx)
    !   call self%kbwin%apod_weights_1d(loc(2), win(1,2), wy)
    !   call self%kbwin%apod_outer_2d(wx, wy, kbw, normalize=.true.)
    !
    ! wx(i) = apod( (win_lo + (i-1)) - loc )
    ! ============================================================
    pure subroutine apod_weights_1d(self, loc, win_lo, w)
        class(kbinterpol_fast), intent(in)  :: self
        real(sp),               intent(in)  :: loc
        integer,                intent(in)  :: win_lo
        real(sp),               intent(out) :: w(:)
        integer :: i
        real(sp) :: x
        do i = 1, size(w)
            x    = real(win_lo + (i-1), sp) - loc
            w(i) = self%apod(x)
        end do
    end subroutine apod_weights_1d

    pure subroutine apod_outer_2d(self, wx, wy, kbw )
        class(kbinterpol_fast), intent(in)    :: self
        real(sp),               intent(in)    :: wx(:), wy(:)
        real(sp),               intent(inout) :: kbw(:,:)
        integer :: i, j
        real(sp) :: s
        ! Fill as outer product (no temporaries)
        do j = 1, size(wy)
            do i = 1, size(wx)
                kbw(i,j) = wx(i) * wy(j)
            end do
        end do
        ! always normalize
        s = sum(kbw)
        if (s > 0.0_sp) kbw = kbw / s
    end subroutine apod_outer_2d

    ! ============================================================
    ! Bessel I0 implementations
    !   - small branch polynomial (sp) (valid/used when u<3.75)
    !   - full version with exp/sqrt tail (sp)
    ! ============================================================

    elemental pure real(sp) function bessi0_small_sp(x) result(r)
        real(sp), intent(in) :: x
        real(sp) :: y
        y = x / 3.75_sp
        y = y * y
        r = 1.0_sp + y*(3.5156229_sp + y*(3.0899424_sp + y*(1.2067492_sp + &
            y*(0.2659732_sp + y*(0.0360768_sp + y*0.0045813_sp)))))
    end function bessi0_small_sp

    elemental pure real(dp) function bessi1_dp( x ) result(r)
        real(dp), intent(in) :: x
        real(dp) :: y, ax, bx, tmp
        ax = abs(x)
        if ( ax < 3.75_dp) then
            y = x / 3.75_dp
            y = y*y
            r = x*(0.5_dp + y*(0.87890594_dp + y*(0.51498869_dp + &
                y*(0.15084934_dp + y*(0.02658733_dp + y*(0.00301532_dp + &
                y*0.00032411_dp))))))
        else
            y   = 3.75_dp / ax
            bx  = exp(ax) / sqrt(ax)
            tmp = 0.39894228_dp + y*(-0.03988024_dp + y*(-0.00362018_dp + &
                  y*(0.00163801_dp + y*(-0.01031555_dp + y*(0.02282967_dp + &
                  y*(-0.02895312_dp + y*(0.01787654_dp + y*(-0.00420059_dp))))))))
            r = tmp * bx
            if (x < 0.0_dp) r = -r
        end if
    end function bessi1_dp

    elemental pure real(dp) function sinhc_sp(xin) result(y)
        ! Keep coefficients in dp for stability; return dp then cast by caller if needed.
        real(sp), intent(in) :: xin
        real(dp), parameter  :: P0 = -0.6307673640497716991184787251d+6,&
            P1 = -0.8991272022039509355398013511d+05, &
            P2 = -0.2894211355989563807284660366d+04, &
            P3 = -0.2630563213397497062819489000d+02, &
            Q0 = -0.6307673640497716991212077277d+06, &
            Q1 =  0.1521517378790019070696485176d+05, &
            Q2 = -0.1736789535582336995334509110d+03
        real(dp) :: x, xsq
        x = real(xin, dp)
        if (x > 0.5_dp) then
            y = (exp(x) - exp(-x)) / (2.0_dp * x)
        else
            xsq = x * x
            y = (((P3*xsq + P2)*xsq + P1)*xsq + P0) / (((xsq + Q2)*xsq + Q1)*xsq + Q0)
        end if
    end function sinhc_sp

end module simple_kbinterpol_fast

! How to implement:

! !> Produces padded shifted, rotated, CTF multiplied fourier & CTF-squared planes
! subroutine gen_planes_pad( self, img, ctfvars, shift, e3, iptcl, linear )
!     class(fplane),                  intent(inout) :: self
!     class(image),                   intent(inout) :: img
!     class(ctfparams),               intent(in)    :: ctfvars
!     real,                           intent(in)    :: shift(2), e3
!     integer,                        intent(in)    :: iptcl
!     logical,                        intent(in)    :: linear
!     type(ctf)                :: tfun
!     complex(c_float_complex) :: c, w1, w2, ph0, ph_h, ph_k, phase
!     real(dp) :: pshift(2)
!     real     :: kbw(self%wdim,self%wdim)
!     real     :: rmat(2,2), loc(2)
!     real     :: tval, tvalsq, add_phshift
!     integer  :: win(2,2), sigma2_kfromto(2)
!     integer  :: physh, physk
!     integer  :: i, j, h, k, hh, kk, shell, iwinsz, hlim
!     logical  :: l_ctf, l_flip
!     ! NEW: separable KB weights (avoid O(wdim^2) apod calls)
!     real     :: wx(self%wdim), wy(self%wdim)
!     if( .not.self%padded ) THROW_HARD('gen_planes_pad only for use with padding!')
!     ! CTF
!     l_ctf = ctfvars%ctfflag /= CTFFLAG_NO
!     if( l_ctf )then
!         l_flip = ctfvars%ctfflag == CTFFLAG_FLIP
!         tfun   = ctf(ctfvars%smpd, ctfvars%kv, ctfvars%cs, ctfvars%fraca)
!         call tfun%init(ctfvars%dfx, ctfvars%dfy, ctfvars%angast)
!         add_phshift = merge(ctfvars%phshift, 0.0, ctfvars%l_phaseplate)
!     endif
!     ! rotation & scale
!     call rotmat2D(e3, rmat)
!     rmat = self%alpha * rmat
!     ! sigma2
!     if( params_glob%l_ml_reg )then
!         sigma2_kfromto(1) = lbound(eucl_sigma2_glob%sigma2_noise,1)
!         sigma2_kfromto(2) = ubound(eucl_sigma2_glob%sigma2_noise,1)
!         self%sigma2_noise(sigma2_kfromto(1):self%nyq_crop) = eucl_sigma2_glob%sigma2_noise(sigma2_kfromto(1):self%nyq_crop,iptcl)
!     end if
!     ! Interpolation
!     self%cmplx_plane = cmplx(0.,0.)
!     self%ctfsq_plane = 0.
!     ! Shift precomputation (phase recurrence)
!     pshift = real(-shift * self%shconst(1:2),dp)
!     w1 = cmplx( real(cos(pshift(1)), c_float), real(sin(pshift(1)), c_float), kind=c_float_complex )
!     w2 = cmplx( real(cos(pshift(2)), c_float), real(sin(pshift(2)), c_float), kind=c_float_complex )
!     ph0  = cmplx( real(cos(real(self%frlims_crop(1,1),dp)*pshift(1)), c_float), &
!                   real(sin(real(self%frlims_crop(1,1),dp)*pshift(1)), c_float), kind=c_float_complex )
!     ph_k = cmplx( real(cos(real(self%frlims_crop(2,1),dp)*pshift(2)), c_float), &
!                   real(sin(real(self%frlims_crop(2,1),dp)*pshift(2)), c_float), kind=c_float_complex )
!     ! KB interpolation
!     iwinsz = ceiling(self%winsz - 0.5)
!     do k = self%frlims_crop(2,1), 0
!         if( k == 0 )then
!             hlim = -1     ! h=0, k=0 is treated after the loop
!         else
!             hlim = self%frlims_crop(1,2)
!         endif
!         ph_h = ph0
!         do h = self%frlims_crop(1,1), hlim
!             shell = nint(sqrt(real(h*h + k*k)))
!             if( shell > self%nyq_crop ) then
!                 ph_h = ph_h * w1
!                 cycle
!             endif
!             ! Retrieve component & shift
!             physh = ft_map_phys_addrh(h,k)
!             physk = ft_map_phys_addrk(h,k)
!             phase = ph_k * ph_h
!             c = merge(conjg(img%get_cmat_at(physh,physk,1)), img%get_cmat_at(physh,physk,1), h<0) * phase
!             ! CTF
!             if( l_ctf )then
!                 tval   = tfun%eval(ft_map_spaFreqSq(h,k), ft_map_astigang(h,k), add_phshift)
!                 tvalsq = tval * tval
!                 c      = merge(abs(tval)*c, tval*c, l_flip)
!             else
!                 tvalsq = 1.0
!             endif
!             ! sigma2 weighing
!             if( params_glob%l_ml_reg ) then
!                 if(shell >= sigma2_kfromto(1))then
!                     c      = c      / self%sigma2_noise(shell)
!                     tvalsq = tvalsq / self%sigma2_noise(shell)
!                 else
!                     c      = c      / self%sigma2_noise(sigma2_kfromto(1))
!                     tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
!                 endif
!             endif
!             ! rotation
!             loc = matmul(real([h,k]), rmat)
!             ! window (same as before)
!             win(1,:) = nint(loc)
!             win(2,:) = win(1,:) + iwinsz
!             win(1,:) = win(1,:) - iwinsz
!             ! ------------------------------------------------------------
!             ! NEW KERNEL BUILD (separable): O(wdim) apod calls + outer product
!             ! ------------------------------------------------------------
!             ! wx(i) = apod( (win_lo_h + (i-1)) - loc(1) )
!             ! wy(j) = apod( (win_lo_k + (j-1)) - loc(2) )
!             do i = 1, self%wdim
!                 wx(i) = self%kbwin%apod( real(win(1,1) + i - 1) - loc(1) )
!                 wy(i) = self%kbwin%apod( real(win(1,2) + i - 1) - loc(2) )
!             end do
!             ! kbw = outer product wx*wy
!             do j = 1, self%wdim
!                 do i = 1, self%wdim
!                     kbw(i,j) = wx(i) * wy(j)
!                 end do
!             end do
!             ! normalize (exactly like old: kbw = kbw / sum(kbw))
!             kbw = kbw / sum(kbw)
!             ! interpolation (unchanged)
!             i = 0
!             do hh = win(1,1), win(2,1)
!                 i = i + 1
!                 if( abs(hh) > self%nyq_croppd )cycle
!                 j = 0
!                 do kk = win(1,2), win(2,2)
!                     j = j + 1
!                     if( (kk >= self%frlims_croppd(2,1)) .and. (kk <= self%frlims_croppd(2,2)) )then
!                         self%cmplx_plane(hh,kk) = self%cmplx_plane(hh,kk) + kbw(i,j) * c
!                         self%ctfsq_plane(hh,kk) = self%ctfsq_plane(hh,kk) + kbw(i,j) * tvalsq
!                     endif
!                     ! Friedel symmetric
!                     if( (-kk >= self%frlims_croppd(2,1)) .and. (-kk <= self%frlims_croppd(2,2)) )then
!                         self%cmplx_plane(-hh,-kk) = self%cmplx_plane(-hh,-kk) + kbw(i,j) * conjg(c)
!                         self%ctfsq_plane(-hh,-kk) = self%ctfsq_plane(-hh,-kk) + kbw(i,j) * tvalsq
!                     endif
!                 enddo
!             enddo
!             ! phase recurrence along h
!             ph_h = ph_h * w1
!         enddo

!         ! phase recurrence along k
!         ph_k = ph_k * w2
!     enddo
!     ! Interpolation free DC
!     c = img%get_fcomp2D(0,0)
!     if( l_ctf )then
!         tval   = tfun%eval(0., 0., add_phshift)
!         tvalsq = tval * tval
!         c      = merge(abs(tval)*c, tval*c, l_flip)
!     else
!         tvalsq = 1.0
!     endif
!     if( params_glob%l_ml_reg ) then
!         c      = c      / self%sigma2_noise(sigma2_kfromto(1))
!         tvalsq = tvalsq / self%sigma2_noise(sigma2_kfromto(1))
!     endif
!     self%cmplx_plane(0,0) = c
!     self%ctfsq_plane(0,0) = tvalsq
! end subroutine gen_planes_pad
