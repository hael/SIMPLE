!@descr: the polar shift map: one reference/particle pair's shift-dependent cross term over a disk of shifts, from the Fourier-Bessel expansion of the polar shift matrix
! The polar shift matrix exp(i(ax sx + ay sy)) on the half ring, written in polar shift coordinates
! s = rho (cos a, sin a), is exp(i z_k sin(theta_p - a)) with z_k = 2 pi k rho / box, so the cross term
! C(s) = sum_k w_k sum_p Re[ Q_k(p) exp(i(ax sx + ay sy)) ] equals, exactly (Jacobi-Anger),
! C(rho, a) = Re sum_n exp(-i n a) sum_k w_k J_n(z_k) H_k[n], H_k[n] = sum_p Q_k(p) exp(i n theta_p).
! The caller forms Q (cross_spectrum: the pair at one grid rotation) and the weights; the object holds the tables.
module simple_polar_shiftmap
use simple_defs,  only: sp, dp, DTWOPI
use simple_error, only: simple_exception
implicit none

public :: polar_shiftmap
private
#include "simple_local_flags.inc"

real(dp), parameter :: SHIFTMAP_JTOL = 1.d-8 !< a Bessel order whose value at the largest argument is below this is dropped

type :: polar_shiftmap
    private
    integer                  :: pftsz      = 0 !< angular samples over a half turn (the stored half ring)
    integer                  :: nrots      = 0 !< 2 * pftsz
    integer                  :: kfromto(2) = 0 !< shell range of the map
    integer                  :: box        = 0 !< logical box: shifts are in its pixels
    integer                  :: nrho       = 0 !< shift magnitudes of the grid
    integer                  :: nalpha     = 0 !< shift directions of the grid
    integer                  :: nharm_max  = 0 !< highest harmonic order kept at any magnitude
    real(dp),    allocatable :: rho(:)         !< the magnitudes (pixels)
    real(dp),    allocatable :: alpha(:)       !< the directions (radians)
    integer,     allocatable :: nharm(:)       !< highest harmonic order kept at each magnitude
    real(dp),    allocatable :: jn(:,:,:)      !< J_n(2 pi k rho / box): (0:nharm_max, kfromto(1):kfromto(2), nrho)
    complex(dp), allocatable :: tw(:,:)        !< exp(i n theta_p): (pftsz, 0:nharm_max)
    complex(dp), allocatable :: ea(:,:)        !< exp(-i n alpha_j): (nalpha, 0:nharm_max)
    logical                  :: exists = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: get_nrho
    procedure :: get_nalpha
    procedure :: get_nharm_max
    procedure :: get_nharm
    procedure :: get_rho
    procedure :: get_alpha
    procedure :: get_shift
    procedure :: cross_spectrum
    procedure :: coeffs
    procedure :: eval
    procedure :: map
    procedure :: peak
    procedure :: shell_weights
    procedure :: eval_rotcoeffs
    procedure :: map_rotcoeffs
end type polar_shiftmap

contains

    !> the tables for a polar geometry (half ring pftsz, shells kfromto, box) and a shift grid: magnitudes rho
    !! (pixels, any order, zero allowed) by nalpha directions (uniform over the turn)
    subroutine new( self, pftsz, kfromto, box, rho, nalpha )
        class(polar_shiftmap), intent(inout) :: self
        integer,               intent(in)    :: pftsz, kfromto(2), box, nalpha
        real(dp),              intent(in)    :: rho(:)
        real(dp) :: z, theta
        integer  :: ir, k, n, p, j
        call self%kill
        if( pftsz < 1 )                                 THROW_HARD('pftsz must be positive; new')
        if( kfromto(1) < 1 .or. kfromto(2) < kfromto(1) ) THROW_HARD('invalid kfromto; new')
        if( box < 2 .or. mod(box,2) /= 0 )              THROW_HARD('box must be even and positive; new')
        if( size(rho) < 1 )                             THROW_HARD('the shift grid needs a magnitude; new')
        if( any(rho < 0.d0) )                           THROW_HARD('shift magnitudes must be nonnegative; new')
        if( nalpha < 1 )                                THROW_HARD('the shift grid needs a direction; new')
        self%pftsz   = pftsz
        self%nrots   = 2 * pftsz
        self%kfromto = kfromto
        self%box     = box
        self%nrho    = size(rho)
        self%nalpha  = nalpha
        allocate(self%rho(self%nrho), source=rho)
        allocate(self%nharm(self%nrho), source=0)
        do ir = 1, self%nrho
            z = DTWOPI * real(kfromto(2),dp) * self%rho(ir) / real(box,dp)
            self%nharm(ir) = harmonic_order(z)
        enddo
        self%nharm_max = maxval(self%nharm)
        allocate(self%jn(0:self%nharm_max, kfromto(1):kfromto(2), self%nrho), source=0.d0)
        do ir = 1, self%nrho
            do k = kfromto(1), kfromto(2)
                z = DTWOPI * real(k,dp) * self%rho(ir) / real(box,dp)
                self%jn(:,k,ir) = bessel_jn(0, self%nharm_max, z)
            enddo
        enddo
        allocate(self%tw(pftsz, 0:self%nharm_max), self%ea(nalpha, 0:self%nharm_max))
        do p = 1, pftsz
            theta = real(p-1,dp) * DTWOPI / real(self%nrots,dp)
            do n = 0, self%nharm_max
                self%tw(p,n) = cmplx(cos(real(n,dp)*theta), sin(real(n,dp)*theta), kind=dp)
            enddo
        enddo
        allocate(self%alpha(nalpha))
        do j = 1, nalpha
            self%alpha(j) = real(j-1,dp) * DTWOPI / real(nalpha,dp)
            do n = 0, self%nharm_max
                self%ea(j,n) = cmplx(cos(real(n,dp)*self%alpha(j)), -sin(real(n,dp)*self%alpha(j)), kind=dp)
            enddo
        enddo
        self%exists = .true.
    end subroutine new

    !> the highest order kept at argument z: J_n(z) is positive and decreasing in n past z, so the orders from
    !! the first one below SHIFTMAP_JTOL on are dropped
    integer function harmonic_order( z )
        real(dp), intent(in) :: z
        integer :: n
        n = max(0, ceiling(z))
        do while( abs(bessel_jn(n, z)) >= SHIFTMAP_JTOL )
            n = n + 1
        enddo
        harmonic_order = max(0, n - 1)
    end function harmonic_order

    subroutine kill( self )
        class(polar_shiftmap), intent(inout) :: self
        if( self%exists )then
            deallocate(self%rho, self%alpha, self%nharm, self%jn, self%tw, self%ea)
            self%pftsz     = 0
            self%nrots     = 0
            self%kfromto   = 0
            self%box       = 0
            self%nrho      = 0
            self%nalpha    = 0
            self%nharm_max = 0
            self%exists    = .false.
        endif
    end subroutine kill

    pure integer function get_nrho( self )
        class(polar_shiftmap), intent(in) :: self
        get_nrho = self%nrho
    end function get_nrho

    pure integer function get_nalpha( self )
        class(polar_shiftmap), intent(in) :: self
        get_nalpha = self%nalpha
    end function get_nalpha

    pure integer function get_nharm_max( self )
        class(polar_shiftmap), intent(in) :: self
        get_nharm_max = self%nharm_max
    end function get_nharm_max

    pure integer function get_nharm( self, irho )
        class(polar_shiftmap), intent(in) :: self
        integer,               intent(in) :: irho
        get_nharm = self%nharm(irho)
    end function get_nharm

    pure real(dp) function get_rho( self, irho )
        class(polar_shiftmap), intent(in) :: self
        integer,               intent(in) :: irho
        get_rho = self%rho(irho)
    end function get_rho

    pure real(dp) function get_alpha( self, ialpha )
        class(polar_shiftmap), intent(in) :: self
        integer,               intent(in) :: ialpha
        get_alpha = self%alpha(ialpha)
    end function get_alpha

    !> the shift of a grid point, in pixels, in the frame the shift matrix acts in (the reference's own)
    pure function get_shift( self, irho, ialpha ) result( s )
        class(polar_shiftmap), intent(in) :: self
        integer,               intent(in) :: irho, ialpha
        real(dp) :: s(2)
        s(1) = self%rho(irho) * cos(self%alpha(ialpha))
        s(2) = self%rho(irho) * sin(self%alpha(ialpha))
    end function get_shift

    !> Q of the pair (reference, particle, |CTF|) at grid rotation irot, in the reference's own indexing: the
    !! reference is shifted in its frame and then rotated, so the shift matrix multiplies entry rp of the
    !! reference, which the rotation lands on ring position q = rp + irot - 1 (modulo the turn) of the particle;
    !! positions past the stored half are the Friedel mates (conjugates) of the first half
    subroutine cross_spectrum( self, ref, ptcl, ctf, irot, Q )
        class(polar_shiftmap), intent(in)  :: self
        complex(sp),           intent(in)  :: ref(self%pftsz,  self%kfromto(1):self%kfromto(2))
        complex(sp),           intent(in)  :: ptcl(self%pftsz, self%kfromto(1):self%kfromto(2))
        real(sp),              intent(in)  :: ctf(self%pftsz,  self%kfromto(1):self%kfromto(2))
        integer,               intent(in)  :: irot
        complex(dp),           intent(out) :: Q(self%pftsz, self%kfromto(1):self%kfromto(2))
        integer :: rp, q_full, q_half, k
        if( irot < 1 .or. irot > self%nrots ) THROW_HARD('rotation index out of range; cross_spectrum')
        do rp = 1, self%pftsz
            q_full = modulo(rp + irot - 2, self%nrots) + 1
            if( q_full <= self%pftsz )then
                q_half = q_full
                do k = self%kfromto(1), self%kfromto(2)
                    Q(rp,k) = cmplx(ref(rp,k),kind=dp) * real(ctf(q_half,k),dp) * conjg(cmplx(ptcl(q_half,k),kind=dp))
                enddo
            else
                q_half = q_full - self%pftsz
                do k = self%kfromto(1), self%kfromto(2)
                    Q(rp,k) = cmplx(ref(rp,k),kind=dp) * real(ctf(q_half,k),dp) * cmplx(ptcl(q_half,k),kind=dp)
                enddo
            endif
        enddo
    end subroutine cross_spectrum

    !> the Fourier coefficients of the map in the shift direction, a_n(rho) = sum_k w_k J_n(z_k(rho)) H_k[n],
    !! for |n| <= nharm(rho) (zero beyond): A(-nharm_max:nharm_max, nrho)
    subroutine coeffs( self, Q, w, A )
        class(polar_shiftmap), intent(in)  :: self
        complex(dp),           intent(in)  :: Q(self%pftsz, self%kfromto(1):self%kfromto(2))
        real(dp),              intent(in)  :: w(self%kfromto(1):self%kfromto(2))
        complex(dp),           intent(out) :: A(-self%nharm_max:self%nharm_max, self%nrho)
        complex(dp), allocatable :: hp(:,:), hm(:,:)
        complex(dp) :: ap, am
        integer     :: k, n, ir
        allocate(hp(self%kfromto(1):self%kfromto(2), 0:self%nharm_max),&
                &hm(self%kfromto(1):self%kfromto(2), 0:self%nharm_max))
        ! the half-ring harmonics of every shell, positive and negative orders
        do n = 0, self%nharm_max
            do k = self%kfromto(1), self%kfromto(2)
                hp(k,n) = sum(Q(:,k) *       self%tw(:,n))
                hm(k,n) = sum(Q(:,k) * conjg(self%tw(:,n)))
            enddo
        enddo
        ! the Bessel-weighted shell sums; J_{-n} = (-1)^n J_n
        A = cmplx(0.d0, 0.d0, kind=dp)
        do ir = 1, self%nrho
            do n = 0, self%nharm(ir)
                ap = cmplx(0.d0, 0.d0, kind=dp)
                am = cmplx(0.d0, 0.d0, kind=dp)
                do k = self%kfromto(1), self%kfromto(2)
                    ap = ap + w(k) * self%jn(n,k,ir) * hp(k,n)
                    am = am + w(k) * self%jn(n,k,ir) * hm(k,n)
                enddo
                A(n,ir) = ap
                if( n > 0 )then
                    if( mod(n,2) == 0 )then
                        A(-n,ir) =  am
                    else
                        A(-n,ir) = -am
                    endif
                endif
            enddo
        enddo
        deallocate(hp, hm)
    end subroutine coeffs

    !> the cross term at magnitude rho(irho) and any direction alpha (radians)
    real(dp) function eval( self, A, irho, alpha )
        class(polar_shiftmap), intent(in) :: self
        complex(dp),           intent(in) :: A(-self%nharm_max:self%nharm_max, self%nrho)
        integer,               intent(in) :: irho
        real(dp),              intent(in) :: alpha
        complex(dp) :: e, en
        integer     :: n
        e    = cmplx(cos(alpha), -sin(alpha), kind=dp) ! exp(-i alpha)
        en   = cmplx(1.d0, 0.d0, kind=dp)
        eval = real(A(0,irho), dp)
        do n = 1, self%nharm(irho)
            en   = en * e
            eval = eval + real(A(n,irho) * en + A(-n,irho) * conjg(en), dp)
        enddo
    end function eval

    !> the cross term on the whole grid: cmap(nalpha, nrho)
    subroutine map( self, A, cmap )
        class(polar_shiftmap), intent(in)  :: self
        complex(dp),           intent(in)  :: A(-self%nharm_max:self%nharm_max, self%nrho)
        real(dp),              intent(out) :: cmap(self%nalpha, self%nrho)
        integer :: ir, n, j
        do ir = 1, self%nrho
            cmap(:,ir) = real(A(0,ir), dp)
            do n = 1, self%nharm(ir)
                do j = 1, self%nalpha
                    cmap(j,ir) = cmap(j,ir) + real(A(n,ir) * self%ea(j,n) + A(-n,ir) * conjg(self%ea(j,n)), dp)
                enddo
            enddo
        enddo
    end subroutine map

    !> the grid point of the largest cross term (the first on ties)
    subroutine peak( self, cmap, irho, ialpha, val )
        class(polar_shiftmap), intent(in)  :: self
        real(dp),              intent(in)  :: cmap(self%nalpha, self%nrho)
        integer,               intent(out) :: irho, ialpha
        real(dp),              intent(out) :: val
        integer :: loc(2)
        loc    = maxloc(cmap)
        ialpha = loc(1)
        irho   = loc(2)
        val    = cmap(ialpha, irho)
    end subroutine peak

    !> the Bessel-weighted shell weights of order n at magnitude rho(irho): wj(k) = w(k) J_n(2 pi k rho / box)
    pure subroutine shell_weights( self, w, n, irho, wj )
        class(polar_shiftmap), intent(in)  :: self
        real(dp),              intent(in)  :: w(self%kfromto(1):self%kfromto(2))
        integer,               intent(in)  :: n, irho
        real(dp),              intent(out) :: wj(self%kfromto(1):self%kfromto(2))
        wj = w * self%jn(n,:,irho)
    end subroutine shell_weights

    !> the cross term from the one-sided coefficients of the all-rotation kernel at one rotation, a(n) = a_n(irot, rho)
    !! for n = 0..nharm(irho) (full-ring harmonics: the negative orders are the conjugates, the full ring is twice the
    !! half ring), at magnitude rho(irho) and any direction alpha: C = Re(a_0)/2 + sum_n Re(exp(-i n alpha) a_n)
    pure real(dp) function eval_rotcoeffs( self, a, irho, alpha )
        class(polar_shiftmap), intent(in) :: self
        complex(sp),           intent(in) :: a(0:self%nharm_max)
        integer,               intent(in) :: irho
        real(dp),              intent(in) :: alpha
        complex(dp) :: e, en
        integer     :: n
        e  = cmplx(cos(alpha), -sin(alpha), kind=dp)
        en = cmplx(1.d0, 0.d0, kind=dp)
        eval_rotcoeffs = 0.5d0 * real(a(0), dp)
        do n = 1, self%nharm(irho)
            en = en * e
            eval_rotcoeffs = eval_rotcoeffs + real(cmplx(a(n),kind=dp) * en, dp)
        enddo
    end function eval_rotcoeffs

    !> the cross term at magnitude rho(irho) on the grid's directions for every rotation, from the one-sided
    !! coefficients a(irot, n) of the all-rotation kernel: cmap(irot, ialpha)
    pure subroutine map_rotcoeffs( self, a, irho, cmap )
        class(polar_shiftmap), intent(in)  :: self
        complex(sp),           intent(in)  :: a(self%nrots, 0:self%nharm_max)
        integer,               intent(in)  :: irho
        real(dp),              intent(out) :: cmap(self%nrots, self%nalpha)
        integer :: n, j
        do j = 1, self%nalpha
            cmap(:,j) = 0.5d0 * real(a(:,0), dp)
        enddo
        do n = 1, self%nharm(irho)
            do j = 1, self%nalpha
                cmap(:,j) = cmap(:,j) + real(cmplx(a(:,n),kind=dp) * self%ea(j,n), dp)
            enddo
        enddo
    end subroutine map_rotcoeffs

end module simple_polar_shiftmap
