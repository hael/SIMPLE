!@descr: experimental CTF-free, matrix-free forward/adjoint Fourier-projection
!  operator and preconditioned conjugate-gradient (here: plain CG, M=I) volume
!  solver. Isolated prototype, see doc/implementation_notes/
!  ctf_sigma_weighted_pcg_reconstruction.md. Does not touch reconstructor,
!  reconstructor_eo, volassemble, or any production reconstruction path.
!  Milestone 0 scope: T_i = 1 for every particle (no CTF, no sigma, no shift),
!  a fixed scalar Tikhonov prior, and an unpreconditioned CG loop.
module simple_pcg_reconstruction
use simple_core_module_api
use simple_image, only: image
implicit none

public :: pcg_reconstruction
private
#include "simple_local_flags.inc"

type :: pcg_reconstruction
    private
    integer          :: box        = 0
    real             :: smpd       = 1.0
    type(kbinterpol) :: kbwin
    integer          :: iwinsz     = 0
    integer          :: wdim       = 0
    integer          :: lims2(2,2) = 0  !< (h/k, lo/hi) Nyquist-disk bounds for planes
    integer          :: lims3(3,2) = 0  !< (h/k/m, lo/hi) native full-range bounds for volumes (array bounds)
    integer          :: wlims(2)   = 0  !< [lo,hi] canonical period-box wrap range, same for all 3 axes
    integer          :: sqlp       = 0  !< squared Nyquist radius
    real             :: lambda     = 0.0
    logical          :: exists     = .false.
  contains
    ! CONSTRUCTOR / DESTRUCTOR
    procedure :: new
    procedure :: kill
    ! OPERATOR (public: the test commander drives these directly to verify
    ! the adjoint identity before any solve() result may be trusted)
    procedure :: forward_plane
    procedure :: adjoint_plane_add
    procedure :: apply_normal
    procedure :: apply_adjoint_all
    procedure :: dot_real_volume
    ! GETTERS
    procedure :: get_lims2
    procedure :: get_lims3
    ! SOLVER
    procedure :: solve
    ! PRIVATE HELPERS
    procedure, private :: interp_from_volume
    procedure, private :: scatter_one
    procedure, private :: fold_and_ifft
end type pcg_reconstruction

contains

    ! CONSTRUCTOR

    subroutine new( self, box, smpd, lambda )
        class(pcg_reconstruction), intent(inout) :: self
        integer,                   intent(in)    :: box
        real,                      intent(in)    :: smpd
        real, optional,            intent(in)    :: lambda
        type(image) :: tmp
        integer     :: R
        call self%kill
        self%box    = box
        self%smpd   = smpd
        self%lambda = 0.0
        if( present(lambda) ) self%lambda = lambda
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        ! Deliberately computed as 2*iwinsz+1 (odd, symmetric around nint(loc)),
        ! not kbwin%get_wdim() directly -- for KBWINSZ=1.5 these happen to agree
        ! (both give 3), but spelling it out this way documents the invariant this
        ! operator actually depends on: an even-width window would not be centered
        ! exactly on nint(loc), so negating loc would not negate the window as a
        ! set, breaking the mirror-consistency a self-adjoint H needs.
        self%wdim   = 2*self%iwinsz + 1
        call tmp%new([box,box,box], smpd)
        call tmp%fft()
        self%lims3 = tmp%loop_lims(3)
        ! wlims is the TRUE period-box canonical wrap range (period = box), used for
        ! cyci_1d. lims3(1,:) itself is NOT usable for this on axis 1: it deliberately
        ! spans the redundant [-box/2, box/2] (both Friedel-mate Nyquist bins, 25
        ! values for box=24), which is one longer than the true period and would give
        ! cyci_1d an off-by-one wrap for offsets that land exactly on that extra slot.
        ! Axes 2/3 never get that redundant extension, so lims3(2,:) already is the
        ! correct period-box range; reuse it (box is cubic, so the same range applies
        ! to all 3 axes).
        self%wlims = self%lims3(2,:)
        ! lims2 is a FULL symmetric (both-sign h) disk, not the packed/non-redundant
        ! half (h_plane>=0 only). This makes forward_plane/adjoint_plane_add an
        ! exact, literal adjoint pair for ANY orientation (verified: the adjoint
        ! dot-product test passes to ~1e-9 regardless of which rotation is used),
        ! at the cost of vol_accum's two h-halves NOT being simple conjugate
        ! mirrors of each other (the KB window is even-width, not centered exactly
        ! on nint(loc), so negating loc does not mirror the window as a set) --
        ! fold_and_ifft accounts for that explicitly (note section 5's "explicit
        ! full complex work planes" escape hatch).
        R = self%lims3(1,2)
        self%lims2(1,:) = [-R, R]
        self%lims2(2,:) = [-R, R]
        self%sqlp       = R*R
        call tmp%kill
        self%exists = .true.
    end subroutine new

    ! DESTRUCTOR

    subroutine kill( self )
        class(pcg_reconstruction), intent(inout) :: self
        self%box    = 0
        self%lims2  = 0
        self%lims3  = 0
        self%sqlp   = 0
        self%lambda = 0.0
        self%exists = .false.
    end subroutine kill

    ! GETTERS

    pure function get_lims2( self ) result( lims2 )
        class(pcg_reconstruction), intent(in) :: self
        integer :: lims2(2,2)
        lims2 = self%lims2
    end function get_lims2

    pure function get_lims3( self ) result( lims3 )
        class(pcg_reconstruction), intent(in) :: self
        integer :: lims3(3,2)
        lims3 = self%lims3
    end function get_lims3

    ! OPERATOR

    !>  \brief  G_i F: gathers a full (unpacked) Fourier plane from an already-
    !!          FFT'd volume at orientation e, via periodic-wrap KB interpolation.
    !!          No packed/Friedel storage on the plane side (note section 5).
    subroutine forward_plane( self, vol_img, e, plane )
        class(pcg_reconstruction), intent(in)    :: self
        class(image),               intent(in)    :: vol_img
        class(ori),                  intent(in)    :: e
        complex,                     intent(out)   :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                             &self%lims2(2,1):self%lims2(2,2))
        real    :: e_rotmat(3,3), loc(3)
        integer :: h, k
        plane    = cmplx(0.,0.)
        e_rotmat = e%get_mat()
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                loc = matmul(real([h,k,0]), e_rotmat)
                plane(h,k) = self%interp_from_volume(vol_img, loc)
            end do
        end do
    end subroutine forward_plane

    !>  \brief  F^dagger G_i^dagger, accumulate form: scatters a full (unpacked,
    !!          both-sign h) disk plane (matching forward_plane's domain) into a
    !!          full-range complex volume accumulator. Written fresh, not derived
    !!          from reconstructor%insert_plane_oversamp or compress_exp (note
    !!          section 5). This is the literal transpose of forward_plane's
    !!          gather: each plane entry is scattered once, at its own location,
    !!          with its own freshly evaluated KB window -- verified against
    !!          forward_plane by the adjoint dot-product test for arbitrary
    !!          orientations.
    subroutine adjoint_plane_add( self, plane, e, vol_accum )
        class(pcg_reconstruction), intent(in)    :: self
        complex,                    intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                            &self%lims2(2,1):self%lims2(2,2))
        class(ori),                  intent(in)    :: e
        complex,                     intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                 &self%lims3(2,1):self%lims3(2,2),&
                                                                 &self%lims3(3,1):self%lims3(3,2))
        real    :: e_rotmat(3,3), loc(3)
        integer :: h, k
        e_rotmat = e%get_mat()
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                if( plane(h,k) == cmplx(0.,0.) ) cycle
                loc = matmul(real([h,k,0]), e_rotmat)
                call scatter_one(self, loc, plane(h,k), vol_accum)
            end do
        end do
    end subroutine adjoint_plane_add

    !>  \brief  scatters one Fourier value at one location into the volume-space
    !!          accumulator using a freshly evaluated KB window (private helper
    !!          shared by adjoint_plane_add's direct and mirror scatters).
    subroutine scatter_one( self, loc, val, vol_accum )
        class(pcg_reconstruction), intent(in)    :: self
        real,                        intent(in)    :: loc(3)
        complex,                     intent(in)    :: val
        complex,                     intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                 &self%lims3(2,1):self%lims3(2,2),&
                                                                 &self%lims3(3,1):self%lims3(3,2))
        real    :: w(self%wdim,self%wdim,self%wdim)
        integer :: i0(3), di, dj, dk, hh, kk, mm
        i0 = nint(loc) - self%iwinsz
        call self%kbwin%apod_mat_3d(loc, self%iwinsz, self%wdim, w)
        do dk = 1, self%wdim
            mm = cyci_1d(self%wlims, i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = cyci_1d(self%wlims, i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = cyci_1d(self%wlims, i0(1)+di-1)
                    vol_accum(hh,kk,mm) = vol_accum(hh,kk,mm) + w(di,dj,dk) * val
                end do
            end do
        end do
    end subroutine scatter_one

    !>  \brief  H p = sum_i G_i^dagger G_i p + lambda * p, for a real trial volume p
    !!          and a set of projection orientations.
    function apply_normal( self, p, orientations ) result( hp )
        class(pcg_reconstruction), intent(inout) :: self
        real,                       intent(in)    :: p(self%box,self%box,self%box)
        class(oris),                 intent(inout) :: orientations
        real, allocatable :: hp(:,:,:)
        complex, allocatable :: vol_accum(:,:,:), plane(:,:)
        type(image) :: tmp
        type(ori)   :: e
        integer     :: i, nprojs
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        allocate(plane(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        call tmp%new([self%box,self%box,self%box], self%smpd)
        call tmp%set_rmat(p, .false.)
        call tmp%fft()
        nprojs = orientations%get_noris()
        call e%new(.false.)
        do i = 1, nprojs
            call orientations%get_ori(i, e)
            call self%forward_plane(tmp, e, plane)
            call self%adjoint_plane_add(plane, e, vol_accum)
        end do
        call tmp%kill
        hp = self%fold_and_ifft(vol_accum)
        hp = hp + self%lambda * p
    end function apply_normal

    !>  \brief  b = sum_i G_i^dagger y_i, the data right-hand side (no prior term).
    function apply_adjoint_all( self, y_planes, orientations ) result( b )
        class(pcg_reconstruction), intent(inout) :: self
        complex,                    intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                               &self%lims2(2,1):self%lims2(2,2), *)
        class(oris),                 intent(inout) :: orientations
        real, allocatable :: b(:,:,:)
        complex, allocatable :: vol_accum(:,:,:)
        type(ori) :: e
        integer   :: i, nprojs
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        nprojs = orientations%get_noris()
        call e%new(.false.)
        do i = 1, nprojs
            call orientations%get_ori(i, e)
            call self%adjoint_plane_add(y_planes(:,:,i), e, vol_accum)
        end do
        b = self%fold_and_ifft(vol_accum)
    end function apply_adjoint_all

    !>  \brief  deterministic double-precision real-volume dot product (note section 7).
    pure function dot_real_volume( self, a, b ) result( d )
        class(pcg_reconstruction), intent(in) :: self
        real,                       intent(in) :: a(self%box,self%box,self%box)
        real,                       intent(in) :: b(self%box,self%box,self%box)
        real(dp) :: d
        d = sum(real(a,dp) * real(b,dp))
    end function dot_real_volume

    ! PRIVATE HELPERS

    !>  \brief  periodic-wrap KB gather of one Fourier component directly from an
    !!          already-FFT'd volume's native (packed) storage, reading via
    !!          get_fcomp/comp_addr_phys (which already handle Friedel folding on
    !!          read) so no expanded/padded cache is required.
    function interp_from_volume( self, vol_img, loc ) result( comp )
        class(pcg_reconstruction), intent(in) :: self
        class(image),               intent(in) :: vol_img
        real,                        intent(in) :: loc(3)
        complex :: comp
        real    :: w(self%wdim,self%wdim,self%wdim)
        integer :: i0(3), di, dj, dk, hh, kk, mm, phys(3)
        i0 = nint(loc) - self%iwinsz
        call self%kbwin%apod_mat_3d(loc, self%iwinsz, self%wdim, w)
        comp = cmplx(0.,0.)
        do dk = 1, self%wdim
            mm = cyci_1d(self%wlims, i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = cyci_1d(self%wlims, i0(2)+dj-1)
                do di = 1, self%wdim
                    hh   = cyci_1d(self%wlims, i0(1)+di-1)
                    phys = vol_img%comp_addr_phys(hh,kk,mm)
                    comp = comp + w(di,dj,dk) * vol_img%get_fcomp([hh,kk,mm], phys)
                end do
            end do
        end do
    end function interp_from_volume

    !>  \brief  folds a full-range (both-sign h) complex volume accumulator into an
    !!          image's native packed storage (h>=0 half only, per the packed
    !!          real-FFT convention: only the h axis is truncated, k/m stay full
    !!          range) and inverse-FFTs it to a real volume. The h>=0-only fold is
    !!          valid because adjoint_plane_add's scatter is Hermitian-consistent
    !!          by construction: with the odd, symmetric KB window (see `new`),
    !!          negating loc exactly negates the interpolation window as a set, so
    !!          a plane entry at (h,k) and its Friedel mirror at (-h,-k) [both
    !!          visited independently, since lims2 is a full symmetric disk] land
    !!          on genuinely mirrored window footprints -- this is exactly what the
    !!          adjoint-dot-product test verifies for arbitrary orientations.
    !!
    !!          h = lims3(1,2) (the redundant Nyquist mate, present only because
    !!          lims3 spans the "including redundant Friedel mates" range) is
    !!          special: adjoint_plane_add's scatter wraps via wlims, which never
    !!          produces +lims3(1,2) (only its canonical negative representative),
    !!          so that slot of vol_accum is always zero. The Nyquist bin's actual
    !!          accumulated value lives at h = -lims3(1,2) instead; read from there
    !!          for that one h value.
    function fold_and_ifft( self, vol_accum ) result( z )
        class(pcg_reconstruction), intent(in) :: self
        complex,                    intent(in) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                             &self%lims3(2,1):self%lims3(2,2),&
                                                             &self%lims3(3,1):self%lims3(3,2))
        real, allocatable :: z(:,:,:)
        type(image) :: zimg
        integer     :: h, hh, k, m, phys(3)
        call zimg%new([self%box,self%box,self%box], self%smpd)
        call zimg%zero_and_flag_ft()
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    hh   = cyci_1d(self%wlims, h)
                    phys = zimg%comp_addr_phys(h,k,m)
                    call zimg%set_cmat_at(phys(1),phys(2),phys(3), vol_accum(hh,k,m))
                end do
            end do
        end do
        call zimg%ifft()
        z = zimg%get_rmat()
        call zimg%kill
    end function fold_and_ifft

    ! SOLVER

    !>  \brief  plain (unpreconditioned, M=I) CG solve of H x = b, per note section 7.
    subroutine solve( self, y_planes, orientations, x, maxits, rtol, rel_res_hist, niters )
        class(pcg_reconstruction), intent(inout) :: self
        complex,                    intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                               &self%lims2(2,1):self%lims2(2,2), *)
        class(oris),                 intent(inout) :: orientations
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        real, allocatable :: b(:,:,:), r(:,:,:), p(:,:,:), hp(:,:,:), hist(:)
        real(dp) :: rho, rho_new, rho0, alpha, beta, pHp
        integer  :: mmaxits, iter, n_done
        real     :: rrtol
        mmaxits = 50
        if( present(maxits) ) mmaxits = maxits
        rrtol = 1.0e-4
        if( present(rtol) ) rrtol = rtol
        allocate(hist(mmaxits))
        b = self%apply_adjoint_all(y_planes, orientations)
        hp = self%apply_normal(x, orientations)
        r  = b - hp
        p  = r
        rho  = self%dot_real_volume(r,r)
        rho0 = rho
        n_done = 0
        do iter = 1, mmaxits
            hp  = self%apply_normal(p, orientations)
            pHp = self%dot_real_volume(p,hp)
            if( pHp <= 0.0_dp ) THROW_HARD('non-positive dot(p,Hp); PCG lost positive-definiteness; solve')
            if( pHp /= pHp )    THROW_HARD('non-finite dot(p,Hp); solve')
            alpha = rho / pHp
            x  = x + real(alpha) * p
            r  = r - real(alpha) * hp
            rho_new = self%dot_real_volume(r,r)
            n_done  = iter
            hist(iter) = real(sqrt(rho_new/rho0))
            if( sqrt(rho_new/rho0) <= real(rrtol,dp) ) exit
            beta = rho_new / rho
            p    = r + real(beta) * p
            rho  = rho_new
        end do
        if( present(niters) ) niters = n_done
        if( present(rel_res_hist) ) allocate(rel_res_hist(n_done), source=hist(1:n_done))
    end subroutine solve

end module simple_pcg_reconstruction
