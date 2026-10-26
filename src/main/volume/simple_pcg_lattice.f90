!@descr: padded-lattice geometry, solve support and KB plane gather/scatter (flat or Ewald paraboloid) of the PCG and wave operators
! (reconstructor_pcg and flex_pcg_t extend it; the wave operators sample the Ewald paraboloid with it)
module simple_pcg_lattice
use simple_core_module_api,   only: cyci_1d, fdim, KBALPHA, KBWINSZ, kbinterpol, OSMPL_PAD_FAC, simple_exception
use simple_image,             only: image
use simple_ftiter,            only: ftiter
use simple_gridding,          only: kb_stencil_envelope_1d
use simple_cartesian_fourier, only: gather_packed_window
use, intrinsic :: iso_c_binding, only: c_float, c_float_complex
!$ use omp_lib,               only: omp_in_parallel
implicit none
private
#include "simple_local_flags.inc"

public :: pcg_lattice
public :: plane_loc, rim_free, win_wraps, win_below_fold, scatter_window, scatter_window_real, scatter_window_pair
public :: scatter_window_nowrap, scatter_window_pair_nowrap, scatter_window_cmplx_nowrap

!> A solve whose unknown lives on the native box while every Fourier operation runs on the padf-times
!! padded lattice (centre-pad in, centre-crop out): the KB window, the OpenMP colouring stride of the
!! scatter, the array and period-box wrap limits, the radius below which a window cannot wrap, the wrap
!! table, the deposition envelope, and the solve support (soft window, hard mask window > 0). The PCG
!! operators extend it. A plane is the full both-sign-h disk of native frequencies (h,k), h^2+k^2 <= R^2,
!! R = box/2; its point (h,k) sits at padf*[h,k,kz] R_i on the padded lattice, kz = curv*(h^2+k^2): the
!! central section for curv = 0, the Ewald paraboloid of a wave model otherwise (ewald_curv). The gather
!! reads a plane through the KB window, the scatter is its literal transpose, coloured over h for threads
!! with a serial pass over the windows that wrap, so it is reproducible at any thread count.
type :: pcg_lattice
    integer          :: box        = 0   !< native box: the solver unknown lives here
    integer          :: boxpd      = 0   !< padf*box: every Fourier operation happens here
    integer          :: padf       = 1   !< oversampling factor, OSMPL_PAD_FAC
    integer          :: pad_off    = 0   !< centred pad/crop offset, (boxpd-box)/2
    integer          :: Rnat       = 0   !< native Nyquist radius, box/2
    integer          :: lims2(2,2) = 0   !< (h/k, lo/hi) plane bounds, full symmetric disk
    integer          :: sqlp       = 0   !< squared native Nyquist radius
    real             :: padsc      = 1.0 !< padf**3, undoes fwd_ft's 1/product(ldim)
    real             :: smpd       = 1.0
    type(kbinterpol) :: kbwin
    integer          :: iwinsz     = 0
    integer          :: wdim       = 0
    integer          :: stride     = 0   !< OpenMP colouring stride of the scatter
    integer          :: lims3(3,2) = 0   !< (h/k/m, lo/hi) volume array bounds of the padded lattice
    integer          :: cdim(3)    = 0   !< complex array shape of the padded lattice
    integer          :: wlims(2)   = 0   !< [lo,hi] canonical period-box wrap range
    integer          :: sq_rim     = 0   !< below this h^2+k^2 a KB window cannot wrap
    integer, allocatable :: wrap(:)      !< cyci_1d over every index a KB window can reach
    real,    allocatable :: dep1d(:)     !< 1D deposition envelope of the padded lattice
    real,    allocatable :: mask(:,:,:)  !< hard solve support P (window > 0)
    real,    allocatable :: window(:,:,:)!< soft output window; the shipped map is window*u
    logical              :: l_mask = .false.
  contains
    procedure :: new_lattice
    procedure :: set_window_sphere
    procedure :: set_window_volume
    procedure :: mask_mul
    procedure :: native_deposition_envelope
    procedure :: divide_deposition_plane
    procedure :: ewald_curv
    procedure :: gather_plane
    procedure :: scatter_plane
    procedure :: gather_cube
    procedure :: scatter_cube
    procedure :: fold_lattice_accum
    procedure :: kill_lattice
end type pcg_lattice

contains

    !> The geometry of a solve at box: KB window, colouring stride, padded box, limits, rim radius,
    !! wrap table and deposition envelope. winsz is the KB window's half-width (KBWINSZ, three taps,
    !! unless a caller needs a more accurate gather)
    subroutine new_lattice( self, box, smpd, winsz )
        class(pcg_lattice), intent(inout) :: self
        integer,            intent(in)    :: box
        real,               intent(in)    :: smpd
        real,    optional,  intent(in)    :: winsz
        type(ftiter) :: fit
        integer      :: lo, hi, i
        real         :: rlim
        call self%kill_lattice
        self%box    = box
        self%smpd   = smpd
        if( present(winsz) )then
            self%kbwin = kbinterpol(winsz, KBALPHA)
        else
            self%kbwin = kbinterpol(KBWINSZ, KBALPHA)
        endif
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        ! odd width centred on nint(loc): negating loc negates the window as a set,
        ! which the h>=0 folds rely on
        self%wdim   = 2*self%iwinsz + 1
        ! oversampling: the unknown lives on the native box, every Fourier operation
        ! on the padf-times padded lattice, as in Fourier gridding
        self%padf    = OSMPL_PAD_FAC
        ! same-colour h-lines stay a full window apart along one axis after rotation: plane points
        ! padf*stride apart in Euclidean norm are at least that / sqrt(3) apart in max-norm
        self%stride  = ceiling(sqrt(3.0) * real(self%wdim) / real(self%padf))
        self%boxpd   = self%padf * box
        self%pad_off = (self%boxpd - box) / 2
        self%Rnat    = box / 2
        ! full symmetric (both-sign h) native disk: gather_plane/scatter_plane are then an exact
        ! adjoint pair for any orientation; fwd_ft divides by product(ldim), padsc restores the
        ! native scale at coincident frequencies
        self%lims2(1,:) = [-self%Rnat, self%Rnat]
        self%lims2(2,:) = [-self%Rnat, self%Rnat]
        self%sqlp       = self%Rnat * self%Rnat
        self%padsc      = real(self%padf)**3
        ! the padded image's logical limits and packed shape, without the image
        call fit%new([self%boxpd,self%boxpd,self%boxpd], smpd)
        self%lims3 = fit%loop_lims(3)
        self%cdim  = [fdim(self%boxpd), self%boxpd, self%boxpd]
        ! true period-box wrap range: lims3(1,:) spans both Friedel Nyquist mates
        ! (one longer than the period), axes 2/3 do not
        self%wlims = self%lims3(2,:)
        ! squared plane radius below which a KB window cannot reach the wrap boundary;
        ! |loc| = padf*sqrt(h^2+k^2) is rotation-independent, so (h,k) decides (conservative)
        rlim = real(min(self%wlims(2) - self%iwinsz, -self%wlims(1) - self%iwinsz)) - 0.5
        self%sq_rim = max(0, int((rlim / real(self%padf))**2) - 1)
        lo = self%wlims(1) - self%iwinsz - 1
        hi = self%wlims(2) + self%iwinsz + 1
        allocate(self%wrap(lo:hi))
        do i = lo, hi
            self%wrap(i) = cyci_1d(self%wlims, i)
        end do
        ! scattering through the KB window multiplies the real-space kernel by the window's transform
        call kb_stencil_envelope_1d(self%kbwin, self%boxpd, self%dep1d)
    end subroutine new_lattice

    !> Spherical support as a CONSTRAINT ON THE SOLVE: the window is image%mask3D_soft of radius mskrad
    !! (pixels) on a unit volume (backgr=0.), the soft mask the gridding restoration applies; the solve runs
    !! on the domain window > 0 and the shipped map is window*u. It removes the solvent, where
    !! deapodization amplifies hardest, and shrinks the problem. No consumer masks a second time.
    !! mskrad <= 0 clears the support.
    subroutine set_window_sphere( self, mskrad )
        class(pcg_lattice), intent(inout) :: self
        real,               intent(in)    :: mskrad
        type(image) :: mimg
        real, allocatable :: ones(:,:,:)
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        self%l_mask = .false.
        if( mskrad <= 0.0 ) return
        allocate(ones(self%box,self%box,self%box), source=1.0)
        call mimg%new([self%box,self%box,self%box], self%smpd)
        call mimg%set_rmat(ones, .false.)
        call mimg%mask3D_soft(mskrad, backgr=0.)
        self%window = mimg%get_rmat()
        call mimg%kill
        deallocate(ones)
        call install_support(self)
    end subroutine set_window_sphere

    !> A caller-supplied real-space [0,1] window on the native box (clipped) and its hard support
    subroutine set_window_volume( self, mskvol )
        class(pcg_lattice), intent(inout) :: self
        class(image),       intent(in)    :: mskvol
        integer :: mdim(3)
        mdim = mskvol%get_ldim()
        if( any(mdim /= self%box) ) THROW_HARD('support window dimensions differ from the solve box; set_window_volume')
        if( mskvol%is_ft() ) THROW_HARD('support window must be in real space; set_window_volume')
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        self%window = mskvol%get_rmat()
        self%window = min(1.0, max(0.0, self%window))
        if( .not. any(self%window > 0.0) ) THROW_HARD('support window is empty; set_window_volume')
        call install_support(self)
    end subroutine set_window_volume

    !> The solve domain is window > 0: P^2 = P, so the projections in operator, right-hand side and
    !! preconditioner are exact, u is the estimate on that domain and the shipped map window*u is one
    !! estimate times one soft window, exactly what the gridding restoration ships
    subroutine install_support( self )
        class(pcg_lattice), intent(inout) :: self
        if( allocated(self%mask) ) deallocate(self%mask)
        allocate(self%mask(self%box,self%box,self%box), source=merge(1.0, 0.0, self%window > 0.0))
        self%l_mask = .true.
    end subroutine install_support

    !> v = P v on the native box (identity without a support)
    pure subroutine mask_mul( self, v )
        class(pcg_lattice), intent(in)    :: self
        real,               intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_mask ) return
        v = v * self%mask
    end subroutine mask_mul

    !> The deposition envelope of the padded lattice cropped to the native box, unnormalized
    function native_deposition_envelope( self ) result( env )
        class(pcg_lattice), intent(in) :: self
        real, allocatable :: env(:,:,:)
        integer :: i, j, k, o
        o = self%pad_off
        allocate(env(self%box,self%box,self%box))
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    env(i,j,k) = self%dep1d(o+i)*self%dep1d(o+j)*self%dep1d(o+k)
                end do
            end do
        end do
        !$omp end parallel do
    end function native_deposition_envelope

    !> Divide the deposition envelope out of plane k of a real-space kernel on the padded lattice
    !! (zero where the envelope vanishes); r is 1-based and at least boxpd on every axis (an image
    !! buffer with its FFT padding rows is fine; only 1..boxpd is touched). Serial, so callers thread
    !! over planes or kernels.
    subroutine divide_deposition_plane( self, r, k )
        class(pcg_lattice), intent(in)    :: self
        real,               intent(inout) :: r(:,:,:)
        integer,            intent(in)    :: k
        real, parameter :: EPS_D = 1.0e-8
        real    :: depval
        integer :: i, j
        do j = 1, self%boxpd
            do i = 1, self%boxpd
                depval = self%dep1d(i) * self%dep1d(j) * self%dep1d(k)
                if( abs(depval) > EPS_D )then
                    r(i,j,k) = r(i,j,k) / depval
                else
                    r(i,j,k) = 0.0
                endif
            end do
        end do
    end subroutine divide_deposition_plane

    !> The curvature coefficient of the Ewald paraboloid, kz = curv*(h^2+k^2) in native Fourier voxels:
    !! kappa*wl/(2*box), kappa = -1, +1, 0 for ewald=neg, pos, no and wl the wavelength in pixels
    !! (ctfvars%wl). Refused when the paraboloid would carry a window past the wrap table: the padded
    !! point must leave the Nyquist sphere by less than half a voxel, so that it rounds onto the lattice's
    !! Nyquist index at most. gather_plane and scatter_plane take curv from here.
    real function ewald_curv( self, kappa, wl ) result( curv )
        class(pcg_lattice), intent(in) :: self
        integer,            intent(in) :: kappa
        real,               intent(in) :: wl
        real :: kzmax, excess
        if( abs(kappa) > 1 ) THROW_HARD('kappa must be -1, 0 or 1; ewald_curv')
        curv   = real(kappa) * wl / (2.0 * real(self%box))
        kzmax  = abs(curv) * real(self%sqlp)
        excess = real(self%padf) * (sqrt(real(self%sqlp) + kzmax**2) - sqrt(real(self%sqlp)))
        ! half a voxel less a margin for the single-precision rounding of loc
        if( excess >= 0.49 ) THROW_HARD('the Ewald paraboloid leaves the padded lattice; ewald_curv')
    end function ewald_curv

    !> The padded-lattice location of plane point (h,k) under rotation rot: padf*[h,k,kz] rot,
    !! kz = curv*(h^2+k^2); curv = 0 is the central section, evaluated as it always was. Module-level,
    !! like the window helpers below: no dispatch in the ~1e8 calls of an accumulation
    pure function plane_loc( self, h, k, rot, curv ) result( loc )
        class(pcg_lattice), intent(in) :: self
        integer,            intent(in) :: h, k
        real,               intent(in) :: rot(3,3), curv
        real :: loc(3)
        if( curv == 0.0 )then
            loc = real(self%padf) * matmul(real([h,k,0]), rot)
        else
            loc = real(self%padf) * matmul([real(h), real(k), curv * real(h*h + k*k)], rot)
        endif
    end function plane_loc

    !> Can the KB window of plane point (h,k) provably not wrap? |loc| = padf*sqrt(h^2+k^2+kz^2) for
    !! every rotation, so the plane point alone decides (conservative, see sq_rim)
    pure logical function rim_free( self, h, k, curv )
        class(pcg_lattice), intent(in) :: self
        integer,            intent(in) :: h, k
        real,               intent(in) :: curv
        real :: kz
        if( curv == 0.0 )then
            rim_free = h*h + k*k <= self%sq_rim
        else
            kz       = curv * real(h*h + k*k)
            rim_free = real(h*h + k*k) + kz*kz <= real(self%sq_rim)
        endif
    end function rim_free

    !> G_i: gathers the plane of rotation rot (and curvature curv) from cmat, the padded lattice's packed
    !! transform; padsc restores the native scale. Threaded over the plane
    subroutine gather_plane( self, cmat, rot, curv, plane )
        class(pcg_lattice), intent(in)  :: self
        complex,            intent(in)  :: cmat(:,:,:)
        real,               intent(in)  :: rot(3,3), curv
        complex,            intent(out) :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                &self%lims2(2,1):self%lims2(2,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, i0(3)
        plane = cmplx(0.,0.)
        ! threaded over the plane unless a caller's particle loop already is
        !$omp parallel do collapse(2) default(shared) private(h,k,loc,i0,w) &
        !$omp schedule(static) proc_bind(close) if(.not. omp_in_parallel())
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                loc = plane_loc(self, h, k, rot, curv)
                i0  = nint(loc) - self%iwinsz
                call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                plane(h,k) = self%padsc*gather_packed_window(cmat,lbound(self%wrap,1),self%wrap,i0,w)
            end do
        end do
        !$omp end parallel do
    end subroutine gather_plane

    !> G_i^dagger, accumulate form: the literal transpose of gather_plane, replicated over the rotations
    !! rots(:,:,g) (R_i S_g for a point group), into the full-range accumulator. h-strided colouring keeps
    !! the windows of one sweep disjoint (rows a stride apart stay a window apart under any rotation and
    !! curvature, since |loc_1 - loc_2| >= padf |h_1 - h_2|), the windows that wrap go in a serial pass
    subroutine scatter_plane( self, plane, rots, curv, vol_accum )
        class(pcg_lattice), intent(in)    :: self
        complex,            intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                  &self%lims2(2,1):self%lims2(2,2))
        real,               intent(in)    :: rots(:,:,:), curv
        complex,            intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                      &self%lims3(2,1):self%lims3(2,2),&
                                                      &self%lims3(3,1):self%lims3(3,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, l, g, i0(3)
        !$omp parallel default(shared) private(h,k,l,g,loc,i0,w) proc_bind(close)
        do g = 1, size(rots,3)
            do l = 0, self%stride-1
                !$omp do schedule(static,1)
                do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp ) cycle
                        if( plane(h,k) == cmplx(0.,0.) ) cycle
                        loc = plane_loc(self, h, k, rots(:,:,g), curv)
                        i0  = nint(loc) - self%iwinsz
                        if( win_wraps(self, i0) ) cycle   ! rim: deferred to the serial pass
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_window_cmplx_nowrap(self, i0, w, self%padsc * plane(h,k), vol_accum)
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        call scatter_rim(self, plane, rots, curv, vol_accum)
    end subroutine scatter_plane

    !> the serial pass of a plane over the wrapping rim (see win_wraps): confined to the outermost shell,
    !! reproducible unlike an atomic; each row visits only its two segments outside the free middle |k| < kb
    subroutine scatter_rim( self, plane, rots, curv, vol_accum )
        class(pcg_lattice), intent(in)    :: self
        complex,            intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                  &self%lims2(2,1):self%lims2(2,2))
        real,               intent(in)    :: rots(:,:,:), curv
        complex,            intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                      &self%lims3(2,1):self%lims3(2,2),&
                                                      &self%lims3(3,1):self%lims3(3,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, g, i0(3), kb, iseg, klo, khi
        do g = 1, size(rots,3)
            do h = self%lims2(1,1), self%lims2(1,2)
                kb = rim_row_start(self, h, curv)
                do iseg = 1, 2
                    if( iseg == 1 )then
                        klo = self%lims2(2,1)
                        khi = min(-kb, self%lims2(2,2))
                    else
                        klo = max(kb, -kb + 1, self%lims2(2,1))
                        khi = self%lims2(2,2)
                    endif
                    do k = klo, khi
                        if( h*h + k*k > self%sqlp ) cycle
                        if( rim_free(self, h, k, curv) ) cycle
                        if( plane(h,k) == cmplx(0.,0.) ) cycle
                        loc = plane_loc(self, h, k, rots(:,:,g), curv)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. win_wraps(self, i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_window(self, i0, w, self%padsc * plane(h,k), vol_accum)
                    end do
                end do
            end do
        end do
    end subroutine scatter_rim

    !> the first k >= 0 of row h that rim_free does not exempt: |k| < kb is free
    pure integer function rim_row_start( self, h, curv ) result( kb )
        class(pcg_lattice), intent(in) :: self
        integer,            intent(in) :: h
        real,               intent(in) :: curv
        ! rim_free needs h^2 + k^2 <= sq_rim, the curvature only shrinks the free disk
        kb = isqrt(self%sq_rim - h*h)
        do while( kb >= 0 )
            if( rim_free(self, h, kb, curv) ) exit
            kb = kb - 1
        end do
        kb = kb + 1
    end function rim_row_start

    !> floor(sqrt(n)), -1 for n < 0
    pure integer function isqrt( n ) result( r )
        integer, intent(in) :: n
        if( n < 0 )then
            r = -1
            return
        endif
        r = int(sqrt(real(n)))
        do while( (r + 1) * (r + 1) <= n )
            r = r + 1
        end do
        do while( r * r > n )
            r = r - 1
        end do
    end function isqrt

    !> The native spectrum of the volume rotated by rot, the cube cube(h,k,l) = padsc * gather(padf [h,k,l] rot)
    !! on the native sphere h^2+k^2+l^2 <= R^2 with |l| < R (the plane disk lims2 times l = -(R-1):(R-1); a
    !! z sampled at N points has one Nyquist value, which the sphere holds only at the origin of the plane, so
    !! it is left out), zero elsewhere: the rotated cube of multislice, the planes l = 0 its central section.
    !! Threaded over the cube unless a caller's particle loop already is
    subroutine gather_cube( self, cmat, rot, cube )
        class(pcg_lattice), intent(in)  :: self
        complex,            intent(in)  :: cmat(:,:,:)
        real,               intent(in)  :: rot(3,3)
        complex,            intent(out) :: cube(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),&
                                               &-(self%lims2(1,2)-1):self%lims2(1,2)-1)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, l, i0(3), rl
        rl   = self%lims2(1,2) - 1
        cube = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,l,loc,i0,w) &
        !$omp schedule(static) proc_bind(close) if(.not. omp_in_parallel())
        do l = -rl, rl
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    if( h*h + k*k + l*l > self%sqlp ) cycle
                    loc = real(self%padf) * matmul(real([h,k,l]), rot)
                    i0  = nint(loc) - self%iwinsz
                    call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                    cube(h,k,l) = self%padsc*gather_packed_window(cmat,lbound(self%wrap,1),self%wrap,i0,w)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine gather_cube

    !> The literal transpose of gather_cube, replicated over rots(:,:,g), accumulated into the full-range
    !! accumulator: the coloured scatter of scatter_plane (points of different h stay padf |h_1 - h_2| apart
    !! under any rotation), the windows that wrap in a serial pass
    subroutine scatter_cube( self, cube, rots, vol_accum )
        class(pcg_lattice), intent(in)    :: self
        complex,            intent(in)    :: cube(self%lims2(1,1):self%lims2(1,2),self%lims2(2,1):self%lims2(2,2),&
                                                 &-(self%lims2(1,2)-1):self%lims2(1,2)-1)
        real,               intent(in)    :: rots(:,:,:)
        complex,            intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                      &self%lims3(2,1):self%lims3(2,2),&
                                                      &self%lims3(3,1):self%lims3(3,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, l, m, g, i0(3), rl, hb, iseg, hlo, hhi
        rl = self%lims2(1,2) - 1
        !$omp parallel default(shared) private(h,k,l,m,g,loc,i0,w) proc_bind(close)
        do g = 1, size(rots,3)
            do m = 0, self%stride-1
                !$omp do schedule(static,1)
                do h = self%lims2(1,1)+m, self%lims2(1,2), self%stride
                    do l = -rl, rl
                        do k = self%lims2(2,1), self%lims2(2,2)
                            if( h*h + k*k + l*l > self%sqlp ) cycle
                            if( cube(h,k,l) == cmplx(0.,0.) ) cycle
                            loc = real(self%padf) * matmul(real([h,k,l]), rots(:,:,g))
                            i0  = nint(loc) - self%iwinsz
                            if( win_wraps(self, i0) ) cycle   ! rim: deferred to the serial pass
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_window_cmplx_nowrap(self, i0, w, self%padsc * cube(h,k,l), vol_accum)
                        end do
                    end do
                end do
                !$omp end do
            end do
        end do
        !$omp end parallel
        ! the wrapping rim, serial (|loc| = padf |[h,k,l]| decides, as rim_free does for a flat plane);
        ! each row visits only its two segments outside the free middle |h| <= hb
        do g = 1, size(rots,3)
            do l = -rl, rl
                do k = self%lims2(2,1), self%lims2(2,2)
                    hb = isqrt(self%sq_rim - k*k - l*l)
                    do iseg = 1, 2
                        if( iseg == 1 )then
                            hlo = self%lims2(1,1)
                            hhi = min(-hb - 1, self%lims2(1,2))
                        else
                            hlo = max(hb + 1, -hb, self%lims2(1,1))
                            hhi = self%lims2(1,2)
                        endif
                        do h = hlo, hhi
                            if( h*h + k*k + l*l > self%sqlp ) cycle
                            if( h*h + k*k + l*l <= self%sq_rim ) cycle
                            if( cube(h,k,l) == cmplx(0.,0.) ) cycle
                            loc = real(self%padf) * matmul(real([h,k,l]), rots(:,:,g))
                            i0  = nint(loc) - self%iwinsz
                            if( .not. win_wraps(self, i0) ) cycle
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_window(self, i0, w, self%padsc * cube(h,k,l), vol_accum)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine scatter_cube

    !> The real volume of a full-range accumulator: laid on the half-complex lattice of the work image
    !! wimg (boxpd^3, the caller's, so that its plans persist), inverse transformed and cropped to the
    !! native box. l_hermitian = .false. takes the h >= 0 half as it is, exact for the accumulator of
    !! Hermitian planes (the PCG backend's, the central section of real images); .true. lays
    !! (A(q) + conj(A(-q)))/2 instead, which makes the result the real part of the full complex inverse
    !! transform, as the two-cap planes of a curved wave model need. A caller that extends the lattice
    !! passes its parent component as self and its own work image, disjoint subobjects
    function fold_lattice_accum( self, vol_accum, wimg, l_hermitian ) result( z )
        class(pcg_lattice), intent(in)    :: self
        complex,            intent(in)    :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                      &self%lims3(2,1):self%lims3(2,2),&
                                                      &self%lims3(3,1):self%lims3(3,2))
        class(image),       intent(inout) :: wimg
        logical,            intent(in)    :: l_hermitian
        real, allocatable :: z(:,:,:)
        complex(c_float_complex), pointer :: cptr(:,:,:)
        real(c_float),            pointer :: rptr(:,:,:)
        integer :: h, hh, k, m, pk, pm, n2, n3, o, j
        ! every packed element (h >= 0, all k and m) is written in place through the image's pointer
        ! (get_cmat_ptr states the contract): a set_cmat_at call per element and the copies of
        ! get_rmat/center_crop were most of the fold; physical addresses as comp_addr_phys
        call wimg%get_cmat_ptr(cptr)
        call wimg%set_ft(.true.)
        n2 = self%boxpd
        n3 = self%boxpd
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,pk,pm) &
        !$omp schedule(static) proc_bind(close)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                pm = m + 1; if( m < 0 ) pm = pm + n3
                pk = k + 1; if( k < 0 ) pk = pk + n2
                if( l_hermitian )then
                    do h = 0, self%lims3(1,2)
                        hh = self%wrap(h)
                        cptr(h+1,pk,pm) = 0.5 * (vol_accum(hh,k,m) + &
                            &conjg(vol_accum(self%wrap(-h),self%wrap(-k),self%wrap(-m))))
                    end do
                else
                    do h = 0, self%lims3(1,2)
                        cptr(h+1,pk,pm) = vol_accum(self%wrap(h),k,m)
                    end do
                endif
            end do
        end do
        !$omp end parallel do
        call wimg%ifft()
        ! back to the native lattice: the unknown never leaves box^3
        call wimg%get_rmat_ptr(rptr)
        o = self%pad_off
        allocate(z(self%box,self%box,self%box))
        !$omp parallel do collapse(2) default(shared) private(j,k) schedule(static) proc_bind(close)
        do k = 1, self%box
            do j = 1, self%box
                z(:,j,k) = rptr(o+1:o+self%box,o+j,o+k)
            end do
        end do
        !$omp end parallel do
    end function fold_lattice_accum

    subroutine kill_lattice( self )
        class(pcg_lattice), intent(inout) :: self
        if( allocated(self%wrap)   ) deallocate(self%wrap)
        if( allocated(self%dep1d)  ) deallocate(self%dep1d)
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        self%l_mask = .false.
        self%box = 0; self%boxpd = 0; self%padf = 1; self%pad_off = 0; self%Rnat = 0
        self%lims2 = 0; self%sqlp = 0; self%padsc = 1.0
        self%iwinsz = 0; self%wdim = 0; self%stride = 0; self%lims3 = 0; self%cdim = 0
        self%wlims = 0; self%sq_rim = 0
    end subroutine kill_lattice

    !> does the window straddle the periodic wrap boundary? The h-strided colouring
    !! separates same-colour h-lines by more than a window only in unwrapped
    !! coordinates; after folding two such windows can overlap at the Nyquist rim
    !! and the sum becomes thread-order dependent. Callers scatter the interior in
    !! the coloured parallel pass and the wrapping rim serially. Module-level, not
    !! type-bound (no dispatch, ~1e8 calls per accumulation), like the scatters
    pure logical function win_wraps( self, i0 )
        class(pcg_lattice),       intent(in) :: self
        integer,                  intent(in) :: i0(3)
        win_wraps = any(i0 < self%wlims(1)) .or. any(i0 + self%wdim - 1 > self%wlims(2))
    end function win_wraps

    !> does the window lie wholly in the h < 0 half, clear of the wrapped Nyquist mate at wlims(1)?
    !! The accumulator of Hermitian planes is read as its h >= 0 half (fold_lattice_accum with
    !! l_hermitian = .false.), so such a window deposits nothing that is ever read: the Friedel
    !! mate of its plane point deposits the conjugate at the mirrored window. Not for the
    !! two-cap planes of a wave model, whose fold reads both halves
    pure logical function win_below_fold( self, i0 )
        class(pcg_lattice),       intent(in) :: self
        integer,                  intent(in) :: i0(3)
        win_below_fold = i0(1) + self%wdim - 1 < 0 .and. i0(1) > self%wlims(1)
    end function win_below_fold

    !> KB-weighted scatter of one complex value into the full-range accumulator
    pure subroutine scatter_window( self, i0, w, val, vol_accum )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,                  intent(in)    :: val
        complex,                  intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, hh, kk, mm
        do dk = 1, self%wdim
            mm = self%wrap(i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = self%wrap(i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = self%wrap(i0(1)+di-1)
                    vol_accum(hh,kk,mm) = vol_accum(hh,kk,mm) + w(di,dj,dk) * val
                end do
            end do
        end do
    end subroutine scatter_window

    !> real-valued counterpart for the density accumulator
    pure subroutine scatter_window_real( self, i0, w, val, vol_accum )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim), val
        real,                     intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, hh, kk, mm
        do dk = 1, self%wdim
            mm = self%wrap(i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = self%wrap(i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = self%wrap(i0(1)+di-1)
                    vol_accum(hh,kk,mm) = vol_accum(hh,kk,mm) + w(di,dj,dk) * val
                end do
            end do
        end do
    end subroutine scatter_window_real

    !> B and D through one wrapped KB-window traversal
    pure subroutine scatter_window_pair( self, i0, w, bval, dval, bacc, dacc )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim), dval
        complex,                  intent(in)    :: bval
        complex,                  intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                       &self%lims3(2,1):self%lims3(2,2),&
                                                       &self%lims3(3,1):self%lims3(3,2))
        real,                     intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
                                                       &self%lims3(2,1):self%lims3(2,2),&
                                                       &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, hh, kk, mm
        do dk = 1, self%wdim
            mm = self%wrap(i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = self%wrap(i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = self%wrap(i0(1)+di-1)
                    bacc(hh,kk,mm) = bacc(hh,kk,mm) + w(di,dj,dk) * bval
                    dacc(hh,kk,mm) = dacc(hh,kk,mm) + w(di,dj,dk) * dval
                end do
            end do
        end do
    end subroutine scatter_window_pair

    !> interior-only scatter for windows proven not to wrap: self%wrap is the
    !! identity there, so direct indexing is bit-identical and the inner run contiguous
    pure subroutine scatter_window_nowrap( self, i0, w, val, vol_accum )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        real,                     intent(in)    :: val
        real,                     intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, h0, kk, mm
        h0 = i0(1)
        do dk = 1, self%wdim
            mm = i0(3) + dk - 1
            do dj = 1, self%wdim
                kk = i0(2) + dj - 1
                do di = 1, self%wdim
                    vol_accum(h0+di-1,kk,mm) = vol_accum(h0+di-1,kk,mm) + w(di,dj,dk) * val
                end do
            end do
        end do
    end subroutine scatter_window_nowrap

    !> B and D through one non-wrapping, contiguous KB-window traversal
    pure subroutine scatter_window_pair_nowrap( self, i0, w, bval, dval, bacc, dacc )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim), dval
        complex,                  intent(in)    :: bval
        complex,                  intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                       &self%lims3(2,1):self%lims3(2,2),&
                                                       &self%lims3(3,1):self%lims3(3,2))
        real,                     intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
                                                       &self%lims3(2,1):self%lims3(2,2),&
                                                       &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, h0, kk, mm
        h0 = i0(1)
        do dk = 1, self%wdim
            mm = i0(3) + dk - 1
            do dj = 1, self%wdim
                kk = i0(2) + dj - 1
                do di = 1, self%wdim
                    bacc(h0+di-1,kk,mm) = bacc(h0+di-1,kk,mm) + w(di,dj,dk) * bval
                    dacc(h0+di-1,kk,mm) = dacc(h0+di-1,kk,mm) + w(di,dj,dk) * dval
                end do
            end do
        end do
    end subroutine scatter_window_pair_nowrap

    !> complex counterpart of scatter_window_nowrap, for the RHS scatter
    pure subroutine scatter_window_cmplx_nowrap( self, i0, w, val, vol_accum )
        class(pcg_lattice),       intent(in)    :: self
        integer,                  intent(in)    :: i0(3)
        real,                     intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,                  intent(in)    :: val
        complex,                  intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        integer :: di, dj, dk, h0, kk, mm
        h0 = i0(1)
        do dk = 1, self%wdim
            mm = i0(3) + dk - 1
            do dj = 1, self%wdim
                kk = i0(2) + dj - 1
                do di = 1, self%wdim
                    vol_accum(h0+di-1,kk,mm) = vol_accum(h0+di-1,kk,mm) + w(di,dj,dk) * val
                end do
            end do
        end do
    end subroutine scatter_window_cmplx_nowrap

end module simple_pcg_lattice
