!@descr: padded-lattice geometry and solve support shared by the PCG reconstructors (reconstructor_pcg, flex_pcg_t)
module simple_pcg_lattice
use simple_core_module_api, only: cyci_1d, KBALPHA, KBWINSZ, kbinterpol, OSMPL_PAD_FAC, simple_exception
use simple_image,           only: image
use simple_gridding,        only: kb_stencil_envelope_1d
implicit none
private
#include "simple_local_flags.inc"

public :: pcg_lattice

!> A solve whose unknown lives on the native box while every Fourier operation runs on the padf-times
!! padded lattice (centre-pad in, centre-crop out): the KB window, the OpenMP colouring stride of the
!! scatter, the array and period-box wrap limits, the radius below which a window cannot wrap, the wrap
!! table, the deposition envelope, and the solve support (soft window, hard mask window > 0). The PCG
!! operators extend it.
type :: pcg_lattice
    integer          :: box        = 0   !< native box: the solver unknown lives here
    integer          :: boxpd      = 0   !< padf*box: every Fourier operation happens here
    integer          :: padf       = 1   !< oversampling factor, OSMPL_PAD_FAC
    integer          :: pad_off    = 0   !< centred pad/crop offset, (boxpd-box)/2
    integer          :: Rnat       = 0   !< native Nyquist radius, box/2
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
    procedure :: kill_lattice
end type pcg_lattice

contains

    !> The geometry of a solve at box: KB window, colouring stride, padded box, limits, rim radius,
    !! wrap table and deposition envelope. fft_nthreads only sizes the plans of the probe image.
    subroutine new_lattice( self, box, smpd, fft_nthreads )
        class(pcg_lattice), intent(inout) :: self
        integer,            intent(in)    :: box
        real,               intent(in)    :: smpd
        integer, optional,  intent(in)    :: fft_nthreads
        type(image) :: tmp
        integer     :: lo, hi, i, nthr
        real        :: rlim
        call self%kill_lattice
        self%box    = box
        self%smpd   = smpd
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        ! odd width centred on nint(loc): negating loc negates the window as a set,
        ! which the h>=0 folds rely on
        self%wdim   = 2*self%iwinsz + 1
        ! same-colour h-lines stay a full window apart along one axis after rotation
        ! (max-norm >= Euclidean norm / sqrt(3)), as in reconstructor%insert_plane_oversamp
        self%stride = ceiling(sqrt(3.0) * real(self%wdim))
        ! oversampling: the unknown lives on the native box, every Fourier operation
        ! on the padf-times padded lattice, as in Fourier gridding
        self%padf    = OSMPL_PAD_FAC
        self%boxpd   = self%padf * box
        self%pad_off = (self%boxpd - box) / 2
        self%Rnat    = box / 2
        nthr = 1
        if( present(fft_nthreads) ) nthr = max(1, fft_nthreads)
        call tmp%new([self%boxpd,self%boxpd,self%boxpd], smpd, wthreads=nthr > 1, fft_nthreads=nthr)
        call tmp%fft()
        self%lims3 = tmp%loop_lims(3)
        self%cdim  = tmp%get_array_shape()
        call tmp%kill
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

    subroutine kill_lattice( self )
        class(pcg_lattice), intent(inout) :: self
        if( allocated(self%wrap)   ) deallocate(self%wrap)
        if( allocated(self%dep1d)  ) deallocate(self%dep1d)
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        self%l_mask = .false.
        self%box = 0; self%boxpd = 0; self%padf = 1; self%pad_off = 0; self%Rnat = 0
        self%iwinsz = 0; self%wdim = 0; self%stride = 0; self%lims3 = 0; self%cdim = 0
        self%wlims = 0; self%sq_rim = 0
    end subroutine kill_lattice

end module simple_pcg_lattice
