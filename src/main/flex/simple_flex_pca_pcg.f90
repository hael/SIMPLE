!@descr: flex_pca coupled M-step on the PCG operator (rec_backend=pcg): pair-weighted Gram kernels and the
!  right-hand sides on the 2x lattice, the per-voxel coupled solve kept as the (floored) block-Jacobi
!  preconditioner, the support inside the solve, CG to a relative-residual tolerance
!  (doc/implementation_notes/flex_pca_envelope_support.md, 3.4), mirroring reconstructor_pcg.
!
!  Conventions. The CG variable is the physical basis volume u on the native lattice. The right-hand side
!  b = S^H y is the exact adjoint of nonuniform sampling: the pair-weighted plane samples deposited through
!  the KB window at DOUBLED coordinates on the 2x lattice, inverse-transformed, cropped and divided by the
!  2x deposition envelope (the 1x KB deposit of the gridding path carries the native gather envelope E and a
!  10-20% interpolation error; the 2x deposit is exact to ~1%). The operator is the bare Toeplitz Gram
!  T = S^H S of the same samples through the doubled-coordinate kernel (scale padf**3 relative to the
!  gridding density under SIMPLE's forward-normalised FFT), band-limited to the Nyquist ball and bracketed by
!  the hard support P. The preconditioner is the per-voxel coupled divide on the gridding density rho
!  (the Fourier-diagonal Jacobi inverse of T), with a shell-relative floor so that it stays bounded on the
!  unsampled voxels every real-space operation (the support) leaks into. Every solve with a positive budget
!  starts from ZERO (the production reconstruct3D_pcg convention), so maxits_pcg buys the same number of
!  corrections here as it does there; budget 0 is the explicit gridding mode and ships
!  solve_coupled_basis_exp untouched. put_back writes E*u so the unchanged tail (inverse KB envelope,
!  FSC-Wiener, band limit, soft mask) ships u.
module simple_flex_pca_pcg
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$ use omp_lib, only: omp_get_thread_num, omp_get_max_threads, omp_in_parallel
use simple_core_module_api
use simple_image,             only: image
use simple_reconstructor,     only: reconstructor
use simple_kbinterpol,        only: kbinterpol
use simple_gridding,          only: kb_stencil_envelope_1d
use simple_math_ft,           only: cyci_1d
use simple_math,              only: ceil_div, floor_div
use simple_cartesian_fourier, only: center_embed_real3d, center_crop_real3d
use simple_flex_reconstructor_latent_ops, only: solve_coupled_basis_exp, pair_index, projected_model_kfromto
use simple_flex_pca_util,     only: cov_env_int
use simple_rnd,               only: gasdev
use simple_parameters,        only: parameters
use simple_imghead,           only: find_ldim_nptcls
implicit none

public :: flex_pcg_t, flex_pcg_outcome_t, test_flex_pcg_operator
public :: flex_pcg_support_volume, flex_env_init, flex_env_active, flex_window_apply, flex_window_apply_rec
public :: flex_pcg_install_window
public :: flex_complement_ncomp, flex_complement_active, flex_window_apply_blk, flex_window_apply_rec_blk
public :: flex_mskfile_set
private
#include "simple_local_flags.inc"

!> kernel scale relative to the gridding density: the native-lattice normal operator is (1/N^3) T under
!! the forward-normalised FFT, the 2x-lattice kernel path returns (1/(2N)^3) T, so the ratio is padf**3
real,    parameter :: FLEX_PCG_KSCALE = real(OSMPL_PAD_FAC)**3
!> diminishing-returns stop on the relative update dx/x (production PCG_XTOL)
real,    parameter :: FLEX_PCG_XTOL   = 1.5e-2
!> recompute b - Hx from scratch every so many iterations (production RESID_REPLACE)
integer, parameter :: FLEX_PCG_RESID_REPLACE = 25
!> shell-relative floor of the preconditioner density (production RHO_FLOOR_FRAC)
real,    parameter :: FLEX_PCG_RHO_FLOOR_FRAC = 1.0e-2
!> relative ridge of the per-voxel coupled solve (as COUPLED_MSTEP_RIDGE_REL in the gridding solve)
real(dp), parameter :: FLEX_PCG_RIDGE_REL = 1.0d-8
!> default Tikhonov term relative to the low-band mean density (production PCG_LAMBDA); env SIMPLE_COV_PCG_LAMBDA
real,     parameter :: FLEX_PCG_LAMBDA_REL_DEFAULT = 1.0e-3
!> resampled support values below this are outside the hard domain (the note's 4.2 item 4; no larger
!! than the solver's PCG_SUPPORT_DIV_MIN = 0.1)
real,     parameter :: FLEX_PCG_SUPPORT_FLOOR = 0.05

!> the envelope support at the covariance box (pcg_mskfile, step 2 of the note), loaded once per process
type(image), save :: flex_env_img
logical,     save :: l_flex_env = .false., l_flex_env_checked = .false.
!> complement block (SIMPLE_COV_COMPLEMENT=k_out, default 0): k_out extra basis components supported on the
!! sphere MINUS the envelope, fitted jointly with the envelope block so the heterogeneity outside the
!! region (micelle, belt) is explained by its own components instead of leaking into the region's
!! latents through the projection overlap (RECOVAR's --use-complement-mask, in the EM). Block 1 = the
!! envelope window, block 2 = the complement window; the two are disjoint up to the soft edges.
type(image), save :: flex_cenv_img
integer,     save :: flex_cenv_ncomp = -1     !< -1: env not read yet
character(len=*), parameter, public :: FLEX_PCG_STOP_INDEFINITE = 'indefinite'

type :: flex_pcg_outcome_t
    character(len=24) :: stop_reason = 'not_started'
    integer :: iteration_count = 0, requested_maxits = 0
    real    :: initial_rel_residual = 0., final_rel_residual = 0., final_rel_update = 0.
    real    :: rhs_norm = 0., start_norm = 0.
    real    :: start_corr = 0., start_scale = 0.   !< corr(b, B x0) and <b,Bx0>/<Bx0,Bx0> of the warm start
    integer :: iters_to_1e2 = 0                    !< first iteration with rel resid <= 1e-2 (0: never)
    logical :: converged = .false., cold_restart_used = .false.
    real    :: seconds = 0.
end type flex_pcg_outcome_t

type :: flex_pcg_t
    private
    integer :: box = 0, boxpd = 0, padf = 1, Rnat = 0, ncomp = 0, npairs = 0
    ! persistent apply_operator scratch: this runs once per CG iteration, so allocating it per
    ! call made the solve mmap- and page-fault-bound rather than arithmetic-bound
    complex, allocatable :: op_cq(:,:,:,:)  !< per-component transforms on the padded lattice
    real    :: smpd = 1.0
    type(kbinterpol) :: kbwin
    integer :: iwinsz = 0, wdim = 0, stride = 0
    integer :: lims3(3,2) = 0, wlims(2) = 0, sq_rim = 0, cdim(3) = 0
    integer, allocatable :: wrap(:)
    ! ---- band-ball index lists (memory): the accumulators and the packed kernels live only on the
    ! lattice points the band can reach. nband is the plane band in native units (set_band); the
    ! expanded list holds every expanded-lattice point within R = pf*(nband+1/2) + sqrt(3)*iwinsz + 1
    ! of a wrap image of the origin, the physical list every physical point the fold reads into.
    ! Slot maps are 0 outside the lists. eorder/efold replay fold_accum's (m,k,h>=0) loop order, so
    ! the fold is bitwise the fold of the dense arrays.
    integer :: nband = 0, nexp = 0, npk = 0, nfold = 0
    integer, allocatable :: emap(:,:,:)    !< (expanded index triple) -> expanded slot
    integer, allocatable :: pmap(:,:,:)    !< (physical i,j,k) -> physical slot
    integer, allocatable :: eijk(:,:)      !< (3,nexp) expanded slot -> index triple
    integer, allocatable :: pijk(:,:)      !< (3,npk)  physical slot -> (i,j,k)
    integer, allocatable :: eorder(:), efold(:) !< (nfold) fold pairs: kpk(:,efold(t)) += kacc(:,eorder(t))
    logical :: l_lists = .false.
    real,    allocatable :: env(:,:,:)      !< native-lattice KB gather envelope E, unity at the centre
    real,    allocatable :: dep(:,:,:)      !< deposition envelope of the 2x lattice (kernel build)
    real,    allocatable :: dep_c(:,:,:)    !< the same cropped to the native box (RHS deapodization)
    real,    allocatable :: mask(:,:,:)     !< hard solve support (window > 0)
    real,    allocatable :: window(:,:,:)   !< soft window (the tail's mask3D_soft multiplies by it)
    real,    allocatable :: khat(:,:,:,:)   !< (npairs, cdim) finalized pair kernels
    real,    allocatable :: ridge(:,:)      !< (ncomp, nshell) optional Fourier shell ridge in the operator
    real,    allocatable :: rhofl(:,:)      !< (ncomp, 0:Rnat) shell-relative preconditioner floors
    real,    allocatable :: lam(:)          !< (ncomp) Tikhonov term of the operator and preconditioner
    real    :: lam_rel = FLEX_PCG_LAMBDA_REL_DEFAULT !< lam = lam_rel x low-band mean of the diagonal density
    logical :: l_mask = .false., l_kernel = .false., l_ridge = .false., l_floor = .false.
    type(image) :: wimg                     !< persistent boxpd^3 work image (keeps its plans)
    type(image) :: nimg                     !< persistent box^3 work image (band limit, ridge, precond)
    ! Per-thread transform pool. apply_operator's component transforms are INDEPENDENT, so they run
    ! one-per-thread with SINGLE-THREADED FFTW plans -- the idiom image%construct_thread_safe_tmp_imgs
    ! already uses (simple_image_core.f90:126, each pooled image built with wthreads=.false.).
    ! Threading inside a transform instead was measured to be a net harm: every stack sample landed
    ! in os_sem_down inside libfftw3f_threads. Width is min(nthr, ncomp) -- more threads than
    ! components would idle.
    type(image), allocatable :: wpool(:)    !< per-thread boxpd^3 images
    type(image), allocatable :: npool(:)    !< per-thread box^3 images (band limit)
    real,    allocatable :: tp_work(:,:,:,:)!< (box,box,box,nthr_pool)
    real,    allocatable :: tp_emb(:,:,:,:) !< (boxpd,boxpd,boxpd,nthr_pool)
    complex, allocatable :: tp_acc(:,:,:,:) !< (cdim,nthr_pool)
    integer :: nthr_pool = 0
    complex, allocatable :: pc_cq(:,:,:,:) !< persistent apply_precond accumulator
    complex, allocatable :: op_hq(:,:,:,:) !< coupled-multiply output, one plane set per component
    integer, allocatable :: pidx(:,:)      !< pair_index(min(q,r),max(q,r)) lookup, built once
    logical :: exists = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: is_ready
    procedure :: set_window_sphere
    procedure :: set_window_volume
    procedure :: set_band
    procedure :: get_npk
    procedure :: get_nexp
    procedure :: get_pijk
    procedure :: list_signature
    procedure :: alloc_accum
    procedure :: alloc_packed
    procedure :: alloc_rhs_accum
    procedure :: alloc_rhs_packed
    procedure :: accumulate
    procedure :: accumulate_rhs
    procedure :: fold_accum
    procedure :: fold_rhs
    procedure :: finalize
    procedure :: finalize_rhs
    procedure :: set_ridge
    procedure :: clear_ridge
    procedure :: set_lambda_relative
    procedure :: solve
    procedure :: cg_core
    procedure :: apply_operator
    procedure :: apply_precond
    procedure :: get_cdim
    procedure :: get_npairs
    procedure :: bytes_accum
    procedure :: bytes_packed
    procedure :: bytes_rhs_accum
    procedure :: bytes_rhs_packed
    procedure, private :: win_wraps
    procedure, private :: build_lists
    procedure, private :: require_lists
    procedure, private :: dot_all
    procedure, private :: bandlimit
    procedure, private :: bandlimit_img
    procedure, private :: ensure_pool
    procedure, private :: prep_floor
end type flex_pcg_t

contains

    ! ---------------- lifecycle ----------------

    !> geometry of the solve at the covariance box: the native lattice holds the unknown, the padf-times
    !! padded lattice the kernels and the right-hand sides; tables as in reconstructor_pcg%new
    subroutine new( self, box, smpd, ncomp )
        class(flex_pcg_t), intent(inout) :: self
        integer,           intent(in)    :: box, ncomp
        real,              intent(in)    :: smpd
        type(image)       :: tmp
        real, allocatable :: env1d(:), dep1d(:)
        real    :: rlim
        integer :: lo, hi, i, j, k, off
        call self%kill
        self%box    = box
        self%smpd   = smpd
        self%ncomp  = ncomp
        self%npairs = (ncomp*(ncomp+1))/2
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        self%wdim   = 2*self%iwinsz + 1
        self%stride = ceiling(sqrt(3.0) * real(self%wdim))
        self%padf   = OSMPL_PAD_FAC
        self%boxpd  = self%padf * box
        self%Rnat   = box / 2
        call tmp%new([self%boxpd,self%boxpd,self%boxpd], smpd)
        call tmp%fft()
        self%lims3 = tmp%loop_lims(3)
        self%cdim  = tmp%get_array_shape()
        call tmp%kill
        self%wlims = self%lims3(2,:)
        rlim = real(min(self%wlims(2) - self%iwinsz, -self%wlims(1) - self%iwinsz)) - 0.5
        self%sq_rim = max(0, int((rlim / real(self%padf))**2) - 1)
        lo = self%wlims(1) - self%iwinsz - 1
        hi = self%wlims(2) + self%iwinsz + 1
        allocate(self%wrap(lo:hi))
        do i = lo, hi
            self%wrap(i) = cyci_1d(self%wlims, i)
        end do
        ! the gather envelope of the native lattice (its reciprocal is the tail's gridding correction)
        call kb_stencil_envelope_1d(self%kbwin, box, env1d)
        allocate(self%env(box,box,box))
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    self%env(i,j,k) = env1d(i)*env1d(j)*env1d(k)
                end do
            end do
        end do
        ! the deposition envelope of the doubled-coordinate scatter on the 2x lattice, and its native crop
        call kb_stencil_envelope_1d(self%kbwin, self%boxpd, dep1d)
        allocate(self%dep(self%boxpd,self%boxpd,self%boxpd))
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, self%boxpd
            do j = 1, self%boxpd
                do i = 1, self%boxpd
                    self%dep(i,j,k) = dep1d(i)*dep1d(j)*dep1d(k)
                end do
            end do
        end do
        !$omp end parallel do
        off = (self%boxpd - box)/2
        allocate(self%dep_c(box,box,box), source=self%dep(off+1:off+box, off+1:off+box, off+1:off+box))
        deallocate(env1d, dep1d)
        call self%wimg%new([self%boxpd,self%boxpd,self%boxpd], smpd)
        call self%nimg%new([box,box,box], smpd)
        ! pair_index LUT: the coupled multiply evaluated pair_index(min(q,r),max(q,r)) once per
        ! (component, component, voxel), which at rank 10 and cdim 129x256x256 is 845 million
        ! out-of-line calls per operator application
        allocate(self%pidx(ncomp,ncomp))
        do j = 1, ncomp
            do i = 1, ncomp
                self%pidx(i,j) = pair_index(min(i,j), max(i,j))
            end do
        end do
        self%lam_rel = FLEX_PCG_LAMBDA_REL_DEFAULT
        call env_real('SIMPLE_COV_PCG_LAMBDA', self%lam_rel)
        self%exists = .true.
    end subroutine new

    !> real-valued environment override (unset or unreadable leaves the default)
    subroutine env_real( name, val )
        character(len=*), intent(in)    :: name
        real,             intent(inout) :: val
        character(len=64) :: buf
        integer :: lenv, stat, ios
        real    :: v
        call get_environment_variable(name, buf, lenv, stat)
        if( stat /= 0 .or. lenv < 1 ) return
        read(buf(1:lenv), *, iostat=ios) v
        if( ios /= 0 ) return
        if( ieee_is_finite(v) .and. v >= 0.0 ) val = v
    end subroutine env_real

    !> Tikhonov term relative to the low-band mean of each component's diagonal density (set before
    !! the solve; 0 disables it)
    subroutine set_lambda_relative( self, lam_rel )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: lam_rel
        self%lam_rel = max(0.0, lam_rel)
        self%l_floor = .false.   ! lam is derived with the floors
    end subroutine set_lambda_relative

    subroutine kill( self )
        class(flex_pcg_t), intent(inout) :: self
        integer :: i_pool
        if( allocated(self%wrap)   ) deallocate(self%wrap)
        if( allocated(self%env)    ) deallocate(self%env)
        if( allocated(self%dep)    ) deallocate(self%dep)
        if( allocated(self%dep_c)  ) deallocate(self%dep_c)
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        if( allocated(self%khat)   ) deallocate(self%khat)
        if( allocated(self%ridge)  ) deallocate(self%ridge)
        if( allocated(self%rhofl)  ) deallocate(self%rhofl)
        if( allocated(self%lam)    ) deallocate(self%lam)
        if( allocated(self%op_cq)  ) deallocate(self%op_cq)
        if( allocated(self%tp_work) ) deallocate(self%tp_work)
        if( allocated(self%tp_emb)  ) deallocate(self%tp_emb)
        if( allocated(self%tp_acc)  ) deallocate(self%tp_acc)
        if( allocated(self%pc_cq)   ) deallocate(self%pc_cq)
        if( allocated(self%op_hq)   ) deallocate(self%op_hq)
        if( allocated(self%pidx)    ) deallocate(self%pidx)
        if( allocated(self%emap)    ) deallocate(self%emap)
        if( allocated(self%pmap)    ) deallocate(self%pmap)
        if( allocated(self%eijk)    ) deallocate(self%eijk)
        if( allocated(self%pijk)    ) deallocate(self%pijk)
        if( allocated(self%eorder)  ) deallocate(self%eorder)
        if( allocated(self%efold)   ) deallocate(self%efold)
        self%nband = 0; self%nexp = 0; self%npk = 0; self%nfold = 0; self%l_lists = .false.
        if( allocated(self%wpool) )then
            do i_pool = 1, size(self%wpool)
                call self%wpool(i_pool)%kill
            end do
            deallocate(self%wpool)
        endif
        if( allocated(self%npool) )then
            do i_pool = 1, size(self%npool)
                call self%npool(i_pool)%kill
            end do
            deallocate(self%npool)
        endif
        self%nthr_pool = 0
        if( self%exists )then
            call self%wimg%kill
            call self%nimg%kill
        endif
        self%l_mask = .false.; self%l_kernel = .false.; self%l_ridge = .false.; self%l_floor = .false.
        self%exists = .false.
    end subroutine kill

    pure logical function is_ready( self )
        class(flex_pcg_t), intent(in) :: self
        is_ready = self%exists .and. self%l_kernel
    end function is_ready

    pure function get_cdim( self ) result( cdim )
        class(flex_pcg_t), intent(in) :: self
        integer :: cdim(3)
        cdim = self%cdim
    end function get_cdim

    pure integer function get_npairs( self )
        class(flex_pcg_t), intent(in) :: self
        get_npairs = self%npairs
    end function get_npairs

    !> bytes of one full-range kernel accumulator (all pairs)
    pure real(dp) function bytes_accum( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_accum = 4.d0 * real(self%npairs,dp) * real(self%nexp,dp)
    end function bytes_accum

    !> bytes of one packed (transport) kernel set
    pure real(dp) function bytes_packed( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_packed = 4.d0 * real(self%npairs,dp) * real(self%npk,dp)
    end function bytes_packed

    !> bytes of one full-range right-hand-side accumulator (all components, complex)
    pure real(dp) function bytes_rhs_accum( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_rhs_accum = 8.d0 * real(self%ncomp,dp) * real(self%nexp,dp)
    end function bytes_rhs_accum

    !> bytes of one packed (transport) right-hand-side set
    pure real(dp) function bytes_rhs_packed( self )
        class(flex_pcg_t), intent(in) :: self
        bytes_rhs_packed = 8.d0 * real(self%ncomp,dp) * real(self%npk,dp)
    end function bytes_rhs_packed

    ! ---------------- support ----------------

    !> spherical window at the solve box: the soft sphere the gridding tail applies, its hard domain
    !! window > 0 as the solve support (reconstructor_pcg%set_mask)
    subroutine set_window_sphere( self, mskrad )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: mskrad
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
        write(logfhandle,'(A,F7.2,A,F6.3,A,F6.3,A,F6.3)') '>>> FLEX_PCA PCG support: sphere radius ', mskrad, &
            &' px, window mean ', sum(self%window)/real(size(self%window)), ', hard support fraction ', &
            &sum(self%mask)/real(size(self%mask)), ', nominal sphere fraction ', &
            &(4.0/3.0)*PI*mskrad**3/real(self%box)**3
    end subroutine set_window_sphere

    !> caller-supplied [0,1] window at the solve box (clipped), the same contract as set_window_sphere
    subroutine set_window_volume( self, mskvol )
        class(flex_pcg_t), intent(inout) :: self
        class(image),      intent(in)    :: mskvol
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

    subroutine install_support( self )
        class(flex_pcg_t), intent(inout) :: self
        if( allocated(self%mask) ) deallocate(self%mask)
        allocate(self%mask(self%box,self%box,self%box), source=merge(1.0, 0.0, self%window > 0.0))
        self%l_mask = .true.
    end subroutine install_support

    ! ---------------- accumulation ----------------

    !> full-range accumulator of the doubled-coordinate kernel scatter, one row per pair;
    !! 1-based on every axis (index = logical - lims3(:,1) + 1): the scatter sees it through an
    !! assumed-shape dummy and fold_accum through an allocatable one, and the two must agree
    subroutine alloc_accum( self, kacc )
        class(flex_pcg_t),  intent(in)  :: self
        real, allocatable,  intent(out) :: kacc(:,:)
        call self%require_lists('alloc_accum')
        allocate(kacc(self%npairs, self%nexp), source=0.0)
    end subroutine alloc_accum

    !> packed (h >= 0, real) kernel set: the transport and reduction form
    subroutine alloc_packed( self, kpk )
        class(flex_pcg_t),  intent(in)  :: self
        real, allocatable,  intent(out) :: kpk(:,:)
        call self%require_lists('alloc_packed')
        allocate(kpk(self%npairs, self%npk), source=0.0)
    end subroutine alloc_packed

    !> full-range accumulator of the doubled-coordinate right-hand-side scatter, one row per component
    subroutine alloc_rhs_accum( self, racc )
        class(flex_pcg_t),     intent(in)  :: self
        complex, allocatable,  intent(out) :: racc(:,:)
        call self%require_lists('alloc_rhs_accum')
        allocate(racc(self%ncomp, self%nexp), source=cmplx(0.,0.))
    end subroutine alloc_rhs_accum

    !> packed (h >= 0, complex) right-hand-side set: the transport and reduction form
    subroutine alloc_rhs_packed( self, rpk )
        class(flex_pcg_t),     intent(in)  :: self
        complex, allocatable,  intent(out) :: rpk(:,:)
        call self%require_lists('alloc_rhs_packed')
        allocate(rpk(self%ncomp, self%npk), source=cmplx(0.,0.))
    end subroutine alloc_rhs_packed

    !> The plane band in native units. Builds the band-ball index lists; every accumulate / packed
    !! allocation needs them. A band of Rnat (the default when nothing is known) covers the whole
    !! lattice reachable by the stencil.
    subroutine set_band( self, nband )
        class(flex_pcg_t), intent(inout) :: self
        integer,           intent(in)    :: nband
        integer :: nb
        nb = max(1, min(self%Rnat, nband))
        if( self%l_lists .and. nb == self%nband ) return
        self%nband = nb
        call self%build_lists
    end subroutine set_band

    subroutine require_lists( self, who )
        class(flex_pcg_t), intent(in) :: self
        character(len=*),  intent(in) :: who
        if( .not. self%l_lists ) THROW_HARD(who//': band lists not built; call set_band first')
    end subroutine require_lists

    pure integer function get_npk( self )
        class(flex_pcg_t), intent(in) :: self
        get_npk = self%npk
    end function get_npk

    pure integer function get_nexp( self )
        class(flex_pcg_t), intent(in) :: self
        get_nexp = self%nexp
    end function get_nexp

    !> physical (i,j,k) of packed slot t
    pure function get_pijk( self, t ) result( ijk )
        class(flex_pcg_t), intent(in) :: self
        integer,           intent(in) :: t
        integer :: ijk(3)
        ijk = self%pijk(:,t)
    end function get_pijk

    !> what two processes must agree on for their packed arrays to be slot-compatible
    pure function list_signature( self ) result( sig )
        class(flex_pcg_t), intent(in) :: self
        integer :: sig(6)
        sig = [self%box, self%boxpd, self%ncomp, self%nband, self%nexp, self%npk]
    end function list_signature

    !> Expanded list: every expanded-lattice index within R of a wrap image of the origin, R the
    !! farthest a stencil point can land for a plane point inside the band (|loc| < pf*(nband+1/2),
    !! nint rounding 1/2, stencil half-diagonal sqrt(3)*iwinsz, +1 margin). Physical list and fold
    !! pairs: fold_accum's (m, k, h >= 0) loop, reading kacc at (wrap(h), k, m) into the physical
    !! address of (h,k,m), replayed here once and stored in the same order.
    subroutine build_lists( self )
        class(flex_pcg_t), intent(inout) :: self
        real    :: R, R2
        integer :: n1, n2, n3, ih, ik, im, h, k, m, hh, hm, km, mm, ne, np, nf, phys(3), es, ps, pass
        if( allocated(self%emap)   ) deallocate(self%emap)
        if( allocated(self%pmap)   ) deallocate(self%pmap)
        if( allocated(self%eijk)   ) deallocate(self%eijk)
        if( allocated(self%pijk)   ) deallocate(self%pijk)
        if( allocated(self%eorder) ) deallocate(self%eorder)
        if( allocated(self%efold)  ) deallocate(self%efold)
        R  = real(self%padf) * (real(self%nband) + 0.5) + sqrt(3.0) * real(self%iwinsz) + 1.0
        R2 = R * R
        n1 = self%lims3(1,2) - self%lims3(1,1) + 1
        n2 = self%lims3(2,2) - self%lims3(2,1) + 1
        n3 = self%lims3(3,2) - self%lims3(3,1) + 1
        allocate(self%emap(n1,n2,n3), source=0)
        do pass = 1, 2
            ne = 0
            do im = 1, n3
                m  = im - 1 + self%lims3(3,1)
                mm = min(abs(m), abs(m - self%boxpd), abs(m + self%boxpd))
                do ik = 1, n2
                    k  = ik - 1 + self%lims3(2,1)
                    km = min(abs(k), abs(k - self%boxpd), abs(k + self%boxpd))
                    do ih = 1, n1
                        h  = ih - 1 + self%lims3(1,1)
                        hm = min(abs(h), abs(h - self%boxpd), abs(h + self%boxpd))
                        if( real(hm*hm + km*km + mm*mm) > R2 ) cycle
                        ne = ne + 1
                        if( pass == 2 )then
                            self%emap(ih,ik,im) = ne
                            self%eijk(:,ne)     = [ih,ik,im]
                        endif
                    end do
                end do
            end do
            if( pass == 1 ) allocate(self%eijk(3,max(1,ne)), source=0)
        end do
        self%nexp = ne
        allocate(self%pmap(self%cdim(1),self%cdim(2),self%cdim(3)), source=0)
        do pass = 1, 2
            np = 0; nf = 0
            do m = self%lims3(3,1), self%lims3(3,2)
                do k = self%lims3(2,1), self%lims3(2,2)
                    do h = 0, -self%lims3(1,1)
                        hh = self%wrap(h)
                        ih = hh - self%lims3(1,1) + 1
                        ik = k  - self%lims3(2,1) + 1
                        im = m  - self%lims3(3,1) + 1
                        es = self%emap(ih,ik,im)
                        if( es == 0 ) cycle
                        phys = self%wimg%comp_addr_phys(h,k,m)
                        ps   = self%pmap(phys(1),phys(2),phys(3))
                        if( ps == 0 )then
                            np = np + 1
                            if( pass == 2 )then
                                self%pmap(phys(1),phys(2),phys(3)) = np
                                self%pijk(:,np) = phys
                            endif
                            ps = np
                        endif
                        nf = nf + 1
                        if( pass == 2 )then
                            self%eorder(nf) = es
                            self%efold(nf)  = ps
                        endif
                    end do
                end do
            end do
            if( pass == 1 )then
                allocate(self%pijk(3,max(1,np)), self%eorder(max(1,nf)), self%efold(max(1,nf)), source=0)
            else
                exit
            endif
            ! pass 1 marked nothing in pmap (ps stays 0 for every new point), so the count is exact
            ! only if pmap is reset between passes
            self%pmap = 0
        end do
        self%npk = np; self%nfold = nf
        self%l_lists = .true.
        write(logfhandle,'(A,I0,A,F6.2,A,I0,A,F6.2,A,I0,A)') '>>> FLEX_PCA PCG band lists: nband=', self%nband, &
            &'  expanded ', 100.0*real(self%nexp)/real(n1)/real(n2)/real(n3), ' % (', self%nexp, &
            &' points), physical ', 100.0*real(self%npk)/real(product(self%cdim)), ' % (', self%npk, ' points)'
        call flush(logfhandle)
    end subroutine build_lists

    !> pair-weighted Gram kernels of one batch: |T_i|^2 M_i(q,r) scattered at padf*R_i*[h,k,0] through the
    !! KB window, the same native-subset plane samples and Friedel storage as the coupled insert
    !! (insert_planes_oversamp_coupled_batch_scaled); the flex plane lattice is already padf-oversampled
    subroutine accumulate( self, kacc, se, orientations, fpls, density_scales, valid, nrecords )
        class(flex_pcg_t),   intent(in)    :: self
        real,                intent(inout) :: kacc(:,:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: density_scales(:,:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(ori) :: o_sym
        real,    allocatable :: rotmats(:,:,:,:), dpack(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3), ctfsq_raw
        integer :: i, isym, nsym, l, h, k, hp, kp, q, r, pf, i0(3), nyq_eff, fpllims_pd(3,2)
        integer :: h_sq, k_max_h, k_lo, k_hi
        if( nrecords < 1 ) return
        if( size(kacc,1) /= self%npairs ) THROW_HARD('pair kernel accumulator has the wrong leading extent; accumulate')
        call self%require_lists('accumulate')
        if( size(kacc,2) /= self%nexp ) THROW_HARD('pair kernel accumulator is not on the band list; accumulate')
        nsym = se%get_nsym()
        pf   = self%padf
        allocate(rotmats(3,3,nsym,nrecords), source=0.)
        allocate(dpack(self%npairs,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            if( .not. allocated(fpls(i)%ctfsq_plane) ) THROW_HARD('ctfsq plane does not exist; accumulate')
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd     = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = self%Rnat
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq/pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            do r = 1, self%ncomp
                do q = 1, r
                    dpack(pair_index(q,r),i) = real(density_scales(q,r,i))
                end do
            end do
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,isym,l,h,k,hp,kp,h_sq,k_max_h,k_lo,k_hi,ctfsq_raw,loc,i0,w,rot) &
        !$omp proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            do isym = 1, nsym
                rot = rotmats(:,:,isym,i)
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), self%stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo = max(fpllims(2,1,i),-k_max_h)
                        k_hi = min(fpllims(2,2,i), k_max_h)
                        hp   = h*pf
                        do k = k_lo, k_hi
                            kp = k*pf
                            if( kp <= 0 )then
                                ctfsq_raw = fpls(i)%ctfsq_plane(hp,kp)
                            else
                                ctfsq_raw = fpls(i)%ctfsq_plane(-hp,-kp)
                            endif
                            if( ctfsq_raw <= TINY ) cycle
                            loc = real(pf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( self%win_wraps(i0) ) cycle   ! rim: serial pass below
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_pairs_nowrap(self, i0, w, dpack(:,i), ctfsq_raw, kacc)
                        end do
                    end do
                    !$omp end do
                end do
                ! wrapping rim, serialized: the colouring does not survive folding
                !$omp single
                do h = fpllims(1,1,i), fpllims(1,2,i)
                    h_sq = h*h
                    if( h_sq > nyq_disks(i) ) cycle
                    k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                    k_lo = max(fpllims(2,1,i),-k_max_h)
                    k_hi = min(fpllims(2,2,i), k_max_h)
                    hp   = h*pf
                    do k = k_lo, k_hi
                        if( h_sq + k*k <= self%sq_rim ) cycle   ! provably cannot wrap
                        kp = k*pf
                        if( kp <= 0 )then
                            ctfsq_raw = fpls(i)%ctfsq_plane(hp,kp)
                        else
                            ctfsq_raw = fpls(i)%ctfsq_plane(-hp,-kp)
                        endif
                        if( ctfsq_raw <= TINY ) cycle
                        loc = real(pf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. self%win_wraps(i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_pairs_wrap(self, i0, w, dpack(:,i), ctfsq_raw, kacc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, dpack, fpllims, nyq_disks)
    end subroutine accumulate

    !> right-hand sides of one batch on the 2x lattice: the data scales E[z_iq] times the transfer-
    !! conjugated plane samples (with the plane padding factor, as the coupled insert deposits them),
    !! scattered at padf*R_i*[h,k,0] through the KB window; the same samples and Friedel storage as
    !! accumulate, so the kernels and the right-hand sides are one sample set
    subroutine accumulate_rhs( self, racc, se, orientations, fpls, data_scales, valid, nrecords )
        class(flex_pcg_t),   intent(in)    :: self
        complex,             intent(inout) :: racc(:,:)
        class(sym),          intent(inout) :: se
        type(ori),           intent(inout) :: orientations(:)
        type(fplane_type),   intent(in)    :: fpls(:)
        real(dp),            intent(in)    :: data_scales(:,:)
        logical,             intent(in)    :: valid(:)
        integer,             intent(in)    :: nrecords
        type(ori) :: o_sym
        real,    allocatable :: rotmats(:,:,:,:), dsc(:,:)
        integer, allocatable :: fpllims(:,:,:), nyq_disks(:)
        complex :: cmplx_raw, vals(self%ncomp)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3), pf2
        integer :: i, isym, nsym, l, h, k, hp, kp, pf, i0(3), nyq_eff, fpllims_pd(3,2)
        integer :: h_sq, k_max_h, k_lo, k_hi
        if( nrecords < 1 ) return
        if( size(racc,1) /= self%ncomp ) THROW_HARD('rhs accumulator has the wrong leading extent; accumulate_rhs')
        call self%require_lists('accumulate_rhs')
        if( size(racc,2) /= self%nexp ) THROW_HARD('rhs accumulator is not on the band list; accumulate_rhs')
        nsym = se%get_nsym()
        pf   = self%padf
        pf2  = real(pf*pf)
        allocate(rotmats(3,3,nsym,nrecords), source=0.)
        allocate(dsc(self%ncomp,nrecords), source=0.)
        allocate(fpllims(3,2,nrecords), nyq_disks(nrecords), source=0)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            if( .not. allocated(fpls(i)%transfer_plane) ) THROW_HARD('transfer plane does not exist; accumulate_rhs')
            rotmats(:,:,1,i) = orientations(i)%get_mat()
            do isym = 2, nsym
                call se%apply(orientations(i), isym, o_sym)
                rotmats(:,:,isym,i) = o_sym%get_mat()
            end do
            fpllims_pd     = fpls(i)%frlims
            fpllims(:,:,i) = fpllims_pd
            fpllims(1,1,i) = ceil_div (fpllims_pd(1,1), pf)
            fpllims(1,2,i) = floor_div(fpllims_pd(1,2), pf)
            fpllims(2,1,i) = ceil_div (fpllims_pd(2,1), pf)
            fpllims(2,2,i) = floor_div(fpllims_pd(2,2), pf)
            nyq_eff = self%Rnat
            if( fpls(i)%nyq > 0 ) nyq_eff = min(nyq_eff, max(1, fpls(i)%nyq/pf))
            nyq_disks(i) = nyq_eff * (nyq_eff + 1)
            dsc(:,i) = real(data_scales(1:self%ncomp,i))
        end do
        call o_sym%kill
        !$omp parallel default(shared) private(i,isym,l,h,k,hp,kp,h_sq,k_max_h,k_lo,k_hi,cmplx_raw,vals,loc,i0,w,rot) &
        !$omp proc_bind(close)
        do i = 1, nrecords
            if( .not. valid(i) ) cycle
            do isym = 1, nsym
                rot = rotmats(:,:,isym,i)
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = fpllims(1,1,i)+l, fpllims(1,2,i), self%stride
                        h_sq = h*h
                        if( h_sq > nyq_disks(i) ) cycle
                        k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                        k_lo = max(fpllims(2,1,i),-k_max_h)
                        k_hi = min(fpllims(2,2,i), k_max_h)
                        hp   = h*pf
                        do k = k_lo, k_hi
                            kp = k*pf
                            if( kp <= 0 )then
                                cmplx_raw = conjg(fpls(i)%transfer_plane(hp,kp))*fpls(i)%cmplx_plane(hp,kp)
                            else
                                cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-hp,-kp))*fpls(i)%cmplx_plane(-hp,-kp))
                            endif
                            if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw)) <= TINY ) cycle
                            loc = real(pf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( self%win_wraps(i0) ) cycle   ! rim: serial pass below
                            vals = dsc(:,i) * (pf2 * cmplx_raw)
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_rhs_nowrap(self, i0, w, vals, racc)
                        end do
                    end do
                    !$omp end do
                end do
                !$omp single
                do h = fpllims(1,1,i), fpllims(1,2,i)
                    h_sq = h*h
                    if( h_sq > nyq_disks(i) ) cycle
                    k_max_h = int(sqrt(real(nyq_disks(i)-h_sq)))
                    k_lo = max(fpllims(2,1,i),-k_max_h)
                    k_hi = min(fpllims(2,2,i), k_max_h)
                    hp   = h*pf
                    do k = k_lo, k_hi
                        if( h_sq + k*k <= self%sq_rim ) cycle
                        kp = k*pf
                        if( kp <= 0 )then
                            cmplx_raw = conjg(fpls(i)%transfer_plane(hp,kp))*fpls(i)%cmplx_plane(hp,kp)
                        else
                            cmplx_raw = conjg(conjg(fpls(i)%transfer_plane(-hp,-kp))*fpls(i)%cmplx_plane(-hp,-kp))
                        endif
                        if( abs(real(cmplx_raw))+abs(aimag(cmplx_raw)) <= TINY ) cycle
                        loc = real(pf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. self%win_wraps(i0) ) cycle
                        vals = dsc(:,i) * (pf2 * cmplx_raw)
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_rhs_wrap(self, i0, w, vals, racc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        deallocate(rotmats, dsc, fpllims, nyq_disks)
    end subroutine accumulate_rhs

    pure logical function win_wraps( self, i0 )
        class(flex_pcg_t), intent(in) :: self
        integer,           intent(in) :: i0(3)
        win_wraps = any(i0 < self%wlims(1)) .or. any(i0 + self%wdim - 1 > self%wlims(2))
    end function win_wraps

    !> all pairs of one sample through one non-wrapping KB window; the pair index is the leading
    !! dimension so the innermost update is a contiguous vector operation
    !> stencil deposit into the band-list accumulator (nowrap: the stencil lies inside wlims)
    subroutine scatter_pairs_nowrap( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        real    :: wv
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_pairs_nowrap')
                    wv = w(di,dj,dk) * val
                    kacc(:,es) = kacc(:,es) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_nowrap
    pure subroutine scatter_pairs_nowrap_dense( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        real    :: wv
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    wv = w(di,dj,dk) * val
                    kacc(:,ih,ik,im) = kacc(:,ih,ik,im) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_nowrap_dense

    subroutine scatter_pairs_wrap( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        real    :: wv
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_pairs_wrap')
                    wv = w(di,dj,dk) * val
                    kacc(:,es) = kacc(:,es) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_wrap
    pure subroutine scatter_pairs_wrap_dense( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        real    :: wv
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    wv = w(di,dj,dk) * val
                    kacc(:,ih,ik,im) = kacc(:,ih,ik,im) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_wrap_dense

    subroutine scatter_rhs_nowrap( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_rhs_nowrap')
                    racc(:,es) = racc(:,es) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_nowrap
    pure subroutine scatter_rhs_nowrap_dense( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    racc(:,ih,ik,im) = racc(:,ih,ik,im) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_nowrap_dense

    subroutine scatter_rhs_wrap( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:)
        integer :: di, dj, dk, ih, ik, im, es
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    es = self%emap(ih,ik,im)
                    if( es == 0 ) THROW_HARD('stencil point outside the band list (plane beyond nband); scatter_rhs_wrap')
                    racc(:,es) = racc(:,es) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_wrap
    pure subroutine scatter_rhs_wrap_dense( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    racc(:,ih,ik,im) = racc(:,ih,ik,im) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_wrap_dense

    !> fold a full-range accumulator into the packed (h >= 0) real set, adding (the accumulator is freed);
    !! h runs to the 2x lattice's Nyquist plane (its accumulator column is the wrapped -Nyquist one)
    !> fold the band-list accumulator into the packed kernels: the same (m, k, h >= 0) visits as the
    !! dense fold, replayed from eorder/efold, so the sums are bitwise those of the dense arrays
    subroutine fold_accum( self, kacc, kpk )
        class(flex_pcg_t),  intent(inout) :: self
        real, allocatable,  intent(inout) :: kacc(:,:)
        real,               intent(inout) :: kpk(:,:)
        integer :: t
        if( .not. allocated(kacc) ) return
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; fold_accum')
        if( size(kacc,2) /= self%nexp .or. size(kpk,2) /= self%npk ) THROW_HARD('arrays are not on the band lists; fold_accum')
        do t = 1, self%nfold
            kpk(:,self%efold(t)) = kpk(:,self%efold(t)) + kacc(:,self%eorder(t))
        end do
        deallocate(kacc)
    end subroutine fold_accum
    subroutine fold_accum_dense( self, kacc, kpk )
        type(flex_pcg_t),   intent(inout) :: self
        real, allocatable,  intent(inout) :: kacc(:,:,:,:)
        real,               intent(inout) :: kpk(:,:,:,:)
        integer :: h, hh, k, m, phys(3), ih, ik, im
        if( .not. allocated(kacc) ) return
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; fold_accum')
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,phys,ih,ik,im) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, -self%lims3(1,1)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    ih   = hh - self%lims3(1,1) + 1
                    ik   = k  - self%lims3(2,1) + 1
                    im   = m  - self%lims3(3,1) + 1
                    kpk(:,phys(1),phys(2),phys(3)) = kpk(:,phys(1),phys(2),phys(3)) + kacc(:,ih,ik,im)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(kacc)
    end subroutine fold_accum_dense

    !> the same fold for the complex right-hand-side accumulator
    subroutine fold_rhs( self, racc, rpk )
        class(flex_pcg_t),     intent(inout) :: self
        complex, allocatable,  intent(inout) :: racc(:,:)
        complex,               intent(inout) :: rpk(:,:)
        integer :: t
        if( .not. allocated(racc) ) return
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; fold_rhs')
        if( size(racc,2) /= self%nexp .or. size(rpk,2) /= self%npk ) THROW_HARD('arrays are not on the band lists; fold_rhs')
        do t = 1, self%nfold
            rpk(:,self%efold(t)) = rpk(:,self%efold(t)) + racc(:,self%eorder(t))
        end do
        deallocate(racc)
    end subroutine fold_rhs
    subroutine fold_rhs_dense( self, racc, rpk )
        type(flex_pcg_t),      intent(inout) :: self
        complex, allocatable,  intent(inout) :: racc(:,:,:,:)
        complex,               intent(inout) :: rpk(:,:,:,:)
        integer :: h, hh, k, m, phys(3), ih, ik, im
        if( .not. allocated(racc) ) return
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; fold_rhs')
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,phys,ih,ik,im) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, -self%lims3(1,1)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    ih   = hh - self%lims3(1,1) + 1
                    ik   = k  - self%lims3(2,1) + 1
                    im   = m  - self%lims3(3,1) + 1
                    rpk(:,phys(1),phys(2),phys(3)) = rpk(:,phys(1),phys(2),phys(3)) + racc(:,ih,ik,im)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(racc)
    end subroutine fold_rhs_dense

    !> packed sums -> operator kernels: inverse transform on the 2x lattice, deposition envelope divided
    !! out, forward transform, real part, scaled to the gridding density (kernel_scale)
    subroutine finalize( self, kpk )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: kpk(:,:)
        real, parameter :: EPS_D = 1.0e-8
        real,    pointer     :: rp(:,:,:)
        integer, parameter :: IBLK = 16
        real,    allocatable :: kt(:,:,:,:)   !< (cdim, npairs): each pair's kernel, written contiguously by its thread
        integer :: ipair, i, j, k, tid, t, ib, i2
        integer(timer_int_kind) :: t_fold
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; finalize')
        if( size(kpk,2) /= self%npk ) THROW_HARD('packed kernel set is not on the band list; finalize')
        t_fold = tic()
        ! khat is pair-leading for the solve's per-voxel k x k multiply. Written pair by pair from the
        ! pair loop, every value landed in its own cache line shared between threads (a fold was memory
        ! bound at ~85 s for 528 pairs on 48 threads). Each thread now writes its pairs as contiguous
        ! slabs and one blocked transpose builds khat: same values, different store order.
        if( allocated(self%khat) )then
            if( size(self%khat,1) /= self%npairs .or. size(self%khat,2) /= self%cdim(1) ) deallocate(self%khat)
        endif
        if( .not. allocated(self%khat) ) allocate(self%khat(self%npairs, self%cdim(1), self%cdim(2), self%cdim(3)))
        allocate(kt(self%cdim(1), self%cdim(2), self%cdim(3), self%npairs))
        ! Scratch is allocated ONCE. The old loop did `tker = get_rmat()` and `ctmp = get_cmat()`
        ! per pair, reallocating a full padded volume 2 x npairs times per M-step (npairs = 55 at
        ! ncomp = 10, boxpd = 2*box_crop), which is what put ~26% of the master's CPU time in the
        ! kernel. Everything else here is memory-bound rather than arithmetic: kpk and khat are
        ! pair-LEADING -- contiguous for the solve's per-voxel k x k matrix, but strided by npairs
        ! in this loop, so each element of the gather and the scatter is its own cache line. Both
        ! are threaded now; the FFT plans are already threaded at this size.
        call self%ensure_pool()
        ! Parallel OVER PAIRS. Each pair needs two transforms on the padded lattice and they are
        ! independent, so with npairs = 55 at rank 10 this is 110 transforms per call and the call
        ! happens once per half set -- comparable to the whole 8-iteration solve. It was a serial
        ! loop with only the inner voxel passes threaded, which left the dominant term on one core.
        ! NOTE: this cost is NOT inside solve()'s reported seconds; it is master-tail time.
        !$omp parallel do default(shared) private(ipair,i,j,k,t,tid,rp) schedule(static)
        do ipair = 1, self%npairs
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_acc(:,:,:,tid) = cmplx(0.,0.)
            do t = 1, self%npk
                self%tp_acc(self%pijk(1,t),self%pijk(2,t),self%pijk(3,t),tid) = cmplx(kpk(ipair,t), 0.)
            end do
            call self%wpool(tid)%set_cmat(self%tp_acc(:,:,:,tid))
            call self%wpool(tid)%ifft()
            ! divide the deposition envelope out IN PLACE on the image's own buffer: ifft leaves
            ! the image in real space, so the get_rmat / set_rmat round trip is pure overhead
            call self%wpool(tid)%get_rmat_ptr(rp)
            do k = 1, self%boxpd
                do j = 1, self%boxpd
                    do i = 1, self%boxpd
                        if( abs(self%dep(i,j,k)) > EPS_D )then
                            rp(i,j,k) = rp(i,j,k) / self%dep(i,j,k)
                        else
                            rp(i,j,k) = 0.0
                        endif
                    end do
                end do
            end do
            call self%wpool(tid)%fft()
            call self%wpool(tid)%get_cmat_sub(self%tp_acc(:,:,:,tid))
            do k = 1, self%cdim(3)
                do j = 1, self%cdim(2)
                    do i = 1, self%cdim(1)
                        kt(i,j,k,ipair) = real(self%tp_acc(i,j,k,tid)) * FLEX_PCG_KSCALE
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        ! blocked transpose (cdim, npairs) -> (npairs, cdim): per (j,k) row, IBLK consecutive i share
        ! one cache line of every pair's slab, and the npairs values of one voxel are written contiguously
        !$omp parallel do default(shared) private(k,j,ib,i2,ipair) schedule(static) collapse(2)
        do k = 1, self%cdim(3)
            do j = 1, self%cdim(2)
                do ib = 1, self%cdim(1), IBLK
                    do i2 = ib, min(ib+IBLK-1, self%cdim(1))
                        do ipair = 1, self%npairs
                            self%khat(ipair,i2,j,k) = kt(i2,j,k,ipair)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(kt)
        self%l_kernel = .true.
        ! This cost sits OUTSIDE solve()'s reported seconds, which is why the master tail looked
        ! unaccounted for: 2*npairs padded transforms per call, once per half set.
        write(logfhandle,'(A,I0,A,I0,A,F8.1)') '>>> FLEX_PCA PCG KERNEL FOLD: npairs=', &
            &self%npairs, '  threads=', self%nthr_pool, '  seconds=', real(toc(t_fold))
        call flush(logfhandle)
    end subroutine finalize

    !> packed right-hand sides -> b = P Pi S^H y on the native lattice: inverse transform on the 2x
    !! lattice, central crop, deposition envelope divided out, Nyquist-ball band limit, support
    subroutine finalize_rhs( self, rpk, b )
        class(flex_pcg_t), intent(inout) :: self
        complex,           intent(in)    :: rpk(:,:)
        real,              intent(out)   :: b(:,:,:,:)
        real, pointer :: rp(:,:,:)
        integer :: q, i, j, k, tid, off, t
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; finalize_rhs')
        if( size(rpk,2) /= self%npk ) THROW_HARD('packed rhs set is not on the band list; finalize_rhs')
        call self%ensure_pool()
        off = (self%boxpd - self%box)/2
        ! parallel over components: independent transforms, per-thread images and scratch. The crop
        ! is inlined against the image buffer (offset (boxpd-box)/2, matching center_crop_real3d)
        ! and the band limit takes the thread's own native image.
        !$omp parallel do default(shared) private(q,i,j,k,t,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_acc(:,:,:,tid) = cmplx(0.,0.)
            do t = 1, self%npk
                self%tp_acc(self%pijk(1,t),self%pijk(2,t),self%pijk(3,t),tid) = rpk(q,t)
            end do
            call self%wpool(tid)%set_cmat(self%tp_acc(:,:,:,tid))
            call self%wpool(tid)%ifft()
            call self%wpool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = &
                &rp(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box) / self%dep_c
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            b(:,:,:,q) = self%tp_work(:,:,:,tid)
        end do
        !$omp end parallel do
    end subroutine finalize_rhs

    subroutine set_ridge( self, ridge )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: ridge(:,:)
        if( size(ridge,1) < self%ncomp ) THROW_HARD('ridge rank mismatch; set_ridge')
        if( allocated(self%ridge) ) deallocate(self%ridge)
        allocate(self%ridge(self%ncomp, size(ridge,2)), source=ridge(1:self%ncomp,:))
        self%l_ridge = .true.
    end subroutine set_ridge

    subroutine clear_ridge( self )
        class(flex_pcg_t), intent(inout) :: self
        if( allocated(self%ridge) ) deallocate(self%ridge)
        self%l_ridge = .false.
    end subroutine clear_ridge

    ! ---------------- operator, preconditioner, solve ----------------

    !> Nyquist-ball band limit on the native lattice (the shell rule of solve_coupled_basis_exp)
    !> allocate the per-thread transform pool once, OUTSIDE any parallel region (FFTW plan
    !! creation is not thread safe; image%new guards it with omp critical, but building the pool
    !! serially is both correct and cheaper). Every pooled image takes wthreads=.false. so that
    !! concurrent transforms never touch FFTW's thread pool.
    subroutine ensure_pool( self )
        class(flex_pcg_t), intent(inout) :: self
        integer :: nthr, t
        if( self%nthr_pool > 0 ) return
        nthr = 1
        !$ nthr = omp_get_max_threads()
        ! width = nthr, not min(nthr,ncomp): finalize folds npairs = ncomp*(ncomp+1)/2 kernels
        ! (55 at rank 10) and can use every thread, while apply_operator uses only the first ncomp
        nthr = max(1, nthr)
        allocate(self%wpool(nthr), self%npool(nthr))
        do t = 1, nthr
            call self%wpool(t)%new([self%boxpd,self%boxpd,self%boxpd], self%smpd, .false.)
            call self%npool(t)%new([self%box,self%box,self%box],       self%smpd, .false.)
        end do
        allocate(self%tp_work(self%box,self%box,self%box,nthr))
        allocate(self%tp_emb(self%boxpd,self%boxpd,self%boxpd,nthr))
        allocate(self%tp_acc(self%cdim(1),self%cdim(2),self%cdim(3),nthr))
        self%nthr_pool = nthr
    end subroutine ensure_pool

    !> band limit against a CALLER-SUPPLIED image, so it can run inside a parallel region over
    !! components. No internal OpenMP: the parallelism is the component loop outside. The shared-image
    !! variant below keeps its own threading for serial callers.
    subroutine bandlimit_img( self, work, img )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(inout) :: work(:,:,:)
        class(image),      intent(inout) :: img
        real, pointer :: rp(:,:,:)
        integer :: lims(3,2), h, k, m, phys(3)
        call img%set_rmat(work, .false.)
        call img%fft()
        lims = img%loop_lims(2)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    if( nint(sqrt(real(h*h + k*k + m*m))) > self%Rnat )then
                        phys = img%comp_addr_phys(h,k,m)
                        call img%set_cmat_at(phys(1),phys(2),phys(3), cmplx(0.,0.))
                    endif
                end do
            end do
        end do
        call img%ifft()
        call img%get_rmat_ptr(rp)
        work = rp(1:self%box,1:self%box,1:self%box)
    end subroutine bandlimit_img

    subroutine bandlimit( self, work )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(inout) :: work(:,:,:)
        real, pointer :: rp(:,:,:)
        integer :: lims(3,2), h, k, m, phys(3)
        call self%nimg%set_rmat(work, .false.)
        call self%nimg%fft()
        lims = self%nimg%loop_lims(2)
        !$omp parallel do collapse(2) default(shared) private(h,k,m,phys) schedule(static)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    if( nint(sqrt(real(h*h + k*k + m*m))) > self%Rnat )then
                        phys = self%nimg%comp_addr_phys(h,k,m)
                        call self%nimg%set_cmat_at(phys(1),phys(2),phys(3), cmplx(0.,0.))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        call self%nimg%ifft()
        call self%nimg%get_rmat_ptr(rp)
        work = rp(1:self%box,1:self%box,1:self%box)
    end subroutine bandlimit

    !> B x = P Pi T Pi P x (+ shell ridge): pair kernels on the 2x lattice, band-limited to the Nyquist
    !! ball, bracketed by the hard support
    subroutine apply_operator( self, x, hx )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: x(:,:,:,:)
        real,              intent(out)   :: hx(:,:,:,:)
        real,    pointer :: rp(:,:,:)
        complex :: acc
        integer :: q, r, i, j, k, sh, h, kk, m, phys(3), nsh, off, tid
        if( .not. self%l_kernel ) THROW_HARD('kernels are not finalized; apply_operator')
        ! Scratch is persistent. This routine runs once per CG iteration and used to allocate and
        ! free, per call, cq (ncomp x cdim complex = 681 MB at box_crop=128 / ncomp=10), acc, work,
        ! a get_cmat result (134 MB) per component and a get_rmat + center_crop pair (75 MB) per
        ! component: roughly 2.8 GB of mmap traffic per operator application. Measured on cnga1 at
        ! box_crop=128 the solve ran at 1.05 of 16 cores with 31% of CPU time in the kernel, i.e.
        ! it was bound by page faults, not arithmetic. center_embed/center_crop are inlined here
        ! against a persistent buffer for the same reason; both use offset = (boxpd-box)/2 and zero
        ! elsewhere, matching simple_cartesian_fourier exactly.
        if( .not. allocated(self%op_cq) ) &
            &allocate(self%op_cq(self%cdim(1),self%cdim(2),self%cdim(3),self%ncomp))
        off = (self%boxpd - self%box)/2
        call self%ensure_pool()
        ! Both component loops are parallel OVER COMPONENTS, each thread owning a padded image, a
        ! native image and its scratch. The transforms dominate: per CG iteration this is 2*ncomp
        ! on the boxpd lattice plus 4*ncomp on the native lattice (bandlimit does two each, twice
        ! per component), all of them independent. The previous code ran them serially through one
        ! shared self%wimg while threading only the k x k voxel loop, which left the dominant term
        ! single-threaded. Threading inside FFTW instead was measured to be a net harm.
        !$omp parallel do default(shared) private(q,tid) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_work(:,:,:,tid) = x(:,:,:,q)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            self%tp_emb(:,:,:,tid) = 0.0
            self%tp_emb(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box, tid) = &
                &self%tp_work(:,:,:,tid)
            call self%wpool(tid)%set_rmat(self%tp_emb(:,:,:,tid), .false.)
            call self%wpool(tid)%fft()
            call self%wpool(tid)%get_cmat_sub(self%op_cq(:,:,:,q))
        end do
        !$omp end parallel do
        ! Coupled multiply, VOXEL-OUTER. Component-outer made each of the ncomp threads stream the
        ! whole khat array independently -- 1.86 GB at rank 10, so ~18.6 GB of reads per call -- and
        ! used only ncomp of nthr threads. Voxel-outer reads khat's contiguous npairs triangle once
        ! per voxel, produces every component's output from it, and collapses over the lattice so
        ! all threads work. Bitwise-exact: each output still sums r ascending.
        if( .not. allocated(self%op_hq) ) &
            &allocate(self%op_hq(self%cdim(1),self%cdim(2),self%cdim(3),self%ncomp))
        !$omp parallel do collapse(3) default(shared) private(i,j,k,q,r,acc) schedule(static)
        do k = 1, self%cdim(3)
            do j = 1, self%cdim(2)
                do i = 1, self%cdim(1)
                    do q = 1, self%ncomp
                        acc = cmplx(0.,0.)
                        do r = 1, self%ncomp
                            acc = acc + self%khat(self%pidx(q,r),i,j,k) * self%op_cq(i,j,k,r)
                        end do
                        self%op_hq(i,j,k,q) = acc
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(q,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            call self%wpool(tid)%set_cmat(self%op_hq(:,:,:,q))
            call self%wpool(tid)%ifft()
            call self%wpool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = &
                &rp(off+1:off+self%box, off+1:off+self%box, off+1:off+self%box)
            call self%bandlimit_img(self%tp_work(:,:,:,tid), self%npool(tid))
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            hx(:,:,:,q) = self%tp_work(:,:,:,tid)
            ! Tikhonov term on the support (lam is derived with the floors; zero before prep_floor)
            if( allocated(self%lam) )then
                if( self%lam(q) > 0.0 )then
                    self%tp_work(:,:,:,tid) = x(:,:,:,q)
                    if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                    hx(:,:,:,q) = hx(:,:,:,q) + self%lam(q) * self%tp_work(:,:,:,tid)
                endif
            endif
        end do
        !$omp end parallel do
        ! Fourier shell ridge on the native lattice (the cross-FSC precision, optional)
        if( self%l_ridge )then
            nsh = size(self%ridge,2)
            !$omp parallel do default(shared) private(q,tid,h,kk,m,sh,phys,rp) schedule(static)
            do q = 1, self%ncomp
                tid = 1
                !$ tid = omp_get_thread_num() + 1
                self%tp_work(:,:,:,tid) = x(:,:,:,q)
                if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                call self%npool(tid)%set_rmat(self%tp_work(:,:,:,tid), .false.)
                call self%npool(tid)%fft()
                block
                    integer :: lims(3,2)
                    complex :: cval
                    lims = self%npool(tid)%loop_lims(2)
                    ! no nested OpenMP: the parallelism is the component loop outside
                    do m = lims(3,1), lims(3,2)
                        do kk = lims(2,1), lims(2,2)
                            do h = lims(1,1), lims(1,2)
                                sh   = nint(sqrt(real(h*h + kk*kk + m*m)))
                                phys = self%npool(tid)%comp_addr_phys(h,kk,m)
                                if( sh < 1 .or. sh > nsh .or. sh > self%Rnat )then
                                    cval = cmplx(0.,0.)
                                else
                                    cval = self%npool(tid)%get_cmat_at(phys(1),phys(2),phys(3)) &
                                        &* self%ridge(q,sh)
                                endif
                                call self%npool(tid)%set_cmat_at(phys(1),phys(2),phys(3), cval)
                            end do
                        end do
                    end do
                end block
                call self%npool(tid)%ifft()
                call self%npool(tid)%get_rmat_ptr(rp)
                self%tp_work(:,:,:,tid) = rp(1:self%box,1:self%box,1:self%box)
                if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
                hx(:,:,:,q) = hx(:,:,:,q) + self%tp_work(:,:,:,tid)
            end do
            !$omp end parallel do
        endif
    end subroutine apply_operator

    !> shell-relative floors of the preconditioner density: a fraction of the shell mean of each
    !! component's diagonal pair, so the per-voxel divide stays bounded on unsampled voxels
    subroutine prep_floor( self, rho, rho_lb )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: rho(:,:,:,:)
        integer,           intent(in)    :: rho_lb(3)
        real(dp), allocatable :: ssum(:,:)
        integer,  allocatable :: scnt(:)
        integer :: lims(3,2), h, k, m, sh, q, ih, ik, im
        if( size(rho,1) /= self%npairs ) THROW_HARD('preconditioner density must be the full packed triangle; prep_floor')
        if( allocated(self%rhofl) ) deallocate(self%rhofl)
        allocate(self%rhofl(self%ncomp, 0:self%Rnat), source=0.0)
        allocate(ssum(self%ncomp, 0:self%Rnat), source=0.0_dp)
        allocate(scnt(0:self%Rnat), source=0)
        lims = self%nimg%loop_lims(2)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > self%Rnat ) cycle
                    ih = h - rho_lb(1) + 1; ik = k - rho_lb(2) + 1; im = m - rho_lb(3) + 1
                    scnt(sh) = scnt(sh) + 1
                    do q = 1, self%ncomp
                        ssum(q,sh) = ssum(q,sh) + real(max(0.0, rho(pair_index(q,q),ih,ik,im)),dp)
                    end do
                end do
            end do
        end do
        do sh = 0, self%Rnat
            if( scnt(sh) < 1 ) cycle
            self%rhofl(:,sh) = FLEX_PCG_RHO_FLOOR_FRAC * real(ssum(:,sh) / real(scnt(sh),dp))
        end do
        ! Tikhonov term: lam_rel x the mean diagonal density over the low band (shells 1..max(4,Rnat/4)),
        ! the analogue of reconstructor_pcg's data scale
        if( allocated(self%lam) ) deallocate(self%lam)
        allocate(self%lam(self%ncomp), source=0.0)
        if( self%lam_rel > 0.0 )then
            block
                integer  :: shi
                real(dp) :: lsum(self%ncomp)
                integer  :: lcnt
                shi  = max(4, self%Rnat/4)
                lsum = 0.0_dp; lcnt = 0
                do sh = 1, min(shi, self%Rnat)
                    lsum = lsum + ssum(:,sh)
                    lcnt = lcnt + scnt(sh)
                end do
                if( lcnt > 0 ) self%lam = self%lam_rel * real(lsum / real(lcnt,dp))
            end block
        endif
        deallocate(ssum, scnt)
        self%l_floor = .true.
    end subroutine prep_floor

    !> z = P G P r: the per-voxel coupled divide on the gridding density (the Fourier-diagonal Jacobi
    !! inverse of T), diagonal pairs floored per shell, bracketed by the support. rho is the packed
    !! expanded-lattice density of the coupled insert with lower bounds rho_lb on the lattice axes.
    subroutine apply_precond( self, r, rho, rho_lb, z )
        class(flex_pcg_t), intent(inout) :: self
        real,              intent(in)    :: r(:,:,:,:)
        real,              intent(in)    :: rho(:,:,:,:)
        integer,           intent(in)    :: rho_lb(3)
        real,              intent(out)   :: z(:,:,:,:)
        real,    pointer     :: rp(:,:,:)
        complex(dp) :: rhs(self%ncomp), sol(self%ncomp)
        real(dp)    :: amat(self%ncomp,self%ncomp), diag_sum, ridge, denom
        integer :: lims(3,2), h, k, m, sh, q, rr, ih, ik, im, phys(3), flag, tid
        if( .not. self%l_floor ) call self%prep_floor(rho, rho_lb)
        call self%ensure_pool()
        if( .not. allocated(self%pc_cq) )then
            block
                integer :: ndim(3)
                ndim = self%nimg%get_array_shape()
                allocate(self%pc_cq(ndim(1),ndim(2),ndim(3),self%ncomp))
            end block
        endif
        ! forward transforms are independent per component: one per thread on the pool, and the
        ! accumulator is persistent (it was ~42 MB allocated and freed on every CG iteration)
        !$omp parallel do default(shared) private(q,tid) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            self%tp_work(:,:,:,tid) = r(:,:,:,q)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            call self%npool(tid)%set_rmat(self%tp_work(:,:,:,tid), .false.)
            call self%npool(tid)%fft()
            call self%npool(tid)%get_cmat_sub(self%pc_cq(:,:,:,q))
        end do
        !$omp end parallel do
        lims = self%nimg%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) schedule(static) &
        !$omp private(h,k,m,sh,q,rr,ih,ik,im,phys,amat,rhs,sol,diag_sum,ridge,denom,flag)
        do m = lims(3,1), lims(3,2)
            do k = lims(2,1), lims(2,2)
                do h = lims(1,1), lims(1,2)
                    phys = self%nimg%comp_addr_phys(h,k,m)
                    sh   = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > self%Rnat )then
                        self%pc_cq(phys(1),phys(2),phys(3),:) = cmplx(0.,0.)
                        cycle
                    endif
                    ih = h - rho_lb(1) + 1; ik = k - rho_lb(2) + 1; im = m - rho_lb(3) + 1
                    diag_sum = 0.d0
                    do q = 1, self%ncomp
                        do rr = q, self%ncomp
                            amat(q,rr) = real(rho(pair_index(q,rr),ih,ik,im), dp)
                            amat(rr,q) = amat(q,rr)
                        end do
                        amat(q,q)  = max(amat(q,q), real(self%rhofl(q,sh),dp)) + real(self%lam(q),dp)
                        diag_sum   = diag_sum + amat(q,q)
                        rhs(q)     = cmplx(self%pc_cq(phys(1),phys(2),phys(3),q), kind=dp)
                    end do
                    ridge = FLEX_PCG_RIDGE_REL * diag_sum / real(self%ncomp, dp)
                    do q = 1, self%ncomp
                        amat(q,q) = amat(q,q) + ridge
                    end do
                    call chol_solve_complex(amat, rhs, sol, self%ncomp, flag)
                    if( flag /= 0 )then
                        do q = 1, self%ncomp
                            denom = max(abs(amat(q,q)), ridge)
                            if( denom > DTINY )then
                                sol(q) = rhs(q) / denom
                            else
                                sol(q) = (0.d0, 0.d0)
                            endif
                        end do
                    endif
                    do q = 1, self%ncomp
                        self%pc_cq(phys(1),phys(2),phys(3),q) = cmplx(real(sol(q), sp), real(aimag(sol(q)), sp))
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        !$omp parallel do default(shared) private(q,tid,rp) schedule(static)
        do q = 1, self%ncomp
            tid = 1
            !$ tid = omp_get_thread_num() + 1
            call self%npool(tid)%set_cmat(self%pc_cq(:,:,:,q))
            call self%npool(tid)%ifft()
            call self%npool(tid)%get_rmat_ptr(rp)
            self%tp_work(:,:,:,tid) = rp(1:self%box,1:self%box,1:self%box)
            if( self%l_mask ) self%tp_work(:,:,:,tid) = self%tp_work(:,:,:,tid) * self%mask
            z(:,:,:,q) = self%tp_work(:,:,:,tid)
        end do
        !$omp end parallel do
    end subroutine apply_precond

    !> real SPD Cholesky solve with a complex right-hand side (as solve_real_spd_complex of the insert)
    pure subroutine chol_solve_complex( amat_in, rhs, sol, n, flag )
        integer,     intent(in)  :: n
        real(dp),    intent(in)  :: amat_in(n,n)
        complex(dp), intent(in)  :: rhs(n)
        complex(dp), intent(out) :: sol(n)
        integer,     intent(out) :: flag
        real(dp) :: chol(n,n), yr(n), yi(n), xr(n), xi(n)
        real(dp) :: sumr, sumi, sumv, tol
        integer  :: i, j, l
        flag = 0
        sol  = (0.d0, 0.d0)
        chol = 0.d0
        tol  = max(DTINY, epsilon(1.d0) * max(1.d0, maxval(abs(amat_in))))
        do j = 1, n
            sumv = amat_in(j,j)
            do l = 1, j - 1
                sumv = sumv - chol(j,l) * chol(j,l)
            end do
            if( sumv <= tol )then
                flag = 1
                return
            endif
            chol(j,j) = sqrt(sumv)
            do i = j + 1, n
                sumv = amat_in(i,j)
                do l = 1, j - 1
                    sumv = sumv - chol(i,l) * chol(j,l)
                end do
                chol(i,j) = sumv / chol(j,j)
            end do
        end do
        do i = 1, n
            sumr = real(rhs(i), dp)
            sumi = aimag(rhs(i))
            do l = 1, i - 1
                sumr = sumr - chol(i,l) * yr(l)
                sumi = sumi - chol(i,l) * yi(l)
            end do
            yr(i) = sumr / chol(i,i)
            yi(i) = sumi / chol(i,i)
        end do
        do i = n, 1, -1
            sumr = yr(i)
            sumi = yi(i)
            do l = i + 1, n
                sumr = sumr - chol(l,i) * xr(l)
                sumi = sumi - chol(l,i) * xi(l)
            end do
            xr(i) = sumr / chol(i,i)
            xi(i) = sumi / chol(i,i)
        end do
        do i = 1, n
            sol(i) = cmplx(xr(i), xi(i), kind=dp)
        end do
    end subroutine chol_solve_complex

    pure function dot_all( self, a, b ) result( d )
        class(flex_pcg_t), intent(in) :: self
        real,              intent(in) :: a(:,:,:,:), b(:,:,:,:)
        real(dp) :: d
        d = sum(real(a,dp) * real(b,dp))
    end function dot_all

    !> preconditioned CG on B x = b from the initial guess x (warm start), true-residual bookkeeping,
    !! periodic residual replacement, dx/x diminishing-returns stop, one cold restart on indefiniteness
    subroutine cg_core( self, b, x, rho, rho_lb, maxits, rtol, outcome, tag )
        class(flex_pcg_t),        intent(inout) :: self
        real,                     intent(in)    :: b(:,:,:,:)
        real,                     intent(inout) :: x(:,:,:,:)
        real,                     intent(in)    :: rho(:,:,:,:)
        integer,                  intent(in)    :: rho_lb(3)
        integer,                  intent(in)    :: maxits
        real,                     intent(in)    :: rtol
        type(flex_pcg_outcome_t), intent(inout) :: outcome
        character(len=*),         intent(in)    :: tag
        real, allocatable :: r(:,:,:,:), p(:,:,:,:), hp(:,:,:,:), z(:,:,:,:)
        real(dp) :: rhoz, rho_new, alpha, beta, pHp, bnorm, rnorm, xnorm, dxnorm, dxx, hnorm, bh
        integer  :: iter, n_done, attempt, verbose
        logical  :: stop_rtol, stop_xtol, l_warm
        if( maxits > 0 .and. .not. self%l_kernel ) THROW_HARD('kernels are not finalized; cg_core')
        verbose = 0
        call cov_env_int('SIMPLE_COV_PCG_VERBOSE', verbose)
        call self%prep_floor(rho, rho_lb)
        allocate(r, mold=b); allocate(p, mold=b); allocate(hp, mold=b); allocate(z, mold=b)
        bnorm = sqrt(self%dot_all(b,b))
        outcome%rhs_norm   = real(bnorm)
        outcome%start_norm = real(sqrt(self%dot_all(x,x)))
        if( bnorm <= 0.0_dp ) THROW_HARD('zero right-hand side; nothing to solve')
        outcome%stop_reason = 'maxits'
        if( rtol <= 0.0 ) outcome%stop_reason = 'fixed_iterations'
        l_warm = outcome%start_norm > 0.0
        if( verbose > 0 ) call audit_terms
        rnorm  = 0.0_dp; dxx = 0.0_dp; n_done = 0
        do attempt = 1, 2
            if( l_warm )then
                call self%apply_operator(x, hp)
                r = b - hp
            else
                ! B*0 = 0, so the zero start's residual IS b and the operator application is pure
                ! waste -- one of nine at maxits=8, ~11% of a cold solve. The production solver
                ! states the same thing: "free: the zero start's residual is b itself, no operator
                ! application" (simple_reconstructor_pcg.f90, solve_accum).
                hp = 0.0
                r  = b
            endif
            rnorm = sqrt(self%dot_all(r,r))
            if( attempt == 1 )then
                outcome%initial_rel_residual = real(rnorm / bnorm)
                ! how the warm start relates to the data: corr(b, B x0) and the least-squares scale
                hnorm = sqrt(self%dot_all(hp,hp))
                bh    = self%dot_all(b,hp)
                outcome%start_corr  = real(bh / max(bnorm*hnorm, 1.0e-30_dp))
                outcome%start_scale = real(bh / max(hnorm*hnorm, 1.0e-30_dp))
                ! Apply that least-squares scale to the warm start. It is FREE: hp = B x0 is
                ! already in hand from the operator application above, so rescaling x and hp and
                ! recomputing the residual costs three vector passes and no transforms. Measured on
                ! cnga1 at box_crop=128, warm-starting from the gridding solution: corr(b,Bx0)=0.77
                ! but scale=0.52, i.e. the gridding iterate points the right way and is ~2x too
                ! large. Unscaled it gives an initial relative residual of 0.959; at the optimal
                ! scale that becomes sqrt(1-corr^2) = 0.634, which is where the COLD solve arrives
                ! only after 8 iterations and ~520 s. Skipped when the scale is not a finite
                ! positive number, and the attempt=1,2 cold restart below still guards curvature.
                if( l_warm .and. ieee_is_finite(real(outcome%start_scale,dp)) .and. &
                    &outcome%start_scale > 0.0 )then
                    x     = outcome%start_scale * x
                    hp    = outcome%start_scale * hp
                    r     = b - hp
                    rnorm = sqrt(self%dot_all(r,r))
                    outcome%start_norm           = real(sqrt(self%dot_all(x,x)))
                    outcome%initial_rel_residual = real(rnorm / bnorm)
                endif
                if( verbose > 0 ) write(logfhandle,'(A,A,A,ES10.3,A,F8.4,A,ES10.3,A,ES10.3,A,ES10.3)') '>>> ', tag, &
                    &' start: rel resid=', outcome%initial_rel_residual, '  corr(b,Bx0)=', outcome%start_corr, &
                    &'  scale <b,Bx0>/<Bx0,Bx0>=', outcome%start_scale, '  lam_rel=', self%lam_rel, &
                    &'  lam(1)=', self%lam(1)
            endif
            call self%apply_precond(r, rho, rho_lb, z)
            p    = z
            rhoz = self%dot_all(r,z)
            if( rhoz <= 0.0_dp ) THROW_HARD('non-positive initial dot(r,z); the preconditioner is not positive definite')
            n_done = 0
            dxx    = 0.0_dp
            do iter = 1, maxits
                call self%apply_operator(p, hp)
                pHp = self%dot_all(p,hp)
                if( .not. ieee_is_finite(pHp) .or. pHp <= 0.0_dp )then
                    if( l_warm .and. attempt == 1 )then
                        write(logfhandle,'(A,A,ES11.3,A,I0,A)') '>>> ', tag, real(pHp), &
                            &' curvature at iteration ', iter, ': warm start discarded, cold restart'
                        outcome%cold_restart_used = .true.
                        x = 0.0
                        l_warm = .false.
                        exit
                    endif
                    outcome%stop_reason = FLEX_PCG_STOP_INDEFINITE
                    outcome%converged   = .false.
                    n_done = iter - 1
                    exit
                endif
                alpha = rhoz / pHp
                x = x + real(alpha) * p
                r = r - real(alpha) * hp
                if( mod(iter, FLEX_PCG_RESID_REPLACE) == 0 )then
                    call self%apply_operator(x, hp)
                    r = b - hp
                endif
                n_done = iter
                rnorm  = sqrt(self%dot_all(r,r))
                xnorm  = sqrt(self%dot_all(x,x))
                dxnorm = abs(alpha) * sqrt(self%dot_all(p,p))
                dxx    = dxnorm / max(xnorm, epsilon(1.0_dp))
                if( outcome%iters_to_1e2 == 0 .and. rnorm / bnorm <= 1.0e-2_dp ) outcome%iters_to_1e2 = iter
                stop_rtol = rtol > 0.0 .and. rnorm / bnorm <= real(rtol,dp)
                stop_xtol = rtol > 0.0 .and. dxx <= real(FLEX_PCG_XTOL,dp)
                if( verbose > 0 ) write(logfhandle,'(A,A,A,I4,A,ES10.3,A,ES10.3,A,ES10.3)') '>>> ', tag, ' it', iter, &
                    &'  rel resid=', real(rnorm/bnorm), '  dx/x=', real(dxx), '  alpha=', real(alpha)
                if( stop_rtol )then
                    outcome%stop_reason = 'rtol'; outcome%converged = .true.; exit
                else if( stop_xtol )then
                    outcome%stop_reason = 'xtol'; outcome%converged = .true.; exit
                endif
                if( iter == maxits ) exit
                call self%apply_precond(r, rho, rho_lb, z)
                rho_new = self%dot_all(r,z)
                beta    = rho_new / rhoz
                p       = z + real(beta) * p
                rhoz    = rho_new
            end do
            if( trim(outcome%stop_reason) == FLEX_PCG_STOP_INDEFINITE ) exit
            if( .not. l_warm .and. attempt == 1 .and. outcome%cold_restart_used ) cycle
            exit
        end do
        outcome%iteration_count    = n_done
        outcome%final_rel_residual = real(rnorm / bnorm)
        outcome%final_rel_update   = real(dxx)
        ! a vanishing curvature means the search direction has reached the operator's null space (the
        ! unregularised masked problem is only positive semi-definite): the iterate so far is the answer
        if( trim(outcome%stop_reason) == FLEX_PCG_STOP_INDEFINITE ) write(logfhandle,'(A,A,A,I0,A)') '>>> ', tag, &
            &': curvature vanished after ', n_done, ' iterations (null space reached); returning the current iterate'
        if( .not. all(ieee_is_finite(x)) ) THROW_HARD(tag//': non-finite solution')
        deallocate(r, p, hp, z)

    contains

        !> How large are the two regularisers inside the operator, relative to the data term, on the
        !! right-hand side as a probe? lam is derived from the GRIDDING density while the data term is
        !! the doubled-lattice Toeplitz Gram, and the cross-FSC ridge arrives on the gridding scale too
        !! (add_invtausq2rho_coupled adds the same numbers to rho). This states the ratios instead of
        !! leaving the relative scaling to be inferred from the source.
        subroutine audit_terms
            real, allocatable :: h0(:,:,:,:), h1(:,:,:,:), lam_save(:)
            logical  :: ridge_save
            real(dp) :: n0, nl, nr
            allocate(h0, mold=b); allocate(h1, mold=b)
            ridge_save = self%l_ridge
            if( allocated(self%lam) ) lam_save = self%lam
            self%l_ridge = .false.
            if( allocated(self%lam) ) self%lam = 0.0
            call self%apply_operator(b, h0)                       ! data term alone
            n0 = sqrt(self%dot_all(h0,h0))
            if( allocated(lam_save) ) self%lam = lam_save
            call self%apply_operator(b, h1)                       ! data + Tikhonov
            h1 = h1 - h0
            nl = sqrt(self%dot_all(h1,h1))
            h1 = h1 + h0                                          ! restore data + Tikhonov
            self%l_ridge = ridge_save
            nr = 0.0_dp
            if( self%l_ridge )then
                call self%apply_operator(b, h0)                   ! data + Tikhonov + FSC prior
                h0 = h0 - h1
                nr = sqrt(self%dot_all(h0,h0))
            endif
            write(logfhandle,'(A,A,A,ES10.3,A,ES10.3,A,L1)') '>>> ', tag, &
                &' operator terms on the rhs probe: Tikhonov/data=', real(nl / max(n0,1.0e-30_dp)), &
                &'  FSC-prior/data=', real(nr / max(n0,1.0e-30_dp)), '  prior_installed=', self%l_ridge
            deallocate(h0, h1)
        end subroutine audit_terms

    end subroutine cg_core

    !> the M-step solve on one half. Budget 0 is the explicit gridding mode: the solution of the 1x
    !! numerators (solve_coupled_basis_exp, deapodized) is delivered untouched. Any positive budget runs
    !! CG on the 2x right-hand sides FROM ZERO, the production reconstruct3D_pcg convention. The gridding
    !! solution is still formed there, but only as the reference the diagnostic log compares against: it is
    !! this operator's block-Jacobi preconditioner applied to the numerators, so a cold CG's own first
    !! iterate already lands essentially on it, while starting there would make a budget of n mean
    !! "gridding plus n" and would carry an unprojected object into a projected solve.
    !! On return the reconstructors carry E*u on the expanded lattice for the unchanged tail.
    subroutine solve( self, Y, rho, rpk, maxits, rtol, outcome, tag )
        class(flex_pcg_t),     intent(inout) :: self
        type(reconstructor),   intent(inout) :: Y(:)
        real,                  intent(in)    :: rho(:,:,:,:)
        complex,               intent(in)    :: rpk(:,:)
        integer,               intent(in)    :: maxits
        real,                  intent(in)    :: rtol
        type(flex_pcg_outcome_t), intent(out) :: outcome
        character(len=*), optional, intent(in) :: tag
        real, allocatable :: b(:,:,:,:), x(:,:,:,:), xgrid(:,:,:,:)
        real, pointer     :: rp(:,:,:)
        integer  :: q, rho_lb(3), verbose, iwarm
        real(dp) :: gnorm, dnorm, gdot
        integer(timer_int_kind) :: t0
        character(len=:), allocatable :: ttag
        ttag = 'FLEX_PCA PCG MSTEP'
        if( present(tag) ) ttag = trim(tag)
        if( size(Y) /= self%ncomp ) THROW_HARD('reconstructor count differs from the solver rank; solve')
        t0 = tic()
        outcome%requested_maxits = maxits
        rho_lb = lbound(Y(1)%cmat_exp)
        verbose = 0
        call cov_env_int('SIMPLE_COV_PCG_VERBOSE', verbose)
        allocate(x(self%box,self%box,self%box,self%ncomp))
        ! The Y reconstructors are built by init_basis_reconstructor through the 3-argument
        ! image%new (simple_flex_pca_em_fit.f90:641), and image%new DEFAULTS wthreads to .true.
        ! (simple_image_core.f90:38) with plan_nthreads = nthr_glob. At box_crop >= 100 that means
        ! every Y(q) carries a 16-thread FFTW plan, and solve then runs 2*ncomp of those transforms
        ! in SERIAL loops (the ifft below and put_back's fft). That is the semaphore-bound
        ! configuration measured harmful on this path: with threaded plans driven one transform at a
        ! time, the master sat in os_sem_down inside libfftw3f_threads. Turning it off here rather
        ! than at the constructor confines the change to the PCG backend, leaving the gridding
        ! path's bit-identical-eigenvalue policy gate untouched. set_wthreads early-returns when the
        ! flag already matches, so only the first solve pays the plan rebuild.
        do q = 1, self%ncomp
            call Y(q)%set_wthreads(.false.)
        end do
        ! the gridding solution, physical units (the tail's inverse envelope undone here)
        call solve_coupled_basis_exp(Y, rho, self%ncomp)
        do q = 1, self%ncomp
            Y(q)%rho_exp = 1.0
            call Y(q)%compress_exp
            call Y(q)%ifft
            call Y(q)%get_rmat_ptr(rp)
            x(:,:,:,q) = rp(1:self%box,1:self%box,1:self%box) / self%env
        end do
        ! budget 0 ships the gridding solution untouched (the tail's soft mask is inside the hard support
        ! anyway, but the half-set FSC that sets the Wiener filter is taken on the unmasked halves and the
        ! support would inflate it); with iterations the solution lives on the support by construction
        if( maxits <= 0 )then
            outcome%stop_reason     = 'block_jacobi'
            outcome%iteration_count = 0
            outcome%start_norm      = real(sqrt(self%dot_all(x,x)))
            call put_back
            outcome%seconds = real(toc(t0))
            deallocate(x)
            return
        endif
        if( .not. self%l_kernel ) THROW_HARD('kernels are not finalized; solve')
        ! cold start (production convention). The zero iterate is trivially on the support, so no masking
        ! of the start is needed; the gridding solution is kept only when the diagnostic is asked for.
        if( verbose > 0 )then
            allocate(xgrid, source=x)
            ! compare on the SUPPORT: the CG solution is zero outside it while the gridding reference is a
            ! full-box object, and with a 2%-occupancy envelope an unmasked comparison reports orthogonality
            ! no matter how close the two solutions are where they are both defined
            if( self%l_mask )then
                do q = 1, self%ncomp
                    xgrid(:,:,:,q) = xgrid(:,:,:,q) * self%mask
                end do
            endif
            gnorm = sqrt(self%dot_all(xgrid,xgrid))
        endif
        ! Start iterate. The cold start is the production convention, and the zero iterate is
        ! trivially on the support while the gridding solution is a full-box object that has to be
        ! masked onto it first. But that gridding (block-Jacobi) solution has already been computed
        ! above and then discarded, cg_core already carries the warm-start machinery (start_corr,
        ! start_scale and the attempt=1,2 cold restart on indefinite curvature), and every box-128
        ! solve so far has stopped on maxits at a relative residual near 0.6 starting from 1.0. When
        ! the solve is truncated that far from convergence the initial iterate IS the budget, so
        ! this is the one lever that cuts the iteration COUNT instead of the cost per iteration.
        ! Opt-in until measured: SIMPLE_COV_PCG_WARM=1 starts from the masked gridding solution.
        iwarm = 1   ! the scaled warm start is the estimator of record (measured 2026-09-14/16); no switch
        if( iwarm > 0 )then
            if( self%l_mask )then
                do q = 1, self%ncomp
                    x(:,:,:,q) = x(:,:,:,q) * self%mask
                end do
            endif
        else
            x = 0.0
        endif
        allocate(b(self%box,self%box,self%box,self%ncomp))
        call self%finalize_rhs(rpk, b)
        call self%cg_core(b, x, rho, rho_lb, maxits, rtol, outcome, ttag)
        if( allocated(xgrid) )then
            gdot  = self%dot_all(x,xgrid)
            xgrid = x - xgrid
            dnorm = sqrt(self%dot_all(xgrid,xgrid))
            write(logfhandle,'(A,A,A,I0,A,ES10.3,A,F8.4)') '>>> ', ttag, ' cold solve vs the gridding &
                &reference: iters=', outcome%iteration_count, '  ||x-xg||/||xg||=', &
                &real(dnorm / max(gnorm, 1.0e-30_dp)), '  corr(x,xg)=', &
                &real(gdot / max(sqrt(self%dot_all(x,x))*gnorm, 1.0e-30_dp))
            deallocate(xgrid)
        endif
        call put_back
        outcome%seconds = real(toc(t0))
        deallocate(b, x)

    contains

        !> E*u back into the reconstructors' expanded lattice: the tail divides by E and ships u
        subroutine put_back
            integer :: qq
            do qq = 1, self%ncomp
                call Y(qq)%set_rmat(x(:,:,:,qq) * self%env, .false.)
                call Y(qq)%fft
                call Y(qq)%expand_exp
            end do
        end subroutine put_back

    end subroutine solve

    ! ---------------- envelope support (pcg_mskfile) ----------------

    !> pcg_mskfile resampled to a target lattice under the note's 4.2 contract: Fourier crop, ringing
    !! diagnostics (min/max/out-of-range fraction before clamping), clamp to [0,1], floor to exactly 0
    !! below FLEX_PCG_SUPPORT_FLOOR, occupancy of the box and of the mskdiam sphere on the native and
    !! the target lattice
    subroutine flex_pcg_support_volume( params, box_t, smpd_t, img, tag )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: box_t
        real,              intent(in)    :: smpd_t
        type(image),       intent(inout) :: img
        character(len=*),  intent(in)    :: tag
        type(image)   :: nat, rc
        real, pointer :: rp(:,:,:)
        real    :: vmin, vmax, frac_out, occ_box_n, occ_sph_n, occ_box_t, occ_sph_t
        integer :: ldim(3), ifoo, n
        if( .not. flex_mskfile_set(params) ) THROW_HARD('pcg_mskfile is not set; flex_pcg_support_volume')
        call find_ldim_nptcls(params%pcg_mskfile, ldim, ifoo)
        if( ldim(1) /= ldim(2) .or. ldim(1) /= ldim(3) ) THROW_HARD('pcg_mskfile must be a cubic volume')
        if( ldim(1) /= params%box ) THROW_HARD('pcg_mskfile box differs from the project box')
        ! native occupancy (the reference for the resampling diagnostic)
        call nat%new(ldim, params%smpd)
        call nat%read(params%pcg_mskfile)
        call nat%get_rmat_ptr(rp)
        call occupancy(rp, ldim(1), 0.5*params%mskdiam/params%smpd, occ_box_n, occ_sph_n)
        call nat%kill
        ! Fourier crop to the target lattice; the result is copied into a fresh image of the target box so
        ! that its buffers are exactly the target's (the in-place clip keeps the native allocation)
        call rc%read_and_crop(params%pcg_mskfile, params%smpd, box_t, smpd_t)
        if( rc%is_ft() ) call rc%ifft
        call rc%get_rmat_ptr(rp)
        n        = box_t**3
        vmin     = minval(rp(1:box_t,1:box_t,1:box_t))
        vmax     = maxval(rp(1:box_t,1:box_t,1:box_t))
        frac_out = real(count(rp(1:box_t,1:box_t,1:box_t) < 0.0 .or. rp(1:box_t,1:box_t,1:box_t) > 1.0)) / real(n)
        rp(1:box_t,1:box_t,1:box_t) = min(1.0, max(0.0, rp(1:box_t,1:box_t,1:box_t)))
        where( rp(1:box_t,1:box_t,1:box_t) < FLEX_PCG_SUPPORT_FLOOR ) rp(1:box_t,1:box_t,1:box_t) = 0.0
        if( .not. any(rp(1:box_t,1:box_t,1:box_t) > 0.0) ) THROW_HARD('pcg_mskfile is empty after resampling')
        call occupancy(rp, box_t, 0.5*params%mskdiam/smpd_t, occ_box_t, occ_sph_t)
        call img%new([box_t,box_t,box_t], smpd_t)
        call img%set_rmat(rp(1:box_t,1:box_t,1:box_t), .false.)
        call rc%kill
        write(logfhandle,'(A,A,A,I0,A,F7.3,A,F6.3,A,F6.3,A,F7.4)') '>>> FLEX_PCA PCG SUPPORT (', trim(tag), &
            &'): pcg_mskfile at box ', box_t, ' smpd ', smpd_t, ': before clamp min=', vmin, ' max=', vmax, &
            &' out-of-[0,1] fraction=', frac_out
        write(logfhandle,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') '>>> FLEX_PCA PCG SUPPORT occupancy (window>0): native box ', &
            &occ_box_n, ' / sphere ', occ_sph_n, '  ->  target box ', occ_box_t, ' / sphere ', occ_sph_t
        call flush(logfhandle)

    contains

        subroutine occupancy( v, nb, rad, f_box, f_sph )
            real,    intent(in)  :: v(:,:,:)
            integer, intent(in)  :: nb
            real,    intent(in)  :: rad
            real,    intent(out) :: f_box, f_sph
            integer :: i, j, k, c, nin, nsph
            real    :: r2, rad2
            c = nb/2 + 1; rad2 = rad*rad
            nin = 0; nsph = 0
            do k = 1, nb
                do j = 1, nb
                    do i = 1, nb
                        r2 = real((i-c)**2 + (j-c)**2 + (k-c)**2)
                        if( r2 <= rad2 )then
                            nsph = nsph + 1
                            if( v(i,j,k) > 0.0 ) nin = nin + 1
                        endif
                    end do
                end do
            end do
            f_box = real(count(v(1:nb,1:nb,1:nb) > 0.0)) / real(nb)**3
            f_sph = 0.0
            if( nsph > 0 ) f_sph = real(nin) / real(nsph)
        end subroutine occupancy

    end subroutine flex_pcg_support_volume

    !> whether pcg_mskfile was given (the string is unallocated when absent)
    logical function flex_mskfile_set( params )
        class(parameters), intent(in) :: params
        flex_mskfile_set = .false.
        if( params%pcg_mskfile%is_allocated() ) flex_mskfile_set = len_trim(params%pcg_mskfile%to_char()) > 0
    end function flex_mskfile_set

    !> idempotent: the envelope at the covariance box when pcg_mskfile is set on the PCG backend
    subroutine flex_env_init( params )
        class(parameters), intent(in) :: params
        if( l_flex_env_checked ) return
        l_flex_env_checked = .true.
        l_flex_env = trim(params%rec_backend) == 'pcg' .and. flex_mskfile_set(params)
        if( .not. l_flex_env ) return
        call flex_pcg_support_volume(params, params%box_crop, params%smpd_crop, flex_env_img, 'covariance box')
    end subroutine flex_env_init

    logical function flex_env_active()
        flex_env_active = l_flex_env_checked .and. l_flex_env
    end function flex_env_active

    !> number of complement (outside-the-envelope) components requested; 0 when off or no envelope
    integer function flex_complement_ncomp( params )
        class(parameters), intent(in) :: params
        character(len=32) :: envval
        integer :: ln, stat, ios, ival
        if( flex_cenv_ncomp < 0 )then
            flex_cenv_ncomp = 0
            call get_environment_variable('SIMPLE_COV_COMPLEMENT', envval, ln, stat)
            if( stat == 0 .and. ln > 0 )then
                read(envval, *, iostat=ios) ival
                if( ios == 0 ) flex_cenv_ncomp = max(0, ival)
            endif
        endif
        call flex_env_init(params)
        flex_complement_ncomp = 0
        if( l_flex_env ) flex_complement_ncomp = flex_cenv_ncomp
    end function flex_complement_ncomp

    logical function flex_complement_active( params )
        class(parameters), intent(in) :: params
        flex_complement_active = flex_complement_ncomp(params) > 0
    end function flex_complement_active

    !> the complement window: the soft sphere (msk_crop) minus the envelope, built once
    subroutine flex_cenv_ensure( params )
        class(parameters), intent(in) :: params
        real, pointer :: cp(:,:,:), wp(:,:,:)
        integer :: b
        real    :: occ
        if( flex_cenv_img%exists() ) return
        b = params%box_crop
        call flex_cenv_img%new([b,b,b], params%smpd_crop)
        call flex_cenv_img%get_rmat_ptr(cp)
        cp = 0.
        cp(1:b,1:b,1:b) = 1.
        if( params%msk_crop > TINY ) call flex_cenv_img%mask3D_soft(params%msk_crop, backgr=0.)
        call flex_cenv_img%get_rmat_ptr(cp)
        call flex_env_img%get_rmat_ptr(wp)
        cp(1:b,1:b,1:b) = cp(1:b,1:b,1:b) * max(0., 1. - wp(1:b,1:b,1:b))
        occ = sum(cp(1:b,1:b,1:b)) / real(b)**3
        write(logfhandle,'(A,I0,A,F7.4,A,F7.4,A)') '>>> FLEX_PCA COMPLEMENT BLOCK: ', flex_cenv_ncomp, &
            &' components on the sphere minus the envelope (occupancy ', occ, ' of the box; envelope ', &
            &sum(wp(1:b,1:b,1:b))/real(b)**3, ')'
        call flush(logfhandle)
    end subroutine flex_cenv_ensure

    !> block-aware window: blk=1 the envelope (or the sphere when no envelope), blk=2 the complement
    subroutine flex_window_apply_blk( img, params, blk )
        type(image),       intent(inout) :: img
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: blk
        real, pointer :: rp(:,:,:), cp(:,:,:)
        integer :: b
        if( blk /= 2 .or. .not. flex_complement_active(params) )then
            call flex_window_apply(img, params)
            return
        endif
        call flex_cenv_ensure(params)
        b = params%box_crop
        if( any(img%get_ldim() /= b) ) THROW_HARD('flex_window_apply_blk: volume is not at the covariance box')
        if( img%is_ft() ) call img%ifft
        call img%get_rmat_ptr(rp)
        call flex_cenv_img%get_rmat_ptr(cp)
        rp(1:b,1:b,1:b) = rp(1:b,1:b,1:b) * cp(1:b,1:b,1:b)
    end subroutine flex_window_apply_blk

    subroutine flex_window_apply_rec_blk( rec, params, blk )
        type(reconstructor), intent(inout) :: rec
        class(parameters),   intent(in)    :: params
        integer,             intent(in)    :: blk
        real, pointer :: rp(:,:,:), cp(:,:,:)
        integer :: b
        if( blk /= 2 .or. .not. flex_complement_active(params) )then
            call flex_window_apply_rec(rec, params)
            return
        endif
        call flex_cenv_ensure(params)
        b = params%box_crop
        if( any(rec%get_ldim() /= b) ) THROW_HARD('flex_window_apply_rec_blk: volume is not at the covariance box')
        if( rec%is_ft() ) call rec%ifft
        call rec%get_rmat_ptr(rp)
        call flex_cenv_img%get_rmat_ptr(cp)
        rp(1:b,1:b,1:b) = rp(1:b,1:b,1:b) * cp(1:b,1:b,1:b)
    end subroutine flex_window_apply_rec_blk

    !> the window every basis-shaped volume is multiplied by: the envelope when it is set, else the
    !! soft sphere of the gridding path (mask3D_soft at msk_crop). Typed entry points on purpose: a
    !! class(image) dummy segfaulted at -O3 on the first call of the M-step tail (gfortran codegen).
    subroutine flex_window_apply( img, params )
        type(image),       intent(inout) :: img   !< at the covariance box
        class(parameters), intent(in)    :: params
        real, pointer :: rp(:,:,:)
        call flex_env_init(params)
        if( l_flex_env )then
            if( any(img%get_ldim() /= params%box_crop) ) THROW_HARD('flex_window_apply: volume is not at the covariance box')
            if( img%is_ft() ) call img%ifft
            call img%get_rmat_ptr(rp)
            call window_product(rp, params%box_crop)
        else
            if( params%msk_crop > TINY ) call img%mask3D_soft(params%msk_crop, backgr=0.)
        endif
    end subroutine flex_window_apply

    !> the same for a reconstructor (the initial basis is built on one)
    subroutine flex_window_apply_rec( rec, params )
        type(reconstructor), intent(inout) :: rec
        class(parameters),   intent(in)    :: params
        real, pointer :: rp(:,:,:)
        call flex_env_init(params)
        if( l_flex_env )then
            if( any(rec%get_ldim() /= params%box_crop) ) THROW_HARD('flex_window_apply_rec: volume is not at the covariance box')
            if( rec%is_ft() ) call rec%ifft
            call rec%get_rmat_ptr(rp)
            call window_product(rp, params%box_crop)
        else
            if( params%msk_crop > TINY ) call rec%mask3D_soft(params%msk_crop, backgr=0.)
        endif
    end subroutine flex_window_apply_rec

    !> elementwise product with the envelope on the logical box (loops: no array temporary)
    subroutine window_product( rp, b )
        real, pointer, intent(inout) :: rp(:,:,:)
        integer,       intent(in)    :: b
        real, pointer :: wp(:,:,:)
        integer :: i, j, k
        if( .not. flex_env_img%exists() ) THROW_HARD('window_product: the envelope window is not loaded')
        call flex_env_img%get_rmat_ptr(wp)
        do k = 1, b
            do j = 1, b
                do i = 1, b
                    rp(i,j,k) = rp(i,j,k) * wp(i,j,k)
                end do
            end do
        end do
    end subroutine window_product

    !> the solver's support at the covariance box: the envelope when set, else the sphere
    subroutine flex_pcg_install_window( op, params )
        type(flex_pcg_t),  intent(inout) :: op
        class(parameters), intent(in)    :: params
        integer :: kfr_pcg(2)
        call flex_env_init(params)
        if( l_flex_env )then
            call op%set_window_volume(flex_env_img)
            write(logfhandle,'(A,F6.3)') '>>> FLEX_PCA PCG support: envelope window (pcg_mskfile), hard support fraction ', &
                &sum(op%mask)/real(size(op%mask))
        else
            if( params%msk_crop > TINY ) call op%set_window_sphere(params%msk_crop)
        endif
        ! the band the planes carry (projected_model_kfromto) plus two shells of margin: the
        ! accumulators and packed kernels of every operator live only on the points this band reaches
        kfr_pcg = projected_model_kfromto(params)
        call op%set_band(kfr_pcg(2) + 2)
    end subroutine flex_pcg_install_window

    ! ---------------- self-test ----------------

    !> Four checks on a random Gaussian test volume at the given box. (A) operator: E Pi T Pi E u through
    !! the doubled-coordinate kernels against the exact nonuniform-DFT Gram (the KB-interpolation Gram is
    !! itself 18-22% off and is not a reference); (B) right-hand side: the 2x deposit of exact samples
    !! against the exact adjoint; (C) box <= 32 only: a full preconditioned CG solve of T u = S^H y from
    !! central-slice samples of the volume, recovery of the volume on a spherical support: the clean,
    !! generously supported baseline by default; with sweep=.true. twelve solves over support size, sample
    !! noise and the Tikhonov term, every clean one held to the baseline criterion. (D) the band-list
    !! kernels and right-hand side against the dense fold, bitwise.
    subroutine test_flex_pcg_operator( box, nsamples, l_pass, passes, sweep )
        integer, intent(in)  :: box, nsamples
        logical, intent(out) :: l_pass
        logical, optional, intent(out) :: passes(4)   !< (A) operator, (B) rhs, (C) solve, (D) band lists
        logical, optional, intent(in)  :: sweep       !< (C) over all twelve settings (default: the baseline)
        type(flex_pcg_t) :: op
        type(flex_pcg_outcome_t) :: out
        real,    allocatable :: kacc(:,:), kpk(:,:), u(:,:,:,:), hu(:,:,:,:), eu4(:,:,:,:)
        real,    allocatable :: kacc4(:,:,:,:), kpk4(:,:,:,:)
        complex, allocatable :: racc4(:,:,:,:), rpk4(:,:,:,:)
        real,    allocatable :: locs(:,:), wts(:), b4(:,:,:,:), x4(:,:,:,:), rho_t(:,:,:,:), ut(:,:,:)
        complex, allocatable :: racc(:,:), rpk(:,:)
        real(dp), allocatable :: be(:,:,:), eu(:,:,:), bx(:,:,:)
        complex(dp), allocatable :: ex(:), ey(:), ez(:), ysmp(:)
        complex(dp) :: fval
        complex     :: cv(1)
        real(dp)    :: twopi_n, cc, scale, err, na, nb
        real    :: w(3,3,3), loc(3), loc2(3), ctr, sig, dx, dy, dz, rotp(3,3), ang(3)
        integer :: i, j, k, s, sgn, i0(3), wdim, iwinsz, nyq, nyqsq, ns, ip, h, kk, rho_lb(3), rho_ub(3)
        integer :: ih, ik, im, di, dj, dk, c, np, isw
        real    :: mfrac, nlev, lamv, lamr(3)
        real(dp) :: yrms
        logical :: pass_a, pass_b, pass_c, pass_d, l_sweep
        integer :: tt, ijk(3), nmiss, nsw
        l_sweep = .false.
        if( present(sweep) ) l_sweep = sweep
        call op%new(box, 1.0, 1)
        call op%set_band(op%Rnat)
        wdim   = 2*ceiling(KBWINSZ - 0.5) + 1
        iwinsz = ceiling(KBWINSZ - 0.5)
        if( wdim /= 3 ) THROW_HARD('test_flex_pcg_operator assumes the 3-tap KB window')
        nyq = box/2
        c   = box/2 + 1
        twopi_n = 2.0_dp * PI / real(box,dp)
        allocate(locs(3,nsamples), wts(nsamples))
        call fixed_seed(20260925)   ! reproducible sample positions, volume modulation and slices
        do s = 1, nsamples
            locs(:,s) = (2.0*[ran3(), ran3(), ran3()] - 1.0) * (0.5*real(nyq))
            wts(s)    = 1.0
        end do
        allocate(u(box,box,box,1), hu(box,box,box,1), eu4(box,box,box,1))
        ctr = real(box)/2.0 + 1.0
        sig = 0.12*real(box)
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    dx = real(i)-ctr; dy = real(j)-ctr; dz = real(k)-ctr
                    u(i,j,k,1) = exp(-(dx*dx+dy*dy+dz*dz)/(2.0*sig*sig)) * (1.0 + 0.3*(ran3()-0.5))
                end do
            end do
        end do
        allocate(ex(box), ey(box), ez(box))
        ! ================= (A) operator against the exact Gram =================
        ! band-limited apodized input Pi(E u), as the operator sees it
        eu4(:,:,:,1) = op%env * u(:,:,:,1)
        call op%bandlimit(eu4(:,:,:,1))
        allocate(eu(box,box,box), be(box,box,box), source=0.0_dp)
        eu = real(eu4(:,:,:,1),dp)
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc = real(sgn) * locs(:,s)
                call exps(loc)
                call sample_exact(eu, fval)
                fval = fval * real(wts(s),dp)
                call adjoint_add(fval, be)
            end do
        end do
        ! E Pi (exact Gram): band-limit the reference the way the operator band-limits its output
        eu4(:,:,:,1) = real(be)
        call op%bandlimit(eu4(:,:,:,1))
        be = real(op%env,dp) * real(eu4(:,:,:,1),dp)
        ! the kernel operator on the same samples and mates
        call op%alloc_accum(kacc)
        allocate(kacc4(op%npairs, op%lims3(1,2)-op%lims3(1,1)+1, op%lims3(2,2)-op%lims3(2,1)+1, &
            &op%lims3(3,2)-op%lims3(3,1)+1), source=0.0)
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc2 = real(op%padf) * real(sgn) * locs(:,s)
                i0   = nint(loc2) - iwinsz
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call scatter_pairs_wrap(op, i0, w, [wts(s)], 1.0, kacc)
                    call scatter_pairs_wrap_dense(op, i0, w, [wts(s)], 1.0, kacc4)
                else
                    call scatter_pairs_nowrap(op, i0, w, [wts(s)], 1.0, kacc)
                    call scatter_pairs_nowrap_dense(op, i0, w, [wts(s)], 1.0, kacc4)
                endif
            end do
        end do
        call op%alloc_packed(kpk)
        call op%fold_accum(kacc, kpk)
        ! (D) band-list accumulation and fold against the dense lattice: bitwise
        allocate(kpk4(op%npairs, op%cdim(1), op%cdim(2), op%cdim(3)), source=0.0)
        call fold_accum_dense(op, kacc4, kpk4)
        nmiss = 0
        do tt = 1, op%npk
            ijk = op%pijk(:,tt)
            if( any(kpk4(:,ijk(1),ijk(2),ijk(3)) /= kpk(:,tt)) ) nmiss = nmiss + 1
            kpk4(:,ijk(1),ijk(2),ijk(3)) = 0.0
        end do
        pass_d = nmiss == 0 .and. all(kpk4 == 0.0)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,L1)') '>>> FLEX PCG TEST (D) band-list kernels vs dense fold: slots=', op%npk, &
            &'  mismatching slots=', nmiss, '  dense mass outside the list=', count(kpk4 /= 0.0), '  pass=', pass_d
        deallocate(kpk4)
        call op%finalize(kpk)
        eu4(:,:,:,1) = op%env * u(:,:,:,1)
        call op%apply_operator(eu4, hu)
        hu(:,:,:,1) = op%env * hu(:,:,:,1)
        call compare3(be, real(hu(:,:,:,1),dp), cc, scale, err, na, nb)
        write(logfhandle,'(A,I0,A,I0,A,F8.5,A,F10.5,A,ES10.3,A,ES10.3)') '>>> FLEX PCG TEST (A) operator box=', box, &
            &' samples=', nsamples, '  kernel vs exact Gram: corr=', real(cc), '  LS scale=', real(scale), &
            &'  rel_resid=', real(err), '  |exact|=', real(na)
        pass_a = abs(scale - 1.0_dp) < 0.05_dp .and. err < 0.1_dp
        ! ================= (B) right-hand side: 2x deposit of exact samples vs exact adjoint =================
        allocate(ysmp(nsamples), source=(0.0_dp,0.0_dp))
        allocate(bx(box,box,box), source=0.0_dp)
        eu = real(u(:,:,:,1),dp)
        do s = 1, nsamples
            loc = locs(:,s)
            call exps(loc)
            call sample_exact(eu, ysmp(s))
        end do
        call op%alloc_rhs_accum(racc)
        allocate(racc4(op%ncomp, op%lims3(1,2)-op%lims3(1,1)+1, op%lims3(2,2)-op%lims3(2,1)+1, &
            &op%lims3(3,2)-op%lims3(3,1)+1), source=cmplx(0.,0.))
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc  = real(sgn) * locs(:,s)
                loc2 = real(op%padf) * loc
                i0   = nint(loc2) - iwinsz
                if( sgn == 1 )then
                    fval = ysmp(s)
                else
                    fval = conjg(ysmp(s))
                endif
                cv(1) = cmplx(real(real(fval)*real(wts(s),dp), sp), real(aimag(fval)*real(wts(s),dp), sp))
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call scatter_rhs_wrap(op, i0, w, cv, racc)
                    call scatter_rhs_wrap_dense(op, i0, w, cv, racc4)
                else
                    call scatter_rhs_nowrap(op, i0, w, cv, racc)
                    call scatter_rhs_nowrap_dense(op, i0, w, cv, racc4)
                endif
                ! exact adjoint of the same sample
                call exps(loc)
                call adjoint_add(fval*real(wts(s),dp), bx)
            end do
        end do
        call op%alloc_rhs_packed(rpk)
        call op%fold_rhs(racc, rpk)
        allocate(rpk4(op%ncomp, op%cdim(1), op%cdim(2), op%cdim(3)), source=cmplx(0.,0.))
        call fold_rhs_dense(op, racc4, rpk4)
        nmiss = 0
        do tt = 1, op%npk
            ijk = op%pijk(:,tt)
            if( any(rpk4(:,ijk(1),ijk(2),ijk(3)) /= rpk(:,tt)) ) nmiss = nmiss + 1
            rpk4(:,ijk(1),ijk(2),ijk(3)) = cmplx(0.,0.)
        end do
        pass_d = pass_d .and. nmiss == 0 .and. all(rpk4 == cmplx(0.,0.))
        write(logfhandle,'(A,I0,A,I0,A,L1)') '>>> FLEX PCG TEST (D) band-list rhs vs dense fold: mismatching slots=', nmiss, &
            &'  dense mass outside the list=', count(rpk4 /= cmplx(0.,0.)), '  pass=', pass_d
        deallocate(rpk4)
        allocate(b4(box,box,box,1))
        call op%finalize_rhs(rpk, b4)
        eu4(:,:,:,1) = real(bx)
        call op%bandlimit(eu4(:,:,:,1))
        bx = real(eu4(:,:,:,1),dp)
        call compare3(bx, real(b4(:,:,:,1),dp), cc, scale, err, na, nb)
        write(logfhandle,'(A,F8.5,A,F10.5,A,ES10.3,A,ES10.3)') '>>> FLEX PCG TEST (B) rhs deposit vs exact adjoint: corr=', &
            &real(cc), '  LS scale=', real(scale), '  rel_resid=', real(err), '  |exact|=', real(na)
        pass_b = abs(scale - 1.0_dp) < 0.05_dp .and. err < 0.1_dp
        deallocate(kpk, rpk, b4, be, bx, ysmp)
        ! ================= (C) preconditioned CG solve from central-slice samples =================
        ! the clean, generously supported baseline (mask 0.40, no noise, no Tikhonov term); the sweep adds
        ! the tighter support, sample noise and the Tikhonov term. Every clean solve must meet the baseline
        ! criterion; the noisy ones characterise what the Tikhonov term is for (unregularised, CG runs
        ! into the null space and the error outside the low-pass grows without bound)
        pass_c = .true.
        if( box <= 32 )then
            np    = 48
            nyqsq = nyq*nyq
            ns = 0
            do ip = 1, np
                do kk = -nyq, nyq
                    do h = -nyq, nyq
                        if( h*h + kk*kk <= nyqsq ) ns = ns + 1
                    end do
                end do
            end do
            deallocate(locs, wts)
            allocate(locs(3,ns), wts(ns))
            wts = 1.0
            s = 0
            do ip = 1, np
                ang(1) = 2.0*PI*ran3(); ang(2) = acos(2.0*ran3()-1.0); ang(3) = 2.0*PI*ran3()
                rotp = euler_rot(ang)
                do kk = -nyq, nyq
                    do h = -nyq, nyq
                        if( h*h + kk*kk > nyqsq ) cycle
                        s = s + 1
                        locs(:,s) = matmul(real([h,kk,0]), rotp)
                    end do
                end do
            end do
            ! exact samples of the volume
            allocate(ysmp(ns), source=(0.0_dp,0.0_dp))
            eu = real(u(:,:,:,1),dp)
            do s = 1, ns
                call exps(locs(:,s))
                call sample_exact(eu, ysmp(s))
            end do
            yrms = sqrt(sum(abs(ysmp)**2) / real(ns,dp))
            ! kernels on the 2x lattice and the 1x gridding density (independent of the data)
            rho_lb = [-(iwinsz+1), -nyq-iwinsz-1, -nyq-iwinsz-1]
            rho_ub = [ nyq+iwinsz+1, nyq+iwinsz+1, nyq+iwinsz+1]
            allocate(rho_t(1, rho_ub(1)-rho_lb(1)+1, rho_ub(2)-rho_lb(2)+1, rho_ub(3)-rho_lb(3)+1), source=0.0)
            call op%alloc_accum(kacc)
            do s = 1, ns
                loc  = locs(:,s)
                loc2 = real(op%padf) * loc
                i0   = nint(loc2) - iwinsz
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call scatter_pairs_wrap(op, i0, w, [1.0], 1.0, kacc)
                else
                    call scatter_pairs_nowrap(op, i0, w, [1.0], 1.0, kacc)
                endif
                i0 = nint(loc) - iwinsz
                call op%kbwin%apod_mat_3d_fast(loc, iwinsz, wdim, w)
                do dk = 1, wdim
                    im = i0(3) + dk - 1 - rho_lb(3) + 1
                    do dj = 1, wdim
                        ik = i0(2) + dj - 1 - rho_lb(2) + 1
                        do di = 1, wdim
                            ih = i0(1) + di - 1 - rho_lb(1) + 1
                            rho_t(1,ih,ik,im) = rho_t(1,ih,ik,im) + w(di,dj,dk)
                        end do
                    end do
                end do
            end do
            call op%alloc_packed(kpk)
            call op%fold_accum(kacc, kpk)
            call op%finalize(kpk)
            allocate(b4(box,box,box,1), x4(box,box,box,1), source=0.0)
            allocate(ut(box,box,box), source=0.0)
            nsw = merge(12, 1, l_sweep)
            do isw = 1, nsw
                mfrac = merge(0.4, 0.25, mod((isw-1)/6, 2) == 0)
                nlev  = merge(0.0, 0.5,  mod((isw-1)/3, 2) == 0)
                lamr  = [0.0, 1.0e-3, 1.0e-2]
                lamv  = lamr(mod(isw-1, 3) + 1)
                call op%set_window_sphere(mfrac*real(box))
                call op%set_lambda_relative(lamv)
                ! right-hand sides of the (noisy) samples
                call op%alloc_rhs_accum(racc)
                do s = 1, ns
                    loc  = locs(:,s)
                    loc2 = real(op%padf) * loc
                    i0   = nint(loc2) - iwinsz
                    call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                    fval = ysmp(s)
                    if( nlev > 0.0 ) fval = fval + real(nlev,dp)*yrms*cmplx(gasdev(), gasdev(), kind=dp)/sqrt(2.0_dp)
                    cv(1) = cmplx(real(real(fval), sp), real(aimag(fval), sp))
                    if( op%win_wraps(i0) )then
                        call scatter_rhs_wrap(op, i0, w, cv, racc)
                    else
                        call scatter_rhs_nowrap(op, i0, w, cv, racc)
                    endif
                end do
                call op%alloc_rhs_packed(rpk)
                call op%fold_rhs(racc, rpk)
                call op%finalize_rhs(rpk, b4)
                x4 = 0.0
                call op%cg_core(b4, x4, rho_t, rho_lb, 60, 1.0e-4, out, 'FLEX PCG TEST (C)')
                ut = u(:,:,:,1) * op%mask
                err = sqrt(sum((real(x4(:,:,:,1),dp) - real(ut,dp))**2)) / sqrt(sum(real(ut,dp)**2))
                call lowpass(x4(:,:,:,1), nyq/2)
                call lowpass(ut, nyq/2)
                cc = sqrt(sum((real(x4(:,:,:,1),dp) - real(ut,dp))**2)) / sqrt(sum(real(ut,dp)**2))
                write(logfhandle,'(A,F5.2,A,F4.2,A,ES8.1,A,I3,A,I3,A,ES9.2,A,ES9.2,A,ES9.2)') &
                    &'>>> FLEX PCG TEST (C) sweep: mask=', mfrac, ' noise=', nlev, ' lam=', lamv, &
                    &'  iters=', out%iteration_count, '  its_to_1e-2=', out%iters_to_1e2, &
                    &'  final resid=', out%final_rel_residual, '  err full=', real(err), '  err inner=', real(cc)
                if( nlev == 0.0 ) pass_c = pass_c .and. cc < 0.05_dp .and. out%final_rel_residual < 1.0e-2
                deallocate(rpk)
            end do
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PCG TEST (C) slices=', np, ' samples=', ns
            deallocate(kpk, b4, x4, rho_t, ut, ysmp)
        endif
        l_pass = pass_a .and. pass_b .and. pass_c .and. pass_d
        if( present(passes) ) passes = [pass_a, pass_b, pass_c, pass_d]
        if( l_pass )then
            write(logfhandle,'(A)') '    PASS: operator, right-hand side, solve and band lists within tolerance'
        else
            write(logfhandle,'(A,4L2)') '    FAIL (operator, rhs, solve, band lists): ', pass_a, pass_b, pass_c, pass_d
        endif
        call op%kill
        deallocate(locs, wts, u, hu, eu4, eu, ex, ey, ez)

    contains

        subroutine fixed_seed( base_seed )
            integer, intent(in) :: base_seed
            integer, allocatable :: seed(:)
            integer :: ii, nn
            call random_seed(size=nn)
            allocate(seed(nn))
            do ii = 1, nn
                seed(ii) = modulo(base_seed + 104729 * (ii - 1), huge(0) - 1) + 1
            enddo
            call random_seed(put=seed)
        end subroutine fixed_seed

        !> separable exponentials of one sample position
        subroutine exps( p )
            real, intent(in) :: p(3)
            integer :: ii
            do ii = 1, box
                ex(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(1),dp)*real(ii-c,dp), dp))
                ey(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(2),dp)*real(ii-c,dp), dp))
                ez(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(3),dp)*real(ii-c,dp), dp))
            end do
        end subroutine exps

        !> F = (1/N^3) sum_n v(n) e^{-i...} for the current exponentials, one axis at a time: the x sum
        !! for every (y,z) column is one real (2,N) x (N,N^2) product (library matmul, also at -O0), then
        !! y and z; (C) takes 38256 central-slice samples at box 32, 1.25e9 voxel terms. v is the box^3
        !! volume by sequence association
        subroutine sample_exact( v, f )
            real(dp),    intent(in)  :: v(box,box*box)
            complex(dp), intent(out) :: f
            real(dp) :: exri(2,box), t(2,box*box)
            integer  :: kk2, j0
            exri(1,:) = real(ex, dp)
            exri(2,:) = aimag(ex)
            t = matmul(exri, v)
            f = (0.0_dp, 0.0_dp)
            do kk2 = 1, box
                j0 = (kk2 - 1)*box
                f = f + ez(kk2) * sum(ey * cmplx(t(1,j0+1:j0+box), t(2,j0+1:j0+box), dp))
            end do
            f = f / real(box,dp)**3
        end subroutine sample_exact

        !> acc(n) += Re[ f e^{+i...} ] for the current exponentials
        subroutine adjoint_add( f, acc )
            complex(dp), intent(in)    :: f
            real(dp),    intent(inout) :: acc(:,:,:)
            integer :: ii, jj, kk2
            !$omp parallel do default(shared) private(ii,jj,kk2) schedule(static)
            do kk2 = 1, box
                do jj = 1, box
                    do ii = 1, box
                        acc(ii,jj,kk2) = acc(ii,jj,kk2) + real(f * conjg(ex(ii)*ey(jj)*ez(kk2)), dp)
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine adjoint_add

        !> spherical low-pass of a native-lattice volume to shell rmax
        subroutine lowpass( v, rmax )
            real,    intent(inout) :: v(:,:,:)
            integer, intent(in)    :: rmax
            real, pointer :: rp(:,:,:)
            integer :: lims(3,2), hh, kh, mh, ph(3)
            call op%nimg%set_rmat(v, .false.)
            call op%nimg%fft()
            lims = op%nimg%loop_lims(2)
            do mh = lims(3,1), lims(3,2)
                do kh = lims(2,1), lims(2,2)
                    do hh = lims(1,1), lims(1,2)
                        if( nint(sqrt(real(hh*hh + kh*kh + mh*mh))) > rmax )then
                            ph = op%nimg%comp_addr_phys(hh,kh,mh)
                            call op%nimg%set_cmat_at(ph(1),ph(2),ph(3), cmplx(0.,0.))
                        endif
                    end do
                end do
            end do
            call op%nimg%ifft()
            call op%nimg%get_rmat_ptr(rp)
            v = rp(1:box,1:box,1:box)
        end subroutine lowpass

        !> ZYZ Euler rotation matrix
        function euler_rot( a ) result( R )
            real, intent(in) :: a(3)
            real :: R(3,3), ca, sa, cb, sb, cg, sg
            ca = cos(a(1)); sa = sin(a(1)); cb = cos(a(2)); sb = sin(a(2)); cg = cos(a(3)); sg = sin(a(3))
            R(1,1) =  ca*cb*cg - sa*sg; R(1,2) =  sa*cb*cg + ca*sg; R(1,3) = -sb*cg
            R(2,1) = -ca*cb*sg - sa*cg; R(2,2) = -sa*cb*sg + ca*cg; R(2,3) =  sb*sg
            R(3,1) =  ca*sb;            R(3,2) =  sa*sb;            R(3,3) =  cb
        end function euler_rot

        !> correlation, least-squares scale of b onto a, relative residual, norms
        subroutine compare3( a, bb, cc_, scale_, err_, na_, nb_ )
            real(dp), intent(in)  :: a(:,:,:), bb(:,:,:)
            real(dp), intent(out) :: cc_, scale_, err_, na_, nb_
            real(dp) :: num, den
            na_ = sqrt(sum(a**2))
            nb_ = sqrt(sum(bb**2))
            num = sum(a*bb)
            den = sum(bb**2)
            cc_ = num / max(na_*nb_, 1.0e-30_dp)
            scale_ = 1.0_dp
            if( den > 0.0_dp ) scale_ = num/den
            err_ = sqrt(sum((a - scale_*bb)**2)) / max(na_, 1.0e-30_dp)
        end subroutine compare3

    end subroutine test_flex_pcg_operator

end module simple_flex_pca_pcg
