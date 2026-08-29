!@descr: CTF/sigma-weighted Fourier-projection operator and
!  preconditioned-conjugate-gradient volume solver used by the reconstruct3D
!  PCG backend. See doc/policies/reconstruct3D_pcg_policy.md.
!
!  Per-particle data is cached once by prep_particles (not re-derived every CG
!  iteration), the particle loop is OpenMP-parallel throughout, and the solve
!  is preconditioned with the sampling-density diagonal. An optional kernelized
!  (Toeplitz/Gram) normal operator (pcgop=kernel) makes per-iteration cost
!  independent of particle count; the matrix-free operator remains the exact
!  reference.
module simple_reconstructor_pcg
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image, only: image
use simple_ctf,   only: ctf
use simple_gridding, only: kb_stencil_envelope_1d
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
implicit none

public :: reconstructor_pcg, pcg_solver_outcome, pcg_fourier_workspace
public :: PCG_OP_MATRIXFREE, PCG_OP_KERNEL
public :: SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT
public :: SHIFT_LM_NO_RELIABLE_UPDATE, SHIFT_LM_STEP_BOUND_REJECTED
public :: SHIFT_LM_INVALID_NUMERICS, SHIFT_LM_ITERATION_LIMIT
public :: POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT
public :: POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED
public :: POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
public :: right_increment_rotation
public :: pcg_raw_accum_compatible
private
#include "simple_local_flags.inc"

integer, parameter :: PCG_OP_MATRIXFREE = 0 !< reference operator: exact, cost ~ nptcls per iteration
integer, parameter :: PCG_OP_KERNEL     = 1 !< kernelized Toeplitz operator: cost independent of nptcls
integer, parameter :: PCG_RAW_ACCUM_VERSION = 1
integer, parameter :: PCG_RAW_PROV_LEN = 256
character(len=16), parameter :: PCG_RAW_ACCUM_MAGIC = 'SIMPLE_PCG_RAW01'

type :: pcg_solver_outcome
    character(len=24) :: stop_reason         = 'not_started'
    integer           :: iteration_count     = 0
    integer           :: requested_maxits    = 0
    real              :: initial_rel_residual = 0.0
    real              :: final_rel_residual   = 0.0
    real              :: final_rel_update     = 0.0
    logical           :: converged            = .false.
    real, allocatable  :: rel_residual_history(:)
    real, allocatable  :: rel_update_history(:)
    real, allocatable  :: preconditioned_residual_history(:)
    real, allocatable  :: iteration_seconds(:)
end type pcg_solver_outcome
integer, parameter :: SHIFT_LM_ACCEPTED_IMPROVEMENT = 1
integer, parameter :: SHIFT_LM_FINITE_NO_IMPROVEMENT = 2
integer, parameter :: SHIFT_LM_NO_RELIABLE_UPDATE = 3
integer, parameter :: SHIFT_LM_STEP_BOUND_REJECTED = 4
integer, parameter :: SHIFT_LM_INVALID_NUMERICS = 5
integer, parameter :: SHIFT_LM_ITERATION_LIMIT = 6
integer, parameter :: POSE_LM_ACCEPTED_IMPROVEMENT = SHIFT_LM_ACCEPTED_IMPROVEMENT
integer, parameter :: POSE_LM_FINITE_NO_IMPROVEMENT = SHIFT_LM_FINITE_NO_IMPROVEMENT
integer, parameter :: POSE_LM_NO_RELIABLE_UPDATE = SHIFT_LM_NO_RELIABLE_UPDATE
integer, parameter :: POSE_LM_STEP_BOUND_REJECTED = SHIFT_LM_STEP_BOUND_REJECTED
integer, parameter :: POSE_LM_INVALID_NUMERICS = SHIFT_LM_INVALID_NUMERICS
integer, parameter :: POSE_LM_ITERATION_LIMIT = SHIFT_LM_ITERATION_LIMIT
real(dp), parameter :: POSE_NUMERIC_FLOOR = epsilon(1._dp)**2

!> Read-only snapshot of the padded Fourier volume for repeated local samples.
!! One snapshot is shared across particles and derivative directions while the
!! real-space volume is fixed. Rebuild it after every set_volume call.
type :: pcg_fourier_workspace
    private
    integer          :: box     = 0
    integer          :: boxpd   = 0
    integer          :: padf    = 1
    integer          :: iwinsz  = 0
    integer          :: wdim    = 0
    integer          :: lims2(2,2) = 0
    integer          :: sqhp    = 0
    integer          :: sqlp    = 0
    real             :: padsc   = 1.0
    type(kbinterpol) :: kbwin
    integer, allocatable :: wrap(:)
    complex, allocatable :: cmat(:,:,:)
    logical :: exists = .false.
  contains
    procedure :: kill => kill_fourier_workspace
    procedure :: get_lims2 => get_fourier_workspace_lims2
    procedure :: set_shell_range => set_fourier_workspace_shell_range
    procedure :: sample_with_grad => sample_fourier_with_grad
    procedure :: shift_residual
    procedure :: shift_jvp
    procedure :: shift_jhz
    procedure :: shift_objective_gradient
    procedure :: refine_shift_lm
    procedure :: rotation_jvp
    procedure :: pose_objective_gradient
    procedure :: refine_pose_lm
    procedure, private :: shift_normal_terms
    procedure, private :: pose_normal_terms
    procedure :: count_stencil_switches
end type pcg_fourier_workspace

type :: reconstructor_pcg
    private
    integer          :: box        = 0   !< native box: the solver unknown lives here
    integer          :: boxpd      = 0   !< padf*box: every Fourier operation happens here
    integer          :: padf       = 1   !< oversampling factor, OSMPL_PAD_FAC
    integer          :: pad_off    = 0   !< centred pad/crop offset, (boxpd-box)/2
    integer          :: Rnat       = 0   !< native Nyquist radius, box/2
    real             :: padsc      = 1.0 !< padf**3, undoes fwd_ft's 1/product(ldim)
    real             :: smpd       = 1.0
    type(kbinterpol) :: kbwin
    integer          :: iwinsz     = 0
    integer          :: wdim       = 0
    integer          :: stride     = 0  !< OpenMP colouring stride for the scatter, see apply_normal
    integer          :: lims2(2,2) = 0  !< (h/k, lo/hi) plane bounds, full symmetric disk
    integer          :: lims3(3,2) = 0  !< (h/k/m, lo/hi) volume array bounds
    integer          :: wlims(2)   = 0  !< [lo,hi] canonical period-box wrap range
    integer          :: sqlp       = 0  !< squared Nyquist radius
    integer          :: sq_rim     = 0  !< below this h^2+k^2 a KB window cannot wrap, see new
    real             :: lambda     = 0.0 !< effective absolute coefficient used by apply_normal
    real             :: lambda_rel = 0.0 !< coefficient relative to the weighted data-operator scale
    real             :: data_scale = 0.0 !< deterministic scale derived from raw data-only D
    logical          :: l_lambda_relative = .false.
    ! ---- per-particle inputs, cached once by prep_particles ----
    integer                       :: nptcls = 0
    integer                       :: nsym   = 1   !< point-group order; 1 = c1 (no replication)
    logical                       :: l_use_ctf = .false.
    real,            allocatable  :: rotmats(:,:,:)  !< (3,3,nptcls)
    real,            allocatable  :: symmats(:,:,:)  !< (3,3,nsym) point-group operators, symmats(:,:,1)=I
    type(ctfparams), allocatable  :: ctfparms(:)     !< (nptcls)
    real,            allocatable  :: shifts(:,:)     !< (2,nptcls), pixels
    real,            allocatable  :: sig2(:,:)       !< (0:R,nptcls) per-particle noise power
    ! ---- lookup tables / work buffers ----
    integer, allocatable :: wrap(:)                  !< precomputed cyci_1d over the whole reachable range
    ! (h,k)-only quantities absT2_plane and build_transfer would otherwise
    ! recompute per particle; built once over the fixed lims2 disk, shared by
    ! all particles. ~258 kB each at box 256. See build_hk_luts.
    real,    allocatable :: spafreqsq_lut(:,:)       !< spatial frequency squared
    real,    allocatable :: ang_lut(:,:)             !< atan2(k,h) astigmatism angle
    integer, allocatable :: shell_lut(:,:)           !< resolution shell index, capped at Rnat
    real,    allocatable :: env(:,:,:)               !< measured KB instrument envelope, see build_env
    real,    allocatable :: invenv(:,:,:)            !< its guarded reciprocal, for deapodization
    logical              :: l_deapod = .true.        !< correct the KB roll-off, see deapod_mul
    real,    allocatable :: mask(:,:,:)              !< soft spherical support constraint, see set_mask
    logical              :: l_mask = .false.
    type(image)          :: wimg                     !< persistent box^3 work image (keeps its FFTW plans)
    logical              :: wimg_exists = .false.
    ! ---- streaming accumulation over particle batches ----
    ! These two Fourier accumulators are the ONLY particle-dependent state the
    ! solve needs; precond, Khat and the RHS are all derived from them. Keeping
    ! them here rather than handing them to callers is what lets a commander
    ! stream batches without ever seeing the padded lattice.
    real,    allocatable :: acc_work(:,:,:)          !< sum_i G_i^dagger |T_i|^2, full-range
    complex, allocatable :: b_work(:,:,:)            !< RHS accumulator, full-range
    real,    allocatable :: b_rhs(:,:,:)             !< folded/deapodized/masked RHS, box^3
    logical              :: l_accum = .false.
    logical              :: l_rhs   = .false.
    integer              :: reduction_next_part = 1  !< fixed raw-artifact association order
    integer              :: reduction_nparts    = 0
    integer              :: reduction_state     = 0
    integer              :: reduction_eo        = -1
    ! ---- sampling-density preconditioner ----
    real,    allocatable :: precond(:,:,:)           !< 1/(rho+lambda) on wimg's cmat layout
    logical              :: l_precond = .false.
    ! ---- kernelized operator ----
    integer              :: op_mode = PCG_OP_MATRIXFREE
    real,    allocatable :: Khat(:,:,:)              !< Gram kernel on wimg's cmat layout, real
    logical              :: l_kernel = .false.
    ! ---- optional FSC/SSNR quadratic prior ----
    ! Requested before raw D is finalized; its absolute scale is then derived
    ! from D, keeping the prior in the same calibrated operator convention.
    real, allocatable :: ml_fsc(:)                    !< independent-half FSC, shells 1:R
    real, allocatable :: ml_prior(:,:,:)              !< calibrated padded Fourier diagonal
    real              :: ml_tau = 1.0                 !< established ML regularization fudge factor
    real              :: ml_hp  = 100.0               !< low-frequency no-prior limit in Angstrom
    logical           :: l_ml_prior_requested = .false.
    logical           :: l_ml_prior = .false.
    ! ---- optional direct NU-evidence replay precision (pcg_priors.md S5) ----
    ! Q_NU = C (sum_b B_b^T W_b B_b) C, with C the native-box mean-centering
    ! projector (C = I - 11^T/N, symmetric idempotent -- it makes a constant
    ! field an EXACT null mode: padded-DC exclusion alone is not enough,
    ! because a native constant becomes a box window after zero-padding and
    ! carries non-DC padded content), B_b the disjoint radial Fourier band
    ! masks on the padded lattice, and W_b = [p(1-a_b)]^2 the graded spatial
    ! lack-of-evidence weight for band b from the frozen compact NU evidence
    ! state. Each B_b = crop o proj o pad is a contraction and the M_b are
    ! disjoint, so sum_b ||B_b x||^2 <= ||x||^2 and with W in [0,1] the
    ! operator norm is bounded by 1 -- the declared normalization. The bank is
    ! NOT a tight frame: the pad/crop sandwich makes sum_b B_b^T B_b < I with
    ! cross-band leakage, so refining/merging the band partition changes the
    ! operator (measured, not assumed invariant; see the Gate A partition
    ! test). Applied matrix-free in the deapodized domain; mode-exclusive with
    ! the FSC/SSNR P_tau (R10), asserted in both attachment orders (the
    ! binary-envelope solvent precision Q_s was removed 2026-08-27;
    ! experiment record in pcg_priors.md S4). The absolute scale is
    ! lambda_nu = rel * data_scale,
    ! derived alongside the ridge lambda. The preconditioner uses the declared
    ! nonnegative approximation: the support-mean band weight fused as a shell
    ! diagonal, mirroring the ML-prior fusion.
    real,    allocatable :: nu_band_w(:,:,:,:)        !< graded spatial weights, box^3 x nbands
    real,    allocatable :: nu_band_limits(:)         !< coarse-to-fine band low-pass limits (A)
    real,    allocatable :: nu_band_wmean(:)          !< support-mean weight per band (precond approx)
    integer(kind=1), allocatable :: nu_band_idx(:,:,:) !< padded-lattice band index, 0 = DC/unvisited
    real              :: nu_lambda_rel = 0.0          !< strength relative to the data scale
    real              :: lambda_nu = 0.0              !< effective absolute strength
    logical           :: l_nu_prior = .false.
    ! ---- per-phase profiling, accumulated over a solve. Exists to answer one
    !      question before any further optimization: of the seconds an iteration
    !      costs, how many are the particle loop (which the kernelized operator
    !      removes) and how many are FFT + bulk cmat traffic on the padf-times
    !      padded lattice (which it does NOT)? Guessing that split wrong means
    !      optimizing the wrong half. ----
    logical  :: l_profile = .false.
    real(dp) :: t_setvol  = 0.0_dp  !< pad + forward FFT of the iterate
    real(dp) :: t_cmatcp  = 0.0_dp  !< get_cmat/set_cmat bulk copies
    real(dp) :: t_ploop   = 0.0_dp  !< the particle loop proper
    real(dp) :: t_fold    = 0.0_dp  !< fold + inverse FFT + crop
    real(dp) :: t_khat    = 0.0_dp  !< kernel pointwise multiply
    real(dp) :: t_prec    = 0.0_dp  !< apply_precond, whole call
    real(dp) :: t_fin_rhs = 0.0_dp  !< RHS fold, deapodization and support
    real(dp) :: t_fin_rho = 0.0_dp  !< deterministic rho shell statistics
    real(dp) :: t_fin_fold = 0.0_dp !< fused reciprocal and packed-Khat pass
    real(dp) :: t_fin_dep = 0.0_dp  !< deposition-envelope construction
    real(dp) :: t_fin_kernel = 0.0_dp !< kernel correction, FFT and calibration
    logical              :: exists = .false.
  contains
    ! CONSTRUCTOR / DESTRUCTOR
    procedure :: new
    procedure :: kill
    ! SETUP
    procedure :: prep_particles
    procedure :: set_sym
    procedure :: build_precond
    procedure :: build_kernel
    procedure :: build_operators
    ! STREAMING SETUP (batch-at-a-time; see begin_accum)
    procedure :: begin_accum
    procedure :: begin_reduction
    procedure :: accumulate_batch
    procedure :: end_accum
    procedure :: write_raw_accum
    procedure :: add_raw_accum
    procedure :: add_raw_accum_weighted
    procedure :: scale_raw_accum
    procedure :: compare_raw_accum
    procedure, private :: accumulate_rhs_density
    procedure, private :: accumulate_absT2
    procedure, private :: accumulate_rhs
    procedure, private :: precond_from_accum
    procedure, private :: fold_accum_to_khat
    procedure, private :: finalize_density_accum
    procedure, private :: finalize_khat
    procedure :: set_deapod
    procedure :: set_mask
    procedure :: set_lambda_relative
    procedure :: set_ml_prior
    procedure, private :: build_env
    procedure, private :: build_kb_envelope_1d
    procedure, private :: build_hk_luts
    procedure, private :: deapod_mul
    procedure, private :: mask_mul
    procedure, private :: calibrate_kernel
    procedure :: measure_kernel_scale
    procedure :: set_op_mode
    ! LOW-LEVEL OPERATOR (public: the test commanders drive these directly to
    ! verify the adjoint identity before any solve() result may be trusted)
    procedure :: set_volume
    procedure :: begin_fourier_workspace
    procedure :: forward_plane
    procedure :: fourier_dot
    procedure :: adjoint_plane_add
    procedure :: build_transfer
    procedure :: whiten_observation
    procedure :: extract_native_plane
    procedure :: dot_real_volume
    ! HIGH-LEVEL OPERATOR (uses the cached per-particle state)
    procedure :: apply_normal
    procedure :: apply_normal_matrixfree
    procedure :: apply_normal_kernel
    procedure :: apply_adjoint_all
    procedure :: assert_prior_attachment_mode
    procedure :: set_nu_prior
    procedure :: apply_nu_precision
    procedure :: get_nu_prior_stats
    ! GETTERS
    procedure :: get_lims2
    procedure :: get_lims3
    procedure :: get_nptcls
    procedure :: get_env
    procedure :: get_invenv
    procedure :: get_rhs
    procedure :: get_raw_accum
    procedure :: get_ml_prior
    procedure :: get_ml_prior_stats
    procedure :: get_data_scale
    procedure :: get_effective_lambda
    procedure :: get_effective_nu_lambda
    ! SOLVER
    procedure :: solve
    procedure :: solve_accum
    procedure, private :: solve_core
    ! PROFILING
    procedure :: reset_profile
    procedure :: report_profile
    procedure, private :: reset_finalize_profile
    procedure :: report_finalize_profile
    ! PRIVATE HELPERS
    procedure, private :: absT2_plane
    procedure, private :: prepare_fused_planes
    procedure, private :: transfer_plane_cmplx
    procedure, private :: fold_and_ifft
    procedure, private :: apply_precond
    procedure, private :: ensure_wimg
    procedure, private :: pad_vol
    procedure, private :: crop_vol
    procedure, private :: update_lambda_from_density
    procedure, private :: build_ml_prior_from_density
    procedure, private :: apply_fourier_diagonal
    procedure, private :: ensure_nu_band_index
    procedure, private :: nu_precond_shell_diag
end type reconstructor_pcg

contains

    !>  \brief  Applies a right tangent-space increment R <- R exp([omega]x).
    pure function right_increment_rotation( rotmat, omega ) result(updated_rotmat)
        real(dp), intent(in) :: rotmat(3,3), omega(3)
        real(dp) :: updated_rotmat(3,3), skew(3,3), exp_skew(3,3)
        real(dp) :: identity(3,3), theta2, theta4, sinc_theta, cosc_theta

        identity = 0._dp
        identity(1,1) = 1._dp
        identity(2,2) = 1._dp
        identity(3,3) = 1._dp
        ! [omega]x u = omega x u.
        skew = reshape([0._dp,omega(3),-omega(2), &
            &-omega(3),0._dp,omega(1),omega(2),-omega(1),0._dp],[3,3])
        theta2 = dot_product(omega,omega)
        if( theta2 < 1.e-8_dp )then
            theta4 = theta2*theta2
            sinc_theta = 1._dp-theta2/6._dp+theta4/120._dp
            cosc_theta = 0.5_dp-theta2/24._dp+theta4/720._dp
        else
            sinc_theta = sin(sqrt(theta2))/sqrt(theta2)
            cosc_theta = (1._dp-cos(sqrt(theta2)))/theta2
        endif
        exp_skew = identity+sinc_theta*skew+cosc_theta*matmul(skew,skew)
        updated_rotmat = matmul(rotmat,exp_skew)
    end function right_increment_rotation

    ! CONSTRUCTOR

    subroutine new( self, box, smpd, lambda )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                   intent(in)    :: box
        real,                      intent(in)    :: smpd
        real, optional,            intent(in)    :: lambda
        type(image) :: tmp
        integer     :: R, lo, hi, i
        real        :: rlim
        call self%kill
        self%box    = box
        self%smpd   = smpd
        self%lambda = 0.0
        if( present(lambda) ) self%lambda = lambda
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        ! Deliberately computed as 2*iwinsz+1 (odd, symmetric around nint(loc)),
        ! not kbwin%get_wdim() directly -- for KBWINSZ=1.5 these agree (both 3),
        ! but spelling it out documents the invariant this operator depends on:
        ! an even-width window would not be centred on nint(loc), so negating
        ! loc would not negate the window as a set, breaking the mirror
        ! consistency fold_and_ifft's h>=0-only fold requires.
        self%wdim   = 2*self%iwinsz + 1
        ! After rotation, source h-lines separated by only wdim can still update
        ! the same 3-D interpolation cell. sqrt(3)*wdim guarantees separation by
        ! a full window width along at least one target axis, because in 3-D the
        ! max-norm is at least the Euclidean norm over sqrt(3). Same scheme as
        ! reconstructor%insert_plane_oversamp.
        self%stride = ceiling(sqrt(3.0) * real(self%wdim))
        ! OVERSAMPLING. The unknown stays on the native box, but every Fourier
        ! operation runs on a padf-times padded lattice: the volume is centre-
        ! padded going in and centre-cropped coming out. This is the same trick
        ! Fourier gridding uses (projector%interp_fcomp_oversamp works at
        ! loc*PAD_FAC and rescales by PAD_FAC**3). Without it, KB interpolation
        ! on the native lattice is only percent-accurate, the roll-off envelope
        ! swings ~30x across the box, and the Gram kernel cannot match the
        ! operator -- three symptoms of one cause.
        self%padf    = OSMPL_PAD_FAC
        self%boxpd   = self%padf * box
        self%pad_off = (self%boxpd - box) / 2
        self%Rnat    = box / 2
        ! fwd_ft divides by product(ldim), so the padded spectrum is 1/padf**3
        ! of the native one at coincident frequencies; this restores the scale.
        self%padsc   = real(self%padf)**3
        call tmp%new([self%boxpd,self%boxpd,self%boxpd], smpd)
        call tmp%fft()
        self%lims3 = tmp%loop_lims(3)
        ! wlims is the TRUE period-box canonical wrap range (period = box).
        ! lims3(1,:) is NOT usable for this on axis 1: it spans the redundant
        ! [-box/2, box/2] (both Friedel-mate Nyquist bins), one longer than the
        ! true period, which would give cyci_1d an off-by-one wrap. Axes 2/3
        ! never get that extension, so lims3(2,:) is the correct range; the box
        ! is cubic so it applies to all three axes.
        self%wlims = self%lims3(2,:)
        ! lims2 is a FULL symmetric (both-sign h) disk, not the packed half.
        ! This makes forward_plane/adjoint_plane_add an exact adjoint pair for
        ! ANY orientation, at the cost of 2x the plane work -- which no longer
        ! matters per iteration once the kernelized operator is in use.
        ! lims2 stays NATIVE: the plane is sampled at native Fourier indices,
        ! which map onto the padded lattice as padf*loc.
        R = self%Rnat
        self%lims2(1,:) = [-R, R]
        self%lims2(2,:) = [-R, R]
        self%sqlp       = R*R
        ! Squared plane radius below which a KB window provably CANNOT reach the
        ! periodic wrap boundary, so the serial rim pass can reject it without
        ! computing anything. The window spans nint(loc) +/- iwinsz and wraps
        ! only if some component leaves [wlims(1), wlims(2)]; since
        ! |nint(loc_j)| <= |loc| + 0.5 and |loc| = padf*sqrt(h^2+k^2) is
        ! rotation-INDEPENDENT, the test collapses to one on (h,k) alone.
        ! Conservative by construction, with an extra -1 of margin: it only ever
        ! rejects points that cannot wrap, never ones that can.
        rlim = real(min(self%wlims(2) - self%iwinsz, -self%wlims(1) - self%iwinsz)) - 0.5
        self%sq_rim = max(0, int((rlim / real(self%padf))**2) - 1)
        call tmp%kill
        ! precomputed cyci_1d over every index the KB window can reach
        lo = self%wlims(1) - self%iwinsz - 1
        hi = self%wlims(2) + self%iwinsz + 1
        allocate(self%wrap(lo:hi))
        do i = lo, hi
            self%wrap(i) = cyci_1d(self%wlims, i)
        end do
        call self%build_hk_luts
        call self%build_env
        ! Default point group is c1: a single identity operator. set_sym
        ! replaces this when the caller requests replication.
        self%nsym = 1
        if( allocated(self%symmats) ) deallocate(self%symmats)
        allocate(self%symmats(3,3,1), source=0.0)
        self%symmats(1,1,1) = 1.0; self%symmats(2,2,1) = 1.0; self%symmats(3,3,1) = 1.0
        self%exists = .true.
    end subroutine new

    !>  \brief  precomputes the (h,k)-only quantities absT2_plane and
    !!          build_transfer would otherwise recompute per particle: spatial
    !!          frequency squared, the astigmatism angle (atan2) and the
    !!          resolution shell (sqrt). The lims2 disk is fixed for the life of
    !!          the object, so these are particle-independent. These three
    !!          expressions are bit-identical to the inlined ones they replace.
    !!          The CTF form that consumes them is NOT: absT2_plane and
    !!          build_transfer factor half_wl2_cs = 0.5*wl*wl*Cs out of the
    !!          per-pixel product, which reassociates. So the residual trace is
    !!          NOT a valid regression signal across that change -- the operator
    !!          stages of test=pcg_recon are (see note 5.4).
    subroutine build_hk_luts( self )
        class(reconstructor_pcg), intent(inout) :: self
        integer :: h, k, R
        R = self%Rnat
        allocate(self%spafreqsq_lut(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(self%ang_lut(      self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(self%shell_lut(    self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                self%spafreqsq_lut(h,k) = (real(h)/real(self%box))**2 + (real(k)/real(self%box))**2
                self%ang_lut(h,k)       = atan2(real(k), real(h))
                self%shell_lut(h,k)     = min(nint(sqrt(real(h*h+k*k))), R)
            end do
        end do
    end subroutine build_hk_luts

    !>  \brief  Separable real-space Kaiser-Bessel instrument envelope.
    !!
    !!          KB interpolation is not a free lunch: gathering
    !!          sum_k w(k-loc) xhat(k) equals the Fourier transform of
    !!          (w_hat . x) sampled at loc, i.e. it silently multiplies the
    !!          volume by the KB instrument function in real space. That is
    !!          exactly why production divides it out after gridding, see
    !!          prep3D_inv_instrfun4mul in simple_gridding.f90.
    !!
    !!          The matrix-free operator applies it TWICE (once gathering,
    !!          once scattering), so H = E T E with E = diag(w_hat) and T the
    !!          bare shift-invariant part. apply_normal_kernel must therefore
    !!          bracket its convolution with the same envelope, or it models a
    !!          different operator entirely -- the envelope falls from 1 at the
    !!          box centre to roughly 0.19 at the edge per axis, so squared it
    !!          varies by a factor of ~30 across the volume.
    !!
    !!          The envelope is the exact discrete transform of the normalized
    !!          three-tap weights returned by apod_mat_3d_fast at the origin.
    !!          It is not kbinterpol%instr, which is the continuous transform and
    !!          differs materially at the box edge. Separability makes the 3-D
    !!          impulse IFFT exactly the outer product of one 1-D cosine sum, so
    !!          no padded 3-D accumulator or FFT is needed.
    subroutine build_env( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, parameter :: EPS_DIV = 1.0e-8
        real, allocatable :: env1d(:)
        real    :: ctrval
        integer :: c, i, j, k, o
        if( allocated(self%env)    ) deallocate(self%env)
        if( allocated(self%invenv) ) deallocate(self%invenv)
        call self%build_kb_envelope_1d(self%boxpd, env1d)
        allocate(self%env(self%box,self%box,self%box))
        o = self%pad_off
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    self%env(i,j,k) = env1d(o+i) * env1d(o+j) * env1d(o+k)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(env1d)
        ! normalize so the envelope is unity at the box centre
        c      = self%box/2 + 1
        ctrval = self%env(c,c,c)
        if( abs(ctrval) < EPS_DIV ) ctrval = 1.0
        self%env = self%env / ctrval
        allocate(self%invenv(self%box,self%box,self%box), source=1.0)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    ! guarded reciprocal, same convention as prep3D_inv_instrfun4mul
                    if( abs(self%env(i,j,k)) < EPS_DIV )then
                        self%invenv(i,j,k) = 0.0
                    else
                        self%invenv(i,j,k) = 1.0 / self%env(i,j,k)
                    endif
                end do
            end do
        end do
    end subroutine build_env

    !> Exact 1-D inverse transform of the discrete, normalized KB stencil at
    !! the Fourier origin. The origin stencil is symmetric and real, so the
    !! transform is the real cosine sum below. Spatial index n/2+1 is the
    !! centered image origin used by image%get_rmat.
    subroutine build_kb_envelope_1d( self, n, env1d )
        class(reconstructor_pcg), intent(in) :: self
        integer,                    intent(in) :: n
        real, allocatable,          intent(out) :: env1d(:)
        ! shared with the gridding backend (simple_gridding): same normalized
        ! origin stencil, period n (the padded lattice here)
        call kb_stencil_envelope_1d(self%kbwin, n, env1d)
    end subroutine build_kb_envelope_1d

    pure function get_env( self ) result( env )
        class(reconstructor_pcg), intent(in) :: self
        real :: env(self%box,self%box,self%box)
        env = self%env
    end function get_env

    !>  \brief  Copy of the right-hand side currently held, whether it was built
    !!          by solve() from whole planes or by end_accum from batches.
    !!          Exists so a test can localize a disagreement between the two
    !!          accumulation paths to the RHS or to the operator, rather than
    !!          only observing that the solutions differ.
    subroutine get_rhs( self, b )
        class(reconstructor_pcg), intent(in)  :: self
        real, allocatable,         intent(out) :: b(:,:,:)
        if( .not. self%l_rhs ) THROW_HARD('no right-hand side has been built; get_rhs')
        allocate(b(self%box,self%box,self%box), source=self%b_rhs)
    end subroutine get_rhs

    !> Copy the open, unfinalized accumulator-domain statistics. This is a
    !! test/diagnostic boundary only: production distribution persists the
    !! same arrays through write_raw_accum and never exposes them to policy.
    subroutine get_raw_accum( self, b, d )
        class(reconstructor_pcg), intent(in) :: self
        complex, allocatable, intent(out) :: b(:,:,:)
        real,    allocatable, intent(out) :: d(:,:,:)
        if( .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open')
        allocate(b(self%lims3(1,1):self%lims3(1,2), &
                   &self%lims3(2,1):self%lims3(2,2), &
                   &self%lims3(3,1):self%lims3(3,2)), source=self%b_work)
        allocate(d(self%lims3(1,1):self%lims3(1,2), &
                   &self%lims3(2,1):self%lims3(2,2), &
                   &self%lims3(3,1):self%lims3(3,2)), source=self%acc_work)
    end subroutine get_raw_accum

    subroutine get_ml_prior( self, prior )
        class(reconstructor_pcg), intent(in) :: self
        real, allocatable, intent(out) :: prior(:,:,:)
        if( .not. self%l_ml_prior ) THROW_HARD('PCG ML prior has not been built')
        allocate(prior, source=self%ml_prior)
    end subroutine get_ml_prior

    !> Summarize the calibrated FSC/SSNR prior without copying its padded
    !! lattice. Ratios compare P_tau with the calibrated data-only Khat over
    !! exactly the bins where the prior is positive. Absolute Khat is used in
    !! the L1 denominator because the finite-support Toeplitz approximation can
    !! contain small negative bins even though the underlying operator is PSD.
    subroutine get_ml_prior_stats( self, npositive, positive_min, positive_max, &
            &prior_to_khat_l1, prior_to_khat_rms )
        class(reconstructor_pcg), intent(in) :: self
        integer, intent(out) :: npositive
        real,    intent(out) :: positive_min, positive_max, prior_to_khat_l1, prior_to_khat_rms
        real(dp) :: prior_l1, khat_l1, prior_sq, khat_sq
        real     :: pval, kval
        integer  :: i, j, k
        if( .not. self%l_ml_prior ) THROW_HARD('PCG ML prior has not been built')
        if( .not. self%l_kernel ) THROW_HARD('PCG ML prior statistics require finalized Khat')
        npositive   = 0
        positive_min = huge(1.0)
        positive_max = 0.0
        prior_l1 = 0.0_dp
        khat_l1  = 0.0_dp
        prior_sq = 0.0_dp
        khat_sq  = 0.0_dp
        !$omp parallel do collapse(3) default(shared) private(i,j,k,pval,kval) &
        !$omp reduction(+:npositive,prior_l1,khat_l1,prior_sq,khat_sq) &
        !$omp reduction(min:positive_min) reduction(max:positive_max) schedule(static)
        do k = 1, size(self%ml_prior,3)
            do j = 1, size(self%ml_prior,2)
                do i = 1, size(self%ml_prior,1)
                    pval = self%ml_prior(i,j,k)
                    if( pval <= 0.0 ) cycle
                    kval = self%Khat(i,j,k)
                    npositive    = npositive + 1
                    positive_min = min(positive_min, pval)
                    positive_max = max(positive_max, pval)
                    prior_l1 = prior_l1 + real(pval,dp)
                    khat_l1  = khat_l1  + abs(real(kval,dp))
                    prior_sq = prior_sq + real(pval,dp)**2
                    khat_sq  = khat_sq  + real(kval,dp)**2
                enddo
            enddo
        enddo
        !$omp end parallel do
        if( npositive < 1 ) THROW_HARD('PCG ML prior statistics found no positive bins')
        prior_to_khat_l1  = real(prior_l1 / max(khat_l1, DTINY))
        prior_to_khat_rms = real(sqrt(prior_sq / max(khat_sq, DTINY)))
        if( .not. ieee_is_finite(prior_to_khat_l1) .or. .not. ieee_is_finite(prior_to_khat_rms) )then
            THROW_HARD('PCG ML prior-to-kernel statistics are not finite')
        endif
    end subroutine get_ml_prior_stats

    pure function get_invenv( self ) result( invenv )
        class(reconstructor_pcg), intent(in) :: self
        real :: invenv(self%box,self%box,self%box)
        invenv = self%invenv
    end function get_invenv

    subroutine set_deapod( self, l_deapod )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                    intent(in)    :: l_deapod
        self%l_deapod = l_deapod
    end subroutine set_deapod

    !> Interpret lambda as a dimensionless coefficient relative to the weighted
    !! data normal operator. The effective absolute coefficient is derived only
    !! after raw D has been reduced or trailing-blended; workers therefore keep
    !! publishing raw (B,D), with no regularization embedded in the artifact.
    subroutine set_lambda_relative( self, lambda_rel )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: lambda_rel
        if( .not. ieee_is_finite(lambda_rel) .or. lambda_rel < 0.0 )then
            THROW_HARD('relative PCG lambda must be finite and non-negative')
        endif
        self%lambda_rel = lambda_rel
        self%data_scale = 0.0
        self%lambda     = 0.0
        self%l_lambda_relative = .true.
    end subroutine set_lambda_relative

    !> Request the established isotropic FSC/SSNR prior. The FSC determines
    !! only the shell-wise relative strength here; its absolute scale is built
    !! later from the finalized raw data density D on the master.
    subroutine set_ml_prior( self, fsc, tau, hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: fsc(:), tau, hp
        ! replay precisions are mode-exclusive (pcg_priors.md R10)
        if( self%l_nu_prior ) THROW_HARD('P_tau and Q_NU are mutually exclusive; NU replay precision is attached')
        if( size(fsc) < 1 ) THROW_HARD('PCG ML prior requires a non-empty FSC')
        if( .not. ieee_is_finite(tau) .or. tau <= 0.0 ) THROW_HARD('PCG ML tau must be finite and positive')
        if( .not. ieee_is_finite(hp) .or. hp <= 0.0 ) THROW_HARD('PCG ML high-pass limit must be finite and positive')
        if( allocated(self%ml_fsc) ) deallocate(self%ml_fsc)
        if( allocated(self%ml_prior) ) deallocate(self%ml_prior)
        allocate(self%ml_fsc(size(fsc)), source=fsc)
        self%ml_tau = tau
        self%ml_hp  = hp
        self%l_ml_prior_requested = .true.
        self%l_ml_prior = .false.
    end subroutine set_ml_prior

    !>  \brief  Hard contract on the operating mode priors attach in. All
    !!          regularizers beyond the plain lambda ridge -- the ML shell
    !!          diagonal today, the planned real-space priors (see
    !!          doc/implementation_notes/pcg_priors.md S3) -- are derived and
    !!          validated in exactly one operator configuration: the kernelized
    !!          normal operator with deapodization ON, where the iterate lives
    !!          in the deapodized (physical-density) domain and the Fourier
    !!          diagonal can be fused with Khat. The other configurations exist
    !!          only as test oracles; attaching a prior there would change its
    !!          meaning silently (the iterate would carry the KB envelope).
    !!          Call this at the solve site whenever a prior is attached.
    subroutine assert_prior_attachment_mode( self )
        class(reconstructor_pcg), intent(in) :: self
        if( self%op_mode /= PCG_OP_KERNEL )then
            THROW_HARD('PCG priors attach in kernel operator mode only')
        endif
        if( .not. self%l_deapod )then
            THROW_HARD('PCG priors attach in deapodized mode only')
        endif
    end subroutine assert_prior_attachment_mode

    !>  \brief  Install the direct NU-evidence replay precision
    !!          Q_NU = C (sum_b B_b^T W_b B_b) C (pcg_priors.md S5.3), with C
    !!          the native mean-centering projector supplying the exact
    !!          constant null mode. band_w holds
    !!          the per-voxel LACK-of-evidence weight for each band (1 - a_b,
    !!          in [0,1], monotone non-decreasing coarse to fine because band
    !!          support is nested coarse-to-fine); band_limits are the
    !!          coarse-to-fine band low-pass boundaries in Angstrom. The stored
    !!          weight is [p*w]^2, graded at the support boundary exactly like
    !!          the retired solvent weight was. Mode exclusion (R10) is
    !!          enforced here and in set_ml_prior. The effective strength
    !!          starts at the relative value (unit data scale, for direct
    !!          operator tests) and is rescaled by data_scale when raw D is
    !!          finalized, exactly like the relative ridge lambda.
    subroutine set_nu_prior( self, band_w, band_limits, lambda_rel )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: band_w(:,:,:,:)
        real,                       intent(in)    :: band_limits(:)
        real,                       intent(in)    :: lambda_rel
        integer :: b, nb, n_supp
        if( self%l_ml_prior_requested ) &
            &THROW_HARD('P_tau and Q_NU are mutually exclusive; FSC/SSNR ML prior is attached')
        if( .not. self%l_mask ) THROW_HARD('NU replay precision requires the support mask; call set_mask first')
        if( .not. ieee_is_finite(lambda_rel) .or. lambda_rel < 0.0 )then
            THROW_HARD('relative NU prior strength must be finite and non-negative')
        endif
        nb = size(band_w,4)
        if( nb < 1 ) THROW_HARD('NU replay precision requires at least one detail band')
        if( size(band_limits) /= nb ) THROW_HARD('NU band limits do not match the band weight count')
        if( any(shape(band_w(:,:,:,1)) /= self%box) ) THROW_HARD('NU band weights do not match the solve box')
        if( .not. all(ieee_is_finite(band_w)) ) THROW_HARD('NU band weights contain non-finite values')
        if( minval(band_w) < 0.0 .or. maxval(band_w) > 1.0 ) &
            &THROW_HARD('NU band weights outside [0,1]; clip before installing')
        do b = 2, nb
            if( .not.(band_limits(b) > 0.0) .or. band_limits(b) >= band_limits(b-1) ) &
                &THROW_HARD('NU band limits must be strictly decreasing and positive (coarse to fine)')
        enddo
        if( .not.(band_limits(1) > 0.0) ) THROW_HARD('NU band limits must be positive')
        if( allocated(self%nu_band_w)      ) deallocate(self%nu_band_w)
        if( allocated(self%nu_band_limits) ) deallocate(self%nu_band_limits)
        if( allocated(self%nu_band_wmean)  ) deallocate(self%nu_band_wmean)
        if( allocated(self%nu_band_idx)    ) deallocate(self%nu_band_idx)
        allocate(self%nu_band_w(self%box,self%box,self%box,nb))
        allocate(self%nu_band_limits(nb), source=band_limits)
        allocate(self%nu_band_wmean(nb),  source=0.0)
        n_supp = count(self%mask > 0.0)
        if( n_supp < 1 ) THROW_HARD('NU replay precision has an empty support')
        do b = 1, nb
            self%nu_band_w(:,:,:,b) = (self%mask * band_w(:,:,:,b))**2
            self%nu_band_wmean(b)   = real(sum(real(self%nu_band_w(:,:,:,b),dp), &
                &mask=self%mask > 0.0) / real(n_supp,dp))
        enddo
        self%nu_lambda_rel = lambda_rel
        self%lambda_nu     = lambda_rel
        self%l_nu_prior    = .true.
    end subroutine set_nu_prior

    !>  \brief  Padded-lattice band index for the disjoint radial Fourier
    !!          masks B_b. Built lazily on the wimg cmat layout, exactly the
    !!          lattice traversal build_ml_prior_from_density uses. Padded DC
    !!          and any physically unaddressed points stay 0 and are excluded
    !!          from every band. NOTE: padded-DC exclusion alone does NOT null
    !!          a native constant (zero-padding turns it into a box window
    !!          with non-DC padded content); the exact constant null mode is
    !!          supplied by the explicit native mean-centering in
    !!          apply_nu_precision. Band 1 covers detail coarser than or equal
    !!          to limits(1); band b covers [limits(b), limits(b-1)); the
    !!          finest band also absorbs everything finer than limits(nb),
    !!          where no candidate evidence exists.
    subroutine ensure_nu_band_index( self )
        class(reconstructor_pcg), intent(inout) :: self
        integer :: cdim(3), h, k, m, phys(3), shpd, b, nb, band
        real    :: res
        if( allocated(self%nu_band_idx) ) return
        if( .not. allocated(self%nu_band_limits) ) THROW_HARD('set_nu_prior has not been called; ensure_nu_band_index')
        call self%ensure_wimg
        cdim = self%wimg%get_array_shape()
        allocate(self%nu_band_idx(cdim(1),cdim(2),cdim(3)), source=0_1)
        nb = size(self%nu_band_limits)
        !$omp parallel do collapse(2) default(shared) private(h,k,m,phys,shpd,b,band,res) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    shpd = nint(sqrt(real(h*h + k*k + m*m)))
                    if( shpd < 1 ) cycle   ! DC excluded: constant fields are unpenalized
                    res  = real(self%box) * self%smpd * real(self%padf) / real(shpd)
                    band = nb
                    do b = 1, nb
                        if( res >= self%nu_band_limits(b) )then
                            band = b
                            exit
                        endif
                    enddo
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    self%nu_band_idx(phys(1),phys(2),phys(3)) = int(band,kind=1)
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine ensure_nu_band_index

    !>  \brief  Q_NU x = C sum_b B_b^T (W_b .* (B_b (C x))) at unit strength;
    !!          callers scale by their lambda. C = I - 11^T/N is the native
    !!          mean-centering projector: symmetric idempotent, applied on
    !!          both sides so Q_NU stays symmetric PSD and a constant field is
    !!          an EXACT null mode (padded-DC exclusion alone cannot provide
    !!          this: a native constant zero-pads to a box window with non-DC
    !!          padded content). Each B_b is crop o IFFT o M_b o FFT o pad
    !!          with M_b a real 0/1 radial mask, hence symmetric (restriction
    !!          is the adjoint of zero-extension and the diagonal is real), so
    !!          every summand is symmetric PSD. The per-band results are
    !!          accumulated in real space so only two padded complex
    !!          workspaces are live at a time.
    function apply_nu_precision( self, x ) result( qx )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: x(self%box,self%box,self%box)
        real,    allocatable :: qx(:,:,:), xb(:,:,:), xc(:,:,:)
        complex, allocatable :: cmat0(:,:,:), cband(:,:,:)
        real(dp) :: mean_dp
        integer :: cdim(3), b, nb, i, j, k
        if( .not. self%l_nu_prior ) THROW_HARD('set_nu_prior has not been called; apply_nu_precision')
        call self%ensure_wimg
        call self%ensure_nu_band_index
        cdim = self%wimg%get_array_shape()
        ! left application of C: remove the native mean before the band bank
        mean_dp = sum(real(x,dp)) / real(self%box,dp)**3
        allocate(xc(self%box,self%box,self%box), source=x)
        xc = xc - real(mean_dp)
        call self%wimg%set_rmat(self%pad_vol(xc), .false.)
        call self%wimg%fft()
        cmat0 = self%wimg%get_cmat()
        allocate(qx(self%box,self%box,self%box), source=0.0)
        allocate(cband(cdim(1),cdim(2),cdim(3)))
        nb = size(self%nu_band_w,4)
        do b = 1, nb
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
            do k = 1, cdim(3)
                do j = 1, cdim(2)
                    do i = 1, cdim(1)
                        if( self%nu_band_idx(i,j,k) == int(b,kind=1) )then
                            cband(i,j,k) = cmat0(i,j,k)
                        else
                            cband(i,j,k) = cmplx(0.0,0.0)
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
            call self%wimg%set_cmat(cband)
            call self%wimg%ifft()
            xb = self%crop_vol(self%wimg%get_rmat())
            xb = self%nu_band_w(:,:,:,b) * xb
            call self%wimg%set_rmat(self%pad_vol(xb), .false.)
            call self%wimg%fft()
            cband = self%wimg%get_cmat()
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
            do k = 1, cdim(3)
                do j = 1, cdim(2)
                    do i = 1, cdim(1)
                        if( self%nu_band_idx(i,j,k) /= int(b,kind=1) ) cband(i,j,k) = cmplx(0.0,0.0)
                    enddo
                enddo
            enddo
            !$omp end parallel do
            call self%wimg%set_cmat(cband)
            call self%wimg%ifft()
            qx = qx + self%crop_vol(self%wimg%get_rmat())
            deallocate(xb)
        enddo
        ! right application of C: re-center the output (C is symmetric, so
        ! this completes C Q~ C and restores the exact adjoint identity)
        mean_dp = sum(real(qx,dp)) / real(self%box,dp)**3
        qx = qx - real(mean_dp)
        deallocate(cmat0, cband, xc)
    end function apply_nu_precision

    !>  \brief  NU replay diagnostics of a final map: the penalty
    !!          R_NU = (lambda_nu/2) x^T Q_NU x and its per-band energies.
    subroutine get_nu_prior_stats( self, x, penalty )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)  :: x(self%box,self%box,self%box)
        real,                       intent(out) :: penalty
        real, allocatable :: qx(:,:,:)
        if( .not. self%l_nu_prior ) THROW_HARD('set_nu_prior has not been called; get_nu_prior_stats')
        qx = self%apply_nu_precision(x)
        penalty = real(0.5_dp * real(self%lambda_nu,dp) * sum(real(x,dp) * real(qx,dp)))
        deallocate(qx)
    end subroutine get_nu_prior_stats

    !>  \brief  Declared nonnegative preconditioner approximation of Q_NU: the
    !!          support-mean band weight as a shell diagonal in raw rho units
    !!          (the /padsc**2 mirrors the ML-prior fusion). Exact spatial
    !!          structure cannot be represented in the Fourier-shell
    !!          preconditioner; the shell mean preserves PSD and captures the
    !!          bandwise average stiffening.
    pure real function nu_precond_shell_diag( self, sh_padded ) result( d )
        class(reconstructor_pcg), intent(in) :: self
        integer,                    intent(in) :: sh_padded
        integer :: b, nb, band
        real    :: res
        d = 0.0
        if( .not. self%l_nu_prior ) return
        if( self%lambda_nu <= 0.0 ) return
        if( sh_padded < 1 ) return
        res  = real(self%box) * self%smpd * real(self%padf) / real(sh_padded)
        nb   = size(self%nu_band_limits)
        band = nb
        do b = 1, nb
            if( res >= self%nu_band_limits(b) )then
                band = b
                exit
            endif
        enddo
        d = self%lambda_nu * self%nu_band_wmean(band) / self%padsc**2
    end function nu_precond_shell_diag

    pure real function get_effective_nu_lambda( self )
        class(reconstructor_pcg), intent(in) :: self
        get_effective_nu_lambda = self%lambda_nu
    end function get_effective_nu_lambda

    !>  \brief  Multiplies a real volume by E^-1, the inverse KB instrument
    !!          envelope -- the deapodization / roll-off correction.
    !!
    !!          WHY THIS MATTERS FOR REAL DATA. The gather makes the operator
    !!          A~ = A.E rather than the true A = D.S (see build_env). For
    !!          SYNTHETIC tests that is invisible, because the observations are
    !!          generated with the same operator, so the envelope cancels -- an
    !!          inverse crime. Real particle images carry no such envelope, so
    !!          fitting them with A~ returns E^-1 x_true instead of x_true, i.e.
    !!          a reconstruction inflated toward the box edges.
    !!
    !!          Applying E^-1 on both sides of the normal operator (and once to
    !!          the right-hand side) turns A~ back into A, so the solve targets
    !!          x_true directly and the prior applies to x_true rather than to a
    !!          rescaled surrogate -- which a post-hoc output correction could
    !!          not achieve.
    !!
    !!          CAVEAT, stated plainly: this correction is large here (~5x per
    !!          axis at the box edge, so ~140x at a corner) ONLY because this
    !!          operator interpolates on the native grid with no oversampling.
    !!          Production grids into a 2x padded box and therefore evaluates
    !!          the envelope over half the range, where the correction is mild
    !!          (~1.5x per axis). Adding oversampling is the proper fix and
    !!          would make this step benign.
    pure subroutine deapod_mul( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                       intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_deapod ) return
        v = v * self%invenv
    end subroutine deapod_mul

    !>  \brief  Constrains the solution to a soft spherical support.
    !!
    !!          This is a CONSTRAINT ON THE SOLVE, not a cosmetic post-hoc mask.
    !!          Writing x = P u and minimizing ||A P u - y||^2 gives the normal
    !!          equations (P H P) u = P b, so P is applied on BOTH sides of the
    !!          operator and once to the right-hand side; apply_normal and solve
    !!          do exactly that. P H P is still symmetric positive semi-definite,
    !!          and its null space -- everything outside the support -- is simply
    !!          never entered when starting from x = 0, the same argument that
    !!          makes the singular preconditioner safe (see build_precond).
    !!
    !!          Two things are gained. The obvious one: solvent outside the
    !!          particle is where E^-1 deapodization amplifies hardest (~3.4x at
    !!          a box corner), which is what puts artefacts at the box sides when
    !!          the display threshold is raised. The less obvious and more
    !!          valuable one: at mskdiam 180 in a 256 box the sphere is only ~18%
    !!          of the volume, so five sixths of the unknowns being solved for
    !!          were solvent the data says almost nothing about. Removing them
    !!          shrinks the problem and improves conditioning.
    !!
    !!          The edge profile is production's own cosedge_r2_3d, obtained by
    !!          running a unit volume through image%mask3D_soft, so the roll-off
    !!          matches what the rest of SIMPLE uses. backgr=0. is passed to stop
    !!          it subtracting a background from what is meant to be a pure
    !!          window function.
    subroutine set_mask( self, mskrad )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: mskrad
        type(image)       :: mimg
        real, allocatable :: ones(:,:,:)
        if( allocated(self%mask) ) deallocate(self%mask)
        self%l_mask = .false.
        if( mskrad <= 0.0 ) return
        allocate(ones(self%box,self%box,self%box), source=1.0)
        call mimg%new([self%box,self%box,self%box], self%smpd)
        call mimg%set_rmat(ones, .false.)
        call mimg%mask3D_soft(mskrad, backgr=0.)
        self%mask = mimg%get_rmat()
        call mimg%kill
        deallocate(ones)
        self%l_mask = .true.
    end subroutine set_mask

    pure subroutine mask_mul( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                     intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_mask ) return
        v = v * self%mask
    end subroutine mask_mul

    subroutine ensure_wimg( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%wimg_exists ) return
        call self%wimg%new([self%boxpd,self%boxpd,self%boxpd], self%smpd)
        self%wimg_exists = .true.
    end subroutine ensure_wimg

    !>  \brief  centre-embed a native box^3 volume into the padded lattice.
    !!          pad_vol and crop_vol are exact adjoints of one another, which is
    !!          what keeps the oversampled operator self-adjoint.
    function pad_vol( self, v ) result( vp )
        class(reconstructor_pcg), intent(in) :: self
        real,                     intent(in) :: v(self%box,self%box,self%box)
        real, allocatable :: vp(:,:,:)
        integer :: o
        o = self%pad_off
        allocate(vp(self%boxpd,self%boxpd,self%boxpd), source=0.0)
        vp(o+1:o+self%box, o+1:o+self%box, o+1:o+self%box) = v
    end function pad_vol

    !>  \brief  centre-crop the padded lattice back to the native box.
    function crop_vol( self, vp ) result( v )
        class(reconstructor_pcg), intent(in) :: self
        real,                     intent(in) :: vp(self%boxpd,self%boxpd,self%boxpd)
        real, allocatable :: v(:,:,:)
        integer :: o
        o = self%pad_off
        allocate(v(self%box,self%box,self%box), &
            &source=vp(o+1:o+self%box, o+1:o+self%box, o+1:o+self%box))
    end function crop_vol

    ! DESTRUCTOR

    subroutine kill( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%wimg_exists ) call self%wimg%kill
        if( allocated(self%rotmats)  ) deallocate(self%rotmats)
        if( allocated(self%ctfparms) ) deallocate(self%ctfparms)
        if( allocated(self%shifts)   ) deallocate(self%shifts)
        if( allocated(self%sig2)     ) deallocate(self%sig2)
        if( allocated(self%symmats)  ) deallocate(self%symmats)
        if( allocated(self%wrap)     ) deallocate(self%wrap)
        if( allocated(self%spafreqsq_lut) ) deallocate(self%spafreqsq_lut)
        if( allocated(self%ang_lut)  ) deallocate(self%ang_lut)
        if( allocated(self%shell_lut)) deallocate(self%shell_lut)
        if( allocated(self%env)      ) deallocate(self%env)
        if( allocated(self%invenv)   ) deallocate(self%invenv)
        if( allocated(self%mask)     ) deallocate(self%mask)
        if( allocated(self%precond)  ) deallocate(self%precond)
        if( allocated(self%Khat)     ) deallocate(self%Khat)
        if( allocated(self%ml_fsc)   ) deallocate(self%ml_fsc)
        if( allocated(self%ml_prior) ) deallocate(self%ml_prior)
        if( allocated(self%nu_band_w)      ) deallocate(self%nu_band_w)
        if( allocated(self%nu_band_limits) ) deallocate(self%nu_band_limits)
        if( allocated(self%nu_band_wmean)  ) deallocate(self%nu_band_wmean)
        if( allocated(self%nu_band_idx)    ) deallocate(self%nu_band_idx)
        if( allocated(self%acc_work) ) deallocate(self%acc_work)
        if( allocated(self%b_work)   ) deallocate(self%b_work)
        if( allocated(self%b_rhs)    ) deallocate(self%b_rhs)
        self%l_accum = .false.
        self%l_rhs   = .false.
        self%box    = 0
        self%lims2  = 0
        self%lims3  = 0
        self%sqlp   = 0
        self%sq_rim = 0
        self%lambda = 0.0
        self%lambda_rel = 0.0
        self%data_scale = 0.0
        self%l_lambda_relative = .false.
        self%nptcls = 0
        self%nsym   = 1
        self%l_use_ctf   = .false.
        self%l_mask      = .false.
        self%l_precond   = .false.
        self%l_kernel    = .false.
        self%ml_tau      = 1.0
        self%ml_hp       = 100.0
        self%l_ml_prior_requested = .false.
        self%l_ml_prior  = .false.
        self%nu_lambda_rel      = 0.0
        self%lambda_nu          = 0.0
        self%l_nu_prior         = .false.
        self%wimg_exists = .false.
        self%op_mode     = PCG_OP_MATRIXFREE
        call self%reset_profile(.false.)
        call self%reset_finalize_profile
        self%exists      = .false.
    end subroutine kill

    ! SETUP

    !>  \brief  caches everything per-particle that does NOT depend on the CG
    !!          iterate: rotation matrices, CTF parameters, shifts and noise
    !!          spectra. Caching avoids re-deriving them inside every
    !!          apply_normal call -- an ori deep copy, string-keyed hash lookups
    !!          and two euler2m evaluations per particle per iteration. All of
    !!          it is small (order 100 kB for 5000 particles) and constant for
    !!          the whole solve. Caching it is also what lets the particle loop
    !!          be shared cleanly across OpenMP threads.
    subroutine prep_particles( self, orientations, use_ctf, sig2 )
        class(reconstructor_pcg), intent(inout) :: self
        class(oris),                intent(inout) :: orientations
        logical,          optional, intent(in)    :: use_ctf
        real,             optional, intent(in)    :: sig2(0:,:)
        type(ori) :: e
        integer   :: i, R
        R = self%lims2(1,2)
        self%nptcls = orientations%get_noris()
        self%l_use_ctf = .false.
        if( present(use_ctf) ) self%l_use_ctf = use_ctf
        if( allocated(self%rotmats)  ) deallocate(self%rotmats)
        if( allocated(self%ctfparms) ) deallocate(self%ctfparms)
        if( allocated(self%shifts)   ) deallocate(self%shifts)
        if( allocated(self%sig2)     ) deallocate(self%sig2)
        allocate(self%rotmats(3,3,self%nptcls), self%shifts(2,self%nptcls), &
            &self%ctfparms(self%nptcls))
        allocate(self%sig2(0:R,self%nptcls), source=1.0)
        call e%new(.false.)
        do i = 1, self%nptcls
            call orientations%get_ori(i, e)
            self%rotmats(:,:,i) = e%get_mat()
            if( self%l_use_ctf )then
                self%ctfparms(i) = e%get_ctfvars()
                self%shifts(:,i) = e%get_2Dshift()
            else
                self%shifts(:,i) = 0.
            endif
        end do
        if( present(sig2) )then
            if( size(sig2,2) /= self%nptcls ) THROW_HARD('sig2 second dimension must be nptcls; prep_particles')
            self%sig2(0:min(R,ubound(sig2,1)),:) = sig2(0:min(R,ubound(sig2,1)),:)
        endif
        call self%ensure_wimg
    end subroutine prep_particles

    !>  \brief  caches the point-group operators for coordinate replication.
    !!
    !!          Symmetry is applied at scatter time by replicating each plane
    !!          pixel at all M orientations R_i . S_g (2.3 in the pcg note):
    !!          absT2_plane is still evaluated once per particle, only the KB
    !!          weight and the 27-tap scatter are replicated. symmats(:,:,1) is
    !!          the identity, so g=1 reproduces the c1 pass exactly. Call after
    !!          new (which installs the c1 default) and before begin_accum.
    subroutine set_sym( self, pgrpsyms )
        class(reconstructor_pcg), intent(inout) :: self
        class(sym),               intent(in)    :: pgrpsyms
        integer :: g
        if( allocated(self%symmats) ) deallocate(self%symmats)
        self%nsym = pgrpsyms%get_nsym()
        allocate(self%symmats(3,3,self%nsym), source=0.0)
        do g = 1, self%nsym
            call pgrpsyms%get_sym_rmat(g, self%symmats(:,:,g))
        end do
    end subroutine set_sym

    subroutine set_op_mode( self, op_mode )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                    intent(in)    :: op_mode
        if( op_mode == PCG_OP_KERNEL .and. .not. self%l_kernel )then
            THROW_HARD('kernelized operator requested but build_kernel has not been called; set_op_mode')
        endif
        self%op_mode = op_mode
    end subroutine set_op_mode

    ! LOW-LEVEL OPERATOR

    !>  \brief  G_i F: gathers a full (unpacked) Fourier plane from an already-
    !!          FFT'd volume at orientation e, via periodic-wrap KB
    !!          interpolation. Retained with the ori-based signature because the
    !!          test commanders drive it directly for the adjoint dot-product
    !!          gate; the solver path uses the cached rotation matrices instead.
    !>  \brief  Loads a native box^3 volume into the operator, centre-padded
    !!          and transformed onto the oversampled lattice. Callers hand over
    !!          a plain real volume and never see the padding -- which is the
    !!          point: the padded lattice is an implementation detail, and
    !!          letting it leak into callers is exactly how convention bugs get
    !!          in.
    subroutine set_volume( self, v )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: v(self%box,self%box,self%box)
        call self%ensure_wimg
        call self%wimg%set_rmat(self%pad_vol(v), .false.)
        call self%wimg%fft()
    end subroutine set_volume

    !>  \brief  Snapshots the current padded Fourier volume once for repeated
    !!          value-and-gradient samples while the real-space volume is fixed.
    subroutine begin_fourier_workspace( self, workspace )
        class(reconstructor_pcg), intent(inout) :: self
        type(pcg_fourier_workspace), intent(inout) :: workspace
        call workspace%kill
        call self%ensure_wimg
        workspace%box    = self%box
        workspace%boxpd  = self%boxpd
        workspace%padf   = self%padf
        workspace%iwinsz = self%iwinsz
        workspace%wdim   = self%wdim
        workspace%lims2  = self%lims2
        workspace%sqhp   = 0
        workspace%sqlp   = self%sqlp
        workspace%padsc  = self%padsc
        workspace%kbwin  = self%kbwin
        ! Preserve the negative physical indices of the periodic wrap table.
        allocate(workspace%wrap(lbound(self%wrap,1):ubound(self%wrap,1)), source=self%wrap)
        ! Copy the padded Fourier lattice once per fixed-volume workspace.
        workspace%cmat   = self%wimg%get_cmat()
        workspace%exists = .true.
    end subroutine begin_fourier_workspace

    !> Release a fixed-volume Fourier snapshot and reset all sampling metadata.
    subroutine kill_fourier_workspace( self )
        class(pcg_fourier_workspace), intent(inout) :: self
        if( allocated(self%wrap) ) deallocate(self%wrap)
        if( allocated(self%cmat) ) deallocate(self%cmat)
        self%box    = 0
        self%boxpd  = 0
        self%padf   = 1
        self%iwinsz = 0
        self%wdim   = 0
        self%lims2  = 0
        self%sqhp   = 0
        self%sqlp   = 0
        self%padsc  = 1.0
        self%exists = .false.
    end subroutine kill_fourier_workspace

    !> Return the native packed-plane bounds used by observations and transfers.
    pure function get_fourier_workspace_lims2( self ) result(lims2)
        class(pcg_fourier_workspace), intent(in) :: self
        integer :: lims2(2,2)
        lims2 = self%lims2
    end function get_fourier_workspace_lims2

    !> Restrict the local objective to an inclusive native Fourier-shell range.
    subroutine set_fourier_workspace_shell_range( self, kfromto )
        class(pcg_fourier_workspace), intent(inout) :: self
        integer, intent(in) :: kfromto(2)
        integer :: khi, klo
        if( .not. self%exists ) error stop 'set_shell_range called on an empty Fourier workspace'
        klo = max(0,kfromto(1))
        khi = min(self%box/2,kfromto(2))
        if( khi < klo ) error stop 'set_shell_range requires an ordered nonempty range'
        self%sqhp = klo*klo
        self%sqlp = khi*khi
    end subroutine set_fourier_workspace_shell_range

    !>  \brief  Samples the packed Fourier snapshot and its three fixed-cell
    !!          spatial derivatives at one oversampled-lattice coordinate.
    pure subroutine sample_fourier_with_grad( self, loc, value, dvalue_dloc, switch_margin )
        class(pcg_fourier_workspace), intent(in)  :: self
        real(sp),                     intent(in)  :: loc(3)
        complex,                      intent(out) :: value, dvalue_dloc(3)
        real(sp),                     intent(out) :: switch_margin(3)
        real(sp) :: w(self%wdim,self%wdim,self%wdim)
        real(sp) :: dw(self%wdim,self%wdim,self%wdim,3)
        integer  :: i0(3)
        if( .not. self%exists ) error stop 'sample_with_grad called on an empty Fourier workspace'
        ! Build w and dw/dloc on the same fixed interpolation stencil.
        call self%kbwin%apod_mat_3d_fast_grad(loc, self%iwinsz, self%wdim, i0, switch_margin, w, dw)
        if( any(i0 < lbound(self%wrap,1)) .or. &
            &any(i0 + self%wdim - 1 > ubound(self%wrap,1)) )then
            error stop 'sample_with_grad location lies outside the periodic wrap table'
        endif
        call gather_window_grad(self, i0, w, dw, value, dvalue_dloc)
        ! Apply native Fourier scaling to the value and all derivatives.
        value       = self%padsc * value
        dvalue_dloc = self%padsc * dvalue_dloc
    end subroutine sample_fourier_with_grad

    !>  \brief  Fixed-volume, CTF-free, unit-noise residual for a shifted
    !!          Fourier projection: r = S(t) G(R)V - y.
    subroutine shift_residual( self, rotmat, shift, observed, residual, objective )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in)  :: observed(self%lims2(1,1):self%lims2(1,2),&
                                          &self%lims2(2,1):self%lims2(2,2))
        complex, intent(out) :: residual(self%lims2(1,1):self%lims2(1,2),&
                                          &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective
        complex :: value, dvalue_dloc(3), phase
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg
        integer :: h, k
        if( .not. self%exists ) error stop 'shift_residual called on an empty Fourier workspace'
        residual = cmplx(0.,0.)
        objective = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                residual(h,k) = phase * value - observed(h,k)
                objective = objective + 0.5_dp * real(conjg(cmplx(residual(h,k),kind=dp)) * &
                    &cmplx(residual(h,k),kind=dp),dp)
            enddo
        enddo
    end subroutine shift_residual

    !>  \brief  Directional derivative of the CTF-free, unit-noise residual
    !!          with respect to the two real image-shift parameters.
    subroutine shift_jvp( self, rotmat, shift, direction, jv )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2), direction(2)
        complex, intent(out) :: jv(self%lims2(1,1):self%lims2(1,2),&
                                    &self%lims2(2,1):self%lims2(2,2))
        complex :: value, dvalue_dloc(3), phase
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, directional_frequency
        integer :: h, k
        if( .not. self%exists ) error stop 'shift_jvp called on an empty Fourier workspace'
        jv = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                directional_frequency = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*direction(1) + real(k,dp)*direction(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                jv(h,k) = cmplx(cmplx(0._dp,directional_frequency,kind=dp) * &
                    &cmplx(phase*value,kind=dp),kind=sp)
            enddo
        enddo
    end subroutine shift_jvp

    !>  \brief  Real-parameter adjoint of the two shift-Jacobian columns.
    subroutine shift_jhz( self, rotmat, shift, z, jhz )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: z(self%lims2(1,1):self%lims2(1,2),&
                                  &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: jhz(2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: jacobian_value
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, frequency
        integer :: axis, h, k
        if( .not. self%exists ) error stop 'shift_jhz called on an empty Fourier workspace'
        jhz = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                do axis = 1, 2
                    frequency = 2._dp * real(PI,dp) * real(merge(h,k,axis==1),dp) / real(self%box,dp)
                    jacobian_value = cmplx(0._dp,frequency,kind=dp) * &
                        &cmplx(phase*value,kind=dp)
                    jhz(axis) = jhz(axis) + real(conjg(jacobian_value) * &
                        &cmplx(z(h,k),kind=dp),dp)
                enddo
            enddo
        enddo
    end subroutine shift_jhz

    !>  \brief  Fused shift objective, gradient and Gauss-Newton block. This
    !!          avoids materializing derivative planes inside local refinement.
    subroutine shift_normal_terms( self, rotmat, shift, observed, objective, gradient, hessian, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(2), hessian(2,2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: model, residual, jacobian(2)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, frequency(2)
        integer :: axis, h, jaxis, k
        if( .not. self%exists ) error stop 'shift_normal_terms called on an empty Fourier workspace'
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k < self%sqhp .or. h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf,sp) * real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp * real(PI,dp) * &
                    &(real(h,dp)*shift(1) + real(k,dp)*shift(2)) / real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                model = cmplx(phase*value,kind=dp)
                if( present(transfer) ) model = model * cmplx(transfer(h,k),kind=dp)
                residual = model - cmplx(observed(h,k),kind=dp)
                frequency = 2._dp * real(PI,dp) * real([h,k],dp) / real(self%box,dp)
                jacobian = cmplx(0._dp,frequency,kind=dp) * model
                objective = objective + 0.5_dp*real(conjg(residual)*residual,dp)
                do axis = 1, 2
                    gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual,dp)
                    do jaxis = 1, 2
                        hessian(axis,jaxis) = hessian(axis,jaxis) + &
                            &real(conjg(jacobian(axis))*jacobian(jaxis),dp)
                    enddo
                enddo
            enddo
        enddo
    end subroutine shift_normal_terms

    !> Return the shift objective and gradient without exposing the computed
    !! two-by-two Gauss-Newton block.
    subroutine shift_objective_gradient( self, rotmat, shift, observed, objective, gradient, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2), &
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(2)
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2), &
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp) :: hessian(2,2)
        if( present(transfer) )then
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian,transfer)
        else
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian)
        endif
    end subroutine shift_objective_gradient

    !>  \brief  Directional derivative of the Fourier gather for a right
    !!          tangent-space rotation. The transfer plane includes CTF and
    !!          whitening when it is supplied by the production caller.
    subroutine rotation_jvp( self, rotmat, shift, direction, jv, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2), direction(3)
        complex, intent(out) :: jv(self%lims2(1,1):self%lims2(1,2),&
                                    &self%lims2(2,1):self%lims2(2,2))
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: derivative
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: args, dloc(3)
        integer :: h, k

        if( .not. self%exists ) error stop 'rotation_jvp called on an empty Fourier workspace'
        jv = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < self%sqhp .or. h*h+k*k > self%sqlp ) cycle
                loc = real(self%padf,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
                ! For row-vector gathers, d(loc)/d(epsilon) = loc x direction.
                dloc = [real(loc(2),dp)*direction(3)-real(loc(3),dp)*direction(2), &
                    &real(loc(3),dp)*direction(1)-real(loc(1),dp)*direction(3), &
                    &real(loc(1),dp)*direction(2)-real(loc(2),dp)*direction(1)]
                args = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(self%box,dp)
                phase = cmplx(cos(args),sin(args),kind=sp)
                derivative = cmplx(phase,kind=dp)*sum(cmplx(dvalue_dloc,kind=dp)*dloc)
                if( present(transfer) ) derivative = derivative*cmplx(transfer(h,k),kind=dp)
                jv(h,k) = cmplx(derivative,kind=sp)
            enddo
        enddo
    end subroutine rotation_jvp

    !>  \brief  Fused objective, five-vector gradient and Gauss-Newton block
    !!          for three right-rotation coordinates and two pixel shifts.
    subroutine pose_normal_terms( self, rotmat, shift, observed, objective, gradient, &
        &hessian, min_switch_margin, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(5), hessian(5,5), min_switch_margin
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: weighted_phase, model, residual, jacobian(5)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, dloc(3,3), frequency(2)
        integer :: axis, h, jaxis, k

        if( .not. self%exists ) error stop 'pose_normal_terms called on an empty Fourier workspace'
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        min_switch_margin = huge(0._dp)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < self%sqhp .or. h*h+k*k > self%sqlp ) cycle
                loc = real(self%padf,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
                call self%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
                min_switch_margin = min(min_switch_margin,real(minval(switch_margin),dp))
                ! Columns are loc x e1, loc x e2 and loc x e3.
                dloc(:,1) = [0._dp,real(loc(3),dp),-real(loc(2),dp)]
                dloc(:,2) = [-real(loc(3),dp),0._dp,real(loc(1),dp)]
                dloc(:,3) = [real(loc(2),dp),-real(loc(1),dp),0._dp]
                arg = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(self%box,dp)
                phase = cmplx(cos(arg),sin(arg),kind=sp)
                weighted_phase = cmplx(phase,kind=dp)
                if( present(transfer) ) weighted_phase = weighted_phase*cmplx(transfer(h,k),kind=dp)
                model = weighted_phase*cmplx(value,kind=dp)
                residual = model-cmplx(observed(h,k),kind=dp)
                do axis = 1, 3
                    jacobian(axis) = weighted_phase*sum(cmplx(dvalue_dloc,kind=dp)*dloc(:,axis))
                enddo
                frequency = 2._dp*real(PI,dp)*real([h,k],dp)/real(self%box,dp)
                jacobian(4:5) = cmplx(0._dp,frequency,kind=dp)*model
                objective = objective+0.5_dp*real(conjg(residual)*residual,dp)
                do axis = 1, 5
                    gradient(axis) = gradient(axis)+real(conjg(jacobian(axis))*residual,dp)
                    do jaxis = 1, 5
                        hessian(axis,jaxis) = hessian(axis,jaxis)+ &
                            &real(conjg(jacobian(axis))*jacobian(jaxis),dp)
                    enddo
                enddo
            enddo
        enddo
        if( min_switch_margin == huge(0._dp) ) min_switch_margin = 0._dp
    end subroutine pose_normal_terms

    !> Return the joint pose objective and five-vector gradient without exposing
    !! the computed Gauss-Newton block or stencil margin.
    subroutine pose_objective_gradient( self, rotmat, shift, observed, objective, gradient, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(out) :: objective, gradient(5)
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp) :: hessian(5,5), min_switch_margin

        if( present(transfer) )then
            call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian, &
                &min_switch_margin,transfer)
        else
            call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian,min_switch_margin)
        endif
    end subroutine pose_objective_gradient

    !>  \brief  Counts active Fourier samples whose nearest-grid interpolation
    !!          stencil changes between two rotation matrices.
    function count_stencil_switches( self, rotmat, trial_rotmat ) result(nswitches)
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3), trial_rotmat(3,3)
        integer :: nswitches
        real(dp) :: loc(3), trial_loc(3)
        integer :: h, k

        if( .not. self%exists ) error stop 'count_stencil_switches called on an empty Fourier workspace'
        nswitches = 0
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h+k*k < self%sqhp .or. h*h+k*k > self%sqlp ) cycle
                loc = real(self%padf,dp)*matmul(real([h,k,0],dp),rotmat)
                trial_loc = real(self%padf,dp)*matmul(real([h,k,0],dp),trial_rotmat)
                if( any(nint(loc) /= nint(trial_loc)) ) nswitches = nswitches+1
            enddo
        enddo
    end function count_stencil_switches

    !>  \brief  Damped two-parameter Gauss-Newton refinement. Only accepted,
    !!          fully recomputed objective values are appended to the trace.
    subroutine refine_shift_lm( self, rotmat, observed, shift, max_iterations, &
        &accepted_objectives, naccepted, status, nattempted, max_trial_step, transfer )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(in) :: rotmat(3,3)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(inout) :: shift(2)
        integer, intent(in) :: max_iterations
        real(dp), intent(out) :: accepted_objectives(0:)
        integer, intent(out) :: naccepted
        integer, intent(out), optional :: status, nattempted
        real(dp), intent(out), optional :: max_trial_step
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp) :: gradient(2), hessian(2,2), trial_gradient(2), trial_hessian(2,2)
        real(dp) :: solve_matrix(2,2), diagonal(2), direction(2), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, det, predicted, actual, ratio, maxdiag
        real(dp) :: discriminant, lambda_max, lambda_min, step_norm, relative_reduction
        integer :: axis, iteration, outcome, attempted
        logical :: bounded_trial
        if( ubound(accepted_objectives,1) < max_iterations )then
            error stop 'refine_shift_lm objective trace is shorter than max_iterations+1'
        endif
        if( present(transfer) )then
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian,transfer)
        else
            call self%shift_normal_terms(rotmat,shift,observed,objective,gradient,hessian)
        endif
        accepted_objectives = huge(0._dp)
        accepted_objectives(0) = objective
        naccepted = 0
        attempted = 0
        bounded_trial = .false.
        outcome = SHIFT_LM_ITERATION_LIMIT
        if( present(max_trial_step) ) max_trial_step = 0._dp
        if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian)) )then
            outcome = SHIFT_LM_INVALID_NUMERICS
            if( present(status) ) status = outcome
            if( present(nattempted) ) nattempted = attempted
            return
        endif
        mu = 1.e-3_dp
        do iteration = 1, max_iterations
            maxdiag = max(maxval([(hessian(axis,axis),axis=1,2)]),1._dp)
            ! Eigenvalues diagnose whether both shift directions are observable.
            discriminant = sqrt(max(0._dp,(hessian(1,1)-hessian(2,2))**2 + &
                &4._dp*hessian(1,2)*hessian(2,1)))
            lambda_max = 0.5_dp*(hessian(1,1)+hessian(2,2)+discriminant)
            lambda_min = 0.5_dp*(hessian(1,1)+hessian(2,2)-discriminant)
            if( lambda_max <= sqrt(epsilon(1._dp))*maxdiag .or. &
                &lambda_min <= sqrt(epsilon(1._dp))*max(lambda_max,1._dp) )then
                outcome = SHIFT_LM_NO_RELIABLE_UPDATE
                exit
            endif
            if( sqrt(dot_product(gradient,gradient)) < 1.e-8_dp )then
                outcome = merge(SHIFT_LM_ACCEPTED_IMPROVEMENT,SHIFT_LM_FINITE_NO_IMPROVEMENT,naccepted>0)
                exit
            endif
            do axis = 1, 2
                diagonal(axis) = max(hessian(axis,axis),sqrt(epsilon(1._dp))*maxdiag,epsilon(1._dp))
            enddo
            solve_matrix = hessian
            solve_matrix(1,1) = solve_matrix(1,1) + mu*diagonal(1)
            solve_matrix(2,2) = solve_matrix(2,2) + mu*diagonal(2)
            det = solve_matrix(1,1)*solve_matrix(2,2) - solve_matrix(1,2)*solve_matrix(2,1)
            if( abs(det) <= epsilon(1._dp)*maxdiag*maxdiag )then
                outcome = SHIFT_LM_NO_RELIABLE_UPDATE
                exit
            endif
            direction(1) = (-solve_matrix(2,2)*gradient(1) + solve_matrix(1,2)*gradient(2)) / det
            direction(2) = ( solve_matrix(2,1)*gradient(1) - solve_matrix(1,1)*gradient(2)) / det
            if( any(.not. ieee_is_finite(direction)) )then
                outcome = SHIFT_LM_INVALID_NUMERICS
                exit
            endif
            ! Shift coordinates are pixels; cap every trial displacement at one pixel.
            step_norm = sqrt(dot_product(direction,direction))
            if( step_norm > 1._dp )then
                direction = direction/step_norm
                bounded_trial = .true.
            endif
            step_norm = min(step_norm,1._dp)
            if( present(max_trial_step) ) max_trial_step = max(max_trial_step,step_norm)
            predicted = -dot_product(gradient,direction) - 0.5_dp * &
                &dot_product(direction,matmul(hessian,direction))
            if( .not. ieee_is_finite(predicted) )then
                outcome = SHIFT_LM_INVALID_NUMERICS
                exit
            elseif( predicted <= 0._dp )then
                mu = 4._dp * mu
                cycle
            endif
            trial_shift = shift + direction
            if( present(transfer) )then
                call self%shift_normal_terms(rotmat,trial_shift,observed,trial_objective,&
                    &trial_gradient,trial_hessian,transfer)
            else
                call self%shift_normal_terms(rotmat,trial_shift,observed,trial_objective,&
                    &trial_gradient,trial_hessian)
            endif
            attempted = attempted + 1
            if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian)) )then
                mu = 4._dp * mu
                outcome = SHIFT_LM_INVALID_NUMERICS
                cycle
            endif
            actual = objective - trial_objective
            ratio = actual / predicted
            if( actual > 0._dp .and. ratio >= 0.25_dp )then
                relative_reduction = actual/max(abs(objective),1._dp)
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                naccepted = naccepted + 1
                accepted_objectives(naccepted) = objective
                if( ratio > 0.75_dp ) mu = max(mu/2._dp,epsilon(1._dp))
                outcome = SHIFT_LM_ACCEPTED_IMPROVEMENT
                if( step_norm < 1.e-8_dp .or. relative_reduction < 1.e-10_dp ) exit
            else
                mu = 4._dp * mu
            endif
        enddo
        if( outcome == SHIFT_LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
            &outcome = SHIFT_LM_STEP_BOUND_REJECTED
        if( present(status) ) status = outcome
        if( present(nattempted) ) nattempted = attempted
    end subroutine refine_shift_lm

    !>  \brief  Scaled, bounded five-parameter LM refinement for a right
    !!          rotation increment and two image shifts.
    subroutine refine_pose_lm( self, rotmat, observed, shift, rotation_scale, max_iterations, &
        &accepted_objectives, naccepted, status, nattempted, max_rotation_step, &
        &max_shift_step, nstencil_switches, transfer, accepted_rotmats, accepted_shifts, &
        &active_parameters, anchor_rotmat, anchor_shift, max_total_rotation, max_total_shift )
        class(pcg_fourier_workspace), intent(in) :: self
        real(dp), intent(inout) :: rotmat(3,3), shift(2)
        complex, intent(in) :: observed(self%lims2(1,1):self%lims2(1,2),&
                                         &self%lims2(2,1):self%lims2(2,2))
        real(dp), intent(in) :: rotation_scale
        integer, intent(in) :: max_iterations
        real(dp), intent(out) :: accepted_objectives(0:)
        integer, intent(out) :: naccepted, status, nattempted
        real(dp), intent(out) :: max_rotation_step, max_shift_step
        integer, intent(out) :: nstencil_switches
        complex, optional, intent(in) :: transfer(self%lims2(1,1):self%lims2(1,2),&
                                                   &self%lims2(2,1):self%lims2(2,2))
        real(dp), optional, intent(out) :: accepted_rotmats(:,:,0:), accepted_shifts(:,0:)
        logical, optional, intent(in) :: active_parameters(5)
        real(dp), optional, intent(in) :: anchor_rotmat(3,3), anchor_shift(2)
        real(dp), optional, intent(in) :: max_total_rotation, max_total_shift
        real(dp) :: gradient(5), hessian(5,5), trial_gradient(5), trial_hessian(5,5)
        real(dp) :: scaled_gradient(5), scaled_hessian(5,5), solve_matrix(5,5)
        real(dp) :: diagonal(5), coordinate_scale(5), scaled_direction(5), direction(5)
        real(dp) :: trial_rotmat(3,3), trial_shift(2), objective, trial_objective
        real(dp) :: mu, predicted, actual, ratio, rotation_norm, shift_norm, hessian_scale
        real(dp) :: relative_reduction, min_switch_margin, trial_switch_margin
        real(dp) :: cumulative_rotation, cumulative_shift, sine_half
        integer :: axis, jaxis, iteration, trial_switches
        logical :: active(5), bounded_trial, cumulative_guard, reliable

        if( max_iterations < 1 ) error stop 'refine_pose_lm requires at least one LM iteration'
        if( rotation_scale <= 0._dp .or. .not. ieee_is_finite(rotation_scale) ) &
            &error stop 'refine_pose_lm requires a positive finite rotation scale'
        if( ubound(accepted_objectives,1) < max_iterations ) &
            &error stop 'refine_pose_lm objective trace is shorter than max_iterations+1'
        if( present(accepted_rotmats) )then
            if( size(accepted_rotmats,1) /= 3 .or. size(accepted_rotmats,2) /= 3 .or. &
                &ubound(accepted_rotmats,3) < max_iterations ) &
                &error stop 'refine_pose_lm rotation trace has invalid dimensions'
        endif
        if( present(accepted_shifts) )then
            if( size(accepted_shifts,1) /= 2 .or. ubound(accepted_shifts,2) < max_iterations ) &
                &error stop 'refine_pose_lm shift trace has invalid dimensions'
        endif
        active = .true.
        if( present(active_parameters) ) active = active_parameters
        if( .not. any(active) ) error stop 'refine_pose_lm requires one active parameter'
        cumulative_guard = present(anchor_rotmat) .and. present(anchor_shift) .and. &
            &present(max_total_rotation) .and. present(max_total_shift)
        if( cumulative_guard .neqv. (present(anchor_rotmat) .or. present(anchor_shift) .or. &
            &present(max_total_rotation) .or. present(max_total_shift)) ) &
            &error stop 'refine_pose_lm cumulative guard requires all four arguments'
        if( cumulative_guard )then
            if( max_total_rotation <= 0._dp .or. max_total_shift <= 0._dp ) &
                &error stop 'refine_pose_lm cumulative bounds must be positive'
        endif
        if( present(transfer) )then
            call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian, &
                &min_switch_margin,transfer)
        else
            call self%pose_normal_terms(rotmat,shift,observed,objective,gradient,hessian,min_switch_margin)
        endif
        accepted_objectives = huge(0._dp)
        accepted_objectives(0) = objective
        if( present(accepted_rotmats) ) accepted_rotmats(:,:,0) = rotmat
        if( present(accepted_shifts) ) accepted_shifts(:,0) = shift
        naccepted = 0
        nattempted = 0
        nstencil_switches = 0
        max_rotation_step = 0._dp
        max_shift_step = 0._dp
        bounded_trial = .false.
        status = POSE_LM_ITERATION_LIMIT
        if( .not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian)) )then
            status = POSE_LM_INVALID_NUMERICS
            return
        endif

        ! Dimensionless variables balance radians against pixel shifts.
        coordinate_scale = [rotation_scale,rotation_scale,rotation_scale,1._dp,1._dp]
        do axis = 1, 5
            scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
            do jaxis = 1, 5
                scaled_hessian(axis,jaxis) = &
                    &coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
            enddo
        enddo
        call apply_pose_parameter_mask(scaled_gradient,scaled_hessian,active)
        mu = 1.e-3_dp
        do iteration = 1, max_iterations
            ! Damping must not hide an unidentifiable five-parameter block.
            call solve_pose_cholesky(scaled_hessian,-scaled_gradient,scaled_direction,reliable)
            if( .not. reliable )then
                status = merge(POSE_LM_ACCEPTED_IMPROVEMENT,POSE_LM_NO_RELIABLE_UPDATE,naccepted>0)
                exit
            endif
            if( sqrt(dot_product(scaled_gradient,scaled_gradient)) < 1.e-8_dp )then
                status = merge(POSE_LM_ACCEPTED_IMPROVEMENT,POSE_LM_FINITE_NO_IMPROVEMENT,naccepted>0)
                exit
            endif
            hessian_scale = max(maxval(abs(scaled_hessian)),POSE_NUMERIC_FLOOR)
            do axis = 1, 5
                diagonal(axis) = max(scaled_hessian(axis,axis), &
                    &sqrt(epsilon(1._dp))*hessian_scale,POSE_NUMERIC_FLOOR)
            enddo
            solve_matrix = scaled_hessian
            do axis = 1, 5
                solve_matrix(axis,axis) = solve_matrix(axis,axis)+mu*diagonal(axis)
            enddo
            call solve_pose_cholesky(solve_matrix,-scaled_gradient,scaled_direction,reliable)
            if( .not. reliable )then
                status = POSE_LM_NO_RELIABLE_UPDATE
                exit
            endif
            direction = coordinate_scale*scaled_direction
            if( any(.not. ieee_is_finite(direction)) )then
                status = POSE_LM_INVALID_NUMERICS
                exit
            endif
            ! Bound rotations in radians and shifts in pixels independently.
            rotation_norm = sqrt(dot_product(direction(1:3),direction(1:3)))
            if( rotation_norm > rotation_scale )then
                direction(1:3) = direction(1:3)*(rotation_scale/rotation_norm)
                bounded_trial = .true.
            endif
            shift_norm = sqrt(dot_product(direction(4:5),direction(4:5)))
            if( shift_norm > 1._dp )then
                direction(4:5) = direction(4:5)/shift_norm
                bounded_trial = .true.
            endif
            rotation_norm = sqrt(dot_product(direction(1:3),direction(1:3)))
            shift_norm = sqrt(dot_product(direction(4:5),direction(4:5)))
            max_rotation_step = max(max_rotation_step,rotation_norm)
            max_shift_step = max(max_shift_step,shift_norm)
            ! Quadratic LM model: predicted = -g^T d - 1/2 d^T H d.
            predicted = -dot_product(gradient,direction)-0.5_dp* &
                &dot_product(direction,matmul(hessian,direction))
            if( .not. ieee_is_finite(predicted) )then
                status = POSE_LM_INVALID_NUMERICS
                exit
            elseif( predicted <= 0._dp )then
                mu = 4._dp*mu
                cycle
            endif
            trial_rotmat = right_increment_rotation(rotmat,direction(1:3))
            trial_shift = shift+direction(4:5)
            nattempted = nattempted+1
            if( cumulative_guard )then
                sine_half = sqrt(sum((trial_rotmat-anchor_rotmat)**2))/(2._dp*sqrt(2._dp))
                cumulative_rotation = 2._dp*asin(max(0._dp,min(1._dp,sine_half)))
                cumulative_shift = sqrt(sum((trial_shift-anchor_shift)**2))
                if( cumulative_rotation > max_total_rotation+10._dp*epsilon(1._dp) .or. &
                    &cumulative_shift > max_total_shift+10._dp*epsilon(1._dp) )then
                    mu = 4._dp*mu
                    bounded_trial = .true.
                    cycle
                endif
            endif
            trial_switches = self%count_stencil_switches(rotmat,trial_rotmat)
            nstencil_switches = nstencil_switches+trial_switches
            if( present(transfer) )then
                call self%pose_normal_terms(trial_rotmat,trial_shift,observed,trial_objective, &
                    &trial_gradient,trial_hessian,trial_switch_margin,transfer)
            else
                call self%pose_normal_terms(trial_rotmat,trial_shift,observed,trial_objective, &
                    &trial_gradient,trial_hessian,trial_switch_margin)
            endif
            if( .not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian)) )then
                mu = 4._dp*mu
                status = POSE_LM_INVALID_NUMERICS
                cycle
            endif
            ! Gain ratio compares the recomputed reduction with the local model.
            actual = objective-trial_objective
            ratio = actual/predicted
            if( actual > 0._dp .and. ratio >= 0.25_dp )then
                relative_reduction = actual/max(abs(objective),1._dp)
                rotmat = trial_rotmat
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                do axis = 1, 5
                    scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
                    do jaxis = 1, 5
                        scaled_hessian(axis,jaxis) = &
                            &coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
                    enddo
                enddo
                call apply_pose_parameter_mask(scaled_gradient,scaled_hessian,active)
                naccepted = naccepted+1
                accepted_objectives(naccepted) = objective
                if( present(accepted_rotmats) ) accepted_rotmats(:,:,naccepted) = rotmat
                if( present(accepted_shifts) ) accepted_shifts(:,naccepted) = shift
                if( ratio > 0.75_dp ) mu = max(mu/2._dp,epsilon(1._dp))
                status = POSE_LM_ACCEPTED_IMPROVEMENT
                if( max(rotation_norm,shift_norm) < 1.e-8_dp .or. relative_reduction < 1.e-10_dp ) exit
            else
                mu = 4._dp*mu
            endif
        enddo
        if( status == POSE_LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial ) &
            &status = POSE_LM_STEP_BOUND_REJECTED
    end subroutine refine_pose_lm

    !> Freeze inactive pose coordinates while retaining one five-vector LM path.
    pure subroutine apply_pose_parameter_mask( gradient, hessian, active )
        real(dp), intent(inout) :: gradient(5), hessian(5,5)
        logical, intent(in) :: active(5)
        integer :: axis

        do axis = 1, 5
            if( active(axis) ) cycle
            gradient(axis) = 0._dp
            hessian(axis,:) = 0._dp
            hessian(:,axis) = 0._dp
            hessian(axis,axis) = 1._dp
        enddo
    end subroutine apply_pose_parameter_mask

    !>  \brief  Cholesky solve with a relative pivot test for a 5-by-5
    !!          symmetric positive-definite pose block.
    pure subroutine solve_pose_cholesky( matrix, rhs, solution, reliable )
        real(dp), intent(in) :: matrix(5,5), rhs(5)
        real(dp), intent(out) :: solution(5)
        logical, intent(out) :: reliable
        real(dp) :: lower(5,5), intermediate(5), pivot, pivot_floor, matrix_scale
        integer :: i, j

        solution = 0._dp
        intermediate = 0._dp
        lower = 0._dp
        reliable = .false.
        if( any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs)) ) return
        matrix_scale = maxval(abs(matrix))
        if( matrix_scale <= POSE_NUMERIC_FLOOR ) return
        pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
        do i = 1, 5
            do j = 1, i-1
                lower(i,j) = (matrix(i,j)-dot_product(lower(i,1:j-1),lower(j,1:j-1)))/lower(j,j)
            enddo
            pivot = matrix(i,i)-dot_product(lower(i,1:i-1),lower(i,1:i-1))
            if( .not. ieee_is_finite(pivot) .or. pivot <= pivot_floor ) return
            lower(i,i) = sqrt(pivot)
        enddo
        do i = 1, 5
            intermediate(i) = (rhs(i)-dot_product(lower(i,1:i-1),intermediate(1:i-1)))/lower(i,i)
        enddo
        do i = 5, 1, -1
            solution(i) = (intermediate(i)-dot_product(lower(i+1:5,i),solution(i+1:5)))/lower(i,i)
        enddo
        reliable = all(ieee_is_finite(solution))
    end subroutine solve_pose_cholesky

    !>  \brief  G_i F: gathers a full (unpacked) Fourier plane at orientation e
    !!          from the volume most recently passed to set_volume. Plane
    !!          coordinates are NATIVE (h,k); they map onto the oversampled
    !!          lattice as padf*loc, and padsc undoes fwd_ft's 1/product(ldim)
    !!          so the result carries native scale.
    subroutine forward_plane( self, e, plane )
        class(reconstructor_pcg), intent(inout) :: self
        class(ori),                 intent(in)    :: e
        complex,                    intent(out)   :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                            &self%lims2(2,1):self%lims2(2,2))
        complex, allocatable :: cmat(:,:,:)
        real    :: e_rotmat(3,3), loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, i0(3)
        call self%ensure_wimg
        cmat     = self%wimg%get_cmat()
        e_rotmat = e%get_mat()
        plane    = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) private(h,k,loc,i0,w) &
        !$omp schedule(static) proc_bind(close)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                loc = real(self%padf) * matmul(real([h,k,0]), e_rotmat)
                i0  = nint(loc) - self%iwinsz
                call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                plane(h,k) = self%padsc * gather_window(self, cmat, i0, w)
            end do
        end do
        !$omp end parallel do
    end subroutine forward_plane

    !>  \brief  <x_hat, accum> over the oversampled Fourier lattice, for the
    !!          volume currently held by set_volume. Exists so the adjoint
    !!          dot-product test can form <x, G^dagger q> without needing to
    !!          know the lattice or the Friedel packing at all.
    function fourier_dot( self, accum ) result( d )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        complex, allocatable :: cmat(:,:,:)
        real(dp) :: d
        complex  :: xv
        integer  :: h, k, m, ph, pk, pm, ny, nz
        call self%ensure_wimg
        cmat = self%wimg%get_cmat()
        ny   = self%boxpd
        nz   = self%boxpd
        d    = 0.0_dp
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = self%lims3(1,1), self%lims3(1,2)
                    if( h >= 0 )then
                        ph = h + 1
                        pk = k + 1; if( k < 0 ) pk = pk + ny
                        pm = m + 1; if( m < 0 ) pm = pm + nz
                        xv = cmat(ph,pk,pm)
                    else
                        ph = -h + 1
                        pk = -k + 1; if( -k < 0 ) pk = pk + ny
                        pm = -m + 1; if( -m < 0 ) pm = pm + nz
                        xv = conjg(cmat(ph,pk,pm))
                    endif
                    d = d + real(conjg(xv)*accum(h,k,m), dp)
                end do
            end do
        end do
    end function fourier_dot

    !>  \brief  F^dagger G_i^dagger, accumulate form: the literal transpose of
    !!          forward_plane's gather. Each plane entry is scattered once, at
    !!          its own location, with its own KB window -- verified against
    !!          forward_plane by the adjoint dot-product test for arbitrary
    !!          orientations. Written fresh, not derived from
    !!          reconstructor%insert_plane_oversamp or compress_exp.
    subroutine adjoint_plane_add( self, plane, e, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        complex,                    intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                            &self%lims2(2,1):self%lims2(2,2))
        class(ori),                 intent(in)    :: e
        complex,                    intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                &self%lims3(2,1):self%lims3(2,2),&
                                                                &self%lims3(3,1):self%lims3(3,2))
        real    :: e_rotmat(3,3)
        e_rotmat = e%get_mat()
        call scatter_plane(self, plane, e_rotmat, vol_accum)
    end subroutine adjoint_plane_add

    !>  \brief  T_i(h,k) = C_i(h,k) * S_i(h,k) / sqrt(sigma2_i(shell)), the full
    !!          complex per-particle transfer. Needed only
    !!          for the right-hand side: see absT2_plane for why apply_normal
    !!          does not use this.
    function build_transfer( self, ctfparms, shift, sig2arr ) result( T )
        class(reconstructor_pcg), intent(in) :: self
        type(ctfparams),           intent(in) :: ctfparms
        real,                      intent(in) :: shift(2)
        real,            optional, intent(in) :: sig2arr(0:)
        complex   :: T(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real      :: cval, args, sw, sum_df, diff_df, angast, wl, half_wl2_cs, accc, phc, cterm, df, phsh, s2
        integer   :: h, k, shell
        logical   :: l_ctf, l_flip
        ! ctfflag, exactly as image%gen_fplane4rec reads it: CTFFLAG_NO means the
        ! images carry no CTF to model, and CTFFLAG_FLIP means they are ALREADY
        ! phase-flipped, so the signed CTF would reintroduce the very phase the
        ! flip removed. Ignoring the flag is invisible on
        ! synthetic fixtures, which are always generated with a signed CTF, and
        ! destroys a real reconstruction on any phase-flipped project.
        l_ctf  = ctfparms%ctfflag /= CTFFLAG_NO
        l_flip = ctfparms%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
            ! see absT2_plane: flat, call-free CTF form of ft_map_ctf_kernel,
            ! inlined to run over the both-sign-h disk via the LUTs.
            ctfvals     = tfun%get_ctfvars(ctfparms%phshift)
            wl          = ctfvals%wl
            half_wl2_cs = 0.5 * wl * wl * ctfvals%cs
            sum_df      = ctfvals%dfx + ctfvals%dfy
            diff_df     = ctfvals%dfx - ctfvals%dfy
            angast      = ctfvals%angast
            accc        = ctfvals%amp_contr_const
            phc         = ctfvals%phshift
        endif
        T = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(h,k,cval,args,shell,sw,cterm,df,phsh,s2) schedule(static) proc_bind(close)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                cval = 1.0
                if( l_ctf )then
                    s2    = self%spafreqsq_lut(h,k)
                    cterm = cos( 2.0 * (self%ang_lut(h,k) - angast) )
                    df    = 0.5 * ( sum_df + cterm * diff_df )
                    phsh  = PI * wl * s2 * (df - half_wl2_cs * s2)
                    cval  = sin( phsh + phc + accc )
                    if( l_flip ) cval = abs(cval)
                endif
                ! SIGN. apply_adjoint_all forms conjg(T)*y, and production applies
                ! exp(i*(-shift)*shconst*(h,k)) to y (gen_fplane4rec: pshift =
                ! -shift*shconst). So T must carry the POSITIVE phase for its
                ! conjugate to reproduce production's correction. A negated phase
                ! here would shift every particle by twice its own displacement in
                ! the wrong direction -- invisible on synthetic
                ! data, where the same build_transfer generates the observations.
                args       = 2.0*PI * (real(h)*shift(1) + real(k)*shift(2)) / real(self%box)
                sw        = 1.0
                if( present(sig2arr) )then
                    shell = min(self%shell_lut(h,k), ubound(sig2arr,1))
                    sw    = 1.0 / sqrt(sig2arr(shell))
                endif
                T(h,k) = cval * cmplx(cos(args), sin(args)) * sw
            end do
        end do
        !$omp end parallel do
    end function build_transfer

    !>  \brief  Applies the production PCG inverse-noise amplitude to one raw
    !!          Fourier observation: y_w = y / sqrt(sigma2(shell)).
    function whiten_observation( self, plane, sig2arr ) result( whitened )
        class(reconstructor_pcg), intent(in) :: self
        complex, intent(in) :: plane(self%lims2(1,1):self%lims2(1,2), &
                                      &self%lims2(2,1):self%lims2(2,2))
        real, intent(in) :: sig2arr(0:)
        complex :: whitened(self%lims2(1,1):self%lims2(1,2), &
                             &self%lims2(2,1):self%lims2(2,2))
        real :: sigma2
        integer :: h, k, shell
        whitened = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                shell = min(self%shell_lut(h,k),ubound(sig2arr,1))
                sigma2 = sig2arr(shell)
                if( .not. ieee_is_finite(sigma2) .or. sigma2 <= 0.0 )then
                    error stop 'whiten_observation requires finite positive sigma2'
                endif
                whitened(h,k) = plane(h,k) / sqrt(sigma2)
            enddo
        enddo
    end function whiten_observation

    !>  \brief  2D, non-interpolated analog of the gather: reads a real
    !!          particle image's own Fourier plane directly into the lims2 disk
    !!          (the image's own box IS the native grid -- no KB window, no
    !!          wrap). img2d must already be FFT'd. This is how y_planes is
    !!          built for REAL particles; forward_plane only builds SYNTHETIC
    !!          ones by projecting a known volume.
    function extract_native_plane( self, img2d ) result( plane )
        class(reconstructor_pcg), intent(in) :: self
        class(image),               intent(in) :: img2d
        complex :: plane(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2))
        integer :: h, k, phys(3)
        plane = cmplx(0.,0.)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                phys = img2d%comp_addr_phys(h,k,0)
                plane(h,k) = img2d%get_fcomp([h,k,0], phys)
            end do
        end do
    end function extract_native_plane

    !>  \brief  deterministic double-precision real-volume dot product.
    pure function dot_real_volume( self, a, b ) result( d )
        class(reconstructor_pcg), intent(in) :: self
        real,                      intent(in) :: a(self%box,self%box,self%box)
        real,                      intent(in) :: b(self%box,self%box,self%box)
        real(dp) :: d
        d = sum(real(a,dp) * real(b,dp))
    end function dot_real_volume

    ! HIGH-LEVEL OPERATOR

    !>  \brief  The operator the SOLVER sees: P H P when a support constraint is
    !!          set (see set_mask), plain H otherwise. The two concrete operators
    !!          below are deliberately left unmasked -- the test commanders drive
    !!          them directly and compare them against one another, and a mask
    !!          would silently hide any disagreement outside its support.
    function apply_normal( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: p(self%box,self%box,self%box)
        real, allocatable :: hp(:,:,:), pm(:,:,:)
        if( self%l_mask )then
            allocate(pm(self%box,self%box,self%box), source=p)
            call self%mask_mul(pm)
        endif
        select case( self%op_mode )
            case( PCG_OP_KERNEL )
                if( self%l_mask )then
                    hp = self%apply_normal_kernel(pm)
                else
                    hp = self%apply_normal_kernel(p)
                endif
            case default
                if( self%l_mask )then
                    hp = self%apply_normal_matrixfree(pm)
                else
                    hp = self%apply_normal_matrixfree(p)
                endif
        end select
        call self%mask_mul(hp)
    end function apply_normal

    !>  \brief  H p = sum_i G_i^dagger |T_i|^2 G_i p + lambda*p -- the exact
    !!          reference operator.
    !!
    !!          Gather, weight and scatter are FUSED into one pass over the
    !!          plane so the rotated coordinate, the wrapped window indices and
    !!          the KB weights are computed once and used for both directions,
    !!          rather than twice.
    !!
    !!          OpenMP: the particle loop is walked in lockstep by all threads,
    !!          with the plane's h loop worksharing inside it. The scatter uses
    !!          the h-strided colouring scheme (see `stride` in new) so threads
    !!          write disjoint 27-voxel footprints into the single shared
    !!          accumulator -- chosen over per-thread accumulators because at
    !!          box 256 those cost ~135 MB per thread.
    function apply_normal_matrixfree( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: p(self%box,self%box,self%box)
        real,    allocatable :: hp(:,:,:), pd(:,:,:)
        complex, allocatable :: vol_accum(:,:,:), cmat(:,:,:)
        real,    allocatable :: absT2(:,:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3)
        complex :: comp
        integer :: i, g, h, k, l, i0(3)
        integer(timer_int_kind) :: tp
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; apply_normal_matrixfree')
        call self%ensure_wimg
        ! A = A~ . E^-1 : deapodize going in, and again coming out of the
        ! adjoint, so the operator is the true D.S rather than D.S.E
        allocate(pd(self%box,self%box,self%box), source=p)
        call self%deapod_mul(pd)
        if( self%l_profile ) tp = tic()
        call self%set_volume(pd)
        if( self%l_profile ) self%t_setvol = self%t_setvol + real(toc(tp),dp)
        ! NOTE: get_cmat returns a COPY. On the padded lattice that is
        ! (boxpd/2+1)*boxpd^2 complex -- 539 MB at box 256 -- allocated and
        ! streamed on every call. t_cmatcp is what separates that traffic from
        ! the transforms themselves; image offers get_rmat_ptr but no cmat
        ! equivalent, so removing it would mean adding one.
        if( self%l_profile ) tp = tic()
        cmat = self%wimg%get_cmat()
        if( self%l_profile ) self%t_cmatcp = self%t_cmatcp + real(toc(tp),dp)
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        allocate(absT2(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        if( self%l_profile ) tp = tic()
        !$omp parallel default(shared) private(i,g,h,k,l,loc,i0,w,comp,rot) proc_bind(close)
        do i = 1, self%nptcls
            ! per-particle real weight |T_i|^2, shared across threads. The
            ! worksharing inside absT2_plane is orphaned and binds to THIS
            ! region, so the CTF evaluation is spread across the team; its
            ! trailing barrier is what makes absT2 safe to read below.
            call self%absT2_plane(i, absT2)
            ! coordinate replication: gather and scatter at every R_i.S_g, with g
            ! outside the h-strided colour sweep so the reference matches the
            ! kernel built by accumulate_absT2. symmats(:,:,1)=I gives the c1 op.
            do g = 1, self%nsym
                rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
                ! fused gather -> weight -> scatter, h-strided for scatter safety
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                        do k = self%lims2(2,1), self%lims2(2,2)
                            if( h*h + k*k > self%sqlp ) cycle
                            loc  = real(self%padf) * matmul(real([h,k,0]), rot)
                            i0   = nint(loc) - self%iwinsz
                            if( win_wraps(self, i0) ) cycle  ! rim: serial pass below
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            comp = self%padsc * gather_window(self, cmat, i0, w)
                            comp = comp * absT2(h,k) * self%padsc
                            call scatter_window(self, i0, w, comp, vol_accum)
                        end do
                    end do
                    !$omp end do
                end do
                ! wrapping rim, serialized: the colouring's separation guarantee does
                ! not survive folding, see win_wraps
                !$omp single
                do h = self%lims2(1,1), self%lims2(1,2)
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp   ) cycle
                        if( h*h + k*k <= self%sq_rim ) cycle   ! provably cannot wrap
                        loc  = real(self%padf) * matmul(real([h,k,0]), rot)
                        i0   = nint(loc) - self%iwinsz
                        if( .not. win_wraps(self, i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        comp = self%padsc * gather_window(self, cmat, i0, w)
                        comp = comp * absT2(h,k) * self%padsc
                        call scatter_window(self, i0, w, comp, vol_accum)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        if( self%l_profile ) self%t_ploop = self%t_ploop + real(toc(tp),dp)
        hp = self%fold_and_ifft(vol_accum)
        call self%deapod_mul(hp)
        if( self%l_ml_prior ) hp = hp + self%apply_fourier_diagonal(p, self%ml_prior)
        if( self%l_nu_prior ) hp = hp + self%lambda_nu * self%apply_nu_precision(p)
        hp = hp + self%lambda * p
    end function apply_normal_matrixfree

    !>  \brief  kernelized (Toeplitz/Gram) normal operator:
    !!          H_data p = crop(IFFT(Khat * FFT(pad(p)))), cost O(box^3 log box)
    !!          and INDEPENDENT of particle count. See build_kernel for how Khat
    !!          is constructed and why it is not a literal impulse
    !!          response.
    function apply_normal_kernel( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: p(self%box,self%box,self%box)
        real,    allocatable :: hp(:,:,:), work(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        real    :: kv
        integer :: cdim(3), i, j, k
        integer(timer_int_kind) :: tp
        if( .not. self%l_kernel ) THROW_HARD('build_kernel has not been called; apply_normal_kernel')
        call self%ensure_wimg
        allocate(work(self%box,self%box,self%box), source=p)
        ! Khat represents the bare Toeplitz operator T. With deapodization ON
        ! the matrix-free operator is E^-1(E T E)E^-1 = T as well, so nothing
        ! more is needed. With it OFF the matrix-free operator is E T E, so the
        ! envelope has to be reinstated on both sides for the two to agree.
        if( .not. self%l_deapod ) work = work * self%env
        if( self%l_profile ) tp = tic()
        call self%wimg%set_rmat(self%pad_vol(work), .false.)
        call self%wimg%fft()
        if( self%l_profile ) self%t_setvol = self%t_setvol + real(toc(tp),dp)
        if( self%l_profile ) tp = tic()
        cmat = self%wimg%get_cmat()
        if( self%l_profile ) self%t_cmatcp = self%t_cmatcp + real(toc(tp),dp)
        cdim = self%wimg%get_array_shape()
        if( self%l_profile ) tp = tic()
        !$omp parallel do collapse(3) default(shared) private(i,j,k,kv) &
        !$omp schedule(static) proc_bind(close)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    kv = self%Khat(i,j,k)
                    if( self%l_ml_prior  .and. self%l_deapod ) kv = kv + self%ml_prior(i,j,k)
                    cmat(i,j,k) = cmat(i,j,k) * kv
                end do
            end do
        end do
        !$omp end parallel do
        if( self%l_profile ) self%t_khat = self%t_khat + real(toc(tp),dp)
        if( self%l_profile ) tp = tic()
        call self%wimg%set_cmat(cmat)
        call self%wimg%ifft()
        hp = self%crop_vol(self%wimg%get_rmat())
        if( self%l_profile ) self%t_fold = self%t_fold + real(toc(tp),dp)
        if( .not. self%l_deapod ) hp = hp * self%env
        if( self%l_ml_prior .and. .not. self%l_deapod )then
            hp = hp + self%apply_fourier_diagonal(p, self%ml_prior)
        endif
        ! band precision on the same (deapodized) domain as the iterate;
        ! attachment mode is enforced upstream, see assert_prior_attachment_mode
        if( self%l_nu_prior ) hp = hp + self%lambda_nu * self%apply_nu_precision(p)
        hp = hp + self%lambda * p
    end function apply_normal_kernel

    !> Apply C^T F^-1 diag(d) F C on the same padded lattice as Khat. This is
    !! used directly by the matrix-free oracle and by the deapodization-off
    !! diagnostic path; production kernel solves fuse d with Khat above.
    function apply_fourier_diagonal( self, p, diag ) result( q )
        class(reconstructor_pcg), intent(inout) :: self
        real, intent(in) :: p(self%box,self%box,self%box)
        real, intent(in) :: diag(:,:,:)
        real, allocatable :: q(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        integer :: cdim(3), i, j, k
        call self%ensure_wimg
        cdim = self%wimg%get_array_shape()
        if( any(shape(diag) /= cdim) ) THROW_HARD('PCG Fourier diagonal shape mismatch')
        call self%wimg%set_rmat(self%pad_vol(p), .false.)
        call self%wimg%fft()
        cmat = self%wimg%get_cmat()
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    cmat(i,j,k) = cmat(i,j,k) * diag(i,j,k)
                enddo
            enddo
        enddo
        !$omp end parallel do
        call self%wimg%set_cmat(cmat)
        call self%wimg%ifft()
        q = self%crop_vol(self%wimg%get_rmat())
    end function apply_fourier_diagonal

    !>  \brief  b = sum_i G_i^dagger(conjg(T_i) * y_i / sqrt(sigma2_i)), the
    !!          data right-hand side (no prior term). Unlike H, the RHS DOES
    !!          need the full complex T_i including the shift phase.
    function apply_adjoint_all( self, y_planes ) result( b )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                               &self%lims2(2,1):self%lims2(2,2), *)
        real,    allocatable :: b(:,:,:)
        complex, allocatable :: vol_accum(:,:,:), weighted(:,:), T(:,:)
        integer :: i, h, k, shell, R
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; apply_adjoint_all')
        R = self%lims2(1,2)
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        allocate(weighted(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(T(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        do i = 1, self%nptcls
            if( self%l_use_ctf )then
                call self%transfer_plane_cmplx(i, T)
                !$omp parallel do collapse(2) default(shared) private(h,k,shell) &
                !$omp schedule(static) proc_bind(close)
                do k = self%lims2(2,1), self%lims2(2,2)
                    do h = self%lims2(1,1), self%lims2(1,2)
                        if( h*h + k*k > self%sqlp )then
                            weighted(h,k) = cmplx(0.,0.)
                        else
                            shell = min(nint(sqrt(real(h*h+k*k))), R)
                            weighted(h,k) = conjg(T(h,k)) * y_planes(h,k,i) / sqrt(self%sig2(shell,i))
                        endif
                    end do
                end do
                !$omp end parallel do
            else
                weighted = y_planes(:,:,i)
            endif
            call scatter_plane(self, weighted, self%rotmats(:,:,i), vol_accum)
        end do
        b = self%fold_and_ifft(vol_accum)
        ! A^dagger = E^-1 A~^dagger : the RHS gets the correction once, the
        ! normal operator twice (once per application of the gather/scatter)
        call self%deapod_mul(b)
    end function apply_adjoint_all

    !>  \brief  Batch form of apply_adjoint_all's loop body: scatters one batch
    !!          of observed planes into a Fourier accumulator WITHOUT folding.
    !!
    !!          Folding is deliberately not done here. fold_and_ifft is an
    !!          inverse FFT followed by a centre-crop, and neither commutes with
    !!          summation over batches once the deapodization is applied -- so
    !!          the fold happens exactly once, in end_accum, over the completed
    !!          accumulator.
    subroutine accumulate_rhs( self, y_batch, nb, ifrom, acc )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                    intent(in)    :: nb, ifrom
        complex,                    intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                              &self%lims2(2,1):self%lims2(2,2), nb)
        complex,                    intent(inout) :: acc(self%lims3(1,1):self%lims3(1,2),&
                                                          &self%lims3(2,1):self%lims3(2,2),&
                                                          &self%lims3(3,1):self%lims3(3,2))
        complex, allocatable :: weighted(:,:), T(:,:)
        integer :: ib, i, h, k, shell, R
        if( nb < 1 ) return
        R = self%lims2(1,2)
        allocate(weighted(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(T(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        do ib = 1, nb
            i = ifrom + ib - 1
            if( self%l_use_ctf )then
                call self%transfer_plane_cmplx(i, T)
                !$omp parallel do collapse(2) default(shared) private(h,k,shell) &
                !$omp schedule(static) proc_bind(close)
                do k = self%lims2(2,1), self%lims2(2,2)
                    do h = self%lims2(1,1), self%lims2(1,2)
                        if( h*h + k*k > self%sqlp )then
                            weighted(h,k) = cmplx(0.,0.)
                        else
                            shell = min(nint(sqrt(real(h*h+k*k))), R)
                            weighted(h,k) = conjg(T(h,k)) * y_batch(h,k,ib) / sqrt(self%sig2(shell,i))
                        endif
                    end do
                end do
                !$omp end parallel do
            else
                weighted = y_batch(:,:,ib)
            endif
            call scatter_plane(self, weighted, self%rotmats(:,:,i), acc)
        end do
    end subroutine accumulate_rhs

    ! STREAMING SETUP
    !
    ! begin_accum / accumulate_batch / end_accum replace the pattern of holding
    ! every particle's Fourier plane resident for the whole run. The images are
    ! needed for one thing only -- forming the RHS -- and for one pass only, so
    ! a caller reads a batch, hands it over, and discards it.
    !
    ! Peak memory during accumulation is one complex RHS and one real density
    ! accumulator plus the work image; that floor is constant in nptcls. See
    ! doc/policies/reconstruct3D_pcg_policy.md.

    subroutine begin_accum( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; begin_accum')
        call self%begin_reduction
    end subroutine begin_accum

    !> Allocate zero raw accumulators without particle metadata. Distributed
    !! masters use this entry point before adding worker artifacts; workers use
    !! begin_accum after prep_particles as before.
    subroutine begin_reduction( self )
        class(reconstructor_pcg), intent(inout) :: self
        call self%ensure_wimg
        if( allocated(self%acc_work) ) deallocate(self%acc_work)
        if( allocated(self%b_work)   ) deallocate(self%b_work)
        if( allocated(self%b_rhs)    ) deallocate(self%b_rhs)
        allocate(self%acc_work(self%lims3(1,1):self%lims3(1,2),&
                              &self%lims3(2,1):self%lims3(2,2),&
                              &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        allocate(self%b_work(self%lims3(1,1):self%lims3(1,2),&
                            &self%lims3(2,1):self%lims3(2,2),&
                            &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        self%l_accum = .true.
        self%l_rhs   = .false.
        self%reduction_next_part = 1
        self%reduction_nparts    = 0
        self%reduction_state     = 0
        self%reduction_eo        = -1
    end subroutine begin_reduction

    !>  \brief  Accumulates one batch into BOTH accumulators. y_batch holds the
    !!          observed planes for particles ifrom .. ifrom+nb-1 of the
    !!          selection passed to prep_particles.
    subroutine accumulate_batch( self, y_batch, nb, ifrom )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                    intent(in)    :: nb, ifrom
        complex,                    intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                              &self%lims2(2,1):self%lims2(2,2), nb)
        complex, allocatable :: bacc(:,:,:)
        real,    allocatable :: dacc(:,:,:)
        if( .not. self%l_accum ) THROW_HARD('begin_accum has not been called; accumulate_batch')
        if( nb < 1 ) return
        if( ifrom < 1 .or. ifrom + nb - 1 > self%nptcls )then
            THROW_HARD('batch range outside the particle selection; accumulate_batch')
        endif
        ! move_alloc detaches the accumulator from `self` before it is passed as
        ! a dummy. Handing self%b_work straight to a procedure that also takes
        ! `self` as intent(inout) makes the two aliases of each other, which the
        ! standard permits a compiler to assume cannot happen -- it may then
        ! cache or reorder the writes. move_alloc is an O(1) descriptor swap, so
        ! removing the hazard costs nothing and beats betting on the optimizer.
        call move_alloc(self%b_work, bacc)
        call move_alloc(self%acc_work, dacc)
        call self%accumulate_rhs_density(y_batch, nb, ifrom, bacc, dacc)
        call move_alloc(bacc, self%b_work)
        call move_alloc(dacc, self%acc_work)
    end subroutine accumulate_batch

    !> Atomically publish one worker's raw full-range B and D accumulators.
    !! No folding, deapodization, shell flooring, kernel finalization or solve
    !! is permitted before this write. The embedded header is the manifest for
    !! this single-file artifact; promotion from .tmp happens only after close.
    subroutine write_raw_accum( self, fname, state, eo, part, nparts, nptcls, provenance )
        class(reconstructor_pcg), intent(in) :: self
        class(string),             intent(in) :: fname
        integer,                   intent(in) :: state, eo, part, nparts, nptcls
        character(len=*),          intent(in) :: provenance
        type(string) :: tmpfname
        character(len=PCG_RAW_PROV_LEN) :: prov_fixed
        integer :: funit, ierr, m
        integer(int64) :: file_size
        if( state < 1 .or. eo < 0 .or. eo > 1 ) THROW_HARD('invalid raw PCG state or half')
        if( part < 1 .or. nparts < part ) THROW_HARD('invalid raw PCG part index')
        if( nptcls < 0 ) THROW_HARD('invalid raw PCG particle count')
        if( nptcls > 0 .and. .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open')
        prov_fixed = ' '
        if( len_trim(provenance) > 0 )then
            prov_fixed(1:min(len_trim(provenance),PCG_RAW_PROV_LEN)) = &
                &provenance(1:min(len_trim(provenance),PCG_RAW_PROV_LEN))
        endif
        tmpfname = fname//'.tmp'
        call del_file(tmpfname)
        call fopen(funit, file=tmpfname, status='REPLACE', action='WRITE', &
            &access='STREAM', iostat=ierr)
        call fileiochk('write_raw_accum opening temporary artifact', ierr)
        write(funit, iostat=ierr) PCG_RAW_ACCUM_MAGIC, PCG_RAW_ACCUM_VERSION
        call fileiochk('write_raw_accum writing magic', ierr)
        write(funit, iostat=ierr) state, eo, part, nparts, nptcls
        call fileiochk('write_raw_accum writing identity', ierr)
        write(funit, iostat=ierr) self%box, self%boxpd, self%padf, self%lims3, self%smpd
        call fileiochk('write_raw_accum writing geometry', ierr)
        write(funit, iostat=ierr) prov_fixed
        call fileiochk('write_raw_accum writing provenance', ierr)
        if( nptcls > 0 )then
            do m = self%lims3(3,1), self%lims3(3,2)
                write(funit, iostat=ierr) self%b_work(:,:,m)
                if( ierr /= 0 ) exit
            enddo
            call fileiochk('write_raw_accum writing B', ierr)
            do m = self%lims3(3,1), self%lims3(3,2)
                write(funit, iostat=ierr) self%acc_work(:,:,m)
                if( ierr /= 0 ) exit
            enddo
            call fileiochk('write_raw_accum writing D', ierr)
        endif
        call fclose(funit)
        file_size = -1_int64
        inquire(file=tmpfname%to_char(), size=file_size)
        if( file_size <= 0_int64 ) THROW_HARD('raw PCG artifact write produced an empty file')
        call simple_rename(tmpfname, fname, overwrite=.true.)
        call tmpfname%kill
    end subroutine write_raw_accum

    !> Add one raw worker artifact to the master's open reduction. Calls must
    !! arrive in ascending part order; this makes the floating-point association
    !! order explicit and reproducible. Slices are streamed to keep temporary
    !! memory negligible relative to the full-range master accumulators.
    subroutine add_raw_accum( self, fname, state, eo, part, nparts, provenance, nptcls )
        class(reconstructor_pcg), intent(inout) :: self
        class(string),             intent(in)    :: fname
        integer,                   intent(in)    :: state, eo, part, nparts
        character(len=*),          intent(in)    :: provenance
        integer,                   intent(out)   :: nptcls
        character(len=16) :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        complex, allocatable :: bslice(:,:)
        real,    allocatable :: dslice(:,:)
        integer :: funit, ierr, m, version, state_file, eo_file, part_file
        integer :: nparts_file, box_file, boxpd_file, padf_file, lims_file(3,2)
        real    :: smpd_file
        integer(int64) :: file_size, stream_pos
        if( .not. self%l_accum ) THROW_HARD('master PCG reduction is not open')
        if( part /= self%reduction_next_part ) THROW_HARD('raw PCG parts must be reduced in ascending order')
        if( .not. file_exists(fname) ) THROW_HARD('raw PCG accumulator artifact is missing')
        prov_expected = ' '
        if( len_trim(provenance) > 0 )then
            prov_expected(1:min(len_trim(provenance),PCG_RAW_PROV_LEN)) = &
                &provenance(1:min(len_trim(provenance),PCG_RAW_PROV_LEN))
        endif
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('add_raw_accum opening artifact', ierr)
        read(funit, iostat=ierr) magic, version
        call fileiochk('add_raw_accum reading magic', ierr)
        read(funit, iostat=ierr) state_file, eo_file, part_file, nparts_file, nptcls
        call fileiochk('add_raw_accum reading identity', ierr)
        read(funit, iostat=ierr) box_file, boxpd_file, padf_file, lims_file, smpd_file
        call fileiochk('add_raw_accum reading geometry', ierr)
        read(funit, iostat=ierr) prov_file
        call fileiochk('add_raw_accum reading provenance', ierr)
        if( magic /= PCG_RAW_ACCUM_MAGIC .or. version /= PCG_RAW_ACCUM_VERSION )then
            THROW_HARD('raw PCG accumulator format mismatch')
        endif
        if( state_file /= state .or. eo_file /= eo .or. part_file /= part )then
            THROW_HARD('raw PCG accumulator identity mismatch')
        endif
        if( nparts_file /= nparts .or. nptcls < 0 ) THROW_HARD('raw PCG accumulator partition mismatch')
        if( box_file /= self%box .or. boxpd_file /= self%boxpd .or. padf_file /= self%padf )then
            THROW_HARD('raw PCG accumulator box mismatch')
        endif
        if( any(lims_file /= self%lims3) ) THROW_HARD('raw PCG accumulator lattice mismatch')
        if( abs(smpd_file-self%smpd) > max(1.0e-6,1.0e-6*abs(self%smpd)) )then
            THROW_HARD('raw PCG accumulator sampling mismatch')
        endif
        if( prov_file /= prov_expected ) THROW_HARD('raw PCG accumulator provenance mismatch')
        if( self%reduction_nparts > 0 )then
            if( self%reduction_nparts /= nparts .or. self%reduction_state /= state .or. &
                &self%reduction_eo /= eo ) THROW_HARD('raw PCG reduction set mismatch')
        endif
        if( nptcls > 0 )then
            allocate(bslice(self%lims3(1,1):self%lims3(1,2), self%lims3(2,1):self%lims3(2,2)))
            allocate(dslice(self%lims3(1,1):self%lims3(1,2), self%lims3(2,1):self%lims3(2,2)))
            do m = self%lims3(3,1), self%lims3(3,2)
                read(funit, iostat=ierr) bslice
                if( ierr /= 0 ) exit
                self%b_work(:,:,m) = self%b_work(:,:,m) + bslice
            enddo
            call fileiochk('add_raw_accum reading B', ierr)
            do m = self%lims3(3,1), self%lims3(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
                self%acc_work(:,:,m) = self%acc_work(:,:,m) + dslice
            enddo
            call fileiochk('add_raw_accum reading D', ierr)
            deallocate(bslice, dslice)
        endif
        file_size  = -1_int64
        stream_pos = -1_int64
        inquire(unit=funit, size=file_size, pos=stream_pos)
        call fclose(funit)
        if( stream_pos-1_int64 /= file_size ) THROW_HARD('raw PCG accumulator has trailing or missing bytes')
        self%reduction_nparts    = nparts
        self%reduction_state     = state
        self%reduction_eo        = eo
        self%reduction_next_part = part + 1
    end subroutine add_raw_accum

    !> Scale an open raw accumulator in the sufficient-statistics domain.
    !! This is the only valid place to apply the u/f and (1-u) continuation
    !! weights: neither finalized kernels nor reconstructed maps are additive.
    subroutine scale_raw_accum( self, weight )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: weight
        if( .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open; scale_raw_accum')
        if( .not. ieee_is_finite(weight) .or. weight < 0.0 ) THROW_HARD('invalid raw PCG accumulator weight')
        self%b_work   = weight * self%b_work
        self%acc_work = weight * self%acc_work
    end subroutine scale_raw_accum

    !> Relative L2 differences between two open raw accumulator chains. This
    !! avoids materializing diagnostic copies of the padded B and D lattices.
    subroutine compare_raw_accum( self, other, b_relerr, d_relerr )
        class(reconstructor_pcg), intent(in) :: self, other
        real,                     intent(out) :: b_relerr, d_relerr
        real(dp) :: bnum, bden, dnum, dden
        integer  :: h, k, m
        if( .not. self%l_accum .or. .not. other%l_accum ) THROW_HARD('raw PCG comparison requires open accumulators')
        if( self%box /= other%box .or. self%boxpd /= other%boxpd .or. any(self%lims3 /= other%lims3) )then
            THROW_HARD('raw PCG comparison geometry mismatch')
        endif
        bnum = 0.0_dp
        bden = 0.0_dp
        dnum = 0.0_dp
        dden = 0.0_dp
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = self%lims3(1,1), self%lims3(1,2)
                    bnum = bnum + real(abs(self%b_work(h,k,m)-other%b_work(h,k,m))**2,dp)
                    bden = bden + real(abs(other%b_work(h,k,m))**2,dp)
                    dnum = dnum + real((self%acc_work(h,k,m)-other%acc_work(h,k,m))**2,dp)
                    dden = dden + real(other%acc_work(h,k,m)**2,dp)
                enddo
            enddo
        enddo
        b_relerr = real(sqrt(bnum / max(1.0_dp,bden)))
        d_relerr = real(sqrt(dnum / max(1.0_dp,dden)))
    end subroutine compare_raw_accum

    !> Non-destructive header compatibility check for a persisted raw
    !! accumulator artifact against the current reconstruction geometry and
    !! continuation identity. Constant-FOV crop continuity applies: a chain
    !! persisted at a SMALLER crop with the same physical field of view is
    !! compatible (add_raw_accum_weighted embeds it by index-aligned
    !! zero-extension). A larger-crop chain, a field-of-view change, an
    !! identity/provenance change, or any read failure counts as incompatible:
    !! callers must discard such a chain and re-seed (bootstrap) instead of
    !! reducing it.
    logical function pcg_raw_accum_compatible( fname, box, smpd, provenance ) result( l_compatible )
        class(string),    intent(in) :: fname
        integer,          intent(in) :: box
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: provenance
        character(len=16) :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        integer :: funit, ierr, version, state_file, eo_file, part_file, nparts_file, nptcls_file
        integer :: box_file, boxpd_file, padf_file, lims_file(3,2)
        real    :: smpd_file, fov_file, fov_cur
        l_compatible = .false.
        if( .not. file_exists(fname) ) return
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
        if( ierr /= 0 ) return
        read(funit, iostat=ierr) magic, version
        if( ierr == 0 ) read(funit, iostat=ierr) state_file, eo_file, part_file, nparts_file, nptcls_file
        if( ierr == 0 ) read(funit, iostat=ierr) box_file, boxpd_file, padf_file, lims_file, smpd_file
        if( ierr == 0 ) read(funit, iostat=ierr) prov_file
        call fclose(funit)
        if( ierr /= 0 ) return
        if( magic /= PCG_RAW_ACCUM_MAGIC .or. version /= PCG_RAW_ACCUM_VERSION ) return
        if( box_file > box .or. padf_file /= OSMPL_PAD_FAC .or. boxpd_file /= padf_file*box_file ) return
        fov_file = real(box_file) * smpd_file
        fov_cur  = real(box) * smpd
        if( abs(fov_file-fov_cur) > 1.0e-4*fov_cur ) return
        prov_expected = ' '
        if( len_trim(provenance) > 0 )then
            prov_expected(1:min(len_trim(provenance),PCG_RAW_PROV_LEN)) = &
                &provenance(1:min(len_trim(provenance),PCG_RAW_PROV_LEN))
        endif
        if( prov_file /= prov_expected ) return
        l_compatible = .true.
    end function pcg_raw_accum_compatible

    !> Add one complete raw artifact with an explicit continuation weight.
    !! Unlike add_raw_accum this routine does not participate in worker-part
    !! ordering. It is for deterministic algebra on already reduced current
    !! and previous chains after the worker reduction has completed.
    !!
    !! Constant-FOV crop continuity: under box*smpd == box_crop*smpd_crop the
    !! padded lattices of consecutive crops share their frequency step, so a
    !! SMALLER previous grid is an index-aligned central subset of the current
    !! one and its (B,D) statistics embed exactly by zero-extension -- the same
    !! autoscale ramp the gridding trailing chain performs. The one
    !! approximation is the old lattice's wrap rim: KB windows that wrapped
    !! around the old period carry aliased mass in the outermost old shells,
    !! at/beyond the producing stage's matching band and decaying as (1-u)^k.
    !! Larger-than-current grids cannot be restricted and are rejected.
    subroutine add_raw_accum_weighted( self, fname, state, eo, part, nparts, provenance, weight, nptcls )
        class(reconstructor_pcg), intent(inout) :: self
        class(string),            intent(in)    :: fname
        integer,                  intent(in)    :: state, eo, part, nparts
        character(len=*),         intent(in)    :: provenance
        real,                     intent(in)    :: weight
        integer,                  intent(out)   :: nptcls
        character(len=16) :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        complex, allocatable :: bslice(:,:)
        real,    allocatable :: dslice(:,:)
        integer :: funit, ierr, m, version, state_file, eo_file, part_file
        integer :: nparts_file, box_file, boxpd_file, padf_file, lims_file(3,2)
        real    :: smpd_file, fov_file, fov_self
        integer(int64) :: file_size, stream_pos
        if( .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open; add_raw_accum_weighted')
        if( .not. ieee_is_finite(weight) .or. weight < 0.0 ) THROW_HARD('invalid weighted raw PCG contribution')
        if( .not. file_exists(fname) ) THROW_HARD('weighted raw PCG accumulator artifact is missing')
        prov_expected = ' '
        if( len_trim(provenance) > 0 )then
            prov_expected(1:min(len_trim(provenance),PCG_RAW_PROV_LEN)) = &
                &provenance(1:min(len_trim(provenance),PCG_RAW_PROV_LEN))
        endif
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('add_raw_accum_weighted opening artifact', ierr)
        read(funit, iostat=ierr) magic, version
        call fileiochk('add_raw_accum_weighted reading magic', ierr)
        read(funit, iostat=ierr) state_file, eo_file, part_file, nparts_file, nptcls
        call fileiochk('add_raw_accum_weighted reading identity', ierr)
        read(funit, iostat=ierr) box_file, boxpd_file, padf_file, lims_file, smpd_file
        call fileiochk('add_raw_accum_weighted reading geometry', ierr)
        read(funit, iostat=ierr) prov_file
        call fileiochk('add_raw_accum_weighted reading provenance', ierr)
        if( magic /= PCG_RAW_ACCUM_MAGIC .or. version /= PCG_RAW_ACCUM_VERSION ) THROW_HARD('weighted raw PCG format mismatch')
        if( state_file /= state .or. eo_file /= eo .or. part_file /= part ) THROW_HARD('weighted raw PCG identity mismatch')
        if( nparts_file /= nparts .or. nptcls < 0 ) THROW_HARD('weighted raw PCG accumulator partition mismatch')
        if( box_file > self%box .or. padf_file /= self%padf .or. boxpd_file /= padf_file*box_file )then
            THROW_HARD('weighted raw box mismatch')
        endif
        fov_file = real(box_file) * smpd_file
        fov_self = real(self%box) * self%smpd
        if( abs(fov_file-fov_self) > 1.0e-4*fov_self ) THROW_HARD('weighted raw PCG field-of-view mismatch')
        if( any(lims_file(:,1) < self%lims3(:,1)) .or. any(lims_file(:,2) > self%lims3(:,2)) )then
            THROW_HARD('weighted raw PCG accumulator lattice is not nested in the current one')
        endif
        if( prov_file /= prov_expected ) THROW_HARD('weighted raw PCG accumulator provenance mismatch')
        if( box_file /= self%box .and. nptcls > 0 .and. weight > 0.0 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> PCG TRAIL: EMBEDDING PREVIOUS-CROP CHAIN (box ', &
                &box_file, ' -> ', self%box, ', CONSTANT-FOV ZERO-EXTENSION)'
        endif
        if( nptcls > 0 .and. weight > 0.0 )then
            allocate(bslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            allocate(dslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) bslice
                if( ierr /= 0 ) exit
                self%b_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) = &
                    &self%b_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) + weight * bslice
            enddo
            call fileiochk('add_raw_accum_weighted reading B', ierr)
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
                self%acc_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) = &
                    &self%acc_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) + weight * dslice
            enddo
            call fileiochk('add_raw_accum_weighted reading D', ierr)
            deallocate(bslice, dslice)
        else if( nptcls > 0 )then
            allocate(bslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            allocate(dslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) bslice
                if( ierr /= 0 ) exit
            enddo
            call fileiochk('add_raw_accum_weighted skipping B', ierr)
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
            enddo
            call fileiochk('add_raw_accum_weighted skipping D', ierr)
            deallocate(bslice, dslice)
        endif
        file_size  = -1_int64
        stream_pos = -1_int64
        inquire(unit=funit, size=file_size, pos=stream_pos)
        call fclose(funit)
        if( stream_pos-1_int64 /= file_size ) THROW_HARD('weighted raw PCG accumulator has trailing or missing bytes')
    end subroutine add_raw_accum_weighted

    !> Accumulate the weighted RHS B and Gram/sampling-density precursor D in
    !! one particle traversal. Both scatters use the same orientation, rotated
    !! coordinate and KB window, so evaluating that geometry twice is pure
    !! overhead. Keeping both updates in the same coloured OpenMP region also
    !! removes the per-particle parallel-region startup from accumulate_rhs.
    !! The update order within each accumulator is unchanged from the two
    !! standalone routines, which remain available to the low-level tests.
    subroutine accumulate_rhs_density( self, y_batch, nb, ifrom, bacc, dacc )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                    intent(in)    :: nb, ifrom
        complex,                    intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                              &self%lims2(2,1):self%lims2(2,2), nb)
        complex,                    intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                           &self%lims3(2,1):self%lims3(2,2),&
                                                           &self%lims3(3,1):self%lims3(3,2))
        real,                       intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
                                                        &self%lims3(2,1):self%lims3(2,2),&
                                                        &self%lims3(3,1):self%lims3(3,2))
        complex, allocatable :: weighted(:,:)
        real,    allocatable :: absT2(:,:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3)
        integer :: ib, i, g, h, k, l, i0(3)
        if( nb < 1 ) return
        allocate(weighted(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        allocate(absT2(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        !$omp parallel default(shared) private(ib,i,g,h,k,l,loc,i0,w,rot) proc_bind(close)
        do ib = 1, nb
            i = ifrom + ib - 1
            call self%prepare_fused_planes(i, y_batch(:,:,ib), weighted, absT2)
            do g = 1, self%nsym
                rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                        do k = self%lims2(2,1), self%lims2(2,2)
                            if( h*h + k*k > self%sqlp ) cycle
                            if( absT2(h,k) == 0. .and. weighted(h,k) == cmplx(0.,0.) ) cycle
                            loc = real(self%padf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( win_wraps(self, i0) ) cycle
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_window_pair_nowrap(self, i0, w, self%padsc * weighted(h,k), &
                                &absT2(h,k), bacc, dacc)
                        end do
                    end do
                    !$omp end do
                end do
                !$omp single
                do h = self%lims2(1,1), self%lims2(1,2)
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp   ) cycle
                        if( h*h + k*k <= self%sq_rim ) cycle
                        if( absT2(h,k) == 0. .and. weighted(h,k) == cmplx(0.,0.) ) cycle
                        loc = real(self%padf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. win_wraps(self, i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_window_pair(self, i0, w, self%padsc * weighted(h,k), &
                            &absT2(h,k), bacc, dacc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        deallocate(weighted, absT2)
    end subroutine accumulate_rhs_density

    !>  \brief  Closes accumulation: folds the RHS once, derives preconditioner
    !!          and (optionally) the Gram kernel, and releases the accumulators.
    subroutine end_accum( self, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                    intent(in)    :: l_kernel
        complex, allocatable :: bwork(:,:,:)
        real,    allocatable :: dwork(:,:,:)
        integer(timer_int_kind) :: tp
        if( .not. self%l_accum ) THROW_HARD('begin_accum has not been called; end_accum')
        if( self%reduction_nparts > 0 )then
            if( self%reduction_next_part /= self%reduction_nparts + 1 )then
                THROW_HARD('raw PCG reduction is incomplete')
            endif
        endif
        call self%reset_finalize_profile
        ! RHS first, and free its accumulator before the kernel work starts:
        ! finalize_khat is the allocation-heavy half and should not overlap it.
        ! move_alloc for the same aliasing reason as accumulate_batch.
        tp = tic()
        call move_alloc(self%b_work, bwork)
        self%b_rhs = self%fold_and_ifft(bwork)
        deallocate(bwork)
        call self%deapod_mul(self%b_rhs)
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        self%t_fin_rhs = real(toc(tp),dp)
        call move_alloc(self%acc_work, dwork)
        call self%finalize_density_accum(dwork, .true., l_kernel)
        deallocate(dwork)
        if( l_kernel ) call self%finalize_khat
        self%l_accum = .false.
    end subroutine end_accum

    ! SETUP: PRECONDITIONER AND KERNEL

    !>  \brief  sampling-density preconditioner M(k) = rho(k) + floor(sh), where
    !!          rho = sum_i G_i^dagger |T_i|^2 is exactly the gridding sampling
    !!          density -- the scatter without the gather. Without it the solver
    !!          would run with M = I, i.e. plain CG despite the name. With
    !!          heterogeneous CTFs H is severely ill-conditioned at CTF zeros,
    !!          so this is the other half of the performance problem: it cuts
    !!          iteration count, independently of per-iteration cost.
    !!
    !!          THE FLOOR MUST BE RELATIVE TO rho, NOT AN ABSOLUTE CONSTANT.
    !!          Using 1/(rho + lambda) with the solver's own lambda looks
    !!          harmless and is not: rho spans many orders of magnitude across
    !!          the oversampled lattice (it falls off as 1/sh^2 and is genuinely
    !!          zero in the gaps between rotated planes), so wherever rho ~ 0 the
    !!          preconditioner returns 1/lambda while well-sampled voxels get
    !!          1/rho -- a relative amplification of the least-constrained modes
    !!          by six to nine orders of magnitude. PCG then spends every
    !!          iteration chasing directions the data barely determines: the
    !!          residual falls for a few iterations, turns around, and the map
    !!          fills with noise. Production never hits this because
    !!          reconstructor%sampl_dens_correct leaves a voxel UNCHANGED when
    !!          rho is below 1e-6 (image%div_cmat_at_1) instead of dividing, and
    !!          because gridding is a single pass with no feedback.
    !!
    !!          Flooring at a fixed fraction of the shell-mean rho caps the
    !!          within-shell dynamic range at 1/RHO_FLOOR_FRAC while staying
    !!          symmetric positive definite, and is invariant to the operator's
    !!          overall scale -- which matters because rho here is accumulated
    !!          without the padsc factors that apply_normal carries.
    subroutine build_precond( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_precond')
        call self%ensure_wimg
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        call self%accumulate_absT2(acc)
        call self%precond_from_accum(acc)
    end subroutine build_precond

    !>  \brief  The ONE pass over the particles that both the preconditioner and
    !!          the Gram kernel are derived from.
    !!
    !!          These were two separate passes until it became apparent they are
    !!          the identical computation: build_precond went through
    !!          scatter_plane and build_kernel had its own inlined copy, but both
    !!          evaluate the same |T_i|^2, at the same rotated coordinate
    !!          padf*R_i*[h,k,0], through the same KB stencil. The accumulators
    !!          differed only by the constant padsc that scatter_plane applies --
    !!          and BOTH consumers are invariant to overall scale: the
    !!          preconditioner's floor is a fraction of the shell mean (see
    !!          precond_from_accum) and the kernel is least-squares rescaled by
    !!          calibrate_kernel. So the fusion is exact, not an approximation.
    !!
    !!          The particle loop is inside a single parallel region rather than
    !!          calling scatter_plane per particle, which opened and closed one
    !!          region per particle -- 5000 of them on a typical run.
    !!
    !!          ifrom/ito restrict the pass to a particle range, which is what
    !!          lets the caller stream batches instead of holding every plane
    !!          resident. Note this routine needs NO image data -- only the
    !!          scalars prep_particles cached -- so the range exists purely to
    !!          keep it in step with the RHS accumulation, which does.
    subroutine accumulate_absT2( self, acc, ifrom, ito )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(inout) :: acc(self%lims3(1,1):self%lims3(1,2),&
                                                        &self%lims3(2,1):self%lims3(2,2),&
                                                        &self%lims3(3,1):self%lims3(3,2))
        integer,          optional, intent(in)    :: ifrom, ito
        real,    allocatable :: absT2(:,:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3)
        integer :: i, g, h, k, l, i0(3), ii_from, ii_to
        ii_from = 1
        ii_to   = self%nptcls
        if( present(ifrom) ) ii_from = ifrom
        if( present(ito)   ) ii_to   = ito
        if( ii_to < ii_from ) return
        allocate(absT2(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        !$omp parallel default(shared) private(i,g,h,k,l,loc,i0,w,rot) proc_bind(close)
        do i = ii_from, ii_to
            call self%absT2_plane(i, absT2)
            ! coordinate replication: scatter |T_i|^2 at every symmetry-related
            ! orientation R_i.S_g. g is OUTSIDE the h-strided colour sweep so the
            ! scatter's separation guarantee still holds per orientation (2.7).
            ! symmats(:,:,1)=I, so nsym=1 reproduces the c1 pass bit-for-bit.
            do g = 1, self%nsym
                rot = matmul(self%rotmats(:,:,i), self%symmats(:,:,g))
                do l = 0, self%stride-1
                    !$omp do schedule(static,1)
                    do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                        do k = self%lims2(2,1), self%lims2(2,2)
                            if( h*h + k*k > self%sqlp ) cycle
                            if( absT2(h,k) == 0. ) cycle
                            loc = real(self%padf) * matmul(real([h,k,0]), rot)
                            i0  = nint(loc) - self%iwinsz
                            if( win_wraps(self, i0) ) cycle  ! rim: serial pass below
                            call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                            call scatter_window_nowrap(self, i0, w, absT2(h,k), acc)
                        end do
                    end do
                    !$omp end do
                end do
                ! wrapping rim, serialized: see win_wraps
                !$omp single
                do h = self%lims2(1,1), self%lims2(1,2)
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp   ) cycle
                        if( h*h + k*k <= self%sq_rim ) cycle   ! provably cannot wrap
                        if( absT2(h,k) == 0. ) cycle
                        loc = real(self%padf) * matmul(real([h,k,0]), rot)
                        i0  = nint(loc) - self%iwinsz
                        if( .not. win_wraps(self, i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        call scatter_window_real(self, i0, w, absT2(h,k), acc)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
    end subroutine accumulate_absT2

    subroutine precond_from_accum( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
                                                               &self%lims3(2,1):self%lims3(2,2),&
                                                               &self%lims3(3,1):self%lims3(3,2))
        call self%finalize_density_accum(rho_accum, .true., .false.)
    end subroutine precond_from_accum

    !> Set the absolute Tikhonov coefficient from a homogeneous scale of the
    !! raw data-only density. The first six native Fourier shells are common to
    !! full and cropped representations of the same physical box, and avoid
    !! making lambda depend on whether higher-resolution shells are present.
    !! Some Euclidean/simulated-data weights deliberately suppress that entire
    !! low-frequency band. In that case extend the same origin-centred support
    !! one native shell at a time until it contains data; do not confuse an
    !! intentionally empty low-frequency band with an empty accumulator.
    !!
    !! Zero bins remain in the selected support average: excluding them would
    !! make the scale nonlinear when fractional updates cover different voxels.
    !! Raw D is converted to the calibrated normal-operator convention by the
    !! same padsc**2 factor as Khat.
    subroutine update_lambda_from_density( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
                                                               &self%lims3(2,1):self%lims3(2,2),&
                                                               &self%lims3(3,1):self%lims3(3,2))
        integer, parameter :: DATA_SCALE_NATIVE_SHELLS = 6
        real(dp), allocatable :: shell_sum(:)
        integer,  allocatable :: shell_count(:)
        real(dp) :: dsum
        real     :: dval
        integer  :: h, k, m, hh, nscale, base_shell, scale_shell, shell, rsq
        allocate(shell_sum(0:self%Rnat), source=0.0_dp)
        allocate(shell_count(0:self%Rnat), source=0)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    rsq = h*h + k*k + m*m
                    if( rsq > (self%padf*self%Rnat)**2 ) cycle
                    shell = ceiling(sqrt(real(rsq)) / real(self%padf))
                    shell = min(shell, self%Rnat)
                    hh   = self%wrap(h)
                    dval = rho_accum(hh,k,m)
                    if( .not. ieee_is_finite(dval) .or. dval < 0.0 )then
                        THROW_HARD('invalid value in raw PCG density accumulator')
                    endif
                    shell_sum(shell)   = shell_sum(shell) + real(dval,dp)
                    shell_count(shell) = shell_count(shell) + 1
                end do
            end do
        end do
        base_shell  = min(DATA_SCALE_NATIVE_SHELLS, self%Rnat)
        scale_shell = base_shell
        dsum        = sum(shell_sum(0:scale_shell))
        nscale      = sum(shell_count(0:scale_shell))
        do while( dsum <= 0.0_dp .and. scale_shell < self%Rnat )
            scale_shell = scale_shell + 1
            dsum   = dsum   + shell_sum(scale_shell)
            nscale = nscale + shell_count(scale_shell)
        enddo
        deallocate(shell_sum, shell_count)
        if( nscale < 1 .or. dsum <= 0.0_dp ) THROW_HARD('cannot derive PCG data scale from empty D')
        if( scale_shell > base_shell )then
            write(logfhandle,'(A,I0)') &
                &'>>> PCG DATA SCALE: LOW BAND EMPTY; EXTENDED THROUGH NATIVE SHELL ', scale_shell
        endif
        self%data_scale = real(dsum / real(nscale,dp)) * self%padsc**2
        if( .not. ieee_is_finite(self%data_scale) .or. self%data_scale <= 0.0 )then
            THROW_HARD('invalid PCG data scale derived from D')
        endif
        if( self%l_lambda_relative ) self%lambda = self%lambda_rel * self%data_scale
        if( self%l_nu_prior        ) self%lambda_nu = self%nu_lambda_rel * self%data_scale
    end subroutine update_lambda_from_density

    !> Build P_tau from the independent-half FSC and raw data-only D.
    !! The padded Toeplitz lattice samples Fourier space padf times more finely
    !! than the native FSC, so each padded radius is assigned to its nearest
    !! native shell. The shell-mean D supplies 1/sigma2; multiplying the final
    !! diagonal by padsc**2 puts it in the calibrated Khat convention.
    subroutine build_ml_prior_from_density( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        real, intent(in) :: rho_accum(self%lims3(1,1):self%lims3(1,2), &
                                      &self%lims3(2,1):self%lims3(2,2), &
                                      &self%lims3(3,1):self%lims3(3,2))
        real(dp), allocatable :: shsum(:), shsum_thr(:,:)
        integer,  allocatable :: shcnt(:), shcnt_thr(:,:)
        integer :: cdim(3), h, k, m, hh, phys(3), shpd, sh, sz
        integer :: nthr, ithr, reslim_ind
        real    :: rval, cc, ssnr, prior_raw
        if( .not. self%l_ml_prior_requested ) return
        if( .not. allocated(self%ml_fsc) ) THROW_HARD('PCG ML FSC is not allocated')
        sz = min(size(self%ml_fsc), self%Rnat)
        if( sz < 1 ) THROW_HARD('PCG ML FSC has no usable shells')
        nthr = 1
        !$ nthr = omp_get_max_threads()
        allocate(shsum(0:sz), source=0.0_dp)
        allocate(shcnt(0:sz), source=0)
        allocate(shsum_thr(0:sz,nthr), source=0.0_dp)
        allocate(shcnt_thr(0:sz,nthr), source=0)
        !$omp parallel default(shared) private(h,k,m,hh,shpd,sh,rval,ithr)
        ithr = 1
        !$ ithr = omp_get_thread_num() + 1
        !$omp do collapse(2) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    shpd = nint(sqrt(real(h*h + k*k + m*m)))
                    sh   = nint(real(shpd) / real(self%padf))
                    if( sh < 0 .or. sh > sz ) cycle
                    hh   = self%wrap(h)
                    rval = rho_accum(hh,k,m)
                    if( rval <= 0.0 ) cycle
                    shsum_thr(sh,ithr) = shsum_thr(sh,ithr) + real(rval,dp)
                    shcnt_thr(sh,ithr) = shcnt_thr(sh,ithr) + 1
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        do ithr = 1, nthr
            do sh = 0, sz
                shsum(sh) = shsum(sh) + shsum_thr(sh,ithr)
                shcnt(sh) = shcnt(sh) + shcnt_thr(sh,ithr)
            enddo
        enddo
        deallocate(shsum_thr, shcnt_thr)
        cdim = self%wimg%get_array_shape()
        if( allocated(self%ml_prior) ) deallocate(self%ml_prior)
        allocate(self%ml_prior(cdim(1),cdim(2),cdim(3)), source=0.0)
        reslim_ind = max(6, calc_fourier_index(self%ml_hp, self%box, self%smpd))
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(h,k,m,phys,shpd,sh,cc,ssnr,prior_raw) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    shpd = nint(sqrt(real(h*h + k*k + m*m)))
                    sh   = nint(real(shpd) / real(self%padf))
                    if( sh < reslim_ind .or. sh > sz ) cycle
                    if( shcnt(sh) < 1 .or. shsum(sh) <= 1.0e-10_dp ) cycle
                    cc        = min(0.999, max(0.001, self%ml_fsc(sh)))
                    ssnr      = cc / (1.0 - cc)
                    prior_raw = real(shsum(sh) / real(shcnt(sh),dp)) / (self%ml_tau * ssnr)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    self%ml_prior(phys(1),phys(2),phys(3)) = prior_raw * self%padsc**2
                enddo
            enddo
        enddo
        !$omp end parallel do
        deallocate(shsum, shcnt)
        self%l_ml_prior = maxval(self%ml_prior) > 0.0
        if( .not. self%l_ml_prior ) THROW_HARD('PCG ML prior contains no positive bins')
    end subroutine build_ml_prior_from_density

    !> Fold the real density accumulator once into the packed Fourier layout.
    !! The preconditioner and Khat are two views of the same D accumulator, so
    !! producing them in one sphere-limited parallel pass avoids a second full
    !! padded-lattice traversal. Shell sums use one private column per OpenMP
    !! thread and are merged in thread-number order to keep the result stable.
    subroutine finalize_density_accum( self, rho_accum, l_precond, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
                                                               &self%lims3(2,1):self%lims3(2,2),&
                                                               &self%lims3(3,1):self%lims3(3,2))
        logical,                    intent(in)    :: l_precond, l_kernel
        real, parameter :: RHO_FLOOR_FRAC = 1.0e-2
        real(dp), allocatable :: shsum(:), shsum_thr(:,:)
        integer,  allocatable :: shcnt(:), shcnt_thr(:,:)
        real,     allocatable :: shfloor(:)
        integer :: h, k, cdim(3), m, hh, phys(3), sh, rho_lim
        integer :: khat_lim, work_lim, nthr, ithr
        real    :: denom, rval
        integer(timer_int_kind) :: tp
        if( .not. l_precond .and. .not. l_kernel ) return
        call self%update_lambda_from_density(rho_accum)
        cdim = self%wimg%get_array_shape()
        call self%build_ml_prior_from_density(rho_accum)
        ! Data only ever reaches |loc| <= padf*Rnat, so beyond that radius rho is
        ! identically zero: those modes are completely unconstrained. Zero them,
        ! exactly as reconstructor%sampl_dens_correct does for sh > sh_lim. A
        ! singular (PSD) preconditioner is the right tool here -- it keeps the
        ! Krylov space out of the null space entirely, so starting from x = 0
        ! those modes are never excited at all.
        rho_lim = self%padf * self%Rnat
        ! D can extend beyond the source sphere only by nearest-grid rounding
        ! and the finite KB stencil. This conservative bound is exact: points
        ! outside it cannot have received an accumulation update.
        khat_lim = rho_lim + ceiling(sqrt(3.0) * (real(self%iwinsz) + 0.5))
        work_lim = rho_lim
        if( l_kernel ) work_lim = max(work_lim, khat_lim)
        if( l_precond )then
            if( allocated(self%precond) ) deallocate(self%precond)
            allocate(self%precond(cdim(1),cdim(2),cdim(3)), source=0.0)
        endif
        if( l_kernel )then
            if( allocated(self%Khat) ) deallocate(self%Khat)
            allocate(self%Khat(cdim(1),cdim(2),cdim(3)), source=0.0)
        endif
        ! PASS 1: mean rho per shell, over the sampled voxels only. Averaging in
        ! the zeros would drag the floor down exactly where sampling is sparse,
        ! which is the regime the floor exists to protect.
        if( l_precond )then
            tp = tic()
            nthr = 1
            !$ nthr = omp_get_max_threads()
            allocate(shsum(0:rho_lim), source=0.0_dp)
            allocate(shcnt(0:rho_lim), source=0)
            allocate(shsum_thr(0:rho_lim,nthr), source=0.0_dp)
            allocate(shcnt_thr(0:rho_lim,nthr), source=0)
            !$omp parallel default(shared) private(h,k,m,hh,sh,rval,ithr)
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            !$omp do collapse(2) schedule(static)
            do m = self%lims3(3,1), self%lims3(3,2)
                do k = self%lims3(2,1), self%lims3(2,2)
                    do h = 0, self%lims3(1,2)
                        sh = nint(sqrt(real(h*h + k*k + m*m)))
                        if( sh > rho_lim ) cycle
                        hh   = self%wrap(h)
                        rval = rho_accum(hh,k,m)
                        if( rval <= 0.0 ) cycle
                        shsum_thr(sh,ithr) = shsum_thr(sh,ithr) + real(rval,dp)
                        shcnt_thr(sh,ithr) = shcnt_thr(sh,ithr) + 1
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel
            do ithr = 1, nthr
                do sh = 0, rho_lim
                    shsum(sh) = shsum(sh) + shsum_thr(sh,ithr)
                    shcnt(sh) = shcnt(sh) + shcnt_thr(sh,ithr)
                end do
            end do
            deallocate(shsum_thr, shcnt_thr)
            allocate(shfloor(0:rho_lim), source=0.0)
            do sh = 0, rho_lim
                if( shcnt(sh) > 0 )then
                    shfloor(sh) = RHO_FLOOR_FRAC * real(shsum(sh) / real(shcnt(sh),dp))
                endif
            end do
            deallocate(shsum, shcnt)
            self%t_fin_rho = real(toc(tp),dp)
        endif
        ! PASS 2: guarded reciprocal and packed Khat. An empty shell leaves the
        ! reciprocal at zero, i.e. treats those modes as unconstrained.
        tp = tic()
        !$omp parallel do collapse(2) default(shared) private(h,k,m,hh,phys,sh,denom) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > work_lim ) cycle
                    hh    = self%wrap(h)
                    phys  = self%wimg%comp_addr_phys(h,k,m)
                    if( l_kernel .and. sh <= khat_lim )then
                        self%Khat(phys(1),phys(2),phys(3)) = rho_accum(hh,k,m)
                    endif
                    if( l_precond )then
                        if( sh <= rho_lim )then
                            denom = max(rho_accum(hh,k,m), 0.0) + shfloor(sh)
                            if( self%l_ml_prior )then
                                denom = denom + self%ml_prior(phys(1),phys(2),phys(3)) / self%padsc**2
                            endif
                            if( self%l_nu_prior )then
                                denom = denom + self%nu_precond_shell_diag(sh)
                            endif
                            if( denom > 0.0 )then
                                self%precond(phys(1),phys(2),phys(3)) = 1.0 / denom
                            endif
                        endif
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        self%t_fin_fold = real(toc(tp),dp)
        if( l_precond )then
            self%l_precond = .true.
            deallocate(shfloor)
        endif
    end subroutine finalize_density_accum

    !>  \brief  Builds the Gram kernel.
    !!
    !!          The kernel is NOT the impulse response of the
    !!          matrix-free operator, h = H_data(delta_at_origin). That does not
    !!          work: the KB weights are normalized to sum 1 (see
    !!          kbinterpol%apod_mat_3d), so G applied to a constant Fourier
    !!          volume returns that constant, which makes the row sums of
    !!          M = sum_i G_i^dagger|T_i|^2 G_i exactly rho -- the gridding
    !!          density. A same-size kernel built that way reproduces gridding
    !!          exactly and PCG would converge in one step to the gridding
    !!          answer, gaining nothing.
    !!
    !!          The correct construction is the standard NUFFT Gram kernel:
    !!          scatter the weights |T_i(p)|^2 onto a 2x OVERSAMPLED Fourier
    !!          grid at DOUBLED coordinates. The extra resolution is precisely
    !!          what resolves sub-pixel frequency offsets, so the resulting
    !!          operator differs from gridding -- the oversampling is the point,
    !!          not an optional safety margin. Since
    !!              t = IFFT_2N(Khat)  and  H p = crop(IFFT_2N(FFT_2N(pad p) * Khat)),
    !!          the real-space kernel t is never needed explicitly: the
    !!          scattered array IS the multiplier.
    !!
    !!          Khat is real (real weights scattered with real KB weights) and
    !!          symmetric (lims2 is a full symmetric disk, so every plane point
    !!          appears with its negation), which is what makes the packed
    !!          half-grid representation exact.
    !!
    !!          Memory: the padded work image is (2*box)^3, about 1.1 GB at
    !!          box 256. That is the price of the kernelized operator.
    subroutine build_kernel( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_kernel')
        call self%ensure_wimg
        ! With oversampling in place the operator ALREADY works on a padf*box
        ! lattice, which is exactly the 2x grid the Gram kernel needs for a
        ! linear (non-wrapping) convolution of a box-supported volume. The two
        ! lattices coincide, so no separate padded image is required.
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        call self%accumulate_absT2(acc)
        call self%fold_accum_to_khat(acc)
        ! freed BEFORE finalize_khat, which is the memory-heavy half
        deallocate(acc)
        call self%finalize_khat
    end subroutine build_kernel

    !>  \brief  Builds preconditioner and kernel from a SINGLE particle pass.
    !!          This is what the solver path uses; build_precond and build_kernel
    !!          remain separately callable for the test commanders, which drive
    !!          one without the other.
    subroutine build_operators( self, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                    intent(in)    :: l_kernel
        real, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_operators')
        call self%ensure_wimg
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        call self%accumulate_absT2(acc)
        call self%finalize_density_accum(acc, .true., l_kernel)
        deallocate(acc)
        if( l_kernel ) call self%finalize_khat
    end subroutine build_operators

    subroutine fold_accum_to_khat( self, kacc )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: kacc(self%lims3(1,1):self%lims3(1,2),&
                                                          &self%lims3(2,1):self%lims3(2,2),&
                                                          &self%lims3(3,1):self%lims3(3,2))
        call self%finalize_density_accum(kacc, .false., .true.)
    end subroutine fold_accum_to_khat

    subroutine finalize_khat( self )
        class(reconstructor_pcg), intent(inout) :: self
        real,    parameter   :: EPS_D = 1.0e-8
        complex, allocatable :: ctmp(:,:,:)
        real,    allocatable :: tker(:,:,:), dep1d(:)
        real :: depval
        integer :: i, j, k, cdim(3)
        integer(timer_int_kind) :: tp
        cdim = self%wimg%get_array_shape()
        ! ---- divide out the DEPOSITION envelope ----
        ! Laying |T|^2 down through the KB window convolves the kernel spectrum
        ! with that window, which multiplies the real-space kernel by the
        ! window's transform. That is a SECOND envelope, distinct from the
        ! gather's, and left in it shows up as a spurious spatial taper. It is
        ! obtained from the exact discrete transform of the same normalized
        ! origin stencil used for deposition. Its separability eliminates a
        ! full complex accumulator and padded inverse FFT.
        tp = tic()
        call self%build_kb_envelope_1d(self%boxpd, dep1d)
        self%t_fin_dep = real(toc(tp),dp)
        tp = tic()
        ctmp = cmplx(self%Khat, 0.)
        call self%wimg%set_cmat(ctmp)
        call self%wimg%ifft()
        tker = self%wimg%get_rmat()
        !$omp parallel do collapse(3) default(shared) private(i,j,k,depval) schedule(static)
        do k = 1, self%boxpd
            do j = 1, self%boxpd
                do i = 1, self%boxpd
                    depval = dep1d(i) * dep1d(j) * dep1d(k)
                    if( abs(depval) > EPS_D )then
                        tker(i,j,k) = tker(i,j,k) / depval
                    else
                        tker(i,j,k) = 0.0
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        call self%wimg%set_rmat(tker, .false.)
        call self%wimg%fft()
        ctmp      = self%wimg%get_cmat()
        self%Khat = real(ctmp)
        deallocate(ctmp, tker, dep1d)
        self%l_kernel = .true.
        call self%calibrate_kernel
        self%t_fin_kernel = real(toc(tp),dp)
    end subroutine finalize_khat

    !>  \brief  One-off least-squares calibration of the kernel's overall scale
    !!          against the matrix-free reference.
    !!
    !!          The scale depends on the DFT normalization convention and on the
    !!          KB deposition used to build Khat. Rather than hard-code a
    !!          convention-dependent constant that would silently rot if the
    !!          image class changed its FFT scaling, it is measured once on a
    !!          deterministic probe.
    !!
    !!          This does NOT paper over a wrong operator: fitting a single
    !!          scalar removes only the scale degree of freedom, so whatever
    !!          error the equivalence test reports afterwards is pure SHAPE
    !!          mismatch -- which is exactly the quantity of interest. The
    !!          fitted factor is reported precisely so a value far from unity
    !!          shows up as the convention problem it would be.
    !>  \brief  Scales Khat to the matrix-free operator. The factor is
    !!          ANALYTIC, not fitted.
    !!
    !!          It used to be obtained by least squares against
    !!          apply_normal_matrixfree on a probe volume -- a full particle
    !!          pass, about 10 s of a 27 s setup, to determine one scalar. The
    !!          fitted value was 6.43e1 on the synthetic fixture and 6.398380e1
    !!          on a 5000-particle real data set: two independent data sets
    !!          within 0.5% of padsc**2 = (padf**3)**2 = 64. The fit was
    !!          measuring a constant.
    !!
    !!          Removing it is not merely an optimization. The trailing/
    !!          fractional-update scheme derives Khat from a stored accumulator
    !!          with no particles resident, so it CANNOT fit anything; an
    !!          analytic factor is the only form that survives that path. See
    !!          doc/policies/reconstruct3D_pcg_policy.md section 5.
    !!
    !!          measure_kernel_scale below still performs the fit, so the test
    !!          suite can assert the constant has not drifted.
    subroutine calibrate_kernel( self )
        class(reconstructor_pcg), intent(inout) :: self
        self%Khat = self%Khat * self%padsc**2
    end subroutine calibrate_kernel

    !>  \brief  Least-squares scale of the kernelized operator against the
    !!          matrix-free reference, on a deterministic probe. Returns 1.0 when
    !!          the analytic factor in calibrate_kernel is correct, so a test can
    !!          assert on the deviation directly. Costs one matrix-free particle
    !!          pass, which is why it is NOT on the setup path.
    function measure_kernel_scale( self ) result( scale )
        class(reconstructor_pcg), intent(inout) :: self
        real,    allocatable :: probe(:,:,:), hm(:,:,:), hk(:,:,:)
        real(dp) :: num, den
        real     :: lam_save, ctr, sig, dx, dy, dz, scale
        integer  :: i, j, k
        logical  :: l_ml_save, l_nu_save
        if( .not. self%l_kernel ) THROW_HARD('build_kernel has not been called; measure_kernel_scale')
        lam_save    = self%lambda
        l_ml_save   = self%l_ml_prior
        l_nu_save   = self%l_nu_prior
        self%lambda = 0.0   ! compare the DATA term only
        self%l_ml_prior = .false.
        self%l_nu_prior = .false.
        allocate(probe(self%box,self%box,self%box))
        ctr = real(self%box)/2.0 + 0.5
        sig = 0.15 * real(self%box)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    dx = real(i)-ctr; dy = real(j)-ctr; dz = real(k)-ctr
                    probe(i,j,k) = exp(-(dx*dx+dy*dy+dz*dz)/(2.0*sig*sig))
                end do
            end do
        end do
        hm  = self%apply_normal_matrixfree(probe)
        hk  = self%apply_normal_kernel(probe)
        num = sum(real(hm,dp)*real(hk,dp))
        den = sum(real(hk,dp)*real(hk,dp))
        scale = 1.0
        if( den > 0.0_dp ) scale = real(num/den)
        self%lambda = lam_save
        self%l_ml_prior = l_ml_save
        self%l_nu_prior = l_nu_save
    end function measure_kernel_scale

    ! PRIVATE HELPERS

    !>  \brief  |T_i|^2 = |C_i|^2 / sigma2_i over the lims2 disk.
    !!
    !!          apply_normal computes conjg(T)*(T*plane) = |T|^2 * plane, and
    !!          the shift factor exp(-i*2*pi*f.t) has unit modulus, so the SHIFT
    !!          CANCELS EXACTLY in the normal operator. H therefore depends only on this real, iteration-
    !!          invariant quantity: no shift phase, no complex multiply, and
    !!          half the transcendentals of a full complex path. The full
    !!          complex transfer is still needed for the RHS -- see
    !!          transfer_plane_cmplx.
    !!          The worksharing loops here are ORPHANED: when called from inside
    !!          apply_normal_matrixfree's or build_kernel's parallel region they
    !!          bind to it and the CTF evaluation is spread across the team;
    !!          when called from build_precond outside any parallel region they
    !!          simply run serially. Both are correct. Wrapping this in
    !!          !$omp single instead would serialize the CTF work and cap the
    !!          achievable speedup by Amdahl's law. Every element is written
    !!          inside a worksharing loop, so there is no benign-looking
    !!          all-threads-write race on the shared buffer.
    subroutine absT2_plane( self, iptcl, absT2 )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                    intent(in)  :: iptcl
        real,                       intent(out) :: absT2(self%lims2(1,1):self%lims2(1,2),&
                                                          &self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real      :: cval, sum_df, diff_df, angast, wl, half_wl2_cs, accc, phc, cterm, df, phsh, s2
        integer   :: h, k, shell
        logical   :: l_ctf, l_flip
        if( .not. self%l_use_ctf )then
            !$omp do collapse(2) schedule(static)
            do k = self%lims2(2,1), self%lims2(2,2)
                do h = self%lims2(1,1), self%lims2(1,2)
                    absT2(h,k) = 1.0
                end do
            end do
            !$omp end do
            return
        endif
        ! same ctfflag semantics as build_transfer; note |CTF|^2 == CTF^2, so the
        ! FLIP case only matters for the complex T, not here -- but the NO case
        ! does, and applying a CTF the images never saw would be wrong.
        l_ctf  = self%ctfparms(iptcl)%ctfflag /= CTFFLAG_NO
        l_flip = self%ctfparms(iptcl)%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(self%ctfparms(iptcl)%smpd, self%ctfparms(iptcl)%kv, &
                &self%ctfparms(iptcl)%cs, self%ctfparms(iptcl)%fraca)
            call tfun%init(self%ctfparms(iptcl)%dfx, self%ctfparms(iptcl)%dfy, self%ctfparms(iptcl)%angast)
            ! hoist the per-particle CTF constants and evaluate with the flat,
            ! call-free transcendental form of simple_math_ctf::ft_map_ctf_kernel.
            ! Its memoized (h,k) maps cover only the h>=0 half, so the kernel is
            ! inlined here to run over this full both-sign-h disk via the LUTs.
            ctfvals     = tfun%get_ctfvars(self%ctfparms(iptcl)%phshift)
            wl          = ctfvals%wl
            half_wl2_cs = 0.5 * wl * wl * ctfvals%cs
            sum_df      = ctfvals%dfx + ctfvals%dfy
            diff_df     = ctfvals%dfx - ctfvals%dfy
            angast      = ctfvals%angast
            accc        = ctfvals%amp_contr_const
            phc         = ctfvals%phshift
        endif
        !$omp do collapse(2) schedule(static) private(cval,shell,cterm,df,phsh,s2)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp )then
                    absT2(h,k) = 0.
                else
                    cval = 1.0
                    if( l_ctf )then
                        s2    = self%spafreqsq_lut(h,k)
                        cterm = cos( 2.0 * (self%ang_lut(h,k) - angast) )
                        df    = 0.5 * ( sum_df + cterm * diff_df )
                        phsh  = PI * wl * s2 * (df - half_wl2_cs * s2)
                        cval  = sin( phsh + phc + accc )
                        if( l_flip ) cval = abs(cval)
                    endif
                    shell      = self%shell_lut(h,k)
                    absT2(h,k) = cval * cval / self%sig2(shell,iptcl)
                endif
            end do
        end do
        !$omp end do
    end subroutine absT2_plane

    !>  \brief  full complex T_i for the cached particle iptcl (RHS only).
    subroutine transfer_plane_cmplx( self, iptcl, T )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                    intent(in)  :: iptcl
        complex,                    intent(out) :: T(self%lims2(1,1):self%lims2(1,2),&
                                                      &self%lims2(2,1):self%lims2(2,2))
        T = self%build_transfer(self%ctfparms(iptcl), self%shifts(:,iptcl), self%sig2(:,iptcl))
    end subroutine transfer_plane_cmplx

    !> Build the two plane values consumed by fused streaming accumulation.
    !! Every thread in the persistent region calls this routine; the orphaned
    !! omp-do partitions the plane between them. The expressions match
    !! build_transfer and absT2_plane, but the CTF is evaluated only once.
    subroutine prepare_fused_planes( self, iptcl, y_plane, weighted, absT2 )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                    intent(in)  :: iptcl
        complex,                    intent(in)  :: y_plane(self%lims2(1,1):self%lims2(1,2),&
                                                            &self%lims2(2,1):self%lims2(2,2))
        complex,                    intent(out) :: weighted(self%lims2(1,1):self%lims2(1,2),&
                                                             &self%lims2(2,1):self%lims2(2,2))
        real,                       intent(out) :: absT2(self%lims2(1,1):self%lims2(1,2),&
                                                         &self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        complex :: tval
        real    :: cval, arg, sw, sum_df, diff_df, angast, wl, half_wl2_cs
        real    :: accc, phc, cterm, df, phsh, s2
        integer :: h, k, shell
        logical :: l_ctf, l_flip
        l_ctf  = self%ctfparms(iptcl)%ctfflag /= CTFFLAG_NO
        l_flip = self%ctfparms(iptcl)%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(self%ctfparms(iptcl)%smpd, self%ctfparms(iptcl)%kv, &
                &self%ctfparms(iptcl)%cs, self%ctfparms(iptcl)%fraca)
            call tfun%init(self%ctfparms(iptcl)%dfx, self%ctfparms(iptcl)%dfy, &
                &self%ctfparms(iptcl)%angast)
            ctfvals     = tfun%get_ctfvars(self%ctfparms(iptcl)%phshift)
            wl          = ctfvals%wl
            half_wl2_cs = 0.5 * wl * wl * ctfvals%cs
            sum_df      = ctfvals%dfx + ctfvals%dfy
            diff_df     = ctfvals%dfx - ctfvals%dfy
            angast      = ctfvals%angast
            accc        = ctfvals%amp_contr_const
            phc         = ctfvals%phshift
        endif
        !$omp do collapse(2) schedule(static) &
        !$omp private(h,k,cval,arg,shell,sw,cterm,df,phsh,s2,tval)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp )then
                    weighted(h,k) = cmplx(0.,0.)
                    absT2(h,k)    = 0.
                    cycle
                endif
                if( .not. self%l_use_ctf )then
                    weighted(h,k) = y_plane(h,k)
                    absT2(h,k)    = 1.0
                    cycle
                endif
                cval = 1.0
                if( l_ctf )then
                    s2    = self%spafreqsq_lut(h,k)
                    cterm = cos(2.0 * (self%ang_lut(h,k) - angast))
                    df    = 0.5 * (sum_df + cterm * diff_df)
                    phsh  = PI * wl * s2 * (df - half_wl2_cs * s2)
                    cval  = sin(phsh + phc + accc)
                    if( l_flip ) cval = abs(cval)
                endif
                shell = self%shell_lut(h,k)
                arg   = 2.0 * PI * (real(h) * self%shifts(1,iptcl) + &
                    &real(k) * self%shifts(2,iptcl)) / real(self%box)
                sw    = 1.0 / sqrt(self%sig2(shell,iptcl))
                tval  = cval * cmplx(cos(arg), sin(arg)) * sw
                weighted(h,k) = conjg(tval) * y_plane(h,k) * sw
                absT2(h,k)    = cval * cval / self%sig2(shell,iptcl)
            end do
        end do
        !$omp end do
    end subroutine prepare_fused_planes

    !>  \brief  folds a full-range (both-sign h) complex volume accumulator into
    !!          the image's packed storage (h>=0 half) and inverse-FFTs it.
    !!
    !!          h = lims3(1,2) (the redundant Nyquist mate, present only because
    !!          lims3 spans the "including redundant Friedel mates" range) is
    !!          special: the scatter wraps via wlims and so never produces
    !!          +lims3(1,2), only its canonical negative representative. The
    !!          Nyquist bin's accumulated value lives at h = -lims3(1,2); the
    !!          wrap table sends it there.
    function fold_and_ifft( self, vol_accum ) result( z )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                &self%lims3(2,1):self%lims3(2,2),&
                                                                &self%lims3(3,1):self%lims3(3,2))
        real, allocatable :: z(:,:,:)
        integer :: h, hh, k, m, phys(3)
        integer(timer_int_kind) :: tp
        call self%ensure_wimg
        if( self%l_profile ) tp = tic()
        call self%wimg%zero_and_flag_ft()
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,phys) &
        !$omp schedule(static) proc_bind(close)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    call self%wimg%set_cmat_at(phys(1),phys(2),phys(3), vol_accum(hh,k,m))
                end do
            end do
        end do
        !$omp end parallel do
        call self%wimg%ifft()
        ! back to the native lattice: the unknown never leaves box^3
        z = self%crop_vol(self%wimg%get_rmat())
        if( self%l_profile ) self%t_fold = self%t_fold + real(toc(tp),dp)
    end function fold_and_ifft

    !>  \brief  z = M^-1 r via FFT, guarded elementwise divide, inverse FFT.
    !!
    !!          THE ENVELOPE BELONGS IN THE PRECONDITIONER TOO. rho is the
    !!          Fourier diagonal of the BARE operator T = A~^dagger W A~, but the
    !!          operator actually being solved is H = E^-1 T E^-1 (deapodization
    !!          on both sides, see deapod_mul). Inverting that gives
    !!          H^-1 = E T^-1 E, so the same real-space envelope has to bracket
    !!          the Fourier-domain divide -- multiplying by env, not invenv,
    !!          because this is the INVERSE of the deapodization sandwich.
    !!          Leaving it out means preconditioning with something that differs
    !!          from the true diagonal by a factor of E^2, which varies smoothly
    !!          but substantially across the box; that mismatch is the most
    !!          likely source of the non-monotone bump seen at iteration 2 with
    !!          both operators.
    function apply_precond( self, r ) result( z )
        class(reconstructor_pcg), intent(inout) :: self
        real,                       intent(in)    :: r(self%box,self%box,self%box)
        real,    allocatable :: z(:,:,:), rw(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        integer :: cdim(3), i, j, k
        integer(timer_int_kind) :: tp
        if( .not. self%l_precond )then
            allocate(z(self%box,self%box,self%box), source=r)
            call self%mask_mul(z)
            return
        endif
        call self%ensure_wimg
        if( self%l_profile ) tp = tic()
        allocate(rw(self%box,self%box,self%box), source=r)
        if( self%l_deapod ) rw = rw * self%env
        call self%wimg%set_rmat(self%pad_vol(rw), .false.)
        call self%wimg%fft()
        cmat = self%wimg%get_cmat()
        cdim = self%wimg%get_array_shape()
        !$omp parallel do collapse(3) default(shared) private(i,j,k) &
        !$omp schedule(static) proc_bind(close)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    cmat(i,j,k) = cmat(i,j,k) * self%precond(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        call self%wimg%set_cmat(cmat)
        call self%wimg%ifft()
        z = self%crop_vol(self%wimg%get_rmat())
        if( self%l_deapod ) z = z * self%env
        ! keeping z inside the support keeps the whole Krylov space there; M^-1
        ! is a Fourier diagonal and would otherwise leak the search directions
        ! back out into the solvent that P H P has just been set up to ignore
        call self%mask_mul(z)
        if( self%l_profile ) self%t_prec = self%t_prec + real(toc(tp),dp)
    end function apply_precond

    ! MODULE-LEVEL KERNELS (not type-bound: kept free of polymorphic dispatch so
    ! gfortran can inline them into the hot loops)

    !>  \brief  KB-weighted gather of one Fourier component from an image's
    !!          packed cmat, with the physical addressing and Friedel folding
    !!          inlined. A direct path through class(image)%comp_addr_phys
    !!          and %get_fcomp here would be true indirect calls into a submodule, 27
    !!          times per plane pixel, which also blocked vectorization.
    pure complex function gather_window( self, cmat, i0, w ) result( comp )
        class(reconstructor_pcg), intent(in) :: self
        complex,                    intent(in) :: cmat(:,:,:)
        integer,                    intent(in) :: i0(3)
        real,                       intent(in) :: w(self%wdim,self%wdim,self%wdim)
        integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
        ny   = self%boxpd
        nz   = self%boxpd
        comp = cmplx(0.,0.)
        do dk = 1, self%wdim
            mm = self%wrap(i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = self%wrap(i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = self%wrap(i0(1)+di-1)
                    if( hh >= 0 )then
                        ph = hh + 1
                        pk = kk + 1; if( kk < 0 ) pk = pk + ny
                        pm = mm + 1; if( mm < 0 ) pm = pm + nz
                        comp = comp + w(di,dj,dk) * cmat(ph,pk,pm)
                    else
                        ph = -hh + 1
                        pk = -kk + 1; if( -kk < 0 ) pk = pk + ny
                        pm = -mm + 1; if( -mm < 0 ) pm = pm + nz
                        comp = comp + w(di,dj,dk) * conjg(cmat(ph,pk,pm))
                    endif
                end do
            end do
        end do
    end function gather_window

    !>  \brief  One packed/Friedel traversal for a KB value and three spatial
    !!          derivatives. Kept non-type-bound so the 27-tap loop can inline.
    pure subroutine gather_window_grad( self, i0, w, dw, value, dvalue_dloc )
        type(pcg_fourier_workspace), intent(in)  :: self
        integer,                     intent(in)  :: i0(3)
        real(sp),                    intent(in)  :: w(self%wdim,self%wdim,self%wdim)
        real(sp),                    intent(in)  :: dw(self%wdim,self%wdim,self%wdim,3)
        complex,                     intent(out) :: value, dvalue_dloc(3)
        complex :: fcomp
        integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
        ny           = self%boxpd
        nz           = self%boxpd
        value        = cmplx(0.,0.)
        dvalue_dloc  = cmplx(0.,0.)
        do dk = 1, self%wdim
            mm = self%wrap(i0(3)+dk-1)
            do dj = 1, self%wdim
                kk = self%wrap(i0(2)+dj-1)
                do di = 1, self%wdim
                    hh = self%wrap(i0(1)+di-1)
                    if( hh >= 0 )then
                        ph = hh + 1
                        pk = kk + 1; if( kk < 0 ) pk = pk + ny
                        pm = mm + 1; if( mm < 0 ) pm = pm + nz
                        fcomp = self%cmat(ph,pk,pm)
                    else
                        ph = -hh + 1
                        pk = -kk + 1; if( -kk < 0 ) pk = pk + ny
                        pm = -mm + 1; if( -mm < 0 ) pm = pm + nz
                        fcomp = conjg(self%cmat(ph,pk,pm))
                    endif
                    ! V(loc) = sum_j w_j(loc) V_j.
                    value = value + w(di,dj,dk) * fcomp
                    ! dV/dloc_a = sum_j (dw_j/dloc_a) V_j.
                    dvalue_dloc = dvalue_dloc + dw(di,dj,dk,:) * fcomp
                end do
            end do
        end do
    end subroutine gather_window_grad

    !>  \brief  Does this window straddle the periodic wrap boundary?
    !!
    !!          THE COLOURING SCHEME IS ONLY VALID FOR WINDOWS THAT DO NOT WRAP.
    !!          The h-strided colouring guarantees that two h-lines of the same
    !!          colour map at least padf*stride apart, which exceeds the window
    !!          width -- but that separation is computed on UNWRAPPED
    !!          coordinates. The accumulator is periodic, so a coordinate just
    !!          past wlims(2) folds back to wlims(1), and two points that were
    !!          most of a period apart can land within a voxel of each other.
    !!          Two threads in the same colour sweep then write the same voxel.
    !!
    !!          Concretely at box 24: |loc| <= padf*R = 24 while
    !!          wlims = [-24,23], and h = -12 and h = +12 differ by 24 = 4*stride
    !!          so they share a colour. Their windows can overlap after folding.
    !!          The same holds at every box size -- it is confined to the Nyquist
    !!          rim, but it makes the result depend on thread scheduling, which
    !!          is how it was found (a batched and a monolithic accumulation of
    !!          identical data disagreed under threads and agreed on one).
    !!
    !!          Callers therefore split the plane: the interior scatters in the
    !!          parallel coloured pass, the wrapping rim in a serial pass.
    !!          Module-level and NOT type-bound, for the reason stated at the
    !!          head of this section: it is evaluated once per plane point per
    !!          particle -- of order 1e8 times per accumulation -- and a
    !!          type-bound call on a class(...) object goes through a dispatch
    !!          the compiler will not inline, which is exactly why
    !!          gather_window and scatter_window live here too.
    pure logical function win_wraps( self, i0 )
        class(reconstructor_pcg), intent(in) :: self
        integer,                   intent(in) :: i0(3)
        win_wraps = any(i0 < self%wlims(1)) .or. any(i0 + self%wdim - 1 > self%wlims(2))
    end function win_wraps

    !>  \brief  KB-weighted scatter of one value into the full-range accumulator.
    pure subroutine scatter_window( self, i0, w, val, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,                    intent(in)    :: val
        complex,                    intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
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

    !> Real-valued counterpart used by the sampling-density accumulator.
    pure subroutine scatter_window_real( self, i0, w, val, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim), val
        real,                       intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
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

    !> Update B and D through one wrapped KB-window traversal.
    pure subroutine scatter_window_pair( self, i0, w, bval, dval, bacc, dacc )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim), dval
        complex,                    intent(in)    :: bval
        complex,                    intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                           &self%lims3(2,1):self%lims3(2,2),&
                                                           &self%lims3(3,1):self%lims3(3,2))
        real,                       intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
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

    !>  \brief  Interior-only scatter for windows the caller has already proven
    !!          cannot wrap (win_wraps == .false.). There the period-box lookup
    !!          self%wrap is the identity over the whole window span, so indexing
    !!          vol_accum directly is bit-identical -- and it makes the inner run
    !!          contiguous (stride-1) and vectorizable, unlike the gathered
    !!          self%wrap form.
    pure subroutine scatter_window_nowrap( self, i0, w, val, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        real,                       intent(in)    :: val
        real,                       intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
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

    !> Update B and D through one non-wrapping, contiguous KB-window traversal.
    pure subroutine scatter_window_pair_nowrap( self, i0, w, bval, dval, bacc, dacc )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim), dval
        complex,                    intent(in)    :: bval
        complex,                    intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                           &self%lims3(2,1):self%lims3(2,2),&
                                                           &self%lims3(3,1):self%lims3(3,2))
        real,                       intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
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

    !>  \brief  Complex-valued counterpart of scatter_window_nowrap, for the RHS
    !!          scatter where the transfer sample is complex. Same interior-only
    !!          contract: the caller has ruled out wrapping, so direct indexing
    !!          is bit-identical to the self%wrap form and the inner run is
    !!          contiguous.
    pure subroutine scatter_window_cmplx_nowrap( self, i0, w, val, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        integer,                    intent(in)    :: i0(3)
        real,                       intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,                    intent(in)    :: val
        complex,                    intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
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

    !>  \brief  scatter a whole plane, h-strided so it is safe to call from
    !!          inside an OpenMP parallel region (used by the non-fused paths:
    !!          adjoint_plane_add, apply_adjoint_all, build_precond).
    subroutine scatter_plane( self, plane, rot, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        complex,                    intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                            &self%lims2(2,1):self%lims2(2,2))
        real,                       intent(in)    :: rot(3,3)
        complex,                    intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                &self%lims3(2,1):self%lims3(2,2),&
                                                                &self%lims3(3,1):self%lims3(3,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot_g(3,3)
        integer :: h, k, l, g, i0(3)
        ! coordinate replication over the point group: scatter the plane at every
        ! R_i.S_g. g is outside the h-strided colour sweep so the scatter stays
        ! race-free per orientation (2.7); symmats(:,:,1)=I gives the c1 pass.
        !$omp parallel default(shared) private(h,k,l,g,loc,i0,w,rot_g) proc_bind(close)
        do g = 1, self%nsym
            rot_g = matmul(rot, self%symmats(:,:,g))
            do l = 0, self%stride-1
                !$omp do schedule(static,1)
                do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp ) cycle
                        if( plane(h,k) == cmplx(0.,0.) ) cycle
                        loc = real(self%padf) * matmul(real([h,k,0]), rot_g)
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
        ! Serial pass over the wrapping rim, where the colouring's separation
        ! guarantee does not survive folding (see win_wraps). Confined to the
        ! outermost shell, so the cost is negligible, and being serial it is also
        ! reproducible -- which an atomic would not be, since the summation order
        ! would still vary.
        do g = 1, self%nsym
            rot_g = matmul(rot, self%symmats(:,:,g))
            do h = self%lims2(1,1), self%lims2(1,2)
                do k = self%lims2(2,1), self%lims2(2,2)
                    if( h*h + k*k > self%sqlp   ) cycle
                    if( h*h + k*k <= self%sq_rim ) cycle       ! provably cannot wrap
                    if( plane(h,k) == cmplx(0.,0.) ) cycle
                    loc = real(self%padf) * matmul(real([h,k,0]), rot_g)
                    i0  = nint(loc) - self%iwinsz
                    if( .not. win_wraps(self, i0) ) cycle
                    call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                    call scatter_window(self, i0, w, self%padsc * plane(h,k), vol_accum)
                end do
            end do
        end do
    end subroutine scatter_plane

    ! GETTERS

    pure function get_lims2( self ) result( lims2 )
        class(reconstructor_pcg), intent(in) :: self
        integer :: lims2(2,2)
        lims2 = self%lims2
    end function get_lims2

    pure function get_lims3( self ) result( lims3 )
        class(reconstructor_pcg), intent(in) :: self
        integer :: lims3(3,2)
        lims3 = self%lims3
    end function get_lims3

    pure integer function get_nptcls( self )
        class(reconstructor_pcg), intent(in) :: self
        get_nptcls = self%nptcls
    end function get_nptcls

    pure real function get_data_scale( self )
        class(reconstructor_pcg), intent(in) :: self
        get_data_scale = self%data_scale
    end function get_data_scale

    pure real function get_effective_lambda( self )
        class(reconstructor_pcg), intent(in) :: self
        get_effective_lambda = self%lambda
    end function get_effective_lambda

    ! SOLVER

    !>  \brief  preconditioned CG solve of H x = b. With
    !!          build_precond called this is genuinely preconditioned; without
    !!          it, M = I and this degenerates to plain CG.
    !!
    !!          This form takes every observed plane at once and is retained for
    !!          the test commanders, whose fixtures are small. Production-sized
    !!          callers stream batches through begin_accum/accumulate_batch/
    !!          end_accum and then call solve_accum, which never materializes
    !!          y_planes at all.
    subroutine solve( self, y_planes, x, maxits, rtol, rel_res_hist, niters, outcome )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                               &self%lims2(2,1):self%lims2(2,2), *)
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out) :: outcome
        if( allocated(self%b_rhs) ) deallocate(self%b_rhs)
        self%b_rhs = self%apply_adjoint_all(y_planes)
        ! b' = P b, completing the (P H P) u = P b normal equations
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters, outcome)
    end subroutine solve

    !>  \brief  Solve against the RHS built by end_accum. Same solver, no
    !!          observed planes resident.
    subroutine solve_accum( self, x, maxits, rtol, rel_res_hist, niters, outcome )
        class(reconstructor_pcg), intent(inout) :: self
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out) :: outcome
        if( .not. self%l_rhs ) THROW_HARD('end_accum has not been called; solve_accum')
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters, outcome)
    end subroutine solve_accum

    !>  \brief  The solver proper. Reads the RHS from self%b_rhs rather than
    !!          taking it as an argument: passing a component of `self` as a
    !!          separate dummy while `self` is intent(inout) is an aliasing
    !!          hazard the standard lets a compiler exploit, and the alternative
    !!          (copying it) would cost 67 MB per solve at box 256 for nothing.
    subroutine solve_core( self, x, maxits, rtol, rel_res_hist, niters, outcome )
        class(reconstructor_pcg), intent(inout) :: self
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out) :: outcome
        ! WHICH RESIDUAL IS REPORTED MATTERS. This used to report
        ! ||r||_Minv / ||r0||_Minv, the PRECONDITIONED norm, which in PCG is not
        ! monotone -- only ||e||_H is -- and with a singular preconditioner it
        ! wanders freely. It duly oscillated, which reads exactly like a solver
        ! that has lost conjugacy while the solve was in fact converging
        ! perfectly well. The headline number and the stopping test are now the
        ! true relative residual ||r||_2 / ||b||_2, which is what a reader
        ! actually wants and costs two dot products against 1.5 s of FFTs. The
        ! M-norm is retained in the solver outcome for the diagnostic file,
        ! because a large gap between it and the true residual says M is a poor
        ! model of H.
        !
        ! RESIDUAL REPLACEMENT. CG propagates the residual by the recurrence
        ! r <- r - alpha*Hp rather than recomputing b - Hx, which is what makes
        ! an iteration cost one operator application instead of two; the two can
        ! drift apart in finite precision. Measured here, they agree to six
        ! significant figures at box 256, so this is kept only as a periodic
        ! audit at a long interval rather than a correction -- reporting ||r||_2
        ! from the recurrence is only legitimate because this check backs it up.
        integer, parameter :: RESID_REPLACE = 25
        ! diminishing-returns stop: relative model-update tolerance dx/x. On
        ! noisy real data |r|/|b| plateaus above rtol while dx/x keeps falling,
        ! so this is what actually terminates a real solve. Internal default;
        ! raise to stop sooner, lower to iterate longer. Suppressed by rtol <= 0,
        ! which is the caller's way of saying "run exactly maxits iterations, no
        ! early exits" -- test=pcg_recon stage 7 depends on that, since comparing
        ! two solves stopped by a data-dependent criterion tests nothing.
        real, parameter :: PCG_XTOL = 1.5e-2
        real, allocatable :: r(:,:,:), p(:,:,:), hp(:,:,:), z(:,:,:), hist(:)
        real, allocatable :: update_hist(:), mnorm_hist(:), iteration_times(:)
        real, allocatable :: rtr(:,:,:)
        real(dp) :: rho, rho_new, rho0, alpha, beta, pHp
        real(dp) :: bnorm, rnorm, xnorm, dxnorm, mnorm, dxx
        integer  :: mmaxits, iter, n_done
        real     :: rrtol
        logical  :: stop_rtol, stop_xtol
        type(pcg_solver_outcome) :: result
        integer(timer_int_kind) :: t_it
        mmaxits = 50
        if( present(maxits) ) mmaxits = maxits
        if( mmaxits < 1 ) THROW_HARD('maxits must be at least 1; solve')
        rrtol = 1.0e-4
        if( present(rtol) ) rrtol = rtol
        if( .not. ieee_is_finite(rrtol) ) THROW_HARD('rtol must be finite; solve')
        result%requested_maxits = mmaxits
        if( rrtol <= 0.0 ) result%stop_reason = 'fixed_iterations'
        if( rrtol > 0.0 )  result%stop_reason = 'maxits'
        allocate(hist(mmaxits), update_hist(mmaxits), iteration_times(mmaxits))
        allocate(mnorm_hist(mmaxits), source=-1.0)
        ! profile the ITERATIONS only: forming the RHS is a one-off setup cost
        ! and folding it in would flatter whichever phase it happens to share.
        call self%reset_profile
        if( all(x == 0.0) )then
            ! zero initialization is the documented baseline;
            ! skip a full operator application that is known to return zero
            allocate(hp(self%box,self%box,self%box), source=0.0)
        else
            hp = self%apply_normal(x)
        endif
        r  = self%b_rhs - hp
        z  = self%apply_precond(r)
        p  = z
        rho  = self%dot_real_volume(r,z)
        rho0 = rho
        if( rho0 <= 0.0_dp ) THROW_HARD('non-positive initial dot(r,z); preconditioner is not positive definite; solve')
        bnorm = sqrt(self%dot_real_volume(self%b_rhs,self%b_rhs))
        if( bnorm <= 0.0_dp ) THROW_HARD('zero right-hand side; nothing to reconstruct; solve')
        rnorm = sqrt(self%dot_real_volume(r,r))
        result%initial_rel_residual = real(rnorm / bnorm)
        n_done = 0
        dxx = 0.0_dp
        do iter = 1, mmaxits
            t_it = tic()
            hp  = self%apply_normal(p)
            pHp = self%dot_real_volume(p,hp)
            if( pHp <= 0.0_dp ) THROW_HARD('non-positive dot(p,Hp); PCG lost positive-definiteness; solve')
            if( pHp /= pHp )    THROW_HARD('non-finite dot(p,Hp); solve')
            alpha = rho / pHp
            x  = x + real(alpha) * p
            r  = r - real(alpha) * hp
            if( mod(iter, RESID_REPLACE) == 0 )then
                rtr  = self%b_rhs - self%apply_normal(x)
                r = rtr
            endif
            n_done  = iter
            ! headline: true relative residual. dx/x says how much the map is
            ! still moving, which for a reconstruction is often the more
            ! practical stopping signal than any residual level.
            rnorm      = sqrt(self%dot_real_volume(r,r))
            xnorm      = sqrt(self%dot_real_volume(x,x))
            dxnorm     = abs(alpha) * sqrt(self%dot_real_volume(p,p))
            dxx        = dxnorm / max(xnorm, epsilon(1.0_dp))
            hist(iter) = real(rnorm / bnorm)
            update_hist(iter) = real(dxx)
            stop_rtol  = rrtol > 0.0 .and. rnorm / bnorm <= real(rrtol,dp)
            stop_xtol  = rrtol > 0.0 .and. dxx <= real(PCG_XTOL,dp)
            if( stop_rtol .or. stop_xtol .or. iter == mmaxits )then
                iteration_times(iter) = real(toc(t_it))
                if( stop_rtol )then
                    result%stop_reason = 'rtol'
                    result%converged   = .true.
                    exit
                endif
                if( stop_xtol )then
                    result%stop_reason = 'xtol'
                    result%converged   = .true.
                    exit
                endif
                exit
            endif
            z       = self%apply_precond(r)
            rho_new = self%dot_real_volume(r,z)
            mnorm   = sqrt(abs(rho_new)/rho0)
            mnorm_hist(iter) = real(mnorm)
            iteration_times(iter) = real(toc(t_it))
            beta = rho_new / rho
            p    = z + real(beta) * p
            rho  = rho_new
        end do
        ! x = P u: the CG variable u is unconstrained outside the support (P H P
        ! annihilates it there, so it never influenced the residual), but it does
        ! accumulate arbitrary values via the preconditioner. This is the step
        ! that makes the returned volume the constrained solution.
        call self%mask_mul(x)
        result%iteration_count  = n_done
        result%final_rel_update = real(dxx)
        if( n_done > 0 ) result%final_rel_residual = hist(n_done)
        if( n_done > 0 )then
            allocate(result%rel_residual_history(n_done), source=hist(1:n_done))
            allocate(result%rel_update_history(n_done), source=update_hist(1:n_done))
            allocate(result%preconditioned_residual_history(n_done), source=mnorm_hist(1:n_done))
            allocate(result%iteration_seconds(n_done), source=iteration_times(1:n_done))
        endif
        if( present(niters) ) niters = n_done
        if( present(rel_res_hist) ) allocate(rel_res_hist(n_done), source=hist(1:n_done))
        if( present(outcome) ) outcome = result
        self%l_profile = .false.
    end subroutine solve_core

    ! PROFILING

    subroutine reset_finalize_profile( self )
        class(reconstructor_pcg), intent(inout) :: self
        self%t_fin_rhs    = 0.0_dp
        self%t_fin_rho    = 0.0_dp
        self%t_fin_fold   = 0.0_dp
        self%t_fin_dep    = 0.0_dp
        self%t_fin_kernel = 0.0_dp
    end subroutine reset_finalize_profile

    subroutine report_finalize_profile( self, funit )
        class(reconstructor_pcg), intent(in) :: self
        integer, optional,          intent(in) :: funit
        real(dp) :: total
        integer  :: out_unit
        out_unit = logfhandle
        if( present(funit) ) out_unit = funit
        total = self%t_fin_rhs + self%t_fin_rho + self%t_fin_fold + &
            &self%t_fin_dep + self%t_fin_kernel
        write(out_unit,'(a)') '>>> PCG ACCUMULATOR FINALIZATION (seconds)'
        write(out_unit,'(a,f9.3)') '    RHS fold + deapod + support : ', self%t_fin_rhs
        write(out_unit,'(a,f9.3)') '    rho shell statistics        : ', self%t_fin_rho
        write(out_unit,'(a,f9.3)') '    fused reciprocal + Khat pack: ', self%t_fin_fold
        write(out_unit,'(a,f9.3)') '    deposition envelope setup   : ', self%t_fin_dep
        write(out_unit,'(a,f9.3)') '    kernel correction + FFT     : ', self%t_fin_kernel
        write(out_unit,'(a,f9.3)') '    ---- accounted subtotal     : ', total
    end subroutine report_finalize_profile

    subroutine reset_profile( self, l_on )
        class(reconstructor_pcg), intent(inout) :: self
        logical, optional,          intent(in)    :: l_on
        self%t_setvol = 0.0_dp
        self%t_cmatcp = 0.0_dp
        self%t_ploop  = 0.0_dp
        self%t_fold   = 0.0_dp
        self%t_khat   = 0.0_dp
        self%t_prec   = 0.0_dp
        self%l_profile = .true.
        if( present(l_on) ) self%l_profile = l_on
    end subroutine reset_profile

    !>  \brief  Per-iteration breakdown of where the solve's time actually goes.
    !!
    !!          Read it as one question: how much of an iteration is the PARTICLE
    !!          LOOP (t_ploop -- the only part the kernelized operator removes)
    !!          and how much is FFT plus bulk traffic on the padded lattice
    !!          (t_setvol + t_cmatcp + t_fold + t_prec -- which switching
    !!          operator does NOT remove, and which only a Fourier-domain
    !!          formulation of the solve eliminates)?
    subroutine report_profile( self, niters, funit )
        class(reconstructor_pcg), intent(in) :: self
        integer,                    intent(in) :: niters
        integer, optional,          intent(in) :: funit
        real(dp) :: rn, tot, ffts
        integer  :: out_unit
        if( niters < 1 ) return
        out_unit = logfhandle
        if( present(funit) ) out_unit = funit
        rn   = real(niters,dp)
        ffts = self%t_setvol + self%t_cmatcp + self%t_fold + self%t_prec
        tot  = ffts + self%t_ploop + self%t_khat
        write(out_unit,'(a)')    '>>> PCG PROFILE (seconds per iteration)'
        write(out_unit,'(a,f9.3)') '    pad + fwd FFT of iterate     : ', self%t_setvol / rn
        write(out_unit,'(a,f9.3)') '    get_cmat/set_cmat copies     : ', self%t_cmatcp / rn
        write(out_unit,'(a,f9.3)') '    particle loop                : ', self%t_ploop  / rn
        write(out_unit,'(a,f9.3)') '    kernel pointwise multiply    : ', self%t_khat   / rn
        write(out_unit,'(a,f9.3)') '    fold + inv FFT + crop        : ', self%t_fold   / rn
        write(out_unit,'(a,f9.3)') '    apply_precond (2 FFT + copy) : ', self%t_prec   / rn
        write(out_unit,'(a,f9.3)') '    ---- accounted subtotal      : ', tot / rn
        if( tot > 0.0_dp )then
            write(out_unit,'(a,f7.1,a)') '    particle loop is ', &
                &100.0_dp * self%t_ploop / tot, '% of accounted time'
            write(out_unit,'(a,f7.1,a)') '    FFT + lattice traffic is ', &
                &100.0_dp * ffts / tot, '% of accounted time'
        endif
    end subroutine report_profile

end module simple_reconstructor_pcg
