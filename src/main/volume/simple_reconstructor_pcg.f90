!@descr: CTF/sigma-weighted Fourier-projection operator and preconditioned
!  conjugate-gradient volume solver of the reconstruct3D PCG backend, see
!  doc/policies/reconstruct3D_pcg_policy.md. Per-particle data is cached once
!  (prep_particles), the particle loops are OpenMP-parallel, and the optional
!  kernelized (Toeplitz) normal operator makes the per-iteration cost independent
!  of the particle count; the matrix-free operator is the exact reference.
module simple_reconstructor_pcg
use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image, only: image
use simple_ctf,   only: ctf
use simple_cartesian_fourier, only: center_embed_real3d, center_crop_real3d, &
    &extract_native_fourier_plane, gather_packed_window
use simple_gridding, only: kb_stencil_envelope_1d, kb_stencil_centered_crop_inv_envelope_1d
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
implicit none

public :: reconstructor_pcg, pcg_solver_outcome
public :: PCG_OP_MATRIXFREE, PCG_OP_KERNEL
public :: pcg_raw_accum_compatible
private
#include "simple_local_flags.inc"

integer, parameter :: PCG_OP_MATRIXFREE = 0 !< reference operator: exact, cost ~ nptcls per iteration
integer, parameter :: PCG_OP_KERNEL     = 1 !< kernelized Toeplitz operator: cost independent of nptcls
integer, parameter :: PCG_RAW_ACCUM_VERSION = 1
integer, parameter :: PCG_RAW_PROV_LEN = 256
real,    parameter :: PCG_SUPPORT_DIV_MIN = 0.1 !< window floor of the x/window warm-start conversion, see window_div
logical, parameter :: PCG_HARD_SOLVE_SUPPORT = .true. !< solve on the hard domain window > 0 and window the output; .false. = soft P H P, see install_support
character(len=16), parameter :: PCG_RAW_ACCUM_MAGIC = 'SIMPLE_PCG_RAW01'

!> stop_reason when dot(p,Hp) is non-positive or non-finite: the iterate is
!! returned as it stood before that step and the caller decides (restart or fail)
character(len=*), parameter, public :: PCG_STOP_INDEFINITE = 'indefinite'

type :: pcg_solver_outcome
    character(len=24) :: stop_reason          = 'not_started'
    integer           :: iteration_count      = 0
    integer           :: requested_maxits     = 0
    real              :: initial_rel_residual = 0.0
    real              :: final_rel_residual   = 0.0
    real              :: final_rel_update     = 0.0
    real(dp)          :: failure_curvature    = 0.0_dp      !< curvature of this attempt's indefinite stop
    integer           :: failure_iteration    = 0           !< iteration of this attempt's indefinite stop
    real(dp)          :: restart_trigger_curvature = 0.0_dp !< warm-attempt curvature that caused a cold restart
    integer           :: restart_trigger_iteration = 0      !< warm-attempt iteration that caused a cold restart
    logical           :: cold_restart_used    = .false.     !< outcome is from the one permitted cold retry
    logical           :: start_rejected       = .false.     !< the nonzero start was worse than zero and was discarded before iterating
    real              :: rejected_start_initial = 0.0       !< initial relative residual of the discarded start
    logical           :: converged            = .false.
    real, allocatable  :: rel_residual_history(:)
    real, allocatable  :: rel_update_history(:)
    real, allocatable  :: preconditioned_residual_history(:)
    real, allocatable  :: iteration_seconds(:)
end type pcg_solver_outcome

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
    ! (h,k)-only lookup tables over the fixed lims2 disk, see build_hk_luts
    real,    allocatable :: spafreqsq_lut(:,:)       !< spatial frequency squared
    real,    allocatable :: ang_lut(:,:)             !< atan2(k,h) astigmatism angle
    integer, allocatable :: shell_lut(:,:)           !< resolution shell index, capped at Rnat
    real,    allocatable :: env(:,:,:)               !< measured KB instrument envelope, see build_env
    real,    allocatable :: invenv(:,:,:)            !< its guarded reciprocal, for deapodization
    logical              :: l_deapod = .true.        !< correct the KB roll-off, see deapod_mul
    real,    allocatable :: mask(:,:,:)              !< solve support P (hard: window > 0), see install_support
    real,    allocatable :: window(:,:,:)            !< soft output window; the shipped map is window*u, see set_mask
    logical              :: l_mask = .false.
    type(image)          :: wimg                     !< persistent box^3 work image (keeps its FFTW plans)
    logical              :: wimg_exists = .false.
    integer              :: fft_nthreads = 1         !< threads owned by each persistent FFTW plan
    ! ---- streaming accumulation over particle batches ----
    ! the only particle-dependent state the solve needs; precond, Khat and the RHS derive from them
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
    ! requested before D is finalized; the absolute scale is derived from D
    real, allocatable :: ml_fsc(:)                    !< independent-half FSC, shells 1:R
    real, allocatable :: ml_prior(:,:,:)              !< calibrated padded Fourier diagonal
    real              :: ml_tau = 1.0                 !< established ML regularization fudge factor
    real              :: ml_hp  = 100.0               !< low-frequency no-prior limit in Angstrom
    logical           :: l_ml_prior_requested = .false.
    logical           :: l_ml_prior = .false.
    ! ---- per-phase profiling over a solve: particle loop vs FFT + lattice traffic ----
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
    logical  :: exists = .false.
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
    ! STREAMING SETUP (batch at a time, see begin_accum)
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
    procedure, private :: finalize_density_accum
    procedure, private :: finalize_khat
    procedure :: set_deapod
    procedure :: set_mask
    procedure :: set_mask_volume
    procedure :: set_lambda_relative
    procedure :: set_ml_prior
    procedure, private :: build_env
    procedure, private :: build_hk_luts
    procedure, private :: deapod_mul
    procedure, private :: install_support
    procedure, private :: mask_mul
    procedure, private :: window_mul
    procedure, private :: window_div
    procedure, private :: calibrate_kernel
    procedure :: measure_kernel_scale
    procedure :: set_op_mode
    ! LOW-LEVEL OPERATOR (public: the test commanders verify the adjoint identity)
    procedure :: set_volume
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
    procedure, private :: update_lambda_from_density
    procedure, private :: build_ml_prior_from_density
    procedure, private :: apply_fourier_diagonal
end type reconstructor_pcg

contains

    !> wall-clock helpers without module-global state, safe for concurrent even/odd solves
    integer(timer_int_kind) function pcg_tic() result(tstart)
        call system_clock(count=tstart)
    end function pcg_tic

    real(dp) function pcg_toc( tstart ) result(seconds)
        integer(timer_int_kind), intent(in) :: tstart
        integer(timer_int_kind) :: tend, rate
        call system_clock(count=tend, count_rate=rate)
        seconds = real(tend-tstart,dp) / real(rate,dp)
    end function pcg_toc

    ! CONSTRUCTOR

    subroutine new( self, box, smpd, lambda, fft_nthreads )
        class(reconstructor_pcg),           intent(inout) :: self
        integer,                            intent(in)    :: box
        real,                               intent(in)    :: smpd
        real,                     optional, intent(in)    :: lambda
        integer,                  optional, intent(in)    :: fft_nthreads
        type(image) :: tmp
        integer     :: R, lo, hi, i
        real        :: rlim
        call self%kill
        self%box    = box
        self%smpd   = smpd
        self%lambda = 0.0
        if( present(lambda) ) self%lambda = lambda
        self%fft_nthreads = nthr_glob
        if( present(fft_nthreads) ) self%fft_nthreads = max(1, fft_nthreads)
        self%kbwin  = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        ! odd width centred on nint(loc): negating loc negates the window as a set,
        ! which the h>=0 fold in fold_and_ifft relies on
        self%wdim   = 2*self%iwinsz + 1
        ! same-colour h-lines stay a full window apart along one axis after rotation
        ! (max-norm >= Euclidean norm / sqrt(3)), as in reconstructor%insert_plane_oversamp
        self%stride = ceiling(sqrt(3.0) * real(self%wdim))
        ! oversampling: the unknown lives on the native box, every Fourier operation
        ! on the padf-times padded lattice (centre-pad in, centre-crop out), as in
        ! Fourier gridding; native-lattice KB interpolation is only percent-accurate
        self%padf    = OSMPL_PAD_FAC
        self%boxpd   = self%padf * box
        self%pad_off = (self%boxpd - box) / 2
        self%Rnat    = box / 2
        ! fwd_ft divides by product(ldim): restores the native scale at coincident frequencies
        self%padsc   = real(self%padf)**3
        call tmp%new([self%boxpd,self%boxpd,self%boxpd], smpd, &
            &wthreads=self%fft_nthreads > 1, fft_nthreads=self%fft_nthreads)
        call tmp%fft()
        self%lims3 = tmp%loop_lims(3)
        ! true period-box wrap range: lims3(1,:) spans both Friedel Nyquist mates
        ! (one longer than the period), axes 2/3 do not
        self%wlims = self%lims3(2,:)
        ! full symmetric (both-sign h) NATIVE disk: forward_plane/adjoint_plane_add
        ! are then an exact adjoint pair for any orientation; padded loc = padf*loc
        R = self%Rnat
        self%lims2(1,:) = [-R, R]
        self%lims2(2,:) = [-R, R]
        self%sqlp       = R*R
        ! squared plane radius below which a KB window cannot reach the wrap boundary;
        ! |loc| = padf*sqrt(h^2+k^2) is rotation-independent, so (h,k) decides (conservative)
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
        ! c1 default: a single identity operator, replaced by set_sym
        self%nsym = 1
        if( allocated(self%symmats) ) deallocate(self%symmats)
        allocate(self%symmats(3,3,1), source=0.0)
        self%symmats(1,1,1) = 1.0; self%symmats(2,2,1) = 1.0; self%symmats(3,3,1) = 1.0
        self%exists = .true.
    end subroutine new

    !> (h,k)-only quantities of the fixed lims2 disk (spatial frequency squared,
    !! astigmatism angle, resolution shell), built once instead of per particle
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

    !> Separable real-space KB instrument envelope: the gather multiplies the volume
    !! by the window's transform (why gridding deapodizes). The matrix-free operator
    !! applies it twice, H = E T E, so the kernel must bracket its convolution with
    !! the same envelope (~0.19 per axis at the box edge). Exact discrete transform
    !! of the normalized origin stencil (not kbinterpol%instr), separable by construction.
    subroutine build_env( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, parameter :: EPS_DIV = 1.0e-8
        real, allocatable :: env1d(:), env1d_padded(:), inv1d(:)
        real    :: ctrval
        integer :: c, i, j, k
        if( allocated(self%env)    ) deallocate(self%env)
        if( allocated(self%invenv) ) deallocate(self%invenv)
        call kb_stencil_envelope_1d(self%kbwin,self%boxpd,env1d_padded)
        call kb_stencil_centered_crop_inv_envelope_1d(self%kbwin,self%boxpd,self%box,inv1d)
        allocate(env1d(self%box),source=env1d_padded(self%pad_off+1:self%pad_off+self%box))
        allocate(self%env(self%box,self%box,self%box))
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    self%env(i,j,k) = env1d(i)*env1d(j)*env1d(k)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(env1d,env1d_padded)
        ! normalize so the envelope is unity at the box centre
        c      = self%box/2 + 1
        ctrval = self%env(c,c,c)
        if( abs(ctrval) < EPS_DIV ) ctrval = 1.0
        self%env = self%env / ctrval
        allocate(self%invenv(self%box,self%box,self%box), source=1.0)
        do k = 1, self%box
            do j = 1, self%box
                do i = 1, self%box
                    if( abs(self%env(i,j,k)) < EPS_DIV )then
                        self%invenv(i,j,k) = 0.0
                    else
                        self%invenv(i,j,k) = inv1d(i)*inv1d(j)*inv1d(k)
                    endif
                end do
            end do
        end do
        deallocate(inv1d)
    end subroutine build_env

    pure function get_env( self ) result( env )
        class(reconstructor_pcg), intent(in) :: self
        real :: env(self%box,self%box,self%box)
        env = self%env
    end function get_env

    !> copy of the RHS currently held (solve or end_accum), lets a test localize a disagreement
    subroutine get_rhs( self, b )
        class(reconstructor_pcg), intent(in)  :: self
        real, allocatable,        intent(out) :: b(:,:,:)
        if( .not. self%l_rhs ) THROW_HARD('no right-hand side has been built; get_rhs')
        allocate(b(self%box,self%box,self%box), source=self%b_rhs)
    end subroutine get_rhs

    !> copy of the open raw accumulators; test/diagnostic boundary only
    subroutine get_raw_accum( self, b, d )
        class(reconstructor_pcg), intent(in)  :: self
        complex, allocatable,     intent(out) :: b(:,:,:)
        real, allocatable,        intent(out) :: d(:,:,:)
        if( .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open')
        allocate(b(self%lims3(1,1):self%lims3(1,2), &
                   &self%lims3(2,1):self%lims3(2,2), &
                   &self%lims3(3,1):self%lims3(3,2)), source=self%b_work)
        allocate(d(self%lims3(1,1):self%lims3(1,2), &
                   &self%lims3(2,1):self%lims3(2,2), &
                   &self%lims3(3,1):self%lims3(3,2)), source=self%acc_work)
    end subroutine get_raw_accum

    subroutine get_ml_prior( self, prior )
        class(reconstructor_pcg), intent(in)  :: self
        real, allocatable,        intent(out) :: prior(:,:,:)
        if( .not. self%l_ml_prior ) THROW_HARD('PCG ML prior has not been built')
        allocate(prior, source=self%ml_prior)
    end subroutine get_ml_prior

    !> P_tau summary over its positive bins relative to the calibrated data-only Khat
    !! (|Khat| in the L1 denominator: the finite-support Toeplitz can hold small negative bins)
    subroutine get_ml_prior_stats( self, npositive, positive_min, positive_max, &
            &prior_to_khat_l1, prior_to_khat_rms )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                  intent(out) :: npositive
        real,                     intent(out) :: positive_min, positive_max, prior_to_khat_l1, prior_to_khat_rms
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
                end do
            end do
        end do
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
        logical,                  intent(in)    :: l_deapod
        self%l_deapod = l_deapod
    end subroutine set_deapod

    !> lambda as a dimensionless coefficient of the weighted data operator; the absolute
    !! value is derived once raw D is reduced, so workers publish unregularized (B,D)
    subroutine set_lambda_relative( self, lambda_rel )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: lambda_rel
        if( .not. ieee_is_finite(lambda_rel) .or. lambda_rel < 0.0 )then
            THROW_HARD('relative PCG lambda must be finite and non-negative')
        endif
        self%lambda_rel = lambda_rel
        self%data_scale = 0.0
        self%lambda     = 0.0
        self%l_lambda_relative = .true.
    end subroutine set_lambda_relative

    !> request the isotropic FSC/SSNR prior: the FSC sets the shell-wise relative
    !! strength, the absolute scale comes from the finalized D on the master
    subroutine set_ml_prior( self, fsc, tau, hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: fsc(:), tau, hp
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

    !> multiplies by E^-1, the inverse KB envelope (deapodization). Real images carry
    !! no envelope, so fitting them with A E returns E^-1 x; applying E^-1 on both
    !! sides of the normal operator and once to the RHS makes the solve target x
    !! itself, so the prior acts on x rather than on a rescaled surrogate
    pure subroutine deapod_mul( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                     intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_deapod ) return
        v = v * self%invenv
    end subroutine deapod_mul

    !> spherical support as a CONSTRAINT ON THE SOLVE: the window is image%mask3D_soft
    !! on a unit volume (backgr=0.), the same soft mask the gridding restoration
    !! applies; the solve runs on the domain window > 0 (install_support) and the
    !! shipped map is window*u. Removes the solvent, where deapodization amplifies
    !! hardest, and shrinks the problem (~18% of the box at mskdiam 180 in a 256 box).
    !! No consumer masks a second time; warm starts return to u through window_div
    subroutine set_mask( self, mskrad )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: mskrad
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
        call self%install_support
    end subroutine set_mask

    !> caller-supplied real-space [0,1] volume as the support P (clipped); same
    !! contract as set_mask. Experimental focused support (pcg_priors_history.md dev item 5)
    subroutine set_mask_volume( self, mskvol )
        class(reconstructor_pcg), intent(inout) :: self
        class(image),             intent(in)    :: mskvol
        integer :: mdim(3)
        mdim = mskvol%get_ldim()
        if( any(mdim /= self%box) )then
            THROW_HARD('support mask volume dimensions differ from the solve box; set_mask_volume')
        endif
        if( mskvol%is_ft() ) THROW_HARD('support mask volume must be in real space; set_mask_volume')
        if( allocated(self%mask)   ) deallocate(self%mask)
        if( allocated(self%window) ) deallocate(self%window)
        self%window = mskvol%get_rmat()
        self%window = min(1.0, max(0.0, self%window))
        if( .not. any(self%window > 0.0) ) THROW_HARD('support mask volume is empty; set_mask_volume')
        call self%install_support
    end subroutine set_mask_volume

    !> solve support from the window. With PCG_HARD_SOLVE_SUPPORT the domain is
    !! window > 0: P^2 = P, so the projections in operator, RHS and preconditioner
    !! are exact, u is the estimate on that domain and the shipped map window*u is
    !! one estimate times one soft window, exactly what the gridding restoration
    !! ships. The soft alternative (P = window) leaves the band a solver-state
    !! dependent mixture of P*u and u that no windowed estimate reproduces
    subroutine install_support( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( allocated(self%mask) ) deallocate(self%mask)
        if( PCG_HARD_SOLVE_SUPPORT )then
            allocate(self%mask(self%box,self%box,self%box), source=merge(1.0, 0.0, self%window > 0.0))
        else
            allocate(self%mask(self%box,self%box,self%box), source=self%window)
        endif
        self%l_mask = .true.
    end subroutine install_support

    pure subroutine mask_mul( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                     intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_mask ) return
        v = v * self%mask
    end subroutine mask_mul

    !> x = window*u, the shipped map
    pure subroutine window_mul( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                     intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_mask ) return
        v = v * self%window
    end subroutine window_mul

    !> output-space initial guess (x = window*u, every shipped half map) -> CG
    !! variable u = x/window where window >= PCG_SUPPORT_DIV_MIN, zero elsewhere.
    !! Projecting with the window instead squared the soft edge on every
    !! warm-started iteration and compounded over a stage. Exact for windowed input;
    !! the floor caps the amplification (1/floor) of content not proportional to
    !! the window (resampling ringing, foreign maps) to the outermost band voxels
    pure subroutine window_div( self, v )
        class(reconstructor_pcg), intent(in)    :: self
        real,                     intent(inout) :: v(self%box,self%box,self%box)
        if( .not. self%l_mask ) return
        where( self%window >= PCG_SUPPORT_DIV_MIN )
            v = v / self%window
        elsewhere
            v = 0.0
        end where
    end subroutine window_div

    subroutine ensure_wimg( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%wimg_exists ) return
        call self%wimg%new([self%boxpd,self%boxpd,self%boxpd], self%smpd, &
            &wthreads=self%fft_nthreads > 1, fft_nthreads=self%fft_nthreads)
        self%wimg_exists = .true.
    end subroutine ensure_wimg

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
        if( allocated(self%window)   ) deallocate(self%window)
        if( allocated(self%precond)  ) deallocate(self%precond)
        if( allocated(self%Khat)     ) deallocate(self%Khat)
        if( allocated(self%ml_fsc)   ) deallocate(self%ml_fsc)
        if( allocated(self%ml_prior) ) deallocate(self%ml_prior)
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
        self%wimg_exists = .false.
        self%op_mode     = PCG_OP_MATRIXFREE
        call self%reset_profile(.false.)
        call self%reset_finalize_profile
        self%exists      = .false.
    end subroutine kill

    ! SETUP

    !> caches the per-particle quantities that do not depend on the iterate
    !! (rotation matrices, CTF parameters, shifts, noise spectra); small, constant
    !! over the solve, and what lets the particle loop be shared across threads
    subroutine prep_particles( self, orientations, use_ctf, sig2 )
        class(reconstructor_pcg),           intent(inout) :: self
        class(oris),                        intent(inout) :: orientations
        logical,                  optional, intent(in)    :: use_ctf
        real,                     optional, intent(in)    :: sig2(0:,:)
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

    !> point-group operators for coordinate replication: each plane pixel is
    !! scattered at all R_i.S_g; symmats(:,:,1) = I reproduces c1. Call after new
    !! and before begin_accum
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
        integer,                  intent(in)    :: op_mode
        if( op_mode == PCG_OP_KERNEL .and. .not. self%l_kernel )then
            THROW_HARD('kernelized operator requested but build_kernel has not been called; set_op_mode')
        endif
        self%op_mode = op_mode
    end subroutine set_op_mode

    ! LOW-LEVEL OPERATOR

    !> loads a native box^3 volume, centre-padded and transformed onto the
    !! oversampled lattice; the padding never leaks to callers
    subroutine set_volume( self, v )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: v(self%box,self%box,self%box)
        call self%ensure_wimg
        call self%wimg%set_rmat(center_embed_real3d(v,self%boxpd), .false.)
        call self%wimg%fft()
    end subroutine set_volume

    !> G_i F: gathers a full (unpacked) Fourier plane at orientation e from the
    !! volume held by set_volume; native (h,k) -> padded padf*loc, padsc restores
    !! the native scale. Ori-based signature kept for the adjoint test
    subroutine forward_plane( self, e, plane )
        class(reconstructor_pcg), intent(inout) :: self
        class(ori),               intent(in)    :: e
        complex,                  intent(out)   :: plane(self%lims2(1,1):self%lims2(1,2),&
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
                plane(h,k) = self%padsc*gather_packed_window(cmat,lbound(self%wrap,1),self%wrap,i0,w)
            end do
        end do
        !$omp end parallel do
    end subroutine forward_plane

    !> <x_hat, accum> over the oversampled lattice for the volume held by
    !! set_volume; lets the adjoint test form <x, G^dagger q> without lattice knowledge
    function fourier_dot( self, accum ) result( d )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                  intent(in)    :: accum(self%lims3(1,1):self%lims3(1,2),&
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

    !> F^dagger G_i^dagger, accumulate form: the literal transpose of the gather,
    !! verified against forward_plane by the adjoint dot-product test
    subroutine adjoint_plane_add( self, plane, e, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        complex,                  intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                        &self%lims2(2,1):self%lims2(2,2))
        class(ori),               intent(in)    :: e
        complex,                  intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        real :: e_rotmat(3,3)
        e_rotmat = e%get_mat()
        call scatter_plane(self, plane, e_rotmat, vol_accum)
    end subroutine adjoint_plane_add

    !> T_i(h,k) = C_i(h,k) S_i(h,k) / sqrt(sigma2_i(shell)), the full complex
    !! transfer; needed for the RHS only (absT2_plane serves the normal operator)
    function build_transfer( self, ctfparms, shift, sig2arr ) result( T )
        class(reconstructor_pcg),           intent(in) :: self
        type(ctfparams),                    intent(in) :: ctfparms
        real,                               intent(in) :: shift(2)
        real,                     optional, intent(in) :: sig2arr(0:)
        complex       :: T(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real          :: cval, args, sw, sum_df, diff_df, angast, wl, half_wl2_cs, accc, phc, cterm, df, phsh, s2
        integer       :: h, k, shell
        logical       :: l_ctf, l_flip
        ! ctfflag as image%gen_fplane4rec reads it: NO = no CTF to model, FLIP =
        ! already phase-flipped (a signed CTF would reintroduce the removed phase)
        l_ctf  = ctfparms%ctfflag /= CTFFLAG_NO
        l_flip = ctfparms%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
            ! flat, call-free CTF form of ft_map_ctf_kernel over the both-sign-h disk
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
                ! POSITIVE shift phase: apply_adjoint_all applies conjg(T), which
                ! reproduces production's exp(-i shift) correction (gen_fplane4rec)
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

    !> y_w = y / sqrt(sigma2(shell)), the production inverse-noise amplitude
    function whiten_observation( self, plane, sig2arr ) result( whitened )
        class(reconstructor_pcg), intent(in) :: self
        complex,                  intent(in) :: plane(self%lims2(1,1):self%lims2(1,2), &
                                                     &self%lims2(2,1):self%lims2(2,2))
        real,                     intent(in) :: sig2arr(0:)
        complex :: whitened(self%lims2(1,1):self%lims2(1,2), &
                           &self%lims2(2,1):self%lims2(2,2))
        real    :: sigma2
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
            end do
        end do
    end function whiten_observation

    !> reads an FFT'd particle image's own Fourier plane into the lims2 disk (no
    !! KB window, no wrap); how y_planes is built for real particles
    function extract_native_plane( self, img2d ) result( plane )
        class(reconstructor_pcg), intent(in) :: self
        class(image),             intent(in) :: img2d
        complex :: plane(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2))
        plane = extract_native_fourier_plane(img2d,self%lims2,self%sqlp)
    end function extract_native_plane

    !> deterministic double-precision dot product
    pure function dot_real_volume( self, a, b ) result( d )
        class(reconstructor_pcg), intent(in) :: self
        real,                     intent(in) :: a(self%box,self%box,self%box)
        real,                     intent(in) :: b(self%box,self%box,self%box)
        real(dp) :: d
        d = sum(real(a,dp) * real(b,dp))
    end function dot_real_volume

    ! HIGH-LEVEL OPERATOR

    !> the operator the solver sees: P H P with a support (set_mask), H otherwise.
    !! The two concrete operators stay unmasked so the tests compare them fully
    function apply_normal( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: p(self%box,self%box,self%box)
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

    !> H p = sum_i G_i^dagger |T_i|^2 G_i p + lambda p, the exact reference operator.
    !! Gather, weight and scatter are fused per plane point. All threads walk the
    !! particle loop in lockstep with the h loop workshared; the h-strided colouring
    !! (stride, see new) keeps the scatter footprints disjoint in one shared accumulator
    function apply_normal_matrixfree( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: p(self%box,self%box,self%box)
        real,    allocatable :: hp(:,:,:), pd(:,:,:)
        complex, allocatable :: vol_accum(:,:,:), cmat(:,:,:)
        real,    allocatable :: absT2(:,:)
        real                    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3)
        complex                 :: comp
        integer                 :: i, g, h, k, l, i0(3)
        integer(timer_int_kind) :: tp
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; apply_normal_matrixfree')
        call self%ensure_wimg
        ! A = A~ E^-1: deapodize in and out of the adjoint, so the operator is D.S not D.S.E
        allocate(pd(self%box,self%box,self%box), source=p)
        call self%deapod_mul(pd)
        if( self%l_profile ) tp = pcg_tic()
        call self%set_volume(pd)
        if( self%l_profile ) self%t_setvol = self%t_setvol + pcg_toc(tp)
        ! get_cmat returns a COPY (539 MB at box 256); t_cmatcp isolates that traffic
        if( self%l_profile ) tp = pcg_tic()
        cmat = self%wimg%get_cmat()
        if( self%l_profile ) self%t_cmatcp = self%t_cmatcp + pcg_toc(tp)
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        allocate(absT2(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        if( self%l_profile ) tp = pcg_tic()
        !$omp parallel default(shared) private(i,g,h,k,l,loc,i0,w,comp,rot) proc_bind(close)
        do i = 1, self%nptcls
            ! orphaned worksharing inside absT2_plane binds to this region; its
            ! trailing barrier makes absT2 safe to read below
            call self%absT2_plane(i, absT2)
            ! replication at every R_i.S_g, g outside the colour sweep (matches accumulate_absT2)
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
                            comp = self%padsc*gather_packed_window(cmat,lbound(self%wrap,1),self%wrap,i0,w)
                            comp = comp * absT2(h,k) * self%padsc
                            call scatter_window(self, i0, w, comp, vol_accum)
                        end do
                    end do
                    !$omp end do
                end do
                ! wrapping rim, serialized: the colouring does not survive folding (win_wraps)
                !$omp single
                do h = self%lims2(1,1), self%lims2(1,2)
                    do k = self%lims2(2,1), self%lims2(2,2)
                        if( h*h + k*k > self%sqlp   ) cycle
                        if( h*h + k*k <= self%sq_rim ) cycle   ! provably cannot wrap
                        loc  = real(self%padf) * matmul(real([h,k,0]), rot)
                        i0   = nint(loc) - self%iwinsz
                        if( .not. win_wraps(self, i0) ) cycle
                        call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                        comp = self%padsc*gather_packed_window(cmat,lbound(self%wrap,1),self%wrap,i0,w)
                        comp = comp * absT2(h,k) * self%padsc
                        call scatter_window(self, i0, w, comp, vol_accum)
                    end do
                end do
                !$omp end single
            end do
        end do
        !$omp end parallel
        if( self%l_profile ) self%t_ploop = self%t_ploop + pcg_toc(tp)
        hp = self%fold_and_ifft(vol_accum)
        call self%deapod_mul(hp)
        if( self%l_ml_prior ) hp = hp + self%apply_fourier_diagonal(p, self%ml_prior)
        hp = hp + self%lambda * p
    end function apply_normal_matrixfree

    !> kernelized (Toeplitz/Gram) operator H_data p = crop(IFFT(Khat FFT(pad p))),
    !! O(box^3 log box) and independent of the particle count; see build_kernel
    function apply_normal_kernel( self, p ) result( hp )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: p(self%box,self%box,self%box)
        real,    allocatable :: hp(:,:,:), work(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        real                    :: kv
        integer                 :: cdim(3), i, j, k
        integer(timer_int_kind) :: tp
        if( .not. self%l_kernel ) THROW_HARD('build_kernel has not been called; apply_normal_kernel')
        call self%ensure_wimg
        allocate(work(self%box,self%box,self%box), source=p)
        ! Khat is the bare Toeplitz T; with deapodization off the matrix-free
        ! operator is E T E, so the envelope is reinstated on both sides
        if( .not. self%l_deapod ) work = work * self%env
        if( self%l_profile ) tp = pcg_tic()
        call self%wimg%set_rmat(center_embed_real3d(work,self%boxpd), .false.)
        call self%wimg%fft()
        if( self%l_profile ) self%t_setvol = self%t_setvol + pcg_toc(tp)
        if( self%l_profile ) tp = pcg_tic()
        cmat = self%wimg%get_cmat()
        if( self%l_profile ) self%t_cmatcp = self%t_cmatcp + pcg_toc(tp)
        cdim = self%wimg%get_array_shape()
        if( self%l_profile ) tp = pcg_tic()
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
        if( self%l_profile ) self%t_khat = self%t_khat + pcg_toc(tp)
        if( self%l_profile ) tp = pcg_tic()
        call self%wimg%set_cmat(cmat)
        call self%wimg%ifft()
        hp = center_crop_real3d(self%wimg%get_rmat(),self%box)
        if( self%l_profile ) self%t_fold = self%t_fold + pcg_toc(tp)
        if( .not. self%l_deapod ) hp = hp * self%env
        if( self%l_ml_prior .and. .not. self%l_deapod )then
            hp = hp + self%apply_fourier_diagonal(p, self%ml_prior)
        endif
        ! band precision on the deapodized domain; attachment mode enforced upstream
        hp = hp + self%lambda * p
    end function apply_normal_kernel

    !> C^T F^-1 diag(d) F C on the Khat lattice; the matrix-free oracle and the
    !! deapodization-off path use it, kernel solves fuse d with Khat
    function apply_fourier_diagonal( self, p, diag ) result( q )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: p(self%box,self%box,self%box)
        real,                     intent(in)    :: diag(:,:,:)
        real,    allocatable :: q(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        integer :: cdim(3), i, j, k
        call self%ensure_wimg
        cdim = self%wimg%get_array_shape()
        if( any(shape(diag) /= cdim) ) THROW_HARD('PCG Fourier diagonal shape mismatch')
        call self%wimg%set_rmat(center_embed_real3d(p,self%boxpd), .false.)
        call self%wimg%fft()
        cmat = self%wimg%get_cmat()
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    cmat(i,j,k) = cmat(i,j,k) * diag(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do
        call self%wimg%set_cmat(cmat)
        call self%wimg%ifft()
        q = center_crop_real3d(self%wimg%get_rmat(),self%box)
    end function apply_fourier_diagonal

    !> b = sum_i G_i^dagger(conjg(T_i) y_i / sqrt(sigma2_i)), the data RHS; unlike
    !! H it needs the full complex T_i including the shift phase
    function apply_adjoint_all( self, y_planes ) result( b )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                  intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
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
        ! A^dagger = E^-1 A~^dagger: the RHS gets the correction once, H twice
        call self%deapod_mul(b)
    end function apply_adjoint_all


    ! STREAMING SETUP: begin_accum / accumulate_batch / end_accum form the RHS and
    ! the density in one pass over particle batches, so no plane stays resident;
    ! peak memory is one complex and one real full-range accumulator plus the work image

    subroutine begin_accum( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; begin_accum')
        call self%begin_reduction
    end subroutine begin_accum

    !> zero raw accumulators without particle metadata: masters call this before
    !! adding worker artifacts, workers use begin_accum after prep_particles
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

    !> accumulates one batch (particles ifrom .. ifrom+nb-1 of the selection) into both accumulators
    subroutine accumulate_batch( self, y_batch, nb, ifrom )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                  intent(in)    :: nb, ifrom
        complex,                  intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                          &self%lims2(2,1):self%lims2(2,2), nb)
        complex, allocatable :: bacc(:,:,:)
        real,    allocatable :: dacc(:,:,:)
        if( .not. self%l_accum ) THROW_HARD('begin_accum has not been called; accumulate_batch')
        if( nb < 1 ) return
        if( ifrom < 1 .or. ifrom + nb - 1 > self%nptcls )then
            THROW_HARD('batch range outside the particle selection; accumulate_batch')
        endif
        ! move_alloc detaches the accumulators from self before they are passed
        ! next to intent(inout) self (an aliasing the standard lets compilers exploit)
        call move_alloc(self%b_work, bacc)
        call move_alloc(self%acc_work, dacc)
        call self%accumulate_rhs_density(y_batch, nb, ifrom, bacc, dacc)
        call move_alloc(bacc, self%b_work)
        call move_alloc(dacc, self%acc_work)
    end subroutine accumulate_batch

    !> atomically publishes one worker's raw full-range B and D (no folding,
    !! deapodization, flooring or solve before this); header = manifest, .tmp promoted after close
    subroutine write_raw_accum( self, fname, state, eo, part, nparts, nptcls, provenance )
        class(reconstructor_pcg), intent(in) :: self
        class(string),            intent(in) :: fname
        integer,                  intent(in) :: state, eo, part, nparts, nptcls
        character(len=*),         intent(in) :: provenance
        type(string)                    :: tmpfname
        character(len=PCG_RAW_PROV_LEN) :: prov_fixed
        integer                         :: funit, ierr, m
        integer(int64)                  :: file_size
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
            end do
            call fileiochk('write_raw_accum writing B', ierr)
            do m = self%lims3(3,1), self%lims3(3,2)
                write(funit, iostat=ierr) self%acc_work(:,:,m)
                if( ierr /= 0 ) exit
            end do
            call fileiochk('write_raw_accum writing D', ierr)
        endif
        call fclose(funit)
        file_size = -1_int64
        inquire(file=tmpfname%to_char(), size=file_size)
        if( file_size <= 0_int64 ) THROW_HARD('raw PCG artifact write produced an empty file')
        call simple_rename(tmpfname, fname, overwrite=.true.)
        call tmpfname%kill
    end subroutine write_raw_accum

    !> adds one raw worker artifact to the open reduction; parts arrive in ascending
    !! order (reproducible association), slices are streamed
    subroutine add_raw_accum( self, fname, state, eo, part, nparts, provenance, nptcls )
        class(reconstructor_pcg), intent(inout) :: self
        class(string),            intent(in)    :: fname
        integer,                  intent(in)    :: state, eo, part, nparts
        character(len=*),         intent(in)    :: provenance
        integer,                  intent(out)   :: nptcls
        character(len=16)               :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        complex, allocatable :: bslice(:,:)
        real,    allocatable :: dslice(:,:)
        integer        :: funit, ierr, m, version, state_file, eo_file, part_file
        integer        :: nparts_file, box_file, boxpd_file, padf_file, lims_file(3,2)
        real           :: smpd_file
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
            end do
            call fileiochk('add_raw_accum reading B', ierr)
            do m = self%lims3(3,1), self%lims3(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
                self%acc_work(:,:,m) = self%acc_work(:,:,m) + dslice
            end do
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

    !> scales an open raw accumulator: the only valid place for the u/f and (1-u)
    !! continuation weights (finalized kernels and maps are not additive)
    subroutine scale_raw_accum( self, weight )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: weight
        if( .not. self%l_accum ) THROW_HARD('raw PCG accumulator is not open; scale_raw_accum')
        if( .not. ieee_is_finite(weight) .or. weight < 0.0 ) THROW_HARD('invalid raw PCG accumulator weight')
        self%b_work   = weight * self%b_work
        self%acc_work = weight * self%acc_work
    end subroutine scale_raw_accum

    !> relative L2 differences between two open raw accumulators, without copies
    subroutine compare_raw_accum( self, other, b_relerr, d_relerr )
        class(reconstructor_pcg), intent(in)  :: self, other
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
                end do
            end do
        end do
        b_relerr = real(sqrt(bnum / max(1.0_dp,bden)))
        d_relerr = real(sqrt(dnum / max(1.0_dp,dden)))
    end subroutine compare_raw_accum

    !> header compatibility of a persisted raw artifact with the current geometry and
    !! identity. Constant-FOV crop continuity: a SMALLER crop with the same field of
    !! view is compatible (zero-extended by add_raw_accum_weighted); a larger crop, a
    !! FOV change, an identity change or a read failure is not (discard and re-seed)
    logical function pcg_raw_accum_compatible( fname, box, smpd, provenance ) result( l_compatible )
        class(string),    intent(in) :: fname
        integer,          intent(in) :: box
        real,             intent(in) :: smpd
        character(len=*), intent(in) :: provenance
        character(len=16)               :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        integer                         :: funit, ierr, version, state_file, eo_file, part_file
        integer                         :: nparts_file, nptcls_file
        integer                         :: box_file, boxpd_file, padf_file, lims_file(3,2)
        real                            :: smpd_file, fov_file, fov_cur
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

    !> adds one complete raw artifact with an explicit weight; deterministic algebra
    !! on reduced chains, outside the worker-part ordering. Under constant FOV the
    !! padded lattices of consecutive crops share their frequency step, so a smaller
    !! previous grid embeds exactly by index-aligned zero-extension (the old wrap rim
    !! carries aliased mass at/beyond the producing band, decaying as (1-u)^k)
    subroutine add_raw_accum_weighted( self, fname, state, eo, part, nparts, provenance, weight, nptcls )
        class(reconstructor_pcg), intent(inout) :: self
        class(string),            intent(in)    :: fname
        integer,                  intent(in)    :: state, eo, part, nparts
        character(len=*),         intent(in)    :: provenance
        real,                     intent(in)    :: weight
        integer,                  intent(out)   :: nptcls
        character(len=16)               :: magic
        character(len=PCG_RAW_PROV_LEN) :: prov_file, prov_expected
        complex, allocatable :: bslice(:,:)
        real,    allocatable :: dslice(:,:)
        integer        :: funit, ierr, m, version, state_file, eo_file, part_file
        integer        :: nparts_file, box_file, boxpd_file, padf_file, lims_file(3,2)
        real           :: smpd_file, fov_file, fov_self
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
        if( magic /= PCG_RAW_ACCUM_MAGIC .or. version /= PCG_RAW_ACCUM_VERSION )then
            THROW_HARD('weighted raw PCG format mismatch')
        endif
        if( state_file /= state .or. eo_file /= eo .or. part_file /= part )then
            THROW_HARD('weighted raw PCG identity mismatch')
        endif
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
            end do
            call fileiochk('add_raw_accum_weighted reading B', ierr)
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
                self%acc_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) = &
                    &self%acc_work(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2), m) + weight * dslice
            end do
            call fileiochk('add_raw_accum_weighted reading D', ierr)
            deallocate(bslice, dslice)
        else if( nptcls > 0 )then
            allocate(bslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            allocate(dslice(lims_file(1,1):lims_file(1,2), lims_file(2,1):lims_file(2,2)))
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) bslice
                if( ierr /= 0 ) exit
            end do
            call fileiochk('add_raw_accum_weighted skipping B', ierr)
            do m = lims_file(3,1), lims_file(3,2)
                read(funit, iostat=ierr) dslice
                if( ierr /= 0 ) exit
            end do
            call fileiochk('add_raw_accum_weighted skipping D', ierr)
            deallocate(bslice, dslice)
        endif
        file_size  = -1_int64
        stream_pos = -1_int64
        inquire(unit=funit, size=file_size, pos=stream_pos)
        call fclose(funit)
        if( stream_pos-1_int64 /= file_size ) THROW_HARD('weighted raw PCG accumulator has trailing or missing bytes')
    end subroutine add_raw_accum_weighted

    !> accumulates the weighted RHS B and the density D in one traversal: same
    !! orientation, coordinate and KB window for both scatters
    subroutine accumulate_rhs_density( self, y_batch, nb, ifrom, bacc, dacc )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                  intent(in)    :: nb, ifrom
        complex,                  intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                          &self%lims2(2,1):self%lims2(2,2), nb)
        complex,                  intent(inout) :: bacc(self%lims3(1,1):self%lims3(1,2),&
                                                       &self%lims3(2,1):self%lims3(2,2),&
                                                       &self%lims3(3,1):self%lims3(3,2))
        real,                     intent(inout) :: dacc(self%lims3(1,1):self%lims3(1,2),&
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

    !> closes accumulation: folds the RHS once, derives preconditioner and optionally Khat
    subroutine end_accum( self, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                  intent(in)    :: l_kernel
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
        ! RHS first and freed before the allocation-heavy kernel work; move_alloc as in accumulate_batch
        tp = pcg_tic()
        call move_alloc(self%b_work, bwork)
        self%b_rhs = self%fold_and_ifft(bwork)
        deallocate(bwork)
        call self%deapod_mul(self%b_rhs)
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        self%t_fin_rhs = pcg_toc(tp)
        call move_alloc(self%acc_work, dwork)
        call self%finalize_density_accum(dwork, .true., l_kernel)
        deallocate(dwork)
        if( l_kernel ) call self%finalize_khat
        self%l_accum = .false.
    end subroutine end_accum

    ! SETUP: PRECONDITIONER AND KERNEL

    !> sampling-density preconditioner M(k) = rho(k) + floor(shell), rho =
    !! sum_i G_i^dagger |T_i|^2 (the gridding density); cuts the iteration count that
    !! heterogeneous CTFs would otherwise inflate. The floor is a fraction of the
    !! shell-mean rho (RHO_FLOOR_FRAC), never an absolute constant: rho spans many
    !! orders of magnitude and is zero between rotated planes, so 1/(rho+lambda)
    !! would amplify the least-constrained modes by orders of magnitude and fill the
    !! map with noise. Scale-invariant (rho lacks the padsc factors of apply_normal)
    subroutine build_precond( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_precond')
        call self%ensure_wimg
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        call self%accumulate_absT2(acc)
        call self%finalize_density_accum(acc, .true., .false.)
    end subroutine build_precond

    !> the one particle pass both the preconditioner and the Gram kernel derive
    !! from: |T_i|^2 scattered at padf*R_i*[h,k,0] through the KB stencil. Both
    !! consumers are scale-invariant (shell-relative floor, analytic calibration),
    !! so sharing the accumulator is exact. One parallel region for all particles;
    !! ifrom/ito keep a batch caller in step with the RHS (no image data needed)
    subroutine accumulate_absT2( self, acc, ifrom, ito )
        class(reconstructor_pcg),           intent(inout) :: self
        real,                               intent(inout) :: acc(self%lims3(1,1):self%lims3(1,2),&
                                                                &self%lims3(2,1):self%lims3(2,2),&
                                                                &self%lims3(3,1):self%lims3(3,2))
        integer,                  optional, intent(in)    :: ifrom, ito
        real, allocatable :: absT2(:,:)
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
            ! replication at every R_i.S_g, g outside the colour sweep; symmats(:,:,1)=I gives c1
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

    !> absolute Tikhonov coefficient from a homogeneous scale of the raw density:
    !! the mean over the first six native shells (common to full and cropped
    !! representations), extended one shell at a time when that band is empty
    !! (some Euclidean/simulated weights suppress it). Zero bins stay in the average
    !! so fractional updates remain linear; padsc**2 converts to the Khat convention
    subroutine update_lambda_from_density( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
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
        end do
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
    end subroutine update_lambda_from_density

    !> P_tau from the independent-half FSC and the raw density: each padded radius
    !! maps to its nearest native shell, the shell-mean D supplies 1/sigma2, and
    !! padsc**2 puts the diagonal in the calibrated Khat convention
    subroutine build_ml_prior_from_density( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2), &
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
                end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        do ithr = 1, nthr
            do sh = 0, sz
                shsum(sh) = shsum(sh) + shsum_thr(sh,ithr)
                shcnt(sh) = shcnt(sh) + shcnt_thr(sh,ithr)
            end do
        end do
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
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(shsum, shcnt)
        self%l_ml_prior = maxval(self%ml_prior) > 0.0
        if( .not. self%l_ml_prior ) THROW_HARD('PCG ML prior contains no positive bins')
    end subroutine build_ml_prior_from_density

    !> folds the real density accumulator once into the packed layout: preconditioner
    !! and Khat are two views of D, produced in one sphere-limited pass; per-thread
    !! shell sums merged in thread order for a stable result
    subroutine finalize_density_accum( self, rho_accum, l_precond, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        logical,                  intent(in)    :: l_precond, l_kernel
        real, parameter :: RHO_FLOOR_FRAC = 1.0e-2
        real(dp), allocatable :: shsum(:), shsum_thr(:,:)
        integer,  allocatable :: shcnt(:), shcnt_thr(:,:)
        real,     allocatable :: shfloor(:)
        integer                 :: h, k, cdim(3), m, hh, phys(3), sh, rho_lim
        integer                 :: khat_lim, work_lim, nthr, ithr
        real                    :: denom, rval
        integer(timer_int_kind) :: tp
        if( .not. l_precond .and. .not. l_kernel ) return
        call self%update_lambda_from_density(rho_accum)
        cdim = self%wimg%get_array_shape()
        call self%build_ml_prior_from_density(rho_accum)
        ! data reaches |loc| <= padf*Rnat only; beyond it rho is zero and the modes
        ! unconstrained. A singular (PSD) preconditioner keeps the Krylov space out of them
        rho_lim = self%padf * self%Rnat
        ! D extends past the source sphere only by rounding and the finite stencil (exact bound)
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
        ! PASS 1: shell-mean rho over the sampled voxels only (zeros would drag the floor down)
        if( l_precond )then
            tp = pcg_tic()
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
            self%t_fin_rho = pcg_toc(tp)
        endif
        ! PASS 2: guarded reciprocal and packed Khat; an empty shell stays unconstrained
        tp = pcg_tic()
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
                            if( denom > 0.0 )then
                                self%precond(phys(1),phys(2),phys(3)) = 1.0 / denom
                            endif
                        endif
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        self%t_fin_fold = pcg_toc(tp)
        if( l_precond )then
            self%l_precond = .true.
            deallocate(shfloor)
        endif
    end subroutine finalize_density_accum

    !> Gram kernel. Not the impulse response of the matrix-free operator: the KB
    !! weights sum to 1, so that reproduces the gridding density and PCG would
    !! converge to gridding in one step. Standard NUFFT construction instead:
    !! |T_i|^2 scattered onto the 2x oversampled grid at doubled coordinates, the
    !! scattered array itself being the multiplier (H p = crop(IFFT(Khat FFT(pad p)))).
    !! Khat is real and symmetric (full symmetric lims2 disk), so the packed
    !! half-grid is exact. Memory: (2 box)^3, ~1.1 GB at box 256
    subroutine build_kernel( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_kernel')
        call self%ensure_wimg
        ! the padf*box lattice of the operator is the 2x grid the linear convolution needs
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=0.0)
        call self%accumulate_absT2(acc)
        call self%finalize_density_accum(acc, .false., .true.)
        ! freed before finalize_khat, the memory-heavy half
        deallocate(acc)
        call self%finalize_khat
    end subroutine build_kernel

    !> preconditioner and kernel from a single particle pass (the solver path);
    !! build_precond and build_kernel stay separately callable for the tests
    subroutine build_operators( self, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                  intent(in)    :: l_kernel
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

    subroutine finalize_khat( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, parameter :: EPS_D = 1.0e-8
        complex, allocatable :: ctmp(:,:,:)
        real,    allocatable :: tker(:,:,:), dep1d(:)
        real                    :: depval
        integer                 :: i, j, k, cdim(3)
        integer(timer_int_kind) :: tp
        cdim = self%wimg%get_array_shape()
        ! divide out the DEPOSITION envelope: scattering |T|^2 through the KB window
        ! multiplies the real-space kernel by the window's transform, a second
        ! envelope distinct from the gather's; separable exact stencil transform
        tp = pcg_tic()
        call kb_stencil_envelope_1d(self%kbwin,self%boxpd,dep1d)
        self%t_fin_dep = pcg_toc(tp)
        tp = pcg_tic()
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
        self%t_fin_kernel = pcg_toc(tp)
    end subroutine finalize_khat

    !> analytic scale of Khat relative to the matrix-free operator: padsc**2. A
    !! least-squares fit measured 64.3 (synthetic) and 63.98 (real data), i.e. this
    !! constant; the trailing path derives Khat from a stored accumulator with no
    !! particles resident, so only an analytic factor survives (policy section 5).
    !! measure_kernel_scale still fits it for the tests
    subroutine calibrate_kernel( self )
        class(reconstructor_pcg), intent(inout) :: self
        self%Khat = self%Khat * self%padsc**2
    end subroutine calibrate_kernel

    !> least-squares scale of the kernel against the matrix-free reference on a
    !! Gaussian probe; 1.0 when calibrate_kernel is right. One particle pass, tests only
    function measure_kernel_scale( self ) result( scale )
        class(reconstructor_pcg), intent(inout) :: self
        real, allocatable :: probe(:,:,:), hm(:,:,:), hk(:,:,:)
        real(dp) :: num, den
        real     :: lam_save, ctr, sig, dx, dy, dz, scale
        integer  :: i, j, k
        logical  :: l_ml_save
        if( .not. self%l_kernel ) THROW_HARD('build_kernel has not been called; measure_kernel_scale')
        lam_save    = self%lambda
        l_ml_save   = self%l_ml_prior
        self%lambda = 0.0   ! compare the DATA term only
        self%l_ml_prior = .false.
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
    end function measure_kernel_scale

    ! PRIVATE HELPERS

    !> |T_i|^2 = |C_i|^2 / sigma2_i over the lims2 disk. The unit-modulus shift
    !! phase cancels in conjg(T) T, so the normal operator needs only this real,
    !! iteration-invariant weight (the RHS needs the complex T, transfer_plane_cmplx).
    !! The worksharing loops are ORPHANED: inside a parallel region the CTF work is
    !! spread across the team, outside one it runs serially; every element is
    !! written inside a workshared loop
    subroutine absT2_plane( self, iptcl, absT2 )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                  intent(in)  :: iptcl
        real,                     intent(out) :: absT2(self%lims2(1,1):self%lims2(1,2),&
                                                      &self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        real          :: cval, sum_df, diff_df, angast, wl, half_wl2_cs, accc, phc, cterm, df, phsh, s2
        integer       :: h, k, shell
        logical       :: l_ctf, l_flip
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
        ! same ctfflag semantics as build_transfer; |CTF|^2 = CTF^2, so only NO matters here
        l_ctf  = self%ctfparms(iptcl)%ctfflag /= CTFFLAG_NO
        l_flip = self%ctfparms(iptcl)%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(self%ctfparms(iptcl)%smpd, self%ctfparms(iptcl)%kv, &
                &self%ctfparms(iptcl)%cs, self%ctfparms(iptcl)%fraca)
            call tfun%init(self%ctfparms(iptcl)%dfx, self%ctfparms(iptcl)%dfy, self%ctfparms(iptcl)%angast)
            ! flat, call-free form of ft_map_ctf_kernel over the both-sign-h disk (LUTs)
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

    !> full complex T_i of cached particle iptcl (RHS only)
    subroutine transfer_plane_cmplx( self, iptcl, T )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                  intent(in)  :: iptcl
        complex,                  intent(out) :: T(self%lims2(1,1):self%lims2(1,2),&
                                                  &self%lims2(2,1):self%lims2(2,2))
        T = self%build_transfer(self%ctfparms(iptcl), self%shifts(:,iptcl), self%sig2(:,iptcl))
    end subroutine transfer_plane_cmplx

    !> the two plane values of the fused accumulation, CTF evaluated once; every
    !! thread of the persistent region calls this and the orphaned omp-do partitions the plane
    subroutine prepare_fused_planes( self, iptcl, y_plane, weighted, absT2 )
        class(reconstructor_pcg), intent(in)  :: self
        integer,                  intent(in)  :: iptcl
        complex,                  intent(in)  :: y_plane(self%lims2(1,1):self%lims2(1,2),&
                                                        &self%lims2(2,1):self%lims2(2,2))
        complex,                  intent(out) :: weighted(self%lims2(1,1):self%lims2(1,2),&
                                                         &self%lims2(2,1):self%lims2(2,2))
        real,                     intent(out) :: absT2(self%lims2(1,1):self%lims2(1,2),&
                                                      &self%lims2(2,1):self%lims2(2,2))
        type(ctf)     :: tfun
        type(ctfvars) :: ctfvals
        complex       :: tval
        real          :: cval, arg, sw, sum_df, diff_df, angast, wl, half_wl2_cs
        real          :: accc, phc, cterm, df, phsh, s2
        integer       :: h, k, shell
        logical       :: l_ctf, l_flip
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

    !> folds a full-range (both-sign h) accumulator into the packed h>=0 storage and
    !! inverse-FFTs it. The redundant Nyquist mate h = lims3(1,2) is never produced by
    !! the wrapping scatter; its value lives at -lims3(1,2), where the wrap table sends it
    function fold_and_ifft( self, vol_accum ) result( z )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                  intent(in)    :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        real, allocatable :: z(:,:,:)
        integer                 :: h, hh, k, m, phys(3)
        integer(timer_int_kind) :: tp
        call self%ensure_wimg
        if( self%l_profile ) tp = pcg_tic()
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
        z = center_crop_real3d(self%wimg%get_rmat(),self%box)
        if( self%l_profile ) self%t_fold = self%t_fold + pcg_toc(tp)
    end function fold_and_ifft

    !> z = M^-1 r via FFT, diagonal multiply, inverse FFT. rho is the diagonal of
    !! the bare operator T while the solve targets H = E^-1 T E^-1, so H^-1 =
    !! E T^-1 E: the real-space envelope brackets the Fourier divide (env, not
    !! invenv); omitting it mis-preconditions by E^2
    function apply_precond( self, r ) result( z )
        class(reconstructor_pcg), intent(inout) :: self
        real,                     intent(in)    :: r(self%box,self%box,self%box)
        real,    allocatable :: z(:,:,:), rw(:,:,:)
        complex, allocatable :: cmat(:,:,:)
        integer                 :: cdim(3), i, j, k
        integer(timer_int_kind) :: tp
        if( .not. self%l_precond )then
            allocate(z(self%box,self%box,self%box), source=r)
            call self%mask_mul(z)
            return
        endif
        call self%ensure_wimg
        if( self%l_profile ) tp = pcg_tic()
        allocate(rw(self%box,self%box,self%box), source=r)
        if( self%l_deapod ) rw = rw * self%env
        call self%wimg%set_rmat(center_embed_real3d(rw,self%boxpd), .false.)
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
        z = center_crop_real3d(self%wimg%get_rmat(),self%box)
        if( self%l_deapod ) z = z * self%env
        ! z inside the support keeps the whole Krylov space there (M^-1 is a Fourier diagonal)
        call self%mask_mul(z)
        if( self%l_profile ) self%t_prec = self%t_prec + pcg_toc(tp)
    end function apply_precond

    !> does the window straddle the periodic wrap boundary? The h-strided colouring
    !! separates same-colour h-lines by more than a window only in unwrapped
    !! coordinates; after folding two such windows can overlap at the Nyquist rim
    !! and the sum becomes thread-order dependent. Callers scatter the interior in
    !! the coloured parallel pass and the wrapping rim serially. Module-level, not
    !! type-bound, so it inlines (~1e8 calls per accumulation), like the scatters
    pure logical function win_wraps( self, i0 )
        class(reconstructor_pcg), intent(in) :: self
        integer,                  intent(in) :: i0(3)
        win_wraps = any(i0 < self%wlims(1)) .or. any(i0 + self%wdim - 1 > self%wlims(2))
    end function win_wraps

    !> KB-weighted scatter of one complex value into the full-range accumulator
    pure subroutine scatter_window( self, i0, w, val, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
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
        class(reconstructor_pcg), intent(in)    :: self
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
        class(reconstructor_pcg), intent(in)    :: self
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
        class(reconstructor_pcg), intent(in)    :: self
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
        class(reconstructor_pcg), intent(in)    :: self
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
        class(reconstructor_pcg), intent(in)    :: self
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

    !> scatters a whole plane, h-strided so it is safe inside a parallel region
    !! (non-fused paths: adjoint_plane_add, apply_adjoint_all)
    subroutine scatter_plane( self, plane, rot, vol_accum )
        class(reconstructor_pcg), intent(in)    :: self
        complex,                  intent(in)    :: plane(self%lims2(1,1):self%lims2(1,2),&
                                                        &self%lims2(2,1):self%lims2(2,2))
        real,                     intent(in)    :: rot(3,3)
        complex,                  intent(inout) :: vol_accum(self%lims3(1,1):self%lims3(1,2),&
                                                            &self%lims3(2,1):self%lims3(2,2),&
                                                            &self%lims3(3,1):self%lims3(3,2))
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot_g(3,3)
        integer :: h, k, l, g, i0(3)
        ! replication at every R_i.S_g, g outside the colour sweep; symmats(:,:,1)=I gives c1
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
        ! serial pass over the wrapping rim (see win_wraps): confined to the outermost
        ! shell, so cheap, and reproducible unlike an atomic
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

    !> preconditioned CG solve of H x = b from every observed plane at once; kept
    !! for the test commanders. Production streams batches through begin_accum /
    !! accumulate_batch / end_accum and calls solve_accum
    subroutine solve( self, y_planes, x, maxits, rtol, rel_res_hist, niters, outcome )
        class(reconstructor_pcg),           intent(inout) :: self
        complex,                            intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                                     &self%lims2(2,1):self%lims2(2,2), *)
        real,                               intent(inout) :: x(self%box,self%box,self%box)
        integer,                  optional, intent(in)    :: maxits
        real,                     optional, intent(in)    :: rtol
        real, allocatable,        optional, intent(out)   :: rel_res_hist(:)
        integer,                  optional, intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out)   :: outcome
        if( allocated(self%b_rhs) ) deallocate(self%b_rhs)
        self%b_rhs = self%apply_adjoint_all(y_planes)
        ! b' = P b, completing the (P H P) u = P b normal equations
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        ! output-space initial guess -> CG variable (see window_div)
        call self%window_div(x)
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters, outcome)
    end subroutine solve

    !> solves against the RHS built by end_accum, no observed planes resident
    !! start_max_rel_resid: a nonzero start whose initial relative residual
    !! exceeds it is discarded for a zero start before the first iteration
    !! (free: the zero start's residual is b itself, no operator application)
    subroutine solve_accum( self, x, maxits, rtol, rel_res_hist, niters, outcome, start_max_rel_resid )
        class(reconstructor_pcg),           intent(inout) :: self
        real,                               intent(inout) :: x(self%box,self%box,self%box)
        integer,                  optional, intent(in)    :: maxits
        real,                     optional, intent(in)    :: rtol
        real, allocatable,        optional, intent(out)   :: rel_res_hist(:)
        integer,                  optional, intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out)   :: outcome
        real,                     optional, intent(in)    :: start_max_rel_resid
        if( .not. self%l_rhs ) THROW_HARD('end_accum has not been called; solve_accum')
        ! the initial guess arrives in the output space (a shipped half map is
        ! window*u): convert it to u rather than projecting again, the exit window defines the output
        call self%window_div(x)
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters, outcome, start_max_rel_resid)
    end subroutine solve_accum

    !> the solver proper. Reads the RHS from self%b_rhs: passing a component of
    !! intent(inout) self as a separate dummy is an aliasing hazard, copying it costs 67 MB at box 256
    subroutine solve_core( self, x, maxits, rtol, rel_res_hist, niters, outcome, start_max_rel_resid )
        class(reconstructor_pcg),           intent(inout) :: self
        real,                               intent(inout) :: x(self%box,self%box,self%box)
        integer,                  optional, intent(in)    :: maxits
        real,                     optional, intent(in)    :: rtol
        real, allocatable,        optional, intent(out)   :: rel_res_hist(:)
        integer,                  optional, intent(out)   :: niters
        type(pcg_solver_outcome), optional, intent(out)   :: outcome
        real,                     optional, intent(in)    :: start_max_rel_resid
        ! the reported and tested residual is the true ||r||_2/||b||_2 (the
        ! preconditioned M-norm is not monotone in PCG and wanders with a singular M;
        ! kept in the outcome as a diagnostic of how well M models H). The recurrence
        ! residual is audited against b - Hx every RESID_REPLACE iterations
        integer, parameter :: RESID_REPLACE = 25
        ! diminishing-returns stop on the relative update dx/x, which is what ends a
        ! real solve (|r|/|b| plateaus above rtol on noisy data); rtol <= 0 disables
        ! both early exits (exactly maxits iterations, needed for solver comparisons)
        real, parameter :: PCG_XTOL = 1.5e-2
        real, allocatable :: r(:,:,:), p(:,:,:), hp(:,:,:), z(:,:,:), hist(:)
        real, allocatable :: update_hist(:), mnorm_hist(:), iteration_times(:)
        real(dp)                 :: rho, rho_new, rho0, alpha, beta, pHp
        real(dp)                 :: bnorm, rnorm, xnorm, dxnorm, mnorm, dxx
        integer                  :: mmaxits, iter, n_done
        real                     :: rrtol
        logical                  :: stop_rtol, stop_xtol
        type(pcg_solver_outcome) :: result
        integer(timer_int_kind)  :: t_it
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
        ! profile the iterations only; forming the RHS is a one-off setup cost
        call self%reset_profile
        bnorm = sqrt(self%dot_real_volume(self%b_rhs,self%b_rhs))
        if( bnorm <= 0.0_dp ) THROW_HARD('zero right-hand side; nothing to reconstruct; solve')
        if( all(x == 0.0) )then
            ! zero initialization: skip the operator application known to return zero
            allocate(hp(self%box,self%box,self%box), source=0.0)
        else
            hp = self%apply_normal(x)
        endif
        r  = self%b_rhs - hp
        rnorm = sqrt(self%dot_real_volume(r,r))
        result%initial_rel_residual = real(rnorm / bnorm)
        if( present(start_max_rel_resid) )then
            if( result%initial_rel_residual > start_max_rel_resid .and. any(x /= 0.0) )then
                ! the start is worse than nothing: discard it before the first
                ! iteration; the zero start's residual is b, no operator applied
                result%start_rejected         = .true.
                result%rejected_start_initial = result%initial_rel_residual
                x     = 0.0
                r     = self%b_rhs
                rnorm = bnorm
                result%initial_rel_residual = 1.0
            endif
        endif
        z  = self%apply_precond(r)
        p  = z
        rho  = self%dot_real_volume(r,z)
        rho0 = rho
        if( rho0 <= 0.0_dp ) THROW_HARD('non-positive initial dot(r,z); preconditioner is not positive definite; solve')
        n_done = 0
        dxx = 0.0_dp
        do iter = 1, mmaxits
            t_it = pcg_tic()
            hp  = self%apply_normal(p)
            pHp = self%dot_real_volume(p,hp)
            if( .not. ieee_is_finite(pHp) .or. pHp <= 0.0_dp )then
                ! lost positive-definiteness: hand the decision to the caller with the previous iterate
                if( .not. present(outcome) )then
                    THROW_HARD('non-positive/non-finite dot(p,Hp); PCG lost positive-definiteness; solve')
                endif
                result%stop_reason       = PCG_STOP_INDEFINITE
                result%failure_curvature = real(pHp)
                result%failure_iteration = iter
                result%converged         = .false.
                exit
            endif
            alpha = rho / pHp
            x  = x + real(alpha) * p
            r  = r - real(alpha) * hp
            if( mod(iter, RESID_REPLACE) == 0 ) r = self%b_rhs - self%apply_normal(x)
            n_done  = iter
            ! headline: true relative residual; dx/x says how much the map still moves
            rnorm      = sqrt(self%dot_real_volume(r,r))
            xnorm      = sqrt(self%dot_real_volume(x,x))
            dxnorm     = abs(alpha) * sqrt(self%dot_real_volume(p,p))
            dxx        = dxnorm / max(xnorm, epsilon(1.0_dp))
            hist(iter) = real(rnorm / bnorm)
            update_hist(iter) = real(dxx)
            stop_rtol  = rrtol > 0.0 .and. rnorm / bnorm <= real(rrtol,dp)
            stop_xtol  = rrtol > 0.0 .and. dxx <= real(PCG_XTOL,dp)
            if( stop_rtol .or. stop_xtol .or. iter == mmaxits )then
                iteration_times(iter) = real(pcg_toc(t_it))
                if( stop_rtol )then
                    result%stop_reason = 'rtol'
                    result%converged   = .true.
                else if( stop_xtol )then
                    result%stop_reason = 'xtol'
                    result%converged   = .true.
                endif
                exit
            endif
            z       = self%apply_precond(r)
            rho_new = self%dot_real_volume(r,z)
            mnorm   = sqrt(abs(rho_new)/rho0)
            mnorm_hist(iter) = real(mnorm)
            iteration_times(iter) = real(pcg_toc(t_it))
            beta = rho_new / rho
            p    = z + real(beta) * p
            rho  = rho_new
        end do
        ! x = window*u: u lives on the solve domain (P H P and the projected
        ! preconditioner never leave it); the window makes the output the shipped map
        call self%window_mul(x)
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
        class(reconstructor_pcg),           intent(in) :: self
        integer,                  optional, intent(in) :: funit
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
        class(reconstructor_pcg),           intent(inout) :: self
        logical,                  optional, intent(in)    :: l_on
        self%t_setvol = 0.0_dp
        self%t_cmatcp = 0.0_dp
        self%t_ploop  = 0.0_dp
        self%t_fold   = 0.0_dp
        self%t_khat   = 0.0_dp
        self%t_prec   = 0.0_dp
        self%l_profile = .true.
        if( present(l_on) ) self%l_profile = l_on
    end subroutine reset_profile

    !> per-iteration time split: particle loop (what the kernelized operator removes)
    !! vs FFT and lattice traffic (what only a Fourier-domain solve would remove)
    subroutine report_profile( self, niters, funit )
        class(reconstructor_pcg),           intent(in) :: self
        integer,                            intent(in) :: niters
        integer,                  optional, intent(in) :: funit
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
