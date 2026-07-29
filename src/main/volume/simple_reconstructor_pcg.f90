!@descr: experimental CTF/sigma-weighted, matrix-free Fourier-projection
!  operator and preconditioned-conjugate-gradient volume solver. Isolated
!  prototype, see doc/policies/reconstruct3D_pcg_policy.md. Does not touch
!  reconstructor, reconstructor_eo, volassemble, or any production
!  reconstruction path.
!
!  Per-particle data is cached once by prep_particles (not re-derived every CG
!  iteration), the particle loop is OpenMP-parallel throughout, and the solve
!  is preconditioned with the sampling-density diagonal. An optional kernelized
!  (Toeplitz/Gram) normal operator (pcgop=kernel) makes per-iteration cost
!  independent of particle count; the matrix-free operator remains the exact
!  reference.
module simple_reconstructor_pcg
use simple_core_module_api
use simple_image, only: image
use simple_ctf,   only: ctf
implicit none

public :: reconstructor_pcg, PCG_OP_MATRIXFREE, PCG_OP_KERNEL
private
#include "simple_local_flags.inc"

integer, parameter :: PCG_OP_MATRIXFREE = 0 !< reference operator: exact, cost ~ nptcls per iteration
integer, parameter :: PCG_OP_KERNEL     = 1 !< kernelized Toeplitz operator: cost independent of nptcls

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
    real             :: lambda     = 0.0
    ! ---- per-particle inputs, cached once by prep_particles ----
    integer                       :: nptcls = 0
    logical                       :: l_use_ctf = .false.
    real,            allocatable  :: rotmats(:,:,:)  !< (3,3,nptcls)
    type(ctfparams), allocatable  :: ctfparms(:)     !< (nptcls)
    real,            allocatable  :: shifts(:,:)     !< (2,nptcls), pixels
    real,            allocatable  :: sig2(:,:)       !< (0:R,nptcls) per-particle noise power
    ! ---- lookup tables / work buffers ----
    integer, allocatable :: wrap(:)                  !< precomputed cyci_1d over the whole reachable range
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
    complex, allocatable :: acc_work(:,:,:)          !< sum_i G_i^dagger |T_i|^2, full-range
    complex, allocatable :: b_work(:,:,:)            !< RHS accumulator, full-range
    real,    allocatable :: b_rhs(:,:,:)             !< folded/deapodized/masked RHS, box^3
    logical              :: l_accum = .false.
    logical              :: l_rhs   = .false.
    ! ---- sampling-density preconditioner ----
    real,    allocatable :: precond(:,:,:)           !< 1/(rho+lambda) on wimg's cmat layout
    logical              :: l_precond = .false.
    ! ---- kernelized operator ----
    integer              :: op_mode = PCG_OP_MATRIXFREE
    real,    allocatable :: Khat(:,:,:)              !< Gram kernel on wimg's cmat layout, real
    logical              :: l_kernel = .false.
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
    logical              :: exists = .false.
  contains
    ! CONSTRUCTOR / DESTRUCTOR
    procedure :: new
    procedure :: kill
    ! SETUP
    procedure :: prep_particles
    procedure :: build_precond
    procedure :: build_kernel
    procedure :: build_operators
    ! STREAMING SETUP (batch-at-a-time; see begin_accum)
    procedure :: begin_accum
    procedure :: accumulate_batch
    procedure :: end_accum
    procedure, private :: accumulate_absT2
    procedure, private :: accumulate_rhs
    procedure, private :: precond_from_accum
    procedure, private :: fold_accum_to_khat
    procedure, private :: finalize_khat
    procedure :: set_deapod
    procedure :: set_mask
    procedure, private :: build_env
    procedure, private :: deapod_mul
    procedure, private :: mask_mul
    procedure, private :: calibrate_kernel
    procedure :: measure_kernel_scale
    procedure :: set_op_mode
    ! LOW-LEVEL OPERATOR (public: the test commanders drive these directly to
    ! verify the adjoint identity before any solve() result may be trusted)
    procedure :: set_volume
    procedure :: forward_plane
    procedure :: fourier_dot
    procedure :: adjoint_plane_add
    procedure :: build_transfer
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
    ! SOLVER
    procedure :: solve
    procedure :: solve_accum
    procedure, private :: solve_core
    ! PROFILING
    procedure :: reset_profile
    procedure :: report_profile
    ! PRIVATE HELPERS
    procedure, private :: absT2_plane
    procedure, private :: transfer_plane_cmplx
    procedure, private :: fold_and_ifft
    procedure, private :: apply_precond
    procedure, private :: ensure_wimg
    procedure, private :: pad_vol
    procedure, private :: crop_vol
end type reconstructor_pcg

contains

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
        call self%build_env
        self%exists = .true.
    end subroutine new

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
    !!          Normalized by instr(0) because apod_mat_3d normalizes the KB
    !!          weights to sum 1, which makes the envelope unity at the centre.
    !!          The envelope is MEASURED, not derived. `kbinterpol%instr` is the
    !!          CONTINUOUS Fourier transform of the KB kernel, but this operator
    !!          interpolates with three discrete, renormalized weight samples --
    !!          the two disagree by about a factor of two at the box edge, which
    !!          is far too much to guess at. Instead the impulse response is
    !!          obtained by scattering a unit value through the operator's own
    !!          KB stencil and running it back through fold_and_ifft, so every
    !!          packing, Nyquist and normalization convention matches the real
    !!          operator by construction rather than by derivation.
    subroutine build_env( self )
        class(reconstructor_pcg), intent(inout) :: self
        real, parameter :: EPS_DIV = 1.0e-8
        complex, allocatable :: vol_accum(:,:,:)
        real    :: w(self%wdim,self%wdim,self%wdim), ctrval
        integer :: i0(3), c, i, j, k
        call self%ensure_wimg
        if( allocated(self%env)    ) deallocate(self%env)
        if( allocated(self%invenv) ) deallocate(self%invenv)
        allocate(vol_accum(self%lims3(1,1):self%lims3(1,2),&
                          &self%lims3(2,1):self%lims3(2,2),&
                          &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        ! deposit the KB stencil at the Fourier origin; its inverse transform
        ! IS the real-space envelope the gather imposes
        call self%kbwin%apod_mat_3d_fast([0.,0.,0.], self%iwinsz, self%wdim, w)
        i0 = -self%iwinsz
        call scatter_window(self, i0, w, cmplx(1.,0.), vol_accum)
        ! fold_and_ifft already crops to the native box, so this measures the
        ! envelope exactly where it is applied. With oversampling it now spans
        ! only +/-(box/2)/boxpd of the padded Nyquist, so the roll-off is mild
        ! (~1.5x per axis) instead of the ~5x per axis of the native-lattice
        ! operator -- which is what made the old deapodization so violent.
        self%env = self%fold_and_ifft(vol_accum)
        deallocate(vol_accum)
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
        if( allocated(self%wrap)     ) deallocate(self%wrap)
        if( allocated(self%env)      ) deallocate(self%env)
        if( allocated(self%invenv)   ) deallocate(self%invenv)
        if( allocated(self%mask)     ) deallocate(self%mask)
        if( allocated(self%precond)  ) deallocate(self%precond)
        if( allocated(self%Khat)     ) deallocate(self%Khat)
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
        self%nptcls = 0
        self%l_use_ctf   = .false.
        self%l_mask      = .false.
        self%l_precond   = .false.
        self%l_kernel    = .false.
        self%wimg_exists = .false.
        self%op_mode     = PCG_OP_MATRIXFREE
        call self%reset_profile(.false.)
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
        type(ctf) :: tfun
        real      :: spafreqsq, ang, cval, arg, sw
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
        endif
        T = cmplx(0.,0.)
        !$omp parallel do collapse(2) default(shared) &
        !$omp private(h,k,spafreqsq,ang,cval,arg,shell,sw) schedule(static) proc_bind(close)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp ) cycle
                cval = 1.0
                if( l_ctf )then
                    spafreqsq = (real(h)/real(self%box))**2 + (real(k)/real(self%box))**2
                    ang       = atan2(real(k), real(h))
                    cval      = tfun%eval(spafreqsq, ang, ctfparms%phshift)
                    if( l_flip ) cval = abs(cval)
                endif
                ! SIGN. apply_adjoint_all forms conjg(T)*y, and production applies
                ! exp(i*(-shift)*shconst*(h,k)) to y (gen_fplane4rec: pshift =
                ! -shift*shconst). So T must carry the POSITIVE phase for its
                ! conjugate to reproduce production's correction. A negated phase
                ! here would shift every particle by twice its own displacement in
                ! the wrong direction -- invisible on synthetic
                ! data, where the same build_transfer generates the observations.
                arg       = 2.0*PI * (real(h)*shift(1) + real(k)*shift(2)) / real(self%box)
                sw        = 1.0
                if( present(sig2arr) )then
                    shell = min(nint(sqrt(real(h*h+k*k))), ubound(sig2arr,1))
                    sw    = 1.0 / sqrt(sig2arr(shell))
                endif
                T(h,k) = cval * cmplx(cos(arg), sin(arg)) * sw
            end do
        end do
        !$omp end parallel do
    end function build_transfer

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
        integer :: i, h, k, l, i0(3)
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
        !$omp parallel default(shared) private(i,h,k,l,loc,i0,w,comp,rot) proc_bind(close)
        do i = 1, self%nptcls
            rot = self%rotmats(:,:,i)
            ! per-particle real weight |T_i|^2, shared across threads. The
            ! worksharing inside absT2_plane is orphaned and binds to THIS
            ! region, so the CTF evaluation is spread across the team; its
            ! trailing barrier is what makes absT2 safe to read below.
            call self%absT2_plane(i, absT2)
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
        !$omp end parallel
        if( self%l_profile ) self%t_ploop = self%t_ploop + real(toc(tp),dp)
        hp = self%fold_and_ifft(vol_accum)
        call self%deapod_mul(hp)
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
        !$omp parallel do collapse(3) default(shared) private(i,j,k) &
        !$omp schedule(static) proc_bind(close)
        do k = 1, cdim(3)
            do j = 1, cdim(2)
                do i = 1, cdim(1)
                    cmat(i,j,k) = cmat(i,j,k) * self%Khat(i,j,k)
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
        hp = hp + self%lambda * p
    end function apply_normal_kernel

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
    ! Peak memory during accumulation is the two full-range accumulators
    ! (2 x 1.08 GB at box 256) plus the work image; that floor is constant in
    ! nptcls, which is the whole point. See
    ! doc/implementation_notes/pcg_reconstruction_as_gridding_replacement.md 5.1.

    subroutine begin_accum( self )
        class(reconstructor_pcg), intent(inout) :: self
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; begin_accum')
        call self%ensure_wimg
        if( allocated(self%acc_work) ) deallocate(self%acc_work)
        if( allocated(self%b_work)   ) deallocate(self%b_work)
        if( allocated(self%b_rhs)    ) deallocate(self%b_rhs)
        allocate(self%acc_work(self%lims3(1,1):self%lims3(1,2),&
                              &self%lims3(2,1):self%lims3(2,2),&
                              &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        allocate(self%b_work(self%lims3(1,1):self%lims3(1,2),&
                            &self%lims3(2,1):self%lims3(2,2),&
                            &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        self%l_accum = .true.
        self%l_rhs   = .false.
    end subroutine begin_accum

    !>  \brief  Accumulates one batch into BOTH accumulators. y_batch holds the
    !!          observed planes for particles ifrom .. ifrom+nb-1 of the
    !!          selection passed to prep_particles.
    subroutine accumulate_batch( self, y_batch, nb, ifrom )
        class(reconstructor_pcg), intent(inout) :: self
        integer,                    intent(in)    :: nb, ifrom
        complex,                    intent(in)    :: y_batch(self%lims2(1,1):self%lims2(1,2),&
                                                              &self%lims2(2,1):self%lims2(2,2), nb)
        complex, allocatable :: work(:,:,:)
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
        call move_alloc(self%b_work, work)
        call self%accumulate_rhs(y_batch, nb, ifrom, work)
        call move_alloc(work, self%b_work)
        call move_alloc(self%acc_work, work)
        call self%accumulate_absT2(work, ifrom, ifrom + nb - 1)
        call move_alloc(work, self%acc_work)
    end subroutine accumulate_batch

    !>  \brief  Closes accumulation: folds the RHS once, derives preconditioner
    !!          and (optionally) the Gram kernel, and releases the accumulators.
    subroutine end_accum( self, l_kernel )
        class(reconstructor_pcg), intent(inout) :: self
        logical,                    intent(in)    :: l_kernel
        complex, allocatable :: work(:,:,:)
        if( .not. self%l_accum ) THROW_HARD('begin_accum has not been called; end_accum')
        ! RHS first, and free its accumulator before the kernel work starts:
        ! finalize_khat is the allocation-heavy half and should not overlap it.
        ! move_alloc for the same aliasing reason as accumulate_batch.
        call move_alloc(self%b_work, work)
        self%b_rhs = self%fold_and_ifft(work)
        deallocate(work)
        call self%deapod_mul(self%b_rhs)
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        call move_alloc(self%acc_work, work)
        call self%precond_from_accum(work)
        if( l_kernel ) call self%fold_accum_to_khat(work)
        deallocate(work)
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
        complex, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_precond')
        call self%ensure_wimg
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
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
        complex,                    intent(inout) :: acc(self%lims3(1,1):self%lims3(1,2),&
                                                          &self%lims3(2,1):self%lims3(2,2),&
                                                          &self%lims3(3,1):self%lims3(3,2))
        integer,          optional, intent(in)    :: ifrom, ito
        real,    allocatable :: absT2(:,:)
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim), rot(3,3)
        integer :: i, h, k, l, i0(3), ii_from, ii_to
        ii_from = 1
        ii_to   = self%nptcls
        if( present(ifrom) ) ii_from = ifrom
        if( present(ito)   ) ii_to   = ito
        if( ii_to < ii_from ) return
        allocate(absT2(self%lims2(1,1):self%lims2(1,2), self%lims2(2,1):self%lims2(2,2)))
        !$omp parallel default(shared) private(i,h,k,l,loc,i0,w,rot) proc_bind(close)
        do i = ii_from, ii_to
            rot = self%rotmats(:,:,i)
            call self%absT2_plane(i, absT2)
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
                        call scatter_window(self, i0, w, cmplx(absT2(h,k),0.), acc)
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
                    call scatter_window(self, i0, w, cmplx(absT2(h,k),0.), acc)
                end do
            end do
            !$omp end single
        end do
        !$omp end parallel
    end subroutine accumulate_absT2

    subroutine precond_from_accum( self, rho_accum )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: rho_accum(self%lims3(1,1):self%lims3(1,2),&
                                                                &self%lims3(2,1):self%lims3(2,2),&
                                                                &self%lims3(3,1):self%lims3(3,2))
        real, parameter :: RHO_FLOOR_FRAC = 1.0e-2
        real(dp), allocatable :: shsum(:)
        integer,  allocatable :: shcnt(:)
        real,     allocatable :: shfloor(:)
        integer :: h, k, cdim(3), m, hh, phys(3), sh, sh_lim
        real    :: denom, rval
        ! fold onto the packed cmat layout and invert with a guard
        cdim = self%wimg%get_array_shape()
        if( allocated(self%precond) ) deallocate(self%precond)
        allocate(self%precond(cdim(1),cdim(2),cdim(3)), source=0.0)
        ! Data only ever reaches |loc| <= padf*Rnat, so beyond that radius rho is
        ! identically zero: those modes are completely unconstrained. Zero them,
        ! exactly as reconstructor%sampl_dens_correct does for sh > sh_lim. A
        ! singular (PSD) preconditioner is the right tool here -- it keeps the
        ! Krylov space out of the null space entirely, so starting from x = 0
        ! those modes are never excited at all.
        sh_lim = self%padf * self%Rnat
        ! PASS 1: mean rho per shell, over the sampled voxels only. Averaging in
        ! the zeros would drag the floor down exactly where sampling is sparse,
        ! which is the regime the floor exists to protect.
        allocate(shsum(0:sh_lim), source=0.0_dp)
        allocate(shcnt(0:sh_lim), source=0)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > sh_lim ) cycle
                    rval = real(rho_accum(self%wrap(h),k,m))
                    if( rval <= 0.0 ) cycle
                    shsum(sh) = shsum(sh) + real(rval,dp)
                    shcnt(sh) = shcnt(sh) + 1
                end do
            end do
        end do
        allocate(shfloor(0:sh_lim), source=0.0)
        do sh = 0, sh_lim
            if( shcnt(sh) > 0 ) shfloor(sh) = RHO_FLOOR_FRAC * real(shsum(sh) / real(shcnt(sh),dp))
        end do
        ! PASS 2: guarded reciprocal. An empty shell (shfloor == 0) leaves its
        ! voxels at the initialized 0.0, i.e. treated as unconstrained.
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    sh = nint(sqrt(real(h*h + k*k + m*m)))
                    if( sh > sh_lim ) cycle
                    if( shfloor(sh) <= 0.0 ) cycle
                    hh    = self%wrap(h)
                    denom = max(real(rho_accum(hh,k,m)), 0.0) + shfloor(sh)
                    phys  = self%wimg%comp_addr_phys(h,k,m)
                    self%precond(phys(1),phys(2),phys(3)) = 1.0 / denom
                end do
            end do
        end do
        self%l_precond = .true.
    end subroutine precond_from_accum

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
        complex, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_kernel')
        call self%ensure_wimg
        ! With oversampling in place the operator ALREADY works on a padf*box
        ! lattice, which is exactly the 2x grid the Gram kernel needs for a
        ! linear (non-wrapping) convolution of a box-supported volume. The two
        ! lattices coincide, so no separate padded image is required.
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
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
        complex, allocatable :: acc(:,:,:)
        if( self%nptcls < 1 ) THROW_HARD('prep_particles has not been called; build_operators')
        call self%ensure_wimg
        allocate(acc(self%lims3(1,1):self%lims3(1,2),&
                    &self%lims3(2,1):self%lims3(2,2),&
                    &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        call self%accumulate_absT2(acc)
        call self%precond_from_accum(acc)
        if( l_kernel ) call self%fold_accum_to_khat(acc)
        deallocate(acc)
        if( l_kernel ) call self%finalize_khat
    end subroutine build_operators

    subroutine fold_accum_to_khat( self, kacc )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: kacc(self%lims3(1,1):self%lims3(1,2),&
                                                           &self%lims3(2,1):self%lims3(2,2),&
                                                           &self%lims3(3,1):self%lims3(3,2))
        integer :: h, k, m, hh, phys(3), cdim(3)
        cdim = self%wimg%get_array_shape()
        if( allocated(self%Khat) ) deallocate(self%Khat)
        allocate(self%Khat(cdim(1),cdim(2),cdim(3)), source=0.0)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    self%Khat(phys(1),phys(2),phys(3)) = real(kacc(hh,k,m))
                end do
            end do
        end do
    end subroutine fold_accum_to_khat

    subroutine finalize_khat( self )
        class(reconstructor_pcg), intent(inout) :: self
        real,    parameter   :: EPS_D = 1.0e-8
        complex, allocatable :: dacc(:,:,:), ctmp(:,:,:)
        real,    allocatable :: tker(:,:,:), dep(:,:,:)
        real    :: wz(self%wdim,self%wdim,self%wdim), depc
        integer :: h, k, i0z(3), cdim(3), m, phys(3), hh, cc
        cdim = self%wimg%get_array_shape()
        ! ---- divide out the DEPOSITION envelope ----
        ! Laying |T|^2 down through the KB window convolves the kernel spectrum
        ! with that window, which multiplies the real-space kernel by the
        ! window's transform. That is a SECOND envelope, distinct from the
        ! gather's, and left in it shows up as a spurious spatial taper. It is
        ! measured rather than derived, for the same reason build_env measures
        ! its own: deposit a unit spike at the origin through the operator's own
        ! stencil and transform back, so conventions match by construction.
        allocate(dacc(self%lims3(1,1):self%lims3(1,2),&
                     &self%lims3(2,1):self%lims3(2,2),&
                     &self%lims3(3,1):self%lims3(3,2)), source=cmplx(0.,0.))
        call self%kbwin%apod_mat_3d_fast([0.,0.,0.], self%iwinsz, self%wdim, wz)
        i0z = -self%iwinsz
        call scatter_window(self, i0z, wz, cmplx(1.,0.), dacc)
        allocate(ctmp(cdim(1),cdim(2),cdim(3)), source=cmplx(0.,0.))
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, self%lims3(1,2)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    ctmp(phys(1),phys(2),phys(3)) = dacc(hh,k,m)
                end do
            end do
        end do
        deallocate(dacc)
        call self%wimg%set_cmat(ctmp)
        call self%wimg%ifft()
        dep  = self%wimg%get_rmat()
        cc   = self%boxpd/2 + 1
        depc = dep(cc,cc,cc)
        if( abs(depc) < EPS_D ) depc = 1.0
        dep = dep / depc
        ctmp = cmplx(self%Khat, 0.)
        call self%wimg%set_cmat(ctmp)
        call self%wimg%ifft()
        tker = self%wimg%get_rmat()
        where( abs(dep) > EPS_D )
            tker = tker / dep
        elsewhere
            tker = 0.
        end where
        call self%wimg%set_rmat(tker, .false.)
        call self%wimg%fft()
        ctmp      = self%wimg%get_cmat()
        self%Khat = real(ctmp)
        deallocate(ctmp, tker, dep)
        self%l_kernel = .true.
        call self%calibrate_kernel
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
    !!          doc/implementation_notes/pcg_reconstruction_as_gridding_replacement.md
    !!          section 1.1.
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
        if( .not. self%l_kernel ) THROW_HARD('build_kernel has not been called; measure_kernel_scale')
        lam_save    = self%lambda
        self%lambda = 0.0   ! compare the DATA term only
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
        type(ctf) :: tfun
        real      :: spafreqsq, ang, cval
        integer   :: h, k, shell, R
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
        R = self%lims2(1,2)
        ! same ctfflag semantics as build_transfer; note |CTF|^2 == CTF^2, so the
        ! FLIP case only matters for the complex T, not here -- but the NO case
        ! does, and applying a CTF the images never saw would be wrong.
        l_ctf  = self%ctfparms(iptcl)%ctfflag /= CTFFLAG_NO
        l_flip = self%ctfparms(iptcl)%ctfflag == CTFFLAG_FLIP
        if( l_ctf )then
            tfun = ctf(self%ctfparms(iptcl)%smpd, self%ctfparms(iptcl)%kv, &
                &self%ctfparms(iptcl)%cs, self%ctfparms(iptcl)%fraca)
            call tfun%init(self%ctfparms(iptcl)%dfx, self%ctfparms(iptcl)%dfy, self%ctfparms(iptcl)%angast)
        endif
        !$omp do collapse(2) schedule(static) private(spafreqsq,ang,cval,shell)
        do k = self%lims2(2,1), self%lims2(2,2)
            do h = self%lims2(1,1), self%lims2(1,2)
                if( h*h + k*k > self%sqlp )then
                    absT2(h,k) = 0.
                else
                    cval = 1.0
                    if( l_ctf )then
                        spafreqsq = (real(h)/real(self%box))**2 + (real(k)/real(self%box))**2
                        ang       = atan2(real(k), real(h))
                        cval      = tfun%eval(spafreqsq, ang, self%ctfparms(iptcl)%phshift)
                        if( l_flip ) cval = abs(cval)
                    endif
                    shell      = min(nint(sqrt(real(h*h+k*k))), R)
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
        real    :: loc(3), w(self%wdim,self%wdim,self%wdim)
        integer :: h, k, l, i0(3)
        !$omp parallel default(shared) private(h,k,l,loc,i0,w) proc_bind(close)
        do l = 0, self%stride-1
            !$omp do schedule(static,1)
            do h = self%lims2(1,1)+l, self%lims2(1,2), self%stride
                do k = self%lims2(2,1), self%lims2(2,2)
                    if( h*h + k*k > self%sqlp ) cycle
                    if( plane(h,k) == cmplx(0.,0.) ) cycle
                    loc = real(self%padf) * matmul(real([h,k,0]), rot)
                    i0  = nint(loc) - self%iwinsz
                    if( win_wraps(self, i0) ) cycle   ! rim: deferred to the serial pass
                    call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                    call scatter_window(self, i0, w, self%padsc * plane(h,k), vol_accum)
                end do
            end do
            !$omp end do
        end do
        !$omp end parallel
        ! Serial pass over the wrapping rim, where the colouring's separation
        ! guarantee does not survive folding (see win_wraps). Confined to the
        ! outermost shell, so the cost is negligible, and being serial it is also
        ! reproducible -- which an atomic would not be, since the summation order
        ! would still vary.
        do h = self%lims2(1,1), self%lims2(1,2)
            do k = self%lims2(2,1), self%lims2(2,2)
                if( h*h + k*k > self%sqlp   ) cycle
                if( h*h + k*k <= self%sq_rim ) cycle       ! provably cannot wrap
                if( plane(h,k) == cmplx(0.,0.) ) cycle
                loc = real(self%padf) * matmul(real([h,k,0]), rot)
                i0  = nint(loc) - self%iwinsz
                if( .not. win_wraps(self, i0) ) cycle
                call self%kbwin%apod_mat_3d_fast(loc, self%iwinsz, self%wdim, w)
                call scatter_window(self, i0, w, self%padsc * plane(h,k), vol_accum)
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
    subroutine solve( self, y_planes, x, maxits, rtol, rel_res_hist, niters )
        class(reconstructor_pcg), intent(inout) :: self
        complex,                    intent(in)    :: y_planes(self%lims2(1,1):self%lims2(1,2),&
                                                               &self%lims2(2,1):self%lims2(2,2), *)
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        if( allocated(self%b_rhs) ) deallocate(self%b_rhs)
        self%b_rhs = self%apply_adjoint_all(y_planes)
        ! b' = P b, completing the (P H P) u = P b normal equations
        call self%mask_mul(self%b_rhs)
        self%l_rhs = .true.
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters)
    end subroutine solve

    !>  \brief  Solve against the RHS built by end_accum. Same solver, no
    !!          observed planes resident.
    subroutine solve_accum( self, x, maxits, rtol, rel_res_hist, niters )
        class(reconstructor_pcg), intent(inout) :: self
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        if( .not. self%l_rhs ) THROW_HARD('end_accum has not been called; solve_accum')
        call self%solve_core(x, maxits, rtol, rel_res_hist, niters)
    end subroutine solve_accum

    !>  \brief  The solver proper. Reads the RHS from self%b_rhs rather than
    !!          taking it as an argument: passing a component of `self` as a
    !!          separate dummy while `self` is intent(inout) is an aliasing
    !!          hazard the standard lets a compiler exploit, and the alternative
    !!          (copying it) would cost 67 MB per solve at box 256 for nothing.
    subroutine solve_core( self, x, maxits, rtol, rel_res_hist, niters )
        class(reconstructor_pcg), intent(inout) :: self
        real,                        intent(inout) :: x(self%box,self%box,self%box)
        integer,          optional,  intent(in)    :: maxits
        real,             optional,  intent(in)    :: rtol
        real, allocatable, optional, intent(out)   :: rel_res_hist(:)
        integer,          optional,  intent(out)   :: niters
        ! WHICH RESIDUAL IS REPORTED MATTERS. This used to report
        ! ||r||_Minv / ||r0||_Minv, the PRECONDITIONED norm, which in PCG is not
        ! monotone -- only ||e||_H is -- and with a singular preconditioner it
        ! wanders freely. It duly oscillated, which reads exactly like a solver
        ! that has lost conjugacy while the solve was in fact converging
        ! perfectly well. The headline number and the stopping test are now the
        ! true relative residual ||r||_2 / ||b||_2, which is what a reader
        ! actually wants and costs two dot products against 1.5 s of FFTs. The
        ! M-norm is still logged, because it is the honest diagnostic for the
        ! PRECONDITIONER: a large gap between the two says M is a poor model of H.
        !
        ! RESIDUAL REPLACEMENT. CG propagates the residual by the recurrence
        ! r <- r - alpha*Hp rather than recomputing b - Hx, which is what makes
        ! an iteration cost one operator application instead of two; the two can
        ! drift apart in finite precision. Measured here, they agree to six
        ! significant figures at box 256, so this is kept only as a periodic
        ! audit at a long interval rather than a correction -- reporting ||r||_2
        ! from the recurrence is only legitimate because this check backs it up.
        integer, parameter :: RESID_REPLACE = 25
        real, allocatable :: r(:,:,:), p(:,:,:), hp(:,:,:), z(:,:,:), hist(:)
        real, allocatable :: rtr(:,:,:)
        real(dp) :: rho, rho_new, rho0, alpha, beta, pHp, nrec, ntru
        real(dp) :: bnorm, rnorm, xnorm, dxnorm, mnorm
        integer  :: mmaxits, iter, n_done
        real     :: rrtol
        integer(timer_int_kind) :: t_it
        real(timer_int_kind)    :: rt_it
        mmaxits = 50
        if( present(maxits) ) mmaxits = maxits
        rrtol = 1.0e-4
        if( present(rtol) ) rrtol = rtol
        allocate(hist(mmaxits))
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
        n_done = 0
        write(logfhandle,'(a,i0,a,es12.5)') '>>> PCG: starting, maxits = ', mmaxits, ', rtol = ', rrtol
        call flush(logfhandle)
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
                nrec = sqrt(self%dot_real_volume(r,r))
                rtr  = self%b_rhs - self%apply_normal(x)
                ntru = sqrt(self%dot_real_volume(rtr,rtr))
                write(logfhandle,'(a,i4,a,es12.5,a,es12.5)') '>>> PCG residual replacement at iter ', &
                    &iter, ':  recurrence = ', real(nrec), '   true = ', real(ntru)
                call flush(logfhandle)
                r = rtr
            endif
            z  = self%apply_precond(r)
            rho_new = self%dot_real_volume(r,z)
            n_done  = iter
            ! headline: true relative residual. dx/x says how much the map is
            ! still moving, which for a reconstruction is often the more
            ! practical stopping signal than any residual level.
            rnorm      = sqrt(self%dot_real_volume(r,r))
            xnorm      = sqrt(self%dot_real_volume(x,x))
            dxnorm     = abs(alpha) * sqrt(self%dot_real_volume(p,p))
            mnorm      = sqrt(abs(rho_new)/rho0)
            hist(iter) = real(rnorm / bnorm)
            rt_it      = toc(t_it)
            ! Per-iteration progress. Without this the solver is silent between
            ! setup and completion, which on a real data set is indistinguishable
            ! from a hang.
            write(logfhandle,'(a,i4,a,es11.4,a,es10.3,a,es10.3,a,f7.2,a)') '>>> PCG iter ', iter, &
                &'  |r|/|b| = ', hist(iter), '  dx/x = ', real(dxnorm/max(xnorm,epsilon(1.0_dp))), &
                &'  (M-norm ', real(mnorm), ')  (', rt_it, ' s)'
            call flush(logfhandle)
            if( rnorm / bnorm <= real(rrtol,dp) ) exit
            beta = rho_new / rho
            p    = z + real(beta) * p
            rho  = rho_new
        end do
        ! x = P u: the CG variable u is unconstrained outside the support (P H P
        ! annihilates it there, so it never influenced the residual), but it does
        ! accumulate arbitrary values via the preconditioner. This is the step
        ! that makes the returned volume the constrained solution.
        call self%mask_mul(x)
        if( present(niters) ) niters = n_done
        if( present(rel_res_hist) ) allocate(rel_res_hist(n_done), source=hist(1:n_done))
        self%l_profile = .false.
        if( n_done > 0 ) call self%report_profile(n_done)
    end subroutine solve_core

    ! PROFILING

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
    subroutine report_profile( self, niters )
        class(reconstructor_pcg), intent(in) :: self
        integer,                    intent(in) :: niters
        real(dp) :: rn, tot, ffts
        if( niters < 1 ) return
        rn   = real(niters,dp)
        ffts = self%t_setvol + self%t_cmatcp + self%t_fold + self%t_prec
        tot  = ffts + self%t_ploop + self%t_khat
        write(logfhandle,'(a)')    '>>> PCG PROFILE (seconds per iteration)'
        write(logfhandle,'(a,f9.3)') '    pad + fwd FFT of iterate     : ', self%t_setvol / rn
        write(logfhandle,'(a,f9.3)') '    get_cmat/set_cmat copies     : ', self%t_cmatcp / rn
        write(logfhandle,'(a,f9.3)') '    particle loop                : ', self%t_ploop  / rn
        write(logfhandle,'(a,f9.3)') '    kernel pointwise multiply    : ', self%t_khat   / rn
        write(logfhandle,'(a,f9.3)') '    fold + inv FFT + crop        : ', self%t_fold   / rn
        write(logfhandle,'(a,f9.3)') '    apply_precond (2 FFT + copy) : ', self%t_prec   / rn
        write(logfhandle,'(a,f9.3)') '    ---- accounted subtotal      : ', tot / rn
        if( tot > 0.0_dp )then
            write(logfhandle,'(a,f7.1,a)') '    particle loop is ', &
                &100.0_dp * self%t_ploop / tot, '% of accounted time'
            write(logfhandle,'(a,f7.1,a)') '    FFT + lattice traffic is ', &
                &100.0_dp * ffts / tot, '% of accounted time'
        endif
        call flush(logfhandle)
    end subroutine report_profile

end module simple_reconstructor_pcg
