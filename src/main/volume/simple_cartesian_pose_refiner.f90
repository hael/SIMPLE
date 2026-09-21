!@descr: Fixed-volume Cartesian Fourier particle-pose refinement numerics
module simple_cartesian_pose_refiner
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
use simple_core_module_api
use simple_image, only: image
use simple_ctf, only: ctf
use simple_cartesian_fourier, only: center_embed_real3d, gather_packed_window_grad
use simple_gridding, only: kb_stencil_centered_crop_inv_envelope_1d
implicit none
private

#include "simple_local_flags.inc"

! Types
public :: cartesian_pose_refiner, cartesian_pose_data, shift_lm_config, pose_lm_config
public :: pose_lm_result, pose_lm_diagnostics

! Cartesian objective identifiers. These are deliberately distinct from
! SIMPLE's polar/PFTC OBJFUN_* contract.
public :: POSE_CONT_OBJECTIVE_CART_NCC, POSE_CONT_OBJECTIVE_CART_EUCLID

! LM status codes
public :: LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT
public :: LM_NO_RELIABLE_UPDATE, LM_STEP_BOUND_REJECTED
public :: LM_INVALID_NUMERICS, LM_ITERATION_LIMIT

! Public procedures
public :: right_increment_rotation

! LM result codes shared by the shift-only and five-parameter solvers
integer, parameter :: LM_ACCEPTED_IMPROVEMENT = 1
integer, parameter :: LM_FINITE_NO_IMPROVEMENT = 2
integer, parameter :: LM_NO_RELIABLE_UPDATE = 3
integer, parameter :: LM_STEP_BOUND_REJECTED = 4
integer, parameter :: LM_INVALID_NUMERICS = 5
integer, parameter :: LM_ITERATION_LIMIT = 6

integer, parameter :: POSE_CONT_OBJECTIVE_CART_NCC = 0
integer, parameter :: POSE_CONT_OBJECTIVE_CART_EUCLID = 1

! Internal numerical constants
real(dp), parameter :: POSE_NUMERIC_FLOOR = epsilon(1._dp)**2
real(dp), parameter :: LM_INITIAL_DAMPING = 1.e-3_dp
real(dp), parameter :: LM_INITIAL_REJECTION_MULTIPLIER = 4._dp
real(dp), parameter :: LM_MAX_DAMPING = 1._dp/epsilon(1._dp)

!> One shift-free, noise-whitened particle observation for a fixed reference.
type :: cartesian_pose_data
    private
    complex, allocatable :: observed(:, :), transfer(:, :)
    real, allocatable :: sigma2(:)
    integer :: shell_range(2) = 0
    logical :: valid = .false.
contains
    procedure :: is_valid => pose_data_is_valid
    procedure :: get_shell_range => pose_data_get_shell_range
end type cartesian_pose_data

type :: cartesian_pose_refiner
    private
    integer          :: box     = 0
    integer          :: boxpd   = 0
    integer          :: padf    = 1
    integer          :: iwinsz  = 0
    integer          :: wdim    = 0
    integer          :: lims2(2,2) = 0
    real             :: padsc   = 1.0
    type(kbinterpol) :: kbwin
    integer, allocatable :: wrap(:)
    complex, allocatable :: cmat(:, :, :)
    logical :: exists = .false.
contains
    procedure :: new_inverse_envelope_reference => new_pose_refiner_inverse_envelope_reference
    procedure :: new_physical_reference => new_pose_refiner_physical_reference
    procedure :: prepare_particle => prepare_pose_particle
    procedure :: prepared_objective_gradient => prepared_pose_objective_gradient
    procedure :: prepared_sigma_contribution => prepared_pose_sigma_contribution
    procedure :: refine_prepared_pose_lm
    procedure :: kill => kill_fourier_workspace
    procedure :: predict_unweighted => predict_unweighted_pose
    ! Retain the efficient dedicated two-parameter LM path.
    procedure :: refine_shift_lm
    procedure, private :: sample_with_grad => sample_fourier_with_grad
    procedure, private :: pose_objective_gradient
    procedure, private :: refine_pose_lm
    procedure, private :: shift_normal_terms
    procedure, private :: pose_normal_terms
end type cartesian_pose_refiner

!> Configuration and bounds for one dedicated two-parameter shift solve.
type :: shift_lm_config
    real(dp) :: shift_step_bound = 1._dp
    integer  :: max_iterations = 40
    integer  :: objective = POSE_CONT_OBJECTIVE_CART_EUCLID
end type shift_lm_config

!> Configuration and bounds for one five-parameter LM solve.
!! Particle observations and transfer weights remain owned by cartesian_pose_data.
type :: pose_lm_config
    real(dp) :: rotation_scale = 1._dp
    real(dp) :: shift_step_bound = 1._dp
    integer  :: max_iterations = 40
    integer  :: objective = POSE_CONT_OBJECTIVE_CART_EUCLID
    logical  :: active_parameters(5) = .true.
    logical  :: use_cumulative_guard = .false.
    real(dp) :: anchor_rotmat(3, 3) = reshape([1._dp, 0._dp, 0._dp, &
        &0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 1._dp], [3, 3])
    real(dp) :: anchor_shift(2) = 0._dp
    real(dp) :: max_total_rotation = 0._dp
    real(dp) :: max_total_shift = 0._dp
end type pose_lm_config

!> Minimal production result from one LM solve.
type :: pose_lm_result
    integer :: status = LM_ITERATION_LIMIT
    integer :: niterations = 0
end type pose_lm_result

!> Optional scalar diagnostics from one LM solve.
type :: pose_lm_diagnostics
    integer :: nattempted = 0
    integer :: naccepted = 0
    integer :: nbound_hits = 0
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step = 0._dp
contains
    procedure :: reset => reset_pose_lm_diagnostics
end type pose_lm_diagnostics

contains

    subroutine reset_pose_lm_diagnostics(self)
        class(pose_lm_diagnostics), intent(inout) :: self

        self%nattempted = 0
        self%naccepted = 0
        self%nbound_hits = 0
        self%max_rotation_step = 0._dp
        self%max_shift_step = 0._dp
    end subroutine reset_pose_lm_diagnostics

    pure logical function pose_data_is_valid(self) result(valid)
        class(cartesian_pose_data), intent(in) :: self
        valid = self%valid
    end function pose_data_is_valid

    !> requested range ∩ available sigma2 shells ∩ Cartesian Nyquist limit
    pure function pose_data_get_shell_range(self) result(shell_range)
        class(cartesian_pose_data), intent(in) :: self
        integer :: shell_range(2)
        shell_range = self%shell_range
    end function pose_data_get_shell_range

    !> Evaluate SIMPLE's CTF object at one signed full-disk Fourier coordinate.
    !! Unlike the memoized hot-loop kernel, this route owns no process-global
    !! Fourier maps and is therefore valid in standalone adapter tests.
    real function pose_cont_ctf_value(tfun, h, k, box, phshift, phase_flip) result(cval)
        type(ctf), intent(in) :: tfun
        integer, intent(in) :: h, k, box
        real, intent(in) :: phshift
        logical, intent(in) :: phase_flip
        real :: angle, spatial_frequency_squared

        spatial_frequency_squared = (real(h)*real(h) + real(k)*real(k))/(real(box)*real(box))
        angle = 0.
        if (h /= 0 .or. k /= 0) angle = atan2(real(k), real(h))
        cval = tfun%eval_canonical(spatial_frequency_squared, angle, phshift)
        if (phase_flip) cval = abs(cval)
    end function pose_cont_ctf_value

    !> Construct an immutable Cartesian Fourier reference from a physical volume.
    !! Apply the inverse Kaiser-Bessel envelope exactly once before padding and FFT.
    subroutine new_pose_refiner_inverse_envelope_reference(self, volume)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:, :, :)
        call load_pose_reference(self, volume, .true.)
    end subroutine new_pose_refiner_inverse_envelope_reference

    !> Construct an immutable Cartesian Fourier reference from a prepared physical volume.
    !! Preserve the supplied amplitudes: pad and FFT without applying an inverse
    !! Kaiser-Bessel envelope. This matches the executed refine3D reference boundary.
    subroutine new_pose_refiner_physical_reference(self, volume)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:, :, :)
        call load_pose_reference(self, volume, .false.)
    end subroutine new_pose_refiner_physical_reference

    subroutine load_pose_reference(self, volume, apply_inverse_envelope)
        class(cartesian_pose_refiner), intent(inout) :: self
        real, intent(in) :: volume(:, :, :)
        logical, intent(in) :: apply_inverse_envelope
        type(image) :: padded_image
        real, allocatable :: prepared(:, :, :), inv1d(:)
        integer :: box, lims3(3, 2), wlims(2), lo, hi, i, j, k

        call self%kill
        box = size(volume, 1)
        if (box < 2 .or. mod(box, 2) /= 0 .or. size(volume, 2) /= box .or. size(volume, 3) /= box) &
            &error stop 'cartesian pose reference requires an even cubic volume'
        self%box = box
        self%padf = OSMPL_PAD_FAC
        self%boxpd = self%padf*box
        self%padsc = real(self%padf)**3
        self%kbwin = kbinterpol(KBWINSZ, KBALPHA)
        self%iwinsz = ceiling(self%kbwin%get_winsz() - 0.5)
        self%wdim = 2*self%iwinsz + 1
        self%lims2(1, :) = [-box/2, box/2]
        self%lims2(2, :) = [-box/2, box/2]
        allocate (prepared, source=volume)
        if (apply_inverse_envelope) then
            call kb_stencil_centered_crop_inv_envelope_1d(self%kbwin, self%boxpd, box, inv1d)
            do k = 1, box
                do j = 1, box
                    do i = 1, box
                        prepared(i, j, k) = prepared(i, j, k)*inv1d(i)*inv1d(j)*inv1d(k)
                    end do
                end do
            end do
            deallocate (inv1d)
        end if
        call padded_image%new([self%boxpd, self%boxpd, self%boxpd], 1.0)
        call padded_image%set_rmat(center_embed_real3d(prepared, self%boxpd), .false.)
        call padded_image%fft()
        lims3 = padded_image%loop_lims(3)
        wlims = lims3(2, :)
        lo = wlims(1) - self%iwinsz - 1
        hi = wlims(2) + self%iwinsz + 1
        allocate (self%wrap(lo:hi))
        do i = lo, hi
            self%wrap(i) = cyci_1d(wlims, i)
        end do
        self%cmat = padded_image%get_cmat()
        self%exists = .true.
        call padded_image%kill
        deallocate (prepared)
    end subroutine load_pose_reference

    !> Prepare one particle for repeated Cartesian pose evaluations.
    !! Store the whitened observation Y/sqrt(sigma2) and the shift-free transfer
    !! C/sqrt(sigma2) over the valid requested shells. The objective applies the
    !! candidate shift phase itself, so including a shift here would apply it twice.
    !! Invalid or unavailable shell variances produce data%valid=.false. without aborting.
    subroutine prepare_pose_particle(self, raw_observed, ctfparms, sigma2, requested_range, data)
        class(cartesian_pose_refiner), intent(in) :: self
        complex, intent(in) :: raw_observed(self%lims2(1, 1):self%lims2(1, 2), &
            &self%lims2(2, 1):self%lims2(2, 2))
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: requested_range(2)
        type(cartesian_pose_data), intent(out) :: data
        type(ctf) :: tfun
        type(ctfvars) :: ctfvals
        real :: cval, sigma
        integer :: h, k, shell, lower_shell, upper_shell, radius_squared
        logical :: use_ctf, phase_flip

        data%valid = .false.
        lower_shell = max(0, requested_range(1))
        upper_shell = min(requested_range(2), ubound(sigma2, 1), self%box/2)
        data%shell_range = [lower_shell, upper_shell]
        if (upper_shell < lower_shell) return
        do shell = lower_shell, upper_shell
            if (.not. ieee_is_finite(sigma2(shell)) .or. sigma2(shell) <= 0.0) return
        end do
        allocate (data%observed(self%lims2(1, 1):self%lims2(1, 2), &
            &self%lims2(2, 1):self%lims2(2, 2)), source=cmplx(0., 0.))
        allocate (data%transfer(self%lims2(1, 1):self%lims2(1, 2), &
            &self%lims2(2, 1):self%lims2(2, 2)), source=cmplx(0., 0.))
        allocate (data%sigma2(lower_shell:upper_shell), source=sigma2(lower_shell:upper_shell))
        use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
        phase_flip = ctfparms%ctfflag == CTFFLAG_FLIP
        if (use_ctf) then
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
            ctfvals = tfun%get_ctfvars(ctfparms%phshift)
        end if
        do k = self%lims2(2, 1), self%lims2(2, 2)
            do h = self%lims2(1, 1), self%lims2(1, 2)
                radius_squared = h*h + k*k
                if (radius_squared < lower_shell*lower_shell .or. &
                    &radius_squared > upper_shell*upper_shell) cycle
                shell = nint(sqrt(real(radius_squared)))
                sigma = sigma2(shell)
                cval = 1.0
                if (use_ctf) then
                    ! Use SIMPLE's CTF object without depending on process-global
                    ! Fourier maps initialized by the production matcher.
                    cval = pose_cont_ctf_value(tfun, h, k, self%box, ctfvals%phshift, phase_flip)
                end if
                data%transfer(h, k) = cval/sqrt(sigma)
                data%observed(h, k) = raw_observed(h, k)/sqrt(sigma)
            end do
        end do
        data%valid = .true.
    end subroutine prepare_pose_particle

    !> Evaluate the selected Cartesian objective and its five derivatives.
    !! Reuse the immutable reference grid and prepared particle so each LM trial changes
    !! only the three-component rotation increment and two native-pixel shifts.
    subroutine prepared_pose_objective_gradient(self, rotmat, shift, data, objective, gradient, &
        &objective_kind)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        real(dp), intent(out) :: objective, gradient(5)
        integer, optional, intent(in) :: objective_kind
        integer :: active_objective

        if (.not. data%valid) error stop 'prepared pose objective requires valid particle data'
        active_objective = POSE_CONT_OBJECTIVE_CART_EUCLID
        if (present(objective_kind)) active_objective = objective_kind
        call self%pose_objective_gradient(rotmat, shift, data%observed, objective, gradient, &
            &data%transfer, data%shell_range, active_objective)
    end subroutine prepared_pose_objective_gradient

    !> Recompute refine3D accounting at one valid terminal Cartesian pose.
    !! Undo particle whitening, evaluate the full-disk residual in native Fourier
    !! coordinates, and return per-shell residual, reference, and particle powers.
    !! This is intentionally separate from LM acceptance: both accepted and valid
    !! rejected terminal poses require accounting consistent with their final pose.
    subroutine prepared_pose_sigma_contribution(self, rotmat, shift, data, sigma_contrib, &
        &ref_pow, ptcl_pow, v)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        real, allocatable, intent(out) :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real, intent(out) :: v
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: model, raw_observed, residual
        real(dp), allocatable :: sigma_sum(:), ref_sum(:), ptcl_sum(:)
        real(dp) :: arg, root_sigma, vnum, vden
        real(sp) :: loc(3), switch_margin(3)
        integer, allocatable :: counts(:)
        integer :: h, k, shell, radius_squared, lower_shell, upper_shell

        if (.not. self%exists) error stop 'prepared sigma contribution requires a Fourier workspace'
        v = -1.
        if (.not. data%valid) return
        lower_shell = data%shell_range(1)
        upper_shell = data%shell_range(2)
        allocate (sigma_sum(lower_shell:upper_shell), source=0._dp)
        allocate (ref_sum(lower_shell:upper_shell), source=0._dp)
        allocate (ptcl_sum(lower_shell:upper_shell), source=0._dp)
        allocate (counts(lower_shell:upper_shell), source=0)
        vnum = 0._dp
        vden = 0._dp
        do k = self%lims2(2, 1), self%lims2(2, 2)
            do h = self%lims2(1, 1), self%lims2(1, 2)
                radius_squared = h*h + k*k
                if (radius_squared < lower_shell*lower_shell .or. &
                    &radius_squared > upper_shell*upper_shell) cycle
                shell = nint(sqrt(real(radius_squared)))
                loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                phase = cmplx(cos(arg), sin(arg), kind=sp)
                root_sigma = sqrt(real(data%sigma2(shell), dp))
                model = cmplx(phase, kind=dp)*cmplx(data%transfer(h, k), kind=dp)* &
                    &cmplx(value, kind=dp)*root_sigma
                raw_observed = cmplx(data%observed(h, k), kind=dp)*root_sigma
                residual = raw_observed - model
                sigma_sum(shell) = sigma_sum(shell) + real(conjg(residual)*residual, dp)
                ref_sum(shell) = ref_sum(shell) + real(conjg(model)*model, dp)
                ptcl_sum(shell) = ptcl_sum(shell) + real(conjg(raw_observed)*raw_observed, dp)
                counts(shell) = counts(shell) + 1
                vnum = vnum + real(conjg(residual)*residual, dp)/real(data%sigma2(shell), dp)
                vden = vden + real(conjg(raw_observed)*raw_observed, dp)/real(data%sigma2(shell), dp)
            end do
        end do
        if (any(counts == 0)) error stop 'prepared sigma contribution found an empty active shell'
        allocate (sigma_contrib(lower_shell:upper_shell))
        allocate (ref_pow(lower_shell:upper_shell))
        allocate (ptcl_pow(lower_shell:upper_shell))
        sigma_contrib = real(sigma_sum/(2._dp*real(counts, dp)), sp)
        ref_pow = real(ref_sum/real(counts, dp), sp)
        ptcl_pow = real(ptcl_sum/real(counts, dp), sp)
        if (vden > 0._dp) then
            v = real(vnum/vden, sp)
        else
            v = -1.
        end if
    end subroutine prepared_pose_sigma_contribution

    !> Refine one pose against an already prepared Cartesian particle.
    !! This is the production entry point: it couples the particle's immutable
    !! observation, transfer, and shell range to the bounded five-parameter LM.
    !! rotmat and shift change only through accepted LM transactions; result is
    !! always returned, while diagnostics is optional operational evidence.
    subroutine refine_prepared_pose_lm(self, rotmat, shift, data, config, result, diagnostics)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(inout) :: rotmat(3, 3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        type(pose_lm_config), intent(in) :: config
        type(pose_lm_result), intent(out) :: result
        type(pose_lm_diagnostics), optional, intent(out) :: diagnostics

        if (.not. data%valid) error stop 'prepared pose LM requires valid particle data'
        ! Forward the immutable prepared particle and grouped solver policy.
        call self%refine_pose_lm(rotmat, shift, data, config, result, diagnostics)
    end subroutine refine_prepared_pose_lm

    !> Apply a local tangent-space rotation increment on the right.
    !! Compute R_new = R exp([omega]x), matching the rotation derivative used by
    !! the pose Jacobian. omega is a three-component rotation vector in radians.
    pure function right_increment_rotation(rotmat, omega) result(updated_rotmat)
        real(dp), intent(in) :: rotmat(3, 3), omega(3)
        real(dp) :: updated_rotmat(3, 3), skew(3, 3), exp_skew(3, 3)
        real(dp) :: identity(3, 3), theta2, theta4, sinc_theta, cosc_theta

        identity = 0._dp
        identity(1, 1) = 1._dp
        identity(2, 2) = 1._dp
        identity(3, 3) = 1._dp
        ! [omega]x u = omega x u.
        skew = reshape([0._dp, omega(3), -omega(2), &
            &-omega(3), 0._dp, omega(1), omega(2), -omega(1), 0._dp], [3, 3])
        theta2 = dot_product(omega, omega)
        if (theta2 < 1.e-8_dp) then
            ! Taylor forms avoid cancellation as the rotation angle approaches zero.
            theta4 = theta2*theta2
            sinc_theta = 1._dp - theta2/6._dp + theta4/120._dp
            cosc_theta = 0.5_dp - theta2/24._dp + theta4/720._dp
        else
            sinc_theta = sin(sqrt(theta2))/sqrt(theta2)
            cosc_theta = (1._dp - cos(sqrt(theta2)))/theta2
        end if
        ! Rodrigues' formula evaluates the SO(3) exponential map.
        exp_skew = identity + sinc_theta*skew + cosc_theta*matmul(skew, skew)
        ! Right multiplication keeps omega in the current particle-pose frame.
        updated_rotmat = matmul(rotmat, exp_skew)
    end function right_increment_rotation

    subroutine kill_fourier_workspace(self)
        class(cartesian_pose_refiner), intent(inout) :: self
        if (allocated(self%wrap)) deallocate (self%wrap)
        if (allocated(self%cmat)) deallocate (self%cmat)
        self%box = 0
        self%boxpd = 0
        self%padf = 1
        self%iwinsz = 0
        self%wdim = 0
        self%lims2 = 0
        self%padsc = 1.0
        self%exists = .false.
    end subroutine kill_fourier_workspace

    !>  \brief  Samples the packed Fourier snapshot and its three fixed-cell
    !!          spatial derivatives at one oversampled-lattice coordinate.
    pure subroutine sample_fourier_with_grad(self, loc, value, dvalue_dloc, switch_margin)
        class(cartesian_pose_refiner), intent(in) :: self
        real(sp), intent(in) :: loc(3)
        complex, intent(out) :: value, dvalue_dloc(3)
        real(sp), intent(out) :: switch_margin(3)
        real(sp) :: w(self%wdim, self%wdim, self%wdim)
        real(sp) :: dw(self%wdim, self%wdim, self%wdim, 3)
        integer :: i0(3)
        if (.not. self%exists) error stop 'sample_with_grad called on an empty Fourier workspace'
        ! Build w and dw/dloc on the same fixed interpolation stencil.
        call self%kbwin%apod_mat_3d_fast_grad(loc, self%iwinsz, self%wdim, i0, switch_margin, w, dw)
        if (any(i0 < lbound(self%wrap, 1)) .or. &
            &any(i0 + self%wdim - 1 > ubound(self%wrap, 1))) then
            error stop 'sample_with_grad location lies outside the periodic wrap table'
        end if
        call gather_packed_window_grad(self%cmat, lbound(self%wrap, 1), self%wrap, &
            &i0, w, dw, value, dvalue_dloc)
        ! Apply native Fourier scaling to the value and all derivatives.
        value = self%padsc*value
        dvalue_dloc = self%padsc*dvalue_dloc
    end subroutine sample_fourier_with_grad

    !> Form the unweighted full-disk Cartesian prediction S(t) G(R)V.
    !! This is the reusable forward-model boundary; particle preparation applies
    !! CTF and shell whitening separately.
    subroutine predict_unweighted_pose(self, rotmat, shift, prediction)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        complex, intent(out) :: prediction(self%lims2(1, 1):self%lims2(1, 2),&
                                            &self%lims2(2, 1):self%lims2(2, 2))
        complex :: value, dvalue_dloc(3), phase
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg
        integer :: h, k
        if (.not. self%exists) error stop 'predict_unweighted called on an empty Fourier workspace'
        prediction = cmplx(0., 0.)
        do k = self%lims2(2, 1), self%lims2(2, 2)
            do h = self%lims2(1, 1), self%lims2(1, 2)
                if (h*h + k*k > (self%box/2)**2) cycle
                loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp*real(PI, dp)* &
                    &(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                phase = cmplx(cos(arg), sin(arg), kind=sp)
                prediction(h, k) = phase*value
            end do
        end do
    end subroutine predict_unweighted_pose

    !>  \brief  Fused Cartesian Euclidean/NCC shift objective, gradient and
    !!          two-by-two Gauss-Newton block.
    !!          Keep this dedicated path: it avoids derivative planes and a masked
    !!          five-by-five solve when only the two image shifts are active.
    subroutine shift_normal_terms(self, rotmat, shift, observed, objective, gradient, hessian, &
        &objective_kind, transfer, shell_range)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        complex, intent(in) :: observed(self%lims2(1, 1):self%lims2(1, 2),&
                                         &self%lims2(2, 1):self%lims2(2, 2))
        complex, optional, intent(in) :: transfer(self%lims2(1, 1):self%lims2(1, 2),&
                                                   &self%lims2(2, 1):self%lims2(2, 2))
        integer, optional, intent(in) :: shell_range(2)
        integer, intent(in) :: objective_kind
        real(dp), intent(out) :: objective, gradient(2), hessian(2, 2)
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: model, particle, residual, jacobian(2)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: arg, frequency(2), particle_power, prediction_power, cross_real
        real(dp) :: model_gradient(2)
        integer :: axis, h, jaxis, k, active_sqhp, active_sqlp
        if (.not. self%exists) error stop 'shift_normal_terms called on an empty Fourier workspace'
        active_sqhp = 0
        active_sqlp = (self%box/2)**2
        if (present(shell_range)) then
            if (shell_range(1) < 0 .or. shell_range(2) > self%box/2 .or. &
                &shell_range(2) < shell_range(1)) &
                &error stop 'shift_normal_terms shell range lies outside the native Fourier disk'
            active_sqhp = shell_range(1)*shell_range(1)
            active_sqlp = shell_range(2)*shell_range(2)
        end if
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        particle_power = 0._dp
        prediction_power = 0._dp
        cross_real = 0._dp
        model_gradient = 0._dp
        do k = self%lims2(2, 1), self%lims2(2, 2)
            do h = self%lims2(1, 1), self%lims2(1, 2)
                if (h*h + k*k < active_sqhp .or. h*h + k*k > active_sqlp) cycle
                loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                arg = 2._dp*real(PI, dp)* &
                    &(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                phase = cmplx(cos(arg), sin(arg), kind=sp)
                model = cmplx(phase*value, kind=dp)
                if (present(transfer)) model = model*cmplx(transfer(h, k), kind=dp)
                particle = cmplx(observed(h, k), kind=dp)
                residual = model - particle
                frequency = 2._dp*real(PI, dp)*real([h, k], dp)/real(self%box, dp)
                jacobian = cmplx(0._dp, frequency, kind=dp)*model
                select case (objective_kind)
                case (POSE_CONT_OBJECTIVE_CART_EUCLID)
                    objective = objective + 0.5_dp*real(conjg(residual)*residual, dp)
                    do axis = 1, 2
                        gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual, dp)
                        do jaxis = 1, 2
                            hessian(axis, jaxis) = hessian(axis, jaxis) + &
                                &real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                        end do
                    end do
                case (POSE_CONT_OBJECTIVE_CART_NCC)
                    particle_power = particle_power + real(conjg(particle)*particle, dp)
                    prediction_power = prediction_power + real(conjg(model)*model, dp)
                    cross_real = cross_real + real(conjg(particle)*model, dp)
                    do axis = 1, 2
                        gradient(axis) = gradient(axis) + real(conjg(particle)*jacobian(axis), dp)
                        model_gradient(axis) = model_gradient(axis) + &
                            &real(conjg(model)*jacobian(axis), dp)
                        do jaxis = 1, 2
                            hessian(axis, jaxis) = hessian(axis, jaxis) + &
                                &real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                        end do
                    end do
                case default
                    error stop 'shift_normal_terms received an unsupported objective'
                end select
            end do
        end do
        if (objective_kind == POSE_CONT_OBJECTIVE_CART_NCC) then
            call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, &
                &model_gradient, objective, gradient, hessian)
        end if
    end subroutine shift_normal_terms

    !>  \brief  Fused Cartesian Euclidean/NCC objective, five-vector gradient
    !!          and Gauss-Newton block for three rotations and two pixel shifts.
    subroutine pose_normal_terms(self, rotmat, shift, observed, objective, gradient, &
        &hessian, min_switch_margin, objective_kind, transfer, shell_range)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        complex, intent(in) :: observed(self%lims2(1, 1):self%lims2(1, 2),&
                                         &self%lims2(2, 1):self%lims2(2, 2))
        real(dp), intent(out) :: objective, gradient(5), hessian(5, 5), min_switch_margin
        complex, optional, intent(in) :: transfer(self%lims2(1, 1):self%lims2(1, 2),&
                                                   &self%lims2(2, 1):self%lims2(2, 2))
        integer, optional, intent(in) :: shell_range(2)
        integer, intent(in) :: objective_kind
        complex :: value, dvalue_dloc(3), phase
        complex(dp) :: weighted_phase, model, particle, residual, jacobian(5)
        real(sp) :: loc(3), switch_margin(3)
        real(dp) :: l_arg, dloc(3, 3), frequency(2)
        real(dp) :: particle_power, prediction_power, cross_real
        real(dp) :: model_gradient(5)
        integer :: active_sqhp, active_sqlp, axis, h, jaxis, k

        if (.not. self%exists) error stop 'pose_normal_terms called on an empty Fourier workspace'
        objective = 0._dp
        gradient = 0._dp
        hessian = 0._dp
        particle_power = 0._dp
        prediction_power = 0._dp
        cross_real = 0._dp
        model_gradient = 0._dp
        min_switch_margin = huge(0._dp)
        active_sqhp = 0
        active_sqlp = (self%box/2)**2
        if (present(shell_range)) then
            if (shell_range(1) < 0 .or. shell_range(2) > self%box/2 .or. &
                &shell_range(2) < shell_range(1)) &
                &error stop 'pose shell range lies outside the native Fourier disk'
            active_sqhp = shell_range(1)*shell_range(1)
            active_sqlp = shell_range(2)*shell_range(2)
        end if
        do k = self%lims2(2, 1), self%lims2(2, 2)
            do h = self%lims2(1, 1), self%lims2(1, 2)
                if (h*h + k*k < active_sqhp .or. h*h + k*k > active_sqlp) cycle
                loc = real(self%padf, sp)*real(matmul(real([h, k, 0], dp), rotmat), sp)
                call self%sample_with_grad(loc, value, dvalue_dloc, switch_margin)
                min_switch_margin = min(min_switch_margin, real(minval(switch_margin), dp))
                ! Columns are loc x e1, loc x e2 and loc x e3.
                dloc(:, 1) = [0._dp, real(loc(3), dp), -real(loc(2), dp)]
                dloc(:, 2) = [-real(loc(3), dp), 0._dp, real(loc(1), dp)]
                dloc(:, 3) = [real(loc(2), dp), -real(loc(1), dp), 0._dp]
                l_arg = 2._dp*real(PI, dp)*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(self%box, dp)
                phase = cmplx(cos(l_arg), sin(l_arg), kind=sp)
                weighted_phase = cmplx(phase, kind=dp)
                if (present(transfer)) weighted_phase = weighted_phase*cmplx(transfer(h, k), kind=dp)
                model = weighted_phase*cmplx(value, kind=dp)
                particle = cmplx(observed(h, k), kind=dp)
                residual = model - particle
                do axis = 1, 3
                    jacobian(axis) = weighted_phase*sum(cmplx(dvalue_dloc, kind=dp)*dloc(:, axis))
                end do
                frequency = 2._dp*real(PI, dp)*real([h, k], dp)/real(self%box, dp)
                jacobian(4:5) = cmplx(0._dp, frequency, kind=dp)*model
                select case (objective_kind)
                case (POSE_CONT_OBJECTIVE_CART_EUCLID)
                    objective = objective + 0.5_dp*real(conjg(residual)*residual, dp)
                    do axis = 1, 5
                        gradient(axis) = gradient(axis) + real(conjg(jacobian(axis))*residual, dp)
                        do jaxis = 1, 5
                            hessian(axis, jaxis) = hessian(axis, jaxis) + &
                                &real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                        end do
                    end do
                case (POSE_CONT_OBJECTIVE_CART_NCC)
                    particle_power = particle_power + real(conjg(particle)*particle, dp)
                    prediction_power = prediction_power + real(conjg(model)*model, dp)
                    cross_real = cross_real + real(conjg(particle)*model, dp)
                    do axis = 1, 5
                        gradient(axis) = gradient(axis) + real(conjg(particle)*jacobian(axis), dp)
                        model_gradient(axis) = model_gradient(axis) + &
                            &real(conjg(model)*jacobian(axis), dp)
                        do jaxis = 1, 5
                            hessian(axis, jaxis) = hessian(axis, jaxis) + &
                                &real(conjg(jacobian(axis))*jacobian(jaxis), dp)
                        end do
                    end do
                case default
                    error stop 'pose_normal_terms received an unsupported objective'
                end select
            end do
        end do
        if (objective_kind == POSE_CONT_OBJECTIVE_CART_NCC) then
            call finalize_cc_normal_terms(particle_power, prediction_power, cross_real, &
                &model_gradient, objective, gradient, hessian)
        end if
        if (min_switch_margin == huge(0._dp)) min_switch_margin = 0._dp
    end subroutine pose_normal_terms

    !> Convert correlation sufficient statistics into 1-correlation gradient
    !! and the Gauss-Newton block of the normalized model vector.
    subroutine finalize_cc_normal_terms(particle_power, prediction_power, cross_real, &
        &model_gradient, objective, gradient, hessian)
        real(dp), intent(in) :: particle_power, prediction_power, cross_real
        real(dp), intent(in) :: model_gradient(:)
        real(dp), intent(out) :: objective
        real(dp), intent(inout) :: gradient(:), hessian(:, :)
        real(dp) :: normalizer
        integer :: axis, jaxis

        if (particle_power <= POSE_NUMERIC_FLOOR .or. prediction_power <= POSE_NUMERIC_FLOOR) then
            objective = ieee_value(0._dp, ieee_quiet_nan)
            gradient = objective
            hessian = objective
            return
        end if
        normalizer = sqrt(particle_power*prediction_power)
        objective = 1._dp - cross_real/normalizer
        gradient = -gradient/normalizer + &
            &cross_real*model_gradient/(normalizer*prediction_power)
        do axis = 1, size(gradient)
            do jaxis = 1, size(gradient)
                hessian(axis, jaxis) = hessian(axis, jaxis)/prediction_power - &
                    &model_gradient(axis)*model_gradient(jaxis)/(prediction_power*prediction_power)
            end do
        end do
    end subroutine finalize_cc_normal_terms

    !> Return the joint pose objective and five-vector gradient without exposing
    !! the computed Gauss-Newton block or stencil margin.
    subroutine pose_objective_gradient(self, rotmat, shift, observed, objective, gradient, &
        &transfer, shell_range, objective_kind)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3), shift(2)
        complex, intent(in) :: observed(self%lims2(1, 1):self%lims2(1, 2),&
                                         &self%lims2(2, 1):self%lims2(2, 2))
        real(dp), intent(out) :: objective, gradient(5)
        complex, optional, intent(in) :: transfer(self%lims2(1, 1):self%lims2(1, 2),&
                                                   &self%lims2(2, 1):self%lims2(2, 2))
        integer, optional, intent(in) :: shell_range(2)
        integer, optional, intent(in) :: objective_kind
        real(dp) :: hessian(5, 5), min_switch_margin
        integer :: active_objective

        active_objective = POSE_CONT_OBJECTIVE_CART_EUCLID
        if (present(objective_kind)) active_objective = objective_kind
        call self%pose_normal_terms(rotmat, shift, observed, objective, gradient, hessian, &
            &min_switch_margin, active_objective, transfer, shell_range)
    end subroutine pose_objective_gradient

    !> Damped two-parameter Gauss-Newton refinement of one prepared particle.
    subroutine refine_shift_lm(self, rotmat, shift, data, config, result, diagnostics)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(in) :: rotmat(3, 3)
        real(dp), intent(inout) :: shift(2)
        type(cartesian_pose_data), intent(in) :: data
        type(shift_lm_config), intent(in) :: config
        type(pose_lm_result), intent(out) :: result
        type(pose_lm_diagnostics), optional, intent(out) :: diagnostics
        real(dp) :: gradient(2), hessian(2, 2), trial_gradient(2), trial_hessian(2, 2)
        real(dp) :: solve_matrix(2, 2), diagonal(2), direction(2), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, rejection_multiplier
        real(dp) :: det, predicted, actual, ratio, maxdiag
        real(dp) :: discriminant, lambda_max, lambda_min, step_norm, relative_reduction
        integer :: axis, iteration, naccepted
        logical :: bounded_trial
        result = pose_lm_result()
        if (present(diagnostics)) call diagnostics%reset
        if (.not. data%valid) then
            result%status = LM_INVALID_NUMERICS
            return
        end if
        if (config%max_iterations < 1) error stop 'refine_shift_lm requires at least one iteration'
        if (config%shift_step_bound <= 0._dp .or. .not. ieee_is_finite(config%shift_step_bound)) &
            &error stop 'refine_shift_lm shift step bound must be positive and finite'
        if (config%objective /= POSE_CONT_OBJECTIVE_CART_EUCLID .and. &
            &config%objective /= POSE_CONT_OBJECTIVE_CART_NCC) &
            &error stop 'refine_shift_lm received an unsupported objective'
        call self%shift_normal_terms(rotmat, shift, data%observed, objective, gradient, hessian, &
            &config%objective, data%transfer, data%shell_range)
        naccepted = 0
        bounded_trial = .false.
        if (.not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian))) then
            result%status = LM_INVALID_NUMERICS
            return
        end if
        mu = LM_INITIAL_DAMPING
        rejection_multiplier = LM_INITIAL_REJECTION_MULTIPLIER
        do iteration = 1, config%max_iterations
            result%niterations = iteration
            maxdiag = max(maxval([(hessian(axis, axis), axis=1, 2)]), 1._dp)
            ! Eigenvalues diagnose whether both shift directions are observable.
            discriminant = sqrt(max(0._dp, (hessian(1, 1) - hessian(2, 2))**2 + &
                &4._dp*hessian(1, 2)*hessian(2, 1)))
            lambda_max = 0.5_dp*(hessian(1, 1) + hessian(2, 2) + discriminant)
            lambda_min = 0.5_dp*(hessian(1, 1) + hessian(2, 2) - discriminant)
            if (lambda_max <= sqrt(epsilon(1._dp))*maxdiag .or. &
                &lambda_min <= sqrt(epsilon(1._dp))*max(lambda_max, 1._dp)) then
                result%status = LM_NO_RELIABLE_UPDATE
                exit
            end if
            if (sqrt(dot_product(gradient, gradient)) < 1.e-8_dp) then
                result%status = merge(LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, naccepted > 0)
                exit
            end if
            do axis = 1, 2
                diagonal(axis) = max(hessian(axis, axis), sqrt(epsilon(1._dp))*maxdiag, epsilon(1._dp))
            end do
            solve_matrix = hessian
            solve_matrix(1, 1) = solve_matrix(1, 1) + mu*diagonal(1)
            solve_matrix(2, 2) = solve_matrix(2, 2) + mu*diagonal(2)
            det = solve_matrix(1, 1)*solve_matrix(2, 2) - solve_matrix(1, 2)*solve_matrix(2, 1)
            if (abs(det) <= epsilon(1._dp)*maxdiag*maxdiag) then
                result%status = LM_NO_RELIABLE_UPDATE
                exit
            end if
            direction(1) = (-solve_matrix(2, 2)*gradient(1) + solve_matrix(1, 2)*gradient(2))/det
            direction(2) = (solve_matrix(2, 1)*gradient(1) - solve_matrix(1, 1)*gradient(2))/det
            if (any(.not. ieee_is_finite(direction))) then
                result%status = LM_INVALID_NUMERICS
                exit
            end if
            ! Shift coordinates are pixels; cap every trial displacement at the configured radius.
            step_norm = sqrt(dot_product(direction, direction))
            if (step_norm > config%shift_step_bound) then
                direction = direction*(config%shift_step_bound/step_norm)
                bounded_trial = .true.
                if (present(diagnostics)) diagnostics%nbound_hits = diagnostics%nbound_hits + 1
            end if
            step_norm = min(step_norm, config%shift_step_bound)
            if (present(diagnostics)) diagnostics%max_shift_step = max(diagnostics%max_shift_step, step_norm)
            predicted = -dot_product(gradient, direction) - 0.5_dp* &
                &dot_product(direction, matmul(hessian, direction))
            if (.not. ieee_is_finite(predicted)) then
                result%status = LM_INVALID_NUMERICS
                exit
            elseif (predicted <= 0._dp) then
                call increase_lm_damping(mu, rejection_multiplier)
                cycle
            end if
            trial_shift = shift + direction
            if (present(diagnostics)) diagnostics%nattempted = diagnostics%nattempted + 1
            call self%shift_normal_terms(rotmat, trial_shift, data%observed, trial_objective,&
                &trial_gradient, trial_hessian, config%objective, data%transfer, data%shell_range)
            if (.not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian))) then
                call increase_lm_damping(mu, rejection_multiplier)
                result%status = LM_INVALID_NUMERICS
                cycle
            end if
            actual = objective - trial_objective
            ratio = actual/predicted
            if (actual > 0._dp .and. ratio >= 0.25_dp) then
                relative_reduction = actual/max(abs(objective), POSE_NUMERIC_FLOOR)
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                naccepted = naccepted + 1
                if (present(diagnostics)) diagnostics%naccepted = naccepted
                if (ratio > 0.75_dp) mu = max(mu/2._dp, epsilon(1._dp))
                rejection_multiplier = LM_INITIAL_REJECTION_MULTIPLIER
                result%status = LM_ACCEPTED_IMPROVEMENT
                if (step_norm < 1.e-8_dp .or. relative_reduction < 1.e-10_dp) exit
            else
                call increase_lm_damping(mu, rejection_multiplier)
            end if
        end do
        if (result%status == LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial) &
            &result%status = LM_STEP_BOUND_REJECTED
    end subroutine refine_shift_lm

    !>  \brief  Scaled, bounded five-parameter LM refinement for a right
    !!          rotation increment and two image shifts.
    subroutine refine_pose_lm(self, rotmat, shift, data, config, result, diagnostics)
        class(cartesian_pose_refiner), intent(in) :: self
        real(dp), intent(inout) :: rotmat(3, 3), shift(2)
        type(cartesian_pose_data), intent(in) :: data
        type(pose_lm_config), intent(in) :: config
        type(pose_lm_result), intent(out) :: result
        type(pose_lm_diagnostics), optional, intent(out) :: diagnostics
        real(dp) :: gradient(5), hessian(5, 5), trial_gradient(5), trial_hessian(5, 5)
        real(dp) :: scaled_gradient(5), scaled_hessian(5, 5), solve_matrix(5, 5)
        real(dp) :: diagonal(5), scaled_direction(5), direction(5)
        real(dp) :: trial_rotmat(3, 3), trial_shift(2)
        real(dp) :: objective, trial_objective, mu, rejection_multiplier
        real(dp) :: predicted, actual, ratio, rotation_norm, shift_norm
        real(dp) :: relative_reduction, min_switch_margin, trial_switch_margin
        real(dp) :: cumulative_rotation, cumulative_shift, sine_half
        integer :: iteration, naccepted
        logical :: active(5), bounded_trial, bounded_step, cumulative_guard
        logical :: accept_trial, identifiable, reliable, stationary

        if (config%max_iterations < 1) error stop 'refine_pose_lm requires at least one LM iteration'
        if (config%rotation_scale <= 0._dp .or. .not. ieee_is_finite(config%rotation_scale)) &
            &error stop 'refine_pose_lm requires a positive finite rotation scale'
        result%status = LM_ITERATION_LIMIT
        result%niterations = 0
        naccepted = 0
        if (present(diagnostics)) call diagnostics%reset
        active = config%active_parameters
        if (.not. any(active)) error stop 'refine_pose_lm requires one active parameter'
        cumulative_guard = config%use_cumulative_guard
        if (cumulative_guard) then
            if (config%max_total_rotation <= 0._dp .or. config%max_total_shift <= 0._dp .or. &
                &.not. ieee_is_finite(config%max_total_rotation) .or. &
                &.not. ieee_is_finite(config%max_total_shift)) &
                &error stop 'refine_pose_lm cumulative bounds must be positive and finite'
            if (any(.not. ieee_is_finite(config%anchor_rotmat)) .or. &
                &any(.not. ieee_is_finite(config%anchor_shift))) &
                &error stop 'refine_pose_lm cumulative anchor must be finite'
        end if
        if (config%shift_step_bound <= 0._dp .or. .not. ieee_is_finite(config%shift_step_bound)) &
            &error stop 'refine_pose_lm shift step bound must be positive and finite'
        if (config%objective /= POSE_CONT_OBJECTIVE_CART_EUCLID .and. &
            &config%objective /= POSE_CONT_OBJECTIVE_CART_NCC) &
            &error stop 'refine_pose_lm received an unsupported objective'
        ! Evaluate the seed pose once. This full Fourier-disk traversal builds
        ! the objective, five-vector gradient, and 5-by-5 Gauss-Newton matrix.
        call self%pose_normal_terms(rotmat, shift, data%observed, objective, gradient, hessian, &
            &min_switch_margin, config%objective, data%transfer, data%shell_range)
        bounded_trial = .false.
        if (.not. ieee_is_finite(objective) .or. any(.not. ieee_is_finite(gradient)) .or. &
            &any(.not. ieee_is_finite(hessian))) then
            result%status = LM_INVALID_NUMERICS
            return
        end if

        mu = LM_INITIAL_DAMPING
        rejection_multiplier = LM_INITIAL_REJECTION_MULTIPLIER
        do iteration = 1, config%max_iterations
            result%niterations = iteration
            ! Form and solve the small damped LM system. This 5-by-5 algebra is
            ! inexpensive compared with evaluating the Fourier plane.
            call build_pose_lm_system(gradient, hessian, config%rotation_scale, mu, active, scaled_gradient, &
                &scaled_hessian, diagonal, solve_matrix, scaled_direction, direction, identifiable, &
                &stationary, reliable, bounded_step, config%shift_step_bound)
            ! Stop when the local system is stationary, unidentifiable, or
            ! numerically unreliable; retain any improvement already accepted.
            if (.not. identifiable) then
                result%status = merge(LM_ACCEPTED_IMPROVEMENT, LM_NO_RELIABLE_UPDATE, naccepted > 0)
                exit
            end if
            if (stationary) then
                result%status = merge(LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, naccepted > 0)
                exit
            end if
            if (.not. reliable) then
                result%status = LM_NO_RELIABLE_UPDATE
                exit
            end if
            if (any(.not. ieee_is_finite(direction))) then
                result%status = LM_INVALID_NUMERICS
                exit
            end if
            bounded_trial = bounded_trial .or. bounded_step
            if (bounded_step .and. present(diagnostics)) diagnostics%nbound_hits = diagnostics%nbound_hits + 1
            rotation_norm = sqrt(dot_product(direction(1:3), direction(1:3)))
            shift_norm = sqrt(dot_product(direction(4:5), direction(4:5)))
            if (present(diagnostics)) then
                diagnostics%max_rotation_step = max(diagnostics%max_rotation_step, rotation_norm)
                diagnostics%max_shift_step = max(diagnostics%max_shift_step, shift_norm)
            end if
            ! Predict the reduction from the local quadratic LM model:
            ! predicted = -g^T d - 1/2 d^T H d.
            predicted = -dot_product(gradient, direction) - 0.5_dp* &
                &dot_product(direction, matmul(hessian, direction))
            if (.not. ieee_is_finite(predicted)) then
                result%status = LM_INVALID_NUMERICS
                exit
            elseif (predicted <= 0._dp) then
                call increase_lm_damping(mu, rejection_multiplier)
                cycle
            end if
            ! Apply the proposed three-component SO(3) increment and the two
            ! image-shift increments without changing the accepted pose yet.
            trial_rotmat = right_increment_rotation(rotmat, direction(1:3))
            trial_shift = shift + direction(4:5)
            if (present(diagnostics)) diagnostics%nattempted = diagnostics%nattempted + 1
            if (cumulative_guard) then
                ! Reject proposals that leave the transaction's rotation or
                ! shift capture basin; increase damping and try a smaller step.
                sine_half = sqrt(sum((trial_rotmat - config%anchor_rotmat)**2))/(2._dp*sqrt(2._dp))
                cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
                cumulative_shift = sqrt(sum((trial_shift - config%anchor_shift)**2))
                if (cumulative_rotation > config%max_total_rotation + 10._dp*epsilon(1._dp) .or. &
                    &cumulative_shift > config%max_total_shift + 10._dp*epsilon(1._dp)) then
                    call increase_lm_damping(mu, rejection_multiplier)
                    if (.not. bounded_step .and. present(diagnostics)) &
                        &diagnostics%nbound_hits = diagnostics%nbound_hits + 1
                    bounded_trial = .true.
                    cycle
                end if
            end if
            ! Recompute the trial objective, gradient, and normal matrix from
            ! the full active Fourier disk. This repeated O(N^2) interpolation
            ! and accumulation is the dominant cost of each LM iteration.
            call self%pose_normal_terms(trial_rotmat, trial_shift, data%observed, trial_objective, &
                &trial_gradient, trial_hessian, trial_switch_margin, config%objective, &
                &data%transfer, data%shell_range)
            if (.not. ieee_is_finite(trial_objective) .or. any(.not. ieee_is_finite(trial_gradient)) .or. &
                &any(.not. ieee_is_finite(trial_hessian))) then
                call increase_lm_damping(mu, rejection_multiplier)
                result%status = LM_INVALID_NUMERICS
                cycle
            end if
            ! Compare the measured reduction with the quadratic prediction.
            ! Accept a trustworthy improvement; otherwise increase damping and
            ! retry from the last accepted pose.
            actual = objective - trial_objective
            ratio = actual/predicted
            accept_trial = actual > 0._dp .and. ratio >= 0.25_dp
            if (accept_trial) then
                relative_reduction = actual/max(abs(objective), POSE_NUMERIC_FLOOR)
                rotmat = trial_rotmat
                shift = trial_shift
                objective = trial_objective
                gradient = trial_gradient
                hessian = trial_hessian
                naccepted = naccepted + 1
                if (present(diagnostics)) diagnostics%naccepted = naccepted
                if (ratio > 0.75_dp) mu = max(mu/2._dp, epsilon(1._dp))
                rejection_multiplier = LM_INITIAL_REJECTION_MULTIPLIER
                result%status = LM_ACCEPTED_IMPROVEMENT
                ! Stop early after an accepted but negligible step or objective
                ! reduction; otherwise continue from this accepted endpoint.
                if (max(rotation_norm, shift_norm) < 1.e-8_dp .or. relative_reduction < 1.e-10_dp) exit
            else
                call increase_lm_damping(mu, rejection_multiplier)
            end if
        end do
        if (result%status == LM_ITERATION_LIMIT .and. naccepted == 0 .and. bounded_trial) &
            &result%status = LM_STEP_BOUND_REJECTED
    end subroutine refine_pose_lm

    !> Increase damping aggressively across consecutive rejected proposals.
    !! The owning solver resets the multiplier after an accepted proposal.
    pure subroutine increase_lm_damping(mu, rejection_multiplier)
        real(dp), intent(inout) :: mu, rejection_multiplier

        mu = min(mu*rejection_multiplier, LM_MAX_DAMPING)
        rejection_multiplier = min(2._dp*rejection_multiplier, LM_MAX_DAMPING)
    end subroutine increase_lm_damping

    !> Construct one scaled, damped, and independently bounded pose proposal.
    pure subroutine build_pose_lm_system(gradient, hessian, rotation_scale, mu, active, &
        &scaled_gradient, scaled_hessian, damping_diagonal, solve_matrix, scaled_step, &
        &physical_step, identifiable, stationary, reliable, bounded, shift_step_bound)
        real(dp), intent(in) :: gradient(5), hessian(5, 5), rotation_scale, mu
        logical, intent(in) :: active(5)
        real(dp), intent(out) :: scaled_gradient(5), scaled_hessian(5, 5)
        real(dp), intent(out) :: damping_diagonal(5), solve_matrix(5, 5)
        real(dp), intent(out) :: scaled_step(5), physical_step(5)
        logical, intent(out) :: identifiable, stationary, reliable, bounded
        real(dp), optional, intent(in) :: shift_step_bound
        real(dp) :: coordinate_scale(5), ignored_step(5), hessian_scale
        real(dp) :: rotation_norm, shift_norm, active_shift_bound
        integer :: axis, jaxis

        coordinate_scale = [rotation_scale, rotation_scale, rotation_scale, 1._dp, 1._dp]
        do axis = 1, 5
            scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
            do jaxis = 1, 5
                scaled_hessian(axis, jaxis) = &
                    &coordinate_scale(axis)*hessian(axis, jaxis)*coordinate_scale(jaxis)
            end do
        end do
        call apply_pose_parameter_mask(scaled_gradient, scaled_hessian, active)
        call solve_pose_cholesky(scaled_hessian, -scaled_gradient, ignored_step, identifiable)
        stationary = sqrt(dot_product(scaled_gradient, scaled_gradient)) < 1.e-8_dp
        damping_diagonal = 0._dp
        solve_matrix = scaled_hessian
        scaled_step = 0._dp
        physical_step = 0._dp
        reliable = .false.
        bounded = .false.
        if (.not. identifiable .or. stationary) return
        hessian_scale = max(maxval(abs(scaled_hessian)), POSE_NUMERIC_FLOOR)
        do axis = 1, 5
            damping_diagonal(axis) = max(scaled_hessian(axis, axis), &
                &sqrt(epsilon(1._dp))*hessian_scale, POSE_NUMERIC_FLOOR)
            solve_matrix(axis, axis) = solve_matrix(axis, axis) + mu*damping_diagonal(axis)
        end do
        call solve_pose_cholesky(solve_matrix, -scaled_gradient, scaled_step, reliable)
        if (.not. reliable) return
        physical_step = coordinate_scale*scaled_step
        active_shift_bound = 1._dp
        if (present(shift_step_bound)) active_shift_bound = shift_step_bound
        rotation_norm = sqrt(dot_product(physical_step(1:3), physical_step(1:3)))
        if (rotation_norm > rotation_scale) then
            physical_step(1:3) = physical_step(1:3)*(rotation_scale/rotation_norm)
            bounded = .true.
        end if
        shift_norm = sqrt(dot_product(physical_step(4:5), physical_step(4:5)))
        if (shift_norm > active_shift_bound) then
            physical_step(4:5) = physical_step(4:5)*(active_shift_bound/shift_norm)
            bounded = .true.
        end if
    end subroutine build_pose_lm_system

    !> Freeze inactive pose coordinates while retaining one five-vector LM path.
    pure subroutine apply_pose_parameter_mask(gradient, hessian, active)
        real(dp), intent(inout) :: gradient(5), hessian(5, 5)
        logical, intent(in) :: active(5)
        integer :: axis

        do axis = 1, 5
            if (active(axis)) cycle
            gradient(axis) = 0._dp
            hessian(axis, :) = 0._dp
            hessian(:, axis) = 0._dp
            hessian(axis, axis) = 1._dp
        end do
    end subroutine apply_pose_parameter_mask

    !>  \brief  Cholesky solve with a relative pivot test for a 5-by-5
    !!          symmetric positive-definite pose block.
    pure subroutine solve_pose_cholesky(matrix, rhs, solution, reliable)
        real(dp), intent(in) :: matrix(5, 5), rhs(5)
        real(dp), intent(out) :: solution(5)
        logical, intent(out) :: reliable
        real(dp) :: lower(5, 5), intermediate(5), pivot, pivot_floor, matrix_scale
        integer :: i, j

        solution = 0._dp
        intermediate = 0._dp
        lower = 0._dp
        reliable = .false.
        if (any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs))) return
        matrix_scale = maxval(abs(matrix))
        if (matrix_scale <= POSE_NUMERIC_FLOOR) return
        pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
        do i = 1, 5
            do j = 1, i - 1
                lower(i, j) = (matrix(i, j) - dot_product(lower(i, 1:j - 1), lower(j, 1:j - 1)))/lower(j, j)
            end do
            pivot = matrix(i, i) - dot_product(lower(i, 1:i - 1), lower(i, 1:i - 1))
            if (.not. ieee_is_finite(pivot) .or. pivot <= pivot_floor) return
            lower(i, i) = sqrt(pivot)
        end do
        do i = 1, 5
            intermediate(i) = (rhs(i) - dot_product(lower(i, 1:i - 1), intermediate(1:i - 1)))/lower(i, i)
        end do
        do i = 5, 1, -1
            solution(i) = (intermediate(i) - dot_product(lower(i + 1:5, i), solution(i + 1:5)))/lower(i, i)
        end do
        reliable = all(ieee_is_finite(solution))
    end subroutine solve_pose_cholesky

end module simple_cartesian_pose_refiner
