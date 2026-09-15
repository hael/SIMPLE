!@descr: Reference, particle-data, and transaction adapters for refine3D pose_cont
module simple_pose_cont_refine3D_adapter
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api, only: CTFFLAG_FLIP, CTFFLAG_NO, CTFFLAG_YES, &
    &PI, TINY, ctfparams, del_file, dp, file_exists, find_img_smpd, find_ldim_nptcls, &
    &simple_exception, string
use simple_image, only: image
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &shift_lm_config, pose_lm_config, pose_lm_result, pose_lm_diagnostics, &
    &LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, LM_NO_RELIABLE_UPDATE, &
    &LM_STEP_BOUND_REJECTED, LM_INVALID_NUMERICS, LM_ITERATION_LIMIT
use simple_refine3D_fnames, only: refine3D_pose_cont_ref_fname
implicit none
private

#include "simple_local_flags.inc"

! Types
public :: cartesian_pose_data, pose_cont_reference_workspace
public :: pose_cont_pose, pose_cont_limits, pose_cont_config
public :: pose_cont_stage_result, pose_cont_transaction_result

! Refinement routes
public :: POSE_CONT_ROUTE_SHIFT_THEN_JOINT, POSE_CONT_ROUTE_JOINT

! Status codes
public :: POSE_CONT_NOT_ATTEMPTED, POSE_CONT_INVALID_PREPARATION
public :: LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, LM_NO_RELIABLE_UPDATE
public :: LM_STEP_BOUND_REJECTED, LM_INVALID_NUMERICS, LM_ITERATION_LIMIT

! Reference and particle adapters
public :: write_pose_cont_reference_artifact, remove_pose_cont_reference_artifacts
public :: prepare_pose_cont_observation
public :: shift_native_to_crop, shift_crop_to_native

integer, parameter  :: POSE_CONT_NOT_ATTEMPTED       = 0
integer, parameter  :: POSE_CONT_INVALID_PREPARATION = -1
integer, parameter  :: POSE_CONT_ROUTE_SHIFT_THEN_JOINT = 1
integer, parameter  :: POSE_CONT_ROUTE_JOINT            = 2

!> Rotation and shift coordinates for an input, staged, or returned pose.
type :: pose_cont_pose
    real(dp) :: rotmat(3, 3) = 0._dp
    real(dp) :: shift(2) = 0._dp
end type pose_cont_pose

!> Cropped-box pixel bounds applied by the shift-only and joint LM stages.
type :: pose_cont_limits
    real(dp) :: shift_step_bound = 0._dp
    real(dp) :: max_total_shift = 0._dp
end type pose_cont_limits

!> Algorithm choices for one pose_cont transaction.
!! Additional refine3D-facing pose_cont policy belongs here rather than in the
!! numerical limits or the Cartesian LM owner.
type :: pose_cont_config
    integer :: route = POSE_CONT_ROUTE_SHIFT_THEN_JOINT
    integer :: max_iterations = 40       !< maximum iterations in each enabled LM stage
    real(dp) :: rotation_scale = 0.1_dp  !< radians per joint-LM proposal
    real(dp) :: max_total_rotation = 15._dp*real(PI,dp)/180._dp !< radians from the seed
end type pose_cont_config

!> Result and accounting shared by the shift-only and joint LM stages.
type :: pose_cont_stage_result
    integer  :: status = POSE_CONT_NOT_ATTEMPTED
    integer  :: attempts = 0
    integer  :: accepts = 0
    integer  :: bound_hits = 0
    integer  :: stencil_switches = 0
    real(dp) :: objective_before = -1._dp
    real(dp) :: objective_after = -1._dp
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step = 0._dp
end type pose_cont_stage_result

!> Transactional result for either supported local LM route.
type :: pose_cont_transaction_result
    ! Terminal transaction decision and returned pose.
    integer :: status = POSE_CONT_NOT_ATTEMPTED
    type(pose_cont_pose) :: pose
    real(dp) :: objective_before = -1._dp
    real(dp) :: objective_after = -1._dp
    real(dp) :: cumulative_rotation = 0._dp
    real(dp) :: cumulative_shift = 0._dp

    ! Aggregate accounting across every LM stage executed by the transaction.
    integer  :: attempts = 0
    integer  :: accepts = 0
    integer  :: bound_hits = 0
    integer  :: stencil_switches = 0
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step = 0._dp

    ! Shift-only stage result and endpoint supplied to the joint stage.
    type(pose_cont_stage_result) :: shift_stage
    type(pose_cont_pose)         :: shift_endpoint

    ! Joint five-parameter stage result.
    type(pose_cont_stage_result) :: joint_stage
end type pose_cont_transaction_result

type :: pose_cont_reference_slot
    private
    type(cartesian_pose_refiner) :: refiner
    logical :: exists = .false.
end type pose_cont_reference_slot

!> Read-only Cartesian grids for all state/half references in one matcher pass.
type :: pose_cont_reference_workspace
    private
    type(pose_cont_reference_slot), allocatable :: even(:), odd(:)
    integer :: box = 0
    logical :: exists = .false.
contains
    ! Production reference lifecycle and particle operations.
    procedure :: new_from_artifacts => new_pose_cont_reference_workspace
    procedure :: kill => kill_pose_cont_reference_workspace
    procedure :: is_ready => pose_cont_reference_workspace_is_ready
    procedure :: prepare_particle => prepare_pose_cont_particle
    procedure :: refine_particle => refine_pose_cont_particle
end type pose_cont_reference_workspace

contains

    ! ========================================================================
    ! Production reference-artifact and workspace lifecycle
    ! ========================================================================

    !> Persist the exact processed real-space reference at the pre-PFTC boundary.
    subroutine write_pose_cont_reference_artifact(refvol, state, half)
        class(image), intent(inout) :: refvol
        integer, intent(in) :: state
        character(len=*), intent(in) :: half
        type(string) :: fname
        integer :: ldim(3)

        if (refvol%is_ft()) THROW_HARD('pose_cont reference artifact must be in real space')
        if (state < 1) THROW_HARD('pose_cont reference artifact requires a positive state')
        ldim = refvol%get_ldim()
        if (ldim(1) < 2 .or. mod(ldim(1), 2) /= 0 .or. any(ldim /= ldim(1))) &
            &THROW_HARD('pose_cont reference artifact must be an even cubic volume')
        call validate_half(half)
        fname = refine3D_pose_cont_ref_fname(state, half)
        call refvol%write(fname, del_if_exists=.true.)
        call fname%kill
    end subroutine write_pose_cont_reference_artifact

    subroutine remove_pose_cont_reference_artifacts(nstates)
        integer, intent(in) :: nstates
        type(string) :: fname
        integer :: state

        do state = 1, nstates
            fname = refine3D_pose_cont_ref_fname(state, 'even')
            if (file_exists(fname)) call del_file(fname)
            fname = refine3D_pose_cont_ref_fname(state, 'odd')
            if (file_exists(fname)) call del_file(fname)
        end do
        call fname%kill
    end subroutine remove_pose_cont_reference_artifacts

    !> Build one immutable Cartesian grid per state/half from validated artifacts.
    subroutine new_pose_cont_reference_workspace(self, nstates, box, smpd)
        class(pose_cont_reference_workspace), intent(inout) :: self
        integer, intent(in) :: nstates, box
        real, intent(in) :: smpd
        integer :: state

        call self%kill
        if (nstates < 1) THROW_HARD('pose_cont reference workspace requires at least one state')
        if (box < 2 .or. mod(box, 2) /= 0) THROW_HARD('pose_cont reference workspace requires an even box')
        if (smpd <= TINY .or. .not. ieee_is_finite(smpd)) &
            &THROW_HARD('pose_cont reference workspace requires positive finite sampling')
        allocate (self%even(nstates), self%odd(nstates))
        do state = 1, nstates
            call load_reference_slot(self%even(state), state, .true., box, smpd)
            call load_reference_slot(self%odd(state), state, .false., box, smpd)
        end do
        self%box = box
        self%exists = .true.
    end subroutine new_pose_cont_reference_workspace

    subroutine load_reference_slot(slot, state, even, box, smpd)
        type(pose_cont_reference_slot), intent(inout) :: slot
        integer, intent(in) :: state, box
        logical, intent(in) :: even
        real, intent(in) :: smpd
        type(image) :: refvol
        type(string) :: fname
        real, allocatable :: volume(:, :, :)
        real :: artifact_smpd
        integer :: ldim(3), nsections

        fname = refine3D_pose_cont_ref_fname(state, merge('even', 'odd ', even))
        if (.not. file_exists(fname)) THROW_HARD('missing pose_cont reference artifact: '//fname%to_char())
        call find_ldim_nptcls(fname, ldim, nsections)
        if (any(ldim /= [box, box, box]) .or. nsections /= box) &
            &THROW_HARD('incompatible pose_cont reference dimensions: '//fname%to_char())
        artifact_smpd = find_img_smpd(fname)
        if (.not. ieee_is_finite(artifact_smpd) .or. &
            &abs(artifact_smpd - smpd) > 10.*epsilon(smpd)*max(abs(smpd), 1.)) &
            &THROW_HARD('incompatible pose_cont reference sampling: '//fname%to_char())
        call refvol%new([box, box, box], smpd)
        call refvol%read(fname)
        volume = refvol%get_rmat()
        ! The artifact already represents refine3D's prepared physical map.
        call slot%refiner%new_physical_reference(volume)
        slot%exists = .true.
        call refvol%kill
        call fname%kill
        deallocate (volume)
    end subroutine load_reference_slot

    subroutine kill_pose_cont_reference_workspace(self)
        class(pose_cont_reference_workspace), intent(inout) :: self
        integer :: state

        if (allocated(self%even)) then
            do state = 1, size(self%even)
                call self%even(state)%refiner%kill
            end do
            deallocate (self%even)
        end if
        if (allocated(self%odd)) then
            do state = 1, size(self%odd)
                call self%odd(state)%refiner%kill
            end do
            deallocate (self%odd)
        end if
        self%box = 0
        self%exists = .false.
    end subroutine kill_pose_cont_reference_workspace

    pure logical function pose_cont_reference_workspace_is_ready(self, state, even) result(ready)
        class(pose_cont_reference_workspace), intent(in) :: self
        integer, intent(in) :: state
        logical, intent(in) :: even

        ready = self%exists
        if (.not. ready) return
        ready = state >= 1 .and. state <= size(self%even)
        if (.not. ready) return
        if (even) then
            ready = self%even(state)%exists
        else
            ready = self%odd(state)%exists
        end if
    end function pose_cont_reference_workspace_is_ready

    ! ========================================================================
    ! Production particle preparation and refinement adapters
    ! ========================================================================

    subroutine prepare_pose_cont_particle(self, state, even, observed, ctfparms, sigma2, &
        &requested_range, data)
        class(pose_cont_reference_workspace), intent(in) :: self
        integer, intent(in) :: state
        logical, intent(in) :: even
        complex, intent(in) :: observed(-self%box/2:self%box/2, &
            &-self%box/2:self%box/2)
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: requested_range(2)
        type(cartesian_pose_data), intent(out) :: data

        if (.not. self%is_ready(state, even)) THROW_HARD('pose_cont reference slot is not ready')
        if (even) then
            call self%even(state)%refiner%prepare_particle(observed, ctfparms, sigma2, &
                &requested_range, data)
        else
            call self%odd(state)%refiner%prepare_particle(observed, ctfparms, sigma2, &
                &requested_range, data)
        end if
    end subroutine prepare_pose_cont_particle

    subroutine refine_pose_cont_particle(self, state, even, seed, data, config, limits, result)
        class(pose_cont_reference_workspace), intent(in) :: self
        integer, intent(in) :: state
        logical, intent(in) :: even
        type(pose_cont_pose), intent(in) :: seed
        type(cartesian_pose_data), intent(in) :: data
        type(pose_cont_config), intent(in) :: config
        type(pose_cont_limits), intent(in) :: limits
        type(pose_cont_transaction_result), intent(out) :: result

        if (.not. self%is_ready(state, even)) THROW_HARD('pose_cont reference slot is not ready')
        if (even) then
            call run_pose_cont_transaction(self%even(state)%refiner, seed, data, config, limits, result)
        else
            call run_pose_cont_transaction(self%odd(state)%refiner, seed, data, config, limits, result)
        end if
    end subroutine refine_pose_cont_particle

    ! ========================================================================
    ! Production pose transaction
    ! ========================================================================

    !> Run the configured local LM route, then commit only when the complete
    !! transaction improves upon the original seed.
    subroutine run_pose_cont_transaction(refiner, seed, data, config, limits, result)
        class(cartesian_pose_refiner), intent(in) :: refiner
        type(pose_cont_pose), intent(in) :: seed
        type(cartesian_pose_data), intent(in) :: data
        type(pose_cont_config), intent(in) :: config
        type(pose_cont_limits), intent(in) :: limits
        type(pose_cont_transaction_result), intent(out) :: result
        real(dp) :: gradient(5)
        real(dp) :: shift_objective_after, joint_objective_after, sine_half
        type(pose_cont_pose) :: staged_pose
        type(shift_lm_config) :: shift_config
        type(pose_lm_config) :: joint_config
        type(pose_lm_result) :: lm_result
        type(pose_lm_diagnostics) :: diagnostics

        ! Initialize every returned pose to the seed so all early exits roll back.
        result = pose_cont_transaction_result()
        result%pose = seed
        result%shift_endpoint = seed

        ! Reject invalid prepared particle data before evaluating either solver.
        if (.not. data%is_valid()) then
            result%status = POSE_CONT_INVALID_PREPARATION
            return
        end if

        ! Bounds are caller policy and must be meaningful before solver setup.
        if (limits%shift_step_bound <= 0._dp .or. limits%max_total_shift <= 0._dp .or. &
            &.not. ieee_is_finite(limits%shift_step_bound) .or. &
            &.not. ieee_is_finite(limits%max_total_shift)) &
            &error stop 'pose_cont transaction requires positive finite shift bounds'

        ! The seed objective is the single acceptance baseline for the transaction.
        call refiner%prepared_objective_gradient(seed%rotmat, seed%shift, data, &
            &result%objective_before, gradient)
        if (.not. ieee_is_finite(result%objective_before)) then
            result%status = LM_INVALID_NUMERICS
            return
        end if

        staged_pose = seed
        shift_objective_after = result%objective_before

        select case (config%route)
        case (POSE_CONT_ROUTE_SHIFT_THEN_JOINT)
            ! Stage 1: refine translation only, holding the seed rotation fixed.
            shift_config = shift_lm_config(shift_step_bound=limits%shift_step_bound, &
                &max_iterations=config%max_iterations)
            call refiner%refine_shift_lm(staged_pose%rotmat, staged_pose%shift, data, &
                &shift_config, lm_result, diagnostics)
            call refiner%prepared_objective_gradient(staged_pose%rotmat, staged_pose%shift, data, &
                &shift_objective_after, gradient)
            call set_stage_result(result%shift_stage, lm_result, diagnostics, &
                &result%objective_before, shift_objective_after)
            result%shift_endpoint = staged_pose

            ! Enforce the caller's cumulative working-grid shift bound explicitly.
            if (sqrt(sum((staged_pose%shift - seed%shift)**2)) > &
                &limits%max_total_shift + 10._dp*epsilon(1._dp)) then
                result%shift_stage%status = LM_STEP_BOUND_REJECTED
                result%shift_stage%bound_hits = result%shift_stage%bound_hits + 1
            end if
            call add_stage_accounting(result, result%shift_stage)
            if (aborts_pose_cont_route(result%shift_stage%status)) then
                ! A failed shift stage cannot supply a trustworthy joint-stage seed.
                result%status = result%shift_stage%status
                result%objective_after = result%objective_before
                return
            end if
        case (POSE_CONT_ROUTE_JOINT)
            ! Direct joint LM starts from the original seed. The shift stage
            ! remains explicitly not attempted in the returned accounting.
        case default
            error stop 'pose_cont transaction received an invalid route'
        end select

        ! Joint stage: refine three rotations and two shifts from the route's
        ! current endpoint (the shift result or the original seed).
        ! Both cumulative guards remain anchored at the original transaction seed.
        joint_config = pose_lm_config(rotation_scale=config%rotation_scale, &
            &shift_step_bound=limits%shift_step_bound, max_iterations=config%max_iterations)
        joint_config%use_cumulative_guard = .true.
        joint_config%anchor_rotmat = seed%rotmat
        joint_config%anchor_shift = seed%shift
        joint_config%max_total_rotation = config%max_total_rotation
        joint_config%max_total_shift = limits%max_total_shift
        call refiner%refine_prepared_pose_lm(staged_pose%rotmat, staged_pose%shift, data, &
            &joint_config, lm_result, diagnostics)
        call refiner%prepared_objective_gradient(staged_pose%rotmat, staged_pose%shift, data, &
            &joint_objective_after, gradient)
        call set_stage_result(result%joint_stage, lm_result, diagnostics, &
            &shift_objective_after, joint_objective_after)
        call add_stage_accounting(result, result%joint_stage)
        if (aborts_pose_cont_route(result%joint_stage%status)) then
            ! Preserve the original seed whenever the joint stage is unreliable.
            result%status = result%joint_stage%status
            result%objective_after = result%objective_before
            return
        end if

        ! Commit the staged pose only if its final objective beats the seed.
        ! Otherwise select_pose_cont_terminal returns the seed unchanged.
        result%objective_after = result%joint_stage%objective_after
        call select_pose_cont_terminal(seed, staged_pose, result%objective_before, &
            &result%objective_after, result%status, result%pose)
        if (result%status == LM_INVALID_NUMERICS) &
            &result%objective_after = result%objective_before
        if (result%status /= LM_ACCEPTED_IMPROVEMENT) return

        ! Report cumulative motion only for a committed improving transaction.
        sine_half = sqrt(sum((result%pose%rotmat - seed%rotmat)**2))/(2._dp*sqrt(2._dp))
        result%cumulative_rotation = 2._dp*asin(max(0._dp, min(1._dp, sine_half)))
        result%cumulative_shift = sqrt(sum((result%pose%shift - seed%shift)**2))
    end subroutine run_pose_cont_transaction

    ! ========================================================================
    ! Internal production transaction helpers
    ! ========================================================================

    pure subroutine set_stage_result(stage, lm_result, diagnostics, objective_before, objective_after)
        type(pose_cont_stage_result), intent(out) :: stage
        type(pose_lm_result), intent(in) :: lm_result
        type(pose_lm_diagnostics), intent(in) :: diagnostics
        real(dp), intent(in) :: objective_before, objective_after

        stage%status = lm_result%status
        stage%attempts = diagnostics%nattempted
        stage%accepts = diagnostics%naccepted
        stage%bound_hits = diagnostics%nbound_hits
        stage%stencil_switches = diagnostics%nstencil_switches
        stage%objective_before = objective_before
        stage%objective_after = objective_after
        stage%max_rotation_step = diagnostics%max_rotation_step
        stage%max_shift_step = diagnostics%max_shift_step
    end subroutine set_stage_result

    pure subroutine add_stage_accounting(result, stage)
        type(pose_cont_transaction_result), intent(inout) :: result
        type(pose_cont_stage_result), intent(in) :: stage
        result%accepts = result%accepts + stage%accepts
        result%attempts = result%attempts + stage%attempts
        result%bound_hits = result%bound_hits + stage%bound_hits
        result%stencil_switches = result%stencil_switches + stage%stencil_switches
        result%max_rotation_step = max(result%max_rotation_step, stage%max_rotation_step)
        result%max_shift_step = max(result%max_shift_step, stage%max_shift_step)
    end subroutine add_stage_accounting

    pure logical function aborts_pose_cont_route(status) result(aborts)
        integer, intent(in) :: status
        select case (status)
        case (LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT)
            aborts = .false.
        case default
            aborts = .true.
        end select
    end function aborts_pose_cont_route

    pure subroutine select_pose_cont_terminal(input_pose, staged_pose, objective_before, &
        &objective_after, status, output_pose)
        type(pose_cont_pose), intent(in) :: input_pose, staged_pose
        real(dp), intent(in) :: objective_before, objective_after
        integer, intent(out) :: status
        type(pose_cont_pose), intent(out) :: output_pose
        logical :: accepted
        if (.not. ieee_is_finite(objective_before) .or. .not. ieee_is_finite(objective_after)) then
            status = LM_INVALID_NUMERICS
            output_pose = input_pose
            return
        end if
        accepted = objective_after < objective_before
        status = merge(LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, accepted)
        if (accepted) then
            output_pose = staged_pose
        else
            output_pose = input_pose
        end if
    end subroutine select_pose_cont_terminal

    ! ========================================================================
    ! Production observation preparation
    ! ========================================================================

    !> Prepare an unshifted, unflipped, masked full-disk Cartesian observation.
    subroutine prepare_pose_cont_observation(raw_img, noise_mask, work_img, mskrad, smpd_crop, &
        &ctfparms_in, observed, ctfparms_out)
        class(image), intent(inout) :: raw_img, work_img
        logical, intent(in) :: noise_mask(:, :, :)
        real, intent(in) :: mskrad, smpd_crop
        type(ctfparams), intent(in) :: ctfparms_in
        complex, allocatable, intent(out) :: observed(:, :)
        type(ctfparams), intent(out) :: ctfparms_out
        integer :: ldim(3)

        ldim = work_img%get_ldim()
        if (ldim(3) /= 1 .or. ldim(1) < 2 .or. mod(ldim(1), 2) /= 0 .or. ldim(2) /= ldim(1)) &
            &THROW_HARD('pose_cont observation workspace must be an even square image')
        if (any(shape(noise_mask) /= raw_img%get_ldim())) &
            &THROW_HARD('pose_cont observation noise mask has incompatible dimensions')
        if (smpd_crop <= TINY .or. .not. ieee_is_finite(smpd_crop)) &
            &THROW_HARD('pose_cont observation requires positive finite sampling')
        select case (ctfparms_in%ctfflag)
        case (CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP)
        case default
            THROW_HARD('unsupported CTF flag for pose_cont observation')
        end select
        call raw_img%norm_noise_fft_clip_shift(noise_mask, work_img, [0., 0.])
        call work_img%ifft_mask_fft(mskrad)
        allocate (observed(-ldim(1)/2:ldim(1)/2, -ldim(2)/2:ldim(2)/2))
        observed = work_img%expand_ft()
        ctfparms_out = ctfparms_in
        ctfparms_out%smpd = smpd_crop
    end subroutine prepare_pose_cont_observation

    ! ========================================================================
    ! Production coordinate and input-validation helpers
    ! ========================================================================

    pure function shift_native_to_crop(shift_native, box, box_crop) result(shift_crop)
        real, intent(in) :: shift_native(2)
        integer, intent(in) :: box, box_crop
        real :: shift_crop(2)
        if (box < 2 .or. box_crop < 2 .or. mod(box, 2) /= 0 .or. mod(box_crop, 2) /= 0) &
            &error stop 'pose_cont shift conversion requires positive even boxes'
        shift_crop = shift_native*real(box_crop)/real(box)
    end function shift_native_to_crop

    pure function shift_crop_to_native(shift_crop, box, box_crop) result(shift_native)
        real, intent(in) :: shift_crop(2)
        integer, intent(in) :: box, box_crop
        real :: shift_native(2)
        if (box < 2 .or. box_crop < 2 .or. mod(box, 2) /= 0 .or. mod(box_crop, 2) /= 0) &
            &error stop 'pose_cont shift conversion requires positive even boxes'
        shift_native = shift_crop*real(box)/real(box_crop)
    end function shift_crop_to_native

    !> Validate the half-set label before deriving an artifact filename.
    subroutine validate_half(half)
        character(len=*), intent(in) :: half
        if (trim(half) /= 'even' .and. trim(half) /= 'odd') &
            &THROW_HARD('pose_cont reference half must be even or odd')
    end subroutine validate_half

end module simple_pose_cont_refine3D_adapter
