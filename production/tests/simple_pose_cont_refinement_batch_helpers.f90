! Test-only batch accounting; this is not part of the production pose API.
module pose_cont_refinement_batch_helpers
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, &
    &SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT, &
    &SHIFT_LM_NO_RELIABLE_UPDATE, SHIFT_LM_STEP_BOUND_REJECTED, &
    &SHIFT_LM_INVALID_NUMERICS, SHIFT_LM_ITERATION_LIMIT, &
    &POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT, &
    &POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED, &
    &POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
implicit none
private

public :: pose_refinement_summary
public :: refine_fixed_volume_shifts, refine_fixed_volume_poses

!> Test-only aggregate of terminal outcomes, objectives and stencil telemetry.
type :: pose_refinement_summary
    integer :: nparticles = 0
    integer :: nimproved = 0
    integer :: nunchanged = 0
    integer :: nunreliable = 0
    integer :: nstep_bound = 0
    integer :: ninvalid = 0
    integer :: niteration_limit = 0
    integer :: naccepted_steps = 0
    integer :: nattempted_steps = 0
    real(dp) :: objective_before = 0._dp
    real(dp) :: objective_after = 0._dp
    real(dp) :: max_trial_step = 0._dp
    real(dp) :: max_rotation_step = 0._dp
    real(dp) :: max_shift_step = 0._dp
    integer :: nstencil_switches = 0
end type pose_refinement_summary

contains

!> Exercise shift refinement over a test batch and verify transactional rollback.
subroutine refine_fixed_volume_shifts(workspace, rotmats, observed, shifts, &
    &max_iterations, statuses, summary, transfers)
    type(cartesian_pose_refiner), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:)
    complex, intent(in) :: observed(:,:,:)
    real(dp), intent(inout) :: shifts(:,:)
    integer, intent(in) :: max_iterations
    integer, intent(out) :: statuses(:)
    type(pose_refinement_summary), intent(out) :: summary
    complex, optional, intent(in) :: transfers(:,:,:)
    real(dp), allocatable :: accepted_objectives(:)
    real(dp) :: input_shift(2), max_trial_step
    integer :: iparticle, naccepted, nattempted, status

    call validate_batch_shapes(workspace,rotmats,observed,shifts,statuses,max_iterations,transfers)
    summary = pose_refinement_summary()
    summary%nparticles = size(shifts,2)
    allocate(accepted_objectives(0:max_iterations))

    do iparticle = 1, summary%nparticles
        input_shift = shifts(:,iparticle)
        if( present(transfers) )then
            call workspace%refine_shift_lm(rotmats(:,:,iparticle),observed(:,:,iparticle), &
                &shifts(:,iparticle),max_iterations,accepted_objectives,naccepted,status, &
                &nattempted,max_trial_step,transfers(:,:,iparticle))
        else
            call workspace%refine_shift_lm(rotmats(:,:,iparticle),observed(:,:,iparticle), &
                &shifts(:,iparticle),max_iterations,accepted_objectives,naccepted,status, &
                &nattempted,max_trial_step)
        endif
        statuses(iparticle) = status
        summary%naccepted_steps = summary%naccepted_steps+naccepted
        summary%nattempted_steps = summary%nattempted_steps+nattempted
        summary%max_trial_step = max(summary%max_trial_step,max_trial_step)
        if( ieee_is_finite(accepted_objectives(0)) ) &
            &summary%objective_before = summary%objective_before+accepted_objectives(0)

        select case(status)
        case(SHIFT_LM_ACCEPTED_IMPROVEMENT)
            if( naccepted < 1 ) error stop 'accepted shift refinement has no accepted LM step'
            if( .not. ieee_is_finite(accepted_objectives(naccepted)) ) &
                &error stop 'accepted shift refinement has a non-finite objective'
            if( accepted_objectives(naccepted) >= accepted_objectives(0) ) &
                &error stop 'accepted shift refinement did not reduce the objective'
            summary%nimproved = summary%nimproved+1
            summary%objective_after = summary%objective_after+accepted_objectives(naccepted)
        case(SHIFT_LM_FINITE_NO_IMPROVEMENT)
            shifts(:,iparticle) = input_shift
            summary%nunchanged = summary%nunchanged+1
            summary%objective_after = summary%objective_after+accepted_objectives(0)
        case(SHIFT_LM_NO_RELIABLE_UPDATE)
            shifts(:,iparticle) = input_shift
            summary%nunreliable = summary%nunreliable+1
            summary%objective_after = summary%objective_after+accepted_objectives(0)
        case(SHIFT_LM_STEP_BOUND_REJECTED)
            shifts(:,iparticle) = input_shift
            summary%nstep_bound = summary%nstep_bound+1
            summary%objective_after = summary%objective_after+accepted_objectives(0)
        case(SHIFT_LM_INVALID_NUMERICS)
            shifts(:,iparticle) = input_shift
            summary%ninvalid = summary%ninvalid+1
            if( ieee_is_finite(accepted_objectives(0)) ) &
                &summary%objective_after = summary%objective_after+accepted_objectives(0)
        case(SHIFT_LM_ITERATION_LIMIT)
            shifts(:,iparticle) = input_shift
            summary%niteration_limit = summary%niteration_limit+1
            summary%objective_after = summary%objective_after+accepted_objectives(0)
        case default
            error stop 'shift refinement returned an unknown LM outcome'
        end select
    enddo
    deallocate(accepted_objectives)
end subroutine refine_fixed_volume_shifts

!> Exercise joint-pose refinement over independent particles in a fixed reference.
subroutine refine_fixed_volume_poses(workspace, rotmats, observed, shifts, rotation_scale, &
    &max_iterations, statuses, summary, transfers)
    type(cartesian_pose_refiner), intent(in) :: workspace
    real(dp), intent(inout) :: rotmats(:,:,:), shifts(:,:)
    complex, intent(in) :: observed(:,:,:)
    real(dp), intent(in) :: rotation_scale
    integer, intent(in) :: max_iterations
    integer, intent(out) :: statuses(:)
    type(pose_refinement_summary), intent(out) :: summary
    complex, optional, intent(in) :: transfers(:,:,:)
    real(dp), allocatable :: objectives_before(:), objectives_after(:)
    real(dp), allocatable :: max_rotation_steps(:), max_shift_steps(:)
    integer, allocatable :: accepted_steps(:), attempted_steps(:), stencil_switches(:)
    integer :: iparticle

    call validate_batch_shapes(workspace,rotmats,observed,shifts,statuses,max_iterations,transfers)
    if( rotation_scale <= 0._dp .or. .not. ieee_is_finite(rotation_scale) ) &
        &error stop 'pose refinement requires a positive finite rotation scale'
    summary = pose_refinement_summary()
    summary%nparticles = size(shifts,2)
    allocate(objectives_before(summary%nparticles),objectives_after(summary%nparticles),source=0._dp)
    allocate(max_rotation_steps(summary%nparticles),max_shift_steps(summary%nparticles),source=0._dp)
    allocate(accepted_steps(summary%nparticles),attempted_steps(summary%nparticles),source=0)
    allocate(stencil_switches(summary%nparticles),source=0)

    ! The reference is read-only; each particle owns its pose and result slots.
    !$omp parallel do default(shared) private(iparticle) schedule(static) proc_bind(close)
    do iparticle = 1, summary%nparticles
        call refine_one_pose(iparticle)
    enddo
    !$omp end parallel do

    ! Reduce in particle order so the test evidence remains deterministic.
    do iparticle = 1, summary%nparticles
        summary%naccepted_steps = summary%naccepted_steps+accepted_steps(iparticle)
        summary%nattempted_steps = summary%nattempted_steps+attempted_steps(iparticle)
        summary%max_rotation_step = max(summary%max_rotation_step,max_rotation_steps(iparticle))
        summary%max_shift_step = max(summary%max_shift_step,max_shift_steps(iparticle))
        summary%nstencil_switches = summary%nstencil_switches+stencil_switches(iparticle)
        summary%objective_before = summary%objective_before+objectives_before(iparticle)
        summary%objective_after = summary%objective_after+objectives_after(iparticle)
        select case(statuses(iparticle))
        case(POSE_LM_ACCEPTED_IMPROVEMENT)
            summary%nimproved = summary%nimproved+1
        case(POSE_LM_FINITE_NO_IMPROVEMENT)
            summary%nunchanged = summary%nunchanged+1
        case(POSE_LM_NO_RELIABLE_UPDATE)
            summary%nunreliable = summary%nunreliable+1
        case(POSE_LM_STEP_BOUND_REJECTED)
            summary%nstep_bound = summary%nstep_bound+1
        case(POSE_LM_INVALID_NUMERICS)
            summary%ninvalid = summary%ninvalid+1
        case(POSE_LM_ITERATION_LIMIT)
            summary%niteration_limit = summary%niteration_limit+1
        case default
            error stop 'pose refinement returned an unknown LM outcome'
        end select
    enddo
    deallocate(objectives_before,objectives_after,max_rotation_steps,max_shift_steps)
    deallocate(accepted_steps,attempted_steps,stencil_switches)

contains

    !> Solve one independent particle and retain or restore its complete pose.
    subroutine refine_one_pose(index)
        integer, intent(in) :: index
        real(dp) :: accepted_objectives(0:max_iterations)
        real(dp) :: input_rotmat(3,3), input_shift(2)
        integer :: naccepted, nattempted, status, nswitches

        input_rotmat = rotmats(:,:,index)
        input_shift = shifts(:,index)
        if( present(transfers) )then
            call workspace%refine_pose_lm(rotmats(:,:,index),observed(:,:,index), &
                &shifts(:,index),rotation_scale,max_iterations,accepted_objectives,naccepted, &
                &status,nattempted,max_rotation_steps(index),max_shift_steps(index),nswitches, &
                &transfers(:,:,index))
        else
            call workspace%refine_pose_lm(rotmats(:,:,index),observed(:,:,index), &
                &shifts(:,index),rotation_scale,max_iterations,accepted_objectives,naccepted, &
                &status,nattempted,max_rotation_steps(index),max_shift_steps(index),nswitches)
        endif
        statuses(index) = status
        accepted_steps(index) = naccepted
        attempted_steps(index) = nattempted
        stencil_switches(index) = nswitches
        if( ieee_is_finite(accepted_objectives(0)) )then
            objectives_before(index) = accepted_objectives(0)
            objectives_after(index) = accepted_objectives(0)
        endif

        if( status == POSE_LM_ACCEPTED_IMPROVEMENT )then
            if( naccepted < 1 ) error stop 'accepted pose refinement has no accepted LM step'
            if( .not. ieee_is_finite(accepted_objectives(naccepted)) ) &
                &error stop 'accepted pose refinement has a non-finite objective'
            if( accepted_objectives(naccepted) >= accepted_objectives(0) ) &
                &error stop 'accepted pose refinement did not reduce the objective'
            objectives_after(index) = accepted_objectives(naccepted)
        else
            rotmats(:,:,index) = input_rotmat
            shifts(:,index) = input_shift
        endif
    end subroutine refine_one_pose

end subroutine refine_fixed_volume_poses

!> Validate the shapes used only by the regression batch harness.
subroutine validate_batch_shapes(workspace, rotmats, observed, shifts, statuses, max_iterations, transfers)
    type(cartesian_pose_refiner), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:), shifts(:,:)
    complex, intent(in) :: observed(:,:,:)
    integer, intent(in) :: statuses(:), max_iterations
    complex, optional, intent(in) :: transfers(:,:,:)
    integer :: lims2(2,2), nparticles

    nparticles = size(shifts,2)
    lims2 = workspace%get_lims2()
    if( max_iterations < 1 ) error stop 'pose-refinement batch requires at least one LM iteration'
    if( size(shifts,1) /= 2 ) error stop 'pose-refinement batch requires two shift coordinates'
    if( size(rotmats,1) /= 3 .or. size(rotmats,2) /= 3 ) &
        &error stop 'pose-refinement batch requires 3-by-3 rotation matrices'
    if( size(rotmats,3) /= nparticles ) error stop 'pose-refinement batch rotation count mismatch'
    if( size(observed,1) /= lims2(1,2)-lims2(1,1)+1 .or. &
        &size(observed,2) /= lims2(2,2)-lims2(2,1)+1 ) &
        &error stop 'pose-refinement batch Fourier-plane shape mismatch'
    if( size(observed,3) /= nparticles ) error stop 'pose-refinement batch observation count mismatch'
    if( present(transfers) )then
        if( any(shape(transfers) /= shape(observed)) ) &
            &error stop 'pose-refinement batch transfer-plane shape mismatch'
    endif
    if( size(statuses) /= nparticles ) error stop 'pose-refinement batch status count mismatch'
end subroutine validate_batch_shapes

end module pose_cont_refinement_batch_helpers
