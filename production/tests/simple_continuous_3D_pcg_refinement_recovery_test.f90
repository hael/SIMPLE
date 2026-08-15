module continuous_3D_pcg_refinement_recovery_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp
use simple_core_module_api, only: euler2m
use simple_pcg_pose_polisher, only: pcg_pose_polish_summary, polish_fixed_volume_poses
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, &
    &POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT, POSE_LM_NO_RELIABLE_UPDATE
implicit none
private
public :: run_pose_recovery

integer, parameter :: MAX_RECOVERY_ITERATIONS = 20
real(dp), parameter :: ROTATION_SCALE = 0.10_dp
real(dp), parameter :: ROTATION_RECOVERY_TOL = 8.e-4_dp
real(dp), parameter :: SHIFT_RECOVERY_TOL = 2.e-3_dp

contains

subroutine run_pose_recovery()
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    type(pcg_pose_polish_summary) :: batch_summary
    real, allocatable :: phantom(:,:,:)
    complex, allocatable :: observed(:,:), zero_plane(:,:), batch_observed(:,:,:), transfers(:,:,:)
    real(dp) :: truth_rotmat(3,3), estimate_rotmat(3,3), initial_rotmat(3,3)
    real(dp) :: truth_shift(2), estimate_shift(2), initial_shift(2)
    real(dp) :: accepted_objectives(0:MAX_RECOVERY_ITERATIONS)
    real(dp) :: initial_rotation_error, final_rotation_error, initial_shift_error, final_shift_error
    real(dp) :: max_rotation_step, max_shift_step, ignored_objective
    real(dp) :: recovery_max_rotation_step, recovery_max_shift_step
    real(dp) :: batch_rotmats(3,3,2), batch_shifts(2,2), weak_rotmat(3,3), weak_shift(2)
    integer :: lims2(2,2), naccepted, nattempted, nstencil_switches, status, statuses(2)
    integer :: recovery_naccepted, recovery_nattempted, recovery_nswitches

    call build_truth_volume(phantom)
    call pcgop%new(TRUTH_VOLUME_BOX,1._sp)
    call pcgop%set_volume(phantom)
    call pcgop%begin_fourier_workspace(workspace)
    call workspace%set_shell_range([2,TRUTH_VOLUME_BOX/2])
    lims2 = pcgop%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(zero_plane,mold=observed)
    zero_plane = cmplx(0.,0.)

    truth_rotmat = real(euler2m([19.,37.,28.]),dp)
    truth_shift = [0.31_dp,-0.24_dp]
    call workspace%shift_residual(truth_rotmat,truth_shift,zero_plane,observed,ignored_objective)
    ! The injected estimate is built independently from Euler angles.
    estimate_rotmat = real(euler2m([20.,36.2,28.7]),dp)
    estimate_shift = [-0.08_dp,0.06_dp]
    initial_rotation_error = rotation_distance(estimate_rotmat,truth_rotmat)
    initial_shift_error = sqrt(sum((estimate_shift-truth_shift)**2))
    call workspace%refine_pose_lm(estimate_rotmat,observed,estimate_shift,ROTATION_SCALE, &
        &MAX_RECOVERY_ITERATIONS,accepted_objectives,naccepted,status,nattempted, &
        &max_rotation_step,max_shift_step,nstencil_switches)
    final_rotation_error = rotation_distance(estimate_rotmat,truth_rotmat)
    final_shift_error = sqrt(sum((estimate_shift-truth_shift)**2))
    call assert_true(status == POSE_LM_ACCEPTED_IMPROVEMENT .and. naccepted > 0, &
        &'joint pose recovery accepted no improvement')
    call assert_true(nattempted >= naccepted, &
        &'joint pose recovery reported fewer trials than accepted steps')
    call assert_true(all(accepted_objectives(1:naccepted) < accepted_objectives(0:naccepted-1)), &
        &'an accepted joint pose step did not lower the fully recomputed objective')
    call assert_true(max_rotation_step <= ROTATION_SCALE+epsilon(1._dp), &
        &'joint pose recovery exceeded its rotation-step bound')
    call assert_true(max_shift_step <= 1._dp+epsilon(1._dp), &
        &'joint pose recovery exceeded its shift-step bound')
    call assert_true(final_rotation_error < initial_rotation_error .and. &
        &final_rotation_error < ROTATION_RECOVERY_TOL, &
        &'known injected rotation was not recovered')
    call assert_true(final_shift_error < initial_shift_error .and. final_shift_error < SHIFT_RECOVERY_TOL, &
        &'known injected shift was not recovered by the joint system')
    call assert_true(all(ieee_is_finite(accepted_objectives(0:naccepted))), &
        &'joint pose recovery produced a non-finite accepted objective')
    recovery_max_rotation_step = max_rotation_step
    recovery_max_shift_step = max_shift_step
    recovery_naccepted = naccepted
    recovery_nattempted = nattempted
    recovery_nswitches = nstencil_switches

    ! At the truth pose, a reliable particle must remain unchanged.
    estimate_rotmat = truth_rotmat
    estimate_shift = truth_shift
    call workspace%refine_pose_lm(estimate_rotmat,observed,estimate_shift,ROTATION_SCALE, &
        &MAX_RECOVERY_ITERATIONS,accepted_objectives,naccepted,status,nattempted, &
        &max_rotation_step,max_shift_step,nstencil_switches)
    call assert_true(status == POSE_LM_FINITE_NO_IMPROVEMENT .and. naccepted == 0, &
        &'an exact joint pose did not return finite_no_improvement')

    ! Batch rollback is atomic when a transfer removes every observable mode.
    allocate(batch_observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),2))
    allocate(transfers(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),2))
    initial_rotmat = real(euler2m([18.2,37.6,27.5]),dp)
    initial_shift = [0.02_dp,-0.03_dp]
    weak_rotmat = real(euler2m([20.1,37.4,27.8]),dp)
    weak_shift = [-0.2_dp,0.15_dp]
    batch_rotmats(:,:,1) = initial_rotmat
    batch_rotmats(:,:,2) = weak_rotmat
    batch_shifts(:,1) = initial_shift
    batch_shifts(:,2) = weak_shift
    batch_observed(:,:,1) = observed
    batch_observed(:,:,2) = cmplx(0.,0.)
    transfers(:,:,1) = cmplx(1.,0.)
    transfers(:,:,2) = cmplx(0.,0.)
    call polish_fixed_volume_poses(workspace,batch_rotmats,batch_observed,batch_shifts, &
        &ROTATION_SCALE,MAX_RECOVERY_ITERATIONS,statuses,batch_summary,transfers)
    call assert_true(statuses(1) == POSE_LM_ACCEPTED_IMPROVEMENT, &
        &'observable batch particle did not improve')
    call assert_true(statuses(2) == POSE_LM_NO_RELIABLE_UPDATE, &
        &'unobservable batch particle did not return no_reliable_update')
    call assert_true(maxval(abs(batch_rotmats(:,:,2)-weak_rotmat)) == 0._dp .and. &
        &maxval(abs(batch_shifts(:,2)-weak_shift)) == 0._dp, &
        &'unobservable batch particle changed despite atomic rollback')
    call assert_true(batch_summary%nimproved == 1 .and. batch_summary%nunreliable == 1, &
        &'joint pose batch summary misclassified its particles')

    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_PCG_POSE initial/final rotation/shift errors: ', &
        &initial_rotation_error,final_rotation_error,initial_shift_error,final_shift_error
    write(*,'(a,2(es14.6,1x),3(i0,1x))') 'CONTINUOUS_3D_PCG_POSE max rotation/shift, accepted/tried/switches: ', &
        &recovery_max_rotation_step,recovery_max_shift_step,recovery_naccepted,recovery_nattempted,recovery_nswitches
    write(*,'(a)') 'CONTINUOUS_3D_PCG_POSE_RECOVERY: PASS'
    deallocate(phantom,observed,zero_plane,batch_observed,transfers)
    call workspace%kill
    call pcgop%kill
end subroutine run_pose_recovery

pure function rotation_distance(left,right) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3)
    real(dp) :: distance, cosine
    cosine = 0.5_dp*(trace3(matmul(transpose(left),right))-1._dp)
    distance = acos(max(-1._dp,min(1._dp,cosine)))
end function rotation_distance

pure function trace3(matrix) result(trace)
    real(dp), intent(in) :: matrix(3,3)
    real(dp) :: trace
    trace = matrix(1,1)+matrix(2,2)+matrix(3,3)
end function trace3

end module continuous_3D_pcg_refinement_recovery_test
