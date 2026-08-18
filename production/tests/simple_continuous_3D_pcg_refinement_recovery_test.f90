module continuous_3D_pcg_refinement_recovery_test
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
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
integer, parameter :: PARALLEL_EQUIVALENCE_PARTICLES = 8
real(dp), parameter :: ROTATION_SCALE = 0.10_dp
real(dp), parameter :: ROTATION_RECOVERY_TOL = 8.e-4_dp
real(dp), parameter :: SHIFT_RECOVERY_TOL = 2.e-3_dp

contains

!> Recover one injected joint rotation and shift from independently generated data.
subroutine run_pose_recovery()
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    type(pcg_pose_polish_summary) :: batch_summary, serial_summary, parallel_summary
    real, allocatable :: phantom(:,:,:)
    complex, allocatable :: observed(:,:), zero_plane(:,:), batch_observed(:,:,:), transfers(:,:,:)
    complex, allocatable :: equivalence_observed(:,:,:), equivalence_transfers(:,:,:)
    real(dp) :: truth_rotmat(3,3), estimate_rotmat(3,3), initial_rotmat(3,3)
    real(dp) :: truth_shift(2), estimate_shift(2), initial_shift(2)
    real(dp) :: accepted_objectives(0:MAX_RECOVERY_ITERATIONS)
    real(dp) :: initial_rotation_error, final_rotation_error, initial_shift_error, final_shift_error
    real(dp) :: max_rotation_step, max_shift_step, ignored_objective
    real(dp) :: recovery_max_rotation_step, recovery_max_shift_step
    real(dp) :: batch_rotmats(3,3,2), batch_shifts(2,2), weak_rotmat(3,3), weak_shift(2)
    real(dp) :: serial_rotmats(3,3,PARALLEL_EQUIVALENCE_PARTICLES)
    real(dp) :: parallel_rotmats(3,3,PARALLEL_EQUIVALENCE_PARTICLES)
    real(dp) :: serial_shifts(2,PARALLEL_EQUIVALENCE_PARTICLES)
    real(dp) :: parallel_shifts(2,PARALLEL_EQUIVALENCE_PARTICLES)
    integer :: lims2(2,2), naccepted, nattempted, nstencil_switches, status, statuses(2)
    integer :: recovery_naccepted, recovery_nattempted, recovery_nswitches
    integer :: serial_statuses(PARALLEL_EQUIVALENCE_PARTICLES)
    integer :: parallel_statuses(PARALLEL_EQUIVALENCE_PARTICLES)
    integer :: iparticle, saved_nthreads, parallel_nthreads

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

    ! The OpenMP path must preserve the serial pose transaction and summary.
    allocate(equivalence_observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2), &
        &PARALLEL_EQUIVALENCE_PARTICLES))
    allocate(equivalence_transfers(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2), &
        &PARALLEL_EQUIVALENCE_PARTICLES))
    do iparticle = 1, PARALLEL_EQUIVALENCE_PARTICLES
        if( mod(iparticle,2) == 1 )then
            serial_rotmats(:,:,iparticle) = initial_rotmat
            serial_shifts(:,iparticle) = initial_shift
            equivalence_observed(:,:,iparticle) = observed
            equivalence_transfers(:,:,iparticle) = cmplx(1.,0.)
        else
            serial_rotmats(:,:,iparticle) = weak_rotmat
            serial_shifts(:,iparticle) = weak_shift
            equivalence_observed(:,:,iparticle) = cmplx(0.,0.)
            equivalence_transfers(:,:,iparticle) = cmplx(0.,0.)
        endif
    enddo
    parallel_rotmats = serial_rotmats
    parallel_shifts = serial_shifts
    saved_nthreads = 1
    !$ saved_nthreads = omp_get_max_threads()
    parallel_nthreads = max(1,min(4,saved_nthreads))
    !$ call omp_set_num_threads(1)
    call polish_fixed_volume_poses(workspace,serial_rotmats,equivalence_observed,serial_shifts, &
        &ROTATION_SCALE,MAX_RECOVERY_ITERATIONS,serial_statuses,serial_summary,equivalence_transfers)
    !$ call omp_set_num_threads(parallel_nthreads)
    call polish_fixed_volume_poses(workspace,parallel_rotmats,equivalence_observed,parallel_shifts, &
        &ROTATION_SCALE,MAX_RECOVERY_ITERATIONS,parallel_statuses,parallel_summary,equivalence_transfers)
    !$ call omp_set_num_threads(saved_nthreads)
    call assert_true(all(parallel_statuses == serial_statuses), &
        &'parallel pose polish changed a terminal outcome')
    call assert_true(maxval(abs(parallel_rotmats-serial_rotmats)) == 0._dp .and. &
        &maxval(abs(parallel_shifts-serial_shifts)) == 0._dp, &
        &'parallel pose polish changed a retained pose')
    call assert_equivalent_summaries(serial_summary,parallel_summary)

    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_PCG_POSE initial/final rotation/shift errors: ', &
        &initial_rotation_error,final_rotation_error,initial_shift_error,final_shift_error
    write(*,'(a,2(es14.6,1x),3(i0,1x))') 'CONTINUOUS_3D_PCG_POSE max rotation/shift, accepted/tried/switches: ', &
        &recovery_max_rotation_step,recovery_max_shift_step,recovery_naccepted,recovery_nattempted,recovery_nswitches
    write(*,'(a,2(i0,1x))') 'CONTINUOUS_3D_PCG_POSE serial/parallel threads: ',1,parallel_nthreads
    write(*,'(a)') 'CONTINUOUS_3D_PCG_POSE_RECOVERY: PASS'
    deallocate(phantom,observed,zero_plane,batch_observed,transfers)
    deallocate(equivalence_observed,equivalence_transfers)
    call workspace%kill
    call pcgop%kill
end subroutine run_pose_recovery

!> Require identical deterministic reductions from serial and parallel batches.
subroutine assert_equivalent_summaries(serial_summary,parallel_summary)
    type(pcg_pose_polish_summary), intent(in) :: serial_summary, parallel_summary
    logical :: same_counts, same_reals

    same_counts = serial_summary%nparticles == parallel_summary%nparticles .and. &
        &serial_summary%nimproved == parallel_summary%nimproved .and. &
        &serial_summary%nunchanged == parallel_summary%nunchanged .and. &
        &serial_summary%nunreliable == parallel_summary%nunreliable .and. &
        &serial_summary%nstep_bound == parallel_summary%nstep_bound .and. &
        &serial_summary%ninvalid == parallel_summary%ninvalid .and. &
        &serial_summary%niteration_limit == parallel_summary%niteration_limit .and. &
        &serial_summary%naccepted_steps == parallel_summary%naccepted_steps .and. &
        &serial_summary%nattempted_steps == parallel_summary%nattempted_steps .and. &
        &serial_summary%nstencil_switches == parallel_summary%nstencil_switches
    same_reals = serial_summary%objective_before == parallel_summary%objective_before .and. &
        &serial_summary%objective_after == parallel_summary%objective_after .and. &
        &serial_summary%max_trial_step == parallel_summary%max_trial_step .and. &
        &serial_summary%max_rotation_step == parallel_summary%max_rotation_step .and. &
        &serial_summary%max_shift_step == parallel_summary%max_shift_step
    call assert_true(same_counts .and. same_reals, &
        &'parallel pose polish changed its deterministic summary')
end subroutine assert_equivalent_summaries

!> Return the geodesic angle between two proper rotation matrices.
pure function rotation_distance(left,right) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3)
    real(dp) :: distance, cosine
    cosine = 0.5_dp*(trace3(matmul(transpose(left),right))-1._dp)
    distance = acos(max(-1._dp,min(1._dp,cosine)))
end function rotation_distance

!> Return the trace used by the rotation-distance identity.
pure function trace3(matrix) result(trace)
    real(dp), intent(in) :: matrix(3,3)
    real(dp) :: trace
    trace = matrix(1,1)+matrix(2,2)+matrix(3,3)
end function trace3

end module continuous_3D_pcg_refinement_recovery_test
