module pose_cont_refinement_solver_test
use pose_cont_refinement_test_helpers, only: assert_true, build_test_volume, &
    &prepare_unweighted_particle, rotation_distance, TEST_BOX
use simple_defs, only: dp
use simple_core_module_api, only: euler2m
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &shift_lm_config, pose_lm_config, pose_lm_result, pose_lm_diagnostics, &
    &POSE_CONT_OBJECTIVE_CART_NCC, &
    &LM_ACCEPTED_IMPROVEMENT, LM_FINITE_NO_IMPROVEMENT, LM_NO_RELIABLE_UPDATE, &
    &LM_STEP_BOUND_REJECTED
use simple_type_defs, only: ctfparams, CTFFLAG_NO
implicit none
private
public :: run_pose_cont_solver

real(dp), parameter :: ROTATION_TOL = 8.e-4_dp
real(dp), parameter :: SHIFT_TOL = 2.e-3_dp

contains

    subroutine run_pose_cont_solver()
        call test_shift_solver()
        call test_joint_solver()
        call test_ncc_solver()
        call test_invalid_and_unobservable_inputs()
        write (*, '(a)') 'POSE_CONT_REFINEMENT_SOLVER: PASS'
    end subroutine run_pose_cont_solver

    subroutine test_shift_solver()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        type(shift_lm_config) :: config
        type(pose_lm_result) :: result
        type(pose_lm_diagnostics) :: diagnostics
        real, allocatable :: volume(:, :, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3, 3), truth_shift(2), shift(2)

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift = [0.31_dp, -0.24_dp]
        call workspace%predict_unweighted(rotation, truth_shift, observed)
        call prepare_unweighted_particle(workspace, observed, data)
        config = shift_lm_config(shift_step_bound=1._dp, max_iterations=20)
        shift = [-0.15_dp, 0.12_dp]
        call workspace%refine_shift_lm(rotation, shift, data, config, result, diagnostics)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &sqrt(sum((shift - truth_shift)**2)) < SHIFT_TOL, &
            &'shift-only LM did not recover a known shift')
        call assert_true(diagnostics%naccepted > 0 .and. &
            &diagnostics%max_shift_step <= config%shift_step_bound + epsilon(1._dp), &
            &'shift-only LM violated its accepted-step contract')
        shift = truth_shift
        call workspace%refine_shift_lm(rotation, shift, data, config, result, diagnostics)
        call assert_true(result%status == LM_FINITE_NO_IMPROVEMENT .and. &
            &all(shift == truth_shift), 'exact shift was not retained')
        call workspace%kill
    end subroutine test_shift_solver

    subroutine test_joint_solver()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        type(pose_lm_config) :: config
        type(pose_lm_result) :: result
        type(pose_lm_diagnostics) :: diagnostics
        real, allocatable :: volume(:, :, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: truth_rotation(3, 3), rotation(3, 3), frozen_rotation(3, 3)
        real(dp) :: truth_shift(2), shift(2), frozen_shift(2), objective_before, objective_after
        real(dp) :: gradient(5)

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        truth_rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift = [0.31_dp, -0.24_dp]
        call workspace%predict_unweighted(truth_rotation, truth_shift, observed)
        call prepare_unweighted_particle(workspace, observed, data)
        config = pose_lm_config(rotation_scale=0.10_dp, shift_step_bound=1._dp, max_iterations=20)
        rotation = real(euler2m([20., 36.2, 28.7]), dp)
        shift = [-0.08_dp, 0.06_dp]
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_before, gradient)
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_after, gradient)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. objective_after < objective_before, &
            &'joint LM did not accept an objective-reducing pose')
        call assert_true(rotation_distance(rotation, truth_rotation) < ROTATION_TOL .and. &
            &sqrt(sum((shift - truth_shift)**2)) < SHIFT_TOL, &
            &'joint LM did not recover the known five-parameter pose')
        call assert_true(diagnostics%max_rotation_step <= config%rotation_scale + epsilon(1._dp) .and. &
            &diagnostics%max_shift_step <= config%shift_step_bound + epsilon(1._dp), &
            &'joint LM exceeded a configured proposal bound')

        rotation = truth_rotation
        shift = truth_shift
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call assert_true(result%status == LM_FINITE_NO_IMPROVEMENT .and. &
            &all(rotation == truth_rotation) .and. all(shift == truth_shift), &
            &'exact joint pose was not retained')

        rotation = real(euler2m([20., 36.2, 28.7]), dp)
        shift = [-0.08_dp, 0.06_dp]
        frozen_rotation = rotation
        frozen_shift = shift
        config%active_parameters = [.false., .false., .false., .true., .true.]
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_before, gradient)
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_after, gradient)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &objective_after < objective_before .and. sqrt(sum((shift - frozen_shift)**2)) > 1.e-10_dp .and. &
            &all(rotation == frozen_rotation), 'shift-only joint solve did not improve only active shifts')

        rotation = real(euler2m([20., 36.2, 28.7]), dp)
        shift = [-0.08_dp, 0.06_dp]
        frozen_rotation = rotation
        frozen_shift = shift
        config%active_parameters = [.true., .true., .true., .false., .false.]
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_before, gradient)
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_after, gradient)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &objective_after < objective_before .and. rotation_distance(rotation, frozen_rotation) > 1.e-10_dp .and. &
            &all(shift == frozen_shift), 'rotation-only joint solve did not improve only active rotations')

        rotation = real(euler2m([20., 36.2, 28.7]), dp)
        shift = [-0.08_dp, 0.06_dp]
        frozen_rotation = rotation
        frozen_shift = shift
        config%active_parameters = .true.
        config%use_cumulative_guard = .true.
        config%anchor_rotmat = rotation
        config%anchor_shift = shift
        config%max_total_rotation = 1.e-12_dp
        config%max_total_shift = 1.e-12_dp
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call assert_true(result%status == LM_STEP_BOUND_REJECTED .and. &
            &all(rotation == frozen_rotation) .and. all(shift == frozen_shift), &
            &'cumulative-bound rejection changed the complete input pose')
        call workspace%kill
    end subroutine test_joint_solver

    subroutine test_ncc_solver()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        type(pose_lm_config) :: config
        type(pose_lm_result) :: result
        type(pose_lm_diagnostics) :: diagnostics
        real, allocatable :: volume(:, :, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: truth_rotation(3, 3), rotation(3, 3), truth_shift(2), shift(2)
        real(dp) :: objective_before, objective_after, gradient(5)

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        truth_rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift = [0.31_dp, -0.24_dp]
        call workspace%predict_unweighted(truth_rotation, truth_shift, observed)
        call prepare_unweighted_particle(workspace, 2.*observed, data)
        config = pose_lm_config(rotation_scale=0.10_dp, shift_step_bound=1._dp, &
            &max_iterations=20, objective=POSE_CONT_OBJECTIVE_CART_NCC)
        rotation = real(euler2m([20., 36.2, 28.7]), dp)
        shift = [-0.08_dp, 0.06_dp]
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_before, gradient, &
            &POSE_CONT_OBJECTIVE_CART_NCC)
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective_after, gradient, &
            &POSE_CONT_OBJECTIVE_CART_NCC)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &objective_after < objective_before, &
            &'Cartesian NCC LM did not accept a correlation-improving pose')
        call assert_true(rotation_distance(rotation, truth_rotation) < ROTATION_TOL .and. &
            &sqrt(sum((shift - truth_shift)**2)) < SHIFT_TOL, &
            &'Cartesian NCC LM did not recover the known gain-scaled pose')
        call workspace%kill()
    end subroutine test_ncc_solver

    subroutine test_invalid_and_unobservable_inputs()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        type(pose_lm_config) :: config
        type(pose_lm_result) :: result
        type(pose_lm_diagnostics) :: diagnostics
        type(ctfparams) :: no_ctf
        real :: sigma2(0:TEST_BOX/2), zero_volume(TEST_BOX, TEST_BOX, TEST_BOX)
        complex :: zero_plane(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3, 3), original_rotation(3, 3), shift(2), original_shift(2)

        zero_volume = 0.
        zero_plane = cmplx(0., 0.)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = 1.
        sigma2(2) = -1.
        call workspace%new_physical_reference(zero_volume)
        call workspace%prepare_particle(zero_plane, no_ctf, sigma2, [2, TEST_BOX/2], data)
        call assert_true(.not. data%is_valid(), 'invalid sigma data was accepted')
        sigma2 = 1.
        call workspace%prepare_particle(zero_plane, no_ctf, sigma2, [2, TEST_BOX/2], data)
        rotation = real(euler2m([19., 37., 28.]), dp)
        shift = [0.2_dp, -0.1_dp]
        original_rotation = rotation
        original_shift = shift
        config = pose_lm_config(rotation_scale=0.10_dp, max_iterations=10)
        call workspace%refine_prepared_pose_lm(rotation, shift, data, config, result, diagnostics)
        call assert_true(result%status == LM_NO_RELIABLE_UPDATE .and. &
            &all(rotation == original_rotation) .and. all(shift == original_shift), &
            &'unobservable joint solve changed the input pose')
        call workspace%kill
    end subroutine test_invalid_and_unobservable_inputs

end module pose_cont_refinement_solver_test
