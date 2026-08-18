module continuous_3D_pcg_refinement_reference_bias_test
use continuous_3D_pcg_refinement_forward_path_test, only: build_forward_path_observations
use continuous_3D_pcg_refinement_halfset_support, only: build_disjoint_half_orientations, &
    &reconstruct_half_fixed, reconstruct_half_trajectory, HALFSET_LAMBDA, HALFSET_MASK_RADIUS, &
    &HALFSET_SMPD, centered_array_correlation
use continuous_3D_pcg_refinement_matched_window_test, only: matched_reference_metrics, &
    &evaluate_matched_reference_batch, gauge_corrected_errors
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp
use simple_ori, only: ori
use simple_oris, only: oris
use simple_reconstructor_pcg, only: reconstructor_pcg, right_increment_rotation
implicit none
private
public :: run_reference_bias_diagnostic

integer, parameter :: NPARTICLES = 48
integer, parameter :: NPROBES = 8
integer, parameter :: NITERATIONS = 3
integer, parameter :: NLAMBDAS = 3
integer, parameter :: NARMS = 22
integer, parameter :: ITERATION_COUNTS(NITERATIONS) = [2,4,8]
real, parameter :: LAMBDAS(NLAMBDAS) = [1.,100.,2000.]
integer, parameter :: SHELL_LOW(2) = [2,4]
integer, parameter :: SHELL_MID(2) = [2,8]
integer, parameter :: SHELL_FULL(2) = [2,BOX/2]
real(dp), parameter :: ROTATION_TOLERANCE = 0.001745_dp
real(dp), parameter :: SHIFT_TOLERANCE = 0.05_dp
real(dp), parameter :: RECOVERY_ROTATION_INITIAL = sqrt(0.010_dp**2+0.008_dp**2+0.006_dp**2)
real(dp), parameter :: RECOVERY_SHIFT_INITIAL = 0.25_dp
character(len=32), parameter :: ARM_NAMES(NARMS) = [character(len=32) :: &
    &'truth_full', 'own_i2_full', 'own_i4_full', 'own_i8_full', &
    &'cross_i2_full', 'cross_i4_full', 'cross_i8_full', 'holdout_i2_full', &
    &'own_i2_low', 'own_i2_mid', 'cross_i2_low', 'cross_i2_mid', &
    &'holdout_i2_low', 'holdout_i2_mid', 'own_i2_unmasked', 'cross_i2_unmasked', &
    &'own_lambda1_full', 'cross_lambda1_full', 'own_lambda100_full', &
    &'cross_lambda100_full', 'own_lambda2000_full', 'cross_lambda2000_full']

contains

!> Separate reconstructed-reference, ownership, iteration, frequency, support, and lambda bias.
subroutine run_reference_bias_diagnostic()
    type(oris) :: even_oris, odd_oris, holdout_oris
    type(ori) :: orientation
    type(reconstructor_pcg) :: sampler
    type(matched_reference_metrics) :: metrics(NARMS)
    complex, allocatable :: even_stages(:,:,:,:), odd_stages(:,:,:,:)
    complex, allocatable :: probe_observed(:,:,:), holdout_planes(:,:,:)
    real, allocatable :: truth_volume(:,:,:)
    real, allocatable :: own_trajectory(:,:,:,:), cross_trajectory(:,:,:,:)
    real, allocatable :: own_unmasked(:,:,:,:), cross_unmasked(:,:,:,:)
    real, allocatable :: holdout_volume(:,:,:), lambda_own(:,:,:,:), lambda_cross(:,:,:,:)
    real, allocatable :: temporary_volume(:,:,:)
    real(dp) :: exact_rotmats(3,3,NPROBES), data_residual
    integer, allocatable :: even_ids(:), odd_ids(:), niters(:)
    integer :: arm, i, ilambda, ignored_niters

    call build_disjoint_half_orientations(NPARTICLES,even_oris,odd_oris,even_ids,odd_ids)
    call build_truth_volume(truth_volume)
    call sampler%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call build_forward_path_observations(sampler,even_oris,truth_volume,even_stages)
    call build_forward_path_observations(sampler,odd_oris,truth_volume,odd_stages)
    allocate(probe_observed,mold=even_stages(:,:,1:NPROBES,3))
    probe_observed = even_stages(:,:,1:NPROBES,3)
    call orientation%new(.false.)
    do i = 1, NPROBES
        call even_oris%get_ori(i,orientation)
        exact_rotmats(:,:,i) = real(orientation%get_mat(),dp)
    enddo
    call test_gauge_removal(exact_rotmats)

    call reconstruct_half_trajectory(even_oris,even_stages(:,:,:,3),ITERATION_COUNTS, &
        &HALFSET_MASK_RADIUS,own_trajectory,niters)
    call reconstruct_half_trajectory(odd_oris,odd_stages(:,:,:,3),ITERATION_COUNTS, &
        &HALFSET_MASK_RADIUS,cross_trajectory,niters)
    call reconstruct_half_trajectory(even_oris,even_stages(:,:,:,3),[2],0.,own_unmasked,niters)
    call reconstruct_half_trajectory(odd_oris,odd_stages(:,:,:,3),[2],0.,cross_unmasked,niters)

    ! Hold out every probe observation while retaining the remaining even-half angular coverage.
    call holdout_oris%new(NPARTICLES-NPROBES,.false.)
    allocate(holdout_planes(-BOX/2:BOX/2,-BOX/2:BOX/2,NPARTICLES-NPROBES))
    do i = NPROBES+1, NPARTICLES
        call even_oris%get_ori(i,orientation)
        call holdout_oris%set_ori(i-NPROBES,orientation)
        holdout_planes(:,:,i-NPROBES) = even_stages(:,:,i,3)
    enddo
    call reconstruct_half_fixed(holdout_oris,holdout_planes,HALFSET_LAMBDA,2, &
        &HALFSET_MASK_RADIUS,holdout_volume,ignored_niters,data_residual)

    allocate(lambda_own(BOX,BOX,BOX,NLAMBDAS),lambda_cross(BOX,BOX,BOX,NLAMBDAS))
    do ilambda = 1, NLAMBDAS
        call reconstruct_half_fixed(even_oris,even_stages(:,:,:,3),LAMBDAS(ilambda),2, &
            &HALFSET_MASK_RADIUS,temporary_volume,ignored_niters,data_residual)
        lambda_own(:,:,:,ilambda) = temporary_volume
        deallocate(temporary_volume)
        call reconstruct_half_fixed(odd_oris,odd_stages(:,:,:,3),LAMBDAS(ilambda),2, &
            &HALFSET_MASK_RADIUS,temporary_volume,ignored_niters,data_residual)
        lambda_cross(:,:,:,ilambda) = temporary_volume
        deallocate(temporary_volume)
    enddo

    arm = 0
    call run_arm(truth_volume,SHELL_FULL)
    do i = 1, NITERATIONS
        call run_arm(own_trajectory(:,:,:,i),SHELL_FULL)
    enddo
    do i = 1, NITERATIONS
        call run_arm(cross_trajectory(:,:,:,i),SHELL_FULL)
    enddo
    call run_arm(holdout_volume,SHELL_FULL)
    call run_arm(own_trajectory(:,:,:,1),SHELL_LOW)
    call run_arm(own_trajectory(:,:,:,1),SHELL_MID)
    call run_arm(cross_trajectory(:,:,:,1),SHELL_LOW)
    call run_arm(cross_trajectory(:,:,:,1),SHELL_MID)
    call run_arm(holdout_volume,SHELL_LOW)
    call run_arm(holdout_volume,SHELL_MID)
    call run_arm(own_unmasked(:,:,:,1),SHELL_FULL)
    call run_arm(cross_unmasked(:,:,:,1),SHELL_FULL)
    do ilambda = 1, NLAMBDAS
        call run_arm(lambda_own(:,:,:,ilambda),SHELL_FULL)
        call run_arm(lambda_cross(:,:,:,ilambda),SHELL_FULL)
    enddo

    call assert_true(arm == NARMS, &
        &'reference-bias matrix did not execute every declared arm')
    call print_reference_bias_decision(metrics)
    do i = 1, NARMS
        call assert_metrics_finite(metrics(i),ARM_NAMES(i))
    enddo
    call assert_true(stationary(metrics(1)), &
        &'matched truth control regressed in the reference-bias diagnostic')
    call assert_true(metrics(1)%recovery_rotation_rms < RECOVERY_ROTATION_INITIAL .and. &
        &metrics(1)%recovery_shift_rms < RECOVERY_SHIFT_INITIAL, &
        &'matched truth control lost its local recovery basin')
    write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS: EVIDENCE COMPLETE'

    call orientation%kill()
    call holdout_oris%kill()
    call even_oris%kill()
    call odd_oris%kill()
    call sampler%kill()

contains

    !> Execute and print one predeclared reference, shell, and reconstruction-policy arm.
    subroutine run_arm(reference_volume,shell_range)
        real, intent(in) :: reference_volume(BOX,BOX,BOX)
        integer, intent(in) :: shell_range(2)

        arm = arm+1
        call evaluate_matched_reference_batch(reference_volume,exact_rotmats,probe_observed, &
            &shell_range,metrics(arm))
        metrics(arm)%reference_truth_correlation = &
            &centered_array_correlation(reference_volume,truth_volume)
        call print_reference_arm(ARM_NAMES(arm),shell_range,metrics(arm))
    end subroutine run_arm
end subroutine run_reference_bias_diagnostic

!> Print fit, stationarity, gauge, recovery, and objective evidence for one arm.
subroutine print_reference_arm(name,shell_range,metrics)
    character(len=*), intent(in) :: name
    integer, intent(in) :: shell_range(2)
    type(matched_reference_metrics), intent(in) :: metrics

    write(*,'(a,1x,a,2(1x,i0),6(1x,es14.6))') &
        &'CONTINUOUS_3D_REFERENCE_BIAS fit shells/residual/scale/objective/rot-grad/shift-grad/truth-corr', &
        &trim(name),shell_range,metrics%fitted_residual, &
        &metrics%amplitude_scale,metrics%exact_objective,metrics%rotation_gradient_rms, &
        &metrics%shift_gradient_rms,metrics%reference_truth_correlation
    write(*,'(a,1x,a,4(1x,es14.6),1x,i0)') &
        &'CONTINUOUS_3D_REFERENCE_BIAS exact raw-rot/raw-shift/gauge-rot/gauge-shift/accepted', &
        &trim(name),metrics%exact_rotation_rms,metrics%exact_shift_rms, &
        &metrics%gauge_rotation_rms,metrics%gauge_shift_rms,metrics%exact_accepted
    write(*,'(a,1x,a,6(1x,es14.6),1x,i0)') &
        &'CONTINUOUS_3D_REFERENCE_BIAS recovery rot/shift/exact-obj-before/after/rec-obj-before/after/accepted', &
        &trim(name),metrics%recovery_rotation_rms,metrics%recovery_shift_rms, &
        &metrics%exact_objective_before,metrics%exact_objective_after, &
        &metrics%recovery_objective_before,metrics%recovery_objective_after, &
        &metrics%recovery_accepted
end subroutine print_reference_arm

!> Print independent contribution flags before the primary diagnosis label.
subroutine print_reference_bias_decision(metrics)
    type(matched_reference_metrics), intent(in) :: metrics(NARMS)
    real(dp) :: base_error, cross_error, holdout_error, iteration_best
    real(dp) :: frequency_best, gauge_best, lambda_best, unmasked_best
    logical :: gauge_component, leave_in_component, frequency_component
    logical :: iteration_component, lambda_component, support_component

    base_error = combined_error(metrics(2))
    cross_error = combined_error(metrics(5))
    holdout_error = combined_error(metrics(8))
    gauge_best = min(gauge_error(metrics(2)),gauge_error(metrics(5)),gauge_error(metrics(8)))
    iteration_best = min(combined_error(metrics(3)),combined_error(metrics(4)), &
        &combined_error(metrics(6)),combined_error(metrics(7)))
    frequency_best = min(combined_error(metrics(9)),combined_error(metrics(10)), &
        &combined_error(metrics(11)),combined_error(metrics(12)), &
        &combined_error(metrics(13)),combined_error(metrics(14)))
    unmasked_best = min(combined_error(metrics(15)),combined_error(metrics(16)))
    lambda_best = minval([combined_error(metrics(17)),combined_error(metrics(18)), &
        &combined_error(metrics(19)),combined_error(metrics(20)), &
        &combined_error(metrics(21)),combined_error(metrics(22))])
    gauge_component = gauge_best < 0.5_dp*min(base_error,cross_error,holdout_error)
    leave_in_component = min(cross_error,holdout_error) < 0.75_dp*base_error
    frequency_component = frequency_best < 0.75_dp*min(base_error,cross_error,holdout_error)
    iteration_component = iteration_best < 0.75_dp*min(base_error,cross_error)
    lambda_component = lambda_best < 0.75_dp*min(base_error,cross_error)
    support_component = unmasked_best < 0.75_dp*min(base_error,cross_error)
    write(*,'(a,6(1x,l1))') &
        &'CONTINUOUS_3D_REFERENCE_BIAS components gauge/leave-in/frequency/iteration/lambda/support', &
        &gauge_component,leave_in_component,frequency_component,iteration_component, &
        &lambda_component,support_component
    if( .not. stationary(metrics(1)) )then
        write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS DIAGNOSIS: MATCHED_OPERATOR_REGRESSION'
    else if( stationary(metrics(2)) )then
        write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS DIAGNOSIS: OWN_REFERENCE_STATIONARY'
    else if( stationary(metrics(5)) .or. stationary(metrics(8)) )then
        write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS DIAGNOSIS: LEAVE_IN_REFERENCE_BIAS'
    else if( gauge_component )then
        write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS DIAGNOSIS: GLOBAL_RIGID_GAUGE_COMPONENT'
    else
        write(*,'(a)') 'CONTINUOUS_3D_REFERENCE_BIAS DIAGNOSIS: RECONSTRUCTED_REFERENCE_BIAS'
    endif
end subroutine print_reference_bias_decision

!> Verify the reported gauge residual removes a known common rotation and translation.
subroutine test_gauge_removal(reference)
    real(dp), intent(in) :: reference(3,3,NPROBES)
    real(dp) :: estimate(3,3,NPROBES), shifts(2,NPROBES)
    real(dp) :: identity(3,3), gauge(3,3), translation(3)
    real(dp) :: rotation_residual, shift_residual
    integer :: i

    identity = 0._dp
    identity(1,1) = 1._dp
    identity(2,2) = 1._dp
    identity(3,3) = 1._dp
    gauge = right_increment_rotation(identity,[0.006_dp,-0.004_dp,0.003_dp])
    translation = [0.12_dp,-0.08_dp,0.05_dp]
    do i = 1, NPROBES
        estimate(:,:,i) = matmul(gauge,reference(:,:,i))
        shifts(:,i) = -matmul(reference(1:2,:,i),translation)
    enddo
    call gauge_corrected_errors(reference,estimate,shifts,rotation_residual,shift_residual)
    call assert_true(rotation_residual < 1.e-3_dp .and. shift_residual < 1.e-10_dp, &
        &'reference-bias rigid-gauge control did not remove a known transform')
    write(*,'(a,2(1x,es14.6))') &
        &'CONTINUOUS_3D_REFERENCE_BIAS gauge control rotation/shift residual', &
        &rotation_residual,shift_residual
end subroutine test_gauge_removal

!> Combine rotation and shift after normalization by their frozen stationarity tolerances.
pure real(dp) function combined_error(metrics) result(error)
    type(matched_reference_metrics), intent(in) :: metrics

    error = sqrt((metrics%exact_rotation_rms/ROTATION_TOLERANCE)**2+ &
        &(metrics%exact_shift_rms/SHIFT_TOLERANCE)**2)
end function combined_error

!> Combine gauge-corrected rotation and shift with the same dimensionless scaling.
pure real(dp) function gauge_error(metrics) result(error)
    type(matched_reference_metrics), intent(in) :: metrics

    error = sqrt((metrics%gauge_rotation_rms/ROTATION_TOLERANCE)**2+ &
        &(metrics%gauge_shift_rms/SHIFT_TOLERANCE)**2)
end function gauge_error

!> Apply the frozen exact-pose stationarity tolerances to one diagnostic arm.
pure logical function stationary(metrics) result(is_stationary)
    type(matched_reference_metrics), intent(in) :: metrics

    is_stationary = metrics%exact_rotation_rms <= ROTATION_TOLERANCE .and. &
        &metrics%exact_shift_rms <= SHIFT_TOLERANCE
end function stationary

!> Require finite collect-all evidence without converting an open hypothesis into a gate.
subroutine assert_metrics_finite(metrics,name)
    type(matched_reference_metrics), intent(in) :: metrics
    character(len=*), intent(in) :: name

    call assert_true(all(ieee_is_finite([metrics%reference_truth_correlation, &
        &metrics%fitted_residual,metrics%amplitude_scale, &
        &metrics%exact_objective,metrics%rotation_gradient_rms,metrics%shift_gradient_rms, &
        &metrics%exact_rotation_rms,metrics%exact_shift_rms,metrics%gauge_rotation_rms, &
        &metrics%gauge_shift_rms,metrics%recovery_rotation_rms,metrics%recovery_shift_rms, &
        &metrics%exact_objective_before,metrics%exact_objective_after, &
        &metrics%recovery_objective_before,metrics%recovery_objective_after])), &
        &trim(name)//' reference-bias arm produced non-finite evidence')
end subroutine assert_metrics_finite

end module continuous_3D_pcg_refinement_reference_bias_test
