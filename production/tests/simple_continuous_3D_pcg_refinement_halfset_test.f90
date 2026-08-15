module continuous_3D_pcg_refinement_halfset_test
use continuous_3D_pcg_refinement_halfset_gridding, only: reconstruct_half_conventionally
use continuous_3D_pcg_refinement_halfset_matrix, only: run_extended_halfset_matrix
use continuous_3D_pcg_refinement_halfset_support, only: apply_fixed_support, &
    &array_l2_norm, build_disjoint_half_orientations, build_independent_observations, calculate_fsc, &
    &centered_array_correlation, HALFSET_LAMBDA, HALFSET_MASK_RADIUS, HALFSET_SMPD, reconstruct_half, &
    &reconstruct_half_trajectory
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp
use simple_oris, only: oris
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_halfset_fsc

integer, parameter :: NPROJS_PER_HALF = 24
integer, parameter :: EVEN_NOISE_SEED = 20260815
integer, parameter :: ODD_NOISE_SEED = 20260816
integer, parameter :: N_LOW_SHELLS = 3
integer, parameter :: N_HIGH_SHELLS = 3
integer, parameter :: N_TRAJECTORY_POINTS = 6
integer, parameter :: TRAJECTORY_ITERATIONS(N_TRAJECTORY_POINTS) = [1, 2, 4, 8, 16, 40]
real, parameter :: OBSERVATION_SNR = 0.5
real(dp), parameter :: SNR_RELATIVE_TOL = 0.12_dp
real(dp), parameter :: NOISE_CORRELATION_TOL = 0.08_dp
real(dp), parameter :: CLEAN_CORRELATION_MIN = 0.85_dp
real(dp), parameter :: LOW_FSC_MIN = 0.40_dp
real(dp), parameter :: FSC_LOW_HIGH_GAP_MIN = 0.10_dp
real(dp), parameter :: FSC_BOUND_SLACK = 1.e-5_dp

contains

subroutine run_halfset_fsc()
    type(oris) :: even_oris, odd_oris
    type(reconstructor_pcg) :: sampler
    complex, allocatable :: even_clean_planes(:,:,:), even_planes(:,:,:), even_planes_snapshot(:,:,:)
    complex, allocatable :: odd_clean_planes(:,:,:), odd_planes(:,:,:)
    real, allocatable :: even_clean_images(:,:,:), even_images(:,:,:)
    real, allocatable :: even_clean_volume(:,:,:), odd_clean_volume(:,:,:)
    real, allocatable :: grid_even_clean_volume(:,:,:), grid_even_volume(:,:,:), grid_fsc(:)
    real, allocatable :: even_noise(:,:,:), even_volume(:,:,:), fsc(:)
    real, allocatable :: grid_odd_clean_volume(:,:,:), grid_odd_volume(:,:,:)
    real, allocatable :: odd_clean_images(:,:,:), odd_images(:,:,:)
    real, allocatable :: odd_noise(:,:,:), odd_volume(:,:,:), truth(:,:,:)
    real, allocatable :: trajectory_clean_even(:,:,:,:), trajectory_clean_odd(:,:,:,:)
    real, allocatable :: trajectory_noisy_even(:,:,:,:), trajectory_noisy_odd(:,:,:,:)
    real, allocatable :: trajectory_fsc(:,:)
    integer, allocatable :: even_ids(:), odd_ids(:)
    integer, allocatable :: trajectory_clean_even_niters(:), trajectory_clean_odd_niters(:)
    integer, allocatable :: trajectory_noisy_even_niters(:), trajectory_noisy_odd_niters(:)
    real(dp) :: even_clean_correlation, even_map_correlation, even_snr, high_fsc, low_fsc
    real(dp) :: grid_even_clean_correlation, grid_even_correlation, grid_high_fsc, grid_low_fsc
    real(dp) :: grid_odd_clean_correlation, grid_odd_correlation
    real(dp) :: noise_correlation, odd_clean_correlation, odd_map_correlation, odd_snr
    real(dp) :: trajectory_correlations(4,N_TRAJECTORY_POINTS)
    real(dp) :: trajectory_norms(4,N_TRAJECTORY_POINTS)
    integer :: even_clean_niters, even_niters, i, odd_clean_niters, odd_niters

    call build_disjoint_half_orientations(NPROJS_PER_HALF, even_oris, odd_oris, even_ids, odd_ids)
    call assert_true(all(mod(even_ids,2) == 0) .and. all(mod(odd_ids,2) == 1), &
        &'half-set source ownership does not follow the fixed parity split')
    do i = 1, NPROJS_PER_HALF
        call assert_true(all(even_ids(i) /= odd_ids), &
            &'even and odd projection ownership overlaps')
    enddo

    call sampler%new(BOX, HALFSET_SMPD, HALFSET_LAMBDA)
    call build_independent_observations(sampler, even_oris, EVEN_NOISE_SEED, OBSERVATION_SNR, &
        &even_planes, even_clean_planes, even_clean_images, even_images, even_noise, even_snr)
    allocate(even_planes_snapshot(lbound(even_planes,1):ubound(even_planes,1), &
        &lbound(even_planes,2):ubound(even_planes,2),size(even_planes,3)), source=even_planes)
    call build_independent_observations(sampler, odd_oris, ODD_NOISE_SEED, OBSERVATION_SNR, &
        &odd_planes, odd_clean_planes, odd_clean_images, odd_images, odd_noise, odd_snr)
    call assert_true(all(even_planes == even_planes_snapshot), &
        &'building the odd observations mutated the even-half data')
    call sampler%kill()

    noise_correlation = centered_array_correlation(even_noise, odd_noise)
    call assert_true(relative_error(even_snr, real(OBSERVATION_SNR,dp)) < SNR_RELATIVE_TOL, &
        &'even-half observations do not realize the requested SNR')
    call assert_true(relative_error(odd_snr, real(OBSERVATION_SNR,dp)) < SNR_RELATIVE_TOL, &
        &'odd-half observations do not realize the requested SNR')
    call assert_true(abs(noise_correlation) < NOISE_CORRELATION_TOL, &
        &'even and odd observation noise is excessively correlated')

    ! Each half owns a separate PCG operator, observations, and solution volume.
    call reconstruct_half(even_oris, even_clean_planes, even_clean_volume, even_clean_niters)
    call reconstruct_half(odd_oris, odd_clean_planes, odd_clean_volume, odd_clean_niters)
    call reconstruct_half(even_oris, even_planes, even_volume, even_niters)
    call reconstruct_half(odd_oris, odd_planes, odd_volume, odd_niters)
    call reconstruct_half_trajectory(even_oris, even_clean_planes, TRAJECTORY_ITERATIONS, HALFSET_MASK_RADIUS, &
        &trajectory_clean_even, trajectory_clean_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_clean_planes, TRAJECTORY_ITERATIONS, HALFSET_MASK_RADIUS, &
        &trajectory_clean_odd, trajectory_clean_odd_niters)
    call reconstruct_half_trajectory(even_oris, even_planes, TRAJECTORY_ITERATIONS, HALFSET_MASK_RADIUS, &
        &trajectory_noisy_even, trajectory_noisy_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_planes, TRAJECTORY_ITERATIONS, HALFSET_MASK_RADIUS, &
        &trajectory_noisy_odd, trajectory_noisy_odd_niters)
    ! The conventional gridding control consumes the same real-space observations.
    call reconstruct_half_conventionally(even_oris, even_clean_images, grid_even_clean_volume)
    call reconstruct_half_conventionally(odd_oris, odd_clean_images, grid_odd_clean_volume)
    call reconstruct_half_conventionally(even_oris, even_images, grid_even_volume)
    call reconstruct_half_conventionally(odd_oris, odd_images, grid_odd_volume)
    call assert_true(even_clean_niters > 0 .and. odd_clean_niters > 0, &
        &'a noiseless half-set PCG control ran zero iterations')
    call assert_true(even_niters > 0 .and. odd_niters > 0, &
        &'a half-set PCG solve ran zero iterations')
    call assert_true(all(ieee_is_finite(even_clean_volume)) .and. &
        &all(ieee_is_finite(odd_clean_volume)), &
        &'a noiseless half-set control contains a non-finite value')
    call assert_true(all(ieee_is_finite(even_volume)) .and. all(ieee_is_finite(odd_volume)), &
        &'a reconstructed half volume contains a non-finite value')
    call assert_true(any(even_volume /= odd_volume), &
        &'the independently reconstructed half maps are identical')
    call assert_true(all(trajectory_clean_even(:,:,:,N_TRAJECTORY_POINTS) == even_clean_volume) .and. &
        &all(trajectory_clean_odd(:,:,:,N_TRAJECTORY_POINTS) == odd_clean_volume) .and. &
        &all(trajectory_noisy_even(:,:,:,N_TRAJECTORY_POINTS) == even_volume) .and. &
        &all(trajectory_noisy_odd(:,:,:,N_TRAJECTORY_POINTS) == odd_volume), &
        &'the forced 40-iteration trajectory does not reproduce the baseline PCG solve')

    call build_truth_volume(truth)
    ! The constrained solve x = P u must be assessed against P x_truth.
    call apply_fixed_support(truth)
    call apply_fixed_support(grid_even_clean_volume)
    call apply_fixed_support(grid_odd_clean_volume)
    call apply_fixed_support(grid_even_volume)
    call apply_fixed_support(grid_odd_volume)
    even_clean_correlation = centered_array_correlation(even_clean_volume, truth)
    odd_clean_correlation = centered_array_correlation(odd_clean_volume, truth)
    even_map_correlation = centered_array_correlation(even_volume, truth)
    odd_map_correlation = centered_array_correlation(odd_volume, truth)
    grid_even_clean_correlation = centered_array_correlation(grid_even_clean_volume, truth)
    grid_odd_clean_correlation = centered_array_correlation(grid_odd_clean_volume, truth)
    grid_even_correlation = centered_array_correlation(grid_even_volume, truth)
    grid_odd_correlation = centered_array_correlation(grid_odd_volume, truth)
    do i = 1, N_TRAJECTORY_POINTS
        trajectory_correlations(1,i) = centered_array_correlation(trajectory_clean_even(:,:,:,i), truth)
        trajectory_correlations(2,i) = centered_array_correlation(trajectory_clean_odd(:,:,:,i), truth)
        trajectory_correlations(3,i) = centered_array_correlation(trajectory_noisy_even(:,:,:,i), truth)
        trajectory_correlations(4,i) = centered_array_correlation(trajectory_noisy_odd(:,:,:,i), truth)
        trajectory_norms(:,i) = [array_l2_norm(trajectory_clean_even(:,:,:,i)), &
            &array_l2_norm(trajectory_clean_odd(:,:,:,i)), array_l2_norm(trajectory_noisy_even(:,:,:,i)), &
            &array_l2_norm(trajectory_noisy_odd(:,:,:,i))]
        call calculate_fsc(trajectory_noisy_even(:,:,:,i), trajectory_noisy_odd(:,:,:,i), fsc)
        if( .not. allocated(trajectory_fsc) ) allocate(trajectory_fsc(size(fsc),N_TRAJECTORY_POINTS))
        trajectory_fsc(:,i) = fsc
    enddo

    call calculate_fsc(even_volume, odd_volume, fsc)
    call calculate_fsc(grid_even_volume, grid_odd_volume, grid_fsc)
    call assert_true(size(fsc) >= N_LOW_SHELLS + N_HIGH_SHELLS, &
        &'FSC curve has too few shells')
    call assert_true(size(grid_fsc) >= N_LOW_SHELLS + N_HIGH_SHELLS, &
        &'the conventional gridding FSC curve has too few shells')
    low_fsc = sum(real(fsc(1:N_LOW_SHELLS),dp)) / real(N_LOW_SHELLS,dp)
    high_fsc = sum(real(fsc(size(fsc)-N_HIGH_SHELLS+1:size(fsc)),dp)) / real(N_HIGH_SHELLS,dp)
    grid_low_fsc = sum(real(grid_fsc(1:N_LOW_SHELLS),dp)) / real(N_LOW_SHELLS,dp)
    grid_high_fsc = sum(real(grid_fsc(size(grid_fsc)-N_HIGH_SHELLS+1:size(grid_fsc)),dp)) / &
        &real(N_HIGH_SHELLS,dp)

    write(*,'(a,2(i0,1x))') 'CONTINUOUS_3D_PCG_HALFSET even/odd projections: ', &
        &NPROJS_PER_HALF, NPROJS_PER_HALF
    write(*,'(a,f6.2)') 'CONTINUOUS_3D_PCG_HALFSET fixed support radius: ', HALFSET_MASK_RADIUS
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_HALFSET realized even/odd SNR: ', even_snr, odd_snr
    write(*,'(a,es14.6)') 'CONTINUOUS_3D_PCG_HALFSET noise correlation: ', noise_correlation
    write(*,'(a,2(i0,1x))') 'CONTINUOUS_3D_PCG_HALFSET PCG even/odd iterations: ', even_niters, odd_niters
    write(*,'(a,2(i0,1x))') 'CONTINUOUS_3D_PCG_HALFSET clean PCG even/odd iterations: ', &
        &even_clean_niters, odd_clean_niters
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_HALFSET clean supported-truth correlations: ', &
        &even_clean_correlation, odd_clean_correlation
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_HALFSET supported-truth correlations: ', &
        &even_map_correlation, odd_map_correlation
    do i = 1, N_TRAJECTORY_POINTS
        write(*,'(a,i0,a,4(es14.6,1x))') 'CONTINUOUS_3D_PCG_TRAJECTORY iterations ', &
            &TRAJECTORY_ITERATIONS(i), ' clean-even/clean-odd/noisy-even/noisy-odd: ', &
            &trajectory_correlations(:,i)
        write(*,'(a,i0,a,4(es14.6,1x))') 'CONTINUOUS_3D_PCG_TRAJECTORY solution norms iterations ', &
            &TRAJECTORY_ITERATIONS(i), ' clean-even/clean-odd/noisy-even/noisy-odd: ', trajectory_norms(:,i)
        call print_trajectory_fsc(TRAJECTORY_ITERATIONS(i), trajectory_fsc(:,i))
    enddo
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_GRID_HALFSET clean supported-truth correlations: ', &
        &grid_even_clean_correlation, grid_odd_clean_correlation
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_GRID_HALFSET supported-truth correlations: ', &
        &grid_even_correlation, grid_odd_correlation
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_GRID_HALFSET low/high FSC means: ', &
        &grid_low_fsc, grid_high_fsc
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_HALFSET low/high FSC means: ', low_fsc, high_fsc
    do i = 1, size(fsc)
        write(*,'(a,i0,a,es14.6)') 'CONTINUOUS_3D_PCG_HALFSET FSC shell ', i, ': ', fsc(i)
    enddo
    do i = 1, size(grid_fsc)
        write(*,'(a,i0,a,es14.6)') 'CONTINUOUS_3D_GRID_HALFSET FSC shell ', i, ': ', grid_fsc(i)
    enddo

    call assert_true(all(ieee_is_finite(grid_even_clean_volume)) .and. &
        &all(ieee_is_finite(grid_odd_clean_volume)) .and. all(ieee_is_finite(grid_even_volume)) .and. &
        &all(ieee_is_finite(grid_odd_volume)), &
        &'the conventional gridding control contains a non-finite value')
    call assert_true(any(grid_even_volume /= grid_odd_volume), &
        &'the conventional independently reconstructed half maps are identical')
    call assert_true(all(ieee_is_finite(grid_fsc)), &
        &'the conventional gridding FSC contains a non-finite value')
    call assert_true(all(trajectory_clean_even_niters == TRAJECTORY_ITERATIONS) .and. &
        &all(trajectory_clean_odd_niters == TRAJECTORY_ITERATIONS) .and. &
        &all(trajectory_noisy_even_niters == TRAJECTORY_ITERATIONS) .and. &
        &all(trajectory_noisy_odd_niters == TRAJECTORY_ITERATIONS), &
        &'a fixed-count PCG trajectory solve did not run its requested iterations')
    call assert_true(all(ieee_is_finite(trajectory_correlations)) .and. &
        &all(ieee_is_finite(trajectory_norms)) .and. all(ieee_is_finite(trajectory_fsc)), &
        &'a fixed-count PCG trajectory diagnostic is non-finite')

    ! One extended run tests angular coverage, iteration regularization, lambda, and gridding.
    call run_extended_halfset_matrix()

    ! Stage 0 gates noisy interpretation on strong clean recovery, not an assumed noisy target.
    call assert_true(maxval(trajectory_correlations(1,:)) > CLEAN_CORRELATION_MIN .and. &
        &maxval(trajectory_correlations(2,:)) > CLEAN_CORRELATION_MIN, &
        &'the noiseless PCG controls do not strongly recover supported truth')
    call assert_true(all(ieee_is_finite(fsc)), &
        &'FSC curve contains a non-finite value')
    call assert_true(all(real(fsc,dp) >= -1._dp-FSC_BOUND_SLACK) .and. &
        &all(real(fsc,dp) <= 1._dp+FSC_BOUND_SLACK), &
        &'FSC curve lies outside [-1,1]')
    call assert_true(low_fsc > LOW_FSC_MIN, &
        &'independent half maps lack low-frequency agreement')
    call assert_true(low_fsc-high_fsc > FSC_LOW_HIGH_GAP_MIN, &
        &'synthetic FSC does not degrade from low to high frequency')
    write(*,'(a)') 'CONTINUOUS_3D_PCG_HALFSET_FSC: PASS'

    call even_oris%kill()
    call odd_oris%kill()
end subroutine run_halfset_fsc

subroutine print_trajectory_fsc(iteration, fsc)
    integer, intent(in) :: iteration
    real, intent(in) :: fsc(:)
    integer :: shell

    do shell = 1, size(fsc)
        write(*,'(a,i0,a,i0,a,es14.6)') 'CONTINUOUS_3D_PCG_TRAJECTORY FSC iterations ', iteration, &
            &' shell ', shell, ': ', fsc(shell)
    enddo
end subroutine print_trajectory_fsc

pure real(dp) function relative_error(actual, expected) result(error)
    real(dp), intent(in) :: actual, expected
    error = abs(actual-expected) / max(abs(expected), tiny(expected))
end function relative_error

end module continuous_3D_pcg_refinement_halfset_test
