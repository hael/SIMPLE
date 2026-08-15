module continuous_3D_pcg_refinement_halfset_matrix
use continuous_3D_pcg_refinement_halfset_gridding, only: reconstruct_half_conventionally
use continuous_3D_pcg_refinement_halfset_support, only: apply_fixed_support, &
    &array_l2_norm, build_disjoint_half_orientations, build_independent_observations, &
    &calculate_fsc, centered_array_correlation, HALFSET_LAMBDA, HALFSET_MASK_RADIUS, HALFSET_SMPD, &
    &reconstruct_half_pair_fixed, reconstruct_half_trajectory, relative_forward_residual, write_mrc_volume
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp
use simple_oris, only: oris
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_syslib, only: simple_mkdir
implicit none
private

integer, parameter :: NHALF = 48
integer, parameter :: EVEN_NOISE_SEED = 20260815
integer, parameter :: ODD_NOISE_SEED = 20260816
integer, parameter :: NITERS = 6
integer, parameter :: ITERATIONS(NITERS) = [1, 2, 4, 8, 16, 40]
integer, parameter :: NSUPPORT_ITERS = 3
integer, parameter :: SUPPORT_ITERATIONS(NSUPPORT_ITERS) = [4, 8, 40]
integer, parameter :: NLAMBDA = 13
integer, parameter :: N_LOW_SHELLS = 3
integer, parameter :: N_HIGH_SHELLS = 3
real, parameter :: LAMBDAS(NLAMBDA) = [1.e-3, 1.e-2, 1.e-1, 1., 10., 100., &
    &1.e3, 2.e3, 3.e3, 5.e3, 1.e4, 1.e5, 1.e6]
real, parameter :: OBSERVATION_SNR = 0.5
real(dp), parameter :: CLEAN_CORRELATION_MIN = 0.85_dp

public :: run_extended_halfset_matrix

contains

subroutine run_extended_halfset_matrix()
    type(oris) :: even_oris, odd_oris
    type(reconstructor_pcg) :: sampler
    complex, allocatable :: even_clean_planes(:,:,:), even_planes(:,:,:)
    complex, allocatable :: odd_clean_planes(:,:,:), odd_planes(:,:,:)
    real, allocatable :: even_clean_images(:,:,:), even_images(:,:,:), even_noise(:,:,:)
    real, allocatable :: odd_clean_images(:,:,:), odd_images(:,:,:), odd_noise(:,:,:)
    real, allocatable :: clean_even_traj(:,:,:,:), clean_odd_traj(:,:,:,:)
    real, allocatable :: noisy_even_traj(:,:,:,:), noisy_odd_traj(:,:,:,:)
    real, allocatable :: unmasked_clean_even(:,:,:,:), unmasked_clean_odd(:,:,:,:)
    real, allocatable :: unmasked_noisy_even(:,:,:,:), unmasked_noisy_odd(:,:,:,:)
    real, allocatable :: clean_even(:,:,:), clean_odd(:,:,:), noisy_even(:,:,:), noisy_odd(:,:,:)
    real, allocatable :: grid_clean_even(:,:,:), grid_clean_odd(:,:,:)
    real, allocatable :: grid_noisy_even(:,:,:), grid_noisy_odd(:,:,:)
    real, allocatable :: supported_truth(:,:,:), fsc(:), trajectory_fsc_shells(:,:)
    integer, allocatable :: even_ids(:), odd_ids(:), clean_even_niters(:), clean_odd_niters(:)
    integer, allocatable :: noisy_even_niters(:), noisy_odd_niters(:)
    integer, allocatable :: unmasked_clean_even_niters(:), unmasked_clean_odd_niters(:)
    integer, allocatable :: unmasked_noisy_even_niters(:), unmasked_noisy_odd_niters(:)
    real(dp) :: trajectory_corr(4,NITERS), lambda_corr(4,NLAMBDA), grid_corr(4)
    real(dp) :: trajectory_norm(4,NITERS)
    real(dp) :: lambda_norm(4,NLAMBDA), lambda_data_residual(4,NLAMBDA)
    real(dp) :: lambda_raw_error(4,NLAMBDA), lambda_scaled_error(4,NLAMBDA)
    real(dp) :: unmasked_corr(4,NSUPPORT_ITERS)
    real(dp) :: trajectory_fsc(2,NITERS), unmasked_fsc(2,NSUPPORT_ITERS)
    real(dp) :: lambda_fsc(2,NLAMBDA), grid_fsc(2)
    real(dp) :: grid_norm(4), grid_data_residual(4), grid_raw_error(4), grid_scaled_error(4)
    real(dp) :: even_snr, odd_snr, noise_correlation
    real(dp) :: even_pair_residuals(2), odd_pair_residuals(2)
    integer :: best_noisy_lambda_index(2), lambda_niters(4,NLAMBDA)
    integer :: even_pair_niters(2), odd_pair_niters(2), i
    character(len=64) :: volume_dir

    call create_volume_directory(volume_dir)
    call build_disjoint_half_orientations(NHALF, even_oris, odd_oris, even_ids, odd_ids)
    call sampler%new(BOX, HALFSET_SMPD, HALFSET_LAMBDA)
    call build_independent_observations(sampler, even_oris, EVEN_NOISE_SEED, OBSERVATION_SNR, &
        &even_planes, even_clean_planes, even_clean_images, even_images, even_noise, even_snr)
    call build_independent_observations(sampler, odd_oris, ODD_NOISE_SEED, OBSERVATION_SNR, &
        &odd_planes, odd_clean_planes, odd_clean_images, odd_images, odd_noise, odd_snr)
    call sampler%kill()
    call build_truth_volume(supported_truth)
    call apply_fixed_support(supported_truth)
    call write_matrix_volume(supported_truth, volume_dir, 'supported_truth.mrc')
    noise_correlation = centered_array_correlation(even_noise, odd_noise)

    call reconstruct_half_trajectory(even_oris, even_clean_planes, ITERATIONS, HALFSET_MASK_RADIUS, &
        &clean_even_traj, clean_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_clean_planes, ITERATIONS, HALFSET_MASK_RADIUS, &
        &clean_odd_traj, clean_odd_niters)
    call reconstruct_half_trajectory(even_oris, even_planes, ITERATIONS, HALFSET_MASK_RADIUS, &
        &noisy_even_traj, noisy_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_planes, ITERATIONS, HALFSET_MASK_RADIUS, &
        &noisy_odd_traj, noisy_odd_niters)
    do i = 1, NITERS
        trajectory_corr(1,i) = centered_array_correlation(clean_even_traj(:,:,:,i), supported_truth)
        trajectory_corr(2,i) = centered_array_correlation(clean_odd_traj(:,:,:,i), supported_truth)
        trajectory_corr(3,i) = centered_array_correlation(noisy_even_traj(:,:,:,i), supported_truth)
        trajectory_corr(4,i) = centered_array_correlation(noisy_odd_traj(:,:,:,i), supported_truth)
        trajectory_norm(:,i) = [array_l2_norm(clean_even_traj(:,:,:,i)), &
            &array_l2_norm(clean_odd_traj(:,:,:,i)), array_l2_norm(noisy_even_traj(:,:,:,i)), &
            &array_l2_norm(noisy_odd_traj(:,:,:,i))]
        call calculate_fsc(noisy_even_traj(:,:,:,i), noisy_odd_traj(:,:,:,i), fsc)
        if( .not. allocated(trajectory_fsc_shells) ) allocate(trajectory_fsc_shells(size(fsc),NITERS))
        trajectory_fsc_shells(:,i) = fsc
        call summarize_fsc(fsc, trajectory_fsc(:,i))
    enddo
    call write_matrix_volume(noisy_even_traj(:,:,:,3), volume_dir, 'iter004_noisy_even.mrc')
    call write_matrix_volume(noisy_odd_traj(:,:,:,3),  volume_dir, 'iter004_noisy_odd.mrc')
    call write_matrix_volume(0.5*(noisy_even_traj(:,:,:,3)+noisy_odd_traj(:,:,:,3)), &
        &volume_dir, 'iter004_noisy_average.mrc')
    call write_matrix_volume(noisy_even_traj(:,:,:,4), volume_dir, 'iter008_noisy_even.mrc')
    call write_matrix_volume(noisy_odd_traj(:,:,:,4),  volume_dir, 'iter008_noisy_odd.mrc')
    call write_matrix_volume(0.5*(noisy_even_traj(:,:,:,4)+noisy_odd_traj(:,:,:,4)), &
        &volume_dir, 'iter008_noisy_average.mrc')
    call write_matrix_volume(noisy_even_traj(:,:,:,6), volume_dir, 'iter040_noisy_even.mrc')
    call write_matrix_volume(noisy_odd_traj(:,:,:,6),  volume_dir, 'iter040_noisy_odd.mrc')
    call write_matrix_volume(0.5*(noisy_even_traj(:,:,:,6)+noisy_odd_traj(:,:,:,6)), &
        &volume_dir, 'iter040_noisy_average.mrc')

    call reconstruct_half_trajectory(even_oris, even_clean_planes, SUPPORT_ITERATIONS, 0., &
        &unmasked_clean_even, unmasked_clean_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_clean_planes, SUPPORT_ITERATIONS, 0., &
        &unmasked_clean_odd, unmasked_clean_odd_niters)
    call reconstruct_half_trajectory(even_oris, even_planes, SUPPORT_ITERATIONS, 0., &
        &unmasked_noisy_even, unmasked_noisy_even_niters)
    call reconstruct_half_trajectory(odd_oris, odd_planes, SUPPORT_ITERATIONS, 0., &
        &unmasked_noisy_odd, unmasked_noisy_odd_niters)
    do i = 1, NSUPPORT_ITERS
        unmasked_corr(1,i) = centered_array_correlation(unmasked_clean_even(:,:,:,i), supported_truth)
        unmasked_corr(2,i) = centered_array_correlation(unmasked_clean_odd(:,:,:,i), supported_truth)
        unmasked_corr(3,i) = centered_array_correlation(unmasked_noisy_even(:,:,:,i), supported_truth)
        unmasked_corr(4,i) = centered_array_correlation(unmasked_noisy_odd(:,:,:,i), supported_truth)
        call calculate_fsc(unmasked_noisy_even(:,:,:,i), unmasked_noisy_odd(:,:,:,i), fsc)
        call summarize_fsc(fsc, unmasked_fsc(:,i))
    enddo

    do i = 1, NLAMBDA
        call reconstruct_half_pair_fixed(even_oris, even_clean_planes, even_planes, LAMBDAS(i), 40, &
            &clean_even, noisy_even, even_pair_niters, even_pair_residuals)
        call reconstruct_half_pair_fixed(odd_oris, odd_clean_planes, odd_planes, LAMBDAS(i), 40, &
            &clean_odd, noisy_odd, odd_pair_niters, odd_pair_residuals)
        lambda_niters(:,i) = [even_pair_niters(1), odd_pair_niters(1), &
            &even_pair_niters(2), odd_pair_niters(2)]
        lambda_corr(1,i) = centered_array_correlation(clean_even, supported_truth)
        lambda_corr(2,i) = centered_array_correlation(clean_odd, supported_truth)
        lambda_corr(3,i) = centered_array_correlation(noisy_even, supported_truth)
        lambda_corr(4,i) = centered_array_correlation(noisy_odd, supported_truth)
        lambda_norm(:,i) = [array_l2_norm(clean_even), array_l2_norm(clean_odd), &
            &array_l2_norm(noisy_even), array_l2_norm(noisy_odd)]
        lambda_data_residual(:,i) = [even_pair_residuals(1), odd_pair_residuals(1), &
            &even_pair_residuals(2), odd_pair_residuals(2)]
        lambda_raw_error(:,i) = [relative_l2_error(clean_even,supported_truth), &
            &relative_l2_error(clean_odd,supported_truth), relative_l2_error(noisy_even,supported_truth), &
            &relative_l2_error(noisy_odd,supported_truth)]
        lambda_scaled_error(:,i) = [scaled_relative_l2_error(clean_even,supported_truth), &
            &scaled_relative_l2_error(clean_odd,supported_truth), &
            &scaled_relative_l2_error(noisy_even,supported_truth), &
            &scaled_relative_l2_error(noisy_odd,supported_truth)]
        call calculate_fsc(noisy_even, noisy_odd, fsc)
        call summarize_fsc(fsc, lambda_fsc(:,i))
        select case(i)
            case(4)
                call write_matrix_volume(noisy_even, volume_dir, 'lambda001_noisy_even.mrc')
                call write_matrix_volume(noisy_odd,  volume_dir, 'lambda001_noisy_odd.mrc')
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda001_noisy_average.mrc')
            case(5)
                call write_matrix_volume(noisy_even, volume_dir, 'lambda010_noisy_even.mrc')
                call write_matrix_volume(noisy_odd,  volume_dir, 'lambda010_noisy_odd.mrc')
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda010_noisy_average.mrc')
            case(6)
                call write_matrix_volume(noisy_even, volume_dir, 'lambda100_noisy_even.mrc')
                call write_matrix_volume(noisy_odd,  volume_dir, 'lambda100_noisy_odd.mrc')
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda100_noisy_average.mrc')
            case(7)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda1e3_noisy_average.mrc')
            case(8)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda2e3_noisy_average.mrc')
            case(9)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda3e3_noisy_average.mrc')
            case(10)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda5e3_noisy_average.mrc')
            case(11)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda1e4_noisy_average.mrc')
            case(12)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda1e5_noisy_average.mrc')
            case(13)
                call write_matrix_volume(0.5*(noisy_even+noisy_odd), volume_dir, 'lambda1e6_noisy_average.mrc')
        end select
        call assert_true(all(even_pair_niters == 40) .and. all(odd_pair_niters == 40), &
            &'a lambda-sweep solve did not run exactly 40 iterations')
    enddo

    call reconstruct_half_conventionally(even_oris, even_clean_images, grid_clean_even)
    call reconstruct_half_conventionally(odd_oris, odd_clean_images, grid_clean_odd)
    call reconstruct_half_conventionally(even_oris, even_images, grid_noisy_even)
    call reconstruct_half_conventionally(odd_oris, odd_images, grid_noisy_odd)
    call apply_fixed_support(grid_clean_even)
    call apply_fixed_support(grid_clean_odd)
    call apply_fixed_support(grid_noisy_even)
    call apply_fixed_support(grid_noisy_odd)
    call write_matrix_volume(grid_noisy_even, volume_dir, 'conventional_noisy_even.mrc')
    call write_matrix_volume(grid_noisy_odd,  volume_dir, 'conventional_noisy_odd.mrc')
    call write_matrix_volume(0.5*(grid_noisy_even+grid_noisy_odd), volume_dir, 'conventional_noisy_average.mrc')
    grid_corr = [centered_array_correlation(grid_clean_even, supported_truth), &
        &centered_array_correlation(grid_clean_odd, supported_truth), &
        &centered_array_correlation(grid_noisy_even, supported_truth), &
        &centered_array_correlation(grid_noisy_odd, supported_truth)]
    grid_norm = [array_l2_norm(grid_clean_even), array_l2_norm(grid_clean_odd), &
        &array_l2_norm(grid_noisy_even), array_l2_norm(grid_noisy_odd)]
    grid_raw_error = [relative_l2_error(grid_clean_even,supported_truth), &
        &relative_l2_error(grid_clean_odd,supported_truth), relative_l2_error(grid_noisy_even,supported_truth), &
        &relative_l2_error(grid_noisy_odd,supported_truth)]
    grid_scaled_error = [scaled_relative_l2_error(grid_clean_even,supported_truth), &
        &scaled_relative_l2_error(grid_clean_odd,supported_truth), &
        &scaled_relative_l2_error(grid_noisy_even,supported_truth), &
        &scaled_relative_l2_error(grid_noisy_odd,supported_truth)]
    call sampler%new(BOX, HALFSET_SMPD)
    call relative_forward_residual(sampler, even_oris, even_clean_planes, grid_clean_even, grid_data_residual(1))
    call relative_forward_residual(sampler, odd_oris, odd_clean_planes, grid_clean_odd, grid_data_residual(2))
    call relative_forward_residual(sampler, even_oris, even_planes, grid_noisy_even, grid_data_residual(3))
    call relative_forward_residual(sampler, odd_oris, odd_planes, grid_noisy_odd, grid_data_residual(4))
    call sampler%kill()
    call calculate_fsc(grid_noisy_even, grid_noisy_odd, fsc)
    call summarize_fsc(fsc, grid_fsc)
    best_noisy_lambda_index(1) = minloc(lambda_raw_error(3,:),dim=1)
    best_noisy_lambda_index(2) = minloc(lambda_raw_error(4,:),dim=1)

    call assert_true(all(clean_even_niters == ITERATIONS) .and. all(clean_odd_niters == ITERATIONS) .and. &
        &all(noisy_even_niters == ITERATIONS) .and. all(noisy_odd_niters == ITERATIONS), &
        &'a 48-view trajectory did not run its requested iterations')
    call assert_true(all(unmasked_clean_even_niters == SUPPORT_ITERATIONS) .and. &
        &all(unmasked_clean_odd_niters == SUPPORT_ITERATIONS) .and. &
        &all(unmasked_noisy_even_niters == SUPPORT_ITERATIONS) .and. &
        &all(unmasked_noisy_odd_niters == SUPPORT_ITERATIONS), &
        &'an unmasked 48-view trajectory did not run its requested iterations')
    call assert_true(all(ieee_is_finite(trajectory_corr)) .and. all(ieee_is_finite(trajectory_norm)) .and. &
        &all(ieee_is_finite(trajectory_fsc)) .and. all(ieee_is_finite(trajectory_fsc_shells)) .and. &
        &all(ieee_is_finite(lambda_corr)) .and. all(ieee_is_finite(lambda_fsc)) .and. &
        &all(ieee_is_finite(lambda_norm)) .and. all(ieee_is_finite(lambda_data_residual)) .and. &
        &all(ieee_is_finite(lambda_raw_error)) .and. all(ieee_is_finite(lambda_scaled_error)) .and. &
        &all(ieee_is_finite(unmasked_corr)) .and. all(ieee_is_finite(unmasked_fsc)) .and. &
        &all(ieee_is_finite(grid_corr)) .and. all(ieee_is_finite(grid_fsc)) .and. &
        &all(ieee_is_finite(grid_norm)) .and. all(ieee_is_finite(grid_data_residual)) .and. &
        &all(ieee_is_finite(grid_raw_error)) .and. all(ieee_is_finite(grid_scaled_error)), &
        &'the comprehensive half-set matrix contains a non-finite correlation')
    write(*,'(a,i0)') 'CONTINUOUS_3D_MATRIX projections per half: ', NHALF
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_MATRIX realized even/odd SNR: ', even_snr, odd_snr
    write(*,'(a,es14.6)') 'CONTINUOUS_3D_MATRIX noise correlation: ', noise_correlation
    do i = 1, NITERS
        write(*,'(a,i0,a,6(es14.6,1x))') 'CONTINUOUS_3D_MATRIX trajectory iterations ', ITERATIONS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd/FSC-low/FSC-high: ', &
            &trajectory_corr(:,i), trajectory_fsc(:,i)
        write(*,'(a,i0,a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX solution norms iterations ', ITERATIONS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', trajectory_norm(:,i)
        call print_fsc_shells('CONTINUOUS_3D_MATRIX', ITERATIONS(i), trajectory_fsc_shells(:,i))
    enddo
    do i = 1, NSUPPORT_ITERS
        write(*,'(a,i0,a,6(es14.6,1x))') 'CONTINUOUS_3D_MATRIX unmasked iterations ', &
            &SUPPORT_ITERATIONS(i), ' clean-even/clean-odd/noisy-even/noisy-odd/FSC-low/FSC-high: ', &
            &unmasked_corr(:,i), unmasked_fsc(:,i)
    enddo
    do i = 1, NLAMBDA
        write(*,'(a,es12.4,a,6(es14.6,1x))') 'CONTINUOUS_3D_MATRIX lambda ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd/FSC-low/FSC-high: ', &
            &lambda_corr(:,i), lambda_fsc(:,i)
        write(*,'(a,es12.4,a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX lambda norms ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', lambda_norm(:,i)
        write(*,'(a,es12.4,a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX lambda data residuals ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', lambda_data_residual(:,i)
        write(*,'(a,es12.4,a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX lambda raw L2 errors ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', lambda_raw_error(:,i)
        write(*,'(a,es12.4,a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX lambda scaled L2 errors ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', lambda_scaled_error(:,i)
        write(*,'(a,es12.4,a,4(i0,1x))') 'CONTINUOUS_3D_MATRIX lambda iterations ', LAMBDAS(i), &
            &' clean-even/clean-odd/noisy-even/noisy-odd: ', lambda_niters(:,i)
    enddo
    write(*,'(a,6(es14.6,1x))') &
        &'CONTINUOUS_3D_MATRIX conventional clean-even/clean-odd/noisy-even/noisy-odd/FSC-low/FSC-high: ', &
        &grid_corr, grid_fsc
    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX conventional norms: ', grid_norm
    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX conventional data residuals: ', grid_data_residual
    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX conventional raw L2 errors: ', grid_raw_error
    write(*,'(a,4(es14.6,1x))') 'CONTINUOUS_3D_MATRIX conventional scaled L2 errors: ', grid_scaled_error
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_MATRIX best noisy lambda even/odd: ', &
        &LAMBDAS(best_noisy_lambda_index)
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_MATRIX best noisy raw L2 even/odd: ', &
        &lambda_raw_error(3,best_noisy_lambda_index(1)), lambda_raw_error(4,best_noisy_lambda_index(2))

    call assert_true(abs(even_snr-OBSERVATION_SNR)/OBSERVATION_SNR < 0.12_dp .and. &
        &abs(odd_snr-OBSERVATION_SNR)/OBSERVATION_SNR < 0.12_dp, &
        &'a 48-view half set missed its requested SNR')
    call assert_true(abs(noise_correlation) < 0.08_dp, &
        &'the 48-view half-set noises are correlated')
    ! Clean recovery is the gate; noisy correlations below target remain diagnostic evidence.
    call assert_true(maxval(trajectory_corr(1,:)) > CLEAN_CORRELATION_MIN .and. &
        &maxval(trajectory_corr(2,:)) > CLEAN_CORRELATION_MIN, &
        &'the 48-view noiseless PCG controls do not strongly recover supported truth')
    call assert_true(interior_minimum(lambda_raw_error(3,:)) .and. interior_minimum(lambda_raw_error(4,:)), &
        &'the noisy raw-L2 lambda optimum is not bracketed by the expanded sweep')
    call assert_true(lambda_raw_error(3,best_noisy_lambda_index(1)) < grid_raw_error(3) .and. &
        &lambda_raw_error(4,best_noisy_lambda_index(2)) < grid_raw_error(4), &
        &'the best scale-sensitive PCG reconstruction does not beat conventional gridding')

    call even_oris%kill()
    call odd_oris%kill()
end subroutine run_extended_halfset_matrix

subroutine summarize_fsc(fsc, means)
    real, intent(in) :: fsc(:)
    real(dp), intent(out) :: means(2)

    call assert_true(size(fsc) >= N_LOW_SHELLS+N_HIGH_SHELLS, &
        &'a comprehensive-matrix FSC curve has too few shells')
    means(1) = sum(real(fsc(1:N_LOW_SHELLS),dp))/real(N_LOW_SHELLS,dp)
    means(2) = sum(real(fsc(size(fsc)-N_HIGH_SHELLS+1:size(fsc)),dp))/real(N_HIGH_SHELLS,dp)
end subroutine summarize_fsc

subroutine print_fsc_shells(prefix, iteration, fsc)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: iteration
    real, intent(in) :: fsc(:)
    integer :: shell

    do shell = 1, size(fsc)
        write(*,'(a,a,i0,a,i0,a,es14.6)') trim(prefix), ' FSC iterations ', iteration, &
            &' shell ', shell, ': ', fsc(shell)
    enddo
end subroutine print_fsc_shells

pure real(dp) function relative_l2_error(volume, truth) result(error)
    real, intent(in) :: volume(:,:,:), truth(:,:,:)
    real(dp) :: denominator

    denominator = sum(real(truth,dp)**2)
    error = sqrt(sum((real(volume,dp)-real(truth,dp))**2)/max(denominator,tiny(denominator)))
end function relative_l2_error

pure real(dp) function scaled_relative_l2_error(volume, truth) result(error)
    real, intent(in) :: volume(:,:,:), truth(:,:,:)
    real(dp) :: denominator, scale, truth_norm_sq

    denominator = sum(real(volume,dp)**2)
    scale = sum(real(volume,dp)*real(truth,dp))/max(denominator,tiny(denominator))
    truth_norm_sq = sum(real(truth,dp)**2)
    error = sqrt(sum((scale*real(volume,dp)-real(truth,dp))**2)/max(truth_norm_sq,tiny(truth_norm_sq)))
end function scaled_relative_l2_error

pure logical function interior_minimum(values) result(is_interior)
    real(dp), intent(in) :: values(:)
    integer :: location(1)

    location = minloc(values)
    is_interior = location(1) > 1 .and. location(1) < size(values)
end function interior_minimum

subroutine create_volume_directory(volume_dir)
    character(len=*), intent(out) :: volume_dir
    character(len=8) :: date
    character(len=10) :: time

    call date_and_time(date=date, time=time)
    volume_dir = 'continuous_3D_matrix_volumes_'//date//'_'//time(1:6)//time(8:10)
    call simple_mkdir(trim(volume_dir), verbose=.false.)
    write(*,'(a,a)') 'CONTINUOUS_3D_MATRIX volume directory: ', trim(volume_dir)
end subroutine create_volume_directory

subroutine write_matrix_volume(volume, volume_dir, basename)
    real, intent(in) :: volume(:,:,:)
    character(len=*), intent(in) :: volume_dir, basename

    call write_mrc_volume(volume, trim(volume_dir)//'/'//trim(basename))
end subroutine write_matrix_volume

end module continuous_3D_pcg_refinement_halfset_matrix
