module continuous_3D_pcg_refinement_noise_volume_test
use continuous_3D_pcg_refinement_noise_test_support, only: centered_correlation, &
    &dp => NOISE_DP, population_mean, population_variance, relative_error
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &set_deterministic_seed, BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_image, only: image
implicit none
private
public :: run_added_volume_noise

integer, parameter :: VOLUME_NOISE_SEED = 20260813
real, parameter :: SMPD = 1.5
real, parameter :: VOLUME_SNR = 0.5
real(dp), parameter :: MOMENT_TOLERANCE = 0.05_dp
real(dp), parameter :: SNR_TOLERANCE = 0.05_dp
real(dp), parameter :: CORRELATION_TOLERANCE = 0.05_dp

contains

!> Verify seeded additive volume noise without modifying the clean fixture.
subroutine run_added_volume_noise()
    type(image) :: clean_image, noisy_even_image, noisy_odd_image, replay_image
    real, allocatable :: clean(:,:,:), noisy_even(:,:,:), noisy_odd(:,:,:), replay(:,:,:)
    real, allocatable :: noise_even(:,:,:), noise_odd(:,:,:)
    real(dp) :: correlation, even_mean_ratio, even_snr, odd_mean_ratio, odd_snr
    real(dp) :: signal_variance

    call build_truth_volume(clean)
    call clean_image%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
    call clean_image%set_rmat(clean, .false.)
    call noisy_even_image%copy(clean_image)
    call noisy_odd_image%copy(clean_image)
    call replay_image%copy(clean_image)
    call set_deterministic_seed(VOLUME_NOISE_SEED)
    call noisy_even_image%add_gauran(VOLUME_SNR)
    call noisy_odd_image%add_gauran(VOLUME_SNR)
    call set_deterministic_seed(VOLUME_NOISE_SEED)
    call replay_image%add_gauran(VOLUME_SNR)

    noisy_even = noisy_even_image%get_rmat()
    noisy_odd = noisy_odd_image%get_rmat()
    replay = replay_image%get_rmat()
    noise_even = noisy_even - clean
    noise_odd = noisy_odd - clean
    signal_variance = population_variance(clean)
    even_snr = signal_variance / population_variance(noise_even)
    odd_snr = signal_variance / population_variance(noise_odd)
    even_mean_ratio = abs(population_mean(noise_even)) / sqrt(population_variance(noise_even))
    odd_mean_ratio = abs(population_mean(noise_odd)) / sqrt(population_variance(noise_odd))
    correlation = centered_correlation(noise_even, noise_odd)

    call assert_true(all(ieee_is_finite(noisy_even)) .and. all(ieee_is_finite(noisy_odd)), &
        &'add_gauran produced a non-finite volume voxel')
    call assert_true(any(noise_even /= 0.) .and. any(noise_odd /= 0.), &
        &'add_gauran did not change the clean truth volume')
    call assert_true(all(noisy_even == replay), &
        &'add_gauran did not reproduce its draw after resetting the deterministic seed')
    call assert_true(any(noise_even /= noise_odd), &
        &'independent volume-noise draws are bit-identical')
    call assert_true(relative_error(even_snr, real(VOLUME_SNR,dp)) < SNR_TOLERANCE, &
        &'even volume-noise draw does not realize the requested SNR')
    call assert_true(relative_error(odd_snr, real(VOLUME_SNR,dp)) < SNR_TOLERANCE, &
        &'odd volume-noise draw does not realize the requested SNR')
    call assert_true(even_mean_ratio < MOMENT_TOLERANCE .and. &
        &odd_mean_ratio < MOMENT_TOLERANCE, &
        &'volume-noise draw has a non-negligible standardized mean')
    call assert_true(abs(correlation) < CORRELATION_TOLERANCE, &
        &'independent volume-noise draws are excessively correlated')

    write(*,'(a,i0,a,f6.3)') 'CONTINUOUS_3D_PCG_NOISE volume seed: ', &
        &VOLUME_NOISE_SEED, '; requested SNR: ', VOLUME_SNR
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_NOISE volume realized even/odd SNR: ', &
        &even_snr, odd_snr
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_NOISE volume mean ratios/correlation: ', &
        &even_mean_ratio, odd_mean_ratio, correlation

    call clean_image%kill()
    call noisy_even_image%kill()
    call noisy_odd_image%kill()
    call replay_image%kill()
end subroutine run_added_volume_noise

end module continuous_3D_pcg_refinement_noise_volume_test
