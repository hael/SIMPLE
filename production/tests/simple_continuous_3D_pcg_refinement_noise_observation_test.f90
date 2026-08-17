module continuous_3D_pcg_refinement_noise_observation_test
use continuous_3D_pcg_refinement_noise_test_support, only: centered_correlation, &
    &dp => NOISE_DP, population_mean, population_variance, relative_error
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &set_deterministic_seed, BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_image, only: image
use simple_ori, only: ori
use simple_projector, only: projector
implicit none
private
public :: run_observation_noise

integer, parameter :: NOBSERVATIONS = 8
integer, parameter :: OBSERVATION_NOISE_SEED = 20260814
real, parameter :: SMPD = 1.5
real, parameter :: OBSERVATION_SNR = 1.0
real(dp), parameter :: MEAN_TOLERANCE = 0.08_dp
real(dp), parameter :: SNR_TOLERANCE = 0.10_dp
real(dp), parameter :: CORRELATION_TOLERANCE = 0.10_dp
real, parameter :: EULERS(3,NOBSERVATIONS) = reshape([&
    &  0.,  0.,   0., &
    & 30., 45.,  10., &
    & 60., 70.,  20., &
    & 90., 55.,  40., &
    &120., 80.,  60., &
    &180., 40.,  90., &
    &240., 65., 120., &
    &300., 85., 150.], [3,NOBSERVATIONS])

contains

!> Verify independent Fourier-observation noise at a requested signal-to-noise ratio.
subroutine run_observation_noise()
    type(image) :: clean_image, noisy_even_image, noisy_odd_image
    real, allocatable :: clean(:,:,:), noisy_even(:,:,:), noisy_odd(:,:,:)
    real, allocatable :: noise_even(:,:,:), noise_odd(:,:,:)
    real(dp) :: correlation, even_mean_ratio, even_snr, odd_mean_ratio, odd_snr
    real(dp) :: signal_variance

    call build_clean_observations(clean)
    call clean_image%new([BOX,BOX,NOBSERVATIONS], SMPD, wthreads=.false.)
    call clean_image%set_rmat(clean, .false.)
    call noisy_even_image%copy(clean_image)
    call noisy_odd_image%copy(clean_image)
    call set_deterministic_seed(OBSERVATION_NOISE_SEED)
    call noisy_even_image%add_gauran(OBSERVATION_SNR)
    call noisy_odd_image%add_gauran(OBSERVATION_SNR)

    noisy_even = noisy_even_image%get_rmat()
    noisy_odd = noisy_odd_image%get_rmat()
    noise_even = noisy_even - clean
    noise_odd = noisy_odd - clean
    signal_variance = population_variance(clean)
    even_snr = signal_variance / population_variance(noise_even)
    odd_snr = signal_variance / population_variance(noise_odd)
    even_mean_ratio = abs(population_mean(noise_even)) / sqrt(population_variance(noise_even))
    odd_mean_ratio = abs(population_mean(noise_odd)) / sqrt(population_variance(noise_odd))
    correlation = centered_correlation(noise_even, noise_odd)

    call assert_true(signal_variance > tiny(signal_variance), &
        &'clean synthetic observations have zero variance')
    call assert_true(all(ieee_is_finite(noisy_even)) .and. all(ieee_is_finite(noisy_odd)), &
        &'noisy synthetic observation contains a non-finite value')
    call assert_true(any(noise_even /= noise_odd), &
        &'independent observation-noise draws are bit-identical')
    call assert_true(relative_error(even_snr, real(OBSERVATION_SNR,dp)) < SNR_TOLERANCE, &
        &'even observation-noise draw does not realize the requested SNR')
    call assert_true(relative_error(odd_snr, real(OBSERVATION_SNR,dp)) < SNR_TOLERANCE, &
        &'odd observation-noise draw does not realize the requested SNR')
    call assert_true(even_mean_ratio < MEAN_TOLERANCE .and. odd_mean_ratio < MEAN_TOLERANCE, &
        &'observation-noise draw has a non-negligible standardized mean')
    call assert_true(abs(correlation) < CORRELATION_TOLERANCE, &
        &'independent observation-noise draws are excessively correlated')

    write(*,'(a,i0,a,i0,a,f6.3)') 'CONTINUOUS_3D_PCG_NOISE observation seed: ', &
        &OBSERVATION_NOISE_SEED, '; projections: ', NOBSERVATIONS, '; requested SNR: ', OBSERVATION_SNR
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_NOISE observation realized even/odd SNR: ', &
        &even_snr, odd_snr
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_NOISE observation mean ratios/correlation: ', &
        &even_mean_ratio, odd_mean_ratio, correlation

    call clean_image%kill()
    call noisy_even_image%kill()
    call noisy_odd_image%kill()
end subroutine run_observation_noise

!> Generate clean projection observations from the independent truth projector.
subroutine build_clean_observations(observations)
    real, allocatable, intent(out) :: observations(:,:,:)
    type(image) :: projection
    type(ori) :: orientation
    type(projector) :: truth_projector
    real, allocatable :: phantom(:,:,:), projection_values(:,:,:)
    integer :: iobs

    call build_truth_volume(phantom)
    call truth_projector%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
    call truth_projector%set_rmat(phantom, .false.)
    call truth_projector%fft()
    call truth_projector%expand_cmat(BOX)
    call projection%new([BOX,BOX,1], SMPD, wthreads=.false.)
    call orientation%new(is_ptcl=.false.)
    allocate(observations(BOX,BOX,NOBSERVATIONS), source=0.)
    allocate(projection_values(BOX,BOX,1))
    do iobs = 1, NOBSERVATIONS
        call orientation%set_euler(EULERS(:,iobs))
        call truth_projector%fproject_serial(orientation, projection)
        call projection%ifft()
        call projection%get_rmat_sub(projection_values)
        observations(:,:,iobs) = projection_values(:,:,1)
    enddo

    call orientation%kill()
    call projection%kill()
    call truth_projector%kill_expanded()
    call truth_projector%kill()
end subroutine build_clean_observations

end module continuous_3D_pcg_refinement_noise_observation_test
