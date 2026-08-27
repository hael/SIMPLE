module continuous_3D_pcg_reconstruction_noise_gauran_test
use continuous_3D_pcg_reconstruction_noise_test_support, only: dp => NOISE_DP, &
    &population_mean, population_variance, relative_error
use continuous_3D_pcg_reconstruction_test_helpers, only: assert_true, &
    &set_deterministic_seed, BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_image, only: image
implicit none
private
public :: run_gauran_replacement

integer, parameter :: UNIT_NOISE_SEED = 20260812
real, parameter :: SMPD = 1.5
real(dp), parameter :: UNIT_MOMENT_TOLERANCE = 0.05_dp

contains

!> Verify deterministic Gaussian noise statistics used by the test fixtures.
subroutine run_gauran_replacement()
    type(image) :: noise_image
    real, allocatable :: samples(:,:,:)
    real(dp) :: sample_mean, sample_variance

    call noise_image%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
    call set_deterministic_seed(UNIT_NOISE_SEED)
    call noise_image%gauran(0., 1.)
    samples = noise_image%get_rmat()
    sample_mean = population_mean(samples)
    sample_variance = population_variance(samples)

    call assert_true(all(ieee_is_finite(samples)), 'gauran produced a non-finite voxel')
    call assert_true(abs(sample_mean) < UNIT_MOMENT_TOLERANCE, &
        &'unit Gaussian sample mean is outside its predeclared tolerance')
    call assert_true(relative_error(sample_variance, 1._dp) < UNIT_MOMENT_TOLERANCE, &
        &'unit Gaussian sample variance is outside its predeclared tolerance')

    write(*,'(a,i0)') 'CONTINUOUS_3D_PCG_NOISE gauran seed: ', UNIT_NOISE_SEED
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_NOISE gauran mean/variance: ', &
        &sample_mean, sample_variance
    call noise_image%kill()
end subroutine run_gauran_replacement

end module continuous_3D_pcg_reconstruction_noise_gauran_test
