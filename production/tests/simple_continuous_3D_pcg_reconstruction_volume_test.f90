module continuous_3D_pcg_reconstruction_volume_test
use continuous_3D_pcg_reconstruction_test_helpers, only: assert_int_equal, assert_true, &
    &build_truth_volume, BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use iso_fortran_env, only: real64
implicit none
private
public :: run_volume_fixture

integer, parameter :: NBLOBS = 4
integer, parameter :: dp = real64

contains

!> Verify that the deterministic asymmetric truth-volume fixture is reproducible.
subroutine run_volume_fixture()
    real, allocatable :: phantom(:,:,:), repeated_phantom(:,:,:)
    real(dp) :: mean_value, relative_reflection_error, relative_rotation_error
    real(dp) :: sum_squares, total, variance

    call build_truth_volume(phantom)
    call build_truth_volume(repeated_phantom)

    call assert_int_equal(size(phantom, 1), BOX, 'truth-volume x dimension')
    call assert_int_equal(size(phantom, 2), BOX, 'truth-volume y dimension')
    call assert_int_equal(size(phantom, 3), BOX, 'truth-volume z dimension')
    call assert_true(all(ieee_is_finite(phantom)), 'truth volume contains a non-finite voxel')
    call assert_true(minval(phantom) >= 0., 'Gaussian truth volume contains a negative voxel')
    call assert_true(all(phantom == repeated_phantom), &
        &'truth volume is not bit-identical when regenerated from fixed constants')

    total = sum(real(phantom, dp))
    sum_squares = sum(real(phantom, dp)**2)
    mean_value = total / real(size(phantom), dp)
    variance = sum((real(phantom, dp) - mean_value)**2) / real(size(phantom), dp)
    call assert_true(sum_squares > 1._dp, 'truth volume has negligible squared norm')
    call assert_true(variance > 1.e-4_dp, 'truth volume has negligible variance')

    ! These loose fingerprints detect accidental changes to the established
    ! PCG fixture while allowing harmless compiler/libm rounding differences.
    call assert_true(abs(total - 460.9012_dp) < 1.e-2_dp, &
        &'truth-volume sum differs from the established PCG fixture')
    call assert_true(abs(sum_squares - 127.6201_dp) < 1.e-2_dp, &
        &'truth-volume squared norm differs from the established PCG fixture')
    call assert_true(abs(real(maxval(phantom), dp) - 0.910925_dp) < 1.e-4_dp, &
        &'truth-volume maximum differs from the established PCG fixture')

    relative_reflection_error = relative_l2_difference(phantom, phantom(BOX:1:-1,:,:))
    relative_rotation_error = relative_l2_difference(phantom, phantom(BOX:1:-1,BOX:1:-1,:))
    call assert_true(relative_reflection_error > 0.5_dp, &
        &'truth volume is insufficiently asymmetric under x reflection')
    call assert_true(relative_rotation_error > 0.5_dp, &
        &'truth volume is insufficiently asymmetric under a 180-degree z rotation')

    write(*,'(a,i0,a,i0)') 'CONTINUOUS_3D_PCG_VOLUME dimensions: ', BOX, '^3; blobs: ', NBLOBS
    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_PCG_VOLUME sum/sumsq/variance: ', &
        &total, sum_squares, variance
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_VOLUME reflection/rotation relative L2: ', &
        &relative_reflection_error, relative_rotation_error
    write(*,'(a)') 'CONTINUOUS_3D_PCG_VOLUME: PASS'
end subroutine run_volume_fixture

!> Return the relative L2 difference between two volume fixtures.
pure real(dp) function relative_l2_difference(first, second) result(relative_error)
    real, intent(in) :: first(:,:,:), second(:,:,:)
    real(dp) :: denominator

    denominator = sum(real(first, dp)**2)
    if( denominator <= 0._dp )then
        relative_error = huge(relative_error)
    else
        relative_error = sqrt(sum((real(first, dp) - real(second, dp))**2) / denominator)
    endif
end function relative_l2_difference

end module continuous_3D_pcg_reconstruction_volume_test
