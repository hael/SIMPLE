module continuous_3D_pcg_refinement_test_helpers
use iso_fortran_env, only: int64
implicit none
private

integer, parameter, public :: CASE_SKIP_EXIT_STATUS = 77
integer, parameter, public :: TRUTH_VOLUME_BOX = 24

integer, parameter :: TRUTH_VOLUME_BLOBS = 4
real, parameter :: TRUTH_VOLUME_CENTRES(3,TRUTH_VOLUME_BLOBS) = reshape([&
    &-5.0, -3.0,  2.0, &
    & 4.0,  5.0, -3.0, &
    & 0.0, -6.0, -5.0, &
    & 3.0, -2.0,  6.0], [3,TRUTH_VOLUME_BLOBS])
real, parameter :: TRUTH_VOLUME_SIGMAS(TRUTH_VOLUME_BLOBS) = [2.0, 2.5, 1.8, 2.2]
real, parameter :: TRUTH_VOLUME_AMPLITUDES(TRUTH_VOLUME_BLOBS) = [1.0, 0.8, 0.6, 0.5]

public :: assert_true
public :: assert_int_equal
public :: assert_real_close
public :: build_truth_volume
public :: set_deterministic_seed
public :: skip_unimplemented_case

contains

!> Stop the selected child case when a logical contract is false.
subroutine assert_true(condition, message)
    logical,          intent(in) :: condition
    character(len=*), intent(in) :: message

    if( .not. condition ) error stop trim(message)
end subroutine assert_true

!> Stop the selected child case when two integer values differ.
subroutine assert_int_equal(actual, expected, message)
    integer,          intent(in) :: actual, expected
    character(len=*), intent(in) :: message
    character(len=256) :: detail

    if( actual == expected ) return
    write(detail,'(a,i0,a,i0)') trim(message)//': actual=', actual, ', expected=', expected
    error stop trim(detail)
end subroutine assert_int_equal

!> Stop the selected child case when an absolute real tolerance is exceeded.
subroutine assert_real_close(actual, expected, tolerance, message)
    use ieee_arithmetic, only: ieee_is_finite
    real,             intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    character(len=256) :: detail

    if( tolerance < 0. ) error stop 'assert_real_close requires a nonnegative tolerance'
    if( ieee_is_finite(actual) .and. ieee_is_finite(expected) )then
        if( abs(actual - expected) <= tolerance ) return
    endif
    write(detail,'(a,es14.6,a,es14.6,a,es14.6)') trim(message)//': actual=', actual, &
        &', expected=', expected, ', tolerance=', tolerance
    error stop trim(detail)
end subroutine assert_real_close

!> Expand one scalar seed into the compiler's full deterministic random seed.
subroutine set_deterministic_seed(base_seed)
    integer, intent(in) :: base_seed
    integer, allocatable :: seed(:)
    integer(int64) :: candidate, modulus
    integer :: i, seed_size

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    modulus = int(huge(0), int64) - 1_int64
    do i = 1, seed_size
        candidate = int(base_seed, int64) + 104729_int64 * int(i - 1, int64)
        seed(i) = int(modulo(candidate, modulus)) + 1
    enddo
    call random_seed(put=seed)
    deallocate(seed)
end subroutine set_deterministic_seed

!> Build the asymmetric deterministic 3-D truth used by all component tests.
subroutine build_truth_volume(volume)
    real, allocatable, intent(out) :: volume(:,:,:)
    real :: centre, dx, dy, dz
    integer :: blob, i, j, k

    allocate(volume(TRUTH_VOLUME_BOX,TRUTH_VOLUME_BOX,TRUTH_VOLUME_BOX), source=0.)
    centre = real(TRUTH_VOLUME_BOX) / 2. + 0.5
    do k = 1, TRUTH_VOLUME_BOX
        do j = 1, TRUTH_VOLUME_BOX
            do i = 1, TRUTH_VOLUME_BOX
                do blob = 1, TRUTH_VOLUME_BLOBS
                    dx = real(i) - centre - TRUTH_VOLUME_CENTRES(1,blob)
                    dy = real(j) - centre - TRUTH_VOLUME_CENTRES(2,blob)
                    dz = real(k) - centre - TRUTH_VOLUME_CENTRES(3,blob)
                    volume(i,j,k) = volume(i,j,k) + TRUTH_VOLUME_AMPLITUDES(blob) * &
                        &exp(-(dx*dx + dy*dy + dz*dz) / (2. * TRUTH_VOLUME_SIGMAS(blob)**2))
                enddo
            enddo
        enddo
    enddo
end subroutine build_truth_volume

!> Mark an explicitly deferred child case with the suite's skip exit status.
subroutine skip_unimplemented_case(case_name, planned_stage)
    character(len=*), intent(in) :: case_name, planned_stage

    write(*,'(a)') 'CONTINUOUS_3D_PCG_CASE ['//trim(case_name)//']: SKIP'
    write(*,'(a)') '  Planned implementation: '//trim(planned_stage)
    stop CASE_SKIP_EXIT_STATUS
end subroutine skip_unimplemented_case

end module continuous_3D_pcg_refinement_test_helpers
