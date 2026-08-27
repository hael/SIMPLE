! Test-only fixtures discovered through the existing pose-suite filename glob.
module pose_cont_refinement_test_helpers
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
public :: run_current_executable_case

contains

!> Run one case in a child process so module state cannot leak between groups.
subroutine run_current_executable_case(case_name, groups_run, groups_passed, groups_skipped, failures)
    character(len=*), intent(in) :: case_name
    integer, intent(inout) :: groups_run, groups_passed, groups_skipped, failures
    character(len=:), allocatable :: executable, command
    character(len=1024) :: command_message
    integer :: executable_length, argument_status, command_status, exit_status

    call get_command_argument(0,length=executable_length,status=argument_status)
    if( argument_status /= 0 .or. executable_length < 1 ) &
        &error stop 'could not determine test executable path'
    allocate(character(len=executable_length) :: executable)
    call get_command_argument(0,executable,status=argument_status)
    if( argument_status /= 0 ) error stop 'could not read test executable path'
#if defined(_WIN32)
    command = 'call "'//executable//'" "case='//trim(case_name)//'"'
#else
    command = "'"//executable//"' 'case="//trim(case_name)//"'"
#endif
    write(*,'(/,a)') '>>> TEST ['//trim(case_name)//'] START'
    groups_run = groups_run+1
    call execute_command_line(command,wait=.true.,exitstat=exit_status,cmdstat=command_status, &
        &cmdmsg=command_message)
    if( command_status /= 0 )then
        failures = failures+1
        write(*,'(a,i0,a)') '>>> TEST ['//trim(case_name)//'] FAIL (cmdstat=',command_status,')'
        if( len_trim(command_message) > 0 ) write(*,'(a)') trim(command_message)
    else if( exit_status == 0 )then
        groups_passed = groups_passed+1
        write(*,'(a)') '>>> TEST ['//trim(case_name)//'] PASS'
    else if( exit_status == CASE_SKIP_EXIT_STATUS )then
        groups_skipped = groups_skipped+1
        write(*,'(a)') '>>> TEST ['//trim(case_name)//'] SKIP'
    else
        failures = failures+1
        write(*,'(a,i0,a)') '>>> TEST ['//trim(case_name)//'] FAIL (exitstat=',exit_status,')'
    endif
end subroutine run_current_executable_case

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

end module pose_cont_refinement_test_helpers
