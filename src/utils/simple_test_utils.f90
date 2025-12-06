module simple_test_utils
use simple_string, only: string
use simple_defs
implicit none
private
public  :: assert_true, assert_int, assert_real, assert_char, assert_string_eq, assert_double, assert_false
public  :: tests_run, tests_failed, report_summary
integer :: tests_run    = 0
integer :: tests_failed = 0

contains

    subroutine assert_true(cond, msg)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        tests_run = tests_run + 1
        if (.not. cond) then
            tests_failed = tests_failed + 1
            write(*,*) 'FAIL: ', trim(msg)
        end if
    end subroutine assert_true

    subroutine assert_int(expected, actual, msg)
        integer, intent(in) :: expected, actual
        character(*), intent(in) :: msg
        tests_run = tests_run + 1
        if (expected /= actual) then
            tests_failed = tests_failed + 1
            write(*,*) 'FAIL: ', trim(msg), ' expected=', expected, ' actual=', actual
        end if
    end subroutine assert_int

    subroutine assert_real(expected, actual, tol, msg)
        real, intent(in) :: expected, actual, tol
        character(*), intent(in) :: msg
        tests_run = tests_run + 1
        if (abs(expected - actual) > tol) then
            tests_failed = tests_failed + 1
            write(*,*) 'FAIL: ', trim(msg), ' expected=', expected, ' actual=', actual
        end if
    end subroutine assert_real

    subroutine assert_char(expected, actual, msg)
        character(*), intent(in) :: expected, actual
        character(*), intent(in) :: msg
        tests_run = tests_run + 1
        if (trim(expected) /= trim(actual)) then
            tests_failed = tests_failed + 1
            write(*,*) 'FAIL: ', trim(msg), ' expected="', trim(expected), '" actual="', trim(actual), '"'
        end if
    end subroutine assert_char

    subroutine assert_string_eq(expected, actual, msg)
        character(len=*), intent(in) :: expected
        type(string),     intent(in) :: actual
        character(len=*), intent(in) :: msg
        character(len=:), allocatable :: act
        if (actual%is_allocated()) then
            act = actual%to_char()
        else
            act = ''   ! treat unallocated as empty
        end if
        call assert_char(expected, act, msg)
    end subroutine assert_string_eq

    ! dp-specific version: assert_double
    subroutine assert_double(expected, actual, msg, ulp_tol)
        real(dp), intent(in)           :: expected
        real(dp), intent(in)           :: actual
        character(len=*), intent(in)   :: msg
        real(dp), intent(in), optional :: ulp_tol
        real(dp) :: tol
        tests_run = tests_run + 1
        ! interpret ulp_tol as a multiplier of machine epsilon
        if (present(ulp_tol)) then
            tol = ulp_tol * epsilon(1.0_dp)
        else
            tol = 10.0_dp * epsilon(1.0_dp)
        end if
        if (abs(expected - actual) > tol) then
            tests_failed = tests_failed + 1
            write(*,'(A,ES24.16,A,ES24.16)') 'ASSERT_DOUBLE FAILED: ', expected, ' /= ', actual
            write(*,'(A,ES14.6)') '  |Δ| = ', abs(expected - actual)
            write(*,'(A,ES14.6)') '  tol = ', tol
            write(*,'(A)') '  '//trim(msg)
        end if
    end subroutine assert_double

    subroutine assert_false(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in), optional :: message
        tests_run = tests_run + 1
        if (condition) then
            tests_failed = tests_failed + 1
            print *, 'Assertion (assert_false) failed!'
            if (present(message)) print *, '  Message: ', message
        end if
    end subroutine assert_false

    subroutine report_summary()
        if (tests_failed == 0) then
            write(*,*) 'ALL TESTS PASSED (', tests_run, ' tests).'
        else
            write(*,*) tests_failed, ' TEST(S) FAILED out of ', tests_run
        end if
    end subroutine report_summary

end module simple_test_utils
