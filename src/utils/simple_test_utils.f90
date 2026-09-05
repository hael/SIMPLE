!@descr: reusable assertion, suite tracking, and reporting utilities for tests
module simple_test_utils
use, intrinsic :: iso_fortran_env, only: output_unit
use simple_string, only: string
use simple_defs,   only: dp, STDLEN
implicit none
private

public :: assert_true, assert_int, assert_real, assert_char, assert_string_eq, assert_double, assert_false
public :: begin_test_suite, end_test_suite, reset_test_report, report_summary
public :: tests_run, tests_failed

type :: test_suite_result
    character(len=STDLEN) :: name = ''
    integer :: checks   = 0
    integer :: failures = 0
    logical :: completed = .false.
end type test_suite_result

type :: test_failure
    character(len=STDLEN) :: suite   = ''
    character(len=STDLEN) :: message = ''
end type test_failure

integer :: tests_run       = 0
integer :: tests_failed    = 0
integer :: active_suite    = 0
character(len=STDLEN) :: report_path = ''
type(test_suite_result), allocatable :: suites(:)
type(test_failure),      allocatable :: failures(:)

contains

    subroutine reset_test_report( filename )
        character(len=*), optional, intent(in) :: filename
        tests_run    = 0
        tests_failed = 0
        active_suite = 0
        report_path  = ''
        if( allocated(suites)   ) deallocate(suites)
        if( allocated(failures) ) deallocate(failures)
        if( present(filename) ) report_path = trim(filename)
        call write_report_file
    end subroutine reset_test_report

    subroutine begin_test_suite( name )
        character(len=*), intent(in) :: name
        type(test_suite_result) :: suite
        if( active_suite > 0 ) call end_test_suite
        suite%name = trim(name)
        if( allocated(suites) )then
            suites = [suites, suite]
        else
            allocate(suites(1))
            suites(1) = suite
        endif
        active_suite = size(suites)
        write(output_unit,'(A)') '---- TEST SUITE: '//trim(suite%name)//' ----'
        call write_report_file
    end subroutine begin_test_suite

    subroutine end_test_suite
        character(len=STDLEN) :: name
        integer :: checks, nfailed
        if( active_suite == 0 ) return
        name    = suites(active_suite)%name
        checks  = suites(active_suite)%checks
        nfailed = suites(active_suite)%failures
        if( nfailed == 0 )then
            if( checks == 0 )then
                write(output_unit,'(A)') 'PASS: '//trim(name)//' completed'
            else
                write(output_unit,'(A,I0,A)') 'PASS: '//trim(name)//' (', checks, ' checks)'
            endif
        else
            write(output_unit,'(A,I0,A,I0,A)') 'FAIL: '//trim(name)//' (', nfailed, ' of ', checks, ' checks failed)'
        endif
        suites(active_suite)%completed = .true.
        active_suite = 0
        call write_report_file
    end subroutine end_test_suite

    subroutine assert_true( cond, msg )
        logical,          intent(in) :: cond
        character(len=*), intent(in) :: msg
        call record_assertion(cond, msg)
    end subroutine assert_true

    subroutine assert_int( expected, actual, msg )
        integer,          intent(in) :: expected, actual
        character(len=*), intent(in) :: msg
        character(len=STDLEN) :: detail
        if( expected == actual )then
            call record_assertion(.true., msg)
        else
            write(detail,'(A,I0,A,I0)') trim(msg)//' expected=', expected, ' actual=', actual
            call record_assertion(.false., detail)
        endif
    end subroutine assert_int

    subroutine assert_real( expected, actual, tol, msg )
        real,             intent(in) :: expected, actual, tol
        character(len=*), intent(in) :: msg
        character(len=STDLEN) :: detail
        if( abs(expected - actual) <= tol )then
            call record_assertion(.true., msg)
        else
            write(detail,'(A,G0,A,G0,A,G0)') trim(msg)//' expected=', expected, ' actual=', actual, ' tolerance=', tol
            call record_assertion(.false., detail)
        endif
    end subroutine assert_real

    subroutine assert_char( expected, actual, msg )
        character(len=*), intent(in) :: expected, actual, msg
        character(len=STDLEN) :: detail
        if( trim(expected) == trim(actual) )then
            call record_assertion(.true., msg)
        else
            detail = trim(msg)//' expected="'//trim(expected)//'" actual="'//trim(actual)//'"'
            call record_assertion(.false., detail)
        endif
    end subroutine assert_char

    subroutine assert_string_eq( expected, actual, msg )
        character(len=*), intent(in) :: expected
        type(string),     intent(in) :: actual
        character(len=*), intent(in) :: msg
        character(len=:), allocatable :: actual_char
        if( actual%is_allocated() )then
            actual_char = actual%to_char()
        else
            actual_char = ''
        endif
        call assert_char(expected, actual_char, msg)
    end subroutine assert_string_eq

    subroutine assert_double( expected, actual, msg, ulp_tol )
        real(dp),         intent(in) :: expected, actual
        character(len=*), intent(in) :: msg
        real(dp), optional, intent(in) :: ulp_tol
        real(dp) :: tol
        character(len=STDLEN) :: detail
        if( present(ulp_tol) )then
            tol = ulp_tol * epsilon(1.0_dp)
        else
            tol = 10.0_dp * epsilon(1.0_dp)
        endif
        if( abs(expected - actual) <= tol )then
            call record_assertion(.true., msg)
        else
            write(detail,'(A,ES24.16,A,ES24.16,A,ES14.6,A,ES14.6)') trim(msg)//' expected=', expected, &
                &' actual=', actual, ' difference=', abs(expected - actual), ' tolerance=', tol
            call record_assertion(.false., detail)
        endif
    end subroutine assert_double

    subroutine assert_false( condition, message )
        logical,                    intent(in) :: condition
        character(len=*), optional, intent(in) :: message
        if( present(message) )then
            call record_assertion(.not. condition, message)
        else
            call record_assertion(.not. condition, 'assert_false condition was true')
        endif
    end subroutine assert_false

    subroutine record_assertion( passed, message )
        logical,          intent(in) :: passed
        character(len=*), intent(in) :: message
        type(test_failure) :: failure
        tests_run = tests_run + 1
        if( active_suite > 0 ) suites(active_suite)%checks = suites(active_suite)%checks + 1
        if( passed ) return
        tests_failed = tests_failed + 1
        failure%message = trim(message)
        if( active_suite > 0 )then
            suites(active_suite)%failures = suites(active_suite)%failures + 1
            failure%suite = suites(active_suite)%name
        else
            failure%suite = 'unscoped'
        endif
        if( allocated(failures) )then
            failures = [failures, failure]
        else
            allocate(failures(1))
            failures(1) = failure
        endif
        write(output_unit,'(A)') 'FAIL: '//trim(failure%message)
        call write_report_file
    end subroutine record_assertion

    subroutine report_summary( filename, failed )
        character(len=*), optional, intent(in)  :: filename
        logical,          optional, intent(out) :: failed
        integer :: report_unit, io_stat
        logical :: report_is_open
        if( active_suite > 0 ) call end_test_suite
        if( present(filename) ) report_path = trim(filename)
        call write_summary(output_unit, list_groups=.false.)
        report_is_open = .false.
        if( len_trim(report_path) > 0 )then
            open(newunit=report_unit, file=trim(report_path), status='replace', action='write', &
                &iostat=io_stat)
            report_is_open = io_stat == 0
            if( report_is_open )then
                call write_summary(report_unit, list_groups=.true.)
                close(report_unit)
                report_is_open = .false.
                write(output_unit,'(A)') 'Test report: '//trim(report_path)
            else
                write(output_unit,'(A)') 'WARNING: could not write test report: '//trim(report_path)
            endif
        endif
        if( present(failed) ) failed = tests_failed > 0
    end subroutine report_summary

    subroutine write_report_file
        integer :: report_unit, io_stat
        logical :: report_is_open
        if( len_trim(report_path) == 0 ) return
        report_is_open = .false.
        open(newunit=report_unit, file=trim(report_path), status='replace', action='write', iostat=io_stat)
        report_is_open = io_stat == 0
        if( report_is_open )then
            call write_summary(report_unit, list_groups=.true.)
            close(report_unit)
            report_is_open = .false.
        endif
    end subroutine write_report_file

    subroutine write_summary( unit, list_groups )
        integer, intent(in) :: unit
        logical, intent(in) :: list_groups
        integer :: i, nsuites, passed_suites, failed_suites, incomplete_suites
        character(len=10) :: suite_status
        nsuites = 0
        passed_suites = 0
        failed_suites = 0
        incomplete_suites = 0
        if( allocated(suites) )then
            nsuites = size(suites)
            passed_suites = count(suites%completed .and. suites%failures == 0)
            failed_suites = count(suites%completed .and. suites%failures > 0)
            incomplete_suites = nsuites - passed_suites - failed_suites
        endif
        write(unit,'(A)') ''
        write(unit,'(A)') '========== SIMPLE TEST SUMMARY =========='
        if( active_suite > 0 )then
            write(unit,'(A)') 'RESULT: INCOMPLETE'
            write(unit,'(A)') 'Active suite: '//trim(suites(active_suite)%name)
        else if( tests_failed == 0 )then
            write(unit,'(A)') 'RESULT: PASS'
        else
            write(unit,'(A)') 'RESULT: FAIL'
        endif
        write(unit,'(A,I0,A,I0,A,I0,A,I0,A)') 'Test groups: ', nsuites, ' run; ', &
            &passed_suites, ' passed; ', failed_suites, ' failed; ', incomplete_suites, ' incomplete.'
        write(unit,'(A,I0,A,I0,A,I0,A)') 'Assertions: ', tests_run, ' evaluated; ', &
            &tests_run - tests_failed, ' passed; ', tests_failed, ' failed.'
        if( list_groups .and. allocated(suites) )then
            write(unit,'(A)') 'Test group details:'
            do i = 1, size(suites)
                if( .not. suites(i)%completed )then
                    suite_status = 'INCOMPLETE'
                else if( suites(i)%failures > 0 )then
                    suite_status = 'FAIL'
                else
                    suite_status = 'PASS'
                endif
                write(unit,'(A,I0,A,I0,A)') '  '//trim(suite_status)//': '//trim(suites(i)%name)// &
                    &' - ', suites(i)%checks, ' assertions; ', suites(i)%failures, ' failed.'
            enddo
        endif
        if( allocated(failures) )then
            write(unit,'(A)') 'Failures:'
            do i = 1, size(failures)
                write(unit,'(A)') '  ['//trim(failures(i)%suite)//'] '//trim(failures(i)%message)
            enddo
        endif
        write(unit,'(A)') '========================================='
    end subroutine write_summary

end module simple_test_utils
