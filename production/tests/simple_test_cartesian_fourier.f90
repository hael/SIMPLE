program simple_test_cartesian_fourier
use cartesian_fourier_kb_test, only: run_kb_derivative
use cartesian_fourier_neutral_extract_test, only: run_neutral_extract
use cartesian_fourier_test_helpers, only: run_current_executable_case
implicit none
character(len=256) :: argument, selected_case
integer :: case_count, failures, groups_passed, groups_run, groups_skipped, iarg, separator

selected_case = ''
case_count = 0
do iarg = 1, command_argument_count()
    call get_command_argument(iarg,argument)
    separator = index(argument,'=')
    if( separator > 0 )then
        if( trim(argument(:separator-1)) == 'case' )then
            case_count = case_count+1
            selected_case = trim(argument(separator+1:))
        endif
    endif
enddo
if( case_count > 1 ) error stop 'neutral Cartesian suite accepts only one case= argument'
if( case_count == 1 .and. len_trim(selected_case) == 0 ) &
    &error stop 'neutral Cartesian suite requires a nonempty case= value'
if( case_count == 0 )then
    failures = 0
    groups_passed = 0
    groups_run = 0
    groups_skipped = 0
    write(*,'(a)') 'Neutral Cartesian Fourier test suite: START'
    call run_current_executable_case('kb_derivative',groups_run,groups_passed,groups_skipped,failures)
    call run_current_executable_case('neutral_extract',groups_run,groups_passed,groups_skipped,failures)
    write(*,'(/,a)') 'Neutral Cartesian Fourier test suite: SUMMARY'
    write(*,'(a,i0)') '  Test groups scheduled: ',2
    write(*,'(a,i0)') '  Test groups run:       ',groups_run
    write(*,'(a,i0)') '  Test groups passed:    ',groups_passed
    write(*,'(a,i0)') '  Test groups skipped:   ',groups_skipped
    write(*,'(a,i0)') '  Test groups failed:    ',failures
    if( failures /= 0 ) error stop 1
    write(*,'(a)') 'Neutral Cartesian Fourier test suite: PASS'
else
    select case(trim(selected_case))
    case('kb_derivative')
        call run_kb_derivative()
    case('neutral_extract')
        call run_neutral_extract()
    case default
        error stop 'unknown neutral Cartesian Fourier case: '//trim(selected_case)
    end select
endif
end program simple_test_cartesian_fourier
