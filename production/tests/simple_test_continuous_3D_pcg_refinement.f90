program simple_test_continuous_3D_pcg_refinement
use continuous_3D_pcg_refinement_test_helpers, only: CASE_SKIP_EXIT_STATUS
use continuous_3D_pcg_refinement_scaffold_test, only: run_scaffold_contract
use continuous_3D_pcg_refinement_volume_test, only: run_volume_fixture
use continuous_3D_pcg_refinement_noise_test, only: run_volume_noise
use continuous_3D_pcg_refinement_halfset_test, only: run_halfset_fsc
use continuous_3D_pcg_refinement_shift_test, only: run_shift_gradient
use continuous_3D_pcg_refinement_shift_polish_test, only: run_shift_polish
use continuous_3D_pcg_refinement_kb_test, only: run_kb_derivative
use continuous_3D_pcg_refinement_rotation_test, only: run_rotation_gradient
use continuous_3D_pcg_refinement_recovery_test, only: run_pose_recovery
use continuous_3D_pcg_refinement_pose_contract_test, only: run_pose_contract
use continuous_3D_pcg_refinement_fixed_reference_test, only: run_fixed_reference_diagnostic
use continuous_3D_pcg_refinement_forward_path_test, only: run_forward_path_diagnostic
use continuous_3D_pcg_refinement_matched_window_test, only: run_matched_window_diagnostic
use continuous_3D_pcg_refinement_reference_bias_test, only: run_reference_bias_diagnostic
use continuous_3D_pcg_refinement_operator_contract_test, only: run_operator_contract_diagnostic
implicit none

character(len=256) :: selected_case
integer :: case_argument_count

call find_selected_case(selected_case, case_argument_count)
if( case_argument_count > 1 ) &
    &error stop 'continuous 3D PCG suite accepts only one case= argument'
if( case_argument_count == 1 .and. len_trim(selected_case) == 0 ) &
    &error stop 'continuous 3D PCG suite requires a nonempty case= value'

if( case_argument_count == 1 )then
    call run_case(selected_case)
else
    call run_suite()
endif

contains

!> Parse one optional case key and reject ambiguous child selection.
subroutine find_selected_case(case_name, occurrences)
    character(len=*), intent(out) :: case_name
    integer,          intent(out) :: occurrences
    character(len=4096) :: argument
    integer :: iarg, argument_status, separator

    case_name = ''
    occurrences = 0
    do iarg = 1, command_argument_count()
        call get_command_argument(iarg, argument, status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read continuous 3D PCG suite argument'
        separator = index(argument, '=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= 'case' ) cycle
        occurrences = occurrences + 1
        case_name = trim(argument(separator+1:))
    enddo
end subroutine find_selected_case

!> Report whether one command argument is a case selector.
logical function is_case_argument(argument) result(is_case)
    character(len=*), intent(in) :: argument
    integer :: separator

    separator = index(argument, '=')
    is_case = separator > 0
    if( is_case ) is_case = trim(argument(:separator-1)) == 'case'
end function is_case_argument

!> Schedule each scientific case in a separate child process and summarize outcomes.
subroutine run_suite()
    character(len=*), parameter :: labels(*) = [character(len=24) :: &
        &'scaffold', 'volume_fixture', 'volume_noise', 'halfset_fsc', &
        &'shift_gradient', 'shift_polish', 'kb_derivative', 'rotation_gradient', 'pose_recovery', &
        &'pose_contract']
    integer :: failures, groups_passed, groups_run, groups_skipped, i

    failures = 0
    groups_passed = 0
    groups_run = 0
    groups_skipped = 0
    write(*,'(a)') 'Continuous 3D PCG refinement test suite: START'
    do i = 1, size(labels)
        call run_case_in_child(trim(labels(i)), groups_run, groups_passed, &
            &groups_skipped, failures)
    enddo

    write(*,'(/,a)') 'Continuous 3D PCG refinement test suite: SUMMARY'
    write(*,'(a,i0)') '  Test groups scheduled: ', size(labels)
    write(*,'(a,i0)') '  Test groups run:       ', groups_run
    write(*,'(a,i0)') '  Test groups passed:    ', groups_passed
    write(*,'(a,i0)') '  Test groups skipped:   ', groups_skipped
    write(*,'(a,i0)') '  Test groups failed:    ', failures
    if( failures /= 0 )then
        write(*,'(a)') 'Continuous 3D PCG refinement test suite: FAIL'
        error stop 1
    endif
    write(*,'(a)') 'Continuous 3D PCG refinement test suite: PASS'
end subroutine run_suite

!> Run one named case in an isolated child so one failure does not hide later cases.
subroutine run_case_in_child(label, groups_run, groups_passed, groups_skipped, failures)
    character(len=*), intent(in) :: label
    integer, intent(inout) :: groups_run, groups_passed, groups_skipped, failures
    character(len=:), allocatable :: child_command, executable
    character(len=4096) :: parent_argument
    character(len=1024) :: command_message
    integer :: argument_status, command_status, executable_len, exit_status, iarg

    call get_command_argument(0, length=executable_len, status=argument_status)
    if( argument_status /= 0 .or. executable_len < 1 ) &
        &error stop 'could not determine continuous 3D PCG suite executable path'
    allocate(character(len=executable_len) :: executable)
    call get_command_argument(0, executable, status=argument_status)
    if( argument_status /= 0 ) &
        &error stop 'could not read continuous 3D PCG suite executable path'

    child_command = executable_command(executable)//' '//quoted_shell_token('case='//trim(label))
    do iarg = 1, command_argument_count()
        call get_command_argument(iarg, parent_argument, status=argument_status)
        if( argument_status /= 0 ) error stop 'could not forward continuous 3D PCG suite argument'
        if( is_case_argument(parent_argument) ) cycle
        child_command = child_command//' '//quoted_shell_token(trim(parent_argument))
    enddo

    groups_run = groups_run + 1
    command_status = -1
    exit_status = -1
    command_message = ''
    write(*,'(/,a)') '>>> TEST ['//trim(label)//'] START'
    call execute_command_line(child_command, wait=.true., exitstat=exit_status, &
        &cmdstat=command_status, cmdmsg=command_message)
    if( command_status /= 0 )then
        failures = failures + 1
        write(*,'(a,i0,a)') '>>> TEST ['//trim(label)//'] FAIL (cmdstat=', command_status, ')'
        if( len_trim(command_message) > 0 ) write(*,'(a)') '    '//trim(command_message)
    else if( exit_status == 0 )then
        groups_passed = groups_passed + 1
        write(*,'(a)') '>>> TEST ['//trim(label)//'] PASS'
    else if( exit_status == CASE_SKIP_EXIT_STATUS )then
        groups_skipped = groups_skipped + 1
        write(*,'(a)') '>>> TEST ['//trim(label)//'] SKIP'
    else
        failures = failures + 1
        write(*,'(a,i0,a)') '>>> TEST ['//trim(label)//'] FAIL (exitstat=', exit_status, ')'
    endif
end subroutine run_case_in_child

!> Quote the current executable path for safe child-process invocation.
function executable_command(executable) result(command)
    character(len=*), intent(in) :: executable
    character(len=:), allocatable :: command

#if defined(_WIN32)
    command = 'call '//quoted_shell_token(executable)
#else
    command = quoted_shell_token(executable)
#endif
end function executable_command

!> Escape one shell token without changing its command-line value.
function quoted_shell_token(value) result(quoted)
    character(len=*), intent(in) :: value
    character(len=:), allocatable :: quoted

#if defined(_WIN32)
    if( index(value, '"') > 0 ) &
        &error stop 'continuous 3D PCG suite cannot forward an argument containing a double quote'
    quoted = '"'//trim(value)//'"'
#else
    if( index(value, achar(39)) > 0 ) &
        &error stop 'continuous 3D PCG suite cannot forward an argument containing a single quote'
    quoted = achar(39)//trim(value)//achar(39)
#endif
end function quoted_shell_token

!> Dispatch one case key to its owning child test module.
subroutine run_case(label)
    character(len=*), intent(in) :: label

    write(*,'(a)') '>>> Continuous 3D PCG refinement child case: '//trim(label)
    select case(trim(label))
    case('scaffold')
        call run_scaffold_contract()
    case('volume_fixture')
        call run_volume_fixture()
    case('volume_noise')
        call run_volume_noise()
    case('halfset_fsc')
        call run_halfset_fsc()
    case('shift_gradient')
        call run_shift_gradient()
    case('shift_polish')
        call run_shift_polish()
    case('kb_derivative')
        call run_kb_derivative()
    case('rotation_gradient')
        call run_rotation_gradient()
    case('pose_recovery')
        call run_pose_recovery()
    case('pose_contract')
        call run_pose_contract()
    case('fixed_reference')
        call run_fixed_reference_diagnostic()
    case('forward_path')
        call run_forward_path_diagnostic()
    case('matched_window')
        call run_matched_window_diagnostic()
    case('reference_bias')
        call run_reference_bias_diagnostic()
    case('operator_contract')
        call run_operator_contract_diagnostic()
    case default
        error stop 'unknown continuous 3D PCG refinement case: '//trim(label)
    end select
end subroutine run_case

end program simple_test_continuous_3D_pcg_refinement
