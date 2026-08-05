!@descr: focused dispatch and command-line contract test for the development abinitio2d_sgd commander
program simple_test_abinitio2d_sgd_dispatch
use simple_test_utils, only: assert_char, assert_true, report_summary, tests_failed
use simple_cmdline, only: cmdline
use simple_string, only: string
use simple_commanders_abinitio2d_sgd, only: prepare_abinitio2d_sgd
use simple_exec_cluster2D, only: exec_cluster2D_commander
implicit none

character(len=64)   :: test_case
character(len=1024) :: executable
integer :: argument_length, argument_status

call get_command_argument(1, test_case, argument_length, argument_status)
if( argument_status == 0 .and. argument_length > 5 )then
    if( test_case(:5) == 'case=' )then
        call run_invalid_child(test_case(6:argument_length))
        stop
    endif
endif

call get_command_argument(0, executable, argument_length, argument_status)
if( argument_status /= 0 .or. argument_length < 1 ) error stop 'could not determine dispatch-test executable path'

call test_valid_preparation()
call expect_invalid_child('objfun_cc', 'abinitio2d_sgd requires objfun=euclid')
call expect_invalid_child('inpl_cont_yes', 'abinitio2d_sgd does not support inpl_cont=yes')
call expect_invalid_child('stage4_off', 'abinitio2d_sgd requires sgd_stage4_mode=on or alternate')

call report_summary()
if( tests_failed /= 0 ) error stop 1

contains

    subroutine test_valid_preparation()
        type(cmdline) :: cline
        call cline%set('prg', 'abinitio2d_sgd')
        call prepare_abinitio2d_sgd(cline)
        call assert_cline_char(cline, 'objfun', 'euclid', 'wrapper defaults to Euclidean objective')
        call assert_cline_char(cline, 'inpl_cont', 'no', 'wrapper disables continuous in-plane path')
        call assert_cline_char(cline, 'sgd_stage4_mode', 'on', 'wrapper defaults to streamed stage policy')
        call assert_cline_char(cline, 'sgd', 'yes', 'wrapper enables internal SGD compatibility flag')
        call assert_cline_char(cline, 'sgd_path', 'stream', 'wrapper selects streaming assignment path')
        call assert_cline_char(cline, 'prg', 'abinitio2D', 'wrapper delegates with standard workflow label')
        call cline%kill()
    end subroutine test_valid_preparation

    subroutine expect_invalid_child(label, expected_text)
        character(len=*), intent(in) :: label, expected_text
        character(len=2048) :: command, logfile
        integer :: command_status, exit_status

        logfile = 'simple_test_abinitio2d_sgd_dispatch_'//trim(label)//'.log'
        command = '"'//trim(executable(:argument_length))//'" case='//trim(label)// &
            &' > "'//trim(logfile)//'" 2>&1'
        call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=command_status)
        call assert_true(command_status == 0, trim(label)//' child command starts')
        call assert_true(exit_status /= 0, trim(label)//' is rejected before workflow execution')
        call assert_true(file_contains(logfile, expected_text), trim(label)//' reports the expected wrapper validation')
    end subroutine expect_invalid_child

    subroutine run_invalid_child(label)
        character(len=*), intent(in) :: label
        type(cmdline) :: cline
        logical :: l_silent, l_did_execute

        call cline%set('prg', 'abinitio2d_sgd')
        select case(trim(label))
        case('objfun_cc')
            call cline%set('objfun', 'cc')
        case('inpl_cont_yes')
            call cline%set('inpl_cont', 'yes')
        case('stage4_off')
            call cline%set('sgd_stage4_mode', 'off')
        case default
            error stop 'unknown abinitio2d_sgd dispatch child case'
        end select
        l_did_execute = .false.
        call exec_cluster2D_commander('abinitio2d_sgd', cline, l_silent, l_did_execute)
        if( .not. l_did_execute ) error stop 'abinitio2d_sgd route was not registered'
        error stop 'invalid abinitio2d_sgd command reached workflow execution'
    end subroutine run_invalid_child

    subroutine assert_cline_char(cline, key, expected, message)
        class(cmdline), intent(in) :: cline
        character(len=*), intent(in) :: key, expected, message
        type(string) :: value
        value = cline%get_carg(key)
        call assert_char(expected, value%to_char(), message)
        call value%kill()
    end subroutine assert_cline_char

    logical function file_contains(filename, expected_text) result(found)
        character(len=*), intent(in) :: filename, expected_text
        character(len=1024) :: line
        integer :: file_unit, ios

        found = .false.
        open(newunit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(file_unit, '(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, trim(expected_text)) > 0 )then
                found = .true.
                exit
            endif
        enddo
        close(file_unit)
    end function file_contains

end program simple_test_abinitio2d_sgd_dispatch
