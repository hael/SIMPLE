!@ descr: module defining the user interfaces for highlevel programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(ui_program), target :: simple_test_mini_stream
type(ui_program), target :: simple_test_simulated_workflow

contains

    subroutine construct_highlevel_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_simple_test_mini_stream(prgtab)
        call new_simple_test_simulated_workflow(prgtab)
    end subroutine construct_highlevel_programs

    subroutine print_highlevel_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('HIGH-LEVEL:', C_UNDERLINED)
        write(logfhandle,'(A)') simple_test_mini_stream%name%to_char()
        write(logfhandle,'(A)') simple_test_simulated_workflow%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_highlevel_programs

    subroutine new_simple_test_mini_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_mini_stream', simple_test_mini_stream, prgtab)
    end subroutine new_simple_test_mini_stream

    subroutine new_simple_test_simulated_workflow( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! add to ui_hash
        call add_ui_program('simple_test_simulated_workflow', simple_test_simulated_workflow, prgtab)
    end subroutine new_simple_test_simulated_workflow

end module simple_test_ui_highlevel
