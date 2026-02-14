!@descr: "tests" UI api (concrete implementation)
module simple_ui_api_tests
use simple_ui_api_modules
implicit none

type(ui_program), target :: test_sim_workflow

contains

    subroutine construct_tests_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_test_sim_workflow(prgtab)
    end subroutine construct_tests_programs

    subroutine print_tests_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('TESTS:', C_UNDERLINED)
        write(logfhandle,'(A)') test_sim_workflow%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_tests_programs

    subroutine new_test_sim_workflow( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call test_sim_workflow%new(&
        &'test_sim_workflow', &         ! name
        &'Test simulation workflow',&   ! descr_short
        &'is test',&                    ! descr long
        &'simple_test_exec',&           ! executable
        &.false.)                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('test_sim_workflow', test_sim_workflow, prgtab)
    end subroutine new_test_sim_workflow

end module simple_ui_api_tests
