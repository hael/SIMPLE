!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(ui_program), target :: mini_stream
type(ui_program), target :: simulated_workflow

contains

    subroutine construct_test_highlevel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mini_stream(tsttab)
        call new_simulated_workflow(tsttab)
    end subroutine construct_test_highlevel_programs

    subroutine print_test_highlevel_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('HIGH-LEVEL:', C_UNDERLINED)
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') simulated_workflow%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_highlevel_programs

    subroutine new_mini_stream( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mini_stream%new(&
        &'mini_stream',&                       ! name
        &'mini_stream',&                       ! descr_short
        &'is a test program for mini_stream',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mini_stream%add_input(UI_IO, )
        ! parameter input/output
        !call mini_stream%add_input(UI_IMG, )
        ! alternative inputs
        !call mini_stream%add_input(UI_PARM, )
        ! search controls
        !call mini_stream%add_input(UI_SRCH, )
        ! filter controls
        !call mini_stream%add_input(UI_FILT, )
        ! mask controls
        !call mini_stream%add_input(UI_MASK, )
        ! computer controls
        !call mini_stream%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mini_stream', mini_stream, tsttab)
    end subroutine new_mini_stream

    subroutine new_simulated_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulated_workflow%new(&
        &'simulated_workflow',&                          ! name
        &'self-contained simulated workflow test',&      ! descr_short
        &'is a self-contained simulated workflow test',&
        &'simple_test_exec',&                            ! executable
        &.false.)                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call simulated_workflow%add_input(UI_IO, )
        ! parameter input/output
        !call simulated_workflow%add_input(UI_IMG, )
        ! alternative inputs
        !call simulated_workflow%add_input(UI_PARM, )
        ! search controls
        !call simulated_workflow%add_input(UI_SRCH, )
        ! filter controls
        !call simulated_workflow%add_input(UI_FILT, )
        ! mask controls
        !call simulated_workflow%add_input(UI_MASK, )
        ! computer controls
        call simulated_workflow%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('simulated_workflow', simulated_workflow, tsttab)
    end subroutine new_simulated_workflow

end module simple_test_ui_highlevel
