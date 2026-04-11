!@descr: module defining the user interfaces for single test programs in the simple_test_exec suite
module simple_test_ui_single
use simple_ui_modules
implicit none

type(ui_program), target :: single_workflow 

contains

    subroutine construct_test_single_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_single_workflow(tsttab)
    end subroutine construct_test_single_programs

    subroutine print_test_single_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('SINGLE:', C_UNDERLINED)
        write(logfhandle,'(A)') single_workflow%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_single_programs

    subroutine new_single_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call single_workflow%new(&
        &'single_workflow',&                            ! name
        &'single workflow',&                           ! descr_short
        &'is a test program for single workflow',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call single_workflow%add_input(UI_IO, )
        ! parameter input/output
        !call single_workflow%add_input(UI_IMG, )
        ! alternative inputs
        !call single_workflow%add_input(UI_PARM, )
        ! search controls
        !call single_workflow%add_input(UI_SRCH, )
        ! filter controls
        !call single_workflow%add_input(UI_FILT, )
        ! mask controls
        !call single_workflow%add_input(UI_MASK, )
        ! computer controls
        !call single_workflow%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('single_workflow', single_workflow, tsttab)
    end subroutine new_single_workflow

end module simple_test_ui_single
