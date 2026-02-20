!@descr: module defining the user interfaces for per class test programs in the simple_test_exec suite
module simple_test_ui_class
use simple_ui_modules
implicit none

type(ui_program), target :: units

contains

    subroutine construct_test_class_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_units(tsttab)
    end subroutine construct_test_class_programs

    subroutine print_test_class_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('CLASS:', C_UNDERLINED)
        write(logfhandle,'(A)') units%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_class_programs

    subroutine new_units( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call units%new(&
        &'units',&                      ! name
        &'units ',&                     ! descr_short
        &'is a test program for units',&
        &'simple_test_exec',&           ! executable
        &.false.)                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call units%add_input(UI_IO, )
        ! parameter input/output
        !call units%add_input(UI_IMG, )
        ! alternative inputs
        !call units%add_input(UI_PARM, )
        ! search controls
        !call units%add_input(UI_SRCH, )
        ! filter controls
        !call units%add_input(UI_FILT, )
        ! mask controls
        !call units%add_input(UI_MASK, )
        ! computer controls
        !call units%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('units', units, tsttab)
    end subroutine new_units

end module simple_test_ui_class
