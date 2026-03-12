!@descr: module defining the user interfaces for per class test programs in the simple_test_exec suite
module simple_test_ui_class
use simple_ui_modules
implicit none

type(ui_program), target :: units, strategy2D

contains

    subroutine construct_test_class_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_units(tsttab)
        call new_strategy2D(tsttab)
    end subroutine construct_test_class_programs

    subroutine print_test_class_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('CLASS:', C_UNDERLINED)
        write(logfhandle,'(A)') strategy2D%name%to_char()
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

    subroutine new_strategy2D( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call strategy2D%new(&
        &'strategy2D',&                      ! name
        &'strategy2D ',&                     ! descr_short
        &'is a test program for strategy2D',&
        &'simple_test_exec',&                ! executable
        &.false.)                            ! requires sp_project
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
        call add_ui_program('strategy2D', strategy2D, tsttab)
    end subroutine new_strategy2D

end module simple_test_ui_class
