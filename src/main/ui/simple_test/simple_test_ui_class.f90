!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_class
use simple_ui_modules
implicit none

type(ui_program), target :: image

contains

    subroutine construct_test_class_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_image(tsttab)
    end subroutine construct_test_class_programs

    subroutine print_test_class_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('CLASS:', C_UNDERLINED)
        write(logfhandle,'(A)') image%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_class_programs

    subroutine new_image( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call image%new(&
        &'image',&                       ! name
        &'image',&                       ! descr_short
        &'is a test program for image',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call image%add_input(UI_IO, )
        ! parameter input/output
        !call image%add_input(UI_IMG, )
        ! alternative inputs
        !call image%add_input(UI_PARM, )
        ! search controls
        !call image%add_input(UI_SRCH, )
        ! filter controls
        !call image%add_input(UI_FILT, )
        ! mask controls
        !call image%add_input(UI_MASK, )
        ! computer controls
        !call image%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('image', image, tsttab)
    end subroutine new_image

end module simple_test_ui_class
