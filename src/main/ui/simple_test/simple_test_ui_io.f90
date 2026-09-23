!@descr: user interfaces for the input/output test programs run by hand on user data (simple_test_exec)
! The hermetic I/O tests are the stack I/O (unit_core), binoris and STAR sub-suites (unit_project).
module simple_test_ui_io
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('io', 'Input/Output', 50)
type(ui_program), target :: mrc2jpeg
type(ui_program), target :: mrc_validate

contains

    subroutine construct_test_io_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mrc2jpeg(tsttab)
        call new_mrc_validate(tsttab)
    end subroutine construct_test_io_programs

    subroutine print_test_io_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('INPUT/OUTPUT:', C_UNDERLINED)
        write(logfhandle,'(A)') mrc2jpeg%name%to_char()
        write(logfhandle,'(A)') mrc_validate%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_io_programs

    subroutine new_mrc2jpeg( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mrc2jpeg%new(&
        &'mrc2jpeg',&                          ! name
        &'mrc2jpeg ',&                         ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mrc2jpeg%add_input(UI_IO, )
        ! parameter input/output
        !call mrc2jpeg%add_input(UI_IMG, )
        ! <no additional inputs>
        !call mrc2jpeg%add_input(UI_PARM, )
        ! search controls
        !call mrc2jpeg%add_input(UI_SRCH, )
        ! filter controls
        !call mrc2jpeg%add_input(UI_FILT, )
        ! mask controls
        !call mrc2jpeg%add_input(UI_MASK, )
        ! computer controls
        !call mrc2jpeg%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mrc2jpeg', mrc2jpeg, tsttab, UI_CATEGORY)
    end subroutine new_mrc2jpeg

    subroutine new_mrc_validate( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mrc_validate%new(&
        &'mrc_validate',&                       ! name
        &'mrc_validate ',&                      ! summary
        &'is a test program for ',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call mrc_validate%add_input(UI_IMG, 'vol', 'file', 'Input volume', &
            &'MRC volume to validate', 'volume.mrc', .true., '')
        ! parameter input/output
        call mrc_validate%add_input(UI_PARM, 'smpd', 'real', 'Sampling distance', &
            &'Sampling distance in Angstrom per voxel', 'e.g. 1.3', .true., '')
        ! alternative inputs
        !call mrc_validate%add_input(UI_PARM, )
        ! <no additional inputs>
        ! search controls
        !call mrc_validate%add_input(UI_SRCH, )
        ! filter controls
        !call mrc_validate%add_input(UI_FILT, )
        ! mask controls
        !call mrc_validate%add_input(UI_MASK, )
        ! computer controls
        !call mrc_validate%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mrc_validate', mrc_validate, tsttab, UI_CATEGORY)
    end subroutine new_mrc_validate

end module simple_test_ui_io
