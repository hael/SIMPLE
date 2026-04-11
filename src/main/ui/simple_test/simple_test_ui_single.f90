!@descr: module defining the user interfaces for single test programs in the simple_test_exec suite
module simple_test_ui_single
use simple_ui_modules
implicit none

type(ui_program), target :: atoms_stats
type(ui_program), target :: detect_atoms
type(ui_program), target :: simulate_nanoprticle
type(ui_program), target :: single_workflow 

contains

    subroutine construct_test_single_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_atoms_stats(tsttab)
        call new_detect_atoms(tsttab)
        call new_simulate_nanoprticle(tsttab)
        call new_single_workflow(tsttab)
    end subroutine construct_test_single_programs

    subroutine print_test_single_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('SINGLE:', C_UNDERLINED)
        write(logfhandle,'(A)') atoms_stats%name%to_char()
        write(logfhandle,'(A)') detect_atoms%name%to_char()
        write(logfhandle,'(A)') simulate_nanoprticle%name%to_char()
        write(logfhandle,'(A)') single_workflow%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_single_programs
 
    subroutine new_atoms_stats( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call atoms_stats%new(&
        &'atoms_stats',&                            ! name
        &'test program for atom stats',&
        &'is a test program for atom stats',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call atoms_stats%add_input(UI_IO, )
        ! parameter input/output
        !call atoms_stats%add_input(UI_IMG, )
        ! alternative inputs
        !call atoms_stats%add_input(UI_PARM, )
        ! search controls
        !call atoms_stats%add_input(UI_SRCH, )
        ! filter controls
        !call atoms_stats%add_input(UI_FILT, )
        ! mask controls
        !call atoms_stats%add_input(UI_MASK, )
        ! computer controls
        !call atoms_stats%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('atoms_stats', atoms_stats, tsttab)
    end subroutine new_atoms_stats

    subroutine new_detect_atoms( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call detect_atoms%new(&
        &'detect_atoms',&                            ! name
        &'test program for atom detection',&
        &'is a test program for atom detection',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call detect_atoms%add_input(UI_IO, )
        ! parameter input/output
        !call detect_atoms%add_input(UI_IMG, )
        ! alternative inputs
        !call detect_atoms%add_input(UI_PARM, )
        ! search controls
        !call detect_atoms%add_input(UI_SRCH, )
        ! filter controls
        !call detect_atoms%add_input(UI_FILT, )
        ! mask controls
        !call detect_atoms%add_input(UI_MASK, )
        ! computer controls
        !call detect_atoms%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('detect_atoms', detect_atoms, tsttab)
    end subroutine new_detect_atoms

    subroutine new_simulate_nanoprticle( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulate_nanoprticle%new(&
        &'simulate_nanoprticle',&                            ! name
        &'test program for simulating nanoparticle',&
        &'is a test program for simulating nanoparticle',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call simulate_nanoprticle%add_input(UI_IO, )
        ! parameter input/output
        !call simulate_nanoprticle%add_input(UI_IMG, )
        ! alternative inputs
        !call simulate_nanoprticle%add_input(UI_PARM, )
        ! search controls
        !call simulate_nanoprticle%add_input(UI_SRCH, )
        ! filter controls
        !call simulate_nanoprticle%add_input(UI_FILT, )
        ! mask controls
        !call simulate_nanoprticle%add_input(UI_MASK, )
        ! computer controls
        !call simulate_nanoprticle%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('simulate_nanoprticle', simulate_nanoprticle, tsttab)
    end subroutine new_simulate_nanoprticle

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
