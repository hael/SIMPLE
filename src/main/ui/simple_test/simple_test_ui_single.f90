!@descr: module defining the user interfaces for single test programs in the simple_test_exec suite
module simple_test_ui_single
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('single', 'SINGLE', 110)
type(ui_program), target :: atoms_stats
type(ui_program), target :: detect_calpha_molecules
type(ui_program), target :: single_workflow 

contains

    subroutine construct_test_single_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_atoms_stats(tsttab)
        call new_detect_calpha_molecules(tsttab)
        call new_single_workflow(tsttab)
    end subroutine construct_test_single_programs

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
        call atoms_stats%add_input(UI_PARM, smpd,    required_override=.true.)
        ! <no additional inputs>
        !call atoms_stats%add_input(UI_PARM, )
        ! search controls
        !call atoms_stats%add_input(UI_SRCH, )
        ! filter controls
        call atoms_stats%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        !call atoms_stats%add_input(UI_FILT, )
        ! mask controls
        !call atoms_stats%add_input(UI_MASK, )
        ! computer controls
        !call atoms_stats%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('atoms_stats', atoms_stats, tsttab, UI_CATEGORY)
    end subroutine new_atoms_stats

    subroutine new_detect_calpha_molecules(tsttab)
        class(ui_hash), intent(inout) :: tsttab
        call detect_calpha_molecules%new(&
        &'detect_calpha_molecules',&
        &'Evaluate C-alpha detection on the built-in 6VXX and 1JYX molecular models',&
        &'Simulates maps from the hard-coded structures and reports uniquely matched C-alpha counts, recall, and precision.',&
        &'simple_test_exec',&
        &.false.)
        call detect_calpha_molecules%add_input(UI_PARM, smpd, required_override=.false.)
        call detect_calpha_molecules%add_input(UI_SRCH, 'angstep', 'num', 'Angular spacing', &
        &'Approximate spacing of the coarse SO(3) orientation grid in degrees', &
        &'degrees{45}', .false., 45.0)
        call detect_calpha_molecules%add_input(UI_SRCH, 'thres', 'num', 'Minimum score', &
        &'Minimum weighted normalized correlation score accepted as a candidate', &
        &'correlation score{0.25}', .false., 0.25)
        call detect_calpha_molecules%add_input(UI_COMP, nthr)
        call add_ui_program('detect_calpha_molecules', detect_calpha_molecules, tsttab, UI_CATEGORY)
    end subroutine new_detect_calpha_molecules

    subroutine new_single_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call single_workflow%new(&
        &'single_workflow',&                            ! name
        &'single workflow',&                           ! summary
        &'is a test program for single workflow',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call single_workflow%add_input(UI_IO, )
        ! parameter input/output
        call single_workflow%add_input(UI_PARM, smpd,    required_override=.true.)
        !call single_workflow%add_input(UI_IMG, )
        ! <no additional inputs>
        !call single_workflow%add_input(UI_PARM, )
        ! search controls
        !call single_workflow%add_input(UI_SRCH, )
        ! filter controls
        call single_workflow%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        !call single_workflow%add_input(UI_FILT, )
        ! mask controls
        !call single_workflow%add_input(UI_MASK, )
        ! computer controls
        !call single_workflow%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('single_workflow', single_workflow, tsttab, UI_CATEGORY)
    end subroutine new_single_workflow

end module simple_test_ui_single
