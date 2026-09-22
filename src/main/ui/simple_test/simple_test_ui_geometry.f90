!@descr: module defining the user interfaces for geometry test programs in the simple_test_exec suite
module simple_test_ui_geometry
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('geometry', 'Geometry', 30)
type(ui_program), target :: angres

contains

    subroutine construct_test_geometry_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_angres(tsttab)
    end subroutine construct_test_geometry_programs

    subroutine new_angres( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call angres%new(&
        &'angres',&                            ! name
        &'angular resolution of the projection-direction spiral',& ! summary
        &'checks find_angres on spirals of 500-4000 directions against recorded values',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call angres%add_input(UI_IO, )
        ! parameter input/output
        !call angres%add_input(UI_IMG, )
        ! <no additional inputs>
        !call angres%add_input(UI_PARM, )
        ! search controls
        !call angres%add_input(UI_SRCH, )
        ! filter controls
        !call angres%add_input(UI_FILT, )
        ! mask controls
        !call angres%add_input(UI_MASK, )
        ! computer controls
        !call angres%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('angres', angres, tsttab, UI_CATEGORY)
    end subroutine new_angres

end module simple_test_ui_geometry
