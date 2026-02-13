!@descr: "tests" UI api (concrete implementation)
module simple_ui_api_tests
use simple_ui_api_modules
implicit none

type(ui_program), target :: test_sim_workflow

contains

    subroutine register_simple_ui_tests(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('test_sim_workflow', test_sim_workflow, prgtab)
    end subroutine register_simple_ui_tests

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_test_sim_workflow
        ! PROGRAM SPECIFICATION
        call test_sim_workflow%new(&
        &'test_sim_workflow', &                                                                       ! name
        &'Test simulation workflow',&                                                                 ! descr_short
        &'is test',&                                                                                  ! descr long
        &'simple_test_exec',&                                                                         ! executable
        &.false.)                                                                                     ! requires sp_project
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
    end subroutine new_test_sim_workflow


end module simple_ui_api_tests
