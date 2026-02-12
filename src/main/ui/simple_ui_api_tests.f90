!@descr: "tests" UI api (concrete implementation)
module simple_ui_api_tests
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: test_sim_workflow

contains

    subroutine register_simple_ui_tests(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('test_sim_workflow', test_sim_workflow, prgtab)
    end subroutine register_simple_ui_tests

end module simple_ui_api_tests
