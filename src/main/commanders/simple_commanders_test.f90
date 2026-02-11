!@descr: for simulating noise, particles, movies, atoms, etc.
module simple_commanders_test
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_sim_workflow
  contains
    procedure :: execute      => exec_test_sim_workflow
end type commander_test_sim_workflow

contains

    subroutine exec_test_sim_workflow( self, cline )
        class(commander_test_sim_workflow), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        call simple_end('**** SIMPLE_TEST_SIM_WORKFLOW NORMAL STOP ****')
    end subroutine exec_test_sim_workflow

end module simple_commanders_test
