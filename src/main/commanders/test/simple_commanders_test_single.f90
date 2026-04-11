!@descr: tests for single 
module simple_commanders_test_single
use simple_commanders_api
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_single_workflow
  contains
    procedure :: execute      => exec_test_single_workflow
end type commander_test_single_workflow

contains

subroutine exec_test_single_workflow( self, cline )
    class(commander_test_single_workflow), intent(inout) :: self
    class(cmdline),                        intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SINGLE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_single_workflow

end module simple_commanders_test_single
