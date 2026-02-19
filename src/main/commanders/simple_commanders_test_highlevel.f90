!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_mini_stream
  contains
    procedure :: execute      => exec_test_mini_stream
end type commander_test_mini_stream

type, extends(commander_base) :: commander_test_simulated_workflow
  contains
    procedure :: execute      => exec_test_simulated_workflow
end type commander_test_simulated_workflow

contains

subroutine exec_test_mini_stream( self, cline )
    class(commander_test_mini_stream),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MINI_STREAM_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mini_stream

subroutine exec_test_simulated_workflow( self, cline )
    class(commander_test_simulated_workflow),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SIMULATED_WORKFLOW_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_simulated_workflow

end module simple_commanders_test_highlevel
