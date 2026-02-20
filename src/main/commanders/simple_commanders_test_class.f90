!@descr: for all class tests
module simple_commanders_test_class
use simple_commanders_api
use simple_stream_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_image
  contains
    procedure :: execute      => exec_test_image
end type commander_test_image

contains

subroutine exec_test_image( self, cline )
    class(commander_test_image),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_IMAGE NORMAL STOP ****')
end subroutine exec_test_image

end module simple_commanders_test_class
