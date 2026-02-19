!@descr: for all numerics tests
module simple_commanders_test_numerics
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_eigh_test
  contains
    procedure :: execute      => exec_test_eigh_test
end type commander_test_eigh_test

type, extends(commander_base) :: commander_test_kbinterpol_fast
  contains
    procedure :: execute      => exec_test_kbinterpol_fast
end type commander_test_kbinterpol_fast

type, extends(commander_base) :: commander_test_maxnloc_test
  contains
    procedure :: execute      => exec_test_maxnloc_test
end type commander_test_maxnloc_test

type, extends(commander_base) :: commander_test_neigh
  contains
    procedure :: execute      => exec_test_neigh
end type commander_test_neigh

contains

subroutine exec_test_eigh_test( self, cline )
    class(commander_test_eigh_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_EIGH_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_eigh_test

subroutine exec_test_kbinterpol_fast( self, cline )
    class(commander_test_kbinterpol_fast),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_KBINTERPOL_FAST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_kbinterpol_fast

subroutine exec_test_maxnloc_test( self, cline )
    class(commander_test_maxnloc_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MAXNLOC_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_maxnloc_test

subroutine exec_test_neigh( self, cline )
    class(commander_test_neigh),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_NEIGH_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_neigh

end module simple_commanders_test_numerics
