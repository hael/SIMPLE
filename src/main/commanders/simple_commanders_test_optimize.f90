!@descr: for all optimize tests
module simple_commanders_test_optimize
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_lbfgsb
  contains
    procedure :: execute      => exec_test_lbfgsb
end type commander_test_lbfgsb

type, extends(commander_base) :: commander_test_lbfgsb_cosine
  contains
    procedure :: execute      => exec_test_lbfgsb_cosine
end type commander_test_lbfgsb_cosine

type, extends(commander_base) :: commander_test_lplims
  contains
    procedure :: execute      => exec_test_lplims
end type commander_test_lplims

type, extends(commander_base) :: commander_test_lpstages_test
  contains
    procedure :: execute      => exec_test_lpstages_test
end type commander_test_lpstages_test

type, extends(commander_base) :: commander_test_opt_lp
  contains
    procedure :: execute      => exec_test_opt_lp
end type commander_test_opt_lp

type, extends(commander_base) :: commander_test_tree_srch
  contains
    procedure :: execute      => exec_test_tree_srch
end type commander_test_tree_srch

contains

subroutine exec_test_lbfgsb( self, cline )
    class(commander_test_lbfgsb),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_LBFGSB_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lbfgsb

subroutine exec_test_lbfgsb_cosine( self, cline )
    class(commander_test_lbfgsb_cosine),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_LBFGSB_COSINE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lbfgsb_cosine

subroutine exec_test_lplims( self, cline )
    class(commander_test_lplims),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_LPLIMS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lplims

subroutine exec_test_lpstages_test( self, cline )
    class(commander_test_lpstages_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_LPSTAGES_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_lpstages_test

subroutine exec_test_opt_lp( self, cline )
    class(commander_test_opt_lp),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_OPT_LP_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_opt_lp

subroutine exec_test_tree_srch( self, cline )
    class(commander_test_tree_srch),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_TREE_SRCH_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_tree_srch

end module simple_commanders_test_optimize
