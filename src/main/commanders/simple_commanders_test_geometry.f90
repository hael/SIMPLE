!@descr: for all geometry tests
module simple_commanders_test_geometry
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_angres
  contains
    procedure :: execute      => exec_test_angres
end type commander_test_angres

type, extends(commander_base) :: commander_test_ori_test
  contains
    procedure :: execute      => exec_test_ori_test
end type commander_test_ori_test

type, extends(commander_base) :: commander_test_oris_test
  contains
    procedure :: execute      => exec_test_oris_test
end type commander_test_oris_test

type, extends(commander_base) :: commander_test_sym_test
  contains
    procedure :: execute      => exec_test_sym_test
end type commander_test_sym_test

type, extends(commander_base) :: commander_test_uniform_euler
  contains
    procedure :: execute      => exec_test_uniform_euler
end type commander_test_uniform_euler

type, extends(commander_base) :: commander_test_uniform_rot
  contains
    procedure :: execute      => exec_test_uniform_rot
end type commander_test_uniform_rot

contains

subroutine exec_test_angres( self, cline )
    class(commander_test_angres), intent(inout) :: self
    class(cmdline),               intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ANGRES_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_angres

subroutine exec_test_ori_test( self, cline )
    class(commander_test_ori_test), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ORI_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ori_test

subroutine exec_test_oris_test( self, cline )
    class(commander_test_oris_test), intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ORIS_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_oris_test

subroutine exec_test_sym_test( self, cline )
    class(commander_test_sym_test), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SYM_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_sym_test

subroutine exec_test_uniform_euler( self, cline )
    class(commander_test_uniform_euler), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_UNIFORM_EULER_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_uniform_euler

subroutine exec_test_uniform_rot( self, cline )
    class(commander_test_uniform_rot),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_UNIFORM_ROT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_uniform_rot

end module simple_commanders_test_geometry
