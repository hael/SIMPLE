!@descr: for all parallel tests
module simple_commanders_test_parallel
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_coarrays
  contains
    procedure :: execute      => exec_test_coarrays
end type commander_test_coarrays

type, extends(commander_base) :: commander_test_openacc
  contains
    procedure :: execute      => exec_test_openacc
end type commander_test_openacc

type, extends(commander_base) :: commander_test_openmp
  contains
    procedure :: execute      => exec_test_openmp
end type commander_test_openmp

type, extends(commander_base) :: commander_test_simd
  contains
    procedure :: execute      => exec_test_simd
end type commander_test_simd

contains

subroutine exec_test_coarrays( self, cline )
    class(commander_test_coarrays),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_COARRAYS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_coarrays

subroutine exec_test_openacc( self, cline )
    class(commander_test_openacc),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_OPENACC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_openacc

subroutine exec_test_openmp( self, cline )
    class(commander_test_openmp),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_OPENMP_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_openmp

subroutine exec_test_simd( self, cline )
    class(commander_test_simd),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SIMD_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_simd

end module simple_commanders_test_parallel
