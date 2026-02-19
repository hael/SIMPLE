!@descr: for all io tests
module simple_commanders_test_io
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_imgfile
  contains
    procedure :: execute      => exec_test_imgfile
end type commander_test_imgfile

type, extends(commander_base) :: commander_test_inside_write
  contains
    procedure :: execute      => exec_test_inside_write
end type commander_test_inside_write

type, extends(commander_base) :: commander_test_io
  contains
    procedure :: execute      => exec_test_io
end type commander_test_io

type, extends(commander_base) :: commander_test_io_parallel
  contains
    procedure :: execute      => exec_test_io_parallel
end type commander_test_io_parallel

type, extends(commander_base) :: commander_test_mrc2jpeg
  contains
    procedure :: execute      => exec_test_mrc2jpeg
end type commander_test_mrc2jpeg

type, extends(commander_base) :: commander_test_mrc_validation
  contains
    procedure :: execute      => exec_test_mrc_validation
end type commander_test_mrc_validation

type, extends(commander_base) :: commander_test_stack_io
  contains
    procedure :: execute      => exec_test_stack_io
end type commander_test_stack_io

type, extends(commander_base) :: commander_test_star_export
  contains
    procedure :: execute      => exec_test_star_export
end type commander_test_star_export

type, extends(commander_base) :: commander_test_starfile_test
  contains
    procedure :: execute      => exec_test_starfile_test
end type commander_test_starfile_test

contains

subroutine exec_test_imgfile( self, cline )
    class(commander_test_imgfile),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_IMGFILE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_imgfile

subroutine exec_test_inside_write( self, cline )
    class(commander_test_inside_write),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_INSIDE_WRITE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_inside_write

subroutine exec_test_io( self, cline )
    class(commander_test_io),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_IO_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_io

subroutine exec_test_io_parallel( self, cline )
    class(commander_test_io_parallel),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_IO_PARALLEL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_io_parallel

subroutine exec_test_mrc2jpeg( self, cline )
    class(commander_test_mrc2jpeg),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MRC2JPEG_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mrc2jpeg

subroutine exec_test_mrc_validation( self, cline )
    class(commander_test_mrc_validation),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_MRC_VALIDATION_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mrc_validation

subroutine exec_test_stack_io( self, cline )
    class(commander_test_stack_io),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_STACK_IO_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_stack_io

subroutine exec_test_star_export( self, cline )
    class(commander_test_star_export),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_STAR_EXPORT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_star_export

subroutine exec_test_starfile_test( self, cline )
    class(commander_test_starfile_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_STARFILE_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_starfile_test

end module simple_commanders_test_io
