!@descr: for all utils tests
module simple_commanders_test_utils
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_ansi_colors
  contains
    procedure :: execute      => exec_test_ansi_colors
end type commander_test_ansi_colors

type, extends(commander_base) :: commander_test_binoris_test
  contains
    procedure :: execute      => exec_test_binoris_test
end type commander_test_binoris_test

type, extends(commander_base) :: commander_test_binoris_io_test
  contains
    procedure :: execute      => exec_test_binoris_io_test
end type commander_test_binoris_io_test

type, extends(commander_base) :: commander_test_cmdline
  contains
    procedure :: execute      => exec_test_cmdline
end type commander_test_cmdline

type, extends(commander_base) :: commander_test_install
  contains
    procedure :: execute      => exec_test_install
end type commander_test_install

type, extends(commander_base) :: commander_test_nice
  contains
    procedure :: execute      => exec_test_nice
end type commander_test_nice

type, extends(commander_base) :: commander_test_serialize
  contains
    procedure :: execute      => exec_test_serialize
end type commander_test_serialize

type, extends(commander_base) :: commander_test_stringmatch
  contains
    procedure :: execute      => exec_test_stringmatch
end type commander_test_stringmatch

type, extends(commander_base) :: commander_test_units
  contains
    procedure :: execute      => exec_test_units
end type commander_test_units

contains

subroutine exec_test_ansi_colors( self, cline )
    class(commander_test_ansi_colors),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_ANSI_COLORS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_ansi_colors

subroutine exec_test_binoris_test( self, cline )
    class(commander_test_binoris_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_BINORIS_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_binoris_test

subroutine exec_test_binoris_io_test( self, cline )
    class(commander_test_binoris_io_test),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_BINORIS_IO_TEST_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_binoris_io_test

subroutine exec_test_cmdline( self, cline )
    class(commander_test_cmdline),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_CMDLINE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_cmdline

subroutine exec_test_install( self, cline )
    class(commander_test_install),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_INSTALL_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_install

subroutine exec_test_nice( self, cline )
    class(commander_test_nice),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_NICE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_nice

subroutine exec_test_serialize( self, cline )
    class(commander_test_serialize),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SERIALIZE_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_serialize

subroutine exec_test_stringmatch( self, cline )
    class(commander_test_stringmatch),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_STRINGMATCH_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_stringmatch

subroutine exec_test_units( self, cline )
    class(commander_test_units),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_UNITS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_units

end module simple_commanders_test_utils
