!@descr: tests for SIMPLE_stream workflows
module simple_commanders_test_stream
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_abinitio2D_stream
  contains
    procedure :: execute => exec_test_abinitio2D_stream
end type commander_test_abinitio2D_stream

type, extends(commander_base) :: commander_test_assign_optics
  contains
    procedure :: execute => exec_test_assign_optics
end type commander_test_assign_optics

type, extends(commander_base) :: commander_test_gen_pickrefs
  contains
    procedure :: execute => exec_test_gen_pickrefs
end type commander_test_gen_pickrefs

type, extends(commander_base) :: commander_test_master
  contains
    procedure :: execute => exec_test_master
end type commander_test_master

type, extends(commander_base) :: commander_test_pick_extract
  contains
    procedure :: execute => exec_test_pick_extract
end type commander_test_pick_extract

type, extends(commander_base) :: commander_test_preproc
  contains
    procedure :: execute => exec_test_preproc
end type commander_test_preproc

type, extends(commander_base) :: commander_test_sieve_cavgs
  contains
    procedure :: execute => exec_test_sieve_cavgs
end type commander_test_sieve_cavgs

contains

subroutine exec_test_abinitio2D_stream( self, cline )
    class(commander_test_abinitio2D_stream), intent(inout) :: self
    class(cmdline),                         intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_ABINITIO2D_STREAM: DUMMY'
end subroutine exec_test_abinitio2D_stream

subroutine exec_test_assign_optics( self, cline )
    class(commander_test_assign_optics), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_ASSIGN_OPTICS: DUMMY'
end subroutine exec_test_assign_optics

subroutine exec_test_gen_pickrefs( self, cline )
    class(commander_test_gen_pickrefs), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_GEN_PICKREFS: DUMMY'
end subroutine exec_test_gen_pickrefs

subroutine exec_test_master( self, cline )
    class(commander_test_master), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_MASTER: DUMMY'
end subroutine exec_test_master

subroutine exec_test_pick_extract( self, cline )
    class(commander_test_pick_extract), intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_PICK_EXTRACT: DUMMY'
end subroutine exec_test_pick_extract

subroutine exec_test_preproc( self, cline )
    class(commander_test_preproc), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_PREPROC: DUMMY'
end subroutine exec_test_preproc

subroutine exec_test_sieve_cavgs( self, cline )
    class(commander_test_sieve_cavgs), intent(inout) :: self
    class(cmdline),                    intent(inout) :: cline
    write(logfhandle,'(a)') '>>> TEST_SIEVE_CAVGS: DUMMY'
end subroutine exec_test_sieve_cavgs

end module simple_commanders_test_stream
