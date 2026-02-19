!@descr: for all network tests
module simple_commanders_test_network
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_socket_client
  contains
    procedure :: execute      => exec_test_socket_client
end type commander_test_socket_client

type, extends(commander_base) :: commander_test_socket_comm_distr
  contains
    procedure :: execute      => exec_test_socket_comm_distr
end type commander_test_socket_comm_distr

type, extends(commander_base) :: commander_test_socket_io
  contains
    procedure :: execute      => exec_test_socket_io
end type commander_test_socket_io

type, extends(commander_base) :: commander_test_socket_server
  contains
    procedure :: execute      => exec_test_socket_server
end type commander_test_socket_server

contains

subroutine exec_test_socket_client( self, cline )
    class(commander_test_socket_client),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SOCKET_CLIENT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_client

subroutine exec_test_socket_comm_distr( self, cline )
    class(commander_test_socket_comm_distr),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SOCKET_COMM_DISTR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_comm_distr

subroutine exec_test_socket_io( self, cline )
    class(commander_test_socket_io),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SOCKET_IO_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_io

subroutine exec_test_socket_server( self, cline )
    class(commander_test_socket_server),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    call simple_end('**** SIMPLE_TEST_SOCKET_SERVER_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_server

end module simple_commanders_test_network
