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
    use simple_socket_comm
    class(commander_test_socket_client), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    type(simple_socket) :: socket
    write(*,*) "Socket client test"
    call socket%open
    call socket%send("TEST MESSAGE FROM THE CLIENT")
    call socket%close
    write(*,*) "Sent message"
    call simple_end('**** SIMPLE_TEST_SOCKET_CLIENT_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_client

subroutine exec_test_socket_comm_distr( self, cline )
    use, intrinsic :: iso_c_binding
    use simple_distr_comm
    class(commander_test_socket_comm_distr),    intent(inout) :: self
    class(cmdline),                             intent(inout) :: cline
    type(distr_comm) :: comm
    write(*,*) "START"
    call comm%init()
    call sleep(10)
    call simple_end('**** SIMPLE_TEST_SOCKET_COMM_DISTR_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_comm_distr

subroutine exec_test_socket_io( self, cline )
    use simple_socket_comm
    use simple_string_utils
    use unix, only : c_pthread_t, c_pthread_create 
    use, intrinsic :: iso_c_binding
    class(commander_test_socket_io),    intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    type(c_pthread_t)   :: server_thread
    type(simple_socket) :: client_socket
    integer             :: server_rc, i
    ! start server thread
    server_rc = c_pthread_create(server_thread, c_null_ptr, c_funloc(socket_server), c_null_ptr)
    ! send a client message every 5 seconds for 5 iterations
    do i=1,5
       call sleep(10)
       call client_socket%open
       call client_socket%send("TEST MESSAGE FROM THE CLIENT. ITERATION : "//int2str(i))
       call client_socket%close
    end do
    call simple_end('**** SIMPLE_TEST_SOCKET_IO_WORKFLOW NORMAL STOP ****')

    contains

       subroutine socket_server() bind(c)
          type(simple_socket) :: socket
          integer             :: fd
          write(*,*) "Starting socket server thread"
          call socket%open()
          call socket%bind_any()
          call socket%listen()
          write(*,*) "Socket server listening"
          do
             call socket%accept(fd)
             call socket%read(fd)
             call socket%close(fd)
          end do
          call socket%close
       end subroutine socket_server

end subroutine exec_test_socket_io

subroutine exec_test_socket_server( self, cline )
    use simple_socket_comm
    use json_kinds
    use json_module
    use unix, only : c_pthread_t, c_pthread_mutex_t 
    use unix, only : c_pthread_create, c_pthread_join
    use unix, only : c_pthread_mutex_init, c_pthread_mutex_destroy
    use unix, only : c_pthread_mutex_lock, c_pthread_mutex_unlock
    use, intrinsic :: iso_c_binding
    class(commander_test_socket_server), intent(inout) :: self
    class(cmdline),                      intent(inout) :: cline
    type(simple_socket)                    :: socket
    integer                                :: fd
    write(*,*) "Socket server test. Waiting for client message"
    call socket%open
    call socket%set_options
    call socket%bind_any
    call socket%listen
    do
       call socket%accept(fd)
       call socket%read(fd)
       call socket%close(fd)
    end do
    call socket%close
    call simple_end('**** SIMPLE_TEST_SOCKET_SERVER_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_socket_server

end module simple_commanders_test_network
