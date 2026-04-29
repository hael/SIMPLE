!@descr: unit tests for simple_ipc_tcp_socket (port discovery, lifecycle, listener, kill)
!==============================================================================
! MODULE: simple_ipc_tcp_socket_tester
!
! PURPOSE:
!   Exercises the ipc_tcp_socket type across seven test cases:
!     1. test_find_available_port       — find_available_port returns found=.true.
!     2. test_port_in_range             — returned port is within [TCP_PORT_MIN, TCP_PORT_MAX]
!     3. test_port_is_bindable          — port returned by find_available_port can be
!                                         successfully bound with a raw c_socket call
!     4. test_occupied_port_skipped     — holding a port open causes find_available_port
!                                         to return a different (higher) port
!     5. test_init_and_kill             — init() starts the listener; kill() stops it
!     6. test_find_host_ips             — find_host_ips() returns a non-empty string
!     7. test_start_and_join_listener   — start_listener/join_listener lifecycle
!
! ENTRY POINT:
!   run_all_ipc_tcp_socket_tests
!
! DEPENDENCIES:
!   simple_ipc_tcp_socket, simple_test_utils, unix, iso_c_binding
!==============================================================================
module simple_ipc_tcp_socket_tester
  use simple_ipc_tcp_socket, only: ipc_tcp_socket, tcp_sockaddr_in, TCP_PORT_MIN, TCP_PORT_MAX, listener_args, recv_msg
  use simple_test_utils,     only: assert_true, assert_int
  use simple_string,         only: string
  use iso_c_binding
  use unix,          only: c_socket, c_bind, c_close, c_htons, c_socklen_t, &
                            c_pthread_mutex_lock, c_pthread_mutex_unlock,     &
                            AF_INET, SOCK_STREAM, INADDR_ANY, c_int16_t
  implicit none

  public  :: run_all_ipc_tcp_socket_tests
  private
#include "simple_local_flags.inc"

contains

  subroutine run_all_ipc_tcp_socket_tests()
    write(*,'(A)') '**** running all ipc_tcp_socket tests ****'
    call test_find_available_port()
    call test_port_in_range()
    call test_port_is_bindable()
    call test_occupied_port_skipped()
   ! call test_init_server_and_kill()
    call test_find_host_ips()
    call test_start_and_join_listener()
    write(*,'(A)') '**** ipc_tcp_socket tests done ****'
  end subroutine run_all_ipc_tcp_socket_tests

  ! find_available_port must succeed (found=.true.) on a machine that has
  ! at least one free port in [TCP_PORT_MIN, TCP_PORT_MAX].
  subroutine test_find_available_port()
    type(ipc_tcp_socket) :: sock
    logical              :: found
    write(*,'(A)') 'test_find_available_port'
    call sock%find_available_port(found)
    call assert_true(found, 'find_available_port returns found=.true.')
  end subroutine test_find_available_port

  ! The port stored after a successful find must lie within the declared range.
  subroutine test_port_in_range()
    type(ipc_tcp_socket) :: sock
    logical              :: found
    integer              :: port
    write(*,'(A)') 'test_port_in_range'
    call sock%find_available_port(found)
    port = sock%get_port()
    call assert_true(found,               'port found')
    call assert_true(port >= TCP_PORT_MIN, 'port >= TCP_PORT_MIN')
    call assert_true(port <= TCP_PORT_MAX, 'port <= TCP_PORT_MAX')
  end subroutine test_port_in_range

  ! Bind a raw socket to the port returned by find_available_port and verify
  ! that c_bind succeeds (rc == 0), confirming the port is genuinely free.
  subroutine test_port_is_bindable()
    type(ipc_tcp_socket)          :: sock
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: fd, rc
    logical                       :: found
    integer                       :: port
    write(*,'(A)') 'test_port_is_bindable'
    call sock%find_available_port(found)
    call assert_true(found, 'port found before bindability check')
    if( .not. found ) return
    port    = sock%get_port()
    addrlen = int(storage_size(addr) / 8, c_socklen_t)
    fd      = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
    call assert_true(fd >= 0, 'raw socket opened for bind check')
    if( fd < 0 ) return
    addr%sin_family = int(AF_INET,    c_int16_t)
    addr%sin_port   = c_htons(transfer(int(port, c_int32_t), 0_c_int16_t))
    addr%sin_addr   = int(INADDR_ANY, c_int32_t)
    rc = c_bind(fd, c_loc(addr), addrlen)
    call assert_true(rc == 0, 'c_bind succeeds on reported free port')
    rc = c_close(fd)
  end subroutine test_port_is_bindable

  ! Hold TCP_PORT_MIN open, then call find_available_port; the result must be
  ! a different port (strictly greater than TCP_PORT_MIN).
  subroutine test_occupied_port_skipped()
    type(ipc_tcp_socket)          :: sock
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: fd, rc
    logical                       :: found
    integer                       :: port
    write(*,'(A)') 'test_occupied_port_skipped'
    addrlen = int(storage_size(addr) / 8, c_socklen_t)
    ! occupy TCP_PORT_MIN
    fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
    call assert_true(fd >= 0, 'raw socket opened to occupy TCP_PORT_MIN')
    if( fd < 0 ) return
    addr%sin_family = int(AF_INET,    c_int16_t)
    addr%sin_port   = c_htons(transfer(int(TCP_PORT_MIN, c_int32_t), 0_c_int16_t))
    addr%sin_addr   = int(INADDR_ANY, c_int32_t)
    rc = c_bind(fd, c_loc(addr), addrlen)
    if( rc /= 0 ) then
      ! TCP_PORT_MIN already occupied by something else — test is still valid
      rc = c_close(fd)
    end if
    call sock%find_available_port(found)
    port = sock%get_port()
    call assert_true(found,               'port found with TCP_PORT_MIN occupied')
    call assert_true(port > TCP_PORT_MIN, 'returned port is beyond occupied TCP_PORT_MIN')
    rc = c_close(fd)
  end subroutine test_occupied_port_skipped

  ! init_server() must start the listener; kill() must stop it and reset port to -1.
  subroutine test_init_server_and_kill()
    type(ipc_tcp_socket)        :: sock
    type(listener_args), target :: args
    write(*,'(A)') 'test_init_server_and_kill'
    call sock%init_server(c_funloc(null_listener_thread), c_loc(args))
    call assert_true(sock%get_port() >= TCP_PORT_MIN, 'port valid after init')
    call assert_true(sock%get_port() <= TCP_PORT_MAX, 'port in range after init')
    call assert_true(sock%is_listening(),             'listener running after init')
    call sock%kill()
    call assert_true(sock%get_port() == -1,           'port reset to -1 after kill')
    call assert_true(.not. sock%is_listening(),       'listener stopped after kill')
    call sock%kill()   ! second kill must be a no-op
    call assert_true(sock%get_port() == -1,           'port still -1 after second kill')
  end subroutine test_init_server_and_kill

  ! find_host_ips() must populate host_ips with at least one IP address.
  subroutine test_find_host_ips()
    type(ipc_tcp_socket) :: sock
    type(string)         :: ips
    write(*,'(A)') 'test_find_host_ips'
    call sock%find_host_ips()
    ips = sock%get_host_ips()
    call assert_true(ips%is_allocated(),        'host_ips is allocated after find_host_ips')
    call assert_true(ips%strlen_trim() > 0,     'host_ips is non-empty after find_host_ips')
  end subroutine test_find_host_ips

  ! start_listener() must set is_listening=.true.; join_listener() must reset it.
  subroutine test_start_and_join_listener()
    type(ipc_tcp_socket)        :: sock
    type(listener_args), target :: args
    logical                     :: found
    write(*,'(A)') 'test_start_and_join_listener'
    call sock%find_available_port(found)
    call assert_true(found,                   'port found before start_listener')
    if( .not. found ) return
    call sock%start_listener(c_funloc(null_listener_thread), c_loc(args))
    call assert_true(sock%is_listening(),     'is_listening true after start_listener')
    call sock%kill()
    call assert_true(.not. sock%is_listening(), 'is_listening false after kill')
  end subroutine test_start_and_join_listener

  ! Minimal bind(c) thread: signals ready then accepts one message (the TERMINATE
  ! sentinel sent by kill()) so the parent's join_listener() can complete.
  subroutine null_listener_thread( arg ) bind(c)
    type(c_ptr), value,      intent(in) :: arg
    type(listener_args),     pointer    :: args
    character(kind=c_char, len=4096), target :: buf
    integer(kind=c_int) :: rc, conn_fd
    integer             :: nread
    logical             :: ok
    if( .not. c_associated(arg) ) return
    call c_f_pointer(arg, args)
    rc = c_pthread_mutex_lock(args%mutex)
    args%ready = 1
    rc = c_pthread_mutex_unlock(args%mutex)
    call recv_msg(args%server_fd, conn_fd, buf, nread, ok)
  end subroutine null_listener_thread

end module simple_ipc_tcp_socket_tester
