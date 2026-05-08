!==============================================================================
! MODULE: simple_ipc_tcp_socket_tester
!@descr: unit tests for IPC TCP client/helpers/server split modules
!
! PURPOSE:
!   Unit tests for the split IPC TCP modules:
!     - simple_ipc_tcp_socket_helpers
!     - simple_ipc_tcp_socket_server
!     - simple_ipc_tcp_socket_client
!
! ENTRYPOINT:
!   run_all_ipc_tcp_socket_tests
!==============================================================================
module simple_ipc_tcp_socket_tester
  use iso_c_binding
  use unix, only: c_pthread_mutex_init, c_pthread_mutex_destroy, c_pthread_mutex_lock, c_pthread_mutex_unlock, &
                  c_socket, c_connect, c_send, c_close, c_htons, c_inet_addr, &
                  AF_INET, SOCK_STREAM, c_socklen_t
  use simple_core_module_api, only: string
  use simple_test_utils, only: assert_true, assert_false, assert_int
  use simple_ipc_tcp_socket_helpers, only: TCP_MAX_MSG, POLLIN, POLLOUT, POLLERR, POLLHUP, POLLNVAL, &
                                           c_pollfd, poll_fds, fd_is_healthy, accept_connection
  use simple_ipc_tcp_socket_server, only: ipc_tcp_socket_server, listener_args, recv_msg, repl_msg, &
                                          tcp_sockaddr_in, TCP_PORT_MIN, TCP_PORT_MAX
  use simple_ipc_tcp_socket_client, only: ipc_tcp_socket_client, TCP_MAX_RETRIES, TCP_TIMEOUT_MS

  implicit none

  integer, parameter :: TEST_UNREACHABLE_PORT = 65000
  integer, parameter :: HELLO_BUF_LEN         = 12
  integer, parameter :: NO_LISTENER_BUF_LEN   = 11

  private
  public :: run_all_ipc_tcp_socket_tests

contains

  subroutine run_all_ipc_tcp_socket_tests()
    write(*,'(A)') '**** running all IPC TCP socket tests ****'
    ! helper module coverage
    call test_helpers_constants()
    call test_helpers_poll_fds_bounds_and_zero()
    call test_helpers_fd_is_healthy_invalid_fd()
    call test_helpers_accept_connection_invalid_fd()
#if !defined(_WIN32)
    ! server module coverage
    call test_server_find_ips_and_port_range()
#if !defined(__FreeBSD__)
  ! listener thread tests are problematic on macOS; skip them
  call test_server_new_kill_lifecycle()
  call test_server_start_join_lifecycle()
  ! client module and integration coverage
  call test_client_server_roundtrip()
  call test_client_send_recv_without_server_fails()
#endif
#endif
    write(*,'(A)') '**** IPC TCP socket tests done ****'
  end subroutine run_all_ipc_tcp_socket_tests

  subroutine test_helpers_constants()
    write(*,'(A)') 'test_helpers_constants'
    call assert_int(1460, TCP_MAX_MSG, 'TCP_MAX_MSG should be MTU payload (1460)')
    call assert_true(POLLIN   > 0_c_short, 'POLLIN should be positive')
    call assert_true(POLLOUT  > 0_c_short, 'POLLOUT should be positive')
    call assert_true(POLLERR  > 0_c_short, 'POLLERR should be positive')
    call assert_true(POLLHUP  > 0_c_short, 'POLLHUP should be positive')
    call assert_true(POLLNVAL > 0_c_short, 'POLLNVAL should be positive')
  end subroutine test_helpers_constants

  subroutine test_helpers_poll_fds_bounds_and_zero()
    type(c_pollfd) :: fds(2)
    integer        :: nready
    write(*,'(A)') 'test_helpers_poll_fds_bounds_and_zero'
    call poll_fds(fds, -1, 0, nready)
    call assert_int(-1, nready, 'poll_fds should reject n < 0')
    call poll_fds(fds, 3, 0, nready)
    call assert_int(-1, nready, 'poll_fds should reject n > size(fds)')
    call poll_fds(fds, 0, 0, nready)
    call assert_int(0, nready, 'poll_fds should return 0 for n == 0')
  end subroutine test_helpers_poll_fds_bounds_and_zero

  subroutine test_helpers_fd_is_healthy_invalid_fd()
    logical :: ok
    write(*,'(A)') 'test_helpers_fd_is_healthy_invalid_fd'
    ok = fd_is_healthy(-1_c_int)
    call assert_false(ok, 'fd_is_healthy should be false for invalid fd')
  end subroutine test_helpers_fd_is_healthy_invalid_fd

  subroutine test_helpers_accept_connection_invalid_fd()
    integer(c_int) :: conn_fd
    write(*,'(A)') 'test_helpers_accept_connection_invalid_fd'
    conn_fd = accept_connection(-1_c_int)
    call assert_true(conn_fd < 0_c_int, 'accept_connection should fail for invalid listener fd')
  end subroutine test_helpers_accept_connection_invalid_fd

  subroutine test_server_find_ips_and_port_range()
    type(ipc_tcp_socket_server) :: server
    type(string)                :: ips
    write(*,'(A)') 'test_server_find_ips_and_port_range'
    call server%find_server_ips()
    ips = server%get_server_ips()
    if( ips%is_allocated() ) then
      call assert_false(index(ips%to_char(), '127.') > 0, 'find_server_ips should skip loopback 127.* entries')
    else
      call assert_true(.true., 'find_server_ips may be empty in minimal/container environments')
    end if
    call assert_int(-1, server%get_port(), 'new server object should have port=-1 before start')
    call assert_true(TCP_PORT_MIN < TCP_PORT_MAX, 'server port range constants should be valid')
  end subroutine test_server_find_ips_and_port_range

  subroutine test_server_new_kill_lifecycle()
    type(ipc_tcp_socket_server)  :: server
    type(listener_args), target  :: largs
    write(*,'(A)') 'test_server_new_kill_lifecycle'
    call init_listener_args(largs)
    call server%new(c_funloc(listener_ready_once), c_loc(largs))
    call assert_true(server%is_listening(), 'server should report listening after new()')
    call assert_port_in_range(server%get_port(), 'server bound port should be inside configured range')
    call server%kill()
    call assert_false(server%is_listening(), 'server should not report listening after kill()')
    call assert_int(-1, server%get_port(), 'server port should reset to -1 after kill()')
    call destroy_listener_args(largs)
  end subroutine test_server_new_kill_lifecycle

  subroutine test_server_start_join_lifecycle()
    type(ipc_tcp_socket_server)  :: server
    type(listener_args), target  :: largs
    integer :: port
    write(*,'(A)') 'test_server_start_join_lifecycle'
    call init_listener_args(largs)
    call server%start_listener(c_funloc(listener_ready_once), c_loc(largs))
    call assert_true(server%is_listening(), 'server should report listening after start_listener()')
    port = server%get_port()
    call assert_port_in_range(port, 'start_listener() should bind inside configured port range')
    call poke_server_port(port)
    call server%join_listener()
    call assert_false(server%is_listening(), 'server should not report listening after join_listener()')
    call assert_int(-1, server%get_port(), 'server port should reset to -1 after join_listener()')
    call destroy_listener_args(largs)
  end subroutine test_server_start_join_lifecycle

  subroutine test_client_server_roundtrip()
    type(ipc_tcp_socket_server)  :: server
    type(listener_args), target  :: largs
    type(ipc_tcp_socket_client)  :: client
    character(kind=c_char), target :: rcv_buffer(TCP_MAX_MSG)
    character(kind=c_char, len=:), allocatable, target :: snd_buffer
    integer :: nread
    logical :: sent
    write(*,'(A)') 'test_client_server_roundtrip'
    call init_listener_args(largs)
    call server%new(c_funloc(listener_echo_ok_once), c_loc(largs))
    call assert_true(server%is_listening(), 'server should be listening before client roundtrip')
    call client%new(string('127.0.0.1'), server%get_port())
    allocate(character(kind=c_char, len=HELLO_BUF_LEN) :: snd_buffer)
    snd_buffer = 'hello-server'
    call client%send_recv_msg(snd_buffer, rcv_buffer, sent, nread)
    call assert_true(sent, 'client send_recv_msg should succeed against local echo listener')
    call assert_true(nread > 0, 'client send_recv_msg should read a non-empty reply')
    call client%kill()
    call server%kill()
    call destroy_listener_args(largs)
  end subroutine test_client_server_roundtrip

  subroutine test_client_send_recv_without_server_fails()
    type(ipc_tcp_socket_client)  :: client
    character(kind=c_char), target :: rcv_buffer(TCP_MAX_MSG)
    character(kind=c_char, len=:), allocatable, target :: snd_buffer
    integer :: nread
    logical :: sent
    write(*,'(A)') 'test_client_send_recv_without_server_fails'
    call assert_true(TCP_MAX_RETRIES > 0, 'client retry constant should be positive')
    call assert_true(TCP_TIMEOUT_MS  > 0, 'client timeout constant should be positive')
    call client%new(string('127.0.0.1'), TEST_UNREACHABLE_PORT)
    allocate(character(kind=c_char, len=NO_LISTENER_BUF_LEN) :: snd_buffer)
    snd_buffer = 'no-listener'
    call client%send_recv_msg(snd_buffer, rcv_buffer, sent, nread)
    call assert_false(sent, 'client send_recv_msg should fail when no listener is running')
    call assert_int(0, nread, 'nread should remain 0 when send_recv_msg fails')
    call client%kill()
  end subroutine test_client_send_recv_without_server_fails

  subroutine init_listener_args(largs)
    type(listener_args), intent(inout) :: largs
    integer(c_int) :: rc
    write(*,'(A)') 'init_listener_args'
    largs%fd       = -1_c_int
    largs%ready    = 0_c_int
    largs%data_ptr = c_null_ptr
    rc = c_pthread_mutex_init(largs%mutex, c_null_ptr)
    call assert_int(0, int(rc), 'listener mutex init should succeed')
  end subroutine init_listener_args

  subroutine destroy_listener_args(largs)
    type(listener_args), intent(inout) :: largs
    integer(c_int) :: rc
    write(*,'(A)') 'destroy_listener_args'
    rc = c_pthread_mutex_destroy(largs%mutex)
    call assert_int(0, int(rc), 'listener mutex destroy should succeed')
  end subroutine destroy_listener_args

  subroutine assert_port_in_range(port, msg)
    integer,          intent(in) :: port
    character(len=*), intent(in) :: msg
    write(*,'(A)') 'assert_port_in_range'
    call assert_true(port >= TCP_PORT_MIN .and. port <= TCP_PORT_MAX, msg)
  end subroutine assert_port_in_range

  subroutine poke_server_port(port)
    integer, intent(in) :: port
    type(tcp_sockaddr_in), target :: addr
    integer(c_int), target :: ping
    integer(c_int) :: fd
    integer(c_int) :: rc
    integer(c_socklen_t) :: addrlen
    write(*,'(A)') 'poke_server_port'
    fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
    call assert_true(fd >= 0_c_int, 'poke helper should open socket')
    if( fd < 0_c_int ) return
    addrlen         = int(storage_size(addr) / 8, c_socklen_t)
    addr%sin_family = int(AF_INET, c_int16_t)
    addr%sin_port   = c_htons(transfer(int(port, c_int32_t), 0_c_int16_t))
    addr%sin_addr   = int(c_inet_addr('127.0.0.1' // c_null_char), c_int32_t)
    rc = c_connect(fd, c_loc(addr), addrlen)
    call assert_int(0, int(rc), 'poke helper should connect to listener')
    if( rc == 0_c_int ) then
      ping = 1_c_int
      rc   = int(c_send(fd, c_loc(ping), int(storage_size(ping)/8, c_size_t), 0_c_int), c_int)
      call assert_true(rc >= 0_c_int, 'poke helper should send one wake payload')
    end if
    rc = c_close(fd)
    call assert_int(0, int(rc), 'poke helper should close socket')
  end subroutine poke_server_port

  subroutine listener_ready_once(arg_ptr) bind(c)
    type(c_ptr), value :: arg_ptr
    type(listener_args), pointer :: args
    character(kind=c_char, len=TCP_MAX_MSG), target :: msg
    integer(c_int) :: conn_fd
    integer :: nread
    logical :: ok
    integer(c_int) :: rc
    write(*,'(A)') 'listener_ready_once'
    if( .not. c_associated(arg_ptr) ) return
    call c_f_pointer(arg_ptr, args)
    rc = c_pthread_mutex_lock(args%mutex)
    args%ready = 1_c_int
    rc = c_pthread_mutex_unlock(args%mutex)
    msg     = ''
    conn_fd = -1_c_int
    call recv_msg(args%fd, conn_fd, msg, nread, ok)
  end subroutine listener_ready_once

  subroutine listener_echo_ok_once(arg_ptr) bind(c)
    type(c_ptr), value :: arg_ptr
    type(listener_args), pointer :: args
    character(kind=c_char, len=TCP_MAX_MSG), target :: in_msg
    character(kind=c_char, len=TCP_MAX_MSG), target :: out_msg
    integer(c_int) :: conn_fd
    integer :: nread, nwrite
    logical :: ok
    integer(c_int) :: rc
    write(*,'(A)') 'listener_echo_ok_once'
    if( .not. c_associated(arg_ptr) ) return
    call c_f_pointer(arg_ptr, args)
    rc = c_pthread_mutex_lock(args%mutex)
    args%ready = 1_c_int
    rc = c_pthread_mutex_unlock(args%mutex)
    in_msg  = ''
    out_msg = ''
    out_msg(1:2) = 'OK'
    conn_fd = -1_c_int
    call recv_msg(args%fd, conn_fd, in_msg, nread, ok, close_after_read=.false.)
    if( ok .and. conn_fd >= 0_c_int ) then
      call repl_msg(conn_fd, out_msg, nwrite, ok)
    end if
  end subroutine listener_echo_ok_once

end module simple_ipc_tcp_socket_tester
