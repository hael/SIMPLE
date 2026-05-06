!==============================================================================
! MODULE: simple_ipc_tcp_socket_server
!
! PURPOSE:
!   TCP listener wrapper used by SIMPLE components. It owns the server socket,
!   starts a listener pthread, and provides thin request/reply helpers for
!   accepted connections.
!
! DESIGN NOTES:
!   - listener_args is passed by pointer to the listener pthread.
!   - The mutex in listener_args is owned by the caller of this module.
!   - kill() unblocks accept() using a localhost sentinel connection, then joins.
!==============================================================================
module simple_ipc_tcp_socket_server
  use iso_c_binding
  use unix,          only: c_pthread_t, c_pthread_mutex_t,                    &
                            c_pthread_create, c_pthread_join,                  &
                            c_pthread_mutex_lock, c_pthread_mutex_unlock,      &
                            c_socket, c_bind, c_listen, c_accept,             &
                            c_connect, c_send, c_setsockopt, c_close, c_read, &
                            c_usleep, c_popen, c_fgets, c_pclose,              &
                            c_htons, c_inet_addr,                              &
                            AF_INET, SOCK_STREAM, INADDR_ANY,                  &
                            SOL_SOCKET, SO_REUSEADDR, c_socklen_t
  use simple_error,                  only: simple_exception
  use simple_core_module_api,        only: logfhandle, string
  use simple_ipc_tcp_socket_helpers, only: c_pollfd, poll_fds

  implicit none

  integer, parameter, public :: TCP_PORT_MIN = 39000
  integer, parameter, public :: TCP_PORT_MAX = 39999

  public :: ipc_tcp_socket_server
  public :: tcp_sockaddr_in
  public :: listener_args
  public :: c_pollfd
  public :: recv_msg
  public :: repl_msg
  public :: poll_fds

  ! Correct C layout for struct sockaddr_in (16 bytes on LP64 Linux/BSD).
  ! unix_netdb's c_sockaddr_in uses c_signed_char for sin_family (1 byte)
  ! instead of the correct uint16_t (2 bytes), corrupting the struct layout.
  type, bind(c) :: tcp_sockaddr_in
    integer(kind=c_int16_t) :: sin_family = 0_c_int16_t
    integer(kind=c_int16_t) :: sin_port   = 0_c_int16_t
    integer(kind=c_int32_t) :: sin_addr   = 0_c_int32_t
    character(kind=c_char)  :: sin_zero(8)
  end type tcp_sockaddr_in

  private
#include "simple_local_flags.inc"

  integer, parameter :: LISTEN_BACKLOG = 5

  ! POD struct passed by pointer to the listener pthread.
  type, bind(c) :: listener_args
    integer(kind=c_int)     :: fd    = -1
    integer(kind=c_int)     :: ready = 0         ! set to 1 once listen() succeeds
    type(c_pthread_mutex_t) :: mutex
    type(c_ptr)             :: data_ptr = c_null_ptr  ! optional opaque payload for listener
  end type listener_args

  type :: ipc_tcp_socket_server
    private
    integer               :: port        = -1
    integer(kind=c_int)   :: fd          = -1
    type(string)          :: server_ips
    type(c_pthread_t)     :: thread
    type(listener_args), pointer :: listener_args => null()
    logical               :: thread_started = .false.
    logical               :: listening   = .false.
  contains
    procedure :: new
    procedure :: kill
    procedure :: start_listener
    procedure :: join_listener
    procedure :: find_server_ips
    procedure :: get_port
    procedure :: get_server_ips
    procedure :: is_listening
  end type ipc_tcp_socket_server

  contains

  !> Reinitialise server object and start a fresh listener thread.
  subroutine new( self, thread_funloc, thread_args_ptr )
    class(ipc_tcp_socket_server), intent(inout) :: self
    type(c_ptr), value,    intent(in)    :: thread_funloc
    type(c_ptr), value,    intent(in)    :: thread_args_ptr
    call self%kill()
    call self%find_server_ips()
    call self%start_listener(thread_funloc, thread_args_ptr)
  end subroutine new

  !> Stop listener thread and release socket/IP state.
  subroutine kill(self)
    class(ipc_tcp_socket_server),    intent(inout) :: self
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: fd, rc
    integer,               target :: kill_msg
    if( self%port < 0 ) return
    ! unblock the listener's accept() by connecting and sending the KILL sentinel
    fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
    if( fd >= 0 ) then
      addrlen         = int(storage_size(addr) / 8, c_socklen_t)
      addr%sin_family = int(AF_INET, c_int16_t)
      addr%sin_port   = c_htons(transfer(int(self%port, c_int32_t), 0_c_int16_t))
      addr%sin_addr   = int(c_inet_addr('127.0.0.1' // c_null_char), c_int32_t)
      rc = c_connect(fd, c_loc(addr), addrlen)
      if( rc == 0 ) then
        kill_msg = 1  ! sentinel value expected by listener-side dispatcher
        rc  = int(c_send(fd, c_loc(kill_msg), int(storage_size(kill_msg)/8, c_size_t), 0_c_int))
      end if
      rc = c_close(fd)
    end if
    call self%join_listener()
    self%port = -1
    call self%server_ips%kill()
  end subroutine kill

  !> Bind an available port, call listen(), and launch the listener pthread.
  subroutine start_listener( self, thread_funloc, thread_args_ptr )
    class(ipc_tcp_socket_server), intent(inout) :: self
    type(c_ptr), value,         intent(in) :: thread_funloc
    type(c_ptr), value,         intent(in) :: thread_args_ptr
    type(listener_args),  pointer :: args
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: rc
    integer(kind=c_int), target   :: optval
    integer                       :: i, iport
    if( .not. c_associated(thread_args_ptr) ) return
    if( .not. c_associated(thread_funloc)   ) return
    call c_f_pointer(thread_args_ptr, args)
    ! Mutex is initialised by the owner of listener_args before calling start_listener.
    self%port = -1
    addrlen   = int(storage_size(addr) / 8, c_socklen_t)
    do iport = TCP_PORT_MIN, TCP_PORT_MAX
      ! open the server socket
      self%fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
      if( self%fd < 0 ) THROW_HARD('start_listener: c_socket failed')
      args%fd = self%fd
      optval = 1
      rc = c_setsockopt(self%fd, SOL_SOCKET, SO_REUSEADDR, c_loc(optval), &
            int(storage_size(optval)/8, c_socklen_t))
      addr%sin_family = int(AF_INET,    c_int16_t)
      addr%sin_port   = c_htons(transfer(int(iport, c_int32_t), 0_c_int16_t))
      addr%sin_addr   = int(INADDR_ANY, c_int32_t)
      rc = c_bind(self%fd, c_loc(addr), addrlen)
      if( rc == 0 ) then
        self%port = iport
        exit
      end if
      rc = c_close(self%fd)
      self%fd = -1
    end do
    if( self%port < 0 ) THROW_HARD('start_listener: failed to find available TCP port')
    write(logfhandle,'(A,I0)')'>>> IPC_TCP_SOCKET available port found: ', self%port
    rc = c_listen(self%fd, int(LISTEN_BACKLOG, c_int))
    if( rc /= 0 ) THROW_HARD('start_listener: c_listen failed')
    args%ready = 0
    args%fd    = self%fd
    ! spawn listener thread, passing args by pointer
    rc = c_pthread_create(self%thread, c_null_ptr, thread_funloc, thread_args_ptr)
    if( rc /= 0 ) THROW_HARD('start_listener: c_pthread_create failed')
    self%thread_started = .true.
    ! wait up to 2 s for the thread to signal it is ready
    do i = 1, 20
      rc = c_pthread_mutex_lock(args%mutex)
      if( args%ready /= 0 ) then
        self%listening = .true.
        rc = c_pthread_mutex_unlock(args%mutex)
        exit
      end if
      rc = c_pthread_mutex_unlock(args%mutex)
      rc = c_usleep(100000)  ! 100 ms
    end do
    if( self%listening ) then
      write(logfhandle,'(A,I0)')'>>> IPC_TCP_SOCKET listener running on port ', self%port
    else
      write(logfhandle,'(A)')'>>> IPC_TCP_SOCKET WARNING: listener ready flag not set'
    end if
  end subroutine start_listener

  !> Join listener thread (if started) and close server fd.
  subroutine join_listener( self )
    class(ipc_tcp_socket_server), intent(inout) :: self
    type(c_ptr)         :: retval
    integer(kind=c_int) :: rc
    if( .not. self%thread_started ) return
    rc = c_pthread_join(self%thread, retval)
    self%thread_started = .false.
    if( self%fd >= 0 ) then
      rc = c_close(self%fd)
      self%fd = -1
    end if
    self%port = -1
    ! The mutex is owned by the caller (e.g. persistent_worker_server) and must be
    ! destroyed by the owner after join_listener returns, not here.
    self%listening = .false.
  end subroutine join_listener

  !> Return .true. only after the listener thread has signaled ready.
  logical function is_listening( self )
    class(ipc_tcp_socket_server), intent(in) :: self
    is_listening = self%listening
  end function is_listening

  !> Return bound TCP listen port, or -1 when not active.
  integer function get_port( self )
    class(ipc_tcp_socket_server), intent(in) :: self
    get_port = self%port
  end function get_port

  !> Return comma-separated server IP list discovered by find_server_ips().
  function get_server_ips( self ) result( ips )
    class(ipc_tcp_socket_server), intent(in) :: self
    type(string)                              :: ips
    ips = self%server_ips
  end function get_server_ips

  !> Run 'hostname -I', collect all non-loopback tokens, and store them as a
  !> comma-separated list in self%server_ips.
  subroutine find_server_ips( self )
    class(ipc_tcp_socket_server), intent(inout) :: self
    character(kind=c_char, len=16) :: cmd, mode
    character(kind=c_char)         :: buf(1024)
    type(c_ptr)                    :: pipe, ret
    integer(kind=c_int)            :: rc
    integer                        :: i, token_start, buflen, tok_len
    character(len=1024)            :: fbuf
    character(len=64)              :: token
    logical                        :: first
    cmd  = 'hostname -I' // c_null_char
    mode = 'r'           // c_null_char
    self%server_ips = string('')
    pipe = c_popen(cmd, mode)
    if( .not. c_associated(pipe) ) return
    ret = c_fgets(buf(1), int(size(buf), c_int), pipe)
    rc  = c_pclose(pipe)
    if( .not. c_associated(ret) ) return
    ! convert null-terminated C buffer to a Fortran character variable
    fbuf   = ''
    buflen = 0
    do i = 1, size(buf)
      if( buf(i) == c_null_char ) exit
      fbuf(i:i) = buf(i)
      buflen    = i
    end do
    ! tokenise on whitespace; skip loopback (127.*) addresses
    first       = .true.
    token_start = 1
    do i = 1, buflen + 1
      if( i > buflen          .or. fbuf(i:i) == ' '     &
                              .or. fbuf(i:i) == char(10) &
                              .or. fbuf(i:i) == char(13) ) then
        tok_len = i - token_start
        if( tok_len > 0 ) then
          token = fbuf(token_start : token_start + tok_len - 1)
          if( token(1:4) /= '127.' ) then
            if( first ) then
              self%server_ips = string(trim(token))
              first         = .false.
            else
              self%server_ips = string(self%server_ips%to_char() // ',' // trim(token))
            end if
          end if
        end if
        token_start = i + 1
      end if
    end do
    write(logfhandle,'(A,A)') '>>> IPC_TCP_SOCKET host IPs: ', self%server_ips%to_char()
  end subroutine find_server_ips

  !> Accept one connection on fd, read into buf, close the connection.
  !> ok=.false. if accept() failed — caller should exit the accept loop.
  !> Accept a single connection and read up to len(buf) bytes.
  !> Optionally closes the accepted fd before return (default: close).
  subroutine recv_msg( fd, conn_fd, buf, nread, ok, close_after_read )
    integer(kind=c_int),           intent(in)         :: fd
    character(kind=c_char, len=*), intent(inout), target :: buf
    integer(kind=c_int),           intent(inout)         :: conn_fd
    integer,                       intent(out)        :: nread
    logical,                       intent(out)        :: ok
    logical, optional,             intent(in)         :: close_after_read
    integer(kind=c_int)    :: rc
    integer(kind=c_size_t) :: nr
    logical                 :: l_close_conn
    l_close_conn = .true.
    if( present(close_after_read) ) l_close_conn = close_after_read
    ok    = .false.
    nread = 0
    conn_fd = c_accept(fd, c_null_ptr, int(0, c_socklen_t))
    if( conn_fd < 0 ) return
    nr    = c_read(conn_fd, c_loc(buf), int(len(buf), c_size_t))
    if( l_close_conn ) rc = c_close(conn_fd)
    if( nr > int(len(buf), c_size_t) ) return  ! c_read error (ssize_t -1 wraps to HUGE)
    nread = int(nr)
    ok    = .true.
  end subroutine recv_msg

  !> Send len(buf) bytes to an accepted connection fd.
  !> Optionally closes the fd after write (default: close).
  subroutine repl_msg( conn_fd, buf, nread, ok, close_after_write )
    integer(kind=c_int),           intent(in)         :: conn_fd
    character(kind=c_char, len=*), intent(inout), target :: buf
    integer,                       intent(out)        :: nread
    logical,                       intent(out)        :: ok
    logical, optional,             intent(in)         :: close_after_write
    integer(kind=c_int)                               :: rc
    logical                                           :: l_close_conn
    l_close_conn = .true.
    if( present(close_after_write) ) l_close_conn = close_after_write
    ok    = .false.
    rc    = int(c_send(conn_fd, c_loc(buf), int(len(buf), c_size_t), 0_c_int))
    nread = rc
    ok    = (rc >= 0)
    if( l_close_conn ) rc = c_close(conn_fd)
  end subroutine repl_msg


end module simple_ipc_tcp_socket_server