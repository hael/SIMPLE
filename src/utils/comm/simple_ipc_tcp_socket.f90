module simple_ipc_tcp_socket
  use iso_c_binding
  use unix,          only: c_pthread_t, c_pthread_mutex_t,                    &
                            c_pthread_create, c_pthread_join,                  &
                            c_pthread_mutex_init, c_pthread_mutex_destroy,     &
                            c_pthread_mutex_lock, c_pthread_mutex_unlock,      &
                            c_socket, c_bind, c_listen, c_accept,             &
                            c_connect, c_send, c_setsockopt, c_close, c_read, &
                            c_usleep, c_useconds_t, c_popen, c_fgets, c_pclose, &
                            c_htons, c_inet_addr,                              &
                            AF_INET, SOCK_STREAM, INADDR_ANY,                  &
                            SOL_SOCKET, SO_REUSEADDR, c_socklen_t
  use simple_error,           only: simple_exception
  use simple_core_module_api, only: logfhandle, string

  implicit none

  integer, parameter, public :: TCP_PORT_MIN = 39000
  integer, parameter, public :: TCP_PORT_MAX = 39999
  integer, parameter, public :: TCP_MAX_MSG  = 1460 ! 1500 MTU - 20 IP header - 20 TCP header

  public :: ipc_tcp_socket
  public :: tcp_sockaddr_in
  public :: listener_args
  public :: recv_msg
  public :: repl_msg

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

  integer, parameter :: LISTEN_BACKLOG  = 5
  integer(kind=c_int), parameter :: SO_SNDTIMEO = 21  ! Linux SO_SNDTIMEO
  integer(kind=c_int), parameter :: SO_RCVTIMEO = 20  ! Linux SO_RCVTIMEO

  ! POSIX struct timeval used with SO_SNDTIMEO via setsockopt.
  type, bind(c) :: c_timeval
    integer(kind=c_long) :: tv_sec  = 0_c_long
    integer(kind=c_long) :: tv_usec = 0_c_long
  end type c_timeval

  ! POD struct passed by pointer to the listener pthread.
  type, bind(c) :: listener_args
    integer(kind=c_int)          :: port      = 0
    integer(kind=c_int)          :: server_fd = -1
    integer(kind=c_int)          :: ready     = 0   ! set to 1 once listen() succeeds
    type(c_pthread_mutex_t)      :: mutex
    type(c_ptr)                  :: data_ptr  = c_null_ptr  ! data_ptr can be used to point to a Fortran-derived type instance, allowing
  end type listener_args

  type :: ipc_tcp_socket
    private
    integer               :: port      = -1
    integer(kind=c_int)   :: server_fd = -1
    type(string)          :: host_ips
    type(string)          :: server_ip
    type(c_pthread_t)     :: thread
    type(c_ptr)           :: thread_argsloc
    type(listener_args), pointer :: listener_args => null()
    logical               :: listening = .false.
  contains
    procedure :: init_server
    procedure :: init_client
    procedure :: kill
    procedure :: start_listener
    procedure :: join_listener
    procedure :: find_available_port
    procedure :: find_available_server
    procedure :: find_host_ips
    procedure :: send_recv_msg
    procedure :: get_port
    procedure :: get_host_ips
    procedure :: is_listening
  end type ipc_tcp_socket

  contains

  subroutine init_server( self, thread_funloc, thread_args_ptr )
    class(ipc_tcp_socket), intent(inout) :: self
    type(c_ptr), value,    intent(in)    :: thread_funloc
    type(c_ptr), value,    intent(in)    :: thread_args_ptr
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: rc
    integer(kind=c_int), target   :: optval
    integer                       :: i
    logical :: found
    call self%find_host_ips()
    call self%find_available_port(found)
    if( .not.found ) THROW_HARD('failed to find available TCP port')
    write(logfhandle,'(A,I0)')'>>> IPC_TCP_SOCKET available port found: ', self%port
    call self%start_listener(thread_funloc, thread_args_ptr)
  end subroutine init_server

  subroutine init_client(self, server_list, port)
    class(ipc_tcp_socket), intent(inout) :: self
    type(string),          intent(in)    :: server_list
    integer,               intent(in)    :: port
    type(string) :: found_ip
    logical      :: found
    self%port = port
    call self%find_available_server(server_list, found)
  end subroutine init_client

  subroutine kill(self)
    class(ipc_tcp_socket),    intent(inout) :: self
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
        kill_msg = 1 ! needs to be WORKER_TERMINATE_MSG
        rc  = int(c_send(fd, c_loc(kill_msg), int(1, c_size_t), 0_c_int))
      end if
      rc = c_close(fd)
    end if
    call self%join_listener()
    self%port = -1
    call self%host_ips%kill()
  end subroutine kill

  !> Bind the found port, start listen(), and launch the listener pthread.
  subroutine start_listener( self, thread_funloc, thread_args_ptr )
    class(ipc_tcp_socket), intent(inout) :: self
    type(c_ptr), value,         intent(in) :: thread_funloc
    type(c_ptr), value,         intent(in) :: thread_args_ptr
    type(listener_args),  pointer :: args
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: rc
    integer(kind=c_int), target   :: optval
    integer                       :: i
    if( .not. c_associated(thread_args_ptr) ) return
    if( .not. c_associated(thread_funloc)   ) return
    call c_f_pointer(thread_args_ptr, args)
    ! Mutex is initialised by the owner of listener_args before calling start_listener.
    args%port     = self%port
   ! self%listener_args%data_ptr = thread_data_ptr
    ! open the server socket
    self%server_fd   = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
    if( self%server_fd < 0 ) THROW_HARD('start_listener: c_socket failed')
    args%server_fd = self%server_fd
    ! SO_REUSEADDR so the port can be rebound immediately after find_available_port
    optval = 1
    rc = c_setsockopt(self%server_fd, SOL_SOCKET, SO_REUSEADDR, &
                      c_loc(optval), int(storage_size(optval)/8, c_socklen_t))
    addrlen         = int(storage_size(addr) / 8, c_socklen_t)
    addr%sin_family = int(AF_INET, c_int16_t)
    addr%sin_port   = c_htons(transfer(int(self%port, c_int32_t), 0_c_int16_t))
    addr%sin_addr   = int(INADDR_ANY, c_int32_t)
    rc = c_bind(self%server_fd, c_loc(addr), addrlen)
    if( rc /= 0 ) THROW_HARD('start_listener: c_bind failed')
    rc = c_listen(self%server_fd, int(LISTEN_BACKLOG, c_int))
    if( rc /= 0 ) THROW_HARD('start_listener: c_listen failed')
    ! spawn listener thread, passing args by pointer
    rc = c_pthread_create(self%thread, c_null_ptr, thread_funloc, thread_args_ptr)
    if( rc /= 0 ) THROW_HARD('start_listener: c_pthread_create failed')
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

  !> Join the listener thread and release resources.
  subroutine join_listener( self )
    class(ipc_tcp_socket), intent(inout) :: self
    type(listener_args), pointer :: args
    type(c_ptr)         :: retval
    integer(kind=c_int) :: rc
    if( .not. self%listening ) return
    rc = c_pthread_join(self%thread, retval)
    if( self%server_fd >= 0 ) then
      rc = c_close(self%server_fd)
      self%server_fd = -1
    end if
    ! The mutex is owned by the caller (e.g. persistent_worker_server) and must be
    ! destroyed by the owner after join_listener returns, not here.
    self%listening = .false.
  end subroutine join_listener

  logical function is_listening( self )
    class(ipc_tcp_socket), intent(in) :: self
    is_listening = self%listening
  end function is_listening

  integer function get_port( self )
    class(ipc_tcp_socket), intent(in) :: self
    get_port = self%port
  end function get_port

  function get_host_ips( self ) result( ips )
    class(ipc_tcp_socket), intent(in) :: self
    type(string)                      :: ips
    ips = self%host_ips
  end function get_host_ips

  !> Run 'hostname -I', collect all space-delimited tokens that are not
  !> loopback (127.*), and store them comma-separated in self%host_ips.
  subroutine find_host_ips( self )
    class(ipc_tcp_socket), intent(inout) :: self
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
    self%host_ips = string('')
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
              self%host_ips = string(trim(token))
              first         = .false.
            else
              self%host_ips = string(self%host_ips%to_char() // ',' // trim(token))
            end if
          end if
        end if
        token_start = i + 1
      end if
    end do
    write(logfhandle,'(A,A)')'>>> IPC_TCP_SOCKET host IPs: ', self%host_ips%to_char()
  end subroutine find_host_ips

  !> Scan TCP_PORT_MIN..TCP_PORT_MAX and return the first port that can be
  !> bound (i.e. not already in use). Sets found=.false. if none are free.
  subroutine find_available_port( self, found )
    class(ipc_tcp_socket), intent(inout) :: self
    logical,               intent(out)   :: found
    type(tcp_sockaddr_in), target :: addr
    integer(kind=c_socklen_t)     :: addrlen
    integer(kind=c_int)           :: fd, rc
    integer                       :: iport
    found      = .false.
    self%port  = -1
    addrlen    = int(storage_size(addr) / 8, c_socklen_t)
    do iport = TCP_PORT_MIN, TCP_PORT_MAX
      fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
      if( fd < 0 ) cycle
      addr%sin_family = int(AF_INET,    c_int16_t)
      addr%sin_port   = c_htons(transfer(int(iport, c_int32_t), 0_c_int16_t))
      addr%sin_addr   = int(INADDR_ANY, c_int32_t)
      rc = c_bind(fd, c_loc(addr), addrlen)
      if( c_close(fd) /= 0 ) continue  ! best-effort close
      if( rc == 0 ) then
        self%port = iport
        found     = .true.
        return
      end if
    end do
  end subroutine find_available_port

  !> Given a comma-separated list of IP addresses, attempt a TCP connection on
  !> the given port to each one and return the first that is reachable (found=.true.).
  subroutine find_available_server( self, server_ips, found )
    class(ipc_tcp_socket), intent(inout) :: self
    class(string),         intent(in)    :: server_ips
    logical,               intent(out)   :: found
    type(tcp_sockaddr_in), target          :: addr
    integer(kind=c_socklen_t)              :: addrlen
    integer(kind=c_int)                    :: fd, rc_connect, rc_close
    character(kind=c_char, len=65), target :: ip_cstr
    character(len=512)                     :: iplist
    character(len=64)                      :: token
    integer                                :: i, tok_start, iplen, tok_len, tlen
    found          = .false.
    self%server_ip = string('')
    if( self%port < 0 ) return
    iplist    = server_ips%to_char()
    iplen     = len_trim(iplist)
    addrlen   = int(storage_size(addr) / 8, c_socklen_t)
    tok_start = 1
    do i = 1, iplen + 1
      if( i > iplen .or. iplist(i:i) == ',' ) then
        tok_len = i - tok_start
        if( tok_len > 0 ) then
          token = adjustl(iplist(tok_start : tok_start + tok_len - 1))
          tlen  = len_trim(token)
          if( tlen > 0 .and. tlen <= 64 ) then
            ip_cstr              = ''
            ip_cstr(1:tlen)      = token(1:tlen)
            ip_cstr(tlen+1:tlen+1) = c_null_char
            fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
            if( fd >= 0 ) then
              addr%sin_family = int(AF_INET, c_int16_t)
              addr%sin_port   = c_htons(transfer(int(self%port, c_int32_t), 0_c_int16_t))
              addr%sin_addr   = int(c_inet_addr(ip_cstr), c_int32_t)
              rc_connect = c_connect(fd, c_loc(addr), addrlen)
              rc_close   = c_close(fd)
              if( rc_connect == 0 ) then
                self%server_ip = string(token(1:tlen))
                found    = .true.
                write(logfhandle,'(A,A,A,I0)')'>>> IPC_TCP_SOCKET available server found: ', &
                  token(1:tlen), ' on port ', self%port
                return
              end if
            end if
          end if
        end if
        tok_start = i + 1
      end if
    end do
  end subroutine find_available_server

  subroutine send_recv_msg( self, buffer, timeout_ms, max_retries, sent, reply, nread )
    class(ipc_tcp_socket),                intent(in)    :: self
    character(len=:), allocatable, target, intent(in)   :: buffer
    integer,                               intent(in)   :: timeout_ms
    integer,                               intent(in)   :: max_retries
    logical,                               intent(out)  :: sent
    character(kind=c_char),        target, intent(inout) :: reply(TCP_MAX_MSG)
    integer,                               intent(out),  optional :: nread
    type(tcp_sockaddr_in),          target :: addr
    type(c_timeval),                target :: tv
    integer(kind=c_socklen_t)              :: addrlen, tv_sz
    integer(kind=c_int)                    :: fd, rc
    integer(kind=c_size_t)                 :: nr
    character(len=65), target              :: ip_cstr
    integer                                :: itry, iplen
    sent = .false.
    if( present(nread) ) nread = 0
    if( self%port < 0 ) return
    if( .not. self%server_ip%is_allocated() ) return
    iplen = self%server_ip%strlen_trim()
    if( iplen == 0 .or. iplen > 64 ) return
    ! build null-terminated IP string
    ip_cstr                  = ''
    ip_cstr(1:iplen)         = self%server_ip%to_char()
    ip_cstr(iplen+1:iplen+1) = c_null_char
    addrlen = int(storage_size(addr) / 8, c_socklen_t)
    tv_sz   = int(storage_size(tv)   / 8, c_socklen_t)
    ! convert timeout_ms -> timeval
    tv%tv_sec  = int(timeout_ms / 1000,            c_long)
    tv%tv_usec = int(mod(timeout_ms, 1000) * 1000, c_long)
    do itry = 1, max_retries
      fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
      if( fd < 0 ) then
        rc = c_usleep(int(timeout_ms * 1000, c_useconds_t))
        cycle
      end if
      ! set per-socket send/receive timeout
      rc = c_setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, c_loc(tv), tv_sz)
      rc = c_setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, c_loc(tv), tv_sz)
      addr%sin_family = int(AF_INET, c_int16_t)
      addr%sin_port   = c_htons(transfer(int(self%port, c_int32_t), 0_c_int16_t))
      addr%sin_addr   = int(c_inet_addr(ip_cstr), c_int32_t)
      rc = c_connect(fd, c_loc(addr), addrlen)
      if( rc == 0 ) then
        rc = int(c_send(fd, c_loc(buffer), int(len(buffer), c_size_t), 0_c_int))
        if( rc >= 0 ) then
          sent = .true.
          nr = c_read(fd, c_loc(reply), int(size(reply), c_size_t))
          if( present(nread) ) nread = int(nr)
          rc = c_close(fd)
          return
        end if
        rc = c_close(fd)
      else
        rc = c_close(fd)
      end if
      if( itry < max_retries ) rc = c_usleep(int(timeout_ms * 1000, c_useconds_t))
    end do
  end subroutine send_recv_msg

  !> Accept one connection on server_fd, read into buf, close the connection.
  !> ok=.false. if accept() failed — caller should exit the accept loop.
  subroutine recv_msg( server_fd, conn_fd, buf, nread, ok, close_after_read )
    integer(kind=c_int),           intent(in)         :: server_fd
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
    conn_fd = c_accept(server_fd, c_null_ptr, int(0, c_socklen_t))
    if( conn_fd < 0 ) return
    nr    = c_read(conn_fd, c_loc(buf), int(len(buf), c_size_t))
    if( l_close_conn ) rc = c_close(conn_fd)
    nread = int(nr)
    ok    = .true.
  end subroutine recv_msg

  subroutine repl_msg( conn_fd, buf, nread, ok )
    integer(kind=c_int),           intent(in)         :: conn_fd
    character(kind=c_char, len=*), intent(inout), target :: buf
    integer,                       intent(out)        :: nread
    logical,                       intent(out)        :: ok
    integer(kind=c_int)    :: rc
    integer(kind=c_size_t) :: nr
    ok    = .false.
    rc    = int(c_send(conn_fd, c_loc(buf), int(len(buf), c_size_t), 0_c_int))
    rc    = c_close(conn_fd)
    nread = int(nr)
    ok    = .true.
  end subroutine repl_msg



end module simple_ipc_tcp_socket