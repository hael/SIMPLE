!@descr: low-level POSIX socket helpers shared by the IPC TCP client and server (poll, accept, liveness)
module simple_ipc_tcp_socket_helpers
  use iso_c_binding
  use unix,                   only: c_accept, c_null_ptr, c_socklen_t
  implicit none

  integer, parameter, public :: TCP_MAX_MSG  = 1460  ! 1500 MTU - 20 IP header - 20 TCP header

  ! POSIX poll event flags
  integer(kind=c_short), parameter, public :: POLLIN   =  1_c_short
  integer(kind=c_short), parameter, public :: POLLOUT  =  4_c_short
  integer(kind=c_short), parameter, public :: POLLERR  =  8_c_short  ! Linux POLLERR
  integer(kind=c_short), parameter, public :: POLLHUP  = 16_c_short  ! Linux POLLHUP
  integer(kind=c_short), parameter, public :: POLLNVAL = 32_c_short  ! Linux POLLNVAL (invalid fd)

  public :: c_pollfd
  public :: fd_is_healthy
  public :: poll_fds
  public :: accept_connection

  ! POSIX struct pollfd — one entry per fd passed to poll(2).
  type, bind(c) :: c_pollfd
    integer(kind=c_int)   :: fd      = -1
    integer(kind=c_short) :: events  = 0_c_short   ! events to watch (e.g. POLLIN)
    integer(kind=c_short) :: revents = 0_c_short   ! events returned by poll()
  end type c_pollfd

  private

  integer(kind=c_int), parameter :: MSG_PEEK     =  2_c_int  ! POSIX MSG_PEEK
  integer(kind=c_int), parameter :: MSG_DONTWAIT = 64_c_int  ! Linux MSG_DONTWAIT

  ! Direct C bindings for POSIX poll(2) and recv(2).
  interface

    function c_poll(fds, nfds, timeout) result(rc) bind(c, name='poll')
      import :: c_int, c_long, c_pollfd
      type(c_pollfd),       intent(inout) :: fds(*)
      integer(kind=c_long), value         :: nfds     ! nfds_t = unsigned long on Linux
      integer(kind=c_int),  value         :: timeout  ! milliseconds; -1 = block, 0 = nowait
      integer(kind=c_int)                 :: rc
    end function c_poll

    function c_recv(sockfd, buf, len, flags) result(nr) bind(c, name='recv')
      import :: c_int, c_ptr, c_size_t
      integer(kind=c_int),    value :: sockfd
      type(c_ptr),            value :: buf
      integer(kind=c_size_t), value :: len
      integer(kind=c_int),    value :: flags
      integer(kind=c_size_t)        :: nr   ! ssize_t; EOF=0, error wraps to large positive
    end function c_recv

  end interface

  contains

  !> Non-blocking liveness check for an open socket fd.
  !>
  !> Returns .true. only when:
  !>   - poll(0) reports no ERR/HUP/NVAL state, and
  !>   - if readable, recv(MSG_PEEK|MSG_DONTWAIT) confirms data (not EOF/error).
  !>
  !> This routine never consumes payload bytes from the socket.
  logical function fd_is_healthy( fd )
    integer(kind=c_int), intent(in) :: fd
    type(c_pollfd)                   :: pfd(1)
    character(kind=c_char), target   :: peek_buf(1)
    integer(kind=c_int)              :: rc
    integer(kind=c_size_t)           :: nr
    fd_is_healthy = .false.
    if( fd < 0 ) return
    pfd(1)%fd      = fd
    pfd(1)%events  = POLLIN
    pfd(1)%revents = 0_c_short
    rc = c_poll(pfd, 1_c_long, 0_c_int)
    if( rc < 0 )                                       return   ! poll error
    if( iand(pfd(1)%revents, POLLNVAL) /= 0_c_short )  return   ! fd not open
    if( iand(pfd(1)%revents, POLLERR)  /= 0_c_short )  return   ! socket error
    if( iand(pfd(1)%revents, POLLHUP)  /= 0_c_short )  return   ! peer closed (RST/FIN)
    if( iand(pfd(1)%revents, POLLIN)   /= 0_c_short ) then
      ! Readable: could be data or EOF — peek one byte to distinguish.
      ! MSG_DONTWAIT ensures we never block; EOF gives nr==0.
      ! Errors from c_recv (ssize_t -1) wrap into a large c_size_t value.
      nr = c_recv(fd, c_loc(peek_buf), 1_c_size_t, ior(MSG_PEEK, MSG_DONTWAIT))
      if( nr == 0_c_size_t .or. nr > 1_c_size_t ) return  ! EOF or recv error
    end if
    fd_is_healthy = .true.
  end function fd_is_healthy

  !> poll(2) fds(1:n) (caller sets %fd) for POLLIN; nready = number ready, < 0 on poll error or
  !> n outside 0..size(fds). timeout_ms: -1 blocks, 0 returns at once. Callers then test
  !> iand(fds(i)%revents, POLLIN).
  subroutine poll_fds( fds, n, timeout_ms, nready )
    type(c_pollfd), intent(inout) :: fds(:)
    integer,        intent(in)    :: n
    integer,        intent(in)    :: timeout_ms
    integer,        intent(out)   :: nready
    integer(kind=c_int) :: rc
    integer             :: i
    if( n < 0 .or. n > size(fds) ) then
      nready = -1
      return
    end if
    if( n == 0 ) then
      nready = 0
      return
    end if
    do i = 1, n
      fds(i)%events  = POLLIN
      fds(i)%revents = 0_c_short
    end do
    rc     = c_poll(fds, int(n, c_long), int(timeout_ms, c_int))
    nready = int(rc)
  end subroutine poll_fds

  !> Accept one incoming connection on \\p fd.
  !> Returns the accepted connection fd, or < 0 on failure.
  function accept_connection( fd ) result( conn_fd )
    integer(kind=c_int), intent(in) :: fd
    integer(kind=c_int)             :: conn_fd
    conn_fd = c_accept(fd, c_null_ptr, int(0, c_socklen_t))
  end function accept_connection

end module simple_ipc_tcp_socket_helpers
