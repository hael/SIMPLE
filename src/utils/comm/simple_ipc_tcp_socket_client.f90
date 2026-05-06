!==============================================================================
! MODULE: simple_ipc_tcp_socket_client
!
! PURPOSE:
!   Lightweight TCP client used by SIMPLE components to send a request buffer
!   and read a reply from any reachable server IP in a comma-separated list.
!
! DESIGN NOTES:
!   - Reuses an existing healthy socket when possible.
!   - Applies SO_KEEPALIVE and send/receive timeouts on new sockets.
!   - Handles partial sends by looping until the full payload is written.
!   - Handles fragmented replies by accumulating reads into rcv_buffer.
!==============================================================================
module simple_ipc_tcp_socket_client
  use iso_c_binding
  use unix,          only: c_socket, c_connect, c_send, c_setsockopt, c_close, c_read, &
                            c_usleep, c_useconds_t,                            &
                            c_htons, c_inet_addr,                              &
                            AF_INET, SOCK_STREAM,                              &
                            SOL_SOCKET, SO_KEEPALIVE, c_socklen_t
  use simple_core_module_api, only: logfhandle, string
  use simple_ipc_tcp_socket_helpers, only: c_pollfd, fd_is_healthy, TCP_MAX_MSG

  implicit none

  integer, parameter, public  :: TCP_MAX_RETRIES = 5
  integer, parameter, public  :: TCP_TIMEOUT_MS  = 1000    ! 1 second timeout for connect/send/recv operations
  logical, parameter, private :: DEBUG           = .false. ! Debug flag

  public :: ipc_tcp_socket_client
  public :: tcp_sockaddr_in
  public :: c_pollfd

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

  integer(kind=c_int),   parameter :: SO_SNDTIMEO    = 21    ! Linux SO_SNDTIMEO
  integer(kind=c_int),   parameter :: SO_RCVTIMEO    = 20    ! Linux SO_RCVTIMEO

  ! POSIX struct timeval used with SO_SNDTIMEO via setsockopt.
  type, bind(c) :: c_timeval
    integer(kind=c_long) :: tv_sec  = 0_c_long
    integer(kind=c_long) :: tv_usec = 0_c_long
  end type c_timeval

  type :: ipc_tcp_socket_client
    private
    integer               :: port       = -1
    integer(kind=c_int)   :: fd         = -1
    type(string)          :: server_ips
    type(string)          :: server_ip
  contains
    procedure :: new
    procedure :: connect
    procedure :: kill
    procedure :: send_recv_msg
  end type ipc_tcp_socket_client

  contains

  !> Reset state and set the target server list and port.
  subroutine new(self, server_list, port)
    class(ipc_tcp_socket_client), intent(inout) :: self
    type(string),          intent(in)    :: server_list
    integer,               intent(in)    :: port
    call self%kill()  ! close any existing connection and reset state
    self%port       = port
    self%server_ips = server_list
  end subroutine new

  !> Close an open socket (if any) and reset object fields.
  subroutine kill(self)
    class(ipc_tcp_socket_client), intent(inout) :: self
    integer(kind=c_int)                         :: rc
    if( self%fd >= 0 ) then
      rc = c_close(self%fd)
      if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: c_close failed in kill(), rc=', rc
    end if
    self%fd         = -1
    self%port       = -1
    self%server_ip  = string('')
    self%server_ips = string('')
  end subroutine kill

  !> Try each IP in self%server_ips on self%port until one connects.
  !> succeeds.  Reuses self%fd if it is still healthy.  SO_KEEPALIVE keeps dead connections from hanging send/recv calls indefinitely, but
  !> per-socket send/receive timeouts are always set on new connections.
  !> self%fd remains -1 if every IP is unreachable.
  subroutine connect( self )
    class(ipc_tcp_socket_client), intent(inout) :: self
    type(tcp_sockaddr_in),               target :: addr
    character(kind=c_char, len=65),      target :: ip_cstr
    integer(kind=c_int),                 target :: optval
    type(c_timeval),                     target :: tv
    integer(kind=c_socklen_t)            :: addrlen, tv_sz
    integer(kind=c_int)                  :: fd, rc
    character(len=512)                   :: iplist
    character(len=64)                    :: token
    integer                              :: i, tok_start, iplen, tok_len, tlen
    if( self%port < 0 )                        return
    if( .not. self%server_ips%is_allocated() ) return
    ! reuse the existing connection if it is still open and healthy
    if( self%fd >= 0 ) then
      if( fd_is_healthy(self%fd) ) then
        if( DEBUG ) write(logfhandle,'(A,A,A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: reusing healthy existing connection to server: ', &
          self%server_ip%to_char(), ' on port ', self%port
        return
      end if
      if( DEBUG ) write(logfhandle,'(A)') '>>> IPC_TCP_SOCKET_CLIENT: existing connection is not healthy, closing and retrying'
      rc = c_close(self%fd)
      self%fd = -1
    end if
    self%server_ip = string('')
    iplist    = self%server_ips%to_char()
    iplen     = len_trim(iplist)
    addrlen   = int(storage_size(addr) / 8, c_socklen_t)
    tv_sz     = int(storage_size(tv)   / 8, c_socklen_t)
    tok_start = 1
    ! convert timeout_ms -> timeval
    tv%tv_sec  = int(TCP_TIMEOUT_MS / 1000,            c_long)
    tv%tv_usec = int(mod(TCP_TIMEOUT_MS, 1000) * 1000, c_long)
    do i = 1, iplen + 1
      if( i > iplen .or. iplist(i:i) == ',' ) then
        tok_len = i - tok_start
        if( tok_len > 0 ) then
          token = adjustl(iplist(tok_start : tok_start + tok_len - 1))
          tlen  = len_trim(token)
          if( tlen > 0 .and. tlen <= 64 ) then
            ! build null-terminated IP string for C calls
            ip_cstr              = ''
            ip_cstr(1:tlen)      = token(1:tlen)
            ip_cstr(tlen+1:tlen+1) = c_null_char
            fd = c_socket(AF_INET, SOCK_STREAM, 0_c_int)
            if( fd >= 0 ) then
              optval = 1_c_int
              rc = c_setsockopt(fd, SOL_SOCKET, SO_KEEPALIVE, c_loc(optval), int(storage_size(optval)/8, c_socklen_t))
              if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: c_setsockopt(SO_KEEPALIVE) failed, rc=', rc
              rc = c_setsockopt(fd, SOL_SOCKET, SO_SNDTIMEO, c_loc(tv), tv_sz)
              if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: c_setsockopt(SO_SNDTIMEO) failed, rc=', rc
              rc = c_setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, c_loc(tv), tv_sz)
              if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: c_setsockopt(SO_RCVTIMEO) failed, rc=', rc
              addr%sin_family = int(AF_INET, c_int16_t)
              addr%sin_port   = c_htons(transfer(int(self%port, c_int32_t), 0_c_int16_t))
              addr%sin_addr   = int(c_inet_addr(ip_cstr), c_int32_t)
              rc = c_connect(fd, c_loc(addr), addrlen)
              if( rc == 0 ) then
                self%fd = fd
                self%server_ip = string(token(1:tlen))
                if( DEBUG ) write(logfhandle,'(A,A,A,I0)')'>>> IPC_TCP_SOCKET connected to server: ', &
                  self%server_ip%to_char(), ' on port ', self%port
                return
              end if
              rc = c_close(fd)
              if( rc /= 0 ) write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: c_close failed after failed connect, rc=', rc
            end if
          end if
        end if
        tok_start = i + 1
      end if
    end do
    write(logfhandle,'(A,I0)') '>>> IPC_TCP_SOCKET_CLIENT: failed to connect to any server IP on port ', self%port
  end subroutine connect

  ! !> Send a request and read a reply with retry/reconnect semantics.
  ! !! sent=.true. only when at least one reply byte has been read.
  ! subroutine send_recv_msg( self, snd_buffer, rcv_buffer, sent, nread )
  !   class(ipc_tcp_socket_client),          intent(inout) :: self
  !   character(kind=c_char),        target, intent(inout) :: rcv_buffer(TCP_MAX_MSG)
  !   character(len=:), allocatable, target, intent(in)    :: snd_buffer
  !   logical,                               intent(out)   :: sent
  !   integer,                               intent(out)   :: nread
  !   integer(kind=c_int)                         :: rc
  !   integer(kind=c_size_t)                      :: nr
  !   character(kind=c_char), allocatable, target :: snd_bytes(:)
  !   integer                                     :: itry, i, nsend, sent_total, nrecv_total, remaining
  !   sent  = .false.
  !   nread = 0
  !   if( len(snd_buffer) < 1 ) return

  !   ! Convert Fortran CHARACTER to contiguous c_char bytes for pointer-safe
  !   ! partial-send retries via c_loc(snd_bytes(k)).
  !   allocate(snd_bytes(len(snd_buffer)))
  !   do i = 1, len(snd_buffer)
  !     snd_bytes(i) = snd_buffer(i:i)
  !   end do

  !   do itry = 1, TCP_MAX_RETRIES
  !     call self%connect()  ! ensure we have a connection before trying to send
  !     if( self%fd >= 0 ) then
  !       ! Keep sending until the full request payload is written.
  !       sent_total = 0
  !       do while( sent_total < size(snd_bytes) )
  !         remaining = size(snd_bytes) - sent_total
  !         rc = int(c_send(self%fd, c_loc(snd_bytes(sent_total + 1)), int(remaining, c_size_t), 0_c_int))
  !         if( rc <= 0 ) exit
  !         nsend      = rc
  !         sent_total = sent_total + nsend
  !       end do

  !       if( sent_total == size(snd_bytes) ) then
  !         ! Read response bytes. Keep appending until EOF/timeout/error after
  !         ! at least one byte, or until the receive buffer is full.
  !         nrecv_total = 0
  !         do while( nrecv_total < size(rcv_buffer) )
  !           remaining = size(rcv_buffer) - nrecv_total
  !           nr = c_read(self%fd, c_loc(rcv_buffer(nrecv_total + 1)), int(remaining, c_size_t))
  !           if( nr > int(TCP_MAX_MSG, c_size_t) ) then  ! c_read error (ssize_t -1 wraps to HUGE)
  !             if( nrecv_total > 0 ) exit
  !             nrecv_total = 0
  !             exit
  !           end if
  !           if( nr == 0_c_size_t ) exit  ! EOF: done if we already got payload bytes
  !           nrecv_total = nrecv_total + int(nr)
  !         end do

  !         if( nrecv_total > 0 ) then
  !           sent  = .true.
  !           nread = nrecv_total
  !           return
  !         end if
  !       end if

  !       ! send or recv failed — discard the broken socket so connect() reconnects
  !       rc = c_close(self%fd)
  !       self%fd = -1
  !     end if
  !     if( itry < TCP_MAX_RETRIES ) rc = c_usleep(int(int(TCP_TIMEOUT_MS, c_long) * 1000_c_long, c_useconds_t))
  !   end do
  !   write(logfhandle,'(A)') '>>> IPC_TCP_SOCKET_CLIENT: send_recv_msg failed after retries'
  ! end subroutine send_recv_msg

  subroutine send_recv_msg( self, snd_buffer, rcv_buffer, sent, nread )
    class(ipc_tcp_socket_client),          intent(inout) :: self
    character(kind=c_char),        target, intent(inout) :: rcv_buffer(TCP_MAX_MSG)
    character(len=:), allocatable, target, intent(in)    :: snd_buffer
    logical,                               intent(out)   :: sent
    integer,                               intent(out)   :: nread
    integer(kind=c_int)                                  :: rc
    integer(kind=c_size_t)                               :: nr
    integer                                              :: itry
    sent  = .false.
    nread = 0
    do itry = 1, TCP_MAX_RETRIES
      call self%connect()  ! ensure we have a connection before trying to send
      if( self%fd >= 0 ) then
        rc = int(c_send(self%fd, c_loc(snd_buffer), int(len(snd_buffer), c_size_t), 0_c_int))
        if( rc >= 0 ) then
          nr = c_read(self%fd, c_loc(rcv_buffer), int(size(rcv_buffer), c_size_t))
          if( nr <= int(TCP_MAX_MSG, c_size_t) ) then  ! error wraps to HUGE(c_size_t)
            sent  = .true.
            nread = int(nr)
            return
          end if
        end if
        ! send or recv failed — discard the broken socket so connect() reconnects
        rc = c_close(self%fd)
        self%fd = -1
      end if
      if( itry < TCP_MAX_RETRIES ) rc = c_usleep(int(int(TCP_TIMEOUT_MS, c_long) * 1000_c_long, c_useconds_t))
    end do
    write(logfhandle,'(A)') '>>> IPC_TCP_SOCKET_CLIENT: send_recv_msg failed after retries'
  end subroutine send_recv_msg

end module simple_ipc_tcp_socket_client