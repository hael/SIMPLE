!@descr: POSIX message-queue wrapper (ipc_mq) for inter-process communication
!==============================================================================
! MODULE: simple_ipc_mq
!
! PURPOSE:
!   Provides a Fortran object-oriented wrapper around POSIX mqueue(7).
!   An ipc_mq instance owns a single named queue identified by the
!   process PID, supporting blocking send, timed send, non-blocking
!   receive, and timed receive. Each overloaded operation accepts
!   either a string object or a raw allocatable character buffer.
!
! TYPES:
!   ipc_mq_attr — mirrors the POSIX mq_attr struct (flags, capacities,
!                 current fill level).
!   ipc_mq      — owns an open queue descriptor; call new() to create,
!                 kill() to close and unlink.
!
! PLATFORM:
!   POSIX message queues are Linux-only. On all other platforms every
!   procedure compiles to a no-op stub so that calling code need not
!   be guarded by #ifdef.
!
! PARAMETERS (hard-coded):
!   MQ_PERM — queue permission bits: owner rw, group r, world r (0644)
!
! DEPENDENCIES:
!   unix, simple_core_module_api
!==============================================================================
module simple_ipc_mq
  use unix,                    only: c_mq_attr, c_mq_open, c_mq_close, c_mq_unlink,  &
                                    c_mq_send, c_mq_timedsend,                       &
                                    c_mq_receive, c_mq_timedreceive,                 &
                                    c_mq_getattr, c_timespec, c_time,                &
                                    O_CREAT, O_RDWR, ENOENT,                         &
                                    c_pid_t, c_mqd_t, c_mode_t,                      &
                                    c_getpid, c_perror
  use, intrinsic :: iso_c_binding, only: c_null_char, c_null_ptr, c_loc, c_long, c_size_t
  use simple_error,           only: simple_exception
  use simple_core_module_api, only: string, int2str, XLONGSTRLEN, logfhandle

  implicit none

  integer, private, parameter :: MQ_PERM = int(o'0644') ! queue permission bits (octal)

  public  :: ipc_mq
  public  :: ipc_mq_attr
  private
#include "simple_local_flags.inc"

  ! Mirrors the POSIX mq_attr struct; populated by update_attributes().
  type :: ipc_mq_attr
    integer(kind=c_long) :: mq_flags   = 0  ! queue flags
    integer(kind=c_long) :: mq_maxmsg  = 0  ! maximum number of messages
    integer(kind=c_long) :: mq_msgsize = 0  ! maximum message size (bytes)
    integer(kind=c_long) :: mq_curmsgs = 0  ! messages currently queued
  end type ipc_mq_attr

  ! Owns a single POSIX message queue. The queue name is derived from the
  ! user-supplied label and the owning process PID so it is unique per process.
  type :: ipc_mq
    private
    type(ipc_mq_attr)     :: attr
    type(string)          :: name
    type(string)          :: mq_path
    integer(kind=c_pid_t) :: mq_ppid
    integer(kind=c_mqd_t) :: mqfd
    integer               :: max_msgsize = 0
    logical               :: active      = .false.
  contains
    procedure :: new
    procedure :: unlink
    procedure :: kill
    procedure :: update_attributes
    procedure :: get_attributes
    procedure :: print_attributes
    procedure :: send_1
    procedure :: send_2
    procedure :: send_timed_1
    procedure :: send_timed_2
    procedure :: receive_1
    procedure :: receive_2
    procedure :: receive_timed_1
    procedure :: receive_timed_2
    procedure :: get_queue_length
    procedure :: is_active
    generic   :: send          => send_1, send_2
    generic   :: send_timed    => send_timed_1, send_timed_2
    generic   :: receive       => receive_1, receive_2
    generic   :: receive_timed => receive_timed_1, receive_timed_2
  end type ipc_mq

contains

#if defined(__linux__)

  ! Create and open a new POSIX message queue. The queue path is
  ! /<name>.<pid>. If max_msgsize is given the queue is configured with
  ! that message size limit and a capacity of 10 messages; otherwise the
  ! system default attributes are used. Any pre-existing queue with the
  ! same path is unlinked first.
  subroutine new( self, name, max_msgsize )
    class(ipc_mq),          intent(inout) :: self
    type(string), optional, intent(in)    :: name
    integer,      optional, intent(in)    :: max_msgsize
    type(c_mq_attr), target               :: attr
    if( present(name)        ) self%name        = name
    if( present(max_msgsize) ) self%max_msgsize = max_msgsize
    self%mq_ppid = c_getpid()
    self%mq_path = '/' // self%name%to_char() // '.' // int2str(self%mq_ppid)
    call self%unlink()
    if( self%max_msgsize > 0 ) then
      attr%mq_msgsize = self%max_msgsize
      attr%mq_maxmsg  = 10
      self%mqfd = c_mq_open(name  = self%mq_path%to_char() // c_null_char, &
                            oflag = ior(O_CREAT, O_RDWR),                  &
                            mode  = int(MQ_PERM, kind=c_mode_t),           &
                            attr  = c_loc(attr))
    else
      self%mqfd = c_mq_open(name  = self%mq_path%to_char() // c_null_char, &
                            oflag = ior(O_CREAT, O_RDWR),                  &
                            mode  = int(MQ_PERM, kind=c_mode_t),           &
                            attr  = c_null_ptr)
    end if
    if( self%mqfd < 0 ) then
      call c_perror('mq_open()' // c_null_char)
      THROW_HARD('failed to create message queue')
    end if
    self%active = .true.
    call self%update_attributes()
  end subroutine new

  ! Unlink (delete) the queue from the system. ENOENT is silently ignored
  ! since it means the queue was already removed.
  subroutine unlink( self )
    class(ipc_mq), intent(inout) :: self
    integer                      :: stat, err
    if( .not. self%active ) return
    stat = c_mq_unlink(self%mq_path%to_char() // c_null_char)
    if( stat < 0 ) then
      err = ierrno()
      if( err /= ENOENT ) then
        call c_perror('mq_unlink()' // c_null_char)
        THROW_HARD('failed to unlink message queue')
      end if
    end if
  end subroutine unlink

  ! Close the queue descriptor and unlink the queue, then mark inactive.
  subroutine kill( self )
    class(ipc_mq), intent(inout) :: self
    integer                      :: stat
    if( .not. self%active ) return
    stat = c_mq_close(self%mqfd)
    if( stat < 0 ) then
      call c_perror('mq_close()' // c_null_char)
      THROW_HARD('failed to close message queue')
    end if
    call self%unlink()
    self%active = .false.
  end subroutine kill

  ! Refresh self%attr from the kernel.
  subroutine update_attributes( self )
    class(ipc_mq), intent(inout) :: self
    type(c_mq_attr)              :: mq_attr
    integer                      :: stat
    stat = c_mq_getattr(self%mqfd, mq_attr)
    if( stat < 0 ) then
      call c_perror('mq_getattr()' // c_null_char)
      THROW_HARD('failed to get message queue attributes')
    end if
    self%attr%mq_flags   = mq_attr%mq_flags
    self%attr%mq_maxmsg  = mq_attr%mq_maxmsg
    self%attr%mq_msgsize = mq_attr%mq_msgsize
    self%attr%mq_curmsgs = mq_attr%mq_curmsgs
  end subroutine update_attributes

  ! Send a string message (blocking, priority 1).
  subroutine send_1( self, msg )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(in)    :: msg
    integer                      :: stat
    if( .not. self%active ) return
    stat = c_mq_send(self%mqfd, msg%to_char(), int(msg%strlen_trim(), kind=c_size_t), 1)
    if( stat < 0 ) then
      call c_perror('mq_send()' // c_null_char)
      THROW_HARD('failed to send message queue message')
    end if
    call self%update_attributes()
  end subroutine send_1

  ! Send a raw character buffer (blocking, priority 1).
  subroutine send_2( self, buffer )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer                                      :: stat
    if( .not. self%active ) return
    stat = c_mq_send(self%mqfd, buffer, int(len(buffer), kind=c_size_t), 1)
    if( stat < 0 ) then
      call c_perror('mq_send()' // c_null_char)
      THROW_HARD('failed to send message queue message')
    end if
    call self%update_attributes()
  end subroutine send_2

  ! Timed send of a string; waits at most timeout_s seconds for queue space.
  ! Returns .true. if the message was enqueued before the deadline.
  function send_timed_1( self, msg, timeout_s ) result( sent )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(in)    :: msg
    integer,       intent(in)    :: timeout_s
    type(c_timespec), target     :: timeout
    logical                      :: sent
    integer                      :: stat
    sent = .false.
    if( .not. self%active ) return
    timeout%tv_nsec = 0
    timeout%tv_sec  = int(c_time(0_c_long)) + timeout_s
    stat = c_mq_timedsend(self%mqfd, msg%to_char(), int(msg%strlen_trim(), kind=c_size_t), 1, timeout)
    if( stat < 0 ) then
      call c_perror('mq_timedsend()' // c_null_char)
      THROW_WARN('failed to timed-send message queue message')
    else
      sent = .true.
    end if
    call self%update_attributes()
  end function send_timed_1

  ! Timed send of a raw character buffer; waits at most timeout_s seconds
  ! for queue space. Returns .true. if the message was enqueued before the deadline.
  function send_timed_2( self, buffer, timeout_s ) result( sent )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer,                       intent(in)    :: timeout_s
    type(c_timespec), target                     :: timeout
    logical                                      :: sent
    integer                                      :: stat
    sent = .false.
    if( .not. self%active ) return
    timeout%tv_nsec = 0
    timeout%tv_sec  = int(c_time(0_c_long)) + timeout_s
    stat = c_mq_timedsend(self%mqfd, buffer, int(len(buffer), kind=c_size_t), 1, timeout)
    if( stat < 0 ) then
      call c_perror('mq_timedsend()' // c_null_char)
      THROW_WARN('failed to timed-send message queue message')
    else
      sent = .true.
    end if
    call self%update_attributes()
  end function send_timed_2

  ! Non-blocking receive into a string. Returns .true. if a message was read.
  function receive_1( self, msg ) result( received )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(inout) :: msg
    character(len=XLONGSTRLEN)   :: buf
    logical                      :: received
    integer                      :: sz, prio
    received = .false.
    if( .not. self%active ) return
    call msg%kill()
    call self%update_attributes()
    if( self%attr%mq_curmsgs > 0 ) then
      sz = c_mq_receive(self%mqfd, buf, len(buf, kind=c_size_t), prio)
      if( sz < 0 ) then
        call c_perror('mq_receive()' // c_null_char)
        THROW_HARD('failed to receive message queue message')
      else if( sz > 0 ) then
        msg      = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_1

  ! Non-blocking receive into an allocatable character buffer.
  ! Returns .true. if a message was read.
  function receive_2( self, buffer ) result( received )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    character(len=XLONGSTRLEN)                   :: buf
    logical                                      :: received
    integer                                      :: sz, prio
    received = .false.
    if( .not. self%active ) return
    if( allocated(buffer) ) deallocate(buffer)
    call self%update_attributes()
    if( self%attr%mq_curmsgs > 0 ) then
      sz = c_mq_receive(self%mqfd, buf, len(buf, kind=c_size_t), prio)
      if( sz < 0 ) then
        call c_perror('mq_receive()' // c_null_char)
        THROW_HARD('failed to receive message queue message')
      else if( sz > 0 ) then
        allocate(character(len=sz) :: buffer)
        buffer   = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_2

  ! Timed receive into a string; waits at most timeout_s seconds.
  ! Returns .true. if a message was read before the deadline.
  function receive_timed_1( self, msg, timeout_s ) result( received )
    class(ipc_mq),    intent(inout) :: self
    class(string),    intent(inout) :: msg
    integer,          intent(in)    :: timeout_s
    type(c_timespec), target        :: timeout
    character(len=XLONGSTRLEN)      :: buf
    logical                         :: received
    integer                         :: sz, prio
    received = .false.
    if( .not. self%active ) return
    call msg%kill()
    call self%update_attributes()
    timeout%tv_nsec = 0
    timeout%tv_sec  = int(c_time(0_c_long)) + timeout_s
    if( self%attr%mq_curmsgs > 0 ) then
      sz = c_mq_timedreceive(self%mqfd, buf, len(buf, kind=c_size_t), prio, timeout)
      if( sz < 0 ) then
        call c_perror('mq_timedreceive()' // c_null_char)
        THROW_WARN('failed to receive timed message queue message')
      else if( sz > 0 ) then
        msg      = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_timed_1

  ! Timed receive into an allocatable character buffer; waits at most
  ! timeout_s seconds. Returns .true. if a message was read before the deadline.
  function receive_timed_2( self, buffer, timeout_s ) result( received )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer,                       intent(in)    :: timeout_s
    type(c_timespec),              target        :: timeout
    character(len=XLONGSTRLEN)                   :: buf
    logical                                      :: received
    integer                                      :: sz, prio
    received = .false.
    if( .not. self%active ) return
    if( allocated(buffer) ) deallocate(buffer)
    timeout%tv_nsec = 0
    timeout%tv_sec  = int(c_time(0_c_long)) + timeout_s
    call self%update_attributes()
    if( self%attr%mq_curmsgs > 0 ) then
      sz = c_mq_timedreceive(self%mqfd, buf, len(buf, kind=c_size_t), prio, timeout)
      if( sz < 0 ) then
        call c_perror('mq_timedreceive()' // c_null_char)
        THROW_WARN('failed to receive timed message queue message')
      else if( sz > 0 ) then
        allocate(character(len=sz) :: buffer)
        buffer   = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_timed_2

#else
  ! Non-Linux stubs — POSIX message queues are unavailable on this platform.

  subroutine new( self, name, max_msgsize )
    class(ipc_mq),          intent(inout) :: self
    type(string), optional, intent(in)    :: name
    integer,      optional, intent(in)    :: max_msgsize
    self%active = .true.
    THROW_WARN('message queues are not available on this system')
  end subroutine new

  subroutine unlink( self )
    class(ipc_mq), intent(inout) :: self
  end subroutine unlink

  subroutine kill( self )
    class(ipc_mq), intent(inout) :: self
  end subroutine kill

  subroutine update_attributes( self )
    class(ipc_mq), intent(inout) :: self
  end subroutine update_attributes

  subroutine send_1( self, msg )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(in)    :: msg
  end subroutine send_1

  subroutine send_2( self, buffer )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
  end subroutine send_2

  function send_timed_1( self, msg, timeout_s ) result( sent )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(in)    :: msg
    integer,       intent(in)    :: timeout_s
    logical                      :: sent
    sent = .false.
  end function send_timed_1

  function send_timed_2( self, buffer, timeout_s ) result( sent )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer,                       intent(in)    :: timeout_s
    logical                                      :: sent
    sent = .false.
  end function send_timed_2

  function receive_1( self, msg ) result( received )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(inout) :: msg
    logical                      :: received
    received = .false.
    call msg%kill()
  end function receive_1

  function receive_2( self, buffer ) result( received )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    logical                                      :: received
    received = .false.
    if( allocated(buffer) ) deallocate(buffer)
  end function receive_2

  function receive_timed_1( self, msg, timeout_s ) result( received )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(inout) :: msg
    integer,       intent(in)    :: timeout_s
    logical                      :: received
    received = .false.
    call msg%kill()
  end function receive_timed_1

  function receive_timed_2( self, buffer, timeout_s ) result( received )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer,                       intent(in)    :: timeout_s
    logical                                      :: received
    received = .false.
    if( allocated(buffer) ) deallocate(buffer)
  end function receive_timed_2

#endif

  ! Copy the cached attributes into attr.
  subroutine get_attributes( self, attr )
    class(ipc_mq),     intent(in)  :: self
    type(ipc_mq_attr), intent(out) :: attr
    attr = self%attr
  end subroutine get_attributes

  ! Write queue attributes to logfhandle.
  subroutine print_attributes( self )
    class(ipc_mq), intent(in) :: self
    if( .not. self%active ) return
    write(logfhandle, '(A,I0)') 'MQ Flags            :', self%attr%mq_flags
    write(logfhandle, '(A,I0)') 'MQ Max Messages     :', self%attr%mq_maxmsg
    write(logfhandle, '(A,I0)') 'MQ Max Message Size :', self%attr%mq_msgsize
    write(logfhandle, '(A,I0)') 'MQ Current Messages :', self%attr%mq_curmsgs
  end subroutine print_attributes

  ! Return the number of messages currently in the queue.
  function get_queue_length( self ) result( queue_length )
    class(ipc_mq),       intent(in) :: self
    integer(kind=c_long)            :: queue_length
    queue_length = self%attr%mq_curmsgs
  end function get_queue_length

  ! Return .true. if the queue is open and ready.
  function is_active( self ) result( active )
    class(ipc_mq), intent(in) :: self
    logical                   :: active
    active = self%active
  end function is_active

end module simple_ipc_mq
