!@descr: various unix message queue utilities
module simple_ipc_mq
use unix
use simple_core_module_api
implicit none

integer, private, parameter :: MQ_PERM = int(o'0644') ! MQ permissions (octal).

public :: ipc_mq
public :: ipc_mq_attr

private
#include "simple_local_flags.inc"

type :: ipc_mq_attr
  integer(kind=c_long) :: mq_flags   = 0
  integer(kind=c_long) :: mq_maxmsg  = 0
  integer(kind=c_long) :: mq_msgsize = 0
  integer(kind=c_long) :: mq_curmsgs = 0
end type ipc_mq_attr

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
  procedure :: receive_1
  procedure :: receive_2
  procedure :: get_queue_length
  procedure :: is_active
  generic   :: send    => send_1, send_2
  generic   :: receive => receive_1, receive_2
end type ipc_mq

contains


#ifdef __linux__
! message queues only work on linux 

  subroutine new( self, name, max_msgsize )
    class(ipc_mq),             intent(inout) :: self 
    type(string),    optional, intent(in)    :: name
    integer,         optional, intent(in)    :: max_msgsize
    type(c_mq_attr), target                  :: attr
    if( present( name )        ) self%name        = name
    if( present( max_msgsize ) ) self%max_msgsize = max_msgsize
    self%mq_ppid = c_getpid()
    self%mq_path = '/' // self%name%to_char() // '.' // int2str( self%mq_ppid )
    ! Unlink, if MQ already exists.
    call self%unlink()
    if( self%max_msgsize > 0 ) then
      attr%mq_msgsize = self%max_msgsize
      attr%mq_maxmsg  = 10
      self%mqfd = c_mq_open(name  = self%mq_path%to_char() // c_null_char, &
                            oflag = ior( O_CREAT, O_RDWR ), &
                            mode  = int( MQ_PERM, kind=c_mode_t ), &
                            attr  = c_loc(attr) &
                           )
    else
      self%mqfd = c_mq_open(name  = self%mq_path%to_char() // c_null_char, &
                            oflag = ior( O_CREAT, O_RDWR ), &
                            mode  = int( MQ_PERM, kind=c_mode_t ), &
                            attr  = c_null_ptr &
                           )
    end if
    if( self%mqfd < 0 ) then
      call c_perror('mq_open()' // c_null_char )
      THROW_HARD('failed to create message queue')
    end if
    self%active = .true.
    call self%update_attributes()
  end subroutine new

  subroutine unlink( self )
    class(ipc_mq), intent(inout) :: self
    integer                      :: stat, err
    if( .not.self%active ) return
    ! message queues only work on linux
    stat = c_mq_unlink(self%mq_path%to_char() // c_null_char)
    if( stat < 0 ) then
      err = ierrno()
      if( err /= ENOENT ) then
        call c_perror('mq_unlink()' // c_null_char)
        THROW_HARD('failed to unlink message queue')
      end if
      return
    end if
  end subroutine unlink

  subroutine kill( self )
    class(ipc_mq), intent(inout) :: self
    integer                      :: stat
    if( .not. self%active ) return
    ! message queues only work on linux
    stat = c_mq_close(self%mqfd)
    if( stat < 0 ) then
      call c_perror('mq_close()' // c_null_char)
      THROW_HARD('failed to close message queue')
    endif
    call self%unlink()
    self%active = .false.
  end subroutine kill

  subroutine update_attributes( self )
    class(ipc_mq), intent(inout) :: self
    type(ipc_mq_attr)            :: attr
    integer                      :: stat 
    type(c_mq_attr)              :: mq_attr  
    stat = c_mq_getattr(self%mqfd, mq_attr)
    if( stat < 0 ) then
      call c_perror('mq_getattr()' // c_null_char)
      THROW_HARD('failed to get message queue attributes')
    endif
    attr%mq_flags   = mq_attr%mq_flags
    attr%mq_maxmsg  = mq_attr%mq_maxmsg
    attr%mq_msgsize = mq_attr%mq_msgsize
    attr%mq_curmsgs = mq_attr%mq_curmsgs
    self%attr = attr
  end subroutine update_attributes

  subroutine send_1( self, msg )
    class(ipc_mq), intent(inout) :: self
    class(string), intent( in )  :: msg
    integer                      :: stat
    if( .not. self%active ) return
    ! message queues only work on linux
    stat = c_mq_send(self%mqfd, msg%to_char(), int(msg%strlen_trim(), kind=c_size_t), 1)
    if(stat < 0) then
        call c_perror('mq_send()' // c_null_char)
        THROW_HARD('failed to send message queue message')
    end if
    call self%update_attributes()
  end subroutine send_1

  subroutine send_2( self, buffer )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer                                      :: stat
    if( .not. self%active ) return
    ! message queues only work on linux
    stat = c_mq_send(self%mqfd, buffer, int(len(buffer), kind=c_size_t), 1)
    if( stat < 0 ) then
        call c_perror('mq_send()' // c_null_char)
        THROW_HARD('failed to send message queue message')
    end if
    call self%update_attributes()
  end subroutine send_2

function receive_1( self, msg ) result( received )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(inout) :: msg
    character( len=XLONGSTRLEN ) :: buf 
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
        msg = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_1

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
      if(sz < 0) then
          call c_perror('mq_receive()' // c_null_char)
          THROW_HARD( 'failed to receive message queue message')
      else if(sz > 0) then
        allocate(character(len=sz) :: buffer)
        buffer = trim(buf(:sz))
        received = .true.
      end if
    end if
    call self%update_attributes()
  end function receive_2

#else
! always inactive on mac -> no message queues

  subroutine new( self, name, max_msgsize )
    class(ipc_mq),             intent(inout) :: self 
    type(string),    optional, intent(in)    :: name
    integer,         optional, intent(in)    :: max_msgsize
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
    class(string), intent( in )  :: msg
    integer                      :: stat
  end subroutine send_1

  subroutine send_2( self, buffer )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    integer                                      :: stat
  end subroutine send_2

  function receive_1( self, msg ) result( received )
    class(ipc_mq), intent(inout) :: self
    class(string), intent(inout) :: msg
    character( len=XLONGSTRLEN ) :: buf 
    logical                      :: received
    integer                      :: sz, prio
    received = .false.
    call msg%kill()
  end function receive_1

  function receive_2( self, buffer ) result( received )
    class(ipc_mq),                 intent(inout) :: self
    character(len=:), allocatable, intent(inout) :: buffer
    character(len=XLONGSTRLEN)                   :: buf
    logical                                      :: received
    integer                                      :: sz, prio
    received = .false.
    if( allocated(buffer) ) deallocate(buffer)
  end function receive_2

#endif

  subroutine get_attributes( self, attr )
    class(ipc_mq),     intent(inout) :: self
    type(ipc_mq_attr), intent(out)   :: attr
    attr = self%attr
  end subroutine get_attributes

  subroutine print_attributes( self )
    class(ipc_mq), intent(inout) :: self
    if( .not. self%active ) return 
    write(logfhandle, '(A,I0)') 'MQ Flags            :', self%attr%mq_flags
    write(logfhandle, '(A,I0)') 'MQ Max Messages     :', self%attr%mq_maxmsg
    write(logfhandle, '(A,I0)') 'MQ Max Message Size :', self%attr%mq_msgsize
    write(logfhandle, '(A,I0)') 'MQ Current Messages :', self%attr%mq_curmsgs
  end subroutine print_attributes

  function get_queue_length( self ) result( queue_length )
    class(ipc_mq),       intent(inout) :: self 
    integer(kind=c_long)               :: queue_length
    queue_length = self%attr%mq_curmsgs
  end function get_queue_length

  function is_active( self ) result( active )
    class(ipc_mq), intent(inout) :: self 
    logical                      :: active
    active = self%active
  end function is_active

end module simple_ipc_mq