module simple_socket_comm
use, intrinsic :: iso_c_binding
implicit none

public :: simple_socket
private

integer(kind=c_int),     parameter :: SOL_SOCKET       = 1
integer(kind=c_int),     parameter :: SO_RCVTIMEO      = 20
integer(kind=c_int),     parameter :: SO_SNDBUF        = 7
integer(kind=c_int),     parameter :: SO_RCVBUF        = 8
integer(kind=c_int),     parameter :: AF_INET          = 2
integer(kind=c_int),     parameter :: AF_INET6         = 10
integer(kind=c_int),     parameter :: AF_UNIX          = 1
integer(kind=c_int),     parameter :: SOCK_STREAM      = 1
integer(kind=c_int),     parameter :: INADDR_ANY       = 0
integer(kind=c_int16_t), parameter :: PORT             = 8099
integer,                 parameter :: TIMEOUT_SECONDS  = 1
integer,                 parameter :: TIMEOUT_USECONDS = 0
integer,                 parameter :: BUFFER_LENGTH    = 1024000

type, bind(c) :: timeval
    integer(kind=c_long) :: seconds
    integer(kind=c_long) :: useconds
end type

type, bind(c) :: in_addr
    integer(kind=c_int32_t) :: s_addr
end type

type, bind(c) :: sockaddr_in
    integer(kind=c_short)   :: sin_family
    integer(kind=c_int16_t) :: sin_port
    type(in_addr)           :: sin_addr
end type

type :: simple_socket
    private
    integer                    :: sock_fd, port
    type(sockaddr_in), pointer :: address
contains
    procedure          :: open
    procedure          :: set_options
    procedure          :: bind_any
    procedure          :: listen
    procedure          :: accept
    procedure, private :: read_1, read_2
    generic            :: read => read_1, read_2
    procedure, private :: send_1
    generic            :: send => send_1
    procedure, private :: close_1, close_2
    generic            :: close => close_1, close_2
end type simple_socket

interface

    function c_socket(domain, type, protocol) bind(c, name="socket")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)          :: c_socket
        integer(kind=c_int),   value :: domain, type, protocol
    end function c_socket

    function c_setsockopt(sockfd, level, optname, optval, optlen) bind(c, name="setsockopt")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)          :: c_setsockopt
        integer(kind=c_int),   value :: sockfd, level, optname, optlen
        type(c_ptr),           value :: optval
    end function c_setsockopt

    function c_bind(sockfd, addrval, addrlen) bind(c, name="bind")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_bind
        integer(kind=c_int),    value :: sockfd
        integer(kind=c_size_t), value :: addrlen
        type(c_ptr),            value :: addrval
    end function c_bind

    function c_listen(sockfd, backlog) bind(c, name="listen")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_listen
        integer(kind=c_int),    value :: sockfd, backlog
    end function c_listen

    function c_accept(sockfd, addrval, addrlen) bind(c, name="accept")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_accept
        integer(kind=c_int),    value :: sockfd
        type(c_ptr),            value :: addrval, addrlen
    end function c_accept

    function c_connect(sockfd, addrval, addrlen) bind(c, name="connect")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_connect
        integer(kind=c_int),    value :: sockfd
        integer(kind=c_size_t), value :: addrlen
        type(c_ptr),            value :: addrval
    end function c_connect

    function c_close(sockfd) bind(c, name="close")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_close
        integer(kind=c_int),    value :: sockfd
    end function c_close

    function c_send(sockfd, bufval, buflen, flags) bind(c, name="send")
        use, intrinsic :: iso_c_binding
        integer(kind=c_size_t)        :: c_send
        integer(kind=c_int),    value :: sockfd, flags
        integer(kind=c_size_t), value :: buflen
        type(c_ptr),            value :: bufval
    end function c_send

    function c_read(sockfd, bufval, buflen) bind(c, name="read")
        use, intrinsic :: iso_c_binding
        integer(kind=c_size_t)        :: c_read
        integer(kind=c_int),    value :: sockfd
        integer(kind=c_size_t), value :: buflen
        type(c_ptr),            value :: bufval
    end function c_read

    function c_htons(i) bind(c, name="htons")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int16_t)        :: c_htons
        integer(kind=c_int16_t), value :: i
    end function c_htons

    function c_inet_pton(af, src, dst) bind(c, name="inet_pton")
        use, intrinsic :: iso_c_binding
        integer(kind=c_int)           :: c_inet_pton
        integer(kind=c_int),    value :: af
        type(c_ptr),            value :: src, dst
    end function c_inet_pton

    function c_errno(i) bind(c, name="perror")
        use, intrinsic :: iso_c_binding
        character(kind=c_char), value :: i
        integer(kind=c_int) :: c_errno
    end function c_errno

end interface

contains

    subroutine open(this, protocol, port, ip)
        class(simple_socket),               intent(inout) :: this
        integer(kind=c_int),      optional, intent(in)    :: protocol
        integer(kind=c_int16_t),  optional, intent(in)    :: port
        character(len=*), target, optional, intent(in)    :: ip
        integer(kind=c_int)             :: l_protocol
        integer(kind=c_int16_t)         :: l_port
        integer(kind=c_int32_t), target :: l_ip
        integer                         :: rc
        l_protocol = 2    !AF_INET
        l_port     = 8099 !PORT
        l_ip       = 0    !INADDR_ANY
        if(present(protocol)) l_protocol = protocol
        if(present(port))     l_port     = port
        if(present(ip)) then
            rc = c_inet_pton(l_protocol, c_loc(ip), c_loc(l_ip))
            if(rc < 1) then
                write(*,*) "failed to convert address to binary"
                write(*,*) c_errno(c_null_char)
                return
            end if
        end if
        allocate(this%address)
        this%address%sin_family      = l_protocol;
        this%address%sin_addr%s_addr = l_ip;
        this%address%sin_port        = c_htons(l_port);
        this%sock_fd = c_socket(l_protocol, SOCK_STREAM, 0)
        if(this%sock_fd < 0) write(*,*) "failed to open socket"
    end subroutine open

    subroutine set_options(this)
        class(simple_socket), intent(inout) :: this
        type(timeval), target :: timeout
        integer(kind=c_int)   :: timeout_size
        integer               :: rc
        timeout%useconds = TIMEOUT_USECONDS
        timeout%seconds  = TIMEOUT_SECONDS
        timeout_size     = 2 * c_long
        rc = c_setsockopt(this%sock_fd, SOL_SOCKET, SO_RCVTIMEO, c_loc(timeout), timeout_size)
        if(rc < 0) write(*,*) "failed to set socket options"
    end subroutine set_options

    subroutine bind_any(this)
        class(simple_socket), intent(inout) :: this
        integer :: rc
        rc = c_bind(this%sock_fd, c_loc(this%address), c_sizeof(this%address) * 2)
        if(rc < 0) then
            write(*,*) "failed to bind socket"
            write(*,*) c_errno(c_null_char)
        endif
    end subroutine bind_any

    subroutine listen(this)
        class(simple_socket), intent(inout) :: this
        integer :: rc
        rc = c_listen(this%sock_fd, 3)
        if(rc < 0) write(*,*) "failed to listen on socket"
    end subroutine listen

    subroutine accept(this, fd)
        class(simple_socket),   intent(inout) :: this
        integer,                intent(inout) :: fd
        integer(kind=c_size_t), target        :: addrlen
        addrlen = c_sizeof(this%address) * 2
        fd = c_accept(this%sock_fd, c_loc(this%address), c_loc(addrlen))
        if(fd < 0) then
            write(*,*) "failed to accept on socket"
            write(*,*) c_errno(c_null_char)
        endif
    end subroutine accept

    subroutine send_1(this, msg)
        class(simple_socket),      intent(inout) :: this
        character (len=*), target, intent(in)    :: msg
        integer(kind=c_size_t),        target    :: msglen
        integer :: rc
        msglen = len(trim(msg))
        rc = c_connect(this%sock_fd, c_loc(this%address), c_sizeof(this%address) * 2)
        if(rc < 0) then
            write(*,*) "failed to connect to socket"
            write(*,*) c_errno(c_null_char)
            return
        endif
        rc = c_send(this%sock_fd, c_loc(msg), msglen, 0)
        if(rc < 0) then
            write(*,*) "failed to send message"
            write(*,*) c_errno(c_null_char)
        endif
    end subroutine send_1

    subroutine read_1(this, ans_str_loc)
        class(simple_socket),          intent(inout) :: this
        character(len=:), allocatable, intent(inout) :: ans_str_loc
        integer(kind=c_size_t),       target :: buflen, buf2len
        character(len=BUFFER_LENGTH), target :: buffer, buffer2
        integer(kind=c_size_t) :: nbytes, msglen
        logical                :: msgstart
        integer                :: n_trials
        msgstart = .false.
        buflen   = BUFFER_LENGTH
        buf2len  = 0
        n_trials = 0
        do while (.true.)
            nbytes = c_read(this%sock_fd, c_loc(buffer), buflen)
            if(buffer(1:5) == "#N1C3" .and. nbytes >= (6 + c_size_t)) then
                msgstart = .true.
                msglen = transfer(buffer(6:5 + c_size_t), 1)
                buffer2 = buffer(6 + c_size_t:nbytes)
                buf2len = nbytes - (5 + c_size_t)
            else if(nbytes > 0 .and. msgstart) then
                buffer2 = buffer2(:buf2len) // buffer(:nbytes)
                buf2len = buf2len + nbytes
            end if
            n_trials = n_trials + 1
            if(buf2len == msglen) exit
            if(n_trials > 5) exit
        end do
        write(*,*) "Received1 : ", buffer2(:msglen), msglen, buf2len, n_trials
        ans_str_loc = buffer2(:msglen)
    end subroutine read_1

    subroutine read_2(this, fd)
        class(simple_socket),  intent(inout) :: this
        integer,               intent(in)    :: fd
        integer(kind=c_size_t),       target :: buflen
        character(len=BUFFER_LENGTH), target :: buffer
        integer(kind=c_size_t) :: nbytes
        buflen = BUFFER_LENGTH
        nbytes = c_read(fd, c_loc(buffer), buflen)
        write(*,*) "Received2 :", buffer(:nbytes)
    end subroutine read_2

    subroutine close_1(this)
        class(simple_socket), intent(inout) :: this
        integer :: rc
        rc = c_close(this%sock_fd)
        if(rc < 0) then
            write(*,*) "failed to close socket"
            write(*,*) c_errno(c_null_char)
        endif
    end subroutine close_1

    subroutine close_2(this, fd)
        class(simple_socket), intent(inout) :: this
        integer,              intent(in)    :: fd
        integer :: rc
        rc = c_close(fd)
        if(rc < 0) then
            write(*,*) "failed to close socket"
            write(*,*) c_errno(c_null_char)
        endif
    end subroutine close_2

end module simple_socket_comm