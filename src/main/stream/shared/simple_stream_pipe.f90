!@descr: length-prefixed framed messages over one non-blocking pipe pair between a stream stage and the master
!==============================================================================
! MODULE: simple_stream_pipe
!
! PURPOSE:
!   One implementation of the framing the stream stages and the master use on
!   their pipes, in place of the send_to_*_in_pipe / receive_from_*_out_pipe
!   copies in every stage. A frame is a C int byte count followed by that many
!   payload bytes.
!
! SEMANTICS:
!   - send: retries immediately on EINTR and with a short sleep on
!     EAGAIN/EWOULDBLOCK. A frame is dropped whole (nothing written within the
!     retry budget) or delivered whole: once any byte is in the pipe the rest
!     follows, because a partial frame desynchronises the reader for good.
!   - receive: a frame already assembled from earlier reads is returned before
!     the pipe is touched, and a read that would block still lets buffered
!     bytes complete a frame. A drain loop therefore sees every queued message,
!     which the per-stage copies did not (they returned on EAGAIN first).
!   - A frame length outside 1..max_frame_bytes (a writer that gave up
!     mid-frame) is reported and the buffered bytes are dropped, so the reader
!     resynchronises on a later frame instead of stopping the process.
!   - discard drops whatever the pipe and the buffer hold, for a reader whose
!     writer has been replaced (a restarted stage).
!   - The descriptors belong to the master (simple_stream_state); kill forgets
!     them and never closes them.
!
! TESTS:
!   simple_stream_pipe_tester
!==============================================================================
module simple_stream_pipe
use, intrinsic :: iso_c_binding, only: c_char, c_int, c_size_t, c_loc, c_sizeof
use unix,                        only: c_read, c_write, c_usleep, EAGAIN, EWOULDBLOCK, EINTR
use simple_error,                only: simple_exception
use simple_gui_metadata_base,    only: gui_metadata_base
implicit none

public :: stream_pipe
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_SEND_RETRIES = 200    ! EAGAIN retries before an unsent frame is dropped
integer, parameter :: RETRY_SLEEP_US   = 10000  ! pause between EAGAIN retries
integer, parameter :: READ_CHUNK_BYTES = 65536  ! one read; a pipe holds no more than this by default

type :: stream_pipe
    private
    integer                       :: fd_read         = -1
    integer                       :: fd_write        = -1
    integer                       :: max_frame_bytes = 0
    integer                       :: expected_len    = -1 ! payload bytes of the frame being assembled; -1 before its header
    character(len=:), allocatable :: pending              ! bytes read and not yet returned
    character(len=:), allocatable :: label                ! names the channel in diagnostics
contains
    procedure          :: new
    procedure          :: send
    procedure          :: send_meta
    procedure          :: receive
    procedure          :: discard
    procedure          :: get_pending_bytes
    procedure          :: kill
    procedure, private :: extract_frame
    procedure, private :: consume
end type stream_pipe

contains

    !> @p fd_read and @p fd_write are the ends this process uses (-1 for none);
    !! @p max_frame_bytes bounds the payload of a received frame.
    subroutine new( self, fd_read, fd_write, max_frame_bytes, label )
        class(stream_pipe), intent(inout) :: self
        integer,            intent(in)    :: fd_read, fd_write, max_frame_bytes
        character(len=*),   intent(in)    :: label
        call self%kill()
        if( max_frame_bytes <= 0 ) THROW_HARD('max_frame_bytes must be positive')
        self%fd_read         = fd_read
        self%fd_write        = fd_write
        self%max_frame_bytes = max_frame_bytes
        self%label           = label
    end subroutine new

    !> Sends @p buffer as one frame; a no-op without a write end or with an empty buffer.
    subroutine send( self, buffer )
        class(stream_pipe), intent(inout) :: self
        character(len=*),   intent(in)    :: buffer
        character(kind=c_char), allocatable, target :: frame(:)
        character(len=:),       allocatable         :: warning
        integer(c_int)    :: msg_len
        integer(c_size_t) :: nwritten
        integer           :: header_bytes, nbytes, nframe, sent, nretries, err_no, rc
        if( self%fd_write < 0 ) return
        nbytes = len(buffer)
        if( nbytes <= 0 ) return
        msg_len      = int(nbytes, c_int)
        header_bytes = int(c_sizeof(msg_len))
        nframe       = header_bytes + nbytes
        allocate(frame(nframe))
        frame(1:header_bytes) = transfer(msg_len, frame(1:header_bytes))
        frame(header_bytes+1:) = transfer(buffer, frame(1:nbytes))
        sent     = 0
        nretries = 0
        do while( sent < nframe )
            nwritten = c_write(self%fd_write, c_loc(frame(sent+1)), int(nframe-sent, c_size_t))
            if( nwritten > 0 )then
                sent     = sent + int(nwritten)
                nretries = 0
                cycle
            endif
            err_no = ierrno()
            if( err_no == int(EINTR) ) cycle
            if( err_no == int(EAGAIN) .or. err_no == int(EWOULDBLOCK) )then
                nretries = nretries + 1
                if( sent == 0 .and. nretries > MAX_SEND_RETRIES )then
                    warning = 'dropped a frame on the '//self%label//' pipe: the reader is not draining it'
                    THROW_WARN(warning)
                    exit
                endif
                rc = c_usleep(RETRY_SLEEP_US)
                cycle
            endif
            warning = 'failed to write to the '//self%label//' pipe'
            THROW_WARN(warning)
            exit
        enddo
        deallocate(frame)
    end subroutine send

    !> Sends one GUI metadata object as a frame, once it holds data.
    subroutine send_meta( self, meta )
        class(stream_pipe),       intent(inout) :: self
        class(gui_metadata_base), intent(in)    :: meta
        character(len=:), allocatable :: buffer
        if( .not. meta%assigned() ) return
        call meta%serialise(buffer)
        call self%send(buffer)
    end subroutine send_meta

    !> Returns .true. with the next complete frame in @p buffer; reads what the pipe holds when none is buffered.
    function receive( self, buffer ) result( got )
        class(stream_pipe),            intent(inout) :: self
        character(len=:), allocatable, intent(inout) :: buffer
        logical :: got
        character(kind=c_char), allocatable, target :: raw(:)
        character(len=:),       allocatable         :: chunk
        integer(c_size_t) :: nread
        integer           :: n
        got = .false.
        if( allocated(buffer) ) deallocate(buffer)
        if( self%fd_read < 0 ) return
        got = self%extract_frame(buffer)
        if( got ) return
        allocate(raw(READ_CHUNK_BYTES))
        nread = c_read(self%fd_read, c_loc(raw(1)), int(READ_CHUNK_BYTES, c_size_t))
        ! nread < 0 with EAGAIN/EWOULDBLOCK only means no new bytes: buffered ones may still complete a frame
        if( nread > 0 )then
            n = int(nread)
            allocate(character(len=n) :: chunk)
            chunk = transfer(raw(1:n), chunk)
            if( allocated(self%pending) )then
                self%pending = self%pending//chunk
            else
                self%pending = chunk
            endif
        endif
        deallocate(raw)
        got = self%extract_frame(buffer)
    end function receive

    !> Drops the bytes buffered and those waiting in the pipe, and any partial frame.
    subroutine discard( self )
        class(stream_pipe), intent(inout) :: self
        character(kind=c_char), allocatable, target :: raw(:)
        integer(c_size_t) :: nread
        self%expected_len = -1
        if( allocated(self%pending) ) deallocate(self%pending)
        if( self%fd_read < 0 ) return
        allocate(raw(READ_CHUNK_BYTES))
        do
            nread = c_read(self%fd_read, c_loc(raw(1)), int(READ_CHUNK_BYTES, c_size_t))
            if( nread <= 0 ) exit
        enddo
        deallocate(raw)
    end subroutine discard

    !> Bytes read and not yet returned as a frame.
    integer function get_pending_bytes( self )
        class(stream_pipe), intent(in) :: self
        get_pending_bytes = 0
        if( allocated(self%pending) ) get_pending_bytes = len(self%pending)
    end function get_pending_bytes

    !> Forgets the descriptors and any partial frame.
    subroutine kill( self )
        class(stream_pipe), intent(inout) :: self
        self%fd_read         = -1
        self%fd_write        = -1
        self%max_frame_bytes = 0
        self%expected_len    = -1
        if( allocated(self%pending) ) deallocate(self%pending)
        if( allocated(self%label)   ) deallocate(self%label)
    end subroutine kill

    ! Moves the next complete frame from the pending bytes into buffer.
    function extract_frame( self, buffer ) result( got )
        class(stream_pipe),            intent(inout) :: self
        character(len=:), allocatable, intent(inout) :: buffer
        logical :: got
        character(len=:), allocatable :: error_message
        integer(c_int) :: msg_len
        integer        :: header_bytes
        got = .false.
        if( .not. allocated(self%pending) ) return
        if( self%expected_len < 0 )then
            header_bytes = int(c_sizeof(msg_len))
            if( len(self%pending) < header_bytes ) return
            msg_len           = transfer(self%pending(1:header_bytes), msg_len)
            self%expected_len = int(msg_len)
            if( self%expected_len <= 0 .or. self%expected_len > self%max_frame_bytes )then
                ! a writer that gave up mid-frame: drop what is buffered and resynchronise
                error_message = 'invalid frame length read from the '//self%label//' pipe; dropping buffered bytes'
                THROW_WARN(error_message)
                deallocate(self%pending)
                self%expected_len = -1
                return
            endif
            call self%consume(header_bytes)
            if( .not. allocated(self%pending) ) return
        endif
        if( len(self%pending) < self%expected_len ) return
        buffer = self%pending(1:self%expected_len)
        call self%consume(self%expected_len)
        self%expected_len = -1
        got = .true.
    end function extract_frame

    ! Drops the first nbytes pending bytes.
    subroutine consume( self, nbytes )
        class(stream_pipe), intent(inout) :: self
        integer,            intent(in)    :: nbytes
        if( len(self%pending) == nbytes )then
            deallocate(self%pending)
        else
            self%pending = self%pending(nbytes+1:)
        endif
    end subroutine consume

end module simple_stream_pipe
