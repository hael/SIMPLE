!@descr: unit tests for stream_pipe: frames round-tripped through a real non-blocking pipe in one process
module simple_stream_pipe_tester
use, intrinsic :: iso_c_binding, only: c_char, c_int, c_size_t, c_loc, c_sizeof
use unix,                        only: c_pipe, c_close, c_write, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_stream_pipe,          only: stream_pipe
implicit none
private
public :: run_all_stream_pipe_tests

integer, parameter :: MAX_FRAME = 1024

contains

    subroutine run_all_stream_pipe_tests()
        write(*,'(A)') '**** running all stream_pipe tests ****'
        call test_single_frame_round_trip()
        call test_queued_frames_all_delivered()
        call test_frame_split_across_writes()
        call test_empty_and_closed_ends()
        call test_bad_length_resyncs()
        call test_discard()
    end subroutine run_all_stream_pipe_tests

    subroutine test_single_frame_round_trip()
        type(stream_pipe)             :: pipe
        character(len=:), allocatable :: buffer
        integer(c_int)                :: fds(2)
        write(*,'(A)') 'test_single_frame_round_trip'
        call open_loopback(fds, pipe)
        call pipe%send('hello stream')
        call assert_true(pipe%receive(buffer), 'a sent frame is received')
        if( allocated(buffer) ) call assert_char('hello stream', buffer, 'the payload survives the round trip')
        call assert_false(pipe%receive(buffer), 'nothing is received after the only frame')
        call close_loopback(fds, pipe)
    end subroutine test_single_frame_round_trip

    !> two frames written before the first receive arrive in one read; the second must
    !! still be delivered by the next receive, although the pipe is then empty (EAGAIN)
    subroutine test_queued_frames_all_delivered()
        type(stream_pipe)             :: pipe
        character(len=:), allocatable :: buffer
        integer(c_int)                :: fds(2)
        write(*,'(A)') 'test_queued_frames_all_delivered'
        call open_loopback(fds, pipe)
        call pipe%send('first')
        call pipe%send('second')
        call assert_true(pipe%receive(buffer), 'the first queued frame is received')
        if( allocated(buffer) ) call assert_char('first', buffer, 'frames keep their order (first)')
        call assert_true(pipe%receive(buffer), 'the second queued frame is received from the buffer')
        if( allocated(buffer) ) call assert_char('second', buffer, 'frames keep their order (second)')
        call assert_false(pipe%receive(buffer), 'nothing is received after the queued frames')
        call close_loopback(fds, pipe)
    end subroutine test_queued_frames_all_delivered

    !> a frame whose header and payload arrive in pieces is returned once, when complete
    subroutine test_frame_split_across_writes()
        character(len=*), parameter :: PAYLOAD = 'split frame'
        type(stream_pipe)                           :: pipe
        character(len=:),       allocatable         :: buffer
        character(kind=c_char), allocatable, target :: frame(:)
        integer(c_int)    :: fds(2), msg_len
        integer(c_size_t) :: nwritten
        integer           :: header_bytes, nframe, nrest
        write(*,'(A)') 'test_frame_split_across_writes'
        call open_loopback(fds, pipe)
        msg_len      = int(len(PAYLOAD), c_int)
        header_bytes = int(c_sizeof(msg_len))
        nframe       = header_bytes + len(PAYLOAD)
        allocate(frame(nframe))
        frame(1:header_bytes)  = transfer(msg_len, frame(1:header_bytes))
        frame(header_bytes+1:) = transfer(PAYLOAD, frame(1:len(PAYLOAD)))
        ! two bytes of the header
        nwritten = c_write(fds(2), c_loc(frame(1)), 2_c_size_t)
        call assert_int(2, int(nwritten), 'the first piece is written')
        call assert_false(pipe%receive(buffer), 'half a header is not a frame')
        ! the rest of the header and two payload bytes
        nwritten = c_write(fds(2), c_loc(frame(3)), int(header_bytes, c_size_t))
        call assert_int(header_bytes, int(nwritten), 'the second piece is written')
        call assert_false(pipe%receive(buffer), 'a header without all its payload is not a frame')
        ! the remaining payload
        nrest    = nframe - 2 - header_bytes
        nwritten = c_write(fds(2), c_loc(frame(3+header_bytes)), int(nrest, c_size_t))
        call assert_int(nrest, int(nwritten), 'the last piece is written')
        call assert_true(pipe%receive(buffer), 'the frame is received once every byte has arrived')
        if( allocated(buffer) ) call assert_char(PAYLOAD, buffer, 'the reassembled payload is intact')
        deallocate(frame)
        call close_loopback(fds, pipe)
    end subroutine test_frame_split_across_writes

    !> a header announcing more than max_frame_bytes (a writer that gave up mid-frame) drops the
    !! buffered bytes with a warning; the next whole frame is received
    subroutine test_bad_length_resyncs()
        type(stream_pipe)                           :: pipe
        character(len=:),       allocatable         :: buffer
        character(kind=c_char), allocatable, target :: header(:)
        integer(c_int)    :: fds(2), bad_len
        integer(c_size_t) :: nwritten
        integer           :: header_bytes
        write(*,'(A)') 'test_bad_length_resyncs'
        call open_loopback(fds, pipe)
        bad_len      = int(MAX_FRAME + 1, c_int)
        header_bytes = int(c_sizeof(bad_len))
        allocate(header(header_bytes))
        header   = transfer(bad_len, header)
        nwritten = c_write(fds(2), c_loc(header(1)), int(header_bytes, c_size_t))
        call assert_int(header_bytes, int(nwritten), 'the bad header is written')
        call assert_false(pipe%receive(buffer), 'a bad length is not a frame')
        call assert_int(0, pipe%get_pending_bytes(), 'and the buffered bytes are dropped')
        call pipe%send('after the bad frame')
        call assert_true(pipe%receive(buffer), 'the next frame is received')
        if( allocated(buffer) ) call assert_char('after the bad frame', buffer, 'intact')
        deallocate(header)
        call close_loopback(fds, pipe)
    end subroutine test_bad_length_resyncs

    !> discard drops the frames waiting in the pipe and a partial one in the buffer
    subroutine test_discard()
        type(stream_pipe)                           :: pipe
        character(len=:),       allocatable         :: buffer
        character(kind=c_char), allocatable, target :: header(:)
        integer(c_int)    :: fds(2), msg_len
        integer(c_size_t) :: nwritten
        write(*,'(A)') 'test_discard'
        call open_loopback(fds, pipe)
        call pipe%send('stale 1')
        call pipe%send('stale 2')
        ! a header without its payload, left in the buffer by a receive
        msg_len = 5_c_int
        allocate(header(int(c_sizeof(msg_len))))
        header   = transfer(msg_len, header)
        call assert_true(pipe%receive(buffer), 'one stale frame is received')
        nwritten = c_write(fds(2), c_loc(header(1)), int(size(header), c_size_t))
        call pipe%discard()
        call assert_int(0, pipe%get_pending_bytes(), 'nothing is left buffered')
        call assert_false(pipe%receive(buffer),      'nor in the pipe')
        call pipe%send('fresh')
        call assert_true(pipe%receive(buffer), 'a frame sent after the discard is received')
        if( allocated(buffer) ) call assert_char('fresh', buffer, 'whole')
        deallocate(header)
        call close_loopback(fds, pipe)
    end subroutine test_discard

    subroutine test_empty_and_closed_ends()
        type(stream_pipe)             :: pipe
        character(len=:), allocatable :: buffer
        integer(c_int)                :: fds(2)
        write(*,'(A)') 'test_empty_and_closed_ends'
        call open_loopback(fds, pipe)
        call assert_false(pipe%receive(buffer), 'an empty pipe yields no frame')
        call assert_false(allocated(buffer),    'the buffer stays unallocated without a frame')
        call close_loopback(fds, pipe)
        call pipe%new(-1, -1, MAX_FRAME, 'closed')
        call pipe%send('ignored')
        call assert_false(pipe%receive(buffer), 'a missing read end yields no frame')
        call pipe%kill()
        call pipe%kill() ! idempotence
    end subroutine test_empty_and_closed_ends

    ! a pipe with a non-blocking read end, and a stream_pipe reading one end and writing the other
    subroutine open_loopback( fds, pipe )
        integer(c_int),    intent(inout) :: fds(2)
        type(stream_pipe), intent(inout) :: pipe
        integer(c_int) :: rc, flags
        fds   = -1
        rc    = c_pipe(fds)
        call assert_int(0, int(rc), 'pipe() succeeds')
        flags = c_fcntl(fds(1), F_GETFL, 0_c_int)
        rc    = c_fcntl(fds(1), F_SETFL, ior(flags, O_NONBLOCK))
        call assert_int(0, int(rc), 'the read end is made non-blocking')
        call pipe%new(int(fds(1)), int(fds(2)), MAX_FRAME, 'test')
    end subroutine open_loopback

    subroutine close_loopback( fds, pipe )
        integer(c_int),    intent(inout) :: fds(2)
        type(stream_pipe), intent(inout) :: pipe
        integer(c_int) :: rc
        call pipe%kill()
        rc = c_close(fds(1))
        call assert_int(0, int(rc), 'the read end closes')
        rc = c_close(fds(2))
        call assert_int(0, int(rc), 'the write end closes')
        fds = -1
    end subroutine close_loopback

end module simple_stream_pipe_tester
