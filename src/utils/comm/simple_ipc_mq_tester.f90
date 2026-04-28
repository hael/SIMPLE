!@descr: unit tests for simple_ipc_mq (lifecycle, I/O, timed receive, throughput)
!==============================================================================
! MODULE: simple_ipc_mq_tester
!
! PURPOSE:
!   Exercises the ipc_mq type across eleven test cases:
!     1.  test_create_and_kill          — basic lifecycle (new/kill/is_active)
!     2.  test_create_and_kill_maxmsg   — lifecycle with explicit max_msgsize
!     3.  test_read_write               — send/receive round-trip (string overload)
!     4.  test_read_write_maxmsg        — round-trip with a size-matched queue
!     5.  test_read_write_buffer        — send/receive round-trip (char-buffer overload)
!     6.  test_read_timed_write         — timed receive, string overload, 10 s deadline
!     7.  test_read_timed_write_buffer  — timed receive, char-buffer overload, 10 s deadline
!     8.  test_send_timed_write         — timed send, string overload, 10 s deadline
!     9.  test_send_timed_write_buffer  — timed send, char-buffer overload, 10 s deadline
!     10. test_print_attributes         — smoke test: print_attributes does not crash
!     11. test_performance              — fill/drain cycle; warns if throughput
!                                         falls below MQ_SEND/RECEIVE_RATE_WARN
!
! ENTRY POINT:
!   run_all_ipc_mq_tests — all tests are Linux-only; the entry point is a
!                          no-op on other platforms.
!
! PARAMETERS (hard-coded):
!   MQ_SEND_RATE_WARN    — minimum acceptable send    rate (msg/s)  (10 000)
!   MQ_RECEIVE_RATE_WARN — minimum acceptable receive rate (msg/s)  (10 000)
!
! DEPENDENCIES:
!   simple_ipc_mq, simple_string, simple_string_utils,
!   simple_test_utils, simple_timer, simple_error
!==============================================================================
module simple_ipc_mq_tester
  use simple_ipc_mq,      only: ipc_mq, ipc_mq_attr
  use simple_string,      only: string
  use simple_string_utils,only: int2str
  use simple_test_utils,  only: assert_true, assert_int, assert_char
  use simple_timer,       only: timer_int_kind, tic, toc
  use simple_error,       only: simple_exception

  implicit none

  integer, parameter :: MQ_SEND_RATE_WARN    = 10000  ! msg/s below which a send    warning is issued
  integer, parameter :: MQ_RECEIVE_RATE_WARN = 10000  ! msg/s below which a receive warning is issued

  public  :: run_all_ipc_mq_tests
  private
#include "simple_local_flags.inc"

contains

  ! Run all ipc_mq unit tests. Linux-only; compiled to a no-op elsewhere.
  subroutine run_all_ipc_mq_tests()
#ifdef __linux__
    if( .not. is_linux_runtime() ) return
    write(*,'(A)') '**** running all ipc mq tests ****'
    call test_create_and_kill()
    call test_create_and_kill_maxmsg()
    call test_read_write()
    call test_read_write_maxmsg()
    call test_read_write_buffer()
    call test_read_timed_write()
    call test_read_timed_write_buffer()
    call test_send_timed_write()
    call test_send_timed_write_buffer()
    call test_print_attributes()
    call test_performance()
#endif
  end subroutine run_all_ipc_mq_tests

  logical function is_linux_runtime()
    character(len=64) :: envval
    integer           :: status, length

    call get_environment_variable('OS', value=envval, length=length, status=status)
    if( status == 0 .and. length > 0 )then
      is_linux_runtime = index(adjustl(envval(:length)), 'Windows_NT') == 0
    else
      is_linux_runtime = .true.
    endif
  end function is_linux_runtime

  ! Verify that new() activates the queue and kill() deactivates it.
  subroutine test_create_and_kill()
    type(ipc_mq) :: mq
    write(*,'(A)') 'test_create_and_kill'
    call mq%new(name=string('TEST_CREATE_AND_KILL'))
    call assert_true(mq%is_active(),        'message queue is active')
    call mq%kill()
    call assert_true(.not. mq%is_active(),  'message queue is inactive')
  end subroutine test_create_and_kill

  ! Same as test_create_and_kill but with an explicit max_msgsize.
  subroutine test_create_and_kill_maxmsg()
    type(ipc_mq) :: mq
    write(*,'(A)') 'test_create_and_kill_maxmsg'
    call mq%new(name=string('TEST_CREATE_AND_KILL_MAXMSG'), max_msgsize=1024)
    call assert_true(mq%is_active(),        'message queue is active')
    call mq%kill()
    call assert_true(.not. mq%is_active(),  'message queue is inactive')
  end subroutine test_create_and_kill_maxmsg

  ! Send one string message and verify it is received intact.
  subroutine test_read_write()
    type(ipc_mq) :: mq
    type(string) :: msg_in, msg_out
    write(*,'(A)') 'test_read_write'
    msg_out = 'TEST MESSAGE'
    call mq%new(name=string('TEST_CREATE_READ_WRITE'))
    call assert_true(mq%is_active(),                                       'message queue is active')
    call mq%send(msg=msg_out)
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive(msg=msg_in),                               'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(),                  'sent and received messages match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_read_write

  ! Send one string message on a size-matched queue; also checks that the
  ! queue attribute mq_msgsize matches the message length.
  subroutine test_read_write_maxmsg()
    type(ipc_mq)      :: mq
    type(string)      :: msg_in, msg_out
    type(ipc_mq_attr) :: mq_attr
    write(*,'(A)') 'test_read_write_maxmsg'
    msg_out = 'TEST MESSAGE'
    call mq%new(name=string('TEST_CREATE_READ_WRITE'), max_msgsize=msg_out%strlen())
    call assert_true(mq%is_active(),                                       'message queue is active')
    call mq%get_attributes(mq_attr)
    call assert_true(mq_attr%mq_msgsize == msg_out%strlen(),               'message queue has correct max msgsize')
    call mq%send(msg=msg_out)
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive(msg=msg_in),                               'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(),                  'sent and received messages match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_read_write_maxmsg

  ! Send one message and retrieve it via receive_timed with a 10-second deadline.
  subroutine test_read_timed_write()
    type(ipc_mq)      :: mq
    type(string)      :: msg_in, msg_out
    type(ipc_mq_attr) :: mq_attr
    write(*,'(A)') 'test_read_timed_write'
    msg_out = 'TEST MESSAGE'
    call mq%new(name=string('TEST_CREATE_READ_TIMED_WRITE'), max_msgsize=msg_out%strlen())
    call assert_true(mq%is_active(),                                       'message queue is active')
    call mq%get_attributes(mq_attr)
    call assert_true(mq_attr%mq_msgsize == msg_out%strlen(),               'message queue has correct max msgsize')
    call mq%send(msg=msg_out)
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive_timed(msg=msg_in, timeout_s=10),           'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(),                  'sent and received messages match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_read_timed_write

  ! Send one message via the character-buffer overload (send_2) and verify it
  ! is received intact via the character-buffer overload (receive_2).
  subroutine test_read_write_buffer()
    type(ipc_mq)                  :: mq
    character(len=:), allocatable :: buf_out, buf_in
    write(*,'(A)') 'test_read_write_buffer'
    buf_out = 'TEST BUFFER MESSAGE'
    call mq%new(name=string('TEST_READ_WRITE_BUFFER'))
    call assert_true(mq%is_active(),                                       'message queue is active')
    call mq%send(buffer=buf_out)
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive(buffer=buf_in),                            'message received from queue')
    call assert_char(buf_in, buf_out,                                      'sent and received buffers match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_read_write_buffer

  ! Send one message and retrieve it via receive_timed with a 10-second deadline
  ! using the character-buffer overload (receive_timed_2).
  subroutine test_read_timed_write_buffer()
    type(ipc_mq)                  :: mq
    character(len=:), allocatable :: buf_out, buf_in
    write(*,'(A)') 'test_read_timed_write_buffer'
    buf_out = 'TEST BUFFER MESSAGE'
    call mq%new(name=string('TEST_READ_TIMED_WRITE_BUFFER'), max_msgsize=len(buf_out))
    call assert_true(mq%is_active(),                                       'message queue is active')
    call mq%send(buffer=buf_out)
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive_timed(buffer=buf_in, timeout_s=10),        'message received from queue')
    call assert_char(buf_in, buf_out,                                      'sent and received buffers match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_read_timed_write_buffer

  ! Send one string message via send_timed with a 10-second deadline and
  ! verify it is received intact.
  subroutine test_send_timed_write()
    type(ipc_mq) :: mq
    type(string) :: msg_in, msg_out
    write(*,'(A)') 'test_send_timed_write'
    msg_out = 'TEST MESSAGE'
    call mq%new(name=string('TEST_SEND_TIMED_WRITE'), max_msgsize=msg_out%strlen())
    call assert_true(mq%is_active(),                                       'message queue is active')
    call assert_true(mq%send_timed(msg=msg_out, timeout_s=10),             'timed send succeeds')
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive(msg=msg_in),                               'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(),                  'sent and received messages match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_send_timed_write

  ! Send one character-buffer message via send_timed with a 10-second deadline
  ! and verify it is received intact.
  subroutine test_send_timed_write_buffer()
    type(ipc_mq)                  :: mq
    character(len=:), allocatable :: buf_out, buf_in
    write(*,'(A)') 'test_send_timed_write_buffer'
    buf_out = 'TEST BUFFER MESSAGE'
    call mq%new(name=string('TEST_SEND_TIMED_WRITE_BUFFER'), max_msgsize=len(buf_out))
    call assert_true(mq%is_active(),                                       'message queue is active')
    call assert_true(mq%send_timed(buffer=buf_out, timeout_s=10),          'timed send succeeds')
    call assert_int( int(mq%get_queue_length(), kind=4), 1,                'message sent to queue')
    call assert_true(mq%receive(buffer=buf_in),                            'message received from queue')
    call assert_char(buf_in, buf_out,                                      'sent and received buffers match')
    call mq%kill()
    call assert_true(.not. mq%is_active(),                                 'message queue is inactive')
  end subroutine test_send_timed_write_buffer

  ! Smoke test: verify that print_attributes completes without error on an
  ! active queue.
  subroutine test_print_attributes()
    type(ipc_mq) :: mq
    write(*,'(A)') 'test_print_attributes'
    call mq%new(name=string('TEST_PRINT_ATTRIBUTES'))
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%print_attributes()
    call mq%kill()
    call assert_true(.not. mq%is_active(), 'message queue is inactive')
  end subroutine test_print_attributes

  ! Fill the queue to capacity and drain it, timing each phase independently.
  ! Issues a warning if either send or receive rate falls below the threshold.
  subroutine test_performance()
    type(ipc_mq)            :: mq
    type(string)            :: msg_in, msg_out
    type(ipc_mq_attr)       :: mq_attr
    integer(timer_int_kind) :: t_start
    integer                 :: i, send_rate, receive_rate
    write(*,'(A)') 'test_performance'
    msg_out = 'MEAGRrTJShyqG4XuiyJ1E4zhJXECzM4PY1zXAG5keC8KfPPAYWaNDCkPQqeWR3c6yT3AxJmm1Kk7iCepfRYnDW2yZbURcQci2FSj6QAchnHXfaAV1EKee1g1TtZK5SUUjbzG3azwYAgChHY2z34xZh1KGxbpZ28A4zTerXFqyWjy3LjFeeYF6ZP1UXJcznDLSgBxCwNNU1tNLdkupM2LM6xr5UW9jrwEDeXrRTuVanh4DGAA27pFmRZJSqPECqvUC4WbSN7mUUnV5hBxcYkZ6cvyfc0trbxeLE7YCGchujN2EhhXycMY9SZ5gEZVvg8wEFxxTQKJww27tBBKbhLaWPGqgX68uhfp7XaAA9MjMpZWbcKUZxugj7jKZZevW0Bhat1JUFX9C3JjAkgGL80MpemD9w9WakCXcTZCPjdeQGVeHr91QXXfPikgNPphmWgj5zKqqK2k5DxU4tyakGBiU4A3KxG9QLiFytewceJ1Hz6FqPD2LDx7qhKxwYrwVcUgxTBxBdhhHFFe0rwF30An0tcHWgqWrz2ju6mbj3740xjmSfu5kGh2CTHC1fmKkTZu72nfBA83UBBzghnGbc0SrRAeK4WyijxQZfy35dfHvyvZ49v97b3GNmYkwzGyqveNZh6UzRYYS6cPhDTztxXuWTVCVe2FL0p5HpBmemezrxJwVjZip7SefJAwn1Qt6Mxj1ZDcZPQUQFGtBeBZa6JFe4FAg8LiwadBFMakLAqqFj0XNttu1M8gh7wffVWGjnCk581HYDve8TbZS5qhjNW7KJB3JFF9YXpXWKdkZRzLBzT0biMHQEQBJb6aUmcfdXkpucbqjtfLSfN6NipHxfzKyYKivKWZ3jPHaNE2AaC0wKvjRnL57Un2yZHqZMLgUnPjBz07yQnpBbWFEuxgahheewpJYJv4mGp8AnT4EK1v3wUE1S02iWgfuHqewvrrMMKu0UFPv8zAyRCxzBid6c3ipHRpgSRacmJgwuMeHpLHdGUkY5pVXDgpSbbJmBHZx6Yn'
    call mq%new(name=string('TEST_PERFORMANCE'))
    call assert_true(mq%is_active(),        'message queue is active')
    call mq%get_attributes(mq_attr)
    call assert_true(mq_attr%mq_maxmsg > 0, 'message queue has non-zero capacity')
    ! Fill phase
    t_start = tic()
    do i = 1, mq_attr%mq_maxmsg
      call mq%send(msg=msg_out)
      call assert_int(int(mq%get_queue_length(), kind=4), i, 'message sent to queue')
    end do
    send_rate = ceiling(mq_attr%mq_maxmsg / toc(t_start))
    if( send_rate < MQ_SEND_RATE_WARN ) &
      THROW_WARN('message queue send rate seems low; ' // int2str(send_rate) // ' msg/s — please investigate')
    ! Drain phase
    t_start = tic()
    do i = 1, mq_attr%mq_maxmsg
      call assert_true(mq%receive(msg=msg_in),                                              'message received from queue')
      call assert_int( int(mq%get_queue_length(), kind=4), int(mq_attr%mq_maxmsg, kind=4) - i, 'message removed from queue')
      call assert_char(msg_in%to_char(), msg_out%to_char(),                                 'sent and received messages match')
    end do
    receive_rate = ceiling(mq_attr%mq_maxmsg / toc(t_start))
    if( receive_rate < MQ_RECEIVE_RATE_WARN ) &
      THROW_WARN('message queue receive rate seems low; ' // int2str(receive_rate) // ' msg/s — please investigate')
    call mq%kill()
    call assert_true(.not. mq%is_active(), 'message queue is inactive')
  end subroutine test_performance

end module simple_ipc_mq_tester
