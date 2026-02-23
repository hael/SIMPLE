!@descr: unit tests for ipc mq module
module simple_ipc_mq_tester
use simple_ipc_mq
use simple_string
use simple_string_utils
use simple_test_utils
use simple_timer
use simple_error

implicit none

integer,  parameter :: MQ_SEND_RATE_WARN    = 10000
integer,  parameter :: MQ_RECEIVE_RATE_WARN = 10000

public :: run_all_ipc_mq_tests
private
#include "simple_local_flags.inc"

contains

  subroutine run_all_ipc_mq_tests()
#ifdef __linux__
    ! message queues only work on linux
    write(*,'(A)') '**** running all ipc mq tests ****'
    call test_create_and_kill()
    call test_create_and_kill_maxmsg()
    call test_read_write()
    call test_read_write_maxmsg()
    call test_performance()
#endif
  end subroutine run_all_ipc_mq_tests

  !---------------- basic lifecycle ----------------

  subroutine test_create_and_kill()
    type(ipc_mq) :: mq
    write(*,'(A)') 'test_create_and_kill'
    call mq%new(name=string('TEST_CREATE_AND_KILL'))
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%kill()
    call assert_true(.not.mq%is_active(), 'message queue is inactive')
  end subroutine test_create_and_kill

  subroutine test_create_and_kill_maxmsg()
    type(ipc_mq) :: mq
    write(*,'(A)') 'test_create_and_kill_maxmsg'
    call mq%new(name=string('TEST_CREATE_AND_KILL_MAXMSG'), max_msgsize=1024)
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%kill()
    call assert_true(.not.mq%is_active(), 'message queue is inactive')
  end subroutine test_create_and_kill_maxmsg

  !---------------- I/O ----------------

  subroutine test_read_write()
    type(ipc_mq) :: mq
    type(string) :: msg_in, msg_out
    write(*,'(A)') 'test_read_write'
    msg_out  = 'TEST MESSAGE'
    call mq%new(name=string('TEST_CREATE_READ_WRITE'))
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%send(msg=msg_out)
    call assert_int(int(mq%get_queue_length(), kind=4), 1,  'message sent to queue')
    call assert_true(mq%receive(msg=msg_in), 'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(), 'sent and received messages match')
    call mq%kill()
    call assert_true(.not.mq%is_active(), 'message queue is inactive')
  end subroutine test_read_write

  subroutine test_read_write_maxmsg()
    type(ipc_mq)      :: mq
    type(string)      :: msg_in, msg_out
    type(ipc_mq_attr) :: mq_attr
    write(*,'(A)') 'test_read_write_maxmsg'
    msg_out  = 'TEST MESSAGE'
    call mq%new(name=string('TEST_CREATE_READ_WRITE'), max_msgsize=msg_out%strlen())
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%get_attributes(mq_attr)
    call assert_true(mq_attr%mq_msgsize == msg_out%strlen(), 'message queue has correct max msgsize')
    call mq%send(msg=msg_out)
    call assert_int(int(mq%get_queue_length(), kind=4), 1,  'message sent to queue')
    call assert_true(mq%receive(msg=msg_in), 'message received from queue')
    call assert_char(msg_in%to_char(), msg_out%to_char(), 'sent and received messages match')
    call mq%kill()
    call assert_true(.not.mq%is_active(), 'message queue is inactive')
  end subroutine test_read_write_maxmsg

  !---------------- Performance ---------------

  subroutine test_performance()
    type(ipc_mq)            :: mq
    type(string)            :: msg_in, msg_out
    type(ipc_mq_attr)       :: mq_attr
    integer(timer_int_kind) :: t_start
    integer                 :: i, send_rate, receive_rate
    write(*,'(A)') 'test_performance'
    msg_out  = 'MEAGRrTJShyqG4XuiyJ1E4zhJXECzM4PY1zXAG5keC8KfPPAYWaNDCkPQqeWR3c6yT3AxJmm1Kk7iCepfRYnDW2yZbURcQci2FSj6QAchnHXfaAV1EKee1g1TtZK5SUUjbzG3azwYAgChHY2z34xZh1KGxbpZ28A4zTerXFqyWjy3LjFeeYF6ZP1UXJcznDLSgBxCwNNU1tNLdkupM2LM6xr5UW9jrwEDeXrRTuVanh4DGAA27pFmRZJSqPECqvUC4WbSN7mUUnV5hBxcYkZ6cvyfc0trbxeLE7YCGchujN2EhhXycMY9SZ5gEZVvg8wEFxxTQKJww27tBBKbhLaWPGqgX68uhfp7XaAA9MjMpZWbcKUZxugj7jKZZevW0Bhat1JUFX9C3JjAkgGL80MpemD9w9WakCXcTZCPjdeQGVeHr91QXXfPikgNPphmWgj5zKqqK2k5DxU4tyakGBiU4A3KxG9QLiFytewceJ1Hz6FqPD2LDx7qhKxwYrwVcUgxTBxBdhhHFFe0rwF30An0tcHWgqWrz2ju6mbj3740xjmSfu5kGh2CTHC1fmKkTZu72nfBA83UBBzghnGbc0SrRAeK4WyijxQZfy35dfHvyvZ49v97b3GNmYkwzGyqveNZh6UzRYYS6cPhDTztxXuWTVCVe2FL0p5HpBmemezrxJwVjZip7SefJAwn1Qt6Mxj1ZDcZPQUQFGtBeBZa6JFe4FAg8LiwadBFMakLAqqFj0XNttu1M8gh7wffVWGjnCk581HYDve8TbZS5qhjNW7KJB3JFF9YXpXWKdkZRzLBzT0biMHQEQBJb6aUmcfdXkpucbqjtfLSfN6NipHxfzKyYKivKWZ3jPHaNE2AaC0wKvjRnL57Un2yZHqZMLgUnPjBz07yQnpBbWFEuxgahheewpJYJv4mGp8AnT4EK1v3wUE1S02iWgfuHqewvrrMMKu0UFPv8zAyRCxzBid6c3ipHRpgSRacmJgwuMeHpLHdGUkY5pVXDgpSbbJmBHZx6Yn'      
    call mq%new(name=string('TEST_PERFORMANCE'))
    call assert_true(mq%is_active(), 'message queue is active')
    call mq%get_attributes(mq_attr)
    call assert_true(mq_attr%mq_maxmsg > 0, 'message queue has non 0 length')
    t_start = tic()
    do i=1, mq_attr%mq_maxmsg
      call mq%send(msg=msg_out)
      call assert_int(int(mq%get_queue_length(), kind=4), i,  'message sent to queue')
    end do
    send_rate = ceiling(mq_attr%mq_maxmsg / toc(t_start))
    if( send_rate < MQ_SEND_RATE_WARN ) THROW_WARN('message queue send rate seems low; ' //  int2str(send_rate) // ' msg/second. please investigate')
    t_start = tic()
    do i=1, mq_attr%mq_maxmsg
      call assert_true(mq%receive(msg=msg_in), 'message received from queue')
      call assert_int(int(mq%get_queue_length(), kind=4), int(mq_attr%mq_maxmsg, kind=4) - i,  'message removed from queue')
      call assert_char(msg_in%to_char(), msg_out%to_char(), 'sent and received messages match')
    end do
    receive_rate = ceiling(mq_attr%mq_maxmsg / toc(t_start))
    if( receive_rate < MQ_RECEIVE_RATE_WARN ) THROW_WARN('message queue receive rate seems low; ' //  int2str(receive_rate) // ' msg/second. please investigate')
    call mq%kill()
    call assert_true(.not.mq%is_active(), 'message queue is inactive')
  end subroutine test_performance

end module simple_ipc_mq_tester