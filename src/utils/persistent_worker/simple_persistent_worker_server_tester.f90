!@descr: unit tests for simple_persistent_worker_server
!==============================================================================
! MODULE: simple_persistent_worker_server_tester
!
! PURPOSE:
!   Exercises persistent_worker_server lifecycle and queue semantics:
!     1. default state after declaration
!     2. queue_task behavior before initialisation
!     3. new()/set_debug()/kill() lifecycle
!     4. queue_task routing across high/norm/low priorities
!     5. unknown priority fallback to normal queue
!
! ENTRY POINT:
!   run_all_persistent_worker_server_tests
!
! DEPENDENCIES:
!   simple_persistent_worker_server, simple_test_utils,
!   simple_persistent_worker_message_task, simple_string
!==============================================================================
module simple_persistent_worker_server_tester
  use simple_persistent_worker_server,       only: persistent_worker_server
  use simple_test_utils,                     only: assert_true, assert_int
  use simple_persistent_worker_message_task, only: qsys_persistent_worker_message_task
  use simple_string,                         only: string
  implicit none

  public  :: run_all_persistent_worker_server_tests
  private
#include "simple_local_flags.inc"

contains

  subroutine run_all_persistent_worker_server_tests()
    write(*,'(A)') '**** running all persistent worker server tests ****'
    call test_server_default_state()
    call test_queue_task_not_initialised()
    call test_server_new_set_debug_and_kill()
    call test_queue_task_priority_routing()
    call test_queue_task_unknown_priority_defaults_to_norm()
    call test_queue_task_norm_queue_full()
    write(*,'(A)') '**** persistent worker server tests done ****'
  end subroutine run_all_persistent_worker_server_tests

  subroutine test_server_default_state()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_server_default_state'
    call assert_int(0, server%get_port(), 'default get_port() == 0')
    call assert_true(.not. server%is_running(), 'default is_running() is false')
    call assert_int(0, server%nthr_workers, 'default nthr_workers == 0')
    call assert_int(0, server%job_count, 'default job_count == 0')
    call assert_true(.not. associated(server%worker_data), 'default worker_data is not associated')
    call assert_true(.not. associated(server%listener_args), 'default listener_args is not associated')
  end subroutine test_server_default_state

  subroutine test_queue_task_not_initialised()
    type(persistent_worker_server)           :: server
    type(qsys_persistent_worker_message_task):: task
    logical                                  :: queued
    write(*,'(A)') 'test_queue_task_not_initialised'
    call task%new()
    queued = server%queue_task(task, string('norm'))
    call assert_true(.not. queued, 'queue_task returns .false. when server is not initialised')
    call assert_int(0, task%job_id, 'task job_id remains 0 when queue_task fails early')
  end subroutine test_queue_task_not_initialised

  subroutine test_server_new_set_debug_and_kill()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_server_new_set_debug_and_kill'

    server%l_debug = .false.
    call server%new(4)

    call assert_true(associated(server%worker_data), 'new() associates worker_data')
    call assert_true(associated(server%listener_args), 'new() associates listener_args')
    call assert_int(4, server%nthr_workers, 'new() stores nthr_workers')

    if( associated(server%worker_data) ) then
      call assert_true(.not. server%worker_data%l_debug, 'new() propagates l_debug to worker_data')
    end if

    call server%set_debug(.true.)
    call assert_true(server%l_debug, 'set_debug() updates server l_debug')
    if( associated(server%worker_data) ) then
      call assert_true(server%worker_data%l_debug, 'set_debug() updates worker_data l_debug')
    end if

    call server%kill()
    call assert_int(0, server%get_port(), 'kill() resets get_port() to 0')
    call assert_true(.not. server%is_running(), 'kill() leaves server not running')
    call assert_true(.not. associated(server%worker_data), 'kill() deassociates worker_data')
    call assert_true(.not. associated(server%listener_args), 'kill() deassociates listener_args')
  end subroutine test_server_new_set_debug_and_kill

  subroutine test_queue_task_priority_routing()
    type(persistent_worker_server)            :: server
    type(qsys_persistent_worker_message_task) :: task_high, task_norm, task_low
    logical                                   :: queued_high, queued_norm, queued_low
    write(*,'(A)') 'test_queue_task_priority_routing'

    call server%new(2)
    if( .not. associated(server%worker_data) ) then
      call assert_true(.false., 'new() must associate worker_data for queue tests')
      return
    end if

    call task_high%new()
    task_high%script_path = '/tmp/high.sh'
    queued_high = server%queue_task(task_high, string('high'))
    call assert_true(queued_high, 'queue_task(high) succeeds')
    call assert_int(1, task_high%job_id, 'first queued task gets job_id 1')
    call assert_int(1, server%worker_data%tasks_priority_high(1)%job_id, 'high-priority task routed to high queue slot 1')

    call task_norm%new()
    task_norm%script_path = '/tmp/norm.sh'
    queued_norm = server%queue_task(task_norm, string('norm'))
    call assert_true(queued_norm, 'queue_task(norm) succeeds')
    call assert_int(2, task_norm%job_id, 'second queued task gets job_id 2')
    call assert_int(2, server%worker_data%tasks_priority_norm(1)%job_id, 'normal-priority task routed to norm queue slot 1')

    call task_low%new()
    task_low%script_path = '/tmp/low.sh'
    queued_low = server%queue_task(task_low, string('low'))
    call assert_true(queued_low, 'queue_task(low) succeeds')
    call assert_int(3, task_low%job_id, 'third queued task gets job_id 3')
    call assert_int(3, server%worker_data%tasks_priority_low(1)%job_id, 'low-priority task routed to low queue slot 1')

    call server%kill()
  end subroutine test_queue_task_priority_routing

  subroutine test_queue_task_unknown_priority_defaults_to_norm()
    type(persistent_worker_server)            :: server
    type(qsys_persistent_worker_message_task) :: task
    logical                                   :: queued
    write(*,'(A)') 'test_queue_task_unknown_priority_defaults_to_norm'

    call server%new(1)
    if( .not. associated(server%worker_data) ) then
      call assert_true(.false., 'new() must associate worker_data for unknown-priority test')
      return
    end if

    call task%new()
    task%script_path = '/tmp/default_norm.sh'
    queued = server%queue_task(task, string('unexpected_priority'))
    call assert_true(queued, 'queue_task(unknown_priority) succeeds')
    call assert_int(1, task%job_id, 'unknown-priority task gets job_id 1')
    call assert_int(1, server%worker_data%tasks_priority_norm(1)%job_id, 'unknown priority falls back to normal queue')

    call server%kill()
  end subroutine test_queue_task_unknown_priority_defaults_to_norm

  subroutine test_queue_task_norm_queue_full()
    type(persistent_worker_server)            :: server
    type(qsys_persistent_worker_message_task) :: task
    logical                                   :: queued
    integer                                   :: i, qsize
    write(*,'(A)') 'test_queue_task_norm_queue_full'

    call server%new(1)
    if( .not. associated(server%worker_data) ) then
      call assert_true(.false., 'new() must associate worker_data for queue-full test')
      return
    end if

    qsize = size(server%worker_data%tasks_priority_norm)

    ! Fill the normal-priority queue to capacity.
    do i = 1, qsize
      call task%new()
      task%script_path = '/tmp/norm_full_fill.sh'
      queued = server%queue_task(task, string('norm'))
      call assert_true(queued, 'queue_task(norm) succeeds while filling queue')
      call assert_int(i, task%job_id, 'job_id increments while filling normal queue')
    end do

    ! One extra insertion must fail because the normal-priority queue is full.
    call task%new()
    task%script_path = '/tmp/norm_full_overflow.sh'
    queued = server%queue_task(task, string('norm'))
    call assert_true(.not. queued, 'queue_task(norm) returns .false. when normal queue is full')
    call assert_int(qsize + 1, task%job_id, 'job_id still increments even when queue insertion fails')

    call server%kill()
  end subroutine test_queue_task_norm_queue_full

end module simple_persistent_worker_server_tester
