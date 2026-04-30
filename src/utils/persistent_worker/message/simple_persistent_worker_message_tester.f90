!@descr: unit tests for persistent-worker wire message modules
!==============================================================================
! MODULE: simple_persistent_worker_message_tester
!
! PURPOSE:
!   Exercises all persistent-worker wire message modules:
!     1. simple_persistent_worker_message_types
!     2. simple_persistent_worker_message_base
!     3. simple_persistent_worker_message_heartbeat
!     4. simple_persistent_worker_message_task
!     5. simple_persistent_worker_message_status
!     6. simple_persistent_worker_message_terminate
!
! ENTRY POINT:
!   run_all_persistent_worker_message_tests
!
! DEPENDENCIES:
!   simple_test_utils and all simple_persistent_worker_message_* modules
!==============================================================================
module simple_persistent_worker_message_tester
  use simple_test_utils, only: assert_true, assert_int
  use simple_persistent_worker_message_types,     only: WORKER_TERMINATE_MSG, WORKER_HEARTBEAT_MSG, WORKER_TASK_MSG, WORKER_STATUS_MSG
  use simple_persistent_worker_message_base,      only: qsys_persistent_worker_message_base
  use simple_persistent_worker_message_heartbeat, only: qsys_persistent_worker_message_heartbeat
  use simple_persistent_worker_message_task,      only: qsys_persistent_worker_message_task
  use simple_persistent_worker_message_status,    only: qsys_persistent_worker_message_status
  use simple_persistent_worker_message_terminate, only: qsys_persistent_worker_message_terminate
  implicit none

  public  :: run_all_persistent_worker_message_tests
  private
#include "simple_local_flags.inc"

contains

  subroutine run_all_persistent_worker_message_tests()
    write(*,'(A)') '**** running all persistent worker message tests ****'
    call test_message_type_enum_values()
    call test_base_message_new_kill_and_serialise()
    call test_heartbeat_new_kill_and_serialise()
    call test_task_new_kill_and_serialise()
    call test_status_new_kill_and_serialise()
    call test_terminate_new_kill_and_serialise()
    write(*,'(A)') '**** persistent worker message tests done ****'
  end subroutine run_all_persistent_worker_message_tests

  subroutine test_message_type_enum_values()
    write(*,'(A)') 'test_message_type_enum_values'
    call assert_int(1, WORKER_TERMINATE_MSG, 'WORKER_TERMINATE_MSG == 1')
    call assert_int(2, WORKER_HEARTBEAT_MSG, 'WORKER_HEARTBEAT_MSG == 2')
    call assert_int(3, WORKER_TASK_MSG,      'WORKER_TASK_MSG == 3')
    call assert_int(4, WORKER_STATUS_MSG,    'WORKER_STATUS_MSG == 4')
  end subroutine test_message_type_enum_values

  subroutine test_base_message_new_kill_and_serialise()
    type(qsys_persistent_worker_message_base) :: msg, decoded
    character(len=:), allocatable             :: buffer
    write(*,'(A)') 'test_base_message_new_kill_and_serialise'

    msg%msg_type = 99
    call msg%new()   ! intentional no-op in base class
    call assert_int(99, msg%msg_type, 'base new() is no-op and preserves msg_type')

    call msg%kill()  ! intentional no-op in base class
    call assert_int(99, msg%msg_type, 'base kill() is no-op and preserves msg_type')

    call msg%serialise(buffer)
    call assert_int(int(sizeof(msg), kind=4), len(buffer), 'base serialise buffer length equals sizeof(base)')
    decoded = transfer(buffer, decoded)
    call assert_int(msg%msg_type, decoded%msg_type, 'base serialise/transfer roundtrip preserves msg_type')
  end subroutine test_base_message_new_kill_and_serialise

  subroutine test_heartbeat_new_kill_and_serialise()
    type(qsys_persistent_worker_message_heartbeat) :: msg, decoded
    character(len=:), allocatable                  :: buffer
    write(*,'(A)') 'test_heartbeat_new_kill_and_serialise'

    call msg%new()
    call assert_int(WORKER_HEARTBEAT_MSG, msg%msg_type, 'heartbeat new() sets msg_type')
    call assert_int(0, msg%worker_id,                'heartbeat new() resets worker_id')
    call assert_int(0, msg%heartbeat_time,           'heartbeat new() resets heartbeat_time')
    call assert_int(0, msg%nthr_used,                'heartbeat new() resets nthr_used')
    call assert_int(0, msg%nthr_total,               'heartbeat new() resets nthr_total')
    call assert_true(len_trim(msg%worker_uid) == 0,  'heartbeat new() resets worker_uid')

    msg%worker_id      = 7
    msg%heartbeat_time = 123456789
    msg%nthr_used      = 3
    msg%nthr_total     = 8
    msg%worker_uid     = 'nodeA_4242'
    call msg%serialise(buffer)
    call assert_int(int(sizeof(msg), kind=4), len(buffer), 'heartbeat serialise buffer length equals sizeof(heartbeat)')
    decoded = transfer(buffer, decoded)
    call assert_int(WORKER_HEARTBEAT_MSG, decoded%msg_type, 'heartbeat roundtrip msg_type')
    call assert_int(7, decoded%worker_id, 'heartbeat roundtrip worker_id')
    call assert_int(123456789, decoded%heartbeat_time, 'heartbeat roundtrip heartbeat_time')
    call assert_int(3, decoded%nthr_used, 'heartbeat roundtrip nthr_used')
    call assert_int(8, decoded%nthr_total, 'heartbeat roundtrip nthr_total')
    call assert_true(trim(decoded%worker_uid) == 'nodeA_4242', 'heartbeat roundtrip worker_uid')

    call msg%kill()
    call assert_int(0, msg%msg_type,               'heartbeat kill() resets msg_type')
    call assert_int(0, msg%worker_id,              'heartbeat kill() resets worker_id')
    call assert_int(0, msg%heartbeat_time,         'heartbeat kill() resets heartbeat_time')
    call assert_int(0, msg%nthr_used,              'heartbeat kill() resets nthr_used')
    call assert_int(0, msg%nthr_total,             'heartbeat kill() resets nthr_total')
    call assert_true(len_trim(msg%worker_uid) == 0,'heartbeat kill() resets worker_uid')
  end subroutine test_heartbeat_new_kill_and_serialise

  subroutine test_task_new_kill_and_serialise()
    type(qsys_persistent_worker_message_task) :: msg, decoded
    character(len=:), allocatable             :: buffer
    write(*,'(A)') 'test_task_new_kill_and_serialise'

    call msg%new()
    call assert_int(WORKER_TASK_MSG, msg%msg_type, 'task new() sets msg_type')
    call assert_int(0, msg%job_id,                 'task new() resets job_id')
    call assert_int(0, msg%queue_time,             'task new() resets queue_time')
    call assert_int(0, msg%start_time,             'task new() resets start_time')
    call assert_int(0, msg%end_time,               'task new() resets end_time')
    call assert_int(0, msg%exit_code,              'task new() resets exit_code')
    call assert_int(0, msg%nthr,                   'task new() resets nthr')
    call assert_true(.not. msg%submitted,          'task new() resets submitted')
    call assert_true(len_trim(msg%script_path) == 0, 'task new() resets script_path')

    msg%job_id     = 21
    msg%queue_time = 100
    msg%start_time = 101
    msg%end_time   = 102
    msg%exit_code  = 0
    msg%nthr       = 4
    msg%submitted  = .true.
    msg%script_path = '/tmp/task_0021.sh'
    call msg%serialise(buffer)
    call assert_int(int(sizeof(msg), kind=4), len(buffer), 'task serialise buffer length equals sizeof(task)')
    decoded = transfer(buffer, decoded)
    call assert_int(WORKER_TASK_MSG, decoded%msg_type, 'task roundtrip msg_type')
    call assert_int(21, decoded%job_id, 'task roundtrip job_id')
    call assert_int(100, decoded%queue_time, 'task roundtrip queue_time')
    call assert_int(101, decoded%start_time, 'task roundtrip start_time')
    call assert_int(102, decoded%end_time, 'task roundtrip end_time')
    call assert_int(0, decoded%exit_code, 'task roundtrip exit_code')
    call assert_int(4, decoded%nthr, 'task roundtrip nthr')
    call assert_true(decoded%submitted, 'task roundtrip submitted')
    call assert_true(trim(decoded%script_path) == '/tmp/task_0021.sh', 'task roundtrip script_path')

    call msg%kill()
    call assert_int(0, msg%msg_type,                'task kill() resets msg_type')
    call assert_int(0, msg%job_id,                  'task kill() resets job_id')
    call assert_int(0, msg%queue_time,              'task kill() resets queue_time')
    call assert_int(0, msg%start_time,              'task kill() resets start_time')
    call assert_int(0, msg%end_time,                'task kill() resets end_time')
    call assert_int(0, msg%exit_code,               'task kill() resets exit_code')
    call assert_int(0, msg%nthr,                    'task kill() resets nthr')
    call assert_true(.not. msg%submitted,           'task kill() resets submitted')
    call assert_true(len_trim(msg%script_path) == 0,'task kill() resets script_path')
  end subroutine test_task_new_kill_and_serialise

  subroutine test_status_new_kill_and_serialise()
    type(qsys_persistent_worker_message_status) :: msg, decoded
    character(len=:), allocatable               :: buffer
    write(*,'(A)') 'test_status_new_kill_and_serialise'

    call msg%new()
    call assert_int(WORKER_STATUS_MSG, msg%msg_type, 'status new() sets msg_type')
    call assert_int(0, msg%status,                   'status new() resets status')
    call assert_true(len_trim(msg%message) == 0,     'status new() resets message')

    msg%status  = 2
    msg%message = 'busy'
    call msg%serialise(buffer)
    call assert_int(int(sizeof(msg), kind=4), len(buffer), 'status serialise buffer length equals sizeof(status)')
    decoded = transfer(buffer, decoded)
    call assert_int(WORKER_STATUS_MSG, decoded%msg_type, 'status roundtrip msg_type')
    call assert_int(2, decoded%status, 'status roundtrip status')
    call assert_true(trim(decoded%message) == 'busy', 'status roundtrip message')

    call msg%kill()
    call assert_int(0, msg%msg_type,             'status kill() resets msg_type')
    call assert_int(0, msg%status,               'status kill() resets status')
    call assert_true(len_trim(msg%message) == 0, 'status kill() resets message')
  end subroutine test_status_new_kill_and_serialise

  subroutine test_terminate_new_kill_and_serialise()
    type(qsys_persistent_worker_message_terminate) :: msg, decoded
    character(len=:), allocatable                  :: buffer
    write(*,'(A)') 'test_terminate_new_kill_and_serialise'

    call msg%new()
    call assert_int(WORKER_TERMINATE_MSG, msg%msg_type, 'terminate new() sets msg_type')
    call assert_int(0, msg%terminate_time,              'terminate new() resets terminate_time')
    call assert_true(len_trim(msg%reason) == 0,         'terminate new() resets reason')

    msg%terminate_time = 555
    msg%reason         = 'shutdown'
    call msg%serialise(buffer)
    call assert_int(int(sizeof(msg), kind=4), len(buffer), 'terminate serialise buffer length equals sizeof(terminate)')
    decoded = transfer(buffer, decoded)
    call assert_int(WORKER_TERMINATE_MSG, decoded%msg_type, 'terminate roundtrip msg_type')
    call assert_int(555, decoded%terminate_time, 'terminate roundtrip terminate_time')
    call assert_true(trim(decoded%reason) == 'shutdown', 'terminate roundtrip reason')

    call msg%kill()
    call assert_int(0, msg%msg_type,                 'terminate kill() resets msg_type')
    call assert_int(0, msg%terminate_time,           'terminate kill() resets terminate_time')
    call assert_true(len_trim(msg%reason) == 0,      'terminate kill() resets reason')
  end subroutine test_terminate_new_kill_and_serialise

end module simple_persistent_worker_message_tester
