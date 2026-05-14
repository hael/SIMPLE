!@descr: unit tests for simple_persistent_worker_server
!==============================================================================
! MODULE: simple_persistent_worker_server_tester
!
! PURPOSE:
!   Validates the public contract of simple_persistent_worker_server:
!     - module-level defaults and constants
!     - server lifecycle (new/kill/is_running/get_port/get_host_ips)
!     - invalid-input handling in new()
!     - idempotent new() behavior while already running
!     - queue_task() request/response semantics for priorities
!     - queue saturation behavior (eventual rejection)
!
! ENTRY POINT:
!   run_all_persistent_worker_server_tests
!==============================================================================
module simple_persistent_worker_server_tester
  use simple_persistent_worker_server,       only: persistent_worker, persistent_worker_server, TCP_BUFSZ
  use simple_test_utils,                     only: assert_true, assert_false, assert_int, assert_string_eq
  use simple_persistent_worker_message_task, only: qsys_persistent_worker_message_task
  use simple_string,                         only: string
  implicit none

  integer, parameter :: TEST_NWORKERS           = 1
  integer, parameter :: TEST_NTHR_WORKERS       = 1
  integer, parameter :: MAX_QUEUE_FILL_ATTEMPTS = 1200

  public  :: run_all_persistent_worker_server_tests
  private
#include "simple_local_flags.inc"

contains

  subroutine run_all_persistent_worker_server_tests()
#if defined(_WIN32) || defined(__FreeBSD__)
    write(*,'(A)') '**** skipping persistent worker server tests on this platform ****'
#else
    write(*,'(A)') '**** running all persistent worker server tests ****'

    call test_module_defaults_and_constants()
    call test_server_default_state()
    call test_get_host_ips_default_empty()
    call test_new_invalid_args_do_not_start()
    call test_new_client_only_valid_address_parse()
    call test_new_client_only_invalid_address_rejected()
    call test_new_and_kill_lifecycle()
    call test_new_is_idempotent_while_running()
    call test_queue_task_priority_requests()
    call test_queue_task_eventually_rejects_when_full()

    write(*,'(A)') '**** persistent worker server tests done ****'
#endif
  end subroutine run_all_persistent_worker_server_tests

  subroutine test_module_defaults_and_constants()
    write(*,'(A)') 'test_module_defaults_and_constants'
    ! These constants are part of the server-worker wire protocol; their
    ! values must be fixed and consistent with the worker-side definitions.
    call assert_int(1460, TCP_BUFSZ, 'TCP_BUFSZ should be 1460 bytes')
    call assert_int(0, persistent_worker%n_workers, 'persistent_worker singleton n_workers default should be 0')
    call assert_int(0, persistent_worker%nthr_per_worker, 'persistent_worker singleton nthr_per_worker default should be 0')
    call assert_false(associated(persistent_worker%server), 'persistent_worker singleton server pointer default should be null')
  end subroutine test_module_defaults_and_constants

  subroutine test_server_default_state()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_server_default_state'
    call assert_int(0, server%get_port(), 'default get_port() should be 0')
    call assert_false(server%is_running(), 'default is_running() should be false')
    call assert_int(0, server%n_workers, 'default n_workers should be 0')
    call assert_int(0, server%nthr_workers, 'default nthr_workers should be 0')
    call assert_int(0, server%job_count, 'default job_count should be 0')
    call assert_false(associated(server%worker_data), 'default worker_data should not be associated')
    call assert_false(associated(server%listener_args), 'default listener_args should not be associated')
  end subroutine test_server_default_state

  subroutine test_get_host_ips_default_empty()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_get_host_ips_default_empty'
    call assert_string_eq('', server%get_host_ips(), 'default get_host_ips() should be empty')
  end subroutine test_get_host_ips_default_empty

  subroutine test_new_invalid_args_do_not_start()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_new_invalid_args_do_not_start'
    ! Invalid arguments to new() should not start the server or allocate resources.
    call server%new(0, TEST_NTHR_WORKERS)
    call assert_false(server%is_running(), 'new(n_workers=0,...) should not start the server')
    call assert_int(0, server%get_port(), 'port should remain 0 when new() rejects n_workers=0')
    call assert_false(associated(server%worker_data), 'worker_data should remain unassociated for invalid n_workers')
    call assert_false(associated(server%listener_args), 'listener_args should remain unassociated for invalid n_workers')
    call server%new(TEST_NWORKERS, 0)
    call assert_false(server%is_running(), 'new(...,nthr_workers=0) should not start the server')
    call assert_int(0, server%get_port(), 'port should remain 0 when new() rejects nthr_workers=0')
    call assert_false(associated(server%worker_data), 'worker_data should remain unassociated for invalid nthr_workers')
    call assert_false(associated(server%listener_args), 'listener_args should remain unassociated for invalid nthr_workers')
    call server%kill()
  end subroutine test_new_invalid_args_do_not_start

  subroutine test_new_and_kill_lifecycle()
    type(persistent_worker_server) :: server
    integer                        :: port
    write(*,'(A)') 'test_new_and_kill_lifecycle'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS)
    call assert_true(server%is_running(), 'new() should start the listener thread')
    port = server%get_port()
    call assert_true(port > 0, 'new() should publish a positive listen port')
    call assert_true(associated(server%worker_data), 'new() should associate worker_data')
    call assert_true(associated(server%listener_args), 'new() should associate listener_args')
    call assert_int(TEST_NWORKERS, server%n_workers, 'new() should persist requested n_workers')
    call assert_int(TEST_NTHR_WORKERS, server%nthr_workers, 'new() should persist requested nthr_workers')
    call server%kill()
    call assert_false(server%is_running(), 'kill() should stop the listener thread')
    call assert_int(0, server%get_port(), 'kill() should reset port to 0')
    call assert_false(associated(server%worker_data), 'kill() should deassociate worker_data')
    call assert_false(associated(server%listener_args), 'kill() should deassociate listener_args')
    call assert_int(0, server%nthr_workers, 'kill() should reset nthr_workers to 0')
    call assert_int(0, server%job_count, 'kill() should reset job_count to 0')
  end subroutine test_new_and_kill_lifecycle

  subroutine test_new_client_only_valid_address_parse()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_new_client_only_valid_address_parse'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS, string('127.0.0.1:34567'))
    call assert_int(34567, server%get_port(), 'client_only new() should parse and publish port')
    call assert_string_eq('127.0.0.1', server%get_host_ips(), 'client_only new() should parse and publish host ips')
    call server%kill()
  end subroutine test_new_client_only_valid_address_parse

  subroutine test_new_client_only_invalid_address_rejected()
    type(persistent_worker_server) :: server
    write(*,'(A)') 'test_new_client_only_invalid_address_rejected'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS, string('invalid_address_without_port'))
    call assert_int(0, server%get_port(), 'client_only new() should reject invalid address without port separator')
    call assert_string_eq('', server%get_host_ips(), 'client_only new() should clear host ips on invalid address')
    call server%kill()
  end subroutine test_new_client_only_invalid_address_rejected

  subroutine test_new_is_idempotent_while_running()
    type(persistent_worker_server) :: server
    integer                        :: first_port
    write(*,'(A)') 'test_new_is_idempotent_while_running'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS)
    call assert_true(server%is_running(), 'server should be running before idempotency check')
    first_port = server%get_port()
    call server%new(TEST_NWORKERS + 1, TEST_NTHR_WORKERS + 1)
    call assert_true(server%is_running(), 'second new() call should leave server running')
    call assert_int(first_port, server%get_port(), 'second new() call should not rebind or change port')
    call assert_int(TEST_NWORKERS, server%n_workers, 'second new() call should keep original n_workers')
    call assert_int(TEST_NTHR_WORKERS, server%nthr_workers, 'second new() call should keep original nthr_workers')
    call server%kill()
  end subroutine test_new_is_idempotent_while_running

  subroutine test_queue_task_priority_requests()
    type(persistent_worker_server)            :: server
    type(qsys_persistent_worker_message_task) :: task
    logical                                   :: queued
    write(*,'(A)') 'test_queue_task_priority_requests'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS)
    call assert_true(server%is_running(), 'server must be running for queue_task priority tests')
    call build_test_task(task, '/tmp/pw_high.sh', 1)
    queued = server%queue_task(task, string('high'))
    call assert_true(queued, 'queue_task(high) should be accepted while queue has capacity')
    call build_test_task(task, '/tmp/pw_norm.sh', 1)
    queued = server%queue_task(task, string('norm'))
    call assert_true(queued, 'queue_task(norm) should be accepted while queue has capacity')
    call build_test_task(task, '/tmp/pw_low.sh', 1)
    queued = server%queue_task(task, string('low'))
    call assert_true(queued, 'queue_task(low) should be accepted while queue has capacity')
    call build_test_task(task, '/tmp/pw_unknown.sh', 1)
    queued = server%queue_task(task, string('unexpected_priority'))
    call assert_true(queued, 'queue_task(unknown priority) should be accepted while queue has capacity')
    call server%kill()
  end subroutine test_queue_task_priority_requests

  subroutine test_queue_task_eventually_rejects_when_full()
    type(persistent_worker_server)            :: server
    type(qsys_persistent_worker_message_task) :: task
    logical                                   :: queued
    integer                                   :: i, success_count, first_failure_index
    write(*,'(A)') 'test_queue_task_eventually_rejects_when_full'
    call server%new(TEST_NWORKERS, TEST_NTHR_WORKERS)
    call assert_true(server%is_running(), 'server must be running for queue saturation test')
    success_count        = 0
    first_failure_index  = 0
    do i = 1, MAX_QUEUE_FILL_ATTEMPTS
      call build_test_task(task, '/tmp/pw_fill.sh', 1)
      queued = server%queue_task(task, string('norm'))
      if( queued ) then
        success_count = success_count + 1
      else
        first_failure_index = i
        exit
      end if
    end do
    call assert_true(success_count > 0, 'queue saturation test should accept at least one task before rejection')
    call assert_true(first_failure_index > 0, 'queue saturation test should eventually observe queue_task rejection')
    if( first_failure_index > 0 ) then
      call assert_int(success_count + 1, first_failure_index, 'first failure should follow the last accepted insertion')
    end if
    call server%kill()
  end subroutine test_queue_task_eventually_rejects_when_full

  subroutine build_test_task(task, script_path, nthr)
    type(qsys_persistent_worker_message_task), intent(inout) :: task
    character(len=*),                          intent(in)    :: script_path
    integer,                                   intent(in)    :: nthr
    call task%new()
    task%script_path = script_path
    task%nthr        = nthr
  end subroutine build_test_task

end module simple_persistent_worker_server_tester
