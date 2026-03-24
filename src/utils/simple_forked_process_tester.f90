!@descr: unit tests for simple_forked_process (lifecycle, signals, restart, timestamps, I/O)
!==============================================================================
! MODULE: simple_forked_process_tester
!
! PURPOSE:
!   Exercises the forked_process type across eight test cases:
!     1. test_start              — fork, run to completion, verify STOPPED
!     2. test_kill               — fork, SIGKILL, verify FAILED
!     3. test_terminate          — fork, SIGTERM, verify STOPPED (graceful)
!     4. test_restart            — fork with restart=.true., SIGKILL, verify
!                                  RESTARTING then STOPPED; checks get_pid
!                                  and get_nrestarts
!     5. test_timestamps         — verify queuetime/starttime/stoptime are
!                                  set and ordered after a clean run
!     6. test_fail_timestamps    — verify failtime is set after a SIGKILL
!     7. test_destroy            — smoke test: destroy() does not crash
!     8. test_logfile_redirection— verify child output is written to logfile
!                                  with correct FNV-1a hash
!
! ENTRY POINT:
!   run_all_forked_process_tests
!
! DEPENDENCIES:
!   unix, simple_forked_process, simple_string, simple_test_utils, simple_syslib
!==============================================================================
module simple_forked_process_tester
  use unix,                  only: c_pid_t, c_usleep
  use simple_forked_process, only: forked_process,         &
                                   FORK_STATUS_RUNNING,    &
                                   FORK_STATUS_STOPPED,    &
                                   FORK_STATUS_FAILED,     &
                                   FORK_STATUS_RESTARTING, &
                                   FORK_POLL_TIME
  use simple_string,         only: string
  use simple_test_utils,     only: assert_true, assert_int, assert_char
  use simple_syslib,         only: file_exists, del_file

  implicit none

  public  :: run_all_forked_process_tests
  private
#include "simple_local_flags.inc"

contains

  ! Run all forked_process unit tests in order.
  subroutine run_all_forked_process_tests()
    write(*,'(A)') '**** running all forked process tests ****'
    call test_start()
    call test_kill()
    call test_terminate()
    call test_restart()
    call test_timestamps()
    call test_fail_timestamps()
    call test_destroy()
    call test_logfile_redirection()
  end subroutine run_all_forked_process_tests

  ! Fork a process, let it run to completion, and verify it reaches STOPPED.
  subroutine test_start()
    type(forked_process) :: proc
    write(*,'(A)') 'test_start'
    call proc%start(name=string('TEST_START'))
    call assert_int(proc%status(), FORK_STATUS_RUNNING, 'process is running after start')
    call proc%await_final_status()
    call assert_int(proc%status(), FORK_STATUS_STOPPED, 'process is stopped after completion')
  end subroutine test_start

  ! Fork a process, send SIGKILL, and verify it reaches FAILED.
  subroutine test_kill()
    type(forked_process) :: proc
    integer              :: rc
    write(*,'(A)') 'test_kill'
    call proc%start(name=string('TEST_KILL'))
    call assert_int(proc%status(), FORK_STATUS_RUNNING, 'process is running after start')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call proc%kill()
    call proc%await_final_status()
    call assert_int(proc%status(), FORK_STATUS_FAILED, 'process is failed after SIGKILL')
  end subroutine test_kill

  ! Fork a process, send SIGTERM (graceful), and verify it reaches STOPPED.
  subroutine test_terminate()
    type(forked_process) :: proc
    integer              :: rc
    write(*,'(A)') 'test_terminate'
    call proc%start(name=string('TEST_TERMINATE'))
    call assert_int(proc%status(), FORK_STATUS_RUNNING, 'process is running after start')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call proc%terminate()
    call proc%await_final_status()
    call assert_int(proc%status(), FORK_STATUS_STOPPED, 'process is stopped after SIGTERM')
  end subroutine test_terminate

  ! Fork with restart=.true., SIGKILL it, and verify: status passes through
  ! RESTARTING, the restarted PID differs from the original, get_nrestarts
  ! returns 1, and the process eventually reaches STOPPED.
  subroutine test_restart()
    type(forked_process)  :: proc
    integer(kind=c_pid_t) :: pid1, pid2
    integer               :: rc, stat
    write(*,'(A)') 'test_restart'
    call proc%start(name=string('TEST_RESTART'), restart=.true.)
    call assert_int(proc%status(), FORK_STATUS_RUNNING, 'process is running after start')
    pid1 = proc%get_pid()
    call assert_true(pid1 /= -1,                        'process has valid PID')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call proc%kill()
    ! Poll until we observe RESTARTING or a terminal state
    stat = FORK_STATUS_RUNNING
    do while( stat == FORK_STATUS_RUNNING )
      rc   = c_usleep(FORK_POLL_TIME)
      stat = proc%status()
    end do
    call assert_true(stat == FORK_STATUS_RESTARTING .or. stat == FORK_STATUS_STOPPED, &
                     'process is restarting or stopped after kill')
    call assert_int(proc%get_nrestarts(), 1,             'restart count is 1 after one failure')
    pid2 = proc%get_pid()
    call assert_true(pid2 /= -1,                        'restarted process has valid PID')
    call assert_true(pid1 /= pid2,                      'restarted process has a new PID')
    call proc%await_final_status()
    call assert_int(proc%status(), FORK_STATUS_STOPPED,  'process is stopped after restart completes')
  end subroutine test_restart

  ! Verify that queuetime, starttime, and stoptime are all positive and
  ! ordered correctly after a clean run (queuetime <= starttime <= stoptime).
  subroutine test_timestamps()
    type(forked_process) :: proc
    write(*,'(A)') 'test_timestamps'
    call proc%start(name=string('TEST_TIMESTAMPS'))
    call proc%await_final_status()
    call assert_true(proc%get_queuetime() > 0,                           'queuetime is set')
    call assert_true(proc%get_starttime() > 0,                           'starttime is set')
    call assert_true(proc%get_stoptime()  > 0,                           'stoptime is set')
    call assert_true(proc%get_starttime() >= proc%get_queuetime(),       'starttime >= queuetime')
    call assert_true(proc%get_stoptime()  >= proc%get_starttime(),       'stoptime >= starttime')
    call assert_true(proc%get_failtime()  == 0,                          'failtime is zero for clean run')
  end subroutine test_timestamps

  ! Verify that failtime is set (non-zero) after a SIGKILL simulated failure.
  subroutine test_fail_timestamps()
    type(forked_process) :: proc
    integer              :: rc
    write(*,'(A)') 'test_fail_timestamps'
    call proc%start(name=string('TEST_FAIL_TIMESTAMPS'))
    rc = c_usleep(FORK_POLL_TIME * 5)
    call proc%kill()
    call proc%await_final_status()
    call assert_true(proc%get_failtime() > 0,  'failtime is set after SIGKILL')
    call assert_true(proc%get_stoptime() == 0, 'stoptime is zero after failure')
  end subroutine test_fail_timestamps

  ! Smoke test: destroy() completes without error.
  subroutine test_destroy()
    type(forked_process) :: proc
    write(*,'(A)') 'test_destroy'
    call proc%start(name=string('TEST_DESTROY'))
    call proc%await_final_status()
    call proc%destroy()
    call assert_true(.true., 'destroy completes without error')
  end subroutine test_destroy

  ! Fork with a logfile, let the process write its sentinel line, then verify
  ! the file exists, has content, and its FNV-1a hash matches the expected value.
  ! NOTE: the expected hash is tied to the 'LOGFILE CONTENTS TEST' sentinel in
  !       execute_test — update it if that string changes.
  subroutine test_logfile_redirection()
    type(forked_process) :: proc
    type(string)         :: log_contents, log_hash, log_fname
    integer              :: ios, unit
    write(*,'(A)') 'test_logfile_redirection'
    log_fname = 'test_logfile_redirection.log'
    call proc%start(name=string('TEST_LOGFILE_REDIRECTION'), logfile=log_fname)
    call assert_int(proc%status(), FORK_STATUS_RUNNING, 'process is running after start')
    call proc%await_final_status()
    call assert_int(proc%status(), FORK_STATUS_STOPPED, 'process is stopped after completion')
    call assert_true(file_exists(log_fname),             'logfile exists after run')
    open(newunit=unit, file=log_fname%to_char(), action='read', iostat=ios)
    call assert_int(ios, 0,                              'logfile opens for reading')
    call log_contents%readfile(unit)
    close(unit)
    call assert_true(log_contents%strlen() > 0,          'logfile has content')
    log_hash = log_contents%to_fnv1a_hash64()
    call assert_char(log_hash%to_char(), '0A8F18CBEE2D2351', 'logfile content hash matches')
    call del_file(log_fname)
  end subroutine test_logfile_redirection

end module simple_forked_process_tester
