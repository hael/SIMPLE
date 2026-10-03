!@descr: unit tests for simple_forked_process (lifecycle, signals, restart, timestamps, I/O)
! Children run the default execute_test, which exits 0 on SIGTERM: terminate() must end in
! STOPPED, kill() (SIGKILL) in FAILED. test_logfile_redirection hashes execute_test's sentinel
! line. test_fork_with_running_monitor forks under a running memory monitor. Skipped on Windows.
module simple_forked_process_tester
  use unix,                  only: c_pid_t, c_usleep
  use simple_forked_process, only: forked_process,         &
                                   FORK_STATUS_RUNNING,    &
                                   FORK_STATUS_STOPPED,    &
                                   FORK_STATUS_FAILED,     &
                                   FORK_STATUS_RESTARTING, &
                                   FORK_POLL_TIME
  use simple_cmdline,        only: cmdline
  use simple_memory_monitor, only: mem_monitor_init, mem_monitor_finish, mem_monitor_is_enabled
  use simple_string,         only: string
  use simple_string_utils,   only: int2str
  use simple_test_utils,     only: assert_true, assert_int, assert_char
  use simple_syslib,         only: file_exists, del_file, get_process_id

  implicit none

  public  :: run_all_forked_process_tests
  private
#include "simple_local_flags.inc"

contains

  ! Run all forked_process unit tests in order.
  subroutine run_all_forked_process_tests()
#if !defined(_WIN32)
    write(*,'(A)') '**** running all forked process tests ****'
    call test_start()
    call test_kill()
    call test_terminate()
    call test_restart()
    call test_timestamps()
    call test_fail_timestamps()
    call test_destroy()
    call test_logfile_redirection()
    call test_fork_with_running_monitor()
#endif
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

  ! Fork while this process's memory monitor runs, as the stream master restarts a stage: the
  ! child starts its own monitor (its command line asks for one), runs and stops, and the
  ! parent's monitor still runs. A child stopping the inherited monitor would join a sampler
  ! thread it does not have and never return, so the wait has a deadline.
  subroutine test_fork_with_running_monitor()
    integer, parameter    :: MAX_POLLS = 200 ! 20 s; execute_test sleeps 2 s
    type(forked_process)  :: proc
    type(cmdline)         :: cline
    type(string)          :: log_fname, child_telemetry
    integer(kind=c_pid_t) :: pid
    integer               :: rc, ipoll, stat
    logical               :: l_own_monitor
    write(*,'(A)') 'test_fork_with_running_monitor'
    call cline%set('memreport', 'yes')
    ! the runner's own monitor, when it runs with memreport=yes, serves as well
    l_own_monitor = .not. mem_monitor_is_enabled()
    if( l_own_monitor ) call mem_monitor_init(cline, 'forked process tester')
    if( .not. mem_monitor_is_enabled() )then
      write(*,'(A)') 'memory telemetry is unsupported here; test skipped'
      return
    endif
    log_fname = 'test_fork_with_running_monitor.log'
    call proc%start(name=string('TEST_FORK_WITH_RUNNING_MONITOR'), logfile=log_fname, cline=cline)
    pid   = proc%get_pid()
    stat  = proc%status()
    ipoll = 0
    do while( stat == FORK_STATUS_RUNNING .and. ipoll < MAX_POLLS )
      rc    = c_usleep(FORK_POLL_TIME)
      stat  = proc%status()
      ipoll = ipoll + 1
    end do
    if( stat == FORK_STATUS_RUNNING )then
      call proc%kill()
      call proc%await_final_status()
    endif
    call assert_int(stat, FORK_STATUS_STOPPED, 'child forked under a running monitor stops by itself')
    child_telemetry = 'memory_usage_'//int2str(pid)//'.csv'
    call assert_true(file_exists(child_telemetry), 'child writes its own memory telemetry')
    call assert_true(mem_monitor_is_enabled(),     'parent monitor still runs after the fork')
    call del_file(child_telemetry)
    call del_file(log_fname)
    if( l_own_monitor )then
      call mem_monitor_finish()
      call del_file('memory_usage_'//int2str(get_process_id())//'.csv')
    endif
    call cline%kill()
  end subroutine test_fork_with_running_monitor

end module simple_forked_process_tester
