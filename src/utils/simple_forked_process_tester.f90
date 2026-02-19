!@descr: various functions for testing unix forking utilities
module simple_forked_process_tester
use unix
use simple_forked_process
use simple_string
use simple_test_utils
use simple_syslib,         only: file_exists, del_file
implicit none

public :: run_all_forked_process_tests
private
#include "simple_local_flags.inc"

contains

  subroutine run_all_forked_process_tests()
    write(*,'(A)') '**** running all forked process tests ****'
    call test_start()
    call test_kill()
    call test_restart()
    call test_terminate()
    call test_logfile_redirection()
#ifdef __linux__
    ! message queues only work on linux
    call test_ipc_mq_communication()
#endif
  end subroutine run_all_forked_process_tests
    
  !---------------- basic lifecycle ----------------

  subroutine test_start()
    type(forked_process) :: fork_proc
    write(*,'(A)') 'test_start'
    call fork_proc%start(name=string('TEST_START'))
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_STOPPED, 'forked process terminated')
  end subroutine test_start
    
  subroutine test_kill()
    type(forked_process) :: fork_proc
    integer              :: rc
    write(*,'(A)') 'test_kill'
    call fork_proc%start(name=string('TEST_KILL'))
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call fork_proc%kill()
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_FAILED, 'forked process killed')
  end subroutine test_kill

  subroutine test_restart()
    type(forked_process)  :: fork_proc
    integer(kind=c_pid_t) :: pid1, pid2
    integer               :: rc
    write(*,'(A)') 'test_restart'
    call fork_proc%start(name=string('TEST_RESTART'), restart=.true.)
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    pid1 = fork_proc%get_pid()
    call assert_true(pid1 /= -1, 'forked process has valid pid')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call fork_proc%kill()
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_STOPPED, 'forked process stopped')
    pid2 = fork_proc%get_pid()
    call assert_true(pid2 /= -1, 'restarted forked process has valid pid')
    call assert_true(pid1 /= pid2, 'restarted forked process has different pid to original')
  end subroutine test_restart

  subroutine test_terminate()
    type(forked_process) :: fork_proc
    integer              :: rc
    write(*,'(A)') 'test_terminate'
    call fork_proc%start(name=string('TEST_TERMINATE'))
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    rc = c_usleep(FORK_POLL_TIME * 5)
    call fork_proc%terminate()
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_STOPPED, 'forked process terminated')
  end subroutine test_terminate

  !---------------- I/O ----------------

  subroutine test_logfile_redirection()
    type(forked_process) :: fork_proc
    type(string)         :: log_contents, log_hash, log_fname
    integer              :: rc
    write(*,'(A)') 'test_logfile_redirection'
    log_fname = 'test_logfile_redirection.log'
    call fork_proc%start(name=string('TEST_LOGFILE_REDIRECTION'), logfile=log_fname)
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_STOPPED, 'forked process terminated')
    call assert_true(file_exists(log_fname), 'log file exists')
    call del_file(log_fname)
  end subroutine test_logfile_redirection

  !---------------- IPC ----------------

  subroutine test_ipc_mq_communication()
    type(forked_process) :: fork_proc
    type(string)         :: mq_name, msg_out
    write(*,'(A)') 'test_ipc_mq_communication'
    call fork_proc%start(name=string('TEST_IPC_MQ'), mq_name=string('TEST_IPC_MQ'))
    call assert_int(fork_proc%status(), FORK_STATUS_RUNNING, 'forked process running')
    call fork_proc%await_final_status()
    call assert_int(fork_proc%status(), FORK_STATUS_STOPPED, 'forked process terminated')
    call assert_true(fork_proc%receive(msg=msg_out), 'received message from child')
    call fork_proc%destroy()
  end subroutine test_ipc_mq_communication

end module simple_forked_process_tester