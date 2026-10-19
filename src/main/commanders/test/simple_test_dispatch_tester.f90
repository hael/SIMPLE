!@descr: platform tests of the test dispatcher (simple_test_dispatch): failing, hanging and worker-hanging children beside passing ones
! Probes that pass, fail, hang, ignore SIGTERM, wait on a worker hanging in its own session, or
! leave a worker behind (in their group, or in a session of their own) must each end as they
! should, within the deadlines and grace periods, with no process of theirs left (none in the
! group, none carrying the marker). An interrupted dispatcher must stop every execution and report
! the rest not run.
module simple_test_dispatch_tester
use, intrinsic :: iso_c_binding,   only: c_int
use, intrinsic :: iso_fortran_env, only: int64
use unix,                 only: c_kill, c_getpid, c_pid_t, c_usleep, c_signal, c_funptr, SIGINT, SIGTERM
use simple_cmdline,       only: cmdline
use simple_syslib,        only: simple_mkdir, simple_rmdir, simple_getcwd
use simple_string,        only: string
use simple_test_utils,    only: assert_true, assert_int
use simple_test_dispatch, only: test_dispatcher, test_execution, EXEC_PASSED, EXEC_FAILED, EXEC_TIMEDOUT, EXEC_NOT_RUN, &
    &EXEC_INTERRUPTED, proc_available
implicit none

public :: run_all_test_dispatch_tests
private

integer, parameter :: DEADLINE = 2 !< seconds of the hanging probes
integer, parameter :: GRACE    = 1 !< seconds between SIGTERM and SIGKILL
integer(c_pid_t)   :: dispatcher_pid = 0

contains

    subroutine run_all_test_dispatch_tests()
#if !defined(_WIN32)
        write(*,'(A)') '**** running all test dispatch tests ****'
        ! the probes that leave processes outside the session are found and stopped through /proc
        if( .not. proc_available() )then
            write(*,'(A)') 'SKIP: test dispatch tests need /proc (Linux)'
            return
        endif
        call test_outcomes()
        call test_interrupt()
#endif
    end subroutine run_all_test_dispatch_tests

    !> every probe ends as it should, in time, and leaves no process behind
    subroutine test_outcomes()
        character(len=*), parameter :: ROOT = 'tmp_test_dispatch'
        character(len=16), parameter :: PROBES(7) = [character(len=16) :: 'pass', 'fail', 'hang', &
            &'stubborn', 'worker_hang', 'orphan', 'escape']
        integer, parameter :: EXPECTED(7) = [EXEC_PASSED, EXEC_FAILED, EXEC_TIMEDOUT, EXEC_TIMEDOUT, EXEC_TIMEDOUT, &
            &EXEC_PASSED, EXEC_PASSED]
        type(test_dispatcher) :: dispatcher
        type(test_execution)  :: exec
        type(string) :: cwd
        integer(int64) :: t0, t1, rate
        integer      :: i
        logical      :: interrupted
        write(*,'(A)') 'test_outcomes'
        call simple_getcwd(cwd)
        call simple_mkdir(ROOT)
        ! a budget of two cores, one per probe: they run two by two, the hanging ones last
        call dispatcher%new(2, GRACE)
        do i = 1, size(PROBES)
            exec          = test_execution()
            exec%label    = PROBES(i)
            exec%prg      = PROBES(i)
            exec%dir      = cwd%to_char()//'/'//ROOT//'/'//trim(PROBES(i))
            exec%ncpu     = 1
            exec%deadline = DEADLINE
            call dispatcher%add(exec)
        enddo
        call system_clock(t0, rate)
        call dispatcher%run(probe_route, interrupted)
        call system_clock(t1)
        call assert_true(.not. interrupted, 'the dispatcher ran to its end')
        do i = 1, size(PROBES)
            exec = dispatcher%get_execution(i)
            call assert_int(EXPECTED(i), exec%state, 'the outcome of the probe '//trim(PROBES(i)))
            call assert_true(group_gone(exec%dir), 'no process of the group of '//trim(PROBES(i))//' is left')
            call assert_int(0, left_after_wait(dispatcher, i), 'no process of '//trim(PROBES(i))//' is left')
        enddo
        exec = dispatcher%get_execution(2)
        call assert_int(3, exec%exit_code, 'a failed child keeps its exit status')
        exec = dispatcher%get_execution(4)
        call assert_int(-9, exec%exit_code, 'a child that ignores SIGTERM ends by SIGKILL')
        ! four pairs at most: the hanging probes need the deadline and the grace period at most once each
        call assert_true(real(t1 - t0) / real(rate) < real(4 * (DEADLINE + GRACE) + 3), &
            &'the run ends within its deadlines and grace periods')
        call dispatcher%kill
        call simple_rmdir(ROOT)
    end subroutine test_outcomes

    !> an interrupted dispatcher stops every child group and starts nothing more
    subroutine test_interrupt()
        character(len=*), parameter :: ROOT = 'tmp_test_dispatch_interrupt'
        character(len=16), parameter :: PROBES(3) = [character(len=16) :: 'worker_hang', 'interrupt', 'pass']
        type(test_dispatcher) :: dispatcher
        type(test_execution)  :: exec
        type(string) :: cwd
        integer      :: i
        logical      :: interrupted
        write(*,'(A)') 'test_interrupt'
        call simple_getcwd(cwd)
        call simple_mkdir(ROOT)
        dispatcher_pid = c_getpid()
        call dispatcher%new(2, GRACE)
        do i = 1, size(PROBES)
            exec          = test_execution()
            exec%label    = PROBES(i)
            exec%prg      = PROBES(i)
            exec%dir      = cwd%to_char()//'/'//ROOT//'/'//trim(PROBES(i))
            exec%ncpu     = 1
            exec%deadline = 60
            call dispatcher%add(exec)
        enddo
        call dispatcher%run(probe_route, interrupted)
        call assert_true(interrupted, 'the dispatcher reports the interrupt')
        exec = dispatcher%get_execution(1)
        call assert_int(EXEC_INTERRUPTED, exec%state, 'a child stopped by the interrupt is interrupted, not failed')
        call assert_true(group_gone(exec%dir), 'and its group is gone')
        call assert_int(0, left_after_wait(dispatcher, 1), 'and so is its worker')
        exec = dispatcher%get_execution(3)
        call assert_int(EXEC_NOT_RUN, exec%state, 'an execution not yet started is not run')
        call dispatcher%kill
        call simple_rmdir(ROOT)
    end subroutine test_interrupt

    !> the probes, run in the forked child (its run directory is the working directory)
    subroutine probe_route( prg, cline )
        character(len=*), intent(in)    :: prg
        class(cmdline),   intent(inout) :: cline
        type(c_funptr)  :: prev
        integer(c_int)  :: rc
        integer         :: funit
        open(newunit=funit, file='pgid.txt', action='write', status='replace')
        write(funit,'(i0)') c_getpid()
        close(funit)
        select case(prg)
            case('pass')
                write(*,'(a)') 'probe passes'
            case('fail')
                error stop 3
            case('hang')
                call hang
            case('stubborn')
                prev = c_signal(SIGTERM, transfer(1_8, prev)) ! SIG_IGN
                call hang
            case('worker_hang')
                call execute_command_line('setsid sleep 600')
            case('orphan')
                call execute_command_line('sleep 600 &')
            case('escape')
                call execute_command_line('setsid sleep 600 &')
            case('interrupt')
                rc = c_kill(dispatcher_pid, SIGINT)
                call hang
        end select
    end subroutine probe_route

    subroutine hang
        integer(c_int) :: rc
        do
            rc = c_usleep(100000)
        enddo
    end subroutine hang

    !> the processes of execution i still alive, allowing those just killed two seconds to end
    integer function left_after_wait( dispatcher, i )
        class(test_dispatcher), intent(in) :: dispatcher
        integer,                intent(in) :: i
        integer(c_int) :: rc
        integer :: itry
        do itry = 1, 20
            left_after_wait = dispatcher%processes_left(i)
            if( left_after_wait == 0 ) return
            rc = c_usleep(100000)
        enddo
    end function left_after_wait

    !> no process is left in the group the probe of dir led (its pid, written to pgid.txt)
    logical function group_gone( dir )
        character(len=*), intent(in) :: dir
        integer(c_int) :: rc
        integer :: funit, ios, pgid, itry
        group_gone = .false.
        open(newunit=funit, file=trim(dir)//'/pgid.txt', action='read', status='old', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) pgid
        close(funit)
        if( ios /= 0 ) return
        ! signal 0 tests for existence: it fails when the group has no process left; a process
        ! killed a moment ago may still wait to be reaped, so allow it two seconds
        do itry = 1, 20
            group_gone = c_kill(-int(pgid, c_pid_t), 0_c_int) /= 0
            if( group_gone ) return
            rc = c_usleep(100000)
        enddo
    end function group_gone

end module simple_test_dispatch_tester
