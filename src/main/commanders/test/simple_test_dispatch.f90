!@descr: the test dispatcher: executions as forked children in their own sessions and directories, packed by CPU reservation, with deadlines
! Each child leads its own session in its own directory and runs the test program in the process.
! The processes of an execution are its session and every process carrying its marker
! SIMPLE_TEST_EXECUTION (the local queue starts workers in sessions of their own); deadlines,
! interrupts and the end of an execution signal all of them, and their CPU use and memory are
! sampled from /proc. The /proc-based parts (markers, samples, the processes left) are Linux-only;
! elsewhere the dispatcher says so and signals the session alone. The contract is in
! doc/policies/test_environment_policy.md, section 1.5.
module simple_test_dispatch
use, intrinsic :: iso_c_binding,   only: c_int, c_long, c_char, c_null_char, c_ptr, c_associated, c_funloc, &
                                         &c_funptr, c_null_funptr, c_size_t
use, intrinsic :: iso_fortran_env, only: int64
use unix,                  only: c_pid_t, c_usleep, c_signal, c_getpid, SIGINT, SIGTERM, SIGKILL, SIGHUP
use simple_defs,           only: logfhandle
use simple_error,          only: simple_exception
use simple_forked_process, only: forked_process, FORK_STATUS_RUNNING, sig_ign, is_sig_ign
use simple_cmdline,        only: cmdline
use simple_string,         only: string
implicit none

public :: test_dispatcher, test_execution, test_route, next_to_start, set_test_environment, proc_available
public :: EXEC_QUEUED, EXEC_RUNNING, EXEC_PASSED, EXEC_FAILED, EXEC_TIMEDOUT, EXEC_NOT_RUN, EXEC_INTERRUPTED
public :: exec_state_name
private
#include "simple_local_flags.inc"

integer, parameter :: EXEC_QUEUED = 0, EXEC_RUNNING = 1, EXEC_PASSED = 2, EXEC_FAILED = 3, EXEC_TIMEDOUT = 4
integer, parameter :: EXEC_NOT_RUN = 5 !< never started: the dispatcher was interrupted first
integer, parameter :: EXEC_INTERRUPTED = 6 !< stopped because the dispatcher was interrupted (no verdict)
integer, parameter :: POLL_US      = 200000 !< poll interval of the children (microseconds)
real,    parameter :: SAMPLE_S     = 1.     !< interval of the CPU and memory samples (seconds)
character(len=*), parameter :: MARKER = 'SIMPLE_TEST_EXECUTION' !< <dispatcher pid>.<execution index>

abstract interface
    !> runs the test program prg under cline, in the calling process
    subroutine test_route( prg, cline )
        import :: cmdline
        character(len=*), intent(in)    :: prg
        class(cmdline),   intent(inout) :: cline
    end subroutine test_route
end interface

!> one test execution: what to run and, once run, how it ended
type :: test_execution
    character(len=64)   :: label       = ''  !< names its directory and its row of the report
    character(len=32)   :: prg         = ''  !< simple_test_exec test=<prg>
    character(len=256)  :: args        = ''  !< the other keys of its command line
    character(len=1024) :: dir         = ''  !< its run directory, absolute
    integer :: seed_offset = 0
    integer :: nthr        = 1               !< OMP_NUM_THREADS of the child
    integer :: ncpu        = 1               !< peak CPU reservation
    real    :: mem_gb      = 0.              !< peak resident memory reservation (GiB)
    real    :: disk_gb     = 0.              !< what its run directory grows to (GiB)
    integer :: deadline    = 0               !< seconds; 0 is none
    integer :: state       = EXEC_QUEUED
    integer :: pid         = 0               !< of the child, which leads its session
    integer :: exit_code   = 0               !< exit status, or minus the signal that ended it
    real    :: wall        = 0.              !< seconds
    real    :: peak_cpu    = -1.             !< cores, the largest sampled use of its session
    real    :: peak_rss_gb = -1.             !< GiB, the largest sampled resident memory of its session
end type test_execution

!> a forked child that runs one execution
type, extends(forked_process) :: dispatch_child
    procedure(test_route), pointer, nopass :: route => null()
    character(len=32)   :: prg  = ''
    character(len=1024) :: dir  = ''
    character(len=32)   :: mark = ''  !< the value of its marker
    integer :: seed_offset = 0
    integer :: nthr        = 1
  contains
    procedure :: execute => execute_child
end type dispatch_child

type :: test_dispatcher
    private
    type(test_execution), allocatable :: execs(:)
    integer :: nexecs = 0
    integer :: budget = 1  !< cores
    integer :: grace  = 30 !< seconds between SIGTERM and SIGKILL
    logical :: keep_outputs = .true. !< .false.: a passed execution's subdirectories are deleted, its disk freed
  contains
    procedure :: new
    procedure :: add
    procedure :: run
    procedure :: get_nexecs
    procedure :: get_execution
    procedure :: processes_left
    procedure :: kill
end type test_dispatcher

logical, volatile :: l_interrupted = .false.

contains

    !> keep_outputs=.false.: every subdirectory of a passed execution's run directory (its outputs) is
    !! deleted when it ends, and its disk returns to the budget; run.log and metrics.tsv stay
    subroutine new( self, budget, grace, keep_outputs )
        class(test_dispatcher), intent(inout) :: self
        integer,                intent(in)    :: budget, grace
        logical, optional,      intent(in)    :: keep_outputs
        call self%kill
        self%budget = max(1, budget)
        self%grace  = max(0, grace)
        self%keep_outputs = .true.
        if( present(keep_outputs) ) self%keep_outputs = keep_outputs
        allocate(self%execs(16))
    end subroutine new

    !> queues one execution; dir is created when it starts
    subroutine add( self, exec )
        class(test_dispatcher), intent(inout) :: self
        type(test_execution),   intent(in)    :: exec
        type(test_execution), allocatable :: tmp(:)
        type(cmdline) :: cline
        ! a malformed command line stops here, before any child runs
        call command_line_of(exec, cline)
        call cline%kill
        if( self%nexecs == size(self%execs) )then
            allocate(tmp(2 * size(self%execs)))
            tmp(1:self%nexecs) = self%execs(1:self%nexecs)
            call move_alloc(tmp, self%execs)
        endif
        self%nexecs = self%nexecs + 1
        self%execs(self%nexecs)       = exec
        self%execs(self%nexecs)%state = EXEC_QUEUED
    end subroutine add

    integer function get_nexecs( self )
        class(test_dispatcher), intent(in) :: self
        get_nexecs = self%nexecs
    end function get_nexecs

    type(test_execution) function get_execution( self, i )
        class(test_dispatcher), intent(in) :: self
        integer,                intent(in) :: i
        get_execution = self%execs(i)
    end function get_execution

    !> how many processes of execution i are alive
    integer function processes_left( self, i )
        class(test_dispatcher), intent(in) :: self
        integer,                intent(in) :: i
        integer,        allocatable :: pids(:), sids(:), marks(:)
        integer(int64), allocatable :: ticks(:), rss(:)
        call scan_proc(pids, sids, ticks, rss, marks)
        processes_left = count(member(sids, marks, self%execs(i)%pid, i))
    end function processes_left

    subroutine kill( self )
        class(test_dispatcher), intent(inout) :: self
        if( allocated(self%execs) ) deallocate(self%execs)
        self%nexecs = 0
    end subroutine kill

    !> Runs every queued execution through route and returns when each has passed, failed or
    !! timed out. interrupted: SIGINT or SIGTERM stopped the dispatcher, which stopped every
    !! child group first; the executions never started are EXEC_NOT_RUN.
    subroutine run( self, route, interrupted )
        use simple_syslib, only: simple_mkdir
        class(test_dispatcher), intent(inout) :: self
        procedure(test_route)                 :: route
        logical,                intent(out)   :: interrupted
        type(dispatch_child), allocatable :: children(:)
        integer(int64), allocatable :: t_start(:), t_term(:)
        logical,        allocatable :: queued(:), l_term(:), l_kill(:), l_deadline(:), l_stopped(:), l_reaped_kill(:)
        integer,        allocatable :: ncpus(:), prev_pid(:)
        real,           allocatable :: mems(:), disks(:)
        integer(int64), allocatable :: prev_ticks(:)
        type(cmdline)  :: cline
        integer(int64) :: t_now, t_last_sample, rate
        integer        :: i, inext, nrunning, in_use, rc
        real           :: mem_budget, disk_budget, mem_in_use, disk_started
        type(c_funptr) :: prev_int, prev_term, prev_hup, prev
        l_interrupted = .false.
        interrupted   = .false.
        if( self%nexecs == 0 ) return
        ! a disposition inherited as ignored (nohup, a background job) stays ignored
        prev_int  = c_signal(SIGINT, sig_ign())
        if( .not. is_sig_ign(prev_int) ) prev = c_signal(SIGINT, c_funloc(on_interrupt))
        prev_hup  = c_signal(SIGHUP, sig_ign())
        if( .not. is_sig_ign(prev_hup) ) prev = c_signal(SIGHUP, c_funloc(on_interrupt))
        prev_term = c_signal(SIGTERM, c_funloc(on_interrupt))
        allocate(children(self%nexecs), t_start(self%nexecs), t_term(self%nexecs), queued(self%nexecs), &
            &l_term(self%nexecs), l_kill(self%nexecs), l_deadline(self%nexecs), l_stopped(self%nexecs), &
            &l_reaped_kill(self%nexecs), ncpus(self%nexecs), mems(self%nexecs), disks(self%nexecs))
        ncpus  = self%execs(1:self%nexecs)%ncpu
        mems   = self%execs(1:self%nexecs)%mem_gb
        disks  = self%execs(1:self%nexecs)%disk_gb
        queued = .true.
        l_term = .false.
        l_kill = .false.
        l_deadline    = .false.
        l_stopped     = .false.
        l_reaped_kill = .true.
        ! the memory free now and the disk of the run directories: what the executions may take
        call resource_budgets(self%execs(1)%dir, mem_budget, disk_budget)
        disk_started = 0.
        if( .not. proc_available() ) write(logfhandle,'(a)') '>>> DISPATCH: no /proc: the CPU and memory '//&
            &'samples, the markers of processes outside the session and the processes left are not available'
        call system_clock(t_last_sample, rate)
        ! the baseline every later sample of the CPU use differs from
        call sample_sessions(SAMPLE_S)
        do
            ! start what fits
            if( .not. l_interrupted )then
                do
                    nrunning   = count(self%execs(1:self%nexecs)%state == EXEC_RUNNING)
                    in_use     = sum(ncpus, mask=self%execs(1:self%nexecs)%state == EXEC_RUNNING)
                    mem_in_use = sum(mems,  mask=self%execs(1:self%nexecs)%state == EXEC_RUNNING)
                    inext      = next_to_start(ncpus, mems, disks, queued, self%budget - in_use, &
                        &mem_budget - mem_in_use, disk_budget - disk_started, nrunning)
                    if( inext == 0 ) exit
                    call start_execution(inext)
                enddo
            endif
            ! poll the running children
            call system_clock(t_now)
            do i = 1, self%nexecs
                if( self%execs(i)%state /= EXEC_RUNNING ) cycle
                ! an ended child is not reaped before its group is signalled: until then its pid and
                ! its group cannot belong to another process
                if( children(i)%has_exited() )then
                    call signal_execution(i, SIGKILL)
                    l_reaped_kill(i) = .false.
                endif
                if( children(i)%status() /= FORK_STATUS_RUNNING )then
                    call finish_execution(i)
                    cycle
                endif
                if( self%execs(i)%deadline > 0 .and. .not. l_term(i) )then
                    if( real(t_now - t_start(i)) / real(rate) > real(self%execs(i)%deadline) )then
                        write(logfhandle,'(a)') '>>> DISPATCH: '//trim(self%execs(i)%label)//' passed its deadline; SIGTERM'
                        call signal_execution(i, SIGTERM)
                        l_term(i)     = .true.
                        l_deadline(i) = .true.
                        t_term(i)     = t_now
                    endif
                endif
                if( l_term(i) .and. .not. l_kill(i) )then
                    if( real(t_now - t_term(i)) / real(rate) > real(self%grace) )then
                        write(logfhandle,'(a)') '>>> DISPATCH: '//trim(self%execs(i)%label)//' outlived the grace period; SIGKILL'
                        call signal_execution(i, SIGKILL)
                        l_kill(i) = .true.
                    endif
                endif
            enddo
            if( real(t_now - t_last_sample) / real(rate) >= SAMPLE_S )then
                call sample_sessions(real(t_now - t_last_sample) / real(rate))
                t_last_sample = t_now
            endif
            ! an interrupted dispatcher stops every child group, then waits for them
            if( l_interrupted .and. .not. interrupted )then
                interrupted = .true.
                write(logfhandle,'(a)') '>>> DISPATCH: interrupted; stopping every execution'
                do i = 1, self%nexecs
                    if( queued(i) )then
                        queued(i) = .false.
                        self%execs(i)%state = EXEC_NOT_RUN
                    endif
                    if( self%execs(i)%state == EXEC_RUNNING .and. .not. l_term(i) )then
                        call signal_execution(i, SIGTERM)
                        l_term(i)    = .true.
                        l_stopped(i) = .true.
                        t_term(i)    = t_now
                    endif
                enddo
            endif
            if( .not. any(queued) .and. .not. any(self%execs(1:self%nexecs)%state == EXEC_RUNNING) ) exit
            rc = c_usleep(POLL_US)
        enddo
        prev = c_signal(SIGINT,  prev_int)
        prev = c_signal(SIGTERM, prev_term)
        prev = c_signal(SIGHUP,  prev_hup)
        do i = 1, self%nexecs
            call children(i)%destroy
        enddo

      contains

        subroutine start_execution( iexec )
            use simple_syslib, only: del_file
            integer, intent(in) :: iexec
            queued(iexec) = .false.
            disk_started  = disk_started + disks(iexec)
            call simple_mkdir(trim(self%execs(iexec)%dir))
            ! a table left in a reused directory is not this execution's
            call del_file(trim(self%execs(iexec)%dir)//'/metrics.tsv')
            call command_line_of(self%execs(iexec), cline)
            children(iexec)%route       => route
            children(iexec)%prg         = self%execs(iexec)%prg
            children(iexec)%dir         = self%execs(iexec)%dir
            children(iexec)%seed_offset = self%execs(iexec)%seed_offset
            children(iexec)%nthr        = self%execs(iexec)%nthr
            write(children(iexec)%mark,'(i0,a,i0)') c_getpid(), '.', iexec
            write(logfhandle,'(a)') '>>> DISPATCH: start '//trim(self%execs(iexec)%label)
            call system_clock(t_start(iexec))
            call children(iexec)%start(name=string(trim(self%execs(iexec)%label)), cline=cline, new_session=.true.)
            self%execs(iexec)%pid   = int(children(iexec)%get_pid())
            self%execs(iexec)%state = EXEC_RUNNING
        end subroutine start_execution

        subroutine finish_execution( iexec )
            integer, intent(in) :: iexec
            integer(int64) :: t_end
            call system_clock(t_end)
            ! nothing the child started outlives it; once the child is reaped (no non-reaping wait
            ! outside Linux) its group may be another process's, so only its marked processes
            if( l_reaped_kill(iexec) ) call signal_execution(iexec, SIGKILL, l_group=.false.)
            self%execs(iexec)%wall      = real(t_end - t_start(iexec)) / real(rate)
            self%execs(iexec)%exit_code = children(iexec)%get_exit_code()
            if( l_deadline(iexec) )then
                self%execs(iexec)%state = EXEC_TIMEDOUT
            else if( self%execs(iexec)%exit_code == 0 )then
                self%execs(iexec)%state = EXEC_PASSED
            else if( l_stopped(iexec) )then
                self%execs(iexec)%state = EXEC_INTERRUPTED
            else
                self%execs(iexec)%state = EXEC_FAILED
            endif
            write(logfhandle,'(a,a,a,i0,a,f10.1,a)') '>>> DISPATCH: end ', trim(self%execs(iexec)%label), ': '//&
                &trim(exec_state_name(self%execs(iexec)%state))//', exit code ', self%execs(iexec)%exit_code, &
                &', ', self%execs(iexec)%wall, ' s'
            ! the outputs of a failed or timed-out execution stay, the evidence of what went wrong
            if( .not. self%keep_outputs .and. self%execs(iexec)%state == EXEC_PASSED ) call delete_outputs(iexec)
        end subroutine finish_execution

        !> the subdirectories of an ended execution's run directory; its share of the disk budget returns
        !! when every one of them is gone
        subroutine delete_outputs( iexec )
            use simple_syslib, only: simple_list_dirs, simple_rmdir
            integer, intent(in) :: iexec
            type(string), allocatable :: subdirs(:)
            integer :: j, status
            logical :: l_all
            subdirs = simple_list_dirs(trim(self%execs(iexec)%dir))
            l_all   = .true.
            do j = 1, size(subdirs)
                call simple_rmdir(trim(self%execs(iexec)%dir)//'/'//subdirs(j)%to_char(), status)
                l_all = l_all .and. status == 0
            enddo
            if( l_all ) disk_started = disk_started - disks(iexec)
        end subroutine delete_outputs

        !> sig to every process of execution iexec, its group first (which it may still be making)
        !! unless l_group is .false.
        subroutine signal_execution( iexec, sig, l_group )
            use unix, only: c_kill
            integer,           intent(in) :: iexec
            integer(c_int),    intent(in) :: sig
            logical, optional, intent(in) :: l_group
            integer,        allocatable :: pids(:), sids(:), marks(:)
            integer(int64), allocatable :: ticks(:), rss(:)
            integer :: ip
            logical :: l_grp
            l_grp = .true.
            if( present(l_group) ) l_grp = l_group
            if( l_grp ) call children(iexec)%signal_group(sig)
            call scan_proc(pids, sids, ticks, rss, marks)
            do ip = 1, size(pids)
                if( member(sids(ip), marks(ip), self%execs(iexec)%pid, iexec) ) rc = c_kill(int(pids(ip), c_pid_t), sig)
            enddo
        end subroutine signal_execution

        !> the CPU use (cores over the last interval) and resident memory of each running execution
        subroutine sample_sessions( dt )
            real, intent(in) :: dt
            integer, allocatable :: pids(:), sids(:), marks(:)
            integer(int64), allocatable :: ticks(:), rss(:)
            integer(int64) :: dticks, mem
            integer :: iexec, ip, j
            call scan_proc(pids, sids, ticks, rss, marks)
            if( size(pids) == 0 ) return
            do iexec = 1, self%nexecs
                if( self%execs(iexec)%state /= EXEC_RUNNING ) cycle
                dticks = 0_int64
                mem    = 0_int64
                do ip = 1, size(pids)
                    if( .not. member(sids(ip), marks(ip), self%execs(iexec)%pid, iexec) ) cycle
                    mem = mem + rss(ip)
                    j = 0
                    if( allocated(prev_pid) ) j = findloc(prev_pid, pids(ip), dim=1)
                    if( j > 0 )then
                        dticks = dticks + max(0_int64, ticks(ip) - prev_ticks(j))
                    else
                        dticks = dticks + ticks(ip)
                    endif
                enddo
                ! without a baseline (no /proc) the CPU use is not measured
                if( allocated(prev_pid) ) self%execs(iexec)%peak_cpu = max(self%execs(iexec)%peak_cpu, &
                    &real(dticks) / real(clock_ticks()) / dt)
                self%execs(iexec)%peak_rss_gb = max(self%execs(iexec)%peak_rss_gb, real(mem) / 1024.**3)
            enddo
            prev_pid   = pids
            prev_ticks = ticks
        end subroutine sample_sessions

    end subroutine run

    !> the command line a child runs: test= and prg= the program, then the key=value arguments
    subroutine command_line_of( exec, cline )
        use simple_string_utils, only: parsestr
        type(test_execution), intent(in)    :: exec
        type(cmdline),        intent(inout) :: cline
        character(len=256) :: tokens(64), args
        integer :: ntok, itok, pos
        call cline%kill
        call cline%set('test', trim(exec%prg))
        call cline%set('prg',  trim(exec%prg))
        if( len_trim(exec%args) == 0 ) return
        args = exec%args
        call parsestr(args, ' ', tokens, ntok)
        do itok = 1, ntok
            pos = index(tokens(itok), '=')
            if( pos < 2 ) THROW_HARD('not a key=value argument: '//trim(tokens(itok)))
            call cline%set_from_text(tokens(itok)(:pos-1), trim(tokens(itok)(pos+1:)))
        enddo
    end subroutine command_line_of

    !> a process of execution i: in the session its child leads, or carrying its marker
    elemental logical function member( sid, mark, pid, i )
        integer, intent(in) :: sid, mark, pid, i
        member = (pid > 0 .and. sid == pid) .or. mark == i
    end function member

    !> The next queued execution to start: the largest reservation queued (the first in queue order
    !! among equals) when its cores, memory and disk fit in what is free, or whenever nothing runs
    !! (a reservation above a budget runs alone); 0 when it may not start now. Smaller executions
    !! wait for every larger one to start, so they cannot starve it.
    pure integer function next_to_start( ncpu, mem, disk, queued, free_cpu, free_mem, free_disk, nrunning )
        integer, intent(in) :: ncpu(:), free_cpu, nrunning
        real,    intent(in) :: mem(:), disk(:), free_mem, free_disk
        logical, intent(in) :: queued(:)
        integer :: l
        next_to_start = 0
        if( .not. any(queued) ) return
        l = maxloc(ncpu, mask=queued, dim=1)
        if( nrunning == 0 )then
            next_to_start = l
        else if( ncpu(l) <= free_cpu .and. mem(l) <= free_mem .and. disk(l) <= free_disk )then
            next_to_start = l
        endif
    end function next_to_start

    !> what the executions may take together: 90 % of the memory available now (MemAvailable of
    !! /proc/meminfo) and the space free on the file system of their directories (GiB); a budget
    !! that cannot be read is unlimited
    subroutine resource_budgets( dir, mem_budget, disk_budget )
        use simple_syslib, only: fs_avail_bytes
        character(len=*), intent(in)  :: dir
        real,             intent(out) :: mem_budget, disk_budget
        character(len=256) :: line
        integer(int64)     :: kib, avail
        integer :: funit, ios, slash
        mem_budget  = huge(1.)
        disk_budget = huge(1.)
        open(newunit=funit, file='/proc/meminfo', action='read', status='old', iostat=ios)
        if( ios == 0 )then
            do
                read(funit,'(a)',iostat=ios) line
                if( ios /= 0 ) exit
                if( index(line, 'MemAvailable:') == 1 )then
                    read(line(14:),*,iostat=ios) kib
                    if( ios == 0 ) mem_budget = 0.9 * real(kib) / 1024.**2
                    exit
                endif
            enddo
            close(funit)
        endif
        ! the directory may not exist yet: its parent does
        slash = index(trim(dir), '/', back=.true.)
        if( slash > 1 )then
            avail = fs_avail_bytes(dir(:slash-1))
            if( avail > 0_int64 ) disk_budget = real(avail) / 1024.**3
        endif
    end subroutine resource_budgets

    !> whether /proc, which the markers, samples and the processes left are read from, exists
    logical function proc_available()
        inquire(file='/proc/self/stat', exist=proc_available)
    end function proc_available

    function exec_state_name( state ) result( name )
        integer, intent(in) :: state
        character(len=:), allocatable :: name
        select case(state)
            case(EXEC_QUEUED);   name = 'queued'
            case(EXEC_RUNNING);  name = 'running'
            case(EXEC_PASSED);   name = 'passed'
            case(EXEC_FAILED);   name = 'failed'
            case(EXEC_TIMEDOUT); name = 'timed_out'
            case(EXEC_NOT_RUN);  name = 'not_run'
            case(EXEC_INTERRUPTED); name = 'interrupted'
            case default;        name = 'unknown'
        end select
    end function exec_state_name

    !> In the forked child: lead the session (done by start), die with the dispatcher, enter the
    !! run directory, write stdout and stderr to run.log there, end at once on SIGTERM (the gate
    !! table and the log are flushed row by row), set the threads and the seed offset, run the test.
    subroutine execute_child( self, cline )
        !$ use omp_lib, only: omp_set_num_threads
        use unix,              only: c_open, c_dup2, c_close, O_WRONLY, O_CREAT, O_TRUNC, c_mode_t
        use simple_syslib,     only: simple_chdir
        use simple_test_utils, only: set_test_seed_offset, fix_test_seed
        class(dispatch_child), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        character(len=16) :: str
        type(c_funptr)    :: prev
        integer(c_int)    :: fd, rc
        integer           :: status
#if defined(__linux__)
        interface
            function c_prctl( option, arg2 ) bind(c, name='prctl') result( ierr )
                import :: c_int, c_long
                integer(c_int),  value :: option
                integer(c_long), value :: arg2
                integer(c_int)         :: ierr
            end function c_prctl
        end interface
        integer(c_int), parameter :: PR_SET_PDEATHSIG = 1
        rc = c_prctl(PR_SET_PDEATHSIG, int(SIGKILL, c_long))
#endif
        prev = c_signal(SIGTERM, c_funloc(on_child_term))
        call simple_chdir(trim(self%dir), status)
        if( status /= 0 ) error stop 'test dispatcher: cannot enter the run directory'
        fd = c_open(trim(self%dir)//'/run.log'//c_null_char, ior(O_WRONLY, ior(O_CREAT, O_TRUNC)), &
            &int(o'644', c_mode_t))
        if( fd < 0 ) error stop 'test dispatcher: cannot open run.log'
        rc = c_dup2(fd, 1_c_int)
        rc = c_dup2(fd, 2_c_int)
        rc = c_close(fd)
        write(str,'(i0)') self%nthr
        call put_env('OMP_NUM_THREADS', trim(str))
        call put_env(MARKER, trim(self%mark))
        !$ call omp_set_num_threads(self%nthr)
        call set_test_seed_offset(self%seed_offset)
        call fix_test_seed
        call self%route(trim(self%prg), cline)
    end subroutine execute_child

    !> every process of /proc: pid, session, CPU ticks (user and system), resident bytes, and the
    !! execution index its marker names when this process dispatched it (0 otherwise)
    subroutine scan_proc( pids, sids, ticks, rss, marks )
        use unix, only: c_opendir, c_closedir, f_readdir, c_dirent
        integer,        allocatable, intent(out) :: pids(:), sids(:), marks(:)
        integer(int64), allocatable, intent(out) :: ticks(:), rss(:)
        type(c_dirent), pointer :: entry
        type(c_ptr)       :: dirp
        character(len=32) :: dname
        integer           :: pid, sid, ios, n, i, nmax, self_pid
        integer(int64)    :: tck, mem
        nmax = 4096
        allocate(pids(nmax), sids(nmax), ticks(nmax), rss(nmax), marks(nmax))
        n        = 0
        self_pid = int(c_getpid())
        dirp     = c_opendir('/proc'//c_null_char)
        if( c_associated(dirp) )then
            do
                entry => f_readdir(dirp)
                if( .not. associated(entry) ) exit
                dname = ''
                do i = 0, min(31, size(entry%d_name) - 1)
                    if( entry%d_name(i) == c_null_char ) exit
                    dname(i+1:i+1) = entry%d_name(i)
                enddo
                if( verify(trim(dname), '0123456789') /= 0 .or. len_trim(dname) == 0 ) cycle
                read(dname,*,iostat=ios) pid
                if( ios /= 0 ) cycle
                if( .not. read_stat(pid, sid, tck, mem) ) cycle
                if( n == nmax )then
                    nmax  = 2 * nmax
                    pids  = [pids,  pids]
                    sids  = [sids,  sids]
                    ticks = [ticks, ticks]
                    rss   = [rss,   rss]
                    marks = [marks, marks]
                endif
                n = n + 1
                pids(n)  = pid
                sids(n)  = sid
                ticks(n) = tck
                rss(n)   = mem
                marks(n) = read_mark(pid, self_pid)
            enddo
            ios = c_closedir(dirp)
        endif
        pids  = pids(1:n)
        sids  = sids(1:n)
        ticks = ticks(1:n)
        rss   = rss(1:n)
        marks = marks(1:n)
    end subroutine scan_proc

    !> the execution index in the marker of process pid when dispatcher_pid set it, else 0; the
    !! environment of a process is NUL-separated, and only another user's is unreadable
    integer function read_mark( pid, dispatcher_pid )
        use, intrinsic :: iso_c_binding, only: c_loc
        use unix, only: c_open, c_read, c_close, O_RDONLY, c_mode_t
        integer, intent(in) :: pid, dispatcher_pid
        character(kind=c_char, len=4096), target :: buf
        character(len=:), allocatable :: env
        character(len=48)  :: fname, prefix
        integer(c_int)     :: fd, rc
        integer(c_size_t)  :: n
        integer :: ios, pos, pend, idx
        read_mark = 0
        write(fname,'(a,i0,a)') '/proc/', pid, '/environ'
        write(prefix,'(a,i0,a)') MARKER//'=', dispatcher_pid, '.'
        fd = c_open(trim(fname)//c_null_char, O_RDONLY, 0_c_mode_t)
        if( fd < 0 ) return
        ! to its end, however long: a /proc file reports no size
        env = ''
        do
            n = c_read(fd, c_loc(buf), int(len(buf), c_size_t))
            if( n <= 0_c_size_t .or. n > int(len(buf), c_size_t) ) exit
            env = env//buf(1:int(n))
        enddo
        rc = c_close(fd)
        ! the variable starts the environment or follows a NUL
        pos = index(achar(0)//env, achar(0)//trim(prefix))
        if( pos == 0 ) return
        pos  = pos + len_trim(prefix) - 1
        pend = index(env(pos+1:), achar(0))
        if( pend == 0 ) pend = len(env) - pos + 1
        read(env(pos+1:pos+pend-1),*,iostat=ios) idx
        if( ios == 0 ) read_mark = idx
    end function read_mark

    !> session, CPU ticks and resident bytes of one process, from /proc/<pid>/stat
    logical function read_stat( pid, sid, tck, mem )
        integer,        intent(in)  :: pid
        integer,        intent(out) :: sid
        integer(int64), intent(out) :: tck, mem
        character(len=1024) :: line
        character(len=32)   :: fname
        character(len=1)    :: state
        integer(int64)      :: f(21)
        integer :: funit, ios, pos
        read_stat = .false.
        sid = 0
        tck = 0_int64
        mem = 0_int64
        write(fname,'(a,i0,a)') '/proc/', pid, '/stat'
        open(newunit=funit, file=trim(fname), action='read', status='old', iostat=ios)
        if( ios /= 0 ) return
        read(funit,'(a)',iostat=ios) line
        close(funit)
        if( ios /= 0 ) return
        ! the command name, in parentheses, may hold blanks: the fields start after its last ')'
        pos = index(line, ')', back=.true.)
        if( pos == 0 ) return
        ! fields 3 to 24: state, ppid, pgrp, session, ..., utime (14), stime (15), ..., rss (24); not
        ! further, rsslim (25) may exceed a signed 64-bit integer
        read(line(pos+1:),*,iostat=ios) state, f
        if( ios /= 0 ) return
        ! a zombie has ended; it only waits to be reaped
        if( state == 'Z' .or. state == 'X' ) return
        sid = int(f(3))
        tck = f(11) + f(12)
        mem = f(21) * page_bytes()
        read_stat = .true.
    end function read_stat

    integer(int64) function clock_ticks()
        interface
            function c_sysconf( name ) bind(c, name='sysconf')
                import :: c_int, c_long
                integer(c_int), value :: name
                integer(c_long)       :: c_sysconf
            end function c_sysconf
        end interface
        integer(c_int), parameter :: SC_CLK_TCK = 2 ! Linux
        clock_ticks = max(1_int64, int(c_sysconf(SC_CLK_TCK), int64))
    end function clock_ticks

    integer(int64) function page_bytes()
        interface
            function c_sysconf( name ) bind(c, name='sysconf')
                import :: c_int, c_long
                integer(c_int), value :: name
                integer(c_long)       :: c_sysconf
            end function c_sysconf
        end interface
        integer(c_int), parameter :: SC_PAGESIZE = 30 ! Linux
        page_bytes = max(1_int64, int(c_sysconf(SC_PAGESIZE), int64))
    end function page_bytes

    !> The environment every child and the workers it starts inherit, as CTest gives it: SIMPLE_PATH
    !! the build this executable belongs to (the parent of its bin/, read from /proc/self/exe; else
    !! the SIMPLE_PATH it was given), its bin/ and scripts/ first on PATH, the local queue.
    subroutine set_test_environment( simple_path )
        character(len=:), allocatable, intent(out) :: simple_path
        character(len=4096) :: exe, path
        integer :: n, status
        logical :: l_exists
        simple_path = ''
        inquire(file='/proc/self/exe', exist=l_exists)
        if( l_exists )then
            n = readlink_self(exe)
            if( n > 4 )then
                ! .../bin/simple_test_exec
                exe = exe(:index(exe(:n), '/', back=.true.) - 1)
                if( exe(max(1,len_trim(exe)-3):len_trim(exe)) == '/bin' ) simple_path = exe(:len_trim(exe)-4)
            endif
        endif
        if( len(simple_path) == 0 )then
            call get_environment_variable('SIMPLE_PATH', path, status=status)
            if( status /= 0 ) THROW_HARD('the test dispatcher needs SIMPLE_PATH or an installed simple_test_exec')
            simple_path = trim(path)
        endif
        inquire(file=simple_path//'/bin/simple_private_exec', exist=l_exists)
        if( .not. l_exists ) THROW_HARD('no simple_private_exec in '//simple_path//'/bin')
        call get_environment_variable('PATH', path)
        call put_env('SIMPLE_PATH', simple_path)
        call put_env('PATH', simple_path//'/bin:'//simple_path//'/scripts:'//trim(path))
        call put_env('SIMPLE_QSYS', 'local')
        call put_env('SIMPLE_EMAIL', 'test@localhost')
        ! for the programs the children start: the dispatcher's runtime, which the forked children
        ! share, read its environment before this point
        call put_env('GFORTRAN_UNBUFFERED_ALL', 'y')
    end subroutine set_test_environment

    !> the target of /proc/self/exe in path; its length, or -1
    integer function readlink_self( path )
        character(len=*), intent(out) :: path
        interface
            function c_readlink( pathname, buf, bufsiz ) bind(c, name='readlink') result( n )
                import :: c_char, c_long, c_size_t
                character(kind=c_char), intent(in)  :: pathname(*)
                character(kind=c_char), intent(out) :: buf(*)
                integer(c_size_t), value            :: bufsiz
                integer(c_long)                     :: n
            end function c_readlink
        end interface
        character(kind=c_char) :: buf(len(path))
        integer :: i
        path = ''
        readlink_self = int(c_readlink('/proc/self/exe'//c_null_char, buf, int(len(path), c_size_t)))
        do i = 1, max(0, readlink_self)
            path(i:i) = buf(i)
        enddo
    end function readlink_self

    subroutine put_env( name, val )
        character(len=*), intent(in) :: name, val
        interface
            function c_setenv( name, value, overwrite ) bind(C, name='setenv') result( ierr )
                import :: c_char, c_int
                character(kind=c_char), intent(in) :: name(*), value(*)
                integer(c_int), value              :: overwrite
                integer(c_int)                     :: ierr
            end function c_setenv
        end interface
        if( c_setenv(trim(name)//c_null_char, trim(val)//c_null_char, 1_c_int) /= 0 )&
            &error stop 'test dispatcher: setenv failed'
    end subroutine put_env

    !> a child past its deadline ends at once, with the status of the signal: _exit is
    !! async-signal-safe (a flush or the Fortran exit is not); its gate rows and log lines were
    !! flushed as they were written
    subroutine on_child_term( signum ) bind(c)
        integer(c_int), value :: signum
        interface
            subroutine c_exit_now( status ) bind(c, name='_exit')
                import :: c_int
                integer(c_int), value :: status
            end subroutine c_exit_now
        end interface
        call c_exit_now(128_c_int + signum)
    end subroutine on_child_term

    subroutine on_interrupt( signum ) bind(c)
        integer(c_int), value :: signum
        l_interrupted = signum > 0
    end subroutine on_interrupt

end module simple_test_dispatch
