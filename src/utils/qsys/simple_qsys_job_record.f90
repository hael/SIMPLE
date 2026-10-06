!@descr: the record a queued job's script keeps of itself, to cancel the job, and directories set aside when a job may still use them
!==============================================================================
! MODULE: simple_qsys_job_record
!
! PURPOSE:
!   An asynchronous job's script (simple_qsys_ctrl, generate_script_2 with an
!   exit-status file) records the job when it starts, in
!   <exit-status file>.job: its pid, its host, and its scheduler and job id
!   ("slurm <id>", "lsf <id>", "pbs <id>" or "none 0"). cancel_queued_job uses
!   that record: scancel, bkill or qdel for a scheduler job; for a job on this
!   host, a SIGTERM to its process group when the script leads one (a local job
!   is started under setsid, simple_qsys_local), which reaches the part jobs a
!   distributed program starts, otherwise to the script and its children (a
!   persistent-worker task). A job still queued (no record yet), one on another
!   host without a scheduler id, and the jobs a program queues itself are not
!   cancelled.
!
!   Liveness (follow-up plan, decision 29): query_job asks the scheduler, or the
!   host, whether a recorded job without an exit status still exists, and
!   check_job_lost does so every JOB_LIVENESS_S; a job gone at two checks in a
!   row (one miss can be a scheduler's completing state) has its exit status
!   written as JOB_LOST_EXIT_CODE, so every caller takes its failure path.
!
!   cancel_unfinished_jobs cancels every such job recorded under a directory,
!   for a process that ended without cancelling its own (a stream stage the
!   master had to kill).
!
!   A directory whose job recorded itself and wrote no exit status may still
!   be in use: fresh_job_dir and set_aside_dir move such a directory aside to
!   <dir>_unfinished<k> before new work starts there. A job that keeps running
!   writes into the moved directory.
!
!   Kept apart from simple_qsys_async_job, which builds on simple_qsys_env, so
!   that the queue controller can use it too.
!==============================================================================
module simple_qsys_job_record
use simple_defs,         only: logfhandle, XLONGSTRLEN
use simple_timer,        only: simple_gettime
use simple_defs_fname,   only: JOB_INFO_EXT
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: file_exists
use simple_syslib,       only: simple_mkdir, dir_exists, simple_rename, exec_cmdline, get_process_id
implicit none

public :: cancel_queued_job, cancel_unfinished_jobs, job_left_unfinished, fresh_job_dir, set_aside_dir
public :: query_job, check_job_lost, JOB_ALIVE, JOB_GONE, JOB_UNKNOWN, JOB_LIVENESS_S, JOB_LOST_EXIT_CODE
private

integer, parameter :: JOB_UNKNOWN        = 0   ! no record yet, a status written, or a local job of another host
integer, parameter :: JOB_ALIVE          = 1
integer, parameter :: JOB_GONE           = 2   ! the scheduler, or the host, no longer has it
integer, parameter :: JOB_LIVENESS_S     = 300 ! between two liveness checks of a job without an exit status
integer, parameter :: JOB_LOST_EXIT_CODE = 254 ! the status written for a job that vanished without one

contains

    !> Cancels the job that recorded itself next to @p exit_code_fname; .true. when a cancel was
    !! sent. A job that has written its exit status, or has not recorded itself yet, is left alone.
    logical function cancel_queued_job( exit_code_fname ) result( l_sent )
        class(string), intent(in) :: exit_code_fname
        character(len=256) :: host, line, sched, schedid
        type(string)       :: info_fname, cmd
        integer            :: funit, ios, pid
        l_sent = .false.
        if( file_exists(exit_code_fname) ) return
        info_fname = exit_code_fname//JOB_INFO_EXT
        if( .not. file_exists(info_fname) ) return
        open(newunit=funit, file=info_fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) pid
        if( ios == 0 ) read(funit,'(A)',iostat=ios) host
        if( ios == 0 ) read(funit,'(A)',iostat=ios) line
        close(funit)
        if( ios /= 0 ) return
        read(line,*,iostat=ios) sched, schedid
        if( ios /= 0 ) return
        select case(trim(sched))
            case('slurm')
                cmd = 'scancel '//trim(schedid)
            case('lsf')
                cmd = 'bkill '//trim(schedid)
            case('pbs')
                cmd = 'qdel '//trim(schedid)
            case DEFAULT
                ! a pid means nothing on another host; the process group when the script leads one
                cmd = 'if [ "$(hostname)" = "'//trim(host)//'" ]; then kill -TERM -- -'//int2str(pid)//' 2> /dev/null || '//&
                    &'{ pkill -TERM -P '//int2str(pid)//'; kill -TERM '//int2str(pid)//'; }; fi'
        end select
        call exec_cmdline(cmd%to_char()//' > /dev/null 2>&1', suppress_errors=.true.)
        l_sent = .true.
    end function cancel_queued_job

    !> Cancels every job recorded under @p dir, at any depth, that wrote no exit status
    !! (cancel_queued_job); @p ncancelled is the number of cancels sent.
    subroutine cancel_unfinished_jobs( dir, ncancelled )
        class(string), intent(in)  :: dir
        integer,       intent(out) :: ncancelled
        character(len=XLONGSTRLEN) :: line
        type(string) :: listfile, exit_code_fname
        integer      :: funit, ios, n
        ncancelled = 0
        if( .not. dir_exists(dir) ) return
        listfile = '__simple_job_records_'//int2str(get_process_id())//'__'
        call exec_cmdline('find "'//dir%to_char()//'" -type f -name "*'//JOB_INFO_EXT//'" > '//listfile%to_char()//&
            &' 2> /dev/null', suppress_errors=.true.)
        if( .not. file_exists(listfile) ) return
        open(newunit=funit, file=listfile%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            n = len_trim(line) - len(JOB_INFO_EXT)
            if( n < 1 ) cycle
            exit_code_fname = line(1:n)
            if( cancel_queued_job(exit_code_fname) ) ncancelled = ncancelled + 1
        enddo
        close(funit, status='delete')
    end subroutine cancel_unfinished_jobs

    !> Whether the job that recorded itself next to @p exit_code_fname still exists: JOB_ALIVE,
    !! JOB_GONE (its scheduler, or this host for a local job, no longer has it), or JOB_UNKNOWN (no
    !! record yet, an exit status written, or a local job of another host).
    integer function query_job( exit_code_fname ) result( state )
        class(string), intent(in) :: exit_code_fname
        character(len=256) :: host, line, sched, schedid
        type(string)       :: info_fname, cmd
        integer            :: funit, ios, pid, exitstat
        state = JOB_UNKNOWN
        if( file_exists(exit_code_fname) ) return
        info_fname = exit_code_fname//JOB_INFO_EXT
        if( .not. file_exists(info_fname) ) return
        open(newunit=funit, file=info_fname%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        read(funit,*,iostat=ios) pid
        if( ios == 0 ) read(funit,'(A)',iostat=ios) host
        if( ios == 0 ) read(funit,'(A)',iostat=ios) line
        close(funit)
        if( ios /= 0 ) return
        read(line,*,iostat=ios) sched, schedid
        if( ios /= 0 ) return
        ! each command exits 0 while the job exists, 1 once it is gone, 2 when it cannot tell (the
        ! scheduler unreachable, its command missing): a scheduler's answer that the id is unknown,
        ! or that the job has ended, is gone; any other failure is not
        select case(trim(sched))
            case('slurm')
                cmd = 'out=$(squeue -h -j '//trim(schedid)//' -o %T 2>&1) || '//&
                    &'{ echo "$out" | grep -qi "invalid job id" && exit 1; exit 2; }; [ -n "$out" ] || exit 1'
            case('lsf')
                cmd = 'out=$(bjobs '//trim(schedid)//' 2>&1); echo "$out" | grep -qi "not found" && exit 1; '//&
                    &'echo "$out" | grep -qE "PEND|RUN|SUSP" && exit 0; echo "$out" | grep -qE "DONE|EXIT" && exit 1; exit 2'
            case('pbs')
                cmd = 'out=$(qstat '//trim(schedid)//' 2>&1) && exit 0; '//&
                    &'echo "$out" | grep -qiE "unknown job|job has finished" && exit 1; exit 2'
            case DEFAULT
                cmd = 'if [ "$(hostname)" = "'//trim(host)//'" ]; then kill -0 '//int2str(pid)//' 2> /dev/null || exit 1; '//&
                    &'else exit 2; fi'
        end select
        call exec_cmdline('sh -c '''//cmd%to_char()//'''', suppress_errors=.true., exitstat=exitstat)
        select case(exitstat)
            case(0)
                state = JOB_ALIVE
            case(1)
                state = JOB_GONE
        end select
        ! the job may have ended between the checks: an exit status written meanwhile is its own
        if( file_exists(exit_code_fname) ) state = JOB_UNKNOWN
    end function query_job

    !> The liveness check of a job polled for its exit status @p exit_code_fname: at most every
    !! @p interval_s (JOB_LIVENESS_S) seconds since @p last_check, query_job; @p nmisses counts the
    !! checks in a row that found it gone. At the second, the job is lost: JOB_LOST_EXIT_CODE is
    !! written as its exit status (by temporary and rename) and .true. returned.
    logical function check_job_lost( exit_code_fname, last_check, nmisses, interval_s ) result( l_lost )
        class(string),     intent(in)    :: exit_code_fname
        integer,           intent(inout) :: last_check, nmisses
        integer, optional, intent(in)    :: interval_s
        integer :: tnow, interval, funit, ios
        l_lost   = .false.
        interval = JOB_LIVENESS_S
        if( present(interval_s) ) interval = interval_s
        tnow = simple_gettime()
        if( last_check > 0 .and. tnow - last_check < interval ) return
        last_check = tnow
        select case(query_job(exit_code_fname))
            case(JOB_ALIVE)
                nmisses = 0
            case(JOB_GONE)
                nmisses = nmisses + 1
            case DEFAULT
                return
        end select
        if( nmisses < 2 ) return
        open(newunit=funit, file=exit_code_fname%to_char()//'.tmp', status='replace', action='write', iostat=ios)
        if( ios /= 0 ) return
        write(funit,'(I0)') JOB_LOST_EXIT_CODE
        close(funit)
        call simple_rename(exit_code_fname//'.tmp', exit_code_fname)
        write(logfhandle,'(A,A,A)') '>>> WARNING: THE JOB OF ', exit_code_fname%to_char(),&
            &' IS GONE WITHOUT AN EXIT STATUS (KILLED BY ITS SCHEDULER, OR ITS NODE LOST); COUNTED AS FAILED'
        l_lost = .true.
    end function check_job_lost

    !> .true. when the job of @p exit_code_fname recorded itself and wrote no exit status: it may
    !! still be running (or queued again), or was killed without a status.
    logical function job_left_unfinished( exit_code_fname )
        class(string), intent(in) :: exit_code_fname
        job_left_unfinished = file_exists(exit_code_fname//JOB_INFO_EXT) .and. .not. file_exists(exit_code_fname)
    end function job_left_unfinished

    !> Before a job @p label runs in @p dir: a directory whose job of that label was left
    !! unfinished is set aside (set_aside_dir), and @p dir is made afresh, so the new job never
    !! shares a directory with one that may still run.
    subroutine fresh_job_dir( dir, label )
        class(string),    intent(in) :: dir
        character(len=*), intent(in) :: label
        if( dir_exists(dir) )then
            if( job_left_unfinished(dir//'/EXIT_CODE_'//label) ) call set_aside_dir(dir)
        endif
        call simple_mkdir(dir)
    end subroutine fresh_job_dir

    !> Moves @p dir to <dir>_unfinished<k>, the first k free; @p aside_out is the new name (empty
    !! when @p dir does not exist).
    subroutine set_aside_dir( dir, aside_out )
        class(string),           intent(in)  :: dir
        type(string), optional,  intent(out) :: aside_out
        type(string) :: aside
        integer      :: k
        if( present(aside_out) ) aside_out = ''
        if( .not. dir_exists(dir) ) return
        k = 1
        do
            aside = dir//'_unfinished'//int2str(k)
            if( .not. dir_exists(aside) ) exit
            k = k + 1
        enddo
        write(logfhandle,'(A,A,A,A)') '>>> MOVING ASIDE ', dir%to_char(), ', WHICH A JOB MAY STILL USE, TO ', aside%to_char()
        call simple_rename(dir, aside)
        if( present(aside_out) ) aside_out = aside
    end subroutine set_aside_dir

end module simple_qsys_job_record
