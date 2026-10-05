!@descr: the record a queued job's script keeps of itself, to cancel the job, and directories set aside when a job may still use them
!==============================================================================
! MODULE: simple_qsys_job_record
!
! PURPOSE:
!   An asynchronous job's script (simple_qsys_ctrl, generate_script_2 with an
!   exit-status file) records the job when it starts, in
!   <exit-status file>.job: its pid, its host, and its scheduler and job id
!   ("slurm <id>", "lsf <id>", "pbs <id>" or "none 0"). cancel_queued_job uses
!   that record: scancel, bkill or qdel for a scheduler job; a SIGTERM to the
!   script and its children when it runs on this host (a local job, or a
!   persistent-worker task on this host). A job still queued (no record yet),
!   one on another host without a scheduler id, and the jobs a program queues
!   itself are not cancelled.
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
use simple_defs,         only: logfhandle
use simple_defs_fname,   only: JOB_INFO_EXT
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_fileio,       only: file_exists
use simple_syslib,       only: simple_mkdir, dir_exists, simple_rename, exec_cmdline
implicit none

public :: cancel_queued_job, job_left_unfinished, fresh_job_dir, set_aside_dir
private

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
                ! a pid means nothing on another host
                cmd = 'if [ "$(hostname)" = "'//trim(host)//'" ]; then pkill -TERM -P '//int2str(pid)//&
                    &'; kill -TERM '//int2str(pid)//'; fi'
        end select
        call exec_cmdline(cmd%to_char()//' > /dev/null 2>&1', suppress_errors=.true.)
        l_sent = .true.
    end function cancel_queued_job

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
