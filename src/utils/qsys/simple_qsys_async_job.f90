!@descr: one program run submitted asynchronously through a qsys environment, with its completion and failure
!==============================================================================
! MODULE: simple_qsys_async_job
!
! PURPOSE:
!   Wraps qsys_env%exec_simple_prg_in_queue_async for a caller that starts a
!   program in a directory and polls it later. The script writes the
!   program's exit status to a file in that directory, so status() can tell a
!   finished run (DONE) from a crashed one (FAILED). Polling a completion
!   marker the program writes on success, as callers did before, cannot see a
!   crash: the marker never appears and the caller waits forever.
!
!   The job runs in its own directory (created if needed); the script and the
!   log are written there as ./distr_<label> and simple_log_<label>.
!==============================================================================
module simple_qsys_async_job
use simple_defs,     only: CWD_GLOB
use simple_string,   only: string
use simple_fileio,   only: del_file, file_exists, read_exit_code, simple_chdir, simple_getcwd
use simple_syslib,   only: simple_mkdir
use simple_cmdline,  only: cmdline
use simple_qsys_env, only: qsys_env
implicit none

public :: qsys_async_job
public :: ASYNC_JOB_IDLE, ASYNC_JOB_RUNNING, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
private

integer, parameter :: ASYNC_JOB_IDLE    = 0 ! not started
integer, parameter :: ASYNC_JOB_RUNNING = 1 ! submitted, no exit status yet
integer, parameter :: ASYNC_JOB_DONE    = 2 ! exited with status 0
integer, parameter :: ASYNC_JOB_FAILED  = 3 ! exited with a non-zero status

type :: qsys_async_job
    private
    type(string) :: dir             ! absolute job directory
    type(string) :: label
    type(string) :: exit_code_fname ! absolute
    logical      :: l_started = .false.
contains
    procedure :: start
    procedure :: status
    procedure :: get_dir
    procedure :: get_log
    procedure :: kill
end type qsys_async_job

contains

    !> Submits @p cline through @p qenv to run in @p dir; @p label names the script, the log and
    !! the exit-status file. A stale exit-status file from an earlier run is removed first.
    subroutine start( self, qenv, cline, dir, label, exec_bin )
        class(qsys_async_job),   intent(inout) :: self
        class(qsys_env),         intent(inout) :: qenv
        class(cmdline),          intent(in)    :: cline
        class(string),           intent(in)    :: dir
        character(len=*),        intent(in)    :: label
        class(string), optional, intent(in)    :: exec_bin
        type(string) :: cwd, cwd_job
        call self%kill()
        call simple_getcwd(cwd)
        call simple_mkdir(dir)
        call simple_chdir(dir)
        call simple_getcwd(cwd_job)
        CWD_GLOB             = cwd_job%to_char()
        self%dir             = cwd_job
        self%label           = label
        self%exit_code_fname = cwd_job//'/EXIT_CODE_'//label
        if( file_exists(self%exit_code_fname) ) call del_file(self%exit_code_fname)
        call qenv%exec_simple_prg_in_queue_async(cline, string('./distr_'//label), string('simple_log_'//label),&
            &exec_bin=exec_bin, exit_code_fname=self%exit_code_fname)
        call simple_chdir(cwd)
        CWD_GLOB       = cwd%to_char()
        self%l_started = .true.
    end subroutine start

    !> ASYNC_JOB_IDLE, _RUNNING, _DONE or _FAILED. An exit-status file that cannot be read yet
    !! (being written) counts as running.
    integer function status( self )
        class(qsys_async_job), intent(in) :: self
        integer :: exit_code
        logical :: err
        status = ASYNC_JOB_IDLE
        if( .not. self%l_started ) return
        status = ASYNC_JOB_RUNNING
        if( .not. file_exists(self%exit_code_fname) ) return
        call read_exit_code(self%exit_code_fname, exit_code, err)
        if( err ) return
        if( exit_code == 0 )then
            status = ASYNC_JOB_DONE
        else
            status = ASYNC_JOB_FAILED
        endif
    end function status

    !> The absolute directory the job runs in.
    function get_dir( self ) result( dir )
        class(qsys_async_job), intent(in) :: self
        type(string) :: dir
        dir = self%dir
    end function get_dir

    !> The absolute path of the job's log.
    function get_log( self ) result( fname )
        class(qsys_async_job), intent(in) :: self
        type(string) :: fname
        fname = self%dir//'/simple_log_'//self%label
    end function get_log

    subroutine kill( self )
        class(qsys_async_job), intent(inout) :: self
        call self%dir%kill
        call self%label%kill
        call self%exit_code_fname%kill
        self%l_started = .false.
    end subroutine kill

end module simple_qsys_async_job
