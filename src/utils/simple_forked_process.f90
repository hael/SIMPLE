!@descr: POSIX fork-based child-process manager with timestamps, auto-restart, and status polling
!==============================================================================
! MODULE: simple_forked_process
!
! PURPOSE:
!   Provides the forked_process type, which wraps a POSIX fork()/waitpid()
!   lifecycle. The parent retains a handle to the child, can poll its status,
!   send signals, and optionally restart it on failure up to FORK_MAX_RESTARTS
!   times. Unix timestamps are recorded at queue, start, stop, and fail events.
!
! TYPES:
!   forked_process — owns a single child PID; bind a commander-style execute()
!                    procedure via type extension or direct override before
!                    calling start().
!
! STATUS CODES (public parameters):
!   FORK_STATUS_FAILED     (-1) — child exited non-zero; restart exhausted
!                                  or disabled
!   FORK_STATUS_RUNNING    ( 0) — child is still running
!   FORK_STATUS_STOPPED    ( 1) — child exited cleanly (exit code 0)
!   FORK_STATUS_RESTARTING ( 2) — child failed and is being restarted
!   FORK_STATUS_SKIPPED    ( 3) — process was skipped (never started)
!
! PARAMETERS (hard-coded):
!   FORK_POLL_TIME    — usleep interval for status polling (µs)    (100 000)
!   FORK_MAX_RESTARTS — maximum automatic restarts before giving up     (10)
!
! DEPENDENCIES:
!   unix, simple_string, simple_syslib, simple_string_utils, simple_cmdline
!==============================================================================
module simple_forked_process
  use unix,                only: c_pid_t, c_int, c_long, c_null_char, &
                                  c_fork, c_kill, c_exit, c_time,     &
                                  c_waitpid, c_usleep, c_perror,      &
                                  SIGTERM, SIGKILL, EXIT_SUCCESS,     &
                                  WNOHANG
  use simple_defs,         only: logfhandle
  use simple_error,        only: simple_exception
  use simple_fileio,       only: fclose                  
  use simple_string,       only: string
  use simple_syslib,       only: file_exists
  use simple_cmdline,      only: cmdline
  use simple_string_utils, only: int2str
  
  implicit none

  integer, public,  parameter :: FORK_STATUS_FAILED     = -1
  integer, public,  parameter :: FORK_STATUS_RUNNING    =  0
  integer, public,  parameter :: FORK_STATUS_STOPPED    =  1
  integer, public,  parameter :: FORK_STATUS_RESTARTING =  2
  integer, public,  parameter :: FORK_STATUS_SKIPPED    =  3
  integer, public,  parameter :: FORK_POLL_TIME         = 100000 ! poll interval (µs)
  integer, private, parameter :: FORK_MAX_RESTARTS      = 10     ! max auto-restarts

  public  :: forked_process
  private
#include "simple_local_flags.inc"

  type :: forked_process
    private
    type(cmdline)         :: cline
    type(string)          :: name
    type(string)          :: description
    type(string)          :: logfile
    integer(kind=c_pid_t) :: pid        = -1
    integer               :: n_restarts = 0
    integer               :: queuetime  = 0  ! Unix timestamp: when process was queued
    integer               :: starttime  = 0  ! Unix timestamp: when process last started
    integer               :: stoptime   = 0  ! Unix timestamp: when process stopped cleanly
    integer               :: failtime   = 0  ! Unix timestamp: when process last failed
    logical               :: running    = .false.
    logical               :: failed     = .false.
    logical               :: stopped    = .false.
    logical               :: skipped    = .false.
    logical               :: restart    = .false. ! enable auto-restart on failure
  contains
    procedure :: execute => execute_test
    procedure :: start
    procedure :: terminate
    procedure :: kill
    procedure :: destroy
    procedure :: skip
    procedure :: status
    procedure :: await_final_status
    procedure :: get_pid
    procedure :: get_nrestarts
    procedure :: get_queuetime
    procedure :: get_starttime
    procedure :: get_stoptime
    procedure :: get_failtime
  end type forked_process

contains

  ! Fork a child process and begin execution. Optionally accept a new cline,
  ! name, logfile, and restart flag. In the child, redirect logfhandle if a
  ! logfile is given, call self%execute(), then exit. In the parent, record
  ! timestamps and build a human-readable description.
  subroutine start( self, restart, name, logfile, cline )
    class(forked_process),           intent(inout) :: self
    logical,               optional, intent(in)    :: restart
    type(string),          optional, intent(in)    :: name, logfile
    type(cmdline),         optional, intent(in)    :: cline
    integer(kind=c_int)                            :: ios
    if( present(restart) ) self%restart = restart
    if( present(logfile) ) self%logfile = logfile
    if( present(name)    ) self%name    = name
    if( present(cline)   ) self%cline   = cline
#if defined(_WIN32)
      self%pid      = -1
      self%running  = .false.
      self%failed   = .false.
      self%stopped  = .false.
      self%skipped  = .true.
      return
#endif
    self%skipped = .false.
    self%pid = c_fork()
    if( self%pid < 0 ) then
      ! Fork failed — terminal error.
      call c_perror('fork()' // c_null_char)
      THROW_HARD('Failed to fork process')
    else if( self%pid == 0 ) then
      ! Child process: optionally redirect log output, execute, then exit.
      if( .not. self%logfile%is_blank() ) then
        if( file_exists(self%logfile%to_char()) ) then
          open(UNIT=logfhandle, FILE=self%logfile%to_char(), IOSTAT=ios, &
               ACTION='WRITE', STATUS='OLD',  POSITION='APPEND')
        else
          open(UNIT=logfhandle, FILE=self%logfile%to_char(), IOSTAT=ios, &
               ACTION='WRITE', STATUS='NEW',  POSITION='APPEND')
        end if
        if( ios /= 0 ) THROW_HARD('Failed to open logfile')
      end if
      call self%execute(self%cline)
      if( .not. self%logfile%is_blank() ) call fclose(logfhandle)
      call c_exit(0)
    else
      ! Parent process: record running state, timestamps, and description.
      self%stoptime  = 0
      self%starttime = int(c_time(0_c_long))
      self%running   = .true.
      self%failed    = .false.
      self%stopped   = .false.
      if( self%queuetime == 0 ) self%queuetime = int(c_time(0_c_long))
      if( self%name%is_blank() ) then
        self%description = 'PID ' // int2str(self%pid)
      else
        self%description = self%name%to_char() // ' (PID:' // int2str(self%pid) // ')'
      end if
    end if
  end subroutine start

  ! Send SIGTERM to the child, requesting a graceful shutdown.
  subroutine terminate( self )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: rc
    if( self%pid < 0 ) return
    rc = c_kill(self%pid, SIGTERM)
    if( rc /= 0 ) THROW_HARD('Failed to send SIGTERM to forked child')
  end subroutine terminate

  ! Send SIGKILL to the child, forcing immediate termination.
  subroutine kill( self )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: rc
    if( self%pid < 0 ) return
    rc = c_kill(self%pid, SIGKILL)
    if( rc /= 0 ) THROW_HARD('Failed to send SIGKILL to forked child')
  end subroutine kill

  ! Mark the process as skipped, which will cause status() to return
  ! FORK_STATUS_SKIPPED and prevent future restarts.
  subroutine skip( self )
    class(forked_process), intent(inout) :: self
    self%skipped  = .true.
    self%restart  = .false.
  end subroutine skip

  ! Default execute implementation used for testing. Installs a SIGTERM
  ! handler that flushes the log and exits cleanly, writes a sentinel line
  ! to logfhandle, then sleeps for 20 poll intervals.
  ! NOTE: changing the sentinel line will break the hash check in
  !       test_logfile_redirection.
  subroutine execute_test( self, cline )
    class(forked_process), intent(inout) :: self
    class(cmdline),        intent(inout) :: cline
    integer                              :: rc
    call signal(SIGTERM, sigterm_handler)
    if( .not. self%logfile%is_blank() ) write(logfhandle, '(A)') 'LOGFILE CONTENTS TEST'
    rc = c_usleep(FORK_POLL_TIME * 20)
  contains
    subroutine sigterm_handler()
      call flush(logfhandle)
      call exit(EXIT_SUCCESS)
    end subroutine sigterm_handler
  end subroutine execute_test

  ! Block until the child reaches a terminal state (STOPPED or FAILED),
  ! sleeping FORK_POLL_TIME µs between status checks.
  subroutine await_final_status( self )
    class(forked_process), intent(inout) :: self
    integer                              :: rc
    do
      select case( self%status() )
        case( FORK_STATUS_RUNNING, FORK_STATUS_RESTARTING )
          rc = c_usleep(FORK_POLL_TIME)
        case( FORK_STATUS_STOPPED, FORK_STATUS_FAILED, FORK_STATUS_SKIPPED )
          exit
        case default
          THROW_HARD('Unknown fork status')
      end select
    end do
  end subroutine await_final_status

  ! No-op destructor placeholder.
  subroutine destroy( self )
    class(forked_process), intent(inout) :: self
  end subroutine destroy

  ! Non-blocking status poll. Uses waitpid(WNOHANG) to check whether the
  ! child has exited. Records stop/fail timestamps. On failure, auto-restarts
  ! the child up to FORK_MAX_RESTARTS times if self%restart is set.
  function status( self ) result( status_code )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: options, stat_loc, rc
    integer                              :: status_code
    options     = WNOHANG
    status_code = FORK_STATUS_RUNNING
    if( self%skipped )then
      status_code = FORK_STATUS_SKIPPED
      return
    endif
    if( self%running ) then
      rc = c_waitpid(self%pid, stat_loc, options)
      if( rc == self%pid ) then
        self%running = .false.
        if( stat_loc == 0 ) then
          self%stopped = .true.
          self%failed  = .false.
        else
          self%stopped = .false.
          self%failed  = .true.
        end if
      end if
    end if
    if( self%stopped ) then
      self%n_restarts = 0
      status_code     = FORK_STATUS_STOPPED
      if( self%stoptime == 0 ) self%stoptime = int(c_time(0_c_long))
    end if
    if( self%failed ) then
      self%failtime = int(c_time(0_c_long))
      if( self%restart ) then
        if( self%n_restarts > FORK_MAX_RESTARTS ) then
          THROW_WARN('max restarts reached for forked process ' // self%description%to_char())
          status_code = FORK_STATUS_FAILED
        else
          call self%start()
          status_code = FORK_STATUS_RESTARTING
          self%n_restarts = self%n_restarts + 1
        end if
      else
        status_code = FORK_STATUS_FAILED
      end if
    end if
  end function status

  ! Return the child's PID.
  function get_pid( self ) result( pid )
    class(forked_process), intent(in) :: self
    integer(kind=c_pid_t)             :: pid
    pid = self%pid
  end function get_pid

  ! Return the number of times the child has been restarted.
  function get_nrestarts( self ) result( n_restarts )
    class(forked_process), intent(in) :: self
    integer                           :: n_restarts
    n_restarts = self%n_restarts
  end function get_nrestarts

  ! Return the Unix timestamp at which the process was first queued.
  function get_queuetime( self ) result( queuetime )
    class(forked_process), intent(in) :: self
    integer                           :: queuetime
    queuetime = self%queuetime
  end function get_queuetime

  ! Return the Unix timestamp of the most recent start() call.
  function get_starttime( self ) result( starttime )
    class(forked_process), intent(in) :: self
    integer                           :: starttime
    starttime = self%starttime
  end function get_starttime

  ! Return the Unix timestamp at which the child stopped cleanly (0 if not yet).
  function get_stoptime( self ) result( stoptime )
    class(forked_process), intent(in) :: self
    integer                           :: stoptime
    stoptime = self%stoptime
  end function get_stoptime

  ! Return the Unix timestamp of the most recent failure (0 if never failed).
  function get_failtime( self ) result( failtime )
    class(forked_process), intent(in) :: self
    integer                           :: failtime
    failtime = self%failtime
  end function get_failtime

end module simple_forked_process
