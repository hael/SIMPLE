!@descr: POSIX fork-based child-process manager with timestamps and status polling
! Extend and override execute(cline); start() forks and runs it in the child (exit 0).
! status() polls waitpid(WNOHANG); a non-zero wait status (including signals) is FAILED.
! A failed child is not restarted here: the owner calls start() again.
! terminate() sends SIGTERM, kill() SIGKILL. The parent's buffered output is flushed before
! the fork, so a child does not write it again; the child ignores SIGINT, so a terminal's
! Ctrl-C reaches the parent, which stops its children in order. A child started with
! new_session leads its own session and process group, which signal_group() reaches whole,
! the processes the child started included.
module simple_forked_process
  use, intrinsic :: iso_c_binding,   only: c_intptr_t
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use unix,                  only: c_pid_t, c_int, c_long, c_null_char, &
                                  c_fork, c_kill, c_exit, c_time, c_setsid, &
                                  c_waitpid, c_usleep, c_perror,      &
                                  SIGTERM, SIGKILL, SIGINT, SIGHUP, EXIT_SUCCESS, &
                                  WNOHANG, c_signal, c_funptr, c_null_funptr
  use simple_defs,           only: logfhandle
  use simple_error,          only: simple_exception
  use simple_fileio,         only: fclose                  
  use simple_string,         only: string
  use simple_syslib,         only: file_exists
  use simple_cmdline,        only: cmdline
  use simple_memory_monitor, only: mem_monitor_init, mem_monitor_finish
  
  implicit none

  integer, public,  parameter :: FORK_STATUS_FAILED     = -1
  integer, public,  parameter :: FORK_STATUS_RUNNING    =  0
  integer, public,  parameter :: FORK_STATUS_STOPPED    =  1
  integer, public,  parameter :: FORK_STATUS_SKIPPED    =  3
  integer, public,  parameter :: FORK_POLL_TIME         = 100000 ! poll interval (µs)

  public  :: forked_process, sig_ign, is_sig_ign
  private
#include "simple_local_flags.inc"

  type :: forked_process
    private
    type(cmdline), allocatable :: cline ! what the child runs execute with (set_cline, or start's)
    type(string)          :: name
    type(string)          :: logfile
    integer(kind=c_pid_t) :: pid        = -1
    integer               :: queuetime  = 0  ! Unix timestamp: when process was queued
    integer               :: starttime  = 0  ! Unix timestamp: when process last started
    integer               :: stoptime   = 0  ! Unix timestamp: when process stopped cleanly
    integer               :: failtime   = 0  ! Unix timestamp: when process last failed
    logical               :: running    = .false.
    logical               :: failed     = .false.
    logical               :: stopped    = .false.
    logical               :: skipped    = .false.
    logical               :: new_session = .false. ! the child leads its own session
    integer(kind=c_int)   :: wstatus    = 0  ! wait status of the reaped child
  contains
    procedure :: execute => execute_test
    procedure :: start
    procedure :: terminate
    procedure :: kill => kill_forked_process
    procedure :: signal_group
    procedure :: get_exit_code
    procedure :: destroy
    procedure :: set_cline
    procedure :: skip
    procedure :: status
    procedure :: has_exited
    procedure :: await_final_status
    procedure :: get_pid
    procedure :: get_queuetime
    procedure :: get_starttime
    procedure :: get_stoptime
    procedure :: get_failtime
  end type forked_process

contains

  ! Fork a child process and begin execution. Optionally accept a new cline,
  ! name and logfile. In the child, lead a new session if new_session, redirect
  ! logfhandle if a logfile is given, call self%execute(), then exit. In the
  ! parent, record timestamps.
  subroutine start( self, name, logfile, cline, new_session )
    class(forked_process),           intent(inout) :: self
    type(string),          optional, intent(in)    :: name, logfile
    type(cmdline),         optional, intent(in)    :: cline
    logical,               optional, intent(in)    :: new_session
    integer(kind=c_int)                            :: ios
    type(c_funptr)                                 :: prev_handler
    if( present(logfile) ) self%logfile = logfile
    if( present(name)    ) self%name    = name
    if( present(new_session) ) self%new_session = new_session
    ! cmdline's defined assignment needs an allocated left-hand side
    if( .not. allocated(self%cline) ) allocate(self%cline)
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
    ! what the parent has buffered would otherwise be written by the child too
    flush(logfhandle,   iostat=ios)
    flush(output_unit,  iostat=ios)
    flush(error_unit,   iostat=ios)
    self%pid = c_fork()
    if( self%pid < 0 ) then
      ! Fork failed — terminal error.
      call c_perror('fork()' // c_null_char)
      THROW_HARD('Failed to fork process')
    else if( self%pid == 0 ) then
      ! Child process: optionally lead a session and redirect log output, execute, then exit.
      if( self%new_session )then
        if( c_setsid() < 0 ) call c_perror('setsid()' // c_null_char)
      end if
      ! Default SIGTERM first: a handler inherited from the parent acts on the parent's
      ! state and threads, which the child does not have; execute() installs its own.
      ! SIGINT is ignored, and stays ignored in what the child execs: a terminal's Ctrl-C
      ! goes to the whole process group, and the parent stops its children in order. SIGHUP
      ! takes its default too, so the parent's handler does not act in the child.
      prev_handler = c_signal(SIGTERM, c_null_funptr)
      prev_handler = c_signal(SIGHUP,  c_null_funptr)
      prev_handler = c_signal(SIGINT,  sig_ign())
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
      ! init memory monitor
      call mem_monitor_init(self%cline, 'simple_stream: ' // self%name%to_char())
      call self%execute(self%cline)
      if( .not. self%logfile%is_blank() ) call fclose(logfhandle)
      call mem_monitor_finish()
      call c_exit(0)
    else
      ! Parent process: record running state and timestamps.
      self%stoptime  = 0
      self%starttime = int(c_time(0_c_long))
      self%running   = .true.
      self%failed    = .false.
      self%stopped   = .false.
      if( self%queuetime == 0 ) self%queuetime = int(c_time(0_c_long))
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
  subroutine kill_forked_process( self )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: rc
    if( self%pid <= 0 ) return
    rc = c_kill(self%pid, SIGKILL)
    if( rc /= 0 ) THROW_HARD('Failed to send SIGKILL to forked child')
  end subroutine kill_forked_process

  ! Send sig to the child's process group (its session, with new_session), or to the child
  ! alone while it has not yet made the group; never to a reaped child, whose pid may be
  ! another process's by now. A group that is gone is not an error.
  subroutine signal_group( self, sig )
    class(forked_process), intent(in) :: self
    integer(kind=c_int),   intent(in) :: sig
    integer(kind=c_int)               :: rc
    if( self%pid <= 0 ) return
    rc = -1
    if( self%new_session ) rc = c_kill(-self%pid, sig)
    if( rc /= 0 .and. self%running ) rc = c_kill(self%pid, sig)
  end subroutine signal_group

  ! The exit status of the reaped child, or minus the signal that ended it.
  integer function get_exit_code( self )
    class(forked_process), intent(in) :: self
    integer :: termsig
    termsig = iand(int(self%wstatus), 127)
    if( termsig == 0 )then
      get_exit_code = iand(ishft(int(self%wstatus), -8), 255)
    else
      get_exit_code = -termsig
    endif
  end function get_exit_code

  ! Mark the process as skipped, which will cause status() to return
  ! FORK_STATUS_SKIPPED.
  subroutine skip( self )
    class(forked_process), intent(inout) :: self
    self%skipped  = .true.
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
        case( FORK_STATUS_RUNNING )
          rc = c_usleep(FORK_POLL_TIME)
        case( FORK_STATUS_STOPPED, FORK_STATUS_FAILED, FORK_STATUS_SKIPPED )
          exit
        case default
          THROW_HARD('Unknown fork status')
      end select
    end do
  end subroutine await_final_status

  ! Releases the command line; the process is not touched.
  subroutine destroy( self )
    class(forked_process), intent(inout) :: self
    if( allocated(self%cline) )then
      call self%cline%kill()
      deallocate(self%cline)
    end if
  end subroutine destroy

  ! The command line every later start runs execute with (one copy, kept here).
  subroutine set_cline( self, cline )
    class(forked_process), intent(inout) :: self
    class(cmdline),        intent(in)    :: cline
    if( .not. allocated(self%cline) ) allocate(self%cline)
    self%cline = cline
  end subroutine set_cline

  ! Non-blocking status poll. Uses waitpid(WNOHANG) to check whether the
  ! child has exited. Records stop/fail timestamps.
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
        self%wstatus = stat_loc
        if( stat_loc == 0 ) then
          self%stopped = .true.
          self%failed  = .false.
        else
          self%stopped  = .false.
          self%failed   = .true.
          ! when the failure is seen, not on every later poll
          self%failtime = int(c_time(0_c_long))
        end if
      end if
    end if
    if( self%stopped ) then
      status_code = FORK_STATUS_STOPPED
      if( self%stoptime == 0 ) self%stoptime = int(c_time(0_c_long))
    end if
    if( self%failed ) status_code = FORK_STATUS_FAILED
  end function status

  ! Whether the running child has ended, without reaping it (waitid with WNOWAIT, Linux): until
  ! status() reaps it, its pid and process group stay its own, so they can still be signalled
  ! without reaching a process that reused the pid. Elsewhere .false.: the owner reaps first.
  logical function has_exited( self )
    class(forked_process), intent(in) :: self
#if defined(__linux__)
    interface
      function c_waitid( idtype, id, info, options ) bind(c, name='waitid') result( rc )
        import :: c_int, c_pid_t
        integer(c_int),   value         :: idtype
        integer(c_pid_t), value         :: id
        integer(c_int),   intent(inout) :: info(32) ! siginfo_t, 128 bytes
        integer(c_int),   value         :: options
        integer(c_int)                  :: rc
      end function c_waitid
    end interface
    integer(c_int), parameter :: P_PID = 1, WEXITED = 4, WNOWAIT = int(z'01000000', c_int)
    integer(c_int) :: info(32)
#endif
    has_exited = .false.
    if( .not. self%running .or. self%pid <= 0 ) return
#if defined(__linux__)
    ! si_signo stays 0 when no child is waitable yet
    info = 0_c_int
    if( c_waitid(P_PID, self%pid, info, ior(WEXITED, ior(WNOWAIT, WNOHANG))) /= 0 ) return
    has_exited = info(1) /= 0_c_int
#endif
  end function has_exited

  ! Return the child's PID.
  function get_pid( self ) result( pid )
    class(forked_process), intent(in) :: self
    integer(kind=c_pid_t)             :: pid
    pid = self%pid
  end function get_pid

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

  ! SIG_IGN, which the unix bindings do not define: the handler value 1 on Linux and macOS.
  function sig_ign() result( handler )
    type(c_funptr) :: handler
    handler = transfer(1_c_intptr_t, handler)
  end function sig_ign

  ! Whether a disposition signal() returned is SIG_IGN.
  logical function is_sig_ign( handler )
    type(c_funptr), intent(in) :: handler
    is_sig_ign = transfer(handler, 0_c_intptr_t) == 1_c_intptr_t
  end function is_sig_ign

end module simple_forked_process
