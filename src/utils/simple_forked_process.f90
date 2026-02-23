!@descr: various unix forking utilities
module simple_forked_process
use unix
use simple_core_module_api
use simple_commander_base
use simple_cmdline
use simple_fileio
use simple_test_utils
use simple_ipc_mq
implicit none

integer, public,  parameter :: FORK_STATUS_FAILED     = -1
integer, public,  parameter :: FORK_STATUS_RUNNING    = 0
integer, public,  parameter :: FORK_STATUS_STOPPED    = 1
integer, public,  parameter :: FORK_STATUS_RESTARTING = 2
integer, public,  parameter :: FORK_POLL_TIME         = 100000 ! poll time in microseconds
integer, private, parameter :: FORK_MAX_RESTARTS      = 10 

public :: forked_process

private
#include "simple_local_flags.inc"

type :: forked_process
  private
  type(cmdline)         :: cline
  type(string)          :: name
  type(string)          :: description
  type(string)          :: logfile
  integer(kind=c_pid_t) :: pid          = -1
  integer               :: n_restarts   = 0
  logical               :: running      = .false.
  logical               :: failed       = .false.
  logical               :: stopped      = .false.
  logical               :: restart      = .false. ! triggers auto restart upon fail

contains
  procedure :: execute => execute_test
  procedure :: start
  procedure :: terminate
  procedure :: kill
  procedure :: destroy
  procedure :: status
  procedure :: await_final_status
  procedure :: get_pid
end type forked_process

contains

  subroutine start( self, restart, name, logfile)
    class(forked_process),           intent(inout) :: self
    type(string),          optional, intent(in)    :: name, logfile
    logical,               optional, intent(in)    :: restart
    integer(kind=c_int)                            :: ios
    if( present(restart) ) self%restart = restart
    if( present(logfile) ) self%logfile = logfile
    if( present(name)    ) self%name    = name
    self%pid = c_fork()
    if( self%pid < 0 ) then
      ! Fork failed. Terminal
      call c_perror('fork()' // c_null_char)
      THROW_HARD( 'Failed to fork process' )
    else if ( self%pid == 0 ) then
      ! Child process.
      if( .not. self%logfile%is_blank() ) then
        ! redirect logfhandle to logfile
        if( file_exists(self%logfile%to_char()) ) then
          open(UNIT=logfhandle, FILE=self%logfile%to_char(), IOSTAT=ios, ACTION='WRITE', STATUS='OLD', POSITION='APPEND')
        else
          open(UNIT=logfhandle, FILE=self%logfile%to_char(), IOSTAT=ios, ACTION='WRITE', STATUS='NEW', POSITION='APPEND')
        end if
        if( ios .ne. 0 ) THROW_HARD('Failed to open logfile')
      end if
      ! execute 
      call self%execute(self%cline)
      ! close logfile
      if( .not. self%logfile%is_blank() ) call fclose(logfhandle)
      ! exit cleanly
      call c_exit(0)
    else
      ! Parent process.
      self%running = .true.
      self%failed  = .false.
      self%stopped = .false.
      if( self%name%is_blank() ) then
        self%description = 'PID ' // int2str(self%pid)
      else
        self%description = self%name%to_char() // '( PID:' // int2str(self%pid) // ' )'
      end if
     end if
  end subroutine start

  subroutine terminate( self )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: sig = SIGTERM
    integer(kind=c_int)                  :: rc
    rc = c_kill(self%pid, sig)
    if( rc /= 0 ) then
      THROW_HARD('Failed to send kill signal to forked child')
    end if 
  end subroutine terminate

  subroutine kill( self )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: rc
    rc = c_kill(self%pid, SIGKILL)
    if( rc /= 0 ) then
      THROW_HARD('Failed to send kill signal to forked child')
    end if 
  end subroutine kill

  subroutine execute_test( self, cline )
    class(forked_process), intent(inout) :: self
    class(cmdline),        intent(inout) :: cline
    type(string)                         :: mq_msg
    integer                              :: rc
    call signal(SIGTERM, sigterm_handler)
    ! note: changing this output will require updating the log hash comparison in test_logfile_redirection
    if( .not. self%logfile%is_blank() ) write(logfhandle, '(A)') "LOGFILE CONTENTS TEST"
    rc = c_usleep(FORK_POLL_TIME * 20)

    contains

    subroutine sigterm_handler()
      call flush(logfhandle)
      call exit(EXIT_SUCCESS)
    end subroutine sigterm_handler

  end subroutine execute_test

  subroutine await_final_status( self )
    class(forked_process), intent(inout) :: self
    integer                              :: rc
    do while(.true.)
      select case(self%status())
        case(FORK_STATUS_RUNNING)
          rc = c_usleep( FORK_POLL_TIME )
          cycle
        case(FORK_STATUS_FAILED)
          exit
        case(FORK_STATUS_STOPPED)
          exit
        case(FORK_STATUS_RESTARTING)
          cycle
        case DEFAULT
          THROW_HARD('Unknown fork status')
      end select
    end do
  end subroutine await_final_status

  subroutine destroy( self )
    class(forked_process), intent(inout) :: self
  end subroutine destroy

  function status( self ) result( status_code )
    class(forked_process), intent(inout) :: self
    integer(kind=c_int)                  :: options = 1 !WNOHANG
    integer(kind=c_int)                  :: stat_loc, rc
    integer                              :: status_code
    if( self%running ) then
      status_code = FORK_STATUS_RUNNING
      rc = c_waitpid(self%pid, stat_loc, options)
      if( rc == self%pid ) then
        self%running = .false.
        if( stat_loc == 0 ) then
          self%stopped = .true.
          self%failed  = .false.
        else
          ! non-zero exit code
          self%stopped = .false.
          self%failed  = .true.
        end if
      end if
    end if
    if( self%stopped ) then
      self%n_restarts = 0
      status_code = FORK_STATUS_STOPPED
    end if
    if( self%failed ) then
      if( self%restart ) then
        self%n_restarts = self%n_restarts + 1
        if( self%n_restarts > FORK_MAX_RESTARTS ) then
          THROW_WARN('max restarts reached for forked thread ' // self%description%to_char())
          status_code = FORK_STATUS_FAILED
        else
          call self%start()
          status_code = FORK_STATUS_RESTARTING
        end if
      else
        status_code = FORK_STATUS_FAILED
      end if
    end if
  end function status

  function get_pid( self ) result( pid )
    class(forked_process), intent(inout) :: self
    integer(kind=c_pid_t)                :: pid
    pid = self%pid
  end function get_pid

end module simple_forked_process