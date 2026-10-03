!@descr: one stream stage as the master runs it: its forked process, command line and master-side pipes
!==============================================================================
! MODULE: simple_stream_master_stage
!
! PURPOSE:
!   The master keeps one stream_master_stage per stage (simple_stream_master_stage_ids
!   gives the ids). It starts, stops and restarts the stage's forked process,
!   reads the stage's GUI messages and sends it the GUI updates, through two
!   stream_pipes on the stage's pipes. One fork type serves every stage: its
!   execute closes the other stages' pipe ends and runs the commander the
!   master gave the stage, so this module knows commanders only as
!   commander_base.
!
! LIFECYCLE:
!   new(id, commander, cline, l_updates, max_frame_bytes) -> start / request_stop /
!   force_stop / discard_pipes + start ... -> kill
!==============================================================================
module simple_stream_master_stage
use simple_defs,                    only: logfhandle
use simple_error,                   only: simple_exception
use simple_string,                  only: string
use simple_cmdline,                 only: cmdline
use simple_commander_base,          only: commander_base
use simple_forked_process,          only: forked_process, FORK_STATUS_RUNNING
use simple_stream_pipe,             only: stream_pipe
use simple_stream_master_stage_ids, only: stage_job_name, stage_label, open_stage_pipes, close_stage_pipes,&
                                         &close_other_pipe_ends, master_fds, stage_fds
implicit none

public :: stream_master_stage, stream_master_stage_fork
private
#include "simple_local_flags.inc"

!> The forked process of one stage, running the stage's commander.
type, extends(forked_process) :: stream_master_stage_fork
    integer                            :: id = 0
    class(commander_base), allocatable :: commander
contains
    procedure :: execute => execute_stage
end type stream_master_stage_fork

type :: stream_master_stage
    integer                        :: id        = 0
    type(stream_master_stage_fork) :: fork
    type(cmdline)                  :: cline
    type(stream_pipe)              :: reader              ! the stage's GUI messages
    type(stream_pipe)              :: writer              ! the GUI updates, to a stage that reads them
    logical                        :: l_updates = .false.
    logical                        :: l_exists  = .false.
contains
    procedure :: new
    procedure :: start
    procedure :: skip
    procedure :: is_running
    procedure :: request_stop
    procedure :: force_stop
    procedure :: discard_pipes
    procedure :: send_update
    procedure :: get_label
    procedure :: kill
end type stream_master_stage

contains

    !---------------- the fork ----------------

    ! In the forked process: only this stage's pipe ends stay open, then its commander runs.
    subroutine execute_stage( self, cline )
        class(stream_master_stage_fork), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        call close_other_pipe_ends(self%id)
        if( .not. allocated(self%commander) ) THROW_HARD('no commander for stream stage '//stage_label(self%id))
        call self%commander%execute(cline)
    end subroutine execute_stage

    !---------------- lifecycle ----------------

    !> Stage @p id, run by @p commander with command line @p cline: its pipes are made, and the
    !! master's ends wrapped; @p l_updates when the stage reads the GUI updates.
    subroutine new( self, id, commander, cline, l_updates, max_frame_bytes )
        class(stream_master_stage), intent(inout) :: self
        integer,                    intent(in)    :: id, max_frame_bytes
        class(commander_base),      intent(in)    :: commander
        class(cmdline),             intent(in)    :: cline
        logical,                    intent(in)    :: l_updates
        integer :: fd_read, fd_write
        call self%kill()
        self%id      = id
        self%fork%id = id
        allocate(self%fork%commander, source=commander)
        self%cline   = cline
        self%l_updates = l_updates
        call open_stage_pipes(id)
        call master_fds(id, fd_read, fd_write)
        call self%reader%new(fd_read, -1, max_frame_bytes, stage_job_name(id))
        call self%writer%new(-1, fd_write, max_frame_bytes, stage_job_name(id))
        self%l_exists = .true.
    end subroutine new

    !> Forks the stage, first or again. A stage that fails stays failed: restarts are the GUI's
    !! call, so the fork never restarts by itself (forked_process would, from inside any status
    !! poll, and would undo a force_stop).
    subroutine start( self )
        class(stream_master_stage), intent(inout) :: self
        call self%fork%start(name=string(stage_job_name(self%id)), logfile=string(stage_job_name(self%id)//'.log'),&
            &cline=self%cline, restart=.false.)
    end subroutine start

    !> The stage is not run (its output exists already).
    subroutine skip( self )
        class(stream_master_stage), intent(inout) :: self
        call self%fork%skip()
    end subroutine skip

    logical function is_running( self )
        class(stream_master_stage), intent(inout) :: self
        is_running = self%fork%status() == FORK_STATUS_RUNNING
    end function is_running

    !> Asks a running stage to stop (SIGTERM); it finishes its pass and finalises.
    subroutine request_stop( self )
        class(stream_master_stage), intent(inout) :: self
        if( self%is_running() ) call self%fork%terminate()
    end subroutine request_stop

    !> Ends a running stage at once (SIGKILL), for one that did not stop when asked.
    subroutine force_stop( self )
        class(stream_master_stage), intent(inout) :: self
        if( .not. self%is_running() ) return
        write(logfhandle,'(A)') '>>> '//stage_label(self%id)//' DID NOT STOP WHEN ASKED; KILLING IT'
        call self%fork%kill()
    end subroutine force_stop

    !> Before a restart: drops what the stopped stage left in its pipes, both its unread messages
    !! and the updates it never read, so the new process and the master start on a frame boundary.
    !! The caller holds the lock of the thread that reads the stage's messages.
    subroutine discard_pipes( self, max_frame_bytes )
        class(stream_master_stage), intent(inout) :: self
        integer,                  intent(in)    :: max_frame_bytes
        type(stream_pipe) :: stale_updates
        integer           :: fd_write, fd_read
        call self%reader%discard()
        call stage_fds(self%id, fd_write, fd_read)
        call stale_updates%new(fd_read, -1, max_frame_bytes, 'stale updates')
        call stale_updates%discard()
        call stale_updates%kill()
    end subroutine discard_pipes

    !> Sends a serialised GUI update to the stage, when it reads them and runs: the master keeps
    !! every pipe end open, so a stopped stage's pipe would only fill up.
    subroutine send_update( self, buffer )
        class(stream_master_stage), intent(inout) :: self
        character(len=*),         intent(in)    :: buffer
        if( .not. self%l_updates ) return
        if( .not. self%is_running() ) return
        call self%writer%send(buffer)
    end subroutine send_update

    function get_label( self ) result( label )
        class(stream_master_stage), intent(in) :: self
        character(len=:), allocatable :: label
        label = stage_label(self%id)
    end function get_label

    !> Closes the master's pipe ends of the stage; the process is not touched.
    subroutine kill( self )
        class(stream_master_stage), intent(inout) :: self
        if( .not. self%l_exists ) return
        call self%reader%kill()
        call self%writer%kill()
        call close_stage_pipes(self%id)
        call self%cline%kill()
        if( allocated(self%fork%commander) ) deallocate(self%fork%commander)
        self%l_updates = .false.
        self%l_exists  = .false.
    end subroutine kill

end module simple_stream_master_stage
