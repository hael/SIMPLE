!@descr: one stream stage as the master runs it: its forked process, command line and master-side pipes
!==============================================================================
! MODULE: simple_stream_master_stage
!
! PURPOSE:
!   The master keeps one stream_master_stage per stage (simple_stream_master_stage_ids
!   gives the ids). It starts, stops and restarts the stage's forked process,
!   reads the stage's GUI messages and sends it the GUI updates, through two
!   stream_pipes on the stage's pipes. NICE answers every heartbeat with its
!   whole state, so an update goes to a stage only when it differs from the
!   last one the stage was sent since it started, and never to a stage asked
!   to stop. The update writer abandons a frame still part-written after
!   UPDATE_PARTIAL_RETRIES (a stage that stopped reading must not hold the
!   master); the channel is mended by the discard before the next start.
!   One fork type serves every stage: its
!   execute closes the other stages' pipe ends and runs the commander the
!   master gave the stage, so this module knows commanders only as
!   commander_base. fork_gui_status reads a forked process into the record
!   the GUI heartbeat takes.
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
use simple_forked_process,          only: forked_process, FORK_STATUS_RUNNING, FORK_STATUS_SKIPPED,&
                                         &FORK_STATUS_FAILED, FORK_STATUS_STOPPED
use simple_gui_assembler,           only: gui_stage_status, GUI_STAGE_STATUS_UNKNOWN, GUI_STAGE_STATUS_RUNNING,&
                                         &GUI_STAGE_STATUS_FAILED, GUI_STAGE_STATUS_FINISHED, GUI_STAGE_STATUS_SKIPPED
use simple_qsys_job_record,         only: cancel_unfinished_jobs
use simple_stream_pipe,             only: stream_pipe
use simple_stream_master_stage_ids, only: stage_job_name, stage_label, open_stage_pipes, close_stage_pipes,&
                                         &close_other_pipe_ends, master_fds, stage_fds
implicit none

public :: stream_master_stage, stream_master_stage_fork, fork_gui_status
private
#include "simple_local_flags.inc"

integer, parameter :: UPDATE_PARTIAL_RETRIES = 200 ! about 2 s of EAGAIN retries before a part-written update is abandoned

!> The forked process of one stage, running the stage's commander.
type, extends(forked_process) :: stream_master_stage_fork
    integer                            :: id = 0
    class(commander_base), allocatable :: commander
contains
    procedure :: execute => execute_stage
end type stream_master_stage_fork

type :: stream_master_stage
    integer                        :: id        = 0
    type(stream_master_stage_fork) :: fork                ! holds the stage's command line (set_cline)
    type(stream_pipe)              :: reader              ! the stage's GUI messages
    type(stream_pipe)              :: writer              ! the GUI updates, to a stage that reads them
    character(len=:), allocatable  :: last_update         ! the last update sent since the stage started
    logical                        :: l_updates = .false.
    logical                        :: l_stop_requested = .false.
    logical                        :: l_exists  = .false.
contains
    procedure :: new
    procedure :: start
    procedure :: skip
    procedure :: is_skipped
    procedure :: is_running
    procedure :: request_stop
    procedure :: force_stop
    procedure :: discard_pipes
    procedure :: send_update
    procedure :: is_new_update
    procedure :: get_label
    procedure :: kill
end type stream_master_stage

contains

    !---------------- the fork ----------------

    ! In the forked process: only this stage's pipe ends stay open, the master's persistent-worker
    ! server and its warm-up registrants are forgotten (the stage reaches the server as a client,
    ! through worker_server), then its commander runs.
    subroutine execute_stage( self, cline )
        use simple_persistent_worker_server, only: forget_inherited_persistent_worker
        use simple_qsys_env,                 only: forget_warmup_envs
        class(stream_master_stage_fork), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        call close_other_pipe_ends(self%id)
        call forget_inherited_persistent_worker()
        call forget_warmup_envs()
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
        call self%fork%set_cline(cline)
        self%l_updates = l_updates
        call open_stage_pipes(id)
        call master_fds(id, fd_read, fd_write)
        call self%reader%new(fd_read, -1, max_frame_bytes, stage_job_name(id))
        call self%writer%new(-1, fd_write, max_frame_bytes, stage_job_name(id))
        call self%writer%limit_partial_frames(UPDATE_PARTIAL_RETRIES)
        self%l_exists = .true.
    end subroutine new

    !> Forks the stage, first or again. A stage that fails stays failed: restarts are the GUI's
    !! call (forked_process never restarts a child by itself). A skipped stage is never started:
    !! its output is the user's earlier run.
    subroutine start( self )
        class(stream_master_stage), intent(inout) :: self
        if( self%is_skipped() )then
            write(logfhandle,'(A)') '>>> '//stage_label(self%id)//' IS SKIPPED AND IS NOT STARTED'
            return
        endif
        if( allocated(self%last_update) ) deallocate(self%last_update)
        self%l_stop_requested = .false.
        call self%fork%start(name=string(stage_job_name(self%id)), logfile=string(stage_job_name(self%id)//'.log'))
    end subroutine start

    !> The stage is not run (its output exists already).
    subroutine skip( self )
        class(stream_master_stage), intent(inout) :: self
        call self%fork%skip()
    end subroutine skip

    !> .true. for a stage that is not run (skip); it stays so.
    logical function is_skipped( self )
        class(stream_master_stage), intent(inout) :: self
        is_skipped = self%fork%status() == FORK_STATUS_SKIPPED
    end function is_skipped

    logical function is_running( self )
        class(stream_master_stage), intent(inout) :: self
        is_running = self%fork%status() == FORK_STATUS_RUNNING
    end function is_running

    !> Asks a running stage to stop (SIGTERM); it finishes its pass and finalises. It is sent no
    !! more updates until it is started again.
    subroutine request_stop( self )
        class(stream_master_stage), intent(inout) :: self
        self%l_stop_requested = .true.
        if( self%is_running() ) call self%fork%terminate()
    end subroutine request_stop

    !> Ends a running stage at once (SIGKILL), for one that did not stop when asked. A killed stage
    !! cannot cancel its jobs, so the master cancels those recorded in the stage's folder that
    !! wrote no exit status.
    subroutine force_stop( self )
        class(stream_master_stage), intent(inout) :: self
        integer :: ncancelled
        if( .not. self%is_running() ) return
        write(logfhandle,'(A)') '>>> '//stage_label(self%id)//' DID NOT STOP WHEN ASKED; KILLING IT'
        call self%fork%kill()
        call cancel_unfinished_jobs(string(stage_job_name(self%id)), ncancelled)
        if( ncancelled > 0 ) write(logfhandle,'(A,I6,A)') '>>> CANCELLED ', ncancelled,&
            &' JOBS OF '//stage_label(self%id)
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
        call self%writer%discard() ! mends a channel broken by an abandoned update
        call stage_fds(self%id, fd_write, fd_read)
        call stale_updates%new(fd_read, -1, max_frame_bytes, 'stale updates')
        call stale_updates%discard()
        call stale_updates%kill()
    end subroutine discard_pipes

    !> Sends a serialised GUI update to the stage, when it reads them, runs, has not been asked to
    !! stop, and the update differs from the last one it was sent: the master keeps every pipe end
    !! open, so a stopped stage's pipe would only fill up, and NICE repeats its whole state in
    !! every answer.
    subroutine send_update( self, buffer )
        class(stream_master_stage), intent(inout) :: self
        character(len=*),         intent(in)    :: buffer
        if( .not. self%l_updates ) return
        if( self%l_stop_requested ) return
        if( .not. self%is_running() ) return
        if( .not. self%is_new_update(buffer) ) return
        call self%writer%send(buffer)
        self%last_update = buffer
    end subroutine send_update

    !> .true. when @p buffer differs from the last update sent to the stage since it started.
    logical function is_new_update( self, buffer )
        class(stream_master_stage), intent(in) :: self
        character(len=*),           intent(in) :: buffer
        is_new_update = .true.
        if( .not. allocated(self%last_update) ) return
        if( len(self%last_update) /= len(buffer) ) return
        is_new_update = self%last_update /= buffer
    end function is_new_update

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
        call self%fork%destroy()
        if( allocated(self%fork%commander) ) deallocate(self%fork%commander)
        if( allocated(self%last_update)    ) deallocate(self%last_update)
        self%l_updates        = .false.
        self%l_stop_requested = .false.
        self%l_exists         = .false.
    end subroutine kill

    !---------------- the GUI heartbeat ----------------

    !> A forked process as the GUI heartbeat reports it: its pid, its times and its status.
    function fork_gui_status( fork ) result( stage_status )
        class(forked_process), intent(inout) :: fork
        type(gui_stage_status) :: stage_status
        stage_status%pid       = int(fork%get_pid())
        stage_status%queuetime = fork%get_queuetime()
        stage_status%starttime = fork%get_starttime()
        stage_status%failtime  = fork%get_failtime()
        stage_status%stoptime  = fork%get_stoptime()
        select case( fork%status() )
            case( FORK_STATUS_RUNNING )
                stage_status%status = GUI_STAGE_STATUS_RUNNING
            case( FORK_STATUS_FAILED )
                stage_status%status = GUI_STAGE_STATUS_FAILED
            case( FORK_STATUS_STOPPED )
                stage_status%status = GUI_STAGE_STATUS_FINISHED
            case( FORK_STATUS_SKIPPED )
                stage_status%status = GUI_STAGE_STATUS_SKIPPED
            case DEFAULT
                stage_status%status = GUI_STAGE_STATUS_UNKNOWN
        end select
    end function fork_gui_status

end module simple_stream_master_stage
