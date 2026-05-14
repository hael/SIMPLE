!==============================================================================
! MODULE: simple_persistent_worker_message_task
!
! PURPOSE:
!   Defines the task wire message exchanged between the persistent-worker
!   server and workers.  The same payload is used for:
!     - queued task requests (master -> listener)
!     - dispatched tasks     (listener -> worker)
!
!   When a persistent worker sends a heartbeat and the server has a queued job that fits
!   the worker's available thread capacity, the server replies with a task
!   message.  The task message carries:
!     - job_id       — unique job counter (>0 when active)
!     - queue_time   — UNIX time of enqueue
!     - start_time   — UNIX time worker started execution
!     - end_time     — UNIX time of completion (0 = pending)
!     - exit_code    — script exit status after completion
!     - nthr         — thread slots required
!     - submitted    — .true. once dispatched to a persistent worker
!     - priority     — .true. when caller requests priority handling
!     - script_path  — absolute path to the bash script
!
! DESIGN CONTRACT:
!   serialise_qsys_persistent_worker_message_task inlines the three-statement
!   TRANSFER body (deallocate / allocate(sizeof(self)) / transfer) rather
!   than delegating to the base serialise procedure.  This is mandatory:
!   sizeof() is resolved against the declared type of the dummy argument, so
!   calling the base procedure would allocate only sizeof(qsys_persistent_worker_message_base)
!   bytes and silently truncate the task-specific fields from the buffer.
!   See simple_persistent_worker_message_base for the full design contract.
!
! USAGE:
!   type(qsys_persistent_worker_message_task) :: task
!   call task%new()
!   task%job_id       = my_job_id
!   task%queue_time   = int(time())
!   task%start_time   = 0
!   task%end_time     = 0
!   task%exit_code    = 0
!   task%nthr         = required_threads
!   task%submitted    = .false.
!   task%priority     = .false.
!   task%script_path  = script_path
!   call task%serialise(buffer)
!
! DEPENDENCIES:
!   simple_persistent_worker_message_base  — base type and serialise contract
!   simple_persistent_worker_message_types — WORKER_TASK_MSG enumerator
!==============================================================================
module simple_persistent_worker_message_task
    use simple_defs,                            only: STDLEN
    use simple_persistent_worker_message_base,  only: qsys_persistent_worker_message_base
    use simple_persistent_worker_message_types, only: WORKER_TASK_MSG
    implicit none

    public  :: qsys_persistent_worker_message_task
    private

    !> Task wire payload for queueing and dispatch.
    !> Carries job identity, lifecycle timestamps, scheduling hints,
    !> required thread count, and executable script path.
    type, extends(qsys_persistent_worker_message_base) :: qsys_persistent_worker_message_task
        integer :: job_id     = 0                  !< unique job counter (>0 when active)
        integer :: queue_time = 0                  !< UNIX time of enqueue
        integer :: start_time = 0                  !< UNIX time worker started execution
        integer :: end_time   = 0                  !< UNIX time of completion (0 = pending)
        integer :: exit_code  = 0                  !< script exit status after completion
        integer :: nthr       = 0                  !< thread slots required
        logical :: submitted  = .false.            !< .true. once dispatched to a persistent worker
        logical :: priority   = .false.            !< .true. if this task should be prioritised over others in the queue
        character(len=STDLEN) :: script_path = ''  !< absolute path to the bash script
    contains
        procedure :: new       => new_qsys_persistent_worker_message_task       !< constructor
        procedure :: kill      => kill_qsys_persistent_worker_message_task      !< destructor
        procedure :: serialise => serialise_qsys_persistent_worker_message_task !< byte-buffer serialiser
    end type qsys_persistent_worker_message_task

contains

    !> Initialise a task message to a clean default state.
    !> Sets msg_type to WORKER_TASK_MSG and clears all payload fields.
    subroutine new_qsys_persistent_worker_message_task( self )
        class(qsys_persistent_worker_message_task), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_TASK_MSG
    end subroutine new_qsys_persistent_worker_message_task

    !> Reset all fields to their default zero / invalid state.
    !> This type owns no dynamic resources; kill() is a pure field reset.
    subroutine kill_qsys_persistent_worker_message_task( self )
        class(qsys_persistent_worker_message_task), intent(inout) :: self
        self%msg_type       = 0
        self%job_id         = 0
        self%queue_time     = 0
        self%start_time     = 0
        self%end_time       = 0
        self%exit_code      = 0
        self%nthr           = 0
        self%submitted      = .false.
        self%priority       = .false.
        self%script_path    = ''
    end subroutine kill_qsys_persistent_worker_message_task

    !> Serialise the full task message into an allocatable byte buffer.
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_persistent_worker_message_task layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    subroutine serialise_qsys_persistent_worker_message_task( self, buffer )
        class(qsys_persistent_worker_message_task), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_persistent_worker_message_task

end module simple_persistent_worker_message_task