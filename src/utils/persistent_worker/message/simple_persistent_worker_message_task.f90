!@descr: task message of the persistent-worker protocol: a queued job request and the job dispatched to a worker
! queue_task() sends it to the listener as WORKER_NEW_TASK_MSG; the listener assigns job_id and later
! dispatches the same record, msg_type unchanged, on a heartbeat with enough free threads.
! queue_time, start_time, end_time and exit_code are never set.
! serialise() override: see simple_persistent_worker_message_base.
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