!@descr: heartbeat message of the persistent-worker protocol: a worker reports liveness and thread load
! Worker -> server each loop: worker_id, worker_uid, heartbeat_time, nthr_used/total (fd set server-side).
! Reply: task, STATUS idle, or TERMINATE.
! serialise() override: see simple_persistent_worker_message_base.
module simple_persistent_worker_message_heartbeat
    use simple_persistent_worker_message_base,  only: qsys_persistent_worker_message_base
    use simple_persistent_worker_message_types, only: WORKER_HEARTBEAT_MSG
    implicit none

    public  :: qsys_persistent_worker_message_heartbeat
    private

    !> Heartbeat wire message sent by a persistent worker to the server.
    !> Carries the worker identity, a liveness timestamp, and thread-load
    !> information used by the server to assign tasks appropriately.
    type, extends(qsys_persistent_worker_message_base) :: qsys_persistent_worker_message_heartbeat
        integer :: worker_id      = 0  !< 1-based worker slot index on the server
        integer :: heartbeat_time = 0  !< UNIX timestamp at time of transmission
        integer :: nthr_used      = 0  !< threads currently executing tasks
        integer :: nthr_total     = 0  !< total thread capacity of this worker
        integer :: fd             = 0  !< server-side connection fd, overwritten on receipt; used to clear the registry entry on disconnect
        character(len=256) :: worker_uid = ''  !< unique worker identifier: <hostname>_<PID>
    contains
        procedure :: new       => new_qsys_persistent_worker_message_heartbeat       !< constructor
        procedure :: kill      => kill_qsys_persistent_worker_message_heartbeat      !< destructor
        procedure :: serialise => serialise_qsys_persistent_worker_message_heartbeat !< byte-buffer serialiser
    end type qsys_persistent_worker_message_heartbeat

contains

    !> Initialise a heartbeat message.
    !> Sets msg_type to WORKER_HEARTBEAT_MSG; all payload fields retain their
    !> default zero values and must be filled by the caller before transmission.
    subroutine new_qsys_persistent_worker_message_heartbeat( self )
        class(qsys_persistent_worker_message_heartbeat), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_HEARTBEAT_MSG
    end subroutine new_qsys_persistent_worker_message_heartbeat

    !> Reset all fields to their default zero / invalid state.
    !> No dynamic resources are held by this type; this is a plain field reset.
    subroutine kill_qsys_persistent_worker_message_heartbeat( self )
        class(qsys_persistent_worker_message_heartbeat), intent(inout) :: self
        self%msg_type       = 0
        self%worker_id      = 0
        self%heartbeat_time = 0
        self%nthr_used      = 0
        self%nthr_total     = 0
        self%fd             = 0
        self%worker_uid     = ''
    end subroutine kill_qsys_persistent_worker_message_heartbeat

    !> Serialise the full heartbeat message into an allocatable byte buffer.
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_persistent_worker_message_heartbeat layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    subroutine serialise_qsys_persistent_worker_message_heartbeat( self, buffer )
        class(qsys_persistent_worker_message_heartbeat), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_persistent_worker_message_heartbeat

end module simple_persistent_worker_message_heartbeat