!@descr: status message of the persistent-worker protocol: the listener's idle or error reply
! Answers a heartbeat with no task (the worker ignores it), acknowledges queue_task (idle = queued,
! error = queue full) and answers unknown message types (error). Codes live in the server module.
! serialise() override: see simple_persistent_worker_message_base.
module simple_persistent_worker_message_status
    use simple_defs,                            only: STDLEN
    use simple_persistent_worker_message_base,  only: qsys_persistent_worker_message_base
    use simple_persistent_worker_message_types, only: WORKER_STATUS_MSG
    implicit none

    public  :: qsys_persistent_worker_message_status
    private

    !> Status wire message sent by the server to a persistent worker.
    !> Carries an idle/error acknowledgment when no task is available.
    type, extends(qsys_persistent_worker_message_base) :: qsys_persistent_worker_message_status
        integer               :: status  = 0   !< status code at time of transmission
        character(len=STDLEN) :: message = ''  !< human-readable status description
    contains
        procedure :: new       => new_qsys_persistent_worker_message_status       !< constructor
        procedure :: kill      => kill_qsys_persistent_worker_message_status      !< destructor
        procedure :: serialise => serialise_qsys_persistent_worker_message_status !< byte-buffer serialiser
    end type qsys_persistent_worker_message_status

contains

    !> Initialise a status message.
    !> Sets msg_type to WORKER_STATUS_MSG; all payload fields retain their
    !> default zero values and must be filled by the caller before transmission.
    subroutine new_qsys_persistent_worker_message_status( self )
        class(qsys_persistent_worker_message_status), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_STATUS_MSG
    end subroutine new_qsys_persistent_worker_message_status

    !> Reset all fields to their default zero / invalid state.
    !> No dynamic resources are held by this type; this is a plain field reset.
    subroutine kill_qsys_persistent_worker_message_status( self )
        class(qsys_persistent_worker_message_status), intent(inout) :: self
        self%msg_type = 0
        self%status   = 0
        self%message  = ''
    end subroutine kill_qsys_persistent_worker_message_status

    !> Serialise the full status message into an allocatable byte buffer.
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_persistent_worker_message_status layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    subroutine serialise_qsys_persistent_worker_message_status( self, buffer )
        class(qsys_persistent_worker_message_status), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_persistent_worker_message_status

end module simple_persistent_worker_message_status