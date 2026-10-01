!@descr: terminate message of the persistent-worker protocol: the server orders a worker to shut down
! Sent in reply to a heartbeat on server shutdown, scale-down, or an out-of-range worker_id/UID clash.
! The worker then cancels its running tasks and exits. The listener's own kill sentinel uses the same code.
! serialise() override: see simple_persistent_worker_message_base.
module simple_persistent_worker_message_terminate
    use simple_defs,                            only: STDLEN
    use simple_persistent_worker_message_base,  only: qsys_persistent_worker_message_base
    use simple_persistent_worker_message_types, only: WORKER_TERMINATE_MSG
    implicit none

    public  :: qsys_persistent_worker_message_terminate
    private

    !> Terminate wire message sent by the server to a persistent worker to command shutdown.
    !> Carries the shutdown timestamp and a human-readable reason; the worker
    !> cancels in-flight tasks and exits upon receipt.
    type, extends(qsys_persistent_worker_message_base) :: qsys_persistent_worker_message_terminate
        integer               :: terminate_time = 0   !< UNIX timestamp at time of transmission
        character(len=STDLEN) :: reason         = ''  !< human-readable terminate reason
    contains
        procedure :: new       => new_qsys_persistent_worker_message_terminate       !< constructor
        procedure :: kill      => kill_qsys_persistent_worker_message_terminate      !< destructor
        procedure :: serialise => serialise_qsys_persistent_worker_message_terminate !< byte-buffer serialiser
    end type qsys_persistent_worker_message_terminate

contains

    !> Initialise a terminate message.
    !> Sets msg_type to WORKER_TERMINATE_MSG; all payload fields retain their
    !> default zero values and must be filled by the caller before transmission.
    subroutine new_qsys_persistent_worker_message_terminate( self )
        class(qsys_persistent_worker_message_terminate), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_TERMINATE_MSG
    end subroutine new_qsys_persistent_worker_message_terminate

    !> Reset all fields to their default zero / invalid state.
    !> No dynamic resources are held by this type; this is a plain field reset.
    subroutine kill_qsys_persistent_worker_message_terminate( self )
        class(qsys_persistent_worker_message_terminate), intent(inout) :: self
        self%msg_type         = 0
        self%terminate_time   = 0
        self%reason           = ''
    end subroutine kill_qsys_persistent_worker_message_terminate

    !> Serialise the full terminate message into an allocatable byte buffer.
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_persistent_worker_message_terminate layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    subroutine serialise_qsys_persistent_worker_message_terminate( self, buffer )
        class(qsys_persistent_worker_message_terminate), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_persistent_worker_message_terminate

end module simple_persistent_worker_message_terminate