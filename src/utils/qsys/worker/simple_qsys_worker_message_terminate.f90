!==============================================================================
! MODULE: simple_qsys_worker_message_terminate
!
! PURPOSE:
!   Provides the concrete terminate message type sent by the queue-system server
!   to worker processes to command an orderly shutdown.
!
!   When the server is stopping it replies to worker heartbeats with a terminate
!   message.  Upon receipt, the worker finishes its current tasks and exits.
!   The terminate message carries:
!     - terminate_time  — UNIX timestamp at the moment of transmission
!     - reason            — human-readable terminate reason
!
! DESIGN CONTRACT:
!   serialise_qsys_worker_message_terminate inlines the three-statement
!   TRANSFER body (deallocate / allocate(sizeof(self)) / transfer) rather
!   than delegating to the base serialise procedure.  This is mandatory:
!   sizeof() is resolved against the declared type of the dummy argument, so
!   calling the base procedure would allocate only sizeof(qsys_worker_message_base)
!   bytes and silently truncate the terminate_time and reason fields from the buffer.
!   See simple_qsys_worker_message_base for the full design contract.
!
! USAGE:
!   type(qsys_worker_message_terminate) :: term
!   call term%new()
!   term%terminate_time = int(time())
!   term%reason         = 'Server shutting down'
!   call term%serialise(buffer)
!
! DEPENDENCIES:
!   simple_qsys_worker_message_base  — base type and serialise contract
!   simple_qsys_worker_message_types — WORKER_TERMINATE_MSG enumerator
!==============================================================================
module simple_qsys_worker_message_terminate
    use simple_defs,                      only: STDLEN
    use simple_qsys_worker_message_base,  only: qsys_worker_message_base
    use simple_qsys_worker_message_types, only: WORKER_TERMINATE_MSG
    implicit none

    public  :: qsys_worker_message_terminate
    private

    !> Terminate wire message sent by the server to a worker to command shutdown.
    !> Carries the shutdown timestamp and a human-readable reason; the worker
    !> finishes in-flight tasks and exits upon receipt.
    type, extends(qsys_worker_message_base) :: qsys_worker_message_terminate
        integer               :: terminate_time = 0   !< UNIX timestamp at time of transmission
        character(len=STDLEN) :: reason         = ''  !< human-readable terminate reason
    contains
        procedure :: new       => new_qsys_worker_message_terminate       !< constructor
        procedure :: kill      => kill_qsys_worker_message_terminate      !< destructor
        procedure :: serialise => serialise_qsys_worker_message_terminate !< byte-buffer serialiser
    end type qsys_worker_message_terminate

contains

    !> Initialise a terminate message.
    !> Sets msg_type to WORKER_TERMINATE_MSG; all payload fields retain their
    !> default zero values and must be filled by the caller before transmission.
    subroutine new_qsys_worker_message_terminate( self )
        class(qsys_worker_message_terminate), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_TERMINATE_MSG
    end subroutine new_qsys_worker_message_terminate

    !> Reset all fields to their default zero / invalid state.
    !> No dynamic resources are held by this type; this is a plain field reset.
    subroutine kill_qsys_worker_message_terminate( self )
        class(qsys_worker_message_terminate), intent(inout) :: self
        self%msg_type         = 0
        self%terminate_time   = 0
        self%reason           = ''
    end subroutine kill_qsys_worker_message_terminate

    !> Serialise the full terminate message into an allocatable byte buffer.
    !>
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_worker_message_terminate layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    !>
    !> WARNING: This body must NOT be replaced with a call to the base serialise
    !> procedure.  sizeof() is resolved against the declared type of the dummy
    !> argument; delegating to the base would yield sizeof(qsys_worker_message_base)
    !> and silently drop terminate_time and reason from the buffer.
    subroutine serialise_qsys_worker_message_terminate( self, buffer )
        class(qsys_worker_message_terminate), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_worker_message_terminate

end module simple_qsys_worker_message_terminate