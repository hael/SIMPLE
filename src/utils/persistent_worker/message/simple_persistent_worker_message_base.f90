!@descr: polymorphic base type of the persistent-worker wire messages
! msg_type leads every message; 0 = uninitialised. Subtypes set it in new().
! Messages travel as their raw memory image. Each subtype overrides serialise() with its own
! copy of the transfer body; do not send extended types through this routine.
module simple_persistent_worker_message_base
    implicit none

    public  :: qsys_persistent_worker_message_base
    private

    !> Base type for all persistent-worker wire messages.
    !> Holds only the leading msg_type discriminator that all messages share.
    type :: qsys_persistent_worker_message_base
        integer :: msg_type = 0  !< wire message type; 0 = uninitialised/invalid
    contains
        procedure :: new         => new_qsys_persistent_worker_message_base  !< default constructor
        procedure :: kill        => kill_qsys_persistent_worker_message_base !< default destructor
        procedure :: serialise   !< serialise to raw byte buffer (must override in subtypes)
    end type qsys_persistent_worker_message_base

contains

    !> Default constructor; intentional no-op.
    !> Concrete subtypes must set msg_type in their own new() override.
    subroutine new_qsys_persistent_worker_message_base( self )
        class(qsys_persistent_worker_message_base), intent(inout) :: self
        ! no-op: msg_type is set by the concrete subtype constructor
    end subroutine new_qsys_persistent_worker_message_base

    !> Default destructor; intentional no-op.
    !> Override in concrete subtypes that hold dynamic resources.
    subroutine kill_qsys_persistent_worker_message_base( self )
        class(qsys_persistent_worker_message_base), intent(inout) :: self
        ! no-op: base type holds no dynamic resources
    end subroutine kill_qsys_persistent_worker_message_base

    !> Copy the raw memory image of self into buffer (sizeof(self) bytes) via TRANSFER.
    !> Subtypes override this with their own copy; do not send extended types through it.
    subroutine serialise( self, buffer )
        class(qsys_persistent_worker_message_base), intent(in)    :: self
        character(len=:), allocatable,   intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise

end module simple_persistent_worker_message_base