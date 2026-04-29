!==============================================================================
! MODULE: simple_qsys_worker_message_base
!
! PURPOSE:
!   Provides the polymorphic base type for all SIMPLE worker wire messages.
!   Every concrete message type (heartbeat, task, status, terminate) extends
!   qsys_worker_message_base and inherits the serialise method.
!
! DESIGN CONTRACT:
!   1.  msg_type MUST be set to the correct WORKER_*_MSG enumerator constant
!       by each concrete subtype constructor or default initialiser.
!       The base default of 0 is intentionally invalid so an uninitialised
!       message is identifiable on the receiving end.
!   2.  serialise() uses sizeof() — a non-standard GCC extension equivalent
!       to storage_size(self)/8 — to determine the byte count of the dynamic
!       type.  Each derived type MUST override serialise() with the same body
!       so that sizeof() is called with the correct declared type, ensuring
!       the full extent of the derived-type layout is transferred.
!       Failure to override will silently truncate larger derived types to
!       the base-type size.
!
! USAGE:
!   type, extends(qsys_worker_message_base) :: my_message
!       integer :: msg_type = MY_MSG_TYPE   ! override default 0
!       ...
!   contains
!       procedure :: serialise => serialise_my_message   ! required override
!   end type my_message
!
! DEPENDENCIES:
!   None (no module-level USE required; sizeof is a compiler built-in).
!==============================================================================
module simple_qsys_worker_message_base
    implicit none

    public  :: qsys_worker_message_base
    private

    !> Abstract base for all worker wire messages.
    !> Holds only the leading msg_type discriminator that all messages share.
    type :: qsys_worker_message_base
        integer :: msg_type = 0  !< wire message type; 0 = uninitialised/invalid
    contains
        procedure :: new         => new_qsys_worker_message_base  !< default constructor
        procedure :: kill        => kill_qsys_worker_message_base !< default destructor
        procedure :: serialise   !< serialise to raw byte buffer (override in subtypes)
    end type qsys_worker_message_base

contains

    !> Default constructor; intentional no-op.
    !> Concrete subtypes must set msg_type in their own new() override.
    subroutine new_qsys_worker_message_base( self )
        class(qsys_worker_message_base), intent(inout) :: self
        ! no-op: msg_type is set by the concrete subtype constructor
    end subroutine new_qsys_worker_message_base

    !> Default destructor; intentional no-op.
    !> Override in concrete subtypes that hold dynamic resources.
    subroutine kill_qsys_worker_message_base( self )
        class(qsys_worker_message_base), intent(inout) :: self
        ! no-op: base type holds no dynamic resources
    end subroutine kill_qsys_worker_message_base

    !> Serialise the message into an allocatable byte buffer using TRANSFER.
    !>
    !> The buffer is (re)allocated to exactly sizeof(self) bytes and filled
    !> with the raw storage representation of self.  Because sizeof() is
    !> evaluated against the declared type of the dummy argument (the base
    !> type), derived-type instances that do NOT override this procedure will
    !> produce a truncated buffer.  Always override in concrete subtypes.
    subroutine serialise( self, buffer )
        class(qsys_worker_message_base), intent(in)    :: self
        character(len=:), allocatable,   intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise

end module simple_qsys_worker_message_base