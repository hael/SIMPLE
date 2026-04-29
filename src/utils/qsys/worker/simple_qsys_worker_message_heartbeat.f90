!==============================================================================
! MODULE: simple_qsys_worker_message_heartbeat
!
! PURPOSE:
!   Provides the concrete heartbeat message type used by SIMPLE worker
!   processes to signal liveness and thread-capacity information to the
!   queue-system server.
!
!   Workers send a heartbeat message at regular intervals (HEARTBEAT_TIMEOUT_MS).
!   The server replies with either a task message or a TERMINATE message.
!   The heartbeat carries:
!     - worker_id       — 1-based slot index on the server side
!     - heartbeat_time  — UNIX timestamp at the moment of transmission
!     - nthr_used       — threads currently executing tasks
!     - nthr_total      — total thread capacity advertised by this worker
!     - worker_uid      — unique worker identifier: <hostname>_<PID>
!
! DESIGN CONTRACT:
!   serialise_qsys_worker_message_heartbeat inlines the three-statement
!   TRANSFER body (deallocate / allocate(sizeof(self)) / transfer) rather
!   than delegating to the base serialise procedure.  This is mandatory:
!   sizeof() is resolved against the declared type of the dummy argument, so
!   calling the base procedure would allocate only sizeof(qsys_worker_message_base)
!   bytes and silently truncate the four heartbeat fields from the buffer.
!   See simple_qsys_worker_message_base for the full design contract.
!
! USAGE:
!   type(qsys_worker_message_heartbeat) :: hb
!   call hb%new()
!   hb%worker_id      = my_id
!   hb%heartbeat_time = int(time())
!   hb%nthr_used      = get_nthr_used()
!   hb%nthr_total     = nthr
!   call hb%serialise(buffer)
!
! DEPENDENCIES:
!   simple_qsys_worker_message_base  — base type and serialise contract
!   simple_qsys_worker_message_types — WORKER_HEARTBEAT_MSG enumerator
!==============================================================================
module simple_qsys_worker_message_heartbeat
    use simple_qsys_worker_message_base,  only: qsys_worker_message_base
    use simple_qsys_worker_message_types, only: WORKER_HEARTBEAT_MSG
    implicit none

    public  :: qsys_worker_message_heartbeat
    private

    !> Heartbeat wire message sent by a worker to the queue-system server.
    !> Carries the worker identity, a liveness timestamp, and thread-load
    !> information used by the server to assign tasks appropriately.
    type, extends(qsys_worker_message_base) :: qsys_worker_message_heartbeat
        integer :: worker_id      = 0  !< 1-based worker slot index on the server
        integer :: heartbeat_time = 0  !< UNIX timestamp at time of transmission
        integer :: nthr_used      = 0  !< threads currently executing tasks
        integer :: nthr_total     = 0  !< total thread capacity of this worker
        character(len=256) :: worker_uid = ''  !< unique worker identifier: <hostname>_<PID>
    contains
        procedure :: new       => new_qsys_worker_message_heartbeat       !< constructor
        procedure :: kill      => kill_qsys_worker_message_heartbeat      !< destructor
        procedure :: serialise => serialise_qsys_worker_message_heartbeat !< byte-buffer serialiser
    end type qsys_worker_message_heartbeat

contains

    !> Initialise a heartbeat message.
    !> Sets msg_type to WORKER_HEARTBEAT_MSG; all payload fields retain their
    !> default zero values and must be filled by the caller before transmission.
    subroutine new_qsys_worker_message_heartbeat( self )
        class(qsys_worker_message_heartbeat), intent(inout) :: self
        call self%kill()
        self%msg_type = WORKER_HEARTBEAT_MSG
    end subroutine new_qsys_worker_message_heartbeat

    !> Reset all fields to their default zero / invalid state.
    !> No dynamic resources are held by this type; this is a plain field reset.
    subroutine kill_qsys_worker_message_heartbeat( self )
        class(qsys_worker_message_heartbeat), intent(inout) :: self
        self%msg_type       = 0
        self%worker_id      = 0
        self%heartbeat_time = 0
        self%nthr_used      = 0
        self%nthr_total     = 0
        self%worker_uid     = ''
    end subroutine kill_qsys_worker_message_heartbeat

    !> Serialise the full heartbeat message into an allocatable byte buffer.
    !>
    !> The buffer is (re)allocated to exactly sizeof(self) bytes — the size of
    !> the complete qsys_worker_message_heartbeat layout — and filled with the
    !> raw storage representation of self via TRANSFER.
    !>
    !> WARNING: This body must NOT be replaced with a call to the base serialise
    !> procedure.  sizeof() is resolved against the declared type of the dummy
    !> argument; delegating to the base would yield sizeof(qsys_worker_message_base)
    !> and silently drop worker_id, heartbeat_time, nthr_used, nthr_total,
    !> and worker_uid from the buffer.
    subroutine serialise_qsys_worker_message_heartbeat( self, buffer )
        class(qsys_worker_message_heartbeat), intent(in)    :: self
        character(len=:), allocatable,        intent(inout) :: buffer
        if( allocated(buffer) ) deallocate(buffer)
        allocate(character(len=sizeof(self)) :: buffer)
        buffer = transfer(self, buffer)
    end subroutine serialise_qsys_worker_message_heartbeat

end module simple_qsys_worker_message_heartbeat