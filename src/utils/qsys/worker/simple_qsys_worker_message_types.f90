!==============================================================================
! MODULE: simple_qsys_worker_message_types
!
! PURPOSE:
!   Defines the C-interoperable wire-protocol enumeration shared by the
!   SIMPLE queue-system server and all worker processes.  The leading
!   integer of every message struct holds one of these constants; the
!   receiver dispatches on this value to determine how to interpret the
!   remainder of the buffer.
!
! WIRE VALUES (fixed; changing any value breaks the protocol):
!   1 — WORKER_TERMINATE_MSG  server→worker  command orderly shutdown
!   2 — WORKER_HEARTBEAT_MSG  worker→server  liveness + thread-load report
!   3 — WORKER_TASK_MSG       server→worker  dispatch script for execution
!   4 — WORKER_STATUS_MSG     server→worker  idle acknowledgment (no task)
!
! DEPENDENCIES:
!   None.
!==============================================================================
module simple_qsys_worker_message_types
    implicit none

    ! Fortran enum constants have no access specifiers; all are implicitly public.
    enum, bind(c)
        enumerator :: WORKER_TERMINATE_MSG = 1  !< server→worker: command orderly shutdown     (value 1)
        enumerator :: WORKER_HEARTBEAT_MSG      !< worker→server: liveness + thread-load info  (value 2)
        enumerator :: WORKER_TASK_MSG           !< server→worker: dispatch script for execution (value 3)
        enumerator :: WORKER_STATUS_MSG         !< server→worker: idle acknowledgment, no task  (value 4)
    end enum

end module simple_qsys_worker_message_types