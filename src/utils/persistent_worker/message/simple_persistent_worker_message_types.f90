!@descr: message-type enumeration of the persistent-worker wire protocol, shared by server and workers
! The leading integer of every message holds one of these values; receivers dispatch on it.
! The values are fixed: changing one breaks the protocol. Code 3 never goes on the wire.
module simple_persistent_worker_message_types
    implicit none

    ! Fortran enum constants have no access specifiers; all are implicitly public.
    enum, bind(c)
        enumerator :: WORKER_TERMINATE_MSG = 1  !< server→worker: shut down; also the kill sentinel sent to the listener (value 1)
        enumerator :: WORKER_HEARTBEAT_MSG      !< worker→server: liveness + thread-load info  (value 2)
        enumerator :: WORKER_TASK_MSG           !< set by task%new() only; queue_task overwrites it (value 3)
        enumerator :: WORKER_NEW_TASK_MSG       !< queue_task→listener, and listener→worker dispatch (value 4)
        enumerator :: WORKER_STATUS_MSG         !< listener reply: idle (no task, or task queued) or error (value 5)
    end enum

end module simple_persistent_worker_message_types