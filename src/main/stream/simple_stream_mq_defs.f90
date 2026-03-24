!@descr: global POSIX message-queue handles for stream master IPC
!==============================================================================
! MODULE: simple_stream_mq_defs
!
! PURPOSE:
!   Declares the two shared ipc_mq handles used for inter-process
!   communication between the stream master and its worker processes.
!   Centralising them here avoids circular USE dependencies and ensures
!   both sides of each channel reference the same object.
!
! VARIABLES:
!   mq_stream_master_in  — inbound queue: messages sent TO the master
!   mq_stream_master_out — outbound queue: messages sent FROM the master
!
! DEPENDENCIES:
!   simple_ipc_mq
!==============================================================================
module simple_stream_mq_defs
  use simple_ipc_mq, only: ipc_mq

  implicit none

  type(ipc_mq), public :: mq_stream_master_in  ! master receives on this queue
  type(ipc_mq), public :: mq_stream_master_out ! master sends on this queue

end module simple_stream_mq_defs
