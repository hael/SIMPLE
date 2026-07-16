!@descr: global stream master pipe descriptors for IPC
!==============================================================================
! MODULE: simple_stream_state
!
! PURPOSE:
!   Declares the shared pipe descriptor arrays used for inter-process
!   communication between the stream master and worker processes.
!   Centralizing them here avoids circular USE dependencies and ensures
!   both sides of each channel reference the same descriptors.
!
! VARIABLES:
!   ipc_pipe_*_in  — stage input pipes: commands sent FROM master TO stage
!   ipc_pipe_*_out — stage output pipes: metadata sent FROM stage TO master
!
! DEPENDENCIES:
!   none
!==============================================================================
module simple_stream_state

  implicit none

  integer,      public :: ipc_pipe_preprocess_in(2)         = [-1, -1]  ! pipe for preprocessing ipc
  integer,      public :: ipc_pipe_preprocess_out(2)        = [-1, -1]  ! pipe for preprocessing ipc
  integer,      public :: ipc_pipe_assign_optics_in(2)      = [-1, -1]  ! pipe for assign_optics ipc
  integer,      public :: ipc_pipe_assign_optics_out(2)     = [-1, -1]  ! pipe for assign_optics ipc
  integer,      public :: ipc_pipe_initial_analysis_in(2)   = [-1, -1]  ! pipe for initial_analysis ipc
  integer,      public :: ipc_pipe_initial_analysis_out(2)  = [-1, -1]  ! pipe for initial_analysis ipc
  integer,      public :: ipc_pipe_refpick_in(2)            = [-1, -1]  ! pipe for refpick ipc
  integer,      public :: ipc_pipe_refpick_out(2)           = [-1, -1]  ! pipe for refpick ipc
  integer,      public :: ipc_pipe_sieve_cavgs_in(2)        = [-1, -1]  ! pipe for sieve_cavgs ipc
  integer,      public :: ipc_pipe_sieve_cavgs_out(2)       = [-1, -1]  ! pipe for sieve_cavgs ipc
  integer,      public :: ipc_pipe_pool2D_in(2)             = [-1, -1]  ! pipe for pool2D ipc
  integer,      public :: ipc_pipe_pool2D_out(2)            = [-1, -1]  ! pipe for pool2D ipc

end module simple_stream_state
