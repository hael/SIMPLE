!@descr: global stream master pipe descriptors for IPC
! Kept here so the master and the forked stages share the same descriptors without circular USE.
! Names are from the master's side: ipc_pipe_*_in carries stage->master metadata,
! ipc_pipe_*_out carries master->stage updates.
module simple_stream_state

  implicit none

  integer,      public :: ipc_pipe_preprocess_in(2)            = [-1, -1]  ! pipe for preprocessing ipc
  integer,      public :: ipc_pipe_preprocess_out(2)           = [-1, -1]  ! pipe for preprocessing ipc
  integer,      public :: ipc_pipe_assign_optics_in(2)         = [-1, -1]  ! pipe for assign_optics ipc
  integer,      public :: ipc_pipe_assign_optics_out(2)        = [-1, -1]  ! pipe for assign_optics ipc
  integer,      public :: ipc_pipe_initial_analysis_in(2)      = [-1, -1]  ! pipe for initial_analysis ipc
  integer,      public :: ipc_pipe_initial_analysis_out(2)     = [-1, -1]  ! pipe for initial_analysis ipc
  integer,      public :: ipc_pipe_refpick_in(2)               = [-1, -1]  ! pipe for refpick ipc
  integer,      public :: ipc_pipe_refpick_out(2)              = [-1, -1]  ! pipe for refpick ipc
  integer,      public :: ipc_pipe_sieve_cavgs_in(2)           = [-1, -1]  ! pipe for sieve_cavgs ipc
  integer,      public :: ipc_pipe_sieve_cavgs_out(2)          = [-1, -1]  ! pipe for sieve_cavgs ipc
  integer,      public :: ipc_pipe_pool2D_in(2)                = [-1, -1]  ! pipe for pool2D ipc
  integer,      public :: ipc_pipe_pool2D_out(2)               = [-1, -1]  ! pipe for pool2D ipc
  integer,      public :: ipc_pipe_abinitio3D_multstate_in(2)  = [-1, -1]  ! pipe for 3D multistate ipc
  integer,      public :: ipc_pipe_abinitio3D_multstate_out(2) = [-1, -1]  ! pipe for 3D multistate ipc

end module simple_stream_state
