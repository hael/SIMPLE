!@descr: settings of the stream commanders: the master (p00), the initial analysis run on its own (p03), multistate 3D (p07)
module simple_defs_stream_commanders
implicit none

! p00 master: what it gives the stages (its choices for a stream run from the GUI); their threads
! and parts are in simple_stream_master_resources
integer, parameter :: MASTER_POOL2D_NCLS            = 150
! p00 master: the master itself
integer, parameter :: MASTER_NCUNITS                = 16    ! computing units of the master's queue environment
integer, parameter :: MASTER_QSYS_NTHR              = 16
integer, parameter :: MASTER_HEARTBEAT_S            = 5     ! between heartbeats
integer, parameter :: MASTER_LAST_HEARTBEAT_TRIES   = 3     ! tries of the last heartbeat, MASTER_LAST_HEARTBEAT_WAIT_S apart
integer, parameter :: MASTER_LAST_HEARTBEAT_WAIT_S  = 1
integer, parameter :: MASTER_MEMLOG_EVERY           = 12    ! heartbeats between memory logs
integer, parameter :: MASTER_OPTICS_STOP_TIMEOUT_S  = 60    ! optics assignment stops first, waited for up to this
integer, parameter :: MASTER_STOP_TIMEOUT_S         = 600   ! a stage still running this long after it was asked to stop is killed
integer, parameter :: MASTER_STARTUP_STOP_TIMEOUT_S = 60    ! after a start-up failure, the stages are killed after this
integer, parameter :: MASTER_LISTENER_SLEEP_US      = 10000 ! the listener's pause between drains

! p03 initial analysis: persistent workers for the stage's jobs
integer, parameter :: P03_NPARTS2D = 8
! p03: the defaults when the master gives none (the stage run on its own; simple_stream_master_resources)
integer, parameter :: P03_NTHR2D   = 16 ! threads of the 2D jobs and of the sieve's chunks
integer, parameter :: P03_NCHUNKS  = 4  ! the sieve's chunks at once
! p03: the 3D route to the picking references (stream fix plan, decision 20)
integer, parameter :: NSTATES_PICKREFS   = 3
integer, parameter :: NSTAGES_PICKREFS   = 3
real,    parameter :: LPSTOP_PICKREFS    = 8.   ! A
integer, parameter :: NSPACE_PICKREFS    = 50   ! reprojections
integer, parameter :: NTHR3D_PICKREFS    = 16   ! threads of the 3D job and the reprojection
integer, parameter :: NRESTARTS_COLLAPSE = 3
real,    parameter :: LPSTART_INI3D      = 100. ! A
real,    parameter :: LPSTOP_INI3D       = 20.  ! A

! p07 multistate 3D
integer, parameter :: P07_NSTATES  = 3 ! states of the reference volumes at most
integer, parameter :: P07_NPARTS3D = 8 ! parts and threads of each 3D job
integer, parameter :: P07_NTHR3D   = 8

end module simple_defs_stream_commanders
