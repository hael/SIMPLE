module simple_stream2D_state
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_qsys_env,           only: qsys_env
use simple_sp_project,         only: sp_project
use simple_stream_chunk,       only: stream_chunk
use simple_starproject,        only: starproject
use simple_starproject_stream, only: starproject_stream
use json_module,               only: json_value
implicit none

!===========================
! 1. Command-lines & queue
!===========================
class(cmdline), pointer :: master_cline => null()
type(cmdline)           :: cline_cluster2D_chunk
type(cmdline)           :: cline_cluster2D_pool

!===========================
! 2. Projects & dimensions
!===========================
type(sp_project), target        :: pool_proj
type(sp_project), allocatable   :: pool_proj_history(:)
type(starproject)               :: starproj
type(scaled_dims)               :: chunk_dims
type(scaled_dims)               :: pool_dims
type(stream_chunk), allocatable :: chunks(:)
type(stream_chunk), allocatable :: converged_chunks(:)

!===========================
! 3. Global control flags
!===========================
logical :: l_abinitio2D      = .false.
logical :: l_no_chunks       = .false.
logical :: l_scaling         = .false.
logical :: l_update_sigmas   = .false.
logical :: l_wfilt           = .false.
logical :: l_stream2D_active = .false.
logical :: l_pool_available  = .false.

!===========================
! 4. Global counters / bookkeeping
!===========================
integer :: pool_iter            = 0
integer :: glob_chunk_id        = 0
integer :: ncls_glob            = 0
integer :: nptcls_per_chunk     = 0
integer :: last_complete_iter   = 0
integer :: numlen               = 0

!===========================
! 5. Resolution / masks
!===========================
real             :: lpstart = 0.0
real             :: lpstop  = 0.0
real             :: lpcen   = 0.0
character(len=6) :: lpthres_type = ""  ! "auto"/"manual"/"off"

!===========================
! 6. GUI / JPEG / stats
!===========================
integer, allocatable      :: pool_jpeg_map(:)
integer, allocatable      :: pool_jpeg_pop(:)
real,    allocatable      :: pool_jpeg_res(:)
type(string)              :: projfile4gui

!===========================
! 7. Snapshot / interactive
!===========================
integer                   :: snapshot_iteration   = 0
integer                   :: snapshot_last_nptcls = 0
integer, allocatable      :: snapshot_selection(:)
type(json_value), pointer :: snapshot_json => null()

!===========================
! 8. Global filenames
!===========================
type(string)              :: refs_glob
type(string)              :: orig_projfile
end module simple_stream2D_state
