!@descr: singleton for common state variables across the stream modules
module simple_stream2D_state
use simple_core_module_api
use simple_cmdline,            only: cmdline
use simple_qsys_env,           only: qsys_env
use simple_sp_project,         only: sp_project
use simple_stream_chunk,       only: stream_chunk
use simple_starproject,        only: starproject
implicit none

!===========================
! 1. Command-lines & queue
!===========================
class(cmdline), pointer :: master_cline => null()
type(cmdline)           :: cline_refine2D_chunk
type(cmdline)           :: cline_refine2D_pool

!===========================
! 2. Projects & dimensions
!===========================
type(sp_project), target        :: pool_proj
! the last POOL_NHISTORY completed iterations, for snapshots: a ring indexed by iteration
! (pool_history_slot); the same iterations' class-average and FRC files are kept on disk
integer, parameter              :: POOL_NHISTORY = 5
type(sp_project)                :: pool_proj_history(POOL_NHISTORY)
integer                         :: pool_history_iter(POOL_NHISTORY) = 0 ! the iteration each slot holds (0: none)
type(starproject)               :: starproj
type(scaled_dims)               :: chunk_dims
type(scaled_dims)               :: pool_dims
type(stream_chunk), allocatable :: chunks(:)
type(stream_chunk), allocatable :: converged_chunks(:)

!===========================
! 3. Global control flags
!===========================
logical :: l_no_chunks       = .false.
logical :: l_scaling         = .false.
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
! 8. Match classes rejection/selection
!===========================
integer, allocatable      :: match_selection(:)
logical                   :: l_match_selection_update = .false.

!===========================
! 9. Global filenames
!===========================
type(string)              :: refs_glob
type(string)              :: orig_projfile

contains

    !> The history slot of pool iteration @p iter (pool_proj_history, pool_history_iter).
    pure integer function pool_history_slot( iter )
        integer, intent(in) :: iter
        pool_history_slot = modulo(iter - 1, POOL_NHISTORY) + 1
    end function pool_history_slot

end module simple_stream2D_state
