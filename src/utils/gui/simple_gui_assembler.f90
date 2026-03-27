!@descr: Assembles GUI metadata objects into a compact JSON document and sends it to the NICE frontend
!==============================================================================
! MODULE: simple_gui_assembler
!
! PURPOSE:
!   Owns a single json-fortran document tree (json_root) that is populated
!   incrementally by the stream pipeline.  Each assemble_* routine updates one
!   named section of the tree; unchanged sections are suppressed via FNV-1a
!   content hashing so only deltas are transmitted.  The assembled document is
!   serialised to a string with to_string() for dispatch over IPC.
!
! TYPES:
!   gui_assembler
!     new()                          — initialise the JSON tree for a job
!     kill()                         — destroy the JSON tree and reset state
!     to_string()                    — serialise the current tree to a string
!     set_stoptime()                 — record the job stop timestamp
!     clear_hashes()                 — reset change-detection hashes
!     is_associated()                — .true. if the JSON root is initialised
!     assemble_stream_heartbeat()    — write process-status section
!     assemble_stream_preprocess()   — write preprocessing section
!     assemble_stream_optics_assignment() — write optics-assignment section
!     assemble_stream_initial_picking()    — write initial-picking section
!     assemble_stream_reference_picking()  — write reference-picking section
!     assemble_stream_opening2D()          — write 2D-classification section
!
! DEPENDENCIES:
!   unix, simple_string, simple_forked_process, simple_gui_metadata_api
!==============================================================================
module simple_gui_assembler
  use unix,                    only: c_time
  use simple_string,           only: string
  use simple_forked_process,   only: forked_process,       &
                                     FORK_STATUS_RUNNING,   &
                                     FORK_STATUS_FAILED,    &
                                     FORK_STATUS_STOPPED,   &
                                     FORK_STATUS_RESTARTING,&
                                     FORK_STATUS_SKIPPED
  use simple_gui_metadata_api
  implicit none

public :: gui_assembler
private
#include "simple_local_flags.inc"

type :: gui_assembler
  private
  type(json_core)           :: json
  type(json_value), pointer :: json_root              ! root of the assembled JSON document
  type(string)              :: preprocess_hash         ! FNV-1a hash of last sent preprocessing section
  type(string)              :: optics_assignment_hash  ! FNV-1a hash of last sent optics-assignment section
  type(string)              :: initial_picking_hash    ! FNV-1a hash of last sent initial-picking section
  type(string)              :: reference_picking_hash  ! FNV-1a hash of last sent reference-picking section
  type(string)              :: opening2D_hash          ! FNV-1a hash of last sent opening2D section
  integer                   :: job_id    = 0           ! pipeline job identifier
  integer                   :: starttime = 0           ! Unix timestamp of job start
  integer                   :: stoptime  = 0           ! Unix timestamp of job stop (0 while running)
  logical                   :: init      = .false.     ! .true. once new() has been called
  real                      :: version   = 1.0         ! JSON schema version sent to the frontend
contains
  procedure :: new
  procedure :: kill
  procedure :: to_string
  procedure :: set_stoptime
  procedure :: clear_hashes
  procedure :: is_associated
  procedure :: assemble_stream_heartbeat
  procedure :: assemble_stream_preprocess
  procedure :: assemble_stream_optics_assignment
  procedure :: assemble_stream_initial_picking
  procedure :: assemble_stream_reference_picking
  procedure :: assemble_stream_opening2D
end type gui_assembler

contains

  ! Initialise the JSON tree for a new job, recording the start timestamp.
  ! Kills any previously initialised state before reinitialising.
  subroutine new( self, jobid )
    class(gui_assembler), intent(inout) :: self
    integer,              intent(in)    :: jobid
    if( self%init ) call self%kill()
    self%init      = .true.
    self%job_id    = jobid
    self%starttime = int(c_time(0_c_long))
    call self%json%initialize(no_whitespace=.true., compact_reals=.true.)
    call self%json%create_object(self%json_root, '') 
    call self%json%add(self%json_root, 'jobid',   self%job_id       )
    call self%json%add(self%json_root, 'version', dble(self%version))
  end subroutine new

  ! Destroy the JSON tree and reset all state to defaults.
  subroutine kill( self )
    class(gui_assembler), intent(inout) :: self
    self%starttime = 0
    self%stoptime  = 0
    self%job_id    = 0
    self%init      = .false.
    call self%json%destroy(self%json_root)
    nullify(self%json_root)
  end subroutine kill

  ! Reset change-detection hashes so every section is retransmitted on next assembly.
  subroutine clear_hashes( self )
    class(gui_assembler), intent(inout) :: self
    call self%preprocess_hash%kill()
    call self%optics_assignment_hash%kill()
    call self%initial_picking_hash%kill()
    call self%reference_picking_hash%kill()
    call self%opening2D_hash%kill()
  end subroutine clear_hashes

  ! Write the stream_heartbeat section: per-process status fields plus a master
  ! aggregate status derived from the union of all child-process states.
  subroutine assemble_stream_heartbeat( self, fork_preprocess, fork_assign_optics, fork_opening2D, fork_reference_picking )
    class(gui_assembler),  intent(inout) :: self
    class(forked_process), intent(inout) :: fork_preprocess, fork_assign_optics, fork_opening2D
    class(forked_process), intent(inout) :: fork_reference_picking
    type(json_value),      pointer       :: json_ptr, json_master_ptr
    integer                              :: n_running, n_failed, n_restarting, n_unknown
    n_running    = 0
    n_failed     = 0
    n_restarting = 0
    n_unknown    = 0
    call self%json%remove_if_present(self%json_root, 'stream_heartbeat')
    call self%json%create_object(json_ptr, 'stream_heartbeat')
    call forked_process_status(string('preprocessing'),    fork_preprocess)
    call forked_process_status(string('assign_optics'), fork_assign_optics)
    call forked_process_status(string('initial_picking'),   fork_opening2D)
    call forked_process_status(string('opening2D'),         fork_opening2D)
    call forked_process_status(string('reference_picking'), fork_reference_picking)
    ! global status
    call self%json%create_object(json_master_ptr, 'master')
    call self%json%add(json_master_ptr, 'timestamp', int(c_time(0_c_long)))
    call self%json%add(json_master_ptr, 'starttime',        self%starttime)
    call self%json%add(json_master_ptr, 'stoptime',          self%stoptime)
    call self%json%add(json_master_ptr, 'pid',                    getpid())
    if( n_unknown > 0 ) then
      call self%json%add(json_master_ptr, 'status', 'unknown')
    else if( n_failed > 0 ) then
      call self%json%add(json_master_ptr, 'status', 'failed')
    else if( n_restarting > 0 ) then
      call self%json%add(json_master_ptr, 'status', 'running')
    else if( n_running == 0 ) then
      call self%json%add(json_master_ptr, 'status', 'finished')
    else
      call self%json%add(json_master_ptr, 'status', 'running')
    endif
    call self%json%add(json_ptr, json_master_ptr)
    call self%json%add(self%json_root, json_ptr)
    nullify(json_ptr)

  contains

    subroutine forked_process_status(my_process_name, my_forked_process)
      class(forked_process), intent(inout) :: my_forked_process
      type(string),          intent(in)    :: my_process_name
      type(json_value),      pointer       :: my_json_proc_ptr
      call self%json%create_object(my_json_proc_ptr, my_process_name%to_char())
      call self%json%add(my_json_proc_ptr, 'timestamp',             int(c_time(0_c_long)))
      call self%json%add(my_json_proc_ptr, 'pid',             my_forked_process%get_pid())
      call self%json%add(my_json_proc_ptr, 'restarts',  my_forked_process%get_nrestarts())
      call self%json%add(my_json_proc_ptr, 'queuetime', my_forked_process%get_queuetime())
      call self%json%add(my_json_proc_ptr, 'starttime', my_forked_process%get_starttime())
      call self%json%add(my_json_proc_ptr, 'failtime',   my_forked_process%get_failtime())
      call self%json%add(my_json_proc_ptr, 'stoptime',   my_forked_process%get_stoptime())
      select case(my_forked_process%status())
        case(FORK_STATUS_RUNNING)
          call self%json%add(my_json_proc_ptr, 'status', 'running')
          n_running = n_running + 1
        case(FORK_STATUS_FAILED)
          call self%json%add(my_json_proc_ptr, 'status', 'failed')
          n_failed = n_failed + 1
        case(FORK_STATUS_STOPPED)
          call self%json%add(my_json_proc_ptr, 'status', 'finished')
        case(FORK_STATUS_RESTARTING)
          call self%json%add(my_json_proc_ptr, 'status', 'restarting')
          n_restarting = n_restarting + 1
        case(FORK_STATUS_SKIPPED)
          call self%json%add(my_json_proc_ptr, 'status', 'skipped')
        case DEFAULT
          call self%json%add(my_json_proc_ptr, 'status', 'unknown')
          n_unknown = n_unknown + 1
      end select
      call self%json%add(json_ptr, my_json_proc_ptr)
    end subroutine forked_process_status

  end subroutine assemble_stream_heartbeat

  ! Write the preprocessing section, including micrographs, histograms, and
  ! timeplots. The whole section is suppressed when its hash matches the
  ! previously sent hash.
  subroutine assemble_stream_preprocess( self, meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_preprocess),              intent(inout) :: meta_preprocess
    type(gui_metadata_micrograph),        allocatable, intent(inout) :: meta_micrographs(:)
    type(gui_metadata_histogram),         allocatable, intent(inout) :: meta_histograms(:)
    type(gui_metadata_timeplot),          allocatable, intent(inout) :: meta_timeplots(:)
    character(kind=CK,len=:),             allocatable                :: buffer
    type(json_value),                     pointer                    :: json_ptr, json_mics_ptr, json_hists_ptr, json_timeplots_ptr
    type(string)                                                     :: str, hash
    integer                                                          :: i_mic, i_hist, i_timeplot
    logical                                                          :: l_add
    call self%json%remove_if_present(self%json_root, 'preprocessing')
    json_ptr => meta_preprocess%jsonise()
    if( .not. associated(json_ptr) ) return
    call self%json%rename(json_ptr, 'preprocessing')
    if( allocated(meta_micrographs) ) then
      l_add = .false.
      call self%json%create_array(json_mics_ptr, 'latest_micrographs')
      do i_mic=1, size(meta_micrographs)
        if( meta_micrographs(i_mic)%assigned() ) then
          l_add = .true.
          call self%json%add(json_mics_ptr, meta_micrographs(i_mic)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_mics_ptr)
    endif
    if( allocated(meta_histograms) ) then
      l_add = .false.
      call self%json%create_object(json_hists_ptr, 'histograms')
      do i_hist=1, size(meta_histograms)
        if( meta_histograms(i_hist)%assigned() ) then
          l_add = .true.
          call self%json%add(json_hists_ptr, meta_histograms(i_hist)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_hists_ptr)
    endif
    if( allocated(meta_timeplots) ) then
      l_add = .false.
      call self%json%create_object(json_timeplots_ptr, 'timeplots')
      do i_timeplot=1, size(meta_timeplots)
        if( meta_timeplots(i_timeplot)%assigned() ) then
          l_add = .true.
          call self%json%add(json_timeplots_ptr, meta_timeplots(i_timeplot)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_timeplots_ptr)
    endif
    call self%json%print_to_string(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%preprocess_hash ) then
      call self%json%add(self%json_root, json_ptr)
      self%preprocess_hash = hash
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr, json_mics_ptr, json_hists_ptr, json_timeplots_ptr)
  end subroutine assemble_stream_preprocess

  ! Write the optics-assignment section, including optics groups. The whole
  ! section is suppressed when its hash matches the previously sent hash.
  subroutine assemble_stream_optics_assignment( self, meta_optics_assignment, meta_optics_groups )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_optics_assignment),       intent(inout) :: meta_optics_assignment
    type(gui_metadata_optics_group),      allocatable, intent(inout) :: meta_optics_groups(:)
    character(kind=CK,len=:),             allocatable                :: buffer
    type(json_value),                     pointer                    :: json_ptr, json_optics_groups_ptr
    type(string)                                                     :: str, hash
    logical                                                          :: l_add
    integer                                                          :: i_group
    call self%json%remove_if_present(self%json_root, 'optics_assignment')
    json_ptr => meta_optics_assignment%jsonise()
    if( .not. associated(json_ptr) ) return
    call self%json%rename(json_ptr, 'optics_assignment')
    if( allocated(meta_optics_groups) ) then
      l_add = .false.
      call self%json%create_array(json_optics_groups_ptr, 'optics_assignments')
      do i_group=1, size(meta_optics_groups)
        if( meta_optics_groups(i_group)%assigned() ) then
          l_add = .true.
          call self%json%add(json_optics_groups_ptr, meta_optics_groups(i_group)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_optics_groups_ptr)
    endif
    call self%json%print_to_string(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%optics_assignment_hash ) then
      call self%json%add(self%json_root, json_ptr)
      self%optics_assignment_hash = hash
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr)
  end subroutine assemble_stream_optics_assignment

  ! Write the initial-picking section, including micrographs. The whole section
  ! is suppressed when its hash matches the previously sent hash.
  subroutine assemble_stream_initial_picking( self, meta_initial_picking, meta_micrographs )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_picking),         intent(inout) :: meta_initial_picking
    type(gui_metadata_micrograph),        allocatable, intent(inout) :: meta_micrographs(:)
    character(kind=CK,len=:),             allocatable                :: buffer
    type(json_value),                     pointer                    :: json_ptr, json_mics_ptr
    type(string)                                                     :: str, hash
    logical                                                          :: l_add
    integer                                                          :: i_mic
    call self%json%remove_if_present(self%json_root, 'initial_picking')
    json_ptr => meta_initial_picking%jsonise()
    if( .not. associated(json_ptr) ) return
    call self%json%rename(json_ptr, 'initial_picking')
    if( allocated(meta_micrographs) ) then
      l_add = .false.
      call self%json%create_array(json_mics_ptr, 'latest_micrographs')
      do i_mic=1, size(meta_micrographs)
        if( meta_micrographs(i_mic)%assigned() ) then
          l_add = .true.
          call self%json%add(json_mics_ptr, meta_micrographs(i_mic)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_mics_ptr)
    endif
    call self%json%print_to_string(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%initial_picking_hash ) then
      call self%json%add(self%json_root, json_ptr)
      self%initial_picking_hash = hash
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr)
  end subroutine assemble_stream_initial_picking

  ! Write the reference-picking section, including micrographs. The whole section
  ! is suppressed when its hash matches the previously sent hash.
  subroutine assemble_stream_reference_picking( self, meta_reference_picking, meta_micrographs, meta_cavgs2D )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_picking),                 intent(inout) :: meta_reference_picking
    type(gui_metadata_micrograph),        allocatable, intent(inout) :: meta_micrographs(:)
    type(gui_metadata_cavg2D),            allocatable, intent(inout) :: meta_cavgs2D(:)
    character(kind=CK,len=:),             allocatable                :: buffer
    type(json_value),                     pointer                    :: json_ptr, json_mics_ptr
    type(string)                                                     :: str, hash
    logical                                                          :: l_add
    integer                                                          :: i_mic
    call self%json%remove_if_present(self%json_root, 'reference_picking')
    json_ptr => meta_reference_picking%jsonise()
    if( .not. associated(json_ptr) ) return
    call self%json%rename(json_ptr, 'reference_picking')
    if( allocated(meta_micrographs) ) then
      l_add = .false.
      call self%json%create_array(json_mics_ptr, 'latest_micrographs')
      do i_mic=1, size(meta_micrographs)
        if( meta_micrographs(i_mic)%assigned() ) then
          l_add = .true.
          call self%json%add(json_mics_ptr, meta_micrographs(i_mic)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_mics_ptr)
    endif
    if( allocated(meta_cavgs2D) ) then
      l_add = .false.
      call self%json%create_array(json_mics_ptr, 'picking_references')
      do i_mic=1, size(meta_cavgs2D)
        if( meta_cavgs2D(i_mic)%assigned() ) then
          l_add = .true.
          call self%json%add(json_mics_ptr, meta_cavgs2D(i_mic)%jsonise())
        endif
      enddo
      if( l_add ) call self%json%add(json_ptr, json_mics_ptr)
    endif
    call self%json%print_to_string(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%reference_picking_hash ) then
      call self%json%add(self%json_root, json_ptr)
      self%reference_picking_hash = hash
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr)
  end subroutine assemble_stream_reference_picking

  ! Write the opening2D (2D classification) section, including any cavgs2D.
  ! The whole section (header + cavgs2D) is suppressed when its hash matches
  ! the previously sent hash.
  subroutine assemble_stream_opening2D( self, meta_opening2D, meta_latest_cavgs2D, meta_final_cavgs2D )
    class(gui_assembler),                   intent(inout) :: self
    type(gui_metadata_cavg2D), allocatable, intent(inout) :: meta_latest_cavgs2D(:), meta_final_cavgs2D(:)
    type(gui_metadata_stream_opening2D),    intent(inout) :: meta_opening2D
    character(kind=CK,len=:),               allocatable   :: buffer
    type(json_value),                       pointer       :: json_ptr, json_cavgs2D_ptr
    type(string)                                          :: str, hash
    logical                                               :: l_add
    integer                                               :: i_cls2D
    call self%json%remove_if_present(self%json_root, 'opening2D')
    json_ptr => meta_opening2D%jsonise()
    if( .not. associated(json_ptr) ) return
    call self%json%rename(json_ptr, 'opening2D')
    l_add = .false.
    if( allocated(meta_final_cavgs2D) ) then
      call self%json%create_array(json_cavgs2D_ptr, 'final_cls2D')
      do i_cls2D=1, size(meta_final_cavgs2D)
        if( meta_final_cavgs2D(i_cls2D)%assigned() ) then
          l_add = .true.
          call self%json%add(json_cavgs2D_ptr, meta_final_cavgs2D(i_cls2D)%jsonise())
        endif
      enddo
    else if( allocated(meta_latest_cavgs2D) ) then
      call self%json%create_array(json_cavgs2D_ptr, 'latest_cls2D')
      do i_cls2D=1, size(meta_latest_cavgs2D)
        if( meta_latest_cavgs2D(i_cls2D)%assigned() ) then
          l_add = .true.
          call self%json%add(json_cavgs2D_ptr, meta_latest_cavgs2D(i_cls2D)%jsonise())
        endif
      enddo
    endif
    if( l_add ) call self%json%add(json_ptr, json_cavgs2D_ptr)
    call self%json%print_to_string(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%opening2D_hash ) then
      call self%json%add(self%json_root, json_ptr)
      self%opening2D_hash = hash
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr)
  end subroutine assemble_stream_opening2D

  ! Record the job stop timestamp (call when the pipeline finishes).
  subroutine set_stoptime( self )
    class(gui_assembler), intent(inout) :: self
    self%stoptime = int(c_time(0_c_long))
  end subroutine set_stoptime

  ! Serialise the current JSON tree to a compact string.
  function to_string( self ) result( str )
    class(gui_assembler),     intent(inout) :: self
    type(string)                            :: str
    character(kind=CK,len=:), allocatable   :: buffer
    call self%json%print_to_string(self%json_root, buffer)
    str = buffer
    if( allocated(buffer) ) deallocate(buffer)
  end function to_string

  ! Return .true. if the JSON root pointer is associated (i.e. new() has been called).
  function is_associated( self ) result( assoc )
    class(gui_assembler), intent(in) :: self
    logical                          :: assoc
    assoc = associated(self%json_root)
  end function is_associated

end module simple_gui_assembler