!@descr: builds the GUI JSON document (stream and batch) from gui_metadata objects
! Each assemble_* replaces its section of json_root. Content sections are dropped when their FNV-1a
! hash matches the last one sent (heartbeats always go). Call clear_hashes() after a failed send.
! The stream heartbeat takes plain stage-status records (gui_stage_status), which the master fills
! from its forked stages: the assembler imports nothing that reaches src/main.
module simple_gui_assembler
  use unix,                                          only: c_time, c_long
  use json_kinds,                                    only: CK
  use json_module,                                   only: json_core, json_value
  use simple_string,                                 only: string
  use simple_gui_metadata_base,                      only: gui_metadata_base
  use simple_gui_metadata_micrograph,                only: gui_metadata_micrograph
  use simple_gui_metadata_histogram,                 only: gui_metadata_histogram
  use simple_gui_metadata_timeplot,                  only: gui_metadata_timeplot
  use simple_gui_metadata_optics_group,              only: gui_metadata_optics_group
  use simple_gui_metadata_cavg2D,                    only: gui_metadata_cavg2D
  use simple_gui_metadata_vol3D,                     only: gui_metadata_vol3D
  use simple_gui_metadata_stream_preprocess,         only: gui_metadata_stream_preprocess
  use simple_gui_metadata_stream_optics_assignment,  only: gui_metadata_stream_optics_assignment
  use simple_gui_metadata_stream_picking,            only: gui_metadata_stream_picking
  use simple_gui_metadata_stream_initial_analysis,   only: gui_metadata_stream_initial_analysis
  use simple_gui_metadata_stream_particle_sieving,   only: gui_metadata_stream_particle_sieving
  use simple_gui_metadata_stream_pool2D,             only: gui_metadata_stream_pool2D
  use simple_gui_metadata_stream_snapshot,           only: gui_metadata_stream_snapshot
  use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
  implicit none

public :: gui_assembler, gui_stage_status
public :: GUI_STAGE_STATUS_UNKNOWN, GUI_STAGE_STATUS_RUNNING, GUI_STAGE_STATUS_FAILED
public :: GUI_STAGE_STATUS_FINISHED, GUI_STAGE_STATUS_SKIPPED
private
#include "simple_local_flags.inc"

! a stage's status in the stream heartbeat
integer, parameter :: GUI_STAGE_STATUS_UNKNOWN  = 0
integer, parameter :: GUI_STAGE_STATUS_RUNNING  = 1
integer, parameter :: GUI_STAGE_STATUS_FAILED   = 2
integer, parameter :: GUI_STAGE_STATUS_FINISHED = 3
integer, parameter :: GUI_STAGE_STATUS_SKIPPED  = 4

! a stage's process in the stream heartbeat: its pid, its times (Unix) and its status
type :: gui_stage_status
  integer :: pid       = 0
  integer :: queuetime = 0
  integer :: starttime = 0
  integer :: failtime  = 0
  integer :: stoptime  = 0
  integer :: status    = GUI_STAGE_STATUS_UNKNOWN
end type gui_stage_status

! the content sections, each with the hash of its last text sent
integer, parameter :: SECTION_PREPROCESS         = 1
integer, parameter :: SECTION_OPTICS_ASSIGNMENT  = 2
integer, parameter :: SECTION_INITIAL_PICKING    = 3
integer, parameter :: SECTION_REFERENCE_PICKING  = 4
integer, parameter :: SECTION_INITIAL_ANALYSIS   = 5
integer, parameter :: SECTION_PARTICLE_SIEVING   = 6
integer, parameter :: SECTION_POOL2D             = 7
integer, parameter :: SECTION_SOLVE3D_MULTISTATE = 8
integer, parameter :: SECTION_PROJECT            = 9
integer, parameter :: NSECTIONS                  = 9

type :: gui_assembler
  private
  type(json_core)           :: json
  type(json_value), pointer :: json_root => null()     ! root of the assembled JSON document
  type(string)              :: hashes(NSECTIONS)       ! FNV-1a hash of each content section last sent
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
  procedure :: assemble_batch_heartbeat
  procedure :: assemble_batch_metadata
  procedure :: assemble_stream_preprocess
  procedure :: assemble_stream_optics_assignment
  procedure :: assemble_stream_initial_picking
  procedure :: assemble_stream_reference_picking
  procedure :: assemble_stream_initial_analysis
  procedure :: assemble_stream_particle_sieving
  procedure :: assemble_stream_pool2D
  procedure :: assemble_stream_solve3D_multistate
  procedure, private :: open_section
  procedure, private :: add_assigned_micrographs
  procedure, private :: add_assigned_histograms
  procedure, private :: add_assigned_timeplots
  procedure, private :: add_assigned_optics_groups
  procedure, private :: add_assigned_cavgs2D
  generic,   private :: add_assigned => add_assigned_micrographs, add_assigned_histograms, add_assigned_timeplots,&
                                        add_assigned_optics_groups, add_assigned_cavgs2D
  procedure, private :: open_list
  procedure, private :: add_item
  procedure, private :: close_list
  procedure, private :: commit_section
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
    call self%clear_hashes()
    call self%json%destroy(self%json_root)
    nullify(self%json_root)
  end subroutine kill

  ! Reset change-detection hashes so every section is retransmitted on next assembly.
  subroutine clear_hashes( self )
    class(gui_assembler), intent(inout) :: self
    integer :: i
    do i = 1,NSECTIONS
      call self%hashes(i)%kill()
    enddo
  end subroutine clear_hashes

  ! Write the stream_heartbeat section: per-process status fields plus a master
  ! aggregate status derived from the union of all child-process states.
  subroutine assemble_stream_heartbeat( self, preprocess, assign_optics, initial_analysis, reference_picking,&
      particle_sieving, pool2D, solve3D_multistate, n_active_persistent_workers )
    class(gui_assembler),   intent(inout) :: self
    type(gui_stage_status), intent(in)    :: preprocess, assign_optics, initial_analysis, reference_picking
    type(gui_stage_status), intent(in)    :: particle_sieving, pool2D, solve3D_multistate
    integer, optional,      intent(in)    :: n_active_persistent_workers
    type(json_value),       pointer       :: json_ptr, json_master_ptr
    integer                               :: n_running, n_failed, n_unknown
    n_running    = 0
    n_failed     = 0
    n_unknown    = 0
    call self%json%remove_if_present(self%json_root, 'stream_heartbeat')
    call self%json%create_object(json_ptr, 'stream_heartbeat')
    call add_stage_status('preprocessing',      preprocess)
    call add_stage_status('assign_optics',      assign_optics)
    ! the initial analysis reports as initial picking and under its NICE key, opening2D
    call add_stage_status('initial_picking',    initial_analysis)
    call add_stage_status('opening2D',          initial_analysis)
    call add_stage_status('reference_picking',  reference_picking)
    call add_stage_status('particle_sieving',   particle_sieving)
    call add_stage_status('pool2D',             pool2D)
    call add_stage_status('solve3D_multistate', solve3D_multistate)
    ! global status
    call self%json%create_object(json_master_ptr, 'master')
    call self%json%add(json_master_ptr, 'timestamp', int(c_time(0_c_long)))
    call self%json%add(json_master_ptr, 'starttime',        self%starttime)
    call self%json%add(json_master_ptr, 'stoptime',          self%stoptime)
    call self%json%add(json_master_ptr, 'pid',                    getpid())
    if( present(n_active_persistent_workers) ) &
      call self%json%add(json_master_ptr, 'active_persistent_workers', n_active_persistent_workers)
    if( n_unknown > 0 ) then
      call self%json%add(json_master_ptr, 'status', 'unknown')
    else if( n_failed > 0 ) then
      if( n_running == 0 ) then
        call self%json%add(json_master_ptr, 'status', 'failed')
      else
        call self%json%add(json_master_ptr, 'status', 'error')
      endif
    else if( n_running == 0 ) then
      call self%json%add(json_master_ptr, 'status', 'finished')
    else
      call self%json%add(json_master_ptr, 'status', 'running')
    endif
    call self%json%add(json_ptr, json_master_ptr)
    call self%json%add(self%json_root, json_ptr)
    nullify(json_ptr)

  contains

    subroutine add_stage_status( name, stage )
      character(len=*),       intent(in) :: name
      type(gui_stage_status), intent(in) :: stage
      type(json_value),       pointer    :: json_stage_ptr
      call self%json%create_object(json_stage_ptr, name)
      call self%json%add(json_stage_ptr, 'timestamp', int(c_time(0_c_long)))
      call self%json%add(json_stage_ptr, 'pid',       stage%pid)
      call self%json%add(json_stage_ptr, 'queuetime', stage%queuetime)
      call self%json%add(json_stage_ptr, 'starttime', stage%starttime)
      call self%json%add(json_stage_ptr, 'failtime',  stage%failtime)
      call self%json%add(json_stage_ptr, 'stoptime',  stage%stoptime)
      select case( stage%status )
        case( GUI_STAGE_STATUS_RUNNING )
          call self%json%add(json_stage_ptr, 'status', 'running')
          n_running = n_running + 1
        case( GUI_STAGE_STATUS_FAILED )
          call self%json%add(json_stage_ptr, 'status', 'failed')
          n_failed = n_failed + 1
        case( GUI_STAGE_STATUS_FINISHED )
          call self%json%add(json_stage_ptr, 'status', 'finished')
        case( GUI_STAGE_STATUS_SKIPPED )
          call self%json%add(json_stage_ptr, 'status', 'skipped')
        case DEFAULT
          call self%json%add(json_stage_ptr, 'status', 'unknown')
          n_unknown = n_unknown + 1
      end select
      call self%json%add(json_ptr, json_stage_ptr)
    end subroutine add_stage_status

  end subroutine assemble_stream_heartbeat

  ! Write the batch_heartbeat section: single top-level status for batch (non-streaming) jobs.
  subroutine assemble_batch_heartbeat( self )
    class(gui_assembler), intent(inout) :: self
    type(json_value),     pointer       :: json_ptr
    call self%json%remove_if_present(self%json_root, 'batch_heartbeat')
    call self%json%create_object(json_ptr, 'batch_heartbeat')
    call self%json%add(json_ptr, 'timestamp', int(c_time(0_c_long)))
    call self%json%add(json_ptr, 'starttime', self%starttime)
    call self%json%add(json_ptr, 'stoptime',  self%stoptime)
    call self%json%add(json_ptr, 'pid',       getpid())
    if( self%stoptime > 0 ) then
      call self%json%add(json_ptr, 'status', 'finished')
    else
      call self%json%add(json_ptr, 'status', 'running')
    endif
    call self%json%add(self%json_root, json_ptr)
    nullify(json_ptr)
  end subroutine assemble_batch_heartbeat

  ! Write the project section from the batch (non-streaming) project metadata, a
  ! gui_metadata_project: only its JSON is needed, so the assembler does not import it.
  subroutine assemble_batch_metadata( self, meta_project )
    class(gui_assembler),     intent(inout) :: self
    class(gui_metadata_base), intent(in)    :: meta_project
    type(json_value),         pointer       :: json_ptr
    if( .not. self%open_section(json_ptr, 'project_metadata', meta_project) ) return
    call self%commit_section(json_ptr, SECTION_PROJECT)
  end subroutine assemble_batch_metadata

  ! Write the preprocessing section, including micrographs, histograms, and timeplots.
  subroutine assemble_stream_preprocess( self, meta_preprocess, meta_micrographs, meta_histograms, meta_timeplots )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_preprocess),              intent(in)    :: meta_preprocess
    type(gui_metadata_micrograph),        allocatable, intent(in)    :: meta_micrographs(:)
    type(gui_metadata_histogram),         allocatable, intent(in)    :: meta_histograms(:)
    type(gui_metadata_timeplot),          allocatable, intent(in)    :: meta_timeplots(:)
    type(json_value),                     pointer                    :: json_ptr
    if( .not. self%open_section(json_ptr, 'preprocessing', meta_preprocess) ) return
    if( allocated(meta_micrographs) ) call self%add_assigned(json_ptr, 'latest_micrographs', meta_micrographs, l_object=.false.)
    if( allocated(meta_histograms)  ) call self%add_assigned(json_ptr, 'histograms',         meta_histograms,  l_object=.true.)
    if( allocated(meta_timeplots)   ) call self%add_assigned(json_ptr, 'timeplots',          meta_timeplots,   l_object=.true.)
    call self%commit_section(json_ptr, SECTION_PREPROCESS)
  end subroutine assemble_stream_preprocess

  ! Write the optics-assignment section, including optics groups.
  subroutine assemble_stream_optics_assignment( self, meta_optics_assignment, meta_optics_groups )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_optics_assignment),       intent(in)    :: meta_optics_assignment
    type(gui_metadata_optics_group),      allocatable, intent(in)    :: meta_optics_groups(:)
    type(json_value),                     pointer                    :: json_ptr
    if( .not. self%open_section(json_ptr, 'optics_assignment', meta_optics_assignment) ) return
    if( allocated(meta_optics_groups) ) call self%add_assigned(json_ptr, 'optics_assignments', meta_optics_groups, l_object=.false.)
    call self%commit_section(json_ptr, SECTION_OPTICS_ASSIGNMENT)
  end subroutine assemble_stream_optics_assignment

  ! Write the initial-picking section, including micrographs.
  subroutine assemble_stream_initial_picking( self, meta_initial_picking, meta_micrographs )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_picking),                 intent(in)    :: meta_initial_picking
    type(gui_metadata_micrograph),        allocatable, intent(in)    :: meta_micrographs(:)
    type(json_value),                     pointer                    :: json_ptr
    if( .not. self%open_section(json_ptr, 'initial_picking', meta_initial_picking) ) return
    if( allocated(meta_micrographs) ) call self%add_assigned(json_ptr, 'latest_micrographs', meta_micrographs, l_object=.false.)
    call self%commit_section(json_ptr, SECTION_INITIAL_PICKING)
  end subroutine assemble_stream_initial_picking

  ! Write the reference-picking section, including micrographs and the picking references.
  subroutine assemble_stream_reference_picking( self, meta_reference_picking, meta_micrographs, meta_cavgs2D )
    class(gui_assembler),                              intent(inout) :: self
    type(gui_metadata_stream_picking),                 intent(in)    :: meta_reference_picking
    type(gui_metadata_micrograph),        allocatable, intent(in)    :: meta_micrographs(:)
    type(gui_metadata_cavg2D),            allocatable, intent(in)    :: meta_cavgs2D(:)
    type(json_value),                     pointer                    :: json_ptr
    if( .not. self%open_section(json_ptr, 'reference_picking', meta_reference_picking) ) return
    if( allocated(meta_micrographs) ) call self%add_assigned(json_ptr, 'latest_micrographs', meta_micrographs, l_object=.false.)
    if( allocated(meta_cavgs2D)     ) call self%add_assigned(json_ptr, 'picking_references', meta_cavgs2D,     l_object=.false.)
    call self%commit_section(json_ptr, SECTION_REFERENCE_PICKING)
  end subroutine assemble_stream_reference_picking

  ! Write the initial-analysis section under its NICE key, opening2D, including the latest and
  ! the selected class averages. meta_vol3D, when present and assigned, holds the metadata for the
  ! single solve3D_cavgs volume chosen for reprojection/picking references and is embedded as a
  ! nested 'volume' object (no FSC/postprocessed fields at this stage).
  subroutine assemble_stream_initial_analysis( self, meta_initial_analysis, meta_latest_cavgs2D, meta_selected_pickrefs, meta_vol3D )
    class(gui_assembler),                        intent(inout) :: self
    type(gui_metadata_stream_initial_analysis),  intent(in)    :: meta_initial_analysis
    type(gui_metadata_cavg2D),      allocatable, intent(in)    :: meta_latest_cavgs2D(:), meta_selected_pickrefs(:)
    type(gui_metadata_vol3D),       optional,    intent(in)    :: meta_vol3D
    type(json_value),               pointer                    :: json_ptr, json_vol3D_ptr
    if( .not. self%open_section(json_ptr, 'opening2D', meta_initial_analysis) ) return
    if( present(meta_vol3D) ) then
      if( meta_vol3D%assigned() ) then
        json_vol3D_ptr => meta_vol3D%jsonise()
        if( associated(json_vol3D_ptr) ) then
          call self%json%rename(json_vol3D_ptr, 'volume')
          call self%json%add(json_ptr, json_vol3D_ptr)
        endif
      endif
    endif
    if( allocated(meta_selected_pickrefs) ) call self%add_assigned(json_ptr, 'selected_pickrefs', meta_selected_pickrefs, l_object=.false.)
    if( allocated(meta_latest_cavgs2D)    ) call self%add_assigned(json_ptr, 'latest_cls2D',      meta_latest_cavgs2D,    l_object=.false.)
    call self%commit_section(json_ptr, SECTION_INITIAL_ANALYSIS)
  end subroutine assemble_stream_initial_analysis

  ! Write the particle-sieving section, including the reference class averages and the
  ! latest class averages.
  subroutine assemble_stream_particle_sieving( self, meta_particle_sieving, meta_latest_cavgs2D, meta_ref_cavgs )
    class(gui_assembler),                        intent(inout) :: self
    type(gui_metadata_stream_particle_sieving),  intent(in)    :: meta_particle_sieving
    type(gui_metadata_cavg2D),      allocatable, intent(in)    :: meta_latest_cavgs2D(:), meta_ref_cavgs(:)
    type(json_value),               pointer                    :: json_ptr
    if( .not. self%open_section(json_ptr, 'particle_sieving', meta_particle_sieving) ) return
    if( allocated(meta_ref_cavgs)      ) call self%add_assigned(json_ptr, 'ref_cls2D',    meta_ref_cavgs,      l_object=.false.)
    if( allocated(meta_latest_cavgs2D) ) call self%add_assigned(json_ptr, 'latest_cls2D', meta_latest_cavgs2D, l_object=.false.)
    call self%commit_section(json_ptr, SECTION_PARTICLE_SIEVING)
  end subroutine assemble_stream_particle_sieving

  ! Write the pool-2D section, including the latest class averages.
  ! meta_snapshot_cavgs2D, when present, holds the cavg2D metadata for the class
  ! averages selected for the current snapshot and is embedded as 'cls2D' nested
  ! under the 'snapshot' object.
  subroutine assemble_stream_pool2D( self, meta_pool2D, meta_latest_cavgs2D, meta_pool2D_snapshot, meta_snapshot_cavgs2D )
    class(gui_assembler),                                intent(inout) :: self
    type(gui_metadata_stream_pool2D),                    intent(in)    :: meta_pool2D
    type(gui_metadata_cavg2D),   allocatable,            intent(in)    :: meta_latest_cavgs2D(:)
    type(gui_metadata_stream_snapshot),           intent(in)    :: meta_pool2D_snapshot
    type(gui_metadata_cavg2D),   allocatable, optional,  intent(in)    :: meta_snapshot_cavgs2D(:)
    type(json_value),            pointer                               :: json_ptr, json_snapshot_ptr
    if( .not. self%open_section(json_ptr, 'pool2D', meta_pool2D) ) return
    if( meta_pool2D_snapshot%assigned() ) then
      json_snapshot_ptr => meta_pool2D_snapshot%jsonise()
      if( associated(json_snapshot_ptr) ) then
        call self%json%rename(json_snapshot_ptr, 'snapshot')
        if( present(meta_snapshot_cavgs2D) ) then
          if( allocated(meta_snapshot_cavgs2D) ) call self%add_assigned(json_snapshot_ptr, 'cls2D', meta_snapshot_cavgs2D, l_object=.false.)
        endif
        call self%json%add(json_ptr, json_snapshot_ptr)
      end if
    end if
    if( allocated(meta_latest_cavgs2D) ) call self%add_assigned(json_ptr, 'latest_cls2D', meta_latest_cavgs2D, l_object=.false.)
    call self%commit_section(json_ptr, SECTION_POOL2D)
  end subroutine assemble_stream_pool2D

  ! Write the multistate solve3D section, with the latest 3D snapshot report, when assigned,
  ! nested as 'snapshot' (as pool2D nests its own).
  ! meta_states_vol3D, when present, holds the vol3D metadata for the current
  ! per-state reconstructed volumes and is embedded as a 'state_volumes' array,
  ! distinct from the lightweight per-state 'states' array already emitted by
  ! meta_solve3D_multistate%jsonise() to avoid a duplicate JSON key.
  ! meta_reprojtiles, when present, holds the individual orthogonal reprojection
  ! tiles (gui_metadata_cavg2D, idx=state) and is nested per-state as a
  ! 'reprojtiles' array inside the matching state_volumes entry; a state without
  ! tiles has none.
  subroutine assemble_stream_solve3D_multistate( self, meta_solve3D_multistate, meta_snapshot, meta_states_vol3D, meta_reprojtiles )
    class(gui_assembler),                                 intent(inout) :: self
    type(gui_metadata_stream_solve3D_multistate),         intent(in)    :: meta_solve3D_multistate
    type(gui_metadata_stream_snapshot),                   intent(in)    :: meta_snapshot
    type(gui_metadata_vol3D),  allocatable, optional,     intent(in)    :: meta_states_vol3D(:)
    type(gui_metadata_cavg2D), allocatable, optional,     intent(in)    :: meta_reprojtiles(:)
    type(json_value),          pointer                                  :: json_ptr, json_states_ptr, json_state_vol_ptr
    type(json_value),          pointer                                  :: json_snapshot_ptr
    logical                                                             :: l_add
    integer                                                             :: i_state
    if( .not. self%open_section(json_ptr, 'solve3D_multistate', meta_solve3D_multistate) ) return
    if( meta_snapshot%assigned() ) then
      json_snapshot_ptr => meta_snapshot%jsonise()
      if( associated(json_snapshot_ptr) ) then
        call self%json%rename(json_snapshot_ptr, 'snapshot')
        call self%json%add(json_ptr, json_snapshot_ptr)
      end if
    end if
    if( present(meta_states_vol3D) ) then
      if( allocated(meta_states_vol3D) ) then
        l_add = .false.
        call self%json%create_array(json_states_ptr, 'state_volumes')
        do i_state = 1, size(meta_states_vol3D)
          if( .not. meta_states_vol3D(i_state)%assigned() ) cycle
          l_add = .true.
          json_state_vol_ptr => meta_states_vol3D(i_state)%jsonise()
          if( present(meta_reprojtiles) ) then
            if( allocated(meta_reprojtiles) ) call add_state_tiles(meta_states_vol3D(i_state)%get_state())
          endif
          call self%json%add(json_states_ptr, json_state_vol_ptr)
        enddo
        if( l_add ) then
          call self%json%add(json_ptr, json_states_ptr)
        else
          call self%json%destroy(json_states_ptr)
        end if
      endif
    endif
    call self%commit_section(json_ptr, SECTION_SOLVE3D_MULTISTATE)

  contains

    ! the assigned tiles of @p state, as the 'reprojtiles' array of its state_volumes entry
    subroutine add_state_tiles( state )
      integer, intent(in)       :: state
      type(json_value), pointer :: json_tiles_ptr
      logical                   :: l_tiles
      integer                   :: i_tile
      l_tiles = .false.
      call self%json%create_array(json_tiles_ptr, 'reprojtiles')
      do i_tile = 1, size(meta_reprojtiles)
        if( .not. meta_reprojtiles(i_tile)%assigned() ) cycle
        if( meta_reprojtiles(i_tile)%get_idx() /= state ) cycle
        l_tiles = .true.
        call self%json%add(json_tiles_ptr, meta_reprojtiles(i_tile)%jsonise())
      enddo
      if( l_tiles ) then
        call self%json%add(json_state_vol_ptr, json_tiles_ptr)
      else
        call self%json%destroy(json_tiles_ptr)
      endif
    end subroutine add_state_tiles

  end subroutine assemble_stream_solve3D_multistate

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
    call self%json%print_to_string_fast(self%json_root, buffer)
    str = buffer
    if( allocated(buffer) ) deallocate(buffer)
  end function to_string

  ! Return .true. if the JSON root pointer is associated (i.e. new() has been called).
  function is_associated( self ) result( assoc )
    class(gui_assembler), intent(in) :: self
    logical                          :: assoc
    assoc = associated(self%json_root)
  end function is_associated

  !---------------- section helpers ----------------

  ! Removes section @p key from the document and opens its replacement, @p meta's JSON renamed
  ! @p key, in @p json_ptr; .false. when @p meta is not assigned, and the section is then left out.
  function open_section( self, json_ptr, key, meta ) result( l_open )
    class(gui_assembler),     intent(inout) :: self
    type(json_value),         pointer       :: json_ptr
    character(len=*),         intent(in)    :: key
    class(gui_metadata_base), intent(in)    :: meta
    logical                                 :: l_open
    call self%json%remove_if_present(self%json_root, key)
    json_ptr => meta%jsonise()
    l_open = associated(json_ptr)
    if( l_open ) call self%json%rename(json_ptr, key)
  end function open_section

  ! add_assigned adds to @p parent the assigned entries of @p items under @p key, as an array, or as
  ! an object of entries named by their own JSON (@p l_object); left out when none is assigned.
  ! One specific per list type: the entries go to add_item one by one.

  subroutine add_assigned_micrographs( self, parent, key, items, l_object )
    class(gui_assembler),          intent(inout) :: self
    type(json_value),              pointer       :: parent
    character(len=*),              intent(in)    :: key
    type(gui_metadata_micrograph), intent(in)    :: items(:)
    logical,                       intent(in)    :: l_object
    type(json_value),              pointer       :: json_items_ptr
    integer                                      :: i
    call self%open_list(json_items_ptr, key, l_object)
    do i = 1, size(items)
      call self%add_item(json_items_ptr, items(i))
    enddo
    call self%close_list(parent, json_items_ptr)
  end subroutine add_assigned_micrographs

  subroutine add_assigned_histograms( self, parent, key, items, l_object )
    class(gui_assembler),         intent(inout) :: self
    type(json_value),             pointer       :: parent
    character(len=*),             intent(in)    :: key
    type(gui_metadata_histogram), intent(in)    :: items(:)
    logical,                      intent(in)    :: l_object
    type(json_value),             pointer       :: json_items_ptr
    integer                                     :: i
    call self%open_list(json_items_ptr, key, l_object)
    do i = 1, size(items)
      call self%add_item(json_items_ptr, items(i))
    enddo
    call self%close_list(parent, json_items_ptr)
  end subroutine add_assigned_histograms

  subroutine add_assigned_timeplots( self, parent, key, items, l_object )
    class(gui_assembler),        intent(inout) :: self
    type(json_value),            pointer       :: parent
    character(len=*),            intent(in)    :: key
    type(gui_metadata_timeplot), intent(in)    :: items(:)
    logical,                     intent(in)    :: l_object
    type(json_value),            pointer       :: json_items_ptr
    integer                                    :: i
    call self%open_list(json_items_ptr, key, l_object)
    do i = 1, size(items)
      call self%add_item(json_items_ptr, items(i))
    enddo
    call self%close_list(parent, json_items_ptr)
  end subroutine add_assigned_timeplots

  subroutine add_assigned_optics_groups( self, parent, key, items, l_object )
    class(gui_assembler),            intent(inout) :: self
    type(json_value),                pointer       :: parent
    character(len=*),                intent(in)    :: key
    type(gui_metadata_optics_group), intent(in)    :: items(:)
    logical,                         intent(in)    :: l_object
    type(json_value),                pointer       :: json_items_ptr
    integer                                        :: i
    call self%open_list(json_items_ptr, key, l_object)
    do i = 1, size(items)
      call self%add_item(json_items_ptr, items(i))
    enddo
    call self%close_list(parent, json_items_ptr)
  end subroutine add_assigned_optics_groups

  subroutine add_assigned_cavgs2D( self, parent, key, items, l_object )
    class(gui_assembler),      intent(inout) :: self
    type(json_value),          pointer       :: parent
    character(len=*),          intent(in)    :: key
    type(gui_metadata_cavg2D), intent(in)    :: items(:)
    logical,                   intent(in)    :: l_object
    type(json_value),          pointer       :: json_items_ptr
    integer                                  :: i
    call self%open_list(json_items_ptr, key, l_object)
    do i = 1, size(items)
      call self%add_item(json_items_ptr, items(i))
    enddo
    call self%close_list(parent, json_items_ptr)
  end subroutine add_assigned_cavgs2D

  ! A new list @p key: an array, or an object (@p l_object)
  subroutine open_list( self, json_items_ptr, key, l_object )
    class(gui_assembler), intent(inout) :: self
    type(json_value),     pointer       :: json_items_ptr
    character(len=*),     intent(in)    :: key
    logical,              intent(in)    :: l_object
    if( l_object ) then
      call self%json%create_object(json_items_ptr, key)
    else
      call self%json%create_array(json_items_ptr, key)
    endif
  end subroutine open_list

  ! @p item's JSON added to the list, when it is assigned. The polymorphic call's pointer result is
  ! taken into a local pointer before it is passed on (simple-modern-fortran skill).
  subroutine add_item( self, json_items_ptr, item )
    class(gui_assembler),     intent(inout) :: self
    type(json_value),         pointer       :: json_items_ptr
    class(gui_metadata_base), intent(in)    :: item
    type(json_value),         pointer       :: json_item_ptr
    if( .not. item%assigned() ) return
    json_item_ptr => item%jsonise()
    call self%json%add(json_items_ptr, json_item_ptr)
  end subroutine add_item

  ! The list added to @p parent when it holds an entry, destroyed otherwise
  subroutine close_list( self, parent, json_items_ptr )
    class(gui_assembler), intent(inout) :: self
    type(json_value),     pointer       :: parent, json_items_ptr
    if( self%json%count(json_items_ptr) > 0 ) then
      call self%json%add(parent, json_items_ptr)
    else
      call self%json%destroy(json_items_ptr)
    endif
  end subroutine close_list

  ! Adds the opened section @p json_ptr to the document when its text differs from the last one
  ! sent for @p section, and destroys it otherwise.
  subroutine commit_section( self, json_ptr, section )
    class(gui_assembler),     intent(inout) :: self
    type(json_value),         pointer       :: json_ptr
    integer,                  intent(in)    :: section
    character(kind=CK,len=:), allocatable   :: buffer
    type(string)                            :: str, hash
    call self%json%print_to_string_fast(json_ptr, buffer)
    str  = buffer
    hash = str%to_fnv1a_hash64()
    if( hash /= self%hashes(section) ) then
      call self%json%add(self%json_root, json_ptr)
      call self%hashes(section)%kill()
      self%hashes(section) = hash
    else
      call self%json%destroy(json_ptr)
    endif
    if( allocated(buffer) ) deallocate(buffer)
    nullify(json_ptr)
  end subroutine commit_section

end module simple_gui_assembler
