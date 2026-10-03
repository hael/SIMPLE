!@descr: GUI metadata for the stream multistate solve3D stage — pipeline stage, particle/state counts, per-state resolution, and user-input flag
! Filled by stream p07: refine_iteration and particles_at_last_refine count solve3D_addon
! passes and the pool size at the last one. last_import_time is stamped when particles_imported
! changes. p07 always sends resolution=0 and never sets user_input; per-state stats follow once
! solve3D is done.
module simple_gui_metadata_stream_solve3D_multistate
  use unix,                     only: c_long, c_time
  use simple_error,             only: simple_exception
  use json_module,              only: json_core, json_value
  use simple_defs,              only: STDLEN
  use simple_string,            only: string
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_solve3D_multistate
  private
#include "simple_local_flags.inc"

  integer, parameter :: MAX_STATES_SOLVE3D_MULTISTATE = 20

  type, extends(gui_metadata_base) :: gui_metadata_stream_solve3D_multistate
    private

    character(len=STDLEN)  :: stage                     = 'unknown'
    integer                :: solve3D_stage             = 0       ! internal solve3D progress: 0=not started, 1=running, 2=complete
    integer                :: refine_iteration          = 0       ! solve3D_addon passes started
    integer                :: nstates                   = 0       ! number of states being resolved
    integer                :: particles_imported        = 0       ! total particles pooled from upstream
    integer                :: particles_at_last_refine  = 0       ! pool size when solve3D or the last addon pass began
    integer                :: last_import_time          = 0       ! Unix timestamp of most recent import event
    logical                :: user_input                = .false. ! .true. once the user has supplied input
    real                   :: resolution                = 0.0     ! overall current low-pass/resolution estimate
    integer                :: state_populations(MAX_STATES_SOLVE3D_MULTISTATE) = 0
    real                   :: state_resolutions(MAX_STATES_SOLVE3D_MULTISTATE) = 0.0
  contains
    procedure :: set
    procedure :: set_user_input
    procedure :: set_state_stats
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_solve3D_multistate

contains

  ! Assign pipeline stage, progress counters, particle counts, and overall
  ! resolution. Updates last_import_time only when the pool size changes.
  subroutine set( self, stage, solve3D_stage, refine_iteration, nstates, particles_imported, particles_at_last_refine, resolution )
    class(gui_metadata_stream_solve3D_multistate), intent(inout) :: self
    type(string), intent(in) :: stage
    integer,      intent(in) :: solve3D_stage, refine_iteration, nstates
    integer,      intent(in) :: particles_imported, particles_at_last_refine
    real,         intent(in) :: resolution
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned                = .true.
    self%stage                     = stage%to_char()
    self%solve3D_stage             = solve3D_stage
    self%refine_iteration          = refine_iteration
    self%nstates                   = nstates
    if( particles_imported /= self%particles_imported ) &
      self%last_import_time        = int(c_time(0_c_long))
    self%particles_imported        = particles_imported
    self%particles_at_last_refine  = particles_at_last_refine
    self%resolution                = resolution
  end subroutine set

  ! Set the user-input flag. May be called independently of set().
  subroutine set_user_input( self, user_input )
    class(gui_metadata_stream_solve3D_multistate), intent(inout) :: self
    logical,                                           intent(in)    :: user_input
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%user_input = user_input
  end subroutine set_user_input

  ! Assign the pooled particle population and resolution estimate for one state.
  subroutine set_state_stats( self, state, population, resolution )
    class(gui_metadata_stream_solve3D_multistate), intent(inout) :: self
    integer,                                          intent(in)    :: state, population
    real,                                             intent(in)    :: resolution
    if( .not.self%l_initialized )                          THROW_HARD('gui metadata object is uninitialised')
    if( state < 1 .or. state > size(self%state_populations) ) THROW_HARD('state is out of range')
    self%l_assigned                    = .true.
    self%state_populations(state)      = population
    self%state_resolutions(state)      = resolution
  end subroutine set_state_stats

  ! Retrieve pipeline stage, progress counters, particle counts, overall
  ! resolution, user-input flag, and the last-import timestamp. Returns
  ! .true. if the object has been assigned.
  function get( self, stage, solve3D_stage, refine_iteration, nstates, particles_imported, particles_at_last_refine, last_import_time, user_input, resolution ) result( l_assigned )
    class(gui_metadata_stream_solve3D_multistate), intent(in)     :: self
    type(string), intent(out) :: stage
    integer,      intent(out) :: solve3D_stage, refine_iteration, nstates
    integer,      intent(out) :: particles_imported, particles_at_last_refine
    integer,      intent(out) :: last_import_time
    logical,      intent(out) :: user_input
    real,         intent(out) :: resolution
    logical                   :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned                = self%l_assigned
    stage                     = trim(self%stage)
    solve3D_stage             = self%solve3D_stage
    refine_iteration          = self%refine_iteration
    nstates                   = self%nstates
    particles_imported        = self%particles_imported
    particles_at_last_refine  = self%particles_at_last_refine
    last_import_time          = self%last_import_time
    user_input                = self%user_input
    resolution                = self%resolution
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_solve3D_multistate), intent(inout) :: self
    type(json_core)                                   :: json
    type(json_value), pointer                         :: json_ptr, json_states_ptr => null(), json_state_ptr => null()
    integer                                            :: i_state
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',                     trim(self%stage)               )
      call json%add(json_ptr, 'solve3D_stage',             self%solve3D_stage             )
      call json%add(json_ptr, 'refine_iteration',          self%refine_iteration          )
      call json%add(json_ptr, 'nstates',                   self%nstates                   )
      call json%add(json_ptr, 'particles_imported',        self%particles_imported        )
      call json%add(json_ptr, 'particles_at_last_refine',  self%particles_at_last_refine  )
      call json%add(json_ptr, 'last_import_time',          self%last_import_time          )
      call json%add(json_ptr, 'user_input',                self%user_input                )
      call json%add(json_ptr, 'resolution',                dble(self%resolution)          )
      if( self%nstates > 0 ) then
        call json%create_array(json_states_ptr, 'states')
        do i_state = 1, self%nstates
          call json%create_object(json_state_ptr, '')
          call json%add(json_state_ptr, 'state',      i_state                              )
          call json%add(json_state_ptr, 'population', self%state_populations(i_state)      )
          call json%add(json_state_ptr, 'resolution', dble(self%state_resolutions(i_state)) )
          call json%add(json_states_ptr, json_state_ptr)
        end do
        call json%add(json_ptr, json_states_ptr)
      end if
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_stream_solve3D_multistate
