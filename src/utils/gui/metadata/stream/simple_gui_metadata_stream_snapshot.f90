!@descr: GUI metadata for a stream snapshot report — id, filename, particle count, timestamp and, for 3D, the states
! Sent by stream p06 (a 2D snapshot of selected classes) and p07 (a 3D snapshot of selected states,
! merged into one) after it writes a snapshot project, or decides not to; the message's tag tells
! which. set() stamps snapshot_time with the current Unix time. states holds a 3D snapshot's
! selected states (up to MAX_STATES_SOLVE3D_MULTISTATE), emitted only when set.
module simple_gui_metadata_stream_snapshot
  use unix,                     only: c_long, c_time
  use simple_error,             only: simple_exception
  use json_module,              only: json_core, json_value
  use simple_defs,              only: LONGSTRLEN
  use simple_string,            only: string
  use simple_gui_metadata_base, only: gui_metadata_base
  use simple_gui_metadata_stream_solve3D_multistate, only: MAX_STATES_SOLVE3D_MULTISTATE

  implicit none

  public :: gui_metadata_stream_snapshot
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_snapshot
    private
    character(len=LONGSTRLEN) :: snapshot_filename = '' ! project file name for the snapshot
    integer                   :: id                = 0
    integer                   :: snapshot_nptcls   = 0  ! number of particles in the snapshot
    integer                   :: snapshot_time     = 0  ! Unix timestamp when the snapshot was written
    integer                   :: nstates           = 0  ! number of selected states (3D)
    integer                   :: states(MAX_STATES_SOLVE3D_MULTISTATE) = 0 ! the selected states (3D)
  contains
    procedure :: kill => kill_override
    procedure :: set
    procedure :: get
    procedure :: get_states
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_snapshot

contains

  ! Assign the snapshot id, filename, particle count and, for a 3D snapshot, its selected states.
  ! Records snapshot_time automatically from the current Unix clock.
  subroutine set( self, id, snapshot_filename, snapshot_nptcls, states )
    class(gui_metadata_stream_snapshot), intent(inout) :: self
    integer,                             intent(in)    :: id
    type(string),                        intent(in)    :: snapshot_filename
    integer,                             intent(in)    :: snapshot_nptcls
    integer, optional,                   intent(in)    :: states(:)
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( snapshot_filename%strlen() > LONGSTRLEN ) THROW_HARD('snapshot_filename exceeds buffer length')
    self%l_assigned        = .true.
    self%id                = id
    self%snapshot_filename = snapshot_filename%to_char()
    self%snapshot_nptcls   = snapshot_nptcls
    self%snapshot_time     = int(c_time(0_c_long))
    self%nstates           = 0
    self%states            = 0
    if( present(states) )then
      if( size(states) > MAX_STATES_SOLVE3D_MULTISTATE ) THROW_HARD('states exceed MAX_STATES_SOLVE3D_MULTISTATE')
      self%nstates             = size(states)
      self%states(1:size(states)) = states
    endif
  end subroutine set

  ! Retrieve the snapshot id, filename, particle count, and timestamp.
  ! Returns .true. if the object has been assigned.
  function get( self, id, snapshot_filename, snapshot_nptcls, snapshot_time ) result( l_assigned )
    class(gui_metadata_stream_snapshot), intent(in)  :: self
    integer,                             intent(out) :: id
    type(string),                        intent(out) :: snapshot_filename
    integer,                             intent(out) :: snapshot_nptcls
    integer,                             intent(out) :: snapshot_time
    logical :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned        = self%l_assigned
    id                = self%id
    snapshot_filename = trim(self%snapshot_filename)
    snapshot_nptcls   = self%snapshot_nptcls
    snapshot_time     = self%snapshot_time
  end function get

  ! The selected states of a 3D snapshot; none for a 2D one.
  function get_states( self ) result( states )
    class(gui_metadata_stream_snapshot), intent(in) :: self
    integer, allocatable :: states(:)
    states = self%states(1:self%nstates)
  end function get_states

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_snapshot), intent(in) :: self
    type(json_core)                                 :: json
    type(json_value),                    pointer    :: json_ptr, json_states_ptr
    integer                                         :: i
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'id',                self%id                     )
      call json%add(json_ptr, 'snapshot_filename', trim(self%snapshot_filename))
      call json%add(json_ptr, 'snapshot_nptcls',   self%snapshot_nptcls        )
      call json%add(json_ptr, 'snapshot_time',     self%snapshot_time          )
      if( self%nstates > 0 ) then
        call json%create_array(json_states_ptr, 'states')
        do i = 1, self%nstates
          call json%add(json_states_ptr, '', self%states(i))
        end do
        call json%add(json_ptr, json_states_ptr)
      end if
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_snapshot), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_snapshot )
        self = gui_metadata_stream_snapshot()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_snapshot
