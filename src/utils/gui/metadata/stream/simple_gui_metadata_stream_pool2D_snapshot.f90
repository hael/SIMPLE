!@descr: GUI metadata for a pool-2D snapshot — id, filename, particle count, and timestamp
!==============================================================================
! MODULE: simple_gui_metadata_stream_pool2D_snapshot
!
! PURPOSE:
!   Extends gui_metadata_base with the fields the pool-2D process sends to
!   the GUI after writing a classification snapshot.  Carries the snapshot
!   id, project filename, the number of particles included, and a Unix timestamp
!   recording when the snapshot was created.
!
! TYPES:
!   gui_metadata_stream_pool2D_snapshot — extends gui_metadata_base
!     set()     — assign id, filename and particle count; records snapshot_time
!                 automatically
!     get()     — retrieve all four fields; returns the l_assigned flag
!     jsonise() — serialise all fields to a json_value tree (base override)
!
! DEPENDENCIES:
!   unix, json_module, simple_defs, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_pool2D_snapshot
  use unix,                     only: c_long, c_time
  use simple_error,             only: simple_exception
  use json_module,              only: json_core, json_value
  use simple_defs,              only: LONGSTRLEN
  use simple_string,            only: string
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_pool2D_snapshot
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_pool2D_snapshot
    private
    character(len=LONGSTRLEN) :: snapshot_filename = '' ! project file name for the snapshot
    integer                   :: id                = 0
    integer                   :: snapshot_nptcls   = 0  ! number of particles in the snapshot
    integer                   :: snapshot_time     = 0  ! Unix timestamp when the snapshot was written
  contains
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_pool2D_snapshot

contains

  ! Assign the snapshot id, filename, and particle count.
  ! Records snapshot_time automatically from the current Unix clock.
  subroutine set( self, id, snapshot_filename, snapshot_nptcls )
    class(gui_metadata_stream_pool2D_snapshot), intent(inout) :: self
    integer,                                    intent(in)    :: id
    type(string),                               intent(in)    :: snapshot_filename
    integer,                                    intent(in)    :: snapshot_nptcls
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( snapshot_filename%strlen() > LONGSTRLEN ) THROW_HARD('snapshot_filename exceeds buffer length')
    self%l_assigned        = .true.
    self%id                = id
    self%snapshot_filename = snapshot_filename%to_char()
    self%snapshot_nptcls   = snapshot_nptcls
    self%snapshot_time     = int(c_time(0_c_long))
  end subroutine set

  ! Retrieve the snapshot id, filename, particle count, and timestamp.
  ! Returns .true. if the object has been assigned.
  function get( self, id, snapshot_filename, snapshot_nptcls, snapshot_time ) result( l_assigned )
    class(gui_metadata_stream_pool2D_snapshot), intent(in)  :: self
    integer,                                    intent(out) :: id
    type(string),                               intent(out) :: snapshot_filename
    integer,                                    intent(out) :: snapshot_nptcls
    integer,                                    intent(out) :: snapshot_time
    logical :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned        = self%l_assigned
    id                = self%id
    snapshot_filename = trim(self%snapshot_filename)
    snapshot_nptcls   = self%snapshot_nptcls
    snapshot_time     = self%snapshot_time
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_pool2D_snapshot), intent(inout) :: self
    type(json_core)                                           :: json
    type(json_value),                           pointer       :: json_ptr
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'id',                self%id                     )
      call json%add(json_ptr, 'snapshot_filename', trim(self%snapshot_filename))
      call json%add(json_ptr, 'snapshot_nptcls',   self%snapshot_nptcls        )
      call json%add(json_ptr, 'snapshot_time',     self%snapshot_time          )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_stream_pool2D_snapshot
