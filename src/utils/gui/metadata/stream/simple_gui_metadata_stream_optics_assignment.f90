!@descr: GUI metadata for the stream optics-assignment stage — micrograph and optics-group assignment counts
! Filled by stream p02. set() stamps last_import_time (Unix time) whenever micrographs_imported
! changes.
module simple_gui_metadata_stream_optics_assignment
  use unix,                     only: c_long, c_time
  use json_kinds
  use json_module,              only: json_core, json_value
  use simple_string,            only: string
  use simple_defs,              only: STDLEN
  use simple_error,             only: simple_exception
  use simple_gui_metadata_base, only: gui_metadata_base

  implicit none

  public :: gui_metadata_stream_optics_assignment
  private
#include "simple_local_flags.inc"

  type, extends(gui_metadata_base) :: gui_metadata_stream_optics_assignment
    private
    character(len=STDLEN) :: stage                    = 'unknown'
    integer               :: micrographs_imported     = 0  ! total micrographs received from import
    integer               :: micrographs_assigned     = 0  ! micrographs successfully placed in an optics group
    integer               :: optics_groups_assigned   = 0  ! number of distinct optics groups created
    integer               :: last_import_time         = 0  ! Unix timestamp of most recent import (stamped by set)
  contains
    procedure :: kill => kill_override
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_optics_assignment

contains

  ! Assign the counts and stage; last_import_time is stamped when micrographs_imported changes.
  subroutine set( self, stage, micrographs_assigned, optics_groups_assigned, micrographs_imported )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    type(string),                                 intent(in)    :: stage
    integer,                                      intent(in)    :: micrographs_assigned, optics_groups_assigned
    integer,                                     intent(in)    :: micrographs_imported
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned               = .true.
    self%stage                    = stage%to_char()
    if( micrographs_imported /= self%micrographs_imported ) &
      self%last_import_time = int(c_time(0_c_long))
    self%micrographs_imported     = micrographs_imported
    self%micrographs_assigned     = micrographs_assigned
    self%optics_groups_assigned   = optics_groups_assigned
  end subroutine set

  ! Retrieve all fields. Returns .true. if the object has been assigned.
  function get( self, stage, micrographs_assigned, optics_groups_assigned, last_import_time, micrographs_imported ) result( l_assigned )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    type(string),                                 intent(out)   :: stage
    integer,                                      intent(out)   :: micrographs_assigned, optics_groups_assigned
    integer,                                      intent(out)   :: last_import_time
    integer,                                      intent(out)   :: micrographs_imported
    logical                                                     :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned               = self%l_assigned
    stage                    = trim(self%stage)
    micrographs_assigned     = self%micrographs_assigned
    optics_groups_assigned   = self%optics_groups_assigned
    last_import_time         = self%last_import_time
    micrographs_imported     = self%micrographs_imported
  end function get

  ! Serialise all fields to a JSON object. Returns a null pointer when
  ! the object has not yet been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    type(json_core)                                             :: json
    type(json_value),                             pointer       :: json_ptr
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, 'stage',                    trim(self%stage)               )
      call json%add(json_ptr, 'micrographs_imported',     self%micrographs_imported      )
      call json%add(json_ptr, 'micrographs_assigned',     self%micrographs_assigned      )
      call json%add(json_ptr, 'optics_groups_assigned',   self%optics_groups_assigned    )
      call json%add(json_ptr, 'last_import_time',         self%last_import_time          )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    select type( self )
      type is( gui_metadata_stream_optics_assignment )
        self = gui_metadata_stream_optics_assignment()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_stream_optics_assignment
