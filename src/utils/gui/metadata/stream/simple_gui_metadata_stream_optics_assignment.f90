!@descr: GUI metadata for the stream optics-assignment stage — micrograph and optics-group assignment counts
!==============================================================================
! MODULE: simple_gui_metadata_stream_optics_assignment
!
! PURPOSE:
!   Extends gui_metadata_base with fields specific to the optics-group
!   assignment stage of the cryo-EM streaming pipeline.  Tracks the number
!   of micrographs assigned to optics groups, the number of distinct optics
!   groups created, and the Unix timestamp of the most recently imported
!   micrograph (supplied by the caller).
!
! TYPES:
!   gui_metadata_stream_optics_assignment — extends gui_metadata_base
!     set()     — assign all fields from caller-supplied values
!     get()     — retrieve all fields; returns the l_assigned flag
!     jsonise() — serialise all fields to a json_value tree (base override)
!
! DEPENDENCIES:
!   json_kinds, json_module, simple_string, simple_defs, simple_error,
!   simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_stream_optics_assignment
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
    integer               :: micrographs_assigned     = 0  ! micrographs successfully placed in an optics group
    integer               :: optics_groups_assigned   = 0  ! number of distinct optics groups created
    integer               :: last_micrograph_imported = 0  ! Unix timestamp of most recent import (caller-supplied)
  contains
    procedure :: set
    procedure :: get
    procedure :: jsonise => jsonise_override
  end type gui_metadata_stream_optics_assignment

contains

  ! Assign all fields from caller-supplied values.
  subroutine set( self, stage, micrographs_assigned, optics_groups_assigned, last_micrograph_imported )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    type(string),                                 intent(in)    :: stage
    integer,                                      intent(in)    :: micrographs_assigned, optics_groups_assigned
    integer,                                      intent(in)    :: last_micrograph_imported
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned               = .true.
    self%stage                    = stage%to_char()
    self%micrographs_assigned     = micrographs_assigned
    self%optics_groups_assigned   = optics_groups_assigned
    self%last_micrograph_imported = last_micrograph_imported
  end subroutine set

  ! Retrieve all fields. Returns .true. if the object has been assigned.
  function get( self, stage, micrographs_assigned, optics_groups_assigned, last_micrograph_imported ) result( l_assigned )
    class(gui_metadata_stream_optics_assignment), intent(inout) :: self
    type(string),                                 intent(out)   :: stage
    integer,                                      intent(out)   :: micrographs_assigned, optics_groups_assigned
    integer,                                      intent(out)   :: last_micrograph_imported
    logical                                                     :: l_assigned
    if( .not. self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned               = self%l_assigned
    stage                    = trim(self%stage)
    micrographs_assigned     = self%micrographs_assigned
    optics_groups_assigned   = self%optics_groups_assigned
    last_micrograph_imported = self%last_micrograph_imported
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
      call json%add(json_ptr, 'micrographs_assigned',     self%micrographs_assigned      )
      call json%add(json_ptr, 'optics_groups_assigned',   self%optics_groups_assigned    )
      call json%add(json_ptr, 'last_micrograph_imported', self%last_micrograph_imported  )
    else
      nullify(json_ptr)
    end if
  end function jsonise_override

end module simple_gui_metadata_stream_optics_assignment
