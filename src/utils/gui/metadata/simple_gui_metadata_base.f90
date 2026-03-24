!@descr: Abstract base type for GUI metadata objects.
!==============================================================================
! MODULE: simple_gui_metadata_base
!
! PURPOSE:
!   Defines gui_metadata_base, the common root for all GUI metadata types.
!   Provides:
!     new/kill    — lifecycle management with an integer subtype tag
!     type        — return the subtype tag set at construction
!     initialized — query initialisation state (safe on uninitialised objects)
!     assigned    — query whether derived-type data has been populated
!     serialise   — binary transfer of base-type fields into a character buffer
!     jsonise     — polymorphic JSON hook; emits {} when assigned, else null
!   Derived types extend gui_metadata_base and override jsonise() to emit
!   their own fields.  Subtype tags are defined in simple_gui_metadata_types.
!
! DEPENDENCIES:
!   json_module, simple_error
!==============================================================================
module simple_gui_metadata_base
use json_module,         only: json_core, json_value
use simple_error,        only: simple_exception

implicit none

public :: gui_metadata_base
private
#include "simple_local_flags.inc"

type :: gui_metadata_base
  integer :: meta_type     = 0       ! integer tag identifying the concrete subtype (see simple_gui_metadata_types)
  logical :: l_initialized = .false. ! set by new(), cleared by kill()
  logical :: l_assigned    = .false. ! set by derived-type setters once data has been populated
contains
  procedure :: new
  procedure :: type
  procedure :: assigned
  procedure :: initialized
  procedure :: jsonise
  procedure :: serialise
  procedure :: kill
end type gui_metadata_base

contains

  !---------------- lifecycle ----------------

  ! Initialise the object with a subtype tag.
  subroutine new( self, meta_type )
    class(gui_metadata_base), intent(inout) :: self
    integer,                  intent(in)    :: meta_type
    self%meta_type     = meta_type
    self%l_initialized = .true.
  end subroutine new

  ! Reset all fields and mark the object as uninitialised.
  subroutine kill( self )
    class(gui_metadata_base), intent(inout) :: self
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%meta_type     = 0
    self%l_assigned    = .false.
    self%l_initialized = .false.
  end subroutine kill

  !---------------- queries ----------------

  ! Return the subtype tag set at construction (see simple_gui_metadata_types).
  function type( self ) result( meta_type )
    class(gui_metadata_base), intent(in) :: self
    integer                              :: meta_type
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    meta_type = self%meta_type
  end function type

  ! Return .true. if new() has been called and kill() has not yet been called.
  function initialized( self ) result( l_initialized )
    class(gui_metadata_base), intent(in) :: self
    logical                              :: l_initialized
    l_initialized = self%l_initialized
  end function initialized

  ! Return .true. if derived-type data has been populated via a setter.
  function assigned( self ) result( l_assigned )
    class(gui_metadata_base), intent(in) :: self
    logical                              :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
  end function assigned

  !---------------- serialisation ----------------

  ! Copy base-type fields into buffer as a raw binary transfer.
  subroutine serialise( self, buffer )
    class(gui_metadata_base),              intent(in)    :: self
    character(len=:),         allocatable, intent(inout) :: buffer
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( allocated(buffer) ) deallocate(buffer)
    allocate(character(len=sizeof(self)) :: buffer)
    buffer = transfer(self, buffer)
  end subroutine serialise

  ! Return a JSON object pointer; emits {} when assigned, null otherwise.
  function jsonise( self ) result( json_ptr )
    class(gui_metadata_base), intent(inout) :: self
    type(json_core)                         :: json
    type(json_value),         pointer       :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
    else
      nullify(json_ptr)
    endif
  end function jsonise

end module simple_gui_metadata_base
