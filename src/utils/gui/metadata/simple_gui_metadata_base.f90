!@descr: gui metadata base structure
module simple_gui_metadata_base
use json_kinds
use json_module
use simple_core_module_api

implicit none

public :: gui_metadata_base
private
#include "simple_local_flags.inc"

type :: gui_metadata_base
  integer               :: meta_type     = 0
  logical               :: l_initialized = .false.
  logical               :: l_assigned    = .false.
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

  subroutine new( self, meta_type )
    class(gui_metadata_base), intent(inout) :: self
    integer,                  intent(in)    :: meta_type 
    self%meta_type     = meta_type
    self%l_initialized = .true.
  end subroutine new

  subroutine kill( self )
    class(gui_metadata_base), intent(inout) :: self
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%meta_type     = 0
    self%l_assigned    = .false.
    self%l_initialized = .false.
  end subroutine kill

  subroutine serialise( self, buffer )
    class(gui_metadata_base),              intent(inout) :: self
    character(len=:),         allocatable, intent(inout) :: buffer
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( allocated(buffer) ) deallocate(buffer)
    allocate(character(len=sizeof(self)) :: buffer)
    buffer = transfer(self, buffer)
  end subroutine serialise

  function jsonise( self ) result( json_ptr )
    class(gui_metadata_base), intent(inout) :: self
    type(json_core)                           :: json
    type(json_value),         pointer         :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) call json%create_object(json_ptr, '')   
  end function jsonise

  function type( self ) result( meta_type )
    class(gui_metadata_base), intent(inout) :: self
    integer                                 :: meta_type
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    meta_type = self%meta_type
  end function type

  function assigned( self ) result( l_assigned )
    class(gui_metadata_base), intent(inout) :: self
    logical                                 :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
  end function assigned

  function initialized( self ) result( l_initialized )
    class(gui_metadata_base), intent(inout) :: self
    logical                                 :: l_initialized
    l_initialized = self%l_initialized
  end function initialized

end module simple_gui_metadata_base