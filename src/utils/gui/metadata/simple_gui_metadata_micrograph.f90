!@descr: gui metadata micrograph structure
module simple_gui_metadata_micrograph
use json_kinds
use json_module
use simple_core_module_api
use simple_gui_metadata_base

implicit none

public :: gui_metadata_micrograph
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_micrograph
  private
  character(len=LONGSTRLEN) :: path   = ''
  real                      :: dfx    = 0.0
  real                      :: dfy    = 0.0
  real                      :: ctfres = 0.0

  contains

  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override

end type gui_metadata_micrograph

contains

  subroutine set( self, path, dfx, dfy, ctfres )
    class(gui_metadata_micrograph), intent(inout) :: self
    type(string),                   intent(in)    :: path 
    real,                           intent(in)    :: dfx, dfy, ctfres
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%path       = path%to_char()
    self%dfx        = dfx
    self%dfy        = dfy
    self%ctfres     = ctfres
  end subroutine set

  function get( self, path, dfx, dfy, ctfres ) result( l_assigned )
    class(gui_metadata_micrograph), intent(inout) :: self
    type(string),                   intent(out)   :: path 
    real,                           intent(out)   :: dfx, dfy, ctfres
    logical                                       :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    path       = trim(self%path)
    dfx        = self%dfx
    dfy        = self%dfy
    ctfres     = self%ctfres
  end function get

  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_micrograph), intent( inout ) :: self
    type(json_core)                                 :: json
    type(json_value),               pointer         :: json_ptr
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, "path",   trim(self%path)  )
      call json%add(json_ptr, "dfx",    dble(self%dfx)   )
      call json%add(json_ptr, "dfy",    dble(self%dfy)   )
      call json%add(json_ptr, "ctfres", dble(self%ctfres))
    endif   
  end function jsonise_override

end module simple_gui_metadata_micrograph