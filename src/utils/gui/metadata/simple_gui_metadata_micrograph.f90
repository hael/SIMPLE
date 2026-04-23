!@descr: GUI metadata type for a single micrograph and its particle coordinates.
!==============================================================================
! MODULE: simple_gui_metadata_micrograph
!
! PURPOSE:
!   Extends gui_metadata_base with micrograph-specific fields:
!     path        — absolute path to the micrograph file
!     dfx/dfy     — defocus values (Angstroms)
!     ctfres      — CTF resolution estimate (Angstroms)
!     i/i_max     — micrograph index within the current batch
!     xdim/ydim   — micrograph dimensions in pixels
!     coordinates — up to 1500 particle box centres (int16 x/y pairs)
!   Provides set/get for scalar fields, set_coordinate/clear_coordinates for
!   the coordinate array, and a jsonise override that emits all fields plus
!   an optional "boxes" array when coordinates are present.
!
! DEPENDENCIES:
!   json_module, simple_defs, simple_string, simple_error, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_micrograph
use json_module,              only: json_core, json_value
use simple_defs,              only: LONGSTRLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_micrograph
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_micrograph
  private
  character(len=LONGSTRLEN) :: path          = ''    ! absolute path to the micrograph file
  real                      :: dfx           = 0.0   ! defocus value along x (Angstroms)
  real                      :: dfy           = 0.0   ! defocus value along y (Angstroms)
  real                      :: ctfres        = 0.0   ! CTF resolution estimate (Angstroms)
  integer                   :: i             = 1     ! index of this micrograph within the current batch
  integer                   :: i_max         = 1     ! total micrographs in the current batch
  integer                   :: xdim          = 0     ! micrograph width in pixels
  integer                   :: ydim          = 0     ! micrograph height in pixels
  integer                   :: n_coordinates = 0     ! number of populated coordinate entries
  integer(kind=2)           :: x_coordinates(1500)   ! particle box centre x (int16, limits transfer size)
  integer(kind=2)           :: y_coordinates(1500)   ! particle box centre y (int16, limits transfer size)
contains
  procedure :: set
  procedure :: set_coordinate
  procedure :: clear_coordinates
  procedure :: get
  procedure :: get_i
  procedure :: get_i_max
  procedure :: jsonise => jsonise_override
end type gui_metadata_micrograph

contains

  !---------------- setters ----------------

  ! Set scalar micrograph fields and mark the object as assigned.
  subroutine set( self, path, dfx, dfy, ctfres, i, i_max )
    class(gui_metadata_micrograph), intent(inout) :: self
    type(string),                   intent(in)    :: path
    real,                           intent(in)    :: dfx, dfy, ctfres
    integer,                        intent(in)    :: i, i_max
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%l_assigned = .true.
    self%path       = path%to_char()
    self%dfx        = dfx
    self%dfy        = dfy
    self%ctfres     = ctfres
    self%i          = i
    self%i_max      = i_max
  end subroutine set

  ! Store a single particle coordinate and update the micrograph dimensions.
  subroutine set_coordinate( self, i_coord, x, y, xdim, ydim )
    class(gui_metadata_micrograph), intent(inout) :: self
    integer,                        intent(in)    :: i_coord, x, y, xdim, ydim
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( i_coord < 1 ) THROW_HARD('i_coord is out of range')
    if( i_coord > size(self%x_coordinates) ) THROW_HARD('i_coord is out of range')
    self%l_assigned = .true.
    self%xdim       = xdim
    self%ydim       = ydim
    if( i_coord > self%n_coordinates ) self%n_coordinates = i_coord
    self%x_coordinates(i_coord) = int2(x)
    self%y_coordinates(i_coord) = int2(y)
  end subroutine set_coordinate

  ! Zero all coordinate fields and reset the coordinate count.
  subroutine clear_coordinates( self )
    class(gui_metadata_micrograph), intent(inout) :: self
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    self%n_coordinates = 0
    self%x_coordinates = 0
    self%y_coordinates = 0
    self%xdim          = 0
    self%ydim          = 0
  end subroutine clear_coordinates

  !---------------- getters ----------------

  ! Return all scalar fields; result is .true. if the object has been assigned.
  function get( self, path, dfx, dfy, ctfres, i, i_max ) result( l_assigned )
    class(gui_metadata_micrograph), intent(in)  :: self
    type(string),                   intent(out) :: path
    real,                           intent(out) :: dfx, dfy, ctfres
    integer,                        intent(out) :: i, i_max
    logical                                     :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    path       = trim(self%path)
    dfx        = self%dfx
    dfy        = self%dfy
    ctfres     = self%ctfres
    i          = self%i
    i_max      = self%i_max
  end function get

  ! Return the micrograph index within the current batch.
  function get_i( self ) result( i )
    class(gui_metadata_micrograph), intent(in) :: self
    integer                                    :: i
    i = self%i
  end function get_i

  ! Return the total number of micrographs in the current batch.
  function get_i_max( self ) result( i_max )
    class(gui_metadata_micrograph), intent(in) :: self
    integer                                    :: i_max
    i_max = self%i_max
  end function get_i_max

  !---------------- serialisation ----------------

  ! Emit all fields as a JSON object; appends a "boxes" array when coordinates
  ! are present.  Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_micrograph), intent(inout) :: self
    type(json_core)                               :: json
    type(json_value),               pointer       :: json_ptr, json_boxes_ptr, json_coords_ptr
    integer                                       :: i_coord
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr, '')
    call json%add(json_ptr, "path",   trim(self%path))
    call json%add(json_ptr, "dfx",    dble(self%dfx))
    call json%add(json_ptr, "dfy",    dble(self%dfy))
    call json%add(json_ptr, "ctfres", dble(self%ctfres))
    call json%add(json_ptr, "i",      self%i)
    call json%add(json_ptr, "i_max",  self%i_max)
    if( self%n_coordinates > 0 ) then
      call json%add(json_ptr, "xdim", self%xdim)
      call json%add(json_ptr, "ydim", self%ydim)
      call json%create_array(json_boxes_ptr, 'boxes')
      do i_coord = 1, self%n_coordinates
        call json%create_object(json_coords_ptr, '')
        call json%add(json_coords_ptr, "x", int(self%x_coordinates(i_coord)))
        call json%add(json_coords_ptr, "y", int(self%y_coordinates(i_coord)))
        call json%add(json_boxes_ptr, json_coords_ptr)
      end do
      call json%add(json_ptr, json_boxes_ptr)
    endif
  end function jsonise_override

end module simple_gui_metadata_micrograph
