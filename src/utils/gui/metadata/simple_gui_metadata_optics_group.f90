!@descr: GUI metadata type for an optics group and its beam-shift scatter plot.
!==============================================================================
! MODULE: simple_gui_metadata_optics_group
!
! PURPOSE:
!   Extends gui_metadata_base with optics-group fields:
!     i/i_max   — optics-group index within the current batch
!     n_shifts  — number of populated beam-shift entries
!     xshifts   — beam-shift x components (Angstroms, up to max_points)
!     yshifts   — beam-shift y components (Angstroms, up to max_points)
!   Provides set/get for all fields, scalar accessors for i, i_max, and
!   max_points, and a jsonise override that emits a "coordinates" array of
!   {x, y} shift objects.
!
! DEPENDENCIES:
!   json_module, simple_core_module_api, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_optics_group
use json_module,              only: json_core, json_value
use simple_error,             only: simple_exception
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_optics_group
private
integer, parameter :: max_points = 100   ! maximum number of beam-shift entries
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_optics_group
  private
  integer :: i        = 1   ! index of this optics group within the current batch
  integer :: i_max    = 1   ! total optics groups in the current batch
  integer :: n_shifts = 0   ! number of populated beam-shift entries
  real    :: xshifts(max_points)   ! beam-shift x components (Angstroms)
  real    :: yshifts(max_points)   ! beam-shift y components (Angstroms)
contains
  procedure :: set
  procedure :: get
  procedure :: get_i
  procedure :: get_i_max
  procedure :: get_max_points
  procedure :: jsonise => jsonise_override
end type gui_metadata_optics_group

contains

  !---------------- setters ----------------

  ! Set all optics-group fields; n_shifts must not exceed max_points or
  ! the size of the supplied shift arrays.
  subroutine set( self, i, i_max, xshifts, yshifts, n_shifts )
    class(gui_metadata_optics_group), intent(inout) :: self
    integer,                          intent(in)    :: i, i_max, n_shifts
    real,                allocatable, intent(in)    :: xshifts(:), yshifts(:)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( n_shifts > max_points     ) THROW_HARD('n_shifts exceeds max_points')
    if( n_shifts > size(xshifts)  ) THROW_HARD('n_shifts exceeds size of xshifts')
    if( n_shifts > size(yshifts)  ) THROW_HARD('n_shifts exceeds size of yshifts')
    self%l_assigned              = .true.
    self%i                       = i
    self%i_max                   = i_max
    self%n_shifts                = n_shifts
    self%xshifts(1:n_shifts)     = xshifts(1:n_shifts)
    self%yshifts(1:n_shifts)     = yshifts(1:n_shifts)
  end subroutine set

  !---------------- getters ----------------

  ! Return all fields; result is .true. if the object has been assigned.
  function get( self, i, i_max, xshifts, yshifts, n_shifts ) result( l_assigned )
    class(gui_metadata_optics_group), intent(in)  :: self
    integer,                          intent(out) :: i, i_max, n_shifts
    real,                allocatable, intent(out) :: xshifts(:), yshifts(:)
    logical                                       :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    i          = self%i
    i_max      = self%i_max
    n_shifts   = self%n_shifts
    allocate(xshifts(self%n_shifts), yshifts(self%n_shifts))
    xshifts = self%xshifts(1:self%n_shifts)
    yshifts = self%yshifts(1:self%n_shifts)
  end function get

  ! Return the optics-group index within the current batch.
  function get_i( self ) result( i )
    class(gui_metadata_optics_group), intent(in) :: self
    integer                                      :: i
    i = self%i
  end function get_i

  ! Return the total number of optics groups in the current batch.
  function get_i_max( self ) result( i_max )
    class(gui_metadata_optics_group), intent(in) :: self
    integer                                      :: i_max
    i_max = self%i_max
  end function get_i_max

  ! Return the maximum number of beam-shift entries supported.
  function get_max_points( self ) result( n )
    class(gui_metadata_optics_group), intent(in) :: self
    integer                                      :: n
    n = max_points
  end function get_max_points

  !---------------- serialisation ----------------

  ! Emit id and a "coordinates" array of {x, y} shift objects as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_optics_group), intent(inout) :: self
    type(json_core)                                 :: json
    type(json_value),                 pointer       :: json_ptr, json_coords_ptr, json_coords_array_ptr
    integer                                         :: i_coord
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr, '')
    call json%add(json_ptr, "id", self%i)
    call json%create_array(json_coords_array_ptr, 'coordinates')
    do i_coord = 1, self%n_shifts
      call json%create_object(json_coords_ptr, '')
      call json%add(json_coords_ptr, "x", dble(self%xshifts(i_coord)))
      call json%add(json_coords_ptr, "y", dble(self%yshifts(i_coord)))
      call json%add(json_coords_array_ptr, json_coords_ptr)
    end do
    call json%add(json_ptr, json_coords_array_ptr)
  end function jsonise_override

end module simple_gui_metadata_optics_group
