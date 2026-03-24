!@descr: GUI metadata type for a time-series plot with one or two data traces.
!==============================================================================
! MODULE: simple_gui_metadata_timeplot
!
! PURPOSE:
!   Extends gui_metadata_base with time-plot fields:
!     name    — display name for the plot
!     labels  — x-axis values (real, up to 512 points)
!     data    — primary y-axis trace (real, same length as labels)
!     data2   — optional secondary y-axis trace (real, same length as labels)
!   Provides set/get for all fields and a jsonise override that emits the
!   name, labels, data, and data2 arrays as a JSON object.
!
! DEPENDENCIES:
!   json_module, simple_core_module_api, simple_gui_metadata_base
!==============================================================================
module simple_gui_metadata_timeplot
use json_kinds
use json_module
use simple_core_module_api
use simple_gui_metadata_base

implicit none

public :: gui_metadata_timeplot
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_timeplot
  private
  character(len=SHORTSTRLEN) :: name     = ''   ! display name for the plot
  real                       :: labels(512)     ! x-axis values
  real                       :: data(512)       ! primary y-axis trace
  real                       :: data2(512)      ! secondary y-axis trace (optional)
  integer                    :: n_labels = 0    ! number of populated points
contains
  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override
end type gui_metadata_timeplot

contains

  !---------------- setters ----------------

  ! Set the plot name, x-axis labels, and one or two data traces.
  ! data2 is optional; omitting it zeros the secondary trace.
  ! All supplied arrays must be the same length and must not exceed 512 points.
  subroutine set( self, name, labels, data, data2 )
    class(gui_metadata_timeplot),              intent(inout) :: self
    type(string),                              intent(in)    :: name
    real,                         allocatable, intent(in)    :: labels(:), data(:)
    real,                         allocatable, intent(in), optional :: data2(:)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( size(labels) /= size(data) ) THROW_HARD('labels and data differ in size')
    if( size(labels) > size(self%labels) ) THROW_HARD('labels exceeds maximum timeplot size')
    self%l_assigned              = .true.
    self%name                    = name%to_char()
    self%n_labels                = size(labels)
    self%labels(1:self%n_labels) = labels
    self%data(1:self%n_labels)   = data
    self%data2                   = 0.0
    if( present(data2) ) then
      if( allocated(data2) ) self%data2(1:self%n_labels) = data2(1:self%n_labels)
    endif
  end subroutine set

  !---------------- getters ----------------

  ! Return the plot name, labels, and both data traces; result is .true. if assigned.
  function get( self, name, labels, data, data2 ) result( l_assigned )
    class(gui_metadata_timeplot),              intent(in)  :: self
    type(string),                              intent(out) :: name
    real,                         allocatable, intent(out) :: labels(:), data(:), data2(:)
    logical                                               :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    name       = trim(self%name)
    allocate(labels(self%n_labels), data(self%n_labels), data2(self%n_labels))
    labels = self%labels(1:self%n_labels)
    data   = self%data(1:self%n_labels)
    data2  = self%data2(1:self%n_labels)
  end function get

  !---------------- serialisation ----------------

  ! Emit name, labels, data, and data2 arrays as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_timeplot), intent(inout) :: self
    type(json_core)                             :: json
    type(json_value),             pointer       :: json_ptr, json_labels_ptr, json_data_ptr, json_data2_ptr
    integer                                     :: i
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr,       trim(self%name))
    call json%create_array(json_labels_ptr, 'labels')
    call json%create_array(json_data_ptr,   'data')
    call json%create_array(json_data2_ptr,  'data2')
    do i = 1, self%n_labels
      call json%add(json_labels_ptr, '', dble(self%labels(i)))
      call json%add(json_data_ptr,   '', dble(self%data(i)))
      call json%add(json_data2_ptr,  '', dble(self%data2(i)))
    end do
    call json%add(json_ptr, json_labels_ptr)
    call json%add(json_ptr, json_data_ptr)
    call json%add(json_ptr, json_data2_ptr)
  end function jsonise_override

end module simple_gui_metadata_timeplot
