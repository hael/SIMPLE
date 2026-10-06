!@descr: GUI metadata type for a time-series plot with one or two data traces.
! Up to MAX_TIMEPLOT_POINTS points; jsonise emits {labels, data, data2} as an object named after the plot, with
! data2 all zeros when set() got no second trace.
module simple_gui_metadata_timeplot
use json_module,              only: json_core, json_value
use simple_defs,              only: SHORTSTRLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_timeplot, MAX_TIMEPLOT_POINTS
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_TIMEPLOT_POINTS = 512 ! points a plot holds

type, extends( gui_metadata_base ) :: gui_metadata_timeplot
  private
  character(len=SHORTSTRLEN) :: name     = ''   ! display name for the plot
  real                       :: labels(MAX_TIMEPLOT_POINTS) = 0. ! x-axis values
  real                       :: data(MAX_TIMEPLOT_POINTS)   = 0. ! primary y-axis trace
  real                       :: data2(MAX_TIMEPLOT_POINTS)  = 0. ! secondary y-axis trace (optional)
  integer                    :: n_labels = 0    ! number of populated points
contains
  procedure :: kill => kill_override
  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override
end type gui_metadata_timeplot

contains

  !---------------- setters ----------------

  ! Set the plot name, x-axis labels, and one or two data traces.
  ! data2 is optional; omitting it zeros the secondary trace.
  ! All supplied arrays must be the same length and must not exceed MAX_TIMEPLOT_POINTS.
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
      if( allocated(data2) ) then
        if( size(data2) < self%n_labels ) THROW_HARD('data2 is smaller than labels')
        self%data2(1:self%n_labels) = data2(1:self%n_labels)
      end if
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
    class(gui_metadata_timeplot), intent(in)    :: self
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

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_timeplot), intent(inout) :: self
    select type( self )
      type is( gui_metadata_timeplot )
        self = gui_metadata_timeplot()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_timeplot
