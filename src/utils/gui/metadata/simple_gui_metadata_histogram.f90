!@descr: GUI metadata type for a labelled histogram.
! Up to MAX_HISTOGRAM_BINS bins of real labels and integer counts; jsonise emits {labels:[...], data:[...]} as an
! object named after the histogram (name is the key, not a field).
module simple_gui_metadata_histogram
use json_module,              only: json_core, json_value
use simple_defs,              only: SHORTSTRLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_gui_metadata_base, only: gui_metadata_base

implicit none

public :: gui_metadata_histogram, MAX_HISTOGRAM_BINS
private
#include "simple_local_flags.inc"

integer, parameter :: MAX_HISTOGRAM_BINS = 512 ! bins a histogram holds

type, extends( gui_metadata_base ) :: gui_metadata_histogram
  private
  character(len=SHORTSTRLEN) :: name     = ''   ! display name for the histogram
  real                       :: labels(MAX_HISTOGRAM_BINS) = 0. ! bin boundary / centre values
  integer                    :: data(MAX_HISTOGRAM_BINS)   = 0  ! bin counts
  integer                    :: n_labels = 0    ! number of populated bins
contains
  procedure :: kill => kill_override
  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override
end type gui_metadata_histogram

contains

  !---------------- setters ----------------

  ! Set the histogram name, labels, and data arrays; arrays must be the same
  ! length and must not exceed MAX_HISTOGRAM_BINS elements.
  subroutine set( self, name, labels, data )
    class(gui_metadata_histogram),              intent(inout) :: self
    type(string),                               intent(in)    :: name
    real,                          allocatable, intent(in)    :: labels(:)
    integer,                       allocatable, intent(in)    :: data(:)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( size(labels) /= size(data) ) THROW_HARD('labels and data differ in size')
    if( size(labels) > size(self%labels) ) THROW_HARD('labels exceeds maximum histogram size')
    self%l_assigned                  = .true.
    self%name                        = name%to_char()
    self%n_labels                    = size(labels)
    self%labels(1:self%n_labels)     = labels
    self%data(1:self%n_labels)       = data
  end subroutine set

  !---------------- getters ----------------

  ! Return the histogram name, labels, and data; result is .true. if assigned.
  function get( self, name, labels, data ) result( l_assigned )
    class(gui_metadata_histogram),              intent(in)  :: self
    type(string),                               intent(out) :: name
    real,                          allocatable, intent(out) :: labels(:)
    integer,                       allocatable, intent(out) :: data(:)
    logical                                                 :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    l_assigned = self%l_assigned
    name       = trim(self%name)
    allocate(labels(self%n_labels), data(self%n_labels))
    labels = self%labels(1:self%n_labels)
    data   = self%data(1:self%n_labels)
  end function get

  !---------------- serialisation ----------------

  ! Emit name, labels array, and data array as a JSON object.
  ! Returns a null pointer when the object has not been assigned.
  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_histogram), intent(in)    :: self
    type(json_core)                              :: json
    type(json_value),              pointer       :: json_ptr, json_labels_ptr, json_data_ptr
    integer                                      :: i
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( .not.self%l_assigned ) then
      nullify(json_ptr)
      return
    endif
    call json%create_object(json_ptr,       trim(self%name))
    call json%create_array(json_labels_ptr, 'labels')
    call json%create_array(json_data_ptr,   'data')
    do i = 1, self%n_labels
      call json%add(json_labels_ptr, '', dble(self%labels(i)))
      call json%add(json_data_ptr,   '',      self%data(i)  )
    end do
    call json%add(json_ptr, json_labels_ptr)
    call json%add(json_ptr, json_data_ptr)
  end function jsonise_override

  ! Resets every field to its default, so a reused object keeps nothing of an earlier message,
  ! and marks the object uninitialised.
  subroutine kill_override( self )
    class(gui_metadata_histogram), intent(inout) :: self
    select type( self )
      type is( gui_metadata_histogram )
        self = gui_metadata_histogram()
      class default
        call self%gui_metadata_base%kill()
    end select
  end subroutine kill_override

end module simple_gui_metadata_histogram
