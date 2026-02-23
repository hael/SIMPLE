!@descr: gui metadata histogram structure
module simple_gui_metadata_histogram
use json_kinds
use json_module
use simple_core_module_api
use simple_gui_metadata_base

implicit none

public :: gui_metadata_histogram
private
#include "simple_local_flags.inc"

type, extends( gui_metadata_base ) :: gui_metadata_histogram
  private
  character(len=SHORTSTRLEN) :: name   = ''
  real                       :: labels(512)
  integer                    :: data(512)
  integer                    :: n_labels

  contains

  procedure :: set
  procedure :: get
  procedure :: jsonise => jsonise_override

end type gui_metadata_histogram

contains

  subroutine set( self, name, labels, data )
    class(gui_metadata_histogram),              intent(inout) :: self
    type(string),                               intent(in)    :: name 
    real,                          allocatable, intent(in)    :: labels(:)
    integer,                       allocatable, intent(in)    :: data(:)
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( size(labels) .ne. size(data) ) THROW_HARD('labels and data differ in size')
    self%l_assigned = .true.
    self%name       = name%to_char()
    self%n_labels   = size(labels)
    self%labels     = transfer(labels, labels)
    self%data       = transfer(data, data)
  end subroutine set

  function get( self, name, labels, data ) result( l_assigned )
    class(gui_metadata_histogram),              intent(inout) :: self
    type(string),                               intent(out)   :: name 
    real,                          allocatable, intent(inout) :: labels(:)
    integer,                       allocatable, intent(inout) :: data(:)
    logical                                                   :: l_assigned
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( allocated(labels) ) deallocate(labels)
    if( allocated(data)   ) deallocate(data)
    l_assigned = self%l_assigned
    name       = trim(self%name)
    allocate(labels(self%n_labels), data(self%n_labels))
    labels = transfer(self%labels, labels, size(labels))
    data   = transfer(self%data,   data,   size(data)  )
  end function get

  function jsonise_override( self ) result( json_ptr )
    class(gui_metadata_histogram), intent(inout) :: self
    type(json_core)                              :: json
    type(json_value),              pointer       :: json_ptr, json_labels_ptr, json_data_ptr
    integer                                      :: i
    if( .not.self%l_initialized ) THROW_HARD('gui metadata object is uninitialised')
    if( self%l_assigned ) then
      call json%create_object(json_ptr, '')
      call json%add(json_ptr, "name",   trim(self%name) )
      call json%create_array(json_labels_ptr, 'labels')
      call json%create_array(json_data_ptr,   'data'  )
      do i=1, self%n_labels
        call json%add(json_labels_ptr, '', dble(self%labels(i)))
        call json%add(json_data_ptr,   '', dble(self%data(i))  )
      enddo
      call json%add(json_ptr, json_labels_ptr)
      call json%add(json_ptr, json_data_ptr  )
    endif
  end function jsonise_override

end module simple_gui_metadata_histogram