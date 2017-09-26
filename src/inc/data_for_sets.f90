
type VECTOR_DATA
    type(SET_DATA) :: data
 end type VECTOR_DATA
public :: SET_DATA
public :: element_isequal, operator(.eq.)

type(VECTOR_DATA), save :: empty_vector_data ! Do not initialise it -
                                             ! it is just a place holder

include 'simple_vector.f90'
