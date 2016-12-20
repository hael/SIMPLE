!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 21st of March 2015.
!
! Name:
! simple_file_defs - basic definitions for file_utils parrallel IO modules.
!
! Description:
! simple_file_defs provides basic definitions for the types and declarations
! used simple_file_utils for the parralel data IO
!*******************************************************************************
!
module simple_file_defs

  use, intrinsic :: iso_c_binding
  use simple_defs

  implicit none

  integer, parameter :: MAXPROCS = 208
  integer, parameter :: recl_kind=selected_int_kind(16)

  !for reals and complex
  integer, parameter :: f_simple = selected_real_kind(6, 37)
  integer, parameter :: f_double = selected_real_kind(15, 307)
  integer, parameter :: f_quadruple = selected_real_kind(33, 4931)

  !for integers to be verified
  integer, parameter :: f_short=selected_int_kind(4)
  integer, parameter :: f_integer=selected_int_kind(8)
  integer, parameter :: f_long=selected_int_kind(16)

  !data structure for the file_utils and handlers
  type fileDetails_Parallel_IO
    character(len=80) :: file
    integer           :: unit
    character(len=80) :: status
    character(len=80) :: position
    character(len=80) :: action
    logical           :: binary
  end type fileDetails_Parallel_IO

contains

end module simple_file_defs
