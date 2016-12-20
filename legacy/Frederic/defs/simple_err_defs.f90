!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 27th of June 2016.
!
! Name:
! simple_err_defs - basic definitions for error handler in simple_glossary
!
! Description:
! simple_err_defs provides basic definitions for the types and declarations
! used in gpu calculations in modules using cuda calls. Using cuda-5.0
!*******************************************************************************
!
module simple_err_defs
  use, intrinsic :: iso_c_binding
  
  use simple_defs
  use simple_cuda_defs

  implicit none

  !return variable on the success failure state of the error
  !general error return code
  integer :: ERR_SUCCESS, ERR_FAIL
  !IO return code
  integer :: INPUT_OUTPUT_ERROR
  integer :: FILE_OPENING_ERROR
  integer :: ERR_OPEN_FILE_SUCCESS
  integer :: ERR_OPEN_FILE_FAIL
  integer :: ERR_CLOSE_FILE_SUCCESS
  integer :: ERR_CLOSE_FILE_FAIL
  integer :: ERR_INQUIRE_SUCCESS
  integer :: ERR_INQUIRE_FAIL
  !generic errors
  integer :: ERR_GENERIC
  !undefined errors
  integer :: ERR_NOT_DEFINED
  !Dynamic memory  errors
  integer :: ERR_ALLOCATE
  integer :: ERR_DEALLOCATE
  integer :: ERR_MEMLIMIT
  integer :: ERR_INVALID_COPY
  !handlers for the file error module
  character(len=*), parameter :: ERRID='Id'
  character(len=*), parameter :: ERRMSG='Message'
  character(len=*), parameter :: ERRACT='Action'
  character(len=*), parameter :: ERR_ADD_INFO='Additional Info'
  
  !Glossary error values
  integer, save, public :: DICT_VALUE_ABSENT

contains

  !TODO: add methods and handlers here

end module simple_err_defs
