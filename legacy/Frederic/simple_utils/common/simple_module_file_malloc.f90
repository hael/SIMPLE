! ============================================================================
! Name        : simple_module_file_malloc
! Author      : Frederic Bonnet
! Version     : 1.0
! Date        : 10th of July 2016
! Description : Module to define a the file malloc
! ============================================================================
!
module simple_module_file_malloc
  use simple_eglossary

  !global parameter of the module telling if the profile has to be activated
  !this parameter can be modified only by dynamic memory module
  integer, parameter :: file_malloc_namelen=32 !length of the character
  logical, save, public :: file_malloc_default_profiling=.true.
  character(len=file_malloc_namelen), save, public :: file_malloc_routine_name=repeat(' ',file_malloc_namelen)

  !TODO: insert the public methods here for handling the file memory
  
end module simple_module_file_malloc


  
