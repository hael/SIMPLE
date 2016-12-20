!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 22nd of March 2015.
!
! Name:
! simple_resources_defs - basic definitions for available resources available
!                         on systems to describe system.
!
! Description:
! simple_resources_defs provides basic definitions for the available resources
! available on the system in addition of systtem query pseudo class.
! This is to be used in the setting of the Parrallel IO and MPI.
! sets the number of nodes avvailable to the user or the number of nodes the
! user would like ot use. The size of the socket boith logical and physical.
! 
!*******************************************************************************
!
module simple_resources_defs
  use, intrinsic :: iso_c_binding
  
  use simple_defs

  implicit none
  !data structure for the file_utils and handlers
  type,bind(c) :: resources_avail
     integer(c_int) :: nnodes           !N of nodes wanted to be used
     integer(c_int) :: size_socket_logi !N of logical cores on each sockets 
     integer(c_int) :: size_socket_phys !N physical cores on each socket
  end type resources_avail

contains

end module simple_resources_defs
