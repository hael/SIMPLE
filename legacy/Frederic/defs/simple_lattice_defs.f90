!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 18th of Jully 2013.
!
! Name:
! physsym_lattice_defs - basic definitions used in all modules.
!
! Description:
! physsym_lattice_defs provides basic definitions for the types and declarations
! and all the lattice QCD modules used through out the code. It is based on
! numerical nrtypes.f90 module.
!*******************************************************************************
!
module simple_lattice_defs

  implicit none

  !inlcude the lattice size
  include 'latticeSize.h'

  !color and directional parameters
  integer,parameter                                   :: nc=3 !The SU(3) color
  integer,parameter                                   :: mu=4 !the direction
  !Dirac indices
  integer,parameter                                   :: nd=4

end module simple_lattice_defs
