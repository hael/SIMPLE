!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 31st of March 2015.
!
! Name:
! initialu - initiales the link variables U in SU(3).
!
! Description:
! subroutine that initialises U to the Identity matrix and attaches it to a 4D
! hypercubic lattice. That is at all point in x,y,z,t and in all directions mu.
!*******************************************************************************
!
      subroutine initialu(ur,ui)
      implicit none
      include 'latticeSize.h'

!     global variable

      integer,parameter                                       :: nc=3
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur,ui

!     local variables

      integer                                                 :: ic

      ur = 0.0d0
      ui = 0.0d0
      do ic = 1, nc
        ur(:,:,:,:,:,ic,ic) = 1.0d0
      end do

      return
      end subroutine initialu
