!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 31st of March 2015.
!
! Name:
! su2random - subroutine that generates a ramdom SU(2) matrix.
!
! Description:
! subroutine that generates a ramdom SU(2) matrix and attaches it to a 4D
! hypercubic lattice. That is at all point in x,y,z,t and in all directions mu.
!*******************************************************************************
!
      subroutine su2random(ursu2,uisu2)
      implicit none
      include 'latticeSize.h'

!     global variables

      integer,parameter                                       :: ncsu2=2
      integer,parameter                                       :: nsigma=ncsu2*ncsu2
      integer,parameter                                       :: mu=4

      double precision,dimension(nx,ny,nz,nt,mu,nsigma)       :: a4vector
      double precision,dimension(nx,ny,nz,nt,mu,ncsu2,ncsu2)  :: ursu2,uisu2

!     local variables

      integer,parameter                                       :: nr=2
      double precision,parameter                              :: pi=3.141592653d0
      double precision,dimension(nx,ny,nz,nt,mu)              :: norm
      double precision,dimension(nx,ny,nz,nt,mu,nr)           :: r,theta
      double precision,dimension(nx*ny*nz*nt*mu*nr)           :: harvestr
      integer,dimension(6)                                    :: shaper=(/nx,ny,nz,nt,mu,nr/)
      integer                                                 ::isigma

!     call one set of random numbers and put them in a big array called r,
!     it can be splited into 2, accessing both sides of the array using parameter nr.
!     the random number is in (0,1) taking the nat log maps it onto (-infty,0]
!     multiplying it by -2 maps it to [0,infty)

      call random_number( harvestr )
      where( harvestr==0.0d0 ) harvestr = 1.0d0
      r = reshape( harvestr,shaper )
      r = sqrt( -2.0d0*log(r) )
      call random_number( harvestr )
      theta = reshape( harvestr,shaper )
      theta = 2.0d0 * pi * theta

!     a4vetor 1 and 2 are independent and gaussian distributed with
!     mean 0 and standard deviation 1. repeat process to get a4vec 3 and 4.
!     It is the cos and sine that normally distributed.
!     a4vector becomes normaly distributed on S^4

      a4vector(:,:,:,:,:,1) = r(:,:,:,:,:,1) * cos( theta(:,:,:,:,:,1) )
      a4vector(:,:,:,:,:,2) = r(:,:,:,:,:,1) * sin( theta(:,:,:,:,:,1) )
      a4vector(:,:,:,:,:,3) = r(:,:,:,:,:,2) * cos( theta(:,:,:,:,:,2) )
      a4vector(:,:,:,:,:,4) = r(:,:,:,:,:,2) * sin( theta(:,:,:,:,:,2) )

      norm = sqrt( sum( a4vector**2,dim=6 ) )

!     it when we divide by the norm= r1**2+r2**2 that the points are brought
!     back onto the surface of the sphere

      do isigma = 1, nsigma
         a4vector(:,:,:,:,:,isigma) = a4vector(:,:,:,:,:,isigma) / norm(:,:,:,:,:)
      end do

!     converting a4vector to an su2 matrix

      ursu2(:,:,:,:,:,1,1) = a4vector(:,:,:,:,:,4)
      ursu2(:,:,:,:,:,2,2) = a4vector(:,:,:,:,:,4)
      ursu2(:,:,:,:,:,1,2) = a4vector(:,:,:,:,:,2)
      ursu2(:,:,:,:,:,2,1) =-a4vector(:,:,:,:,:,2)
      uisu2(:,:,:,:,:,1,1) = a4vector(:,:,:,:,:,3)
      uisu2(:,:,:,:,:,2,2) =-a4vector(:,:,:,:,:,3)
      uisu2(:,:,:,:,:,1,2) = a4vector(:,:,:,:,:,1)
      uisu2(:,:,:,:,:,2,1) = a4vector(:,:,:,:,:,1)

      return
      end subroutine su2random
