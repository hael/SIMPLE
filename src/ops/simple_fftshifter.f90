!>------------------------------------------------------------
!! Supplies 1D and 2D FFT shift procedures ala matlab fftshift
!!
!! Uses a generic interface so that procedures can be called with
!! any variation of complex, real, 1D or 2D arrays
!! 2D can also be called with C_DOUBLE_COMPLEX variables
!!
!! fftshift swaps the left and right sides of an array.  For a
!!  2D array, the top and bottom halves are also swapped.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! Modified Oct 2017 - 3D
!! Michael Eager (michael dot eager at monash dot edu)
!!------------------------------------------------------------
module simple_fftshifter
    include 'simple_lib.f08'
    implicit none
    !> Generic interface to 1D, 2D real, complex, and C-complex fftshift routines
    interface fftshift
        module procedure fftshift2cc,  fftshift3cc, fftshift3c, fftshift3r, fftshift2c, fftshift2r, fftshift1c, fftshift1r
    end interface
    !> Generic interface to 1D, 2D real, complex, and C-complex inverse fftshift routines
    interface ifftshift
        module procedure ifftshift2cc_z,  ifftshift3cc, ifftshift3c, ifftshift3r, ifftshift2cc, ifftshift2c, ifftshift2r, ifftshift1c, ifftshift1r
    end interface
contains
    !! 1D fftshifters
    ! complex
    subroutine fftshift1c(fftimage)
        implicit none
        complex, intent(inout),dimension(:):: fftimage
        integer::nx
        complex,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(ii) = fftimage(i)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift1c
    ! real
    subroutine fftshift1r(fftimage)
        implicit none
        real, intent(inout),dimension(:):: fftimage
        integer::nx
        real,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(ii) = fftimage(i)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift1r

    !! 2D fftshifters
    ! single precision real
    subroutine fftshift2r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        real, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2r

    ! single precision complex
    subroutine fftshift2c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2c

    ! double precision complex
    subroutine fftshift2cc(fftimage)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex(kind=dp), allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j)
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift2cc

    ! double precision complex - pseudo 3D slices
    subroutine fftshift2cc_z(fftimage,fixed_axis)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:,:):: fftimage
        integer, intent(in) :: fixed_axis
        integer::nx,ny
        complex(kind=dp), allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i,j,fixed_axis)
            enddo
        enddo
        fftimage(:,:,fixed_axis) = tmp
        deallocate(tmp)
    end subroutine fftshift2cc_z
    !! 3D FFTSHIFTERS

    ! 3D single precision float
    subroutine fftshift3r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:,:):: fftimage
        integer :: nx,ny,nz
        real, allocatable,dimension(:,:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny,nz))
        do k=1,nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i,j,k)
                enddo
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift3r

    !3D single precision complex
    subroutine fftshift3c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:,:) :: fftimage
        integer :: nx,ny,nz
        complex, allocatable,dimension(:,:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny,nz))
        do k=1,nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i,j,k)
                enddo
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift3c

    ! 3D double precision complex
     subroutine fftshift3cc(fftimage)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:,:):: fftimage
        integer :: nx,ny,nz
        complex(kind=dp), allocatable,dimension(:,:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny,nz))
        do k=1, nz
             kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i,j,k)
                enddo
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine fftshift3cc

    !! INVERSE IFFTSHIFTERS

    ! 1D real
    subroutine ifftshift1r(fftimage)
        implicit none
        real, intent(inout),dimension(:):: fftimage
        integer::nx
        real,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift1r

    ! 1D complex
    subroutine ifftshift1c(fftimage)
        implicit none
        complex, intent(inout),dimension(:):: fftimage
        integer::nx
        complex,allocatable,dimension(:) :: tmp
        integer::i,ii

        nx=size(fftimage)
        allocate(tmp(nx))
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift1c

     ! 2D real single-precision
    subroutine ifftshift2r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        real, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj)
            end do
        end do
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2r

    ! 2D single precision complex
    subroutine ifftshift2c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex, allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj)
            end do
        end do
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2c

    ! 2D double complex
    subroutine ifftshift2cc(fftimage)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:):: fftimage
        integer::nx,ny
        complex(kind=dp), allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx
                tmp(i,j) = fftimage(ii,jj)
            end do
        end do
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift2cc

    ! 2D (3D slices) double complex
    subroutine ifftshift2cc_z(fftimage, fixed_axis)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:,:):: fftimage
        integer, intent(in) :: fixed_axis
        integer::nx,ny
        complex(kind=dp), allocatable,dimension(:,:) :: tmp
        integer::i,j,ii,jj

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        allocate(tmp(nx,ny))
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii,jj,fixed_axis)
            end do
        end do
        fftimage(:,:,fixed_axis) = tmp
        deallocate(tmp)
    end subroutine ifftshift2cc_z

    !! 3D INVERSE FFTSHIFTERS
    ! 3D double precision complex
     subroutine ifftshift3cc(fftimage)
        implicit none
        complex(kind=dp), intent(inout),dimension(:,:,:):: fftimage
        integer :: nx,ny,nz
        complex(kind=dp), allocatable,dimension(:,:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny,nz))
        do k=1, nz
             kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz

            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(i,j,k) = fftimage(ii,jj,kk)
                enddo
            enddo
        enddo
        fftimage = tmp
        deallocate(tmp)
    end subroutine ifftshift3cc

    !3D single precision complex
    subroutine ifftshift3c(fftimage)
        implicit none
        complex, intent(inout),dimension(:,:,:) :: fftimage
        integer :: nx,ny,nz
        complex, allocatable,dimension(:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny))
        do k=1,nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(i,j) = fftimage(ii,jj,kk)
                enddo
            enddo
            fftimage(:,:,kk) = tmp
        enddo
        deallocate(tmp)
    end subroutine ifftshift3c

    ! 3D single precision float
    subroutine ifftshift3r(fftimage)
        implicit none
        real, intent(inout),dimension(:,:,:):: fftimage
        integer :: nx,ny,nz
        real, allocatable,dimension(:,:) :: tmp
        integer :: i,j,k,ii,jj,kk

        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        allocate(tmp(nx,ny))
        do k=1,nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(i,j) = fftimage(ii,jj,kk)
                enddo
            enddo
            fftimage(:,:,kk) = tmp
        enddo
        deallocate(tmp)
    end subroutine ifftshift3r

end module simple_fftshifter
