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
    use simple_defs
    use simple_error, only: allocchk
!    include 'simple_lib.f08'
    implicit none
    !> Generic interface to 1D, 2D real, complex, and C-complex fftshift routines
    interface fftshift
        module procedure fftshift2cc,  fftshift3cc, fftshift3c, fftshift3r, fftshift2c, &
            &fftshift2r, fftshift1c, fftshift1r
    end interface fftshift
    !> Generic interface to 1D, 2D real, complex, and C-complex inverse fftshift routines
    interface ifftshift
        module procedure  ifftshift3cc, ifftshift3c, ifftshift3r, ifftshift2cc, &
            &ifftshift2c, ifftshift2r, ifftshift1c, ifftshift1r
    end interface ifftshift
    private
#include "simple_local_flags.inc"
    public :: fftshift, ifftshift
contains


    !! 1D fftshifters
    ! complex
    subroutine fftshift1c(fftimage)
        complex, intent(inout) :: fftimage(:)
        complex, allocatable   :: tmp(:)
        integer                :: nx, lbx, i, ii
        nx=size(fftimage); lbx = lbound(fftimage,1)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 1D complex tmp failed")
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(ii) = fftimage(i-lbx+1)
        enddo
        fftimage(lbx:nx-lbx+1) = tmp(:nx)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift1c

    ! real 1D fftshift
    subroutine fftshift1r(fftimage)
        real, intent(inout) :: fftimage(:)
        integer :: nx, lbx
        real,allocatable    :: tmp(:)
        integer :: i,ii
        nx=size(fftimage); lbx = lbound(fftimage,1)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 1D real tmp failed")
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx
            tmp(ii) = fftimage(i-lbx+1)
        enddo
        fftimage(lbx:nx-lbx+1) = tmp(:nx)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift1r

    !! 2D fftshifters

    ! single precision real 2D
    subroutine fftshift2r(fftimage)
        real, intent(inout) :: fftimage(:,:)
        integer::nx,ny,lbx,lby
        real, allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 2D real '
        nx=size(fftimage,1); lbx = lbound(fftimage,1)
        ny=size(fftimage,2); lby = lbound(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 2D real tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i-lbx+1,j-lby+1)
            enddo
        enddo
        fftimage(lbx:nx-lbx+1, lby:ny-lby+1) = tmp(:nx,:ny)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift2r

    ! single precision complex 2D
    subroutine fftshift2c(fftimage)
        implicit none
        complex, intent(inout) :: fftimage(:,:)
        integer::nx,ny,lbx,lby
        complex, allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 2D complex '
        nx=size(fftimage,1); lbx = lbound(fftimage,1)
        ny=size(fftimage,2); lby = lbound(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 2D complex tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i-lbx+1,j-lby+1)
            enddo
        enddo
        fftimage(lbx:nx-lbx+1, lby:ny-lby+1) = tmp(:nx,:ny)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift2c

    ! double precision complex 2D
    subroutine fftshift2cc(fftimage)
        complex(kind=dp), intent(inout) :: fftimage(:,:)
        integer::nx,ny,lbx,lby
        complex(kind=dp), allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 2D complex double '
        nx=size(fftimage,1); lbx = lbound(fftimage,1)
        ny=size(fftimage,2); lby = lbound(fftimage,2)
        allocate(tmp(nx,ny),stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 2D double complex tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i-lbx+1,j-lby+1)
            enddo
        enddo
        fftimage(lbx:nx-lbx+1, lby:ny-lby+1) = tmp(:nx,:ny)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift2cc

    ! ! double precision complex - pseudo 3D slices
    subroutine fftshift2cc_z(fixed_axis,fftimage)
        integer, intent(in) :: fixed_axis
        complex(kind=dp), intent(inout) :: fftimage(:,:,:)
        integer::nx,ny,lbx,lby
        complex(kind=dp), allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 2D slice complex'
        nx=size(fftimage,1); lbx = lbound(fftimage,1)
        ny=size(fftimage,2); lby = lbound(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 2D slice complex tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(ii,jj) = fftimage(i-lbx+1,j-lby+1,fixed_axis)
            enddo
        enddo
        fftimage(lbx:nx-lbx+1, lby:ny-lby+1,fixed_axis) = tmp(:nx,:ny)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift2cc_z

    !! 3D FFTSHIFTERS

    ! 3D single precision float
    subroutine fftshift3r(fftimage)
        real, intent(inout) :: fftimage(:,:,:)
        integer :: nx,ny,nz,lbx,lby,lbz
        real, allocatable  :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 3D real'
        nx=size(fftimage,1); lbx = lbound(fftimage,1)
        ny=size(fftimage,2); lby = lbound(fftimage,2)
        nz=size(fftimage,3); lbz = lbound(fftimage,3)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 3D real tmp failed")
        ! omp parallel do collapse(3) default(shared) private(i,j,k,ii,jj,kk)&
        ! omp schedule(static) proc_bind(close)
        do k=1,nz
            do j=1, ny
                do i=1, nx
                    jj = mod(j+(ny+1)/2,ny)
                    if(jj==0) jj = ny
                    kk = mod(k+(nz+1)/2,nz)
                    if(kk==0) kk = nz
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i-lbx+1,j-lby+1,k-lbz+1)
                enddo
            enddo
        enddo
        ! omp end parallel do
        fftimage(lbx:nx-lbx+1,lby:ny-lby+1,lbz:nz-lbz+1) = tmp(:nx,:ny,:nz)
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift3r

    !3D single precision complex
    subroutine fftshift3c(fftimage,lims)
        implicit none
        complex, intent(inout)  :: fftimage(:,:,:)
        integer, intent(inout),optional :: lims(3,2)
        integer :: nx,ny,nz
        complex, allocatable  :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk,lms(3,2)
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 3D complex'
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        VerbosePrint 'In fftshifter::fftshift 3D complex ', nx, ny, nz
        if(.not.present(lims))then
            lms(1,1) = lbound(fftimage,1)
            lms(1,2) = ubound(fftimage,1)
            lms(2,1) = lbound(fftimage,2)
            lms(2,2) = ubound(fftimage,2)
            lms(3,1) = lbound(fftimage,3)
            lms(3,2) = ubound(fftimage,3)
        else
            lms=lims
        endif
        VerbosePrint "fftshift3c CMAT lower/upper bounds ", lms
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 3D complex tmp failed")
        ! omp parallel do default(shared) private(i,j,k,ii,jj,kk)&
        ! omp schedule(static) proc_bind(close) collapse(1)
        do k=1,nz
            do j=1, ny
                do i=1, nx
                    kk = mod(k+(nz+1)/2,nz)
                    if(kk==0) kk = nz

                    jj = mod(j+(ny+1)/2,ny)
                    if(jj==0) jj = ny

                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i-lms(1,1)+1,j-lms(2,1)+1,k-lms(3,1)+1)
                enddo
            enddo
        enddo
        ! omp end parallel do
        fftimage(lms(1,1):lms(1,2),lms(2,1):lms(2,2),lms(3,1):lms(3,2)) = tmp
        if(allocated(tmp))deallocate(tmp)
        VerbosePrint 'In fftshifter::fftshift 3D complex done'
    end subroutine fftshift3c

    ! 3D double precision complex
    subroutine fftshift3cc(fftimage,lims)
        complex(kind=dp), intent(inout) :: fftimage(:,:,:)
        integer, intent(in) :: lims(3,2)
        integer :: nx,ny,nz,lbx,lby,lbz
        complex(kind=dp), allocatable  :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk
        verbose=.true.
        VerbosePrint 'In fftshifter::fftshift 3D complex double'
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::fftshift 3D double complex tmp failed")
        ! omp parallel do collapse(3) default(shared) private(i,j,k,ii,jj,kk)&
        ! omp schedule(static) proc_bind(close)
        do k=1, nz
            do j=1, ny
                do i=1, nx
                    kk = mod(k+(nz+1)/2,nz)
                    if(kk==0) kk = nz
                    jj = mod(j+(ny+1)/2,ny)
                    if(jj==0) jj = ny
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(ii,jj,kk) = fftimage(i-lims(1,1)+1,j-lims(2,1)+1,k-lims(3,1)+1)
                enddo
            enddo
        enddo
        ! omp end parallel do
        fftimage(lims(1,1):lims(1,2),lims(2,1):lims(2,2),lims(3,1):lims(3,2)) = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine fftshift3cc

    !! INVERSE IFFTSHIFTERS

    ! 1D real
    subroutine ifftshift1r(fftimage)
        real, intent(inout) :: fftimage(:)
        integer::nx
        real,allocatable  :: tmp(:)
        integer::i,ii

        nx=size(fftimage)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 1D real tmp failed")
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift1r

    ! 1D complex
    subroutine ifftshift1c(fftimage)
        implicit none
        complex, intent(inout) :: fftimage(:)
        integer::nx
        complex,allocatable  :: tmp(:)
        integer::i,ii

        nx=size(fftimage)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 1D complex tmp failed")
        do i=1, nx
            ii = mod(i+(nx+1)/2,nx)
            if(ii==0) ii = nx

            tmp(i) = fftimage(ii)
        enddo
        fftimage = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift1c

    ! 2D real single-precision
    subroutine ifftshift2r(fftimage)
        real, intent(inout) :: fftimage(:,:)
        integer::nx,ny
        real, allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint "in ifftshift2r "
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 2D real tmp failed")
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
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift2r

    ! 2D single precision complex
    subroutine ifftshift2c(fftimage)
        complex, intent(inout) :: fftimage(:,:)
        integer::nx,ny,lbx,lby
        complex, allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint "in ifftshift2c "
        nx=size(fftimage,1);lbx=lbound(fftimage,1)
        ny=size(fftimage,2);lby=lbound(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 2D complex tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx

                tmp(i,j) = fftimage(ii-lbx+1,jj-lby+1)
            end do
        end do
        fftimage(lbx:nx-lbx+1,lby:ny-lby+1) = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift2c

    ! 2D double complex
    subroutine ifftshift2cc(fftimage)
        complex(kind=dp), intent(inout) :: fftimage(:,:)
        integer :: nx,ny,lbx,lby
        complex(kind=dp), allocatable  :: tmp(:,:)
        integer::i,j,ii,jj
        verbose=.true.
        VerbosePrint "in ifftshift2cc"
        nx=size(fftimage,1);lbx=lbound(fftimage,1)
        ny=size(fftimage,2);lby=lbound(fftimage,2)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 2D double complex tmp failed")
        do j=1, ny
            jj = mod(j+(ny+1)/2,ny)
            if(jj==0) jj = ny
            do i=1, nx
                ii = mod(i+(nx+1)/2,nx)
                if(ii==0) ii = nx
                tmp(i,j) = fftimage(ii-lbx+1,jj-lby+1)
            end do
        end do
        fftimage(lbx:nx-lbx+1,lby:ny-lby+1) = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift2cc

    ! ! 2D (3D slices) double complex
    ! subroutine ifftshift2cc_z(fftimage, fixed_axis)
    !     implicit none
    !     complex(kind=dp), intent(inout) :: fftimage(:,:,:)
    !     integer, intent(in) :: fixed_axis
    !     integer::nx,ny,lbx,lby
    !     complex(kind=dp), allocatable  :: tmp(:,:)
    !     integer::i,j,ii,jj
    !     verbose=.true.
    !     VerbosePrint "in ifftshift2cc_z "
    !     nx=size(fftimage,1);lbx=lbound(fftimage,1)
    !     ny=size(fftimage,2);lby=lbound(fftimage,2)
    !     allocate(tmp(nx,ny),source=fftimage(:,:,1),stat=alloc_stat)
    !     if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 2D slice complex tmp failed")
    !     do j=1, ny
    !         jj = mod(j+(ny+1)/2,ny)
    !         if(jj==0) jj = ny
    !         do i=1, nx
    !             ii = mod(i+(nx+1)/2,nx)
    !             if(ii==0) ii = nx

    !             tmp(i,j) = fftimage(ii-lbx+1,jj-lby+1,fixed_axis)
    !         end do
    !     end do
    !     fftimage(lbx:nx-lbx+1,lby:ny-lby+1,fixed_axis) = tmp
    !     deallocate(tmp)
    ! end subroutine ifftshift2cc_z

    !! 3D INVERSE FFTSHIFTERS
    ! 3D double precision complex
    subroutine ifftshift3cc(fftimage,lims)
        complex(kind=dp), intent(inout) :: fftimage(:,:,:)
        integer, intent(in) :: lims(3,2)
        integer :: nx,ny,nz
        complex(kind=dp), allocatable  :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk
        verbose=.true.
        VerbosePrint "in ifftshift3cc "
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 3D double complex tmp failed")
        do k=1, nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz

            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                !  if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(i,j,k) = fftimage(ii-lims(1,1)+1,jj-lims(2,1)+1,kk-lims(3,1)+1)
                enddo
            enddo
        enddo
        ! omp end parallel do
        fftimage(lims(1,1):lims(1,2),lims(2,1):lims(2,2),lims(3,1):lims(3,2)) = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift3cc

    !3D single precision complex
    subroutine ifftshift3c(fftimage,lims)
        complex, intent(inout)  :: fftimage(:,:,:)
        integer, intent(inout),optional :: lims(3,2)
        integer :: nx,ny,nz,lims_here(3,2)
        complex, allocatable  :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk
        verbose=.true.
        VerbosePrint 'In fftshifter::ifftshift 3D complex'
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        VerbosePrint 'In fftshifter::ifftshift 3D complex size', nx, ny, nz
        if(.not.present(lims))then
            lims_here(1,1) = lbound(fftimage,1)
            lims_here(1,2) = ubound(fftimage,1)
            lims_here(2,1) = lbound(fftimage,2)
            lims_here(2,2) = ubound(fftimage,2)
            lims_here(3,1) = lbound(fftimage,3)
            lims_here(3,2) = ubound(fftimage,3)
        else
            lims_here = lims
        endif
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 3D complex tmp failed")
        do k=1,nz
            kk = mod(k+(nz+1)/2,nz)
            if(kk==0) kk = nz
            do j=1, ny
                jj = mod(j+(ny+1)/2,ny)
                if(jj==0) jj = ny
                do i=1, nx
                    ii = mod(i+(nx+1)/2,nx)
                    if(ii==0) ii = nx
                    tmp(i,j,k) = fftimage(ii-lims_here(1,1)+1,jj-lims_here(2,1)+1,kk-lims_here(3,1)+1)
                enddo
            enddo
        enddo
        ! omp end parallel do
        fftimage(lims_here(1,1):lims_here(1,2),lims_here(2,1):lims_here(2,2),lims_here(3,1):lims_here(3,2)) = tmp
        if(allocated(tmp))deallocate(tmp)
    end subroutine ifftshift3c

    ! 3D single precision float
    subroutine ifftshift3r(fftimage)
        real, intent(inout) :: fftimage(:,:,:)
        integer :: nx,ny,nz
        real, allocatable   :: tmp(:,:,:)
        integer :: i,j,k,ii,jj,kk
        verbose=.true.
        VerbosePrint 'In fftshifter::ifftshift 3D real'
        nx=size(fftimage,1)
        ny=size(fftimage,2)
        nz=size(fftimage,3)
        if(allocated(tmp))deallocate(tmp)
        allocate(tmp(nx,ny,nz),source=fftimage,stat=alloc_stat)
        if(alloc_stat /= 0) call allocchk("simple_fftshifter::ifftshift 3D real tmp failed")
        do k=1,nz
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
        if(allocated(tmp))deallocate(tmp)
        VerbosePrint 'In fftshifter::ifftshift 3D real done'
    end subroutine ifftshift3r

end module simple_fftshifter
