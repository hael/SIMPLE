!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 18th of Jully 2013.
!
! Name:
! matrixGetter - Various utilities and matrix getter for other modules.
!
! Description:
! matrixGetter provides initialisation of matrix to be used in other modules:
!*******************************************************************************
!
module matrixGetter

  use simple_defs
  use simple_lattice_defs

  implicit none

contains
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 1d array of size (n) with random numbers.
!
!*******************************************************************************
! SOURCE

  subroutine getZRandomVec_cpu(n,vecZA)
    implicit none
  
    !global varaibles

    integer                         :: n
    complex(dp)                     :: vecZA(*)

    !local variables

    double precision,dimension(n)   :: harvestr
    double precision,dimension(n)   :: harvesti

    !start of the execution commands

    vecZA(1:n) = 0.0d0         !initializing the matrix

!    call start_timer("getZRandomVec")

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0
    call random_number( harvesti )
    where( harvesti == 0.0d0 ) harvesti = 1.0d0

    vecZA(1:n) = cmplx(harvestr(1:n) , harvesti(1:n) , dp)

!    call stop_timer("getZRandomVec")

    return
  end subroutine getZRandomVec_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
!
!*******************************************************************************
! SOURCE

  subroutine getZRandomMat_cpu(n,m,lda,matZA)
    implicit none
  
    ! global varaibles

    integer                         :: LDA,n,m
    complex(dp)                     :: matZA(LDA,*)

    ! local variables

    double precision,dimension(n*m) :: harvestr
    double precision,dimension(n*m) :: harvesti

    double precision,dimension(n,m) :: re_matZA
    double precision,dimension(n,m) :: im_matZA

    integer                         :: i,j
    integer                         :: size_of_elt

    ! start of the execution commands

    matZA(1:n,1:m) = 0.0d0         !initializing the matrix

!    call start_timer("getZRandomMat")

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0
    call random_number( harvesti )
    where( harvesti == 0.0d0 ) harvesti = 1.0d0

    size_of_elt = 2 * kind(matZA)

    do i=1,n
       do j=1,m
          re_matZA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          im_matZA(i,j) = harvesti( (j - 1)*lda + (i-1) + 1 )
       end do
    end do

    matZA(1:n,1:m) = cmplx(re_matZA(1:n,1:m) , im_matZA(1:n,1:m) , dp)

!    call stop_timer("getZRandomMat")

    return
  end subroutine getZRandomMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
!
!*******************************************************************************
! SOURCE

  subroutine getDRandomMat_cpu(n,m,lda,matDA)
    implicit none
  
    ! global varaibles

    integer                         :: LDA,n,m
    real(dp)                        :: matDA(LDA,*)

    ! local variables

    double precision,dimension(n*m) :: harvestr

    double precision,dimension(n,m) :: re_matDA

    integer                         :: i,j
    integer                         :: size_of_elt

    ! start of the execution commands

    matDA(1:n,1:m) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0

    size_of_elt = 2 * kind(matDA)

    do i=1,n
       do j=1,m
          re_matDA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
       end do
    end do

    matDA(1:n,1:m) = re_matDA(1:n,1:m)

    return
  end subroutine getDRandomMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with fixed numbers.
!
!*******************************************************************************
! SOURCE

  subroutine getDFixedMat_cpu(n,m,lda,matDA,iran)
    implicit none
  
    ! global varaibles

    integer                         :: LDA,n,m,iran
    real(dp)                        :: matDA(LDA,*)

    ! local variables

    double precision,dimension(1)   :: harvestr

    double precision,dimension(n,m) :: re_matDA

    integer                         :: rand_int
    integer                         :: i,j
    integer                         :: size_of_elt

    ! start of the execution commands

    matDA(1:n,1:m) = 0.0d0         !initializing the matrix

    where( harvestr == 0.0d0 ) harvestr = 1.0d0

    rand_int = int(harvestr(1) * 10)

    size_of_elt = 2 * kind(matDA)

    do i=1,n
       do j=1,m
          re_matDA(i,j) = iran * rand_int + (i + j)
       end do
    end do

    matDA(1:n,1:m) = re_matDA(1:n,1:m)

    return
  end subroutine getDFixedMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
!
!*******************************************************************************
! SOURCE

  subroutine getRandomLinks_cpu(ur_rnd,ui_rnd)
    implicit none

    ! global varaibles
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc) :: ur_rnd,ui_rnd

    ! local variables
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc) :: r,theta

    integer,dimension(7)                             :: shapeu=(/nx,ny,nz,nt,mu,nc,nc/)
    double precision,dimension(nx*ny*nz*nt*mu*nc*nc) :: harvestu

    integer                                          :: size_of_elt_ur
    integer                                          :: size_of_elt_ui
    
    ! start of the execution commands

    !initializing the matrix
    call initialu(ur_rnd,ui_rnd)

    !generating the random number!
    call random_number( harvestu )
    where( harvestu==0.0d0 ) harvestu = 1.0d0
    r = reshape( harvestu,shapeu )
    call random_number( harvestu )
    where( harvestu==0.0d0 ) harvestu = 1.0d0
    theta = reshape( harvestu,shapeu )

    ur_rnd = r
    ui_rnd = theta

    size_of_elt_ur = 2 * kind(ur_rnd)
    size_of_elt_ui = 2 * kind(ui_rnd)

    return
  end subroutine getRandomLinks_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
!
!*******************************************************************************
! SOURCE
  subroutine getSU3RandomLinks_cpu(ur_rnd,ui_rnd)
    implicit none

    ! global varaibles
    double precision,dimension(nx,ny,nz,nt,mu,nc,nc)        :: ur_rnd,ui_rnd

    ! local variables

    integer                         :: size_of_elt_ur
    integer                         :: size_of_elt_ui
    
    ! start of the execution commands

    !initializing the matrix
    call initialu(ur_rnd,ui_rnd)

    size_of_elt_ur = 2 * kind(ur_rnd)
    size_of_elt_ui = 2 * kind(ui_rnd)
    
    !getting some uniformally generated SU(3) matrix
    call su3random(ur_rnd,ui_rnd)

    return
  end subroutine getSU3RandomLinks_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! and construct symmetric matrix A such that A^{T} = A
!
!*******************************************************************************
! SOURCE

  subroutine getSymZRandomMat_cpu(n,lda,matSymZA)
    implicit none
  
    ! global varaibles

    integer                                    :: LDA,n
    complex(dp)                                :: matSymZA(LDA,*)

    ! local variables

    double precision,dimension(n*n)           :: harvestr
    double precision,dimension(n*n)           :: harvesti

    double precision,dimension(n,n)           :: re_matSymZA
    double precision,dimension(n,n)           :: im_matSymZA

    integer                                   :: i,j

    ! start of the execution commands

    matSymZA(1:n,1:n) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0
    call random_number( harvesti )
    where( harvesti == 0.0d0 ) harvesti = 1.0d0

    re_matSymZA(1:n,1:n) = 0.0d0 
    im_matSymZA(1:n,1:n) = 0.0d0 

    do i=1,n
       do j=i,n
          re_matSymZA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          im_matSymZA(i,j) = harvesti( (j - 1)*lda + (i-1) + 1 )
          re_matSymZA(j,i) = re_matSymZA(i,j)
          im_matSymZA(j,i) = im_matSymZA(i,j)
       end do
    end do

    matSymZA(1:n,1:n) = cmplx(re_matSymZA(1:n,1:n) , im_matSymZA(1:n,1:n) , dp)

    return
  end subroutine getSymZRandomMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a double precision 2d array of size (n)x(m) with random numbers.
! and construct symmetric matrix A such that A^{T} = A
!
!*******************************************************************************
! SOURCE

  subroutine getSymDRandomMat_cpu(n,lda,matSymDA)
    implicit none
  
    ! global varaibles

    integer                                    :: LDA,n
    double precision                           :: matSymDA(LDA,*)

    ! local variables

    double precision,dimension(n*n)           :: harvestr
    double precision,dimension(n,n)           :: re_matSymDA

    integer                                   :: i,j

    ! start of the execution commands

    matSymDA(1:n,1:n) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0

    re_matSymDA(1:n,1:n) = 0.0d0 

    do i=1,n
       do j=i,n
          re_matSymDA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          re_matSymDA(j,i) = re_matSymDA(i,j)
       end do
    end do

    matSymDA(1:n,1:n) = re_matSymDA(1:n,1:n)

    return
  end subroutine getSymDRandomMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a double precision 2d array of size (n)x(m) with random numbers.
! and construct symmetric matrix A such that A^{T} = A
!
!*******************************************************************************
! SOURCE

  subroutine getSymSRandomMat_cpu(n,lda,matSymSA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n
    real(sp)                                  :: matSymSA(LDA,*)

    ! local variables

    real(sp),dimension(n*n)                   :: harvestr
    real(sp),dimension(n,n)                   :: re_matSymSA

    integer                                   :: i,j

    ! start of the execution commands

    matSymSA(1:n,1:n) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0

    re_matSymSA(1:n,1:n) = 0.0d0 

    do i=1,n
       do j=i,n
          re_matSymSA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          re_matSymSA(j,i) = re_matSymSA(i,j)
       end do
    end do

    matSymSA(1:n,1:n) = re_matSymSA(1:n,1:n)

    return
  end subroutine getSymSRandomMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a double precision 2d array of size (n)x(m) with with a
! sequence number 1 to 9 in row major order (fortran) style.
!
!*******************************************************************************
! SOURCE
  subroutine get1to9RowMajZdpMat_cpu(n,lda,matZA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n
    complex(dp)                               :: matZA(LDA,*)

    !local variables 

    real(dp),allocatable :: re_matZA(:,:)
    real(dp),allocatable :: im_matZA(:,:)

    ! counters

    integer :: i,j

    !start of the excution commands

    allocate(re_matZA(n,n))
    allocate(im_matZA(n,n))

    !  matZA = 0.0d0
    re_matZA = 0.0d0
    im_matZA = 0.0d0

    call getIndentityMatrix_cpu(n,LDA,matZA)

    do i=1,n
       do j=1,n
          re_matZA(i,j) = (i - 1)*lda + (j-1) + 1 !row major
          im_matZA(i,j) = (i - 1)*lda + (j-1) + 1 !row major
       end do
    end do

    matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , dp)

    deallocate(re_matZA)
    deallocate(im_matZA)

    return
  end subroutine get1to9RowMajZdpMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a real(dp) 2d array of size (n)x(m) with with a
! sequence number 1 to 9 in coluumn major order (C) style.
!
!*******************************************************************************
! SOURCE
  subroutine get1to9ColMajZdpMat_cpu(n,lda,matZA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n
    complex(dp)                               :: matZA(LDA,*)

    !local variables 

    real(dp),allocatable :: re_matZA(:,:)
    real(dp),allocatable :: im_matZA(:,:)

    ! counters

    integer :: i,j

    !start of the excution commands

    allocate(re_matZA(n,n))
    allocate(im_matZA(n,n))

    !  matZA = 0.0d0
    re_matZA = 0.0d0
    im_matZA = 0.0d0

    call getIndentityMatrix_cpu(n,LDA,matZA)

    do i=1,n
       do j=1,n
          re_matZA(i,j) = (j - 1)*lda + (i-1) + 1 !column major
          im_matZA(i,j) = (j - 1)*lda + (i-1) + 1 !column major
       end do
    end do

    matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , dp)

    deallocate(re_matZA)
    deallocate(im_matZA)

    return
  end subroutine get1to9ColMajZdpMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a double precision 2d array of size (n)x(m) with with a
! sequence number 1 to 9 in row major order (fortran) style.
!
!*******************************************************************************
! SOURCE
  subroutine get1to9RowMajZspMat_cpu(n,lda,matZA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n
    complex(sp)                               :: matZA(LDA,*)

    !local variables 

    real(sp),allocatable :: re_matZA(:,:)
    real(sp),allocatable :: im_matZA(:,:)

    ! counters

    integer :: i,j

    !start of the excution commands

    allocate(re_matZA(n,n))
    allocate(im_matZA(n,n))

    !  matZA = 0.0d0
    re_matZA = 0.0
    im_matZA = 0.0

    call getIndentityZspMatrix_cpu(n,LDA,matZA)

    do i=1,n
       do j=1,n
          re_matZA(i,j) = (i - 1)*lda + (j-1) + 1 !row major
          im_matZA(i,j) = (i - 1)*lda + (j-1) + 1 !row major
       end do
    end do

    matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , sp)

    deallocate(re_matZA)
    deallocate(im_matZA)

    return
  end subroutine get1to9RowMajZspMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a real(sp) 2d array of size (n)x(m) with with a
! sequence number 1 to 9 in coluumn major order (C) style.
!
!*******************************************************************************
! SOURCE
  subroutine get1to9ColMajZspMat_cpu(n,lda,matZA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n
    complex(sp)                               :: matZA(LDA,*)

    !local variables 

    real(sp),allocatable :: re_matZA(:,:)
    real(sp),allocatable :: im_matZA(:,:)

    ! counters

    integer :: i,j

    !start of the excution commands

    allocate(re_matZA(n,n))
    allocate(im_matZA(n,n))

    !  matZA = 0.0d0
    re_matZA = 0.0
    im_matZA = 0.0

    call getIndentityZspMatrix_cpu(n,LDA,matZA)

    do i=1,n
       do j=1,n
          re_matZA(i,j) = (j - 1)*lda + (i-1) + 1 !column major
          im_matZA(i,j) = (j - 1)*lda + (i-1) + 1 !column major
       end do
    end do

    matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , sp)

    deallocate(re_matZA)
    deallocate(im_matZA)

    return
  end subroutine get1to9ColMajZspMat_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a real(sp) 2d array of size (n)x(m) with with a
! sequence number 1 to 9 in coluumn major order (C) style.
!
!*******************************************************************************
! SOURCE
  subroutine get1to9ColMajDspMat_3D_cpu(n,m,l,lda,matSA)
    implicit none
  
    ! global varaibles

    integer                                   :: LDA,n,m,l
    real(sp)                                  :: matSA(n,m,*)

    !local variables 

    real(sp),allocatable :: re_matSA(:,:,:)

    ! counters

    integer :: i,j,k

    !start of the excution commands

    allocate(re_matSA(n,m,l))

    !  matZA = 0.0d0
    re_matSA = 0.0

    call getIndentityDspMatrix_3D_cpu(n,m,l,LDA,matSA)

    do i=1,n
       do j=1,m
          do k=1,l
             re_matSA(i,j,k) = ((j - 1) + m*(k-1))*n + (i-1) + 1 !column major
          end do
       end do
    end do

    matSA(1:n,1:m,1:l) = re_matSA(1:n,1:m,1:l)

    deallocate(re_matSA)

    return
  end subroutine get1to9ColMajDspMat_3D_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! using the random generator from the GPU magma routines.
!*******************************************************************************
! SOURCE
subroutine getIndentityMatrix_cpu(n,LDA,matZA)
  implicit none

  !global variables

  integer :: LDA,n
  complex(dp) :: matZA(LDA,*)

  !local variables 

  double precision,allocatable :: re_matZA(:,:)
  double precision,allocatable :: im_matZA(:,:)

  integer :: in

  !start of the excution commands

  allocate(re_matZA(n,n))
  allocate(im_matZA(n,n))

!  matZA = 0.0d0
  re_matZA = 0.0d0
  im_matZA = 0.0d0

  forall ( in=1:n ) re_matZA(in,in) = 1.0d0

  matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , dp)

  deallocate(re_matZA)
  deallocate(im_matZA)

  return
end subroutine getIndentityMatrix_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! using the random generator from the GPU magma routines.
!*******************************************************************************
! SOURCE
subroutine getIndentityZspMatrix_cpu(n,LDA,matZA)
  implicit none

  !global variables

  integer :: LDA,n
  complex(sp) :: matZA(LDA,*)

  !local variables 

  real(sp),allocatable :: re_matZA(:,:)
  real(sp),allocatable :: im_matZA(:,:)

  integer :: in

  !start of the excution commands

  allocate(re_matZA(n,n))
  allocate(im_matZA(n,n))

!  matZA = 0.0d0
  re_matZA = 0.0
  im_matZA = 0.0

  forall ( in=1:n ) re_matZA(in,in) = 1.0

  matZA(1:n,1:n) = cmplx(re_matZA(1:n,1:n) , im_matZA(1:n,1:n) , sp)

  deallocate(re_matZA)
  deallocate(im_matZA)

  return
end subroutine getIndentityZspMatrix_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! using the random generator from the GPU magma routines.
!*******************************************************************************
! SOURCE
subroutine getIndentityDspMatrix_3D_cpu(n,m,l,LDA,matDA)
  implicit none

  !global variables

  integer :: LDA,n,m,l
  real(sp) :: matDA(n,m,*)

  !local variables

  real(sp),allocatable :: re_matDA(:,:,:)

  integer :: in

  !start of the excution commands

  allocate(re_matDA(n,m,l))

!  matZA = 0.0d0
  re_matDA = 0.0

  forall ( in=1:n ) re_matDA(in,in,in) = 1.0

  matDA(1:n,1:m,1:l) = re_matDA(1:n,1:m,1:l)

  deallocate(re_matDA)

  return
end subroutine getIndentityDspMatrix_3D_cpu

end module matrixGetter
