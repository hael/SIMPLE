! ============================================================================
! Name        : testing_dgetri_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the functionality of the My_zgetri on CPU using the Lacpacks
!             :
! ============================================================================
!
program testing_dgetri_cpu

  use simple_defs
  use matrixGetter
  use invert_cpu_utils
  use greeting_version
  use simple_timing

  implicit none
  integer, parameter              :: nmx=3001
  integer, parameter              :: start_nmx=1001, istep=1000

  integer                         :: err
  integer                         :: size_of_elt
  integer                         :: lda, n
  integer                         :: ldb,ldc
  integer                         :: k

  real(dp)                        :: alpha,beta

  real(dp),allocatable            :: matSymDE_cpu(:,:)
  real(dp),allocatable            :: invSymDA(:,:)  !inverted sym Z matrix

  !checking the matrix inversion has work correctly
  real(dp),allocatable            :: matSymDF_cpu(:,:)  !matrix multiplication 

  !timing variables

  integer                         :: inmx

  !index variable
  integer                         :: ilda
  integer                         :: jlda

  !start of the execution commands
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  !starting the cuda environment
  write(*,*)'                                                               '
  write(*,*)'***************************************************************'
  write(*,*)'     Now testing the dgetri and dgemm Lapacks (CPU)            '
  write(*,*)'***************************************************************'

  do inmx = start_nmx, nmx, istep
     lda = inmx

     allocate(matSymDE_cpu(inmx,inmx))
     allocate(invSymDA(inmx,inmx))
     allocate(matSymDF_cpu(inmx,inmx))

     matSymDE_cpu = 0.0d0
     invSymDA = 0.0d0
     matSymDF_cpu = 0.0d0
     
     call random_seed(SIZE=k)   !initializing the seed of the random generator

     call getSymDRandomMat_cpu(inmx,lda,matSymDE_cpu)

     !checking the matrix matSymDE_cpu after initialisation
!     write(*,*)
!     write(*,*) "The product of matrix matSymDE_cpu"
!     write(*,*)
!     do ilda = 1,LDA - (inmx-2)
!        do jlda = 1,LDA - (inmx-2)
!           write(*,'(x,i4,x,i4,4x,f16.8,4x,f16.8)') ilda, jlda, matSymDE_cpu(ilda,jlda)
!        end do
!     end do

    ! now using the magma to get iunverse matrix

     size_of_elt = kind(matSymDE_cpu)

     invSymDA = 0.0d0 !initializing the inverted matrix
     call getBlockInverseDMatrix_cpu(inmx,inmx,lda,size_of_elt,inmx,matSymDE_cpu,invSymDA)
 
     !checking the matrix matSymDE_cpu after CUDA experience
!     write(*,*)
!     write(*,*) "The matSymDE_cpu(ilda,jlda),"
!     write(*,*)
!     do ilda = 1,LDA - (inmx-2)
!        do jlda = 1,LDA - (inmx-2)
!           write(*,*) ilda, jlda, matSymDE_cpu(ilda,jlda)
!        end do
!     end do
     
     !matSymDF_cpu = matmul(matSymDE_cpu,invSymDA)   !this is the CPU slow way
     !this is the cublas fast way for matrix multiplication
     alpha = 1.0d0
     beta = 0.0d0
     ldb = lda
     ldc = lda
     matSymDF_cpu = 0.0d0

     !performing the matrix multiplication
     call start_timer_cpu("l_dgemm")
     call dgemm("N", "N", inmx, inmx, inmx, alpha, matSymDE_cpu, &
          lda, invSymDA, ldb, beta, matSymDF_cpu, ldc)
     call stop_timer_cpu("l_dgemm")

     write(*,'(1x,a,1x,i5,1x,a,f20.8)')"When N =",inmx,&
          "--> Double precision Full sum(A * A^{-1}): ",sum(matSymDF_cpu)

     deallocate(matSymDE_cpu)
     deallocate(invSymDA)
     deallocate(matSymDF_cpu)

  end do

  !shutting down the timers
  call stop_Alltimers_cpu()

end program testing_dgetri_cpu



