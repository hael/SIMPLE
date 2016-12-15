! ============================================================================
! Name        : testing_dgetri_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the functionality of the My_zgetri on cuda on GPU
!             :
! ============================================================================
!
program testing_dgetri_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use invert_gpu_utils
  use greeting_version
  use simple_timing

  implicit none
#define devptr_t integer*8
  integer, parameter              :: nmx=13001
  integer, parameter              :: start_nmx=8001, istep=1000

  integer                         :: err
  integer                         :: size_of_elt
  integer                         :: lda, n
  integer                         :: ldb,ldc
  integer                         :: k

  devptr_t                        :: devPtrA_new
  !checking the matrix inversion has work correctly
  devptr_t                        :: devPtrA_invSymZA
  devptr_t                        :: devPtrA_SymZE_gpu
  devptr_t                        :: devPtrA_SymZF_gpu

  real(dp)                        :: alpha,beta

  real(dp),allocatable            :: matSymZE_gpu(:,:)
  real(dp),allocatable            :: invSymZA(:,:)  !inverted sym Z matrix

  real(dp),allocatable            :: matSymZF_gpu(:,:)  !matrix multiplication 

  !timing variables

  integer                         :: inmx

  !index variable
  integer                         :: ilda
  integer                         :: jlda

  !start of the execution commands
  !start of the greeting message
  call hello_gpu_magma()
  call timestamp()
  call start_Alltimers_cpu()

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. 0 ) write(*,*) 'cublas init failed'

  write(*,*)'                                                               '
  write(*,*)'***************************************************************'
  write(*,*)' Now testing the dgetri and dgemm cuBlas and MAGMA_dgetrf (GPU)'
  write(*,*)'***************************************************************'

  do inmx = start_nmx, nmx, istep

     lda = inmx

     allocate(matSymZE_gpu(inmx,inmx))
     allocate(invSymZA(inmx,inmx))

     matSymZE_gpu = 0.0d0

     allocate(matSymZF_gpu(inmx,inmx))

     call random_seed(SIZE=k)       !initializing the seed of the random generator

     call getSymDRandomMat_cpu(inmx,lda,matSymZE_gpu)

     !checking the matrix matSymZE_gpu after initialisation
!     write(*,*)
!     write(*,*) "The product of matrix matSymZE_gpu"
!     write(*,*)
!     do ilda = 1,LDA - (inmx-2)
!        do jlda = 1,LDA - (inmx-2)
!           write(*,'(x,i4,x,i4,4x,f16.8,4x,f16.8)') ilda, jlda, matSymZE_gpu(ilda,jlda)
!        end do
!     end do

    ! now using the magma to get iunverse matrix

     size_of_elt = kind(matSymZE_gpu)
      
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_new)
     if (err .ne. 0 ) call simple_cudblas_stat_return(err)

     err = cublas_set_matrix (inmx, inmx, size_of_elt ,matSymZE_gpu, inmx, devPtrA_new, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     invSymZA = 0.0d0 !initializing the inverted matrix
     call getBlockInverseDMatrix_gpu(inmx,inmx,lda,size_of_elt,inmx,devPtrA_new,invSymZA)
 
     !checking the matrix matSymZE_gpu after CUDA experience
!     write(*,*)
!     write(*,*) "The Real(matSymZE_gpu(ilda,jlda))," , ",aimag(matSymZE_gpu(ilda,jlda))"
!     write(*,*)
!     do ilda = 1,LDA - (inmx-2)
!        do jlda = 1,LDA - (inmx-2)
!           write(*,*) ilda, jlda, matSymZE_gpu(ilda,jlda)
!        end do
!     end do
     !freein memory on the device before doing the matrix multiplication test
     err = cublas_free(devPtrA_new)
     
     !matSymZF_gpu = matmul(matSymZE_gpu,invSymZA)   !this is the CPU slow way
     !this is the cublas fast way for matrix multiplication
     alpha = 1.0d0
     beta = 0.0d0
     ldb = lda
     ldc = lda
     matSymZF_gpu = 0.0d0
     !allocating memory
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_invSymZA)
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_SymZE_gpu)
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_SymZF_gpu)
     !setting up the matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_elt ,invSymZA, inmx, devPtrA_invSymZA, inmx )
     if ( err .ne. 0 ) write(*,*)"Error: ",err," Data upload on the GPU failed"
     err = cublas_set_matrix (inmx, inmx, size_of_elt ,matSymZE_gpu, inmx, devPtrA_SymZE_gpu, inmx )
     if ( err .ne. 0 ) write(*,*)"Error: ",err," Data upload on the GPU failed"
     !performing the matrix multiplication
     call start_timer_cpu("n_cublasDGEMM")
     call cublas_DGEMM("N", "N", inmx, inmx, inmx, alpha, devPtrA_SymZE_gpu, &
          lda, devPtrA_invSymZA, ldb, beta, devPtrA_SymZF_gpu, ldc)
     call stop_timer_cpu("n_cublasDGEMM")
     err = cublas_get_matrix ( inmx, inmx, size_of_elt, devPtrA_SymZF_gpu, inmx, matSymZF_gpu, inmx )

!     write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double Complex precision matrix matSymZF_gpu is: ",sum(matSymZF_gpu)
     write(*,*)
     write(*,'(1x,a,1x,i5,1x,a,f20.8)')"When N =",inmx,&
          "--> Double precision Full sum(A * A^{-1}): ",sum(matSymZF_gpu)
     write(*,*)

     deallocate(matSymZE_gpu)
     deallocate(invSymZA)
     deallocate(matSymZF_gpu)

     err = cublas_free(devPtrA_invSymZA)
     err = cublas_free(devPtrA_SymZE_gpu)
     err = cublas_free(devPtrA_SymZF_gpu)

  end do

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_dgetri_gpu



