! ============================================================================
! Name        : testing_zgetri_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the functionality of the My_zgetri on cuda
!             :
! ============================================================================
!
program testing_zgetri_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use invert_gpu_utils
  use greeting_version
  use simple_timing

  implicit none
#define devptr_t integer*8
  integer, parameter              :: nmx=11500
  integer, parameter              :: start_nmx=5500, istep=1000
!  integer, parameter              :: start_nmx=9500, istep=5
  integer, parameter              :: n_by_m_chkA = 3 !size check matrix 4 vec decomp

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
  !device pointer for the vector to check the matrix inversion
  devptr_t                        :: devPtrA_vecZEr
  devptr_t                        :: devPtrA_vecZEi
  devptr_t                        :: devPtrA_vecZAr
  devptr_t                        :: devPtrA_vecZAi

  complex(dp)                     :: alpha,beta

  real(dp),allocatable            :: vecZEr(:), vecZEi(:) !vectors of matrix
  real(dp),allocatable            :: vecZAr(:), vecZAi(:) !vectors of matrix

  complex(dp),allocatable         :: result_mat(:,:)

  complex(dp),allocatable         :: matSymZE_gpu(:,:)
  complex(dp),allocatable         :: invSymZA(:,:)  !inverted sym Z matrix

  complex(dp),allocatable         :: matSymZF_gpu(:,:)  !matrix multiplication 

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
  write(*,*)' Now testing the zgetri and zgemm cuBlas and MAGMA_zgetrf (GPU)'
  write(*,*)'***************************************************************'

  do inmx = start_nmx, nmx, istep

     lda = inmx

     allocate(matSymZE_gpu(inmx,inmx))
     allocate(invSymZA(inmx,inmx))

     matSymZE_gpu = 0.0d0

     allocate(matSymZF_gpu(inmx,inmx))

     call random_seed(SIZE=k)       !initializing the seed of the random generator

     call getSymZRandomMat_cpu(inmx,lda,matSymZE_gpu)

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

     size_of_elt = 2 * kind(matSymZE_gpu)
      
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_new)
     if (err .ne. 0 ) call simple_cudblas_stat_return(err)

     err = cublas_set_matrix (inmx, inmx, size_of_elt ,matSymZE_gpu, inmx, devPtrA_new, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     invSymZA = 0.0d0 !initializing the inverted matrix
     call getBlockInverseZMatrix_gpu(inmx,inmx,lda,size_of_elt,inmx,devPtrA_new,invSymZA)
 
     !checking the matrix matSymZE_gpu after CUDA experience
!     write(*,*)
!     write(*,*) 'The original matrix after the GPU inversion' 
!     write(*,*) "The Real(matSymZE_gpu(ilda,jlda))," , ",aimag(matSymZE_gpu(ilda,jlda))"
!     write(*,*)
!     do ilda = 1,LDA - (inmx-2)
!        do jlda = 1,LDA - (inmx-2)
!           write(*,*) ilda, jlda, Real(matSymZE_gpu(ilda,jlda))," , ",aimag(matSymZE_gpu(ilda,jlda))
!        end do
!     end do
     !freein memory on the device before doing the matrix multiplication test
     err = simple_cublas_free(devPtrA_new)

     if (inmx < 9501 ) then
        
        write(*,*) 'The size of the matrix is less than: 9501'
        write(*,*) 'Performing a cublas_zgemm for checking inversion'

        !matSymZF_gpu = matmul(matSymZE_gpu,invSymZA)   !this is the CPU slow way
        !this is the cublas fast way for matrix multiplication
        alpha = (1.0d0,0.0d0)
        beta = (0.0d0,0.0d0)
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
        call start_timer_cpu("n_cublasZgemm")
        call cublas_ZGEMM("N", "N", inmx, inmx, inmx, alpha, devPtrA_SymZE_gpu, &
             lda, devPtrA_invSymZA, ldb, beta, devPtrA_SymZF_gpu, ldc)
        err = cublas_get_matrix ( inmx, inmx, size_of_elt, devPtrA_SymZF_gpu, inmx, matSymZF_gpu, inmx )
        call stop_timer_cpu("n_cublasZgemm")

!        write(*,*)
!        write(*,*) "The Real(matSymZF_gpu(ilda,jlda))," , ",aimag(matSymZF_gpu(ilda,jlda))"
!        write(*,*)
!        do ilda = 1,LDA - (inmx-2)
!           do jlda = 1,LDA - (inmx-2)
!              write(*,*) ilda, jlda, Real(matSymZF_gpu(ilda,jlda))," , ",aimag(matSymZF_gpu(ilda,jlda))
!           end do
!        end do
        write(*,*)
        write(*,'(1x,a,1x,i5,1x,a,2f20.8)')"When N =",inmx,&
             "--> Double precision Full sum(A * A^{-1}): ",sum(matSymZF_gpu)
        write(*,*)

        err = simple_cublas_free(devPtrA_invSymZA)
        err = simple_cublas_free(devPtrA_SymZE_gpu)
        err = simple_cublas_free(devPtrA_SymZF_gpu)

     else if ( inmx > 9501 ) then

        write(*,*) 'The size of the matrix exceeds n: ',inmx
        write(*,*) 'Performing a vector decomposition for checking inversion'

        !disecting matrix invSymZA and matSymZE_gpu into the real and imaginary part
        !into vectors

        allocate(result_mat(n_by_m_chkA,n_by_m_chkA))
        result_mat = 0.0d0

        !allocating the ressources for the computation
        allocate(vecZEr(inmx))
        allocate(vecZEi(inmx))
        allocate(vecZAr(inmx))
        allocate(vecZAi(inmx))

        call start_timer_cpu("vector_chk")
        do ilda=1,n_by_m_chkA
           do jlda=1,n_by_m_chkA
              vecZEr(1:inmx) = Real(matSymZE_gpu(ilda,1:inmx))
              vecZEi(1:inmx) = aimag(matSymZE_gpu(ilda,1:inmx))

              vecZAr(1:inmx) = Real(invSymZA(1:inmx,jlda))
              vecZAi(1:inmx) = aimag(invSymZA(1:inmx,jlda))

              result_mat(ilda,jlda) = cmplx( &
                   sum(vecZEr(1:inmx) * vecZAr(1:inmx) - vecZEi(1:inmx) * vecZAi(1:inmx)) , &
                   sum(vecZEr(1:inmx) * vecZAi(1:inmx) + vecZEi(1:inmx) * vecZAr(1:inmx)) , &
                   dp)
           end do
        end do
        call stop_timer_cpu("vector_chk")

!        write(*,*)
!        write(*,*) "The Real(result_mat(ilda,jlda))," , ",aimag(result_mat(ilda,jlda))"
!        write(*,*)
!        do ilda = 1,3
!           do jlda = 1,3
!              write(*,*) ilda, jlda, result_mat(ilda,jlda)
!           end do
!        end do
        write(*,*)
        write(*,'(1x,a,1x,i5,1x,a,2f20.8)')"When N =",inmx,&
             "--> Double precision 3x3 subblok sum(A * A^{-1}): ",sum(result_mat)
        write(*,*)
        !dealocating the ressources on host
        deallocate(result_mat)
        deallocate(vecZEr)
        deallocate(vecZEi)
        deallocate(vecZAr)
        deallocate(vecZAi)

     end if

     deallocate(matSymZE_gpu)
     deallocate(invSymZA)
     deallocate(matSymZF_gpu)

  end do

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_zgetri_gpu



