! ============================================================================
! Name        : matrixMultiplier_mpi.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 22th of October 2014
! Copyright   : Your copyright notice
! Description : Calculate matrix multiplication of nxm mxn dense matrix in MPI
! ============================================================================
program matrixMul_dgemm_mpi_1

  use simple_defs
  use greeting_version
  use simple_mpi_defs
  use matrixGetter
  use simple_random
  use simple_timing
  use simple_matrixHandler_cpu

  implicit none
  include 'matrixMultiplier_mpi.h'
  !The counters
  integer                                    :: im,in

  !start of the execution commands
!  call spacer(1,20,.true.,"lines")

  call timestamp()
  call start_Alltimers_cpu()
  call init_fixed_random_seed(iseed)

  !lets get some random matrix
  k=n
  allocate(matDA(m,k))
  allocate(matDB(k,n))
  allocate(matDC(m,n))
  allocate(matDC_mpi(m,n))
  lda = k
  ldb = n
  ldc = m

  call getDRandomMat_cpu(n,m,lda,matDA)
  call getDRandomMat_cpu(n,m,ldb,matDB)

  call start_timer_cpu("my_dgemm_cpu")
  call my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
  call stop_timer_cpu("my_dgemm_cpu")

  !start of the MPI computation
  call hello_mpi(my_rank,n_proc)
  call simple_mpi_init(mpierr,my_rank,n_proc)
  call MPI_GET_PROCESSOR_NAME(hostname, len, mpierr)
  write(*,*)'Number of proc = ',n_proc,' My rank = ',my_rank, &
       "Length of host = ",len,'Running on = ',hostname

  call start_timer_cpu("my_dgemm_mpi_1")
  call my_dgemm_mpi_1(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC_mpi,ldc)
  call stop_timer_cpu("my_dgemm_mpi_1")

  !end  of the MPI computation
  call simple_mpi_finalize(mpierr)

  !Check if the matrices are the same by taking teh difference
  
  write(*,'(1x,a)')"The difference of the first 5 elements of (matDC-matDC_mpi) = "
  do im=1, 5
     do in = 1, 5
        write(*,'(2x,f16.8,$)')matDC(im,in) - matDC_mpi(im,in)
     end do
     print *, ' '
  end do
  write(*,'(1x,a,2x,f16.8)')"The sum of the difference sum((matDC-matDC_mpi)) = ",&
       sum(matDC(1:m,1:n) - matDC_mpi(1:m,1:n))

!  call bye_mpi()
  !freeing the memory from the allocated matrix A,B and C
  deallocate(matDA)
  deallocate(matDB)
  deallocate(matDC)
  deallocate(matDC_mpi)

  call stop_Alltimers_cpu()
  call timestamp()
!  call spacer(1,20,.true.,"lines")

end program matrixMul_dgemm_mpi_1
