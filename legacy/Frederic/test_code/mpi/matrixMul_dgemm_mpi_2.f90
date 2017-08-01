! ============================================================================
! Name        : matrixMultiplier_mpi.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 22th of October 2014
! Copyright   : Your copyright notice
! Description : Calculate matrix multiplication of nxm mxn dense matrix in MPI
! ============================================================================
program matrixMul_dgemm_mpi_2

  use simple_defs
  use greeting_version
  use simple_mpi_defs
  use matrixGetter
  use simple_random
  use simple_timing
  use simple_matrixHandler_cpu

  implicit none
  include 'matrixMultiplier_mpi.h'
  !local variables
  !The counters
  integer                                    :: im,in,ik

  !start of the execution commands
  call timestamp()
  call start_Alltimers_cpu()
  call init_fixed_random_seed(iseed)
  !lets get some random matrix
  k=n
  allocate(matDA(m,k))
  allocate(matDB(k,n))
  allocate(matDC(m,n))
  allocate(matDA_mpi(m,k))
!  allocate(matDB_mpi(k,n))
  allocate(matDC_mpi(m,n))
  lda = k
  ldb = n
  ldc = m

  call getDRandomMat_cpu(m,n,lda,matDA)
  call getDRandomMat_cpu(m,n,ldb,matDB)

  call start_timer_cpu("my_dgemm_cpu")
  call my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
  call stop_timer_cpu("my_dgemm_cpu")

  write(*,*)'Printing the result from my_dgemm_cpu'
  do im = 1,3
     do in = 1,3
        print *, "matDC(", im, ",", in, ") = ", matDC(im,in)
     end do
  end do

  call hello_mpi(my_rank,n_proc)
  call simple_mpi_init(mpierr,my_rank,n_proc)
  call MPI_GET_PROCESSOR_NAME(hostname, len, mpierr)
  write(*,*)'Number of proc = ',n_proc,' My rank = ',my_rank, &
       "Length of host = ",len,'Running on = ',hostname

  write(*,*)MASTER,FROM_MASTER,FROM_WORKER
  write(*,*)mpierr

  call start_timer_cpu("my_dgemm_mpi_2")
  call my_dgemm_mpi_2(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC_mpi,ldc)
  call stop_timer_cpu("my_dgemm_mpi_2")

  !end  of the MPI computation
  call simple_mpi_finalize(mpierr)
  write(*,*)'Printing the result'

  do im = 1,2
     do in = 1,2
        print *, "matDC_mpi(", im, ",", in, ") = ", matDC_mpi(im,in),"diff = ",matDC_mpi(im,in)-matDC(im,in)
     end do
  end do

  !freeing the memory from the allocated matrix A,B and C
  deallocate(matDA)
  deallocate(matDB)
  deallocate(matDC)
  deallocate(matDA_mpi)
!  deallocate(matDB_mpi)
  deallocate(matDC_mpi)

  call stop_Alltimers_cpu()
  call timestamp()

end program matrixMul_dgemm_mpi_2
