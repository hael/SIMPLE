program matrixMul_dgemm_mpi_3

  use simple_defs
  use greeting_version
  use simple_mpi_defs
  use matrixGetter
  use simple_random
  use simple_timing
  use simple_matrixHandler_cpu
!  use physsym_matrixHandler_cpu

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
!  allocate(matDA_mpi(m,k))
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
  do im = 1,2
     do in = 1,2
        print *, "matDC(", im, ",", in, ") = ", matDC(im,in)
     end do
  end do

  call simple_mpi_init(mpierr,my_rank,n_proc)
  call hello_mpi(my_rank,n_proc)
  call MPI_GET_PROCESSOR_NAME(hostname, len, mpierr)
  write(*,*)'Number of proc = ',n_proc,' My rank = ',my_rank, &
       "Length of host = ",len,'Running on = ',hostname

  call start_timer_cpu("my_dgemm_mpi_3")
  call my_dgemm_mpi_3(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC_mpi,ldc)
  call stop_timer_cpu("my_dgemm_mpi_3")

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
!  deallocate(matDA_mpi)
!  deallocate(matDB_mpi)
  deallocate(matDC_mpi)

  call stop_Alltimers_cpu()
  call timestamp()

  end program matrixMul_dgemm_mpi_3
!*******************************************************************************
! DESCRIPTION
! subroutine to multiply dense mxn matrix using mpi distributed memory
!*******************************************************************************
! SOURCE
  subroutine my_dgemm_mpi_3(my_rank,n_proc,mpierr,m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
    implicit none 
    include 'mpif.h'
    !global variables
    integer                                    :: my_rank      !rank of process
    integer                                    :: n_proc       !number of processes
    integer                                    :: mpierr 
    integer                                    :: lda,ldb,ldc
    integer                                    :: n,m,k
    double precision,intent(in)                :: matDA(lda,*)
    double precision,intent(in)                :: matDB(ldb,*)
    double precision,intent(inout)             :: matDC(ldc,*)
    !local variables
    integer                                    :: nworkers
    !MPI variables
    integer                                    :: status(MPI_STATUS_SIZE)
    integer,parameter                          :: MASTER = 0
    integer,parameter                          :: FROM_MASTER = 1
    integer,parameter                          :: FROM_WORKER = 2
    !The counters
    integer                                    :: im,in,ik
    integer                                    :: i_proc
    integer                                    :: iworkers, isource
    !start of the execution commands

    !initializing the matDC
    matDC(1:m,1:n) = 0.0d0
#ifdef SIMPLE_MPI
    nworkers = n_proc - 1
    !*************************  Master task **************************
    if ( my_rank .eq. MASTER ) then
       !send matrix data to the work tasks
       do iworkers=1,nworkers
          call MPI_SEND(matDA, m*k, MPI_DOUBLE_PRECISION, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
          call MPI_SEND(matDB, k*n, MPI_DOUBLE_PRECISION, iworkers,FROM_MASTER, MPI_COMM_WORLD, mpierr)
       end do
       !getting matDC back from the worker
       do isource = 1, nworkers
          call MPI_RECV(matDC, m*n, MPI_DOUBLE_PRECISION, isource, FROM_WORKER, MPI_COMM_WORLD, status, mpierr)
       end do

    endif
  !*************************  Worker task **************************
    if ( my_rank > MASTER ) then
       !receiving the matDA and matDB from the master
       call MPI_RECV(matDA, m*k, MPI_DOUBLE_PRECISION, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)
       call MPI_RECV(matDB, k*n, MPI_DOUBLE_PRECISION, MASTER, FROM_MASTER, MPI_COMM_WORLD, status, mpierr)
       !performing the product
       do im=1,m
          do in=1,n
             do ik=1,k
                matDC(im,in) = matDC(im,in) + matDA(im,ik) * matDB(ik,in)
             end do
          end do
       end do
       !send result back to master task
       call MPI_SEND(matDC, m*n, MPI_DOUBLE_PRECISION, MASTER, FROM_WORKER, MPI_COMM_WORLD, mpierr)

    end if

#else 
    write(*,*)"***************************WARNING************************************"
    write(*,*)"You need to compile with -DSIMPLE_MPI to acces the mpi computtation  "
    write(*,*)"Proceeding with standard 1-CPU calculation using my_degemm_cpu routine"
    write(*,*)"**********************************************************************"
!    call my_dgemm_cpu(m,n,k,matDA,lda,matDB,ldb,matDC,ldc)
#endif

  return
end subroutine my_dgemm_mpi_3
