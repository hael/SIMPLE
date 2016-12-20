!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 4th of March 2015.
!
! Name:
! simple_mpi_defs - basic definitions for mpi used in all modules.
!
! Description:
! simple_mpi_defs provides basic definitions for the types and declarations
! used in mpi calculations in modules using mpi calls. Using openmpi
!*******************************************************************************
!
module simple_mpi_defs
  !use mpi
  use simple_defs

  implicit none
#if defined (SIMPLE_MPI)
  include 'mpif.h'
#endif
  integer            :: mpierr 
  integer            :: my_rank                 ! rank of process
  integer            :: n_proc                  ! number of processes
  integer            :: len
#if defined (SIMPLE_MPI)
  integer            :: status(MPI_STATUS_SIZE) !status of the communication
  character(MPI_MAX_PROCESSOR_NAME) :: hostname
#endif
  !definitions of general variables
  integer,parameter                 :: MASTER=0
  integer,parameter                 :: FROM_MASTER=1
  integer,parameter                 :: FROM_WORKER=2
contains

  !subroutine to initialise mpi environment
  subroutine simple_mpi_init(mpierr,my_rank,n_proc)
    implicit none
    integer            :: mpierr                  !error status
    integer            :: rc                      !error code return code
    integer            :: my_rank                 ! rank of process
    integer            :: n_proc                  ! number of processes
#ifdef SIMPLE_MPI
    ! start up MPI
    call MPI_Init(mpierr)
    if ( mpierr /= MPI_SUCCESS ) then
       write(*,'(a)')"Error starting MPI program. Terminating."
       call MPI_ABORT(MPI_COMM_WORLD, rc, mpierr)
    end if
    ! find out process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, mpierr)
    ! find out number of processes
    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, mpierr)
#else 
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DSIMPLE_MPI to acces the mpi          "
    write(*,*)"computation using MPI                                           "
    write(*,*)"****************************************************************"
#endif
    return
  end subroutine simple_mpi_init

  !subroutine to finish mpi environment
  subroutine simple_mpi_finalize(mpierr)
    implicit none
    integer            :: mpierr 
#ifdef SIMPLE_MPI
    call MPI_Finalize(mpierr)
#else
    write(*,*)"***************************WARNING******************************"
    write(*,*)"You need to compile with -DSIMPLE_MPI to acces the mpi          "
    write(*,*)"computation using MPI                                           "
    write(*,*)"****************************************************************"
#endif
    return
  end subroutine simple_mpi_finalize

end module simple_mpi_defs
