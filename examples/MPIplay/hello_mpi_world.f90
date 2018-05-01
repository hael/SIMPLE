program hello_mpi_world
use mpi_f08
implicit none
integer :: ierr
call MPI_INIT ( ierr )
print *, "Hello world"
call MPI_FINALIZE ( ierr )
end program hello_mpi_world
