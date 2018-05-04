program simple_test_mpi
include 'simple_lib.f08'
use mpi_f08
implicit none
integer ierr

call MPI_INIT ( ierr )
print *, "Hello world"
call MPI_FINALIZE ( ierr )

stop

end program simple_test_mpi
