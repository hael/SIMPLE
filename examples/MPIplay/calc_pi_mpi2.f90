program calc_pi_mpi1
use mpi_f08
implicit none

integer                     :: i, n, myid, ntasks, ierr, islave
TYPE(MPI_Status)            :: status
integer,          parameter :: master = 0, msgtag1 = 11, msgtag2 = 12 
double precision, parameter :: pi25dt = 3.141592653589793238462643d0
double precision            :: a,h,pi,sum,x,mypi
character(len=32)           :: cmd_arg

! initialization of the MPI env
call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, myid,   ierr)
call MPI_Comm_size(MPI_COMM_WORLD, ntasks, ierr)

! broadcasting of input data (n)
if( myid == 0 )then ! only the master
    call get_command_argument(1, cmd_arg)
    read(cmd_arg,'(I3)') n
endif
call MPI_Bcast(n, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

! parallel calculation of the quadrature (summation)
h = 1.d0 / dble(n) ! stride
sum = 0.d0
!$omp parallel do default(shared) private(i,x) reduction(+:sum) schedule(static) proc_bind(close)
do i=myid + 1, n, ntasks
    x = h * (dble(i) - 0.5d0)
    sum = sum + f(x)
end do
!$omp end parallel do
mypi = h * sum

! reduction of the subtotals
call MPI_Reduce(mypi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, master, MPI_COMM_WORLD, ierr)
if( myid == 0 ) print *, 'pi is approximately: ', pi, ' Error is: ', abs(pi-pi25dt)

! finalize
call MPI_Finalize(ierr)

contains

    double precision function f( a )
        double precision, intent(in) :: a
        f = 4.d0 / (1.d0 + a*a)
    end function f

end program calc_pi_mpi1
