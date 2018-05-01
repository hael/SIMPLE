! Parallel Trapezoidal Rule, first version

! Input: None.
! Output: Estimate of the integral from a to b of f(x) using the trapezoidal rule and n trapezoids.

! Algorithm:

! 1.  Each process calculates "its" interval of integration.
! 2.  Each process estimates the integral of f(x) over its interval using the trapezoidal rule.
! 3a. Each process != 0 sends its integral to 0.
! 3b. Process 0 sums the calculations received from the individual processes and prints the result.
program trapezoidal_mpi1
use mpi_f08
implicit none


integer :: my_rank  ! My process rank
integer :: p        ! # processes
real    :: a        ! left endpoint
real    :: b        ! right endpoint
integer :: n        ! # trapezoids
real    :: h        ! Trapezoid base length
real    :: a_local  ! left enpoint for my process
real    :: b_local  ! right endpoint for my process
integer :: n_local  ! # trapezoids for my calculation
real    :: integral ! integral over [a_local,b_local]
real    :: total    ! total integral
integer :: source   ! process sending integral
integer :: dest     ! all messages go to 0
TYPE(MPI_Status) :: status
integer          :: tag, ierr

! data a, b, n, dest, tag /0.0, 1.0, 1024, 0, 50/

! Let the system do what it needs to start up MPI
call MPI_Init(ierr)

! Get rank of my process
call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

! Find out how many processes are being used
call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

call get_data(a, b, n, my_rank, p)

! get and broadcast data
! call get_and_bcast_data(a, b, n, my_rank)

h       = (b - a) / real(n) ! h is the same for all processes
n_local = n / p             ! so is the number of trapeziods

! lenght of each process' interval of integration = n_local * h
! so my interval starts at

a_local  = a + my_rank * n_local * h
b_local  = a_local + n_local * h
integral = Trap(a_local, b_local, n_local, h)

! add up integrals calculated by each process
if( my_rank == 0 )then
    total = integral ! takes care of process 0:s contribution
    do source=1,p-1
        call MPI_Recv(integral, 1, MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierr)
        total = total + integral ! adding the remaining processes contributions
    end do
else
    call MPI_Send(integral, 1, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
endif

! print the result
if( my_rank ==  0 )then
    print *, 'With n = ', n, ' trapezoids, our estimate of the integral from ', a, ' to ', b, ' = ', total
endif

! shut down mpi
call MPI_Finalize(ierr)

contains

    subroutine get_data(a, b, n, my_rank, p)
        real,    intent(out) :: a, b
        integer, intent(out) :: n
        integer, intent(in)  :: my_rank, p
        TYPE(MPI_Status) :: status
        integer :: source, dest, tag, ierr
        if( my_rank == 0 )then
            print *, 'Enter a, b, and n'
            read *, a, b, n
            do dest=1,p-1
                tag = 0
                call MPI_Send(a, 1, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                tag = 1
                call MPI_Send(b, 1, MPI_REAL, dest, tag, MPI_COMM_WORLD, ierr)
                tag = 2
                call MPI_Send(n, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
            end do
        else
            source = 0
            tag = 0
            call MPI_Recv(a, 1, MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierr)
            tag = 1
            call MPI_Recv(b, 1, MPI_REAL, source, tag, MPI_COMM_WORLD, status, ierr)
            tag = 2
            call MPI_Recv(n, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, status, ierr)
        endif
    end subroutine get_data

    subroutine get_and_bcast_data(a, b, n, my_rank)
        real,    intent(out) :: a, b
        integer, intent(out) :: n
        integer, intent(in)  :: my_rank
        integer :: ierr
        if( my_rank == 0 )then
            print *, 'Enter a, b, and n'
            read *, a, b, n
        endif
        ! int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
        call MPI_Bcast(a, 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(b, 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end subroutine get_and_bcast_data

    real function Trap(a_local, b_local, n_local, h)
        real,    intent(in) :: a_local, b_local, h
        integer, intent(in) :: n_local
        integer :: i
        real    :: x, integral
        integral = (f(a_local) + f(b_local)) / 2.0
        x = a_local
        do i=1,n_local - 1
            x = x + h
            integral = integral + f(x)
        end do
        integral = integral * h
        Trap     = integral
    end function Trap

    real function f(x)
        real, intent(in) :: x
        f = x*x
    end function f

end program trapezoidal_mpi1
