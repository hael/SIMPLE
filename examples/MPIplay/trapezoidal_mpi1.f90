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


integer :: myid     ! My process ID
integer :: ntasks   ! # processes
real    :: a        ! left endpoint
real    :: b        ! right endpoint
integer :: n        ! # trapezoids
real    :: h        ! Trapezoid base length
real    :: a_local  ! left enpoint for my process
real    :: b_local  ! right endpoint for my process
integer :: n_local  ! # trapezoids for my calculation
real    :: integral ! integral over [a_local,b_local]
real    :: total    ! total integral
integer :: dest     ! all messages go to 0
type(MPI_Status)   :: status
type(MPI_Request)  :: request
integer, parameter :: master = 0, msgtag1 = 1, msgtag2 = 2, msgtag3 = 3, msgtag4 = 4
integer            :: ierr, islave

! Let the system do what it needs to start up MPI
call MPI_Init(ierr)

! Get id of my process
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

! Find out how many processes are being used
call MPI_Comm_size(MPI_COMM_WORLD, ntasks, ierr)

! get data
! call get_data(a, b, n, myid, ntasks)

! get and broadcast data
call get_and_bcast_data(a, b, n, myid)

h       = (b - a) / real(n) ! h is the same for all processes
n_local = n / ntasks       

! lenght of each process' interval of integration = n_local * h
! so my interval starts at

a_local  = a + myid * n_local * h
b_local  = a_local + n_local * h
integral = Trap(a_local, b_local, n_local, h)

! add up integrals calculated by each process

! if( myid == 0 )then
!     total = integral ! takes care of process 0:s contribution
!     do islave=1,ntasks - 1
!         call MPI_Recv(integral, 1, MPI_REAL, islave, msgtag4, MPI_COMM_WORLD, status, ierr)
!         total = total + integral ! adding the remaining processes contributions
!     end do
! else
!     call MPI_Isend(integral, 1, MPI_REAL, master, msgtag4, MPI_COMM_WORLD, request, ierr)
! endif
! call MPI_Wait(request,status,ierr)

call MPI_Reduce(integral, total, 1, MPI_REAL, MPI_SUM, master, MPI_COMM_WORLD, ierr)

! print the result
if( myid ==  0 )then
    print *, 'With n = ', n, ' trapezoids, our estimate of the integral from ', a, ' to ', b, ' = ', total
endif

! shut down mpi
call MPI_Finalize(ierr)

contains

    subroutine get_data(a, b, n, myid, p)
        real,    intent(out) :: a, b
        integer, intent(out) :: n
        integer, intent(in)  :: myid, p
        TYPE(MPI_Status) :: status
        integer :: ierr
        if( myid == 0 )then
            print *, 'Enter a, b, and n'
            read *, a, b, n
            do islave=1,ntasks-1
                call MPI_Send(a, 1, MPI_REAL,    islave, msgtag1, MPI_COMM_WORLD, ierr)
                call MPI_Send(b, 1, MPI_REAL,    islave, msgtag2, MPI_COMM_WORLD, ierr)
                call MPI_Send(n, 1, MPI_INTEGER, islave, msgtag3, MPI_COMM_WORLD, ierr)
            end do
        else
            call MPI_Recv(a, 1, MPI_REAL,    master, msgtag1, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(b, 1, MPI_REAL,    master, msgtag2, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(n, 1, MPI_INTEGER, master, msgtag3, MPI_COMM_WORLD, status, ierr)
        endif
    end subroutine get_data

    subroutine get_and_bcast_data(a, b, n, myid)
        real,    intent(out) :: a, b
        integer, intent(out) :: n
        integer, intent(in)  :: myid
        integer :: ierr
        if( myid == 0 )then
            print *, 'Enter a, b, and n'
            read *, a, b, n            
        endif
        ! int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
        call MPI_Bcast(a, 1, MPI_REAL,    master, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(b, 1, MPI_REAL,    master, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(n, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
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
