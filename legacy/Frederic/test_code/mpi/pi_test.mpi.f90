! ============================================================================
! Name        : mpi_pi_example.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 19th of July 2013
! Copyright   : Your copyright notice
! Description : Calculate Pi in MPI
! ============================================================================

subroutine calc_pi(rank, num_procs)
  use mpi
  
  implicit none

  integer, intent(in) :: rank
  integer, intent(in) :: num_procs

  integer          :: i
  integer          :: ierror
  integer          :: num_intervals
  double precision :: h
  double precision :: mypi
  double precision :: pi
  double precision :: sum
  double precision :: x

  ! set number of intervals to calculate
  if (rank == 0) num_intervals = 1000000000

  ! tell other tasks how many intervals
  call MPI_Bcast(num_intervals, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

  ! now everyone does their calculation

  h = 1.0d0 / num_intervals
  sum = 0.0d0

  do i = rank + 1, num_intervals, num_procs
     x = h * (i - 0.5d0);
     sum = sum + (4.0d0 / (1.0d0 + x*x))
  end do

  mypi = h * sum

  ! combine everyone's calculations
  call MPI_Reduce(mypi, pi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
       MPI_COMM_WORLD, ierror)

  if (rank == 0) print *, "PI is approximately ", pi
end subroutine calc_pi

program mpi_pi_example
  use mpi
  use simple_timing

  implicit none

  integer, parameter :: LEN = 100               ! message length

  integer            :: ierror                  ! error code
  integer            :: my_rank                 ! rank of process
  integer            :: num_procs               ! number of processes
  integer            :: source                  ! rank of sender
  integer            :: dest                    ! rank of receiver
  integer            :: tag                     ! tag for messages
  character(len=LEN) :: message                 ! storage for message
  integer            :: status(MPI_STATUS_SIZE) ! return status for receive

  dest = 0
  tag = 0
  call start_Alltimers_cpu()

  ! start up MPI
  call MPI_Init(ierror)

  ! find out process rank
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)

  ! find out number of processes
  call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierror)


  if (my_rank .ne. 0) then
     ! create message
     write (message, *) "Greetings from process ", my_rank
     call MPI_Send(message, LEN, MPI_CHARACTER, &
          dest, tag, MPI_COMM_WORLD, ierror)
  else
     print *, "Num processes: ", num_procs
     do source = 1, num_procs-1
        call MPI_Recv(message, LEN, MPI_CHARACTER, source, tag, &
             MPI_COMM_WORLD, status, ierror)
        print *, "Process 0 received ", message
     end do

     ! now return the compliment
     write (message, *) "Hi, how are you?"
  end if

  call MPI_Bcast(message, LEN, MPI_CHARACTER, dest, MPI_COMM_WORLD, ierror)

  if (my_rank .ne. 0) then
     print *, "Process ", my_rank, " received ", message
  end if

  ! calculate PI
  call start_timer_cpu("calc_pi")
  call calc_pi(my_rank, num_procs)
  call stop_timer_cpu("calc_pi")

  ! shut down MPI
  call MPI_Finalize(ierror)

  call stop_Alltimers_cpu()

  stop
end program mpi_pi_example
