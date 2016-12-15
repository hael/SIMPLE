! primes_mpi.f90
!     Determine the first 1000 primes - using MPI
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program primes_mpi
    use mpi

    implicit none

    integer                  :: rank
    integer                  :: main               = 0
    integer, parameter       :: tag_new_task       = 1 ! Get new task
    integer, parameter       :: tag_results        = 2 ! Transmit results
    integer                  :: error

    integer                  :: number_images

    integer, dimension(2)    :: range
    integer                  :: number_tasks
    integer, dimension(1000) :: prime
    integer                  :: number_primes
    integer, dimension(MPI_STATUS_SIZE) :: status
    logical                  :: new_results
    logical                  :: new_task
    logical                  :: ready

    call mpi_init( error )
    if ( error /= MPI_SUCCESS ) then
        write(*,*) 'Error while starting MPI-based program:' ,error
        stop
    endif

    call mpi_comm_rank( MPI_COMM_WORLD, rank, error )
    call mpi_comm_size( MPI_COMM_WORLD, number_images, error )

    if ( number_images < 2 ) then
        write( *, *) 'We need at least two images to run!'
        call mpi_finalize
        stop
    endif

    !
    ! What we do depends on the rank:
    ! Rank 0 is the main program that hands out the chunks
    ! and gathers the results
    ! All others do the work
    !

    if ( rank == 0 ) then
        ready         = .false.
        new_results   = .false.
        number_tasks  = 0
        number_primes = 0

        ! The main program:
        ! Hand out the chunks and receive the results
        !
        range(2)      = 0

        do while ( .not. ready )

            call handle_communication
        enddo
    else
        !
        ! Worker programs:
        ! Get a task (new range, determine the primes)
        !
        new_task = .false.
        do
            call get_range

            if ( new_task ) then
                call find_primes
                call mpi_send( prime, number_primes, MPI_INTEGER, main, &
                         tag_results, MPI_COMM_WORLD, status, error )
            else
                exit
            endif
        enddo
    endif

    !
    ! Print the results
    !
    if ( rank == 0 ) then
        write(*,'(10i5)') prime
    endif

    call mpi_finalize

    stop

contains

!
! Communicate with the worker images
!
subroutine handle_communication
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer                             :: error
    integer                             :: count
    integer, dimension(100)             :: result
    integer                             :: i
    integer                             :: number_store

    do
        call mpi_recv( result, size(result), MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
                 MPI_COMM_WORLD, status, error )

        if ( error /= MPI_SUCCESS ) then
            write(*,*) 'Some error occurred while receiving: ',error
            call mpi_finalize( error )
            stop
        endif

        !
        ! Send a new task or store the results?
        !
        if ( status(MPI_TAG) == tag_new_task ) then
            range(1) = range(2) + 1
            range(2) = range(2) + 100
            call mpi_send( range, 2, MPI_INTEGER, status(MPI_SOURCE), tag_new_task, &
                     MPI_COMM_WORLD, status, error )
        else
            call mpi_get_count( status, MPI_INTEGER, count, error )

            number_store = min( size(prime) - number_primes, count )
            prime(number_primes+1:number_primes+number_store) = result(1:number_store)
            number_primes = number_primes + number_store

            !
            ! Signal the worker programs to stop
            !
            if ( number_primes == size(prime) ) then
                do i = 1,number_images-1
                    range = -1
                    call mpi_send( range, 2, MPI_INTEGER, i, tag_new_task, &
                             MPI_COMM_WORLD, status, error )
                enddo
                exit
            endif
        endif
    enddo

end subroutine handle_communication

!
! Get the range - if any
!
subroutine get_range
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer                             :: error

    call mpi_send( 1,     1,           MPI_INTEGER, main, tag_new_task, MPI_COMM_WORLD, status, error )
    call mpi_recv( range, size(range), MPI_INTEGER, main, tag_new_task, MPI_COMM_WORLD, status, error )

    if ( error /= MPI_SUCCESS ) then
        write(*,*) 'Some error occurred while receiving: ',error
        call mpi_finalize( error )
        stop
    endif

    !
    ! Do we have a task?
    !
    new_task = range(2) > range(1)

end subroutine get_range

!
! Subroutine to get a task and search for
! primes inside the new range
!
subroutine find_primes

    integer                 :: lower
    integer                 :: upper
    logical                 :: isprime
    logical                 :: got_task
    integer                 :: np
    integer                 :: i
    integer                 :: j
    integer                 :: residue
    integer                 :: maxindex

    np       = 0

    lower = range(1)
    upper = range(2)

    do i = lower,upper
        isprime = .true.
        do j = 2,int(sqrt(real(i)))
            residue = mod(i,j)
            if ( residue == 0 ) then
                isprime = .false.
                exit
            endif
        enddo
        if ( isprime ) then
            np = np + 1
            prime(np) = i
        endif
    enddo

    number_primes = np

end subroutine find_primes

end program primes_mpi
