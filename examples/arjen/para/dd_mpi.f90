! dd_openmp.f90 --
!     Idea:
!     Use domain-decomposition as an example of how to use OpenMP
!
!     We keep it simple:
!     two domains with different sizes, attached at one side
!
!     Check:
!     - Compare the output to the other methods
!     - What happens if the receive is left out? Still closing in on 1?
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program dd_openmp

    use mpi

    implicit none

    real, dimension(:,:), allocatable :: temperature
    real, dimension(:), allocatable   :: temp_interface

    integer                           :: itime
    integer                           :: ibound
    integer                           :: icopy
    integer                           :: side
    integer                           :: xmax
    integer                           :: ymax
    integer                           :: thid

    real                              :: deltt
    real                              :: coeff

    integer                             :: error
    integer                             :: rank
    integer                             :: number
    integer                             :: handle
    integer                             :: tag
    integer, dimension(MPI_STATUS_SIZE) :: status

    call mpi_init( error )
    if ( error /= MPI_SUCCESS ) then
        write(*,*) 'Error while starting MPI-based program:' ,error
        stop
    endif

    call mpi_comm_rank( MPI_COMM_WORLD, rank, error )
    call mpi_comm_size( MPI_COMM_WORLD, number, error )

    !
    ! Do not do anything if we have a rank larger than 1
    !
    if ( rank > 1 ) then
        call mpi_finalize( error )
        stop
    endif

    !
    ! Allocate the arrays we need
    !
    deltt = 0.1
    coeff = 1.0    ! Contains the thermal conductivity and grid cell size

    if ( rank == 0 ) then
        allocate( temperature(20,20) )
        temperature = 0.0

        !
        ! Left boundary value
        !
        temperature(:,1) = 1.0
        !
        ! Right interface
        !
        ibound   = size(temperature,2)
        icopy    = ibound - 1

        tag      = 1 ! Rank from which we receive data or send data to

    else
        allocate( temperature(20,30) )
        temperature = 0.0
        !
        ! Right boundary value
        !
        temperature(:,30) = 1.0
        !
        ! Left interface
        !
        ibound   = 1
        icopy    = ibound + 1

        tag      = 0 ! Rank from which we receive data or send data to
    endif

    !
    ! Allocate the array we need to transfer the temperature
    ! at the interface
    !
    allocate( temp_interface(20) )
    temp_interface = 0.0

    !
    ! From now on: compute
    !
    do itime = 1,10000

        !
        ! Set the Neumann boundary conditions
        !
        !call set_boundaries( temperature )
        side = size(temperature,1)
        temperature(1,:)    = temperature(2,:)
        temperature(side,:) = temperature(side-1,:)

        !
        ! Copy the temperature at the interface from the other image
        !
        temperature(:,ibound) = temp_interface(:)

        !
        ! Determine the new values
        !
        !call set_timestep( temperature, deltt, coeff )

        xmax = size(temperature,1) - 1
        ymax = size(temperature,2) - 1
        temperature(2:xmax,2:ymax) = temperature(2:xmax,2:ymax) + deltt * coeff * &
            ( temperature(1:xmax-1,2:ymax) + temperature(3:xmax+1,2:ymax) + &
              temperature(2:xmax,1:ymax-1) + temperature(2:xmax,3:ymax+1) &
              - 4.0*temperature(2:xmax,2:ymax) )

        !
        ! Copy the values to the other image - do not wait for an answer
        !
        temp_interface(:) = temperature(:,icopy)
        call mpi_isend( temp_interface, size(temp_interface), MPI_REAL, tag, tag, &
                 MPI_COMM_WORLD, handle, error )

        if ( error /= MPI_SUCCESS ) then
            write(*,*) 'Some error occurred while sending: ',error
            call mpi_finalize( error )
            stop
        endif

        !
        ! Receive them from the other side (use rank as the tag!)
        !
        call mpi_recv( temp_interface, size(temp_interface), MPI_REAL, tag, rank, &
                 MPI_COMM_WORLD, status, error )

        if ( error /= MPI_SUCCESS ) then
            write(*,*) 'Some error occurred while receiving: ',error
            call mpi_finalize( error )
            stop
        endif

        !
        ! Make sure all images wait for the next step
        ! - this is implicit in the fact that we have to receive the data first
        !

        write(*,*) itime, rank, temperature(10,10), temp_interface(10)
    enddo

    call mpi_finalize( error )
end program
