! dd_openmp.f90 --
!     Idea:
!     Use domain-decomposition as an example of how to use OpenMP
!
!     We keep it simple:
!     two domains with different sizes, attached at one side
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program dd_openmp

    use omp_lib

    implicit none

    type domain_data
        real, dimension(:,:), allocatable :: temperature
        integer                           :: ibound
        integer                           :: icopy
        integer                           :: todomain
    end type

    type(domain_data), dimension(2), target :: domain

    real, dimension(:,:), allocatable :: temp_interface
    real, dimension(:,:), pointer     :: temperature

    integer                           :: itime
    integer                           :: ibound
    integer                           :: icopy
    integer                           :: todomain
    integer                           :: side
    integer                           :: xmax
    integer                           :: ymax
    integer                           :: thid

    real                              :: deltt
    real                              :: coeff

    !
    ! Allocate the arrays we need
    !
    deltt = 0.1
    coeff = 1.0    ! Contains the thermal conductivity and grid cell size

    call omp_set_num_threads( 2 )

!$omp parallel

    write(*,*) 'In thread' , omp_get_thread_num()

    if ( omp_get_thread_num() == 0 ) then
        write(*,*) 'Thread 1'
        allocate( domain(1)%temperature(20,20) )
        domain(1)%temperature = 0.0

        !
        ! Left boundary value
        !
        domain(1)%temperature(:,1) = 1.0
        !
        ! Right interface
        !
        domain(1)%ibound   = size(domain(1)%temperature,2)
        domain(1)%icopy    = domain(1)%ibound - 1
        domain(1)%todomain = 2
    else
        write(*,*) 'Thread 2'
        allocate( domain(2)%temperature(20,30) )
        domain(2)%temperature = 0.0
        !
        ! Right boundary value
        !
        domain(2)%temperature(:,30) = 1.0
        !
        ! Left interface
        !
        domain(2)%ibound   = 1
        domain(2)%icopy    = domain(2)%ibound + 1
        domain(2)%todomain = 1
    endif
!$omp end parallel

    !
    ! Allocate the array we need to transfer the temperature
    ! at the interface
    !
    allocate( temp_interface(20,2) )
    temp_interface = 0.0

    !
    ! From now on: compute
    !
!$omp parallel private( thid, side, ibound, icopy, xmax, ymax, todomain, &
!$omp                   temperature )

    do itime = 1,10000
        thid = 1 + omp_get_thread_num()

        temperature => domain(thid)%temperature
        ibound      =  domain(thid)%ibound
        icopy       =  domain(thid)%icopy
        todomain    =  domain(thid)%todomain

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
        temperature(:,ibound) = temp_interface(:,thid)

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
        ! Copy the values to the other image
        !
        temp_interface(:,todomain) = temperature(:,icopy)

        !
        ! Make sure all images wait for the next step
        !

        !$omp barrier

        write(*,*) itime, thid, temperature(10,10), temp_interface(10,thid)
    enddo
!$omp end parallel
    stop
end program
