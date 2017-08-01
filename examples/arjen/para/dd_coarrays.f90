! dd_coarrays.f90 --
!     Idea:
!     Use domain-decomposition as an example of how to use coarrays
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
program coarrays

    implicit none

    real, dimension(:,:), allocatable               :: temperature
    real, dimension(:), codimension[:], allocatable :: temp_interface

    integer                                         :: itime
    integer                                         :: ibound
    integer                                         :: icopy
    integer                                         :: toimage
    integer                                         :: side
    integer                                         :: xmax
    integer                                         :: ymax

    real                                            :: deltt
    real                                            :: coeff

    !
    ! Check the number of images - the program only works correctly
    ! with two images!
    !
    if ( num_images() /= 2 ) then
        write(*,*) 'This program requires exactly two images!'
        stop
    endif

    !
    ! Allocate the arrays we need
    !
    sync all

    deltt = 0.1
    coeff = 1.0    ! Contains the thermal conductivity and grid cell size

    if ( this_image() == 1 ) then
        allocate( temperature(20,20) )
        temperature = 0.0

        !
        ! Left boundary value
        !
        temperature(:,1) = 1.0
        !
        ! Right interface
        !
        ibound  = size(temperature,2)
        icopy   = ibound - 1
        toimage = 2
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
        ibound  = 1
        icopy   = ibound + 1
        toimage = 1
    endif

    !
    ! Allocate the one coarray we need to transfer the temperature
    ! at the interface
    !
    allocate( temp_interface(20)[*] )
    temp_interface = 0.0

    !
    ! Wait until both images are done with the initialisation
    !
    sync all

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
        temperature(:,ibound) = temp_interface

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
        temp_interface(:)[toimage] = temperature(:,icopy)

        !
        ! Make sure all images wait for the next step
        !
        sync all
        write(*,*) itime, this_image(), temperature(10,10), temp_interface(10)
    enddo

!contains
!    ...
end program
