! particles.f90 --
!     Example of modelling particles transported in a flow field
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!

! particle_modelling --
!     Basic module providing the framework
!
module particle_modelling
    use points2d3d

    implicit none

    type, extends(point2d) :: basic_particle_type
        real :: mass
    contains
        procedure :: force => force_basic_particle
        procedure, pass(particle) &
                  :: new_position => position_basic_particle
    end type basic_particle_type

    type :: basic_flow_field_type
    contains
        procedure :: flow_velocity => velocity_basic_flow_field
    end type basic_flow_field_type

contains

subroutine position_basic_particle( particle, flow_field, deltt )

    class(basic_particle_type),   intent(inout), target :: particle
    class(basic_flow_field_type), intent(in)            :: flow_field
    real, intent(in)                                    :: deltt

    ! Empty routine
end subroutine position_basic_particle

subroutine force_basic_particle( particle, force_vector )

    class(basic_particle_type),   intent(inout) :: particle
    class(point2d), intent(out)                 :: force_vector

    ! Empty routine
end subroutine force_basic_particle

subroutine velocity_basic_flow_field( flow_field, position, velocity )

    class(basic_flow_field_type), intent(in)    :: flow_field
    class(point2d), intent(in)                  :: position
    class(point2d), intent(out)                 :: velocity

    ! Empty routine
end subroutine velocity_basic_flow_field
end module particle_modelling


! analytical_field_module --
!     Implementation of a simple flow field
!
module analytical_field_module
    use particle_modelling

    implicit none

    type, extends(basic_flow_field_type) :: analytical_field
        real :: x_velocity
        real :: y_velocity
    contains
        procedure :: flow_velocity => velocity_analytical_flow_field
    end type analytical_field

contains

subroutine velocity_analytical_flow_field( flow_field, position, velocity )

    class(analytical_field), intent(in)    :: flow_field
    class(point2d), intent(in)             :: position
    class(point2d), intent(out)            :: velocity

    ! The flow field is uniform, so always return the same vector
    velocity%x = flow_field%x_velocity
    velocity%y = flow_field%y_velocity

end subroutine velocity_analytical_flow_field
end module analytical_field_module


! oil_particles --
!     Naive implementation of an oil model
!
!     Note:
!     We do not use the force routine in this implementation,
!     so we can use the dummy version from the basic particle type.
!
module oil_particles
    use particle_modelling

    implicit none

    type, extends(basic_particle_type) :: oil_particle
        logical :: stuck
    contains
        procedure, pass(particle) &
                  :: new_position => position_oil_particle
    end type oil_particle

contains
subroutine position_oil_particle( particle, flow_field, deltt )
    class(oil_particle), target  :: particle
    real                         :: deltt
    class(basic_flow_field_type) :: flow_field

    class(point2d), pointer      :: position
    type(point2d)                :: velocity
    type(point2d)                :: random_displacement

    real                         :: r

    !
    ! The particle may get stuck to the bottom ...
    !
    call random_number( r )
    if ( r > 0.99 ) then
        particle%stuck = .true.
    endif
    !
    ! If it is stuck to the bottom, no further motion
    !
    if ( particle%stuck ) then
        return
    endif

    !
    ! Else let it be transported with the flow. Add a
    ! random displacement due to mixing and turbulence
    !
    position => particle%point2d

    call flow_field%flow_velocity( position, velocity )
    call random_displacement%random_vector

    position = position + deltt * velocity + &
               random_displacement

end subroutine position_oil_particle
end module oil_particles


! test_oil_particles --
!     Main program to demonstrate the oil particles modelling
!
program test_oil_particles

    use oil_particles
    use analytical_field_module

    implicit none

   ! type(oil_particle), dimension(100000) :: particle
    type(oil_particle), dimension(10)     :: particle
    type(analytical_field)                :: flow_field
    real                                  :: deltt

    integer                               :: time
    integer                               :: number_times
    integer                               :: p

    flow_field = analytical_field( 1.0, 0.1 )

    number_times = 10

    do time = 1,number_times
        write(*,*) 'Time:', time
        do p = 1,size(particle)
            call particle(p)%new_position( flow_field, deltt )
            call particle(p)%print
        enddo
    enddo

end program test_oil_particles
