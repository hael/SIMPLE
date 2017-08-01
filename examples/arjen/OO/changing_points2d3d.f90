! changing_points2d3d.f90 --
!     Module for wrapping the original point2d type
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module new_points2d3d
    use points2d3d, point2d_original       => point2d, &
                    add_vector_2d_original => add_vector_2d

    implicit none

    type, extends(point2d_original) :: point2d
    contains
        procedure :: add_vector => add_vector_2d
    end type

contains

function add_vector_2d( point, vector )
    class(point2d), intent(in)           :: point
    class(point2d_original), intent(in)  :: vector
    class(point2d_original), allocatable :: add_vector_2d

    ! Workaround for gfortran 4.6
    if ( allocated( add_vector_2d ) ) then
        deallocate( add_vector_2d )
    endif

    write(*,*) 'Calling "add_vector"'

    allocate( add_vector_2d )
    add_vector_2d = point%point2d_original%add_vector(vector)

end function add_vector_2d
end module new_points2d3d

module points2d_functionality
    !
    ! Original module
    !
    !use points2d3d

    !
    ! Use the new one
    !
    use new_points2d3d
end module points2d_functionality

program test_new_points2d3d
    use points2d_functionality

    type(point2d) :: point, vector, result

    point  = point2d( 1.0, 2.0 )
    vector = point2d( 0.5, 0.5 )

    result = point%add_vector(vector)

    call result%print

end program test_new_points2d3d
