! abstract_point.f90 --
!     Example of the use of abstract types
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module abstract_points
    implicit none

    type, abstract :: abstract_point
        ! No coordinates, leave that to the extending types
    contains
        procedure(add_vector), deferred :: add
    end type abstract_point

    !
    ! Define what the named interface "add_vector" should look like
    !
    abstract interface
        subroutine add_vector( point, vector )
            import abstract_point
            class(abstract_point), intent(inout) :: point
            class(abstract_point), intent(in)    :: vector
        end subroutine add_vector
    end interface
end module abstract_points

! points2d --
!     Use the abstract points module for a concrete type
!
module points2d
    use abstract_points

    type, extends(abstract_point) :: point2d
        real :: x, y
    contains
        procedure :: add => add_vector_2d
    end type point2d
contains

subroutine add_vector_2d( point, vector )
    class(point2d), intent(inout)        :: point
    class(abstract_point), intent(in)    :: vector

    select type (vector)
        class is (point2d)
            point%x = point%x + vector%x
            point%y = point%y + vector%y
    end select
end subroutine add_vector_2d
end module points2d

! test_points --
!     Use the abstract type
!
program test_points
    use points2d

    class(abstract_point), pointer :: p
    type(point2d), target          :: point
    type(point2d)                  :: vector

    point  = point2d(1.0,2.0)
    vector = point2d(0.5,0.5)

    p => point
    call p%add( vector )

    write(*,*) 'Resulting point:', point

end program test_points
