! prototypes_geom.f90 --
!     Using prototypes to determine the properties of an object
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module polygons

    implicit none

    type polygon_type
        real, dimension(:), allocatable   :: x, y
        !
        ! Initialisaton of a procedure pointer is a Fortran 2008 feature
        ! procedure(compute_value), pointer, pass(polygon) :: area => area_polygon
        procedure(compute_value), pointer :: area
    contains
        procedure :: draw => draw_polygon
    end type polygon_type

    abstract interface
        real function compute_value(polygon)
            import :: polygon_type
            class(polygon_type) :: polygon
        end function compute_value
    end interface

contains

subroutine draw_polygon( polygon )
    class(polygon_type) :: polygon

    integer             :: i

    write(*,*) "Draw polygon:"
    write(*,'(2f10.4)') (polygon%x(i), polygon%y(i) ,i=1,size(polygon%x))

end subroutine draw_polygon

!
! Determine area of polygon - dummy
!
real function area_polygon( polygon )
    class(polygon_type) :: polygon

    area_polygon = -huge(1.0)

end function area_polygon

!
! Alternative for rectangles: simpler method
!
real function area_rectangle( polygon )
    class(polygon_type) :: polygon

    associate( x => polygon%x, y => polygon%y )
        area_rectangle = abs( (x(2)-x(1)) * (y(3) - y(2)) )
    end associate
end function area_rectangle

subroutine new_polygon( polygon, x, y )
    real, dimension(:)               :: x, y
    type(polygon_type)               :: polygon

    allocate( polygon%x(size(x)), polygon%y(size(x)))
    polygon%x = x
    polygon%y = y

    polygon%area => area_polygon ! See remark above
end subroutine new_polygon

!
! Alternative method to construct a rectangle
! Override the default method for computing the area!
!
subroutine new_rectangle( rectangle, x1, y1, width, height )
    real                             :: x1, y1, width, height
    type(polygon_type)               :: rectangle

    allocate( rectangle%x(4), rectangle%y(4) )
    rectangle%x = (/ x1, x1+width, x1+width, x1 /)
    rectangle%y = (/ y1, y1, y1+height, y1+height /)

    rectangle%area => area_rectangle

end subroutine new_rectangle

end module polygons

! test_polygons
!     Demonstrate the use of these polygons
!
program test_polygons
    use polygons

    implicit none

    type(polygon_type) :: arbitrary_polygon
    type(polygon_type) :: rectangle

    call new_polygon(   arbitrary_polygon, (/ 1.0, 2.0, 3.0 /), (/ -1.0, 2.0, -1.0 /) )
    call new_rectangle( rectangle, 0.0, 0.0, 10.0, 3.0 )

    call arbitrary_polygon%draw
    call rectangle%draw
    write(*,*) 'Surface polygon:   ', arbitrary_polygon%area()
    write(*,*) 'Surface rectangle: ', rectangle%area()

end program test_polygons
