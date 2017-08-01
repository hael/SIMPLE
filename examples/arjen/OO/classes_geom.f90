! classes_geom.f90 --
!     Classes of geometrical objects as an illustration
!     of object-oriented programming in Fortran
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module geometrical_objects

    implicit none

    !
    ! General shape
    !
    type, abstract :: shape
        ! No data
    contains
        procedure(get_shape_area), deferred           :: get_area
        !procedure                                    :: size  -- does not work
    end type shape

    abstract interface
        real function get_shape_area( this )
            import                   :: shape
            class(shape), intent(in) :: this
        end function get_shape_area
    end interface

    !
    ! Rectangle
    !
    type, extends(shape) :: rectangle
        real :: width, height
    contains
        procedure :: get_area             => get_rectangle_area
        procedure :: size                 => rectangle_size
    end type rectangle

    !
    ! Square
    ! Note:
    ! square_size must have the same interface as its rectangle parent!
    !
    type, extends(rectangle) :: square
    contains
        procedure :: get_area => get_square_area
        procedure :: size     => square_size
    end type square

contains

!
! The various routines and functions we need
!
real function get_rectangle_area( this )
    class(rectangle), intent(in) :: this

    get_rectangle_area = this%width * this%height

end function get_rectangle_area

subroutine rectangle_size( this, width, height )
    class(rectangle), intent(inout) :: this
    real, intent(in)                :: width
    real, intent(in), optional      :: height

    this%width  = width
    if ( present(height) ) then
        this%height = height
    else
        this%height = width
    endif

end subroutine rectangle_size

subroutine square_size( this, width, height )
    class(square), intent(inout) :: this
    real, intent(in)             :: width
    real, intent(in), optional   :: height  ! Ignored

    this%width  = width
    this%height = 0.0

end subroutine square_size

real function get_square_area( this )
    class(square), intent(in) :: this

    get_square_area = this%width ** 2

end function get_square_area

end module geometrical_objects

!
! Small test program
program test_objects

    use geometrical_objects

    implicit none

    type list_of_objects
        class(shape), pointer :: object
    end type
    type(list_of_objects), dimension(2) :: list

    type(rectangle), target    :: rect
    type(square), target       :: sq


    integer                    :: i

    call rect%size( 1.0, 2.0 )
    call sq%size( 1.5 )

    list(1)%object => rect
    list(2)%object => sq

    do i = 1,size(list)
        write(*,*) 'Area: ', list(i)%object%get_area()
    enddo

end program test_objects
