! random_walk.f90 --
!     Simulate a random walk in two and three dimensions
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
module points2d3d

    implicit none

    type point2d
        real :: x, y
    contains
        procedure               :: print           => print_2d
        procedure               :: random_vector   => random_vector_2d
        procedure               :: add_vector      => add_vector_2d
        procedure, pass(vector) :: scale_by_factor => scale_by_factor_2d
        procedure               :: assign          => assign_2d
        generic                 :: operator(+)     => add_vector
        generic                 :: operator(*)     => scale_by_factor
        generic                 :: assignment(=)   => assign
    end type point2d

    type, extends(point2d) :: point3d
        real :: z
    contains
        procedure               :: print           => print_3d
        procedure               :: random_vector   => random_vector_3d
        procedure               :: add_vector      => add_vector_3d
        procedure, pass(vector) :: scale_by_factor => scale_by_factor_3d
        procedure               :: assign          => assign_3d
    end type point3d

contains

!
! Auxiliary routine to convert from class to point
!
subroutine assign_2d( point1, point2 )
    class(point2d), intent(inout) :: point1
    class(point2d), intent(in)    :: point2

    point1%x = point2%x
    point1%y = point2%y

end subroutine assign_2d

subroutine assign_3d( point1, point2 )
    class(point3d), intent(inout) :: point1
    class(point2d), intent(in)    :: point2

    point1%x = point2%x
    point1%y = point2%y
    point1%z = 0.0

    select type (point2)
        class is (point3d)
            point1%z = point2%z
    end select

end subroutine assign_3d

subroutine print_2d( point )
    class(point2d) :: point

    write(*,'(2f10.4)') point%x, point%y
end subroutine print_2d

subroutine print_3d( point )
    class(point3d) :: point

    write(*,'(3f10.4)') point%x, point%y, point%z
end subroutine print_3d

subroutine random_vector_2d( point )
    class(point2d) :: point

    call random_number( point%x )
    call random_number( point%y )

    point%x = 2.0 * (point%x - 0.5)
    point%y = 2.0 * (point%y - 0.5)

    write(*,*) 'Random: ', point%x, point%y

end subroutine random_vector_2d

subroutine random_vector_3d( point )
    class(point3d) :: point

    call point%point2d%random_vector
    call random_number( point%z )

    point%z = 2.0 * (point%z - 0.5)

end subroutine random_vector_3d

function add_vector_2d( point, vector )
    class(point2d), intent(in)  :: point, vector
    class(point2d), allocatable :: add_vector_2d

    ! Workaround for gfortran 4.6
    if ( allocated( add_vector_2d ) ) then
        deallocate( add_vector_2d )
    endif

    allocate( add_vector_2d )
    add_vector_2d%x = point%x + vector%x
    add_vector_2d%y = point%y + vector%y

    write(*,*) 'Add point: ', point%x, point%y
    write(*,*) 'Add vector: ', vector%x, vector%y
    write(*,*) 'Add result: ', add_vector_2d%x, add_vector_2d%y

end function add_vector_2d

function add_vector_3d( point, vector )
    class(point3d), intent(in)  :: point
    class(point2d), intent(in)  :: vector

    class(point2d), allocatable :: add_vector_3d
    class(point3d), allocatable :: add_result

    allocate( add_result )

    select type (vector)
        class is (point3d)
            !add_result%point2d = point%point2d + vector%point2d
            add_result%point2d = point%point2d%add_vector( vector%point2d )
            add_result%z       = point%z       + vector%z

        class default
            add_result%point2d = point%point2d + vector
            add_result%z       = 0.0
    end select

    call move_alloc( add_result, add_vector_3d )

end function add_vector_3d

function scale_by_factor_2d( factor, vector )
    real, intent(in)            :: factor
    class(point2d), intent(in)  :: vector
    class(point2d), allocatable :: scale_by_factor_2d

    ! Workaround for gfortran 4.6
    if ( allocated( scale_by_factor_2d ) ) then
        deallocate( scale_by_factor_2d )
    endif

    allocate( scale_by_factor_2d )

    scale_by_factor_2d%x = factor * vector%x
    scale_by_factor_2d%y = factor * vector%y

    write(*,*) 'Scale vector: ', vector%x, vector%y
    write(*,*) 'Scale factor: ', factor
    write(*,*) 'Scale result: ', scale_by_factor_2d%x, scale_by_factor_2d%y

end function scale_by_factor_2d

function scale_by_factor_3d( factor, vector )
    real, intent(in)            :: factor
    class(point3d), intent(in)  :: vector

    class(point2d), allocatable :: scale_by_factor_3d
    class(point3d), allocatable :: scale_result

    allocate( scale_result )

    scale_result%point2d = factor * vector%point2d
    scale_result%z       = factor * vector%z

    call move_alloc( scale_result, scale_by_factor_3d )

end function scale_by_factor_3d

end module points2d3d

program random_walk

    use points2d3d   ! Both 2D and 3D points available

    type(point2d), target   :: point_2d
    type(point3d), target   :: point_3d
    type(point2d), target   :: vector_2d
    type(point3d), target   :: vector_3d

    !
    ! A variable of class point2d can point to point_2d but
    ! also to point_3d
    !
    class(point2d), pointer :: point
    class(point2d), pointer :: vector

    integer        :: nsteps = 100
    integer        :: i
    integer        :: trial
    real           :: deltt  = 0.1

    do trial = 1,2

        ! Select what type of point ...

        if ( trial == 1 ) then
            point  => point_2d
            vector => vector_2d
            write(*,*) 'Two-dimensional walk:'
        else
            point  => point_3d
            vector => vector_3d
            write(*,*) 'Three-dimensional walk:'
        endif

        call point%random_vector

!        do i = 1,nsteps
         do i = 1,2!!
            call vector%random_vector

            point = point + deltt * vector

            call point%print
        enddo
        exit!!
    enddo
end program random_walk
