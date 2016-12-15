! integrate_reverse.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Integrate using reverse communication
!
module integration_library

    implicit none

    type integration_parameters
        private
        integer :: state = -1          ! Not-initialised
        integer :: steps
        integer :: i
        real    :: x
        real    :: xmin
        real    :: xmax
        real    :: deltx
        real    :: result
        real    :: sum
    end type integration_parameters

    !
    ! Parameters describing actions
    !
    integer, parameter :: get_value = 1
    integer, parameter :: completed = 2
    integer, parameter :: failure   = 3

contains

subroutine set_parameters( data, xmin, xmax, steps )

    type(integration_parameters) :: data
    real, intent(in)             :: xmin
    real, intent(in)             :: xmax
    integer, intent(in)          :: steps

    if ( steps <= 0 ) then
        return
    endif

    data%xmin   = xmin
    data%xmax   = xmax
    data%steps  = steps

    data%state  = 1
    data%sum    = 0.0
    data%i      = 0
    data%deltx  = (xmax - xmin) / steps

end subroutine set_parameters

subroutine integrate_trapezoid( data, value, result, action, x )

    type(integration_parameters) :: data
    real, intent(in)             :: value
    real, intent(out)            :: result
    integer, intent(out)         :: action
    real, intent(out)            :: x

    result  = 0.0
    if ( data%state == -1 ) then
        action = failure
        return
    endif

    !
    ! We split the computation into steps
    !
    select case ( data%state )
        case ( 1 )
            x          = data%xmin
            action     = get_value
            data%state = 2
        case ( 2 )
            data%result = 0.5 * value
            x           = data%xmax
            action      = get_value
            data%state  = 3
        case ( 3 )
            data%result = data%result + 0.5 * value
            x           = data%xmin + data%deltx
            data%i      = 1
            action      = get_value
            data%state  = 4
        case ( 4 )
            data%result = data%result + value
            if ( data%i < data%steps-1 ) then
                data%i     = data%i    + 1
                x          = data%xmin + data%i * data%deltx
                action     = get_value
                data%state = 4
            else
                result     = data%result * data%deltx
                action     = completed
            endif
        case default
            write(*,*) 'Programming error - unknown state: ', data%state
            stop
    end select

end subroutine integrate_trapezoid

end module integration_library

module functions
    implicit none
contains
real function f( x, a, b )
    real, intent(in) :: x
    real, intent(in) :: a, b  ! Not used

    f = x
end function f
end module functions

program test_integrate

    use integration_library
    use functions

    implicit none

    real    :: xmin, xmax, result, value, x
    real    :: a, b
    integer :: steps
    integer :: action

    type(integration_parameters) :: data

    a     = 1.0
    b     = 2.0

    xmin  = 1.0
    xmax  = 10.0
    steps = 10

    call set_parameters( data, xmin, xmax, steps )

    do
        call integrate_trapezoid( data, value, result, action, x )

        select case ( action )
            case ( get_value )
                value = f(x,a,b)
                write(*,*) x, value

            case ( completed )
                exit

            case ( failure )
                write(*,*) 'Error: invalid arguments'

            case default
                write(*,*) 'Programming error: unknown action - ', action
        end select
    enddo

    write(*,*) 'Computed: ', result
    write(*,*) 'Expected: ', 0.5*(xmin+xmax)*(xmax-xmin)

end program test_integrate
