! integrate_class_v3.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Test the code for the interface design problem - integration
!     routines using Fortran 2003's classes. Version 3
!
module integration_library

    implicit none

    type, abstract :: function_parameters
        procedure(eval), pointer, pass(params) :: feval
    end type function_parameters

    abstract interface
        real function eval(x, params)
            import :: function_parameters
            class(function_parameters) :: params
        end function eval
    end interface

contains

subroutine integrate_trapezoid( params, xmin, xmax, steps, result )

    interface
        real function f( x, params )
            import function_parameters
            real, intent(in)           :: x
            class(function_parameters) :: params
        end function f
    end interface

    class(function_parameters) :: params
    real, intent(in)           :: xmin, xmax
    integer, intent(in)        :: steps
    real, intent(out)          :: result

    integer                    :: i
    real                       :: x
    real                       :: deltx

    if ( steps <= 0 ) then
        result = 0.0
        return
    endif

    deltx = (xmax - xmin) / steps

    result = ( params%feval(xmin) + params%feval(xmax) )/ 2.0

    do i = 2,steps
        x      = xmin + (i - 1) * deltx
        result = result + params%feval(x)
    enddo

    result = result * deltx
end subroutine integrate_trapezoid

end module integration_library

module functions
    use integration_library

    implicit none

    type, extends(function_parameters) :: my_function_parameters
        real :: a
    end type my_function_parameters

contains
real function f( x, params )

    real, intent(in)              :: x
    class(my_function_parameters) :: params

    f = params%a * x

end function f

end module functions

program test_integrate

    use integration_library
    use functions

    implicit none

    type(my_function_parameters) :: params
    real                         :: xmin, xmax, result
    integer                      :: steps

    params%a     =  1.0
    params%feval => f
    xmin         = 1.0
    xmax         = 10.0
    steps        = 10

    call integrate_trapezoid( params, xmin, xmax, steps, result )

    write(*,*) 'Computed: ', result
    write(*,*) 'Expected: ', 0.5*(xmin+xmax)*(xmax-xmin)

end program test_integrate
