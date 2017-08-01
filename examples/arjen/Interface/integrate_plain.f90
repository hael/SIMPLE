! integrate_plain.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Test the code for the interface design problem - integration
!     routines
!
module integration_library

    implicit none

contains

subroutine integrate_trapezoid( f, xmin, xmax, steps, result )

    interface
        real function f( x )
            real, intent(in) :: x
        end function f
    end interface

    real, intent(in)    :: xmin, xmax
    integer, intent(in) :: steps
    real, intent(out)   :: result

    integer             :: i
    real                :: x
    real                :: deltx

    if ( steps <= 0 ) then
        result = 0.0
        return
    endif

    deltx = (xmax - xmin) / steps

    result = ( f(xmin) + f(xmax) )/ 2.0

    do i = 2,steps
        x      = xmin + (i - 1) * deltx
        result = result + f(x)
    enddo

    result = result * deltx
end subroutine integrate_trapezoid

end module integration_library

module functions
    implicit none
contains
real function f( x )
    real, intent(in) :: x

    f = x
end function f
end module functions

program test_integrate

    use integration_library
    use functions

    implicit none

    real    :: xmin, xmax, result
    integer :: steps

    xmin  = 1.0
    xmax  = 10.0
    steps = 10

    call integrate_trapezoid( f, xmin, xmax, steps, result )

    write(*,*) 'Computed: ', result
    write(*,*) 'Expected: ', 0.5*(xmin+xmax)*(xmax-xmin)

end program test_integrate
