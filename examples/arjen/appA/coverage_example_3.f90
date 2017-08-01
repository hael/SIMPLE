! coverage_example_3.f90 --
!     Show the effect of static analysis
!
!     ifort -Qdiag-enable:sc3
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program test_coverage
    implicit none

    real    :: a
    real    :: b

    call setvalue( a, 1 )

    write(*,*) 'Parameter = ', a

contains
subroutine setvalue( param, type )

    real, intent(inout) :: param
    integer, intent(in) :: type

    if ( type == 0 ) then
        param = b * param * exp(-param)
    endif

end subroutine setvalue
end program test_coverage
