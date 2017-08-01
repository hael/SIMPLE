! coverage_example.f90 --
!     Show how to use the gcov utility
!
!     set GCOV_PREFIX_STRIP=1000
!     gfortran --coverage | -ftest-coverage -fprofile-arcs
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

    real :: a = 0.0

    call setvalue( a, 1 )

    write(*,*) 'Parameter = ', a

contains
subroutine setvalue( param, type )

    real, intent(inout) :: param
    integer, intent(in) :: type

    if ( type == 0 ) then
        param = param * exp(-param)
    endif

end subroutine setvalue
end program test_coverage
