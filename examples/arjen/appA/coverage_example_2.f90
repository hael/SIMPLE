! coverage_example_2.f90 --
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

    real    :: a
    integer :: times
    real    :: b

    write(*,*) 'Times?'
    read(*,*) times

    call setvalue( a, 1, times )

    write(*,*) 'Parameter = ', a

contains
subroutine setvalue( param, type, times )

    real, intent(inout) :: param
    integer, intent(in) :: type
    integer, intent(in) :: times

    integer             :: i

    if ( type == 0 ) then
        do i = 1,times
            param = b * param * exp(-param)
        enddo
    endif

end subroutine setvalue
end program test_coverage
