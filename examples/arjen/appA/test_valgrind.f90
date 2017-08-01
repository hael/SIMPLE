! test_valgrind.f90 --
!     Check if valgrind can report some obvious allocation errors
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program test_valgrind

    integer, dimension(:), pointer :: data

    allocate( data(100) )

    data(1) = 1.0
    data(101) = 2.0
    nullify( data )
end program test_valgrind
