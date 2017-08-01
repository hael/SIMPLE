! check_stack.f90 --
!     See what the program does with an ever increasing automatic array
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program check_stack
    implicit none

    integer :: size

    size = 1
    do
        size = size * 2
        write(*,*) 'Size: ', size

        call create_automatic_array( size )
    enddo
contains
subroutine create_automatic_array( size )
    integer               :: size
    real, dimension(size) :: data

    data(1) = 0.0
end subroutine create_automatic_array

end program check_stack
