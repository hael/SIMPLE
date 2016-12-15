! runprogram.f90 --
!     Produce a table for a given function
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
program runprogram
    use iso_fortran_env

    implicit none

    real    :: param, x, y
    integer :: i

    do
        read(*,*) param

        do i = 0,200
            x = i * 0.05
            y = func( param, x )

            write(*,*) x, y
        enddo

        flush( output_unit )

    enddo
contains

real function func( param, x )
    real :: param, x

    func = exp(-x) * cos(param*x)

end function func

end program runprogram
