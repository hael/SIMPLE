! interp.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Straightforward version of interpolation
!
module interpolation

    implicit none

contains

real function interpolate( x, y, xp )

    real, dimension(:), intent(in) :: x
    real, dimension(:), intent(in)  :: y
    real, intent(in)                :: xp

    integer                         :: i
    integer                         :: idx

    !
    ! Search for the interval that contains xp
    !
    idx = size(x) - 1

    do i = 2,size(x)-1
        if ( xp < x(i) ) then
            idx = i - 1
            exit
        endif
    enddo

    interpolate = y(idx) + (xp - x(idx)) * (y(idx+1) - y(idx)) / &
                      (x(idx+1) -  x(idx))

end function interpolate

end module interpolation

program test_interpolation
    use interpolation

    implicit none

    real, dimension(6) :: x = (/ 0.0, 1.0, 10.0, 10.0, 20.0, 120.0 /)
    real, dimension(6) :: y = (/ 0.0, 2.0, 20.0, 20.0, 10.0, 10.0 /)

    integer            :: i
    real               :: xp

    do i = 1,25
        xp = -4.0 + 1.0 * i
        write(*,'(2g12.4)') xp, interpolate( x, y, xp )
    enddo
end program
