! robust_interp.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Robust version of interpolation
!
module interpolation

    implicit none

    type interpolation_data
        logical                         :: useable = .false.
        integer                         :: extrapolation
        real, dimension(:), allocatable :: x
        real, dimension(:), allocatable :: y
    end type interpolation_data

    integer, parameter                  :: extrapolation_none     = 0
    integer, parameter                  :: extrapolation_constant = 1
    integer, parameter                  :: extrapolation_linear   = 2

contains

function interpolation_object( x, y, extrapolation )
    type(interpolation_data)        :: interpolation_object
    real, dimension(:), intent(in)  :: x
    real, dimension(:), intent(in)  :: y
    integer, intent(in)             :: extrapolation

    integer                         :: i
    integer                         :: ierr
    integer                         :: n
    logical                         :: success

    interpolation_object%useable = .false.

    if ( allocated(interpolation_object%x) ) then
        deallocate(interpolation_object%x )
    endif
    if ( allocated(interpolation_object%y) ) then
        deallocate(interpolation_object%y )
    endif

    !
    ! Set the extrapolation method
    !
    interpolation_object%extrapolation = extrapolation_none
    if ( extrapolation == extrapolation_constant .or. &
         extrapolation == extrapolation_linear        ) then
        interpolation_object%extrapolation = extrapolation
    endif

    !
    ! Enough data? If not, simply return
    !
    if ( size(x) < 2 .or. size(y) < size(x) ) then
        return
    endif

    !
    ! Data sorted?
    !
    success = .true.

    do i = 2,size(x)
        if ( x(i) < x(i-1) ) then
            success = .false.
            exit
        endif
    enddo

    if ( .not. success ) then
        return
    endif

    !
    ! Copy the data
    !
    n = size(x)
    allocate( interpolation_object%x(n), &
              interpolation_object%y(n), stat = ierr )

    if ( ierr /= 0 ) then
        return
    endif

    !
    ! We allow array y to be larger than x, so take care of that
    !
    interpolation_object%x(1:n) = x(1:n)
    interpolation_object%y(1:n) = y(1:n)

    interpolation_object%useable = .true.


end function interpolation_object

subroutine interpolate( object, xp, estimate, success )

    type(interpolation_data)        :: object
    real, intent(in)                :: xp
    real, intent(out)               :: estimate
    logical, intent(out)            :: success

    integer                         :: i
    integer                         :: idx
    integer                         :: nd
    real                            :: dx

    estimate = 0.0
    success  = .false.

    if ( .not. object%useable ) then
        return
    endif

    !
    ! Check extrapolation
    !
    nd = size(object%x)

    if ( object%extrapolation == extrapolation_none ) then
        if ( xp < object%x(1)  ) return
        if ( xp > object%x(nd) ) return
    endif
    if ( object%extrapolation == extrapolation_constant ) then
        if ( xp < object%x(1)              ) then
            estimate = object%x(1)
            return
        endif
        if ( xp > object%x(nd) ) then
            estimate = object%x(nd)
            return
        endif
    endif

    !
    ! Search for the interval that contains xp
    ! (Linear extrapolation is taken care of
    ! automatically)
    !
    idx = nd - 1

    do i = 2,nd - 1
        if ( xp < object%x(i) ) then
            idx = i - 1
            exit
        endif
    enddo

    dx = object%x(idx+1) - object%x(idx)

    if ( dx /= 0.0 ) then
        estimate = object%y(idx) + &
            (xp - object%x(idx)) * (object%y(idx+1) - object%y(idx)) / dx
    else
        estimate = 0.5 * (object%y(idx+1) + object%y(idx))
    endif

    success  = .true.

end subroutine interpolate

real function simple_interpolate( x, y, xp )

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

    simple_interpolate = y(idx) + (xp - x(idx)) * (y(idx+1) - y(idx)) / &
                      (x(idx+1) -  x(idx))

end function simple_interpolate

end module interpolation

program test_interpolation
    use interpolation

    implicit none

    real, dimension(6) :: x = (/ 0.0, 1.0, 10.0, 10.0, 20.0, 20.0 /)
    real, dimension(6) :: y = (/ 0.0, 2.0, 20.0, 20.0, 10.0, 10.0 /)

    type(interpolation_data) :: interp

    integer            :: i
    real               :: xp
    real               :: result
    real               :: r
    logical            :: success

    interp = interpolation_object( x, y, extrapolation_constant )

    do i = 1,25
        xp = -4.0 + 1.0 * i
        call interpolate( interp, xp, result, success )
        r = simple_interpolate( x, y,  xp)
        write(*,'(3g12.4,5x,l)') xp, r, result, success
    enddo
end program
