!!****h* Utilities/sort_routines
! FUNCTION
!     The sort_routines module provides several routines to sort arrays
!!****
!
module sort_routines
    implicit none

    !****m* sort_routines/sort
    ! NAME
    !    sort - generic name for the sorting routines
    ! PURPOSE
    !    The sort interface allows the programmer to use a single
    !    generic name
    !****
    !
    interface sort
        module procedure sort_int
        module procedure sort_real
    end interface

contains

!****if* sort_routines/sort_int
! NAME
!    sort_int - sort array of integers
! SYNOPSIS
subroutine sort_int( array )
    integer, dimension(:) :: array
!****
    ... actual code ...
end subroutine sort_int

!****if* sort_routines/sort_real
! NAME
!    sort_int - sort array of reals
! SYNOPSIS
subroutine sort_real( array )
    real, dimension(:) :: array
!****
    ... actual code ...
end subroutine sort_real

end module sort_routines
