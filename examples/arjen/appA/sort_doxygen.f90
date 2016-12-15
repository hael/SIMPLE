!> @file
!! The sort_routines module provides several routines to sort arrays
!<
!
!> @defgroup Sort library
!! The sort_routines module provides several routines to sort arrays
!!
module sort_routines
    implicit none

    !> generic name for the sorting routines
    !! The sort interface allows the programmer to use a single
    !! generic name
    !<
    interface sort
        module procedure sort_int
        module procedure sort_real
    end interface

contains

!>sort array of integers
!<
subroutine sort_int( array )
    integer, dimension(:) :: array   !< Array to be sorted

end subroutine sort_int

!>sort array of integers
!<
subroutine sort_real( array )
    real, dimension(:) :: array      !< Array to be sorted

end subroutine sort_real

end module sort_routines
