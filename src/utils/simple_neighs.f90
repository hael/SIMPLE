module simple_neighs
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"

interface neigh_8
    module procedure neigh_8_1
    module procedure neigh_8_2
    module procedure neigh_8_3
end interface neigh_8

interface neigh_8_3D
    module procedure neigh_8_3D_1
    module procedure neigh_8_3D_2
end interface neigh_8_3D

interface neigh_4_3D
    module procedure neigh_4_3D_1
    module procedure neigh_4_3D_2
end interface neigh_4_3D

contains

    ! Returns 8-neighborhoods of the pixel position px
    ! it returns the pixel INDECES of the 8-neigh in CLOCKWISE order,
    ! starting from any 8-neigh. It doesn't consider the pixel itself.
    subroutine neigh_8_1( ldim, px, neigh_8, nsz )
        integer, intent(in)    :: ldim(3), px(3)
        integer, intent(inout) :: neigh_8(3,8)
        integer, intent(out)   :: nsz
        integer :: i, j
        i = px(1)
        j = px(2) ! 2-dim matrix asssumed
        if ( i-1 < 1 .and. j-1 < 1 ) then
            neigh_8(1:3,1) = [i+1,j,1]
            neigh_8(1:3,2) = [i+1,j+1,1]
            neigh_8(1:3,3) = [i,j+1,1]
            nsz = 3
        else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then
            neigh_8(1:3,1) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j-1,1]
            neigh_8(1:3,3) = [i,j-1,1]
            nsz = 3
        else if (j-1 < 1  .and. i+1 >ldim(1)) then
            neigh_8(1:3,3) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j+1,1]
            neigh_8(1:3,1) = [i,j+1,1]
            nsz = 3
        else if (j+1 > ldim(2) .and. i-1 < 1) then
            neigh_8(1:3,1) = [i,j-1,1]
            neigh_8(1:3,2) = [i+1,j-1,1]
            neigh_8(1:3,3) = [i+1,j,1]
            nsz = 3
        else if( j-1 < 1 ) then
            neigh_8(1:3,5) = [i-1,j,1]
            neigh_8(1:3,4) = [i-1,j+1,1]
            neigh_8(1:3,3) = [i,j+1,1]
            neigh_8(1:3,2) = [i+1,j+1,1]
            neigh_8(1:3,1) = [i+1,j,1]
            nsz = 5
        else if ( j+1 > ldim(2) ) then
            neigh_8(1:3,1) = [i-1,j,1]
            neigh_8(1:3,2) = [i-1,j-1,1]
            neigh_8(1:3,3) = [i,j-1,1]
            neigh_8(1:3,4) = [i+1,j-1,1]
            neigh_8(1:3,5) = [i+1,j,1]
            nsz = 5
        else if ( i-1 < 1 ) then
            neigh_8(1:3,1) = [i,j-1,1]
            neigh_8(1:3,2) = [i+1,j-1,1]
            neigh_8(1:3,3) = [i+1,j,1]
            neigh_8(1:3,4) = [i+1,j+1,1]
            neigh_8(1:3,5) = [i,j+1,1]
            nsz = 5
        else if ( i+1 > ldim(1) ) then
            neigh_8(1:3,1) = [i,j+1,1]
            neigh_8(1:3,2) = [i-1,j+1,1]
            neigh_8(1:3,3) = [i-1,j,1]
            neigh_8(1:3,4) = [i-1,j-1,1]
            neigh_8(1:3,5) = [i,j-1,1]
            nsz = 5
        else
            neigh_8(1:3,1) = [i-1,j-1,1]
            neigh_8(1:3,2) = [i,j-1,1]
            neigh_8(1:3,3) = [i+1,j-1,1]
            neigh_8(1:3,4) = [i+1,j,1]
            neigh_8(1:3,5) = [i+1,j+1,1]
            neigh_8(1:3,6) = [i,j+1,1]
            neigh_8(1:3,7) = [i-1,j+1,1]
            neigh_8(1:3,8) = [i-1,j,1]
            nsz = 8
        endif
    end subroutine neigh_8_1

    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! and the value of the pixel itself (the last one).
    subroutine neigh_8_2( ldim, imat, px, neigh_8, nsz )
        integer, intent(in)    :: ldim(3), imat(ldim(1),ldim(2),ldim(3)), px(3)
        integer, intent(inout) :: neigh_8(9), nsz
        integer :: i, j
        i = px(1)
        j = px(2) ! 2-dim matrix assumed
        ! identify neighborhood
        if( i-1 < 1 .and. j-1 < 1 )then                  ! NW corner
            neigh_8(1) = imat(i+1,j,1)
            neigh_8(2) = imat(i+1,j+1,1)
            neigh_8(3) = imat(i,j+1,1)
            neigh_8(4) = imat(i,j,1)
            nsz = 4
        else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then ! SE corner
            neigh_8(1) = imat(i-1,j,1)
            neigh_8(2) = imat(i-1,j-1,1)
            neigh_8(3) = imat(i,j-1,1)
            neigh_8(4) = imat(i,j,1)
            nsz = 4
        else if (j-1 < 1  .and. i+1 >ldim(1)) then       ! SW corner
            neigh_8(3) = imat(i-1,j,1)
            neigh_8(2) = imat(i-1,j+1,1)
            neigh_8(1) = imat(i,j+1,1)
            neigh_8(4) = imat(i,j,1)
            nsz = 4
        else if (j+1 > ldim(2) .and. i-1 < 1) then       ! NE corner
            neigh_8(1) = imat(i,j-1,1)
            neigh_8(2) = imat(i+1,j-1,1)
            neigh_8(3) = imat(i+1,j,1)
            neigh_8(4) = imat(i,j,1)
            nsz = 4
        else if( j-1 < 1 ) then                          ! N border
            neigh_8(5) = imat(i-1,j,1)
            neigh_8(4) = imat(i-1,j+1,1)
            neigh_8(3) = imat(i,j+1,1)
            neigh_8(2) = imat(i+1,j+1,1)
            neigh_8(1) = imat(i+1,j,1)
            neigh_8(6) = imat(i,j,1)
            nsz = 6
        else if ( j+1 > ldim(2) ) then                   ! S border
            neigh_8(1) = imat(i-1,j,1)
            neigh_8(2) = imat(i-1,j-1,1)
            neigh_8(3) = imat(i,j-1,1)
            neigh_8(4) = imat(i+1,j-1,1)
            neigh_8(5) = imat(i+1,j,1)
            neigh_8(6) = imat(i,j,1)
            nsz = 6
        else if ( i-1 < 1 ) then                         ! W border
            neigh_8(1) = imat(i,j-1,1)
            neigh_8(2) = imat(i+1,j-1,1)
            neigh_8(3) = imat(i+1,j,1)
            neigh_8(4) = imat(i+1,j+1,1)
            neigh_8(5) = imat(i,j+1,1)
            neigh_8(6) = imat(i,j,1)
            nsz = 6
        else if ( i+1 > ldim(1) ) then                   ! E border
            neigh_8(1) = imat(i,j+1,1)
            neigh_8(2) = imat(i-1,j+1,1)
            neigh_8(3) = imat(i-1,j,1)
            neigh_8(4) = imat(i-1,j-1,1)
            neigh_8(5) = imat(i,j-1,1)
            neigh_8(6) = imat(i,j,1)
            nsz = 6
        else                                              ! DEFAULT
            neigh_8(1) = imat(i-1,j-1,1)
            neigh_8(2) = imat(i,j-1,1)
            neigh_8(3) = imat(i+1,j-1,1)
            neigh_8(4) = imat(i+1,j,1)
            neigh_8(5) = imat(i+1,j+1,1)
            neigh_8(6) = imat(i,j+1,1)
            neigh_8(7) = imat(i-1,j+1,1)
            neigh_8(8) = imat(i-1,j,1)
            neigh_8(9) = imat(i,j,1)
            nsz = 9
        endif
    end subroutine neigh_8_2

    ! Returns 8-neighborhoods of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! and the value of the pixel itself (the last one).
    subroutine neigh_8_3( ldim, rmat, px, neigh_8, nsz )
        integer, intent(in)    :: ldim(3), px(3)
        real,    intent(in)    :: rmat(ldim(1),ldim(2),ldim(3))
        real,    intent(inout) :: neigh_8(9)
        integer, intent(out)   :: nsz
        integer :: i, j
        i = px(1)
        j = px(2) ! 2-dim matrix assumed
        ! identify neighborhood
        if( i-1 < 1 .and. j-1 < 1 )then                  ! NW corner
            neigh_8(1) = rmat(i+1,j,1)
            neigh_8(2) = rmat(i+1,j+1,1)
            neigh_8(3) = rmat(i,j+1,1)
            neigh_8(4) = rmat(i,j,1)
            nsz = 4
        else if (j+1 > ldim(2) .and. i+1 > ldim(1)) then ! SE corner
            neigh_8(1) = rmat(i-1,j,1)
            neigh_8(2) = rmat(i-1,j-1,1)
            neigh_8(3) = rmat(i,j-1,1)
            neigh_8(4) = rmat(i,j,1)
            nsz = 4
        else if (j-1 < 1  .and. i+1 >ldim(1)) then       ! SW corner
            neigh_8(3) = rmat(i-1,j,1)
            neigh_8(2) = rmat(i-1,j+1,1)
            neigh_8(1) = rmat(i,j+1,1)
            neigh_8(4) = rmat(i,j,1)
            nsz = 4
        else if (j+1 > ldim(2) .and. i-1 < 1) then       ! NE corner
            neigh_8(1) = rmat(i,j-1,1)
            neigh_8(2) = rmat(i+1,j-1,1)
            neigh_8(3) = rmat(i+1,j,1)
            neigh_8(4) = rmat(i,j,1)
            nsz = 4
        else if( j-1 < 1 ) then                          ! N border
            neigh_8(5) = rmat(i-1,j,1)
            neigh_8(4) = rmat(i-1,j+1,1)
            neigh_8(3) = rmat(i,j+1,1)
            neigh_8(2) = rmat(i+1,j+1,1)
            neigh_8(1) = rmat(i+1,j,1)
            neigh_8(6) = rmat(i,j,1)
            nsz = 6
        else if ( j+1 > ldim(2) ) then                   ! S border
            neigh_8(1) = rmat(i-1,j,1)
            neigh_8(2) = rmat(i-1,j-1,1)
            neigh_8(3) = rmat(i,j-1,1)
            neigh_8(4) = rmat(i+1,j-1,1)
            neigh_8(5) = rmat(i+1,j,1)
            neigh_8(6) = rmat(i,j,1)
            nsz = 6
        else if ( i-1 < 1 ) then                         ! W border
            neigh_8(1) = rmat(i,j-1,1)
            neigh_8(2) = rmat(i+1,j-1,1)
            neigh_8(3) = rmat(i+1,j,1)
            neigh_8(4) = rmat(i+1,j+1,1)
            neigh_8(5) = rmat(i,j+1,1)
            neigh_8(6) = rmat(i,j,1)
            nsz = 6
        else if ( i+1 > ldim(1) ) then                   ! E border
            neigh_8(1) = rmat(i,j+1,1)
            neigh_8(2) = rmat(i-1,j+1,1)
            neigh_8(3) = rmat(i-1,j,1)
            neigh_8(4) = rmat(i-1,j-1,1)
            neigh_8(5) = rmat(i,j-1,1)
            neigh_8(6) = rmat(i,j,1)
            nsz = 6
        else                                              ! DEFAULT
            neigh_8(1) = rmat(i-1,j-1,1)
            neigh_8(2) = rmat(i,j-1,1)
            neigh_8(3) = rmat(i+1,j-1,1)
            neigh_8(4) = rmat(i+1,j,1)
            neigh_8(5) = rmat(i+1,j+1,1)
            neigh_8(6) = rmat(i,j+1,1)
            neigh_8(7) = rmat(i-1,j+1,1)
            neigh_8(8) = rmat(i-1,j,1)
            neigh_8(9) = rmat(i,j,1)
            nsz = 9
        endif
    end subroutine neigh_8_3

    ! Returns 8-neighborhoods (in 3D they are 27) of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! of the first slice, then central slice and finally third slice.
    ! The value of the pixel itself is saved as the last one.
    ! This function is for volumes.
    subroutine neigh_8_3D_1(ldim, imat, px, neigh_8, nsz )
        integer, intent(in)    :: ldim(3), imat(ldim(1),ldim(2),ldim(3)), px(3)
        integer, intent(inout) :: neigh_8(27), nsz
        integer :: i, j, k
        logical :: i_ok, j_ok, k_ok
        i = px(1)
        j = px(2)
        k = px(3)
        i_ok = (i-1 > 0 .and. i+1 <= ldim(1))
        j_ok = (j-1 > 0 .and. j+1 <= ldim(2))
        k_ok = (k-1 > 0 .and. k+1 <= ldim(3))
        if( i-1 < 1 .and. j-1 < 1 .and. k-1 < 1 )then
            neigh_8(1) = imat(i+1,j,k)
            neigh_8(2) = imat(i+1,j+1,k)
            neigh_8(3) = imat(i,j+1,k)
            neigh_8(4) = imat(i+1,j,k+1)
            neigh_8(5) = imat(i+1,j+1,k+1)
            neigh_8(6) = imat(i,j+1,k+1)
            neigh_8(7) = imat(i,j,k+1)
            neigh_8(8) = imat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j+1 > ldim(2) .and. k+1 > ldim(3) )then
            neigh_8(1) = imat(i-1,j,k)
            neigh_8(2) = imat(i-1,j-1,k)
            neigh_8(3) = imat(i,j-1,k)
            neigh_8(4) = imat(i-1,j,k-1)
            neigh_8(5) = imat(i-1,j-1,k-1)
            neigh_8(6) = imat(i,j-1,k-1)
            neigh_8(7) = imat(i,j,k-1)
            neigh_8(8) = imat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j-1 < 1 .and. k+1 > ldim(3) )then
            neigh_8(1) = imat(i+1,j,k-1)
            neigh_8(2) = imat(i+1,j+1,k-1)
            neigh_8(3) = imat(i,j+1,k-1)
            neigh_8(4) = imat(i,j,k-1)
            neigh_8(5) = imat(i+1,j,k)
            neigh_8(6) = imat(i+1,j+1,k)
            neigh_8(7) = imat(i,j+1,k)
            neigh_8(8) = imat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j-1 < 1 .and. k-1 < 1 ) then
            neigh_8(1) = imat(i-1,j,k)
            neigh_8(2) = imat(i-1,j+1,k)
            neigh_8(3) = imat(i,j+1,k)
            neigh_8(4) = imat(i-1,j,k+1)
            neigh_8(5) = imat(i-1,j+1,k+1)
            neigh_8(6) = imat(i,j+1,k+1)
            neigh_8(7) = imat(i,j,k+1)
            neigh_8(8) = imat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > ldim(2) .and. k-1 < 1 ) then
            neigh_8(1) = imat(i+1,j,k)
            neigh_8(2) = imat(i+1,j-1,k)
            neigh_8(3) = imat(i+1,j,k+1)
            neigh_8(4) = imat(i,j,k+1)
            neigh_8(5) = imat(i,j,k)
            neigh_8(6) = imat(i+1,j-1,k+1)
            neigh_8(7) = imat(i,j-1,k+1)
            neigh_8(8) = imat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j-1 < 1 .and. k+1 > ldim(3) ) then
            neigh_8(1) = imat(i-1,j,k)
            neigh_8(2) = imat(i-1,j+1,k)
            neigh_8(3) = imat(i,j+1,k)
            neigh_8(4) = imat(i-1,j,k-1)
            neigh_8(5) = imat(i-1,j+1,k-1)
            neigh_8(6) = imat(i,j+1,k-1)
            neigh_8(7) = imat(i,j,k-1)
            neigh_8(8) = imat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > ldim(2) .and. k+1 > ldim(3) ) then
            neigh_8(1) = imat(i+1,j,k)
            neigh_8(2) = imat(i+1,j-1,k)
            neigh_8(3) = imat(i+1,j,k-1)
            neigh_8(4) = imat(i,j,k-1)
            neigh_8(5) = imat(i,j,k)
            neigh_8(6) = imat(i+1,j-1,k-1)
            neigh_8(7) = imat(i,j-1,k-1)
            neigh_8(8) = imat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j+1 > ldim(2) .and. k-1 < 1 ) then
            neigh_8(1) = imat(i-1,j,k)
            neigh_8(2) = imat(i-1,j-1,k)
            neigh_8(3) = imat(i-1,j,k+1)
            neigh_8(4) = imat(i,j,k+1)
            neigh_8(5) = imat(i,j,k)
            neigh_8(6) = imat(i-1,j-1,k+1)
            neigh_8(7) = imat(i,j-1,k+1)
            neigh_8(8) = imat(i,j-1,k)
            nsz = 8
            return
        elseif( i-1 < 1 .and. k-1 < 1 .and. j_ok )then
            neigh_8(1)  = imat(i,j-1,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i,j+1,k)
            neigh_8(4)  = imat(i+1,j-1,k)
            neigh_8(5)  = imat(i+1,j,k)
            neigh_8(6)  = imat(i+1,j+1,k)
            neigh_8(7)  = imat(i,j-1,k+1)
            neigh_8(8)  = imat(i,j,k+1)
            neigh_8(9)  = imat(i,j+1,k+1)
            neigh_8(10) = imat(i+1,j-1,k+1)
            neigh_8(11) = imat(i+1,j,k+1)
            neigh_8(12) = imat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. k+1 > ldim(3) .and. j_ok )then
            neigh_8(1)  = imat(i,j-1,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i,j+1,k)
            neigh_8(4)  = imat(i+1,j-1,k)
            neigh_8(5)  = imat(i+1,j,k)
            neigh_8(6)  = imat(i+1,j+1,k)
            neigh_8(7)  = imat(i,j-1,k-1)
            neigh_8(8)  = imat(i,j,k-1)
            neigh_8(9)  = imat(i,j+1,k-1)
            neigh_8(10) = imat(i+1,j-1,k-1)
            neigh_8(11) = imat(i+1,j,k-1)
            neigh_8(12) = imat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > ldim(1) .and. k+1 > ldim(3) .and. j_ok )then
            neigh_8(1)  = imat(i,j-1,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i,j+1,k)
            neigh_8(4)  = imat(i-1,j-1,k)
            neigh_8(5)  = imat(i-1,j,k)
            neigh_8(6)  = imat(i-1,j+1,k)
            neigh_8(7)  = imat(i,j-1,k-1)
            neigh_8(8)  = imat(i,j,k-1)
            neigh_8(9)  = imat(i,j+1,k-1)
            neigh_8(10) = imat(i-1,j-1,k-1)
            neigh_8(11) = imat(i-1,j,k-1)
            neigh_8(12) = imat(i-1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > ldim(1) .and. k-1 < 1 .and. j_ok )then
            neigh_8(1)  = imat(i,j-1,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i,j+1,k)
            neigh_8(4)  = imat(i-1,j-1,k)
            neigh_8(5)  = imat(i-1,j,k)
            neigh_8(6)  = imat(i-1,j+1,k)
            neigh_8(7)  = imat(i,j-1,k+1)
            neigh_8(8)  = imat(i,j,k+1)
            neigh_8(9)  = imat(i,j+1,k+1)
            neigh_8(10) = imat(i-1,j-1,k+1)
            neigh_8(11) = imat(i-1,j,k+1)
            neigh_8(12) = imat(i-1,j+1,k+1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k-1 < 1 .and. i_ok )then
            neigh_8(1)  = imat(i-1,j,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i+1,j,k)
            neigh_8(4)  = imat(i-1,j+1,k)
            neigh_8(5)  = imat(i,j+1,k)
            neigh_8(6)  = imat(i+1,j+1,k)
            neigh_8(7)  = imat(i-1,j,k+1)
            neigh_8(8)  = imat(i,j,k+1)
            neigh_8(9)  = imat(i+1,j,k+1)
            neigh_8(10) = imat(i-1,j+1,k+1)
            neigh_8(11) = imat(i,j+1,k+1)
            neigh_8(12) = imat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( j+1 > ldim(2) .and. k-1 < 1 .and. i_ok )then
            neigh_8(1)  = imat(i-1,j,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i+1,j,k)
            neigh_8(4)  = imat(i-1,j-1,k)
            neigh_8(5)  = imat(i,j-1,k)
            neigh_8(6)  = imat(i+1,j-1,k)
            neigh_8(7)  = imat(i-1,j,k+1)
            neigh_8(8)  = imat(i,j,k+1)
            neigh_8(9)  = imat(i+1,j,k+1)
            neigh_8(10) = imat(i-1,j-1,k+1)
            neigh_8(11) = imat(i,j-1,k+1)
            neigh_8(12) = imat(i+1,j-1,k+1)
            nsz = 12
            return
        elseif( j+1 > ldim(2) .and. k+1 > ldim(3) .and. i_ok )then
            neigh_8(1)  = imat(i-1,j,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i+1,j,k)
            neigh_8(4)  = imat(i-1,j-1,k)
            neigh_8(5)  = imat(i,j-1,k)
            neigh_8(6)  = imat(i+1,j-1,k)
            neigh_8(7)  = imat(i-1,j,k-1)
            neigh_8(8)  = imat(i,j,k-1)
            neigh_8(9)  = imat(i+1,j,k-1)
            neigh_8(10) = imat(i-1,j-1,k-1)
            neigh_8(11) = imat(i,j-1,k-1)
            neigh_8(12) = imat(i+1,j-1,k-1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k+1 > ldim(3) .and. i_ok )then
            neigh_8(1)  = imat(i-1,j,k)
            neigh_8(2)  = imat(i,j,k)
            neigh_8(3)  = imat(i+1,j,k)
            neigh_8(4)  = imat(i-1,j+1,k)
            neigh_8(5)  = imat(i,j+1,k)
            neigh_8(6)  = imat(i+1,j+1,k)
            neigh_8(7)  = imat(i-1,j,k-1)
            neigh_8(8)  = imat(i,j,k-1)
            neigh_8(9)  = imat(i+1,j,k-1)
            neigh_8(10) = imat(i-1,j+1,k-1)
            neigh_8(11) = imat(i,j+1,k-1)
            neigh_8(12) = imat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. j-1 < 1 .and. k_ok )then
            neigh_8(1)  = imat(i+1,j,k-1)
            neigh_8(2)  = imat(i+1,j+1,k-1)
            neigh_8(3)  = imat(i,j+1,k-1)
            neigh_8(4)  = imat(i,j,k-1)
            neigh_8(5)  = imat(i+1,j,k)
            neigh_8(6)  = imat(i+1,j+1,k)
            neigh_8(7)  = imat(i,j+1,k)
            neigh_8(8)  = imat(i+1,j,k+1)
            neigh_8(9)  = imat(i+1,j+1,k+1)
            neigh_8(10) = imat(i,j+1,k+1)
            neigh_8(11) = imat(i,j,k+1)
            neigh_8(12) = imat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > ldim(2) .and. i+1 > ldim(1) .and. k_ok ) then
            neigh_8(1)  = imat(i-1,j,k-1)
            neigh_8(2)  = imat(i-1,j-1,k-1)
            neigh_8(3)  = imat(i,j-1,k-1)
            neigh_8(4)  = imat(i,j,k-1)
            neigh_8(5)  = imat(i-1,j,k)
            neigh_8(6)  = imat(i-1,j-1,k)
            neigh_8(7)  = imat(i,j-1,k)
            neigh_8(8)  = imat(i-1,j,k+1)
            neigh_8(9)  = imat(i-1,j-1,k+1)
            neigh_8(10) = imat(i,j-1,k+1)
            neigh_8(11) = imat(i,j,k+1)
            neigh_8(12) = imat(i,j,k)
            nsz = 12
            return
        else if ( j-1 < 1  .and. i+1 >ldim(1) .and. k_ok ) then
            neigh_8(1)  = imat(i,j+1,k-1)
            neigh_8(2)  = imat(i-1,j+1,k-1)
            neigh_8(3)  = imat(i-1,j,k-1)
            neigh_8(4)  = imat(i,j,k-1)
            neigh_8(5)  = imat(i,j+1,k)
            neigh_8(6)  = imat(i-1,j+1,k)
            neigh_8(7)  = imat(i-1,j,k)
            neigh_8(8)  = imat(i,j+1,k+1)
            neigh_8(9)  = imat(i-1,j+1,k+1)
            neigh_8(10) = imat(i-1,j,k+1)
            neigh_8(11) = imat(i,j,k+1)
            neigh_8(12) = imat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > ldim(2) .and. i-1 < 1 .and. k_ok ) then
            neigh_8(1)  = imat(i,j-1,k-1)
            neigh_8(2)  = imat(i+1,j-1,k-1)
            neigh_8(3)  = imat(i+1,j,k-1)
            neigh_8(4)  = imat(i,j,k-1)
            neigh_8(5)  = imat(i,j-1,k)
            neigh_8(6)  = imat(i+1,j-1,k)
            neigh_8(7)  = imat(i+1,j,k)
            neigh_8(8)  = imat(i,j-1,k+1)
            neigh_8(9)  = imat(i+1,j-1,k+1)
            neigh_8(10) = imat(i+1,j,k+1)
            neigh_8(11) = imat(i,j,k+1)
            neigh_8(12) = imat(i,j,k)
            nsz = 12
            return
        else if( j-1 < 1 .and. i_ok .and. k_ok ) then
            neigh_8(1)  = imat(i+1,j,k-1)
            neigh_8(2)  = imat(i+1,j+1,k-1)
            neigh_8(3)  = imat(i,j+1,k-1)
            neigh_8(4)  = imat(i-1,j+1,k-1)
            neigh_8(5)  = imat(i-1,j,k-1)
            neigh_8(6)  = imat(i,j,k-1)
            neigh_8(7)  = imat(i+1,j,k)
            neigh_8(8)  = imat(i+1,j+1,k)
            neigh_8(9)  = imat(i,j+1,k)
            neigh_8(10) = imat(i-1,j+1,k)
            neigh_8(11) = imat(i-1,j,k)
            neigh_8(12) = imat(i+1,j,k+1)
            neigh_8(13) = imat(i+1,j+1,k+1)
            neigh_8(14) = imat(i,j+1,k+1)
            neigh_8(15) = imat(i-1,j+1,k+1)
            neigh_8(16) = imat(i-1,j,k+1)
            neigh_8(17) = imat(i,j,k+1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if ( j+1 > ldim(2) .and. i_ok .and. k_ok ) then
            neigh_8(1)  = imat(i-1,j,k-1)
            neigh_8(2)  = imat(i-1,j-1,k-1)
            neigh_8(3)  = imat(i,j-1,k-1)
            neigh_8(4)  = imat(i+1,j-1,k-1)
            neigh_8(5)  = imat(i+1,j,k-1)
            neigh_8(6)  = imat(i,j,k-1)
            neigh_8(7)  = imat(i-1,j,k)
            neigh_8(8)  = imat(i-1,j-1,k)
            neigh_8(9)  = imat(i,j-1,k)
            neigh_8(10) = imat(i+1,j-1,k)
            neigh_8(11) = imat(i+1,j,k)
            neigh_8(12) = imat(i-1,j,k+1)
            neigh_8(13) = imat(i-1,j-1,k+1)
            neigh_8(14) = imat(i,j-1,k+1)
            neigh_8(15) = imat(i+1,j-1,k+1)
            neigh_8(16) = imat(i+1,j,k+1)
            neigh_8(17) = imat(i,j,k+1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if ( i-1 < 1 .and. j_ok .and. k_ok  ) then
            neigh_8(1)  = imat(i,j-1,k-1)
            neigh_8(2)  = imat(i+1,j-1,k-1)
            neigh_8(3)  = imat(i+1,j,k-1)
            neigh_8(4)  = imat(i+1,j+1,k-1)
            neigh_8(5)  = imat(i,j+1,k-1)
            neigh_8(6)  = imat(i,j,k-1)
            neigh_8(7)  = imat(i,j-1,k)
            neigh_8(8)  = imat(i+1,j-1,k)
            neigh_8(9)  = imat(i+1,j,k)
            neigh_8(10) = imat(i+1,j+1,k)
            neigh_8(11) = imat(i,j+1,k)
            neigh_8(12) = imat(i,j-1,k+1)
            neigh_8(13) = imat(i+1,j-1,k+1)
            neigh_8(14) = imat(i+1,j,k+1)
            neigh_8(15) = imat(i+1,j+1,k+1)
            neigh_8(16) = imat(i,j+1,k+1)
            neigh_8(17) = imat(i,j,k+1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if ( i+1 > ldim(1) .and. j_ok .and. k_ok  ) then
            neigh_8(1)  = imat(i,j+1,k-1)
            neigh_8(2)  = imat(i-1,j+1,k-1)
            neigh_8(3)  = imat(i-1,j,k-1)
            neigh_8(4)  = imat(i-1,j-1,k-1)
            neigh_8(5)  = imat(i,j-1,k-1)
            neigh_8(6)  = imat(i,j,k-1)
            neigh_8(7)  = imat(i,j+1,k)
            neigh_8(8)  = imat(i-1,j+1,k)
            neigh_8(9)  = imat(i-1,j,k)
            neigh_8(10) = imat(i-1,j-1,k)
            neigh_8(11) = imat(i,j-1,k)
            neigh_8(12) = imat(i,j+1,k+1)
            neigh_8(13) = imat(i-1,j+1,k+1)
            neigh_8(14) = imat(i-1,j,k+1)
            neigh_8(15) = imat(i-1,j-1,k+1)
            neigh_8(16) = imat(i,j-1,k+1)
            neigh_8(17) = imat(i,j,k+1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if ( k-1 < 1 .and. i_ok .and. j_ok ) then
            neigh_8(1)  = imat(i-1,j-1,k)
            neigh_8(2)  = imat(i-1,j-1,k+1)
            neigh_8(3)  = imat(i-1,j,k)
            neigh_8(4)  = imat(i-1,j+1,k)
            neigh_8(5)  = imat(i-1,j+1,k+1)
            neigh_8(6)  = imat(i-1,j,k+1)
            neigh_8(7)  = imat(i,j-1,k)
            neigh_8(8)  = imat(i+1,j-1,k)
            neigh_8(9)  = imat(i+1,j,k)
            neigh_8(10) = imat(i+1,j+1,k)
            neigh_8(11) = imat(i,j+1,k)
            neigh_8(12) = imat(i,j-1,k+1)
            neigh_8(13) = imat(i+1,j-1,k+1)
            neigh_8(14) = imat(i+1,j,k+1)
            neigh_8(15) = imat(i+1,j+1,k+1)
            neigh_8(16) = imat(i,j+1,k+1)
            neigh_8(17) = imat(i,j,k+1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if ( k+1 > ldim(3) .and. i_ok .and. j_ok) then
            neigh_8(1)  = imat(i-1,j-1,k)
            neigh_8(2)  = imat(i-1,j-1,k-1)
            neigh_8(3)  = imat(i-1,j,k)
            neigh_8(4)  = imat(i-1,j+1,k)
            neigh_8(5)  = imat(i-1,j+1,k-1)
            neigh_8(6)  = imat(i-1,j,k-1)
            neigh_8(7)  = imat(i,j-1,k)
            neigh_8(8)  = imat(i+1,j-1,k)
            neigh_8(9)  = imat(i+1,j,k)
            neigh_8(10) = imat(i+1,j+1,k)
            neigh_8(11) = imat(i,j+1,k)
            neigh_8(12) = imat(i,j-1,k-1)
            neigh_8(13) = imat(i+1,j-1,k-1)
            neigh_8(14) = imat(i+1,j,k-1)
            neigh_8(15) = imat(i+1,j+1,k-1)
            neigh_8(16) = imat(i,j+1,k-1)
            neigh_8(17) = imat(i,j,k-1)
            neigh_8(18) = imat(i,j,k)
            nsz = 18
            return
        else if(i_ok .and. j_ok .and. k_ok) then
            neigh_8(1)  = imat(i-1,j-1,k-1)
            neigh_8(2)  = imat(i,j-1,k-1)
            neigh_8(3)  = imat(i+1,j-1,k-1)
            neigh_8(4)  = imat(i+1,j,k-1)
            neigh_8(5)  = imat(i+1,j+1,k-1)
            neigh_8(6)  = imat(i,j+1,k-1)
            neigh_8(7)  = imat(i-1,j+1,k-1)
            neigh_8(8)  = imat(i-1,j,k-1)
            neigh_8(9)  = imat(i,j,k-1)
            neigh_8(10) = imat(i-1,j-1,k)
            neigh_8(11) = imat(i,j-1,k)
            neigh_8(12) = imat(i+1,j-1,k)
            neigh_8(13) = imat(i+1,j,k)
            neigh_8(14) = imat(i+1,j+1,k)
            neigh_8(15) = imat(i,j+1,k)
            neigh_8(16) = imat(i-1,j+1,k)
            neigh_8(17) = imat(i-1,j,k)
            neigh_8(18) = imat(i-1,j-1,k+1)
            neigh_8(19) = imat(i,j-1,k+1)
            neigh_8(20) = imat(i+1,j-1,k+1)
            neigh_8(21) = imat(i+1,j,k+1)
            neigh_8(22) = imat(i+1,j+1,k+1)
            neigh_8(23) = imat(i,j+1,k+1)
            neigh_8(24) = imat(i-1,j+1,k+1)
            neigh_8(25) = imat(i-1,j,k+1)
            neigh_8(26) = imat(i,j,k+1)
            neigh_8(27) = imat(i,j,k)
            nsz = 27
            return
        else
            write(logfhandle, *) 'i, j, k =  ', i, j, k
            THROW_HARD('Case not covered!; neigh_8_3D')
        endif
    end subroutine neigh_8_3D_1

    ! Returns 8-neighborhoods (in 3D they are 27) of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a CLOCKWISE order, starting from any 4-neigh
    ! of the first slice, then central slice and finally third slice.
    ! The value of the pixel itself is saved as the last one.
    ! This function is for volumes.
    subroutine neigh_8_3D_2( ldim, rmat, px, neigh_8, nsz )
        integer, intent(in)    :: ldim(3), px(3)
        real,    intent(in)    :: rmat(ldim(1),ldim(2),ldim(3))
        real,    intent(inout) :: neigh_8(27)
        integer, intent(out) :: nsz
        integer :: i, j, k
        logical :: i_ok, j_ok, k_ok
        i = px(1)
        j = px(2)
        k = px(3)
        i_ok = (i-1 > 0 .and. i+1 <= ldim(1))
        j_ok = (j-1 > 0 .and. j+1 <= ldim(2))
        k_ok = (k-1 > 0 .and. k+1 <= ldim(3))
        if( i-1 < 1 .and. j-1 < 1 .and. k-1 < 1 )then
            neigh_8(1) = rmat(i+1,j,k)
            neigh_8(2) = rmat(i+1,j+1,k)
            neigh_8(3) = rmat(i,j+1,k)
            neigh_8(4) = rmat(i+1,j,k+1)
            neigh_8(5) = rmat(i+1,j+1,k+1)
            neigh_8(6) = rmat(i,j+1,k+1)
            neigh_8(7) = rmat(i,j,k+1)
            neigh_8(8) = rmat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j+1 > ldim(2) .and. k+1 > ldim(3) )then
            neigh_8(1) = rmat(i-1,j,k)
            neigh_8(2) = rmat(i-1,j-1,k)
            neigh_8(3) = rmat(i,j-1,k)
            neigh_8(4) = rmat(i-1,j,k-1)
            neigh_8(5) = rmat(i-1,j-1,k-1)
            neigh_8(6) = rmat(i,j-1,k-1)
            neigh_8(7) = rmat(i,j,k-1)
            neigh_8(8) = rmat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j-1 < 1 .and. k+1 > ldim(3) )then
            neigh_8(1) = rmat(i+1,j,k-1)
            neigh_8(2) = rmat(i+1,j+1,k-1)
            neigh_8(3) = rmat(i,j+1,k-1)
            neigh_8(4) = rmat(i,j,k-1)
            neigh_8(5) = rmat(i+1,j,k)
            neigh_8(6) = rmat(i+1,j+1,k)
            neigh_8(7) = rmat(i,j+1,k)
            neigh_8(8) = rmat(i,j,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j-1 < 1 .and. k-1 < 1 ) then
            neigh_8(1) = rmat(i-1,j,k)
            neigh_8(2) = rmat(i-1,j+1,k)
            neigh_8(3) = rmat(i,j+1,k)
            neigh_8(4) = rmat(i-1,j,k+1)
            neigh_8(5) = rmat(i-1,j+1,k+1)
            neigh_8(6) = rmat(i,j+1,k+1)
            neigh_8(7) = rmat(i,j,k+1)
            neigh_8(8) = rmat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > ldim(2) .and. k-1 < 1 ) then
            neigh_8(1) = rmat(i+1,j,k)
            neigh_8(2) = rmat(i+1,j-1,k)
            neigh_8(3) = rmat(i+1,j,k+1)
            neigh_8(4) = rmat(i,j,k+1)
            neigh_8(5) = rmat(i,j,k)
            neigh_8(6) = rmat(i+1,j-1,k+1)
            neigh_8(7) = rmat(i,j-1,k+1)
            neigh_8(8) = rmat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j-1 < 1 .and. k+1 > ldim(3) ) then
            neigh_8(1) = rmat(i-1,j,k)
            neigh_8(2) = rmat(i-1,j+1,k)
            neigh_8(3) = rmat(i,j+1,k)
            neigh_8(4) = rmat(i-1,j,k-1)
            neigh_8(5) = rmat(i-1,j+1,k-1)
            neigh_8(6) = rmat(i,j+1,k-1)
            neigh_8(7) = rmat(i,j,k-1)
            neigh_8(8) = rmat(i,j,k)
            nsz = 8
            return
        elseif(i-1 < 1 .and. j+1 > ldim(2) .and. k+1 > ldim(3) ) then
            neigh_8(1) = rmat(i+1,j,k)
            neigh_8(2) = rmat(i+1,j-1,k)
            neigh_8(3) = rmat(i+1,j,k-1)
            neigh_8(4) = rmat(i,j,k-1)
            neigh_8(5) = rmat(i,j,k)
            neigh_8(6) = rmat(i+1,j-1,k-1)
            neigh_8(7) = rmat(i,j-1,k-1)
            neigh_8(8) = rmat(i,j-1,k)
            nsz = 8
            return
        elseif(i+1 > ldim(1) .and. j+1 > ldim(2) .and. k-1 < 1 ) then
            neigh_8(1) = rmat(i-1,j,k)
            neigh_8(2) = rmat(i-1,j-1,k)
            neigh_8(3) = rmat(i-1,j,k+1)
            neigh_8(4) = rmat(i,j,k+1)
            neigh_8(5) = rmat(i,j,k)
            neigh_8(6) = rmat(i-1,j-1,k+1)
            neigh_8(7) = rmat(i,j-1,k+1)
            neigh_8(8) = rmat(i,j-1,k)
            nsz = 8
            return
        elseif( i-1 < 1 .and. k-1 < 1 .and. j_ok )then
            neigh_8(1)  = rmat(i,j-1,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i,j+1,k)
            neigh_8(4)  = rmat(i+1,j-1,k)
            neigh_8(5)  = rmat(i+1,j,k)
            neigh_8(6)  = rmat(i+1,j+1,k)
            neigh_8(7)  = rmat(i,j-1,k+1)
            neigh_8(8)  = rmat(i,j,k+1)
            neigh_8(9)  = rmat(i,j+1,k+1)
            neigh_8(10) = rmat(i+1,j-1,k+1)
            neigh_8(11) = rmat(i+1,j,k+1)
            neigh_8(12) = rmat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. k+1 > ldim(3) .and. j_ok )then
            neigh_8(1)  = rmat(i,j-1,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i,j+1,k)
            neigh_8(4)  = rmat(i+1,j-1,k)
            neigh_8(5)  = rmat(i+1,j,k)
            neigh_8(6)  = rmat(i+1,j+1,k)
            neigh_8(7)  = rmat(i,j-1,k-1)
            neigh_8(8)  = rmat(i,j,k-1)
            neigh_8(9)  = rmat(i,j+1,k-1)
            neigh_8(10) = rmat(i+1,j-1,k-1)
            neigh_8(11) = rmat(i+1,j,k-1)
            neigh_8(12) = rmat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > ldim(1) .and. k+1 > ldim(3) .and. j_ok )then
            neigh_8(1)  = rmat(i,j-1,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i,j+1,k)
            neigh_8(4)  = rmat(i-1,j-1,k)
            neigh_8(5)  = rmat(i-1,j,k)
            neigh_8(6)  = rmat(i-1,j+1,k)
            neigh_8(7)  = rmat(i,j-1,k-1)
            neigh_8(8)  = rmat(i,j,k-1)
            neigh_8(9)  = rmat(i,j+1,k-1)
            neigh_8(10) = rmat(i-1,j-1,k-1)
            neigh_8(11) = rmat(i-1,j,k-1)
            neigh_8(12) = rmat(i-1,j+1,k-1)
            nsz = 12
            return
        elseif( i+1 > ldim(1) .and. k-1 < 1 .and. j_ok )then
            neigh_8(1)  = rmat(i,j-1,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i,j+1,k)
            neigh_8(4)  = rmat(i-1,j-1,k)
            neigh_8(5)  = rmat(i-1,j,k)
            neigh_8(6)  = rmat(i-1,j+1,k)
            neigh_8(7)  = rmat(i,j-1,k+1)
            neigh_8(8)  = rmat(i,j,k+1)
            neigh_8(9)  = rmat(i,j+1,k+1)
            neigh_8(10) = rmat(i-1,j-1,k+1)
            neigh_8(11) = rmat(i-1,j,k+1)
            neigh_8(12) = rmat(i-1,j+1,k+1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k-1 < 1 .and. i_ok )then
            neigh_8(1)  = rmat(i-1,j,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i+1,j,k)
            neigh_8(4)  = rmat(i-1,j+1,k)
            neigh_8(5)  = rmat(i,j+1,k)
            neigh_8(6)  = rmat(i+1,j+1,k)
            neigh_8(7)  = rmat(i-1,j,k+1)
            neigh_8(8)  = rmat(i,j,k+1)
            neigh_8(9)  = rmat(i+1,j,k+1)
            neigh_8(10) = rmat(i-1,j+1,k+1)
            neigh_8(11) = rmat(i,j+1,k+1)
            neigh_8(12) = rmat(i+1,j+1,k+1)
            nsz = 12
            return
        elseif( j+1 > ldim(2) .and. k-1 < 1 .and. i_ok )then
            neigh_8(1)  = rmat(i-1,j,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i+1,j,k)
            neigh_8(4)  = rmat(i-1,j-1,k)
            neigh_8(5)  = rmat(i,j-1,k)
            neigh_8(6)  = rmat(i+1,j-1,k)
            neigh_8(7)  = rmat(i-1,j,k+1)
            neigh_8(8)  = rmat(i,j,k+1)
            neigh_8(9)  = rmat(i+1,j,k+1)
            neigh_8(10) = rmat(i-1,j-1,k+1)
            neigh_8(11) = rmat(i,j-1,k+1)
            neigh_8(12) = rmat(i+1,j-1,k+1)
            nsz = 12
            return
        elseif( j+1 > ldim(2) .and. k+1 > ldim(3) .and. i_ok )then
            neigh_8(1)  = rmat(i-1,j,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i+1,j,k)
            neigh_8(4)  = rmat(i-1,j-1,k)
            neigh_8(5)  = rmat(i,j-1,k)
            neigh_8(6)  = rmat(i+1,j-1,k)
            neigh_8(7)  = rmat(i-1,j,k-1)
            neigh_8(8)  = rmat(i,j,k-1)
            neigh_8(9)  = rmat(i+1,j,k-1)
            neigh_8(10) = rmat(i-1,j-1,k-1)
            neigh_8(11) = rmat(i,j-1,k-1)
            neigh_8(12) = rmat(i+1,j-1,k-1)
            nsz = 12
            return
        elseif( j-1 < 1 .and. k+1 > ldim(3) .and. i_ok )then
            neigh_8(1)  = rmat(i-1,j,k)
            neigh_8(2)  = rmat(i,j,k)
            neigh_8(3)  = rmat(i+1,j,k)
            neigh_8(4)  = rmat(i-1,j+1,k)
            neigh_8(5)  = rmat(i,j+1,k)
            neigh_8(6)  = rmat(i+1,j+1,k)
            neigh_8(7)  = rmat(i-1,j,k-1)
            neigh_8(8)  = rmat(i,j,k-1)
            neigh_8(9)  = rmat(i+1,j,k-1)
            neigh_8(10) = rmat(i-1,j+1,k-1)
            neigh_8(11) = rmat(i,j+1,k-1)
            neigh_8(12) = rmat(i+1,j+1,k-1)
            nsz = 12
            return
        elseif( i-1 < 1 .and. j-1 < 1 .and. k_ok )then
            neigh_8(1)  = rmat(i+1,j,k-1)
            neigh_8(2)  = rmat(i+1,j+1,k-1)
            neigh_8(3)  = rmat(i,j+1,k-1)
            neigh_8(4)  = rmat(i,j,k-1)
            neigh_8(5)  = rmat(i+1,j,k)
            neigh_8(6)  = rmat(i+1,j+1,k)
            neigh_8(7)  = rmat(i,j+1,k)
            neigh_8(8)  = rmat(i+1,j,k+1)
            neigh_8(9)  = rmat(i+1,j+1,k+1)
            neigh_8(10) = rmat(i,j+1,k+1)
            neigh_8(11) = rmat(i,j,k+1)
            neigh_8(12) = rmat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > ldim(2) .and. i+1 > ldim(1) .and. k_ok ) then
            neigh_8(1)  = rmat(i-1,j,k-1)
            neigh_8(2)  = rmat(i-1,j-1,k-1)
            neigh_8(3)  = rmat(i,j-1,k-1)
            neigh_8(4)  = rmat(i,j,k-1)
            neigh_8(5)  = rmat(i-1,j,k)
            neigh_8(6)  = rmat(i-1,j-1,k)
            neigh_8(7)  = rmat(i,j-1,k)
            neigh_8(8)  = rmat(i-1,j,k+1)
            neigh_8(9)  = rmat(i-1,j-1,k+1)
            neigh_8(10) = rmat(i,j-1,k+1)
            neigh_8(11) = rmat(i,j,k+1)
            neigh_8(12) = rmat(i,j,k)
            nsz = 12
            return
        else if ( j-1 < 1  .and. i+1 >ldim(1) .and. k_ok ) then
            neigh_8(1)  = rmat(i,j+1,k-1)
            neigh_8(2)  = rmat(i-1,j+1,k-1)
            neigh_8(3)  = rmat(i-1,j,k-1)
            neigh_8(4)  = rmat(i,j,k-1)
            neigh_8(5)  = rmat(i,j+1,k)
            neigh_8(6)  = rmat(i-1,j+1,k)
            neigh_8(7)  = rmat(i-1,j,k)
            neigh_8(8)  = rmat(i,j+1,k+1)
            neigh_8(9)  = rmat(i-1,j+1,k+1)
            neigh_8(10) = rmat(i-1,j,k+1)
            neigh_8(11) = rmat(i,j,k+1)
            neigh_8(12) = rmat(i,j,k)
            nsz = 12
            return
        else if ( j+1 > ldim(2) .and. i-1 < 1 .and. k_ok ) then
            neigh_8(1)  = rmat(i,j-1,k-1)
            neigh_8(2)  = rmat(i+1,j-1,k-1)
            neigh_8(3)  = rmat(i+1,j,k-1)
            neigh_8(4)  = rmat(i,j,k-1)
            neigh_8(5)  = rmat(i,j-1,k)
            neigh_8(6)  = rmat(i+1,j-1,k)
            neigh_8(7)  = rmat(i+1,j,k)
            neigh_8(8)  = rmat(i,j-1,k+1)
            neigh_8(9)  = rmat(i+1,j-1,k+1)
            neigh_8(10) = rmat(i+1,j,k+1)
            neigh_8(11) = rmat(i,j,k+1)
            neigh_8(12) = rmat(i,j,k)
            nsz = 12
            return
        else if( j-1 < 1 .and. i_ok .and. k_ok ) then
            neigh_8(1)  = rmat(i+1,j,k-1)
            neigh_8(2)  = rmat(i+1,j+1,k-1)
            neigh_8(3)  = rmat(i,j+1,k-1)
            neigh_8(4)  = rmat(i-1,j+1,k-1)
            neigh_8(5)  = rmat(i-1,j,k-1)
            neigh_8(6)  = rmat(i,j,k-1)
            neigh_8(7)  = rmat(i+1,j,k)
            neigh_8(8)  = rmat(i+1,j+1,k)
            neigh_8(9)  = rmat(i,j+1,k)
            neigh_8(10) = rmat(i-1,j+1,k)
            neigh_8(11) = rmat(i-1,j,k)
            neigh_8(12) = rmat(i+1,j,k+1)
            neigh_8(13) = rmat(i+1,j+1,k+1)
            neigh_8(14) = rmat(i,j+1,k+1)
            neigh_8(15) = rmat(i-1,j+1,k+1)
            neigh_8(16) = rmat(i-1,j,k+1)
            neigh_8(17) = rmat(i,j,k+1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if ( j+1 > ldim(2) .and. i_ok .and. k_ok ) then
            neigh_8(1)  = rmat(i-1,j,k-1)
            neigh_8(2)  = rmat(i-1,j-1,k-1)
            neigh_8(3)  = rmat(i,j-1,k-1)
            neigh_8(4)  = rmat(i+1,j-1,k-1)
            neigh_8(5)  = rmat(i+1,j,k-1)
            neigh_8(6)  = rmat(i,j,k-1)
            neigh_8(7)  = rmat(i-1,j,k)
            neigh_8(8)  = rmat(i-1,j-1,k)
            neigh_8(9)  = rmat(i,j-1,k)
            neigh_8(10) = rmat(i+1,j-1,k)
            neigh_8(11) = rmat(i+1,j,k)
            neigh_8(12) = rmat(i-1,j,k+1)
            neigh_8(13) = rmat(i-1,j-1,k+1)
            neigh_8(14) = rmat(i,j-1,k+1)
            neigh_8(15) = rmat(i+1,j-1,k+1)
            neigh_8(16) = rmat(i+1,j,k+1)
            neigh_8(17) = rmat(i,j,k+1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if ( i-1 < 1 .and. j_ok .and. k_ok  ) then
            neigh_8(1)  = rmat(i,j-1,k-1)
            neigh_8(2)  = rmat(i+1,j-1,k-1)
            neigh_8(3)  = rmat(i+1,j,k-1)
            neigh_8(4)  = rmat(i+1,j+1,k-1)
            neigh_8(5)  = rmat(i,j+1,k-1)
            neigh_8(6)  = rmat(i,j,k-1)
            neigh_8(7)  = rmat(i,j-1,k)
            neigh_8(8)  = rmat(i+1,j-1,k)
            neigh_8(9)  = rmat(i+1,j,k)
            neigh_8(10) = rmat(i+1,j+1,k)
            neigh_8(11) = rmat(i,j+1,k)
            neigh_8(12) = rmat(i,j-1,k+1)
            neigh_8(13) = rmat(i+1,j-1,k+1)
            neigh_8(14) = rmat(i+1,j,k+1)
            neigh_8(15) = rmat(i+1,j+1,k+1)
            neigh_8(16) = rmat(i,j+1,k+1)
            neigh_8(17) = rmat(i,j,k+1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if ( i+1 > ldim(1) .and. j_ok .and. k_ok  ) then
            neigh_8(1)  = rmat(i,j+1,k-1)
            neigh_8(2)  = rmat(i-1,j+1,k-1)
            neigh_8(3)  = rmat(i-1,j,k-1)
            neigh_8(4)  = rmat(i-1,j-1,k-1)
            neigh_8(5)  = rmat(i,j-1,k-1)
            neigh_8(6)  = rmat(i,j,k-1)
            neigh_8(7)  = rmat(i,j+1,k)
            neigh_8(8)  = rmat(i-1,j+1,k)
            neigh_8(9)  = rmat(i-1,j,k)
            neigh_8(10) = rmat(i-1,j-1,k)
            neigh_8(11) = rmat(i,j-1,k)
            neigh_8(12) = rmat(i,j+1,k+1)
            neigh_8(13) = rmat(i-1,j+1,k+1)
            neigh_8(14) = rmat(i-1,j,k+1)
            neigh_8(15) = rmat(i-1,j-1,k+1)
            neigh_8(16) = rmat(i,j-1,k+1)
            neigh_8(17) = rmat(i,j,k+1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if ( k-1 < 1 .and. i_ok .and. j_ok ) then
            neigh_8(1)  = rmat(i-1,j-1,k)
            neigh_8(2)  = rmat(i-1,j-1,k+1)
            neigh_8(3)  = rmat(i-1,j,k)
            neigh_8(4)  = rmat(i-1,j+1,k)
            neigh_8(5)  = rmat(i-1,j+1,k+1)
            neigh_8(6)  = rmat(i-1,j,k+1)
            neigh_8(7)  = rmat(i,j-1,k)
            neigh_8(8)  = rmat(i+1,j-1,k)
            neigh_8(9)  = rmat(i+1,j,k)
            neigh_8(10) = rmat(i+1,j+1,k)
            neigh_8(11) = rmat(i,j+1,k)
            neigh_8(12) = rmat(i,j-1,k+1)
            neigh_8(13) = rmat(i+1,j-1,k+1)
            neigh_8(14) = rmat(i+1,j,k+1)
            neigh_8(15) = rmat(i+1,j+1,k+1)
            neigh_8(16) = rmat(i,j+1,k+1)
            neigh_8(17) = rmat(i,j,k+1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if ( k+1 > ldim(3) .and. i_ok .and. j_ok) then
            neigh_8(1)  = rmat(i-1,j-1,k)
            neigh_8(2)  = rmat(i-1,j-1,k-1)
            neigh_8(3)  = rmat(i-1,j,k)
            neigh_8(4)  = rmat(i-1,j+1,k)
            neigh_8(5)  = rmat(i-1,j+1,k-1)
            neigh_8(6)  = rmat(i-1,j,k-1)
            neigh_8(7)  = rmat(i,j-1,k)
            neigh_8(8)  = rmat(i+1,j-1,k)
            neigh_8(9)  = rmat(i+1,j,k)
            neigh_8(10) = rmat(i+1,j+1,k)
            neigh_8(11) = rmat(i,j+1,k)
            neigh_8(12) = rmat(i,j-1,k-1)
            neigh_8(13) = rmat(i+1,j-1,k-1)
            neigh_8(14) = rmat(i+1,j,k-1)
            neigh_8(15) = rmat(i+1,j+1,k-1)
            neigh_8(16) = rmat(i,j+1,k-1)
            neigh_8(17) = rmat(i,j,k-1)
            neigh_8(18) = rmat(i,j,k)
            nsz = 18
            return
        else if(i_ok .and. j_ok .and. k_ok) then
            neigh_8(1)  = rmat(i-1,j-1,k-1)
            neigh_8(2)  = rmat(i,j-1,k-1)
            neigh_8(3)  = rmat(i+1,j-1,k-1)
            neigh_8(4)  = rmat(i+1,j,k-1)
            neigh_8(5)  = rmat(i+1,j+1,k-1)
            neigh_8(6)  = rmat(i,j+1,k-1)
            neigh_8(7)  = rmat(i-1,j+1,k-1)
            neigh_8(8)  = rmat(i-1,j,k-1)
            neigh_8(9)  = rmat(i,j,k-1)
            neigh_8(10) = rmat(i-1,j-1,k)
            neigh_8(11) = rmat(i,j-1,k)
            neigh_8(12) = rmat(i+1,j-1,k)
            neigh_8(13) = rmat(i+1,j,k)
            neigh_8(14) = rmat(i+1,j+1,k)
            neigh_8(15) = rmat(i,j+1,k)
            neigh_8(16) = rmat(i-1,j+1,k)
            neigh_8(17) = rmat(i-1,j,k)
            neigh_8(18) = rmat(i-1,j-1,k+1)
            neigh_8(19) = rmat(i,j-1,k+1)
            neigh_8(20) = rmat(i+1,j-1,k+1)
            neigh_8(21) = rmat(i+1,j,k+1)
            neigh_8(22) = rmat(i+1,j+1,k+1)
            neigh_8(23) = rmat(i,j+1,k+1)
            neigh_8(24) = rmat(i-1,j+1,k+1)
            neigh_8(25) = rmat(i-1,j,k+1)
            neigh_8(26) = rmat(i,j,k+1)
            neigh_8(27) = rmat(i,j,k)
            nsz = 27
            return
        else
            write(logfhandle, *) 'i, j, k =  ', i, j, k
            THROW_HARD('Case not covered!; neigh_8_3D')
        endif
    end subroutine neigh_8_3D_2

    ! Returns 4-neighborhoods (in 3D they are 6) of the pixel position px in self
    ! it returns the COORDINATES of the 8-neigh in a fixed order.
    ! The value of the pixel itself is NOT saved.
    ! This function is for volumes.
    subroutine neigh_4_3D_1(ldim, px, neigh_4, nsz)
        integer, intent(in)    :: ldim(3), px(3)
        integer, intent(inout) :: neigh_4(3,6)
        integer, intent(out)   :: nsz
        integer :: i, j, k
        i = px(1)
        j = px(2)
        k = px(3)
        if( i == 1 .and. j == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j+1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            nsz = 3
            return
        elseif( i == 1 .and. j == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,4) = [i+1,j,k]
            nsz = 4
            return
        elseif( i == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i+1,j,k]
            nsz = 4
            return
        elseif( j == 1 .and. k == 1) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j+1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
        endif
        if( i+1 == ldim(1) .and. j+1 == ldim(2) .and. k+1 == ldim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,2) = [i,j-1,k]
            neigh_4(1:3,3) = [i-1,j,k]
            nsz = 3
            return
        elseif( i+1 == ldim(1) .and. j+1 == ldim(2)) then
            neigh_4(1:3,1) = [i,j,k+1]
            neigh_4(1:3,2) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
        elseif( i+1 == ldim(1) .and. k+1 == ldim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,3) = [i,j+1,k]
            neigh_4(1:3,3) = [i,j-1,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
        elseif( j+1 == ldim(2) .and. k+1 == ldim(3)) then
            neigh_4(1:3,1) = [i,j,k-1]
            neigh_4(1:3,2) = [i,j-1,k]
            neigh_4(1:3,3) = [i+1,j,k]
            neigh_4(1:3,4) = [i-1,j,k]
            nsz = 4
            return
        endif
        neigh_4(1:3,1) = [i,j,k+1]
        neigh_4(1:3,2) = [i,j,k-1]
        neigh_4(1:3,3) = [i,j+1,k]
        neigh_4(1:3,4) = [i,j-1,k]
        neigh_4(1:3,5) = [i+1,j,k]
        neigh_4(1:3,6) = [i-1,j,k]
        nsz = 6
    end subroutine neigh_4_3D_1

    ! Returns 4-neighborhoods (in 3D they are 6) of the pixel position px in self
    ! it returns the INTENSITY values of the 8-neigh in a fixed order.
    ! The value of the pixel itself is NOT saved.
    ! This function is for volumes.
    subroutine neigh_4_3D_2(ldim, imat, px, neigh_4, nsz )
        integer, intent(in)    :: ldim(3), imat(ldim(1),ldim(2),ldim(3)), px(3)
        integer, intent(inout) :: neigh_4(6), nsz
        integer :: i, j, k
        i = px(1)
        j = px(2)
        k = px(3)
        if(i+1<ldim(1) .and. i-1>0 .and. j+1<ldim(2) .and. j-1>0 .and. k+1<ldim(3) .and. k-1>0) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(2) = imat(i,j,k-1)
            neigh_4(3) = imat(i,j+1,k)
            neigh_4(4) = imat(i,j-1,k)
            neigh_4(5) = imat(i+1,j,k)
            neigh_4(6) = imat(i-1,j,k)
            nsz = 6
            return
        endif
        if( i == 1 .and. j == 1 .and. k == 1) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(2) = imat(i,j+1,k)
            neigh_4(3) = imat(i+1,j,k)
            nsz = 3
            return
        elseif( i == 1 .and. j == 1) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(2) = imat(i,j,k-1)
            neigh_4(3) = imat(i,j+1,k)
            neigh_4(4) = imat(i+1,j,k)
            nsz = 4
            return
        elseif( i == 1 .and. k == 1) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(3) = imat(i,j+1,k)
            neigh_4(3) = imat(i,j-1,k)
            neigh_4(4) = imat(i+1,j,k)
            nsz = 4
            return
        elseif( j == 1 .and. k == 1) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(2) = imat(i,j+1,k)
            neigh_4(3) = imat(i+1,j,k)
            neigh_4(4) = imat(i-1,j,k)
            nsz = 4
            return
        endif
        if( i+1 == ldim(1) .and. j+1 == ldim(2) .and. k+1 == ldim(3)) then
            neigh_4(1) = imat(i,j,k-1)
            neigh_4(2) = imat(i,j-1,k)
            neigh_4(3) = imat(i-1,j,k)
            nsz = 3
            return
        elseif( i+1 == ldim(1) .and. j+1 == ldim(2)) then
            neigh_4(1) = imat(i,j,k+1)
            neigh_4(2) = imat(i,j,k-1)
            neigh_4(3) = imat(i,j-1,k)
            neigh_4(4) = imat(i-1,j,k)
            nsz = 4
            return
        elseif( i+1 == ldim(1) .and. k+1 == ldim(3)) then
            neigh_4(1) = imat(i,j,k-1)
            neigh_4(3) = imat(i,j+1,k)
            neigh_4(3) = imat(i,j-1,k)
            neigh_4(4) = imat(i-1,j,k)
            nsz = 4
            return
        elseif( j+1 == ldim(2) .and. k+1 == ldim(3)) then
            neigh_4(1) = imat(i,j,k-1)
            neigh_4(2) = imat(i,j-1,k)
            neigh_4(3) = imat(i+1,j,k)
            neigh_4(4) = imat(i-1,j,k)
            nsz = 4
            return
        endif
    end subroutine neigh_4_3D_2

end module simple_neighs
