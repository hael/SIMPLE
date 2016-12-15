! timearrays.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Time the impact of various types of arrays
!
module array_types

    implicit none

contains

subroutine plain_arrays( data, sz, mean )
    real, dimension(:) :: data
    integer            :: sz
    real               :: mean

    integer            :: i

    do i = 1,sz
        data(i) = i
    enddo

    mean = sum( data(1:sz) ) / sz

end subroutine plain_arrays

subroutine auto_arrays( sz, mean )
    integer             :: sz
    real, dimension(sz) :: data     ! Local
    real                :: mean

    integer             :: i

    do i = 1,sz
        data(i) = i
    enddo

    mean = sum( data(1:sz) ) / sz

end subroutine auto_arrays

subroutine alloc_arrays( sz, mean )
    real, dimension(:), allocatable :: data     ! Local
    integer                         :: sz
    real                            :: mean

    integer                          :: i

    allocate( data(1:sz) )
    do i = 1,sz
        data(i) = i
    enddo

    mean = sum( data(1:sz) ) / sz

end subroutine alloc_arrays

subroutine pointer_arrays( sz, mean )
    real, dimension(:), pointer     :: data     ! Local
    integer                         :: sz
    real                            :: mean

    integer                         :: i

    allocate( data(1:sz) )
    do i = 1,sz
        data(i) = i
    enddo

    mean = sum( data(1:sz) ) / sz

    deallocate( data )

end subroutine pointer_arrays

end module array_types

program timearrays
    use array_types

    implicit none

    real, dimension(100000), save :: rdata
    real                          :: rmean
    integer                       :: sz
    integer                       :: total_plain
    integer                       :: total_auto
    integer                       :: total_alloc
    integer                       :: total_pointer
    integer                       :: time1, time2
    integer                       :: i
    integer                       :: j

    sz = 1
    do j = 1,4

        total_plain   = 0
        total_auto    = 0
        total_alloc   = 0
        total_pointer = 0

        sz = sz * 10
        call system_clock( time1 )
        do i = 1,1000000
            call plain_arrays( rdata, sz, rmean )
        enddo
        call system_clock( time2 )
        total_plain = total_plain + (time2-time1)

        call system_clock( time1 )
        do i = 1,1000000
            call auto_arrays( sz, rmean )
        enddo
        call system_clock( time2 )
        total_auto = total_auto + (time2-time1)

        call system_clock( time1 )
        do i = 1,1000000
            call alloc_arrays( sz, rmean )
        enddo
        call system_clock( time2 )
        total_alloc = total_alloc + (time2-time1)

        call system_clock( time1 )
        do i = 1,1000000
            call pointer_arrays( sz, rmean )
        enddo
        call system_clock( time2 )
        total_pointer = total_pointer + (time2-time1)

        write(*,*) 'Size: ', sz
        write(*,*) 'Plain:           ', total_plain
        write(*,*) 'Automatic:       ', total_auto
        write(*,*) 'Allocate:        ', total_alloc
        write(*,*) 'Pointer:         ', total_pointer
    enddo

end program
