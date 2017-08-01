! timealloc.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Time the impact of allocation
!
module allocation

    implicit none

contains

subroutine alloc( sz )
    integer :: sz

    real, dimension(:), allocatable :: array

    if ( sz > 0 ) then
        allocate( array(1:sz) )
    endif
end subroutine alloc
end module allocation

program timealloc
    use allocation

    implicit none

    integer :: sz
    integer :: i
    integer :: j
    integer :: time1, time2, time3
    integer :: total_noalloc, total_alloc
    real    :: r

    total_noalloc = 0
    total_alloc   = 0
    do j = 1,10
        call system_clock( time1 )
        do i = 1,1000000
            call random_number( r )
            sz = 100000 * r
            call alloc( -1 )
        enddo
        call system_clock( time2 )
        total_noalloc = total_noalloc + (time2-time1)

        call system_clock( time1 )
        do i = 1,1000000
            call random_number( r )
            sz = 100000 * r
            call alloc( sz )
        enddo
        call system_clock( time2 )
        total_alloc = total_alloc + (time2-time1)
    enddo

    write(*,*) 'No allocation:   ', total_noalloc
    write(*,*) 'With allocation: ', total_alloc

end program
