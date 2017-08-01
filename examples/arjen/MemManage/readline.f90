! readline.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Read complete lines from file
!
module readline_utility
    use iso_fortran_env

    implicit none

contains

subroutine readline( lun, line, success )

    integer, intent(in)                        :: lun
    character(len=:), allocatable, intent(out) :: line
    logical, intent(out)                       :: success

    character(len=0)                           :: newline

    success = .true.

    call readline_piece_by_piece( newline )

contains

recursive subroutine readline_piece_by_piece( newline )
    character(len=*)                :: newline

    character(len=10)               :: piece
    integer                         :: ierr
    integer                         :: sz

    read( lun, '(a)', advance = 'no', size = sz, iostat = ierr ) piece

    if ( ierr /= 0 .and. ierr /= iostat_eor  ) then
        allocate( character(len=len(newline)):: line )
        line = newline
        success = .false.
        return
    endif

    !
    ! Have we gotten to the end of the line or not?
    !
    if ( sz >= len(piece)  ) then
        call readline_piece_by_piece( newline // piece )
    else
        allocate( character(len=len(newline)+sz):: line )
        line = newline // piece(1:sz)
        success = .true.
    endif
end subroutine readline_piece_by_piece
end subroutine readline

end module readline_utility

program test_readline
    use readline_utility

    integer                       :: lun
    logical                       :: success
    character(len=:), allocatable :: line

    lun = 10
    open( lun, file = 'test_readline.inp' )

    do
        call readline( lun, line, success )

        if ( .not. success ) then
            exit
        endif

        write(*,*) len(line), '>', line, '<'
        deallocate( line )
    enddo
end program test_readline
