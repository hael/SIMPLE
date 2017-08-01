! robust_stat.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Basic statistical parameters - robust version
!
program robust_stat

    character(len=80) :: line
    real              :: value
    real              :: sumsq
    real              :: sum
    real              :: var
    integer           :: i
    integer           :: j
    integer           :: nodata
    integer           :: nomissing
    integer           :: noerrors
    integer           :: noempty
    integer           :: nolines
    integer           :: ierr
    logical           :: first_value = .true.

    open( 10, file = 'robust_stat.data', status = 'old', iostat = ierr )

    if ( ierr /= 0 ) then
        write(*,*) 'Error opening file with input data - robust_stat.data'
        write(*,*) 'Check that it exists'
        stop
    endif

    !
    ! One value per line, ? means a missing value ...
    ! (Any other value that can not be read is regarded to be an error,
    ! empty lines are reported but not counted)
    !
    sum       = 0.0
    sumsq     = 0.0
    nodata    = 0
    nomissing = 0
    nolines   = 0
    noerrors  = 0

    do
        read( 10, '(a)', iostat = ierr ) line

        if ( ierr < 0 ) then
            !
            ! End of file
            !
            exit
        elseif ( ierr > 0 ) then
            !
            ! Some reading error occurred - report it
            !
            write(*,*) 'Error reading line no.', nolines+1
            write(*,*) 'Skipping the rest of the file'
            exit
        else
            !
            ! Get rid of tabs and carriage returns
            !
            call cleanup_line( line )
            !
            ! Analyse the contents:
            ! - Empty line?
            ! - Missing value?
            ! - Not a valid number?
            ! - Valid number
            !
            nolines = nolines + 1

            if ( line == ' ' ) then
                noempty = noempty + 1
                cycle
            endif
            if ( adjustl(line) == '?' ) then
                nomissing = nomissing + 1
                cycle
            endif

            read( line, *, iostat = ierr ) value

            if ( ierr /= 0 ) then
                noerrors = noerrors + 1
                cycle
            endif

            !
            ! If the value is out of range, report it and
            ! skip it
            !
            if ( abs(value) > sqrt(huge(value)) .or. &
                 ( abs(value) < sqrt(tiny(value)) .and. abs(value) /= 0.0 ) ) then
                write(*,*) 'Value out of range: ', value, ' - ignoring it!'
                nomissing = nomissing + 1
                cycle
            endif

            !
            ! We do have a valid value
            !
            if ( first_value ) then
                first_value = .false.
                offset     = value
            endif

            sum    = sum    + (value - offset)
            sumsq  = sumsq  + (value - offset) ** 2
            nodata = nodata + 1
        endif
    enddo

    close( 10 )

    !
    ! Report our findings
    !
    write(*,*) 'Outcome:'
    write(*,*) '    Number of lines read:  ', nolines
    write(*,*) '    Number of empty lines: ', noempty
    write(*,*) '    Number of valid data:  ', nodata
    write(*,*) '    Number of missing data:', nomissing
    write(*,*) '    Number of invalid data:', noerrors
    write(*,*) ' '

    if ( nodata > 0 ) then
        vmean = offset + sum / nodata
        write(*,*) '    Mean value:            ', vmean

        if ( nodata > 1 ) then
            stdev  = sqrt( (sumsq - sum ** 2 / nodata) / (nodata-1) )
            write(*,*) '    Standard deviation:    ', stdev
        else
            write(*,*) '    Standard deviation:    too few data'
        endif
    else
        write(*,*) '    Mean value:            too few data'
    endif
contains

subroutine cleanup_line( line )
    character(len=*), intent(inout) :: line

    logical                         :: found
    integer                         :: k
    integer                         :: i
    integer, dimension(3)           :: chars = (/ 9, 10, 13 /)

    found = .true.
    do while ( found )

        found = .false.
        !
        ! Remove any tabs, carriage returns and newlines
        !
        do i = 1,size(chars)
            k = index( line, achar(chars(i)) )
            if ( k > 0 ) then
                found     = .true.
                line(k:k) = ' '
            endif
        enddo
    endif

end subroutine cleanup_line
end program robust_stat
