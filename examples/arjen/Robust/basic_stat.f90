! basic_stat.f90 --
!
!     Example belonging to "Modern Fortran in Practice" by Arjen Markus
!
!     This work is licensed under the Creative Commons Attribution 3.0 Unported License.
!     To view a copy of this license, visit http://creativecommons.org/licenses/by/3.0/
!     or send a letter to:
!     Creative Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
!
!     Basic statistical parameters - straightforward version
!
program basic_stat
    implicit none

    real    :: value
    real    :: sumsq
    real    :: sum
    real    :: stdev
    real    :: stdev2
    real    :: vmean
    real    :: var
    integer :: i
    integer :: j
    integer :: nodata
    integer :: nomissing
    integer :: ierr

    open( 10, file = 'basic_stat.data', status = 'old', iostat = ierr )

    if ( ierr /= 0 ) then
        write(*,*) 'Error opening file with input data - basic_stat.data'
        write(*,*) 'Check that it exists'
        stop
    endif

    !
    ! One value per line, ? means a missing value ...
    ! (As a ? can not be read into a number, we consider each line that causes
    ! a read error to be a missing value)
    !
    sum       = 0.0
    sumsq     = 0.0
    nodata    = 0
    nomissing = 0

    do
        read( 10, *, iostat = ierr ) value

        if ( ierr < 0 ) then
            !
            ! End of file
            !
            exit
        elseif ( ierr > 0 ) then
            !
            ! Missing value
            !
            nomissing = nomissing + 1
            cycle
        endif

        sum    = sum    + value
        sumsq  = sumsq  + value ** 2
        nodata = nodata + 1
    enddo

    close( 10 )

    !
    ! Report our findings
    !
    write(*,*) 'Outcome:'
    write(*,*) '    Number of valid data:  ', nodata
    write(*,*) '    Number of missing data:', nomissing
    write(*,*) ' '

    if ( nodata > 0 ) then
        vmean = sum / nodata
        write(*,*) '    Mean value:            ', vmean

        if ( nodata > 1 ) then
            stdev2 = (sumsq - sum ** 2 / nodata) / (nodata-1)
            write(*,*) stdev2
            stdev  = sqrt( (sumsq - sum ** 2 / nodata) / (nodata-1) )
            write(*,*) '    Standard deviation:    ', stdev
        else
            write(*,*) '    Standard deviation:    too few data'
        endif
    else
        write(*,*) '    Mean value:            too few data'
    endif
end program basic_stat
