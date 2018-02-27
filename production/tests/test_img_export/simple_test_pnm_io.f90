!! Based on John Burkardt's  PPMA and PBMA modules (LPGL Copyright, March 2003- 2007)
module simple_test_pnm_io
    include 'simple_lib.f08'
    use simple_pnm
    implicit none

!    public :: test_pnm_io

contains
    subroutine test_pnm_io
        use simple_pnm
        call simple_timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_IO '
        call ppma_test01 ( )
        call ppma_test02 ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PPMA_IO'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call simple_timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_IO'

        call pgma_test01 ( )
        call pgma_test02 ( )
        call pgma_test03 ( )
        !
        !  Terminate.
        !
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_IO'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call simple_timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_IO portable grayscale image'

        call pgma_test01 ( )
        call pgma_test02 ( )
        !
        !  Terminate.
        !
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_IO'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call simple_timestamp ( )


    end subroutine test_pnm_io

    subroutine ppma_test01 ( )
        integer  , parameter :: ncol = 300
        integer  , parameter :: nrow = 300
        integer   :: b(nrow,ncol)
        character ( len = 80 ) :: file_name = 'test01.ascii.ppm'
        integer   :: g(nrow,ncol)
        integer   :: ierror
        integer   :: r(nrow,ncol)
        write(*,'(a)') ' '
        write(*,'(a)') 'TEST01'
        write(*,'(a)') '  PPMA_EXAMPLE sets up sample PPMA data.'
        write(*,'(a)') '  PPMA_WRITE writes an ASCII PPMA file.'
        call ppma_example ( nrow, ncol, r, g, b )
        call ppma_write ( file_name, nrow, ncol, r, g, b, ierror )
        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'TEST01 - Fatal error!'
            write ( *, '(a,i6)' ) 'PPMA_WRITE returns IERROR = ', ierror
            return
        end if
        write(*,'(a)') ' '
        write(*,'(a)') '  Wrote the header and data for "' &
            // trim ( file_name ) //'".'
        write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
        write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol
        return
    end subroutine ppma_test01
    subroutine ppma_test02 ( )
        integer, allocatable :: r(:,:), g(:,:), b(:,:)
        character(len=80) :: file_name = 'test02.ascii.ppm'
        integer   :: file_unit
        integer   :: i, j, k, ierror, ios, ncol, nrow
        integer   :: rgb_max
        write(*,'(a)') ' '
        write(*,'(a)') 'TEST02'
        write(*,'(a)') '  PPMA_READ_HEADER reads the header.'
        write(*,'(a)') '  PPMA_READ_HEADER reads the data.'

        call ppma_write_test ( file_name )
        write(*,'(a)') ' '
        write(*,'(a)') '  PPMA_WRITE_TEST created some data.'
        call fopen ( file_unit, file = file_name, status = 'old', iostat = ios )

        if ( ios /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'TEST02 - Fatal error!'
            write(*,'(a)') '  Could not open the file.'
            return
        end if

        call ppma_read_header ( file_unit, nrow, ncol, rgb_max, ierror )

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'TEST02 - Fatal error!'
            write(*,'(a)') '  Error while reading the header.'
            return
        end if

        write(*,'(a)') ' '
        write(*,'(a)') '  PPMA_READ_HEADER read the header.'
        write(*,'(a)') ' '
        write ( *, '(a,i6)' ) '  Number of rows of data =    ', nrow
        write ( *, '(a,i6)' ) '  Number of columns of data = ', ncol
        write ( *, '(a,i6)' ) '  Maximum RGB value =         ', rgb_max

        allocate ( r(nrow,ncol), g(nrow,ncol), b(nrow,ncol) )

        call ppma_read_data ( file_unit, nrow, ncol, r, g, b, ierror )

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'TEST02 - Fatal error!'
            write(*,'(a)') '  PPMA_READ_DATA failed.'
            return
        end if

        write(*,'(a)') ' '
        write(*,'(a)') '  PPMA_READ_DATA read the data.'
        write(*,'(a)') ' '
        write(*,'(a)') '  Ten sample values:'
        write(*,'(a)') ' '
        do k = 1, 10
            i = ( ( 10 - k ) * 1 + ( k - 1 ) * nrow ) / ( 10 - 1 )
            j = ( ( 10 - k ) * 1 + ( k - 1 ) * ncol ) / ( 10 - 1 )
            write ( *, '(2i4,4x,3i6)' ) i, j, r(i,j), g(i,j), b(i,j)
        end do

        call ppma_check_data ( nrow, ncol, rgb_max, r, g, b, ierror )

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'TEST02 - Error!'
            write ( *, '(a,i6)' ) '  The data was not accepted by PPMA_CHECK_DATA.'
            return
        end if

        write(*,'(a)') ' '
        write(*,'(a)') '  The data was accepted by PPMA_CHECK_DATA.'

        deallocate ( r ,g, b )
        call fclose( file_unit )
    end subroutine ppma_test02

    subroutine pgma_test01 ( )
        integer , parameter :: ncol = 300
        integer , parameter :: nrow = 300

        character ( len = 80 ) :: file_name = 'pgma_io_prb_01.ascii.pgm'
        integer  g(nrow,ncol)
        integer  ierror

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01'
        write ( *, '(a)' ) '  PGMA_EXAMPLE sets up ASCII PGM data.'
        write ( *, '(a)' ) '  PGMA_WRITE writes an ASCII PGM file.'

        call pgma_example ( nrow, ncol, g )

        call pgma_write ( file_name, nrow, ncol, g, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8)' ) 'PGMA_WRITE returns IERROR = ', ierror
            return
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Wrote the header and data for "' &
            // trim ( file_name ) //'".'
        write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
        write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol

        return
    end subroutine pgma_test01
    subroutine pgma_test02 ( )
        character ( len = 80 ) :: file_name = 'pgma_io_prb_02.ascii.pgm'
        integer file_unit
        integer, allocatable, dimension ( :, : ) :: g
        integer :: i, ierror, ios, j, k, maxg, ncol, nrow

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02'
        write ( *, '(a)' ) '  PGMA_READ reads an ASCII PGM file.'

        call pgma_write_test ( file_name )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  PGMA_WRITE_TEST created some data.'

        call fopen ( file_unit, file = file_name, status = 'old', iostat = ios )

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TEST02 - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file.'
            return
        end if

        call pgma_read_header ( file_unit, nrow, ncol, maxg )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  PGMA_READ_HEADER read the header.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
        write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
        write ( *, '(a,i8)' ) '  Maximum G value =           ', maxg

        allocate ( g(nrow,ncol) )

        call pgma_read_data ( file_unit, nrow, ncol, g )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  PGMA_READ_DATA read the data.'

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Sample data:'
        write ( *, '(a)' ) ' '

        do k = 1, 10
            i = ( ( 10 - k ) * 1 + ( k - 1 ) * nrow ) / ( 10 - 1 )
            j = ( ( 10 - k ) * 1 + ( k - 1 ) * ncol ) / ( 10 - 1 )
            write ( *, '(i4,2x,i4,2x,i6)' ) i, j, g(i,j)
        end do

        call pgma_check_data ( nrow, ncol, maxg, g, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'TEST02'
            write ( *, '(a,i8)' ) '  The data was not accepted by PGMA_CHECK_DATA.'
            return
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The data was accepted by PGMA_CHECK_DATA.'

        deallocate ( g )

        call fclose( file_unit )
    end subroutine pgma_test02

    subroutine pgma_test03 ( )

        integer , parameter :: ncol = 300
        integer , parameter :: ngray = 11
        integer , parameter :: nrow = 300

        character ( len = 80 ) :: file_name = 'pgma_io_prb_03.ascii.pgm'
        integer  g(nrow,ncol)
        real ( kind = 8 ), dimension ( ngray ) :: gray = (/ &
            0.000D+00, 0.291D+00, 0.434D+00, 0.540D+00, 0.629D+00, &
            0.706D+00, 0.774D+00, 0.837D+00, 0.895D+00, 0.949D+00, &
            1.000D+00 /)
        integer  i
        integer  ierror
        integer  j
        integer  k

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST03'
        write ( *, '(a)' ) '  PGMA_WRITE writes an ASCII PGM file.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  In this example, we make a sort of grayscale'
        write ( *, '(a)' ) '  checkerboard.'

        do i = 1, nrow
            do j = 1, ncol
                k = ( i - 1 + j - 1 ) * ngray / min ( nrow, ncol )
                k = 1 + mod ( k, ngray )
                g(i,j) = int ( 255.0D+00 * gray(k) )
            end do
        end do

        call pgma_write ( file_name, nrow, ncol, g, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a,i8)' ) 'PGMA_WRITE returns IERROR = ', ierror
        else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  Wrote the header and data for "' &
                // trim ( file_name ) //'".'
            write ( *, '(a,i8)' ) '  Number of rows of data =    ', nrow
            write ( *, '(a,i8)' ) '  Number of columns of data = ', ncol
        end if

        return
    end subroutine pgma_test03
    subroutine pbma_example ( row_num, col_num, b )
        integer, intent(in) :: col_num, row_num
        integer, intent(inout) :: b(row_num,col_num)
        integer :: i, j
        real (kind=dp) :: r, test, x, xc, y, yc

        xc = real ( col_num, kind = 8 ) / 2.0D+00
        yc = real ( row_num, kind = 8 ) / 2.0D+00
        r = real ( min ( row_num, col_num ), kind = 8 ) / 3.0D+00

        do i = 1, row_num
            y = real ( i, kind = 8 )
            do j = 1, col_num
                x = real ( j, kind = 8 )
                test = r - sqrt ( ( x - xc )**2 + 0.75D+00 * ( y - yc )**2 )
                if ( abs ( test ) <= 3.0D+00 ) then
                    b(i,j) = 1
                else
                    b(i,j) = 0
                end if
            end do
        end do
    end subroutine pbma_example

    subroutine pbma_read_test ( file_in_name )
        character(len=*), intent(in) :: file_in_name
        integer, allocatable, dimension ( :, : ) :: b
        integer :: file_in_unit, ierror, ios
        integer :: col_num, row_num

        call fopen ( file_in_unit, file = file_in_name, status = 'old', &
            iostat = ios )

        if ( ios /= 0 ) then
            ierror = 1
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_READ_TEST - Fatal error!'
            write(*,'(a)') '  Could not open the file.'
            stop
        end if
        !  Read the header.
        call pbma_read_header ( file_in_unit, row_num, col_num )
        !  Allocate the data.
        allocate ( b(row_num,col_num) )
        !  Read the data.
        call pbma_read_data ( file_in_unit, row_num, col_num, b )

        call fclose( file_in_unit )
        !  Check the data.
        call pbma_check_data ( row_num, col_num, b )

        write(*,'(a)') ' '
        write(*,'(a)') 'PBMA_READ_TEST:'
        write(*,'(a)') '  PBMA_CHECK_DATA has approved the data from the file.'

        deallocate ( b )

        return
    end subroutine pbma_read_test

    subroutine pbma_write_test ( file_out_name )
        integer, allocatable, dimension ( :, : ) :: b
        character(len=*) file_out_name
        integer :: col_num
        integer :: row_num

        row_num = 200
        col_num = 200
        !  Allocate memory.
        allocate ( b(row_num,col_num) )
        !  Set the data.
        call pbma_example ( row_num, col_num, b )
        !  Write the data to the file.
        call pbma_write ( file_out_name, row_num, col_num, b )

        deallocate ( b );

        return
    end subroutine pbma_write_test

    subroutine ppma_read_test ( file_in_name, ierror )
        integer, allocatable, dimension ( :, : ) :: b
        character(len=*), intent(in) :: file_in_name
        integer :: file_in_unit
        integer, allocatable, dimension ( :, : ) :: g
        integer, intent(out) :: ierror
        integer :: ios
        integer :: col_num
        integer :: row_num
        integer, allocatable, dimension ( :, : ) :: r
        integer :: rgb_max

        call fopen (file_in_unit, file = file_in_name, status = 'old', &
            iostat = ios )

        if ( ios /= 0 ) then
            ierror = 1
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_TEST - Fatal error!'
            write(*,'(a)') '  Could not open the file.'
            return
        end if
        !  Read the header.
        call ppma_read_header ( file_in_unit, row_num, col_num, rgb_max, ierror )

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_TEST - Fatal error!'
            write(*,'(a)') '  PPMA_READ_HEADER failed.'
            return
        end if
        !  Allocate the data.
        allocate ( r(row_num,col_num) )
        allocate ( g(row_num,col_num) )
        allocate ( b(row_num,col_num) )
        !  Read the data.
        call ppma_read_data ( file_in_unit, row_num, col_num, r, g, b, ierror )

        if ( ierror /= 0 ) then

            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_TEST - Fatal error!'
            write(*,'(a)') '  PPMA_READ_HEADER failed.'
            deallocate ( r, g, b )
            return
        end if

        call fclose( file_in_unit )
        !  Check the data.
        call ppma_check_data ( row_num, col_num, rgb_max, r, g, b, ierror )

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_TEST - Fatal error!'
            write(*,'(a)') '  PPMA_CHECK_DATA did not approve the data.'
        else
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_TEST:'
            write(*,'(a)') '  PPMA_CHECK_DATA has approved the data from the file.'
        end if
        deallocate (r , b, g)

    end subroutine ppma_read_test
    subroutine ppma_example ( row_num, col_num, r, g, b )
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(inout) :: b(row_num,col_num)
        integer, intent(inout) :: g(row_num,col_num)
        integer, intent(inout) :: r(row_num,col_num)
        real (kind=dp) :: f1, f2, f3
        integer :: i, j
        real (kind=dp) :: x, y

        do i = 1, row_num
            y = real ( row_num - i, kind = 8 ) / real ( row_num - 1, kind = 8 )
            do j = 1, col_num
                x = real ( j - 1, kind = 8 ) / real ( col_num - 1, kind = 8 )
                f1 = 4.0D+00 * ( x - 0.5D+00 )**2
                f2 = sin ( 3.14159265D+00 * x )
                f3 = x
                if ( y <= f1 ) then
                    r(i,j) = int ( 255.0D+00 * f1 )
                else
                    r(i,j) = 50
                end if
                if ( y <= f2 ) then
                    g(i,j) = int ( 255.0D+00 * f2 )
                else
                    g(i,j) = 150
                end if
                if ( y <= f3 ) then
                    b(i,j) = int ( 255.0D+00 * f3 )
                else
                    b(i,j) = 250
                end if
            end do
        end do
    end subroutine ppma_example

    subroutine ppma_write_test ( file_out_name )
        integer, allocatable, dimension ( :, : ) :: b
        character(len=*), intent(in) :: file_out_name
        integer, allocatable, dimension ( :, : ) :: g
        integer :: ierror
        integer :: col_num
        integer :: row_num
        integer, allocatable, dimension ( :, : ) :: r

        row_num = 300
        col_num = 300
        !  Allocate memory.
        allocate ( r(row_num,col_num) )
        allocate ( g(row_num,col_num) )
        allocate ( b(row_num,col_num) )
        !  Set the data.
        call ppma_example ( row_num, col_num, r, g, b )
        !  Write the data to the file.
        call ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )

        deallocate ( r );
        deallocate ( g );
        deallocate ( b );

        if ( ierror /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_WRITE_TEST - Fatal error!'
            write(*,'(a)') '  PPMA_WRITE failed.'
        end if

        return
    end subroutine ppma_write_test

    subroutine pgma_example ( row_num, col_num, g )
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(inout) :: g(row_num,col_num)
        integer :: i, j
        integer, parameter :: periods = 3
        real ( kind = 8 ), parameter :: pi = 3.14159265D+00
        real ( kind = 8 ) x
        real ( kind = 8 ) y

        do i = 1, row_num
            y = 2.0D+00 * real ( i - 1, kind = 8 ) &
                / real ( row_num - 1, kind = 8 ) - 1.0D+00
            do j = 1, col_num
                x = 2.0D+00 * pi * real ( periods * ( j - 1 ), kind = 8 ) &
                    / real ( col_num - 1, kind = 8 )
                g(i,j) = int ( 20.0D+00 * ( sin ( x ) - y + 2.0D+00 ) )
            end do
        end do

        return
    end subroutine pgma_example

    !! PGMA_WRITE_TEST tests the ASCII PGM write routines.
    subroutine pgma_write_test ( file_out_name )
        character(len=*), intent(in) :: file_out_name
        integer, allocatable, dimension ( :, : ) :: g
        integer ierror
        integer col_num
        integer row_num

        row_num = 300
        col_num = 300
        !  Allocate memory.
        allocate ( g(row_num,col_num) )
        !  Set the data.
        call pgma_example ( row_num, col_num, g )
        !  Write the data to the file.
        call pgma_write ( file_out_name, row_num, col_num, g, ierror )

        deallocate ( g );

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_WRITE_TEST - Fatal error!'
            write ( *, '(a)' ) '  PGMA_WRITE failed.'
            stop
        end if
    end subroutine pgma_write_test

    subroutine pgma_read_test ( file_in_name, ierror )
        character(len=*), intent(in) :: file_in_name
        integer, intent(out) :: ierror
        integer, allocatable :: g(:,:)
        integer :: file_in_unit, ios, g_max, col_num, row_num

        call fopen ( file_in_unit, file = file_in_name, status = 'old', &
            iostat = ios )
        if ( ios /= 0 ) then
            ierror = 1
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_TEST - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file.'
            stop
        end if
        !  Read the header.
        call pgma_read_header ( file_in_unit, row_num, col_num, g_max )
        !  Allocate the data.
        allocate ( g(row_num,col_num) )
        !  Read the data.
        call pgma_read_data ( file_in_unit, row_num, col_num, g )
        call fclose (file_in_unit )
        !  Check the data.
        call pgma_check_data ( row_num, col_num, g_max, g, ierror )

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_TEST - Warning!'
            write ( *, '(a)' ) '  PGMA_CHECK_DATA did not approve the data.'
        else
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_TEST:'
            write ( *, '(a)' ) '  PGMA_CHECK_DATA has approved the data from the file.'
        end if
        deallocate ( g )
    end subroutine pgma_read_test



end module simple_test_pnm_io
