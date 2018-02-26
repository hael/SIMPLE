
!! Portable aNy Map image formats  -- Netpbm formats pbm, ppm, bmp

!! Based on John Burkardt's  PPMA and PBMA modules (LPGL Copyright, March 2003)
!! Modified by Michael Eager 2018, Monash University
module simple_pnm
    include 'simple_lib.f08'
    use, intrinsic :: iso_c_binding
    implicit none


!     public :: pbma_read_data, pbma_check_data, pbma_read_header
!     public :: pbma_write, pbma_write_data, pbma_write_header

!     public :: ppma_read_data, ppma_read_header, ppma_check_data
!     public :: ppma_write, ppma_write_data, ppma_write_header

!     public :: pgma_check_data, pgma_read_data, pgma_read_header
!     public :: pgma_write, pgma_write_data, pgma_write_header
!     private

#include "simple_local_flags.inc"




    !      type PNM
    !          type(Pixel), allocatable :: buffer(:)
    ! !         type(Pixel), allocatable :: ext_buffer(:)
    !          integer :: width ,height
    !      contains
    !          procedure :: read
    !          procedure :: write
    !          procedure :: new
    !          procedure :: getHeight
    !          procedure :: getWidth
    !          procedure :: getPixelAt
    !          procedure :: setPixelAt_Int
    !          procedure :: setPixelAt_Real
    !          procedure :: kill
    !      end type PNM



    ! !     !     !
    ! !     !     !  KRONROD is provided by the C library, and so the following
    ! !     !     !  INTERFACE block must be set up to describe how data is to
    ! !     !     !  be passed.
    ! !     !     !
    ! !     interface
    ! !         !         subroutine read_pnm (  ) bind ( c )
    ! !         !             use iso_c_binding
    ! !         !             integer ( c_int ), VALUE :: n
    ! !         !             real ( c_double ), VALUE :: eps
    ! !         !             real ( c_double ) :: x(*)
    ! !         !             real ( c_double ) :: w1(*)
    ! !         !             real ( c_double ) :: w2(*)
    ! !         !         end subroutine  read_pnm
    ! !     end interface

  contains


    !         integer function getWidth(self)
    !             class(PNM), intent(inout) :: self
    !             getWidth = self%width
    !         end function getWidth
    !         integer function getHeight(self)
    !             class(PNM), intent(inout) :: self
    !             getHeight=self%height
    !         end function getHeight
    !         function buffer() result(data)
    !             class(PNM) intent(inout) :: self
    !             integer,allocatable :: data(:)
    !         end function buffer

    !     function getPixelAt(self, x, y) result(px)
    !         class(PNM), intent(inout) :: self
    !          integer, intent(in) :: x,y
    !         type(Pixel) :: px
    !         px = self%buffer(x + self%width*y )
    !     end function getPixelAt

    !     subroutine setPixelAt_Real(self, x, y, val)
    !         class(PNM), intent(inout) :: self
    !         integer, intent(in) :: x,y
    !         real, intent(inout) :: val
    !         type(Pixel) :: px
    !         integer :: rawval
    !         rawval = px%set_grey(val)
    !         self%buffer(x + self%width*y )= px
    !     end subroutine setPixelAt_Real

    !     subroutine setPixelAt_Int(self, x, y, val)
    !         class(PNM), intent(inout) :: self
    !         integer, intent(in) :: x,y
    !         integer, intent(inout) :: val
    !         type(Pixel) :: px
    !         integer :: rawval
    !         rawval = px%set_blue(val)
    !         rawval = px%set_red(val)
    !         rawval = px%set_green(val)
    !         self%buffer(x + self%width*y )= px
    !     end subroutine setPixelAt_Int

    !      subroutine new (self, width,  height, buf, file_name)
    !          class(PNM), intent(inout):: self
    !          integer, intent(in) :: width, height
    !          real, intent(inout),optional :: buf(:,:)
    !          character(len=*), intent(inout),optional:: file_name
    !          integer :: i,j
    !          self%width = width
    !          self%height = height
    !         if (allocated(self%buffer)) deallocate(self%buffer)
    !         if(present(buf))then
    !             if( size(buf,1)*size(buf,2) /= width * height) then
    !                 stop 'simple_pnm_io:: new input buffer not the same size as width|height inputs '
    !             end if
    !             allocate( self%buffer(width*height))
    !             do i=1, self%width
    !                 do j=1, self%height
    !                     call self%setPixelAt_Real(i,j,buf(i,j))
    !                 end do
    !             end do

    !         else
    !             allocate( self%buffer(width*height) )
    !         end if
    !      end subroutine new


    !     subroutine kill (self)
    !         class(PNM),intent(inout) :: self

    !         deallocate(self%buffer)
    !     end subroutine kill

subroutine getint ( done, ierror, inunit, ival, string )

!*****************************************************************************80
!
!! GETINT reads an integer from a file.
!
!  Discussion:
!
!    The file, or at least the part read by GETINT, is assumed to
!    contain nothing but integers.  These integers may be separated
!    by spaces, or appear on separate lines.  Comments, which begin
!    with "#" and extend to the end of the line, may appear anywhere.
!
!    Each time GETINT is called, it tries to read the next integer
!    it can find.  It remembers where it was in the current line
!    of text.
!
!    The user should open a text file on FORTRAN unit INUNIT,
!    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
!    care of reading in a new STRING as necessary, and extracting
!    as many integers as possible from the line of text before
!    reading in the next line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 October 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, logical DONE.
!
!    On input, if this is the first call, or the user has changed
!    STRING, then set DONE = TRUE.
!
!    On output, if there is no more data to be read from STRING,
!    then DONE is TRUE.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred while trying to read the integer.
!
!    Input, integer ( kind = 4 ) INUNIT, the FORTRAN unit from which to read.
!
!    Output, integer ( kind = 4 ) IVAL, the integer that was read.
!
!    Input/output, character ( len = * ) STRING, the text of the most recently
!    read line of the file.
!
  logical done
  integer  i
  integer ierror
  integer inunit
  integer ios
  integer ival
  integer last
  character ( len = * ) string
  character ( len = 80 ) word

  do

    call word_next_rd ( string, word, done )

    if ( .not. done ) then
      exit
    end if

    read ( inunit, '(a)', iostat = ios ) string

    if ( ios /= 0 ) then
      ierror = 1
      return
    end if

    i = index ( string, '#' )
    if ( i /= 0 ) then
      string(i:) = ' '
    end if

  end do

  call str2int ( word, ierror, ival )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GETINT - Fatal error!'
    write ( *, '(a)' ) '  Error trying to convert string to integer.'
    stop
  end if

  return
end subroutine getint

subroutine word_next_rd ( line, word, done )

!*****************************************************************************80
!
!! WORD_NEXT_RD "reads" words from a string, one at a time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) LINE, a string, presumably containing
!    words separated by spaces.
!
!    Output, character ( len = * ) WORD.
!    If DONE is FALSE,
!      WORD contains the "next" word read from LINE.
!    Else
!      WORD is blank.
!
!    Input/output, logical DONE.
!    On input, on the first call, or with a fresh value of LINE,
!      set DONE to TRUE.
!    Else
!      leave it at the output value of the previous call.
!    On output, if a new nonblank word was extracted from LINE
!      DONE is FALSE
!    ELSE
!      DONE is TRUE.
!    If DONE is TRUE, then you need to provide a new LINE of data.
!
!  Local Parameters:
!
!    NEXT is the next location in LINE that should be searched.
!
  implicit none

  logical done
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) lenl
  character ( len = * ) line
  integer ( kind = 4 ), save :: next = 1
  character ( len = 1 ), parameter :: TAB = char(9)
  character ( len = * ) word

  lenl = len_trim ( line )

  if ( done ) then
    next = 1
    done = .false.
  end if
!
!  Beginning at index NEXT, search LINE for the next nonblank.
!
  ilo = next

  do
!
!  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
!
    if ( lenl < ilo ) then
      word = ' '
      done = .true.
      next = lenl + 1
      return
    end if
!
!  ...If the current character is blank, skip to the next one.
!
    if ( line(ilo:ilo) /= ' ' .and. line(ilo:ilo) /= TAB ) then
      exit
    end if

    ilo = ilo + 1

  end do
!
!  To get here, ILO must be the index of the nonblank starting
!  character of the next word.
!
!  Now search for the LAST nonblank character.
!
  next = ilo + 1

  do

    if ( lenl < next ) then
      word = line(ilo:next-1)
      return
    end if

    if ( line(next:next) == ' ' .or. line(next:next) == TAB ) then
      exit
    end if

    next = next + 1

  end do

  word = line(ilo:next-1)

  return
end subroutine word_next_rd

    subroutine pbma_read_data ( file_in_unit, row_num, col_num, b )
        integer, intent(in) :: col_num, row_num
        integer, intent(in) :: file_in_unit
        integer, intent(inout) :: b(row_num,col_num)
        character  c
        integer :: i, ierror, ios, j, k, k_max
        character ( len = 80 ) :: string

        ierror = 0
        k = 0
        k_max = 0
        string = ' '

        do i = 1, row_num
            do j = 1, col_num
                do
                    if ( k_max <= k ) then
                        read ( file_in_unit, '(a)', iostat = ios ) string
                        if ( ios /= 0 ) then
                            ierror = 1
                            write(*,'(a)') ' '
                            write(*,'(a)') 'PBMA_READ_DATA - Fatal error!'
                            write(*,'(a)') '  Problem reading data.'
                            stop
                        end if

                        k = 0
                        k_max = len_trim ( string )
                        if ( k_max <= 0 ) then
                            cycle
                        end if
                    end if
                    k = k + 1
                    c = string(k:k)
                    if ( c == '0' ) then
                        b(i,j) = 0
                        exit
                    else if ( c == '1' ) then
                        b(i,j) = 1
                        exit
                    end if
                end do
            end do
        end do
    end subroutine pbma_read_data

    subroutine pbma_read_header ( file_in_unit, row_num, col_num )
        integer :: file_in_unit
        integer :: ierror
        integer :: ios
        character ( len = 2 ) magic
        integer :: col_num
        integer :: row_num
        logical done
        character ( len = 255 ) :: string
        !  Read the first line of data, which must begin with the magic number.
        read ( file_in_unit, '(a)', iostat = ios ) magic

        if ( ios /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  End or error while reading file.'
            ierror = 2
            stop
        end if

        if ( .not. stringsAreEqual ( magic, 'P1' ) ) then
            ierror = 3
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_READ_HEADER - Fatal error.'
            write(*,'(a)') '  First two bytes are not magic number "P1".'
            write(*,'(a)') '  First two bytes are: "' // magic // '".'
            stop
        end if
        !  Now search for COL_NUM and ROW_NUM.
        done = .true.
        string = ' '

        call getint ( done, ierror, file_in_unit, col_num, string )

        if ( ierror /= 0 ) then
            call fclose( file_in_unit )
            ierror = 4
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  Problem reading COL_NUM.'
            stop
        end if

        call getint ( done, ierror, file_in_unit, row_num, string )

        if ( ierror /= 0 ) then
            ierror = 4
            call fclose( file_in_unit )
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  Problem reading ROW_NUM.'
            stop
        end if
    end subroutine pbma_read_header

    subroutine pbma_write ( file_out_name, row_num, col_num, b )

        !*****************************************************************************80
        !! PBMA_WRITE writes an ASCII PBM file.
        !  Example:
        !    P1
        !    # feep.pbma created by PBMA_IO(PBMA_WRITE).
        !    24 7
        !    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        !    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
        !    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
        !    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
        !    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
        !    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
        !    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        !  Licensing:
        !    This code is distributed under the GNU LGPL license.
        !  Modified:
        !    04 June 2010
        !  Author:
        !    John Burkardt
        !  Parameters:
        !    Input, character(len=*) FILE_OUT_NAME, the name of the file.
        !    Input, integer :: ROW_NUM, COL_NUM, the number of rows
        !    and columns of data.
        !    Input, integer :: B(ROW_NUM,COL_NUM), the bit value of each
        !    pixel.  These should be 0 or 1.
        implicit none

        integer :: col_num
        integer :: row_num

        integer :: b(row_num,col_num)
        character(len=*) file_out_name
        integer :: file_out_unit
        integer :: ierror
        integer :: ios

        ierror = 0
        !  Open the file.
        call fopen ( file_out_unit, file = file_out_name, status = 'replace', &
            form = 'formatted', access = 'sequential', iostat = ios )

        if ( ios /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_WRITE - Fatal error!'
            write(*,'(a)') '  Could not open the file.'
            stop
        end if
        !  Write the header.
        call pbma_write_header ( file_out_name, file_out_unit, row_num, col_num )
        !  Write the data.
        call pbma_write_data ( file_out_unit, row_num, col_num, b )
        !  Close the file.
        call fclose( file_out_unit )

        return
    end subroutine pbma_write
    subroutine pbma_write_data ( file_out_unit, row_num, col_num, b )
        integer, intent(in) :: col_num, row_num
        integer, intent(in) :: b(row_num,col_num)
        integer, intent(in) :: file_out_unit
        integer :: i
        integer :: jhi
        integer :: jlo
        !  Write the header.
        do i = 1, row_num
            do jlo = 1, col_num, 60
                jhi = min ( jlo + 59, col_num )
                write ( file_out_unit, '(60i1)' ) b(i,jlo:jhi)
            end do
        end do

        return
    end subroutine pbma_write_data

    subroutine pbma_write_header ( file_out_name, file_out_unit, row_num, col_num )
        character(len=*), intent(in):: file_out_name
        integer, intent(in) :: file_out_unit
        integer, intent(in) :: col_num, row_num
        character ( len = 2 ) :: magic = 'P1'

        !  Write the header.
        write ( file_out_unit, '(a2)' ) magic
        write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
            // ' created by PBMA_IO::PBMA_WRITE.F90.'
        write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num

        return
    end subroutine pbma_write_header

    subroutine pbma_check_data ( row_num, col_num, b )
        integer, intent(in) :: col_num, row_num
        integer, intent(in) :: b(row_num,col_num)
        integer :: ierror

        ierror = 0
        if ( minval ( b(1:row_num,1:col_num) ) < 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_CHECK_DATA - Fatal error!'
            write(*,'(a)') '  At least one bit value is below 0.'
            ierror = 1
            stop
        end if
        if ( 1 < maxval ( b(1:row_num,1:col_num) ) ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PBMA_CHECK_DATA - Fatal error!'
            write(*,'(a)') '  At least one bit value exceeds 1.'
            ierror = 1
            stop
        end if
    end subroutine pbma_check_data

    subroutine ppma_read_data ( file_in_unit, row_num, col_num, r, g, b, ierror )
        integer :: col_num
        integer :: row_num

        integer :: b(row_num,col_num)
        logical done
        integer :: file_in_unit
        integer :: g(row_num,col_num)
        integer :: i
        integer :: ierror
        integer :: j
        integer :: r(row_num,col_num)
        character ( len = 80 ) string

        ierror = 0
        done = .true.
        string = ' '

        do i = 1, row_num
            do j = 1, col_num

                call getint ( done, ierror, file_in_unit, r(i,j), string )

                if ( ierror /= 0 ) then
                    ierror = 5
                    call fclose( file_in_unit )
                    write(*,'(a)') ' '
                    write(*,'(a)') 'PPMA_READ_DATA - Fatal error!'
                    write(*,'(a)') '  Problem reading R data.'
                    return
                end if

                call getint ( done, ierror, file_in_unit, g(i,j), string )

                if ( ierror /= 0 ) then
                    ierror = 5
                    call fclose( file_in_unit )
                    write(*,'(a)') ' '
                    write(*,'(a)') 'PPMA_READ_DATA - Fatal error!'
                    write(*,'(a)') '  Problem reading G data.'
                    return
                end if

                call getint ( done, ierror, file_in_unit, b(i,j), string )

                if ( ierror /= 0 ) then
                    ierror = 5
                    call fclose( file_in_unit )
                    write(*,'(a)') ' '
                    write(*,'(a)') 'PPMA_READ_DATA - Fatal error!'
                    write(*,'(a)') '  Problem reading B data.'
                    return
                end if

            end do
        end do
    end subroutine ppma_read_data
    subroutine ppma_read_header ( file_in_unit, row_num, col_num, rgb_max, ierror )
        integer, intent(in) :: file_in_unit
        integer, intent(out) :: ierror
        integer, intent(out) :: col_num
        integer, intent(out) :: row_num
        integer, intent(out) :: rgb_max
        integer :: ios
        character ( len = 2 ) magic
        logical  done
        character ( len = 80 ) string
        !  Read the first line of data, which must begin with the magic number.
        read ( file_in_unit, '(a)', iostat = ios ) magic

        if ( ios /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  End or error while reading file.'
            ierror = 2
            return
        end if

        if ( .not. stringsAreEqual ( magic, 'P3' ) ) then
            ierror = 3
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_HEADER - Fatal error.'
            write(*,'(a)') '  First two bytes are not magic number "P3".'
            write(*,'(a)') '  First two bytes are: "' // magic // '".'
            return
        end if
        !  Now search for COL_NUM, ROW_NUM, and RGB_MAX.
        done = .true.
        string = ' '

        call getint ( done, ierror, file_in_unit, col_num, string )

        if ( ierror /= 0 ) then
            call fclose( file_in_unit )
            ierror = 4
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  Problem reading COL_NUM.'
            return
        end if

        call getint ( done, ierror, file_in_unit, row_num, string )

        if ( ierror /= 0 ) then
            ierror = 4
            call fclose( file_in_unit )
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  Problem reading ROW_NUM.'
            return
        end if

        call getint ( done, ierror, file_in_unit, rgb_max, string )

        if ( ierror /= 0 ) then
            ierror = 4
            call fclose( file_in_unit )
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_READ_HEADER - Fatal error!'
            write(*,'(a)') '  Problem reading RGB_MAX.'
            return
        end if

        return
    end subroutine ppma_read_header

    subroutine ppma_check_data ( row_num, col_num, rgb_max, r, g, b, ierror )
        integer, intent(out) :: ierror
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: rgb_max
        integer, intent(in) ::  b(row_num,col_num)
        integer, intent(in)  :: g(row_num,col_num)
        integer, intent(in)  :: r(row_num,col_num)

        ierror = 0
        !  Make sure no color is negative.
        if ( minval ( r(1:row_num,1:col_num) ) < 0 .or. &
            minval ( g(1:row_num,1:col_num) ) < 0 .or. &
            minval ( b(1:row_num,1:col_num) ) < 0 ) then
            ierror = 1
            return
        end if
        !  Make sure no color is greater than RGB_MAX.
        if ( rgb_max < maxval ( r(1:row_num,1:col_num) ) .or. &
            rgb_max < maxval ( g(1:row_num,1:col_num) ) .or. &
            rgb_max < maxval ( b(1:row_num,1:col_num) ) ) then
            ierror = 1
            return
        end if

    end subroutine ppma_check_data

    subroutine ppma_write ( file_out_name, row_num, col_num, r, g, b, ierror )
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: b(row_num,col_num)
        integer, intent(in) :: g(row_num,col_num)
        integer, intent(in) :: r(row_num,col_num)
        character (len=*), intent(in) :: file_out_name
        integer, intent(inout) :: ierror
        integer :: ios, file_out_unit, rgb_max

        ierror = 0
        !  Compute the maximum color value.
        rgb_max = max ( &
            maxval ( r(1:row_num,1:col_num) ), &
            maxval ( g(1:row_num,1:col_num) ), &
            maxval ( b(1:row_num,1:col_num) ) )
        !  Open the file.


        call fopen ( file_out_unit, file = file_out_name, status = 'replace', &
            form = 'formatted', access = 'sequential', iostat = ios )
        if ( ios /= 0 ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_WRITE - Fatal error!'
            write(*,'(a)') '  Could not open the file.'
            ierror = 2
            return
        end if
        !  Write the header.
        call ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
            rgb_max, ierror )
        !  Write the data.
        call ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )
        !  Close the file.
        call fclose( file_out_unit )
        !  Report
        if ( debug ) then
            write(*,'(a)') ' '
            write(*,'(a)') 'PPMA_WRITE - Note:'
            write(*,'(a)') '  The data was checked and written.'
            write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
            write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
            write ( *, '(a,i8)' ) '  Maximum RGB value RGB_MAX =      ', rgb_max
        end if

    end subroutine ppma_write
    subroutine ppma_write_data ( file_out_unit, row_num, col_num, r, g, b, ierror )
        integer, intent(in) :: file_out_unit
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: b(row_num,col_num)
        integer, intent(in) :: g(row_num,col_num)
        integer, intent(in) :: r(row_num,col_num)
        integer, intent(inout), optional :: ierror
        integer :: i,j, jhi, jlo

        ierror = 0
        !  Write the header.
        do i = 1, row_num
            do jlo = 1, col_num, 4
                jhi = min ( jlo + 3, col_num )
                write ( file_out_unit, '(12i5)' ) ( r(i,j), g(i,j), b(i,j), j = jlo, jhi )
            end do
        end do

        return
    end subroutine ppma_write_data
    subroutine ppma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
        rgb_max, ierror )
        character(len=*), intent(in) :: file_out_name
        integer, intent(in)  :: file_out_unit
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: rgb_max
        integer, intent(inout), optional :: ierror

        character ( len = 2 ) :: magic = 'P3'
        ierror = 0
        !  Write the header.
        write ( file_out_unit, '(a2)' ) magic
        write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
            // ' created by PPMA_IO::PPMA_WRITE.F90.'
        write ( file_out_unit, '(i5,2x,i5)' ) col_num, row_num
        write ( file_out_unit, '(i5)' ) rgb_max

        return
    end subroutine ppma_write_header


    subroutine pgma_check_data ( row_num, col_num, g_max, g, ierror )
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: g(row_num,col_num)
        integer, intent(out) :: ierror
        integer, intent(in) :: g_max

        ierror = 0

        if ( minval ( g(1:row_num,1:col_num) ) < 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
            write ( *, '(a)' ) '  At least one gray value is below 0.'
            ierror = 1
            stop
        end if

        if ( g_max < maxval ( g(1:row_num,1:col_num) ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_CHECK_DATA - Fatal error!'
            write ( *, '(a)' ) '  At least one gray value exceeds G_MAX.'
            write ( *, '(a,i12)' ) '  G_MAX = ', g_max
            ierror = 1
            stop
        end if
    end subroutine pgma_check_data


    subroutine pgma_read_data ( file_in_unit, row_num, col_num, g )
        integer, intent(in) :: col_num, row_num
        integer, intent(in) :: file_in_unit
        integer, intent(inout) :: g(row_num,col_num)

        logical done
        integer  :: i, ierror, j
        character(len=80) :: string

        ierror = 0
        done = .true.
        string = ' '

        do i = 1, row_num
            do j = 1, col_num

                call getint ( done, ierror, file_in_unit, g(i,j), string )

                if ( ierror /= 0 ) then
                    call fclose( file_in_unit )
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'PGMA_READ_DATA - Fatal error!'
                    write ( *, '(a)' ) '  Problem reading G data.'
                    stop
                end if

            end do
        end do
    end subroutine pgma_read_data

    subroutine pgma_read_header( file_in_unit, row_num, col_num, g_max )
        integer, intent(in) :: file_in_unit
        integer, intent(out) :: col_num
        integer, intent(out) :: row_num
        integer, intent(out) :: g_max
        integer :: ios, ierror
        character ( len = 2 ) magic
        logical  done
        character ( len = 80 ) string
        !  Read the first line of data, which must begin with the magic number.
        read ( file_in_unit, '(a)', iostat = ios ) magic

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
            write ( *, '(a)' ) '  End or error while reading file.'
            ierror = 2
            stop
        end if

        if ( .not. stringsAreEqual ( magic, 'P2' ) ) then
            ierror = 3
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error.'
            write ( *, '(a)' ) '  First two bytes are not magic number "P2".'
            write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
            stop
        end if
        !  Now search for COL_NUM, ROW_NUM, and G_MAX.
        done = .true.
        string = ' '

        call getint ( done, ierror, file_in_unit, col_num, string )

        if ( ierror /= 0 ) then
            call fclose( file_in_unit )
            ierror = 4
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
            write ( *, '(a)' ) '  Problem reading COL_NUM.'
            stop
        end if
        call getint ( done, ierror, file_in_unit, row_num, string )
        if ( ierror /= 0 ) then
            ierror = 4
            call fclose( file_in_unit )
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
            write ( *, '(a)' ) '  Problem reading ROW_NUM.'
            stop
        end if
        call getint ( done, ierror, file_in_unit, g_max, string )
        if ( ierror /= 0 ) then
            ierror = 4
            call fclose( file_in_unit )
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
            write ( *, '(a)' ) '  Problem reading G_MAX.'
            stop
        end if

    end subroutine pgma_read_header

    subroutine pgma_write ( file_out_name, row_num, col_num, g, ierror )
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: g(row_num,col_num)
        character (len=*), intent(in) :: file_out_name
        integer, intent(inout) :: ierror
        integer :: ios, g_max, file_out_unit

        ierror = 0
        !  Compute the maximum color value.
        g_max = maxval ( g(1:row_num,1:col_num) )
        !  Open the file.

        call fopen ( file_out_unit, file = file_out_name, status = 'replace', &
            form = 'formatted', access = 'sequential', iostat = ios )

        if ( ios /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_WRITE - Fatal error!'
            write ( *, '(a)' ) '  Could not open the file.'
            ierror = 2
            stop
        end if
        !  Write the header.
        call pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
            g_max, ierror )
        !  Write the data.
        call pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )
        !  Close the file.
        call fclose( file_out_unit )
        !  Report
        if ( debug ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_WRITE - Note:'
            write ( *, '(a)' ) '  The data was checked and written.'
            write ( *, '(a,i8)' ) '  Number of data rows ROW_NUM =    ', row_num
            write ( *, '(a,i8)' ) '  Number of data columns COL_NUM = ', col_num
            write ( *, '(a,i8)' ) '  Maximum gray value G_MAX =       ', g_max
        end if
    end subroutine pgma_write

    subroutine pgma_write_data ( file_out_unit, row_num, col_num, g, ierror )
        integer, intent(in) :: file_out_unit
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num
        integer, intent(in) :: g(row_num,col_num)
        integer, intent(inout), optional :: ierror
        integer :: i,j, jhi, jlo

        ierror = 0
        do i = 1, row_num
            do jlo = 1, col_num, 12
                jhi = min ( jlo + 11, col_num )
                write ( file_out_unit, '(12i5)' ) g(i,jlo:jhi)
            end do
        end do
    end subroutine pgma_write_data

    subroutine pgma_write_header ( file_out_name, file_out_unit, row_num, col_num, &
        g_max, ierror )
        character(len=*), intent(in) ::  file_out_name
        integer, intent(in) :: file_out_unit
        integer, intent(in) :: col_num
        integer, intent(in) :: row_num,g_max
        integer, intent(inout) :: ierror
        character ( len = 2 ) :: magic = 'P2'

        ierror = 0
        !  Write the header.
        write ( file_out_unit, '(a2)' ) magic
        write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) &
            // ' created by PGMA_IO::PGMA_WRITE.'
        write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num
        write ( file_out_unit, '(i8)' ) g_max

    end subroutine pgma_write_header

end module simple_pnm





! module pnmpixel
!     implicit none
!     public :: Pixel, BytesPerPixel
!     integer, parameter :: BytesPerPixel=4, MaxBitDepth=256
!     type Pixel
!         private
!         integer(1) :: data(4)
!     contains
!         procedure :: get_red
!         procedure :: set_red
!         procedure :: get_green
!         procedure :: set_green
!         procedure :: get_blue
!         procedure :: set_blue
!         procedure :: get_grey
!         procedure :: set_grey
!         procedure :: get_alpha
!         procedure :: set_alpha
!         !    generic :: red => get_red, set_red
!     end type Pixel
! contains
!     integer function get_red (self)
!         class(Pixel), intent(inout) :: self
!         get_red=self%data(1)
!     end function get_red
!     integer function set_red (self,arg)
!         class(Pixel), intent(inout) :: self
!         integer, intent(inout) :: arg
!         if (arg < 0 ) arg = 0
!         if(arg >= MaxBitDepth) arg = MaxBitDepth-1
!         self%data(1) = arg
!         set_red=self%data(1)
!     end function set_red
!     integer function get_green (self)
!         class(Pixel), intent(inout) :: self
!         get_green=self%data(2)
!     end function get_green
!     integer function set_green (self,arg)
!         class(Pixel), intent(inout) :: self
!         integer, intent(inout) :: arg
!         if (arg < 0 ) arg = 0
!         if(arg >= MaxBitDepth) arg = MaxBitDepth-1
!         self%data(2) = arg
!         set_green=self%data(2)
!     end function set_green
!     integer function get_blue (self)
!         class(Pixel), intent(inout) :: self
!         get_blue=self%data(3)
!     end function get_blue
!     integer function set_blue (self,arg)
!         class(Pixel), intent(inout) :: self
!         integer, intent(inout) :: arg
!         if (arg < 0 ) arg = 0
!         if(arg >= MaxBitDepth) arg =MaxBitDepth-1
!         self%data(3) = arg
!         set_blue = self%data(3)
!     end function set_blue
!     real function get_grey (self)
!         class(Pixel), intent(inout) :: self
!         get_grey = REAL( SUM(self%data(1:3)) ) / (REAL(MaxBitDepth-1) * 3.0)
!     end function get_grey
!     integer function get_alpha (self)
!         class(Pixel), intent(inout) :: self
!         get_alpha=self%data(4)
!     end function get_alpha
!     integer function set_alpha (self,arg)
!         class(Pixel), intent(inout) :: self
!         real, intent(inout) :: arg
!         if(arg<1.0) arg= arg * MaxBitDepth
!         self%data(4) = INT(arg)
!         set_alpha=self%data(4)
!     end function set_alpha
!     integer function set_grey (self, arg)
!         class(Pixel), intent(inout) :: self
!         real, intent(inout) :: arg
!         if (arg < 0.0) arg = 0.
!         if (arg > 1.0) arg = 1.0
!         self%data(1:3)  = NINT(arg *(MaxBitDepth-1))
!         set_grey=self%data(1)
!     end function set_grey

! end module pnmpixel
