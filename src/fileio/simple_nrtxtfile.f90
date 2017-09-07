! deals with text files of numbers
module simple_nrtxtfile
use simple_defs
use simple_strings ! use all in there
use simple_fileio
implicit none

public :: nrtxtfile
private

integer, parameter :: OPEN_TO_READ  = 1
integer, parameter :: OPEN_TO_WRITE = 2

type :: nrtxtfile
    private
    integer               :: funit              !< FILE*
    character(len=STDLEN) :: fname              !< filename
    integer               :: recs_per_line = 0  !< recordings per line
    integer               :: ndatalines    = 0  !< # lines
    integer               :: access_type        !< read/write
contains
    procedure          :: new
    procedure          :: readNextDataLine
    procedure, private :: writeDataLineReal
    procedure, private :: writeDataLineInt
    generic            :: write => writeDataLineReal, writeDataLineInt
    procedure          :: writeCommentLine
    procedure          :: get_nrecs_per_line
    procedure          :: get_ndatalines
    procedure          :: kill
end type

contains

    subroutine new( self, fname, access_type, wanted_recs_per_line )
        class(nrtxtfile),  intent(inout) :: self
        character(len=*),  intent(in)    :: fname       !< filename
        integer,           intent(in)    :: access_type !< Either OPEN_TO_READ or OPEN_TO_WRITE
        integer, optional, intent(in)    :: wanted_recs_per_line
        character(len=LINE_MAX_LEN)      :: buffer      !< will hold a line from the file
        integer                          :: recs_on_curr_line
        integer                          :: tot_nr_of_recs, tmpunit
        integer                          :: ios         !< ios is negative if an end of record condition is encountered or if
                                                        !! an endfile condition was detected. It is positive if an error was
                                                        !! detected. ios is zero otherwise.
        character(len=512)               :: io_msg
        ! always destruct first
        call self%kill
        tot_nr_of_recs   = 0;  recs_on_curr_line = 0
        self%fname       = fname
        self%access_type = access_type
        ! if we are opening to read, work out the file details..
        if( self%access_type .eq. OPEN_TO_READ )then
            if(.not.fopen(tmpunit, file=self%fname, iostat=ios, status='old'))&
                 call fileio_errmsg("simple_nrtxtfile::new; Error when opening file for reading ",ios)
            if( ios .ne. 0 )then
                if(.not.fclose(self%funit,ios)) call fileio_errmsg("nrtxt close failed ",ios)
                write(*,'(a)') 'simple_nrtxtfile::new; Error when opening file for reading: ', trim(self%fname)
                stop
            endif
            self%funit = tmpunit
            do
                ! work out records per line, and number_of_lines
                read(self%funit, '(a)', iostat=ios) buffer
                if ( ios == 0 ) then
                    if( strIsComment(buffer) .or. strIsBlank(buffer) )then
                        ! don't do anything
                    else
                        self%ndatalines = self%ndatalines+1
                        ! work out how many records on this line..
                        recs_on_curr_line = cntRecsPerLine(buffer)
                        tot_nr_of_recs = tot_nr_of_recs + recs_on_curr_line
                    endif
                else
                    exit ! we found the end...
                endif
            enddo
            ! if this file is correct, then the value left in recs_on_curr_line
            ! should be the same as the tot_nr_of_recs / number_of_lines.
            ! This is possibly a bit of a nasty way to check things, and may need to
            ! be changed.
            if( recs_on_curr_line * self%ndatalines .eq. tot_nr_of_recs )then
                ! All seems to be ok..
                self%recs_per_line = recs_on_curr_line
                rewind(self%funit)
            else
                ! Something went wrong..
                stop 'nrtxtfile::new; Not all lines contain the same number of records?'
            endif
        else if (self%access_type .eq. OPEN_TO_WRITE) then

            if(.not.fopen(self%funit, file=self%fname, iostat=ios, status='replace',iomsg=io_msg))&
                  call fileio_errmsg("nrtxt::new Error when opening file "&
                          //trim(self%fname)//' ; '//trim(io_msg),ios)
            if( ios .ne. 0 )then
                if(.not.fclose(self%funit,ios))then
                     call fileio_errmsg("nrtxt::new Error when closing file "&
                          //trim(self%fname)//' ; '//trim(io_msg),ios)
                 else
                     write(*,'(a)') 'simple_nrtxtfile::new; Error when opening file for writing: '&
                          //trim(self%fname)//' ; '//trim(io_msg)
                 end if
                stop
            endif
            if (present(wanted_recs_per_line)) then
                self%recs_per_line = wanted_recs_per_line
            else
                self%recs_per_line = 1
            endif
            self%ndatalines = 0
        endif
    end subroutine new

    subroutine readNextDataLine( self, read_data )
        class(nrtxtfile), intent(inout) :: self
        real,             intent(inout) :: read_data(:) !< output data
        character(len=LINE_MAX_LEN) :: buffer
        integer                     :: ios
        character(len=256)          :: io_message
        ! Check we are open to read
        if( self%access_type .ne. OPEN_TO_READ )then
            stop 'simple_nrtxtfile::readNextDataLine; File is not OPEN_TO_READ'
        endif
        ! Check the passed array is big enough
        if( size(read_data) .lt. self%recs_per_line )then
            stop 'simple_nrtxtfile::readNextDataLine; Supplied array is smaller than records per line'
        endif
        ! read the next line data line
        do
            read(self%funit, '(a)', iostat=ios, iomsg=io_message) buffer
            if( ios .ne. 0 )then
                call fileio_errmsg('simple_nrtxtfile::readNextDataLine; Encountered iostat error: '//trim(io_message),ios)
            endif
            if( .not. strIsComment(buffer) .and. .not. strIsBlank(buffer) )then
                buffer = trim(adjustl(buffer))
                read(buffer, *) read_data(1:self%recs_per_line)
                exit
            endif
        enddo
    end subroutine readNextDataLine

    subroutine writeDataLineReal( self, data_to_write )
        class(nrtxtfile), intent(inout) :: self
        real, intent(in)                :: data_to_write(:) !< input data
        integer                         :: record_counter
        if( self%access_type .ne. OPEN_TO_WRITE )then
            stop 'simple_nrtxtfile::writeDataLineReal; File is not OPEN_TO_WRITE'
        endif
        if( size(data_to_write) .lt. self%recs_per_line )then
             call simple_stop( 'simple_nrtxtfile::writeDataLineReal; Supplied array is smaller than records per line')
        endif
        do record_counter = 1, self%recs_per_line
            write(self%funit, '(g14.7,a)', advance='no') data_to_write(record_counter), ' '
        enddo
        ! finish the line
        write(self%funit,*)
        ! increase the number of data lines
        self%ndatalines = self%ndatalines+1
    end subroutine writeDataLineReal

    subroutine writeDataLineInt( self, data_to_write )
        class(nrtxtfile), intent(inout) :: self
        integer, intent(inout)          :: data_to_write(:) !< input data
        integer                         :: record_counter
        ! Check we are open to write
        if (self%access_type .ne. OPEN_TO_WRITE) then
            call simple_stop('simple_nrtxtfile::writeDataLineInt; File is not OPEN_TO_WRITE')
        endif
        ! Check size
        if( size(data_to_write) .lt. self%recs_per_line )then
            stop 'simple_nrtxtfile::writeDataLineInt; Supplied array is smaller than records per line'
        endif
        ! write out the line..
        do record_counter = 1,self%recs_per_line
            write(self%funit, '(g14.7,a)', advance='no') real(data_to_write(record_counter)), ' '
        enddo
        ! finish the line
        write(self%funit,*)
        ! increase the number of data lines
        self%ndatalines = self%ndatalines + 1
    end subroutine writeDataLineInt

    subroutine writeCommentLine( self, comment_to_write )
        class(nrtxtfile), intent(inout) :: self
        character(len=*), intent(in)    :: comment_to_write !< input data
        ! Check we are open to write
        if (self%access_type .ne. OPEN_TO_WRITE) then
            stop 'simple_nrtxtfile::writeCommentLine; File is not OPEN_TO_WRITE'
        endif
        ! write out the line..
        write(self%funit, '(2a)') '# ', trim(adjustl(comment_to_write))
    end subroutine writeCommentLine

    pure integer function get_nrecs_per_line( self )
        class(nrtxtfile), intent(in) :: self
        get_nrecs_per_line = self%recs_per_line
    end function get_nrecs_per_line

    pure integer function get_ndatalines( self )
        class(nrtxtfile), intent(in) :: self
        get_ndatalines = self%ndatalines
    end function get_ndatalines

    subroutine kill(self)
        class(nrtxtfile), intent(inout) :: self
        integer :: ios
        if( is_open(self%funit) )then
            if(.not.fclose(self%funit,iostat=ios))&
             call fileio_errmsg("nrtxtfile; kill; Failed to close ",ios)
        end if
        self%recs_per_line = 0
        self%ndatalines = 0
    end subroutine kill

end module simple_nrtxtfile
