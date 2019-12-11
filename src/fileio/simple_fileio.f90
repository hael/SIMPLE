! generic fileio module
module simple_fileio
use simple_defs
use simple_strings, only: upperCase,stringsAreEqual, strIsBlank, int2str,int2str_pad,cpStr
use simple_error,   only: allocchk, simple_exception, simple_error_check
use simple_syslib,  only: file_exists, is_open, is_file_open, is_io, simple_abspath,&
&exec_cmdline, del_file, syslib_copy_file
implicit none

public :: fileiochk, fopen, fclose, nlines,  filelength, funit_size, is_funit_open, get_open_funits
public :: add2fbody, swap_suffix, get_fbody, fname_new_ext, fname2ext, fname2iter, basename,  get_fpath
public :: make_dirnames, make_filenames, filepath, del_files, fname2format, read_filetable, write_filetable
public :: write_singlelineoftext, arr2file, file2rarr, simple_copy_file, make_relativepath
private
#include "simple_local_flags.inc"

interface fclose
    module procedure fclose_1
    module procedure fclose_2
end interface fclose

interface filepath
    module procedure filepath_1
    module procedure filepath_2
end interface filepath

integer, parameter :: MAX_UNIT_NUMBER = 1000

contains

    !> \brief  is for checking file IO status
    subroutine fileiochk( message, iostat , die)
        character(len=*),  intent(in)    :: message  !< error message
        integer, optional, intent(inout) :: iostat   !< error status
        logical, optional, intent(in)    :: die      !< do you want to terminate or not
        logical :: die_this
        integer :: iostat_this
        die_this=.true.
        iostat_this=2
        if( present(die) ) die_this=die
        if( present(iostat) ) iostat_this=iostat
        if( iostat_this /= 0 ) write(logfhandle,'(a)') message
        if (iostat_this == -1)then
            write(logfhandle,'(a)') "fileio: EOF reached "
        else if (iostat_this == -2) then
            write(logfhandle,'(a)') "fileio: End-of-record reached "
        else if( iostat_this /= 0 ) then
            if( die_this ) THROW_HARD('I/O')
        endif
    end subroutine fileiochk

    !> FOPEN enforce F2008 style open
    !!
    !! Usage: if(.not.fopen(fnr, fname, STATUS='REPLACE', action='WRITE', iostat=file_stat))&
    !!        call fileiochk('In: commander_rec :: volassemble', file_stat )
    !!
    subroutine fopen(funit, file, status, action, iostat, access, form, recl, async, pad,&
        &decimal, round, delim, blank, convert, iomsg, position, errmsg)
        integer,                    intent(inout) :: funit
        character(len=*),           intent(in)    :: file
        integer,          optional, intent(inout) :: iostat
        integer,          optional, intent(inout) :: recl
        character(len=*), optional, intent(in)    :: status, access, async, action, &
            &blank, pad, form, decimal, round, delim, convert, iomsg, position, errmsg
        integer               :: iostat_this,  recl_this
        character(len=STDLEN) :: filename,iomsg_this
        character(len=30)     :: async_this, access_this, action_this, status_this,&
            &blank_this, pad_this, decimal_this, delim_this, form_this, round_this,&
            &position_this, errmsg_this
        ! check to see if filename is empty
        write(filename,'(A)') trim(adjustl(file))
        if ( strIsBlank(filename) )then
            write(logfhandle,*) 'simple_system::fopen filename blank'
            if(present(iomsg))  write(logfhandle,*) trim(adjustl(iomsg))
            if(present(errmsg)) write(logfhandle,*) "Message: ", trim(adjustl(errmsg))
            return
        end if
        errmsg_this="In simple_fileio::fopen "
        if (present(errmsg)) errmsg_this=errmsg
        if (.not. (present(iostat) .or. present(form) .or. present(recl) .or.&
            & present(async) .or. present(pad) .or. present(action) .or. present(status)&
            & .or. present(position) .or. present(access) .or. present(decimal) .or. &
            present(round) .or. present(delim) .or. present(blank) ) )then
            open(NEWUNIT=funit, FILE=trim(adjustl(filename)),IOSTAT=iostat_this)
            call fileiochk(trim(adjustl(errmsg_this))//" fopen basic open "//trim(filename), iostat_this,.false.)
            if(is_io(funit)) THROW_HARD( "newunit returned "//int2str(funit))
            return ! if ok
        end if
        ! Optional args
        ! CONVERT ignored in file open argument list
        if(present(iostat))iostat_this=iostat
        !! Default
        write(action_this,'(A)') 'READWRITE'
        write(status_this,'(A)') 'UNKNOWN'
        write(position_this,'(A)') 'APPEND'
        if (present(status))then
            if (stringsAreEqual(status, 'OLD',.false.))     write(status_this,'(A)')  upperCase(status)
            if (stringsAreEqual(status, 'SCRATCH',.false.)) write(status_this,'(A)')  upperCase(status)
            if (stringsAreEqual(status, 'REPLACE',.false.)) write(status_this ,'(A)') upperCase(status)
            if (stringsAreEqual(status, 'NEW',.false.))     write( status_this,'(A)')  upperCase(status)
        end if
        ! ACTION: READ, WRITE, or READWRITE (default).
        if (present(action))then
            if (stringsAreEqual(action, 'WRITE',.false.)) write(action_this ,'(A)') upperCase(action)
            if (stringsAreEqual(action, 'READ',.false.))  write(action_this ,'(A)') upperCase(action)
        end if
        if ( (stringsAreEqual(status_this, 'NEW',.false.))  .and. &
            (stringsAreEqual(action_this, 'READ',.false.))  .and. &
            (.not. file_exists(filename) ) )then
            write(logfhandle,*) "::fopen incompatible status=NEW and action=READ ", trim(filename)," does not exist"
            return ! false
        end if
        if(present(position)) then
            if ( (stringsAreEqual(status_this, 'OLD',.false.))  .and. &
                (stringsAreEqual(position, 'APPEND',.false.))  .and. &
                (.not. file_exists(filename) ) )then
                write(logfhandle,*) "::fopen incompatible status=OLD and position=APPEND  when ",&
                    trim(filename)," does not exist"
                write( status_this,'(A)')  upperCase('NEW')
            end if
        end if
        ! access: DIRECT (random access) or SEQUENTIAL  or STREAM (F2003)
        write(access_this ,'(A)') 'SEQUENTIAL'
        if (present(access))then
            if (stringsAreEqual(access, 'DIRECT',.false.))  write(access_this ,'(A)') upperCase(access)
            if (stringsAreEqual(access, 'STREAM',.false.))then
                write(access_this ,'(A)') upperCase(access)
            endif
        end if
        recl_this=-1
        if(present(recl)) recl_this=recl
        !! Common file open
        if (.not.( present(form) .or. present(async) .or. present(pad) .or. &
            &present(decimal) .or. present(round) .or. present(delim) .or. present(blank) ) )then
            if (stringsAreEqual(access_this, 'DIRECT',.false.) .and. (recl_this > 0) ) then

                if (present(action))then
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,RECL=recl_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,RECL=recl_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this)
                    end if
                else ! no action
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,RECL=recl_this,&
                            &STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this,RECL=recl_this)
                    end if
                end if
            else
                if (present(action))then
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this)
                    end if
                else ! no action
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this)
                    end if
                end if
            end if
            ! call fileiochk(trim(adjustl(errmsg_this))//" fopen common open "//trim(filename), iostat_this,.false.)
            if(present(iostat))iostat=iostat_this
            if(funit/=0 .and. is_io(funit)) THROW_HARD("newunit returned "//int2str(funit))
            return
        end if
        recl_this=-1
        if(present(recl)) recl_this=recl
        write(pad_this,'(A)') 'YES'
        if (present(pad))then
            if (stringsAreEqual(pad, 'NO',.false.))  write(pad_this ,'(A)') upperCase(pad)
        end if
        write(async_this,'(A)')'NO'
        if (present(async))then
            if (stringsAreEqual(async, 'YES',.false.))  write(async_this ,'(A)')upperCase(async)
        end if
        write(blank_this,'(A)')'NULL'
        if (present(blank))then
            if (stringsAreEqual(blank, 'ZERO',.false.)) write( blank_this ,'(A)') upperCase(blank)
        end if
        write(decimal_this,'(A)')'POINT'
        if (present(decimal))then
            if (stringsAreEqual(decimal, 'COMMA',.false.)) write( decimal_this ,'(A)') upperCase(decimal)
        end if
        write(delim_this,'(A)')'NONE'
        if (present(delim))then
            if (stringsAreEqual(delim, 'APOSTROPHE',.false.)) write( delim_this ,'(A)') upperCase(delim)
            if (stringsAreEqual(delim, 'QUOTE',.false.))  write(delim_this ,'(A)') upperCase(delim)
        end if
        write(form_this,'(A)')'FORMATTED'
        if (present(form))then
            if (stringsAreEqual(form, 'UNFORMATTED',.false.))  write(form_this ,'(A)') upperCase(form)
            if (stringsAreEqual(form, 'BINARY',.false.))  write(form_this ,'(A)') upperCase(form)
        end if
        write(round_this,'(A)')'PROCESSOR_DEFINED'
        if (present(round))then
            if (stringsAreEqual(round, 'SUPPRESS',.false.)) write( round_this ,'(A)') upperCase(round)
            if (stringsAreEqual(round, 'PLUS',.false.))  write(round_this ,'(A)') upperCase(round)
            if (stringsAreEqual(round, 'UNDEFINED',.false.)) write( round_this ,'(A)') upperCase(round)
        end if
        if(present(iomsg)) iomsg_this=iomsg
        ! execute open under specific conditions
        if (stringsAreEqual(form_this, 'FORMATTED',.false.)) then
            open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                &ACTION=action_this,STATUS=status_this,ACCESS=access_this,&
                &BLANK=blank_this,FORM='FORMATTED', ROUND=round_this,&
                &IOMSG=iomsg_this)
        else
            if (stringsAreEqual(access_this, 'DIRECT',.false.))then
                open( NEWUNIT=funit, FILE=filename, IOSTAT=iostat_this, &
                    &ACTION=action_this, STATUS=status_this,&
                    &ACCESS=access_this, FORM=form_this, RECL=recl_this,&
                    &IOMSG=iomsg_this)
            else
                if (recl_this == -1)then
                    open( NEWUNIT=funit, FILE=filename, IOSTAT=iostat_this, &
                        &ACTION=action_this, STATUS=status_this,&
                        &ACCESS=access_this, FORM=form_this,&
                        &IOMSG=iomsg_this)
                else
                    open( NEWUNIT=funit, FILE=filename, IOSTAT=iostat_this, &
                        &ACTION=action_this, STATUS=status_this,&
                        &ACCESS=access_this, RECL=recl_this, FORM=form_this,&
                        &IOMSG=iomsg_this)
                end if
            end if
        end if
        call fileiochk(trim(adjustl(errmsg_this))//" fopen opening "//trim(filename), iostat_this, .false.)
        if(is_io(funit)) THROW_HARD("newunit returned "//int2str(funit))
        if(present(iostat))iostat=iostat_this
        if(present(recl))recl=recl_this
    end subroutine fopen

    !> FCLOSE replacement for intrinsic close
    !!
    !! Usage: call fclose( fnr,file_stat,errmsg='In: <src>::<func> failed ' )
    !!
    subroutine fclose_1 (funit,iostat,errmsg)
        integer,          intent(in)           :: funit
        integer,          intent(inout)        :: iostat
        character(len=*), intent(in), optional :: errmsg
        character(len=STDLEN) :: msg_this
        msg_this="SIMPLE_FILEIO::fclose_1 failed closing unit "//int2str(funit)
        if (present(errmsg)) write(msg_this,'(A)') trim(adjustl(errmsg))

        if (is_open(funit)) then
            CLOSE (funit,IOSTAT=iostat)
            call simple_error_check(iostat, trim(msg_this))
        end if
    end subroutine fclose_1

    subroutine fclose_2 (funit,errmsg,dispose,status)
        integer,          intent(in)           :: funit
        character(len=*), intent(in), optional :: status,dispose,errmsg
        character(len=30)     :: status_this
        character(len=STDLEN) :: msg_this
        integer               :: iostat
        if(is_io(funit)) then
            return !true
        end if
        !! status or dispose: 'KEEP' or 'DELETE', unless file was opened with status=SCRATCH
        write(status_this,'(A)')'KEEP'
        if (present(dispose))then
            if (stringsAreEqual(dispose, 'DELETE',.false.))  write(status_this ,'(A)') upperCase(dispose)
        end if
        if (present(status))then
            if (stringsAreEqual(status, 'DELETE',.false.))  write(status_this ,'(A)') upperCase(status)
        end if
        msg_this=" no message"
        if( present(errmsg) ) write(msg_this,'(A)') errmsg
        if( is_open(funit) ) close(funit,IOSTAT=iostat,STATUS=status_this)
    end subroutine fclose_2

    !> \brief return the number of lines in a textfile
    function nlines( fname ) result( n )
        use simple_strings
        character(len=*), intent(in)  :: fname !< input filename
        character(len=:), allocatable :: tfile
        integer          :: n, funit, ios,io_status
        character(len=1) :: junk
        if( file_exists(fname) )then
            tfile=fname
            call fopen(funit, tfile, status='unknown', action='read', iostat=io_status)
            call fileiochk(":nlines error opening file "//trim(tfile), io_status)
            n = 0
            do
                read(funit,*,IOSTAT=ios) junk
                if(ios /= 0)then
                    exit
                else
                    if(.not. strIsComment(junk))  n = n + 1
                endif
            end do
            call fclose_1( funit, io_status ,errmsg=" Error closing file in ::nlines "//trim(tfile))
        else
            n = 0
        endif
    end function nlines

    !> \brief  return the size of a binary file
    function filelength( fname ) result( filesz )
        character(len=*), intent(in) :: fname !< input filename
        integer                      :: filesz, funit, ios, cnt,recl
        character(len=1)             :: junk
        if(  file_exists(fname) )then
            recl=1
            call fopen(funit, fname, ACTION='read', IOSTAT=ios,&
                ACCESS='direct', form='unformatted', recl=recl)
            call fileiochk('simple_fileio :: nlines opening '//trim(fname), ios)
            cnt = 0
            filesz = 0
            do
                cnt = cnt+1
                read(funit,rec=cnt,IOSTAT=ios) junk
                if(ios /= 0)then
                    exit
                else
                    filesz = filesz+1
                endif
            end do
            call fclose_1( funit, ios ,errmsg=" Error closing file in ::filelength "//trim(fname))
        else
            filesz = 0
        endif
    end function filelength

    !> \brief  return file size in bytes
    function funit_size(unit) result(sz)
        integer, intent(in) :: unit !< input file unit
        integer(kind=8)     :: sz
        inquire(unit=unit,size=sz)
    end function funit_size

    !> \brief  return file opened state
    logical function is_funit_open(funit)
        integer, intent(in) :: funit !< input file unit
        inquire(unit=funit, OPENED=is_funit_open)
    end function is_funit_open

    ! computes an array of integers made of all currently opened units.
    ! Output : nbunits : number of opened units, units ( iunit ) : unit number for the opened unit
    ! #iunit with 1<= iunit <= nbunits
    subroutine get_open_funits( nbunits , units )
        integer, intent ( out ) :: nbunits
        integer, allocatable    :: units(:)
        integer :: iunit
        logical :: lopen
        integer :: step
        ! Loop over the steps.
        ! Step #1 : count the number of opened units
        ! Step #2 : store the number of opened units
        do step=1,2
            nbunits=0
            do iunit=1, MAX_UNIT_NUMBER
                if( iunit /= 5 .and. iunit /= 6 .and. iunit /= 9 ) then
                    inquire( UNIT = iunit, opened = lopen )
                    if( lopen )then
                        if( step == 1 )then ! count the number of opened units
                            nbunits = nbunits + 1
                        else                 ! store the number of opened units
                            nbunits = nbunits + 1
                            units ( nbunits ) = iunit
                        endif
                    end if
                end if
            end do
            ! At the end of step #1, allocate the array
            if( step == 1) allocate( units(nbunits) )
        enddo
    end subroutine get_open_funits

    !> \brief  is for adding to filebody
    function add2fbody( fname, suffix, str ) result( newname )
        character(len=*), intent(in)  :: fname, suffix, str
        character(len=:), allocatable :: newname
        integer :: pos
        pos = index(fname, suffix) ! position of suffix
        allocate(newname, source=fname(:pos-1)//trim(str)//trim(suffix))
    end function add2fbody

    !> \brief  is for substituting suffices
    function swap_suffix( fname, suffix, old_suffix ) result( newname )
        character(len=*), intent(in)  :: fname, suffix, old_suffix
        character(len=:), allocatable :: newname
        integer :: pos
        pos = index(fname, old_suffix) ! position of old_suffix
        allocate(newname, source=fname(:pos-1)//trim(suffix))
    end function swap_suffix

    !> \brief  is for extracting the body of a file
    function get_fbody( fname, suffix, separator ) result( fbody )
        character(len=*), intent(in) :: fname, suffix !< file extension
        character(len=STDLEN)        :: fbody
        logical,            optional :: separator
        integer :: pos
        logical :: l_separator
        l_separator = .true.
        if(present(separator))l_separator = separator
        if( l_separator )then
            pos = index(fname, '.'//suffix) ! position of suffix
        else
            pos = index(fname, suffix) ! position of suffix
        endif
        fbody = fname(:pos-1)
    end function get_fbody

    !> \brief  is for putting a new extension on filename
    !! \param fname Input filename
    function fname_new_ext( fname, suffix ) result( new_fname )
        character(len=*), intent(in)  :: fname, suffix !< filename and new file extension
        character(len=STDLEN)         :: fbody, new_fname
        character(len=:), allocatable :: ext
        ext   = fname2ext(trim(fname))
        fbody = get_fbody(trim(fname), ext)
        new_fname = trim(fbody)//'.'//trim(suffix)
    end function fname_new_ext

    !>  \brief  Return the 3-letter extension of a fname if present (without the period)
    pure function fname2ext( fname )
        character(len=*), intent(in)  :: fname    !< filename
        character(len=:), allocatable :: fname2ext
        integer :: length, pos
        length = len_trim(fname)
        pos = scan(fname(1:length),'.',back=.true.)
        if( pos == 0 )then
            allocate(fname2ext, source='   ')
        else
            allocate(fname2ext, source=trim(fname(pos+1:length)))
        endif
    end function fname2ext

    pure integer function fname2iter( fname )
        use simple_strings, only: map_str_nrs, str2int
        character(len=*), intent(in)  :: fname
        character(len=:), allocatable :: iter_num_ext, nrstr
        logical,          allocatable :: lnrs(:)
        integer :: ind, i, istat
        ind = index(fname, 'iter')
        if(ind == 0)then
            fname2iter = 0
            return
        endif
        iter_num_ext = fname(ind:)
        lnrs         = map_str_nrs(iter_num_ext)
        nrstr        = ''
        do i=1,len_trim(iter_num_ext)
            if( lnrs(i) ) nrstr = nrstr//iter_num_ext(i:i)
        end do
        call str2int(nrstr, istat, fname2iter)
    end function fname2iter

    !> strip directory from filenames (POSIX basename)
    pure function basename( fname ) result( new_fname)
        character(len=*), intent(in)  :: fname     !< abs filename
        character(len=:), allocatable :: new_fname
        integer :: length, pos
        length = len_trim(fname)
        pos = scan(fname(1:length),PATH_SEPARATOR,back=.true.)
        if( pos == 0 )then
            allocate(new_fname, source=trim(fname))
        else
            allocate(new_fname, source=trim(fname(pos+1:length)))
        endif
    end function basename

    !> get path from file name  (POSIX dirname)
    pure function get_fpath( fname ) result( path )
        character(len=*), intent(in)  :: fname !< abs filename
        character(len=:), allocatable :: path
        integer :: pos
        pos = scan(fname,PATH_SEPARATOR,back=.true.)
        if( pos == 0 )then
            allocate(path, source=PATH_HERE)
        else
            allocate(path, source=trim(fname(:pos)))
        endif
    end function get_fpath

    !>  \brief  returns numbered names (body) with 0-padded integer strings
    function make_dirnames( body, n, numlen ) result( names )
        use simple_strings, only: int2str, int2str_pad
        character(len=*),  intent(in) :: body
        integer,           intent(in) :: n
        integer, optional, intent(in) :: numlen
        character(len=STDLEN), allocatable :: names(:)
        integer :: nnumlen, i
        nnumlen = len(int2str(n))
        if( present(numlen) ) nnumlen = numlen
        allocate(names(n))
        do i=1,n
            names(i) = trim(body)//int2str_pad(i, nnumlen)
        end do
    end function make_dirnames

    !>  \brief  returns numbered file-names with 0-padded integer strings
    function make_filenames( body, n, ext, numlen, suffix ) result( names )
        use simple_strings, only: int2str, int2str_pad
        character(len=*),           intent(in) :: body, ext
        integer,                    intent(in) :: n
        integer,          optional, intent(in) :: numlen
        character(len=*), optional, intent(in) :: suffix
        character(len=STDLEN), allocatable     :: names(:)
        integer :: nnumlen, i
        logical :: suffix_present
        nnumlen = len(int2str(n))
        if( present(numlen) ) nnumlen = numlen
        suffix_present = present(suffix)
        allocate(names(n))
        do i=1,n
            if( suffix_present )then
                names(i) = trim(body)//int2str_pad(i,nnumlen)//trim(suffix)//ext
            else
                names(i) = trim(body)//int2str_pad(i,nnumlen)//ext
            endif
        end do
    end function make_filenames

    !> concatenate strings together to with '/' to create a filename
    !! input args are restricted to STDLEN after trimming
    !! for allocatable character results
    !! e.g. fname = filepath('dir1',trim(dirfixed), diralloc, file//ext )
    function filepath_1(p1, p2, p3, p4) result( fname )
        character(len=*),  intent(in)           :: p1
        character(len=*),  intent(in)           :: p2
        character(len=*),  intent(in), optional :: p3
        character(len=*),  intent(in), optional :: p4
        character(len=:), allocatable     :: fname
        character(len=STDLEN)      :: s1,s2,s3,s4
        integer :: endpos1,endpos2,endpos3,endpos4,startpos2,startpos3,startpos4
        startpos2=1;startpos3=1;startpos4=1
        s1 = trim(adjustl(p1))
        endpos1 = len_trim(s1)
        ! Exceptions
        if(endpos1<1) THROW_HARD("first arg too small")
        if(s1(endpos1:endpos1)==PATH_SEPARATOR) endpos1 = endpos1-1
        if(endpos1<1) THROW_HARD("first arg cannot be /")
        s2 = trim(adjustl(p2))
        endpos2 = len_trim(s2)
        if(endpos2<1) THROW_HARD("second arg too small")
        if(s2(endpos2:endpos2)==PATH_SEPARATOR) endpos2 = endpos2-1
        if(endpos2==0) THROW_HARD("second arg cannot be / ")
        if(s2(1:1)==PATH_SEPARATOR) startpos2 = 2
        if(endpos2-startpos2==0) THROW_HARD("second arg cannot be "//trim(s2))
        if(present(p3))then
            s3 = trim(adjustl(p3))
            endpos3 = len_trim(s3)
            if(endpos3<1) THROW_HARD("third arg too small")
            if(s3(endpos3:endpos3)==PATH_SEPARATOR) endpos3 = endpos3-1
            if(endpos3<1) THROW_HARD("third arg cannot be /")
            if(s3(1:1)==PATH_SEPARATOR) startpos3 = 2
            if(endpos3-startpos3==0) THROW_HARD("third arg cannot be "//trim(s3))

            if(present(p4))then
                s4 = trim(adjustl(p4))
                endpos4 = len_trim(s4)
                if(endpos4==0) THROW_HARD("fourth arg too small")
                if(s4(endpos4:endpos4)==PATH_SEPARATOR) endpos4 = endpos4-1
                if(endpos4==0) THROW_HARD("fourth arg  cannot be /")
                if(s4(1:1)==PATH_SEPARATOR) startpos4 = 2
                if(endpos4-startpos4==0) THROW_HARD("fourth arg cannot be "//trim(s4))
                allocate(fname,source=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos2:endpos2)//&
                    &PATH_SEPARATOR//s3(startpos3:endpos3)//PATH_SEPARATOR//s4(startpos4:endpos4))
            else
                !! Concat arg1 arg2 and arg3
                allocate(fname,source=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos2:endpos2)//&
                    &PATH_SEPARATOR//s3(startpos3:endpos3))
            endif
        else
            !! Concat arg1 and arg2
            allocate(fname,source=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos2:endpos2))
        endif
    end function filepath_1

    !> concatenate strings together to create a filename similar to filepath_1 above
    !! for non-allocatable character results
    !! Third and Fourth args are optional
    function filepath_2(p1, p2, p3, p4, nonalloc ) result(fname)
        character(len=*),  intent(in)           :: p1
        character(len=*),  intent(in)           :: p2
        character(len=*),  intent(in), optional :: p3
        character(len=*),  intent(in), optional :: p4
        logical :: nonalloc
        character(len=STDLEN)      :: s1,s2,s3,s4
        character(len=LONGSTRLEN)  ::fname
        integer :: endpos1,endpos2,endpos3,endpos4,startpos2,startpos3,startpos4
        startpos2=1;startpos3=1;startpos4=1
        s1 = trim(adjustl(p1))
        endpos1 = len_trim(s1)
        if(endpos1<1) THROW_HARD("first arg too small")
        if(s1(endpos1:endpos1)==PATH_SEPARATOR) endpos1 = endpos1-1
        if(endpos1<1) THROW_HARD("first arg cannot be / ")
        s2 = trim(adjustl(p2))
        endpos2 = len_trim(s2)
        if(endpos2<1) THROW_HARD("second arg too small")
        if(s2(endpos2:endpos2)==PATH_SEPARATOR) endpos2 = endpos2-1
        if(endpos2==0) THROW_HARD("second arg cannot be / ")
        if(s2(1:1)==PATH_SEPARATOR) startpos2 = 2
        if(endpos2-startpos2==0) THROW_HARD("second arg cannot be "//trim(s2))
        if(present(p3))then
            s3 = trim(adjustl(p3))
            endpos3 = len_trim(s3)
            if(endpos3<1) THROW_HARD("third arg too small")
            if(s3(endpos3:endpos3)==PATH_SEPARATOR) endpos3 = endpos3-1
            if(endpos3<1) THROW_HARD("third arg cannot be /")
            if(s3(1:1)==PATH_SEPARATOR) startpos3 = 2
            if(endpos3-startpos3==0) THROW_HARD("third arg cannot be "//trim(s3))

            if(present(p4))then
                s4 = trim(adjustl(p4))
                endpos4 = len_trim(s4)
                if(endpos4==0) THROW_HARD("fourth arg too small")
                if(s4(endpos4:endpos4)==PATH_SEPARATOR) endpos4 = endpos4-1
                if(endpos4==0) THROW_HARD("fourth arg  cannot be /")
                if(s4(1:1)==PATH_SEPARATOR) startpos4 = 2
                if(endpos4-startpos4==0) THROW_HARD("fourth arg cannot be "//trim(s4))
                !! concatenate four pathnames
                fname=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos2:endpos2)//&
                    &PATH_SEPARATOR//s3(startpos3:endpos3)//PATH_SEPARATOR//s4(startpos4:endpos4)
            else
                !! concatenate three pathnames
                fname=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos2:endpos2)//&
                    &PATH_SEPARATOR//s3(startpos3:endpos3)
            endif
        else
            !! concatenate two pathnames
            fname=s1(1:endpos1)//PATH_SEPARATOR//s2(startpos4:endpos2)
        endif
    end function filepath_2

    !> \brief  is for deleting consecutively numbered files with padded number strings
    subroutine del_files( body, n, ext, numlen, suffix )
        character(len=*),           intent(in) :: body !< input filename body
        integer,                    intent(in) :: n    !< total num for del, formatted as body[n].ext
        character(len=*), optional, intent(in) :: ext  !< input filename extension
        integer,          optional, intent(in) :: numlen !< number length
        character(len=*), optional, intent(in) :: suffix !< file suffix
        character(len=STDLEN), allocatable :: names(:)
        integer :: ifile
        if( present(ext) )then
            names = make_filenames( body, n, ext, numlen=numlen, suffix=suffix )
        else
            names = make_dirnames( body, n, numlen=numlen)
        endif
        do ifile=1,n
            if( file_exists(names(ifile)) ) call del_file(names(ifile))
        end do
    end subroutine del_files

    !>  \brief Return a one letter code for the file format designated by the extension in the fname
    !!         if .mrc: M
    !!         if .spi: S
    !!         if .img: I
    !!         if .hed: I
    !!         else: N
    pure function fname2format( fname )
        character(len=*), intent(in)  :: fname        !< input filename
        character(len=1)              :: fname2format
        character(len=:), allocatable :: extension
        extension = fname2ext(fname)
        select case(extension)
        case ('img','hed')
            fname2format = 'I'
        case ('mrc','map','ctf','mrcs')
            fname2format = 'M'
        case ('spi')
            fname2format = 'S'
        case('bin','raw','sbin')
            fname2format = 'B'
        case('dbin')
            fname2format = 'D'
        case('txt', 'asc', 'box','dat')
            fname2format = 'T'
        case('pdb')
            fname2format = 'P'
        case('simple')
            fname2format = 'O'
        case('star')
            fname2format = 'R'
#ifdef USING_TIFF
        case('tif','tiff')
            fname2format = 'J'
#endif
        case DEFAULT
            fname2format = 'N'
        end select
    end function fname2format

    !>  \brief  reads a filetable into an array
    subroutine read_filetable( filetable, filenames )
        character(len=*),                       intent(in)  :: filetable    !< input table filename
        character(len=LONGSTRLEN), allocatable, intent(out) :: filenames(:) !< array of filenames
        integer :: nl, funit, iline,io_stat
        nl = nlines(trim(filetable))
        call fopen(funit,filetable,'old','unknown',io_stat)
        call fileiochk("read_filetable failed to open file "//trim(filetable),io_stat )
        allocate( filenames(nl), stat=alloc_stat )
        if(alloc_stat /= 0) call allocchk ('In: read_filetable; simple_fileio  ', alloc_stat)
        do iline=1,nl
            read(funit,'(a1024)') filenames(iline)
        end do
        call fclose_1(funit,io_stat)
        call fileiochk("read_filetable failed to close",io_stat)
    end subroutine read_filetable

    !>  \brief  writes a filetable array to a text file
    subroutine write_filetable( filetable, filenames )
        character(len=*),          intent(in)  :: filetable    !< output table filename
        character(len=LONGSTRLEN), intent(in)  :: filenames(:) !< array of filenames
        integer :: nl, funit, iline, io_stat
        nl = size(filenames)
        call fopen(funit,filetable, 'replace', 'unknown', io_stat)
        call fileiochk("write_filetable failed to open file "//filetable,io_stat )
        do iline=1,nl
            write(funit,'(a)') trim(filenames(iline))
        end do
        call fclose_1(funit,io_stat)
        call fileiochk("write_filetable failed to close",io_stat)
    end subroutine write_filetable

    !>  \brief  (over)writes a file with a single file of text
    subroutine write_singlelineoftext( filename, text )
        character(len=*), intent(in)  :: filename, text
        integer :: funit, io_stat
        call fopen(funit,filename, 'replace', 'unknown', io_stat)
        call fileiochk("write_singlelineoftext failed to open file: "//filename,io_stat )
        write(funit,'(A)')trim(text)
        call fclose_1(funit,io_stat)
        call fileiochk("write_singlelineoftext failed to close: "//filename,io_stat)
    end subroutine write_singlelineoftext

    !> \brief  for converting a real array 2 file
    subroutine arr2file( arr, fnam )
        real,             intent(in) :: arr(:)    !< array of filenames
        character(len=*), intent(in) :: fnam      !< input table filename
        real    :: rval
        integer :: recsz, i, funit,io_stat
        inquire(iolength=recsz)rval
        rval = size(arr)
        funit=-1
        call fopen(funit,fnam,'replace','unknown', iostat=io_stat,access='direct',form='unformatted',recl=recsz)
        call fileiochk("arr2file_1 fopen failed "//trim(fnam),io_stat)
        write(funit, rec=1, iostat=io_stat) rval
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        call fclose_1(funit,io_stat, errmsg="arr2file_1 fclose_1 failed "//trim(fnam))
    end subroutine arr2file

    !> \brief  for converting a file generated by arr2file back to an array
    function file2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam  !< input table filename
        real, allocatable            :: arr(:) !< array of filenames
        real    :: rval
        integer :: recsz, i, n, funit,io_stat
        if( file_exists(trim(fnam)) )then
            inquire(iolength=recsz) rval
            call fopen(funit,fnam,'old','unknown', io_stat,'direct','unformatted',recl=recsz)
            call fileiochk("file2rarr fopen failed "//trim(fnam),io_stat)
            read(funit, rec=1,iostat=io_stat) rval
            n = nint(rval)
            allocate( arr(n), stat=alloc_stat )
            if(alloc_stat /= 0) call allocchk('In: file2arr; simple_fileio ', alloc_stat)
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            call fclose_1(funit,io_stat,errmsg="file2rarr fclose_1 failed "//trim(fnam))
        else
            THROW_HARD(trim(fnam)//' does not exist')
        endif
    end function file2rarr

    subroutine simple_copy_file(fname1, fname2, status)
        character(len=*),  intent(in)  :: fname1, fname2 !< input filenames
        integer, optional, intent(out) :: status
        call syslib_copy_file(fname1, fname2, status)
    end subroutine simple_copy_file

    ! builds a relative path with respect to working directory
    ! given absolute paths of directory and filename returns file relative path
    ! with respect to working directory (minus execution directory, eg 1_xxx)
    ! if the file lies outside the working directory then the absolute file path is kept
    subroutine make_relativepath(cwd, fname, newfname)
        character(len=*),          intent(in)  :: cwd, fname
        character(len=LONGSTRLEN), intent(out) :: newfname
        character(len=:), allocatable :: fname_here
        character(LONGSTRLEN)         :: cwd_here, projdir
        integer                       :: l_cwd, l_fname, l, slashpos_left, slashpos_right, ipos
        ! get absolute filename if necessary
        if( fname(1:1).eq.'/' )then
            ! was already absolute
            fname_here = trim(fname)
            if( .not.file_exists(fname_here) )THROW_HARD('File does not exist: '//trim(fname_here))
        else
            fname_here = simple_abspath(fname, errmsg='simple_fileio::make_relativepath: '//trim(fname))
        endif
        ! remove final '/' from cwd for safety
        l_cwd = len_trim(cwd)
        if(cwd(l_cwd:l_cwd)=='/')then
            cwd_here = cwd(1:l_cwd-1)
            l_cwd = l_cwd-1
        else
            cwd_here = cwd(1:l_cwd)
        endif
        ! removes execution directory from cwd
        slashpos_left  = scan(cwd_here,'/',back=.false.)
        slashpos_right = scan(cwd_here,'/',back=.true.)
        if(slashpos_right > slashpos_left)then
            projdir = cwd_here(1:slashpos_right-1)
        else
            THROW_HARD('Incorrect directory structure; simple_fileio :: make_relativepath')
        endif
        ! builds path
        ipos = index(trim(fname_here),trim(projdir))
        if( ipos == 1 )then
            ! file inside project directory, use relative path
            l = len_trim(projdir)
            l_fname = len_trim(fname_here)
            if(fname_here(l+1:l+1) == '/')then
                newfname = '..'//fname_here(l+1:l_fname)
            else
                newfname = '../'//fname_here(l+1:l_fname)
            endif
        else
            ! file outside of project directory or weird structure, use absolute path
            newfname = trim(fname_here)
        endif
    end subroutine make_relativepath

end module simple_fileio
