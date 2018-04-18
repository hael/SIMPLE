! generic fileio module
module simple_fileio
use simple_defs
use, intrinsic :: iso_fortran_env, only: stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT, stdin=>INPUT_UNIT
use simple_strings, only: upperCase,stringsAreEqual, strIsBlank, int2str,int2str_pad,cpStr
use simple_error,   only: allocchk, simple_stop, simple_error_check
use simple_syslib, only: file_exists, is_open, is_file_open, is_io,  &
    &exec_cmdline, del_file, simple_list_files, simple_glob_list_tofile
implicit none

interface arr2file
    module procedure arr2file_1
!    module procedure arr2file_2
end interface arr2file

interface fclose
    module procedure fclose_1
    module procedure fclose_2
end interface fclose

integer, parameter :: MAX_UNIT_NUMBER = 1000

contains

    !> \brief  is for checking file IO status
    subroutine fileiochk( message, iostat , die)
        character(len=*), intent(in)              :: message  !< error message
        integer,          intent(inout), optional :: iostat   !< error status
        logical,          intent(in),    optional :: die      !< do you want to terminate or not
        logical :: die_this
        integer :: iostat_this
        die_this=.true.
        iostat_this=2
        if(present(die)) die_this=die
        if (present(iostat)) iostat_this=iostat
        if( iostat_this /= 0 ) write(stderr,'(a)') message
        if (iostat_this == -1)then
            write(stderr,'(a)') "fileio: EOF reached (PGI version)"
        else if (iostat_this == -2) then
            write(stderr,'(a)') "fileio: End-of-record reached (PGI version)"
        else if( iostat_this /= 0 ) then
            call simple_error_check(iostat_this,'File I/O Error#'//int2str(iostat_this))
            if(die_this)call simple_stop('File I/O stop ',__FILENAME__,__LINE__)
        endif
    end subroutine fileiochk

    !> FOPEN enforce F2008 style open so that PGI/Intel behave correctly
    !!
    !! Usage: if(.not.fopen(fnr, fname, STATUS='REPLACE', action='WRITE', iostat=file_stat))&
    !!        call fileiochk('In: commander_rec :: eo_volassemble', file_stat )
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
            print *, 'simple_system::fopen filename blank'
            if(present(iomsg))  print *, trim(adjustl(iomsg))
            if(present(errmsg))  print *, "Message: ", trim(adjustl(errmsg))
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
            if(is_io(funit)) call simple_stop( "simple_fileio::fopen newunit returned "//int2str(funit) ,&
                __FILENAME__,__LINE__)
            return ! if ok
        end if
        ! Optional args
        if(present(convert)) then
            ! print *, 'CONVERT ignored in file open argument list'
        end if
        if(present(iostat))iostat_this=iostat
        !! Default
        write(action_this,'(A)') 'READWRITE'
        write(status_this,'(A)') 'UNKNOWN'
        write(position_this,'(A)') 'APPEND'
        if (present(status))then
            if (stringsAreEqual(status, 'OLD',.false.))  write(status_this,'(A)')  upperCase(status)
            if (stringsAreEqual(status, 'SCRATCH',.false.)) write(status_this,'(A)')  upperCase(status)
            if (stringsAreEqual(status, 'REPLACE',.false.)) write(status_this ,'(A)') upperCase(status)
            if (stringsAreEqual(status, 'NEW',.false.)) write( status_this,'(A)')  upperCase(status)
        end if
        ! ACTION: READ, WRITE, or READWRITE (default).
        if (present(action))then
            if (stringsAreEqual(action, 'WRITE',.false.))  write(action_this ,'(A)') upperCase(action)
            if (stringsAreEqual(action, 'READ',.false.))  write(action_this ,'(A)') upperCase(action)
        end if
        if ( (stringsAreEqual(status_this, 'NEW',.false.))  .and. &
            (stringsAreEqual(action_this, 'READ',.false.))  .and. &
            (.not. file_exists(filename) ) )then
            print *, "::fopen incompatible status=NEW and action=READ ", trim(filename)," does not exist"
            return ! false
        end if
        if(present(position)) then
            if ( (stringsAreEqual(status_this, 'OLD',.false.))  .and. &
                (stringsAreEqual(position, 'APPEND',.false.))  .and. &
                (.not. file_exists(filename) ) )then
                print *, "::fopen incompatible status=OLD and position=APPEND  when ",&
                    trim(filename)," does not exist"
                write( status_this,'(A)')  upperCase('NEW')
            end if
        end if
        ! access: DIRECT (random access) or SEQUENTIAL  or STREAM (F2003)
        write(access_this ,'(A)') 'SEQUENTIAL'
        if (present(access))then
            if (stringsAreEqual(access, 'DIRECT',.false.))  write(access_this ,'(A)') upperCase(access)
            if (stringsAreEqual(access, 'STREAM',.false.))then
#ifdef PGI
               ! print *,"** Cannot 'stream' in current PGI version, using DIRECT"
                write(access_this,'(A)') upperCase(access)
                if (stringsAreEqual(status_this,'NEW',.false.)) write( status_this,'(A)')  'UNKNOWN'
#else
                write(access_this ,'(A)') upperCase(access)
#endif
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
            call fileiochk(trim(adjustl(errmsg_this))//" fopen common open "//trim(filename), iostat_this,.false.)
            if(present(iostat))iostat=iostat_this
            if(funit/=0 .and. is_io(funit)) call simple_stop( "::fopen newunit returned "//int2str(funit),&
                __FILENAME__,__LINE__)
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
        if(is_io(funit)) call simple_stop( "::fopen newunit returned "//int2str(funit) ,__FILENAME__,__LINE__)
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
        if(present(errmsg)) write(msg_this,'(A)') errmsg
        if (is_open(funit)) then
            CLOSE (funit,IOSTAT=iostat,STATUS=status_this)
            call simple_error_check(iostat, "SIMPLE_FILEIO::fclose_1 failed closing unit "//int2str(funit)//&
                &" ; "//trim(msg_this))
        end if
    end subroutine fclose_2


    !> \brief return the number of lines in a textfile
    function nlines( fname ) result( n )
        character(len=*), intent(in)    :: fname !< input filename
        character(len=:), allocatable   :: tfile
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
                    n = n + 1
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
        integer, intent(in)          :: unit !< input file unit
        integer(kind=8)              :: sz
        inquire(unit=unit,size=sz)
    end function funit_size

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


    !> \brief  is for making a file-table (to be able to commander execute programs that depend on them)
    ! subroutine make_filetable( tabname, n, body, ext, numlen, suffix )
    !     character(len=*),           intent(in) :: tabname !< file-table (string)
    !     integer,                    intent(in) :: n       !< total num of files
    !     character(len=*),           intent(in) :: body    !< filename body
    !     character(len=*),           intent(in) :: ext     !< filename extension
    !     integer,          optional, intent(in) :: numlen  !< number length
    !     character(len=*), optional, intent(in) :: suffix  !< file suffix
    !     character(len=STDLEN), allocatable :: names(:)
    !     integer :: ifile, fnr, ios
    !     names = make_filenames( body, n, ext, numlen=numlen, suffix=suffix )
    !     call fopen(fnr, file=tabname, status='replace', action='write', iostat=ios)
    !     call fileiochk('simple_fileio :: make_filetable '//trim(tabname), ios)
    !     do ifile=1,n
    !         write(fnr,'(a)') trim(names(ifile))
    !     end do
    !     call fclose_1( fnr, ios )
    !     call fileiochk(" Error closing file in ::make_filetable ",ios)
    ! end subroutine make_filetable

    !> \brief  return file size in bytes
    ! function file_size(fname) result(sz)
    !     character(len=*), intent(in) :: fname !< input filename
    !     integer(kind=8)              :: sz
    !     inquire(file=trim(adjustl(fname)),size=sz)
    ! end function file_size

    !> \brief  is for checking file kind
    !> \param fname,suffix string args to check suffix
    ! function file_kind( fname, suffix ) result( yep )
    !     character(len=*), intent(in) :: fname, suffix
    !     integer :: pos
    !     logical :: yep
    !     pos = index(fname, suffix) ! position of suffix
    !     if( pos == 0 )then
    !         yep = .false.
    !     else
    !         yep = .true.
    !     endif
    ! end function file_kind

    !> \brief  is for adding to filebody
    function add2fbody( fname, suffix, str ) result( newname )
        character(len=*), intent(in)  :: fname, suffix, str
        character(len=:), allocatable :: newname
        integer :: pos
        pos = index(fname, suffix) ! position of suffix
        allocate(newname, source=fname(:pos-1)//trim(str)//trim(suffix))
    end function add2fbody

    ! function add2fbody_and_new_ext( fname, suffix, str, new_ext ) result( newname )
    !     character(len=*), intent(in)  :: fname, suffix, str, new_ext
    !     character(len=:), allocatable :: newname
    !     integer :: pos
    !     pos = index(fname, suffix) ! position of suffix
    !     allocate(newname, source=fname(:pos-1)//trim(str)//trim(new_ext))
    ! end function add2fbody_and_new_ext

    !> \brief  is for deleting from fbody
    ! function del_from_fbody( fname, suffix, str ) result( newname )
    !     character(len=*), intent(in)  :: fname, suffix, str
    !     character(len=:), allocatable :: newname
    !     integer :: pos
    !     pos = index(fname, str) ! position of str
    !     allocate(newname, source=fname(:pos-1)//trim(suffix))
    ! end function del_from_fbody

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

    !! Used once in simple_picker
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

    !> strip directory and suffix from filenames
    pure function basename( fname ) result( new_fname)
        character(len=*), intent(in)  :: fname     !< abs filename
        character(len=:), allocatable :: new_fname
        integer :: length, pos
        length = len_trim(fname)
        pos = scan(fname(1:length),'/',back=.true.)
        if( pos == 0 )then
            allocate(new_fname, source=trim(fname))
        else
            allocate(new_fname, source=trim(fname(pos+1:length)))
        endif
    end function

    !! Unused
    !> strip last component from file name
    pure function dirname( fname ) result( abspath )
        character(len=*), intent(in)  :: fname !< abs filename
        character(len=:), allocatable :: abspath !< abs file path
        integer :: length, pos
        length = len_trim(fname)
        pos = scan(fname(1:length),'/',back=.true.)
        allocate(abspath, source=trim(fname(1:pos)))
    end function

    !>  \brief  returns the integer number identifier of a filename
    ! subroutine fname2ind( str, ivar )
    !     use simple_strings, only: map_str_nrs, str2int
    !     character(len=*), intent(in)  :: str    !< abs filename
    !     integer,          intent(out) :: ivar   !< file index number
    !     logical, allocatable          :: pos(:)
    !     character(len=:), allocatable :: str_copy
    !     integer :: j, lstr, io_stat, nrrange(2)
    !     lstr = len(str);  nrrange = 0
    !     pos = map_str_nrs(str)
    !     if( any(pos) )then
    !         do j=lstr,1,-1
    !             if( pos(j) )then
    !                 nrrange(1) = j
    !                 nrrange(2) = j
    !                 do while( pos(nrrange(1)) )
    !                     nrrange(1) = nrrange(1)-1
    !                 end do
    !                 nrrange(1) = nrrange(1)+1
    !                 exit
    !             endif
    !         end do
    !         allocate(str_copy, source=str(nrrange(1):nrrange(2)))
    !         call str2int(str_copy, io_stat, ivar)
    !     else
    !         allocate(str_copy, source='1')
    !         call str2int(str_copy, io_stat, ivar)
    !     endif
    !     if( allocated(pos)      ) deallocate(pos)
    !     if( allocated(str_copy) ) deallocate(str_copy)
    ! end subroutine fname2ind

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
            case ('mrc','map','st','ctf','mrcs')
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
            case DEFAULT
                fname2format = 'N'
        end select
    end function fname2format

    !! Unused
    !>  \brief  to check if same file format
    !! \param fname1,fname2 input filenames
    ! pure logical function same_format( fname1, fname2 )
    !     character(len=*), intent(in) :: fname1, fname2
    !     character(len=1) :: form1, form2
    !     form1 = fname2format(fname1)
    !     form2 = fname2format(fname2)
    !     same_format = form1 == form2
    ! end function same_format

    !>  \brief  reads a filetable into an array
    subroutine read_filetable( filetable, filenames )
        character(len=*),                   intent(in)  :: filetable    !< input table filename
        character(len=STDLEN), allocatable, intent(out) :: filenames(:) !< array of filenames
        integer :: nl, funit, iline,io_stat
        nl = nlines(filetable)
        call fopen(funit,filetable,'old','unknown',io_stat)
        call fileiochk("read_filetable failed to open file "//trim(filetable),io_stat )
        allocate( filenames(nl), stat=alloc_stat )
        if(alloc_stat /= 0) call allocchk ('In: read_filetable; simple_fileio  ', alloc_stat)
        do iline=1,nl
            read(funit,'(a256)') filenames(iline)
        end do
        call fclose_1(funit,io_stat)
        call fileiochk("read_filetable failed to close",io_stat)
    end subroutine read_filetable

    !>  \brief  writes a filetable array to a text file
    subroutine write_filetable( filetable, filenames )
        character(len=*),      intent(in)  :: filetable  !< output table filename
        character(len=STDLEN), intent(in)  :: filenames(:)!< array of filenames
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

    !!Unused
    !> \brief  for converting a file generated by txtfile2arr back to an array
    ! function txtfile2rarr( fnam ) result( arr )
    !     character(len=*), intent(in) :: fnam    !< input table filename
    !     real, allocatable :: arr(:)             !< array of filenames
    !     integer :: i, n, funit, io_stat
    !     if( file_exists(trim(fnam)) )then
    !         n = nlines(fnam)
    !         allocate( arr(n), stat=alloc_stat )
    !         if(alloc_stat /= 0) call allocchk('In: txtfile2rarr; simple_fileio  ', alloc_stat)
    !         call fopen(funit,fnam,'old','unknown',io_stat)
    !         call fileiochk("txtfile2rarr failed to open  "//trim(fnam), io_stat)
    !         do i=1,n
    !             read(funit,*) arr(i)
    !         end do
    !         call fclose_1(funit,io_stat,errmsg="txtfile2rarr failed to close "//trim(fnam))
    !     else
    !         write(*,*) fnam
    !         stop 'file does not exist; txtfile2rarr; simple_fileio      '
    !     endif
    ! end function txtfile2rarr

    ! !! Unused
    ! !> \brief  merging two text files into a single array
    ! !! \param file1,file2 input filenames for merging
    ! function merge_txtfiles( file1, file2 )  result( arr )
    !     character(len=*),      intent(in)  :: file1, file2
    !     character(len=STDLEN), allocatable :: arr(:)
    !     integer :: n1, n2, cnt, funit, i, io_stat
    !     logical :: here(2)
    !     here(1) = file_exists(trim(file1))
    !     here(2) = file_exists(trim(file2))
    !     n1 = 0
    !     n2 = 0
    !     if( file_exists(trim(file1)) ) n1 = nlines(file1)
    !     if( file_exists(trim(file2)) ) n2 = nlines(file2)
    !     allocate( arr(n1+n2), stat=alloc_stat )
    !     if(alloc_stat /= 0) call allocchk('In: merge_txtfiles; simple_fileio  ', alloc_stat)
    !     if( here(1) )then
    !         call fopen(funit,file1,'old','unknown',io_stat)
    !         call fileiochk("merge_txtfiles failed "//trim(file1), io_stat)
    !         cnt = 0
    !         do i=1,n1
    !             cnt = cnt+1
    !             read(funit,*) arr(cnt)
    !         end do
    !         call fclose_1(funit,io_stat)
    !         call fileiochk("merge_txtfiles failed to close "//trim(file1), io_stat)
    !         if( .not. here(2) ) return
    !     else
    !         call fopen(funit,file2,'old','unknown',io_stat)
    !         call fileiochk("merge_txtfiles failed to open "//trim(file2), io_stat)
    !         cnt = 0
    !         do i=1,n2
    !             cnt = cnt+1
    !             read(funit,*) arr(cnt)
    !         end do
    !         call fclose_1(funit,io_stat)
    !         call fileiochk("merge_txtfiles failed to close "//trim(file2), io_stat)
    !         return
    !     endif
    !     call fopen(funit,file2,'old','unknown',io_stat)
    !     call fileiochk("merge_txtfiles failed to open "//trim(file2), io_stat)
    !     do i=1,n2
    !         cnt = cnt+1
    !         read(funit,*, iostat=io_stat) arr(cnt)
    !     end do
    !     call fclose_1(funit,io_stat,errmsg="merge_txtfiles failed to close "//trim(file2))
    ! end function merge_txtfiles

    ! !! Unused
    ! !> \brief  for converting a file generated by file2arr back to an array
    ! function file2iarr( fnam ) result( arr )
    !     character(len=*), intent(in)  :: fnam             !< input table filename
    !     integer,          allocatable :: arr(:)           !< array of filenames
    !     integer :: recsz, i, n, funit, ival, io_stat
    !     if( file_exists(trim(fnam)) )then
    !         inquire(iolength=recsz) ival
    !         call fopen(funit,fnam,'OLD','unknown', io_stat,'direct','unformatted',recsz)
    !         call fileiochk("file2iarr fopen failed "//trim(fnam),io_stat)
    !         read(funit, rec=1) n
    !         allocate( arr(n), stat=alloc_stat )
    !         if(alloc_stat /= 0) call allocchk('In: file2iarr; simple_fileio  ', alloc_stat)
    !         do i=1,n
    !             read(funit, rec=i+1) arr(i)
    !         end do
    !         call fclose_1(funit,io_stat, errmsg="file2iarr failed to close "//trim(fnam))
    !     else
    !         write(*,*) fnam
    !         call simple_stop('file does not exist; file2iarr; simple_fileio      ', __FILENAME__,__LINE__)
    !     endif
    ! end function file2iarr

    !! Used once in reconstructor_eo
    !> \brief  for converting a real array 2 file
    subroutine arr2file_1( arr, fnam )
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
    end subroutine arr2file_1

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
            call simple_stop(trim(fnam)//' does not exist; file2rarr; simple_fileio ', __FILENAME__,__LINE__)
        endif
    end function file2rarr

    !> \brief  for converting an integer array 2 file
    ! subroutine arr2file_2( arr, fnam )
    !     integer,          intent(in) :: arr(:) !< array of data
    !     character(len=*), intent(in) :: fnam   !< output filename
    !     integer :: recsz, i, funit, ival, io_stat
    !     inquire(iolength=recsz) ival
    !     ival = size(arr)
    !     call fopen(funit,fnam,'replace','unknown', io_stat,'direct','unformatted',recl=recsz)
    !     call fileiochk("arr2file_2 fopen failed "//trim(fnam),io_stat)
    !     write(funit, rec=1) ival
    !     do i=1,size(arr)
    !         write(funit, rec=i+1) arr(i)
    !     end do
    !     call fclose_1(funit,io_stat, errmsg="arr2file_2 fclose failed "//trim(fnam))
    ! end subroutine arr2file_2

    !> \brief  for converting a real 2D array 2 file
    ! subroutine arr2D2file( arr, fnam )
    !     real,             intent(in) :: arr(:,:) !< array of data
    !     character(len=*), intent(in) :: fnam     !< output filename
    !     real    :: dim1, dim2
    !     integer :: funit, io_stat
    !     dim1 = real(size(arr,dim=1))
    !     dim2 = real(size(arr,dim=2))
    !     call fopen(funit,fnam,'replace','write', io_stat, 'STREAM')
    !     call fileiochk("simple_fileio::arr2D2file fopen failed "//trim(fnam),io_stat)
    !     write(unit=funit,pos=1,iostat=io_stat) dim1
    !     call fileiochk('simple_fileio::arr2D2file: writing stream startbyte 1 to: '//trim(fnam), io_stat)
    !     write(unit=funit,pos=5,iostat=io_stat) dim2
    !     call fileiochk('simple_fileio::arr2D2file: writing stream startbyte 5 to: '//trim(fnam), io_stat)
    !     write(unit=funit,pos=9,iostat=io_stat) arr(:,:)
    !     call fileiochk('simple_fileio::arr2D2file: writing stream startbyte 9 to: '//trim(fnam), io_stat)
    !     call fclose_1(funit,io_stat, errmsg=" arr2D2file Error closing file # "//int2str(funit))
    ! end subroutine arr2D2file

    !> \brief  for converting a real 2D array 2 file
    ! function file2arr2D( fname ) result( arr )
    !     character(len=*), intent(in) :: fname   !< input filename
    !     real, allocatable :: arr(:,:)          !< array of data
    !     real    :: dim1r, dim2r
    !     integer :: dim1, dim2, funit, io_stat
    !     call fopen(funit,fname,'old','read',io_stat,'STREAM')
    !     call fileiochk("simple_fileio::file2arr2D fopen failed "//trim(fname),io_stat)
    !     read(unit=funit,pos=1,iostat=io_stat) dim1r
    !     call fileiochk("simple_fileio::file2arr2D  reading stream startbyte 1 from: "// trim(fname), io_stat)
    !     read(unit=funit,pos=5,iostat=io_stat) dim2r
    !     call fileiochk("simple_fileio::file2arr2D  reading stream startbyte 5 from: "// trim(fname), io_stat)
    !     dim1 = nint(dim1r)
    !     dim2 = nint(dim2r)
    !     if( allocated(arr) ) deallocate(arr)
    !     allocate( arr(dim1,dim2), stat=alloc_stat )
    !     if(alloc_stat /= 0) call allocchk('In: simple_fileio:: file22Darr arr ', alloc_stat)
    !     read(unit=funit,pos=9,iostat=io_stat) arr(:,:)
    !     call fileiochk("simple_fileio::file2arr2D  reading stream startbyte 9 from: "// trim(fname), io_stat)
    !     call fclose_1(funit,io_stat, errmsg="Error closing file "//trim(fname))
    ! end function file2arr2D

    !> \brief  for converting a real array 2 file
    ! subroutine arr2txtfile( arr, fname )
    !     real,             intent(in) :: arr(:) !< array of data
    !     character(len=*), intent(in) :: fname !< output filename
    !     integer :: i, funit, io_stat
    !     call fopen(funit, fname,'REPLACE', 'write', io_stat)
    !     call fileiochk("simple_fileio ::arr2txtfile, tried to open file "//trim(fname), io_stat )
    !     do i=1,size(arr)
    !         write(funit,*) arr(i)
    !     end do
    !     call fclose_1(funit,io_stat, errmsg="Error closing file "//trim(fname))
    ! end subroutine arr2txtfile

    ! FILE-HANDLING JIFFYS

    !> \brief  for reading raw images using stream access
    ! subroutine read_raw_image( fname, mat, first_byte )
    !     character(len=*), intent(in)  :: fname
    !     double precision, intent(out) :: mat(:,:,:)
    !     integer,          intent(in)  :: first_byte
    !     integer :: filnum, io_stat
    !     character(len=100) :: io_message

    !     call fopen(filnum, fname, 'OLD', 'READ', io_stat, 'STREAM', convert='NATIVE')
    !     call fileiochk("Error opening file "//trim(fname) , io_stat)
    !     read(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
    !     ! Check the read was successful
    !     call fileiochk('simple_fileio::read_raw_image; reading '//trim(fname)//&
    !         ' IO error message was: '// trim(io_message),io_stat)
    !     call fclose_1(filnum, io_stat,errmsg="Error closing file "//trim(fname))
    ! end subroutine read_raw_image

    !> \brief  for writing raw images using stream access
    ! subroutine write_raw_image( fname, mat, first_byte )
    !     character(len=*), intent(in) :: fname
    !     real,             intent(in) :: mat(:,:,:)
    !     integer,          intent(in) :: first_byte
    !     integer :: filnum, io_stat
    !     character(len=100) :: io_message
    !     call fopen(filnum,fname, 'REPLACE', 'WRITE', io_stat, 'STREAM')
    !     call fileiochk("Error opening file "//trim(fname), io_stat )
    !     write(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
    !     ! Check the write was successful
    !     call fileiochk('simple_fileio::write_raw_image; writing '//trim(fname)//&
    !         ' IO error message was: '// trim(io_message),io_stat)
    !     call fclose_1(filnum, io_stat,errmsg="Error closing file "//trim(fname))
    ! end subroutine write_raw_image

    !! Unused
    !> From flibs file_list
    ! Copyright (c) 2008, Arjen Markus
    ! subroutine file_list_old( dir, list , suppress_errors, outfile)
    !     character(len=*), intent(in)            :: dir
    !     character(len=*), pointer, dimension(:) :: list
    !     logical, intent(in), optional           :: suppress_errors
    !     character(len=*), intent(inout), optional  :: outfile
    !     character(len=200)                      :: cmd,redirect,tmpfile
    !     character(len=1)                        :: line
    !     integer                                 :: luntmp
    !     integer                                 :: i
    !     integer                                 :: ierr
    !     integer                                 :: count
    !     redirect=" "
    !     if (present(suppress_errors))redirect="2>/dev/null "
    !     if (.not. present(outfile))then
    !         open( newunit = luntmp, status = 'scratch' )
    !         inquire( luntmp, name = tmpfile ) ! Hope this is okay
    !         close( luntmp )
    !     else
    !         tmpfile = trim(adjustl(outfile))
    !         if( file_exists(tmpfile) ) call del_file(tmpfile)
    !     end if
    !     cmd = 'ls -tr ' // ' ' // trim(dir) // ' ' // trim(redirect) // tmpfile
    !     call exec_cmdline( cmd )
    !     open( newunit = luntmp, file = tmpfile )
    !     !
    !     ! First count the number of files, then allocate and fill the array
    !     !
    !     do
    !         read( luntmp, '(a)', iostat = ierr ) line
    !         if ( ierr == 0 ) then
    !             count = count + 1
    !         else
    !             exit
    !         end if
    !     end do
    !     rewind( luntmp )
    !     allocate( list(count) )
    !     do i = 1,count
    !         read( luntmp, '(a)' ) list(i)
    !     end do
    !     close( luntmp, status = 'delete' )
    ! end subroutine file_list_old

    subroutine simple_copy_file(fname1, fname2, status)
        character(len=*), intent(in)           :: fname1, fname2 !< input filenames
        integer, intent(out), optional :: status
        character(len=STDLEN) :: cmd
        status=0
        cmd = 'cp '//trim(fname1)//'  '//trim(fname2)
        call exec_cmdline(cmd)
    end subroutine simple_copy_file

    subroutine ls_mrcfiletab( dir, filetabname )
        character(len=*),intent(in)  :: dir, filetabname
         character(len=STDLEN) :: cmd
         cmd = 'ls -tr '//trim(dir)//'/*.mrc*'//' > '//trim(filetabname)
         call exec_cmdline(cmd)

        ! integer :: stat
        ! stat = simple_glob_list_tofile(glob=trim(dir//'/*.mrc*'), outfile=trim(filetabname), tr=.true.)
        ! if(stat/=0) call fileiochk("ls_mrcfiletab failed "//trim(dir))
    end subroutine ls_mrcfiletab

    ! subroutine ls_fbody_mrcfiletab( fbody, filetabname )
    !     character(len=*),intent(in)  :: fbody, filetabname
    !     character(len=STDLEN) :: cmd
    !     cmd = 'ls -tr '//trim(fbody)//'*.mrc*'//' > '//trim(filetabname)
    !     print *,trim(cmd)
    !     call exec_cmdline(cmd)

    !     ! integer :: stat
    !     ! stat = simple_glob_list_tofile(glob=trim(fbody//'*.mrc*'), outfile=trim(filetabname), tr=.true.)
    !     ! if(stat/=0) call fileiochk("ls_fbody_mrcfiletab failed "//trim(fbody))
    ! end subroutine ls_fbody_mrcfiletab

    ! subroutine ls_filetab( fbody, ext, filetabname )
    !     character(len=*), intent(in)  :: fbody, ext, filetabname
    !     character(len=STDLEN) :: cmd
    !     cmd = 'ls -tr '//trim(fbody)//'*'//trim(ext)//' > '//trim(filetabname)
    !     call exec_cmdline(cmd)

    !     ! integer :: stat
    !     ! stat = simple_glob_list_tofile(glob=trim(fbody)//'*.'//trim(ext), outfile=trim(filetabname), tr=.true.)
    !     ! if(stat/=0) call fileiochk("ls_filetab failed "//trim(fbody))
    ! end subroutine ls_filetab

    ! subroutine sys_del_files( fbody, ext )
    !     character(len=*),      intent(in)  :: fbody, ext
    !     character(len=STDLEN), allocatable :: fnames(:)
    !     character(len=STDLEN), parameter   :: ftab = 'ftab_from_sys_del_files.txt'
    !     integer :: i, last
    !     call ls_filetab(fbody, ext, ftab) ! filetable written to disc
    !     call read_filetable(ftab, fnames) ! filetable read back in
    !     last = size(fnames)
    !     do i=1,last
    !         call del_file(fnames(i))
    !     end do
    !     call del_file(ftab)
    !     deallocate(fnames)
    ! end subroutine sys_del_files

    ! function get_last_fname( fbody, ext ) result( fname )
    !     character(len=*),      intent(in)  :: fbody, ext
    !     character(len=STDLEN), allocatable :: fnames(:)
    !     character(len=STDLEN), parameter   :: ftab = 'ftab_from_sys_find_last_fname.txt'
    !     character(len=STDLEN) :: fname
    !     integer :: last
    !     call ls_filetab(fbody, ext, ftab) ! filetable written to disc
    !     call read_filetable(ftab, fnames) ! filetable read back in
    !     last = size(fnames)
    !     fname = fnames(last)
    !     call del_file(ftab)
    !     deallocate(fnames)
    ! end function get_last_fname

    ! subroutine merge_docs( docnames, fname_merged )
    !     character(len=STDLEN), intent(in) :: docnames(:)
    !     character(len=*),      intent(in) :: fname_merged
    !     character(len=STDLEN) :: cmd
    !     integer :: ndocs, idoc
    !     call del_file(fname_merged)
    !     ndocs = size(docnames)
    !     do idoc=1,ndocs
    !         cmd = 'cat '//trim(docnames(idoc))//' >> '//trim(fname_merged)
    !         call exec_cmdline(cmd)
    !     end do
    ! end subroutine merge_docs

end module simple_fileio
