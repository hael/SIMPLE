! generic fileio       class

! #define FCLOSE(X,Y)  if(.not.file_close( X, Y )) &
! #define IOERRMSG(Y,Z) write(,'(/,A,/,":",I0,": IOSTAT ",I0," ",A)') __FILENAME__,__LINE__, Y,Z

module simple_fileio
    use simple_defs
    use simple_syslib
!use ISO_C_BINDING
use, intrinsic :: iso_fortran_env, only: stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT, stdin=>INPUT_UNIT
!use, intrinsic :: iso_fortran_env
use simple_strings, only: upperCase,stringsAreEqual, strIsBlank, int2str,int2str_pad
implicit none

interface arr2file
    module procedure arr2file_1
    module procedure arr2file_2
end interface arr2file


! interface
!     !> Declare the interface for POSIX fsync function
!     function fsync (fd) bind(c,name="fsync")
!     use iso_c_binding, only: c_int
!         integer(c_int), value :: fd
!         integer(c_int) :: fsync
!     end function fsync
! end interface
! private
! #include "simple_local_flags.inc"

!public :: alloc_errchk, err_msg, errcheck
!public :: fopen, fclose, nlines,

contains

    !> \brief  is for checking file IO status
    subroutine fileio_errmsg( message, iostat , die)
        use simple_syslib, only: simple_error_check
        character(len=*), intent(in) :: message  !< error message
        integer, intent(inout)          :: iostat !< error status
        logical, intent(in), optional   :: die    !< do you want to terminate or not
        integer :: errno
        logical :: this_die
        this_die=.true.
        if(present(die)) this_die=die
        if (iostat == -1)then
            write(stderr,'(a)') "fileio: EOF reached (PGI version)"
        else if (iostat == -2) then
            write(stderr,'(a)') "fileio: End-of-record reached (PGI version)"
        else if( iostat /= 0 ) then
            write(stderr,'(a,1x,I0 )') 'ERROR: File I/O failure, IOS# ',iostat
            call simple_error_check(iostat)
            if(this_die)stop
        endif
    end subroutine fileio_errmsg

    !> FOPEN enforce F2008 style open so that PGI/Intel behave correctly
    function fopen (funit,file,status,action,iostat,access,form,recl,async,pad,&
         decimal,round,delim,blank,convert,iomsg,position)
        integer, intent(inout) :: funit
        character(len=*), intent(in) :: file
        integer, intent(inout), optional :: iostat
        integer, intent(inout), optional :: recl
        character(len=*), intent(in), optional :: access, async, action, &
             &status, blank, pad, form, decimal, round, delim, convert, &
             &iomsg, position
        integer  :: iostat_this,  recl_this
        character(len=STDLEN) :: filename,iomsg_this
        character(len=30) :: async_this, access_this, action_this, status_this,&
             &blank_this, pad_this, decimal_this, delim_this, form_this , &
             &round_this,position_this
        logical :: fopen,is_open
        fopen=.false.
        ! check to see if filename is empty
        write(filename,'(A)') trim(adjustl(file))
        if ( strIsBlank(filename) )then
            print *, 'simple_system::fopen filename blank'
            if(present(iomsg))  print *, trim(adjustl(iomsg))
            return
        end if

        if (.not. (present(iostat) ) )then
            open(NEWUNIT=funit, FILE=trim(adjustl(filename)),IOSTAT=iostat_this)
            call simple_error_check(iostat_this,"fileio::fopen basic open "//trim(filename))
            fopen = .true.
            return ! true
        end if
        ! Optional args
        if(present(convert)) then
!            print *, 'CONVERT ignored in file open argument list'

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
        ! ACTION: READ, WRITE, or READWRITE.
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
                print *, "::fopen incompatible status=OLD and position=APPEND  when ", trim(filename)," does not exist"
                write( status_this,'(A)')  upperCase('NEW')
            end if

        end if
        ! access: DIRECT (random access) or SEQUENTIAL
        write(access_this ,'(A)') 'SEQUENTIAL'
        if (present(access))then
            if (stringsAreEqual(access, 'DIRECT',.false.))  write(access_this ,'(A)') upperCase(access)
            if (stringsAreEqual(access, 'STREAM',.false.))then
#ifdef PGI
                print *,"** Cannot 'stream' in current PGI version, using DIRECT"
                write(access_this,'(A)')'DIRECT'
#else
                write(access_this ,'(A)') upperCase(access)
#endif
            endif
        end if
        !! Common file open
        if (.not.( present(form) .or. present(recl) .or. present(async) .or. present(pad) .or. &
             &present(decimal) .or. present(round) .or. present(delim) .or. present(blank) ) )then
            if (present(position))then
                !! Appending to file
                open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                     &ACTION=action_this,STATUS=status_this,ACCESS=access_this,POSITION=position_this)
            else
                open( NEWUNIT=funit,FILE=filename,IOSTAT=iostat_this,&
                     &ACTION=action_this,STATUS=status_this,ACCESS=access_this)
            end if
            call simple_error_check(iostat_this,"::fopen common open ACTION/STATUS/ACCESS")
            fopen=.true.
            if(present(iostat))iostat=iostat_this
            return
        end if

        ! err_this=2000
       ! if(present(err))err_this=err
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
                     &ACCESS=access_this, FORM=form_this,&
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
        call simple_error_check(iostat_this,"fileio::fopen extra open ")
        if(present(iostat))iostat=iostat_this
        if(present(recl))recl=recl_this
        fopen=.true.
        return ! true
        ! Seriously bad hack for PGI system call
!91      print stderr, "fopen: ERR called file " ,trim(filename)," err # ",iostat_this
!        return
    end function fopen


    function fclose (unit,iostat,status,dispose)
        integer, intent(in) :: unit
        integer, intent(inout) :: iostat
        character(len=*), intent(in), optional :: status,dispose
        integer :: iostat_this
        character(len=30) :: status_this,dispose_this
        logical :: fclose
        fclose=.false.

        !! status or dispose: 'KEEP' or 'DELETE', unless file was opened with status=SCRATCH
        write(status_this,'(A)')'KEEP'
        if (present(dispose))then
            if (stringsAreEqual(dispose, 'DELETE',.false.))  write(status_this ,'(A)') upperCase(dispose)
        end if
        if (present(status))then
            if (stringsAreEqual(status, 'DELETE',.false.))  write(status_this ,'(A)') upperCase(status)
        end if
        if (is_open(unit)) then
            CLOSE (unit,IOSTAT=iostat,STATUS=status_this)
            call simple_error_check(iostat, "fclose failed ")
        end if
        fclose=.true.
        return
! 92      print ERROR_UNIT, "fclose failed ", int2str(iostat_this)
!         return
    end function fclose


    !> \brief return the number of lines in a textfile
    function nlines( fname ) result( n )
        character(len=*), intent(in) :: fname !< input filename
        character(len=:), allocatable   :: tfile
        integer          :: n, funit, ios,io_status
        character(len=1) :: junk
        if( file_exists(fname) )then
            tfile=fname
            if(.not.fopen(funit, tfile, status='unknown', action='read', iostat=io_status))then
                call fileio_errmsg(":nlines error opening file "//trim(tfile), io_status)
            end if
            n = 0
            do
                 read(funit,*,IOSTAT=ios) junk
                 if(ios /= 0)then
                     exit
                 else
                     n = n + 1
                 endif
            end do
            if(.not.fclose( funit, io_status )) call fileio_errmsg(" Error closing file in ::nlines ",io_status)
        else
            n = 0
        endif
    end function nlines

    !> \brief  return the size of a binary file
    function filelength( fname ) result( filesz )
        character(len=*), intent(inout) :: fname !< input filename
        integer                      :: filesz, funit, ios, cnt,recl
        character(len=1)             :: junk
        if(  file_exists(fname) )then
            recl=1
            if(.not.fopen(funit, fname, ACTION='read', IOSTAT=ios,&
            ACCESS='direct', form='unformatted', recl=recl))then
                call fileio_errmsg('simple_fileio       :: nlines', ios)
            end if
            !open( unit=funit, file=trim(fname), action='read',&
            !access='direct', form='unformatted', recl=1 )
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
            if(.not.fclose( funit, ios )) &
                 call fileio_errmsg(" Error closing file in ::filelength ",ios)
        else
            filesz = 0
        endif
    end function filelength

    !> \brief  return the record size of a binary file
    ! function reclength( fname, nentries ) result( recsz )
    !     character(len=*), intent(in) :: fname     !< input filename
    !     integer,          intent(in) :: nentries  !< total num of entries
    !     integer                      :: recsz
    !     recsz = filelength(fname)/nentries
    ! end function reclength

    !> \brief  return file size in bytes
    function funit_size(unit) result(sz)
        integer, intent(in)          :: unit !< input file unit
        integer(kind=8)              :: sz
        inquire(unit=unit,size=sz)
    end function funit_size

    !> \brief  return file size in bytes
    function file_size(fname) result(sz)
        character(len=*), intent(in) :: fname !< input filename
        integer(kind=8)              :: sz
        inquire(file=trim(adjustl(fname)),size=sz)
    end function file_size

    !> \brief  is for deleting a file
    subroutine del_file( file )
        character(len=*), intent(in) :: file !< input filename
        integer :: fnr, ios

        if( file_exists(file) )then
            !call cpStr(file,tfile)
            if(.not.fopen(fnr,file,STATUS='OLD',IOSTAT=ios))&
                call fileio_errmsg( "del_file failed to open file designated for deletion", ios)
            if( ios == 0 )then
                if(.not.fclose(fnr, status='delete',iostat=ios))&
                     call fileio_errmsg("::fclose failed in del_file ",ios)
            end if
        endif
    end subroutine del_file

    !> \brief  is for deleting consequtively numbered files with padded number strings
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

    !> \brief  is for making a file-table (to be able to commander execute programs that depend on them)
    subroutine make_filetable( tabname, n, body, ext, numlen, suffix )
        character(len=*),           intent(in) :: tabname !< file-table (string)
        integer,                    intent(in) :: n       !< total num of files
        character(len=*),           intent(in) :: body    !< filename body
        character(len=*),           intent(in) :: ext     !< filename extension
        integer,          optional, intent(in) :: numlen  !< number length
        character(len=*), optional, intent(in) :: suffix  !< file suffix
        character(len=STDLEN), allocatable :: names(:)
        integer :: ifile, fnr, ios

        names = make_filenames( body, n, ext, numlen=numlen, suffix=suffix )
        if(.not.fopen(fnr, file=tabname, status='replace', action='write', iostat=ios))then
            call fileio_errmsg('simple_fileio       :: make_filetable', ios)
        end if
        do ifile=1,n
            write(fnr,'(a)') trim(names(ifile))
        end do
        if(.not.fclose( fnr, ios )) call fileio_errmsg(" Error closing file in ::make_filetable ",ios)
    end subroutine make_filetable

    !> \brief  is for making a file-table (to be able to commander execute programs that depend on them)
    !! \param tab1,tab2,tab3,tab4 input filenames
    subroutine make_multitab_filetable( tabname, tab1, tab2, tab3, tab4 )
        character(len=*),                intent(in) :: tabname
        character(len=STDLEN),           intent(in) :: tab1(:), tab2(:)
        character(len=STDLEN), optional, intent(in) :: tab3(:), tab4(:)
        integer :: ntabs, n, fnr, ifile, ios

        ntabs = 2
        if( present(tab3) ) ntabs = 3
        if( present(tab4) ) ntabs = 4
        n = size(tab1)

        if(.not.fopen(fnr,tabname, status='replace', action='write', iostat=ios))&
            call fileio_errmsg('simple_fileio       :: make_multitab_filetable', ios)

        do ifile=1,n
            if( ntabs == 2 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))
            if( ntabs == 3 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))&
            &//' '//trim(tab3(ifile))
            if( ntabs == 4 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))&
            &//' '//trim(tab3(ifile))//' '//trim(tab4(ifile))
        end do
        if(.not.fclose( fnr, ios )) call fileio_errmsg(" Error closing file in ::make_multitab_filetable ",ios)
    end subroutine make_multitab_filetable

    !> \brief  is for checking file kind
    !> \param fname,suffix string args to check suffix
    function file_kind( fname, suffix ) result( yep )
        character(len=*), intent(in) :: fname, suffix
        integer :: pos
        logical :: yep
        pos = index(fname, suffix) ! position of suffix
        if( pos == 0 )then
            yep = .false.
        else
            yep = .true.
        endif
    end function file_kind



    !>  Wrapper for POSIX system call stat
    subroutine sys_stat( filename, status, buffer, doprint )
#if defined(INTEL)
        use ifport
        use ifposix
#endif
        character(len=*),     intent(in) :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        logical, optional,    intent(in)    :: doprint
        logical :: l_print = .true.
#if defined(GNU)
        allocate(buffer(13), source=0)
        call stat(trim(adjustl(filename)), buffer, status)
#elif defined(PGI)
        integer,allocatable :: statb(:)
        integer             :: stato
#include "simple_local_flags.inc"
        include 'lib3f.h'
        status =  stat(trim(adjustl(filename)), buffer)
!        DebugPrint 'fileio       sys_stat PGI stato ', status
!        DebugPrint 'fileio       sys_stat PGI size of buffer ', size(statb)

#elif defined(INTEL)
        ! use ifport
        ! integer(4) statarray(12), istat
        ! ISTAT = STAT (1, statarray)
        ! if (.NOT. istat) then
        !     print *, statarray
        ! end if

        integer(4) :: ierror
        integer(8) :: jhandle
        allocate(buffer(13), source=0)
        call pxfstat (trim(adjustl(filename)), len_trim(adjustl(filename)), jhandle, ierror)
        call PXFSTRUCTCREATE('stat', jhandle, ierror)
        if(ierror.EQ.0) then
            if(ierror.EQ.0) then
                CALL PXFINTGET (jhandle,'st_dev',buffer(1), ierror)    ! Device ID
                CALL PXFINTGET (jhandle,'st_ino',buffer(2), ierror)    ! Inode number
                CALL PXFINTGET (jhandle,'st_mode' , buffer(3), ierror) !  File mode
                CALL PXFINTGET (jhandle,'st_nlink' ,buffer(4), ierror) ! Number of links
                CALL PXFINTGET (jhandle,'st_uid' ,buffer(5), ierror)   ! Owner’s uid
                CALL PXFINTGET (jhandle,'st_gid' ,buffer(6), ierror)   ! Owner’s gid
                buffer(7)=0 ! ID of device containing directory entry for file (0 if not available)
                CALL PXFINTGET (jhandle,'st_size',buffer(8), ierror)   ! File size (bytes)
                CALL PXFINTGET (jhandle,'st_atime',buffer(9), ierror)  ! Last access time
                CALL PXFINTGET (jhandle,'st_mtime',buffer(10), ierror) ! Last modification time
                CALL PXFINTGET (jhandle,'st_ctime',buffer(11), ierror) ! Last file status change time
                buffer(12)=0 ! Preferred I/O block size (-1 if not available)
                buffer(13)=0 ! Number of blocks allocated (-1 if not available)
                call PXFSTRUCTFREE(jhandle,ierror)
            else
                print *, 'Filehandling sys_stat PXFSTAT failed, file ', trim(adjustl(filename)),' error ', ierror
            end if
            if (ierror.NE.0)then
                print *, 'Filehandling sys_stat PXFINTGET failed, file ', trim(adjustl(filename)),' error ', ierror
            end if
        else
            call PXFSTRUCTFREE(jhandle,ierror)
            stop 'Filehandling sys_stat  failed - cannot create structure for jhandle1'
        end if
        status=ierror
#endif
        if( present(doprint) )l_print = doprint
        if( l_print )then
            write(*,*) 'command: stat ', trim(adjustl(filename))
            write(*,*) 'status of execution: ', status
        endif
    end subroutine sys_stat

    !> \brief  is for adding to filebody
    function add2fbody( fname, suffix, str ) result( newname )
        character(len=*), intent(in)  :: fname, suffix, str
        character(len=:), allocatable :: newname
        integer :: pos
        pos = index(fname, suffix) ! position of suffix
        allocate(newname, source=fname(:pos-1)//trim(str)//trim(suffix))
    end function add2fbody

    !> \brief  is for extracting the body of a file
    function get_fbody( fname, suffix ) result( fbody )
        character(len=*), intent(in) :: fname, suffix !< file extension
        character(len=STDLEN)        :: fbody
        integer :: pos
        pos = index(fname, '.'//suffix) ! position of suffix
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

    pure function remove_abspath( fname ) result( new_fname)
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
    end function remove_abspath

    pure function extract_abspath( fname ) result( abspath )
        character(len=*), intent(in)  :: fname !< abs filename
        character(len=:), allocatable :: abspath !< abs file path
        integer :: length, pos
        length = len_trim(fname)
        pos = scan(fname(1:length),'/',back=.true.)
        allocate(abspath, source=trim(fname(1:pos)))
    end function extract_abspath

    !>  \brief  returns the integer number identifier of a filename
    subroutine fname2ind( str, ivar )
        use simple_strings, only: map_str_nrs, str2int
        character(len=*), intent(in)  :: str    !< abs filename
        integer,          intent(out) :: ivar   !< file index number
        logical, allocatable          :: pos(:)
        character(len=:), allocatable :: str_copy
        integer :: j, lstr, io_stat, nrrange(2)
        lstr = len(str);  nrrange = 0
        pos = map_str_nrs(str)
        if( any(pos) )then
            do j=lstr,1,-1
                if( pos(j) )then
                    nrrange(1) = j
                    nrrange(2) = j
                    do while( pos(nrrange(1)) )
                        nrrange(1) = nrrange(1)-1
                    end do
                    nrrange(1) = nrrange(1)+1
                    exit
                endif
            end do
            allocate(str_copy, source=str(nrrange(1):nrrange(2)))
            call str2int(str_copy, io_stat, ivar)
        else
            allocate(str_copy, source='1')
            call str2int(str_copy, io_stat, ivar)
        endif
        if( allocated(pos)      ) deallocate(pos)
        if( allocated(str_copy) ) deallocate(str_copy)
    end subroutine fname2ind

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
        character(len=STDLEN), allocatable :: names(:)
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
            case ('txt', 'asc', 'box','dat')
                fname2format = 'T'
            case DEFAULT
                fname2format = 'N'
        end select
    end function fname2format

    !>  \brief  to check if same file format
    !! \param fname1,fname2 input filenames
    pure logical function same_format( fname1, fname2 )
        character(len=*), intent(in) :: fname1, fname2
        character(len=1) :: form1, form2
        form1 = fname2format(fname1)
        form2 = fname2format(fname2)
        same_format = form1 == form2
    end function same_format

    !>  \brief  reads a filetable into an array
    subroutine read_filetable( filetable, filenames )
        character(len=*),                   intent(in)  :: filetable    !< input table filename
        character(len=STDLEN), allocatable, intent(out) :: filenames(:) !< array of filenames
        integer :: nl, funit, alloc_stat, iline,io_stat

        nl    = nlines(filetable)
        if(.not.fopen(funit,filetable,'old','unknown',io_stat))&
            call fileio_errmsg("read_filetable failed to open file "//filetable,io_stat )
        allocate( filenames(nl), stat=alloc_stat )
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: read_filetable; simple_fileio      '
            stop
        endif
        do iline=1,nl
            read(funit,'(a256)') filenames(iline)
        end do
        if(.not.fclose(funit,io_stat)) call fileio_errmsg("read_filetable failed to close",io_stat)
    end subroutine read_filetable

    !>  \brief  writes a filetable array to a text file
    subroutine write_filetable( filetable, filenames )
        character(len=*),      intent(in)  :: filetable  !< output table filename
        character(len=STDLEN), intent(in) :: filenames(:)!< array of filenames
        integer :: nl, funit, iline, io_stat
        nl = size(filenames)
        if(.not.fopen(funit,filetable, 'replace', 'unknown', io_stat))then
            call fileio_errmsg("write_filetable failed to open file "//filetable,io_stat )
        end if
        do iline=1,nl
            write(funit,'(a)') trim(filenames(iline))
        end do
        if(.not.fclose(funit,io_stat)) call fileio_errmsg("write_filetable failed to close",io_stat)
    end subroutine write_filetable

    !> \brief  for converting a file generated by txtfile2arr back to an array
    function txtfile2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam    !< input table filename
        real, allocatable :: arr(:)             !< array of filenames
        integer :: i, n, alloc_stat, funit, io_stat
       
        if( file_exists(trim(fnam)) )then
            n = nlines(fnam)
            allocate( arr(n), stat=alloc_stat )
            if( alloc_stat /= 0 ) then
                write(*,'(a)') 'ERROR: Allocation failure!'
                write(*,'(a)') 'In: file2arr; simple_fileio      '
                stop
            endif
            if(.not.fopen(funit,fnam,'old','unknown',io_stat))&
                 call fileio_errmsg("txtfile2rarr failed to open  "//trim(fnam), io_stat)
            do i=1,n
                read(funit,*) arr(i)
            end do
            if(.not.fclose(funit,io_stat)) call fileio_errmsg("txtfile2rarr failed to close  "//trim(fnam), io_stat)
        else
            write(*,*) fnam
            stop 'file does not exist; txtfile2rarr; simple_fileio      '
        endif
    end function txtfile2rarr

    !> \brief  merging two text files into a single array
    !! \param file1,file2 input filenames for merging
    function merge_txtfiles( file1, file2 )  result( arr )
        character(len=*), intent(in) :: file1, file2
        character(len=STDLEN), allocatable :: arr(:)
        integer :: n1, n2, alloc_stat, cnt, funit, i, io_stat
        logical :: here(2)
        here(1) = file_exists(trim(file1))
        here(2) = file_exists(trim(file2))
        n1 = 0
        n2 = 0
        if( file_exists(trim(file1)) ) n1 = nlines(file1)
        if( file_exists(trim(file2)) ) n2 = nlines(file2)
        allocate( arr(n1+n2), stat=alloc_stat )
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: merge_txtfiles; simple_fileio      '
            stop
        endif
        if( here(1) )then
            if(.not.fopen(funit,file1,'old','unknown',io_stat)) &
                 call fileio_errmsg("merge_txtfiles failed "//trim(file1), io_stat)
            cnt = 0
            do i=1,n1
                cnt = cnt+1
                read(funit,*) arr(cnt)
            end do
            if(.not.fclose(funit,io_stat)) &
                 call fileio_errmsg("merge_txtfiles failed to close "//trim(file1), io_stat)
            if( .not. here(2) ) return
        else
            if(.not.fopen(funit,file2,'old','unknown',io_stat))&
                 call fileio_errmsg("merge_txtfiles failed to open "//trim(file2), io_stat)
            cnt = 0
            do i=1,n2
                cnt = cnt+1
                read(funit,*) arr(cnt)
            end do
            if(.not.fclose(funit,io_stat))&
                 call fileio_errmsg("merge_txtfiles failed to close "//trim(file2), io_stat)
            return
        endif
        if(.not.fopen(funit,file2,'old','unknown',io_stat))&
             call fileio_errmsg("merge_txtfiles failed to open "//trim(file2), io_stat)
        do i=1,n2
            cnt = cnt+1
            read(funit,*) arr(cnt)
        end do
        if(.not.fclose(funit,io_stat))&
             call fileio_errmsg("merge_txtfiles failed to close "//trim(file2), io_stat)
    end function merge_txtfiles

    !> \brief  for converting a file generated by file2arr back to an array
    function file2iarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam             !< input table filename
        integer, allocatable :: arr(:)                   !< array of filenames
        integer :: recsz, i, n, alloc_stat, funit, ival, io_stat
        if( file_exists(trim(fnam)) )then
            inquire(iolength=recsz) ival
            if(.not.fopen(funit,fnam,'OLD','unknown', io_stat,'direct','unformatted',recsz))&
                call fileio_errmsg("file2iarr fopen failed "//trim(fnam),io_stat)
            read(funit, rec=1) n
            allocate( arr(n), stat=alloc_stat )
            if( alloc_stat /= 0 ) then
                write(*,'(a)') 'ERROR: Allocation failure!'
                write(*,'(a)') 'In: file2arr; simple_fileio      '
                stop
            endif
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            if(.not.fclose(funit,io_stat))&
                 call fileio_errmsg("file2iarr failed to close "//trim(fnam), io_stat)
        else
            write(*,*) fnam
            stop 'file does not exist; file2iarr; simple_fileio      '
        endif
    end function file2iarr

    !> \brief  for converting a real array 2 file
    subroutine arr2file_1( arr, fnam )
        real,             intent(in) :: arr(:)    !< array of filenames
        character(len=*), intent(in) :: fnam      !< input table filename
        real    :: rval
        integer :: recsz, i, funit,io_stat
        inquire(iolength=recsz) rval
        rval = size(arr)
        if(.not.fopen(funit,fnam,'replace','unknown', io_stat,'direct','unformatted',recl=recsz))then
            call fileio_errmsg("arr2file_1 fopen failed "//trim(fnam),io_stat)
        endif
        write(funit, rec=1) rval
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        if(.not.fclose(funit,io_stat))&
             call fileio_errmsg("arr2file_1 fclose failed "//trim(fnam),io_stat)
    end subroutine arr2file_1

    !> \brief  for converting a file generated by arr2file back to an array
    function file2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam  !< input table filename
        real, allocatable            :: arr(:) !< array of filenames
        real    :: rval
        integer :: recsz, i, n, alloc_stat, funit,io_stat
        if( file_exists(trim(fnam)) )then
            inquire(iolength=recsz) rval
            if(.not.fopen(funit,fnam,'old','unknown', io_stat,'direct','unformatted',recl=recsz))then
                call fileio_errmsg("file2rarr fopen failed "//trim(fnam),io_stat)
            endif
            read(funit, rec=1) rval
            n = nint(rval)
            allocate( arr(n), stat=alloc_stat )
            if( alloc_stat /= 0 ) then
                write(*,'(a)') 'ERROR: Allocation failure!'
                write(*,'(a)') 'In: file2arr; simple_fileio      '
                stop
            endif
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            if(.not.fclose(funit,io_stat))&
                call fileio_errmsg("file2rarr fclose failed "//trim(fnam),io_stat)
        else
            write(*,*) fnam, ' does not exist; file2rarr; simple_fileio      '
            stop
        endif
    end function file2rarr

    !> \brief  for converting an integer array 2 file
    subroutine arr2file_2( arr, fnam )
        integer,          intent(in) :: arr(:)!< array of data
        character(len=*), intent(in) :: fnam !< output filename
        integer :: recsz, i, funit, ival, io_stat
        inquire(iolength=recsz) ival
        ival = size(arr)
        if(.not.fopen(funit,fnam,'replace','unknown', io_stat,'direct','unformatted',recl=recsz))then
            call fileio_errmsg("arr2file_2 fopen failed "//trim(fnam),io_stat)
        endif
        write(funit, rec=1) ival
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        if(.not.fclose(funit,io_stat)) &
             call fileio_errmsg("arr2file_2 fopen failed "//trim(fnam),io_stat)
    end subroutine arr2file_2

    !> \brief  for converting a real 2D array 2 file
    subroutine arr2D2file( arr, fnam )
        real,             intent(in) :: arr(:,:) !< array of data
        character(len=*), intent(in) :: fnam     !< output filename
        real    :: dim1, dim2
        integer :: funit, io_stat
        dim1 = real(size(arr,dim=1))
        dim2 = real(size(arr,dim=2))
        if(.not.fopen(funit,fnam,'replace','write', io_stat, 'STREAM'))then
            call fileio_errmsg("arr2D2file fopen failed "//trim(fnam),io_stat)
        endif
        write(unit=funit,pos=1,iostat=io_stat) dim1
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when writing stream startbyte 1 to: ', trim(fnam)
            stop 'I/O error; arr2D2file; simple_fileio      '
        endif
        write(unit=funit,pos=5,iostat=io_stat) dim2
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when writing stream startbyte 5 to: ', trim(fnam)
            stop 'I/O error; arr2D2file; simple_fileio      '
        endif
        write(unit=funit,pos=9,iostat=io_stat) arr(:,:)
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when writing stream startbyte 9 to: ', trim(fnam)
            stop 'I/O error; arr2D2file; simple_fileio      '
        endif
        if(.not.fclose(funit,io_stat)) &
             call fileio_errmsg("Error closing file # "//int2str(funit), io_stat )
    end subroutine arr2D2file

    !> \brief  for converting a real 2D array 2 file
    function file2arr2D( fname ) result( arr )
        character(len=*), intent(in) :: fname   !< input filename
        real, allocatable :: arr(:,:)          !< array of data
        real    :: dim1r, dim2r
        integer :: dim1, dim2, funit, io_stat, alloc_stat

        if(.not.fopen(funit,fname,'old','read',io_stat,'STREAM'))then
            call fileio_errmsg("file2arr2D fopen failed "//trim(fname),io_stat)
        endif
        read(unit=funit,pos=1,iostat=io_stat) dim1r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 1 from: ', trim(fname)
            stop 'I/O error; file22Darr; simple_fileio      '
        endif
        read(unit=funit,pos=5,iostat=io_stat) dim2r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 5 from: ', trim(fname)
            stop 'I/O error; file22Darr; simple_fileio      '
        endif
        dim1 = nint(dim1r)
        dim2 = nint(dim2r)
        if( allocated(arr) ) deallocate(arr)
        allocate( arr(dim1,dim2), stat=alloc_stat )
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: simple_fileio       :: file22Darr'
            stop
        endif
        read(unit=funit,pos=9,iostat=io_stat) arr(:,:)
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 9 from: ', trim(fname)
            stop 'I/O error; file22Darr; simple_fileio      '
        endif
        if(.not.fclose(funit,io_stat)) call fileio_errmsg("Error closing file "//trim(fname),io_stat)
    end function file2arr2D

    !> \brief  for converting a real array 2 file
    subroutine arr2txtfile( arr, fname )
        real,             intent(in) :: arr(:) !< array of data
        character(len=*), intent(in) :: fname !< output filename
        integer :: i, funit, io_stat

        if(.not.fopen(funit, fname,'REPLACE', 'write', io_stat))then
            call fileio_errmsg("simple_fileio       ::arr2txtfile, tried to open file "//trim(fname), io_stat )
        end if
        do i=1,size(arr)
            write(funit,*) arr(i)
        end do
       if(.not.fclose(funit,io_stat)) call fileio_errmsg("Error closing file "//trim(fname),io_stat)
    end subroutine arr2txtfile



    ! FILE-HANDLING JIFFYS


    !> \brief  for reading raw images using stream access
    subroutine read_raw_image( fname, mat, first_byte )
        character(len=*), intent(in)  :: fname
        double precision, intent(out) :: mat(:,:,:)
        integer, intent(in)           :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message

        if(.not.fopen(filnum, fname, 'OLD', 'READ', io_stat, 'STREAM', convert='NATIVE'))then
            call fileio_errmsg("Error opening file "//trim(fname) , io_stat)
        end if
        read(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the read was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        if(.not.fclose(filnum, io_stat)) call fileio_errmsg("Error closing file "//trim(fname),io_stat)
    end subroutine read_raw_image

    !> \brief  for writing raw images using stream access
    subroutine write_raw_image( fname, mat, first_byte )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: mat(:,:,:)
        integer,          intent(in) :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message

        if(.not.fopen(filnum,fname, 'REPLACE', 'WRITE', io_stat, 'STREAM'))then
            call fileio_errmsg("Error opening file "//trim(fname), io_stat )
        end if
        write(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        if(.not.fclose(filnum, io_stat))  call fileio_errmsg("Error closing file "//trim(fname),io_stat)
    end subroutine write_raw_image



    subroutine ls_mrcfiletab( dir, filetabname )
        character(len=*),intent(in)  :: dir, filetabname
        character(len=STDLEN) :: cmd
        cmd = 'ls -tr '//trim(dir)//'/*.mrc*'//' > '//trim(filetabname)
        call exec_cmdline(cmd)
    end subroutine ls_mrcfiletab

    subroutine ls_filetab( fbody, ext, filetabname )
        character(len=*), intent(in)  :: fbody, ext, filetabname
        character(len=STDLEN) :: cmd
        cmd = 'ls -tr '//trim(fbody)//'*'//trim(ext)//' > '//trim(filetabname)
        call exec_cmdline(cmd)
    end subroutine ls_filetab

    subroutine sys_del_files( fbody, ext )
        character(len=*),      intent(in)  :: fbody, ext
        character(len=STDLEN), allocatable :: fnames(:)
        character(len=STDLEN), parameter   :: ftab = 'ftab_from_sys_del_files.txt'
        integer :: i, last
        call ls_filetab(fbody, ext, ftab) ! filetable written to disc
        call read_filetable(ftab, fnames)      ! filetable read back in
        last = size(fnames)
        do i=1,last
            call del_file(fnames(i))
        end do
        call del_file(ftab)
        deallocate(fnames)
    end subroutine sys_del_files

    function get_last_fname( fbody, ext ) result( fname )
        character(len=*),      intent(in)  :: fbody, ext
        character(len=STDLEN), allocatable :: fnames(:)
        character(len=STDLEN), parameter   :: ftab = 'ftab_from_sys_find_last_fname.txt'
        character(len=STDLEN) :: fname
        integer :: last
        call ls_filetab(fbody, ext, ftab) ! filetable written to disc
        call read_filetable(ftab, fnames)      ! filetable read back in
        last = size(fnames)                    
        fname = fnames(last)
        call del_file(ftab)
        deallocate(fnames)
    end function get_last_fname

    subroutine merge_docs( docnames, fname_merged )
        character(len=STDLEN), intent(in) :: docnames(:)
        character(len=*),      intent(in) :: fname_merged
        character(len=STDLEN) :: cmd
        integer :: ndocs, idoc
        call del_file(fname_merged)
        ndocs = size(docnames)
        do idoc=1,ndocs
            cmd = 'cat '//trim(docnames(idoc))//' >> '//trim(fname_merged)
            call exec_cmdline(cmd)
        end do
    end subroutine merge_docs


end module simple_fileio
