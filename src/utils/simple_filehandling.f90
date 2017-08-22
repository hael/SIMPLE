! generic fileio       class
module simple_fileio      
use simple_defs
implicit none

interface arr2file
    module procedure arr2file_1
    module procedure arr2file_2
end interface arr2file

interface
    !> Declare the interface for POSIX fsync function
    function fsync (fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
        integer(c_int), value :: fd
        integer(c_int) :: fsync
    end function fsync
end interface

contains

    !> \brief return the number of lines in a textfile
    function nlines( fname ) result( n )
        character(len=*), intent(in) :: fname !< input filename
        integer          :: n, funit, ios
        logical          :: here
        character(len=1) :: junk
        inquire(FILE=fname, EXIST=here)
        if( here )then
            funit = get_fileunit( )
            open( unit=funit, file=trim(fname) )
            n = 0
            do
                read(funit,*,IOSTAT=ios) junk
                if(ios /= 0)then
                    exit
                else
                    n = n + 1
                endif
            end do
            close( unit=funit )
        else
            n = 0
        endif
    end function nlines

    !> \brief  return the size of a binary file
    function filelength( fname ) result( filesz )
        character(len=*), intent(in) :: fname !< input filename
        integer                      :: filesz, funit, ios, cnt
        logical                      :: here
        character(len=1)             :: junk
        inquire(FILE=fname, EXIST=here)
        if( here )then
            funit = get_fileunit( )
            open( unit=funit, file=trim(fname), action='read',&
            access='direct', form='unformatted', recl=1 )
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
            close( unit=funit )
        else
            filesz = 0
        endif
    end function filelength

    !> \brief  return the record size of a binary file
    function reclength( fname, nentries ) result( recsz )
        character(len=*), intent(in) :: fname     !< input filename
        integer,          intent(in) :: nentries  !< total num of entries
        integer                      :: recsz
        recsz = filelength(fname)/nentries
    end function reclength

    !> \brief  return file size in bytes
    function file_size(fname) result(sz)
        character(len=*), intent(in) :: fname !< input filename
        integer(kind=8)              :: sz
        inquire(file=trim(adjustl(fname)),size=sz)
    end function file_size

    !> \brief  is for deleting a file
    subroutine del_file( file )
        character(len=*), intent(in) :: file !< input filename
        integer :: fnr, file_status
        if( file_exists(file) )then
            fnr = get_fileunit()
            open(unit=fnr, iostat=file_status, file=trim(file), status='old')
            if( file_status == 0 ) close(fnr, status='delete')
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
        integer :: ifile, fnr, file_status
        names = make_filenames( body, n, ext, numlen=numlen, suffix=suffix )
        fnr = get_fileunit()
        open(unit=fnr, status='replace', action='write', file=tabname, iostat=file_status)
        call fopen_err('simple_fileio       :: make_filetable', file_status)
        do ifile=1,n
            write(fnr,'(a)') trim(names(ifile))
        end do
        close(unit=fnr)
    end subroutine make_filetable

    !> \brief  is for making a file-table (to be able to commander execute programs that depend on them)
    !! \param tab1,tab2,tab3,tab4 input filenames
    subroutine make_multitab_filetable( tabname, tab1, tab2, tab3, tab4 )
        use simple_strings, only: int2str, int2str_pad
        character(len=*),                intent(in) :: tabname
        character(len=STDLEN),           intent(in) :: tab1(:), tab2(:)
        character(len=STDLEN), optional, intent(in) :: tab3(:), tab4(:)
        integer :: ntabs, n, fnr, ifile, file_status
        ntabs = 2
        if( present(tab3) ) ntabs = 3
        if( present(tab4) ) ntabs = 4
        n = size(tab1)
        fnr = get_fileunit()
        open(unit=fnr, status='replace', action='write', file=tabname, iostat=file_status)
        call fopen_err('simple_fileio       :: make_filetable', file_status)
        do ifile=1,n
            if( ntabs == 2 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))
            if( ntabs == 3 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))&
            &//' '//trim(tab3(ifile))
            if( ntabs == 4 ) write(fnr,'(a)') trim(tab1(ifile))//' '//trim(tab2(ifile))&
            &//' '//trim(tab3(ifile))//' '//trim(tab4(ifile))
        end do
        close(unit=fnr)
    end subroutine make_multitab_filetable

    !> \brief  is for getting a free fileunit
    function get_fileunit( ) result( iunit )
        integer :: i, ios, iunit
        logical :: lopen
        iunit = 0
        do i=1,99
            if (i /= 5 .and. i /= 6 .and. i /= 9)then
                inquire(unit=i, opened=lopen, iostat=ios)
                if( ios == 0 ) then
                    if( .not. lopen ) then
                        iunit = i
                        return
                    end if
                end if
            end if
        end do
    end function get_fileunit

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

    !> \brief  is for checking file open status
    subroutine fopen_err( message, file_status )
        character(len=*), intent(in) :: message  !< error message
        integer, intent(in)          :: file_status!< error status
        if( file_status /= 0 ) then
            write(*,'(a)') 'ERROR: File opening failure!'
            write(*,'(a)') message
            stop
        endif
    end subroutine fopen_err

    !>  \brief  check if a file exists on disk
    logical function file_exists(fname)
        character(len=*), intent(in) :: fname
        inquire(file=trim(adjustl(fname)), exist=file_exists)
    end function file_exists

    !>  \brief  check whether a IO unit is currently opened
    logical function is_open( unit_number )
        integer, intent(in)   :: unit_number
        integer               :: io_status
        character(len=STDLEN) :: io_message
        io_status = 0
        inquire(unit=unit_number, opened=is_open,iostat=io_status,iomsg=io_message)
        if (io_status .ne. 0) then
            print *, 'is_open: IO error ', io_status, ': ', trim(adjustl(io_message))
            stop 'IO error; is_open; simple_fileio      '
        endif
    end function is_open

    !>  \brief  check whether a file is currently opened
    logical function is_file_open( fname )
        character(len=*), intent(in)  :: fname
        integer               :: io_status
        character(len=STDLEN) :: io_message
        io_status = 0
        inquire(file=fname, opened=is_file_open,iostat=io_status,iomsg=io_message)
        if (io_status .ne. 0) then
            print *, 'is_open: IO error ', io_status, ': ', trim(adjustl(io_message))
            stop 'IO error; is_file_open; simple_fileio      '
        endif
    end function is_file_open

    !>  Wrapper for system call
    subroutine sys_stat( filename, status, buffer, doprint )
#if defined(INTEL)
        use ifport
        use ifposix
        character(len=*),     intent(inout) :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        logical, optional,    intent(in)    :: doprint
        logical    :: l_print = .true.
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
#elif defined(PGI)
        include 'lib3f.h'
        character(len=*),     intent(inout) :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        logical, optional,    intent(in)    :: doprint
        integer,allocatable :: statb(:)
        integer             :: stato
        logical             :: l_print = .true.
        stato =  stat(trim(adjustl(filename)), statb)
        status = stato
        buffer = statb
#ifdef _DEBUG
        print *, 'fileio       sys_stat PGI stato ', stato
        print *, 'fileio       sys_stat PGI size of statb ', size(statb)
#endif

#elif defined(GNU)
        character(len=*),     intent(inout) :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        logical, optional,    intent(in)    :: doprint
        logical :: l_print = .true.
        allocate(buffer(13), source=0)
        call stat(trim(adjustl(filename)), buffer, status)
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
        integer :: nl, funit, alloc_stat, iline
        nl    = nlines(filetable)
        funit = get_fileunit()
        open(unit=funit, status='old', file=filetable)
        allocate( filenames(nl), stat=alloc_stat )
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: read_filetable; simple_fileio      '
            stop
        endif
        do iline=1,nl
            read(funit,'(a256)') filenames(iline)
        end do
        close(funit)
    end subroutine read_filetable

    !>  \brief  writes a filetable array to a text file
    subroutine write_filetable( filetable, filenames )
        character(len=*),      intent(in)  :: filetable  !< output table filename
        character(len=STDLEN), intent(in) :: filenames(:)!< array of filenames
        integer :: nl, funit, iline
        nl = size(filenames)
        funit = get_fileunit()
        open(unit=funit, status='replace', file=filetable)
        do iline=1,nl
            write(funit,'(a)') trim(filenames(iline))
        end do
        close(funit)
    end subroutine write_filetable

    !> \brief  for converting a file generated by txtfile2arr back to an array
    function txtfile2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam    !< input table filename
        real, allocatable :: arr(:)             !< array of filenames
        integer :: i, n, alloc_stat, funit
        logical :: here
        inquire(FILE=fnam, EXIST=here)
        if( here )then
            funit = get_fileunit()
            n = nlines(fnam)
            allocate( arr(n), stat=alloc_stat )
            if( alloc_stat /= 0 ) then
                write(*,'(a)') 'ERROR: Allocation failure!'
                write(*,'(a)') 'In: file2arr; simple_fileio      '
                stop
            endif
            open(unit=funit, status='old', file=fnam)
            do i=1,n
                read(funit,*) arr(i)
            end do
            close(funit)
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
        integer :: n1, n2, alloc_stat, cnt, funit, i
        logical :: here1, here2
        inquire(FILE=file1, EXIST=here1)
        inquire(FILE=file2, EXIST=here2)
        n1 = 0
        n2 = 0
        if( here1 ) n1 = nlines(file1)
        if( here2 ) n2 = nlines(file2)
        allocate( arr(n1+n2), stat=alloc_stat )
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') 'In: merge_txtfiles; simple_fileio      '
            stop
        endif
        if( here1 )then
            funit = get_fileunit()
            open(unit=funit, status='old', file=file1)
            cnt = 0
            do i=1,n1
                cnt = cnt+1
                read(funit,*) arr(cnt)
            end do
            close(funit)
            if( .not. here2 ) return
        else
            funit = get_fileunit()
            open(unit=funit, status='old', file=file2)
            cnt = 0
            do i=1,n2
                cnt = cnt+1
                read(funit,*) arr(cnt)
            end do
            close(funit)
            return
        endif
        open(unit=funit, status='old', file=file2)
        do i=1,n2
            cnt = cnt+1
            read(funit,*) arr(cnt)
        end do
        close(funit)
    end function merge_txtfiles

    !> \brief  for converting a file generated by file2arr back to an array
    function file2iarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam             !< input table filename 
        integer, allocatable :: arr(:)                   !< array of filenames   
        integer :: recsz, i, n, alloc_stat, funit, ival
        logical :: here
        inquire(FILE=fnam, EXIST=here)
        if( here )then
            inquire(iolength=recsz) ival
            funit = get_fileunit()
            open(unit=funit, status='old', file=fnam, access='direct', form='unformatted', recl=recsz)
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
            close(funit)
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
        integer :: recsz, i, funit
        inquire(iolength=recsz) rval
        rval = size(arr)
        funit = get_fileunit()
        open(unit=funit, status='replace', file=fnam, access='direct', form='unformatted', recl=recsz)
        write(funit, rec=1) rval
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        close(funit)
    end subroutine arr2file_1

    !> \brief  for converting a file generated by arr2file back to an array
    function file2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam  !< input table filename 
        real, allocatable            :: arr(:) !< array of filenames
        real    :: rval
        integer :: recsz, i, n, alloc_stat, funit
        logical :: here
        inquire(FILE=fnam, EXIST=here)
        if( here )then
            inquire(iolength=recsz) rval
            funit = get_fileunit()
            open(unit=funit, status='old', file=fnam, access='direct', form='unformatted', recl=recsz)
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
            close(funit)
        else
            write(*,*) fnam, ' does not exist; file2rarr; simple_fileio      '
            stop
        endif
    end function file2rarr

    !> \brief  for converting an integer array 2 file
    subroutine arr2file_2( arr, fnam )
        integer,          intent(in) :: arr(:)!< array of data
        character(len=*), intent(in) :: fnam !< output filename 
        integer :: recsz, i, funit, ival
        inquire(iolength=recsz) ival
        ival = size(arr)
        funit = get_fileunit()
        open(unit=funit, status='replace', file=fnam, access='direct', form='unformatted', recl=recsz)
        write(funit, rec=1) ival
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        close(funit)
    end subroutine arr2file_2

    !> \brief  for converting a real 2D array 2 file
    subroutine arr2D2file( arr, fnam )
        real,             intent(in) :: arr(:,:) !< array of data
        character(len=*), intent(in) :: fnam     !< output filename 
        real    :: dim1, dim2
        integer :: funit, io_stat
        dim1 = real(size(arr,dim=1))
        dim2 = real(size(arr,dim=2))
        funit = get_fileunit()
        open(unit=funit,access='STREAM',file=fnam,action='write',status='replace')
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
        close(funit)
    end subroutine arr2D2file

    !> \brief  for converting a real 2D array 2 file
    function file2arr2D( fnam ) result( arr )
        character(len=*), intent(in) :: fnam   !< input filename 
        real, allocatable :: arr(:,:)          !< array of data    
        real    :: dim1r, dim2r
        integer :: dim1, dim2, funit, io_stat, alloc_stat
        funit = get_fileunit()
        open(unit=funit,access='STREAM',file=fnam,action='read',status='old')
        read(unit=funit,pos=1,iostat=io_stat) dim1r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 1 from: ', trim(fnam)
            stop 'I/O error; file22Darr; simple_fileio      '
        endif
        read(unit=funit,pos=5,iostat=io_stat) dim2r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 5 from: ', trim(fnam)
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
            io_stat, ' when reading stream startbyte 9 from: ', trim(fnam)
            stop 'I/O error; file22Darr; simple_fileio      '
        endif
        close(funit)
    end function file2arr2D

    !> \brief  for converting a real array 2 file
    subroutine arr2txtfile( arr, fnam )
        real,             intent(in) :: arr(:) !< array of data  
        character(len=*), intent(in) :: fnam !< output filename 
        integer :: i, funit
        funit = get_fileunit()
        open(unit=funit, status='replace', file=fnam)
        do i=1,size(arr)
            write(funit,*) arr(i)
        end do
        close(funit)
    end subroutine arr2txtfile

end module simple_fileio      
