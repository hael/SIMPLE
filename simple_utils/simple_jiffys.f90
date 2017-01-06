!==Class simple_jiffys
!
! simple_jiffys provides jiffys. The code is distributed with the hope that it will be useful, 
! but _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification is regulated by the GNU General 
! Public License. *Author:* Hans Elmlund, 2011-08-18.
! 
!==Changes are documented below
!
module simple_jiffys
use simple_defs ! singleton
use, intrinsic :: iso_c_binding
implicit none

logical, parameter :: debug=.false.

character(len=*), parameter :: LOWER_CASE_LETTERS = 'abcdefghijklmnopqrstuvwxyz'
character(len=*), parameter :: UPPER_CASE_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(len=*), parameter :: INTEGERS = '1234567890'
character(kind=c_char,len=*), parameter :: NEW_LINES_C = C_FORM_FEED // &
C_NEW_LINE // C_CARRIAGE_RETURN // C_VERTICAL_TAB
character(kind=c_char,len=*), parameter :: BLANK_C_CHARACTERS = C_NULL_CHAR // C_HORIZONTAL_TAB
character(len=*), parameter :: BLANK_CHARACTERS =  ' '//BLANK_C_CHARACTERS
character(len=*), parameter :: BLANKS_AND_NEW_LINES = BLANK_CHARACTERS // NEW_LINES_C
character(len=*), parameter :: INTEGER_DIGITS = '10'   ! Maximum number of digits expected when reading an integer

private :: arr2file_1, arr2file_2, get_mrcfile_info, get_spifile_info, debug, LOWER_CASE_LETTERS,&
UPPER_CASE_LETTERS, INTEGERS, NEW_LINES_C, BLANK_C_CHARACTERS, BLANK_CHARACTERS, BLANKS_AND_NEW_LINES,&
INTEGER_DIGITS, assert_eq_2, assert_eq_3, assert_eq_4, assert_eq_n
public

interface arr2file
    module procedure arr2file_1
    module procedure arr2file_2
end interface arr2file

interface assert_eq
    module procedure assert_eq_2,assert_eq_3,assert_eq_4,assert_eq_n
end interface assert_eq

interface swap
    module procedure swap_i,swap_r,swap_rv,swap_c, swap_cv,swap_cm,&
    masked_swap_rs,masked_swap_rv,masked_swap_rm
end interface swap

contains
    
    ! FILE FUNCTIONS
    
    !> \brief  is for counting the number of lines in a textfile
    function nlines( fname ) result( n )
        character(len=*), intent(in) :: fname
        integer                      :: n, funit, ios
        logical                      :: here
        character(len=1)             :: junk
        inquire(FILE=fname, EXIST=here)
        if( here )then
            funit = get_fileunit( )
            open( unit=funit, file=fname )
            n = 0
            do
                 read(funit,*,IOSTAT=ios) junk
                 if(ios /= 0)then
                     exit
                 else
                     n = n+1
                 endif
            end do
            close( unit=funit )
        else
            n = 0
        endif
    end function nlines
    
    !> \brief  is returning the size of a binary file
    function filelength( fname ) result( filesz )
        character(len=*), intent(in) :: fname
        integer                      :: filesz, funit, ios, cnt
        logical                      :: here
        character(len=1)             :: junk
        inquire(FILE=fname, EXIST=here)
        if( here )then
            funit = get_fileunit( )
            open( unit=funit, file=fname, action='read',&
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
    
    !> \brief  is returning the record size of a binary file
    function reclength( fname, nentries ) result( recsz )
        character(len=*), intent(in) :: fname
        integer, intent(in)          :: nentries
        integer                      :: recsz
        recsz = filelength(fname)/nentries
    end function reclength
    
    !> \brief  is for deleting binary file
    subroutine del_binfile( file )
        character(len=*), intent(in) :: file
        open(29, file=file, form='unformatted')
        close (29, status='delete')
    end subroutine del_binfile
    
    !> \brief  is for deleting text file
    subroutine del_txtfile( file )
        character(len=*), intent(in) :: file
        open (29, file=file, form='formatted')
        close (29, status='delete')
    end subroutine del_txtfile
    
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
    
    !> \brief  is for checking file open stat
    subroutine fopen_err( message, file_stat )
        character(len=*), intent(in) :: message
        integer, intent(in)          :: file_stat
        if( file_stat /= 0 ) then
            write(*,'(a)') 'ERROR: File opening failure!'
            write(*,'(a)') message
            stop
        endif
    end subroutine fopen_err
    
    !> \brief  is for raising command line exception
    subroutine cmdline_err( cmdstat, cmdlen, arg, pos )
        integer, intent(in)          :: cmdstat, cmdlen, pos
        character(len=*), intent(in) :: arg
        if( cmdstat == -1 )then
            write(*,*) 'ERROR! while parsing the command line; simple_exec'
            write(*,*) 'The string length of argument: ', arg, 'is: ', cmdlen
            write(*,*) 'which likely exeeds the lenght limit STDLEN'
            write(*,*) 'Create a symbolic link with shorter name in the cwd'
            stop
        endif
        if( arg(:pos-1) .ne. 'prg' )then
            write(*,'(a)') 'ERROR!'
            write(*,'(a)') 'prg=simple_program required to be first on command line'
            write(*,'(a)') 'Please, refer to the manual for a comprehensive '
            write(*,'(a)') 'list of all programs and their specific documentation'
            stop
        endif
    end subroutine cmdline_err
    
    !>  \brief  Check a file exists on disk
    logical function file_exists(fname)
        character(len=*), intent(in) :: fname
        inquire(file=trim(adjustl(fname)), exist=file_exists)
    end function file_exists
    
    !> \brief   Find file size in bytes
    function file_size(fname) result(sz)
        character(len=*), intent(in) :: fname
        integer(kind=8)              :: sz
        inquire(file=trim(adjustl(fname)),size=sz)
    end function file_size
    
    !>  \brief  Check whether a IO unit is current open
    logical function is_open( unit_number )
        integer, intent(in)   :: unit_number
        integer               :: io_status
        character(len=STDLEN) :: io_message
        io_status = 0
        inquire(unit=unit_number, opened=is_open,iostat=io_status,iomsg=io_message)
        if (io_status .ne. 0) then
            print *, 'is_open: IO error ', io_status, ': ', trim(adjustl(io_message))
            stop 'IO error; is_open; simple_smpd'
        endif
    end function is_open
    
    !> \brief  for converting a real array 2 file
    subroutine arr2file_1( arr, fnam )
        real, intent(in)             :: arr(:)
        character(len=*), intent(in) :: fnam
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
        character(len=*), intent(in) :: fnam
        real, allocatable            :: arr(:)
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
            call alloc_err("In: file2arr; simple_jiffys", alloc_stat)
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            close(funit)
        else
            write(*,*) fnam, ' does not exist; file2rarr; simple_jiffys'
            stop 
        endif
    end function file2rarr
    
    !> \brief  for converting an integer array 2 file
    subroutine arr2file_2( arr, fnam )
        integer, intent(in)          :: arr(:)
        character(len=*), intent(in) :: fnam
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
        real, intent(in)             :: arr(:,:)
        character(len=*), intent(in) :: fnam
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
            stop 'I/O error; arr2D2file; simple_jiffys'
        endif
        write(unit=funit,pos=5,iostat=io_stat) dim2
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when writing stream startbyte 5 to: ', trim(fnam)
            stop 'I/O error; arr2D2file; simple_jiffys'
        endif
        write(unit=funit,pos=9,iostat=io_stat) arr(:,:)
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when writing stream startbyte 9 to: ', trim(fnam)
            stop 'I/O error; arr2D2file; simple_jiffys'
        endif
        close(funit)
    end subroutine arr2D2file
    
    !> \brief  for converting a real 2D array 2 file
    function file2arr2D( fnam ) result( arr )
        character(len=*), intent(in) :: fnam
        real, allocatable            :: arr(:,:)
        real    :: dim1r, dim2r
        integer :: dim1, dim2, funit, io_stat, alloc_stat
        funit = get_fileunit()
        open(unit=funit,access='STREAM',file=fnam,action='read',status='old')
        read(unit=funit,pos=1,iostat=io_stat) dim1r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 1 from: ', trim(fnam)
            stop 'I/O error; file22Darr; simple_jiffys'
        endif
        read(unit=funit,pos=5,iostat=io_stat) dim2r
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 5 from: ', trim(fnam)
            stop 'I/O error; file22Darr; simple_jiffys'
        endif
        dim1 = nint(dim1r)
        dim2 = nint(dim2r)
        if( allocated(arr) ) deallocate(arr)
        allocate( arr(dim1,dim2), stat=alloc_stat )
        call alloc_err("In: simple_jiffys :: file22Darr", alloc_stat )
        read(unit=funit,pos=9,iostat=io_stat) arr(:,:)
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(arr2D2file): I/O error ',&
            io_stat, ' when reading stream startbyte 9 from: ', trim(fnam)
            stop 'I/O error; file22Darr; simple_jiffys'
        endif
        close(funit)
    end function file2arr2D
    
    !> \brief  for converting a real array 2 file
    subroutine arr2txtfile( arr, fnam )
        real, intent(in)             :: arr(:)
        character(len=*), intent(in) :: fnam
        integer :: i, funit
        funit = get_fileunit()
        open(unit=funit, status='replace', file=fnam)
        do i=1,size(arr)
            write(funit,*) arr(i)
        end do
        close(funit)
    end subroutine arr2txtfile

    !> \brief  for reading raw images using stream access
    subroutine read_raw_image( fname, mat, first_byte )
        character(len=*), intent(in)  :: fname
        double precision, intent(out) :: mat(:,:,:)
        integer, intent(in)           :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message
        filnum = get_fileunit()
        open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
        read(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the read was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        close(filnum)
    end subroutine read_raw_image
    
    !> \brief  for writing raw images using stream access
    subroutine write_raw_image( fname, mat, first_byte )
        character(len=*), intent(in) :: fname
        real, intent(out)            :: mat(:,:,:)
        integer, intent(in)          :: first_byte
        integer :: filnum, io_stat
        character(len=100) :: io_message
        filnum = get_fileunit()
        open(unit=filnum, status='REPLACE', action='WRITE', file=fname, access='STREAM')
        write(unit=filnum,pos=first_byte,iostat=io_stat,iomsg=io_message) mat
        ! Check the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(rwSlices): I/O error ', io_stat, ' when reading from: ', fname
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_raw_image; simple_jiffys'
        endif
        close(filnum)
    end subroutine write_raw_image
    
    !> \brief  for converting a file generated by txtfile2arr back to an array
    function txtfile2rarr( fnam ) result( arr )
        character(len=*), intent(in) :: fnam
        real, allocatable            :: arr(:)
        integer :: i, n, alloc_stat, funit
        logical :: here
        inquire(FILE=fnam, EXIST=here)
        if( here )then
            funit = get_fileunit()
            n = nlines(fnam)
            allocate( arr(n), stat=alloc_stat )
            call alloc_err("In: file2arr; simple_jiffys", alloc_stat)
            open(unit=funit, status='old', file=fnam)
            do i=1,n
                read(funit,*) arr(i)
            end do
            close(funit)
        else
            write(*,*) fnam
            stop 'file does not exist; txtfile2rarr; simple_jiffys'
        endif
    end function txtfile2rarr
    
    !> \brief  merging two text files into a single array
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
        call alloc_err('In: merge_txtfiles; simple_jiffys', alloc_stat)
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
        character(len=*), intent(in) :: fnam
        integer, allocatable         :: arr(:)
        integer :: recsz, i, n, alloc_stat, funit, ival
        logical :: here
        inquire(FILE=fnam, EXIST=here)
        if( here )then
            inquire(iolength=recsz) ival
            funit = get_fileunit()
            open(unit=funit, status='old', file=fnam, access='direct', form='unformatted', recl=recsz)
            read(funit, rec=1) n
            allocate( arr(n), stat=alloc_stat )
            call alloc_err("In: file2arr; simple_jiffys", alloc_stat)
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            close(funit)
        else
            write(*,*) fnam
            stop 'file does not exist; file2iarr; simple_jiffys'
        endif
    end function file2iarr
    
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
        character(len=*), intent(in) :: fname, suffix
        character(len=STDLEN)        :: fbody
        integer :: pos
        pos = index(fname, suffix) ! position of suffix
        fbody = fname(:pos-1)
    end function get_fbody
    
    !> \brief  is for putting a new extension on filename
    function fname_new_ext( fname, suffix ) result( new_fname )
        character(len=*), intent(in) :: fname, suffix
        character(len=STDLEN)        :: new_fname
        new_fname = trim(get_fbody(trim(fname), fname2ext(trim(fname))))//suffix
    end function fname_new_ext
    
    !>  \brief  Return the 3-letter extension of a fname if present (without the period)
    pure function fname2ext( fname )
        character(len=*), intent(in)  :: fname
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
        character(len=*), intent(in)  :: fname
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
    
    !>  \brief  returns the integer number identifier of a filename
    subroutine fname2ind( str, ivar )
        character(len=*), intent(in)  :: str
        integer, intent(out)          :: ivar
        logical, allocatable          :: pos(:)
        character(len=:), allocatable :: str_copy
        integer :: j, lstr, io_stat, nrrange(2)
        lstr = len(str)
        pos = map_str_nrs(str)
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
        deallocate(pos,str_copy)
    end subroutine fname2ind
    
    !>  \brief Return a one letter code for the file format designated by the extension in the fname
    !!         if .mrc: M
    !!         if .spi: S
    !!         if .img: I
    !!         if .hed: I
    !!         else: N
    pure function fname2format( fname )
        character(len=*), intent(in)  :: fname
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
    pure logical function same_format( fname1, fname2 )
        character(len=*), intent(in) :: fname1, fname2
        character(len=1) :: form1, form2
        form1 = fname2format(fname1)
        form2 = fname2format(fname2)
        same_format = form1 == form2
    end function same_format

    !>  \brief  is for finding logical dimension and numbe rof particles in stack
    subroutine find_ldim_nptcls( fname, ldim, nptcls, doprint, formatchar, endconv )
        character(len=*), intent(in)                         :: fname      !< filename
        integer, intent(out)                                 :: ldim(3)    !< logical dimension
        integer, intent(out)                                 :: nptcls     !< number of particles
        logical, intent(in), optional                        :: doprint    !< do print or not
        character(len=1), intent(in), optional               :: formatchar !< input format
        character(len=:), intent(out), allocatable, optional :: endconv    !< endian conversion
        integer                       :: mode, iform, maxim
        real                          :: smpd
        character(len=:), allocatable :: conv
        character(len=1)              :: form
        logical                       :: ddoprint
        ddoprint = .false. 
        if( present(doprint) ) ddoprint = doprint
        if( present(formatchar) )then
            form = formatchar
        else
            form = fname2format(fname)
        endif
        nptcls = 0
        select case (form)
            case('M','F')
                call get_mrcfile_info(fname, ldim, form, smpd, ddoprint )
                nptcls = ldim(3)
            case('S')
                call get_spifile_info(fname, ldim, iform, maxim, smpd, conv, ddoprint)
                nptcls = maxim
            case DEFAULT
                write(*,*) 'fname: ', fname
                write(*,*) 'format descriptor: ', fname2format(fname)
                stop 'File format not supported; find_ldim_nptcls; simple_procimgfile'
        end select
        if( present(endconv) )then
            if( allocated(endconv) ) deallocate(endconv)
            select case (form)
                case('M','F')
                    allocate(endconv, source='NATIVE')
                case('S')
                    allocate(endconv, source=conv)
            end select
        endif
    end subroutine find_ldim_nptcls
    
    !>  \brief is for gettign a part of the info in a MRC image header
    subroutine get_mrcfile_info( fname, ldim, form, smpd, doprint )
        use simple_imghead, only: ImgHead, MrcImgHead, MrcFeiImgHead
        character(len=*), intent(in) :: fname
        character(len=1), intent(in) :: form
        integer, intent(out)         :: ldim(3)
        real, intent(out)            :: smpd
        logical, intent(in)          :: doprint
        class(imghead), allocatable  :: hed
        integer :: filnum
        ldim = 0
        smpd = 0.
        if( file_exists(fname) )then
            select case(form)
                case('M')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    filnum = get_fileunit()
                    open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
                    call hed%read(filnum)
                    close(filnum)
                    ldim = hed%getDims()
                    smpd = hed%getPixSz()
                    if( doprint )then
                        write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', ldim(1), ldim(2), ldim(3)
                        write(*,'(a,1x,f15.8)')  'Pixel size: ', smpd
                    endif
                case('F')
                    allocate(MrcImgHead :: hed)
                    call hed%new
                    filnum = get_fileunit()
                    open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='BIG_ENDIAN')
                    call hed%read(filnum)
                    close(filnum)
                    if( doprint ) call hed%print
                case DEFAULT
                    write(*,*) 'The inputted file is not an MRC file; get_mrcfile_info; simple_jiffys'
                    write(*,*) fname
                    stop
            end select
        else
            write(*,*) 'The below file does not exists; get_mrcfile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif
    end subroutine get_mrcfile_info
    
    !>  \brief is for gettign a part of the info in a SPIDER image header
    subroutine get_spifile_info( fname, ldim, iform, maxim, smpd, conv, doprint )
        character(len=*), intent(in)               :: fname
        integer, intent(out)                       :: ldim(3), maxim, iform
        real, intent(out)                          :: smpd
        character(len=:), allocatable, intent(out) :: conv
        logical, intent(in)                        :: doprint
        real    :: spihed(40)
        integer :: filnum, cnt, i
        if( file_exists(fname) )then
            if( fname2format(fname) .eq. 'S' )then
                if( allocated(conv) ) deallocate(conv)
                filnum = get_fileunit()
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='NATIVE')
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='NATIVE')
                    call print_spihed
                    return
                endif
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='BIG_ENDIAN') ! 
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='BIG_ENDIAN')
                    call print_spihed
                    return
                endif
                open(unit=filnum, status='OLD', action='READ', file=fname, access='STREAM', convert='LITTLE_ENDIAN') ! 
                call read_spihed
                close(filnum)
                if( .not. any(ldim < 1) )then
                    allocate(conv, source='LITTLE_ENDIAN')
                    call print_spihed
                    return
                endif
            else
                write(*,*) 'The inputted file is not a SPIDER file; get_spifile_info; simple_jiffys'
                write(*,*) fname
                stop
            endif
        else
            write(*,*) 'The below file does not exists; get_spifile_info; simple_jiffys'
            write(*,*) fname
            stop
        endif
        
        contains
            
            subroutine read_spihed
                cnt = 0
                do i=1,40*4,4
                    cnt = cnt+1
                    read(unit=filnum ,pos=i) spihed(cnt)
                end do
                ldim  = int([spihed(12), spihed(2), spihed(1)])
                iform = int(spihed(5))
                maxim = int(spihed(26))
                smpd  = spihed(38)
            end subroutine
            
            subroutine print_spihed
                if( doprint )then
                    write(*,'(a,3(i0,1x))') 'Number of columns, rows, sections: ', int(spihed(12)), int(spihed(2)), int(spihed(1))
                    write(*,'(a,1x,i3)')    'Iform descriptor: ', int(spihed(5))
                    write(*,'(a,1x,f7.0)')  'The number of the highest image currently used in the stack: ', spihed(26)
                    write(*,'(a,1x,f7.3)')  'Pixel size: ', spihed(38)
                endif
            end subroutine
            
    end subroutine get_spifile_info
    
    !>  \brief  reads a filetable into an array
    subroutine read_filetable( filetable, filenames )
        character(len=*), intent(in)                    :: filetable
        character(len=STDLEN), allocatable, intent(out) :: filenames(:)
        integer :: nl, funit, alloc_stat, iline
        nl    = nlines(filetable)
        funit = get_fileunit()
        open(unit=funit, status='old', file=filetable)
        allocate( filenames(nl), stat=alloc_stat )
        call alloc_err('In: read_filetable; simple_jiffys', alloc_stat)
        do iline=1,nl
            read(funit,'(a256)') filenames(iline)
        end do
        close(funit)
    end subroutine read_filetable
    
    ! ERROR CHECKS

    !> \brief  is for checking allocation 
    subroutine alloc_err( message, alloc_stat )
        character(len=*), intent(in) :: message
        integer, intent(in)          :: alloc_stat
        if( alloc_stat /= 0 ) then
            write(*,'(a)') 'ERROR: Allocation failure!'
            write(*,'(a)') message
            stop
        endif
    end subroutine alloc_err
    
    ! STRING MANIPULATIONS
    
    !>  \brief  Convert string to all upper case
    !!  Adapted from
    !!  Figure 3.5B, pg 80, "Upgrading to Fortran 90", by Cooper Redwine,
    !!  1995 Springer-Verlag, New York.
    function upperCase( str )
        character(len=*),   intent(in)  ::  str
        character(len=len(str))  ::  upperCase
        integer :: i, n
        upperCase = str
        do i=1,len(str)
            n = index(LOWER_CASE_LETTERS,upperCase(i:i))
            if( n .ne. 0 ) upperCase(i:i) = UPPER_CASE_LETTERS(n:n)
        enddo
    end function upperCase
    
    !>  \brief  turn integer variable into character variable
    function int2str(intg) result(string)
        integer, intent(in)           :: intg
        character(len=:), allocatable :: string
        allocate(character(len=ndigs_int(intg)) :: string)
        write(string,'(i0)') intg
    end function int2str
    
    !>  \brief  turns a string into an integer variable
    subroutine str2int( string, io_stat, ivar )
        character(len=*), intent(in) :: string
        integer, intent(out)         :: io_stat, ivar
        read(string,*, iostat=io_stat) ivar
    end subroutine str2int
    
    !>  \brief  turn integer variable into zero padded character variable
    function int2str_pad(intg, numlen) result(string)
        integer, intent(in)           :: intg, numlen
        character(len=:), allocatable :: string, str_tmp, str_tmp2
        integer :: str_tmp_len
        str_tmp = int2str(intg)
        if( len(str_tmp) > numlen )then
            stop 'length of number > than desired length; simple_jiffys::in2str_zero_padded'
        else if( len(str_tmp) == numlen )then
            allocate(string, source=str_tmp)
        else
            do while( len(str_tmp) < numlen )
                allocate(str_tmp2, source='0'//str_tmp)
                if( allocated(string) ) deallocate(string)
                allocate(string, source=str_tmp2)
                deallocate(str_tmp)
                allocate(str_tmp, source=str_tmp2)
                deallocate(str_tmp2)
            end do
        endif
        deallocate(str_tmp)
    end function int2str_pad
    
    !>  \brief  number of digits in an integer
    pure integer function ndigs_int(intg)
        integer, intent(in) :: intg
        if (intg .eq. 0) then
            ndigs_int = 1
        else
            ndigs_int = int(log10(real(intg)))+1
        endif
    end function ndigs_int
    
    !> \brief copies the trimmed and justified version of the instr into outstr
    subroutine cpStr(instr, outstr)
        character(len=*), intent(in) :: instr
        character(len=:), allocatable, intent(inout) :: outstr
        outstr = trim(adjustl(instr))
    end subroutine cpStr
    
    !>  \brief  maps out the number characters in a string
    !!          .true. means is number, .false. means is not
    function map_str_nrs( str ) result( map )
        character(len=*), intent(in)  :: str
        logical, allocatable :: map(:)
        integer              :: lstr, i, j, ind
        lstr = len(str)
        allocate(map(lstr))
        map = .false.
        do j=1,lstr
            do i=0,9
                ind = index(str(j:j), int2str(i))
                if( ind == 0 )then
                    cycle
                else
                    map(j) = .true.
                    exit
                endif
            end do
         end do
    end function map_str_nrs
    
    !>  \brief works out whether a character string is a comment line
    logical function strIsComment( line )
        character(len=*), intent(in)  ::  line
        integer ::  pos1
        pos1 = firstNonBlank(line)
        ! if we didn't find a non-blank character, then the line must be blank!
        if( pos1 .eq. 0 )then
            strIsComment = .false.
        else if( scan(line(pos1:),'#!c;', back=.false.) .eq. 1 )then
            ! the first non-blank character is a comment character. if that's a c, we need to check that the following character is a space.
            if( stringsAreEqual(line(pos1:pos1),'c') )then
                if( stringsAreEqual(line(pos1+1:pos1+1),' ') )then
                    strIsComment = .true.
                else
                    strIsComment = .false.
                endif
            else
                strIsComment = .true.
            endif
        else
            strIsComment = .false.
        endif
    end function strIsComment

    !>  \brief  check for presence of substring
    logical function str_has_substr( str, substr )
        character(len=*), intent(in) :: str, substr
        str_has_substr = .not. (index(str, substr) == 0)
    end function str_has_substr

    !>  \brief works out whether a character string is blank
    pure logical function strIsBlank(line)
        character(len=*),   intent(in)  ::  line
        if( trim(line) .eq. '' )then
             strIsBlank = .true.
        else if( firstNonBlank(line) .eq. 0 )then
             strIsBlank = .true.
        else
             strIsBlank = .false.
        endif
    end function strIsBlank

    !>  \brief  Find the first non-blank character in a string and return its position.
    pure integer function firstNonBlank( string, back )
        character(len=*), intent(in)  :: string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last non-blank character in the string
        logical ::  bback
        ! reverse?
        bback = .false.
        if (present(back)) bback = back
        firstNonBlank = verify(string,blanks_and_new_lines,back=bback)
    end function firstNonBlank

    !>  \brief  find the first blank character in a string and return its position.
    pure integer function firstBlank(string, back)
        character(len=*), intent(in)  :: string
        logical, optional, intent(in) :: back !<  if .true., we'll look for the last blank character in the string
        logical ::  bback
        ! reverse?
        if (present(back)) then
            bback = back
        else
            bback = .false.
        endif
        firstBlank = scan(string,blanks_and_new_lines,back=bback)
    end function firstBlank

    !>  \brief  test whether two strings are equal, ignoring blank characters
    logical function stringsAreEqual(instr1, instr2, case_sensitive)
        character(len=*), intent(in)  :: instr1
        character(len=*), intent(in)  :: instr2
        logical, optional, intent(in) :: case_sensitive
        integer :: first_non_blank_1, first_non_blank_2
        integer :: last_non_blank_1, last_non_blank_2
        logical :: ccase_sensitive
        character(len=:), allocatable ::  local_str1, local_str2
        ! Will we be case sensitive?
        ccase_sensitive = .true.
        if( present(case_sensitive) ) ccase_sensitive = case_sensitive
        ! Copy to local string, with or without capitalization
        if( ccase_sensitive )then
            call cpStr(instr1,local_str1)
            call cpStr(instr2,local_str2)
        else
            local_str1 = upperCase(instr1)
            local_str2 = upperCase(instr2)
        endif
        ! Find positions of first and last non-blank characters
        first_non_blank_1 = verify(local_str1,blank_characters)
        last_non_blank_1  = verify(local_str1,blank_characters,back=.true.)
        first_non_blank_2 = verify(local_str2,blank_characters)
        last_non_blank_2  = verify(local_str2,blank_characters,back=.true.)
        if( first_non_blank_1 .eq. 0 .and. first_non_blank_2 .eq. 0 )then
            ! both strings are blank
            stringsAreEqual = .true.
        else if( first_non_blank_1 .eq. 0 .or. first_non_blank_2 .eq. 0 )then
            ! one string is blank, the other isn't
            stringsAreEqual = .false.
        else if( index(local_str1(first_non_blank_1:last_non_blank_1), &
                       local_str2(first_non_blank_1:last_non_blank_2)) .eq. 1  &
                 .and. last_non_blank_1-first_non_blank_1 .eq. last_non_blank_2-first_non_blank_2 )then
            ! neither string is blank, and the strings match
            stringsAreEqual = .true.
        else
            ! all other cases
            stringsAreEqual = .false.
        endif
    end function stringsAreEqual
    
    !>  \brief  count the number of records on a line
    function cntRecsPerLine( line, separators ) result( nrecs )
        character(len=*)            :: line        !<  line to be split
        character(len=*), optional  :: separators  !<  characters which separate words, if not present, default is blank characters (space, tabs...)
        character(len=line_max_len) :: buffer
        integer :: nrecs, pos1, pos2
        if (strIsBlank(line) .or. strIsComment(line)) then
            nrecs = 0
        else
            buffer = trim(line)
            nrecs = 0
            do
                if (present(separators)) then
                    ! find the first non-separator character
                    pos1 = verify(buffer,separators,back=.false.)
                    if( pos1 .ne. 0 )then
                        ! find the last non-separator character after that
                        pos2 = pos1+scan(buffer(pos1:),separators,back=.false.)-2
                    endif
                else
                    ! find the first non-blank character
                    pos1 = firstNonBlank(buffer)
                    if (pos1 .ne. 0) then
                        ! find the last non-blank character after that
                        pos2 = pos1+firstBlank(buffer(pos1:))-2
                    endif
                endif
                if( pos1 .ne. 0 .and. pos2 .ne. 0 )then ! if we found a word
                    nrecs = nrecs + 1
                endif
                if( pos2+1 .ge. len_trim(buffer) )then ! if we reached end of the line
                    exit
                else
                    buffer = buffer(pos2+1:)
                endif
            enddo
        endif
    end function cntRecsPerLine

    ! PRETTY PROGRESS & ENDING
    
    !> \brief is for printing a progress bar
    subroutine progress(i, imax)
        integer, intent(in) :: i, imax
        integer :: k, curr, prev
        character(len=1) :: back
        back = char(8)
        if( i>= imax )then
            ! delete the bar and the percentage
            write(6,'(256a1)', advance='no') (back, k =1,59)
            ! write new bar and percentage
            write(6,'(2x,1a4,2x,1a1,256a1)', advance='no') '100%','|', ('=', k=1,50)
            write(6,'(a)') '| done.'
            flush 6            
        else
            prev = nint( 100.*real(i-1)/real(imax) )
            curr = nint( 100.*real(i)/real(imax) )
            if( curr>prev )then
                ! delete the bar and the percentage
                write(6,'(256a1)', advance='no') (back, k =1,(50*i/imax)+9)
                ! write new bar and percentage
                write(6,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') curr,'%','|', ('=', k=1,50*i/imax)
                flush 6
            endif
        endif
    end subroutine progress

    !> \brief  is for pretty ending
    subroutine simple_end( str )
        character(len=*), intent(in) :: str
        write(*,'(A)') "       _______ _____ _______  _____         _______"
        write(*,'(A)') "       |______   |   |  |  | |_____] |      |______"
        write(*,'(A)') "       ______| __|__ |  |  | |       |_____ |______"      
        write(*,'(A)') " "                                                            
        write(*,'(A)') " _)_ ( _   _     ) o  _             _   _   _   o  _   _"  
        write(*,'(A)') " (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"   
        write(*,'(A)') "         (_                  (\   (_         _)      (_"
        write(*,'(A)') ""
        write(*,'(A)') str
    end subroutine simple_end
  
    !> \brief  for pretty haloween ending
    subroutine haloween_end( str )
        character(len=*), intent(in) :: str
        write(*,'(A)') " #"     
        write(*,'(A)') " ##"
        write(*,'(A)') " ###"   
        write(*,'(A)') "  ####"   
        write(*,'(A)') "   #####              _______ _____ _______  _____         _______"
        write(*,'(A)') "   #######            |______   |   |  |  | |_____] |      |______"
        write(*,'(A)') "    #######           ______| __|__ |  |  | |       |_____ |______"
        write(*,'(A)') "    ########"
        write(*,'(A)') "    #########    _)_ ( _   _     ) o  _             _   _   _   o  _   _"
        write(*,'(A)') "    ##########   (_   ) ) )_)   (  ( ) ) (_( \)    )_) ) ) (_(  ( ) ) )_)"
        write(*,'(A)') "    ##########            (_                  (\   (_         _)      (_"
        write(*,'(A)') "   ###########"
        write(*,'(A)') " ##############                                                    #"          
        write(*,'(A)') "###############                                                    #"
        write(*,'(A)') " ##############                                                   ##" 
        write(*,'(A)') "   #############                                                 ##"
        write(*,'(A)') "    #############                                              ###"
        write(*,'(A)') "    ###############                                         #####"
        write(*,'(A)') "     ###############                                    ########"
        write(*,'(A)') "     ################                                ###########"
        write(*,'(A)') "     ################                              ############"
        write(*,'(A)') "     ################                           #############"
        write(*,'(A)') "    #################       #                 ###############"
        write(*,'(A)') "   ###################     ##    #           ################"
        write(*,'(A)') "  #####################   ###   ##          #################"
        write(*,'(A)') " #######################  ########         ###################"
        write(*,'(A)') "   ######################  #######        #####################"
        write(*,'(A)') "     #############################       #######################"
        write(*,'(A)') "        #########################       #####################"
        write(*,'(A)') "           ################################################"
        write(*,'(A)') "            ##############################################"
        write(*,'(A)') "              ###########################################"
        write(*,'(A)') "               #########################################"
        write(*,'(A)') "                #######################################"
        write(*,'(A)') "                 ######################################"
        write(*,'(A)') "                 ######################################"
        write(*,'(A)') "                 ###########################      ######"
        write(*,'(A)') "                 ####  ###################           ###"
        write(*,'(A)') "                 ###    ################              ##"
        write(*,'(A)') "                 ##     ##  ##########                 #"
        write(*,'(A)') "                 #      #   #### ####"
        write(*,'(A)') "                            ##    ###"
        write(*,'(A)') "                            #     ##"
        write(*,'(A)') "                                  #"
        write(*,'(A)') str
    end subroutine haloween_end
    
    !> \brief  Lexographical sort. On input, strArr is a one-dimensional array of character strings to be 
    !!         sorted in ascending lexical order. On output, strArr is the sorted array. The characters of 
    !!         the elements of the string array are not modified. If blanks or punctuation characters are 
    !!         to be ignored, this needs to be taken care of before calling.
    subroutine lexSort( strArr, CaseSens ) 
        character(len=*),  intent(inout) :: strArr(:)
        logical, optional, intent(in)    :: CaseSens
        integer, allocatable :: indexarray(:)
        integer              :: low, high, alloc_stat, k
        logical              :: LCaseSens
        LCaseSens = .false.
        if( present(CaseSens) ) LCaseSens = CaseSens
        low  = 1 
        high = size(strArr) 
        allocate(indexarray(high), stat=alloc_stat)
        call alloc_err('In: lexSort; simple_jiffys', alloc_stat) 
        indexarray = (/(k,k=low,high)/) 
        call quicksort(low, high) 
        strArr = strArr(indexarray) 
        deallocate(indexarray)
    
      contains
      
          recursive subroutine quicksort( low, high )
              integer, intent (in) :: low, high 
              integer :: pivotlocation 
              if( low < high )then 
                  call partition(low, high, pivotlocation) 
                  call quicksort(low, pivotlocation-1) 
                  call quicksort(pivotlocation+1, high) 
              end if 
          end subroutine quicksort
      
          subroutine partition( low, high, pivotlocation )
              integer, intent (in) :: low, high 
              integer, intent(out) :: pivotlocation 
              integer :: k, lastsmall 
              call swap(indexarray(low), indexarray((low+high)/2)) 
              lastsmall = low 
              do k=low+1,high 
                  if(strComp(strArr(indexarray(k)), strArr(indexarray(low))))then 
                      lastsmall = lastsmall + 1 
                      call swap(indexarray(lastsmall), indexarray(k)) 
                  end if 
              end do 
              call swap(indexarray(low), indexarray(lastsmall)) 
              pivotlocation = lastsmall 
          end subroutine partition
      
          subroutine swap( m, n ) 
              integer, intent (inout) :: m, n 
              integer :: temp 
              temp = m 
              m = n 
              n = temp 
          end subroutine swap
      
          function strComp( p, q ) result( lexLess ) 
              character(len=*), intent(in) :: p, q 
              logical :: lexLess 
              integer :: kq, k 
              if( LCaseSens )then 
                  lexLess = p < q 
              else 
                  kq = 1 
                  do k=1,max(len_trim(p),len_trim(q)) 
                      if( UpperCase(p(k:k)) == UpperCase(q(k:k)) )then 
                          cycle 
                      else 
                          kq = k 
                          exit 
                      end if 
                  end do 
                  lexLess = UpperCase(p(kq:kq)) < UpperCase(q(kq:kq)) 
              end if 
          end function strComp
      
          function UpperCase( letter ) result( L ) 
              character(len=*), intent (in) :: letter 
              character(len=1)              :: L 
              character(len=26), parameter  :: Lower = "abcdefghijklmnopqrstuvwxyz"
              character(len=26), parameter  :: Upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ" 
              integer :: k 
              k = index(Lower, letter) 
              if (k > 0) then 
                  L = Upper(k:k) 
              else 
                  L = letter 
              end if 
          end function UpperCase
        
    end subroutine lexSort
    
    ! ASSERTIONS AND SWAPS FROM NR
    
    function assert_eq_2(n1,n2,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2
        integer :: assert_eq_2
        if (n1 == n2) then
            assert_eq_2=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_2; simple_jiffys'
        end if
    end function assert_eq_2

    function assert_eq_3(n1,n2,n3,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3
        integer :: assert_eq_3
        if (n1 == n2 .and. n2 == n3) then
            assert_eq_3=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_3; simple_jiffys'
        end if
    end function assert_eq_3
          
    function assert_eq_4(n1,n2,n3,n4,string)
        character(len=*), intent(in) :: string
        integer, intent(in) :: n1,n2,n3,n4
        integer :: assert_eq_4
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
            assert_eq_4=n1
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_4; simple_jiffys'
        end if
    end function assert_eq_4
    
    function assert_eq_n(nn,string)
        character(len=*), intent(in) :: string
        integer, dimension(:), intent(in) :: nn
        integer :: assert_eq_n
        if (all(nn(2:) == nn(1))) then
            assert_eq_n=nn(1)
        else
            write (*,*) 'nrerror: an assert_eq failed with this tag: ', string
            stop 'program terminated by assert_eq_n; simple_jiffys'
        end if
    end function assert_eq_n
    
    subroutine swap_i(a,b)
          integer, intent(inout) :: a,b
          integer :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_i

    subroutine swap_r(a,b)
          real, intent(inout) :: a,b
          real :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_r

    subroutine swap_rv(a,b)
          real, dimension(:), intent(inout) :: a,b
          real, dimension(size(a)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_rv

    subroutine swap_c(a,b)
          complex, intent(inout) :: a,b
          complex :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_c

    subroutine swap_cv(a,b)
          complex, dimension(:), intent(inout) :: a,b
          complex, dimension(size(a)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_cv

    subroutine swap_cm(a,b)
          complex, dimension(:,:), intent(inout)  :: a,b
          complex, dimension(size(a,1),size(a,2)) :: dum
          dum=a
          a=b
          b=dum
    end subroutine swap_cm

    subroutine masked_swap_rs(a,b,mask)
          real,         intent(inout) :: a,b
          logical(lgt), intent(in)    :: mask
          real :: swp
          if (mask) then
              swp=a
              a=b
              b=swp
          end if
    end subroutine masked_swap_rs

    subroutine masked_swap_rv(a,b,mask)
          real,         dimension(:), intent(inout) :: a,b
          logical(lgt), dimension(:), intent(in)    :: mask
          real, dimension(size(a)) :: swp
          where (mask)
              swp=a
              a=b
              b=swp
          end where
    end subroutine masked_swap_rv

    subroutine masked_swap_rm(a,b,mask)
          real,         dimension(:,:), intent(inout) :: a,b
          logical(lgt), dimension(:,:), intent(in)    :: mask
          real, dimension(size(a,1),size(a,2)) :: swp
          where (mask)
              swp=a
              a=b
              b=swp
          end where
    end subroutine masked_swap_rm

end module simple_jiffys
