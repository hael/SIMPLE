module simple_fileio
use simple_defs
use simple_string
use simple_string_utils
use simple_error
use simple_syslib
implicit none
#include "simple_local_flags.inc"

interface arr2file
    module procedure arr2file_sp
    module procedure arr2file_dp
end interface arr2file

interface arr2txtfile
    module procedure arr2txtfile_1
    module procedure arr2txtfile_2
end interface arr2txtfile

interface del_files
    module procedure del_files_1
    module procedure del_files_2
end interface del_files

interface add2fbody
    module procedure add2fbody_1
    module procedure add2fbody_2
    module procedure add2fbody_3
end interface add2fbody

interface fname_new_ext
    module procedure fname_new_ext_1
    module procedure fname_new_ext_2
end interface fname_new_ext

interface filepath
    module procedure filepath_1
    module procedure filepath_2
    module procedure filepath_3
    module procedure filepath_4
end interface filepath

interface append2basename
    module procedure append2basename_1
    module procedure append2basename_2
end interface append2basename

interface swap_suffix
    module procedure swap_suffix_1
    module procedure swap_suffix_2
end interface swap_suffix

interface get_fbody
    module procedure get_fbody_1
    module procedure get_fbody_2
end interface get_fbody


contains

    subroutine fileiochk( message, iostat , die )
        character(len=*),  intent(in) :: message  !< error message
        integer,           intent(in) :: iostat   !< error status
        logical, optional, intent(in) :: die      !< do you want to terminate or not
        logical :: die_this
        if( iostat == 0 ) return
        die_this=.true.
        if( present(die) ) die_this=die
        write(logfhandle,'(a)') message
        if( iostat == -1 )then
            write(logfhandle,'(a)') "fileio: EOF reached "
        else if( iostat == -2 )then
            write(logfhandle,'(a)') "fileio: End-of-record reached "
        else
            if( die_this ) THROW_HARD('I/O')
        endif
    end subroutine fileiochk

    subroutine fopen(funit, file, status, action, iostat, access, form, recl, async, pad,&
        &decimal, round, delim, blank, convert, iomsg, position, errmsg)
        integer,                    intent(inout) :: funit
        class(string),              intent(in)    :: file
        integer,          optional, intent(inout) :: iostat
        integer,          optional, intent(inout) :: recl
        character(len=*), optional, intent(in)    :: status, access, async, action, &
            &blank, pad, form, decimal, round, delim, convert, iomsg, position, errmsg
        integer               :: iostat_this,  recl_this
        character(len=STDLEN) :: iomsg_this
        character(len=30)     :: async_this, access_this, action_this, status_this,&
            &blank_this, pad_this, decimal_this, delim_this, form_this, round_this,&
            &position_this, errmsg_this
        ! check to see if file is empty
        if( file%is_blank() )then
            write(logfhandle,*) 'simple_system::fopen filename blank'
            if(present(iomsg))  write(logfhandle,*) trim(adjustl(iomsg))
            if(present(errmsg)) write(logfhandle,*) "Message: ", trim(adjustl(errmsg))
            return
        endif
        errmsg_this="In simple_fileio::fopen "
        if (present(errmsg)) errmsg_this=errmsg
        if (.not. (present(iostat) .or. present(form) .or. present(recl) .or.&
            & present(async) .or. present(pad) .or. present(action) .or. present(status)&
            & .or. present(position) .or. present(access) .or. present(decimal) .or. &
            present(round) .or. present(delim) .or. present(blank) ) )then
            open(NEWUNIT=funit, FILE=file%to_char(),IOSTAT=iostat_this)
            call fileiochk(trim(adjustl(errmsg_this))//" fopen basic open "//file%to_char(), iostat_this,.false.)
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
            if (strings_are_equal(status, 'OLD',    .false.)) write(status_this, '(A)') upperCase(status)
            if (strings_are_equal(status, 'SCRATCH',.false.)) write(status_this, '(A)') upperCase(status)
            if (strings_are_equal(status, 'REPLACE',.false.)) write(status_this ,'(A)') upperCase(status)
            if (strings_are_equal(status, 'NEW',    .false.)) write(status_this, '(A)') upperCase(status)
        end if
        ! ACTION: READ, WRITE, or READWRITE (default).
        if (present(action))then
            if (strings_are_equal(action, 'WRITE',  .false.)) write(action_this ,'(A)') upperCase(action)
            if (strings_are_equal(action, 'READ',   .false.)) write(action_this ,'(A)') upperCase(action)
        end if
        if ( (strings_are_equal(status_this, 'NEW',.false.))  .and. &
            (strings_are_equal(action_this, 'READ',.false.))  .and. &
            (.not. file_exists(file) ) )then
            write(logfhandle,*) "::fopen incompatible status=NEW and action=READ ", file%to_char()," does not exist"
            return ! false
        end if
        if(present(position)) then
            if ( (strings_are_equal(status_this, 'OLD',.false.))  .and. &
                (strings_are_equal(position, 'APPEND',.false.))  .and. &
                (.not. file_exists(file) ) )then
                write(logfhandle,*) "::fopen incompatible status=OLD and position=APPEND  when ",&
                    file%to_char()," does not exist"
                write( status_this,'(A)')  upperCase('NEW')
            end if
        end if
        ! access: DIRECT (random access) or SEQUENTIAL  or STREAM (F2003)
        write(access_this ,'(A)') 'SEQUENTIAL'
        if (present(access))then
            if (strings_are_equal(access, 'DIRECT',.false.))  write(access_this ,'(A)') upperCase(access)
            if (strings_are_equal(access, 'STREAM',.false.))then
                write(access_this ,'(A)') upperCase(access)
            endif
        end if
        recl_this=-1
        if(present(recl)) recl_this=recl
        !! Common file open
        if (.not.( present(form) .or. present(async) .or. present(pad) .or. &
            &present(decimal) .or. present(round) .or. present(delim) .or. present(blank) ) )then
            if (strings_are_equal(access_this, 'DIRECT',.false.) .and. (recl_this > 0) ) then

                if (present(action))then
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,RECL=recl_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,RECL=recl_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this)
                    end if
                else ! no action
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,RECL=recl_this,&
                            &STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this,RECL=recl_this)
                    end if
                end if
            else
                if (present(action))then
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                            &ACTION=action_this,STATUS=status_this,ACCESS=access_this)
                    end if
                else ! no action
                    if (present(position))then
                        !! Appending to file
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this,POSITION=position_this)
                    else
                        open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                            &STATUS=status_this,ACCESS=access_this)
                    end if
                end if
            end if
            if(present(iostat))iostat=iostat_this
            if(funit/=0 .and. is_io(funit)) THROW_HARD("newunit returned "//int2str(funit))
            return
        end if
        recl_this=-1
        if(present(recl)) recl_this=recl
        write(pad_this,'(A)') 'YES'
        if (present(pad))then
            if (strings_are_equal(pad, 'NO',.false.))  write(pad_this ,'(A)') upperCase(pad)
        end if
        write(async_this,'(A)')'NO'
        if (present(async))then
            if (strings_are_equal(async, 'YES',.false.))  write(async_this ,'(A)')upperCase(async)
        end if
        write(blank_this,'(A)')'NULL'
        if (present(blank))then
            if (strings_are_equal(blank, 'ZERO',.false.)) write( blank_this ,'(A)') upperCase(blank)
        end if
        write(decimal_this,'(A)')'POINT'
        if (present(decimal))then
            if (strings_are_equal(decimal, 'COMMA',.false.)) write( decimal_this ,'(A)') upperCase(decimal)
        end if
        write(delim_this,'(A)')'NONE'
        if (present(delim))then
            if (strings_are_equal(delim, 'APOSTROPHE',.false.)) write( delim_this ,'(A)') upperCase(delim)
            if (strings_are_equal(delim, 'QUOTE',.false.))  write(delim_this ,'(A)') upperCase(delim)
        end if
        write(form_this,'(A)')'FORMATTED'
        if (present(form))then
            if (strings_are_equal(form, 'UNFORMATTED',.false.))  write(form_this ,'(A)') upperCase(form)
            if (strings_are_equal(form, 'BINARY',.false.))  write(form_this ,'(A)') upperCase(form)
        end if
        write(round_this,'(A)')'PROCESSOR_DEFINED'
        if (present(round))then
            if (strings_are_equal(round, 'SUPPRESS',.false.)) write( round_this ,'(A)') upperCase(round)
            if (strings_are_equal(round, 'PLUS',.false.))  write(round_this ,'(A)') upperCase(round)
            if (strings_are_equal(round, 'UNDEFINED',.false.)) write( round_this ,'(A)') upperCase(round)
        end if
        if(present(iomsg)) iomsg_this=iomsg
        ! execute open under specific conditions
        if (strings_are_equal(form_this, 'FORMATTED',.false.)) then
            open( NEWUNIT=funit,FILE=file%to_char(),IOSTAT=iostat_this,&
                &ACTION=action_this,STATUS=status_this,ACCESS=access_this,&
                &BLANK=blank_this,FORM='FORMATTED', ROUND=round_this,&
                &IOMSG=iomsg_this)
        else
            if (strings_are_equal(access_this, 'DIRECT',.false.))then
                open( NEWUNIT=funit, FILE=file%to_char(), IOSTAT=iostat_this, &
                    &ACTION=action_this, STATUS=status_this,&
                    &ACCESS=access_this, FORM=form_this, RECL=recl_this,&
                    &IOMSG=iomsg_this)
            else
                if (recl_this == -1)then
                    open( NEWUNIT=funit, FILE=file%to_char(), IOSTAT=iostat_this, &
                        &ACTION=action_this, STATUS=status_this,&
                        &ACCESS=access_this, FORM=form_this,&
                        &IOMSG=iomsg_this)
                else
                    open( NEWUNIT=funit, FILE=file%to_char(), IOSTAT=iostat_this, &
                        &ACTION=action_this, STATUS=status_this,&
                        &ACCESS=access_this, RECL=recl_this, FORM=form_this,&
                        &IOMSG=iomsg_this)
                end if
            end if
        end if
        call fileiochk(trim(adjustl(errmsg_this))//" fopen opening "//file%to_char(), iostat_this, .false.)
        if(is_io(funit)) THROW_HARD("newunit returned "//int2str(funit))
        if(present(iostat))iostat=iostat_this
        if(present(recl))recl=recl_this
    end subroutine fopen

    subroutine fclose( funit )
        integer, intent(in) :: funit
        if( is_open(funit) ) close(funit)
    end subroutine fclose

    subroutine wait_for_closure( fname )
        class(string), intent(in) :: fname
        logical :: exists, closed
        integer :: wait_time
        wait_time = 0
        do
            if( wait_time == 60 )then
                call simple_exception('been waiting for a minute for file: '//fname%to_char(), 'simple_fileio.f90', __LINE__,l_stop=.false.)
                wait_time = 0
                flush(OUTPUT_UNIT)
            endif
            exists = file_exists(fname)
            closed = .false.
            if( exists )closed = .not. is_file_open(fname)
            if( exists .and. closed )exit
            call sleep(1)
            wait_time = wait_time + 1
        enddo
    end subroutine wait_for_closure

    function nlines( fname ) result( n )
        class(string), intent(in) :: fname !< input filename
        integer          :: n, funit, ios,io_status
        character(len=1) :: junk
        if( file_exists(fname) )then
            call fopen(funit, fname, status='unknown', action='read', iostat=io_status)
            call fileiochk(":nlines error opening file "//fname%to_char(), io_status)
            n = 0
            do
                read(funit,*,IOSTAT=ios) junk
                if(ios /= 0)then
                    exit
                else
                    if(.not. str_is_comment(junk))  n = n + 1
                endif
            end do
            call fclose(funit)
        else
            n = 0
        endif
    end function nlines

    function filelength( fname ) result( filesz )
        class(string), intent(in) :: fname !< input filename
        integer                   :: filesz, funit, ios, cnt,recl
        character(len=1)          :: junk
        if(  file_exists(fname) )then
            recl=1
            call fopen(funit, fname, ACTION='read', IOSTAT=ios,&
                ACCESS='direct', form='unformatted', recl=recl)
            call fileiochk('simple_fileio :: nlines opening '//fname%to_char(), ios)
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
            call fclose(funit)
        else
            filesz = 0
        endif
    end function filelength

    function funit_size(unit) result(sz)
        integer, intent(in) :: unit !< input file unit
        integer(kind=8)     :: sz
        inquire(unit=unit,size=sz)
    end function funit_size

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
        integer, parameter      :: MAX_UNIT_NUMBER = 1000
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

    function add2fbody_1( fname, suffix, str ) result( newname )
        class(*),      intent(in)  :: fname
        class(string), intent(in) ::suffix, str
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp = cast_str_types(fname) 
        pos = index(fname_tmp, suffix%to_char()) ! position of suffix
        newname = fname_tmp(:pos-1)//str%to_char()//suffix%to_char()
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function add2fbody_1

    function add2fbody_2( fname, suffix, str ) result( newname )
        class(*),         intent(in)  :: fname
        character(len=*), intent(in)  :: suffix, str
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp = cast_str_types(fname) 
        pos = index(fname_tmp, trim(suffix)) ! position of suffix
        newname = fname_tmp(:pos-1)//trim(str)//trim(suffix)
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function add2fbody_2

    function add2fbody_3( fname, suffix, str ) result( newname )
        class(*),         intent(in)  :: fname
        class(string),    intent(in)  :: suffix
        character(len=*), intent(in)  :: str
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp  = cast_str_types(fname)
        pos = index(fname_tmp, suffix%to_char()) ! position of suffix
        newname = fname_tmp(:pos-1)//trim(str)//suffix%to_char()
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function add2fbody_3

    function append2basename_1( fname, suffix ) result( newname )
        class(string), intent(in) :: fname, suffix 
        type(string) :: bname, newname, ext
        bname   = basename(fname)
        ext     = fname2ext(bname)
        newname = add2fbody(bname, string('.'//ext%to_char()), suffix)
        call bname%kill
        call ext%kill
    end function append2basename_1

    function append2basename_2( fname, suffix ) result( newname )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: suffix 
        type(string) :: bname, newname, ext
        bname   = basename(fname)
        ext     = fname2ext(bname)
        newname = add2fbody(bname, string('.'//ext%to_char()), string(suffix))
        call bname%kill
        call ext%kill
    end function append2basename_2

    function rm_from_fbody( fname, suffix, str ) result( newname )
        class(*),      intent(in)  :: fname 
        class(string), intent(in)  :: suffix, str 
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp = cast_str_types(fname)
        pos       = index(fname_tmp, str%to_char()) ! position of str
        newname   = fname_tmp(:pos-1)//suffix%to_char()
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function rm_from_fbody

    function swap_suffix_1( fname, suffix, old_suffix ) result( newname )
        class(*),      intent(in)  :: fname
        class(string), intent(in)  :: suffix, old_suffix
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp = cast_str_types(fname)
        pos = index(fname_tmp, old_suffix%to_char()) ! position of old_suffix
        newname = fname_tmp(:pos-1)//suffix%to_char()
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function swap_suffix_1

    function swap_suffix_2( fname, suffix, old_suffix ) result( newname )
        class(*),         intent(in)  :: fname
        character(len=*), intent(in)  :: suffix, old_suffix
        character(len=:), allocatable :: fname_tmp
        type(string) :: newname
        integer :: pos
        fname_tmp = cast_str_types(fname)
        pos = index(fname_tmp, trim(adjustl(old_suffix))) ! position of old_suffix
        newname = fname_tmp(:pos-1)//trim(adjustl(suffix))
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function swap_suffix_2

    function get_fbody_1( fname, suffix, separator ) result( fbody )
        class(*),          intent(in) :: fname
        class(string),     intent(in) :: suffix
        logical, optional, intent(in) :: separator
        character(len=:), allocatable :: fname_tmp
        type(string) :: fbody
        integer      :: pos
        logical      :: l_separator
        l_separator = .true.
        if(present(separator))l_separator = separator
        fname_tmp = cast_str_types(fname)
        if( l_separator )then
            pos = index(fname_tmp, '.'//suffix%to_char()) ! position of suffix
        else
            pos = index(fname_tmp, suffix%to_char())      ! position of suffix
        endif
        fbody = fname_tmp(:pos-1)
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function get_fbody_1

    function get_fbody_2( fname, suffix, separator ) result( fbody )
        class(*),          intent(in) :: fname
        character(len=*),  intent(in) :: suffix
        logical, optional, intent(in) :: separator
        character(len=:), allocatable :: fname_tmp
        type(string) :: fbody
        integer      :: pos
        logical      :: l_separator
        l_separator = .true.
        if(present(separator))l_separator = separator
        fname_tmp = cast_str_types(fname)
        if( l_separator )then
            pos = index(fname_tmp, '.'//trim(adjustl(suffix))) ! position of suffix
        else
            pos = index(fname_tmp, trim(adjustl(suffix)))      ! position of suffix
        endif
        fbody = fname_tmp(:pos-1)
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function get_fbody_2

    function fname_new_ext_1( fname, suffix ) result( new_fname )
        class(string), intent(in) :: fname, suffix
        type(string) :: ext, fbody, new_fname
        ext       = fname2ext(fname)
        fbody     = get_fbody(fname, ext)
        new_fname = fbody%to_char()//'.'//suffix%to_char()
    end function fname_new_ext_1

    function fname_new_ext_2( fname, suffix ) result( new_fname )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: suffix
        type(string) :: ext, fbody, new_fname
        ext       = fname2ext(fname)
        fbody     = get_fbody(fname, ext)
        new_fname = fbody%to_char()//'.'//trim(adjustl(suffix))
    end function fname_new_ext_2

    function fname2ext( fname ) result( ext )
        class(string),    intent(in)  :: fname
        character(len=:), allocatable :: fname_tmp
        type(string) :: ext
        integer      :: length, pos
        length    = fname%strlen_trim()
        fname_tmp = fname%to_char()
        pos = scan(fname_tmp(1:length),'.',back=.true.)
        if( pos == 0 )then
            ext = ''
        else
            ext = trim(fname_tmp(pos+1:length))
        endif
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function fname2ext

    integer function fname2iter( fname )
        class(string),    intent(in)  :: fname
        character(len=:), allocatable :: iter_num_ext, nrstr, fname_tmp
        logical,          allocatable :: lnrs(:)
        integer :: ind, i
        fname_tmp = fname%to_char()
        ind = index(fname_tmp, 'iter')
        if(ind == 0)then
            fname2iter = 0
            return
        endif
        iter_num_ext = fname_tmp(ind:)
        lnrs         = map_str_nrs(iter_num_ext)
        nrstr        = ''
        do i=1,len_trim(iter_num_ext)
            if( lnrs(i) ) nrstr = nrstr//iter_num_ext(i:i)
        end do
        fname2iter = str2int(nrstr)
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function fname2iter

    ! thread-safe version of basename
    function basename( fname ) result( new_fname)
        class(string),    intent(in)  :: fname
        character(len=:), allocatable :: fname_tmp
        type(string) :: new_fname
        integer      :: length, pos
        length    = fname%strlen_trim()
        fname_tmp = fname%to_char()
        pos = scan(fname_tmp(1:length),PATH_SEPARATOR,back=.true.)
        if( pos == 0 )then
            new_fname = trim(fname_tmp)
        else
            new_fname = trim(fname_tmp(pos+1:length))
        endif
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function basename

    function stemname( fname ) result( new_fname)
        class(string),    intent(in)  :: fname
        character(len=:), allocatable :: fname_tmp
        type(string) :: new_fname
        integer :: length, pos
        length    = fname%strlen_trim()
        fname_tmp = fname%to_char()
        pos       = scan(fname_tmp(1:length),PATH_SEPARATOR,back=.true.)
        if(pos == length)then !< case with trailling slash
            pos = scan(fname_tmp(1:length-1),PATH_SEPARATOR,back=.true.)
        end if
        if( pos == 0 )then
            new_fname = trim(fname_tmp)
        else
            new_fname = trim(fname_tmp(1:pos - 1))
        endif
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function stemname

    function get_fpath( fname ) result( path )
        class(string),    intent(in)  :: fname
        character(len=:), allocatable :: fname_tmp
        type(string) :: path
        integer      :: pos
        fname_tmp = fname%to_char()
        pos       = scan(fname_tmp,PATH_SEPARATOR,back=.true.)
        if( pos == 0 )then
            path = PATH_HERE
        else
            path = trim(fname_tmp(:pos))
        endif
        if( allocated(fname_tmp) ) deallocate(fname_tmp)
    end function get_fpath

    !>  \brief  returns numbered names (body) with 0-padded integer strings
    function make_dirnames( body, n, numlen ) result( names )
        character(len=*),  intent(in) :: body
        integer,           intent(in) :: n
        integer, optional, intent(in) :: numlen
        type(string), allocatable :: names(:)
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
        character(len=*),           intent(in) :: body, ext
        integer,                    intent(in) :: n
        integer,          optional, intent(in) :: numlen
        character(len=*), optional, intent(in) :: suffix
        type(string), allocatable :: names(:)
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
    function filepath_1(s1, s2, s3) result( fname )
        class(string),           intent(in) :: s1, s2
        class(string), optional, intent(in) :: s3
        type(string) :: fname
        if( present(s3) )then
            fname = s1%to_char()//PATH_SEPARATOR//s2%to_char()//PATH_SEPARATOR//s3%to_char()
        else
            fname = s1%to_char()//PATH_SEPARATOR//s2%to_char()
        endif
    end function filepath_1

    function filepath_2(s1, s2) result( fname )
        class(string),    intent(in) :: s1
        character(len=*), intent(in) :: s2
        type(string) :: fname
        fname = s1%to_char()//PATH_SEPARATOR//trim(adjustl(s2))
    end function filepath_2

    function filepath_3(s1, s2) result( fname )
        character(len=*), intent(in) :: s1
        character(len=*), intent(in) :: s2
        type(string) :: fname
        fname = trim(adjustl(s1))//PATH_SEPARATOR//trim(adjustl(s2))
    end function filepath_3

    function filepath_4(s1, s2) result( fname )
        character(len=*), intent(in) :: s1
        class(string),    intent(in) :: s2
        type(string) :: fname
        fname = trim(adjustl(s1))//PATH_SEPARATOR//s2%to_char()
    end function filepath_4

    !> \brief  is for deleting consecutively numbered files with padded number strings
    subroutine del_files_1( body, n, ext, numlen, suffix )
        character(len=*),           intent(in) :: body !< input filename body
        integer,                    intent(in) :: n    !< total num for del, formatted as body[n].ext
        character(len=*), optional, intent(in) :: ext  !< input filename extension
        integer,          optional, intent(in) :: numlen !< number length
        character(len=*), optional, intent(in) :: suffix !< file suffix
        type(string), allocatable :: names(:)
        integer :: ifile
        if( present(ext) )then
            names = make_filenames( body, n, ext, numlen=numlen, suffix=suffix )
        else
            names = make_dirnames( body, n, numlen=numlen)
        endif
        do ifile=1,n
            if( file_exists(names(ifile)) ) call del_file(names(ifile))
        end do
    end subroutine del_files_1

    subroutine del_files_2( list )
        class(string), intent(in) :: list(:)
        integer :: ifile, n
        n = size(list)
        do ifile=1,n
            if( file_exists(list(ifile)) ) call del_file(list(ifile))
        end do
    end subroutine del_files_2

    subroutine move_files_in_cwd( dir, files_that_stay )
        class(string),           intent(in) :: dir
        class(string), optional, intent(in) :: files_that_stay(:)
        type(string), allocatable :: file_list(:)
        integer :: n_stay, n_move, i, j
        logical :: l_move
        n_stay = 0
        if( present(files_that_stay) )then
            n_stay = size(files_that_stay)
        endif
        call simple_list_files('*', file_list)
        n_move = 0
        if( allocated(file_list) ) n_move = size(file_list)
        if( n_move > 0 ) call simple_mkdir(dir)
        if( n_stay > 0 .and. n_move > 0 )then
            do i = 1, n_move
                l_move = .true.
                do j = 1, n_stay
                    if( file_list(i)%has_substr(files_that_stay(j)) )then
                        l_move = .false.
                        exit
                    endif
                end do
                if( l_move ) call simple_rename(file_list(i), dir//file_list(i))
            end do
        else if( n_move > 0 )then
            do i = 1, n_move
                call simple_rename(file_list(i), dir//file_list(i))
            end do
        endif
        if( n_move > 0 ) call file_list%kill
    end subroutine move_files_in_cwd

    subroutine move_files2dir( dir, file_list )
        class(string), intent(in) :: dir, file_list(:)
        integer :: n_move, i
        n_move = size(file_list)
        if( n_move > 0 )then
            do i = 1, n_move
                call simple_rename(file_list(i), dir//file_list(i))
            end do
        endif
    end subroutine move_files2dir

    !>  \brief Return a one letter code for the file format designated by the extension in the fname
    !!         if .mrc: M
    !!         if .spi: S
    !!         if .img: I
    !!         if .hed: I
    !!         else: N
    function fname2format( fname )
        class(string),    intent(in)  :: fname
        character(len=1)              :: fname2format
        type(string) :: extension
        extension = fname2ext(fname)
        select case(extension%to_char())
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
            case('eer')
                fname2format = 'K'
            case('gain')
                fname2format = 'L'
#endif
            case DEFAULT
                fname2format = 'N'
        end select
    end function fname2format

    subroutine read_filetable( filetable, filenames )
        class(string),             intent(in)    :: filetable    !< input table filename
        type(string), allocatable, intent(inout) :: filenames(:) !< array of filenames
        integer :: nl, funit, iline, io_stat
        nl = nlines(filetable)
        if( nl == 0 ) return
        call fopen(funit, filetable, 'old', 'unknown', io_stat)
        call fileiochk("read_filetable failed to open file "//filetable%to_char(), io_stat)
        if( allocated(filenames) )then
            call filenames%kill
            deallocate(filenames)
        endif
        allocate(filenames(nl))
        do iline = 1, nl
            call filenames(iline)%readline(funit, io_stat)
        end do
        call fclose(funit)
    end subroutine read_filetable

    !>  \brief  writes a filetable array to a text file
    subroutine write_filetable( filetable, filenames )
        class(string),  intent(in) :: filetable    !< output table filename
        class(string),  intent(in) :: filenames(:) !< array of filenames
        type(string) :: str_line
        integer      :: nl, funit, iline, io_stat
        nl = size(filenames)
        call fopen(funit,filetable, 'replace', 'unknown', io_stat)
        call fileiochk("write_filetable failed to open file "//filetable%to_char(),io_stat )
        do iline=1,nl
            str_line = simple_abspath(filenames(iline))
            call str_line%writeline(funit, io_stat)
            if( io_stat /= 0 ) THROW_HARD('writing of line: '//int2str(iline)//' failed')
            call str_line%kill
        end do
        call fclose(funit)
        call fileiochk("write_filetable failed to close",io_stat)
    end subroutine write_filetable

    !>  \brief  (over)writes a file with a single file of text
    subroutine write_singlelineoftext( filename, text )
        class(string), intent(in)  :: filename, text
        integer :: funit, io_stat
        call fopen(funit, filename, 'replace', 'unknown', io_stat)
        call fileiochk("write_singlelineoftext failed to open file: "//filename%to_char(),io_stat )
        write(funit,'(A)') text%to_char()
        call fclose(funit)
    end subroutine write_singlelineoftext

    !>  \brief  read exit code from file. Format is a single integer on a one line file
    subroutine read_exit_code( filename, exit_code, err )
        class(string), intent(in)  :: filename
        integer,       intent(out) :: exit_code
        logical,       intent(out) :: err        ! error dealing with file, unrelated to exit code
        integer :: nl, funit, io_stat
        exit_code = -1
        err       = .true.
        nl = nlines(filename)
        if( nl /= 1 ) return
        call fopen(funit, filename, 'old', 'unknown', io_stat)
        call fileiochk("read_exit_code failed to open file "//filename%to_char(), io_stat,die=.false.)
        err = io_stat /= 0
        if( .not.err )then
            read(funit,'(i32)') exit_code
        endif
        call fclose(funit)
    end subroutine read_exit_code

    !> \brief  for converting a real array 2 file
    subroutine arr2file_sp( arr, fnam )
        real,          intent(in) :: arr(:)
        class(string), intent(in) :: fnam      !< input table filename
        real    :: rval
        integer :: recsz, i, funit,io_stat
        inquire(iolength=recsz)rval
        rval = size(arr)
        call fopen(funit,fnam,'replace','unknown', iostat=io_stat,access='direct',form='unformatted',recl=recsz)
        call fileiochk("arr2file_sp fopen failed "//fnam%to_char(),io_stat)
        write(funit, rec=1, iostat=io_stat) rval
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        call fclose(funit)
    end subroutine arr2file_sp

    !> \brief  for converting a real array 2 file
    subroutine arr2file_dp( arr, fnam )
        real(kind=dp), intent(in) :: arr(:)
        class(string), intent(in) :: fnam      !< input table filename
        real(dp) :: rval
        integer  :: recsz, i, funit,io_stat
        inquire(iolength=recsz)rval
        rval = size(arr)
        call fopen(funit,fnam,'replace','unknown', iostat=io_stat,access='direct',form='unformatted',recl=recsz)
        call fileiochk("arr2file_dp fopen failed "//fnam%to_char(),io_stat)
        write(funit, rec=1, iostat=io_stat) rval
        do i=1,size(arr)
            write(funit, rec=i+1) arr(i)
        end do
        call fclose(funit)
    end subroutine arr2file_dp

    !> \brief  for writing real array
    subroutine arr2txtfile_1( arr, fnam )
        real,          intent(in) :: arr(:)
        class(string), intent(in) :: fnam
        integer :: i,funit,io_stat
        call fopen(funit,fnam,'replace', 'unknown', iostat=io_stat, form='formatted')
        call fileiochk("arr2txtfile fopen failed "//fnam%to_char(),io_stat)
        do i=1,size(arr)
            write(funit,*) arr(i)
        end do
        call fclose(funit)
    end subroutine arr2txtfile_1

    !> \brief  for writing integer array
    subroutine arr2txtfile_2( arr, fnam )
        integer,       intent(in) :: arr(:)
        class(string), intent(in) :: fnam
        integer :: i,funit,io_stat
        call fopen(funit,fnam,'replace', 'unknown', iostat=io_stat, form='formatted')
        call fileiochk("arr2txtfile fopen failed "//fnam%to_char(),io_stat)
        do i=1,size(arr)
            write(funit,*) arr(i)
        end do
        call fclose(funit)
    end subroutine arr2txtfile_2

    !> \brief  for converting a file generated by arr2file_sp back to an array
    function file2rarr( fnam ) result( arr )
        class(string), intent(in) :: fnam  !< input table filename
        real, allocatable         :: arr(:)
        real    :: rval
        integer :: recsz, i, n, funit,io_stat
        if( file_exists(fnam) )then
            inquire(iolength=recsz) rval
            call fopen(funit, fnam,'old','unknown', io_stat,'direct','unformatted',recl=recsz)
            call fileiochk("file2rarr fopen failed "//fnam%to_char(),io_stat)
            read(funit, rec=1,iostat=io_stat) rval
            n = nint(rval)
            allocate( arr(n) )
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            call fclose(funit)
        else
            THROW_HARD(fnam%to_char()//' does not exist')
        endif
    end function file2rarr

    !> \brief  for converting a file generated by arr2file_dp back to an array
    function file2drarr( fnam ) result( arr )
        class(string), intent(in) :: fnam  !< input table filename
        real(dp), allocatable     :: arr(:)
        real(dp) :: rval
        integer  :: recsz, i, n, funit,io_stat
        if( file_exists(fnam) )then
            inquire(iolength=recsz) rval
            call fopen(funit,fnam,'old','unknown', io_stat,'direct','unformatted',recl=recsz)
            call fileiochk("file2drarr fopen failed "//fnam%to_char(),io_stat)
            read(funit, rec=1,iostat=io_stat) rval
            n = nint(rval)
            allocate( arr(n) )
            do i=1,n
                read(funit, rec=i+1) arr(i)
            end do
            call fclose(funit)
        else
            THROW_HARD(fnam%to_char()//' does not exist')
        endif
    end function file2drarr

    !> \brief  for converting a real matrix to file
    ! Can be imported into python with:
    ! import numpy as np
    ! f    = open('somename','rb')
    ! dims = np.fromfile(f, dtype=np.int32, count=2)
    ! A    = np.fromfile(f, dtype=np.float32, count=np.prod(dims)).reshape(dims, order='F')
    ! f.close()
    subroutine rmat2file( mat, fname )
        real,          intent(in) :: mat(:,:)    !< Input matrix
        class(string), intent(in) :: fname       !< Output filename
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("mat2file fopen failed: "//fname%to_char(),io_stat)
        write(unit=funit,pos=1) shape(mat)
        write(unit=funit,pos=(2*sizeof(funit)+1)) mat
        call fclose(funit)
    end subroutine rmat2file

    !> \brief  for converting a logical matrix to file
    subroutine lmat2file( mat, fname )
        logical,       intent(in) :: mat(:,:)    !< Input matrix
        class(string), intent(in) :: fname       !< Output filename
        integer :: funit,io_stat
        call fopen(funit, fname, access='STREAM', action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk("lmat2file fopen failed: "//fname%to_char(),io_stat)
        write(unit=funit,pos=1) shape(mat)
        write(unit=funit,pos=(2*sizeof(funit)+1)) mat
        call fclose(funit)
    end subroutine lmat2file

    !> \brief  for converting a file generated by rmat2file back to a matrix
    subroutine file2rmat( fname, rmat )
        class(string),     intent(in)    :: fname       !< input filename
        real, allocatable, intent(inout) :: rmat(:,:)   !< output matrix
        integer :: dims(2), funit,io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('file2rmat; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        if( allocated(rmat) ) deallocate(rmat)
        allocate(rmat(dims(1),dims(2)))
        read(unit=funit,pos=(sizeof(dims)+1)) rmat
        call fclose(funit)
    end subroutine file2rmat

    !> \brief  for converting a file generated by lmat2file back to a logical matrix
    subroutine file2lmat( fname, mat )
        class(string),        intent(in)    :: fname    !< input filename
        logical, allocatable, intent(inout) :: mat(:,:) !< output matrix
        integer :: dims(2), funit,io_stat
        if( .not.file_exists(fname) ) THROW_HARD(fname%to_char()//' does not exist')
        call fopen(funit, fname, access='STREAM', action='READ', status='OLD', iostat=io_stat)
        call fileiochk('file2lmat; fopen failed: '//fname%to_char(), io_stat)
        read(unit=funit,pos=1) dims
        if( allocated(mat) ) deallocate(mat)
        allocate(mat(dims(1),dims(2)))
        read(unit=funit,pos=(sizeof(dims)+1)) mat
        call fclose(funit)
    end subroutine file2lmat

    subroutine simple_copy_file(fname1, fname2, status)
        class(string),     intent(in)  :: fname1, fname2 !< input filenames
        integer, optional, intent(out) :: status
        integer(dp),       parameter   :: MAXBUFSZ = nint(1e8) ! 100 MB max buffer size
        character(len=1),  allocatable :: byte_buffer(:)
        integer(dp) :: sz, nchunks, leftover, bufsz, bytepos, in, out, ichunk
        integer     :: ioerr
        if( present(status) )status = 0
        if( .not. file_exists(fname1) ) THROW_HARD('file '//fname1%to_char()//' does not exist')
        ! we need to inquire size before opening file as stream access
        ! does not allow inquire of size from file unit
        inquire(file=fname1%to_char(),size=sz)
        open(newunit=in, file=fname1%to_char(), status="old", action="read", access="stream", iostat=ioerr)
        if( ioerr /= 0 )then
            write(logfhandle,*) "In simple_copy_file, failed to open input file ", fname1%to_char()
            call simple_error_check(ioerr,"simple_copy_file input file not opened")
            if( present(status) ) status = ioerr
            return
        endif
        if( sz <= MAXBUFSZ )then
            nchunks  = 1
            leftover = 0
            bufsz    = sz
        else
            nchunks  = sz / MAXBUFSZ
            leftover = sz - nchunks * MAXBUFSZ
            bufsz    = MAXBUFSZ
        endif
        ! allocate raw byte buffer
        allocate(byte_buffer(bufsz))
        ! process output file
        open(newunit=out, file=fname2%to_char(), status="replace", action="write", access="stream", iostat=ioerr)
        if( ioerr /= 0 )then
            write(logfhandle,*)"In simple_copy_file, failed to open output file ", fname2%to_char()
            call simple_error_check(ioerr,"simple_copy_file output file not opened")
            if( present(status) ) status = ioerr
            return
        endif
        bytepos = 1
        do ichunk=1,nchunks
            read(in,   pos=bytepos, iostat=ioerr) byte_buffer
            if( ioerr /= 0 ) THROW_HARD("failed to read byte buffer")
            write(out, pos=bytepos, iostat=ioerr) byte_buffer
            if( ioerr /= 0 ) THROW_HARD("failed to write byte buffer")
            bytepos = bytepos + bufsz
        end do
        ! take care of leftover
        if( leftover > 0 )then
            read(in,   pos=bytepos, iostat=ioerr) byte_buffer(:leftover)
            if( ioerr /= 0 ) THROW_HARD("failed to read byte buffer")
            write(out, pos=bytepos, iostat=ioerr) byte_buffer(:leftover)
            if( ioerr /= 0 ) THROW_HARD("failed to write byte buffer")
        endif
        close(in)
        close(out)
    end subroutine simple_copy_file

    function get_relative_path ( path, root, trimlength ) result ( newpath )
        class(string),     intent(in)    :: path, root
        integer, optional, intent(inout) :: trimlength
        integer      :: l_trimlength
        type(string) :: newpath
        l_trimlength = 0
        if(present(trimlength)) l_trimlength = trimlength 
        if(l_trimlength .eq. 0) l_trimlength = path%substr_ind(root)
        if(l_trimlength .gt. 0) then
            newpath = string(path%to_char([l_trimlength,path%strlen_trim()]))
        else
            newpath = path
        end if
        if(present(trimlength)) trimlength = l_trimlength
    end function get_relative_path

end module simple_fileio


