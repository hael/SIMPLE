module simple_syslib
use simple_defs
use simple_error
use simple_string
use simple_string_utils
use, intrinsic :: iso_fortran_env
use, intrinsic :: iso_c_binding
implicit none
! local version of throw_hard to enable public feature
#define THROW_ERROR(msg) call simple_exception(msg, __FILENAME__ , __LINE__)

!> glibc interface CONFORMING TO POSIX.1-2001, POSIX.1-2008, SVr4, 4.3BSD.
interface

    function rmdir(dirname) bind(C, name="rmdir")
        use, intrinsic :: iso_c_binding
        integer(c_int) :: rmdir
        character(c_char),dimension(*),intent(in)  ::  dirname
    end function rmdir

    !! mkdir() attempts to create a directory named pathname. mkdir returns zero
    !! on success, or -1 if an error occurred (in which case, errno is set
    !! appropriately). If errno equals EEXIST pathname already exists (not
    !! necessarily as a directory). This includes the case where pathname is a
    !! symbolic link, dangling or not.
    function mkdir(path,mode) bind(c,name="mkdir")
        use, intrinsic :: iso_c_binding
        integer(c_int) :: mkdir
        character(kind=c_char,len=1),dimension(*),intent(in) :: path
        integer(c_int16_t), value :: mode
    end function mkdir

    !! symlink() creates a symbolic link named linkpath to target. On success,
    !! zero is returned. On error, -1 is returned, and errno is set
    !! appropriately.
    function symlink(target_path, link_path) bind(c,name="symlink")
        use, intrinsic :: iso_c_binding
        integer(c_int) :: symlink
        character(kind=c_char,len=1),dimension(*),intent(in) :: target_path
        character(kind=c_char,len=1),dimension(*),intent(in) :: link_path
    end function symlink

    ! !!  sync() causes all buffered modifications to file metadata and data to be
    ! !!  written to the underlying filesystems.
    function fsync (fd) bind(c,name="fsync")
      use iso_c_binding, only: c_int
        integer(c_int), value :: fd
        integer(c_int) :: fsync
    end function fsync

    ! For float parsing only!
    function sscanf(str, fmt, val) bind(C, name="sscanf")
        use iso_c_binding, only : c_int, c_char, c_float
        integer(kind=c_int) :: sscanf
        character(kind=c_char,len=1), dimension(*),intent(in)  :: str, fmt
        real(kind=c_float),                        intent(out) :: val
    end function sscanf

end interface

!> SIMPLE_POSIX.c commands
interface

    function isdir(dirname, str_len) bind(C, name="isdir")
        import
        integer(c_int) :: isdir
        character(c_char),dimension(*),intent(in)  ::  dirname
        integer(c_int), intent(in) :: str_len
    end function isdir

    function makedir(dirname, str_len) bind(C, name="makedir")
        import
        integer(c_int) :: makedir
        character(c_char),dimension(*),intent(in)  ::  dirname
       integer(c_int), intent(in) :: str_len
    end function makedir

    function removedir(dirname,str_len, count) bind(C, name="remove_directory")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: removedir
        character(c_char),dimension(*),intent(in)  ::  dirname
        integer(c_int), intent(in) :: str_len
        integer(c_int), intent(inout) :: count
    end function removedir

    function list_dirs(path, str_len, list_fout, str_len_fout, count) bind(c,name="list_dirs")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: list_dirs                                 !> return success
        character(kind=c_char,len=1),dimension(*),intent(in):: path !> input pathname
        integer(c_int), intent(in)    :: str_len                    !> input pathname string length
        character(kind=c_char,len=1),dimension(*),intent(in):: list_fout !> output list string file name
        integer(c_int), intent(in)    :: str_len_fout               !> output list file name string length
        integer(c_int), intent(inout) :: count                      !> return number of elements in results
    end function list_dirs

    function wait_pid(pid) bind(c,name="wait_pid")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: wait_pid                                    !> return PID of forked process
        integer(c_int), intent(in) :: pid
    end function wait_pid

    function touch(filename, len) bind(c,name="touch")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: touch                                       !> return success of touch
        character(kind=c_char,len=1),dimension(*),intent(in) :: filename
        integer(c_int), intent(in) :: len
    end function touch

    function get_absolute_pathname(infile, inlen, outfile, outlen) bind(c,name="get_absolute_pathname")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: get_absolute_pathname                             !> return status
        character(kind=c_char,len=1),dimension(*),intent(in)    :: infile   !> input pathname
        integer(c_int), intent(in)  :: inlen                                !> input pathname string length
        character(kind=c_char,len=1),dimension(*),intent(inout) :: outfile  !> output pathname
        integer(c_int), intent(out) :: outlen                               !> output pathname string length
    end function get_absolute_pathname

    function get_sysinfo(HWM, totRAM, shRAM, bufRAM, peakBuf) bind(c,name="get_sysinfo")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: get_sysinfo
        integer(c_long), intent(inout) :: HWM                !> high-water mark
        integer(c_long), intent(inout) :: totRAM             !> total RAM usage
        integer(c_long), intent(inout) :: shRAM              !> shared RAM usage
        integer(c_long), intent(inout) :: bufRAM             !> this process's buffered RAM
        integer(c_long), intent(inout) :: peakBuf            !> this process's peak RAM usage
    end function get_sysinfo

    function regexp_match(source,src_len, regex,rgx_len) bind(C,name="regexp_match")
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: regexp_match
        character(kind=c_char,len=1),dimension(*),intent(in)    :: source  !> input string
        integer(c_int), intent(in) :: src_len
        character(kind=c_char,len=1),dimension(*),intent(in)    :: regex   !> input RE string
        integer(c_int), intent(in) :: rgx_len
    end function

end interface

contains

    subroutine exec_cmdline( cmdline, waitflag, suppress_errors, exitstat)
        class(*),          intent(in)  :: cmdline
        logical, optional, intent(in)  :: waitflag, suppress_errors
        integer, optional, intent(out) :: exitstat
        character(len=:), allocatable  :: cmdstr
        character(len=100) ::errmsg
        integer ::  cstat, exec_stat
        logical :: l_doprint, wwait, l_suppress_errors
        l_doprint = .false.
        wwait     = .true.
        l_suppress_errors = .false.
        if( present(waitflag)        ) wwait = waitflag
        if( present(suppress_errors) ) l_suppress_errors = suppress_errors
        select type(cmdline)
            type is (string)
                if( l_suppress_errors )then
                    allocate(cmdstr, source=cmdline%to_char()//' '//SUPPRESS_MSG)
                else
                    allocate(cmdstr, source=cmdline%to_char())
                endif
            type is (character(*))
                if( l_suppress_errors )then
                    allocate(cmdstr, source=trim(adjustl(cmdline))//' '//SUPPRESS_MSG)
                else
                    allocate(cmdstr, source=trim(adjustl(cmdline)))
                endif
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        call execute_command_line(cmdstr, wait=wwait, exitstat=exec_stat, cmdstat=cstat, cmdmsg=errmsg)
        if( .not. l_suppress_errors ) call raise_sys_error( cmdstr, exec_stat, cstat, errmsg )
        if( l_doprint )then
            write(logfhandle,*) 'command            : ', cmdstr
            write(logfhandle,*) 'status of execution: ', exec_stat
        endif
        if(present(exitstat))exitstat=exec_stat

        contains

            !>  Handles error from system call
            subroutine raise_sys_error( cmd, exit_status, cmdstat, cmdmsg )
                integer,          intent(in) :: exit_status, cmdstat
                character(len=*), intent(in) :: cmd
                character(len=*), intent(in) :: cmdmsg
                logical :: err
                err = .false.
                if( exit_status /= 0 )then
                    write(logfhandle,*)'System error', exit_status,' for command: ', trim(adjustl(cmd))
                    err = .true.
                endif
                if( cmdstat /= 0 )then
                    write(logfhandle,*)cmdmsg
                    call simple_error_check(cmdstat,' command could not be executed: '//trim(adjustl(cmd)))
                    write(logfhandle,*)'cmdstat = ',cmdstat,' command could not be executed: ', trim(adjustl(cmd))
                    err = .true.
                endif
            end subroutine raise_sys_error

    end subroutine exec_cmdline

    function simple_getenv( name , status )  result( envval )
        character(len=*),  intent(in)  :: name
        integer,           intent(out) :: status
        type(string)               :: envval
        character(len=XLONGSTRLEN) :: retval
        integer                    :: length
        call get_environment_variable( trim(adjustl(name)), value=retval, length=length, status=status)
        if( status == -1 ) write(logfhandle,*) 'value string too short; simple_syslib :: simple_getenv_2'
        if( status ==  1 )then
            write(logfhandle,*) 'environment variable: ', trim(adjustl(name)), ' is not defined; simple_syslib :: simple_getenv_'
            envval = NIL
            return
        endif
        if( status ==  2 ) write(logfhandle,*) 'environment variables not supported by system; simple_syslib :: simple_getenv'
        if( length ==  0 .or. status /= 0 )then
            envval = NIL
            return
        end if
        envval = trim(adjustl(retval))
    end function simple_getenv

    !> \brief Touch file, create file if necessary
    subroutine simple_touch( fname, status )
        class(*),          intent(in)  :: fname
        integer, optional, intent(out) :: status
        integer :: iostat
        select type(fname)
            type is (string)
                iostat  = touch(fname%to_char(), fname%strlen_trim())
            type is (character(*))
                iostat  = touch(fname, len_trim(fname))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        if(iostat/=0)then
            call simple_error_check(iostat, "simple_touch")
        endif
        if(present(status))status=iostat
    end subroutine simple_touch

    !> \brief Soft link file
    subroutine syslib_symlink( f1, f2, status)
        class(*),          intent(in)  :: f1, f2
        integer, optional, intent(out) :: status
        integer :: iostat
         select type(f1)
            type is (string)
                select type(f2)
                    type is (string)
                        iostat = symlink(f1%to_char()//achar(0), f2%to_char()//achar(0))
                    type is (character(*))
                        iostat = symlink(f1%to_char()//achar(0), trim(adjustl(f2))//achar(0))
                    class default
                        call simple_exception('Incompatible type, f2 ', __FILENAME__ , __LINE__)
                end select
            type is (character(*))
                 select type(f2)
                    type is (string)
                        iostat = symlink(trim(adjustl(f1))//achar(0), f2%to_char()//achar(0))
                    type is (character(*))
                        iostat = symlink(trim(adjustl(f1))//achar(0), trim(adjustl(f2))//achar(0))
                    class default
                        call simple_exception('Incompatible type, f2 ', __FILENAME__ , __LINE__)
                end select
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        if(iostat/=0)then
            call simple_error_check(iostat, "In syslib_symlink")
        endif
        if(present(status))status=iostat
    end subroutine syslib_symlink

    !> \brief  Rename or move file
    subroutine simple_rename( filein, fileout, overwrite )
        class(*),          intent(in) :: filein, fileout !< input filename
        logical, optional, intent(in) :: overwrite       !< default true
        character(kind=c_char, len=:), allocatable :: f1, f2
        character(len=:),              allocatable :: msg
        logical :: force_overwrite
        integer :: file_status
        if( .not. file_exists(filein) ) THROW_ERROR("designated input file does not exist")
        force_overwrite = .true.
        if( present(overwrite) ) force_overwrite=overwrite
        if( file_exists(fileout) .and. (force_overwrite) ) call del_file(fileout)
        select type(filein)
            type is(string)
                allocate(f1, source=filein%to_char()//achar(0))
            type is(character(*))
                allocate(f1, source=trim(adjustl(filein))//achar(0))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        select type(fileout)
            type is(string)
                allocate(f2, source=fileout%to_char()//achar(0))
            type is(character(*))
                allocate(f2, source=trim(adjustl(fileout))//achar(0))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        file_status = rename(trim(f1), trim(f2))
        if( file_status /= 0 )then
            allocate(msg,source="simple_rename failed to rename file "//f1)
            call simple_error_check(file_status, trim(msg))
            deallocate(msg)
        endif
        deallocate(f1,f2)
    end subroutine simple_rename

    function simple_chmod( pathname, mode ) result( status )
        class(*),         intent(in) :: pathname
        character(len=*), intent(in) :: mode
        integer :: status
        select type(pathname)
            type is(string)
                status = chmod(pathname%to_char(), mode) !! intrinsic GNU
                if( status/=0 ) call simple_error_check(status,"simple_syslib::simple_chmod chmod failed "//pathname%to_char())
            type is(character(*))
                status = chmod(trim(adjustl(pathname)), mode) !! intrinsic GNU
                if( status/=0 ) call simple_error_check(status,"simple_syslib::simple_chmod chmod failed "//trim(adjustl(pathname)))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
    end function simple_chmod

    !>  Wrapper for POSIX system call stat
    subroutine simple_file_stat( filename, status, buffer)
        class(*),             intent(in)    :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        allocate(buffer(13), source=0)
        select type(filename)
            type is(string)
                status = stat(filename%to_char(), buffer)
            type is(character(*))
                status = stat(trim(adjustl(filename)), buffer)
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
    end subroutine simple_file_stat

    logical function is_io(unit)
        integer, intent(in) :: unit
        is_io=.false.
        if (unit == ERROR_UNIT .or. unit == OUTPUT_UNIT .or. unit == INPUT_UNIT) is_io= .true.
    end function is_io

    !>  \brief  check whether a IO unit is currently opened
    logical function is_open( unit_number )
        integer, intent(in)   :: unit_number
        integer               :: io_status
        character(len=STDLEN) :: io_message
        io_status = 0
        is_open=.false.
        inquire(unit=unit_number, opened=is_open,iostat=io_status,iomsg=io_message)
        if(is_iostat_eor(io_status) .or. is_iostat_end(io_status)) return
        if (io_status .ne. 0) then
            write(logfhandle,*) 'is_open: I/O error ', io_status, ': ', trim(adjustl(io_message))
            THROW_ERROR('I/O')
        endif
    end function is_open

    !>  \brief  check if a file exists on disk
    !! return logical true=dir exists, false=dir does not exist
    logical function dir_exists( dname )
        class(*), intent(in) :: dname
        integer :: status
        character(kind=c_char, len=:), allocatable :: d1
        dir_exists=.false.
        select type(dname)
            type is(string)
                allocate(d1,source=dname%to_char()//achar(0))
            type is(character(*))
                allocate(d1,source=trim(adjustl(dname))//achar(0))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        status = isdir(trim(d1), len_trim(d1))
        deallocate(d1)
        if (status == 1) dir_exists = .true.
    end function dir_exists

    !>  \brief  check if a file exists on disk
    !! return logical true=FILE exists, false=FILE does not exist
    logical function file_exists( fname )
        class(*), intent(in) :: fname
        select type(fname)
            type is(string)
                inquire(file=fname%to_char(),      exist = file_exists)
            type is(character(*))
                inquire(file=trim(adjustl(fname)), exist = file_exists)
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
    end function file_exists

    !>  \brief  check whether a file is currently opened
    logical function is_file_open( fname )
        class(string), intent(in) :: fname
        integer                   :: io_status
        character(len=STDLEN)     :: io_message
        io_status = 0
        inquire(file=fname%to_char(), opened=is_file_open,iostat=io_status,iomsg=io_message)
        if (io_status .ne. 0) then
            THROW_ERROR('I/O '//trim(adjustl(io_message)))
        endif
    end function is_file_open

    !> \brief  Get current working directory
    subroutine simple_getcwd( cwd )
        class(string), intent(inout) :: cwd
        character(len=XLONGSTRLEN)    :: cwd_tmp
        integer :: io_status
        io_status = getcwd(cwd_tmp)
        if(io_status /= 0)then
            call simple_error_check(io_status, "syslib:: simple_getcwd failed to get path "//trim(cwd_tmp))
        else
            cwd = trim(adjustl(cwd_tmp))
        endif
    end subroutine simple_getcwd

    !> \brief  Change working directory
    !! return optional status 0=success
    subroutine simple_chdir( newd, status )
        class(*),          intent(in)  :: newd   !< target pathname
        integer, optional, intent(out) :: status
        type(string) :: targetdir
        integer      :: io_status
        logical      :: dir_e
        if(present(status)) status = 1
        targetdir = simple_abspath(newd, check_exists=.true.)
        inquire(file=targetdir%to_char(), EXIST=dir_e, IOSTAT=io_status)
        if(dir_e) then
            io_status = chdir(targetdir%to_char())
            if(io_status /= 0)then
                select case (io_status)
                case (2)  ! ENOENT
                    write (*,*)'The directory '//targetdir%to_char()//' does not exist'
                case (20) ! ENOTDIR
                    write (*,*) targetdir%to_char()//' is not a directory'
                case default
                    write (*,*)'Error with code ', io_status
                end select
                call simple_error_check(io_status, &
                    "syslib:: simple_chdir failed to change path "//targetdir%to_char())
            endif
        else
            THROW_ERROR("directory does not exist")
        endif
        if(present(status)) status = io_status
        call targetdir%kill
    end subroutine simple_chdir

    !> \brief  Make directory
    subroutine simple_mkdir( dir, verbose )
        class(*),          intent(in) :: dir
        logical, optional, intent(in) :: verbose
        logical, parameter :: ignore_here = .false.
        character(kind=c_char, len=:), allocatable :: path
        character(len=:), allocatable :: tmpdir
        type(string) :: tmpdir_str, path_str
        integer      :: io_status, lenstr
        logical      :: l_verbose
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        select type(dir)
            type is(string)
                tmpdir = dir%to_char()
            type is(character(*))
                tmpdir = trim(adjustl(dir))
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        lenstr = len_trim(tmpdir)
        if( lenstr==0 ) then
            if( l_verbose )write(logfhandle,*)"syslib:: simple_mkdir arg empty "//trim(tmpdir)
            return
        else if( (lenstr<=2) .and. (tmpdir(1:1)=='/' .or. tmpdir(1:1)=='.') )then
            ! ignore '/' '.' './' '..'
            if( l_verbose )write(logfhandle,*)"syslib:: simple_mkdir arg special char: "//trim(tmpdir)
        endif
        io_status = 0
        tmpdir_str = trim(adjustl(tmpdir))
        if(.not. dir_exists(tmpdir_str)) then
            ! prepare path for C function
            allocate(path, source=trim(tmpdir)//c_null_char)
            io_status = makedir(trim(adjustl(path)), len_trim(tmpdir))
            path_str = trim(adjustl(path))
            if(.not. dir_exists(path_str)) then
                if( l_verbose )write(logfhandle,*)" syslib:: simple_mkdir failed to create "//trim(path)
                if(.not. ignore_here)then
                    if(io_status /= 0) call simple_error_check(io_status, &
                        "syslib:: simple_mkdir failed to create "//trim(path))
                endif
            endif
            deallocate(path)
        end if
        call tmpdir_str%kill
        call path_str%kill
    end subroutine simple_mkdir

    !> \brief  Remove directory
    !! return status 0=success for directory exists or directory created
    !! return error status for other removedir results
    subroutine simple_rmdir( d, status )
        class(*),          intent(in)  :: d
        integer, optional, intent(out) :: status
        character(kind=c_char,len=:), allocatable :: path
        integer :: io_status, err, length, count
        logical :: dir_e
        dir_e = .false.
        select type(d)
            type is(string)
                inquire(file=d%to_char(), exist=dir_e)
                if(dir_e) then
                    allocate(path, source=d%to_char()//c_null_char)
                endif
            type is(character(*))
                inquire(file=trim(adjustl(d)), exist=dir_e)
                if(dir_e) then
                    allocate(path, source=trim(adjustl(d))//c_null_char)
                endif
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        io_status=0
        count=0
        if(dir_e)then
            length = len_trim(adjustl(path))
            io_status = removedir(trim(adjustl(path)), length, count)
            if(io_status /= 0)then
                err = int(IERRNO(), kind=4 )
                call simple_error_check(io_status, "syslib:: simple_rmdir failed to remove "//trim(adjustl(path)))
            endif
            deallocate(path)
        end if
        if(present(status)) status = io_status
    end subroutine simple_rmdir

    !> ensure C-strings get converted to fortran-style strings
    subroutine syslib_c2fortran_string( str, len )
        character(len=*), intent(inout) :: str
        integer, intent(out), optional :: len
        integer :: l
        l = index(str, char(0))
        if(present(len)) len = l-1
        if(l>0) str(l:)=' '
    end subroutine syslib_c2fortran_string

    function find_next_int_dir_prefix( dir2list, last_prev_dir ) result( next_int_dir_prefix )
        class(string),           intent(in)  :: dir2list
        class(string), optional, intent(out) :: last_prev_dir
        character(len=:), allocatable :: str
        type(string),     allocatable :: dirs(:)
        logical,          allocatable :: nrmap(:)
        integer,          allocatable :: dirinds(:)
        integer :: i, j, last_nr_ind
        integer :: next_int_dir_prefix, ndirs, loc(1)
        dirs = simple_list_dirs(dir2list)
        last_nr_ind = 1
        if( allocated(dirs) )then
            ndirs = size(dirs)
        else
            next_int_dir_prefix = 1
            return
        endif
        allocate(dirinds(ndirs), source=0)
        do i=1,ndirs
            str = dirs(i)%to_char()
            if( char_is_a_number(str(1:1)) )then
                nrmap = map_str_nrs(trim(str))
                do j=1,size(nrmap)
                    if( nrmap(j) )then
                        last_nr_ind = j
                    else
                        exit
                    endif
                enddo
                dirinds(i) = str2int(str(1:last_nr_ind))
            endif
        end do
        if( any(dirinds > 0) )then
            loc = maxloc(dirinds)
            next_int_dir_prefix = dirinds(loc(1)) + 1
            if( present(last_prev_dir) ) last_prev_dir = dirs(loc(1))
        else
            next_int_dir_prefix = 1
        endif
    end function find_next_int_dir_prefix

    function simple_list_dirs( path, status ) result( list )
        class(*),          intent(in)  :: path
        integer, optional, intent(out) :: status
        type(string),                 allocatable :: list(:)
        character(kind=c_char,len=:), allocatable :: pathhere
        character(len=XLONGSTRLEN) :: buffer
        type(string) :: list_fname
        integer      :: stat, i,num_dirs, luntmp, pid
        select type(path)
            type is(string)
                allocate(pathhere, source=path%to_char()//c_null_char)
            type is(character(*))
                allocate(pathhere, source=trim(adjustl(path))//c_null_char)
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        pid        = getpid()
        list_fname =  '__simple_dirlist_'//int2str(part_glob)//'_'//int2str(pid)//'__'
        stat       = list_dirs(trim(pathhere),len_trim(pathhere), list_fname%to_char(), list_fname%strlen_trim(), num_dirs)
        if( stat /= 0 )THROW_ERROR("failed to process list_dirs "//trim(pathhere))
        open(newunit=luntmp, file=list_fname%to_char())
        allocate(list(num_dirs))
        do i = 1,num_dirs
            read(luntmp, '(a)') buffer
            list(i) = trim(adjustl(buffer))
        enddo
        close(luntmp)
        call del_file(list_fname)
        deallocate(pathhere)
        if(present(status)) status= stat
    end function simple_list_dirs

    subroutine simple_list_files( pattern, list )
        character(len=*),          intent(in)    :: pattern
        type(string), allocatable, intent(inout) :: list(:)
        character(len=XLONGSTRLEN) :: buffer
        type(string)               :: cmd, tmpfile
        character(len=1)           :: junk
        integer                    :: sz, funit, ios, i, nlines, pid
        pid     = getpid()
        tmpfile = '__simple_filelist_'//int2str(pid)//'__'
        cmd     = 'ls -1f '//trim(pattern)//' > '//tmpfile%to_char()
        call exec_cmdline( cmd, suppress_errors=.true.)
        inquire(file=tmpfile%to_char(), size=sz)
        if( allocated(list) ) deallocate(list)
        if( sz > 0 )then
            open(newunit=funit, file=tmpfile%to_char())
            nlines = 0
            do
                read(funit,*,iostat=ios) junk
                if(ios /= 0)then
                    exit
                else
                    nlines = nlines + 1
                endif
            end do
            rewind(funit)
            allocate( list(nlines) )
            do i=1,nlines
                read(funit, '(a)') buffer
                list(i) = trim(adjustl(buffer))
            enddo
            close(funit, status='delete')
        else
            open(newunit=funit, file=tmpfile%to_char())
            close(funit, status='delete')
        endif
    end subroutine simple_list_files

    subroutine simple_list_files_regexp( dir, regexp, list, chronological )
        class(string),             intent(in)    :: dir
        character(len=*),          intent(in)    :: regexp
        type(string), allocatable, intent(inout) :: list(:)
        logical, optional,         intent(in)    :: chronological
        type(string)                  :: cmd, tmpfile
        character(len=1)              :: junk
        character(len=XLONGSTRLEN)    :: buffer
        integer :: sz, funit, ios, i, nlines
        logical :: l_chrono
        if( len_trim(adjustl(regexp)) == 0) return
        l_chrono = .false.
        if( present(chronological) ) l_chrono = chronological
        tmpfile = '__simple_filelist_'//int2str(part_glob)//'__'
        ! builds command
        if( l_chrono )then
            cmd = 'ls -1f -rt '//dir%to_char()//' | grep -E '''//adjustl(trim(regexp))//''' > '//tmpfile%to_char()
        else
            cmd = 'ls -1f '//dir%to_char()//' | grep -E '''//adjustl(trim(regexp))//''' > '//tmpfile%to_char()
        endif
        call exec_cmdline(cmd, suppress_errors=.true.)
        inquire(file=tmpfile%to_char(), size=sz)
        if( allocated(list) ) deallocate(list)
        if( sz > 0 )then
            open(newunit=funit, file=tmpfile%to_char())
            nlines = 0
            do
                read(funit,*,iostat=ios) junk
                if(ios /= 0)then
                    exit
                else
                    nlines = nlines + 1
                endif
            end do
            rewind(funit)
            allocate( list(nlines) )
            do i=1,nlines
                read(funit, '(a)') buffer
                list(i) = dir%to_char()//'/'//trim(adjustl(buffer))
            enddo
            close(funit, status='delete')
        else
            open(newunit=funit, file=tmpfile%to_char())
            close(funit, status='delete')
        endif
    end subroutine simple_list_files_regexp

    !> \brief  is for deleting a file
    subroutine del_file( file )
        class(*), intent(in) :: file !< input filename
        integer :: fnr, file_status
        if( file_exists(file) )then
            select type(file)
                type is(string)
                    open(newunit=fnr,file=file%to_char(),STATUS='OLD',IOSTAT=file_status)
                    if( file_status == 0 )then
                        close(fnr, status='delete',IOSTAT=file_status)
                        if(file_status /=0) THROW_ERROR("failed to close file "//file%to_char())
                    end if
                type is(character(*))
                     open(newunit=fnr,file=trim(adjustl(file)),STATUS='OLD',IOSTAT=file_status)
                    if( file_status == 0 )then
                        close(fnr, status='delete',IOSTAT=file_status)
                        if(file_status /=0) THROW_ERROR("failed to close file "//trim(adjustl(file)))
                    end if
                class default
                    call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
            end select
        endif
    end subroutine del_file

    integer(4) function get_process_id( )
        get_process_id = getpid()
    end function get_process_id

    function simple_abspath( infile, status, check_exists ) result( absolute_name )
        class(*),          intent(in)  :: infile
        integer, optional, intent(out) :: status
        logical, optional, intent(in)  :: check_exists
        character(kind=c_char,len=XLONGSTRLEN) :: infilename_c
        character(kind=c_char,len=XLONGSTRLEN) :: outfilename_c
        type(string) :: absolute_name
        integer      :: lengthin,  lengthout, status_here
        logical      :: check_exists_here
        check_exists_here = .true.
        if( present(check_exists) )check_exists_here = check_exists
        select type(infile)
            type is(string)
                if( check_exists_here )then
                    if( .not.file_exists(infile) )then
                        write(logfhandle,*)' cwd: '//trim(CWD_GLOB)
                        THROW_ERROR('input file: '//infile%to_char()//' does not exist')
                    endif
                endif
                lengthin     = infile%strlen_trim()
                infilename_c = infile%to_char()//achar(0)
            type is(character(*))
                if( check_exists_here )then
                    if( .not.file_exists(infile) )then
                        write(logfhandle,*)' cwd: '//trim(CWD_GLOB)
                        THROW_ERROR('input file: '//trim(infile)//' does not exist')
                    endif
                endif
                lengthin     = len_trim(adjustl(infile))
                infilename_c = trim(adjustl(infile))//achar(0)
            class default
                call simple_exception('Unsupported type', __FILENAME__ , __LINE__)
        end select
        status_here = get_absolute_pathname(infilename_c, lengthin, outfilename_c, lengthout)
        call syslib_c2fortran_string(outfilename_c)
        if( lengthout > 1)then
            absolute_name = trim(outfilename_c(1:lengthout))
        else
            absolute_name = infile
        end if
        if(present(status))status = status_here
    end function simple_abspath

    subroutine print_slurm_env()
        character(len=255) :: env_value
        character(len=63), allocatable :: env_vars(:);
        integer :: i, stat, len
        len = 63
        allocate (env_vars(7))
        env_vars = [character(len=63) :: "slurm_jobid", "slurm_job_user", "slurm_job_cpus_per_node", "slurm_mem_per_cpu", "slurmd_nodename", "slurm_job_account", "slurm_submit_dir"]
        call get_environment_variable("slurm_jobid", env_value, len, stat)
        if( stat .eq. 0 )then
            write(logfhandle,*) ""
            write(logfhandle,*) "##### simple slurm env #####"
            write(logfhandle,*) ""
            do i = 1, 7
                call get_environment_variable(trim(env_vars(i)), env_value, len, stat)
                if(stat .eq. 0) then
                    write(logfhandle,*) trim(env_vars(i)), achar(9), " : ", achar(9), trim(env_value)
                    endif
                end do
                write(logfhandle,*) ""
                write(logfhandle,*) "############################"
                write(logfhandle,*) ""
        endif
        deallocate(env_vars)
    end subroutine print_slurm_env

end module simple_syslib
