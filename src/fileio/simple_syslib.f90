!!
!! System library functions and error checking
!!
!! Error routines:   simple_stop, allocchk, simple_error_check, raise_sys_error
!! System routines:  exec_cmdline, simple_sleep, simple_isenv simple_getcwd, simple_chdir simple_chmod
!! File routines:    simple_file_stat, wait_for_closure is_io is_open file_exists is_file_open
!! Memory routines:  simple_mem_usage, simple_dump_mem_usage

!! Function:  get_sys_error simple_getenv  get_lunit cpu_usage get_process_id

!! New interface: get_file_list, list_dirs, subprocess, glob_list_tofile, show_dir_content_recursive
!! New OS calls:  simple_list_dirs, simple_list_files, simple_rmdir, simple_del_files, exec_subprocess

module simple_syslib
    use simple_defs
    use iso_c_binding
#if defined(GNU)
    use, intrinsic :: iso_fortran_env, only: &
        &stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT,&
        &IOSTAT_END, IOSTAT_EOR, COMPILER_VERSION, COMPILER_OPTIONS
#elif defined(INTEL) || defined(PGI)
    use, intrinsic :: iso_fortran_env, only: &
        &stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT,&
        &IOSTAT_END, IOSTAT_EOR
#endif
#if defined(INTEL)
    use ifport, killpid=>kill, intel_ran=>ran
    use ifcore
#endif
    use simple_error
    implicit none
#if defined(PGI)
    !! include lib3f.h without kill, bessel fns, etc

    !     Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
    !
    ! NVIDIA CORPORATION and its licensors retain all intellectual property
    ! and proprietary rights in and to this software, related documentation
    ! and any modifications thereto.  Any use, reproduction, disclosure or
    ! distribution of this software and related documentation without an express
    ! license agreement from NVIDIA CORPORATION is strictly prohibited.
    !


    ! Interfaces for lib3f routines.
    ! Version 3.4
    interface

        subroutine abort()
        end subroutine abort

        integer function access(fil,mod)
            character*(*), intent(in) :: fil, mod
        end function access

        integer function chdir(path)
            character*(*), intent(in) :: path
        end function chdir

        integer function chmod(nam,mode)
            character*(*), intent(in) :: nam
            integer, intent(in) :: mode
        end function chmod

        pure character*(24) function ctime(stime)
            integer, intent(in) :: stime
        end function ctime

        pure subroutine date(str)
            character*(*), intent(out) :: str
        end subroutine date
        pure real function etime(tarray)
            real, intent(out) :: tarray(2)
        end function etime

        pure real function dtime(tarray)
            real, intent(out) :: tarray(2)
        end function dtime

        subroutine exit(s)
            integer, intent(in) :: s
        end subroutine exit

        pure subroutine fdate(str)
            character*(*), intent(out) :: str
        end subroutine fdate

        integer function fgetc(lu,ch)
            integer, intent(in) :: lu
            character*(*), intent(out) :: ch
        end function fgetc

        subroutine flush(lu)
            integer, intent(in) :: lu
        end subroutine flush

        integer function fork()
        end function fork

        integer function fputc(lu,ch)
            integer, intent(in) :: lu
            character*(*), intent(in) :: ch
        end function fputc

        subroutine free(p)
            integer, intent(in) :: p
        end subroutine free

        integer function fseek(lu,offset,from)
            integer, intent(in) :: lu, offset, from
        end function fseek

        integer function ftell(lu)
            integer, intent(in) :: lu
        end function ftell

        subroutine getarg(i,c)
            integer, intent(in) :: i
            character*(*), intent(out) :: c
        end subroutine getarg

        integer function iargc()
        end function iargc

        integer function getc(ch)
            character*(*), intent(out) :: ch
        end function getc

        integer function getcwd(dir)
            character*(*), intent(out) :: dir
        end function getcwd

        subroutine getenv(en,ev)
            character*(*), intent(in) :: en
            character*(*), intent(out) :: ev
        end subroutine getenv

        integer function getfd(lu)
            integer, intent(in) :: lu
        end function getfd

        integer function getgid()
        end function getgid

        subroutine getlog(name)
            character*(*), intent(out) :: name
        end subroutine getlog

        integer function getpid()
        end function getpid

        integer function getuid()
        end function getuid

        pure subroutine gmtime(stime,tarray)
            integer, intent(in) :: stime
            integer, intent(out) :: tarray(9)
        end subroutine gmtime

        integer function hostnm(nm)
            character*(*), intent(out) :: nm
        end function hostnm

        pure subroutine idate(date_array)
            integer, intent(out) :: date_array(3)
        end subroutine idate

        ! subroutine ioinit(cc,bz,ap,pf,vb)
        ! logical, intent(in) :: cc, bz, ap, vb
        ! character*(*), intent(in) :: pf
        ! end subroutine

        logical function isatty(lu)
            integer, intent(in) :: lu
        end function isatty

        pure subroutine itime(iarray)
            integer, intent(out) :: iarray(3)
        end subroutine itime

        ! integer function kill(pid,sig)
        ! integer, intent(in) :: pid, sig
        ! end function

        integer function link(n1,n2)
            character*(*), intent(in) :: n1, n2
        end function link

        pure integer function lnblnk(a1)
            character*(*), intent(in) :: a1
        end function lnblnk

        pure integer function loc(a)
            integer, intent(in) :: a
        end function loc

        pure subroutine ltime(stime,tarray)
            integer, intent(in) :: stime
            integer, intent(out) :: tarray(9)
        end subroutine ltime

        integer function malloc(n)
            integer, intent(in) :: n
        end function malloc

        pure integer function mclock()
        end function mclock

        integer function outstr(ch)
            character*(*), intent(in) :: ch
        end function outstr

        integer function putc(ch)
            character*(*), intent(in) :: ch
        end function putc

        integer function putenv(str)
            character*(*), intent(in) :: str
        end function putenv

        pure double precision function dflmin()
        end function dflmin

        pure double precision function dflmax()
        end function dflmax

        pure double precision function dffrac()
        end function dffrac

        pure integer function inmax()
        end function inmax

        integer function rename(from,to)
            character*(*), intent(in) :: from, to
        end function rename

        pure integer function rindex(a1,a2)
            character*(*), intent(in) :: a1, a2
        end function rindex

        integer function signal(sig,proc,flag)
            integer, intent(in) :: sig, flag
            external proc
        end function signal

        subroutine sleep(itime)
            integer, intent(in) :: itime
        end subroutine sleep

        integer function stat(nm,statb)
            character*(*), intent(in) :: nm
            integer, intent(out) :: statb(*)
        end function stat

        integer function lstat(nm,statb)
            character*(*), intent(in) :: nm
            integer, intent(out) :: statb(*)
        end function lstat

        integer function fstat(lu,statb)
            integer, intent(in) :: lu
            integer, intent(out) :: statb(*)
        end function fstat

        integer function stime(tp)
            integer, intent(in) :: tp
        end function stime

        integer function symlnk(n1,n2)
            character*(*), intent(in) :: n1, n2
        end function symlnk

        integer function system(str)
            character*(*), intent(in) :: str
        end function system

        pure integer function time()
        end function time

        pure integer function times(buf)
            integer, intent(out) :: buf(*)
        end function times

        character*(100) function ttynam(lu)
            integer, intent(in) :: lu
        end function ttynam

        integer function unlink(fil)
            character*(*), intent(in) :: fil
        end function unlink

        integer function wait(st)
            integer, intent(out) :: st
        end function wait

        subroutine pxffileno(lu,fd,err)
            integer, intent(in) :: lu
            integer, intent(out) :: fd,err
        end subroutine pxffileno

    end interface
#endif
    private :: raise_sys_error
    ! private


    !> libc interface
    interface
        ! rmdir    CONFORMING TO POSIX.1-2001, POSIX.1-2008, SVr4, 4.3BSD.
        ! On  success,  zero is returned.  On error, -1 is returned, and errno is
        ! set appropriately.
        function rmdir(dirname) bind(C, name="rmdir")
            use, intrinsic :: iso_c_binding
            integer(c_int) :: rmdir
            character(c_char),dimension(*),intent(in)  ::  dirname
        end function rmdir

        function mkdir(path,mode) bind(c,name="mkdir")
            use, intrinsic :: iso_c_binding
            integer(c_int) :: mkdir
            character(kind=c_char,len=1),dimension(*),intent(in) :: path
            integer(c_int16_t), value :: mode
        end function mkdir
    end interface
    !> SIMPLE_POSIX.c commands
    interface
         function makedir(dirname) bind(C, name="makedir")
             import
             integer(c_int) :: makedir
             character(c_char),dimension(*),intent(in)  ::  dirname
         end function makedir

         function removedir(dirname) bind(C, name="removedir")
             import
             integer(c_int) :: removedir
             character(c_char),dimension(*),intent(in)  ::  dirname
         end function

         function get_file_list(path,  ext, count) bind(c,name="get_file_list")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: get_file_list                           !> return success
            character(kind=c_char,len=1),dimension(*),intent(in)   :: path
            character(kind=c_char,len=1),dimension(*),intent(in)   :: ext
            integer(c_int), intent(inout) :: count                    !> number of elements in results
        end function get_file_list

        function get_file_list_modified(path, ext, count, flag) bind(c,name="get_file_list_modified")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: get_file_list_modified                  !> return success
            character(kind=c_char,len=1),dimension(*),intent(in)   :: path
            character(kind=c_char,len=1),dimension(3),intent(in)   :: ext
            integer(c_int), intent(inout) :: count                    !> number of elements in results
            integer(c_int), intent(in), value :: flag                 !> 1st bit reverse, 2nd bit alphanumeric sort or modified time
        end function

        function glob_file_list(av, count, flag) bind(c,name="glob_file_list")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: glob_file_list                           !> return success
            character(kind=c_char,len=1),dimension(*),intent(in)    :: av  !> glob string
            integer(c_int), intent(inout) :: count                     !> number of elements in results
            integer(c_int), intent(in)    :: flag                      !> flag 1=time-modified reverse
        end function glob_file_list

        function glob_rm_all(av, count) bind(c,name="glob_rm_all")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: glob_rm_all                              !> return success
            character(kind=c_char,len=1),dimension(*),intent(in):: av  !> glob string
            integer(c_int), intent(inout) :: count                     !> number of elements in results
        end function glob_rm_all

        function list_dirs(path,  count) bind(c,name="list_dirs")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: list_dirs                                 !> return success
            character(kind=c_char,len=1),dimension(*),intent(in)   :: path
            type(c_ptr) :: file_list_ptr
            integer(c_int), intent(inout) :: count                      !> return number of elements in results
        end function

        subroutine show_dir_content_recursive(path) bind(c,name="show_dir_content_recursive")
            use, intrinsic :: iso_c_binding
            implicit none
            character(kind=c_char,len=1),dimension(*),intent(in) :: path
        end subroutine

        function subprocess(cmd,args) bind(c,name="subprocess")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: subprocess                                  !> return PID of forked process
            character(kind=c_char,len=1),dimension(*),intent(in) :: cmd   !> executable path
            character(kind=c_char,len=1),dimension(*),intent(in) :: args  !> arguments
        end function

        function fcopy(file1, file2) bind(c,name="fcopy")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: fcopy                                       !> return success of fcopy
            character(kind=c_char,len=1),dimension(*),intent(in) :: file1
            character(kind=c_char,len=1),dimension(*),intent(in) :: file2
        end function fcopy

        subroutine free_file_list(p, n) bind(c, name='free_file_list')
            use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_char
            implicit none
            type(c_ptr), intent(in), value :: p
            integer(c_int), intent(in), value :: n
        end subroutine free_file_list

        function get_absolute_pathname(infile, outfile, outfile_length) bind(c,name="get_absolute_pathname")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: get_absolute_pathname
            character(kind=c_char,len=1),dimension(*),intent(in) :: infile
            type(c_ptr) :: outfile   !character(kind=c_char,len=1),dimension(*),intent(in) :: outfile
            integer(c_int) :: outfile_length  !> string lengths
        end function get_absolute_pathname
        function get_sysinfo(HWM, totRAM, shRAM, bufRAM, peakBuf) bind(c,name="get_sysinfo")
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_int) :: get_sysinfo
            integer(c_long), intent(inout) :: HWM, totRAM, shRAM, bufRAM, peakBuf
        end function get_sysinfo
    end interface
contains


    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline, waitflag , suppress_errors)
        character(len=*),  intent(in) :: cmdline
        logical, optional, intent(in) :: waitflag, suppress_errors
        character(len=STDLEN) :: cmsg, cmdmsg
        character(13) :: suppress_msg="2>/dev/null"
        integer ::  cstat, exec_stat
        logical :: l_doprint, wwait
        l_doprint = .false.
        wwait     = .true.
        if( present(waitflag) ) wwait = waitflag
        cmsg = trim(adjustl(cmdline))
        if( present(suppress_errors) )   cmsg = trim(adjustl(cmsg // suppress_msg))
#if defined(PGI)
        ! include 'lib3f.h'  ! PGI declares kill,wait here
        exec_stat = system(trim(adjustl(cmsg)))
        ! #elif defined(INTEL)
        !        exec_stat = system(trim(adjustl(cmdline)))
#else
        !! GNU and INTEL
        call execute_command_line( trim(adjustl(cmsg)), wait=wwait, exitstat=exec_stat, cmdstat=cstat, cmdmsg=cmdmsg)
        call raise_sys_error( cmsg, exec_stat, cstat, cmdmsg )
#endif
        if( l_doprint )then
            write(*,*) 'command: ', trim(adjustl(cmsg))
            write(*,*) 'status of execution: ', exec_stat
        endif
    end subroutine exec_cmdline


    !>  Wrapper for simple_posix's subprocess : this uses system fork & execp
    subroutine exec_subprocess( cmdline, pid)
        character(len=*),  intent(in)      :: cmdline
        integer, intent(out)               :: pid
        character(len=STDLEN)              :: tmp
        character(len=STDLEN), allocatable :: cmdexec, cmdargs
        integer                            :: pos, cmdlen
        tmp = trim(adjustl(cmdline))
        pos = INDEX(trim(tmp),' ')
        cmdlen = LEN_TRIM(tmp)
        if(pos == 0 .or. pos > cmdlen) call simple_stop("exec_subprocess called with invalid command "//trim(cmdline))
        allocate(cmdexec, source=tmp(1:pos-1)//c_null_char)
        allocate(cmdargs, source=tmp(pos:)//c_null_char)
        pid = subprocess( cmdexec , cmdargs )
        deallocate(cmdexec, cmdargs)
    end subroutine exec_subprocess

    !>  Handles error from system call
    subroutine raise_sys_error( cmd, exitstat, cmdstat, cmdmsg )
        integer,               intent(in) :: exitstat, cmdstat
        character(len=*),      intent(in) :: cmd
        character(len=STDLEN), intent(in) :: cmdmsg
        logical :: err
        err = .false.
        if( exitstat /= 0 )then
            write(*,*)'System error', exitstat,' for command: ', trim(adjustl(cmd))
            err = .true.
        endif
        if( cmdstat /= 0 )then
            call simple_error_check()
            write(*,*)'cmdstat = ',cmdstat,' command could not be executed: ', trim(adjustl(cmd))
            err = .true.
        endif
        ! if( err ) write(*,*) trim(adjustl(cmdmsg))
    end subroutine raise_sys_error

    !> isenv; return 0 if environment variable is present
    logical function simple_isenv( name )
        character(len=*), intent(in) :: name
        character(len=STDLEN)        :: varval
        integer                      :: length, status
        simple_isenv=.false.
        status=1
#if defined(PGI)
        call getenv( trim(adjustl(name)), varval)
        if(len_trim(varval) /= 0)   status=0
#else
        !! Intel and GNU F2003 included
        call get_environment_variable( trim(adjustl(name)), status=status)
#endif
        if(status==0) simple_isenv=.true.
    end function simple_isenv

    !> simple_getenv gets the environment variable string and status
    function simple_getenv( name , retval, allowfail)  result( status )
        character(len=*), intent(in)       :: name
        character(len=STDLEN), intent(out) :: retval
        logical, intent(in),optional       :: allowfail
        integer                            :: length, status
#if defined(PGI)
        call getenv( trim(adjustl(name)), retval)
#else
        !! Intel and GNU F2003 included
        call get_environment_variable( trim(name), value=retval, length=length, status=status)
        if( status == -1 ) write(*,*) 'value string too short; simple_syslib :: simple_getenv'
        if( status ==  1 )then
            write(*,*) 'environment variable: ', trim(name), ' is not defined; simple_syslib :: simple_getenv'
            retval = 'undefined'
            return
        endif
        if( status ==  2 ) write(*,*) 'environment variables not supported by system; simple_syslib :: simple_getenv'
        if( length ==  0 .or. status /= 0 )then
            retval=""
            return
        end if
#endif
    end function simple_getenv

    subroutine simple_sleep( secs )
        integer, intent(in) :: secs
        integer(dp)         :: longtime
#if defined(INTEL)
        integer             :: msecs
        msecs = 1000*INT(secs)
        call sleepqq(msecs)  !! milliseconds
#else
        call sleep(INT(secs)) !! intrinsic
#endif
    end subroutine simple_sleep


    !! SYSTEM FILE OPERATIONS
    subroutine simple_copy_file(fname1, fname2)
        character(len=*), intent(in)           :: fname1, fname2 !< input filenames
        ! character(len=1) c
        ! integer :: isrc, idst, ierr, irec, i
        ! open(newunit=isrc, file=fname1, access='direct', status='old', action='read', iostat=ierr, recl=1)
        ! open(newunit=idst, file=fname2, access='direct', status='replace', action='write', iostat=ierr)
        ! irec = 1
        ! do
        !     read(unit=isrc, rec=irec, iostat=ierr) c
        !     if (ierr.ne.0) exit
        !     write(unit=idst, rec=i) c
        !     irec = irec + 1
        ! end do
        ! close(isrc)
        ! close(idst)
        integer :: status
        status = fcopy(trim(fname1), trim(fname2))
        if (status/=0)&
            call simple_error_check(status,"simple_syslib::simple_copy_file failed "//trim(fname1)//" "//trim(fname2))
    end subroutine simple_copy_file



    !> \brief  Rename or move file
    function simple_rename( filein, fileout , overwrite) result(file_status)
        character(len=*), intent(in)  :: filein, fileout !< input filename
        logical, intent(in), optional :: overwrite      !< default true
        integer                       :: file_status
        logical                       :: force_overwrite
        force_overwrite=.true.
        if(present(overwrite)) force_overwrite=overwrite
        if( file_exists(filein) )then
            if( file_exists(fileout) .and. (.not. force_overwrite) )&
                call simple_stop("simple_rename failed to rename file,  designated output filename already exists "//trim(fileout))
            file_status = rename(filein, fileout)
            if(file_status /= 0) call simple_error_check(file_status,"simple_rename failed to rename file "//trim(filein))
        else
            call simple_stop( "simple_fileio::simple_rename, designated input filename doesn't exist "//trim(filein))
        end if

    end function simple_rename


    function simple_chmod(pathname, mode ) result( status )
        character(len=*), intent(in) :: pathname, mode
        integer :: status, imode
#if defined(INTEL)
        ! Function method
        status = chmod(pathname, mode)
#elif defined(PGI)
        ! integer :: imode
         imode=0 ! convert symbolic to octal
         if ( index(mode, 'x') /=0) imode=o'110'
         if ( index(mode, 'w') /=0) imode=o'220'
        if ( index(mode, 'r') /=0) imode=o'440'
        !status = system("chmod "//trim(adjustl(mode))//" "//trim(adjustl(pathname)))
        status = chmod(pathname, imode)
#else
        imode = INT(o'000') ! convert symbolic to octal
        if ( index(mode, 'x') /=0) imode=IOR(imode,INT(o'111'))
        if ( index(mode, 'w') /=0) imode=IOR(imode,INT(o'222'))
        if ( index(mode, 'r') /=0) imode=IOR(imode,INT(o'444'))
        status = chmod(pathname, mode) !! intrinsic GNU
#endif
        if(status/=0)&
            call simple_error_check(status,"simple_syslib::simple_chmod chmod failed "//trim(pathname))
    end function simple_chmod


    !>  Wrapper for POSIX system call stat
    subroutine simple_file_stat( filename, status, buffer, doprint )
#if defined(INTEL)
        use ifposix
#endif
        character(len=*),     intent(in)    :: filename
        integer,              intent(inout) :: status
        integer, allocatable, intent(inout) :: buffer(:)  !< POSIX stat struct
        logical, optional,    intent(in)    :: doprint
        logical :: l_print = .true., currently_opened=.false.
        integer :: funit
        character(len=STDLEN) :: io_message

#if defined(GNU)
        allocate(buffer(13), source=0)
        status = stat(trim(adjustl(filename)), buffer)

#elif defined(PGI)
        include 'lib3f.h'
        status =  stat(trim(adjustl(filename)), buffer)
        !        DebugPrint 'fileio       sys_stat PGI stato ', status
        !        DebugPrint 'fileio       sys_stat PGI size of buffer ', size(statb)
#elif defined(INTEL)

        inquire(file=trim(adjustl(filename)), opened=currently_opened, iostat=status)
        if(status /= 0) call simple_error_check(status,"simple_syslib::simple_sys_stat inquire failed "//trim(filename))
        if(.not.currently_opened) open(newunit=funit,file=trim(adjustl(filename)),status='old')
        !allocate(buffer(13), source=0)
        status = STAT (trim(adjustl(filename)) , buffer)
        if (.NOT. status) then
            call simple_error_check(status, "In simple_syslib::simple_file_stat "//trim(filename))
            print *, buffer
        end if
        if(.not.currently_opened) close(funit)
        ! integer(4) :: ierror, fsize
        ! integer(jhandle_size) :: jhandle
        ! fsize=len_trim(adjustl(filename))
        ! allocate(buffer(13), source=0)
        ! call PXFSTRUCTCREATE('stat', jhandle, ierror)
        ! if(ierror.EQ.0) then
        !    call pxfstat (trim(adjustl(filename)), fsize, jhandle, ierror)
        !    if(ierror.EQ.0) then
        !         CALL PXFINTGET (jhandle,'st_dev',buffer(1), ierror)    ! Device ID
        !         CALL PXFINTGET (jhandle,'st_ino',buffer(2), ierror)    ! Inode number
        !         CALL PXFINTGET (jhandle,'st_mode' , buffer(3), ierror) !  File mode
        !         CALL PXFINTGET (jhandle,'st_nlink' ,buffer(4), ierror) ! Number of links
        !         CALL PXFINTGET (jhandle,'st_uid' ,buffer(5), ierror)   ! Owner's uid
        !         CALL PXFINTGET (jhandle,'st_gid' ,buffer(6), ierror)   ! Owner's gid
        !         buffer(7)=0 ! ID of device containing directory entry for file (0 if not available)
        !         CALL PXFINTGET (jhandle,'st_size',buffer(8), ierror)   ! File size (bytes)
        !         CALL PXFINTGET (jhandle,'st_atime',buffer(9), ierror)  ! Last access time
        !         CALL PXFINTGET (jhandle,'st_mtime',buffer(10), ierror) ! Last modification time
        !         CALL PXFINTGET (jhandle,'st_ctime',buffer(11), ierror) ! Last file status change time
        !         buffer(12)=0 ! Preferred I/O block size (-1 if not available)
        !         buffer(13)=0 ! Number of blocks allocated (-1 if not available)
        !         call PXFSTRUCTFREE(jhandle,ierror)
        !     else
        !         print *, 'Filehandling sys_stat PXFSTAT failed, file ', trim(adjustl(filename)),' error ', ierror
        !     end if
        !     if (ierror.NE.0)then
        !         print *, 'Filehandling sys_stat PXFINTGET failed, file ', trim(adjustl(filename)),' error ', ierror
        !     end if
        ! else
        !     call PXFSTRUCTFREE(jhandle,ierror)
        !     stop 'Filehandling sys_stat  failed - cannot create structure for jhandle1'
        ! end if
        ! status=ierror
#endif
        if( present(doprint) )l_print = doprint
        if( l_print )then
            write(*,*) 'command: stat ', trim(adjustl(filename))
            write(*,*) 'status of execution: ', status
        endif
    end subroutine simple_file_stat

    logical function is_io(unit)
        integer, intent(in) :: unit
        is_io=.false.
        if (unit == stderr .or. unit == stdout .or. unit == stdin) is_io= .true.
    end function is_io

    !>  \brief  check whether a IO unit is currently opened
    logical function is_open( unit_number )
        integer, intent(in)   :: unit_number
        integer               :: io_status
        character(len=STDLEN) :: io_message
        io_status = 0
        inquire(unit=unit_number, opened=is_open,iostat=io_status,iomsg=io_message)
        if (io_status .ne. 0) then
            print *, 'is_open: IO error ', io_status, ': ', trim(adjustl(io_message))
            call simple_stop ('IO error; is_open; simple_fileio      ')
        endif
    end function is_open

    !>  \brief  check if a file exists on disk
    logical function file_exists(fname)
        character(len=*), intent(in) :: fname
        inquire(file=trim(adjustl(fname)), exist = file_exists)
    end function file_exists

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

    !>  \brief  waits for file to be closed
    subroutine wait_for_closure( fname )
        character(len=*), intent(in)  :: fname
        logical :: exists, closed
        integer :: wait_time
        wait_time = 0
        do
            if( wait_time == 60 )then
                write(*,'(A,A)')'>>> WARNING: been waiting for a minute for file: ',trim(adjustl(fname))
                wait_time = 0
                flush(stdout)
            endif
            exists = file_exists(fname)
            closed = .false.
            if( exists )closed = .not. is_file_open(fname)
            if( exists .and. closed )exit
            call simple_sleep(1)
            wait_time = wait_time + 1
        enddo
    end subroutine wait_for_closure

    !> \brief  Get current working directory
    subroutine simple_getcwd( cwd )
        character(len=*), intent(inout) :: cwd   !< output pathname
        integer :: io_status
#if defined(INTEL)
        io_status = GETCWD(cwd)
#else
        io_status = getcwd(cwd)
#endif
        if(io_status /= 0) call simple_error_check(io_status, &
            "syslib:: simple_getcwd failed to get path "//trim(cwd))
    end subroutine simple_getcwd

    !> \brief  Change working directory
    subroutine simple_chdir( newd , oldd, status )
        character(len=*), intent(in)            :: newd   !< output pathname
        character(len=*), intent(out), optional :: oldd
        integer, intent(out), optional          :: status
        character(len=STDLEN)                   :: olddir
        integer :: io_status
        logical :: dir_e, qq
        if(present(status)) status = 1
        if(present(oldd))then
            call simple_getcwd(olddir)
            oldd = trim(olddir)
        endif

        inquire(file=newd, exist=dir_e)
        if(dir_e) then
#if defined(INTEL)
            qq =  changedirqq(d)
            io_status = INT(qq)
#else
            io_status = chdir(newd)
#endif

            if(io_status /= 0) call simple_error_check(io_status, &
                "syslib:: simple_chdir failed to change path "//trim(newd))
        else
            call simple_stop("syslib:: simple_chdir directory does not exist ")
        endif
        if(present(status)) status = io_status
    end subroutine simple_chdir

    !> Make directory
    subroutine simple_mkdir_old( path )
        character(len=*), intent(in) :: path
        integer :: iostat
        logical :: dir_e
        inquire( file=trim(adjustl(path)), exist=dir_e , iostat=iostat)
        if (.not. dir_e ) then
            call exec_cmdline('mkdir -p '//trim(adjustl(path))//' | true')
        end if
    end subroutine simple_mkdir_old

    !> \brief  Make directory -- fail when ignore is false
    subroutine simple_mkdir( d , ignore, status)
        use iso_c_binding
        character(len=*), intent(in)            :: d
        logical,          intent(in), optional  :: ignore
        integer,          intent(out), optional :: status
        character(len=:), allocatable :: path
        character(len=STDLEN) :: subdirs(5)
        character(len=STDLEN) :: tmppath
        integer :: nargs, iarg, ival, io_stat
        integer :: io_status, idir, ndirs
        logical :: dir_e, ignore_here, qq
        ignore_here = .false.
        inquire(file=trim(d), exist=dir_e)
        if(.not. dir_e) then
            if(index(trim(d), '/') == 0)then
                ndirs=1
                subdirs(1) = trim(d)
            else
                call parsestr(trim(d),'/',subdirs,ndirs)
            endif
            do idir=1, ndirs
                call nsplit_str(trim(d),'/', tmppath, idir)

#if defined(INTEL)
            allocate(path, source=trim(tmppath))
            qq =  makedirqq(path)
            io_status = INT(qq)
#elif defined(GNU)
            allocate(path, source=trim(tmppath)//c_null_char)
            !o_status = makedir(path) !trim(d)//c_null_char)
            io_status= mkdir(trim(path), int(o'777',c_int16_t))
#elif defined(PGI)
            allocate(path, source=trim(tmppath))
            io_status = mkdir(path, int(o'777',c_int16_t))
#endif
            deallocate(path)
            if(.not. ignore_here)then
                if(io_status /= 0) call simple_error_check(io_status, &
                    "syslib:: simple_mkdir failed to create "//trim(d))
            endif
            enddo
        else
            if(global_verbose) print *," Directory ", d, " already exists, simple_mkdir ignoring request"
        end if
        if(present(status)) status = io_status
    end subroutine simple_mkdir

    !> \brief  Remove directory
    subroutine simple_rmdir( d , status)
        character(len=*), intent(in)            :: d
        integer,          intent(out), optional :: status
        character(len=:), allocatable           :: path
        integer                                 :: io_status
        logical                                 :: dir_e, qq
        inquire(file=d, exist=dir_e)
        if(dir_e) then
#if defined(INTEL)
            allocate(path, source=trim(d))
            qq =  deldirqq(d)
            io_status = INT(qq)
#elif defined(GNU)
            allocate(path, source=trim(d)//c_null_char)
            io_status = rmdir(path)
#elif defined(PGI)
            allocate(path, source=trim(d))
            io_status = rmdir(d)
#endif
            deallocate(path)
            if(io_status /= 0) call simple_error_check(io_status, &
                "syslib:: simple_rmdir failed to remove "//trim(d))
        else
            print *," Directory ", d, " does not exists, simple_rmdir ignoring request"
        end if
        if(present(status)) status = io_status
    end subroutine simple_rmdir

    !> file list command based on flibs-0.9
    subroutine simple_filetmp_list( pattern, list )
        character(len=*), intent(in)           :: pattern
        character(len=*), intent(out), pointer :: list(:)
        character(len=STDLEN)                  :: tmpfile
        character(len=STDLEN)                  :: cmd
        character(len=1)                       :: line
        character(len=len(pattern))            :: prefix
        integer                                :: luntmp
        integer                                :: i, ierr, count
#if _DEBUG
        print*,"simple_filetmp_list creating tmp file __simple_filelist__"
#endif

        tmpfile = '__simple_filelist__'
        cmd = trim(pattern) // ' ' // trim(redirect) // trim(tmpfile) &
            // ' ' // suppress_msg

        call exec_cmdline( cmd )
#if _DEBUG
        print*,"simple_filetmp_list reading __simple_filelist__"
#endif
        open(newunit = luntmp, file = tmpfile)
        !
        ! First count the number of files, then allocate and fill the array
        !
        count = 0
        do
            read( luntmp, '(a)', iostat = ierr ) line
            if ( ierr == 0 ) then
                count = count + 1
            else
                exit
            endif
        enddo
        rewind( luntmp )
        allocate( list(count) )
        do i = 1,count
            read( luntmp, '(a)' ) list(i)
        enddo
        close( luntmp, status = 'delete' )
#if _DEBUG
        print*,"simple_filetmp_list done "
#endif
    end subroutine simple_filetmp_list

    function simple_list_dirs(path, outfile, status) result(list)
        character(len=*), intent(in)            :: path
        character(len=*), intent(in), optional  :: outfile
        integer,          intent(out), optional :: status
        character(len=STDLEN),pointer           :: list(:)
        character(len=STDLEN)                   :: cur
        character(len=:), allocatable           :: pathhere
        integer                                 :: stat, i,num_dirs, luntmp

        allocate(pathhere, source=trim(path)//c_null_char)
        stat = list_dirs(trim(path), num_dirs)
        if(stat/=0)call simple_stop("simple_syslib::simple_list_dirs failed to process list_dirs "//trim(pathhere))
        if(present(outfile))then
            call simple_copy_file('__simple_filelist__', trim(outfile))
        endif
        open(newunit = luntmp, file = '__simple_filelist__')
        allocate( list(num_dirs) )
        do i = 1,num_dirs
            read( luntmp, '(a)' ) list(i)
        enddo
        close( luntmp, status = 'delete' )

        deallocate(pathhere)
        if(present(status)) status= stat
    end function simple_list_dirs

    !> File list -- wrapper for C functions get_file_list and glob_file_list
    !!  Path method
    !>    list = simple_list_files('./')
    !!    list = simple_list_files('./', ext='mrc', outfile='filetab.txt', tr=.true.)
    !!  Glob method:
    !!    list = simple_list_file(glob='dir/tmp*.txt')
    !! Glob function searches for all the pathnames matching pattern according to the rules
    !! used by the shell (see man 7 glob.  No tilde expansion or parameter substitution is done
    function simple_list_files(path, ext, glob, outfile, tr, status) result(list)
        character(len=*), intent(in), optional               :: path
        character(len=*), intent(in), optional               :: ext
        character(len=*), intent(in), optional               :: glob
        character(len=*), intent(in), optional               :: outfile
        logical,          intent(in), optional               :: tr  !> "ls -tr " reverse time-modified flag
        integer,          intent(out), optional              :: status
       ! type(c_ptr)                                         :: file_list_ptr
       ! character(kind=c_char,len=1), pointer, dimension(:) :: tmp
        character(len=STDLEN), allocatable :: list(:)
        character(len=STDLEN)              :: cur
        character(len=1)                   :: sep='/'
        character(len=STDLEN)              :: pathhere, thisglob, thisext !> pass these to C routines
        integer                            :: i,stat, luntmp
        integer(c_int)                     :: time_sorted_flag,num_files

        time_sorted_flag = 0
        global_debug=.true.
        pathhere= ""; thisglob=""; thisext=""
        if(present(glob)) then
            thisglob =trim(glob)//achar(0)
            if (global_debug) print *," In simple_syslib::simple_list_files glob:", trim(thisglob),":"
        else
            thisglob =trim("*")//achar(0)
        end if
        if(present(path))then
            if (global_debug) print *," In simple_syslib::simple_list_files path:", trim(path),":"
            if (len_trim(path)==0)then
                pathhere=trim("./")//achar(0)
            elseif (index(trim(path),sep)==len_trim(path)) then
                if (global_debug) &
                print *," In simple_syslib::simple_list_files path is not empty and contains a slash"
                pathhere =trim(path)//achar(0)
            else
                if (global_debug) print *," In simple_syslib::simple_list_files appending a separator to path"
                pathhere =trim(path)//trim(sep)//achar(0)
            end if
            ! ext only has effect in path method
            if(present(ext)) then
                thisext=trim(ext)//achar(0)
            else
                thisext = ""//achar(0)
            end if
            if (global_debug) print *," In simple_syslib::simple_list_files pathhere:", trim(pathhere),":"
        end if
        !! Check both path and glob methods
        if(present(path) .and. present(glob)) then
            !deallocate(thisglob)
            thisglob =trim(pathhere)//trim(thisglob)//achar(0)
            if (global_debug) print *," In simple_syslib::simple_list_files glob:", trim(thisglob),":"
        else if(present(path) .and. (.not. present(glob))) then
            ! path already set above
        else if( (.not.present(path)) .and. present(glob)) then
            ! glob already set above
        else
            call simple_stop("simple_syslib::simple_list_files Error must have either path or glob in args")
        end if
        time_sorted_flag = 0
        if(present(tr)) then
            if(tr) time_sorted_flag = 1 ! .eqv. .true.)
        end if
        !! Arguments all parsed -- Calling get_file_list or glob_file_list
        if(present(glob)) then
            !! GLOB takes precedence over get_file_list
            if(global_debug) print *, ' Calling  glob_file_list thisglob:', trim(thisglob),":"
            stat = glob_file_list(trim(thisglob), num_files, time_sorted_flag)
        else
            if(global_debug) print *, ' Calling  get_file_list pathhere:', trim(pathhere),":  ext:",trim(thisext),":"
            if(time_sorted_flag == 1) then
                stat = get_file_list_modified(trim(pathhere),trim(thisext), num_files, 3)
            else
                stat = get_file_list_modified(trim(pathhere),trim(thisext), num_files, 0)
            end if
        end if
        if(stat/=0)call simple_stop("simple_syslib::simple_list_files failed to process file list "//trim(pathhere))
        if(global_debug) print *, ' In simple_syslib::simple_list_files  num_files : ',  num_files

        if(present(outfile))then
            call simple_copy_file('__simple_filelist__', trim(outfile))
        endif
        open(newunit = luntmp, file = '__simple_filelist__')
        allocate( list(num_files) )
        do i = 1,num_files
            read( luntmp, '(a)' ) list(i)
        enddo
        close( luntmp, status = 'delete' )
        if(present(status)) status= stat
    end function simple_list_files


    !> Glob list : Emulate ls "glob" > outfile
    !!    list = glob_list_tofile(glob='dir/*tmp*.txt')
    function simple_glob_list_tofile(glob, outfile, tr) result(status)
        character(len=*), intent(in)           :: glob
        character(len=*), intent(in)           :: outfile
        logical,          intent(in), optional :: tr                   !> "ls -tr " reverse time-modified flag

        character(len=STDLEN), allocatable                  :: list(:)
        character(len=:), allocatable                       :: thisglob
        character(1) :: sep='/'
        integer      :: status
        type(c_ptr)  :: file_list_ptr
        integer      :: i,num_files, luntmp, time_sorted_flag
        time_sorted_flag = 0 !   call simple_chdir(path, cur)
        ! call simple_filetmp_list(trim(ls_command)//' --almost-all -1', list)
        if(len(glob)==0) then
            allocate(thisglob, source=trim(glob)//achar(0))
        else
            allocate(thisglob, source='*'//achar(0))
        end if
        time_sorted_flag = 0
        if(present(tr))then
            if(tr) time_sorted_flag = 1
        end if
        if(global_debug) print *, 'Calling  glob_file_list ', trim(thisglob)
        status = glob_file_list(trim(thisglob), num_files, time_sorted_flag)
        if(status/=0)call simple_stop("simple_syslib::simple_glob_list_files failed to process file list "//trim(thisglob))
        status= simple_rename('__simple_filelist__', trim(outfile))
        if(status/=0) call simple_stop("simple_syslib::simple_glob_list_files failed to copy tmpfile to "//trim(outfile))
        deallocate(thisglob)
    end function simple_glob_list_tofile


    !> \brief  is for deleting a file
    subroutine del_file( file )
        character(len=*), intent(in) :: file !< input filename
        integer :: fnr, file_status
        if( file_exists(file) )then
            !call cpStr(file,tfile)
            open(newunit=fnr,file=file,STATUS='OLD',IOSTAT=file_status)
            if( file_status == 0 )then
                close(fnr, status='delete',IOSTAT=file_status)
                if(file_status /=0) call simple_stop("simple_syslib::del_file failed to close file "//trim(file))
            end if
        endif
    end subroutine del_file

    !> generic deletion of files using POSIX glob, emulate rm -f glob
    function simple_del_files(glob, dlist, status) result(glob_elems)
        character(len=*), intent(in), optional   :: glob
        character(len=*), intent(out), optional  :: dlist(:)
        integer,          intent(out), optional  :: status
        type(c_ptr)                                :: listptr
        character(kind=c_char,len=STDLEN), pointer :: list(:)
        character(len=STDLEN)                      :: cur
        character(len=:), allocatable              :: thisglob
        integer                                    :: i, glob_elems,iostatus, luntmp

        if(present(glob))then
            thisglob=trim(glob)//c_null_char
        else
            thisglob='*'//c_null_char
        endif
        !! glob must be protected by c_null char
        iostatus =  glob_file_list(trim(thisglob), glob_elems, 0)  ! simple_posix.c
        if(status/=0) call simple_stop("simple_syslib::simple_del_files glob failed")
        open(newunit = luntmp, file = '__simple_filelist__')
        allocate( list(glob_elems) )
        do i = 1,glob_elems
            read( luntmp, '(a)' ) list(i)
        enddo
        close( luntmp, status = 'delete' )

        if ( glob_elems > 0) then
            do i=1,glob_elems
               ! if(global_debug) print*, 'Deleting ', list(i)
                call del_file(list(i))
            enddo
            do i=1,glob_elems
                if(file_exists(list(i))) call simple_stop(" simple_del_files failed to delete "//trim(list(i)))
            enddo
        else
            print *,"simple_syslib::simple_del_files no files matching ", trim(thisglob)
        endif
       ! call simple_chdir(cur)
        if(present(dlist)) dlist = list
        if(present(status))status=iostatus
        deallocate(list)
        deallocate(thisglob)
    end function simple_del_files

    !> forced deletion of dirs and files using POSIX glob -- rm -rf glob
    !! num_deleted = simple_rm_force("tmpdir*/tmp*.ext", dlist=list_of_deleted_elements, status=stat)
    function simple_rm_force(glob, dlist, status) result(glob_elems)
        character(len=*), intent(in), optional     :: glob
        character(len=*), intent(out), optional    :: dlist(:)
        integer,          intent(out), optional    :: status
        character(kind=c_char,len=STDLEN), pointer :: list(:)
        character(len=STDLEN)                      :: cur
        character(len=:), allocatable              :: thisglob
        type(c_ptr)                                :: listptr
        integer                                    :: i, glob_elems,iostatus, luntmp

        if(present(glob))then
            allocate(thisglob, source=trim(glob)//c_null_char)
        else
            allocate(thisglob, source='*'//c_null_char)
        endif
        !! glob must be protected by c_null char
        iostatus =  glob_rm_all(trim(thisglob), glob_elems)  ! simple_posix.c
        if(status/=0) call simple_stop("simple_syslib::simple_del_files glob failed")
        open(newunit = luntmp, file = '__simple_filelist__')
        allocate( list(glob_elems) )
        do i = 1,glob_elems
            read( luntmp, '(a)' ) list(i)
        enddo
        close( luntmp, status = 'delete' )

        if ( glob_elems > 0) then
            do i=1, glob_elems
                if(global_debug) print*, 'Checking deleted file or dir ', list(i)
                if(file_exists(list(i))) then
                    print *, " failed to delete "//trim(list(i))
                    call simple_stop(" simple_syslib::simple_rm_force ")
                end if
            enddo
        else
            print *,"simple_syslib::simple_del_files no files matching ", trim(thisglob)
        endif
       ! call simple_chdir(cur)
        if(present(dlist)) dlist = list
        if(present(status))status=iostatus
        deallocate(list)
        deallocate(thisglob)
    end function simple_rm_force


    ! !>  \brief  get logical unit of file
    ! integer function get_lunit( fname )
    !     character(len=*), intent(in) :: fname

    !     integer               :: io_status
    !     character(len=STDLEN) :: io_message
    !     io_status = 0
    !     inquire(file=trim(adjustl(fname)),unit=get_lunit,iostat=io_status,iomsg=io_message)
    !     if (io_status .ne. 0) then
    !         print *, 'is_open: IO error ', io_status, ': ', trim(adjustl(io_message))
    !         call simple_stop ('IO error; is_open; simple_fileio      ')
    !     endif
    ! end function get_lunit

    !> simple_timestamp prints time stamp (based on John Burkardt's website code)
    subroutine simple_timestamp ( )
        character(len= 8) :: ampm
        integer (kind=sp) :: d
        integer (kind=sp) :: h
        integer (kind=sp) :: m
        integer (kind=sp) :: mm
        character (len=9 ), parameter, dimension(12) :: month = (/ &
            'January  ', 'February ', 'March    ', 'April    ', &
            'May      ', 'June     ', 'July     ', 'August   ', &
            'September', 'October  ', 'November ', 'December ' /)
        integer    :: n, s, y, values(8)

        call date_and_time(values=values)
        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)
        if ( h < 12 ) then
            ampm = 'AM'
        else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
                ampm = 'Noon'
            else
                ampm = 'PM'
            end if
        else
            h = h - 12
            if ( h < 12 ) then
                ampm = 'PM'
            else if ( h == 12 ) then
                if ( n == 0 .and. s == 0 ) then
                    ampm = 'Midnight'
                else
                    ampm = 'AM'
                end if
            end if
        end if
        write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
            d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
    end subroutine simple_timestamp

    function cpu_usage ()
        real    :: cpu_usage
        integer :: ios, i
        integer :: unit,oldidle, oldsum, sumtimes = 0
        real    :: percent = 0.
        character(len = 4) lineID ! 'cpu '
        integer, dimension(9) :: times = 0
        cpu_usage=0.0
#ifdef LINUX
        write(*, *) 'CPU Usage'
        open(newunit=unit, file = '/proc/stat', status = 'old', action = 'read', iostat = ios)
        if (ios /= 0) then
            print *, 'Error opening /proc/stat'
            stop
        else
            read(unit, fmt = *, iostat = ios) lineID, (times(i), i = 1, 9)
            if (ios /= 0) then
                print *, 'Error reading /proc/stat'
                stop
            end if
            close(unit, iostat = ios)
            if (ios /= 0) then
                print *, 'Error closing /proc/stat'
                stop
            end if
            if (lineID /= 'cpu ') then
                print *, 'Error reading /proc/stat'
                stop
            end if
            sumtimes = sum(times)
            percent = (1. - real((times(4) - oldidle)) / real((sumtimes - oldsum))) * 100.
            write(*, fmt = '(F6.2,A2)') percent, '%'
            oldidle = times(4)
            oldsum = sumtimes

        end if
        cpu_usage=percent
#else
        write(*, *) 'CPU Usage not available'
#endif
    end function cpu_usage

    !! SYSTEM INFO ROUTINES

    integer function get_process_id( )
        get_process_id = getpid()
    end function get_process_id

    subroutine print_compiler_info(file_unit)
        use simple_strings, only: int2str
#ifdef INTEL
        character*56  :: str
#endif

#ifdef GNU
        character(*), parameter :: compilation_cmd = compiler_options()
        character(*), parameter :: compiler_ver = compiler_version()
#endif

        integer, intent (in), optional :: file_unit
        integer  :: file_unit_op
        integer  :: status


        if (present(file_unit)) then
            file_unit_op = file_unit
        else
            file_unit_op = stdout
        end if
#if defined(GNU)
        write( file_unit_op, '(A,A,A,A)' ) &
            ' This file was compiled by ', trim(adjustl(compiler_ver)), &
            ' using the options ', trim(adjustl(compilation_cmd))
#endif
#ifdef INTEL
        status = for_ifcore_version( str )
        write( file_unit_op, '(A,A)' ) &
            ' Intel IFCORE version ', trim(adjustl(str))
        status = for_ifport_version( str )
        write( file_unit_op, '(A,A)' ) &
            ' Intel IFPORT version ', trim(adjustl(str))
#endif

    end subroutine print_compiler_info


#if defined(INTEL)

    ! pure subroutine rt_init (arg)
    !     real, intent(in) :: arg
    !     io_status = for_rtl_finish_ ( )

    ! end subroutine rt_init
    ! pure subroutine rt_shutdown (arg)
    !     real, intent(in) :: arg
    !     io_status = for_rtl_finish_ ( )

    ! end subroutine rt_shutdown

    ! IFCORE checks
    function fastmem_policy() result (policy)
        integer(4) :: old_policy, policy
        print *,"Print the current Results of initial for_set_fastmem_policy(FOR_K_FASTMEM_INFO):"
        old_policy = for_set_fastmem_policy(FOR_K_FASTMEM_INFO)
        select case (old_policy)
        case (FOR_K_FASTMEM_NORETRY)
            print *,"    Issue a Severe error if FASTMEM is not available."
        case (FOR_K_FASTMEM_RETRY_WARN)
            print *,"    Issue a Warning if FASTMEM is not available, then use the default memory allocator."
        case (FOR_K_FASTMEM_RETRY)
            print *,"    If FASTMEM is not available, then silently use the default memory allocator."
        end select

        print *,"Set the FASTMEM policy to RETRY_WARN:"
        policy = for_set_fastmem_policy(FOR_K_FASTMEM_RETRY_WARN)
    end function fastmem_policy

    ! integer function hbw_availability
    !     use ifcore
    !     !        integer(4) :: hbw_availability
    !     hbw_availability = for_get_hbw_availability()
    !     print *,"Results of for_get_hbw_availability():"
    !     select case (hbw_availability)
    !     case (FOR_K_HBW_AVAILABLE)
    !         print *,"    HBM is available."
    !     case (FOR_K_HBW_NO_ROUTINES)
    !         print *,"    The libmemkind routines needed for hbw memory allocation are not linked-in."
    !     case (FOR_K_HBW_NOT_AVAILABLE)
    !         print *,"    HBM is NOT available on this node."
    !     end select
    ! end function hbw_availability

    ! function get_hbw_size (partition, total, free) result(istatus)
    !     use ifcore
    !     use, intrinsic :: iso_c_binding
    !     integer(kind=4),         intent(in)  :: partition
    !     integer(kind=int_ptr_kind()), intent(inout) :: total,free
    !     integer(kind=C_INT) :: part, istatus
    !     integer(kind=C_SIZE_T) :: tot, freemem
    !     istatus = for_get_hbw_size(part,tot,freemem )
    !     call simple_error_check(istatus, "syslib::get_hbw_size ")
    !     total = tot
    !     free = freemem
    ! end function get_hbw_size

#endif

    ! Suggestion from https://stackoverflow.com/a/30241280
    subroutine simple_mem_usage(valueRSS,valuePeak,valueSize,valueHWM)
        implicit none
        integer(kind=8), intent(out) :: valueRSS
        integer(kind=8), intent(out), optional :: valuePeak
        integer(kind=8), intent(out), optional :: valueSize
        integer(kind=8), intent(out), optional :: valueHWM

        character(len=200) :: filename=' '
        character(len=80)  :: line
        character(len=8)   :: pid_char=' '
        integer            :: pid,unit
        logical            :: ifxst, debug
        debug=.false.
        valueRSS=-1    ! return negative number if not found

        !--- get process ID

        pid=getpid()
        write(pid_char,'(I8)') pid
        filename='/proc/'//trim(adjustl(pid_char))//'/status'
        if(debug) print *,'simple_mem_usage:debug:  Fetching ', trim(filename)
        !--- read system file

        inquire (file=trim(filename),exist=ifxst)
        if (.not.ifxst) then
            write (*,*) 'system file does not exist'
            return
        endif

        open(newunit=unit, file=filename, action='read')
        ! the order of the following do loops is dependent on cat /proc/PID/status listing
        if(present(valuePeak))then
            do
                read (unit,'(a)',end=110) line
                if (line(1:7).eq.'VmPeak:') then
                    read (line(8:),*) valuePeak
                    if(debug) print *,'simple_mem_usage:debug:  Peak ', valuePeak
                    exit
                endif
            enddo
110         continue
        endif
        if(present(valueSize))then
            do
                read (unit,'(a)',end=120) line
                if (line(1:7).eq.'VmSize:') then
                    read (line(8:),*) valueSize
                    if(debug) print *,'simple_mem_usage:debug:  VM Size ', valueSize
                    exit
                endif
            enddo
120         continue
        endif
        if(present(valueHWM))then
            do
                read (unit,'(a)',end=130) line
                if (line(1:6).eq.'VmHWM:') then
                    read (line(7:),*) valueHWM
                    if(debug) print *,'simple_mem_usage:debug:  peak RAM ', valueHWM
                    exit
                endif
            enddo
130         continue
        endif
        do
            read (unit,'(a)',end=140) line
            if (line(1:6).eq.'VmRSS:') then
                read (line(7:),*) valueRSS
                if(debug) print *,'simple_mem_usage:debug: RSS ', valueRSS
                exit
            endif
        enddo
140     continue
        close(unit)
        return
    end subroutine simple_mem_usage

    subroutine simple_dump_mem_usage(dump_file)
        character(len=*), intent(inout), optional :: dump_file
        character(len=200)    :: filename=' '
        character(len=80)     :: line
        character(len=8)      :: pid_char=' '
        character(len=STDLEN) :: command
        integer               :: pid,unit
        logical               :: ifxst

        !--- get process ID --  make thread safe

        !  omp critical dump_mem
        pid=getpid()
        write(pid_char,'(I8)') pid
        filename='/proc/'//trim(adjustl(pid_char))//'/status'
        command = 'grep -E "^(VmPeak|VmSize|VmHWM|VmRSS):"<'//trim(filename)//'|awk {a[NR-1]=$2}END{print a[0],a[1],a[2],a[3]} '
!!         | awk {a[NR-1]=$2} END{print a[0],a[1],a[2],a[3]}
        if(present(dump_file)) command = trim(command)//' >> '//trim(dump_file)
        call exec_cmdline(trim(command))
        !  omp end critical
    end subroutine simple_dump_mem_usage


    function simple_full_path (infile) result(canonical_name)
        character(kind=c_char, len=*), intent(in)  :: infile
        character(len=STDLEN) :: canonical_name
        integer(1), dimension(:), pointer          :: iptr
        type(c_ptr) :: cstring
        character, pointer                            :: fstring(:)
        integer                                       :: slen, i


        integer :: lengthin, status, lengthout
        character(len=LINE_MAX_LEN), target :: fstr
        character(len=STDLEN) :: infile_c
        lengthin = len_trim(infile)
         print *, " address fstr ", loc(fstr)
        cstring = c_loc(fstr)
        infile_c = trim(infile)//achar(0)
        print *, " address cstring ", loc(cstring)
        status = get_absolute_pathname(infile_c, cstring, lengthout )
        print *, " address cstring ", loc(cstring)
        print *, " address fstr ", loc(fstr)
     !   print *," fstr ", trim(adjustl(fstr))
       if( lengthout > 1)then
           call c_f_pointer(cstring, fstring, [STDLEN] )
           write(*,'(a)') fstring
!           allocate(canonical_name)
           !           canonical_name = transfer(strptr,c_char)
           canonical_name= trim(adjustl(fstr))
           ! do i= 1,lengthout
           !     canonical_name(i) = fstring(i)
           ! enddo
        else
            canonical_name =infile
        end if

        print *, 'input file/path name: ', infile
        print *, 'absolute file/path name: ', canonical_name
    end function simple_full_path

    subroutine set_absolute_str_array(n, cstring) bind(C)
        use iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_loc, c_null_char
        implicit none
        integer(kind=c_int),               intent(in) :: n
        type(c_ptr), dimension(n), target, intent(in) :: cstring
        character, pointer                            :: fstring(:)
        integer                                       :: slen, i

        do i = 1, n
            call c_f_pointer(cstring(i), fstring, [4])
            write(*,*) fstring
        end do

    end subroutine set_absolute_str_array

end module simple_syslib
