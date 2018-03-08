!!
!! System library functions and error checking
!!

module simple_syslib
    use simple_defs
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
        end function

        pure real function dtime(tarray)
        real, intent(out) :: tarray(2)
        end function

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

        subroutine gerror(str)
            character*(*), intent(out) :: str
        end subroutine gerror

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

        integer function ierrno()
        end function ierrno

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

        subroutine perror(str)
            character*(*), intent(in) :: str
        end subroutine perror

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

contains

    !! ERROR Routines

    subroutine simple_stop (msg,f,l)
        character(len=*),      intent(in) :: msg
        character(len=*),      intent(in), optional :: f !< filename of caller
        integer,               intent(in), optional :: l !< line number from calling file
        integer   :: last_err
        if(present(f).and.present(l))&
            write(stderr,'("Stopping in file ",/,A,/," at line ",I0)') f,l
        write(stderr,'("Message: ",/,a,/)') trim(msg)
        last_err = get_sys_error()
        stop
    end subroutine simple_stop

#ifdef PGI
    subroutine simple_cuda_stop (msg,f,l)
        use cudafor
        character(len=*),      intent(in) :: msg
        character(len=*),      intent(in), optional :: f !< filename of caller
        integer,               intent(in), optional :: l !< line number from calling file
        integer   :: last_err
        last_err = cudaGetLastError()
        if ( last_err .ne. 0)then
           write(*,"('CUDA Error ',a,', File: ',a,', Line: ',i0)")msg, f, l
           write(*,"('Last CUDA Error :',i0)") last_err
           write(*,"(4x,a)")cudaGetErrorString( last_err )
           call simple_stop("In simple_syslib::simple_cuda_stop ",f,l)
       end if

    end subroutine simple_cuda_stop
#endif
    function get_sys_error () result(err)
        integer :: err
        CHARACTER(len=100) :: msg
        err = int( IERRNO(), kind=4 ) !!  EXTERNAL;  no implicit type in INTEL

        if( err < 0)then !! PGI likes to use negative error numbers
            !#ifdef PGI
            !            msg = gerror()
            !#else
            call gerror(msg) !! EXTERNAL;
            !#endif
            write(stderr,'("SIMPLE_SYSLIB::SYSERROR NEG ",I0)') err
            write(stderr,*) trim(msg)
        else if (err /= 0) then
#ifdef GNU
            write(msg,'("Last detected error (SIMPLE_SYSLIB::SYSERROR) ",I0,":")') err
            call perror(trim(adjustl(msg)))

#elif defined(INTEL)
            ! err= getlasterror()
            ! using IERRNO result from above
            select case (err)
            case (1) ! EPERM
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 1 Insufficient permission for operation"
            case (2) ! ENOENT
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 2 No such file or directory"
            case (3) ! ESRCH
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 3 No such process"
            case (4) ! EIO
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 4 INTR error"
            case (5) ! EIO
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 5 I/O error"
            case (6) ! ENXIO
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 5 NXIO error"
            case (7) ! E2BIG)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 7 Argument list too long"
            case (8) ! ENOEXEC)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 8 File is not executable"
            case (9)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 9 Bad file"
            case (10)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 8 Bad child"
            case (11)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 9 AGAIN error"
            case ( ENOMEM)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 12 Not enough resources"
            case ( EACCES)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 13  Access denied; the file's permission setting does not allow the specified access. Permission denied"
            case (14)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 14 FAULT"
            case (15)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 15 NOTBLK       "
            case (16) ! ENOEXEC)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 16 BUSY         "
            case (17) ! ENOEXEC)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 17 EXIST        "
            case ( EXDEV)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 18 Cross-device link"
            case (19) ! ENOEXEC)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 19 ERR$NODEV"
            case ( ENOTDIR)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 20 Not a directory"
            case (21)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 21 Is a directory"
            case (  EINVAL)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 22 Invalid argument"
            case (23)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 23 NFILE       "
            case (24)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 24 MFILE       "
            case (25)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 25 NOTTY       "
            case (26)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 26 TXTBSY      "
            case (27)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 27 FBIG        "
            case (28)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 28 NOSPC       "
            case (29)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 29 SPIPE       "
            case (30)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 30 ROFS        "
            case (31)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 31 MLINK       "
            case (32)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 32 PIPE        "
            case (33)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 33 DOM         "
            case (34)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 34 RANGE       "
            case (35)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 35 UCLEAN      "
            case (36)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 36 DEADLOCK    "
            case (38)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 38 NAMETOOLONG "
            case (39)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 39 NOLCK       "
            case (40)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 40 NOSYS       "
            case (41)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 41 NOTEMPTY    "
            case (42)
                write(stderr,'(A)')"SIMPLE_SYSLIB::SYSERROR 42 ILSEQ       "
            end select
            write(msg,'("SIMPLE SYSERROR ",I0)') err
            call perror(msg)
#else
            !! PGI
            ! msg = gerror()
            ! write(stderr,'("SIMPLE_SYSLIB::SYSERROR ",I0)',advance='no') err
            ! write(stderr,*) trim(msg)
            write(msg,'("Last detected error (SIMPLE_SYSLIB::SYSERROR) ",I0,":")') err
            call perror(trim(adjustl(msg)))
#endif
        end if
    end function get_sys_error

    !> \brief  is for checking allocation
    subroutine allocchk( message, alloc_err, file, line, iomsg )
        character(len=*), intent(in)           :: message
        integer,          intent(in), optional :: alloc_err
        character(len=*), intent(in), optional :: file !< filename of caller
        integer,          intent(in), optional :: line !< line number from calling file
        character(len=*), intent(in), optional :: iomsg !< IO message
        integer                                :: syserr, alloc_status
        alloc_status=alloc_stat    !! global variable from simple_defs
        if(present(alloc_err))alloc_status=alloc_err
        if (alloc_status/=0)then
            write(stderr,'(a)') 'ERROR: Allocation failure!'
            call simple_error_check(alloc_status)
            if(present(iomsg))&
                write(stderr,'("IO Message ",A)') trim(adjustl(iomsg))
            if(present(file).and.present(line))&
                write(stderr,'("Stopping in file ",/,A,/," at line ",I0)') file,line
            call simple_stop(message)
        endif
    end subroutine allocchk

    subroutine simple_error_check(io_stat, msg)
        integer,          intent(in), optional :: io_stat
        character(len=*), intent(in), optional :: msg
        integer :: io_stat_this, last_sys_error
        character(len=STDLEN ) :: newmsg
        if (present(io_stat))then
            io_stat_this = io_stat
        else
            io_stat_this = get_sys_error()
        end if

#ifdef PGI
        if(io_stat_this < 0)then
            if (IS_IOSTAT_END(io_stat_this))then
                write(stderr,'(a)')"fclose EOF reached (PGI version)"
                io_stat_this=0
            else if (IS_IOSTAT_EOR(io_stat_this)) then
                write(stderr,'(a)')"fclose End-of-record reached (PGI version)"
                io_stat_this=0
            end if
        end if
#else
        !! Intel and GNU
        if (io_stat_this==IOSTAT_END)  write(*,'(a,1x,I0 )') 'ERRCHECK: EOF reached, end-of-file reached IOS# ', io_stat_this
        if (io_stat_this==IOSTAT_EOR)  write(*,'(a,1x,I0 )') 'ERRCHECK: EOR reached, read was short, IOS# ', io_stat_this

#endif
        if( io_stat_this /= 0 ) then
            write(stderr,'(a,1x,I0 )') 'ERROR: File I/O failure, IOS# ', io_stat_this
            if(present(msg)) write(stderr,'(a)') trim(adjustl(msg))
            !! do not stop yet -- let the fopen/fclose call to finish
            !! stop
        endif

    end subroutine simple_error_check



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
    function simple_isenv( name )  result( status )
        character(len=*), intent(in)  :: name
        character(len=STDLEN) :: varval
        integer :: length, status
        status=1
#if defined(PGI)
        call getenv( trim(adjustl(name)), varval)
        if(len_trim(varval) /= 0)   status=0
#else
        !! Intel and GNU F2003 included
        call get_environment_variable( trim(adjustl(name)), status=status)
#endif
    end function simple_isenv

    !> simple_getenv gets the environment variable string and status
    function simple_getenv( name , retval, allowfail)  result( status )
        character(len=*), intent(in)  :: name
        character(len=STDLEN), intent(out)    :: retval
        logical, intent(in),optional :: allowfail
        integer :: length, status
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
#if defined(INTEL)
        integer  :: msecs
        msecs = 1000*secs
        call sleepqq(msecs)  !! milliseconds
#else
        call sleep(secs) !! intrinsic
#endif
    end subroutine simple_sleep


    !! SYSTEM FILE OPERATIONS

    function simple_chmod(pathname, mode ) result( status )
        character(len=*), intent(in) :: pathname, mode
        integer :: status, imode
#if defined(INTEL)
        ! Function method
        status = chmod(pathname, mode)
#elif defined(PGI)
        ! integer :: imode
        ! imode=0 ! convert symbolic to octal
        ! if ( index(mode, 'x') /=0) imode=o'110'
        ! if ( index(mode, 'w') /=0) imode=o'220'
        ! if ( index(mode, 'r') /=0) imode=o'440'
        status = system("chmod "//trim(adjustl(mode))//" "//trim(adjustl(pathname)))
#else
        imode = INT(o'000') ! convert symbolic to octal
        if ( index(mode, 'x') /=0) imode=IOR(imode,INT(o'111'))
        if ( index(mode, 'w') /=0) imode=IOR(imode,INT(o'222'))
        if ( index(mode, 'r') /=0) imode=IOR(imode,INT(o'444'))
        status = chmod(pathname, mode) !! intrinsic
#endif
        if(status/=0)&
            call simple_error_check(status,"simple_syslib::simple_chmod chmod failed "//trim(pathname))
    end function simple_chmod


    !>  Wrapper for POSIX system call stat
    subroutine simple_file_stat( filename, status, buffer, doprint )
#if defined(INTEL)
        use ifposix
#endif
        character(len=*),     intent(in) :: filename
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
#include "simple_local_flags.inc"
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
        io_status = getcwd(cwd)
        if(io_status /= 0) call simple_error_check(io_status, &
             "syslib:: simple_getcwd failed to get path "//trim(cwd))
    end subroutine simple_getcwd

    !> \brief  Change working directory
    subroutine simple_chdir( newd )
        character(len=*), intent(in) :: newd   !< output pathname
        integer :: io_status
        io_status = chdir(newd)
        if(io_status /= 0) call simple_error_check(io_status, &
            "syslib:: simple_getcwd failed to change path "//trim(newd))
    end subroutine simple_chdir

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
        real :: cpu_usage

        integer :: ios, i
        integer :: unit,oldidle, oldsum, sumtimes = 0
        real :: percent = 0.
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
    CHARACTER(*), PARAMETER :: compilation_cmd = COMPILER_OPTIONS()
    CHARACTER(*), PARAMETER :: compiler_ver = COMPILER_VERSION()
#endif

    integer , intent (in), optional :: file_unit
    integer  :: file_unit_op
    integer       :: status


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
        integer, intent(out) :: valueRSS
        integer, intent(out), optional :: valuePeak
        integer, intent(out), optional :: valueSize
        integer, intent(out), optional :: valueHWM

        character(len=200):: filename=' '
        character(len=80) :: line
        character(len=8)  :: pid_char=' '
        integer :: pid,unit
        logical :: ifxst

        valueRSS=-1    ! return negative number if not found

        !--- get process ID

        pid=getpid()
        write(pid_char,'(I8)') pid
        filename='/proc/'//trim(adjustl(pid_char))//'/status'

        !--- read system file

        inquire (file=filename,exist=ifxst)
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
                endif
            enddo
110         continue
        endif
        if(present(valueSize))then
            do
                read (unit,'(a)',end=120) line
                if (line(1:7).eq.'VmSize:') then
                    read (line(8:),*) valueSize
                endif
            enddo
120         continue
        endif
        if(present(valueHWM))then
            do
                read (unit,'(a)',end=130) line
                if (line(1:6).eq.'VmHWM:') then
                    read (line(7:),*) valueHWM
                endif
            enddo
130         continue
        endif
        do
            read (unit,'(a)',end=140) line
            if (line(1:6).eq.'VmRSS:') then
                read (line(7:),*) valueRSS
                exit
            endif
        enddo
140     continue
        close(unit)

        return
    end subroutine simple_mem_usage

    subroutine simple_dump_mem_usage(dump_file)
        character(len=*), intent(inout), optional :: dump_file
        character(len=200):: filename=' '
        character(len=80) :: line
        character(len=8)  :: pid_char=' '
        character(len=STDLEN) :: command
        integer :: pid,unit
        logical :: ifxst

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

end module simple_syslib
