!!
!! System library functions and error checking
!!
#include "simple_lib.f08"
module simple_syslib
    use simple_defs
    use simple_strings, only: cpStr
    use, intrinsic :: iso_fortran_env, only: &
        &stderr=>ERROR_UNIT, stdout=>OUTPUT_UNIT,&
        &IOSTAT_END, IOSTAT_EOR
    implicit none

#if defined(INTEL)
    ! interface
    !     integer function ierrno()
    !     end function ierrno
    !     pure character*(24) function ctime(stime)
    !         integer, intent(in) :: stime
    !     end function ctime
    !     pure integer function time()
    !     end function time
    ! end interface
#endif
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

#ifdef PGI

    ! include 'lib3f.h'
    ! public :: exec_cmdline, simple_getenv, ls_mrcfiletab, ls_filetab, &
    !    sys_del_files, sys_get_last_fname, simple_sleep, &
    !     sys_merge_docs

#endif

contains

    subroutine simple_stop (msg,f,l)
        character(len=*),      intent(in) :: msg
        character(len=*),      intent(in), optional :: f !< filename of caller
        integer,               intent(in), optional :: l !< line number from calling file
        if(present(f).and.present(l))&
            write(stderr,'("Stopping in file ",/,A,/," at line ",I0)') f,l
        write(stderr,'("Message: ",/,a,/)') trim(msg)
        write(stderr,'(I0)') get_sys_error()
        stop
    end subroutine simple_stop

    function get_sys_error () result(err)
#ifdef INTEL
        use ifport
        use ifcore
#endif
        integer :: err
        CHARACTER(len=100) :: msg
        err = INT( IERRNO(), kind=4 ) !! intrinsic,  no implicit type in INTEL

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
            write(msg,'("SIMPLE_SYSLIB::SYSERROR ",I0)') err
            call perror(msg)

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

            write(msg,'("SIMPLE_SYSLIB::SYSERROR ",I0)') err
            call perror(msg)
#endif
        end if
    end function get_sys_error



    !> \brief  is for checking allocation
    subroutine alloc_errchk( message, alloc_status, file,line, iomsg )
        character(len=*), intent(in)           :: message
        integer,          intent(in)           :: alloc_status
        character(len=*), intent(in), optional :: file !< filename of caller
        integer,          intent(in), optional :: line !< line number from calling file
        character(len=*), intent(in), optional :: iomsg !< IO message
        integer                                :: syserr
        
        if (alloc_status/=0)then
            write(stderr,'(a)') 'ERROR: Allocation failure!'
            call simple_error_check(alloc_status)
            if(present(iomsg))&
                write(stderr,'("IO Message ",A)') trim(adjustl(iomsg))
            if(present(file).and.present(line))&
                write(stderr,'("Stopping in file ",/,A,/," at line ",I0)') file,line
            call simple_stop(message)
        endif
    end subroutine alloc_errchk



    subroutine simple_error_check (iostat, msg)
#ifdef INTEL
        use ifport
        use ifcore
#endif
        integer,          intent(in), optional :: iostat
        character(len=*), intent(in), optional :: msg
        integer :: io_stat
        character(len=STDLEN ) :: newmsg
        if (present(iostat))then
            io_stat = iostat
        else
            io_stat = get_sys_error()
        end if

#ifdef PGI
        if(io_stat < 0)then
            if (IS_IOSTAT_END(io_stat))then
                write(stderr,'(a)')"fclose EOF reached (PGI version)"
                iostat=0
            else if (IS_IOSTAT_EOR(io_stat)) then
                write(stderr,'(a)')"fclose End-of-record reached (PGI version)"
                iostat=0
            end if
        end if
#else
        !! Intel and GNU
        if (io_stat==IOSTAT_END)  write(*,'(a,1x,I0 )') 'ERRCHECK: EOF reached, end-of-file reached IOS# ', io_stat
        if (io_stat==IOSTAT_EOR)  write(*,'(a,1x,I0 )') 'ERRCHECK: EOR reached, read was short, IOS# ', io_stat

#endif
        if( io_stat /= 0 ) then
            write(stderr,'(a,1x,I0 )') 'ERROR: File I/O failure, IOS# ', io_stat
            if(present(msg)) write(stderr,'(a)') trim(adjustl(msg))
            !! do not stop yet -- let the fopen/fclose call to finish
            !! stop
        endif

    end subroutine simple_error_check


    subroutine print_compiler_info(file_unit)
      use simple_strings, only: int2str
#ifdef GNU
        use, intrinsic :: iso_fortran_env, only: compiler_version, &
            &compiler_options
#endif
#ifdef INTEL
        use ifport
        use ifcore
        integer       :: res
        character*56  :: str
#endif
        integer , intent (in), optional :: file_unit
        integer  :: file_unit_op

        if (present(file_unit)) then
            file_unit_op = file_unit
        else
            file_unit_op = stdout
        end if
#ifdef GNU
        write( file_unit_op, '(/4a/)' ) &
            ' This file was compiled by ', COMPILER_VERSION(), &
            ' using the options ', COMPILER_OPTIONS()
#endif
#ifdef INTEL
        res = for_ifcore_version( str )
        if (res == 0 ) call simple_stop("print_compiler_info for_ifcore_version returned "//int2str(res))
        write( file_unit_op, '(/2a/)' ) " Intel IFCORE version ", trim(adjustl(str))
        res = for_ifport_version( str )
        if (res == 0 ) call simple_stop("print_compiler_info for_ifport_version returned "//int2str(res))
        write( file_unit_op, '(/2a/)' ) ' Intel IFPORT version ', trim(adjustl(str))
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
        use ifcore
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


    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline, waitflag )
#if defined(INTEL)
        use ifport
#endif
        character(len=*),  intent(in) :: cmdline
        logical, optional, intent(in) :: waitflag
        character(len=STDLEN) :: cmsg
        integer ::  cstat, exec_stat
        logical :: l_doprint = .true., wwait = .true.
        wwait = .true.
        if( present(waitflag) ) wwait = waitflag

#if defined(PGI)
        ! include 'lib3f.h'  ! PGI declares kill,wait here
        exec_stat = system(trim(adjustl(cmdline)))

        ! #elif defined(INTEL)
        !        exec_stat = system(trim(adjustl(cmdline)))

#else
        !! GNU
        call execute_command_line( trim(adjustl(cmdline)), wait=wwait, exitstat=exec_stat, cmdstat=cstat, cmdmsg=cmsg)
        call raise_sys_error( cmdline, exec_stat, cstat, cmsg )
#endif

        if( l_doprint )then
            write(*,*) 'command: ', trim(adjustl(cmdline))
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

    character(len=STDLEN) function simple_getenv( name ) ! result( varval )
        character(len=*), intent(in)  :: name
        character(len=STDLEN)         :: value
        character(len=:), allocatable :: varval
        integer :: length, status

#if defined(PGI)
        call getenv( trim(name), value)
#else
        !! Intel and GNU F2003 included
        call get_environment_variable( trim(name), value=value, length=length, status=status)
        if( status == -1 ) write(*,*) 'value string too short; simple_syslib :: simple_getenv'
        if( status ==  1 ) write(*,*) 'environment variable: ', trim(name), ' is not defined; simple_syslib :: simple_getenv'
        if( status ==  2 ) write(*,*) 'environment variables not supported by system; simple_syslib :: simple_getenv'
        if( length ==  0 .or. status /= 0 ) return
#endif
        
        write(simple_getenv,'(A)') value
        !        call alloc_errchk("In syslib::simple_getenv ", alloc_stat)
    end function simple_getenv

    subroutine simple_sleep( secs )
#if defined(INTEL)
        use ifport
#endif
        integer, intent(in) :: secs
#if defined(INTEL)
        integer  :: msecs
        msecs = 1000*secs
        call sleepqq(msecs)  !! milliseconds 
#else
        call sleep(secs) !! intrinsic
#endif        
    end subroutine simple_sleep

    integer function simple_chmod(pathname, mode  )
#if defined(INTEL)
        use ifport
#endif
        character(len=*), intent(in) :: pathname, mode
#if 1
        
        simple_chmod = chmod(pathname, mode)
        call simple_error_check()
#else
        call chmod(pathname, mode, status=simple_chmod) !! intrinsic
#endif        
    end function simple_chmod


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


end module simple_syslib

