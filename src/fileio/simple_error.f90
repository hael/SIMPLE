module simple_error
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
    interface
        subroutine gerror(str)
            character*(*), intent(out) :: str
        end subroutine gerror
        integer function ierrno()
        end function ierrno
        subroutine perror(str)
            character*(*), intent(in) :: str
        end subroutine perror
    end interface
#endif


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

    function get_sys_error () result(err)
        integer :: err
        character(len=100) :: msg
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

            write(msg,'("Last detected error (SIMPLE_SYSLIB::SYSERROR) ",I0,":")') err
            call perror(trim(adjustl(msg)))

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

end module simple_error
