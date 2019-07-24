module simple_error
use simple_defs
use, intrinsic :: iso_fortran_env
implicit none

contains

    subroutine simple_exception( msg, file, line, l_stop )
        character(len=*),  intent(in) :: msg, file
        integer,           intent(in) :: line
        logical, optional, intent(in) :: l_stop
        logical :: ll_stop
        ll_stop = .true.
        if( present(l_stop) ) ll_stop = l_stop
        if( ll_stop )then
            write(logfhandle,'(A)', advance='no') 'ERROR! '//trim(msg)
        else
            write(logfhandle,'(A)', advance='no') 'WARNING! '//trim(msg)
        endif
        if( l_distr_exec_glob )then
            write(logfhandle,'(A,I5)') '; '//trim(file)//'; line: ', line, '; part: ', part_glob, ' of distributed execution'
        else
            write(logfhandle,'(A,I5)') '; '//trim(file)//'; line: ', line
        endif
        if( ll_stop )then
#if defined(GNU) && defined(_DEBUG)
                    call backtrace()
#endif
            call exit(EXIT_FAILURE)
        endif
    end subroutine simple_exception

    !> \brief  is for checking allocation
    subroutine allocchk( message, alloc_err, file, line, iomsg )
        character(len=*), intent(in)           :: message
        integer,          intent(in), optional :: alloc_err
        character(len=*), intent(in), optional :: file !< filename of caller
        integer,          intent(in), optional :: line !< line number from calling file
        character(len=*), intent(in), optional :: iomsg !< IO message
        integer                                :: alloc_status
        alloc_status=alloc_stat    !! global variable from simple_defs
        if(present(alloc_err))alloc_status=alloc_err
        if (alloc_status/=0)then
            write(logfhandle,'(a)') 'ERROR: Allocation failure!'
            call simple_error_check(alloc_status)
            if(present(iomsg))&
                write(logfhandle,'("IO Message ",A)') trim(adjustl(iomsg))
            if(present(file).and.present(line))&
                write(logfhandle,'("Stopping in file ",/,A,/," at line ",I0)') file,line
            call simple_exception(message,__FILENAME__,__LINE__)
        endif
    end subroutine allocchk

    subroutine simple_error_check(io_stat, msg)
        integer,                    intent(in) :: io_stat
        character(len=*), optional, intent(in) :: msg

        !! Intel and GNU
        if (io_stat==IOSTAT_END)  write(logfhandle,'(a,1x,I0 )') 'ERRCHECK: EOF reached, end-of-file reached IOS# ', io_stat
        if (io_stat==IOSTAT_EOR)  write(logfhandle,'(a,1x,I0 )') 'ERRCHECK: EOR reached, read was short, IOS# ', io_stat
        if( io_stat /= 0 ) then
            write(logfhandle,'(a,1x,I0 )') 'ERROR: IOS# ', io_stat
            if(present(msg)) write(logfhandle,'(a)') trim(adjustl(msg))
            !! do not stop yet -- let the fopen/fclose call to finish
            !! stop
        endif
    end subroutine simple_error_check

end module simple_error
