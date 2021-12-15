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
            call backtrace()
            stop
        endif
    end subroutine simple_exception

    subroutine simple_error_check(io_stat, msg)
        integer,                    intent(in) :: io_stat
        character(len=*), optional, intent(in) :: msg
        if (io_stat==IOSTAT_END)  write(logfhandle,'(a,1x,I0 )') 'ERRCHECK: EOF reached, end-of-file reached IOS# ', io_stat
        if (io_stat==IOSTAT_EOR)  write(logfhandle,'(a,1x,I0 )') 'ERRCHECK: EOR reached, read was short, IOS# ', io_stat
        if( io_stat /= 0 ) then
            write(logfhandle,'(a,1x,I0 )') 'ERROR: IOS# ', io_stat
            if(present(msg)) write(logfhandle,'(a)') trim(adjustl(msg))
        endif
    end subroutine simple_error_check

end module simple_error
