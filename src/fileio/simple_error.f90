!@descr: exception handling
module simple_error
use simple_defs
use, intrinsic :: iso_fortran_env
implicit none

! simple_exception( msg, file, line )         hard error, never returns (THROW_HARD)
! simple_exception( msg, file, line, l_stop ) warning when l_stop is false (THROW_WARN)
! The three-argument form resolves to a NORETURN specific, so the optimizer treats the
! code after THROW_HARD as unreachable and does not report variables left unset on the
! error branch (-Wmaybe-uninitialized in Release builds).
interface simple_exception
    module procedure simple_exception_hard
    module procedure simple_exception_opt
end interface simple_exception

contains

    subroutine simple_exception_hard( msg, file, line )
!GCC$ ATTRIBUTES NORETURN :: simple_exception_hard
        character(len=*), intent(in) :: msg, file
        integer,          intent(in) :: line
        call simple_exception_opt(msg, file, line, .true.)
        error stop 1 ! not reached: simple_exception_opt stops when l_stop is true
    end subroutine simple_exception_hard

    subroutine simple_exception_opt( msg, file, line, l_stop )
        character(len=*), intent(in) :: msg, file
        integer,          intent(in) :: line
        logical,          intent(in) :: l_stop
        if( l_stop )then
            write(logfhandle,'(A)', advance='no') 'ERROR! '//trim(msg)
        else
            write(logfhandle,'(A)', advance='no') 'WARNING! '//trim(msg)
        endif
        if( l_distr_worker_glob )then
            write(logfhandle,'(A,I5)') '; '//trim(file)//'; line: ', line, '; part: ', part_glob, ' of distributed execution'
        else
            write(logfhandle,'(A,I5)') '; '//trim(file)//'; line: ', line
        endif
        if( l_stop )then
            ! best-effort removal of any on-disk cache this process owns; the hook
            ! nullifies itself on entry, so a throw during cleanup cannot recurse here
            if( associated(cache_cleanup_glob) ) call cache_cleanup_glob()
            call backtrace()
            error stop 1
        endif
    end subroutine simple_exception_opt

    subroutine simple_error_check(io_stat, msg)
        integer,                    intent(in) :: io_stat
        character(len=*), optional, intent(in) :: msg
        if( io_stat /= 0 .and. io_stat /= IOSTAT_END .and. io_stat /= IOSTAT_EOR ) then
            write(logfhandle,'(a,1x,I0 )') 'ERROR: IOS# ', io_stat
            if(present(msg)) write(logfhandle,'(a)') trim(adjustl(msg))
        endif
    end subroutine simple_error_check

end module simple_error