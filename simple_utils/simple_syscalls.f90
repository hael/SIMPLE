!==Module simple_syscalls
!
! simple_syscalls is a little module for calculating the relative and actual CPU-time.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution or modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund, 2009-10-01.
! 
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-10-01
!
module simple_syscalls
use simple_jiffys ! singleton
use simple_defs   ! singleton
implicit none

private :: raise_sys_error

contains

    !> is the fortran 90 variant of the classic dtime
    real function dtime( time )
        real                  :: time(2)
        double precision,save :: last_time = 0
        double precision      :: this_time
        intrinsic cpu_time
        call cpu_time(this_time)
        time(1) = real(this_time-last_time)
        time(2) = 0.
        dtime = time(1)
        last_time = this_time
    end function dtime

    !> is the fortran 90 variant of the classic etime
    real function etime( time )
        real :: time(2)
        call cpu_time(etime)
        time(1) = etime
        time(2) = 0
    end function etime

    !> is for getting the actual cputime
    function getabscpu( lprint ) result( actual )
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: actual
        actual = etime( tarray )
        if( lprint )then
            write(*,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        endif
    end function getabscpu

    !> is for getting the relative cpu-time
    function getdiffcpu( lprint ) result( delta )   
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: delta
        delta = dtime( tarray )
        if( lprint )then
            write(*,'(A,F9.2)') 'Relative cpu-time:', delta
        endif
    end function getdiffcpu

    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline )
        character(len=*), intent(in) :: cmdline
        integer               :: estat, cstat
        character(len=STDLEN) :: cmsg
        call execute_command_line( trim(adjustl(cmdline)), exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
        call raise_sys_error( cmdline, estat, cstat, cmsg )
    end subroutine exec_cmdline

    !>  Handles error from system call
    subroutine raise_sys_error( cmd, exitstat, cmdstat, cmdmsg )
        integer,               intent(in) :: exitstat, cmdstat
        character(len=STDLEN), intent(in) :: cmd, cmdmsg
        logical :: err
        err = .false.
        if( exitstat /= 0 )then
            write(*,*)'System error', exitstat,' for command: ', trim(adjustl(cmd))
            err = .true.
        endif 
        if( cmdstat /= 0 )then
            write(*,*)'cmdstat /= 0, command could not be executed: ', trim(adjustl(cmd))
            err = .true.
        endif
        if( err ) write(*,*) trim(cmdmsg)
    end subroutine raise_sys_error

    function sys_get_env_var( name ) result( varval )
        character(len=*), intent(in)  :: name
        character(len=STDLEN)         :: value
        character(len=:), allocatable :: varval
        integer :: length, status
        call get_environment_variable( trim(name), value=value, length=length, status=status)
        if( status == -1 ) write(*,*) 'value string too short; simple_syscalls :: get_env_var'
        if( status ==  1 ) write(*,*) 'environment variable: ', trim(name), ' is not defined; simple_syscalls :: get_env_var'
        if( status ==  2 ) write(*,*) 'environment variables not supported by system; simple_syscalls :: get_env_var'
        if( length ==  0 .or. status /= 0 ) return
        allocate(varval, source=trim(value))
    end function sys_get_env_var

end module simple_syscalls
