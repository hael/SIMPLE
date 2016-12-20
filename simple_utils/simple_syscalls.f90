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

    real function dtime( time )
    ! is the fortran 90 variant of the classical dtime
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

    real function etime( time )
    ! is the fortran 90 variant of the classical etime
        real :: time(2)
        call cpu_time(etime)
        time(1) = etime
        time(2) = 0
    end function

    function getabscpu( lprint ) result( actual )
    ! is for getting the actual cputime
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: actual
        actual = etime( tarray )
        if( lprint )then
            write(*,'(A,2X,F9.2)') 'Actual cpu-time:', actual
        endif
    end function

    function getdiffcpu( lprint ) result( delta )
    ! is for getting the relative cpu-time
        logical, intent(in) :: lprint
        real                :: tarray(2)
        real                :: delta
        delta = dtime( tarray )
        if( lprint )then
            write(*,'(A,F9.2)') 'Relative cpu-time:', delta
        endif
    end function

    !>  Wrapper for system call
    subroutine exec_cmdline( cmdline )
        character(len=*), intent(inout) :: cmdline
        integer               :: estat, cstat
        character(len=STDLEN) :: cmsg
        call execute_command_line( trim(cmdline), exitstat=estat, cmdstat=cstat, cmdmsg=cmsg)
        call raise_sys_error( cmdline, estat, cstat, cmsg )
    end subroutine exec_cmdline

    !>  Handles error from system call
    subroutine raise_sys_error( cmd, exitstat, cmdstat, cmdmsg )
        integer,               intent(in) :: exitstat, cmdstat
        character(len=STDLEN), intent(in) :: cmd, cmdmsg
        logical :: dostop
        dostop = .false.
        if( exitstat /= 0 )then
            print *, 'System error', exitstat,' for command: ', cmd
            write(*,*)'System error', exitstat,' for command: ', cmd
            !dostop = .true.
        endif 
        if( cmdstat /= 0 )then
            print *, 'Command could not be executed: ', cmd
            write(*,*)'Command could not be executed: ', cmd
            !dostop = .true.
        endif 
        if( dostop )stop
    end subroutine raise_sys_error

    !> interface to unix mkdir
    subroutine sys_mkdir( dir, args )
        character(len=*),           intent(in) :: dir
        character(len=*), optional, intent(in) :: args
        character(len=STDLEN) :: cmd
        if( .not. file_exists( trim(dir) ) )then
            if( present(args) )then
                cmd = 'mkdir '//trim(args)//' '//trim(dir)
            else
                cmd = 'mkdir '//trim(dir)
            endif
            call exec_cmdline( cmd )
        endif
    end subroutine sys_mkdir

    !> interface to unix mkdir
    subroutine sys_cp( source, dest, args )
        character(len=*),           intent(in) :: source, dest
        character(len=*), optional, intent(in) :: args
        character(len=STDLEN) :: cmd
        if( file_exists( source ) )then
            if( present(args) )then
                cmd = 'cp '//trim(args)//' '//trim(source)//' '//trim(dest)
            else
                cmd = 'cp '//' '//trim(source)//' '//trim(dest)
            endif
            call exec_cmdline( cmd )
        endif
    end subroutine sys_cp

end module simple_syscalls
