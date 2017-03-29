!==Module simple_restart
!
! This is a small module for handling tasks related to restart
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Authors:* Cyril Reboul & Hans Elmlund 2017
!
module simple_restart
use simple_defs 
use simple_strings, only: str_has_substr, parse, split, int2str_pad, str2int
use simple_cmdline, only: cmdline
implicit none

public  :: check_restart, parse_restart
private
integer, parameter :: MAXNKEYS=100

contains

    !> detemines whether this is a restart scenario 
    subroutine check_restart( line, is_restart )
        character(len=STDLEN), intent(inout) :: line
        logical,               intent(inout) :: is_restart
        character(len=STDLEN) :: args(MAXNKEYS)
        integer               :: nargs
        is_restart = .false.
        call parse( line, ' ', args, nargs )
        if( nargs /= 3 )return
        if( str_has_substr(args(2), 'prg=').and.str_has_substr(args(3), 'startit=') )then
            is_restart = .true.
        else if( str_has_substr(args(2), 'prg=').and.str_has_substr(args(3), 'fname=') )then
            is_restart = .true.
        endif
    end subroutine check_restart

    !> parses restart files and produces corresponding cmdline object
    subroutine parse_restart( prg, line, cline, keys_required, keys_optional )
        character(len=*),           intent(in)    :: prg
        character(len=STDLEN),      intent(inout) :: line
        type(cmdline),              intent(inout) :: cline
        character(len=*), optional, intent(inout) :: keys_required(:), keys_optional(:)
        character(len=STDLEN) :: args(MAXNKEYS), arg, fname
        integer               :: nargs, iostat, iter
        call parse( line, ' ', args, nargs )
        if( nargs /= 3 )stop 'Invalid command line for restart; simple_restart::parse_restart 1'
        call split( args(2), '=', arg )
        if( trim(arg).ne.'prg' )stop 'Invalid command line for restart; simple_restart::parse_restart 2'
        if( trim(args(2)).ne.trim(prg) )stop 'PRG key arguments are incompatible; simple_restart::parse_restart 3'
        call split( args(3), '=', arg )
        select case( trim(arg) )
            case('startit')
                ! restart from an iteration
                call str2int( args(3), iostat, iter )
                if( iostat.ne.0 )stop 'invalid argument to startit; simple_restart::parse_restart'
                fname = 'prime3D_restart_iter' // int2str_pad( iter, 3 ) // '.txt'
            case('fname')
                ! restart from file
                fname = trim( adjustl(args(3)) )
            case DEFAULT
                print *,'Invalid restart key: ', trim(arg)
                stop
        end select
        call cline%read( fname, keys_required, keys_optional )
        if( cline%defined('startit') )then
            iter = nint( cline%get_rarg('startit') )
            call cline%set( 'startit', real(iter+1) )
        else
            call cline%set( 'startit', 1. )
        endif
    end subroutine parse_restart

end module simple_restart
