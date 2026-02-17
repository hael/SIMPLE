module single_exec_nano3D
use simple_string,            only: string
use simple_cmdline,           only: cmdline
use single_commanders_nano3D, only: commander_refine3D_nano, commander_autorefine3D_nano
use simple_exec_helpers,      only: restarted_exec
implicit none

public :: exec_nano3D_commander
private

type(commander_autorefine3D_nano) :: xautorefine3D_nano
type(commander_refine3D_nano)     :: xrefine3D_nano

contains

    subroutine exec_nano3D_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'autorefine3D_nano' )
                if( cline%defined('nrestarts') )then
                    call restarted_exec(cline, string('autorefine3D_nano'), string('single_exec'))
                else
                    call xautorefine3D_nano%execute(cline)
                endif
            case( 'refine3D_nano' )
                call xrefine3D_nano%execute(cline)
            case default    
                l_did_execute = .false.
        end select
    end subroutine exec_nano3D_commander

end module single_exec_nano3D
