!@descr: execution of orientation commanders
module simple_exec_ori
use simple_cmdline, only: cmdline
use simple_commanders_ori, only: commander_make_oris, commander_oriops, commander_oristats, commander_vizoris
implicit none

public :: exec_ori_commander
private

type(commander_make_oris) :: xmake_oris
type(commander_oriops)    :: xoriops
type(commander_oristats)  :: xoristats
type(commander_vizoris)   :: xvizoris   

contains

    subroutine exec_ori_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'make_oris' )
                call xmake_oris%execute(cline)
            case( 'oriops' )
                call xoriops%execute(cline)
            case( 'oristats' )
                call xoristats%execute(cline)
            case( 'vizoris' )
                call xvizoris%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_ori_commander

end module simple_exec_ori
