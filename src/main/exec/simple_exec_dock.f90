!@descr: execution of map docking commanders
module simple_exec_dock
use simple_cmdline,           only: cmdline
use simple_commanders_volops, only: commander_dock_volpair, commander_volanalyze
implicit none

public :: exec_dock_commander
private

type(commander_dock_volpair) :: xdock_volpair
type(commander_volanalyze)   :: xvolanalyze

contains

    subroutine exec_dock_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'dock_volpair' )
                call xdock_volpair%execute(cline)
            case( 'volanalyze' )
                call xvolanalyze%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_dock_commander

end module simple_exec_dock
