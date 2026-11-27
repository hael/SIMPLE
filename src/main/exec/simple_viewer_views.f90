!@descr: execution of the simple_viewer commanders
module simple_viewer_views
use simple_cmdline,             only: cmdline
use simple_commanders_select2D, only: commander_select2D
use simple_commanders_pick_positions, only: commander_pick_positions
implicit none

public :: exec_viewer_commander
private

type(commander_select2D) :: xselect2D
type(commander_pick_positions) :: xpick_positions

contains

    subroutine exec_viewer_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'select2D' )
                call xselect2D%execute(cline)
            case( 'pick_positions' )
                call xpick_positions%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_viewer_commander

end module simple_viewer_views
