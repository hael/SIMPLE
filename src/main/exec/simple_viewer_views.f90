!@descr: execution of the simple_viewer commanders (select2D)
module simple_viewer_views
use simple_cmdline,             only: cmdline
use simple_commanders_select2D, only: commander_select2D
implicit none

public :: exec_viewer_commander
private

type(commander_select2D) :: xselect2D

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
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_viewer_commander

end module simple_viewer_views
