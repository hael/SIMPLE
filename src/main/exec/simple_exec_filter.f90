!@descr: execution of filtering commanders
module simple_exec_filter
use simple_cmdline,             only: cmdline
use simple_commanders_resolest, only: commander_uniform_filter2D, commander_uniform_filter3D
use simple_commanders_imgops,   only: commander_filter
implicit none

public :: exec_filter_commander
private

type(commander_filter)           :: xfilter
type(commander_uniform_filter2D) :: xuniform_filter2D
type(commander_uniform_filter3D) :: xuniform_filter3D

contains

    subroutine exec_filter_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'filter' )
                call xfilter%execute(cline)
            case( 'uniform_filter2D' )
                call xuniform_filter2D%execute(cline)
            case( 'uniform_filter3D' )
                call xuniform_filter3D%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_filter_commander

end module simple_exec_filter
