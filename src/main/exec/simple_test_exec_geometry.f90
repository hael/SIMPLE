!@descr: execution of test geometry processing commanders
module simple_test_exec_geometry
use simple_cmdline,                  only: cmdline
use simple_commanders_test_geometry, only: commander_test_angres
implicit none

public :: exec_test_geometry_commander
private

type(commander_test_angres) :: xangres

contains

    subroutine exec_test_geometry_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'angres' )
                call xangres%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_geometry_commander

end module simple_test_exec_geometry
