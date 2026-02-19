!@descr: execution of test highlevel processing commanders
module simple_test_exec_highlevel
use simple_cmdline,                   only: cmdline
use simple_commanders_test_highlevel, only: commander_test_mini_stream, commander_test_simulated_workflow
implicit none

public :: exec_highlevel_commander
private

type(commander_test_mini_stream)        :: xmini_stream
type(commander_test_simulated_workflow) :: xsimulated_workflow

contains

    subroutine exec_highlevel_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'mini_stream' )
                call xmini_stream%execute(cline)
            case( 'simulated_workflow' )
                call xsimulated_workflow%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_highlevel_commander

end module simple_test_exec_highlevel
