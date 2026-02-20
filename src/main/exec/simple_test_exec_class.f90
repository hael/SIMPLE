!@descr: execution of test class processing commanders
module simple_test_exec_class
use simple_cmdline,               only: cmdline
use simple_commanders_test_class, only: commander_test_units
implicit none

public :: exec_test_class_commander
private

type(commander_test_units) :: xunits

contains

    subroutine exec_test_class_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'units' )
                call xunits%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_class_commander

end module simple_test_exec_class
