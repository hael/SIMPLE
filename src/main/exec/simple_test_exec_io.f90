!@descr: execution of test input/output processing commanders (manual, user-data cases)
module simple_test_exec_io
use simple_cmdline,            only: cmdline
use simple_commanders_test_io, only: commander_test_mrc2jpeg, commander_test_mrc_validate
implicit none

public :: exec_test_io_commander
private

type(commander_test_mrc2jpeg)       :: xmrc2jpeg
type(commander_test_mrc_validate)   :: xmrc_validate

contains

    subroutine exec_test_io_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'mrc2jpeg' )
                call xmrc2jpeg%execute(cline)
            case( 'mrc_validate' )
                call xmrc_validate%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_io_commander

end module simple_test_exec_io
