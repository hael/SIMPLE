!@descr: execution of test parallel processing commanders
module simple_test_exec_parallel
use simple_cmdline,                  only: cmdline
use simple_commanders_test_parallel, only: commander_test_coarrays
implicit none

public :: exec_test_parallel_commander
private

type(commander_test_coarrays)           :: xcoarrays

contains

    subroutine exec_test_parallel_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'coarrays' )
                call xcoarrays%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_parallel_commander

end module simple_test_exec_parallel
