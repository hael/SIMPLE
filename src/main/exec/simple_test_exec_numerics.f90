!@descr: execution of test numerics processing commanders
module simple_test_exec_numerics
use simple_cmdline,                  only: cmdline
use simple_commanders_test_numerics, only: commander_test_eigh_test, commander_test_kbinterpol_fast, &
                                           commander_test_maxnloc_test, commander_test_neigh
implicit none

public :: exec_test_numerics_commander
private

type(commander_test_eigh_test)       :: xeigh_test
type(commander_test_kbinterpol_fast) :: xkbinterpol_fast
type(commander_test_maxnloc_test)    :: xmaxnloc_test
type(commander_test_neigh)           :: xneigh

contains

    subroutine exec_test_numerics_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'eigh_test' )
                call xeigh_test%execute(cline)
            case( 'kbinterpol_fast' )
                call xkbinterpol_fast%execute(cline)
            case( 'maxnloc_test' )
                call xmaxnloc_test%execute(cline)
            case( 'neigh' )
                call xneigh%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_numerics_commander

end module simple_test_exec_numerics
