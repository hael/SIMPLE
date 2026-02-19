!@descr: execution of test optimize processing commanders
module simple_test_exec_optimize
use simple_cmdline,                  only: cmdline
use simple_commanders_test_optimize, only: commander_test_lbfgsb, commander_test_lbfgsb_cosine, &
                                           commander_test_lplims, commander_test_lpstages_test, &
                                           commander_test_opt_lp, commander_test_tree_srch
implicit none

public :: exec_test_optimize_commander
private

type(commander_test_lbfgsb)        :: xlbfgsb
type(commander_test_lbfgsb_cosine) :: xlbfgsb_cosine
type(commander_test_lplims)        :: xlplims
type(commander_test_lpstages_test) :: xlpstages_test
type(commander_test_opt_lp)        :: xopt_lp
type(commander_test_tree_srch)     :: xtree_srch

contains

    subroutine exec_test_optimize_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'lbfgsb' )
                call xlbfgsb%execute(cline)
            case( 'lbfgsb_cosine' )
                call xlbfgsb_cosine%execute(cline)
            case( 'lplims' )
                call xlplims%execute(cline)
            case( 'lpstages_test' )
                call xlpstages_test%execute(cline)
            case( 'opt_lp' )
                call xopt_lp%execute(cline)
            case( 'tree_srch' )
                call xtree_srch%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_optimize_commander

end module simple_test_exec_optimize
