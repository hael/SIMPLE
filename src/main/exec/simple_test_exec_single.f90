!@descr: execution of test single processing commanders
module simple_test_exec_single
use simple_cmdline,                only: cmdline
use simple_commanders_test_single, only: commander_test_atoms_stats, commander_test_detect_atoms, commander_test_simulate_nanoprticle, commander_test_single_workflow
implicit none

public :: exec_test_single_commander
private

type(commander_test_atoms_stats)          :: xatoms_stats
type(commander_test_detect_atoms)         :: xdetect_atoms
type(commander_test_simulate_nanoprticle) :: xsimulate_nanoprticle
type(commander_test_single_workflow)      :: xsingle_workflow

contains

    subroutine exec_test_single_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'atoms_stats' )
                call xatoms_stats%execute(cline)
            case( 'detect_atoms' )
                call xdetect_atoms%execute(cline)
            case( 'simulate_nanoprticle' )
                call xsimulate_nanoprticle%execute(cline)
            case( 'single_workflow' )
                call xsingle_workflow%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_single_commander

end module simple_test_exec_single
