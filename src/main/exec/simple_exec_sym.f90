!@descr: execution of symmetry-related commanders
module simple_exec_sym
use simple_cmdline,           only: cmdline
use simple_commanders_volops, only: commander_symaxis_search, commander_symmetry_test, commander_symmetrize_map
implicit none

public :: exec_sym_commander
private

type(commander_symaxis_search) :: xsymsrch
type(commander_symmetrize_map) :: xsymmetrize_map
type(commander_symmetry_test)  :: xsymtst

contains

    subroutine exec_sym_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'symaxis_search' )
                call xsymsrch%execute(cline)
            case( 'symmetrize_map' )
                call xsymmetrize_map%execute(cline)
            case( 'symmetry_test' )
                call xsymtst%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_sym_commander

end module simple_exec_sym
