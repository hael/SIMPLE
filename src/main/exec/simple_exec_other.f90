!@descr: execution of other commanders
module simple_exec_other
use simple_cmdline,                 only: cmdline
use simple_commanders_atoms,        only: commander_cif2pdb
use simple_commanders_distr,        only: commander_split
use simple_commanders_misc,         only: commander_fractionate_movies_distr
use simple_commanders_project_ptcl, only: commander_split_stack
implicit none

public :: exec_other_commander
private

type(commander_cif2pdb)                  :: xcif2pdb
type(commander_fractionate_movies_distr) :: xfractionate_movies
type(commander_split)                    :: xsplit
type(commander_split_stack)              :: xsplit_stack

contains

    subroutine exec_other_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case('cif2pdb')
                call xcif2pdb%execute(cline)
            case('fractionate_movies')
                call xfractionate_movies%execute(cline)
            case( 'split' )
                call xsplit%execute(cline)
            case( 'split_stack' )
                call xsplit_stack%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_other_commander

end module simple_exec_other
