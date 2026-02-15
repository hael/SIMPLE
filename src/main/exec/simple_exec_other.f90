!@descr: execution of other commanders
module simple_exec_other
use simple_cmdline,                 only: cmdline
use simple_commanders_distr,        only: commander_split
use simple_commanders_project_ptcl, only: commander_split_stack
implicit none

public :: exec_other_commander
private

type(commander_split)       :: xsplit
type(commander_split_stack) :: xsplit_stack

contains

    subroutine exec_other_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'split' )
                call xsplit%execute(cline)
            case( 'split_stack' )
                call xsplit_stack%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_other_commander

end module simple_exec_other
