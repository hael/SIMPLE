module simple_exec_res
use simple_cmdline,             only: cmdline
use simple_commanders_resolest, only: commander_fsc, commander_clin_fsc
implicit none

public :: exec_res_commander
private

type(commander_fsc)      :: xfsc
type(commander_clin_fsc) :: xclin_fsc

contains

    subroutine exec_res_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'fsc' )
                call xfsc%execute(cline)
            case( 'clin_fsc' )
                call xclin_fsc%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_res_commander

end module simple_exec_res
