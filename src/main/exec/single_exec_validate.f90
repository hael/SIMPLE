module single_exec_validate
use simple_cmdline,                 only: cmdline
use single_commanders_experimental, only: commander_cavgsproc_nano, commander_cavgseoproc_nano, commander_ptclsproc_nano
implicit none

public :: exec_validate_commander
private

type(commander_cavgseoproc_nano) :: xcavgseoproc
type(commander_cavgsproc_nano)   :: xcavgsproc
type(commander_ptclsproc_nano)   :: xptclsproc

contains

    subroutine exec_validate_commander(which, cline, l_silent, l_did_execute)
        character(len=*), intent(in)    :: which
        class(cmdline),   intent(inout) :: cline
        logical,          intent(inout) :: l_did_execute
        logical,          intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'cavgseoproc_nano' )
                call xcavgseoproc%execute(cline)
            case( 'cavgsproc_nano' )
                call xcavgsproc%execute(cline)
            case( 'ptclsproc_nano' )
                call xptclsproc%execute(cline)
            case default    
                l_did_execute = .false.
        end select
    end subroutine exec_validate_commander

end module single_exec_validate
