module simple_exec_volume
use simple_cmdline,           only: cmdline
use simple_commanders_volops, only: commander_centervol, commander_reproject, commander_volops
implicit none

public :: exec_volume_commander
private

type(commander_centervol) :: xcenter
type(commander_reproject) :: xreproject
type(commander_volops)    :: xvolops

contains

    subroutine exec_volume_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'center' )
                call xcenter%execute(cline)
            case( 'reproject' )
                call xreproject%execute(cline)
            case( 'volops' )
                call xvolops%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_volume_commander

end module simple_exec_volume
