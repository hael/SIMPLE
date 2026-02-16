!@descr: execution of validation commanders
module simple_exec_validate
use simple_cmdline,             only: cmdline
use simple_commanders_validate, only: commander_mini_stream, commander_check_refpick
use simple_commanders_atoms,    only: commander_map2model_fsc, commander_model_validate
implicit none

public :: exec_validate_commander
private

type(commander_check_refpick)  :: xcheck_refpick
type(commander_mini_stream)    :: xmini_stream
type(commander_model_validate) :: xmodel_validate

contains

    subroutine exec_validate_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'check_refpick' )
                call xcheck_refpick%execute(cline)
            case( 'mini_stream' )
                call xmini_stream%execute(cline)
            case( 'model_validate' )
                call xmodel_validate%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_validate_commander

end module simple_exec_validate
