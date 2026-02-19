!@descr: execution of masking commanders
module simple_exec_mask
use simple_cmdline,         only: cmdline
use simple_commanders_mask, only: commander_auto_spher_mask, commander_mask, commander_automask2D
implicit none

public :: exec_mask_commander
private

type(commander_auto_spher_mask) :: xauto_spher_mask
type(commander_automask2D)      :: xautomask2D
type(commander_mask)            :: xmask

contains

    subroutine exec_mask_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'auto_spher_mask' )
                call xauto_spher_mask%execute(cline)
            case( 'automask2D' )
                call xautomask2D%execute(cline)
            case( 'mask' )
                call xmask%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_mask_commander

end module simple_exec_mask
