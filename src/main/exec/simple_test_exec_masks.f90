!@descr: execution of test masks processing commanders
module simple_test_exec_masks
use simple_cmdline,               only: cmdline
use simple_commanders_test_masks, only: commander_test_nano_mask, commander_test_score_volume_shape
implicit none

public :: exec_test_masks_commander
private

type(commander_test_nano_mask)               :: xnano_mask
type(commander_test_score_volume_shape)      :: xscore_volume_shape

contains

    subroutine exec_test_masks_commander( which, cline, l_silent, l_did_execute )
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'nano_mask' )
                call xnano_mask%execute(cline)
            case( 'score_volume_shape' )
                call xscore_volume_shape%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_test_masks_commander

end module simple_test_exec_masks
