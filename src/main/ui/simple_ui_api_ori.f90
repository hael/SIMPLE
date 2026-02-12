!@descr: "ori" UI api (concrete implementation)
module simple_ui_api_ori
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: make_oris
type(ui_program), target :: orisops
type(ui_program), target :: oristats
type(ui_program), target :: vizoris

contains

    subroutine register_ui_ori(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('make_oris', make_oris, prgtab)
        call add_ui_program('orisops',   orisops,   prgtab)
        call add_ui_program('oristats',  oristats,  prgtab)
        call add_ui_program('vizoris',   vizoris,   prgtab)
    end subroutine register_ui_ori

end module simple_ui_api_ori
