!@descr: "mask" UI api (concrete implementation)
module simple_ui_api_mask
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: auto_spher_mask
type(ui_program), target :: automask2D
type(ui_program), target :: mask

contains

    subroutine register_ui_mask(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('auto_spher_mask', auto_spher_mask, prgtab)
        call add_ui_program('automask2D',      automask2D,      prgtab)
        call add_ui_program('mask',            mask,            prgtab)
    end subroutine register_ui_mask

end module simple_ui_api_mask
