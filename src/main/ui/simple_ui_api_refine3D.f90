!@descr: "refine3D" UI api (concrete implementation)
module simple_ui_api_refine3D
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: refine3D
type(ui_program), target :: refine3D_auto
type(ui_program), target :: reconstruct3D
type(ui_program), target :: postprocess
type(ui_program), target :: automask

contains

    subroutine register_simple_ui_refine3D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('refine3D',      refine3D,      prgtab)
        call add_ui_program('refine3D_auto', refine3D_auto, prgtab)
        call add_ui_program('reconstruct3D', reconstruct3D, prgtab)
        call add_ui_program('postprocess',   postprocess,   prgtab)
        call add_ui_program('automask',      automask,      prgtab)
    end subroutine register_simple_ui_refine3D

end module simple_ui_api_refine3D
