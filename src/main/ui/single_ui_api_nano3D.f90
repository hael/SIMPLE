!@descr: "single_ui_api_nano3D" UI api (concrete implementation)
module single_ui_api_nano3D
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: autorefine3D_nano
type(ui_program), target :: refine3D_nano

contains

    subroutine register_single_ui_nano3D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('autorefine3D_nano', autorefine3D_nano, prgtab)
        call add_ui_program('refine3D_nano', refine3D_nano, prgtab)
    end subroutine register_single_ui_nano3D

end module single_ui_api_nano3D
