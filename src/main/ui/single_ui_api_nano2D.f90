!@descr: "single_ui_api_nano2D" UI api (concrete implementation)
module single_ui_api_nano2D
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: analysis2D_nano
type(ui_program), target :: center2D_nano
type(ui_program), target :: cluster2D_nano
type(ui_program), target :: estimate_diam

contains

    subroutine register_single_ui_nano2D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('analysis2D_nano', analysis2D_nano, prgtab)
        call add_ui_program('center2D_nano', center2D_nano, prgtab)
        call add_ui_program('cluster2D_nano', cluster2D_nano, prgtab)
        call add_ui_program('estimate_diam', estimate_diam, prgtab)
    end subroutine register_single_ui_nano2D

end module single_ui_api_nano2D
