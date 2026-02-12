!@descr: "filter" UI api (concrete implementation)
module simple_ui_api_filter
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: filter
type(ui_program), target :: uniform_filter2D
type(ui_program), target :: uniform_filter3D

contains

    subroutine register_simple_ui_filter(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('filter',           filter,           prgtab)
        call add_ui_program('uniform_filter2D', uniform_filter2D, prgtab)
        call add_ui_program('uniform_filter3D', uniform_filter3D, prgtab)
    end subroutine register_simple_ui_filter

end module simple_ui_api_filter
