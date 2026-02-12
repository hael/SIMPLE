!@descr: "dock" UI api (concrete implementation)
module simple_ui_api_dock
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: dock_volpair
type(ui_program), target :: volanalyze

contains

    subroutine register_simple_ui_dock(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('dock_volpair', dock_volpair, prgtab)
        call add_ui_program('volanalyze',   volanalyze,   prgtab)
    end subroutine register_simple_ui_dock

end module simple_ui_api_dock
