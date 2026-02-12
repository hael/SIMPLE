!@descr: "resolution" UI api (concrete implementation)
module simple_ui_api_resolution
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: fsc
type(ui_program), target :: clin_fsc

contains

    subroutine register_ui_resolution(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('fsc',      fsc,      prgtab)
        call add_ui_program('clin_fsc', clin_fsc, prgtab)
    end subroutine register_ui_resolution

end module simple_ui_api_resolution
