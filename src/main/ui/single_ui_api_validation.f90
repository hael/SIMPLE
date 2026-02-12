!@descr: "single_ui_api_validation" UI api (concrete implementation)
module single_ui_api_validation
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: cavgseoproc_nano
type(ui_program), target :: cavgsproc_nano
type(ui_program), target :: ptclsproc_nano

contains

    subroutine register_single_ui_validation(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('cavgseoproc_nano', cavgseoproc_nano, prgtab)
        call add_ui_program('cavgsproc_nano', cavgsproc_nano, prgtab)
        call add_ui_program('ptclsproc_nano', ptclsproc_nano, prgtab)
    end subroutine register_single_ui_validation

end module single_ui_api_validation
