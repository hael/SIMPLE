!@descr: "validation" UI api (concrete implementation)
module simple_ui_api_validation
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none
public :: register_ui_validation

type(ui_program), target :: model_validation
type(ui_program), target :: mini_stream
type(ui_program), target :: check_refpick

contains

    subroutine register_ui_validation(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('model_validation', model_validation, prgtab)
        call add_ui_program('mini_stream',      mini_stream,      prgtab)
        call add_ui_program('check_refpick',    check_refpick,    prgtab)
    end subroutine register_ui_validation

end module simple_ui_api_validation
