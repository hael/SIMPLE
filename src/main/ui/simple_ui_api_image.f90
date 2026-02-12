!@descr: "image" UI api (concrete implementation)
module simple_ui_api_image
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: binarize
type(ui_program), target :: convert
type(ui_program), target :: ctf_phaseflip
type(ui_program), target :: ctfops
type(ui_program), target :: scale
type(ui_program), target :: select_
type(ui_program), target :: stack
type(ui_program), target :: stackops

contains

    subroutine register_simple_ui_image(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('binarize',      binarize,      prgtab)
        call add_ui_program('convert',       convert,       prgtab)
        call add_ui_program('ctf_phaseflip', ctf_phaseflip, prgtab)
        call add_ui_program('ctfops',        ctfops,        prgtab)
        call add_ui_program('scale',         scale,         prgtab)
        call add_ui_program('select_',       select_,       prgtab)
        call add_ui_program('stack',         stack,         prgtab)
        call add_ui_program('stackops',      stackops,      prgtab)
    end subroutine register_simple_ui_image

end module simple_ui_api_image
