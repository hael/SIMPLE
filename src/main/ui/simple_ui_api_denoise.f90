!@descr: "denoise" UI api (concrete implementation)
module simple_ui_api_denoise
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: icm2D
type(ui_program), target :: icm3D
type(ui_program), target :: ppca_denoise
type(ui_program), target :: ppca_denoise_classes
type(ui_program), target :: ppca_volvar

contains

    subroutine register_ui_denoise(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('icm2D',                icm2D,                prgtab)
        call add_ui_program('icm3D',                icm3D,                prgtab)
        call add_ui_program('ppca_denoise',         ppca_denoise,         prgtab)
        call add_ui_program('ppca_denoise_classes', ppca_denoise_classes, prgtab)
        call add_ui_program('ppca_volvar',          ppca_volvar,          prgtab)
    end subroutine register_ui_denoise

end module simple_ui_api_denoise
