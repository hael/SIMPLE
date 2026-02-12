module simple_ui_api_abinitio3D
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: abinitio3D
type(ui_program), target :: abinitio3D_cavgs
type(ui_program), target :: estimate_lpstages
type(ui_program), target :: multivol_assign
type(ui_program), target :: noisevol

contains

    subroutine register_simple_ui_abinitio3D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('abinitio3D',        abinitio3D,        prgtab)
        call add_ui_program('abinitio3D_cavgs',  abinitio3D_cavgs,  prgtab)
        call add_ui_program('estimate_lpstages', estimate_lpstages, prgtab)
        call add_ui_program('multivol_assign',   multivol_assign,   prgtab)
        call add_ui_program('noisevol',          noisevol,          prgtab)
    end subroutine register_simple_ui_abinitio3D

end module simple_ui_api_abinitio3D
