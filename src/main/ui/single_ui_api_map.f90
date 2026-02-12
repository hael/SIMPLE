!@descr: "single_ui_api_map" UI api (concrete implementation)
module single_ui_api_map
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: conv_atom_denoise
type(ui_program), target :: tsegmaps_core_finder

contains

    subroutine register_single_ui_map(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('conv_atom_denoise', conv_atom_denoise, prgtab)
        call add_ui_program('tsegmaps_core_finder', tsegmaps_core_finder, prgtab)
    end subroutine register_single_ui_map

end module single_ui_api_map
