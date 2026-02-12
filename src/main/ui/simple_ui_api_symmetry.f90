!@descr: "symmetry" UI api (concrete implementation)
module simple_ui_api_symmetry
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: symaxis_search
type(ui_program), target :: symmetrize_map
type(ui_program), target :: symmetry_test

contains

    subroutine register_ui_symmetry(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('symaxis_search', symaxis_search, prgtab)
        call add_ui_program('symmetrize_map', symmetrize_map, prgtab)
        call add_ui_program('symmetry_test',  symmetry_test,  prgtab)
    end subroutine register_ui_symmetry

end module simple_ui_api_symmetry
