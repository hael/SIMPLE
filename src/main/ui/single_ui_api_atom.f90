!@descr: "single_ui_api_atom" UI api (concrete implementation)
module single_ui_api_atom
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: atoms_register
type(ui_program), target :: atoms_rmsd
type(ui_program), target :: atoms_stats
type(ui_program), target :: core_atoms_analysis
type(ui_program), target :: crys_score
type(ui_program), target :: detect_atoms
type(ui_program), target :: simulate_atoms

contains

    subroutine register_single_ui_atom(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('atoms_register', atoms_register, prgtab)
        call add_ui_program('atoms_rmsd', atoms_rmsd, prgtab)
        call add_ui_program('atoms_stats', atoms_stats, prgtab)
        call add_ui_program('core_atoms_analysis', core_atoms_analysis, prgtab)
        call add_ui_program('crys_score', crys_score, prgtab)
        call add_ui_program('detect_atoms', detect_atoms, prgtab)
        call add_ui_program('simulate_atoms', simulate_atoms, prgtab)
    end subroutine register_single_ui_atom

end module single_ui_api_atom
