module single_ui_api
use simple_ui_program, only: ui_program
implicit none
type(ui_program), target :: analysis2D_nano
type(ui_program), target :: atoms_register
type(ui_program), target :: atoms_stats
type(ui_program), target :: autorefine3D_nano
type(ui_program), target :: cavgseoproc_nano
type(ui_program), target :: cavgsproc_nano
type(ui_program), target :: center2D_nano
type(ui_program), target :: cluster2D_nano
type(ui_program), target :: conv_atom_denoise
type(ui_program), target :: crys_score
type(ui_program), target :: denoise_trajectory
type(ui_program), target :: detect_atoms
type(ui_program), target :: estimate_diam
type(ui_program), target :: extract_substk
type(ui_program), target :: graphene_subtr
type(ui_program), target :: ptclsproc_nano
type(ui_program), target :: refine3D_nano
type(ui_program), target :: simulate_atoms
type(ui_program), target :: tseries_atoms_rmsd
type(ui_program), target :: tseries_core_atoms_analysis
type(ui_program), target :: tseries_import
type(ui_program), target :: tseries_import_particles
type(ui_program), target :: tseries_make_pickavg
type(ui_program), target :: tseries_motion_correct
type(ui_program), target :: track_particles
type(ui_program), target :: trajectory_core_finder
type(ui_program), target :: trajectory_make_projavgs
type(ui_program), target :: trajectory_reconstruct3D
type(ui_program), target :: trajectory_swap_stack
end module single_ui_api
