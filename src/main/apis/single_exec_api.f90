!@descr: Aggregated public API for single_exec
module single_exec_api
use simple_core_module_api
use simple_cmdline,                 only: cmdline, cmdline_err
use simple_commanders_atoms,        only: commander_pdb2mrc, commander_conv_atom_denoise, commander_atoms_stats, commander_atoms_register,&
commander_crys_score, commander_atoms_rmsd, commander_core_atoms_analysis, commander_detect_atoms
use simple_commanders_cavgs,        only: commander_map_cavgs_selection
use simple_commanders_cluster2D,    only: commander_ppca_denoise_classes
use simple_commanders_imgproc,      only: commander_estimate_diam
use simple_commanders_ori,         only: commander_vizoris
use simple_commanders_project_core, only: commander_new_project, commander_update_project, commander_print_project_info, commander_print_project_field,&
commander_extract_subproj
use simple_commanders_project_ptcl, only: commander_import_particles, commander_prune_project_distr
use simple_commanders_sim,          only: commander_simulate_atoms
use simple_exec_helpers,            only: script_exec, update_job_descriptions_in_project, restarted_exec
use simple_jiffys,                  only: simple_print_git_version, simple_print_timer
use simple_user_interface,          only: make_user_interface, list_single_prgs_in_ui
use single_commanders_experimental, only: commander_cavgsproc_nano, commander_cavgseoproc_nano, commander_ptclsproc_nano,&
commander_trajectory_make_projavgs, commander_tsegmaps_core_finder
use single_commanders_nano2D,       only: commander_analysis2D_nano, commander_center2D_nano, commander_cluster2D_nano
use single_commanders_nano3D,       only: commander_refine3D_nano, commander_autorefine3D_nano,&
commander_trajectory_reconstruct3D_distr, commander_trajectory_reconstruct3D_distr
use single_commanders_trajectory,   only: commander_import_trajectory, commander_track_particles_distr, commander_graphene_subtr,&
commander_trajectory_denoise, commander_extract_substk, commander_trajectory_swap_stack
use single_commanders_tseries,      only: commander_tseries_import, commander_tseries_make_pickavg, commander_tseries_motion_correct_distr
end module single_exec_api
