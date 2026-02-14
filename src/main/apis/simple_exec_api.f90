!@descr: Aggregated public API for simple_exec
!! 
!! This module centralizes all commander types used by the simple_exec
!! front-end. It:
!!  - USEs the various simple_commanders_* modules
!!
!! simple_exec should depend only on this module (plus simple_exec_helpers),
!! rather than on individual simple_commanders_* modules. To add a new
!! command:
!!  1. USE its defining simple_commanders_* module here.
!!  2. Add the corresponding commander_* symbol to the PUBLIC list.
!!  3. Declare a type(commander_...) variable in simple_exec and wire it
!!     into the SELECT CASE(prg) dispatch.
module simple_exec_api
use simple_core_module_api
use simple_exec_helpers,   only: script_exec, restarted_exec, update_job_descriptions_in_project
use simple_jiffys,         only: simple_print_git_version, simple_print_timer
use simple_user_interface, only: make_user_interface, list_simple_prgs_in_ui, list_simple_test_prgs_in_ui
use iso_fortran_env,       only: output_unit
use simple_cmdline,        only: cmdline, cmdline_err

! imgproc commanders, standard image processing routines
use simple_commanders_imgproc, only: commander_ctfops, commander_ctf_phaseflip

use simple_commanders_cluster2D, only: commander_ppca_denoise_classes

! stkops commanders, image stack operations
use simple_commanders_stkops, only: commander_convert, commander_stack, commander_stackops

! imgops commanders, standard image operation routines
use simple_commanders_imgops, only: commander_binarize, commander_filter, commander_normalize, commander_ppca_denoise,&
commander_scale

! resoltest commanders, for resolution estimation and regularization testing
use simple_commanders_resolest, only: commander_fsc, commander_clin_fsc,&
commander_uniform_filter2D, commander_uniform_filter3D, commander_icm2D, commander_icm3D

! volops commanders, operations on volumes
use simple_commanders_volops, only: commander_symaxis_search, commander_symmetry_test,&
commander_symmetrize_map, commander_dock_volpair, commander_centervol, commander_reproject,&
commander_volanalyze, commander_volops

! mask commanders, masking and envelope maskign routines
use simple_commanders_mask, only: commander_auto_spher_mask, commander_mask, commander_automask2D

! validate commanders, for validation of parts of the pipelined stream process
use simple_commanders_validate, only: commander_mini_stream, commander_check_refpick

! oris commanders, construction and operations on oris (per-particle parameters)
use simple_commanders_oris, only: commander_make_oris, commander_orisops, commander_oristats, commander_vizoris

! check commanders, check number of and dimensions of images
use simple_commanders_checks, only: commander_info_image, commander_info_stktab

! misc commanders, miscallenous commanders (mostly printing)
use simple_commanders_misc, only: commander_print_fsc, commander_print_magic_boxes, commander_print_dose_weights

! sim commanders, for simulation of nosie, particls & movies
use simple_commanders_sim, only: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie

! atoms commanders, routines involving PDB files
use simple_commanders_atoms, only: commander_map2model_fsc, commander_pdb2mrc, commander_model_validation

! distr commanders, support routines for distributed execution
use simple_commanders_distr, only: commander_split

! test commanders, for testing purposes
use simple_commanders_test, only: commander_test_sim_workflow

end module simple_exec_api
