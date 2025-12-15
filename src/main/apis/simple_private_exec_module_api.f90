!> Aggregated public API for simple_private_exec.
!! 
!! This module centralizes all commander types used by the simple_private_exec
!! front-end. It:
!!  - USEs the various simple_commanders_* modules
!!
!! simple_private_exec should depend only on this module (plus simple_exec_helpers),
!! rather than on individual simple_commanders_* modules. To add a new
!! command:
!!  1. USE its defining simple_commanders_* module here.
!!  2. Add the corresponding commander_* symbol to the PUBLIC list.
!!  3. Declare a type(commander_...) variable in simple_exec and wire it
!!     into the SELECT CASE(prg) dispatch.
module simple_private_exec_module_api
use simple_core_module_api
use simple_cmdline,        only: cmdline, cmdline_err
use simple_jiffys,         only: simple_print_timer
use simple_private_prgs,   only: make_private_user_interface
use simple_symanalyzer,    only: print_subgroups
use simple_syslib,         only: print_slurm_env
use simple_user_interface, only: make_user_interface, print_ui_json, write_ui_json, print_stream_ui_json

! pick commanders, picking routines
use simple_commanders_pick, only: commander_extract, commander_reextract, commander_pick_extract,&
commander_pick, commander_make_pickrefs

! preprocess commanders,  pre-processing routines
use simple_commanders_preprocess, only: commander_preprocess, commander_motion_correct,&
commander_gen_pspecs_and_thumbs, commander_ctf_estimate

! cluster2D commanders, for simultanous 2D alignment and clustering of single-particle images
use simple_commanders_cluster2D, only: commander_make_cavgs, commander_cluster2D, commander_cluster2D_distr,&
commander_cavgassemble, commander_prob_tab2D

! cavgs commanders, for operations on class averages
use simple_commanders_cavgs, only: commander_rank_cavgs, commander_shape_rank_cavgs

! project commanders, operations on projects (sp_project) and associated files
use simple_commanders_project, only: commander_export_cavgs, commander_print_project_vals, commander_prune_project,&
commander_scale_project_distr

! refine3D commanders, low-level methods for refine3D for ab initio 3D reconstruction and 3D refinement
use simple_commanders_refine3D, only: commander_refine3D, commander_check_3Dconv, commander_prob_tab,&
commander_prob_align

! euclid commanders, for obtaining signal statistics for noise normalized Euclidean distance functiuon evalutaion
use simple_commanders_euclid, only: commander_calc_pspec_distr, commander_calc_pspec, commander_calc_pspec_assemble,&
commander_calc_group_sigmas

! rec commanders, 3D reconstrunction from aligned particles with Fourier griddding
use simple_commanders_rec, only: commander_volassemble, commander_reconstruct3D

! check commanders, check number of and dimensions of images
use simple_commanders_checks, only: commander_check_box, commander_check_nptcls, commander_check_stoch_update,&
commander_check_update_frac

! volops commanders, operations on volumes
use simple_commanders_volops, only: commander_postprocess

! mask commanders, masking and envelope maskign routines
use simple_commanders_mask, only: commander_automask

! imgops commanders, standard image operation routines
use simple_commanders_imgops, only: commander_scale, commander_binarize

! misc commanders, miscallenous commanders (mostly printing)
use simple_commanders_misc, only: commander_kstest, commander_pearsn

! oris commanders, construction and operations on oris (per-particle parameters)
use simple_commanders_oris, only: commander_rotmats2oris

! tseries commanders, methods operating on time-series data obtained with GLC-EM
use simple_commanders_tseries, only: commander_tseries_track_particles, commander_tseries_motion_correct

! distr commanders, support routines for distributed execution
use simple_commanders_distr, only: commander_split
end module simple_private_exec_module_api
