!> Aggregated public API for simple_exec.
!! 
!! This module centralizes all commander types used by the simple_exec
!! front-end. It:
!!  - USEs the various simple_commanders_* modules, and
!!  - re-exports the commander_* types via PUBLIC declarations.
!!
!! simple_exec should depend only on this module (plus simple_exec_helpers),
!! rather than on individual simple_commanders_* modules. To add a new
!! command:
!!  1. USE its defining simple_commanders_* module here.
!!  2. Add the corresponding commander_* symbol to the PUBLIC list.
!!  3. Declare a type(commander_...) variable in simple_exec and wire it
!!     into the SELECT CASE(prg) dispatch.
module simple_exec_module_api

! DEFINE API

! project commanders, operations on projects (sp_project) and associated files
use simple_commanders_project, only: commander_new_project, commander_update_project, commander_print_project_info,&
commander_print_project_field, commander_zero_project_shifts, commander_import_movies, commander_import_boxes,&
commander_import_particles, commander_import_cavgs, commander_replace_project_field, commander_selection,&
commander_merge_projects, commander_extract_subproj,&
commander_sample_classes, commander_write_mic_filetab, commander_prune_project_distr

! RELION commanders, export utility for RELION
use simple_commanders_relion, only: commander_export_relion

! starproject commanders, for star-file mirroring of SIMPLE project files
use simple_commanders_starproject, only: commander_import_starproject, commander_export_starproject,&
commander_assign_optics_groups

! preprocess commanders,  pre-processing routines
use simple_commanders_preprocess, only: commander_preprocess_distr, commander_motion_correct_distr,&
commander_gen_pspecs_and_thumbs_distr, commander_ctf_estimate_distr, commander_pick_distr,&
commander_extract_distr, commander_reextract_distr, commander_map_cavgs_selection

! cluster2D commanders, for simultanous 2D alignment and clustering of single-particle images
use simple_commanders_cluster2D, only: commander_make_cavgs_distr, commander_cluster2D_autoscale, commander_write_classes,&
commander_ppca_denoise_classes

! ab initio 2D commanders, for simultanous 2D alignment and clustering of single-particle images
use simple_commanders_abinitio2D, only: commander_abinitio2D

! stream_cluster2D commanders, for testing sieving in an offline fashion
use simple_stream_cluster2D_subsets, only: stream_cluster2D_subsets

! cavgs commanders, for operations on class averages
use simple_commanders_cavgs, only: commander_cluster_cavgs, commander_select_clusters, commander_match_cavgs

! imgproc commanders, standard image processing routines
use simple_commanders_imgproc, only: commander_cluster_stack, &
commander_match_stacks, commander_binarize, commander_convert, commander_ctfops, commander_ctf_phaseflip,&
commander_filter, commander_normalize, commander_ppca_denoise, commander_scale, commander_stack, commander_stackops

! resoltest commanders, for resolution estimation and regularization testing
use simple_commanders_resolest, only: commander_estimate_lpstages, commander_fsc, commander_clin_fsc,&
commander_uniform_filter2D, commander_uniform_filter3D, commander_icm2D, commander_icm3D

! volops commanders, operations on volumes
use simple_commanders_volops, only: commander_noisevol, commander_symaxis_search, commander_symmetry_test,&
commander_symmetrize_map, commander_dock_volpair, commander_postprocess, commander_centervol, commander_reproject,&
commander_volanalyze, commander_volops

! ab initio commanders, ab initio 3D reconstruction and multi-particle analysis workflows
use simple_commanders_abinitio, only: commander_abinitio3D_cavgs, commander_abinitio3D, commander_multivol_assign

! refine3D commanders, low-level methods for refine3D for ab initio 3D reconstruction and 3D refinement
use simple_commanders_refine3D, only: commander_refine3D_distr, commander_refine3D_auto

! rec commanders, 3D reconstrunction from aligned particles with Fourier griddding
use simple_commanders_rec, only: commander_reconstruct3D_distr

! mask commanders, masking and envelope maskign routines
use simple_commanders_mask, only: commander_automask, commander_auto_spher_mask, commander_mask, commander_automask2D

! validate commanders, for validation of parts of the pipelined stream process
use simple_commanders_validate, only: commander_mini_stream, commander_check_refpick

! oris commanders, construction and operations on oris (per-particle parameters)
use simple_commanders_oris, only: commander_make_oris, commander_orisops, commander_oristats, commander_vizoris

! check commanders, checck number of and dimensions of images
use simple_commanders_checks, only: commander_info_image, commander_info_stktab

! misc commanders, miscallenous commanders (mostly printing)
use simple_commanders_misc, only: commander_print_fsc, commander_print_magic_boxes, commander_print_dose_weights

! sim commanders, for simulation of nosie, particls & movies
use simple_commanders_sim, only: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie

! atoms commanders, routines involving PDB files
use simple_commanders_atoms, only: commander_map2model_fsc, commander_pdb2mrc

! distr commanders, support routines for distributed execution
use simple_commanders_distr, only: commander_split

! PUBLICIZE API

! project commanders, operations on projects (sp_project) and associated files
public :: commander_new_project, commander_update_project, commander_print_project_info,&
commander_print_project_field, commander_zero_project_shifts, commander_import_movies, commander_import_boxes,&
commander_import_particles, commander_import_cavgs, commander_replace_project_field, commander_selection,&
commander_merge_projects, commander_extract_subproj,&
commander_sample_classes, commander_write_mic_filetab, commander_prune_project_distr

! RELION commanders, export utility for RELION
public :: commander_export_relion

! starproject commanders, for star-file mirroring of SIMPLE project files
public :: commander_export_starproject, commander_assign_optics_groups

! preprocess commanders,  pre-processing routines
public :: commander_preprocess_distr, commander_motion_correct_distr,&
commander_gen_pspecs_and_thumbs_distr, commander_ctf_estimate_distr, commander_pick_distr,&
commander_extract_distr, commander_reextract_distr, commander_map_cavgs_selection 

! cluster2D commanders, for simultanous 2D alignment and clustering of single-particle images
public :: commander_make_cavgs_distr, commander_cluster2D_autoscale, commander_write_classes,&
commander_ppca_denoise_classes

! ab initio 2D commanders, for simultanous 2D alignment and clustering of single-particle images
public :: commander_abinitio2D

! stream_cluster2D commanders, for testing sieving in an offline fashion
public :: stream_cluster2D_subsets

! cavgs commanders, for operations on class averages
public :: commander_cluster_cavgs, commander_select_clusters, commander_match_cavgs

! imgproc commanders, standard image processing routines
public :: commander_cluster_stack,&
commander_match_stacks, commander_binarize, commander_convert, commander_ctfops, commander_ctf_phaseflip,&
commander_filter, commander_normalize, commander_ppca_denoise, commander_scale, commander_stack, commander_stackops

! resoltest commanders, for resolution estimation and regularization testing
public :: commander_estimate_lpstages, commander_fsc, commander_clin_fsc,&
commander_uniform_filter2D, commander_uniform_filter3D, commander_icm2D, commander_icm3D

! volops commanders, operations on volumes
public :: commander_noisevol, commander_symaxis_search, commander_symmetry_test,&
commander_symmetrize_map, commander_dock_volpair, commander_postprocess, commander_centervol, commander_reproject,&
commander_volanalyze, commander_volops

! ab initio commanders, ab initio 3D reconstruction and multi-particle analysis workflows
public :: commander_abinitio3D_cavgs, commander_abinitio3D, commander_multivol_assign

! refine3D commanders, low-level methods for refine3D for ab initio 3D reconstruction and 3D refinement
public :: commander_refine3D_distr, commander_refine3D_auto

! rec commanders, 3D reconstrunction from aligned particles with Fourier griddding
public :: commander_reconstruct3D_distr

! mask commanders, masking and envelope maskign routines
public :: commander_automask, commander_auto_spher_mask, commander_mask, commander_automask2D

! validate commanders, for validation of parts of the pipelined stream process
public :: commander_mini_stream, commander_check_refpick

! oris commanders, construction and operations on oris (per-particle parameters)
public :: commander_make_oris, commander_orisops, commander_oristats, commander_vizoris

! check commanders, checck number of and dimensions of images
public :: commander_info_image, commander_info_stktab

! misc commanders, miscallenous commanders (mostly printing)
public :: commander_print_fsc, commander_print_magic_boxes, commander_print_dose_weights

! sim commanders, for simulation of nosie, particls & movies
public :: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie

! atoms commanders, routines involving PDB files
public :: commander_map2model_fsc, commander_pdb2mrc

! distr commanders, support routines for distributed execution
public :: commander_split

end module simple_exec_module_api
