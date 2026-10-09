!@descr: the unit-test suites: the build's fast gate (test=unit_<area>), its umbrella (test=units) and the platform-tier forked-process suite
module simple_commanders_test_class
use simple_commanders_api
use simple_test_utils,                       only: begin_test_suite, end_test_suite, reset_test_report, report_summary, &
    &set_fixed_seed
! core library tester modules
use simple_ansi_ctrls_tester,                only: run_all_ansi_ctrls_tests
use simple_string_tester,                    only: run_all_string_tests
use simple_syslib_tester,                    only: run_all_syslib_tests
use simple_fileio_tester,                    only: run_all_fileio_tests
use simple_stack_io_tester,                  only: run_all_stack_io_tests
use simple_class_sample_io_tester,           only: run_all_class_sample_io_tests
use simple_chash_tester,                     only: run_all_chash_tests
use simple_vrefhash_tester,                  only: run_all_vrefhash_tests
use simple_hash_tester,                      only: run_all_hash_tests
use simple_linked_list_tester,               only: run_all_list_tests
use simple_cmdline_tester,                   only: run_all_cmdline_tests
use simple_rec_list_tester,                  only: run_all_rec_list_tests
use simple_ori_tester,                       only: run_all_ori_tests
use simple_oris_tester,                      only: run_all_oris_tests
use simple_sym_tester,                       only: run_all_sym_tests
use simple_stat_tester,                      only: run_all_stat_tests
use simple_linalg_tester,                    only: run_all_linalg_tests
use simple_kbinterpol_tester,                only: run_all_kbinterpol_tests
use simple_srch_sort_loc_tester,             only: run_all_srch_sort_loc_tests
use simple_decay_funs_tester,                only: run_all_decay_funs_tests
use simple_pca_tester,                       only: run_all_pca_tests
use simple_opt_tester,                       only: run_all_opt_tests
use simple_lpstages_tester,                  only: run_all_lpstages_tests
use simple_image_msk_tester,                 only: run_all_mask_tests, run_all_nano_mask_tests, run_all_image_bin_tests
use simple_segmentation_tester,              only: run_all_segmentation_tests
use simple_accum_blend_tester,               only: run_all_accum_blend_tests
use simple_ctf_tester,                       only: run_all_ctf_tests
use simple_starfile_tester,                  only: run_all_starfile_tests
use simple_starproject_tester,               only: run_all_starproject_tests
use simple_starproject_stream_tester,        only: run_all_starproject_stream_tests
use simple_binoris_tester,                   only: run_all_binoris_tests
use simple_sp_project_tester,                only: run_all_sp_project_tests
use simple_project_merge_tester,             only: run_all_project_merge_tests
use simple_state_weight_set_tester,          only: run_all_state_weight_set_tests
use simple_class_compatibility_tester,       only: run_all_class_compatibility_tests
use simple_ptcl_sieve_tester,                only: run_all_ptcl_sieve_tests
use simple_mic_import_tester,                only: run_all_mic_import_tests
use simple_mic_selection_tester,             only: run_all_mic_selection_tests
use simple_optics_groups_tester,             only: run_all_optics_groups_tests
use simple_optics_maps_tester,               only: run_all_optics_maps_tests
use simple_motion_gain_tester,               only: run_all_motion_gain_tests
use simple_motion_model_tester,              only: run_all_motion_model_tests
use simple_particle_extractor_tester,        only: run_all_particle_extractor_tests
use simple_gui_metadata_tester,              only: run_all_gui_metadata_tests
use simple_gui_assembler_tester,             only: run_all_gui_assembler_tests
use simple_ui_hash_tester,                   only: run_all_ui_hash_tests
use simple_ui_visibility_tester,             only: run_all_ui_visibility_tests
use simple_rnd_tester,                       only: run_all_rnd_tests
use simple_diff_map_graphs_tester,           only: run_all_diff_map_graphs_tests
use simple_qsys_ctrl_tester,                 only: run_all_qsys_ctrl_tests
use simple_qsys_env_tester,                  only: run_all_qsys_env_tests
use simple_cavg_quality_relations_tester,    only: run_all_cavg_quality_relations_tests
use simple_sigma2_state_tester,              only: run_all_sigma2_state_tests
use simple_sigma2_bootstrap_tester,          only: run_all_sigma2_bootstrap_tests
use simple_eul_prob_tab2D_tester,            only: run_all_eul_prob_tab2D_tests
use simple_cavg_sums_tester,                 only: run_all_cavg_sums_tests
use simple_trail_chain_manifest_tester,      only: run_all_trail_chain_manifest_tests
use simple_gauran_tester,                    only: run_all_gauran_tests
use simple_rec3D_strategy_tester,            only: run_all_rec3D_strategy_tests
use simple_frozen_accum_tester,              only: run_all_frozen_accum_tests
use simple_reconstructor_tester,             only: run_all_reconstructor_tests
use simple_rec3D_service_tester,             only: run_all_rec3D_service_lib_tests
use simple_solve3D_manifest_tester,          only: run_all_solve3D_manifest_tests
use simple_project_superset_tester,          only: run_all_project_superset_tests
use simple_solve3D_addon_report_tester,      only: run_all_solve3D_addon_report_tests
use simple_volpair_metrics_tester,           only: run_all_volpair_metrics_tests
use simple_nu_filter_tester,                 only: run_all_nu_filter_tests
use simple_solve3D_addon_dock_tester,        only: run_all_solve3D_addon_dock_tests
use simple_pcg_halfset_tester,               only: run_all_pcg_halfset_tests
use simple_pftc_inplane_tester,              only: run_all_pftc_inplane_tests
use simple_strategy3D_inplane_tester,        only: run_all_strategy3D_inplane_tests
use simple_strategy3D_cont_tester,           only: run_all_strategy3D_cont_tests
use simple_convergence_tester,               only: run_all_convergence_tests
use simple_strategy3D_cont_1jyx_tester,      only: run_all_strategy3D_cont_1jyx_tests
use simple_cartesian_fourier_tester,         only: run_all_cartesian_fourier_tests
use simple_cartft_calc_tester,               only: run_all_cartft_calc_tests
use simple_cartft_pose_opt_tester,           only: run_all_cartft_pose_opt_tests
use simple_flex_pca_tester,                  only: run_all_flex_pca_tests, run_all_flex_pca_lib_tests
use simple_flex_pca_application_tester,      only: run_all_flex_pca_application_tests
use simple_flex_pcg_tester,                  only: run_all_flex_pcg_tests, run_all_flex_pcg_lib_tests, &
    &run_all_flex_pcg_sweep_tests
use simple_ipc_tcp_socket_tester,            only: run_all_ipc_tcp_socket_tests
use simple_http_post_tester,                 only: run_all_http_post_tests
use simple_persistent_worker_message_tester, only: run_all_persistent_worker_message_tests
use simple_persistent_worker_server_tester,  only: run_all_persistent_worker_server_tests
use simple_stream_pipe_tester,               only: run_all_stream_pipe_tests
use simple_forked_process_tester,            only: run_all_forked_process_tests
use simple_imghead_tester,                   only: run_all_imghead_tests
use simple_image_tester,                     only: run_all_image_tests
use simple_jpg_tester,                       only: run_all_mrc2jpeg_tests
use simple_mrc_validate_tester,              only: run_all_mrc_validate_tests
use simple_volume_shape_tester,              only: run_all_volume_shape_tests
use simple_ftiter_tester,                    only: run_all_ftiter_tests
use simple_ftexp_shsrch_tester,              only: run_all_ftexp_shsrch_tests
use simple_bspline_smoother_tester,          only: run_all_bspline_smoother_tests
use simple_online_var_tester,                only: run_all_online_var_tests
use simple_aff_prop_tester,                  only: run_all_aff_prop_tests
use simple_hac_tester,                       only: run_all_hac_tests
use simple_kmeans_tester,                    only: run_all_kmeans_tests
use simple_kcenter_tester,                   only: run_all_kcenter_tests
use simple_gmm_tester,                       only: run_all_gmm_tests
use simple_xd_gmm_tester,                    only: run_all_xd_gmm_tests
use simple_atoms_tester,                     only: run_all_atoms_tests
use simple_cif2mrc_tester,                   only: run_all_cif2mrc_tests
use simple_calpha_finder_tester,             only: run_all_calpha_finder_tests
use simple_nano_species_tester,              only: run_all_nano_species_tests
use simple_pdb2mrc_tester,                   only: run_all_pdb2mrc_tests
use simple_image_serialize_tester,           only: run_all_image_serialize_tests
use simple_cavg_registration_tester,         only: run_all_cavg_registration_tests
use simple_polarft_corr_tester,              only: run_all_polarft_corr_tests
use simple_openmp_offload_tester,            only: run_openmp_offload_tests
use simple_stream_tester,                    only: run_all_stream_optics_tests, run_all_stream_pickrefs_tests, &
    &run_all_stream_pick_extract_tests
use simple_stream_chain_tester,              only: run_all_stream_chain_sieve_tests, run_all_stream_chain_movie_tests
use simple_stream_stage_preprocess_tester,   only: run_all_stream_stage_preprocess_tests
use simple_stream_stage_optics_tester,       only: run_all_stream_stage_optics_tests
use simple_stream_stage_initial_analysis_tester, only: run_all_stream_stage_initial_analysis_tests
use simple_stream_stage_refpick_tester,      only: run_all_stream_stage_refpick_tests
use simple_stream_stage_sieve_tester,        only: run_all_stream_stage_sieve_tests
use simple_stream_stage_pool2D_tester,       only: run_all_stream_stage_pool2D_tests
use simple_stream_pool2D_tester,             only: run_all_stream_pool2D_tests
use simple_stream_stage_solve3D_tester,   only: run_all_stream_stage_solve3D_tests
use simple_stream_master_tester,             only: run_all_stream_master_tests, run_stream_heartbeat_tests
use simple_stream_job_sets_tester,           only: run_all_stream_job_sets_tests
use simple_stream_meta_plots_tester,         only: run_all_stream_meta_plots_tests
use simple_stream_watcher_tester,            only: run_all_stream_watcher_tests
use simple_ui,                               only: validate_ui_json
implicit none
#include "simple_local_flags.inc"

! Fast gate: the fifteen test=unit_<area> suites, one CTest entry each, label `fast`
! (doc/refactoring_notes/completed/uniform_test_environment_refactoring.md, section 5.1): assertions
! through simple_test_utils, localhost only, no downloads or user data, one OpenMP thread.
! suite=<name> runs one sub-suite (suite=list names them); test=units runs every area (not the gate);
! lib_<area> are nightly library suites (section 5.2.1); forked_process and openmp_offload are label
! `platform`. SIMPLE_UNIT_ORDER=reverse runs a table backwards: a differing result is a state leak.


type, extends(commander_base) :: commander_test_units
  contains
    procedure :: execute      => exec_test_units
end type commander_test_units

type, extends(commander_base) :: commander_test_unit_core
  contains
    procedure :: execute      => exec_test_unit_core
end type commander_test_unit_core

type, extends(commander_base) :: commander_test_unit_ori
  contains
    procedure :: execute      => exec_test_unit_ori
end type commander_test_unit_ori

type, extends(commander_base) :: commander_test_unit_image
  contains
    procedure :: execute      => exec_test_unit_image
end type commander_test_unit_image

type, extends(commander_base) :: commander_test_unit_numerics
  contains
    procedure :: execute      => exec_test_unit_numerics
end type commander_test_unit_numerics

type, extends(commander_base) :: commander_test_unit_project
  contains
    procedure :: execute      => exec_test_unit_project
end type commander_test_unit_project

type, extends(commander_base) :: commander_test_unit_motion
  contains
    procedure :: execute      => exec_test_unit_motion
end type commander_test_unit_motion

type, extends(commander_base) :: commander_test_unit_ui
  contains
    procedure :: execute      => exec_test_unit_ui
end type commander_test_unit_ui

type, extends(commander_base) :: commander_test_unit_ipc
  contains
    procedure :: execute      => exec_test_unit_ipc
end type commander_test_unit_ipc

type, extends(commander_base) :: commander_test_unit_reconstruction
  contains
    procedure :: execute      => exec_test_unit_reconstruction
end type commander_test_unit_reconstruction

type, extends(commander_base) :: commander_test_lib_reconstruction
  contains
    procedure :: execute      => exec_test_lib_reconstruction
end type commander_test_lib_reconstruction

type, extends(commander_base) :: commander_test_unit_pftc_align2D3D
  contains
    procedure :: execute      => exec_test_unit_pftc_align2D3D
end type commander_test_unit_pftc_align2D3D

type, extends(commander_base) :: commander_test_unit_cart_align3D
  contains
    procedure :: execute      => exec_test_unit_cart_align3D
end type commander_test_unit_cart_align3D

type, extends(commander_base) :: commander_test_lib_cart_align3D
  contains
    procedure :: execute      => exec_test_lib_cart_align3D
end type commander_test_lib_cart_align3D

type, extends(commander_base) :: commander_test_unit_single
  contains
    procedure :: execute      => exec_test_unit_single
end type commander_test_unit_single

type, extends(commander_base) :: commander_test_lib_single
  contains
    procedure :: execute      => exec_test_lib_single
end type commander_test_lib_single

type, extends(commander_base) :: commander_test_unit_stream
  contains
    procedure :: execute      => exec_test_unit_stream
end type commander_test_unit_stream

type, extends(commander_base) :: commander_test_lib_stream
  contains
    procedure :: execute      => exec_test_lib_stream
end type commander_test_lib_stream

type, extends(commander_base) :: commander_test_unit_parallel
  contains
    procedure :: execute      => exec_test_unit_parallel
end type commander_test_unit_parallel

type, extends(commander_base) :: commander_test_unit_heterogeneity
  contains
    procedure :: execute      => exec_test_unit_heterogeneity
end type commander_test_unit_heterogeneity

type, extends(commander_base) :: commander_test_lib_heterogeneity
  contains
    procedure :: execute      => exec_test_lib_heterogeneity
end type commander_test_lib_heterogeneity

type, extends(commander_base) :: commander_test_openmp_offload
  contains
    procedure :: execute      => exec_test_openmp_offload
end type commander_test_openmp_offload

type, extends(commander_base) :: commander_test_forked_process
  contains
    procedure :: execute      => exec_test_forked_process
end type commander_test_forked_process

! a sub-suite: a name and a no-argument procedure that asserts through simple_test_utils
abstract interface
    subroutine no_arg_test()
    end subroutine no_arg_test
end interface

type :: unit_suite
    character(len=32) :: name = ''
    procedure(no_arg_test), pointer, nopass :: run => null()
end type unit_suite

integer, parameter :: MAX_SUITES = 128   ! must hold every sub-suite of `units`; add_suite stops on overflow

contains

    ! ---- area tables ---------------------------------------------------------
    ! One function per area returns its sub-suites in table order. The umbrella
    ! (test=units) is the concatenation of all thirteen.

    subroutine suites_core( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'ANSI formatting',      run_all_ansi_ctrls_tests)
        call add_suite(s, n, 'string',               run_all_string_tests)
        call add_suite(s, n, 'syslib',               run_all_syslib_tests)
        call add_suite(s, n, 'fileio',               run_all_fileio_tests)
        call add_suite(s, n, 'stack I/O',            run_all_stack_io_tests)
        call add_suite(s, n, 'class sample I/O',     run_all_class_sample_io_tests)
        call add_suite(s, n, 'character hash',       run_all_chash_tests)
        call add_suite(s, n, 'hash',                 run_all_hash_tests)
        call add_suite(s, n, 'value-reference hash', run_all_vrefhash_tests)
        call add_suite(s, n, 'linked list',          run_all_list_tests)
        call add_suite(s, n, 'record list',          run_all_rec_list_tests)
        call add_suite(s, n, 'command line',         run_all_cmdline_tests)
    end subroutine suites_core

    subroutine suites_ori( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'orientation',            run_all_ori_tests)
        call add_suite(s, n, 'orientation collection', run_all_oris_tests)
        call add_suite(s, n, 'symmetry',               run_all_sym_tests)
        call add_suite(s, n, 'Euler shift',            test_euler_shift)
    end subroutine suites_ori

    subroutine suites_image( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'image',                run_all_image_tests)
        call add_suite(s, n, 'mrc2jpeg',             run_all_mrc2jpeg_tests)
        call add_suite(s, n, 'mrc validate',          run_all_mrc_validate_tests)
        call add_suite(s, n, 'image header',         run_all_imghead_tests)
        call add_suite(s, n, 'Fourier iterator',     run_all_ftiter_tests)
        call add_suite(s, n, 'B-spline smoother',    run_all_bspline_smoother_tests)
        call add_suite(s, n, 'masks',                run_all_mask_tests)
        call add_suite(s, n, 'nano mask',            run_all_nano_mask_tests)
        call add_suite(s, n, 'volume shape',         run_all_volume_shape_tests)
        call add_suite(s, n, 'binary image',         run_all_image_bin_tests)
        call add_suite(s, n, 'segmentation',         run_all_segmentation_tests)
        call add_suite(s, n, 'trailing-reconstruction blend', run_all_accum_blend_tests)
        call add_suite(s, n, 'CTF',                  run_all_ctf_tests)
        call add_suite(s, n, 'image serialisation',  run_all_image_serialize_tests)
    end subroutine suites_image

    subroutine suites_numerics( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'online variance',         run_all_online_var_tests)
        call add_suite(s, n, 'random draws',            run_all_rnd_tests)
        call add_suite(s, n, 'affinity propagation',    run_all_aff_prop_tests)
        call add_suite(s, n, 'hierarchical clustering', run_all_hac_tests)
        call add_suite(s, n, 'k-means',                 run_all_kmeans_tests)
        call add_suite(s, n, 'k-center',                run_all_kcenter_tests)
        call add_suite(s, n, 'Gaussian mixture',        run_all_gmm_tests)
        call add_suite(s, n, 'extreme deconvolution',   run_all_xd_gmm_tests)
        call add_suite(s, n, 'statistics',              run_all_stat_tests)
        call add_suite(s, n, 'linear algebra',          run_all_linalg_tests)
        call add_suite(s, n, 'Kaiser-Bessel kernel',    run_all_kbinterpol_tests)
        call add_suite(s, n, 'search, sort, locate',    run_all_srch_sort_loc_tests)
        call add_suite(s, n, 'decay schedules',         run_all_decay_funs_tests)
        call add_suite(s, n, 'PCA',                     run_all_pca_tests)
        call add_suite(s, n, 'cavg quality relations',  run_all_cavg_quality_relations_tests)
        call add_suite(s, n, 'diffusion-map graphs',    run_all_diff_map_graphs_tests)
        call add_suite(s, n, 'optimisers',              run_all_opt_tests)
        call add_suite(s, n, 'low-pass stages',         run_all_lpstages_tests)
        ! motion-correction shift search on expanded Fourier transforms (an optimiser, not an image test)
        call add_suite(s, n, 'shift search',            run_all_ftexp_shsrch_tests)
    end subroutine suites_numerics

    subroutine suites_project( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'STAR file',               run_all_starfile_tests)
        call add_suite(s, n, 'STAR project',            run_all_starproject_tests)
        call add_suite(s, n, 'STAR stream export',      run_all_starproject_stream_tests)
        call add_suite(s, n, 'binoris',                 run_all_binoris_tests)
        call add_suite(s, n, 'project records',         run_all_sp_project_tests)
        call add_suite(s, n, 'project merge',           run_all_project_merge_tests)
        call add_suite(s, n, 'class compatibility',     run_all_class_compatibility_tests)
        call add_suite(s, n, 'particle sieve',          run_all_ptcl_sieve_tests)
        call add_suite(s, n, 'micrograph import',       run_all_mic_import_tests)
        call add_suite(s, n, 'micrograph selection',    run_all_mic_selection_tests)
        call add_suite(s, n, 'optics groups',           run_all_optics_groups_tests)
        call add_suite(s, n, 'optics maps',             run_all_optics_maps_tests)
        call add_suite(s, n, 'solve3D manifest',        run_all_solve3D_manifest_tests)
        call add_suite(s, n, 'project superset',        run_all_project_superset_tests)
        call add_suite(s, n, 'solve3D addon report', run_all_solve3D_addon_report_tests)
        call add_suite(s, n, 'state weight set',        run_all_state_weight_set_tests)
    end subroutine suites_project

    subroutine suites_motion( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'motion gain',             run_all_motion_gain_tests)
        call add_suite(s, n, 'motion model',            run_all_motion_model_tests)
        call add_suite(s, n, 'particle extractor',      run_all_particle_extractor_tests)
    end subroutine suites_motion

    subroutine suites_ui( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'UI JSON',       suite_ui_json)
        call add_suite(s, n, 'GUI metadata',  run_all_gui_metadata_tests)
        call add_suite(s, n, 'GUI assembler', run_all_gui_assembler_tests)
        call add_suite(s, n, 'UI hash',       run_all_ui_hash_tests)
        call add_suite(s, n, 'UI visibility', run_all_ui_visibility_tests)
    end subroutine suites_ui

    subroutine suites_ipc( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        ! localhost only, bounded; `forked process` is deliberately not here
        call add_suite(s, n, 'IPC TCP socket',            run_all_ipc_tcp_socket_tests)
        call add_suite(s, n, 'HTTP POST',                 run_all_http_post_tests)
        call add_suite(s, n, 'persistent worker message', run_all_persistent_worker_message_tests)
        call add_suite(s, n, 'persistent worker server',  run_all_persistent_worker_server_tests)
        call add_suite(s, n, 'stream pipe',               run_all_stream_pipe_tests)
    end subroutine suites_ipc

    subroutine suites_reconstruction( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'rec3D backend',             run_all_rec3D_strategy_tests)
        call add_suite(s, n, 'observation noise',         run_all_gauran_tests)
        call add_suite(s, n, 'class-average carry-over',  run_all_cavg_sums_tests)
        call add_suite(s, n, 'trailing chain identity',   run_all_trail_chain_manifest_tests)
        call add_suite(s, n, 'frozen accumulator',        run_all_frozen_accum_tests)
        call add_suite(s, n, 'weighted insertion',        run_all_reconstructor_tests)
        call add_suite(s, n, 'volume pair metrics',       run_all_volpair_metrics_tests)
        call add_suite(s, n, 'nonuniform filtering',      run_all_nu_filter_tests)
    end subroutine suites_reconstruction

    !> registration on the polar Fourier transform, shared by the 2D and 3D searches
    subroutine suites_pftc_align2D3D( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'polar correlation',        run_all_polarft_corr_tests)
        call add_suite(s, n, 'continuous in-plane',      run_all_pftc_inplane_tests)
        call add_suite(s, n, 'refine3D in-plane state',  run_all_strategy3D_inplane_tests)
        call add_suite(s, n, '2D probability table I/O', run_all_eul_prob_tab2D_tests)
        call add_suite(s, n, 'sigma2 state',             run_all_sigma2_state_tests)
        call add_suite(s, n, 'sigma2 bootstrap',         run_all_sigma2_bootstrap_tests)
        call add_suite(s, n, 'cavg registration',          run_all_cavg_registration_tests)
    end subroutine suites_pftc_align2D3D

    !> the Cartesian (continuous) 3D registration: the neutral Cartesian Fourier layer (shared
    !! with PCG), the Cartesian calculator and pose optimizer, the Cartesian 3D strategy, and
    !! the pose statistics of the convergence report
    subroutine suites_cart_align3D( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'Cartesian Fourier', run_all_cartesian_fourier_tests)
        call add_suite(s, n, 'cart calculator',   run_all_cartft_calc_tests)
        call add_suite(s, n, 'pose optimizer',    run_all_cartft_pose_opt_tests)
        call add_suite(s, n, 'pose strategy',     run_all_strategy3D_cont_tests)
        call add_suite(s, n, 'pose statistics',   run_all_convergence_tests)
    end subroutine suites_cart_align3D

    !> nightly: 1 000 simulated 1JYX particles through the Cartesian 3D strategy under both
    !! objectives at lp 8 A and 4 A, minutes
    subroutine suites_lib_cart_align3D( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'pose 1JYX recovery', run_all_strategy3D_cont_1jyx_tests)
    end subroutine suites_lib_cart_align3D

    !> SINGLE (nanoparticles, atomic models): the atoms module, the C-alpha candidate search and the
    !! species-discovery numerics
    subroutine suites_single( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'atoms',          run_all_atoms_tests)
        call add_suite(s, n, 'cif2mrc',        run_all_cif2mrc_tests)
        call add_suite(s, n, 'C-alpha finder', run_all_calpha_finder_tests)
        call add_suite(s, n, 'species',        run_all_nano_species_tests)
    end subroutine suites_single

    !> nightly pdb2mrc coverage of the built-in 6VXX and 1JYX models
    subroutine suites_lib_single( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'pdb2mrc', run_all_pdb2mrc_tests)
    end subroutine suites_lib_single

    !> STREAM stages (src/main/stream/stages), step by step without a queue or waits:
    !! preprocessing, optics assignment, initial analysis, reference picking, particle sieving,
    !! pool 2D, solve 3D, the job sets preprocessing and reference picking run on, the GUI plots
    !! preprocessing fills, the master's GUI commands, metadata store and stage pipes, and the
    !! watcher's history
    subroutine suites_stream( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'preprocessing',      run_all_stream_stage_preprocess_tests)
        call add_suite(s, n, 'optics assignment',  run_all_stream_stage_optics_tests)
        call add_suite(s, n, 'initial analysis',   run_all_stream_stage_initial_analysis_tests)
        call add_suite(s, n, 'reference picking',  run_all_stream_stage_refpick_tests)
        call add_suite(s, n, 'particle sieving',   run_all_stream_stage_sieve_tests)
        call add_suite(s, n, 'pool 2D',            run_all_stream_stage_pool2D_tests)
        call add_suite(s, n, 'pool 2D object',     run_all_stream_pool2D_tests)
        call add_suite(s, n, 'solve 3D',        run_all_stream_stage_solve3D_tests)
        call add_suite(s, n, 'stream master',      run_all_stream_master_tests)
        call add_suite(s, n, 'job sets',           run_all_stream_job_sets_tests)
        call add_suite(s, n, 'meta plots',         run_all_stream_meta_plots_tests)
        call add_suite(s, n, 'stream watcher',     run_all_stream_watcher_tests)
    end subroutine suites_stream

    !> nightly: the stream stages that run in-process, with the arguments the stream gives them
    !! (Ruben's stream tests; optics assignment waits a minute for the stream watcher);
    !! doc/refactoring_notes/planned/stream_area_tests_handover.md says what they should pin beyond counts
    subroutine suites_lib_stream( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'optics assignment',  run_all_stream_optics_tests)
        call add_suite(s, n, 'picking references', run_all_stream_pickrefs_tests)
        call add_suite(s, n, 'pick and extract',   run_all_stream_pick_extract_tests)
        ! decision 19 of the 5 October 2026 fix plan: the stages chained without the master; they
        ! submit jobs to the local queue (installed executables) and take tens of minutes
        call add_suite(s, n, 'sieve to 3D',        run_all_stream_chain_sieve_tests)
        call add_suite(s, n, 'movies to 3D',       run_all_stream_chain_movie_tests)
    end subroutine suites_lib_stream

    !> distributed execution: the job controller and the queue-system environment
    subroutine suites_parallel( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'qsys control',     run_all_qsys_ctrl_tests)
        call add_suite(s, n, 'qsys environment', run_all_qsys_env_tests)
    end subroutine suites_parallel

    !> heterogeneity analysis (flex_pca): latent model, state weights, deconvolution and
    !! the PCG M-step operator
    subroutine suites_heterogeneity( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'flex PCA',          run_all_flex_pca_tests)
        call add_suite(s, n, 'flex PCG operator', run_all_flex_pcg_tests)
    end subroutine suites_heterogeneity

    !> nightly: application-level two-state phantom, deconvolution of 20000 particles at realistic
    !! noise, the PCG M-step operator at box 64, and the twelve-setting PCG solve sweep at box 32
    subroutine suites_lib_heterogeneity( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'flex PCA two-state phantom', run_all_flex_pca_application_tests)
        call add_suite(s, n, 'flex PCA deconvolution 20k', run_all_flex_pca_lib_tests)
        call add_suite(s, n, 'flex PCG operator 64',       run_all_flex_pcg_lib_tests)
        call add_suite(s, n, 'flex PCG solve sweep',       run_all_flex_pcg_sweep_tests)
    end subroutine suites_lib_heterogeneity

    !> nightly library suite: minutes, full boxes allowed, same assertions and runner
    subroutine suites_lib_reconstruction( s, n )
        type(unit_suite), intent(inout) :: s(:)
        integer,          intent(inout) :: n
        call add_suite(s, n, 'PCG half-set',         run_all_pcg_halfset_tests)
        call add_suite(s, n, 'addon report docking', run_all_solve3D_addon_dock_tests)
        call add_suite(s, n, 'fractional reconstruction', run_all_rec3D_service_lib_tests)
    end subroutine suites_lib_reconstruction

    subroutine add_suite( s, n, name, proc )
        type(unit_suite),      intent(inout) :: s(:)
        integer,               intent(inout) :: n
        character(len=*),      intent(in)    :: name
        procedure(no_arg_test)               :: proc
        if( n >= size(s) ) THROW_HARD('too many unit sub-suites; raise MAX_SUITES')
        n = n + 1
        s(n)%name = name
        s(n)%run  => proc
    end subroutine add_suite

    ! ---- commanders ----------------------------------------------------------

    subroutine exec_test_units( self, cline )
        class(commander_test_units), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_core(s, n)
        call suites_ori(s, n)
        call suites_image(s, n)
        call suites_numerics(s, n)
        call suites_project(s, n)
        call suites_motion(s, n)
        call suites_ui(s, n)
        call suites_ipc(s, n)
        call suites_reconstruction(s, n)
        call suites_pftc_align2D3D(s, n)
        call suites_cart_align3D(s, n)
        call suites_heterogeneity(s, n)
        call suites_parallel(s, n)
        call suites_single(s, n)
        call suites_stream(s, n)
        call run_unit_suites('units', cline, s(1:n))
    end subroutine exec_test_units

    subroutine exec_test_unit_core( self, cline )
        class(commander_test_unit_core), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_core(s, n)
        call run_unit_suites('unit_core', cline, s(1:n))
    end subroutine exec_test_unit_core

    subroutine exec_test_unit_ori( self, cline )
        class(commander_test_unit_ori), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ori(s, n)
        call run_unit_suites('unit_ori', cline, s(1:n))
    end subroutine exec_test_unit_ori

    subroutine exec_test_unit_image( self, cline )
        class(commander_test_unit_image), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_image(s, n)
        call run_unit_suites('unit_image', cline, s(1:n))
    end subroutine exec_test_unit_image

    subroutine exec_test_unit_numerics( self, cline )
        class(commander_test_unit_numerics), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_numerics(s, n)
        call run_unit_suites('unit_numerics', cline, s(1:n))
    end subroutine exec_test_unit_numerics

    subroutine exec_test_unit_project( self, cline )
        class(commander_test_unit_project), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_project(s, n)
        call run_unit_suites('unit_project', cline, s(1:n))
    end subroutine exec_test_unit_project

    subroutine exec_test_unit_motion( self, cline )
        class(commander_test_unit_motion), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_motion(s, n)
        call run_unit_suites('unit_motion', cline, s(1:n))
    end subroutine exec_test_unit_motion

    subroutine exec_test_unit_ui( self, cline )
        class(commander_test_unit_ui), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ui(s, n)
        call run_unit_suites('unit_ui', cline, s(1:n))
    end subroutine exec_test_unit_ui

    subroutine exec_test_unit_ipc( self, cline )
        class(commander_test_unit_ipc), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_ipc(s, n)
        call run_unit_suites('unit_ipc', cline, s(1:n))
    end subroutine exec_test_unit_ipc

    subroutine exec_test_unit_reconstruction( self, cline )
        class(commander_test_unit_reconstruction), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_reconstruction(s, n)
        call run_unit_suites('unit_reconstruction', cline, s(1:n))
    end subroutine exec_test_unit_reconstruction

    subroutine exec_test_unit_pftc_align2D3D( self, cline )
        class(commander_test_unit_pftc_align2D3D), intent(inout) :: self
        class(cmdline),                                   intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_pftc_align2D3D(s, n)
        call run_unit_suites('unit_pftc_align2D3D', cline, s(1:n))
    end subroutine exec_test_unit_pftc_align2D3D

    subroutine exec_test_unit_cart_align3D( self, cline )
        class(commander_test_unit_cart_align3D), intent(inout) :: self
        class(cmdline),                                 intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_cart_align3D(s, n)
        call run_unit_suites('unit_cart_align3D', cline, s(1:n))
    end subroutine exec_test_unit_cart_align3D

    subroutine exec_test_lib_cart_align3D( self, cline )
        class(commander_test_lib_cart_align3D), intent(inout) :: self
        class(cmdline),                                intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_cart_align3D(s, n)
        call run_unit_suites('lib_cart_align3D', cline, s(1:n))
    end subroutine exec_test_lib_cart_align3D

    subroutine exec_test_unit_single( self, cline )
        class(commander_test_unit_single), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_single(s, n)
        call run_unit_suites('unit_single', cline, s(1:n))
    end subroutine exec_test_unit_single

    subroutine exec_test_lib_single( self, cline )
        class(commander_test_lib_single), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_single(s, n)
        call run_unit_suites('lib_single', cline, s(1:n))
    end subroutine exec_test_lib_single

    subroutine exec_test_unit_stream( self, cline )
        class(commander_test_unit_stream), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_stream(s, n)
        call run_unit_suites('unit_stream', cline, s(1:n))
    end subroutine exec_test_unit_stream

    subroutine exec_test_lib_stream( self, cline )
        class(commander_test_lib_stream), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_stream(s, n)
        call run_unit_suites('lib_stream', cline, s(1:n))
    end subroutine exec_test_lib_stream

    subroutine exec_test_unit_parallel( self, cline )
        class(commander_test_unit_parallel), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_parallel(s, n)
        call run_unit_suites('unit_parallel', cline, s(1:n))
    end subroutine exec_test_unit_parallel

    subroutine exec_test_unit_heterogeneity( self, cline )
        class(commander_test_unit_heterogeneity), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_heterogeneity(s, n)
        call run_unit_suites('unit_heterogeneity', cline, s(1:n))
    end subroutine exec_test_unit_heterogeneity

    subroutine exec_test_lib_heterogeneity( self, cline )
        class(commander_test_lib_heterogeneity), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_heterogeneity(s, n)
        call run_unit_suites('lib_heterogeneity', cline, s(1:n))
    end subroutine exec_test_lib_heterogeneity

    !> OpenMP target offload on a device (Cyril's): setup, persistence, async, cuFFT against
    !! FFTW, cuBLAS, the KB device forms; nthr= and device= on the command line; label `platform`,
    !! registered with USE_OPENMP_OFFLOAD; a failed check stops through THROW_HARD
    subroutine exec_test_openmp_offload( self, cline )
        class(commander_test_openmp_offload), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        call run_openmp_offload_tests(cline)
        call simple_end('**** SIMPLE_TEST_OPENMP_OFFLOAD NORMAL STOP ****')
    end subroutine exec_test_openmp_offload

    subroutine exec_test_lib_reconstruction( self, cline )
        class(commander_test_lib_reconstruction), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(unit_suite) :: s(MAX_SUITES)
        integer :: n
        n = 0
        call suites_lib_reconstruction(s, n)
        call run_unit_suites('lib_reconstruction', cline, s(1:n))
    end subroutine exec_test_lib_reconstruction

    subroutine exec_test_forked_process( self, cline )
        class(commander_test_forked_process), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(unit_suite) :: s(2)
        integer :: n
        n = 0
        call add_suite(s, n, 'forked process',   run_all_forked_process_tests)
        call add_suite(s, n, 'stream heartbeat', run_stream_heartbeat_tests)
        call run_unit_suites('forked_process', cline, s(1:n))
    end subroutine exec_test_forked_process

    ! ---- the runner ----------------------------------------------------------

    !> Runs the sub-suites of one area in this process, in its own dated
    !! directory, accumulating failures through simple_test_utils; exits
    !! non-zero if any check failed. `suite=<name>` (from the command line)
    !! runs one sub-suite, and `suite=list` prints the accepted names without
    !! starting a test run; SIMPLE_UNIT_ORDER=reverse walks the table backwards.
    subroutine run_unit_suites( label, cline, suites )
        character(len=*), intent(in)    :: label
        class(cmdline),   intent(inout) :: cline
        type(unit_suite), intent(in)    :: suites(:)
        character(8)          :: datestr
        character(len=STDLEN) :: folder
        type(string)          :: original_cwd, report_file, only_suite
        character(len=32)     :: order_env
        logical               :: test_failed, l_reverse
        integer               :: i, isuite, nrun, iostat
        only_suite = ''
        if( cline%defined('suite') )then
            only_suite = cline%get_carg('suite')
            only_suite = suite_id(only_suite%to_char())
        endif
        if( only_suite == 'list' )then
            write(logfhandle,'(A)') 'Available suites for '//trim(label)//':'
            do i = 1, size(suites)
                write(logfhandle,'(A)') '  '//trim(suite_id(suites(i)%name))
            end do
            return
        endif
        call date_and_time(date=datestr)
        folder = 'SIMPLE_TEST_'//trim(label)//'_'//datestr
        call simple_getcwd(original_cwd)
        report_file = original_cwd//'/'//trim(folder)//'/simple_test_'//trim(label)//'_report.txt'
        call simple_mkdir(folder)
        call simple_chdir(folder)
        call reset_test_report(report_file%to_char())
        ! an optional developer switch: read directly so that an unset variable is silent
        call get_environment_variable('SIMPLE_UNIT_ORDER', value=order_env, status=iostat)
        l_reverse = iostat == 0 .and. trim(order_env) == 'reverse'
        nrun = 0
        do i = 1, size(suites)
            isuite = i
            if( l_reverse ) isuite = size(suites) + 1 - i
            if( only_suite%strlen_trim() > 0 )then
                if( .not. (only_suite == suite_id(suites(isuite)%name)) ) cycle
            endif
            ! every sub-suite starts from one fixed seed, so the full run, a focused run
            ! (suite=<name>) and a reversed run draw the same numbers in it; the seed also restarts
            ! the SIMPLE_SEED count of seed_rnd, which production code calls (parameters%new)
            call set_fixed_seed(20260923)
            call begin_test_suite(trim(suites(isuite)%name))
            call suites(isuite)%run()
            call end_test_suite
            nrun = nrun + 1
        end do
        call report_summary(failed=test_failed)
        call simple_chdir(original_cwd%to_char())
        if( nrun == 0 ) THROW_HARD('no sub-suite '//only_suite%to_char()//' in '//trim(label)//'; use suite=list')
        if( test_failed ) error stop 1
        call simple_end('**** SIMPLE_TEST_'//trim(label)//' NORMAL STOP ****')
    end subroutine run_unit_suites

    !> command-line spelling of a sub-suite name: lowercase, blanks and hyphens as underscores,
    !! commas and slashes dropped ('search, sort, locate' -> search_sort_locate, 'stack I/O' ->
    !! stack_io). Only the name up to its last non-blank counts: the stored names are character(32),
    !! and their padding must not become underscores, or suite= matches nothing;
    !! scripts/check_test_registry.py applies the same rule to the lists in the test UI
    function suite_id( name ) result( id )
        character(len=*), intent(in) :: name
        character(len=len(name))     :: id, lname
        integer :: i, j
        lname = lowercase(adjustl(name))
        id    = ''
        j     = 0
        do i = 1, len_trim(lname)
            select case( lname(i:i) )
                case( ',', '/' )
                    cycle
                case( ' ', '-' )
                    j = j + 1
                    id(j:j) = '_'
                case default
                    j = j + 1
                    id(j:j) = lname(i:i)
            end select
        end do
    end function suite_id

    ! ---- wrappers for test procedures that take arguments -----------------------

    subroutine suite_ui_json
        write(logfhandle,'(a)') 'VALIDATING UI JSON FILE:'
        call validate_ui_json
        write(logfhandle,'(a)') 'PASSED UI JSON FILE TEST'
    end subroutine suite_ui_json

    ! ---- local sub-suites (formerly contained in exec_test_units) -----------------

    subroutine test_euler_shift
        type(ori) :: o
        integer   :: i
        real      :: euls(3), euls_shifted(3)
        logical   :: doshift
        call o%new(is_ptcl=.false.)
        do i=1,100000
            euls(1) = ran3()*800.-400.
            euls(2) = ran3()*500-250.
            euls(3) = ran3()*800.-400.
            call o%set_euler(euls)
            euls_shifted = o%get_euler()
            doshift = .false.
            if( euls_shifted(1) < 0. .or. euls_shifted(1) > 360. ) doshift = .true.
            if( euls_shifted(2) < 0. .or. euls_shifted(2) > 180. ) doshift = .true.
            if( euls_shifted(3) < 0. .or. euls_shifted(3) > 360. ) doshift = .true.
            if( doshift ) THROW_HARD('euler shifting does not work!')
        end do
    end subroutine test_euler_shift

end module simple_commanders_test_class
