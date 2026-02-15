!@descr: Aggregated public API for simple_exec
module simple_exec_api
use simple_core_module_api
use simple_exec_helpers,    only: script_exec, restarted_exec, update_job_descriptions_in_project
use simple_jiffys,          only: simple_print_git_version, simple_print_timer
use simple_ui,              only: make_ui, list_simple_prgs_in_ui, list_simple_test_prgs_in_ui
use iso_fortran_env,        only: output_unit
use simple_cmdline,         only: cmdline, cmdline_err
use simple_exec_project,    only: exec_project_commander
use simple_exec_preproc,    only: exec_preproc_commander
use simple_exec_cluster2D,  only: exec_cluster2D_commander
use simple_exec_cavgproc,   only: exec_cavgproc_commander
use simple_exec_abinitio3D, only: exec_abinitio3D_commander
use simple_exec_refine3D,   only: exec_refine3D_commander
use simple_exec_denoise,    only: exec_denoise_commander
use simple_exec_filter,     only: exec_filter_commander
use simple_exec_image,      only: exec_image_commander
use simple_exec_mask,       only: exec_mask_commander
use simple_exec_ori,        only: exec_ori_commander
use simple_exec_print,      only: exec_print_commander
use simple_exec_sim,        only: exec_sim_commander
use simple_exec_volume,     only: exec_volume_commander
use simple_exec_res,        only: exec_res_commander
use simple_exec_dock,       only: exec_dock_commander
use simple_exec_validate,   only: exec_validate_commander
use simple_exec_other,      only: exec_other_commander
end module simple_exec_api
