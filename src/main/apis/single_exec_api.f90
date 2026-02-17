!@descr: Aggregated public API for single_exec
module single_exec_api
use simple_core_module_api
use simple_exec_helpers,    only: script_exec, restarted_exec, update_job_descriptions_in_project
use simple_jiffys,          only: simple_print_git_version, simple_print_timer
use simple_ui,              only: make_ui, list_single_prgs_in_ui
use iso_fortran_env,        only: output_unit
use simple_cmdline,         only: cmdline, cmdline_err
use simple_ui,              only: make_ui, list_single_prgs_in_ui
use single_exec_tseries,    only: exec_tseries_commander
use single_exec_trajectory, only: exec_trajectory_commander
use single_exec_nano2D,     only: exec_nano2D_commander
use single_exec_nano3D,     only: exec_nano3D_commander
use single_exec_map,        only: exec_map_commander
use single_exec_atom,       only: exec_atom_commander
use single_exec_validate,   only: exec_validate_commander
end module single_exec_api

