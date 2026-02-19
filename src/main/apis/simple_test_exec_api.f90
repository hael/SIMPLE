!@descr: Aggregated public API for simple_test_exec
module simple_test_exec_api
use simple_core_module_api
use simple_exec_helpers,        only: script_exec, restarted_exec, update_job_descriptions_in_project
use simple_jiffys,              only: simple_print_git_version, simple_print_timer
use simple_ui,                  only: make_test_ui, list_simple_test_prgs_in_ui
use iso_fortran_env,            only: output_unit
use simple_cmdline,             only: cmdline, cmdline_err
use simple_test_exec_fft,       only: exec_test_fft_commander
use simple_test_exec_geometry,  only: exec_test_geometry_commander
use simple_test_exec_highlevel, only: exec_test_highlevel_commander
use simple_test_exec_io,        only: exec_test_io_commander
use simple_test_exec_masks,     only: exec_test_masks_commander
use simple_test_exec_network,   only: exec_test_network_commander
use simple_test_exec_numerics,  only: exec_test_numerics_commander
use simple_test_exec_optimize,  only: exec_test_optimize_commander
use simple_test_exec_parallel,  only: exec_test_parallel_commander
use simple_test_exec_stats,     only: exec_test_stats_commander
use simple_test_exec_utils,     only: exec_test_utils_commander
end module simple_test_exec_api
