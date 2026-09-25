!@descr: Aggregated public API for simple_test_exec
module simple_test_exec_api
use simple_core_module_api
use simple_exec_helpers,        only: script_exec, restarted_exec, update_job_descriptions_in_project
use simple_jiffys,              only: simple_print_git_version, simple_print_timer
use simple_ui,                  only: make_test_ui, list_simple_test_prgs_in_ui, get_test_prg_ptr
use simple_ui_program,          only: ui_program
use iso_fortran_env,            only: output_unit
use simple_cmdline,             only: cmdline, cmdline_err
use simple_test_exec_class,     only: exec_test_class_commander
use simple_test_exec_highlevel, only: exec_test_highlevel_commander
use simple_test_exec_single,    only: exec_test_single_commander
use simple_test_exec_stream,    only: exec_test_stream_commander
end module simple_test_exec_api
