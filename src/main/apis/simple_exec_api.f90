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

! validate commanders, for validation of parts of the pipelined stream process
use simple_commanders_validate, only: commander_mini_stream, commander_check_refpick

! sim commanders, for simulation of nosie, particls & movies
use simple_commanders_sim, only: commander_simulate_noise, commander_simulate_particles, commander_simulate_movie

! atoms commanders, routines involving PDB files
use simple_commanders_atoms, only: commander_map2model_fsc, commander_model_validation

! distr commanders, support routines for distributed execution
use simple_commanders_distr, only: commander_split

! test commanders, for testing purposes
use simple_commanders_test, only: commander_test_sim_workflow

end module simple_exec_api
