!@descr: runs select2D, the GTK 3 window for selecting a project's class averages (simple_viewer projfile=<p>.simple)
program simple_viewer
use iso_fortran_env,     only: output_unit
use simple_defs,         only: logfhandle, SIMPLE_GIT_VERSION
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_timer,        only: timer_int_kind, tic, toc
use simple_syslib,       only: is_open
use simple_fileio,       only: fclose
use simple_jiffys,       only: simple_print_git_version, simple_print_timer
use simple_cmdline,      only: cmdline
use simple_ui,           only: make_ui
use simple_viewer_views, only: exec_viewer_commander
implicit none
#include "simple_local_flags.inc"
type(cmdline)           :: cline
type(string)            :: prg
integer(timer_int_kind) :: t0
real(timer_int_kind)    :: rt_exec
logical                 :: l_silent, l_did_execute
t0 = tic()
call make_ui
! no prg= needed: cmdline%parse sets the one program the registry holds for simple_viewer
call cline%parse
prg = cline%get_carg('prg')
l_silent      = .false.
l_did_execute = .false.
call exec_viewer_commander(prg%to_char(), cline, l_silent, l_did_execute)
if( .not. l_did_execute ) THROW_HARD('Program "'//prg%to_char()//'" not recognized by simple_viewer')
! no job description is appended to the input project: a save records the job in the project it writes
if( logfhandle .ne. output_unit )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
if( .not. l_silent )then
    call simple_print_git_version(SIMPLE_GIT_VERSION)
    rt_exec = toc(t0)
    call simple_print_timer(rt_exec)
endif
end program simple_viewer
