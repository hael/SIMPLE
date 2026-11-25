!@descr: runs the simple_viewer programs: select2D, the GTK 3 window for selecting a project's class averages (simple_viewer prg=select2D projfile=<p>.simple)
program simple_viewer
use iso_fortran_env,     only: output_unit
use simple_defs,         only: logfhandle, SIMPLE_GIT_VERSION, STDLEN, XLONGSTRLEN
use simple_error,        only: simple_exception
use simple_string,       only: string
use simple_string_utils, only: str_has_substr
use simple_timer,        only: timer_int_kind, tic, toc
use simple_syslib,       only: is_open
use simple_fileio,       only: fclose
use simple_jiffys,       only: simple_print_git_version, simple_print_timer
use simple_cmdline,      only: cmdline, cmdline_err
use simple_ui,           only: make_ui, get_prg_ptr, list_viewer_prgs_in_ui
use simple_ui_program,   only: ui_program
use simple_viewer_views, only: exec_viewer_commander
implicit none
#include "simple_local_flags.inc"
character(len=STDLEN)      :: xarg
character(len=XLONGSTRLEN) :: entire_line
type(cmdline)              :: cline
type(string)               :: prg
type(ui_program), pointer  :: ptr2prg
integer                    :: cmdstat, cmdlen, pos
integer(timer_int_kind)    :: t0
real(timer_int_kind)       :: rt_exec
logical                    :: l_silent, l_did_execute
t0 = tic()
! the program name, as simple_exec reads it
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err(cmdstat, cmdlen, xarg, pos)
prg = trim(xarg(pos+1:))
call make_ui
if( str_has_substr(entire_line, 'prg=list') )then
    call list_viewer_prgs_in_ui
    stop
endif
! validate program name before parsing (gives a friendly error for typos)
ptr2prg => null()
call get_prg_ptr(prg, ptr2prg)
if( .not. associated(ptr2prg) )then
    THROW_HARD('Program "'//prg%to_char()//'" not recognized. Use prg=list to see available programs.')
endif
call cline%parse
prg = cline%get_carg('prg')
l_silent      = .false.
l_did_execute = .false.
call exec_viewer_commander(prg%to_char(), cline, l_silent, l_did_execute)
if( .not. l_did_execute ) THROW_HARD('Program "'//prg%to_char()//'" not recognized by simple_viewer. Use prg=list to see available programs.')
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
