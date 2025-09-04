! executes the shared-memory parallelised programs in SIMPLE_STREAM
program simple_stream
include 'simple_lib.f08'
use simple_user_interface, only: make_user_interface, list_stream_prgs_in_ui
use simple_cmdline,        only: cmdline, cmdline_err
use simple_exec_helpers
use simple_commander_stream
use simple_commander_stream2D

implicit none
#include "simple_local_flags.inc"

! PROGRAMS
type(commander_stream_preprocess)           :: xpreprocess
type(commander_stream_pick_extract)         :: xpick_extract
type(commander_stream_gen_picking_refs)     :: xgen_picking_refs
type(commander_stream_assign_optics)        :: xassign_optics
type(commander_stream_sieve_cavgs)          :: xsieve_cavgs
type(commander_stream_cluster2D)            :: xcluster2D_stream
type(commander_stream_abinitio2D)           :: xabinitio2D_stream

! OTHER DECLARATIONS
character(len=STDLEN)                       :: xarg, prg, entire_line
type(cmdline)                               :: cline
integer                                     :: cmdstat, cmdlen, pos
integer(timer_int_kind)                     :: t0
real(timer_int_kind)                        :: rt_exec

! start timer
t0 = tic()
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_stream_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
! generate script for queue submission?
call script_exec(cline, trim(prg), 'simple_stream')

select case(trim(prg))
    case( 'preproc' )
        call xpreprocess%execute(cline)
    case( 'pick_extract' )
        call xpick_extract%execute(cline)
    case( 'gen_picking_refs' )
        call xgen_picking_refs%execute(cline)
    case( 'assign_optics' )
        call xassign_optics%execute(cline)
    case( 'sieve_cavgs' )
        call xsieve_cavgs%execute(cline)
    case( 'cluster2D_stream' )
        call xcluster2D_stream%execute(cline)
    case( 'abinitio2D_stream' )
        call xabinitio2D_stream%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
! close log file
if( logfhandle .ne. OUTPUT_UNIT )then
    if( is_open(logfhandle) ) call fclose(logfhandle)
endif
call simple_print_git_version('dc6dd592')
! end timer and print
rt_exec = toc(t0)
call simple_print_timer(rt_exec)
end program simple_stream
