! QUANT = Quantitative Unsupervised analysis of Atomic-resolution NanoparTicle 3D charge density maps
program quant_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_spproj_hlev,    only: update_job_descriptions_in_project
use simple_commander_quant
use simple_user_interface
implicit none
#include "simple_local_flags.inc"
type(detect_atoms_commander)           :: xdetect_atoms
type(nano_softmask_commander)          :: xnano_softmask
type(atoms_stats_commander)            :: xatoms_stats
type(tseries_atoms_analysis_commander) :: xtseries_atoms_analysis
type(dock_coords_commander)            :: xdock_coords
! OTHER DECLARATIONS
character(len=STDLEN) :: xarg, prg, entire_line
type(cmdline)         :: cline
integer               :: cmdstat, cmdlen, pos
! parse command-line
call get_command_argument(1, xarg, cmdlen, cmdstat)
call get_command(entire_line)
pos = index(xarg, '=') ! position of '='
call cmdline_err( cmdstat, cmdlen, xarg, pos )
prg = xarg(pos+1:)     ! this is the program name
! make UI
call make_user_interface
if( str_has_substr(entire_line, 'prg=list') )then
    call list_quant_prgs_in_ui
    stop
endif
! parse command line into cline object
call cline%parse
select case(prg)
    case( 'detect_atoms' )
        call cline%set('mkdir', 'no')
        call xdetect_atoms%execute(cline)
    case( 'atoms_stats' )
        call cline%set('mkdir', 'yes')
        call xatoms_stats%execute(cline)
    case( 'tseries_atoms_analysis' )
        call xtseries_atoms_analysis%execute(cline)
    case( 'nano_softmask' )
        call xnano_softmask%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
end program quant_exec
