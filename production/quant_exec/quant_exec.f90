program quant_exec
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline, cmdline_err
use simple_spproj_hlev,    only: update_job_descriptions_in_project
use simple_commander_quant
use simple_user_interface
implicit none
#include "simple_local_flags.inc"
! QUANT PROGRAMS
type(detect_atoms_commander)           :: xdetect_atoms
type(atoms_rmsd_commander)             :: xatoms_rmsd
type(radial_dependent_stats_commander) :: xradial_dependent_stats
type(atom_cluster_analysis_commander)  :: xatom_cluster_analysis
type(nano_softmask_commander)          :: xnano_softmask
type(geometry_analysis_commander)      :: xgeometry_analysis
type(radial_sym_test_commander)        :: xrad_sym_tst
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
    ! QUANT PROGRAMS
    case( 'detect_atoms' )
        call xdetect_atoms%execute(cline)
    case( 'radial_dependent_stats' )
        call xradial_dependent_stats%execute(cline)
    case( 'atom_cluster_analysis' )
        call xatom_cluster_analysis%execute(cline)
    case( 'nano_softmask' )
        call xnano_softmask%execute(cline)
    case('geometry_analysis')
          call xgeometry_analysis%execute(cline)
    case('atoms_rmsd')
        call xatoms_rmsd%execute(cline)
    case('radial_sym_test')
        call xrad_sym_tst%execute(cline)
    case DEFAULT
        THROW_HARD('prg='//trim(prg)//' is unsupported')
end select
call update_job_descriptions_in_project( cline )
end program quant_exec
