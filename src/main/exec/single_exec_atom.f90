!@descr:
module single_exec_atom
use simple_cmdline,          only: cmdline
use simple_commanders_atoms, only: commander_conv_atom_denoise, commander_atoms_stats, commander_atoms_register,&
commander_crys_score, commander_atoms_rmsd, commander_core_atoms_analysis, commander_detect_atoms
use simple_commanders_sim,   only: commander_simulate_atoms
implicit none

public :: exec_atom_commander
private

type(commander_atoms_register)      :: xatoms_register
type(commander_atoms_rmsd)          :: xatoms_rmsd
type(commander_atoms_stats)         :: xatoms_stats
type(commander_core_atoms_analysis) :: xcore_atoms_analysis
type(commander_crys_score)          :: xcrys_score
type(commander_detect_atoms)        :: xdetect_atoms
type(commander_simulate_atoms)      :: xsimulate_atoms

contains

    subroutine exec_atom_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(inout) :: l_did_execute
        logical,             intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'atoms_register' )
                call cline%set('mkdir', 'no')
                call xatoms_register%execute(cline)
            case( 'atoms_rmsd' )
                call xatoms_rmsd%execute(cline)
            case( 'atoms_stats' )
                call cline%set('mkdir', 'yes')
                call xatoms_stats%execute(cline)
            case( 'core_atoms_analysis' )
                call xcore_atoms_analysis%execute(cline)
            case( 'crys_score' )
                call cline%set('mkdir', 'no')
                call xcrys_score%execute(cline)
            case( 'detect_atoms' )
                call cline%set('mkdir', 'no')
                call xdetect_atoms%execute(cline)
            case( 'simulate_atoms' )
                call cline%set('mkdir', 'no')
                call xsimulate_atoms%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_atom_commander

end module single_exec_atom