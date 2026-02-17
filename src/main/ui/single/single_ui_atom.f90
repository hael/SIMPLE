!@descr: module defining the user interfaces for atom-related programs in the single_exec suite
module single_ui_atom
use simple_ui_modules
implicit none

type(ui_program), target :: atoms_register
type(ui_program), target :: atoms_rmsd
type(ui_program), target :: atoms_stats
type(ui_program), target :: core_atoms_analysis
type(ui_program), target :: crys_score
type(ui_program), target :: detect_atoms
type(ui_program), target :: simulate_atoms

contains

    subroutine construct_single_atom_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_atoms_register(prgtab)
        call new_atoms_rmsd(prgtab)
        call new_atoms_stats(prgtab)
        call new_core_atoms_analysis(prgtab)
        call new_crys_score(prgtab)
        call new_detect_atoms(prgtab)
        call new_simulate_atoms(prgtab)
    end subroutine construct_single_atom_programs

    subroutine print_single_atom_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('ATOM ANALYSIS PROGRAMS:', C_UNDERLINED)
        write(logfhandle,'(A)') atoms_register%name%to_char()
        write(logfhandle,'(A)') atoms_rmsd%name%to_char()
        write(logfhandle,'(A)') atoms_stats%name%to_char()
        write(logfhandle,'(A)') core_atoms_analysis%name%to_char()
        write(logfhandle,'(A)') crys_score%name%to_char()
        write(logfhandle,'(A)') detect_atoms%name%to_char()
        write(logfhandle,'(A)') simulate_atoms%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_atom_programs

    subroutine new_atoms_register( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_register%new(&
        &'atoms_register',&                                                                           ! name
        &'Registration of two nanoparticles',&                                                        ! descr_short
        &'is a program that registers two nanoparticles given the maps and the atom position maps.',& ! descr long
        &'single_exec',&                                                                              ! executable
        &.false., gui_advanced=.false.)                                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call atoms_register%add_input(UI_IMG, 'fname', 'file', 'PDB file list', 'PDB file list', 'e.g. pdb_files.txt', .true., '')
        ! parameter input/output
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call atoms_register%add_input(UI_COMP, nthr)
        call add_ui_program('atoms_register', atoms_register, prgtab) 
    end subroutine new_atoms_register

    subroutine new_atoms_rmsd( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_rmsd%new(&
        &'atoms_rmsd',&                                                               ! name
        &'Analysis of results obtained with trajectory_reconstruct3D and detect_atoms',& ! descr_short
        &'is a program that analyzes atomic time-series coordinates',&                ! descr long
        &'single_exec',&                                                              ! executable
        &.false., gui_advanced=.false.)                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call atoms_rmsd%add_input(UI_PARM, smpd)
        call atoms_rmsd%add_input(UI_PARM, 'pdbfiles',  'file', 'txt', 'List of PDB format coords files',  'List of input coords files in PDB format', .true., '')
        call atoms_rmsd%add_input(UI_PARM, 'frac_diam', 'num',  'Fraction of atomic diameter', 'Fraction of atomic diameter used for thresholding{0.5}', '{0.5}', .false., 0.5)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call atoms_rmsd%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        call add_ui_program('atoms_rmsd', atoms_rmsd, prgtab) 
    end subroutine new_atoms_rmsd

    subroutine new_atoms_stats( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_stats%new(&
        &'atoms_stats',&                                                                              ! name
        &'Statistical test for radial dependent symmetry',&                                           ! descr_short
        &'is a program that generates statistics at different radii and across the whole nano map.',& ! descr long
        &'single_exec',&                                                                              ! executable
        &.false., gui_advanced=.false.)                                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call atoms_stats%add_input(UI_IMG, 'vol1', 'file', 'Raw volume', 'Raw volume of grey valued pixel intensities', &
        & 'input volume e.g. vol.mrc', .true., '')
        call atoms_stats%add_input(UI_IMG, 'vol2', 'file', 'Connected components volume', 'Connected components volume produced by detect atoms', &
        & 'input volume e.g. *CC.mrc', .true., '')
        call atoms_stats%add_input(UI_IMG, 'vol3', 'file', 'Volume', 'Nanoparticle volume to use for lattice fitting', &
        & 'input volume 4 lattice fit e.g. vol3.mrc', .false., '')
        ! parameter input/output
        call atoms_stats%add_input(UI_PARM, smpd)
        call atoms_stats%add_input(UI_PARM, 'pdbfile', 'file', 'PDB', 'Input coords file in PDB format', 'Input coords file in PDB format', .true., '')
        call atoms_stats%add_input(UI_PARM, 'pdbfile2', 'file', 'PDB', 'subset coords for stats calc', 'subset coords file in PDB format for stats calc', .false., '')
        call atoms_stats%add_input(UI_PARM, 'rmsd_file','file', 'bin', 'per-atom e/o rmsd:s', 'per-atom e/o rmsd:s from CS model building', .false., '')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call atoms_stats%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        ! mask controls
        ! <empty>
        ! computer controls
        call atoms_stats%add_input(UI_COMP, nthr)
        call add_ui_program('atoms_stats', atoms_stats, prgtab) 
    end subroutine new_atoms_stats

    subroutine new_core_atoms_analysis( prgtab )
        class(ui_hash), intent(inout) :: prgtab        
        ! PROGRAM SPECIFICATION
        call core_atoms_analysis%new(&
        &'core_atoms_analysis',&                                                      ! name
        &'Analysis of results obtianed with trajectory_reconstruct3D and detect_atoms',& ! descr_short
        &'is a program that analysis atomic time-series coordinates',&                ! descr long
        &'single_exec',&                                                              ! executable
        &.false., gui_advanced=.false.)                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call core_atoms_analysis%add_input(UI_PARM, smpd)
        call core_atoms_analysis%add_input(UI_PARM, 'pdbfiles',  'file', 'txt', 'List of PDB format coords files',  'List of input coords files in PDB format', .true., '')
        call core_atoms_analysis%add_input(UI_PARM, 'frac_diam', 'num',  'Fraction of atomic diameter', 'Fraction of atomic diameter used for thresholding{0.5}', '{0.5}', .false., 0.5)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call core_atoms_analysis%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '')
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        call add_ui_program('core_atoms_analysis', core_atoms_analysis, prgtab) 
    end subroutine new_core_atoms_analysis

    subroutine new_crys_score( prgtab )
        class(ui_hash), intent(inout) :: prgtab        
        ! PROGRAM SPECIFICATION
        call crys_score%new(&
        &'crys_score',&                                                 ! name
        &'Computing crystal score',&                                    ! descr_short
        &'is a program that computes crystal score.',&                  ! descr long
        &'single_exec',&                                                ! executable
        &.false., gui_advanced=.false.)                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call crys_score%add_input(UI_IMG, 'fname',   'file', 'Core PDB folder list', 'Core PDB folder list', 'e.g. core_pdbs.txt', .true., '')
        call crys_score%add_input(UI_IMG, 'pdbfile', 'file', 'PDB input coordinates file to estimate moldiam', 'Input coordinates file in PDB format', 'PDB file e.g. startvol_ATMS.pdb', .true., 'startvol_ATMS.pdb')
        ! parameter input/output
        call crys_score%add_input(UI_PARM, smpd)
        call crys_score%add_input(UI_PARM, box)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call crys_score%add_input(UI_FILT, element)
        ! mask controls
        ! <empty>
        ! computer controls
        call crys_score%add_input(UI_COMP, nthr)
        call add_ui_program('crys_score', crys_score, prgtab) 
    end subroutine new_crys_score

    subroutine new_detect_atoms( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call detect_atoms%new(&
        &'detect_atoms', &                                      ! name
        &'Detect atoms in atomic-resolution nanoparticle map',& ! descr_short
        &'is a program for identifying atoms in atomic-resolution nanoparticle maps and generating bin and connected-components map',& ! descr long
        &'single_exec',&                                        ! executable
        &.false., gui_advanced=.false.)                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call detect_atoms%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Nanoparticle volume to analyse', &
        & 'input volume e.g. vol.mrc', .true., '')
        call detect_atoms%add_input(UI_IMG, 'vol2', 'file', 'Volume', 'Nanoparticle volume to use for lattice fitting', &
        & 'input volume 4 lattice fit e.g. vol2.mrc', .false., '')
        ! parameter input/output
        call detect_atoms%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call detect_atoms%add_input(UI_FILT, element)
        ! mask controls
        ! <empty>
        ! computer controls
        call detect_atoms%add_input(UI_COMP, nthr)
        call add_ui_program('detect_atoms', detect_atoms, prgtab) 
    end subroutine new_detect_atoms

    subroutine new_simulate_atoms( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_atoms%new(&
        &'simulate_atoms',&                                              ! name
        &'Simulate atoms or FCC lattice density',&                       ! descr_short
        &'is a program for simulation of atoms or FCC lattice density',& ! descr_long
        &'single_exec',&                                                 ! executable
        &.false., gui_advanced=.false.)                                  ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_atoms%add_input(UI_IMG, 'pdbfile', 'file', 'PDB', 'Input coordinates file in PDB format', 'Input coordinates file', .false., '')
        call simulate_atoms%add_input(UI_IMG, outvol)
        ! parameter input/output
        call simulate_atoms%add_input(UI_PARM, smpd)
        call simulate_atoms%add_input(UI_PARM, box)
        call simulate_atoms%add_input(UI_PARM, element, required_override=.false.)
        call simulate_atoms%add_input(UI_PARM, moldiam)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_atoms%add_input(UI_COMP, nthr)
        call add_ui_program('simulate_atoms', simulate_atoms, prgtab) 
    end subroutine new_simulate_atoms


end module single_ui_atom
