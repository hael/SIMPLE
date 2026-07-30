!@descr: module defining the user interfaces for atom-related programs in the single_exec suite
module single_ui_atom
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('atom', 'Atom Analysis', 50)
type(ui_program), target :: atoms_register
type(ui_program), target :: atoms_rmsd
type(ui_program), target :: atoms_stats
type(ui_program), target :: core_atoms_analysis
type(ui_program), target :: crys_score
type(ui_program), target :: detect_atoms
type(ui_program), target :: simulate_nanoparticle

contains

    subroutine construct_single_atom_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_atoms_register(prgtab)
        call new_atoms_rmsd(prgtab)
        call new_atoms_stats(prgtab)
        call new_core_atoms_analysis(prgtab)
        call new_crys_score(prgtab)
        call new_detect_atoms(prgtab)
        call new_simulate_nanoparticle(prgtab)
    end subroutine construct_single_atom_programs

subroutine new_atoms_register( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_register%new(&
        &'atoms_register',&                                                                           ! name
        &'Align two nanoparticle maps using their detected atom positions',&                          ! summary
        &'is a program that registers two nanoparticles given the maps and the atom position maps.',& ! descr long
        &'single_exec',&                                                                              ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Register Nanoparticle Maps') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call atoms_register%add_input(UI_FILE, 'fname', 'file', 'PDB file list', 'PDB file list', 'e.g. pdb_files.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call atoms_register%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('atoms_register', atoms_register, prgtab, UI_CATEGORY)
    end subroutine new_atoms_register

    subroutine new_atoms_rmsd( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_rmsd%new(&
        &'atoms_rmsd',&                                                               ! name
        &'Measure coordinate changes in atomic nanoparticle time-series data',&           ! summary
        &'is a program that analyzes atomic time-series coordinates',&                ! descr long
        &'single_exec',&                                                              ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Analyze Atomic Trajectories') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call atoms_rmsd%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call atoms_rmsd%add_input(UI_FILE, 'pdbfiles',  'file', 'txt', 'List of PDB format coords files',  'e.g. pdb_files.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call atoms_rmsd%add_input(UI_PARM, 'frac_diam', 'num',  'Fraction of atomic diameter', 'Fraction of atomic diameter used for thresholding{0.5}', '{0.5}', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call atoms_rmsd%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('atoms_rmsd', atoms_rmsd, prgtab, UI_CATEGORY)
    end subroutine new_atoms_rmsd

    subroutine new_atoms_stats( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call atoms_stats%new(&
        &'atoms_stats',&                                                                              ! name
        &'Compute radial statistics to assess symmetry in a nanoparticle map',&                      ! summary
        &'is a program that generates statistics at different radii and across the whole nano map.',& ! descr long
        &'single_exec',&                                                                              ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Analyze Radial Symmetry') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call atoms_stats%add_input(UI_IMG, 'vol1', 'file', 'Raw volume', 'Raw volume of grey valued pixel intensities', &
        & 'input volume e.g. vol.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call atoms_stats%add_input(UI_IMG, 'vol2', 'file', 'Connected components volume', 'Connected components volume produced by detect atoms', &
        & 'input volume e.g. *CC.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call atoms_stats%add_input(UI_IMG, 'vol3', 'file', 'Volume', 'Nanoparticle volume to use for lattice fitting', &
        & 'input volume 4 lattice fit e.g. vol3.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call atoms_stats%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call atoms_stats%add_input(UI_FILE, 'pdbfile', 'file', 'PDB', 'Input coords file in PDB format', 'e.g. atoms.pdb', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call atoms_stats%add_input(UI_FILE, 'pdbfile2', 'file', 'PDB', 'subset coords for stats calc', 'e.g. subset.pdb', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call atoms_stats%add_input(UI_FILE, 'rmsd_file','file', 'bin', 'per-atom e/o rmsd:s', 'e.g. atom_rmsd.bin', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call atoms_stats%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! <empty>
        ! computer controls
        call atoms_stats%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('atoms_stats', atoms_stats, prgtab, UI_CATEGORY)
    end subroutine new_atoms_stats

    subroutine new_core_atoms_analysis( prgtab )
        class(ui_hash), intent(inout) :: prgtab        
        ! PROGRAM SPECIFICATION
        call core_atoms_analysis%new(&
        &'core_atoms_analysis',&                                                      ! name
        &'Analysis of results obtianed with trajectory_reconstruct3D and detect_atoms',& ! summary
        &'is a program that analysis atomic time-series coordinates',&                ! descr long
        &'single_exec',&                                                              ! executable
        &.false., visibility=UI_VIS_DEVELOPER)                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call core_atoms_analysis%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call core_atoms_analysis%add_input(UI_FILE, 'pdbfiles',  'file', 'txt', 'List of PDB format coords files',  'e.g. pdb_files.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call core_atoms_analysis%add_input(UI_PARM, 'frac_diam', 'num',  'Fraction of atomic diameter', 'Fraction of atomic diameter used for thresholding{0.5}', '{0.5}', .false., 0.5, &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call core_atoms_analysis%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('core_atoms_analysis', core_atoms_analysis, prgtab, UI_CATEGORY)
    end subroutine new_core_atoms_analysis

    subroutine new_crys_score( prgtab )
        class(ui_hash), intent(inout) :: prgtab        
        ! PROGRAM SPECIFICATION
        call crys_score%new(&
        &'crys_score',&                                                 ! name
        &'Calculate a crystal-quality score for a nanoparticle density map',& ! summary
        &'is a program that computes crystal score.',&                  ! descr long
        &'single_exec',&                                                ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Score Crystal Structure') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call crys_score%add_input(UI_FILE, 'fname', 'file', 'PDB file list', 'PDB file list', 'e.g. np_pdbs.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call crys_score%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call crys_score%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call crys_score%add_input(UI_FILT, element, &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! <empty>
        ! computer controls
        call crys_score%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('crys_score', crys_score, prgtab, UI_CATEGORY)
    end subroutine new_crys_score

    subroutine new_detect_atoms( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call detect_atoms%new(&
        &'detect_atoms', &                                      ! name
        &'Identify atoms and connected regions in an atomic-resolution nanoparticle map',& ! summary
        &'is a program for identifying atoms in atomic-resolution nanoparticle maps and generating bin and connected-components map',& ! descr long
        &'single_exec',&                                        ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Detect Atoms') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call detect_atoms%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Nanoparticle volume to analyse', &
        & 'input volume e.g. vol.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call detect_atoms%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call detect_atoms%add_input(UI_FILT, element, &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call detect_atoms%add_input(UI_MASK, mskdiam, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call detect_atoms%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('detect_atoms', detect_atoms, prgtab, UI_CATEGORY)
    end subroutine new_detect_atoms

    subroutine new_simulate_nanoparticle( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_nanoparticle%new(&
        &'simulate_nanoparticle',&                                              ! name
        &'Simulate nanoparticle for lattice density',&                          ! summary
        &'is a program for simulation of nanoparticle for lattice density',&    ! help
        &'single_exec',&                                                        ! executable
        &.false., visibility=UI_VIS_ADVANCED)                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_nanoparticle%add_input(UI_FILE, 'pdbfile', 'file', 'PDB', 'Input coordinates file in PDB format', 'Input coordinates file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call simulate_nanoparticle%add_input(UI_FILE, 'pdbout', 'file', 'Output PDB', &
        &'Output coordinates file for simulated atoms', 'output coordinates e.g. simatms.pdb', .false., 'simatms.pdb', &
        &visibility=UI_VIS_ADVANCED)
        call simulate_nanoparticle%add_input(UI_IMG, outvol, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call simulate_nanoparticle%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call simulate_nanoparticle%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        call simulate_nanoparticle%add_input(UI_PARM, element, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_nanoparticle%add_input(UI_PARM, moldiam, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_nanoparticle%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('simulate_nanoparticle', simulate_nanoparticle, prgtab, UI_CATEGORY)
    end subroutine new_simulate_nanoparticle


end module single_ui_atom
