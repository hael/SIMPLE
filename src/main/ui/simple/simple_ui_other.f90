!@descr: module defining the user interfaces for miscellaneous programs in the simple_exec suite
module simple_ui_other
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('other', 'Other Utilities', 190)
type(ui_program), target :: cif2pdb
type(ui_program), target :: sigma2_convert

contains

    subroutine construct_other_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_cif2pdb(prgtab)
        call new_sigma2_convert(prgtab)
    end subroutine construct_other_programs

    subroutine new_cif2pdb( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call cif2pdb%new(&
        &'cif2pdb',&                                       ! name
        &'Convert PDBx/mmCIF atomic coordinates into PDB format',& ! summary
        &'is a program for converting PDBx/mmCIF to PDB',& ! help
        &'simple_exec',&                                   ! executable
        &.false., &
        &visibility=UI_VIS_STANDARD, display_name='Convert mmCIF to PDB') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cif2pdb%add_input(UI_FILE, 'ciffile', 'file', 'PDBx/mmCIF input coordinates file', 'Input coordinates file in PDBx/mmCIF format', 'PDBx/mmCIF file e.g. molecule.cif', .true., 'molecule.cif', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        ! computer controls
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>               
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('cif2pdb', cif2pdb, prgtab, UI_CATEGORY)
    end subroutine new_cif2pdb

    subroutine new_sigma2_convert( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call sigma2_convert%new(&
        &'sigma2_convert',&
        &'Import or export canonical sigma2 state explicitly',&
        &'converts complete legacy particle files or grouped STAR seeds into canonical state, or exports canonical groups to STAR',&
        &'simple_exec',&
        &.true., visibility=UI_VIS_DEVELOPER)
        call sigma2_convert%add_input(UI_PARM, 'sigma_action', 'multi', 'Conversion action', &
        &'Conversion action(star_import|parts_import|star_export)', '', .true., 'star_export', &
        &choices=ui_choices([character(len=12) :: 'star_import', 'parts_import', 'star_export']), &
        &visibility=UI_VIS_STANDARD)
        call sigma2_convert%add_input(UI_PARM, 'oritype', 'multi', 'Particle orientation field', &
        &'Particle orientation field(ptcl2D|ptcl3D){ptcl3D}', '', .false., 'ptcl3D', &
        &choices=ui_choices([character(len=6) :: 'ptcl2D', 'ptcl3D']), visibility=UI_VIS_DEVELOPER)
        call sigma2_convert%add_input(UI_FILE, 'infile', 'file', 'Input path', &
        &'Grouped STAR input, or the complete legacy part-file prefix ending before the padded part number', &
        &'e.g. sigma2_groups.star', .false., '', visibility=UI_VIS_DEVELOPER)
        call sigma2_convert%add_input(UI_FILE, 'outfile', 'file', 'Output path', &
        &'Canonical state output for import, or STAR output for export', 'e.g. sigma2_state.bin', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call sigma2_convert%add_input(UI_SRCH, sigma_est, visibility=UI_VIS_DEVELOPER)
        call sigma2_convert%add_input(UI_COMP, nparts, required_override=.false., visibility=UI_VIS_DEVELOPER)
        call add_ui_program('sigma2_convert', sigma2_convert, prgtab, UI_CATEGORY)
    end subroutine new_sigma2_convert

end module simple_ui_other
