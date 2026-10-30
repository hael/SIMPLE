!@descr: module defining the user interfaces for map post-processing programs in the simple_exec suite
module simple_ui_postprocess
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('postprocess', 'Post-processing', 69)
type(ui_program), target :: postprocess
type(ui_program), target :: postprocess_nu
type(ui_program), target :: postprocess_pdb
type(ui_program), target :: enhance

contains

    subroutine construct_postprocess_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_postprocess(prgtab)
        call new_postprocess_nu(prgtab)
        call new_postprocess_pdb(prgtab)
        call new_enhance(prgtab)
    end subroutine construct_postprocess_programs

    subroutine new_postprocess( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',&                                                      ! name
        &'Filter and sharpen a reconstructed density map for interpretation',& ! summary
        &'is a program for map post-processing. Use program volops to estimate the B-factor with the Guinier plot',& ! help
        &'simple_exec',&                                                      ! executable
        &.true., &
        &visibility=UI_VIS_STANDARD, display_name='Post-process Density Map') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess%add_input(UI_IMG, 'vol1', 'file', 'Volume override', 'Volume override for selected state', &
            &'input volume e.g. vol_state01.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call postprocess%add_input(UI_PARM, 'state', 'num', 'State to postprocess', 'State to postprocess{1}', 'Input state{1}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call postprocess%add_input(UI_PARM, 'imgkind', 'str', 'Volume kind', 'Project output volume kind; unfil = the unfiltered (base) pair average beside the volume, '//&
            &'solvent = the solvent-prior pair average (pcg_solvent=yes){vol}', &
            &'project output kind: vol, vol_cavg, unfil or solvent', .false., 'vol', &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false., 20., &
        &visibility=UI_VIS_ADVANCED)
        call postprocess%add_input(UI_FILT, 'fsc', 'file', 'FSC file', 'Binary FSC file setting the FSC=0.143 cutoff when no '//&
            &'_even_unfil/_odd_unfil pair sits beside the volume', &
            &'e.g. fsc_state01.bin file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call postprocess%add_input(UI_FILT, bfac, &
        &visibility=UI_VIS_ADVANCED)
        call postprocess%add_input(UI_FILT, mirr, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call postprocess%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call postprocess%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('postprocess', postprocess, prgtab, UI_CATEGORY)
    end subroutine new_postprocess

    subroutine new_postprocess_nu( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call postprocess_nu%new(&
        &'postprocess_nu',&                                    ! name
        &'Nonuniform evidence-bounded postprocessing of even/odd half-maps',& ! summary
        &'is a program for NU-evidence local sharpening of a 3D reconstruction: local amplitude restoration in which &
        &both the confidence field and the target spectrum derive from cross-half NU evidence (model-free LocScale &
        &analogue). Operates on the project like postprocess: the state volume and its unregularized even/odd pair &
        &are taken from the out segment; when the regularized pair is present the refinement filter competition is run &
        &first and its references and local-resolution map are written. Outputs carry the _pproc_nu suffix and are &
        &display/interpretation maps only, never inputs to FSC correction or resolution estimation. Isolated from the &
        &standard postprocess program (global B-factor), which is unchanged',& ! help
        &'simple_exec',&                                       ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED, display_name='Nonuniform Postprocessing')                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess_nu%add_input(UI_IMG, outvol, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call postprocess_nu%add_input(UI_PARM, 'state', 'num', 'State to postprocess', 'State to postprocess{1}', 'Input state{1}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call postprocess_nu%add_input(UI_FILT, mirr, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        ! mask controls
        call postprocess_nu%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call postprocess_nu%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('postprocess_nu', postprocess_nu, prgtab, UI_CATEGORY)
    end subroutine new_postprocess_nu

    subroutine new_enhance( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call enhance%new(&
        &'enhance',&
        &'Enhance unprocessed half-maps with local agreement filtering and bounded sharpening',&
        &'is a deterministic, model-free half-map postprocessor. It estimates outer-shell noise and a conservative '//&
        &'molecular envelope, derives local FSC support in overlapping tapered windows, applies FSC-derived Wiener '//&
        &'suppression and bounded local Guinier sharpening, and overlap-adds the windows. It writes outvol plus '//&
        &'matching _confidence, _locres and _envelope maps. All outputs are interpretation products; only the input '//&
        &'half-maps retain resolution authority',&
        &'simple_exec',&
        &.false., visibility=UI_VIS_STANDARD, display_name='Enhance Half-maps')
        call enhance%add_input(UI_IMG, 'vol_even', 'file', 'Even half-map', &
            &'Unprocessed even half-map with independent noise', 'even_half.mrc file', .true., '', &
            &visibility=UI_VIS_STANDARD)
        call enhance%add_input(UI_IMG, 'vol_odd', 'file', 'Odd half-map', &
            &'Unprocessed odd half-map with independent noise', 'odd_half.mrc file', .true., '', &
            &visibility=UI_VIS_STANDARD)
        call enhance%add_input(UI_IMG, outvol, required_override=.true., &
            &visibility=UI_VIS_STANDARD)
        call enhance%add_input(UI_PARM, smpd, required_override=.false., &
            &visibility=UI_VIS_ADVANCED)
        call enhance%add_input(UI_COMP, nthr, &
            &visibility=UI_VIS_STANDARD)
        call add_ui_program('enhance', enhance, prgtab, UI_CATEGORY)
    end subroutine new_enhance

    subroutine new_postprocess_pdb( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call postprocess_pdb%new(&
        &'postprocess_pdb',&
        &'Post-process a density map inside a docked atomic-model envelope',&
        &'docks PDB-derived density to the input map over 100-20 Angstroms, builds a soft envelope with a '//&
        &'six-pixel cosine edge, and applies half-map FSC weighting, low-pass filtering and global B-factor '//&
        &'sharpening. Requires matching <map>_even and <map>_odd half-maps beside the input volume. Their unmasked '//&
        &'FSC determines the filter and FSC=0.143 cutoff. Writes <map>_lp, <map>_pproc, the simulated and docked '//&
        &'model densities, and a binary FSC curve.',&
        &'simple_exec',&
        &.false., visibility=UI_VIS_STANDARD, display_name='Post-process with Atomic-Model Envelope')
        call postprocess_pdb%add_input(UI_IMG, 'vol1', 'file', 'Volume', &
        &'Input density map with matching _even and _odd half-maps in the same directory', &
        &'e.g. map.mrc', .true., '', visibility=UI_VIS_STANDARD)
        call postprocess_pdb%add_input(UI_FILE, 'pdbfile', 'file', 'Atomic model', &
        &'PDB coordinates used to generate and dock the masking envelope', &
        &'e.g. model.pdb', .true., '', visibility=UI_VIS_STANDARD)
        call postprocess_pdb%add_input(UI_PARM, smpd, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_PARM, center_pdb, visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_FILE, 'fsc', 'file', 'FSC output', &
        &'Output binary FSC curve, recalculated from the half-maps; overwrites this path when supplied. '//&
        &'Otherwise writes <map>_fsc.bin in the execution directory', &
        &'e.g. map_fsc.bin', .false., '', visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_FILT, bfac, &
        &help_override='B-factor in Angstroms^2; when omitted, estimated by a Guinier fit if the FSC cutoff is '//&
        &'better than 5 Angstroms, otherwise zero', visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_PARM, mirr, &
        &help_override='Write an additional x-mirrored postprocessed map; mirr=no disables this output', &
        &visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_MASK, mskdiam, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call postprocess_pdb%add_input(UI_COMP, nthr, visibility=UI_VIS_STANDARD)
        call add_ui_program('postprocess_pdb', postprocess_pdb, prgtab, UI_CATEGORY)
    end subroutine new_postprocess_pdb

end module simple_ui_postprocess
