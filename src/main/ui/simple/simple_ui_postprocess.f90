!@descr: module defining the user interfaces for map post-processing programs in the simple_exec suite
module simple_ui_postprocess
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('postprocess', 'Post-processing', 69)
type(ui_program), target :: postprocess
type(ui_program), target :: postprocess_nu

contains

    subroutine construct_postprocess_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_postprocess(prgtab)
        call new_postprocess_nu(prgtab)
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
        call postprocess%add_input(UI_PARM, 'imgkind', 'str', 'Volume kind', 'Project output volume kind{vol}', &
            &'project output kind: vol or vol_cavg', .false., 'vol', &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false., 20., &
        &visibility=UI_VIS_ADVANCED)
        call postprocess%add_input(UI_FILT, 'fsc', 'file', 'FSC file', 'Binary FSC file for optimal filtering', &
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
        &analogue). Feed the UNREGULARIZED even/odd half pair; outputs carry the _nu_sharp suffix and are &
        &display/interpretation maps only, never inputs to FSC correction or resolution estimation. Isolated from the &
        &standard postprocess program (global B-factor), which is unchanged',& ! help
        &'simple_exec',&                                       ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call postprocess_nu%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Unregularized odd half-map',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call postprocess_nu%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Unregularized even half-map', 'vol2.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call postprocess_nu%add_input(UI_IMG, outvol, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call postprocess_nu%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! filter controls
        call postprocess_nu%add_input(UI_FILT, nu_refine, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call postprocess_nu%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call postprocess_nu%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('postprocess_nu', postprocess_nu, prgtab, UI_CATEGORY)
    end subroutine new_postprocess_nu

end module simple_ui_postprocess
