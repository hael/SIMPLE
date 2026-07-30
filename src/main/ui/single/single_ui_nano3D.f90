!@descr: module defining the user interfaces for 3D analysis of nanoparticles in the single_exec suite
module single_ui_nano3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('nano3d', '3D Reconstruction', 40)
type(ui_program), target :: autorefine3D_nano
type(ui_program), target :: refine3D_nano
type(ui_program), target :: abinitio3D_nano

contains

    subroutine construct_single_nano3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio3D_nano(prgtab)
        call new_autorefine3D_nano(prgtab)
        call new_refine3D_nano(prgtab)
    end subroutine construct_single_nano3D_programs

subroutine new_abinitio3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D_nano%new(&
        &'abinitio3D_nano',&                                                                               ! name
        &'Generate an initial 3D model for a nanoparticle data set',&    ! summary
        &'is a wrapper around abinitio3D that applies nanoparticle-oriented defaults while allowing overrides',& ! help
        &'single_exec',&                                                                                   ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Create Initial Nanoparticle Model') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! search controls
        call abinitio3D_nano%add_input(UI_SRCH, nsample, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call abinitio3D_nano%add_input(UI_FILT, hp, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_nano%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', 'Starting low-pass limit', 'low-pass limit in Angstroms{3}', .false., 3., &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_nano%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit', 'low-pass limit in Angstroms{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call abinitio3D_nano%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call abinitio3D_nano%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_nano%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('abinitio3D_nano', abinitio3D_nano, prgtab, UI_CATEGORY)
    end subroutine new_abinitio3D_nano

    subroutine new_autorefine3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call autorefine3D_nano%new(&
        &'autorefine3D_nano',&                                                            ! name
        &'Automatically refine a metallic nanoparticle structure from particle images',&  ! summary
        &'is a distributed workflow for automated 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! help
        &'single_exec',&                                                                  ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Automated Nanoparticle Refinement') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call autorefine3D_nano%add_input(UI_IMG, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call autorefine3D_nano%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call autorefine3D_nano%add_input(UI_PARM, element, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call autorefine3D_nano%add_input(UI_SRCH, nspace, &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_SRCH, 'maxits', 'num', 'Max iterations', 'Maximum number of iterations', 'Max # iterations{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of particle assignments for nanoparticle-refinement convergence', 'overlap fraction', .false., .98, &
        &visibility=UI_VIS_DEVELOPER)
        call autorefine3D_nano%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call autorefine3D_nano%add_input(UI_SRCH, nrestarts, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call autorefine3D_nano%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call autorefine3D_nano%add_input(UI_FILT, 'lp', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms{1.5}', .true., 1.5, &
        &visibility=UI_VIS_STANDARD)
        call autorefine3D_nano%add_input(UI_FILT, icm, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call autorefine3D_nano%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        call autorefine3D_nano%add_input(UI_MASK, 'mskdiam_detect', 'num', 'Detect-atoms mask diameter', 'Optional mask diameter in Angstroms passed only to detect_atoms', 'mask diameter in Angstroms{0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call autorefine3D_nano%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        call autorefine3D_nano%add_input(UI_COMP, script, &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('autorefine3D_nano', autorefine3D_nano, prgtab, UI_CATEGORY)
    end subroutine new_autorefine3D_nano

    subroutine new_refine3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D_nano%new(&
        &'refine3D_nano',&                                                                                                    ! name
        &'Refine a metallic nanoparticle structure from particle images',&                                                     ! summary
        &'is a distributed workflow for 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! help
        &'single_exec',&                                                                                                      ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Refine Nanoparticle Structure') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D_nano%add_input(UI_IMG, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_IMG, 'vol_odd',  'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_IMG, 'vol_even', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call refine3D_nano%add_input(UI_SRCH, nspace, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, maxits, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, update_frac, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call refine3D_nano%add_input(UI_SRCH, sigma_est, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', &
        &'Refinement mode(snhc|shc|neigh|shc_neigh|prob|prob_state|prob_neigh){neigh}', '', .false., 'neigh', &
        &choices=ui_choices([character(len=10) :: 'snhc', 'shc', 'neigh', 'shc_neigh', 'prob', 'prob_state', 'prob_neigh']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call refine3D_nano%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_FILT, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms{1.0}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_FILT, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', &
        &'Low-pass limit used to limit the resolution during refinement', 'low-pass limit in Angstroms{0.5}', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_FILT, icm, &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', &
        &'Strength of noise reduction', '(0.01-3.0){0.1}', .false., 0.1, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D_nano%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call refine3D_nano%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_nano%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine3D_nano', refine3D_nano, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_nano

end module single_ui_nano3D
