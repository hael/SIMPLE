!@descr: module defining the user interfaces for 3D analysis of nanoparticles in the single_exec suite
module single_ui_nano3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('nano3d', '3D Reconstruction', 40)
type(ui_program), target :: autorefine3D_nano
type(ui_program), target :: refine3D_nano
type(ui_program), target :: solve3D_nano
type(ui_program), target :: import_tilt_series
type(ui_program), target :: refine_tomo

contains

    subroutine construct_single_nano3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_solve3D_nano(prgtab)
        call new_autorefine3D_nano(prgtab)
        call new_refine3D_nano(prgtab)
        call new_import_tilt_series(prgtab)
        call new_refine_tomo(prgtab)
    end subroutine construct_single_nano3D_programs

subroutine new_solve3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve3D_nano%new(&
        &'solve3D_nano',&                                                                                  ! name
        &'De novo 3D map determination for a nanoparticle data set',&    ! summary
        &'is a wrapper around solve3D that applies nanoparticle-oriented defaults while allowing overrides',& ! help
        &'single_exec',&                                                                                   ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='De Novo Nanoparticle Map') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! search controls
        call solve3D_nano%add_input(UI_SRCH, nsample, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call solve3D_nano%add_input(UI_FILT, hp, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_nano%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', 'Starting low-pass limit', 'low-pass limit in Angstroms{3}', .false., 3., &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_nano%add_input(UI_FILT, 'lpstop',     'num', 'Final low-pass limit', 'Final low-pass limit', 'low-pass limit in Angstroms{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call solve3D_nano%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call solve3D_nano%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_nano%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('solve3D_nano', solve3D_nano, prgtab, UI_CATEGORY)
    end subroutine new_solve3D_nano

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
        call refine3D_nano%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', &
        &'Refinement mode; cont refines the stored poses continuously(shc|neigh|prob|prob_neigh|cont){neigh}', &
        &'', .false., 'neigh', &
        &choices=ui_choices([character(len=10) :: 'shc', 'neigh', 'prob', 'prob_neigh', 'cont']), &
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

    subroutine new_import_tilt_series( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_tilt_series%new(&
        &'import_tilt_series',&                                                                        ! name
        &'Add an atomic electron tomography tilt series and its angles to a SIMPLE project',&           ! summary
        &'is a program for importing a tilt series with per-projection angles in the SIMPLE, single-axis tilt or '//&
        &'GENFIRE/RESIRE convention, ready for refine_tomo',&                                          ! help
        &'single_exec',&                                                                               ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Import Tilt Series')                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call import_tilt_series%add_input(UI_IMG, 'stk', 'file', 'Tilt series', &
        &'Stack of projection images, image x along the first array dimension', 'e.g. projs.mrcs', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call import_tilt_series%add_input(UI_FILE, 'angtab', 'file', 'Projection angles', &
        &'Per-projection angles in the convention of angconv, one projection per line in stack order', &
        &'e.g. angles.txt', .true., '', visibility=UI_VIS_STANDARD)
        call import_tilt_series%add_input(UI_PARM, 'angconv', 'multi', 'Angle convention', &
        &'Convention of angtab: simple, an oritab (e1 e2 e3 x y); tilt, one single-axis tilt angle about image x per line, '//&
        &'optionally followed by x y; genfire, phi theta psi of GENFIRE/RESIRE per line, optionally followed by x y'//&
        &'(simple|tilt|genfire){simple}', '', .true., 'simple', &
        &choices=ui_choices([character(len=7) :: 'simple', 'tilt', 'genfire']), visibility=UI_VIS_STANDARD)
        call import_tilt_series%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! search controls
        call import_tilt_series%add_input(UI_SRCH, nspace, &
        &visibility=UI_VIS_ADVANCED)
        ! add to ui_hash
        call add_ui_program('import_tilt_series', import_tilt_series, prgtab, UI_CATEGORY)
    end subroutine new_import_tilt_series

    subroutine new_refine_tomo( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine_tomo%new(&
        &'refine_tomo',&                                                                                ! name
        &'Refine the projection poses of an atomic electron tomography tilt series',&                   ! summary
        &'is a workflow that reconstructs by gridding, masks the density, reconstructs by masked PCG and refines '//&
        &'every projection''s rotation and shift continuously in a short pass at a fixed low-pass limit',& ! help
        &'single_exec',&                                                                                ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Refine Tilt Series')                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! search controls
        call refine_tomo%add_input(UI_SRCH, 'trs', 'num', 'Shift bound', &
        &'Largest total shift a projection may move from its imported shift', 'in pixels{3}', .false., 3., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'athres_cont', 'num', 'Rotation bound', &
        &'Largest total rotation a projection may move from its imported pose', 'in degrees{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'maxits', 'num', 'Max iterations per stage', &
        &'Maximum number of refinement iterations in each low-pass stage', 'Max # iterations{2}', .false., 2., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'minits', 'num', 'Min iterations per stage', &
        &'Minimum number of refinement iterations in each low-pass stage', 'Min # iterations{0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'cont_conv_rot', 'num', 'Convergence rotation', &
        &'A projection is stable when it moved by at most this rotation since its last pose', 'in degrees{0.5}', &
        &.false., 0.5, visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'cont_conv_shift', 'num', 'Convergence shift', &
        &'A projection is stable when it moved by at most this shift since its last pose', 'in Angstroms{1.0}', &
        &.false., 1.0, visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_SRCH, 'cont_conv_frac', 'num', 'Convergence fraction', &
        &'A stage converges when this fraction of the projections is stable', '(0,1]{0.9}', .false., 0.9, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call refine_tomo%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Low-pass limit of the first refinement stage', 'low-pass limit in Angstroms{2}', .false., 2., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Low-pass limit of the last refinement stage; no stage goes beyond the half-map FSC 0.143 resolution', &
        &'low-pass limit in Angstroms{2}', .false., 2., visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_FILT, 'nstages', 'num', 'Number of stages', &
        &'Number of low-pass stages, geometrically spaced from lpstart to lpstop', '# stages{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        call refine_tomo%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG iterations', &
        &'Iterations of every masked PCG reconstruction', '# iterations{5}', .false., 5., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine_tomo%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        call refine_tomo%add_input(UI_MASK, 'amsklp', 'num', 'Envelope low-pass limit', &
        &'Low-pass limit of the map the density envelope (the PCG support) is drawn from', &
        &'low-pass limit in Angstroms{6}', .false., 6., visibility=UI_VIS_ADVANCED)
        ! computer controls
        call refine_tomo%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine_tomo', refine_tomo, prgtab, UI_CATEGORY)
    end subroutine new_refine_tomo

end module single_ui_nano3D
