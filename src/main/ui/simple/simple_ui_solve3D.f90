!@descr: module defining the user interfaces for the solve3D de novo 3D map determination programs in the simple_exec suite
module simple_ui_solve3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('solve3d', 'De Novo 3D Map Determination', 50)
type(ui_program), target :: solve3D
type(ui_program), target :: solve3D_cavgs
type(ui_program), target :: solve3D_addon
type(ui_program), target :: estimate_lpstages
type(ui_program), target :: noisevol

contains

    subroutine construct_solve3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_solve3D(prgtab)
        call new_solve3D_cavgs(prgtab)
        call new_solve3D_addon(prgtab)
        call new_estimate_lpstages(prgtab)
        call new_noisevol(prgtab)
    end subroutine construct_solve3D_programs

    subroutine new_solve3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve3D%new(&
        &'solve3D',&                                                                       ! name
        &'De novo 3D map determination from particle images',&                             ! summary
        &'is a distributed workflow for de novo map determination from particles that '//&
        &'couples ab initio 3D reconstruction with initial 3D refinement',&                ! help
        &'simple_exec',&                                                                   ! executable
        &.true.,&                                                                          ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='De Novo 3D Map Determination')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call solve3D%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', 'Starting reference volume &
        & for particle matching', 'input starting volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call solve3D%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend from stage 3 onward; stages 1 and 2 always use gridding(gridding|pcg){gridding}', &
        &'', .false., 'gridding', group="search", &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_PARM, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations from stage 3 onward; the cold original-sampling final reconstruction uses at least 5', &
        &'iterations{2}', .false., 2., group="search", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call solve3D%add_input(UI_PARM, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., group="search", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call solve3D%add_input(UI_PARM, 'pcg_solvent', 'binary', 'PCG soft solvent prior', &
        &'Soft solvent prior on the PCG base solve: a real-space ridge pulling solvent toward zero, solvent identified '//&
        &'per half from a prior-free solve of that half (smoothed absolute density, Otsu, logistic weight), then the '//&
        &'same cold solve again with the ridge; half-independent, so the pair stays gold standard; the support is untouched; '//&
        &'active in solve3D from stage 7 (one stage after NU filtering starts)(yes|no){no}', '', .false., 'no', group="search", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call solve3D%add_input(UI_PARM, 'pcg_solvent_lambda', 'num', 'PCG solvent prior strength', &
        &'Ridge coefficient of the solvent prior relative to the data scale (1 = as strong as the low-band data term); '//&
        &'not given = estimated per state and iteration by cross-validation of the prior-free half pair', 'coefficient{auto}', &
        &.false., 1.0, group="search", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call solve3D%add_input(UI_PARM, 'pcg_solvent_check', 'binary', 'PCG solvent prior strength check', &
        &'Validation of the automatic solvent-prior strength: re-solve the whole strength grid for real and print the re-solve objective and residuals beside the closed-form estimate; eight extra pair solves per state and iteration, no effect on the result(yes|no){no}', '', .false., 'no', group="search", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call solve3D%add_input(UI_PARM, 'cavg_ini', 'binary', '3D initialization on class averages', '3D initialization on class averages(yes|no){no}','', .false., 'no', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_PARM, 'cavg_ini_ext', 'binary', 'External class-average 3D initialization', &
            &'Use existing ptcl3D orientations and state assignments from a prior solve3D_cavgs run; skips the symmetry-search stage(yes|no){no}','', .false., 'no', group="model", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call solve3D%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){no}','', .false., 'no', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_SRCH, pgrp, group="model", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_SRCH, pgrp_start, group="model", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_SRCH, nsample, group="search", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_SRCH, 'balance', 'multi', 'Fractional-update sampling units', &
        &'Units every fractional-update sample is balanced over: none draws globally from the lowest update-count '//&
        &'particles; class gives every selected 2D class the same share; cavg first gives every group of similar '//&
        &'class averages the same share, then every class inside a group, so a preferred view spread over many '//&
        &'classes no longer dominates the sample and classes of one view keep equal footing(none|class|cavg){cavg}', &
        &'', .false., 'cavg', group="search", choices=ui_choices([character(len=5) :: 'none', 'class', 'cavg']), &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_SRCH, 'nclust', 'num', 'Number of class-average groups', &
        &'Number of groups of similar class averages formed with balance=cavg, by average linkage on their aligned '//&
        &'correlation; with no more selected classes than this, every class is its own group{20}', '# groups{20}', &
        &.false., 20., group="search", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('balance', [character(len=4) :: 'cavg']))
        call solve3D%add_input(UI_SRCH, 'nstages', 'num', 'Last solve3D stage to run',&
            &'Last solve3D stage to run; default is 5 for nstates>1 and 8 otherwise; &
            &a multi-state run writes final volumes at its last stage',&
            &'last stage', .false., 8., group="search", visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_SRCH, nstates, group="search", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_SRCH, 'state', 'num', 'Continuation state label', &
            &'State label to select from an existing multi-state solve3D project and continue as a single-state stage-5 search', &
            &'state label', .false., 1., group="search", visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of particle assignments for solve3D stage convergence', 'overlap fraction', .false., .95, &
        &group="search", visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call solve3D%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'lpstart',     'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms',  .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'lpstop',     'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms; default is 6 for nstates>1 &
            &and 8 otherwise',    .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, lp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'force_lp_range', 'binary', 'Force low-pass range', &
            &'Use lpstart/lpstop directly for solve3D low-pass stages instead of class-FRC-derived limits(yes|no){no}','', .false., 'no', group="filter", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call solve3D%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
            &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}; nonuniform_lpset promotes the &
            &NU frontier into an explicit merged-reference LP-set matching run','', .false., 'nonuniform', &
            &group="filter", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']))
        call solve3D%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, envmsklp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'lpstart_ini3D',     'num', 'Starting low-pass limit ini3D', 'Starting low-pass limit ini3D',&
            &'low-pass limit for the initial stage of ini3D in Angstroms',  .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D%add_input(UI_FILT, 'lpstop_ini3D',     'num', 'Final low-pass limit ini3D', 'Final low-pass limit ini3D',&
            &'low-pass limit for the final stage of ini3D in Angstroms',    .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call solve3D%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_MASK, 'automsk', 'multi', 'Refinement envelope mode', &
            &'Use the density envelope, or prefer the lag-one NU-evidence envelope with density fallback, '//&
            &'from the staged automasking point(yes|nu|no){no}', &
            &'', .false., 'no', group="mask", visibility=UI_VIS_STANDARD, &
        &choices=ui_choices([character(len=3) :: 'yes', 'nu', 'no']))
        ! computer controls
        call solve3D%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_COMP, nthr,                                 group="compute", visibility=UI_VIS_STANDARD)
        call solve3D%add_input(UI_COMP, 'nthr_ini3D', 'num', 'Number of threads for ini3D phase, give 0 if unsure', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .false., 0., group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('solve3D', solve3D, prgtab, UI_CATEGORY)
    end subroutine new_solve3D

    !> Grow a completed solve3D solution with the particles a superset
    !! project adds: the frozen particles contribute their signal, unsearched,
    !! to every reconstruction; the others are searched against the union from
    !! stage 3 to the base run's last stage. Every setting that describes the
    !! solution comes from the frozen project's run manifest; the command line
    !! carries only compute effort, convergence and diagnostics.
    subroutine new_solve3D_addon( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve3D_addon%new(&
        &'solve3D_addon',&                                                                     ! name
        &'Extend a solve3D solution with the particles of a superset project',&                ! summary
        &'is a distributed workflow that searches the particles a superset project adds '//&
        &'against a frozen solve3D solution whose particles contribute unsearched; '//&
        &'when the run completes, its project replaces the superset project file',&            ! help
        &'simple_exec',&                                                                       ! executable
        &.true., visibility=UI_VIS_ADVANCED, display_name='Extend De Novo 3D Map')             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call solve3D_addon%add_input(UI_FILE, projfile, &
        &label_override       = 'Superset project', &
        &help_override        = 'SIMPLE project file whose particles are a superset of the frozen set in projfile_frozen; '//&
        &'the particles it adds are searched and, when the run completes, it is replaced by the run''s project', &
        &required_override    = .true., visibility=UI_VIS_STANDARD)
        call solve3D_addon%add_input(UI_FILE, 'projfile_frozen', 'file', 'Frozen solution project', &
        &'Project of a completed solve3D run (its run directory copy) whose particles are frozen', &
        &'e.g. 1_solve3D/myproject.simple', .true., '')
        ! parameter input/output
        call solve3D_addon%add_input(UI_PARM, 'addon_diag', 'binary', 'Cohort-only diagnostic map', &
        &'Also reconstruct the searched particles alone, without the frozen term, into addon_diag/(yes|no){no}', &
        &'', .false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call solve3D_addon%add_input(UI_PARM, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations (PCG base runs); default is the base run''s value', &
        &'iterations', .false., 2., group="search", visibility=UI_VIS_ADVANCED)
        call solve3D_addon%add_input(UI_PARM, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system (PCG base runs); default is the base run''s value', &
        &'iterations', .false., 0., group="search", visibility=UI_VIS_ADVANCED)
        call solve3D_addon%add_input(UI_PARM, 'pcg_solvent_check', 'binary', 'PCG solvent prior strength check', &
        &'Validation of the automatic solvent-prior strength when the base run used the solvent prior(yes|no){no}', &
        &'', .false., 'no', group="search", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        ! search controls
        call solve3D_addon%add_input(UI_SRCH, nsample, group="search", visibility=UI_VIS_STANDARD)
        call solve3D_addon%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of the searched particles'' assignments for early stopping in stage 3{0.95}', &
        &'overlap fraction', .false., .95, group="search", visibility=UI_VIS_ADVANCED)
        ! overridable settings of the frozen run: how the added particles are sampled and masked
        call solve3D_addon%add_input(UI_SRCH, 'balance', 'multi', 'Fractional-update sampling units', &
        &'Units every fractional-update sample of the searched particles is balanced over (see solve3D); '//&
        &'default is the frozen run''s value; none needs no 2D solution in the superset project, class and cavg '//&
        &'need its 2D classes(none|class|cavg)', &
        &'', .false., 'cavg', group="search", choices=ui_choices([character(len=5) :: 'none', 'class', 'cavg']), &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_addon%add_input(UI_SRCH, 'nclust', 'num', 'Number of class-average groups', &
        &'Number of groups of similar class averages formed with balance=cavg; default is the frozen run''s value', &
        &'# groups', .false., 20., group="search", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('balance', [character(len=4) :: 'cavg']))
        ! mask controls
        call solve3D_addon%add_input(UI_MASK, mskdiam, required_override=.false., &
        &help_override='Mask diameter in A of the searched particles and the union maps; default is the frozen run''s '//&
        &'value (the frozen accumulators are mask-free, so it may differ from the frozen run''s)', &
        &group="mask", visibility=UI_VIS_ADVANCED)
        ! computer controls
        call solve3D_addon%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call solve3D_addon%add_input(UI_COMP, nthr,                                 group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('solve3D_addon', solve3D_addon, prgtab, UI_CATEGORY)
    end subroutine new_solve3D_addon

    subroutine new_solve3D_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call solve3D_cavgs%new(&
        &'solve3D_cavgs',&                                                                      ! name
        &'De novo 3D map determination from 2D class averages',&                                ! summary
        &'is a distributed workflow for de novo map determination from class averages that '//&
        &'couples ab initio 3D reconstruction with initial 3D refinement',&                     ! help
        &'simple_exec',&                                                                        ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='De Novo 3D Map from Class Averages') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call solve3D_cavgs%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_cavgs%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call solve3D_cavgs%add_input(UI_SRCH, pgrp_start, &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_cavgs%add_input(UI_SRCH, nstates, group="search", visibility=UI_VIS_STANDARD)
        call solve3D_cavgs%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of class-average assignments for solve3D stage convergence', 'overlap fraction', .false., .95, &
        &group="search", visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call solve3D_cavgs%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_cavgs%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_cavgs%add_input(UI_FILT, 'lpstart',     'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms', .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call solve3D_cavgs%add_input(UI_FILT, 'lpstop',     'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call solve3D_cavgs%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        ! computer controls
        call solve3D_cavgs%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call solve3D_cavgs%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('solve3D_cavgs', solve3D_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_solve3D_cavgs

    subroutine new_estimate_lpstages( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call estimate_lpstages%new(&
        &'estimate_lpstages',&                                                                                             ! name
        &'Estimation of low-pass limits, shift boundaries, and downscaling parameters for solve3D',&                       ! summary
        &'is a program for estimation of low-pass limits, shift boundaries, and downscaling parameters for solve3D',&      ! help
        &'simple_exec',&                                                                                                   ! executable
        &.true., &
        &display_name='Estimate solve3D Stages') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call estimate_lpstages%add_input(UI_FILE, projfile, &
        &visibility=UI_VIS_STANDARD)
        call estimate_lpstages%add_input(UI_PARM, 'nstages', 'num', 'Number of low-pass limit stages', 'Number of low-pass limit stages', '# stages', .true., 8., &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('estimate_lpstages', estimate_lpstages, prgtab, UI_CATEGORY)
    end subroutine new_estimate_lpstages

    subroutine new_noisevol( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call noisevol%new(&
        &'noisevol',&                         ! name
        &'Generate one or more white-noise volumes',& ! summary
        &'is a program for generating noise volume(s)',&
        &'simple_exec',&                      ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED, display_name='Generate Noise Volumes')                             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call noisevol%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call noisevol%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        call noisevol%add_input(UI_PARM, 'nstates', 'num', 'Number states', 'Number states', '# states', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('noisevol', noisevol, prgtab, UI_CATEGORY)
    end subroutine new_noisevol

end module simple_ui_solve3D
