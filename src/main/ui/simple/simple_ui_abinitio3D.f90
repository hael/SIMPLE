!@descr: module defining the user interfaces for ab initio 3D reconstruction programs in the simple_exec suite
module simple_ui_abinitio3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('abinitio3d', 'Ab Initio 3D Reconstruction', 50)
type(ui_program), target :: abinitio3D
type(ui_program), target :: abinitio3D_cavgs
type(ui_program), target :: estimate_lpstages
type(ui_program), target :: noisevol

contains

    subroutine construct_abinitio3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio3D(prgtab)
        call new_abinitio3D_cavgs(prgtab)
        call new_estimate_lpstages(prgtab)
        call new_noisevol(prgtab)
    end subroutine construct_abinitio3D_programs

subroutine new_abinitio3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D%new(&
        &'abinitio3D',&                                                                    ! name
        &'Generate an initial 3D model from particle images',&                             ! summary
        &'is a distributed workflow for generating an ab initio 3D model from particles',& ! help
        &'simple_exec',&                                                                   ! executable
        &.true.,&                                                                          ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Create Initial 3D Model')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call abinitio3D%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', 'Starting reference volume &
        & for particle matching', 'input starting volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call abinitio3D%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){no}','', .false., 'no', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, pgrp, group="model", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_SRCH, pgrp_start, group="model", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, 'cavg_ini', 'binary', '3D initialization on class averages', '3D initialization on class averages(yes|no){no}','', .false., 'no', group="model", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, 'cavg_ini_ext', 'binary', 'External class-average 3D initialization', &
            &'Use existing ptcl3D orientations and state assignments from a prior abinitio3D_cavgs run; skips the symmetry-search stage(yes|no){no}','', .false., 'no', group="model", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call abinitio3D%add_input(UI_SRCH, nsample, group="search", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_SRCH, 'nstages', 'num', 'Last ab initio stage to run',&
            &'Last abinitio3D stage to run; default is 5 for multivol_mode=independent and 8 otherwise; &
            &independent mode writes final volumes at its last stage',&
            &'last stage', .false., 8., group="search", visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, nstates, group="search", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_SRCH, 'state', 'num', 'Continuation state label', &
            &'State label to select from an existing multi-state abinitio3D project and continue as a single-state stage-5 search', &
            &'state label', .false., 1., group="search", visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, 'multivol_mode', 'multi', 'Multi-volume ab initio mode', 'Multi-volume ab initio mode(single|independent|docked){single}','', .false., 'single', &
        &choices=ui_choices([character(len=11) :: 'single', 'independent', 'docked']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, 'split_stage', 'num', 'Docked-state split stage', &
        &'Stage at which docked multi-volume ab-initio reconstruction splits into state-specific volumes', &
        &'stage number', .false., 6., group="search", visibility=UI_VIS_DEVELOPER)
        call abinitio3D%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of particle assignments for ab-initio convergence', 'overlap fraction', .false., .95, &
        &group="search", visibility=UI_VIS_DEVELOPER)
        call abinitio3D%add_input(UI_SRCH, objfun_den, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, objfun_den_w, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_SRCH, 'projrec', 'binary', 'Projection-direction reconstruction', &
            &'Assemble raw 2D Fourier numerator/CTF-squared sums by projection direction before compact 3D reconstruction(yes|no){no}','', .false., 'no', group="search", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        ! filter controls
        call abinitio3D%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms',  .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms; default is 6 for multivol_mode=independent &
            &and 8 otherwise',    .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, lp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'force_lp_range', 'binary', 'Force low-pass range', &
            &'Use lpstart/lpstop directly for abinitio3D low-pass stages instead of class-FRC-derived limits(yes|no){no}','', .false., 'no', group="filter", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call abinitio3D%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
            &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}; nonuniform_lpset promotes the &
            &NU frontier into an explicit merged-reference LP-set matching run','', .false., 'nonuniform', &
            &group="filter", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']))
        call abinitio3D%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, conical_fsc, group="filter", visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'lpstart_ini3D',  'num', 'Starting low-pass limit ini3D', 'Starting low-pass limit ini3D',&
            &'low-pass limit for the initial stage of ini3D in Angstroms',  .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D%add_input(UI_FILT, 'lpstop_ini3D',  'num', 'Final low-pass limit ini3D', 'Final low-pass limit ini3D',&
            &'low-pass limit for the final stage of ini3D in Angstroms',    .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call abinitio3D%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_MASK, 'automsk', 'multi', 'Perform envelope masking', &
            &'Whether to generate/apply an envelope mask from the staged automasking point(yes|tight|no){no}','', .false., 'no', group="mask", visibility=UI_VIS_STANDARD, &
        &choices=ui_choices([character(len=5) :: 'yes', 'tight', 'no']))
        call abinitio3D%add_input(UI_MASK, envref, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call abinitio3D%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_COMP, nthr,                              group="compute", visibility=UI_VIS_STANDARD)
        call abinitio3D%add_input(UI_COMP, 'nthr_ini3D', 'num', 'Number of threads for ini3D phase, give 0 if unsure', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .false., 0., group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('abinitio3D', abinitio3D, prgtab, UI_CATEGORY)
    end subroutine new_abinitio3D

    subroutine new_abinitio3D_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D_cavgs%new(&
        &'abinitio3D_cavgs',&                                                                   ! name
        &'Generate an initial 3D model from 2D class averages',&                               ! summary
        &'is a distributed workflow for generating an ab initio 3D model from class averages',& ! help
        &'simple_exec',&                                                                        ! executable
        &.true., visibility=UI_VIS_STANDARD, display_name='Create Initial 3D Model from Classes') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call abinitio3D_cavgs%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call abinitio3D_cavgs%add_input(UI_SRCH, pgrp_start, &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_SRCH, nstates, group="search", visibility=UI_VIS_STANDARD)
        call abinitio3D_cavgs%add_input(UI_SRCH, 'multivol_mode', 'multi', 'Multi-volume class-average ab initio mode', &
            &'Multi-volume class-average ab initio mode(single|independent|docked){single}','', .false., 'single', &
        &choices=ui_choices([character(len=11) :: 'single', 'independent', 'docked']), &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_SRCH, 'split_stage', 'num', 'Docked-state split stage', &
        &'Stage at which docked multi-volume class-average ab-initio reconstruction splits into state-specific volumes', &
        &'stage number', .false., 6., group="search", visibility=UI_VIS_DEVELOPER)
        call abinitio3D_cavgs%add_input(UI_SRCH, 'overlap', 'num', 'Convergence overlap target', &
        &'Required overlap of class-average assignments for ab-initio convergence', 'overlap fraction', .false., .95, &
        &group="search", visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call abinitio3D_cavgs%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms', .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call abinitio3D_cavgs%add_input(UI_FILT, conical_fsc, group="filter", visibility=UI_VIS_ADVANCED)
        ! mask controls
        call abinitio3D_cavgs%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        ! computer controls
        call abinitio3D_cavgs%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call abinitio3D_cavgs%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('abinitio3D_cavgs', abinitio3D_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_abinitio3D_cavgs

    subroutine new_estimate_lpstages( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call estimate_lpstages%new(&
        &'estimate_lpstages',&                                                                                             ! name
        &'Estimation of low-pass limits, shift boundaries, and downscaling parameters for ab initio 3D',&                  ! summary
        &'is a program for estimation of low-pass limits, shift boundaries, and downscaling parameters for ab initio 3D',& ! help
        &'simple_exec',&                                                                                                   ! executable
        &.true.)                                                                                                           ! requires sp_project
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
        &visibility=UI_VIS_ADVANCED)                             ! requires sp_project
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

end module simple_ui_abinitio3D
