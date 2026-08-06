!@descr: module defining the user interfaces for 3D refinement programs in the simple_exec suite
module simple_ui_refine3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('refine3d', 'Refine 3D Workflows', 60)
type(ui_program), target :: automask
type(ui_program), target :: postprocess
type(ui_program), target :: reconstruct3D
type(ui_program), target :: bootstrap_rec3D
type(ui_program), target :: refine3D
type(ui_program), target :: refine3D_auto
type(ui_program), target :: refine3D_multi

contains

    subroutine construct_refine3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_automask(prgtab)
        call new_postprocess(prgtab)
        call new_reconstruct3D(prgtab)
        call new_bootstrap_rec3D(prgtab)
        call new_refine3D(prgtab)
        call new_refine3D_auto(prgtab)
        call new_refine3D_multi(prgtab)
    end subroutine construct_refine3D_programs

subroutine new_automask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call automask%new(&
        &'automask',&                                    ! name
        &'Create a spherical mask from the estimated particle diameter',& ! summary
        &'is a program for automated envelope masking',& ! help
        &'simple_exec',&                                 ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call automask%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call automask%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 8., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call automask%add_input(UI_MASK, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        call automask%add_input(UI_MASK, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for envelope mask generation', 'Volume threshold, give 0 if unknown', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call automask%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6., &
        &visibility=UI_VIS_ADVANCED)
        call automask%add_input(UI_MASK, automsk, &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call automask%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('automask', automask, prgtab, UI_CATEGORY)
    end subroutine new_automask

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

    subroutine new_reconstruct3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call reconstruct3D%new(&
        &'reconstruct3D',&                                               ! name
        &'3D reconstruction from oriented particles',&                   ! summary
        &'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,&
        & given input orientations and state assignments. The algorithm is based on direct Fourier inversion&
        & with a Kaiser-Bessel (KB) interpolation kernel',&
        &'simple_exec',&                                                 ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call reconstruct3D%add_input(UI_SRCH, trs, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call reconstruct3D%add_input(UI_SRCH, ptcl_src, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_SRCH, 'projrec', 'binary', 'Projection-direction reconstruction',&
        &'Assemble raw 2D Fourier numerator/CTF-squared sums by projection direction before compact 3D reconstruction(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call reconstruct3D%add_input(UI_FILT, envfsc, &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map',&
        &'Postprocess reconstructed volumes using the generated FSC curves','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call reconstruct3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call reconstruct3D%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call reconstruct3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('reconstruct3D', reconstruct3D, prgtab, UI_CATEGORY)
    end subroutine new_reconstruct3D

    subroutine new_bootstrap_rec3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call bootstrap_rec3D%new(&
        &'bootstrap_rec3D',&                                             ! name
        &'bootstrap ML-regularized 3D reconstruction',&                  ! summary
        &'generates an unregularized even/odd reconstruction, estimates weighted global sigma2 curves from the half-map difference,&
        & and reruns reconstruct3D with Euclidean ML regularization',&
        &'simple_exec',&                                                 ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call bootstrap_rec3D%add_input(UI_PARM, 'which_iter', 'num', 'Sigma iteration index',&
        &'Iteration index used for the generated sigma2_groups file{1}', 'iteration{1}', .false., 1.0, &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_FILE, 'outfile', 'file', 'Resolution output prefix',&
        &'Optional FSC/resolution text output prefix; state tags are appended', 'e.g. resolution',&
        &.false., 'RESOLUTION_FINAL.txt', &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call bootstrap_rec3D%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        call bootstrap_rec3D%add_input(UI_SRCH, nstates, &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call bootstrap_rec3D%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map',&
        &'Postprocess ML-regularized reconstructed volumes using the generated FSC curves','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_FILT, combine_eo, &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        call bootstrap_rec3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call bootstrap_rec3D%add_input(UI_COMP, nparts, required_override=.false., &
        &visibility=UI_VIS_DEVELOPER)
        call bootstrap_rec3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('bootstrap_rec3D', bootstrap_rec3D, prgtab, UI_CATEGORY)
    end subroutine new_bootstrap_rec3D

    subroutine new_refine3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',&                                                                               ! name
        &'Refine a 3D structure from particle images by projection matching',& ! summary
        &'is a distributed workflow for 3D refinement based on probabilistic projection matching',& ! help
        &'simple_exec',&                                                                            ! executable
        &.true.,&                                                                                   ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Refine 3D Structure')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%add_input(UI_IMG, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. vol.mrc', .false., 'vol1.mrc', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call refine3D%add_input(UI_SRCH, nspace, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, trs, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}','', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, maxits, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, update_frac, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, pgrp, group="search", visibility=UI_VIS_STANDARD)
        call refine3D%add_input(UI_SRCH, nstates, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun_den, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun_den_w, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'projrec', 'binary', 'Projection-direction reconstruction',&
        &'Assemble raw 2D Fourier numerator/CTF-squared sums by projection direction before compact 3D reconstruction(yes|no){no}','', .false., 'no', group="search", visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call refine3D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc|shc|neigh|shc_neigh|prob|prob_state|prob_neigh){shc}','',&
        &.false., 'shc', group="search", &
        &choices=ui_choices([character(len=10) :: 'snhc', 'shc', 'neigh', 'shc_neigh', 'prob', 'prob_state', 'prob_neigh']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'prob_neigh_mode', 'multi', 'Prob-neigh neighborhood mode', &
        &'Prob-neigh neighborhood mode(state|geom|shc|snhc){state}','', .false., 'state', &
        &group="search", &
        &choices=ui_choices([character(len=5) :: 'state', 'geom', 'shc', 'snhc']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}','', .false.,&
        &'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, sigma_est, group="search", &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call refine3D%add_input(UI_FILT, hp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, lplim_crit, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, lp_backgr, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(none|uniform|fsc|nonuniform|nonuniform_lpset){none}','', .false., 'none', group="filter", &
        &choices=ui_choices([character(len=16) :: 'none', 'uniform', 'fsc', 'nonuniform', 'nonuniform_lpset']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, ml_reg, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, conical_fsc, group="filter", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, nu_refine, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, combine_eo, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call refine3D%add_input(UI_MASK, automsk, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_MASK, nu_msk_sig, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call refine3D%add_input(UI_COMP, nparts, required_override=.false., group="compute", visibility=UI_VIS_STANDARD)
        call refine3D%add_input(UI_COMP, nthr,                              group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine3D', refine3D, prgtab, UI_CATEGORY)
    end subroutine new_refine3D

    subroutine new_refine3D_auto( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D_auto%new(&
        &'refine3D_auto',&                                                                          ! name
        &'Automatically refine a single 3D structure from particle images',&                         ! summary
        &'is an automated workflow for single-state 3D refinement based on probabilistic projection matching',& ! help
        &'simple_exec',&                                                                            ! executable
        &.true.,&                                                                                   ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Automated 3D Refinement')
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D_auto%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', 'Starting reference volume &
        & for particle matching', 'input starting volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call refine3D_auto%add_input(UI_SRCH, maxits,      required_override=.false., group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'minits', 'num', 'Minimum automatic iterations', &
        &'Minimum number of iterations used when automatic refinement calculates its iteration budget', &
        &'minimum iterations', .false., 10., group="search", visibility=UI_VIS_DEVELOPER)
        call refine3D_auto%add_input(UI_SRCH, 'nsample', 'num', 'Projection samples', &
        &'Number of projection samples used by automatic refinement', 'number of samples', .false., 25000., &
        &group="search", visibility=UI_VIS_DEVELOPER)
        call refine3D_auto%add_input(UI_SRCH, pgrp,                                  group="search", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, sigma_est, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', .false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}','', .false.,&
        &'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call refine3D_auto%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}','', .false., 'nonuniform', group="filter", &
        &choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'nu_refine', 'binary', 'NU resolution expansion refinement', &
        & 'Allow one high-resolution nonuniform-filter bank expansion per refinement iteration(yes|no){yes}','', .false., 'yes', group="filter", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, combine_eo, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'res_target', 'num', 'Resolution target (in A)',&
        & 'Resolution target in Angstroms', 'Resolution target in Angstroms', .false., 3., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D_auto%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_MASK, automsk, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_MASK, nu_msk_sig, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call refine3D_auto%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine3D_auto', refine3D_auto, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_auto

    subroutine new_refine3D_multi( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D_multi%new(&
        &'refine3D_multi',&                                                                         ! name
        &'Refine multiple 3D states from particle images or project state labels',&                 ! summary
        &'is an automated workflow for multi-state 3D refinement from project state labels or command-line nstates',& ! help
        &'simple_exec',&                                                                            ! executable
        &.true.,&                                                                                   ! requires sp_project
        &visibility=UI_VIS_STANDARD, display_name='Multi-state 3D Refinement')
        ! search controls
        call refine3D_multi%add_input(UI_SRCH, maxits,      required_override=.false., group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, nstates,     required_override=.false., group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'nsample', 'num', 'Particle sample target', &
        &'Particles sampled per iteration; set 0 to derive the automatic target from the number of states', &
        &'particles (0=automatic)', .false., 0., group="search", visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call refine3D_multi%add_input(UI_SRCH, sigma_est, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'multivol_mode', 'multi', 'Multi-volume refinement mode', &
        &'Multi-volume refinement mode(input_oris_refine|input_oris_fixed){input_oris_refine}','', .false., 'input_oris_refine', &
        &group="search", &
        &choices=ui_choices([character(len=17) :: 'input_oris_refine', 'input_oris_fixed']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'prob_neigh_mode', 'multi', 'Prob-neigh neighborhood mode', &
        &'Prob-neigh neighborhood mode(state|geom){geom}','', .false., 'geom', &
        &group="search", &
        &choices=ui_choices([character(len=5) :: 'state', 'geom']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, pgrp,                                  group="search", visibility=UI_VIS_STANDARD)
        call refine3D_multi%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', .false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', &
        &'Continue previous refinement(yes|no){no}','', .false.,&
        &'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_SRCH, 'overlap', 'num', 'State-overlap convergence target', &
        &'Required overlap of state assignments for convergence in the probabilistic-neighborhood stage', &
        &'overlap fraction', .false., .99, group="search", visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call refine3D_multi%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(fsc|nonuniform_lpset|none){nonuniform_lpset}','', .false., 'nonuniform_lpset', group="filter", &
        &choices=ui_choices([character(len=16) :: 'fsc', 'nonuniform_lpset', 'none']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_FILT, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', &
        &'Low-pass limit used to limit the resolution during refinement', 'low-pass limit in Angstroms', .false., 6., &
        &group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_FILT, ml_reg, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D_multi%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call refine3D_multi%add_input(UI_MASK, 'automsk', 'multi', 'Perform envelope masking', &
        &'Generate/apply the NU-evidence envelope; requires filt_mode=nonuniform_lpset(yes|tight|no){no}','', .false., 'no', group="mask", &
        &choices=ui_choices([character(len=5) :: 'yes', 'tight', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_multi%add_input(UI_MASK, nu_msk_sig, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call refine3D_multi%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call refine3D_multi%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine3D_multi', refine3D_multi, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_multi

end module simple_ui_refine3D
