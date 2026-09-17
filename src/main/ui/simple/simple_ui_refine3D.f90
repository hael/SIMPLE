!@descr: module defining the user interfaces for 3D refinement programs in the simple_exec suite
module simple_ui_refine3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('refine3d', 'Refine 3D Workflows', 60)
type(ui_program), target :: refine3D
type(ui_program), target :: refine3D_auto

contains

    subroutine construct_refine3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_refine3D(prgtab)
        call new_refine3D_auto(prgtab)
    end subroutine construct_refine3D_programs

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
        call refine3D%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend for per-iteration half-map assembly(gridding|pcg){gridding}', &
        &'', .false., 'gridding', group="search", &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_PARM, 'box_crop', 'num', 'Refinement box', &
        &'Even Fourier-cropped refinement box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0, group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_PARM, 'euclid_diag', 'binary', 'Euclid scale diagnostics', &
        &'Per-iteration report of the reference/particle amplitude ratio per band and the euclid objective quantiles(yes|no){no}','', .false., 'no', visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call refine3D%add_input(UI_PARM, 'projrec', 'binary', 'Projection-direction reconstruction',&
        &'Assemble raw 2D Fourier numerator/CTF-squared sums by projection direction before compact 3D reconstruction(yes|no){no}','', .false., 'no', visibility=UI_VIS_ADVANCED, &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
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
        call refine3D%add_input(UI_SRCH, 'pose_cont', 'binary', &
        &'Experimental five-parameter pose refinement', &
        &'Run transactional Cartesian LM after the established matcher; mutually exclusive with '// &
        &'refine=pose_cont(yes|no){no}', '', &
        &.false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'pose_cont_route', 'multi', &
        &'Continuous pose LM route', &
        &'LM route used by pose_cont=yes or '// &
        &'refine=pose_cont(shift_then_joint|joint){shift_then_joint}', '', &
        &.false., 'shift_then_joint', group="search", &
        &choices=ui_choices([character(len=16) :: 'shift_then_joint', 'joint']), &
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
        call refine3D%add_input(UI_SRCH, 'inpl_cont', 'binary', &
        &'Continuous in-plane refinement', &
        &'Joint continuous Euclidean in-plane and shift refinement(yes|no){yes}', '', &
        &.false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun_den, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun_den_w, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', &
        &'Refinement mode(snhc|shc|neigh|shc_neigh|prob|prob_state|prob_neigh|pose_cont){shc}','',&
        &.false., 'shc', group="search", &
        &choices=ui_choices([character(len=10) :: 'snhc', 'shc', 'neigh', 'shc_neigh', &
        &'prob', 'prob_state', 'prob_neigh', 'pose_cont']), &
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
        call refine3D%add_input(UI_FILT, envmsklp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'amsklp', 'num', 'NU evidence envelope smoothing limit',&
        & 'Low-pass limit for NU evidence envelope generation in Angstroms', &
        &'low-pass limit in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, ml_reg, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations; independent of refine3D outer maxits', 'iterations{2}', &
        &.false., 2., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D%add_input(UI_FILT, 'pcg_mskfile', 'file', 'PCG support-constraint mask volume', &
        &'Real-space [0,1] mask volume installed as the hard support constraint of every PCG solve (the projected '//&
        &'system P H P; experimental focused/support mode); spherical mskdiam support when absent', &
        &'e.g. focusmask.mrc', .false., '', group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D%add_input(UI_FILT, conical_fsc, group="filter", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_FILT, combine_eo, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call refine3D%add_input(UI_MASK, automsk_refine3D, group="mask", &
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
        call refine3D_auto%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend for per-iteration half-map assembly, forwarded to the refine3D child and the '//&
        &'bootstrap/final reconstructions(gridding|pcg){gridding}', &
        &'', .false., 'gridding', group="search", &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), visibility=UI_VIS_ADVANCED)
        ! search controls
        call refine3D_auto%add_input(UI_SRCH, maxits,      required_override=.false., group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'minits', 'num', 'Minimum automatic iterations', &
        &'Minimum number of iterations used when automatic refinement calculates its iteration budget', &
        &'minimum iterations', .false., 10., group="search", visibility=UI_VIS_DEVELOPER)
        call refine3D_auto%add_input(UI_SRCH, 'nsample', 'num', 'Projection samples', &
        &'Number of projection samples used by automatic refinement', 'number of samples', .false., 25000., &
        &group="search", visibility=UI_VIS_DEVELOPER)
        call refine3D_auto%add_input(UI_SRCH, 'ref_pose_init', 'multi', 'External-reference pose initialization', &
        &'When vol1 has independent provenance, run one fixed-reference CC pose-initialization pass at 15 Angstroms &
        &before Euclidean refinement(cc|none){none}', '', .false., 'none', group='search', &
        &choices=ui_choices([character(len=4) :: 'cc', 'none']), visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, pgrp,                                  group="search", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_SRCH, ptcl_src, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, sigma_est, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'inpl_cont', 'binary', &
        &'Continuous in-plane refinement', &
        &'Joint continuous Euclidean in-plane and shift refinement(yes|no){yes}', '', &
        &.false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', .false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'regpass', 'binary', 'Global registration pass', &
        &'One global (refine=prob) registration pass of all particles against the masked startup references, &
        &band-limited at the FSC=regpass_fsc resolution of the startup pair, before the neighbourhood iterations(yes|no){yes}', &
        &'', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_SRCH, 'regpass_fsc', 'num', 'Registration-pass FSC criterion', &
        &'FSC value of the startup pair whose resolution band-limits the registration pass', &
        &'FSC value in (0,1){0.143}', .false., 0.143, group="search", &
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
        call refine3D_auto%add_input(UI_FILT, 'amsklp', 'num', 'NU evidence envelope smoothing limit',&
        & 'Low-pass limit for NU evidence envelope generation in Angstroms', &
        &'low-pass limit in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}','', .false., 'nonuniform', group="filter", &
        &choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations during refinement; the cold original-sampling final reconstruction uses at least 5', &
        &'iterations{2}', &
        &.false., 2., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D_auto%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D_auto%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, envmsklp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, combine_eo, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_FILT, 'res_target', 'num', 'Resolution target (in A)',&
        & 'Resolution target in Angstroms', 'Resolution target in Angstroms', .false., 3., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call refine3D_auto%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_MASK, automsk_refine3D, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D_auto%add_input(UI_MASK, nu_msk_sig, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call refine3D_auto%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call refine3D_auto%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('refine3D_auto', refine3D_auto, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_auto




end module simple_ui_refine3D
