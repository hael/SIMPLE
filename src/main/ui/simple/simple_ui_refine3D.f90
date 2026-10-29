!@descr: module defining the user interfaces for 3D refinement programs in the simple_exec suite
module simple_ui_refine3D
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('refine3d', 'Refine 3D Workflows', 60)
type(ui_program), target :: refine3D
type(ui_program), target :: refine3D_auto
type(ui_program), target :: refine3D_local

contains

    subroutine construct_refine3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_refine3D(prgtab)
        call new_refine3D_auto(prgtab)
        call new_refine3D_local(prgtab)
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
        &visibility=UI_VIS_ADVANCED, display_name='Refine 3D Structure')
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
        call refine3D%add_input(UI_SRCH, 'topk', 'num', 'Poses per particle in the reconstruction', &
        &'Probabilistic modes (refine=prob|prob_neigh): the poses a particle enters the reconstruction at, its committed '//&
        &'draw and up to topk-1 extras sampled without replacement from the refined entries of its row, each with its '//&
        &'share of the tempered posterior over the sampled set (an extra below 1% is not inserted); 1 is the committed '//&
        &'pose alone, the search unchanged; gridding backend{1}', 'poses{1}', .false., 1., group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'shift_prior', 'binary', 'Shift prior in the probability tables', &
        &'The tables'' shift refinements and seed estimates carry a Gaussian prior on the shift''s distance from the '//&
        &'seed, with the previous iteration''s shift variance as its width and the particle''s temperature converting '//&
        &'it to loss units; the row holds the MAP loss, so a refinement does not wander onto a noise optimum far '//&
        &'from its seed(yes|no){yes}', '', .false., 'yes', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'pose_cont', 'binary', &
        &'Continuous pose polish', &
        &'Follow every discrete iteration by a continuous Cartesian pose pass over the same '// &
        &'particles; not with refine=cont(yes|no){no}', '', &
        &.false., 'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'athres_cont', 'num', 'Continuous pose rotation bound', &
        &'Halfwidth of the total rotation a continuous Cartesian pose solve (refine=cont, pose_cont=yes) may '// &
        &'move a particle from its seed: the capture range for wrongly assigned orientations', &
        &'in degrees{10}', .false., 10., group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'cont_conv_rot', 'num', 'Continuous pose convergence rotation', &
        &'A refine=cont particle is stable when it moved by at most this rotation since its last pose', &
        &'in degrees{0.5}', .false., 0.5, group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'cont_conv_shift', 'num', 'Continuous pose convergence shift', &
        &'A refine=cont particle is stable when it moved by at most this shift since its last pose', &
        &'in Angstroms{1.0}', .false., 1.0, group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'cont_conv_frac', 'num', 'Continuous pose convergence fraction', &
        &'refine=cont converges when this fraction of every state is stable, over iterations that together '// &
        &'sampled this fraction of the particles', '(0,1]{0.9}', .false., 0.9, group="search", visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, maxits, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, update_frac, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, pgrp, group="search", visibility=UI_VIS_STANDARD)
        call refine3D%add_input(UI_SRCH, nstates, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, objfun, group="search", &
        &visibility=UI_VIS_ADVANCED)
        call refine3D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', &
        &'Refinement mode; cont refines the stored poses continuously(shc|neigh|prob|prob_neigh|cont){shc}','',&
        &.false., 'shc', group="search", &
        &choices=ui_choices([character(len=10) :: 'shc', 'neigh', &
        &'prob', 'prob_neigh', 'cont']), &
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
        call refine3D%add_input(UI_FILT, 'pcg_solvent', 'binary', 'PCG soft solvent prior', &
        &'Soft solvent prior on the PCG base solve: a real-space ridge pulling solvent toward zero, solvent identified '//&
        &'per half from a prior-free solve of that half (smoothed absolute density, Otsu, logistic weight), then the '//&
        &'same cold solve again with the ridge; half-independent, so the pair stays gold standard; the support is untouched(yes|no){no}', '', .false., 'no', group="filter", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call refine3D%add_input(UI_FILT, 'pcg_solvent_lambda', 'num', 'PCG solvent prior strength', &
        &'Ridge coefficient of the solvent prior relative to the data scale (1 = as strong as the low-band data term); '//&
        &'not given = estimated per state and iteration by cross-validation of the prior-free half pair', 'coefficient{auto}', &
        &.false., 1.0, group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call refine3D%add_input(UI_FILT, 'pcg_solvent_check', 'binary', 'PCG solvent prior strength check', &
        &'Validation of the automatic solvent-prior strength: re-solve the whole strength grid for real and print the re-solve objective and residuals beside the closed-form estimate; eight extra pair solves per state and iteration, no effect on the result(yes|no){no}', '', .false., 'no', group="filter", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call refine3D%add_input(UI_FILT, 'pcg_mskfile', 'file', 'PCG support-constraint mask volume', &
        &'Real-space [0,1] mask volume installed as the hard support constraint of every PCG solve (the projected '//&
        &'system P H P; experimental focused/support mode); spherical mskdiam support when absent', &
        &'e.g. focusmask.mrc', .false., '', group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
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
        call add_refine3D_workflow_inputs(refine3D_auto, .false.)
        call add_ui_program('refine3D_auto', refine3D_auto, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_auto

    subroutine new_refine3D_local( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call refine3D_local%new(&
        &'refine3D_local', &
        &'Continue single-state 3D refinement from stored particle poses', &
        &'uses bounded continuous Cartesian pose refinement with existing poses and even/odd references', &
        &'simple_exec', .true., visibility=UI_VIS_STANDARD, display_name='Local 3D Refinement')
        call add_refine3D_workflow_inputs(refine3D_local, .true.)
        call add_ui_program('refine3D_local', refine3D_local, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_local

    subroutine add_refine3D_workflow_inputs( prg, l_local )
        type(ui_program), intent(inout) :: prg
        logical, intent(in) :: l_local
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call prg%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', 'Starting reference volume &
        & for particle matching', 'input starting volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call prg%add_input(UI_PARM, 'rec_backend', 'multi', 'Reconstruction backend', &
        &'Reconstruction backend for per-iteration half-map assembly, forwarded to the refine3D child and the '//&
        &'bootstrap/final reconstructions(gridding|pcg){gridding}', &
        &'', .false., 'gridding', group="search", &
        &choices=ui_choices([character(len=8) :: 'gridding', 'pcg']), visibility=UI_VIS_ADVANCED)
        call prg%add_input(UI_PARM, 'state', 'num', 'State to refine', &
        &'Refine this state alone, on a copy of the project that holds it as its only state; the project '//&
        &'itself is not written. With m_estimator=flex the copy holds every particle with a positive weight '//&
        &'for the state, otherwise the particles labelled with it{all particles}', 'state index{none}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call prg%add_input(UI_PARM, 'm_estimator', 'multi', 'Particle weighting', &
        &'flex weighs every particle of the refined state by its frozen weight in the project''s state weight set '//&
        &'(written by flex_pca) in every reconstruction; requires state. no reconstructs from the hard state '//&
        &'labels(no|flex){no}', '', .false., 'no', &
        &choices=ui_choices([character(len=4) :: 'no', 'flex']), visibility=UI_VIS_ADVANCED)
        ! search controls
        call prg%add_input(UI_SRCH, maxits,      required_override=.false., group="search", &
        &visibility=UI_VIS_ADVANCED)
        call prg%add_input(UI_SRCH, 'minits', 'num', 'Minimum automatic iterations', &
        &'Minimum number of iterations used when automatic refinement calculates its iteration budget', &
        &'minimum iterations', .false., 10., group="search", visibility=UI_VIS_DEVELOPER)
        call prg%add_input(UI_SRCH, 'nsample', 'num', 'Particle samples', &
        &'Number of particles sampled per refinement iteration', 'number of samples', .false., 25000., &
        &group="search", visibility=UI_VIS_DEVELOPER)
        if( .not. l_local )then
            call prg%add_input(UI_SRCH, 'ref_pose_init', 'multi', 'External-reference pose initialization', &
            &'When vol1 has independent provenance, run one fixed-reference CC pose-initialization pass at 15 Angstroms &
            &before Euclidean refinement(cc|none){none}', '', .false., 'none', group='search', &
            &choices=ui_choices([character(len=4) :: 'cc', 'none']), visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'pose_cont', 'multi', 'Continuous pose refinement', &
            &'Run a continuous Cartesian polishing pass after each discrete iteration(no|yes){no}', '', .false., 'no', &
            &group='search', choices=ui_choices([character(len=3) :: 'no', 'yes']), &
            &visibility=UI_VIS_ADVANCED)
        endif
        call prg%add_input(UI_SRCH, 'athres_cont', 'num', 'Continuous pose rotation bound', &
        &'Halfwidth of the total rotation a continuous Cartesian pose solve (local refinement or pose_cont=yes) may '// &
        &'move a particle from its seed: the capture range for wrongly assigned orientations', &
        &'in degrees{10}', .false., 10., group='search', visibility=UI_VIS_ADVANCED)
        if( l_local )then
            call prg%add_input(UI_SRCH, trs, group='search', visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'cont_conv_rot', 'num', 'Continuous pose convergence rotation', &
            &'A continuous pose (refine3D_local) is stable when it moved by at most this rotation since its last pose', &
            &'in degrees{0.5}', .false., 0.5, group='search', visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'cont_conv_shift', 'num', 'Continuous pose convergence shift', &
            &'A continuous pose (refine3D_local) is stable when it moved by at most this shift since its last pose', &
            &'in Angstroms{1.0}', .false., 1.0, group='search', visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'cont_conv_frac', 'num', 'Continuous pose convergence fraction', &
            &'Continuous refinement converges when this fraction of every state is stable, over iterations that '// &
            &'together sampled this fraction of the particles', '(0,1]{0.9}', .false., 0.9, group='search', &
            &visibility=UI_VIS_ADVANCED)
        endif
        call prg%add_input(UI_SRCH, pgrp,                                  group="search", visibility=UI_VIS_STANDARD)
        if( .not. l_local )then
            call prg%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
            &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', .false., 'no', group="search", &
            &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
            &visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'regpass', 'binary', 'Global registration pass', &
            &'One global (refine=prob) registration pass of all particles against the masked startup references, &
            &band-limited at the FSC=regpass_fsc resolution of the startup pair, before the neighbourhood iterations(yes|no){yes}', &
            &'', .false., 'yes', group="search", &
            &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
            &visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_SRCH, 'regpass_fsc', 'num', 'Registration-pass FSC criterion', &
            &'FSC value of the startup pair whose resolution band-limits the registration pass', &
            &'FSC value in (0,1){0.143}', .false., 0.143, group="search", &
            &visibility=UI_VIS_ADVANCED)
        endif
        if( l_local )then
            call prg%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
            &for accelerated computation(yes|no){no}','', .false., 'no', group="search", &
            &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
            &visibility=UI_VIS_ADVANCED)
        else
            call prg%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
            &for accelerated computation(yes|no){yes}','', .false., 'yes', group="search", &
            &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
            &visibility=UI_VIS_ADVANCED)
        endif
        call prg%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}','', .false.,&
        &'no', group="search", &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call prg%add_input(UI_FILT, 'amsklp', 'num', 'NU evidence envelope smoothing limit',&
        & 'Low-pass limit for NU evidence envelope generation in Angstroms', &
        &'low-pass limit in Angstroms', .false., 8., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        if( l_local )then
            call prg%add_input(UI_FILT, 'lpstop', 'num', 'Matching resolution cap', &
            &'Cap the matching band while keeping separate even/odd references', 'resolution in Angstroms', &
            &.false., 1., group='filter', visibility=UI_VIS_ADVANCED)
            call prg%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
            &'Half-specific filtering mode(none|nonuniform){nonuniform}', '', .false., 'nonuniform', group='filter', &
            &choices=ui_choices([character(len=10) :: 'none', 'nonuniform']), visibility=UI_VIS_ADVANCED)
        else
            call prg%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
            &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}','', .false., 'nonuniform', group="filter", &
            &choices=ui_choices([character(len=16) :: 'none', 'nonuniform', 'nonuniform_lpset']), &
            &visibility=UI_VIS_ADVANCED)
        endif
        call prg%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Maximum kernel PCG iterations during refinement; the cold original-sampling final reconstruction uses at least 5', &
        &'iterations{2}', &
        &.false., 2., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call prg%add_input(UI_FILT, 'maxits_ml', 'num', 'Regularized-solve PCG iterations', &
        &'Coupled PCG iterations of the ML-regularized system from the closed-form Wiener start; 0 = closed form only', 'iterations{0}', &
        &.false., 0., group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call prg%add_input(UI_FILT, 'pcg_solvent', 'binary', 'PCG soft solvent prior', &
        &'Soft solvent prior on the PCG base solve: a real-space ridge pulling solvent toward zero, solvent identified '//&
        &'per half from a prior-free solve of that half (smoothed absolute density, Otsu, logistic weight), then the '//&
        &'same cold solve again with the ridge; half-independent, so the pair stays gold standard; the support is untouched(yes|no){no}', '', .false., 'no', group="filter", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('rec_backend', [character(len=3) :: 'pcg']))
        call prg%add_input(UI_FILT, 'pcg_solvent_lambda', 'num', 'PCG solvent prior strength', &
        &'Ridge coefficient of the solvent prior relative to the data scale (1 = as strong as the low-band data term); '//&
        &'not given = estimated per state and iteration by cross-validation of the prior-free half pair', 'coefficient{auto}', &
        &.false., 1.0, group="filter", visibility=UI_VIS_ADVANCED, &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call prg%add_input(UI_FILT, 'pcg_solvent_check', 'binary', 'PCG solvent prior strength check', &
        &'Validation of the automatic solvent-prior strength: re-solve the whole strength grid for real and print the re-solve objective and residuals beside the closed-form estimate; eight extra pair solves per state and iteration, no effect on the result(yes|no){no}', '', .false., 'no', group="filter", &
        &visibility=UI_VIS_ADVANCED, choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &activation=ui_activation_equals_any('pcg_solvent', [character(len=3) :: 'yes']))
        call prg%add_input(UI_FILT, envfsc, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        call prg%add_input(UI_FILT, envmsklp, group="filter", &
        &visibility=UI_VIS_ADVANCED)
        if( .not. l_local )then
            call prg%add_input(UI_FILT, combine_eo, group="filter", &
            &visibility=UI_VIS_ADVANCED)
        endif
        call prg%add_input(UI_FILT, 'res_target', 'num', 'Resolution target (in A)',&
        & 'Resolution target in Angstroms', 'Resolution target in Angstroms', .false., 3., group="filter", &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call prg%add_input(UI_MASK, mskdiam, group="mask", visibility=UI_VIS_STANDARD)
        call prg%add_input(UI_MASK, automsk_refine3D, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        call prg%add_input(UI_MASK, nu_msk_sig, group="mask", &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call prg%add_input(UI_COMP, nparts, group="compute", visibility=UI_VIS_STANDARD)
        call prg%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
    end subroutine add_refine3D_workflow_inputs




end module simple_ui_refine3D
