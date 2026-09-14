!@descr: module defining the user interfaces for heterogeneity-analysis programs in the simple_exec suite
module simple_ui_heterogeneity
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = &
    category_descriptor('heterogeneity', 'Heterogeneity Analysis', 65)
type(ui_program), target :: flex_pca
type(ui_program), target :: refine3D_states
type(ui_program), target :: classify3D_refs
type(ui_program), target :: ptcl3D_state_consensus

contains

    subroutine construct_heterogeneity_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_flex_pca(prgtab)
        call new_refine3D_states(prgtab)
        call new_classify3D_refs(prgtab)
        call new_ptcl3D_state_consensus(prgtab)
    end subroutine construct_heterogeneity_programs

    subroutine new_flex_pca( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call flex_pca%new(&
        &'flex_pca',&
        &'Projection-aware covariance heterogeneity states',&
        &'fits a sigma-whitened low-rank 3D covariance factor directly from fixed-pose particles, performs MAP latent inference, and reconstructs combined/even/odd kernel states without a diffusion graph',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_STANDARD, display_name='Covariance Heterogeneity States')
        call flex_pca%add_input(UI_IMG, 'vol1', 'file', &
            'Consensus mean volume', 'Fixed mean subtracted from every particle; when omitted, the &
            &project consensus map (out segment, state 1) at native sampling is used', &
            'e.g. vol1.mrc (consensus mean)', .false., '', &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'neigs', 'num', &
            'Covariance components (default 16)', 'Number of fitted low-rank covariance factors; capped at 48', &
            '# components', .false., 16.0, &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'npreimages', 'num', &
            'Max state volumes (default 16)', &
            'Upper bound on the kernel-regression targets in latent space; with the default state_axis=0 &
            &these are diffusion k-centers over all retained components. The two-gate merge collapses &
            &indistinct states, so the recovered count is <= this', &
            'max # states 3-32', .false., 16.0, &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'min_state_frac', 'num', &
            'Minimum state population fraction (default 0 = off)', &
            'Population floor: every delivered state must hold at least this fraction of the embedded &
            &particles. Under-populated clusters are dropped, the targets re-placed on the retained &
            &particles with the count raised by the deficit until npreimages states qualify, and the &
            &dropped or unassigned particles receive a random label among the delivered states, which &
            &are then reconstructed from their hard labels. Incompatible with preimage_auto and the merge', &
            'fraction of particles 0-1', .false., 0.0, &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'preimage_auto', 'binary', &
            'Determine the state count automatically (default no)', &
            'Raises the state ceiling to 32 (unless npreimages is given) and enables the two-gate merge, &
            &so the delivered state count is recovered from the data rather than requested(yes|no){no}', &
            '', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_FILT, 'niter', 'num', &
            'Covariance fit iterations (default 5)', 'Alternating projection/backprojection covariance-factor iterations', &
            '# iterations 1-20', .false., 5.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'state_axis', 'num', &
            'Latent axis for state targets (default 0 = diffusion k-center)', &
            'With 0 the state targets are diffusion k-centers over ALL retained covariance components, &
            &which covers a continuous reaction coordinate and branched compositional states with the &
            &same constants. A negative value places them along a density-spread path instead. &
            &A positive value places them along that single component, &
            &which discards the other components and tends to concentrate the particles on one state. &
            &SIMPLE_COV_KMEANS=1 recovers the former k-means placement', &
            'component index', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'nkern', 'num', &
            'Latent components entering the state kernel (default 0 = all)', &
            'Decouples the state stage from neigs. neigs sets how many eigenvolumes are estimated; &
            &nkern sets how many of them define "nearby" for target placement and kernel weighting. &
            &Components past the first few are usually dominated by fitting noise, and each one still &
            &contributes to the Mahalanobis distance, so leaving them in lets noise directions decide &
            &which particles pool at a target and drives every state map toward the consensus. Rank &
            &components by observed spread over posterior variance and keep those above ~1.5', &
            'leading components, 0=all', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'infile', 'file', &
            'Cached embedding to resume from', &
            'Path to a flex_pca_embedding.bin written by an earlier run. Skips the covariance &
            &basis fit and the per-particle embedding (~77% of the runtime) and re-enters at the &
            &state-weighting stage, so a different npreimages/min_neff/state_axis can be tried &
            &without refitting. The cache is tied to the particle selection it was built from', &
            'e.g. flex_pca_embedding.bin', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'nbins', 'num', &
            'Cross-validated bandwidth bins (default 1 = off)', &
            'Number of kernel bandwidth bins swept per state. With 1 the &
            &chi-squared floor is used directly, which is the most local and noisiest choice. &
            &With >1 each bin is reconstructed on both halfsets and the bandwidth maximising &
            &even/odd agreement is kept per state. Costs 2*nbins extra reconstruction passes', &
            '# bins, 1=off', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! min_neff is not a flex_pca input: on the default path the GMM replaces the kernel weights and
        ! bandwidth, so it cannot change the maps. Reachable as SIMPLE_COV_MIN_NEFF for the opt-out paths.
        call flex_pca%add_input(UI_FILT, 'heldout', 'binary', &
            'Cross-halfset (held-out) embedding', &
            'Fit the covariance basis on one halfset and embed the other, then swap, so no particle is &
            &projected onto a basis estimated from it; removes in-sample bias and reports the halfset &
            &subspace principal angles. Costs two covariance estimations', &
            '(yes|no){no}', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, 'column_separation', 'num', &
            'Minimum grid separation between columns (default 2)', &
            'Selected frequencies closer than this are suppressed; also decorrelates the column noise', &
            'grid units', .false., 2.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_SRCH, 'n_probe_iters', 'num', &
            'Probe subspace-iteration refinements (default 5)', &
            'EM / probe subspace iterations refining the column basis. Probe volumes aggregate the whole &
            &slice instead of one Fourier voxel, which is the main lever on per-particle latent quality. &
            &An upper bound rather than a fixed count: the loop stops early once the mean principal-angle &
            &cosine between successive bases reaches 0.97. Set 0 to disable', &
            '# iterations', .false., 5.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_FILT, lp, required_override=.false., &
            label_override='Low-pass limit (derived: 2.5*smpd_crop)', &
            group="regularization", visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_PARM, 'box_crop', 'num', &
            'Working box size (default 64)', 'Even low-resolution box used for covariance fitting and the latent embedding', &
            'pixels', .false., 64.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'box_rec', 'num', &
            'State-map reconstruction box (default: native project box)', &
            'Even box for the delivered state maps; decoupled from box_crop so the maps are not &
            &limited to the covariance Nyquist. The commander resolves it to the native project box, &
            &so the maps come out at the native sampling; it falls back to box_crop only when the &
            &project geometry cannot be read. Capped at the native box', &
            'pixels', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_PARM, 'oritype', 'str', &
            'Particle orientation segment', 'Fixed to ptcl3D', 'ptcl3D', .false., 'ptcl3D', &
        &visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_SRCH, sigma_est, visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_MASK, mskdiam, required_override=.false., &
            group="mask", visibility=UI_VIS_STANDARD)
        call flex_pca%add_input(UI_COMP, nparts, required_override=.false., &
            group="compute", visibility=UI_VIS_ADVANCED)
        call flex_pca%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call add_ui_program('flex_pca', flex_pca, prgtab, UI_CATEGORY)
    end subroutine new_flex_pca

    subroutine new_refine3D_states( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call refine3D_states%new(&
        &'refine3D_states',&
        &'Refine conformational states from a same-lineage particle orientation scaffold',&
        &'refines conformational states from existing particle poses. The pose policy fixes the projection direction, searches a local neighborhood, or permits global matching',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='Conformational State Refinement')
        call refine3D_states%add_input(UI_SRCH, maxits, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, nstates, required_override=.false., group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'flex', 'binary', 'Initialize states with flex PCA', &
        &'Run flex_pca to derive the initial particle states and state volumes; the default for state=0/1 input, &
        &skipped automatically when the project already carries multi-state labels; flex=no selects stochastic &
        &state initialization(yes|no){yes}', '', &
        &.false., 'yes', group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'min_state_frac', 'num', 'Minimum flex state population fraction', &
        &'With flex=yes, every initial state must hold at least this fraction of the particles; under-populated &
        &flex clusters are dropped and their particles randomized over the delivered states{0.1}', &
        &'fraction of particles', .false., 0.1, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'nsample', 'num', 'Particle sample target', &
        &'Particles sampled per iteration; set 0 to derive the automatic target from the number of states', &
        &'particles (0=automatic)', .false., 0., group='search', visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call refine3D_states%add_input(UI_SRCH, 'sticky_class_sampling', 'binary', 'Reuse one sampled cohort', &
        &'Keep a projection-balanced stochastic cohort fixed across frequency stages(yes|no){no}', '', &
        &.false., 'no', group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, sigma_est, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'pose_policy', 'multi', 'Pose-search policy', &
        &'Pose-search policy(fixed|local|global){global}. Fixed keeps the projection direction and optimizes the in-plane angle and x/y translations; local uses the current geometric neighborhood; global permits full probabilistic matching', &
        &'fixed, local, or global', .false., 'global', group='search', &
        &choices=ui_choices([character(len=6) :: 'fixed', 'local', 'global']), visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_SRCH, 'local_ang_bound', 'num', 'Local projection-angle bound', &
        &'Advanced override for the automatically derived local projection-direction neighborhood in degrees', &
        &'degrees (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, 'local_inpl_bound', 'num', 'Local in-plane angle bound', &
        &'Advanced upper bound for the automatically derived probabilistic in-plane neighborhood in degrees', &
        &'degrees (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, 'local_shift_bound', 'num', 'Local shift bound', &
        &'Advanced override for the automatically derived local translational half-width in pixels', &
        &'pixels (-1=automatic)', .false., -1., group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_SRCH, 'inpl_cont', 'binary', 'Continuous in-plane refinement', &
        &'Joint continuous Euclidean in-plane and shift refinement(yes|no){yes}', '', .false., 'yes', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_SRCH, ptcl_src, group='search', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) by their center of gravity and map shifts back to the particles(yes|no){no}', '', &
        &.false., 'no', group='search', choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', &
        &'Continue previous refinement(yes|no){no}', '', .false., 'no', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_SRCH, 'overlap', 'num', 'State-overlap convergence target', &
        &'Required overlap of state assignments for convergence during frequency marching', 'overlap fraction', &
        &.false., .99, group='search', visibility=UI_VIS_DEVELOPER)
        call refine3D_states%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(fsc|nonuniform_lpset|none){nonuniform_lpset}', '', .false., 'nonuniform_lpset', group='filter', &
        &choices=ui_choices([character(len=16) :: 'fsc', 'nonuniform_lpset', 'none']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Starting low-pass limit for the common frequency schedule across states', 'low-pass limit in Angstroms', &
        &.false., 10., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Final low-pass limit for the common frequency schedule across states', 'low-pass limit in Angstroms', &
        &.false., 6., group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, ml_reg, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, envfsc, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_FILT, envmsklp, group='filter', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_MASK, 'automsk', 'binary', 'Perform envelope masking', &
        &'Generate/apply the NU-evidence envelope; requires filt_mode=nonuniform_lpset(yes|no){no}', '', .false., 'no', group='mask', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_MASK, nu_msk_sig, group='mask', visibility=UI_VIS_ADVANCED)
        call refine3D_states%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call refine3D_states%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('refine3D_states', refine3D_states, prgtab, UI_CATEGORY)
    end subroutine new_refine3D_states

    subroutine new_classify3D_refs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call classify3D_refs%new(&
        &'classify3D_refs',&
        &'Classify particles against supplied references and reconstruct the resulting state maps',&
        &'competitively classifies particles against a complete set of references that may have independent provenance, then reconstructs maps from the resulting hard assignments',&
        &'simple_exec',&
        &.true.,&
        &visibility=UI_VIS_STANDARD, display_name='Reference-guided 3D Classification')
        call classify3D_refs%add_input(UI_IMG, 'vol1', 'file', 'Reference volume', &
        &'First member of the complete vol1..volN fixed-reference set', 'input volume e.g. vol1.mrc', .false., 'vol1.mrc', &
        &visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, 'maxits', 'num', 'Maximum iterations', &
        &'Total number of frequency-marched classification iterations; must be >= 1', 'iterations{50}', &
        &.false., 50., group='search', visibility=UI_VIS_ADVANCED, preserve_default=.true.)
        call classify3D_refs%add_input(UI_SRCH, nstates, required_override=.true., group='search', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, 'nsample', 'num', 'Particle sample target per state', &
        &'Particles sampled per iteration per state', 'particles{10000}', .false., 10000., group='search', &
        &visibility=UI_VIS_DEVELOPER, preserve_default=.true.)
        call classify3D_refs%add_input(UI_SRCH, 'inpl_cont', 'binary', 'Continuous in-plane refinement', &
        &'Joint continuous Euclidean in-plane and shift refinement(yes|no){yes}', '', .false., 'yes', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_SRCH, pgrp, group='search', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_SRCH, ptcl_src, group='search', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', &
        &'Center reference volume(s) and map shifts back to particles(yes|no){no}', '', .false., 'no', group='search', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
        &'Filtering mode(nonuniform_lpset|none){nonuniform_lpset}', '', .false., 'nonuniform_lpset', group='filter', &
        &choices=ui_choices([character(len=16) :: 'nonuniform_lpset', 'none']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'lpstart', 'num', 'Starting low-pass limit', &
        &'Starting low-pass limit for frequency-marched classification', 'low-pass limit in Angstroms', .false., 10., &
        &group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', &
        &'Final low-pass limit for frequency-marched classification', 'low-pass limit in Angstroms', .false., 6., &
        &group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, ml_reg, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, envfsc, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_FILT, envmsklp, group='filter', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_MASK, mskdiam, group='mask', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_MASK, 'automsk', 'binary', 'Perform envelope masking', &
        &'Generate/apply the NU-evidence envelope; requires filt_mode=nonuniform_lpset(yes|no){no}', '', .false., 'no', group='mask', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_MASK, nu_msk_sig, group='mask', visibility=UI_VIS_ADVANCED)
        call classify3D_refs%add_input(UI_COMP, nparts, group='compute', visibility=UI_VIS_STANDARD)
        call classify3D_refs%add_input(UI_COMP, nthr, group='compute', visibility=UI_VIS_STANDARD)
        call add_ui_program('classify3D_refs', classify3D_refs, prgtab, UI_CATEGORY)
    end subroutine new_classify3D_refs

    subroutine new_ptcl3D_state_consensus( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ptcl3D_state_consensus%new(&
        &'ptcl3D_state_consensus', &                                    ! name
        &'Build consensus particle-state assignments across projects',& ! summary
        &'is a program that builds a consensus particle state assignment from a file table of SIMPLE projects &
        &and writes it to the target project ptcl3D field', &           ! help
        &'simple_exec',&                                                ! executable
        &.false.)                                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ptcl3D_state_consensus%add_input(UI_FILE, projtab,&
        &help_override        = 'Text file listing SIMPLE project files (*.simple) containing ptcl3D state assignments',&
        &placeholder_override = 'e.g. projtab.txt',&
        &required_override          = .true.,&
        &group="data", visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call ptcl3D_state_consensus%add_input(UI_FILE, projfile,&
        &help_override        = 'Target SIMPLE project file that receives the consensus ptcl3D state assignment',&
        &required_override          = .true.,&
        &group="data", visibility=UI_VIS_STANDARD)
        call ptcl3D_state_consensus%add_input(UI_SRCH, nstates,&
        &help_override        = 'Number of state labels to match; inferred from projtab when omitted',&
        &required_override          = .false.,&
        &group="state", visibility=UI_VIS_ADVANCED)
        call ptcl3D_state_consensus%add_input(UI_PARM, prune,&
        &group="data", visibility=UI_VIS_ADVANCED)
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
        call add_ui_program('ptcl3D_state_consensus', ptcl3D_state_consensus, prgtab, UI_CATEGORY)
    end subroutine new_ptcl3D_state_consensus

end module simple_ui_heterogeneity

