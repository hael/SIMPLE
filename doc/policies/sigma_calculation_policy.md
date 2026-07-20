# Sigma Calculation Policy

## Scope

This document defines the current policy for Euclidean sigma calculation and use in SIMPLE.

It applies to workflows that run with `objfun='euclid'`, especially:

- `cluster2D`
- `abinitio2D`
- `refine3D`
- reconstruction paths that consume Euclidean sigma spectra

## Public policy

Sigma calculation is a Euclidean-objective feature.

- If `objfun='euclid'`, SIMPLE computes and persists sigma spectra.
- If `objfun!='euclid'`, the Euclidean sigma path is inactive.
- Sigma application during restoration/reconstruction is gated separately by `ml_reg='yes'`.

Two gates must never be conflated:

- **Sigma calculation/update gate:** `cc_objfun == OBJFUN_EUCLID`
- **ML restoration/reconstruction consumption gate:** `params%l_ml_reg`

Probabilistic refinement modes are not an exception. `prob_align`,
`prob_align2D`, `prob_tab`, `prob_tab2D`, and `prob_tab_neigh` consume the
current grouped sigma spectra for Euclidean scoring, but they do not produce the
final residual sigma update for the iteration. The final sigma update is always
the matcher pass after the probabilistic assignment has been applied to the
particle orientation/state records. Never guard matcher-side sigma updates with
`.not. l_prob_align_mode`, and never treat a probabilistic pre-alignment pass as
a replacement for post-assignment residual sigma calculation.

The effective reconstruction gate is `params%l_ml_reg`, which is only true when:

- `ml_reg='yes'`
- `cc_objfun == OBJFUN_EUCLID`

This is enforced in `src/main/params/simple_parameters_phases.f90` during derived-parameter validation.

## Ownership and persistence

The runtime owner of sigma state is `builder%esig`, the `euclid_sigma2` object held by the builder. That object owns:

- `sigma2_noise`: group-expanded per-particle spectra used by alignment/reconstruction
- `sigma2_part`: per-particle spectra written during the current search/update pass
- `sigma2_groups`: group spectra loaded from the STAR representation

See [src/main/simple_builder.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_builder.f90:30) and [src/main/simple_euclid_sigma2.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_euclid_sigma2.f90:15).

Persisted sigma data has two forms:

- partition-local binary files: `sigma2_noise_partN.dat`
- grouped STAR files: `sigma2_it_<iter>.star`

The binary format is defined in [src/fileio/simple_sigma2_binfile.f90](/Users/elmlundho/src/SIMPLE/src/fileio/simple_sigma2_binfile.f90:1).

Probabilistic commanders own candidate-table generation and assignment
artifacts only. They may load sigma through the normal matcher preparation path,
but ownership of durable sigma update remains with:

- matcher/search code: per-particle residual spectra after assignment
- Euclid commanders/strategies: bootstrap and group consolidation
- reconstruction/restoration code: ML consumption of grouped spectra

## Grouping policy

`sigma_est` controls grouping:

- `group`: group by particle `stkind`
- `global`: use one shared group per even/odd half-set

This is normalized into `params%l_sigma_glob` in `src/main/params/simple_parameters_phases.f90`.

When group STAR files are loaded, the selected group spectrum is expanded back onto each active particle according to:

- particle `eo`
- particle `stkind` unless `l_sigma_glob`

See [src/main/simple_euclid_sigma2.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_euclid_sigma2.f90:171).

## Native-grid and staged-refinement policy

Persistent sigma spectra live on the original particle-image grid, not on a
stage's cropped working grid:

```text
native shell k  <->  params%box, params%smpd
storage range   =    1 : fdim(params%box)-1
```

`calc_pspec`, the matcher residual update, and `euclid_sigma2` all use this
native shell convention. `box_crop` and `smpd_crop` only determine how much of
that stored model a stage consumes. A crop change between refine3D stages does
not change the meaning or size of a stored sigma record and is not a reason to
recalculate or rewrite it.

The native records already contain shells outside an early stage's active
matching band. When a later stage activates those shells, the normal matcher
residual update replaces their values. No stage-handoff extrapolation or
interpolation is required.

Reconstruction already owns the only required grid adaptation.
`image%gen_fplane4rec` calls `simple_math_ft::upsample_sigma2` to interpolate
the selected native sigma range onto the padded reconstruction Fourier plane.
This is an existing parallel consumer path and must not be duplicated in stage
initialization. See
[src/main/image/simple_image_ctf.f90](/Users/elmlundho/src/SIMPLE/src/main/image/simple_image_ctf.f90:203)
and
[src/utils/math/simple_math_ft.f90](/Users/elmlundho/src/SIMPLE/src/utils/math/simple_math_ft.f90:1).

A change to the original `params%box` or `params%smpd` is a different native
particle representation and is outside this reuse contract.

## Calculation lifecycle

### 1. Bootstrap spectra

`calc_pspec` is the initial image-power-spectrum bootstrap. For refine3D it runs
on the first sigma-producing stage (`startit <= 1`), after enforcing even/odd
partitioning. A later staged run (`startit > 1`) reuses the existing
`sigma2_noise_partN.dat` files instead. The normal first
`calc_group_sigmas` call then reads those files and creates the grouped model
for the new iteration.

This reuse needs no public command-line mode, controller handoff metadata,
binary-format extension, file conversion, or new interpolation path. The staged
controller already advances `startit` and retains the native-grid partition
files in the working directory. Existing readers remain responsible for file
validation. A missing or unusable file fails through that established path; it
must not silently trigger another expensive `calc_pspec` bootstrap. See
[src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:523).

`abinitio2D` sets `sigma_est='global'` and, unless overridden, uses the default
Euclidean objective with `ml_reg='yes'`. It runs `calc_pspec` before the first
`cluster2D` stage, so the first 2D search stage has grouped sigma spectra
available. See
[src/main/commanders/simple/simple_commanders_abinitio2D.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_abinitio2D.f90:38) and
[src/main/commanders/simple/simple_commanders_abinitio2D.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_abinitio2D.f90:276).

`abinitio3D` sets `objfun='euclid'` and `sigma_est='global'`, then delegates
stage execution to `refine3D`. Its first stage bootstraps sigma and later stages
reuse the residual partition files written by the preceding matcher. See
[src/main/commanders/simple/simple_commanders_abinitio.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_abinitio.f90:444) and
[src/main/simple_abinitio_utils.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_abinitio_utils.f90:283).

`calc_pspec` bootstraps two global noise spectra, one for each even/odd half-set,
from masked particle images. It selects a random, even/odd-stratified sample of
at most 25,000 active particles across the complete run, so distributed workers
only evaluate their intersection with that shared selection. The assembler uses
the existing global-estimate scaling and propagation path. In the bootstrap
assembly step:

- sampled spectra are scaled and averaged globally per even/odd half-set
- the average image power spectrum is subtracted
- negative values are repaired by neighbor propagation
- the requested grouped STAR layout is written with identical initial global
  curves in every group, so existing `sigma_est='group'` readers remain valid
- the two global curves are propagated to every active record in partition-local
  `sigma2_noise_partN.dat` files

This is the established `sigma_est='global'` propagation behavior, now used for
all initial bootstraps. It is not a permanent replacement for group estimation:
the normal orientation-conditioned matcher residual update and later group
consolidation specialize the initially shared curves.

See [src/main/strategies/parallelization/simple_calc_pspec_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_calc_pspec_strategy.f90:130) and [src/main/commanders/simple/simple_commanders_euclid_distr.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_euclid_distr.f90:15).

### 2. Per-particle update during search

During Euclidean 2D/3D matching, sigma is recalculated after the particle orientation/shift update.

- 2D uses `refkind='class'`
- 3D uses `refkind='proj'`

The residual is formed in polar Fourier space after shift, rotation, and CTF
application, and the per-shell squared residual energy is normalized by
`2 * pftsz`. See
[src/main/pftc/simple_polarft_corr.f90](/Users/elmlundho/src/SIMPLE/src/main/pftc/simple_polarft_corr.f90:801).

Call sites:

- [src/main/strategies/search/simple_strategy2D_matcher.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/search/simple_strategy2D_matcher.f90:153)
- [src/main/strategies/search/simple_strategy3D_matcher.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/search/simple_strategy3D_matcher.f90:148)

At the end of the search pass, the updated per-particle spectra are written back to the partition binary file by `build%esig%write_sigma2`.

This update must run in probabilistic and non-probabilistic modes. In
probabilistic modes, the sequence is:

1. consolidate/load the current grouped sigma spectra for the iteration
2. run `prob_align*` to sample particles, generate probability tables, and write
   the assignment artifact
3. run the normal matcher, which consumes the assignment, updates
   orientation/state/shift records, recalculates per-particle residual sigmas,
   and writes `sigma2_noise_partN.dat`
4. consolidate partition sigmas into the next grouped STAR handoff

The important bug-prevention rule is that step 3 is mandatory for `prob`,
`prob_state`, and `prob_neigh`. The probabilistic table is an accelerated
assignment mechanism, not the sigma update mechanism.

The assignment-only side of this split is visible in
[src/main/commanders/simple/simple_commanders_prob.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_prob.f90:172) for 3D and
[src/main/commanders/simple/simple_commanders_prob.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_prob.f90:363) for 2D.

### 3. Iteration consolidation

The policy is that grouped STAR files are the durable handoff artifact between iterations and runs.

`sigma2_it_<N>.star` is the grouped model consumed by matcher iteration `N`.
`sigma2_group_iter` encodes this convention for both workflows: consolidation
before matching writes the current iteration number, while consolidation after
matching writes the next iteration number.

Current implementation details differ slightly by workflow:

- `cluster2D` consolidates after the iteration and writes `which_iter + 1` for
  the next iteration. See
  [src/main/strategies/parallelization/simple_cluster2D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_cluster2D_strategy.f90:285).
- shared-memory `refine3D` consolidates and writes `which_iter` at iteration entry.
  See [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:442) ..
- distributed `refine3D` follows the same policy: current-iteration
  consolidation before `prob_align*`/worker scheduling. See
  [src/main/strategies/parallelization/simple_refine3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_refine3D_strategy.f90:834).

The authoritative consolidation implementation is `exec_calc_group_sigmas` in [src/main/commanders/simple/simple_commanders_euclid.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_euclid.f90:41).

That implementation:

- reads all partition `sigma2_noise_partN.dat` files
- groups by `eo` and either `stkind` or one global bucket
- computes arithmetic group averages over active particles; each active
  particle contributes one spectrum to its even/odd and `stkind` (or global)
  bucket
- writes `sigma2_it_<which_iter>.star`

This is intentionally different from bootstrap `calc_pspec_assemble`, which
derives spectra from average image-power subtraction and repairs negative
values before writing its grouped result.

## Workflow summaries

### `abinitio2D`

`abinitio2D` is a staged `cluster2D` driver. It sets `sigma_est='global'` and,
unless the user overrides the objective, uses ML-regularized Euclidean 2D
output. The first stage bootstraps sigma with `calc_pspec`; each `cluster2D`
iteration then updates particle assignments, recalculates per-particle residual
sigma in `cluster2D_exec`, and consolidates grouped sigma for the next
iteration. Final class-average generation is ML-regularized and records the
final grouped sigma STAR file in the project output. See
[src/main/commanders/simple/simple_commanders_abinitio2D.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_abinitio2D.f90:347) and
[src/main/commanders/simple/simple_commanders_abinitio2D.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_abinitio2D.f90:358).

When `abinitio2D` uses a probabilistic `refine='prob*'` mode, `prob_align2D` owns sampling and
probabilistic class assignment, but `cluster2D_exec` still owns the residual
sigma update after that assignment is consumed.

### `abinitio3D`

`abinitio3D` sets `objfun='euclid'` and `sigma_est='global'`, then executes
staged `refine3D` runs. The stage controller uses `snhc_smpl` in stage 1,
`shc_smpl` in stage 2, `prob` in middle stages, and `prob_neigh` in later
stages; fixed-orientation multivolume mode uses `prob_state`. See
[src/main/simple_abinitio_controller.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_abinitio_controller.f90:220).

ML consumption is staged separately from sigma calculation. Standard
`abinitio3D` disables `ml_reg` for stages 1-2 and enables it for stages 3-8,
while the `abinitio3D_cavgs` route forces `ml_reg='no'` for refinement stages.
Sigma calculation still follows `objfun='euclid'`; `ml_reg` only decides whether
the grouped spectra are applied during reconstruction/restoration. See
[src/main/simple_abinitio_controller.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_abinitio_controller.f90:216) and
[src/main/simple_abinitio_controller.f90](/Users/elmlundho/src/SIMPLE/src/main/simple_abinitio_controller.f90:262).

For every `abinitio3D` stage with a probabilistic refinement mode, the same
critical invariant applies: `prob_align*` chooses/writes assignments, and the
following `refine3D_exec` pass must recalculate residual sigma after those
assignments have been applied.

## Consumption policy

Grouped sigma spectra are loaded into `builder%esig` and expanded per particle.
They are consumed in two distinct places:

- Euclidean objective scoring during 2D/3D matching whenever
  `cc_objfun == OBJFUN_EUCLID`
- restoration/reconstruction only when `params%l_ml_reg` is true

For matching, sigma enters the polar Fourier objective functions as
`sigma2_noise`. See
[src/main/pftc/simple_polarft_corr.f90](/Users/elmlundho/src/SIMPLE/src/main/pftc/simple_polarft_corr.f90:58).

For 3D reconstruction, `params%l_ml_reg` causes the sigma groups to be loaded
before `calc_3Drec`. See
[src/main/commanders/simple/simple_commanders_rec.f90](/Users/elmlundho/src/SIMPLE/src/main/commanders/simple/simple_commanders_rec.f90:67) and
[src/main/strategies/parallelization/simple_rec3D_strategy.f90](/Users/elmlundho/src/SIMPLE/src/main/strategies/parallelization/simple_rec3D_strategy.f90:149).

For 2D class restoration, the same grouped sigma state is loaded before Fourier-plane generation. See [src/main/class/simple_classaverager_restore.f90](/Users/elmlundho/src/SIMPLE/src/main/class/simple_classaverager_restore.f90:79).

The Fourier-plane policy is:

- upsample sigma from the cropped shell grid to the padded reconstruction grid
- divide both the complex Fourier sample and its CTF power term by sigma

See [src/main/image/simple_image_ctf.f90](/Users/elmlundho/src/SIMPLE/src/main/image/simple_image_ctf.f90:203) and [src/main/image/simple_image_ctf.f90](/Users/elmlundho/src/SIMPLE/src/main/image/simple_image_ctf.f90:263).

## Current architectural boundary

The current boundary is acceptable and should be preserved:

- `simple_euclid_sigma2`: sigma object model and STAR/binary translation
- matcher/search code: per-particle sigma updates
- probabilistic commanders: candidate probability tables and assignment artifacts only
- Euclid commanders/strategies: bootstrap and iteration-level consolidation
- refine3D stage initialization: bootstrap once, then retain existing native-grid
  partition files
- reconstruction/restoration code: sigma consumption when `l_ml_reg`

The main rule to preserve is simple: sigma remains on the native particle grid.
Grouped STAR files are the normal cross-iteration model, while the newer
matcher-written partition files are retained across refine3D stage boundaries
so the next stage can consolidate them without rerunning `calc_pspec`.

## Forward-looking evolution

The most plausible future simplification is to make all Euclidean sigma
handling follow the current `sigma_est='global'` policy and remove the
`sigma_est` option. Bootstrap and steady-state consolidation would then produce
one spectrum per even/odd half-set rather than separate `stkind` spectra, and
matching and reconstruction would use those global curves for every active
particle.

This would not remove the matcher-side per-particle residual calculation. Those
residuals remain the observations from which the next global even/odd spectra
are estimated.

This direction is deliberately deferred. Before removing the grouped mode,
`sigma_est='global'` must be validated in high-resolution refinement against
the current grouped policy. The comparison must cover representative homogeneous
and mixed particle populations and demonstrate comparable:

- FSC and reported resolution;
- high-frequency stability and map quality;
- orientation, shift, and convergence behavior;
- robustness when stacks or particle groups have different noise levels.

Until that evidence exists, `sigma_est='group'` and `sigma_est='global'` remain
supported current behavior. Removal of the option is not part of the present
stage-reuse refactoring.
