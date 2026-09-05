# Staged Ab Initio 3D Reconstruction

## Problem

Infer particle orientations, translations, optional state labels, and 3D maps
without a trusted starting model. The [refine3D](refine3d.md) alternation
converges to whatever local optimum is nearest its starting point, and from
random orientations that optimum is usually a featureless blob or a
symmetric artifact. Ab initio 3D is a continuation scheme with the same three
control knobs as [ab initio 2D](abinitio2d.md): admitted bandwidth,
stochasticity of the search, and number of particles per iteration, plus two
that are specific to 3D: the angular sampling density and the point-group
symmetry imposed.

## Initialization

Particle orientations are drawn uniformly at random; with independent
multi-state initialization, state labels are also drawn uniformly. The starting
volume is deliberately uninformative: Gaussian noise of standard deviation
`5 * box/box_crop` plus a soft Gaussian sphere of radius `box/4`, written as
the map and both half maps. The first iteration therefore aligns particles to
a nearly isotropic reference and the reconstruction from those poses is the
first real signal. Alternatives are a reconstruction from the random
orientations (`rndstart`), class averages (`abinitio3D_cavgs`), or an external
volume after a correlation-based pose-initialization pass.

## Continuation schedule

The particle workflow has up to eight stages. Each stage is one `refine3D`
run with its own bandwidth, angular sampling, crop, iteration cap, search
mode, and sampling target:

| stage | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| max iterations | 20 | 20 | 17 | 15 | 12 | 12 | 12 | 25 |
| `nspace` | 500 | 1000 | 1000 | 1000 | 2500 | 2500 | 5000 | 5000 |
| search | prob_neigh (shc) | prob_neigh (shc) | prob | prob | prob | prob_neigh (state) | prob_neigh (state) | prob_neigh (state) |
| shifts | 0 | 0 | on | on | on | on | on | on |
| ML regularization | no | no | yes | yes | yes | yes | yes | yes |
| symmetry | start group | start group | start group, then search | target | target | target | target | target |

**Bandwidth.** The ladder is defined on FRC values, not resolutions. From
the class-average FRC curve `frc(k)` that produced the input, the two end
criteria are

```text
crit_1 = max(frc(k at lp_start), 0.65),    crit_n = max(frc(k at lp_final), 0.03),
```

stepped linearly across stages, and each stage's low-pass is the resolution at
which the FRC first drops below that stage's criterion, never finer than
`lp_final`. `lp_start` lies in `[10, 20]` A, `lp_final` is the median of the
three best class-average resolutions clamped to `[4.5, 6]` A. If the FRC is
unusable the ladder degrades to a linear one in resolution. Stage 1 is never
finer than 20 A. The stage-1 reference is band-limited by a Gaussian at its
low-pass (as in 2D) rather than by an FSC filter.

**Sampling grid.** Each stage crops to `smpd_target = max(2, lp/3)` at the
ends and interpolates the box in between, always preserving
`box * smpd = box_crop * smpd_crop`; the shift bound is `min(8, max(2, 12 A /
smpd_crop))` pixels. When a single downscaled particle cache is used, every
eligible stage adopts the final crop while the low-pass still marches.

**Search stochasticity.** Stages 1 and 2 use the stochastic direct
neighborhood (first-improvement over a random subset of directions);
stages 3 to 5 use the full probabilistic table with balanced assignment;
stages 6 to 8 restrict the table to pooled coarse-peak neighborhoods so the
late search is local and exhaustive within it. The coarse subspace is
rescaled with `nspace` so its granularity stays constant. Independent
multi-state runs use plain first-improvement hill climbing for the first two
stages.

**Particle subset.** A fixed target `nsample` (default 10 000 per state)
gives `update_frac = nsample * nstates / N_active`, capped at 0.9 so fractional
updates stay on; if the target exceeds 90 percent of the active set the run
switches to full updates and disables trailing. From stage 5 (4 for independent
multi-state) the subset is redrawn stochastically each iteration and only
the best-scoring half of each direction bin is eligible (`frac_best = 0.5`,
rising to 0.85 in the last two stages), which is a soft outlier rejection.

**Regularization.** ML regularization starts at stage 3, once the halves are
independent enough for the FSC-derived signal prior to mean something. From
stage 6 the reference bandwidth may become spatially varying
([nonuniform filtering](nonuniform_filtering.md)); stage 8 adds the density
envelope for FSC masking and, when requested, the lag-by-one reference mask.
NU high-resolution shell extension is not used in ab initio.

## Symmetry

If the target point group `pgrp` differs from the starting group `pgrp_start`
(usually C1), stages up to 3 refine in the starting group. At stage 3 the
symmetry axis is found by maximizing the mutual consistency of the
symmetry mates of the volume's band-limited polar Fourier samples:

```text
T_g = samples rotated by M_g R_axis,     S = sum_g T_g,
score(R_axis) = (1/n_sym) sum_g Re<S, T_g> / sqrt(||S||^2 ||T_g||^2),
```

over 1000 spiral directions on the northern hemisphere at a low-pass of at
least 12 A, followed by simplex refinement of the 10 best from three restarts.
The map is symmetrized by averaging its mates in the found frame, the
orientations are transformed into that frame and reduced to the asymmetric
unit, and a symmetric reconstruction becomes the next reference. In
multi-state runs each state finds and applies its own axis.

## Multi-state modes

- `single`: one map.
- `independent`: `nstates` maps refined from random labels and independent
  random starts. Its default five-stage schedule stops at 6 A, before the
  local, NU, trailing, and automask stages, because its purpose is inspection.
- `docked`: one shared map through stage 5, then a balanced cohort of
  particles (capped at 100 000) is randomly split into states and refined with
  shared geometric neighborhoods, so the states begin from one consensus
  geometry and diverge only in density.

## Final reconstruction

The FSC during the run is diagnostic only: the halves share initialization
and the ab initio search is not gold-standard. At the end, a fresh
all-particle reconstruction is made at the original sampling with no staged
search or trailing controls, after a fill-in pass has given every active
particle at least one post-schedule pose. Docked mode first verifies that
every active particle received a post-split update.

## Rationale

- Starting from noise with random poses means the first reconstruction is
  the projection-slice average of the data itself, and the coarse bandwidth
  keeps that average from locking onto high-frequency coincidences.
- Angular sampling follows bandwidth: at 20 A a 4-degree grid oversamples
  the data, so `nspace = 500` is enough and the search is cheap when it is
  most stochastic.
- Searching in C1 first and imposing symmetry only after the map has a
  shape avoids the classic failure where a symmetric start point is an
  unbreakable fixed point of the alternation.
- Balanced probabilistic assignment in the middle stages does for
  projection directions what it does for 2D classes: it prevents a few
  directions from absorbing the data before the map is informative.

## Implementation

- Stage schedule and policies: `src/main/simple_abinitio_controller.f90`;
  bandwidth ladder in `src/main/simple_abinitio_utils.f90` and
  `src/utils/filter/simple_estimate_ssnr.f90`.
- Orchestration: `src/main/commanders/simple/simple_commanders_abinitio.f90`.
- Symmetry axis search: `src/main/volume/simple_volpft_symsrch.f90`,
  `src/main/simple_symanalyzer.f90`.
- Base estimator: `src/main/commanders/simple/simple_commanders_refine3D.f90`.
- Policy: `doc/policies/abinitio3D_policy.md`.
