# Staged Ab Initio 2D Classification

## Problem

Construct stable 2D class averages and assignments from random
initialization. The alternating estimator in
[Cluster2D](cluster2d_class_averaging.md) converges to a local optimum, and
from a random start the nearest local optimum is poor: a few classes capture
most particles, high-frequency noise is fitted before low-frequency shape, and
shifts absorb misassignment. Ab initio 2D is a continuation scheme that steers
the alternation toward a good optimum by controlling three things over time:
the admitted bandwidth, the amount of stochasticity in the search, and the
number of particles updated per iteration.

## Continuation schedule

The run is a sequence of `cluster2D` stages (six by default, five when the
even/odd final stage is disabled), grouped into two phases: stages 1 to 4 are
the exploratory phase, stages 5 and 6 the consolidation phase.

**Bandwidth.** The low-pass limit starts coarse and halves its distance to the
target each stage:

```text
lp_start = clamp(mskdiam/12, 15, 20) A
lp_stop  = clamp(mskdiam/22,  6,  8) A
lp_1 = lp_start,   lp_i = lp_{i-1} - (lp_{i-1} - lp_stop)/2,   lp_{n-1} = lp_stop.
```

The defaults scale with the mask diameter because a larger particle has more
low-frequency structure to classify on. The final stage does not fix `lp` at
all: it lets the band follow the class FRCs, so the last iterations match at
whatever resolution the classes have actually reached.

**Search stochasticity.** Stages 1 and 2 use stochastic neighborhood hill
climbing with sampled in-plane and class choices (`snhc_smpl`), which visits
only a random fraction of classes per particle and accepts random improvements.
From stage 3 the search switches to the probabilistic mode: a global
class-likelihood table is built and particles are assigned to classes by a
balanced stochastic assignment ([sampling](sampling_and_fractional_updates.md)).
In `prob_snhc`, the intermediate stages keep that table sparse (only a
neighborhood of classes is scored per particle) and the last stage densifies it
so every particle's previous class is an explicit candidate and class-overlap
statistics are meaningful.

**Particle subset.** With `N_active` particles and a per-iteration target
`nsample` (default 200 000),

```text
update_frac = min(1, nsample / N_active),
```

active when below 0.99. Stage 1 uses one fixed random subset for all its
iterations (sticky), so that the random references are overwritten by a
consistent population; stage 2 begins carrying class sums forward
fractionally; from stage 3 the subset is redrawn every iteration, biased toward
particles with the fewest updates so coverage evens out. Stage iteration
counts are stretched when `update_frac` is small so that each particle is
expected to be visited at least once per phase:

```text
nits_coverage = clamp( ceil((1/update_frac - 1)/3), 1, 10 ).
```

**Regularization.** Stage 1 references are Gaussian-noise images (or random
particles, or averages of random labels), band-limited by a Gaussian at `lp_1`
rather than by an FRC-derived filter, with no shifts, no centering, and no ML
regularization. Any FRC between random halves would be meaningless there. From
stage 2 onward class averages are restored as MAP estimates with the
FRC-derived signal prior, and shifts are allowed up to
`trs = min(5, max(2, 12 A / smpd_crop))` pixels.

**Sampling grid.** The whole run works on one downscaled copy of the data,
targeting 2.67 A per pixel with a minimum box of 88 pixels; the field of view
is preserved. Nothing in the schedule needs finer sampling than `lp_stop`.

## Candidate sampling within a stage

In the probabilistic stages the per-particle candidate distribution over
classes is

```text
p_j = exp[-(d_j - d_min)] / sum_l exp[-(d_l - d_min)],
```

where `d_j` is the whitened Euclidean loss for class `j` at the seed shift
(for the correlation objective, `d = 1 - clamp(cc, 0, 1)`). Because the losses
are already noise-normalized log-likelihoods, no temperature is needed; the
`d_min` shift is only overflow protection. One candidate is drawn, its shift
is refined, and the matcher commits it as a hard assignment. The class
average is never a soft weighted sum over all candidates.

## Coverage and final products

Fractional updates leave some particles with stale poses. After any sampled
stage, the run therefore ends with a dense, greedy, all-particle pass with
fractional updates off. Every active particle's class, angle, and shift is
refreshed before the final class averages are restored, so the published
averages and FRCs describe the assignments that accompany them.

The run publishes merged, even, and odd class stacks, per-class FRCs, ranked
class products, and the per-particle assignments, scores, and noise spectra.

## Rationale

- Coarse-to-fine bandwidth is a homotopy: at 20 A the objective is smooth and
  has few minima, and the solution found there is a good starting point for the
  sharper objective at 15 A, and so on.
- Stochastic acceptance early and greedy acceptance late is annealing. The
  switch to balanced probabilistic assignment at stage 3 addresses the
  class-collapse failure mode directly: a class cannot take a second particle
  until it wins another draw, so populations stay comparable.
- Random noise references are deliberately uninformative, so the first stage's
  structure comes entirely from the data.
- Fractional updates make the cost per iteration independent of dataset size;
  the fill-in pass restores the guarantee that the final products are
  self-consistent.

## Implementation

- Stage schedule: `src/main/simple_abinitio2D_controller.f90`.
- Orchestration: `src/main/commanders/simple/simple_commanders_abinitio2D.f90`.
- Iteration execution: `src/main/strategies/parallelization/simple_cluster2D_strategy.f90`.
- Class-likelihood table: `src/main/simple_eul_prob_tab2D.f90`.
- Policy: `doc/policies/abinitio2D_policy.md`.
