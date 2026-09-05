# Online Streaming Pipeline

## Problem

Process an acquisition while movies are still arriving, producing corrected
micrographs, CTF estimates, particle picks, and a global 2D classification
that improves as data accumulate, with the same estimators the batch
workflow uses. The online setting changes three things about the algorithms:
the dataset size is unknown and growing, so schedules cannot be planned in
advance; contamination must be rejected before it enters a global estimate
it would otherwise bias; and the cost per unit of new data must stay bounded.

## Stage graph

Data flow through six stages, each consuming completed outputs of the one
before:

```text
movies -> preprocessing -> optics groups -> initial 2D analysis
       -> reference picking and extraction -> particle sieve -> global pooled 2D.
```

Each stage runs the batch estimator on a bounded set of inputs and publishes
its result only when that set is finished; a stage never reads a partially
computed upstream result, and every input is consumed exactly once, so no
datum is counted twice in any downstream sum.

## Preprocessing and optics

Movies are grouped into small sets and each set is passed through
[motion correction](motion_correction.md),
[CTF estimation](ctf_estimation.md), and the segmentation
[picker](particle_picking.md). Micrographs failing the CTF-resolution
(10 A in streaming), ice-fraction (1.0), or astigmatism gates are excluded.
Accepted micrographs are assigned to optics groups from their acquisition
metadata; the group is the unit over which the noise model
and beam-tilt parameters are shared downstream.

## Bootstrap: initial 2D analysis

The segmentation picks are classified once with a small
[ab initio 2D](abinitio2d.md) run, `ncls = clamp(N/100, 10, 100)` classes at a
final low-pass of 8 A. Its selected class averages become the reference bank
for the reference picker. This is what turns a reference-free opening into a
reference-guided stream, and it is the only stage whose output is consumed
as a model rather than as data.

## Reference picking and extraction

With references available, each new micrograph is picked by exhaustive
Pearson matching against the bank (12 in-plane rotations per reference in
streaming), and particles are extracted at the box size implied by the
consensus diameter. Reference revisions replace the bank; micrographs picked
under an older bank are not repicked.

## Particle sieve

New particles enter a two-tier chunked classification whose purpose is to
reject junk before it reaches the global pool:

| | coarse chunk | fine chunk |
|---|---|---|
| particles per chunk | about 5 000 | about 10 000 |
| classes | 100 | 100 |
| sampled particles per iteration | 2 000 | 10 000 |
| final low-pass | 15 A | 10 A |

A coarse chunk closes when its unassigned particles exceed the threshold (the
micrograph that crosses it is included, keeping slices contiguous); fine
chunks merge the survivors of rejection-complete coarse chunks. Each chunk
runs the ab initio 2D schedule on a 128-pixel crop.

**Class rejection.** For every class average, 14 features are computed:
log population, log resolution, centering, log foreground and background
local variance and their ratio, an FRC proxy, center-to-edge SNR, connected
component area fraction, presence, band-pass variance in the 100 to 40 A band,
and a diffuse-signal score. Features are standardized per chunk by median and
MAD and clipped at four sigma. Hard gates remove classes with population below
0.35 percent of the chunk, bad pixels, no connected component, or density
extending beyond the mask. The remaining classes are scored by a logistic
model with pairwise feature interactions,

```text
eta = b_0 + sum_i b_i z_i + sum_{i<j} g_ij z_i z_j + rho * neighbor signal,
P(accept) = sigmoid(eta / T),
```

and accepted when `P >= 0.5`. Particles of rejected classes are removed from
the stream. The model rejects on shape statistics, not on resolution alone,
so early low-resolution chunks are not emptied.

## Global pooled 2D

Accepted particles are imported into a growing pool and classified with the
sampled stochastic-neighborhood [Cluster2D](cluster2d_class_averaging.md)
estimator, `ncls = 200` by default. Because the dataset grows, the
coarse-to-fine schedule is expressed in the pool's own iteration count rather
than in stages:

```text
lp(it) = lp_stop + (lp_start - lp_stop) (20 - it) / 20,     it < 20,
```

with Gaussian-limited references and no shifts for the first 5 iterations,
then a 5-pixel shift bound; the extremal annealing runs over the same 20
iterations. From iteration 20 the low-pass follows the class FRCs. The
first iteration waits until at least `max(20 ncls, 500 x arrival rate)`
particles are present, and the pool pauses when no new particles have
arrived, resuming when the count has grown by a rate-dependent margin.
Beyond 500 000 particles the update becomes fractional,
`update_frac = (N_old - N_previously_classified) / N_old`, so the cost per
iteration is bounded by the newly arrived data plus a fixed sample.

Once the pool is at its class limit, has been at Nyquist for five
consecutive iterations, and a larger box would materially improve Nyquist,
the crop is enlarged and `lp_stop` lowered accordingly. Snapshot
requests copy the current class products without disturbing the live
estimate.

## Rationale

- Running the unchanged batch estimators on bounded sets keeps every
  streaming result reproducible offline; the only online-specific logic is
  scheduling and rejection.
- Rejecting at the chunk level, with a model trained on class-average
  features, removes contamination when it is cheap to identify (in a small,
  homogeneous chunk) rather than after it has been absorbed into a
  200-class pool.
- Expressing the pool's continuation schedule in iterations rather than
  stages is what lets the same coarse-to-fine logic work when the dataset
  has no known final size.

## Implementation

- Stages: `src/main/stream/simple_stream_p01_preprocess_new.f90` through
  `src/main/stream/simple_stream_p06_pool2D_new.f90`,
  `src/main/stream/simple_stream_p03_initial_analysis.f90`.
- Sieve and chunking: `src/main/sieve/simple_ptcl_sieve.f90`.
- Class-quality features and model: `src/main/cavg_quality/simple_cavg_quality_*.f90`.
- Pool schedule: `src/main/stream/simple_stream_pool2D_utils.f90`.
- Diameter consensus: `src/main/stream/simple_mini_stream_utils.f90`.
- Policies: `doc/policies/sieving_and_rejection/`.
