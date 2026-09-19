# Nonuniform Filtering

## Problem

A single FSC-derived low-pass cutoff filters the whole map at the resolution
of its *average* voxel. A well-ordered core is then over-smoothed and a
flexible periphery is left with noise. Nonuniform filtering estimates a
spatially varying cutoff `c(v)` from the data themselves, so that each voxel
is filtered at the resolution at which the two independent half maps still
agree there.

The input is an unfiltered even/odd half pair on one grid. The output is a
label field `c(v)` over a candidate bank of cutoffs, the filtered halves and
merged map, and a local-resolution map.

## Model

Let `E` and `O` be the raw halves and `E_c`, `O_c` their low-pass versions at
candidate cutoff `c` from the static bank

```text
20, 15, 12, 10, 8, 6, 5, 4 A,
```

capped at FSC0.143/1.5 of the raw pair (only members coarser than that are
retained, never fewer than two); when ML regularization is active the
regularized half pair is one more member beside the finest retained one,
at the same Potts coordinate, so the unary alone decides between them.

If the true local resolution at voxel `v` is `c`, then `E_c(v)` predicts
`O(v)` up to noise, and vice versa; filtering finer than `c` lets through
noise that the other half does not share, and filtering coarser loses
signal that it does. The evidence for candidate `c` at `v` is therefore the
cross-half prediction error

```text
r1_c(v) = [E(v) - O_c(v)] / sigma(|v|),
r2_c(v) = [E_c(v) - O(v)] / sigma(|v|),
C_c(v)  = H(r1_c(v)) + H(r2_c(v)),
```

with the Huber loss `H(r) = r^2/2` for `|r| <= 1.345`, else
`1.345(|r| - 0.6725)`. The transition 1.345 gives 95 percent efficiency at
the Gaussian; the linear tail stops isolated large residuals (a stray
strong voxel) from dominating.

**Whitening.** The noise scale `sigma(r)` is a function of real-space radius,
estimated once from the raw difference `E - O`: in each radial shell,

```text
sigma_j = 1.4826 * median | (E-O)_j - median(E-O)_j |,
```

with gap filling, radial smoothing, and linear interpolation between shell
centers. A radial profile rather than a global scalar is required because
the gridding correction and any tapered solve support make reconstruction
noise grow toward the periphery; a global scale would put central and
peripheral residuals in different Huber regimes.

**Support.** All statistics are evaluated inside the sphere of diameter
`mskdiam`. The sphere, rather than a density mask, is used because it keeps
the candidate objective comparable everywhere and retains a solvent
population that later serves as a noise reference.

## Algorithm

1. **Candidate costs.** For each `c`, compute `C_c(v)` on the support and
   smooth it with a mask-normalized tent kernel of radius `min(1.5 c, 30 A)`.
   Smoothing at a scale proportional to `c` matches the spatial extent over
   which a cutoff of `c` is meaningful; voxels outside the support never
   influence in-support costs.
2. **Initial labels.** `c(v) = argmin_c C_c^smoothed(v)`.
3. **Ordered-label Potts smoothing.** Minimize over the 26-neighbor lattice

   ```text
   E(c) = sum_v C_{c(v)}(v) + sum_{v~w} phi( |rank(c(v)) - rank(c(w))| ) / deg(v),
   ```

   where `rank` is the index in the retained bank and `phi` is zero for
   adjacent ranks and a linear-quadratic hinge for larger jumps. The prior is
   *ordered*: a 20-to-15 A boundary costs less than a 20-to-4 A jump. Degree
   normalization keeps the regularization strength uniform at the sphere
   boundary. Iterated conditional modes with eight-color sweeps (no two
   concurrently updated voxels are neighbors) and tie preservation.
4. **Synthesis.** At each voxel take `E_{c(v)}(v)` and `O_{c(v)}(v)`; the
   merged map is their average; the local-resolution map stores `c(v)` in
   Angstrom inside the support and zero outside or beyond Nyquist.

When ML regularization is active, the regularized half pair (the
closed-form voxelwise Wiener shrinkage of the raw pair) joins the bank as
one more member beside the finest retained rung, at that rung's prior
coordinate, and competes for every voxel like a hard rung: its cost is the
cross-half prediction error of the regularized halves. Where it wins, the
map keeps the estimator's own high-resolution content; where a hard rung
wins, the estimator over-reached there. Its label resolution is the raw
pair's FSC=0.143, and it joins the bank the moment that is at or beyond the
ladder's finest rung; within the ladder the rungs compete alone (the cut at
`fsc/1.5` always keeps a rung finer than the FSC).

## High-resolution extension (retired 2026-09-16)

A sequential shell walk (`nu_refine=yes`) used to challenge the frontier of
the finest populated label with the next unrepresented Fourier shell. Its
acceptance criterion was equivalent to a local FSC above 0.5, conservative
for the map and no better than the FSC=0.143 extent for the matching band,
and it never let the regularized pair compete. A generated dense ladder
replaced it for two days (2026-09-16 to 18) and regressed PfCRT; the static
ladder with the regularized pair beside the finest retained rung is the one
competition for every workflow.

## Handoff to matching

The FSC and the NU filter answer different questions: the FSC reports the
average resolution, the NU field reports where the map is better than
average. The matching low-pass for the next iteration is the finest member
of the bank: the finest rung of the ladder cut at `fsc/1.5`, or the
regularized pair once its FSC=0.143 is at or beyond the ladder's finest
rung, bounded by any explicit `lp` and by `lpstop`. The same rule serves
every workflow, and which labels won voxels decides the filter, never the
band. In plain `nonuniform` mode the even and odd NU halves stay separate
references; in `nonuniform_lpset` the merged NU map is used with a single
band. No further low-pass is applied on top of an NU reference.

The band leads the FSC by 1.25-1.5x within the ladder, which is what
carries a climb. With the band at the previous FSC=0.143 crossing the
alignment never sees the shells beyond the crossing where the reference
still holds signal; the crossing then moves at most about one shell per
iteration and the orientation assignment crawls or freezes (PfCRT: 0/4 and
1/9 restarts converged, 2026-09-16/18). That is underfitting: the
sigma2-weighted objective, the stochastic assignment and the
evidence-limited NU reference guard against fitting noise, and the runs
with a 1.5x lead (PfCRT 10/10, +0.47 A of FSC per stage-6 iteration
against +0.19 A) show no sign of it. The uncapped July 2026 ladder led by
about 2x and froze 2 of 10 restarts, so the lead is bounded by the cut.

The reproducibility envelope derived from the same candidate costs is a
separate estimator: [NU-evidence envelope masking](nu_evidence_envelope_mask.md).

## Rationale

- Cross-half prediction error is a direct, model-free test of local
  resolution: it needs no assumption about the signal, only that the two
  halves have independent noise.
- Comparing candidates with one shared whitening profile makes the
  minimum-cost label meaningful; comparing candidates each smoothed at its
  own scale would bias the boundary.
- The ordered Potts prior encodes that resolution varies continuously in
  space, which is what distinguishes local resolution from voxelwise noise.

## Implementation

- Bank, costs, labels: `src/main/nu_filt/simple_nu_filter*.f90`; the per-state
  driver: `src/main/volume/simple_nu_state_filter.f90`.
- Noise scale and Huber objective: `src/main/image/simple_image_calc.f90`.
- Integration into volume assembly:
  `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Matching bandwidth handoff:
  `src/main/strategies/search/simple_matcher_refvol_utils.f90`.
- Policy: `doc/policies/nonuniform_filtering_policy.md`.
