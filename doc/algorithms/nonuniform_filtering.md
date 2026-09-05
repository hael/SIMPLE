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
candidate cutoff `c` from the bank

```text
20, 15, 12, 10, 8, 6, 5, 4 A.
```

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

When ML regularization is active, the regularized half pair may
conservatively replace the finest member of the bank, since it is already the
best available estimate at that cutoff.

## High-resolution extension (`nu_refine=yes`)

The bank's finest member bounds what can be selected. Voxels currently at the
finest populated label form a frontier; the next unrepresented Fourier shell
is offered as a challenger on that frontier only. It is accepted if at least
5 percent of the tested voxels (and an absolute minimum count) prefer it.
Accepted shells may continue outward; the process stops at the first
rejected, unsupported, or off-grid shell. The complete label field then
receives the same ordered smoothing, and the accepted depth is carried to the
next iteration so extension is monotone across a refinement.

## Handoff to matching

The FSC and the NU filter answer different questions: the FSC reports the
average resolution, the NU field reports where the map is better than
average. After filtering, the finest cutoff selected anywhere in the support
becomes the matching low-pass for the next iteration, bounded by any explicit
`lp` and by `lpstop`, so that particles are aligned against all the signal the
reference actually contains. In plain `nonuniform` mode the even and odd NU
halves stay separate references; in `nonuniform_lpset` the merged NU map is
used with a single band. No further low-pass is applied on top of an NU
reference.

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

- Bank, costs, labels, extension: `src/main/nu_filt/simple_nu_filter*.f90`.
- Noise scale and Huber objective: `src/main/image/simple_image_calc.f90`.
- Integration into volume assembly:
  `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Matching bandwidth handoff:
  `src/main/strategies/search/simple_matcher_refvol_utils.f90`.
- Policy: `doc/policies/nonuniform_filtering_policy.md`.
