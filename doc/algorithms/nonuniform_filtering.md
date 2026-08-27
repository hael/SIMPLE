# Nonuniform Filtering

## Problem

Construct reference maps whose local bandwidth follows reproducible half-map
signal rather than one global FSC cutoff. The result is a spatially varying
low-pass field, filtered even/odd/merged maps, and a matching-bandwidth handoff.

## Support and inputs

The filter consumes unfiltered even and odd volumes on the same grid. A sphere
derived from `mskdiam` is the immutable objective support. It is deliberately
not replaced by a density mask or evidence envelope: the broad sphere keeps the
candidate objective comparable and retains a solvent population for evidence
calibration.

When ML regularization is active in static mode, a compatible auxiliary half
pair may conservatively replace the finest discrete bank member. High-resolution
extension instead owns its own challenger path and does not use that auxiliary
replacement.

## Candidate bank and whitening

The static low-pass bank is

```text
20, 15, 12, 10, 8, 6, 5, 4 A.
```

For candidate `c`, let `E_c` and `O_c` be the low-pass-filtered halves. The raw
even/odd difference estimates a radial noise profile `sigma(r)` by shell-wise
Gaussian-scaled MAD, with gap filling, smoothing, and voxelwise radial
interpolation. Radial whitening is required because deapodization and tapered
solve support make reconstruction noise spatially non-stationary.

The cross-half residuals at voxel `v` are

```text
r1_c(v) = [E(v)-O_c(v)] / sigma(radius(v))
r2_c(v) = [E_c(v)-O(v)] / sigma(radius(v)).
```

Their unary cost is the sum of Huber losses with transition `1.345`. Each
candidate cost is smoothed over the sphere with a mask-normalized tent kernel
whose radius is `min(1.5*c,30 A)`. Values outside support cannot influence an
in-support cost.

## Local label field

The initial label at each supported voxel minimizes the smoothed unary. SIMPLE
then applies an ordered-label Potts prior on the 26-neighbor lattice. Neighbor
penalties operate on retained-bank coordinates, tolerate adjacent resolution
steps, penalize larger jumps with a linear-quadratic hinge, normalize by the
number of supported neighbors, and preserve the current label on numerical
ties. Eight-color sweeps make concurrent updates non-neighboring.

This is not an unordered categorical smoother: a 20-to-15 A boundary is less
severe than a 20-to-4 A jump.

## Output synthesis

At each voxel, the selected label chooses the corresponding candidate value in
the even and odd banks. The merged filtered map is the average of the filtered
halves. A local-resolution map stores the selected Angstrom value inside
support and zero outside support or beyond Nyquist.

The base reconstructed halves remain primary data products. `_nu_filt` and
`_nu_locres` are derived references and diagnostics.

## High-resolution extension

With `nu_refine=yes`, voxels currently assigned to the finest populated label
form a frontier. The next unrepresented Fourier shell is challenged only on
that frontier. It is accepted when at least 5 percent and an absolute minimum
number of tested voxels prefer it. Accepted contiguous shells may continue;
the process stops at the first rejected, unsupported, or off-grid shell.

The active extension bank is thinned and capped for memory. After accepted
steps, the complete label field receives the same ordered-label cleanup. The
accepted depth is persisted for the next iteration.

## Matching handoff

FSC reporting and NU filtering have different bandwidth roles. FSC does not
truncate the candidate bank. After filtering, the finest cutoff selected
anywhere inside the NU sphere becomes the project-level matching `lp` for the
next iteration, bounded by explicit user `lp` and `lpstop` policy.

Plain `nonuniform` preserves independent even/odd matching and prefers the
corresponding NU halves. `nonuniform_lpset` promotes the selected limit and uses
the merged NU reference under LP-set topology. Ordinary low-pass filtering is
not applied again on top of an NU reference.

The reproducibility envelope derived from the candidate bank is a separate
algorithm: see [NU-evidence envelope masking](nu_evidence_envelope_mask.md).

## Implementation

- Filter bank and labels: `src/main/nu_filt/simple_nu_filter*.f90`.
- Volume lifecycle: `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Matching-reference selection:
  `src/main/strategies/search/simple_matcher_refvol_utils.f90`.
- Policy: `doc/policies/nonuniform_filtering_policy.md`.

