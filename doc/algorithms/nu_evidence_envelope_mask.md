# NU-Evidence Envelope Masking

## Problem

Derive a soft molecular envelope from the same cross-half prediction errors
that drive [nonuniform filtering](nonuniform_filtering.md), so that the mask
selects voxels where the two half maps *agree* rather than voxels that are
merely dense. The estimator is opt-in (`nu_envmsk=yes`).

The envelope measures reproducible, locally ordered signal. It is not a density
threshold and is not the map of selected NU low-pass labels. This distinction
allows the method to reject strong but cross-half-inconsistent density, such as
disordered solvent or a detergent belt.

## Support and Candidate Bank

The even and odd half maps must have the same dimensions and sampling distance.
`mskdiam` defines a centered spherical support. This sphere is used for the NU
objective, evidence statistics, and envelope segmentation; neither a density
automask nor the resulting evidence envelope can replace it.

The static low-pass bank is

```text
20, 15, 12, 10, 8, 6, 5, 4 A.
```

Let `E` and `O` be the raw even and odd maps and let `E_c` and `O_c`
be the maps filtered with candidate `c`.

## Squared-Error Objective at One Noise Level

First, SIMPLE estimates one candidate-independent noise level from the raw
even-minus-odd values of the observed voxels inside the sphere (voxels where
both half maps are exactly zero are unobserved and left out), the
Gaussian-scaled MAD

\[
\sigma_0 = 1.4826\,\operatorname{median}
\left| (E-O)-\operatorname{median}(E-O) \right|.
\]

A median over the whole support is not inflated by the minority of voxels
where the halves disagree over density, which a locally estimated level would
be. If the MAD is degenerate the RMS supplies the fallback, and an all-zero
pair takes the unit level. For every candidate and supported voxel `v`, the
cross-half residuals are

\[
r_{1,c}(v)=\frac{E(v)-O_c(v)}{\sigma_0}, \qquad
r_{2,c}(v)=\frac{E_c(v)-O(v)}{\sigma_0},
\]

and the candidate cost is the plain squared error

\[
C_c(v)=\tfrac12\left[r_{1,c}(v)^2+r_{2,c}(v)^2\right].
\]

Reconstruction noise is Gaussian, so no robust loss is needed; the large
residuals come from signal that a too-coarse candidate removed and are
charged in full. The shared level makes candidates comparable at each voxel
and drops out of the per-voxel competition altogether; it enters the
evidence only through the margins below, which are thresholded across voxels.

## Evidence Margin

Before the candidate-specific smoothing used to select NU filter labels,
SIMPLE records the raw coarsest cost and the best raw cost:

\[
B(v)=C_{20\,\mathrm{A}}(v), \qquad
M(v)=\min_c C_c(v).
\]

The absolute improvement is

\[
D(v)=\max(0,B(v)-M(v)).
\]

`D` is embedded in the spherical support and smoothed once with a
mask-normalized 3D tent kernel. The smoothing radius is

\[
r_{\mathrm{smooth}}=
\min(1.5\,\texttt{amsklp},30\,\mathrm{A}),
\]

subject to a minimum of one sampling interval. Smoothing the difference once
is equivalent to smoothing baseline and best terms identically. It avoids the
boundary bias that would result from comparing costs smoothed at their
candidate-dependent NU scales.

The evidence value is the smoothed absolute improvement:

\[
e(v)=\widetilde D(v).
\]

An absolute margin is used, rather than a baseline-to-best ratio, because the
costs are already noise-normalized and candidate-independent; a ratio
would let a high-contrast core outvote weak but ordered density.

The margin does depend on the noise level, unlike the per-voxel competition:
where the noise is higher, as toward the box edge under the gridding
correction, the same true improvement yields a smaller margin. One global
level therefore costs some recall of peripheral density against a
spatially varying one; the measured trade-off and the alternative (dividing
the margins by the reconstruction's known geometric noise factor) are
recorded in the refactoring note `nu_euclidean_loss.md`.

## Robust Evidence Score

SIMPLE estimates the no-evidence population from the observed values inside
the sphere (voxels where both half maps are exactly zero are excluded):

\[
\mu_0=\operatorname{median}(e), \qquad
s_0=1.4826\,\operatorname{median}|e-\mu_0|.
\]

For `nu_msk_sig = k`, the evidence threshold and normalized score are

\[
t=\mu_0+k s_0, \qquad q(v)=\frac{e(v)-t}{s_0},
\]

with a numerical floor on `s_0`. Positive `q` favors signal and negative `q`
favors solvent.

This null estimate assumes solvent occupies most of the spherical support.
Because `e` is a best-of-bank statistic, solvent values are not exactly zero:
finite noise can make one candidate win by chance. The null therefore applies
to the exact candidate bank and smoothing configuration being evaluated.

An envelope-constrained base pair (a PCG solve on a density support) leaves
the far solvent unobserved, so the remaining support is not solvent-dominated.
`set_nu_evidence_null_shell` then designates the null geometrically: \(\mu_0\)
and \(s_0\) are taken on the dilation ring of the density envelope (dilated
minus core, where the base support has full weight), labels are free only on
the observed density envelope, and voxels outside it are fixed solvent.

## Binary MRF Segmentation

The initial binary field labels voxels with positive score as signal. SIMPLE
then minimizes a two-label Markov random field by iterative conditional modes.
For voxel `v`, let `d(v)` be the number of supported voxels in its
26-neighbor neighborhood and let `n_s(v)` be the number currently labeled as
signal. The two local energies are

\[
E_{\mathrm{signal}}(v)=-q(v)+\beta\frac{d(v)-n_s(v)}{d(v)},
\]

\[
E_{\mathrm{solvent}}(v)=q(v)+\beta\frac{n_s(v)}{d(v)},
\]

where the production policy fixes \(\beta=1\). The voxel takes the lower-energy
label. Degree normalization prevents voxels at the spherical boundary from
receiving a different effective regularization strength.

Updates use an eight-color 3D schedule, so voxels updated concurrently are not
26-neighbors. Iteration stops when no label changes or after six sweeps. This
prior regularizes boundary area but does not enforce connectivity.

## Topology, Morphology, and Softening

The binary MRF result is converted to the final soft envelope as follows:

1. find all 26-connected signal components;
2. retain every component whose size is at least 0.1 times the largest
   component size;
3. fill enclosed background cavities;
4. grow the binary field by 1 A; and
5. apply a cosine soft edge of 6 A.

Both lengths are converted to voxels from the map's sampling distance (at
least one voxel), so the physical finish is approximately constant across
samplings.

Keeping components relative to the largest, rather than keeping only one
component, permits separated ordered domains to survive when their linker is
not reproducible enough to pass the evidence threshold.

## Parameters and Defaults

| Parameter | Default | Effect |
|---|---:|---|
| `nu_envmsk` | `no` | Enable evidence-map and envelope generation |
| `mskdiam` | required | Diameter in Angstrom of spherical NU/evidence support |
| `amsklp` | `8` A | Physical evidence scale; sets the margin-smoothing radius |
| `nu_msk_sig` | `3.0` | Threshold in Gaussian-scaled MADs above the evidence median |

`nu_msk_sig` and `amsklp` are the only public envelope-shape tuning parameters.
The remaining algorithm parameters are fixed but retain explicit roles:

| Fixed choice | Value | Effect |
|---|---:|---|
| Scale-free evidence | no | A baseline-to-best ratio can prevent weak but ordered density from being outvoted by a high-contrast core; production uses the absolute squared-error cost margin |
| Density weight | 0.0 | Positive weight retains strong but poorly ordered density; zero keeps the mask evidence-only |
| MRF beta | 1.0 | Boundary smoothness; higher values give smoother boundaries |
| Minimum component fraction | 0.1 | Smallest connected component kept relative to the largest |
| Binary growth | 1 A | Expands the accepted binary support before softening |
| Cosine edge | 6 A | Softens the molecular-envelope boundary |

## Outputs and Interpretation

Two products accompany the NU-filtered maps:

- `*_nu_evidence.mrc`: the smoothed absolute evidence margin;
- `*_nu_envmask.mrc`: the component-filtered, hole-filled, dilated, soft
  envelope mask.

The evidence map is the primary diagnostic. A useful run should show a distinct
ordered-signal population rather than a continuously varying solvent field. If
the reported signal fraction exceeds 50%, the whole-support median/MAD estimate
cannot safely be interpreted as a solvent null; increase `mskdiam`, tighten the
threshold, or use a different null model before consuming the envelope.

The evidence envelope is selected from cross-half agreement, which makes it
unsuitable for FSC solvent correction: the FSC would then be computed inside a
region chosen for half-map agreement, biasing it upward. FSC masking uses an
independent density automask, and the NU objective keeps the spherical
`mskdiam` support.

## Implementation

- Command entry point and output naming:
  `src/main/commanders/simple/simple_commanders_resolest.f90`
- Candidate bank and raw evidence accumulation:
  `src/main/nu_filt/simple_nu_filter_bank.f90`
- Evidence calculation and binary MRF:
  `src/main/nu_filt/simple_nu_filter_envmask.f90`
- Noise level and squared-error objective:
  `src/main/image/simple_image_calc.f90`
- Component filtering, hole filling, dilation, and soft edge:
  `src/main/image/simple_image_msk.f90`
- Synthetic regression: the neutral fixture (a sphere of band-limited common
  signal in a noise support) is `test_evidence_envelope` of
  `src/main/nu_filt/simple_nu_filter_tester.f90` (`unit_reconstruction`); it
  reports the envelope's recall of true density and its solvent
  false-positive rate; the filter is exercised end to end through
  `simple_exec prg=nu_filt3D`

Historical design constraints: [superseded NU-evidence envelope masking note](../implementation_notes/rejected/nu_evidence_envelope_masking.md) (its reference-masking design was retired; the envelope itself is live).
