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

## Radially Whitened Huber Objective

First, SIMPLE estimates one candidate-independent radial noise profile from the
raw even-minus-odd values inside the sphere. Supported voxels are grouped by
real-space radius, and each radial shell receives the Gaussian-scaled MAD

\[
\sigma_j = 1.4826\,\operatorname{median}
\left| (E-O)_j-\operatorname{median}(E-O)_j \right|.
\]

Empty or numerically degenerate shells are filled from valid neighbors and the
profile is smoothed radially. If no shell has a valid MAD, a global MAD, then
RMS, supplies the fallback. Linear interpolation between shell centers gives
`sigma(r(v))` at voxel `v`.

This whitening is candidate-independent but spatially varying. It accounts for
radial noise changes introduced by reconstruction deapodization and tapered
solve support; a single global scale would put central and peripheral residuals
in different Huber regimes. For every candidate and supported voxel `v`, the
cross-half residuals are

\[
r_{1,c}(v)=\frac{E(v)-O_c(v)}{\sigma(r(v))}, \qquad
r_{2,c}(v)=\frac{E_c(v)-O(v)}{\sigma(r(v))}.
\]

The candidate cost is

\[
C_c(v)=H_{1.345}(r_{1,c}(v))+H_{1.345}(r_{2,c}(v)),
\]

where

\[
H_\delta(r)=
\begin{cases}
\tfrac12 r^2, & |r|\le\delta,\\
\delta\left(|r|-\tfrac12\delta\right), & |r|>\delta.
\end{cases}
\]

The common radial profile makes candidates comparable at each voxel. The Huber
loss remains quadratic near the expected local noise level but prevents
isolated large residuals from dominating the evidence.

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
Huber costs are already noise-normalized and candidate-independent; a ratio
would let a high-contrast core outvote weak but ordered density.

## Robust Evidence Score

SIMPLE estimates the no-evidence population from all values inside the sphere:

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
| Scale-free evidence | no | A baseline-to-best ratio can prevent weak but ordered density from being outvoted by a high-contrast core; production uses the absolute Huber-cost margin |
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
- Noise scale and Huber objective:
  `src/main/image/simple_image_calc.f90`
- Component filtering, hole filling, dilation, and soft edge:
  `src/main/image/simple_image_msk.f90`
- Synthetic regression:
  `production/tests/simple_test_nu_envmask.f90`

Design constraints: [nu_evidence_envelope_masking.md](../implementation_notes/nu_evidence_envelope_masking.md).
