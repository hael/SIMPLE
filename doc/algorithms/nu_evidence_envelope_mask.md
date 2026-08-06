# NU-Evidence Envelope Mask in `nu_filt3D`

## Purpose

`nu_filt3D` can derive a soft molecular envelope from the same cross-half
prediction errors used by the nonuniform filter. Enable it with:

```text
nu_envmsk=yes
```

The envelope measures reproducible, locally ordered signal. It is not a density
threshold and is not the map of selected NU low-pass labels. This distinction
allows the method to reject strong but cross-half-inconsistent density, such as
disordered solvent or a detergent belt.

## Support and Candidate Bank

The even and odd half maps must have the same dimensions and sampling distance.
`mskdiam` defines a centered spherical support. This sphere is used for the NU
objective, evidence statistics, and envelope segmentation; neither a density
automask nor the resulting evidence envelope can replace it.

The standalone command evaluates the static low-pass bank

```text
20, 15, 12, 10, 8, 6, 5, 4 A.
```

Let `E` and `O` be the raw even and odd maps and let `E_c` and `O_c`
be the maps filtered with candidate `c`.

## Noise-Normalized Huber Objective

First, SIMPLE estimates one candidate-independent half-map noise scale from the
raw even-minus-odd values inside the sphere:

\[
\sigma = 1.4826\,\operatorname{median}
\left| (E-O)-\operatorname{median}(E-O) \right|.
\]

An RMS fallback is used if this MAD estimate is numerically zero. For every
candidate and supported voxel `v`, the cross-half residuals are

\[
r_{1,c}(v)=\frac{E(v)-O_c(v)}{\sigma}, \qquad
r_{2,c}(v)=\frac{E_c(v)-O(v)}{\sigma}.
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

The common noise scale makes candidates comparable. The Huber loss remains
quadratic near the expected noise level but prevents isolated large residuals
from dominating the evidence.

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

By default the evidence value is the smoothed absolute improvement:

\[
e(v)=\widetilde D(v).
\]

With `nu_msk_rel=yes`, SIMPLE also smooths `B` at the same scale and uses a
scale-free cost-improvement ratio. Define

\[
f=0.1\,\operatorname{median}_v \widetilde B(v),
\]

with a positive numerical floor, then

\[
B_f(v)=\max(\widetilde B(v),f),\qquad
M_f(v)=\max(B_f(v)-\widetilde D(v),f),
\]

and

\[
e(v)=\frac{B_f(v)}{M_f(v)}-1.
\]

The floor prevents near-zero costs from producing arbitrarily large ratios.
Unlike a bounded fractional reduction, this statistic remains compatible with
the unbounded robust threshold used by the segmentation.

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

## Optional Density Rescue Term

When `nu_msk_dens` is nonzero, `nu_filt3D` forms the even/odd average, filters
it to `amsklp`, and robustly standardizes its density inside the sphere:

\[
z_\rho(v)=\frac{\rho(v)-\operatorname{median}(\rho)}
{1.4826\,\operatorname{median}|\rho-\operatorname{median}(\rho)|}.
\]

The segmentation score becomes

\[
q'(v)=q(v)+\texttt{nu_msk_dens}\,z_\rho(v).
\]

The default weight is zero. A positive value can retain strong but poorly
ordered peripheral density, but it weakens the method's ability to distinguish
ordered molecular signal from disordered density.

## Binary MRF Segmentation

The initial binary field labels voxels with positive score as signal. SIMPLE
then minimizes a two-label Markov random field by iterative conditional modes.
For voxel `v`, let `d(v)` be the number of supported voxels in its
26-neighbor neighborhood and let `n_s(v)` be the number currently labeled as
signal. The two local energies are

\[
E_{\mathrm{signal}}(v)=-q'(v)+\beta\frac{d(v)-n_s(v)}{d(v)},
\]

\[
E_{\mathrm{solvent}}(v)=q'(v)+\beta\frac{n_s(v)}{d(v)},
\]

where `nu_msk_beta` supplies \(\beta\). The voxel takes the lower-energy label.
Degree normalization prevents voxels at the spherical boundary from receiving
a different effective regularization strength.

Updates use an eight-color 3D schedule, so voxels updated concurrently are not
26-neighbors. Iteration stops when no label changes or after six sweeps. This
prior regularizes boundary area but does not enforce connectivity.

## Topology, Morphology, and Softening

The binary MRF result is converted to the final soft envelope as follows:

1. find all 26-connected signal components;
2. retain every component whose size is at least `nu_msk_minvol` times the
   largest component size;
3. fill enclosed background cavities;
4. grow the binary field by `binwidth` voxel layers; and
5. apply a cosine soft edge of `edge` voxels.

Keeping components relative to the largest, rather than keeping only one
component, permits separated ordered domains to survive when their linker is
not reproducible enough to pass the evidence threshold.

## Parameters and Defaults

| Parameter | Default | Effect |
|---|---:|---|
| `nu_envmsk` | `no` | Enable evidence-map and envelope generation |
| `mskdiam` | required | Diameter in Angstrom of spherical NU/evidence support |
| `amsklp` | `8` A | Evidence smoothing scale and optional density low-pass |
| `nu_msk_sig` | `3.0` | Threshold in Gaussian-scaled MADs above the evidence median |
| `nu_msk_beta` | `1.0` | Binary MRF boundary regularization |
| `nu_msk_dens` | `0.0` | Optional robust density-score weight |
| `nu_msk_rel` | `no` | Use the scale-free baseline-to-best cost ratio |
| `nu_msk_minvol` | `0.1` | Minimum component size relative to the largest |
| `binwidth` | `1` voxel | Binary dilation before softening |
| `edge` | `6` voxels | Cosine soft-edge width |

## Outputs and Interpretation

With envelope generation enabled, `nu_filt3D` writes two additional products
beside its NU-filtered maps:

- `*_nu_evidence.mrc`: the smoothed absolute evidence margin, or the
  dimensionless cost-improvement ratio when `nu_msk_rel=yes`;
- `*_nu_envmask.mrc`: the component-filtered, hole-filled, dilated, soft
  envelope mask.

The evidence map is the primary diagnostic. A useful run should show a distinct
ordered-signal population rather than a continuously varying solvent field. If
the reported signal fraction exceeds 50%, the whole-support median/MAD estimate
cannot safely be interpreted as a solvent null; increase `mskdiam`, tighten the
threshold, or use a different null model before consuming the envelope.

The evidence envelope is selected from cross-half agreement. It must not be
used for FSC solvent correction, and it does not define the NU objective
support. These roles remain assigned to an independent density automask and the
spherical `mskdiam` support, respectively.

## Implementation Locations

- Command lifecycle and output naming:
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

Design constraints and deferred workflow integration are discussed separately
in [`NU-Evidence Envelope Masking`](../implementation_notes/nu_evidence_envelope_masking.md).
