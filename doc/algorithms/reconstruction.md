# Cartesian 3D Reconstruction

## Problem

Given particle images with fixed poses, shifts, CTFs, state labels (or
per-state weights, see the weighted variant below), half-set labels, and
per-shell noise powers, estimate the even and odd 3D volumes of each state. This is the linear inverse problem of tomography with a
frequency-dependent, sign-changing transfer function. Poses are treated as
known; their estimation is the concern of [refine3D](refine3d.md).

## Forward model

For volume `x` with Fourier transform `F`, particle `i` observes

```text
y_i = C_i S_i G_i F + n_i,       n_i ~ N(0, diag(sigma2_i(k))),
```

where `G_i` gathers the central section at orientation `R_i` from the 3D
Fourier grid with the Kaiser-Bessel (KB) kernel, `S_i` is the shift phase
ramp, and `C_i` is the CTF. The KB kernel is

```text
apod(u) = (1/W) I_0( beta sqrt(1 - (2u/W)^2) ),   |u| <= W/2,
```

with window `W = 3` grid points and `beta` chosen for oversampling factor
`alpha = 2`; the 3D Fourier grid is `2x` the native box so that the kernel's
real-space envelope falls outside the reconstructed volume. Particles are
grouped by state and half before accumulation, and point-group symmetry is
imposed by inserting every symmetry-related orientation of each particle.

## Gridding backend (default)

The weighted least-squares solution of the forward model, ignoring the
coupling between neighboring grid points, is the ratio of two accumulators:

```text
B(q) = sum_i sum_{g in sym} G_{ig}^T [ C_i S_i^* y_i / sigma2_i ]
D(q) = sum_i sum_{g in sym} G_{ig}^T [ |C_i|^2 / sigma2_i ].
```

`B` is the CTF-weighted backprojection of the whitened data, `D` the
CTF-squared sampling density. Each particle plane is deposited into the `3^3`
grid neighborhood of every polar sample with separable KB weights normalized
to unit sum per axis. Without ML regularization the `1/sigma2` factors are
absent. Only the `k <= 0` half is stored; Friedel symmetry supplies the rest.

Restoration:

1. **Density division.** `F(q) = B(q) / D(q)` inside the resolution limit,
   zero beyond. The division is bare: with ML regularization the
   Wiener term has already been added to `D` ([refine3D](refine3d.md)), and
   without it the density of a real dataset is never small enough inside the
   limit to need a floor. The one exception is a state reconstructed from
   fractional weights (below).
2. **Inverse FFT** to the padded real-space grid.
3. **Gridding correction.** Multiplication by the reciprocal of the KB
   kernel's real-space envelope, computed as the discrete Fourier transform of
   the normalized stencil so that it matches exactly what deposition did.
4. **Crop** to the native box.

The merged map is restored from `B_even + B_odd` and `D_even + D_odd` (with ML
regularization, after each half has received its own prior), not by averaging
restored halves; the FSC is computed on the halves before the gridding
correction, which is a pure amplitude envelope and cancels in the ratio.

## PCG backend (opt-in)

`rec_backend=pcg` solves the same model without the decoupling
approximation. With `K_i = sigma_i^{-1} C_i S_i G_i`, the normal equations are

```text
( sum_i K_i^T K_i + Lambda ) F = sum_i K_i^T sigma_i^{-1} y_i,
```

i.e. `(H + Lambda) F = B`, where `H` couples grid points through the KB
kernel overlap and `Lambda` is an optional prior. Because every `K_i^T K_i`
is a shift-invariant deposit of `|C_i|^2 / sigma2_i` at one orientation, `H`
is a convolution operator on the padded lattice whose kernel is the Fourier
transform of `D` deposited with the *squared* stencil. Applying `H` is then
one padded FFT, a pointwise product, and an inverse FFT, so the cost of an
iteration is independent of the number of particles: the particle loop runs
once, to build `B` and `D`, and never again.

The solve is left-preconditioned conjugate gradients with the reciprocal of
`D` plus a shell-relative floor as preconditioner (an absolute floor would
amplify unsampled modes by orders of magnitude). A soft spherical support
`P` may be imposed as `(P H P) u = P b`, which keeps the operator symmetric
positive semidefinite. Modes outside the reachable Fourier sphere stay zero.
Convergence is reported as the relative residual `||b - HF||/||b||` and the
relative update `||dF||/||F||`; the iteration count is kept short because the
Toeplitz kernel is an approximation of `H` whose positive-curvature region is
limited.

Two solves are made from one accumulation: a base solve producing
unfiltered halves, whose FSC remains the resolution authority, and an ML
replay that adds the FSC-derived `1/tau2` shell prior and warm-starts from the
previous same-half solution. PCG maps receive no gridding correction and no
second density division. Trailing reconstruction keeps its own chain of raw
`(B, D)` accumulators and blends it with the population rule of
[sampling](sampling_and_fractional_updates.md) before any prior is applied.
Shared-memory execution runs the same route as distributed execution, in one
process, as the only worker and then the master, so update fractions and
trailing chains behave the same in both modes.

## Weighted state reconstruction

`reconstruct3D m_estimator=flex` reconstructs each state from per-particle
weights instead of hard labels. The weights come from the project's state
weight set, which `flex_pca` publishes. Every selected particle enters each
state it has weight for, scaled by that weight. With `state=X` only state `X`
is reconstructed, from every particle with a positive weight for it; without
`m_estimator=flex`, `state=X` uses the particles labelled `X`. The state label
stays the selection flag, so a particle labelled 0 is not reconstructed; FLEX
labels each weighted particle with its largest-weight state. Both backends
take weights, in shared memory and distributed.

The state weight set (`src/main/project/simple_state_weight_set.f90`) is one
binary file per state, holding that state's weight in `[0, 1]` and a
hard-label flag for every particle, plus a text manifest. The manifest is
published last, by temporary name and rename, and the project's `out` segment
points at it. It records the generation, the kind, the producer, the particle
count, the particle layout digest, each file's size and checksum, and per
state the applied mass `sum w`, the effective sample size (ESS)
`(sum w)^2 / sum w^2` and the hard population. A set opens only when all of
this validates against the current project, so a half-written generation is
never read. The kind is `partition` when every weighted particle's weights sum
to one over the states, and `kernel` otherwise. The generation and the layout
digest together are the set's identity.

In gridding, `insert_plane_oversamp` takes an optional weight `w` that scales
the data term and the density term alike, so `B` and `D` above both carry the
factor `w_i(s)`. A particle is a member of state `s` when its weight is above
zero. Membership lists are built per state and half, and one half-map
reconstructor is held at a time, so memory does not grow with the number of
states. A particle is read once per state it weighs into. Partial files and
their summation are unchanged.

In PCG, a particle is a member of state `s` when its weight exceeds
`PCG_WEIGHT_THRESHOLD` (0.01, in `simple_reconstructor_pcg.f90`). Its noise
spectrum is divided by its weight, `sigma2_i / w_i(s)`, which multiplies its
contribution to both `B` and `D` by `w_i(s)`; the threshold keeps that
division well conditioned. The raw accumulator header records the applied
mass beside the particle count, and the identity of the weight set. The
master refuses a part built under another set.

Fractional weights make the density small and erratic in poorly occupied
regions. The gridding assembly (volassemble) therefore applies a
shell-relative floor `D(q) >= mean_shell(D) / 1000` to each half's density
before restoration, and to the merged density before the merged map is
restored, but only when the state's weights are fractional: some weight lies
strictly between 0 and 1. The floor is never applied on hard runs or with 0/1
weights. FLEX has no reconstruction code of its own; its trial maps come from
a weight table it owns rather than a published set, and they ask for the
floor through volassemble's internal `rho_floor=yes` key. PCG applies no extra
floor, since its preconditioner already floors the density shell by shell.

Weights of 0 and 1 equal to the hard labels reproduce the hard maps byte for
byte on both backends: membership equals the labels, a weight of one leaves
every product unchanged, and no floor is applied.

Under a weight set, population counts become masses. The state update counts
of the population rule, the population a trailing seed represents and the
mass a trailing chain represents are all applied masses: sums of the weights
as the backend applies them (after its membership threshold) over the same
particles that were counted before. The realized update fraction is the ratio
of two such masses. Gates
that decide whether a state is populated use the ESS instead of a head count:
the assembly's state populations, the PCG nonuniform-filter low-pass handoff
and the refinement's probability table. Post-processing skips a state without
weight. The weights are frozen for a trailing chain. Each chain manifest, for
gridding and PCG alike, records the identity of the weight set it was built
under, and a chain read under another identity is discarded and re-seeded,
never blended. The frozen accumulators of `solve3D_addon` record the same
identity and refuse another.

`refine3D_auto state=X m_estimator=flex` refines one state with these weights
in a work project. Its iterations and its final native reconstruction take the
same weighted path, and the final map is registered with its hard population,
applied mass and ESS as separate fields.

## Complexity

Accumulation is `O(N_particles x N_sym x N_polar x 27)`; restoration is one
3D FFT on the padded grid. Under a weight set, `N_particles` is the summed
state membership, since a particle counts once per state it weighs into.
Distributed execution reduces per-partition accumulators by summation, which
is exact.

## Rationale

- Depositing whitened, CTF-multiplied data and dividing by CTF-squared
  density is the least-squares estimator for a frequency-dependent transfer
  function: where the CTF is near zero for one defocus, particles at other
  defoci fill the gap in `D`, and no particle is ever divided by its own CTF.
- Kaiser-Bessel interpolation has near-optimal energy concentration for a
  given window, and its instrument function is known analytically, so the
  interpolation blur is removed exactly rather than approximately.
- Even/odd halves reconstructed from disjoint particles are the only way to
  obtain a resolution estimate whose noise terms are independent.
- Scaling a particle's data and density by its weight is the same as dividing
  its noise variance by that weight. The weighted map is therefore the
  least-squares estimator in which a particle of weight `w` counts as `w` of
  a particle, and both backends implement the same model.

## Implementation

- Gridding accumulator and restoration: `src/main/volume/simple_reconstructor.f90`;
  kernel in `src/main/interp/simple_kbinterpol.f90`, envelope in
  `src/main/interp/simple_gridding.f90`.
- Plane preparation (CTF, shift, whitening): `src/main/image/simple_image_ctf.f90`.
- Even/odd pair restoration (volassemble, with the fractional-weight density
  floor): `src/main/commanders/simple/simple_commanders_rec_distr.f90`.
- Reconstruction service (weight source, backend, dispatch):
  `src/main/strategies/parallelization/simple_rec3D_service.f90`; weighted
  membership and insertion in `src/main/strategies/search/simple_matcher_3Drec.f90`.
- State weight set: `src/main/project/simple_state_weight_set.f90`; file
  format in `src/fileio/simple_state_weights_file.f90`.
- Half-map diagnostics (FSC, cFAR): `src/main/volume/simple_halfmap_diagnostics.f90`.
- PCG operator and solver: `src/main/volume/simple_reconstructor_pcg.f90`;
  orchestration in `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`.
- Policies: `doc/policies/KB_Interpolation_Policy.md`,
  `doc/policies/3D/reconstruct3D_pcg_policy.md`.
