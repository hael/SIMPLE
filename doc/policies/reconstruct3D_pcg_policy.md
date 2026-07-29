# `reconstruct3D_pcg` policy

This document is the contract of the code that is executed today. It describes
the experimental CTF- and sigma-weighted preconditioned-conjugate-gradient
(PCG) reconstruction path: what it does, the design invariants it depends on,
and what it deliberately does not yet support. Forward-looking design work is
collected in section 10.

`reconstruct3D_pcg` is an isolated prototype. It shares no code with the
production Fourier-gridding path (`reconstructor`, `reconstructor_eo`,
`volassemble`, the online matcher), never writes to the project, and does not
add a mode switch to any production program. This is a hard boundary, kept so
that the comparison against production reconstruction stays honest and so that
a bug in this path cannot affect production output.

## 1. Scope and fixed inputs

`reconstruct3D_pcg` reads a project and reconstructs exactly one volume for one
state by solving a CTF- and noise-weighted least-squares problem with PCG. For
fixed orientations, states, and shifts it answers one question: does solving
the weighted normal equations exactly improve a controlled reconstruction
enough to justify its cost, relative to the diagonal gridding approximation?

Every per-particle input is read-only and taken from the same metadata the
production reconstructor uses:

- orientation and in-plane shift are inputs, never optimized;
- exactly one state is reconstructed (`state`, default 1);
- particles are selected by `sample4rec` (state and, where available,
  `updatecnt` gated), then filtered to the requested state;
- CTF parameters, phase shift, and `ctfflag` come from the project;
- per-particle Euclidean `sigma2` is read from the existing estimation
  machinery (section 7);
- `mskdiam`, if given, constrains the solve to a soft spherical support
  (section 4).

The command requires `nparts=1` and `pgrp=c1`; both are enforced with a hard
error rather than silently ignored. No even/odd split, symmetry, or distributed
execution is performed (section 9).

Nothing is written back to the project. Output goes to a new numbered execution
directory: the reconstructed map
(`reconstruct3D_pcg_state<NN>.mrc`) and a solver diagnostics table
(`reconstruct3D_pcg_diagnostics.txt`) with iteration count, per-iteration
relative residual, and timing.

## 2. Where the code lives

| File | Contents |
| --- | --- |
| `src/main/volume/simple_reconstructor_pcg.f90` | `reconstructor_pcg` operator/solver type |
| `src/main/commanders/simple/simple_commanders_reconstruct3D_pcg.f90` | `reconstruct3D_pcg` command |
| `src/main/ui/simple/simple_ui_volume.f90` | UI record |
| `src/main/exec/simple_exec_volume.f90` | exec routing |
| `src/main/params/simple_parameters*.f90` | `pcgop`, `rtol` parameters |
| `src/fileio/simple_sigma2_files.f90` | shared, builder-free sigma2 discovery/loading |
| `src/main/commanders/test/simple_commanders_test_highlevel.f90` | operator/solver tests |

## 3. Statistical model

Let `x` be a real volume and `F` the 3-D Fourier transform. For particle `i`:
`G_i` extracts an oriented Cartesian Fourier plane by Kaiser-Bessel (KB)
interpolation, `S_i` is the shift phase, `C_i` is the complex CTF (including the
phase-flip convention), and `N_i` is the diagonal noise covariance from the
particle's `sigma2` spectrum. The whitened observation operator is
`K_i = N_i^{-1/2} C_i S_i G_i F`, and the solver targets the MAP normal
equations

```text
H x = b,    H = sum_i K_i^dagger K_i + Lambda,    b = sum_i K_i^dagger N_i^{-1/2} y_i
```

with `Lambda = lambda*I` a fixed positive prior held constant for the whole
solve. No current-volume-dependent mask, FSC estimate, filter, or nonlinear
clipping appears inside the operator; that linearity is what PCG relies on.

The normal operator carries only the real weight `|T_i|^2 = |C_i|^2 / sigma2_i`:
the shift phase is unit-modulus and cancels between the forward and adjoint
passes of `K_i^dagger K_i`. The full complex transfer
`T_i = C_i S_i / sqrt(sigma2_i)` is needed only for the right-hand side `b`,
whose adjoint uses `conjg(T_i)` (required for shifted particles and phase
shifts, even when the CTF itself is real).

## 4. Numerical design and invariants

These are the load-bearing conventions of the operator. Changing any of them
requires re-passing the operator tests in section 8.

**Full-symmetric-disk planes.** `forward_plane`/`adjoint_plane_add` always use
an unpacked, both-sign-`h` disk, never the packed half-plane. This makes them
an exact adjoint pair for any orientation at the cost of 2x plane work. The
packed half-plane's Nyquist-bin bookkeeping is an orientation-dependent silent
energy-loss trap that a single-orientation adjoint test does not catch; the
full disk removes that surface permanently. The KB window is used at its
natural odd width (`2*iwinsz+1`, 3 points for `KBWINSZ=1.5`); an even width
would not be centred on `nint(loc)` and would break the mirror consistency the
full-disk fold depends on.

**Oversampling.** The unknown lives on the native box, but every Fourier
operation runs on a `OSMPL_PAD_FAC`-times padded lattice: the volume is
centre-padded going in and centre-cropped coming out, with `pad_vol`/`crop_vol`
as exact adjoints. This is the same arrangement Fourier gridding uses. Without
it, KB interpolation is only percent-accurate, the roll-off envelope swings
~30x across the box, and the Gram kernel cannot be made to match the operator.
Because the padded lattice is already `2*box`, it *is* the grid the linear
(non-wrapping) Gram convolution needs, so the kernel and the operator share one
lattice.

**Two KB envelopes, both measured.** KB interpolation multiplies the volume by
an instrument envelope in real space. Two arise and both are handled:

| Envelope | Origin | Handling |
| --- | --- | --- |
| Gather | reading the volume through the KB window | `build_env`, divided out by `deapod_mul` |
| Deposition | laying `|T|^2` onto the kernel grid | measured and divided out in the kernel build |

Both are measured by scattering a unit spike through the operator's own stencil
and transforming back, not derived from `kbinterpol%instr` (the continuous
transform, which disagrees with the three discrete renormalized weights by
about a factor of two at the box edge). The matrix-free operator applies the
gather envelope twice, so it is `A = A_env^{-1}(E T E)A_env^{-1}` with `E` the
gather envelope; `deapod_mul` brackets the operator to recover the true `T`.
For real particles this correction is not optional: synthetic data generated by
the operator carries the same envelope and cancels it (an inverse crime), but
real particles do not, so an uncorrected solve returns `E^{-1} x_true`.

**Soft spherical support.** When `mskdiam` is set, the solve is constrained to a
soft sphere by writing `x = P u` and minimising `||A P u - y||^2`, giving
`(P H P) u = P b`. `P` is applied on both sides of the operator and once to the
RHS. `P H P` stays symmetric positive semidefinite and its null space (outside
the support) is never entered from `x = 0`. Besides removing edge artefacts
where deapodization amplifies hardest, this shrinks the problem: at `mskdiam`
180 in a 256 box the sphere is ~18% of the volume, so most unknowns are solvent
the data barely constrains. The edge profile is production's own
`cosedge_r2_3d`.

## 5. Preconditioner and kernelized operator

**Preconditioner.** `M(k) = rho(k) + floor(shell)`, where
`rho = sum_i G_i^dagger |T_i|^2 G_i` is the gridding sampling density. The floor
is a fixed fraction of the shell-mean `rho`, *not* an absolute constant: `rho`
spans many orders of magnitude and is genuinely zero between rotated planes, so
an absolute floor would amplify the least-constrained modes by six to nine
orders of magnitude and PCG would fill the map with noise. Modes beyond
`padf*Rnat` are unconstrained and zeroed, making `M` singular but positive
semidefinite, which keeps the Krylov space out of the null space entirely. The
envelope belongs in `M` too: the operator being solved is `E^{-1} T E^{-1}`, so
`M^{-1}` brackets its Fourier divide with `E` (not `E^{-1}`).

**Kernelized (Toeplitz/Gram) operator.** Opt-in via `pcgop=kernel` (the
default). It replaces the per-iteration particle loop with one padded FFT, a
pointwise multiply by a precomputed real kernel `Khat`, and an inverse FFT, so
per-iteration cost is independent of particle count (~7x faster per iteration).
`Khat` is built with the standard NUFFT Gram construction: scatter `|T_i|^2`
onto the `2x` oversampled grid at doubled coordinates. The oversampling is the
mechanism, not a safety margin — it resolves the sub-pixel frequency offsets
that distinguish the operator from gridding. (A literal impulse-response kernel
does not work: KB weights are normalised to sum 1, so it reduces exactly to
gridding and PCG would converge in one step to the gridding answer.) The kernel
scale is calibrated once by least squares against the matrix-free operator on a
deterministic probe.

The kernelized operator is shift-invariant and the true operator is not, so it
remains an approximation (~3.4% interior error against the matrix-free
reference). The matrix-free path (`pcgop=matrixfree`) is the exact reference and
the fallback. The kernel should not be trusted as equivalent until it is shown
to give equivalent reconstructions on real data, not merely a similar operator.

## 6. PCG solver

Standard left-preconditioned CG of `H x = b` (`P H P` under a support
constraint). Zero initialisation is the baseline. All real-volume dot products
use a deterministic double-precision reduction. The solve fails explicitly on
non-finite or non-positive `dot(p,Hp)` and on a non-positive-definite
preconditioner or zero RHS, rather than continuing past a broken assumption.

The reported headline and stopping test are the **true** relative residual
`||r||_2 / ||b||_2`, not the preconditioned M-norm (which is not monotone under
a singular preconditioner and reads like a diverging solver while the solve is
in fact converging). The M-norm is still logged as the honest diagnostic for
the preconditioner. The recurrence residual is periodically audited against a
recomputed `b - Hx`; they agree to six significant figures, which is what
licenses reporting the cheap recurrence norm.

## 7. Sigma2 handling

Sigma2 weighting is gated on `objfun`, matching `reconstruct3D`: `objfun=cc`
runs unweighted (`sigma2 = 1`), `objfun=euclid` requires sigma2 files and
weights the fit by them. A missing sigma2 file under `objfun=euclid` is a hard
error, never a silent unweighted fallback — a silent fallback would quietly
change which objective is minimised.

Sigma2 is per-particle-per-shell, read via `euclid_sigma2` over the group STAR
file (as `flex_analysis` does), then upsampled to the operator's shell range.
The discovery/carry-over/loading logic lives once in
`src/fileio/simple_sigma2_files.f90`, called by both `flex_analysis` and
`reconstruct3D_pcg`. That module is deliberately builder-free: `builder`
depends on `euclid_sigma2`, so the sigma2 side must not depend back on
`builder`; callers pass `pftc`, `esig`, and orientations explicitly.

## 8. Tests

All operator/solver tests are in
`simple_commanders_test_highlevel.f90` and run in memory with no project I/O.
Each is a release gate for the ones that build on it; a visually plausible map
is not evidence that the CTF/sigma adjoint is correct.

| Test | Gate |
| --- | --- |
| `test=pcg_recon` | single gate, eight fail-fast stages: (1) adjoint dot-product with `T_i=1`; (2) the same with nonzero shift, astigmatic CTF and sigma2, isolating `build_transfer`; (3) normal-operator symmetry and positive-definiteness; (4) heterogeneous phantom recovery; (5) kernelized-vs-matrix-free equivalence, all-voxel and interior; (6) kernel shift-invariance, CTF-dependence, and the preconditioner; (7) streaming batch accumulation reproducing the monolithic solve; (8) deapodization against envelope-free data — the one stage that avoids the inverse crime |

The CTF-free and hetero tests generate their observations with `forward_plane`,
so the gather envelope cancels; they gate operator *algebra* only. Envelope
correctness for real particles is gated solely by stage 8 of `pcg_recon`.

## 9. Deliberate non-goals (current path)

The following are not implemented and are hard-errored or absent, not silently
approximated:

- no even/odd half-set split (section 10);
- no point-group symmetry — `pgrp=c1` is enforced (section 10);
- no `nparts>1` / distributed execution;
- no orientation, shift, state, sigma, FSC, or resolution write-back;
- no orientation search, probabilistic assignment, online partial
  reconstruction, fractional updates, or multi-state reconstruction;
- no adaptive regulariser, automask, nonuniform filtering, or FSC feedback
  inside PCG;
- no MPI/distributed reduction and no GPU/offload path;
- no reuse of SPIDER BP-CG real-space code (the design is Fourier
  central-section; the architectures do not transfer — there is no licensing
  barrier, SIMPLE is GPL-3.0 and SPIDER is GPL-2.0-or-later).

## 10. Deferred and future work

**Next: even/odd cross-validation and state assignment.** Getting the path
usable in real workflows requires independent even/odd half-set solves so an
FSC can be measured and the result compared against the frozen gridding
reference (visual map, FSC, local anisotropy indicators, objective trajectory,
wall time, peak memory). The command can already be run twice with explicit
even/odd particle selections; the missing pieces are half-set ownership in the
selection/solve path, the FSC comparison against gridding, and multi-state
assignment so more than one volume can be reconstructed per run. Half-set
independence is a hard requirement: no cross-half fitted volume may be used to
influence a particle's parameters.

**Deferred: point-group symmetry.** Real workflows need symmetry, which the
current `pgrp=c1` guard forbids. Symmetry must be handled either by expanding
the particle set over the asymmetric unit or by symmetrising inside the
operator; the choice interacts with the kernelized operator's shift-invariance
assumption and needs its own design pass.

**Longer term.** The validated operator provides one consistent weighted
likelihood for both reconstruction and continuous parameter refinement (poses,
shifts, restricted CTF/phase), which would be a separate command/strategy with
its own design review — block-coordinate or variable-projection refinement over
a fixed-volume PCG solve, staged one parameter family at a time, with
trust-region acceptance, parameter bounds/priors, and half-set safeguards
against the alignment-vs-blur and CTF-vs-sharpening degeneracies. Minibatch
stochastic optimisation of the same likelihood is a still-later option for very
large data sets, with PCG remaining the reference volume solver. None of this
is in scope for the current fixed-pose command.

## 11. Reused vs. reimplemented

`reconstruct3D_pcg` reuses established data preparation and conventions but
implements its own tested adjoint pair rather than repurposing production
gridding routines:

- **CTF/sigma physics** is evaluated directly via `ctf`/`ctfparams`
  (`simple_ctf.f90`) per pixel, the same physics `image%apply_ctf` uses.
  `image%gen_fplane4rec` is *not* reused: it requires a module-global memoised
  cache (which conflicts with this path's no-module-global rule) and uses a
  packed, padded-box plane convention that clashes with the full-disk native
  representation above.
- **Particle I/O** uses production's batched pattern
  (`prepimgbatch`/`discrete_read_imgbatch`), with `norm_noise` +
  `taper_edges_particle` + `fft` applied per particle before plane extraction —
  the same normalise/taper/transform steps production fuses into
  `norm_noise_taper_edge_pad_fft`.
- **The adjoint is written fresh**, not derived from
  `reconstructor%compress_exp` or `insert_plane_oversamp`, because those are
  production gridding/storage conversions, not established linear adjoints.
