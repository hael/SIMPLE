# `reconstruct3D` PCG backend policy

Contract of the code that runs today: the opt-in CTF- and sigma-weighted
preconditioned-conjugate-gradient (PCG) reconstruction path. Regularization
research (the direct NU-evidence replay precision, Wilson priors) and the
record of removed prior experiments -- including the retired binary-envelope
solvent prior, removed 2026-08-27 -- live in
`doc/implementation_notes/pcg_priors.md`.

The production `reconstruct3D` command accepts the selector
`rec_backend=gridding|pcg`, with `gridding` unchanged as the default. Its `pcg`
branch runs independent kernel accumulation and solves for the even and odd
halfsets, writes the standard halfmaps and FSC, and forms the merged map by
averaging those two dense halfmaps. It supports both shared-memory execution and
distributed worker accumulation followed by a master solve. The former
`reconstruct3D_pcg` command has been retired; there is one production workflow
and one backend selector.

The production backend currently supports full-box fixed-pose
reconstruction from raw or denoised particle sources, multiple populated
states, CTF-aware `cc` or sigma-weighted `euclid`, point-group replication,
standard halfmap/FSC/merged-map names, project output registration, and the
usual final postprocessing handoff. FSC summaries include both the 0.5 and
0.143 resolutions and the cFAR diagnostic, using the same spherical or
density-envelope mask as the corresponding FSC. Each selected particle is read
once into its state/half accumulation and is never reread during PCG iterations.

In distributed execution, each worker atomically publishes one versioned raw
artifact per `(state,half,part)`. It contains the unfolded complex `B`, real
`D`, particle count, geometry, and provenance. The master validates and adds
these artifacts in ascending part-number order. Only after the complete
state/half reduction does it fold the RHS, calculate rho floors, finalize
`Khat`, or solve. Empty partitions publish a valid header-only artifact so the
association order and completeness check do not depend on particle balance.

Box cropping is supported under the constant-field-of-view contract
(`box*smpd == box_crop*smpd_crop`, enforced at entry). It deliberately rejects
`projrec=yes`, fractional/trailing reconstruction, and `conical_fsc=yes`
regularization. Those cases require additional equivalence tests and must not
silently fall back to gridding or matrix-free PCG. (The designed
fractional/trailing algebra — raw `(B,D)` chains blended as
`(u/f) current + (1-u) previous` at full mass, priors applied only after the
blend — is recorded here for when that guard lifts.)

## 1. Production scope and fixed inputs

For each populated state and halfset, the backend solves a CTF- and
noise-weighted fixed-pose least-squares problem with PCG. Orientations, in-plane
shifts, state assignments, half assignments, CTF parameters, phase shifts and
`ctfflag` are read from the project and are never optimized by reconstruction.
The Euclidean objective uses per-particle `sigma2`; correlation is unweighted.
`mskdiam` supplies the soft solve support and `pgrp` is applied by coordinate
replication.

Shared execution accumulates and solves both halfsets in-process. Distributed
workers publish raw accumulators and the master reduces, finalizes and solves
them. Standard project volumes, halfmaps, FSC, cFAR, resolution fields and
postprocessing handoffs are written by the ordinary `reconstruct3D` workflow.
Per-half diagnostics record stop reason, requested/completed iteration count,
convergence flag, initial/final true relative residual, final relative update,
per-iteration residuals and timings.

Particle images are streamed through a bounded batch. The backend calls
`begin_accum`, repeats `accumulate_batch`, closes with `end_accum`, and solves
with `solve_accum`; it never materializes all observed particle planes at once.
The two particle-dependent accumulated quantities are the weighted RHS and the
`|T_i|^2` Gram/sampling-density precursor. The RHS, preconditioner, and kernel
are derived from those quantities, so the PCG iteration itself performs no
particle-image I/O. Peak image/plane memory is therefore constant in particle
count.

Streaming accumulation fuses those two updates. For each particle it evaluates
the CTF, rotated coordinate, and KB interpolation window once, then updates
both accumulators in the same coloured OpenMP traversal. The monolithic RHS
and density routines remain the monolithic test oracle; the streaming-versus-
monolithic gate covers both single-batch and multi-batch fused accumulation.
It gates `B` and the `D`-derived kernel directly at `1e-6` relative error. The
20-step volume comparison has a separate `5e-4` bound because the guarded
reciprocal preconditioner and Krylov recurrence amplify single-precision
accumulation roundoff; passing the volume gate cannot excuse failure of either
direct accumulator gate.

`B` is complex and `D` is real. Finalization computes deterministic per-thread
rho shell statistics, then produces the reciprocal preconditioner and packed
`Khat` together in one OpenMP pass limited to the reachable Fourier sphere.

## 2. Where the code lives

| File | Contents |
| --- | --- |
| `src/main/volume/simple_reconstructor_pcg.f90` | `reconstructor_pcg` operator/solver type |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | shared worker/master PCG orchestration |
| `src/defs/simple_refine3D_fnames.f90` | raw `(state,half,part)` artifact names |
| `src/main/commanders/simple/simple_commanders_rec.f90` | `reconstruct3D` commander and backend default |
| `src/main/ui/simple/simple_ui_refine3D.f90` | `reconstruct3D` UI and backend selector |
| `src/main/strategies/parallelization/simple_rec3D_strategy.f90` | backend dispatch |
| `src/main/params/simple_parameters*.f90` | `pcgop`, `rtol` parameters |
| `src/fileio/simple_sigma2_files.f90` | shared, builder-free sigma2 discovery/loading |
| `src/main/commanders/test/simple_commanders_test_highlevel.f90` | operator/solver tests |

## 3. Statistical model

For volume `x` and 3-D transform `F`: `G_i` extracts an oriented Cartesian
Fourier plane by Kaiser-Bessel (KB) interpolation, `S_i` is the shift phase,
`C_i` the complex CTF (including the phase-flip convention), `N_i` the diagonal
noise covariance from the particle's `sigma2`. With
`K_i = N_i^{-1/2} C_i S_i G_i F`, the solver targets

```text
H x = b,    H = sum_i K_i^dagger K_i + Lambda,    b = sum_i K_i^dagger N_i^{-1/2} y_i
```

`Lambda = lambda*I` is a fixed positive prior, constant for one solve, with
the absolute coefficient `1e-3` (`PCG_LAMBDA`). The relative-lambda CLI
(`pcg_lambda_rel`) was removed as unused; the internal `set_lambda_relative`
mechanism and the deterministic linear fixed-band data scale `s_data(D)`
remain for tests, diagnostics, and future prior anchoring
(`pcg_priors.md`). No current-volume-dependent mask, FSC, filter or
nonlinear clipping appears inside the operator — that linearity is what PCG
relies on.

**The normal operator carries only the real weight `|T_i|^2 = |C_i|^2/sigma2_i`.**
The shift phase is unit-modulus and cancels between the forward and adjoint
passes of `K_i^dagger K_i`. The full complex `T_i = C_i S_i / sqrt(sigma2_i)` is
needed only for `b`, whose adjoint uses `conjg(T_i)` — required for shifted
particles and phase shifts even when the CTF itself is real.

### Output amplitude convention and scale continuity

Solved maps are written at the **data-quotient convention** — the plain
weighted least-squares solution, with no box-size scaling — identical to the
gridding backend after the legacy division/multiplication pair was retired
(`doc/implementation_notes/drop_legacy_box_division.md`). Deapodization is
applied inside the solver; PCG maps must never receive the gridding
correction or a second sampling-density correction. Scale continuity is a
contract: the euclid/sigma2 equilibrium survives any *stable* reference
amplitude but not an abrupt scale change between consecutive refinement
iterations, so any backend handoff, prior, or convention change must
preserve the amplitude scale seen by matching (the historical
gridding-to-PCG handoff crash was exactly such a jump).

### ML two-map contract and warm starts

Refinement solves are two phases from one particle accumulation. The *base*
solve (`H_data + lambda I`, cold start) produces the `_unfil` half pair;
FSC/cFAR and resolution metadata come from that pair. The *ML replay*
(`H_data + P_tau + lambda I`, where `P_tau` is the FSC/SSNR shell-diagonal
precision applied in both the operator and the preconditioner) deterministically
replays kernel finalization from the persisted raw `(B,D)` and produces the
standard maps. The ML replay warm-starts from the previous refinement
iteration's ML half map when one exists on disk — strictly the same half
(gold-standard independence), constant-FOV `read_and_crop` across crop
changes, support re-masked after resampling, the first-iteration noise
`startvol` excluded by name — and otherwise from the base solution. Neither
`P_tau` nor lambda is ever accumulated into raw `B` or `D`.

### Beyond-band diagnostic

`report_beyond_band_excess` (module-level in the strategy, both execution
paths) compares the RMS of shells beyond the matching band with the band-edge
shell and logs `>>> PCG BEYOND-BAND EXCESS` at ratio >= 10. It is the
regression signal for solver defects that park energy above the matched band,
where a later stage transition would expose them to euclid matching. The
structural mitigation is ML-replay convergence (warm start plus adequate
iterations), not spectral smoothing (see the removed-experiment record in
`pcg_priors.md`).

### Backend regression gate

`test=rec3D_backends` reconstructs one fixed particle set with both backends
and hard-fails on gated shell-amplitude, FSC, and radial-flatness criteria
(band capped at `lp`); its ground-truth mode adds map-to-truth FSC and radial
LS-profile flatness. Mutations restoring the legacy box factor or omitting
deapodization must fail it; it is the standing acceptance harness for
convention and deapodization changes.

## 4. Numerical invariants

Load-bearing conventions. Changing any of them requires re-passing §9.

**Full-symmetric-disk planes.** `forward_plane`/`adjoint_plane_add` always use
an unpacked, both-sign-`h` disk, never the packed half. This makes them an exact
adjoint pair for *any* orientation, at 2x plane work. The packed half-plane's
Nyquist-bin bookkeeping is an orientation-dependent silent energy-loss trap that
a single-orientation adjoint test does not catch; the full disk removes that
surface permanently. The KB window is used at its natural odd width
(`2*iwinsz+1`); an even width would not be centred on `nint(loc)` and would
break the mirror consistency the fold depends on.

**Oversampling.** The unknown lives on the native box, but every Fourier
operation runs on an `OSMPL_PAD_FAC`-times padded lattice, centre-padded in and
centre-cropped out, `pad_vol`/`crop_vol` exact adjoints — the same arrangement
gridding uses. Without it KB interpolation is only percent-accurate, the
roll-off envelope swings ~30x across the box, and the Gram kernel cannot match
the operator. Since the padded lattice is already `2*box`, it *is* the grid the
linear (non-wrapping) Gram convolution needs, so kernel and operator share one
lattice.

**Two KB envelopes, both exact discrete transforms.**

| Envelope | Origin | Handling |
| --- | --- | --- |
| Gather | reading the volume through the KB window | `build_env`, divided out by `deapod_mul` |
| Deposition | laying `|T|^2` onto the kernel grid | divided out in the kernel build |

Both use the separable cosine transform of the operator's normalized discrete
KB stencil. This is algebraically identical to scattering a unit spike and
running the 3-D FFT, without allocating another padded volume. It is not
`kbinterpol%instr`: that continuous transform disagrees with the three discrete
renormalized weights by ~2x at the box edge. The matrix-free operator applies
the gather envelope twice, so it is `A = A_env^{-1}(E T E)A_env^{-1}`;
`deapod_mul` brackets it to recover the true `T`. Not optional for real
particles: synthetic data from the operator carries the same envelope and
cancels it (an inverse crime), real particles do not, so an uncorrected solve
returns `E^{-1} x_true`.

**Soft spherical support.** With `mskdiam` set, the solve is constrained by
writing `x = P u` and minimising `||A P u - y||^2`, giving `(P H P) u = P b` —
`P` on both sides of the operator and once on the RHS. `P H P` stays symmetric
positive semidefinite and its null space is never entered from `x = 0`. Besides
removing edge artefacts where deapodization amplifies hardest, it shrinks the
problem: at `mskdiam` 180 in a 256 box the sphere is ~18% of the volume. Edge
profile is production's `cosedge_r2_3d`.

**Scatter reproducibility.** Every scatter into a Fourier accumulator is
parallelized by h-strided colouring, whose separation guarantee holds on
*unwrapped* coordinates only. The accumulator is periodic, so windows that reach
the wrap boundary can collide after folding even in the same colour. All scatter
sites therefore split: the parallel colour sweep skips any window `win_wraps`
flags, and a serial pass handles the rim (pre-filtered by `sq_rim`, a radius
below which wrapping is provably impossible). Serial, not `atomic` — an atomic
would fix the race but leave summation order thread-dependent, and bitwise
run-to-run reproducibility is what makes a residual trace usable as a regression
signal. The rim is a thin outer shell, so the serial cost is small.

## 5. Preconditioner and kernelized operator

**Preconditioner.** `M(k) = rho(k) + floor(shell)`, with
`rho = sum_i G_i^dagger |T_i|^2 G_i` the sampling density. The floor is a fixed
fraction of the shell-mean `rho`, **not** an absolute constant: `rho` spans many
orders of magnitude and is genuinely zero between rotated planes, so an absolute
floor amplifies the least-constrained modes by six to nine orders of magnitude
and PCG fills the map with noise. Modes beyond `padf*Rnat` are unconstrained and
zeroed, making `M` singular but positive semidefinite, which keeps the Krylov
space out of the null space. The envelope belongs in `M` too: the operator being
solved is `E^{-1} T E^{-1}`, so `M^{-1}` brackets its Fourier divide with `E`,
not `E^{-1}`.

**Kernelized (Toeplitz/Gram) operator.** `pcgop=kernel`, the default. Replaces
the per-iteration particle loop with one padded FFT, a pointwise multiply by a
precomputed real `Khat`, and an inverse FFT — per-iteration cost independent of
particle count (~7x faster per iteration). `Khat` uses the standard NUFFT Gram
construction: scatter `|T_i|^2` onto the 2x oversampled grid at doubled
coordinates. The oversampling is the mechanism, not a safety margin; it resolves
the sub-pixel frequency offsets that distinguish the operator from gridding. (A
literal impulse-response kernel does not work: KB weights sum to 1, so it
reduces exactly to gridding and PCG converges in one step to the gridding
answer.) Overall scale is the analytic constant `padsc^2`;
`measure_kernel_scale` is a test-only check that it is right.

The kernelized operator is shift-invariant and the true operator is not, so it
stays an approximation (~3.4% interior error). The matrix-free path is the exact
numerical reference for unit tests and small diagnostic fixtures, **not a
feasible real-workflow backend and not a production fallback**. Its
per-iteration particle loop becomes
prohibitive at experimental particle counts and under symmetry replication.

The approximation is not reliably positive definite under late over-iteration:
on the deterministic fixture the residual bottoms out and then curvature fails
after the useful map has already saturated. Production therefore defaults to
two iterations and rejects `maxits>8`. Low-level tests may use more iterations
when diagnosing the approximation boundary.

`pcgop=kernel` is therefore the only candidate for production workflows. It
must be validated against matrix-free on deterministic, small enough fixtures
where both can run. Agreement of a residual trace alone is insufficient: a
uniform operator-scale error cancels from the CG recurrence while changing map
amplitude. Kernel validation must compare the operator action and scale
directly, then compare fixed-iteration solutions built from the same RHS. Real
workflow acceptance additionally requires independently reconstructed halfmaps
and FSC; a similar operator is not enough.

## 6. Symmetry

Point-group symmetry is applied by **coordinate replication** inside the
operator: each plane pixel is gathered and scattered at all `M` orientations
`R_i . S_g`. Replication is applied to `H` (`accumulate_absT2`,
`apply_normal_matrixfree`) and to `b` (`scatter_plane`, reached from
`apply_adjoint_all` and `accumulate_rhs`), so the system solved is consistently
symmetrized.

**Composition order is production's.** `matmul(R_i, S_g)` in the row-vector
convention `loc = matmul([h,k,0], rot)` is exactly what `sym%apply` produces and
what `reconstructor%insert_fplane` scatters at. `test=pcg_recon` stage 9 asserts
this by composing the reference path through `sym%apply` rather than repeating
the operator's own expression.

`symmats(:,:,1)` is exactly the identity for every supported group, so `pgrp=c1`
(`nsym=1`) is bit-identical to the pre-symmetry path.

Two consequences to keep in view:

- **The `g` loop sits outside the h-strided colour sweep.** Inside it, `M`
  replicas would write `M` footprints per colour and separation would have to
  hold for every pair `(g, g')`, which it does not.
- **Cost scales linearly in `nsym`.** Fine for c2/d2; icosahedral (60) multiplies
  the whole accumulation, including the serial rim pass, by 60. The
  lattice-exact permutation path — available only for groups whose operators are
  signed permutations (`c2`, `c4`, `d2`, `d4`, `t`, `o`) — is the way out and is
  not implemented.

Symmetry-as-constraint and data replication are the same estimator, so there is
no SNR difference between them; symmetry's gain is fewer effective unknowns.

## 7. PCG solver

Standard left-preconditioned CG of `H x = b` (`P H P` under a support
constraint). Zero initialisation. All real-volume dot products use a
deterministic double-precision reduction. The solve fails explicitly on
non-finite or non-positive `dot(p,Hp)`, a non-positive-definite preconditioner,
or a zero RHS, rather than continuing past a broken assumption.

The headline and stopping test are the **true** relative residual
`||r||_2 / ||b||_2`, not the preconditioned M-norm — which is not monotone under
a singular preconditioner and reads like a diverging solver while the solve is
converging. The M-norm is still logged as the diagnostic for the
*preconditioner*: a large gap between the two says `M` models `H` poorly. The
recurrence residual is periodically audited against a recomputed `b - Hx`; they
agree to six significant figures, which licenses reporting the cheap recurrence
norm.

A second criterion, `dx/x <= PCG_XTOL`, exits on diminishing returns: on noisy
real data `||r||/||b||` plateaus above `rtol` while the map keeps settling, so
`dx/x` is what actually terminates a real solve. It is **suppressed by
`rtol <= 0`**, the caller's way of demanding exactly `maxits` iterations —
stage 7 depends on that, since comparing two solves stopped by a data-dependent
criterion asserts nothing.

Successful completion returns a typed `pcg_solver_outcome` identifying `rtol`,
`xtol`, intentional `fixed_iterations`, or exhausted `maxits`, with the
requested/completed iteration count and initial/final residual/update values.
Broken curvature, a non-finite state, a non-positive preconditioner, and a zero
RHS are hard errors. Conversion of those failures into recoverable workflow
outcomes remains future production work.

The final fixed iteration does not apply the preconditioner because no next
search direction will consume it. The terminal residual and update diagnostics
remain unchanged; only the unused M-norm is omitted.

## 8. Sigma2 handling

Gated on `objfun`, matching `reconstruct3D`: `objfun=cc` runs unweighted
(`sigma2 = 1`), `objfun=euclid` requires sigma2 files and weights the fit by
them. A missing sigma2 file under `objfun=euclid` is a hard error, never a silent
unweighted fallback — that would quietly change which objective is minimised.

Sigma2 is per-particle-per-shell, read via `euclid_sigma2` over the group STAR
file (as `flex_analysis` does), then upsampled to the operator's shell range.
Discovery/carry-over/loading lives once in `src/fileio/simple_sigma2_files.f90`,
called by both `flex_analysis` and the `reconstruct3D` PCG strategy. That module is
deliberately builder-free: `builder` depends on `euclid_sigma2`, so the sigma2
side must not depend back on `builder`; callers pass `pftc`, `esig` and
orientations explicitly.

## 9. Tests

`test=pcg_recon` in `simple_commanders_test_highlevel.f90` — one gate, nine
fail-fast stages, in memory with no project I/O. Each stage gates the ones after
it; a visually plausible map is not evidence that the CTF/sigma adjoint is
correct.

| # | Stage |
| --- | --- |
| 1 | adjoint dot-product identity, `T_i = 1` |
| 2 | same with nonzero shift, astigmatic CTF and sigma2 — isolates `build_transfer` |
| 3 | normal-operator symmetry and positive-definiteness |
| 4 | heterogeneous phantom recovery |
| 5 | kernelized-vs-matrix-free operator, scale/energy, and fixed-iteration solution baseline |
| 6 | kernel shift-invariance, CTF-dependence, and the preconditioner |
| 7 | streaming batches and serialized fixed-order raw reduction reproduce monolithic accumulation |
| 8 | deapodization against envelope-free data — the one stage without an inverse crime |
| 9 | symmetry replication equals a c1 build of the symmetry-expanded particle set |

Stages 1-7 generate observations with `forward_plane`, so the gather envelope
cancels; they gate operator *algebra* only. Envelope correctness for real
particles is gated solely by stage 8.

Matrix-free is the oracle for kernel development. Any change to kernel
construction, scale, masks, symmetry, resolution limiting, weights,
preconditioning, or box conversion must add or extend a deterministic
kernel-versus-matrix-free gate. Production-sized runs exercise kernel only;
they do not establish correctness by comparing kernel output with another
kernel output.

Not covered, and worth remembering before trusting a change: the backend's own
particle I/O loop, and replicated symmetry through the matrix-free operator
(stage 9 compares kernel to kernel).

## 10. Current backend exclusions

Not implemented, and hard-errored or absent rather than silently approximated:

- no orientation search or pose optimization inside reconstruction;
- no online partial reconstruction or fractional/trailing updates;
- no adaptive regulariser, automask, nonuniform filtering or FSC feedback inside
  PCG;
- no box cropping, projection-direction compression, conical-FSC regularization,
  or GPU/offload path;
- no reuse of SPIDER BP-CG real-space code — the design is Fourier
  central-section and the architectures do not transfer. There is no licensing
  barrier: SIMPLE is GPL-3.0 and SPIDER GPL-2.0-or-later.

## 11. Reused vs. reimplemented

Established data preparation and conventions are reused; the adjoint pair is
implemented and tested here rather than repurposed from production gridding.

- **CTF/sigma physics** evaluated directly via `ctf`/`ctfparams`
  (`simple_ctf.f90`), the same physics `image%apply_ctf` uses. The per-pixel
  evaluation is the flat `ft_map_ctf_kernel` form with per-particle constants
  hoisted, *inlined rather than called*: the library routine's memoized `(h,k)`
  maps span the `h >= 0` half only, while `lims2` is a full both-sign-`h` disk.
  `image%gen_fplane4rec` is not reused — it requires a module-global memoised
  cache (against this path's no-module-global rule) and a packed, padded-box
  plane convention that clashes with the full-disk native representation.
- **Particle I/O** uses production's batched pattern
  (`prepimgbatch`/`discrete_read_imgbatch`) with `norm_noise` +
  `taper_edges_particle` + `fft` per particle before plane extraction — the same
  steps production fuses into `norm_noise_taper_edge_pad_fft`. Batches stream
  into the accumulators and are discarded, so peak memory is constant in
  `nptcls`.
- **The adjoint is written fresh**, not derived from
  `reconstructor%compress_exp` or `insert_plane_oversamp`: those are production
  gridding/storage conversions, not established linear adjoints.

## 12. Execution-path identity and performance rules

**Shared-memory and distributed execution are two parallelizations of one
algorithm.** Output conventions, warm starts, and diagnostics are implemented
once at module level and invoked identically from both entry points; nothing
may depend on which path produced a map.

Durable performance rules (from the retired production-readiness note):

- one logical particle read per accumulation phase (direct source reads for
  standalone `reconstruct3D`; the validated downscaled cache for cache-enabled
  refinement); particle residency bounded by `MAXIMGBATCHSZ`; no particle
  plane cache — kernel iterations are data-free after `(B,D)`;
- the preconditioner uses the padded Toeplitz lattice; do not trade the
  pad/crop geometry for a native-grid FFT optimization;
- process one state/half at a time; do not overlap the largest accumulation
  and solve scratch allocations; keep `D` real; keep the fused
  reciprocal/`Khat` packing over the reachable Fourier sphere and the exact
  separable discrete KB-envelope construction;
- kernel failure never retries with matrix-free (too slow, hides defects,
  impossible after distributed reduction); fail structurally instead;
- symmetry cost scales with group order (coordinate replication); validate
  C1/C2/D2 and treat high-order groups as a measured-budget contract; the
  lattice-exact permutation path for signed-permutation groups is the
  eventual optimization and is not implemented;
- use the per-phase timing diagnostics, not total wall time, to choose the
  next optimization.
