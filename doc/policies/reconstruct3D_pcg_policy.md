# `reconstruct3D_pcg` policy

Contract of the code that runs today: the experimental CTF- and sigma-weighted
preconditioned-conjugate-gradient (PCG) reconstruction path. Deferred design
work is in section 11.

`reconstruct3D_pcg` is an isolated prototype: it shares no code with the
production gridding path (`reconstructor`, `reconstructor_eo`, `volassemble`,
the online matcher), never writes to the project, and adds no mode switch to any
production program. That boundary is hard — it keeps the comparison against
production honest and keeps a bug here out of production output.

## 1. Scope and fixed inputs

Reads a project and reconstructs one volume for one state by solving a CTF- and
noise-weighted least-squares problem with PCG. For fixed orientations, states
and shifts it answers one question: does solving the weighted normal equations
exactly beat the diagonal gridding approximation by enough to justify the cost?

Every per-particle input is read-only, from the same metadata production uses:

- orientation and in-plane shift are inputs, never optimized;
- one state per run (`state`, default 1);
- particles selected by `sample4rec` (state and, where available, `updatecnt`
  gated), then filtered to the requested state;
- CTF parameters, phase shift and `ctfflag` from the project;
- per-particle Euclidean `sigma2` from the existing estimation machinery (§8);
- `mskdiam`, if given, constrains the solve to a soft spherical support (§4);
- `pgrp` is applied by coordinate replication (§6).

`nparts=1` is enforced with a hard error. No even/odd split and no distributed
execution (§10).

Nothing is written back to the project. Output goes to a new numbered execution
directory: `reconstruct3D_pcg_state<NN>.mrc` and a diagnostics table
`reconstruct3D_pcg_diagnostics.txt` (iteration count, per-iteration relative
residual, timing).

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

For volume `x` and 3-D transform `F`: `G_i` extracts an oriented Cartesian
Fourier plane by Kaiser-Bessel (KB) interpolation, `S_i` is the shift phase,
`C_i` the complex CTF (including the phase-flip convention), `N_i` the diagonal
noise covariance from the particle's `sigma2`. With
`K_i = N_i^{-1/2} C_i S_i G_i F`, the solver targets

```text
H x = b,    H = sum_i K_i^dagger K_i + Lambda,    b = sum_i K_i^dagger N_i^{-1/2} y_i
```

`Lambda = lambda*I` is a fixed positive prior, constant for the whole solve. No
current-volume-dependent mask, FSC, filter or nonlinear clipping appears inside
the operator — that linearity is what PCG relies on.

**The normal operator carries only the real weight `|T_i|^2 = |C_i|^2/sigma2_i`.**
The shift phase is unit-modulus and cancels between the forward and adjoint
passes of `K_i^dagger K_i`. The full complex `T_i = C_i S_i / sqrt(sigma2_i)` is
needed only for `b`, whose adjoint uses `conjg(T_i)` — required for shifted
particles and phase shifts even when the CTF itself is real.

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

**Two KB envelopes, both measured.**

| Envelope | Origin | Handling |
| --- | --- | --- |
| Gather | reading the volume through the KB window | `build_env`, divided out by `deapod_mul` |
| Deposition | laying `|T|^2` onto the kernel grid | measured and divided out in the kernel build |

Both are measured by scattering a unit spike through the operator's own stencil
and transforming back, not derived from `kbinterpol%instr` — the continuous
transform disagrees with the three discrete renormalized weights by ~2x at the
box edge. The matrix-free operator applies the gather envelope twice, so it is
`A = A_env^{-1}(E T E)A_env^{-1}`; `deapod_mul` brackets it to recover the true
`T`. Not optional for real particles: synthetic data from the operator carries
the same envelope and cancels it (an inverse crime), real particles do not, so
an uncorrected solve returns `E^{-1} x_true`.

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
stays an approximation (~3.4% interior error). `pcgop=matrixfree` is the exact
reference and the fallback. The kernel is not equivalent until it gives
equivalent *reconstructions* on real data — a similar operator is not enough.

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
  not implemented (§11).

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

## 8. Sigma2 handling

Gated on `objfun`, matching `reconstruct3D`: `objfun=cc` runs unweighted
(`sigma2 = 1`), `objfun=euclid` requires sigma2 files and weights the fit by
them. A missing sigma2 file under `objfun=euclid` is a hard error, never a silent
unweighted fallback — that would quietly change which objective is minimised.

Sigma2 is per-particle-per-shell, read via `euclid_sigma2` over the group STAR
file (as `flex_analysis` does), then upsampled to the operator's shell range.
Discovery/carry-over/loading lives once in `src/fileio/simple_sigma2_files.f90`,
called by both `flex_analysis` and `reconstruct3D_pcg`. That module is
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
| 5 | kernelized-vs-matrix-free equivalence, all-voxel and interior |
| 6 | kernel shift-invariance, CTF-dependence, and the preconditioner |
| 7 | streaming batch accumulation reproduces the monolithic solve |
| 8 | deapodization against envelope-free data — the one stage without an inverse crime |
| 9 | symmetry replication equals a c1 build of the symmetry-expanded particle set |

Stages 1-7 generate observations with `forward_plane`, so the gather envelope
cancels; they gate operator *algebra* only. Envelope correctness for real
particles is gated solely by stage 8.

Not covered, and worth remembering before trusting a change: the command's own
particle I/O loop, and replicated symmetry through the matrix-free operator
(stage 9 compares kernel to kernel).

## 10. Deliberate non-goals

Not implemented, and hard-errored or absent rather than silently approximated:

- no even/odd half-set split (§11);
- no `nparts>1` / distributed execution;
- no orientation, shift, state, sigma, FSC or resolution write-back;
- no orientation search, probabilistic assignment, online partial
  reconstruction, fractional updates or multi-state reconstruction;
- no adaptive regulariser, automask, nonuniform filtering or FSC feedback inside
  PCG;
- no MPI/distributed reduction, no GPU/offload path;
- no reuse of SPIDER BP-CG real-space code — the design is Fourier
  central-section and the architectures do not transfer. There is no licensing
  barrier: SIMPLE is GPL-3.0 and SPIDER GPL-2.0-or-later.

## 11. Deferred work

**Next: even/odd cross-validation and state assignment.** Independent half-set
solves, so an FSC can be measured and compared against the frozen gridding
reference (map, FSC, local anisotropy, objective trajectory, wall time, peak
memory). The command can already be run twice with explicit even/odd selections;
missing are half-set ownership in the selection/solve path, the FSC comparison,
and multi-state assignment. Half-set independence is a hard requirement: no
cross-half fitted volume may influence a particle's parameters.

**Symmetry, remaining.** The lattice-exact permutation path gated against
replication (§6), and a per-solve asymmetry diagnostic `||x - Pi x|| / ||x||`.

**Longer term.** The operator gives one consistent weighted likelihood for both
reconstruction and continuous parameter refinement (poses, shifts, restricted
CTF/phase). That is a separate command with its own design review — see
`doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md` — and
none of it is in scope for the fixed-pose command.

## 12. Reused vs. reimplemented

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
