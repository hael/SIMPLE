# CTF- and sigma-weighted PCG reconstruction: isolated prototype

## Status

**Implementation in progress.** Milestones 0-3 are built (CTF-free operator;
heterogeneous CTF/shift/sigma; a real-project command; then parallelization,
the section 8 preconditioner and the section 8.1 kernelized operator). See
section 0 for what's implemented, what's deferred, and where the
implementation deviated from this note's original plan — including two places
where the note's own recipe turned out to be wrong (section 4's description of
`transfer_plane`, and section 8.1's impulse-response kernel) and one where its
stated rationale was wrong (section 10's licensing claim).

This note proposes an experimental reconstruction command that solves a
CTF- and noise-weighted least-squares problem with preconditioned conjugate
gradients (PCG).  It is deliberately independent of the production Fourier
gridding/`volassemble` path.  Its first purpose is to establish whether a
more exact iterative data model gives a measurable improvement for a set of
*already aligned* particles.

It is not a replacement for the current reconstruction path and must not
change its numerical result, artifacts, or performance characteristics.

## 0. Implementation status

Code lives in:

- `src/main/volume/simple_pcg_reconstruction.f90` -- the `pcg_reconstruction`
  operator/solver type (milestones 0-1).
- `src/main/commanders/simple/simple_commanders_pcg_recon.f90` -- the
  `reconstruct3D_pcg` production command (milestone 2).
- Tests, all in `src/main/commanders/test/simple_commanders_test_highlevel.f90`:
  - `test=pcg_recon_ctf_free`, `test=pcg_recon_ctf_hetero` -- operator algebra
    (adjointness, symmetry, self-consistent recovery). These generate their own
    observations, so they deliberately run with deapodization off; see "Why the
    synthetic tests could not see the modelling bias" below.
  - `test=pcg_recon_deapod` -- the roll-off correction, against envelope-free
    observations. The only test here that avoids the inverse crime.
  - `test=pcg_recon_kernel` -- the section 8.1 equivalence gate for the
    kernelized operator, plus the section 8 preconditioner.

### Milestone 0 -- CTF-free operator (section 9, stages 2-4)

Implements `forward_plane`/`adjoint_plane_add` (`G_i F` and its adjoint),
`apply_normal`/`apply_adjoint_all`, and a plain (unpreconditioned) CG `solve`.
`T_i = 1` for every particle (no CTF, no shift, no sigma). Validated by
`test=pcg_recon_ctf_free`: adjoint dot-product test (stage 2), normal-operator
symmetry/positive-definiteness test (stage 3), and a no-CTF/no-noise
synthetic-phantom recovery test (stage 4) -- all in-memory, no project I/O.

### Milestone 1 -- heterogeneous CTF/shift/sigma (section 9, stage 5)

Adds `build_transfer` (`T_i = C_i S_i / sqrt(sigma2_i)`), and optional
`use_ctf`/`sig2arr` arguments on `apply_normal`/`apply_adjoint_all`/`solve`
(default off, reproducing milestone 0 exactly). Validated by
`test=pcg_recon_ctf_hetero`: the same three-stage sequence, now with nonzero
shifts, several distinct defocus/astigmatism groups, and a shared
noise-power profile -- this is what actually exercises section 5's release
gate ("must be run with nonzero shifts and a complex transfer"), which
milestone 0's `T_i=1` case could not.

### Milestone 2 -- real-project command (section 9, stage 6, partial)

`reconstruct3D_pcg` reads a project, selects particles via `sample4rec`
filtered to one state, reads per-particle orientation/CTF/shift straight from
the project (read-only), reads particle images in production's own batched
I/O pattern (`prepimgbatch`/`discrete_read_imgbatch`), and solves. Matches
section 2's fixed-input contract: single state, orientations/shifts as
inputs, `nparts=1`, no even/odd, no symmetry, no distributed execution,
writes to a new experimental directory, never writes back to the project,
does not enter through `commander_volassemble`.

**Not yet done, from stage 6's full ask:** FSC comparison against a frozen
gridding reference, local anisotropy indicators, and peak-memory reporting.
Wall time is recorded. This is deliberately deferred, not dropped -- getting
the command running on real data and producing a solver diagnostics table
came first; the comparison step is the next piece of work, once run against
real project data.

### Milestone 3 -- performance, convergence, and correctness of the model

Milestone 2 ran but was unusable: roughly 57 s per CG iteration at 5000
particles and box 128, single-threaded, so hours per reconstruction. Three
causes, all now addressed:

1. **Nothing was cached.** Every CG iteration re-derived each particle's `ori`
   (a deep copy), CTF parameters (string-keyed hash lookups), shift, and
   rotation matrix (`euler2m`, twice per particle), then re-evaluated the CTF
   per Fourier pixel. All of it is constant for the whole solve.
   `prep_particles` now caches it once.
2. **No OpenMP anywhere.** The gather and plane-level loops now parallelize
   directly; the scatter uses the same h-strided colouring scheme as
   `reconstructor%insert_plane_oversamp`, writing safely into one shared
   accumulator rather than per-thread copies.
3. **`solve` was plain CG despite the name** -- `M = I`. The section 8
   preconditioner is now implemented (`build_precond`), which matters most
   exactly where this method should win: with heterogeneous CTFs, `H` is
   severely ill-conditioned at CTF zeros.

Also: gather and scatter are **fused** (KB weights computed once per plane
pixel instead of twice); `apod_mat_3d_fast` replaces `apod_mat_3d`; packed
addressing and Friedel folding are inlined instead of going through
`class(image)` indirect calls 27 times per pixel; and work images persist so
FFTW plans are not rebuilt every iteration.

#### Oversampling -- one cause behind three symptoms

Milestone 2's operator interpolated on the **native** lattice, with no
oversampling. That single omission produced three separate-looking failures:

- the KB roll-off envelope swung ~30x across the box, so the deapodization
  correction was violent (~140x at a corner);
- KB interpolation was only percent-accurate rather than the ~1e-3 that
  `KBALPHA = 2` is designed to deliver;
- the section 8.1 Gram kernel could not be made to agree with the operator.

The operator now pads by `OSMPL_PAD_FAC` for all Fourier work while the
**unknown stays on the native lattice** -- centre-padded going in, centre-
cropped coming out, with `pad_vol`/`crop_vol` as exact adjoints. This is the
same arrangement Fourier gridding uses (`interp_fcomp_oversamp` works at
`loc*PAD_FAC` and rescales by `PAD_FAC**3`). The envelope correction is now
mild (~1.5x per axis) and the kernel agrees with the reference.

Because the padded lattice is already `2*box`, it *is* the grid a linear
(non-wrapping) Gram convolution needs, so the kernel and the operator share
one lattice and no separate padded image exists.

#### Two distinct envelopes

KB interpolation is not neutral: gathering `sum_k w(k-loc) xhat(k)` equals
`FT[what . x]` sampled at `loc`, i.e. it multiplies the volume by an envelope
in real space. Two of them appear, and both had to be handled:

| Envelope | Origin | Handling |
| --- | --- | --- |
| Gather | reading the volume through the KB window | `build_env`, divided out by `deapod_mul` |
| Deposition | laying `abs(T)**2` onto the kernel grid | measured and divided out in `build_kernel` |

Both are **measured**, not derived. `kbinterpol%instr` is the *continuous*
transform of the KB kernel, while this operator uses three discrete
renormalized weight samples; the two disagreed by about a factor of two at the
box edge. Measuring the impulse response through the operator's own stencil
makes every packing, Nyquist and normalization convention match by
construction.

#### Why the synthetic tests could not see the modelling bias

`pcg_recon_ctf_free` and `pcg_recon_ctf_hetero` generate their observations
with `forward_plane`, so the gather envelope appears identically in the data
and the model and cancels -- an inverse crime. That is what makes them a clean
gate on operator *algebra*, and also why they scored 0.9998 correlation while
the forward model was still wrong for real particles, which carry no such
envelope. `pcg_recon_deapod` is the one test that avoids this, manufacturing
envelope-free data via the identity `forward_plane(x/env) = FT[x]`.

#### Status of the kernelized operator

`test=pcg_recon_kernel` passes: 3.4% interior error against the matrix-free
reference (10.1% including boundaries, where shift-invariance is expected to
break down), the kernel is bit-identically invariant to particle shifts, and
it tracks CTF changes. The scale calibration factor lands at `padf**6`, which
is exactly the `padsc**2` the matrix-free operator carries and the kernel does
not -- i.e. the constant is understood rather than absorbed.

3.4% is good enough for a fast approximate operator but is **not** exact, so
the matrix-free path remains the default and the reference; the kernel is
opt-in via `pcgop=kernel`. It should not become the default until it is shown
to give equivalent reconstructions on real data, not merely a similar operator.

**Not started:** stage 7 (thread-count reproducibility), the FSC/anisotropy/
peak-memory parts of stage 6, and everything in section 12 (continuous joint
refinement).

### Where implementation deviated from this note, and why

1. **`image%gen_fplane4rec` (section 4) was not reused.** Tracing its actual
   implementation (`simple_image_ctf.f90`) surfaced two problems: (a) it
   requires a one-time module-global cache (`simple_memoize_ft_maps`), which
   conflicts directly with section 6's "must not use a module global for
   ... FFT plans, or particle planes"; (b) its packed plane convention
   (`k<=0` half-plane over a *padded*, `OSMPL_PAD_FAC`-scaled box) differs
   from the full-symmetric-disk, native-box representation this
   implementation uses (point 3 below), which would have reintroduced the
   packed-half-plane indexing risk section 5 explicitly says to avoid.
   `build_transfer` instead evaluates `C_i`/`S_i` directly via `ctf`/
   `ctfparams` (`src/main/ctf/simple_ctf.f90`) -- the same underlying CTF
   physics and reciprocal-pixel convention `image%apply_ctf` uses, just
   invoked per-pixel rather than through the `fplane_type` wrapper.

2. **A factual correction to this note's section 4.** It describes
   `fplane_type%transfer_plane` (observation mode) as containing
   `CTF*shift/sqrt(sigma2)`. Tracing the actual code shows the shift is *not*
   in `transfer_plane` -- it's folded into `cmplx_plane` instead
   (`cmplx_plane = S_i*y_i/sqrt(sigma2_i)`, `transfer_plane = C_i/sqrt(sigma2_i)`
   only). This doesn't affect the implementation (which never called
   `gen_fplane4rec`, per point 1), but is worth fixing here since a future
   reader relying on this note's description of `transfer_plane` would be
   misled.

3. **The full-symmetric-disk plane representation (section 5's "may use ...
   full complex work planes" escape hatch) became the permanent design, not
   just a validation aid.** `forward_plane`/`adjoint_plane_add` always use an
   unpacked, both-sign-`h` disk, never the packed half-plane. This was
   load-bearing, not cosmetic: getting the packed convention's Nyquist-bin
   bookkeeping right (which slot the redundant Friedel mate's energy
   actually lands in after `adjoint_plane_add`'s scatter) was the hardest
   bug in milestone 0 -- an orientation-dependent, silent energy-loss bug
   that a single-orientation adjoint test did not catch; it only surfaced
   once tested across multiple, sufficiently different rotations. Keeping
   the full-disk representation permanently avoids reintroducing that
   surface for a performance gain that hasn't been shown to matter yet.

4. **`forward_plane`, `adjoint_plane_add`, `apply_normal`, `apply_adjoint_all`,
   `build_transfer`, and `extract_native_plane` are public**, not the private
   helpers section 6's sketch shows. The test commanders (and, for
   `extract_native_plane`, the `reconstruct3D_pcg` commander) need to drive
   the operator's individual pieces directly to validate the adjoint
   identity *before* trusting `solve()` -- consistent with section 9's own
   "do not proceed past a failed earlier stage" discipline, just implemented
   as public API rather than a private-module self-test.

5. **Sigma2 is per-particle, read via `euclid_sigma2` as flex_analysis does.**
   Milestone 2 had used the lightweight `sigma2_binfile` and averaged the
   selection into one shared profile, on the assumption that the full
   `euclid_sigma2` path required standing up search/polar-FT machinery. That
   assumption was wrong: `euclid_sigma2%new` takes a `polarft_calc` only to
   perform a pointer assignment (`polarft_core%assign_sigma2_noise`), so an
   unconstructed `build%pftc` is fine and no strategy3D toolbox is needed.
   Milestone 3 therefore carries sigma files over from prior runs, picks the
   highest-iteration group star file, forces global sigma dispatch (reused
   star files are typically single-group), and reads genuine
   per-particle-per-shell spectra.

   That logic was initially duplicated from flex_analysis. It now lives once,
   in `src/fileio/simple_sigma2_files.f90`, which both flex_analysis and
   `reconstruct3D_pcg` call. The module is deliberately builder-free --
   `builder` depends on `euclid_sigma2`, so anything on the sigma2 side must
   not depend back on `builder`; callers pass `pftc`, `esig` and the
   orientations explicitly.

6. **`objfun`-gated sigma2, not file-presence-gated.** Not specified in this
   note at all -- added to match `reconstruct3D`'s own convention (`objfun=cc`
   needs no sigmas, `objfun=euclid` requires them) per direct instruction,
   rather than the note's original silent-fallback framing. `reconstruct3D_pcg`
   now requires (`THROW_HARD`s if absent) a sigma2 file when `objfun=euclid`,
   instead of quietly reconstructing unweighted if the file happens to be
   missing.

7. **The KB window is used at its natural odd width** (3 points for
   `KBWINSZ=1.5`), not blindly assumed to be even. This isn't a deviation
   from anything the note specifies, but is worth recording: an even-width
   window would not be centered on `nint(loc)`, breaking the
   mirror-consistency the full-disk adjoint construction (point 3) depends
   on. `KBWINSZ=1.5` already gives an odd window in this codebase, so no
   change to KB parameters was needed -- `simple_pcg_reconstruction.f90`
   just computes the width as `2*iwinsz+1` explicitly, to document the
   invariant instead of trusting it silently.

8. **Section 8.1's impulse-response kernel does not work, and was replaced.**
   The note says to build the Gram kernel as `h = H_data(delta_at_origin)`.
   That reduces *exactly* to gridding, and therefore gains nothing at all.
   The reason: `kbinterpol%apod_mat_3d` normalizes the KB weights to sum to 1,
   so `G` applied to a constant Fourier volume returns that constant, which
   makes the row sums of `M = sum_i G_i^dagger |T_i|^2 G_i` precisely `rho`,
   the gridding sampling density. A same-size kernel built that way gives
   `H = F^dagger diag(rho) F + Lambda`, whose exact solution is the gridding
   reconstruction -- PCG would converge in one step to the answer it was
   supposed to improve on.

   The implementation instead uses the standard NUFFT Gram construction:
   scatter the weights `|T_i(p)|^2` onto a **2x oversampled** Fourier grid at
   **doubled** coordinates. The oversampling is not a safety margin, it is the
   entire mechanism -- it is what resolves sub-pixel frequency offsets that the
   native grid cannot represent, and hence what makes the operator differ from
   gridding at all. Since `t = IFFT_2N(Khat)`, the real-space kernel is never
   formed explicitly; the scattered array *is* the Fourier multiplier.

   This remains an approximation (it is shift-invariant, and the true operator
   is not), which is why it is gated on the equivalence test and is not the
   default.

9. **Section 10's licensing rationale is wrong.** It forbids SPIDER-derived
   code "because SPIDER is GPL-licensed". SIMPLE is **GPL-3.0** (`LICENSE`)
   and SPIDER is **GPL-2.0-or-later** (see the header of
   `SPIDER/src/bpcg.f`), and "or later" means SPIDER's code may legitimately
   be used under GPL-3.0. There is no licensing barrier.

   The *practical* conclusion is unchanged, for a different reason: SPIDER's
   BP CG is real-space, ray-driven, and restricted to a sphere via run-length
   voxel spans, whereas this design is Fourier central-section. The
   architectures do not transfer. Two of its ideas remain worth revisiting --
   merging duplicate viewing directions with integer multiplicity weights (so
   per-iteration cost scales with distinct orientations rather than particle
   count), and its derivative-stencil regularizers with roughly 15-20
   iterations -- but neither requires copying code.

Everything else in section 10's non-goals list still holds: no changes to
`reconstructor`/`reconstructor_eo`/`volassemble`/online matcher, no new
production `refine3D` argument, no MPI/GPU, and no adaptive regularizer/
automask/nonuniform filtering/FSC feedback inside PCG.

## 1. Question being tested

The current Cartesian reconstruction is efficient because it approximates the
normal equations as diagonal on the gridded 3-D Fourier lattice.  Per-particle
data contribute a numerator proportional to `CTF* y / sigma2` and a density
proportional to `CTF**2 / sigma2`, after which the volume is obtained by a
pointwise division plus regularization.

The proposed prototype instead applies the complete forward model and its
adjoint at every iteration.  The scientific question is therefore narrow:

> For fixed orientations, states, and shifts, does solving the weighted normal
> equations improve a controlled reconstruction enough to justify its cost?

This is an operator-validation and reconstruction experiment.  It does **not**
include orientation search, probabilistic assignment, online partial
reconstruction, fractional updates, multi-state reconstruction, symmetry,
FSC-driven prior refresh, automasking, or nonuniform postprocessing.

## 2. Fixed-input prototype contract

The initial command should have a new, unambiguous name, for example
`reconstruct3D_pcg`.  It reads a project and reconstructs exactly one volume
from the selected particles using their current project-carried assignments:

- orientation and in-plane shift are inputs, not optimized;
- only one state is selected;
- particles are selected by the existing project/segment mechanisms;
- the CTF, phase shift, `ctfflag`, and per-particle Euclidean `sigma2` are
  taken from the same metadata/estimation machinery used by reconstruction;
- no project orientation, shift, state, sigma, FSC, or resolution metadata is
  written back;
- output is a new experimental execution directory containing the map,
  optional density/preconditioner map, and a solver diagnostics table.

For the first milestone, require `nparts=1`, no even/odd split, no symmetry,
and no distributed execution.  The command may be run twice with explicit
even/odd particle selections later, but even/odd ownership is not part of the
first numerical kernel.

The command belongs in a dedicated experimental commander/strategy route.  It
must not enter through `commander_volassemble`, reuse its output filenames, or
add a mode switch to the production reconstructor.  This preserves an honest
comparison with the present Fourier-gridding implementation.

## 3. Statistical model

Let `x` be a real-space volume and let `F` be the 3-D Fourier transform.  For
particle `i`, define:

- `G_i`: extract an oriented Cartesian Fourier plane from `F x` by the chosen
  Kaiser--Bessel interpolation;
- `S_i`: complex phase factor for the project-carried 2-D shift;
- `C_i`: complex CTF transfer, including the configured phase-flip convention;
- `N_i`: diagonal noise covariance whose radial entries are the particle's
  `sigma2` spectrum;
- `y_i`: the normalized, padded particle Fourier plane before CTF correction.

The whitened observation operator is

\[
K_i = N_i^{-1/2} C_i S_i G_i F.
\]

The prototype solves the MAP normal equations

\[
H x = b,
\qquad
H = \sum_i K_i^\dagger K_i + \Lambda,
\qquad
b = \sum_i K_i^\dagger N_i^{-1/2}y_i.
\]

`Lambda` is a fixed, positive semidefinite prior for the whole solve.  The
first version should use the same radial `1/tau2` concept as the production
reconstructor, calculated once before PCG begins.  A conical/directional prior
is a later extension, provided it too remains fixed during a solve.

No current-volume-dependent mask, FSC estimate, filter, tau update, or
nonlinear clipping may appear inside `apply_normal`.  Any such operation would
break the linear Hermitian positive-definite system on which PCG relies.

### 3.1 Per-particle operations

For one trial direction `p`, the data part of `H p` is:

```text
v       = FFT3(p)
q_i     = G_i(v)
m_i     = T_i * q_i                       ! T_i = C_i S_i / sqrt(sigma2_i)
a_i     = conjg(T_i) * m_i
z       = FFT3_inverse(sum_i G_i_dagger(a_i))
H p     = z + Lambda(p)
```

The right-hand side is the same adjoint path with
`m_i = y_i / sqrt(sigma2_i)`:

```text
b = FFT3_inverse(sum_i G_i_dagger(conjg(T_i) * y_i / sqrt(sigma2_i)))
```

The adjoint must use `conjg(T_i)`, not `T_i`, even though many conventional
CTFs are real.  This is necessary for shifted particles, phase shifts, and
future complex transfer functions.

## 4. Existing SIMPLE contracts to reuse

The prototype should reuse established data preparation rather than recreate
CTF or sigma rules.

`image%gen_fplane4rec(..., observation_model=.true.)` is the intended source
of particle planes.  In that mode it creates a whitened observation
`y/sqrt(sigma2)` and separately stores the forward complex transfer
`CTF*shift/sqrt(sigma2)` in `fplane_type%transfer_plane`.  This is the correct
representation for the PCG model.  It also handles the current CTF flag and
phase-flip convention.

This differs intentionally from the ordinary reconstruction mode, where
`fplane_type%cmplx_plane` already contains `conjg(CTF*shift)*y/sigma2` and
`ctfsq_plane` contains `abs(CTF*shift)**2/sigma2`.  Do not use a conventionally
preweighted plane as a PCG observation and apply the transfer again: doing so
would double-count the CTF and sigma factors.

The following current components are therefore dependencies, not templates to
copy:

| Need | Existing owner | Prototype use |
| --- | --- | --- |
| Particle normalization, padding, FFT, CTF/sigma plane generation | `image%gen_fplane4rec` | Generate observation-model planes once per batch. |
| Particle CTF and assigned shift | `sp_project` / current reconstruction preparation | Obtain exactly the current assignment metadata. |
| Fourier-plane geometry and KB interpolation conventions | `fplane_type`, `kbinterpol`, `reconstructor` | Reuse coordinate/index conventions; establish a new tested adjoint pair. |
| `sigma2` ownership | `builder%esig` | Read the existing fixed spectra; do not estimate sigma in this command. |
| Fixed radial prior | current `reconstructor` prior logic | Reuse the scientific definition, but expose it as a fixed linear operator. |

## 5. Numerical representation and adjoint requirement

The solver unknown is a real volume.  Each normal-operator application may
temporarily use Fourier storage, but both the forward and adjoint must refer to
the same canonical grid and the same FFT normalization.

This deliberately avoids calling `reconstructor%compress_exp` as the adjoint:
that routine is a production gridding/storage conversion, not an established
linear adjoint.  Instead, introduce a private experimental operator with two
matched methods:

```fortran
call op%forward(volume_real, particle_index, model_plane)
call op%adjoint_add(residual_plane, particle_index, fourier_accumulator)
```

`forward` is `G_i F` and `adjoint_add` is `F^dagger G_i^dagger` after the
caller has applied `conjg(transfer_plane)`.  They must use identical:

- Fourier-coordinate and Friedel-symmetry conventions;
- padded/native sampling conversion;
- KB stencil support, apodization, and normalization;
- particle plane inclusion mask and Nyquist rule;
- treatment of the `k = 0` plane and self-conjugate frequencies.

The first implementation may use explicit full complex work planes and a full
complex 3-D accumulation buffer if that makes the inner product unambiguous.
That costs memory but removes packed-half-plane multiplicity ambiguity during
validation.  Compact storage and fused kernels are optimization work only
after the adjoint test passes.

The required complex inner-product test is

\[
\frac{|\langle G_i F x, q\rangle - \langle x, F^\dagger G_i^\dagger q\rangle|}
     {\max(1,|\langle G_i F x, q\rangle|,|\langle x,F^\dagger G_i^\dagger q\rangle|)}
\leq 10^{-5}
\]

in single precision on small deterministic fixtures.  The test must be run
with nonzero shifts and a complex transfer.  It is a release gate for all
later PCG tests.

## 6. Proposed isolated module boundary

Add a new numerical module in `src/main/volume/`, tentatively
`simple_pcg_reconstruction.f90`.  It should be private by default and export
only a focused solver type, for example `pcg_reconstruction`.

Suggested responsibilities:

```text
pcg_reconstruction
  new(params, project, selection, options)
  prepare_observations(images_or_stack_batch)
  build_rhs_and_preconditioner()
  solve(result_volume, diagnostics)
  kill()

private helpers
  apply_normal(p, hp)
  apply_prior(p, prior_p)
  apply_preconditioner(r, z)
  forward_plane(volume, orientation, plane)
  adjoint_plane_add(plane, orientation, fourier_sum)
  dot_real_volume(a, b)
```

The solver owns only numerical work buffers and explicit immutable references
to the parsed `parameters`, project/builder context, and selected particles.
It must not use a module global for `parameters`, `builder`, sigma, FFT plans,
or particle planes.  The commander owns file and execution-directory policy;
the solver owns no user-facing filenames.

For the first implementation, particle processing should be batched.  Prepare
each batch once, use it for the forward and adjoint passes of that operator
application, then release it.  Reopen neither stack nor per-particle files in
an inner PCG loop.  A later cache policy can trade RAM for iteration time.

## 7. PCG algorithm and stopping rules

Use standard left-preconditioned CG for the fixed system.  With `M` the
positive diagonal approximation to `H`:

```text
x = 0                                 or a documented supplied initial map
r = b - H x
z = M^-1 r
p = z
rho = dot(r, z)

for iter = 1:maxit
    hp    = H p
    alpha = rho / dot(p, hp)
    x     = x + alpha * p
    r     = r - alpha * hp
    z     = M^-1 r
    rho_new = dot(r, z)
    p     = z + (rho_new / rho) * p
    rho   = rho_new
end for
```

All real-volume dots use a deterministic double-precision reduction.  Fail
explicitly on non-finite values, non-positive `dot(p,hp)`, or a zero/invalid
preconditioner denominator.  Do not silently continue after a broken
Hermitian/positive-definite assumption.

Initial stopping criteria:

- relative preconditioned residual `sqrt(rho/rho0) <= rtol`;
- maximum iterations `maxit`;
- optional stagnation guard based on no meaningful residual decrease for a
  fixed number of iterations.

Record, per iteration: objective data term, fixed-prior term, total objective,
relative residual, `alpha`, `beta`, `dot(p,hp)`, and wall time.  The objective
must decrease to rounding tolerance for the deterministic test cases.

## 8. Preconditioner

The initial preconditioner should be the present gridding approximation:

\[
M(\mathbf{k}) = \rho(\mathbf{k}) + \lambda(\mathbf{k}),
\]

where `rho` is independently accumulated from
`abs(transfer_plane)**2` through the same KB adjoint stencil and `lambda` is
the fixed prior on the canonical Fourier lattice.  Applying `M^-1` consists
of FFT, elementwise guarded division, and inverse FFT.

This is not merely convenient: it tests whether PCG can correct the
off-diagonal interpolation/support effects while retaining the speed benefits
of the gridding solution.  The direct gridded reconstruction can also be used
as an optional initial `x`, but zero initialization is the required baseline.

### 8.1 Phase-2 option: kernelized normal operator

After the matrix-free operator has passed every phase-1 validation gate, an
optional acceleration may replace its repeated particle passes.  The key
observation is that, for fixed orientations and fixed particle transfer
functions, the data normal operator is approximately shift invariant in the
interior of the reconstruction box:

\[
H_{\mathrm{data}} = \sum_i K_i^\dagger K_i.
\]

The project-carried shift is a unit-modulus Fourier phase in `T_i`.  It cancels
between the forward and adjoint paths in `K_i^dagger K_i`; the corresponding
normal weight is therefore

\[
|T_i|^2 = |C_i|^2 / \sigma_i^2.
\]

For fixed orientations, CTFs, plane limits, and sigma spectra, construct a
3-D point-spread/normal kernel once and apply it as a **linear**, zero-padded
3-D convolution:

```text
build once:
  h_data = matrix_free_data_normal(delta_at_volume_origin)
  Hhat   = FFT3(padded(h_data))

each PCG iteration:
  hp_data = crop(IFFT3(FFT3(padded(p)) * Hhat))
  hp      = hp_data + Lambda(p)
```

**Correction (milestone 3): the impulse-response construction above does not
work.** It reduces exactly to gridding. The KB weights are normalized to sum
to 1, so `G` applied to a constant Fourier volume returns that constant, which
makes the row sums of `M` precisely the gridding density `rho`; the resulting
operator is `F^dagger diag(rho) F + Lambda`, whose exact solution *is* the
gridding reconstruction. What is implemented instead is the standard NUFFT
Gram construction: scatter `abs(T_i)**2` onto a **2x oversampled** Fourier
grid at **doubled** coordinates. The oversampling is the mechanism, not a
margin -- it resolves the sub-pixel frequency offsets the native grid cannot
represent, which is precisely what distinguishes the operator from gridding.
The scattered array is itself the Fourier multiplier, so the real-space kernel
`h_data` is never formed. See section 0, deviation 8.

This is a normal-operator acceleration, not a new reconstruction model.  It
does not use the current `rho` division as its operator, and it does not
replace the fixed prior with a data-density heuristic.  The matrix-free
operator remains the reference implementation.

#### Boundary and validity limits

The convolution property is not automatic for a finite reconstruction:

- finite Fourier-plane masks and frequency cutoffs;
- finite Kaiser--Bessel support and interpolation truncation;
- real-space reconstruction support and output cropping; and
- any spatially varying particle-selection rule

can make the normal operator non-Toeplitz, especially near boundaries.  The
kernelized path must use enough zero padding for **linear**, never circular,
convolution.  It must not be enabled if the compared matrix-free operator has
an error above the agreed tolerance.

Even if the interior is convolutional, retain PCG.  Output cropping and a
finite-support volume make the actual system Toeplitz rather than simply
diagonalizable by a same-size FFT; the padded convolution only makes each
matrix-vector product inexpensive.

#### Required phase-2 equivalence tests

For deterministic small fixtures, compare `H_data p` from the kernelized and
matrix-free paths for random real vectors `p`:

\[
\frac{\|H_{\mathrm{kernel}}p-H_{\mathrm{matrix-free}}p\|_2}
     {\max(1,\|H_{\mathrm{matrix-free}}p\|_2)} \leq \epsilon_{\mathrm{kernel}}.
\]

Report this error separately for all voxels and for an interior region whose
margin is at least the interpolation-kernel support.  Set
`epsilon_kernel` from measured single-precision roundoff in the phase-1
operator test; do not choose it after observing a desired reconstruction.

The phase-2 test matrix must include heterogeneous, astigmatic CTFs and
per-particle sigma spectra.  It must additionally verify that changing only
particle shifts leaves the normal kernel unchanged while changing a CTF or
sigma spectrum rebuilds it.  Compare full PCG convergence, objective, output
map, wall time, and peak memory against the matrix-free reference before
claiming an acceleration.

## 9. Validation sequence

Implement and pass these in order.  Do not proceed past a failed earlier
stage.

1. **Transfer convention fixture.** For several Fourier pixels, verify the
   observation-model plane stores `y/sqrt(sigma2)` and the transfer stores
   `CTF*shift/sqrt(sigma2)`, including nonzero shift, CTF phase shift, and
   phase-flip mode.
2. **Adjoint-dot-product test.** Validate `forward`/`adjoint_add` as specified
   above using random seeded small volumes/planes.
3. **Normal-operator test.** Check `dot(p,Hq) == dot(Hp,q)` to tolerance and
   `dot(p,Hp) > 0` for nonzero `p` with a positive prior.
4. **No-CTF/noise synthetic reconstruction.** Generate noiseless projections
   from a known small volume.  PCG must recover the expected solution and
   converge monotonically.
5. **Heterogeneous CTF and sigma synthetic reconstruction.** Add distinct
   CTFs, shifts, and noise spectra.  Compare recovered-map error and weighted
   objective against direct gridding under the identical fixed inputs.
6. **Real aligned-particle smoke test.** Compare the gridding and PCG output
   on a small, frozen selection: visual map, FSC to a frozen reference, local
   anisotropy indicators, objective trajectory, wall time, and peak memory.
7. **Reproducibility.** Run the same fixture with one and several OpenMP
   threads.  Permit documented floating-point reduction variation only; require
   the same convergence decision and scientifically indistinguishable map.
8. **Kernelized normal operator (phase 2 only).** Run the equivalence tests in
   section 8.1 before allowing the FFT-convolution operator in a solver run.
   Compare its PCG result and diagnostics with the matrix-free reference, then
   quantify the speed and memory trade-off.

The experiment must report a failure if the operator test fails.  A visually
plausible reconstruction is not evidence that the CTF/sigma adjoint is correct.

## 10. Deliberate non-goals for phase 1

- no modification of `reconstructor`, `reconstructor_eo`, `volassemble`, or
  online matcher partial reconstruction;
- no new production `refine3D` argument or mode switch;
- no MPI/distributed reduction;
- no GPU/offload path;
- no adaptive regularizer, automask, nonuniform filtering, or FSC feedback
  inside PCG;
- no attempt to reproduce SPIDER BP-CG real-space interpolation.

**Correction (milestone 3).** This section previously also forbade "code
copied or translated from SPIDER" on the grounds that SPIDER is GPL-licensed.
That rationale is wrong: SIMPLE is GPL-3.0 and SPIDER is GPL-2.0-*or-later*,
so SPIDER code could legitimately be used here. The design is nonetheless an
independent implementation, because SPIDER's BP CG is real-space and
ray-driven while this is Fourier central-section -- the architectures do not
transfer. See section 0, deviation 9.

## 11. Acceptance decision after phase 1

Keep the prototype separate unless all of the following are true:

1. adjoint and normal-operator tests pass robustly;
2. CTF/sigma handling agrees with the frozen SIMPLE convention;
3. the solver gives a reproducible benefit on at least one scientifically
   meaningful difficult case, not only a lower algebraic residual;
4. its cost is quantified relative to direct gridding and is acceptable for a
   clearly defined use case.

If those conditions hold, a separate design note can address even/odd support,
distributed batches, a stable command-line contract, and whether any mature
operator primitive should be shared with production reconstruction.  The
kernelized normal operator in section 8.1 remains optional until it has passed
its own equivalence tests; it is not part of the phase-1 acceptance decision.

## 12. Forward evolution: continuous joint refinement

The longer-term value of the validated operator is that it provides one
consistent weighted likelihood for both reconstruction and continuous
parameter refinement.  With volume `x`, per-particle pose parameters
`theta_i`, and CTF/phase parameters `eta_i`, the target is

\[
\min_{x,\theta,\eta}
\sum_i \left\|
N_i^{-1/2}\left[y_i-C_i(\eta_i)S_i(\theta_i)P_i(\theta_i)x\right]
\right\|^2
+ R_x(x) + R_\theta(\theta) + R_\eta(\eta).
\]

For fixed poses and microscope parameters, the PCG system specified in this
note is the conditional MAP volume solve.  For a fixed volume, the same
forward model supplies the residual and derivatives needed to update poses,
shifts, phase, or selected CTF parameters.  This avoids optimizing alignment
under one interpolation/CTF/noise convention and reconstructing a volume under
another.

### 12.1 Intended optimization shape

Do not attempt fully simultaneous unconstrained optimization first.  The
future route should be block-coordinate or variable-projection refinement:

```text
repeat until outer convergence:
  1. hold theta and eta fixed; solve x with the validated PCG volume solver
  2. hold x fixed; update a constrained block of theta and/or eta
  3. accept/reject or damp the block update using the same weighted objective
```

The existing discrete/probabilistic search machinery remains valuable as an
initializer and candidate generator.  Continuous updates are a local
refinement layer after a credible discrete assignment, not a replacement for
global pose search.

### 12.2 Staged parameter expansion

Add one parameter family at a time, with its own synthetic derivative and
recovery tests:

1. **Continuous 2-D shifts.** These are the first extension because their
   Fourier derivatives are analytic.  If `S_i(f)=exp(-2*pi*i*f dot t_i)`,
   then `dS_i/dt_x=-2*pi*i*f_x*S_i` and likewise for `t_y`.
2. **Local rotations.** Parameterize a small orientation increment in the
   tangent space of `SO(3)`, not three unconstrained Euler-angle increments.
   Differentiate the oriented Fourier-coordinate interpolation consistently
   with the forward operator.
3. **Restricted CTF/phase refinement.** Begin with strongly regularized global
   or optics-group parameters, such as phase shift or a defocus offset.  Only
   then consider per-particle defocus/astigmatism corrections, with explicit
   priors and bounds.

For each block, use a trust-region Gauss--Newton or damped L-BFGS update, not
an undamped Newton step.  The implementation must expose the residual,
Jacobian-vector product, and adjoint-Jacobian-vector product as tested
operator methods; it must not form a dense particle-parameter Hessian.

### 12.3 Identifiability and validation safeguards

Joint refinement has genuine degeneracies: alignment error can resemble volume
blur, CTF error can resemble sharpening or a changed volume prior, and an
over-flexible per-particle CTF model can fit noise.  The future design must
therefore include:

- fixed or tightly controlled volume priors during each outer step;
- parameter-specific bounds, Gaussian priors, and resolution schedules;
- trust-region acceptance based on an independently recomputed objective;
- synthetic recovery tests with known shifts, rotations, phase shifts, and
  defocus perturbations;
- independent even/odd half-set volume solves and parameter updates, with no
  cross-half fitted volume used to improve a particle's parameters; and
- held-out or cross-validation diagnostics before accepting additional
  parameter freedom.

The initial `reconstruct3D_pcg` command remains intentionally fixed-pose.
Continuous joint refinement is a separate future command/strategy and needs a
new design review once the phase-1 operator, PCG solver, and half-set contracts
are established.

### 12.4 Later option: stochastic-gradient refinement

The same likelihood also permits minibatch stochastic optimization.  For a
batch `B` sampled from the fixed particle selection, an unbiased volume-gradient
estimator is

\[
g_B(x)=\frac{N}{|B|}\sum_{i\in B}K_i^\dagger(K_i x-y_i^\mathrm{white})
+\Lambda x,
\]

when particles are sampled uniformly and the CTF/sigma weighting remains
inside each `K_i`.  A preconditioned update is

\[
x \leftarrow x-\gamma M^{-1}g_B(x),
\]

where `M` is the fixed gridding-density-plus-prior preconditioner from section
8.  The same batch residual provides gradients for the active pose or CTF
parameter block.

This is attractive for very large data sets and eventual continuous joint
refinement because each batch follows the existing single-read particle I/O
pattern.  It is not a phase-1 replacement for deterministic PCG.  PCG remains
the reference volume solver and the periodic global correction/re-anchoring
step against which a stochastic path is judged.

An initial stochastic design should use stratified minibatches that represent
the two half sets, orientation coverage, and optics/CTF groups.  It should use
the interpretable Fourier preconditioner above with conservative momentum or a
stochastic Gauss--Newton variant, rather than beginning with an unconstrained
adaptive optimizer.  CTF zeros and sparsely represented orientations otherwise
produce high-variance high-frequency updates.

Required safeguards are:

- retain a fixed prior within each scheduled stochastic epoch;
- periodically evaluate the exact full-data objective and gradient norm;
- periodically run deterministic PCG or a full-data normal-equation check;
- schedule step sizes and any momentum explicitly, with rollback on an
  independently evaluated objective increase;
- maintain completely independent even/odd batches and fitted volumes; and
- compare against deterministic PCG on every synthetic fixture before using
  stochastic updates for continuous pose or CTF refinement.

Stochastic refinement is therefore a phase-3 research option: useful once the
operator is proven and the cost of full global passes becomes the limiting
factor, but not a substitute for establishing the deterministic problem first.
