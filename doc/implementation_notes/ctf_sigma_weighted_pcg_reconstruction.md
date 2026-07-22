# CTF- and sigma-weighted PCG reconstruction: isolated prototype

## Status

**Proposal for review. No implementation is implied by this note.**

This note proposes an experimental reconstruction command that solves a
CTF- and noise-weighted least-squares problem with preconditioned conjugate
gradients (PCG).  It is deliberately independent of the production Fourier
gridding/`volassemble` path.  Its first purpose is to establish whether a
more exact iterative data model gives a measurable improvement for a set of
*already aligned* particles.

It is not a replacement for the current reconstruction path and must not
change its numerical result, artifacts, or performance characteristics.

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

The initial, safest kernel construction is the impulse response shown above.
It reuses the already verified matrix-free `K^dagger K` path and cannot change
the CTF/sigma convention.  A later analytic construction may sum the rotated
Kaiser--Bessel projection autocorrelation weighted by `abs(T_i)**2`, but only
after it agrees with the impulse kernel and the matrix-free operator.

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
- no attempt to reproduce SPIDER BP-CG real-space interpolation;
- no code copied or translated from SPIDER.

The last point is intentional: SPIDER is GPL-licensed.  This design is based
on the standard weighted least-squares/PCG formulation and SIMPLE's existing
data-model contracts, with a new implementation and independent tests.

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
