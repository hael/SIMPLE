# Ghostbuster's forward and adjoint operators in SIMPLE

Status: plan, 2026-10-23, rulings of the same day applied (section "Rulings").
Nothing implemented. Sources read: specter (the package that contains
Ghostbuster) at `5786651` in `~/src/specter`, SIMPLE at `630c2604a`. specter
was read, not run. The work is carried out by an unattended run on the Dell
(section 7).

## Summary

Ghostbuster predicts a cryo-EM image from a 3D electrostatic potential with
wave optics: the beam picks up a phase in the specimen, travels through the
objective lens, and the detector records the squared modulus of the wave. Its
reconstruction runs that prediction backwards with automatic differentiation.
This plan brings both directions into SIMPLE as Fortran, first as a standalone
program, `reconstruct3D_wave`, that reconstructs half maps from a project with
known poses and CTF. Refinement workflows come later and build on it: Ewald
sphere correction, CTF refinement and pose refinement all need this forward
operator and its derivative.

The plan is not a translation. Three facts make a SIMPLE-native design both
shorter and faster than the Python original.

1. For the single-scattering models (projection, first Born, Rytov) the exit
   wave depends on the volume only through one complex 2D spectrum: the
   volume's Fourier transform sampled on a paraboloid, the Ewald surface.
   specter samples it with a non-uniform FFT library. SIMPLE's PCG
   reconstructor already samples the volume's transform at arbitrary 3D
   points with a Kaiser-Bessel window and has the exact transpose of that
   sampling. Adding the paraboloid height to the sampling point is a one-line
   change to a tested, threaded, deterministic kernel.
2. Everything after that spectrum is 2D: a pointwise exponential, a complex
   lens transfer, a squared modulus and a Fourier crop. SIMPLE got a 2D
   complex FFT engine for this on 2026-10-23 (`simulate_scattering`, commit
   `f4cac4773`). The derivative of each step is short and is written by hand,
   with dot-product and finite-difference tests in place of autograd.
3. Ghostbuster's reconstruction (Adam epochs, then a preconditioned nonlinear
   conjugate-gradient solve with a Wiener prior from the half-map FSC) needs
   a preconditioner, a prior, a support mask, half sets and an FSC. SIMPLE
   has all of them in its PCG backend. Only the solver loop is new, and it is
   written as a generic engine that later refinement can reuse.

What it buys, stated plainly. For cryo-EM proteins the expected gain is the
Ewald-sphere curvature, which limits large particles, low voltages and high
resolution. The nonlinear terms (Rytov, multislice) matter for strong
scatterers such as the metal nanoparticles of SINGLE. Every claim is tested
against simulated data with a known answer.

## Rulings

Hans, 2026-10-23:

1. Tomography is out of scope for now: no tilt-series import changes, no
   `refine_tomo` changes, no tomographic datasets.
2. Absolute scale as proposed: the gain from the dose (section 3.4), the gain
   fitted from the data reported as a diagnostic only.
3. No integration yet. A standalone application first: a new shared-memory
   program, `reconstruct3D_wave`. `reconstruct3D`, `refine3D` and their
   strategies are not touched.
4. Alignment stays on SIMPLE's flat references; working implementation first,
   integration later.
5. `simulate_scattering` moves onto the shared operator with corrected
   conventions; Ruben is told.
6. This is the foundation for Ewald-sphere correction and CTF refinement,
   which come later.
7. Names: `wave_model=projection|first_born|rytov|multislice`,
   `ewald=neg|pos|no`; `rec_backend=wave` is kept for the later integration.
8. Multislice stays in this run, as its last implementation phase, for
   single-particle sets of strong scatterers.

## Terms

- **Potential** `V`: the specimen's electrostatic potential on the native
  box `N^3`, voxel size `a` (SIMPLE's `smpd`), in volts.
- **Padded lattice**: the `2N` grid on which SIMPLE's PCG backend does every
  Fourier operation (`pcg_lattice`, `OSMPL_PAD_FAC = 2`). The unknown stays on
  the native box.
- **Gather and scatter**: the Kaiser-Bessel (KB) interpolation that reads the
  volume's transform at a 3D point, and its exact transpose that deposits a
  value there (`forward_plane` and `scatter_plane` in
  `simple_reconstructor_pcg.f90`).
- **Ewald paraboloid**: the transform of the volume seen by an image is not
  the central plane `kz = 0` but `kz = -lambda |k|^2 / 2` (wavelength
  `lambda`), for specter's default "negative" curvature sign.
- **Surface spectrum** `S_i(h,k)`: the volume's transform sampled on the
  posed paraboloid of particle `i`, on the full disk of native image
  frequencies `(h,k)`, both signs of `h`.
- **Exit wave** `psi`: the complex electron wave leaving the specimen.
- **Lens transfer** `T`: the complex objective-lens factor `exp(i gamma)`,
  `gamma` the aberration phase.
- **Forward operator**: `V -> predicted image`. **Linearisation**: the change
  of the predicted image for a small change `dV`. **Adjoint**: the transpose
  of the linearisation applied to an image residual; summed over particles it
  is the gradient of the image misfit, the "backprojection".
- **Weak-phase model**: SIMPLE's current image model, image = CTF (contrast
  transfer function) times the central section, which is the first-order
  limit of every model above.
- **PCG backend**: SIMPLE's reconstruction by preconditioned conjugate
  gradients (`rec_backend=pcg`), which solves the weak-phase least-squares
  problem with the gather and scatter above.
- **FSC, SSNR, ML prior**: the Fourier shell correlation of the two half
  maps; the spectral signal-to-noise ratio `FSC / (1 - FSC)` it implies; the
  per-shell regularisation `1/tau^2` built from it (maximum-likelihood, as in
  RELION).
- **SINGLE**: SIMPLE's programs for metal nanoparticles imaged by
  high-resolution TEM.

## 1. What Ghostbuster computes

`Reconstructor.forward` runs the simulator's chain on the volume being
reconstructed (`ghostbuster/_reconstructor.py`):

1. **Scatter** (`imaging/_specimen.py`, `scattering/_batch.py`). Under
   `rotate_mode="surface"` (`scattering/_surface.py`) the slice sum
   `sum_z F2[V_z](k) F_z(k)` (`F_z` the Fresnel propagator from slice `z` to
   the exit plane) is evaluated as the volume's transform on the posed
   paraboloid, with one type-2 non-uniform FFT of the unrotated volume and a
   phase. Projection samples the flat plane. With `c = i sigma a (sqrt(1 -
   alpha^2) + i alpha)` (`sigma` the interaction constant, `alpha` amplitude
   contrast): projection and Rytov give `psi = exp(IFFT[c S])`, first Born
   `psi = 1 + IFFT[c S]`. Multislice alternates slice transmission
   `exp(i sigma a V_z)` and one-slice Fresnel propagation, and needs the
   rotated volume (`rotate_mode="real"`, trilinear `grid_sample`).
2. **Oversample**: the exponential and the intensity are formed on a grid
   twice as fine (`oversample=2`, the default), since both double the
   bandwidth.
3. **Lens**: `psi_T = IFFT[T FFT[psi]]` with `T = exp(-i chi)` (defocus,
   astigmatism, Cs, phase shift, beam tilt, trefoil, B-factor). The defocus
   refers to the specimen midplane, not the exit face.
4. **Detect**: `I = |psi_T|^2`, Fourier-cropped to the pixel grid, times
   `dose * a^2` for expected counts. No noise, no camera, no coherence
   envelopes in the inverse.

The loss is the squared error in counts, per-particle contrast scale applied.
The reconstruction default is Rytov (`configs/reconstruct.toml`), the cheapest
model with curvature.

The backward direction has two levels. The adjoint is PyTorch's autograd; for
the surface sampling it is a type-1 non-uniform FFT (`nufft.py`). The
reconstruction (`ghostbuster/_converge.py`) runs a few Adam epochs per half,
then minimises

```text
1/2 sum_i || y_i - p_i(V) ||^2 + 1/2 sum_k lambda(k) |V(k)|^2,
lambda(k) = W(k) / SSNR(k),   SSNR = FSC / (1 - FSC) of the Adam half maps
```

by preconditioned Polak-Ribiere conjugate gradients, preconditioner
`1 / (W(k) + lambda(k))`, `W(k)` the data curvature measured per radial band
with random probes, step length from a Gauss-Newton estimate refined on 5% of
the particles. Four iterations from a five-epoch Adam map converge. The map is
then unfiltered by `1 + 1/SSNR`. specter also has a tomogram reconstructor
(Adam only, multislice through slices of the tilted slab); it is out of scope
here.

## 2. What SIMPLE already has

- **A matched gather/scatter pair on the padded lattice**
  (`simple_reconstructor_pcg.f90`, `simple_pcg_lattice.f90`). `forward_plane`
  gathers at `loc = padf * [h,k,0] R`, `scatter_plane` is its literal
  transpose, threaded by h-strided colouring with a serial pass over the
  wrapping rim, bitwise reproducible. `test=pcg_recon` proves the adjoint
  identity to `1e-5`. Deapodisation, symmetry replication over `R S_g`, the
  solve support and the full both-sign-`h` disk are in place.
- **A density, a preconditioner and priors**: `accumulate_absT2`,
  `finalize_density_accum` with the unmeasured-voxel rule
  (`PCG_RHO_SAMPLED_FRAC`), the ML prior from the half-map FSC
  (`build_ml_prior_from_density`), sigma2 whitening, even/odd halves,
  automask, half-map diagnostics.
- **A generic linear PCG engine** (`opt/simple_pcg_solver.f90`): an abstract
  operator with deferred procedures, an adapter per client, the engine owning
  only its vectors. L-BFGS-B exists but its workspace at `n = N^3` is several
  GB. There is no nonlinear CG.
- **The particle read path of reconstruction**: `ptcl_read_session`,
  `prep_rec_observation` and the particle cache (`cache=yes`), used by the PCG
  worker's accumulation (`accumulate_worker_state_half` in
  `simple_rec3D_pcg_strategy.f90`).
- **A 2D complex wave engine** (`image/simple_scattering_simulator.f90`,
  `simulate_scattering`, Ruben, 2026-10-23): owned c2c FFTW buffers, the
  interaction constant, the electron wavelength, first Born and multislice,
  an objective lens with B-factor, Poisson noise. It rotates the whole volume
  per particle (`rotvol`) and does one 2D FFT per slice. It has no adjoint.
- **The CTF phase** (`simple_ctf.f90`): `evalPhSh = pi lambda s^2 (df - 1/2
  lambda^2 s^2 Cs)`, which is `-chi` in Kirkland's (and specter's) sign, so
  specter's `T = exp(-i chi)` is `exp(+i evalPhSh)` in SIMPLE. `evalPhSh` is
  private today.

What is missing: any curvature in the gather and scatter, a complex lens
transfer in reconstruction, the derivative of the wave chain, and a nonlinear
solver.

Two findings on `simulate_scattering`, confirmed against specter in the first
phase. It conjugates both the Fresnel propagator (`exp(+i pi lambda d k^2)`)
and the lens transfer (`exp(+i chi)`) relative to Kirkland, while keeping the
transmission `exp(+i sigma V)`. That is specter's image of the potential with
its sign reversed: a positive potential is bright at underfocus instead of
dark, and the multiple-scattering terms are those of `-V`. specter's own
propagator had the conjugate sign until 2026-10-16, when it was found to walk
each slice the wrong way in depth. Second, its defocus refers to the exit
face; specter and CTF fits refer to the midplane. Its tests (zero and constant
potentials) cannot see either. A numpy check of the first: a weak Gaussian
phase object at 1 um underfocus, 300 kV, changes the centre intensity by
-0.026 under Kirkland's signs (dark) and by +0.027 under the simulator's
(bright).

## 3. Design

### 3.1 The surface spectrum

The gather and scatter move from `reconstructor_pcg` into `pcg_lattice` and
take a curvature `kappa`:

```text
loc = padf * [h, k, kz(h,k)] R_i,      kz(h,k) = kappa * wl * (h^2 + k^2) / (2N)
```

with `wl = lambda / a` (SIMPLE's `ctfvars%wl`) and `kappa = -1, +1, 0` for
`ewald=neg`, `pos` and `no`. At 300 kV, 1 A pixels and `N = 256` the
paraboloid leaves the plane by 0.63 Fourier voxels at Nyquist; at 100 kV,
0.8 A and `N = 400`, by 2.3. With `kappa = 0` the code is the current PCG
code, bitwise. The colouring argument holds unchanged (the map
`(h,k) -> loc` is a rotation of a graph over `(h,k)`); the rim prefilter
`sq_rim` must include `kz`.

`S_i(h,k) = gather(loc) * exp(+2 pi i (h s_x + k s_y) / N)`, the shift on the
model side as `build_transfer` puts it, so the particle preparation of the
PCG backend (`prep_rec_observation`) is reused unchanged.

Symmetry is replication over `R_i S_g` in the gather and the scatter, as in
PCG: the model sees the symmetrised volume, which is what specter's solver
does with its exact point-group average.

Defocus reference: sampling the paraboloid about the box centre makes the
nominal defocus the midplane's with no extra phase. specter's propagation to
the exit face and its midplane defocus shift cancel to exactly this for even
boxes (numpy, first Born on a random 16^3 volume: equal to `2e-15`).

### 3.2 Image formation on the padded grid

The padded lattice is `2N` per side, which is specter's default image
oversampling, so the fine grid comes for free. For particle `i`:

```text
u   = IFFT_2N[ embed(S_i) ]                       complex; S_i zero beyond the native disk
psi = exp(c u)        projection (kappa = 0), Rytov
psi = 1 + c u         first Born
psi_T = IFFT_2N[ T * FFT_2N[psi] ]                T on the fine-grid frequencies
I_hat = crop_N( FFT_2N[ |psi_T|^2 ] )             native disk, the detector's band limit
p_i   = -g * I_hat  (k /= 0)                      in the units of SIMPLE's particles
```

`c = i sigma a (sqrt(1 - alpha^2) + i alpha)` with `alpha = fraca`. The minus
sign and the gain `g` map relative intensity to SIMPLE's contrast-inverted,
noise-normalised particles (section 3.4). The residual is formed on the
native disk in Fourier space, whitened by the particle's sigma2, exactly where
PCG forms its right-hand side.

The weak-phase limit fixes every convention. With `V -> eps V`,
`p_i / eps -> 2 g sigma a * CTF * S_i(flat)` where `CTF = sin(evalPhSh +
phshift + amp_contr_const)` is SIMPLE's CTF: the new operator collapses onto
`forward_plane` times `build_transfer` (numpy, with amplitude contrast 0.1
and Cs: holds to second order in the phase). This identity is the first test
written, and it ties orientation, shift, defocus sign, astigmatism angle,
amplitude contrast and scale to code that refinement has used for years. It
also shows that the amplitude contrast of the wave models is `fraca`, the
same quantity the project already carries per particle.

The oversampling factor is a constant (2, as specter) until the Dell run has
measured the error and the speed of 1 against 2; a parameter is added only if
the measurement makes 1 worth offering.

### 3.3 Lens transfer

`T = exp(i (evalPhSh + phshift)) * exp(-B s^2 / 4)`, the phase from the `ctf`
type itself (a public `eval_phase` that returns what `evalPhSh` computes, the
phase a required argument as `phase_shift_ctf_policy.md` demands). Amplitude
contrast lives in `c`, not in `T`. Beam tilt and trefoil are not modelled:
SIMPLE's project has no fields for them yet, and they come with CTF
refinement.

### 3.4 Units and the scale link

The nonlinear models need the potential in volts and the images as relative
intensity (vacuum = 1); the weak-phase model does not. Two conversions:

- **Images**: SIMPLE particles are contrast-inverted and normalised to unit
  background noise. Following specter's `image_units="normalized"`, relative
  intensity is `1 - x / sqrt(D a^2)` with `D` the dose in e/A^2 (SIMPLE's
  `total_dose`), so `g = sqrt(D a^2)`. The gain the data prefer (a one-scalar
  fit of the quadratic term) is reported, not used.
- **Maps**: by the weak-phase identity, a SIMPLE map at the data-quotient
  scale is `2 g sigma a` times the potential in volts, up to SIMPLE's FFT
  normalisation, which the test pins. The solver works in volts; the program
  writes its maps at SIMPLE's scale, so they compare directly with
  `reconstruct3D` maps and a later integration sees no amplitude jump (the
  scale-continuity contract of `reconstruct3D_pcg_policy.md`).

### 3.5 Linearisation and adjoint

Each step above is a small operator with a hand-written linearisation and
transpose: the KB gather (transpose: the scatter that exists), the embedding
and crop (transposes of each other up to the FFT scale), the c2c FFTs, the
pointwise exponential (`d psi = c psi du`), the lens multiply (`conj(T)`), the
squared modulus (`dI = 2 Re(conj(psi_T) d psi_T)`). The adjoint of a particle
is a 2D chain ending in one scatter of a complex plane at the paraboloid
points, the same deposit the PCG right-hand side makes. The gradient is the
real part of the accumulated lattice transformed back, as in
`apply_adjoint_all`.

### 3.6 Reconstruction: one solver for every model

Minimise, per state and half set,

```text
L(V) = 1/2 sum_i sum_(h,k) |y_i - p_i(V)|^2 / sigma2_i(shell)  +  1/2 sum_q lambda(q) |V(q)|^2
```

on the support, with:

- **Start**: the half's flat PCG map, solved in the same process with
  `reconstructor_pcg`'s public API from the same particle pass, converted to
  volts.
- **Prior**: SIMPLE's ML prior from the FSC of the two start maps, built by
  the existing PCG code. It is the same quantity as specter's
  `W(k)/SSNR(k)`.
- **Preconditioner**: `1 / (rho + lambda)`, `rho` the PCG density of the
  start solve. It is the diagonal of the linearised data term, which specter
  measures with random probes.
- **Iteration**: preconditioned Polak-Ribiere conjugate gradients. The step
  length uses the exact data slope and the Gauss-Newton curvature `||J d||^2`
  from the linearisation on a particle subset, refined by a parabola through
  the loss on the same subset, as specter does. Default four iterations
  (`maxits`).
- **Output**: the regularised map and the unfiltered half map (`1 + 1/SSNR`,
  exact for a Wiener solve), and the FSC of the unfiltered pair.

For first Born with a linearised intensity the objective is quadratic, and
this is linear CG on the exact two-cap least-squares problem up to the step
length taken from a subset. For Rytov, projection and multislice it is the
nonlinear solve. One loop serves all of them.

### 3.7 Multislice

Multislice is not one sample of the transform. The rotated cube is formed
from the lattice: its spectrum on native plane `l` is the gather at
`padf [h,k,l] R`, then one 3D inverse FFT. Its transpose is the scatter, so
the pair is exactly adjoint, which neither `rotvol` nor trilinear resampling
gives. The recursion runs on the `2N` fine grid with Kirkland's signs,
transmission `exp(c v_j)` and propagator `exp(-i pi lambda a |k|^2)`, the exit
wave referred to the midplane. The adjoint is the reverse sweep, `conj(P)` and
`conj(t_j)`, with the wave stored at checkpoints every `ceil(sqrt(N))` slices
and recomputed in between. Cost per particle at `N = 256`: of order a second
per core, which suits nanoparticle boxes and small sets, not 100 000-particle
cryo-EM sets. Rytov serves those.

### 3.8 Parallel execution and performance

Particles are processed in batches: each thread owns its 2D work fields and
c2c plans (made once, in the FFTW planner critical section, as
`scattering_simulator` does) and runs the 2D chain of one particle; the
batch's gradient planes are then scattered with the coloured scatter, the
`accumulate_batch` pattern of PCG. Single precision fields, double precision
losses and dot products.

Order of magnitude at `N = 256`: a Rytov gradient is about eight 512^2
complex FFTs plus a gather and a scatter, some tens of milliseconds per
particle per core, so a pass over 100 000 particles on the Dell's 24 cores
takes minutes and four iterations per half ten to thirty times a PCG
accumulation. First Born with a linearised intensity needs no FFT at all
(`p = A S(k) + B conj(S(-k))` on the native disk) and costs about two PCG
right-hand sides per pass.

### 3.9 What maps onto what

| Ghostbuster / specter | SIMPLE |
| --- | --- |
| Non-uniform FFT of the unrotated volume on the posed surface, and its type-1 adjoint | KB gather and scatter of the PCG lattice with the paraboloid height |
| Volume rotation for multislice (`grid_sample`) | Lattice gather of every plane, one 3D FFT; scatter for the adjoint |
| PyTorch FFTs, complex tensors | FFTW c2c fields from `scattering_simulator`, `image` r2c |
| `Aberration` (`exp(-i chi)`) | `ctf` phase, `exp(+i evalPhSh)`, B-factor |
| Wavelength, interaction constant | One home in `simple_math_ctf` |
| Autograd | Hand-written linearisation and transpose per step, tested |
| Adam epochs | None: start from the flat PCG map |
| Measured curvature preconditioner | PCG density |
| Wiener prior from Adam half maps | ML prior from the FSC of the start maps |
| Gold-standard halves, FSC, auto mask | Even/odd, half-map diagnostics, `automask3D` / `mskdiam` |
| Exact point-group average | Symmetry replication in gather and scatter |
| Lightning DDP over GPUs | OpenMP in one process; distributed execution at integration |
| `.cs` / `.star` input | The SIMPLE project; specter data enter through RELION star files |

Not ported: Adam and Lightning, the kinematic model, the linear "ctf" model
(SIMPLE's existing model is that), anisotropic magnification, per-particle
contrast scale, the absorption map, inelastic absorption, beam tilt and
trefoil, the tomogram reconstructor, and the simulator-only parts (ice,
crowding, damage, detector MTF, coincidence loss).

## 4. Architecture: how the code fits into SIMPLE

The rule set is `doc/policies/layering_and_project_ownership_policy.md`,
`doc/policies/compile_time_policy.md`, `doc/policies/ui_layer_policy.md` and
the `simple-modern-fortran` skill. This section says where every piece goes,
what it may import, who owns what, and which existing code it changes.

### 4.1 Where each piece lives

From the bottom of the import graph to the top:

| Module | Layer | What it holds | New imports | Used by |
| --- | --- | --- | --- | --- |
| `src/utils/math/simple_math_ctf.f90` | stateless utility | gains `electron_wavelength(kv)` and `electron_interaction_constant(kv)`, the one home of both | none | `ctf`, `wave_imaging`, the simulator |
| `src/utils/math/simple_rnd.f90` | stateless utility | gains the Poisson sampler that sits in the simulator today | none | the simulator |
| `src/main/ctf/simple_ctf.f90` | domain | `ctf` gains a public `eval_phase` (what the private `evalPhSh` computes); its private `kV2wl` calls the shared wavelength | `simple_math_ctf` | `wave_imaging` |
| `src/main/volume/simple_pcg_lattice.f90` | domain | gains the gather, the coloured scatter and the rim test, moved from `reconstructor_pcg`, with the curvature `kappa` | none beyond its own | `reconstructor_pcg`, `wave_projector`, `wave_multislice` |
| `src/main/image/simple_wave_field.f90` (new) | domain | type `wave_field`: a 2D complex field on an `n x n` grid with its own c2c FFTW plans, forward and backward transforms, embedding of a native-disk spectrum into the `2N` grid and the band-limited crop back, both with their transposes | `simple_fftw3` | `wave_imaging`, `wave_multislice` |
| `src/main/image/simple_wave_imaging.f90` (new) | domain | type `wave_imaging`: the 2D chain of one image for projection, first Born and Rytov (exponential, lens transfer, intensity, crop), its linearisation and its adjoint; owns one thread's work fields | `wave_field`, `ctf`, `simple_math_ctf` | `wave_projector`, the simulator |
| `src/main/volume/simple_wave_projector.f90` (new) | domain | type `wave_projector extends pcg_lattice`: the volume's padded transform, the particles' cached geometry (rotation, shift, CTF parameters), one `wave_imaging` per thread, the lattice gradient accumulator; for a batch of observation planes: loss, gradient accumulation, `||J d||^2`, loss changes along a direction, and the forward image alone | `pcg_lattice`, `wave_imaging`, `ori`/`oris`, `sym` | the strategy, the simulator |
| `src/main/volume/simple_wave_multislice.f90` (new, last phase) | domain | type `wave_multislice extends pcg_lattice`: rotated cube from the lattice, recursion, reverse sweep with checkpoints, the same batch interface as `wave_projector` | `pcg_lattice`, `wave_field` | the strategy, the simulator |
| `src/main/opt/simple_nlcg_solver.f90` (new) | numerical engine | abstract `nlcg_objective` (deferred: size, loss and gradient, curvature along a direction, loss changes along it, preconditioner, dot, projection onto the model's constraints) and `nlcg_solve`, preconditioned Polak-Ribiere with the step rule of section 3.6, on the pattern of `simple_pcg_solver.f90` | none of SIMPLE's domain | the strategy; later CTF and pose refinement |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | strategy | gains the observation-batch frame of reconstruction (selection of a state's half, batched reading through `ptcl_read_session`, `prep_rec_observation` or the pass's records, native Fourier planes), extracted from the PCG worker so that the PCG worker and the wave objective share one copy | none | PCG worker, wave strategy |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` (new) | strategy | `exec_rec3D_wave(params, build)`: the program body; and `rec3D_wave_objective extends nlcg_objective`, whose deferred procedures run the particle passes on the caller's `params` and `build` | the volume and opt modules above, `reconstructor_pcg`, `simple_matcher_ptcl_io`, `simple_halfmap_diagnostics`, sigma2 | the commander |
| `src/main/commanders/simple/simple_commanders_rec.f90` | commander | `commander_rec3D_wave`: normalizes `cline`, one `params%new` and builder, calls `exec_rec3D_wave`, ends | the strategy | the router |
| `src/main/exec/simple_exec_refine3D.f90` | exec | the `reconstruct3D_wave` case beside `reconstruct3D` | the commander | `simple_exec` |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90` | ui | the program descriptor, category `reconstruct3d` | none | the UI registry |
| `src/main/params/simple_parameters*.f90` | params | the keys of section 4.5 and their checks | none | everything above the domain |
| `src/main/commanders/simple/simple_commanders_sim.f90`, `src/main/ui/simple/simple_ui_sim.f90` | commander, ui | `simulate_scattering` on `wave_projector` and `wave_multislice` | the volume modules | `simple_exec` |
| `src/main/image/simple_scattering_simulator.f90` | domain | deleted once the simulator runs on the shared pieces (the multislice phase); its FFT buffers become `wave_field`, its physics `wave_imaging` and `wave_multislice` | | |

Testers sit beside each module (`simple_<module>_tester.f90`).

### 4.2 Rules every piece follows

- **Imports point downward.** Utilities, then `ctf`, then `image`, then
  `volume`; `opt` imports no domain module; the strategy imports all of them;
  the commander imports the strategy. No domain module imports a strategy
  module: that is why the particle passes live in the strategy, not in
  `wave_projector`, and why this work does not add to the known deviation of
  domain modules using `ptcl_read_session` (layering policy, section 8).
- **Kinds of step.** The program body is a strategy routine on
  `(params, build)`; the commander is a thin `execute(cline)` (layering
  policy, section 3). Numerics live in domain modules. Nothing calls a
  commander from a strategy.
- **One copy of each frame.** The observation-batch frame is extracted from
  the PCG worker, not copied; the gather and scatter are moved, not copied;
  the wavelength has one home. The PCG worker's results are byte-identical
  after the extraction (`simple_test_exec test=identity` on a
  `reconstruct3D rec_backend=pcg` case).
- **Encapsulated classes.** Each new stateful module holds one type named
  after it, with private components, `new`/`kill` (symmetric, idempotent) and
  behaviour as type-bound procedures (`simple-modern-fortran`). Large
  components are `allocatable`, allocated in `new` and released in `kill`.
  No type holds a copy of `parameters`; it copies the fields it reads
  (`box_crop`, `smpd_crop`, `kv`, ...).
- **Threads.** Per-thread objects (`wave_imaging`, its FFT plans) sit in
  allocatable arrays sized by `omp_get_max_threads()` at the point of use,
  never by a stored thread count (`thread_budget_and_safety.md`). FFTW plans
  are made and destroyed in the `fftw_planner` critical section. The scatter
  is the existing coloured scatter with its serial rim pass; losses and dot
  products are summed per particle in double precision and reduced in
  particle order. Results are therefore identical at any thread count, which
  the tests check at 1, 11, 16 and 24 threads, repeated.
- **Errors** go through `THROW_HARD`, built in a string variable when long
  (the macro's argument is never continued across lines). The simulator's
  `error stop` calls go with the move.
- **Compile time.** `use ..., only:` with the names used; no umbrella module
  in new code and nothing added to an umbrella; no new CMake target (the
  library globs its sources); testers declare large objects
  `class(T), allocatable`.
- **Release-4 rules.** No compatibility shims: `simulate_scattering`'s
  `scattering_model` key becomes `wave_model` and its `absorption` key is
  retired in favour of `fraca`, with no alias. No `goto`, no reliance on
  integer overflow, no short-circuit on an allocatable's element. Dead code
  found on the way is deleted.

### 4.3 The program `reconstruct3D_wave`

```text
simple_exec prg=reconstruct3D_wave projfile=<project> wave_model=rytov ewald=neg \
    mskdiam=<A> [mskfile=<mask>] [pgrp=<group>] [box_crop=<px>] [maxits=4] [nthr=<n>]
```

- **Input**: the project's `ptcl3D` field (poses, shifts, states, even/odd
  labels), its stacks and per-particle CTF (`dfx`, `dfy`, `angast`,
  `phshift`, `kv`, `cs`, `fraca`), `total_dose`, and the sigma2 of the
  particles when `objfun=euclid` (as PCG reads them).
- **Steps**: per state, the two flat PCG start solves, the ML prior from their
  FSC, then the wave solve of each half; the unfiltered pair's FSC and the
  half-map diagnostics; the maps written at SIMPLE's scale.
- **Project contract**: reads the project once (the builder); writes the
  maps and FSC files and registers them in the `out` segment
  (`add_vol2os_out`, `add_fsc2os_out`); one write of the `out` segment at the
  end. The particle field is not written.
- **Particle reads**: one pass per solver iteration plus the start pass, all
  through the shared observation-batch frame, so `cache=yes` serves the
  repeated passes from SIMPLE's particle cache. No new cache.
- **Execution**: shared memory only. `nparts` is refused with a message
  naming the later integration.

### 4.4 Existing code that changes, and what each change must keep

| Change | Must keep |
| --- | --- |
| Gather and scatter moved into `pcg_lattice` with `kappa` | Every PCG result bitwise at `kappa = 0`: `pcg_recon`, `lib_reconstruction`, `unit_reconstruction` unchanged, and the identity check of a `reconstruct3D rec_backend=pcg` case |
| Observation-batch frame extracted from the PCG worker | The same identity check |
| `ctf` exposes `eval_phase`; `kV2wl` calls the shared wavelength | The CTF suite unchanged (same formula, so bitwise) |
| `simulate_scattering` on the shared operator | Intended change of output (conventions corrected); its tests rewritten to be sign- and defocus-sensitive |
| `flex_pcg_t` (extends `pcg_lattice`) | Unchanged; its suites pass |

### 4.5 Parameters and UI

New keys, declared in `simple_parameters.f90` and parsed and checked there:
`wave_model` (`projection|first_born|rytov|multislice`, default `rytov`) and
`ewald` (`neg|pos|no`, default `neg`). Existing keys reused: `total_dose`,
`maxits`, `mskdiam`, `mskfile`, `pgrp`, `box_crop`, `objfun`, `nthr`, `cache`.
`simulate_scattering`'s `scattering_model` becomes `wave_model` (its values
`first_born` and `multislice` keep their spelling); its `absorption` key is
retired, since the amplitude contrast of the wave models is `fraca`; its
`noise_sdev` stays. The `reconstruct3D_wave` descriptor goes into
`simple_ui_reconstruct3D.f90`, shared input metadata into
`simple_ui_params_common.f90`. `scripts/check_test_registry.py` and the UI
visibility tester pass.

### 4.6 Tests

Unit testers beside their modules, registered as sub-suites with `add_suite`
in `simple_commanders_test_class.f90` and listed in the `suite=` help of
`simple_test_ui_class.f90`: the fast ones under `unit_reconstruction` (and
`unit_numerics` for the solver engine on analytic objectives), the
ground-truth ones under `lib_reconstruction`. No new CTest entry, so
`SIMPLE_CTEST_BUDGET` does not change. A high-level test of the program is
decided at integration.

### 4.7 What waits for integration

`rec_backend=wave` in `reconstruct3D` and `refine3D` (through
`simple_rec3D_service`), distributed execution (workers publishing raw
gradient and curvature sums, the master stepping), the NU filter, trailing
chains, state weight sets and fractional updates, Ewald correction by double
insertion in the PCG and gridding backends, CTF refinement, pose refinement,
and tomography. The design keeps the way open: the objective's particle pass
is the piece a worker would run, the solver needs only gradient sums and
subset curvatures, and the engine in `opt` is generic.

## 5. Validation

### 5.1 SIMPLE's own tests

Fast sub-suites (box 16 to 32, well under a second each):

- Adjoint identity `<J dV, r> = <dV, J^T r>` for every linear step, for the
  whole linearisation at a random point, and for the multislice sweep:
  `1e-5`.
- Finite-difference check of the loss gradient for each model, double
  precision accumulation: `1e-3`.
- Weak-phase identity against `forward_plane` and `build_transfer` (3.2).
- Curvature off: the surface gather with `kappa = 0` is `forward_plane`,
  bitwise.
- A point scatterer at height `z` images as the flat model at defocus
  `df + z` for `ewald=neg` and `df - z` for `pos`: the size and sign of the
  curvature, with no inverse crime.
- Model hierarchy: Rytov with `kappa = 0` is projection; first Born is the
  linearisation of Rytov; multislice approaches first Born as the potential
  goes to zero.
- The solver engine on analytic objectives (a quadratic, where it must equal
  linear PCG, and a nonconvex one with a known minimum).
- specter golden cases: a few 16^3 and 24^3 phantoms with known poses, CTF
  and both curvature signs, imaged by specter in the first phase and embedded
  in the tester. specter's own committed golden files cannot be used as they
  are: they include beam tilt and trefoil.
- Library tier: a phantom imaged by multislice, solved with Rytov and first
  Born, map against truth, beating the flat PCG map on a curvature-limited
  geometry.

### 5.2 specter as the oracle

specter has its own suite (111 test files, golden outputs, parity tests
between its simulator and its reconstructor). On the Dell it serves as the
oracle: installed outside SIMPLE, its tests for the parts used here run first,
then an export script in the run's scratch writes reference cases. Nothing in
SIMPLE calls Python; the references are files.

Conventions cross the boundary through RELION star files, read and written by
each package's own code: SIMPLE's `export_relion` gives specter its poses and
CTF (`star_path`), specter's simulated stacks come back with a star file.

Parity measured: images, exit waves and gradients (specter's autograd against
SIMPLE's adjoint, for fixed residual images) at box 64 and 128, random poses,
shifts and astigmatism, all four models, both curvature signs. The tolerance
is the KB lattice error, about `1e-3` relative (specter's non-uniform FFT runs
at `1e-6`); it is measured at zero curvature first, where both sides are exact
by another route.

### 5.3 Ground-truth reconstructions

Datasets simulated by specter with multislice (so no model reconstructs its
own data), truth poses, Poisson noise:

- **Cryo-EM, ordinary**: a protein of about 150 A at 300 kV, 1 A pixels,
  6 000 particles. Expected: all models agree with the flat PCG map.
- **Cryo-EM, curvature-limited**: a particle of about 300 A at 100 kV, where
  the periphery loses signal from about 5 A (`sqrt(2 lambda D)`). Expected:
  first Born and Rytov beat flat PCG beyond that; the wrong curvature sign
  loses to the right one.
- **Strong scatterer**: a 2 to 3 nm Pt nanoparticle, single-particle set,
  300 kV, HRTEM defocus. Expected: multislice recovers atoms the linear map
  misplaces or misses.

Scored against truth: FSC and real-space correlation in a mask for the
protein sets (more than 2 000 particles), atom recall, precision and position
error for the nanoparticle (SINGLE's atom detection). specter's own converged
reconstruction of each set is in the comparison.

## 6. Risks

- **Accuracy of the KB lattice.** About `1e-3` against specter's `1e-6`. It
  sets the parity tolerance; it is far below the noise of any particle.
- **Scale of real data.** Dose-weighted sums do not have Poisson background
  noise, so `g` from the dose is approximate. It only scales the nonlinear
  terms, which are small for proteins; the fitted gain shows how far off it is.
- **Local minima for strong scatterers.** specter reports that a specimen with
  17 rad peak phase converged from a 15 A low-passed start and not from zero.
  The start here is a linear map; a low-pass schedule is the fallback.
- **Repeated particle passes.** Four to five reads per half; the particle
  cache makes them cheap when `cache=yes`, otherwise they cost what a PCG
  accumulation costs each time.
- **Thread safety of the curved scatter.** The colouring still separates
  windows; the rim prefilter changes. Covered by the repeated thread tests.

## 7. Execution

### 7.1 Common rules

- Build incrementally in the run's build tree (`make -C build -j`, then
  `scripts/run_fast_gate.sh build` and `make -C build install`); run the
  touched area's `simple_test_exec test=unit_<area>` and the library suites
  this plan names. Start long jobs in the background.
- New OpenMP code is tested at 1, 11, 16 and 24 threads, each twice, for
  identical results.
- Bugs found are fixed in the phase that finds them. No floor or tolerance
  is relaxed unless this plan or a ruling says so; tolerances come from the
  measurements of Phase 0 and Phase 2, recorded with their evidence.
- Timing on the Dell only, on the Release build, nothing else heavy running,
  twice each way.
- Plan text, progress entries and the report read on their own: SIMPLE terms
  explained, acronyms spelled out, no item IDs.

### 7.2 Phases

0. **Oracle and conventions.** specter installed outside the repository on the
   Dell; its tests for scattering, aberrations, the surface sampling and the
   reconstructor pass. The SIMPLE-to-specter convention table through star
   files, checked with an off-centre two-atom phantom at random poses. The
   reference cases of 5.1 and 5.2 and the three datasets of 5.3 written, with
   their truth. `simulate_scattering` compared with specter on the same
   potential and pose, settling the two findings of section 2.
   Exit: the oracle's tests pass; the convention table, the references and the
   datasets are on disk with a manifest; the two findings confirmed or
   refuted with numbers.
1. **Foundations.** Gather and scatter in `pcg_lattice` with `kappa`; the
   observation-batch frame in `simple_matcher_ptcl_io`; `eval_phase` and the
   shared wavelength and interaction constant; `wave_field`; the Poisson
   sampler in `simple_rnd`; the solver engine in `opt`. Exit: the identity
   obligations of 4.4 met; the new testers pass (adjoints of embed and crop,
   FFT round trips, solver engine on analytic objectives); fast gate green.
2. **Single-scatter operator.** `wave_imaging` and `wave_projector`: forward,
   linearisation and adjoint for projection, first Born and Rytov, batched,
   threaded. The fast tests of 5.1. The parity campaign of 5.2. Thread
   identity. Cost per particle against a PCG accumulation, and the error and
   speed of oversampling 1 against 2. Exit: every test of 5.1 for these models
   passes; parity within the tolerance measured at zero curvature; results
   identical across thread counts; the cost table recorded.
3. **The program.** `reconstruct3D_wave` end to end (parameters, UI, router,
   commander, strategy with the objective over particle passes);
   `simulate_scattering`'s single-scatter models moved onto `wave_projector`
   (it gains projection and Rytov), its remaining multislice path corrected to
   Kirkland's signs and the midplane defocus, its keys renamed and retired, its
   tests rewritten to be sign- and defocus-sensitive. Ground truth on
   the two cryo-EM datasets of 5.3 against gridding, flat PCG and specter.
   Exit: the program runs on both datasets; the comparisons recorded; the
   curvature-limited set shows the expected ordering or the reason it does
   not is found and fixed.
4. **Multislice.** `wave_multislice` with its tests (adjoint of the sweep,
   finite differences, weak limit to first Born, specter parity);
   `wave_model=multislice` in the program; `simulate_scattering`'s multislice on
   `wave_multislice` and `simple_scattering_simulator.f90` deleted; ground
   truth on the nanoparticle set. Exit: tests pass; the nanoparticle
   comparison recorded.
5. **Documentation and report.** A policy for the wave model
   (`doc/policies/3D/wave_model_policy.md`: conventions table, scale link,
   architecture of section 4, tests); user-guide pages for
   `reconstruct3D_wave` and the updated `simulate_scattering`; rows in the
   test-environment policy; Ruben's note updated; the report beside this plan.
   Exit: the documents read on their own; `check_test_registry.py` and
   `check_descr.py` pass.

### File table

A phase edits only the files listed for it; testers of a listed module are
included. A file not listed gets a row with a one-line reason before it is
edited.

| File | Purpose | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/ghostbuster_operators.md` | this plan: progress, measured tolerances | each |
| `src/utils/math/simple_math_ctf.f90` | electron wavelength and interaction constant | 1 |
| `src/utils/math/simple_rnd.f90` | Poisson sampler | 1 |
| `src/main/ctf/simple_ctf.f90` | public `eval_phase`; `kV2wl` on the shared wavelength | 1 |
| `src/main/volume/simple_pcg_lattice.f90` | gather and scatter with curvature | 1, 2, 4 |
| `src/main/volume/simple_reconstructor_pcg.f90` | calls the moved gather and scatter | 1 |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | shared observation-batch frame | 1, 3 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | PCG worker on the shared frame | 1 |
| `src/main/image/simple_wave_field.f90` | new: 2D complex fields and plans | 1, 2, 4 |
| `src/main/opt/simple_nlcg_solver.f90` | new: nonlinear CG engine | 1, 3 |
| `src/main/image/simple_wave_imaging.f90` | new: the 2D chain, linearisation, adjoint | 2, 3 |
| `src/main/volume/simple_wave_projector.f90` | new: single-scatter operator | 2, 3, 4 |
| `src/main/volume/simple_wave_multislice.f90` | new: multislice operator | 4 |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | new: program body and objective | 3, 4 |
| `src/main/commanders/simple/simple_commanders_rec.f90` | `commander_rec3D_wave` | 3 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `simulate_scattering` on the shared operator | 3, 4 |
| `src/main/image/simple_scattering_simulator.f90` | reduced in 1 and 3, deleted in 4 | 1, 3, 4 |
| `src/main/exec/simple_exec_refine3D.f90` | router case | 3 |
| `src/main/exec/simple_exec_sim.f90` | if the simulator's commander changes name | 3 |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90` | program descriptor | 3 |
| `src/main/ui/simple/simple_ui_sim.f90` | `simulate_scattering` keys | 3 |
| `src/main/ui/simple_ui_params_common.f90` | shared metadata of `wave_model`, `ewald` | 3 |
| `src/main/params/simple_parameters.f90` | new keys; renamed and retired keys | 3 |
| `src/main/params/simple_parameters_parse.f90` | their registration | 3 |
| `src/main/params/simple_parameters_phases.f90` | their checks | 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration | 1, 2, 3, 4 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` help | 1, 2, 3, 4 |
| `doc/user_guide/simulate_scattering.md` | updated keys and conventions | 3, 5 |
| `doc/user_guide/reconstruct3D_wave.md` | new user-guide page | 5 |
| `doc/user_guide/user_guide.md` | index entry | 5 |
| `doc/policies/3D/wave_model_policy.md` | new policy | 5 |
| `doc/policies/test_environment_policy.md` | rows for the new sub-suites | 5 |
| `doc/implementation_notes/planned/ghostbuster_scattering_simulator.md` | Ruben's note: conventions corrected, simulator on the shared operator | 5 |
| `doc/implementation_notes/planned/ghostbuster_operators_report.md` | the run's report | 5 |

## 8. Progress

(none yet)
