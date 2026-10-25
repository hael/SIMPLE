# Ghostbuster's forward and adjoint operators in SIMPLE

Status: implemented by the unattended run on the Dell (section 7), Phases 0
to 5, 2026-10-23 to 2026-10-24, uncommitted against base `5111068bf`; the
progress and evidence are in section 8, the summary in
`ghostbuster_operators_report.md` beside this plan, the contract of the code in
`doc/policies/3D/wave_model_policy.md`. Written 2026-10-23, rulings of the same
day applied (section "Rulings"), from specter (the package that contains
Ghostbuster) at `5786651` in `~/src/specter` and SIMPLE at `630c2604a`.

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

Phase 0 measured both by running `simulate_scattering` (base `5111068bf`) and
specter on the same potential and pose (identity pose, 300 kV, Cs 2.7 mm, 1 um
nominal underfocus, a 96-voxel box at 1 A, no noise, no B-factor, both on the
volume's own grid). In plain words: the simulator draws the image of the
negated potential, and its defocus refers to the exit face, half a box
downstream of the midplane that specter (and a CTF fit) refer to. With no
amplitude contrast, the simulator's image of `V` correlates with specter's
image of `-V` at 0.986 (myoglobin, first Born), 0.984 (myoglobin, multislice)
and 1.0000 (a Gaussian phase object), and with specter's image of `V` at -0.93
and -0.90; the phase object (a Gaussian of 3 A standard deviation and 1 V
peak, weaker than the numpy check's) is brighter at its centre by 0.0026 in
the simulator and darker by 0.0026 in specter. The best match lies at a
specter defocus 48.6 to 51.4 A below the simulator's nominal one, half the box
depth being 48 A: the simulator refers its defocus to the exit face. One
refinement of the first finding: the simulator's absorption keeps Kirkland's
sign. With amplitude contrast 0.07 on both sides, specter's image of `-V`
matches only at 0.952 to 0.959 (myoglobin), while specter with the phase of
`-V` and the absorption of `+V` (amplitude contrast `-0.07` on `-V`) matches
at 0.983 to 0.984 (myoglobin) and 1.0000 (phase object), at the same
half-depth offset. So only the phase is inverted; the corrected simulator
flips the propagator and the lens, not the absorption. Script and numbers:
`scratch/keep/oracle/compare_simulate_scattering.py`,
`scratch/keep/oracle/simulate_scattering_vs_specter/report.json` under the run
directory `/home/elmlundho/agent_runs/wave_operators`.

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

Measured in Phase 2: the cancellation holds for specter's negative curvature
sign; for its positive sign specter refers the defocus one slice off the box
centre (an even box's centre slice is `nz/2 + 1` slices from that sign's exit
face, `nz/2` for the other), so specter's positive-sign images equal SIMPLE's
at one voxel less defocus (images within 1e-4 with that shift, 3e-3 to 9e-3
without). SIMPLE keeps the symmetric reference about the box centre for both
signs. The direction is pinned by a point scatterer: with `ewald=neg` a point
`z0` voxels along `+z` images as the flat model at defocus `df + z0` (the test
of section 5.1), and the parity with specter (section 8, Phase 2) confirms
`kappa = -1` for specter's negative sign: the opposite sign is off by 0.2 to
0.7 in images and gradients.

### 3.2 Image formation on the padded grid

The padded lattice is `2N` per side, which is specter's default image
oversampling, so the fine grid comes for free. For particle `i`:

```text
u   = IFFT_2N[ embed(S_i) ]                       complex; S_i zero beyond the native disk
psi = exp(c u)        projection (kappa = 0), Rytov
psi = 1 + c u         first Born
psi_T = IFFT_2N[ T * FFT_2N[psi] ]                T on the fine-grid frequencies
I_hat = crop_N( FFT_2N[ |psi_T|^2 ] )             native disk, the detector's band limit
p_i   = -g * I_hat / N^2  (k /= 0)                in the units of SIMPLE's particles (1/N^2: the unnormalised FFT)
```

`c = i sigma a (sqrt(1 - alpha^2) + i alpha)` with `alpha = fraca`. The minus
sign and the gain `g` map relative intensity to SIMPLE's contrast-inverted,
noise-normalised particles (section 3.4). The residual is formed on the
native disk in Fourier space, whitened by the particle's sigma2, exactly where
PCG forms its right-hand side.

Three details measured in Phase 2 belong to this chain. Rytov takes its
exponential where the wave leaves the specimen, as specter does and as the
Rytov approximation is defined: on the curved models the field is carried
half a box downstream to the exit face (Kirkland's propagator
`P = exp(-i pi lambda (N a / 2) |k|^2)` on the native disk), exponentiated
there, and carried back to the midplane for the lens (`P^-1` folded into
`T`); for first Born the two cancel exactly, and without curvature the model
has no depth, so Rytov with `kappa = 0` stays projection. Taking the
exponential at the midplane instead differed from specter by 2% (a protein of
0.7 rad) to 17% (a Pt particle of 2.6 rad) in images; at the exit face the
exit waves agree to `8e-6`. The phase shift of a phase plate acts on the
scattered waves only: `T(0) = 1` for the unscattered beam (a phase common to
all frequencies would vanish from the intensity). And SIMPLE's wavelength is
the CTF's `12.26 / sqrt(V + 0.9784e-6 V^2)`, 3.4e-4 short of the exact
relativistic value specter uses; the defocus phase `pi lambda df |k|^2`
reaches some 300 rad at Nyquist for 2 um at 300 kV, so the two packages'
lens phases differ by about 0.1 rad there, and by about 0.5 rad at a Pt
particle's Nyquist through Cs. Parity is measured with specter given SIMPLE's
constants; the constant itself is SIMPLE's, shared with the CTF by design.

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
  fit of the quadratic term) is reported, not used. `a` is the particles'
  native pixel size: the reconstruction normalises each particle to unit
  background noise at its native box before any Fourier crop, and a crop
  keeps the values, so `g = sqrt(D) smpd` also at `box_crop` (Phase 3). The
  fit: with `l` the linear part of the prediction (the chain linearised at
  zero along the solution) and `q = p - l` its nonlinear part, the data fix
  `g V` through `l`, so `y - l = (g / g_data) q` and the fitted gain is
  `g sum w|q|^2 / sum w Re(conj(y - l) q)`. It is meaningful for Rytov and
  projection, whose nonlinear terms are those of a phase object; first Born
  lacks the `-phi^2/2` of the exponential, and its fit is not a gain.
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

Measured and decided in Phase 3. The PCG backend's normal operator, density
and prior carry a factor `boxpd^3` (its adjoint folds the scattered lattice
without the `1/boxpd^3` of the transform's transpose; its own solution does not
see the factor), so the prior in volts is `kappa^2 P / boxpd^3` with
`kappa = 2 g sigma a N` the scale of a PCG map in volts; the projector's tester
checks the link to `2e-6`. The preconditioner is PCG's own (a constant scale of
it leaves the iterates unchanged). The start is the half's base PCG solve
(`maxits_pcg`) shrunk by the ML prior in closed form (`maxits_ml` coupled
iterations), the PCG master's regularised map; it is written beside the wave
maps as `recvol_stateNN_<half>_flat.mrc`, the flat model under the same prior
and support. The step subset is every k-th particle of the half, 5% of it
with at least 256 and at most 1000 particles (its planes stay in memory for
the solve). The Wiener filter is undone isotropically, by the PCG operator's
own shrink factor averaged over each Fourier shell: undone voxel by voxel it
amplified the least-sampled voxels by up to `P / rho`, and the crop back to
the native box spread them over every shell (the first run's unfiltered pair
crossed an FSC of 0.143 at 17 A instead of 2 A).

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

As built in Phase 4. The cube is the native spectrum of the rotated potential
on the Nyquist sphere, `C(h,k,l)` gathered at `padf [h,k,l] R_i S_g` (with
`|l| < N/2`: the one Nyquist value a z sampled at `N` points has lies on the
sphere only at the plane's origin), times `N^3`, the shift phase and the
symmetry mean; one transform along `l` gives the sections' disk spectra
`s_j = (1/N) sum_l C exp(+2 pi i l z_j / N)`, `z_j = j - (N/2 + 1)` from the
box centre, so no real-space rotated cube is formed. Each section enters the
twice-finer grid by the disk embedding of the single-scatter chain, its
transmission is `exp(c u_j)`, the propagator acts between sections, and the
exit wave is carried back to the box centre (`N/2 - 1` sections for the beam
entering through the first section, `N/2` through the last) before
`wave_imaging`'s lens and detector, which take a given exit wave through the
same code as their own models. The wave projector runs it as the model
`multislice` under its one batch frame; its gradient passes go in chunks of
one particle per thread, each chunk's cubes scattered in particle order. At box
160 one particle costs about 0.9 s forward, 1.2 s for the adjoint sweep and
0.8 s for its cube gather on one core, and each cube scatter 0.3 s on four
threads.

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
| `src/main/volume/simple_pcg_lattice.f90` | domain | gains the gather, the coloured scatter and the rim test, moved from `reconstructor_pcg`, with the curvature `kappa` | `simple_cartesian_fourier` (the packed-window gather the moved code uses) | `reconstructor_pcg`, `wave_projector`, `wave_multislice` |
| `src/main/image/simple_wave_field.f90` (new) | domain | type `wave_field`: a 2D complex field on an `n x n` grid with its own c2c FFTW plans, forward and backward transforms, embedding of a native-disk spectrum into the `2N` grid and the band-limited crop back, both with their transposes | `simple_fftw3` | `wave_imaging`, `wave_multislice` |
| `src/main/image/simple_wave_imaging.f90` (new) | domain | type `wave_imaging`: the 2D chain of one image for projection, first Born and Rytov (exponential, lens transfer, intensity, crop), its linearisation and its adjoint; owns one thread's work fields | `wave_field`, `ctf`, `simple_math_ctf` | `wave_projector`, the simulator |
| `src/main/volume/simple_wave_projector.f90` (new) | domain | type `wave_projector extends pcg_lattice`: the volume's padded transform, the particles' cached geometry (rotation, shift, CTF parameters), one `wave_imaging` per thread, the lattice gradient accumulator; for a batch of observation planes: loss, gradient accumulation, `||J d||^2`, loss changes along a direction, and the forward image alone | `pcg_lattice`, `wave_imaging`, `ori`/`oris`, `sym`, `image` (the padded work image), `ctf` (the wavelength of the paraboloid), `simple_gridding` (the KB deapodisation), `simple_cartesian_fourier` (the centred embedding) | the strategy, the simulator |
| `src/main/volume/simple_wave_multislice.f90` (new, last phase) | domain | type `wave_multislice`: one thread's multislice of a particle from its rotated cube's spectrum (slices, recursion, reverse sweep with checkpoints), its linearisation and adjoint; the cube itself is gathered from and scattered to the lattice by `pcg_lattice`, and `wave_projector` runs it under its batch frame as the model `multislice` (Phase 4: a second batch frame would copy the first) | `wave_field`, `wave_imaging` (the lens and detection of the exit wave), `simple_fftw3` (the transform along the beam) | `wave_projector` |
| `src/main/opt/simple_nlcg_solver.f90` (new) | numerical engine | abstract `nlcg_objective` (deferred: size, loss and gradient, curvature along a direction, loss changes along it, preconditioner, dot, projection onto the model's constraints) and `nlcg_solve`, preconditioned Polak-Ribiere with the step rule of section 3.6, on the pattern of `simple_pcg_solver.f90` | none of SIMPLE's domain | the strategy; later CTF and pose refinement |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | strategy | gains the observation-batch frame of reconstruction (selection of a state's half, batched reading through `ptcl_read_session`, `prep_rec_observation` or the pass's records, native Fourier planes), extracted from the PCG worker so that the PCG worker and the wave objective share one copy | `simple_cartesian_fourier`, `simple_euclid_sigma2` (what the extracted frame calls) | PCG worker, wave strategy |
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
  are made and destroyed in the unnamed critical section, the one lock every
  FFTW planner call of SIMPLE takes (`image` takes it; FFTW's planner is not
  thread-safe, and a named section would be a second, independent lock;
  corrected in Phase 1, together with the scattering simulator, which used a
  named one). The scatter
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
simple_exec prg=reconstruct3D_wave projfile=<project> wave_model=rytov ewald=neg total_dose=<e/A^2> \
    mskdiam=<A> pgrp=<group> [pcg_mskfile=<mask>] [box_crop=<px>] [objfun=cc|euclid] [maxits=4] [nthr=<n>]
```

- **Input**: the project's `ptcl3D` field (poses, shifts, states, even/odd
  labels), its stacks and per-particle CTF (`dfx`, `dfy`, `angast`,
  `phshift`, `kv`, `cs`, `fraca`), `total_dose`, and the sigma2 of the
  particles when `objfun=euclid` (as PCG reads them: the project's canonical
  noise spectra, which an earlier refinement or reconstruction registered;
  without them the program stops). The halves come from the project's one
  maker, `ensure_field_halves`, which labels a field only when it has none.
  The support of the start and wave solves is `pcg_mskfile` when given (the
  plan's `mskfile`: no such key exists, and this one already names the PCG
  solve support), else the `mskdiam` sphere; `automsk` other than `no` is
  refused, since the density support needs a reference the program does not
  have.
- **Steps**: per state, the two flat PCG start solves, the ML prior from their
  FSC, then the wave solve of each half; the unfiltered pair's FSC and the
  half-map diagnostics; the maps written at SIMPLE's scale.
- **Project contract**: reads the project once (the builder); writes the
  maps and FSC files and registers them in the `out` segment
  (`add_vol2os_out`, `add_fsc2os_out`); one write of the `out` segment at the
  end. The particle field is not written.
- **Particle reads**: one pass per solver iteration plus the start pass, all
  through the shared observation-batch frame. SIMPLE's particle disk cache
  (`cache=yes`) is two-dimensional only; the repeated passes are served by
  the existing in-memory pass buffer instead (`fill_pass_buffer`), filled once
  per state with the state's particles at the native box, so the stacks are
  read once. No new cache.
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
`ewald` (`neg|pos|no`, default `neg`, and `no` for `projection`, which samples
the central section and refuses a curvature). Existing keys reused:
`total_dose` (required), `maxits` (default 4), `mskdiam`, `pcg_mskfile`,
`pgrp` (required, as for `reconstruct3D`), `box_crop`, `objfun`, `maxits_pcg`,
`maxits_ml`, `nthr`; `nparts`, `automsk` other than `no` and `pcg_solvent` are
refused.
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

Conventions cross the boundary through RELION star files. specter writes and
reads them with its own code; on SIMPLE's side Phase 0 found that
`export_relion` writes no poses and no origins at all (only image and
micrograph names, optics, CTF and coordinates), and that SIMPLE's star import
and its three star exports (the particle table of `simple_starfile.f90`, the
star project of `simple_starproject.f90` and the stream's export in
`simple_starproject_stream.f90`) copy RELION's `rlnOriginXAngst/YAngst`, in
angstroms, into the particle shift `x, y`, in pixels, and back, without the
pixel size. The sign is right; the units are wrong whenever the pixel size is
not 1 A. This run does not fix either (both lie outside its subsystem and its
files, and its data do not depend on them); they are reported to the
maintainer. The run's data therefore carry SIMPLE's own orientation table
beside each star file, written by the conversion below, and specter receives
SIMPLE poses through the same conversion in the oracle scripts. specter's own
star writer has a defect of its own: it stores its phase shift, in radians, in
`rlnPhaseShift`, which RELION and specter's reader take in degrees; the oracle
converts it when it rewrites specter's star files for SIMPLE. SIMPLE's star
import also needs RELION 3.1's optics block (specter writes a single unnamed
block, on which the import stops with an empty particle table) and resolves
image names against its working directory; the oracle rewrites specter's star
files accordingly.

The convention table, measured in Phase 0. A two-blob phantom off the box
centre (box 64, 1.3 A pixels) was imaged by specter at ten random RELION poses
and origins, imported with SIMPLE's star import and projected by SIMPLE's
`simulate_particles` from the imported orientations: under the conversion in
the table below the two blob positions agree within 0.019 pixel in every image
(as imported, with angstroms read as pixels: 1.5 pixels; with the opposite
sign: 9.9 pixels). The reverse direction, SIMPLE's random orientations and
shifts converted to RELION and imaged by specter, agrees within 0.016 pixel.
For the CTF, a weak noise potential (box 128, 1 A, 40 random poses, astigmatic
underfocus 1.4 and 0.9 um at 30 degrees, Cs 2.7 mm, amplitude contrast 0.1)
was imaged by specter (first Born on the Ewald surface) and by
`simulate_particles` from the imported terms: the first zero of the power
spectrum lies at 0.0611 and 0.0687 1/A along x and y for specter and at 0.0610
and 0.0687 for SIMPLE; the cross spectrum of the two image sets has the sign
of contrast inversion at 99.9% of the frequencies where `|CTF| > 0.3`, against
41% with the defocus sign flipped and 52 to 53% with the astigmatism angle
turned by 90 degrees or negated.

| Quantity | SIMPLE (`ori`, `ctf`) | RELION star | specter |
| --- | --- | --- | --- |
| Euler angles | `e1, e2, e3` in degrees, `R = Rz(e3) Ry(e2) Rz(e1)` with SIMPLE's passive `rotmat` (the transpose of specter's matrix of the same angles), the projection samples `[h,k,0] R` | `rlnAngleRot, rlnAngleTilt, rlnAnglePsi` = `e1, e2, e3` | quaternion `roma.euler_to_unitquat("ZYZ", [rot, tilt, psi])`; the template is sampled at `R r` |
| In-plane shift | `x, y` in pixels; model-side phase `exp(+2 pi i (h x + k y) / N)` | `rlnOriginXAngst, rlnOriginYAngst` = `x * smpd, y * smpd` | translation in A = the origin; the density moves by minus the translation |
| Defocus | `dfx, dfy` in um, positive underfocus (for the wave models, at the midplane) | `rlnDefocusU, rlnDefocusV` = `1e4 dfx, 1e4 dfy` in A | `dfu, dfv` in A, positive underfocus, the midplane |
| Astigmatism angle | `angast` in degrees from `+x` | `rlnDefocusAngle` = `angast` | `dfang` = `angast` (azimuth measured from `+y` with the opposite form of the cosine term, which comes to the same) |
| Phase shift | `phshift` in radians | `rlnPhaseShift` in degrees | `phaseshift` in radians |
| Spherical aberration | `cs` in mm | `rlnSphericalAberration` in mm | `cs` in A (`1e7` per mm) |
| Amplitude contrast | `fraca` | `rlnAmplitudeContrast` | `alpha` in `c = i sigma a (sqrt(1 - alpha^2) + i alpha)` |
| Voltage, pixel size | `kv`, `smpd` | `rlnVoltage`, `rlnImagePixelSize` | `voltage`, `voxel_size` |
| Transfer function | `sin(-chi + phshift + atan(fraca / sqrt(1 - fraca^2)))`; the code's `evalPhSh` is `-chi + phshift` | | wave lens `exp(-i chi)`, Kirkland's `chi = pi lambda k^2 (Cs lambda^2 k^2 / 2 - df)` minus the phase shift |
| Image polarity | particles contrast-inverted, unit background noise | | counts `D a^2 |psi|^2` (dose `D` in e/A^2, pixel size `a`, no detector), protein dark; this run's stacks `x = -(counts - D a^2) / sqrt(D a^2)`, the inverse of Ghostbuster's normalised units with contrast flip |
| Volume axes | MRC, `x` fastest, centre voxel `N/2 + 1` (1-based) | | `(Z, Y, X)`, centre voxel `N // 2` (0-based): the same voxel |
| Ewald curvature | `ewald=neg` (the key of the rulings, section 3.1) | | `ews_curvature_sign="negative"`: the beam enters at `z` index 0, the volume's first MRC section |

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

The three sets as Phase 0 generated them, under
`scratch/keep/oracle/datasets/` in the run directory (each with
`truth_potential.mrc` in volts, `particles.mrcs`, specter's star file, a RELION
3.1 star file SIMPLE's import reads, `simple_oris.txt` in SIMPLE's conventions
and `manifest.json`). Common to all: specter multislice on the twice-finer
grid, Ewald curvature negative, defocus at the midplane,
Poisson noise, amplitude contrast 0.07, Cs and voltage as listed, nothing else
(no ice, neighbours, detector MTF, coincidence loss, dose, chromatic or spatial
envelope, B-factor). The protein and nanoparticle sets propagate in a frame
padded to twice the box, so no wave wraps at the frame edge; the ribosome set
propagates on its box, as specter's reconstructor does, because padding made
it cost 7 s a particle (eight hours) against 1.8 s without; its Fresnel spread
over the box depth is at most 12 A at the band edge, inside the 24 A margin
around the particle. Each particle's potential is rendered from its rotated and
shifted atoms, so no interpolated volume enters the data, and the truth is the
same renderer's unrotated potential. Origins are uniform within 3 A (1 A for
the nanoparticle), astigmatism up to 300 A (50 A) at a uniform angle. The
stacks are contrast-inverted and normalised by the dose itself,
`x = -(counts - D a^2) / sqrt(D a^2)`, the exact inverse of specter's
"normalized" image units, so the gain `g = sqrt(D a^2)` of section 3.4 holds
exactly for these data.

| Set | Structure | Voltage, pixel, box | Particles | Dose, Cs, underfocus | Generation |
| --- | --- | --- | --- | --- | --- |
| `protein_300kV` | RuBisCO (PDB 1BXN, specter's test file), 36 288 atoms, about 120 A across (135 A at its widest) | 300 kV, 1.0 A, 192 | 6 000 | 40 e/A^2, 2.7 mm, 0.8 to 2.5 um | 1 h 44 min |
| `ribosome_100kV` | human 80S ribosome (PDB 6QZP, downloaded), 219 463 atoms, 284 to 300 A across (336 A at its widest) | 100 kV, 1.2 A, 320 | 4 000 | 40 e/A^2, 2.7 mm, 0.8 to 2.0 um | 2 h 1 min |
| `pt_np_300kV` | Pt nanoparticle EFC1_time11s of `/data/NanoX/zenodo`, 469 atoms, 2.7 nm | 300 kV, 0.358 A (as the NanoX data), 160 | 2 000 | 1 000 e/A^2, 1.2 mm, 400 to 800 A (about Scherzer) | 21 min |

The curvature-limited set in numbers: at 100 kV the wavelength is 0.0370 A,
so for the ribosome's 300 A (336 A at its widest) `sqrt(2 lambda D)` is 4.7 A
(5.0 A), and the Ewald paraboloid leaves the central plane by 1.2 Fourier
voxels at Nyquist (2.4 A). The data reach well past that limit: SIMPLE's
reconstruction with the truth poses (below) has half maps crossing a Fourier
shell correlation of 0.5 at 2.78 A and 0.143 at 2.40 A, Nyquist, so the band
from 4.7 A to 2.4 A, where the curvature acts, holds signal. The flat map
nevertheless follows the truth to a correlation of 0.95 at 4.7 A and 0.51 at
2.4 A (real-space correlation 0.81): the loss is gradual in the shell average,
because most of the mass lies near the particle's centre plane; the wave models
are to recover the periphery, which the comparison of Phase 3 measures in a
mask.

Each set was checked with SIMPLE before anything depends on it: the stack and
`simple_oris.txt` imported with `import_particles`, reconstructed by
`reconstruct3D` (gridding) with the truth poses, the half maps compared with
each other and their mean with the truth (`scratch/validate/<set>/report.json`).
The protein's map reaches Nyquist (2.0 A) both ways: the half maps cross a
Fourier shell correlation of 0.143 at 2.02 A, the map stays above 0.5 against
the truth to 2.02 A, and its real-space correlation with the truth in a sphere
is 0.78. The nanoparticle's map agrees with the truth to a Fourier shell
correlation of 0.98 to 0.998 down to 1.4 A (real-space correlation 0.70 in a
sphere), turns negative between 1.17 and 1.0 A and recovers to 0.9 between
0.94 and 0.83 A, while the half maps agree to Nyquist (0.72 A): the flat model
fails in the band of the particle's own Ewald limit, `sqrt(2 lambda D)` = 0.98 A
for its 24 A, which is what the wave models are to repair.

## 6. Risks

- **Accuracy of the KB lattice.** About `1e-3` against specter's `1e-6`. It
  sets the parity tolerance; it is far below the noise of any particle.
  Measured in Phase 2: the PCG lattice's three-tap window interpolates a
  section to 0.7 to 1% of the image contrast (1.3 to 1.8% of a gradient;
  corrected in Phase 5 to the maxima of section 8, which read 0.6 to 0.9% and
  1.2 to 1.8% here before),
  ten times the estimate; a seven-tap window brings the operator to 4e-5 to
  9e-4. Both are far below the noise of a particle; the wave operator keeps
  the PCG lattice's window, and `new_lattice` takes a window size for a
  caller that needs more.
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
| `src/main/volume/simple_reconstructor_pcg.f90` | calls the moved gather and scatter (1); its fold of the scatter accumulator moves into `pcg_lattice`, which the wave operator needs with a Hermitian symmetrisation for the two-cap planes of the curved models (2); public application of its ML prior and preconditioner and the inverse of its Wiener shrink, which the wave solve uses in volts (3) | 1, 2, 3 |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | shared observation-batch frame | 1, 3 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90` | PCG worker on the shared frame (1); its solve-support and half-map diagnostics routines made public for the wave strategy, unchanged (3) | 1, 3 |
| `src/main/image/simple_wave_field.f90` | new: 2D complex fields and plans | 1, 2, 4 |
| `src/main/opt/simple_nlcg_solver.f90` | new: nonlinear CG engine | 1, 3 |
| `src/main/image/simple_wave_imaging.f90` | new: the 2D chain, linearisation, adjoint (2); the lens and detection of a given exit wave, with their linearisation and adjoint, for multislice, and the detector step of `simulate_scattering` (4) | 2, 3, 4 |
| `src/main/volume/simple_wave_projector.f90` | new: single-scatter operator | 2, 3, 4 |
| `src/main/volume/simple_wave_projector_golden_tester.f90` | new: specter's golden images as data for the projector's tester (box 16 and 24) | 2 |
| `src/main/volume/simple_wave_multislice.f90` | new: multislice operator | 4 |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | new: program body and objective | 3, 4 |
| `src/main/commanders/simple/simple_commanders_rec.f90` | `commander_rec3D_wave` | 3 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `simulate_scattering` on the shared operator | 3, 4 |
| `src/main/image/simple_scattering_simulator.f90` | reduced in 1 and 3, deleted in 4 | 1, 3, 4 |
| `src/main/exec/simple_exec_refine3D.f90` | router case | 3 |
| `src/main/exec/simple_exec_sim.f90` | if the simulator's commander changes name | 3 |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90` | program descriptor (3); `wave_model=multislice` in its help (4) | 3, 4 |
| `src/main/ui/simple_ui_visibility_tester.f90` | the count of programs in the reconstruct3d category | 3 |
| `src/main/ui/simple/simple_ui_sim.f90` | `simulate_scattering` keys (3); its description once multislice runs on the lattice (4) | 3, 4 |
| `src/main/ui/simple_ui_params_common.f90` | shared metadata of `wave_model`, `ewald` | 3 |
| `src/main/params/simple_parameters.f90` | new keys; renamed and retired keys | 3 |
| `src/main/params/simple_parameters_parse.f90` | their registration | 3 |
| `src/main/params/simple_parameters_phases.f90` | their checks | 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration | 1, 2, 3, 4 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `suite=` help | 1, 2, 3, 4 |
| `doc/user_guide/simulate_scattering.md` | updated keys and conventions (3); multislice on the lattice (4) | 3, 4, 5 |
| `doc/user_guide/reconstruct3D_wave.md` | new user-guide page | 5 |
| `doc/user_guide/user_guide.md` | index entry | 5 |
| `doc/policies/3D/wave_model_policy.md` | new policy | 5 |
| `doc/policies/test_environment_policy.md` | rows for the new sub-suites | 5 |
| `doc/policies/3D/reconstruct3D_pcg_policy.md` | the gather, scatter, wrap test and fold it names now live in `pcg_lattice` (moved in Phase 1 and 2); a pointer to the wave-model policy | 5 |
| `doc/implementation_notes/planned/ghostbuster_scattering_simulator.md` | Ruben's note: conventions corrected, simulator on the shared operator | 5 |
| `doc/implementation_notes/planned/ghostbuster_operators_report.md` | the run's report | 5 |

## 8. Progress

### Phase 0, oracle and conventions (2026-10-23 to 2026-10-24, the Dell)

What it means: specter runs on the Dell and passes its own tests for
everything SIMPLE will reproduce; the way poses, shifts and CTF terms cross
between SIMPLE and specter is measured and written down (section 5.2); the
reference cases and the three ground-truth sets are on disk with their truth
and a manifest; and both findings about `simulate_scattering` are confirmed
with numbers (section 2), with one refinement: only its phase is inverted, its
absorption has the right sign. Two defects of SIMPLE's star interchange were
found on the way and are left for the maintainer (section 5.2). Paths below
are under the run directory `/home/elmlundho/agent_runs/wave_operators`.

Files changed in the repository: this plan only (sections 2, 5.2, 5.3, 8).

The oracle. specter `5786651` was copied from `~/src/specter` to
`scratch/keep/specter` and installed there with uv 0.13.5 (its cache and Python
under `scratch/keep/uv`; Python 3.11, PyTorch's CUDA 12.1 wheels running on the
CPU). Its seven test files for the parts SIMPLE reproduces (scattering, surface
scattering, Ghostbuster, forward-model parity, aberration functions, the
converged solver, FFTs) ran with eight workers: 196 passed, 12 skipped (eleven
need a graphics processor, one needs two), 1 failed
(`scratch/specter_tests.log`). The failure,
`test_multislice_supersample_matches_fine_grid_reference[3-negative]`, demands
bit-for-bit equality between multislice at oversampling 3 and a spelled-out
reference; on this processor they differ by 1.9e-6 on a wave of modulus 1
(single-precision rounding of the 96-point transforms), while the
oversampling-2 cases, the only ones this run uses, are bitwise equal
(`scratch/check_ms3.py`). The oracle is accepted with that recorded.

SIMPLE at the base was built in Release (`./compile_clean.sh`, gfortran 15.2.1,
system FFTW, as the maintainer's checkout is configured;
`scratch/build_release.log`). Its fast gate fails one suite at the base,
deterministically, unrelated to this run: `unit_pftc_align2D3D`, check "a class
tie draws one of the tied classes" (`simple_eul_prob_tab2D_tester.f90`, code
last changed by the base's commit `0147279c4`, "prob_beta posterior temperature
factor"). It is outside this run's files and is not fixed here; later phases
will see the same failure until the maintainer fixes it.

Exit criteria and how they were met:

- The oracle's tests pass: as above, with the one oversampling-3 rounding
  failure explained.
- The convention table is on disk and in this plan: section 5.2, measured with
  `scratch/keep/oracle/bridge_positions.py` (blob positions within 0.019 pixel
  from specter to SIMPLE, 0.016 pixel back) and `bridge_ctf.py` (first CTF
  zeros equal within 0.0001 1/A along both axes, cross-spectrum sign right at
  99.9% of frequencies, wrong defocus sign and angle detected), reports in
  `scratch/keep/oracle/bridge/`.
- The references are on disk with a manifest: `scratch/keep/oracle/references/`
  (14 small cases for the fast testers, 25 parity cases at box 64, 96 and 128,
  every model and both curvature signs, multislice also at quarter-turn
  poses), described in `scratch/keep/oracle/MANIFEST.md`. specter computed them
  in double precision; every image recomputes from its files to 8e-7 relative
  or better, and every stored gradient agrees with a central finite difference
  of specter's own images to 1.7e-5 relative or better
  (`references/verify_golden.json`, `verify_parity.json`).
- The datasets are on disk with their truth and manifests: section 5.3
  (`scratch/keep/oracle/datasets/`, generation logs `scratch/gen_datasets.log`
  and `scratch/gen_ribosome.log`). Each was checked by a SIMPLE reconstruction
  from its stack and orientation table with the truth poses
  (`scratch/validate/<set>/report.json`): the protein and ribosome maps reach
  Nyquist (half-map correlation 0.143 at 2.02 and 2.40 A) and follow the truth
  with real-space correlations of 0.78 and 0.81; the nanoparticle's map follows
  the truth to 1.4 A and fails in the band of its own Ewald limit, near 1 A.
- The two findings on `simulate_scattering` are confirmed with numbers: section
  2 (`scratch/keep/oracle/simulate_scattering_vs_specter/report.json`).

Decisions taken in this phase, with their reasons:

- The ground-truth sets render every particle from its rotated atoms rather
  than rotating one volume, so that specter's trilinear rotation does not blur
  the data; this costs 0.3 s per particle at box 192 on 24 cores.
- The stacks are normalised by the dose, not by each image's background
  statistics as specter's simulator does, so that the gain of the scale link is
  exact for these data and a later test can tell a scale error from noise.
- Each set ships SIMPLE's own orientation table beside its star files, since
  SIMPLE's star import reads origins in angstroms as pixels.
- The ribosome set propagates without FFT padding (the other two with it):
  with padding it ran at 7 s a particle, about eight hours, against 1.8 s
  without (2 h 1 min in all); specter's reconstructor propagates unpadded as
  well, and the wave's lateral spread over the box depth stays inside the
  margin around the particle.
- specter's reconstruction of the sets (section 5.3, "in the comparison") is
  left to the phase that compares: Ghostbuster reads only CryoSPARC `.cs`
  files, and on the processor a gold-standard Rytov run of the protein set
  is estimated at about ten hours.
- The comparison of `simulate_scattering` with amplitude contrast needed
  specter at a negative amplitude contrast (phase of `-V`, absorption of `+V`);
  specter refuses it in its settings check, which the comparison script
  relaxes for that variant only.
- specter's star writer stores the phase shift in radians where RELION wants
  degrees; the oracle's RELION 3.1 rewrite converts it. The datasets carry no
  phase shift; the small reference cases carry 0.3 rad and their RELION 3.1
  files say 17.19 degrees.

### Phase 1, foundations (2026-10-24, the Dell)

What it means: the pieces every later phase builds on are in place and tested,
and SIMPLE's existing reconstruction is untouched in its results. The
Kaiser-Bessel (KB) gather and scatter of planes now live in the padded lattice
type and take the Ewald curvature; the reconstruction's particle frame (which
particles, their orientations and CTF at the cropped box, their noise spectrum,
their observations as Fourier planes) is one shared copy; the wavelength and
the interaction constant have one home; the CTF exposes its phase; a 2D
complex field type and a nonlinear conjugate-gradient engine exist with their
tests. A reconstruction with the PCG backend (preconditioned conjugate
gradients) gives byte-identical maps before and after. Paths below are under
the run directory `/home/elmlundho/agent_runs/wave_operators`.

Files changed: `simple_math_ctf.f90` (`electron_wavelength`,
`electron_interaction_constant`), `simple_ctf.f90` (`eval_phase`, the former
private `evalPhSh`, public; `kV2wl` calls the shared wavelength),
`simple_rnd.f90` and its tester (`poisson_rnd`, moved from the simulator),
`simple_scattering_simulator.f90` and its tester (on the shared wavelength,
interaction constant and sampler), `simple_pcg_lattice.f90` (the plane disk
`lims2`, `sqlp`, `padsc`; `ewald_curv`, `plane_loc`, `rim_free`,
`gather_plane`, `scatter_plane` and the window helpers, moved),
`simple_reconstructor_pcg.f90` (calls them), `simple_matcher_ptcl_io.f90`
(`rec_state_half`, `rec_obs_oris`, `rec_obs_inv_sig2`, `read_rec_obs_batch`),
`simple_rec3D_pcg_strategy.f90` (the PCG worker on that frame), the new
`simple_wave_field.f90`, `simple_nlcg_solver.f90` and the testers of these and
of the lattice, `simple_commanders_test_class.f90` and
`simple_test_ui_class.f90` (three sub-suites registered), and this plan.

Exit criteria and how they were met:

- The identity obligations of section 4.4. A fixed case (400 noisy,
  CTF-modulated, shifted projections of myoglobin's potential, box 64,
  `scratch/keep/identity/`, reused by Phase 3) was reconstructed by
  `reconstruct3D rec_backend=pcg` at one thread with the base build (kept as
  `scratch/keep/build_base`) and with this phase's, for point group c1 with the
  correlation objective, d2, and c1 with `objfun=euclid` (the halves' noise
  spectra, which exercises the moved noise-spectrum path). Every map, FSC
  (Fourier shell correlation) and sigma2 file is byte-identical, and the
  diagnostics text files are identical once their wall times are removed
  (`compare_runs.sh`; `scratch/identity_final.log`). `simple_test_exec
  test=identity` was run on the same directories and records the same
  outcome; it cannot be the criterion here, since the registry entry it needs
  (`simulate_particles`, the only one whose declared patterns select this
  case's maps) also declares patterns this case does not produce and the
  diagnostics carry timings. The 6 000-particle protein set of Phase 0
  reconstructed at 24 threads also gives identical maps on both builds.
  `pcg_recon` reports all 42 gate metrics identical, and `lib_reconstruction`,
  `unit_reconstruction`, `unit_image` (with the CTF and simulator suites) and
  `unit_numerics` pass on both builds with identical numbers
  (`scratch/suites_base/`, `scratch/suites_new2/`).
- The new testers pass. `wave field` (unit_reconstruction): the forward
  transform equals a direct DFT (discrete Fourier transform) to 1e-5, round
  trips to 1e-6, the transposes of the transforms, the embedding and the crop
  to 1e-5. `PCG lattice` (unit_reconstruction): the scatter is the transpose
  of the gather on the paraboloid with point-group replication to 1e-5; the
  curved gather sits at `padf [h, k, kz] R` to 3.1e-7 relative (identity and
  random rotation); a point 5 voxels along the beam carries the paraboloid's
  phase, sign and size, to 6.1e-3 rad (bound 0.05 rad); gather and scatter are
  bitwise identical at 1, 11, 16 and 24 threads, each twice. `nonlinear
  conjugate gradient` (unit_numerics): on a quadratic the engine reproduces
  linear PCG's iterates (to 1e-4 after six iterations), with and without the
  parabola refinement, solves the system to 1e-4 and, with a support
  projection, the restricted system; on a nonconvex least squares with a
  saddle it reaches the known minimum to 1e-4 with a loss that never rises,
  absorbs steps ten times too long by halving and refuses a step it cannot
  bring downhill. `random draws` gained the Poisson sampler's mean and
  variance on both of its branches. Each sub-suite runs in 0.02 s or less.
  `scripts/check_test_registry.py` passes.
- The fast gate (`scripts/run_fast_gate.sh build`, `scratch/fast_gate_p1.log`):
  the source checks (descriptors, test registry, dependency graph) pass and 14
  of the 15 CTest entries pass, in 4.5 s; the one failure is the base's own,
  `unit_pftc_align2D3D`'s check "a class tie draws one of the tied classes",
  identical to the base build's failure recorded in Phase 0 and in code this
  run does not touch. The gate is green apart from it.

Measured on the way, for the phases that follow:

- The KB lattice's interpolation error is larger than the "about 1e-3" of
  section 6 for content off the box centre: for a point 5 padded voxels
  (a tenth of the padded box) from the centre it reaches 0.8% between lattice
  points, 0.1% near them (the first form of the curvature test compared
  complex values against 5e-3 and measured 1.03e-2; the error was the same for
  both curvature signs and followed the fractional position on the lattice, so
  the test was split into the location, exact to rounding, and the phase).
  Phase 2's parity tolerance has to be measured with this in mind.
- PCG wall time on the protein set at 24 threads, base against this phase,
  alternated five times on an idle machine: the accumulation took 6.72 s
  against 6.85 s on average (spread about 0.5 s per build), the whole run
  25.4 s against 25.7 s; no slowdown beyond the noise
  (`scratch/timing_pcg/times.txt`).

Decisions:

- The window helpers, `plane_loc` and `rim_free` are module procedures of
  `pcg_lattice`, not type-bound: `flex_pcg_t` (which extends the lattice and
  must stay unchanged) has a private type-bound `win_wraps` of its own, which a
  public binding of the same name would clash with, and module procedures
  avoid dispatch in the hot loops.
- `ewald_curv` refuses a paraboloid that would take a padded point more than
  0.49 voxel past the Nyquist sphere, so that no window reaches past the wrap
  table.
- The FFTW planner lock is the unnamed critical section, the one SIMPLE's
  `image` takes (section 4.2 named another, which would have been a second,
  independent lock); the scattering simulator is corrected to it as well.
- A review of the phase's diff by a second reader found no change to the PCG
  path at zero curvature; its findings on the solver engine (the gradient
  tolerance must use the projected gradient; a step that still raises the
  loss after backtracking must not be taken; backtracking restarts from the
  stored point) were fixed and are now tested.
- Left for Phase 5: `reconstruct3D_pcg_policy.md` names `win_wraps`, `sq_rim`
  and `scatter_plane`, which now live in `pcg_lattice`.
- Not fixed, as recorded in Phase 0: the base's fast-gate failure in
  `unit_pftc_align2D3D`.

### Phase 2, the single-scatter operator (2026-10-24, the Dell)

What it means: SIMPLE can now predict a particle image from a potential under
the projection, first Born and Rytov models, on the Ewald paraboloid of
either curvature sign, and has the exact derivative of that prediction with
respect to the potential, threaded and reproducible. Against specter it
agrees as closely as the PCG lattice's interpolation allows: within 1% of the
image contrast and 2% of a gradient with SIMPLE's three-tap
Kaiser-Bessel (KB) window, and to 4e-5 to 1e-3 with a seven-tap window, where
the operator's own numerics show. A gradient pass costs 1.7 to 3 times a PCG
accumulation of the same particles, not the ten to thirty times section 3.8
feared. Paths below are under `/home/elmlundho/agent_runs/wave_operators`.

Files changed: the new `simple_wave_imaging.f90` (the 2D chain of one
particle, its linearisation and adjoint), `simple_wave_projector.f90` (the
operator over a particle set on the PCG lattice) and their testers, the new
`simple_wave_projector_golden_tester.f90` (specter's golden images as data);
`simple_pcg_lattice.f90` (the fold of a scatter accumulator, moved from
`reconstructor_pcg` with a Hermitian option for the two-cap planes of the
curved models; an optional KB window size), `simple_reconstructor_pcg.f90`
(its fold now the lattice's, same arithmetic), `simple_wave_field.f90` (an
embedding onto a grid of the native box itself), the wave field's tester,
`simple_commanders_test_class.f90` and `simple_test_ui_class.f90` (sub-suites
`wave imaging` and `wave projector` of `unit_reconstruction`), this plan.

Exit criteria and how they were met:

- Every test of section 5.1 for these models passes (`unit_reconstruction`,
  sub-suites `wave imaging`, `wave projector`, `wave field`, `PCG lattice`;
  0.3 s together):
  - the transpose of the linearisation, at a random point of about half a
    radian, for the 2D chain and for the whole volume-to-image operator with
    point-group replication, every model and both curvature signs: 2e-8 to
    2e-6 (bound 1e-5);
  - the gradient of the misfit against a fourth-order central difference on
    potentials of about half a radian (the nonlinear regime): 8e-6 to 1.1e-4
    (bound 1e-3). Second-order differences of the single-precision chain
    are limited by rounding at small steps (scanned and recorded in the
    journal), hence the fourth-order stencil at 3e-2 of the point's scale;
  - the weak-phase identity: at zero potential the linearisation is
    `2 g sigma a N` times `forward_plane` times `build_transfer` of the PCG
    backend, the CTF with its shift, defocus and astigmatism, to 1.8e-6
    (bound 1e-5); in the 2D chain, `-2 sigma a CTF` to 4e-6 at every
    frequency but the origin, where only absorption acts (it carries no
    misfit);
  - zero curvature is the flat gather bitwise (the same routine), and Rytov
    without curvature is projection bitwise;
  - first Born is the linearisation of Rytov at zero (to 1e-6);
  - a point scatterer 4 voxels along the beam images, with `ewald=neg`, as the
    flat first Born at defocus `df + 4 A` to 8.5e-3 (bound 2e-2, the
    interpolation error Phase 1 measured for a point off the centre), and at
    `df - 4 A` 14 times worse: the size and sign of the curvature;
  - the line losses equal the misfit of the moved volume (1e-5);
  - specter's golden images (box 16 at 300 kV and box 24 at 100 kV, five
    models and three, two particles each, embedded as data): 7.3e-3 to
    1.03e-2 of the interior image contrast with the three-tap window (bound
    1.10e-2, below), 3.6e-5 to 5.6e-5 with seven taps (bound 6.1e-4).
- Parity within the tolerance measured at zero curvature. The tolerance is
  the KB lattice's interpolation error, measured on the projection cases,
  where specter's non-uniform transform and SIMPLE's lattice are both exact
  by another route: 8.2e-3 of the image contrast and 1.3e-2 of a gradient
  over the parity cases (8.77e-3 on the golden projection cases), ten times
  the 1e-3 of section 6. The bound of the golden test is 1.25 times the
  largest of these, the margin covering the spread over poses. The curved
  models sit at the same level (table below). The one excess is Rytov on the
  Pt particle (9.9e-3 images, 1.8e-2 gradients): the native-disk design
  (section 3.2, `S` zero beyond the disk) leaves out the content specter's
  native square holds beyond Nyquist (13% of the field's energy for one of
  the four particles, the transform of a box-filling strong scatterer sampled
  off the lattice), and the exponential of a 2.6-radian field mixes it back
  into the disk; inside the disk the two fields agree to 1e-5. First Born and
  projection of the same particle are unaffected (to 2e-4 with seven taps).
- Results identical across thread counts: predictions, misfits, gradients,
  curvatures, line losses, linearised planes and adjoint gradients are
  bitwise equal at 1, 11, 16 and 24 threads,
  each run twice (`wave projector`), and the lattice's curved gather and
  scatter likewise (`PCG lattice`).
- The cost table (below).
- The identity obligations still hold after the fold moved: the Phase 1 case
  (`scratch/keep/identity`, c1, d2, c1 with the noise-spectrum objective) is
  byte-identical to the base build, `pcg_recon` reports identical metrics,
  `lib_reconstruction` an identical log (`scratch/final_checks_p2.txt`,
  `scratch/identity_p2.log`). The fast gate passes apart from the base's own
  `unit_pftc_align2D3D` (Phase 0); `check_test_registry.py` passes.

The parity campaign (`scratch/keep/parity/`: `parity_driver.f90`, a program
linked against the run's build and outside the repository, `parity.py`,
`run_campaign.sh`, the tables `parity_*.tsv`). Every single-scatter reference
case of Phase 0 (golden box 16 and 24; parity box 64 myoglobin at 300 kV, box
128 bacteriorhodopsin at 100 kV, box 96 Pt particle at 300 kV; random poses,
shifts, astigmatism; eight particles, four for Pt). The potential is
band-limited by a raised cosine between 0.70 and 0.85 of the native Nyquist
radius and the stored residual likewise, so that both packages see the same
input and none where their discretisations differ (SIMPLE's native disk with
both Nyquist signs against specter's native square with one); specter
recomputes its images, exit waves and gradient of `1/2 sum r * counts` in
double precision; the errors are relative to the image contrast `I - 1`,
the scattered exit wave `psi - 1` and the gradient, on the interior
frequencies below `R - 2` (`R` the native Nyquist radius). specter's
positive curvature sign is compared at one voxel less SIMPLE defocus (section
3.1). Largest error per model and curvature sign over all cases:

| Model | KB window | image | exit wave | gradient |
| --- | --- | --- | --- | --- |
| projection | three taps (SIMPLE's) | 8.2e-3 | 7.3e-3 | 1.3e-2 |
| first Born, Ewald negative | three taps | 8.4e-3 | 7.9e-3 | 1.4e-2 |
| first Born, Ewald positive | three taps | 9.9e-3 | 9.0e-3 | 1.4e-2 |
| Rytov, Ewald negative | three taps | 9.9e-3 | 8.1e-3 | 1.8e-2 |
| Rytov, Ewald positive | three taps | 8.4e-3 | 8.8e-3 | 1.6e-2 |
| projection | seven taps | 4.9e-4 | 2.3e-4 | 5.2e-3 |
| first Born, Ewald negative | seven taps | 3.0e-4 | 1.1e-4 | 5.9e-3 |
| first Born, Ewald positive | seven taps | 9.1e-4 | 6.4e-3 | 6.0e-3 |
| Rytov, Ewald negative | seven taps | 7.8e-3 (Pt; 7.2e-4 without) | 4.7e-3 | 1.3e-2 |
| Rytov, Ewald positive | seven taps | 3.1e-3 (Pt; 6.8e-4 without) | 6.8e-3 | 9.8e-3 |

The positive sign's exit waves carry specter's one-slice offset of the
reference plane (3e-3 to 6e-3). Controls: with specter's own wavelength and
interaction constant (3.4e-4 off SIMPLE's) the images differ by 4 to 5% and
gradients by 8% on these sets, through the lens phase at high frequency; with
the curvature sign reversed by 0.2 to 0.7; without the positive-sign defocus
compensation by up to 1.2e-2. The exit-face comparison of Rytov (section
3.2): 7.6e-6 on bacteriorhodopsin with seven taps.

Oversampling 1 against 2 (seven-tap window, so that the lattice does not mask
it): on the weak golden phantoms both agree with specter to 1e-4; on the
protein and Pt sets oversampling 1 errs by 0.7 to 2.7% where oversampling 2
errs by 3e-4 to 8e-3, since the exponential and the intensity double the
bandwidth. It is 1.8 times faster at box 256. The oversampling stays 2 and
no parameter is offered (section 3.2); `osf` remains an argument of
`wave_projector%new` for such measurements, and its transpose is tested.

Cost per particle on 24 cores (the Release build, the machine otherwise idle,
two runs within 5%; `scratch/keep/parity/timing_driver.f90`, `timing.log`):

| Box | Oversampling | PCG accumulation | predict (projection, first Born, Rytov) | gradient | gradient / PCG |
| --- | --- | --- | --- | --- | --- |
| 128 | 1 | 0.31 ms | 0.17-0.19 ms | 0.42-0.45 ms | 1.3-1.4 |
| 128 | 2 | 0.32 ms | 0.45-0.58 ms | 0.86-0.94 ms | 2.7-3.0 |
| 256 | 1 | 2.25 ms | 1.0-1.1 ms | 2.2 ms | 1.0 |
| 256 | 2 | 2.28 ms | 2.0-2.2 ms | 3.7-4.0 ms | 1.7-1.8 |

Decisions and findings, with their reasons:

- Rytov exponentiates at the exit face (section 3.2): specter's model and the
  definition of the approximation; the midplane differed by 2 to 17%.
- The phase shift leaves the unscattered beam alone, `T(0) = 1`: applied to
  every frequency it would be a common phase and vanish from the intensity
  (found by the weak-phase test, 0.2 off before).
- SIMPLE's wavelength stays the CTF's (the plan's one home, bitwise
  identical to the CTF); it is 3.4e-4 short of the exact relativistic value,
  which moves the lens phase by about 0.1 rad at Nyquist for 2 um of
  defocus at 300 kV. This concerns SIMPLE's CTF as a whole and is left for
  the maintainer.
- The wave operator keeps the PCG lattice's three-tap window: its error is
  far below a particle's noise and it matches the PCG start map's
  interpolation; `new_lattice` takes an optional window size, used here to
  separate the lattice's error from the operator's and tested by the golden
  cases at seven taps.
- specter's positive curvature sign is one voxel of defocus off the
  symmetric reference (section 3.1); SIMPLE keeps the symmetric one.
- The fold moved into the lattice is passed the extension's own work image
  with the parent component as the object (`self%pcg_lattice%...`), so the
  two arguments are disjoint; `flex_pcg_t` already binds `fold_accum`, hence
  the name `fold_lattice_accum`.
- The finite-difference tests needed a fourth-order stencil on half-radian
  potentials: in single precision a central difference on weak potentials is
  limited by rounding (the first form failed by up to 3e-2 for that reason;
  the derivative was right, as the adjoint identities and the step scan
  showed).
- The point-scatterer test first used a 0.05 V point, whose intensity change
  (6e-5) is below single-precision rounding of the intensity; it uses 50 V.
- A review of the phase by a second reader asked for: the golden data module
  named as a tester, so that builds without tests leave it out; specific
  module imports instead of the umbrella module, and the model constants
  imported from `simple_wave_imaging` alone; the lattice's gather not
  opening a nested team inside the operator's particle loop, and the
  operator refusing to run inside a caller's parallel region; input checks
  (phase-flipped particles, amplitude contrast, batch and particle ranges,
  an even box) before the parallel regions; the batch planes passed without
  copies; the lens geometry tabulated once per object. All were made; the
  tests and the parity tables were rerun, the latter moving by 1e-6 at most.

### Phase 3, the program (2026-10-24, the Dell)

What it means: SIMPLE now has the program `reconstruct3D_wave`, which
reconstructs half maps from a project with known poses and contrast transfer
function (CTF) under the projection, first Born or Rytov model, on either
curvature sign of the Ewald sphere, and `simulate_scattering` draws its images
with the same operator and corrected conventions. On the curvature-limited
ribosome at 100 kV, Rytov with the data's curvature sign raises the agreement
of the map with the truth beyond the curvature limit (the Fourier shell
correlation, FSC, of the map against the truth between 3 and 2.4 A) from 0.66
for the flat model to 0.81, and the wrong sign drops it to 0.21; specter's own
converged reconstruction reaches 0.83 there and agrees with SIMPLE's in every
other band. On the ordinary protein at 300 kV all models agree, Rytov slightly
ahead (0.762 against 0.733 between 3 and 2 A) and level with specter (0.765).
Paths below are under `/home/elmlundho/agent_runs/wave_operators`.

Files changed: the new `simple_rec3D_wave_strategy.f90` (the program body and
the objective the nonlinear conjugate-gradient engine minimises) and its
library tester; `simple_rec3D_pcg_strategy.f90` (the accumulation of a
state's half, extracted from the PCG worker so that the worker and the wave
program share it; its solve-support and half-map diagnostics routines made
public); `simple_reconstructor_pcg.f90` (the Wiener shrink factor shared by
the closed form and its isotropic inverse; public application of the ML
prior and of the preconditioner); `simple_wave_projector.f90` and its tester
(the terms of the gain fit; new tests below); `simple_commanders_rec.f90`,
`simple_exec_refine3D.f90`, `simple_ui_reconstruct3D.f90`,
`simple_ui_visibility_tester.f90` (the program); `simple_parameters.f90`,
`simple_parameters_parse.f90`, `simple_parameters_phases.f90`,
`simple_ui_params_common.f90` (keys `wave_model` and `ewald`;
`scattering_model` and `absorption` removed, no alias);
`simple_commanders_sim.f90`, `simple_ui_sim.f90`,
`simple_scattering_simulator.f90` and its tester, `doc/user_guide/simulate_scattering.md`
(`simulate_scattering` on the shared operator); `simple_commanders_test_class.f90`
and `simple_test_ui_class.f90` (the library sub-suite); this plan.

Exit criteria and how they were met:

- The program runs on both datasets: every run of
  `scratch/keep/gt3/run_gt3.sh` ended normally (logs and `runs.txt` in
  `scratch/gt3/<set>/gt3/`), with objective `cc` (the sets carry no noise
  spectra) and the same even/odd halves as specter's (`scratch/keep/gt3/make_cs.py`
  writes specter's CryoSPARC file and SIMPLE's orientation table from the
  set's truth). On the ribosome (box 320, 4 000
  particles) at 24 threads with nothing else running, twice each and
  alternated (`scratch/keep/gt3/timing.sh`, `timing.txt`): `reconstruct3D`
  with PCG 79 and 90 s, `reconstruct3D_wave` (Rytov, four iterations per
  half) 554 and 577 s, about 6.5 times; of the latter, the flat start takes
  about 100 s, each half's solve about 220 s. The protein's wave runs took
  4 to 6.5 min at 12 threads beside specter.
- The comparisons of section 5.3 are recorded: gridding and flat PCG from
  `reconstruct3D`, the wave program's flat start (the regularised PCG map
  under the same prior and support), first Born and Rytov at both curvature
  signs and without curvature, and specter's gold-standard converged
  reconstruction (`specter reconstruct particle`, Rytov, negative sign,
  default solver: five Adam epochs and four converged iterations; protein
  41 min, ribosome 2 h 16 min at 12 threads), scored against the truth
  potential by `scratch/keep/gt3/score.py` (`scores_<set>.json`): the FSC
  of the half-map mean against the truth averaged over resolution bands
  (insensitive to a map's filtering), the real-space correlation in a soft
  sphere of the particle's extent and in its outer half (the periphery, where
  the curvature acts most), the half-map FSC. The wave rows score the
  unfiltered pair for the real-space correlations (the regularised maps are
  Wiener-filtered; the flat start has no unfiltered twin, hence its lower
  real-space values).

  The protein (`protein_300kV`, curvature limit 2.17 A, Nyquist 2.0 A):

| Map | truth FSC 20-8A | truth FSC 8-5A | truth FSC 5-4A | truth FSC 4-3A | truth FSC 3-2A | real-space corr. (sphere) | (periphery) | half-map FSC 0.143 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gridding (reconstruct3D) | 0.9929 | 0.9874 | 0.9879 | 0.9571 | 0.7353 | 0.808 | 0.808 | 2.02 A |
| flat PCG (reconstruct3D) | 0.9925 | 0.9869 | 0.9883 | 0.9572 | 0.7297 | 0.806 | 0.805 | 2.02 A |
| flat start (same prior) | 0.9925 | 0.9869 | 0.9883 | 0.9574 | 0.7331 | 0.755 | 0.755 | 2.02 A |
| first Born, ewald=neg | 0.9913 | 0.9897 | 0.9883 | 0.9575 | 0.7465 | 0.813 | 0.814 | 2.02 A |
| first Born, ewald=pos | 0.9913 | 0.9897 | 0.9883 | 0.9552 | 0.7077 | 0.801 | 0.802 | 2.02 A |
| first Born, ewald=no | 0.9912 | 0.9896 | 0.9883 | 0.9568 | 0.7340 | 0.809 | 0.810 | 2.02 A |
| Rytov, ewald=neg | 0.9988 | 0.9919 | 0.9877 | 0.9607 | 0.7623 | 0.817 | 0.816 | 2.02 A |
| Rytov, ewald=pos | 0.9988 | 0.9919 | 0.9875 | 0.9574 | 0.7072 | 0.801 | 0.798 | 2.02 A |
| projection (Rytov, ewald=no) | 0.9988 | 0.9919 | 0.9877 | 0.9601 | 0.7486 | 0.812 | 0.811 | 2.02 A |
| specter, gold converged, unfiltered | 0.9988 | 0.9922 | 0.9878 | 0.9606 | 0.7648 | 0.816 | 0.815 | 2.02 A |
| specter, gold converged, filtered | 0.9988 | 0.9922 | 0.9878 | 0.9607 | 0.7642 | 0.788 | 0.788 | 2.02 A |

  The ribosome (`ribosome_100kV`, curvature limit 4.7 A, Nyquist 2.4 A):

| Map | truth FSC 20-8A | truth FSC 8-5A | truth FSC 5-4A | truth FSC 4-3A | truth FSC 3-2.4A | real-space corr. (sphere) | (periphery) | half-map FSC 0.143 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gridding (reconstruct3D) | 0.9751 | 0.9575 | 0.9485 | 0.9042 | 0.6658 | 0.831 | 0.822 | 2.42 A |
| flat PCG (reconstruct3D) | 0.9743 | 0.9572 | 0.9485 | 0.9015 | 0.6504 | 0.828 | 0.821 | 2.42 A |
| flat start (same prior) | 0.9743 | 0.9573 | 0.9485 | 0.9023 | 0.6596 | 0.808 | 0.804 | 2.42 A |
| first Born, ewald=neg | 0.9795 | 0.9660 | 0.9591 | 0.9277 | 0.7789 | 0.845 | 0.845 | 2.42 A |
| first Born, ewald=pos | 0.9799 | 0.9651 | 0.9453 | 0.8199 | 0.3584 | 0.774 | 0.755 | 2.42 A |
| first Born, ewald=no | 0.9783 | 0.9640 | 0.9556 | 0.9094 | 0.6669 | 0.825 | 0.821 | 2.42 A |
| Rytov, ewald=neg | 0.9984 | 0.9927 | 0.9868 | 0.9559 | 0.8100 | 0.876 | 0.874 | 2.42 A |
| Rytov, ewald=pos | 0.9984 | 0.9901 | 0.9653 | 0.7947 | 0.2063 | 0.775 | 0.743 | 2.42 A |
| projection (Rytov, ewald=no) | 0.9985 | 0.9926 | 0.9851 | 0.9411 | 0.7066 | 0.860 | 0.852 | 2.42 A |
| specter, gold converged, unfiltered | 0.9982 | 0.9922 | 0.9862 | 0.9559 | 0.8328 | 0.879 | 0.877 | 2.42 A |
| specter, gold converged, filtered | 0.9982 | 0.9922 | 0.9862 | 0.9559 | 0.8328 | 0.878 | 0.877 | 2.42 A |

  The gain the data prefer (section 3.4), against the gain from the dose:
  protein 6.19 to 6.22 (Rytov), 6.16 to 6.18 (projection) against 6.32;
  ribosome 7.16 to 7.20 (Rytov), 7.03 to 7.06 (projection) against 7.59.
  First Born's fit is not a gain (section 3.4). The shortfall is real and
  found: SIMPLE normalises each particle to the spread of its background
  outside the mask, and in these stacks that background carries the
  particle's delocalised signal (at 2 um underfocus): its spread is 1.098 of
  the true noise for the ribosome and 1.035 for the protein
  (`scratch/keep/gt3/bg_sdev.py`), so the particles are 9% and 3.5% smaller
  than the dose implies. With the dose corrected for it (gain 6.91 for the
  ribosome) the fitted gain agrees to 1% (6.84, 6.87) and the map's scale in
  volts (0.791 of the truth by least squares in the sphere) equals specter's
  (0.797). The scale link of real data therefore needs a noise estimate free
  of signal; the fitted gain shows when it is off.
- The curvature-limited set shows the expected ordering: beyond the
  curvature limit Rytov (0.810) and first Born (0.779) with the data's sign
  beat the flat maps (0.650 to 0.666) and lose with the wrong sign (0.206,
  0.358); without curvature Rytov is projection (0.707), which gains from the
  nonlinear term alone; first Born, which lacks it, stays flat (0.667). In
  the bands below the limit the nonlinear models gain as well (20 to 8 A:
  0.998 against 0.975), because the ribosome's phase at 100 kV is no longer
  weak.

Tests added or rewritten in this phase (each passes; `scratch/tests_p3/`,
`scratch/suites_p3/`):

- `scattering simulator` (`unit_image`), rewritten to be sign- and
  defocus-sensitive: a weak random sheet of potential in the midplane images
  as SIMPLE's weak-phase model `-2 sigma a CTF (*) u` to 2.1e-3 of its
  contrast (bound 5e-3, twice the expected second order and rounding); the
  same sheet 16 A downstream images at defocus `df + 16 A` to 2.0e-3 and at
  `df - 16 A` with an error of 0.27; with the beam reversed (`ewald=pos`) at
  `df - 16 A` to 2.5e-3; the earlier simulator's conjugated lens and
  propagator, or its exit-face defocus, fail these. The analytic limits and
  the detector's Poisson counts are kept.
- `wave projector` (`unit_reconstruction`): the scale link of the program, the
  Gauss-Newton curvature at zero against `kappa^2 / boxpd^3` times PCG's
  normal operator, to 2.1e-6 (bound 1e-4); the vacuum, relative intensity 1
  at zero potential, every model; multislice of a weak potential against
  first Born on the paraboloid of the same sign, 1.3e-3 and 1.8e-3 at the
  identity pose (bound 1e-2), and through `rotvol` at a general pose a
  correlation of 0.9994 against 0.84 and 0.74 for the transposed pose; the
  gain terms join the thread-identity test.
- `wave reconstruction` (`lib_reconstruction`, the library tier of section
  5.1): a phantom of 150 Gaussian atoms 100 A across, box 128 at 1 A, 300
  particles imaged by multislice at 100 kV with detector noise; between the
  curvature limit (2.7 A) and Nyquist the truth FSC is 0.775 for the flat
  start, 0.814 for Rytov with the data's sign and 0.578 for the other. It
  takes 189 s at 8 threads on the idle machine, most of it the multislice
  simulation (`rotvol`'s gather over the padded lattice, about 0.4 s a
  particle), which the multislice phase replaces.

Other checks: the PCG identity case of Phase 1 is byte-identical on this
phase's build, after the extraction of the accumulation frame
(`scratch/identity_p3.log`: c1, d2 and c1 with the euclid objective, 0
different); `reconstruct3D_wave` gives byte-identical maps and FSC files at
1, 11, 16 and 24 threads, each twice (`scratch/keep/gt3/thread_identity.sh`,
`scratch/thread_identity_p3.log`); `unit_image` (934 assertions),
`unit_reconstruction` (418), `unit_numerics` (1365), `lib_reconstruction` (86)
and `pcg_recon` pass; the fast gate passes 14 of 15 entries, the one failure
the base's own `unit_pftc_align2D3D` (`scratch/fast_gate_p3.log`);
`check_test_registry.py` and `check_descr.py` pass.

Decisions taken in this phase, with their reasons:

- specter reconstructs with `rotate_mode=surface` (its exact non-uniform
  sampling of the Ewald paraboloid, the model SIMPLE's operator implements):
  its default real-space rotation costs 0.51 s a particle against 0.012 s at
  box 192 (`scratch/keep/gt3/time_specter.py`), about 8 h for the protein and
  33 h for the ribosome, where its converged solver's batch of 128 would also
  exceed the Dell's memory.
- The program's design decisions of sections 3.4, 3.6 and 4.3 (measured and
  decided here): the prior in volts carries `boxpd^3`; the Wiener filter is
  undone per shell; the flat start is written beside the wave maps; the
  repeated particle passes use the pass buffer (the particle cache is
  two-dimensional only); `pcg_mskfile` is the support key (no `mskfile` key
  exists); halves from `ensure_field_halves`; `objfun=euclid` consumes the
  project's canonical noise spectra.
- `winsz` (an existing key, default the lattice's three taps) sets the wave
  operator's window: measured on the ribosome, seven taps change no band
  (0.8102 against 0.8101 between 3 and 2.4 A) at 1.4 times the run time
  (745 s against 549 s, 24 threads), so the default stays; the key serves a caller who wants the operator's numerics
  alone.
- `ewald` defaults to `no` for `wave_model=projection` (the central section by
  definition) and to `neg` otherwise; projection with a curvature is refused
  (both programs).
- `simulate_scattering`'s multislice keeps the rotated cube of `rotvol` until
  the multislice phase builds it from the lattice, now from one padded
  transform per run (byte-identical images, half the work). `rotvol`'s
  Fourier gridding does not correct the Kaiser-Bessel window's roll-off: a
  weak potential's multislice image through it is 8 to 10% weaker off the
  box centre than without it (correlation 0.9995). Its lens is SIMPLE's CTF
  with the unscattered beam untouched, its propagator Kirkland's; its midplane
  is section `N/2 + 1`, `N/2 - 1` sections upstream of the exit with the beam
  entering through the first section and `N/2` through the last (the even
  box's asymmetry, found by the reversed-beam test); particle shifts are
  refused there. The single-scatter models apply the particles' shifts and
  phase shifts; the output orientation table records the CTF parameters used.
- On the native grid the astigmatic lens phase of the Nyquist row and column
  is ambiguous (`+N/2` and `-N/2` are one frequency); tests band-limit to
  the disk.

Open, for the maintainer:

- The last band of the ribosome stays 0.02 below specter's (0.810 against
  0.833 between 3 and 2.4 A; 0.733 against 0.786 at 2.43 A) although the
  half-map FSCs agree (0.686 against 0.681 at 2.43 A): an error common to both
  halves. Not the window (seven taps), not convergence (ten iterations give
  0.8104), not the scale link (corrected gain, 0.8097). What remains lies in
  the observation preparation the program shares with PCG, which this phase
  may not change: the particles' edge taper, at a defocus whose
  high-resolution delocalisation (about 220 A at 2.5 A for 1.5 um at 100 kV)
  exceeds the box margin, and the background-based noise normalisation.
- The noise normalisation's background spread (1.035 to 1.098 here) breaks the
  absolute scale link of the nonlinear models whenever the box leaves little
  signal-free margin; the fitted gain detects it.

### Phase 4, multislice (2026-10-24, the Dell)

What it means: SIMPLE's wave operator now has multislice, the model of
multiple scattering, with its exact derivative, and `reconstruct3D_wave` and
`simulate_scattering` both use it; the old scattering simulator is gone. Its
images agree with specter's where specter's own rotation is exact (to 7e-5
for bacteriorhodopsin, 2e-3 for myoglobin), its derivative with specter's
once specter's rotation artefact is removed (2.5 to 5%), and it tends to
first Born on the Ewald paraboloid for weak potentials. On the 469-atom Pt
nanoparticle, the nonlinear models (Rytov and multislice alike) repair part
of the band where the linear map fails, around the particle's own Ewald limit
near 1 A, and place the atoms more accurately; every map, the linear ones
included, finds all 469 atoms at this dose. specter's own multislice
reconstruction matches SIMPLE's over the last band; its converged Rytov
repairs more of the band near the limit, and the reason is found: the
program's prior, built from the flat maps' agreement, is strongest exactly
there.
Paths below are under `/home/elmlundho/agent_runs/wave_operators`.

Files changed: the new `simple_wave_multislice.f90` (the multislice chain of
one particle: the rotated spectrum cube to its sections by one transform along
the beam, transmission `exp(c v_j)` and Kirkland's propagator on the twice-finer
grid, the lens of `wave_imaging`, and its exact derivative and adjoint, a
reverse sweep with checkpoints every `ceil(sqrt N)` sections) and its tester;
`simple_pcg_lattice.f90` and its tester (`gather_cube` and `scatter_cube`, the
particle's rotated cube read from and written to the padded lattice, the second
the exact transpose of the first); `simple_wave_imaging.f90` and its tester (the
model constant `WAVE_MSLICE`; lens and detection entry points for a given exit
wave, sharing the arithmetic of the two-dimensional chain; the detector's
Poisson and Gaussian noise, moved here from the simulator);
`simple_wave_projector.f90` and its tester (multislice as the fourth model: one
chain per thread, gradient passes in chunks of the thread count whose cubes are
scattered in particle order); `simple_rec3D_wave_strategy.f90` and its library
tester (`wave_model=multislice`; the tester's data are now made by the shared
operator's multislice); `simple_commanders_sim.f90`, `simple_ui_sim.f90`,
`doc/user_guide/simulate_scattering.md` (`simulate_scattering`'s multislice on
the shared operator; particle shifts apply to every model);
`simple_ui_reconstruct3D.f90` (help text); `simple_scattering_simulator.f90`
and its tester deleted; `simple_commanders_test_class.f90` and
`simple_test_ui_class.f90` (suite `wave multislice` in `unit_reconstruction`;
`scattering simulator` removed from `unit_image`); this plan (file table;
the row of `wave_multislice` in 4.1; "As built in Phase 4" in 3.7).

Exit criteria and how they were met:

- `wave_multislice` and its tests pass (`scratch/suites_p4/`, the suite
  `wave multislice`):
  - the adjoint of the sweep: `<J x, y> = <x, J* y>` on random cubes to 2.9e-7
    and 8.0e-7 for the two beam directions (bound 1e-5);
  - finite differences: the adjoint's gradient of a misfit against a
    fourth-order difference to 6.7e-4 and 4.3e-4 (bound 1e-3);
  - the weak limit: a weak sheet of potential (peak phase 5e-3 rad) in the
    midplane images as SIMPLE's weak-phase model to 1.2e-3 of its contrast; 16 A
    downstream it images at defocus `df + 16 A` to 1.5e-3 and at `df - 16 A` with
    an error of 0.12; with the beam reversed at `df - 16 A` to 2.5e-3 (bound
    5e-3, the expected second order); zero and constant potentials give the
    vacuum and the analytic absorption `exp(-2 sigma a N alpha V0)`;
  - in `wave projector` (`unit_reconstruction`): the projector's multislice
    against first Born on the paraboloid of the same sign for a weak potential,
    7.4e-4 and 9.8e-4 (bound 2.2e-2); a smooth off-centre blob at 100 kV, where the
    two curvature signs differ: 2.3e-3 against the same sign and 2.1e-2 against
    the opposite one (criterion: the opposite at least four times the same; the
    signs differ by only about 2% at box 24, so a flip swaps the errors); the
    adjoint and gradient tests cover multislice at both signs, and the vacuum
    and thread-identity tests cover it as well;
  - the lattice's tester: the cube gather's adjoint to 1e-5, its central plane
    bit-identical to the flat section, thread identity.
- specter parity (`scratch/keep/parity/`, `parity_osf2_simpleconst_only_multislice.tsv`,
  `_w3_` with seven taps, `run_multislice*.log`). specter runs multislice only
  with its real-space (trilinear) rotation. At quarter-turn poses, where that
  rotation is exact, the images agree to 6.5e-5 (bacteriorhodopsin) and 1.8e-3
  (myoglobin, equal to specter's own difference between its two rotation modes
  there, 1.9e-3), the exit waves to 7e-5 and 1.2e-3. At general poses they
  differ by 9 to 22%, which is specter's interpolation: its two rotation modes
  differ from each other by 12.6% (myoglobin), 12.7% (bacteriorhodopsin) and 19%
  (the Pt particle) on first Born (`specter_rotation_check.py`, `.log`). The
  gradients first differed by 30%: specter's real-space-rotation gradient
  differs from its own surface-mode gradient by 31.5% on first Born and
  decorrelates beyond 0.75 of Nyquist, while SIMPLE's gradients are exact
  (white-direction finite differences: first Born to 6e-6, multislice 4.6e-4;
  `scratch/p4_debug/grad_dbg.f90`). Corrected by specter's real-minus-surface
  difference (`multislice_gradient_check.py`, `.log`): 2.5e-2 for myoglobin
  (1.4e-2 below 0.7 of Nyquist) and 5.2e-2 for bacteriorhodopsin (3.1e-2),
  against multislice-versus-first-Born differences of 6% and 12%.
- `wave_model=multislice` in the program: `reconstruct3D_wave` runs it on the
  Pt set (below); its maps are byte-identical at 1, 11, 16 and 24 threads, twice
  each (300 Pt particles cropped to box 64, 8 files per run, 0 differences;
  `scratch/keep/gt4/thread_identity_ms.sh`, `scratch/thread_identity_p4.log`).
  Cost at box 160: about 0.9 s forward, 1.2 s adjoint and 0.8 s gather per
  particle on one thread; 90 min for the 2 000-particle set at 12 threads (Rytov
  85 s).
- simulate_scattering's multislice is the shared operator's, and
  `simple_scattering_simulator.f90` is deleted: the four models run through
  `wave_projector`; its tests moved to the operator's suites (above); the
  library test `wave reconstruction` now draws its multislice data through it
  and passes (truth FSC beyond the curvature limit: flat start 0.798, Rytov
  with the data's sign 0.840, the other sign 0.599; 291 s with specter running
  beside it).
- Ground truth on the nanoparticle set: below.

The nanoparticle (`pt_np_300kV`: a 469-atom Pt particle, box 160 at 0.358 A,
2 000 particles, 1 000 per half, at 300 kV; its own Ewald limit `sqrt(2 lambda D)` about 0.98 A;
scratch/gt3/pt_np_300kV/gt3, `scratch/keep/gt4/run_gt4.sh`, mask diameter 40 A,
12 threads beside specter). Truth FSC of the half-map mean near the band where
the flat model fails, over the last band, the real-space correlation, and the
atoms: SINGLE's atom detection (`single_exec prg=detect_atoms smpd=0.358
element=Pt mskdiam=34 atom_thres=no`, the same for every map), scored against
the truth coordinates one-to-one within 1.0 A (`scratch/keep/gt4/atoms.sh`,
`score_atoms.py`, `atoms.txt`; on the truth map itself 469 of 469, root mean
square error 0.024 A). Wave rows score the unfiltered pair; atoms are detected
on the program's merged map. The table is `scratch/keep/gt4/table_pt.md`
(`summarize_pt.py`); the coarse bands of section 5.3 are in `table_pt_bands.md`
(`scratch/keep/gt3/summarize.py pt_np_300kV`).

| Map | 1.27 A | 1.22 A | 1.17 A | 1.12 A | 1.08 A | 1.04 A | 1.00 A | 3 to 0.72 A | real-space corr. | atoms found (of 469), precision | rms position error |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gridding (reconstruct3D) | 0.835 | 0.626 | 0.033 | -0.087 | -0.053 | 0.048 | 0.366 | 0.750 | 0.735 | 469, 1.000 | 0.035 A |
| flat PCG (reconstruct3D) | 0.832 | 0.621 | 0.017 | -0.089 | -0.056 | 0.046 | 0.361 | 0.746 | 0.732 | 469, 1.000 | 0.036 A |
| flat start (same prior) | 0.839 | 0.642 | 0.002 | -0.249 | -0.158 | -0.049 | 0.404 | 0.732 | 0.703 | - | - |
| first Born, ewald=neg | 0.783 | 0.418 | 0.091 | 0.041 | -0.009 | 0.002 | 0.214 | 0.747 | 0.703 | 469, 1.000 | 0.029 A |
| Rytov, ewald=neg | 0.933 | 0.698 | 0.308 | 0.202 | 0.068 | 0.179 | 0.529 | 0.805 | 0.753 | 469, 1.000 | 0.029 A |
| multislice, ewald=neg | 0.937 | 0.730 | 0.380 | 0.182 | 0.068 | 0.164 | 0.498 | 0.801 | 0.749 | 469, 1.000 | 0.030 A |
| multislice, ewald=pos | 0.922 | 0.694 | 0.367 | 0.176 | 0.066 | 0.149 | 0.426 | 0.734 | 0.721 | 469, 1.000 | 0.030 A |
| Rytov, prior 10x weaker (tau=10) | 0.943 | 0.774 | 0.457 | 0.495 | 0.235 | 0.341 | 0.642 | 0.836 | 0.760 | 469, 1.000 | 0.029 A |
| multislice, prior 10x weaker (tau=10) | 0.945 | 0.802 | 0.581 | 0.469 | 0.223 | 0.314 | 0.624 | 0.834 | 0.756 | 469, 1.000 | 0.029 A |
| specter Rytov, gold converged | 0.946 | 0.793 | 0.502 | 0.626 | 0.397 | 0.441 | 0.677 | 0.859 | 0.757 | 469, 1.000 | 0.027 A |
| specter multislice, Adam epochs only | -0.612 | -0.899 | -0.974 | -0.956 | -0.935 | -0.946 | -0.937 | 0.089 | 0.473 | 444, 1.000 | 0.099 A |
| specter multislice, gold converged | 0.851 | 0.677 | 0.260 | 0.315 | 0.144 | 0.194 | 0.486 | 0.801 | 0.760 | 469, 1.000 | 0.027 A |

Fitted gains against 11.32 from the dose: multislice 10.29, 10.39 (`neg`),
10.19, 10.31 (`pos`); Rytov 10.75, 10.83 (the shortfall of Phase 3's
background normalisation).

What the table shows: every map finds all 469 atoms with no false ones; the
nonlinear models place them better (0.029 to 0.030 A against 0.035 to
0.036 A for the flat maps), and Rytov and multislice with the data's sign
repair part of the band around 1.17 A (0.31 and 0.38 against 0.03 flat), with
the wrong sign only in the last band (0.734 against 0.801). The band from 1.12
to 1.04 A stays nearly empty in the program (0.07 to 0.20) while specter's
converged Rytov fills it (0.40 to 0.63). The cause is found and is the program's
prior, not the operator: the maximum-likelihood prior is built from the flat
start pair's FSC, which collapses exactly where the flat model fails, so the
prior is strongest there; more iterations do not help (12 give the same as 4:
the steps stall at the solver's minimum fraction). With the prior ten times
weaker (`tau=10`) Rytov reaches 0.50, 0.24, 0.34 there and multislice 0.47,
0.22, 0.31 (last band 0.836 and 0.834 against specter's 0.859); specter's
prior comes from its own wave-model maps. Multislice and Rytov agree closely on
this particle at 300 kV.

specter's multislice (`scratch/keep/gt4/run_specter.sh pt_np_300kV 16
multislice real`, log `specter_pt_multislice.log`; its only rotation for
multislice is the real-space one, 8 h 28 min at 16 threads) agrees with
SIMPLE's after its converged stage: the last band 0.801 for both, atoms 469 of
469 at 0.027 A against 0.030 A; specter is ahead from 1.12 to 1.04 A (0.32,
0.14, 0.19 against 0.18, 0.07, 0.16) and behind from 1.27 to 1.17 A (0.85,
0.68, 0.26 against 0.94, 0.73, 0.38). Its Adam stage alone (five epochs) is no reference: those
maps anticorrelate with the truth beyond 1.3 A (FSC near -0.95) while their
halves agree to 0.74 A, with the atoms unshifted (mean offset below 0.003 A,
`scratch/p4_debug/atom_shift.py`); the converged solver removes it. Scored half
against half (each map alone), SIMPLE's multislice and specter's give the same
last band (0.767 both).

Other checks: the PCG identity case of Phase 1 is byte-identical on this
phase's build (`scratch/identity_p4.log`: c1, d2, c1 with the euclid
objective, 0 different); `unit_image` (925 assertions), `unit_reconstruction`
(440), `unit_numerics` (1365), `lib_reconstruction` (86) and `pcg_recon` pass;
the fast gate passes 14 of 15 entries, the one failure the base's own
`unit_pftc_align2D3D` (`scratch/fast_gate_p4.log`); `check_test_registry.py`
and `check_descr.py` pass (`scratch/final_checks_p4.txt`,
`scratch/check_descr_p4.log`); `git diff --summary` shows no mode changes.

Decisions taken in this phase, with their reasons:

- The chain is a per-thread type in `simple_wave_multislice.f90`, run by
  `wave_projector` as its fourth model, not a second extension of the lattice
  with its own batch interface (section 4.1 said so): the batch frame, the
  scale link, the gain terms and the thread-ordered scatter are the projector's,
  and a second frame would copy them. The cube gather and scatter belong to the
  lattice.
- The constant is `WAVE_MSLICE`: Fortran is case-insensitive, and
  `WAVE_MULTISLICE` would collide with the type `wave_multislice`.
- Multislice requires a curvature sign (`ewald=neg` or `pos`; both programs
  refuse `no`): the beam direction is the sign.
- The orientation test's criterion is a factor of four between the opposite and
  the same sign (my first factor of ten was not attainable at box 24, where the
  two signs differ by about 2%); recorded before it was relied on.
- specter's multislice parity uses its real-space rotation, its only mode for
  multislice; parity at general poses is therefore bounded by specter's
  interpolation, measured separately.
- The ML prior is left as the plan fixes it (from the start pair) because
  changing it would change Phase 3's results; see Open.

Open, for the maintainer:

- The prior near the Ewald limit: `tau=10` shows that a weaker prior there
  recovers most of the band on the nanoparticle. Rebuilding it from the wave
  pair's FSC after a solve under the present prior would not: that pair's
  half-map FSC from 1.12 to 1.04 A is 0.00 to 0.05 (the prior holds it down),
  against 0.13 to 0.42 under `tau=10` and 0.29 to 0.76 for specter's Rytov
  (whose prior comes from its own wave-model maps). A prior from the pair of a
  weakly regularised first wave solve is the candidate; a later phase or the
  maintainer decides.
- For Phase 5 (its files): the index text in `doc/user_guide/user_guide.md`
  and Ruben's `ghostbuster_scattering_simulator.md` still describe the deleted
  simulator.

### Phase 5, documentation and report (2026-10-24, the Dell)

What it means: the work is documented for the people who will use and maintain
it. A policy states the contract of the wave operators and of
`reconstruct3D_wave` (conventions, scale link, architecture, tests, known
limits); a user-guide page explains the program; `simulate_scattering`'s page,
the user-guide index and the test-environment policy are current; the PCG
backend's policy names the routines where they now live; Ruben's note records
what changed in his simulator and why; and the run's report sits beside this
plan. Documents only: no source file changed in this phase.

Files changed: the new `doc/policies/3D/wave_model_policy.md`,
`doc/user_guide/reconstruct3D_wave.md` and
`doc/implementation_notes/planned/ghostbuster_operators_report.md`;
`doc/user_guide/user_guide.md` (index entries), `doc/user_guide/simulate_scattering.md`
(pointers to the policy and the program), `doc/policies/test_environment_policy.md`
(the sub-suites `nonlinear conjugate gradient`, `PCG lattice`, `wave field`,
`wave imaging`, `wave multislice`, `wave projector` and `wave reconstruction`
in their area suites), `doc/policies/3D/reconstruct3D_pcg_policy.md` (added to
the file table in this phase: the gather, scatter, wrap test and fold it names
moved into `pcg_lattice` in Phases 1 and 2, the item Phase 1 left for this
phase), `doc/implementation_notes/planned/ghostbuster_scattering_simulator.md`
(a dated section on the corrected conventions and the move onto the shared
operator, the contract made current, the first implementation kept as its
record), and this plan (status line; file table; the Phase 4 entry's particle
count: the Pt set holds 2 000 particles, 1 000 per half, which the logs' `N`
counts).

Exit criteria and how they were met:

- The documents read on their own: each defines its terms and acronyms on first
  use, names no item IDs, and points to the plan only for evidence. A second
  reader checked them against this plan and the code (findings and their
  fixes in `state/journal_phase_5.md`).
- `scripts/check_test_registry.py` and `scripts/check_descr.py` pass
  (`scratch/checks_p5.log`).

The second reader's findings, all fixed in the documents: the user guide gave
`objfun=cc` as the program's default, which is `euclid` (SIMPLE's default; the
example now passes `cc`, and the default is left open for the maintainer in
the report, since this phase changes no code); the pass buffer can fall back
to reading the stacks every pass; the three-tap lattice error is 0.7 to 1% of
the image contrast at its maxima (section 6 corrected too); the lens formula
counted the phase shift twice; the solver test's description; the refusals; the
prediction's `1/N^2` (section 3.2 corrected too); terms defined on first use;
the old simulator's first-Born constant-potential test has no direct successor
(the weak-phase identity and the model hierarchy cover the linear term).
