# Wave-model policy

Contract of the code that runs today: SIMPLE's wave operators, which predict a
particle image from a three-dimensional electrostatic potential with wave optics,
their hand-written derivatives, the program `reconstruct3D_wave` that
reconstructs half maps with them, and `simulate_scattering`, which draws data
with the same operators. The design came from Ghostbuster, the reconstruction
inside the specter package, and was rebuilt on SIMPLE's own machinery; the plan
with its measurements is
`doc/implementation_notes/planned/ghostbuster_operators.md` (sections 3 and 4,
and section 8 for the evidence of every number below), the run's summary
`doc/implementation_notes/planned/ghostbuster_operators_report.md`. The user
guides are `doc/user_guide/reconstruct3D_wave.md` and
`doc/user_guide/simulate_scattering.md`.

Terms used throughout:

- **Potential** `V`: the specimen's electrostatic potential on the native box
  `N^3`, voxel size `a` (`smpd`), in volts.
- **Padded lattice**: the `2N` grid on which SIMPLE's PCG (preconditioned
  conjugate gradients) reconstruction does every Fourier operation
  (`pcg_lattice`). The unknown stays on the native box.
- **Gather and scatter**: the Kaiser-Bessel (KB) interpolation that reads the
  volume's transform at a 3D point, and its exact transpose.
- **Ewald paraboloid**: the surface of the volume's transform that an image of a
  thick specimen samples, `kz = kappa lambda |k|^2 / 2` rather than the central
  plane `kz = 0`.
- **Exit wave** `psi`: the complex electron wave leaving the specimen. **Lens
  transfer** `T`: the objective lens's complex factor.
- **Weak-phase model**: SIMPLE's image model elsewhere, image = CTF (contrast
  transfer function) times the central section, which is the first-order limit
  of every wave model.

## 1. Scope

Fixed by the maintainer (rulings of 2026-10-23, in the plan):

- `reconstruct3D_wave` is a standalone shared-memory program. `reconstruct3D`,
  `refine3D`, the reconstruction service and their strategies do not use the
  wave operators; `rec_backend=wave` is reserved for that later integration.
  `nparts` is refused.
- Alignment stays on SIMPLE's flat references. Poses and CTF are inputs.
- No tomography.
- `simulate_scattering` runs on the same operators, with Kirkland's sign
  conventions and the defocus at the specimen midplane.
- Keys: `wave_model=projection|first_born|rytov|multislice`,
  `ewald=neg|pos|no`. No aliases of retired keys.

Not modelled: beam tilt and trefoil (the project has no fields for them),
anisotropic magnification, per-particle contrast scale, inelastic absorption,
detector response, ice and coherence envelopes. Specter itself is never a
dependency: nothing in the repository calls Python, and no test needs specter;
its outputs enter the tests as data.

## 2. Where the code lives

From the bottom of the import graph to the top. Imports point downward; no
domain module imports a strategy.

| File | Contents |
| --- | --- |
| `src/utils/math/simple_math_ctf.f90` | `electron_wavelength(kv)` and `electron_interaction_constant(kv)`, the one home of both (the CTF's wavelength) |
| `src/utils/math/simple_rnd.f90` | `poisson_rnd`, the detector's Poisson sampler |
| `src/main/ctf/simple_ctf.f90` | public `eval_phase`, the CTF's aberration phase including the phase shift (a required argument), used for the lens transfer |
| `src/main/volume/simple_pcg_lattice.f90` | the padded lattice shared with the PCG backend: `gather_plane` and `scatter_plane` with the curvature (`ewald_curv`, `plane_loc`, `rim_free`, `win_wraps`), `fold_lattice_accum` (with the Hermitian option for the two-cap planes of the curved models), and `gather_cube` and `scatter_cube`, a particle's rotated spectrum cube for multislice |
| `src/main/image/simple_wave_field.f90` | `wave_field`: a 2D complex field with its own c2c FFTW plans, the native-disk embedding into the twice-finer grid and the crop back, each with its transpose |
| `src/main/image/simple_wave_imaging.f90` | `wave_imaging`: one thread's 2D chain of one particle (exponential or linear exit wave, lens, intensity, crop) with its linearisation and adjoint; the lens and detection of a given exit wave for multislice; `detect_intensity`, the dose and detector noise of `simulate_scattering`; the model constants `WAVE_PROJECTION`, `WAVE_FIRST_BORN`, `WAVE_RYTOV`, `WAVE_MSLICE` |
| `src/main/volume/simple_wave_multislice.f90` | `wave_multislice`: one thread's multislice of a particle from its rotated cube's spectrum, with its linearisation and the adjoint sweep |
| `src/main/volume/simple_wave_projector.f90` | `wave_projector extends pcg_lattice`: the operator over a particle set (the volume's padded transform, the particles' geometry and optics, one `wave_imaging` and one `wave_multislice` per thread); loss and gradient, curvature along a direction, line losses, linearised planes, predictions, exit waves, the terms of the fitted gain |
| `src/main/opt/simple_nlcg_solver.f90` | `nlcg_objective` (abstract) and `nlcg_solve`: preconditioned Polak-Ribiere nonlinear conjugate gradients with a Gauss-Newton step refined by a parabola, generic for later CTF and pose refinement |
| `src/main/strategies/search/simple_matcher_ptcl_io.f90` | the observation-batch frame of reconstruction (`rec_state_half`, `rec_obs_oris`, `rec_obs_inv_sig2`, `read_rec_obs_batch`), one copy shared by the PCG worker and the wave objective |
| `src/main/strategies/parallelization/simple_rec3D_wave_strategy.f90` | `exec_rec3D_wave(params, build)`, the program body, and `rec3D_wave_objective extends nlcg_objective`, whose passes run on the caller's `params` and `build` |
| `src/main/commanders/simple/simple_commanders_rec.f90` | `commander_rec3D_wave`: defaults and refusals on `cline`, one `params` and builder, the strategy call |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `simulate_scattering` on `wave_projector` |
| `src/main/ui/simple/simple_ui_reconstruct3D.f90`, `simple_ui_sim.f90`, `simple_ui_params_common.f90` | the descriptors and the shared metadata of `wave_model` and `ewald` |

Each new stateful module holds one type with private components and symmetric,
idempotent `new`/`kill`. Testers sit beside their modules.

## 3. Image formation

For particle `i` with rotation `R_i`, shift `s` and CTF parameters:

```text
S_i(h,k) = gather(padf [h, k, kz] R_i) exp(+2 pi i (h s_x + k s_y) / N)    native disk, both signs of h
kz       = kappa wl (h^2 + k^2) / (2N),  wl = lambda / a,  kappa = -1 (ewald=neg), +1 (pos), 0 (no)
u        = IFFT_2N[ embed(S_i) ]
psi      = exp(c u)       projection (kappa = 0); Rytov (curved, see below)
psi      = 1 + c u        first Born
psi_T    = IFFT_2N[ T FFT_2N[psi] ]
I_hat    = crop_N( FFT_2N[ |psi_T|^2 ] )                                    native disk
p_i      = -g I_hat / N^2  (k /= 0)                                         SIMPLE's particle units (1/N^2: the unnormalised transform)
```

with `c = i sigma a (sqrt(1 - alpha^2) + i alpha)`, `sigma` the interaction
constant and `alpha = fraca`, so the amplitude contrast of the wave models is
the project's per-particle `fraca`. The rules:

- **Native disk.** `S_i` is zero beyond the native Nyquist radius; the residual
  is formed on the native disk in Fourier space, whitened by the particle's
  sigma2 when `objfun=euclid`, where the PCG backend forms its right-hand side.
- **Oversampling 2.** The exponential and the intensity are formed on the
  twice-finer grid, which is the padded lattice's size. Oversampling 1 errs by
  0.7 to 2.7% on the protein and Pt sets against 3e-4 to 8e-3 for 2; no key
  offers it.
- **Rytov at the exit face.** On the curved models the field is carried half a
  box downstream with Kirkland's propagator
  `P = exp(-i pi lambda (N a / 2) |k|^2)`, exponentiated there, and carried back
  to the midplane for the lens (`P^-1` folded into `T`). Without curvature
  Rytov is projection, bitwise.
- **Lens.** `T = exp(i eval_phase(phshift)) exp(-bfac |k|^2 / 4)`, with
  `T(0) = 1`: a phase plate's shift acts on the scattered waves only. SIMPLE's
  `eval_phase` is `-chi` in Kirkland's sign, so this is Kirkland's (and
  specter's) `exp(-i chi)`. Amplitude contrast lives in `c`, not in `T`.
- **Wavelength.** The CTF's, `12.26 / sqrt(V + 0.9784e-6 V^2)`, 3.4e-4 short of
  the exact relativistic value; shared with the CTF by design.
- **Multislice.** The particle's cube `C(h,k,l)` is gathered at
  `padf [h,k,l] R_i S_g` on the Nyquist sphere (`|l| < N/2`), times `N^3`, the
  shift phase and the symmetry mean; one transform along `l` gives the
  sections' spectra `s_j = (1/N) sum_l C exp(+2 pi i l z_j / N)`,
  `z_j = j - (N/2 + 1)`. Each section is embedded in the fine grid, transmitted
  by `exp(c u_j)` and propagated by `exp(-i pi lambda a |k|^2)` to the next; the
  exit wave is carried back to the box centre (`N/2 - 1` sections when the beam
  enters through the first section, `ewald=neg`; `N/2` through the last,
  `ewald=pos`) and passes `wave_imaging`'s lens and detector. No real-space
  rotated cube is formed: the cube gather and scatter are exact transposes, the
  same lattice and window as the reconstruction. The adjoint is the reverse
  sweep with `conj(P)` and `conj(t_j)`, checkpoints every `ceil(sqrt(N))`
  sections. Multislice needs a curvature sign (the beam direction); `ewald=no`
  is refused.

**The weak-phase identity** ties every convention to code refinement has used
for years: at zero potential the linearisation is `2 g sigma a N` times the PCG
backend's `forward_plane` times its `build_transfer` (CTF with shift, defocus,
astigmatism and amplitude contrast). The `wave projector` and `wave imaging`
tests hold it to `1e-5`.

## 4. Conventions

Measured against specter and RELION in the run's first phase (plan section
5.2); the wave operators and `simulate_scattering` follow the SIMPLE column.

| Quantity | SIMPLE | specter |
| --- | --- | --- |
| Euler angles | `e1, e2, e3` in degrees, `R = Rz(e3) Ry(e2) Rz(e1)` (SIMPLE's passive `rotmat`, the transpose of specter's), sections sampled at `[h,k,0] R`; RELION `rlnAngleRot, Tilt, Psi` | quaternion of `("ZYZ", [rot, tilt, psi])` |
| Shift | `x, y` in pixels, model-side phase `exp(+2 pi i (h x + k y) / N)`; RELION's origin is `x smpd, y smpd` in angstroms | translation in angstroms equal to the RELION origin |
| Defocus | `dfx, dfy` in micrometres, positive underfocus, at the specimen midplane (the box centre) | `dfu, dfv` in angstroms, the midplane |
| Astigmatism angle | `angast` in degrees from `+x` | `dfang`, the same value |
| Phase shift | `phshift` in radians, on the scattered waves only | radians (its star writer stores radians where RELION wants degrees) |
| Spherical aberration | `cs` in mm | angstroms |
| Lens | `exp(+i eval_phase)` | `exp(-i chi)`, the same |
| Propagator, transmission | Kirkland's: `exp(-i pi lambda d |k|^2)`, `exp(c V)` | the same |
| Image polarity | particles contrast-inverted, unit background noise; a positive potential is dark at underfocus in relative intensity | counts `D a^2 |psi|^2` |
| Volume axes | MRC, `x` fastest, centre voxel `N/2 + 1` | `(Z, Y, X)`, centre `N // 2` (0-based): the same voxel |
| Curvature | `ewald=neg`: a point `z0` voxels along `+z` images as the flat model at defocus `df + z0`; the beam enters through the first MRC section | `ews_curvature_sign="negative"` |

specter's positive sign refers its defocus one section off the box centre
(`N/2 + 1` sections from that sign's exit face, `N/2` for the other); SIMPLE
keeps the symmetric reference for both signs, so specter's positive-sign images
equal SIMPLE's at one voxel less defocus.

## 5. Units and the scale link

- **Images.** SIMPLE's particles are contrast-inverted and normalised to unit
  background noise; relative intensity is `1 - x / g` with the gain
  `g = sqrt(total_dose) smpd` at the particles' native pixel size (the
  normalisation precedes any Fourier crop, so `box_crop` keeps it).
  `total_dose` is required.
- **Maps.** By the weak-phase identity a SIMPLE map is `2 g sigma a N` times
  the potential in volts. The solver works in volts; the program writes its maps
  at SIMPLE's scale, so they compare directly with `reconstruct3D` maps.
- **Prior in volts.** The PCG backend's normal operator, density and ML
  (maximum-likelihood) prior, the per-shell regularisation it builds from the
  half maps' FSC (Fourier shell correlation),
  carry a factor `boxpd^3` (`boxpd = 2N`, the padded box); the prior in volts is `kappa_s^2 P / boxpd^3` with
  `kappa_s = 2 g sigma a N` (the `wave projector` test holds the link to
  `1e-4`).
- **The fitted gain** (`gain_fitted` in the diagnostics) is the one-scalar fit
  of the nonlinear part of the prediction to the data. It is a diagnostic, never
  used. It is meaningful for Rytov, projection and multislice; first Born's fit
  is not a gain. It falls short of the dose's gain when the particles' noise
  normalisation sees delocalised signal in the background (3.5 to 9% on the
  run's sets); it is the signal that the scale link of real data is off.

## 6. `reconstruct3D_wave`

Per state, in shared memory:

1. The halves come from the project's one maker, `ensure_field_halves`. The
   state's particles are read once into the in-memory pass buffer
   (`fill_pass_buffer`) when it fits the memory budget, and every later pass
   copies from it; otherwise (the log reports the buffer's status) every pass
   reads the stacks.
2. **Start**: each half's flat PCG solve (`maxits_pcg`), the ML prior from the
   FSC of the pair
   (`tau`, `hp`), the regularised flat map in closed form (`maxits_ml` coupled
   iterations on top). It is written as `recvol_stateNN_<half>_flat.mrc`.
3. **Wave solve** of each half from its start, by `nlcg_solve` with `maxits`
   (default 4) iterations: preconditioner the PCG backend's, step length from
   the exact slope and the Gauss-Newton curvature `||J d||^2` on a subset (every
   k-th particle of the half, 5% with at least 256 and at most 1000 particles),
   refined by a parabola through the subset's loss. Support: `pcg_mskfile` when
   given, else the `mskdiam` sphere.
4. **Outputs**: the regularised halves `recvol_stateNN_even.mrc`, `_odd.mrc`,
   their mean `recvol_stateNN.mrc`; the unfiltered halves `_even_unfil.mrc`,
   `_odd_unfil.mrc` (the Wiener filter undone isotropically, by the PCG
   operator's own shrink factor averaged over each shell); the FSC of the
   unfiltered pair (`fsc_stateNN.bin`, `RESOLUTION_STATENN`); per half
   `reconstruct3D_wave_stateNN_<half>_diagnostics.txt` (model, gains, solver
   history, passes).
5. **Project**: read once by the builder; the merged map and the FSC are
   registered in the `out` segment, written once at the end. The particle field
   is not written.

Refused: `nparts`; `automsk` other than `no`; `pcg_solvent=yes`; a missing or
non-positive `total_dose`; `objfun=euclid` without the project's noise spectra;
`wave_model=projection` with a curvature and `multislice` without one; a state
with particles in one half only; particles of a state at different voltages;
a project in which no state has particles.

Defaults: `wave_model=rytov`, `ewald=neg` (`no` for projection), `maxits=4`,
`winsz=1.5` (the PCG backend's three-tap KB window; `winsz=3` gives seven taps,
ten times smaller interpolation errors, at 1.4 times the run time, and changed
no band of the run's ribosome set).

## 7. Numerical invariants

- **Zero curvature is the PCG backend.** `gather_plane` and `scatter_plane` at
  `kappa = 0` are the PCG code, bitwise; every PCG result is byte-identical to
  the code before the move (the identity case of the plan's Phase 1, rerun in
  every phase).
- **Curvature bound.** `ewald_curv` refuses a paraboloid that would take a
  padded point more than 0.49 voxel past the Nyquist sphere, so no window
  reaches past the wrap table. The rim prefilter of the scatter includes `kz`.
- **Thread determinism.** Results are identical at any thread count: per-thread
  `wave_imaging`, `wave_multislice` and FFT plans in arrays sized by
  `omp_get_max_threads()` at the point of use; FFTW plans made and destroyed in
  the unnamed critical section (the one lock of every FFTW planner call in
  SIMPLE); scatters through the coloured h-strided pass with the serial rim
  pass; losses and dot products in double precision per particle, reduced in
  particle order; multislice gradient passes in chunks of one particle per
  thread, each chunk's cubes scattered in particle order. The projector refuses
  to run inside a caller's parallel region, and the lattice's gather opens no
  nested team inside the particle loop.
- **Single precision fields, double precision sums.**
- **Lattice accuracy.** The three-tap window interpolates a section to 0.7 to
  1% of the image contrast (1.3 to 1.8% of a gradient); seven taps reach 4e-5
  to 9e-4. Both lie far below a particle's noise; parity with specter is
  bounded by them.

## 8. Tests

Sub-suites (one tester module each; no CTest entry of their own, so
`SIMPLE_CTEST_BUDGET` is unchanged):

| Sub-suite (suite) | What it holds |
| --- | --- |
| `PCG lattice` (`unit_reconstruction`) | the curved scatter is the transpose of the gather with point-group replication (`1e-5`); the curved point lies at `padf [h,k,kz] R`; a point along the beam carries the paraboloid's phase and sign (0.05 rad); the cube gather's transpose (`1e-5`), its plane `l = 0` bitwise the flat gather; gather and scatter of planes and cubes bitwise identical at 1, 11, 16 and 24 threads, twice each |
| `wave field` (`unit_reconstruction`) | the forward transform is the direct DFT (`1e-5`), round trips (`1e-6`), the transposes of the transforms, the embedding and the crop (`1e-5`) |
| `wave imaging` (`unit_reconstruction`) | the transpose of the 2D linearisation (`1e-5`); the gradient against a fourth-order finite difference on half-radian fields (`1e-3`); the weak-phase limit `-2 sigma a CTF` (`1e-5`); first Born is Rytov's linearisation (`1e-6`); the image of a given exit wave is the chain's, bitwise, and its adjoint; the detector's counts are nonnegative integers |
| `wave multislice` (`unit_reconstruction`) | the sweep's adjoint (`1e-5`) and gradient (`1e-3`) at both beam directions; a weak sheet in the midplane images as the weak-phase model, 16 A downstream at `df + 16 A` and not at `df - 16 A`, with the beam reversed at `df - 16 A` (`5e-3`); zero and constant potentials give the vacuum and the analytic absorption |
| `wave projector` (`unit_reconstruction`) | the operator's transpose (`1e-5`) and gradient (`1e-3`) for every model, both signs, with point-group replication; the weak-phase identity against the PCG backend (`1e-5`); Rytov without curvature is projection bitwise; the vacuum; a point scatterer at height `z0` images at `df + z0` for `ewald=neg` (`2e-2`) and at least ten times worse at `df - z0`; line losses (`1e-5`); the scale link (`1e-4`); specter's golden images at box 16 and 24 (1.10e-2 of the contrast with three taps, 6.1e-4 with seven); multislice tends to first Born of its own sign (`2.2e-2`) and a smooth off-centre blob tells the signs apart (the opposite sign at least four times the error of the same); every output bitwise identical at 1, 11, 16 and 24 threads, twice each, for Rytov and multislice |
| `nonlinear conjugate gradient` (`unit_numerics`) | on a quadratic the engine reproduces linear PCG's iterates (`1e-4`), with and without the parabola refinement of the step, and run to convergence solves the system and, with a support projection, the restricted one (`1e-4`); on a nonconvex least squares it reaches the known minimum with a loss that never rises, halves overlong steps and refuses a step it cannot bring downhill |
| `wave reconstruction` (`lib_reconstruction`) | a phantom of 150 Gaussian atoms 100 A across, box 128 at 1 A, 300 particles imaged by multislice at 100 kV with detector noise; between the curvature limit (2.7 A) and Nyquist the truth FSC of Rytov with the data's sign beats the flat start and the opposite sign |

The tolerances are the plan's measured ones (section 8 of the plan, with their
evidence); they are not loosened to pass. A change to the operators reruns
`unit_reconstruction`, `unit_numerics`, `lib_reconstruction` and `pcg_recon`,
the PCG identity case, and the thread checks at 1, 11, 16 and 24 threads.

## 9. Known limits

- **The prior near the Ewald limit.** The ML prior is built from the flat
  start pair's FSC, which collapses exactly in the band where the flat model
  fails, so the prior is strongest there. On the run's Pt nanoparticle the wave
  maps repair only part of that band (truth FSC 0.07 to 0.20 between 1.12 and
  1.04 A against 0.40 to 0.63 for specter's converged Rytov); a ten times
  weaker prior (`tau=10`) recovers most of it. A prior from a weakly
  regularised first wave solve is the candidate; the program keeps the plan's
  prior until the maintainer decides.
- **Strong scatterers filling the box.** The native-disk design drops the
  content specter's native square holds beyond Nyquist; for a box-filling Pt
  particle under Rytov (2.6 rad) this moved images by 1% against specter.
- **Cost of multislice.** At box 160 about 0.9 s forward, 1.2 s adjoint and
  0.8 s cube gather per particle on one core: suited to nanoparticle boxes and
  small sets. A Rytov gradient pass costs 1.7 to 3 times a PCG accumulation.
- **Scale of real data.** Dose-weighted sums do not have Poisson background
  noise, and background-based normalisation sees delocalised signal; the
  fitted gain shows by how much the gain from the dose is off.

## 10. What waits for integration

`rec_backend=wave` in `reconstruct3D` and `refine3D` through the
reconstruction service; distributed execution (workers publishing gradient and
curvature sums, the master stepping); the nonuniform (NU) filter, trailing chains, state
weight sets and fractional updates; Ewald correction in the PCG and gridding
backends; CTF and pose refinement; tomography. The objective's particle pass is
the piece a worker would run, and `nlcg_solve` needs only gradient sums and
subset curvatures.
