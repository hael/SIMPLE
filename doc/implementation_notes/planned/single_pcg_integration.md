# The atomicity prior on the PCG backend, for SINGLE and atomic electron tomography

Implementation note. First version 2026-09-15; revised 2026-10-15 at master
`e2673a5d7`. Design for review; nothing here is implemented. Companions:
`pcg_backend_overview.md` (the solver), `reconstruct3D_pcg_policy.md`,
`pcg_decision_log.md`, `pcg_priors_history.md` (the priors tried before),
`aet_codes_and_simple.md` (the FePt work) and
`doc/refactoring_notes/planned/refine_tomo.md`.

## 0. What the revision changes

- **No ICM step.** The first version began by reproducing the ICM filter of
  `refine3D_nano` inside the solver. ICM is little used any more and never on
  denoised data, so that step and its Laplacian term are gone. Version 1 has
  one prior and one strength.
- **Tomography first.** The prior is tried first on atomic electron
  tomography (AET) through `refine_tomo`, then on SINGLE particle data. A tilt
  series is fixed and small, a run takes minutes, and the FePt data come with
  a published atomic model; simulated tilt series with exact atoms are
  planned with Peter Ercius (NCEM).
- **The solver has changed.** Unmeasured Fourier voxels are now excluded from
  the preconditioner (`PCG_RHO_SAMPLED_FRAC`), a real-space solvent ridge
  exists with its mean folded into the preconditioner, the data operator's
  real-space diagonal (`realspace_diag`) is available, and a density mask can
  be the solve support (`pcg_mskfile`). Section 5 builds on these.
- **The lessons of the rejected priors are rules here** (section 2), above
  all how a strength must be set so that it carries across data sets.

## 1. Why now: what the FePt tilt series showed

The FePt nanoparticle of Yang et al. (2017), 68 projections from −65° to
+64°, reconstructed in SIMPLE and scored against their published model
(23,196 atoms, 28% Fe) with the same atom detector:

| Map | Level-1 recall | Pt above background | Fe above background | Atom elongation |
| --- | --- | --- | --- | --- |
| SIMPLE gridding | 0.25 | 63% | 4% | 2.00 |
| SIMPLE PCG, 5 iterations | 0.35 | 66% | 2% | 2.03 |
| `refine_tomo` (refined poses, density mask, PCG) | 0.42 | 74% | 3% | 2.07 |
| Their GENFIRE map | 0.41 | 90% | 7% | 1.89 |

"Above background" is the share of atoms whose smoothed density at the true
position exceeds the 99th percentile of the density between atoms;
elongation is the ratio of the longest to the shortest axis of the mean Pt
atom. What these numbers say:

- **Detection is matched, the background is not.** Least squares with a
  density support plus pose refinement detects as many atoms as GENFIRE.
  Their map is still cleaner between atoms.
- **Positivity is not the reason.** A GENFIRE-style iteration with positivity
  on our data made the atoms as round as theirs (1.90) and the map free of
  negative values, but lifted no more atoms above the background, and recall
  fell from 0.30 to 0.24.
- **Every map smears atoms along the beam.** The missing 50° wedge doubles
  the atom length in that direction, in their map too.
- **Fe is barely visible anywhere.** Almost no Fe atom is detected, so
  species discovery splits the Pt in two instead of finding Fe.

The atomicity prior addresses all three directly: it pulls the space between
atoms to zero, it extrapolates into the missing wedge with atom-shaped
density instead of nothing, and it returns a least-squares amplitude for
every atom, which is what species discovery needs (section 3.3).

## 2. Lessons from the earlier priors

From `pcg_priors_history.md`, `pcg_decision_log.md`, the NU-evidence
sharpening note and LocScale (`~/src/locscale`, v2.4.1):

- **Strength anchored to the low-band data scale does not carry over.** The
  biharmonic smoothness prior, the binary solvent prior and the NU-evidence
  precision all set their strength relative to `data_scale`, the mean of the
  data diagonal over the six lowest shells. Their effect then depended on how
  steeply the data diagonal falls across the band, a property of each data
  set: the same coefficient was inert on PfCRT and strong on msp1. The NU
  prior worked (it beat everything on the 1WCM truth phantom) and was removed
  for that reason, after a controller that steered it by the FSC of
  non-independent halves wrecked ten runs.
- **A spectral target is not the operative ingredient.** The Wilson prior (a
  shell-diagonal precision from a reference spectrum) helped at both ends of
  the band but was beaten everywhere by local evidence. The atomic kernels
  carry the Wilson envelope themselves: a kernel with width B has power
  `|f(s)|² exp(−B s²/2)`, and the positions add the structure factor. No
  separate spectrum is needed.
- **The solvent ridge shows the mechanics.** It is a real-space diagonal term
  `λ_s (1 − w(r))` with a protein weight built per half from that half's own
  prior-free map. The closed-form regularized solve cannot represent it, so it
  needs a real solve, and without its mean in the preconditioner two
  iterations stalled. The FSC is never computed from a prior'd pair.
- **LocScale.** Model-based amplitude references sharpen well. Its hybrid
  mode fills unmodelled density with pseudo-atoms, the same split as atoms
  plus a free residual. It smooths refined B factors over 25 Å, which argues
  for a smooth B field over free per-atom B. Randomly placed pseudo-atoms
  impose a liquid-like structure factor that real proteins do not have, so a
  kernel prior needs real atoms or a real model; nanoparticles and AET
  specimens have them. Its network mode warns about hallucinated density.
  Its Wilson completion is no longer called in v2.4.1.
- **Rules every prior obeys.** A symmetric positive semi-definite term; a
  linear operator frozen during a solve; nothing accumulated into the raw
  accumulators; attached only in deapodized kernel mode; strength zero
  reproduces the prior-free result bit for bit; half independence kept or its
  loss declared; resolution reported from the prior-free pair only; judged by
  truth, by held-out data or against a model, never by the FSC of prior'd
  halves; one variable per experiment, thresholds written down before a run.

## 3. The prior

### 3.1 Definition

The map is modelled as `x = G c + e`. Column `g_i` of `G` is the density
kernel of atom `i` at its model position, `c` are free per-atom amplitudes,
and `e` is what the atoms do not explain. Only `e` is penalized:

    Q(x) = λ_A · xᵀ (I − P_G) x,        P_G = G (GᵀG)⁻¹ Gᵀ

`P_G` is the orthogonal projector onto the span of the atom kernels, so
`I − P_G` is symmetric and idempotent and the term is already in
normal-equation form. The regularized solve is

    (P H P + λ_A (I − P_G)) u = P b

on the solve support `P`, solved by the existing CG.

### 3.2 Properties

- **Scale-free.** `P_G` does not change when a column of `G` or the map is
  rescaled. The model contributes positions, species and kernel shape, never
  amplitudes, so the simulated and reconstructed maps never need a common
  scale.
- **Vacuum.** Where no kernel reaches, `P_G x = 0` and the term is
  `λ_A ‖x‖²`: the space between atoms is pulled to zero by the same operator.
  This is the useful half of positivity, kept linear.
- **The missing wedge.** The kernels have full Fourier support. Atom
  amplitudes set by the measured data therefore come with the matching
  density in the unmeasured wedge, while `λ_A` suppresses whatever else would
  fill it. Atoms should come out rounder without a nonlinear constraint.
- **One parameter spans the range.** `λ_A = 0` is the present solve. As
  `λ_A → ∞` the solution is exactly `G ĉ`, atoms with least-squares
  amplitudes. No separate atomic-basis solver is needed.

### 3.3 Per-atom amplitudes

Species discovery now classifies an aperture sum on a filtered map. On the
synthetic particles that sum underestimates the Pt/Ni ratio by 23% and the
Pt/Al ratio by 43%.

The replacement is the amplitude fitted to the data, `ĉ = argmin ‖A G c − b‖²`.
It solves `Gᵀ H G c = Gᵀ r` with the reconstructor's own normal operator `H`
and right-hand side `r`, by a small CG in atom space. This is the strong-prior
limit of the regularized solve.

The amplitude measured on a map, `(GᵀG)⁻¹ Gᵀ x`, does not work. In the numpy
prototype on FePt (2026-10-15, `~/agent_runs/runs/atm_prior_mac/numpy_proto/RESULTS.md`)
it separated atoms from false peaks no better than the peak height. A
two-class mixture on it did not split Fe from Pt at all.

The amplitude fitted to the projections did both:

- **Atoms against false peaks:** recall at 90% precision rose from 0.75 to
  0.85.
- **Fe against Pt:** the unsupervised two-class split called 76 to 80% of
  the Fe atoms Fe at about 70% precision.

Refining each atom's position against the data (Gauss-Newton steps, with the
amplitudes refitted after each step) raised the Fe recall further. It is
tested on its own (section 9, step 4).

### 3.4 Half independence

The atom positions come from the merged base map, so they are a shared
frozen field. Both halves are then pulled toward the same atoms, and the FSC
of the regularized pair is inflated. As for the solvent ridge: resolution is
reported from the prior-free base pair only.

## 4. Kernels and the one fitted width

- **TEM data (SINGLE).** The five-Gaussian electron scattering factor of each
  element, as `atoms%convolve` uses, broadened by `B_i = B_min + B_data`.
- **Annular dark-field STEM (AET).** Intensity scales roughly as `Z^1.7`,
  and less predictably for low-angle detectors. The prior ignores amplitudes,
  so only the kernel shape matters: a Gaussian per species, whose width comes
  from the fit below. The FePt refinement of Yang et al. used one B per
  species, which is a precedent for a per-species width.
- **The width fit.** On the base map `x_b` with the model simulated at
  `B_min`, the per-shell regression slope `c_s = Re⟨m̂_s, x̂_s⟩ / ‖m̂_s‖²` is
  unbiased under map noise. Fit `log c_s = log K − B_data s²/4` over shells
  where the model correlates with the map at 0.5 or more; discard `K`; clamp
  `B_data ≥ 0`. One `B_data` in version 1. A smooth B field, the LocScale
  lesson, comes later; free per-atom B does not.
- **Footprints and Gram matrix.** Kernels are evaluated within `cutoff`
  (4 to 8 pixels) and clipped to the solve support, so `P_G` lives on the
  solve domain. `GᵀG` is sparse (only overlapping footprints couple) and
  positive definite as long as atoms are separated. It is solved by an inner
  CG or a banded Cholesky built once per model. When kernels overlap strongly
  (wide B at low resolution), the Gram solve gets a small ridge `ε I`.

## 5. Operator and preconditioner

- **Operator.** `apply_normal` adds `λ_A (p − P_G p)` next to `add_ridge`,
  after the support mask. `P_G p` is a gather, a Gram solve and a scatter, all
  in real space; no FFT is added. Cost per CG iteration: about
  `N_atoms × footprint` operations, small next to the FFTs.
- **Preconditioner.** The denominator of the measured voxels gains `λ_A`
  (in the kernel's units, divided by `padsc²` as for the ML prior). The
  unmeasured voxels keep the present zero. The numpy prototype on FePt
  (2026-10-15, `~/agent_runs/runs/atm_prior_mac/numpy_proto/RESULTS.md`)
  tried both. Giving the unmeasured voxels the gain `1/λ_A` stalled CG at
  small strengths and gave worse maps at every strength. With the zero, CG
  converged at every strength, and the wedge still filled through the support
  coupling: Pt atoms went from twice as long along the beam as across to 1.4
  to 1.8 times.
- **Iteration budget.** Not the production two iterations. The prior acts
  through the coupled solve, so the budget is set by convergence on FePt and
  on the synthetic tilt series (an `rtol`, with `maxits` as a cap; 20 to 30
  expected). The kernel operator's behaviour when over-iterated is watched
  with the existing indefinite-stop guard and checked against the
  matrix-free operator on the small fixtures.
- **What is untouched.** The raw accumulators, the ML prior and its closed
  form (not used on these paths), the solvent ridge.

## 6. Strength

The earlier priors failed on their anchor, not their form. This one is set
so that it does not depend on how steeply the data fall across the band.

- **Anchor.** `λ_A = λ_A,rel × realspace_diag`. `realspace_diag` is the mean
  of the data diagonal over all native shells, the data precision per voxel
  in real space. A real-space prior competes with that, not with the low-band
  mean.
- **Estimate.** `λ_A,rel` is a Wiener ratio from moments of the base map:

      N   = noise variance per voxel
      S   = max(mean(x_b²) over the support − N, ε N)       signal power
      E   = mean((P_G x_b)²) over the support               power the atoms explain
      E'  = max(E − N · N_atoms / |V|, 0)                   minus the noise in span(G)
      C²  = min(E'/S, 0.98)                                 explained fraction of the signal
      λ_A,rel = (N/S) / (1 − C²)

  `N` comes from the half-map difference where the halves are independent
  (an AET tilt series split by interleaved tilts: `var((x_even − x_odd)/2)`
  inside the support), and from the empty region of the box otherwise. The
  ratio is the sigma-A weight of crystallographic map combination, built from
  the map's own moments; no constant is tuned.
- **Check by held-out data, never by FSC.** For a tilt series, each half's
  regularized map predicts the other half's projections, which are other
  tilts. For particles, each half's map predicts the other half's images
  (the cross-validation of the first version). The score is the prediction
  error over the run's band, on a grid of `{0, ¼, ½, 1, 2, 4} × λ_A,rel`.
- **Acceptance.** The estimate becomes the default only if, on every test
  set of section 9, steps 3 and 5, its held-out prediction error is within
  1% of the best on the grid. Otherwise `lambda_atm` stays an explicit input.
  The test is on the error, not on the strength, because the optimum is
  flat. In the numpy prototype on FePt (every fourth tilt held out), the
  estimate came out at λ_rel 0.5 to 0.7 and the held-out optimum at 0.1.
  The error differed by 1% across that range. Every regularized map
  predicted the held-out tilts better than the prior-free map at its best
  iteration count (0.117 against 0.121).

## 7. Where it sits

### 7.1 Tomography (`refine_tomo`)

Version 1 changes only the final reconstruction:

1. The masked PCG base map from the refined poses, as now.
2. `detect_atoms` on the base map with the element list, in process, giving
   the model, with the atoms species discovery recovers from the residual.
   The model errs on inclusion: a spurious kernel costs little (its
   amplitude goes to zero), a missing one pulls a real atom toward zero.
   Surface atoms are the ones at risk, being the weakest and partly
   occupied on real maps.
3. `B_data`, the Gram matrix and `λ_A` from the base map.
4. The regularized solve per half.

Products: `tomo_final.mrc` (regularized), `tomo_final_base.mrc`,
`tomo_final_res.mrc` (`x_b − G ĉ`, the residual on the base map) and a
sidecar `_atmprior.txt` (`B_data`, `N`, `S`, `C²`, `λ_A,rel`, `N_atoms`, CG
statistics). Peaks of the residual above `3√N` are missing-atom candidates;
atoms with `ĉ_i` far below their species' median are removal candidates.
Using the regularized map as the reference of the pose refinement comes
later, after the drift found on 2026-10-15 is understood.

### 7.2 SINGLE particles (`refine3D_nano`, `autorefine3D_nano`)

The loop is an EM iteration. `detect_atoms` (the M-step) builds the model
after each refinement stage. The next stage's reconstructions use the
regularized solve (the E-step) as the alignment reference, with the base
solve kept for resolution and for the residual. Stage 1, without a model,
runs prior-free. ICM is not combined with the prior.

## 8. Code map

- **New `src/main/volume/simple_pcg_atomic_prior.f90`, type
  `pcg_atomic_prior`.**
  - State: positions in crop voxels, species, widths, radial kernel tables,
    packed footprints, the Gram matrix and its factor, `λ_A`, the support.
  - Procedures: `new`/`kill`, `build`, `gather`, `scatter`, `solve_gram`,
    `project`, `apply` (`λ_A (x − P_G x)`), `amplitudes` (fitted to the
    data through the reconstructor's normal operator, section 3.3),
    `residual`, `fit_width`, `estimate_lambda`, `stats`.
- **`simple_reconstructor_pcg.f90`.**
  - `set_atomic_prior`.
  - `apply_normal` adds `prior%apply` in both operator modes.
  - `finalize_density_accum` adds `λ_A` to the denominator of the measured
    voxels; the unmeasured ones keep zero.
  - `get_atomic_prior_stats`.
- **`simple_rec3D_pcg_strategy.f90`.** The base solve, the prior build and
  the regularized solve per half; the held-out scan when asked; the products
  and the sidecar.
- **Parameters.** `atm_prior` (yes/no), `lambda_atm` (relative; 0 uses the
  estimate), `bfac_data` (override of the fit), `atm_cv` (the held-out scan).
  `atm_prior=yes` requires `rec_backend=pcg` and a model (`pdbfile`, or the
  in-process detection in `refine_tomo`).
- **`single_commanders_tomo.f90` and `single_commanders_nano3D.f90`.** The
  model handoff, the parameters and the products.
- **`atoms%convolve`.** Its per-element tables become a kernel builder shared
  by the simulator and the prior.

## 9. Order of work and the gates

1. **The module.** Unit tests on a synthetic box:
   - `P_G` is a projector (`P_G² = P_G`, symmetric, `P_G G = G`);
   - `apply` is symmetric to 1e-6;
   - `P_G` does not change when columns are rescaled;
   - the Gram solve is accurate;
   - in the reconstructor, strength zero reproduces the prior-free solve bit
     for bit, and the kernel and matrix-free operators agree with the prior
     on.
2. **Synthetic tilt series.** The atomic phantom of the PCG sparse-sampling
   test, with a sparse and a dense tilt range.
   - With the true model: rounder atoms and a cleaner background than
     prior-free, at the same detection.
   - With a damaged model (10% of atoms removed, half of them surface
     atoms, 10% spurious added, positions jittered 0.2 Å, width off by a
     factor two, model map scaled by 10): the missing atoms are the top
     residual peaks; the spurious ones get `ĉ ≈ 0` and no density; `B_data`
     does not depend on the scale; and detection on the regularized map
     beats the prior-free map.
   - Scored against the phantom's atoms, surface atoms (fewer than 12
     neighbours) apart from the rest, over several noise realizations; not
     by FSC.
3. **FePt through `refine_tomo`.** The `λ_A` grid around the estimate, with
   the thresholds written down before the run:
   - level-1 recall and precision, Fe recall, recall of surface atoms
     (fewer than 12 neighbours in the published model);
   - background separation and atom elongation;
   - held-out tilt prediction error;
   - distance to the published model.
4. **Amplitudes for species.** `ĉ` against the aperture intensities on FePt
   and on the synthetic species cases of `species_discovery`, by
   coordination number.
5. **Across data sets.** At least these four, so that a calibration that
   holds only on FePt is caught:
   - the multislice simulations with Peter Ercius (NCEM), with exact atoms;
   - a single-species AET set (the tungsten tip of Xu et al., 2015);
   - an AET Pt/Ni set if obtainable;
   - the SINGLE Pt particle data.

   The estimate's acceptance rule (section 6) is judged here.
6. **SINGLE particles** in `refine3D_nano` and `autorefine3D_nano`, judged
   by held-out image prediction and the atom statistics of `detect_atoms`.
7. **Later.** A smooth B field; a per-atom trust weight (if surface or
   low-Z atoms fail the hard inclusion); non-negative amplitudes; the prior
   in the pose-refinement references.

## 10. Related methods

- **Crystallography.** Atomicity has long been used to extend phases:
  Sayre's equation (1952) and the squaring methods, and density modification
  (solvent flattening, histogram matching). In cryo-EM this was taken up by
  Terwilliger et al. (2020). ARP/wARP (Perrakis et al., 1999) alternates
  placing free atoms and recomputing the map, the same loop as section 7.2.
  Charge flipping (Oszlányi & Sütő, 2004) is a dual-space iteration that
  needs no prior model.
- **Radio interferometry.** CLEAN (Högbom, 1974) deconvolves an incompletely
  sampled Fourier plane by subtracting point sources and keeping a residual
  map, the same idea as the residual of section 7.1. Self-calibration
  (Cornwell & Wilkinson, 1981) alternates the source model with the
  instrument calibration. It cannot recover the absolute position or flux
  scale, much as our pose refinement drifted in the tilt scale and frame
  rotation.
- **Materials STEM.** Statistical atom counting from annular dark-field
  intensities (StatSTEM, De Backer et al., 2016) fits a Gaussian per atomic
  column and counts atoms with a mixture model, using integrated intensities
  as section 3.3 does. Discrete tomography (DART, Batenburg & Sijbers, 2011)
  and atom counting gave 3D atomic positions of nanoparticles (Van Aert et
  al., 2011). De Backer et al. (2022, npj Comput. Mater. 8:216) build a 3D
  Pt particle from a single projection: per-column atom counts with their
  mixture-model probabilities, a neighbour-count prior and an embedded-atom
  energy, searched by a genetic algorithm. They judge the result by the
  share of surface atoms placed correctly against simulated truth over 30
  noise realizations, the score adopted in section 9. Their time series
  shows surface atoms moving under the beam while the atom count stays
  fixed, consistent with the partly occupied outer sites seen on the SINGLE
  Pt maps.
- **Sparse recovery.** Off-the-grid point-source recovery (Candès &
  Fernandez-Granda, 2014) and the atomic norm (Chandrasekaran et al., 2012)
  formalize when a sum of point sources can be recovered from limited
  Fourier data.
- **Cryo-EM maps.** Gaussian-mixture representations of maps (Chen & Ludtke,
  2021) and denoiser priors inside refinement (RELION Blush, Kimanius et al.,
  2024). The latter show the gain of strong priors and the risk of density
  the data do not support.

## 11. Status (2026-10-16)

Implemented for `reconstruct3D` (`atm_prior=yes pdbfile= [lambda_atm= bfac_data= maxits_atm= rtol_atm=]`):
`simple_pcg_atomic_prior` with its unit tests (`unit_reconstruction`, sub-suite `atomic_prior`), the hook in
`simple_reconstructor_pcg` (sections 5 and 8) and the plumbing in `simple_rec3D_pcg_strategy`. Not yet done:
`refine_tomo` and `refine3D_nano` integration, the held-out scan inside the program (`atm_cv`), the TEM
five-Gaussian kernel and the `atoms%convolve` kernel builder.

Choices made in the implementation:

- **Kernel.** One Gaussian for all species, `B = B_min + B_data`, `B_min` a Gaussian of one voxel.
- **Width fit.** It simulates its own model with a Gaussian of 1.5 voxels and a 5-sigma footprint. The narrowest
  kernel aliases and a truncated footprint adds a high-frequency plateau; both bent the per-shell slope by 6%
  in the unit test.
- **Noise of the estimate.** The solve is per half, so the estimate uses the half-map noise `2 var((xe-xo)/2)`.
  The merged ratio of section 6 is logged beside it and is half of it.
- **Gram solve.** Jacobi-preconditioned CG in double precision with a ridge of `1e-7` of each diagonal entry.
- **Budget.** 30 iterations by default. On FePt the CG residual is 1.5e-3 at 30 and 1.2e-3 at 60, and the atom
  scores move by 1 to 2% between the two.
- **Half maps.** The prior'd halves are no noise reference: species discovery on them overcalls. Use the
  prior-free `_base` halves.

FePt with the refined poses of `2_refine_tomo` (30 iterations, model = `detect_atoms` with discovery on the base map):

| Map | Level-1 recall (precision) | Surface recall | Elongation | Data error, 68 tilts | Held-out error |
| --- | --- | --- | --- | --- | --- |
| Prior-free, 5 iterations | 0.416 (0.983) | 0.144 | 2.07 | 0.0690 | 0.0862 |
| Prior-free, best (10 iterations) | | | | | 0.0852 |
| Prior-free, 30 iterations | 0.372 (0.972) | 0.162 | 2.08 | 0.0691 | 0.0890 |
| Prior, 1/4 x estimate | 0.457 (0.985) | 0.196 | 1.93 | 0.0658 | 0.0858 |
| Prior, estimate (0.415) | 0.442 (0.986) | 0.182 | 1.80 | 0.0674 | 0.0862 |

The held-out errors come from separate maps of 51 tilts that predict the other 17. Unlike the numpy prototype,
the prior does not predict held-out tilts better than the prior-free map at its best iteration count; it matches
it. On this set the estimate is within 1% of the best strength on the grid (section 6).
