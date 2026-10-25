# The atomicity prior on the PCG backend: open work

Status (updated 2026-10-25): the prior works inside `reconstruct3D`, and `autorefine3D_nano`
can bootstrap it from each round's `detect_atoms` model. It is not yet used by `refine_tomo` or
standalone `refine3D_nano`. Design, rationale, the rules every prior obeys, and the
record of what was built are in
[`../completed/single_pcg_integration_completed.md`](../completed/single_pcg_integration_completed.md)
(read its sections 3 to 6 for the maths and section 11 for the implementation choices).

Terms. PCG is the preconditioned conjugate-gradient reconstruction backend. The atomicity prior
models the map as `x = G c + e`, where each column of `G` is the density kernel of one atom at its
model position, `c` are free per-atom amplitudes, and only `e` is penalized, with strength
`lambda_A` times `x^T (I - P_G) x` (`P_G` projects onto the span of the atom kernels). AET is atomic
electron tomography. `refine_tomo` is the SINGLE tilt-series workflow. SINGLE is the nanoparticle
workflow of SIMPLE. A "held-out" check predicts the images or tilts that a map was not built from.

The learning so far. On FePt (68 tilts) with refined poses, the prior beats the prior-free map on
atom recall (0.44 against 0.42 at best prior-free), surface recall (0.18 against 0.14) and atom
roundness (elongation 1.80 against 2.07), but it only matches, not beats, the prior-free map on
held-out tilt prediction (0.0862 against 0.0852 at its best iteration count). The estimated
strength was within 1 percent of the best on the grid. Per-atom amplitudes fitted to the data
separated atoms from false peaks and Fe from Pt in a numpy prototype; that is not yet in the
program.

## 1. What is implemented (do not redo)

`reconstruct3D atm_prior=yes pdbfile= [lambda_atm= bfac_data= maxits_atm= rtol_atm=]` with
`rec_backend=pcg`: the module `src/main/volume/simple_pcg_atomic_prior.f90`, the operator and
preconditioner hooks in `simple_reconstructor_pcg.f90`, the plumbing in
`simple_rec3D_pcg_strategy.f90`, and the `atomic_prior` sub-suite of `unit_reconstruction`.
Choices fixed there: one Gaussian kernel for all species, `B = B_min + B_data`; the width fit uses
its own 1.5-voxel, 5-sigma model; the strength estimate uses the half-map noise; Jacobi-PCG Gram
solve in double precision with a ridge of 1e-7 of each diagonal entry; 30 iterations by default;
the prior'd halves are never a noise reference (use the prior-free `_base` halves).

## 2. Open work

1. **`refine_tomo` integration (first, because AET is the first target).** Version 1 changes only
   the final reconstruction: (1) the masked PCG base map from the refined poses, as now; (2)
   `detect_atoms` on it, in process, with the element list, giving the model including atoms
   recovered from the residual (the model errs toward inclusion: a spurious kernel costs little,
   a missing one pulls a real atom to zero); (3) `B_data`, the Gram matrix and `lambda_A` from the
   base map; (4) the regularized solve per half. Products: `tomo_final.mrc` (regularized),
   `tomo_final_base.mrc`, `tomo_final_res.mrc` (base map minus `G c`), and a sidecar
   `_atmprior.txt` (`B_data`, `N`, `S`, `C^2`, `lambda_A,rel`, `N_atoms`, CG statistics). Residual
   peaks above `3 sqrt(N)` are missing-atom candidates; atoms with `c_i` far below their species'
   median are removal candidates. Using the regularized map as the pose-refinement reference waits
   until the drift found on 2026-10-15 is understood. Files: `single_commanders_tomo.f90`, the
   strategy above.
2. **The held-out scan inside the program (`atm_cv`).** A new parameter. For a tilt series each
   half's regularized map predicts the other half's tilts; for particles it predicts the other
   half's images. Score the prediction error over the run's band on the grid
   `{0, 1/4, 1/2, 1, 2, 4} x lambda_A,rel`. This is how the estimate is checked, never by the FSC
   of prior'd halves.
3. **The TEM kernel.** The five-Gaussian electron scattering factor of each element, broadened by
   `B_i = B_min + B_data`, as `atoms%convolve` uses. Make `atoms%convolve`'s per-element tables a
   kernel builder shared by the simulator and the prior. The prior currently uses one Gaussian for
   all species, which suits annular dark-field STEM (AET) but not TEM data (SINGLE).
4. **SINGLE particles, partly implemented.** `autorefine3D_nano atm_prior=yes rec_backend=pcg`
   now runs stage 1 prior-free, builds a model with `detect_atoms`, and passes that PDB to the next
   stage's atomic-prior reconstruction. The published map is the regularized solve; resolution,
   the next atom detection, and species discovery use the prior-free `_base` pair. Every stage
   retains `_base`, `_res`, and `_atmprior.txt` beside the existing products. ICM is rejected.
   An element list such as `Pt,Ni` is passed only to `detect_atoms`, so the child refinement does
   not parse a species selector. With discovery enabled, the next prior uses `_species.pdb`, which
   includes residual-recovered weak sites, and all species products are retained per round. The
   same-stem species table now supplies the fitted per-atom B, occupancy and species-posterior
   expectation of the class/environment model intensity. Their Gaussian model defines the prior's
   one-dimensional unpenalized subspace; only
   its global scale stays free and is therefore set by the particle data. A PDB without this table
   keeps the previous position-only, free-amplitude prior. Species-dependent TEM kernels still need
   item 3. Standalone `refine3D_nano`, the full `single_workflow` acceptance trajectory and the
   held-out image check remain open.
5. **Per-atom amplitudes for species discovery.** Write the data-fitted amplitude
   `c = argmin || A G c - b ||^2` (solves `G^T H G c = G^T r` with the reconstructor's normal
   operator, by a small CG in atom space) and let `discover_species` use it instead of the aperture
   sum on the filtered map, which underestimates the Pt/Ni ratio by 23 percent and Pt/Al by 43
   percent on the synthetic particles. The amplitude measured on the map, `(G^T G)^-1 G^T x`,
   does not work (numpy prototype); do not use it. See
   [`species_discovery_unfinished.md`](species_discovery_unfinished.md).
6. **Acceptance runs of the design not recorded as done.** See section 3.
7. **Later, not now:** a smooth B field (not free per-atom B); a per-atom trust weight if surface
   or low-Z atoms fail the hard inclusion; non-negative amplitudes; the prior inside the
   pose-refinement references.
8. **Site visit demonstration.** The unsupervised pipeline run on the 2020/2021 SINGLE data sets,
   models compared with the curated ones. Needs items 3, 4, 5 and the image hold-out run of item 2.

## 3. Acceptance checks still to run

Written down before the runs, from the design:

- **Synthetic tilt series** (the atomic phantom of the PCG sparse-sampling test, sparse and dense
  tilt range). With the true model: rounder atoms and a cleaner background than prior-free at the
  same detection. With a damaged model (10 percent of atoms removed, half of them surface atoms;
  10 percent spurious added; positions jittered 0.2 Å; width off by 2; model map scaled by 10):
  the missing atoms are the top residual peaks, spurious atoms get `c` near 0 and no density,
  `B_data` does not depend on the scale, and detection on the regularized map beats prior-free.
  Score against the phantom's atoms, surface atoms (fewer than 12 neighbours) apart, over several
  noise realizations; not by FSC.
- **Amplitudes for species:** `c` against the aperture intensities on FePt and on the synthetic
  species cases, by coordination number.
- **Across data sets, at least four,** so that a calibration that holds only on FePt is caught:
  multislice simulations with Peter Ercius (NCEM) with exact atoms; a single-species AET set (the
  tungsten tip of Xu et al., 2015); an AET Pt/Ni set if obtainable; the SINGLE Pt particle data.
- **Default strength rule.** The estimated `lambda_A,rel` becomes the default only if, on every
  test set above, its held-out prediction error is within 1 percent of the best on the grid.
  Otherwise `lambda_atm` stays an explicit input. The test is on the error, not the strength,
  because the optimum is flat. Met so far on FePt only.
- **SINGLE particles:** held-out image prediction and atom statistics (item 4).

## 4. Rulings and rules that bind

- No ICM step; version 1 has one prior and one strength.
- Tomography first, then SINGLE.
- Rules every prior obeys: a symmetric positive semi-definite term; a linear operator frozen during
  a solve; nothing accumulated into the raw accumulators; attached only in deapodized kernel mode;
  strength zero reproduces the prior-free result bit for bit; half independence kept or its loss
  declared; resolution reported from the prior-free pair only; judged by truth, held-out data or
  against a model, never by the FSC of prior'd halves; one variable per experiment, thresholds
  written down before a run.
- Strength anchor: `lambda_A = lambda_A,rel x realspace_diag`, with `lambda_A,rel` the Wiener ratio
  of the completed note's section 6. Unmeasured voxels keep a zero in the preconditioner (giving
  them `1/lambda_A` stalled CG).
- Amplitudes of the model never enter the prior (scale-free), so simulated and reconstructed maps
  need no common scale.
- Related companion changes: the rule in the species design that atoms with amplitude below 0.2
  of the median are removal candidates would flag a whole light species; reference it to the
  noise instead.
