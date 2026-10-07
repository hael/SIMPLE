# NU filtering: switch the candidate cost to the plain Euclidean loss

Date: 2026-10-06. Status: completed 2026-10-06. The code change is built and
its unit tests pass; the documents and the skill describe the new cost;
section 4 is decided as option 1. The refine3D and PCG comparison runs of
Phase 2 were dropped by decision (see Progress); the neutral fixture lives on
as a unit test that reports the envelope recall on every build.

## 1. What is there now

The nonuniform (NU) filter chooses, voxel by voxel, among a bank of low-passed
versions of the half maps by the cross-half prediction error

```text
r1_c(v) = [E(v) - O_c(v)] / sigma(|v|),    r2_c(v) = [E_c(v) - O(v)] / sigma(|v|),
C_c(v)  = H(r1_c(v)) + H(r2_c(v)),
```

where `H` is the Huber loss with transition 1.345 (introduced 2026-05-21) and
`sigma(|v|)` is a radial noise profile: the Gaussian-scaled MAD of the
real-space difference `E - O` in 4-pixel shells of distance from the box
centre (`image%nu_objective_noise_profile`, introduced 2026-08-26). The costs
are smoothed with a tent kernel, an ordered-label Potts prior regularizes the
label field, and the same cost and profile feed the NU evidence envelope and
the PCG solvent-prior strength estimate.

## 2. Why it changes

Under the Huber loss the noise scale decides the winner: a residual is
quadratic or linear depending on how it compares with `1.345 * sigma`, so the
scale has to be right voxel by voxel, and that is what the radial profile
tries to supply. It cannot supply it without damage:

- **The profile hides the disagreement NU is meant to find.** In real data
  `E - O` is larger over density than over solvent: heterogeneity, alignment
  errors and other genuine disagreement between the halves add to the noise.
  A profile estimated from `E - O` therefore inflates `sigma` where the
  disagreement is, shrinks the residuals there, and pushes exactly those
  regions toward finer labels. Radius is also the wrong partition: a compact,
  centred particle puts density in the inner shells and solvent in the outer
  ones, and flexible peripheral parts inflate their own shells.
- **The Huber tail discounts the wrong residuals.** Reconstruction noise is a
  sum over very many Fourier coefficients and is Gaussian; there are no heavy
  noise tails to guard against. The large residuals come from real signal that
  a too-coarse candidate removed, and the linear tail under-penalizes exactly
  that over-smoothing. Isolated voxels are already limited by the tent
  smoothing of the costs.
- **The plain squared error needs no local scale.** With
  `C_c(v) = (r1_c(v)^2 + r2_c(v)^2) / 2` and residuals divided by any
  candidate-independent `sigma(v)`, the winner at a voxel does not depend on
  `sigma(v)`. Noise that is higher in some region (the gridding correction
  raises it toward the box edge) makes the finer candidates costlier there,
  which is the correct response for that map. This is the criterion of
  cryoSPARC's non-uniform refinement (Punjani, Zhang and Fleet, Nature Methods
  2020): cross-validated squared error between one half map, filtered, and
  the other, over a local window.

## 3. Design

- **Cost.** `C_c(v) = (r1_c(v)^2 + r2_c(v)^2) / 2` with the residuals divided
  by one global level `sigma_0`, so that the cost equals the present Huber
  cost wherever that cost is in its quadratic regime and the Potts prior keeps
  its calibration. `sigma_0` affects only the balance between the data term
  and the Potts prior, never the per-voxel ranking.
- **Noise level.** `sigma_0` is the Gaussian-scaled MAD of `E - O` over the
  observed voxels of the NU support (voxels where both halves are exactly zero
  stay excluded). A median over the whole support is not inflated by a
  minority of disagreeing density, unlike a median over a radial shell.
- **Removed.** The Huber loss and its constant, `nu_objective_noise_profile`,
  the cached profile and its provenance keys (`whitening_shells`,
  `whitening_checksum`, `whitening_rmax`).
- **Unchanged.** The static bank and its cut at FSC 0.143 / 1.5, the
  regularized member, the tent smoothing radii, the ordered Potts prior, the
  synthesis, and the handoff to matching.
- **Solvent-prior strength (PCG).** `estimate_solvent_prior_lambda` and its
  second caller choose the prior strength by the mean cross-half prediction
  cost; they take the same squared-error cost with `sigma_0`.
- **Evidence null.** The lower-quartile centring of the null margin stays; it
  was introduced for an offset of the Huber cost and is re-measured, not
  assumed, under the squared error.

## 4. Open decision: the evidence envelope

The evidence envelope thresholds cost margins (`C_zero - min C_signal`)
across voxels, so unlike the per-voxel competition it does depend on the
local noise level. The radial profile was introduced for it: with one global
scale, the gridding correction's noise rise at the periphery cost recall of
true density (0.48 against 0.61 with the profile, neutral fixture, commit
`53f173787`). Options:

1. Global `sigma_0` only, and measure the recall (simplest; may fall back
   toward 0.48).
2. Divide the margins by the known geometric factor `g(v)^2` of the
   reconstruction: the gridding correction for gridding maps, the solve-support
   weight for PCG maps. Exact and not density-dependent, so it hides no
   disagreement, but it needs the factor passed from volume assembly to NU.

Proposed: option 1 first; option 2 only if Phase 2 shows the recall loss.

## 5. Phases

- **Phase 0, baseline.** On a build of the base state: the neutral fixture
  (envelope recall of true density, component count, NU label histogram,
  in-band ratio), the NU unit tests, one `refine3D` with nonuniform filtering
  on beta-galactosidase, and one PCG reconstruction with the solvent prior.
  Record every number with its command line.
- **Phase 1, implementation and tests.** The design of section 3. Tests: on a
  phantom with stationary noise the per-voxel ranking is unchanged when
  `sigma_0` is multiplied by any constant; on a phantom with a region of
  known coarser local resolution and extra half-to-half disagreement, the
  region receives its coarser label; the existing NU, evidence and PCG
  solvent-prior tests pass with their floors re-derived where a floor was set
  in Huber units, each change recorded with its reason.
- **Phase 2, comparison and documents.** Repeat Phase 0. Report label
  histograms, envelope recall and component counts side by side; `refine3D`
  final FSC 0.143 within the larger of the run-to-run spread and 3 %; the
  solvent-prior strength before and after. Then take the section 4 decision
  on the evidence numbers, and update
  `doc/policies/NU/nonuniform_filtering_policy.md`,
  `doc/algorithms/nonuniform_filtering.md`,
  `doc/algorithms/nu_evidence_envelope_mask.md` and the skill
  `simple-nonuniform-regularization`.

### File table

| File | Change | Phases |
|---|---|---|
| `doc/refactoring_notes/completed/nu_euclidean_loss.md` | plan and progress | each |
| `src/main/image/simple_image_calc.f90`, `simple_image.f90` | squared-error `nu_objective`, global `sigma_0`, profile routine removed | 1 |
| `src/main/nu_filt/simple_nu_filter_bank.f90`, `simple_nu_filter_evidence.f90`, `simple_nu_filter_vars.f90`, `simple_nu_filter_state.f90`, `simple_nu_filter_envmask.f90` | callers, cached level, provenance, null centring | 1 2 |
| `src/main/volume/simple_pcg_solvent_sidecar.f90` | solvent-prior callers | 1 |
| unit testers of the files above | tests | 1 |
| `doc/policies/NU/nonuniform_filtering_policy.md`, `doc/algorithms/nonuniform_filtering.md`, `doc/algorithms/nu_evidence_envelope_mask.md` | documents | 2 |
| `.github/skills/simple-nonuniform-regularization/` | skill | 2 |

## Progress

### Phase 1, 2026-10-06: implemented, not yet compiled

The code change of section 3 is in the working tree as an uncommitted diff;
compilation and the test run are with Hans. Phase 0 (the baseline on the base
state) has not been run; it needs a build of commit `9b77876a3` before this
diff is applied.

What changed:

- `image%nu_objective` computes `C_c(v) = (r1^2 + r2^2) / 2` with both
  residuals divided by one scalar `noise_scale`; the loss is formed in double
  precision and capped at `1/epsilon` as before; a non-finite residual still
  costs the cap. The Huber loss, its constant and the radial profile routine
  (`nu_objective_noise_profile`) are gone.
- `image%nu_objective_noise_scale` existed already (no callers since the
  profile replaced it) and now is `sigma_0`: the Gaussian-scaled MAD of
  `E - O` over the voxels of the support where the halves are not both exactly
  zero, with the profile's finite check, its RMS fallback and the unit level
  for an all-zero pair.
- `simple_nu_filter`: `setup_nu_dmats` computes `sigma_0` once and caches it
  (`nu_noise_scale_cached` replaces the cached profile and its radius); the
  evidence null candidate uses the same level. The three provenance keys
  `whitening_shells`, `whitening_checksum`, `whitening_rmax` are replaced by one
  key `noise_scale`, and `NU_EVIDENCE_ALGORITHM` is bumped to `nu_evidence_v2`
  so that a frozen evidence state built under the old cost can never share an
  identity with one built under the new one. Log lines and comments naming
  the Huber cost or the whitening profile were rewritten.
- `simple_pcg_solvent_sidecar`: `estimate_solvent_prior_lambda` and
  `solvent_prior_cross_half_objective` take `sigma_0` from
  `nu_objective_noise_scale` on the prior-free pair and call the squared-error
  cost, as section 3 specifies.
- The lower-quartile centring of the evidence null is unchanged, to be
  re-measured in Phase 2.
- Section 4 follows the proposal: option 1, global `sigma_0` only; the
  envelope recall is measured in Phase 2 before option 2 is considered.

Tests. The plan assumed existing NU, evidence and PCG solvent-prior unit tests
whose floors would be re-derived; there are none (the standalone NU test
programs were deleted on 2026-09-23 because they asserted nothing), so nothing
was re-derived. A new tester, `src/main/nu_filt/simple_nu_filter_tester.f90`,
is registered as the sub-suite `nonuniform filtering` of `unit_reconstruction`
(suite table, test UI help and the policy's suite table updated;
`scripts/check_test_registry.py` and `scripts/check_descr.py` pass). It
checks:

- the unary against its closed form, voxel by voxel, and that one noise level
  of residual in each term costs exactly one;
- `sigma_0` against the known standard deviation of `E - O` (within 3 %),
  unchanged when half of the support is set to exact zero in both halves, and
  equal to one for an all-zero pair;
- on a phantom under stationary noise (a white Gaussian field smoothed with a
  Gaussian of 0.7 px, unit standard deviation, plus independent noise of
  standard deviation 0.5 per half; four candidates low-passed at 12, 8, 6 and
  4.5 Å): multiplying `sigma_0` by 3.7 or 0.2 scales every unary by the
  inverse square and changes the winner at no voxel whose top-two margin
  exceeds single-precision rounding;
- on the two-resolution phantom through `setup_nu_dmats` and
  `optimize_nu_cutoff_finds` (box 64 at 2 Å, support sphere of 100 Å: the
  signal is the 0.7 px field left of the box centre and the same white field
  smoothed with 2.5 px right of it, both of unit standard deviation; noise
  0.5 per half everywhere and an extra independent 0.7 per half on the coarse
  side): leaving 8 px either side of the split out, the median selected
  low-pass is at or finer than 6.5 Å on the fine side and at or coarser than
  9.5 Å on the coarse side, each holding for at least 90 % of the side's
  voxels.

The thresholds of the last test were chosen on a numpy emulation of the
competition (Butterworth order 8 bank at the ladder cutoffs, the squared-error
unary, tent smoothing at the candidate radii with mask normalisation, the
coarse-to-fine like-for-like selection and the coarsen-only ordered-label ICM
with beta from the mean top-two gap). On three seeds the emulation selected a
median of 4.9 Å on the fine side and 11.6 to 14.2 Å on the coarse side with
every voxel of each side on the asserted side of the threshold, before and
after the Potts pass; the asserted margins are therefore wide. The emulation
is throwaway scratch, not part of the repository.

Build and fast gate, 2026-10-06 (Hans): the diff compiles and
`unit_reconstruction` passes with the new sub-suite in 3.4 s for the whole
area suite, so the full competition at box 64 fits the budget. The gate's one
failure, `unit_stream` / initial analysis / `test_rebuild_init_mics`, was a
pre-existing race in `stream_watcher%watch` unrelated to this change (the
settle check read the clock before the forked directory listing, and an
access to a project by another process across a second boundary made its
age −1, which a settle time of −1 rejects); fixed in
`src/main/stream/simple_stream_watcher.f90` as a separate hunk. The whole
gate then passed.

### Phase 2 and close-out, 2026-10-06

Hans's decision: finish with the documents, the skill and a resurrected
neutral fixture; the refine3D (beta-galactosidase) and PCG solvent-prior
comparison runs of Phase 2, and with them the Phase 0 baseline on the base
state, are dropped. The comparison numbers the plan asked for (label
histograms, final FSC within 3 % or the run-to-run spread, the solvent-prior
strength before and after) are therefore not recorded; a refinement with
nonuniform filtering on real data is the first place a regression would show,
and the labels it selects are visible in the `_nu_locres` map and the
`>>> NU candidate label assignments` table under `NU_DEV_OUTPUT`.

The neutral fixture (the deleted `simple_test_nu_envmask` program: a sphere of
radius 15 px of band-limited common signal in a noise support of radius 22 px,
box 48 at 2 Å, uniform per-half noise of amplitude 1, two fixed linear
congruential streams) is now `test_evidence_envelope` of the tester. It builds
the compact evidence state (valid, null fraction within the readiness bounds,
nine candidates, four bands, provenance naming `nu_evidence_v2` and carrying
`noise_scale`), checks that the coarse-band support is higher inside the
molecule and the explicit null wins preferentially in the solvent, that the
evidence margin at 4 Å separates the two by more than a factor of two, and
segments the envelope at 3 null MADs with the binary prior; it prints and
asserts the envelope's recall of true density (floor 0.90) and its solvent
false-positive rate (ceiling 0.15). On the same pair zeroed outside a hard
radius of 18 px it checks that `sigma_0` moves by less than 5 % (the zero/zero
exclusion), that the observed fraction equals the hard support's share of the
sphere, and that every unobserved voxel is frozen at the null with maximal
uncertainty and no band support. The numpy emulation of the envelope path
(the margin at 4 Å, the median/MAD null over the observed support, the
degree-normalised ICM) gave recall 1.00 at a false-positive rate of 0.07 on
this fixture, and `sigma_0` of 0.437 against 0.430 on the hard-supported pair
(the true standard deviation of `E - O` is 0.408); the floors have the
corresponding margin. The build's own numbers are on the test's output line
`envelope recall of true density`.

Section 4 is decided as option 1: one global `sigma_0`, nothing passed from
volume assembly to NU. The fixture's recall under the squared error is at the
ceiling, so on it the global level costs nothing; the plan's concern was the
periphery of real maps under the gridding correction, which this fixture does
not model (its noise is stationary). Option 2, dividing the margins by the
reconstruction's known geometric noise factor, stays open as the remedy if
real-data evidence envelopes prove tight at the periphery; its design is in
section 4 and needs the factor handed from assembly to `setup_nu_dmats`.

The fixture's first run on Hans's debug build (2026-10-06, bounds-checked)
stopped in `regularize_evidence_labels` on
`if( l_constrain .and. nu_solvent_lmask(i,j,k) ) cycle`: Fortran evaluates
both operands, and without a solvent envelope the array is unallocated. A
release build had been reading past an unallocated descriptor on every
evidence build without an envelope. The test is nested now (`if( l_constrain )`
then the array test), and a sweep of `src/` found no other
`allocated(x) .and. x(...)`-style hazard. The rest of the fixture is
unmeasured until the next build.

Documents updated: `doc/policies/NU/nonuniform_filtering_policy.md` (section 6
and the evidence paragraphs), `doc/algorithms/nonuniform_filtering.md` (model
and rationale), `doc/algorithms/nu_evidence_envelope_mask.md` (objective and
margin sections, the fixed-choice table, the tests pointer), the
`simple-nonuniform-regularization` skill and its reference map, and a status
note on the planned `nu_evidence_mask_cross_validated_threshold.md`, whose
formula was written for the Huber cost. The dated history documents under
`doc/implementation_notes/completed/` keep their Huber-era wording.
