# Probabilistic in-plane search: mathematical and performance review

**Latest review:** the 2026-11-03 end-to-end mathematical audit is in section 15.
Sections 12-14 retain the previous integration reviews; sections 1-11 retain the original review.

**Reviewed:** 2026-10-30, `master` at `a2c1b1f86`.
**Scope:** commits `519c67d2a` and `a2c1b1f86`, the proposed integration into solve2D/solve3D, and the relevant existing probability-table and reconstruction contracts.
**Status:** review findings and recommendations; no implementation changes made or recommendations adopted by this review.

The implementation/decision record is [polar_shift_search.md](polar_shift_search.md). This document is the separately requested review record, not a second implementation plan. The related history is in [pose_marginalization.md](pose_marginalization.md), especially sections 11, 20.8–20.13 and 21.

## Assessment

The Fourier–Bessel approach is mathematically sound and worth pursuing. It provides an efficient way to evaluate the existing discrete polar objective over translations and in-plane rotations without Cartesian regridding. The proposed posterior construction needs corrections before it represents a joint posterior, and the recorded timings do not yet establish that integration will be fast.

The distinction matters: correct evaluation of a likelihood at every pose does not by itself establish correct probability masses, marginalization, sampling, or reconstruction updates.

The new machinery is not connected to the production search at the reviewed revision. Findings about Phase 2 therefore concern the proposed integration, while findings about the existing table sampler describe a pre-existing approximation rather than a regression introduced by these two commits.

| Priority | Finding | Required resolution |
| --- | --- | --- |
| P1 | Equal weights on the polar translation grid use the wrong integration measure | Include translation cell areas and handle the origin explicitly |
| P1 | Drawing a pose per class, then selecting a class from that pose's loss, is not a joint posterior draw | Select classes/projections using integrated in-plane evidence |
| P1 | L-BFGS applied after a posterior draw changes its distribution | Separate optimization from posterior sampling |
| P1 | Screening candidates at the current shift can exclude the hypothesis needed to recover a large offset | Establish translation-aware candidate coverage and an exploration mechanism |
| P1 | The Phase 1 calculator submodule is absent from the reviewed tree | Restore the implementation and reproduce its tests and timings |
| P2 | Total-cost estimates omit the number of mapped references per particle | Benchmark actual candidate counts and complete iterations |
| P2 | Prior, frame, precision and calibration contracts need stronger validation | Test the probability model as well as the correlation identity |

P1 means a prerequisite for claiming a correct, reproducible joint-posterior implementation. P2 means a material qualification or validation requirement.

## 1. What is mathematically sound

The Jacobi–Anger expansion separates the shift phase into Bessel factors and angular harmonics. On the calculator's conventions, the shifted reference is rotated into the particle frame, so the cross spectrum must be indexed in the reference frame. The mechanics class follows that convention.

The documented all-rotation identity,

```text
a_n(irot, rho) = (1 / nrots) sum_m exp(-2 pi i m (irot-1) / nrots)
                sum_k w_k J_n(z_k(rho)) Rhat_k[m-n] conj(Yhat_k[m]),
```

is consistent with the discrete Fourier conventions, with `Y` the CTF-modulated particle. Friedel extension gives the stated symmetry of the coefficients and the factor of one half in the zero-order contribution to the half-ring objective.

An independent Python/SciPy calculation checked this identity against direct phase-ramp summation at 1,440 probes, using an angularly varying CTF. The maximum absolute cross-term error was `1.57e-6`, or `4.89e-10` relative to the Cauchy–Schwarz scale. This supports the algebra; it is not a test of the absent production submodule.

The following choices are also appropriate:

- Preserve the existing objective's shell weights, including `k/sigma2` for Euclidean matching.
- Reuse memoized angular spectra rather than computing shifted-reference forward transforms at every translation.
- Build geometry and Bessel tables once, and keep mutable evaluation work per thread.
- Use one numerical implementation for 2D and 3D.

“Exact” should be qualified as identity with the discrete polar objective, subject to harmonic truncation and floating-point error. Exact evaluation of a truncated trigonometric polynomial does not establish an exact continuous maximum or an accurate posterior integral on a coarse grid.

## 2. Correct the polar probability measure

**P1 — plan ruling 5, the flattened-grid draw.**

The translation prior is a density in Cartesian shift coordinates. Under `s = rho (cos alpha, sin alpha)`, the area element is

```text
ds_x ds_y = rho d_rho d_alpha.
```

Uniform radial and angular spacing therefore does not give equal-area cells. Applying `draw_posterior_row` directly to their losses overweights small radii. For radial cell boundaries `rho_minus`, `rho_plus`, the translation cell area is

```text
Delta_A = 0.5 * (rho_plus^2 - rho_minus^2) * Delta_alpha.
```

At `rho=0`, all directions describe the same shift. Represent its central disk once, or divide its area across the duplicate entries. Merely multiplying by `rho` and assigning zero weight to the origin is not a satisfactory finite-cell treatment either.

For a flat likelihood, a Gaussian prior with per-axis `sigma=10 A`, truncated to a `30 A` disk, gives these mean radii:

| Integration measure | Mean radius |
| --- | ---: |
| Cartesian area, equivalently correctly weighted polar cells | 12.30 A |
| Equal weights on uniform polar coordinates, in the fine-grid limit | 7.91 A |

This is a probability-model error, not a small discretization effect. Changing the number of angular samples per ring makes explicit area weights even more important.

## 3. Select the class using evidence, not a sampled loss

**P1 — Phase 2a's instruction to keep the class draw unchanged.**

The current machinery draws an in-plane rotation for each reference and stores that rotation's loss. The table then exponentiates those stored losses to select the class or projection. The new plan extends the first draw to translations while retaining the second operation.

That construction is not equivalent to drawing from a joint posterior. Let `f_c(j)` be the unnormalized mass of in-plane pose `j` for class `c`. The required class weight is `Z_c = sum_j f_c(j)`, not `f_c(j_c)` for a conditional draw `j_c`.

A finite counterexample with equal class priors and equal pose-cell sizes is:

```text
Class A: in-plane likelihoods [1.000, 0.010]
Class B: in-plane likelihoods [0.505, 0.505]
```

Both classes have evidence `1.01`, so the correct class probability is `P(A)=0.5`. Exact enumeration of the draw-then-rescore procedure gives `P(A)=0.658065`. The sharper class is favored despite equal evidence.

Relevant production sites are `gen_posterior_val` in [simple_polarft_corr.f90](../../../src/main/pftc/simple_polarft_corr.f90), `ref_assign` in [simple_eul_prob_tab.f90](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90), and `assign_posterior` in [simple_eul_prob_tab2D.f90](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90). Section 11 of the earlier pose-marginalization note already identifies the row as an approximate direction sampler.

### A coherent and inexpensive replacement

For particle `i`, class/projection candidate `c`, and in-plane cell `j`, construct

```text
ell_cj = -beta * T_i * L_cj
         - |u_cj|^2 / (2 sigma^2)
         + log(pi_c)
         + log(Delta_V_cj)

log Z_c = logsumexp_j(ell_cj).
```

Here `u_cj` is the absolute particle shift relative to the pick, `pi_c` is the intended class/direction prior, and `Delta_V_cj` includes the relevant angular and translation integration measures. Count each prior or measure factor once. Define `T_i` as the untempered `posterior_temperature`; apply `beta` once.

Draw `c` proportional to `Z_c`, and draw `j` conditionally within `c`. The calculation can retain one conditional pose sample and one `log Z_c` per candidate while evaluating its map. Once the candidate is selected using evidence, its saved conditional sample is already the correct one. A second expensive map evaluation is unnecessary.

Keep the evidence and the sampled pose's loss in separate fields. They serve different purposes, and evidence must not be exponentiated again as if it were a raw normalized loss. If the reference set is truncated, this is a posterior conditional on that retained set; its omitted mass remains a separate approximation to assess.

## 4. Optimization after sampling is a different algorithm

**P1 — the Phase 2a map-to-L-BFGS handoff.**

A posterior draw followed by deterministic L-BFGS is no longer distributed according to that posterior. In a single Gaussian basin, every optimized sample reaches the same mode, eliminating all within-basin variance. In a multimodal objective, the result is a distribution over optimized basins induced by their attraction regions.

This affects both reconstruction and the interpretation of sigma2 residuals. The residual at a posterior sample has a posterior-expectation interpretation; the residual after optimizing that sample does not inherit it.

For probabilistic operation, use optimization to locate basins or construct proposals, and sample on an adaptively refined grid or through a transition that preserves the target posterior. Keep deterministic polishing as an explicitly different estimator or comparison arm.

The statement that a posterior narrower than the grid becomes “deterministic in effect” is not a sufficient accuracy argument. Such a posterior can lie between grid points, giving an incorrect evidence, class ranking and apparent certainty. Refine where probability mass concentrates and check convergence under further refinement.

## 5. Candidate screening must permit large-shift recovery

**P1 — Decisions fixed and ruling 12.**

Running the expensive map on a bounded candidate set is a reasonable performance decision. However, selecting that set at the current shift can exclude the correct class precisely because the particle is miscentred. A correct physical displacement is common to the particle, but its inference is coupled to the reference, orientation and reference-centering convention. Searching an incorrect candidate does not guarantee recovery of that displacement for the next iteration.

Retain a bounded budget, but establish a mechanism that covers this failure mode: broad translation-aware screening at low resolution, exploration candidates, or another proposal scheme with a defined treatment of selection bias. An MCMC alternative would need proposal corrections and mixing measurements; a few transitions are not automatically independent posterior draws.

On a small particle subset, compare the selected candidates with exhaustive joint evaluation. Measure retained posterior mass and recovery of the true class/projection across offset magnitude and SNR. Do not use fixed-shift ranking alone as evidence that the discarded candidates are unimportant.

## 6. Prior, frame and likelihood-accuracy requirements

The planned untempered prior is a coherent choice, but the existing `set_shift_prior` in [simple_pftc_shsrch_grad.f90](../../../src/main/pftc/simple_pftc_shsrch_grad.f90) encodes a penalty in normalized loss units using `1 / (2 sigma^2 T)`. Exponentiating that combined loss at `beta T` tempers the prior too. Reusing that path unchanged would violate the new contract. Assemble the likelihood and prior separately, including at very small `beta`.

The map's optimizer-frame shift, stored particle-frame shift, and absolute shift relative to the pick must be distinguished. Existing assignment paths add increments to the stored shift. The absolute prior must evaluate the resulting absolute displacement, with the correct rotation/frame conversion, rather than the increment alone. Define the disk's support in that same frame.

The calculator builds `argtransf` using native `box` dimensions in [simple_polarft_core.f90](../../../src/main/pftc/simple_polarft_core.f90). Convert physical sigma, radius and grid spacing into those pixels. Validate cropped and uncropped paths and reference recentring with nonzero incoming shifts; fixtures with no crop and zero incoming shifts cannot establish this contract.

The proposed `sigma=10 A` is a modeling choice, not a value inferred from the stated mean radial displacement of `20 A`. An isotropic two-dimensional Gaussian has `E[rho] = sigma sqrt(pi/2)`, so that mean would imply `sigma approximately 16 A` under that model. A `3 sigma` disk excludes approximately 1.11% of a two-dimensional Gaussian's mass. View-dependent biases may require a different model; recovered shift histograms alone are not unbiased measurements of picking errors.

Measure numerical error in the units the sampler uses: `beta T_i Delta_L`, alongside errors in evidence and posterior mass. Agreement of `exp(-L)` to a fixed absolute tolerance is not by itself a uniform guarantee on posterior odds. The Bessel cutoff controls individual omitted coefficients, not directly the total log-likelihood error.

## 7. Source completeness and reproducibility

**P1 — the Phase 1 submodule is missing at the reviewed revision.**

The plan reports an implementation in `src/main/pftc/simple_polarft_shiftmap.f90`, but that file is absent from both HEAD and the checkout. [simple_polarft_calc.f90](../../../src/main/pftc/simple_polarft_calc.f90) declares its bindings; [simple_polar_shiftmap_tester.f90](../../../src/main/pftc/simple_polar_shiftmap_tester.f90) calls them; and `polarft_calc%kill` calls `kill_shiftmap` unconditionally.

The reviewed source therefore lacks the definitions needed to link these calls. This is a source inspection finding, not an observed linker diagnostic: a full build was not attempted. The recorded Debug/Release passes and Phase 1 timings cannot be reproduced from the reviewed commit until the submodule is restored.

The existing kernel fixture also uses `ctf=no`, identical reference halves and no crop. After restoring the implementation, extend identity coverage to angularly varying CTF, distinct half references and noise spectra, cropped geometry, and the prior/frame cases above. The independent calculation in section 1 does not replace those production tests.

## 8. Performance assessment

The plan's recorded timings are useful kernel measurements, but not evidence that a complete integrated iteration will be fast. They are historical measurements quoted from the plan, not reproduced in this review.

| Fixture | Existing 5 x 5 coarse search | Map at the 5-pixel window | Map at the 30 A disk |
| --- | ---: | ---: | ---: |
| Box 64, lp 8 A | 0.178 ms/pair | 0.096 ms/pair | 1.093 ms/pair |
| Box 64, lp 2.7 A | 0.479 ms/pair | 0.938 ms/pair | Not measured |

The fine-resolution window case is slower than the existing coarse search. These grids also differ in coverage and sampling, so equal-accuracy and equal-coverage comparisons are still required.

The note's roughly eight-second iteration estimate assumes one mapped reference per particle. In 3D, `calc_num2sample` in [simple_eul_prob_tab_utils.f90](../../../src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90) uses

```text
num_smpl = min(num_all, max(1, int(athres * num_all / 180))).
```

For 1,000 directions and a 10-degree bound, this is 55 references, not a count bounded numerically by `prob_athres`. Using the recorded `1.093 ms/pair`, 200,000 particles, and ideal scaling across 24 threads:

| Mapped references per particle | Map time alone |
| ---: | ---: |
| 1 | 9.1 s |
| 10 | 91.1 s |
| 55 | 501 s, approximately 8.3 min |

These are optimistic extrapolations of total map work, not measured incremental overhead. Integration replaces some existing work, fractional sampling changes the particle count, and thread scaling may be imperfect. A complete comparison must account for all three.

### Optimization priorities

1. **Replace direct direction synthesis with batched FFTs.** `map_rotcoeffs` in [simple_polar_shiftmap.f90](../../../src/main/pftc/simple_polar_shiftmap.f90) loops over harmonics and directions for every rotation and radius. Its cost is approximately `Nrot * Nalpha * sum_r Nharm(r)`. FFT synthesis is the natural first target. Choose transform lengths that preserve the retained harmonics; do not silently truncate them to fit a shorter transform.
2. **Use appropriate sampling per ring and refine significant cells.** Inner rings need fewer angular samples. Preserve cell areas and validate evidence under refinement; a grid adequate for initializing a local maximum is not necessarily adequate for posterior integration.
3. **Reuse spectra and weighted tables.** Amortize particle/reference spectrum preparation and Bessel-weighted shell tables where their dependencies permit it. Avoid pair-loop allocation, and benchmark cache/memory effects before duplicating large memoized arrays.
4. **Reduce maps while they are produced.** Accumulate stable evidence and retain conditional samples by block, rather than keeping a full map for every candidate. The evidence correction in section 3 should require reductions over values already evaluated, not another search.
5. **Measure complete stage costs.** Include actual candidate counts, representative centring masks, matching bands, CTF, crop, per-thread memory, and both one-thread and full-team execution. Report search time and iteration time, including saved local-optimization work.

A speed claim should compare equal posterior accuracy and recovery quality, not only time per grid point.

## 9. Qualifications for the surrounding marginalization work

The draw-centred `topk` cone is useful local averaging, but is not generally exact Rao–Blackwellization of the original draw. Its retained support depends on the randomly committed pose, and the current insertion path uses untempered weights while the committed draw uses `prob_beta`.

Even with equal temperatures and no top-K truncation, overlapping draw-centred neighborhoods can change the expected insertion. For probabilities `[0.2, 0.5, 0.3]` and retained neighborhoods `{1,2}`, `{1,2,3}`, `{2,3}`, averaging the posterior within the drawn neighborhood gives expected insertion weights `[0.157143, 0.580357, 0.262500]`, not the original probabilities. Its empirical benefit can stand without claiming exact preservation of the estimator's expectation.

One correct posterior draw does reproduce posterior-weighted linear accumulators in expectation for a fixed model. That does not imply the restored map, a ratio of accumulators with subsequent processing, has exactly the same expectation. The earlier note already records this distinction.

RELION's formulation marginalizes class and pose and accumulates posterior-weighted reconstruction and noise statistics, with adaptive coarse-to-fine sampling to reduce cost. A corrected SIMPLE sampler can provide a Monte Carlo approximation to those statistics without duplicating RELION's complete implementation. See [Scheres, 2012, sections 2.1 and 2.3](https://pmc.ncbi.nlm.nih.gov/articles/PMC3690530/).

Keep known-pose simulations as machinery/calibration tests. Finished-run poses are useful real-data proxies, but do not establish truth or eliminate all dependence on the earlier inference. For the discrete highest-density quantile proposed in section 21 of the pose-marginalization note, uniformity requires suitable randomization within probability atoms/ties, or comparison with the appropriate discrete null distribution.

## 10. Acceptance evidence still needed

| Check | What it establishes |
| --- | --- |
| Flat-likelihood prior recovery at several radial/angular resolutions | Correct polar cell masses, origin handling and untempered prior |
| Equal-evidence class counterexample | Correct hierarchical joint sampling rather than draw-then-rescore bias |
| Single Gaussian and multimodal in-plane posteriors | Correct within-basin variance and basin probabilities |
| Grid-refinement comparison of evidence, moments and credible mass | Posterior integration accuracy, including narrow off-grid peaks |
| Cropped/native and recentred-frame equivalence | Physical shift units, absolute-prior anchor and increment conversion |
| Direct-sum kernel identities with angular CTF and distinct halves | Objective and half-selection correctness in production |
| Exhaustive small-subset comparison | Candidate retention and large-offset recovery |
| Known-pose simulations over offset, SNR and stage | Recovery accuracy and credible-set coverage |
| Repeated full-team execution | Thread safety beyond the one-thread fast gate |
| Complete-stage timing at actual candidate budgets | Practical speed and memory cost |

## 11. Verification performed in this review

- Confirmed the reviewed revision and a clean starting worktree.
- Inspected the mechanics class, calculator interfaces and lifecycle, testers, current table draws, shift prior, candidate-budget formula, and weighted-insertion consumers.
- Ran an independent Fourier–Bessel/direct-sum calculation with fixed seed `20261030`, `box=64`, `pftsz=48`, shells `2:11`, radii `[0, 1.5, 4, 9, 15]`, all 96 rotations, and three random directions per radius/rotation. Angular CTF values were nonuniform; shell weights were `k / (0.5 + 0.1 k)`. Results are in section 1.
- Exactly enumerated the two-class sampling counterexample and the draw-centred cone example; numerically integrated the flat-likelihood radial prior examples.
- Passed `python3 scripts/check_test_registry.py . --verbose`.
- Passed `python3 scripts/check_flex_dag.py . --verbose`.
- Passed `python3 scripts/check_descr.py .`.

No Fortran build, area suite, fast gate, Linux/BOX test, end-to-end recovery run or performance benchmark was executed in this review. This checkout had no `build/`, and the required Phase 1 submodule was absent. Historical results quoted from the implementation note remain historical evidence rather than newly verified results.

## 12. Integration re-review, 2026-10-31

**Scope:** `04f6f0cb7` through `4114f9f4a`, plus the 26 files with
uncommitted changes at review start. This covers the stored-frame prior
correction, retirement of the n_eff controller, tempering changes,
weighted insertion and its cone, the optimized polar kernel, both 2D
fills, the dense 3D fill, the stage controllers and centring/mask handoff.
Unrelated trajectory, nanoparticle and hardware-planner developments are
outside this shift-search review. Late infrastructure additions present
at the final build snapshot (22:36 local) were also inspected. Further 2D
implementation edits arrived during the final gate; they remain outside
this closed review. No production source was edited by this review.

### Improvements verified in the current source

The missing calculator submodule is present. The kernel's identity checks
now exercise angular CTF, distinct reference halves and noise spectra, and
report errors in likelihood units. The sampler includes polar cell areas,
normalizes the absolute prior separately from the tempered likelihood,
and retains the origin's finite cell area. Refined map candidates carry
integrated evidence for the reference draw, and the map's committed pose
is not passed through the old in-plane optimizer. The original shift prior
now follows the stored shift at each angle, including the joint angle
derivative. The old n_eff solve and process-wide beta feedback are gone.

The remaining issues are chiefly at the interfaces between these pieces.

### R1 — P1: the extra-pose conversion double-rotates window-search entries

`simple_strategy3D_utils::extra_pose` unconditionally calls
`shift_ref2ptcl`. However, `refine_best_refs` and `refine_best_neighbors`
call `minimize(sh_rot=.true.)` and write its returned shift directly to
the table. That return already applies exactly `shift_ref2ptcl`:
`cxy(2:) = matmul(cxy(2:), rotmat)`. The passed prior-frame tests also
establish that this is the stored particle-frame shift.

Consequently `topk>1`, gridding and `shift_sigma=0` now rotate an extra's
increment twice. A native table increment (2,-1) at 90 degrees becomes
(1,2), rather than remaining (2,-1). The map writes optimizer-frame cells,
so its entries do need the conversion. One table now has two possible
frame contracts, and the new helper assumes the map's for both. Normalize
entries at the producer or distinguish the two consumers explicitly.

The new `test_extra_pose_frame` does not catch this: it fabricates a
reference-frame input and compares two consumers supplied with the same
assumption, without passing through the window-search producer. Add a
producer-to-consumer check for each mode, with nonzero shifts and angles.

### R2 — P1: the map's cone test compares shifts from different frames

The dense map fill records `xcell` in each candidate reference's frame.
`select_extras` passes these vectors directly to `prob_insert_cone%within`,
whose translation term assumes a common stored particle frame. Thus the
cone can admit poses beyond its displacement bound or reject nearby ones.

Independent numeric counterexample, C1: equal reference-frame increments
(20,0) px, distinct directions two degrees apart in theta and in-plane
angles 0 and 10 degrees, mask diameter 100 A, sampling 1.3 A/px and cone
resolution 10 A. The geodesic rotation is 10.19754 degrees. The present
test sees zero translation difference and accepts the 8.89903 A angular
displacement. After conversion at each in-plane angle the translation
difference is 4.53210 A; the required bound is 13.43113 A and rejects the
pair. Convert the shifts before both cone eligibility and insertion; the
previous stored shift is common to the particle and cancels in their
difference.

### R3 — P1: weighted insertion treats marginal evidence as a pose loss

With the map, the row is `D_c = -log Z_c(beta)/(beta T)`. The unchanged
`top_posterior_extras` receives that row at temperature T, giving weights
proportional to `Z_c(beta)^(1/beta)`. Those are neither the posterior
weights of the retained cell poses at T nor the integrated evidence at T.
Even at beta=1, a reference's integrated evidence is not the likelihood
of its particular sampled cell.

For example, two retained cells can have equal pointwise likelihood and
prior while their references' marginal evidences differ by a factor of
ten because one reference has more supported cells. At beta=0.3 the
current insertion ratio becomes 2154.43 rather than equal pose weights.
The fill obtains `loss_cell` but discards it. Keep the sampled pose's loss
and prior separately from the marginal reference score, and define and
test the insertion distribution for the map path. The existing window
path's row-loss interpretation cannot be reused unchanged.

### R4 — P2: a small requested disk can produce cells outside its radius

`setup_shiftmaps::build_map` sets `nrho=max(2,floor(radius/step)+1)` and
then `rho(i)=(i-1)*step`. If radius<step, the outer point is step, outside
both the requested radius and potentially `trs`. For box=64,
kstop=4, radius=1 px, it builds magnitudes [0,4] px; kstop=8 gives [0,2].
Such a map is accepted by the sampler without a later radius clamp.
Use a one-radius origin grid or an endpoint inside the requested disk;
test the setup API for radius below, equal to and just above a grid step.
Current kernel tests construct their radii directly and miss this path.

### R5 — P2: the new scope guard breaks the scheduled continuous polish

`prep_strategy3D` rejects `l_shift_map` whenever refine is not `prob`,
before its early return for a Cartesian pass. Both shared-memory and
distributed iteration strategies copy the parent's command into the
scheduled `pose_cont=yes` pass and change refine to `cont` without
removing `shift_sigma`. Thus `refine3D refine=prob shift_sigma=10
pose_cont=yes` reaches the new guard in its own internally scheduled
polish and aborts. Strip the table-only key from that pass, or explicitly
reject the combination before starting the iteration if preserving the
strict posterior draw is the chosen contract. Do not launch workers with
a combination the new guard necessarily refuses.

### R6 — P2: map-mode posterior diagnostics silently write default rows

With `prob_post_diag=yes` the dense fill allocates `post_batch` and writes
it after processing every batch. `process_particle_with_shift` returns
through `process_particle_with_map` before `posterior_diag`, and the map
routine never fills these records. Their default particle index, mode,
temperature and statistics are therefore written as zero. Implement the
map diagnostic, or refuse this combination explicitly. Also distinguish
the map's marginal-evidence row from pose-loss diagnostics and from an
attempt to infer the beta=1 posterior by rescaling an evidence computed
at another beta.

### Coverage and scientific limits

The complete current fast gate passes, but the 2D table tester covers
part I/O, and the 3D table tester uses synthetic rows. Neither exercises
the map fill through weighted insertion or the shift_sigma controller
handoff. The 65 polar-map checks validate the kernel and sampler on their
own; they cannot establish the frame or score semantics of these
downstream consumers. No new registered tests were added by this review.

The 2D screen evaluates its candidate classes over translations. The 3D
screen still starts from only the best `4*projs_ns` references at the
stored shift. It therefore cannot rescue a correct large-offset direction
outside that pool. This is a remaining approximation and validation gap,
not evidence that such an exclusion occurred in the supplied data. The
exhaustive-subset retained-mass test from section 5 remains outstanding,
as does posterior grid-refinement validation. The shift-stability exit
also fires outside the centring stage; the plan already records that its
sampled-motion statistic need not measure posterior convergence.

### Validation observed in this re-review

- Mac Release incremental build (`make -C build -j8`): PASS, no compiler
  warnings observed.
- Full fast gate: PASS, 15/15 in 35.14 s, including registry, descriptor
  and FLEX-DAG checks. Polar alignment: 1390 assertions; pose statistics:
  135 assertions. The sandbox denied the CPU-count sysctl, so CTest ran
  serially; the kernel's explicit four-thread identity check still ran.
- `make -C build install`: PASS.
- Independent numeric counterexamples for R1/R2/R3 and the radius formula
  in R4 reproduce the values above; these are not new compiled integration
  tests. R5/R6 are source-traced findings, not observed runtime failures.
- `simulated_workflow_1jxy`: PASS, 696.67 s. This includes solve2D and both
  shared-memory and two-part distributed solve3D. Truth correlations were
  0.90010 and 0.89485, with truth FSC=0.143 at 4.30345 and 4.62222 A;
  both wrong-truth controls passed. This is the map-off workflow and does
  not exercise the problematic map/topk combinations.
- Final incremental build of the 22:36 local snapshot: PASS. It includes
  the then-unused `nscreen`, `trs_window`, `mout` and mass-outside
  definitions. A hash manifest of 845 Fortran sources remained unchanged
  through that build. One `-Wuninitialized` warning was emitted at the
  untouched `simple_commanders_cavgs.f90:760` (`diams = pack(...)`).
- Final fast gate on that build: PASS, 15/15, 34.97 s. During this gate,
  further live changes appeared in `simple_commanders_prob`,
  `simple_eul_prob_tab2D` and `simple_eul_prob_tab_utils`, wiring the new
  2D mass-outside/screen-pool development. Those ongoing edits are not
  covered by this review or the passing results. An extra install attempt
  was canceled during its consequent rebuild; the earlier install passed.
  The completed workflow uses the earlier build. The 2D screen-rank cost
  instrumentation was inspected; the six findings above remain present.
- No Linux/BOX run, real-data rerun, commit or push by this review.

## 13. Follow-up review, 2026-11-01

**Scope:** `master` at `c980026e5`, including all 19 files dirty at review
start. That batch was committed as `f3cfb911c` during the review; all 845
Fortran source hashes still match the reviewed/built snapshot. Rechecked
the previous six findings against the current
producers, consumers and tests, and reviewed the new coarse-to-fine 2D
centring stage, its drain, convergence and part-file changes. The source
snapshot contains 845 Fortran files. No production source or tests were
edited by this review.

**Stage contract clarified by the maintainer, 2026-11-01:** stage 1 is
deliberately shift-free because introducing shifts there causes class
collapse. Stage 2 owns the current centring budget of **five iterations**
(global iterations 6-10 under the default schedule). S1 concerns initializing
the centring state on entry to stage 2; it does not request shifts in stage 1.
S3 concerns use of that five-iteration budget. Extending centring through
stage 3 is a proposal for consideration, recorded in the living plan's
Open decision 4; it has not been implemented or validated by this review.

### Disposition of the previous findings

- R1/R2: corrected. The map producer converts its cell into particle axes;
  the window producer already returns those axes. The committed and extra
  consumers add the increment directly, and the cone receives a common
  frame. The new producer/consumer test exercises the optimizer itself.
- R3: the original evidence-to-point-loss error is removed, but the new map
  insertion loss omits the absolute prior; see S2 below.
- R4: corrected. A disk below one radial step now uses its own rim. The
  setup API is covered below, at and above the grid step.
- R5: corrected in both execution modes. The continuous child command and
  job description remove `shift_sigma`.
- R6: corrected by an explicit refusal of map-mode posterior diagnostics.

### S1 — P1: centring state is not initialized on entry to stage 2

Both initializers in `simple_refine2D_strategy` reset `shincr` to -1 only
inside the `startit == 1` block (shared line 179, distributed line 335).
However, the solve2D controller numbers iterations across stages:
`set_refine2D_stage_iteration_policy` takes the preceding `endit + 1` and
passes it as the child's `startit`. After the deliberate five-iteration,
shift-free stage 1, stage 2 starts at global iteration 6. Its centring state
must be initialized at that stage boundary. The new named slot initially
contains zero, which
`fill_tab_range` interprets as a particle already done with the disk.

Observed on the 400-particle five-view fixture, with current build binaries:

| route | first centring iteration | particles searched | disk searches |
|---|---:|---:|---:|
| solve2D, shared memory | 6 | 400 | 0 |
| solve2D, two parts | 6 | 200 + 200 | 0 + 0 |
| standalone refine2D centring invocation, startit=1 | 1 | 400 | 400 |

The standalone invocation is a control for the centring initializer,
independent of solve2D's deliberately shift-free stage 1.

The solve2D runs later search the disk only for particles whose preceding
window move exceeds `trs_window`: shared iterations 7/8 had 141/103 disk
particles; distributed had 134/101. Thus this does not disable every later
disk search, but a badly offset particle trapped in a small window move
can avoid the global search throughout. Reset the drain state when the
master enters the centring stage, independently of global iteration
numbering; distinguish that boundary from a deliberate continuation.
Persist the distributed reset before its consumers run.

Validation needed: enter solve2D stage 2 after a nonzero stage-1 `endit`
and assert every first-time sampled particle takes the disk, in both
execution modes. The current convergence and named-slot tests do not
exercise stage initialization.

### S2 — P1: map insertion weights drop the absolute shift prior

The map fill now stores `loss_cell` in `prob_refined_entry%loss`
(`simple_eul_prob_tab`, line 540), and `select_extras` correctly uses that
field instead of marginal evidence. But `shiftmap_draw` defines this as
the data term alone. `top_posterior_extras` receives no prior, so the new
weights are proportional to `exp(-T L_data)` rather than the point
posterior `exp(-T L_data - |u_prev + du|^2/(2 sigma_px^2))`.

For equal data losses, sigma=2 px and absolute stored shifts of 0 and
4 px, the current insertion weights are 0.5/0.5. The posterior weights
are 0.880797/0.119203. Both shifts are inside the 3-sigma disk and can
be inside the insertion cone. Drawing a cell under the prior does not
make its separately renormalized point weights equal to this posterior.
The window path's returned loss already includes its configured prior
factor; the map's raw data loss does not.

Keep the reference evidence separate, as now, but carry a physical
posterior pose loss (including the absolute prior, scaled by the physical
coefficient T) for insertion, or explicitly redefine and validate a
likelihood-only insertion estimator. The present policy calls these
posterior weights and claims one distribution on both paths. This matters
by default now that topk is 3. The example is a closed-form counterexample,
not a measured real-data map failure.

### S3 — P2: class convergence shortens the five-iteration centring budget

`check_conv2D` computes ordinary class-overlap convergence first, then the
new drain code only sets `converged` to true when its streak reaches two.
It never clears an already true class verdict while the disk is still
needed (`simple_convergence`, lines 237-249). Lowering the stage's minits
therefore permits the ordinary 0.8 class-overlap threshold to terminate it
before the documented disk criterion or the iteration budget.

Observed in the same probes: shared stage 2 ended at global iteration 8 with
class overlap 0.825 and 11.50% still needing the disk; distributed ended
at 8 with 6.25% still needing it. Neither printed `DISK EXHAUSTED`.
Stage 2 had five available centring iterations, globally numbered 6-10
(EXTR_LIM_LOCAL=20, PHASES(1)=4). These runs used **three of those five**,
iterations 6-8. Global iteration 10 is the endpoint, not a ten-iteration
centring budget. Both fractions exceed DISK_NEED_FRAC=0.02.

Keep the five-iteration stage-2 budget explicit. If early stopping on the
drain is retained, use that verdict as the centring stage's convergence
condition, leaving the owning strategy's budget/coverage checks in force. Add a
convergence case with stable classes but a non-drained shift field. The
existing test deliberately gives the classes zero overlap and misses it.

### S4 — P2: the hybrid fallback discards a locally improved pose

The translated-box polish calls `minimize` without `xy_in`
(`simple_eul_prob_tab2D`, line 442). That API initializes its acceptance
baseline at zero shift (`simple_pftc_shsrch_grad`, lines 581-599), even
when zero is outside the translated search box. It returns `irot=0` when
the local optimum does not beat that origin baseline. The new caller
interprets this as nothing beating the coarse seed and commits the
unpolished coarse point (lines 453-459).

For full-band losses origin=0.1, coarse seed=0.3, local optimum=0.2,
the local optimization improved the seed but does not beat the origin.
The caller stores 0.3, losing the valid 0.2 result. This can occur when
the half-band coarse map and full-band likelihood disagree. It is also
the case where the original full-band origin row was better than the
coarse point now written over it.

Make the seed and success contract explicit, preserve the best valid pose
(including the origin when it wins), and distinguish invalid optimization
from non-improvement against the origin. If the polish must compare against
the coarse seed instead, expose that baseline deliberately: simply passing
`xy_in` also disables this API's coarse lattice. The same fallback can force
the coarse point when the optimizer actually returns the origin, since an
equal-to-baseline result is not an improvement. Test translated boxes with
the above ordering of three independent losses and with the origin on a bound.
This finding is source-traced, not an observed fixture failure.

### S5 — P2: the cost record omits the drained particles' window searches

The new `grad_shsrch_win` objects perform the window path, but the fill
passes only `grad_shsrch_obj` to `record_fill_profile2D` before killing
both sets (`simple_eul_prob_tab2D`, lines 485-489). Their minimize,
objective, gradient, grid and alternation counters are lost.

The shared solve2D probe's iteration-6 record reports 400 refined entries
but **zero shift searches/objective evaluations/gradient evaluations**:
all work used the uncounted window objects. Every particle had at least
one refinement minimization, plus any seed search. Later mixed iterations
undercount the window subset. Wall timings remain recorded, but these
counters cannot support the cost comparison. Reduce both objects' profiles
before destruction without counting the particle/region totals twice.

### Additional inactive API defect

The old `temp_mass` diagnostic option remains in `shiftmap_draw` although
the new 2D workflow no longer calls it. With `mass_out` and `temp_mass`,
the diagnostic overwrites `s` with the rescaled-temperature normalizer;
the subsequent stochastic draw uses that new normalizer with the original
tempered cell weights. For two equal-area, zero-prior cells with log
weights 0 and -1, temp=1 and temp_mass=2, the first-cell probability
changes from 0.731059 to 0.880797 merely by requesting the diagnostic.
The existing mass-out test uses a greedy draw without `temp_mass`.

There is no current production caller of `temp_mass`, so this does not
explain the observed solve2D behavior. Remove the retired diagnostic API
or preserve the original draw normalizer and test diagnostic invariance
under a fixed RNG seed. This is a latent kernel issue, separate from the
five active findings above.

### Validation and limits

- Existing-tree Mac Release build (`make -C build -j8`): PASS; no compiler
  warnings observed. Fortran source hashes unchanged through build and gate.
- Full fast gate: PASS, 15/15, 34.30 s, including registry, descriptor and
  FLEX-DAG checks. Polar suite: 1369 assertions, including 73 shift-map,
  72 refine3D in-plane-state and 42 2D table-I/O assertions.
- `make -C build install`: PASS.
- Manual probes: shared solve2D, two-part distributed solve2D, and
  standalone refine2D all finish normally. They establish the disk counts,
  convergence behavior and missing profile counts above; they are not
  registered tests and do not establish a map-quality regression rate.
  Files/logs: `/private/tmp/simple_shift_rereview_20261101/`.
- S2/S4 and the inactive API counterexamples are independent numeric/source
  analyses. No new compiled regression tests were added by this review.
- `simulated_workflow_1jxy`: PASS, 756.90 s, including shared-memory and
  two-part distributed solve3D. Truth correlations 0.90010/0.89422; truth
  FSC=0.143 at 4.30345/4.45714 A. Both wrong-truth controls and the
  two-part worker check pass. The output confirms topk=3 with extra poses
  actually used. It exercises the window path, not the map-insertion prior
  in S2. All source hashes still match after completion.
- The new origin-posterior screen quantile does not bound the posterior
  over translated candidates. Its retained mass and coarse-grid convergence
  remain scientific validation gaps. Earlier reported benchmark results
  in the living plan were not rerun or claimed as this review's results.
- No Linux/BOX run, production fix, commit or push by this review.

## 14. Current state and remaining work (2026-11-02)

Reviewed clean `master` at `69c6c9d78c531dea3dfdc515ee89dea2c4935b43`, including the
3D/table consolidation, the parameter default of `shift_sigma=10` and dynamic
centring. Unrelated XML, molecule-fixture and other development changes are
outside this review. Production code and testers were not edited. The current
status/next actions are maintained in the existing living plan, not a new plan.

The original `prob_neff=3` feedback policy has been retired; the current
`prob_neff_readout` reports entropy at the configured draw beta and at beta 1
without changing it. The findings below concern the current shift-search
integration, not a recurrence of that feedback controller.

The dynamic design addresses a real scheduling problem: fractional sampling
can introduce a particle after any fixed centring stage. The master now resets
`shsig=-1` even in the first shift-free stage, and later tables and particle
preparation consult the same per-particle state. Disk particles use the box
mask; the sparse 2D and neighbourhood 3D tables now have disk paths too. The
stored width chooses the domain, rather than becoming a new data-derived prior.
Those ownership and lifecycle choices are sound.

Stage 1 of solve2D remains deliberately shift-free. Its old stage-2 centring
budget was five iterations, globally 6-10; the proposed extension through
stage 3 is now superseded by ruling 19. No finding here requests shifts in
stage 1. Solve3D's first two stages likewise retain their no-shift policy.

### Previous findings at this revision

| Earlier item | Current disposition |
|---|---|
| R1/R2, increment frames and insertion cone | The native-increment producer/consumer contract remains coherent; no recurrence found |
| S1, missing first centring reset | Corrected by the workflow reset before the shift-free stage, with distributed persistence; observed first disk counts below |
| S2, insertion omitted absolute prior | Corrected in dense 3D; neighbourhood insertion still differs, T2 below |
| S3, class convergence bypassed drain | Superseded: the fixed centring stage and drain are retired |
| S4, coarse fallback overwrote a better origin | The caller now retains the origin row when the coarse fallback is worse; the seed/baseline contract still merits the Phase-4 tests below |
| S5, window work absent from profile | Both search-object sets now reach the profile |
| Inactive `temp_mass` normalizer defect | Retired API removed; new physical-moment accumulators are separate from the draw normalizer |

### T1 — P1: continuous polish projects a disk seed into the smaller stage window

`simple_strategy2D_srch.f90:300-301` and
`simple_strategy3D_srch.f90:475-476` still give the joint optimizer bounds
`[-trs,trs]` about zero. The dynamic table can now return a valid disk
increment outside that window in any shifting stage. Both committed poses
and the 3D extras pass through this continuation when continuous refinement
is enabled, as it is by default. The table's wider disk object does not set
the consumer's joint limits.

This also breaks the monotonicity/score contract, rather than merely limiting
the gain from continuous refinement. `minimize_joint` sets `x=xy_in` and
records its baseline from the optimizer's first objective evaluation
(`simple_pftc_shsrch_grad.f90:752-765`). L-BFGS-B projects an out-of-bounds
start before requesting that evaluation (`simple_opt_lbfgsb.f90:137-142`,
`984-997`). Consequently the baseline is the clipped pose. On a valid
non-improvement or a bound-pinned result, the function returns the original
`xy_in` with the clipped pose's score (`simple_pftc_shsrch_grad.f90:832-839`).
An accepted interior result is likewise compared against the clipped pose,
not the authoritative incoming one.

An independent one-dimensional example makes the failure explicit: let the
incoming disk coordinate be 12, `trs=5`, and the loss be
`L(x)=0.1+0.01*(x-12)^2`. The optimizer first evaluates 5, the constrained
minimum, with loss 0.59. The no-improvement return keeps coordinate 12 but
attaches `exp(-0.59)` instead of its correct `exp(-0.1)`. This is a
source-traced counterexample, not a measured real-data failure rate.

Give the continuation a domain that contains the authoritative disk seed
and evaluate that actual incoming pose for acceptance. Preserve a consistent
pose/score on every fallback. Test disk seeds beyond `trs` through the real
2D, 3D and extra-pose consumers, including non-improvement, boundary and
interior outcomes. The current map-initialised-polish fixture does not cover
this table-to-continuous handoff outside the ordinary window.

### T2 — P1: neighbourhood insertion still uses a different posterior loss

The dense table now stores a separate physical pose loss,
`L_data + |u_prev+du|^2/(2*sigma_px^2*T)`, in
`prob_refined_entry%loss` (`simple_eul_prob_tab.f90:584-588`). This is the
correction requested by S2. In contrast, the neighbourhood disk path records
only the polished data distance (`simple_eul_prob_tab_neigh.f90:1547-1551`).
`gather_state_row` copies candidate `%dist` to the insertion input, and
`select_extras` passes it directly to `top_posterior_extras` (lines 947 and
899-901). There is no corresponding physical pose-prior term. The window
case can instead include the legacy relative prior described in T3.

At least for first centring, where the pick prior is explicitly active, the
two routes therefore cannot describe the same insertion posterior. Conditional
on two poses surviving the same insertion cone, equal data losses with
`sigma_px=10` and absolute shift magnitudes 20 and 22 give dense weights
0.603483/0.396517 but neighbourhood disk weights 0.5/0.5. Drawing/selecting a
coarse cell under the prior does not supply this separately normalized
point-pose weight.

Carry the intended physical insertion loss separately from the candidate
selection distance on the neighbourhood route, and compare dense/neighbourhood
weights on the same externally specified poses. The written policy also needs
one explicit distinction: the coarse initializer applies the pick prior only
at first centring, while dense insertion currently applies it on every mapped
iteration. The dynamic note says both "first centring only" and "the pick
prior ... in the 3D insertion's loss". These may be deliberately different
estimators; this review does not assume the dense insertion must drop its
prior after the first iteration. State the initializer, polish and insertion
targets separately before extending the width estimator.

### T3 — P2: mapped sparse/window polish retains the legacy increment prior

The dense mapped fills bypass the old `set_shift_prior` setup. Sparse 2D
(`simple_eul_prob_tab2D.f90:623-624`) and both neighbourhood 3D routes
(`simple_eul_prob_tab_neigh.f90:284-285`, `533-534`) set it whenever
`l_shift_prior` and a shift search are active, including map mode. Their
window branch uses that same object in `minimize_in_cell`; only the disk
branch switches to the separate object without the legacy prior.

The window map explicitly has zero pick-prior weight, but its polish then
optimizes a data term plus the old population-increment penalty. Because
`minimize_in_cell` calls `minimize` without `xy_in`, that relative penalty is
centred on zero increment. Thus switching dense to sparse/neighbourhood, or
switching a particle from disk to window, changes the full-band objective as
well as its domain. It can suppress a window correction or change a
candidate's distance relative to the raw origin row. A coarse fallback is
re-scored by `gen_corr_for_rot_8` without that penalty, so polished and
fallback values can differ in meaning too.

Make the mapped polish objective explicit and apply it consistently. Under
the written data-only window policy, skip/clear this legacy prior in map
mode. If a separate local regularizer is desired, specify and test it on all
table paths instead. Add a dense/sparse/neighbourhood comparison with an
off-zero optimum and `shift_prior=yes`; the existing kernel identity checks
do not exercise these per-particle setup calls. The continuous consumer's
own start-centred prior also needs an explicit disposition in the target
definition above, rather than silently treating it as the absolute pick prior.

### Remaining development: recommended order and qualifications

1. **Integration corrections and dynamic-state tests first.** The proposed
   mixed sentinel/mover/drained fill regression is absent. Test the exact
   `3*shsig > trs` boundary, sentinel persistence through no-shift stages,
   bound inflation, unchanged state on unsampled particles, first arrivals
   after stage 2/3, restart versus fresh-workflow resets, stream injection and
   cropped-record handoff. Cover dense/sparse/neighbourhood and shared/parts
   consumers, including agreement of the preparation mask with the table's
   domain. Direct tests of physical-temperature moments and their frame
   transformations are also missing. The new default makes these completion
   gates more pressing than another search mode or ladder rung.

2. **Laplace width is useful, but the proposed 2D stencil is conditional.**
   The observed near-zero coarse widths mean mass concentrated in one grid
   cell; they do not establish subpixel accuracy. The proposed full-band
   curvature can resolve that scale. However, `gen_corr_for_rot_8` takes an
   integer angle, whereas the proposal measures after continuous in-plane
   polish. Use the actual committed fractional-angle objective; PFTC already
   exposes fractional-angle evaluators. Put the numerical method in the
   owning PFTC/search domain and let the strategies orchestrate it.

   A 2 x 2 shift Hessian at fixed angle gives conditional covariance. The
   current coarse moments marginalize over in-plane rotation. With appreciable
   shift-angle coupling these are different search-width criteria. For a
   joint Hessian partitioned into shift and angle blocks, the marginal shift
   covariance is `(Hss - Hsa * Haa^-1 * Has)^-1`. For the positive-definite
   example `H=[[1,0,.99],[0,1,0],[.99,0,1]]`, the largest conditional shift
   standard deviation is 1; the marginal one is 7.0888. A window of 5 would
   be selected by one rule and rejected by the other. Either marginalize the
   coupled local model or explicitly adopt and calibrate a conditional rule.

   Specify the active physical-temperature loss/prior, stencil scale,
   positive-definiteness and near-flat/negative-curvature fallback. Unresolved
   curvature should keep coverage conservative, not mark a particle centred
   with zero width. Retain the coarse multimodality flag; a sharp wrong mode
   can have narrow local curvature. The proposed largest-eigenvalue width is
   smaller than the current square-root trace (by sqrt(2) for an isotropic
   2D Gaussian), so it changes the convention rather than being uniformly
   more conservative. Validate on known Gaussian/coupled/multimodal examples,
   grid-phase/stencil sweeps and low-SNR fixtures before interpreting width
   histograms on data as calibration.

3. **Measure translation-aware coverage and close Phase 3.** A pool containing
   99% of posterior mass at the stored shift need not contain the reference
   preferred after translation. Coarse width moments describe only evaluated
   hypotheses; their median/90th percentile is not the retained-mass
   measurement. Compare pool results to exhaustive translated evaluations on
   a tractable subset: lost mass, excluded best modes, grid convergence,
   class/view dependence and late first-time samples. Preserve the explicitly
   chosen neighbourhood/subspace scope in that comparison. In particular, a
   new class/direction can have a better, separated peak outside the window
   while the in-window peak is narrow and its move stays interior; neither
   local width nor bound inflation certifies coverage in that case. Test this
   before accepting the dynamic note's claim that class/direction changes
   need no guard. Complete paired
   seeded controls, same-particle per-class shift distributions versus the
   external benchmark, the centring-frame/mask ruling, and re-extraction with
   downstream 3D. The maintainer's reported beta-gal/2D outcomes are evidence
   recorded in the plan; they were not rerun by this review.

4. **Phase 4 is complete as table integration, not as lattice retirement.**
   Every table now gets a map before local polish. Nevertheless,
   `minimize_in_cell` calls `minimize` without `xy_in`, with `coarse_init=true`
   on the table objects, so a 5 x 5 all-angle lattice still runs inside the
   translated box. A subcell seed and a clear best-valid-pose/baseline
   contract remain work if the aim is to replace that lattice. Test origin
   versus coarse versus local full-band loss orderings, then measure its
   fraction of full-iteration cost. The current workflow's first dense 3D
   record has 7743 searches and 193575 grid evaluations, exactly 25 per
   search. Its refinement timer includes both maps and local polish, so it
   does not isolate the lattice's cost. The dated inventory predates the
   neighbourhood integration; classic non-table modes still use their
   existing search. Do not count deliberate no-shift stages as missing map
   integration or expand unsupported modes merely to complete an inventory.

5. **Phase 5 should follow profiling and numerical parity.** Reusing harmonic
   coefficients for Fourier–Bessel-domain polish is plausible, but the
   retained coefficients need bounded per-thread/pair ownership, not storage
   for every particle/reference pair. Both current initializer geometries
   stop at `kscreen=max(klo+1,khi/2)`; reusing their coefficients alone would
   change the full-band polish objective. Budget the additional full-band
   harmonics and truncation requirements on the selected pairs. Specify
   whether Phase 5 replaces the table's shift polish or also the later joint
   continuous-angle solve: coefficients at integer rotations alone do not
   replace the latter's fractional-angle objective. Prove the same normalized objective
   under angular CTF, noise weighting, distinct halves and the chosen prior;
   compare analytic derivatives as well as matched-start optima. Treat the
   polar origin and angular seam explicitly, using stable Cartesian
   derivatives or limiting formulas at zero radius. Verify thread identity
   and unchanged draw/RNG behavior where promised, and measure total time and
   memory. Removing the remaining local lattice may be more valuable than
   accelerating one polish evaluation. Add an intermediate ladder radius
   only when calibrated widths and measured cost justify it.

### Validation and limits at the reviewed snapshot

- Existing-tree Mac Release build, `make -C build -j8`: PASS; no compiler
  warnings observed. Install: PASS. All 822 Fortran source hashes remained
  unchanged through build and gate.
- Fast gate: PASS, 15/15, 36.73 s wall, including registry, descriptor and
  FLEX checks. Polar suite: 1399 assertions, including 72 shift-map and 98
  refine3D in-plane-state assertions. The sandbox denied the gate's CPU-count
  `sysctl`; it ran one CTest worker with its usual OpenMP checks.
- `simulated_workflow_1jxy`: PASS, 878.95 s, including shared-memory and
  two-part distributed solve3D under the current default map path. Truth
  correlations 0.898208/0.900111; truth FSC=0.143 at 4.23051/4.37895 A.
  Both wrong-truth controls and the two-part worker check pass. The output
  confirms topk=3 with extras actually used. All 822 source hashes still
  match after completion. This is one fixed-seed workflow, not a repeated
  random-start quality comparison or a regression test for T1-T3.
- Manual 400-particle five-view fixture, shared/full-sampling solve2D: normal
  stop, 19.3 s. No disk search in stage 1. First shifting iteration: 400/400
  disk searches, all not yet centred. At entry to stage 3, 5/400 still need
  the disk; subsequent stages can reopen it. The first table benchmark now
  counts the searches rather than losing the window object's work.
- Same input, two-part solve2D with `nsample=100`: normal stop, 33.6 s;
  worker completion and both part records observed. First shifting iteration:
  100/100 disk searches, all not yet centred. Stage-3 entry still has 33
  sampled particles marked not yet centred, and stage-6 entry has 12; both
  are included in the disk set. This establishes later-stage dispatch, not correct shifts
  against truth. These are manual observations, not registered regressions.
  The uncontrolled random starts finish with three populated classes in the
  full run and one in the partial run, so they do not establish class-quality
  parity or a regression rate. Logs and independent arithmetic:
  `/private/tmp/simple_shift_review_20261102/`.
- T1-T3 are source-traced integration findings; the numeric examples are
  independent arithmetic, not new compiled regression tests or attributed
  real-data map failures. No new real-data or Linux/BOX run, production fix,
  commit or push by this review.

## 15. End-to-end mathematical closeout audit (2026-11-03)

The maintainer requested the current closing diff and the entire methodology,
with performance excluded. The initial snapshot was `15a7b0774` plus the
closing diff, including the classic searches, cell-start polish, rulings
21-22 and PCG weighted insertion. During review that same source was committed
as `22703d0ec`; all 823 Fortran hashes still match the initial snapshot. The
intervening multi-state posterior proposal is separate, unimplemented work
and is outside this audit. This review changes no production code or testers.

The closing changes resolve the previous T1-T3 findings under their adopted
contracts: the joint domain contains the incoming disk pose; both 3D tables
use `insert_pose_loss`; and the legacy increment prior is disabled throughout
map-mode polish. Ruling 22 deliberately changes the prior to an increment
about the current origin. This audit does not repeat the old absolute-pick
requirement. The explicit zero keys also correct the former controls, and the
PCG extension has the right weighted least-squares algebra.

Four live issues remain. The first two break numerical invariants across the
pipeline. These findings mostly predate the closing diff; they matter to an
end-to-end closeout even though a diff-only review would miss them.

### M1 — P1: the map, FFT score row and continuous polish use different losses

The Fourier-Bessel map and direct fixed-pose residual retain the angular
Nyquist component. The Euclidean FFT row drops it:
`simple_polarft_corr.f90:209,221` stop at `pftsz`, and
`gen_euclid_crvec` passes `pftsz`, rather than `pftsz+1`, to the shell sum
(lines 292,296). The continuous Euclidean evaluator also stops at `pftsz`
(line 1302). The existing cc continuous evaluator already includes the bin.

The continuous-angle development deliberately preserved the legacy FFT
semantics (`continuous_inplane_rotation_development_plan.md`, Phase 0 and
the interpolation contract). This finding does not treat that preservation
as an accidental change. It concerns the combined method's use of the
full-bin map/direct residual alongside that baseline without one common
energy convention or a bound on their discrepancy in likelihood units.

This reaches production through origin/pool scoring, the alternating
polisher's all-angle callback and returned score, table losses/insertion
weights, and the continuous consumer. The alternating optimizer's fixed-angle
objective/gradient uses the direct residual, so even that polish combines
different objective conventions. This prevents a general same-pose likelihood
identity and a monotonicity claim against the map/direct objective.

An isolated compiled probe against the current library used the closed-form
half-ring pattern `R(p)=X(p)=(-1)^(p-1)`, with an even half-ring length of 36,
unit CTF and noise variance, and zero shift. The expected loss is zero at
rotation 1 and four at rotation 2. Observed:

| Pose | Direct residual | Shift map | FFT row | Continuous evaluator |
|---|---:|---:|---:|---:|
| Rotation 1, zero shift | 0 | 0 | 2 | 2 |
| Rotation 2, zero shift | 4 | 4 | 2 | — |

The effect is not confined to a pure endpoint pattern. With
`R=X=1+0.05*(-1)^(p-1)`, shells 2-31, box 64 and unit variance, the true loss
is zero while the row/continuous loss is 0.00498736. At the calculator's own
physical temperature that is **3.88759 nats**, despite the small normalized
loss discrepancy. The map gives approximately zero (single-precision error
of -3.58e-7). These are controlled polar-array counterexamples, not a claim
about the endpoint content of the maintainer's data or its failure rate.

Use one angular-spectrum convention in map, score rows, fixed-angle
objectives/gradients, continuous evaluator and pose weights. Including the
endpoint consistently is one option; a deliberately projected objective must
project all terms consistently and be specified as such. Add same-pose tests
with nonzero endpoint content, at zero/nonzero shifts and grid/fractional
angles, checking energy differences in nats and derivatives. The existing
kernel test compares against the direct residual, and the older bank test
explicitly tolerates the dropped bin; neither pins this whole invariant.

### M2 — P1: record availability changes a disk particle's reconstruction observation

`prepimg4align` normalizes a disk particle against `build%lmsk_box`.
`polarize_batch_particles3D` keeps that normalized transform's unshifted crop
as the reconstruction record (`simple_matcher_ptcl_batch.f90:290-296`).
`get_rec_observation` retrieves it without renormalizing. With no usable
record, both reconstruction backends instead call `prep_rec_observation`
with `build%lmsk`, the ordinary mask (`simple_matcher_ptcl_io.f90:676` for
PCG; `simple_matcher_3Drec.f90:568` for gridding).

The masks determine which background pixels estimate the mean and variance,
so this changes the observation's scale and offset, not just its spatial
support. It is particularly relevant to an off-centred particle whose signal
falls into the ordinary mask's supposed background. Record availability is
controlled by the memory budget and partition size. Consequently otherwise
identical scientific settings can produce different reconstruction inputs.
The normalization factor depends on each particle's content, so a single
final-map rescaling cannot generally remove the difference.

A compiled probe used a deterministic sinusoidal background plus an
off-centred Gaussian signal in a 64-pixel image, cropped to 32. It exercised
the actual record writer/reader and fallback preparation helpers. The
ordinary/box background standard deviations were 2.32706 and 1.00209.
The reread observation's demeaned amplitude was **0.430623** times the record's,
matching the independent ratio of those standard deviations; the relative
L2 difference was 0.591298. A control using the same box mask in the fallback
agreed with the record exactly. This is a controlled observation test, not a
measured map-quality difference.

Give reconstruction records and rereads one normalization policy. Either
produce records under the reconstruction's invariant policy or carry enough
normalization information to reproduce the same observation on reread.
Selecting a mask from the current `shsig` after assignment alone is
insufficient: that state may already have changed since the record was made.
Test disk/window observations with records available, absent and partially
missing, through both backends and across cropping. This is a numerical
correctness issue, independent of whether the buffer is faster.

### M3 — P2: the harmonic convolution indexes outside the ring for wider searches

The convolution needs the spectral index `(m-n) modulo nrots` for every
Bessel order n. `gen_shiftmap` implements this with two slices that handle
only `n <= nrots` (`simple_polarft_shiftmap.f90:274-275`). The public radius
is bounded by the box, but neither setup nor construction enforces that
harmonic-order limit. For example, box 128, kstop 32, rho 39 and pftsz 36
require order 87 while the ring has 72 positions.

Reproduced in the current Debug executable on an imported 248-pixel,
1.3 A/px stack, with a standalone `refine2D` invocation:

```text
ncls=5 cls_init=ptcl mskdiam=20 lp=2.6 box_crop=248
nthr=1 mkdir=no maxits=3 refine=prob trs=4 shift_sigma=26
```

At its first disk search it aborts at line 274 with index 0 outside
`1:72`. This is a wider-prior/small-mask boundary case, not an observed
failure of default sigma 10 on the supplied benchmark sets.

Wrap the spectral index for arbitrary n, or reject an explicitly documented
unsupported range before entering the fill. Do not silently cap the Bessel
order: the discarded coefficients can still matter. Angular sampling of
the translated object is a separate accuracy condition; exact evaluation
of the discrete ring does not remove its sampling/aliasing limits. Add a
case with required order greater than one ring length against a direct sum.

### M4 — P2: an excluded pose can erase all eligible insertion weights

`top_posterior_extras` takes `dmin` over every valid row, including entries
outside the eligibility mask (`simple_eul_prob_tab_utils.f90:1129-1132`).
It then zeros weights beyond 80 nats relative to that minimum before
renormalizing the kept set. If an excluded pose is much better, every
eligible weight can become zero. The fallback then gives the committed
pose weight one, even when an eligible extra has comparable conditional
posterior weight.

The current-library probe used distances `[0,0.9,0.901]`, temperature 100,
mask `[false,true,true]`, committed index 2 and topk 2. It returns one pose
with weight 1. Removing only the excluded first row returns two poses with
weights **0.524980/0.475020**, as the closed-form conditional distribution
requires. Thus changing a hypothesis that cannot be inserted changes the
distribution over an unchanged eligible set. The dense and neighbourhood
routes can both reach this helper, including rows with unrefined distances.

Use the minimum over the committed pose and eligible extras for this local
normalization. Retain a regression that changes an excluded row without
changing any kept-set weights. This is a numerical tail case; no frequency
of occurrence on real runs is established by this review.

### What is mathematically sound

Write P for the half-ring length, N=2P, and
`W_X=sum_k (k/sigma_k^2) sum_p |X_kp|^2`. On the direct objective,

```text
L(theta,s) = sum_k (k/sigma_k^2) sum_p |X_kp - CTF_kp R_theta,s,kp|^2 / W_X
T = pi W_X / (2P)
T L = sum_k (pi k/P) sum_p |residual_kp|^2 / (2 sigma_k^2).
```

The temperature therefore restores the Gaussian residual energy under the
stated polar quadrature and per-real-component noise convention. Masking
and polar sampling still make this a diagonal-noise approximation; that is
distinct from a normalization or factor-of-two error.

For the full Friedel ring, the cross term follows by expanding the phase
ramp with Jacobi-Anger and using the DFT of the reference and CTF-modulated
particle. The shifted spectral index is `m-n`; the full-ring cross sum is
twice the real half-ring sum, giving the half-weight on order zero and one
copy of each positive order. Reference power correctly retains its
rotation-dependent angular CTF term. The kernel respects the half-specific
reference and sigma2 spectrum. An independent direct-DFT/Bessel calculation
with an angular CTF checked 200 poses: maximum cross error relative to the
Cauchy-Schwarz scale 1.79e-10, normalized-loss error 1.61e-10 and energy
error 5.10e-8 nats. This confirms the identity within truncation, inside the
supported index range; M1 concerns its consumers.

Polar cell areas are the correct measure for integration, stochastic draws
and the coarse moments. The origin's area is divided among its duplicate
direction entries. The new isotropic increment prior is invariant under
rotation into particle axes and under a common change of the current origin.
Its weight `1/(2*(shift_sigma/smpd)^2)` has the correct pixel units.
`insert_pose_loss = L + |du|^2/(2*sigma_px^2*T)` gives precisely
`exp(-T L) exp(-|du|^2/(2 sigma_px^2))`; the prior is not accidentally
tempered in those physical insertion weights. Turning the old increment
prior off under the key avoids adding another penalty during polish.

The increment is converted from optimizer/reference axes to particle axes
once, then added to the stored shift. The revised joint domain includes
that incoming pose, fixing T1. The cell-start polisher now retains the better
of its cell and local result and compares against the stored-shift row,
fixing the earlier fallback contract independently of M1's objective issue.
Stage 1 of solve2D and solve3D's first two stages remain shift-free by design.

The cone's rotation is minimized over point-group equivalents, and
`r*alpha + |delta_shift| <= resolution` is a sufficient rim-displacement
bound: `2r sin(alpha/2) <= r alpha`, followed by the triangle inequality.
Its physical scale is heuristic, but the stated geometric bound is sound.
The table conditions on the committed state and renormalizes after the
cone/topk/weight-floor truncation, subject to M4.

For fixed pose responsibilities omega_ij, weighted reconstruction minimizes

```text
Q(V) = sum_i b_i sum_j omega_ij ||Sigma_i^-1/2 (y_i - A_ij V)||^2,
B    = sum_i,j b_i omega_ij A_ij* Sigma_i^-1 y_i,
H    = sum_i,j b_i omega_ij A_ij* Sigma_i^-1 A_ij.
```

There are no cross-pose terms in this expected-residual objective. The new
PCG entries implement it correctly: the same unshifted observation at each
pose, that pose's transfer/shift, and its responsibility applied once in
both B and density. They do not square the responsibility or treat poses
as additional population mass. The committed shift is already cropped;
extra stored shifts are multiplied by `box_crop/box` once. This preserves
the phase `2*pi*h*shift/box` and agrees with the gridding pose convention.
Since each particle's weights sum to one, fractional/trailing population
accounting continues to count particles. M2 concerns the observation fed
into this otherwise-correct algebra.

### Statistical interpretation and deliberate approximations

These points qualify the method; they are not additional code defects or
requests to undo the maintainer's deferred work.

- **The initializer selects a cell mass.** Its greedy score is
  `-beta*T*L - prior + log(cell_area)`. That is correct for the most probable
  discrete polar cell. A Cartesian-density MAP would omit the area term.
  The distinction matters for a deterministic seed: an isotropic Gaussian
  with standard deviation 3 px on radial steps of 4 px has its density
  mode at zero but its largest per-direction polar-cell mass at radius 4.
  This is an initializer heuristic; it is not a proof of a continuous
  global minimum. A density-MAP seed is a possible modeling choice, rather
  than a correction required for stochastic sampling.
- **The pool is conditional on the current shift.** Retaining 99% of its
  posterior there does not retain 99% after translation. Full-band polish
  then replaces in-plane integration by a locally optimized point per
  selected reference. The direction draw therefore approximates a
  distribution over those constructed candidates. It is not the earlier
  fully marginalized joint pose draw. Mode heights alone omit basin volumes
  (a Laplace marginal would contain a Hessian-determinant factor).
- **Continuous polish transports frozen weights.** Extras and committed
  poses take the same continuous route, which resolves their old asymmetry.
  But insertion weights are computed at the table poses and retained after
  their individual data-only refinements. This is a transported finite
  mixture; it need not equal the physical point posterior evaluated at the
  final poses. If that latter interpretation is required, re-score data and
  increment-prior terms after polish and renormalize. Frozen table weights
  are explicitly the current ruling, so this is a modeling qualification.
- **The moving prior is an iterative regularizer.** Ruling 22 makes it a
  prior conditional on the current origin, with an empirical picking-width
  scale. It is equivariant to a consistent origin change and avoids the old
  pull back from 2D-aligned origins. It does not define one fixed posterior
  about the original pick across all iterations. Initializer, data-only
  polish and physical insertion are intentionally different targets.
- **Coarse width controls coverage heuristically.** The square-root trace
  is the RMS shift deviation around the coarse posterior's mean, after
  rotating each cell into common particle axes. It measures only evaluated
  cells of the selected reference. One occupied coarse cell can give zero
  width without subpixel certainty; an excluded or separated outside-window
  mode is invisible. A narrow interior mode after a class/direction change
  need not trigger either width or bound inflation. The Laplace proposal
  can resolve local curvature but cannot certify absent global mass. Its
  conditional-versus-angle-marginal covariance choice remains relevant.
  Three sigmas also needs a dimensional convention: an isotropic 2D
  Gaussian has 98.89% mass inside a radius of three per-axis sigmas.

The dynamic-state tests, stochastic 2D comparison and Laplace estimator were
explicitly deferred by the maintainer. That scope decision is respected.
For mathematical reporting, describe the current method as adaptive
coarse-to-fine pose optimization with probabilistic assignment and a local
weighted reconstruction mixture, with the approximations above stated.

### Validation and evidence

- Existing-tree Mac gfortran-16 **Debug** build: PASS, no compiler warnings
  observed; install PASS. Fast gate: PASS, 15/15, 89.5 s, warning-only budget
  overrun. Polar suite: 1429 assertions; reconstruction: 546, including the
  new PCG entry checks. The CPU-count sysctl is sandbox-denied; the gate
  uses one CTest worker with its usual OpenMP team checks.
- PCG fractional-update gate: PASS, 6.90 s. It validates population/chain
  behavior; it is not a new topk-with-CTF-and-cropping oracle.
- The full Debug `simulated_workflow_1jxy` was started in the background
  alongside that gate and is still running at this review's handoff. No
  current workflow pass or map-quality result is claimed. Its live CTest
  output is `/private/tmp/simple_polar_math_review_20261103/workflows.log`;
  the metrics are `build/test_runs/simulated_workflow_1jxy/metrics.tsv`.
- Independent Python/standard-library DFT/Bessel, Gaussian, insertion and
  weighted-normal-equation calculations are in
  `/private/tmp/simple_polar_math_review_20261103/math_oracles.py` and JSON.
- Compiled scratch probes linked the unchanged current library, in an
  isolated local clone's worktree, without adding a CMake target or touching
  repository testers. `math_probe.f90`, `math_probe_weak_nyquist.f90` and
  `observation_probe.f90`, with their logs, are under
  `/private/tmp/simple_polar_math_review_20261103/probe_checkout/`.
  They reproduce M1, M2 and M4 through production APIs. The standalone
  commander reproduction of M3 is under `harmonic_bounds/` in the same root.
- No new real-data or Linux/BOX run, production fix, commit or push by this
  review. The maintainer's reported 10/10 PorF result and membrane-protein
  improvements are not this review's measurements. These counterexamples
  establish mathematical/numerical failures, not their frequency on those
  data sets.
