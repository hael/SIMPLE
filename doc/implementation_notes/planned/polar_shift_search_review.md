# Probabilistic in-plane search: mathematical and performance review

**Latest review:** the 2026-11-01 follow-up is in section 13. Section 12
retains the previous integration review; sections 1-11 retain the original review.

**Reviewed:** 2026-10-30, `master` at `a2c1b1f86`.
**Scope:** commits `519c67d2a` and `a2c1b1f86`, the proposed integration into solve2D/solve3D, and the relevant existing probability-table and reconstruction contracts.
**Status:** review findings and recommendations; no implementation changes made or recommendations adopted by this review.

The living implementation plan remains [polar_shift_search.md](polar_shift_search.md). This document is the separately requested review record, not a second implementation plan. The related history is in [pose_marginalization.md](pose_marginalization.md), especially sections 11, 20.8–20.13 and 21.

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
