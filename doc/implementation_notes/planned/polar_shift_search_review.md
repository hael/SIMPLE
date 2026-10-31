# Probabilistic in-plane search: mathematical and performance review

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
