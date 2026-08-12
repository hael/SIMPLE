# PCG reconstruction production-readiness plan

## Status

Approved implementation plan. Phase 0 and the shared-memory `reconstruct3D`
kernel backend are implemented; runtime validation is in progress. Implemented behavior is defined in
`doc/policies/reconstruct3D_pcg_policy.md`; this note tracks the remaining
production gap, validation plan, and later research.

## 1. Bottom line

The fixed-pose PCG command is a credible experiment, not a production backend.

1. **Production PCG must use `pcgop=kernel`.** Matrix-free cannot scale to real
   experimental workflows.
2. **Matrix-free remains the unit-test oracle.** It validates kernel numerics on
   small fixtures; it is never a workflow backend or fallback.
3. **Productionize the unregularized fixed-pose route first.** Priors and pose
   refinement must not hide base operator or integration errors.

The first release should be opt-in. Gridding remains default until kernel PCG
shows a reproducible advantage at matched compute.

Review decisions already settled:

- the merged map combines the independently solved halfmaps; there is no
  additional full-data PCG solve;
- the public selector is `rec_backend=gridding|pcg`;
- the first integration milestone is `reconstruct3D`; `refine3D` follows;
- kernel failure never retries with matrix-free because that route is too slow;
- trailing remains in the accumulator domain.

## 2. Production gap

| Area | Required outcome |
| --- | --- |
| regression | Linux/macOS CI, debug/release, thread-count baselines |
| command I/O | deterministic project/stack fixture, including a short final batch |
| gold standard | independent even/odd kernel solves and normal FSC artifacts |
| kernel validation | matrix-free unit oracles plus kernel-vs-gridding real halfmaps |
| distribution | reduce additive statistics before kernel/preconditioner finalization |
| workflow parity | crop/resolution, states, weights, fractional updates, box growth |
| assembly | standard `volassemble`, postprocessing, project, and restart artifacts |
| operations | typed validation, structured outcomes, atomic publication, resource diagnostics |

Continuous pose/shift refinement remains separate in
`doc/implementation_notes/continuous_3D_refinement_on_pcg_operator.md`.

## 3. Production rules

1. **Keep ownership boundaries.** PCG owns the fixed-pose operator, additive
   data statistics, kernel/preconditioner, solve, and diagnostics. Existing
   workflow code owns particle selection, halves/states, FSC, trailing,
   filtering, postprocessing, project output, and iteration handoff.

2. **Reduce raw sufficient statistics.** Each `(state,half,part)` writes:

   ```text
   B = weighted RHS accumulator
   D = |CTF|^2/sigma2 Gram/sampling-density accumulator
   ```

   Reduce `(B,D)` before deapodization, support, flooring, kernel finalization,
   `lambda`, or solving. Never average finalized partition kernels.

3. **Use restart-safe artifacts.** A versioned manifest written last records
   state, half, part, geometry, particle count, symmetry, objective, and
   operator version. Missing or mismatched components fail before solving.
   Initial box-256 compact target: about 269 MB for `D_folded` plus 67 MB for
   `B_raw` per artifact. Do not persist `Khat`, preconditioners, or particle
   planes.

4. **Preserve half independence.** Even and odd inputs, statistics, solves, and
   halfmaps remain separate through FSC. The merged volume is formed by
   combining the two dense halfmaps after FSC. Do not run a third full-data PCG
   solve.

5. **Blend trailing statistics, not volumes.** Store unregularized `(B,D)`.
   Mirror the current gridding `volassemble` contract: keep separate even/odd
   RHS and sampling/Gram components as one manifest-validated chain, seed it at
   full sampling mass, scale current partials by `u/f`, scale the previous chain
   by `1-u`, and restore halfmaps/FSC only after the blend. Apply `lambda` only
   after blending. A previous solution may warm start PCG but is not
   authoritative.

6. **Expose one production backend.** A selector such as
   `rec_backend=gridding|pcg` maps `pcg` to kernel. Do not expose matrix-free in
   production workflow configuration.

7. **Fail rather than switch operators.** Kernel curvature or non-finite errors
   return a structured failure. Never retry with matrix-free; that would hide a
   defect, is too slow for the target workflow, and cannot work after
   distributed statistic reduction.

8. **Reuse the existing particle cache, not a PCG plane cache.** A standalone
   `reconstruct3D` reads every selected source particle once, so constructing a
   new cache for that pass would add I/O. When PCG is later used by `refine3D`,
   consume the existing validated downscaled particle cache through the same
   explicit, uniform cache contract as matcher reconstruction. In both cases,
   discard batches after `(B,D)` accumulation; PCG iterations reread no
   particles.

9. **Treat symmetry support as a measured cost contract.** Coordinate
   replication is mathematically general but accumulation cost grows roughly
   with group order (for example C2: 2, D2: 4, icosahedral: 60 operators).
   Validate C1 and C2/D2 first, add a high-order performance case, and report
   setup progress/work estimates. Do not advertise a high-order group as
   production-supported until it meets the declared runtime budget.

## 4. Implementation sequence

| Phase | Work | Exit condition |
| --- | --- | --- |
| 0. Freeze baseline | register existing test in CI; establish kernel-oracle tolerances; remeasure box-256 time/RSS; fix help drift | nine stages pass reproducibly with retained logs |
| 1. Harden command | typed inputs; structured stop reasons; on-disk kernel fixture; negative tests; atomic output | real command I/O is covered and failures leave no valid-looking output |
| 2. Shared `reconstruct3D` | separate even/odd solves and FSC; combine halfmaps; standard outputs; C1/C2/D2 gates; add `rec_backend` dispatch | **Implemented; runtime fixture and C1/C2/D2 acceptance remain** |
| 3. Distributed `reconstruct3D` | versioned `(B,D)` reduction; sequential state/half solves; restart/provenance tests; standard `volassemble` handoff | N-part and monolithic fixed-iteration results agree and the reconstruct3D milestone is complete |
| 4. Refinement semantics | states, crop/resolution, weights on both `B` and `D`, fractional updates, box growth, accumulator trailing, lambda scaling, existing cache consumption | ordering, partitioning, updates, restarts, and cached/uncached runs preserve the objective |
| 5. Integrate `refine3D` | route shared/distributed refinement reconstruction through the selector while retaining normal assembly/project handoff | supported opt-in refine3D workflows pass production-usable gates |
| 6. Promotion | benchmark matched kernel PCG and gridding | promote only with a material reproducible advantage |

The current shared implementation intentionally rejects box cropping,
fractional/trailing reconstruction, `projrec=yes`, `conical_fsc=yes`
regularization, and all distributed part execution. The standard cFAR
diagnostic is supported and uses the same mask as the FSC. These guards keep
the first runtime validation on the fixed-pose full-box kernel path; they are
removed only with the phase-specific equivalence tests above.

Trailing should first be tested in a standalone harness covering `uf=1`,
convergence, round ordering, warm starts, box growth, artifact round trips, and
restart equivalence.

## 5. Validation and release gates

### 5.1 Matrix-free oracle tests

Both operators receive identical RHS, metadata, symmetry, masks, limits, and
probe vectors. Every kernel-affecting change compares:

1. `H_kernel p` with `H_matrixfree p`, including interior and maximum error;
2. bilinear symmetry and nonnegative energy;
3. operator scale and reconstruction amplitude;
4. fixed-iteration solutions and a final true residual evaluated matrix-free;
5. recovery from independently generated known-truth observations.

Residual traces alone cannot detect a uniform operator-scale error. Coverage
must grow to include shifts, astigmatic CTF, CC and sigma weighting,
unequal/zero weights, support, resolution limits, multiple boxes, C1, C2/D2,
and one non-lattice-exact group.

### 5.2 Workflow and real-data tests

The on-disk workflow fixture always runs kernel. It covers project discovery,
image preprocessing, CTF/shift/sigma metadata, halves, states, short batching,
manifests, outputs, FSC shape, and diagnostics.

For the first milestone, diagnostics remain a human-readable key/value record.
They include solver stop reason, iteration count, initial/final true residual,
final relative update, state/half and particle counts, kernel/operator version,
objective and regularization settings, timing by phase, peak RSS when available,
and input/output artifact provenance. Add JSON or STAR only when a workflow or
test requires machine ingestion.

Real-data tests compare independent kernel halfmaps with gridding for:

- heterogeneous CTF and moderate/large particle count;
- low particle count;
- axial-view-heavy C2/D2 data;
- a high-order symmetry performance/limit case.

Freeze particles, metadata, geometry, postprocessing, compute layout, and
hardware. Compare FSC, map differences/amplitude, anisotropy, time by phase,
peak RSS, disk use, solver outcome, and repeatability.

### 5.3 Production-usable gate

All must pass:

- deterministic kernel-vs-matrix-free unit tests;
- independent halfmaps/FSC and real-data kernel-vs-gridding acceptance;
- monolithic/partitioned and shared/distributed parity;
- supported state, crop, resolution, weight, fractional, and trailing paths;
- standard assembly/project artifacts and safe restart/failure behavior;
- explicit convergence, iteration-limit, and failure outcomes;
- declared memory/disk budgets;
- unchanged default gridding behavior.

## 6. Performance and memory rules

- One logical particle read per accumulation phase: direct source-stack reads
  for one-off `reconstruct3D`, or reads from the already validated workflow
  cache for cache-enabled `refine3D`.
- Particle-image residency bounded by `MAXIMGBATCHSZ`.
- No particle-plane cache: kernel iterations are data-free after `(B,D)`.
- Process one state/half at a time; reuse FFT plans where safe.
- Do not overlap the largest accumulation and solve scratch allocations.
- Measure I/O, preprocessing, scatters, reduction, kernel setup, solves,
  assembly, and peak worker/master RSS.
- The shared implementation reports metadata/operator preparation, particle
  I/O/preprocessing, fused `(B,D)` accumulation, accumulator finalization,
  solve, map output, and FSC/cFAR time independently for each state/half. Use
  these timings—not total wall time alone—to choose the next optimization.
- Production streaming updates `B` and `D` in one interpolation traversal and
  one persistent OpenMP region per batch. Preserve the monolithic two-path
  implementation as the numerical oracle.
- Treat the historical box-256 peak of about 3 GB as a rebaseline target, not
  current evidence.
- Optimize sigma2 storage only after measurement and without changing weights.

## 7. Later research

Do not add priors until the unregularized kernel route passes Phase 2.

| Idea | Priority | Recommendation |
| --- | --- | --- |
| Tikhonov lambda scaling | P0 | define against effective data mass before distribution |
| solvent-flatness quadratic prior | P1 | first new-prior experiment after independent halfmaps |
| soft state weights | P1 / Phase 4 | land after hard-state parity; weight `B` and `D` |
| smoothness quadratic prior | P1/P2 | controlled solvent-flatness comparator |
| symmetry projection/permutation | P2 | pursue after distributed profiling |
| total variation | P3 | separate nonlinear/proximal solver project |
| non-negativity | P3 | defer until density sign/baseline is validated |
| wavelet L1 | P3 | defer behind TV |
| joint pose/shift refinement | P3 | separate program, outside productionization |

Any prior must preserve half independence, positive-semidefinite CG structure,
conditioning, data-mass scaling, and separate reporting of fit and prior energy.
Evaluate it against an unregularized independent-half control.

## 8. Remaining decisions

1. Define “material advantage” before promotion benchmarks.
2. Set the measured runtime/memory threshold for advertising high-order
   symmetry groups as production-supported.

## 9. First patch series status

Completed: typed selector/UI, backend and execution-mode dispatch, CI gate,
strengthened kernel-oracle checks, structured solver outcomes, and the shared
independent-halfmap implementation. Still required: the deterministic on-disk
kernel fixture and C1/C2/D2 runtime acceptance.

This series does not change `refine3D`, enable distributed PCG, or make PCG the
default.
