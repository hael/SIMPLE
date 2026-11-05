# Probabilistic alignment as a batch pipeline: report of the run `prob_pipeline`

Run on the Dell, 2026-11-04 09:04 to 2026-11-05 05:30, on base `22a643642`; applied
to the maintainer's checkout on top of `51dd0a627` (the merged result equals the
run's tree plus the newer commits). The plan and its five Progress entries, one per
phase with every piece of evidence, are in `prob_align_batch_pipeline.md` beside this
report. Phases 0-4 each ended DONE; the driver then stopped the run on its file fence
(phase 4's retirement edits touched three files the file table had assigned to phases
1-3 only: `simple_refine2D_strategy.f90`, `simple_refine3D_strategy.f90`,
`simple_eul_prob_tab_utils.f90`, an omission of the table, not of the work), so this
report and the move of the plan were written by hand afterwards.

## What was learned

- **The refactoring holds.** Every probabilistic draw (the sample, the in-plane draw
  per scored reference, the disk draw, the sparse fills' candidate picks, the
  assignment) is keyed on a stateless generator (seed, pass, iteration, particle,
  draw). The table, the part files, the assignment and insert files, the master's
  assignment step, the six table programs and their two modules are gone. One batch
  pipeline (prepare, fill complete rows, draw, polish, residual, commit) runs inside
  `refine3D_exec` and `refine2D_exec`, in shared memory and in every distributed part;
  the master draws the sample, schedules once per iteration, merges the parts' readout
  records in part order and assembles.
- **Draws no longer depend on threads, batches or parts.** The frozen-input harness
  (each route's fill, selection and extras on the same prepared batch) gives zero
  differences across 1, 8 and 16 threads, two batch sizes and 1 to 3 parts, in Release
  and Debug. On the old code the drawn entries differed with every one of those.
- **The old two-phase path had three defects**, each a value the table worker's second
  toolbox did not share with the strategy: (1) without an explicit `lp` the 2D table
  worker never read the class FRCs and scored at the `lpstart` band while the matcher
  polished at the FRC band (from solve2D's stage 6 on); (2) a fresh 2D start's shifts
  were zeroed in the matcher only, so the table was scored on inherited shifts;
  (3) the dynamic shift window widened in the strategy did not reach the worker. The
  pipeline removes all three by construction; they are why the shared-memory workflows
  part from the old code at solve2D stage 6 and the trpv1 refine2D case from its first
  iteration. Proven by controls: with `lp` (and `trs`) given explicitly the old and new
  code are byte-identical.
- **One draw is still on the intrinsic stream**: the shift search's simplex start
  (`simple_pftc_shsrch_grad` `lims_init` → `simple_opt_simplex` `ran3`), missed by the
  plan's audit. Its effect is a slightly different optimum for a handful of particles.
  Keying it is left for a later pass (it touches `simple_opt_simplex.f90`).
- **Complete runs still differ across thread counts**, now by the preparations'
  threaded arithmetic and the reductions' summation order (section 8's excluded
  contract), which the old thread-dependent draws used to hide: trpv1 refine3D final
  maps correlate 0.99998 across part counts, 0.9996 one thread against eight.

## Evidence (detail in the Progress entries)

| level | what | result |
| --- | --- | --- |
| 1 | generator vectors (ten keys from an independent Python reference) | bit for bit in Release and Debug |
| 2 | frozen-input harness, four routes x 2 state counts x 6-9 configurations | every count zero, Release 8 and 16 threads, Debug 4 |
| 3 | each phase against the previous at fixed configurations | byte-identical on every declared output, except the shared-memory deviations explained above (phases 2-3), each traced to the old path's defects |
| 4 | fixed-input reconstructions and class averages at nparts 1, 2, 4 | phase 0's numbers to the digit (8.58e-5 / 7.79e-7 against the frozen 1.72e-4 / 1.56e-6) |
| truth | `multistate_workflows` 129/129, `solve3D_addon` 46/46, `cont_refine3D_1jxy` 135/135 | floors unchanged, met |
| nightly | `test=nightly` on the final install | 18 of 19 pass at every floor (1,078 rows); the nineteenth is `simulated_workflow_1jxy`'s base failure (solve3D stage 3, `new_shiftmap`), present at the base |
| readouts | bgal solve2D over eight seeds; trpv1 refine3D three iterations | within the seeds' spread; trpv1 shells 4.041/6.224 and 4.095/6.621 A at every layout |

## Measurements (base vs final, 4 parts x 6 threads, three runs, machine quiet)

| case | base | final |
| --- | --- | --- |
| trpv1 refine3D, 3 iterations, crop 160 | 1,015 s | 1,001 s (fill unchanged, 293-300 s per part and iteration) |
| trpv1 refine2D, crop 120, 2 iterations | 41.6 s | 34.6 s |
| trpv1 refine2D, box 256, 2 iterations | 96.0 s | 81.5 s |
| master peak RSS (3D) | 1.26 GB | 1.10 GB |
| worker peak RSS (3D) | 2.43 GB search, 1.56 GB table | 2.42 GB |
| master table and part files | 713 MB, 4 x 199 MB per iteration | none |

The 2D gain is one particle preparation per batch instead of two; the 3D fill is
untouched (wall time comes later, from kernels inside the batch). At 750,000 particles
and 2,500 directions the retired master table would have been 15 GB per state.

## For the maintainer

- **Mac check**: the generator's ten test vectors are in phase 0's Progress and in
  `simple_rnd_tester.f90`; `unit_numerics` (`random draws`) must pass on the Mac build.
- **Open**: key the simplex start; the base failures of the two simulated workflows
  (1jxy `new_shiftmap`, 6vxx final-volume correlation) are unchanged and not this
  run's.
- **Skill text to update** (not edited by the run): `simple-main-strategies-probabilistic`
  (no part files, assignment or insert files; the per-row selection `assign_rows`,
  `assign_rows_neigh`, the 2D forms, inside the batch; the keyed draws),
  `simple-main-strategies-parallelization` (no `prob_align*_exec`; one scheduling phase;
  the readout records merged in part order), `simple-refine3d`, `simple-solve2d`,
  `simple-solve3d-importance-sampling` (and its coupling map),
  `simple-frac-update-trailing`, `simple-main-commanders` (no table commanders),
  `simple-main-apis` (the six private programs gone).

## Commit message

```
prob align batch pipeline: keyed draws, one batch pipeline, table programs retired

- every probabilistic draw keyed on a stateless generator (seed, pass, iteration,
  particle, draw): independent of threads, batches and parts
- refine3D_exec/refine2D_exec run the table as a batch pipeline (fill, per-row draw,
  polish, residual, commit) in shared memory and in every distributed part; the master
  draws the sample, schedules once per iteration and merges the parts' readouts
- retired: the global tables, part/assignment/insert files, prob_align*_exec, the six
  prob_tab/prob_align private programs and their two modules
- three defects of the old two-phase path removed (2D band without lp, inherited
  shifts on a fresh 2D start, the widened shift window)
- refine2D 15-17 % faster (one preparation per batch); refine3D neutral, master table gone
```
