# Single-state speed, days 1 and 2 and the 2D table: report of the run `prob_speed`

Report, 2026-10-15, of the unattended run `prob_speed` on the test machine
`dell` (Intel Core Ultra 9 285, 24 cores, Oracle Linux 8.10). It carried out
days 1 and 2 of `single_state_speed_week.md` (beside this report) and the
same changes in the 2D probability table, in five phases. Nothing is
committed: the change is the working tree against `1bbbdcc24`. The plan's
section "Progress" holds the detail (command lines, logs, per-run numbers);
paths below are under `/home/elmlundho/agent_runs/prob_speed`.

## What changed, phase by phase

The diff of each phase is `review/phase_N.diff`; `review/full_after_N.diff`
is the whole tree after phase N.

- **Phase 0, baseline** (`review/phase_0.diff`). `simulated_workflow_1jxy`
  runs `solve3D` a second time from the same `solve2D` project with
  `nparts=2` (two local worker processes), checked against the same truth.
  That step exposed a defect: the final reconstruction of a distributed
  one-thread `solve3D` left `nthr` off its workers' command line, and the
  workers stopped. Fixed in `simple_final_rec.f90`. The baseline cases were
  run on the unchanged code.
- **Phase 1, day 1** (`review/phase_1.diff`).
  - The `prob_neigh` pre-step driver writes its own benchmark report
    (`PROB_ALIGN_NEIGH_BENCH.txt`).
  - The master's job loops (`simple_qsys_ctrl.f90`) refresh the done flags
    before testing them and poll every 0.1 s on the local queue (1 s stays
    for batch queues and the persistent worker).
  - The redundant project-segment write before the distributed pre-step is
    gone; the guard flag of the read before it was always true and is
    deleted.
  - Every 3D table fill with a shift search runs with dynamic OpenMP
    scheduling.
  - The summary script totals the tables' thread-seconds.
- **Phase 2, day 2** (`review/phase_2.diff`). The developer key `prob_fast`
  is retired (parameters, parser, `refine3D` and `solve3D`-family
  interfaces; a command line passing it is refused). Its particle-frame paths
  are the only paths in the `prob` table, and the `prob_neigh` table gets the
  same path: the shift seed is converted once to the particle frame, the
  particle is shifted once and every reference, in the coarse and the fine
  pass, is scored against it with the bank kernel (which includes the angular
  Nyquist bin); refinement starts from the seed converted at each reference's
  rotation; the assignment records the seed as it is.
- **Phase 3, the 2D table** (`review/phase_3.diff`). The 2D table part
  reports the same fill profile as the 3D ones; both 2D fills (the plain one
  of `refine=prob`, the `prob_snhc` one) use the particle-frame seed and
  dynamic scheduling. After the stop and the rulings: a root seed per
  simulated suite, and the fix of the map-symmetrization crash (below).
- **Phase 4, close** (`review/phase_4.diff`, written by the driver): the
  living documents (the per-suite seeds in `test_environment_policy.md`,
  the symmetry search's band in `solve3D_policy.md`; the policies of
  `refine3D`, `solve3D_chunks` and `solve2D` were brought up to date in the
  phases that changed their subjects), the final repeats, this report.

## Evidence against the baseline

Beta-galactosidase, 5,513 particles, box 256; `solve3D` from the common
`solve3D_cavgs` result with `nsample=2500 cavg_ini_ext=yes`; distributed
runs `nparts=4 nthr=6`, shared memory `nthr=24`. Release builds, one case at
a time on an otherwise idle machine. Phase 0 is the mean of its runs (two,
three for A and A-shm where a third before-run was made in Phase 1); Phase 4
is one final run of the same command line (3D without `prob_fast`).

| case | Phase 0 | final | change | FSC 0.143, Phase 0 / final |
| --- | --- | --- | --- | --- |
| A, `solve3D` distributed | 659.6 s | 478.7 s | -27 % | 4.030 / 4.030 Å |
| A-shm, `solve3D` shared memory | 436.6 s | 386.0 s | -12 % | 4.185-4.239 / 4.080 Å |
| E, `solve2D` distributed | 149.0 s | 80.0 s | -46 % | score 0.387 / 0.387 |
| E-shm, `solve2D` shared memory | 98.8 s | 95.1 s | -4 % | score 0.387 / 0.387 |

Where the time went, from the phase-by-phase comparisons (each at least
twice each way, interleaved):

- Job polling: A -115 s, E -61 s (the master's time around the slowest part
  in A fell from about 150 s to 28 s).
- Dynamic scheduling of the 3D fills: A -23 s; load balance of the table
  parts 0.75 to 0.96 (`prob`) and 0.78 to 0.94 (`prob_neigh`).
- Particle-frame shift in `prob_neigh`: A -41 s, A-shm -35 s; its score loop
  2.8 to 2.9 times fewer thread-seconds (926-948 to 319-340 in A; 299 in the
  final run).
- 2D table (particle-frame shift and dynamic scheduling): E -4.6 s, E-shm
  -6.7 s; 2D score loop 2.7 times fewer thread-seconds, load balance 0.98.
- A-shm's total varies with the number of iterations a run needs to converge
  (33 to 38), which blurs its per-item gains; its table time fell
  consistently. E-shm's single final run (95.1 s) lies at the top of the
  Phase 3 runs of the same build (87.7-91.4 s); the Phase 3 comparison, three
  runs each way, measured 7 % less than the build before it.

Acceptance of the final runs (plan, "Execution"): A's FSC 0.143 resolution
equals Phase 0; A-shm's 4.080 Å lies within the Phase 0 spread (4.030-4.239
Å over three runs) of its mean; E and E-shm keep the mean score 0.387 and
class resolutions of 10.03 and 9.84 Å (limits 10.19 and 10.09 Å).

Simulated tests on the final build all pass: `simulated_workflow_1jxy` (both
steps; whole-volume correlation 0.8863 and 0.8837), `simulated_workflow_6vxx`
(0.8763), `solve3D_addon` and `cont_refine3D_1jxy`. The identical changes of
Phase 1 left the one-thread suites byte-identical, the distributed step
included. The Debug build passes the fast gate (15 suites).

## Decisions taken under delegation

All are recorded in the plan's progress section of the phase named.

- The `nthr` fix of the final reconstruction (Phase 0), added to the file
  table: the distributed step could not pass without it.
- The byte comparison ignores two files' run-specific content: the run
  identifier in `solve3D_manifest.txt` and the layout digest in the header of
  `sigma2_state.bin`, which hash absolute paths (Phase 0).
- The polling interval: 0.1 s, local queue only (Phase 1).
- The project write: removing the write made the read guard visibly dead;
  deleted under the dead-code rule (Phase 1).
- The `prob_align_neigh` report as a family of its own,
  `PROB_ALIGN_NEIGH` (Phase 1).
- The A-shm acceptance band from three Phase 0 runs (a third before-run in
  Phase 1 gave 4.030 Å, inside the band once counted) (Phase 1).
- `record_seed_shift` moved to the shared table routines and the rotating
  `materialize_seed_shift` deleted once no table used it (Phase 3); the
  multi-state `prob_neigh` `shc` mode, which also seeds, uses the
  particle-frame seed too (Phase 2).
- The map-symmetrization fix's form: the search's high-pass is
  `max(hp, low-pass)`; a first form (an octave below the low-pass) put the
  origin into the band and was replaced (Phase 3).

## Deviations from the plan

- Files added to the file table: `simple_final_rec.f90` (Phase 0),
  `simple_eul_prob_tab_tester.f90` (Phases 2 and 3), the import changes in
  the 3D tables (Phase 3), `simple_solve3D_utils.f90` and
  `simple_volpft_corrcalc.f90` (Phase 3, ruling), `solve3D_policy.md`
  (Phase 4), the workflow tester extended to Phase 3 (ruling).
- Phase 3 stopped once: `simulated_workflow_6vxx` failed its truth check
  after the 2D change. The maintainer ruled the change stays and gave the
  suite its own seed (20261003); 1jxy keeps 20260923 and stayed
  byte-identical.
- The map-symmetrization crash was fixed in this run (ruling). Cause: on
  data so poor that the stage ladder plans the symmetry-search stage at
  124.8 Å, that low-pass met the default high-pass of half the box (124.15 Å)
  and left no Fourier shell; every grid score was undefined, the best index
  stayed 0, and reading orientation 0 stopped the run. Fix: the band always
  holds the low-pass shell (`simple_solve3D_utils.f90`), and the polar volume
  correlator refuses a band without a shell (`simple_volpft_corrcalc.f90`).
  The two crashing seeds now run through.

## Left open for the maintainer

- PfCRT (the plan's baseline dataset) was not run; the plan leaves that run
  to the maintainer.
- `solve3D_addon`'s cohort pose excess has come closer to its floor (1.0
  degrees in Phase 1, 2.6 in Phase 2, 4.5-4.6 in the Phase 3 runs and 3.7 in
  the final one, floor 5); it passes, but the margin is thin. The test runs
  on eight threads, so its draws vary from run to run.
- `simulated_workflow_6vxx` remains sensitive to the draw at 80 particles
  (two of four seeds missed the map in the diagnosis); more particles would
  make it a sturdier gate.
- The summary script still warns on every run of a key reported twice
  (VOLASSEMBLE at iteration 1) and a missing strategy record at the closing
  sigma iteration; harmless, not changed.
- Findings for day 3: the `prob` table now dominates the pre-step
  (score loop about 1,950 thread-seconds in A against 300 for `prob_neigh`);
  its shift refinement and the reference memo traffic are the next targets,
  as the plan has them.
