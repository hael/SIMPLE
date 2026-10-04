# Review report: trailing reconstruction without the finished-halfmap blend

Date: 2026-10-04. Unattended run on the Oracle Linux test machine (Dell Pro
Micro, 24 cores), three phases. Nothing is committed: the change is the
working tree against base commit `d32aa9d16`. The phase diffs are in
`/home/elmlundho/agent_runs/trailing_halfmap/review/phase_0.diff`,
`phase_1.diff` and `phase_2.diff`. The full record, including every command
line and per-iteration table, is section 12 of the plan beside this file
(`trailing_reconstruction_without_halfmap_blend.md`).

## Terms

- **Fractional update.** With `update_frac` below 1, each 3D refinement
  iteration aligns and reconstructs only a sample of the particles. `f` is the
  sample's share of the particles updated so far.
- **Trailing reconstruction** (`trail_rec=yes`). Blends the current sample
  with what came before. "What came before" is kept on disk as an
  **accumulator chain**: unregularized Fourier sums and sampling densities,
  one per state and half set.
- **Finished half maps.** The even and odd half maps after restoration
  (density correction, regularization, deapodization).
- **FSC.** The Fourier shell correlation between the two half maps. Its 0.143
  crossing, in ångström (Å, smaller is better), is the resolution used here.
- **Gridding and PCG.** The two 3D reconstruction backends. Gridding is the
  default; PCG solves the reconstruction by the preconditioned
  conjugate-gradient method.
- **solve3D and refine3D.** solve3D determines a map de novo in stages.
  refine3D refines an existing one. volassemble is the gridding step that
  assembles the parts' partial reconstructions into maps.

## What changed and why

The maintainer ruled that "no trailing should ever happen on finished
halfmaps". When a blend was due but no chain existed yet, both backends:

- seeded the chain correctly;
- for that one iteration, took the FSC from the previous iteration's half
  maps;
- blended the restored half maps with those previous maps.

Now that iteration seeds the chain from the current sample at full mass
(sums times `1/f`) and ships the current sample's own map. It reads no
previous map and runs no extra reconstruction. From the next iteration on,
every blend uses the chain, as before.

- **Phase 0 (`phase_0.diff`, the plan only).** Debug build and fast gate,
  then a baseline of the old code on beta-galactosidase (5,513 particles):
  - case A: solve3D with gridding, three runs;
  - case B: solve3D with PCG, three runs;
  - case C: refine3D with `trail_rec=yes` and `update_frac=0.2`;
  - case D: a two-state attempt;
  - the two simulated workflows.
- **Phase 1 (`phase_1.diff`).**
  - In `simple_commanders_rec_distr.f90` (volassemble) and
    `simple_rec3D_pcg_strategy.f90` (PCG master), the previous-map FSC, the
    finished-map blend and the PCG bootstrap output are removed. The callers
    in `simple_refine3D_strategy.f90` and `simple_rec3D_strategy.f90` lose
    the matching argument.
  - The support-provenance kind `mixed` goes from
    `simple_halfmap_diagnostics.f90`.
  - The new rule for a state with almost no sample is described under
    decisions below.
  - The tester `simple_accum_blend_tester.f90` gains a chain-start case:
    60 checks, against 40 before.
  - The trailing skill, its contract reference and the fractional-update and
    PCG policies are updated.
- **Phase 2 (`phase_2.diff`, documents only).**
  - Cases A to D repeated with the Phase 0 command lines.
  - Item C6 removed from the release 4 inventory, and the "What remains"
    paragraph of the release 4 report updated.
  - The plan moved to `completed/`, and this report written.

## Evidence against the baseline

Final map, FSC 0.143 (Å). Acceptance (plan 11.3) as read in section 12: a
Phase 2 run may be worse than the Phase 0 mean by at most the larger of the
Phase 0 spread and 3 % of that mean.

| Case | Phase 0 (old code) | Bound | Phase 2 (new code) |
| --- | --- | --- | --- |
| A, gridding | 3.93, 3.89, 3.89 | 4.02 | 3.80, 3.93, 3.98 |
| B, PCG | 3.93, 5.18, 4.29 | 5.72 | 3.93, 3.89, 3.93 |
| C, refine3D | 4.13; 4.24, 4.24 (old code, rerun in Phase 2) | 4.33 | 4.30, 4.24, 4.24 |

- **Final maps and resources.** All cases pass. Wall time and peak memory
  are unchanged: case A about 28 minutes and 3.1 GiB; case B about 40 minutes
  and 6 to 6.6 GiB.
- **Log lines.** Every start-without-chain iteration logs "SEEDED FULL-MASS
  TRAILING CHAIN ...; THIS ITERATION USES THE CURRENT SAMPLE ONLY" and no
  blend; all later iterations blend from the chain. No run of the new code
  prints the old legacy line.
- **The switch-on iteration.** Its printed FSC 0.143 drops, because it now
  describes the 1,000-particle map actually shipped rather than the previous
  maps:
  - solve3D: 5.53 Å became 7.6 to 8.2 Å;
  - case C: 3.89 Å became 4.53 to 4.60 Å.

  One or two iterations later the values match the baseline, and the
  orientation and shift changes match throughout.
- **Gridding check (Phase 1).** On identical partial reconstructions, without
  `vol1`, volassemble with `trail_rec=yes` and no chain gives the same maps
  as `trail_rec=no` (the current sample alone), to a relative 4e-7.
- **No-sample check (Phase 1).** It passes on both backends.

## Decisions taken under delegation

Each is recorded in the plan where noted.

1. **Run configuration (section 12, Phase 0).**
   - Long runs use an out-of-tree Release build: Debug is too slow for
     solve3D.
   - `nsample=1000` in cases A and B, because the default forces full
     sampling on 5,513 particles and trailing never switches on.
   - `startit=2` in case C, because `refine=prob` with `startit=1` clears the
     update history and no blend would be due.
   - A third PCG run, because the first two diverged.
2. **The case D outcome is recorded rather than forced (section 12,
   Phase 0).** Case D could not be produced: refine3D always moves some
   sampled particles into the small state, so it is never left unsampled.
3. **No-sample rule (section 6.1 and the Phase 1 notes).** This covers a
   state with at least one sampled particle but `f` below 0.001:
   - with a chain, the chain carries unchanged and the iteration is
     restored from it;
   - without a chain, the previous volume is carried forward, provided the
     run directory holds it with its half maps and FSC;
   - with neither, the run stops with an error that says to raise
     `update_frac`.

   Seeding from the sample in that last case was tried. It ships a map of
   about one particle; postprocess crashed on it and PCG refused a half set
   with no particles. Reading `vol<state>` would bring back the removed
   dependence. Both backends apply the same rule.
4. **File-table additions in Phase 1, recorded in the table.**
   `simple_halfmap_diagnostics.f90`, because `mixed` had no other writer.
   One stale comment in `simple_refine3D_strategy.f90`.
5. **Gridding check as an integration script (Phase 1 notes).** No
   volassemble unit-test harness exists, so the check of plan section 8 runs
   volassemble by hand on real partial reconstructions
   (`scratch/gridding_check.sh`) instead of adding a test file.
6. **PCG log line (section 6.2 and the Phase 1 notes).** The PCG master now
   logs the same seed line as gridding, as plan section 9 expects; before,
   it logged none.
7. **Acceptance reading and the case C spread (section 12, Phase 2).** The
   Phase 0 spread of case C (one run) was completed with two runs of the old
   code in Phase 2. They came from an out-of-tree build of the base commit
   via `git archive`, deleted afterwards. No tolerance was changed.

## Deviations from the plan

- **Section 8.** The gridding check is a script, not a test in the
  repository (decision 5).
- **Section 6.1.** The no-sample rule adds an error for the one case where
  there is nothing to restore or carry (decision 3).
- **Section 11.1.** Case D, which the plan marks optional, was not obtained
  (decision 2).

## Open for the maintainer

- **Section 10 questions.** Marking the seeding iteration in the convergence
  output would explain the printed FSC drop, which is cosmetic. Nothing
  measured argues for seeding solve3D's chain earlier: the final maps are
  unchanged.
- **Unchanged and pre-existing behaviour.**
  - A carried-forward state of a multi-state gridding run gets resolution 0
    in the project, as a dropped state already does.
  - solve3D_addon with PCG and a near-empty cohort sample still blends by
    the population rule, where gridding carries the chain.
  - The NU (nonuniform) filter still accepts an evidence source,
    `previous_shipped`, that nothing produces.
- **No CTest coverage of trailing.** Neither simulated high-level workflow
  reaches trailing: about 80 particles force full sampling. No CTest entry
  exercises trailing end to end.
- **Flaky simulated workflow.** `simulated_workflow_6vxx` failed once on the
  unchanged base (a failed de novo solution) and passed on the rerun.
- **Code indexes.** `doc/code_overview/fortran-indexes` was not regenerated,
  as instructed.
