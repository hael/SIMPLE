# Report: retire `prob_state`; flex PCA is the only particle state initializer

Run `prob_state_flex`, 2026-10-06 to 2026-10-07, on dell (24 cores, GCC 15.2).
The plan, with its rulings and a detailed Progress section, is
[prob_state_retirement_flex_split.md](prob_state_retirement_flex_split.md)
beside this report. Nothing is committed: the change is the working tree
against commit `6914834d3`. The per-phase diffs are
`/home/elmlundho/agent_runs/prob_state_flex/review/phase_N.diff`; logs and
scripts are under `/home/elmlundho/agent_runs/prob_state_flex/scratch`.

## Outcome in one paragraph

`refine=prob_state` (keep each particle's projection direction, score it
against every state, draw the state) no longer exists. `refine3D_states`
splits a single-state project into states only with projection-aware flex
principal component analysis (flex PCA, `flex_pca`), which needs
`nstates >= 3`; its `flex` key, its `prob_state` initialization phase and
`pose_policy=fixed` are gone, and `pose_policy` is `local|global`. Following
your ruling of 2026-10-06, docked multi-state mode (one consensus model up to
a split stage, then states) is removed from `solve3D` and `solve3D_cavgs`
instead of being rebuilt around flex: `nstates` alone now selects single-state
or independent multi-state refinement. On beta-galactosidase (5,513
particles), the single-state `solve3D` reaches 3.80 Å (baseline 3.84 Å) and
`refine3D_states` 4.03 Å for its best state (baseline 4.03 Å); both smoke
tests of the multi-state mode complete with every state populated.

## What changed, phase by phase

- **Phase 0, baseline** (`review/phase_0.diff`, plan only). Release build of
  the unchanged base; six beta-galactosidase runs: the common class-average
  start S (`solve3D_cavgs`), the single-state `solve3D` from S, docked
  `solve3D nstates=3` twice (A1, A2), `refine3D_states` (B) and docked
  `solve3D_cavgs nstates=3` (C). The A logs showed that docked particle
  `solve3D` never ran `prob_state` (it split a cohort at random and handed
  over to `refine3D_states`), which led to your ruling. The ruling, the exact
  command lines and the acceptance limits are in the plan.
- **Phase 1, implementation** (`review/phase_1.diff`, 47 files). Removed
  `prob_state` from the stage controller, the refine3D strategy, matcher,
  allocator and probabilistic search, the probability-table commanders and
  the state-only table of `simple_eul_prob_tab.f90`, the parameters and the
  UI. Reworked `refine3D_states` (flex-only initialization, no init phase,
  `local|global`, greedy final pass). Removed docked mode: the split
  checkpoint module, cohort pass, handoff, `split_stage`,
  `sticky_class_sampling`, the `sampled_only` sampler path, the docked
  manifest fields; the controller and commanders decide by `nstates`; both
  workflows refuse `multivol_mode`. Updated policies, algorithm chapters, the
  how-to, the two skills, the memory benchmark script and two unit tests.
- **Phase 2, comparison runs and this report** (`review/phase_2.diff`): four
  runs on the new build, the plan moved to `completed/`, this report.

## Decisions taken under delegation

All are recorded in the plan (Rulings section or Progress).

1. Phase 0: case C, already running when the ruling arrived, was finished as
   context; no new baselines were run (ruling item 4).
2. Phase 0: the ruling's smoke test `solve3D nstates=3 cavg_ini_ext=yes`
   "from S" cannot run, because S has one state and that route requires the
   particles to carry the requested states (verified on the base build). It
   was run from the result of the `solve3D_cavgs nstates=3` smoke test.
3. Phase 1: `solve3D` and `solve3D_cavgs` refuse `multivol_mode` with an
   explicit message, because the parser accepts every known key for every
   program and a silently ignored mode would mislead.
4. Phase 1: the solve3D manifest keeps refusing unknown fields (its unit test
   makes that the release-4 rule), so manifests written before this change no
   longer read; only `solve3D_addon` replays of older runs are affected.
   `.simple` project files are untouched.
5. Phase 1: the `startit` offset of `refine3D_states`, used only by the
   docked handoff, was dropped; the multi-state minimum-population constant
   was renamed `PROB_MULTISTATE_MIN_POP`; a helper whose name contained
   `state_assign` was renamed so the plan's text check holds literally.
6. Phase 1: the generated Fortran indexes under `doc/code_overview/` were not
   regenerated.

## Deviations from the plan

- The ruling replaced the plan's flex split inside `solve3D` and the balanced
  class-average partition of `solve3D_cavgs`; neither was built, and the
  plan's two planned tests (balanced partition, docked `nstates=2` stop) were
  replaced by command-line checks: `refine3D_states nstates=2` on
  single-state input stops with "state=0/1 input requires nstates >= 3 for
  flex_pca initialization"; `pose_policy=fixed`, `flex=no`, `split_stage` and
  `solve3D multivol_mode` are refused.
- Phase 2 validation followed ruling item 4 instead of the plan's A/B/C
  repeat. The plan's title and "What changes" sections still describe the
  original flex-split design; they are marked as replaced where the ruling
  applies.
- Thirteen files beyond the plan's original table were edited, each added to the file table
  with its reason before editing.

## Evidence against the baseline

| Check | Baseline (Phase 0) | After (Phase 2) |
|---|---|---|
| single-state `solve3D`, final FSC 0.143 (Fourier shell correlation) | 3.84 Å | 3.80 Å (limit 3.96 Å) |
| `refine3D_states` B, best-state FSC 0.143 | 4.03 Å, 2 states (flex dropped one) | 4.03 Å, 2 states (limit 4.15 Å) |
| `solve3D_cavgs nstates=3` | (docked C: 9/19/22 class averages) | 30/5/15 class averages, all populated |
| `solve3D nstates=3 cavg_ini_ext=yes` | not runnable from S | 4,189/471/853 particles, all populated |
| Debug and Release builds | no compiler warnings, linker "executable stack" warnings only | same |
| fast test gate | 14 of 14 | 14 of 14 (Debug and Release) |

Also passing: `lib_heterogeneity`, `scripts/check_test_registry.py`, the
plan's text checks (no `prob_state`, `input_oris_fixed`,
`fill_tab_state_only` or `state_assign` in `src/`). An independent review of
the Phase 1 diff found no regression; it noted that `solve3D cavg_ini=yes
nstates>1`, which used to fail, now works.

## Open for the maintainer

1. The high-level test `solve3D_addon` fails at its first `solve3D` with
   "balance=cavg needs selected class averages" (its simulated particles have
   no class averages). It fails the same way on the unchanged base, so it
   predates this run.
2. `refine3D refine=prob_state` (like any unknown `refine` value, in the base
   too) is refused by the distributed workers while the master keeps waiting.
   Refusing it up front would need a list of refine modes kept beside the
   matcher.
3. In `refine3D_states`, the sampling target counts particles that flex
   dropped to state 0 (B: 5,513 instead of about 4,200), so the realized
   update fraction is higher than reported; present in the baseline too.
4. `multivol_mode` is now read only by the convergence rule of
   `refine3D_states`; it can be removed together with a replacement for that
   test (kept in this run, as the ruling asks).
5. Older solve3D manifests no longer read (decision 4 above).
6. Inputs that reproduce the Phase 2 runs remain in
   `/home/elmlundho/agent_runs/prob_state_flex/scratch/keep/` (S and the
   single-state `solve3D` result, 4.5 GB); delete them when no longer needed.
