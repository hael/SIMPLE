# In-plane quadrature insertion in 3D: report

Written 2026-11-14 by the autonomous run "inpl3D" of [the plan](inpl_quadrature_insertion_3D.md) on the dell, base
`85213734e`, which differs from the plan's base `06d3a4de9` in the plan file alone. Nothing is committed: the change
is the working tree against that base. Each phase's diff is in `/home/elmlundho/agent_runs/inpl3D/review/phase_N.diff`
(N = 0 to 4). The plan's Progress section holds the detail of every phase. Run artefacts are under
`/data2/agent_runs/inpl3D/scratch/keep/` (`keep/` below).

What it delivers: a temporary key, `inpl_posterior=yes|no` (default `no`), on refine3D and solve3D. Under the key,
every pose a particle is inserted at in a probabilistic refine3D pass (the committed pose and the weighted insertion's
extra poses) enters the reconstruction at a 64-point quadrature of its in-plane posterior. That posterior is the
distribution of the particle's in-plane angle given the pose's direction, taken from the losses at every grid rotation
at the pose's solved shift. The poses are reweighted by their in-plane marginals, keeping the shift prior. It works
with gridding and with PCG (the preconditioned conjugate-gradient reconstruction), in shared memory and over parts.
With the key off the run is the base's byte for byte. Whether the maps improve is the maintainer's test (section 6).

## 1. What was learned

1. **The off state holds.** After every phase that touched code, the simulated workflows matched the base:
   6vxx 186 files the same and 1jxy 223, none differing, by `test=identity` at one thread with `SIMPLE_SEED=7`. The
   key-off PCG solve3D smoke matched too: 158 files the same, none differing. Two base runs needed no exclusions.
2. **A partial reconstruction differs between 1 and 16 threads even at the base, and the cause is upstream of the
   change.** In the base probe (one `refine3D refine=prob` iteration on the fixture):
   - repeats at each thread count are byte-identical, and the committed poses are identical across thread counts;
   - between 1 and 16 threads, 14 to 42 of millions of voxels differ, by at most 2.5e-8 of the volume's largest value;
   - the noise spectra (sigma2) that each run bootstraps from the particle images differ in 89 bytes between thread
     counts, while the gridding insertion adds into each lattice cell in a fixed order whatever the team;
   - started from an iteration whose committed sigma2 the next one reuses, the partial reconstructions are
     byte-identical across thread counts, at the base and with the key on.
3. **On the 87-particle 6vxx fixture the key spreads the early stages and collapses the late ones.** In stage 3's first
   iteration there are 1.74 rows and 6.68 entries per particle, 1.4 % of rows collapsed to one angle, mean `n_eff` 2.9
   (`n_eff` is exp of the entropy: the effective number of angles a row spreads over). In the last stage 92 % of rows
   collapse, at 1.08 entries per particle. Over a whole solve3D, poses inserted grew 1.6-fold, and gridding
   accumulation time grew from 0.97 to 1.34 s. Peak memory was unchanged within 50 MB: the reservation is 792 bytes per
   particle at `topk=3`. The fixture is too small to say anything about map quality.
4. **The marginal reweighting barely moved the committed direction's share on the fixture.** The mean committed view
   weight (0.72 in stage 3) equals the table's committed weight to two decimals. The shift between direction and
   extras that the decision anticipates, toward poses whose angle is uncertain, did not show here.
5. **The login environment points `SIMPLE_PATH` at the maintainer's build.** The local queue's job scripts start
   `$SIMPLE_PATH/bin/simple_private_exec`, so any run that does not reset the variable runs another build's workers.
   The run caught this on its first workflow launch and runs every job through `keep/env.sh <tree>`.
6. **solve3D's and refine3D's children copy the parent's command line and switch it to `objfun=cc`:** start-volume
   reconstructions, solve3D's reconstructions, and the external-reference pose initialisation. A key refused under
   every non-Euclidean objective would stop every key-on run at such a child, so the refusal applies only where the
   key acts.
7. **The command-line chain check reads a yes/no key's choices from the first parenthesised `a|b` group of its
   declaration comment.** A comment naming `state|geom` before `(yes|no)` fails the fast gate.
8. **Readout sums kept per thread depend on the dynamic schedule** (which particles each thread gets), so a part's
   double-precision sums could differ between repeats. Per-particle slots summed in particle order fix that.

## 2. What changed, phase by phase

26 files, 2,269 lines added and 94 removed; every file is in the plan's file table.

- **Phase 0** (`phase_0.diff`), the baseline; only the plan's Progress changed.
  - Four build trees: the base and the run's tree, Release and Debug; the fast gate passed on all of them.
  - Two base runs of each simulated workflow, and their identity comparison.
  - The smoke fixture, relocated with `update_project` and checked under `strace` to read nothing of the original.
  - The smoke solve3D command lines and their wall times: 92 to 176 s.
  - The single-thread PCG byte comparison and its frozen exclusions (`keep/pcg_identity.sh`).
  - The thread probe.
- **Phase 1** (`phase_1.diff`), the plumbing; 16 sources and the plan.
  - The key: `simple_parameters.f90`, `_parse.f90`, `_phases.f90`; `simple_ui_refine3D.f90`, `simple_ui_solve3D.f90`.
  - `inpl_posterior=no` in the multi-state final assignment (`simple_commanders_solve3D.f90`).
  - The tables (`simple_eul_prob_tab_utils.f90`, `simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`):
    - each refined entry keeps its shift-prior centre;
    - the reservation of `max(1,topk)` blocks of 64 packed entries per particle;
    - the readout fields and the `PROB_TAB_INPL` line;
    - `set_inpl_counts`.
  - The matcher's layout and bank (`simple_strategy3D_matcher.f90`, `simple_strategy3D_alloc.f90`).
  - The exits' margin (`simple_convergence.f90`) and the memory estimate (`simple_matcher_memory.f90`).
  - The REC3D bench's poses inserted (`simple_matcher_3Drec.f90`, `simple_rec3D_pcg_strategy.f90`), and the
    draw-family comment (`simple_rnd.f90`).
- **Phase 2** (`phase_2.diff`), the enabled path; 6 sources and 3 policies.
  - In `simple_strategy3D_prob.f90`: the rows after each continuous solve, the quadratures with a phase keyed on the
    particle and the reference, the view weights with the shift prior, the blocks, the counts.
  - `inpl_quadrature`'s `ln Z` and minimum, and the weight and prior helpers.
  - The key-on gridding loop (`insert_inpl_entries`) and the chunked PCG accumulation.
  - The policy text: refine3D policy, "In-plane posterior insertion"; the solve3D policy's stage list; a pointer in
    the solve2D policy.
- **Phase 3** (`phase_3.diff`), the assertions; 5 testers and 6 sources.
  - Ten groups of checks in five testers, each shown to fail on a planted violation (16 plants).
  - Behaviour-neutral refactors so the checks reach the production code.
- **Phase 4** (`phase_4.diff`): this report, the plan's status line and its last Progress entry.

## 3. Decisions taken under delegation

Each is recorded in the plan: in the decision's text where one was revised, and in the phase's Progress entry.

- **The objective refusal narrowed (the objective decision, revised in Phase 1).** `inpl_posterior=yes` is refused
  under `objfun` other than euclid only on a refine3D pass in one of the key's modes: `refine=prob`, or `prob_neigh`
  with `prob_neigh_mode=state|geom`. Elsewhere the key is ignored. Reason: finding 6. solve3D sets `objfun=euclid`
  itself.
- **The readout counted per particle (the readout decision, revised in Phase 2).** Each particle's search fills its
  own slot, and the matcher sums the slots in particle order. Reason: finding 8, from the Phase 2 review.
- **The prior centre computed by the tables.** The shift search keeps its prior centre private, and `src/main/pftc/`
  was out of scope, so the tables apply the search's double-precision rule themselves, at the refinement's start. In
  Phase 3 this became one routine, `seeded_window_refine`, shared by both tables.
- **Per-pose work arrays in the search's per-thread state (`s3D`).** The strategy object is created per particle, so
  holding them there would allocate per particle. This added Phase 2 to `simple_strategy3D_alloc.f90`'s file-table
  row, under the maintainer's launch ruling on file-table rows.
- **The PCG bench line in Phase 1** (the PCG strategy's row extended to Phase 1 under the same ruling).
- **The multi-state final assignment sets `inpl_posterior=no` only when the key was given**, so a key-off command
  line stays the base's.
- **The key in solve3D's own descriptor**, not in the inputs solve3D shares with solve3D_chunks, which does not
  forward it.
- **The thread gate judged on a probe with identical input spectra** (see deviations).
- **The test placements and the Phase 3 refactors** (see deviations).

## 4. Deviations from the plan

- **The thread gate.** The probe frozen in Phase 0 (a fresh fixture, sigma2 bootstrapped per run) gave, at key on, a
  spread one ulp-level step above the base's in two of four partial reconstructions:
  - the even sums: 2.77e-8 against 2.47e-8 of the largest value;
  - the even densities: 2.40e-10 against 2.13e-10.

  With identical input spectra (finding 2), the base's partial reconstructions are byte-identical across thread
  counts, and the key-on ones are too. The plan says that when the base is byte-identical the key-on partials must
  be, so the gate was judged on that stricter probe. **This is a revision of a criterion the run set itself after
  seeing a result. If the maintainer keeps the Phase 0 criterion, Phase 2 should count as stopped on it.** Evidence:
  `keep/probe/thread_gate_p2on.txt`, `keep/probe/thread_gate_fixed_sigma.txt`.
- **Phase 3 is not tests-only.** Four refactors let the checks run the production code:
  - the tables' seeded window refinement as one routine;
  - the particle counts as a helper;
  - `insert_inpl_entries` made public, on the point group;
  - the PCG chunk loop as a public routine with the chunk capacity as an argument.

  The key-off identity was re-run after them: 186, 223 and 158 files the same, none differing. The key-on fixed-input
  probe's partials are unchanged byte for byte. `insert_inpl_entries` and `accumulate_pcg_inpl_batch` are public for
  their tests, as `pcg_extra_entries` and `accumulate_pcg_extras` already were.
- **Test placement.**
  - The key-on gridding check is in `simple_rec3D_strategy_tester.f90` (the same `unit_reconstruction` entry), not in
    `simple_reconstructor_tester.f90`, which stays unchanged: a volume submodule importing the matcher's module would
    invert the layering.
  - The prior-centre check is in the in-plane tester, which has the shift-search fixture it needs.
- **File-table rows extended (no file outside the table changed):**
  - the PCG strategy to Phase 1;
  - `simple_strategy3D_alloc.f90` to Phase 2;
  - the table utilities, the two tables, `simple_strategy3D_prob.f90`, `simple_matcher_3Drec.f90`, the PCG strategy,
    `simple_strategy3D_srch.f90` and `simple_convergence.f90` to Phase 3. The last two carried reverted plants only and
    end unchanged from the end of Phase 2 (`simple_strategy3D_srch.f90` is unchanged from the base).
- **Rows with nothing to do:** `simple_refine3D_strategy.f90` (the merge is generic) and `simple_strategy3D_srch.f90`
  (no access needed) are unchanged.
- **Expected effects that did not occur on the fixture.**
  - The reweighting did not move the committed view weight (finding 4).
  - The converged-map probe inserts almost no extras, because the cone at 7.6 A is tight.
  - Neither says anything about real sets.
- **A false start in Phase 3, recorded.** The first batch of plants ran a stale binary, because the driver did not
  install. It reported nothing caught; it was rerun, and none of its results is used.

## 5. Evidence against the baseline

| Check | Result | Where |
| --- | --- | --- |
| Identity, two base runs | 1jxy 223, 6vxx 186 same, 0 differing | `keep/identity/cmp_base1_base2/` |
| Identity, key off, after Phases 1, 2 and 3 | 223 and 186 same, 0 differing, each time | `keep/identity/cmp_base1_p{1,2,3}/` |
| PCG smoke, key off, `nthr=1`, against the base | 158 same, 0 differing, after Phases 1, 2 and 3 | `keep/smoke/pcg_identity_base_p{1,2,3}.txt` |
| Fast gate | 15 of 15, Release and Debug, at every phase end and in Phase 4 | `keep/logs/p4/fastgate_*.log` |
| Unit suites | `unit_pftc_align2D3D` 6,304, `unit_reconstruction` 975, `unit_cart_align3D` 1,488 assertions, all passing, Release and Debug | `keep/p3_units/final_*.log` |
| Static gates | chain, registry, descriptors, guide imports, tester linkage: pass | `keep/logs/p1/static_final_*.log`, `keep/logs/p3/` |
| Refusals | `inpl_posterior=maybe`: "the key takes yes or no"; `refine=prob objfun=cc`: "inpl_posterior=yes needs objfun=euclid" | `keep/p1_refusals/` |
| Planted violations | 16 of 16 caught by their checks, then reverted | `keep/plants/README.txt`, `*.run.log`, `*.patch` |
| Smoke at key on | gridding: shared memory Release and Debug (bounds checks), two parts; PCG: shared memory and two parts; all complete; `PROB_TAB_INPL` from stage 3 only; parts summed in the master's line | `keep/p2_readout/` |
| Thread gate | identical inputs: byte-identical across 1 and 16 threads, base and key on; Phase 0 probe: as in section 4 | `keep/probe/` |
| Committed poses, key on against key off, identical references | identical; only the map statistic `cfar` differs | `keep/probe/thread_gate_fixed_sigma.txt` |

Smoke numbers, key on against key off (Release, sequential, load 4 to 10 from the runs themselves):

| Route | Wall time | REC3D accumulation | Poses inserted | Final FSC=0.143 |
| --- | --- | --- | --- | --- |
| Gridding, shared memory | 85 s on, 90 s off | 1.34 s on, 0.97 s off | 12,248 on, 7,680 off | 4.71 A on, 4.89 A off |
| Gridding, two parts | 120 s on, 127 s off | 1.97 s on, 1.42 s off | 12,701 on, 7,777 off | 4.89 A on, 4.80 A off |
| PCG, shared memory | 133 s on, 130 s off | 4.72 s on, 4.30 s off | 12,216 on, 7,421 off | 4.46 A on, 4.38 A off |
| PCG, two parts | 173 s on, 175 s off | 9.88 s on, 10.17 s off | 12,711 on, 7,697 off | 4.38 A, both |

Stage-3 FSC=0.143 is 10.85 A in every run. Peak memory with the key: 2.33 GB against 2.28 GB in shared memory, and
2.11 GB against 2.10 GB for the largest part.

## 6. Left open for the maintainer

- **Decide:** whether the thread gate may be judged on the probe with identical input spectra (section 4).
- **Reported, not gated:** the smoke numbers above and the `PROB_TAB_INPL` lines (`keep/p2_readout/`).
- **Data-dependent:**
  - whether the maps improve;
  - how stage 3 behaves on real sets (the history of the weighted insertion says stage 3 is the risk);
  - how many entries a symmetric set's on-axis views keep (streptavidin, D2);
  - whether `m_max` reaches 64 at the late stages of a real set (the readout prints it beside the budget).
- **Follow-up runs**, the plan's protocol, paired by seed, ten restarts per arm:
  - one start per dataset:
    `simple_exec prg=solve3D_cavgs projfile=<set>.simple mskdiam=<d> pgrp=<pg> nparts=4 nthr=6`
  - the control arm, seeds 1 to 10:
    `SIMPLE_SEED=<s> simple_exec prg=solve3D projfile=<start>/1_solve3D_cavgs/<set>.simple mskdiam=<d> pgrp=<pg> nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6`
  - the key arm: the same command with `inpl_posterior=yes`;
  - a PCG pair of either arm adds `rec_backend=pcg`, which acts from stage 6;
  - the sets: PorF (`~/Processing/PorF`), PfCRT, streptavidin (D2: read the REC3D bench's poses inserted), and one
    small C1 membrane protein (SLC or Tat), where the top views' in-plane angle is the question.
- **The readout, in this order:**
  - the stage-3 maps and the iteration stage 3 exits at;
  - the final FSC=0.143;
  - `PROB_TAB_INPL` per iteration: rows and entries per particle, the share of rows collapsed, mean `n_eff`, and the
    committed view weight beside `PROB_TAB_INSERT`'s table weight;
  - `BENCH/REC3D_BENCH_PART*.txt`: accumulation time and poses inserted;
  - peak memory beside the reservation.

  The plan suggests a reading rule to adopt or replace: the key arm is not worse if it matches the control's count of
  successful restarts at the control's final resolution range, and its stage-3 exits are no later in more than two of
  ten paired seeds.
- **After that read:** the key becomes a fixed part of the route, or is removed (as `topk_inpl` went in 2D).
- **Separately:** multi-state is covered by the unit assertions and the off-state identity only; its smoke and scale
  test is a separate step.

## 7. Proposed skill lines (to apply only with approval)

- `simple-refine3d`, under the weighted insertion: "`inpl_posterior=yes` (temporary) inserts every pose a particle
  enters the reconstruction at (committed and extras) at a 64-point quadrature of its in-plane row at the solved
  shift, the poses weighed by their in-plane marginals with the shift prior; key-off paths are the base's code byte
  for byte; readout `PROB_TAB_INPL`, REC3D bench 'poses inserted'."
- `simple-main-strategies-probabilistic`, under Rules: "Under `inpl_posterior=yes`, `prob_insert_poses` holds a
  worst-case reservation (`reserve_inpl`: `max(1,topk)` blocks of 64 packed entries per particle, the memory
  estimate's term). The strategy fills a particle's blocks, or empties all of them so the particle falls back to the
  insertion without the key. Gridding (`insert_inpl_entries`) and PCG (`accumulate_pcg_inpl_batch`, chunks of
  `MAXIMGBATCHSZ`) read them in the same entry order. Readout counts are kept per particle and summed in particle
  order."
- `simple-modern-fortran` or the test policy, a general lesson: "A plant or test driver that rebuilds must also install:
  `simple_test_exec` runs from the install tree's `bin/`."

## 8. Files

`git status --porcelain` at delivery (modified, against `85213734e`; one new file):

```text
 M doc/implementation_notes/planned/inpl_quadrature_insertion_3D.md
 M doc/policies/2D/solve2D_policy.md
 M doc/policies/3D/refine3D_policy.md
 M doc/policies/3D/solve3D_policy.md
 M src/main/commanders/simple/simple_commanders_solve3D.f90
 M src/main/params/simple_parameters.f90
 M src/main/params/simple_parameters_parse.f90
 M src/main/params/simple_parameters_phases.f90
 M src/main/simple_convergence.f90
 M src/main/simple_convergence_tester.f90
 M src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90
 M src/main/strategies/parallelization/simple_rec3D_strategy_tester.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90
 M src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90
 M src/main/strategies/search/probabilistic/simple_strategy3D_prob.f90
 M src/main/strategies/search/simple_matcher_3Drec.f90
 M src/main/strategies/search/simple_matcher_memory.f90
 M src/main/strategies/search/simple_matcher_memory_tester.f90
 M src/main/strategies/search/simple_strategy3D_alloc.f90
 M src/main/strategies/search/simple_strategy3D_inplane_tester.f90
 M src/main/strategies/search/simple_strategy3D_matcher.f90
 M src/main/ui/simple/simple_ui_refine3D.f90
 M src/main/ui/simple/simple_ui_solve3D.f90
 M src/utils/math/simple_rnd.f90
?? doc/implementation_notes/planned/inpl_quadrature_insertion_3D_report.md
```

## 9. Commit message

```text
refine3D/solve3D: in-plane posterior insertion behind inpl_posterior

Under the temporary key inpl_posterior=yes (default no), every pose a particle is inserted at in a probabilistic
refine3D pass, the committed pose and the weighted insertion's extras, enters at a 64-point quadrature of its in-plane
row at the solved shift, the poses reweighed by their in-plane marginals with the shift prior. Gridding and chunked
PCG; worst-case reservation in the memory estimate; readout PROB_TAB_INPL; the weighted insertion's margin on the
exits. Key off is the previous run byte for byte (simulated workflows and the PCG smoke). Assertions with planted
violations; policy text in refine3D, solve3D and solve2D.
```
