# Multi-state ab initio: find what broke, restore what worked, until it separates

Plan of 2026-11-06 for the autonomous run `ms_abinitio` on dell2. Written after two
runs on Afan's three-state benchmark (`ms_posterior`, `ms_posterior2`; the second's
report is `doc/implementation_notes/completed/multistate_posterior_state_draw_report.md`)
and Hans's account of what worked: multi-state ab initio from a random partition of
the particles into states (the 2025 paper, `het_mode=independent`) resolved many
states on real data, and straight stochastic hill climbing with sampled in-plane
angles worked surprisingly often; a single-state solution split afterwards never
worked on real data.

Unlike the last run, this run does not stop at a verdict. Its goal is a de novo
multi-state `solve3D` that separates the benchmark's states, and it keeps working
the candidate list below until that holds, the list is exhausted, or a stop
condition of section 3 is met.

## 1. What is known

**The benchmark** (`/data2/multistate_bench` on dell2, read-only): 80,000 simulated
particles, 128 px at 3.0 A, SNR 0.01, three ribosome assemblies of 24,000 (state 1,
smallest), 48,000 (state 2, medium) and 8,000 (state 3, largest). True labels, poses
and maps are given. cryoSPARC's ab initio recovers the states (the README).

**What the earlier runs measured** (current code unless named):

- With the true poses and the truth maps as references, every particle is assigned its
  true state at 10 A. The state signal is there.
- From random balanced labels on the true poses, on the 20,000-particle subset S
  (6,000 / 12,000 / 2,000), two runs ended at accuracy 0.69 and 0.65: the smallest
  assembly separates, the medium is split over two states, the largest is merged into
  one of them. The joint state draw and the old argmin behave alike here. (Corrected in
  Phase 1: run 2 compared the two rules on the small `labels` route only, never on S. On S,
  label recovery is a stochastic fork for every build tested, see Phase 1's Progress.)
- Three de novo runs on S ended at 0.59, 0.42 and 0.58; the current code's de novo run
  of the full set before the October posterior change ended at 0.45 (not paired).
- In every de novo run the two hill-climbing stages reached 0.42-0.45 (chance is
  0.34), and the first probabilistic iteration re-drew about half the labels and
  erased it (stage 3 at 0.38-0.40). The state posterior at the search's factor 0.3
  spreads over 2.2-2.7 states at 20-23 A, and a quarter to half of each state's mass
  sits on entries the fill did not shift-refine.
- On S the 2D classes reach 22.6 A, so the five multi-state stages ran at 24, 22.6,
  22.6, 20.2 and then 6 A.
- The small simulated tests pass: `multistate_workflows suite=labels` resolves every
  state in 13 of 20 seeds (gated at 6), `suite=shc` in 15 of 20 (gated at 6). They do not
  reproduce the benchmark's failure.

**How the code changed since the method that worked** (the 2025 paper, Acta Cryst. D
81, ic5124; code tag `v3.1-beta`, 2025-05-27):

| What | The paper's code | Now | When |
| --- | --- | --- | --- |
| Early search | stages 1-4 `shc_smpl`: hill climbing over directions with the in-plane angle sampled within an angular neighbourhood (10 degrees to start, adapted) | stages 1-2 `shc` with the best in-plane angle; `shc_smpl` deleted | default changed about 2026-04; deleted 2026-10-06 (dffa21364) |
| Probabilistic draw | sampled within the angular neighbourhood of the best candidates (`prob_athres`), with a reference-balanced frontier (every reference draws its closest free particle) | one draw from the whole row's tempered posterior (`prob_beta` 0.3), no neighbourhood, no frontier | 2026-10-24 (811934ca5) |
| Stages | 8 for every state count | 5 for `nstates>1`, `lpstop` 6 A | 2026-06-08 (32a324ab9; corrected in Phase 0); the stages' content changed 2026-10-16/18 |
| State label | the best state | the joint draw over the row | 2026-11-05 |

Cyril reported that multi-state stopped working on 2026-10-20; the build of 2026-10-15
(1bbbdcc24) still separated his two-state set in 2 of 4 runs.

## 2. Rulings

1. Multi-state ab initio from a random partition of the particles into states. No
   single-state solution split afterwards, no flex initialisation, no docked mode
   (Hans, 2026-11-06).
2. Hard state labels in the reconstructions (2026-11-03). If the evidence says only soft
   state membership can separate the benchmark, the run stops and writes that up.
3. No new parameter of any kind (every field of `simple_parameters` becomes an accepted
   key of every program). A search mode that existed before may come back as a value of
   the existing `refine` key; the solve3D controller chooses the modes per stage and
   state count, as it does now. Experimental variants are patches applied to `git
   archive` exports in scratch, never keys.
4. `nstates=1` stays bit for bit identical at one thread (single-state solve3D was just
   optimised; protect it). A change that cannot be confined to `nstates>1` is a stop.
5. No floor, tolerance or timeout of an existing test is relaxed. New floors come only
   from a profile of the finished code by the floor rule.
6. The benchmark drives the development this time. Decisions are taken on it, on the
   fixtures of section 4, in the order section 5 sets.
7. State 0 is exclusion; a true state in state 0 is a failure. A de novo verdict rests on
   several runs, never one.
8. Stream code and the GUI are not edited; the 2D tables are not changed.
9. No de novo runs with the old builds (Hans, 2026-11-07, on the Phase 0 diagnostics so far). Phase 0
   ends with A25's label recovery on S and HEAD's label recovery on S at the full set's update
   fraction (0.375); the de novo runs of A25 and A1015 in Phase 0's step 3 are dropped, and the
   diagnostics are read by Phase 0's table without them.
10. No `freq_march_aggressive` (Hans, 2026-11-08; the key is to be removed): no test or candidate uses
    it. The ladder hypothesis (C3) is tested with the explicit-limits path that exists, `lpstart=<A>
    lpstop=6 force_lp_range=yes` on the `solve3D` command line (a linear ramp from lpstart to
    lpstop over the stages, taken as given), with an lpstart that does not sit at the 2D classes'
    22-24 A, combined with C1 and without it where the evidence needs it. If the forced ramp is what
    separates the states, the candidate kept is a controller default for `nstates>1` that plans a
    linear ramp instead of the class-FRC ladder (no new key), `nstates=1` unchanged.
11. The balanced subset B judges Phase 2 (2026-11-08, the maintainer's delegate). With the true poses
    S separates in label recovery only by a race over the largest assembly's 2,000 particles (HEAD 2
    of 6, the paper's code 2 of 2), and every de novo candidate on S ended at 0.41-0.47 after stage 4,
    whatever it changed. So candidates are judged by de novo runs on B (8,000 of each assembly; label
    recovery on true poses separates it reliably). B's 2D solution is made first (`solve2D ncls=30
    mskdiam=300 nthr=16` on a copy of PZB, as P2S was made). Before any candidate, one de novo run of
    HEAD unchanged on B: if it separates, HEAD de novo on F follows before anything is built. The
    stopping rule (a de novo run below 0.60 after stage 4 is stopped) and the candidate order stay. S
    stays in the report as a stress case with its counts.
12. The target is SNR 0.05 (Hans, 2026-11-08, after Cyril reported success on real data with the code
    this run started from): "SNR=0.01 is kind of intense. I'd be happy with code that can cope with
    SNR=0.05." Afan's SNR 0.01 benchmark stops driving the development; its findings (the rare-state
    race on S, no regression, the 2D-pinned ladder, HEAD de novo on B at 0.699, C1 and C3 below HEAD)
    are kept as a stress case. A benchmark F05 is simulated at SNR 0.05 by Afan's recipe (80,000
    particles, the three truth maps, populations 0.3 / 0.6 / 0.1), its truth checked as in Phase 0,
    with its project and 2D solution under `scratch/keep/`. HEAD unchanged de novo on F05 twice, then
    true-pose label recovery once. If HEAD separates both runs, no code change is needed (Phase 4
    next; Phase 3's gates apply only to a kept change); one of two -> a third run decides (2 of 3); if
    HEAD fails, the candidates (C4, C2, then C1 and C3 if the evidence favours them) are judged on F05
    with HEAD's runs as the bar. C2 is not run on the SNR 0.01 data.
13. C5 first on F05 (2026-11-09, the maintainer's delegate). Both HEAD runs on F05 fail the same way:
    the greedy class sampling of stages 1-3 searches a fixed sample, the particles outside it (most of
    the smallest assembly) keep their random start labels until stage 4, and the final states' strays
    are those late arrivals, not lower-quality particles. C5, stochastic sampling from stage 1 for
    `nstates>1`, therefore runs first, as "what the evidence names". C4 and C2 follow only if C5 does
    not reach the success criterion, under the same stopping and second-run rules. If C5's first
    run separates, its second and third runs come before anything else.

## 3. Stop conditions

The run continues until the success criterion of section 4 holds. It stops only when:

- the candidate list of section 5 (with what Phase 1 adds) is exhausted without success;
- the only remaining remedy needs a new parameter, soft state membership, a split of a
  single-state solution, or a data-driven initial labelling (rulings 1-3);
- a change cannot keep `nstates=1` identical;
- the `shc` route's profile falls below its floor (6 of 20) on a change the run keeps;
- the machine cannot run what a phase needs.

A stop is written up as a question for Hans with the evidence, never worked around.

## 4. Fixtures, measures and the success criterion

**Data** (read-only, on dell2):
- full set F: `/data2/multistate_bench/full80k/` (PT true poses and labels, PZ zeroed
  poses, P2 the 2D solution of PZ, `ncls=30`), made by the first run with the code of
  2026-11-03, readable by the current build;
- subset S: `/data2/multistate_bench/subset20k/` (PTS, PZS, P2S, the indices and the
  truth in row order);
- balanced subset B: 8,000 of each true state (24,000), drawn in Phase 0 the way S was
  (`make_subset.py`, a new seed), from PT and PZ by `selection prune=yes`;
- the tools at `/data2/multistate_bench/tools` (`acc3.py` and friends).

**Fixtures**, fastest first:
- **Unit and route tests**: the fast gate; `multistate_workflows` routes `labels`
  (4 min), `states`, `cavgs`, `shc`, `hard_shc`, `consensus` (5-20 min each at 8 threads).
- **Label recovery (LR)**: true poses, balanced random labels drawn once per dataset
  (seed recorded) with starting maps from `reconstruct3D nstates=3`, then
  `refine3D refine=<mode> nstates=3 pgrp=c1 mskdiam=300 lp=10 nspace=2500
  update_frac=<the stages' fraction on that set> balance=none maxits=8 objfun=euclid
  ml_reg=yes prob_beta_end=1 nthr=16 print_states=yes`. It isolates the state
  assignment from the pose search (about 1 h on S on dell2). On F the update fraction is
  0.375 (30,000 per iteration, as the stages sample), on S and B 0.9.
- **De novo (DN)**: `simple_exec prg=solve3D pgrp=c1 mskdiam=300 nstates=3 nparts=2
  nthr=8 print_states=yes projfile=<copy of P2S or P2>` (about 6 h on S, 10-14 h on F).

**Measures** per run: accuracy under the best matching of found to true states,
recall per true state, ARI, populations, per stage and per iteration
(`print_states=yes`); the state readout lines; the final maps docked onto the truth.

**Verdict** per run: *separated* when every true state is recovered (its matched found
state holds at least 0.75 of its particles) and accuracy is at least 0.90; *partial*
when accuracy is at least 0.75 or exactly two states are recovered; *failed* otherwise.

**Success criterion of the run**: the code the run ends with separates S de novo in at
least 2 of 3 runs and F de novo in at least 1 of 2, with `nstates=1` identical, the
`shc` profile at or above its floor, and the nightly passing. Read from 2026-11-08 (ruling
11): B de novo in at least 2 of 3 runs in place of S, and F in at least 1 of 2, the other
gates unchanged; S is reported as a stress case. Superseded by ruling 12 (2026-11-08): de novo
on F05 (SNR 0.05) separates in at least 2 of 3 runs, the other gates unchanged for any kept
change; the S and B criteria are recorded as superseded.

## 5. Phases

### Phase 0: builds, data, and why it fails here

1. Build HEAD in `repo/build` (`./compile_clean.sh`). Build two anchors from `git archive`
   exports in `scratch/keep/`: **A25**, tag `v3.1-beta` (2025-05-27, the paper's code:
   `abinitio3D`, `multivol_mode=independent`, `shc_smpl` stages, bounded draw), and
   **A1015**, 1bbbdcc24 (2026-10-15, before the whole-row posterior: bounded draw and
   frontier, `shc` stages 1-2, prob 3-5; five stages for `nstates>1` unless `nstages`, as
   Phase 0 found). Old code under
   CMake 4.2 may need `-DCMAKE_POLICY_VERSION_MINIMUM=3.5`; old code under gfortran
   15.2 may need small compile fixes, made in the export only and listed in the journal.
   If an anchor does not build after reasonable effort, record why and go on without it.
2. The balanced subset B (section 4). Each anchor reads its own projects: import the
   needed sets with the anchor's own build (`import_particles` with the field printouts
   as `oritab`, or whatever that build accepts; check row identity as the first run
   did), and make its own 2D solution where its ab initio needs one.
3. Diagnostics, each on its own on the machine, in this order, recorded per run:
   - LR with HEAD on F and on B (S is known: 0.69, 0.65);
   - LR with A1015 and with A25 on S (their own refine3D, the same settings as far as
     the build takes them; record what differs);
   - DN with A25 on S, twice (the paper's method, end to end); DN with A1015 on S once.
     Dropped by ruling 9 (2026-11-07): no de novo runs with the old builds. In their place, HEAD's
     LR on S at the full set's update fraction (0.375), which separates the subset's size from the
     update fraction.
4. Read the diagnostics by this table, written into Progress before Phase 1:

   | Observation | Reading | Next |
   | --- | --- | --- |
   | An anchor separates (LR or DN) where HEAD fails | a regression since that anchor | Phase 1 |
   | HEAD's LR separates F or B but not S | the subset or the imbalance made S unfair | Phase 2 targets F (and B) |
   | No build separates LR on F | hard labels from a random start do not reach the truth on this data | stop: write up (ruling 2) |
   | Otherwise | the de novo path is the problem, not the state rule | Phase 2 |

Exit: the builds (or the reasons an anchor was skipped), B, the diagnostics with their
verdicts, and the reading.

### Phase 1: locate the regression (only if Phase 0 found one)

Bisect between the anchor that separates and HEAD on the cheapest fixture that shows
the difference (LR on S if it does, else DN on S), by first-parent commits that touch the
multi-state path (`git log --first-parent -- src/main/solve src/main/strategies
src/main/commanders/simple/simple_commanders_solve3D.f90`), one export and build per
probe. The suspects to try first, in order: 811934ca5 (2026-10-24, the whole-row
posterior and the frontier deleted), the multi-state schedule of 2026-10-16/18,
dffa21364 (2026-10-06, `shc_smpl` deleted), e2673a5d7 (2026-10-15, the shift grid in the
multi-state tables). Name the change or changes and show the mechanism (per-iteration
accuracy and the labels' movement), then carry the remedy into Phase 2 as its first
candidate.

Exit: the change named with its evidence, or "no single change" with the probes.

### Phase 2: restore and build, on HEAD, for `nstates>1` only

Candidates, in the order Phase 0/1's evidence sets (default order below); each is a
change of HEAD confined to `nstates>1`, unit-tested, with `nstates=1` identical:

- **C1, the paper's early search**: restore `shc_smpl` (from dffa21364^, rewritten to
  today's search frame: the keyed draws, the polar shift map, the per-thread rules) and
  have the controller run it for `nstates>1` in the stages before the probabilistic
  search, through stage 4 as the paper did.
- **C2, the bounded probabilistic draw**: for `nstates>1`, the draw within the angular
  neighbourhood of the best candidates (`prob_athres`, as before 811934ca5), and the
  reference-balanced frontier if Phase 1 implicates it. The whole-row posterior stays
  for `nstates=1`.
- **C3, the schedule**: the paper's eight stages for `nstates>1`, and a ladder that does
  not sit at the 2D classes' resolution for four stages and then jump (on S: 24, 22.6,
  22.6, 20.2, 6 A). Per ruling 10 the ladder is tested with `lpstart=<A> lpstop=6
  force_lp_range=yes` on the command line (never `freq_march_aggressive`), and kept, if it
  separates, as a linear-ramp default of the controller for `nstates>1`.
- **C4, the per-stage annealing tail**: `prob_beta_end=1` on every multi-state
  probabilistic stage (an existing key).
- **C5, the sampling** (named by HEAD's runs on F05, ruling 13): stochastic class-balanced
  sampling in every stage of a multi-state run, where HEAD samples greedily in stages 1-3.
- What Phase 1 names, first.

For each candidate: the fast gate and the touched unit suites; the single-state
identity (`simulated_workflow_6vxx` at one thread, base export against the tree,
`test=identity`); `labels`, `shc`, `hard_shc` once each; LR on S where the candidate acts
on the probabilistic stages; then DN on S twice. A candidate is kept when it improves the
DN verdicts or the per-stage accuracy without breaking a gate; candidates may be
combined. Keep a table per candidate in Progress: what changed, the gates, LR and DN
results per stage.

Exit: the kept combination meets the S part of the success criterion (2 of 3 DN runs
separated, the third run made here), or every candidate has been tried (then the stop
condition of section 3 applies). From ruling 11 the judge is B, not S: one de novo run of
HEAD unchanged on B first (if it separates, HEAD de novo on F next), then the candidates on
B in their order, 2 of 3 B runs separated for the kept code. From ruling 12 the judge is
F05: HEAD unchanged de novo on F05 twice (and a third if they split), then label recovery
once; no code change if HEAD separates; otherwise the candidates on F05, C5 first (ruling 13).

### Phase 3: validation and robustness

- DN on F twice with the kept code.
- `simple_test_exec test=profile target=multistate_workflows suite=shc nseeds=20` (floor
  6), `suite=labels` (floor 6), the full `multistate_workflows` entry, `test=nightly
  ncunits=16`.

Exit: the success criterion of section 4, or its failure written up with the evidence.

### Phase 4: documents and report

The policy documents (`solve3D_policy.md` section 6 and the stage schedule,
`refine3D_policy.md` for any restored mode or draw), the test policy if a route changes,
proposed skill text in the report, the report beside the completed plan, the plan moved.

## 6. Common rules

- Nothing is committed; git is read-only; anchors and probes are `git archive` exports.
- One timed job at a time on the machine; record the load with every time.
- Unset `SIMPLE_SEED` for every benchmark run; tests use their own seeds.
- New or changed OpenMP code: a tester at 16 threads, repeated, bit for bit with one thread.
- No `goto`; THROW macros never continued with `//&`; no short circuit on allocatables;
  `use ..., only:`; no large type as a plain component.
- Bulky outputs are deleted once their numbers are in the journal; `scratch/keep/runs/`
  keeps every run's labels, logs, statistics and final maps.
- Progress entries and the report read without other documents; name things by what
  they do; the report ends with a commit message.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_ab_initio_recovery.md` | this plan and its progress | each |
| `src/main/strategies/search/simple_strategy3D_shc_smpl.f90` | C1: the restored hill climbing with sampled in-plane angles | 2 |
| `src/main/strategies/search/simple_strategy3D_matcher.f90` | C1: the restored mode's dispatch | 2 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | C1: the sampled in-plane neighbourhood, if the search object needs it | 2 |
| `src/main/strategies/search/simple_strategy3D_alloc.f90` | C1: per-thread state of the restored mode | 2 |
| `src/main/commanders/simple/simple_commanders_refine3D.f90` | C1: `refine3D` accepts `refine=shc_smpl` (its list of the matcher's modes) | 2 |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | C1: the run seed adopted (and exported to the distributed parts) for the keyed draws of `shc_smpl`, as for the probabilistic modes | 2 |
| `src/main/params/simple_parameters.f90` | C1: `shc_smpl` back in the `refine` key's documented values (no new key) | 2 |
| `src/main/params/simple_parameters_phases.f90` | C1: the mode's derived flags | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | C2: the bounded draw and the frontier for several states | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | C2: the same in the neighbourhood table | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | C2: the bounded sampler | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | C2: its tests | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | C1, C3, C4, C5: the stage modes, schedule, tails and sampling for several states | 2 |
| `src/main/solve/simple_solve3D_utils.f90` | C5: the interface of the deleted `solve3D_stoch_sampl_stage` | 2 |
| `src/main/solve/simple_solve3D_reseed_tester.f90` | C5: the emitted stage sampling for one and several states, tested | 2 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | C3: the multi-state stage count | 2 |
| `src/main/simple_convergence.f90` | a stage rule the candidates need | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90` | path checks of the routes when the stage modes change | 2 3 |
| `doc/policies/3D/solve3D_policy.md` | section 6 and the schedule | 4 |
| `doc/policies/3D/refine3D_policy.md` | a restored mode or draw | 4 |
| `doc/policies/importance_sampling_fractional_update_policy.md` | C5: the multi-state sampling switch | 4 |
| `doc/policies/test_environment_policy.md` | a route's description if it changes | 3 4 |
| `doc/implementation_notes/completed/multistate_ab_initio_recovery_report.md` | the report | 4 |
| `doc/implementation_notes/completed/multistate_ab_initio_recovery.md` | the plan, moved at the end | 4 |

Other files a phase needs are added to this table, with the reason, by that phase.

## Progress

(one section per phase, written by the run)

### Phase 0 (2026-11-07, dell2): builds, the balanced subset, the diagnostics and their reading

`R` below is `/home/hael/agent_runs/ms_abinitio`; the journal is `R/state/journal_phase_0.md`,
the scripts are in `R/scratch/keep/scripts/`, every run's labels, logs and statistics are in
`R/scratch/keep/runs/<run>/`, and `R/scratch/keep/runs/summary.tsv` holds one line per run.
Files changed: this plan only. Accuracy is the share of particles whose found state is matched to
their true state under the best one-to-one matching of found to true states; the recall of a true
state is the share of its particles that its matched found state holds; ARI is the adjusted Rand
index (0 for labels unrelated to the truth, 1 for the truth up to renaming). The true states are
the smallest assembly (state 1), the medium (state 2) and the largest (state 3).

**Builds.** All three with Ubuntu's gfortran 15.2 and CMake 4.2.3, the compilers of the
maintainer's own build trees on dell2 (Release, OpenMP, no environment modules):
- HEAD (d4e89630c) in `repo/build` by `./compile_clean.sh`: 4 minutes, the fast gate 15 of 15,
  installed.
- **A25**, tag `v3.1-beta` (bd6f82ddd, 2025-05-27), exported to `R/scratch/keep/src_A25` and
  built in `build_A25` with `-DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DGUI=off -DUSE_LIBTIFF=off
  -DBUILD_TESTS=OFF` and the install prefix set to the build tree (its cache otherwise keeps
  `/usr/local`).
- **A1015**, 1bbbdcc24 (2026-10-15), exported to `src_A1015` and built in `build_A1015` with
  `-DBUILD_TESTS=OFF`.

Neither anchor needed a source patch: both compile and link with gfortran 15.2 as they are.
Only their executables are needed, so their tests were not built.

**What the anchors are** (read from their sources; the plan's section 1 needed two corrections):
- *A1015 runs five multi-state stages, not eight.* Its multi-state `solve3D` already defaulted to
  `nstages=5` and `lpstop=6` A. That default dates from 32a324ab9 (2026-06-08), not from
  2026-10-16/18 as section 1's table says. The 2026-10-16/18 commits changed the stages' content:
  `prob_neigh` at stage 5 for a few days, NU filtering from stage 5, and the stage-2 scan exit.
  So the A1015 de novo run below uses its own five-stage default, as Cyril's runs of that build
  did.
- *A1015's refine3D has no `prob_beta` or `prob_beta_end`* (an unknown key stops a SIMPLE program).
  Its state rule is the greedy reference frontier of `ref_assign`. Every reference points at its
  closest unassigned particle, and the reference with the smallest such distance takes that
  particle. The particle's state is that reference's: in effect the closest reference over all
  states. The direction is a draw among the best frontier references within the `prob_athres`
  window, weighted by `exp(-(d - dmin))`, with no temperature.
- *A25 (the paper's code) runs `abinitio3D`.* For `nstates>1` it selects
  `multivol_mode=independent`. The random start draws random orientations and a balanced
  labelling ('squared_uniform', one of each label per block of `nstates` particles), then noise
  volumes. The schedule:
  - eight stages over 500, 500, 1,000, 1,000, 1,000, 1,000, 2,500 and 2,500 directions;
  - `shc_smpl` in stages 1-4 (hill climbing with the in-plane angle sampled within a
    neighbourhood) and `prob` in stages 5-8;
  - class-balanced sampling from the 2D solution in every stage;
  - `objfun=euclid` with one global noise spectrum (`sigma_est=global`) and `prob_athres=10`.

  Its refine3D has no `prob_beta`, `prob_beta_end` or `print_states`, and writes no
  per-iteration label file. In its runs, a watcher copies the project file each time the program
  rewrites it, and the copies' labels are printed afterwards with A25's own
  `print_project_field`.

**The data.**
- *B, the balanced subset.* `R/scratch/keep/scripts/make_subset_B.py` draws 8,000 particles of
  each true state (24,000 in all) the way S was drawn: `random.Random(20261107).sample` within
  each state. Since the benchmark has only 8,000 of the largest assembly, B holds all of them.
  B shares 5,939 particles with S. The outputs are in `R/scratch/keep/subsetB/`.
- *PTB and PZB* were made with HEAD from copies of PT and PZ by `selection oritype=ptcl3D
  infile=subsetB_flags.txt prune=yes`. HEAD's pruning set every PTB row to state 1 (PT itself
  holds 24,000 / 48,000 / 8,000), so the true labels were put back by `selection
  infile=subsetB_truth.txt`.
- *Row checks.* `check_subset_rows.py` checks each row against the source project's printout:
  defocus, astigmatism, angles, shifts, half, image and particle indices, and for PTB the true
  label. It found 0 mismatches in 24,000 rows.
- *Each anchor's own projects.* `new_project` then `import_particles
  stk=/data2/multistate_bench/particles.mrcs oritab=<PT_field_import.txt | PZ_field_final.txt>
  ctf=yes smpd=3 kv=300 cs=2.7 fraca=0.07`; then S by `selection infile=subset_flags.txt
  prune=yes`, and for the true-pose project the true labels by `selection infile=subset_truth.txt`.
  Both builds read the printouts as they are and both can prune. Every anchor project has 0
  mismatches on every key:
  - the full A25 true-pose project against its printout (80,000 rows);
  - the four anchor subset projects against theirs (20,000 rows each).
- *Projects that need care.* `print_project_field` fails with an I/O error on the read-only
  `/data2` projects: work on copies.

**The label-recovery (LR) starts.** `draw_balanced_labels.py` shuffles a list of balanced labels
with Python's `random.Random(seed)`:
- S: seed 20261106, run 2's start, reproduced byte for byte (accuracy 0.338);
- F: seed 2026110701 (0.335);
- B: seed 2026110702 (0.337).

`lr_prep.sh` applies the labels with `selection` and reconstructs the starting maps with
`reconstruct3D nstates=3 objfun=cc`, each build its own. The anchors' runs on S start from run 2's
labelling, so they are paired with run 2's two HEAD runs on S (0.685 and 0.648).

The LR command, HEAD: `refine3D refine=prob nstates=3 pgrp=c1 mskdiam=300 lp=10 nspace=2500
update_frac=<0.375 on F, 0.9 on S and B> balance=none maxits=8 objfun=euclid ml_reg=yes
prob_beta_end=1 nthr=16 print_states=yes iterstats=yes vol1..3=<the starting maps>`. The anchors
take the same command without `prob_beta_end`, which neither has. A25 also takes no
`print_states`, and `balance=no` (its key is yes|no).

**Diagnostic 1: HEAD's LR on F separates.** `lr_F_head`, 1 h 36 min, 5.2 GB; load about 16, and
16-25 during the first ten minutes while the anchors' starting maps were made. Each iteration
updates 30,000 particles. The sample cycles through the set, so after 8 iterations every particle
has been updated exactly three times.

| iteration | factor | accuracy | recall smallest / medium / largest | state n_eff at the factor / at 1 |
| --- | --- | --- | --- | --- |
| 1 | 0.30 | 0.338 | 0.34 / 0.33 / 0.35 | 2.60 / 1.67 |
| 2 | 0.30 | 0.361 | 0.37 / 0.34 / 0.44 | 2.50 / 1.58 |
| 3 | 0.30 | 0.449 | 0.48 / 0.40 / 0.65 | 2.05 / 1.37 |
| 4 | 0.30 | 0.615 | 0.67 / 0.56 / 0.80 | 1.28 / 1.08 |
| 5 | 0.41 | 0.826 | 0.87 / 0.79 / 0.94 | 1.01 / 1.00 |
| 6 | 0.55 | 0.925 | 0.95 / 0.90 / 1.00 | 1.00 / 1.00 |
| 7 | 0.74 | 0.962 | 0.97 / 0.95 / 1.00 | 1.00 / 1.00 |
| 8 | 1.00 | 0.9985 | 0.995 / 1.000 / 0.9995 | 1.00 / 1.00 |

The run ends separated: ARI 0.995, found populations 7,996 / 48,114 / 23,890 against the true
8,000 / 48,000 / 24,000. Here the rare largest assembly gains its own state first. The unrefined
share of each state block's mass stays at or below 0.036.

**Diagnostic 2: HEAD's LR on B separates.** `lr_B_head`, 1 h 05 min, 4.6 GB, load about 16.
Accuracy per iteration: 0.353, 0.443, 0.683, 0.830, 0.967, 0.996, 0.998 and 0.998. The run ends
with ARI 0.994, recall 0.994 / 1.000 / 1.000 and populations 7,952 / 8,048 / 8,000. The smallest
and largest assemblies organise first (recall 0.83 and 0.88 at iteration 3) while the medium's
state shrinks to 4,295 particles. The medium then recovers, and the set is separated at
iteration 5, before the annealing tail has done much (factor 0.41). State n_eff at 0.3: 1.96,
1.87, 1.26 and 1.14 over iterations 1-4.

**Diagnostic 3: A1015's LR on S is partial and still climbing.** `lr_S_A1015`, 1 h 46 min,
4.7 GB, load about 16. Accuracy per iteration: 0.354, 0.424, 0.554, 0.629, 0.656, 0.705, 0.761
and 0.834. The run ends with ARI 0.63 and recall 0.991 / 0.729 / 0.999. Through iteration 4 it
follows run 2's HEAD runs from the same start (HEAD at iteration 4: 0.629 / 0.667, recall of the
largest 0.20 / 0.30; A1015: 0.629, 0.08). Then the paths part:
- A1015 holds 1,773-1,998 of the 2,000 largest-assembly particles in one state together with
  medium particles (6,969 at iteration 4, 6,609 at 5, 3,256 at 8).
- Over iterations 5-8 that state sheds its medium particles to the pure medium state (4,880 ->
  8,744). The medium's recall rises 0.45 -> 0.73 and the accuracy still gains 0.07 per iteration
  at the end.
- In HEAD's runs the largest assembly's state kept the larger part of the medium (7,931 and
  7,074 at the end) and the partition froze under the annealing tail.

**Diagnostic 4: A25's LR on S separates.** `lr_S_A25`, 1 h 56 min, 18.9 GB, load about 16. Its
labels come from copies of the project file taken after each iteration. The copy after the last
iteration was missed, so those labels were read from the run's own project file. The run starts
from the same labels as A1015 and run 2's HEAD runs.

| iteration | A25 | A1015 | HEAD (run 2, two runs) | recall of the largest: A25 / A1015 / HEAD |
| --- | --- | --- | --- | --- |
| 1 | 0.357 | 0.354 | 0.350, 0.357 | 0.34 / 0.32 / 0.33, 0.33 |
| 2 | 0.396 | 0.424 | 0.383, 0.414 | 0.33 / 0.26 / 0.30, 0.31 |
| 3 | 0.480 | 0.554 | 0.478, 0.553 | 0.32 / 0.19 / 0.26, 0.31 |
| 4 | 0.650 | 0.629 | 0.629, 0.667 | 0.32 / 0.08 / 0.20, 0.30 |
| 5 | 0.727 | 0.656 | 0.691, 0.689 | 0.41 / 0.99 / 0.16, 0.27 |
| 6 | 0.770 | 0.705 | 0.699, 0.683 | 0.63 / 1.00 / 0.09, 0.19 |
| 7 | 0.835 | 0.761 | 0.691, 0.668 | 0.96 / 1.00 / 0.02, 0.08 |
| 8 | **0.927** | 0.834 | 0.685, 0.648 | 1.00 / 1.00 / 0.00, 0.01 |

The largest assembly's recall moves differently in the two anchors:
- *In A1015 it jumps at iteration 5*, where the matching of found to true states changes. By
  iteration 4 the largest assembly already sits almost whole (1,773 of 2,000) in a state with
  6,969 medium particles. From iteration 5 that state is the one matched to it.
- *In A25 the matching never changes.* The largest assembly gathers step by step into the
  smaller of the two medium-holding states: 631, 818, 1,260, 1,915 and 1,993 of its particles
  at iterations 4-8. Over the same iterations that state's medium particles fall: 4,024, 4,015,
  3,811, 3,201, then 1,445.

A25 ends separated: ARI 0.81, recall 0.999 / 0.880 / 0.997. One state holds 5,992 particles of
the smallest assembly and nothing else. A second holds 10,555 of the medium. The third holds 1,993
of the largest together with 1,445 medium particles, down from 3,201 one iteration before, and
the accuracy was still rising by 0.09 per iteration.

Both anchors take the same path. The smallest assembly separates first. The largest stays in a
state with part of the medium, and from iteration 5-6 that state loses its medium particles to
the pure medium state. In HEAD's runs it is the other way round: the largest assembly's state
keeps the larger part of the medium, the remaining medium particles form the third state, and the
largest assembly ends inside a medium state (recall 0.00-0.01). A25 is somewhat faster than A1015
on that path: at iteration 7 it is at 0.835 against 0.761.

Two differences between the runs:
- *No annealing tail in the anchors.* Neither anchor has the key, so neither had the tail of
  HEAD's runs, which from iteration 4 raises the factor to 1. A higher factor makes HEAD's draw
  more like the anchors' closest-state rule, not less.
- *Different sampling of the 18,000 particles per iteration.* A25 samples `balance=no`, the others
  `balance=none` (the same: no class balancing). Each build otherwise draws its sample its own
  way.

**Diagnostic 5, added: HEAD's LR on S at the full set's update fraction is partial and loses the
largest assembly.** `lr_S_head_uf0375`, 25 minutes, 3.6 GB, load about 16. HEAD separated F at an
update fraction of 0.375 and failed S at 0.9, so F and S differed in two ways at once: the set (4
times as many particles per state) and the update fraction. Each iteration at 0.375 reassigns
fewer particles and so moves the maps less. This run takes S, run 2's start and the update
fraction of F: 7,500 particles per iteration, every particle updated three times in the eight
iterations. It replaces the dropped de novo runs (ruling 9).

Accuracy per iteration: 0.342, 0.383, 0.504, 0.632, 0.770, 0.831, 0.855 and 0.884. The run ends
with ARI 0.73 and recall 0.948 / 1.000 / 0.000:
- the smallest and the medium assemblies separate;
- the largest is never gathered: its recall falls from 0.31 at iteration 2 to 0.00, and it ends
  inside the medium state (1,995 of its 2,000 particles, with all 12,000 medium ones);
- the third found state drains to 221 particles, and its map's FSC=0.5 resolution falls to 192 A.

It counts as partial only because two true states are recovered. So on S HEAD loses the rare
assembly at either update fraction: the update fraction is not what fails on S.

**The diagnostics in one table** (LR = label recovery: true poses, a balanced random start, eight
iterations, `nspace=2500`, `lp=10`):

| run | build | set | update fraction | accuracy | ARI | recall smallest / medium / largest | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `lr_F_head` | HEAD | F | 0.375 | 0.9985 | 0.995 | 0.995 / 1.000 / 1.000 | separated |
| `lr_B_head` | HEAD | B | 0.9 | 0.998 | 0.994 | 0.994 / 1.000 / 1.000 | separated |
| run 2, two runs | HEAD's code | S | 0.9 | 0.685, 0.648 | 0.48, 0.48 | 0.96-0.98 / 0.59-0.66 / 0.00-0.01 | failed |
| `lr_S_head_uf0375` | HEAD | S | 0.375 | 0.884 | 0.73 | 0.948 / 1.000 / 0.000 | partial (the largest lost) |
| `lr_S_A1015` | A1015 | S | 0.9 | 0.834 | 0.63 | 0.991 / 0.729 / 0.999 | partial, still rising |
| `lr_S_A25` | A25 | S | 0.9 | 0.927 | 0.81 | 0.999 / 0.880 / 0.997 | separated, still rising |

The de novo runs of A25 and A1015 were not made (ruling 9). Neither were the anchors' 2D
solutions, which only those runs needed.

**The reading** (Phase 0's table):

| Observation | Holds? | Reading | Next |
| --- | --- | --- | --- |
| An anchor separates (LR or DN) where HEAD fails | **Yes.** A25 separates S in LR (0.927) from the start on which HEAD failed twice (0.685, 0.648). A1015 takes the same path and is at 0.834 and rising after eight iterations. | A regression since A25. A1015 already takes the anchors' path (the largest assembly's state sheds its medium particles), so the change that matters probably lies between 1bbbdcc24 and HEAD. | Phase 1 |
| HEAD's LR separates F or B but not S | **Yes.** F 0.9985 and B 0.998 against S 0.65-0.88. | The subset is harder. S holds 2,000 particles of the largest assembly; F and B hold 8,000. The imbalance alone does not fail (F is 30 / 60 / 10 % and separates), and neither does the update fraction (diagnostic 5). But S is not unfair: the paper's code separates it, so S is the fixture that tells the builds apart. | Phase 2 targets F and B as well as S |
| No build separates LR on F | No: HEAD separates F. | | no stop |
| Otherwise | | | |

What this means for the next phases:
- *Hard labels from a random start can reach the truth on this data* (ruling 2). HEAD does it on
  F and B, A25 on S. Nothing here calls for soft state membership.
- *The state rule on true poses has regressed for a rare state with few particles.* On the same
  start, the closest-state rule of the anchors keeps the largest assembly gathered with part of
  the medium, and that state then sheds the medium. HEAD's draw from the tempered posterior over
  the whole row sends the largest assembly into a medium state. The partition freezes there, under
  the annealing tail at 0.9 and also at 0.375.
- *Two differences are not yet separated:* the state rule itself, and the anchors' lack of an
  annealing tail. As section 1 records, run 2 found the old closest-state rule and the joint draw
  alike on S; Phase 1's probes settle it.
- *Bisecting.* For Phase 1 the cheapest fixture that shows the difference is LR on S from run 2's
  start, and A1015 can serve as the good end:
  - it is within a month of HEAD, against 17 months for A25;
  - it shows the anchors' signature by iteration 5-6, so a probe can stop after iteration 6
    (phase notes).

  The signature is whether the largest assembly ends gathered in one state, with that state
  losing medium particles, or spread into a medium state with its recall falling towards 0. The
  suspects after 1bbbdcc24 in section 5's order:
  - 811934ca5 (2026-10-24: the whole-row posterior, the frontier deleted);
  - the schedule of 2026-10-16/18 (it does not act on LR);
  - e2673a5d7 (2026-10-15: the shift grid in the multi-state tables, after 1bbbdcc24);
  - the joint draw of 2026-11-05.

  dffa21364 (`shc_smpl` deleted) precedes A1015 and does not act on LR. It remains a de novo
  candidate (C1).
- *De novo is a separate question.* LR on true poses separating F does not mean de novo
  separates F: the first run's de novo run of F, on older code, failed at 0.45. Phase 2 judges
  the candidates on de novo runs, as planned.

**What later phases can reuse** (all under `R/scratch/keep/`):
- the anchors' exports and build trees (`src_A25`, `build_A25`, `src_A1015`, `build_A1015`);
- their projects (`projects_A25`, `projects_A1015`: PT, PZ, PTS, PZS). Probes after 1bbbdcc24
  should read A1015's, since later code reads earlier projects;
- B (`subsetB/`, `projects/PTB`, `projects/PZB`);
- the label starts (`labels/`) and the LR start directories with their starting maps
  (`starts/S_A1015`, `S_A25`, `S_head`, `F_head`, `B_head`);
- the scripts:
  - `lr_prep.sh`, `lr_run_head.sh`, `lr_run_anchor.sh`;
  - `anchor_import.sh`, `anchor_subset.sh`;
  - `watch_project.sh` and `snap_labels.sh` (labels of builds without `print_states`);
  - `archive_run.sh`;
  - `tools/p4_iter_table.py --set S|B|F` (per-iteration accuracy beside the readout).

**Exit items.**
- *The builds.* HEAD, A25 and A1015 were built, none needing a source patch; both anchors run.
- *B.* Drawn, its projects made, its rows checked.
- *The diagnostics with their verdicts.* LR on F and B with HEAD (both separated), LR on S with
  A1015 (partial) and A25 (separated), and the added HEAD LR on S at 0.375 (partial, the largest
  lost). The de novo runs of the anchors were dropped by ruling 9.
- *The reading.* Above: a regression since A25, with A1015 already on the anchors' path, so
  Phase 1 follows. Phase 2 also judges on F and B.

Evidence: `R/scratch/keep/runs/{lr_F_head, lr_B_head, lr_S_A1015, lr_S_A25,
lr_S_head_uf0375}/` (logs, per-iteration labels, `iter_table.tsv`, statistics, final maps);
`R/scratch/keep/runs/summary.tsv`; build logs `R/scratch/logs/build_{head,A25,A1015}.log`.

### Phase 1 (2026-11-07/08, dell2): looking for the regression in label recovery on S

`R` is `/home/hael/agent_runs/ms_abinitio`; the journal is `R/state/journal_phase_1.md`. Files changed:
this plan only. The section 1 sentence that the joint draw and the old argmin "behave alike" on S is
corrected in place. Every run below is label recovery (LR) on S from run 2's random start (seed
20261106), with the anchors' command: `refine3D refine=prob nstates=3 pgrp=c1 mskdiam=300 lp=10
nspace=2500 update_frac=0.9 balance=none objfun=euclid ml_reg=yes nthr=16 print_states=yes` with the
starting maps, and no annealing tail. Each run's labels, log, statistics and per-iteration tables
are in `R/scratch/keep/runs/<run>/`, with one line per run in `R/scratch/keep/runs/summary.tsv`.
Each probe is a `git archive` export built in `R/scratch/keep/build_<name>` without tests
(`probe_build.sh`).

**The signature.** Accuracy after four iterations does not tell the builds apart: it was 0.63-0.67
in A25, A1015 and HEAD alike. The difference lies in what happens to the largest assembly, and it
shows from iteration 5:
- By then the smallest assembly has its own state. The medium assembly is spread over a big state
  (about 7,000 of its 12,000 particles) and a smaller one (about 4,000). The largest assembly's
  2,000 particles are split between those two.
- *Good path* (the anchors in Phase 0): the largest assembly gathers in one of the two states, and
  that state then loses its medium particles to the other. Its map becomes the largest's, and the
  accuracy climbs to 0.83-0.96.
- *Bad path* (HEAD in run 2): the largest assembly ends inside the big medium state, whose map stays
  a medium map. Medium particles keep flowing into that state, the largest's recall falls towards
  0, and the accuracy stays at 0.65-0.80.

`fork_class.py` names the path at the last iteration from the state that holds most of the largest
assembly:
- GOOD: that state holds fewer than 5,000 medium particles, and fewer than two iterations before;
- BAD: it holds 6,500 or more.

`label_flows.py` counts the medium particles moving between states from one iteration to the next.
On the bad path (run 2's HEAD run, iterations 4-7), 366, 242 and 137 move into the largest's state
per iteration, against 93, 162 and 231 out of it. On the good path (A1015, the same iterations),
101, 42 and 22 move in, against 513, 879 and 1,094 out.

**The probes, in the order run** (the start is A1015's project and maps unless marked "HEAD's
start"; a run stopped early was stopped once its path was clear, as the phase notes allow):

| run | build | path | at the last iteration |
| --- | --- | --- | --- |
| control | HEAD, no annealing tail, HEAD's start | BAD | 0.798 at 8; the largest's state gains medium, 8,095 -> 10,029 |
| HEAD's parent | d592c90e9 (2026-11-06, the closest-state rule), HEAD's start | BAD | stopped at 6: 0.726; medium in the largest's state 7,281 -> 8,123 |
| 811934ca5 | 2026-10-24: the posterior draw on every table, the frontier deleted | GOOD | 0.945 at 8, separated |
| 9c2fcdaf9 | 2026-10-31 | BAD | stopped at 6: 0.740 |
| 707dfa130 | 2026-10-27 | GOOD (slow) | 0.773 at 8; the largest gathers from iteration 7 |
| 16b0cfc49 | 2026-10-29: after the shift prior's default and the weighted insertion | GOOD | stopped at 6: 0.785 |
| 6ddc4a0d2 | 2026-10-30 | GOOD | stopped at 6: 0.781 |

HEAD's parent, the last build before the joint draw (its label is the closest reference over all
states), fails. So the joint state draw of 2026-11-05 is not the regression. The control fails without the tail, so
the annealing tail is not either. 811934ca5 separates, so the posterior draw and the frontier's
deletion are not.

The bisection then closed on 6ddc4a0d2 (good) and 9c2fcdaf9 (bad). Seven first-parent commits
separate them, and none touches what label recovery runs:
- a new polar shift map module and its tester (519c67d2a);
- new procedures and components of the polar Fourier calculator that `refine3D` does not call yet
  (a2c1b1f86);
- plan notes;
- a half-set fix in `make_cavgs` and nanoparticle code (653c5cd21);
- a code review that touches documents only;
- the module file the build was missing (9c2fcdaf9).

The commits a2c1b1f86..00dc4ec55 do not build as committed (an undefined `gen_shiftmap`: the
module file came only with 9c2fcdaf9), so that window was not probed. The two ends therefore run the
same label-recovery code and took opposite paths.

**Label recovery on S is a stochastic fork.** Replicates of the same build from the same start go
either way. Every run seeds itself from `/dev/urandom`; `SIMPLE_SEED` is unset as the plan's common
rules require.

| build | runs on the good path / runs |
| --- | --- |
| A25 (the paper's code, its own start) | 2 / 2 (Phase 0's run 0.927; the replicate 0.837 after 8 iterations) |
| A1015 (2026-10-15) | 1 / 2 |
| 811934ca5 (2026-10-24) | 1 / 1 |
| 707dfa130 (2026-10-27) | 1 / 1 |
| 16b0cfc49 (2026-10-29) | 1 / 2 |
| 6ddc4a0d2 and 9c2fcdaf9 (the same label-recovery code) | 1 / 3 |
| HEAD's parent | 0 / 1 |
| HEAD | 2 / 6: 0 / 3 from HEAD's start (run 2's two runs with the tail and the control), 2 / 3 from A1015's start |

HEAD separates S in two runs: 0.958 and 0.935 after seven iterations, the largest assembly
gathered whole (1,990 of 2,000) with 691 and 1,178 medium particles left in its state.

The two starts were compared to see whether HEAD's 2 of 3 against 0 of 3 says more than chance:
- the starting half-maps and merged maps correlate at 1.0000;
- every particle record matches in every field;
- every other project segment matches except path metadata;
- the FSC files are byte-identical.

HEAD's start directory also carries provenance and support sidecar files. Only the `refine3D_local`
start reads the first, and only a postprocessing command the second. Plain `refine3D` reads neither,
and the iteration-1 logs from the two starts differ only by the random draws. The difference is
chance (one-sided Fisher exact test, p = 0.2).

**What decides the fork.** At iteration 4 the outcome can mostly be read off. Take the ratio of the
largest assembly's particles in the big medium state to those in the small one:
- the good runs had 0.6-1.4 (811934ca5 1.23, 707dfa130 1.37, 16b0cfc49 1.43, 6ddc4a0d2 1.27,
  HEAD 0.61 and 0.81, A25 1.43 and 1.19);
- the bad runs had 1.6-3.1, with one exception: A1015's second run had 1.38 and still failed;
- A1015's first run took the other route: the big state gathered almost all of the largest assembly
  (1,773 at iteration 4, a fifth of its content) and lost its medium instead.

The run succeeds when one state captures enough of the 2,000 particles early for its map to lean
towards the largest assembly. It fails when they stay diluted at a tenth of a medium state. That
also explains why F and B, with 8,000 particles of the largest assembly, separate reliably (Phase
0): wherever those particles gather, they dominate.

**The reading.**
- *No single change.* The bisection's ends run the same code, and the rate of the good path cannot
  be told apart between the builds with the runs made. Before 2026-11-06 it was 7 of 11 runs (A25 2
  of 2); HEAD and its parent are 2 of 7. That difference is not significant (one-sided Fisher
  exact test, p = 0.17), though it leans the old code's way. The anchor's separation of S in
  Phase 0 and run 2's two failures from one start were draws of the same fork.
- *Phase 0's reading is revised.* Its first row ("an anchor separates where HEAD fails: a
  regression") rested on single runs. With replicates, HEAD separates S in label recovery as the
  anchors do. Its rate cannot be told apart from theirs with these runs, though these runs cannot
  rule out that the older code's rate was higher (2 of 7 against 7 of 11). The table's last row applies: the
  de novo path is the problem, not the state rule. Phase 0's second row still holds: HEAD separates
  F and B reliably, and S is the hard case because of its 2,000-particle rare state.
- *For Phase 2* there is no named remedy to carry. The candidates acting on the de novo path come
  first, in this order:
  - C1, the paper's early search with sampled in-plane angles;
  - C3, the schedule;
  - C4, the per-stage tail.

  C2, the bounded draw, acts on the state assignment, which label recovery does not show to be
  broken; it comes last. Every judgement on S must rest on several runs: one label-recovery or de
  novo run can land on either side of the fork. The de novo verdicts already follow that (ruling 7).
- *A measure for Phase 2.* The fork happens during the first probabilistic iterations. A candidate
  is better on S if it raises the share of runs whose largest-assembly state wins the race. The
  ratio at iteration 4 and `fork_class.py` read that in the label-recovery fixture.

**Exit.** "No single change", with the probes: 6 probe builds, HEAD's control without the tail,
and 7 replicates (HEAD three times, A25, A1015, 16b0cfc49 and 9c2fcdaf9 once each). The joint draw, the annealing tail, the posterior draw with the
frontier's deletion, the shift prior's default and the weighted insertion are each excluded by a
probe on the good path after them or by HEAD's parent on the bad path. The mechanism is shown by
the medium flows and the iteration-4 configuration of 18 runs.

Evidence:
- `R/scratch/keep/runs/` holds `lr_S_head_notail`, `lr_S_pHm1`, `lr_S_p811`, `lr_S_p1027`,
  `lr_S_p1029`, `lr_S_p1030`, `lr_S_p1031` and `lr_S_rep_*`;
- the probe builds (15 exports, 9 of them pre-built for branches the bisection did not take) were deleted
  after the phase. `R/scratch/keep/probe_builds_index.txt` lists their commits, and
  `R/scratch/keep/scripts/probe_build.sh <rev> <name>` rebuilds any of them;
- the tools are in `R/scratch/keep/tools/`: `largest_traj.py`, `label_flows.py`, `fork_class.py` and
  `mrc_corr.py`.

### Phase 2 (2026-11-08, dell2): the candidates

`R` is `/home/hael/agent_runs/ms_abinitio`; the journal is `R/state/journal_phase_2.md`. The
candidates follow the order Phase 1's reading set: C1 (the paper's early search) first.

**C1: hill climbing with sampled in-plane angles in multi-state stages 1-4.** `refine=shc_smpl`
returns as a value of `refine`. Before 2026-10-06 (dffa21364) it was the early search of the
paper's multi-state ab initio. The restored mode is rewritten to today's search frame:
- *The search, as `shc`'s.* The previous pose is prepared, the shift is searched on the previous
  reference, and the references are then visited in random order until the first one scores
  better than the previous. The chosen reference's in-plane angle and shift are polished as in
  `shc`, with the polar shift map when the key is on, and the pose is committed.
- *The in-plane angle is drawn, not taken.* For each reference it is drawn within the
  neighbourhood of its best angles, by the paper-era bounded rule (`sample_bounded_dists`). The
  window is `prob_athres` (10 degrees), capped by the state's mean in-plane distance.
- *Both draws are keyed.* The reference order (`keyed_candidate_order`, the previous reference
  last) and the in-plane draw (`rnd_keyed` on the pass, the particle and the reference) are keyed,
  so a particle's search does not depend on the thread that runs it. The matcher fixes the pass's
  key scope before the particle loop. `refine3D` adopts the run seed for it as for the
  probabilistic modes, so distributed parts read the same seed.
- *The greedy fallback.* A particle never searched, or one on the shift map's disk, scores every
  reference, as under `shc`. There are no random greedy turns.

The `solve3D` controller runs `shc_smpl` in stages 1-4 of a multi-state run from particles, and
`prob` in stage 5 with its annealing tail and the final assignment. Stages 2-4 end by the scan rule
stage 2 already had: the mean share of the search space scanned above 95 %. `solve3D_cavgs` keeps
its own schedule (the mode policy now knows which route asks). `nstates=1` is untouched.

Files:
- `simple_strategy3D_shc_smpl.f90` (new) and the matcher (`simple_strategy3D_matcher.f90`);
- the search specification (`simple_strategy3D_srch.f90`), which now carries the key scope;
- `sample_bounded_dists` (`simple_eul_prob_tab_utils.f90`), with an optional keyed variate;
- the `refine` key's documentation (`simple_parameters.f90`) and the shift map's list of modes
  (`simple_parameters_phases.f90`);
- `refine3D`'s list of accepted modes (`simple_commanders_refine3D.f90`) and the seed adoption
  (`simple_refine3D_strategy.f90`), both added to the file table;
- the stage controller (`simple_solve3D_controller.f90`);
- `multistate_workflows` (`simple_commanders_test_highlevel_msworkflows.f90`): the `shc` route's
  path check, and a new by-hand suite `shc_smpl_threads`.

The new mode is not offered in `refine3D`'s UI choices. It is internal, like `greedy`, `eval` and
`sigma`, so no metadata the GUI reads changes.

Gates (scratch logs under `R/scratch/`):

| gate | result |
| --- | --- |
| build and fast gate | 15 of 15 (`logs/c1_fastgate.log`) |
| single-state identity: `simulated_workflow_6vxx` at one thread, base export against the tree | 186 of 186 declared outputs identical (`identity_c1/identity.log`) |
| `shc_smpl_threads` (new): one `refine3D refine=shc_smpl` iteration over 3,000 particles at one thread and twice at 16 threads, on the `labels` start, `objfun=cc` and the window shift search | 2,405 particles searched by `shc_smpl`, 0 differing records one thread against 16, 0 between the two 16-thread runs (`routes/c1_shc_smpl_threads/`) |
| `labels` | pass, 3 min 20 s |
| `shc` | pass, 17 min 24 s; path `shc_smpl` 500, `shc_smpl` 1,000, `prob` 2,500 directions; every true state resolved, accuracy 0.984 |
| `hard_shc` | pass (its labels reported, not gated), 11 min 44 s; the conformational pair not resolved (accuracy 0.527) |

The thread check runs at `objfun=cc` with the window shift search for two reasons:
- under `euclid` the noise spectra are estimated over the team, whose summation order varies
  with the thread count;
- on the shift map's disk every particle would take the greedy search, so the restored mode would
  not run.

So the check compares only what the restored mode decides. As in run 2, `hard_shc` resolves
nothing at its seed; that route is a robustness profile.

**C1 de novo on S.** The runs use `simple_exec prg=solve3D pgrp=c1 mskdiam=300 nstates=3
nparts=2 nthr=8 print_states=yes`, each on its own copy of the subset's 2D solution P2S, with
`SIMPLE_SEED` unset. Evidence is in `R/scratch/keep/runs/<run>/` (`metrics.txt`, the labels of
every stage and iteration, the log), with a line per run in `summary.tsv`. Accuracy after each
stage:

| run | stage 1 | stage 2 | stage 3 | stage 4 | stage 5 | what happened |
| --- | --- | --- | --- | --- | --- | --- |
| C1, first version | 0.443 | 0.444 | 0.454 | 0.467 | - | stopped in stage 5's first iteration: stages 3-4 had ended after 3 iterations each (below) |
| C1 (`dn_S_c1_1`) | 0.409 | 0.414 | 0.407 | 0.415 | 0.404 at its 6th iteration | stopped: flat through half of stage 5 |
| C1 and the forced ramp (`dn_S_c13_1`) | 0.434 | 0.424 | 0.436 | 0.447 | - | stopped after stage 4 |
| run 2, HEAD (three runs) | 0.39-0.45 | 0.39-0.45 | 0.37-0.40 | 0.39-0.51 | 0.42-0.59 | failed |

*The first version's stage exit.* It also ended stages 3-4 on stage 2's scan rule (the mean share
of the search space scanned above 95 %), and that rule measured the wrong thing there. On entering
stage 3 every particle is uncentred, so the polar shift map sends it to the exhaustive search,
which scans the whole space. The stages therefore ended after 3 iterations each. That rule was
removed for stages 3-4 (stage 2 keeps it), the build and fast gate were redone, and the `shc`
route was run again: it passes, every state resolved, accuracy 0.981.

*Why the runs were stopped.* Stage 5 of `dn_S_c1_1` stayed at 0.404-0.414 for six iterations.
Run 2's stage 5 added at most 0.10, so that run could not have reached even partial.

*The reading.* In all three runs the hill climbing with sampled in-plane angles leaves the labels
at 0.41-0.47 after stage 4, where run 2's stages left them at 0.39-0.51. At 24-20 A it does not
organise the states. With the forced ramp (ruling 10: 20, 16.5, 13, 9.5, 6 A) its stages at 13 and
9.5 A converge within 6 iterations each: the hill climbing freezes the labels once no reference
improves (state overlap above 0.96 and 99 % scanned), still at 0.44.

C1 does not improve the per-stage accuracy, so it is not kept alone. Its second pure run and the
second run with the ramp were not made: three runs show the same stall, and ruling 10 directs the
ramp to be tested without C1. That test is running.

**C1's code is reverted.** C1 is not kept, so its diff is saved as
`R/scratch/keep/variants/c1_shc_smpl.patch` and every file it touched is restored from HEAD. The
two file-table rows added for it stay, with their reason. Its gates above remain the evidence that
the restored mode works as designed.

**C3, the forced linear ramp, without C1** (ruling 10; the base build, i.e. HEAD exported, with
`lpstart=20 lpstop=6 force_lp_range=yes`: `shc` at 20 and 16.5 A, `prob` at 13, 9.5 and 6 A).
The value 20 A is the top of the planner's own stage-1 range (10-20 A): stage 1 starts from random
orientations.

| run | stage 1 (20 A) | stage 2 (16.5 A) | stage 3 (13 A, `prob`) | stage 4 (9.5 A, `prob`) | stage 5 (6 A) | what happened |
| --- | --- | --- | --- | --- | --- | --- |
| `dn_S_c3_1` | 0.444 | 0.440 | 0.433 | 0.469 | 0.466-0.472 over 5 iterations | stopped |
| `dn_S_c3_2` | 0.440 | 0.432 | 0.392 at its 5th iteration | - | - | stopped when ruling 11 made S a stress case |

On S the probabilistic stages at 13 and 9.5 A move the labels little, from 0.43 to 0.47. The first
probabilistic iterations again redraw the labels: 0.44 to 0.39-0.41, as in run 2. A run below 0.60
after stage 4 is stopped from now on: in every run, stage 5 has added at most 0.10, so such a run
cannot reach separated.

**C4, the per-stage tail**, would be tested with the existing key `prob_beta_end=1` at entry
(every stage then carries it). **C2, the bounded draw**, was written:
- for several states, the closest reference's state, and the direction drawn within that state's
  block among its angular window;
- the draw keyed, with the untempered weights exp(-(d - dmin)) of the rule before 811934ca5, and the
  window from the iteration's model.

It changes what ten existing checks of the probability-table suite assert about the joint draw:
temperature limits, ties, and the state of the extras. So it was moved to an experiment patch,
`R/scratch/keep/variants/c2_bounded_draw.patch`, built from an export
(`R/scratch/keep/build_c2`), and the tree was restored. Neither was run on S before ruling 11.

**Ruling 11: B judges Phase 2.** The S runs still queued were stopped. B's 2D solution, P2B, is
made with HEAD (`solve2D ncls=30 mskdiam=300 nthr=16` on a copy of PZB), and one de novo run of HEAD
unchanged on B comes before any candidate.

**B's 2D solution.** P2B was made with HEAD (`solve2D ncls=30 mskdiam=300 nthr=16` on a copy of
PZB): 8 min 54 s, 0.64 GB, in `R/scratch/keep/projects/P2B`.

**HEAD de novo on B fails, at 0.699** (`dn_B_head_1`, the base build unchanged; 7 h 03 min, load
about 16). Its class-FRC ladder is 22.6, 21.3, 20.2 and 19.2 A, then 6 A. Stage 4's 19.2 A was
promoted to 6 A by the FSC=0.5 rule at the stage boundary.

| stage | 1 `shc` | 2 `shc` | 3 `prob` | 4 `prob` | 5 `prob`, with the tail and the final assignment |
| --- | --- | --- | --- | --- | --- |
| accuracy | 0.555 | 0.557 | 0.639 | 0.680 | 0.699 (ARI 0.385; recall 0.606 / 0.649 / 0.841) |

What happens:
- *Stages 1-2.* On B the hill climbing already reaches 0.555. On S no stage got past 0.47.
- *Stage 3.* The first probabilistic iterations redraw part of the labels (0.557 to 0.508). The
  stage then isolates the largest assembly: its state holds 6,613 of its 8,000 particles with 196
  medium ones by iteration 43.
- *Stages 4-5.* The smallest and the medium assemblies stay mixed over the other two states
  (3,658 / 5,162 and 4,341 / 2,642 at iteration 43), and the accuracy creeps from 0.68 to 0.70
  through 28 iterations at 6 A. The draw keeps changing about 15 % of the labels per iteration (state
  overlap 0.85) without progress. The state posterior spreads over 1.35 effective states at the
  search's factor 0.3 and 1.10 at 1, with 88 % of the draws at the state mode. That is a stable
  mixture which the hot draw only churns.

The verdict is failed: one true state is recovered, and the accuracy is below 0.75. So HEAD does
not separate B, and the candidates are judged on B in their order: one run each, and a second run
for a candidate whose first beats this run at some stage or separates. Each run is stopped by the
rule if it is below 0.60 after stage 4 (`stop_rule_watch.sh`).

**The candidates on B.** Each candidate ran once on its own copy of P2B. None beat HEAD's run at
any stage after stage 2, so none had a second run.

| run | stage 1 | stage 2 | stage 3 | stage 4 | what happened |
| --- | --- | --- | --- | --- | --- |
| HEAD (`dn_B_head_1`, above) | 0.555 | 0.557 | 0.639 | 0.680 | ran to the end: 0.699 |
| C1 (`dn_B_c1_1`, the variant build from `c1_shc_smpl.patch`) | 0.545 | 0.561 | 0.587 | 0.605 | stopped after stage 4, 1 h 17 min |
| C3, the forced ramp (`dn_B_c3_1`, 20, 16.5, 13, 9.5 and 6 A) | 0.537 | 0.537 | 0.607 | 0.645 | stopped after stage 4, 1 h 52 min |
| C4, `prob_beta_end=1` (`dn_B_c4_1`) | 0.554 | 0.555 | 0.506 at its 3rd iteration | - | stopped by ruling 12 |

C1 and C3 ended stage 4 above 0.60, so the stopping rule did not apply. They were stopped by hand
because they were below HEAD at stages 3 and 4, and stage 5 had added at most 0.10 in every run,
so neither could reach 0.90.
- *C1.* The hill climbing in stages 3-4 does not isolate the largest assembly, as `prob` does in
  HEAD's stage 3. After stage 4 its recall is 0.572, 0.457 and 0.785.
- *C3.* The ramp gives the same picture as HEAD, one stage later: the largest assembly isolated,
  the smallest and the medium mixed over two states (recall 0.620, 0.595, 0.720).
- *C2* was not run on B: ruling 12 came before its turn.

**Ruling 12: SNR 0.05 judges.** Afan's SNR 0.01 benchmark becomes a stress case: its findings
(S and B above) are reported, not used to decide. C4's B run was in stage 3, not past it, so it was
stopped, and C2 was not started on SNR 0.01 data. The judge is F05, a new simulation by Afan's
recipe at SNR 0.05. First HEAD unchanged de novo on F05, twice, then HEAD's true-pose label
recovery on F05 once.

**F05, the SNR 0.05 benchmark** (`R/scratch/keep/F05`; the command, the seed and the tables in
`make_f05.sh` and `make_f05.out`). It was made with the base build:

`simulate_particles vol1..3=<class1_8c9c, class2_8c96, class3_8c8x> nstates=3 state_pops=0.3,0.6,0.1
nptcls=80000 smpd=3 snr=0.05 ctf=yes kv=300 cs=2.7 fraca=0.07 defocus=1.5 dferr=0.5 astigerr=0.05
sherr=10 pgrp=c1 mskdiam=300`, with `SIMPLE_SEED=20261109` for this command only.

- *The particles.* 128-pixel boxes at 3 A, 24,000, 48,000 and 8,000 per state. Orientations are
  random; shifts lie within 10 pixels; defocus is 1-2 um with at most 0.05 um astigmatism.
- *The truth.* `truth_f05.txt` holds every particle's pose, shift, defocus and state.
  `f05_tables.py` turns it into `oritab_true.txt`, `oritab_zero.txt` (the poses and shifts zeroed)
  and `labels_f05.txt`.
- *The projects.* PT05 has the true poses and labels and PZ05 the zeroed poses with one state.
  PZ05's 2D solution, P2_05, comes from `solve2D ncls=30 mskdiam=300 nthr=16` on a copy.

The deviations from Afan's recipe:
1. The input maps are the README's SIMPLE reconstructions from the SNR 0.01 particles, not the
   PDB-derived maps of his simulation.
2. The noise is set per image: each image's noise SD comes from its own signal variance, so every
   image has SNR 0.05. Afan used one SD for all images, so his SNR is a global mean and varies
   per state with the assembly's mass.
3. The projection is the classic one; `wave_model` is not used (the README names none).
4. The astigmatism is `astigerr=0.05` um, for the README's U-V spread of 500 A.
5. No B-factor is applied, as in Afan's README.

The truth checks out:
- *The maps.* A reconstruction from PT05's true poses and labels (`reconstruct3D nstates=3 pgrp=c1
  objfun=cc`, 38 s) matches its own input map with correlation 0.987, 0.994 and 0.984. The six
  cross terms are 0.31-0.57.
- *The records.* PT05's 80,000 records equal the truth tables in every key: defocus, angles,
  shifts, even/odd half, stack index and state. PZ05's equal the zeroed table.

The label-recovery start on F05 is a balanced random labelling (`draw_balanced_labels.py`, seed
2026110901, accuracy 0.335) with its maps from `reconstruct3D nstates=3 objfun=cc`
(`R/scratch/keep/F05/lr_start_f05`).

**P2_05, F05's 2D solution:** 22 min 47 s. Its class-FRC resolution is 6.1 A for half the
particles, against 19-23 A on the SNR 0.01 sets.

**HEAD de novo on F05, run 1: partial, 0.831** (`dn_F05_head_1`, the base build unchanged; 6 h 51
min, load about 16). The class-FRC ladder is 20.2, 8.9, 6.7, 6.1 and 6 A.

| stage | 1 `shc` (20 iterations) | 2 `shc` (4) | 3 `prob` (17) | 4 `prob` (15) | 5 `prob` (10), then the final assignment |
| --- | --- | --- | --- | --- | --- |
| accuracy | 0.469 | 0.465 | 0.547 | 0.838 | 0.836; final assignment 0.831 (ARI 0.712; recall 0.650 / 0.952 / 0.652) |

What happens:
- *Stages 1-3 never see the smallest assembly.* These stages sample greedily: each 2D class gives
  its best particles by 2D score, the same 30,000 in every iteration. That sample holds 11,470 of
  the largest assembly, 18,529 of the medium one and 1 of the smallest. Its 8,000 particles keep
  their random start labels through stage 3. The rule is HEAD's design: the controller switches
  multi-state runs to stochastic sampling at stage 4 (`STOCH_SAMPL_STAGE_INDEP`).
- *Stage 4 samples stochastically and the labels move.* Accuracy rises from 0.547 to 0.838 in six
  iterations, then stays flat through stage 5 and the annealing tail. The state posterior is
  sharp throughout: 1.00 effective states, 99.9 % of the draws at the mode.
- *The failure is one mixed state.* Found states 1 and 2 are nearly pure (the medium assembly's
  recall is 0.952). Found state 3 holds 5,219 of the smallest assembly, 8,397 of the largest and
  1,766 of the medium.
- *The strays are mostly particles the greedy stages never searched.* Of the largest assembly's
  particles outside the greedy sample, 40 % sit in found state 3, against 22 % of those inside it.
  For the medium assembly the shares are 6.8 % and 1.2 %. Quality does not explain this. Among the
  largest assembly's particles outside the sample, those in found state 3 and those in found
  state 1 have the same 2D score (median 0.386 and 0.388), and the sample's particles score
  barely higher (0.40). The strays are the particles that arrived late.

**HEAD de novo on F05, run 2: partial, 0.824** (`dn_F05_head_2`; 7 h 23 min). The same picture
as run 1: the smallest assembly keeps its start labels until stage 4, and one found state mixes it
with part of the two others (5,248 of the smallest, 7,127 of the largest, 4,142 of the medium).

| stage | 1 | 2 | 3 | 4 | 5, then the final assignment |
| --- | --- | --- | --- | --- | --- |
| accuracy | 0.515 | 0.511 | 0.558 | 0.844 | 0.837; final assignment 0.824 (ARI 0.666; recall 0.703 / 0.913 / 0.656) |

HEAD fails F05 in both runs, so the candidates are judged on F05 with these two runs as the bar
(0.838 and 0.844 after stage 4, 0.831 and 0.824 at the end). Ruling 13 puts C5 first.

**HEAD's label recovery on F05 separates, exactly** (`lr_F05_head`; the true poses, a balanced
random start, seed 2026110901, at 0.335; `update_frac=0.375`, `prob_beta_end=1`; 58 min). Accuracy
per iteration: 0.381, 0.602, 0.807, 0.890, 0.973, then 1.000 from the sixth iteration on. With the
true poses, today's state assignment recovers the three states completely. The de novo failure
therefore lies upstream of it, in what the stages hand it, which is where the greedy sampling acts.

**C5: stochastic sampling from stage 1 for several states.** Variant patch
`R/scratch/keep/variants/c5_stoch_sampl_stage1.patch`, built from an export (`build_c5`); the tree
stays at HEAD until a candidate is kept.
- *The change.* The stage controls set `greedy_sampling=no` in every stage of a multi-state run;
  `nstates=1` keeps its rules.
- *What goes with it.* The multi-state switch stage (`STOCH_SAMPL_STAGE_INDEP`) and its one-caller
  function `solve3D_stoch_sampl_stage` are deleted, with the function's interface in
  `simple_solve3D_utils.f90`.
- *What stays greedy.* The start volumes' initial sample. It marks its particles as updated, so
  stage 1's stochastic draw takes the others first.

In stage 1 every iteration now samples a different set, with about 1,850 of the smallest
assembly's particles (HEAD: 1). The class quotas give the smallest assembly about 6 % of the sample
for its 10 % of the particles.

| run | stage 1 | stage 2 | stage 3 | stage 4 | what happened |
| --- | --- | --- | --- | --- | --- |
| HEAD, runs 1 and 2 (the bar) | 0.469, 0.515 | 0.465, 0.511 | 0.547, 0.558 | 0.838, 0.844 | partial, 0.831 and 0.824 at the end |
| C5, run 1 (`dn_F05_c5_1`) | 0.775 | 0.744 | 0.786 | 0.774 (recall 0.648 / 0.964 / 0.011) | stopped by hand after stage 4; not judged (below) |
| C5, run 2 (`dn_F05_c5_2`) | 0.730 | 0.766 | 0.759 | 0.905 | partial: stage 5 0.918, final assignment 0.920 (ARI 0.821; recall 0.918 / 0.963 / 0.673); 7 h 22 min |
| C5, run 3 (`dn_F05_c5_3`) | 0.825 | 0.771 | 0.874 | 0.887 | **separated**: stage 5 0.902, final assignment 0.9025 (ARI 0.770; recall 0.926 / 0.905 / 0.819); 8 h 23 min |
| C5, run 4 (`dn_F05_c5_4`) | 0.783 | 0.761 | 0.819 | 0.829 | partial: stage 5 0.892, final assignment 0.9007 (ARI 0.799; recall 0.848 / 0.971 / 0.639); 8 h 49 min, the gates beside it |

The early stages change a lot:
- *Stages 1-2.* Stage 1 ends at 0.775. Stage 2 runs its 20 iterations (HEAD's ended after 3-4 on
  the scan rule) and ends with every assembly near 0.7 recall.
- *Stages 3-4.* The largest assembly splits over two states, and the smallest loses its state.
  Its particles end in the state of the largest assembly's bigger part (4,853) and in the medium
  assembly's (3,062).

*Run 1 was stopped by mistake.* It was stopped after stage 4, below HEAD, on the reasoning that
stage 5 adds at most 0.10 and so cannot reach 0.90. That bound came from the SNR 0.01 runs, and
C5's second run disproved it on F05: its stage 4 alone gained 0.146. Run 1 was on the same path,
only more slowly. The largest assembly's split-off part shrank through stage 4 (13,055 to 9,138),
and the accuracy had been rising from its 66th iteration (0.758 to 0.774). Its outcome is
therefore unknown. Run 1 is not judged, and C5's two of three are judged on runs 2, 3 and, if
needed, 4. From then on no F05 run is stopped by hand above 0.60 after stage 4; only the formal
rule applies.

*Run 2 nearly separates.*
- *Stages 2-3.* The largest assembly is split over two states and the smallest has none (recall
  0.000 after stage 3), as in run 1.
- *Stage 4.* The split-off part shrinks (8,695 to 3,955 particles by the 10th iteration), and the
  freed state takes the smallest assembly. Its recall rises by about 0.1 per iteration (0.004 to
  0.607), and the stage gains 0.146.
- *Stage 5.* The recall rises more slowly (+0.06 over 10 iterations) and ends at 0.673, below the
  0.75 a separated verdict asks. 2,617 of the smallest assembly's particles stay in the medium
  assembly's state.
- *Against HEAD.* Run 2 beats HEAD at every stage. As C5's first judged run, it earns a second
  judged run, run 3. If run 3 also fails, C5 is out and C4 and C2 follow (ruling 13). If it
  separates, a run 4 decides, in place of the unjudged run 1.

*Run 3 separates F05 de novo, the first run to do so.* Every recall is above 0.75 from stage 3
on, and the accuracy crosses 0.90 in stage 5. The smallest assembly's state still holds 4,441 of
the medium assembly and 1,765 of the largest. C5 now has one separated run of two judged, so a
fourth run decides.

**C5 in the tree, with its gates** (made while run 4 runs; reverted from HEAD if run 4 fails). The
file table gained `simple_solve3D_utils.f90` (the deleted function's interface), the tester and,
for Phase 4, the sampling policy. The tree's controller and interface equal the variant that ran.
The tester's test of emitted stage commands now checks the sampling of stages 1-5:
- one state: `greedy_sampling=yes`, with `frac_best` 1.0 and 0.5 from stage 5;
- three states: `greedy_sampling=no` and `frac_best` 1.0.

For three states, HEAD's controller would fail those checks in stages 1-3, which sample greedily.
This was read from the code, not run.

| gate | result |
| --- | --- |
| build | no warnings in the changed files |
| `unit_reconstruction suite=solve3D_cavgs_stage_rules` | 81 of 81 (61 before) |
| fast gate | 15 of 15 (`logs/c5_fastgate.log`) |
| single-state identity: `simulated_workflow_6vxx` at one thread, base export against the tree | 186 of 186 declared outputs identical (`identity_c5/identity.log`) |
| `labels` | pass, 9 min 12 s beside a de novo run; accuracy 0.999, every true state resolved |
| `shc` | pass, 58 min 52 s beside a de novo run; accuracy 0.994, every true state resolved |
| `hard_shc` | pass (its labels reported, not gated), 37 min 27 s; accuracy 0.800, the conformational pair merged (C1: 0.527) |

In each route the shuffled-label control does not resolve every state, as intended (as under C1).

*C5 alone: one of three judged runs separates.* Run 4 ends at 0.9007, but the smallest assembly's
recall is 0.639: 2,892 of its particles sit in the medium assembly's state, and its own state holds
3,633 of the largest. Over runs 2-4 every final accuracy is at or above 0.90 (0.920, 0.9025,
0.9007), against HEAD's 0.831 and 0.824. What fails is the smallest assembly's recall: 0.673,
0.819 and 0.639.

C5 is kept. Under the plan's rule it improves both the de novo verdicts (one run separated against
none) and the accuracy at every stage, and it breaks no gate. The success criterion is now sought
by combining the ruled candidates with it, C4 first, then C2. The bar is C5's runs: a second run
follows when a first separates or ends above C5's median final, 0.9025.

**C5 and C4** (`dn_F05_c54_1`: the C5 build with `prob_beta_end=1`, so every multi-state `prob`
stage anneals to 1). Stages 0.800, 0.754, 0.876, 0.886 and 0.896. The final assignment gives
**partial, 0.8973** (ARI 0.761; recall 0.923 / 0.904 / 0.782; 8 h 12 min). Every recall is above
0.75 from stage 3 on, and the accuracy is 0.003 short of 0.90. It is the same picture as C5's
run 3: the smallest assembly's state also holds 4,616 of the medium assembly. The annealing tails
add nothing visible. The run is below C5's median final and does not separate, so there is no
second run. C5 and C2 follow.

**C5 and C2, not judged** (`dn_F05_c52_1`: the C5 build with the bounded draw,
`variants/c2_bounded_draw.patch`). Stages 0.701, 0.765 and 0.867. At iteration 67, early in stage 4,
the accuracy was 0.903 (ARI 0.781), but the smallest assembly's recall was 0.348 (recall 0.919 /
0.987 / 0.348): 5,215 of its 8,000 particles had joined the largest assembly's state. The smallest
assembly had already lost its state in stage 3 (recall 0.273 at iteration 52).

**Stopped 2026-11-11 07:35 at the maintainer's request**, to sync with the current code before
going further. C5 was ported onto a321f2904 (the stage-sampling block moved, the logic was
unchanged): Release build, fast gate 15/15, `unit_reconstruction` 22/22. The single-state
identity and the multi-state routes were not rerun on the port. Evidence for the next step, from
the journal: shc brings the smallest assembly to a recall of 0.70-0.77 by stage 2 in every C5-based
run (0.725, 0.657, 0.768, 0.738, 0.765, 0.701). Stage 3's jump to 6.7 A `prob` then either keeps it
(C5 run 3 0.792, C5+C4 0.781) or destroys it (C5 run 1 0.660, then 0.011 at stage 4; run 2 0.000;
C5+C2 0.27). So the next candidate when the run resumes is C3 on the C5 code (`lpstart=20 lpstop=6
force_lp_range=yes`, `prob` at 13, 9.5 and 6 A). C1 is not favoured. The run's records stay on
dell2 under `~/agent_runs/ms_abinitio/` (`scratch/keep/runs/summary.tsv`, `state/journal_phase_*.md`,
`scratch/keep/variants/`).
