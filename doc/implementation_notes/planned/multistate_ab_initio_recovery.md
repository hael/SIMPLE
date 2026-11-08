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
  one of them. The joint state draw and the old argmin behave alike here.
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
| Stages | 8 for every state count | 5 for `nstates>1`, `lpstop` 6 A | 2026-10-16/18 |
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
`shc` profile at or above its floor, and the nightly passing.

## 5. Phases

### Phase 0: builds, data, and why it fails here

1. Build HEAD in `repo/build` (`./compile_clean.sh`). Build two anchors from `git archive`
   exports in `scratch/keep/`: **A25**, tag `v3.1-beta` (2025-05-27, the paper's code:
   `abinitio3D`, `multivol_mode=independent`, `shc_smpl` stages, bounded draw), and
   **A1015**, 1bbbdcc24 (2026-10-15, before the whole-row posterior: bounded draw and
   frontier, `shc` stages 1-2, prob 3-5, eight stages unless `nstages`). Old code under
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
  22.6, 20.2, 6 A).
- **C4, the per-stage annealing tail**: `prob_beta_end=1` on every multi-state
  probabilistic stage (an existing key).
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
condition of section 3 applies).

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
| `src/main/params/simple_parameters.f90` | C1: `shc_smpl` back in the `refine` key's documented values (no new key) | 2 |
| `src/main/params/simple_parameters_phases.f90` | C1: the mode's derived flags | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | C2: the bounded draw and the frontier for several states | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | C2: the same in the neighbourhood table | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | C2: the bounded sampler | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_tester.f90` | C2: its tests | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | C1, C3, C4: the stage modes, schedule and tails for several states | 2 |
| `src/main/commanders/simple/simple_commanders_solve3D.f90` | C3: the multi-state stage count | 2 |
| `src/main/simple_convergence.f90` | a stage rule the candidates need | 2 |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90` | path checks of the routes when the stage modes change | 2 3 |
| `doc/policies/3D/solve3D_policy.md` | section 6 and the schedule | 4 |
| `doc/policies/3D/refine3D_policy.md` | a restored mode or draw | 4 |
| `doc/policies/test_environment_policy.md` | a route's description if it changes | 3 4 |
| `doc/implementation_notes/completed/multistate_ab_initio_recovery_report.md` | the report | 4 |
| `doc/implementation_notes/completed/multistate_ab_initio_recovery.md` | the plan, moved at the end | 4 |

Other files a phase needs are added to this table, with the reason, by that phase.

## Progress

(one section per phase, written by the run)
