# Multi-state solve3D regression: diagnosis on a simulated two-state set

Planned 2026-10-21 for a Dell run of about two days (arun run `ms_regression`), base 72930a3a4.
Revised the same day after the maintainer's rulings (section "Rulings"). Nothing in this note has
been run yet. The facts in "What we know" were read from the code and the history and checked by a
second reader; the run checks again whatever it relies on.

## What was reported

Cyril (2026-10-20): two runs of de novo `solve3D` with `nstates=4` on his real data ("Sov") no
longer give what he obtained before, which was four states with a fair bit of smearing. He also
has a simulated set of 30 000 particles in two states that he cannot reproduce. The set is in
`~/Processing/simulated_vcrnapol` with the script that makes it (`generate.sh`). He does not know
which build gave his earlier result.

Two failed runs are not yet proof of a regression: on our own three-state test the de novo route
resolves every state in only 13 of 20 runs. The run therefore compares how often builds separate
the states, and *where* in the stage ladder the separation forms or is lost.

## Terms

- **State**: one of the `nstates` structures a multi-state refinement reconstructs. Each particle
  carries a state label; state 0 means the particle is excluded.
- **Stage**: `solve3D` runs a ladder of refine3D executions ("stages") at increasing resolution.
  With `nstates>1` there are five: stages 1-2 use `shc` (stochastic hill climbing, in which states
  and projection directions are searched on an equal footing), stages 3-4 use `prob` (assignment
  from a table of particle-reference scores: projection directions drawn at random, the state
  taken as the best over all states), stage 5 uses `prob_neigh` in state mode (the same over a
  neighbourhood of directions shared by the states). Each stage has a planned low-pass limit (its
  "rung") and a maximum number of iterations.
- **Distributed / shared memory**: `nparts>1` runs a master process plus worker processes per
  iteration; without `nparts` one process does everything. Cyril runs distributed (`nparts=3`).
- **res05**: the resolution at which a state's half-map FSC falls to 0.5.
- **Resolution-limit updates**: the rules that move a stage's limits during a run: stage omission
  (a stage is skipped when res05 is already finer than its rung), FSC=0.5 promotion of a stage's
  limit, the mid-stage exit (a prob stage ends once res05 passes the next rung), and the matching
  limit taken from the NU hand-off.
- **NU filtering and NU hand-off**: nonuniform low-pass filtering of the reconstructed maps during
  assembly, on from stage 5; the hand-off is the limit the filter reports back for matching. The
  matching limit never goes finer than `lpstop`, 6 Å for multi-state.
- **Accuracy**: the fraction of the particles active at the start of `solve3D` whose final state
  equals the true state, taking the better of the two ways to match found states to true states;
  state 0 counts as wrong. Putting every particle in one state scores 0.567 on this set.
- **ARI**: adjusted Rand index between found and true labels (0 for random, 1 for perfect).
- **Signature**: the stage at which accuracy rises above 0.75 or falls back, plus the stages
  omitted or cut short, read from the labels written after every stage and iteration.
- **Anchor**: an older build expected to work, measured beside the current one.

## Rulings

Of 2026-10-21, on the first version of this note:

1. **The speed changes stay.** Nothing obtained in the speed refactoring is given up: the
   multi-state schedule (`prob_neigh` in stage 5, the stage-2 scan exit, the 2 500-direction cap),
   the table and assembly changes. A diagnostic variant may undo one of them to localise a cause;
   a fix may not. A per-run time of at most 30 minutes is expected.
2. **NU filtering at stage 5 stays as it is**, whatever the evidence; it is not varied.
3. **Resolution-limit updates may change for `nstates>1`** and stay exactly as they are for
   `nstates==1`.
4. **Distributed and shared memory are equivalent.** Where they differ, the difference is a defect
   and is fixed.
5. **Stages 1-2 stay as they are**: `shc` with sticky sampling treats states and projection
   directions on an equal footing, which is the approach that worked in the past.
6. **A temperature on the state assignment in stages 3-4** is welcome. Going straight to the
   greedy state assignment in the probabilistic step may be unwise. From stage 5 onwards the
   assignment stays as it is now: projection directions drawn at random, states greedy.
7. **An emptied state is re-seeded by a balanced random relabelling**; a restart of the ladder may
   have to be considered when that happens.

Earlier rulings that hold: de novo multi-state stays, without several starts and selection
(2026-10-18); sampling and relabelling are never driven by the maps or poses (2026-10-04), so no
fix relabels particles against the current maps; state 0 means exclusion and a true state in
state 0 is a failure (2026-10-17); de novo routes are judged over many runs, never one
(2026-10-18); bugs found are fixed in the run (2026-10-17).

Readings of these rulings that the run follows (correct them before launch if they are wrong):

- Ruling 1: a time change that comes only from a multi-state resolution-limit rule (more
  iterations because a stage is no longer omitted) is allowed by ruling 3 and is reported per
  stage; a slower kernel, table or assembly is not allowed.
- Ruling 4 covers the frozen stage-1 references and the empty-state handling (long-standing items
  1 and 2): both are fixed, single-state included, since shared memory already does the right
  thing.
- Ruling 6: the temperature is a parameter (default greedy, so `refine3D` outside `solve3D` is
  unchanged) that `solve3D` sets for stages 3 and 4 only. A particle's state is drawn from a
  softmax over states of its best score in each state, divided by the temperature; projection
  directions are then drawn within the drawn state as now. The run implements it in its own phase,
  chooses the values by a small sweep and keeps it only if it does not lower the separation counts.
- Ruling 7: whenever a state is empty (or below the prob stages' minimum population) at a stage
  start, all active particles are relabelled at random with balanced populations, in both modes
  and for any number of states (today an empty state 1 stops a prob stage, three or more states
  with an empty last state stop, and shared memory relabels only at the first iteration). The run
  counts these events and whether runs recover from them; a restart rule is proposed from that
  evidence, not implemented.

## What we know before any run

### The simulated set

`generate.sh` simulates 17 000 particles of EMD-51275 (labelled state 1) and 13 000 of EMD-51273
(state 2), box 180 at 2.055 Å, SNR 0.1, CTF with 2 ± 1 µm defocus, random poses, shifts up to
3 pixels, B-factor 100 ± 5. It stacks and shuffles them, imports them with `deftab=`, runs
`solve2D ncls=10`, then `solve3D` single-state and with `nstates=2`, each with `nparts=3 nthr=16`.

Checks made on 2026-10-21:

- The two maps are not in a common frame (real-space correlation 0.12 as given). After a rotation
  search they agree closely: FSC 0.94 at 120 Å, 0.83 at 30-50 Å, 0.48 at 13 Å, against 0.99 for
  identical maps at the lowest shell. About 30 % of the density differs. This is a moderate
  difference, visible at low resolution, so 30 000 particles at SNR 0.1 should separate reliably
  with good poses. Found maps are compared with the truth after docking, in both hands.
- `import_particles deftab=` takes only CTF parameters, `state` and `eo` from the file, so neither
  poses nor labels leak into the project (the script sets every state to 1 first). The import
  stores the stack by absolute path, and project row i is stack image i.
- `stackops nran=30000` draws a permutation without repeats and carries each particle's
  parameters with it, so `out.txt` line i describes stack image i.
- `labels.txt` is column 17 of `out.txt`. Which column holds `state=` depends on how SIMPLE writes
  the parameter line, so the run takes the labels by key and compares them with column 17.
- The set must be generated without `SIMPLE_SEED`, as the script does. Every SIMPLE process
  restarts the seed sequence, so with a fixed seed the two `make_oris` and the two
  `simulate_particles` calls would draw the same poses, shifts and noise for the first 13 000
  particles of each state.

### What changed on the multi-state path

The commits since 2026-09-25 that change code a `solve3D nstates=2` run from a 2D project
executes. Stream, GUI and test-only commits are left out.

| date | commit | change | evidence when it landed |
| --- | --- | --- | --- |
| 10-05 | 7c24f2b45, 311b8817a | class-average balancing by default (with 10 classes the same as class balancing); distributed shc sample drawn once on the master | not measured on multi-state |
| 10-13 | 2307ea96f | global sigma2 only | single-state and beta-galactosidase runs |
| 10-13 | cfbc56bf3 | passes that search nothing no longer relabel states; a collapsed project now ends collapsed | not measured on multi-state |
| 10-15 | e2673a5d7 | `prob_fast` retired: the particle-frame shift grid and quadratic step in every table, so also in the multi-state stages 3-4 | its report says multi-state was not run |
| 10-16 | 57a447151 | the multi-state schedule: stage 5 `prob` → `prob_neigh` in state mode at 2 500 directions; stage 2 also ends when more than 95 % of the search is scanned; state assignment rewritten (same rule) | paired 20-run test on the three-state set (2026-10-17): 16 → 11 of 20 runs separate every state; accepted |
| 10-18 | cd13708eb | stage omission and FSC=0.5 promotion driven by the best state's res05 | not run on multi-state |
| 10-18 | 8c8d8b427 | NU filtering from stage 5, also for multi-state; the mid-stage exit set to the next rung | 13 of 20 on the profile of 2026-10-19 (code of 10-18) |
| 10-19 | 0c9542e50, 24d50aabf | stage progression changes; the "no longer gaining" guard of the mid-stage exit made able to fire | unit test only |
| 10-19 | af13bc828 | 2D class sums from cropped particles | equal to rounding, not bitwise; the class FRCs set `solve3D`'s ladder |
| 10-15 to 10-20 | 04e520e78, 43de9b340, aea5b1846, af13bc828, f078cf683, dafd4d046, 72930a3a4 | speed and structure changes: parallel table reads, dynamic loops, reconstruction from the search's cropped records, the 3D table by reference, OpenMP private-variable fixes, one project owner on the master | byte-identical at one thread, including a three-state refine3D case; thread counts above one not covered (a race in the table reads, present from 04e520e78, was fixed in 43de9b340) |

Measured on multi-state at all: the three-state test (3 000 simulated particles of three
beta-galactosidase states, shared memory, no 2D start) from 2026-10-17 to 2026-10-19 only. Cyril's
path differs from it in four ways: it is distributed, it starts from a `solve2D` project (so the
class-average balancing is active), it has two states, and it is ten times larger.

### Candidate causes, and what tells them apart

1. **Frozen stage-1 references in distributed runs** (long-standing item 1). Old, but it matters
   more now that stage 2 is shorter: in a distributed run the states can only start to form in
   stage 2. *Signature:* shared memory separates where distributed does not; accuracy after
   stage 2 low in distributed runs.
2. **Stage 2 shortened** (10-16, a speed item that stays). *Signature:* fewer stage-2 iterations
   and lower accuracy after stage 2 than the anchor.
3. **Resolution-limit updates cutting stages 3-4** (10-18, 10-19). At the start both states hold
   mixtures of nearly the same structure, so the best state can reach FSC=0.5 early and the prob
   stages, where the states split, can be skipped or shortened. *Signature:* the log lines
   `SOLVE3D OMITTING STAGE`, `STAGE RESOLUTION EXIT` and `promoted by FSC=0.5`; no
   `states_stage3.txt` or `states_stage4.txt` (an omitted stage writes none); short stages.
4. **Greedy state assignment in stages 3-4 with the prob table of 10-15.** The state is the argmin
   over all references, and only the best references across all states get a refined shift, so
   the leading state is favoured; a state that falls behind early cannot come back. *Signature:*
   stages 3-4 run their full length but gain less accuracy than in the anchor, or one state's
   population drains during them.
5. **Stage 5** (`prob_neigh` in state mode since 10-16; NU filtering since 10-18, which stays) and
   the matching limit taken from the finest NU hand-off over the populated states (a
   resolution-limit update). *Signature:* separated after stage 4, lost in stage 5.
6. **A defect in the speed and structure changes**, for example a thread race or a state indexed
   wrongly, which one-thread identity cannot show. *Signature:* a crash, NaN, a state whose map
   never changes, or a bisection that lands on one of these commits.
7. **No regression.** The method decides the states early by chance (diagnosis of 2026-10-18).
   *Signature:* the anchor fails as often, and in the same way, as the current build.

### Long-standing behaviour found while reading

Older than the window, present in the anchor too. Both are defects by ruling 4 and are fixed in
Phase 3.

1. **In distributed runs the stage-1 references never change.** The master keeps one `params` for
   the whole stage and builds the references from `params%vols_even/odd`, set once from the
   start volumes; assembly writes the new maps under other names and updates only the command
   line (`simple_refine3D_strategy.f90`, `simple_matcher_refvol_utils.f90`,
   `simple_parameters_phases.f90`). From stage 2 on the start volume carries the name assembly
   rewrites, so only stage 1 is frozen. Shared memory rebuilds `params` from the command line every
   iteration and is not affected. A note of 2026-08-23 (`drop_legacy_box_division.md`, section 8)
   saw the reference power flat over all 20 stage-1 iterations. In a distributed multi-state run
   this means twenty `shc` iterations against two noise volumes; the three-state test (shared
   memory) never saw it. Single-state distributed runs have the same frozen stage 1.
2. **Empty states.** With two states, a distributed stage start relabels every particle at random
   (balanced, `gen_labelling 'uniform'`) when state 2 is empty; an empty state 1 stops a prob stage
   with an error; with three or more states an empty last state is an error; shared memory
   relabels only at the first iteration.
3. **The NU hand-off of a state without an FSC=0.143 crossing** is the finest limit of the bank (a
   9-particle state handed off 2.6 Å in the 2D speed run of 2026-10-19); through the minimum over
   states it sets the matching limit of every state. A multi-state resolution-limit rule, so open
   to change for `nstates>1` (ruling 3).

## Before launch (no Dell time)

1. **Ask Cyril** (text at the end) for his command lines, a few log lines, his core count and a
   check of `labels.txt`. His log lines may already show cause 3.
2. **Copy the data to the Dell**:
   `rsync -a ~/Processing/simulated_vcrnapol/ dell:/data/simulated_vcrnapol/` (or any read-only
   place; the run reads it from `DATA_DIR` in `phases.md`).
3. **Commit this note** (arun copies tracked files only).

## How runs are made and judged

**Builds.** The run's tree (HEAD, 72930a3a4) builds in `repo/build`. Other commits are exported
with `git archive <commit> | tar -x -C scratch/src_<commit>` and built in `scratch/build_<commit>`;
the repository's HEAD never moves. The anchor is 1bbbdcc24 (2026-10-15), which predates the
`prob_fast` retirement, the schedule change and every later change in the table above. Never build
ab2b5f5cd (it does not build; use 72930a3a4). A build from 04e520e78 up to 8c8d8b427 gets
43de9b340's fix of the table-read race applied (its two production files), so that the race does
not stand in for the commit being tested. Builds no longer needed are deleted (the sweep keeps
build trees).

**Thread layout.** One run at a time on the whole machine: distributed `nparts=3 nthr=8`, shared
memory `nthr=24`. Timings are recorded per stage; they are judged only in Phase 5 (ruling 1).

**Projects.** P2: import and `solve2D ncls=10 mskdiam=270` made by the anchor, read by every
newer build (project files of older builds must stay readable). P2H: the same made by HEAD, for
the end-to-end arm. Each `solve3D` run works in its own project directory, launched from it (before
b57538b7b the sigma2 store was looked up in the current directory), sharing the one stack.

**The command** (Cyril's, Dell threads, labels after every stage and iteration):

```
simple_exec prg=solve3D pgrp=c1 mskdiam=270 nparts=3 nthr=8 nstates=2 print_states=yes projfile=<its copy of P2>
```

No `SIMPLE_SEED`: in a distributed run every worker would restart the same random sequence every
iteration, which Cyril's runs do not do. Runs of an arm are independent repetitions.

**Recorded per run:** the `SIMPLE Git Commit` line; wall time per stage; `states_stage<N>.txt`
after every stage, `states_iter<NNN>.txt` after every iteration (both modes at HEAD; shared memory
only at the anchor), the final labels; the log lines `ENTERING REFINE3D STAGE`, `STAGE <n> WITH`,
`promoted by FSC=0.5`, `SOLVE3D OMITTING STAGE`, `STAGE RESOLUTION EXIT`,
`STATE <s> JOINT DISTRIBUTION OVERLAP ... POPULATION`, `STATE OVERLAP`, the per-state resolution
lines, the NU hand-off per state, any message about empty states or relabelling; the final state
maps. Labels, logs and maps of every run go to `scratch/keep/runs/<arm>_<n>/`.

**Metrics:** accuracy, ARI and the confusion matrix (true state against found state 0, 1, 2) of the
final labels and of every stage's and iteration's labels; populations against 17 000 / 13 000;
the number of particles `solve2D` left inactive (excluded from accuracy, reported); each final map
docked to each true map in both hands (`dock_volpair`, the mirror by `volops`): correlation and
FSC=0.5 resolution against the truth; stages omitted or cut short; iterations per stage; relabel
events; the signature.

**Verdict per run:** *separated* at accuracy ≥ 0.90; *partial* from 0.75 to 0.90 (the smearing
Cyril describes); *failed* below 0.75. If the known-reference control of Phase 0 stays below 0.95,
the threshold becomes its accuracy minus 0.05. The threshold is fixed in the journal at the end of
Phase 0 and never moved afterwards.

**Verdict per build.** Eight runs for the main arms, four for probes. Counts of this size separate
large differences only (8 of 8 against 2 of 8, not 6 against 4), so the run decides on the count
and the signature together: the signature (where the states separate or are lost, which stages are
omitted) is close to deterministic for a build. A build is *good* when most of its runs separate,
*bad* when most fail with the same signature; a mixed build is judged by whether its failures share
HEAD's signature.

**Tools.** SIMPLE programs do the SIMPLE work (`simulate_particles`, `import_particles`,
`reconstruct3D`, `refine3D`, `dock_volpair`, `volops`, `fsc`, `print_project_field`). Accuracy and
the tables are computed from text dumps by a small script kept in `scratch/` that only reads
outputs.

**Budget**, at 30 minutes per run (Phase 0 measures it and the journal rescales):

| phase | content | estimate |
| --- | --- | --- |
| 0 | builds, data, truth and known-reference controls, P2 and P2H, first H and A runs | 3 h |
| 1 | H and A to eight runs each, four shared-memory runs, two end-to-end runs | 11 h |
| 2 | diagnostic variants (six runs each) and at most four bisection probes (four runs each) | 12 h |
| 3 | the mode fixes and the regression fix, with their runs and checks | 8 h |
| 4 | the temperature: code, a sweep of three values (four runs each), eight runs at the chosen one | 9 h |
| 5 | Cyril's path as a test, `test=nightly`, the three-state profile, timings, documents, report | 8 h |

If the time runs short, drop in this order: bisection probes beyond two, the second end-to-end run,
shared-memory runs beyond two, the sweep down to two values. A probe aimed at a signature that shows
by stage N may be stopped once `states_stage<N>.txt` is written.

## Phases

### Phase 0: builds, the data, controls, timing

1. Build HEAD and the anchor.
2. Generate the set with HEAD as `generate.sh` does (no `SIMPLE_SEED`), with the Dell's thread
   counts, into `scratch/keep/data`. Check that the first lines of `oris_1.txt` and `oris_2.txt`
   differ. Take the labels by key from `out.txt` and compare them with column 17; if they differ,
   note it for Cyril, whose own evaluation then read the wrong column.
3. Truth control: import `out.txt` with `oritab=` (true poses and states) into a separate truth
   project and reconstruct both states (`reconstruct3D nstates=2`). Each map against its input map
   (FSC, same frame) proves the labels, the simulation and the shift convention; the FSC between
   the two truth maps after docking is the difference a separation must find. Repeat with the
   labels shuffled: the maps a failed separation gives.
4. Known-reference control: `refine3D nstates=2 refine=prob` for two iterations on the truth
   project, with the two truth reconstructions as references, distributed. The accuracy it keeps
   is the ceiling for this data and sets the threshold.
5. P2 with the anchor and P2H with HEAD; count the particles `solve2D` leaves inactive.
6. One HEAD and one anchor two-state run on P2, timed. They are the first runs of arms H and A.
   Write the rescaled plan into the journal.

Exit: the truth maps match their inputs, the labels are confirmed, the threshold is fixed, P2 and
P2H exist, the time per run is known and the rescaled plan is written.

### Phase 1: reproduce, the anchor, the modes

| arm | build | what | runs |
| --- | --- | --- | --- |
| H | HEAD | `solve3D nstates=2` on P2, distributed | 8 |
| A | anchor | the same | 8 |
| HS | HEAD | the same in shared memory | 4 |
| E | HEAD | the same on P2H (Cyril's whole pipeline) | 2 |

Write the side-by-side table of H, A and HS: accuracy after every stage, iterations per stage,
omitted stages and early exits, populations, relabel events, final maps against the truth,
signatures.

Decision, written into the journal with the numbers before Phase 2 starts:

- **A good, H bad**: the regression is confirmed. Phase 2 starts from H's signature.
- **Both bad with the same signature**: measure an older anchor (b818cd7cc of 2026-10-05, with its
  own 2D project and two H runs on that project too). If it also fails, there is no regression in
  the window on this set: Phase 2 is skipped and Phases 3 and 4 go ahead (the mode fixes and the
  temperature address the mechanism either way).
- **H good**: the failure does not reproduce on the Dell. Run `generate.sh` end to end once with
  scaled thread counts, describe what differs from Cyril's setup, and go on with Phases 3 to 5.
- **HS good, H bad** points at the distributed path (cause 1 or 6); **E bad, H good** points at
  the 2D side (compare the two ladders and the class FRCs).
- **Mixed counts**: decide by signature; if the signatures do not separate the builds, Phase 2 runs
  the diagnostic variants only, without bisection.

Exit: the decision and the numbers behind it are in the journal.

### Phase 2: localise

Step 1, diagnostic variants on HEAD, chosen by H's signature, six runs each. They are scratch-only
builds or command lines that localise; none of them is a fix (rulings 1 and 2):

- **Stage-1 references refreshed in distributed runs**: a scratch patch that updates
  `params%vols`, `vols_even` and `vols_odd` of populated states after assembly. For cause 1.
- **No resolution-limit updates for several states**: a scratch patch that turns stage omission,
  FSC=0.5 promotion and the mid-stage exit off when `nstates>1`, only if H's logs show them. For
  cause 3.
- **The schedule before 10-16**: a scratch patch re-derived from 57a447151's hunks (stage 5 back
  to `prob`, no stage-2 scan exit). For causes 2 and 5; diagnostic only.
- **The prob table before 10-15**: the commit pair e2673a5d7 against its parent, if stages 3-4 run
  their full length and still gain less. For cause 4.

Step 2, only if no variant restores the anchor's behaviour: at most four bisection probes over the
first-parent commits between the anchor and HEAD that touch the path, chosen by the signature (a
stage-5 signature, for example, puts the first probe between 57a447151 and 8c8d8b427). In order:
b57538b7b, e2673a5d7, 65bee3073, 57a447151, 13795c1d9, 04e520e78, 189e0a355, fb4e653ac, cd13708eb,
8c8d8b427 (these five with 43de9b340's fix), 43de9b340, 9531dee62, f4d6accff, aea5b1846, af13bc828,
0c9542e50, f078cf683, dafd4d046, 24d50aabf, 497cc974b, 68278fa9a, 72930a3a4.

For the cause found, show where it acts (stage and iteration, from the labels) and how (log lines,
populations per iteration).

Exit: the change or changes responsible are named and the mechanism shown, or "no single cause"
with the evidence.

### Phase 3: the mode fixes and the regression fix

Fixed in this order, each verified before the next:

1. **Distributed stage-1 references** (long-standing item 1, ruling 4): after assembly the master
   updates the reference names of populated states as shared memory does, so stage 1 refines
   against its own reconstructions. Verified by a scratch-only print of the files the master reads
   at each iteration, then four H runs and two single-state distributed runs on P2 (map against the
   truth, before and after).
2. **Empty states** (long-standing item 2, ruling 7): at every stage start, in both modes and for
   any number of states, a state that is empty or below the prob stages' minimum population leads
   to a balanced random relabelling of all active particles, logged. A unit check covers an empty
   first, middle and last state in both modes.
3. **The regression's cause**, as Phase 2 named it, within the rulings: a resolution-limit rule
   may change for `nstates>1` only; NU filtering at stage 5 and the speed changes stay. A defect
   (code that does not do what its commit intended, a state indexed wrongly, a thread error) is
   fixed at its cause. If the only fix found would undo a speed change or touch NU filtering, the
   run stops and describes it.

Verification of the phase: eight H runs and four HS runs on the fixed tree, against Phase 1's H
and A; distributed and shared memory now agree within the counts and in their signatures. For
every change gated on `nstates>1`, single-state stays byte-identical at one thread
(`simple_test_exec test=identity` on the single-state workflow case). Changes 1 and 2 are expected
to change single-state distributed results; they are judged by `test=nightly` in Phase 5.

Exit: the fixes are in the tree with their evidence, or the run has stopped with a written design
question.

### Phase 4: a temperature on the state assignment in stages 3-4

1. Read how the prob table's scores are scaled (`simple_eul_prob_tab.f90`: the scores
   `ref_assign` compares, and how they are normalised) and write it into the journal; the
   temperature is defined on that scale.
2. A parameter (name and default agreed with the parameter conventions; the default is greedy, so
   every other caller is unchanged), read by `assign_state_labels`: with temperature T > 0 a
   particle's state is drawn from the softmax over states of its best score in each state divided
   by T; projection directions are drawn within the drawn state as now. Thread safety: per-particle
   draws from per-thread generators seeded so that the result does not depend on the thread count;
   a tester at 24 threads, repeated, against one thread.
3. `solve3D` sets T for stages 3 and 4 only (stage 4 at half the stage-3 value unless the sweep
   says otherwise); stage 5 stays greedy (ruling 6).
4. Sweep three stage-3 values (for example 0.5, 1 and 2 on the scale found in step 1), four runs
   each; then eight runs at the best value. Kept only if the separation count is not lower than
   Phase 3's and the time per run is unchanged.

Exit: the temperature is in the tree with its values and evidence, or it is reported as not
helping and left out.

### Phase 5: Cyril's path as a test, final checks, documents, report

- The only de novo multi-state test from particles runs shared memory with three states and no 2D
  start. Add the smallest route that covers Cyril's path: distributed (`nparts=2`), from a `solve2D`
  project, two states, as a sub-suite of `multistate_workflows` outside the default CTest entry (its
  core reservation would change), with a `success` row held by `test=profile`. Run
  `scripts/check_test_registry.py`; keep the CTest budget and the policy's counts consistent if an
  entry changes.
- `test=nightly` on the final tree: passes as on the base.
- The three-state test: `simple_test_exec test=profile target=multistate_workflows suite=shc
  nseeds=20 nthr=8 ncunits=24 keep_outputs=no`, alone on the machine, compared by count with the 13
  of 20 of 2026-10-19 (the test code has changed since, so runs cannot be paired). Fewer than 13 is
  reported with the discordant routes; 9 or fewer means the final tree is worse and stops the run.
- Timings (ruling 1): the mean time per stage of the final H runs against Phase 1's, and
  beta-galactosidase case A once on base and final tree (the speed week's single-state case:
  `solve3D pgrp=d2 mskdiam=180 nsample=2500 cavg_ini_ext=yes nparts=4 nthr=6` on
  `bgal_addon.simple`, 5 513 particles, without the retired `prob_fast` key). Kernel, table and
  assembly times may not grow; stage-count changes from ruling 3 are listed.
- Update policy documents where behaviour changed (solve3D, refine3D, the empty-state handling).
- Write the report next to the completed plan. It reads on its own.

## Common rules

- Nothing is committed; git is read-only; other commits are built from `git archive` exports.
- Fix as you go. A fix that changes results waits until Phase 3, after the localisation, so that it
  does not move the comparisons it would confound; a defect that blocks runs (a crash) is fixed
  when found.
- Thresholds are fixed before the runs they judge; never loosened, never rerun until a run passes.
- New or changed OpenMP code gets a tester at 24 threads, repeated, compared bit for bit with one
  thread; no string temporaries inside parallel regions; work arrays sized by the team in use.
- No `goto`; THROW macros never continued with `//&`; no `.and.` short-circuit on allocatable
  arrays.
- Bulky outputs of a run are deleted once its numbers are in the journal; `scratch/keep/runs/`
  keeps the labels, logs and final maps of every run, and `scratch/keep/` the data, P2 and P2H.

### File table

| file | why | phase |
| --- | --- | --- |
| `doc/implementation_notes/planned/multistate_solve3D_regression.md` | this plan and its progress | each |
| `src/main/strategies/parallelization/simple_refine3D_strategy.f90` | stage-1 references and empty states in both modes | 3 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the temperature in `assign_state_labels` | 4 |
| `src/main/solve/simple_solve3D_controller.f90` | the temperature for stages 3-4; multi-state resolution-limit rules if Phase 2 points there | 3 4 |
| `src/main/params/simple_parameters.f90`, `src/main/params/simple_parameters_parse.f90`, `src/main/params/simple_parameters_phases.f90` | the temperature parameter | 4 |
| `src/main/ui/simple/simple_ui_refine3D.f90` | the temperature parameter's metadata | 4 |
| `src/main/commanders/test/simple_commanders_test_highlevel_msworkflows.f90`, `src/main/commanders/test/simple_test_registry.f90`, `doc/policies/test_environment_policy.md` | Cyril's path as a sub-suite | 5 |
| `doc/implementation_notes/completed/multistate_solve3D_regression_report.md` | the report | 5 |

Other files a fix needs are added to this table, with the reason, by the phase that finds the
defect.

## Message to Cyril

> Thanks Cyril, we are on it and will spend a couple of days on the Dell with your simulated set.
> To aim that time well, could you send:
>
> 1. The command lines of the Sov runs (old and new), and the line `**** SIMPLE Git Commit:` near
>    the top of any old log you still have.
> 2. From `solve3D_multi.log`, the output of
>    `grep -E 'Git Commit|ENTERING REFINE3D STAGE|STAGE +[0-9]+ WITH|OMITTING STAGE|RESOLUTION EXIT|JOINT DISTRIBUTION OVERLAP' solve3D_multi.log`
> 3. How many cores the machine has (the script runs three parts of 16 threads).
> 4. A check of the labels: which column of `out.txt` holds `state=` depends on how SIMPLE writes
>    the line, so `grep -o 'state=[0-9]*' out.txt` is safer than column 17. Do the two agree?

## Progress

Nothing run yet.
