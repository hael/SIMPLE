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
  *Corrected in Phase 0 (see Progress):* docked with `dock_volpair`, the two maps agree more
  closely than this, FSC 0.88-0.90 at 30-50 Å, 0.81 at 13 Å and about 0.75 from 7 Å to 4 Å, so the
  states differ less than assumed here.
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
| `src/main/solve/simple_solve3D_controller.f90`, `src/main/solve/simple_solve3D_utils.f90`, `src/main/commanders/simple/simple_commanders_solve3D.f90`, `src/main/solve/simple_solve3D_reseed_tester.f90`, `doc/policies/3D/solve3D_policy.md` | the maintainer's stage-skipping fix (no stage omitted or cut short when `nstates>1`): received uncommitted, built and checked; edited only if it does not build or its unit test fails | 6 |
| `doc/implementation_notes/planned/multistate_solve3D_regression_report.md` | the review report of run `ms_fixcheck` | 6 |
| `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` | the benchmark's middle thread team: half the full team instead of the Mac's five performance cores (test side only) | 7 |
| `doc/implementation_notes/planned/multistate_stage5_prob_check_report.md` | the report of the check on the current code (run `ms_prob5check`) | 7 |

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

### Phase 0: builds, data, controls, first runs (2026-10-21, done)

*Made by the stopped run `ms_regression`; this section is carried over unchanged from that run's
working tree so that Phase 6 below can refer to it.*

Files changed: this note only (this section, and a correction under "The simulated set"). Evidence
lives under `/home/elmlundho/agent_runs/ms_regression/scratch/keep/` (shortened to `keep/` below);
the session journal is `/home/elmlundho/agent_runs/ms_regression/state/journal_phase_0.md`.

**What was learned.** On the first pair of runs the current build (HEAD) fails and the anchor
(1bbbdcc24, 2026-10-15) separates the two states. In both builds nothing separates during stages
1-2 (accuracy about 0.51 after stage 2). The anchor then runs stage 3 in full and all 15 iterations
of stage 4; the states start to split in the fifth iteration of stage 4 and are 0.86 apart by its
end, and stage 5 completes the separation (0.97). HEAD skips stage 3 (stage omission: the best
state's half-map FSC=0.5 resolution was already finer than stage 3's limit), ends stage 4 after
three iterations (the mid-stage resolution exit), and in stage 5 the second state drains from 6 178
to 939 particles. One pair is not a verdict, but this is the signature the plan lists for
resolution-limit updates cutting stages 3-4 short. The two simulated structures are also more alike
than the plan assumed (below), which makes a long probabilistic phase more necessary.

**Builds.** HEAD (the run's tree, 11968c20 = base 72930a3a4 plus this note) in `repo/build` with
`./compile_clean.sh` (Release, fast gate 15 of 15 passed, installed); the anchor exported with
`git archive` and built in `keep/build_1bbbdcc24` (Release, tests on). Both with gfortran 15.2.1
from gcc-toolset-15 and CMake 3.26.5; Python 3.12 first on `PATH` for the gate. The login
environment points at another checkout, so every run sources `keep/tools/env_head.sh` or
`env_anchor.sh`, which also unset `SIMPLE_SEED`. A run of an exported build prints an empty
"SIMPLE Git Commit" line (an export has no hash); such runs are identified by their build directory.

**Data.** `generate.sh` with thread counts capped at 24 and the import/2D/3D block removed,
run with HEAD and no `SIMPLE_SEED` in `keep/data` (5.5 min, `scratch/logs/generate.log`). The first
lines of `oris_1.txt` and `oris_2.txt` differ, so the two simulations drew independent poses.
`out.txt` has 30 000 lines; taken by key, 17 000 are state 1 and 13 000 state 2, and column 17 is
the `state=` field on every line with HEAD's writer, so `labels.txt` agrees with the keyed labels.
The check proposed to Cyril, `grep -o 'state=[0-9]*'`, also matches `mi_state=` and counts every
line twice; the key has to be taken as a whole field (a token starting with `state=`). Keyed labels:
`keep/data/labels_true.txt`.

**Accuracy tool.** `keep/tools/acc.py` (reads outputs only). On the true labels it gives 1.000;
all particles in one state 0.567; shuffled labels 0.508 (the expected chance level
(17²+13²)/30² = 0.509); swapped state numbers 1.000; 100 particles in state 0 count as wrong
(0.997) and are excluded with `--active`. Per-run tables come from `run_metrics.py`; one line per run
goes to `keep/runs/summary.tsv`.

**Truth control.** `import_particles` resets every state to 1, also when `oritab=` gives states,
so the truth project takes its states from `selection oritype=ptcl3D infile=labels_true.txt`
(which writes the file's values as states). With the true poses, shifts and CTF from `out.txt`,
`reconstruct3D nstates=2` gives maps that match their input maps by FSC (Fourier shell correlation,
same frame, no docking): 0.999 at 40, 20 and 10 Å, 0.98 at 6 Å, 0.93-0.94 at 4.2 Å (Nyquist 4.11 Å),
and about 0 against the other input (the inputs are in different frames). So the labels, the
simulation and the shift convention are right. With shuffled labels each map correlates about 0.80
with the state-1 input and 0.58 with the state-2 input in every shell: a 17:13 mixture. The FSC=0.5
resolution of a mixture against a truth stays at Nyquist, so found maps are judged by docked
correlation and FSC values per shell, not by that crossing alone. (`keep/truth/`, comparisons in
`keep/truth/cmp/`.)

**How different the two states are.** Docked onto each other (`dock_volpair` between 40 and
10 Å, the mirror tried with `volops mirr=x`), the truth maps correlate 0.88 in real space in the
direct hand (0.15 mirrored), with FSC 0.94 at 123 Å, 0.88-0.90 at 30-50 Å, 0.85 at 15 Å, 0.81 at
13 Å and about 0.75 from 7 Å to 4 Å. Docking the input maps directly gives the same. The earlier
estimate (0.48 at 13 Å, about 30 % of the density different) was wrong: about a fifth of the signal
differs, at every resolution.

**Known-reference control.** `refine3D nstates=2 refine=prob` for two iterations, distributed
(`nparts=3 nthr=8`), on the truth project with the two truth reconstructions as references and
HEAD's stage-4 settings on this set (low-pass 7.0 Å, images cropped to 160 pixels, 1 000
projection directions, no nonuniform filtering), every particle updated: accuracy 0.996 after the
first iteration, 0.995 after the second (ARI 0.98), 5.5 min (`keep/known_ref/`). This is the
ceiling for this data. It is above 0.95, so the **threshold stays as planned and is now fixed**:
a run has *separated* at final accuracy ≥ 0.90, is *partial* from 0.75 to 0.90 and has *failed*
below 0.75.

**2D starts.** P2: `import_particles deftab=` and `solve2D ncls=10 mskdiam=270 nthr=24` with the
anchor (501 s); P2H: the same with HEAD (364 s). In both `solve2D` leaves every particle active
(`keep/P2`, `keep/P2H`).

**First runs** (`solve3D pgrp=c1 mskdiam=270 nparts=3 nthr=8 nstates=2 print_states=yes` on its own
copy of P2, alone on the machine; labels after every stage and iteration, logs and final maps in
`keep/runs/H_1` and `keep/runs/A_1`):

| stage | H_1 (HEAD) | A_1 (anchor) |
| --- | --- | --- |
| 1 (`shc`, 20.6 Å, crop 88) | 12 iterations, 79 s, accuracy 0.52 | 12 iterations, 115 s, 0.51 |
| 2 (`shc`, 8.2 Å, crop 112) | 20 iterations, 848 s, 0.51 | 20 iterations, 949 s, 0.51 |
| 3 (`prob`, 7.4 Å) | omitted: "planned LP 7.4 Å already reached at FSC=0.5" | 17 iterations, 2 686 s, 0.51 |
| 4 (`prob`, 1 000 directions, crop 160) | 7.0 Å; 3 iterations, 387 s, ended by the resolution exit (FSC=0.5 at 4.6 Å past the next limit, no longer gaining); 0.53, populations 23 822 / 6 178 | 6.0 Å (7.0 promoted by FSC=0.5); 15 iterations, 4 915 s; populations dip to 24 085 / 5 915 at the fourth iteration, then accuracy rises every iteration from 0.50 to 0.86 |
| 5 (2 500 directions, 6.0 Å) | `prob_neigh` with NU filtering, 12 iterations, 1 529 s; state 2 drains to 939 particles; 0.57 | `prob`, 12 iterations, 6 699 s; 0.93 → 0.97 |
| whole run | 2 854 s; **failed**, accuracy 0.573, ARI 0.006 | 15 375 s; **separated**, accuracy 0.971, ARI 0.886, populations 13 414 / 16 586 |

H_1's final map of state 1 docks onto both truth maps at correlation 0.97 (a mixture of the two
structures) and its map of state 2 onto neither (0.13 and 0.23). Each of A_1's maps docks onto its
own truth at 0.99, and the two came out in opposite hands (one matches its truth mirrored, the
other directly), as two independently refined states may. Neither run printed an empty-state or
relabelling message. Both builds warn twice at the start of stage 1 that "your input doc has
multiple states but NSTATES is not given" (`simple_builder.f90`): some program of the distributed
stage-1 iteration runs without `nstates`; Phase 2 or 3 traces it.

**Rescaled plan** (written in the journal before Phase 1): an anchor run takes 4.3 h against 48 min
for HEAD, so eight anchor runs alone would take 34 h. Arm A gets 4 runs in all, judged on final
labels as planned (the decision rests on count and signature together, and the anchor's signature
is the expensive part); arm H keeps 8; shared-memory runs are cut to 2 and end-to-end runs to 1 (the
plan's drop order). Phase 1 then takes about 21 h. HEAD variants without stage omission or the
resolution exit are expected to take 1-1.5 h per run; Phase 2 runs 4 per variant unless the first
runs disagree, and at most two bisection probes, each stopped once the stage that shows the
signature has written its labels.

Exit criteria: the truth maps match their inputs (FSC 0.999 to 10 Å); the labels are confirmed by
key and by column; the threshold is fixed at 0.90 / 0.75 (control 0.995); P2 and P2H exist with no
inactive particle; the time per run is known (HEAD 48 min, anchor 4.3 h); the rescaled plan is
written.

### Phase 6, 2026-10-21: the stage-skipping fix on the two-state set

Run `ms_fixcheck`. This is not one of the planned phases 0-5: it makes one check of a fix the
maintainer wrote after Phase 0. The fix (uncommitted, in the five files the file table lists for
phase 6) stops a multi-state `solve3D` from omitting a stage or ending one early by the mid-stage
resolution exit; single-state runs are unchanged. The question is whether the fixed build now
separates the two states of Cyril's simulated set, as the anchor did and HEAD did not.

Evidence lives under `/home/elmlundho/agent_runs/ms_fixcheck/scratch/` (logs in `logs/`, kept runs
in `keep/runs/<run>/` with `summary.tsv`); the journal is
`/home/elmlundho/agent_runs/ms_fixcheck/state/journal_phase_6.md`.

**Build and unit test.** The tree (base 11968c20 plus the fix) built as received with
`./compile_clean.sh` (Release, tests on, gfortran 15.2.1, the recipe of Phase 0); the fast gate passed
15 of 15 and the build was installed (`logs/build_head_full.log`). `simple_test_exec
test=unit_reconstruction`, which holds the suite "solve3D_cavgs stage rules" with the new
`test_stage_skipping`, passed: 356 of 356 assertions, 19 in that suite (`unit/unit_reconstruction.log`).
No file of the fix needed an edit.

**Runs.** Arm HF: the fixed build, distributed, on P2, with exactly arm H's command (`solve3D
pgrp=c1 mskdiam=270 nparts=3 nthr=8 nstates=2 print_states=yes`, no `SIMPLE_SEED`), one run at a
time alone on the machine, with the tools of Phase 0 adapted to this run's paths. The verdict
threshold of Phase 0 holds: separated at final accuracy ≥ 0.90, partial 0.75 to 0.90, failed below.

After the first two runs (written before the next two were started):

| stage | H_1 (HEAD) | A_1 (anchor) | HF_1 (fixed) | HF_2 (fixed) |
| --- | --- | --- | --- | --- |
| 1 | 0.52, 12 it, 1 min | 0.51, 12 it, 2 min | 0.51, 12 it, 1 min | 0.52, 12 it, 1 min |
| 2 | 0.51, 20 it, 14 min | 0.51, 20 it, 16 min | 0.51, 13 it, 9 min | 0.51, 15 it, 10 min |
| 3 | omitted | 0.51, 17 it, 45 min | 0.51, 17 it, 21 min | 0.51, 17 it, 21 min |
| 4 | 0.53, 3 it, 6 min | 0.86, 15 it, 82 min | 0.52, 15 it, 29 min | 0.51, 15 it, 29 min |
| 5 | 0.57, 12 it, 25 min | 0.97, 12 it, 112 min | 0.57, 12 it, 26 min | 0.67, 12 it, 26 min |
| whole run | 0.573, 48 min, failed | 0.971, 256 min, separated | 0.572, 86 min, failed | 0.674, 88 min, failed |

(Each cell: accuracy after the stage, iterations, wall time.) Both fixed runs did what the fix
promises: no stage omitted, no resolution exit, stages 3 and 4 run to their iteration limits. Neither
separated. HF_1 looks like H_1: the states stay mixed through stage 4 and state 2 drains in stage 5
to 2 934 particles. HF_2 stays mixed through stage 4 too, but its states begin to split in the fifth
iteration of stage 5 and accuracy climbs by about 0.025 per iteration to 0.674 when stage 5's
twelve iterations run out. Unlike the anchor, neither run ran stage 4 at 6.0 Å: the anchor's
stage-boundary FSC=0.5 promotion fired because one of its stage-3 maps resolved to the edge of the
image box, while the fixed runs' stage-3 maps reached only about 11 Å at FSC=0.5, so no promotion
applies under either version of the rule.

**All five runs.** HF_3 and HF_4 on P2, then EF_1: the fixed build on P2H, the 2D project made by
the current build (Cyril's whole pipeline). Same command, same conditions.

| stage | H_1 (HEAD) | A_1 (anchor) | HF_1 | HF_2 | HF_3 | HF_4 | EF_1 (P2H) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.52, 12 it, 1 min | 0.51, 12 it, 2 min | 0.51, 12 it, 1 min | 0.52, 12 it, 1 min | 0.53, 12 it, 1 min | 0.52, 12 it, 1 min | 0.53, 12 it, 1 min |
| 2 | 0.51, 20 it, 14 min | 0.51, 20 it, 16 min | 0.51, 13 it, 9 min | 0.51, 15 it, 10 min | 0.52, 11 it, 7 min | 0.50, 17 it, 12 min | 0.50, 20 it, 14 min |
| 3 | omitted | 0.51, 17 it, 45 min | 0.51, 17 it, 21 min | 0.51, 17 it, 21 min | 0.52, 15 it, 18 min | 0.51, 17 it, 27 min | 0.50, 17 it, 21 min |
| 4 | 0.53, 3 it, 6 min | 0.86, 15 it, 82 min | 0.52, 15 it, 29 min | 0.51, 15 it, 29 min | 0.51, 15 it, 30 min | 0.52, 14 it, 34 min | 0.78, 15 it, 37 min |
| 5 | 0.57, 12 it, 25 min | 0.97, 12 it, 112 min | 0.57, 12 it, 26 min | 0.67, 12 it, 26 min | 0.51, 12 it, 26 min | 0.59, 12 it, 26 min | 0.91, 12 it, 26 min |
| final accuracy, ARI | 0.573, 0.006 | 0.971, 0.886 | 0.572, 0.009 | 0.674, 0.121 | 0.515, 0.000 | 0.589, 0.028 | 0.912, 0.679 |
| populations (17 000 / 13 000 true) | 29 061 / 939 | 13 414 / 16 586 | 27 066 / 2 934 | 18 652 / 11 348 | 11 879 / 18 121 | 21 749 / 8 251 | 19 311 / 10 689 |
| stage 4 limit | 7.0 Å | 6.0 Å (promoted) | 7.0 Å | 7.0 Å | 7.0 Å | 6.0 Å (promoted) | 6.0 Å (6.7 promoted) |
| wall time | 48 min | 256 min | 86 min | 88 min | 83 min | 100 min | 100 min |
| verdict | failed | separated | failed | failed | failed | failed | **separated** |

ARI is the adjusted Rand index (0 random, 1 perfect). A stage that ends before its iteration limit
(stage 2 in most fixed runs, stage 3 in HF_3, stage 4 in HF_4) ended by its own convergence rule; in
stage 2 that is the scan exit of 10-16 (more than 95 % of the search scanned). No run printed
`SOLVE3D OMITTING STAGE` or `STAGE RESOLUTION EXIT`, an empty-state or relabelling message, an error
or a NaN; every run printed the two known "NSTATES is not given" warnings of the power-spectrum step
at stage 1.

**Verdict count.** The fixed build separated 1 of 5 runs: 0 of 4 on P2 and 1 of 1 on P2H. Before the
fix, H_1 failed; the anchor's one run, A_1, separated. The plan's reading of counts this small holds:
they separate large differences only, so the signatures carry the weight.

**Where the states split.** In the two separated runs the split has the same shape. In A_1 state 2
shrinks to 5 915 particles by the fourth iteration of stage 4 (at 6.0 Å) and then accuracy rises
every iteration, from 0.50 to 0.86 by the end of the stage; stage 5 completes it. In EF_1 state 2
shrinks to 2 874 particles by the third iteration of stage 4 (also at 6.0 Å), the split starts in the
fifth and sixth (0.56 → 0.78 by the end of the stage, 0.865 after the first iteration of stage 5,
0.912 at the end). Of the failed fixed runs, HF_2 starts the same kind of split only in the fifth
iteration of stage 5 and runs out of iterations at 0.674 while still gaining about 0.03 per
iteration; HF_4 gains weakly in stage 5 (0.54 → 0.59, flattening); HF_1 loses state 2 in stage 5 as
H_1 did; HF_3 never moves.

**What a run costs.** 83 to 100 minutes (mean 89 on P2, 100 on P2H), against H_1's 48 minutes and
A_1's 4.3 hours. The extra 40 minutes over H_1 are stages 3 and 4 run in full (about 20 and 30
minutes); per iteration the fixed build is as fast as H_1 (stage 3 75-94 s, stage 4 116-150 s,
stage 5 about 130 s, against the anchor's 158, 328 and 558 s). Ruling 1 expected about 30 minutes
per run; H_1 already took 48, and the further 40 minutes come only from the multi-state rule change,
which the plan's reading of ruling 1 allows.

**Unexpected.**

- HF_3 never found the structure. Its final maps match neither truth map (docked correlation at most
  0.39 between 40 and 10 Å, at most 0.29 between 60 and 20 Å, both hands) nor each other (0.12), and
  their projections are streaky blobs (`scratch/diag_HF3/`). Its half maps agreed to 5.4 Å at
  FSC=0.143 in stage 4: a consistent wrong model. Its intermediate maps were not kept (the per-stage
  images were added to the collector from HF_4 on), so the stage where it went wrong is unknown;
  stages 1 and 2, where the de novo model forms, are not touched by the fix.
- After stage 3 the smaller state's map can be a striped blob holding 5 000-6 000 particles whose
  labels barely move (in HF_4's stage 3, 18 particles changed state over the last 14 iterations), and
  recover into a proper structure at 6.0 Å in stage 4. It happened in HF_4, which failed, and in
  EF_1, which separated, so it does not by itself decide the outcome.
- The FSC=0.5 promotion of stage 4 to 6.0 Å fired in the anchor run, HF_4 and EF_1, and not in HF_1
  to HF_3, whose stage-3 maps reached only about 11 Å at FSC=0.5. A read-only check of the code found
  no change in how the half-map FSC is measured since the anchor; the promotion rule itself changed
  on 10-18 (one rung instead of the measured resolution), but at 11 Å neither version promotes. Both
  separated runs had stage 4 at 6.0 Å; HF_4 had it too and did not separate.

**Answer.** The fix does what it says: multi-state stages are no longer omitted or cut short. It
does not restore the anchor's behaviour on its own: on the anchor's 2D project the fixed build
failed 4 of 4 times, and the one separated run came from the current build's 2D project. Where it
separates, it does so like the anchor (in stage 4 at 6.0 Å, completed in stage 5) at 100 minutes
instead of 4.3 hours. What makes the difference between a run that splits in stage 4 and one that
does not is not settled by five runs. Candidates the evidence leaves open: whether stage 4 runs at
6.0 Å (it did in both separated runs, and in one failed run); the stage-2 scan exit (the two
separated runs, and H_1, ran stage 2 for its full 20 iterations; the four fixed P2 runs stopped after
11 to 17); the prob stages' greedy state assignment (stages 3-4 of the failed runs run their full
length without gain, the plan's cause 4); and the 2D project (P2H gave the one separated fixed run,
but one run cannot separate this from chance).

**Phase 6 exit.** The build passed (fast gate 15 of 15) and its unit test passed (356 of 356, with
`test_stage_skipping`); the five runs were made, one at a time alone on the machine, and judged
against the fixed threshold (1 of 5 separated); this section and the report
(`multistate_solve3D_regression_report.md`) are written. Files changed in this phase: this note
only (Phase 0 carried over, this section, the phase-6 rows of the file table) and the report; the
fix's five files are as received.

### Phase 7, 2026-10-22: multi-state prob in stages 3-5 on the two-state set; the table profile; the benchmark on the Dell

Run `ms_prob5check`, on HEAD ee5b68be1 (no stage skipping for multi-state, multi-state `prob` in stages 3-5 with
`prob_neigh` from stage 6, the vectorized per-pair polar kernels). Report:
`multistate_stage5_prob_check_report.md` beside this note. Journal:
`/home/elmlundho/agent_runs/ms_prob5check/state/journal_phase_7.md`. Evidence under
`/home/elmlundho/agent_runs/ms_prob5check/scratch/` (abbreviated `scratch/` below).

Files changed: `src/main/commanders/test/simple_commanders_test_highlevel_matchbench.f90` (the benchmark's middle
thread team: half the full team instead of the Mac's 5, so 12 on the Dell and still 5 on the Mac; test side only),
this note (file table, this entry), the report. No production code.

**Answers.**

1. Separation on P2 (the anchor's 2D project), four runs, same command as before: 1 separated (0.965), 1 partial
   (0.806, still rising when stage 5 ended), 2 failed (0.548, 0.510). Earlier: the build with `prob_neigh` in
   stage 5 0 of 4 on P2 (1 of 5 with its P2H run), the anchor 2 of 4. Four runs per build do not rank the builds.
   With `prob` in stage 5 a split that starts late in stage 4 continues through stage 5; a minority that drains in
   stage 4 (H5_2) or labels frozen after a stage-3 lock-on (H5_3) do not recover.
2. The probability table is 86-93 % of each prob stage's wall time (77-94 % of single iterations). Of its fill's thread time, scoring is
   67 / 57 / 68 % in stages 3 / 4 / 5 (pooled; 51-88 % per run and stage), the shift refinement of the best
   references 31 / 41 / 31 % (11-48 %), the seed 0, the rest 1.3-1.6 %. Scoring is about half of a prob iteration
   (49-60 %); a 3.4-4.1 times faster scoring would make prob iterations about 1.5-1.8 times faster.
3. On the Dell the new euclid kernel is slower than the old (one P-core 15.6 against 14.4 µs per pair, one E-core
   29.2 against 23.7, 0.75-0.96 times the old speed at 12 and 24 threads, 0.80-0.85 in an alternating-order repeat); the cc
   kernel is 1.18 times faster on a P-core, 0.92 on an E-core. One core is limited by arithmetic (cost flat from
   1,000 to 10,000 references), the full machine partly by memory (cost per pair and thread doubles and grows from
   31 to 43 µs with the references). `blocked` gains 1.7-2.1 times at 24 threads (nothing on the Mac),
   `matrix_t` 3.4-4.1 times at 24 and 3.8-4.9 at 12 (3.7 on the Mac).

**Runs** (`scratch/keep/runs/H5_<n>/`; per-stage accuracy, iterations, minutes in the report):

| run | stage 2 | after stage 3 | stage 4 | stage 5 | final | minutes |
|---|---|---|---|---|---|---|
| H5_1 | 13 it (scan exit) | 0.52, no lock-on (poses 23° / 16°) | 7.0 Å; poses right (2.3° / 2.6°), split starts in its last 4 iterations | 0.55 → 0.81 | partial 0.806 | 160 |
| H5_2 | 16 it (scan exit) | 0.51, no lock-on (24° / 18°) | 7.0 Å; poses right, minority drains 11.1k → 7.8k | flat 0.54-0.55, minority 5.2k | failed 0.548 | 141 |
| H5_3 | 20 it | 0.51, lock-on: map 0.97, poses 8.2° / 2.7°, promoted to 6.0 Å | labels frozen at 23.9k / 6.1k after an early swing of about 1k | frozen | failed 0.510 | 140 |
| H5_4 | 14 it (scan exit) | 0.51, half-maps at the crop edge and promoted, but no right map (≤ 0.70, 22° / 15°) | poses right, split from its 8th iteration, 0.50 → 0.58 | 0.71, 0.80, 0.88, 0.92 … 0.965 | separated 0.965 | 165 |

No stage omitted, no early exit, no empty state or relabelling in any run; final pose errors 0.9-1.2°. Time per run
140-165 minutes against 83-100 for the build with `prob_neigh` in stage 5: stage 5 now takes 66-91 minutes (5.5-7.6
per iteration over 5,000 references) against 26.

**Decisions** (journal): Release builds by `compile_clean.sh`; the middle team change; a pinned benchmark series
(the kernel 4.18 does not know the P/E layout, so unpinned one-thread rows mix core kinds), without its old-build
`matrix_t` runs (that layout does not call the changed file); an alternating-order repeat of old against new; a
vectorization report of both kernel versions (the old loop was already AVX2-vectorized on x86); snapshots after
stage 4 as well; inputs read in place from `ms_regression`, tools copied from `ms_stage3diag` and adapted.

**Exit items.**

- Builds pass: HEAD in `repo/build` and the HEAD export with only `simple_polarft_corr.f90` at f3ccb5fbd in
  `scratch/src_oldkernel/build` (tree compared file by file), Release, fast gate 15 of 15 each, again 15 of 15 after
  the test-side change (`scratch/logs/build_*_full.log`, `fastgate_repo_after_edit.log`). BLAS: system OpenBLAS
  0.3.15 (`/usr/lib64/libopenblas.so`).
- Benchmark measured on both builds: `per_pair`, `per_pair_cc`, `replica`, `blocked`, `matrix_t` over the six
  cases at 24, 12 and 1 threads, best of three after a warm-up, `OPENBLAS_NUM_THREADS=1`, alone on the machine
  (`scratch/bench/{head,old}/polar_matching_bench.tsv`, tables `bench_tables.md`); the pinned series
  (`scratch/bench/pin_*`), the repeat (`scratch/bench/rep*`); `suite=check` on HEAD: every layout's 24-thread
  repeats bit-identical to one thread, `matrix_t` within 4.0e-7 relative (`scratch/bench/check_table.md`).
- Four runs made and judged: H5_1-4 above, one at a time on the whole machine, no `SIMPLE_SEED`
  (`scratch/logs/run_H5_<n>.log`).
- Profile extracted: `scratch/keep/runs/H5_<n>/profile.txt`, pooled in `scratch/profile_table.md`
  (`scratch/keep/tools/profile.py`, `profile_table.py`).
- Report written: `multistate_stage5_prob_check_report.md`.
