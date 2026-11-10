# How to write a plan for an agent in SIMPLE

Guidelines for the notes under `doc/implementation_notes/planned/` and
`doc/refactoring_notes/planned/` that a Claude Code session or an unattended
`arun` run carries out. They are not a style preference. Each one is backed by
what happened when a plan in this repository had or lacked it: about 60 plans
written between July and November 2026 were compared with their run records,
their reports and the code that is in production today (section 9 lists them).

The short form:

1. Pin the base commit and check every premise against the code at that commit.
2. Say where the code will live, what it depends on, and what is not part of it.
3. Put the decisions and rulings in one numbered block at the top, with numbers in them.
4. Give every phase an exit the agent can check, and a file table with a phase column.
5. Phase 0 measures the base: build, gate, the production entries twice, the acceptance rule on unchanged code.
6. Tag every floor and tolerance with where its number came from, and freeze it before the run it judges.
7. Identical changes are proven byte for byte at one thread; result-changing changes are judged against a before-profile over seeds.
8. Write the stop conditions and the list of gates no ruling may waive.
9. Nothing is "implemented" until it has been built and run; nothing is a default before a real-data read.
10. End with a report: what was learned, decisions taken under delegation, deviations, evidence, what is left.

## 1. What the record shows

The comparison separates cleanly along a few lines.

| Correlated with clean delivery and code that survived | Correlated with stops, rework, unfinished splits or removal |
| --- | --- |
| A measured Phase 0 baseline and a tolerance derived from the base's own spread (`prob_align_batch_pipeline`: "twice the base's largest observed difference"; `structure2D_and_project_owner`: two base runs per case) | An acceptance rule written without measuring spread: the base itself failed it (`global_sigma2_only` first stop; `particle_io` solve2D rule) |
| Byte identity at one thread with a fixed seed and declared exclusions (`identical_speed_run`, `speed2D_run`, `speed2d3d_run`) | Identity judged on a thread-order-dependent metric with a two-run spread (`multistate_speed_run` Phase 2 stalled) |
| Decisions closed before launch, numbered and quantified (`flex_on_simple_machinery` D1-D12, the stream fix plan's 32-row table) | Open decisions "taken as recommended" under delegation (`pose_cont_refactoring` O4-O8: two reversed at review) |
| A change-type table naming the gate per phase (`cont_lm_acceleration`: identical / endpoints change / rounding changes) | Floors from an emulation or a closed form applied to real maps whose physics the design had not measured (`species_discovery`: four stops in one phase) |
| Pre-registered verdict thresholds with a measured calibration control (`multistate_posterior_state_draw`) | Verdicts on three to ten restarts without a decision rule (`pose_marginalization`: "three of seven against two of ten") |
| A statement of what is reserved for the maintainer, with the artifact to collect (`nano_pipeline_pt`: "the maps are not judged here") | A gate the agent cannot judge left as a gate ("rounder atoms and a cleaner background", `single_pcg_integration`) |
| Code built and run before a phase is marked done | Phases written without a compiler (`flex_on_simple_machinery` 5-6, `prob_search_speedup`, the stream plans of 5-8 October: "nothing has been compiled or run") |
| A new rule backed by a gate script (`check_cli_chain.py`, `check_guide_imports.py`) | A rule in a policy only: the compile-time rules were "widely violated with no check" four weeks later; dated comment lines regrew from 5 to about 90 |

Two things did not predict survival: the length of the plan and the weight of its
process. The heaviest record in the set (a frozen contract with twelve acceptance
criteria, eleven phases and checksummed evidence) and the lightest (implemented the
same day, nothing compiled) both ended in `rejected/` within weeks. What predicted
survival was a stated scientific endpoint, an integration owner and a measured base.

## 2. Before writing

**Pin the base.** Name the commit the plan was written against. Reviews and plans
that cite `file:line` at a pinned commit could be re-verified by the next review;
"the current tree" could not. Routine and module names are the durable anchor; line
numbers must say "at commit X".

**Check the premises against the code.** Run `rg` for every routine, key, file and
constant the plan names, at the pinned commit, and write down what the check changed
(`identical_speed_run` "What the code check changed"; `speed2d3d_run` "What the
reading found" corrected two premises before work started). Plans written from
reports and memory carried wrong facts in four of seven refactoring runs: a wrong
base commit, a constant that differed from the code, a stale claim that a name did
not exist, a fixture size that did not survive x86.

**Name the integration owner.** Say which layer and which existing routine the new
code extends, following `ui -> exec -> commander -> strategy/domain object` and
`doc/policies/layering_and_project_ownership_policy.md`. Both features in the record
that were rewritten wholesale within weeks had landed as a parallel route beside the
owner: a second program, an adapter with eleven public types, a strategy that could
not satisfy the strategy contract. A validation contract that forbade integration
("no production caller") delivered correct numerics with no home, and the home had
to be designed afterwards.

**Name the unproven dependencies.** Write "depends on X (validated: yes/no)". Of the
nine rejected plans that were never implemented, five were written on top of a
component that was itself removed within days or weeks. Nothing in those plans said
so; the banners added later did.

**Check the queue.** Say whether the layer the plan touches is slated for redesign
and whether another run is editing the same files. About half of the identical
optimisations of three October runs lived seventeen days before the batch pipeline
retired the tables they optimised; a plan written against `ref_assign` while a
concurrent run deleted it needed a page-long re-mapping ruling and still left stale
names for the final review. Sequence the runs, or write the plan after the other one
lands.

**Set the principles before the inventory.** For a cleanup, write what is kept on
purpose, what is out of bounds and which policy decides, before listing candidates.
The release-4 cleanup did and every deletion held. The dead-code inventory of
2026-11-02 extrapolated from a compile-time motivation to a deletion list, called
user-run programs and general-purpose APIs dead, and was refused. The rule that came
out is in the `simple-dead-code` skill: grep finds candidates, the compiler probe
proves users, the owner decides for general-purpose APIs.

## 3. The shape of the plan

Write the plan once. Its history goes into Progress rows and the report, not into
rewrites of the plan text. Notes with twelve to sixty-six commits were plans that
absorbed their journal; they needed retrofitted banners ("as written 2026-08-29;
superseded, see the banner") and later a consolidation session. The two arun runs
written in the same weeks have two or three commits each and have not been touched.

In order:

1. **Status banner.** One line: delivered / planned / superseded by X, with the date.
   When a plan supersedes another, say so in both, in the same commit.
2. **Why.** One paragraph with the measured motivation (a run time, a failure count,
   a resolution), not an estimate table. Estimate tables were wrong in direction and
   size every time they were later measured (section 6).
3. **Contract.** The invariants the change keeps, including the off-state invariant:
   with the feature off, which lines of the new code may run? `inpl_cont=no` was not
   the legacy path because "objective selection was used as feature selection"; the
   missing test was the negative test.
4. **Decisions and rulings.** One numbered block, dated, with the quantities in it
   ("max(5·nstates, 10 % of the frozen particles)", "retry once in a fresh folder").
   Never renumber. Later plans cite and replace by number ("replaces fix plan decision
   26"), and the quantities survived in code and policy where the plan around them was
   superseded. Add "Defaults taken without a question" for names and constants the
   maintainer may veto cheaply; none of those defaults caused a rollback.
5. **Not in this plan.** Non-goals and the files out of bounds. Every arun plan that
   had this section kept to it; scope creep in the record is confined to interactive
   notes. A non-goal in a sentence lasted ten days ("cc explicitly out of the first
   implementation"); the enforceable form is the file table plus the stop rule.
6. **Phases, each with "Exit:".** An exit is a predicate the agent can check: a test
   name and its result, a byte comparison, a count, a file that must or must not
   exist (`rg` for six retired symbols). "Phase 1's exit required the tree to hold no
   `sigma_est`" is an exit; "the design is in place" is not. Each planned phase of ten
   to twenty-five lines produced 65 to 430 lines of Progress, so keep phases to one
   deliverable each.
7. **Change-type table.** For each phase: identical / result-changing / tests-only,
   and the gate that follows from it (section 5). The plan that had this ran seven
   stages without a stop and every gain is attributed to one stage.
8. **Common rules.** No new parameter unless ruled (never breached in the record;
   experimental variants are patches on a `git archive` export); no floor, tolerance
   or assertion relaxed unless the plan or a ruling says so; fix as you go, with a file
   table row before the edit; the thread policy; the disk policy.
9. **Stop conditions.** Short and about design, not about failures: a floor that would
   have to move, a file outside the table, a fixed decision the code contradicts, a
   fixture the machine cannot run. Also list the gates no ruling may waive. A ruling
   that reduced a phase's gate to the unit suites shipped a result-changing fix
   untested on the identity cases.
10. **File table.** Path, purpose, phases that may edit it. A folder row suits a
    deletion-heavy phase; a retirement phase owns every file that names the retired
    symbols (the driver killed `prob_pipeline` after Phase 4 on exactly this). A row
    for a file that may not need an edit says "if needed". The report lists every file
    changed outside the table.
11. **Progress.** One entry per phase, written by the run: date, machine, base,
    files changed, each exit item and how it was met, log paths, decisions taken under
    delegation, "what it means, in short". A later entry that supersedes an earlier
    claim says so.

## 4. Phase 0 is a measurement, not a warm-up

Phase 0 runs on the base, before any edit:

- Build Release and Debug, run the fast gate, record the machine (compiler, cores,
  load) and keep a base build for later comparisons.
- Run the production entries the change touches twice with a fixed seed. That shows
  which outputs are deterministic (identity gates are set only there) and what the
  spread is (tolerances are derived from it). The balance-sampling run found this way
  that `solve3D` is not reproducible under a fixed seed, which no earlier plan had a
  mechanism to notice.
- Run the acceptance rule on the unchanged code. In two runs the base failed the
  plan's own rule; both stopped in a later phase instead.
- Record base failures and leave them. Plans that named known base failures in
  advance ("recorded, not judged") classified them without stopping; plans that did
  not stopped on a FLEX segfault and an add-on drift that were not theirs.
- Verify the command lines against the UI's required keys. One run lost its first
  phase to four commands the UI refused.
- Dry-run the fixture by hand. Three stops came from fixtures the gate's unit could
  not process: planted noise windows that do not align with class chunks, Gaussian
  atoms "the detector assumes" on maps with B-factors and negative wells, a selection
  rule and a test bound that cannot both hold.
- State the data preconditions as a Phase 0 stop ("if the labelled files are absent,
  stop and say so"). Data not in hand is the dominant reason work in the record is
  still open.

## 5. Acceptance criteria

**Tag every number.** Each floor, tolerance, band or target carries its origin:
*measured* (Phase 0 spread, a before-profile, an oracle), *literature*, or
*guessed*. In the record, nearly every stop traces to a guessed or emulation-derived
number meeting real data; measured numbers were never contested.

**Freeze before the judged run.** "Tolerances come from the measurements of Phase 0
and Phase 2, recorded with their evidence" is the single most useful sentence in the
record (`ghostbuster_operators`). "Floors only from a profile of the finished code by
the floor rule, never from a first run" is its complement
(`doc/policies/test_environment_policy.md`, section 6.2). A number printed by the code
under test is not a criterion: a wrong gradient passed finite-difference checks whose
tolerance the harness had printed.

**Identical changes.** Byte identity at one thread with a fixed seed, exclusions
declared in Phase 0 and never added because an item changed a file, every case run
twice on the base first, and "never rerun until it passes; an item that cannot be made
identical is restored and recorded". This caught three agent mistakes before any item
was judged in one run, a fixture corruption in another, and three latent defects of
the replaced path in a third. Identity across compilation units is not guaranteed in
Release; allow a ulp bound for code-generation work and demand bytes only for logic.
Byte identity does not prove thread safety; a race was fixed after identity passed.
Add the parallel-region rules of `speed2D_run.md` ("Why thread stress testing") when a
parallel region changes.

**Result-changing changes.** Judge against a before-profile of the base over N seeds,
paired per seed, with the acceptance "inside the profile's spread" or a floor derived
by the floor rule. State N. Three to ten restarts decided too much in the record; the
stage-3 diagnosis computed that its own split "would still turn up about one time in
eleven" and asked for eight to ten per arm. If N is small, label the result "a lead,
not a finding".

**Pre-register verdicts.** For a scientific question, write the decision rule before
the runs: the thresholds, the calibration control that fixes them ("if the
known-reference control keeps less than 0.95, both thresholds drop by the shortfall,
fixed before the first de novo run and never moved"), the restart count, and the
all-against-all comparison if two builds are compared. This turned anecdotes into
clean negatives twice. A fixed absolute bar with no control (3.6 Å in 10 of 10) was
abandoned.

**Bound the coverage metrics.** "The teacher's direction lies in the 99 % set" was
satisfied by a set that covered the sphere; "every block carries a polished entry" was
unachievable without changing single-state results. Write a metric as coverage and
width, and ask at writing time what would pass it trivially.

**Mutation-test the gates.** Plant each kind of violation and show the check fails,
then remove the plant (`check_cli_chain.py` was accepted on nine plants; the
box-division gates were shown to fail on the restored division). A negative control
per gate catches "checks that could not fail", which the nightly triage found in
existing test code.

**Never pass by reseeding, never relax by ruling without a diagnosis.** A seed
changed to make a truth test pass was reverted in three days. An acceptance rule
dropped by ruling after a failure masked a regression that was reverted five days
later; the agent's log had already shown the mechanism. A failing acceptance is a
diagnosis task.

**Reserve the human judgments explicitly.** Where the maintainer's eye decides (a
map, a class-average montage), write "reported, not gated" with the artifact and its
path, and say what the agent must collect. The run whose plan said "the maps are not
judged here" closed cleanly with an intended open item; the runs whose plans left a
qualitative gate ("rounder atoms") could not close.

**No default before a real-data read.** Both scientific removals in `rejected/` were
preceded by a new default set on the day the code landed, with the falsifying real-data
test named and deferred. Fixture-only validation misled on real data twice more
(a purity gain that "did not carry to the bench's class populations"; a teacher too
broad on real sets). Name the real set, where it is, and who reads it.

## 6. Performance plans

- **Estimates are hypotheses.** Write each with the measurement that will test it,
  not as a table of savings. The headline item of the speed week (estimated 400 to
  500 s) gave 8 %; polling, estimated smaller, was the real gain; the worker fold
  (estimated 9.5 %) measured neutral and was dropped; a bandwidth model was falsified
  by its own benchmark.
- **Measure first, with a reorder clause.** "If the numbers contradict the order of
  Phase 3, reorder it and say why"; "item 6 runs only if the loop dominates". Skipping
  on the stated condition is then a success, not a deviation.
- **Protocol.** Identity at one thread; interleaved before/after builds from one
  snapshot of the dirty tree; at least two or three repeats each way, with the A-A
  spread as the noise floor; an idle machine with the load recorded (one run waited
  until 02:35); pinned single-core and unpinned team timings on hybrid CPUs, with six
  or more rounds for the team; thread-seconds and load balance beside wall time, since
  they attribute gains that iteration-count variance hides.
- **Both architectures before a commit.** A kernel committed at 1.30× on the Mac ran
  at 0.80 to 0.85× on the Dell. Any compiler-shaped change is measured on both.
- **State the data's power.** Say what the dataset can and cannot measure (one
  page-cached stack cannot time I/O; a set "without real states" cannot judge a
  multi-state schedule). Where it cannot, rule "regression only" and report counts
  that must change (frozen images read 25,200 → 0) as the primary evidence.
- **Keep rule fixed before data.** "At least 5 % faster at 11 threads on the largest
  case and at most two units in the last place"; "not more than 2 % slower in any of
  four settings".
- **Schedule the real-target measurement.** Every speed plan deferred "PfCRT and msp1
  by the maintainer afterwards"; none of those measurements is in a note. Make it a
  phase with an owner and a date, or say that the gain is unmeasured.
- **Record the losers with their mechanism** ("doubles the working set", "+12 % and
  +18 %"), so they are not re-litigated. The PCG priors history and the decision log
  are the model: hypothesis, fixture, target set before the run, the measurement that
  killed it, one-sentence lesson, removal commit; and a present-tense page of what is
  deliberately not there.

## 7. Running and closing

- **Built and run, or not done.** No phase is marked done on code that has not been
  compiled and run with its exit checked. Phases written without a compiler produced
  build errors, an allocation bug, untested fixes and validation that was never
  recorded afterwards. "Completed" with "nothing has been compiled or run" is a
  directory misuse; the closing rule of every plan is "a recorded run".
- **Validation at close time or never.** Every "user-owned validation" matrix in the
  record is still open; two were archived on 2026-11-04 as "not claimed as passed".
  If the maintainer is to run something, the plan names the run and the date, and the
  close waits for it.
- **Multi-day runs need a sync rule.** A three-day run on a tree that master outran by
  nine commits forced a manual stop, a hand merge and a full re-profile. Name the
  merge point or freeze master for the run.
- **Rulings name their files.** A ruling that removed a file-table row the diff had
  touched stopped the driver. A ruling is applied where it belongs in the plan the
  first time it is used ("Ruling 4 is applied here for the first time").
- **Previous build defective.** Say what to do when a byte-identity exit fails because
  the replaced path was wrong. The agent had to rule for itself twice; its rulings were
  right, but the plan gave no authority.
- **Splits keep their tails.** When a run is split, the remaining phases go to a named
  run or are closed in the same edit. An `_unfinished` note was deleted while three
  documents still link to it; a Phase 4 "final runs and report" was never run and its
  numbers were never reused.
- **Every new rule gets a check in the fast gate**, in the pattern of
  `scripts/check_test_registry.py`, or it erodes within weeks (section 1).
- **File, do not delete.** The notes of removed code move to `rejected/` in the commit
  that removes the code, with a one-line class in the banner: idea never pursued /
  negative result / superseded / rewritten / measured out. The largest removed
  experiment of the summer (2D SGD, about 1,300 lines and nine tests) has no record
  because its notes were deleted.
- **Plan is history, policy is truth.** The plan says which policy is the living
  description and updates it in the same phase. Plans that did this aged well; plans
  that were the only description of a behaviour remain "unverified caveats".

## 8. The report

Written by the last phase as `<plan>_report.md` beside the completed plan, in this
order, because this is what later notes cited and reused:

1. **What was learned**, first, as findings with their evidence.
2. **What changed, phase by phase**, with file counts.
3. **Decisions taken under delegation**, each with its reason and its Progress entry.
4. **Deviations from the plan**, including expected effects that did not occur
   ("multislice recovers atoms the linear map misses" did not arise) and every file
   changed outside the table.
5. **Evidence against the baseline**: gates, suites, identity comparisons, the
   numbers with their spreads.
6. **Left open for the maintainer**, separated into reported-not-gated items,
   data-dependent items and follow-up runs.
7. **A commit message.**

Reports written this way made the negative results reusable: both failed multi-state
runs fed the next plan's rulings.

## 9. Interactive sessions

The same elements, scaled. A session that specifies, codes, builds and tests in one
sitting still pins the base, writes the off-state invariant, runs the production entry
before and after, and freezes any number before the run it judges. The plan is written
once; the session appends dated Progress entries and, when a decision is reversed, the
superseded entry says so in its last sentence (the rule of `pcg_decision_log.md`).
Review documents produced in a session are committed; two plans in the record cite
reviews that exist nowhere.

## 10. Template

```markdown
# <Title>

Status, <date>: planned. Base: `<commit>`. Supersedes: none. Depends on: <X> (validated: yes).

## 1. Why
<one paragraph with the measured motivation>

## 2. Contract
<invariants; the off-state invariant; what stays byte-identical>

## 3. Decisions and rulings
1. (<date>, <who>) ...          <!-- quantities in the text; never renumbered -->
Defaults taken without a question: ...

## 4. Not in this plan
<non-goals; files and folders out of bounds; what another run owns>

## 5. Change types and gates
| Phase | Kind (identical / result-changing / tests-only) | Gate |

## 6. Acceptance
<each number tagged measured | literature | guessed, and where it is frozen>
<identity cases with exclusions; before-profile seeds; verdict rule and control>
<reported, not gated: artifact and path>

## 7. Common rules and stop conditions
<no new parameter; no floor relaxed; fix as you go with a table row first>
<stops: ...; gates no ruling may waive: ...>

## 8. Phases
- Phase 0, baseline. <builds, gate, entries twice, acceptance on base, fixture dry run>
  Exit: ...
- Phase 1, ... Exit: ...

### File table
| File | What | Phases |

## Progress
(one entry per phase, written by the run)
```

## 11. Anti-patterns and the episode behind each

| Pattern | What happened |
| --- | --- |
| Savings estimate table as motivation | Headline item delivered 8 %; the unestimated item won (`single_state_speed_week`) |
| Acceptance rule without a measured spread | The base failed it; first of three stops (`global_sigma2_only`) |
| Floor relaxed by ruling after a failure | Masked a regression reverted five days later (`multistate_speed_run`, 709b9d761) |
| Seed changed to pass a truth test | Reverted in three days (61695ce7a) |
| Phases marked done without a compiler | Name clashes, a double allocation, validation never recorded (`flex_on_simple_machinery` 5-6) |
| Default on the day the code lands | Both scientific removals in `rejected/` (`local_linear_diffusion_map_pre_image`, `nu_evidence_envelope_masking`) |
| Parallel route or program beside the owner | Rewritten wholesale within weeks (`refine3D_pose_cont`, `sieve_ini3D`) |
| Validation contract that forbids integration | Correct numerics with no home for three months (`continuous_3D_pose_end_polishing`) |
| Plan on top of an unvalidated component | Five of nine never-implemented rejections |
| Coverage metric without a size bound | Satisfied by the whole sphere (`learned_pose_state_assignment`) |
| Qualitative gate ("rounder", "cleaner") | Could not close (`single_pcg_integration`) |
| One-machine code-generation measurement | 1.30× on the Mac, 0.80× on the Dell (`polar_kernel_codegen`) |
| Rule in a policy without a gate script | Compile-time rules violated unchecked; dated comments regrew |
| Plan written against code a concurrent run rewrites | A page-long re-mapping ruling, stale names at review (`multistate_posterior_state_draw`) |
| Multi-day run without a sync rule | Manual stop, hand merge, full re-profile (`test_suite_truth_and_robustness`) |
| "User-owned validation" section | Never closed, in every case |
| Deleting the notes of removed code | No record of the largest removed experiment (2D SGD) |
| Process weight instead of an endpoint | The heaviest and the lightest record both rejected within weeks |

## 12. What this rests on

The comparison was made on 2026-11-09 over 205 notes in the two lifecycle trees
(121 implementation notes, 93 refactoring notes; 18 rejected). About 60 were read in
full against their reports, Progress entries and the code at `71c4a41aa`, with `rg`
and `git log -S` for every named artifact: the autonomous refactoring runs
(`prob_align_batch_pipeline`, `pose_cont_refactoring`, `flex_on_simple_machinery`,
`structure2D_and_project_owner`, `global_sigma2_only`, `prob_search_speedup`,
`agent_friendly_code`), the science runs (`multistate_posterior_state_draw`,
`learned_pose_state_assignment`, `polar_shift_search`, `identical_speed_run`,
`cont_lm_acceleration`, `multistate_stage3_diagnosis`, `pose_marginalization`), the
performance runs (`single_state_speed_week`, `multistate_speed_run`, `speed2D_run`,
`speed2d3d_run`, `polar_kernel_codegen`, `polar_matching_layout_bench`,
`particle_io_few_particle_stacks`, `hardware_resource_planning`), the nano and wave
runs (`species_discovery`, `nano_pipeline_pt`, `trajectory_solve3D`,
`atomic_density_model`, `nano_cavgs_rejection`, `ghostbuster_operators`,
`single_pcg_integration`), the interactive-era plans of July to September
(continuous in-plane rotation, continuous pose refinement, `drop_legacy_box_division`,
`pcg_priors_history`, `uniform_test_environment_refactoring`, `abinitio3D_addon`, the
sigma2 and reconstructor refactors), the stream chain of October, the test-suite plans
and handovers, the policy-set audit, all eighteen rejected notes and the three cleanup
records.

Limits: one maintainer and one code base; the survival horizon of the November runs
is days to weeks; the arun format itself is six weeks old, so its long-term survival
is not yet measured. The guidelines should be re-checked against the record when the
next batch of runs has aged.
