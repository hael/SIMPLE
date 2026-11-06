# Agent-friendly code: report of the run `agent_friendly`

2026-11-05, on the dell, unattended. Plan: [agent_friendly_code.md](agent_friendly_code.md) (sections 3, 5 and 7
say what and why, section 8 how the run was scoped and judged). Its Progress section holds the detail of each
phase, with every number and log path. Base: `6565f5ba5`. Logs: `/data2/agent_runs/agent_friendly/scratch/`
(`scratch/` below; `~/agent_runs/agent_friendly` links to the run directory).

## What was learned

- **Checks find defects that guides only describe.** Turning the silent rules of 3.1 into stops found real
  defects along the way:
  - `rot_to_asym` returned an uninitialised record under c1 at the south pole;
  - `eo=-1` particles were summed into the even half;
  - the gridding density files had no header, so a reader of the wrong box read a prefix without a word;
  - `read_segment(..., fromto)` ignored the range for `.simple` files, and its string branch stored rows at the
    wrong index;
  - `export_cls2D` changed the project it exported;
  - icosahedral maps aligned by the symmetry search were in another setting than refinement uses (phase 4, O1).

  Each new check is a tester assertion that fails on the code before it: 134 assertions in 13 suites on the base
  for phase 2's fixes, and 8 on the previous phase's code for phase 3's.
- **A static check of the command-line chain was cheap and broad.** `scripts/check_cli_chain.py`, now in the
  fast gate, reported 196 problems on the base: 93 keys accepted and silently ignored, 94 yes/no keys taking any
  value, 4 user-interface inputs that stopped the program, and 4 private programs never checked. It passes now.
- **Behaviour-preserving work can be proven byte for byte; behaviour-changing work needs attribution.** The
  phase-3 renames and moves left every workflow output identical. One ruled fix of phase 2 (X5, the atomic centre
  at ldim/2) changes the simulated input of every workflow entry, so the entries' numbers moved. Rebuilding the
  tree without that one fix showed that nothing else changed apart from two ruled project keys.
- **One run of a de novo workflow cannot judge quality.** `simulated_workflow_6vxx` does not depend on
  `SIMPLE_SEED`, and one data set shows only "different". The 3.7 judging runs (deferred) are what would show
  better or worse.
- **Premises need measuring.** O4's ruling assumed every L-BFGS-B client passes its tolerance into the line
  search. The shift search does not. Handing it the textbook line-search constants made it worse (271 instead of
  150 evaluations, larger shift errors).
- **A Fortran language server copes with SIMPLE.** `fortls` 3.2.2 indexes all 832 files in 1.6 s at 216 MB.
  It resolves `use` chains and type-bound bindings, reaches submodule bodies in two steps, and gives exact
  references for type-specific bindings. Generic calls resolve to the binding only, and the JSON library gives
  false errors.

## What changed, phase by phase

The driver's `/home/elmlundho/agent_runs/agent_friendly/review/phase_N.diff` files for phases 0 to 5 hold only the
changes to `.github/skills` and the generated Fortran indexes: no Fortran source, no plan, no new file
(`phase_4.diff` is empty). The complete per-phase diffs for those phases were rebuilt from the driver's
phase-start snapshots and HEAD by `scratch/keep/full_diffs/make_full_diffs.py`: `scratch/keep/full_diffs/phase_N.diff`,
with a file list and line counts in `phase_N.stat`. The engine was fixed before phase 6, so `review/phase_6.diff`
and `review/full.diff` are complete. The whole change is `git diff 6565f5ba5` of the working tree plus its
untracked files.

| Phase | Files (+/- lines) | What |
| --- | --- | --- |
| 0, baseline and guides | 19 (+433/-264) | Base built (Release, Debug, a kept base export); gates 15/15; workflow baselines recorded: 1jxy stops in solve3D stage 3 (`new_shiftmap`), 6vxx misses the map (correlation 0.561, FSC 20.8 Å), mini_stream passes. Guides brought up to date with the batch pipeline (`ea4c1c20a`); links to completed notes fixed. |
| 1, command-line chain | 27 (+1383/-444) | `check_cli_chain.py` in the fast gate. R1 (`mul`, `import_type` bound), R2 (four inputs dropped or renamed), R3 (generator stops at the type's end), P6 (`add_yesno`), P3 (private programs table built and complete), P2 (stream launcher validated against the descriptors), X3 (stream aliases gone). Each kind of violation planted and caught. |
| 2, small fixes and 3.1 | 102 (+1659/-296) | All 24 small items of section 7 (V1 already done) and 9 of the 13 rules of 3.1 as stops, each with a tester check. |
| 3, medium items, renames, glossary | 101 (+9519/-8510, most of it regenerated indexes) | X2, C1, V2, N10, O11, V5; G3 in part; 380 numbered specifics renamed in five folders; `doc/glossary.md`; the `fortls` experiment (outside the repository). |
| 4, measurements | the plan only | I1, O1, O4, S4 measured on scratch builds; trial patches kept as diffs in `scratch/logs/phase4_trial_patches/`. |
| 5, the report | the plan, this report, the index | The plan moved to `completed/`; `doc/refactoring_notes/README.md` updated. |
| 6, finishing (ruled after the report) | 72 (`review/phase_6.diff`; list in `scratch/logs/phase6_files.txt`) | O1 and O4 (CTF fit) adopted; X1 and G3 finished in bounds; the particle-row slot stop with its producers (`frame` became slot 56, other keys stripped); the sieve fix (Option 1) and the prefix stop; the stream UI contract regenerated; the frozen-set schema raised; the three known issues fixed (stale Cartesian slot, atomic rename, CTF memo sampling); guides, policies and this report. |

## Phase 6

The maintainer accepted every recommendation of the report ("I accept everything. Let's finish this up.") and
ruled that nothing in bounds stays deferred. The sieve came into bounds; the items that need `ms_posterior2`'s
files go to the follow-up run `agent_friendly2`. The ruling is recorded in the plan's sections 8.1 and 8.3.

What was learned:
- **The new stops caught nothing in the existing workflows**, and the two the code itself needed were already
  ruled: O1 for icosahedral runs, and the slot producers.
  - Without O4, the phase's tree reproduces phase 3's workflow outputs exactly.
  - The Debug build, where the slot stop is active, runs every unit suite and the `library` label without
    tripping it.
  - So the X1 and G3 completions, the slot stripping, the atomic rename, the CTF memo check and the sieve fix
    change no result, only what happens when a rule is broken.
- **The icosahedral axis search is not exact.** After O1 its aligned map scores 0.91 to 0.98 under the builder's
  operators, against 0.495 in the old setting, and the value changes between runs. The O1 test therefore asks for
  the suite's control ceiling, not its exact-target floor.
- **One stop caught a tester's own design.** The first version of the stale-slot test passed on the old code,
  because evaluating an invalid slot already stopped there. Checking each new test on the previous tree found it.

Evidence: fast gate 15/15 in Release and Debug; `test=units` 152 groups, 14,910 assertions, all passed in both;
the `library` label 7/7 in Release (in Debug, with the slot stop active, it tripped nowhere). The new checks fail on the phase-5 tree: 13 assertions over 6 suites. The
workflow entries differ from phase 3 only through O4: whole-volume correlation 1jxy 0.645, 6vxx 0.509; the
defocus errors move by 0.1 nm or less. Details and log paths: the plan's Phase 6 Progress entry.

## Decisions taken under delegation

Each is recorded in the plan's Progress entry of its phase.

- **Phase 1:** 86 declared fields that `parameters` derives are marked "not a command-line key" and left out of
  the accepted keys, rather than bound. Seven accepted-but-unbound keys that programs or the user interface name
  are bound. No parameter was added.
- **Phase 2:** C5 is narrowed to "a partition worker (not `stream=yes`) may not write its input project". Stream
  jobs own their job projects.
- **Phase 2:** stops that can fire inside OpenMP parallel regions (the polar calculator's memo and weights checks,
  `polarize*`, `apply_ctf`) are kept as fatal stops for programming errors, with fixed messages.
- **Phase 2, from the diff review:** I5 was extended to the polar calculator and the class sums; `esig%kill` bumps
  the spectra counter; `score_project_cavgs` also restores the worker identity.
- **Phase 2:** test children that write a project first set one OpenMP thread, so a regression fails the suite
  instead of hanging it.
- **Phase 3:**
  - N10 takes the investigator's design A: the service applies the request's settings and restores the caller's.
    Design B (consumers read the request) would need the out-of-bounds refine3D strategy.
  - The renames stay within the five folders the file table allows.
  - V2's trailing-chain manifest version was raised, so old chains are re-seeded.
- **Phase 4:**
  - I1 judged by solve2D over nine seeds plus the workflow entries;
  - O1 judged by a symmetrized two-blob map and a `c1` against `pgrp=i` reconstruction;
  - O4 judged in three variants for the CTF fit;
  - S4 on synthetic beam-shift patterns, because no data on the machine has per-micrograph beam shifts.
- **Phase 6:**
  - O1 deletes RELION's icosahedral table and the `icorelion` argument, which the change made unused.
  - O4 states the shift search's strict line search at its calls.
  - The slot stop is compiled in Debug builds only (`_DEBUG`).
  - The slot producers follow the ruling: `frame` is read back after a write, so it became a slot. The CTF and
    optics constants are not read back, so they are stripped where they are copied in.
  - G3 changes allocation sizes only, not the parallel regions' code; `roavg` keeps its own forced team.
  - `simple_rename` keeps an ignored `overwrite` argument for one out-of-bounds caller.
  - The O1 test asks for the control ceiling, after four phantom geometries.

## Deviations from the plan

- **Section 8 deferrals and the phase-6 remainder go to `agent_friendly2`** (on `ms_posterior2`'s committed tree):
  - R4/X4 (the `internal:` parameters), C3, C4, the folder moves of 3.4 (R6) and the judging runs of 3.7;
  - the 3.1 check that workers verify the master's sample (the refine3D master);
  - `check_pass_team` and the probability tables' `nthr_glob` sizing (`probabilistic/`);
  - X2's commander import in `simple_refine3D_strategy.f90`;
  - the `overwrite` argument that `simple_commanders_solve3D.f90` still passes to `simple_rename`.
- **3.5 renames:** `simple_cmdline.f90`'s `set_1..6` remain (not in the file table).
- **Weaker checks:**
  - O10's check does not fail on the base, because the child stops on its missing movie there.
  - O6 has no tester check, since the optimiser is a private pointer.
  - The p05 sieve restart case passes on the old tree by construction (it plants the lists).
  - The GPU-offload reconstructor's X1 pass is not compiled on this machine (`USE_OPENMP_OFFLOAD`).
- **File table:** rows were added when a phase needed a file (listed in each Progress entry).

## Evidence against the baseline

- **Gates and suites:**
  - fast gate 15/15 in Release and Debug at the end of every phase;
  - `simple_test_exec test=units` all passed in both builds: 152 groups, 14,910 assertions in phase 6;
  - the `library` label 7/7 in phases 1-3 and 7/7 in phase 6;
  - `check_guide_imports`, `check_test_registry`, `check_descr` and `check_cli_chain` pass.
- **Workflow entries** (`SIMPLE_SEED=20260923`, compared with `test=identity`):

| Phase | Result |
| --- | --- |
| 1 | identical to the base, apart from paths and run identifiers |
| 2 | X5 moves the simulated input by half a voxel, so every output differs. Without X5 the tree is identical to the base except the ruled keys `mcconvention` (O9) and `thumb_den` (O14). Final numbers: 1jxy now completes (correlation 0.641; over two parts 0.638, FSC 35.7 Å); 6vxx correlation 0.519, FSC 49.9 Å (base 0.561, 20.8 Å; both below the 0.80 floor); mini_stream passes |
| 3 | identical to phase 2 except the 8 density files, which differ only by V2's 64-byte header |
| 6 | O4 changes the CTF fits, so everything downstream differs: 1jxy 0.645 (over two parts 0.631, FSC 49.9 Å), 6vxx 0.509; without O4 the tree is identical to phase 3 |

- **New checks on the previous code:**
  - phase 2: 134 assertions in 13 suites fail on the base;
  - phase 3: 8 fail on the phase-2 tree;
  - phase 6: 13 fail on the phase-5 tree.

## Measurements and recommendations (phase 4), and the ruling

| Item | Finding | Ruled 2026-11-05 |
| --- | --- | --- |
| I1, class averages deapodized with the stencil envelope | Amplitudes change as predicted (up to 1.9× at the edge). solve2D over 8 seeds: no gain (median best class 8.4 → 10.6 Å). Workflow entries moved both ways | not adopted; waits for the judging runs of 3.7 |
| O1, one icosahedral setting | Maps aligned by `symaxis_search`/`symmetrize_map` were in RELION's setting, refinement uses SPIDER's. FSC between `c1` and `pgrp=i` reconstructions from exact poses: 0.32-0.70 before, 0.93-0.99 with everything on SPIDER's setting | adopted in phase 6, with a `pgrp=i` library test |
| O4, L-BFGS-B stopping tests | CTF fit: standard line search with default stopping keeps the accuracy, about 10% faster; a literal translation doubles the defocus error. Shift search: the standard line search fails | adopted for the CTF fit in phase 6; the shift search keeps its strict line search, stated at the call |
| S4, batch optics grouping as the stream's | No beam-shift data on the machine. Synthetic: identical on well-separated shots; single linkage chains close positions, the batch clustering fragments jitter | not adopted; waits for 3.7 and a real EPU session |

## Left open for the maintainer

- **The follow-up run `agent_friendly2`:** the items listed under deviations.
- **Recorded, not changed:**
  - the icosahedral `symaxis_search` is neither exact nor deterministic between runs;
  - the motion correction's shift search (`ftexp_shsrch`, `1e-4`) and the quality model (`1e-7`) pass
    unmeasured line-search constants; the opt guide names them.
- **The `fortls` plugin** is recommended as a user-level install; it is kept in `scratch/keep/fortls_plugin/`.
- **Workflow baselines:** 1jxy completes (since X5) and the numbers moved with X5 and O4. Whether that is better
  or worse needs the judging runs of 3.7.

## Commit message

```
Agent-friendly code: checks for silent rules, CLI chain gate, renames

Run agent_friendly of doc/refactoring_notes/completed/agent_friendly_code.md:
- scripts/check_cli_chain.py in the fast gate, and the command-line fixes it found (R1-R3, P2, P3, P6, X3);
- the small fixes of the 2026-11-04 rulings and the rules of 3.1 as stops, each with a tester check that fails
  on the code before it, including the particle-row slot stop (Debug) and the sieve's prefix stop;
- X1, X2, C1, V2 (a header on gridding density files), N10, O11, V5, G3;
- O1 (one icosahedral setting, SPIDER's) and O4 (the CTF fit's line search), as measured;
- the particle sieve keeps its chunk lists, so a restart repairs chunked_mics.txt for the whole run;
- 380 numbered generic specifics renamed by what distinguishes them;
- doc/glossary.md; folder guides and policies updated; measurements of I1 and S4 reported.

See doc/refactoring_notes/completed/agent_friendly_code_report.md.
```
