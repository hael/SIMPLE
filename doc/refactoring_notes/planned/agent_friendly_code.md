# Agent-friendly code

Status, 2026-11-05: delivered with this note: section 2 and the defect fixes of section 5 (1-7).
Every ruling is decided. The batch alignment pipeline has landed (`ea4c1c20a`), so the rest of
section 7 runs as the autonomous run `agent_friendly` on the dell (section 8), except the items that
change files of the run `ms_posterior2` (on dell2 at the same time), which wait for it.

## 1. Why

A coding agent finds code by text search, reads files in slices, keeps nothing between sessions,
and learns that it was wrong only from the compiler, a `THROW_HARD` or a test. Four properties
decide how well it works in a codebase, and they help developers for the same reasons:

1. Locality: what a change touches can be found by search and read in a few files.
2. Explicit contracts: what a routine depends on, and what must happen before it, is visible at
   the call or enforced by the code.
3. Loud failure: a mistake stops the build, the run or a test within seconds.
4. One name per concept, one way to do a thing.

SIMPLE already has much of this: no global singletons (`params_glob` and `build_glob` have no
references), testers next to their code, a fast gate of about 30 s and unit suites of seconds,
static consistency checks in the gate (test registry, `!@descr`, flex DAG), the short-comment and
single-code-path rules, no backwards-compatibility shims, and domain folders that mostly do not
import `parameters` or `builder` (none in `ori`, `opt`, `interp`, `cftc`). The items below address
what remains. They are small steps that neither depend on nor block the
[domain-driven design proposal](next_release/staged_domain_driven_design_refactor.md), whose
value objects (its section 6) are the long-term answer to unit and key mistakes.

## 2. Delivered with this note: folder guides that load themselves

- The `src/main` folder guides (`.github/skills/simple-main-*`) lead with the folder's rules,
  silent ones first, then its call order, what it does not own, its couplings and its tests.
  New guides: `cftc`, `flex`, `sieve`, `sigma2`, and `strategies` split by subfolder;
  `preprocess`, `image_processing` and `quality` are folded into the guides of their owners.
  `apis` has no import (its guide is a map), and `src/main` only points to `simple-main-root`,
  since a guide there would load on every read below it.
- A `CLAUDE.md` in the folder imports the guide, so Claude Code loads it the first time a file in
  that folder or below is read; `src/CLAUDE.md` and `production/CLAUDE.md` import
  `simple-modern-fortran`. Copilot reads the skills as before. Edit the `SKILL.md`, never the
  `CLAUDE.md`.
- `scripts/check_guide_imports.py` runs in the fast gate: Claude Code drops a missing import
  without notice.
- Basis: in the four sessions on record before this change (2026-11-02 and 03), nine subsystems
  with a guide were edited and none had its guide loaded before the first edit.

## 3. Proposed items

Ordered by payoff for the effort.

### 3.1 Silent rules become checks

Each rule a guide marks as silent is a place where neither an agent nor a developer gets feedback.
Where a check is cheap, add it and cut the guide entry to the rule and the name of the check.
Forms: a generation counter for call-order rules (the setter bumps it, the preparation records it,
the user compares once per call, never in an inner loop); a stop instead of a default; a debug-build
assertion where the check costs time. Each converted rule gets a tester that commits the violation
and expects the stop.

Seed list from the guide audit:

| Rule (guide) | Silent today because | Proposed check |
| --- | --- | --- |
| A particle row carries only its named slots (ori, project) | the particle writer stores `ori2prec` alone, so any other key is dropped at the next project write | a debug-build stop when the writer meets a non-slot key |
| A new particle batch is memoised before scoring (pftc) | `reallocate_ptcls` keeps the old transforms when the batch size is unchanged | particle generation counter, compared by the scoring entry points |
| The sigma2 weights are refreshed after the spectra change (pftc, sigma2) | the self-term is rebuilt only by `memoize_sigma2_weights` | spectra generation counter |
| `prep_sigmas_objfun` runs after every `pftc%new` (search) | `new` drops the spectra pointer and the weights stay unset | a Euclidean score without weights stops |
| A pass runs with `omp_get_max_threads() == nthr_glob` (strategies) | a larger team writes out of bounds, unchecked in Release | one check at the start of a pass |
| The master writes its sample before workers are scheduled (parallelization) | workers reproduce the previous iteration's sample | workers compare the sample's iteration with `which_iter` |
| `finalise_ptcl_import` runs once, after import and copy (star) | skipped or repeated, shifts are wrong | a finalised flag; a second finalise or an export of an unfinalised import stops |
| The records marked included stay a prefix of the list (sieve) | otherwise a record is chunked again | an assertion where the list is written |
| `cleanup_nu_filter` runs between states (nu-filt) | setup clears neither the solvent clamp nor the null shell, so the next state uses the previous clamp | setup stops when a clamp is still set |
| A gridding partial or chain file matches the reader's box (volume) | the files have no header, so a larger file is read as a prefix | a box and sampling header, checked by `read_rho` |
| CTF application on an image follows `memoize_ft_maps` for its box (image) | nothing records which box is memoised; `ctfops` skips it (section 5) | the memo records its box; `apply_ctf` stops on a mismatch |
| `polarize*` matches the memoised polar geometry (image, pftc) | box, `pdim` and mode are not compared | the memo records its geometry; `polarize*` stops on a mismatch |
| A model's `feature_weights` and `linear_coefficients` follow the feature bank's order (cavg_quality) | a short vector is zero-filled | stop on a length mismatch |

The batch alignment pipeline refactoring
([plan](../../implementation_notes/completed/prob_align_batch_pipeline.md)) rewrites
`strategies/search`; the rows for that folder wait for it.

### 3.2 A static check of the command-line chain

The params and ui audit found keys that fail without error, or fail at a user's first try. A static
check over the source tree catches all of them, in the pattern of `check_test_registry.py`:

- declared in `parameters` but never bound by the parser: `mul` (read by `oriops` and `volops`,
  so always 1.0) and `import_type` (always `auto`, so the destination checks of
  `import_starproject` never run);
- UI inputs that are not accepted keys, so giving one stops the program: `positive`, `heldout`,
  `rmsd_file`, `stk_traj`;
- the accepted-key generator (`scripts/simple_args_generator.pl`) reads past the type, because
  `contains` is not indented: `cntfile`, `silent`, `final_planned` and the fragment `1.0]` are
  accepted keys;
- private programs in neither the public UI nor `simple_private_prgs.f90` get no required-key
  check, and that table's count is hard-coded;
- a key bound twice: `rec_states` (`simple_parameters_parse.f90`).

Proposal: `scripts/check_cli_chain.py` in the fast gate (every UI input is an accepted key, every
accepted key is bound, every private program has a table entry) and a generator that stops at
`end type`. Fixing today's violations needs the rulings of section 5.

### 3.3 Orientation keys

Particle parameters have fixed slots with named indices (`I_*` in `simple_defs_ori`). Any other
key goes to a hash: `get` of a missing key returns 0.0 (`simple_hash.f90`, `get`) and `set` of a
misspelt key adds a new one (`simple_ori.f90`, `set_1`). About 8,000 calls of the get/set family
name their key as a string literal (630 distinct keys across orientations, command lines and
reports); 671 guard with `isthere`. The particle slots `class`, `eo`, `corr` and others count as
present only when non-zero. A scan of the 630 keys for pairs one edit apart found no live
misspelling.

Ruled (section 5, 12): a missing key keeps returning 0. What remains is the presence rule in the
ori guide, the fixed `get_neven` (section 5, 7) and, optionally, a registry of known literal keys
checked statically as in 3.2, which catches a misspelt key at the gate without changing `get`.

### 3.4 One folder per owner

Guides attach to folders, so a folder that mixes owners gets no guide or one loaded for unrelated
work. `src/main` holds sixteen loose files: builder, cmdline, convergence, the final
reconstruction, the refine3D stage plan, startup, states controller and workflow controller,
external-reference pose initialization, the class-average work project and the simulator, with
their testers. A guide in `src/main` would load on every read below it. Proposal: move them into
named folders, for example `core/` (builder, cmdline) and `refine3D/` (the refine3D modules and
the external-reference initialization), with the rest placed by owner decision. Cost: `git mv`,
a CMake re-configure (`src/CMakeLists.txt` collects `main/*.f90` recursively), guide and doc path
updates; no aliases (release 4).

### 3.5 Descriptive names for numbered specifics

About 100 generic bindings dispatch to specifics numbered `_1`, `_2` (`mirror_ref_pft_1/_2`,
`gen_corr_for_rot_8_1/_2`). At a call site the body that runs depends on the argument types, and a
search lands on the binding. Rename each specific by what distinguishes it (`_by_index`,
`_by_pft`) and keep the generic name. Mechanical, one folder per change; the compiler checks it.

### 3.6 A glossary

One concept under several names costs more than abbreviations do. The even/odd half is `iseven`
(200 uses), `ieo` (146) and the key `'eo'` (98); state 0 means rejected; `_8` marks
double-precision variants and `l_` logical flags; `box`/`box_crop`, `smpd`/`smpd_crop`,
`kfromto`, `trs`, `nspace`, the `ptcl2D`/`ptcl3D` and `cls2D`/`cls3D` segments,
`sampled`/`updatecnt`. Proposal: one page, `doc/glossary.md`, linked from `simple-architecture`
and CLAUDE.md; new code uses its terms. It records today's code names and does not replace the
domain language of the DDD proposal (its section 4).

### 3.7 A judging run per workflow

A change meant to improve results needs a run that can show it. `simulated_workflow` misses the
map in 1 to 2 of 20 runs (profile of 2026-10-19), so one run shows breakage, not improvement.
Proposal: for refine3D, solve2D, solve3D and the stream, one fixed-seed, fixed-data run of a few
minutes with its metric, its noise band and how to read it, in
`doc/policies/test_environment_policy.md`. `test=profile` and `test=identity` provide the
machinery.

### 3.8 A Fortran language server (experiment)

Claude Code connects language servers through plugins (an `.lsp.json` naming the server and the
extensions): diagnostics after each edit, and a read-only `LSP` tool for definitions and
references by symbol, which plain search cannot do for the 254 bindings named
`new`/`kill`/`get`/`set`/`execute`. There is no official Fortran plugin; the candidate is
`fortls` (pip), not installed on the dell. Try it in user space with a local plugin for `.f90`
and measure: resolution of type-bound procedures into submodules, indexing time and memory on
SIMPLE, and false diagnostics from `THROW_*` macros and generated sources. Drop it if it does not
cope.

## 4. Not proposed

- Breaking up `parameters` (662 fields) or the `(params, build)` signature (134 routines). It
  hides what a routine reads, but it is the layering policy and a rewrite would not pay. New
  domain code stays on narrow inputs (the exception today: ten `volume` modules use `parameters`).
- Splitting files by size (sixteen over 2,000 lines). An agent reads in slices; cohesion matters
  more than length.
- Renaming idiomatic `new`/`kill`/`get`/`set`. The lifecycle symmetry is worth more than easier
  search; 3.8 addresses search.

## 5. Rulings

Defects found by the guide audit, decided 2026-11-04:

1. Sieve: `cleanup_chunk` deletes `chunk_mics.txt` (and `consumed_coarse.txt`) after a chunk's
   rejection, so `rebuild_chunked_mics` gives up once a coarse chunk has been rejected and the
   crash repair of the [sieve policy](../../policies/sieving_and_rejection/ptcl_sieve_policy.md)
   section 4.7 no longer holds (a crash in the window duplicates a coarse chunk's particles).
   Confirmed and reproduced; a report went to the sieve's developer, who owns the fix.
2. STAR: RELION's `rlnCoordinateX/Y` is the particle centre, SIMPLE's `xpos/ypos` the corner of its
   box (the EMAN convention). Fixed: `finalise_ptcl_import` takes half the stack's box off the
   imported coordinates once (`coords_to_corners`); `star_project` tests the round trip.
3. `ctfops`: `apply_ctf_imgfile` now memoises the CTF maps at twice the box, as
   `apply_ctf_wpad` pads a real-space image, and forgets them after.
4. Optional arguments are present at some calls and absent at others (now a rule in
   `simple-modern-fortran`): `ft_expanded%shift` takes `self_out` as required, and
   `apply_ctf_imgfile` lost its never-passed `bfac`.
5. `atoms%new_from_pdb` reads the record name, so HETATM survives a read and a write (`atoms`
   tester).
6. Pick: `picksegdiam` converts its default box to pixels of the micrograph
   (`BOX_EXP_FAC*moldiam_max/smpd`), and the multi-diameter Gaussian references of
   `simple_picker_utils` are sized in pixels like the other branches.
7. `get_neven` counts the even rows of a particle field once the field has odd rows, so
   `get_nevenodd() == 0` still means "no halves yet" (`oris` tester).

Decided 2026-11-04 (the work is in section 7):

8. `mul` and `import_type` (3.2): bind them.
9. `positive`, `heldout`, `rmsd_file`, `stk_traj` (3.2): bind each one a program reads, drop the rest.
10. The key generator (3.2): stop at `end type`.
11. The eight parameters documented as `internal:` (`unfil_unmasked`, `cc_emit_sigma`,
   `cohort_sampling`, `scan_exit`, `prob_post_diag`, `cont_route`, `sigma_commit_deferred`,
   `lpexit`) and the five `INPROCESS_KEYS` of `simple_cmdline.f90`: retire them; add none.
12. 3.3: a missing key keeps returning 0.
13. 3.4: `src/main/core/` (builder, cmdline) and `src/main/refine3D/` (the refine3D modules and the
   external-reference initialization); the other loose modules placed when moved.
14. 3.7: `simulated_workflow` 6vxx and 1jxy at fixed seeds for solve2D, solve3D and refine3D, judged
   by map correlation and FSC against the truth inside the band of the 2026-10-19 profile;
   `mini_stream_6vxx` for the stream.
15. 3.8: `fortls` 3.2.2 is installed in user space on the dell (`~/.local/bin`).

The open questions of the guide audit (69 of 70) were answered the same day; the guides state the
answers, and the code work they rule is in section 7.

## 6. Order

1. The small items of section 7, then 3.2 (one script; its first findings are already fixed by
   then). 3.1: one rule at a time, outside `strategies/search` until the batch pipeline lands.
2. 3.8 (an afternoon; `fortls` is installed) and 3.6 (a day).
3. The medium items of section 7, with 3.4 (the folders are named) and 3.5, folder by folder.
4. 3.7 (data and metrics chosen), then the large items of section 7, each with its measurement.
5. 3.3's optional key registry, if 3.2 shows it is worth it.

## 7. Work from the guide rulings

The code changes the 2026-11-04 rulings call for, by size. Each comes with a tester check where
one fits; rows that touch `strategies/search` or the probability tables wait for the batch pipeline.
Ids are those of the rulings (proposal R, guide questions X, P, C, G, S, V, I, O, N).

Small (a few lines in one routine):

- Command line: bind `mul` and `import_type` (R1); per UI input, bind or drop (R2); the generator
  stops at `end type` (R3); validate yes/no keys (P6); every private program in
  `simple_private_prgs.f90`, count computed (P3); `print_ui_stream` under the stream contract (P2);
  drop the stream aliases `initial_analysis` and `opening2D` (X3).
- Project and layering: stop project writes in a distributed worker (C5).
- Orientations: `rot_to_asym` initialises its local from the record (O2); `eo = -1` stops the half
  sums instead of joining the even half (X1); `center_pdbcoord` centres at ldim/2 (X5).
- Optimisers: the clients deallocate their optimiser (O6); an all-zero start is a start, with an
  explicit "no start" (O3); differential evolution uses `npop` and `set_inipop` (O5).
- CTF and motion: `ctf_estimate_fit%kill` keeps the process-wide FT maps (O7); stop
  `extractfrommov=yes` unless `mcconvention=simple` (O9); stop `nframesgrp` outside time series
  (O10).
- Picking: the 8 A erosion also when picks are assigned to bins (O13); `thumb_den` on the segdiam
  route (O14).
- Images: the polarizer records its geometry and `polarize*` checks it (I2); `shconst` on odd boxes
  is 2 pi/n (I5); `cartft_calc%set_sigma2` invalidates prepared slots (I3).
- STAR: `export_cls2D` works on a copy (S2); `export_cls2D`, `export_stardata` and
  `relion_project` write a temporary file and rename it (S3).
- Nano and flex: stop on an element missing from `get_lattice_params` (N3); restore
  `omp_max_active_levels` (N9); graphene subtraction and `score_project_cavgs` take the caller's
  threads and seed (N1, N6).
- Comments and docs: the gridding FSC comments (V1); the routers' extra duties into layering
  policy section 8 (P1); the flex gateway comment (N8).

Medium (a module or a flow):

- `refine_motion_model` reconstructs through the strategy routine, not `commander_rec3D` (X2).
- `read_segment(..., fromto)` honours the range for `.simple` files (C1).
- One builder of the class-average row layout, used by `solve3D_cavgs` too (C3).
- The run state of `simple_solve3D_utils` becomes a type (C4).
- A box and sampling header on gridding partials and chains, checked by `read_rho` (V2).
- `rec3D_request` carries the reconstruction settings `reconstruct_state_halves` now overrides (N10).
- Moves: `pspecs` to `strategies/search` (O11), `simple_dose_series_manifest` to `motion/` (V5), the
  loose `src/main` modules to `core/` and `refine3D/` (R6, 3.4).
- Per-thread search state sized by the team (G3), after the batch pipeline.

Large, or measured before adopting:

- Retire the eight `internal:` parameters and the five `INPROCESS_KEYS` (R4, X4): each needs another
  channel from controller to child (a typed argument of the strategy routine, or project state).
- 2D class averages deapodize with the stencil envelope (I1): compare a solve2D run with the
  current build.
- One icosahedral setting (O1): check the `symaxis_search` to refine3D path with an icosahedral test
  first.
- L-BFGS-B clients set `factr`/`pgtol` (O4): compare shift-search and CTF-fit results before and
  after.
- Batch optics grouping moves to the stream's (S4): compare the groups on one data set with beam tilt.
- The judging runs of 3.7 and the language-server experiment of 3.8.

Deferred until the run `ms_posterior2` has landed, since they change its files
(`doc/implementation_notes/planned/multistate_posterior_state_draw.md`, "File table"): retiring the
`internal:` parameters and the `INPROCESS_KEYS` (R4, X4), one class-average row builder (C3), the run
state of `simple_solve3D_utils` as a type (C4), the folder moves of 3.4 (R6) and the judging runs of
3.7 (the test registry).

## 8. The autonomous run `agent_friendly`

### 8.1 Scope

The run carries out section 7 without the deferred items: the small items, the medium items, the
measurements of the large ones (measured and reported, never adopted), the checks of 3.1 and 3.2,
the renames of 3.5, the glossary of 3.6 and the language-server experiment of 3.8. The rulings are
section 5 and the ids of section 7; the folder guides already state each guide question's ruling.
Two kinds of file are out of bounds whatever a row of the file table seems to allow: the 20 files
of `ms_posterior2`'s file table (among them `simple_refine3D_strategy.f90`,
`simple_prob_align_strategy.f90`, `simple_convergence.f90`, `simple_solve3D_controller.f90`,
`simple_commanders_solve3D.f90`, `simple_commanders_prob.f90`, the files of
`strategies/search/probabilistic/`, `simple_test_registry.f90` and the policies `refine3D_policy.md`,
`solve3D_policy.md`, `refine3D_states_policy.md` and `test_environment_policy.md`), and
`src/main/sieve/`, whose recovery defect its developer is fixing. If an item cannot be done
without one of them, it is recorded as deferred with the reason, and the run goes on.

### 8.2 The machine and how the work is judged

The dell (24 cores: 8 performance, 16 efficiency; 62 GB; GCC 15.2; CMake 3.26; Python 3.12 at
`/usr/bin/python3.12`, `python3` is 3.6). The run builds and tests only in its own copy and its own
build trees. The maintainer works on the machine meanwhile; record `uptime` with every timing.

- Every change comes with the check that shows it: a tester assertion that fails on the base and
  passes now where the change is testable in a unit, otherwise the program or workflow run that
  exercises it, named in Progress.
- Each phase ends with the fast gate (15/15) in a Release and a Debug build, the unit suites of the
  touched areas, and `scripts/check_guide_imports.py`, `check_test_registry.py`, `check_descr.py`.
- Behaviour that should not change does not: phases 2 and 3 run `simulated_workflow` (1jxy and
  6vxx) and `mini_stream_1jxy` with `SIMPLE_SEED` set, on the phase's build and the previous phase's,
  and compare them with `test=identity`; a difference is either explained by a ruled change or a
  stop. Pure renames and moves (3.5, O11, V5) must be byte-identical.
- No floor, tolerance or assertion is relaxed. Failures the base already has are recorded in
  phase 0 and not fixed unless their cause is in this run's files.

### 8.3 Phases

- **Phase 0, baseline and guides.** Build the base (Release and Debug) and record the fast gate, the
  unit suites and the base failures of the workflow entries above. Bring the guides up to date with
  the landed batch pipeline: every statement in the strategies guides (umbrella, `-parallelization`,
  `-search`, `-probabilistic`), pftc, sigma2, class, cftc, commanders, exec, apis, root, volume and
  their `references/` about the probability tables, part files, table programs, particle batches or
  `prob_align` is checked against HEAD and rewritten, and every link to a note that moved to
  `completed/` is fixed. Exit: both builds pass the gate; Progress records the baselines and the
  guide changes.
- **Phase 1, the command-line chain.** R1, R2, R3, P2, P3, P6, X3 and the gate check of 3.2
  (`scripts/check_cli_chain.py`: every UI input an accepted key, every accepted key bound, every
  private program in `simple_private_prgs.f90`, no key bound twice), run by `run_fast_gate.sh`.
  Exit: the check passes and fails on a planted violation of each kind (recorded, then removed).
- **Phase 2, the small fixes and the checks of 3.1.** The small items of section 7 not in phase 1,
  and the rows of 3.1 whose files are in bounds, each with its check. Exit: every item done or
  recorded as deferred with the reason; the identity comparison of 8.2 against phase 1.
- **Phase 3, the medium items, renames and the glossary.** X2, C1, V2, N10, O11, V5, G3, the renames
  of 3.5, the glossary of 3.6 (`doc/glossary.md`, linked from `simple-architecture`), and the
  experiment of 3.8 (a local plugin for `fortls` 3.2.2 in `~/.local/bin`, kept under the run's
  scratch; its findings in Progress, nothing of it in the repository). Exit: as phase 2, renames
  and moves byte-identical.
- **Phase 4, measurements.** I1, O1, O4 and S4, each measured as section 7 says on scratch builds
  and data, with the numbers and a recommendation in Progress. No repository code changes.
- **Phase 5, the report.** `doc/refactoring_notes/completed/agent_friendly_code_report.md`: what was
  learned, the evidence, the measurements and their recommendations, the deferred items with their
  reasons, and a commit message; the refactoring notes index updated; this note moved to
  `doc/refactoring_notes/completed/` with `mv`.

### 8.4 Common rules

- No new parameter, except the fields ruling R2 binds for an input a program already reads.
- Unreferenced is not dead: delete only what a ruling deletes (X3) or what a change itself makes
  unused.
- A change that alters a rule a folder guide states updates the guide's `SKILL.md` in the same
  phase; the importing `CLAUDE.md` files are never edited.
- New tests go into existing suites (a new sub-suite changes the test policy, which is out of
  bounds).
- Release 4: no backwards compatibility except reading earlier `.simple` projects.

### File table

| File | What | Phases |
|---|---|---|
| `doc/refactoring_notes/planned/agent_friendly_code.md` | this plan, its Progress | each |
| `.github/skills/` | the folder and workflow guides | each |
| `params/` | R1, R2, P6, O9, O10 | 1, 2 |
| `ui/` | R2, P2, P3 | 1 |
| `exec/` | P3 | 1 |
| `apis/` | P3 | 1 |
| `production/` | X3, P3 | 1 |
| `scripts/` | R3, the gate check of 3.2, run_fast_gate.sh | 1 |
| `project/` | C5, C1 | 2, 3 |
| `ori/` | O2, X1, renames | 2, 3 |
| `opt/` | O3, O5, O6 | 2 |
| `ctf/` | O6, O7 | 2 |
| `motion/` | O10, V5 (the dose-series manifest moves here) | 2, 3 |
| `preprocess/` | O9, O11 | 2, 3 |
| `pick/` | O13, O14 | 2 |
| `image/` | I2, I5, 3.1 (CTF memo), renames | 2, 3 |
| `image_processing/` | 3.1 (CTF memo) | 2 |
| `cftc/` | I3 | 2 |
| `pftc/` | O6, 3.1 (memo counters), renames | 2, 3 |
| `class/` | X1 | 2 |
| `sigma2/` | 3.1 (weights) | 2 |
| `volume/` | V2, V5, N10 | 3 |
| `nu_filt/` | 3.1 (clamp) | 2 |
| `star/` | S2, S3 | 2 |
| `nano/` | X5, N1 | 2 |
| `flex/` | N9, N10 | 2, 3 |
| `cavg_quality/` | N6, 3.1 (model length) | 2 |
| `quality/` | 3.1 | 2 |
| `src/defs/` | N3 | 2 |
| `src/fileio/` | 3.1 (particle writer), renames | 2, 3 |
| `src/utils/` | renames | 3 |
| `strategies/search/simple_*.f90` | X1 (simple_matcher_3Drec.f90), G3, the pspecs move; not probabilistic/ | 2, 3 |
| `strategies/parallelization/simple_refine_motion_model_strategy.f90` | X2 | 3 |
| `strategies/parallelization/simple_rec3D_service.f90` | X2, N10 | 3 |
| `strategies/parallelization/simple_rec3D_service_tester.f90` | X2, N10 | 3 |
| `strategies/parallelization/simple_volassemble_strategy.f90` | V2 | 3 |
| `strategies/parallelization/simple_pick_strategy.f90` | O14 | 2 |
| `commanders/test/simple_commanders_test_class.f90` | suite tables of new tester routines | 2, 3 |
| `doc/policies/layering_and_project_ownership_policy.md` | X2 (section 8) | 3 |
| `doc/glossary.md` | 3.6 | 3 |
| `doc/code_overview/fortran-indexes/` | regenerated after renames and moves | 3, 5 |
| `doc/refactoring_notes/completed/agent_friendly_code_report.md` | the report | 5 |
| `doc/refactoring_notes/completed/agent_friendly_code.md` | this plan, moved | 5 |
| `doc/refactoring_notes/README.md` | the index | 5 |

## Progress

