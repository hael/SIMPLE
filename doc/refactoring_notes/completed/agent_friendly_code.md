# Agent-friendly code

Status, 2026-11-05: delivered with this note: section 2 and the defect fixes of section 5 (1-7). The run
`agent_friendly` (section 8) carried out section 7 on the dell on 2026-11-05, uncommitted for the maintainer's
review; its record is the Progress section below and its summary the
[report](agent_friendly_code_report.md). Phase 6 (ruled after the report) adopted O1 and O4 and finished the
in-bounds remainder. Left for the follow-up run `agent_friendly2` on `ms_posterior2`'s tree, since they change its
files: R4/X4, C3, C4, the moves of 3.4 (R6), the judging runs of 3.7 and the phase-6 remainder listed in its
Progress entry.

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
[domain-driven design proposal](../planned/next_release/staged_domain_driven_design_refactor.md), whose
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
without one of them, it is recorded as deferred with the reason, and the run goes on. Ruled 2026-11-05 for
phase 6: `src/main/sieve/` is in bounds (this run's fix supersedes the developer report's), and the items that need
one of `ms_posterior2`'s files go to the follow-up run `agent_friendly2`, on `ms_posterior2`'s committed tree.

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

- **Phase 6, finishing (ruled 2026-11-05).** The maintainer accepted every recommendation of the
  report. Within bounds (the 20 files of `ms_posterior2` stay out; their items go to the follow-up
  run `agent_friendly2` on `ms_posterior2`'s tree):
  - O1: one icosahedral setting, SPIDER's, everywhere (the symmetry search, `symmetrize_map`, the
    symmetry analysis), as the phase-4 trial measured; and O4 for the CTF fit only: the standard
    line search with the default stopping tests (the phase-4 variant that kept the accuracy).
    I1 and S4 are not adopted; they wait for the judging runs of 3.7.
  - X1 finished in the iterative (PCG) reconstruction backend, `flex/`, the GPU-offload
    reconstructor and the batch scoring: an `eo=-1` row stops, as in the class sums.
  - G3 finished for the polar calculator's buffers (`pftc/`) and the other in-bounds helpers.
  - The particle-row slot stop of 3.1: for each of the four producers that write other keys, a key
    that is read back after a project write becomes a slot (`references/particle_slots.md` of the
    ori guide), any other key is stripped before the write; then the stop is added.
  - The sieve (now in bounds): `chunk_mics.txt` and `consumed_coarse.txt` on `cleanup_chunk`'s keep
    list, the rule of 3.1 that included records stay a prefix of the list as a stop, both with the
    tester cases the report for the sieve's developer describes, the sieve and restart policies
    updated.
  - `production/stream_ui_contract.json` regenerated so the launcher contract check compares a
    current file; the frozen-set schema version raised so pre-header sets are refused cleanly.
  - The known issues of the report: an invalidated Cartesian slot stops instead of being dropped as a
    degenerate particle; `simple_rename(..., overwrite=.true.)` replaces the target atomically;
    `ft_maps_match` compares the sampling as well as the box.
  - The folder guides and the report brought up to date; the report gains a section for phase 6.
  Exit: every item done with its check, or recorded with the evidence of why it cannot be done
  without an out-of-bounds file; the gates and identity comparison of 8.2 against phase 3's tree,
  differences explained by O1, O4 or a fix of this phase.

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
| `params/` | R1, R2, P6, O9, O10; C5 (the worker's project, recorded at `params%new`) | 1, 2 |
| `ui/` | R2, P2, P3 | 1 |
| `exec/` | P3 | 1 |
| `apis/` | P3 | 1 |
| `production/` | X3, P3; P2 (`validate_complete_ui_json.py`, `CMakeLists.txt`, `stream_launcher_contract.json`) | 1 |
| `scripts/` | R3, the gate check of 3.2, run_fast_gate.sh; `default_audit.py` reads `add_yesno` (P6) | 1 |
| `project/` | C5, C1 | 2, 3 |
| `ori/` | O2, X1, renames | 2, 3 |
| `opt/` | O3, O5, O6 | 2 |
| `ctf/` | O6, O7 | 2 |
| `motion/` | O10, V5 (the dose-series manifest moves here) | 2, 3 |
| `preprocess/` | O9, O11 | 2, 3 |
| `pick/` | O13, O14 | 2 |
| `stream/stages/simple_stream_stage_initial_analysis_tester.f90` | O13 (the bin of an eroded pick, tested on the stage's picker) | 2 |
| `image/` | I2, I5, 3.1 (CTF memo), renames | 2, 3 |
| `image_processing/` | 3.1 (CTF memo) | 2 |
| `cftc/` | I3 | 2 |
| `pftc/` | O6, 3.1 (memo counters), renames | 2, 3 |
| `class/` | X1; G3 (if the class averager's per-thread images are sized by the team) | 2, 3 |
| `sigma2/` | 3.1 (weights) | 2 |
| `volume/` | V1, O6 (`simple_volpft_symsrch.f90` leaks its optimisers too); V2, V5, N10 | 2, 3 |
| `nu_filt/` | 3.1 (clamp) | 2 |
| `star/` | S2, S3 | 2 |
| `nano/` | X5, N1 | 2 |
| `commanders/single/single_commanders_trajectory.f90` | N1 (graphene subtraction takes the caller's parameters) | 2 |
| `flex/` | N9, N8 (the gateway comment, `README.md`), N10 | 2, 3 |
| `cavg_quality/` | N6, 3.1 (model length) | 2 |
| `quality/` | 3.1 | 2 |
| `src/defs/` | N3; C5 (`worker_projfile_glob`) | 2 |
| `src/fileio/` | 3.1 (particle writer), C1 (`simple_binoris.f90`: the string segments' row range), renames | 2, 3 |
| `src/utils/` | N6 (`simple_rnd`: the intrinsic generator's state saved and restored); renames | 2, 3 |
| `strategies/search/simple_*.f90` | X1 (simple_matcher_3Drec.f90), 3.1 (the pass's team, checked in the matchers), I2 (`polarize` callers), G3, the pspecs move; not probabilistic/ | 2, 3 |
| `strategies/parallelization/simple_refine_motion_model_strategy.f90` | X2 | 3 |
| `strategies/parallelization/simple_rec3D_service.f90` | X2, N10 | 3 |
| `strategies/parallelization/simple_rec3D_strategy.f90` | X2 (`run_rec3D`, the reconstruction's strategy entry) | 3 |
| `commanders/simple/simple_commanders_rec.f90` | X2 (`commander_rec3D` becomes a thin wrapper of `run_rec3D`) | 3 |
| `doc/policies/refine_motion_model_policy.md` | X2 (its reconstruction step names the strategy routine) | 3 |
| `strategies/parallelization/simple_rec3D_service_tester.f90` | X2, N10 | 3 |
| `strategies/parallelization/simple_volassemble_strategy.f90` | V2 | 3 |
| `strategies/parallelization/simple_pick_strategy.f90` | O14 | 2 |
| `commanders/test/simple_commanders_test_class.f90` | suite tables of new tester routines | 2, 3 |
| `doc/policies/layering_and_project_ownership_policy.md` | P1 (section 8), C5 (section 5: the worker's refused writes); X2 (section 8) | 2, 3 |
| `doc/glossary.md` | 3.6 | 3 |
| `doc/refactoring_notes/planned/next_release/x_execution_layer_adoption.md` | 3.5 (it names the renamed `generate_script` specifics) | 3 |
| `flex/fit/simple_flex_pca_embed.f90` | 3.5 (the one direct caller of a renamed specific outside the renamed folders: `locate_2`) | 3 |
| `doc/code_overview/fortran-indexes/` | regenerated after renames and moves | 3, 5 |
| `doc/refactoring_notes/completed/agent_friendly_code_report.md` | the report | 5 |
| `doc/refactoring_notes/completed/agent_friendly_code.md` | this plan, moved | 5 |
| `doc/refactoring_notes/README.md` | the index | 5 |
| `.github/skills/` | the guides (phase 6) | 6 |
| `production/` | the stream UI contract | 6 |
| `scripts/` | checks | 6 |
| `src/defs/` | | 6 |
| `src/fileio/` | | 6 |
| `src/utils/` | simple_rename | 6 |
| `volume/` | X1, O1 | 6 |
| `flex/` | X1 | 6 |
| `ori/` | O1, slots | 6 |
| `ctf/` | O4 | 6 |
| `opt/` | O4 | 6 |
| `pftc/` | G3 | 6 |
| `image/` | ft_maps_match, O1 | 6 |
| `cftc/` | the Cartesian slot | 6 |
| `class/` | X1 | 6 |
| `sieve/` | the keep list, the prefix stop | 6 |
| `stream/` | slot producers, the sieve's callers | 6 |
| `star/` | slot producers | 6 |
| `project/` | slot producers | 6 |
| `nano/` | slot producers | 6 |
| `cavg_quality/` | slot producers | 6 |
| `motion/` | slot producers | 6 |
| `preprocess/` | slot producers | 6 |
| `strategies/search/simple_*.f90` | X1 batch scoring, G3; not probabilistic | 6 |
| `commanders/test/simple_commanders_test_class.f90` | suite tables | 6 |
| `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` | the sieve's lists | 6 |
| `doc/policies/stream/restart_policy.md` | the sieve's lists | 6 |
| `doc/refactoring_notes/completed/agent_friendly_code_report.md` | the report | 6 |
| `doc/refactoring_notes/README.md` | the index | 6 |
| `solve/simple_solve3D_manifest.f90` | its `simple_rename` call loses the `overwrite` argument (replacement is always atomic) | 6 |
| `commanders/test/simple_commanders_test_highlevel_backends.f90` | the same `simple_rename` call change | 6 |
| `commanders/simple/simple_commanders_solve3D_chunks.f90` | the same `simple_rename` call change | 6 |
| `strategies/parallelization/simple_rec3D_pcg_strategy.f90` | X1 (the PCG backend's half selection stops on `eo=-1`) | 6 |
| `strategies/parallelization/simple_cls_expansion_strategy.f90` | X1 (its even/odd split stops on `eo=-1`) | 6 |
| `strategies/parallelization/simple_calc_pspec_strategy.f90` | G3 (its per-thread sums sized by the team) | 6 |
| `commanders/simple/simple_commanders_volops_tester.f90` | O1 (`symaxis_search` with `pgrp=i` leaves the map in the builder's setting) | 6 |

## Progress

### Phase 0: baseline and guides (2026-11-05, the dell)

Base: `6565f5ba5`, the tree the run started from. The run's logs are under
`/data2/agent_runs/agent_friendly/scratch` (`~/agent_runs/agent_friendly` links there), written `scratch/` below.

**How SIMPLE builds here.** GCC 15.2.1 of gcc-toolset-15 (gfortran, gcc, g++), CMake 3.26.5, no environment
modules, the system OpenBLAS with OpenMP (`OPENBLAS_CORETYPE=Haswell`), 24 cores. The maintainer's checkout
`~/src/SIMPLE` builds Debug with tests, `NICE=OFF`, `USE_ARCHOPT=ON`, installing into its build tree; the run never
touches it. The run's builds: Release in the repository's `build/` (`./compile_clean.sh`); Debug in
`scratch/build_debug` (the same CMake configuration as `compile_debug.sh`, which would overwrite `build/`, with
`-DPython3_EXECUTABLE=/usr/bin/python3.12`), then the fast gate and `make install`; and the base exported with
`git archive` to `scratch/keep/src_base`, built Release in `scratch/keep/build_base`, kept for later comparisons
(`scratch/keep/summary.tsv` lists every kept artefact).

**Builds and the fast gate** (the fast gate is the 15 unit suites `unit_<area>` that every build runs before it
installs, with the static checks `check_test_registry.py`, `check_descr.py`, `check_flex_dag.py` and
`check_guide_imports.py` first):

| build | compile (wall, `-j24`) | gate | gate wall | slowest suites |
| --- | --- | --- | --- | --- |
| Release | 1 min 44 s (load 0.08 at start) | 15/15 | 6.2 s (29.1 core-seconds) | `unit_reconstruction` 6.2 s, `unit_pftc_align2D3D` 5.0 s, `unit_cart_align3D` 3.9 s |
| Debug | 38 s (load 6.5) | 15/15 | 20.4 s | `unit_reconstruction` 20.4 s, `unit_pftc_align2D3D` 15.7 s, `unit_heterogeneity` 6.9 s |

Both gates were run again after the guide changes below, 15/15 each (`scratch/logs/phase0_final_gates.log`). Logs:
`scratch/logs/phase0_release.log`, `phase0_release_gate.log`, `phase0_debug_gate.log`, `phase0_debug_gate_ctest.log`,
`phase0_base.log`, times in `phase0_times.txt`.

**Workflow entries on the base** (Release, run as CTest runs them, with `SIMPLE_SEED=20260923`; run directories and
logs kept in `scratch/keep/phase0_base_runs/` for the identity comparisons of later phases):

| entry | result | evidence |
| --- | --- | --- |
| `mini_stream_1jxy` (six simulated micrographs through CTF estimation, picking, extraction and `solve2D`) | passes every check: defocus error at most 0.067 µm (ceiling 0.15), 21 of 48 placed particles picked (ratio 0.44, floor 0.25), 2 classes; 51 s at 8 threads | `mini_stream_1jxy.log` |
| `simulated_workflow_1jxy` (movies to `solve2D` and `solve3D`) | stops in `solve3D` stage 3, iteration 41: `new_shiftmap` "the map's last shell must lie in the band" (the polar shift map is asked for a band that its last shell falls outside); the low-pass limit stayed at 124.8 Å through stages 1-2. The base failure recorded on 2026-11-03 and in the batch-pipeline report, unchanged; 4 min 10 s | `simulated_workflow_1jxy.log` lines 2395-2404 |
| `simulated_workflow_6vxx` | exits normally (the map's quality is reported, not gated): the final map missed the truth, whole-volume correlation 0.561 against a minimum of 0.80 (the 2026-11-03 base: 0.54), docking correlation 0.876, FSC 0.143 at 20.8 Å; `solve2D` 4 classes; 12 min 16 s at one thread | `simulated_workflow_6vxx.log` lines 4854-4867 |

These two `simulated_workflow` outcomes are the base's; no item of this run addresses them, so later phases
compare against them rather than fix them. The three entries ran at the same time, on a machine otherwise idle
(load 1-2.5); they were not timed for comparison.

**Guides brought up to date with the batch pipeline** (`ea4c1c20a`: probabilistic alignment, `refine=prob` and
`prob_neigh`, now runs as a batch pipeline inside the matchers `refine3D_exec` and `refine2D_exec`, in shared memory
and in every distributed part, with every random draw keyed on the run's seed, the pass, the iteration, the
particle and the draw; the table programs `prob_tab*` and `prob_align*`, their part files, the master's global
table, its assignment step `ref_assign` and the assignment and insertion files are gone, and the master only draws
the sample, schedules once per iteration and merges the parts' small readout records). Four read-only audits
checked every statement of the guides named in phase 0, and of the workflow guides that describe the same flow,
against the code at HEAD; I made the edits and checked each name the new text uses with `rg`. A fifth read-only
review checked every added line against the code; its eight corrections (among them: the master schedules the
parts twice when the `pose_cont` polish follows; the fill-and-draw stops differ between the dense 3D, 2D and
neighbourhood tables; a part keeps more than the assignment map across batches) are applied.

- `simple-main-strategies-probabilistic`: rewritten. It described the retired design throughout; it now states the
  pipeline's order, the keyed draws and their scope, the iteration's model set once per pass, the alternation of
  fill and draw per batch, the fallbacks as the code has them now (a neighbourhood particle without a candidate is
  re-seeded at its previous orientation), the readout records and their merge, the weighted insertion's per-batch
  extras, and the right test names (`2d_probability_table`, not `2d_probability_table_io`;
  `polar_matching_bench suite=fills`; `unit_numerics suite=random_draws`).
- `simple-main-strategies`: the master's steps (no `prob_align`; the readout merge added), the 2D matcher's module
  pointers (`set_b_p_ptrs2D` no longer exists), the retired 2D table-to-buffer hand-over, the one accepted
  commander call from a strategy is now the sigma2 commit's.
- `simple-main-strategies-parallelization`: description and intro without `prob_align*_exec`; the sample draw in the
  probabilistic modes (`sample_ptcls4prob3D/2D`); a new silent rule (the parts' draws depend on the run seed and
  `which_iter` that the master hands them, so other values draw differently without error); the readout merge and
  the 2D clean-up in the masters' order; the shared-memory file reads without the old `prob_align` re-read.
- `simple-main-strategies-search`: the pass order with the pipeline's steps; the classic modes' thread-dependent
  draws against the keyed `prob*` draws; the empty partition's zero readout record; the retired table workers'
  duplicate preparation removed; the probability-table and `fills` tests.
- `simple-main-pftc`: the keyed in-plane and disk draws and their scope; the `fills` test.
- `simple-main-commanders`, `simple-main-apis`, `simple-main-ui`, `simple-architecture`,
  `simple-main-root/references/solve_controllers.md`: the deleted commanders and programs removed from lists and
  examples; the worker frame now named as the matchers'; the apis guide's out-of-date remark about the private API
  module's header (the header is correct); the private programs in neither key table named as they are now
  (`refine2D_distr`, `calc_pspec_assemble`, `aggregate_chunks`, `scale_project`); `l_prob_beta_end_override`.
- Workflow guides: `simple-refine3d` and its `bayesian-model.md` and `refine3d-code-map.md`, `simple-solve2d`,
  `simple-solve3d-importance-sampling` and its `coupling-map.md`, and `simple-frac-update-trailing`'s
  `frac-update-contract.md`: the pre-step, table workers, part files and assignment file replaced by the pipeline,
  the owner's keyed sample draw and the parts' reproduction of it.
- No finding in `simple-main-sigma2`, `-class`, `-cftc`, `-exec`, `-volume` (and its `lifecycles.md`), or in
  `simple-frac-update-trailing/SKILL.md`.
- Links to notes moved to `completed/`: every `.md` path the guides and their `references/` name was resolved
  against `git ls-files`. Two were broken and are fixed: `simple-main-nano`
  (`implementation_notes/completed/atomic_density_model.md`) and `simple-main-pftc`
  (`refactoring_notes/completed/polar_shift_search_performance_refactoring.md`). The others resolve, apart from
  `doc/code_overview/code_base_map.md` and `test_inventory.md`, which the architecture guide names as generated and
  untracked.

**Stale text found outside this phase's files** (left for the phase that owns the file, or for the maintainer):
`doc/policies/2D/solve2D_policy.md` line 382 names `sample_ptcls4iter2D`, which does not exist (the 2D draws are
`sample_ptcls4prob2D` and `sample_ptcls4update2D`); `doc/policies/3D/refine3D_policy.md` headings at lines 491 and
833 still name `ref_assign` and `ref_assign_neigh` (a file of `ms_posterior2`, out of bounds); a comment in
`strategies/search/simple_matcher_ptcl_batch.f90` line 142 names "the probability-table commanders" (phases 2-3);
the description of `polar_matching_bench` in `ui/simple_test/simple_test_ui_highlevel.f90` line 225 links
`implementation_notes/planned/polar_matching_layout_bench.md`, now in `completed/` (phase 1); comments in
`strategies/search/probabilistic/` (`simple_eul_prob_tab.f90` lines 70, 83, 1290-1292, `_neigh` line 85, `_utils`
line 1480) still speak of the global table and `ref_assign` (out of bounds: `ms_posterior2`'s files). Not
verified, for phase 2 to look at: in shared-memory `refine3D` with `volrec=no` and no polish, nothing appears to
write the particle field before the sigma2 commit's commander reads it from the file
(`publish_before_file_assembly` is false without partial reconstructions; `simple_refine3D_strategy.f90` near line
1345); whether a sigma2 update runs in that configuration was not checked.

Rulings: `spec/rulings.md` adds none beyond section 5 and the ids of section 7 for this phase.

Exit: both builds pass the fast gate (15/15, before and after the guide changes); the baselines and the guide
changes are recorded above. Files changed: this plan and 18 files under `.github/skills/` (listed by `git status`).

### Phase 1: the command-line chain (2026-11-05, the dell)

A command-line key passes through four places: its declaration in `type(parameters)`
(`src/main/params/simple_parameters.f90`), the list of accepted keys that the build generates from the
declarations (`scripts/simple_args_generator.pl`, compiled as `simple_args.f90`), its binding onto the field
(`simple_parameters_parse.f90`) and the program descriptors that offer it in the help and the graphical
interface (`src/main/ui`). The phase makes a gap between them fail the build gate instead of a user's run.

**The gate check** `scripts/check_cli_chain.py` (Python 3.6, standard library; 0.3 s) runs in
`scripts/run_fast_gate.sh` before the unit suites. It runs the generator itself into a temporary directory, so
it and the build agree on the accepted keys, and fails on:

1. a UI input, or a key of the private programs' table, that is no accepted key (giving it stops the program);
2. an accepted key that `parameters` never binds (accepted and silently ignored), unless `parameters` reads it
   from the command line itself or it is the one listed exception (`nrots`, read by `make_pickrefs` from its
   command line; the field of that name is the search grid the builder derives);
3. a key bound twice;
4. a yes/no key bound without `add_yesno`, `add_yesno` on another key, or a yes/no (`binary`) UI input whose
   declaration names no `(yes|no)` choices;
5. a key of the hand-written stream launcher (`print_stream_ui_json`) that is no accepted key;
6. a program of `simple_private_exec` with neither a public descriptor nor an entry in
   `src/main/ui/simple_private_prgs.f90` (its command line was never checked), a table entry no dispatch case
   runs, or a name in the table twice.

On the base tree it reports 196 problems (`scratch/logs/phase1_check_base.log`): the four UI inputs of ruling 9,
`rec_states` bound twice, 93 accepted keys never bound, 94 yes/no keys without validation, 18 yes/no UI inputs
without declared choices and the four private programs without an entry. On this phase's tree it passes. A
violation of each kind planted into a copy of the tree (nine plants: an unbound key, a new declaration without a
binding, a UI input that is no key, a private program without an entry, a table key that is no key, a key bound
twice, a yes/no key bound with `add_char`, a yes/no input whose declaration lost its choices, a stream launcher
key that is no key) fails it, and the unchanged copy passes (`scratch/logs/phase1_plant.sh`, `phase1_plant.log`;
the plants were made in a temporary copy, never in the repository).

**Decision: derived fields are no keys.** The rule "every accepted key is bound" found 111 declared fields that
were accepted but never bound, not the two the audit had named. 86 of them (`cwd`, `exec_dir`, the `l_*` flags,
`box_croppd`, `refs_even`, `vol`, ...) are values `parameters` derives or sets itself, and no code puts any of
them on a command line (searched in `src/`, `production/`, `nice/` and `scripts/`). Their declarations now say
"not a command-line key" in the comment, the convention `inpl_cont` already used, and the generator leaves such
lines out. A user who gives one now gets the unknown-key stop instead of having it ignored. The others are bound:
`mul` and `import_type` (ruling 8: `oriops`/`volops` now apply `mul`, `import_starproject` now runs its
destination checks), and `import_dir`, `xmldir`, `xmlloc`, `test`, `istart`, `reliongroups`, `tiltgroupmax`
(UI inputs or keys programs read from their command lines; their fields were read nowhere, so binding changes
nothing else). No parameter was added.

**Rulings carried out**

- R1 (ruling 8): `mul` and `import_type` bound; the duplicate binding of `rec_states` removed.
- R2 (ruling 9): no program reads `positive` (`automask2D`), `heldout` (`flex_pca`) or `rmsd_file`
  (`atoms_stats`): the inputs are dropped. `trajectory_graphene_subtr` offered `stk_traj` (no key) and
  `stk_backgr` (a key it ignores) while it reads `stk` and `stk2`: its two inputs now name `stk` and `stk2`,
  with the old labels; the shared `stk_traj` descriptor, its placeholder and its file-suffix entry in the
  JSON validator, made unused, are removed.
- R3 (ruling 10): the generator reads `type(parameters)` from its first line to its `contains` (or `end type`),
  drops the commas of array constructors, stops on anything that is not a Fortran name, and writes its temporary
  files into its output directory. The keys `cntfile`, `silent`, `final_planned`, `1.0]` and a second `ext` are
  gone.
- P6: the parameter registry has `add_yesno`, which stops on any value other than `yes` or `no` before the
  field's length would cut it (`len=3` turned `yess` into `yes`). 130 keys use it: the 93 declared with the
  choices `(yes|no)` and 37 that took only yes or no without saying so, which now declare `(yes|no){default}`.
  `scripts/default_audit.py`, which derives the UI's displayed defaults from the bindings, reads `add_yesno`
  (the generated defaults are identical to the base's). Known effect: `pose_cont=only`, retired, now stops with
  the generic yes/no message before the specific hint in `simple_parameters_phases.f90` is reached
  (`refine3D_auto` still prints the hint, which checks the value before `params%new`).
- P3: the private programs' table is built with `add_private_prg` and `push_req`/`push_opt`; the count follows
  (`n_private_prgs` was the literal 17) and a name given twice stops. New entries for the four workers in neither
  table: `refine2D_distr` (requires `projfile`, which the stream's 2D pool always sets), `calc_pspec_assemble`
  (`projfile`, `nparts`; run in-process only), `aggregate_chunks` (`projfile`, `projfile_target`, which its
  commander reads) and `scale_project` (`projfile`).
- P2: every build now checks the stream launcher (`prg=print_ui_stream`, written by hand) against the program
  descriptors (`production/validate_complete_ui_json.py`, given `simple_private_exec` with
  `--launcher-executable`): each process is a `simple_stream` program, each user input it gives a process is
  defined, is an input of that program and carries its default. On the base build this fails: `preproc` was given
  a user input `scale` that the launcher never defines and `preproc` does not take; the reference is removed.
  The opt-in contract targets (`validate_stream_ui_contract`, `update_stream_ui_contract`) now also compare and
  update `production/stream_launcher_contract.json`, created here; a planted change of it is reported
  (`scratch/logs/phase1_plant_launcher_contract.log`). The programs' reference
  `production/stream_ui_contract.json` already differed from the descriptors at the base (the base fails the
  contract target, `scratch/logs/phase1_base_stream_contract.log`); it is left as it was, for the maintainer's
  review, and the validator now compares the two references independently so a stale one does not hide the
  other's diff.
- X3: `simple_stream` dispatches `gen_pickrefs` only; the aliases `initial_analysis` and `opening2D` (which
  never ran, because the command-line parser stops on an unregistered name first) are gone.

**Tests.** `unit_core suite=parameters` gained three checks: the accepted keys (no dummy argument after the type,
no array-constructor fragment, no derived flag; `mul` and `import_type` present), the new bindings survive
`params%new`, and a yes/no key given `yess` or `true` stops a forked child process while `yes` lets it finish. The
same tester compiled against the base library fails 10 of its 16 checks (`scratch/logs/phase1_base_tester_run.log`;
the scratch build was deleted after the run) and passes 16 of 16 now.

**Gates and runs.** Release and Debug: fast gate 15/15 with the new check (`scratch/logs/phase1_release_gate_final.log`,
`phase1_debug_gate_final.log`); `simple_test_exec test=units` passes in both (`phase1_units_build.log`,
`phase1_units_build_debug.log`); the `library` label 7 of 7 (`phase1_library.log`, 262 s). The workflow entries
of 8.2 on the phase's final Release build, `SIMPLE_SEED=20260923`, against the phase 0 base runs
(`simple_test_exec test=identity`, then the differing project files compared as text with the run directory and
`solve3D_run_id` normalised): every declared output but the project files is byte-identical (6vxx 177 files,
1jxy 132), and the project files differ only in their run paths, the `solve3D_run_id` time stamp and the padding
these change (`scratch/logs/phase1_identity_final.log`, `phase1_identity_final_*.txt`). `mini_stream_1jxy`, whose
run directory carries a process id, compared by `diff -r`: images, boxes and STAR files identical, the rest process
ids, paths and timings (`phase1_identity_mini_stream_1jxy.log` and the journal).
`simulated_workflow_6vxx` gives the base's numbers (whole-volume correlation 0.5611, FSC 0.143 at 20.8 Å),
`simulated_workflow_1jxy` stops where the base stops (`new_shiftmap`, solve3D stage 3), `mini_stream_1jxy` passes.

**Guides updated**: `simple-main-params`, `simple-main-ui` and its `references/cli_change_chain.md`,
`simple-main-exec`, `simple-architecture` (the gate's checks).

**Found, for later phases**: `src/main/motion/simple_motion_model_tester.f90` sets `scale_movies` on a command
line that `parameters` never read (the line had no effect before and is no key now; motion/, phase 2). The
yes/no rule covers keys declared or offered as yes/no; a key that takes other values (`automsk`, `ctf`, `mirr`,
`stats`, ...) keeps its own validation where it has one.

Files changed: `scripts/simple_args_generator.pl`, `scripts/check_cli_chain.py` (new), `scripts/run_fast_gate.sh`,
`scripts/default_audit.py`; `src/main/params/simple_parameters.f90`, `_parse.f90`, `_registry.f90`, `_tester.f90`;
`src/main/ui/simple_private_prgs.f90`, `simple_ui.f90`, `simple_ui_param.f90`, `simple_ui_params_common.f90`,
`simple/simple_ui_mask.f90`, `simple/simple_ui_heterogeneity.f90`, `single/single_ui_atom.f90`,
`single/single_ui_trajectory.f90`, `simple_test/simple_test_ui_highlevel.f90` (the `polar_matching_bench` help
named its note under `planned/`; it is in `completed/`, found in phase 0); `production/simple_stream.f90`, `production/validate_complete_ui_json.py`,
`production/CMakeLists.txt`, `production/stream_launcher_contract.json` (new); this plan; the guides above.

### Phase 2: the small fixes and the checks of 3.1 (2026-11-05, the dell)

Two kinds of work: the small code changes the 2026-11-04 rulings asked for (section 7), and the rules of 3.1 that
failed silently and now stop. Every changed rule has a tester check that commits the violation in a forked child
process and expects the stop, or computes the corrected value. Unless an item says otherwise, the check fails on
the base library (the evidence build is under **Tests**).

**Section 7, the small items**

- C5, a worker writing the project (narrowed). A partition worker (`part=` given, not `stream=yes`) now stops when
  it would write its input project: `write`, `write_segment_inside`, `clear_segment_inside` and
  `write_non_data_segments` (`project/simple_sp_project_io.f90`). `params%new` records the worker's project file
  (`worker_projfile_glob`, `src/defs/simple_defs.f90`). It writes part files as before. A stream job (`stream=yes`)
  owns its job project and keeps writing it, and the master is no worker. Layering policy section 5.2 and the
  project guide state the rule as checked.
- O2: `rot_to_asym` (`ori/simple_sym.f90`) starts from the input record, so a pose already in the asymmetric unit
  (under c1, a south-pole record, which no symmetry operation moves) keeps its angles and its other keys instead of
  coming back zeroed.
- X1: a particle with no even/odd half (`eo=-1`) is no longer added to the even half's sums. The new
  `oris%get_half` stops on any value other than 0 or 1. The class averager calls it in a serial pass before its
  parallel region (`class/simple_classaverager_restore.f90`), and 3D gridding calls it in both of its groupings
  (`strategies/search/simple_matcher_3Drec.f90`, where `normalized_eo` and its counters are gone). Not covered,
  because the files are outside this phase: the iterative (PCG) reconstruction backend and `flex/` drop such a row
  without a word. The GPU-offload reconstructor and the scoring of a batch still treat it as odd and even.
- X5: `center_pdbcoord` (`nano/simple_atoms.f90`) puts the atomic centre at voxel ldim/2+1, that is at
  `ldim/2*smpd`. The old `(ldim-1)/2*smpd` put an even box's centre half a voxel off. This moves every simulated
  volume made from a PDB file by half a voxel, so it changes all three workflow entries (below). The `cif2mrc`
  tester's analytical reference now uses the same centre; its tolerances are unchanged.
- O6: the optimiser's clients deallocate the optimiser they allocate in `kill`: `pftc_shsrch_grad`,
  `ctf_estimate_cost` (`kill4Dcont`) and `volpft_symsrch`. No tester check: the optimiser is a private pointer.
- O3: `opt_spec` has `nostart`. An all-zero `x` is a start like any other, and the simplex draws a random start
  only under `nostart=.true.`. `specify` resets `nostart` and `npop`.
- O5: differential evolution uses the `npop` given to `specify`, falling back to its mode's default, and places
  the members of `set_inipop` after the start point. The random members are drawn first, so a run without an
  initial population draws the same stream as before.
- O7: `ctf_estimate_fit%kill` no longer forgets the process-wide Fourier-index maps (`memoize_ft_maps`), which
  other users in the process still hold.
- O9: motion correction records `mcconvention` on the micrograph's row (`motion/simple_motion_correct_iter.f90`).
  Extraction from the movie frames (`extractfrommov=yes`) stops unless the convention is `simple`
  (`preprocess/simple_particle_extractor.f90`).
- O10: `nframesgrp` outside a time series stops in `motion_correct_iter%iterate`. Without the stop, the grouped
  movie would be recorded as the movie. This check does not fail on the base, whose child also stops, on the
  missing movie the test gives it. Also in the motion tester: the `scale_movies` line found in phase 1 is now
  `smpd_downscale=3.0`, with an assertion that the parameters scale movies by 0.5.
- O13: when picks are assigned to their diameter bins, the segmentation picker adds the 8 Å that binarization
  erodes (`2*SMPD_SHRINK1`), as its bin decision already does (`pick/simple_segdiam_bin_picker.f90`). The box file
  still records the measured diameter.
- O14: the segmentation-diameter picker writes the denoised thumbnail (`thumb_den`): `pick_1` takes the thumbnail
  name, and `exec_segdiampick`, the picker iterator and `simple_pick_strategy.f90` pass it on. The strategy stores
  the key only when a name came back.
- I2: `polarize` and `polarize_oversamp` take the caller's polar dimensions (`pdim`) and stop unless they, the box
  and the mode (normal or oversampled) are those of the memo. Memoising inside an OpenMP parallel region stops.
  Callers updated: `simple_polarft_core.f90`, `simple_corrmat.f90`, the graphene subtraction and the Cartesian
  calculator tester. `mem_poldim` moved into the polarizer submodule.
- I5: the Fourier shift constant is `2*pi/n` for an odd box `n` as well; the integer `n/2` gave `2*pi/(n-1)`.
  Fixed in the image (`simple_image_core.f90`, `simple_image_calc.f90`), `simple_ft_expanded.f90`, the polar
  calculator's argument-transfer constant (`simple_polarft_core.f90`) and the class-sum shift
  (`class/simple_cavg_sums.f90`). The last two were found by the diff review; they keep the polar search, the
  class sums and `image%shift` in agreement on odd boxes. The workflow entries use even boxes and do not change.
- I3: `cartft_calc%set_sigma2` (the continuous refinement's Cartesian calculator) unprepares the Euclidean particle
  slots prepared under the previous spectra. Two Cartesian tester routines relied on the old behaviour and now set
  the spectra first.
- S2: `export_cls2D` works on a copy of the class rows, so the project no longer gains `ncls` or loses `stk`.
- S3: `export_cls2D`, `export_opticsgroups`, `export_stardata` (which now refuses a name not ending in `.tmp`) and
  the three RELION writers write `<name>.tmp` and rename it.
- N3: `get_lattice_params` stops on an element it has no lattice for (`src/defs/simple_defs_atoms.f90`).
- N1: the graphene subtraction (`nano/simple_tseries_graphene_subtr.f90`) takes the caller's parameters, so it no
  longer resets the process's threads and seed through a `params%new` of its own. It also sets its particle
  instead of overwriting the reference. The commander passes its parameters
  (`commanders/single/single_commanders_trajectory.f90`).
- N6: `score_project_cavgs` builds its relation parameters under the caller's team. It restores the caller's thread
  count, team, directory, random state (`get_rnd_state` and `put_rnd_state`, new in `src/utils/math/simple_rnd.f90`)
  and, after the review, the worker identity (`part_glob`, `l_distr_worker_glob`, the recorded project).
- N9: the flex PCA basis restores OpenMP's maximum active levels after it changes them; checked in the
  `lib_heterogeneity` entry.
- V1: already done at HEAD (51dd0a627); nothing to change.
- P1: layering policy section 8 has the routers' extra duties: "Routers choose by key" and "Front doors write the
  project". The exec guide points there.
- N8: the flex gateway comment and `src/main/flex/README.md` name the project the run writes: the `mkdir=yes`
  copy, or the caller's under `mkdir=no`.

**The rules of 3.1**

| Rule | Outcome |
| --- | --- |
| A new particle batch is memoised before scoring | done: the calculator counts particle changes (`ptcl_gen`, bumped by every particle setter: `reallocate_ptcls`, `set_ptcl_pft`, `set_ptcl_ctfmat`, `polarize_ptcl_pft`, `set_ptcl_fcomp`, `cp_even_ref2ptcl`, `shift_ptcl`, `create_polar_absctfmats`); `memoize_ptcls` records it; nine scoring entries stop on a stale memo (`check_ptcl_memo`, once per call) |
| The sigma2 weights follow the spectra | done: `euclid_sigma2` counts spectra changes (`spectra_gen`, bumped by `new_cart`, `read_global` and, after the review, `kill`); the calculator points at the counter (`assign_sigma2_noise`); `memoize_sigma2_weights` records it; six Euclidean scorers stop on newer spectra (`check_euclid_weights`) |
| `prep_sigmas_objfun` after every `pftc%new` | done by the same check: a Euclidean score without weights stops (the base crashed with a segmentation fault) |
| A pass runs with `omp_get_max_threads() == nthr_glob` | done: `check_pass_team` at the start of `refine3D_exec` and `refine2D_exec` stops when the team is larger |
| The master writes its sample before workers run | deferred: a worker-side check needs the sample's iteration in the project. That means a new per-iteration project write from the sampling routines, which layering policy section 8 forbids adding, or a change to the refine3D master (`simple_refine3D_strategy.f90`, a file of `ms_posterior2`) |
| `finalise_ptcl_import` once, after import | done: `star_file` has a flag set by the imports. A second finalise, or an export or optics assignment of an unfinalised import, stops |
| The included records stay a prefix (sieve) | deferred: `src/main/sieve/` is out of bounds for this run |
| `cleanup_nu_filter` between states | done: `setup_nu_dmats` stops while the previous pair's solvent clamp or null-shell mask is still set |
| A gridding file matches the reader's box | phase 3 (V2) |
| CTF application follows `memoize_ft_maps` for its box | done: `ft_maps_match` (`image/simple_memoize_ft_maps.f90`); `ctf2img` and `apply_ctf` stop on another box |
| `polarize*` matches the memoised geometry | done (I2) |
| A model's feature vectors follow the feature bank | done: `read_feature_weights` and the linear coefficients stop on a length other than the bank's (`cavg_quality/simple_cavg_quality_model.f90`) |
| A particle row carries only its named slots | deferred: a debug-build stop in the particle writer would stop four producers that write other keys today: `trajectory_solve3D` (`frame`), `import_particles` (CTF keys), the STAR import (optics keys) and the binary-orientation tester. Each needs an owner's decision: make the key a slot, or strip it. |

Several of these stops are reached inside OpenMP parallel regions: the memo and weights checks at the top of the
scorers, `polarize*` in the batch preparation, and `apply_ctf` in the class averager. They are fatal stops for a
programming error, with fixed messages and no shared state written, but several threads may print the message at
once.

**Found by the diff review and not acted on** (for the maintainer):
- I3 marks a stale slot invalid, which its consumers treat as a degenerate particle, without an error. The
  production order (spectra before the batch) is right.
- `simple_rename(..., overwrite=.true.)` deletes the target before renaming (`src/fileio/simple_syslib.f90`), so
  S3's write-and-rename still leaves a moment with no file. No data is lost: the `.tmp` file stays.
- `ft_maps_match` compares the box, not the sampling, as ruled.

**Tests.** The new checks are in the area suites: optimisers (a zero start; `npop` and the initial population),
CTF (`kill` keeps the maps), symmetry (`rot_to_asym` under c1), orientations (`get_half`), STAR project (finalise
once, `export_cls2D` on a copy), project records and parameters (C5), image (odd-box shift constant, CTF maps of
another box), polar correlation (six call-order violations: no weights, newer spectra, an unmemoised batch, and
`polarize` with another box, mode or `pdim`), Cartesian calculator (new spectra unprepare a slot), atoms (centre,
lattice lookup, graphene subtraction keeps threads and seed), `cif2mrc`, nonuniform filtering (the clamp),
class-average quality (feature lengths), particle extractor (O9), motion (O10), random draws (the generator's
state saved and restored), in-plane strategy (the pass's team) and the stream's initial analysis (O13). Each stop
is checked in a forked child process, by its status and its exit code, so a crash does not pass for a stop. A child
that writes the project first sets one OpenMP thread: after a fork, the parent's thread pool is not in the child,
and a write that should have stopped would otherwise hang the suite instead of failing it.

Evidence on the base: the base library with this phase's testers. Four testers were edited to the base's
signatures in the scratch copy (atoms, polar correlation, Cartesian calculator and project records), and the five
that need a new interface (random draws, optimisers, orientations, parameters, in-plane strategy) kept their base
versions. That build fails 13 suites with 134 assertions, every one of them a new check
(`scratch/logs/phase2_base_overlay_units.log`; the scratch build was deleted). The base child of the "no weights"
check died of a segmentation fault (exit code -11). The O10 check does not fail on the base (see O10).

**Gates and runs.** Release and Debug: fast gate 15/15 (`scratch/logs/phase2_release_gate_final.log`,
`phase2_debug_gate_final.log`); `simple_test_exec test=units` 152 groups, 14,862 assertions, all passed in both
(`phase2_units_release_final.log`, `phase2_units_debug_final.log`); `check_guide_imports`, `check_test_registry`,
`check_descr` and `check_cli_chain` pass; no file mode changed. The `library` label passes 7 of 7 (237 s, `phase2_library_final.log`).

The workflow entries of 8.2 (`SIMPLE_SEED=20260923`), against the phase 0 base runs (the phase 1 runs were identical
to them apart from paths and run identifiers):

- On this phase's final Release build, every output differs from the base. X5 moves the simulated volume, and with
  it every simulated movie, by half a voxel before any processing (`scratch/logs/phase2_identity_final_vs_base.log`).
  `simulated_workflow_1jxy` now runs to the end (rc=0, 14:36): whole-volume correlation 0.6410, docking 0.8669; over
  two parts 0.6379 and an FSC (Fourier shell correlation) of 0.143 at 35.66 Å. The base stopped in solve3D stage 3
  (`new_shiftmap`, the map's last shell outside the band). `simulated_workflow_6vxx` gives whole-volume correlation
  0.5185 and an FSC of 0.143 at 49.92 Å, against the base's 0.5611 and 20.8 Å. Both are below the entry's 0.80 floor,
  a missed map as at the base, and the entry reports without failing. The entry's result does not depend on
  `SIMPLE_SEED`: with 20260924, the base and this build each reproduce their 20260923 outputs byte for byte
  (`phase2_seed2_*_6vxx.log`, `phase2_identity_seed2_phase2.txt`). A second seed therefore gives no measure of
  run-to-run spread, and one data set cannot tell whether the moved input made 6vxx worse and 1jxy better or
  only different. `mini_stream_1jxy` passes (pick ratio
  0.417, base 0.438; 2 classes).
- Attribution: the same tree with only X5 reverted, built in scratch (`scratch/logs/phase2_identity_attr_noX5.log`,
  `phase2_attr_noX5_*.log`). `simulated_workflow_6vxx`: 177 files byte-identical to the base, whole-volume
  correlation 0.5611 and FSC 20.80 Å, as at the base. `simulated_workflow_1jxy`: 132 identical, stopping where the
  base stops. The project files differ only in two keys: `mcconvention=simple` on every micrograph row (O9), and
  `thumb_den`, blank at the base, now the denoised thumbnail's path (O14). A micrograph with no picks gets no
  `thumb_den` key; its thumbnail file is written but not referenced. `mini_stream_1jxy`: images, boxes and STAR files
  identical; the rest are paths, process ids, timings and the sigma2 state's layout digest, which hashes the stack
  paths. So apart from X5, only the ruled O9 and O14 change a workflow's output. I5's odd-box constants and the
  other checks change nothing here (the entries use even boxes).
- The final build against an earlier build of this phase, before the review's fixes (I5 in the polar calculator
  and class sums, `esig%kill`, the restored worker identity): identical apart from paths
  (`phase2_identity_final_vs_runA.log`).
- The final runs are kept as phase 3's reference: `scratch/keep/phase2_runs` (line in `scratch/keep/summary.tsv`).

Timings were taken with other jobs running (load 2 to 22): 1jxy 14:36, 6vxx 13:02, mini_stream 1:00.

**Exit items**
- Every item: 24 small items done (V1 found already done). Of the 13 rows of 3.1, nine are done, three deferred
  with their reasons, and one belongs to phase 3.
- Identity comparison of 8.2: done above. The differences are explained by X5, O9 and O14, all ruled.
- Fast gate 15/15 in Release and Debug.
- Unit suites of the touched areas: all of `test=units`, both builds, plus the `library` label.
- `check_guide_imports`, `check_test_registry`, `check_descr`: pass.

**Guides updated**: `simple-main-opt`, `-pftc`, `-ctf`, `-strategies`, `-strategies-search`, `-star`, `-ori` (and
`references/conventions.md`), `-project`, `-class`, `-nano`, `-flex`, `-nu-filt`, `-pick`, `-motion`, `-exec`,
`-image` (and `references/fourier-and-memos.md`), `-cftc`, `-sigma2`, `simple-cluster-cavgs-quality`. Each
converted rule names its check. `simple-cluster-cavgs-quality` and `simple-main-stream` contained committed
merge-conflict markers, which are resolved.

**File table**: rows added for `commanders/single/single_commanders_trajectory.f90` (N1),
`stream/stages/simple_stream_stage_initial_analysis_tester.f90` (O13) and `volume/` (V1, O6); widened: `params/`
and `src/defs/` (C5), `src/utils/` (N6), the layering policy (P1, C5), `flex/` (N8) and `strategies/search/` (the
pass's team, the `polarize` callers).

Files changed: `src/defs/simple_defs.f90`, `simple_defs_atoms.f90`; `src/utils/math/simple_rnd.f90` (+ tester);
`src/main/` `cavg_quality/simple_cavg_quality_model.f90`, `_selection.f90`, `_lattice_tester.f90`;
`cftc/simple_cartft_calc.f90` (+ tester); `class/simple_cavg_sums.f90`, `simple_classaverager_restore.f90`;
`commanders/single/single_commanders_trajectory.f90`; `ctf/simple_ctf_estimate_cost.f90`, `_fit.f90`,
`simple_ctf_tester.f90`; `flex/README.md`, `fit/simple_flex_pca_basis.f90`,
`run/simple_flex_pca_application_tester.f90`, `run/simple_flex_pca_project_gateway.f90`;
`image/simple_ft_expanded.f90`, `simple_image.f90`, `_calc.f90`, `_core.f90`, `_ctf.f90`, `_polar.f90`,
`_tester.f90`, `simple_memoize_ft_maps.f90`; `motion/simple_motion_correct_iter.f90`, `simple_motion_model_tester.f90`;
`nano/simple_atoms.f90` (+ tester), `simple_cif2mrc_tester.f90`, `simple_tseries_graphene_subtr.f90`;
`nu_filt/simple_nu_filter_bank.f90` (+ tester); `opt/simple_opt_de.f90`, `_simplex.f90`, `_spec.f90`, `_tester.f90`;
`ori/simple_oris.f90`, `_getters.f90`, `_tester.f90`, `simple_sym.f90` (+ tester);
`params/simple_parameters_phases.f90`, `_tester.f90`; `pftc/simple_pftc_shsrch_grad.f90`, `simple_polarft_calc.f90`,
`_core.f90`, `_corr.f90`, `_corr_tester.f90`, `_ctf.f90`, `_geom.f90`, `_memo.f90`, `_shiftmap.f90`;
`pick/simple_picker_iter.f90`, `simple_picker_utils.f90`, `simple_picksegdiam.f90`, `simple_segdiam_bin_picker.f90`;
`preprocess/simple_particle_extractor.f90` (+ tester); `project/simple_sp_project_io.f90` (+ tester);
`sigma2/simple_euclid_sigma2.f90`; `star/simple_relion.f90`, `simple_starproject.f90` (+ tester),
`simple_starproject_utils.f90`; `strategies/parallelization/simple_pick_strategy.f90`; `strategies/search/simple_corrmat.f90`,
`simple_matcher_3Drec.f90`, `simple_matcher_ptcl_batch.f90`, `simple_matcher_refvol_utils.f90`,
`simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`, `simple_strategy3D_inplane_tester.f90`;
`stream/stages/simple_stream_stage_initial_analysis_tester.f90`; `volume/simple_volpft_symsrch.f90`;
`doc/policies/layering_and_project_ownership_policy.md`; this plan; the guides above.

### Phase 3: the medium items, renames and the glossary (2026-11-05, the dell)

The medium items of section 7 that are in bounds: X2, C1, V2, N10, O11, V5 and G3. Also the renames of 3.5,
the glossary of 3.6, and the language-server experiment of 3.8. Four read-only investigations, one per item
group, mapped the code first; their designs were adopted with the changes noted below.

**Items**

- X2, the nested reconstruction of `refine_motion_model`. `simple_rec3D_strategy.f90` has a strategy routine,
  `run_rec3D(cline)`. It holds the old body of `commander_rec3D`: the command-line defaults, the strategy
  lifecycle and the stop line. `commander_rec3D` is now a thin wrapper of it.
  `simple_refine_motion_model_strategy.f90` calls `run_rec3D` and no longer imports a commander. The same code
  runs in the same order, so its outputs do not change. `refine_motion_model_policy.md` and the
  parallelization guide name the routine. Exercised by `lib_reconstruction` (suite "fraction dose export",
  `reconstruct_motion_fractions` end to end). Still importing `commander_rec3D` is
  `simple_refine3D_strategy.f90`, a file of `ms_posterior2` and out of bounds. Layering policy section 8 already
  lists it.
- C1, a row range read from a `.simple` file. `read_segment(..., fromto)` now fills only the rows
  `fromto(1):fromto(2)`, at their index in the segment, of a field sized to the whole segment; the other rows
  stay fresh. This holds for `ptcl2D`, `ptcl3D` and `cls3D`, as the text branch already did. A range with any
  other segment now stops, where the text branch ignored it. The binary reader's string segments (`cls3D`)
  ignored the range and stored rows at batch-local indices; they now read from the range's first byte, and a
  read without a range is unchanged. No production caller passes a range today (the one caller, in
  `init_params_spproj_tbox2D`, has no callers itself), so no workflow changes. Test: `unit_project`,
  `test_read_segment_range`.
- V2, the gridding density files (`rho_*`: worker partials, trailing chains, frozen sets). They carry no box,
  so a reader took the receiver's shape and read a larger file as a prefix. They now start with a 64-byte header:
  a magic word, the box, the density's array shape and the sampling (`simple_reconstructor.f90`, `write_rho`).
  `read_rho` and the readers of open files check it first: `check_rho_header` in volume assembly and the frozen
  accumulator, serially, before the four-thread reads. They stop on another box or sampling and on a file
  without the header. The `update_frac` path that pads a smaller previous chain checks the header's box against
  its numerator, and the physical extent within 1%. The trailing-chain manifest's version is 4, so chains written
  before the change are re-seeded rather than stopping. The numerator files are still read in the reader's
  shape: they are MRC files whose header the readers skip. Test: `unit_reconstruction`, `test_rho_header`.
- N10, the settings FLEX's trial reconstructions override. `rec3D_request` now carries them: `nstates`,
  `one_part`, `no_nonuniform`, `no_trail_rec`, `all_states` and `no_automsk`. `rec3D_service%execute` applies them
  to the parameters for the reconstruction and restores the caller's values (and `vols`) at its end.
  `reconstruct_state_halves` sets the request and no longer keeps its own save-and-restore list. The defaults
  mean "the caller's value", so the three strategy callers are unchanged. The investigator's alternative, the
  consumers reading the request instead of `params` (about 50 reads across `calc_3Drec`, the PCG worker and the
  PCG master), was not taken: the PCG master is also called from the out-of-bounds refine3D strategy. The
  command line handed to volume assembly still carries the same settings as keys. Test:
  `lib_reconstruction` "fractional reconstruction", `test_request_settings`: a one-column request under a
  caller with two states and `state=2` reconstructs one state, and the caller's `nstates`, state selection and
  `state` come back. It uses the new fields, so it cannot run on the previous build. `lib_heterogeneity` runs
  the FLEX path end to end.
- O11 and V5, moves with `mv`, module names kept: `preprocess/simple_pspecs.f90` to `strategies/search/`, next to
  its one user, and `volume/simple_dose_series_manifest.f90` with its tester to `motion/`. The tester stays in the
  `unit_reconstruction` suite, so the test UI is unchanged. CMake collects `main/*.f90` recursively; both builds
  were re-configured. Guides: motion, the volume lifecycles reference, strategies-search.
- G3, per-thread search state sized by the team, partly. The pass's own state is now sized by the OpenMP team
  its regions run (`pass_team()`, `simple_matcher_ptcl_io.f90`): the 3D and 2D search singletons (`s3D`, `s2D`),
  the particle images (`alloc_ptcl_imgs`), the reconstruction's image heaps and the 2D references' work images.
  The work-image dummy arguments are assumed-shape. Three things still size by `nthr_glob`: the probability
  tables (`probabilistic/`, out of bounds), the polar calculator's work buffers, and the helpers outside a pass
  (`corrmat`, the 2D utilities). A team larger than `nthr_glob` would still write out of bounds there, so
  `check_pass_team` stays. Test: `test_pass_state_team_sized` (in-plane strategy suite): with `nthr_glob=1`
  and a team of 2, the particle images number 2. No parallel region's arithmetic changed, only the sizes
  (equal to before whenever the team is `nthr_glob`).

**3.5, descriptive names for numbered specifics.** All the numbered specifics of generics in the folders the
file table allows were renamed: `src/fileio`, `src/main/image`, `src/main/ori`, `src/main/pftc` and `src/utils`.
That is 149 generic statements dispatching to 371 numbered specifics, plus 9 numbered procedures named the same
way but outside a generic: 380 names. Examples: `set_1..5` of `oris` became `set_real`, `set_char`, `set_int`,
`set_dp` and `set_str`; `mirror_ref_pft_1/2` became `mirror_ref_pft_by_pft` and `_by_ind`;
`gen_corr_for_rot_8_1/2` became `_noshift` and `_shifted`; `locate_1/2` became `locate_real` and `locate_int`. The
generic names stay. A script (`scratch/logs/phase3_rename.py`, map `scratch/logs/numspec_map.txt`) replaced whole
words in each module's own files only: the parent, its submodules, its tester beside it, and the one outside
caller, `flex/fit/simple_flex_pca_embed.f90` (`locate_2`). That was 1,905 replacements in 55 files. Names that
repeat across modules (`set_1` in `ori`, `oris`, `hash`, `chash`, `cmdline`) were therefore renamed per module.
Three proposed names collided with existing ones and were changed: `append_at`, `same_dims_img` and
`split_string`. The `apply_filter` binding, whose implementation carried the generic's own name, became
`apply_filter_single` by hand. The compiler checks the rest: Release and Debug build without error. Left
numbered: the specifics of `src/main/simple_cmdline.f90` (`set_1..6`; not in the file table) and the
non-generic `open_1` of `simple_discrete_stack_io.f90`. The Fortran indexes under
`doc/code_overview/fortran-indexes/` were regenerated (`scripts/gen_fortran_indexes.pl`), and the planned note
`x_execution_layer_adoption.md` names the new `generate_script` specifics.

**3.6, the glossary.** `doc/glossary.md` lists the recurring concepts under their code names today:
- the workflows, passes, master/worker/part;
- the project segments, `state`, the even/odd half under its three names, `sampled`/`updatecnt`;
- `box`/`box_crop`, `smpd`/`smpd_crop`, `kfromto`, `trs`, `nspace`;
- the search and score terms;
- the code layers;
- the name conventions (`_8`, `l_`, `_glob`, `new`/`kill`, generic and specific).

`simple-architecture` links it.

**3.8, a Fortran language server (experiment; nothing of it is in the repository).** A local plugin in
`scratch/keep/fortls_plugin/`: `.lsp.json` runs `~/.local/bin/fortls` 3.2.2 for `.f90` files, and `fortls.json` is
its configuration. It lists `src/**`, `production` and the generated `build/lib/simple` as sources and
`src/inc` for includes, and preprocesses `.f90`. A stdio client (`probe.py`, logs in `logs/`) measured it, and a
separate `claude --plugin-dir` session on the repository used it (`logs/session.json`).

| Measured | Result |
| --- | --- |
| Indexing | 832 files, 393,000 lines: 1.6 s, 216 MB resident, unchanged after the queries |
| Query time | definition and implementation 20 to 660 ms; references 0.5 to 1.7 s |
| Definition of a call to a type-bound procedure or generic (`img%polarize`, `o%set`, `pftc%new`, `pftc%gen_corr_for_rot_8`) | the binding line in the type, not the specific's body |
| Into a submodule | "implementation" on the call gives the parent's interface block; a second "implementation" on that gives the submodule body (`simple_image_polar.f90`, `simple_sp_project_io.f90`). The session agent found both bodies this way |
| Polymorphic call of a deferred binding (`xautomask%execute`, declared type known) | the concrete commander's `execute` body |
| `use` chains | an imported variable (`nthr_glob`) resolves to `simple_defs.f90`; its 218 references are exact |
| References | exact for a type's own binding (`oris%get_half`: 12 in 4 files; `pftc%new`: 31); for a generic common to many types (`set`) 1,921 hits: by name, not by type |
| Generated sources | `simple_args` resolves only with `build/lib/simple` among the sources; without it the definition is empty and fortls reports "Module simple_args not found" (information) |
| Diagnostics, 80 random files of `src/main` | no error from the `THROW_*` macros (with `src/inc` and preprocessing on); 32 false errors in one file: names of the JSON library imported through `json_module` ("not found in scope"); 96 warnings "variable masks variable in parent scope" (real shadowing by a contained procedure's local) |
| A separate session | the LSP tool was offered (deferred, loaded by search); 13 calls answered five navigation questions in 60 s; the first call met "server is starting" and the retry worked |

Assessment: `fortls` copes with SIMPLE. It indexes the whole tree in seconds with little memory, which matters
for an agent starting cold. It resolves the `use` chains and the type-bound bindings, and it reaches submodule
bodies in two steps. Its references are exact for type-specific bindings, which plain search cannot do for the
254 bindings named `new`, `kill`, `get`, `set` and `execute`. Its weaknesses:
- generics give the binding, not the specific chosen by the arguments;
- references to a generic shared by many types are name-based;
- the third-party JSON module gives false errors;
- the generated sources need the build directory configured.

Recommendation for the maintainer: worth adopting as a user-level plugin, with the build directory among the
sources. The false JSON errors are the one source of noise in diagnostics.

**Evidence on the previous build.** The phase-2 tree (the base sources plus the driver's phase-start snapshot)
was built with this phase's project and frozen-accumulator testers. There, `test_read_segment_range` fails 5 of
its checks: the rows outside the range are read, and a range on the mic segment does not stop.
`test_rho_header` fails 3: a smaller box, another sampling and a file without the header are all read without a
stop (`scratch/logs/phase3_p2tree_unit_project.log`, `phase3_p2tree_unit_reconstruction.log`; the scratch build
was deleted). The N10 and G3 tests use new interfaces and cannot compile there; X2, the moves and the renames are
judged by identity.

**Gates and runs.** Release and Debug: fast gate 15/15 (`scratch/logs/phase3_release_gate1.log`,
`phase3_debug_gate.log`, Debug re-configured for the moves); `simple_test_exec test=units` 152 groups, 14,880
assertions, all passed in both (`phase3_units_release1.log`, `phase3_units_debug.log`); `check_guide_imports`,
`check_test_registry`, `check_descr` and `check_cli_chain` pass; no file mode changed (the three "delete mode"
lines of `git diff --summary` are the moves). The `library` label passes 7 of 7 on the final build (222 s, `phase3_library_final.log`), as it did after X2, C1, V2
and N10 alone (`phase3_library1.log`).

The workflow entries of 8.2 (`SIMPLE_SEED=20260923`) on this phase's Release build, against the phase 2 runs,
by `simple_test_exec test=identity` (`scratch/logs/phase3_identity_p3_vs_p2.log` and the per-entry reports):
- **`simulated_workflow_1jxy`:** 202 files identical. The 8 project files are equal after paths are
  normalised. The 6 density files `rho_recvol_state01_part*_{even,odd}.mrc` of solve3D and its two-part run are
  64 bytes longer, starting with `SMPLRHO1`, and from byte 65 on identical to phase 2's (`cmp -i 0:64`).
- **`simulated_workflow_6vxx`:** 174 identical, 7 projects equal after normalisation, and its 2 density files
  differ by the header only, as for 1jxy.
- **`mini_stream_1jxy`:** images, boxes and STAR files identical; the rest is paths, process ids and timings.

So the renames, the moves, X2, C1, N10 and G3 are byte-identical in the workflows, and V2 changes only the density
files' header. The numbers equal phase 2's:
- 1jxy runs to the end, with whole-volume correlation 0.6410 and, over two parts, 0.6379 and an FSC of 0.143 at
  35.66 Å (13:06);
- 6vxx: correlation 0.5185, FSC 49.92 Å (12:51);
- mini_stream passes (0:59).

These timings ran beside a Debug build and the tests (load 2 to 10). The runs are kept as
`scratch/keep/phase3_runs` (the phase 2 runs, superseded, were deleted; `scratch/keep/summary.tsv`).

**Found by the diff review and not changed** (a read-only reviewer checked every item; no defect):
- Frozen sets (`solve3D_addon`) written before V2 pass their manifest check, whose schema version was not
  raised, and then stop with the "no header" error, where trailing chains are re-seeded through their manifest
  version. This happens only when the binary changes in the middle of a run.
- The binary reader indexes particle rows by their absolute row and string rows relative to the header's first
  row. The two agree whenever that first row is 1, as in every project file.

**Exit items**
- X2, C1, V2, N10, O11, V5 done; G3 done for the state of `strategies/search`, the rest recorded (above);
  renames of 3.5 done in the five folders the table allows; glossary written and linked; experiment run and
  recorded, nothing of it in the repository.
- The identity comparison of 8.2 against phase 2: done above. Renames and moves are byte-identical, and the only
  difference is V2's header.
- Fast gate 15/15 in Release and Debug.
- Unit suites of the touched areas: all of `test=units`, both builds, plus the `library` label.
- `check_guide_imports`, `check_test_registry`, `check_descr`: pass.

**Guides updated**: `simple-architecture` (the glossary), `simple-main-project` (C1), `simple-main-volume` and
`references/lifecycles.md` (V2, V5), `simple-main-flex` (N10), `simple-main-strategies` (G3),
`simple-main-strategies-parallelization` (X2), `simple-main-strategies-search` (O11), `simple-main-motion` (V5).

**File table**: rows added for `simple_rec3D_strategy.f90`, `simple_commanders_rec.f90` and
`refine_motion_model_policy.md` (X2), `x_execution_layer_adoption.md` and `flex/fit/simple_flex_pca_embed.f90`
(3.5); widened: `src/fileio/` (C1), `class/` (G3, unused in the end).

Files changed (101; `scratch/logs/phase3_files.txt`):
- X2: `strategies/parallelization/simple_rec3D_strategy.f90`, `simple_refine_motion_model_strategy.f90`,
  `commanders/simple/simple_commanders_rec.f90`, `doc/policies/refine_motion_model_policy.md`.
- C1: `src/fileio/simple_binoris.f90`, `project/simple_sp_project_io.f90` (+ tester).
- V2: `volume/simple_reconstructor.f90`, `simple_frozen_accum.f90` (+ tester), `simple_trail_chain_manifest.f90`,
  `strategies/parallelization/simple_volassemble_strategy.f90`.
- N10: `strategies/parallelization/simple_rec3D_service.f90` (+ tester), `flex/run/simple_flex_pca_state_service.f90`.
- Moves: `strategies/search/simple_pspecs.f90` (from `preprocess/`); `motion/simple_dose_series_manifest.f90` and
  its tester (from `volume/`).
- G3: `strategies/search/simple_strategy3D_alloc.f90`, `simple_strategy2D_alloc.f90`,
  `simple_matcher_ptcl_batch.f90`, `simple_matcher_ptcl_io.f90`, `simple_matcher_3Drec.f90`,
  `simple_matcher_pftc_prep.f90`, `simple_strategy3D_inplane_tester.f90`.
- Renames: 55 files of `src/fileio`, `src/main/image`, `src/main/ori`, `src/main/pftc` and `src/utils`, plus
  `flex/fit/simple_flex_pca_embed.f90` (`scratch/logs/phase3_rename.log`). Also `simple_qsys_job_record.f90` (a
  comment) and `simple_ft_expanded.f90` (a message).
- Docs: `doc/glossary.md` (new); `doc/code_overview/fortran-indexes/*` (regenerated);
  `doc/refactoring_notes/planned/next_release/x_execution_layer_adoption.md`; this plan; the guides above.

### Phase 4: measurements (2026-11-05, the dell)

The four larger changes of section 7 were measured on scratch builds that carry each trial change: copies of the
phase 3 tree under `scratch/p4/`, built in Release. The repository keeps no trial code, and this entry is the
phase's only repository change. Four read-only investigations mapped each change and its measurement first. The
trial patches are kept as diffs against that tree in `scratch/logs/phase4_trial_patches/`, the scripts in
`scratch/logs/phase4_scripts/`, and the logs and scores in `scratch/logs/phase4_*`. The trial trees and builds were
deleted after the numbers were recorded. Timings were taken with the
maintainer's own jobs and this run's single-threaded workflow runs active (load 7 to 18, noted with each).

**I1, class averages deapodized with the stencil envelope.** Deapodization corrects a reconstruction for the
interpolation kernel's attenuation. The class-sum insertion now splats onto the native `box_crop` lattice. The
class averager (`cavger_restore_cavgs`) still corrects with the continuous kernel's instrument function at the
padded period (`prep2D_inv_instrfun4mul(ldim_crop, ldim_croppd, ...)`), which the old padded gather needed.
- **Trial** (patch `i1.diff`): a new `prep2D_inv_kbenvelope4mul`, the inverse of the three-tap stencil's envelope
  at the `box_crop` period, as the 3D gridding reconstructor uses (`kb_stencil_inv_envelope_1d`), called in its
  place.
- **Run:** solve2D on the extracted particles of the 1jxy workflow run (80 particles, `ncls=4`, `mskdiam=180`, one
  thread), current build against the trial build, at nine seeds (`scratch/logs/phase4_i1_seeds.txt`, script `seeds.sh`, profile `radial.py`; 8.5 s a run).
- **Amplitudes:** they change as predicted. The final class averages' radial amplitude, relative to the centre in
  eight steps to half the box, is 1, 0.92, 0.61, 0.22, 0.14, 0.17, 0.09, 0.08 now and 1, 0.94, 0.69, 0.30, 0.19,
  0.25, 0.18, 0.15 with the trial: up to 1.9 times stronger at the edge.
- **Results:** they do not improve, within the spread of this small set. Over eight seeds the median resolution
  of the best class was 8.4 Å now and 10.6 Å with the trial; the median final resolution estimate 40.8 Å and
  50.7 Å. The mean score (a correlation) was 0.445 and 0.463; the trial raises it through the stronger amplitudes
  alone. Both builds left a class empty in two of eight runs. The workflow entries on the trial build
  (logs `scratch/logs/phase4_i1wf_*.log`) moved both ways, from inputs that otherwise match
  phase 3's:
  - 1jxy: whole-volume correlation 0.626, and over two parts 0.628 at an FSC of 49.9 Å, against 0.641 and 0.638
    at 35.7 Å now;
  - 6vxx: 0.555 at 20.8 Å, against 0.519 at 49.9 Å now.

  One run of each cannot separate a change in quality from the spread of de novo runs.
- **Recommendation:** do not adopt on this evidence. The trial is the consistent choice by the kernel's
  arithmetic, but it changes the averages' amplitude profile without a measured gain. Judge it again with a
  two-dimensional metric against the truth (class averages against projections of the simulated volume) on a
  larger simulated set before adopting.

**O1, one icosahedral setting.** SIMPLE builds icosahedral symmetry in SPIDER's setting (a 5-fold axis on z)
everywhere except the volume symmetry search and analysis. `simple_volpft_symsrch` (behind `symaxis_search`,
`symmetrize_map` and `symmetry_test`) and `simple_symanalyzer` pass `icorelion=.true.` and use RELION's setting (a
2-fold axis on z). refine3D and reconstruct3D apply the builder's SPIDER operators.
- **Trial** (patch `o1a.diff`): the three `icorelion=.true.` removed, so everything uses SPIDER's setting.
- **Run** (scripts `o1_run.sh`, `o1_rec.sh`, `blobs.py`; logs `scratch/logs/phase4_o1_*`):
  1. An icosahedral map: two Gaussian blobs off the centre (`blobs.py`), symmetrized with `pgrp=i` into 120 blobs.
     A first try with a low-passed noise volume was useless: icosahedral averaging made it nearly spherical.
  2. Rotated off axis, then aligned by `symaxis_search pgrp=i`.
  3. Checks: which rotations about z leave the aligned map invariant (shell correlations, `fsc`). Then 1,500
     particles projected from it at known poses, reconstructed with `pgrp=c1` and with `pgrp=i`, and the two
     maps compared.
- **Current build:** the aligned map has its 2-fold on z (correlation 0.92 to 1.0 under a 180° turn), and a 72°
  turn decorrelates it like an arbitrary rotation (0.21, 0.11, -0.05 from the fourth shell). The map is not in the
  setting refine3D uses. The `pgrp=i` reconstruction agrees with the `c1` one only to 0.32 to 0.70 in the middle
  shells: the builder's operators average the map wrongly.
- **Trial:** the 72° turn keeps the map (0.94 to 1.0), the 180° turn decorrelates it, and the two reconstructions
  agree to 0.93 to 0.99 in every shell.
- **Note:** production solve3D with `pgrp=i` runs `symmetrize_map` and then `reconstruct3D pgrp=i`
  (`simple_solve3D_utils.f90`), so icosahedral runs are likely affected today. No test covers `pgrp=i`.
- **Recommendation:** adopt one setting, SPIDER's (the trial: three one-line changes in `volume/`). Making
  everything RELION's would also need the asymmetric unit (`build_eullims`, `within_platonic_asymunit`) rewritten.
  Add the reconstruction check above as a library test of `pgrp=i` in the same change.

**O4, L-BFGS-B clients set `factr`/`pgtol`.** In SIMPLE's L-BFGS-B (limited-memory quasi-Newton with bounds),
`factr` and `pgtol` decide when to stop, while `ftol` and `gtol` are the line search's sufficient-decrease and
curvature constants (standard 1e-3 and 0.9). The investigation corrected the premise for one of the two clients.
The shift search (`pftc_shsrch_grad`) already sets `factr`/`pgtol` and never sets `ftol`/`gtol`, so its line search
runs at `opt_spec`'s defaults, 1e-5 for both: nearly exact. The CTF fit (`ctf_estimate_cost`) passes its
tolerances 1e-1 and 1e-3 as `ftol`/`gtol`, and stops at the default `factr`/`pgtol`.
- **Shift search:**
  - Run: the "continuous in-plane" suite (`unit_pftc_align2D3D`), with the optimiser's evaluation counts printed
    (scratch tester change only, patch `o4base.diff`; logs `scratch/logs/phase4_o4_shift_{base,b}.log`).
  - Trial B (patch `o4b.diff`) sets the standard line-search constants. It made the searches worse:
    - evaluations 150 → 271 over the nine searches;
    - five of nine end in "abnormal termination in line search" instead of converging on the gradient test;
    - shift errors against the truth 0.025/0.036/0.056 px → 0.080/0.036/0.188 px.
  - The strict line search is part of what makes this search work.
- **CTF fit:** 26 simulated micrographs of the 1jxy and 6vxx workflow runs, re-estimated per build and compared
  with the simulation's truth (scripts `ctf_run.sh`, `ctf_score.py`; `scratch/logs/phase4_o4_ctf_timing.txt`, `phase4_o4_ctf_score_*.txt`; patches `o4a/o4b/o4c.diff`):

| Build | mean defocus error | astigmatism error | angle error | CTF correlation | time (two sets) |
| --- | --- | --- | --- | --- | --- |
| current | 5.6 nm | 19.6 nm | 0.79° | 0.1594 | 16.0 s, repeat 16.7 s |
| A: literal (`factr`=0.1/ε, `pgtol`=1e-3, standard line search) | 11.4 nm | 14.7 nm | 1.51° | 0.1561 | 13.2 s |
| B: low accuracy (`factr`=1e12, `pgtol`=1e-3, standard line search) | 5.4 nm | 19.4 nm | 1.12° | 0.1593 | 12.7 s, repeat 13.3 s |
| C: standard line search only | 5.5 nm | 19.6 nm | 0.78° | 0.1594 | 14.3 s |

- **Recommendation:**
  - For the CTF fit, adopt C: standard line-search constants, default stopping. It keeps the accuracy and is
    about 10% faster. B is 20% faster but costs angle accuracy; A loses accuracy.
  - For the shift search, keep its line search strict. Make that explicit (`ftol`/`gtol` = 1e-5 at the call) rather
    than leaving it to the defaults, and revise the opt guide's rule: not every client should hand `ftol`/`gtol`
    back to the line search's standard values.
  - The motion correction's shift search (`ftexp_shsrch`, `ftol=gtol=1e-4`) was not measured.

**S4, batch optics grouping moved to the stream's.** No data set on this machine has per-micrograph beam
shifts. Both paths read them from EPU XML files (`read_epu_beamshift`), and there is none under `/data` or
`/data2`. The two STAR files with `rlnBeamTiltX/Y` carry a single global value from RELION's refinement. The sets
with EPU file names (`/data/SAGA-polii`, FlhB) carry no shifts, and with zero shifts both paths give one group per
EPU `Data_` group. The real-data comparison waits for an EPU session with its XML files (only the metadata, tens of
MB, is needed).
- **Synthetic comparison instead:** a scratch tester (patch `s4.diff`, log `scratch/logs/phase4_s4_tester.log`) ran the batch clustering (`h_clust`) and the
  stream's single linkage (`assign_optics_groups`) on the same shifts, five patterns at thresholds 0.02, 0.05 and
  0.1 (`scratch/logs/phase4_s4_synthetic.txt`, `s4_score.py`, scored by adjusted Rand index against the truth):
  - **Well separated** (a 3×3 grid at spacing 0.15, a hexagon at 0.12): identical and exact up to threshold 0.05.
  - **Stream chains:**
    - a slow drift becomes one group 1.2 wide;
    - a 5×5 grid at spacing 0.08 gives 23 groups at the default 0.05 (two neighbours merged, a group 0.106 wide),
      and one group at 0.1;
    - the hexagon gives 3 groups at 0.1.
  - **Batch fragments:**
    - two jittered positions 0.06 apart become 14 or 8 groups;
    - the drift becomes 174 to 195 groups;
    - the 5×5 grid becomes 35 groups at 0.1, more than the 25 positions, because its one-pass centroid merge is
      not transitive;
    - its groups stay about the threshold wide.
- **Recommendation:** do not move the batch to the stream's grouping on this evidence alone. Single linkage
  matches the hierarchical clustering on well-separated multi-shot patterns, is faster and needs no distance
  matrix (which at 20,000 micrographs in one group is 1.6 GB), but it chains neighbouring positions when their
  spacing approaches the threshold. If it is adopted, add a check that warns when a group is wider than twice the
  threshold, and compare both on a real EPU session first.

**Exit items**
- I1, O1, O4, S4 each measured on scratch builds, with numbers and a recommendation (S4 on synthetic data, with
  the reason).
- No repository code changes: `git status` shows only this plan changed in this phase.
- No floor, tolerance or assertion touched; the trial trees and their builds were deleted after the numbers
  were recorded.

Files changed: this plan (this entry).

### Phase 5: the report (2026-11-05, the dell)

- **The report:** `doc/refactoring_notes/completed/agent_friendly_code_report.md`. It covers what was learned,
  what changed phase by phase, the decisions taken under delegation, the deviations from this plan, the evidence
  against the baseline, the measurements with their recommendations, what is left open for the maintainer, and a
  commit message.
- **The driver's per-phase diffs are incomplete.** `/home/elmlundho/agent_runs/agent_friendly/review/phase_N.diff`
  hold only the changes to `.github/skills` and the generated indexes; `phase_4.diff` is empty. The complete
  per-phase diffs were rebuilt from the driver's phase-start snapshots (`state/phase_N.snapshot.json`,
  `state/phase_N.start/`) and HEAD, by `scratch/keep/full_diffs/make_full_diffs.py`. They are in
  `scratch/keep/full_diffs/phase_N.diff`, with per-file line counts in `phase_N.stat`. Files per phase: 19, 27,
  102, 101 and 1 (the plan).
- **The move:** this plan moved from `planned/` to `completed/` with `mv`. Its one link into `planned/` was
  corrected, and its status line points to the report. In `doc/refactoring_notes/README.md` the row moved from
  open work to the completed records.
- **Exit items:**
  - the report written;
  - the index updated;
  - the note moved;
  - fast gate 15/15 in Release and Debug and the static checks after the move
    (`scratch/logs/phase5_release_gate.log`, `phase5_debug_gate.log`).
- **Files changed:** this plan (moved), the report (new), `doc/refactoring_notes/README.md`.

### Phase 6: finishing (2026-11-05, the dell)

**Ruled 2026-11-05 (Hans, recorded in section 8.1 and 8.3):** every recommendation of the report is accepted and
nothing in bounds is deferred. `src/main/sieve/` came into bounds. The items that need one of `ms_posterior2`'s
files go to the follow-up run `agent_friendly2`. Four read-only investigations (the slot producers, the
remaining X1 and G3 sites, the sieve, the small items) mapped the work first.

**Items**

- **O1, one icosahedral setting.**
  - What changed: SPIDER's setting (a 5-fold axis on z) is now used everywhere. The three `icorelion=.true.`
    calls are gone (`simple_volpft_symsrch.f90`, behind `symaxis_search`; `simple_symanalyzer.f90`, behind
    `symmetrize_map` and `symmetry_test`). With no caller left, `sym%new` lost its optional `icorelion` and the
    RELION operator table `make_i_relion` was deleted, as the change itself made it unused.
  - Test: `test_symaxis_search_i` (library suite "volume programs"). A turned icosahedral phantom aligned by
    `symaxis_search pgrp=i` must be invariant under the builder's operators: its self-correlation under them must
    reach the suite's existing control ceiling (0.8), while the turned input's stays below it. Guide: ori.
  - Why the ceiling and not the 0.95 floor of the D2 case: the icosahedral axis search is not exact. The aligned
    map scored 0.914 in one build and 0.983 in the other, and the two values swapped between runs; the old setting
    scores 0.495.
  - The phantom: 120 Gaussian blobs from two seeds at 28 and 32 Å. Three geometries were tried first: two seeds
    at 24 to 30 Å (control 0.83, too symmetric already), two at 33 and 40 Å (aligned 0.939), and one seed
    (control 0.99).
- **O4, the CTF fit's line search.** The CTF fit (`ctf_estimate_cost`) now passes the standard line-search
  constants (`ftol=1e-3`, `gtol=0.9`) and the default stopping tests: phase 4's variant C, which kept the
  accuracy. The polar shift search states its strict line search at both calls (`ftol=gtol=1e-5`, the defaults it
  ran on). The opt guide states the rule and the measured exception.
- **X1, finished.** An `eo=-1` particle that is used stops (`oris%get_half`) in:
  - the iterative (PCG) reconstruction backend's half selections (`simple_rec3D_pcg_strategy.f90`);
  - the reconstruction's state halves (`rec_state_half`);
  - batch scoring (a serial pass before both parallel regions; the half is now taken from `get_eo == 0`
    instead of `eo <= 0`);
  - the class expansion's split;
  - the flex gateway's selection;
  - the GPU-offload reconstructor (a serial pass; this module builds only with `USE_OPENMP_OFFLOAD`, so it was
    not compiled here).

  Rejected particles (state 0) are not looked at. Test: `test_rec_state_half_stops_on_eo` (in-plane strategy
  suite): a used particle with `eo=-1` stops, and a rejected one is skipped. Guides: search, flex.
- **G3, finished in bounds.** Per-thread state is now sized by the team in use: the polar calculator's
  `heap_vars`, memo workspace and CTF work arrays; the class averager's work slices and images; the image
  filters' thread sums; `corrmat`'s and the 2D utilities' search objects (automatic arrays made allocatable);
  `pspecs`; `calc_pspec`'s thread sums. Only the allocation sizes and loop bounds changed, not the parallel
  regions' code. `roavg` keeps `num_threads(nthr_glob)`, which forces its own team to match its arrays.
  `check_pass_team` stays, because the probability tables (`probabilistic/`, out of bounds) still size by
  `nthr_glob`. The existing forced-team testers (`test_search_threads`, `test_posterior_threads`,
  `test_kernel_threads`, which compare teams bit for bit) pass in both builds. Guides: strategies, pftc.
- **The particle-row slot stop (3.1).** Per the ruling:
  - `frame`, written by `trajectory_solve3D` and read back after a project write by the trajectory test, is now
    slot 56. The record keeps its on-disk width of 64 reals, so projects of earlier releases still read.
  - The keys of the other producers are stripped before the write. `import_particles` (`add_single_stk`,
    `add_stktab_2`) strips the CTF constants and any text-oritab keys. The STAR import (`finalise_ptcl_import`)
    strips the optics constants, after moving them to the stacks.
  - The binary-orientation tester's dropped-key case became the stop test.
  - New: `ori%nonslot_key` and `strip_nonslot_keys`, with `oris` wrappers.
  - The stop: `check_particle_slots` in both particle writers of `simple_binoris`, compiled in Debug builds only
    (`_DEBUG`).
  - Tests: `test_nonslot_key_on_particle_row` (forked: Debug stops, Release writes without the key), the slot
    pins moved to 56 and 57-64, a `frame` round trip, and the ori tester's `frame` pins. The Debug unit suites
    and the Debug `library` label tripped no other producer. Guides: ori, star.
- **The sieve (Option 1 of the report for the sieve's developer).**
  - `cleanup_chunk` keeps `chunk_mics.txt` and `consumed_coarse.txt`, so `rebuild_chunked_mics` repairs
    `chunked_mics.txt` after a crash for the whole run.
  - The rule of 3.1 that the included records stay a prefix of the record list is a stop
    (`check_included_prefix`, at the entry of `generate_chunks_coarse` and in `write_chunked_mics`).
  - Tests: `test_rejected_chunk_keeps_its_lists` (a rejection keeps both lists; the rebuild after it lists the
    next chunk) and `test_included_prefix_stop` (forked), in suite "particle sieve";
    `test_restore_after_rejected_chunk` in the p05 tester, suite "particle sieving": the restart after a rejected
    chunk and a crash restores both sets, with nothing imported again. The give-up case's comment now names the
    chunks it still covers.
  - Policies: `ptcl_sieve_policy.md` sections 4.7, 6.1 and 8.2; `restart_policy.md`'s state table. Guide: sieve.
- **`production/stream_ui_contract.json`** was regenerated (`update_stream_ui_contract`): the descriptor
  reference had been stale since before the run. `validate_stream_ui_contract` now passes; the launcher contract
  is unchanged.
- **The frozen-set schema version** is 4, so a context written before the density header is refused cleanly.
  The tester's fixture is version 4, and version 3 is refused.
- **The known issues of the report:**
  - **Stale Cartesian slot:** a slot unprepared by new spectra is now marked stale, and any use of it stops
    (`ptcl_is_valid`). A degenerate particle keeps its handling. Test: a forked child evaluates a stale slot and
    must stop by its check.
  - **Atomic rename:** `simple_rename` no longer deletes the target first. `rename(2)` replaces it atomically, so
    the file is never absent. Its `overwrite` argument now has no effect; it was removed from the 13 in-bounds
    callers. It stays in the signature only because `simple_commanders_solve3D.f90`, a file of `ms_posterior2`
    that a regex had touched and that was restored from HEAD, still passes it; `agent_friendly2` removes both.
  - **CTF memo:** `ft_maps_match` compares the sampling as well as the box. Test: an `apply_ctf` at another
    sampling stops.

**Rows added to the file table:** `solve/simple_solve3D_manifest.f90`,
`commanders/test/simple_commanders_test_highlevel_backends.f90` and `commanders/simple/simple_commanders_solve3D_chunks.f90`
(their `simple_rename` calls); `strategies/parallelization/simple_rec3D_pcg_strategy.f90` and
`simple_cls_expansion_strategy.f90` (X1); `simple_calc_pspec_strategy.f90` (G3);
`commanders/simple/simple_commanders_volops_tester.f90` (O1).

**The new checks on the previous phase's tree** (the base sources plus the driver's phase-6 start snapshot, with
this phase's testers; logs `scratch/logs/phase6_p5tree_*.log`; the scratch build was deleted). The checks that
fail before the fixes:
- the sieve: 4 (the two lists kept, the rebuild after a rejection, the prefix stop);
- the CTF memo's sampling: 2;
- the frozen-set schema: 3;
- X1 in `rec_state_half`: 1;
- the stale Cartesian slot: 2 (the first version of that child evaluated the slot, which the old code already
  refused, so it now asks `ptcl_is_valid` as the consumers do);
- O1: 1 (the old setting's aligned map scores 0.495).

The p05 restart case passes on the old tree by construction: it plants the lists itself and documents the
restart contract that the kept lists make true. The slot stop, the frame slot, the stale-slot flag and G3 use new
interfaces.

**Gates, suites and runs.**
- Fast gate 15/15 in Release and Debug: `scratch/logs/phase6_release_gate_final.log`, `phase6_debug_gate_final.log`.
- `test=units` 152 groups, 14,910 assertions, all passed in both builds: `phase6_units_release_final.log`,
  `phase6_units_debug_final.log`.
- `check_guide_imports`, `check_test_registry`, `check_descr`, `check_cli_chain` pass; no file mode changed.
- The `library` label in Release: 7/7 (248 s, `phase6_library_release.log`). The first run collided with a
  rebuild of the test binary and is superseded.
- In Debug, where the slot stop is active, the `library` label tripped it nowhere. Its one failure was the O1
  test's first phantom (`phase6_library_debug.log`), and `lib_reconstruction` passes in both builds with the final
  test (`phase6_volops_*.log`).

The workflow entries of 8.2 (`SIMPLE_SEED=20260923`) on this phase's Release build, against the phase 3 runs
(`scratch/keep/phase3_runs`):
- **The phase's build** (logs `scratch/logs/phase6_run_*.log`, `phase6_identity_p6_vs_p3.log`): O4 changes the
  CTF fits, so everything downstream differs.
  - mini_stream's defocus errors move by 0.1 nm or less per micrograph.
  - 1jxy: whole-volume correlation 0.645, and over two parts 0.631 at an FSC of 49.9 Å (phase 3: 0.641, and
    0.638 at 35.7 Å).
  - 6vxx: 0.509 at 49.9 Å (phase 3: 0.519). Both remain below the 0.80 floor, as at the base.
- **Attribution:** the same tree with O4 reverted (`phase6_identity_attr_noO4.log`, `phase6_attr_noO4_*.log`).
  - 1jxy: identical to phase 3 (208 files, and the 8 projects after normalisation).
  - 6vxx: identical (176 files; the 7 projects equal token for token, 4 of them with line-only differences from
    path padding).
  - mini_stream: paths and timings only.

  So the X1 stops, the slot stripping, G3, the rename, the CTF memo check and the sieve fix leave the workflow
  outputs unchanged, and O4 accounts for every difference.

**Open after this phase** (each needs a file of `ms_posterior2`'s table, so it goes to `agent_friendly2`):
- R4/X4, C3, C4, R6 (3.4), the judging runs of 3.7;
- the 3.1 check that workers verify the master's sample (the refine3D master);
- `check_pass_team` and the probability tables' `nthr_glob` sizing (`probabilistic/`);
- X2's remaining commander import in `simple_refine3D_strategy.f90`;
- the `overwrite` argument still passed by `simple_commanders_solve3D.f90`.

Recorded:
- the icosahedral `symaxis_search` is not deterministic run to run, and not exact (0.91 to 0.98 self-correlation
  after alignment);
- I1 and S4 wait for 3.7, by ruling.

**Exit items**
- Every item of 8.3's phase 6 is done with its check, or listed above with the out-of-bounds file it needs.
- The gates and identity comparison of 8.2 against phase 3's tree: done; the differences are explained by O4.
- Guides and the report brought up to date; the report has a Phase 6 section and a new commit message.

Files changed: 72 (`scratch/logs/phase6_files.txt` and `review/phase_6.diff`). The groups:
- the sieve and its testers;
- `simple_binoris`, `simple_defs_ori`, `ori`/`oris` (the slot), the project and STAR producers;
- the X1 and G3 sites;
- `simple_syslib` and the `simple_rename` callers;
- the Cartesian calculator, `simple_memoize_ft_maps`, the image CTF code;
- `simple_frozen_accum`;
- the O1 and O4 sources;
- the testers named above, the volops tester;
- `production/stream_ui_contract.json`;
- the two policies, eleven guides, the refactoring-notes index, this plan and the report.
