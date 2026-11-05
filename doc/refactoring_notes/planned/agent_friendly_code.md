# Agent-friendly code

Status, 2026-11-04: on hold until the batch alignment pipeline refactoring
([plan](../../implementation_notes/planned/prob_align_batch_pipeline.md)) has landed, since both
change the same code. Delivered with this note: section 2 and the defect fixes of section 5 (1-7).
Every ruling is decided; the work they call for is section 7, small items first.

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
([plan](../../implementation_notes/planned/prob_align_batch_pipeline.md)) rewrites
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

