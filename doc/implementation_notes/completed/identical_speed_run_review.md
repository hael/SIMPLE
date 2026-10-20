# Identical speed changes: review against the SIMPLE standards

Review, 2026-10-18, of the uncommitted refactoring made by the run
`speed_mac` (2026-10-17/18) on the maintainer's Mac. It covers every change
in the run's diff (29 files, about 5,000 lines) and judges the architecture
against the standards of the SIMPLE library. The changes this leads to are
proposed separately, in `speed_architecture_proposals.md`. The run's own
account is `identical_speed_run_report.md`; the plan is
`identical_speed_run.md`. Both are part of the diff.

The standards used are `CLAUDE.md`, the skills `simple-modern-fortran`,
`simple-architecture` and the subsystem skills of the touched areas,
`doc/policies/compile_time_policy.md`, `doc/policies/test_environment_policy.md`,
and the maintainer's standing rules (no `flag .and. array(i)` in one
condition, no goto, no integer-overflow tricks, no `//&` inside the THROW
macros, dead code deleted, short comments without history).

Line numbers refer to the changed files in the worktree
`~/src/SIMPLE-wt-speed`. Nothing was compiled or run for this review; it is
a reading of the code, made by one reviewer per area and checked against the
sources.

## Verdict

The refactoring is sound and almost all of it meets the standards. Every
change made on 2026-10-17 left the outputs of the one-thread identity cases
byte-identical on the Mac. The two changes made on 2026-10-18 (volume
assembly handing its project on, and the kept projection labels) were checked
by the unit tests only, by ruling, and still need the test machine.

Three things do not meet the standards:

- **New entry points on commanders.** `prob_align`, volume assembly and
  `commit_sigma2` gained typed entry points beside `execute(cline)`, called
  by the `refine3D` strategy with in-memory projects. SIMPLE has no such
  pattern. The precedent for the same need (the PCG assembly) puts the body
  in a strategy routine that takes `(params, build)`. The new entry points
  work, but they make commanders into library functions.
- **The kept reprojection builder reaches into the builder's internals.** It
  repeats parts of `build_spproj` and `build_general_tbox` by hand. Two of the
  run's own defects (reused killed volumes, and a write into the wrong
  project) came from exactly this.
- **Peak memory on the master grows by one or two full projects** in the
  distributed iteration, because volume assembly copies its whole project out
  before it frees its reconstructors, and the copy lives until the sigma
  commit. On a 750,000-particle project that is about 1 GB per copy (two
  such projects took 1.95 GB in the run's copy timing); the master's peak
  itself was not measured.

Two defects matter: a truncated dense table part file is read without an
error, and the dense file no longer lets the master check which reference a
row belongs to. A handful of smaller findings (one dead public routine, two
dead per-thread buffers, an unchecked "load before evaluate" contract, stale
skills, stale comments) are listed below. None of them affects the identity
results.

## What the refactoring is

Terms used throughout:

- **Identity test**: the same case, at one thread with a fixed seed, on the
  build before and the build after a change; the declared outputs must be
  byte for byte the same.
- **Master and workers**: in distributed `refine3D` the master process
  schedules worker processes per partition of the particles; between
  iterations it reads their results, builds the probability tables'
  assignment, assembles the volumes and filters them.
- **Dense table**: the `prob` table, which holds a distance for every
  particle and every reference. **Neighbourhood table**: the `prob_neigh`
  table, which holds a few candidates per particle. Each worker writes its
  part of the table to a **part file**; the master reads all parts.
- **NU filter**: the nonuniform filter applied in volume assembly.
- **V**: the memory of one padded volume: 16.6 MB at box 160, 109 MB at box
  300.

The run carried out the identical items of the speed week in four areas.

| area | what changed | identity evidence |
| --- | --- | --- |
| probability tables | neighbourhood table read in one pass (counts per particle in the part file); compact record for the dense table (8 bytes per entry instead of 28 in the file and 24 in memory, plus a list of refined shifts); parallel read of the dense parts with a coverage check | one-thread cases, Mac |
| polar kernels | shift matrix converted only on the searched shells; repeated all-rotation objective reused inside the shift search; one double-precision copy of the reference per minimisation; per-thread buffers in the score loop, which removed heap trampolines and per-call allocations | one-thread cases and exact tester checks at 1 and 4 threads, Mac |
| volume assembly and NU filter | timers inside the filter and per state; two unused volume reads removed; the filter's input halves kept in memory; the objective computed right after each candidate; the composition filters the input pair again, so the filter writes no scratch files | one-thread cases, Mac; the composition timed on PfCRT half maps, box 300: 10.5 to 9.9 s |
| the master's project handling | trailing reread skipped; matching low-pass read from one record; reprojection builder kept for the run; `prob_align` given a copy of the master's project; volume assembly given the master's project and handing its own on to `commit_sigma2` and the low-pass refresh; trailing counts taken from memory, which keeps back-filled projection labels | the first four by one-thread cases on the Mac; the last two by the unit tests only |

On the identity cases (87 particles, box 192) wall times moved by about 2%,
the neighbourhood table's refinement by about 9%. The changes target large
runs; their gains are for the test machine to measure.

## The changes by area, with the standards

### Probability tables

**Architecture.**

- Two new types in `simple_eul_prob_tab_utils.f90`: `prob_dense_entry`
  (distance and in-plane index, 8 bytes, no default values, line 35) and
  `prob_refined_entry` (particle, reference, shift, line 42), and a format
  word `PROB_TAB_DENSE_FORMAT = -2` (line 31).
- The master's dense table is now three arrays: `rows(nrefs,nptcls)`, the
  refined entries in particle order, and offsets per particle (a compressed
  sparse row layout). `dense_candidate` rebuilds the full candidate for the
  assignment.
- The dense part file is: header (nrefs, nptcls, -2, number of refined
  entries), particle indices, the seed table, the rows, the refined entries.
  The neighbourhood part file keeps its stream of candidates and gains the
  candidate count of every particle after the particle indices.
- The sparse writer moved into the neighbourhood type, which now overrides
  `new_worker`, `begin_write`, `write_tab` and `read_tabs_to_glob`. The dense
  type writes rows per batch at the particles' positions.
- The master reads the dense parts in an OpenMP loop (line 913). After the
  loop it checks that every particle is in exactly one part and that the parts
  agree on the rotation grid (lines 919-924).
- No reader of the old 3D layout remains. The 2D table keeps its own format.

**Standards.** Imports, lifecycle (`rows`, `refined` and the offsets are
freed in `kill`), dead code (`fill_tab`, the old reader and chunk constant
deleted), nested conditions, THROW macros and comment length all meet the
rules. The new tester checks run in the fast suite `unit_pftc_align2D3D`; no
registration was needed. Deviations, all small:

- The format word is public in the utilities module, but only the dense
  table uses it; it belongs to the dense table as a private parameter.
- `write_rows` and `dense_candidate` are public bindings that the
  neighbourhood type inherits, where they are meaningless. `write_rows` there
  would overwrite the sparse file's header.
- The testers import `omp_lib` without `only:` (tester line 18) and change
  the global thread count with `omp_set_num_threads(4)` instead of a
  `num_threads` clause, which the test policy prefers.
- A stale comment at line 265 still speaks of a candidate buffer per thread,
  and lines 262-267 carry a run measurement (older).

**Defects.**

1. **A truncated dense part file is read silently.** The rows and refined
   entries are read with `iostat` and checked with `fileiochk` (lines 977,
   996, 1001), which only prints on end of file and does not stop
   (`simple_fileio.f90`, lines 71-72). Since the master's table is no longer
   initialised, a short file leaves undefined distances that drive the
   assignment, and the coverage check passes because a particle is counted
   before its row is read. Fix: stop on any nonzero `iostat`, and compare the
   file size with the size the header implies.
2. **The reference check is lost.** A row's index now means "reference
   number among the active states" on the worker that wrote it. The master
   checks only the number of references. The old format stored the full
   reference index, which the master mapped and checked. If master and
   worker disagreed on the active states with the same count, rows would be
   misread silently. Fix: write the active state list after the header and
   compare it on read.
3. **THROW inside the parallel read.** A part that fails stops the program
   from a worker thread; two failing at once race in the error path. Fix: a
   status per part, raised after the loop.
4. **Integer kinds.** File offsets are 64-bit; the refined offsets and counts
   are default integers. About 3e8 refined entries at 750,000 particles and
   7,500 references with the default angular threshold is safe; a much wider
   threshold or more directions would pass 2^31 without a guard. Fix: 64-bit
   offsets or the size guard already used in `new_fixed_candidate_store`.
5. The dense and the neighbourhood readers treat a particle missing from all
   parts differently (stop against a silent default candidate). One rule
   should hold for both.

### Polar kernels

**Architecture.**

- New per-thread state in `heap_vars` (`simple_polarft_calc.f90`, lines 45,
  50-52): a double-precision reference buffer for one shift search and three
  score, distance and cumulative buffers of length `nrots`. They are
  allocated in `new` and freed in `kill` for every thread.
- New state in the shift-search object (`simple_pftc_shsrch_grad.f90`,
  lines 61-64): the last all-rotation objective and its shift. It is active
  only inside `minimize` (lines 402-403, cleared at 484-485) and reused when
  the shift matches bit for bit (line 337).
- New public routines: `load_ref_aln_8` and the two evaluators
  `gen_corr_for_rot_8_aln` and `gen_corr_grad_for_rot_8_aln`, which read the
  loaded reference. The coarse searches became private. A sampler over an
  array of distances (`sample_likelihood_dists`) replaced internal
  procedures in the score loop.
- `gen_corr_grad_only_for_rot_8` and `shift_ref` are deleted; no reference
  remains.

**Standards.** The submodule structure is respected (interfaces in the
parent, bodies and private helpers in `simple_polarft_corr`). Lifecycle,
nested conditions, the bitwise shift comparison (`transfer` to integers,
portable) and comment lengths meet the rules. Deviations:

- `gen_corr_grad_for_rot_8` is now dead in production: its only caller was
  the shift search, and only the new tester calls it (tester line 924). Under
  the dead-code rule it goes; the tester can compare against a fresh load.
- Two older per-thread buffers next to the new ones, `kcorrs` and
  `pft_r1_8` (`simple_polarft_calc.f90`, lines 42 and 48), are allocated and
  freed but never used. They should go with this change; that also returns
  the memory the new buffer adds.
- The tester section header names the run ("exact checks of the identical
  speed changes", tester line 856), which is history in code.
- A bare `use omp_lib` in the in-plane tester (line 21).

**Contracts and risks.**

- The new evaluators read whatever reference the calling thread last
  loaded. Only `minimize` and the tester follow this today. A future caller
  that evaluates without loading gets silent garbage. Fix: keep the loaded
  reference's index and half in `heap_vars` and stop on a mismatch (one
  integer compare per evaluation).
- The loaded state now lives for a whole minimisation instead of one
  evaluation, so a nested parallel region or an untied task inside a search
  would corrupt it. Nothing does that today; the pftc skill should say it.
- The bench line "prob_tab objective evaluations" changed meaning (it no
  longer counts cache hits and now counts the two evaluations of `minimize`
  itself), so it cannot be compared with older reports.

### Volume assembly and the nonuniform filter

**Architecture.**

- Four wall-time accumulators in `simple_nu_filter_vars.f90` (line 71), reset
  at the start of each filter sequence and read by `get_nu_filter_timings`.
  The volume-assembly report gains a section of split timings (filter parts,
  and the restore and filter time of every state) outside its accounted
  total; the summary script shows them.
- The call sequence keeps its shape (set up, optimise, compose, clean up),
  but `nu_filter_vols` now takes a required input pair
  (`simple_nu_filter.f90`, line 471) and filters it again per label. **The
  caller must keep the input pair alive and unchanged from setup to the
  composition.** Every caller does: volume assembly, the PCG master,
  `nu_filt3D`; `postprocess_nu` never composes, and `flex_pca` runs without
  the filter.
- The scratch cache files (`nu_filter_cache_*`), their names, writes and
  deletion, and the unused `nu_filter_vol` are deleted.
- In volume assembly the filter's input halves are kept from the restore
  (`keep_half`, `simple_commanders_rec_distr.f90`, lines 721-729) instead of
  being read back from the files just written. The trailing counts use the
  particle field in memory (decided on 2026-10-18), which keeps back-filled
  projection labels.

**Standards.** Module state sits where the compile-time policy wants it (the
shared companion module), imports and lifecycle meet the rules, the
multi-state memory policy still holds (one state's volumes at a time), and
the policy documents were updated. Deviations:

- The skills are stale (not updated, by design; the run's report proposes
  text). `simple-main-nu-filt` shows `nu_filter_vols` without the input pair
  (an example that would not compile), still speaks of scratch cache files,
  and should state the new rule that the input pair lives until the
  composition. `simple-nonuniform-regularization` and its reference map name
  the cache files.
- The header of `compose_nu_labels_from_volume` has four lines (limit three).
- The planned note `next_release/calpha_nu_evidence_mask.md` still calls the
  deleted `nu_filter_vol`.

**Risks.**

- **Setup memory grows by two volumes** (the two transforms live through the
  candidate loop): +33 MB at box 160, +217 MB at box 300. The run's "same
  peak" measurement was made with `nu_filt3D`, not with volume assembly under
  ML regularisation, where setup is the peak. Two identical savings offset it:
  the auxiliary pair is held three times during setup (the caller's copy, the
  filter's copy and a stash), and the PCG solvent path keeps the base pair
  through a composition that reads the solvent pair.
- The composition and the candidates are filtered by two code paths that
  must stay bit-equal (a per-label Butterworth filter against the stored
  bank). Only the tester guards it. Composing through the same candidate
  routine would make it hold by construction.
- `keep_half` sizes its copy from `params` rather than from the half it
  copies; correct today, fragile if the grid changes.

### The master's project handling

**Architecture.**

- **New entry points.** `commander_prob_align%execute_from_project` and the
  same on `commander_prob_align_neigh` (`simple_commanders_prob.f90`, lines
  23, 29, 267-272), `commander_volassemble%execute_hand_on` (lines 21,
  837-843) and `commander_commit_sigma2%execute_from_project`
  (`simple_commanders_euclid.f90`, lines 20, 53-58). Each and the old
  `execute` call one private body (`run_prob_align`, `run_volassemble`,
  `run_commit_sigma2`), which builds `params` once.
- **The builder** takes an optional source project on `build_spproj` and
  `init_params_and_build_general_tbox` (`simple_builder.f90`, lines 130-139,
  166-205) and copies it instead of reading the file, then sets the project's
  working directory as before.
- **The kept reprojection builder** is an allocatable component of the
  distributed strategy (`simple_refine3D_strategy.f90`, line 106), made on
  the first iteration and released in clean-up. Later iterations read only
  the segments the model needs, remake its volumes, write the particle field
  and release the project (lines 572-609). Between iterations it holds the
  toolboxes, no volumes and no project.
- **Project reads per distributed iteration** (gridding, `prob`, Euclidean
  objective, NU filter): five full reads before, one now (the master's own);
  the particle field once by the kept builder; up to four full copies in
  memory. Writes are unchanged. At 751,729 particles a read took 1.3 s and a
  copy 0.28 s on the Mac, so the saving is about 4 s per iteration at that
  size.

**Standards.** `params` is built once per path and the command line is set
before it; the builder's optional source is a true option (most callers omit
it); the kept builder is allocatable and released, as the compile-time
policy asks; no new module-level imports. Deviations:

- **Layering.** In SIMPLE a commander is a command object reached through
  `execute(cline)`. The new entry points make three commanders into library
  calls with typed arguments, and leave about 520 lines of assembly workflow
  (`run_volassemble`, lines 845-1369) in a commander module, now called from a
  strategy. The repository's own precedent for this need is
  `execute_rec3D_pcg_distributed_master(params, build, ...)` in
  `strategies/parallelization`, called from the strategy with the master's
  own builder. The standard-conforming shape is a strategy or domain routine
  that both `execute(cline)` and the strategy call.
- **Builder internals in the strategy.** The kept builder rereads segments,
  sets the working directory and remakes volumes by hand (lines 584-595),
  repeating `build_spproj` and `build_general_tbox`. Lines 590-591 set the
  working directory that `update_projinfo` has just set (to the same value).
  A builder method (or
  materialising on the master's own builder, as the shared-memory path does)
  would own this.
- **One-record read in a strategy.** `read_segment` ignores its record range
  for project files, so the low-pass refresh opens the binary file itself.
  That puts project-file knowledge in a strategy; the fix belongs in the
  project layer.
- **Heavy local.** `vol_proj` is a plain local `type(sp_project)` with a
  separate flag (line 1319); allocatable would drop the flag and free it
  early.
- **Dead release.** The release of the kept builder in `distr_initialize`
  (line 1085) can never fire, because the strategy is created fresh for every
  run; its comment describes a case that cannot happen.

**Risks.**

- **Peak memory** (above): volume assembly copies its whole project out
  (`simple_commanders_rec_distr.f90`, line 903) before `cleanup_context`
  frees its reconstructors, and the copy lives through post-processing until
  the sigma commit, when the master holds three projects. Volume assembly
  changes only `res`, `res05`, `cfar`, `lp`, `lp_est` and `proj`; the
  refresh needs one value and the commit needs fields the master already has.
  Handing back only those, or copying after the clean-up and freeing early,
  removes the cost.
- **"The copy equals the file" is unchecked.** The review traced every write
  between the master's read and each consumer and found the condition holds
  today. A future write in between would make the copies diverge without any
  error. A single owner of the project (below) removes the question.
- **The final full write carries stale fields**: the master's particle field
  never receives the last assembly's `cfar`, `lp_est` or back-filled `proj`.
  This is as before; adopting the handed-on field would fix it, but it changes
  the final project.
- The PCG backend and `refine=sigma` still do full reads in the refresh and
  the commit; the master's own builder could serve both.

## Standards at a glance

| standard | verdict |
| --- | --- |
| `use ..., only:`; procedure-only imports inside procedures | met in production code; two testers import `omp_lib` bare |
| no umbrella modules, no new CMake targets, no per-file flags | met |
| heavy types allocatable as components | met (the kept builder); one heavy plain local (`vol_proj`) |
| `new`/`kill` symmetry | met in every new type and buffer |
| encapsulation of stateful modules | the NU filter remains module state (unchanged); new public bindings inherited where meaningless (`write_rows`) |
| `optional` only where callers omit it | met |
| `params` built once per commander path | met |
| layering `ui -> exec -> commander -> strategy/domain` | **deviates**: typed entry points on three commanders, called from a strategy |
| orchestration in commanders/strategies, numerics in domain modules | met; one project-file read moved into a strategy |
| no `flag .and. array(i)`; no goto; no overflow tricks; THROW macro rule | met |
| dead code deleted | mostly; `gen_corr_grad_for_rot_8`, `kcorrs`, `pft_r1_8`, the dead builder release remain |
| comments short, no history or measurements | mostly; one tester header names the run, one comment stale |
| tests in tester modules of registered suites | met; testers change the global thread count |
| skills and policy documents current | policies yes; skills `simple-main-nu-filt` and `simple-nonuniform-regularization` stale |
| memory policy for multi-state assembly | met; master peak memory and NU setup memory grow (not measured) |

## What follows from this review

The actions this review leads to (the small fixes that finish the change,
the question of the commander entry points, and the architectural steps for
performance) are in a separate proposal note,
`speed_architecture_proposals.md`, so that they can be ranked apart from
what was done.

## What the test machine still has to show

- Identity on the x86 build, which the Mac's ARM build does not prove: the
  simulated workflow of beta-galactosidase (both steps), the three-state
  `refine3D`, the correlation-objective and the PCG cases at one thread, base
  against final.
- For the 2026-10-18 changes: volume assembly's hand-on byte-identical on
  those cases; for the kept labels, which callers reach the back-filled
  labels and that the only differences are the labels and what follows from
  them.
- Peak memory of the master in a distributed `refine3D` at large particle
  counts, before and after (the handed-on project, the NU setup).
- Wall-time gains on the cases of the runs `prob_speed` and `ms_speed`.

## Files of the diff

Probability tables: `simple_eul_prob_tab.f90`, `simple_eul_prob_tab_neigh.f90`,
`simple_eul_prob_tab_utils.f90`, `simple_eul_prob_tab_tester.f90`,
`simple_commanders_prob.f90`. Polar kernels: `simple_polarft_calc.f90`,
`simple_polarft_core.f90`, `simple_polarft_corr.f90`, `simple_polarft_geom.f90`,
`simple_pftc_shsrch_grad.f90`, `simple_pftc_inplane_tester.f90`. Volume
assembly and the NU filter: `simple_nu_filter.f90`, `_vars`, `_state`, `_bank`,
`_apply`, `_tester`, `simple_nu_state_filter.f90`,
`simple_commanders_rec_distr.f90`, `simple_commanders_resolest.f90`,
`scripts/plot_refine3d_bench.py`. The master: `simple_refine3D_strategy.f90`,
`simple_builder.f90`, `simple_commanders_euclid.f90`. Documents:
`refine3D_policy.md`, `nonuniform_filtering_policy.md`,
`single_state_speed_week_completed.md`, `multistate_speed_run_completed.md`, and the run's plan
and report.
