# Stream initial 3D from the sieve's class averages, plan, 8 October 2026

An optional route for the stream's first 3D. p06 builds its first publication from the sieve's own
2D: the class averages, FRCs, class table and per-particle 2D parameters of the sets it has
imported, combined. p07 then runs `solve3D_cavgs` on those class averages to make the initial
volumes, maps the classes' poses and states to the particles, and refines at the particle level
with `solve3D`. Later publications, and p07's addon and final runs, stay as they are. The route is
off unless a command-line parameter turns it on.

## Decisions

1. **The trigger is today's.** p06 publishes the initial set at the first iteration that leaves
   `NPTCLS_FIRST3D` (100,000) particles selected by the pool, or at iteration 10, whichever
   comes first. The route changes what the publication holds, not when it is made.
2. **solve3D_cavgs gets every accepted class, combined.** These are the classes the sieve
   accepted in the sets of the initial publication. They are concatenated with their even and
   odd halves and their FRCs. That can be several hundred to a few thousand class averages.
3. **The particle-level run starts from the class averages' poses and states.** Each class's 3D
   orientation and state is mapped to its particles (`map2ptcls`). `solve3D` then runs with
   `cavg_ini_ext=yes`: it starts after the symmetry stage and keeps the class-to-state
   assignment.
4. **Fewer states are accepted.** If `solve3D_cavgs` leaves fewer than `nstates` populated
   states, p07 continues with the populated ones, renumbered from 1. One state left gives a
   single-state run.

Defaults taken without a question:
- **The parameter** is `sieve_ini3D` (yes|no, default no), on the stream master's UI entry. The
  master passes it to p06 only. p07 follows the publication: it takes the route when the first
  publication carries `sieve_ini3D=yes` in its out segment, as it follows `pool_final=yes`. So p06
  and p07 cannot disagree, and a p07 restart that takes that publication takes the route again.
  Change the name before implementing if another reads better.
- **A failed `solve3D_cavgs`** (the job fails or leaves no result) falls back to today's first
  run, `solve3D` on the first set, with a warning. Decision 4 covers a collapse, not a failure.
- **The initial publication holds the sieve's selection,** not the pool's: the particles each
  set accepted, with the sieve's 2D parameters. The first set is then taken as today, at most
  `nptcls3D_max` particles, whole stacks in order, and keeps its selection for the session.

## 1. Where things stand

- **The sieve's sets** (p05 `hand_off`) are copies of finished chunk projects. Each holds the
  chunk's own 2D:
  - in `os_out`, the chunk's `cavg` and `frc2D`, at the native box, since a chunk's `solve2D`
    writes its final class averages there;
  - `cls2D`, whose states are the sieve's class rejection;
  - `ptcl2D`, with chunk-local class labels, in-plane angles, shifts and the sieve's selection.

  `cleanup_chunk` keeps a chunk's final class averages and the last iteration's even and odd
  stacks.
- **p06 drops all of it at import.** `append_project_sets` adds the particles without 2D
  parameters, and the pool draws its own classes. Every publication (`build_pool_publication`)
  carries the pool's classes, class averages and FRCs.
- **p07's first set** is the first publication as the pool selects it (`take_first_publication`),
  at most `nptcls3D_max` particles as whole stacks in order. `take_cavgs` copies the
  publication's main class-average stack and FRCs into `quality_selection/<id>/`, not the even
  and odd halves.
- **solve3D_cavgs** reconstructs `nstates` volumes from the class averages of a project's `cavg`
  entry. Its `cls2D` states select the classes, and it needs the `_even` and `_odd` stacks beside
  the main one. It leaves each class's 3D orientation and state in `cls3D` and the volumes as
  `vol_cavg` entries, and reports `final_nstates`.
- **`map2ptcls`** (`simple_sp_project_ptcl`) composes each `cls3D` orientation with its
  particles' 2D orientations, through `ptcl2D`'s class labels. It sets their `ptcl3D` pose,
  `proj`, `state` and `corr`, and their state in both segments.
- **solve3D `cavg_ini_ext=yes`** needs prior `ptcl3D` poses and, with `nstates > 1`, every state
  populated (`validate_cavg_ini_ext_states`). It skips the symmetry search and starts after the
  symmetry stage. `cavg_ini=yes` does not fit this route, because it runs its class-average 3D
  with one state.

## 2. The parameter

- `sieve_ini3D` is declared in `simple_parameters.f90` as `'no'`, registered in
  `simple_parameters_parse.f90`, and added to the stream master's UI entry (`simple_ui_stream`).
- The master's `make_stage_clines` sets it on p06's command line when it is yes. p07 needs no
  parameter (see the defaults).

## 3. p06: the initial publication from the sieve's 2D

1. **Collecting.** While `sieve_ini3D=yes` and nothing has been published, each imported set's
   class averages (main, `_even`, `_odd`) and FRCs are copied into p06's `sieve_cavgs/<set stem>/`,
   and the set is recorded. The sieve's chunk folders are not relied on at publication time.
   A set without class averages (the sieve's empty final set, a chunk with every class rejected)
   adds no classes.
2. **Building** (`build_sieve_publication`, beside `build_pool_publication` in
   `simple_stream_refine2D_utils`). It runs at the same trigger as the first publication. From
   the recorded sets, re-read from the sieve's completed folder, it builds:
   - the sets' micrographs, stacks and particles, stacks renumbered from 1, each particle with its
     image index in its stack (`indstk`, as `build_pool_publication` sets it);
   - the sieve's 2D parameters, with each set's class labels offset by the classes of the sets
     before it, and the sieve's selection;
   - `cls2D`, the sets' class tables concatenated, the sieve's states kept, each `pop` counted
     from the selected particles;
   - `ptcl3D` prepared as today: a copy of `ptcl2D` with the 2D clustering removed and the shifts
     kept;
   - the newest optics map's groups, as today.
3. **The files.** The combined class averages, their halves and the FRCs are written beside the
   publication under today's names (`<id>_cavgs.mrcs` with `_even` and `_odd`, `<id>_frcs.bin`):
   - the stacks are concatenated in set order;
   - the FRCs are a new `class_frcs` filled class by class from each set's file, which needs one
     box and pixel size across the sets (the native ones; a mismatch is refused).

   `os_out` registers them, with the sieve's mask diameter and `sieve_ini3D=yes`. The project is
   written last, under a temporary name renamed into place, as today. `sieve_cavgs/` is then
   removed.
4. **Afterwards.** Later publications are the pool's, as today, and p06 collects nothing more.
   A restarted pool with publications on disk resumes at iteration 25, as today, so the route
   applies to a fresh pool's first publication only.

## 4. p07: solve3D_cavgs, then solve3D

1. **The first set** is taken as today (`take_first_publication`). On this publication its
   selection is the sieve's.
2. **`take_cavgs` copies the even and odd halves too,** so the stage's `cavg` entry has them. This
   costs little, and `solve3D_cavgs` needs them.
3. **A new phase, `PHASE_CAVGS3D`.** When the first publication carries `sieve_ini3D=yes`,
   `next_job` starts `solve3D_cavgs` instead of `solve3D` on the first set, in `solve3D_cavgs/`:
   - its settings are `nstates` (3), `pgrp=c1`, the mask diameter, `nparts3D` and `nthr3D`;
   - its low-pass limits are solve3D_cavgs' own defaults;
   - it runs with `exit_collapse=no`, so a collapse leaves fewer populated states rather than an
     early exit.
4. **Its result** (`finish_cavgs3D`):
   - the job's `cls3D` and its `vol_cavg` entries are read;
   - the populated states are counted, and with fewer than `nstates` they are renumbered from 1
     (decision 4). The stage's states in use become a new `nstates_eff`, used wherever p07 now
     uses `nstates`: status, volumes, FSCs, the snapshot, solve3D and the final refine3D. The
     addon runs replay theirs from the manifest;
   - `map2ptcls` maps the classes' poses and states to the stage's rows, after which only
     first-set rows stay selected. `map2ptcls` sets a state on every row of a class, so the excess
     the cap deselected must be deselected again.
5. **solve3D** then starts as today's first run (the cap and queue, `balance=none`), with
   `cavg_ini_ext=yes` and `nstates=nstates_eff`.
6. **A failed `solve3D_cavgs`** falls back to today's first run, with a warning (see the
   defaults). A failed `solve3D` stops the stage, as today.
7. **The GUI** shows "running solve3D_cavgs" during the new phase, with the solve3D progress of 1.
   The `nstates` it gets is `nstates_eff`.

## 5. Docs

- **`stream_3D_ingestion_policy.md`:** a section on the route, covering the initial publication
  from the sieve's 2D, the cavgs phase, the state compaction and the fallback. Also its scope,
  the §3 item on when publications are made (what the first holds in this route), and tests.
- **`pool2D_policy.md`:** the publication variant and `sieve_cavgs/`.
- **`restart_policy.md`:** p06 (the route only in a fresh pool), p07 (a restart that takes the
  marked publication takes the route), and `sieve_cavgs/` and `solve3D_cavgs/`.
- **`ipc_policy.md` or the stream README:** the master passes `sieve_ini3D` to p06.
- **The `simple-main-stream` skill,** with approval: the optional route.

## 6. Tests

- **`unit_stream` "pool 2D":** `build_sieve_publication` on fixture sets with small class-average
  stacks and FRC files, checking:
  - the class offsets, the combined `cls2D` with states and populations, and `indstk`;
  - the selection is the sieve's, and the marker is set;
  - the combined files hold every class in set order;
  - a set without class averages adds no classes;
  - a box mismatch is refused.
- **`unit_stream` "solve 3D":**
  - `take_cavgs` copies the halves;
  - the route is taken only on a marked first publication;
  - the state compaction (a pure helper): three states with one empty become two, renumbered;
  - after `map2ptcls`, rows the cap deselected stay deselected;
  - `nstates_eff` drives the status and volume loops.
- **`unit_ui`:** the master's new entry.

## 7. Order

1. The parameter and its forwarding. Behaviour is unchanged.
2. p06's collection and the sieve-based publication, with tests. p07 takes it through today's
   path, which ignores the marker.
3. p07's halves in `take_cavgs`, `PHASE_CAVGS3D`, the mapping, `cavg_ini_ext` and `nstates_eff`,
   with tests.
4. The docs.

## 8. Risks

- **Size.** `solve3D_cavgs` on a few thousand class averages takes time and memory. It can run
  distributed (`nparts`), but it was tuned on smaller sets.
- **Mixed classes.** The chunks' classes are not aligned with one another: the same view appears
  in many chunks' class averages. That is what solve3D_cavgs sees, but it is more redundant than
  one 2D run's classes.
- **Missing halves.** The class-average halves depend on what `cleanup_chunk` keeps; their names
  are checked during implementation.
- **Fewer states change the rest of the session.** `nstates_eff` below `nstates` holds for the
  whole session: addon runs and the final refine3D use it.

## Status

Implemented, 8 October 2026: sections 2 to 6, the skill line included.
Nothing has been compiled or run. The plan moves to `completed/` once `unit_stream` and `unit_ui`
pass and a stream has run the route.

Where the implementation differs from the plan:
- **States are counted on the particles, not the classes.** `take_cavgs3D_result` maps first,
  then renumbers the states the selected particles hold. `cavg_ini_ext` needs every requested
  state populated with particles, and the cap can leave a state that only classes hold.
- **The stage's states** are `params%nstates`, set to the states held, rather than a new
  `nstates_eff`. Every use in p07 follows it.
- **A box or FRC mismatch, or a missing file,** falls back to the pool's first publication with a
  warning, rather than stopping p06.
- **The FRCs** keep the sets' own box, which may differ from the class averages'. Their pixel size
  follows from the class averages'.
- **A first set too small for solve3D:** the marker waits. If a later publication is merged
  before solve3D can start, `solve3D_cavgs` runs on the classes then current (the pool's,
  registered with their class averages by `take_cavgs`).
