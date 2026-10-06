# Stream 3D snapshots plan, 6 October 2026

Particle snapshots from one or more selected states of stream multistate 3D (p07), while the stream
runs; state selections once it has finished; and a 2D snapshot told apart from a 3D one throughout.
The SIMPLE side follows the 2D snapshot path (NICE → master → p06 → snapshot report). The NICE side
reuses the batch 3D viewer's state selector and its state-selection job.

## Decisions

1. **A 3D snapshot merges the selected states into one.** Particles of the selected states become
   state 1, keeping their 3D orientations; all others are deselected. The state artifacts
   (volumes, FSCs) leave the project's `out` segment. This is the result the batch
   `selection oritype=ptcl3D states=...` gives, ready for a new refine3D.
2. **The files are those of a 2D snapshot, plus the selected volumes.** Besides the project:
   - the micrograph and ptcl2D STAR files, with the optics map and the optics offset, as for 2D;
   - each selected state's volume copied into the snapshot folder, as files the project does not
     register (it has one merged state).
3. **After the stream has finished, selections use the batch state selection.** NICE launches the
   existing `selection oritype=ptcl3D states=...` job on p07's final project, as the batch 3D
   viewer does. The result is a batch job, not a stream particle set. SIMPLE needs no new program.
4. **NICE names the kinds `snapshot2D` and `snapshot3D`.** A particle set's type says which; NICE
   maps each type to its stage folder (`classification_2D`, `solve3D_multistate`). Existing
   `snapshot` sets become `snapshot2D` through a NICE data migration. Release 4 keeps no
   compatibility.

Defaults taken without a question:
- **The GUI answer key** is `snapshot3D` {`id`, `selection` (state numbers), `filename`}, beside
  `snapshot2D`, with the same rule for the name: a bare `*.simple` name of at most 128 characters.
- **A snapshot takes p07's current result,** the latest run (solve3D, an add-on pass or the final
  refine3D). p07 keeps no history of runs. A request before the first result is answered as not
  written (0 particles, no file), as p06 answers one for an iteration it no longer keeps.
- **The optics table** is the one p07's project holds, which is the newest publication's. The
  offset is p06's: `max(nicedispid - 1, 0) * OPTICS_ID_DELTA`.
- **The copied volumes** are each selected state's raw volume, `vol_state<NN>.mrc` with the
  original state number.

## 1. Where things stand

- **2D snapshots:**
  1. NICE's `snapshot_classification_2D` adds a particle set of type `snapshot` and puts
     `snapshot2D` in `master_update`; `get_master_update` sends it once.
  2. The master parses it into `gui_metadata_stream_update`, which goes to p06.
  3. p06's `write_snapshot` writes `classification_2D/snapshots/<stem>/<name>`, with its class
     averages, FRCs and STAR files.
  4. p06 answers with a `gui_metadata_stream_pool2D_snapshot` and the snapshot's class-average
     tiles, which the assembler nests as `pool2D.snapshot`.
  5. NICE fills the particle set from that report, and links it into batch jobs from
     `classification_2D/snapshots/...`. That path is hard-coded in `stream_views.py`
     (`view_stream_link_particle_set`), `job_builder_views.py` (`_snapshot_source`) and
     `batch_views.py`.
- **p07 reads no GUI updates.** The master sends updates to p01, p03 and p06 only
  (`l_updates` in p00; IPC policy §5.2). p07's project is
  `solve3D_multistate/solve3D_multistate.simple`, written as p07 finishes (`write_stage_project`).
- **Batch 3D:** `_cls3D_state_selector.html` toggles states and "save selection" posts them to
  `batch_cls3D_selection`. `BatchJob.createStateSelection` then launches `selection
  oritype=ptcl3D states=...`, or `state=` for a single state. The stream multistate page includes
  the same selector, but its stream branch has no selection controls.

## 2. SIMPLE

### 2.1 The request

1. **The update type** (`gui_metadata_stream_update`) gets:
   - `snapshot3D_id` and `snapshot3D_filename`, the latter of `MAX_SNAPSHOT2D_FNAME_LEN`
     characters, renamed `MAX_SNAPSHOT_FNAME_LEN` for both;
   - `snapshot3D_selection(MAX_STATES_SOLVE3D_MULTISTATE)` with its length;
   - `set_snapshot3D_update`, `get_snapshot3D_update` and `has_snapshot3D_update`.

   The update type imports `MAX_STATES_SOLVE3D_MULTISTATE`, a metadata module, so the layering
   rule holds.
2. **The master's parser** (GUICMD) reads `snapshot3D` with the same name rule as `snapshot2D`.
   It drops a selection that is empty, longer than `MAX_STATES_SOLVE3D_MULTISTATE`, or holds a
   state outside 1..`MAX_STATES_SOLVE3D_MULTISTATE`, with a warning.
3. **p00** sets `l_updates` for `STAGE_SOLVE3D`, so p07 receives updates. Only a stage that reads
   updates gets them, and the same answer then reaches p06, p07 and the others: each stage acts on
   its own fields, as today.

### 2.2 p07 writes the snapshot

1. **`apply_gui_updates`.** p07 drains its update pipe once per pass, as p06 does, and calls
   `write_snapshot` for a new `snapshot3D` id. Each id is written once (`last_snapshot_id`), and
   the request never interrupts a running job.
2. **`write_snapshot(update)`:**
   1. Without a result yet, answer "not written".
   2. Otherwise, copy the stage project. Set state 1 on the particles of the selected states in
      `ptcl3D` and `ptcl2D` alike, so the ptcl2D STAR matches, and state 0 on the rest.
   3. Remove the state artifacts from `out` (`remove_state_artifacts_from_osout`, as the batch
      selection does).
   4. Write the project to `snapshots/<stem>/<name>`.
   5. Write `<stem>_micrographs.star` and `<stem>_particles.star` with the optics offset.
   6. Copy each selected state's volume as `vol_state<NN>.mrc`.
   7. Report back.
3. **Shared steps.** The optics offset and the two STAR calls are p06's too. A small shared helper,
   which writes a snapshot project's STAR files, could serve both. It would live in
   `src/main/stream/shared`, and p06 would move onto it in the same change.
4. **Restart.** A restarted p07 starts its snapshot ids from zero. NICE sends each request once,
   so no request repeats. `restart_policy.md` gets a row for p07's `snapshots/`.

### 2.3 The report

1. **One snapshot type for both kinds.** `gui_metadata_stream_pool2D_snapshot` becomes
   `gui_metadata_stream_snapshot` (module `simple_gui_metadata_stream_snapshot`). It gains the
   selected states (`MAX_STATES_SOLVE3D_MULTISTATE`, empty for 2D), emitted as `states` when set.
   A new tag, `GUI_METADATA_STREAM_SOLVE3D_SNAPSHOT_TYPE`, marks the 3D report; the 2D report
   keeps its tag.
2. **The master's meta store** keeps the latest 3D snapshot report, which `clear_stage` drops with
   p07's lists.
3. **The assembler** nests the 3D report as `solve3D_multistate.snapshot`, as the 2D one goes in
   `pool2D.snapshot`. The section's hash covers it, as it covers everything else.

### 2.4 After the stream has finished

There is no new SIMPLE program; NICE runs the existing `selection` on p07's final project
(decision 3). That needs two checks:
- the final project's stack and volume paths must resolve from another directory, which means
  absolute paths. If they don't, p07 writes its final project with absolute paths;
- `selection` must accept a project whose `ptcl2D` and `ptcl3D` both carry the stream's rows.

## 3. NICE

1. **Types.** Rename `snapshot` to `snapshot2D` wherever particle sets are made, read or linked:
   - `streamjob.snapshot_classification_2D`;
   - the zoom template's `particle_set.type == "snapshot"`;
   - `view_stream_link_particle_set`'s type set;
   - `_snapshot_source`'s type check.

   One map from type to stage folder replaces the three hard-coded `classification_2D` paths:
   `snapshot2D` → `classification_2D`, `snapshot3D` → `solve3D_multistate`. A data migration
   renames the type in existing jobs' `particle_sets_stats`.
2. **Live 3D snapshots:**
   - `streamjob.snapshot_solve3D(selected_states)` adds a `snapshot3D` particle set (name
     "particle set N", file `snapshot_<id>.simple`) and puts `snapshot3D` in `master_update`;
   - `get_master_update` pops `snapshot3D` after sending, as it does `snapshot2D`;
   - `update_stats` fills the set from `stats_json["solve3D_multistate"]["snapshot"]`: particle
     count, time, file name, states;
   - a view `view_stream_snapshot_solve3D` takes the posted states;
   - in the stream branch of `_cls3D_state_selector.html` (the multistate zoom page), a
     "snapshot" button posts the kept states while the job runs. It works like the batch
     branch's "save selection": tiles toggle, and at least one state stays selected.
3. **After the stream has finished:** the same selector posts to a view that calls
   `BatchJob.createStateSelection(project, workspace, <stream>/solve3D_multistate/solve3D_multistate.simple,
   states)`. The page shows this button instead of "snapshot" once the job is no longer running.
4. **Particle sets of type `snapshot3D`** show their states (with the state tiles the page already
   has) where 2D sets show class averages, and link into batch jobs like 2D sets, from
   `solve3D_multistate/snapshots/...`.
5. **Tests:** the migration, the new views, the type-to-folder map, and `_snapshot_source` for
   both types.

## 4. Docs

- **IPC policy:** §5.2 (p07 reads updates: `snapshot3D`), the §6 table (`snapshot3D` and its
  rules), §6.3 (dropped 3D selections), and §8 tests.
- **`stream_3D_ingestion_policy.md`,** p07's policy: a section on 3D snapshots covering what they
  hold, when they are written, and the "not written" answer.
- **`restart_policy.md`:** p07's `snapshots/`.
- **`pool2D_policy.md`:** the snapshot report's new type name.
- **The `simple-main-stream` skill,** with approval: p07 reads updates.

## 5. Tests (SIMPLE)

- **`unit_stream` "stream master":**
  - `snapshot3D` parsing: valid, an unsafe name, an empty or oversized selection, a state out of
    range;
  - p07 receiving updates;
  - the store keeping the 3D report.
- **`unit_stream` "solve 3D":** `write_snapshot` on a fixture project with three states, checking:
  - the merged state and the deselected particles;
  - `out` without the state artifacts;
  - the copied volumes and the STAR files;
  - the report;
  - before a result, the "not written" answer.
- **`unit_ui`:**
  - the snapshot type with and without states (serialisation and JSON);
  - the multistate section nesting `snapshot`.

## 6. Order

1. SIMPLE 2.1 to 2.3, with their tests and docs. The stream then writes 3D snapshots that NICE does
   not ask for yet; the new key is ignored until NICE sends it.
2. NICE 3.1, the types and folder map with the migration: 2D keeps working under its new type.
3. NICE 3.2 and 3.4, live 3D snapshots.
4. SIMPLE 2.4's checks, then NICE 3.3, selections after the stream has finished.

## Status

Implemented, 6 October 2026: SIMPLE 2.1 to 2.3, NICE 3.1 to 3.5, and the docs of section 4,
including the `simple-main-stream` skill line. Nothing has been compiled or run yet, in SIMPLE or
NICE.
The plan moves to `completed/` once `unit_ui`, `unit_stream` and the NICE tests pass.

Where the implementation differs from the plan:
- **The snapshot's source** is `result_projfile`, the project of p07's latest finished run, which
  `finish_run` records. The stage project is written only as p07 finishes, so it can't be the
  source. `ptcl2D` is merged only when it has as many rows as `ptcl3D`.
- **No shared STAR helper (2.2.3).** p07 calls the two STAR writers directly, as p06 does.
  `OPTICS_ID_DELTA` moved from the p06 stage to `simple_defs_stream`, so both use one offset.
- **`pool2D_policy.md`** needed no change, because it does not name the report type.
- **An existing 2D bug was fixed.** The stage's report overwrote a set's file name with the
  absolute path, so the set could no longer be linked. An unwritten report also left an empty file
  name. `apply_snapshot_report` now keeps the bare name, stores the path as `path`, and drops the
  name when nothing was written. This applies to both kinds.
- **Batch jobs find a stream snapshot's folder** from the set's type in the stream job's
  `particle_sets_stats`. A `snapshot:<job>:<id>` source does not carry the type.
- **2.4's checks are still open:**
  - the stream's final multistate project resolves from the selection job's folder;
  - `selection` accepts a project with both `ptcl2D` and `ptcl3D` rows.

  Both need a finished stream. NICE checks only that the project exists before it launches the
  job.
