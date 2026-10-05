# Stream IPC Policy

How the stream's processes talk: the stages and the master over pipes, the master and the GUI over
HTTP. What a message may contain, what is dropped, and what a change must preserve.

## 1. Scope

- **Pipes and framing:** `src/main/stream/shared/simple_stream_pipe.f90`. The pipe descriptors
  are in `src/main/stream/shared/simple_stream_state.f90`, and the master's stage table in
  `src/main/stream/master/simple_stream_master_stage_ids.f90` and
  `src/main/stream/master/simple_stream_master_stage.f90`.
- **GUI metadata types:** `src/utils/gui/metadata/` (`gui_metadata_base` and its extensions).
- **The master's store, GUI commands and heartbeat:**
  - `src/main/stream/master/simple_stream_master_meta_store.f90`;
  - `src/main/stream/master/simple_stream_master_gui_commands.f90`;
  - `src/utils/gui/simple_gui_assembler.f90`;
  - the heartbeat loop of `src/main/commanders/stream/simple_commanders_stream_p00_master.f90`.

## 2. Processes and pipes

1. The master forks one process per stage. Each stage has two pipes, both made by the master
   before the fork: stage to master (`ipc_pipe_<stage>_in`) and master to stage
   (`ipc_pipe_<stage>_out`).
2. The master keeps every pipe end open for the whole run. A forked stage closes the ends that
   are not its own (`close_other_pipe_ends`) before its commander runs.
3. `stream_pipe` wraps the ends a process uses; it never closes the descriptors (they are the
   master's).
4. A stage is forked only before the master's listener thread starts, or while holding the
   listener's lock (`doc/policies/stream/restart_policy.md`).

## 3. Framing

A frame is a C `int` byte count followed by that many payload bytes.

1. **Send:** retries at once on `EINTR`, and with a 10 ms pause on `EAGAIN` up to 200 times. A
   frame is delivered whole or dropped whole: once a byte is in the pipe, the rest follows.
2. **Receive:** returns a frame already assembled from earlier reads before reading again, so a
   drain loop sees every queued message.
3. **Resync:** a length outside 1..`max_frame_bytes` is reported, and the buffered bytes are
   dropped, so the reader resynchronises on a later frame instead of stopping.
4. **Discard:** drops what the pipe and the buffer hold, for a reader whose writer has been
   replaced (a restarted stage).
5. **Frame bound:** `max_frame_bytes` is `max_metadata_size()`, the largest serialised metadata
   type.

## 4. Stage to master: GUI metadata

1. **A message is one serialised GUI metadata object** whose first bytes are its type
   (`GUI_METADATA_*_TYPE`).
2. **Serialisation is a byte copy** (`transfer`) of the object, so:
   - an object sent over a pipe has no allocated allocatable or pointer component. A type with
     one overrides `serialise` and sends a copy without it (`gui_metadata_vol3D` drops its
     reprojection tiles);
   - sender and receiver must be the same binary, which forked stages are.
3. **Capacities are fixed** by the types (named constants where callers need them:
   `MAX_PICKREFS_SELECTION`, `MAX_SNAPSHOT2D_SELECTION`, `MAX_MIC_COORDINATES`). Callers stay
   within them: the picks drawn on a thumbnail are cut at `MAX_MIC_COORDINATES`. The setters stop
   on overflow as a guard against programming errors; nothing from outside may reach that guard.
4. **A reused metadata object is reset with `kill`**: every type extending `gui_metadata_base`
   overrides `kill` to assign its default-initialised self, so no field of the previous message
   survives into the next. Every component therefore keeps a default initialiser; `new` sets only
   the base flags.
5. **The master keeps the latest message per type** (`stream_master_meta_store`). A frame whose
   length is not its tag's type's (a desynchronised pipe read as a frame) is dropped with a warning
   before it is copied into anything (`frame_fits`).
6. **What it keeps:**
   - a status replaces the previous one;
   - an item of a list (micrograph, optics group, class average, volume, reprojection tile)
     goes to slot `i` of a list of `i_max`, which is remade when `i_max` changes.
7. **Locking:** the master's listener thread drains every stage's pipe into the store while
   holding the metadata lock. The main loop holds the same lock to assemble a heartbeat, and to
   discard a stage's pipes and fork it again.

## 5. Master to stage: updates

1. Each GUI answer becomes one `gui_metadata_stream_update` holding that answer's fields. NICE
   answers every heartbeat with its whole state, so most answers repeat the last one.
2. The update goes to the running stages that read updates, and only when it differs from the
   last update the stage was sent since it was started (`stream_master_stage%send_update`). A
   restarted stage therefore receives the next answer whole:
   - p01 reads the CTF resolution, astigmatism and ice thresholds;
   - p03 reads the picking-reference selection and its cycle;
   - p06 reads the 2D mask diameter and snapshot requests.
3. A stopped stage, and a stage asked to stop, is not sent updates (its pipe would only fill).
4. A stage drains its updates once per pass and applies them in order.
5. **The master's update writer gives up on a part-written frame** after
   `UPDATE_PARTIAL_RETRIES` (about 2 s; `stream_pipe%limit_partial_frames`): a stage that stopped
   reading must not hold the master. The channel is then broken, later updates are dropped, and
   the discard before the stage's next start mends it. The stages' writers to the master have no
   such limit, since the master's listener always drains them.

## 6. Master and GUI

1. **Heartbeat:** every `HEARTBEAT_S` (5 s) the master POSTs the assembler's JSON document to
   `niceserver`.
   - Each content section is sent only when its hash (FNV-1a) differs from the last one sent.
   - A heartbeat the GUI did not accept (no answer, or a status other than 200) clears the
     hashes, so the next heartbeat sends everything.
2. **The answer** (status 200) is a JSON object; text that is not JSON is logged and ignored.
   Recognised keys:

   | Key | Effect |
   |---|---|
   | `terminate` | stop the stream |
   | `terminate_<key>`, `restart_<key>` | stop, or restart when stopped, one stage (`<key>` from `stage_gui_key`); NICE keeps `restart_<key>` until it sees the stage running, so the master acts on it once, again only after the key has left an answer, and never once the stream is stopping |
   | `ctfresthreshold`, `astigthreshold`, `icefracthreshold` | preprocessing thresholds |
   | `pickrefs_selection`, `pickrefs_cycle` | a picking-reference selection of p03's cycle |
   | `mskdiam2D` | the 2D pool's mask diameter |
   | `snapshot2D` {`id`, `iteration`, `selection`, `filename`} | a 2D snapshot request; p06 answers each id once with a snapshot report, where 0 particles and no file mean it was not written (an iteration no longer kept) |

3. **Dropped selections:** a selection larger than the update holds is dropped whole with a
   warning. Acting on part of it would act on classes the user did not choose; the rest of the
   answer is applied.
4. **Nothing in an answer stops the master.**

## 7. Stopping

1. Every stage, and the master, handles SIGTERM (the master also SIGINT) by setting a flag only
   (`simple_stream_sigterm`); the flag is polled between steps.
2. **The master's stop:**
   - optics assignment is asked first, with up to `OPTICS_STOP_TIMEOUT_S` (60 s);
   - then every running stage is asked on each pass;
   - a stage still running after `STOP_TIMEOUT_S` (600 s) is killed (SIGKILL).
3. `forked_process` never restarts a child by itself: a stage restarts only when the GUI asks.

## 8. Change rules

- A new metadata type: no allocatable or pointer component in what it sends (or a `serialise`
  override), and a named capacity for any array a caller fills from outside.
- A new GUI answer key: parsed in `stream_master_gui_commands`, with oversized arrays dropped, and
  added to section 6.
- Any read of a stage pipe on the master's main thread takes the metadata lock.
- Tests:
  - `unit_ipc` "stream pipe": framing, resync, discard, an abandoned part-written frame;
  - `unit_stream` "stream master": stage names and keys; GUI answers, including invalid and
    oversized ones; the store, including a frame of the wrong length; a stage's pipes from both
    sides; an update sent once;
  - `unit_ui` "GUI metadata": every type a stage sends survives serialisation
    (a copy received by transfer holds the same fields).

## 9. Known gaps

- **`gui_metadata_project`** has allocatable components and no `serialise` guard. It is not sent
  over a pipe today.
- **The sieve-reference selection** (`ref_selection`) is still sent by NICE and ignored; the
  update type no longer has a field for it.
- **`increase_nmics`** is still sent by NICE and ignored.
- **The byte-copy format** ties the processes to one binary; a stage started by `exec` would need
  a real wire format.
