# SIMPLE GUI Onboarding Guide

This note is for a new developer starting on the batch side of the SIMPLE GUI.
The immediate goal is to understand the existing NICE GUI, make the batch GUI
more robust, and keep its communication model aligned with the newer stream
work.

## Scope

SIMPLE has three related but distinct GUI concerns:

- NICE: the Django-based web GUI under `nice/`.
- Batch, also called classic in the NICE code: discrete `simple_exec` or
  `single_exec` jobs launched from a workspace.
- Stream: automated streaming jobs launched through `simple_stream`.
- Desktop: a robust local application experience around NICE and the SIMPLE
  binaries.

For the batch GUI work, focus first on the classic path. Use stream as the
reference implementation for more structured job communication, status, and
control.

## First Mental Model

The GUI is not a separate scientific engine. It is a launcher, monitor, and
project browser around the Fortran executables.

```text
Fortran UI definitions
        |
        | simple_private_exec prg=print_ui_json
        v
Django/NICE form rendering
        |
        | create JobClassicModel and job.script
        v
simple_exec or single_exec
        |
        | updates workspace.simple, logs, thumbnails, and GUI status
        v
NICE project/job views
```

The stream path has the same broad idea, but the execution model is a
multi-stage streaming pipeline with richer status and update messages.

## Repository Map

Start with these files and directories:

| Area | Path | Why it matters |
| --- | --- | --- |
| NICE build | `nice/CMakeLists.txt` | Optional Django/Tailwind build, venv setup, migrations, static assets. |
| Local launcher | `nice/nice_local.py.in` | Current local GUI bootstrap; useful but not yet a robust desktop app. |
| Django models | `nice/nice_lite/models.py` | Defines projects, datasets, workspaces, stream jobs, batch jobs, dispatch templates. |
| URL routing | `nice/nice_lite/urls.py` | Separates global, classic, stream, and API endpoints. |
| Classic views | `nice/nice_lite/views_classic.py` | Web entry points for batch workspaces and jobs. |
| Classic job object | `nice/nice_lite/data_structures/jobclassic.py` | Owns batch job creation, selection jobs, status updates, deletion. |
| SIMPLE launcher | `nice/nice_lite/data_structures/simple.py` | Bridges Django to `simple_exec`, `single_exec`, and `simple_stream`. |
| Classic job form | `nice/nice_lite/app_views/newjobview.py` | Builds the new-job form from Fortran-generated UI JSON. |
| Classic job panels | `nice/nice_lite/app_views/jobview.py` | Reads `workspace.simple` and renders project stats, micrographs, classes, and logs. |
| Classic templates | `nice/nice_lite/templates/nice_classic/` | HTML for batch workspace, job, logs, micrograph, and class-average panels. |
| Classic JS | `nice/nice_lite/static/nice_classic/` | Client-side behavior for batch views. |
| Stream templates | `nice/nice_lite/templates/nice_stream/` | More developed stream monitoring UI. |
| Stream JS | `nice/nice_lite/static/nice_stream/` | Useful reference for richer interactive controls. |
| UI definitions | `src/main/ui/` | Fortran-side definitions of programs and input parameters. |
| SIMPLE UI group | `src/main/ui/simple_ui_simple_group.f90` | Registers non-stream `simple_exec` GUI programs. |
| Stream UI group | `src/main/ui/simple_ui_stream_group.f90` | Registers stream GUI programs. |
| UI JSON printer | `src/main/ui/simple_ui.f90` | Implements `print_ui_json` and stream UI JSON output. |
| Private exec | `production/simple_private_exec.f90` | Dispatches `print_ui_json` and `print_ui_stream`. |
| GUI metadata | `src/utils/gui/metadata/` | Typed Fortran metadata objects that serialize to JSON. |
| GUI assembler | `src/utils/gui/simple_gui_assembler.f90` | Newer stream-side JSON assembly with change detection. |
| Stream communicator | `src/utils/comm/simple_stream_communicator.f90` | HTTP communicator used by stream stages. |
| Legacy NICE comm | `src/utils/gui/simple_nice.f90` | Older socket/thread communication object used by several batch commanders. |
| GUI tests | `production/tests/simple_test_gui_metadata.f90` and `production/tests/simple_test_gui_assembler.f90` | Fortran-side GUI metadata tests. |
| Django tests | `nice/nice_lite/test/` | Existing model and data-structure tests. |

Also inspect `nice/nice_lite_dev/`. It appears to contain newer stream-facing
work, including `StreamJob`, master heartbeat fields, and a more explicit
stream master launch path. Before porting ideas, decide whether `nice_lite` or
`nice_lite_dev` is the active development base for your branch.

## Batch GUI Flow

The batch side is called classic in the Django code.

1. User chooses a workspace in the classic UI.
2. User opens a new job page.
3. `NewJobView` asks `SIMPLE` for UI JSON.
4. `SIMPLE.loadUIJSON()` runs:

   ```text
   simple_private_exec prg=print_ui_json
   ```

5. The form is built from the JSON sections for the selected program.
6. `JobClassic.new()` creates a `JobClassicModel`, assigns a job directory, and
   records user arguments.
7. `SIMPLE.dispatch()` writes `job.script`.
8. The dispatch script copies the parent `workspace.simple`, runs
   `simple_exec prg=update_project`, then executes the requested program.
9. The job writes `stdout.log`, `stderr.log`, and an updated `workspace.simple`.
10. Job views read `workspace.simple` using:

   ```text
   simple_exec prg=print_project_info json=yes
   simple_exec prg=print_project_field json=yes
   ```

11. The GUI renders global project stats, micrograph tables, 2D class-average
   panels, plots, histograms, selections, and logs.

Important batch files:

- `nice/nice_lite/data_structures/jobclassic.py`
- `nice/nice_lite/data_structures/simple.py`
- `nice/nice_lite/app_views/newjobview.py`
- `nice/nice_lite/app_views/jobview.py`
- `nice/nice_lite/templates/nice_classic/`
- `src/main/commanders/simple/`
- `src/main/ui/simple/`

## Stream GUI Flow

The stream side exists both as the installed `nice_lite` version and as newer
work in `nice_lite_dev`.

The current installed `nice_lite` stream launcher starts several independent
stream process scripts based on `simple_private_exec prg=print_ui_stream`.
The newer `nice_lite_dev` direction launches a stream master:

```text
simple_stream prg=master
```

That master owns stream subprocesses, assembles GUI metadata, sends a single
versioned JSON payload, and receives update commands from NICE. This is the
communication direction worth learning from.

Important stream files:

- `nice/nice_lite/data_structures/job.py`
- `nice/nice_lite/data_structures/simple.py`
- `nice/nice_lite/views_stream.py`
- `nice/nice_lite/templates/nice_stream/`
- `nice/nice_lite_dev/data_structures/streamjob.py`
- `nice/nice_lite_dev/data_structures/simple.py`
- `src/main/stream/simple_stream_p00_master.f90`
- `src/utils/gui/simple_gui_assembler.f90`
- `src/utils/comm/simple_stream_communicator.f90`
- `src/utils/gui/metadata/stream/`

## Communication Model To Align

The batch and stream GUIs should eventually share the same concepts:

- job id
- schema version
- status
- stage
- timestamp
- process id
- stats payload
- update command payload
- termination request

Stream already points in this direction. Batch currently has lighter-weight
status handling through `JobClassicModel.status`, `heartbeat`, and `update`,
plus older `simple_nice_comm` usage in selected commanders.

A useful target envelope is:

```json
{
  "version": 1,
  "jobid": 123,
  "kind": "batch",
  "heartbeat": {
    "status": "running",
    "stage": "iteration 4",
    "pid": 9999,
    "timestamp": 1779111111
  },
  "stats": {},
  "update": {}
}
```

For stream, `kind` can be `stream` and the heartbeat can contain stage-specific
entries. For batch, it can be a single process heartbeat. The GUI should not
need a completely different vocabulary for each.

## Desktop GUI Direction

The current local entry point is `nice_local`. It syncs NICE into
`~/.nice/nice_4_0`, runs migrations, creates a local user, resets the local
dispatch template, starts Django's development server, and opens a browser.

That is useful for development, but a robust desktop experience needs more:

- deterministic app data location
- safe repeated startup
- dynamic localhost port selection
- clear logs
- supervised server process
- graceful shutdown
- browser/window ownership
- migration and static-asset checks
- packaging of Python runtime expectations
- clear failure messages when SIMPLE binaries or libraries are missing

Treat the first robust desktop milestone as a hardened local NICE runtime. A
native shell such as Tauri, Electron, or pywebview can come later; it should
wrap a stable local runtime rather than hide fragile startup behavior.

## Suggested Reading Order

Day 1:

1. Read this note.
2. Open `nice/nice_lite/models.py`.
3. Compare `JobModel` and `JobClassicModel`.
4. Open `nice/nice_lite/data_structures/simple.py`.
5. Trace `SIMPLE.dispatch()` for batch jobs.
6. Trace `SIMPLEStream.dispatch()` for stream jobs.

Day 2:

1. Open `nice/nice_lite/app_views/newjobview.py`.
2. Run or inspect `simple_private_exec prg=print_ui_json`.
3. Open `src/main/ui/simple_ui.f90`.
4. Open `src/main/ui/simple_ui_program.f90`.
5. Open one program definition under `src/main/ui/simple/`, for example
   `simple_ui_cluster2D.f90`.

Day 3:

1. Open `nice/nice_lite/data_structures/jobclassic.py`.
2. Trace `JobClassic.new()` through job directory creation and parent project
   handling.
3. Open `nice/nice_lite/app_views/jobview.py`.
4. Trace how global stats and field stats are read from `workspace.simple`.
5. Open the classic templates and JavaScript for job, logs, micrographs, and
   2D class panels.

Day 4:

1. Open `nice/nice_lite_dev/data_structures/streamjob.py`.
2. Read `StreamJob.update_stats()`.
3. Open `src/utils/gui/simple_gui_assembler.f90`.
4. Open `src/utils/comm/simple_stream_communicator.f90`.
5. Sketch which ideas should be shared with batch communication.

Day 5:

1. Run the Django tests if the local environment is available.
2. Run the Fortran GUI metadata tests if the build is available.
3. Pick one narrow first ticket from the list below.

## First Good Tickets

These are intentionally small enough to build confidence without changing the
scientific behavior of SIMPLE.

1. Document the difference between `nice_lite` and `nice_lite_dev`, then decide
   which tree should receive active batch GUI work.
2. Add a batch communication contract document and compare it with the stream
   payload shape in `simple_gui_assembler`.
3. Add more explicit heartbeat/status fields to `JobClassicModel`, mirroring
   the useful parts of `JobModel.master_*` from `nice_lite_dev`.
4. Build a small classic API test that posts a batch heartbeat and verifies the
   database status/update response.
5. Harden `nice_local` startup: avoid fragile rename behavior, detect occupied
   ports, write logs, and make repeated launches idempotent.
6. Choose one batch command, such as `import_movies`, and trace its
   `simple_nice_comm` updates from commander to Django.
7. Prototype a generic batch heartbeat using the same status vocabulary as the
   stream master.
8. Add a UI-level smoke test for new batch job form rendering from
   `print_ui_json`.

## Guardrails

- Keep Fortran command-line parameters registered through the typed SIMPLE
  parameter path. Do not invent GUI-only command-line keys.
- Preserve the project-file contract. Batch jobs usually operate on a copied
  `workspace.simple` in the job directory.
- Keep child command lines sparse. Avoid restating defaults unless the child
  behavior actually needs to change.
- Do not let Django become the scientific source of truth. The project file and
  SIMPLE executables remain authoritative for scientific state.
- Prefer a shared GUI communication vocabulary over separate batch and stream
  special cases.
- Add tests around communication and dispatch behavior before changing launch
  semantics.
- Keep desktop packaging separate from scientific workflow changes.

## Open Questions For The Project Lead

- Should active development happen in `nice_lite`, `nice_lite_dev`, or should
  `nice_lite_dev` be merged back first? - Yes. We should move nice_lite_dev to 
  nice_lite once we're ready to start active development. 
- Should the future batch communicator use HTTP like stream, or adapt the older
  `simple_nice_comm` socket/thread model? use HTTP like stream
- Which batch commands are the first-class desktop workflows? import procedures
- What minimum offline/local behavior is required for a desktop release? None - already implemented
- Should the desktop app support local-only execution first, or also configure
  SLURM/LSF submission from day one? Local only to start.
- What is the expected support matrix for macOS and Linux desktop installs? All possible

## Practical Definition Of Done For The First Milestone

A useful first milestone is not a redesigned GUI. It is:

- batch jobs launch reliably from NICE,
- job status transitions are explicit and testable,
- at least one batch command reports structured heartbeat/stage information,
- the GUI can send a small update or terminate command and the job can consume
  it,
- local desktop startup is repeatable and logs failures clearly,
- the stream and batch sides use the same names for common communication
  concepts.

That gives the batch GUI work a stable foundation while keeping it synchronized
with the stream direction.
