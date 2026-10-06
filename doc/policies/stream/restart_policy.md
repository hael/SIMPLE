# Stream Restart Policy

What happens when a stream stage is started again, by the GUI or in a run folder that already
holds its output. What each stage keeps, restores and redoes, and what a change must preserve.

## 1. How a stage is started again

1. **From the GUI**, after the stage has stopped (`restart_<key>` in a heartbeat answer). NICE
   keeps the key in its answers until it sees the stage running, so the master acts on a request
   once, and again only after the key has left an answer; it acts on none once the stream is
   stopping. The master:
   - drops what the stopped stage left in its pipes, both its unread messages and the updates it
     never read, and mends an update channel broken by an abandoned frame (`discard_pipes`);
   - drops the stage's lists from its store (`clear_stage`), so the GUI shows only what the new
     process sends;
   - forks the stage again with the same command line, while holding the listener's lock;
   - sends it the next GUI answer whole, whatever it sent the previous process.
2. **By starting the master again** in a run folder where the stages' folders exist. The master
   forks every stage as at a first start, and each stage finds its folder.
3. **By starting one stage on its own**, in its folder or with `dir_exec=<folder>` (p04, p06).

The master never restarts a stage by itself, and `forked_process` never restarts a child, so a
stage that fails stays failed until the GUI restarts it. A stage the master skipped (p01 with
`dir_preprocess`, p03 with `pickrefs`) stays skipped: a GUI restart of it is ignored and logged,
since its output is the user's earlier run. The master holds no restart state of its own: each
stage decides what it restores.

A master started again with `dir_preprocess` finds the link it made at its first start
(`preprocessing`) and keeps it when it points to the same folder.

## 2. How a stage recognises a restart

A stage recognises a restart by its output folder existing before `params%new` creates it, or by
`dir_exec` where the stage takes one. It reads the command line for this, before the parameters
exist, and then removes `dir_exec` from the command lines it passes on.

## 3. Per stage

| Stage | On restart |
|---|---|
| p01 preprocessing | Recognised by the execution folder existing (`outdir`, or `dir_exec`). A stop cancels the jobs in flight. Its job sets are restored: every movie of the completed sets, accepted or rejected, goes into the watcher's history, import indices continue after the highest given, and the accepted micrographs are imported again. Set numbering continues after the highest set, completed or left unfinished; the job folder with the unfinished sets is set aside (`<folder>_unfinished<k>`) and their movies submitted again. The thresholds the GUI set come back from `gui_thresholds.txt`, and a generated gain reference is reused. The micrograph STAR file is rewritten. A leftover `TERM_STREAM` is removed. |
| p02 optics assignment | The micrograph segment starts empty, and every completed upstream project is imported again (the watcher's history starts empty). The micrographs the newest map in its folder lists come back with that map's group ids (`latest_optics_map_table`), and no map is published until all of them are imported again, or a pass finds nothing more to import (a warning then names how many did not come back): the later stages keep applying the newest map, which is complete, instead of a partial one that would put the rest in group 1. Then the groups are assigned again keeping their ids, and new ones continue after the highest; optics-map ids continue after the newest map. A leftover `TERM_STREAM` is removed. |
| p03 initial analysis | A stop cancels its extraction, 2D and 3D jobs. When the picking references are published, they are sent to the GUI again and the stage is finished at once (`doc/policies/stream/reference_generation_policy.md`). Otherwise the previous run's working folders are cleared (micrograph copies, picks, extractions, the sieve's chunks and hand-offs, the 2D and 3D runs, the selections, the balancing), and the plan starts again from cycle 1 on every completed upstream project. p03 uses no `TERM_STREAM`. |
| p04 reference picking | A stop cancels the jobs in flight. Its job sets are restored: the completed sets are imported again, and the upstream projects they were made from go into the watcher's history, so none is picked twice. Set numbering continues after the highest set, completed or left unfinished; the job folder with the unfinished sets is set aside and their inputs submitted again. A leftover `TERM_STREAM` is removed. The picking templates are made again from the published references when the first new set arrives. |
| p05 particle sieving | A stop cancels the 2D jobs of the running chunks. The sieve restores its chunks from their folders and marker files (`doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`, section 4), including a chunk final ingestion staged; an unfinished chunk whose job recorded itself and wrote no exit status has that job cancelled and its folder set aside (`<folder>_unfinished<k>`), and is made afresh from its project as made. Every set the sieve chunked from (its `chunked_mics.txt`, one line per chunked micrograph) is imported again with those micrographs marked chunked, so the rest of a partly chunked set is still sieved, and goes into the watcher's history; sets imported but not chunked from are imported again. The sieve is made as soon as the upstream folder is attached, with nothing new to import. A leftover `TERM_STREAM` is removed. |
| p06 pool 2D | A stop cancels the running pool iteration. A restart cancels an iteration a crashed stage left running, then removes the previous pool's files from its folder (`cleanup_root_folder`: images, text, STAR, JPEG, binary and data files, the chunk folders, `TERM_STREAM`), and `REFINE2D_FINISHED`, the iteration's exit status and its project as made. The pool starts again from every set the sieve has handed off. Snapshots are kept. Publications for 3D continue their numbering (`restore_export_id`). |
| p07 multistate 3D | A stop cancels the running job. A leftover `TERM_STREAM` is removed and the stage starts again from the newest publication not listed in `rejected_publications.txt`, which is its first set: `solve2D` and `solve3D` run again in their folders. A run folder (`solve2D`, `solve3D`, `solve3D_addon/it_<n>`) whose job recorded itself and wrote no exit status is first moved aside to `<folder>_unfinished<k>` (`fresh_job_dir`), so a run never shares a folder with a job that may still be running. |

## 4. What carries state across a restart

| State | Kept in | Read by |
|---|---|---|
| completed job sets, their numbering and origins | the completed folder (`DIR_STREAM_COMPLETED`) of p01 and p04 | `stream_job_sets%restore` |
| optics-map ids | the maps in p02's folder | p02 (next id); readers take the newest |
| published picking references | `INITIAL_ANALYSIS_PICKREFS` in p03's folder | p03, p04 |
| sieve chunks and their state | chunk folders and their marker files (`SOLVE2D_FINISHED`, `REJECTION_FINISHED`, `COMPLETE`, `REJECTION_FAILED`, `FINAL_INGESTION`), a chunk's job record and exit status (`EXIT_CODE_solve2D`, `.job`), `chunked_mics.txt` | the sieve in p03 and p05 |
| queued jobs | each job's record (`<exit-status file>.job`: pid, host, scheduler and id) and exit status, beside the job | the stage's cancel on stop; a restart's fresh-folder check; the liveness check |
| the sieve's hand-offs, with its empty final set (`sieve_final_c<n>_f<m>.simple`) | p05's completed folder | p06 |
| upstream idle or stopped | `STREAM_IDLE`, `STREAM_FINISHED` in p01's and p04's folders; both removed when the stage starts, the idle marker when new work arrives | p03 and p04 (p01's), p05 (p04's) |
| publications for 3D and their numbering | p06's completed folder (the newest two) | p06 (next id), p07 |
| snapshots | `snapshots/` in p06's folder | the GUI, users |

Anything kept only in memory is lost:
- the pool's history of iterations;
- the 2D pool itself, which starts again;
- the sieve's final-ingestion flag, which its trigger sets again from the markers;
- the GUI metadata the master held, which the stages send again.

## 5. Invariants a change must keep

1. **Hand-offs are atomic.** Every file a stage hands on is written under a temporary name and
   renamed, so a restart never takes in a partial file.
2. **Numbering never goes back:** job sets, optics-map ids and publication ids continue after the
   highest one on disk.
3. **Completed work is not repeated.** A stage that restores its completed inputs (p01, p04, p05)
   puts them in its watcher's history.
4. **Downstream copes with repeats.** A stage that redoes its work by design (p02, p03 before
   publishing, p06, p07) produces outputs that downstream copes with: p07 matches a
   publication's particles by stack and image, so a restarted pool's order changes nothing, and
   p02's new maps are taken by id.
5. **Published picking references are final.**
6. **Each stage's restart is documented** in its module header (RESTART) and in section 3.

## 6. Change rules

- A change of what a stage restores updates its module header and section 3.
- A new marker file, or file that carries state, goes in section 4, and is kept by any cleanup
  that runs on restart.
- Tests: each stage tester in `unit_stream` covers its restart:
  - preprocessing: the previous run's sets imported again;
  - optics assignment: the map ids and the `TERM_STREAM` removal;
  - initial analysis: published references are final;
  - reference picking: the restart history and clearing;
  - pool 2D: the publication numbering, and the removal of the previous iteration's completion,
    exit status and project;
  - particle sieving: the `TERM_STREAM` removal, the restore of a partly chunked set from
    `chunked_mics.txt`, and the sieve made with nothing new to import;
  - multistate 3D: the `TERM_STREAM` removal;
  - the master (`unit_stream` "stream master"): a skipped stage is never started;
  - the sieve (`unit_project` "particle sieve"): its import and the final-ingestion chunk;
  - jobs (`unit_parallel` "qsys control"): the job record, the fresh folder, the multi-job
    script's exit status.
  The chained tests (`lib_stream` "sieve to 3D", "movies to 3D") run the stages from fixtures
  without the master; they do not restart a stage.

## Job lifecycle

Every queued job of the stream (follow-up plan, decisions 28-30):
- runs through `qsys_async_job` (p03's and p07's jobs, p04's `make_pickrefs`, the sieve's chunks,
  the pool iteration) or, for p01's and p04's job sets, the queue controller's streaming
  scheduler; both write the job's record and exit status beside it;
- is checked for liveness while it has no exit status: every `JOB_LIVENESS_S` (5 minutes) its
  scheduler is asked whether its id still exists (`squeue`, `bjobs`, `qstat`), or, for a local job,
  whether its process does on the recorded host (`query_job`). A job gone at two checks in a row
  (one miss can be a scheduler's completing state) gets `JOB_LOST_EXIT_CODE` (254) written as its
  exit status, so the caller's failure path runs (a retry, or the stage's stop);
- is cancelled through its record: `scancel`, `bkill` or `qdel` for a scheduler job; a local job
  runs in a session and process group of its own (`setsid`, where the system has it), and its
  cancel signals the group, reaching the part jobs a distributed program starts; a
  persistent-worker task, the script and its children.

## 7. Known gaps

- **Job sets keep the controller's streaming scheduler.** p01's and p04's sets run through
  `qsys_ctrl`'s streaming scheduler, not `qsys_async_job`, since the stream programs they run
  rely on its per-part numbering and done markers; they share the record, exit status, liveness
  check and cancel (follow-up plan, decision 28, done for every other job).
- **A liveness check needs the scheduler's answer.** A scheduler that cannot be reached, or a
  local job on another host than the stage's, cannot be checked; such a job is waited for.
- **The master forks a restarted stage from a multithreaded process.**
  - The memory monitor is fork-safe, and forks are made under the listener's lock.
  - Other threads (the persistent-worker server's) are still copied in whatever state they are.
  - Proposal: start a stage by fork and `exec` of `simple_stream prg=<stage>`, with the pipe
    descriptors on its command line.
- **p06 starts the pool from scratch.** A restart repeats every iteration. Until the restarted
  pool has classified a stack again, p07 holds its rows deselected (the first set's keep their
  selection). With publications on disk, the restarted pool publishes from iteration 25 only.
- **p07 runs `solve2D` and `solve3D` again** on everything, the newest publication being its
  first set.
- **p03 repeats its plan** when it stopped before publishing references.
