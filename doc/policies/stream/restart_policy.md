# Stream Restart Policy

What happens when a stream stage is started again, by the GUI or in a run folder that already
holds its output. What each stage keeps, restores and redoes, and what a change must preserve.

## 1. How a stage is started again

1. **From the GUI**, after the stage has stopped (`restart_<key>` in a heartbeat answer). The
   master:
   - drops what the stopped stage left in its pipes, both its unread messages and the updates it
     never read (`discard_pipes`);
   - forks the stage again with the same command line, while holding the listener's lock.
2. **By starting the master again** in a run folder where the stages' folders exist. The master
   forks every stage as at a first start, and each stage finds its folder.
3. **By starting one stage on its own**, in its folder or with `dir_exec=<folder>` (p04, p06).

The master never restarts a stage by itself. Stages are forked without `forked_process`
auto-restart, so a stage that fails stays failed until the GUI restarts it. The master holds no
restart state of its own: each stage decides what it restores.

## 2. How a stage recognises a restart

A stage recognises a restart by its output folder existing before `params%new` creates it, or by
`dir_exec` where the stage takes one. It reads the command line for this, before the parameters
exist, and then removes `dir_exec` from the command lines it passes on.

## 3. Per stage

| Stage | On restart |
|---|---|
| p01 preprocessing | Its job sets are restored: the accepted micrographs of the completed sets are imported again, their movies go into the watcher's history, and set numbering continues after the highest completed set. Unfinished sets are dropped and their movies submitted again. The micrograph STAR file is rewritten. A leftover `TERM_STREAM` is removed. |
| p02 optics assignment | The micrograph segment starts empty, and every completed upstream project is imported again (the watcher's history starts empty). The groups are assigned again; optics-map ids continue after the newest map in its folder. A leftover `TERM_STREAM` is removed. |
| p03 initial analysis | When the picking references are published, they are sent to the GUI again and the stage is finished at once (`doc/policies/stream/reference_generation_policy.md`). Otherwise the plan starts again from cycle 1, and every completed upstream project is imported again. The sieve takes up the chunks in its folders, and the highest-numbered `solve3D_cavgs` folder holds a 3D result. p03 uses no `TERM_STREAM`. |
| p04 reference picking | Its job sets are restored: the completed sets are imported again, and the upstream projects they were made from go into the watcher's history, so none is picked twice. Set numbering continues; unfinished sets are dropped and their inputs submitted again. A leftover `TERM_STREAM` is removed. The picking templates are made again from the published references when the first new set arrives. |
| p05 particle sieving | The sieve restores its chunks from their folders and marker files (`doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`, section 4), including a chunk final ingestion staged. The projects it has already chunked (its `imported_projects.txt`) go into the watcher's history; projects imported but not yet chunked are imported again. A leftover `TERM_STREAM` is removed. |
| p06 pool 2D | The previous pool's files are removed from its folder (`cleanup_root_folder`: images, text, STAR, JPEG, binary and data files, the chunk folders, `TERM_STREAM`). The pool starts again from every set the sieve has handed off. Snapshots are kept. Publications for 3D continue their numbering (`restore_export_id`). |
| p07 multistate 3D | A leftover `TERM_STREAM` is removed and the stage starts again from the newest publication; `solve3D` runs again in its folder. |

## 4. What carries state across a restart

| State | Kept in | Read by |
|---|---|---|
| completed job sets, their numbering and origins | the completed folder (`DIR_STREAM_COMPLETED`) of p01 and p04 | `stream_job_sets%restore` |
| optics-map ids | the maps in p02's folder | p02 (next id); readers take the newest |
| published picking references | `OPENING2D_PICKREFS` in p03's folder | p03, p04 |
| sieve chunks and their state | chunk folders and their marker files (`SOLVE2D_FINISHED`, `REJECTION_FINISHED`, `COMPLETE`, `REJECTION_FAILED`, `FINAL_INGESTION`), `imported_projects.txt` | the sieve in p03 and p05 |
| the sieve's hand-offs | p05's completed folder | p06 |
| publications for 3D and their numbering | p06's completed folder (the newest two) | p06 (next id), p07 |
| snapshots | `snapshots/` in p06's folder | the GUI, users |

Anything kept only in memory is lost:
- the pool's history of iterations;
- the 2D pool itself, which starts again;
- the sieve's final-ingestion flag, which its trigger sets again;
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
  - pool 2D: the publication numbering;
  - particle sieving and multistate 3D: the `TERM_STREAM` removal;
  - the sieve (`unit_project` "particle sieve"): its import and the final-ingestion chunk.

## 7. Known gaps

- **The master forks a restarted stage from a multithreaded process.**
  - The memory monitor is fork-safe, and forks are made under the listener's lock.
  - Other threads (the persistent-worker server's) are still copied in whatever state they are.
  - Proposal: start a stage by fork and `exec` of `simple_stream prg=<stage>`, with the pipe
    descriptors on its command line.
- **p06 starts the pool from scratch.** A restart repeats every iteration. Until the restarted
  pool has classified a stack again, p07 holds its rows deselected.
- **p07 runs `solve3D` again** on everything.
- **p03 repeats its plan** when it stopped before publishing references.
