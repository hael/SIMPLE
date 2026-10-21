---
name: simple-main-stream
description: Use when working on SIMPLE's streaming pipeline, including the p00-p07 commanders in src/main/commanders/stream, the stage types in src/main/stream/stages, the master's parts in src/main/stream/master, the shared stage modules in src/main/stream/shared (pipes, sigterm, job sets, GUI senders, micrograph and optics helpers), the 2D pool in src/main/stream/pool2D, and the src/main/stream support layers (watchers, particle-sieve integration, streaming variants of preprocessing and refine2D workflows).
---

# SIMPLE stream pipeline

The commanders live in `src/main/commanders/stream`, the stage types they drive in `src/main/stream/stages`, the master's parts in `src/main/stream/master`, the modules the stages share in `src/main/stream/shared`, and the 2D pool in `src/main/stream/pool2D`; `src/main/stream` itself keeps the watcher and the stream utilities.

## Read First

- `doc/refactoring_notes/planned/stream_refactor.md`: what each module does, how each stage differs from the stage it replaced, and the remaining clean-up
- the policy for the contract you change, in `doc/policies/stream/` (index and the per-stage `pNN_overview.md` diagrams in its README):
  - `reference_generation_policy.md`: p03's routes to the picking references, their publication and precedence, p04's use of them
  - `stream_3D_ingestion_policy.md`: p06's publications for 3D, p07's import of them, its reference volumes, the assignment of each publication to them, and its 3D snapshots
  - `pool2D_policy.md`: the 2D pool's schedule, sampling, resolution, dimensions, snapshots and final project
  - `ipc_policy.md`: the stage-to-master pipes, the GUI metadata wire format, the master's heartbeat and the GUI's answers
  - `restart_policy.md`: what each stage keeps, restores and redoes when started again
- `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` for the particle sieve (p03, p05)

Pipeline stages (one forked process each, launched by p00); commanders in `src/main/commanders/stream`, stage types in `src/main/stream/stages`:

- `simple_commanders_stream_p00_master.f90`: GUI-side master; the stage table (`simple_stream_master_stage`, ids in `simple_stream_master_stage_ids`), a listener thread filling `simple_stream_master_meta_store`, the heartbeat loop that applies `simple_stream_master_gui_commands`, and the stages' threads and parts (`simple_stream_master_resources`: named defaults, environment overrides)
- `simple_commanders_stream_p01_preprocess.f90` to `simple_commanders_stream_p07_solve3D_multistate.f90`: thin commanders, each driving one stage type:
  - p01 `simple_stream_stage_preprocess`: motion correction, CTF estimation, segmentation picking
  - p02 `simple_stream_stage_optics`: optics-group assignment and maps
  - p03 `simple_stream_stage_initial_analysis`: opening solve2D / solve3D_cavgs on segmentation picks; produces picking references
  - p04 `simple_stream_stage_refpick`: reference-based picking and extraction
  - p05 `simple_stream_stage_sieve`: continuous particle sieving via `ptcl_sieve`
  - p06 `simple_stream_stage_pool2D`: global 2D pool classification by the pool object it holds (`stream_pool2D`), snapshots, publication of the classified pool state for 3D
  - p07 `simple_stream_stage_solve3D`: multistate 3D by class-average assignment (`doc/policies/stream/stream_3D_ingestion_policy.md`, sections 4, 9 and 10); its first set is p06's first publication, the sieve's own 2D of the sets handed off, taken whole and keeping its selection for the session; from it p07 makes reference volumes once (solve2D into 150 classes, the pool model's selection, cls_expansion into about 500, `solve3D_cavgs`), checked (even/odd angular distance and state agreement, two populated states), retried, with a fallback on a later publication's pool classes, and kept in `state_refs/`; each later publication is merged into rows that only grow, scored with the pool model (logged, not applied), split by cls_expansion, matched against the references by `assign_cavgs_states` (simple_private_exec, `simple_cavgs_state_match`), rejected by fit, mapped onto its particles, whose poses are then refined locally with the states fixed (`refine3D refine=cont`), in `assign/it_<n>/`; it reads GUI updates and writes 3D snapshots of selected states (merged into state 1, with their volumes); its rules, check, assignment mapping, refinement seeds and job command lines live in `src/main/solve/simple_state_assignment.f90`, shared with the batch program `solve3D_pyramid`, so a change there changes both

Shared pieces in `src/main/stream/shared`: `simple_stream_pipe` (framing), `simple_stream_state` (pipe descriptors), `simple_stream_sigterm`, `simple_stream_gui_senders`, `simple_stream_job_sets`, `simple_optics_maps`, `simple_optics_groups`, `simple_mic_import`, `simple_mic_selection`, `simple_stream_meta_plots`.

The 2D pool in `src/main/stream/pool2D`:

- `simple_stream_pool2D.f90`: the pool as a type (`stream_pool2D`, private components; tester `simple_stream_pool2D_tester`): iterations, history, dimensions and mask, snapshots, publications for 3D, the final project; each iteration is a queued `refine2D_distr` job (simple_private_exec, the batch `commander_refine2D` and its strategy, workers `prg=refine2D`), so changes to batch refine2D reach the pool
- `simple_stream_refine2D_utils.f90`: stateless helpers: folder clean-up, iteration files, set appending, class draws, publication building and naming, snapshot sprite sheets

Supporting layers in `src/main/stream`:

- `simple_stream_watcher.f90`: directory watcher used both for movie import and for stage-to-stage project handoff
- `simple_stream_utils.f90`, `simple_mini_stream_utils.f90`: utilities

## Structure

- Stages exchange data through project files in watched directories (written as `.tmp` and renamed into place), and GUI metadata through `stream_pipe` frames on pipes (stage->master `*_in`, master->stage `*_out`, named in `simple_stream_state`)
- A stage is a type with public components and named steps (`new`, `iterate`, `finished`, `finalize`, `kill`); its tester assembles it step by step without a queue or waits (`unit_stream`)
- p06 holds the pool as an object, empty from `init_params` and started at the first import; the pool's state is private: p06 hands it imported sets (`append_sets`) and reads what the GUI shows from `stats()`; its files have fixed names, so a folder holds one pool

## Particle Sieving

p03 and p05 use `ptcl_sieve` (`src/main/sieve/simple_ptcl_sieve.f90`). Read `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md` before changing sieve lifecycle, chunk tiers, hand-off, or rejection, and keep `new`/`kill` symmetry in the stages that drive it.

## Guardrails

- p00's listener thread holds the metadata lock while draining the stages' readers into the store; any main-thread read or `discard` of a stage pipe must take the same lock
- A forked child copies the master's memory but not its threads: fork a stage only before the listener starts or while holding its lock, and keep C state that a thread owns fork-safe (`pthread_atfork`, as the memory monitor does)
- `forked_process` never restarts a child by itself: a failed stage stays failed and is forked again (`stream_master_stage%start`) only when the GUI asks
- The master returns from `execute`, so the entry point stops the persistent worker server it owns; forked stages exit directly from `forked_process` and never run the entry point's tail
- p03 publishes the picking references once, as `INITIAL_ANALYSIS_PICKREFS`, by a rename (`publish_pickrefs`). A GUI selection pre-empts the 3D route, and published references are final; don't write that file from anywhere else
- Watcher handoff assumes files are only returned once untouched for `report_time`; do not weaken that check, and write files a stage hands on as a temporary and rename them
- GUI metadata is serialised by its bytes (`transfer`): never send an object with an allocated component (see `gui_metadata_vol3D%serialise`)
- The GUI metadata types and the assembler import nothing that reaches `src/main` (IPC policy §8): domain work for the GUI belongs to the caller (`simple_gui_project_builder`, `simple_oris_utils`, `simple_stream_meta_plots`)
- Signal handlers only set a flag (`simple_stream_sigterm`); poll it between steps, never log, join or exit inside a handler
- The pool command line sets `msk_crop` explicitly, so changing `mskdiam` must also recompute `msk_crop` (`set_mask`, as `set_mskdiam` does)
- The pool may run downscaled (its private working `dims` vs its native box/smpd; its command line carries `box_crop`/`smpd_crop`); anything written for downstream use must be rescaled or labelled with the matching smpd
- The `projfile` of a project's `projinfo` may be a bare basename; pass an explicit filename to `write_segment_inside`/`write` from stage code
- Particle fields read by stream selection must have a live writer; check before relying on one

## Working Rule

Do not assume batch semantics apply unchanged to stream code; stage boundaries and artifact cadence matter here.
