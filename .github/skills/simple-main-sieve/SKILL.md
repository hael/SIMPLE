---
name: simple-main-sieve
description: Use when working in SIMPLE's src/main/sieve subsystem, the particle sieve (ptcl_sieve) used by stream initial analysis (p03), stream particle sieving (p05) and the particle_sieving program, including coarse and fine 2D chunks, their solve2D jobs, class-average rejection (sieve quality preset, class compatibility), chunk sentinel files and restart import, retries of failed chunks, final ingestion and the empty final set, hand-off to the completed folder with optics maps, chunked_mics.txt, counters and combine_completed_chunks.
---

# `src/main/sieve`: the particle sieve

Owns `ptcl_sieve`: it cuts imported particle records into coarse chunks (and merges passed coarse
chunks into fine ones), runs a `solve2D` job per chunk, rejects class averages, hands terminal
chunks to a completed folder and counts what passed. `simple_ptcl_sieve_utils` cuts a project into
pre-chunked sub-projects for the `particle_sieving` program. It does not watch or import (the
driving stage: p03, p05, `simple_commanders_sieve`), score class averages (`cavg_quality`,
`class/simple_class_compatibility`), run queues (`qsys_async_job`) or decide what the pool does with
a hand-off (p06). Contract: `doc/policies/sieving_and_rejection/ptcl_sieve_policy.md`; read it
before changing the lifecycle, tiers, hand-off or rejection, and update it in the same change.

## Rules

Silent ones first: nothing reports these.

- Everything is relative to the working directory: `new` adopts every chunk folder already in
  `chunks_coarse/` and `chunks_fine/` (restart import) and `chunked_mics.txt` is written there;
  `submit` changes directory and restores it and `CWD_GLOB`. Make and cycle it in the driving
  stage's folder; a fresh run needs a fresh folder.
- `cycle` order is policy: `collect_and_reject`, `generate_chunks_coarse`, `generate_chunks_fine`
  (not with `single_pass`), `submit` (fine before coarse), `hand_off_final_set`.
- Sentinel files are the state a restart reads (`SOLVE2D_FINISHED`, `REJECTION_FINISHED`,
  `COMPLETE`, `REJECTION_FAILED`, `FINAL_INGESTION`, the job record and exit status). A chunk
  exists once its project is renamed in; what it took is written before. After rejection,
  `cleanup_chunk` deletes every file of the chunk folder not on its keep list: a new file a restart
  needs goes on that list and in the policy (sections 4 and 8).
- Hand-off: a copy in `completedir` written as `<stem>.tmp` and renamed, before `COMPLETE`; optics
  groups go on the copy only (fine chunks merge coarse ones and offset their group ids).
- The end of intake is only the empty `sieve_final_c<n>_f<m>.simple` (`sieve_final=yes`), handed
  on once final ingestion is set and every chunk is complete or failed. No chunk is marked final.
  p05 sets and withdraws final ingestion on p04's markers; p03 sets it once and never withdraws it.
- Rejection reasons are `<tier>_reject:<reason>` without blanks (the orientation reader splits at
  blanks). Both tiers use the `sieve` preset and hard-gate context; another is a policy change that
  needs a validation run.
- A chunk whose job fails is retried once from `chunk_input.simple` in a fresh folder (the old one
  kept as `<folder>_attempt1`); a second failure drops it with its particles counted. The stage
  never stops for a chunk.
- Counters are cumulative across restarts (rebuilt by `new` from the chunks on disk).
- The records marked included stay a prefix of `project_list`: coarse chunking takes the contiguous
  range from the first record not included. `generate_chunks_coarse` and `write_chunked_mics` stop
  otherwise (`check_included_prefix`). A restart pushes the chunked micrographs first (p05
  `restore_imports`).
- `cleanup_chunk` keeps what a chunk took (`chunk_mics.txt`, `consumed_coarse.txt`) with its sentinels:
  p05's restart rebuilds `chunked_mics.txt` from the coarse chunks' lists (`rebuild_chunked_mics`), and
  gives up, with a warning, only on a chunk without its list (`pre_chunked`, or older than the fix).

## Order of calls

`new(params, settings, completedir[, pre_chunked, optics_dir])` (calls `kill` first; `settings`
holds everything the sieve reads, `params` only the queue environment; tier overrides 0 mean the
default) -> { `set_final_ingestion`/`unset_final_ingestion`; `cycle(project_list)`; `get_*` } until
`get_finished` -> optional `combine_completed_chunks` -> on stop `cancel` (the running chunks' jobs),
then `kill`; `kill` alone leaves them running.

## Couplings

- Each chunk is a `solve2D` job: a change to solve2D/refine2D reaches every chunk; the `DEFAULT_*`
  constants are also the p05 commander's and `particle_sieving`'s defaults.
- p06 reads the hand-offs in chunk-id order with the final sets last (`sort_sets`); p05's mask is
  p04's `moldiam.txt`, p03's the box default (`reference_generation_policy.md` 3.1).

## Tests

- `build/bin/simple_test_exec test=unit_project suite=particle_sieve` (part of a 2-3 s suite):
  lifecycle, restart import, counters, `get_finished`, hand-off with optics, retries, the lists kept by
  the cleanup and the rebuild after a rejection, the prefix stop.
- `test=unit_stream suite=particle_sieving` and `suite=initial_analysis`: the drivers, a partly
  chunked set restored from `chunked_mics.txt`, the restore after a rejected chunk and a crash.
- `test=stream_chains suite=sieve_to_3d` runs real chunk jobs (after `make -C build install`, tens
  of minutes). A change to selection moves p05's counts: compare the chain's sieve and pool numbers.

## See also

`doc/policies/stream/pool2D_policy.md` (section 2), `restart_policy.md` (p05);
skills `simple-main-stream`, `simple-cluster-cavgs-quality`, `simple-solve2d`.
