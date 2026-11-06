---
name: simple-main-strategies
description: Use when working in SIMPLE's src/main/strategies subsystem, the execution-policy layer between commanders and domain objects: shared-memory and distributed-master strategies, the master's steps on its own project (volume and class-average assembly, the sigma2 commit, the sample draw, the merge of the parts' readout records), and the search engine (refine3D/refine2D matchers and their probabilistic batch pipeline, per-particle 2D/3D strategies, particle batches, sampling, probability tables). Routes to simple-main-strategies-parallelization, -search and -probabilistic.
---

# `src/main/strategies`: how a workflow runs

Owns the execution policy between commanders and domain objects: shared memory or the distributed
master, the iteration, the master's steps on its own `(params, build)` (`parallelization/`), and the
search engine: matchers, per-particle 2D/3D strategies, particle batches, sampling (`search/`) and the
probability tables (`search/probabilistic/`). It does not own what the data means (`pftc/`, `class/`,
`sigma2/`, `volume/`, `project/`, `ori/`), the command line (`ui/`, `params/`) or the sequencing of whole
programs (`solve/`, `refine3D_states`). The rule set is `doc/policies/layering_and_project_ownership_policy.md`.

## Rules

- While a pass runs, `omp_get_max_threads()` equals `nthr_glob`. Per-thread state is indexed by
  `omp_get_thread_num() + 1` and sized by the team in use (`pass_team()` or `omp_get_max_threads()` read
  before the allocation): the search singletons `s3D`, `s2D`, `alloc_ptcl_imgs`, the reconstruction's image
  heaps, the polar calculator's work buffers (at `pftc%new`), the class averager's work images and the
  helpers (`corrmat`, the 2D utilities, `pspecs`, `calc_pspec`). The probability tables (`probabilistic/`)
  still use `nthr_glob`: a larger team would write out of bounds there, so each pass stops on one
  (`check_pass_team` at the start of `refine3D_exec` and `refine2D_exec`). A
  master step sets both
  (`omp_set_num_threads(nthr)`, `nthr_glob = nthr`); its caller saves both and restores both.
- One pass at a time per process: `s3D`, `s2D`, the 2D matcher's module pointers (`p_ptr`, `b_ptr` of
  `simple_strategy2D_matcher`), the class averager and the particle buffer are module state. The probability
  table is a local of the pass (`refine3D_exec`, `refine2D_exec`); nothing is handed over between passes.
- A part of a stage that the master or another workflow runs is a routine `<x>_exec(params, build, ...)` in
  `parallelization/` with a `run_<x>(cline)` entry; its commander is a thin `execute(cline)` calling
  `run_<x>`. A strategy never calls a commander for a part of a stage (shared-memory refine2D/refine3D
  calling the `commit_sigma2` commander is the one accepted exception) and never imports the commander
  that wraps it (an in-process commander comes in as a `commander_base` dummy). `rg` the verb before writing a step.
- Workers read the project file and write part files; only the owner of the project writes it.

## Subfolder guides

- `parallelization/`: `simple-main-strategies-parallelization` (strategies, the master's iteration and writes).
- `search/`: `simple-main-strategies-search` (matcher pass, batches, per-particle searches).
- `search/probabilistic/`: `simple-main-strategies-probabilistic` (tables, the batch pipeline's fill and per-row
  draw, the keyed draws, readout records).

## Tests

Unit suites are in the subfolder guides. Workflow checks (`build/bin/simple_test_exec test=<entry>`):

- A change that claims identical results: `test=identity` on two builds at one thread
  (`doc/policies/test_environment_policy.md` section 1.5); never judged on multi-threaded runs.
- `simulated_workflow suite=1jxy`: its second `solve3D` runs `nparts=2`, the only registered run of the
  probabilistic refine3D master over several parts (the sample draw, the parts' batch pipeline and the readout
  merge, assembly, sigma2 commit). The map correlation is a robustness profile (1-2 misses in 20 runs is
  normal): one run shows breakage, not improvement; compare `test=profile` runs with `dir1=`/`dir2=`.
- `cont_refine3D_1jxy` (`refine=cont` over two parts, the `pose_cont` polish), `multistate_workflows`
  (`suite=states` or `consensus`; 1 913 s at eight threads for all routes), `flex_pca_blobs` (two workers).
- The distributed refine2D master runs only in the stream's 2D pool (`mini_stream_*`, `stream_chains`;
  `stream2d=yes`, one part in `mini_stream`); check a non-stream 2D master change by hand with `nparts=2`.

## See also

`doc/policies/3D/refine3D_policy.md` sections 3-4, `doc/policies/2D/solve2D_policy.md`,
`doc/implementation_notes/completed/prob_align_batch_pipeline.md` (the batch pipeline that replaced the global table);
skills `simple-refine3d`, `simple-solve2d`, `simple-solve3d-importance-sampling`, `simple-modern-fortran`.
