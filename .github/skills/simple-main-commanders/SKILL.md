---
name: simple-main-commanders
description: Use when working in SIMPLE's src/main/commanders subsystem, including the abstract commander_base interface (execute(cline) only), the SIMPLE, SINGLE and test commander families that own CLI program behavior, command-line normalization before params%new and child command lines, controllers that run whole programs as stages (solve3D, solve2D, refine3D_auto/local/states, solve3D_chunks/pyramid), thin wrappers of strategy routines, worker shells, the unit-suite tables and the high-level test registry of simple_test_exec.
---

# `src/main/commanders`: the command objects of the CLI programs

One `commander_<x>` per program, reached from the `exec/` routers and, for workers, from
`simple_private_exec`. A commander normalizes and validates `cline`, makes the one `params`, calls a
routine and ends; controllers (`exec_solve3D`, `exec_solve2D`, `exec_refine3D_workflow`,
`refine3D_states`, `solve3D_chunks`, `_pyramid`) sequence whole programs. Not here: a part of a stage
another workflow runs (`strategies/parallelization`: `<x>_exec(params, build, ...)` with
`run_<x>(cline)`), stage ladders (`solve/`), numerics (domain folders), `commanders/stream` (skill
`simple-main-stream`). Rules: `doc/policies/layering_and_project_ownership_policy.md`.

## Rules

Silent ones first: nothing reports these.

- A called commander builds its own `params` and builder from the file. Only a controller calls
  one, for a whole program (the shared-memory refine strategies' sigma2 commit through
  `commander_commit_sigma2` is the one ruled exception, policy section 2): write the segments it reads first (`write_segment_inside`),
  read back the ones it wrote after (`read_segment`); your in-memory project is stale across the
  call. `refine3D` and `refine2D` return their last iteration as `endit` on the caller's cline:
  delete it before a stage, number the next stage from it.
- `params%new` changes the process (skill `simple-main-params`): with `mkdir=yes` it enters a new
  `<n>_<prg>/` (copying a required project there), so set `mkdir=no` on every child cline after, or the
  child nests a directory and you read stale files. Each `new`, so each in-process commander call,
  reseeds (one more or less shifts every later draw of a seeded run; the keyed draws of `refine=prob*`
  excepted) and resets the threads.
- `cline%set` takes any key: `params%new` ignores one it does not bind, a worker's parse stops on
  one outside the vocabulary (`THROW_HARD`). Silent in shared memory, fatal distributed. No
  parameter for an internal switch; the in-process handshakes (`INPROCESS_KEYS`,
  `simple_cmdline.f90`) never reach job descriptions and, like the `internal:` keys, are ruled to be
  retired: add none. After `new`, read `params%`, not `cline`.
- One `params%new` per execution path, after the cline is final and validated: `exec_refine3D`
  refuses an unknown `refine=` before `create_refine3D_strategy`, so a master never waits on workers
  that stopped. A strategy-object program makes `params` in the strategy's `initialize`.
- A master step's commander is one call to `run_<x>` (`commander_volassemble`,
  `commander_cavgassemble`, `commander_commit_sigma2`). No typed entry beside `execute`: they were
  deleted and stay deleted. The probabilistic alignment has no master step and no commander: it runs
  inside the matchers (skill `simple-main-strategies-probabilistic`).
- A worker program is a shell: `params` and builder, the body, part finished
  (`commander_refine3D_distr_worker` hands `refine3D_exec` its `(params, build)`).
- The worker frame (`sample4update_reprod`, `split_nobjs_even`, `qsys_declare_part_finished`) is the
  matchers' (`refine3D_exec`, `refine2D_exec`): the three together in a new commander are a third copy
  (policy section 7). Deviations here (policy section 8), not patterns: `postprocess_volume_from_files`
  in `simple_commanders_volops.f90`, called by a strategy (`simple_rec3D_service`); controller routines
  of hundreds of lines.
- New code imports `builder`, `parameters`, `cmdline`, `sp_project` directly, not
  `simple_commanders_api`; shared defaults go in `src/defs/simple_default_clines.f90`.

## Where to start

`simple/simple_commander_base.f90`; a wrapper, `simple_commanders_rec_distr.f90`; a strategy object,
`exec_calc_pspec`; controllers, `exec_solve3D` and `exec_refine3D_workflow`. `single/`: SINGLE
programs (skill `simple-main-nano`).

## Tests

- `test/simple_commanders_test_class.f90` holds the `suites_<area>` tables of `test=unit_<area>`,
  `lib_<area>`, `units`. A high-level entry is a row of `test/simple_test_registry.f90` plus its
  `production/CMakeLists.txt` registration, added together (`scripts/check_test_registry.py`).
- A commander change: `scripts/run_fast_gate.sh build`, then its program's entry:
  `test=simulated_workflow suite=1jxy` (solve2D, solve3D, then solve3D with `nparts=2`; ~1000 s),
  `test=multistate_workflows suite=<cavgs|shc|states|consensus>`, `test=cont_refine3D_1jxy`
  (refine3D_auto/local), `test=solve3D_addon`, `test=flex_pca_blobs` (refine3D_states_cavgs).
- A simulated workflow's map correlation is a robustness profile (1-2 misses in 20 seeds is
  normal): one run shows breakage, not improvement (`test=profile`). A change that must not move
  results: one-thread runs on two builds, then `test=identity target=<entry> dir1= dir2=`.

## See also

`doc/policies/test_environment_policy.md`; skills `simple-main-root` (builder, cmdline, `solve/`),
`simple-main-strategies-parallelization`, `simple-main-params`, `simple-main-exec`.
