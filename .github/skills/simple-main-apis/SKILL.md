---
name: simple-main-apis
description: Use when working in SIMPLE's src/main/apis subsystem, including the per-executable API modules that wire exec routers into simple_exec, single_exec, simple_test_exec and simple_private_exec, and the existing umbrella modules (simple_core_module_api, simple_commanders_api, simple_stream_api, simple_ori_api, simple_pftc_api, simple_pftc_srch_api, simple_gui_metadata_api) that new code must not widen.
---

# `src/main/apis`: import bundles

Two kinds of module, both re-exports. The per-executable API modules are the one import of their
front door: `simple_exec_api` (`production/simple_exec.f90`), `single_exec_api`,
`simple_test_exec_api` and `simple_private_exec_api` (`production/simple_private_exec_driver.f90`);
`production/simple_stream.f90` imports its commanders directly. The umbrellas re-export many modules
to many importers and are not to grow. What a router or commander does is in `exec/` and
`commanders/`.

## Rules

- Do not add an umbrella module or a re-export, and do not add a name to an existing umbrella
  (`doc/policies/compile_time_policy.md`): each module-level `use` is a build-order edge for every
  importer. The umbrellas and their importers: `simple_core_module_api` (file
  `simple_core_api.f90`, about 280 modules), `simple_commanders_api` (about 60: commanders,
  strategies, `solve/`, root), `simple_pftc_srch_api` (about 30, mostly `strategies/search`),
  `simple_ori_api` (the `ori/` modules), `simple_stream_api` (`sieve/simple_ptcl_sieve_utils.f90`),
  `simple_pftc_api` (`pftc/simple_polarft_calc.f90`), `simple_gui_metadata_api` (three testers),
  and `simple_ui_modules` in `ui/`. New code imports what it uses with `use <module>, only:`.
- An umbrella that imports a module whole grows with it: `simple_core_module_api` re-exports all of
  `simple_defs`, `simple_defs_fname`, `simple_defs_conv`, `simple_defs_stream`, `simple_type_defs`,
  `omp_lib` and `iso_c_binding`, so a public name added to one of them reaches every importer, where
  it can clash with a local name. A module an umbrella imports, even with `only:`, compiles before
  all of the umbrella's importers: a module-level `use` added to `simple_exec_helpers`,
  `simple_builder` or `simple_parameters` delays every importer of `simple_commanders_api`.
- `simple_pftc_srch_api` re-exports `polarft_calc` with what the search strategies use beside it
  (its header: to avoid a circular dependency involving `polarft_calc`). `simple_pftc_api` and
  `simple_ori_api` are the import bundles of `polarft_calc` and of `ori/`: what is added there
  compiles before them.
- A new router module's `use` goes into the API module of its executable, not into the
  `production/` file, which holds only the call (skill `simple-main-exec`). A new worker's
  commander goes into `simple_private_exec_api`; the dispatch is `dispatch_private_prg` in
  `production/simple_private_exec_driver.f90` (`simple_private_exec_api`'s header lists the steps).
- Removing a name from an umbrella is a compile error in each module that used it through the
  umbrella: give that module its own `use ..., only:`.

## Tests

The build is the test. Claim a compile-time gain only from a like-for-like
`scripts/profile_build.sh` comparison (CLAUDE.md).

## See also

`doc/policies/compile_time_policy.md`; skills `simple-modern-fortran`, `simple-main-exec`.
