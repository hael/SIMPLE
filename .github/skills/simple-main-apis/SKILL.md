---
name: simple-main-apis
description: Use when working in SIMPLE's src/main/apis subsystem, including the per-executable API modules that wire exec routers into simple_exec, single_exec, simple_test_exec and simple_private_exec, and the existing umbrella modules (simple_core_module_api, simple_commanders_api, simple_stream_api, orientation, polar Fourier and GUI metadata APIs) that new code should not widen.
---

# SIMPLE `src/main/apis`

`apis/` is the import-aggregation layer.

## Read First

- `simple_exec_api.f90`
- `simple_test_exec_api.f90`
- `simple_stream_api.f90`
- `single_exec_api.f90`
- `simple_private_exec_api.f90`
- `simple_commanders_api.f90`
- `simple_core_api.f90` (defines module `simple_core_module_api`)

## Role

- Collects the modules an executable family needs: each `*_exec_api` is
  imported by its `production/` entrypoint (`simple_private_exec_api` by
  `simple_private_exec_driver.f90`). `production/simple_stream.f90` imports the
  stream commanders directly; `simple_stream_api` is a stream-helper umbrella
- Keeps production entrypoints thin
- Often reveals the real dependency surface of a workflow faster than reading the executable alone
- `simple_pftc_srch_api` is a facade that breaks a circular dependency around `polarft_calc`

## Working Rule

- A new exec router goes into the API module of the executable that calls it,
  not into the executable file.
- Do not add new umbrella modules or re-exports, and do not widen the existing
  umbrellas (`simple_core_module_api`, `simple_commanders_api`,
  `simple_stream_api`): each
  module-level import is a build-order edge for every importer. New library code
  imports the modules it uses with `use ..., only:`
  (`doc/policies/compile_time_policy.md`).
