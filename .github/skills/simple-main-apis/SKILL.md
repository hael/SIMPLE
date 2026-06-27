---
name: simple-main-apis
description: Use when working in SIMPLE's src/main/apis subsystem, including public API aggregator modules that bundle imports for executables, commanders, core helpers, GUI metadata, orientation APIs, polar Fourier APIs, and stream/test entrypoints.
---

# SIMPLE `src/main/apis`

`apis/` is the import-aggregation layer.

## Read First

- `simple_exec_api.f90`
- `simple_test_exec_api.f90`
- `simple_stream_api.f90`
- `single_exec_api.f90`
- `simple_commanders_api.f90`
- `simple_core_api.f90`

## Role

- Collects and re-exports the modules an executable family needs
- Keeps production entrypoints thin
- Often reveals the real dependency surface of a workflow faster than reading the executable alone

## Working Rule

If you add a new exec/commander dependency and an entrypoint should see it, update the relevant API aggregator instead of bloating the executable file.
