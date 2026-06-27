---
name: simple-main-root
description: Use when working with the root modules directly under SIMPLE's src/main directory, especially orchestration modules like builder, parameters, command line parsing, convergence, ab initio controllers, probability tables, simulation, symmetry, and other cross-subsystem foundation code.
---

# SIMPLE Main Root

This skill covers files directly in `src/main/`, not the subfolders.

## Key Files

- `simple_builder.f90`
- `simple_cmdline.f90`
- `simple_parameters.f90`
- `simple_convergence.f90`
- `simple_eul_prob_tab.f90`
- `simple_eul_prob_tab2D.f90`
- `simple_eul_prob_tab_neigh.f90`
- `simple_abinitio_controller.f90`
- `simple_abinitio2D_controller.f90`
- `simple_abinitio_utils.f90`

## What Lives Here

- Cross-cutting orchestration state
- CLI parsing and validation
- Global parameter materialization
- Probability-table infrastructure for search workflows
- Builder/toolbox construction
- Shared controllers for major workflows

## When Editing

- Check whether the behavior really belongs at the root or should live in a subsystem folder.
- Be careful with changes here because these modules often affect many executables and workflows.
