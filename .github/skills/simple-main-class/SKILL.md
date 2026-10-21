---
name: simple-main-class
description: Use when working in SIMPLE's src/main/class subsystem, including the class averager (online class sums, restoration, even/odd/merged class averages, part reduction, carried fractional-update sums), class FRCs, the on-disk class sums (cavg_sums), class averages matched against state volumes, and size-based class compatibility, used by refine2D, make_cavgs, the stream and solve3D state assignment.
---

# SIMPLE `src/main/class`

This folder owns class-average products and related statistics.

## Read First

- `simple_classaverager.f90`: the `cavger_*` interface and module state (`cavgs_even`, `cavgs_odd`, `cavgs_merged`)
- `simple_classaverager_restore.f90`: online sums (`cavger_init_online`, `cavger_update_sums`), restoration, reduction of parts, carried sums
- `simple_classaverager_core.f90`: the `stack` and `cavgs_set` work types
- `simple_cavg_sums.f90`: the carried set (`cavg_state.bin`) and the workers' contributions
- `simple_class_frcs.f90`: class FRC bookkeeping
- `simple_cavgs_state_match.f90`: class averages matched against state volumes (`assign_cavgs_states`)
- `simple_class_compatibility.f90`: support model and size-based class rejection (particle sieve, stream p03)
- `doc/policies/2D/solve2D_policy.md` (ownership, carry-over), `class_average_bootstrap_policy.md`, `particle_cache_policy.md`

## Ownership

- Assembly from partial sums is `cavgassemble_exec` in `src/main/strategies/parallelization/simple_cavgassemble_strategy.f90`, run on the caller's `(params, build)` by the distributed refine2D master; `run_cavgassemble` builds its own for the `cavgassemble` commander (`simple_commanders_mkcavgs.f90`, a wrapper)
- The 2D matcher (`simple_strategy2D_matcher.f90`) accumulates the online sums; the refine2D strategy restores them in shared memory (`restore_cavgs_inmem`)
- Only the assembly owner blends and publishes the carried set; workers accumulate the current sample from zero

## Guardrails

- `cavger_update_sums` buckets a batch by (class, half) and runs the buckets under a dynamic schedule, largest first, each in batch order, so the sums do not depend on the thread count; keep it so
- `cavger_init_online(cropped_ptcls=.true.)` means the particles handed in are at `box_crop`; the 2D matcher does this when cached or when `crop_particles2D` holds (`box_crop < box`, `ptcl2D`, outside the stream pool), cropping in memory from the batch's single read
- The Fourier crop round trip (crop, inverse FFT, forward FFT) is exact only to rounding: cropped and full-box sums differ, never expect bitwise identity between them

## Connections

- Depends on `image/`, `interp/` (KB kernel), `ctf/`, `project/`, and `ori/`
- Tests: `unit_reconstruction` ("class-average carry-over"), `unit_pftc_align2D3D` ("cavgs state match"), `unit_project` ("class compatibility")

## Working Rule

When a bug looks like “class products are wrong,” check both the averaging core and the project/orientation state feeding it.
