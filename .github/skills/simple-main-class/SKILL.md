---
name: simple-main-class
description: Use when working in SIMPLE's src/main/class subsystem, including the class averager (online class sums, restoration, even/odd/merged class averages, part reduction, carried fractional-update sums), class FRCs, the on-disk class sums (cavg_sums), class averages matched against state volumes, and size-based class compatibility, used by refine2D, make_cavgs, the stream and solve3D state assignment.
---

# `src/main/class`: class averages and their sums

Owns the class averager (`simple_classaverager.f90`, a module singleton: `cavgs_even/odd/merged`, the
restoration state in `_restore`), its sums on disk (`simple_cavg_sums.f90`: the carried set
`cavg_state.bin`, a worker's contribution `cavg_contrib_fname(part)`), class FRCs (`simple_class_frcs.f90`),
class averages matched to state volumes (`simple_cavgs_state_match.f90`, for `assign_cavgs_states`) and the
size-based compatibility of the particle sieve and the stream's initial analysis
(`simple_class_compatibility.f90`). It does not decide whether an iteration blends, which particles are
summed or who assembles: `refine2D_blends_carryover`, `carryover_needs_full_update`, the 2D matcher
(`strategies/`); assembly is `cavgassemble_exec` (distributed) or `restore_cavgs_inmem` (shared memory).

## Rules

Silent ones first.

- Only the assembly owner reads or writes `cavg_state.bin`: `cavger_commit_carryover` (shared memory) or
  `cavger_assemble_sums_from_parts` (the distributed master). The matcher and the workers start from zero
  (`cavger_init_online`) and write only their current sample (`cavger_write_contribution`). Distributed
  `make_cavgs` also publishes `cavg_state.bin` (through `cavgassemble`), shared memory does not: it is not
  refine2D's starting set; read nothing from it after `make_cavgs`.
- The owner reads the particle field (state, `updatecnt`, class, eo) to weight the blend and recount the
  populations: the workers' documents must be merged into it first.
- The carried sums are shifted once, by each class's centring offset, before the blend. The offset is set
  where the reference is prepared (`cavger_set_center_offset` in `strategies/search/simple_matcher_pftc_prep.f90`)
  and travels in the contribution; workers that disagree get part 1's, with a warning.
- An unusable carried set (missing, other classes, box or sampling) restores from the current sample with
  only a warning; the strategies check first (`carryover_needs_full_update`) and run a full update.
- `cavger_update_sums` buckets a batch by (class, half) and runs the buckets largest first under a dynamic
  schedule, each in batch order, so the sums do not depend on the thread count: keep it so (no tester
  checks it). The master adds contributions in part order: sums differ with `nparts` by rounding.
- An active particle with `eo = -1` (in no half) stops the class sums and the carry-over recount
  (`oris%get_half`); `ensure_field_halves` does not report such a row (it can come from a reactivated
  `selection` of a merged chunk project).
  `cavger_transf_oridat(updated_only=.true.)` drops state-0 and never-updated particles without a message.
- `cavger_init_online(cropped_ptcls=.true.)`: the images handed in are at `box_crop`, noise-normalised at
  the full box and not again; shifts and CTF are scaled to the crop. The caller sizes them (`crop_particles2D`).
  Full-box images are normalised under the mask their alignment chose (`build%obs_on_box`). Cropped and
  full-box sums differ by rounding (crop, inverse and forward FFT): never expect bitwise identity.
- With `ml_reg` the sums use `build%esig`'s plane grid (`prep_fplane_grid`) and restoration the global
  spectra (`cavger_read_euclid_sigma2`).

## Where to start

- Online order: `cavger_new` -> `cavger_read_all` (references) -> `cavger_init_online` -> per batch
  `cavger_transf_oridat`, `cavger_update_sums` -> `cavger_dealloc_online` -> contribution or carry-over commit
  -> `cavger_restore_cavgs(frcs)` -> `cavger_gen2Dclassdoc` -> write -> `cavger_kill`. Offline (`make_cavgs`):
  `cavger_new`, `cavger_read_euclid_sigma2`, `cavger_assemble_sums`, restore.
- The averager reads particles through `ptcl_read_session` (`strategies/search/simple_matcher_ptcl_io.f90`),
  a domain-to-strategy import listed in the layering policy's section 8.

## Tests

- `build/bin/simple_test_exec test=unit_reconstruction suite=class_average_carry_over` (cavg_sums round trip,
  refusals, owner blend, padding; the suite about 9 s), `unit_pftc_align2D3D suite=cavgs_state_match` (2.5 s),
  `unit_project suite=class_compatibility`.
- Workflow: `simulated_workflow` (shared-memory solve2D), `mini_stream_*` and `stream_chains` (the pool's
  distributed refine2D master, carried sums, `cavger_pad_carried_sums`).

## See also

`doc/policies/2D/solve2D_policy.md` (ownership, carry-over), `doc/policies/2D/class_average_bootstrap_policy.md`,
`doc/policies/2D/particle_cache_policy.md`, `doc/policies/sieving_and_rejection/class_compatibility_policy.md`,
`doc/algorithms/refine2d_class_averaging.md`; skills `simple-solve2d`, `simple-cluster-cavgs-quality`.
