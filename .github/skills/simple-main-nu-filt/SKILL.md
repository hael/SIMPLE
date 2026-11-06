---
name: simple-main-nu-filt
description: Use when working in SIMPLE's src/main/nu_filt subsystem, including simple_nu_filter and its submodules, candidate-bank setup, the source pair of nu_filter_vols, like-for-like label selection, mask-packed objective costs, ordered-label Potts smoothing, tent smoothing, NU sub-timers, NU local-resolution output, the NU-evidence envelope and evidence state (nu_filt3D, postprocess_nu sharpening), and cleanup of module-level filter state.
---

# `src/main/nu_filt`: the volume-domain nonuniform filter

Owns `simple_nu_filter` and its submodules: the candidate bank (a static low-pass ladder capped from
the pair's FSC), mask-packed cross-half costs, coarse-to-fine label selection with ordered-label Potts
smoothing, composition of the filtered pair, the local-resolution map, the NU-evidence envelope and
the evidence state. It filters the one even/odd pair it is given. When NU runs, which pair is base,
auxiliary or apply pair, which mask plays which role and what is written belong to the callers:
`nonuniform_filter_state` (`volume/simple_nu_state_filter.f90`, both backends), `exec_nu_filt3D`
(`simple_commanders_resolest.f90`), `simple_commanders_postprocess_nu.f90`. The workflow contract is
skill `simple-nonuniform-regularization`.

## Rules

Silent ones first: nothing reports these.

- One sequence per process at a time, never two pairs concurrently: a rule, not a limitation to lift
  (the state is module-level, `simple_nu_filter_vars`). Order: `setup_nu_dmats` ->
  `set_nu_evidence_null_shell` -> envelope -> `set_nu_solvent_envelope` -> `optimize_nu_cutoff_finds`
  -> `build_nu_evidence_state` -> `nu_filter_vols` -> outputs -> `cleanup_nu_filter`. A clamp set
  after `optimize_nu_cutoff_finds` leaves the labels unclamped, a null shell set after the envelope
  leaves it on the mixture null; other steps out of order stop (`THROW_HARD`).
- `cleanup_nu_filter` after every sequence, on every path: `setup_nu_dmats` stops when the previous
  pair's clamp or null shell is still set (it does not clear them).
- The source pair of `nu_filter_vols` stays alive and unchanged from setup to composition: it is
  filtered again per label, and the composed voxels must be the candidates scored. Keep
  `prepare_nu_candidate_fts`/`filter_nu_candidate` and `compose_nu_labels_from_volume` on one taper,
  transform and Butterworth filter. Another pair (`pcg_solvent=yes`) gets the base pair's labels.
- Label index runs coarse to fine; an auxiliary pair is appended last, only at or beyond the finest
  retained rung, at that rung's Potts coordinate. The level loop, the Potts cap (a voxel only moves
  coarser), the clamp (label 1) and the handoff (`get_nu_filter_bank_finest_lp`, the last label)
  depend on that order.
- The NU domain is the sphere `setup_nu_dmats` builds from `mskdiam`; callers never pass a density
  mask as support. Costs are mask-packed (`nu_mask_vox` order): nothing outside the sphere enters the
  noise level, smoothing or selection; it is composed at the coarsest label. Exact zero/zero voxels
  (a PCG density support) are kept out of the noise level and evidence calibration.
- Units: `lowpass_limits`, `mskdiam`, `fsc_res`, `amsklp`, envelope lengths in Angstrom; `cutoff_finds`
  are Fourier indices of the input box; the filter runs at the pair's box (`box_crop` in assembly).
  Without `fsc_res` the bank is uncapped (`nu_filt3D` only).
- The Potts beta is a real OpenMP reduction: label fields repeat at a fixed thread count, not across
  counts. Compare runs at equal `nthr`.
- `build_nu_evidence_state` needs `evidence_source=` at setup, no auxiliary pair, the setup pair
  itself (fingerprint) and the unary bank, so it runs before `nu_filter_vols`, which releases the
  bank; it stops otherwise. The envelope reads `nu_ev_base`/`nu_ev_best`, which last until
  `cleanup_nu_filter`; the evidence state outlives it.
- `src/utils/filter/simple_tent_smooth.f90` is held bit for bit to a reference by
  `test_tent_smooth_exact` here.
- A new timed phase follows `nu_rt_*` (reset in `init_nu_filter`); a new field also goes through
  `get_nu_filter_timings` and `nu_state_filter_timings` (`volume/`) into the `VOLASSEMBLE` report.

## Where to start

- `simple_nu_filter.f90`: interfaces, tuning constants and their reasons. Bodies: `_bank` (setup,
  selection), `_state` (init, cleanup, packing, clamp), `_potts`, `_apply`, `_stats`, `_envmask`,
  `_evidence`, `_sharpen`. Algorithm, mask roles, envelope policy: [references/algorithm.md](references/algorithm.md).

## Tests

- `build/bin/simple_test_exec test=unit_reconstruction suite=nonuniform_filtering` (about 3 s): closed
  form and noise level, the label field of a two-resolution phantom, composition, the evidence
  envelope, the tent identity, the label field repeated at every processor.
- `simple_exec prg=nu_filt3D` on a refined half-map pair before and after: `_nu_filt`, `_nu_locres`,
  the `>>> NU` lines. In refinement: the `VOLASSEMBLE` NU split, `>>> NU AUXILIARY MEMBER:`,
  `>>> NU MATCHING LOW-PASS HANDOFF:`.
- Behaviour: `simulated_workflow_6vxx`/`_1jxy` (NU from solve3D stage 5; a robustness profile),
  `cont_refine3D_1jxy suite=nu_handoff`.

## See also

`doc/policies/NU/nonuniform_filtering_policy.md` §8-9 and 12, `doc/policies/3D/automasking_policy.md`;
skills `simple-nonuniform-regularization`, `simple-main-volume`.
