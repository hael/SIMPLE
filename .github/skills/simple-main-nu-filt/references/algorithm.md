# NU filter: algorithm, masks and fixed policy

Reference for `src/main/nu_filt`; the rules are in the guide (`../SKILL.md`). The design and its
rationale are `doc/policies/NU/nonuniform_filtering_policy.md` §8-9 and 12; the assembly order and
the ML coupling are skill `simple-nonuniform-regularization`.

## Bank (`setup_nu_dmats`, `init_nu_filter`)

- Static ladder `lowpass_limits` = 20, 15, 12, 10, 8, 6, 5, 4 A (`simple_nu_filter_vars`). With
  `fsc_res` the rungs at or coarser than `fsc_res/NU_BANK_FSC_HEADROOM` (1.5) are kept, at least two;
  without it (`nu_filt3D`) the ladder is uncapped. At most `NU_DMAT_CANDIDATE_CAP` (24) labels.
- The auxiliary pair (the ML-regularised halves, `ml_reg=yes`; in assembly its resolution is the
  state's FSC=0.143) is stashed as a copy and appended as the last label only when its Fourier
  index is at or beyond the finest retained rung; it shares that rung's Potts coordinate. Otherwise
  the rungs compete alone (`>>> NU AUXILIARY MEMBER:` says which).
- Support: the sphere of `mskdiam` (A) on the input grid, packed once (`nu_mask_vox`).
- Each candidate pair is filtered from the transforms of the edge-tapered halves
  (`prepare_nu_candidate_fts`, `filter_nu_candidate`: taper `COSMSKHALFWIDTH`, Butterworth) and
  scored at once: the squared cross-half prediction error at one global noise level
  (`image%nu_objective`, `nu_objective_noise_scale` over the observed support). No candidate volume
  is kept. Stored mask-packed: raw costs `raw_dmats_mask`, costs smoothed at each candidate's own
  radius `dmats_mask` (radius 0.5 x 3 x lowpass, at most 30 A), the smoothing norm `norms_mask`, and
  the raw coarsest and best costs `nu_ev_base`/`nu_ev_best` for the envelope.

## Labels (`optimize_nu_cutoff_finds`, `_potts`)

- Coarse to fine, like for like: at level `l` the incumbent and candidate `l` are both smoothed at
  `l`'s radius from the raw costs; `l` wins only when strictly cheaper. `l`'s own column is the
  setup's `dmats_mask` and its norm the setup's `norms_mask` (`load_nu_smooth_norm`): neither is
  smoothed again.
- The background clamp (armed by `set_nu_solvent_envelope`) sets support voxels outside the envelope
  to label 1 before and after the Potts sweeps; the sweeps re-optimise from the costs and would
  erode a one-time override.
- Ordered-label Potts: 26-neighbourhood, 8 colours walked as packed per-colour voxel lists, at most
  6 sweeps; beta = 2 x the mean gap between the best and second-best cost (a real reduction); a
  voxel moves only to a label no finer than the one it entered with. Adjacent ladder steps are free,
  larger jumps pay a linear-quadratic hinge. It is part of the algorithm, not a user mode.
- The unary bank stays allocated for `build_nu_evidence_state`; `nu_filter_vols` releases it.

## Composition (`nu_filter_vols`, `_apply`)

Each half of the source pair is tapered, transformed and filtered per label; the coarsest filter
seeds the whole grid, support voxels take their own label's filter, auxiliary voxels are copied from
the stashed pair. An apply pair other than the setup pair (`pcg_solvent=yes`) receives the setup
pair's label field.

## Masks

- The sphere of `mskdiam` is the only NU domain; callers never pass a density mask or a logical
  envelope as support.
- The density envelope (`automask3D_stateNN.mrc`, `automask3D` at `envmsklp`) is independent of NU
  support: it feeds `envfsc` and, under `automsk=yes` or as the `automsk=nu` fallback, is the
  background clamp and multiplies the `_nu_filt` references.
- The NU-evidence envelope (`nu_envmask3D_stateNN.mrc`, `write_nu_evidence_envmask`) is a
  diagnostic under `automsk=yes`; under `automsk=nu`, when valid, it is the clamp and reference mask
  and its lag-one file the gridding `envfsc` mask (density fallback). It is never NU support and never
  the PCG solve support. Which mask plays which role is decided in `nonuniform_filter_state`.
- An envelope-constrained base pair (PCG density support) gets a Euclidean null shell
  (`set_nu_evidence_null_shell`: labels free on the density envelope, null statistics on its
  dilation ring).

## Standalone envelope (`nu_filt3D`) and evidence

- Two shape controls only: `nu_msk_sig` (threshold, in null MADs) and `amsklp` (evidence scale, A).
  Fixed policy beside the constants in `simple_nu_filter.f90`: absolute evidence
  (`NU_ENVMASK_RELATIVE`), zero density weight, MRF beta 1, smallest kept component 0.1 of the
  largest; growth 1 A and a 6 A cosine edge, converted to the nearest voxel count at the input
  sampling (at least one).
- The evidence state (`build_nu_evidence_state`, `nu_evidence_state`) is a compact immutable record
  of one setup (selected labels and cutoffs, uncertainty, nested band support at 20, 12, 8, 5 A plus
  appended bands, summary, identity). `assert_nu_evidence_replay_ready` checks it (null fraction
  0.01-0.90); `nu_evidence_sharpen_vol` (`postprocess_nu`) consumes it after `cleanup_nu_filter`.

## Timings and logs

Under `L_BENCH_GLOB` the sequence accumulates `nu_rt_filtering`, `_objective`, `_optimise` (of it
`_opt_levels` and `_opt_potts`) and `_compose`, reset in `init_nu_filter`; `get_nu_filter_timings`
feeds `nu_state_filter_timings` and the `VOLASSEMBLE` report. `set_nu_filter_report` and
`NU_DEV_OUTPUT` persist across sequences.
