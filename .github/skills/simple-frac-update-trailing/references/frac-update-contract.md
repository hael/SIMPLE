# Fractional Update Contract

## Stage Policy

`simple_abinitio_controller.f90` emits the per-stage `refine3D` command line:

- `update_frac` decides the target outer particle subset.
- `trail_rec` turns on the trailing reconstruction consumer path in later stages.
- `nspace` and stage search policy are search/reference concerns; they must not
  silently desynchronize the produced handoff artifact from what the next stage
  will consume.

## Particle Selection

`simple_matcher_smpl_and_lplims.f90` chooses the active subset for
non-probabilistic refinement:

- `sample4update_class` for balanced class-aware sampling.
- `sample4update_cnt` for low-update-count biased sampling.
- `sample4update_all` when fractional update is inactive.
- `sample4update_fillin` for late fill-in behavior.

`simple_commanders_prob.f90` owns the probabilistic pre-step:

- `prob_align` samples once and writes the subset.
- `prob_tab` and `refine3D_exec` reproduce that exact subset with
  `sample4update_reprod`.

`simple_oris_sampling.f90` and `simple_oris_getters.f90` define the handshake:

- `sampled` marks the current sampling round.
- `updatecnt` records persistent update history.
- `get_update_frac` computes the realized fraction from the latest sampled round
  over the active updated pool.

## Current Partial Artifact

`simple_strategy3D_matcher.f90` writes partition-local current updates:

- The current volume reconstruction path: partial even/odd reconstructions and rho files.
- Polar/obsfield alternatives must write their own equivalent current-update artifact.

This is not a final reference yet. It is the current fractional update,
partitioned for assembly.

In the online reconstruction path, those current updates are produced from the same
particle batch read used for matching. Batch construction keeps the raw images
needed by reconstruction; after assignment, reconstruction preparation consumes
those images. Do not replace this with a later full particle pass that re-reads
the stack to lower peak memory unless the performance policy explicitly changes.

## Assembly and Trailing Reconstruction

`simple_commanders_rec_distr.f90::exec_volassemble` owns the current volume assembly path:

- It reduces partial reconstructions.
- It always computes realized per-state fractions `f` (`get_state_update_fracs`);
  the applied map-update weight `u` equals `f` unless a single-state
  `ufrac_trec` override is provided.
- If `trail_rec=yes`, trailing happens in the accumulator domain, mirroring 2D
  class-average restoration. The previous artifact is the per-state chain
  `trailrec_stateNN_{even,odd}.mrc` + `rho_*` files + `trailrec_stateNN.txt`
  manifest: blended, unregularized e/o Fourier sums and sampling densities at
  full-dataset sampling mass `D`:

```text
trailed_sums_eo = (u/f) * current_partial_sums_eo + (1 - u) * chain_sums_eo
trailed_rho_eo  = (u/f) * current_partial_rho_eo  + (1 - u) * chain_rho_eo
```

  The `u/f` scaling of the current partials makes the restored current-map
  coefficient exactly `u` (the historical `ufrac_trec` meaning) and keeps the
  chain at full mass: `(u/f)*(f*D) + (1-u)*D = D`. With no override, `u = f`
  and the current partials enter at full accumulator weight.
- The blended accumulators are written back as the new chain BEFORE restoration,
  because ML regularization adds 1/tau^2 to rho in place.
- A single sampling-density correction after the blend restores the trailed
  halves, so each Fourier component is weighted by its accumulated sampling
  density. The FSC is estimated post-blend from the restored trailed halves and
  describes the artifact written to disk; the merged volume and nonuniform-filter
  inputs are derived from the same blended statistics.
- Bootstrap (chain absent): previous even/odd halfmaps from the command line are
  mandatory. That iteration keeps the legacy behavior — FSC prior from the
  previous halfmaps and a volume-domain blend of the restored halves — while the
  chain is seeded with the current partials scaled by `1/f`, so the stored chain
  carries full-dataset mass. Without that normalization a chain seeded at
  fractional mass `f` would give the next iteration an effective update weight
  of `f / (f + (1-f)*f)` — for `f=0.1`, a requested 10 percent update would act
  like ~53 percent. From the next iteration the chain path is authoritative and
  the halfmap inputs are ignored.
- Full-weight seeding: a stage-boundary full reconstruction writes the chain at
  full-dataset weight when the internal `trail_seed=yes` cline handshake is set
  (`simple_abinitio_utils.f90::calc_rec`, gated on the consuming stage's
  `trail_rec`), so the consuming trailing stage starts from complete statistics
  instead of warming up.
- Artifact-set integrity: the manifest is deleted before and rewritten after
  the four data files, recording per-component byte sizes, a generation
  counter, and provenance (box, sampling, particle population, state layout).
  `validate_trail_chain` accepts a chain only when the manifest parses, the
  provenance matches the current project, every component size matches, and
  the grid is compatible; any failure discards the complete set and re-seeds.
  This rejects mixed-generation remnants from interrupted writes and stale
  chains from other runs sharing a directory.
- Continued runs from another directory carry chains over as complete sets
  only, destination cleaned first and manifest copied last
  (`simple_refine3D_strategy.f90::carry_over_trail_rec_chains`).
- The recurrence, override weighting, and bootstrap-normalization contracts are
  covered by the deterministic `simple_test_exec prg=trail_rec_blend` test.

The chain files become the previous artifact for the next iteration or stage.
Their names deliberately avoid the `recvol_state` stem so partial-reconstruction
globs and cleanup never match them.

## Downsampling Compatibility

`read_gridding_pair_accumulators`, contained in `restore_state_from_parts` in
`simple_commanders_rec_distr.f90`, is the previous-artifact
compatibility mechanism, used for both partial reconstructions and the trailing
accumulator chain:

- It requires previous MRC and rho files for both even and odd halves.
- When `l_update_frac` is active and previous dimensions are smaller than current
  dimensions, it pads with zeros (the abinitio autoscale ramp).
- It rejects previous dimensions larger than the current dimensions. For the
  trailing chain, `trail_chain_available` checks dimensions before the read and
  discards + re-seeds a larger chain instead of failing the run.

This is a producer/reader contract for previous artifacts, not a license to
invent a new fallback representation.

## Minimal Mirror for Alternate Representations

To mirror this for `polar=obsfield`:

- Keep current update production in obsfield partials.
- Extract dense reprojection references from obsfield partials before assembly
  feedback operations.
- Keep previous handoff in dense reprojection references
  (`POLAR_REFS_even/odd`, with merged refs when available), not in sparse
  restored obsfield artifacts.
- Require previous artifacts when `trail_rec=yes`.
- Compute the same realized update fraction.
- Do not blend previous/current updates in the obsfield domain, and do not
  provide obsfield-domain APIs for assembly feedback operations.
- Blend previous and current in the dense reprojection model before FSC/FRC calculation.
- Calculate FSC/FRCs from the dense reprojection model after trailing, not
  directly from obsfields.
- Apply low-resolution even/odd insertion to the dense reprojection model only;
  do not mutate sparse obsfield artifacts or add obsfield-level lowres insertion
  APIs for this policy step.
- Generate the reprojection model to the next consumer's interpolation limit
  before FSC, trailing, and low-resolution insertion.
- After restoring obsfields, extract the dense reprojection model through an
  explicit restored-extraction API. Generic extraction APIs that accept
  `invtau2` and `prior_start` are for raw obsfield extraction only; using them
  on restored fields risks double-normalization.
- Handle Fourier-range compatibility by overlap: crop previous larger reference
  ranges to the current interpolation limit and pad previous smaller ranges with
  zero contribution outside the overlap.
- When `nspace_next > nspace`, remap previous dense references state-by-state
  onto the denser current projection grid, using closest previous projections in
  the same state. Reject previous per-state grids larger than the current one
  unless a deliberate coarsening policy is added.
