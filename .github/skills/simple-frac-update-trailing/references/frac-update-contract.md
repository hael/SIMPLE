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
- It computes `update_frac_trail_rec` from `ufrac_trec` if provided (single-state);
  otherwise from realized per-state fractions (`get_state_update_fracs`).
- If `trail_rec=yes`, trailing happens in the accumulator domain, mirroring 2D
  class-average restoration. The previous artifact is the per-state chain
  `trailrec_stateNN_{even,odd}.mrc` + `rho_*` files: blended, unregularized e/o
  Fourier sums and sampling densities:

```text
trailed_sums_eo = current_partial_sums_eo + (1 - update_frac_trail_rec) * chain_sums_eo
trailed_rho_eo  = current_partial_rho_eo  + (1 - update_frac_trail_rec) * chain_rho_eo
```

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
  chain is seeded with the current sums. From the next iteration the chain path
  is authoritative and the halfmap inputs are ignored.
- Full-weight seeding: a stage-boundary full reconstruction writes the chain at
  full-dataset weight when the internal `trail_seed=yes` cline handshake is set
  (`simple_abinitio_utils.f90::calc_rec`, gated on the consuming stage's
  `trail_rec`), so the consuming trailing stage starts from complete statistics
  instead of warming up.
- Continued runs from another directory carry the chain files over alongside the
  FSCs (`simple_refine3D_strategy.f90::carry_over_trail_rec_chains`).

The chain files become the previous artifact for the next iteration or stage.
Their names deliberately avoid the `recvol_state` stem so partial-reconstruction
globs and cleanup never match them.

## Downsampling Compatibility

`simple_reconstructor_eo%read_eos_parallel_io` is the previous-artifact
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
