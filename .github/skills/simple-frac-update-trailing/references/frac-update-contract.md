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

`simple_commanders_rec_distr.f90::exec_cartesian_assembly` owns the current volume assembly path:

- It reduces partial reconstructions.
- It computes FSC and sampling-density correction.
- If `trail_rec=yes`, it requires previous even/odd volumes from the command line.
- It computes `update_frac_trail_rec` from `ufrac_trec` if provided; otherwise from `get_update_frac`.
- It blends current and previous even/odd volumes:

```text
trailed_even = update_frac_trail_rec * current_even + (1 - update_frac_trail_rec) * previous_even
trailed_odd  = update_frac_trail_rec * current_odd  + (1 - update_frac_trail_rec) * previous_odd
```

The output halfmaps become the previous artifact for the next iteration or stage.

## Downsampling Compatibility

`simple_reconstructor_eo%read_eos_parallel_io` is the current previous-halfmap compatibility mechanism:

- It requires previous MRC and rho files for both even and odd halves.
- When `l_update_frac` is active and previous dimensions are smaller than current
  dimensions, it pads with zeros.
- It rejects previous dimensions larger than the current dimensions.

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
