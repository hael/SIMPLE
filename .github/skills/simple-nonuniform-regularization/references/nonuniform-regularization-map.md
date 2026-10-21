# Nonuniform Regularization Map

## Public Policy

The user-facing selector is `filt_mode`:

- `none`
- `uniform`
- `fsc`
- `nonuniform`
- `nonuniform_lpset`

`filt_mode=nonuniform` activates NU-filtered volume products while preserving
ordinary half-map matching unless another policy changes the bandwidth.
`filt_mode=nonuniform_lpset` activates the same NU products and promotes the
selected NU bandwidth into LP-set matching with merged registration-reference
topology. The public policy is documented in
`doc/policies/NU/nonuniform_filtering_policy.md`.

`src/main/params/simple_parameters.f90` declares the value and
`simple_parameters_phases.f90` derives `l_nonuniform` and
`l_nonuniform_lpset`. Unlike `uniform` and `fsc`, plain `nonuniform` does not
set `l_lpauto`/`l_lpset`; filtering is handled when Cartesian state volumes are
assembled. `automsk=nu` requires an NU mode, and `automsk=tight` is rejected
with NU filtering. 2D workflows reject NU modes.

For staged `solve3D`, `simple_solve3D_controller.f90`
(`set_refine3D_filtering_policy`) enables NU filtering from `NU_FILTER_STAGE`
(or `nu_stage`), promotes staged `nonuniform` to `nonuniform_lpset` (solve3D
has no gold-standard stage), and does not enable NU on the class-average route.

## Execution Path

Owners:

- `volassemble_exec` in
  `src/main/strategies/parallelization/simple_volassemble_strategy.f90`, on
  the caller's `(params, build)`: the distributed `refine3D` master runs it on
  its own builder; `run_volassemble` (behind `commander_volassemble` in
  `simple_commanders_rec_distr.f90`) builds its own for command-line callers
- the PCG master in `simple_rec3D_pcg_strategy.f90`, whose project handoff is
  `filter_pcg_nonuniform_maps` in the same strategy module
- `nonuniform_filter_state` in `src/main/volume/simple_nu_state_filter.f90`,
  shared by both backends: the per-state NU competition

High-level order for each state in gridding assembly:

1. Reduce the partition-local even/odd sums (trailing blend when active).
2. Restore the even/odd pair (`restore_gridding_pair`: deapodized, soft sphere
   at `msk_crop` once), keep the halves as written as the NU inputs (and the
   ML-regularized pair as the auxiliary inputs when `ml_reg=yes`), and write
   the FSC and resolution diagnostics.
3. Restore and write the merged state volume.
4. `nonuniform_filter_state`:
   - under active `automsk`, build the conservative density envelope of the
     base pair (`automask3D` at `envmsklp`); `automsk=nu` writes it as
     `automask3D_stateNN.mrc`
   - `setup_nu_dmats` on the base pair with `mskdiam` and `fsc_res=res0143`,
     the ML pair offered as the auxiliary member
   - under active `automsk`, write the NU-evidence envelope
     (`write_nu_evidence_envmask`) and arm the background clamp
     (`set_nu_solvent_envelope`): the valid evidence envelope under
     `automsk=nu`, else the density envelope
   - `optimize_nu_cutoff_finds`
   - `nu_filter_vols` on the base pair (or the apply pair), then kill the base
     pair
   - multiply the references by the active envelope under active `automsk`
   - write `_nu_filt` even, odd, merged and `_nu_locres`
   - record the handoff (`get_nu_filter_bank_finest_lp`) and `cleanup_nu_filter`
5. `update_project_nu_alignment_lowpass`: the finest valid handoff over the
   populated states becomes the particles' `lp` and `lp_est`.

## Support and Automask Policy

NU support is not configurable. `setup_nu_dmats` accepts `mskdiam` in
Angstrom and constructs spherical support internally. Callers cannot supply a
density automask or arbitrary logical envelope.

Mask artifact compatibility is
`state_mask_is_compatible(mskfile_state, box, smpd, exists, compatible)` in
`src/main/volume/simple_vol_pproc_policy.f90`: all three dimensions match the
box and the sampling distance matches within tolerance. It protects reused
mask artifacts (the lag-one NU-evidence mask of the gridding `envfsc` in
volume assembly, postprocess, the solve3D final reconstruction), not NU
support.

Artifacts:

- `automask3D_stateNN.mrc`: density envelope
- `nu_envmask3D_stateNN.mrc`: NU-evidence envelope

## Numerical Filter Module

Implementation owner:

- `src/main/nu_filt/simple_nu_filter.f90`
- `src/main/nu_filt/simple_nu_filter_*.f90` (working state in
  `simple_nu_filter_vars.f90`)

Required lifecycle:

```fortran
call setup_nu_dmats(vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, fsc_res=res0143)
call optimize_nu_cutoff_finds()
call nu_filter_vols(vol_even_nu, vol_odd_nu, vol_src_even, vol_src_odd)
call cleanup_nu_filter()
```

Static low-pass ladder (`lowpass_limits`):

- `20., 15., 12., 10., 8., 6., 5., 4.` Angstrom

Given `fsc_res` the ladder keeps the rungs at or coarser than
`fsc_res/NU_BANK_FSC_HEADROOM` (1.5), at least two; without it (`nu_filt3D`)
it is uncapped. There is no high-resolution extension.

## Candidate Selection Model

`setup_nu_dmats`:

- validates input even/odd dimensions and sampling
- validates positive `mskdiam` and constructs nonempty spherical support
- builds Butterworth filters for the retained bank
- appends the auxiliary pair as the last label when its resolution is at or
  beyond the finest retained rung (a stashed copy), else leaves it out
- estimates one global noise level `sigma_0` (`image%nu_objective_noise_scale`:
  Gaussian-scaled MAD of the raw even-odd difference over the observed
  support, exact zero/zero voxels excluded)
- filters each candidate pair from the transforms of the edge-tapered halves
  and computes its cost at once: the squared cross-half prediction error
  `((E - O_c)^2 + (E_c - O)^2) / (2 sigma_0^2)` (`image%nu_objective`); the
  auxiliary cost comes from the pair itself; no candidate volume is kept
- stores the raw costs (`raw_dmats_mask`), the costs smoothed at each
  candidate's own radius (`dmats_mask`) and the smoothing norm of that radius
  (`norms_mask`), all mask-packed

`optimize_nu_cutoff_finds`:

- selects coarse to fine: at each level the incumbent and the candidate are
  smoothed at the candidate's radius from the raw costs, and the candidate
  wins only with a strictly lower cost; the candidate's own column and norm
  are the setup's
- applies the background clamp (when armed) before and after the Potts sweeps
- runs ordered-label Potts smoothing, which never moves a voxel finer than the
  label it entered with, except for degenerate exits such as a single label or
  numerical-zero prior strength
- stores the winners in `filtmap`; the unary bank stays allocated for
  `build_nu_evidence_state`

The ordered-label prior uses ladder coordinates; the auxiliary label shares
the finest rung's coordinate. Adjacent coordinate steps are tolerated; larger
jumps are penalized with a linear-quadratic hinge.

`nu_filter_vols`:

- releases the unary bank
- filters each half of the source pair per label (the coarsest seeds the whole
  grid) and scatters the voxels by `filtmap`
- copies auxiliary voxels where the auxiliary label is selected

## ML-Regularized Coupling

When `params%l_ml_reg` is true (`nu_aux_member`), volume assembly uses:

- base input: the `_unfil` even/odd pair
- auxiliary member: the ML-regularized even/odd pair, at the state's
  FSC=0.143 resolution

The auxiliary pair joins the bank beside the finest retained rung once its
resolution is at or beyond it, and competes voxel by voxel under the same
unary, smoothing and prior. It never replaces a rung.

## Outputs and Consumers

Output suffixes:

- `NUFILT_SUFFIX = '_nu_filt'`
- `NULOCRES_SUFFIX = '_nu_locres'`

For each state, `nonuniform_filter_state` writes:

- `vol_state_even_nu_filt.mrc`
- `vol_state_odd_nu_filt.mrc`
- `vol_state_nu_filt.mrc`
- `vol_state_nu_locres.mrc`

Names are built via `add2fbody(..., MRC_EXT, NUFILT_SUFFIX)` and
`NULOCRES_SUFFIX`.

Matcher reference loading is in:

- `src/main/strategies/search/simple_matcher_refvol_utils.f90`

Plain `nonuniform` tries `_nu_filt` even/odd references first, then falls back
to regular even/odd files. `nonuniform_lpset` with active LP-set matching uses
merged registration-reference topology and prefers the merged `_nu_filt`
product. Fallback matters before the first assembled nonuniform products exist.

## Tests and Manual Checks

The unit tests are `src/main/nu_filt/simple_nu_filter_tester.f90`, the
sub-suite `nonuniform filtering` of `unit_reconstruction`: the objective's
closed form, the noise level, the ranking's independence of that level, the
label field on a two-resolution phantom, the per-label composition, the
evidence state and envelope on the neutral fixture (a sphere of band-limited
common signal in a noise support; it reports envelope recall and solvent
false-positive rate), the tent smoothing against a column-by-column reference
bit for bit, and the label field repeated at every processor. Run the filter
end to end with `simple_exec prg=nu_filt3D` on a refined half-map pair.

For code changes, consider:

- `simple_exec prg=nu_filt3D` before and after the change on the same maps
- the `VOLASSEMBLE` bench report's NU split before and after
- checking nonuniform mode still falls back to regular references before
  `_nu_filt` files exist
- checking multi-state automask naming and compatibility
- checking `ml_reg=yes` auxiliary admission (`>>> NU AUXILIARY MEMBER:`) and
  the handoff line (`>>> NU MATCHING LOW-PASS HANDOFF:`)

## Current Limitations

- The module uses module-level state, which makes concurrent use awkward.
- The selected-cutoff map and summary statistics are exposed through module
  state and logging rather than explicit outputs.
