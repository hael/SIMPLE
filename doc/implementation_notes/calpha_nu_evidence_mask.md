# NU-Evidence Resolution Mask for C-alpha Detection

> Status 2026-09-16: the `nu_refine` shell walk referred to below is retired; the evidence state is built from the generated ladder (`nonuniform_filtering_policy.md` section 8). The `nu_refine` control no longer exists.

## Status

Design, revised 2026-09-11 after review (section 11). Not implemented.
Decisions taken in the revision:

- input is the unregularized half pair, `vol1` (odd) and `vol2` (even);
- the workflow is its own commander and `simple_exec` program (section 2);
- the atom response is computed on the NU-filtered merged map (section 4);
- the target support `q_T` is read from the compact NU evidence state, not from
  a second evidence pathway (section 3);
- the objective is written in physical units and the result may consist of a
  bounded number of independently rooted components (section 5);
- the ROI is also the support for absolute B-factor scaling of the map ahead of
  initial model building (section 10);
- the lifecycle template is `postprocess_nu`, validated on PfCRT 2026-09-11:
  run on the `_unfil` half pair it produced the best PfCRT map to date. The
  input must be the unregularized pair; a regularized pair flattens the
  evidence margin (see the `nu_filt3D` warning) and the result is unusable.

This is a model-building region-of-interest mask. It is not a replacement for
the NU filter-field envelope (`nu_envmask3D`), the density automask, the
spherical NU objective support, or the FSC mask. It supports the C-alpha
detection work in [`calpha_finder.md`](calpha_finder.md).

## 1. Problem and contract

Automatic atomic model building first needs the regions of a high-resolution
map that contain reproducible detail. The masking input is an NU resolution
floor `T_NU` in Angstroms, for example 5 A. It is not the resolution at which
C-alpha building is expected to work and not the resolution of the map; it is
a permissive quality boundary that rejects parts of an otherwise finer map
that have degraded too far for the model-building workflow. A voxel meets the
floor when its supported resolution is `T_NU` or finer, so smaller values are
better.

A hard threshold on a local-resolution map is insufficient: it fragments,
discards confidence, and cannot distinguish a short weak link between two
strong regions from a large unsupported region. The contract is instead:

> Find contiguous regions supported as strongly as possible at `T_NU` or
> finer, allowing limited lower-resolution material where the evidence cost is
> justified by connectivity.

The floor is a soft evidence constraint. Connectivity is a topological
property of each returned region; it never forces the solver to build a bridge
through unsupported density. A remote island is dropped when connecting it
costs more than the evidence inside it is worth, and a structure made of
genuinely separate parts is returned as separate components.

## 2. Workflow: `simple_exec prg=nu_roi3D`

The evidence state (`nu_evidence_state`) is in-memory only; nothing serializes
it, and `detect_calpha` searches a single `vol1`. The ROI therefore has to be
built where the evidence exists, in a dedicated commander that owns the whole
NU lifecycle for one half pair.

Registration mirrors `postprocess_nu`: `commander_nu_roi3D` in a new
`simple_commanders_nu_roi3D.f90`, a `ui_program` entry next to `nu_filt3D` /
`postprocess_nu` in `simple_ui_filter.f90`, and a `case('nu_roi3D')` in
`simple_exec_filter.f90`. The commander owns I/O, parameters and lifecycle; a
new domain object (`nu_roi`, model-building domain) owns the evidence
conversion, atom response, region assembly and the scaling of section 10.

Inputs:

- `vol1`, `vol2`: the unregularized (`_unfil`) odd / even half maps, same
  grid. The commander should refuse or at least warn when the evidence margin
  is flat, as `nu_filt3D` does, since that is the signature of a regularized
  pair.
- `smpd`, `mskdiam`: as for `postprocess_nu`.
- `nu_refine`: evidence-gated shell walk, default `yes` as in `postprocess_nu`
  (cost is paid once per map).
- `nu_floor`: `T_NU` in Angstroms, default 5.
- `outvol`: soft ROI, default `vol_nu_roi.mrc`.
- `nthr`.

Lifecycle, in order:

1. `setup_nu_dmats(even, odd, mskdiam, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE)`,
   `optimize_nu_cutoff_finds`, and the accepted shell walk when `nu_refine=yes`
   (`extend_nu_filter_highres_shells`). One evidence identity, exactly the
   `postprocess_nu` lifecycle.
2. Merge `(even + odd) / 2` and NU-filter it while the bank is live:
   `nu_filter_vol(merged, merged_nu)`. This is the atom-response work map.
3. `build_nu_evidence_state(even, odd, state)`, `cleanup_nu_filter`,
   `assert_nu_evidence_replay_ready(state)`, `print_nu_evidence_summary`.
4. `unpack_nu_evidence_state(state, selected_cutoff=, band_support=)`; the
   summary provides `band_limits` and `n_bands`.
5. Atom response on `merged_nu` (section 4), target support and utility
   (section 3), region assembly (section 5), products (section 6).
6. Absolute B-factor scaling inside the ROI (section 10) and the
   model-building map.

Products (stem from `outvol`):

- `<stem>.mrc`: soft ROI in [0,1];
- `<stem>_class.mrc`: 0 outside, 1 connector, 2 core;
- `<stem>_qT.mrc`: target support `q_T`;
- `<stem>_atom.mrc`: normalized atom-support field `a_C`;
- `<stem>_cost.mrc`: accumulated connector cost (bottleneck confidence);
- `<stem>_nufilt.mrc`: the NU-filtered merged map the atom response was
  computed on;
- `<stem>_mb.mrc`: the model-building map, B-scaled inside the ROI
  (section 10), the map `detect_calpha` should search;
- `<stem>_summary.txt`: section 6.

Consumer: `detect_calpha` gains an optional `mskfile` (the soft ROI) that
gates the search and a `_class` volume that ranks candidates. Nothing else in
`detect_calpha` changes.

## 3. Target support from the compact evidence state

`build_nu_evidence_state` already computes, per packed support voxel `v` and
per band `b`,

```text
band_support(v, b) = sum over signal candidates with lp <= band_limit(b)
                     of p(candidate | v)
```

where `p` is the temperature-calibrated posterior over {null, candidates}
(`candidate_measure * exp(-(E - E_min) / temperature)`, normalized). That is
exactly

```text
q_T(v) = P(local detail is supported through T | NU evidence)
```

for `T` in the band ladder: the static 20 / 12 / 8 / 5 A bands, extended under
`nu_refine=yes` at ratio `NU_EVIDENCE_BAND_RATIO = 0.64` (3.2, 2.05 A, ...)
over challenger-accepted shells only, pruned finest-first when a band earns
less than `NU_EVIDENCE_MIN_BAND_SUPPORT`. The competing coarser labels and the
explicit null are integrated out already. The selected label
(`selected_cutoff`) is an argmin, discards confidence and jumps between
neighboring candidates; it is used only to index the atom kernel bank
(section 4) and the pseudo-atom B-factors (section 10), never as the mask
statistic.

`T_NU` is snapped to the nearest retained band in spatial frequency
(`1 / T`), and the snap is logged (`>>> NU ROI: floor 5.00 A snapped to band
5.00 A`). If `T_NU` is finer than the finest retained band the commander stops
with an error naming the finest band; the evidence cannot support the request.
Interpolation between bands in `1 / T` is a later refinement, not needed for
the first slice.

Two properties of the compact posterior are accepted as-is:

- it is Potts-conditioned: `evidence_site_energy` adds `beta * pair_cost`
  against the neighbors' selected labels, so `q_T` already carries the label
  smoothing of the NU field and is a conditional given the MAP neighborhood
  (somewhat over-confident). The sharpening path uses the same posterior. The
  only additional smoothing is one physical Gaussian of sigma `T_NU / 2`
  applied to `q_T` before seeding, so that seeds are resolution-element sized
  rather than voxel sized.
- `uncertainty` (normalized entropy over all candidates including the null) is
  already integrated into `q_T` and does not enter the objective. It is a
  diagnostic only (section 6).

Domain. The state lives on the packed spherical `mskdiam` support
(`n_nu_mask`, `nu_mask_vox`). Everywhere else `q_T = 0`. Solvent-clamped
voxels carry `band_support(:,1) = 1` and nothing finer, so `q_T = 0` for any
`T` finer than the coarsest band. Unobserved voxels (no half-map pair) are the
explicit null with zero band support. All three cases fall out of the same
utility

```text
u_T(v) = log((q_T(v) + eps) / (1 - q_T(v) + eps))
```

as confident exclusion, `u_T = log(eps)`. The ROI is produced on the grid of
`vol1` / `vol2`; if the model-building work map is on another grid (refinement
crop versus full box) the ROI is resampled by the consumer, not here.

## 4. Atom response on the NU-filtered merged map

The nanoparticle path provides the template primitive, not the operator.
`phasecorr_one_atom` is phase correlation (spectrally whitened) against a
carbon convolved at `lp = 2 * smpd` with an `8 * smpd` cutoff; on a protein
map at 3-4 A that is a high-pass noise amplifier and is not reused. What is
reused is `atoms%convolve(vol, cutoff, lp)`, which already builds an
electron-scattering-factor carbon at a requested resolution, and the
normalized cross-correlation machinery of `calpha_finder` (`image%ccf_into`,
three correlations for the weighted data sum, squared-data sum and target-data
sum).

The work map is the NU-filtered merged map `merged_nu` from step 2 of the
lifecycle, for two reasons. First, its bandwidth is known per voxel: the NU
filter has low-passed every voxel at its `selected_cutoff`, so the carbon
kernel can be matched locally instead of to one global estimate. Second, the
NU filter has already removed detail the half pair does not reproduce, so the
atom response cannot reward noise that the evidence rejected; on the
unfiltered merged map it would, and the ~100k score-eligible voxels at
`thres = 0.25` in the 2026-09-10 `detect_calpha` baseline are exactly that
population.

Computation:

1. Kernel bank: one carbon per retained candidate cutoff finer than or equal
   to `T_NU` (`atoms%convolve` with `lp` = that cutoff, cutoff radius `4 * lp`);
   coarser voxels get no atom response. This resolves the matched-filter
   question: the bank is indexed by the evidence, not inferred from the map.
2. For each kernel, weighted NCC of `merged_nu` inside a sphere of radius
   `lp` with a cosine taper, exactly as the C-alpha target is scored. The
   response at voxel `v` is taken from the kernel matching
   `selected_cutoff(v)`.
3. Robust normalization: median and MAD of the response over the NU null
   population (packed support outside the evidence envelope), converted to a
   z-score and clipped at 0 below. The evidence envelope is already computed
   in the lifecycle (`nu_evidence_envelope`), so no new solvent estimate is
   introduced.
4. Atom-support field `a_C`: the clipped z-score pooled (mean) over a sphere
   of radius `T_NU / 2`, so the objective rewards regions populated by
   plausible atomic densities rather than isolated correlation peaks.
   `u_C = log((a_C + eps) / (1 + eps))`, non-positive to zero, so the atom
   term can only be a reward relative to solvent, never a penalty.

Role. The NU floor remains the gate: a strong atom response with no
reproducible support at `T_NU` never establishes core, and `a_C` never enters
the seed selection. The atom term does two things: it raises the reward of
regions that are both reproducible and atomic, and it makes a short
connector through slightly sub-floor evidence that still carries atomic
density cheaper than one through structureless density. The circularity
caveat from the review stands: `a_C` is a weaker relative of the detector the
ROI serves, so its contribution must be established by the NU-only ablation
in section 9, and `alpha = 0` must remain a supported setting. Half-map atom
responses are not computed; their reproducibility is what `q_T` already
measures.

## 5. Connected-region objective

Let `M` be a binary region. Everything is expressed in physical units so that
the balance of the terms does not drift with sampling: volume integrals carry
`smpd^3`, area integrals `smpd^2`, path costs are per Angstrom of path.

```text
maximize  smpd^3 * sum(v in M) [u_T(v) + alpha * u_C(v)]
          - lambda_boundary * smpd^2 * boundary_area(M)

subject to M = union of at most K components, each connected to its own root.
```

Boundary roughness is not a separate term. The morphological open/close at a
physical radius `T_NU / 2` in the assembly (step 6 below) is the boundary
regularizer; it is well defined, sampling-independent and already needed for
the delivered mask.

The per-voxel exclusion cost is bounded by `log(1 / eps)`, so an island with
integrated reward `R` can pay for a connector of at most
`R / (smpd^3 * log(1 / eps))` voxels. The maximum bridge length is therefore
not an independent parameter: it is derived from `eps` and the island reward,
and the commander reports the implied bound in Angstroms and the longest
connector actually accepted. `eps` is the one knob, and it is stated as a
minimum accepted support probability, not as a length.

Solver: evidence-weighted geodesic assembly, per component.

1. Seeds are taken from the finest retained band, not from `q_{T_NU}`: with a
   5 A floor on a 3.5 A map `q_5` is about 1 over the whole protein and does
   not rank. Seed regions are connected components of
   `q_finest (smoothed) >= p_seed`; the root is the seed with the largest
   integrated reward.
2. Every other positive-reward region (`u_T + alpha * u_C > 0`, connected
   components after the same smoothing) is a candidate island.
3. Minimum-cost geodesic path from each island to the current component; the
   edge cost is `max(0, -(u_T + alpha * u_C))` per Angstrom, zero through
   positive-reward voxels.
4. Attach the island with the best `reward - (path + boundary) cost` if that
   is positive; repeat from step 3 with the enlarged component.
5. When no island can be attached, the component is closed. Remove it, and
   start the next component from the best remaining seed while that seed's
   integrated reward exceeds `R_min` and fewer than `K` components exist.
6. Per component: open/close at radius `T_NU / 2` with connectivity
   preserved, then a growth margin of `T_NU / 2` with a cosine edge of width
   `T_NU / 2` for the delivered soft mask.

An exact maximum-weight connected subgraph on the voxel graph is not
practical; a prize-collecting Steiner tree on supervoxels is a possible later
solver. A binary MRF followed by keeping the largest component is not
equivalent: it can discard islands but cannot decide whether a bridge is worth
paying for.

## 6. Mask semantics and outputs

Three classes are retained even though the primary product is one soft mask:

- **core**: `q_finest >= p_seed` after smoothing (directly supported at the
  finest evidenced band);
- **connector**: inside `M` but not core, i.e. supported only through `T_NU`
  or included to preserve connectivity;
- **outside**.

Connectors must not be reported as if they achieved the requested resolution.
The summary (`<stem>_summary.txt`) contains: `T_NU` requested and snapped,
the band ladder, number of components, per-component core and connector
volumes in A^3, fraction of the mask that is core, longest accepted connector
in A, weakest accepted bottleneck (minimum `q_T` along any accepted path), the
implied bridge bound from `eps`, `alpha`, the evidence-state summary (null
fraction, uncertain fraction, supported fraction per band), and the scaling
record of section 10 (target B, fitted B inside the ROI, applied ΔB).

`detect_calpha` searches the full soft ROI on `<stem>_mb.mrc`, ranks core
detections above connector detections, and treats connector detections
conservatively; the `_cost` volume tells it how far a detection sits behind a
weak bottleneck.

## 7. Scientific and workflow boundaries

Mask connectivity is not backbone connectivity. A contiguous region can
contain merged helices, branches, membrane, nucleic acid, ligand, or density
that does not define a traceable C-alpha path. The ROI is a search region and
confidence prior for C-alpha detection, not a trace.

Hard ownership constraints:

- the ROI is never the NU objective support; `setup_nu_dmats` keeps the
  spherical `mskdiam` support and its solvent / null population;
- the ROI is never used for FSC estimation or phase-randomized correction; it
  is selected from half-map reproducibility at the target frequency and is
  circular for that purpose;
- `nu_envmask3D_stateNN.mrc` is not overloaded; it defines the NU filter-field
  background and answers a different question;
- the unmasked half maps remain the primary artifacts; `nu_roi3D` writes only
  the products listed in section 2.

## 8. Parameter policy

Public: `nu_floor` (`T_NU`, A). Everything else is either derived from the
evidence or fixed by policy in the domain object:

- `p_seed`: seed probability on the finest band, policy (start at 0.8);
- `eps`: minimum accepted support probability, policy (start at 0.02), with
  the implied bridge bound reported;
- `alpha`: atom-term weight, policy (start at 0.5 with `alpha = 0` as the
  ablation);
- `R_min`, `K`: minimum integrated reward for a new root and the component
  cap, policy (start at `K = 4`, `R_min` = 5 % of the root's reward);
- `lambda_boundary`: policy, per A^2;
- smoothing sigma, open/close radius, growth margin and edge width are all
  `T_NU / 2`; the atom pooling radius is `T_NU / 2`; the kernel bank and
  cutoff radii follow the retained candidate cutoffs;
- scaling target B (section 10): policy, 0 A^2 by default.

No parameter is expressed in voxels. Calibration of the policy constants is
part of the validation below; they are promoted to public inputs only if the
validation shows they need to vary between maps.

## 9. Validation plan

Mask quality is validated separately from downstream detection quality.

Fixture. The `detect_calpha_molecules` benchmark builds one
uniform-resolution map from 6VXX / 1JYX with no half pair, so no NU evidence
can be built on it. A new fixture generates a half pair from the same
coordinates with spatially varying resolution: a locally varying B-factor
field plus independent noise in the two halves, with controlled cases for a
short degraded gap between two well-resolved regions, a distant well-resolved
island behind a long unsupported gap, and two unrelated nearby chains.
`simple_test_nu_envmask` is the closest existing fixture and is the starting
point. Because the fixture is generated from coordinates with a known
B-factor field, it also validates the scaling of section 10: the fitted B
inside the ROI must recover the generating B.

Mask-level measurements:

- precision and recall of the `T_NU`-qualified region against the known
  local-resolution field;
- precision and recall of the atom response against known atomic positions
  on `merged_nu`;
- core and connector volumes;
- false bridges between the two unrelated chains (must be zero);
- recovery across the short gap; rejection of the distant island;
- number of components returned on the two-chain case (must be two);
- stability across sampling distance, box size, `T_NU`, and band
  granularity (`nu_refine=yes/no`);
- behavior when no seed passes `p_seed` or the null is uncalibrated
  (`assert_nu_evidence_replay_ready` failure): a clean error, no mask.

Model-building measurements: `detect_calpha` with and without the ROI, on the
NU-filtered map and on the B-scaled map, on the fixture and on representative
experimental maps, with the ablations NU-only (`alpha = 0`), combined, and NU
core with atom-assisted connectors. Candidate precision / recall are reported
separately in core, connector and excluded regions. The mask is useful only if
it improves the operating tradeoff or reduces search cost without hiding
recoverable C-alpha sites; the atom term is retained only if the combined
ablation beats NU-only; the scaling is retained only if detection on
`_mb.mrc` beats detection on `_nufilt.mrc`.

## 10. Absolute B-factor scaling for model building

If the ROI works it is also the right support for preparing the map for
initial atomic model building. The idea is LocScale's model-free mode
(`~/src/locscale`, `preprocessing/pipeline.py`): fill the mask with
pseudo-atoms, derive a reference spectrum, scale the map to it. Reading that
code shows which parts carry the information, and which SIMPLE already has.

What LocScale does: mask -> atom count (mask volume / 1.55, or MW / 13.14 amu)
-> pseudo-atoms (element `O`, minimum "bond length", gradient or
random-placement-with-kick solver) -> iterative servalcat ADP refinement
against the half maps -> B shifted so that p(B < 0.01) = 0 -> gemmi
electron-scattering reference map -> windowed local scaling
(`WindowedScaler`: per window, radial |ref| / |target| ratio applied bin by
bin, central voxel kept). The global variant `scale_profiles` is a Wilson-fit
ΔB between the Wilson cutoff and the FSC cutoff, applied as
`exp(ΔB s^2 / 4)`. Packing artefacts of randomly placed atoms are handled by
`measure_debye`, `generate_no_debye_profile` and blending with theoretical
profiles beyond the Wilson cutoff.

Four observations decide the SIMPLE design:

1. With uniform B the atom positions are nearly redundant for global scaling.
   Above the Wilson cutoff a dense pack of identical atoms has the radial
   spectrum `N |f(s)|^2 exp(-B s^2 / 2) S(s)`; the positions contribute only
   the mask-shape term below the Wilson cutoff and the packing term `S(s)`.
   The target profile is the analytic scattering-factor curve at a chosen B,
   which `atoms%convolve` already parameterizes. What this buys is real:
   `nu_evidence_sharpen_vol` fits one Guinier B to the map itself, a relative
   correction, whereas matching to a physical atom ensemble at a chosen B
   gives an absolute target, which is what model-building programs want.
2. The information LocScale adds is spatial and comes from the refined
   per-atom ADPs, not from the atoms. SIMPLE already holds that spatial
   information in the evidence state (`selected_cutoff`, `q_T`), half-map
   based and computed in this lifecycle. Pseudo-atom B-factors can therefore
   be set from the local NU cutoff instead of from servalcat.
3. Never scale beyond the evidence. Bin-by-bin ratio scaling to a noise-free
   reference amplifies whatever the target holds at frequencies without
   signal. The guard is the NU cutoff: scale only inside the local evidenced
   passband and keep the NU-filtered value elsewhere, the same principle as
   the sharpening v2 design after the PfCRT over-sharpening record.
4. Randomly packed atoms have a liquid-like `S(s)` and none of the 4.7 A
   beta / 10 A helix features a real protein has; a raw bin-by-bin ratio
   would divide those out of the map. Use a smooth reference (Wilson /
   Guinier fit) or LocScale's blending.

Design, two slices:

- **Slice 1, global absolute B inside the ROI (first).** No atoms placed.
  Radial amplitude profile of `merged_nu` inside the soft ROI, Guinier fit
  between the Wilson cutoff (from the ROI volume, LocScale's `find_wilson_cutoff`
  method) and the finest evidenced cutoff; the same fit on the analytic
  carbon scattering-factor profile at the target B (policy, 0 A^2); apply
  `exp(ΔB s^2 / 4)` to the ROI-masked map, cosine-tapered at each voxel's
  `selected_cutoff` so nothing beyond the local evidence is amplified. Product
  `<stem>_mb.mrc`; fitted B, target B and ΔB in the summary. This is a small
  addition to `nu_evidence_sharpen_vol`, which is validated on PfCRT: the
  first `_mb.mrc` is simply that sharpened map inside the ROI, and the
  absolute-B target is the refinement on top of it, to be compared against
  it on PfCRT before it replaces it.
- **Slice 2, NU-B pseudo-model with windowed scaling (later).** Pseudo-atoms
  placed on ROI core voxels by Poisson-disk sampling at 1.5 A minimum
  distance, count from the ROI volume (LocScale's empirical 1.55 correction as
  the starting point; the absolute scale cancels in B scaling); per-atom B
  from `selected_cutoff` through a monotone resolution-to-B map calibrated on
  the fixture; reference map from `atoms%convolve`; windowed radial scaling
  as in `WindowedScaler`, with the scale factor forced to 1 beyond the
  window's evidenced cutoff and the reference profile smoothed per point 4.
  This is LocScale without servalcat. Built only if slice 1 leaves a
  measurable spatial residual in the fixture's recovered B field.

The ROI's origin in half-map reproducibility means the scaling is not
circular even though the ROI itself uses the atom response of section 4.

## 11. Review findings (2026-09-11)

The design above was revised against the code; these are the facts that fixed
the decisions.

- `build_nu_evidence_state` already computes `q_T` for the band ladder as
  `band_support`; the compact state suffices, and it is Potts-conditioned
  (section 3).
- The evidence state is in-memory only, with no serializer; the ROI must be
  produced where the state exists, hence the dedicated commander (section 2).
- `uncertainty` is the normalized entropy of the same posterior and is not an
  independent cost (section 3).
- `phasecorr_one_atom` is whitened phase correlation against a Nyquist-sharp
  atom and is not a matched filter; `atoms%convolve(…, lp)` and the
  `calpha_finder` NCC machinery are the reusable parts (section 4).
- The atom response is a weaker relative of the detector the ROI serves; the
  NU floor is the gate and the NU-only ablation is mandatory (sections 4, 9).
- A voxel-sum objective drifts with sampling and `eps` already bounds the
  bridge length; physical units and a derived, reported bound (section 5).
- With a 5 A floor on a finer map `q_{T_NU}` does not rank; seeds come from
  the finest evidenced band (section 5).
- The 6VXX / 1JYX benchmark has no half pair; the fixture in section 9 is a
  prerequisite for any of the validation.
- `postprocess_nu` on the PfCRT `_unfil` pair (2026-09-11) produced the best
  PfCRT map to date, which validates post hoc the three things this design
  depends on: the evidence state built from a post-hoc pair, the shell walk
  (`nu_refine=yes`) on such a pair, and the evidence-localized sharpening.
  Run on a regularized pair the same program is unusable; the input contract
  in section 2 follows from that.

## 12. Decisions taken and remaining

Taken: half-pair input; own commander / program; atom response on the
NU-filtered merged map with a kernel bank indexed by the local cutoff; `q_T`
from the compact state with `T_NU` snapped to the band ladder; seeds from the
finest band; up to `K` independently rooted components; physical-unit
objective with the bridge bound derived from `eps`; no uncertainty term in
the objective; absolute B scaling inside the ROI as the model-building map,
global slice first.

Still open:

1. Program name (`nu_roi3D` is a placeholder) and whether the domain object
   lives in `nano/` next to `calpha_finder` or in a new model-building
   directory.
2. Whether user-supplied seed coordinates should be accepted. Trivial to add
   later; not part of the first slice.
3. Whether the ROI should also enter `detect_calpha` candidate scoring as a
   continuous prior, or only gate the search and rank by class. First slice:
   gate and rank.
4. Interpolation of `q_T` between bands in `1 / T` versus snapping. First
   slice: snap.
5. Pseudo-atom element for slice 2 (carbon here, oxygen in LocScale); the
   scattering-factor shapes differ little and the choice only matters if the
   absolute scale is ever used.
