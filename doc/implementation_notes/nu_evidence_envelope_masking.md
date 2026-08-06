# NU-Evidence Envelope Masking

> Implementation note. The standalone routine exists and is reachable through
> `simple_exec prg=nu_filt3D nu_envmsk=yes`. Workflow integration described
> in section 4 is **not implemented**. Validate the routine first (sections 6-8),
> then implement.
>
> **Support decision:** NU filtering now uses only spherical `mskdiam` support.
> `setup_nu_dmats` constructs the sphere internally, so density-derived and
> NU-evidence envelopes cannot enter the normalized Huber objective domain.
> Dilated-envelope support and collar-based null estimation are deferred unless
> the spherical memory cost proves prohibitive in representative runs.
>
> **Validation snapshot:** the focused NU test executables build and link. On
> the current synthetic fixture, spherical support contains 44,473 voxels and
> the default absolute margin gives recall 1.000 with solvent false-positive
> rate 0.054. The scale-free cost-improvement ratio gives recall 0.952 with
> solvent false-positive rate 0.008 on the same fixture.

The implemented standalone algorithm is described in
[`NU-Evidence Envelope Mask in nu_filt3D`](../algorithms/nu_evidence_envelope_mask.md).

## 1. Summary

The nonuniform filter's unary objective is a per-voxel cross-half prediction
error. In solvent the half maps are uncorrelated at every bandwidth, so no
candidate should beat the coarsest one systematically; finite noise and taking
the minimum over several candidates still produce positive voxelwise margins.
Inside ordered density the objective can have a real minimum at the local SNR
crossover. The per-voxel improvement over the coarsest baseline is therefore an
empirical measure of local orderedness, and it is what this routine segments on:

```
margin(v) = cost(v, coarsest) - min over candidates of cost(v, c)
```

Two properties observed on a detergent-solubilised membrane protein:

- the detergent belt is **excluded cleanly**, because it is real density with no
  cross-half consistency at any bandwidth;
- with symmetric smoothing (section 3) the boundary tracks the protein envelope
  closely without any density term.

That combination is what makes the mask interesting: it discriminates *ordered*
density from *disordered* density, which an amplitude threshold cannot do.

The margin is a best-of-bank statistic. Its solvent distribution depends on
candidate count, candidate correlation, objective smoothing, and any accepted
high-resolution extension. The null must therefore be validated for the exact
bank used to generate an envelope; it is not a universal zero-centered law.

It is also what makes it dangerous in the wrong consumer. Section 2 is the
central policy claim of this note.

## 2. Mask ownership: one envelope per consumer, not one artifact for all

Before the spherical-support refactor, `automask3D_stateNN.mrc` was a single
artifact serving NU support, `envfsc`, and `envref`. The first role has now been
removed. That separation should be preserved for the NU-evidence envelope:

| Consumer | Mask source | Rationale |
|---|---|---|
| NU support mask | spherical `mskdiam` | needs solvent inside it to estimate the noise scale and evidence null |
| Matching reference before reprojection (`envref`) | **NU-evidence envelope** | highest-value use; belt removal is pure gain for alignment |
| FSC solvent correction (`envfsc`) | density-derived (existing ICM/Otsu masker), or none | resolution-derived is circular |
| Derived map for display or deposition support | NU-evidence envelope | preserve unmasked base/half maps as primary artifacts |

### 2.1 Why not the FSC path

Any tight mask inflates FSC by removing uncorrelated solvent, and phase
randomization exists to correct for exactly that. The correction assumes the
mask was chosen **independently of the half-map correlation structure**.

A NU-evidence mask violates that assumption by construction: it is selected to
contain the voxels where the half maps agree at high frequency. Randomizing the
phases destroys the very correlation the selection was based on, so the
randomized-masked curve cannot reproduce the selection bias and the correction
systematically under-corrects. The result is resolution inflation that appears
to be properly solvent-corrected.

A density-derived mask is selected on the local mean rather than on cross-half
agreement. Density and resolution do correlate, so it is not perfectly
independent, but it is far less circular and it is the regime phase
randomization was designed for.

**The existing `automask3D` density masker must therefore be kept, not
retired.** The two masks answer different questions.

### 2.2 Why not the NU support mask

Two reasons, one statistical and one practical.

Statistical: the null is estimated from the support itself. If the support is
the tight envelope, there is no solvent left to estimate it from, and the
envelope erodes monotonically across iterations while looking plausible at every
step.

Practical: `dmats_mask` is `(n_nu_mask, n_candidates)` and support size is the
dominant memory term. Using the envelope as NU support was a deliberate memory
decision in the former policy. A full `mskdiam` sphere is not free — on a 300^3
box an envelope might be ~3M voxels against ~12M for the sphere, roughly 300 MB
against 1.2 GB at the 24-candidate cap, plus full-grid temporaries.

The initial implementation accepts this memory cost and uses only spherical
support. This is the conservative choice for the current normalized Huber
objective and whole-support evidence-null estimator. A dilated previous
envelope remains a possible future memory optimization, but only together with
the null change in section 3.2 and explicit shrinkage guards.

### 2.3 Why `envref` is the right consumer, and its one trap

Solvent-flattening the reference is where a good envelope pays off most: the
detergent belt is noise for alignment, and removing it improves the signal in
every reprojection.

The trap is that `envref` makes the mask self-fulfilling. Anything zeroed in the
reference produces no signal in the reprojections, so particles never align that
region and it never improves — whether or not it was real. For the belt that is
the desired outcome. For a flexible domain that fell below threshold on one
iteration it is a trap that closes permanently.

Periodic regeneration only helps if the mask can *recover* ground. Spherical NU
support makes the omitted region observable to the evidence calculation, but a
temporal shrink/recovery guard is still required before enabling `envref`.

## 3. Required changes to the routine before integration

### 3.1 Symmetric evidence smoothing (done)

`dmats_mask` is smoothed at candidate-dependent radii, `min(1.5 x LP, 30 A)`:
30 A for the 20 A baseline against 6 A for the 4 A member. That is correct for
label selection and wrong for locating a boundary — it blurred the baseline five
times harder than its competitors and eroded the envelope inward on a 30 A
scale, removing any domain thinner than ~60 A entirely.

Fixed by caching the raw coarsest cost and running per-voxel minimum in
`setup_nu_dmats` before `smooth_nu_objective`, then smoothing the *difference*
once at a scale set by `amsklp`. Smoothing the difference is equivalent to
smoothing both terms with the same kernel. The filter's own label selection is
untouched.

### 3.2 Collar-based null estimation (deferred with spherical-only support)

The current null is the median and MAD over the whole support, which assumes
solvent is the **majority**. True for a generous sphere, false for a dilated
envelope where the molecule may be 60-80% of the voxels. The routine warns above
50% but cannot correct for it.

Proposed: estimate the null from the **collar** — voxels inside the dilated
support but outside the previous envelope. That region is solvent by
construction, since it is precisely what the previous iteration rejected, and it
removes the majority assumption entirely.

Failure mode to guard: if one iteration's envelope comes out too tight, the
collar contains real density, the null inflates, the threshold rises, and the
next envelope is tighter still. Guards:

- generous dilation (15-20 A), so the collar stays mostly true solvent even when
  the envelope is somewhat wrong;
- a low quantile within the collar rather than the median;
- a floor on how much the envelope may shrink per iteration.

This work is not required while NU support remains exclusively spherical. If
dilated support is reconsidered, bootstrap from the sphere at `startit`; do not
use a tight density automask as a substitute for a solvent-containing collar.
The temporal shrink/recovery guard remains required before `envref` integration
even with spherical NU support, because reference masking itself can make a
false-negative region self-fulfilling.

### 3.3 `nu_refine` interaction (not done)

Raw evidence is accumulated over the static bank only, inside the candidate loop
in `setup_nu_dmats`. `extend_nu_filter_highres_shell_next` appends candidates
afterwards and does not update it. Harmless for `nu_filt3D`, which never
extends, but `refine3D_auto` enables extension and the accepted shells must feed
the evidence field before envelope generation is integrated into `volassemble`.

## 4. Workflow integration (deferred)

### 4.1 Artifact and lifecycle

Introduce a second per-state artifact rather than overloading the existing one:

- `automask3D_stateNN.mrc` — unchanged, density-derived, consumed by `envfsc`
  and retained as the fallback for `envref`
- `nu_envmask3D_stateNN.mrc` — new, NU-evidence envelope, consumed by `envref`
  and final-map masking

Lifecycle is **lag-by-one**: iteration *N*'s evidence writes the envelope used by
iteration *N+1*, mirroring how `nu_highres_depth_stateNN.txt` already persists
across iterations. This needs no reordering of the assembly sequence and no
second NU pass. The envelope must be derived before `nu_filter_vols`, which
calls `release_nu_filter_unary_storage` and frees `dmats_mask`.

### 4.2 Consumer plumbing

`prepare_matching_reference_mask` in
`src/main/strategies/search/simple_matcher_refvol_utils.f90:349` already reads a
per-state mask, checks it against `box_crop`/`smpd_crop`, and falls back to a
spherical mask when missing or incompatible. That fallback is exactly the
bootstrap behaviour lag-by-one needs at `startit`. Point it at the new filename
with a fallback chain of `nu_envmask3D -> automask3D -> sphere`.

### 4.3 Hard constraints

- the NU-evidence envelope must never reach `envfsc` (section 2.1);
- the NU support must remain the spherical `mskdiam` support (section 2.2);
- envelope regeneration must not be able to shrink without bound (section 3.2).

## 5. Implementation phases

- [x] Phase 0a: make spherical `mskdiam` support an invariant of the NU setup API.
- [x] Phase 0b: align public policy and repository skills with spherical-only support.
- [ ] Phase 0c: validate the standalone routine across the test set (section 7)
      and fix defaults from the parameter study (section 6).
- [ ] Phase 1: add the temporal shrink/recovery guard required by `envref`.
- [ ] Deferred: collar-based null estimation (3.2), only if a future memory
      optimization reintroduces dilated support.
- [ ] Phase 2: feed accepted `nu_refine` extension shells into the raw evidence
      baseline (3.3).
- [ ] Phase 3: `volassemble` writes `nu_envmask3D_stateNN.mrc` under lag-by-one;
      `simple_vol_pproc_policy` gains the second artifact in its plan type.
- [ ] Phase 4: repoint `prepare_matching_reference_mask` at the new artifact
      with the three-step fallback chain.
- [ ] Phase 5: after workflow validation, promote the NU-evidence artifact and
      consumer lifecycle to `doc/policies/automasking_policy.md` and
      `doc/policies/nonuniform_filtering_policy.md`.

## 6. Parameter optimization

### 6.1 The parameter set

| Parameter | Default | Role |
|---|---|---|
| `nu_msk_sig` | 3.0 | threshold, in null MADs above the null median |
| `amsklp` | 8.0 | envelope scale; margin smoothing radius is `min(1.5 x amsklp, 30 A)` |
| `nu_msk_beta` | 1.0 | binary MRF boundary smoothness |
| `nu_msk_dens` | 0.0 | weight of the local density term |
| `nu_msk_rel` | no | scale-free cost-improvement ratio |
| `nu_msk_minvol` | 0.1 | smallest connected component kept, as a fraction of the largest |
| `binwidth` | 1 | dilation layers |
| `edge` | 6 | cosine edge width |

### 6.2 Structure of the space — do not grid all eight

`nu_msk_sig` and `amsklp` are the dominant pair and are **coupled**: larger
`amsklp` pools evidence over a wider support, which raises the margin in weak
regions and is therefore equivalent to loosening the threshold. They trade off
along a ridge and must be swept jointly.

Everything else is secondary:

- `nu_msk_beta` mostly controls speckle and boundary regularity once the evidence
  is right. Weakly coupled to `nu_msk_sig`, since higher beta can partially
  rescue a slightly-too-high threshold by filling in. Sweep third, coarsely.
- `nu_msk_dens` is a **rescue lever, not a default**. It reintroduces the
  detergent belt in proportion to its weight, since density cannot distinguish
  ordered protein from disordered lipid. Engage only if flexible density is
  still being lost after `nu_msk_sig`/`amsklp` are settled. Observed: 1.0 is far
  too high; expect the useful range below 0.4 if it is needed at all.
- `nu_msk_rel` is experimental. Its original bounded fractional reduction could
  produce an empty envelope whenever `median + nu_msk_sig * MAD` exceeded one.
  It now uses the unbounded ratio `baseline / best - 1`, with both costs
  protected by the denominator floor. Keep it off by default until its intended
  low-occupancy use has a dedicated fixture and real-data validation.
- `nu_msk_minvol`, `binwidth`, `edge` are the topology and morphology tail. They
  determine the mask's *finish*, not its *quality*. Freeze them at sensible
  values before sweeping anything else, chosen from the consumer's needs —
  `envref` wants a soft, slightly generous mask.

**Recommended procedure**: freeze the tail, sweep `(nu_msk_sig, amsklp)` as a 2D
grid, then `nu_msk_beta`, then rescue levers only if needed.

### 6.3 Sweep efficiency

The expensive step is `setup_nu_dmats`. Structure the sweep to exploit what each
parameter actually invalidates:

- `amsklp` and `nu_msk_rel` change the **evidence field** (`calc_nu_evidence_margin`);
- `nu_msk_sig`, `nu_msk_beta`, `nu_msk_dens` change only the **segmentation**
  (`calc_nu_evidence_score`, `segment_nu_evidence`);
- `nu_msk_minvol`, `binwidth`, `edge` change only the **tail**
  (`envmask3D_from_lmask`).

So the natural sweep is an outer loop over `amsklp` and an inner loop over
everything else, with `setup_nu_dmats` and `optimize_nu_cutoff_finds` run once
per outer iteration. A dedicated sweep driver in `production/tests/` doing this
in-process is worth the small effort; a shell loop over `nu_filt3D` re-pays the
NU setup on every point.

### 6.4 Objective functions

Ranked by how directly they measure what we want, against how available they are.

**Tier 1 — synthetic, ground truth known.** Extend
`production/tests/simple_test_nu_envmask.f90`. Report Dice/IoU, recall, and
solvent false-positive rate, and add:

- **mean surface distance** to the true envelope. More sensitive than Dice to the
  erosion failure mode, which is what bit us — Dice degrades gracefully under
  uniform erosion while surface distance reports it directly.
- **thin-domain retention**: build a fixture with a deliberately thin
  (30-50 A) lower-resolution appendage and score the fraction retained. This is
  the specific failure the symmetric-smoothing fix targets and it must have a
  regression pin.
- Keep the default absolute-margin fixture as the required regression. Test
  `nu_msk_rel` separately on an explicitly low-occupancy fixture; the present
  general fixture must not assert that both statistics produce the same mask.

**Tier 2 — real data with an atomic model.**

- *Model coverage*: fraction of model atoms falling inside the mask. Target 1.0.
- *Volume ratio*: mask volume divided by expected molecular volume. For protein,
  `V(A^3) ~ 1.21 x MW(Da)`. Masks legitimately exceed this because of hydration
  and the soft edge, so the target is a controlled factor (roughly 1.3-2.0),
  not 1.0.
- Combined scalar: penalise coverage below 1 and volume ratio above target.
  These two together are a strong and cheap precision/recall pair.

**Tier 3 — real data, no model. This is the general case.**

- **Half-set reproducibility.** Split the particles into two disjoint halves,
  reconstruct each into its own gold-standard pair, derive a mask from each
  independently, and measure Dice between the two masks. A mask driven by real
  structure is reproducible; one driven by noise is not. This requires no ground
  truth and directly measures whether the evidence is real rather than fitted.
  **Recommended as the primary objective for real data.**
- *Solvent flatness*: RMS of the map outside the mask compared against RMS in a
  far-corner region. Elevated RMS outside means the mask is too tight.
- *Annotated-feature exclusion*: for a membrane protein, annotate the belt once
  and score exclusion on every sweep point thereafter.

**Tier 4 — downstream, expensive and confounded.** Resolution and map quality
after N `refine3D` iterations with `envref=yes`. Use only as a **final
confirmation** on the one or two settings that survive tiers 1-3. Never as the
sweep objective.

### 6.5 Do not tune on the reported FSC resolution

`envref` plus a resolution-derived mask is precisely the loop that inflates
reported resolution (section 2.1). Optimising parameters against it selects for
the mask that most aggressively retains only high-resolution voxels — the
degenerate solution, which is a tight shell around the best-ordered core. It
will look like a large improvement.

Tier 3 half-set reproducibility is the honest substitute: it rewards masks that
are *real*, not masks that are *tight*.

### 6.6 Prefer a plateau over a peak

Report the objective as a surface over `(nu_msk_sig, amsklp)` and choose the
**centroid of the region within a few percent of the best**, not the argmax. A
setting that is optimal but on a knife edge is worse than a slightly suboptimal
one in a broad basin, because real datasets vary and a knife-edge optimum will
not transfer.

Final defaults should be the setting that is *acceptable on every* test specimen,
not optimal on any one. Record per-dataset optima alongside the consensus so the
spread is visible.

### 6.7 A likely outcome worth planning for

The margin scales with SNR, so a threshold tuned at 3 A with 500k particles may
not hold at 6 A with 20k. If the sweep shows `nu_msk_sig` is strongly dataset- or
SNR-dependent, that is evidence the fixed MAD multiple is the wrong
parameterisation, and the answer is to make the threshold adaptive — Otsu on the
margin histogram rather than a fixed multiple of the null scale. That would
remove the parameter entirely and would mirror what `automsk=tight` already does
for density. Treat this as a live design fork, not a fallback.

## 7. Test set

The set must span the failure modes, not just the easy case:

1. **Detergent-solubilised membrane protein** — the belt case. Already
   validated qualitatively; needs a quantitative baseline.
2. **Flexible multi-domain complex** — the cut-off-domains case, and the
   strongest test of whether `nu_msk_dens` is needed.
3. **Small rigid globular protein** — baseline sanity; should be easy.
4. **Map with an internal cavity or channel** — exercises 3D hole filling.
5. **Low-occupancy or partially-occupied subunit** — the case `nu_msk_rel` was
   built for.
6. **One specimen at several particle counts** (for example full, 1/4, 1/16) —
   the SNR-dependence check underpinning section 6.7.

## 8. What to record per run

The routine already logs null median, null MAD, threshold, envelope scale,
scale-free flag, ICM beta and iteration count, seed and signal voxel counts and
percentages, and component counts kept versus found. Emit these plus the section
6.4 metrics as one CSV row per sweep point so surfaces can be plotted directly.

Watch the signal percentage specifically: above 50% the median/MAD null is not
trustworthy and the routine says so. Spherical geometry does not guarantee this
condition: if it triggers, `mskdiam` is too tight or the whole-support null model
is unsuitable. Automatic workflow use of that envelope must stop rather than
continuing after a warning.

## 9. Current state

Implemented and reachable through `nu_filt3D`:

- `src/main/nu_filt/simple_nu_filter_envmask.f90` — evidence margin, null
  estimation, binary MRF segmentation, diagnostics
- `src/main/image/simple_image_msk.f90` — `envmask3D_from_lmask`, the topology
  and morphology tail
- `src/main/image/simple_image_bin.f90` — 3D hole filling unlocked, corner-seed
  guard added, `find_ccs` corrected to true 26-connectivity
- `production/tests/simple_test_nu_envmask.f90` — synthetic absolute-margin and
  scale-free cost-improvement-ratio regressions
- `production/tests/simple_test_cc_connectivity.f90` — connectivity regression

Outputs `<vol>_nu_evidence.mrc` (raw margin field) and `<vol>_nu_envmask.mrc`
(soft mask). The evidence map is the primary diagnostic: if it does not separate
solvent from density, no threshold or smoothness setting will rescue the mask.

Not yet implemented:

- workflow-owned `nu_envmask3D_stateNN.mrc` artifacts
- lag-by-one `envref` consumption and temporal recovery guards
- evidence updates from accepted `nu_refine` extension shells
- a purpose-built low-occupancy `nu_msk_rel` regression
