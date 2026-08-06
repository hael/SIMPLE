# NU-Evidence Envelope Masking

> Implementation note. The standalone routine is reachable through
> `simple_exec prg=nu_filt3D nu_envmsk=yes`. The workflow integration described
> in section 4 is implemented in `volassemble` and matching-reference
> preparation after real-volume validation by Cyril and Hans.
>
> **Support decision:** NU filtering now uses only spherical `mskdiam` support.
> `setup_nu_dmats` constructs the sphere internally, so density-derived and
> NU-evidence envelopes cannot enter the normalized Huber objective domain.
> Dilated-envelope support and collar-based null estimation are deferred unless
> the spherical memory cost proves prohibitive in representative runs.
>
> **Reference-masking decision:** the standalone `envref` control has been
> removed. Matching references are solvent-flattened with the current
> NU-evidence envelope whenever `automsk /= no`; there is no separate opt-in.
> In refinement workflows `automsk=yes|tight` is valid only while
> `filt_mode=nonuniform|nonuniform_lpset`; other filtering modes must use
> `automsk=no`.
> The lag-by-one envelope is regenerated freely each cycle and is allowed to
> **shrink** as resolution improves. The former monotonic recovery guard and the
> `pct_signal > 50%` regeneration failsafe have both been removed as too
> conservative.
>
> **`envfsc` is an unfinished branch.** `envfsc=yes` is not currently supported:
> the density automask that used to feed it is no longer generated in the loop,
> and the NU-evidence envelope must never reach the FSC (section 2.1). The route
> needs a simple, fast density-based mask generated on the fly before it can be
> re-enabled. Until then `envfsc=yes` exits with a hard error; it never consumes
> a pre-existing mask or silently degrades to the broad-sphere FSC.
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
artifact serving NU support, `envfsc`, and reference masking. The first role has
now been removed. That separation should be preserved for the NU-evidence
envelope:

| Consumer | Mask source | Rationale |
|---|---|---|
| NU support mask | spherical `mskdiam` | needs solvent inside it to estimate the noise scale and evidence null |
| Matching reference before reprojection | **NU-evidence envelope**, applied whenever `automsk /= no` | highest-value use; belt removal is pure gain for alignment |
| FSC solvent correction (`envfsc`) | **unfinished branch; hard error** until a fast density mask can be generated on the fly | resolution-derived is circular; NU envelope must not enter FSC |
| Derived map for display or deposition support | NU-evidence envelope | preserve unmasked base/half maps as primary artifacts |

Reference masking is no longer gated by a separate `envref` control; that
command-line and library variable has been removed. Whenever automasking is
enabled (`automsk = yes|tight`), the matcher solvent-flattens each reference with
the lagged NU-evidence envelope, falling back to the sphere only until the first
envelope exists.

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

The density `automask3D` routine remains available for standalone and legacy
consumers, but it is no longer generated in the refinement loop. `envfsc=yes` is
therefore an **unfinished branch** that exits with a hard error. Re-enabling it
properly requires a small, fast density-based masker selected on local mean and
generated on the fly from the current map rather than loaded from persistent
state. The two masks answer different questions.

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

### 2.3 Why the matching reference is the right consumer, and its one trap

Solvent-flattening the reference is where a good envelope pays off most: the
detergent belt is noise for alignment, and removing it improves the signal in
every reprojection.

The trap is that reference masking makes the mask self-fulfilling. Anything
zeroed in the reference produces no signal in the reprojections, so particles
never align that region and it never improves — whether or not it was real. For
the belt that is the desired outcome. For a flexible domain that fell below
threshold on one iteration it is a trap that could close permanently.

The mitigation is not a shrink guard but the opposite: the envelope is
regenerated from scratch each cycle and is **free to shrink or grow** as the
evidence changes. Spherical NU support keeps every omitted region observable to
the evidence calculation, so a domain that recovers reproducible signal at higher
resolution re-enters the envelope on the next regeneration. A monotonic
grow-only guard was tried and rejected: it froze early, possibly loose, envelopes
and defeated the resolution-driven tightening this feature exists to provide.

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
With spherical NU support the envelope is regenerated from scratch each cycle and
may shrink freely, so no temporal guard is needed on the reference-masking path.

### 3.3 `nu_refine` interaction (done)

Raw evidence is accumulated over the static bank inside `setup_nu_dmats`.
`extend_nu_filter_highres_shell_next` snapshots each challenger's raw,
unsmoothed objective and commits it to the evidence minimum only after the shell
is accepted. The refinement envelope therefore describes the completed accepted
bank.

## 4. Workflow integration

### 4.1 Artifact and lifecycle

Introduce a second per-state artifact rather than overloading the existing one:

- `automask3D_stateNN.mrc` — legacy density-derived artifact, no longer
  generated or consumed by the FSC path in the refinement loop
- `nu_envmask3D_stateNN.mrc` — new, NU-evidence envelope, consumed by
  reference masking (whenever `automsk /= no`) and final-map masking

Lifecycle is **lag-by-one**: iteration *N*'s evidence writes the envelope used by
iteration *N+1*, mirroring how `nu_highres_depth_stateNN.txt` already persists
across iterations. This needs no reordering of the assembly sequence and no
second NU pass. The envelope must be derived before `nu_filter_vols`, which
calls `release_nu_filter_unary_storage` and frees `dmats_mask`.

The envelope is regenerated on the ordinary automask cadence and simply
overwrites the per-state file each time. There is no monotonic recovery guard and
no `pct_signal > 50%` failsafe: the envelope is allowed to shrink as resolution
improves. The over-occupancy condition is still reported in the stats block as a
diagnostic, but it no longer suppresses regeneration.

### 4.2 Consumer plumbing

`prepare_matching_reference_mask` in
`src/main/strategies/search/simple_matcher_refvol_utils.f90:349` reads a
per-state mask, checks it against `box_crop`/`smpd_crop`, and falls back to a
spherical mask when missing or incompatible. It runs whenever `automsk /= no`
(which requires an active NU filtering mode; there is no separate `envref`
control) and follows the fallback chain
`nu_envmask3D -> automask3D -> sphere`, which also provides the bootstrap at
`startit`.

### 4.3 Hard constraints

- the NU-evidence envelope must never reach `envfsc` (section 2.1);
- the NU support must remain the spherical `mskdiam` support (section 2.2);
- the envelope must be free to shrink as resolution improves; no monotonic
  grow-only guard (section 2.3).

## 5. Implementation phases

- [x] Phase 0a: make spherical `mskdiam` support an invariant of the NU setup API.
- [x] Phase 0b: align public policy and repository skills with spherical-only support.
- [x] Phase 0c: validate the standalone routine across representative volumes
      and fix the public defaults from the parameter study (section 6).
- [x] Phase 1: allow free shrink/grow regeneration; no recovery guard, no
      `pct_signal` failsafe (section 2.3, 4.1).
- [ ] Deferred: collar-based null estimation (3.2), only if a future memory
      optimization reintroduces dilated support.
- [x] Phase 2: feed accepted `nu_refine` extension shells into the raw evidence
      baseline (3.3).
- [x] Phase 3: `volassemble` writes `nu_envmask3D_stateNN.mrc` under lag-by-one;
      `simple_vol_pproc_policy` gains the second artifact in its plan type.
- [x] Phase 4: reference masking (any `automsk /= no`) reads the new artifact
      via `prepare_matching_reference_mask` with the three-step fallback chain;
      the `envref` control is removed.
- [x] Phase 5: after workflow validation, promote the NU-evidence artifact and
      consumer lifecycle to `doc/policies/automasking_policy.md` and
      `doc/policies/nonuniform_filtering_policy.md`.
- [ ] Deferred: `envfsc=yes` needs a simple, fast density-based masker generated
      on the fly before it can be re-enabled (section 2.1); until then it throws.

## 6. Parameter optimization

### 6.1 Public parameter set

| Parameter | Default | Role |
|---|---|---|
| `nu_msk_sig` | 3.0 | threshold, in null MADs above the null median |
| `amsklp` | 8 A | physical envelope scale; margin smoothing radius is `min(1.5 x amsklp, 30 A)` |

These are the only two public envelope-shape controls. `nu_envmsk` remains the
feature toggle and `mskdiam` remains the mandatory spherical support diameter;
neither is an envelope tuning constant.

The secondary choices are fixed production policy:

| Constant | Value | Role |
|---|---:|---|
| scale-free evidence | no | A baseline-to-best ratio can prevent weak but well-ordered density from being outvoted by a high-contrast core; production retains the absolute Huber-cost margin |
| density weight | 0.0 | A positive local-density weight retains strong but poorly ordered density; zero keeps segmentation evidence-only |
| MRF beta | 1.0 | Binary MRF boundary smoothness; higher values produce smoother boundaries |
| minimum component fraction | 0.1 | Smallest connected component kept, expressed as a fraction of the largest |
| binary growth | 1 A | slightly expand the accepted support |
| cosine edge | 6 A | soften the final mask boundary |

The last two values translate the previous 1-pixel and 6-pixel defaults under
the requested 1 A/pixel reference sampling. At run time each is converted to
the nearest voxel count from the actual input sampling, with a minimum of one
voxel. This changes units and sampling behavior, not the masking sequence.

### 6.2 Structure of the two-dimensional space

`nu_msk_sig` and `amsklp` are the dominant pair and are **coupled**: larger
`amsklp` pools evidence over a wider support, which raises the margin in weak
regions and is therefore equivalent to loosening the threshold. They trade off
along a ridge and must be swept jointly.

All secondary values are intentionally frozen. This preserves the validated
algorithm and makes `(nu_msk_sig, amsklp)` the complete optimization space.
The internal NU API retains alternate evidence and segmentation options for
diagnostic tests, but they are not standalone `nu_filt3D` controls.

### 6.3 Sweep efficiency

The expensive step is `setup_nu_dmats`. `amsklp` changes the evidence field
(`calc_nu_evidence_margin`), whereas `nu_msk_sig` changes only the segmentation
threshold (`calc_nu_evidence_score`). The natural sweep is therefore an outer
loop over `amsklp` and an inner loop over `nu_msk_sig`, reusing NU setup and
cutoff optimization. A dedicated in-process driver can exploit this; a shell
loop over `nu_filt3D` repeats setup at every point.

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
- Keep the production absolute-margin fixture as the required regression. Any
  internal scale-free diagnostic belongs in a separate, explicitly
  low-occupancy fixture and must not change the standalone interface.

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
after N `refine3D` iterations with reference masking active (`automsk /= no`).
Use only as a **final confirmation** on the one or two settings that survive
tiers 1-3. Never as the sweep objective.

### 6.5 Do not tune on the reported FSC resolution

Reference masking plus a resolution-derived mask is precisely the loop that
inflates reported resolution (section 2.1). Optimising parameters against it
selects for the mask that most aggressively retains only high-resolution voxels —
the degenerate solution, which is a tight shell around the best-ordered core. It
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
2. **Flexible multi-domain complex** — the cut-off-domains case and strongest
   test of whether the two public controls transfer across local resolution.
3. **Small rigid globular protein** — baseline sanity; should be easy.
4. **Map with an internal cavity or channel** — exercises 3D hole filling.
5. **Low-occupancy or partially-occupied subunit** — tests whether the absolute
   margin remains useful at weak occupancy.
6. **One specimen at several particle counts** (for example full, 1/4, 1/16) —
   the SNR-dependence check underpinning section 6.7.

## 8. What to record per run

The routine already logs null median, null MAD, threshold, envelope scale, ICM
beta and iteration count, seed and signal voxel counts and percentages, and
component counts kept versus found. Emit these plus the section 6.4 metrics as
one CSV row per sweep point so surfaces can be plotted directly.

Watch the signal percentage specifically: above 50% the median/MAD null is not
trustworthy and the routine says so. Spherical geometry does not guarantee this
condition: if it triggers, `mskdiam` is too tight or the whole-support null model
is unsuitable. This is now reported as a diagnostic only; it no longer blocks
regeneration, so a persistently high signal percentage is a signal to revisit
`mskdiam` or the null model rather than an automatic stop.

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

Workflow integration now provides `nu_envmask3D_stateNN.mrc`, lag-by-one
reference-mask consumption (whenever `automsk /= no`; the `envref` control has
been removed), free shrink/grow regeneration with no recovery guard or
`pct_signal` failsafe, and accepted-shell evidence updates. `envfsc=yes` remains
an unfinished hard-error branch pending a fast density-based mask generated on
the fly. A purpose-built
low-occupancy internal scale-free-evidence regression remains future work.
