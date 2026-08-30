# NU-evidence local sharpening: model-free LocScale from cross-half evidence

## Status

Proposal (2026-08-27), **SCHEDULED as PRESSING (2026-08-29, user
direction)** — item 1 on the active dev list in `pcg_priors.md`
(Stage 6.6 run records). The original precondition is met: the direct
NU-evidence prior cleared its Gate C/D program and the Stage 6.6
nu_refine evidence-bank extension validated on 1WCM. The motivation is
now empirical, not speculative: on real data (PfCRT and others) a single
isotropic B-factor does not produce acceptable postprocessed maps — only
bgal- and streptavidin-like specimens tolerate it. This remains a
postprocessing experiment that consumes Stage 6 infrastructure, not a
competitor to it. No solver, base-solve, replay, or artifact behavior
changes are proposed here. Note the Wilson-target variant 2.3(c) is
DEAD: the Wilson prior was adjudicated against and removed from the
codebase (2026-08-29, `pcg_priors.md` Stage 7 record); variants 2.3(a)
evidence-derived and 2.3(b) local-B remain the candidates.

## 1. The idea

LocScale sharpens a map by rescaling *local* Fourier amplitudes against a
reference while preserving the observed phases. The pseudo-atomic (or hybrid)
reference plays two roles: a local confidence estimate (where is the map
good, and to what resolution?) and a target spectrum (what should the
amplitudes look like?).

The NU replay evidence state already provides the first role without any
model: `a_b(v)` is a calibrated, graded, per-band local support confidence
derived from noise-whitened cross-half prediction with an explicit calibrated
null (`nu_evidence_state`, built by `build_nu_evidence_state` from the
unregularized base half pair). Compared with a pseudo-atomic reference it
has no model bias, no map-model circularity, and no reference-refinement
cycle; compared with FDR/local-resolution confidence maps it carries an
explicit noise model and a validated null competitor. The proposal is
therefore: **LocScale-style local amplitude scaling in which both the
confidence field and (in the first variants) the target spectrum are derived
from the NU evidence — fully model-free.**

## 2. Construction

### 2.1 The band frame is already built

The PCG reconstructor's `Q_NU` machinery is an analysis/synthesis frame:
disjoint radial band operators `B_b` on the padded lattice with per-voxel
spatial fields applied between analysis and synthesis
(`apply_nu_precision`). The replay applies penalties `W_b = [p(1-a_b)]^2`;
the sharpener applies **gains** on the same stack:

```text
x_sharp = sum_b g_b(v) * (B_b x),   g_b(v) >= 0
```

with the DC/mean component passed through unchanged. Phases are untouched by
construction (real nonnegative gains on band-limited real-space components).
The frozen compact evidence state is shared read-only, exactly as
`pcg_priors.md` §6/§11 mandates for diagnostics — one evidence identity, no
second NU analysis that can disagree with the replay.

### 2.2 Two gain layers

1. **Local Wiener layer (suppression).** `g_b(v)` proportional to the local
   band SSNR/(1+SSNR) surrogate derived from the same cross-half costs that
   produce `a_b`. This upgrades the production NU filter from *binary local
   cutoff selection* to *graded local amplitude weighting*, and it supplies
   the principled shipped-map rolloff for the NU replay's beyond-band
   retention (the Gate C watch item: `Q_NU` suppresses by lack of evidence,
   not by SSNR, so the raw replay map carries unshrunk near/beyond-band
   amplitude that postprocessing owns).
2. **Restoration layer (the LocScale analogue).** Boost damped-but-supported
   bands toward a target spectrum, bounded by evidence: amplification only
   within bands the evidence supports, tapered by confidence and by the
   entropy/uncertainty field at ambiguous boundaries. This is the layer that
   distinguishes the proposal from the existing NU filter, which can only
   attenuate.

### 2.3 Target-spectrum options (model-free ladder)

In increasing ambition; (a) is the first implementation:

- **(a) Self-referential:** match each locus's band spectrum to the mean
  spectrum of the highest-confidence regions of the same map. Fully
  data-driven, no assumptions beyond "the best parts of this map are what
  this molecule's spectrum looks like."
- **(b) Local B-factor fit:** fit a two-parameter local decay
  (scale, B) to the bandwise cross-half agreement per voxel and invert it.
  Compact and robust to the coarseness of the 4-band frame; effectively a
  model-free local B-factor map with calibrated confidence.
- **(c) Wilson expected spectrum:** the `pcg_priors.md` §5.5 Wilson object,
  used here in its gentlest possible role — a sharpening *target* rather
  than a prior. Requires only composition-level assumptions, no atomic
  coordinates. This variant is the natural bridge between the sharpening
  experiment and the Wilson prior stage, and must wait for that stage's
  spectrum-source mechanism.

### 2.4 Band granularity

Four bands (20/12/8/5 A) are coarse for sharpening; LocScale operates in
fine radial shells. Mitigations in order of preference: the local-B fit of
2.3(b), which interpolates smoothly across band boundaries; or a finer gain
bank from the full 8-label candidate bank plus accepted `nu_refine`
extensions — affordable post-hoc because the cost is paid once per map, not
once per CG iteration. Record the frame normalization when refining bands
(the `Q_NU` lesson: the pad/crop band frame is NOT tight, so band-partition
changes change the operator and must be measured, not assumed invariant).

## 3. Discipline and boundaries

- **Half-map independence.** Applied identically to both halves (or to the
  merged map), the sharpener inflates shipped-pair FSC exactly as shared
  regularization does. Same rule as the replay: the unregularized base pair
  keeps sole resolution authority; the sharpened map is a
  display/interpretation product. Never feed the sharpened map or its
  evidence into FSC solvent correction or resolution claims.
- **Not a solver component.** The LocScale risk rows in `pcg_priors.md` §9
  (amplitude target nonlinear in `x`; common targets carrying phases)
  concern in-solve use. As postprocessing none of them apply: the operator
  acts once on a finished estimate, phases are preserved, and no CG
  assumption is involved.
- **Amplification is evidence-bounded by design.** Where a global-B
  sharpener amplifies noise indiscriminately, gains here are capped by band
  support and the calibrated null: no evidence, no amplification. This is
  the feature that justifies the experiment; any implementation that adds an
  uncapped user gain knob has left the design.
- **Ownership.** Evidence: `src/main/nu_filt/` (the frozen compact state and
  its expansion, shared with the replay). Gain synthesis and application:
  volume-domain postprocessing beside the existing NU filter
  (assembly-owned, per the nonuniform filtering policy — matchers and
  search never sharpen). Products are derived (`_nu_sharp` naming beside
  `_nu_filt`/`_nu_locres`), never replacements for base reconstructions.

## 3b. Implementation v1 (2026-08-29, user-directed design decisions)

User-directed: the path is ISOLATED from the standard postprocess commander
(global B-factor + FSC filter, unchanged) behind a dedicated `postprocess_nu`
commander. What was built:

- **Commander** `postprocess_nu`
  (`src/main/commanders/simple/simple_commanders_postprocess_nu.f90`, routed
  in `simple_exec_filter`, UI in `simple_ui_filter` beside `nu_filt3D`).
  Standalone file interface, no sp_project: `vol1` (odd) / `vol2` (even)
  UNREGULARIZED half maps, `smpd`, `mskdiam`, optional `outvol`, `nthr`,
  `nu_refine` (default **yes** here — finer evidence granularity is
  affordable post-hoc, §2.4; the cost is paid once per map, not per CG
  iteration). Evidence lifecycle mirrors the Q_NU replay exactly:
  `setup_nu_dmats -> optimize_nu_cutoff_finds -> [accepted shell walk] ->
  build_nu_evidence_state -> assert_nu_evidence_replay_ready` — one evidence
  identity, no second NU analysis. Outputs: `_nu_sharp` even/odd/merged.
- **Gain synthesis** (`simple_nu_filter_sharpen.f90`, submodule of
  `simple_nu_filter`): `x_sharp = mean(x) + sum_b g_b(v) * B_b(x - mean(x))`
  with disjoint radial bands from the evidence ladder using the identical
  coarse-to-fine first-match partition as `ensure_nu_band_index`, and
  mean-centering standing in for the Q_NU `C` projector (`apply_filter`
  would otherwise ride DC on every band). Per band:
  `g_b = a_b * (1 + a_b * (min(CAP, max(1, t_b / max(e_b(v), eps))) - 1))`
  where `a_b = 1 - band_w` is the calibrated support confidence (§2.2 Wiener
  layer), `e_b(v)` the local band RMS (squared band component low-passed at
  `NU_SHARP_ENERGY_SMOOTH_FAC=2.0` x the band limit, floored at
  `NU_SHARP_ENERGY_EPS_REL=0.05` x in-support band RMS), and the target
  `t_b` the RMS of `e_b` over voxels with `a_b >= NU_SHARP_REF_CONFIDENCE=0.9`
  (§2.3(a) self-referential); restoration disables for a band with fewer
  than `NU_SHARP_MIN_REF_VOXELS=1000` reference voxels or a target below the
  floor. Restoration cap `NU_SHARP_MAX_GAIN=4.0` (amplitude ratio). Gains
  are computed once from the merged pair and applied identically to both
  halves. No user gain knob exists, per the design mandate; the constants
  are recorded design, and changing them is an experiment.
- **Known v1 limitations, recorded up front:** (i) the gain field steps to
  zero at the spherical evidence-support edge (the sharpened map decays to
  its mean outside `mskdiam`) — acceptable for a display product, revisit
  with a soft support taper if edge artifacts show; (ii) the uncertainty/
  entropy taper of §2.2 is not yet in the gains (confidence-only taper);
  (iii) restoration targets are per-band scalars (variant 2.3(a)), the
  local-B interpolation of 2.3(b) remains the upgrade path if band
  coarseness limits.

### v1 outcome (2026-08-29, PfCRT, FAILED -- recorded, superseded by v2)

First real-data run (PfCRT unfil pair, box 300): the evidence envelope was
genuinely good -- flat-zero solvent tracking the molecular shape, and the
shell walk found a real high-resolution core (10 accepted steps to ~3 A).
But the periphery collapsed to coarse-band blobs and the core was
over-sharpened salt-and-pepper noise. Root cause, confirmed by design
review: the v1 Wiener layer used band support confidence `a_b` as the
shrinkage, but confidence is a calibrated SUPPORT PROBABILITY that
saturates to 1 wherever evidence exists -- it carries no SSNR. So the
finest bands of the raw unfiltered input passed at full noise power, and
the ratio-restoration layer then boosted local amplitude dips toward a
noise-level regional RMS target by up to 4x. Confidence != SSNR; v1
conflated them. (Operational note from the same session: the first run
crashed because `setup_nu_dmats` was called without `evidence_source` --
the provenance contract guard worked as intended -- and a run on
mismatched inputs showed that the mskdiam-too-large fallback warning must
be treated as an input error: check the half-map box/smpd headers.)

## 3c. Implementation v2 (2026-08-29): classical shrink-then-sharpen, localized

What the classical pipeline gets right and v1 skipped: shrinkage and
sharpening are always coupled, in that order (the standard postprocess
applies the FSC-derived per-shell Wiener `fsc2optlp` and only sharpens
inside that rolloff); LocScale's boost is implicitly bounded by the
reference's physical amplitude falloff; and the B-factor is a smooth
two-parameter Guinier fit, not a per-band amplitude ratio. v2 is that
recipe with the evidence deciding WHERE each ingredient acts:

1. **Sharpen:** one classical Guinier B-factor from the merged map
   (`guinier_bfac(HPLIM_GUINIER, finest_evidenced_cutoff)`), applied only
   when the finest evidenced cutoff is finer than
   `NU_SHARP_BFAC_FINEST_A = 5 A` -- the standard postprocess gate.
2. **Filter (the local Wiener surrogate):** per-voxel Butterworth
   low-pass at the frozen state's evidenced local cutoff
   (`selected_cutoff`, from the Potts-smoothed label optimization),
   composed production-NU-filter style from the <=24 distinct per-cutoff
   filtered versions. Sharpening therefore never extends beyond the local
   passband.
3. **Solvent:** null-claimed voxels (cutoff 0) and voxels outside the
   spherical support flatten to the map mean -- the evidence
   envelope behavior the v1 run validated.

Output (user-directed, 2026-08-29): the shipped product is a SINGLE
sharpened merged volume, classical-postprocess style. Every v2 operation
is linear, so sharpening the merged map is identical to averaging two
identically-processed halves; per-half `_nu_sharp` outputs were dropped.
The half pair still supplies the evidence and the resolution authority.

A global `fsc2optlp` is deliberately NOT applied: the global FSC averages
over the map and would erase exactly the core detail the walk validated;
the Butterworth rolloff at the calibrated local cutoff is the local
shrinkage surrogate. The v1 restoration layer (ratio boost toward a
regional target) is retired; if v2's single global B proves too blunt
across regions of very different decay, the recorded upgrade path is the
local-B fit 2.3(b) -- per-region Guinier, still inside the local
passband. The v1 `NU_SHARP_*` band-gain constants were removed from the
code with the design.

## 4. Validation plan

The Stage 6 harness already measures everything needed, with ground truth:

1. **Fixture (1WCM phantom, `test=rec3D_backends` infrastructure):**
   truth-FSC and per-shell/radial LS comparison of (i) unsharpened map,
   (ii) global-B sharpening, (iii) post-hoc NU filter (current production),
   (iv) evidence-sharpened map — same inputs, same masks. This is the
   §11.3 "incremental value over post-hoc NU filtering" measurement the
   PCG doc already calls for, applied to sharpening.
2. **Heterogeneous-quality fixture:** the deliberately weak/coarse
   peripheral-domain phantoms from Gate C. The claim that distinguishes
   local from global sharpening is differential restoration: the weak
   domain sharpened to its own supported band, the strong domain to its
   finer one, neither over-amplified.
3. **Real data (bgal, streptavidin):** no truth, so judge by the
   established truth-free observables plus visual/model-map assessment
   against deposited structures; record before/after local spectra in the
   high- and low-confidence regions.
4. **Negative controls:** identity behavior at zero restoration gain;
   solvent regions must not brighten (the null bounds them); a
   deliberately mis-scaled target must be visibly rejected by the
   confidence caps (mutation-style check).

Acceptance thresholds recorded before the runs, per R9.

## 5. Sequencing

1. Land the outstanding Stage 6 items first (bgal ablation, Gate C
   variants, Gate A algebra tests). This note is parked until then.
2. First implementation: gains 2.2 with target 2.3(a), fixture validation
   (items 1–2, 4).
3. Local-B variant 2.3(b) if band coarseness limits (1); real-data pass
   (item 3).
4. Wilson target 2.3(c) only after the Wilson stage exists — do not build a
   parallel spectrum source for sharpening.

## 6. Relationship to the existing roadmap

| Component | Relationship |
| --- | --- |
| `Q_NU` replay (`pcg_priors.md` §5) | Same frozen evidence, same band frame; replay regularizes *in-solve*, sharpening restores *post-hoc*. Complementary, never combined implicitly. |
| Production NU filter | Sharpening generalizes it: graded gains instead of binary local cutoff selection, plus restoration. Long-term the NU filter is the `g`-restoration-off special case. |
| Wilson prior (§5.5) | Supplies target 2.3(c); sharpening is the lowest-risk consumer of the Wilson spectrum and a natural first validation of it. |
| Beyond-band retention watch item (Gate C record) | The Wiener layer is the principled shipped-map rolloff that closes it. |
| LocScale-2.0 [1] | Same operation class; confidence and target here are model-free from cross-half evidence instead of pseudo-atomic references. |

## References

1. A. Bharadwaj, R. de Bruin, and A. J. Jakobi, "Confidence-guided cryo-EM
   map optimisation with LocScale-2.0," *Nature Communications*, vol. 17,
   article 8778, 2026.
2. A. J. Jakobi, M. Wilmanns, and C. Sachse, "Model-based local density
   sharpening of cryo-EM maps," *eLife*, vol. 6, e27131, 2017.
3. A. Singer, "Wilson statistics: derivation, generalization and
   applications to electron cryomicroscopy," *Acta Crystallographica
   Section A*, vol. 77, pp. 472-479, 2021.
4. M. Beckers and A. J. Jakobi, "Confidence maps: statistical inference of
   cryo-EM maps," *Acta Crystallographica Section D*, vol. 76, pp. 332-339,
   2020.
