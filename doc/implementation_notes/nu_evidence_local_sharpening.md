# NU-evidence local sharpening: model-free LocScale from cross-half evidence

## Status

Proposal (2026-08-27), not yet scheduled. Recorded now so the idea is
specified while the Stage 6 NU replay work is fresh; implementation must not
begin before the direct NU-evidence prior clears its Gate C/D program
(`pcg_priors.md` §8, §10) — this is a postprocessing experiment that consumes
Stage 6 infrastructure, not a competitor to it. No solver, base-solve, replay,
or artifact behavior changes are proposed here.

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
