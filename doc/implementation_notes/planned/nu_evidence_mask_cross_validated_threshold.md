# NU-evidence mask: threshold chosen by cross-half prediction error

## Status

**Planning only (2026-10-01). No code is proposed for immediate implementation.**

The idea is supported by a synthetic check in numpy (section 4). Nothing has
been run in SIMPLE or on real half-map pairs. The first step, when wanted, is
a diagnostic table with no change of behaviour (section 6). Code references
are to `c961fb432`.

## 1. Question

The strength of the PCG solvent prior is chosen by cross-validation with the
NU objective (`estimate_solvent_prior_lambda`, decision log 2026-09-22), and
that works. The NU-evidence mask has been hard to get right. Can the same
kind of optimisation choose the mask?

## 2. What transfers from the solvent prior

Not the solver, the selection rule. A mask is an estimator choice. It can be
chosen the way the prior strength is: score what the consumer actually uses
against the other raw half, scan a very small family, take the minimum, and
flag a flat curve or a minimum on the grid edge.

The envelope is built today as follows (`simple_nu_filter_envmask.f90`,
[algorithm note](../../algorithms/nu_evidence_envelope_mask.md)):

```text
margin D = C_coarsest - min_c C_c, smoothed once at amsklp
  -> null median and MAD (whole observed support, or the dilation ring of
     the density envelope for a constrained base pair)
  -> threshold at nu_msk_sig MADs above the null median
  -> binary MRF -> component filter, hole filling, growth, cosine edge
```

The fragile parts are the null model in its two regimes, the validity gates
that go with it (solvent majority, shell sufficiency), and the fixed MAD
multiple. Section 6.7 of
[superseded NU-evidence envelope masking note](../rejected/nu_evidence_envelope_masking.md) already
suspected that the threshold depends on SNR.

## 3. Proposal

Keep the smoothed margin as the statistic that orders voxels, and keep the
morphological finish. Replace the null estimate and `nu_msk_sig` by a scan:

```text
for each threshold t on a small grid:
    m_t  = envelope at threshold t, with its soft edge
    P_e  = m_t * (NU-filtered even half) + (1 - m_t) * (coarsest candidate of the even half)
    P_o  = the same for the odd half
    J(t) = mean over the spherical support of
           H((P_e - O)/sigma(r)) + H((P_o - E)/sigma(r))
choose the minimum of J
```

`P` is the background clamp the envelope is used for: the selected label
inside, the coarsest bank candidate outside. `E` and `O` are the raw base
halves, `sigma(r)` the radial whitening profile and `H` the Huber loss, all
as in `image%nu_objective`.

Points that belong to the proposal:

- **Paired differences.** `J` is dominated by irreducible noise, and the
  differences between thresholds are below 0.1 % of it. Only per-voxel cost
  differences between two thresholds are informative, with a standard error
  from a block jackknife over the support.
- **Tie-break.** Among thresholds within one paired standard error of the
  minimum, take the tightest mask.
- **Flat curve.** When no threshold differs from the others, the objective
  does not define the mask and the density envelope is used. This replaces
  the null-validity gates.
- **Second scalar.** `amsklp` is coupled with the threshold and can be
  scanned jointly on a two-dimensional grid.
- **MRF.** Whether the binary MRF is still needed once the threshold is
  chosen this way is an open test.
- **Cost.** A few tens of objective passes on volumes the NU filter already
  holds. No solves and no new filtering.

The selection is in-sample: the margin and the score use the same pair. That
is acceptable for one or two scalars over about a million voxels, as it is
for the prior strength. It would not be acceptable for per-voxel freedom.

## 4. Synthetic check

Setup: 128^3 voxels at 1.5 A; a pseudo-atom "protein" filling 2.6 % of the
spherical support; a smooth 20 A "micelle" around its lower part; white
noise, independent per half; bank 20, 15, 12, 10, 8, 6, 5, 4 A; whitened
Huber unary; coarse-to-fine like-for-like label selection with Gaussian cost
smoothing. No Potts prior, no bank cap, no radial noise profile. Two to
three noise seeds per level.

| Half-map FSC=0.143 | Fixed 3 MAD: mask size, Dice vs protein | Cross-validated: mask size, Dice vs protein | Threshold chosen |
| --- | --- | --- | --- |
| about 6.7 A | 7.7 % of the sphere, 0.51 | 3.0 %, 0.93 | about 30 MAD |
| about 16 A | 4.6 %, 0.72 | 3.0 %, 0.91 | about 7 MAD |
| 21-24 A | 2.9 %, 0.86 | flat curve, erratic (3 to 30 %) | none |

What the check shows:

- The threshold that is right in MAD units changes several-fold with the
  noise level. A fixed multiple cannot transfer between datasets.
- In the first two rows the chosen mask coincides with the one that
  minimises the error against the true signal. Its gain over the 3-MAD mask
  is 2 to 4 paired standard errors, and the masks within one standard error
  of the minimum span 3.0 to 3.5 % of the sphere.
- In the third row nothing finer than the coarsest rung is reproducible. The
  curve is flat and the minimum is arbitrary, so the flat flag and the
  fallback are required.

Three variants were tried in the same setup and did not give a usable mask:

- **Zero as the background.** With `P = m_t * (NU-filtered half)` and zero
  outside, `J` rises with any masking at every noise level. The objective
  never prefers hard removal to a coarse low-pass, which is the principle of
  2026-09-02 (references are never multiplied by an envelope) found again
  from the data. It can therefore price a multiply-by-envelope mask, as the
  difference in `J` between multiplying and clamping, but it cannot choose
  one.
- **Zero as a bottom rung of the ladder.** A support defined as "any
  candidate beats zero" covered 35 to 66 % of the sphere, with 6 to 49 %
  occupancy in solvent far from the particle. Zero against 20 A is too
  weakly determined.
- **Local cross-half regression coefficient**, thresholded at one half. The
  result is haloed and includes the micelle (Dice 0.3 to 0.5 against the
  protein).

## 5. Limits

- The mask keeps the meaning it has today: structure finer than the coarsest
  rung is reproducible here. A detergent micelle stays outside. The mask
  does not become a valid PCG solve support, and it does not become a safe
  reference multiplier.
- It is still selected on cross-half agreement and must never reach the FSC
  (`automasking_policy.md`).
- Not tested: real noise (radial profile, deapodization), the Potts prior,
  the bank cap, weak or flexible domains, and constrained base pairs with
  zero/zero voxels inside the sphere.

## 6. First step, when wanted

A check-style table in the envelope writer or `nu_filt3D`, in the manner of
`pcg_solvent_check`, with no change of behaviour:

```text
threshold (MAD)   mask % of support   dJ +- SE (paired, against the minimum)
```

printed on real `_unfil` pairs (bgal, exp_gate, msp1, PfCRT and one
particle-count series), with the row nearest 3 MAD marked. It decides
whether real pairs show the interior minimum before anything is designed.

## 7. Related uses of the same scan

- A measured guard for `automsk=nu`: the price in `J` of multiplying the
  references by the envelope, against clamping the background.
- The solvent-prior weight map: the Otsu threshold and the logistic width
  could be scanned in the same way when the strength curve is flat ("the
  weight map, not the strength, is the limit").
