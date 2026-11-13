# The PCG order prior: one solvent and detergent prior for soluble and membrane proteins

## Status

**In the loop, 2026-11-08, untested on real data.** `pcg_order=yes` (with
`rec_backend=pcg`) puts a ridge on what lacks NU evidence of order into the
PCG re-solve of the base pair. Its ordered fraction and strength are chosen per
state and iteration in closed form. `solve3D` starts PCG and the prior with the
gold-standard phase (stage 6). Code: `simple_pcg_order_prior` (`src/main/volume/`),
the PCG master (`simple_rec3D_pcg_strategy`), the operator (`reconstructor_pcg`).
The offline program `simple_exec prg=order_prior` (`simple_order_prior_offline`)
checks the closed-form choice against full solves of a Fourier surrogate of the
operator. Unit sub-suites: `unit_reconstruction suite=order_prior` and
`suite=order_prior_offline`.

The order prior replaces the solvent prior `pcg_solvent` (decision log
2026-09-18 to 09-22) and no longer uses its density weight (section 3).
Sections 9 to 11 record the offline rounds that shaped it. This note takes over the first step of
[the cross-validated NU-evidence threshold](nu_evidence_mask_cross_validated_threshold.md):
the evidence mask is the prior's order weight.

## 1. The problem

`solve3D` with gold-standard matching from stage 6 gives usable maps of small membrane proteins. For the
smallest ones, though, the map does not separate ordered protein from disordered
detergent: the two are connected in the map and in every mask derived from it.
The solvent prior flattened solvent with a ridge outside a density weight `w`, a
logistic on `|x|` smoothed at `max(8 A, 2 x FSC0.143)` and Otsu-thresholded. At
that scale a micelle is as dense as the protein it surrounds, so `w` is about 1
on both and the micelle stays. The goal is an alignment reference without the
micelle. The final reconstruction never sees the prior.

## 2. What separates protein from detergent: order, not density

The NU evidence margin is about zero on a micelle and positive on ordered
protein. The margin is the gain of the best finer candidate over the 20 A
candidate in cross-half prediction, smoothed at `amsklp`
([algorithm](../../algorithms/nu_evidence_envelope_mask.md)). A micelle is
reproducible at 20 A and at nothing finer.

The order weight `w_ord` is a logistic on the margin. Its threshold is the
`(1-f)` quantile over the support, where `f` is the ordered fraction; its width
is the spread of the margin below the threshold. There is no binary MRF and no
morphology, because a prior weight needs no hard boundary. The weight is built
from both halves, like the NU label field, so a pair it shaped never reaches an
FSC.

## 3. The model

```text
penalty = lambda h (1 - w_ord) x^2,     the same weight for both halves
```

`h` is the operator's real-space diagonal, the data precision per voxel (the
mean of `D` over all shells). In these units one `lambda` means the same thing for
every molecule: a fully penalized voxel shrinks by `1/(1 + lambda)`. A strength
relative to the low-band data scale, the first version, acted 0.46 to 11.8 times
as strongly in `h` depending on the shape of `D(k)`, so each molecule needed its
own value (Hans, 2026-11-08: "we cannot design a system where we need to set a
different regularization parameter for different kinds of molecules").

The prior penalizes whatever lacks evidence of order: solvent and disordered
density such as detergent. The density weight `w` of `pcg_solvent` was in the
penalty (`1 - w w_ord`) until 2026-11-08. It leaked into the protein: 12-25 % of
the core's voxels had `1 - w > 0.2` (surface, thin and low-density features near
its Otsu threshold), which a strength that removes the micelle cuts to under
half. The maps looked heavily low-pass filtered (Hans). The order weight alone
puts that penalty above 0.2 on 0.7-7 % of the core and still penalizes the
solvent (mean `w_ord` 0.04-0.19 there). The density weight now only defines the
regions of section 4.

## 4. The choice of `f` and `lambda_rel`

Two regions are fixed from the prior-free pair before the scan:

- **Protein core:** support voxels that are dense (mean of the two halves'
  density weights >= 1/2) and whose NU label is at or finer than 1.5 x FSC=0.143.
- **Peri:** the dense voxels outside the core.

The signal of a region is its cross-half sum of even x odd, in which the halves'
independent noise cancels. For every `f` in `ORDER_FRAC_GRID` (1 % to 64 %) and
`lambda_rel` in `ORDER_LAMBDA_GRID` (0.1 to 20), the selector predicts the
signal each region keeps. Among the settings that keep at least 95 % of the
core's signal (`choose_order_setting`, shared with the offline check):

- with dense density around the core (a micelle, at least 1000 voxels), the most
  peri signal removed;
- with none (a soluble protein), the strongest `lambda_rel`, at its smallest
  fraction: solvent flattening;
- none qualifies, or the core is too small to check: no prior this iteration.

Auto is the default with `pcg_order=yes` (Hans, 2026-11-08, on the auto maps of
section 11: "absolutely fine"; to be revised on the loop tests). An explicit
`pcg_order_lambda`, for testing, is applied as given every iteration, at the
smallest fraction that keeps 95 % of the core (else the largest fraction).

Whole-map cross-half agreement cannot choose the removal. It counts every
reproducible voxel as signal, detergent included, so any removal costs it and it
always prefers to keep the micelle (Hans, 2026-11-07). Two guards failed before
the magnitude guard on the core signal, which removing noise does not lower and
removing density does:

- **Core objective within one paired SE:** it admitted a weak ridge on the core
  itself.
- **Removal within one SE:** solver-level differences of a fraction of a percent
  are many SE.

Neither score says anything about alignment. References were once multiplied by
an envelope and lost pose discrimination on PfCRT (2026-09-02). The ridge
attenuates rather than zeroes, and only the loop can test it.

## 5. In the loop

```text
base solve (prior-free) -> FSC, w per half
  -> select_order_prior (closed form) on the previous iteration's
     <vol>_nu_evidence.mrc and <vol>_nu_locres.mrc
  -> ridge on both operators -> the same cold solve again (the order pair)
  -> NU on the base pair, its label field composed from the prior'd pair
```

- **Closed form, banded:** each half is split into octave bands of Fourier
  shells (`[0,1]`, `[2,3]`, `[4,7]`, ...). Per voxel, band `b` shrinks by
  `D_b / (D_b + lambda h (1 - w_ord))`, where `D_b` is the operator's diagonal
  over the band's shells (`reconstructor_pcg%get_shell_diagonal`). The
  predicted half is the sum of the shrunk bands. Only the dense voxels are
  scanned, and no solve is repeated. One scalar `h` for all shells, the first
  version, missed where the core's signal sits in `D(k)` (section 11).
- **Lag one:** the NU evidence and labels the selector reads come from the
  previous iteration of the same state and box. The NU run of the current
  iteration needs the re-solved pair, and the re-solve needs the prior. The
  NU state filter writes the evidence map only under `pcg_order=yes`, so the
  first iteration with the key, and the first at a new box, run without the
  prior and say so in the log. Without NU filtering the prior never acts; a
  standalone `reconstruct3D` uses it only on evidence an earlier run left in
  its directory.
- **Gold standard and outputs:** the prior-free pair keeps the FSC, the NU
  competition, the evidence and the `_unfil` files. The prior'd pair is written
  as `<half>_order`; its ML replay is the state map, and the NU references are
  composed from it. The weight is written as `pcg_order_weight_stateNN.mrc`,
  overwritten every iteration.
- **Final reconstruction:** prior-free (`simple_final_rec` forwards no
  `pcg_order`).
- **Exclusions:** `atm_prior=yes` and `pcg_order=yes` are exclusive.
  `reconstruct3D_wave` and `refine3D_states_cavgs` refuse the prior.
- **solve3D:** gridding through stage 5, PCG and the prior from stage 6
  (`PCG_START_STAGE`), independently of `gold_std_stage`. The PCG accumulation
  and solves are the cost, so starting them late is the saving.
- **Trailing chains:** they are per backend. The stage-boundary reconstruction
  takes the next stage's backend and seeds the PCG chain at full mass. The prior
  acts after the blend, on the restored pair.

## 6. The offline surrogate of the PCG re-solve

Offline there are no accumulators, so the operator is modelled. The PCG base
operator is close to Fourier-diagonal: CTF-squared, sigma2-whitened sampling
density. For an unregularized reconstruction its diagonal is the inverse of the
half-map noise power, which the pair measures:

```text
sigma_n^2(k) = N <|F(E - O)|^2>_shell / 2     (fft divides by N)
D(k)         = 1 / sigma_n^2(k)
(F^-1 D F + lambda h diag(1 - w_ord)) x = F^-1 D F x0     per half
h            = mean D over all shells
```

The shared CG engine (`pcg_solve`) solves it with a two-region preconditioner,
`sqrt(v) D^-1 sqrt(v) + sqrt(1-v) (D + lambda h)^-1 sqrt(1-v)` with `v = w_ord`,
which is exact inside a region of one kind.

With `D` constant, the operator reduces to the closed form (unit tests). Where
`D(k)` falls with `k`, a ridge acts as a local low-pass rather than a mask. The
re-solve then keeps the low-resolution density the data determine and
suppresses the rest, which the closed form cannot see. Shells beyond the last
one with data are held at zero; such shells occur in a pair Fourier-padded from
a crop box.

## 7. The program

```text
simple_exec prg=order_prior vol1=<odd>_unfil.mrc vol2=<even>_unfil.mrc smpd=<A> mskdiam=<A> \
    [outvol=<name>.mrc] [amsklp=8] [pcg_order_lambda=<x>] nthr=<n>
```

Feed the unregularized pair: an ML-regularized pair flattens both the noise model
and the evidence. `mskdiam` must leave solvent around the micelle, because the
NU noise level and the evidence need it.

The program runs in this order:

1. The density weights (the regions).
2. NU on the pair, capped at FSC/1.5 as in assembly, giving the evidence margin
   and the label field.
3. The closed-form choice and its tables.
4. The same grid solved through the surrogate (warm-started along
   `lambda_rel`), chosen by the same rule.
5. Agreement or disagreement of the two choices.
6. The surrogate pair at the closed-form choice, which is what the loop ships.

Outputs, with `<name>` from `outvol` (default `order_prior`):

| file | content |
| --- | --- |
| `<name>_base.mrc`, `_base_nu_filt.mrc` | the input pair average, and its NU filtering |
| `<name>_nu_evidence.mrc`, `_nu_locres.mrc` | the smoothed evidence margin; the label field's local resolution |
| `<name>_density_weight_even/odd.mrc` | per-half density weights (they define the regions) |
| `<name>_order_weight.mrc` | `w_ord` at the chosen fraction (1 minus it is what the ridge acts on) |
| `<name>_order_even/odd.mrc`, `_order.mrc` | the prior'd pair at the closed-form choice |
| `<name>_order_nu_filt.mrc` | the base pair's label field composed from the prior'd pair: the reference the matcher would see |

## 8. What to read

- **In the loop:** one `PCG ORDER PRIOR` line per state and iteration, giving
  the ordered fraction, `lambda_rel`, the peri signal removed and the core kept,
  or `prior off this iteration` (no evidence yet, or no setting keeps the core).
  Read it next to the `_order` halves and the order weight. Diagnostics added
  2026-11-13 to find why the NU references looked low-passed with
  `automsk=no`:
  - the closed-form tables;
  - the spread of `1 - w_ord` over the core, the dense density around it and
    the rest of the domain;
  - the choice a guard on every octave band to FSC=0.143 would make (the
    loop still uses the total);
  - after the re-solve, a table per octave band: the band's share of the core
    signal, `D_b/h`, `lambda_rel/(D_b/h)` (a voxel with penalty `p` keeps
    `1/(1 + p x` that`)` of the band), and the core kept and density around it
    removed, predicted and measured on `<half>_order` against the base pair.

  `pcg_order_regions_stateNN.mrc` maps the regions (2 core, 1 dense around it).
- **Offline:** the two `CLOSED FORM` tables and the two `operator` tables, then
  `the operator at the closed-form choice removes ... keeps ...` and whether the
  two choices agree. The `BAND GUARD` choice and the per-band table (measured on
  the surrogate's pair) are printed as in the loop. A steep `D(k)` (streptavidin, embb, exp_gate) is where they
  can differ.

## 9. First round: three models by whole-map agreement (2026-11-07)

Three models were scanned, with every setting chosen by whole-map agreement:

- `soluble`: the density ridge alone.
- `ordered`: a ridge outside the order weight.
- `membrane`: the density ridge plus a gradient penalty
  `lambda_d w (1-w_ord) |grad x|^2`, which smooths the micelle rather than
  removing it.

The inputs were gridding `_unfil` pairs on a sphere in
`~/Processing/solvent_prior/<pair>` (`mskdiam` 80, 160, 170, 150, 155, 150 A;
agent-runs 1505-1510, Release, uncommitted tree on 913d50765):

| pair | FSC=0.143 | `D(k)` low band / mean | solvent level | soluble `lambda_rel` | membrane - soluble | ordered core (`f = 0` against the minimum) |
| --- | --- | --- | --- | --- | --- | --- |
| streptav (soluble) | 3.1 A | 11.8 | -0.18, 17 SE | 0.5 | +0.0004 +- 0.0005 | none: collapses to soluble |
| PfCRT | 5.35 A | 1.15 | 2 SE | 5 | -0.033 +- 0.003 | 0.5 % (within 2 SE) |
| PorF | 4.5 A | 0.46 | -0.33, 10 SE | 2 | -0.036 +- 0.003 | 2-4 % (5.6 SE) |
| embb | 3.9 A | 9.3 | +0.25, 17 SE | 1 | -0.009 +- 0.001 | 4 % (4 SE) |
| exp_gate | 3.9 A | 14.8 | +0.09, 14 SE | 1 | -0.006 +- 0.001 | 8 % (6 SE) |
| tat | 3.8 A | 1.12 | +0.33, 14 SE | 0.1 (edge) | -0.108 +- 0.005 | most dense density ordered |

By agreement, `membrane` won on every membrane protein and `ordered` lost on all
six. No prior's NU-composed references beat plain NU references.

Hans's review of the maps:

- streptavidin's prior'd map was excellent: a clean background and a low-passed
  surface, a reference for alignment rather than a final map;
- `ordered` was worst and `membrane` best;
- on PorF and embb the whole micelle was still there in `membrane`;
- on PfCRT the micelle took over the protein: the chosen ordered core was 0.5 %,
  and the protein was smoothed as detergent.

The scoreboard measured the wrong thing (section 4). Smoothing keeps the micelle,
so agreement favoured it, and `ordered` was held to the mildest removal.

Two other findings are kept:

- `D(k)` has no single shape: it is steep, flat or rising depending on the pair.
- Five pairs carry a solvent level that both halves share, and a ridge toward
  zero removes it.

## 10. The single model offline (2026-11-07)

The model was the density ridge plus a separate detergent ridge on
`w (1 - w_ord)`, chosen by full surrogate solves. Agent-runs 1525-1530 (Release,
uncommitted tree on 913d50765), 16 threads; maps in
`~/agent_work/SIMPLE/out_<run>`. The PfCRT pair here is a newer reconstruction
(FSC=0.143 4.2 A) than the one in section 9.

| pair | core, % of support | dense outside the core | chosen `f`, `lambda_m` | signal removed outside the core | core signal kept | run time |
| --- | --- | --- | --- | --- | --- | --- |
| PfCRT | 5.0 | 4.6 % | 16 %, 50 (grid edge) | 75 % | 99.6 % | 4.6 min |
| PorF | 8.5 | 4.1 % | 16 %, 2 | 47 % | 95.7 % | 2.3 min |
| embb | 9.2 | 8.0 % | 16 %, 20 | 95 % | 95.3 % | 7.9 min |
| exp_gate | 11.9 | 6.1 % | 32 % (edge), 50 (edge) | 53 % | 97.3 % | 17.4 min |
| tat | 9.0 | 18.6 % | none: below 32 % no setting keeps the core | | | 2.3 min |
| streptav (soluble) | 21.1 | 0.01 % | 32 %, 1 | 50 % of almost nothing | 96.2 % | 1.2 min |

Hans on PfCRT:

- at `f` 16 %, `lambda_m` 5 under the earlier selection (agent-run 1519):
  "excellent, the best solvent flattening I have ever seen";
- at the stronger choice of run 1526 (`lambda_m` 50): "maps and weights look
  good".

On every membrane pair, the ordered fraction must reach 8-16 % before the order
weight covers the core. Below that the ridge eats the core: under 85 % of its
signal is kept at 4 %.

What carried into section 4:

- **One ridge with one strength:** the integration needs a cheap in-loop
  choice, and a single ridge has one.
- **The net-signal condition:** streptavidin lost 3.8 % of its core with
  nothing to remove.
- **The grid:** extended to 48 % and 64 %, because tat needs `f` above its
  dense fraction of 28 %.

Still open: at 3.8 A the core (labels <= 5.7 A) leaves protein at 6-10 A
outside it, which counts as removable (tat: twice the core).

## 11. The closed form against the operator (2026-11-08)

### Scalar closed form

The first in-loop version shrank each voxel by `h / (h + lambda (1 - w w_ord))`,
with one `h` (the mean of `D` over all shells) and `lambda_rel` from 0.5. The six
pairs were run through `order_prior` (agent-runs 1536-1541, Release,
uncommitted tree on 913d50765, 16 threads):

| pair | `D` low band / mean | closed-form choice: removed, core kept | operator's choice | operator at the closed-form choice: removed, core kept |
| --- | --- | --- | --- | --- |
| PfCRT | 1.28 | 16 %, 1.0: 59.5 %, 95.9 % | 16 %, 0.5 | 52.0 %, 94.0 % |
| PorF | 0.46 | 16 %, 2.0: 49.1 %, 96.8 % | density ridge 0.5 | 59.3 %, 86.8 % |
| embb | 9.32 | density ridge 0.5: core 67.2 % | density ridge 0.5 | 53.7 %, 75.1 % |
| tat | 1.12 | 16 %, 0.5: 51.2 %, 96.3 % | density ridge 0.5 | 59.7 %, 76.4 % |
| streptav | 11.76 | density ridge 0.5: core 67.8 % | density ridge 0.5 | 53.6 %, 81.3 % |
| exp_gate | 14.76 | density ridge 0.5: core 48.0 % | density ridge 0.5 | 52.2 %, 65.3 % |

Two failures:

- **The scalar `h` misses where the core's signal sits.** The core's cross-half
  signal is at low and mid resolution. Where `D` is below `h` there (tat at
  8-17 A, PorF below 40 A), the operator removes much more of the core than
  predicted. Where `D` falls steeply, it removes less.
- **The fallback broke the guard.** The weakest strength was applied even when it
  failed the core guard: streptav, embb and exp_gate would have lost 20-35 % of
  the core. With steep `D` (low band 9-15 x the mean), `lambda_rel` 0.5 is
  already 5-7 x the real-space diagonal.

The fixes: the banded closed form (section 5), the off rule (section 4), and the
grid extended down to 0.05.

### Banded closed form

The same six pairs with the banded closed form, the off rule and the extended
grid (agent-runs 1548-1553, Release run 1547, same tree, 16 threads):

| pair | closed-form choice: removed, core kept | operator's choice: removed, core kept | operator at the closed-form choice | run time |
| --- | --- | --- | --- | --- |
| PfCRT | 16 %, 1.0: 53.8 %, 96.1 % | 16 %, 0.5: 36.4 %, 96.0 % | 52.0 %, 94.0 % | 3.2 min |
| PorF | 16 %, 1.0: 50.2 %, 96.2 % | 16 %, 0.2: 16.3 %, 97.6 % | 46.0 %, 91.2 % | 1.7 min |
| embb | 16 %, 0.05: 22.5 %, 97.1 % | density ridge 0.05: core 95.1 % | 21.1 %, 93.9 % | 4.9 min |
| tat | 16 %, 0.2: 37.7 %, 96.9 % | 16 %, 0.05: 12.4 %, 95.6 % | 36.4 %, 86.4 % | 1.5 min |
| streptav | density ridge 0.05: core 95.9 % | density ridge 0.05: core 97.3 % | 2.9 %, 97.3 % | 0.7 min |
| exp_gate | none | none | | 9.9 min |

The off rule now holds: no pair loses a third of its core. The closed form
still underestimates the core loss on the membrane pairs. Near the guard, the
operator's loss is about 2 x the predicted one on PfCRT, PorF and embb, and
4.5 x on tat; on streptavidin it is 0.7 x. The factor is roughly constant along
`lambda_rel` for each pair.

The bands are not the cause. On tat the density ridge alone (`f = 1`, `lambda_rel`
0.2) keeps 94.2 % under the operator against 97.6 % predicted. A low band's
content spans the protein and its micelle, so a ridge on the micelle lowers it
inside the core too. A per-voxel shrink cannot see this nonlocal effect.

### Order weight alone, strength in `h` (2026-11-08)

The penalty `lambda h (1 - w_ord)` on the six pairs, in auto and at one strength
for all, `pcg_order_lambda=5` (agent-runs 1564-1575, Release run 1563, 16
threads). Under 1 % of the core's voxels now carry a penalty above 0.2 on PfCRT,
PorF and streptavidin, and 3.9 % on tat; it was 12-25 % with the density weight.
The operator at the closed-form choice:

| pair | `lambda` 5: fraction, removed, core kept | auto: closed-form choice | auto: removed, core kept |
| --- | --- | --- | --- |
| PfCRT | 16 %, 70 %, 98.3 % | 8 %, 2.0 | 78 %, 92.6 % |
| PorF | 16 %, 73 %, 91.9 % | 16 %, 20 | 81 %, 88.7 % |
| streptav | 48 %, -, 104 % | 48 %, 20 | -, 106 % |
| embb | 16 %, 73 %, 88.3 % | 16 %, 20 | 92 %, 84.8 % |
| tat | 16 %, 93 %, 64.6 % | 16 %, 5.0 | 93 %, 64.6 % |
| exp_gate | 32 %, 31 %, 90.0 % | 32 %, 20 | 51 %, 85.4 % |

The closed form keeps 95-100 % of every core; the operator keeps less on the
membrane pairs, most on tat. The voxels of the core are no longer penalized, so
the loss is the nonlocal one of section 11: on tat the dense density around the
core is twice the core and probably mostly protein below 6 A. On streptavidin the
operator raises the core's signal (above 100 %): flattening the solvent moves low
band content into the protein. Maps: `~/agent_work/SIMPLE/out_<run>/<pair>_<auto|l5>_order_nu_filt.mrc`.

## 12. With `automsk=yes` (2026-11-13)

`refine3D_local` with `pcg_order=yes automsk=yes` cut PorF in half. What
survived had the best resolution Hans had seen on it. The solve support was the
density envelope of the published reference, which under `ml_reg=yes` is the
ML replay of the prior'd pair. The ridge's removals fell below the envelope
threshold and were solved as zero from then on. The support is now made from the
previous iteration's prior-free pair, for every `automsk` run
(`reconstruct3D_pcg_policy.md`, section 3). The division of labour:

- the envelope, from prior-free data, owns the solvent as a hard support;
- the order prior owns the disordered density inside it (the micelle) as a soft
  ridge, its fraction and evidence null taken inside the envelope;
- the NU references are masked by the current prior-free pair's envelope.

## 13. Next

The question is whether the prior improves the alignment, not what it deletes.
From `gold_std_stage` (default stage 6) the halves align independently, and the final map is
prior-free, so the loop judges the prior. The order weight is shared by the
halves, like the NU label field: the evidence of order is a cross-half quantity,
and without the per-half density weight (section 3) nothing in the penalty is
per half. It is smooth (8 A, logistic), a shared soft solvent mask rather than
strict two-fold independence; accepted for the tests (Hans, 2026-11-08). An FSC
gain only at low resolution would point to it. Strict independence would need
quarter reconstructions, each half's weight from the evidence between its two
quarters.

Test on PfCRT, PorF and streptavidin in auto (`pcg_order=yes` alone), the setting
whose offline maps Hans approved. Each dataset gets a `pcg_order=no` arm on
the same build, since PCG now starts at stage 6. Compare the final prior-free
maps and their FSC, the stage 6-8 references, the per-iteration `PCG ORDER
PRIOR` lines, and the run time.
