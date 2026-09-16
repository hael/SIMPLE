# Retiring the NU shell walk: a dense static ladder with the regularized pair in the competition

## Status

Design, 2026-09-16, revised from the 2026-09-06 proposal of the same file.
**Awaiting Hans's approval; no code is proposed for immediate implementation.**
The 09-06 text is superseded in full (its §1–2 described the ML pair as a
second, separately solved estimator; since 2026-09-14 it is not, see §1.2).

Decision reached in discussion (Hans, 2026-09-16): retire the `nu_refine=yes`
shell walk from `refine3D_auto`; run the static-bank competition that every
`abinitio3D` run has used, with the ladder densified at its fine end, capped
by the gold-standard FSC, and with the regularized pair as a competitor. This
note is the record to think against before implementing. Section 8 lists the
decisions still open.

## 1. Context

### 1.1 What the two modes are today

Both backends run one post-hoc NU competition per state and iteration on the
unregularized base pair: a bank of low-passed candidates, a whitened Huber
cross-half prediction error per voxel and candidate (raw even vs filtered odd
and vice versa), candidate-scale smoothing, coarse-to-fine like-for-like
selection, ordered-label Potts cleanup, synthesis of the `_nu_filt` even/odd
references from the selected labels, the `_nu_locres` map, and the handoff of
the finest selected label (with the 1% signal-voxel floor) as the next
iteration's matching band (`nonuniform_filtering_policy.md` §8–12).

`nu_refine=no` (static-bank mode, every abinitio3D NU stage): the bank is the
ladder `lowpass_limits = [20,15,12,10,8,6,5,4]` Å, truncated at
`fsc/NU_BANK_FSC_HEADROOM` (1.5× finer than the pair's FSC=0.143, never fewer
than two rungs), and with `ml_reg=yes` the regularized pair *replaces the
finest retained rung* when its FSC=0.143 resolution is finer than that rung
(`setup_nu_dmats`, `aux_replacement_idx`; the handoff then reports the
auxiliary's resolution for voxels that selected it).

`nu_refine=yes` (refine3D_auto default since 2026-09-06): full ladder, no FSC
cap (2026-09-11), and a shell walk that challenges one Fourier shell at a time
beyond the finest rung, accepting a shell only on a significant majority of
the frontier (z ≥ 3, 2026-09-13), with frontier bookkeeping, a per-challenge
filtered pair and unary, candidate thinning under `NU_DMAT_CANDIDATE_CAP` (24),
ordered-label cleanup with walked shells sharing the finest rung's Potts
coordinate (2026-09-13), and the walk depth persisted for restarts.

### 1.2 What changed since the 09-06 proposal

- The regularized pair is no longer a separately iterated solve. Since
  2026-09-14 it is the closed form `base(v) · (rho+floor)/(rho+floor+P_τ)`
  on the padded lattice (`reconstructor_pcg%shrink_by_ml_prior`; on gridding
  `add_invtausq2rho` is the same division): a deterministic, voxelwise,
  sampling-aware Wiener filter of the base pair whose shell shrinkage is the
  FSC (`W_k = τ·SSNR/(1+τ·SSNR) = FSC_k` at τ=1). The "ML estimator" of the
  old note is a *filter shape*, not an estimator, and it is the MMSE
  cross-half predictor per shell by construction.
- The walk's stopping criterion is the local FSC=0.5 crossing: adding a
  shell at full weight lowers the cross-half error only where local
  half-map SNR > 1. So the walk is a correct estimator of where a *hard*
  cutoff should stop, and what it leaves on the table is the filter shape
  above that point, not the stopping point (2026-09-15 discussion).
- Record of the walk on real targets: aldolase run 15/16 — first challenge
  rejected every iteration (27% wins, z ≈ −140), contribution nil; aldolase
  38.8k before the majority test — walked to 3.37 Å on coin flips, band
  finer than the FSC, flat FSC, decaying cFAR; PfCRT — the July result
  predates the walk and the 2026-09-11 diagnosis found neither band nor bank
  to be the limiter. Three corrective changes in one week (majority test,
  ladder-only Potts pricing, FSC-gate removal). The static-bank competition
  has carried every abinitio3D run.
- `refine3D_auto` now has a gold-standard, envelope-corrected FSC
  (`envfsc=yes`) to anchor a cap on, which abinitio3D does not.

### 1.3 What the walk buys, and what a dense ladder reproduces

Two things, both at the fine end: hard cutoffs at shell granularity between
the last rung (4 Å) and the resolution limit, and an extension rule driven by
local evidence rather than a global cap. A ladder densified below 6 Å gives
the granularity with the proven selection machinery and none of the frontier
logic. What it cannot do is extend a region past `fsc/1.5` of the global
FSC=0.143 — a core 1.5× better than the global number — and there the
regularized candidate's roll-off is already doing the right thing.

## 2. Decision

For `refine3D_auto`:

1. Retire the shell walk. `nu_refine=yes` is no longer a mode; the static-bank
   competition runs every iteration.
2. Densify the ladder at its fine end, generated from the box, not a fixed
   list (§3.1).
3. Cap the bank at `fsc/NU_BANK_FSC_HEADROOM` from the gold-standard FSC=0.143
   of the base pair, as static mode does now (§3.2).
4. The regularized pair joins the bank as a candidate placed at its own
   FSC=0.143 resolution, rather than replacing the finest rung (§3.3).
5. Handoff, `_nu_locres`, evidence envelope, Potts cleanup and synthesis are
   unchanged in kind; the Potts coordinate is redefined so dense rungs do not
   change the price of a resolution jump (§3.4).

`abinitio3D` uses the same bank (decided 2026-09-16, one implementation):
its NU stages get the dense ladder and the finest-member rule. This does
not change its top end -- the regularized pair already replaces the finest
rung and already sets the handoff when selected -- and it removes the hard
rungs finer than the (non-gold-standard, inflated) FSC=0.143 that the
`fsc/1.5` cap admits today, which is the safer direction there. Note that
the 4.5 A ladder bound guards only the FSC=0.5 promotion of the planned
non-NU stage limits; the NU-stage handoff has carried no ceiling since
2026-09-08 and stays uncapped. Should the unified bank show the abinitio
band running ahead on an inflated FSC, a ceiling on the regularized
candidate's content extent in abinitio3D's NU stages is the one-line guard
to add, on evidence. abinitio3D runs are therefore a regression gate of
this change (§6), not a later adoption.

Soft (FSC-shaped, anchored) synthesis in place of the hard cutoff at the
selected label, discussed 2026-09-15, is **not** part of this change. It is a
separate later step on top of the dense ladder (§7).

## 3. Design

### 3.1 The ladder (decided 2026-09-16)

Coarse rungs as now: 20, 15, 12, 10, 8, 6 Å. Below 6 Å, hard rungs at a
constant shell spacing of the current (cropped) box, nominally
`NU_LADDER_FINE_STEP = 2` shells, from the 6 Å shell to the fine bound of
§3.3 (the base pair's FSC=0.143 rung when `ml_reg=yes`, the `fsc/1.5` cap
otherwise). The list is built once per call from `box`, `smpd` and the
bound; `lowpass_limits` stays as the coarse part.

Bank budget: at most `NU_BANK_MAX_MEMBERS = 16` candidates in total (six
coarse rungs, the fine rungs, and the regularized candidate when present).
When 2-shell spacing would exceed the budget the step is widened to fit
(`step = max(2, ceil(nshells_fine / n_fine_slots))`), so the fine rungs
always span the whole 6 Å-to-bound range at the coarsest spacing the
budget allows, rather than truncating and leaving the regularized
candidate alone to cover the gap. Aldolase (box 192, 1.3 Å; 6 Å = k42,
FSC=0.143 = k68): 26 shells into 9 slots → step 3, rungs at k 45..69;
PfCRT (box 300, 1.035 Å; 6 Å = k52, ~4 Å = k78): step 3; a 3 Å map at
box 300 (k52..k103): step 6. The budget is a memory cap, never a failure:
one packed unary column is `n_mask × 4 B` (aldolase ~1.6 MB, box 300
~6 MB per candidate) and one filtered pair on scratch disk per candidate
(2 × box³ × 4 B; 216 MB per candidate at box 300, ~3.5 GB for 16).
`NU_DMAT_CANDIDATE_CAP` (24, sized for the walk's window) is replaced by
this budget.

### 3.2 The cap (decided: keep 1.5)

`fsc/NU_BANK_FSC_HEADROOM` (1.5) of the gold-standard FSC=0.143 of the base
pair, as in static mode. With `ml_reg=yes` the fine bound of the hard
rungs is the FSC=0.143 rung itself (§3.3), which is always coarser than
the cap, so the cap only acts when `ml_reg=no` (hard rungs to `fsc/1.5`,
today's static behaviour) and in the coarse-resolution regime where it
drops rungs the pair cannot support (never fewer than two).

### 3.3 The regularized pair in the competition (decided 2026-09-16)

Rule: with `ml_reg=yes` the regularized pair is always the FINEST member of
the bank and replaces every hard rung at or finer than the pair's
FSC=0.143 resolution. Hard rungs therefore run from 6 Å up to the last
rung coarser than the FSC=0.143 shell; the regularized candidate follows
as the last member, with its FSC=0.143 resolution as its label resolution
(handoff, `_nu_locres`) and the Potts coordinate of that resolution (§3.4).
This is the current static aux-replacement rule ("replace the finest
retained rung") generalized to a dense ladder, and it keeps the bank to
one candidate per coordinate: no two members compete at the same
resolution, so the like-for-like coarse-to-fine selection needs no new
tie rule.

What is given up relative to "add it beside a hard rung at the same
resolution": a best-resolved core that would prefer full amplitude to the
FSC=0.143 shell (or a finer hard rung up to the cap) over the global
Wiener shrink there. Under the previous ladder it never had that option
either (the 4 Å rung was replaced), and the regularized candidate's
`rho`-aware shrink is the better estimator of that region in expectation;
if the diagnostics of §5 show the finest hard rung winning a large,
coherent population right below the regularized candidate, the "add
beside" variant is the thing to try, as a separate experiment.

What this candidate is: the base pair with the global Wiener roll-off,
voxelwise `rho`-aware. Against the hard rung one step coarser it carries
the shells between local FSC 0.5 and 0.14 at their proper weight instead
of nothing, and shrinks undersampled voxels by their own `rho`. It is
evaluated with the same unary, smoothed at its own scale, and selected
coarse-to-fine like the other members; no extra margin in the first
experiment.

Half independence as in the old note §3.3: a voxel selecting the regularized
candidate takes its even output from the regularized even map and its odd
from the regularized odd; the merged `_nu_filt` remains the average of the
two synthesized halves.

### 3.4 Potts coordinate

The hinge prices coordinate jumps, and with ladder-position coordinates a
dense fine end would price a 6→4 Å transition as several steps instead of
two and over-smooth the fine end. Coordinate = position on the 8-rung
ladder, interpolated in `log(1/resolution)` for the fine rungs and for the
regularized candidate at its FSC=0.143 resolution, so the price of a given
resolution jump is what it is today whatever the spacing the budget chose.

### 3.5 Matching band: the handoff follows the reference content

This is the part of the scheme that must NOT be carried over from
abinitio3D unchanged (Hans, 2026-09-16). abinitio3D is right to be
restrictive: its FSC is not gold-standard, its stage limits are a planned
ladder with FSC=0.5 promotion only at stage boundaries, and its NU handoff
is the finest populated rung of a coarse ladder -- all of which protect a
non-gold-standard search from ratcheting on noise it fitted itself. A
gold-standard high-resolution refinement has the opposite need: every shell
the reference carries with signal should be in the objective, as RELION
does by matching to the FSC=0.143 resolution of its regularized reference.

The principle: the band is set by the *content extent* of the selected
references, not by a label index. A hard rung carries content exactly to
its cutoff. The regularized candidate carries content with the Wiener
weight `W_k = FSC_k` out to where the FSC vanishes, i.e. usefully to the
pair's FSC=0.143 resolution (`W ≈ 0.14`; the prior floor at FSC 0.001 makes
everything beyond effectively zero). Matching beyond a hard rung would
match against an empty reference and is pointless; matching to the
regularized candidate's extent is matching a Wiener-weighted reference,
the standard practice.

Rule: for every label, define its content extent -- the cutoff for a hard
rung, the base pair's FSC=0.143 resolution for the regularized candidate.
The handoff is the finest content extent among labels holding >= 1% of the
signal voxels (at that label or finer; `NU_ALIGN_LP_MIN_SIGNAL_PCT`, the
existing floor, decided kept), capped by a user `lpstop`, with no headroom
beyond it (decided). Since the regularized candidate is the finest member
and is expected to occupy the top of the field wherever the local FSC
curve resembles the global one, the band will normally be the
gold-standard FSC=0.143 resolution -- the same band the FSC-driven non-NU
path uses (`lplim_crit=0.143`), so `filt_mode=none` and
`filt_mode=nonuniform` agree on the band and differ only in the
references -- while regions the competition filtered coarser still enter
the objective only to their own cutoff through the reference itself. If
the regularized candidate holds fewer than 1% of the signal voxels the band
falls back to the finest hard rung, which is the correct reading of "the
data do not support the global roll-off anywhere". The 1% floor is the
same floor every label is subject to; the regularized candidate being
one member deciding the band is no different from the finest hard rung
doing so today, and a 1% coherent population at the global resolution is
exactly the case where the FSC=0.143 band is warranted there.

Consequences to note: the band can now advance every iteration as the
gold-standard FSC advances, coupled to poses through the reference as
in any refinement; the gold-standard split is the guard against the
band/FSC co-adaptation loop, as it is in RELION. `incrreslim` (ten shells
beyond the criterion on the FSC-driven path) is not applied to this
handoff in the first experiment; whether headroom beyond FSC=0.143 buys
anything with a regularized reference that is already ~zero there is an
open question (§8). The registration pass keeps its own band (FSC=0.8 via
`lpstop`) and is unaffected. `_nu_locres` keeps reporting label
resolutions, not content extents: it describes where the map is resolved,
the handoff describes what the matcher may use.

With the soft-synthesis step of §7 every label becomes an anchored
roll-off and the same rule gives band = anchored 0.143 point of the finest
populated label; the rule above is its special case for hard rungs plus
one global roll-off.

### 3.6 Evidence bands and envelope

`NU_EVIDENCE_BAND_LIMITS = [20,12,8,5]` and the geometric extension "over
walked candidates" lose their extension case and stay at the static four
bands. The NU evidence envelope (automsk=yes) is built from the same unary
bank and is unaffected in kind. Both to be checked, not assumed, when the
walk code is removed.

## 4. What is removed

- The walk: `simple_nu_filter_extend.f90` (831 lines: challenge, majority
  test, frontier, thinning, walked-label cleanup), the walk-depth
  persistence (`write_nu_highres_steps_for_state`, `n_highres_steps`), the
  `NU BANK UNCAPPED` path in `init_nu_filter`, `nu_highres_extension_stats`
  and the `NU SHELL WALK` / `accepted shell steps` log lines,
  `NU_HIGHRES_EXTENSION_*` constants, `NU_HIGHRES_EXTENSION_RETAIN_STRIDE`.
- `nu_refine` as a mode: the key is removed (decided); `l_nu_refine`
  branches in `simple_nu_state_filter`, `simple_commanders_rec_distr`,
  `simple_final_rec`, the abinitio controller (which already emits `no`),
  `refine3D_auto` defaults and UI, and the parameter/parse/UI
  registrations. An unknown key is accepted by the parser, so a stale
  `nu_refine=yes` on a command line is silently ignored; the removal is
  logged in the decision log and policies.
- `nu_static_aux_replacement` (which excluded the aux under `nu_refine=yes`)
  collapses to "regularized pair is the finest member whenever `ml_reg=yes`".
- `postprocess_nu` (decided): the same machinery -- dense ladder within the
  budget, regularized pair as finest member when a regularized pair is
  supplied, same handoff-free synthesis -- through the shared
  `nonuniform_filter_state`/`nu_filter_vols` path; its call into
  `extend_nu_filter_highres_shells` goes with the walk.
- Policy text: `nonuniform_filtering_policy.md` §8 (second paragraph), §10
  (High-Resolution Extension, entire), §12 (walk references), §13
  (`nu_refine=yes` default); `refine3D_auto_policy.md` NU/FSC paragraph and
  defaults; `reconstruct3D_pcg_policy.md` nu_refine paragraph;
  `pcg_backend_overview.md` §5 and log-line list; `pcg_decision_log.md`
  entry. `postprocess_nu` (which shares `extend_nu_filter_highres_shells`)
  needs its own decision (§8).

## 5. Diagnostics

Per state and iteration, on the master:

- retained rung list (count, coarsest, finest, cap) — one line;
- selected-label histogram including the regularized candidate as its own
  row (`Source = ML`), with voxel counts and percentages, as the current
  `NU LOW-PASS ASSIGNMENTS` table does for `Base` rows;
- for the regularized candidate: percentage of voxels where it was selected,
  and of those, the median smoothed-cost improvement over the hard rung at
  the same coordinate; and its iteration-to-iteration occupancy (the
  fixed-point oscillation check of the old note §7);
- handoff and its rung; `_nu_locres` unchanged in meaning;
- runtime and peak RSS of the competition against the walk's.

## 6. Validation and gates

A/B on identical inputs, iteration counts and matching settings, current
default (`nu_refine=yes`) against the dense static bank, on:

- aldolase (18.6k, box 192): the null control — at the ceiling, the
  expectation is identical FSC, cFAR and pose jitter, with the regularized
  candidate occupying the top labels;
- PfCRT (box 300): the case where band and reference shape have mattered;
- streptavidin: the small, low-SNR case;
- bgal: the > 4 Å-resolution case where the walk should have mattered most
  and where the dense ladder must reproduce what it gave between 4 Å and
  the map's limit.

Pass: FSC=0.143 and FSC=0.5 equal or better; cFAR trajectory not decaying
faster; inter-iteration pose distance not larger; regularized-candidate
occupancy spatially coherent and stable across iterations; runtime and
memory acceptable on both backends. abinitio3D regression gate on the
same targets from the `~/for_claude` restart inputs (bgal, PfCRT,
streptavidin canonical): stage-by-stage FSC and matching band not worse
than the current ladder, final maps equal by eye, and the NU-stage band
not running ahead of the previous ladder's on an inflated FSC.

## 7. Later steps, explicitly not in this change

- Soft synthesis: apply the global FSC's Wiener profile translated so its
  0.5 crossing sits at the selected rung, instead of the hard cutoff;
  selection stays on hard candidates. With soft references the handoff could
  follow the anchored profile's 0.143 point (`k(v) + (k_0.143 − k_0.5)` of
  the global curve) — the RELION `incr_size` analogue. One experiment each,
  after §6 passes.

## 8. Decisions (Hans, 2026-09-16) and what remains open

Decided:

1. Fine-rung spacing: constant 2 shells, widened to fit the bank budget.
2. `NU_BANK_FSC_HEADROOM` stays 1.5.
3. Regularized candidate: the finest member, replacing every hard rung at
   or finer than the pair's FSC=0.143 resolution (§3.3).
4. `nu_refine` key removed.
5. `postprocess_nu` uses the same machinery (dense ladder, regularized
   pair as finest member).
6. Bank budget `NU_BANK_MAX_MEMBERS = 16`, met by widening the fine step,
   never by failing (§3.1).
7. abinitio3D uses the same dense ladder and finest-member rule -- one
   implementation, abinitio runs as a regression gate (§2, §6).
8. No band headroom beyond the regularized candidate's FSC=0.143 extent.
9. The shared 1% signal-voxel floor applies to the regularized candidate
   as to every label (§3.5).

Open:

- Whether the "add beside a hard rung" variant of §3.3 is worth an
  experiment, decided on the §5 diagnostics after the first run.
