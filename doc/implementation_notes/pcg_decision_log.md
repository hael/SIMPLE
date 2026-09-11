# PCG backend: decision log

One dated paragraph per decision, in order: what was decided, what forced
it, what it superseded. The full record behind each entry is in
`pcg_priors_history.md` (cited by its dated headings and, after
2026-08-29, by "dev item N" of its Active dev list). The present-tense
description of the backend is `pcg_backend_overview.md`; a new decision is
appended here AND the affected overview paragraph is rewritten. Entries
that are no longer in force say so in their last sentence.

---

**2026-08-25 -- Support contract verified; cross-iteration ML warm start
adopted.** Audit confirmed `P H P` on both sides of the normal operator and
one `P` on the right-hand side (`test=pcg_recon` stage 3b asserts symmetry
and positivity). The ML replay was made to warm-start from the previous
iteration's ML half (validated on bgal), because a 2-iteration replay
from the unregularized base could not close the beyond-band gap. History:
Stage 0/1. *Warm starts retired 2026-09-10.*

**2026-08-26 -- Binary-envelope solvent prior `Q_s` on by default at
`pcg_solvent_lambda_rel=0.1`; abinitio3D forces `automsk=yes` from the
first NU stage; suppression readout in the convergence report.** A
real-space penalty on solvent variation (weighted centering, `L^T L`),
anchored to the data scale, calibrated on bgal (active regime ~1e-1,
over-flattening from 3e-1) and the neutral fixture. History: §4, the
DECISION entries of 2026-08-26, Stage 5 records. *Prior removed
2026-08-27; the forced automasking removed 2026-08-28; the readout left
with the prior.*

**2026-08-27 -- Suppression thresholds are not portable; over-flattening
guidance moved to shipped-vs-base FSC inflation.** The streptavidin and
bgal ladders did not overlay (convergence state and data-scale anchoring
differ across box sizes), so the absolute-percentage bounds were retracted
and the inflation of the shipped pair's FSC=0.143 over the base pair's
(> 5%) became the portable over-regularization signal. History: Stage 5
records of 2026-08-27. *Left with the priors.*

**2026-08-27 -- Direct NU-evidence precision `Q_NU` becomes priority 1;
`Q_s` removed.** The NU competition's per-voxel cross-half evidence,
compacted into a frozen evidence state with an explicit calibrated null,
replaces `P_tau` in the replay in NU mode (mode-exclusive; never
additive). The first run failed on a saturated null (a detection
threshold used as a likelihood offset); the same-day redesign to a
lower-quartile null center passed, and the first truth-judged 1WCM
two-way ablation showed `Q_NU` beating `P_tau`. `Q_s` and its plumbing
were removed. History: §3, §5, Stage 6 records. *`Q_NU` removed
2026-09-06.*

**2026-08-28 -- `Q_NU` default-on in NU mode (`pcg_nu_lambda_rel=0.1`
dynamic default); post-hoc NU filtering and forced automasking retired
from the PCG NU path; `Q_NU` supports trailing reconstruction (evidence
pair = FSC pair).** First NU-replay abinitio3D on bgal: the loop worked,
the solvent-era stage-6 overlap collapse was largely cured, cost was
production-viable. History: dated DECISION entries of 2026-08-28.
*Superseded 2026-09-06.*

**2026-08-29 -- Wilson prior `Q_W` implemented, corrected twice on
truth-judged failures, validated, adjudicated against `Q_NU` and REMOVED
the same day.** On 1WCM `Q_NU` dominated at every shell even with `Q_W`
given the ground-truth spectrum. Stage 7 closed. History: Stage 7.1
records and the removal DECISION.

**2026-08-29 -- `automsk=yes` regenerates the evidence envelope under the
`Q_NU` replay; refine3D_auto joins the PCG bypass; adaptive band
granularity (Stage 6.6) reduced to the `nu_refine` shell walk.** The
a-priori frontier-tracked candidate proposal degraded on 1WCM and was
withdrawn; the final form mirrors the gridding challenger (accept a shell
only on unary win-fraction >= 5%), validated on 1WCM. History: DECISIONs
of 2026-08-29, Stage 6.6 records. *The evidence envelope and the
`nu_refine` walk survive; the `Q_NU`-specific parts left with `Q_NU`.*

**2026-08-29/30 -- `Q_NU` inert on PfCRT at the default strength;
auto-lambda controller implemented.** Suppression 3% on PfCRT versus ~35%
on 1WCM and bgal at the same `lambda_rel`: the data-scale anchoring uses
the low-band mean of `D` while the prior competes in the fine shells, so
the effective strength is a dataset property. A one-pole plant model and
a memoryless secant controller on the persisted suppression readout were
added, then (2026-08-31) an AIMD outer loop on the setpoint driven by the
shipped-pair FSC trajectory. History: dev item 2. *Removed 2026-09-06.*

**2026-08-31 -- Final reconstructions keep the PCG backend.** The
original-sampling final reconstructions of abinitio3D and refine3D_auto
had silently dropped to a gridding `objfun=cc` pass; they now run
`bootstrap_rec3D` on the refinement's backend (cc pass for sigma
availability, then the regularized pass). History: dev item 2
(FINAL-RECONSTRUCTION Q_NU POLICY).

**2026-09-01 -- `envfsc=yes` restored as the refine3D_auto default;
solvent constraint via the lowest-bin convention (dev item 4); direct PCG
support constraint `pcg_mskfile` (dev item 5, experimental); refine3D_auto
startup made reconstruct -> mask -> re-reconstruct -> refine (dev item 7).**
The conservative density automask at `envmsklp` returned as the FSC
envelope; the NU filter and prior regained a solvent clamp (background
voxels fixed to the coarsest candidate) without re-multiplying references;
an arbitrary [0,1] support can be installed as `P`; and the PfCRT
refine3D_auto collapse was traced by elimination to a first matching
against unfiltered raw halves. History: dev items 4, 5, 7. *The
"references multiplied by the NU envelope" hypothesis (dev item 6) was
falsified the same day.*

**2026-09-02 -- References are never multiplied by an envelope;
`automsk=yes` means a heavily low-passed background instead.** With the
startup fixed, the PfCRT refine3D_auto collapse tracked the first
consumption of envelope-masked references on both backends
(`rec_backend=gridding` reproduced it; `filt_mode=none automsk=no` was
healthy). The NU evidence envelope omits the detergent by design, so a
reference multiplied by it is protein-only while the images contain the
micelle: a model-data mismatch under the euclid objective. Principle
(user): a reference must never hard-remove density present in the images;
excluded density is down-weighted by fixing the filter-field background
(the envelope complement) to the coarsest bank candidate, cisTEM-style.
The code-review response of the same day removed the post-hoc
`nu_filter_vols` from the PCG path in every mode. History: dev items 8, 9.

**2026-09-06 -- `Q_NU` REMOVED from the PCG backend in its entirety; the
PCG path follows the gridding NU competition exactly.** Together with the
auto-lambda/auto-target controllers, the stats file, the bootstrap
calibration pass, the parameters, UI, readout, operator and
`test=pcg_priors`. Rationale (user): an elegant estimator that ends in a
parameter optimization the competition never needed; the belt on PfCRT
and the msp1 trajectories are the record. The gridding volassemble NU
section was extracted into `simple_nu_state_filter` and is called from
both backends: base pair seeds the bank, `P_tau` pair is the auxiliary
member, `nu_refine=yes` runs the shell walk, finest selected label is the
handoff. History: dev item 2 (DECISION 2026-09-06).

**2026-09-06 -- Solve-support policy correction: no density envelope
anywhere unless `automsk=yes`.** Dev item 5 had made the envelope the
replay support unconditionally in ML mode; every msp1 and streptavidin
run since 2026-09-01 solved the shipped pair on the envelope while the
base pair stayed spherical, so their two FSCs were never comparable. Rule:
`automsk=yes` gates the envelope; under it the replay takes the envelope
and `envfsc=yes` extends it to the base. History: dev item 2
(SOLVE-SUPPORT POLICY CORRECTION). *The envfsc split superseded
2026-09-09.*

**2026-09-06 -- Matching low-pass ceiling in the NU stages; sigma2 loader
carry-over removed; stage-boundary FSC=0.5 promotion.** The fine NU handoff
must not set the abinitio3D matching band by itself (user verdict); the
implicit sibling/parent-directory sigma2 carry-over that had seeded ten
msp1 runs from a foreign run's residuals was removed entirely. History: dev
item 2 (USER VERDICT, LOG EVIDENCE, msp1 STAGE-7 COLLAPSE ROOT CAUSE).

**2026-09-07 -- PfCRT regression (7/7 at 6.0-6.5 A): the ceiling was the
ladder's FINAL value, not the 4.5 A hard bound; NU stages run their full
budget.** Map-level audit by cross-FSC against the July reference showed
the sets parting at stage 6. Fixes: cap = `LPSTOP_BOUNDS(1)`, `minits =
maxits` in NU stages. Not implicated: the backend, the stage-1 floor, the
FSC=0.5 promotion, the sigma bootstrap. History: dev item 2 (PfCRT
REGRESSION 2026-09-07).

**2026-09-07 -- Final-reconstruction stage made standalone
(`bootstrap_rec3D`); canonical sigma2 store reviewed; sequence validated
on msp1 (best msp1 map to date).** Image-power seed, gridding ML bootstrap
map carrying the workflow's filtering, residual sigma pass, PCG final map
with a cold-solve floor of five iterations. History: dev items 10-12.

**2026-09-08 -- `nu_input=gridding|ml` alternatives tried and reverted.**
Feeding the NU competition a gridding-equivalent `M^-1 b` pair made the
PfCRT histograms worse and stage 7 collapsed; the truncated-CG base pair
remains the NU input. History: run record 2026-09-08b.

**2026-09-09 -- One support, applied once: gridding products deapodized
and given the same soft spherical support as PCG; no second mask anywhere;
`evaluate_halfmap_pair` masks nothing of its own; postprocess skips its
post-hoc mask for any volume with the provenance sidecar; hard solve
domain + output window replaces the soft `P H P`.** Audit of every
`mask3D_soft` consumer; the soft formulation had let the band depend on
the CG state where `0 < P < 1`. History: records 2026-09-09b/c.

**2026-09-09 -- `automsk=yes` implies `envfsc=yes` on both backends; the
density envelope constrains both PCG solves.** The 09-06 split (spherical
base under `envfsc=no`) is retired: on PCG the envelope is the support of
base and replay once a prior reconstruction exists, so the FSC pair is
envelope-constrained in the estimator; on gridding the same envelope is
applied post hoc with the phase-randomized correction. The masked-FSC
bias of the constrained pair is acknowledged and REPORTED (`>>> FSC MODE`),
not corrected. Final reconstructions inherit `automsk`. History: dev item
2 (SOLVE-SUPPORT POLICY, ENVFSC COUPLING; CODE REVIEW RESPONSE).

**2026-09-09 -- NU evidence envelope keeps control of the filtering; its
null has two regimes; shared physical dilation minimum.** With the PCG
base pair envelope-constrained, the envelope's mixture null (median/MAD
over the whole support) sat on the zero/zero spike; the null is now
designated on the density envelope's dilation ring for a constrained pair
(Euclidean shell, at full weight of the base support) and stays the
robust mixture estimate for a spherical pair, with the density envelope
as the fallback background. The dilation has a shared minimum of 7.5 A
(the former abinitio3D 7 layers at 1.075 A), an explicit `binwidth` wins
either way, and the ring's evidence occupancy is logged. `pcg_mskfile` is
reported as a constrained support. The vestigial `AMSK_FREQ` regeneration
planner was removed. History: dev item 2 (NU EVIDENCE ENVELOPE UNDER THE
DOUBLE SUPPORT; CODE REVIEW RESPONSE).

**2026-09-10 -- Cross-iteration warm starts RETIRED; base from zero, replay
from the shell-shrunk current base; `maxits_pcg=2` kept.** PfCRT
gridding_vs_pcg (10 vs 9 restarts): every PCG solve was ITS=2 fixed; the
replay's relative residual never fell below 1 and grew to 10^3-10^4 in the
failed runs through the warm start from the previous ML half; the base
pair carried no fine-shell evidence at the stage-6 take-off (0% of the
mask at 7.96 A vs gridding's 6.4%). The user's simulated-data calibration
(two iterations beat gridding, beyond five nothing interpretable changes)
fixes the budget. The first cold-start run: 3.98 A with side chains,
identical to the gridding successes, cost unchanged. `INIT=` added to the
solve summary line. History: dev item 2 (PfCRT REGRESSION gridding_vs_pcg
2026-09-10, RESULT).

**2026-09-10 -- Replay start rejected when worse than zero (streptavidin
canonical 9/10).** With the canonical sigma store the PCG backend missed one
streptavidin run in ten (gridding 10/10; legacy 10/10 on both): the map
had the wrong symmetry axis. The logs of the good and the bad run are
identical in the sigma flow (no PCG-specific branch exists in the
canonical store) and part in stage 3, the first PCG/ML stage and the
symmetry-search stage, where the bad run stalled at 6 A instead of
reaching 4.5 A, so the axis was searched on a weaker map. The one
PCG-specific fragility there: the replay started from the shell-shrunk
base at `INIT` 2-4 and ended at `RESID` 0.4-2, a half-converged transient
that is both the matching reference and the symmetry-search input at
that stage. Guard, inside the solver and free: `solve_core` computes the
start's residual anyway, and when it exceeds 1
(`PCG_START_MAX_REL_RESID`) the start is discarded before the first
iteration and the residual is `b` itself (no operator application; a
first version that re-ran the whole solve was withdrawn as too costly).
Logged per solve. Whether this closes the 9/10 is to be measured; the
store itself is not implicated. History: this log entry.

**2026-09-10 -- Documentation split.** `pcg_priors.md` frozen as
`pcg_priors_history.md`; this log and `pcg_backend_overview.md` created.

**2026-09-11 -- FSC gate removed from the NU bank and the matching low-pass
under `nu_refine=yes`; extension is NU evidence only.** refine3D_auto had
stopped extending in resolution. The 2026-09-08 `fsc/1.5` cap
(`init_nu_filter(fsc_res=)`, `nu_bank_cap_find`) truncated the static
bank and bounded the shell walk (`max_find`) by the base pair's FSC=0.143,
so the finest label -- and with it the matching low-pass handoff -- could
never lead the global FSC by more than about two ladder labels, while the
FSC of the base pair lags the local evidence precisely where the map is
extending; the loop that used to carry refine3D_auto from the FSC to the
evidence-supported resolution (July 2026) was closed. Decision (Hans):
`nonuniform_filter_state` passes no FSC to the bank under `nu_refine=yes`
(full static ladder), the shell walk runs without `max_find`, and the
handoff stays the raw finest selected label; `get_nu_bank_cap_find`
removed. The cap remains in static-bank mode (`nu_refine=no`, the
abinitio3D ladder), where the walk does not run and the PfCRT
regularized-pair pinning it was written for was in any case fixed at the
root by the like-for-like coarse-to-fine selection of the same day. Logged
as `NU BANK UNCAPPED (nu_refine=yes)`. Stale policy text describing a
PCG-only two-shell FSC bound on the walk and a 5%-support handoff with FSC
headroom was corrected at the same time (`nonuniform_filtering_policy.md`
sections 8, 10, 12; `reconstruct3D_pcg_policy.md`;
`refine3D_auto_policy.md`). Compile and the refine3D_auto rerun are Hans's.
