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

**2026-09-11 -- Final and startup PCG reconstructions get a support
reference.** On bgal the refine3D_auto resolution doc reported the gridding
FSC wording for a PCG run: `bootstrap_rec3D` deleted `vol<state>` from the
final `reconstruct3D` command line, so `build_pcg_state_support` found no
reference, the base pair bootstrapped on the sphere and only the replay
took the envelope; the reported FSC was the post-hoc/corrected one. Fix:
the final PCG `reconstruct3D` receives the step-2 gridding bootstrap map
(same name, same sampling, automasked as the refinement) as `vol<state>`;
refine3D_auto's startup reconstruction receives the initial volume as
`vol1`. Both base pairs are now envelope-constrained and `FSC MODE` reads
estimator-constrained, as in every refinement iteration.

**2026-09-13 -- Shell-walk acceptance is a majority test; matching handoff
floor on signal voxels.** Review of the aldolase refine3D_auto log (38.8k
particles, box 192 at 1.3 A, PCG) with the ladder-only Potts pricing in
place: the cleanup no longer erases the walk (`voxels on walked labels 62333
-> 62333`), which exposed the acceptance rule. The challenge compares filters
one shell apart voxel by voxel with no margin (null win rate 50%) against a 5%
threshold, so every shell was accepted until the frontier halved below 32
voxels: depth = log2(frontier/32) = 12 at the bootstrap, 11 at iteration 1,
with four shells accepted while the majority preferred the coarser filter.
The raw-finest handoff then set the matching band from 54 voxels (3.37 A vs
FSC=0.143 3.62 A) and, once the walk died from iteration 2 on (0% wins on
frontiers of 4-424 voxels), from seeded remnants of 4-36 voxels; the band
sat finer than the FSC in every iteration, the FSC never moved in 20
iterations, orientation overlap was 0.998 from iteration 4 and cFAR decayed
0.70 -> 0.55. Changes: (1) `extend_nu_filter_highres` accepts a shell only if
`wins - n/2 >= NU_HIGHRES_EXTENSION_MAJORITY_Z (3) * sqrt(n)/2` on the tested
frontier, seed floor unchanged (`majority_z` in the stats; logged per
challenge with the win count instead of the previously misleading
`extended 0` on rejections). (2) `get_nu_filtmap_finest_selected_lp` takes
`min_signal_pct`: the floor is relative to the signal voxels (mask minus
`nu_solvent_lmask`, `count_nu_solvent_clamped`), and
`record_nu_alignment_lowpass_limit` passes `NU_ALIGN_LP_MIN_SIGNAL_PCT` (1%)
with `min_assigned_pct=0`, logging `NU MATCHING LOW-PASS HANDOFF` on the
master. On the aldolase bootstrap both rules land at 3.67-3.73 A. Compile and
rerun are Hans's; the expected signature is a walk that stops where the win
rate crosses 50%, a handoff at or just finer than the FSC=0.143, and cFAR no
longer decaying.

**2026-09-13 -- Ordered-label Potts prior prices the discrete ladder only.**
Hans's hypothesis: the prior is too conservative for `nu_refine=yes`. In the
code every accepted shell became one more integer coordinate
(`new_coords(i) = real(i)` at the three bank rebuilds and in
`setup_nu_candidate_coords`), so the hinge `(d-1)+(d-1)^2` priced one
Fourier shell like one ladder rung and grew quadratically with the walk
length, while adjacent walked candidates differ little in unary; the
post-extension cleanup (`refine_nu_extension_filtmap_ordered_labels`, run
after every accepted shell) pulled the leading edge back and the next
frontier shrank under the 32-voxel / 5% gate. Change
(`nu_potts_coord_for_label`, `nu_static_ladder_count` in
`simple_nu_filter_bank.f90`): coordinate = ladder position for labels
1..min(8, n_base), the finest ladder position for every walked shell. The
coordinate is no longer a label identity: `nu_effective_base_label_for_candidate`
returns the clamped label directly, and the evidence-state ordering assertion
checks `cutoff_finds` instead of coordinates (the evidence posterior keeps its
own shell-distance continuation). New master log line per cleanup:
`NU post-extension cleanup: voxels on walked labels a -> b`. Discontinuity
statistics use the same coordinates and now count walked transitions as
identical. Compile and rerun are Hans's; the expected signature is walked
populations surviving the cleanup and `accepted shell steps` advancing across
iterations.

**2026-09-14 -- Regularized pair in closed form; replay solve and its
FSC-shrunk start retired; MRES on the summary line.** Trigger: aldolase
`refine3D_auto` (18.6k ptcls, box 192) with every replay start rejected
(`INIT` 3-7) and every from-zero replay at `RESID` 1.0-1.9 after two
iterations against base `RESID` 0.034 -- the same pattern as every PCG run
since the warm-start retirement. The sidecars in
`~/for_claude/PfCRT_regression` show what it is: (1) the preconditioner and
the kernel operator are both Fourier multipliers, so up to the support crop
the preconditioned operator is the identity for base and replay alike, and
in the preconditioned norm they converge alike (`iter2_preconditioned_resid`
0.18-0.31 base, 0.26-0.45 replay); the L2 residual the summary prints
collapses for the base (the first step fits `b` everywhere, noise shells
included) and stays at `~b_k` for the replay wherever the prior dominates
and the solution is ~0 by design -- it measures unfitted noise, not an
unconverged map. (2) `P_tau = rho_mean/(tau*SSNR)` with the FSC clamped at
0.001 is 1000x the density beyond the band; `ml_prior_to_data_khat_l1` is
400-600 in every abinitio stage and any refinement whose FSC reaches zero
inside Nyquist (7 only in Sep07 stage 10, the one replay that ever looked
converged on the L2 number). (3) The FSC-shrunk start `W_k = FSC_k` is the
`P_tau` optimum only for a shell-isotropic operator; `rho` varies within a
shell by orders of magnitude and the base solution at undersampled voxels is
bounded only by the 1% floor, so `(H+P)Wx_base ~ rho_mean(1-FSC)x_base(v)`
there, tens of times `b(v)`: `INIT >> 1` on every anisotropically sampled
dataset, rejection correct, start pointless (the first CG direction from
zero, `M^-1 b`, is already the voxelwise Wiener-shrunk gridding map).
Decision (Hans): ship the closed form. `reconstructor_pcg%shrink_by_ml_prior`
scales the base solution's padded-lattice coefficients by
`1 - P_tau*precond = (rho+floor)/(rho+floor+P_tau)` voxelwise (window
converted to the solve domain and back), on the replayed operator with the
prior installed -- the diagonal-model optimum, without the global step-length
overshoot two CG iterations from zero produce -- and returns the L2 and
preconditioned residuals of the result against the replay system (one
operator application) as a diagnostic of the support coupling it leaves out.
Both execution paths (`regularize_state_half`, the distributed ML job) take
it; `regularized_ml_initial_guess`, `ml_shrinkage_filter` and
`fsc2shrink_filter` are removed; no production solve starts nonzero. The
summary line gains `MRES=` (final preconditioned relative residual; the last
`z` is now formed at exit so the final entry exists in the sidecar as
`final_rel_resid_m`), `KIND=ml` reads `ITS=0 STOP=closed_form`. The
regularized pair costs the raw re-read plus finalization plus ~4 FFT passes
instead of two CG iterations. What to compare on the aldolase rerun: the
shipped pair against the former 2-iteration replay (the map Hans judged
visually good) and against `fsc2optlp` of the `_unfil` pair, by FSC between
them and by eye; `MRES` of the closed form tells how far the diagonal optimum
sits from the coupled one. Compile and runs are Hans's.

**2026-09-16 -- The NU shell walk retired; one generated ladder with the
regularized pair as its finest member for every workflow.** Trigger: the
`nu_refine=yes` walk never let the ML-regularized volume compete for voxels
(it was the finest STATIC member and the walk challenged hard shells beyond
it), and its acceptance rule -- a full-weight shell helps cross-half
prediction only where the half-map SNR exceeds 1 -- is local FSC>0.5, which
is conservative for filtering the map and no better than the FSC-0.143
extent for setting the matching band; on aldolase it pinned the band while
the actual ceiling was the data model (no per-particle CTF or polishing).
Meanwhile every abinitio3D run used the discrete ladder plus auxiliary
competition and it works. Decision (Hans): delete the walk and densify the
static ladder instead, one implementation for abinitio3D, refine3D_auto,
refine3D and postprocess_nu. The bank is generated per box: coarse rungs
20/15/12/10/8/6 A, then fine rungs every `NU_LADDER_FINE_STEP`=2 Fourier
shells up to one shell coarser than the closed-form regularized pair
(`ml_reg=yes`), which is the finest member and competes voxelwise like any
rung; without a regularized pair the rungs are bounded at `fsc/1.5`
(postprocess_nu, `ml_reg=no`); the fine spacing widens (never truncates)
to fit `NU_BANK_MAX_MEMBERS`=16. Potts coordinates and the evidence
candidate masses are positions on the reference ladder 20-4 A interpolated
in log(1/res), so the smoothing prior is scale-free across boxes. The
matching handoff is the content extent of the finest label holding >=1% of
the signal voxels at it or finer: a rung's cutoff, the pair's FSC=0.143
resolution for the regularized member, no headroom (abinitio3D keeps its
own 4.5 A stage limit for the non-NU FSC=0.5 promotion only). Removed:
`simple_nu_filter_extend.f90`, the `nu_refine` key (parameters, phases, UI,
every workflow default), `refine_nonuniform_filter_bank`, the highres-step
sidecars, `nu_static_ladder_count`, the walk statistics and the majority
test, the `NU BANK CAP/UNCAPPED` log lines (now one `NU BANK` line naming the
rungs, the fine step and the finest member). Design record:
`nu_refine_ml_estimator_competition.md` (rewritten). Compile and the
aldolase/PfCRT A/B reruns are Hans's; the expected signature is the
regularized label populating the core, hard rungs the periphery, and the
handoff tracking the FSC=0.143 extent rather than a shell-walk frontier.

**2026-09-16 -- Final reconstruction: automsk and the lag-one reference
forwarded on the direct route too.** bgal refine3D_auto (native box 256 =
registration box, committed sigmas reused, so `calc_final_rec` took the
direct `reconstruct3D` route rather than `bootstrap_rec3D`): the shipped
PCG map ran on the sphere (`rec_final_state01_pcg_support.txt:
solve_support=sphere`, FSC mode "density envelope applied post hoc")
while every refinement iteration had been estimated on the density
envelope. `prep_final_rec_cline` forwarded `automsk` and `vol<state>` only
on the bootstrap route. Now `automsk` is inherited on both routes and the
direct route sets `filt_mode=none` and passes the last refinement volume
as `vol<state>` on PCG, so `build_pcg_state_support` derives the same
lag-one envelope and the doc reads "halves estimated on a support
envelope" like the iterations. The first-run readings of the retired-walk
bank on this bench: abinitio3D 822 s, stage-6 bank 13 rungs + MLreg (14.6%
of the mask), refine3D_auto 1040 s, registration pass moved 1.8% of
directions, FSC0.5 4.30 -> 4.03 A, MLreg ~6.5% at box 256 with the 3.63 A
rung below it at 0.02%; MRES of the closed form 0.86-1.3 at box 256 versus
~0.1 in the cropped stages. Open: the gridding regularized pair won 0% in
the abinitio final's bootstrap bank (PCG wins 6-15% at the same box), to be
checked on the gridding aux path.

**2026-09-16 -- Regularized solve: maxits_ml coupled iterations from the
closed form.** bgal refine3D_auto at the native box 256: the closed form's
preconditioned residual `MRES` read 0.86-1.3 at every iteration and 0.93 in
the abinitio final, against ~0.1 in the cropped abinitio stages, whose band
reached Nyquist (box 140 at 2.33 A/px). At the native box a third of the
shells lie beyond the band with `P_tau` ~1000x the data, and that is where
the diagonal and coupled answers can differ; whether they differ inside the
band (the map) or only beyond it (bookkeeping) the number cannot say. The
regularized member still won 6-7% of the mask, so nothing was blocked.
Decision (Hans): run 2 coupled iterations of the regularized system FROM
the closed form (`maxits_ml`, default 2, `rtol=0`, no start rejection,
indefinite stop falls back to the closed form). Unlike the retired from-zero
replay, the start already has the prior's spectral shape, so CG can only
move toward the coupled solution; ~2 s per half per iteration at box 256.
Diagnostics: `INIT` on the `KIND=ml` line is the closed form's L2 residual,
`MRES` the final one, and a `CLOSED-FORM MRES a -> b` line reports the FSC
between the start and the solved map (0.5/0.143 crossings and the minimum
over the pair's FSC>0.143 band). Reading rule for the rerun: band minimum
~1 with `MRES` down means the iterations only touched beyond-band content
and the closed form alone was adequate; band minimum well below 1 means the
closed form was leaving in-band signal on the table and the iterations
stay. Same bench rerun (abinitio3D + refine3D_auto, PCG) is Hans's.

**2026-09-16 -- maxits_ml result: refinement-invisible; default 0.** Same
bench rerun (7_refine3D_auto, box 256, `maxits_ml=2`): `MRES` of the
regularized half went 1.02 -> 0.29 at the bootstrap, 0.71 -> 0.49 and
0.81 -> 0.67 in the iterations (two iterations do not reach the coupled
solution either); the FSC between the closed form and the solved map stayed
>= 0.965 on every in-band shell and crossed 0.5 only at 2.7-3.1 A, beyond
the 3.59 A band; the regularized member's share of the mask rose 0.4%; and
every iteration's FSC0.5/0.143, handoff and the final map were identical to
run 6 to three decimals. Decision (Hans): default `maxits_ml=0` (closed
form shipped as before); the parameter and the diagnostics stay as the
knob. The diagonal model is adequate for the regularized member's two jobs
(competing for voxels, setting the band); the residual it leaves lives in
the beyond-band shells.

**2026-09-16 -- PfCRT abinitio3D regression: the merged-reference climb
needs the band ahead of FSC=0.143; headroom restored for nonuniform_lpset.**
`~/for_claude/PfCRT_regression/current_broken` (five restarts): with the
regularized member handing off its FSC=0.143 resolution (the 2026-09-16
band rule) the NU stages crept one shell per iteration (stage 6: 9.41 ->
8.39 A over 30 iterations, bank `4 hard rungs 19.4 -> 10.0 A + regularized
pair at 9.41 A`), stage 8 was carried only by its explicit 6.0 A limit and
the final map stalled at 7-8 A. The healthy run (`latest2`, former static
bank) at the same point: FSC=0.143 8.87 A, bank capped at fsc/1.5 = 5.86
A, the regularized pair in the 5.97 A slot winning 5.4% and the 8 A rung
15%, handoff 5.97 A, then 5.0 and 4.5 A as the FSC caught up to 4.5 A.
The staged climb is a ratchet: both halves align to one merged reference,
so its FSC=0.143 is not a gold-standard band, and matching 1.5x ahead of
it is what pulls the FSC up. bgal, Msp1 and streptavidin climb either way
(high SNR); PfCRT does not. Fix: under `filt_mode=nonuniform_lpset` the
regularized member sits and hands off at `max(fsc/NU_BANK_FSC_HEADROOM,
NU_LPSET_BAND_FLOOR=4.5 A)` (the former bank's cap and finest rung), the
hard rungs run up to that shell; under `nonuniform` (gold-standard
refine3D_auto) it stays at FSC=0.143 with no headroom (validated on bgal
today). One implementation, one mode-dependent resolution for the
regularized member (`nu_aux_effective_resolution`). PfCRT rerun is
Hans's; expected signature: stage-6 handoff ~5.9 A at FSC 8.9 A and the
final at 4.5 A as in `latest2`.

**2026-09-17 -- PfCRT `latest3` (10/10 abinitio3D restarts converge):
map quality gap is in refine3D_auto, not abinitio3D; band floor 4 A;
registration pass gated.** abinitio3D with the restored lpset headroom
reaches its band floor in stage 8 on every restart (FSC=0.143 4.44-4.50 A
at the 4.5 A floor; final FSC=0.5 4.70 A) versus `Jul02_very_good`
(band 4.14 A from the July ladder's 4 A rung, FSC=0.5 4.63 A) -- a small,
real cost of the 4.5 A rung, so `NU_LPSET_BAND_FLOOR` is the reference
ladder's 4 A rung (the 2026-09-08 record already noted 4.5 A pinning the
handoff that asked for 4.14 A). The large gap is refine3D_auto: July's
started iteration 1 at 4.09 A from the abinitio poses and reached
3.61/4.03 A by iteration 3; `latest3` ran the registration pass at the
startup pair's FSC=0.8 band, 8.87 A, which moved 33% of the directions
beyond the basin width (17% beyond twice it), dropped iteration 1 to
5.54/7.96 A, and the four-iteration budget ended it at 4.31/6.61 A with
B-factor -137. At that band a small asymmetric membrane protein's
orientations are not discriminable; bgal (large, D2) reassigned 1.7%.
Decision (Hans): keep the pass unconditionally -- the FSC=0.8 band is
too conservative, not the pass itself -- and set `regpass_fsc=0.5`
(the PfCRT startup pair's FSC=0.5 was 7.96 A against 8.87 A at 0.8),
then, agreeing that one shell finer changes nothing for a particle whose
low-resolution band is the micelle, `regpass_fsc=0.143`: the pass runs
at the working band, where a global search confirms a good registration
and re-basins only the misregistered particles (cost: a full-band
`prob_tab` over 5000 directions, 2-3x the coarse pass). A discard gate on
the reassigned fraction was drafted and withdrawn; instead (Hans) the pass
is `refine=greedy`: exhaustive argmax with the incumbent's direction among
the candidates, so a particle moves only to a better pose, where `prob`
samples the assignment from the table and re-basins at random wherever it
is flat. Rerun of `latest3`'s refine3D_auto is Hans's; the numbers to
read are the reassigned fraction at 4.14 A and iteration 1's FSC against
July's 4.09 A.

**2026-09-17 -- automsk=nu review (Hans's implementation).** Kept: the
opt-in mode, the lag-one NU-evidence envelope as the gridding post-hoc FSC
mask and as the reference mask (assembly `_nu_filt` products, matcher
fallback), density fallback whenever the artifact is absent, incompatible,
empty or has an invalid null, `envfsc=yes` derived from any active
automsk, the `support_kind`/`mask_kind` plumbing and `FSC MODE` reading
`backend= support= posthoc_mask= phase_randomization=`. Changed at review:
(1) the NU envelope is never the PCG solve support -- under nu the density
envelope constrains both solves as under yes, because the evidence null of
a constrained pair is designated on the density envelope's dilation ring
at full support weight, outside any evidence envelope, so a NU-supported
solve empties its own null, invalidates the next envelope and falls back to
a density envelope derived from a map that is zero outside the NU support:
an oscillating support with no way back for excluded density; (2) the
bootstrap final-reconstruction route (box change) kept setting
`filt_mode=none` for the shipped map while inheriting `automsk=nu`, which
the new validation rejects -- it now keeps the caller's filt_mode under nu
like the direct route; (3) the 2026-09-02 PfCRT collapse rationale, removed
from the matcher comment and the rewritten automasking policy, is back in
both: nu is not a default and must not become one.

**2026-09-17 -- Potts coordinate back to the label index; PfCRT map-quality
gap traced to the abinitio3D poses.** (1) `postprocess_nu` on the Sep-11
PfCRT halves produced noise with the retired-walk code. Cause found in the
code, log pending: the 2026-09-16 Potts coordinate (position on the
reference ladder in log(1/resolution)) puts the generated fine rungs about a
tenth of a rank apart while `NU_LABEL_SMOOTH_STEP_TOL=1` frees jumps of one
rank, so the smoothing prior was inert across the whole fine end -- the
label field could alternate freely between rungs 4.6-4.1 A voxel by voxel.
In refinement the regularized member wins the core and the FSC hid it; in
the evidence state, which has no regularized member, the fine cutoffs
scatter and the per-voxel sharpening passes noise. Coordinates are the
label index again in both the filter and the evidence competition; the
log-resolution geometry survives as the evidence candidate mass only.
(2) refine3D_auto on `latest3` (greedy pass at 4.44 A: 6.2% reassigned,
no FSC change) plateaus at 4.03/4.50 A from iteration 1 with B -90..-100,
against the Sep-11 run's 3.93-3.98/4.14-4.25 A with B -73. The gap is
present in the BOOTSTRAP, before any refinement: gold-standard
reconstruction of the Sep-17 abinitio poses reads 4.14/4.50 A, of the
Sep-11 abinitio poses 3.93/4.09 A, from the same particles. The abinitio
poses are worse although the abinitio's own merged-reference FSC reads the
same (4.5 A band in both). Candidates, both introduced between the runs:
the inert Potts prior above (the NU-stage references' hard-rung pattern),
and the closed-form regularized pair (MRES 1.3 on PfCRT at box 300; the
Sep-11 references' core came from the two-iteration replay). A/B after the
recompile: rerun abinitio3D + refine3D_auto with the index coordinate; if
still short, `maxits_ml=2` isolates the closed form.

**2026-09-17 -- postprocess_nu restored on the Sep-11 PfCRT halves with the
index Potts coordinate.** Even/odd run at box 300: 11 distinct evidenced
cutoffs with coherent populations (4.33 A 181k, 3.80 A 44k, 3.38/3.05 A
3-5k voxels), four evidence bands, B -49; "PfCRT back to being gorgeous"
(Hans). Confirms the 09-16 log-resolution coordinate as the cause of the
noise (and the prime suspect for the abinitio pose gap; A/B pending). Two
fixes on the way: postprocess_nu's own half-map FSC was computed on
real-space volumes (`image%fsc` reads Fourier coefficients) and read 0.000
A, so the bank ran unbounded to Nyquist -- harmless here (rungs beyond 3 A
won nothing) but wrong; now transformed before the FSC, so the bank is
bounded at fsc/1.5 as designed. And identical inputs (vol1 = vol2, an FSC
of 1 everywhere, the finest cutoff awarded on 1.36M voxels, a sharpened
noise ball) are refused with a message.

**2026-09-18 -- Index Potts coordinate reverted for the filter competition
(PfCRT `latest4`).** With the label index, 3/9 abinitio3D restarts reached
4.6 A (4.85-8.2 A otherwise) against 10/10 at 4.4-4.5 A with the
log-resolution coordinate (`latest3`). Mechanism in the logs: at FSC 8.2 A
the regularized member sits at 5.45 A but the handoff reads `7.96 A (raw
finest label 5.45 A)` for the rest of the stage -- nothing finer than the
7.96 A rung holds 1% of the signal voxels, because the dense ladder puts
four to eight members between 8 A and the regularized pair and on the
index coordinate the core's jump out of the 8 A surround costs several
hinge steps instead of the one it cost on the 8-rung ladder; on PfCRT's
small core (12% envelope occupancy) the prior wins. The log-resolution
coordinate is the scale-free generalization of the old prices and is
validated by latest3, bgal, Msp1 and streptavidin; restored for the filter
competition. The evidence competition keeps the index coordinate
(separate Potts problem, separate beta; validated by the "gorgeous"
postprocess_nu on the Sep-11 halves). Correction of the 09-17 entry: the
first "ball of noise" postprocess_nu run had vol1 = vol2, so the
log-resolution coordinate was never shown to produce a noisy
postprocess_nu with correct inputs; the index change was made on a
misdiagnosis. Reruns are Hans's: abinitio3D (expect latest3 behaviour)
and postprocess_nu on the Sep-11 halves (expect unchanged).

**2026-09-18 -- NU machinery restored to the ed36eb4c static ladder + aux
competition, as the only mechanism.** Decision (Hans): the discrete ladder
`[20,15,12,10,8,6,5,4]` A capped at `fsc/1.5` with the ML-regularized pair
replacing the finest retained label -- abinitio3D's machinery at commit
ed36eb4c, which produced the best PfCRT maps -- is the NU competition for
abinitio3D, refine3D_auto and postprocess_nu; nothing else. The
`src/main/nu_filt` submodules and `simple_nu_state_filter` are restored
from ed36eb4c verbatim except: the extend submodule and its interfaces,
the walk statistics type and the walk-only constants are removed (no code
path can extend the bank), the `nu_refine` control stays removed, the
auxiliary member is used whenever `ml_reg=yes`, and Hans's `automsk=nu`
edits (envelope selection, `write_nu_evidence_envmask(mask_out, l_valid)`)
are kept. Gone with this: the generated dense ladder, `NU_BANK_MAX_MEMBERS`,
`NU_LADDER_*`, the log-resolution and index Potts coordinates (integer
ladder coordinates again), the lpset headroom/floor (`NU_LPSET_BAND_FLOOR`;
the cap provides the headroom as it always did), the content-extent
handoff. postprocess_nu passes the half-map FSC so its bank is capped like
the refinement's. Compile and the PfCRT/bgal reruns are Hans's.

**2026-09-18 -- The auxiliary pair competes with the finest rung instead of
replacing it.** Decision (Hans): with `ml_reg=yes` the ML-regularized pair
is appended as one more bank member beside the finest retained rung, and
the two compete voxel by voxel; there must be no prior penalty for
replacing a finest-rung voxel by the regularized pair -- if it wins the
unary it is included. Coupling to the ordered-label prior: the auxiliary
shares the finest rung's Potts coordinate (a boundary between two filter
shapes of the same terminal resolution is not a resolution discontinuity),
so every existing price is unchanged and the unary alone decides between
them. Side effects: the auxiliary is always included when supplied (the
ed36eb4c rule ignored it unless finer than the finest rung); the handoff
sorts labels by resolution so it needs no change, and if the auxiliary
falls below the 1% signal-voxel floor the band now drops to the finest
rung rather than to the one below it; the auxiliary's own unary margin
diagnostic now includes the finest rung. `setup_nu_dmats` appends the
label (own Fourier index, bwfilters column unused, filtered pair never
cached), `setup_nu_candidate_coords` assigns the shared coordinate.
Untested; the PfCRT abinitio3D + refine3D_auto pair against July's
3.61/3.98 A is the test.

**2026-09-18 -- postprocess_nu operates on the project.** Decision (Hans):
like `postprocess`, `postprocess_nu` takes `projfile` (+ `state`, `mskdiam`,
`nthr`, optional `outvol`), fetches the state volume from the out segment,
its `_even_unfil/_odd_unfil` pair (evidence input) and, when present, its
`_even/_odd` regularized pair, with which it first runs the refinement's
filter competition -- static ladder capped at fsc/1.5, the regularized pair
beside the finest rung -- printing the bank, the assignment table and the
matching handoff and writing the references and local-resolution map;
then the unchanged evidence sharpening. Every product carries
`_pproc_nu` (`<vol>_pproc_nu.mrc` sharpened; `<vol>_ref_pproc_nu.mrc`,
`<vol>_even/_odd_ref_pproc_nu.mrc`, `<vol>_locres_pproc_nu.mrc`). First
run on the Sep-11 PfCRT abinitio3D project: FSC 6.61/3.98 A, bank capped
at 2.65 A (all 8 rungs), regularized pair beside the 4 A rung (both on
shell 78, 3.981 A): rung 6.6% of the mask, regularized pair 3.9%, handoff
3.98 A; "map looks great" (Hans). The
`vol1/vol2` file inputs are gone. This is the cheap test of the auxiliary
competition on any completed run directory.

**2026-09-18 -- solvent prior as an opt-in SOFT real-space ridge on the replay.**
Decision (Hans): "there is solvent left to flatten, no doubt about that",
but "the prior needs to be soft". A first cut refined the replay's hard
support with a Wang-type mask; with the default `maxits_ml=0` that is the
closed-form shrink multiplied by a solvent mask, i.e. exactly the
post-solve masking the backend forbids, and with iterations it is still a
hard cut at a single Otsu boundary. Discarded. What ships:
`pcg_solvent=yes|no` (default `no`; refine3D/refine3D_auto/abinitio3D/
reconstruct3D; requires `rec_backend=pcg`; forwarded like the other
backend keys, stripped for the gridding children) with
`pcg_solvent_lambda` (default 1.0, relative to `data_scale`). On: the
regularized operator gains `lambda_s (1 - w(r))` on its diagonal, `w` the
logistic protein weight built PER HALF from that half's own current base
map (half-independent prior, regularized pair stays gold standard) at
`max(8 A, 2 x res0143)` (Otsu inside the production support, width = the
solvent class's spread); `maxits_ml` defaults to 2 with the prior on
(parameters class, `PCG_SOLVENT_MAXITS_ML`; explicit 0 refused) because
the closed form cannot see a real-space term. Weight volumes written per
half (`pcg_solvent_weight_stateNN_even|odd.mrc`, overwritten every
iteration, never accumulated), even/odd correlation and solvent-fraction
gap logged (the half-independence check). In abinitio3D the key is withheld
through stages 3-7 and passed only to stage 8; shortened workflows never use
the prior. Support, base pair, FSC oracle and NU bank inputs are untouched on
purpose: the resolution claim carries no extra
mask and the prior's effect is confined to the shipped regularized pair.
No low-pass guard on `w`: a ridge modulator, unlike a multiplicative
mask, injects no spectrum, and the statistic is band-limited at 2 x res
anyway. Runs under `automsk=no` too (sphere as support). Provenance
`solvent_prior=soft lambda_rel=<x>`. Off: bit-identical.
Untested; first test is the PfCRT cold-start run with `pcg_solvent=yes`
against its 3.98 A baseline, judged on the map (cavities, detergent
belt, skirt) rather than the FSC, then a `pcg_solvent_lambda` sweep
(0.3, 1, 3).

**2026-09-19 -- solvent prior moved to the base solve (per-half re-solve).**
Observation (Hans, several data sets): stage maps under the replay-ridge
design looked less noise-fitted but the final map looked too low-pass
filtered. Diagnosis: the on/off comparison was confounded with
`maxits_ml` 2 vs 0, and the coupled replay's proper P_tau attenuation
plus the classical postprocess's `2FSC/(1+FSC)` double-attenuates a
regularized map (the closed form under-attenuates, hiding it; the stage
maps ship NU-filtered base rungs, hiding it). Decision (Hans): keep runs
comparable -- a separate cold solve with the prior on the base system,
its pair fed to the NU machinery, `maxits_ml=0` and the closed-form
shrink as standard. What ships: prior-free base pair (operators kept)
-> its FSC sets the smoothing scale, per-half weights from the
prior-free halves -> ridge installed -> both halves solved again from
zero with the same budget on the same accumulators -> that pair is the
base pair everywhere (FSC, NU evidence, replay, B factor). NU caveat
recorded, not acted on: the NU null (dilation ring / robust bulk) is
measured on the prior'd pair; watch `NU NULL SHELL GEOMETRY`, the null
bias median/MAD and the assignment table for collapse. Summary lines
`KIND=pre` (prior-free) then `KIND=base` (prior'd). `postprocess` gained
`fsc_filt=no` (skip the optimal filter, low-pass at FSC=0.143 only) and
`imgkind=unfil` (classical route from the unfiltered pair average) for
the re-reconstruct/re-postprocess matrix on the completed runs.

**2026-09-19 -- abinitio3D robustness: one bank rule for every workflow;
the matching band is the finest member of the bank.** Log analysis of every
PfCRT set (record: `claude/pfcrt_abinitio_band_lead_analysis.md`; note that
the `MATCHING LOW-PASS LIMIT (THIS ITERATION)` readout is the handoff for
the NEXT iteration): with the regularized member handing off its FSC=0.143
(b9cc8e31 `current_broken` 0/4; 4ff3a1f6 `latest5` 1/9, the ed36eb4c
restoration of 09-18 had dropped the 09-16 headroom placement) the band
equals the FSC in 41-42% of stage-6 iterations, the crossing advances
about one shell per iteration and the assignment crawls or freezes; with
the member at `max(fsc/1.5, 4.5)` (6a41a40f `Sep16_10of10`) the band was
pinned in 4% of iterations, +0.47 A per iteration, 10/10 converged, and
the 4.5 A floor held stages 7-8 below their crop Nyquist (4.14/3.88 A).
July's uncapped rung handoff led by 1.55x median (p90 1.88) and froze 2/10
in stage 6. `latest4` (281ca643, index Potts coordinate): the member one
prior step from the finest rung won <= 0.6% of voxels, the handoff fell
back to the rungs, 3/9. Replacement versus competition is not the
separator (replacement at FSC=0.143 failed identically, iteration 86 of
`current_broken`: MLreg 3.9%, band 8.62 = FSC). Stages 1-5 end alike in
Sep16 and latest5 and better than July's. Hans's reading, adopted: the
band at the FSC starves the alignment (underfitting); the objective's
sigma2 weighting, the stochastic assignment and the evidence-limited NU
reference guard against overfitting, so the 09-16 "inflated FSC" story is
withdrawn. Rule (Hans): the ladder cut at fsc/1.5; the ML-regularized pair
carries its FSC=0.143 and joins beside the finest rung the moment that is
at or beyond the ladder's finest rung at the box (`setup_nu_dmats`; within
the ladder the rungs compete alone, since the cut always keeps a rung
between fsc/1.5 and fsc); the finest member of the bank is the handoff
(`get_nu_filter_bank_finest_lp`, `record_nu_alignment_lowpass_limit` and
postprocess_nu), the 1% signal-voxel floor and the raw finest label are
diagnostics on the handoff line. Same code path for abinitio3D,
refine3D_auto, refine3D_states, classify3D_refs and postprocess_nu; the
interim lpset-only cap placement of the same morning was superseded before
any run. Band per FSC on the PfCRT boxes: 10.7 -> 7.96, 8.9 -> 5.97,
7.6 -> 5.97, 6.9 -> 5.01, 6.0 -> 3.98 (or the crop Nyquist 4.44/4.14),
5.0 -> 3.98, 4.03 -> 3.98, 3.93 -> 3.93 (pair admitted). Expected signature
(Hans's rerun): handoff 5.97 A at an 8.9 A FSC in the first stage-6
iteration, `NU AUXILIARY MEMBER ... within the ladder ... the rungs compete
alone` throughout the lpset stages on PfCRT, FSC=0.143 at the stage-6
Nyquist within ~10 iterations, overlap < 0.5 through stage 6, stages 7-8
climbing to 4.1/3.9 A. refine3D_auto on PfCRT: the pair admitted only once
the FSC reaches 3.98 A. Not changed, on record: the stage early-stop
(overlap 0.9/0.95) cannot tell a frozen assignment from convergence.

**2026-09-21 -- solvent prior: estimate on the base pair, apply to the
prior'd pair; isotropic postprocess protocol; postprocess_nu names.**
`~/for_claude/pcg_solvent=yes` (exp_gate, msp1 at ae97341c; the
streptavidin logs there are ec5527df, pre-prior, not comparable): the
`_lp` maps the best yet, the `_pproc` maps noise. Measured: the prior'd
pair was the `_unfil` pair, its FSC (3.64/3.43 A against 3.88/3.91 A
prior-free; even/odd weight correlation 0.971/0.996, no randomization,
stage-8 crossing pinned at the crop Nyquist in every iteration) set the
optimal filter, which stays open to FSC=0.05, and the B-factor
(-121/-115, density-windowed on top of the ridge) sharpened the
closed-form map with nothing closing it. postprocess_nu on that pair
"looks like madness": its evidence null and the NU whitening (the MAD of
even minus odd per radial shell over a support that is 80-93% solvent)
are measured where the prior removed the noise. Decisions (Hans): the
FSC is never computed from a solvent-prior'd pair; the prior must
influence every reference voxel, not only the band beyond the ladder, so
the NU label field is derived on the prior-free base pair and applied to
the prior'd pair (`nonuniform_filter_state` apply pair -> `nu_filter_vols`
per-label composition of the apply halves); the prior-free pair is the
base pair everywhere (FSC, cap, band, competition, evidence, `_unfil`,
start-volume handoff), the prior'd pair is `_even_solvent/_odd_solvent`,
the base of the replay and of the references (both PCG paths:
`reduce_solve_state_pair` solvent outputs, `resolve_base_pair_with_solvent_prior`
into the solvent pair); `automsk=nu` needs no special case (its null is
the base pair's). Postprocess (classical): cutoff = FSC=0.143 of the
reconstruction's FSC file (now always the base pair's curve; the `_unfil`
pair's own FSC computed there only when no file is given), B
from the map being sharpened between `HPLIM_GUINIER` and the cutoff,
sharpen, Butterworth at the cutoff; `fsc_filt` and the density window
retired; `imgkind=unfil|solvent` postprocess the pair averages with the
same cutoff (`pair_stem`). postprocess_nu: evidence from `_unfil`,
sharpening applied to the `_solvent` pair when present; outputs named
like postprocess with `_nu` in the suffix (`_pproc_nu`, `_pproc_nu_mirr`,
`_locres_nu`; the `_ref_pproc_nu` references are gone). Also on record:
the final PCG map is the 5-iteration cold solve; the gridding
reconstruction before it is bootstrap_rec3D's sigma seed, the only NU
competition in the final block (its gridding ML pair beyond the ladder
wins 0% of voxels: open since 09-18, inconsequential for the map). To
watch on the rerun: the B-factor of a solvent-flattened map, and whether
the prior'd references keep what made the exp_gate/msp1 `_lp` maps.
