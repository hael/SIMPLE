# PCG reconstruct3D → euclid crash: investigation handoff

**Status:** TWO distinct defects identified. Defect A (§13): PCG solves can
leave hot spectral content just above the matching band, exposed when a stage
transition extends `kto`; it is diagnosed in the log (`PCG BEYOND-BAND
EXCESS` line, §14) and the solver-level fix is open (§13.5). Defect B (§15,
found on the 2026-08-21 rerun with the actual run directory in hand, fixed in
§15.3): the PCG solution is a uniform ~×200 above gridding's data-quotient
amplitude convention, so the gridding→PCG backend handoff at stage 3 jumps
the reference amplitudes faster than the euclid sigma2 equilibrium can adapt
— this, not defect A, is the original iter-41/42 crash. Root cause: gridding
divides its output by the original box size (`reconstructor_eo%mag_correction`,
a retained historical normalization) and PCG did not; PCG now applies the
same convention on output (§15.3), identically in both execution paths (§16).
Retiring the convention itself is planned in
`drop_legacy_box_division.md`. §§1–12 are the original handoff, kept
verbatim.

**Date:** 2026-08-20
**Workspace:** `/home/elmlundho/src/SIMPLE` (branch `master`)
**Run dir:** `/home/elmlundho/cache/streptav3/6_abinitio3D` (streptavidin, 10335)

---

## 1. Symptom / crash

`abinitio3D` with `rec_backend=pcg` aborts at the start of **stage 3**
(first stage where euclid + ML regularization turn on), inside `prob_align`:

```
ERROR! all probability-table reference distances are invalid;
simple_eul_prob_tab.f90; line: 612
```

Crash occurs at the first or second stage-3 iteration (iter 41 or 42 depending
on run). It is deterministic across repeated runs; the exact iteration and the
reference amplitude vary slightly with stochastic sampling.

Preceding warnings:
```
WARNING: particles with invalid probability-table distances: ~5807
WARNING: particles with flat probability-table reference distances: 1
```

---

## 2. How to reproduce

- Workflow: `abinitio3D` on the streptavidin project, `rec_backend=pcg`.
- Original box 128, cropped box 88, `smpd_crop` 1.5593.
- Crash is at stage 3 (LP 8.6), the first stage using the euclid objective
  (`objfun=euclid`) together with ML regularization (`ML_REG_START_STAGE=3`).
- Gridding backend (`rec_backend=gridding`) does NOT crash — this is
  PCG-specific.
- **PCG + 2 CG iterations is proven to work on other datasets** (user-confirmed),
  so this is a dataset-specific instability, not a universal PCG bug.

---

## 3. Euclid overflow mechanism (confirmed)

The euclid distance uses `v = 1 + crvec/norm`; `euclid_dist_from_crvec` returns
`huge` (→ invalid) when `v > -log(TINY) ≈ 23.03`. The crash fires when *every*
reference distance for some particle is invalid.

`crvec ≈ ft_ctf2 · ft_ref2` where `ft_ref2 = |pfts_refs|²`. In `v`, the per-shell
`sigma2` weighting appears in **both** `crvec` and `norm` and largely cancels, so
the raw `ref²/ptcl²` amplitude ratio is exposed directly in `v`. Therefore the
matcher is only stable when reference reprojection amplitude ≈ particle
amplitude.

Diagnostic evidence (one crash run, stage 3 iter 41; reference = the stage-2
base-PCG volume, ML off):

| quantity | value |
|---|---|
| `max|pfts_refs|` (reference reprojection) | 11.9 (other runs: up to 33.8) |
| `max|ptcl pft|` (particle) | 0.008 |
| ratio `ref/ptcl` | ~1500× (other runs ~4000×) |
| `v` max | 177 (≫ 23 threshold) |
| `sigma2` min/max | 3.5e-6 / 0.87 |

So the reference volume reprojects ~1500–4000× larger in amplitude than the
noise-normalized (`norm_noise`, unit-variance real space) particles. That
amplitude mismatch is the proximate cause of the overflow.

---

## 4. Decisive experiment: PCG from stage 1

`PCG_REC_START_STAGE` was changed `3 → 1` (see §8) so PCG runs from stage 1.
Result: still crashes at stage 3, and the iter-41 reference (a **stage-2 base
PCG** volume, ML **off**) already overflows euclid. This proves:

- The scale error is **NOT** ML-specific (ML doesn't turn on until stage 3;
  the offending reference is a base-only stage-2 map).
- The crash tracks the **euclid onset** (stage 3), not the PCG onset. Stages 1–2
  use `cc` (scale-invariant), so a mis-scaled/corrupt map is invisible there and
  detonates the instant euclid (scale-sensitive) activates.

---

## 5. Scale instrumentation results (the key data)

A temporary `PCG SCALE DIAG` print was added to `reduce_solve_state_half`
(distributed base+ML solve), printing the raw solved-map amplitude before any
ML rescale. Reading the stage-1 vs stage-2 base solves:

| | BOX_CROP | SMPD_CROP | MAX\|X\| | L2 | DATA_SCALE | LAMBDA |
|---|---|---|---|---|---|---|
| **Stage 1 (healthy)** | 88 | 1.5593 | ~7 (stable) | ~350 (stable) | ~1.2–1.7e8 | 1.000e-3 |
| **Stage 2 (ruined)** | 88 | 1.5593 | 6–20 (unstable, growing) | 350–1900 (unstable) | ~5e5 (discrete ~290× drop) | 1.000e-3 |

Critical observations:

1. **`BOX_CROP`/`SMPD_CROP` are constant (88 / 1.5593) across all stages.**
   The earlier "box_crop changes between stages" hypothesis is **refuted** —
   geometry is not the variable.
2. `DATA_SCALE` drops **discretely** ~290× at the stage-1→2 boundary (a switched
   policy, not drift), and the map amplitude simultaneously destabilizes and
   grows.
3. `LAMBDA` is pinned at the absolute default `1e-3` throughout.
4. The **stage-2 map is genuinely corrupt** (user inspected the maps), not merely
   rescaled — `MAX|X|` and `L2` swing wildly (7→20, 350→1900) iteration to
   iteration, which is a bad/unstable solve, not a constant-factor scale error.

### `data_scale` is inert as a cause (important)

`data_scale` exists **only** to support relative-lambda:
`update_lambda_from_density` does `if(l_lambda_relative) lambda = lambda_rel*data_scale`.
Relative lambda is **not enabled** (LAMBDA stays absolute `1e-3`), so `data_scale`
feeds nothing in the solve. Its 290× collapse is a **thermometer, not the fire**.
`data_scale` (and the `l_lambda_relative`/`lambda_rel`/`set_lambda_relative`/
`pcg_lambda_rel` machinery) is a candidate for removal as dead code — but that is
cleanup, orthogonal to the bug.

However: `data_scale` is computed from the raw density `D`, and `Khat`, the RHS
`b`, and the preconditioner `rho` are **all built from that same `D`**. So the
collapse it reports reflects a real change in the accumulated operator that DOES
feed the (fixed-lambda, 2-iteration) solve. The right quantities to instrument
next are `b`, `Khat`, and `rho` conditioning directly — not `data_scale`.

---

## 6. Current best hypothesis

On this dataset, at stage 2 the accumulated data operator `D` (hence `Khat`/`b`/
preconditioner `rho`) changes discretely, and the euclid per-shell `sigma2`
enters with an enormous dynamic range (stage-3 diag showed `sigma2`
3.5e-6 … 0.87, ~250,000× spread; the 3.5e-6 low-shell value spikes `1/sigma2`).
That wrecks the conditioning of `rho`/the preconditioner. The PCG solve is a
fixed **2 CG iterations from x=0** — provably enough on well-conditioned data,
but insufficient for the ill-conditioned stage-2 operator on this dataset — so
it emits a corrupt, over-scaled map. The map then reprojects ~1500–4000× too
large and overflows euclid at stage 3.

This is consistent with every constraint:
- geometry constant (box_crop 88 throughout),
- not ML (base stage-2 map already overflows),
- "more than just sigma values" (it's sigma's dynamic range hitting the
  preconditioner + fixed 2-iter budget → a corrupt map, not a rescaled one),
- 2 iterations works on well-conditioned datasets.

**NOT yet proven.** Needs the isolation in §7.

---

## 7. Recommended next step (decisive isolation)

Reconstruct the **identical stage-2 particles/orientations** with the gridding
backend (`rec_backend=gridding`) and compare the map to the PCG map from the
same inputs:

- Gridding clean **and** PCG blown up from the same `b`/`D` ⇒ inputs are fine;
  fault is the **PCG solve conditioning** (preconditioner/lambda/iteration
  interaction) → fix lives in `simple_reconstructor_pcg` (preconditioner and/or
  make lambda track the operator scale so the solve stays well-conditioned).
- Gridding **also** off ⇒ the inputs (`D`/`b`, i.e. `sigma2`/CTF/orientation
  coverage) are corrupt → fix upstream.

Complementary instrumentation (cheap): extend the DIAG to print RHS `b` L2 and
`rho` min/max (the real solve inputs) instead of / in addition to `data_scale`.

---

## 8. Uncommitted code changes made during this investigation

All changes validated with editor/language-server diagnostics only (no compile).

### Temporary — MUST be removed before any commit

1. **`src/main/pftc/simple_polarft_corr.f90`** — one-shot `EUCLID DIAG` block
   (guarded by `logical :: l_euclid_diag_done`) in
   `gen_prob_likelihood_euclid_val`, printing `norm`, `wsqsum`, `crvec`, `v`,
   `sigma2`, `max|pfts_refs|`, `max|ptcl pft|`, `max|ctfmat|`.
2. **`src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`** —
   `PCG SCALE DIAG` print in `reduce_solve_state_half` (after
   `validate_solved_map`, before `rescale_ml_map_to_base`), printing
   `BOX_CROP`, `SMPD_CROP`, `MAX|X|`, `L2`, `DATA_SCALE`, `LAMBDA`.

### Experimental — revert or reconsider once root cause is fixed

3. **`src/main/simple_abinitio_controller.f90`** — `PCG_REC_START_STAGE`
   changed `3 → 1` (to run PCG from stage 1 for the isolation experiment).
   Revert to `3` unless intentionally kept. `ML_REG_START_STAGE` unchanged (3).
4. **`simple_rec3D_pcg_strategy.f90`** — `rescale_ml_map_to_base` helper (clamps
   ML map to base map whole-volume L2) + ML solve zero-start. This guard **did
   not fix** the crash (wrong metric: whole-volume L2 is solvent-dominated and
   does not reflect low-shell reprojection concentration). Reconsider/remove
   when the real fix lands.

### Likely-keep (separate task, done earlier this session)

5. `inpl_cont` CLI exposure for abinitio2D/3D, abinitio3D_cavgs,
   refine3D_auto/multi (`simple_ui_abinitio3D.f90`, `simple_ui_refine3D.f90`).
6. Logging reductions: one-time rotational-storage notice
   (`simple_syslib.f90` `l_io_rot_notice_emitted`); removed CLUSTER2D/REFINE3D
   CONTINUOUS aggregate logging and gated cache messages to `part==1`
   (`simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90`).

---

## 9. Ruled out

- Empty volume / invalid reprojection bin / KB weight normalization / padsc²
  mismatch (earlier session).
- ML regularization as the cause (base stage-2 map overflows on its own).
- box_crop / smpd_crop change between stages (constant 88 / 1.5593 throughout).
- PCG 2-iteration budget being universally too small (works on other datasets).
- `data_scale` as a direct cause (feeds only the unused relative-lambda path).
- `rescale_ml_map_to_base` whole-volume-L2 guard (did not fire usefully / did
  not fix).

---

## 10. Key code locations

- `src/main/simple_eul_prob_tab.f90:611-612` — the THROW ("all reference
  distances invalid").
- `src/main/pftc/simple_polarft_corr.f90` — euclid kernel
  (`gen_prob_likelihood_euclid_val`, `euclid_dist_from_crvec`) + temp diag.
- `src/main/pftc/simple_polarft_memo.f90` — `memoize_refs`, `ft_ref`/`ft_ref2`
  (no normalization of reference reprojections).
- `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`
  - `solve_state_half` (~224, shared base solve)
  - `regularize_state_half` (~352, shared ML solve)
  - `reduce_solve_state_half` (~1416, distributed base+ML solve — the path the
    failing run uses; workers print `PCG RAW WORKER`, master `PCG DISTRIBUTED`)
  - `rescale_ml_map_to_base` (~1702)
  - worker `accumulate_raw` (~820) — builds raw accumulator; `sig2` from
    `resample_sigma2(..., R, 1.0, sig2)` (scale arg hardcoded 1.0)
  - `prepare_fused_planes` in the PCG op — RHS `b = conjg(T)·y·sw` and
    `absT2 = cval²/sig2`; both carry one power of `1/sig2` (consistent, verified).
- `src/main/volume/simple_reconstructor_pcg.f90`
  - `update_lambda_from_density` (~2720) — computes `data_scale`, sets relative
    lambda (unused path)
  - `calibrate_kernel` (~3149) — `Khat *= padsc²` (analytic, padsc=padf³=8)
  - `fold_and_ifft` / `pad_vol` / `crop_vol` — FFT round-trip (norm-preserving
    embeds; padded-box FFT)
  - `apply_normal_kernel` / `apply_normal_matrixfree` / `apply_precond` — the
    operator and preconditioner
  - `get_data_scale` (~3786), `get_effective_lambda` (~3791)
- `src/main/volume/simple_reconstructor.f90:583` — gridding `sampl_dens_correct`
  (divide by rho) — the data-scale normalization PCG lacks an equivalent of.
- `src/utils/math/simple_math_ft.f90:277` — `resample_sigma2` (linear-interp
  resample; `input_shell_step` scale arg).

---

## 11. Open questions for the next system

1. Gridding vs PCG on identical stage-2 inputs: is the map corrupt from bad
   inputs or a bad solve? (the decisive test, §7)
2. What exactly switches at the stage-1→2 boundary to drop the operator `D`
   ~290× discretely? (objfun/sigma weighting turning on? sampling/update_frac?)
   The change is discrete → a stage-gated policy.
3. Does the euclid `sigma2` dynamic range (3.5e-6 … 0.87) blow up the
   preconditioner `rho`? Instrument `rho` min/max.
4. Should the real fix be: (a) make PCG lambda track the operator scale so the
   2-iter solve stays conditioned; (b) floor/clamp low-`sigma2` spikes feeding
   `D`; (c) a data-scale output normalization matching gridding's
   `sampl_dens_correct`; or a combination?
5. Fold the `data_scale`/relative-lambda dead-code removal into the real fix.

---

## 12. Forensic tooling on the run box

`/tmp/mrcstat.py` (pure-Python `struct`/`array`, numpy-free) prints
dim/min/max/mean/sd/L2 for `recvol_state01*` maps in the run dir. Useful for
comparing merged/ML vs base half-maps and (next) PCG vs gridding outputs.

---

## 13. Follow-up analysis (2026-08-20): the beyond-band spectral defect

Source: full re-read of `ABINITIO3D_OUTPUT_RESTART1` (a restart run with
`PCG_REC_START_STAGE=1`, euclid diag + scale diag active from iteration 1;
stages here are 20 iterations each: stage 1 = iters 1–20, stage 2 = 21–40,
stage 3 = iter 41 → crash).

### 13.1 The measured event timeline (all numbers from the run log)

| iter | event | max\|pfts_refs\| | sigma2 min/max | DATA_SCALE | RESID (2 its) | map L2 |
|---|---|---|---|---|---|---|
| 1–20 (stage 1, kto=14) | healthy | 0.0050–0.0058 | 4.2e-6 / 9.9e-6 | 1.25–1.70e8 | 0.22–0.25 | 309–357 |
| 21 (stage 2 entry, kto→15) | **refs jump ×345**; map & sigma2 unchanged | **1.76–1.91** | 3.4e-6 / 9.5e-6 | 1.68e8 | 0.22 | 353 |
| 22 | sigma2 absorbs shell 15; weights collapse | 1.91 | 3.6e-6 / **4.4e-2** | **5.9e5** | 0.47/0.57 | 760 |
| 23–40 | runaway feedback | 5.4–34.1 | — | ~5e5 | **0.86–1.00** | 350–1946 |
| 41 (stage 3, kto→16) | euclid all-invalid THROW | 11.9 | 3.5e-6 / 0.87 | — | — | — |

Two additional log facts: the euclid `v` max at iter 21 is 448–547 (not the
~120,000 an all-shell ×345 rescale would give), and the FSC=0.143 resolution
jumps from ~8.6 Å (stage 1) to **3.1 Å = Nyquist** by iter 30 with lp=9.1
alignment — the half-maps share correlated garbage across the whole band (no
gold standard in abinitio3D; both halves align to common refs).

### 13.2 What this proves (corrections to §§4–6)

1. **Nothing "switches" in reference prep or the solve at the stage-1→2
   boundary.** The switch is simply `kto: 14 → 15` (lp 9.8 → 9.1). At iter 21
   the volume is bit-identical in scale to iter 20 (SCALE DIAG), sigma2 files
   are unchanged, yet `max|pfts_refs|` jumps ×345. The polar reference model
   simply starts sampling shell 15 — and shell 15 of every stage-1 PCG map is
   ~345x hotter than its in-band (k≤14) content. Confirmed twice: iter-21 refs
   (from the stage-1-final map) and iter-22 refs (from the iter-21 map, solved
   pre-pollution) both read 1.9.
2. **In-band content stays healthy through iter 21.** The freshly estimated
   sigma2 at iter 22 keeps its minimum (3.6e-6, shells ≤14 residuals unchanged)
   while the maximum explodes ×4700 — the residual blow-up is confined to the
   newly matched shell. An all-shell scale error is excluded.
3. **The DATA_SCALE ~290x collapse (§5) is a downstream symptom, one iteration
   later.** calc_group_sigmas at iter 22 estimates sigma2(15) from iter-21
   residuals against the hot shell → 4.4e-2. `resample_sigma2` flat-extends
   sigma2(kto) to Nyquist, so the reconstruction weights of shells 15–44 (~96%
   of Fourier voxels) drop ×4600, and D's aggregate scale drops ×290.
4. **The solve then genuinely stops converging.** Relative residual after the
   fixed 2 iterations: 0.23 (iters ≤21) → 0.9–1.0 (iters 23+). The wild
   MAX|X|/L2 swings of §5 are this non-convergence, not the root cause. The
   sigma2 dynamic-range/preconditioner hypothesis of §6 is right, but about the
   *propagation*, not the origin.
5. **The euclid crash needs two ingredients that only meet at stage 3**: refs
   hot at ≥2 matched shells (15,16) plus enough particles where *every*
   reference is invalid. Stage 2 survives on warnings because sigma2(15)
   rebalances `v` (max 3.8 at iter 22) once estimated.
6. In the **original run** (PCG only from stage 3), stages 1–2 are gridding
   (clean beyond-band, quotient-normalized), so the first exposure happens at
   stage 3 where PCG + ML + euclid all start together — same defect class,
   different first-exposure site (crash at iter 41/42).

### 13.3 The open core question, sharply posed

Why do PCG base solves leave shell-(kto+1…) amplitudes at ~whitened-data scale
(observed ×345 over in-band; ×~1000 over the expected noise-suppressed
quotient)? Constraints established by code reading:

- The accumulation is full-band (`sqlp = Rnat²`; planes to native Nyquist) at
  every stage; there is no k-space data cliff at the alignment band.
- `b` and `D` carry consistent `1/sigma2` powers (`prepare_fused_planes`
  re-verified), so the pointwise quotient `b/D ≈ y/ctf` is sigma2-independent —
  a *diagonal* solve cannot produce the hot shell; neither can from-zero CG
  *undershoot* explain it (weak modes stay near x=0, they cannot overshoot).
- Therefore either (a) the quotient itself is hot at these shells in the
  Cartesian plane convention (per-component noise at shell 15 >> per-component
  signal at shells ≤14, insufficiently suppressed by sqrt(N_eff) on the padded
  lattice), and gridding differs only via its pointwise normalization +
  compressed-lattice averaging; or (b) a non-diagonal effect (preconditioner
  shell-floor at sampling-gap voxels, Khat/mask shell coupling, α global
  scaling) parks energy there. The 1.9 ≈ O(1) whitened scale and
  345 ≈ 1/sqrt(sigma2(14)) coincidence favor a whitened-b-scale mechanism.

### 13.4 Decisive instrumentation (cheap, one run)

Extend the existing `PCG SCALE DIAG` (reduce_solve_state_half) with a per-shell
amplitude profile of the solved map around the band edge, e.g. after
`validate_solved_map`:

    ! shell-resolved |X| profile, native shells 10..20
    call volume%new(...); call volume%set_rmat(x,.false.); call volume%fft()
    ! for sh in 10..20: print sqrt(mean |X(sh)|^2)

plus the same profile for `b` (get_rhs) and `rho` min/mean/max per shell.
Prediction: at stage 1 the |X| profile shows a step of ~×300 between shells 14
and 15. Where the step tracks `b`'s profile → quotient/convention issue (a);
where it appears only in `X` → solver effect (b). The §7 gridding-vs-PCG
same-inputs comparison remains the gold-standard cross-check — compare *shell
profiles*, not whole-map stats (whole-volume L2/max hide the defect; that is
why the §8.4 whole-volume-L2 ML guard could never fire — and in this run it
never did: all failing solves are `KIND=base`, which the guard does not wrap).

### 13.5 Fix candidates (solver level)

The map keeps its full band throughout: resolution extending beyond the
matching limit is the validation of the map, so the fix belongs in the solve.

1. **Warm-start the base solve from the preconditioned quotient** `x0 = P·b`
   (≈ the gridding solution) instead of x0=0, so the 2-iteration budget refines
   a data-consistent spectrum rather than leaving transients; costs one extra
   operator application. (Distinct from the §8.4 ML warm-start, which starts
   from the *base map* for a *different* operator.)
2. **Put the shell-floor into the operator, not only the preconditioner** — a
   shell-relative Tikhonov Λ(sh) = frac·mean_rho(sh) (the absolute λ=1e-3 is
   ~1e-11 of D and regularizes nothing). This damps sampling-gap voxels
   consistently in H and P.
3. **Sigma2 hygiene at stage transitions**: when kto grows, the first estimate
   of the new shell is made against references sampled from a map that was
   never constrained by matching at that shell; once the solver's beyond-band
   content is data-consistent this is the normal euclid band extension.

### 13.6 Answers to §11's open questions

1. (gridding vs PCG on identical inputs) Still worth running, but compare
   shell profiles; the in-band inputs are demonstrably fine.
2. (what switches at stage 1→2) Nothing in policy: `kto` grows 14→15 and
   exposes pre-existing beyond-band junk. The ~290x DATA_SCALE drop happens at
   the *second* stage-2 iteration, via the first sigma2 estimate that includes
   shell 15.
3. (does sigma2 range blow up rho) Yes, but as propagation: post-pollution
   D spans ~4 extra orders across the band edge and the fixed 2-iteration solve
   returns residual ~1.0.
4. (which fix) See §13.5 — (1)/(2) are the structural candidates;
   data-scale-tracking lambda alone would not remove the beyond-band defect at
   stage 1 conditions (sigma2 was flat then; the defect predates any sigma2
   dynamic-range issue).
5. (dead code removal) unchanged, orthogonal.

Restart-specific caveat: this restart's aggressive LP schedule (stage 1 already
at 9.8 Å) was derived from the original run's FSC files, which were themselves
contaminated to near-Nyquist by the same defect — FSC from PCG maps without
gold-standard splitting inherits the correlated beyond-band junk (8.07–8.58 Å
"resolution" during stage 1 with lp 9.8 alignment; 3.1 Å by iter 30).

---

## 14. Defect A diagnostic and cleanup (2026-08-20/21)

`report_beyond_band_excess` (module-level, `simple_rec3D_pcg_strategy.f90`,
called after every solve in both execution paths) compares the RMS of the
shells beyond the matching band with the band-edge shell and logs
`>>> PCG BEYOND-BAND EXCESS ... BEYOND/EDGE RMS RATIO=` when the ratio
reaches 10. The map is not modified. This is the regression signal for the
§13 defect: it would have printed ~3e2 every iteration of the restart run, and
it is silent when the solver behaves. The matching band is
`params%kfromto(2)` (`matched_band_kstop`, 0 = unknown → silent), read the
same way in every execution path.

The 2026-08-21 rerun (default `PCG_REC_START_STAGE=3`, proper sigma2) showed
the beyond-band content of the stage-3 PCG solves to be mild (ratio ≤ 10), so
in that configuration defect A is not active; it was severe only in the
from-stage-1 configuration with bootstrap-scale, near-flat sigma2.

**Cleanup of the original investigation's changes:**

- Removed the temporary `EUCLID DIAG` (simple_polarft_corr.f90) and
  `PCG SCALE DIAG` prints (§8.1–8.2).
- Removed the `rescale_ml_map_to_base` whole-volume-L2 guard and the ML
  zero-start experiment (§8.4); the guard compared the wrong pair of maps
  (ML vs base within one iteration — see §15 for the pair that matters) and
  never fired on the failing base solves; with the right anchor in place the
  original ML warm start is the better-converging choice.
- `PCG_REC_START_STAGE` reverted to 3 (§8.3).

---

## 15. Second defect (2026-08-21 rerun): absolute-scale discontinuity at the gridding→PCG handoff

The workflow was rerun with the default `PCG_REC_START_STAGE=3`, so stages
1–2 used gridding and PCG (base+ML, euclid-weighted) first activated at stage
3. The crash reproduced in the ORIGINAL mode: iteration 41 aligned fine on the
stage-2 gridding references, both PCG solves converged (RESID 0.11 base /
0.07 ml), and iteration 42's prob table threw all-invalid on the first PCG
references. Shell spectra of the written iteration-41 maps are smooth and
falling across the matched band with no excess beyond it (ratio ≤ 10).
Defect A was not the killer in this configuration.

### 15.1 The measurement (run dir 6_abinitio3D, shellspec on the actual maps)

| map | in-band character | L2 |
|---|---|---|
| `recvol_state01_stage01_lp.mrc` (stage-1 final, **gridding**, lp snapshot) | smooth, falling | **0.99** |
| `recvol_state01_iter041_lp.mrc` (iter-41, **PCG**, lp snapshot) | smooth, falling, same shape | **194** |
| `recvol_state01_even_unfil.mrc` (iter-41 PCG base half) | smooth | 206 |
| `recvol_state01_even.mrc` (iter-41 PCG ML half) | smooth | 196 |

Per-shell RMS ratios PCG/gridding over the deeply shared band (shells 0–10):
~150–260, i.e. a **uniform ~×200 amplitude convention difference** between the
PCG solution and gridding's `sampl_dens_correct` data-quotient convention
(base and ML equally — the ML solve is not the cause).

### 15.2 Why ×200 is lethal at the handoff but not in steady state

The euclid objective is scale-sensitive, but the group sigma2 estimates
equilibrate to whatever stable reference amplitude they see (the §13 restart
run aligned happily for 20+ iterations on PCG-scale maps, refs 0.0055 with
sigma2 ~5e-6; the gridding stages align happily at gridding scale). What
euclid cannot survive is an abrupt scale change between consecutive
iterations: references ×200 → `crvec` ×4e4 → `v = 1 + crvec/norm` far above
the 23.03 invalidation threshold for EVERY rotation and reference of every
particle, one iteration before calc_group_sigmas could re-equilibrate. The
gridding→PCG backend switch inside stage 3 is exactly such a jump; the
original run's iter-41/42 crash is this defect, not defect A. (Note the
asymmetry: refs *below* particle scale drive v→1 — flat but usable — which is
why the gridding stages ran fine even though their reference amplitudes are
tiny in the polar convention.)

### 15.3 Fix implemented: PCG adopts the gridding amplitude convention

`pcg_mag_correction`/`apply_output_convention` (module-level,
`simple_rec3D_pcg_strategy.f90`): every solved map, in both execution paths
and for both solve kinds, is divided by `params%box` right after the solve —
exactly what `reconstructor_eo` does to its maps (`mag_correction`). The
solves themselves stay in the solver's native convention: a stored base map
that seeds an ML solve is multiplied back by the same factor before it is
used as the warm start. With this the first PCG map after the gridding stages
lands on the same scale as the map it replaces, the euclid/sigma2
equilibrium is undisturbed, and no per-iteration measurement or rescaling is
needed (the interim `apply_scale_continuity` anchor was removed).

The measured remainder (~1.5 over 128) is not a convention: it reflects the
different alignment states of the two maps compared in §15.1 (the ratio rises
with resolution) plus second-order KB-deconvolution differences, and is
within what sigma2 absorbs between consecutive iterations. The same-inputs
dual-backend test in `drop_legacy_box_division.md` §5.2 pins it down.

### 15.4 Root cause: `mag_correction = box`, and where the deapodization sits

Tracing both chains with SIMPLE's transform convention (`fft` = (1/N)Σ,
`ifft` = Σ, a true pair, so Fourier cropping and zero-padding are
amplitude-neutral): gridding rescales the padded 2D plane by `pf²` back to
the native coefficient convention, inserts with per-axis-normalized KB
weights into numerator and `rho` alike, divides pointwise, inverse-transforms
on the native lattice — and then divides by the ORIGINAL box size
(`simple_reconstructor_eo.f90:119`, `mag_correction = real(params%box)`,
"consistent with the current scheme"; applied at lines 522/532/571/582/601/602
and 770). PCG's factors (`padsc` in `b`, `padsc²` in `Khat`, the 176-lattice
interpolation, the same normalized weights, the deapodization bracket on `H`
and `b`) all cancel, so its solution is the plain data quotient. PCG/gridding
= box = 128, against 160–260 measured on non-identical maps.

The gridding deapodization is not in the reconstructor objects: `volassemble`
(`simple_commanders_rec_distr.f90`, `restore_state_from_parts`) multiplies the
merged volume and the nonuniform-filter source halves by `gridcorr_img =
prep3D_inv_instrfun4mul(ldim, ldim_pd, smpd_crop)`; the half-map files are
not corrected; `filter_pcg_nonuniform_maps` applies no gridcorr to PCG maps
(deapodized inside the solver), so there is no double correction. The
correction itself is wrong for the current reconstructor: it evaluates the
continuous KB instrument function at `1/ldim_croppd` — the envelope of a
window on the 2x padded lattice — while the reconstructor inserts on the
native lattice with a 1.5-native-voxel window, whose envelope (period `box`)
is twice as steep. Gridding maps are therefore under-deapodized by
R_grid = 0.94/0.79/0.62/0.53 at r = 10/20/30/40 (on axis, box 88).

Same-alignment measurement (stage-2 gridding snapshot = iteration 40, vs the
PCG iteration-41 map; both lp snapshots of the merged volumes): L2 ratio
inside r < 24 is 175–177. The box convention alone gives 128; box × gridding's
envelope deficit gives 140–156; the remaining ×1.2–1.26 rises with spatial
frequency and toward the core, consistent with the shift search switching on
at stage 3 (`trs = 0` in stages 1–2). Decomposition of the ~200 first measured
on stage-1 vs stage-3 maps: ×128 convention, ×~1.15 gridding envelope
deficit (streptavidin is small; up to ~1.8 at the periphery of a box-filling
particle), ×~1.2–1.4 map/alignment differences. Details and the fix for the
envelope in `drop_legacy_box_division.md` §2.1 and §5.3.

Why the convention should go rather than be propagated, and the plan for
doing so, are in `drop_legacy_box_division.md`.

---

## 16. Design principle applied on review (2026-08-21)

**Shared-memory and distributed execution must be methodologically
identical.** They are two parallelizations of one algorithm. The output
convention and the diagnostics are implemented once at module level and
invoked the same way from `execute_rec3D_pcg_shared` and
`execute_rec3D_pcg_distributed_master`; no entry point receives special
arguments, and nothing depends on which path produced the maps.
