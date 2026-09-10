# The PCG reconstruction backend: what it is today

Status: current as of 2026-09-10. Present tense only. This page describes
the backend as it runs on master; it does not argue for it. The binding
contracts are the policy documents (`doc/policies/3D/reconstruct3D_pcg_policy.md`,
the PCG sections of `refine3D_policy.md` and `refine3D_auto_policy.md`,
`automasking_policy.md`); when they and this page disagree, the policy is
right and this page needs a fix. How the backend got here, decision by
decision, is `pcg_decision_log.md`; the full experiment record is
`pcg_priors_history.md`. New decisions are appended to the log, and the
affected paragraph here is rewritten, never annotated.

## 1. What it is

`rec_backend=pcg` replaces the gridding density quotient with the weighted
least-squares reconstruction: it solves the normal equations of the
CTF/sigma-weighted Fourier-slice projection operator by preconditioned
conjugate gradients, on a real-space support, with the optional ML prior
inside the operator. The particle pass is done once per iteration into
raw accumulators; the operator is then applied in kernelized (Toeplitz)
form, so the cost of a CG iteration is a handful of FFTs and independent
of the particle count. `pcgop=kernel` is required in production; the
matrix-free operator is the exact reference and exists for tests.

Both backends produce the same kinds of products (even/odd/merged state
volumes, an `_unfil` base pair under `ml_reg=yes`, FSC/cFAR and the
resolution text), at the same amplitude convention (the data-quotient
convention; deapodization is inside the solver), and feed the same
assembly-owned nonuniform (NU) competition afterwards. The PCG backend
carries no prior of its own beyond the FSC/SSNR precision `P_tau`; every
NU-derived prior that was tried inside the solve has been removed
(decision log, 2026-08-27, 08-29, 09-06).

## 2. The two solves

Per state and per half, from one particle accumulation:

- **Base solve** `(H + lambda_0 I) x = b`, `lambda_0 = PCG_LAMBDA = 1e-3`
  (absolute Tikhonov ridge). Its half pair is the `_unfil` pair: the FSC,
  cFAR and resolution authority, and the seed of the NU candidate bank.
- **ML replay** `(H + P_tau + lambda_0 I) x = b`, replayed from the same
  raw accumulators with the FSC/SSNR shell-diagonal precision `P_tau`
  (`tau`, `hp`) in both the operator and the preconditioner. Its pair is
  the shipped map and the NU auxiliary member.

**Starts.** The base solve starts from zero. The replay starts from the
current same-half base solution passed through the closed-form shrinkage
filter (each shell scaled by its FSC-implied Wiener factor, the `P_tau`
optimum in closed form; logged as `PCG ML REGULARIZED INIT`), unless that
start is worse than zero (initial relative residual above 1), in which
case the solver discards it before the first iteration -- at no cost, a
zero start's residual is `b` itself -- and the log says so (`start worse
than zero ... discarded, solved from zero`); this is the normal case in
the low-resolution stages, where the shrinkage start has `INIT` 2-4. There
are no cross-iteration warm starts: nothing from a previous iteration's
half maps enters a solve (2026-09-10).

**Budget.** `maxits_pcg=2`, `rtol=0` (exactly two iterations) in
refinement; the original-sampling final reconstruction uses at least
`FINAL_PCG_MAXITS_FLOOR = 5`. Two iterations from a state-free start is
the calibrated regime (simulated data: two iterations beat gridding,
beyond five the residual moves but nothing interpretable in the map does).
A solve from a nonzero start that loses positive-definiteness is retried
once from zero; a solve from zero that loses it is fatal.

**Preconditioner.** The sampling-density diagonal `1/(rho + floor + P_tau)`
with a shell-relative floor (`RHO_FLOOR_FRAC`).

## 3. The support

The solve is constrained by a real-space support `P`: the system solved
is the projected `P H P u = P b` on the hard domain `window > 0`, and the
shipped map is `window * u` (`PCG_HARD_SOLVE_SUPPORT`). Which window:

- **`automsk=no`**: the soft spherical support at `msk_crop`, the same
  `mask3D_soft` the gridding restoration applies after deapodization.
  Base and replay both.
- **`automsk=yes`**: the conservative density envelope -- `automask3D` of
  the lag-one reference volume at `envmsklp` (20 A), Otsu core, largest
  component, spherical dilation by `binwidth` layers, outward cosine skirt
  of `edge` voxels -- constrains BOTH the base and the replay once a
  prior reconstruction exists. The first reconstruction bootstraps the
  base on the sphere and builds the replay support from its own pair.
- **`pcg_mskfile=<vol>`**: an explicit [0,1] volume, development mode,
  constrains every solve regardless of `automsk` and is reported as the
  state support.

The dilation has a shared physical minimum, `ENVMSKWIDTH_A_MIN = 7.5 A`
(`binwidth = max(binwidth, ceiling(7.5/smpd_crop))` whenever the envelope
is in use; an explicit `binwidth` wins in either direction). The support
is the one place a mask enters the estimate; no PCG product is multiplied
by any mask after the solve (postprocess, `_pproc`, `_mirr` included),
and the support-provenance sidecar `<vol>_pcg_support.txt`
(`solve_support=density|sphere`, `solve_kind=base|regularized|mixed|gridding`)
records what the shipped pair carries.

## 4. The FSC and the envelope

`automsk=yes` implies `envfsc=yes` on both backends (derived in
parameter validation, logged when it overrides a `no`). The FSC is then
obtained in one of three ways, and the way used is REPORTED on every
evaluation as `>>> FSC MODE` in the log and in the resolution text:

- spherical support from the reconstruction, no envelope, no correction
  (`envfsc=no`);
- density envelope applied post hoc to the halves with the
  phase-randomized solvent correction (gridding, `envfsc=yes`);
- halves estimated on a support envelope, no post-hoc mask and no
  correction (PCG under `automsk=yes`, or `pcg_mskfile`). The envelope
  contribution to the FSC is deliberately NOT removed; a common window on
  both halves can contribute correlated power, and the policy is to say
  so rather than hide it.

The `automask3D_stateNN.mrc` artifact is written on the envfsc path for its
other consumers (postprocess of non-PCG products, the abinitio final rec).

## 5. NU filtering and the evidence envelope

The NU competition is assembly-owned and identical on both backends
(`simple_nu_state_filter`): the discrete candidate bank is built from the
base (`_unfil`) pair, the replay pair joins as the finest auxiliary member
(`ml_reg=yes, nu_refine=no`), the shell walk extends the bank
(`nu_refine=yes`, refine3D_auto), and the finest selected label with at
least 5% assigned support is the matching low-pass handoff. The NU
objective always runs on the spherical `mskdiam` support.

Under `automsk=yes` the NU evidence envelope (`nu_envmask3D_stateNN.mrc`,
regenerated every competition from the live evidence) is the mask that
controls the filtering: the filter-field background is its complement, and
background voxels take the coarsest bank candidate, so excluded density
(detergent, disordered belt) reaches the matching references heavily
low-pass filtered rather than removed. References are never multiplied by
an envelope. The evidence envelope's null model has two regimes, keyed on
how the base pair was solved:

- spherical base pair (gridding; PCG bootstrap): robust median +
  `nu_msk_sig` MAD of the margin over the observed support; valid while
  solvent holds the majority;
- envelope-constrained base pair (PCG `automsk=yes`, `pcg_mskfile`): the
  estimator has removed the far solvent, so the null is designated by
  Euclidean geometry -- median/MAD on the density envelope's dilation
  ring at full weight of the base support (`set_nu_evidence_null_shell`);
  labels are free on the observed density envelope and fixed solvent
  outside it, so the evidence envelope is nested inside the density
  envelope; valid while the shell is populated.

If the null is invalid or the envelope empty, the density envelope itself
is armed as the background (`EVIDENCE FALLBACK`), and the provenance
string records which envelope ran.

## 6. Where it runs

- **refine3D / refine3D_auto / refine3D_states**: `rec_backend=pcg`
  selects the distributed master (`execute_rec3D_pcg_distributed_master`):
  workers accumulate and write raw statistics, the master reduces,
  solves both halves concurrently, runs the FSC, the replay and the NU
  competition, and ships. Trailing reconstruction (`trail_rec`) keeps
  accumulator chains; the bootstrap iteration blends with the previous
  pair at the update weight exactly as the gridding bootstrap.
- **abinitio3D**: the PCG backend engages from `PCG_REC_START_STAGE = 3`;
  stages 1-2 are gridding. Stage policy sets `envfsc`/`automsk` at the
  last stage (`AUTOMSK_STAGE = ENVFSC_STAGE = NSTAGES`). In the NU stages
  the NU handoff sets the matching low-pass, capped at the ladder's hard
  fine bound `LPSTOP_BOUNDS(1)` (4.5 A; a coarser explicit `lpstop` is
  retained), with the stage limit promoted to the previous stage's FSC=0.5
  crossing when that is finer; NU stages run their full iteration budget
  (`minits = maxits`).
- **Final reconstruction** (`bootstrap_rec3D`, both workflows): image-power
  sigma seed, a gridding ML bootstrap map carrying the workflow's
  `filt_mode`/`nu_refine`/`automsk` (the residual sigmas depend on the
  regularization of the reference they are scored against), one residual
  sigma pass, then the shipped PCG map at the native box: cold, at least
  five iterations, `filt_mode=none`, `nu_refine=no`, `automsk` inherited
  (so its support matches the refinement's), `postprocess=yes`.
- **Shared-memory** (`nparts=1`) uses `execute_rec3D_pcg_shared`, the same
  policy in one process; shared-memory trailing is not supported.

## 7. Reading a run

Per half and solve, one summary line:

    >>> PCG DISTRIBUTED | STATE= 1 | HALF=even | KIND=base | N=  8416 | ITS= 2 | INIT= 1.000E+00 | RESID= 6.100E-02 | TIME=   1.9 s | STOP=fixed_iterations

`INIT` is the relative residual of the start (1.0 from zero), `RESID` the
final one. Healthy values on PfCRT at box 140-160: base `1 -> 0.06-0.09`,
replay `~0.2 -> ~0.04` once the map is at 4-5 A; in the low-resolution
stages the shrinkage start reads `INIT` 2-4, is rejected, and the replay
runs from zero (`INIT` 1.0). A replay `RESID` that climbs across
iterations is the signature of a bad start.

Other lines to grep: `PCG SOLVE SUPPORT` (which support, and why),
`FSC MODE`, `PCG ML REGULARIZED INIT`, `NU ENVELOPE OCCUPANCY`,
`NU DILATION RING OCCUPANCY` (how much of the density envelope's dilation
ring the evidence labels signal -- the number to consult before tightening
`binwidth`), `NU NULL SHELL GEOMETRY` (envelope/support Dice, ring
retained at full weight), `NU BACKGROUND` (evidence envelope or fallback),
`NU BANK CAP`, `NU LOW-PASS ASSIGNMENTS`, `NU filter promoted matching
low-pass`, `PCG BEYOND-BAND EXCESS` (post-band RMS >= 10x the band-edge
shell; the regression signal for solver defects), `RECONSTRUCTION MASTER
PHASE` (wall time, thread-seconds, peak RSS). Per-half solve diagnostics
(residual history, data scale, effective lambda, prior statistics) go to
`reconstruct3D_pcg_stateNN_<half>_<kind>*.txt` sidecars.

The stage-6 take-off is the thing to check on an abinitio3D run: at the
first NU iteration the base pair should support a few percent of the mask
at ~8 A, the matching low-pass should be promoted below 6 A, and FSC-0.143
should drop from ~9 to ~4.5 A within the stage. If the matching low-pass
sits at the stage ladder value, the base pair is not carrying fine-shell
evidence.

## 8. Cost

At PfCRT scale (16.8k particles, box 140-160, 10 parts x 8 threads) the
PCG master phase is ~16 s per refinement iteration against gridding's
~4 s, and the refine3D worker step is ~12 s against ~8 s because the PCG
workers write the raw accumulators. On a 133-iteration abinitio3D that is
7460 s against 4820 s. Iterations of CG themselves are cheap (~1 s per
half per iteration under the kernel operator); the cost is the
accumulation and the master's reduction, not the solve.

## 9. Tests

`simple_test_exec test=pcg_recon` (operator identities, support contract,
window band regression), `test=pcg_frac_update` (trailing accumulator
arithmetic), `test=rec3D_backends` (gridding vs PCG on a project with
poses, with optional truth). `doc/implementation_notes/automsk_yes_code_review.md`
lists the nine-case validation matrix for the envelope/NU changes of
2026-09-09.

## 10. What is deliberately not there

No solvent prior, no Wilson prior, no NU-evidence prior in the solve, no
auto-lambda; no cross-iteration warm starts; no post-hoc mask on any PCG
product; no envelope multiplication of matching references; no AMSK_FREQ
regeneration cadence. Each was implemented, measured and removed; the
decision log has the date and the evidence for every one.
