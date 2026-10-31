# Polar shift map: a global shift search on the polar Fourier transform

Planned 2026-10-30. One living note: the contract (the mechanics and what is fixed), the plan (phases, exit
criteria, file table), the decisions still open, and the progress record. The mechanics and its unit test are
written; nothing is wired into a search yet.

## Why

Benchmarks against cryoSPARC on the same particles (Hans, 2026-10-30) found origin shifts after 2D and
solve3D of up to 30 Å relative to the picked coordinates, averaging up to 20 Å on some projects. cryoSPARC
finds them by a discrete search over shifts in 0.5 Å bins. SIMPLE cannot: the in-plane search of a pair is a
5 x 5 coarse grid over the window [-trs, trs]^2 (`coarse_search_opt_angle`, `gen_objfun_vals` at each point)
followed by L-BFGS-B in the window, and the window is the stage's `trslim = min(8, max(2, AHELIX_WIDTH /
smpd_crop))` pixels, about +-10 Å. A particle picked 30 Å off centre is outside the window at every iteration,
the L-BFGS basin is about lp/2 wide, and the shift prior is centred on the start, so the pose never recovers
and the particle is scored wrongly against every reference for the whole run.

The Cartesian remedy is to inverse-FFT the cross spectrum and read all shifts at once (the structure of phase
correlation and of cryoSPARC's search). The polar representation has an exact analogue in which the inverse
FFT over shifts becomes a Fourier-Bessel (Hankel) transform; this note develops it on the calculator's own
conventions and plans its use as the coarse initializer of the windowed L-BFGS search.

## The contract

### The identity

The calculator shifts the reference in its own frame and then rotates it (`gen_corr_cc_for_rot_8_2`,
`gen_corrs`: `shmat * pfts_refs`, then the rotation). On the half ring p = 1..pftsz, theta_p = (p-1) 2 pi /
nrots, the shift matrix of s = (sx, sy) (pixels of the box the calculator's `argtransf` was built on,
`ldim(1)`) is exp(i (ax sx + ay sy)) with ax = k sin(theta_p) 2 pi / box and ay = -k cos(theta_p) 2 pi / box
(`simple_polarft_core.f90`). For a fixed grid rotation irot the shift-dependent part of every objective is the
cross term

    C(s) = sum_k w_k sum_p Re[ Q_k(p) exp(i (ax(p,k) sx + ay(p,k) sy)) ],

with Q_k(p) the pair's cross spectrum in the reference's own indexing: the shift matrix multiplies entry rp of
the reference and the rotation lands that entry on ring position q = rp + irot - 1 (modulo the turn) of the
particle, so Q_k(rp) = REF(rp,k) CTF(q,k) conj(X~(q,k)), X~ the Friedel-extended particle (its second half the
conjugate of the first). For objfun=cc, w_k = k and cc = C / sqrt(sum_k k |CTF REF_k|^2 sum_k k |X_k|^2), the
denominator shift-independent; for objfun=euclid, w_k = k / sigma2_k and the loss is const(irot) - 2 C /
wsqsum. In polar shift coordinates s = rho (cos a, sin a) the argument is ax sx + ay sy = z_k sin(theta_p - a)
with z_k = 2 pi k rho / box, and the Jacobi-Anger expansion exp(i z sin phi) = sum_n J_n(z) exp(i n phi)
separates the sum over the ring:

    C(rho, a) = Re sum_n exp(-i n a) a_n(rho),   a_n(rho) = sum_k w_k J_n(z_k(rho)) H_k[n],
    H_k[n]    = sum_p Q_k(p) exp(i n theta_p)      (the half ring's harmonics; J_{-n} = (-1)^n J_n).

Three consequences:

1. It is exact. It is the same discrete sum the L-BFGS evaluates, with no regridding or interpolation; the only
   truncation is in n. The hand-off to the L-BFGS is therefore from the global maximum of the function it then
   refines.
2. The direction is exact. At fixed rho, C is a trigonometric polynomial in a of degree nharm(rho): the
   inverse DFT over n samples it without error, and any direction can be evaluated. Only rho is gridded.
3. rho = 0 keeps the order n = 0 alone and a_0(0) = sum_k w_k sum_p Re Q_k(p) is the present correlation: the
   existing code is the zeroth term, and the unit test.

### The truncation rule

Past n = z, J_n(z) is positive and decreasing in n (checked numerically for z up to 120). The orders kept at a
magnitude rho are 0 .. nharm(rho), nharm(rho) the last order whose Bessel value at the largest argument
z_max = 2 pi kmax rho / box is at least `SHIFTMAP_JTOL` = 1e-8 (`simple_polar_shiftmap`). In scale-free form
z_max = 2 pi rho / lp with rho and lp in the same units. Orders needed (the usual Wiscombe-type estimate
z + 4 z^(1/3) + 2 for comparison):

| z_max | example | kept at 1e-4 | 1e-6 | 1e-8 | 1e-10 | z + 4 z^(1/3) + 2 |
|---|---|---|---|---|---|---|
| 2.4 | 30 Å at lp 80 Å | 8 | 11 | 13 | 15 | 10 |
| 9.4 | 30 Å at lp 20 Å | 18 | 21 | 24 | 27 | 20 |
| 23.6 | 30 Å at lp 8 Å | 35 | 39 | 43 | 47 | 37 |
| 47.1 | 30 Å at lp 4 Å | 61 | 66 | 71 | 76 | 64 |
| 94.2 | 30 Å at Nyquist, 1 Å/px | 111 | 118 | 124 | 130 | 114 |

### The grid

Magnitudes rho_i (pixels, zero included) and nalpha directions uniform over the turn. The L-BFGS basin is
about lp/2 wide, so a magnitude step of lp/4 (in Å, converted at the search sampling) puts the nearest grid
point inside the basin: nrho = 4 rho_max / lp, i.e. 6 magnitudes at lp 20 Å and 15 at lp 8 Å for a 30 Å
radius. The direction count matters only for the readout (the polynomial is exact); 2 nharm + 1 samples
represent it completely, fewer are fine for a peak search followed by L-BFGS.

### Sampling of the ring

The ring must resolve the phase ramp of the largest shift at the outer shell: nrots >= 4 pi kmax rho / box,
i.e. rho <= pftsz box / (2 pi kmax) pixels. With pftsz = pi msk (`magic_pftsz`) that is rho <= msk lp / (2
smpd), always true for picking errors; it is a condition on the discrete objective itself, not on the
expansion, which is an identity of the discrete sum.

### Cost

Per pair at one rotation: the harmonics (a DFT of the half ring per shell and order, 2 K pftsz (N+1), or one
FFT of the zero-padded half ring per shell, K nrots log nrots), then the Bessel sums nrho (2N+1) K and the
evaluation nrho nalpha (2N+1). Against a brute-force grid of pi nrho^2 points at the same step, each a
K pftsz sum, that is one to two orders of magnitude less.

Over all rotations the production form works from the memoised spectra. On the Friedel-extended full ring the
expansion holds with a_{-n} = conj(a_n) (only n >= 0 is needed, and the full-ring sum is twice the half-ring
objective), the memoised half spectra determine the full ones (`ft_ref`, `ft_ptcl_ctf` keep bins 0..pftsz of
the extended ring's FFT, and the extension gives F^[-m] = (-1)^m conj(F^[m])), and for every order n the
rotational correlation of the reference's spectrum displaced by n bins with the particle's is one inverse
transform of length nrots after the Bessel-weighted shell sum:

    a_n(irot, rho) = (1/nrots) sum_m exp(-2 pi i m (irot-1)/nrots) sum_k w_k J_n(z_k(rho)) R^_k[m-n] conj(X^_k[m]).

That costs (N+1) nrho (K nrots multiply-adds + one FFT) per pair and no forward transforms. The present
coarse grid pays, per grid point, K forward FFTs of the shifted reference plus K nrots products and one
inverse FFT; at equal coverage (nalpha ~ 2N+1 directions) the map is roughly 10-20 times cheaper and exact in
the direction. Both forms are verified numerically (Progress, 2026-10-30).

### What is built

`src/main/pftc/simple_polar_shiftmap.f90`: the class `polar_shiftmap`. `new(pftsz, kfromto, box, rho(:),
nalpha)` builds the tables (Bessel values per order, shell and magnitude; the ring's and the directions'
twiddles; the orders kept); `cross_spectrum(ref, ptcl, ctf, irot, Q)` forms Q of a pair at a grid rotation
in the calculator's convention; `coeffs(Q, w, A)` the coefficients a_n(rho) for any shell weights; `eval(A,
irho, alpha)` the cross term at a magnitude and any direction; `map(A, cmap)` the grid; `peak(cmap, irho,
ialpha, val)` its maximum; `get_shift(irho, ialpha)` the shift of a grid point in the frame the shift matrix
acts in. The harmonics are a direct DFT with a twiddle table (the mechanics, not the kernel); all arithmetic is
double precision; the object is read-only after `new`, so one instance serves every thread.

`src/main/pftc/simple_polar_shiftmap_tester.f90`, sub-suite `polar shift map` of `unit_pftc_align2D3D`
(`simple_test_exec test=unit_pftc_align2D3D suite=polar_shift_map`): the map against the closed-form cosine
sum of a pure phase ramp on and off the grid, the peak of that ramp, the zero-magnitude row and the
truncation rule; the map's normalised cross term against `gen_corr_for_rot_8` (objfun=cc, the reference
shifted and rotated by the calculator) at the four rotations that exercise every branch of `rotate_ref_8`
and at 40 random poses, to 1e-9; and peak recovery of a particle that is the reference shifted and rotated by
a known grid pose, exact and at unit SNR, where the calculator's cc of 1 at the true pose certifies that the
tester's construction follows the calculator's conventions.

## Decisions fixed

- The map scores the same objective the L-BFGS refines: the shell weights are the objective's (k for cc,
  k/sigma2 for euclid), no per-coefficient whitening. Per-coefficient phase normalisation amplifies every shell
  whose per-coefficient SNR is below one, which for a single particle is most of them; the peak width is set
  by the band limit, not by whitening.
- The shift prior is added on the grid (it is quadratic in s and a per-rotation factor, `prior_factors`), so
  the grid maximum is the MAP of the same posterior the L-BFGS polishes; that is also what keeps a disk of
  several hundred hypotheses honest at low SNR.
- Where the map runs (ruling 12). On the refined candidates of a particle, where the L-BFGS runs today (<=
  prob_athres references in 3D, the 2D refined set), in the centring stage as in every other: the global search
  of their shifts before the windowed L-BFGS, never against every class. A picking offset is a property of the
  particle, not of the (reference, rotation) hypothesis: the shift found on the committed candidate is the stored
  shift every entry of the next iteration is scored at, so a mis-centred particle is scored wrongly for one
  iteration and then correctly. Phase 1's timing table decided it: the map at the disk costs about 150 times one
  per-class score of the fill.
- The L-BFGS window after the map is the magnitude step about the grid maximum (the basin), not the search
  radius.
- Same mechanics for 2D and 3D (one calculator).

## Rulings (2026-10-30)

Hans's rulings on the maintainer's decisions, from the design discussion of 2026-10-30, with the schedule
facts read from the code the same day.

1. The radius is the prior's: radius = 3 sigma, one quantity, no separate switch. sigma = 10 Å from the
   benchmark (a Gaussian whose 3 sigma is the 30 Å maximum recovered), so the disk is 30 Å; the grid step
   is lp/4 (5 Å at lp 20 Å, six magnitudes); the L-BFGS limits are +-one magnitude step about the grid
   start, the basin, never the disk.
2. The prior is absolute: a Gaussian centred on the pick with sigma 10 Å, applied on the grid as the
   per-rotation factor so grid and L-BFGS score one posterior. The present random-walk prior on the
   increments (`shift_prior_sigma2`) stays in the stages the map does not run in.
3. One objective for both levels (Decisions fixed). The grid is joint in rotation and shift: the map gives
   (irot, rho, alpha) for every rotation at once, no alternation and no rotation at the wrong shift. The
   direction is a trigonometric polynomial, so the sub-grid start is Newton in alpha and a three-point
   parabola in rho at no objective evaluation. Where the map runs, this replaces `grid_fit` and the
   particle-frame seed bank.
4. Superseded by ruling 12. The first form, the map inside the fill against every class at the cheap stage, was
   a cost claim; Phase 1's timing table (Progress, 2026-10-30) put the map at the disk at about 150 times one
   per-class score of today's fill, above 30 after the obvious optimisations.
5. The in-plane choice is a draw over the three in-plane degrees of freedom (ruled 2026-10-30). The first
   ruling, the greedy best (irot, shift) per class with the draw over classes only, was reopened the same
   day: the tempered in-plane draw has been the mechanism that diversifies the search and converges the hard
   particles, and a greedy shift would make that draw conditional on one shift hypothesis. Today
   `gen_posterior_val` draws the rotation from p_r proportional to exp(-T (L_r - L_min)) over all nrots
   rotations at the stored shift, with T = prob_beta x pi wsqsum / (2 pftsz) (`posterior_temperature`,
   `draw_posterior_row`), the class is drawn from the row posterior, and the shift is local and
   deterministic. The map gives the same normalised loss on the whole (irot, rho, alpha) grid,
   L(irot, s) = const(irot) - 2 C(irot, s) / wsqsum (a shift preserves the reference's power per shell), so
   one `draw_posterior_row` over the flattened grid, prior included, is the joint draw; its cost is the
   map's, the draw itself one cumulative sum. Two conditions: the prior is not tempered, p proportional to
   exp(-beta T (L - L_min)) x prior(s), because tempering it at beta 0.1 widens sigma by 1/sqrt(beta) and the
   disk becomes uniform; and the anneal tail (`l_anneal_tail`, prob_beta_end = 1) ends at the physical
   posterior, where a sharpened class average narrows the shift posterior below the grid step and the draw
   is deterministic in effect. Early draws then spread at the 10 Å scale, inside the band of lp 20-30 Å; this
   is the single-sample form of RELION's marginalisation over psi and the translations under its Gaussian
   translational prior. The greedy choice is the same draw at infinite temperature on the shift degrees and
   stays as the comparison arm of Phase 3 behind the switch (ruling 11). To measure: the entropy of a grid of
   some 45,000 points against today's 300 (more diffuse draws at equal beta; the untempered prior bounds
   it), and the fill cost (ruling 4).
6. The class-centring rule and the prior agree on the frame: cavger's centre-of-mass recentring of the
   class averages (members shifted accordingly) and the absolute prior on the particle are expressed in one
   frame, and the automask mass-centre shift and the cavg centring are one thing, not two. Settled after
   reading cavger, before Phase 2.
7. Masks. The centring stage runs with the widest mask, the mskdiam = 0 effect of `derive_image_settings`
   ((box - COSMSKHALFWIDTH) smpd), which also puts the 2D limits at lpstart 20, lpstop 8 and lpcen 30 Å and
   grows pftsz; a user mskdiam is ignored in that stage. At its exit a mask diameter is estimated from the
   class averages (automask2D's diameters, population-weighted, empty and junk classes excluded, a high
   quantile rather than the maximum) and becomes the next stage's parameter through `automask2D_mskdiam`;
   a user mskdiam is honoured from that stage on. No mid-run mask change: the stage boundary builds a new
   calculator anyway.
8. The centring stage exits on the shift statistics in addition to its usual criteria: the median
   per-iteration shift change over the sampled particles below lp/4 for two iterations.
9. Re-extraction (`reextract` with the recovered centres, shifts reset, the prior re-anchored) and Joe's
   data are later questions, after a working implementation.
10. Where the centring stage sits. Shifts are never searched against random or noise-initialised
    references. Today solve2D (six stages) sets trs = 0, center = no and Gaussian-regularised references in
    stage 1 and opens its window at stage 2 under snhc_smpl, one stage before prob at stage 3
    (`simple_solve2D_controller.f90`, `set_refine2D_stage_reference_policy`; the window
    `min(5, max(2, AHELIX_WIDTH / smpd_crop))` cropped pixels, about 12 Å, set in
    `simple_commanders_solve2D.f90`). solve3D (eight stages) sets trs = 0 in stages 1-2 and opens its window
    at stage 3 together with prob (`simple_solve3D_controller.f90`, `build_refine3D_stage_cfg`;
    `min(8, ...)` in `simple_estimate_ssnr.f90`). Neither window reaches the 20-30 Å the benchmark found.
    The centring stage is therefore solve2D's stage 2 and solve3D's stage 3: the map opens the 30 Å disk
    where the 12 Å window opens today, against references one stage old.

11. Under a switch (Hans, 2026-10-30), so the integration point is tested, not guessed. The map is off by
    default and every existing path is unchanged when it is off (identity on the `simulated_workflow`
    entries). One key, the prior's sigma in Å (working name `shift_sigma`; 0 or absent = the present window
    search), turns it on in whatever program carries it: radius 3 sigma, the absolute prior, the joint draw
    in the fill, the map as the coarse initializer at the window's radius elsewhere. refine2D and refine3D
    honour it on their command lines first (Phase 2a), so a stage is tested by handing the key to that stage
    alone; the solve2D/solve3D controllers adopt it as a stage parameter (Phase 2b) once Phase 3's synthetic
    recovery says which stage. The greedy arm of the comparison is a second, test-level key, not a user
    setting.

12. The map runs as the L-BFGS runs today, on the sampled subset of classes or projections (the refined
    candidates), in the centring stage too (Hans, 2026-10-30, on Phase 1's timing). The kernel is optimised first,
    as fast as it can be made, before Phase 2a builds on it; its exactness is kept (the identity tests), its
    speed measured against the Phase 1 kernel built from its commit (`tree`) on the same fixture.

What Phase 3 measures, in order: the fraction of particles that leave rejected classes for good ones; the
shift histogram against the Gaussian cryoSPARC measured, sigma included; class quality; after re-extraction,
the downstream solve3D. The synthetic recovery test against magnitude and SNR runs alongside.

## Open decisions

1. Closed (Hans, 2026-10-30): the 20 Å average is a mean magnitude, so a spread of picking errors, the case
   this fixes, not a systematic origin offset. Not every particle has a large shift: the large ones are
   probably the side views of the membrane protein, which the picker places with a systematic bias. Phase 3
   therefore reports the shift histogram per class as well as overall, since a view-dependent bias shows as
   a per-class offset and the frame question of ruling 6 then matters for those classes.
2. The centring frame (ruling 6), after reading cavger.

## Phases

### Phase 0: build and pass the mechanics (done 2026-10-30, Progress below)

The two new files were written without a compiler. Re-configure (new sources are globbed: `cmake -S . -B
build` or `./compile_debug.sh`), build, fix compile errors without changing the contract (the identity, the
truncation rule, the API), run `simple_test_exec test=unit_pftc_align2D3D suite=polar_shift_map`, then the
area suite and the fast gate; `scripts/check_test_registry.py`, `check_descr.py`, `check_flex_dag.py` (all three
pass on the written tree). Record in the journal and in Progress the orders kept and the maximum errors the
tests print.

Exit: the sub-suite, the area suite and the fast gate pass on the Debug and Release builds.

### Phase 1: the all-rotation kernel from the memo (done 2026-10-30, Progress below)

In a new submodule of the calculator, `src/main/pftc/simple_polarft_shiftmap.f90`: the all-rotation form from
`ft_ref` and `ft_ptcl_ctf` (full spectra by the twisted symmetry), per-thread work arrays sized by the team,
cc and euclid weights, the per-rotation prior factor applied on the grid, the grid maximum returned as (irot,
shift) in the optimizer's frame. The `polar_shiftmap` object holds the tables (built once per calculator from
box, pftsz, kfromto, the radius and lp). Tests in the tester: identity of the kernel with the mechanics class
to 1e-6 relative at random poses (cc and euclid), thread identity at 1 and 4 threads, and a timing against
`coarse_search_opt_angle` at equal coverage on the fixture sizes of `simple_pftc_inplane_tester`.

Exit: identity and timing table in Progress; `unit_pftc_align2D3D` and the gate pass.

### Phase 2: the switch in refine2D/refine3D, then the controllers (after open decision 2)

Phase 2a, behind the switch of ruling 11, no controller change. `pftc_shsrch_grad`: the map replaces
`coarse_search_opt_angle` when the radius exceeds the window; the L-BFGS limits become the magnitude step
about the grid start; the prior per ruling 2. Identity test: with the radius equal to `trs` the new path and
the old path give the same pose on the in-plane tester's fixtures (same objective, finer grid: equal or
better value, same final L-BFGS pose to 0.01 px). The tables (`eul_prob_tab2D`, `eul_prob_tab`): the fill scores every class at the stored shift as today and refines
the candidates; with the key on, each refined candidate's shift search starts from the map's joint draw of ruling
5 over its whole (irot, rho, alpha) grid (greedy on the test key), prior untempered, and the row holds the drawn
pose and its loss; the class draw over the rows is unchanged. Exit 2a adds the per-candidate cost of the map to
the fill's cost record (`record_fill_profile`).

Phase 2b, after Phase 3's synthetic part. solve2D's stage 2 and solve3D's stage 3, or whichever stage the
recovery test names, become the centring stage (ruling 10): the disk of 3 sigma, the absolute prior, the
box-wide mask, the shift-stabilisation exit and, in 2D, the automatic mask diameter handed to the next stage
(rulings 7-8). The later stages run the map at the window's radius on the refined candidates.

Exit 2a: the in-plane, table and strategy testers pass, the Phase 1 timing table extended to the fill, the
identity with the key off. Exit 2b: the controllers' testers pass; the synthetic recovery per stage in
Progress.

### Phase 3: recovery on data

Synthetic: particles simulated with picking errors drawn up to 30 Å (uniform direction, magnitudes 0-30 Å),
solve2D and solve3D with the map on, the recovered shifts against the truth per stage (the histogram, the
fraction within lp/4, the failures). Real: one of the benchmarked sets, SIMPLE's shift distribution against
cryoSPARC's on the same particles. Report in Progress: the answer first (does the map recover the shifts
cryoSPARC finds, at what cost per iteration), then the tables.

Exit: the report; the ruling on masks and re-extraction recorded.

### File disposition

| file | why | phase |
|---|---|---|
| `doc/implementation_notes/planned/polar_shift_search.md` | this note: progress, rulings as they are applied | each |
| `src/main/pftc/simple_polar_shiftmap.f90` | the mechanics class; phase 1 adds the table accessors of the all-rotation kernel | 0 1 |
| `src/main/pftc/simple_polar_shiftmap_tester.f90` | its unit tests; phase 1 adds the kernel's identity and timing | 0 1 2 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration (done); phase 2 if a strategy sub-suite is added | 0 2 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | suite help list (done) | 0 2 |
| `doc/policies/test_environment_policy.md` | the area table (done) | 0 2 |
| `src/main/pftc/simple_polarft_calc.f90` | the kernel's interface, the shift map and its per-thread work as components | 1 2 |
| `src/main/pftc/simple_polarft_shiftmap.f90` | the all-rotation kernel (new submodule) | 1 2 |
| `src/main/pftc/simple_polarft_core.f90` | `kill` releases the shift map | 1 |
| `src/main/pftc/simple_polarft_memo.f90` | only if the memo must keep what the kernel needs | 1 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | the coarse initializer | 2 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | the identity of the new and old coarse paths | 2 |
| `src/main/params/simple_parameters.f90` | the radius / switch | 2 |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | plumbing | 2 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | plumbing | 2 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the joint draw in the fill of the centring stage | 2 |
| `src/main/pftc/simple_polarft_corr.f90` | `gen_posterior_val` over the map's grid | 2 |
| `src/main/solve/simple_solve2D_controller.f90` | stage 2 as the centring stage: mask, prior, exit, mskdiam at exit | 2 |
| `src/main/solve/simple_solve3D_controller.f90` | stage 3: the disk in place of the window | 2 |
| `src/main/commanders/simple/simple_commanders_solve2D.f90` | the stage window | 2 |
| `doc/policies/3D/refine3D_policy.md` | the shift search paragraph | 2 3 |
| `doc/policies/2D/solve2D_policy.md` | the shift search paragraph | 2 3 |

Phase 2's UI files are added to the table with the ruling that names the key.

## Common rules

- Nothing is committed; the maintainer reviews the diff.
- The contract above is fixed: a compile fix may not change the identity, the truncation rule or the API; a
  test tolerance is not relaxed to pass (the mechanics are exact to double precision; a larger error is a bug).
- Fortran rules of `.github/skills/simple-modern-fortran/SKILL.md`: `use ..., only:`; THROW macros on one
  line; no goto; no `.and.` short-circuit on allocatable arrays; large locals `class(T), allocatable`;
  comments short.
- New OpenMP code is stress-tested at the machine's full thread count against one thread.

## Progress

### 2026-10-30: the mechanics and its test (Claude, Cowork sandbox, no compiler)

- Derived the identity on the calculator's conventions (`argtransf`, `gen_shmat*`, `rotate_ref_8`,
  `gen_corr_cc_for_rot_8_2`); the shift acts on the reference before the rotation, so Q is indexed by the
  reference and the particle index follows the rotation.
- Numerical verification in Python (double precision, mpmath Bessel values; box 64, pftsz 50, shells 2-16,
  magnitudes to 15 px): the per-rotation form against a line-by-line emulation of `gen_corr_cc_for_rot_8_2`
  at 60 random (irot, rho, alpha) incl. every rotation branch, max |cc difference| 4e-12; the closed-form ramp,
  2e-12 at scale 6750; rho = 0 equals the unshifted correlation; peak recovery of a shifted, rotated copy on a
  1 px / 5 degree grid: true rotation, shift within 0.46 px. The all-rotation form from the memoised half
  spectra (twisted symmetry, one inverse FFT per order and magnitude) against the same brute force at 200
  probes: 7e-8 absolute at scale 6750. Truncation table above from the same scripts.
- Written: `simple_polar_shiftmap.f90`, `simple_polar_shiftmap_tester.f90`, the sub-suite registration (test
  table, UI help, policy table); `check_test_registry.py`, `check_descr.py` pass.

### 2026-10-30: Phase 0 done (Mac Release build through the build bridge, gfortran 16.2.0)

- `cmake -S . -B build` (new sources), `make -C build -j12`: both new files compile without a warning;
  `make -C build install`.
- `simple_test_exec test=unit_pftc_align2D3D suite=polar_shift_map`: PASS, 28 checks, 0.05 s. Orders kept at
  magnitudes 0, 1, 2.5, 5, 8, 12 px (box 64, kmax 12): 0, 9, 13, 18, 23, 30. Closed form: max error 9.5e-9 on
  the grid at scale 3.1e3, 6.9e-10 off the grid. Against `gen_corr_for_rot_8`: max |cc difference| 3.7e-8.
  Peak recovery: true pose found, exact and at unit SNR (cc 0.703 at the true pose).
- One change on the first run: the tolerance of the comparison with `gen_corr_for_rot_8` was set at 1e-9 and
  the difference is 3.7e-8, because the calculator builds `argtransf` from single-precision angles (`polar`
  is `real(sp)`), not because of the expansion (the double-precision closed form holds to 1e-12 relative).
  The tolerance is 1e-6 with the reason in the tester; the contract is unchanged.
- The area suite: PASS, 10 sub-suites, 1303 assertions; the fast gate: 15/15, 43 s*proc. Debug build not run
  (it would replace the Release tree); to be run by the maintainer or the arun run before Phase 1.

### 2026-10-30: rulings, the stage schedule and the Debug gate (Claude, FWM-S116118, agent kit)

- Debug build of 519c67d2a through the agent kit (`~/agent_kit/agent_cycle.sh SIMPLE debug`, run 2, log in
  `~/agent_work/SIMPLE`): configure, build, parallel install and the fast gate 15/15 in about two minutes on 24
  cores. This closes Phase 0's open Debug item. The kit builds outside the checkout; `build/` is the
  maintainer's.
- Rulings 1-10 recorded from the design discussion; Decisions fixed amended for ruling 4 (the fill form in
  the centring stage); Phase 2 rewritten as the centring stage; the greedy ruling reopened with the
  three-degree draw proposal (open decision 1).
- The stage schedule read from the controllers (ruling 10): shifts open at solve2D stage 2 and solve3D
  stage 3, both windows about 12 Å.
- Hans agreed to the three-degree draw (ruling 5) and asked for a switch (ruling 11); open decision 1
  closed, Phase 2 split into 2a (refine2D/refine3D, the switch) and 2b (the controllers, after the synthetic
  recovery names the stage).
- Open decision 1 closed: mean magnitude (picking error, probably biased side views); Phase 3 adds the
  per-class shift histogram. Go for Phase 1 (Hans).

### 2026-10-30: Phase 1 done (Claude, agent kit runs 5-8, Debug and Release, gfortran 16.2.0)

- The kernel: `src/main/pftc/simple_polarft_shiftmap.f90`, a submodule of the calculator, with the shift map and
  its per-thread work as allocatable components (`new_shiftmap`, `kill_shiftmap`, `gen_shiftmap`,
  `shiftmap_cross`, `shiftmap_objfun_grid` with the optional prior, `shiftmap_peak`, `shiftmap_dims`,
  `shiftmap_shift`). The full-ring spectra come from the memoised half spectra by F[-m] = (-1)^m conj(F[m]);
  the shell products of each order in double precision, the Bessel-weighted shell sums per magnitude cast to
  single and transformed in batches of nk through the memo workspace's batched plan; the reference power per
  rotation from the memoised CTF2 and REF2 spectra as the Euclidean kernels form it. The mechanics class gained
  the table accessors the kernel needs (`shell_weights`, `eval_rotcoeffs`, `map_rotcoeffs`; the one-sided
  full-ring coefficients, C = Re(a_0)/2 + sum_{n>=1} Re(exp(-i n alpha) a_n)). No change to any existing path.
- Identity (box 64, pftsz 48, shells 2-11, 64 poses incl. the four rotation branches, random directions): the
  kernel's cross term against the mechanics class 9.3e-5 at scale 6.0e3 for cc and 7.0e-5 at 4.9e3 for euclid,
  1.5e-8 relative (tolerance 1e-6 of the Cauchy-Schwarz scale; the spectra are single precision); the grid's
  objective against `gen_corr_for_rot_8` 4.0e-8 (cc) and 1.9e-8 (euclid, as exp(-L)); the zero-magnitude row
  against the unshifted objective at every rotation 1.5e-8 (tolerance 1e-6). Four threads against one, eight
  particles: bit for bit. Sub-suite 41 checks, 0.15 s Release; the area suite and the fast gate 15/15 on both
  builds (Debug run 5, Release run 7).
- Timing (Release, box 64, 1.3 A/px, mask 60 A, euclid, one thread, ms per pair; the coarse grid is the 5 x 5
  lattice of `coarse_search_opt_angle` over +-5 px, 25 all-rotation evaluations of a shifted reference; the map's
  directions 2 z_max + 3):

  | lp (A) | nrots | shells | coarse 5x5 | map at the 5 px window (grid) | map at the 30 A disk (grid) |
  |---|---|---|---|---|---|
  | 8.0 | 144 | 9 | 0.178 | 0.096 (4 x 13) | 1.093 (15 x 47) |
  | 2.7 | 144 | 30 | 0.479 | 0.938 (10 x 33) | not a use case |

  Per grid point (every rotation included): coarse 7.1 and 19 us, the map 1.8 and 2.8 us at the window, 1.6 us at
  the disk, 4-7 times cheaper per point rather than the 10-20 estimated in Cost. The direction sums dominate
  (nalpha x nrots x sum over magnitudes of nharm): an inverse FFT over the directions and fewer directions at
  the inner magnitudes (spacing proportional to rho) would take most of it; not done, the exactness was the
  point of this phase.
- The cost claim of ruling 4, read from the table: the map at the disk costs about 1.1 ms per pair at this size
  against about 7 us for one per-class score of today's fill (one all-rotation evaluation at the stored shift),
  150 times; the obvious optimisations leave it above 30 times. Against every class the centring stage's fill
  would therefore cost 100 times today's; against the committed or a few candidate classes per particle it is
  negligible (200,000 particles x 1 ms / 24 threads = 8 s per iteration). Phase 2a's first decision is the
  class set the map runs on in the centring stage: every class, the top-k of the stored-shift posterior, or
  the drawn class alone (the picking offset is the particle's, not the class hypothesis's).
- Rulings 12 (Hans): the map runs on the sampled subset, as the L-BFGS does today; the kernel is optimised
  before Phase 2a. The Phase 1 state is the checkpoint the optimisation is measured against.
