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
- The map runs on the refined candidates of a particle (<= prob_athres references in 3D, the 2D refined set)
  as the coarse initializer of their shift search, before the windowed L-BFGS, not inside the table fill. A
  picking offset is a property of the particle, not of the (reference, rotation) hypothesis: the shift found
  on the committed candidate at one iteration is the stored shift every entry of the next iteration is scored
  at, so a mis-centred particle is scored wrongly for one iteration and then correctly.
- The L-BFGS window after the map is the magnitude step about the grid maximum (the basin), not the search
  radius.
- Same mechanics for 2D and 3D (one calculator).

## Decisions for the maintainer

1. The radius. Options: a new parameter (picking-error scale in Å, converted at the search sampling), or `trs`
   reinterpreted as the radius of the map with the L-BFGS window derived from lp. Recommendation: one
   quantity; tie the radius to the prior's width where the prior is on (radius = 3 sigma), which makes the
   global search and the windowed search one mechanism.
2. The prior's centre. The present prior is centred on the start (a random-walk prior: it restrains each step,
   not the walk). With a global search the picked position is the natural centre (an absolute prior, sigma at
   the picking-error scale), and it is the prior knowledge cryoSPARC's discrete search leans on in its first
   iterations. Recommendation: absolute centre when the map is on; measure both.
3. Masks and extraction. A particle 30 Å off centre inside a mask of radius mskdiam/2 is clipped and no
   correlation operator recovers what the mask removed. The centring stage needs the particle mask radius
   grown by the radius (or box taper only), and once centres are established the particles should be
   re-extracted (`reextract`), as cryoSPARC users do. Decide whether the centring stage is the first stage of
   solve2D/solve3D or a dedicated step before them.
4. The diagnostic before building further: on the benchmarked sets, is the 20 Å average the mean vector or the
   mean magnitude? A mean vector that large is a systematic origin offset (cryoSPARC recentres the volume on its
   centre of mass each iteration; SIMPLE's origin is the box centre) and needs no per-particle search; a
   near-zero mean with a 20 Å spread is picking error, the case this fixes. And SIMPLE's shift histogram on the
   same particles: pile-up at +-trs is the truncation signature.

## Phases

### Phase 0: build and pass the mechanics (done 2026-10-30, Progress below)

The two new files were written without a compiler. Re-configure (new sources are globbed: `cmake -S . -B
build` or `./compile_debug.sh`), build, fix compile errors without changing the contract (the identity, the
truncation rule, the API), run `simple_test_exec test=unit_pftc_align2D3D suite=polar_shift_map`, then the
area suite and the fast gate; `scripts/check_test_registry.py`, `check_descr.py`, `check_flex_dag.py` (all three
pass on the written tree). Record in the journal and in Progress the orders kept and the maximum errors the
tests print.

Exit: the sub-suite, the area suite and the fast gate pass on the Debug and Release builds.

### Phase 1: the all-rotation kernel from the memo

In a new submodule of the calculator, `src/main/pftc/simple_polarft_shiftmap.f90`: the all-rotation form from
`ft_ref` and `ft_ptcl_ctf` (full spectra by the twisted symmetry), per-thread work arrays sized by the team,
cc and euclid weights, the per-rotation prior factor applied on the grid, the grid maximum returned as (irot,
shift) in the optimizer's frame. The `polar_shiftmap` object holds the tables (built once per calculator from
box, pftsz, kfromto, the radius and lp). Tests in the tester: identity of the kernel with the mechanics class
to 1e-6 relative at random poses (cc and euclid), thread identity at 1 and 4 threads, and a timing against
`coarse_search_opt_angle` at equal coverage on the fixture sizes of `simple_pftc_inplane_tester`.

Exit: identity and timing table in Progress; `unit_pftc_align2D3D` and the gate pass.

### Phase 2: the coarse initializer (after rulings 1-2)

`pftc_shsrch_grad`: the map replaces `coarse_search_opt_angle` when the radius exceeds the window; the L-BFGS
limits become the magnitude step about the grid maximum; the prior per ruling 2. Identity test: with the
radius equal to `trs` the new path and the old path give the same pose on the in-plane tester's fixtures
(same objective, finer grid: equal or better value, same final L-BFGS pose to 0.01 px). Strategies: the
switch and the radius plumbed through `parameters` and the 2D/3D search objects; UI entries for the key.

Exit: the in-plane and strategy testers pass; refine2D/refine3D on `simulated_workflow` entries unchanged
when the radius equals the window.

### Phase 3: recovery on data (after rulings 3-4)

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
| `src/main/pftc/simple_polar_shiftmap.f90` | the mechanics class | 0 1 |
| `src/main/pftc/simple_polar_shiftmap_tester.f90` | its unit tests; phase 1 adds the kernel's identity and timing | 0 1 2 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration (done); phase 2 if a strategy sub-suite is added | 0 2 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | suite help list (done) | 0 2 |
| `doc/policies/test_environment_policy.md` | the area table (done) | 0 2 |
| `src/main/pftc/simple_polarft_calc.f90` | the kernel's interface | 1 2 |
| `src/main/pftc/simple_polarft_shiftmap.f90` | the all-rotation kernel (new submodule) | 1 2 |
| `src/main/pftc/simple_polarft_memo.f90` | only if the memo must keep what the kernel needs | 1 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | the coarse initializer | 2 |
| `src/main/pftc/simple_pftc_inplane_tester.f90` | the identity of the new and old coarse paths | 2 |
| `src/main/params/simple_parameters.f90` | the radius / switch | 2 |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | plumbing | 2 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | plumbing | 2 |
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
