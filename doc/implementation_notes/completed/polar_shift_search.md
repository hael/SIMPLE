# Polar shift map: a global shift search on the polar Fourier transform

Started 2026-10-30. One living note: the contract, implementation plan, open decisions and progress record.

## Status: closed 2026-11-03

Completed. The polar shift map runs by default (`shift_sigma=10` A, objfun=euclid) in every polar mode with a
shift search: the four probability tables and, since ruling 20, the classic SHC, greedy and neighbourhood
strategies. Each particle is centred as it needs (ruling 19: the disk when not yet centred in the workflow or
wide, the window otherwise), the map initialises every polish, and the polish starts at the map's cell. Stage 1
of solve2D and solve3D's first two stages stay shift-free by design.

Hans's validation (2026-11-03), after the cell-start polish and the classic strategies: "best porF so far in
10/10. Same quality straight out of solve3D as cryo-sparc after NU refinement"; solve2D re-run on the small
membrane proteins: "It looks amazing". The closing changes (ruling 21) are in the last progress entry.

Left open when the note closed:

| Item | State |
|---|---|
| The Laplace width | Its own planned note: `../planned/polar_shift_laplace_width.md` (Hans's call) |
| The centring frame (open decision 2) | Closed by ruling 22 (2026-11-03): the prior moves with the particle |
| The stochastic 2D classic modes (`snhc`, `snhc_smpl`, the stream's pool) on a fixture that separates them from run-to-run spread | Deferred (Hans: not now) |
| The dynamic-state regressions the review asked for | Deferred to the Laplace width's change (Hans: later) |
| Phase 3's written exit (per-class shift histograms against cryoSPARC's, re-extraction, downstream solve3D) | Hans's data; the 2D runs on the small membrane proteins are in |
| The automatic mask's parity (5.63 against 4.29 A on beta-gal) | Accepted as good enough (Hans: the small beta-gal set says little) |
| The box's mask in a standalone refine2D's two shift-free iterations | Not now (Hans) |

The dated inventories, superseded schedules and progress below are retained as history.

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
- The kernel is the data term only. The sampler assembles the untempered prior and the polar cell measure itself
  (review, sections 2 and 6, adopted 2026-10-31): the shift search's prior in normalised loss units,
  1 / (2 sigma^2 T), is tempered along with the loss when the row is exponentiated at beta T, so the map's grid
  does not carry it. The grid gives the loss (or correlation), the shift of every cell in the optimizer's frame
  (`shiftmap_shift`) and the cell's translation area (`shiftmap_cell_area`: rings between the midpoints of the
  magnitudes, the central disk divided across the directions that all name the origin).
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
    speed measured against the Phase 1 kernel built from its commit (`tree`) on the same fixture. Done 2026-10-31
    (Progress): 2.3 times at the disk, 1.9 at the windows, about ten times cheaper per grid point than the coarse
    lattice.

13. Measured, not hard-coded (Hans, 2026-10-31, on the cost plan): the centring stage's length and the screen's
    breadth follow what the fill measures, not constants. Two quantities, both free in the fill: the posterior mass
    beyond the window's radius under a particle's map (the disk is needed while sampled particles carry mass out
    there; the stage ends when almost none do, within its budget: the exit ruling 8 reached for), and the origin
    rank of the class the screen puts first (the screen's breadth for the next iteration is the measured rank
    quantile of the last, every class at the first iteration; the tolerance, a confidence level, is the only
    parameter). Phase 2c builds the measures into the fill and the cost record first, then the exit and the
    adaptive breadth on them.
14. The map initialises the classic search (Hans, 2026-10-31): the map's best cell over all rotations, refined
    below the grid (Newton in alpha, a parabola in rho), is the L-BFGS start with the bounds one grid step about
    it, in place of the lattice seed and the seed bank; every mode is reviewed for it, solve2D's stages, refine2D's
    and refine3D's modes, solve3D's stages, single and multi-state alike (Phase 4, from an inventory of today's
    shift initialisations). The polish in the Fourier-Bessel domain, on the map's harmonic correlations
    T_n(k, irot) with analytic gradients in rho and alpha, about five times cheaper per evaluation, follows where
    the map is built anyway (Phase 5; not pursued 2026-11-03, the measurement in the entry of that date).

What Phase 3 measures, in order: the fraction of particles that leave rejected classes for good ones; the
shift histogram against the Gaussian cryoSPARC measured, sigma included; class quality; after re-extraction,
the downstream solve3D. The synthetic recovery test against magnitude and SNR runs alongside.

15. `topk=3` by default (Hans, 2026-11-01, after the 3D sweep). The weighted insertion therefore runs under the key,
    which closed R3: every refined entry carries its pose's own loss (`prob_refined_entry%loss`, the row's distance
    under the window search and the drawn cell's data term under the map), and the extras are weighed by those
    losses on both paths (`select_extras`); the refusal of the key with `topk > 1` is gone.

16. The consolidation (Hans, 2026-11-01: "hack after hack"; "go"): the 2D centring stage is coarse-to-fine, the coarse
    map over the disk at a particle's likely classes (its posterior at the origin), the L-BFGS polish at full
    resolution from the coarse optimum, no full-band disk map, rows the polished losses as the window path's (no
    evidence-based draw); a particle whose last move stayed inside the next stage's window takes the window path
    (`shincr`, slot 54), the stage ends when few still need the disk. The rank record, `nscreen`, the mass measure
    and its exit, the screen-posterior bound and the 2D draw are gone; the key reaches the centring stage alone.
    The 3D dense fill keeps the full map for now (its cost on data is the next measurement).

17. solve2D runs the key by default for now (Hans, 2026-11-01, so the tests do not forget it): the commander sets
    `shift_sigma` to `SHIFT_SIGMA_DEFAULT` (10 A) when none is given, like its `refine`; `shift_sigma=0` is the
    control (the window search in every stage). The parameter's default stays 0: a refine2D outside solve2D (the
    stream's pool, a standalone run) keeps the window search, since its sparse table refuses the key. Every solve2D
    caller gets the default, the stream's initial analysis and the pyramid included. Development pauses until
    solve2D is tested on data with the drain reset in place (the earlier bench numbers predate that fix).
    Extended 2026-11-02 (Hans, after the beta-gal 3D test: "shift_sigma=10 by default across the board, including
    solve3D_cavgs"): the parameter's own default is 10 A; the key is on wherever a mode can carry it (a probability
    table with a shift search under objfun=euclid) and silently off elsewhere, where an explicit key is refused;
    the commanders set no default of their own.

18. The weighted insertion's extras take the continuous in-plane route with the committed pose (Hans, 2026-11-01:
    "we ought to inpl_cont-refine all of them, or none of them"; "I do not think it will make a huge difference"):
    each from its own cell and shift through the same local joint solve, the weights the table's, the route status
    the committed pose's; under the key nothing is polished, as before (Phase 5's Fourier-Bessel polish, not pursued).

19. Dynamic centring (Hans, 2026-11-02; the proposal and its rulings in `polar_shift_search_dynamic_centring.md`; built the same day, Progress below):
    the sampling strategies inject particles the centring stage never sees, so the centring becomes a per-particle
    state across every stage and table, 2D and 3D: the sentinel per workflow and on injection, the disk particle
    prepared under the box's mask, the move at the bound and the posterior width combined into one prior width that
    chooses the next box at three standard deviations (the pick prior the only prior in the likelihood), the
    neighbourhood table's disk over the subspace pool, `trs_window` and the drain exits retired.
20. The classic strategies take the table's structure, not a seed (Hans, 2026-11-02): under the key `inpl_srch_first`
    goes, in every state count. A seed at the previous reference is optimal for that one state only; the states are
    not docked and nothing records which projection of each other state was optimal last time, which is why the
    classic multi-state search never seeded. As in the tables, the candidates are scored at the stored shift, the
    best few get a map each (the window, or the disk by the particle's width) and the best cell is polished in its box.
    Built 2026-11-03: the pool is every candidate the mode scored and the previous one, window or disk alike (a
    first build that mapped only SHC's committed class in the window lost the 2D fixture, Progress below).
21. The close (Hans, 2026-11-03): the targets of the three priors are stated once. The coarse draw (the
    initialiser) carries the pick prior at a particle's first centring; every polish (the cell's, the classic
    one, the continuous route) is the data's alone, so the window search's increment prior (`shift_prior`) is
    off under the key (the review's T3); the weighted insertion weighs poses by the data term with the map's prior
    at the physical temperature, in the dense and the neighbourhood tables alike (`insert_pose_loss`, the review's
    T2; the prior's centre is ruling 22's). The PCG backend takes the weighted insertion's poses ("fix it").
    solve2D's terminal greedy pass carries the key ("make the change"). The Laplace width goes to its own planned
    note; the stochastic 2D fixture, the dynamic-state regressions and the shift-free iterations' mask are
    deferred; the automatic mask is good enough for now.
22. The prior moves with the particle (Hans, 2026-11-03, closing open decision 2): "Each particle is shifted to its
    previous rotational origin offset in prepimg4align and then moved with the center of the class average so
    having the prior defined at 0,0 must be holding us back". The map's prior is on the cell's shift, the increment
    from the stored shift (which carries the search's moves and the class centring's), not on the total
    displacement from the pick: in the first centring's draw and in the insertion's loss alike. In 3D the first
    centring no longer pulls a particle from its 2D shift back towards its pick.

The 2026-11-01 follow-up review is retained in `polar_shift_search_review.md`, section 13.
Its stage-entry, drain and cost findings were corrected or superseded by the subsequent consolidation and
dynamic state. The coarse fallback now preserves the better origin row; dense insertion now carries an
absolute prior loss. Section 14 rechecks those dispositions and records the remaining path differences.

The maintainer's stage-1 clarification remains in force: no shift search there because it causes class
collapse. The former centring budget was five iterations in stage 2 (global 6-10), and the proposal to
extend it through stage 3 is superseded by the per-particle state of ruling 19, as Open decision 4 records.

## Open decisions

1. Closed (Hans, 2026-10-30): the 20 Å average is a mean magnitude, so a spread of picking errors, the case
   this fixes, not a systematic origin offset. Not every particle has a large shift: the large ones are
   probably the side views of the membrane protein, which the picker places with a systematic bias. Phase 3
   therefore reports the shift histogram per class as well as overall, since a view-dependent bias shows as
   a per-class offset and the frame question of ruling 6 then matters for those classes.
2. The centring frame (ruling 6), after reading cavger. Closed by ruling 22 (2026-11-03): the prior moves with the
   particle. Hans had asked for more information first. In short: with `center=yes` the 2D search recentres each class average before each
   iteration (`calc_2Dref_offset`: the mass centre by default) and moves every member's stored shift by the same
   offset (`add_shift2class`), while the pick prior is on the total displacement from the pick. For a view the
   picker places off-centre the two disagree; the prior acts at a particle's first centring in 2D and in every
   3D insertion's loss.
3. Closed (2026-10-31, revised by the re-review's R1 the same day): the table has one frame contract, every
   increment native (the particle's frame). The window path's `minimize(sh_rot=.true.)` already returns that; the
   map fills convert the drawn cell at the drawn rotation when they record it; every consumer adds the increment
   unconverted and builds the optimizer-frame shift only for the objective. My morning fix (`extra_pose` converting
   the extras' increments) had assumed the map's frame for the window path's entries too and double-rotated them;
   it is retired and `test_table_increment_axes` (the refine3D in-plane state suite) pins the contract through
   the optimizer itself. The weighted insertion (topk > 1) is refused under the key until R3 is done.
4. Superseded 2026-11-02 by the dynamic centring proposal
   (`polar_shift_search_dynamic_centring.md`, Hans's finding: the sampling
   strategies inject particles the centring stage never sees; the centring
   becomes a per-particle state across the stages, 2D and 3D). The original
   text: extend the faster centring search
   across stages 2 and 3, giving a nominal budget of ten centring iterations
   (five per stage), with stage 1 remaining shift-free. This is a proposal,
   not an adopted schedule. A coherent extension would select the dense
   coarse-to-fine path in both stages, initialize `shincr` on entry to stage
   2 and carry it into stage 3, preserving the drain instead of restarting
   every particle's disk search. Both stages' `trs_window` must describe the
   ordinary window after centring (stage 4), since stage 3 would itself have
   the wide disk range. Keep the frequency ladder explicit and decide the
   mask-estimate/adoption boundary across these stages. Compare the corrected
   five-iteration baseline with this extension using paired seeds: shift
   recovery, class populations/collapse, class quality and total wall time.

## Review findings and their disposition

`polar_shift_search_review.md` (2026-10-30, commit 00dc4ec55) reviewed the mathematics, the posterior
construction and the cost. Its findings and what the plan does with them:

| finding | disposition |
|---|---|
| P1 equal weights on the polar grid are the wrong measure (section 2) | kernel: `shiftmap_cell_area` (2026-10-31); Phase 2a: the sampler's log cell area in every cell's log mass |
| P1 a pose draw per class, then the class from that pose's loss, is not a joint draw (section 3) | Phase 2a requirement: the class is drawn from the evidence log Z_c = logsumexp over the candidate's map of (-beta T L - prior + log pi_c + log area), one conditional pose kept per candidate; evidence and the sampled loss in separate fields |
| P1 L-BFGS after the draw changes the distribution (section 4) | ruled (Hans, 2026-10-31): the existing table paths keep their polish; the new path is strict, the drawn pose is not polished; the coupling to continuous pose search is a later question |
| P1 screening the candidates at the current shift can exclude the right class (section 5) | ruled (Hans, 2026-10-31): a low-resolution map against every class as the screen; exploration candidates kept as the alternative; Phase 3 measures retained mass against the exhaustive joint evaluation on a subset |
| P1 the submodule missing from a2c1b1f86 (section 7) | committed in 9c2fcdaf9 |
| P2 the cost must count the mapped references (section 8) | the kernel's optimisation (ruling 12); Phase 2a's cost record counts candidates; 55 references per particle at a 10 degree bound in 3D is the figure to beat |
| P2 prior, frame, precision, calibration (section 6) | kernel: the data term only, the error reported in the sampler's units (T dL); fixture with angular CTF, distinct halves and spectra (2026-10-31); frame and crop validation in Phase 2a |
| sigma 10 A is a modelling choice; a 20 A mean radius implies sigma about 16 A | sigma is the parameter of the switch (ruling 11), a setting, not a constant |

The integration re-review of 2026-10-31 (the review's section 12, on `04f6f0cb7`..`4114f9f4a` plus the uncommitted
tree) and the dispositions, all applied the same day:

| finding | disposition |
|---|---|
| R1 P1: `extra_pose` double-rotates the window path's entries, which `minimize(sh_rot=.true.)` already converts; the table had two frame contracts | one contract: every table increment is native. The map fills convert the drawn cell at the drawn rotation when they record it (2D dense and sparse, 3D dense); the consumers add it unconverted and score the pose in the optimizer's frame (`shift_ptcl2ref` for the objective only); `extra_pose` retired, the extras go through `candidate_pose`; `test_table_increment_axes` pins the producer's frame through the optimizer (sh_rot true against false converted at the found rotation) and the consumer's addition |
| R2 P1: the cone compares the map's reference-frame cells as if native | follows from R1: the entries the cone sees are native |
| R3 P1: the weighted insertion weighs the evidence row as a pose loss | done 2026-11-01 (ruling 15): every refined entry carries its pose's own loss beside the row; the extras are weighed by those losses on both paths; the refusal is gone |
| R4 P2: a disk below one step reaches outside its radius | `build_map`: below one step the grid is the origin and the rim at the radius; `test_disk_within_radius` (below, at, above a step) |
| R5 P2: the scheduled Cartesian polish inherits the key and hits the scope guard | the polish pass's command and job description drop `shift_sigma` |
| R6 P2: `prob_post_diag=yes` under the map writes default rows | refused explicitly in the 3D fill |
| the 3D screen's pool of 4 projs_ns cannot rescue a direction outside it | the measured breadth of Phase 2c (the rank quantile) goes to the 3D fill with its measures; the first iteration's pool stays bounded in 3D |

The follow-up review of 2026-11-01 (the review's section 13, on `c980026e5` plus the consolidation) and the dispositions,
applied the same day:

| finding | disposition |
|---|---|
| S1 P1: the drain state (`shincr`) was reset only for runs starting at iteration 1; solve2D's centring stage starts at 6, so on entry every particle counted as done with the disk (the fixture: 0 disk searches at the first centring iteration, in both modes) | the refine2D master resets `shincr` at the start of every run under the key, whatever its first iteration's number, and writes the segment in the distributed mode; the master prints the disk set it is about to search (`DISK SEARCHED BY n OF m`); the bench re-measured with the disk searched at the entry |
| S2 P1: the 3D insertion's pose loss was the drawn cell's data term alone, dropping the absolute prior | `prob_refined_entry%loss` under the map is the posterior loss, the data term with the absolute prior on the total displacement at the physical temperature |
| S3 P2: the class-overlap verdict could end the centring stage while the disk was still needed | under the key the drain is the stage's convergence (the class verdict does not end it; the budget does); the convergence tester covers converged classes with a non-drained field |
| S4 P2: when the box's optimum did not beat the optimizer's origin baseline, the coarse point overwrote a better origin row | the coarse optimum, scored at full resolution, is the pose only if it beats the origin's row; else that row stands |
| S5 P2: the drained particles' window objects were not in the cost record | both object sets' counters reach the record (`grad_win`) |
| latent: `temp_mass` in `shiftmap_draw` overwrote the draw's normaliser | the retired diagnostic (rho_out, mass_out, temp_mass) removed from the sampler and its test |
| the origin-posterior screen quantile does not bound the posterior over translated candidates; retained mass and coarse-grid convergence unmeasured | open: the fixture's 95.8 % (its refined count of one) is the only measure so far |

The review's section 14 (2026-11-02) and the dispositions, all applied 2026-11-03:

| finding | disposition |
|---|---|
| T1 the continuous route projects a disk seed into the window and scores the projected pose | the joint domain contains the particle's box and the incoming shift (`shift_box_radius`), 2D and 3D, committed pose and extras; `test_continuous_seed_beyond_window` fails at all five starts with the old bound |
| T2 the neighbourhood insertion lacks the dense table's physical pose-prior loss | one definition, `insert_pose_loss`, in both tables; `test_insert_pose_loss` reproduces the review's 0.603483 / 0.396517 |
| T3 the mapped sparse and window polish keep the legacy increment prior | the increment prior is off under the key (`l_shift_prior`, the continuous route included); `test_key_derivation` |
| the dynamic-state regressions | deferred with the Laplace width (Hans) |
| the Laplace width's qualifications (fractional angle, conditional against marginal, indefinite curvature) | carried into `../planned/polar_shift_laplace_width.md` |

The review's section 15 (2026-11-03, the end-to-end audit) and the dispositions, applied the same evening (Hans: "Go
ahead and do all this"):

| finding | disposition |
|---|---|
| M1 the Euclidean score row and the continuous evaluator drop the angular Nyquist bin the map and the direct residual carry | every Euclidean kernel sums to `pftsz+1` (`gen_raw_euclid_vals_impl`, `gen_euclid_crvec`, `gen_raw_euclid_grad_at_angle`), as the cc kernels, the bank and the batched rows already did; `test_angular_nyquist`; the bank tests' tolerances tightened from 5e-2 to 1e-4 |
| M2 a disk particle's reconstruction observation depends on whether its record is kept (box-mask record, ordinary-mask reread) | the preparation records its mask per particle (`builder%obs_mask`, `note_obs_mask`), every reread in the process reads it (`prep_rec_observation_ptcl`, the gridding full-size path, the 2D class sums from full-size images); `test_rec_observation_mask` |
| M3 Bessel orders beyond the ring's length index outside the ring (a write past the array in Release) | the spectral index wraps modulo the ring (`gen_shiftmap`); `test_orders_beyond_ring` |
| M4 an excluded row can zero every eligible insertion weight | the extras' normaliser runs over the committed and the eligible rows (`top_posterior_extras`); the review's example in `test_top_posterior_extras` |
| the initialiser's cell-mass mode | measured against the density mode (the progress entry) |
| the width's dimensional convention | stated in the Laplace note (the coarse width is an RMS radius, sqrt(2) per-axis sigmas) |

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

### Phase 2: the switch in refine2D/refine3D, then the controllers (after open decision 2; 2a in progress, Progress below)

Phase 2a, behind the switch of ruling 11, no controller change. `pftc_shsrch_grad`: the map replaces
`coarse_search_opt_angle` when the radius exceeds the window; the L-BFGS limits become the magnitude step
about the grid start; the prior per ruling 2. Identity test: with the radius equal to `trs` the new path and
the old path give the same pose on the in-plane tester's fixtures (same objective, finer grid: equal or
better value, same final L-BFGS pose to 0.01 px). The tables (`eul_prob_tab2D`, `eul_prob_tab`): the fill scores every class at the stored shift as today and refines
the candidates; with the key on, each refined candidate gets its map, and from it, in one pass, its evidence
log Z_c = logsumexp over the grid of (-beta T L - |u|^2 / (2 sigma^2) + log area) with the untempered prior on the
absolute displacement from the pick and the cell measure, and one conditional pose draw of ruling 5 (greedy on
the test key). The row holds the evidence, the drawn pose and the drawn pose's loss in separate fields; the class
draw is from the evidence (review, section 3), the candidates' coverage of large offsets is established (review,
section 5), and the drawn pose is not polished (review, section 4; the polish is the deterministic arm). Exit 2a
adds the per-candidate cost of the map to the fill's cost record (`record_fill_profile`).

Phase 2b (done 2026-10-31, Progress: solve2D, then solve3D and the key's scope). solve2D's stage 2 and solve3D's stage 3, or whichever stage the
recovery test names, become the centring stage (ruling 10): the disk of 3 sigma, the absolute prior, the
box-wide mask, the shift-stabilisation exit and, in 2D, the automatic mask diameter handed to the next stage
(rulings 7-8). The later stages run the map at the window's radius on the refined candidates.

Exit 2a: the in-plane, table and strategy testers pass, the Phase 1 timing table extended to the fill, the
identity with the key off. Exit 2b: the controllers' testers pass; the synthetic recovery per stage in
Progress.

### Phase 2c: the adaptive centring stage (ruling 13)

The fill measures, per particle, the posterior mass beyond the window's radius under the map of the row's best
candidate and the origin rank of the screen's choice; both travel in the table's refined entry to the master.
The consumer writes the mass to the particle field; `check_conv2D` ends the centring stage once the fraction of
sampled particles with mass beyond the window above a tolerance falls below a tolerance for two iterations. The
master sets the next iteration's screen breadth from the rank quantile of the last (`nscreen`, an internal
control: every class at the first iteration); the tables screen the origin's best that many. The window's radius
reaches the table from the controller. The cost record prints both distributions. Exit: the fixtures' recovery
unchanged, the beta-gal bench's centring stage ending by itself and its screens shrinking after the first
iteration, the window stages' screens collapsing to about twice the refined count.

### 2026-11-01: the re-review chain green, the automatic mask's estimators, topk=3 and the insertion under the map (Claude, runs 463-502)

- The chain after the review fixes (the frame test rebuilt on a white reference and a rotated, shifted particle, the
  calculator's rotation helper being private): Debug and Release gates 15/15; suites in-plane state 72, polar shift
  map 73, 2D I/O 42, pose statistics 135; the frame contract on the fixtures: 2D map from true references 100 %
  (0.75 px), 3D map3 100 % (0.90 px; the 5 px window 90.5 %, the 16 px window 99 %).
- Bench with the user's mask and the evidence-based stage 1 (key2c): 239 s (stage 1 87 s, the centring stage 73 s),
  best class 4.29 A, 67 % of the particles at 8 A or better, the best set so far, but in 38 populated classes of 90
  (the control 65, the key with the stage-1 shift 57): the evidence draw concentrates the particles; Hans to judge on
  heterogeneous data.
- The automatic mask: the pre-estimate from the unaligned mean is unreliable, 110 A under Otsu and 119 A as the extent
  above 3 sigma of the mean's noise for the 180 A molecule (the rotational average's tail is weak, the class averages'
  edges are what automask2D measures), and stage 1 under that mask cost the run (4.66 A, 38-53 % at 8 A); a second
  defect, the ladder was not re-derived when the estimate was adopted (stages 3-5 ran at 8.2/7.1/6.0 A from the 110 A
  ladder). Now: no pre-estimate (stage 1 under the box's mask), the class-average estimate adopted at the exits of
  stage 1 and of the centring stage with the ladder re-derived each time (`rederive_lplims`). Measured: stage 1 at
  20 A under the box's mask, the estimates 165 and 183 A, the ladder 13.8 then 10.9/9.4/8.0 A; best class 5.63 A,
  55 % at 8 A or better in 42 classes, 190 s, against 4.29 A and 67 % with `mskdiam=180`: the automatic mode works
  but is not at parity, and the gap is stage 1's (20 A under the box's mask against 15 A under the molecule's);
  an automatic run needs the mask before stage 1, which no estimate from unaligned particles gave. Give `mskdiam`
  for the sweep.
- Ruling 15: `topk=3` by default; R3 done (the pose loss beside the row, one insertion distribution). At the new
  default the fixtures hold: 2D map 100 % (0.74 px), 3D map3 100 % (0.88 px; the 5 px window 91 %, 16 px 99 %), the
  insertion active under the map from the second iteration (the cone at 9.5 A; on the fixture almost no second
  refined entry falls inside it). The gates found two tester assumptions of the old default: the convergence tester's
  single-draw tail test pins `topk=1` now, and the 2D part reader's refined-entry byte count follows the widened
  record (20 bytes); the 3D table suite run standalone on the full team exposed a latent sizing bug in
  `select_extras` (per-thread masks sized by `nthr_glob`, not the running team), fixed per the OpenMP rule.
- Hans's bgal run on the pushed-candidate tree collapsed the classes at stage 2: the evidence-only stage 1 (no shift
  assigned) is the cause (Phase 2d, revised above); the mode removed, stage 1 under the key assigns the drawn cell;
  the bench re-measured in the next run (populated classes per iteration, the class quality).
- Hans on TRPV1: the centring stage's first iteration took 141 s against 4.6 s for the window iteration before it
  (every class screened for every particle, the one breadth still hard-coded), the class averages much improved.
  The second iteration cost the same (132 s): the rank quantile's breadth stays near every class while any particle
  is still off. Now each particle's pool at every iteration is measured from its own origin posterior (the classes
  holding SCREEN_QUANTILE = 0.99 of the mass at the draw's temperature, capped by the rank breadth when set, at least
  twice the refined count), both 2D fills, and the full disk maps go to the screen's candidates holding that mass
  of the screen's posterior, at most the refined count. Hans (2026-11-01): "the centring stage is VERY expensive
  when the number of particles scale up"; "a grid/L-BFGS hybrid?"; "hack after hack"; "go" on the consolidation
  (ruling 16). Built the same day: `simple_eul_prob_tab2D` rewritten (the dense fill coarse-to-fine with the drain,
  the sparse fill without the map; the part files back to seed, flag and temperature per particle), `shincr` in
  slot 54, `DISK_NEED_FRAC` the exit's constant, the refine2D master's reset of `shincr`, the key to the centring
  stage alone, `nscreen`/`shmassout`/the profile's rank and mass counts removed, the 2D UI without `shift_draw` and
  with `trs_window` for standalone runs. Validated (runs 562-577): Debug and Release gates 15/15; the 2D fixture from
  true references 95.8 % within lp/4 at a median error of 0.21 px (the full map: 100 % at 0.74 px; the fixture's five
  classes give a refined count of one, so the class is the coarse map's choice alone, and the 4 % are its misses,
  which the polish in a lattice box about the coarse optimum did not change); the bench with `mskdiam=180`: no
  collapse at the centring stage (90 populated through its four iterations, the decline from stage 3 as the
  control's, 61 at the end against the control's 65), best class 4.29 A, 24 % at 6.5 A or better, 57 % at 8 A, the
  centring stage 20.8 s (fills 2.8-3.6 s against 12-13 s with the full map; the disk needed by 13, 19, 9, 4 % of the
  particles), the run 111.5 s against the control's 100.6 s. The 3D dense fill keeps the full map.
- The follow-up review (section 13) applied the same day, the Debug gate 15/15 (run 579): the drain state reset at
  the start of every run under the key (before, only at iteration 1, so solve2D's centring stage, starting at
  iteration 6, searched no disk on entry: the bench numbers above were measured with the disk barely used), the
  drain the centring stage's only exit under the key, the posterior pose loss in the 3D insertion, the origin row
  kept when the polish loses, the window searches in the cost record, the retired mass diagnostic removed. Then
  ruling 17: solve2D defaults the key (`SHIFT_SIGMA_DEFAULT`), `shift_sigma=0` the control; Hans tests on data
  (bgal, TRPV1, a badly picked set; key against control, mask given against estimated) before further development.
  Ruling 18 built the same day (`resolve_extra_poses` with the polish, `refine_extra_pose_continuously`).

### 2026-11-02: the dynamic centring built (Claude, autonomous run, agent kit runs 692-)

Ruling 19 as built, in the plan's steps:

- The 2D particle disk cache removed (`simple_ptcl_cache` and its sites in the matcher, the batch preparation,
  the table worker, the masters, the stream's initial analysis, the state-assignment route, the parameters and the
  UI; the exit hook in `simple_defs`). The cropped-particle path in memory stays (`prepimg4align_cropped`, the
  records). The flex PCA plane cache kept `cache_dir` and got its own location helpers (`plane_cache_dir`,
  `plane_cache_run_token`, `plane_cache_space_ok`), gated by `cache_dir` given instead of `cache=yes`.
- The state: `shsig` (slot 55, `I_SHSIG`) beside `shincr`; the sentinel set by the refine2D and refine3D masters
  at a run's first iteration (every workflow's first stage starts there, whatever the key's state in that run) and
  by the stream's pool for the particles it appends (`append_sets`). The width: `shiftmap_draw` measures the
  posterior's standard deviation at the physical temperature (the mean in the particle's axes from each rotation's
  sum rotated once; the second moment is rotation-invariant), carried by `prob_candidate%sdev`,
  `prob_refined_entry%sdev` (the 2D part's refined entry 24 bytes) and `ptcl_ref%sdev` to the consumers;
  `record_shift_state` (utils) writes the move and `next_shift_width` (the width, or the previous one without a
  width, at least the radius searched when the move reached its bound within half a map step).
- The rule in the four fills (`shift_needs_disk`: the sentinel, or three widths beyond `trs`): the disk of
  `disk_radius_px` (three sigma, bounded by the box) under the prior in the disk's own search object per thread,
  the window of `trs` otherwise; the sparse 2D fill scores every class for a disk particle, then screens; the
  neighbourhood table screens the state's evaluated candidates (the state mode's subspace search in front, the
  geom mode's neighbourhood as it is) with the disk map and polishes the best in the disk object.
- The preparation: `ptcl_on_disk` (2Dprep) reads the same field and rule; `prepimg4align` and the cropped
  preparation take the box's noise mask (`build%lmsk_box`) and soft mask (`msk_box_crop`) for a disk particle,
  the records too.
- Retired: the centring stage (`FIRST_SHIFT_STAGE` keeps stage 2 on the dense table under the key and under the
  box's mask only for an automatic run's mask estimate), the stage-3 raise of `trs` in solve3D, `trs_window`,
  the drain exits of `check_conv2D`/`check_conv3D` and their tests, `DISK_NEED_FRAC`. The drivers print the disk
  set by cause and the widths' median and 90th percentile (`print_disk_set`).
- Debug gates 15/15 at each step (runs 695, 697). The first Release chain (runs 698-747) found two defects: the
  derivation had made the key depend on the run's `refine`, so the solve3D master (default `shc`) carried no key
  and refused an explicit one; and a re-centring under the pick prior pulled centred particles back towards their
  picks (the draw tempers the data term by `prob_beta` 0.3, not the prior), the solve2D fixture falling to 40 %
  within lp/4 against the window's 79.5 %. Fixed (run 749): `l_shift_map` is the key's presence (euclid), which the
  masters pass on, `l_centring` a run that centres (a table with a shift search), which the fills, the preparation
  and the consumers read; the pick prior applies at a particle's first centring only, a re-centring is the data's
  alone; the bound test is a tenth of the radius or half a map step where the map is finer. The second Release
  chain (runs 750-799): the gate 15/15, the in-plane state suite 98/98 and the 3D table suite 496/496 at four
  threads; the refine2D fixture from true references 100 % within lp/4 at 0.19 px mean (the window of 5 px 100 %
  at 0.30 px, of 16 px 100 % at 0.19 px: every polish map-initialised now); the refine3D fixture 100 % at 0.53 px
  (the 5 px window 96.2 %); solve2D from the random start at full sampling, the key 98.8 % at 1.33 px against the
  window's 72.8 % (the stage design gave 96.8 %); at a quarter sampled per iteration (`nsample=100`, the case the
  stage design missed: three quarters of the particles first sampled after what was the centring stage) the key
  61.5 % against 50.0 %; solve3D de novo at a quarter sampled, the key 20.2 % against 13.0 % (the de novo fixture
  is variance-dominated, as recorded). Hans validates on data. (Corrected 2026-11-03: the solve2D and solve3D
  "window" numbers here were not controls; their stages ran the map, the last progress entry.)
- Hans's beta-gal solve3D with the dynamic centring (2026-11-02): every particle on the disk at the workflow's
  first shift iteration (4598 of 4598, the table 26.7 s against 7.8 s after), then 37, 37, 11 and 5-10 movers an
  iteration; no change to the stage's course or the final 4.03 A. The coarse map's width reads 0.00 for nearly
  every particle (a resolution limit: the cell is 3.6 px against a 4.4 px window), so the re-centring is carried by
  the move at the bound; the Laplace width at the polished optimum is the informative alternative, written up in
  the proposal's section 2.6 for Hans's call. The neighbourhood driver's disk print was missing and is added.

### 2026-11-03: the review's section 15, M1-M4 (Claude, agent kit runs 1230-1292)

- M1, the angular Nyquist bin: the Euclidean score rows (`gen_raw_euclid_vals_impl`, `gen_euclid_crvec`) and the
  continuous Euclidean evaluator (`gen_raw_euclid_grad_at_angle`) summed the angular spectrum to `pftsz`, dropping
  the bin the map, the direct residual, the cc kernels, the bank and the batched rows carry; all sum to `pftsz+1`
  now, so one energy convention holds across the map, the rows, the polish's objectives and the continuous route.
  `test_angular_nyquist`: on random spectra the row is the direct residual to 3.7e-8 at every rotation, unshifted
  and shifted, and the continuous evaluator the row to 2.6e-8 at the grid angles; on the alternating half ring
  (all energy at the Nyquist bin) the four give 0 at rotation 1 and 4 at rotation 2 (the row gave 2 and 2 before).
  The bank tests' tolerances went from 5e-2 to 1e-4; one of them compared the cross bank loss with the shifted
  kernel's best over all rotations, which against an unrelated half is another rotation, and now compares at the
  bank's own rotation.
- M2, the record and the reread: the preparation records the mask it normalised a particle under
  (`builder%obs_mask`, `note_obs_mask`, in `prepimg4align` and the 2D cropped path), and every later normalisation
  in the process reads it (`obs_on_box`): both backends' rereads (`prep_rec_observation_ptcl`), the gridding
  full-size path and the 2D class sums from full-size images. Without a preparation in the process (a standalone
  reconstruction) the ordinary mask, as before. `test_rec_observation_mask` (an off-centred Gaussian on a structured
  background): the reread equals the record exactly; the ordinary mask would differ by 25 % of the image's range.
- M3, orders beyond the ring: the spectral index of the harmonic convolution wraps modulo the ring
  (`gen_shiftmap`); `test_orders_beyond_ring` maps orders to 64 on a ring of 32 against the direct objective to
  1.5e-7. Before, the slice wrote past the array: Debug stops on it, Release would corrupt memory silently.
- M4, the insertion's normaliser: `top_posterior_extras` takes its minimum over the committed and the eligible rows;
  the review's case (distances 0, 0.9, 0.901, the first excluded) keeps both poses at 0.524980 / 0.475020 and is
  unchanged when the excluded row moves.
- Builds and tests: Debug and the gate (run 1235); the in-plane state 129/129, the polar shift map 79/79, the polar
  correlation 80/80, the continuous in-plane 168/168, the 3D table 508/508 (runs 1236-1240); Release and the gate
  (run 1241).
- Fixtures after M1-M4 (true references or volume, the key): the 2D table 100 % (0.18 px; 5 px window 0.22 px);
  the 3D table 100 % (0.61 px; 5 px window 96.2 %); 3D `shc` 96.5 %, `neigh` 100 %, 2D `greedy` 78.0 %, 2D `snhc` at
  5 px 35.5 and 60.0 %, all as before within their spreads. solve2D from its random start (four runs each): the
  key 80.0, 39.8, 96.0, 77.2 %, the control 98.8, 97.8, 91.2, 89.0 %; on the build before today's changes the key
  100, 59.8, 75.5, 99.2 % and the control 19.5, 97.5, 97.0, 77.8 %. Its first, shift-free stage collapses classes at
  random (3 to 6 of the 10 requested populated, purity 58-99 %), so this fixture ranks neither the key against the
  control nor one build against the other.
- The initialiser's mode, measured (the review's statistical point): a build whose greedy cell is the density's mode
  (the log cell area left out of the pick, the draw unchanged) against the mass mode, interleaved pairs at a 5 px
  window: 2D `snhc` 57.8, 40.0, 40.0 % against 60.0, 52.8, 60.0 %; `snhc_smpl` 60.0 against 59.8 %; `greedy` 78.0
  against 78.0 %; 3D `shc` 97.2 against 97.0 %; the 3D table 96.2 against 96.2 %. The density mode does not help and
  the stalled SNHC runs are not the area term's; the mass mode stays (the variant was a scratch build, runs 1366-1422).
- Release on the final source and the gate (run 1423).

### 2026-11-03: the prior moves with the particle (ruling 22; Claude, agent kit runs 1135-1228)

- As built: the map's prior is on the cell's shift, the increment from the stored shift the particle is prepared
  at (`shiftmap_prior` and `shiftmap_draw` lose their centre argument and its per-rotation tables: the penalty is
  `prior_w rho^2`, the same at every rotation and direction); the insertion's loss is the data term with the prior
  on the pose's increment (`insert_pose_loss(dist, du, prior_w, temp)`, `map_prior_weight`). Every caller (the four
  fills, the classic map search, both insertions) drops the stored shift it passed only for the prior. The shift-map
  kernel test loses its stored-shift case (the draws centred on the pick); the insertion test's poses are
  increments of 20 and 22 px (the same weights).
- Builds and tests: Debug and the gate (run 1135), the polar shift map suite 71/71, the 3D table 504/504, the
  in-plane state suite 121/121 (runs 1136-1138), Release and the gate (run 1139).
- A/B against the build before this change (snapshot of run 1062's Release, `run@close`), refine3D `prob` from the
  true volume, the key, `trs=5`, starting from the true shifts (the case where the prior at the pick pulls back),
  from the true shifts with 2 px of noise, and from zero: 100 % within lp/4 in every pair, at SNR 0.1 and at SNR
  0.03 (the same geometry simulated anew), medians within 0.01 px. With the true volume as reference the data term
  outweighs the prior by tens of nats wherever it sits, so this fixture cannot show the prior's centre; the
  conditions where it matters (weak references, small particles, the first centring after 2D) are Hans's data.

### 2026-11-03: the close (ruling 21; Claude, agent kit runs 1057-1133)

- **The controls were not controls (a bug since the default of 2026-11-02, `bfef77378`).** The stage
  controllers switched the key off by deleting `shift_sigma` from a child command line, which since the parameter's
  default became 10 A gives the child the key. Every solve2D and solve3D run given `shift_sigma=0` ran the map in
  its shift stages (run 783, the 2026-11-02 solve2D "window" fixture, prints the disk set: 400 of 400 not yet
  centred); standalone refine2D and refine3D controls were right. Fixed by setting the key explicitly, zero where
  it is off: solve2D's stages and the seeded restart's probabilistic pass, solve3D's stages, refine3D's Cartesian
  polish (its command line and its job description), solve2D's terminal greedy pass (which carries the key, ruling
  21). The solve2D fixture at full sampling, a true control now: the control 19.5 % within lp/4 (two classes left,
  purity 40 %, no disk line in its log) against the key's 100 % (purity 100 %, 22 disk lines). Any control Hans ran
  through solve2D or solve3D between 2026-11-02 and this fix had the map on in its shift stages.
- **T3:** under the key the window search's increment prior is off (`l_shift_prior` false when `l_centring`): the
  tables' window polish, the classic polish and the continuous route are the data's alone. `test_key_derivation`
  pins the key's derivation per mode, the prior's switch and `shift_sigma=0`.
- **T2:** `insert_pose_loss` (with `map_prior_weight`) is the one definition of the insertion's pose loss; the
  dense fill calls it, and the neighbourhood's `select_extras` replaces each refined candidate's data term by it
  under the key. `test_insert_pose_loss` gives the review's 0.603483 / 0.396517 at two temperatures.
- **PCG under topk:** the PCG worker takes the weighted insertion's poses: `pcg_extra_entries` lays the extras out
  as operator entries after the committed ones (the committed weight on the committed entry, each extra's on its
  own, the extra's Euler triplet and absolute stored shift at the cropped box) and `accumulate_pcg_extras` adds a
  batch's extras from its particles' observations; the matcher reads the insertion table under either backend.
  `test_pcg_extra_entries` (the rec3D backend suite): extras at the committed pose with weights summing to one
  reproduce the plain accumulation in one batch and in two, a moved extra changes it.
- Builds and tests: Debug and the gate (run 1058), the in-plane state suite 121/121, the 3D table suite 504/504,
  the reconstruction suites 546/546 (runs 1059-1061), Release and the gate (run 1062).
- Fixtures (true volume or references, the key unless said): the sparse 2D table (`prob_snhc`, 5 px window) 100 %
  against its control's 79.5 %; refine3D gridding with `topk=3` 100 %; the neighbourhood table with `topk=3`
  97.8-98.2 %; at 2000 directions, where the cone holds 6-7 eligible poses per particle, the posteriors are sharp
  enough that the extras mostly fall under the insertion's weight floor (1.00-1.04 poses per particle), and PCG with
  `topk=3` matches `topk=1` (83.0 % against 82.5 %). PCG trails gridding on this fixture with or without today's
  changes: at 500 directions and `topk=1`, 22-47 % in four runs of the build before the close, 29.5-58.8 % in three
  after, gridding 100 %; a property of PCG on 400 particles in three iterations, not of the insertion.

### 2026-11-03: Phase 4's classic strategies built, the continuous route's domain fixed (Claude, agent kit runs 802-)

- Ruling 20 as built: `simple_shiftmap_srch` holds the classic searches' map context (the particle's box and
  temperatures, the window's and the disk's polish objects, built per run by `setup_shiftmap_srch` in both
  matchers) and the dense table's screen and polish for a list of candidates scored at the stored shift: the pool
  of the best by that score (`screen_pool_size`, at least twice the polish count), a coarse map each, the best
  cells polished with `minimize_in_cell`. Under the key the classic searches seed nothing (`inpl_srch_first`, in
  any state count), the previous candidate competes at its stored-shift score, `inpl_srch` polishes one,
  `inpl_srch_peaks` the peaks, the written class carries its own mapped shift (2D `assign_ori`), and the commit
  records the shift state (2D `assign_ori`, 3D `extract_peak_ori`). Under `shc` (3D) and in 2D a particle that
  needs the disk is searched greedily (every candidate scored), as the sparse table scores every class for one.
  `l_centring` covers the classic polar modes (2D `snhc`, `snhc_smpl`, `snhc_smpl_many`, `greedy`; 3D `shc`,
  `neigh`, `greedy`, `greedy_inpl`), not `tseries`, `inpl`, `inpl_smpl`, `greedy_smpl` or `cont`.
- The review's T1 fixed: the continuous route's joint domain contains the particle's box (`shift_box_radius`: the
  disk when `shift_needs_disk`, under the key) and the incoming shift, in 2D and 3D, for the committed pose and the
  extras. `test_continuous_seed_beyond_window` (in-plane state suite) polishes a pose found at 5.3 px under
  `trs=2` and checks it comes back unclipped and scored no worse than itself; with the old bound it fails at all
  five starts on the score (the shift kept, the clipped pose's score attached), as the review traced. The table
  fixtures are unchanged by the fix: 2D 100 % at 0.19 px, 3D 100 % at 0.54 px.
- The fixtures, true references or volume, 400 particles whose true shifts have a median of 8.6 px, the key
  (`shift_sigma=10`) against the window (`shift_sigma=0`), within lp/4 (2.31 px):

  | mode | window | key | no key |
  |---|---|---|---|
  | refine3D `shc`, 3 iterations | 5 px | 95.2 %, 96.0 % | 78.8 % |
  | refine3D `neigh`, 3 iterations | 5 px | 100.0 %, 100.0 % | 86.0 % |
  | refine2D `greedy`, 12 iterations | 5 px | 78.5 % (purity 98.5 %) | 43.0 % (67.5 %) |
  | refine2D `snhc`, 12 iterations | 24 px | 79.8 % (100 %) | 20.2 % (67.5 %) |
  | refine2D `snhc` | 10 px | 65.2 % | 39.5 % |
  | refine2D `snhc` | 5 px | 39.2, 39.5, 40.0 % | 57.0, 75.2, 41.5 % |
  | refine2D `snhc_smpl` | 5 px | 60.0, 42.2 % | 55.2 % |
  | refine2D `snhc_smpl` | 10 px | 40.0 % | 70.5 % |

  The 3D modes, 2D greedy and SNHC at a wide window gain clearly. The stochastic 2D modes at small windows are
  not resolved by this fixture: its classic 2D search first searches shifts at iteration 3, against class averages
  rebuilt for two iterations without shifts, so the first centring is poor (43.5 % at iteration 3, against
  84.8 % for the table's first centring from the true references), the classes stay about 80 % pure in every
  SNHC run, key or not, and the window's own repeats spread from 41 to 75 %. One pattern repeats: SNHC at a 5 px
  window, about one coarse map step on this 248 px box, ends at 40 % with the particles stopped (last move
  0.01 px) about 3 px from their view-mates. A first build that mapped only the class SHC committed to (a wrong
  class's window map moves the particle to that class's optimum) was worse; the pool over every scored class and
  the previous one replaced it. The stream's pool (`snhc_smpl`, appended particles greedy against settled class
  averages) is closer to the table's first centring than to this fixture; its measurement is Hans's data.
- The lattice inside the box polish retired: `minimize_in_cell` starts the L-BFGS at the map's cell
  (`minimize(xy_in=cell, prior_origin=.true.)`: the shift prior, when set, stays centred at the stored shift), keeps
  the better of the polish and the cell at its best rotation, and reports whether it beat the stored shift's best
  rotation (the callers' contract). A/B from one source, the lattice build snapshotted (`run@lattice`) against the
  cell start, two repeats interleaved, the table fixtures (true references or volume, 400 particles, the key, 8
  threads, load 9-10 on 24 cores); the repeats agree within 1-2 %:

  | fixture | refinement (thread-s), lattice / cell | fill wall (s), lattice / cell | objective evaluations, lattice / cell | within lp/4 |
  |---|---|---|---|---|
  | 3D, trs=24, 3 iterations | 102.1 / 49.5 | 13.6 / 7.0 | 418 M / 80 M | 100 / 100 % (median 0.60 / 0.61 px) |
  | 3D, trs=5 | 84.6 / 31.1 | 11.4 / 4.6 | 418 M / 80 M | 100 / 100 % |
  | 2D, trs=24, 12 iterations | 14.6 / 6.4 | 2.05 / 1.00 | 60 M / 9.9 M | 100 / 100 % (0.18 / 0.18 px) |
  | 2D, trs=5 | 10.8 / 2.6 | 1.45 / 0.44 | 60 M / 9.9 M | 100 / 100 % |

  The classic fixtures with the cell start: 3D `shc` 95.8 %, `neigh` 99.8 %, 2D `greedy` 78.5 % (the lattice:
  95.2-96.0, 100, 78.5). The saving is larger than the 2026-11-03 Phase 5 entry estimated from the beta-gal profile
  (a quarter of the refinement): on these fixtures the polishes were most of the refinement. `test_map_initialised_polish`
  holds with the cell start (same optimum and rotation as the lattice-initialised search). Beta-gal is not rerun.

### 2026-11-03: Phase 5 not pursued; the polish measured on Hans's beta-gal solve3D run (Claude's reading)

- The cost of the polish, read from `5_solve3D/BENCH/PROB_TAB_BENCH_PART1.txt` (1150 particles, 6 threads, box
  120, shells 2-36, the dynamic centring):

  | iteration | fill wall | iteration wall | refinement thread-s | maps per particle | polishes per particle |
  |---|---|---|---|---|---|
  | 1 (every particle on the disk) | 21.7 s | 34.4 s | 121.6 | 110 | 55 |
  | 2 (the window) | 6.0 s | 14.7 s | 27.7 | 38 | 19 |
  | 11-15 (the window) | 1.9 s | 8.6 s | 2.9 | 4 | 2 |

  The refinement is the maps plus the polishes; the profile cannot separate them because the pool is exactly
  twice the refined count in every block (`npool`'s floor). With the 2026-10-31 like-for-like record at box 120
  (a window map 0.24 ms, a classic minimisation 0.27 ms) the window iterations split about two thirds maps, one
  third polish; the disk iteration is the disk maps.
- What the polish is: 25 grid evaluations per search in every block, i.e. the object's 5 x 5 lattice runs over
  the box because `minimize_in_cell` calls `minimize` without a start (`xy_in` is what turns the lattice off),
  then the L-BFGS (7.8 gradients per search). The objective count agrees: about 34 all-rotation evaluations per
  polish, 25 of them the lattice. The L-BFGS proper is a quarter of the polish, under a tenth of the refinement,
  the refinement 70-80 % of the fill, the fill 40 % of an early iteration and 20 % of a late one: Phase 5's best
  case is about 3 % of an early iteration and 1 % of a late one.
- The band: both maps are built on the coarse band ending at `kscreen` (`setup_shiftmaps`), the polish is what
  brings the shift to the full band. A Fourier-Bessel polish on the map's intermediate works at the map's band;
  at the full band the intermediate doubles the shells and the orders kept, about four coarse maps per
  candidate, more than the polish. Hans: not pursued.
- The levers, in order of return per effort, each measurable on the in-plane tester's identity and the fixtures:
  the lattice inside the box polish (start the L-BFGS at the map's cell with the box as bounds, as ruling 14
  ruled; about three quarters of the polish, a quarter of the refinement; the risk to measure is a full-band
  secondary optimum inside the box, the lattice at half a step the present guard, a 3 x 3 lattice the fallback);
  the pool's floor of twice the refined count (the screen-rank record: almost nothing beyond the refined count
  from the later stages on); the map's own costs (measured FFTW plans, fewer directions on the inner rings).

### 2026-11-01: solve2D validated on data (Hans's runs; the approach adopted)

- With the drain reset in place and the key solve2D's default, Hans ran the consolidated centring stage on the
  export gate set, PfCRT and TRPV1: the classes gain resolution and visual quality on the first two, TRPV1 works
  well, and the 2D harvests more useful particles than before. Hans: "I consider the approach validated and I am
  ready to take the next steps in the dev plan." Phase 3's 2D question (does the map recover what the window
  cannot, at a cost the data sets bear) is answered yes; the per-class shift histograms against cryoSPARC's on
  the membrane-protein set, the mask ruling and re-extraction (Phase 3's exit) remain to be recorded.

### 2026-11-01: the 3D consolidation and Phase 4 in the tables (Claude, autonomous run, agent kit runs 604-)

Hans's bundle for an autonomous run: the 3D side of the consolidation, Phase 4, and the smaller items that fit.

- The 3D dense fill is the 2D centring stage's design (`process_particle_with_map`): rows at the stored shift, the
  screen pool by `screen_pool_size` (moved to the utils module, shared), the coarse map over the box, the polish in
  a box of one map step about the cell (`pftc_shsrch_grad%minimize_in_cell`, the one home of the box polish and its
  origin fallback, used by every table), the refined entries' loss the polished loss with the absolute prior at the
  physical temperature. The drain: `strategy3D_prob` writes `shincr` after the polish (the extras and the committed
  pose now take the continuous in-plane route under the key too, since the rows are polished losses like the
  window's), both refine3D masters reset it at the start of every run under the key, `prob_align_exec` prints the
  disk set, and `check_conv3D` makes the drain the centring stage's convergence (the pose rule, the annealing tail
  and the resolution exit deferred while the disk is needed; `test_shift_map_exit3D`). solve3D: stage 3 is the
  centring stage (`cfg%trs_window` the next stage's half-width). The full map (`SHIFTMAP_FULL`) is no longer built
  by the fills; `shift_draw` and the fill's `maps` count retired. Debug gate 15/15 (run 604). The random-view
  fixture (refine3D refine=prob from the true volume, 3 iterations, Debug tree, run 615): the disk searched by 400,
  303, 16 of 400 particles, the drain 76, 4.0, 0.75 %, 99.8 % of the shifts within lp/4 at 0.54 px mean error (the
  retired full map: 100 % at 0.56 px; the 5 px window: 99.2 % at 0.62 px).
- Phase 4 as built: under the key every table's polish starts from a coarse map over its box. A third map geometry,
  `SHIFTMAP_WINDOW` (the coarse band, radius the stage's window: `trs_window` in the centring stage, `trs`
  elsewhere), is set up beside the disk's screen map (`setup_shiftmaps(sigma, radius_disk, radius_win, ...)`). In
  the centring stage a drained particle takes the window map (no prior) through the same screen and polish as the
  disk; in every later stage every particle does, in the 2D dense and sparse (`prob_neigh`) fills, the 3D dense fill
  and the neighbourhood table's candidates (`refine_best_neighbors`), so the previous reference's seed, its
  shifted-particle bank and the 5 x 5 lattice at every rotation retire wherever the key runs (`shift_sigma=0` keeps
  them). The controllers pass the key to every stage with a shift search (solve2D stages 2-6, solve3D 3-8), the
  window only to the centring stage; `prep_strategy3D` admits `prob_neigh`. The classic (non-table) strategies keep
  their lattice: solve2D's stage 1 and solve3D's stages 1-2 search no shifts, so the key never reaches them in the
  workflows. The pin: `test_map_initialised_polish` (the in-plane state suite) compares the map-initialised polish
  with the lattice-initialised search on the white-reference fixture at five rotations.
- Validated (Release, runs 620-668): the gate 15/15, the in-plane state suite 98/98 and the 3D table suite 496/496
  at four threads; the 2D refine2D fixture from true references 99.8 % within lp/4 at 0.23 px mean (0.18 px median;
  before this run 95.8 % at 0.21 px median: the window map's polish in the later iterations recovers the coarse
  map's class misses); the 3D refine3D fixture from the true volume 100 % at 0.53 px (the 5 px window 86 %, the
  16 px window 100 % at 0.53 px); the solve3D fixture de novo (400 random-view particles, no 2D solution): the key
  78 % within lp/4 and 5.66 A against the window's 56.5 % and 6.58 A, the centring stage drained over six
  iterations (30, 33, 14, 4.5, 0.25, 0 % still needing the disk) and ended by the drain; beta-gal solve3D from the
  bench project's 2D solution (d2, 4463 particles): the key's centring stage drained 13, 2.0, 0.7, 0.2 % over four
  iterations (the disk searched by 4463, 588, 91, 32), both runs end at 4.03 A, the key 1002 s against 815 s
  (stage 3's four disk iterations). One defect found and fixed: the 2D controller gave the key to stage 1, whose
  SNHC strategy has no table and refused it (the solve2D fixture and the beta-gal 2D key run stopped at stage 1);
  the key now starts at the centring stage. The beta-gal 2D control (`shift_sigma=0`): 64 classes, best 4.29 A,
  26 % at 6.5 A or better, 50 % at 8 A, 99.8 s.
- After the fix (Release, runs 669-686): the solve2D fixture from the random start, the key 96.8 % within lp/4 at
  1.62 px mean against the window's 78.0 % at 2.96 px; beta-gal solve2D with the key (the default): no collapse
  (90 populated through the centring stage's five iterations, 57 at the end against the control's 64), best class
  4.03 A against 4.29, 24 % at 6.5 A or better against 26 %, 58 % at 8 A against 50 %, 143 s against 99.8 s. The
  centring stage drained 28, 23, 8.6, 3.1, 1.5 % over its five iterations (the disk searched by 4463, 1254, 1019,
  384, 138) and ended on its budget, one iteration short of the drain's two-below-2 % rule; its table cost 5.7 then
  4.0-4.3 s per part and iteration whatever the disk count, since the drained particles' window map costs about
  what the disk does. The later sparse stages 0.5-0.6 s against 0.4 s; the last, dense stage 3.6-4.1 s against
  1.3-1.4 s (the screen over every class, then the polish): the map-initialised window polish is the price of
  Phase 4 in 2D, +43 s on beta-gal, to weigh against its gain on data. The base tree (HEAD before this run) on the
  solve3D fixture: the window 1.5 % within lp/4 (2142 s), the key 13.2 % (1191 s), against this tree's 56.5 % and
  78 %: the window path's code is the same on both trees, so the de novo fixture's run-to-run variance dominates
  (as recorded for Phase 2d) and the solve3D fixture does not attribute; the refine3D fixture from the true volume
  and the beta-gal 3D pair are the 3D evidence.

### 2026-10-31: Phase 2c measured, Phase 2d's first stage measured, the re-review applied (Claude, runs 383-459)

- Phase 2c on the bench, the field through the particle record (slot 54, `shmassout`) and the mass at the physical
  temperature: the disk-needed fraction (posterior mass beyond the 4.4 px window above 5 %) drains 80, 53, 35, 22,
  15 % over the centring stage's five iterations (the fraction above 50 %: 75, 42, 25, 15, 8 %), so the stage ran
  its budget under either tolerance; the screen breadth went 90, 90, 76, 58 in the centring stage and 18 (twice the
  refined count) in every later stage, whose fills fell from 1.2 to 0.4 s and the dense stage's from 4 to 1.7 s; the
  class quality was unchanged by the length of the stage (a two-iteration run on a zero field gave the same). Which
  statistic ends the stage (the mass above a tolerance, the posterior's mode beyond the window, or either relative
  to the first iteration) is Hans's ruling; all three run the stage to its budget on beta-gal.
- Phase 2d, the first stage on the shift-marginalised evidence (every class mapped over the disk, no shift assigned):
  stage-1 purity on the five-view fixture over four random starts 80, 60, 40, 60 % against 31, 20, 25, 33 % without
  (`nstages=1`); two full fixture runs each inconclusive (the later stages' randomness dominates there); on the bench
  65 % of the particles in classes at 8 A or better, best class 4.35 A, the best set so far, at 88 s for stage 1
  (14 s per iteration: 0.7 ms per map at 15 A against the 0.4 estimated), 238 s for the run. The cost is the open
  item: the evidence over the screen's band, or over a pool after the first iteration, are the candidates.
- The automatic mask: Hans's full-bench pair with the estimator fixed still put the automatic run 1.4 A behind at
  the top; the cause is the ladder and stage 1's mask from the box-wide default; an automatic run now measures its
  first mask from the mean of up to 2000 unaligned particles (built, measured in the next chain).
- The review's section 12 (six findings) applied the same day: the dispositions in the review table above. The frame
  finding (R1) was mine from the morning: `extra_pose` had assumed the map's frame for the window path's native
  entries; the table now has one contract, native increments, converted at the producer.

### 2026-10-31: the automatic mask on the full bench (Hans's runs 4_solve2D, 5_solve2D; Claude's reading)

- With the estimator fixed (183.5 A against 180) the automatic run is still worse: best class 5.63 A against 4.24,
  six classes at 6.7 A or better against twelve, 169 s against 186. From stage 3 on the masks agree within 2 %;
  what differs is before: without a mask on the command line the ladder is derived from the box-wide default
  (stage 1 at 20 A, then 14, 11, 9.5, against 15, 11.5, 9.8, 8.9) and stage 1 runs under the box's mask. The gap
  is visible at the end of stage 2 (best 8.0 A against 7.2). Built: an automatic run measures its first mask from
  the mean of up to 2000 unaligned particles (the molecule blurred by the pick error, an overestimate by about the
  picking spread), which sets stage 1's mask and the ladder; the centring stage's class averages refine it after.

### Phase 2d: the shift-invariant first stage (Hans, 2026-10-31; the "no shift" part reverted 2026-11-01)

Revised 2026-11-01 (Hans: "it collapses into too few classes when it hits stage 2"): stage 1 under the key runs the
dense table over the disk like the centring stage, the class from the evidence and the drawn cell assigned; the
`shift_assign=no` variant below left the particles at their picks, so stage 1's averages were blurred by the pick
error and the centring stage, mapping the disk against them, collapsed the populated classes from 90 to 50 in its
first iteration (38 at the end; 57 with the shift assigned, 65 for the control; the per-iteration counts from the
runs' `clusters2D_iterNNN.star`). The marginalisation RELION applies is in the E-step; its M-step still places the
particle at its posterior shifts, which the drawn cell is here. The mode and its key are gone. The variant with the
drawn cell assigned still collapsed the classes at the centring stage's entry (Hans's run, 2026-11-01), so stage 1
under the key is the plain SNHC window stage again, as in Phase 2b (79 populated at the centring stage's entry, 57 at
the end, against 65 for the control; the pooled screen of 2c gives the same endpoint). Phase 2d is closed without a
shift-invariant stage 1; the measured purity gain on the fixture did not carry to the bench's class populations.
Original text:

Stage 1 is shift-blind: at zero shift against random references a pick off by most of a period at 15 A matches
nothing, and the fixture's class collapses happen there. The shift-invariant objective is the marginal evidence
over the disk, which the map returns in one pass: under the key solve2D's first stage runs the dense table with
every class mapped, draws the class from the evidence and assigns no shift (`shift_assign=no`, the particle stays
at its pick), as RELION and cryoSPARC marginalise shifts in their expectation step. Cost about 0.4 ms per map at
15-20 A, 7 s per iteration on the bench. Measured on the five-view fixture (stage-1 purity over random starts,
`nstages=1`) and on the bench. 3D: the same evidence exists in the dense 3D fill, but the cost scales with the
references scored (nspace 500-1000 in solve3D's early stages) and the first stage is noise-driven by design; the
candidate is solve3D's stage 2 with the low-resolution map over a pool, after the neighbourhood table has the map.

### Phase 3: recovery on data

Synthetic: particles simulated with picking errors drawn up to 30 Å (uniform direction, magnitudes 0-30 Å),
solve2D and solve3D with the map on, the recovered shifts against the truth per stage (the histogram, the
fraction within lp/4, the failures). Real: one of the benchmarked sets, SIMPLE's shift distribution against
cryoSPARC's on the same particles. Report in Progress: the answer first (does the map recover the shifts
cryoSPARC finds, at what cost per iteration), then the tables.

Exit: the report; the ruling on masks and re-extraction recorded.

### Phase 4: the map initialises the classic search (ruling 14; the tables done 2026-11-01, Progress below)

An inventory first: where each mode seeds its shift search today (`estimate_shift_seed`, `seed_start`, the seed
bank, `coarse_search_opt_angle`, `grid_fit`, the shift prior's centre), for refine2D (`snhc_smpl`, `prob`,
`prob_neigh`), refine3D (`shc`, `neigh`, `greedy`, `prob`, `prob_neigh` in its modes, `cont`), solve2D's and
solve3D's stages, single and multi-state. Then `pftc_shsrch_grad%minimize` takes the map's refined best cell as
its start and one grid step as its bounds wherever a map exists for the pair, the lattice seed and the seed bank
retire where the map replaces them, and the in-plane tester's identity pins the optimum against the old path.

#### The inventory of today's shift initialisations (2026-10-31, read from the code)

The optimizer `pftc_shsrch_grad%minimize`: with `coarse_init`, `coarse_search_opt_angle` scores a 5 by 5 lattice
over the bounds (plus or minus `trs`) at every rotation and starts from the best (point, rotation); `grid_fit`
adds a quadratic step from the lattice and `grid_ptcl` scores the lattice from a particle-frame bank (both only
in the multi-state 3D table, `l_grid = nstates > 1`); with `xy_in` the lattice is skipped and the start is the
seed. Then the L-BFGS-B alternations with the in-plane callback (`new_alternating`), or a fixed rotation
(`new_fixed`: the 2D time series, the per-peak polish `grad_shsrch_obj2`). The shift prior, when on, is centred
at the stored shift.

| owner | seed of the shift search | start of each reference's search |
|---|---|---|
| 2D strategies (`snhc`, `snhc_smpl`, `greedy`, `greedy_smpl`, `inpl_smpl`; shifts from the third iteration, `do_inplsrch`) | `inpl_srch_first`: the lattice and L-BFGS at the previous class, converted to the particle frame (`xy_first`) | `xy_in = xy_first` at every class; the lattice when there is no seed |
| 2D `tseries` | none | the lattice at the fixed rotation |
| 2D tables (`prob`, `prob_neigh`) | the previous class's lattice and L-BFGS; the seed bank memoises the particle at it; every class scored at the seed | `seed_start`: the seed in the reference's frame, clamped to `trs`, as `xy_in`; under the key the map about the stored shift replaces all of it |
| 3D strategies (`shc`, `greedy`, `greedy_sub`, `greedy_inpl`, `neigh`) | `inpl_srch_first` at the previous reference, single state only (`l_sh_first = doshift .and. nstates <= 1`) | the references scored at `xy_first`, the final `inpl_srch` from it; multi-state: the lattice at the best reference, no seed |
| 3D `cont` | the Cartesian pose optimizer's own 5-parameter solve (`shift_first`) | not polar |
| 3D dense table (`prob`) | `estimate_shift_seed`: the previous (state, projection)'s lattice and L-BFGS, single state only; multi-state: seed 0 with `grid_fit` and `grid_ptcl` instead; every reference scored at the seed (bank) | `refine_best_refs` from the seed as `xy_in`; under the key the map |
| 3D neighbourhood table (`prob_neigh`) | modes `shc`/`snhc` (solve3D's early stages): the seed in `shc` mode only; modes `geom`/`state` (solve3D's late stages, multi-state): the seed single state only | the candidates' L-BFGS from the seed; no map |
| solve2D | stage 1 `trs = 0`; stages 2-6 the window (the key: stage 2 the disk) | as the stage's mode |
| solve3D | stages 1-2 `trs = 0`; 3-5 `prob`, 6-8 `prob_neigh` in `state` mode, multi-state `state` from stage 6 (`shc` before) | as the stage's mode |

Three facts for Phase 4: multi-state runs have no previous-reference seed anywhere (the 3D search object, the dense
table, the subspace modes), only the lattice; the neighbourhood table's `snhc`, `geom` and `state` modes seed
nothing; and every lattice is a 25-point grid at all rotations, a cost comparable to one map at the window's radius,
without the joint rotation. The map-initialised polish removes the asymmetry: one map per (reference, particle)
gives the joint start for single and multi-state alike, the seed bank and the lattice retire where it runs, and
the neighbourhood table gets the map it lacks.

The classic strategies (ruling 20; built 2026-11-03, Progress below): under the key the previous-reference seed
goes, in every state count, and the strategies take `process_particle_with_map`'s structure
(`simple_shiftmap_srch`). The candidate loop stays as it is, at the stored shift (SHC's first improvement, greedy's
sweep, the neighbourhood); the screen covers every scored candidate and the previous one, the best cells polished
with `minimize_in_cell`, one for the committed pose (`inpl_srch`), the peaks for `inpl_srch_peaks`, then
`record_shift_state`. `l_centring` extends to these modes. Open: the 2D gate `do_inplsrch` (shifts from the third
iteration) keeps a freshly injected particle off the disk for two iterations. The stream's pool runs `snhc_smpl`
(p06; the pool admits `snhc` and `snhc_smpl` only) and is centred this way now.

### Phase 5: the polish in the Fourier-Bessel domain (ruling 14; not pursued, Hans, 2026-11-03)

The design was: where the map is built, keep its intermediate, the per-rotation harmonic correlations
T_n(k, irot), and polish in (rho, alpha) on them, each evaluation a sum over shells and orders instead of
shells by ring samples, the gradient analytic. Not pursued, on the beta-gal 3D profile (the 2026-11-03 entry):
the L-BFGS it would replace is a quarter of the polish and under a tenth of the refinement, so its best case is
a few per cent of an iteration; and both maps are built on the coarse band ending at `kscreen`, so a polish on
their harmonics returns a half-band shift, while the full-band intermediate costs about four coarse maps per
candidate, more than the polish it replaces. The cheaper levers are the lattice inside `minimize_in_cell`
(ruling 14 as ruled), the pool and the map's own costs.

### File disposition

| file | why | phase |
|---|---|---|
| `doc/implementation_notes/completed/polar_shift_search.md` | this note: progress, rulings as they are applied (planned/ until the close, 2026-11-03) | each |
| `src/main/pftc/simple_polar_shiftmap.f90` | the mechanics class; phase 1 adds the table accessors of the all-rotation kernel | 0 1 |
| `src/main/pftc/simple_polar_shiftmap_tester.f90` | its unit tests; phase 1 adds the kernel's identity and timing | 0 1 2 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | sub-suite registration (done); phase 2 if a strategy sub-suite is added | 0 2 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | suite help list (done) | 0 2 |
| `doc/policies/test_environment_policy.md` | the area table (done) | 0 2 |
| `src/main/pftc/simple_polarft_calc.f90` | the kernel's interface, the shift map and its per-thread work as components | 1 2 |
| `src/main/pftc/simple_polarft_shiftmap.f90` | the all-rotation kernel (new submodule) | 1 2 |
| `src/main/pftc/simple_polarft_core.f90` | `kill` releases the shift map | 1 |
| `src/main/pftc/simple_polarft_memo.f90` | only if the memo must keep what the kernel needs | 1 |
| `src/main/pftc/simple_pftc_shsrch_grad.f90` | `minimize_in_cell`: the polish in a box of one map step about a map's cell, with the origin fallback (every table's) | 4 (done 2026-11-01) |
| `src/main/strategies/search/simple_strategy3D_inplane_tester.f90` | `test_map_initialised_polish`: the identity of the map-initialised and lattice-initialised optima | 4 (done 2026-11-01) |
| `src/main/params/simple_parameters.f90` | the radius / switch | 2 |
| `src/main/strategies/search/simple_strategy2D_srch.f90` | plumbing; the classic map path, the continuous route's domain | 2, 20 |
| `src/main/strategies/search/simple_strategy3D_srch.f90` | plumbing; the classic map path, the continuous route's domain | 2, 20 |
| `src/main/strategies/search/simple_shiftmap_srch.f90` | the classic searches' map context: the box, the screen and the polish | 20 |
| `src/main/strategies/search/simple_strategy2D_greedy.f90`, `simple_strategy2D_snhc.f90`, `simple_strategy3D_utils.f90` | every scored class kept under the key; the 3D commit's shift state | 20 |
| `src/main/strategies/search/simple_strategy2D_matcher.f90`, `simple_strategy3D_matcher.f90` | the maps built per run; a disk particle searched greedily | 20 |
| `src/main/solve/simple_solve2D_controller.f90`, `simple_solve3D_controller.f90`, `src/main/strategies/parallelization/simple_refine3D_strategy.f90`, `src/main/commanders/simple/simple_commanders_solve2D.f90` | the key set explicitly where it is off (zero) and on solve2D's terminal pass | 21 |
| `src/main/strategies/parallelization/simple_rec3D_pcg_strategy.f90`, `simple_rec3D_strategy_tester.f90` | the PCG worker's weighted insertion (`pcg_extra_entries`, `accumulate_pcg_extras`) and its test | 21 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90`, `simple_eul_prob_tab.f90`, `simple_eul_prob_tab_utils.f90` | one insertion loss (`insert_pose_loss`) in both 3D tables (T2) | 21 |
| `src/main/params/simple_parameters_phases.f90` | the increment prior off under the key (T3) | 21 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab2D.f90` | the joint draw in the fill of the centring stage | 2 |
| `src/main/pftc/simple_polarft_corr.f90` | `gen_posterior_val` over the map's grid | 2 |
| `src/main/solve/simple_solve2D_controller.f90` | stage 2 as the centring stage (`CENTRING_STAGE`): the table, the disk, the box's mask, minits for the shift exit, the key per stage; the mask estimate at its exit (`solve2D_mskdiam_estimate`) | 2b (done 2026-10-31) |
| `src/main/solve/simple_solve3D_controller.f90` | stage 3 the centring stage (`trs_window`), the key to every stage with a shift search | 2b 4 (done 2026-11-01) |
| `src/main/strategies/search/simple_strategy2D_alloc.f90`, `simple_strategy3D_alloc.f90` | the key refused outside the tables (`prob`, `prob_neigh`) | 2b 4 |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_neigh.f90` | the window map initialises the candidates' polish, no seed under the key | 4 (done 2026-11-01) |
| `src/main/ui/simple/simple_ui_solve3D.f90` | the key on solve3D | 2b |
| `doc/policies/3D/refine3D_policy.md`, `solve3D_policy.md` | the map's paragraph, the key's stage scope (done 2026-10-31) | 2b 3 |
| `src/main/commanders/simple/simple_commanders_solve2D.f90` | the centring stage's trs and mask, mskdiam=0 as the automatic mask, the estimate's adoption | 2b (done 2026-10-31) |
| `src/main/simple_convergence.f90`, `_tester.f90` | `check_conv2D`'s shift exit (ruling 8) and its test | 2b (done 2026-10-31) |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab.f90` | the 3D table consolidated (`process_particle_with_map`: the screen, the box polish, the drain), the window map in the later stages | 2a 2b 4 (done 2026-11-01) |
| `src/main/strategies/search/probabilistic/simple_eul_prob_tab_utils.f90` | `screen_pool_size`, `SCREEN_QUANTILE` shared by the tables; the dynamic centring's rule and state (`shift_needs_disk`, `disk_radius_px`, `next_shift_width`, `record_shift_state`), `sdev` in the records | 4 19 |
| `src/main/strategies/search/simple_matcher_2Dprep.f90`, `simple_matcher_ptcl_batch.f90`, `src/main/simple_builder.f90` | the disk particle's box mask (`ptcl_on_disk`, `lmsk_box`, `msk_box_crop`); the particle cache removed | 19 |
| `src/main/stream/pool2D/simple_stream_pool2D.f90` | the sentinel on the appended particles | 19 |
| `src/main/flex/fit/simple_flex_pca_plane_cache.f90` | its own location helpers, gated by `cache_dir` | 19 |
| `src/main/pftc/simple_polarft_calc.f90`, `simple_polarft_shiftmap.f90` | `SHIFTMAP_WINDOW`, `setup_shiftmaps` with the window's radius; the full map no longer built by the fills | 4 |
| `src/main/ui/simple/simple_ui_refine2D.f90`, `simple_ui_refine3D.f90` | the key on refine2D, refine3D (2a) and solve2D (2b) | 2a 2b |
| `src/main/image/simple_image_msk.f90` | `automask2D_set` writes its mask stacks only when asked (write2disk), as its flag intended; the estimate calls it silently | 2b |
| `doc/policies/2D/solve2D_policy.md` | the centring stage's contract (done 2026-10-31) | 2b 3 |

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

### 2026-10-31: the kernel optimised (Claude, agent kit runs 9-33; baseline `tree phase1` = 9c2fcdaf9)

- What changed. The direction synthesis and the transform over the ring are one batched 2D transform per
  magnitude of the shell sums G_n[m, rho] over (order, ring position), the orders folded modulo nalpha onto the
  direction axis (exact on the sampled directions for any count); since only the real part is used and the real
  part of a forward transform is the forward transform of the Hermitian half X_h = (X + conj(X[-n,-m]))/2, it is
  one complex-to-real transform of half the size (fed conj(X_h): FFTW's c2r is the backward transform). The
  shell products come from two contiguous slices instead of a modulo index; the shell sums are single precision
  like the spectra (the identity is unchanged: the spectra set the floor); the transform's input keeps its zero
  rows across pairs; the per-rotation normalisation is applied once per rotation. `shiftmap_ndirections(nmin)`
  gives the smallest 7-smooth count >= nmin: a prime direction count (47 at the disk) made the transform slower
  than the sums it replaced. Review-driven: the kernel is the data term only (no prior; the sampler assembles the
  untempered prior and the measure), `shiftmap_cell_area` gives the polar cell measure, `set_ptcl_ctfmat` lets a
  test set a per-particle angular CTF, and the grid layout is vals(ialpha, irot, irho).
- Identity (Release, box 64, pftsz 48, shells 2-11, two distinct reference halves, particle 1 odd with an
  angularly varying CTF, euclid with two sigma2 spectra): kernel against the mechanics 1.6e-8 (cc) and 1.3e-8
  (euclid) of the Cauchy-Schwarz scale; the grid against `gen_corr_for_rot_8` 5.2e-8 (cc) and 4.8e-8 (euclid,
  as exp(-L)); in the sampler's units the euclid error is 4.3e-5 nats, T |dL| at the particle's temperature; the
  zero-magnitude row 1.5e-8 and 3.7e-8; a 7-direction grid with 30 orders folded 2.8e-8; four threads bit for
  bit; cell areas sum to the disk. 53 checks, 0.15 s; both gates 15/15 (Debug run 30, Release run 32).
- Timing (Release, the fixture of the 2026-10-30 table, FFT-friendly direction counts, ms per pair):

  | lp (A) | coarse 5x5 | map at the window: before -> after (grid) | map at the disk: before -> after (grid) |
  |---|---|---|---|
  | 8.0 | 0.170 | 0.096 -> 0.051 (4 x 14) | 1.119 -> 0.487 (15 x 48) |
  | 2.7 | 0.474 | 0.923 -> 0.499 (10 x 35) | not a use case |

  Per grid point with every rotation included: 0.7 us at the disk and 1.4 us at the lp 2.7 window against 6.8
  and 19 us for the lattice, ten times cheaper, the low end of the Cost section's estimate. What is left: the
  c2r transform (about half), the shell sums, the symmetrisation pass. Measured plans (FFTW_MEASURE, what a
  production run without a fixed seed gets) would take more off the transform; fewer directions on the inner
  rings and the sampler's choice of counts are the remaining levers, both the sampler's business.
- The first c2r attempt returned the map with both indices negated; the zero-magnitude row caught it (a
  sign-convention bug is a reflection of the rotation axis). A folded-direction test that normalised by a value
  it was checking was replaced by the direct comparison against the calculator.

### 2026-10-31: Phase 2a, the switch and the sampler in the 2D tables (Claude, agent kit runs 34-119)

- The key `shift_sigma` (A, `simple_parameters.f90`, parsed, derived into `l_shift_map` after `l_doshift`, which it
  requires; refine2D's UI). The sampler on the calculator: `shiftmap_prior` (the Gaussian on the absolute
  displacement |u_prev + R(irot)^T x| times the cell's area over every cell, its log normaliser and the origin's
  log mass) and `shiftmap_draw` (the log masses -temp L - prior + log area, their maximum, the evidence as the
  log-sum-exp over the normaliser, one cumulative draw or the greedy cell, the data term at the cell). Tests: at
  temperature zero the draws follow the prior under the polar measure (mean radius 12.19 against the discrete
  prior's 12.28 and the continuous 12.30 of the review), the evidence is 0, the origin's mass is the central
  disk's share, a stored shift centres the draws on the pick; at the particle's temperature the evidence is the
  direct sum over the grid (4e-6) and the greedy draw its maximum. `gen_posterior_evid` gives the in-plane
  marginal of a class at the stored shift as an effective loss.
- The 2D tables, dense and sparse fills alike, with the key on: the seed is the stored shift the sections already
  carry (no L-BFGS seed); every class is scored at it by the posterior in-plane draw and its row holds the
  evidence, the in-plane marginal times the origin's prior mass, as an effective loss at prob_beta x T; the
  refined candidates get their map, their row the disk's evidence -log Z / (prob_beta T) and the drawn cell's
  shift and rotation. The master's class draw is unchanged and is therefore from the evidence (review, section
  3), every class on one scale. `strategy2D_prob` applies the drawn cell as the particle-frame increment
  (shift_ref2ptcl at the drawn rotation), scores the pose with the objective at the pose, and does not polish
  it (ruling, 2026-10-31). With the key off nothing changes: the gate's table, in-plane and strategy tests pass
  on both builds (Debug run 36, Release run 37). The 3D table and the controllers are not touched.
- Recovery on simulated particles (the 6VXX map of the simulate_particles test at 1.3 A/px, box 248, 400
  particles over five fixed views, in-plane angles random, Gaussian shifts of sigma 7 px clipped at 16 px: mean
  8.8 px, 11.5 A, largest 20.8 px; SNR 0.1, no CTF; refine2D refine=prob objfun=euclid ncls=5 lp=12 center=no
  minits=maxits=12 from the five true reprojections as references; recovered against truth in native pixels,
  the basin lp/4 = 2.3 px):

  | variant | within the basin | mean, median, 90th pct error (px) | class purity | iterations to settle |
  |---|---|---|---|---|
  | map, sigma 10 A, trs 5 | 88.8 % | 2.10, 0.73, 11.6 | 89 % | 2 (8.4 px of motion in the first) |
  | window trs 5 | 79.0 % | 2.47, 0.48, 10.8 | 79 % | 4 (5.1, 2.8, 0.9 px) |
  | window trs 16 | 82.5 % | 1.82, 0.25, 7.4 | 82 % | 3 |
  | no recovery (zero shifts) | 5 % | 8.8, 8.4, 15.0 | | |

  Every particle outside the basin is a misclassified one, in all three variants: the shift search never fails
  on a correctly classified particle, and the map's joint search classifies more of them. The map's median error
  is the grid's quantisation (the drawn cell is not polished); the windows' is the L-BFGS's. Misclassification
  grows with the true shift (map: 100, 87, 86, 79 % within the basin at |s| 0-5, 5-10, 10-15, 15-25 px): the
  class ranking at the stored shift excludes the right class for the particles that are far off, the review's
  section 5, which the translation-aware screen (the next 2a step) is for. The runs took 14.1 s (map) against
  15.9 and 15.7 s (windows), twelve iterations on 8 threads. The in-plane angles of the correctly classified
  are recovered to a median 0.35 degrees.
- What the fixture taught: random references with 400 random views over the sphere do not classify in twelve
  iterations (class purity 50 %, every variant at the null), and standalone refine2D calls its convergence
  after two iterations unless `minits` holds it; the recovery statistic needs per-class offsets removed when
  `center=no` (they were below 1 px here). The kit's runs needed the runtime environment CTest gives the tests
  (SIMPLE_QSYS, SIMPLE_PATH, the tree's bin on PATH, now exported by the kit's project file) and `mkdir=no` so
  the programs update the given project in place. The chain is `~/agent_work/SIMPLE/scripts/e2e_shiftmap.sh`,
  the statistics `compare_shifts.py` beside it.
- Left in 2a after this entry: the screen (next entry), the greedy test key, the 3D table, the fill's cost
  record, the in-plane tester's identity of the deterministic coarse path, and the frame/crop validation beyond
  this run (which went through SIMPLE_SCALE's downscaling with the shifts in native pixels and recovered them).

### 2026-10-31: the translation-aware screen (Claude, agent kit runs 121-134)

- The calculator holds two shift maps (`SHIFTMAP_FULL`, `SHIFTMAP_SCREEN`; every map routine takes the index):
  the screen's is built over the lower half of the band (`kscreen = max(kfromto(1)+1, kfromto(2)/2)`) with its
  own grid, magnitudes a quarter of that band's wavelength apart out to the same disk and a direction count
  resolving its ring, so its orders, magnitudes and directions are all about half the full map's; its losses
  and temperature are the band's (`wsqsum_ptcl_band`). Identity: the screen's map is the mechanics over its
  band to 3e-9 of scale (65 checks, both gates).
- In both fills, with the key on, the screen runs against every class scored (greedy, the evidence only) and
  its evidence over the disk ranks the classes the full map then refines; the rows keep the full-band
  quantities (the origin evidence of every class, the disk evidence of the refined ones), so the class draw's
  scale is unchanged. The candidate choice is therefore translation-aware (review, section 5).
- Recovery, the same fixture and references as the previous entry (five views, 400 particles, shifts of mean
  8.8 px, lp 12 A, twelve iterations):

  | variant | within the basin | mean, median, 90th pct error (px) | class populations | settled after |
  |---|---|---|---|---|
  | map + screen, sigma 10 A | 100.0 % | 0.75, 0.74, 1.18 | 80 80 80 80 80 | 1 iteration (8.7 px, then 0.5) |
  | window trs 5 | 79.0 % | 2.47, 0.47, 11.0 | 159 59 49 51 82 | 3-4 |
  | window trs 16 | 82.8 % | 1.81, 0.25, 7.4 | 149 72 50 54 75 | 3 |

  Every particle classified correctly and every shift recovered to the grid's quantisation (0.75 px at a
  2.3 px step: the drawn cell is not polished), after one iteration, in 15.5 s against 16.0 and 15.8 s for the
  windows. Misclassification was the whole of the previous entry's tail; the screen removes it on this fixture.
  The windows are unchanged, as the gate's identity requires.
- Left in 2a after this entry: the greedy test key, the fill's cost record and the 3D table (next entry), the
  in-plane tester's identity of the deterministic coarse path, and Phase 3's measurements at lower SNR and from
  random references (the centring stage's real start).

### 2026-10-31: the 3D table, the greedy key and the cost record (Claude, agent kit runs 136-228)

- `shift_draw=yes|no` (refine2D and refine3D UI, developer level): `no` takes the map's cell of largest mass, the
  comparison arm of Phase 3. The fill's cost record counts the maps and the screens (`maps`, `screens` of
  `prob_fill_profile`, printed by the bench writer). The calculator builds both maps from the prior's width
  (`setup_shiftmaps`), which both tables call.
- The 3D table (`eul_prob_tab`), with the key on: as in 2D, the seed is the stored shift, every reference is
  scored at it with the origin evidence in its row, the screen ranks a pool of the best 4 projs_ns of them by
  the coarse map's evidence, and the projs_ns best get the full map, their row its evidence and the drawn
  cell. `strategy3D_prob` applies the drawn cell without polish, scored by the objective at the pose, converted
  to the native frame as the polish converts it: the first run without that conversion moved every particle by
  the right amount in the wrong direction (17 % within the basin), which is how the frame contract was pinned.
- Recovery in 3D (the random-view stack of 400 particles, shifts uniform to 12 px per axis, mean 8.8 px; refine3D
  refine=prob from the true volume, nspace 500, lp 12 A, two iterations so the window cannot walk the distance):

  | variant | within the basin | mean, median, 90th pct error (px) | motion per iteration (px) | time |
  |---|---|---|---|---|
  | map + screen, sigma 10 A | 100.0 % | 0.88, 0.91, 1.43 | 9.0, 0.4 | 23.4 s |
  | window trs 5 | 91.8 % | 1.00, 0.70, 1.85 | 5.8, 3.1 | 14.9 s |
  | window trs 16 | 99.5 % | 0.60, 0.60, 0.84 | 9.2, 0.3 | 15.3 s |

  The map's error is the grid's quantisation; the windows' is the L-BFGS's. With six iterations every variant
  reaches about 99 % here (the window walks), so the iteration count is the discriminating quantity on this
  fixture, as in 2D. The map run costs 8 s more: the screen over the pool of 220 references and the 55 full
  maps per particle (the counts are in the bench record), the review's cost question made concrete, to be
  reduced with the pool and the screen's band.
- Both gates pass on both builds (Debug run 191, Release runs 192 and 217; the Debug gate of the final tree
  below). The UI registry validation caught two description mistakes: a binary input declares its two choices
  explicitly, and the parser does not read them from the help text.
- Left in 2a: the in-plane tester's identity of the deterministic coarse path; the pool and band of the 3D
  screen against its cost; Phase 3 at lower SNR and from random references.

### 2026-10-31: the random start (Claude, agent kit runs 231-241)

- The same five-view 2D fixture and twelve iterations, but from random references (refine2D's own start: five
  random images), the centring stage's real situation:

  | variant | class purity | classes (view counts) | shift error, mean (px) | within the basin |
  |---|---|---|---|---|
  | map + screen, sigma 10 A | 60 % | views 1+3 merged, views 0+4 merged, view 2 split in three | 5.0 | 1.5 % |
  | window trs 5 | 47 % | mixtures | 8.3 | 0 % |
  | window trs 16 | 52 % | mixtures | 7.4 | 0.2 % |
  | no recovery | | | 8.8 | 5 % |

  The map produces clean merged classes (two views each, 160 particles) and halves the error; the windows
  produce mixtures and recover nothing. With the key the shifts lock after one iteration (9.8 px, then 0.6):
  once the classes have formed, wrong or not, every member's shift is consistent with its class average, and
  nothing in a standalone refine2D from random images separates merged views. This is solve2D's stage-1
  problem (the sticky sampling, the Gaussian-regularised references, the class count above the view count and
  the recentring are its machinery), so the next measurement is through solve2D itself with the key handed to
  its stages: Phase 2b, the controllers.
- Phase 2a is complete in the working tree (17 files, both gates on both builds, Debug run 230).
- The extras' frame (open decision 3) fixed and tested; the 3D map run with topk=3 unchanged at 100 % (run 248).

### 2026-10-31: Phase 2b, solve2D's centring stage (Claude, agent kit runs 258-312)

- Built, under the key: solve2D's stage 2 is the centring stage (`CENTRING_STAGE`, `solve2D_centring_stage`): it
  runs the table of the stages that follow, its `trs` is 3 sigma in native pixels (at least the window) and its mask
  the box's, its `minits` is its first iteration so that `check_conv2D`'s new shift exit (ruling 8: the median shift
  change of the sampled particles below lp/4 for two iterations, `test_shift_map_exit`) can end it; at its exit the
  mask diameter is estimated from its class averages (`solve2D_mskdiam_estimate`: automask2D's diameters of the
  populated classes, their population-weighted 90th percentile, through `automask2D_mskdiam`), logged, and adopted by
  the later stages under `mskdiam=0` (a positive mskdiam is honoured from stage 3 on). The later stages map their
  window: `setup_shiftmaps` takes the disk radius as an argument and both tables pass min(3 sigma, trs). Stage 1 never
  sees the key; the stream's checkpointing refuses it; `cls_init=prev` skips the centring stage. solve3D: the key
  goes to its stages with a shift search and stage 3 (`PROB_REFINE_STAGE`) opens the disk the same way (measured
  below when done). `automask2D_set` now writes its mask stacks only when asked (its `write2disk` was computed and
  ignored). The policy: `doc/policies/2D/solve2D_policy.md`, the centring stage's subsection.
- Gates: Debug (run 259) and Release (run 265) fast gates 15/15; the suites pose statistics 135, UI visibility 380,
  masks 74, polar shift map 65, refine3D in-plane state 101 checks on Debug (runs 260-264).
- The five-view fixture (400 particles, 80 per view, Gaussian shifts of sigma 7 px) through solve2D itself from its
  random start, ncls 10, the six stages, 12-16 s per run on 8 threads (the key costs nothing measurable at this
  size), four runs per variant; recovery = the fraction within lp/4 at lp 12 (2.3 px) after the per-class offset:

  | variant | purity per run | recovery per run | mean error (px) |
  |---|---|---|---|
  | key, mskdiam=0 (the estimate, 187-232 A, for the later stages) | 60, 60, 80, 100 % | 37, 99.5, 100, 100 % | 2.4, 1.2, 1.1, 1.2 |
  | key, mskdiam=180 | 100, 100, 20, 40 % | 100, 100, 49, 60 % | 1.2, 1.1, 4.8, 2.9 |
  | key off, mskdiam=180 (today) | 58, 40, 61, 98 % | 47, 64, 59, 97 % | 3.6, 2.5, 2.1, 1.5 |

- Readings. (a) With the key the recovery is complete whenever the classes are not degenerate, merged views included
  (two views per class centre their members alike: purity 60 %, recovery 99.5 %); the window recovers in proportion
  to the purity. (b) The classification from the random start is unstable in every variant on this fixture: the
  collapse into one or two classes (purity 20-40 %) happens in stage 1, before any shift search (run 296: class
  overlap 1.000 from stage 1's second iteration), and nothing later separates a collapsed class; it struck the key
  twice and the window once in four runs each, so it is the fixture's (80 particles per view, ten classes), not the
  key's. The box-wide first stage of the mask-free mode did not hurt. (c) The shift exit never ended the centring
  stage: the tempered draws keep the median shift change above lp/4 while the stage runs (in the collapsed run the
  increments spread over the whole disk, mean 12 px, max 22); it fires from stage 4 on, once the draws sit on one
  cell (median change 0.00). Ruling 8's statistic measures the draw's spread, not the posterior's convergence; the
  stage ran its five-iteration budget. A statistic of the posterior itself (the drawn cell's mass, or the map's
  entropy averaged over the sampled particles) is the candidate for a real exit; left as it is until Phase 3 says
  the budget is wrong. (d) The estimate of the mask diameter, 187-232 A for particles simulated under a 180 A mask
  (molecular diameters 156-193 A times 1.2), is usable as the later stages' mask.

### 2026-10-31: Phase 2b, solve3D and the key's scope (Claude, agent kit runs 314-357)

- The key through solve3D's stages (stage 3, the first with a shift search, opens the disk: `trs` raised to 3 sigma
  in native pixels; the stages after it map their window) on the random-view fixture (400 particles, shifts of mean
  9.2 px) from solve3D's own random start, key against off, the recovery after the fit of the map's 3D offset
  (`compare_shifts3D.py`): the first run with the key drifted in stages 6-8 (21.8 % within lp/4, mean error 14.6 px,
  map 8.1 A; the window: 100 %, 0.56 px, 5.6 A). Stages 3-5 were healthy (stage 3's first iteration moves the
  particles 4.9 px on average, 18 px at most, the disk's reach, against the window's 2.8 and 5.7; stage 5 ended at
  5.8 A). The greedy cell and a narrow prior (sigma 4 A) drifted the same, so not the draw; in isolation from the
  converged project at lp 7.5 A and trs 8 px (runs 332-338): `refine=prob` with the key 100 %, `prob_neigh` with the
  key 71 % with increments growing to 11.3 px (the window's corner), `prob_neigh` without it 100 %. The cause: the
  neighbourhood table has its own fill (`fill_tab_neigh_range`) without the map, while `strategy3D_prob` under the
  key treated every row's increment as the map's cell (the frame conversion at the drawn rotation, the score at that
  pose). Fixed by scope, not by a second implementation: the neighbourhood fill refuses the key, `prep_strategy3D`
  refuses it for every `refine` but `prob`, `prep_strategy2D_glob` for every 2D `refine` without a table, and the
  controllers hand it only to those stages (solve3D stages 3-5; solve2D's table stages, the terminal greedy pass
  dropping it). Rebuilt (Debug run 341, Release 346, both gates 15/15; the four suites 135/65/101/380): solve3D
  key 100 % (0.56 px, final map 5.5 A) and off 100 % (0.57 px, 6.2 A), 233 against 263 s. This fixture's shifts are
  within the window's reach over solve3D's hundred iterations, so it shows the key harmless, not better; the disk's
  case is Phase 3's data. The map in the neighbourhood table is the item to build when the later stages need it.
- The policies: `doc/policies/3D/refine3D_policy.md` (the map's paragraph beside the shift prior's),
  `solve3D_policy.md` (the key's stage scope), `solve2D_policy.md` (the refusal outside the tables).
- Standalone refine2D/refine3D: the disk is min(3 sigma, trs), so the full disk needs `trs` of at least 3 sigma in
  native pixels (the e2e scripts' map variants now use trs=24 at sigma 10 A and 1.3 A per pixel).

### 2026-10-31: Phase 3, first data: the beta-gal bench (Hans's runs 8_solve2D and 9_solve2D, Claude's reading)

- `solve2D ncls=90 mskdiam=0 nparts=4 shift_sigma=10` on the bench's 4463 particles after the class rejection (box
  256 at 1.275 A, cropped to 120; stage 2 at trs 23.5 px under the box's mask, the estimate 261 A adopted) against
  run 8 on the same particles (taken as the window control, settings to be confirmed) and the original run 2 on all
  5513 particles. The picks are off by 8.3 px on average (10.6 A), 15 % by more than 13 px, 1 % by more than 22 px.
- The shifts: the key run's and the control's agree per particle to 1.2 px median (80 % within 2 px), the same as
  two window runs agree with each other (run 2 against 8: 1.1 px median), and the three distributions are the same
  (mean 8.3-8.5 px, 90th percentile 14.9-15.2 px). On this set the window reaches every shift over its 27
  iterations, so the disk changes nothing and the map reproduces the window's answer.
- The class averages, same particles: key 55 populated classes, best 6.16 A, 7 % of the particles in classes at
  6.5 A or better, 32 % at 8 A or better, population-weighted mean 10.57 A; control 54 classes, 5.63 A, 11 %, 39 %,
  10.98 A. A weaker top end under the key, the mean the same; the 261 A mask against 180 (the class resolution is
  an FRC under the mask) and the later stages' tempered shift draws (increments of 1.5 px on average through stages
  3-5, 0.1 after the anneal; the window's polish does not jitter) are the two candidates, separable by a key run at
  `mskdiam=180` and by `shift_draw=no`. Both runs are far below run 2 (84 classes, 4.29 A, 23 % at 6.5 A), which is
  the particle set or the run, not the key, and is a question for Hans.
- The log: the centring stage's shift exit never fired (median change 2.6-6 px against lp/4 of 2.75 px); it fired in
  stages 3-5 without effect (minits is the budget there) and in stage 6; the table's cost at the disk 10 s per
  iteration against 3 s at the window and 7.5 s dense, on 6 threads; the sparse stages' class overlap near zero
  until the dense stage is the sparse draw's normal reporting (the fixture's window runs show it too).
- Scripts: `~/agent_work/SIMPLE/scripts/shift_stats.py` (the magnitude histogram, per-class offsets), the fields in
  `~/agent_work/SIMPLE/bgal_ptcl2D_*.txt`, `bgal_cls2D_*.txt`.

### 2026-10-31: integration re-review

The follow-up in `polar_shift_search_review.md`, section 12, reviews the
kernel, sampler, controllers and the uncommitted integration batch present
at review start. Its cutoff and the later live 2D edits are recorded there.
Six issues remain open: the new extra-pose conversion double-rotates window
entries (reopening the frame conclusion of Open decision 3); map entries
reach the insertion cone in different reference frames; weighted insertion
uses marginal evidence as a point-pose loss; a disk narrower than a grid
step can exceed its requested radius; the scheduled continuous pass inherits
a now-refused table-only key; and map-mode posterior diagnostics write
default rows. The numerical counterexamples, source paths and validation
are recorded in the review. No production changes were made by that review.

### 2026-10-31: the cost on the beta-gal bench, the screen-rank record, the mask estimate fixed (Claude, runs 369-378)

- Hans: the key doubles the run time on the sweep; the automatic mask hurts (his class docs: key with mskdiam=180
  best class 4.35 A, eight at 5.4 A or better, as the original run; mskdiam=0, the estimate 261 A, best 5.5 A).
- The estimate overshot because it went through `automask2D_mskdiam`, the picking box's rule (diameter + 2 soft
  edges, times 1.2): the population-weighted 90th percentile of the class diameters was 185 A for beta-gal (user's
  mskdiam 180). Fixed: the estimate is the diameter itself, capped at the box; it now prints 183-185 A.
- The cost, same particles (4463, ncls 90, 4 parts x 6 threads, this Mac, load 15-20 from other jobs but within
  5 % of Hans's idle records), per stage in seconds: window 8.6 / 13.0 / 17.5 / 17.5 / 17.8 / 23.3 = 100.6; key with
  mskdiam=180: 8.0 / 72.8 / 19.5 / 20.2 / 8.3 / 33.2 = 165.5; key with mskdiam=0 (lp ladder from 20 A): 8.6 / 58.8 /
  22.5 / 22.7 / 8.8 / 37.6 = 162.5. The centring stage is +60 of the +65 s; the later stages' maps are about
  neutral (stage 6 +10, stage 5 -9 from fewer iterations).
- Inside the centring stage (part 1, 1195 particles, per iteration at lp 11.5 A): fill 11-12 s wall, score loop
  45-55 thread-s (the origin evidence and the screen map of every scanned class, about 70 per particle: 0.65 ms
  each against 5.6 us for a plain score), the full maps 11-24 thread-s (9.6 to 4.3 per particle, 2.2 ms each). At
  the window radius a map costs 0.24 ms, a classic L-BFGS minimisation 0.27 ms (9 evaluations over all rotations,
  11 gradients, from a 25-point lattice seed), a screen 90 us: the screen over every class, not the map, is the
  overhead of the later stages.
- The screen-rank record (new in the 2D cost record: the origin rank of the class the screen puts first, bins
  rank 1 / within the refined count / 2x / 4x / beyond): at the first centring iteration 5.9 % of the particles
  beyond four times the refined count and 17 % beyond twice; 2.4 % beyond 4x at the second iteration, 1.1 % at
  the fifth, 0.3 % from stage 3 on. A pool of the origin's best four times the refined count would halve the
  screens and miss those.
- Levers, with the expected saving on this bench: a cheaper screen (half the shells or rings, every class kept)
  about -30 % of the stage; three centring iterations instead of five, -40 %; the pool, -30 % at the 6 % cost
  above; the user's mask plus the disk in place of the box's, about -20 % (ruling 7); dropping the later stages'
  maps, about -5 s. Awaiting Hans's choice.
- Class quality of the three runs (single runs, the random start dominates): window best 4.24 A, 14 % of the
  particles in classes at 6.5 A or better; key 180: 4.60 A, 23 %; key with the estimate (185 A): 4.41 A, 10 %.
