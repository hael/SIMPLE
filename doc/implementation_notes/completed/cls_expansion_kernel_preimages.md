# cls_expansion: sub-class averages as kernel pre-images (`subcls_avg=kernel`)

    simple_exec prg=cls_expansion projfile=... ncls=<total subclasses> mskdiam=... subcls_avg=kernel nparts=10 nthr=8

Status (2026-10-14): implemented as an optional route; the default stays `subcls_avg=model`
(`flex_cls_expansion_2026_10_01.md`). Unit-tested, checked end to end (shared memory and two parts) on the
three-state phantom of the `flex_pca_blobs` gate. Revised the same day after the first real-data run: the
first version kept a class whole unless the CV gain was significant and capped the representatives by `minpop`
and noise overlap, and delivered 54 of 500 budgeted averages. The route is a generative model (a pre-image per
member), so it now fills the budget with the pre-images of that many members (section 2, steps 4-5; section 7).

## 1. Why

Since 2026-10-12 the sub-class averages are the model's conditional averages,
`P_h + U_ext,h (I - S_h)^-1 (zbar_s,h - zbar_h)` per subclass `s` and half `h`. `solve3D_cavgs` and
`refine3D_states_cavgs` read each of them as two independent rows (even, odd), plan their low-pass ladder
from the registered `frc2D`, and build state maps from the rows. That breaks three of their assumptions:

1. **Shared noise.** Every sibling carries the same half parent `P_h`, i.e. the same noise, plus a
   combination of the same `neigs` noisy basis columns. Siblings at similar poses add the same noise in the
   same place: the 3D solve sees structure that no particle subset holds. A class with fewer than
   `2 * 3 (neigs + 1)` members has no half models at all and uses the joint deviations in both halves: there
   the even-minus-odd images of all siblings are identical (section 5).
2. **Borrowed resolution.** The even/odd FRC of `P_e + D_e` against `P_o + D_o` is essentially the parent's
   (the deviation is small next to the mean). The ladder and every resolution claim of a sub-average, and of
   state maps built from siblings, are the parent's.
3. **Exaggeration.** The members are selected on their posterior means `m`; by the tower property
   `E[z | m in cell] = E[m | m in cell]`, so the deshrink `(I - S)^-1` inflates every deviation, and `P + U z`
   at an inflated `z` extrapolates a variation that need not be linear into an image no object projects to.

Plain hard averages are no fix: the latent is a function of a member's own noise along `span(U)` in the fit
band, so selecting on it averages selected noise ("Einstein from noise").

The argument for hard assignment over soft weights (with many subclasses the soft averages converge to each
other) holds for responsibilities normalised across components. A kernel weight is a responsibility without
that denominator: it does not shrink when representatives are added. That decouples the weight distribution of
a particle from the clustering.

## 2. Method

Per 2D class with at least two subclasses in the budget, after the class model of the default route is fitted
(bands, nuisance columns and all; `flex_cls_expansion_2026_10_01.md` section 1, steps 1-2):

1. **Coefficient split.** The global fit band is split into two interleaved sets of 4x4 Fourier blocks in a
   checkerboard. Blocks, not single coefficients: the class-frame soft mask
   correlates direct neighbours at about 0.2, which would leak the selection noise into the validation. A
   model is fitted on each set (`flex_cls_fit`, the class's nuisance columns, rank `neigs`, lower when a set
   holds fewer than `3 (neigs + 4)` coefficients; a band too narrow for any rank keeps the class whole).
2. **ML latents.** Every member's maximum-likelihood latent `zhat_i` under each of the three models, without
   the latent prior and with the nuisance profiled out, and its covariance `L_i^-1`, the latent block of the
   inverse of the joint normal matrix.
3. **Directions.** Per model, the generalised eigenvectors of the members' latent scatter against their mean
   ML covariance, normalised so that the mean covariance is the identity (noise-whitened), ordered by scatter
   over noise. The in-sample fit makes every direction look like signal in
   the members' own scatter (a rank-16 basis fitted to noise sits at the top of the Marchenko-Pastur spectrum
   of the coefficients, far above the edge for 16 latent dimensions), so the rank is not read from these
   eigenvalues but chosen by the cross-validation.
4. **Kernel and cross-validation.** A member at location `c` weighs `w_j(c) = N(zhat_j; c, L_j^-1 + tau^2 I)`,
   its likelihood of being there, in the leading `r` directions; an uninformative member weighs little
   anywhere. For `r in {1, 2, 4, 8, 16}` (up to the rank of the set models) and
   `tau in {0, 0.25, 0.5, 1, 2, 4, 8, inf}` (noise units; `inf` is the class mean) the members selected on set
   A's latents predict every validation member's set-B coefficients, leaving the member out, as the CTF-Wiener
   kernel mean `sum_j w_j c_j W_j y_j / (sum_j w_j c_j^2 W_j + eps)`, and vice versa (the
   target's own latent noise enters as the identity, the class's mean covariance). The validation noise is
   independent of the selection, so the whitened error per held-out coefficient is 1 plus the pre-image's
   error. The rank `r` with the smallest error at a finite bandwidth is kept, with the CV optimum `tau*` there
   (`inf` included) and the best finite bandwidth `tau_g` (`tau*` itself when finite). The paired per-member
   gain of `tau_g` over `tau = inf` gives a z-score, reported (log, `sep`) and not acted on: below about 3 the
   class shows no reproducible variation at its SNR and its pre-images are near copies of its mean.
5. **Pre-images and representatives.** On the full model's latents in its leading `r` directions. Member
   `i`'s pre-image is the kernel mean of the members with weights `w_ij = N(zhat_j; zhat_i, L_j^-1 + C + tau_g^2 I)`,
   `C` the members' mean covariance (the latent uncertainty of the centre, as in the CV): the generative model,
   one image per member, at the bandwidth the CV chose (`preimage`). The representatives are the
   pre-images of `nsub` members, the class's budget (`select`). Candidates are the
   informative members (an even stride of at most 1000); they are taken greedily by coverage: each step the
   candidate that adds most to `sum_j max_r L_rj`, `L_rj` the members' likelihoods under the candidate's kernel
   (a submodular facility-location objective: the densest region first, then what the accepted ones leave
   uncovered, so a second state comes before the gap between two and a near copy of an accepted pre-image comes
   last), first among the candidates whose pre-image holds a Kish effective population `(sum w)^2 / sum w^2` of
   at least `minpop`, then among the others. The budget is delivered whenever the class has that many members.
   Kept whole: a class with a budget of one, a fit band too narrow for the coefficient split, or latents that
   give no CV (too few informative members, a covariance not positive definite).
6. **Averages.** Each representative's full, even and odd average is the kernel-weighted, CTF-corrected
   mean of all, the even and the odd members with the same weights (`average`). On the fitted
   coefficients the latent is a function of the member's own noise, so the weighted mean is biased towards
   the representative there. Stein's identity removes it: for Gaussian noise
   `E[w(zhat) e_k] = -c_k w (U_k . L^-1 g + N_k . S_nz g)`, `g = W (L_r^-1 + C + tau_g^2 I)^-1 (zr - c)`, `W` the kernel's
   directions (`zr = zhat W`), `S_nz` the nuisance-latent block of the inverse; the sample version is added to
   the numerator. Above the fit band the noise never entered the latent and nothing is corrected. The full
   average is filtered by `2F/(1+F)` of its own even/odd FRC, as on the default route.
7. **Labels and readouts.** A member's hard label is the representative of its largest weight (the
   `ptcl2D` contract is unchanged). Per subclass: `pop` (hard members), `neff` (Kish), `sep` (the parent's CV
   z-score) and `repro`, the even/odd correlation of the deviations from the half parents over the fit band
   (the default route's is over every coefficient; beyond the band the deviation of a small subclass is noise).

What the 3D programs get: `ncls` rows, FRCs of the representatives' own particles, deviations supported by
the members, and noise shared between siblings only through the members both weigh. That sharing is not capped
any more: a budget larger than the class's variation supports gives near copies (the largest overlap is logged
per class; `neff` and the z-score say how much each row and class carry).

## 3. Readouts in the log

    Cls expansion kernel CV: class=<c> whitened error per held-out coefficient by latent rank 1 2 4 8 16 (best bandwidth; -1 not tried): <5 values>
    Cls expansion kernel CV: class=<c> at rank <r>, tau = 0 0.25 0.5 1 2 4 8 (noise units) and the class mean: <8 values>
    Cls expansion kernel: class=<c> members=<n> set rank=<k> kernel rank=<r> CV tau=<tau*> pre-image tau=<tau_g> CV z (pre-images over the class mean)=<z> representatives=<m> of <budget> max overlap=<rho>
    Cls expansion kernel subclass: class=<c> sub=<s> pre-image of member <i> pop=<hard members> neff=<Kish> repro=<even/odd> FRC 0.143=<A> A
    >>> CLS_EXPANSION KERNEL: delivered <total> of the <ncls> budgeted subclasses (classes kept whole: ...)

The last line appears only when a class could not be split (see step 5).

The CV curves are the evidence: a class whose curve falls to the class mean has no reproducible variation at
the class's SNR; the rank curve says how many latent directions carry it.

## 4. Implementation

| layer | file | what |
|---|---|---|
| params | `simple_parameters.f90`, `_parse.f90`, `_phases.f90` | `subcls_avg` (`model|kernel`, default `model`), validated |
| ui | `simple_ui_refine2D.f90` (`cls_expansion`) | `subcls_avg` (advanced), the program description |
| strategy | `simple_cls_expansion_strategy.f90` | `split_class_flex` routes after the class model to its internal `deliver_kernel`: the class's `flex_cls_kernel`, its averages, FRCs, filter and images; a class reports the subclasses it delivered and their `neff`; on this route a class kept whole (step 5) is reported, not an error. The image and resolution helpers are shared with the default route (`flex_planes_to_images`, `flex_res143`) |
| domain | `src/main/flex/simple_flex_cls_kernel.f90` (new) | the class `flex_cls_kernel`, private state: `new` (ML latents, directions, the two set models, the CV), `select` (the representatives), `average`, `preimage` (any member's), `labels`, `neff`, `report` (the log lines), `kill` |
| test | `src/main/flex/simple_flex_cls_kernel_tester.f90`, a white-box submodule like the PCG self-test (suite `flex cls_expansion kernel` in `unit_heterogeneity`, built only with `BUILD_TESTS=ON`) | section 5 |

The distributed merge needed no change: it counts the rows of the part maps, so classes with fewer subclasses
than budgeted merge as they come (checked: two parts reproduce the shared-memory class map exactly).
`scripts/check_flex_dag.py` ignores the new module as it ignores `simple_flex_cls_expansion` (outside the
FLEX-PCA layering).

Cost per class with `N` members, at most 1000 validation members and 1000 candidates: two more model fits on
half the band; for each of up to five ranks and seven bandwidths `2000 N (r^2 + 3 n_q)` for the CV (`n_q`
held-out coefficients; negligible weights are skipped); `1000 N r^2` for the candidates' pre-images and
`1000 N` per representative for the coverage; the averages as on the default route. On the phantom below (60 classes of 60, box 64) the kernel
route took 12 s against 10 s for the default route on two threads.

## 5. Tests and checks

**Unit tests**, `simple_test_exec test=unit_heterogeneity suite=flex_cls_expansion_kernel` (65 checks, about
2 s). Synthetic classes on the half-plane of a 32-pixel box, CTF-like factors per member and coefficient,
complex Gaussian noise of unit power:

- the coefficient sets are disjoint, cover the band, balanced, constant on a block;
- noise-free members: the ML latent is the generating latent and its covariance the inverse information
  (closed form), also with a profiled nuisance column;
- two clusters of latents, a budget of four: four representatives, each at a member's latent, the first two
  at the cluster centres, every member labelled with a representative of its cluster, every one holding
  `minpop`; the budget still filled when no pre-image reaches `minpop`; with one allowed, the denser cluster;
- a homogeneous class weighted around an off-centre location: the plain weighted mean deviates along the
  latent direction by `c cov / (2 cov + C + tau^2)` (0.500 measured, 0.5 expected with `C = 0`; 0.334 and 1/3
  with `C = cov`); the Stein-corrected one by -0.005 and -0.002;
- a homogeneous class: no significant CV gain (z below 3), still four pre-images at a finite bandwidth; with a
  budget of one it is kept whole and delivers its CTF-corrected mean, also as every member's pre-image;
- two states (60/40, rank-2 models, budget two): z at least 3, kernel rank 1, the pre-images at the CV
  optimum, labels at least 90% pure (the best possible is about 96%), the averages on either side of the class
  mean along the state direction;
- the same two states with a budget of six: six pre-images of members, each with at least `minpop` effective
  members (120-260 of 400), both states represented, every representative's average equal to its member's
  pre-image.

Over five further seeds all pass: the Stein residual stays within 0.045 (threshold 0.12), the homogeneous
class's CV z-score at most 0.9 (threshold 3), purity 0.95-0.97.

**Three-state phantom**, first version of the representatives (the z gate, the overlap cap and the `minpop`
floor; not re-run after the revision) (`flex_pca_blobs`' `states_cavgs` fixture: 60 views, each an exact 2D class of 20
particles per state, SNR 1, box 64 at 2 A, one 4 A lobe differs between states; `ncls=180 lp=12`). Sibling
noise: correlation of the even-minus-odd images of two sub-averages of one class inside the mask. Deviations:
`c (A_s - class mean)` against the noise-free truth of the sub-average's members, below 21 A.

| | label purity | sibling noise correlation | sub-average FRC 0.143 | deviation correlation / slope / rel. error |
|---|---|---|---|---|
| `model`, `minpop=20` | 0.685 | 1.000 (identical) | 9.85 A (the parent's) | 0.51 / 0.82 / 1.96 |
| `kernel`, `minpop=10` | 0.994 | -0.004 | 12.8 A (20 particles' own) | 0.64 / 0.96 / 1.32 |

The kernel CV chose rank 2 (three states span two directions) in 59 of 60 classes. With `minpop=20` the
first version kept 53 of 60 classes whole: a state of exactly 20 particles cannot hold a weighted subclass of
Kish population 20. That floor, the z gate and the overlap cap are what made the first real-data run deliver 54
of 500; the revised rule only prefers pre-images that reach `minpop`.

## 6. Open

- **Real-data checks.** On a class split by both routes: the correlation of the even-minus-odd images of two
  siblings inside the mask (near 0 expected here, high on the default route), and `refine3D_states_cavgs` /
  `solve3D_cavgs` on both outputs (fake states, resolution claims). `solve3D_cavgs` cannot run on the phantom
  (its ladder asks for a crop box of 88 at a native box of 64).
- **Redundant rows.** A budget beyond what a class's variation supports now gives near copies with shared
  noise (logged overlap near 1). If `refine3D_states_cavgs` / `solve3D_cavgs` count them as independent
  evidence, a weight per row (its `neff` over the shared part) or an overlap-aware budget is the next step.
- **Calibration.** The kernel and the Stein term take the ML covariances at face value; a sigma2 off by a
  factor scales them. The CV picks rank and bandwidth regardless; the Stein term would be off by the same
  factor.
- **Junk.** The default placement uses the log fit residual as a junk axis; here junk members only lower
  their weights through their covariances. A residual coordinate in the kernel is a candidate.
- **Budget on capacity.** The budget is allocated before the classes are seen; the delivered counts could
  feed a second allocation.
