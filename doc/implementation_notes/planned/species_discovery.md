# SINGLE without a species input: atom detection, Gaussian density model, species discovery and map simulation

Implementation note, 2026-10-07, revised 2026-10-08 after review. Written
against master `108c8c9d`. Companion page: `single_pcg_integration.md`
(section 11 says how the two fit together).

Status. The design is for review and is not approved as a replacement of the
present detection. What is approved is a small, opt-in prototype: the
foundations and the species-free level 1 of section 10 (phases 0 to 2, which
change nothing for an existing run) and the residual recovery and species
call behind an experimental key, producing diagnostics only (phase 3). The
prototype is developed alongside the existing detection machinery, which
keeps working unchanged throughout; section 10.0 states the rules. Phases 0
to 3 ran on 2026-10-08. The second run, phases 4 and 5 (10.5, 10.6), follows
two rulings of the same day: B factors go into the element kernels of the
simulator, and the species of a synthesised particle are given on the
command line as a list; its test renders ground truth with the element
kernels, not with the atoms the detector assumes. Promotion of anything
beyond diagnostics waits for the validation of section 10.7. Section 15
records the reviews and what they changed.

Nothing in SIMPLE was compiled or run. The numbers in section 8 come from a
numpy emulation of the proposed pipeline on synthetic particles
(`emul_species.py`, which needs `emul_detect.py`; both are kept out of the
tree, like the other reference scripts), checked in section 8.8 by a second,
independent implementation (`codex_checks.py`). They show how the algorithm
behaves on idealised input and are not a measurement on data: the generator
produces the Gaussian atoms the detector assumes, so agreement between the two
says the design is coherent, not that it works on a reconstruction.

## 0. What is asked, and what constrains the design

The request: no `element` on the command line. The program finds the atoms,
finds how many kinds of atom there are and which atom is which, and simulates
a map from what it found. The first data are HAADF-STEM images (high-angle
annular dark-field scanning transmission electron microscopy, whose contrast
rises with atomic number) of Pt/Ni particles.

Two further points from the discussion shape the design.

- The B factor rises from the core to the surface, so anything that matches
  the density against a fixed shape has a different error rate at the surface
  than in the core. The density should be modelled by a generic Gaussian fit
  whose parameters can be examined as a function of radius. This is adopted:
  every atom gets its own fitted width, the species call does not use
  matching, and the radial profiles of the fit are a product (sections 1, 3.7,
  3.8, 7).
- An alternative is to specify the species, take the maximum over per-species
  correlation maps for the segmentation, and assign species with a Gaussian
  mixture model that includes the B factors. The mixture model with embedded
  B factors is what sections 3.6 and 3.7 are. The maximum over per-species
  maps does not change the segmentation, and section 5 says why and what
  specifying the species does buy.

"Species" in this note means an intensity class. HAADF contrast orders the
atoms by atomic number, so the data can say that there are two kinds of atom,
which is which, and what their intensity ratio is. They cannot name the
elements without a calibration, and nothing in the density model needs the
names (sections 5 and 12).

"Pt:Ni is 3-6:1" can be read as the atom ratio or as the intensity ratio.
The design does not depend on the reading: the emulation covers atom ratios of
3:1 and 6:1 crossed with intensity ratios of 3:1 and 6:1. Which one is meant
matters for what to expect from the data (Q1).

Constraints taken from the code:

- `refine3D_nano` has no species dependence
  (`single_commanders_nano3D.f90:269-295`). The species enters
  `autorefine3D_nano` only through `detect_atoms` and the `_SIM.mrc` it
  returns as `vol1`, and once upstream, in the `startvol.mrc` of
  `analysis2D_nano` (`single_commanders_nano2D.f90:196-203`).
- The simulated map is an alignment reference only. Trailing reconstruction
  blends accumulators, never volumes
  (`.github/skills/simple-frac-update-trailing/references/frac-update-contract.md`),
  so `_SIM.mrc` puts no density into the reconstruction. The intensities
  measured in the next map come from the particle images. (As read in
  `derive_sampling_settings`, `simple_parameters_phases.f90:504-512`,
  `l_trail_rec` is cleared when `update_frac` is above 0.99, and
  `autorefine3D_nano` sets `trail_rec=yes ufrac_trec=0.5` but not
  `update_frac`. If that reading is right, trailing is off in the default
  SINGLE loop. Not traced beyond that routine.)
- The volumes `refine3D_nano` ships are zero outside a soft sphere of radius
  `msk_crop` (`simple_commanders_rec_distr.f90:461-466`, `:638-682`), and
  `nanoparticle%new` masks again (`simple_nanoparticle.f90:271-278`). Any
  noise reference has to come from inside that sphere and away from the
  atoms.
- `image%phase_corr` divides by one global norm that depends on the map
  (`simple_image_calc.f90:1482-1485`). Correlation maps of different residuals
  are therefore not on one scale until they are standardised.
- `atoms%convolve` reads neither `beta` nor `occupancy` and silently skips
  any `Z` outside its table (`simple_atoms.f90:1235-1236`); `set_name` and
  `set_element` stop on an unknown symbol (`:727-743`, `:890-904`).
- A single-species data set must come out as it does today: one class, no
  atoms added.
- `discard_atoms` deletes an atom by clearing its component in `_BIN`,
  relabelling `_CC` and calling `find_centers`, which deallocates and
  reallocates `atominfo` (`simple_nanoparticle.f90:767-768`, `:1116`,
  `:1151`); `identify_atomic_pos` already has to recompute `valid_corr`
  after it for that reason (`:601-603`). Nothing fitted per atom survives
  pruning unless pruning comes first or works on the atom table. Its size
  test compares the longest extent of the component with half the atom
  radius (`:1132-1139`), so a component that is generated from the model
  can never fail it.
- The regression test of the present path is `single_atoms_stats`
  (`production/CMakeLists.txt:307-308`): it simulates a Pt particle in
  process (`simulate_nanoparticle`, box 160, `moldiam=20`), runs
  `detect_atoms element=Pt` and `atoms_stats`, and requires recall and
  precision of 0.90 or more within 1 A, a root-mean-square position error
  of at most 0.75 A and a correlation of at least 0.80 between `_SIM.mrc`
  and the input. The last recorded result is recall 1.000, precision 1.000,
  0.015 A, 1.000. The sub-suite `nanoparticle atoms` of `lib_single` that
  the first version of this note named no longer exists (moved 2026-09-25).

## 1. The model

    rho(r) = b(r) + sum_i A_i exp(-4 pi^2 |r - r_i|^2 / B_i)

`b` is a smooth background. Each atom is one isotropic Gaussian with its own
position, amplitude and width. Its integrated intensity is

    I_i = A_i (B_i / 4 pi)^(3/2)

The atoms form a mixture of `K` classes. Class `k` has an intensity `I_k`,
and an atom of class `k` has `I_i` near `I_k`. The width `B_i` belongs to the
atom and not to the class: it is free per atom, and it is expected to vary
with radius.

The reason for this split is physical. In an incoherent image the integral of
an atom is set by its scattering power, and displacement, the probe and the
reconstruction only spread it. So `I` carries the species and `B` carries the
environment. A statistic that mixes the two, such as the peak height, a
correlation with a fixed shape, or the amplitude at a fixed width, cannot
tell a light atom from a heavy atom that is smeared. Section 8.4 shows that
this is not a small effect: with widths that grow towards the surface, the
amplitude at a fixed width invents a species in a single-species particle.

Unknowns found from the map: the number of atoms, `r_i`, `B_i`, the number of
classes `K`, the class of each atom with its probability, and the `I_k`.

## 2. Length scales from the data

The element currently supplies an atom radius `r` and, through the lattice
table, a neighbour cutoff. Both follow from one measured number, the
nearest-neighbour distance `d_NN`: the median over atoms of the distance to
the nearest other centre, taken from the centres of the first binarisation.

- `r = 0.4 d_NN`. For Pt this gives 1.106 A against the tabulated 1.10 A, so
  `split_atoms`, `validate_atoms` and `discard_atoms` keep their behaviour.
- Neighbour cutoff `1.267 d_NN`, which is what `find_rMax` returns for FCC Pt
  (3.504 A). This assumes close packing. The general form is the midpoint
  between the first two peaks of the pair-distance histogram plus `0.15 r`;
  not emulated (Q5).
- The contact-score ceiling of `discard_atoms` (12 or 4, from the crystal
  system) is 12 in the prototype, the close-packed value that goes with the
  neighbour cutoff above (Q5), so that the pruning logic is unchanged.
  Deriving it from the data (the median contact score of the inner 30% of
  the atoms) belongs to 10.7.
- Fit and aperture radius `d_NN / 2`.
- Radial coordinate of an atom: its distance from the centre of the atom
  positions (`cendist`, already in `atom_stats`). For a particle that is far
  from spherical the depth below the surface is the better coordinate (Q3).

The distance tests of `split_atom` multiply squared voxel distances by `smpd`
where `smpd^2` is meant (`simple_nanoparticle.f90:914-969`), so the exclusion
radius written as `0.9 r` acts as 0.59 A for Pt at 0.358 A. The species-free
path should state it in Angstroms as a fraction of `d_NN`: 0.21 reproduces
today's effective value, 0.36 the written one. To be settled on the Pt test.

## 3. Pipeline

Nine parts inside `detect_atoms`; 3.9 gives the order they execute in (level
1, its pruning, then the rest). 3.1 and 3.9 are the present flow, with the
Gaussian kernel and the `d_NN` length scales substituted only when `element`
is absent. The rest is new and runs only under `discover_species=yes`.

### 3.1 Level 1: the present flow with a pseudo-atom

`phasecorr_one_atom`, the Otsu threshold ladder, the bisection,
`discard_small_ccs` and `split_atoms` run as now, with one change: the
template and the equal-atom simulation of `t2c` use the Gaussian of section 1
at a reference width `B_ref`, not an element. At 0.358 A the Pt template
correlates at 0.997 with a single Gaussian of sigma 0.42 A (`B` = 13.9 A^2),
which is the starting `B_ref`; with an element, `B_ref` is the Gaussian
fitted to that element's template. After the first fit (3.4) the median `B`
of the strong atoms is compared with `B_ref` and the ratio is reported
(`_species.txt`); re-running the detection at the fitted width is left to
10.7, because in the prototype level 1 and its products must not depend on
the opt-in key. Nothing is carried between calls of `detect_atoms`.

Level 1 finds the strongest class and nothing else. That is the behaviour to
preserve for single-species data, and the reason more levels are needed: the
Otsu gate sits at about a quarter of the strong peak, and in every
multi-species case of section 8 level 1 found the strongest class (all of it
in all cases but two, 395 and 394 of 396 there) and none of the weakest.

### 3.2 Noise reference

The noise region is the set of voxels inside the shipped sphere, at least the
soft-edge width from its rim, and farther than `1.5 d_NN` from every accepted
atom. It is rebuilt after each level. Three numbers come from it:

- the mean and standard deviation of each correlation map, which turn that
  map into a z map and remove the global norm of `phase_corr`;
- `sigma_n`, the standard deviation of the map itself;
- `s_A`, the standard deviation of the matched amplitude at `B_ref`, and
  `s_I`, that of the aperture intensity (3.5), both measured at 400 or more
  phantom sites placed at random in the region.

Measuring the noise of each statistic by applying that statistic to empty
positions makes no assumption about the colour of the noise. The emulation
used white noise, noise band-limited at 1 A, and strongly correlated noise,
with the same thresholds.

If the region is too small, because `mskdiam` hugs the particle, the program
stops and asks for a larger mask. The noise in the region is a lower bound on
the noise inside the particle, where alignment error adds to it.

**Half maps as the noise reference.** With `vol_even` and `vol_odd` given,
`sigma_n` and the mean and standard deviation that standardise each
correlation map come from the half-map difference (its noise is sqrt(2)
times that of a half and twice that of the average, so the full-map value is
the difference's divided by two), measured over the same region. The
difference is free of any density, so atoms that have not been found yet
cannot enter it. That matters: the region excludes 1.5 `d_NN` around every
accepted atom, which removes an undetected surface shell (its atoms sit
within one `d_NN` of the strong atoms) but not an undetected domain made only
of the weak species and thicker than that. In the check of section 8.8 such a
domain inside a mask of radius 14 to 16 A raised the region's standard
deviation by 36 to 40% (5% at 21 A); a robust spread removed only a quarter
to two fifths of that excess, because the matched filter spreads every atom
over a hundred voxels. The
direction of the error is conservative (a higher effective threshold, fewer
detections), and with a single map and no half maps it is accepted and
reported as the ratio of the region's robust to plain spread.

**The detection threshold is calibrated, not fixed.** A threshold on a z map
of local maxima is a multiple-testing problem: the number of noise maxima
above `k` grows with the volume searched and depends on the colour of the
noise, and the per-voxel Gaussian tail does not give it (it overpredicts the
count of maxima six times in the check of section 8.8). The region does. The
number of maxima of a smooth Gaussian field above `u` has the form
`C (u^2 - 1) exp(-u^2 / 2)` (the Euler characteristic density in three
dimensions), with `C` proportional to the volume searched. `C` is fitted to
the counts of region maxima above 2.5, 3.0 and 3.5, scaled by the ratio of
the search volume to the region volume, and `k` is the `u` at which the
expected number of false atoms in the search volume is a target of 0.2 per
particle. Scaling the region's own count by volume predicted the search
volume's count to within a factor 0.8 to 1.0 in the check; the calibrated `k`
came out at 5.0 for white and 1 A band-limited noise and 4.7 for strongly
correlated noise, so the fixed `k` = 5 of the first version of this note was
right for the emulated geometry and would be wrong for a larger search volume
or a different noise colour. The expected false count at the chosen `k` is a
reported diagnostic (section 7).

### 3.3 Residual levels

After level 1, repeat:

1. Subtract the fitted density of all accepted atoms from the map as read.
2. Matched-filter the residual with the same template; standardise to a z
   map with 3.2.
3. Candidates are local maxima of the z map above a threshold `k`
   (non-maximum suppression within `0.35 d_NN`) with at least `NVOX_THRESH` of
   the 27 surrounding voxels above `k`. The position is the centroid of the
   positive residual within `0.35 d_NN`.
4. Accept, fit each new atom on the residual, go to 1.

Two stages. Stage A uses the calibrated threshold `k_A` of 3.2 (5.0 on the
emulated geometry), searches the whole shipped sphere, and accepts any
candidate farther than `0.7 d_NN` from every accepted atom; it repeats until
it adds nothing. Stage B uses `k_B` = `k_A` - 1 and accepts a candidate only
if it passes a neighbour gate: no accepted atom closer than `0.85 d_NN`, and
at least three within `1.15 d_NN`. Its search volume is therefore the
accepted atoms dilated by `1.15 d_NN`. It also repeats until it adds nothing.
At most eight levels per stage. The gate's requirement is the key
`min_nbrs`, the number of already found atoms within `1.15 d_NN` that a
recovered atom needs (default 3; 0 switches the gate off), and every atom
records the stage and level that accepted it and its detection z, so that a
run with the gate off can be compared atom by atom.

Three choices here were forced by the emulation.

- Local maxima, not connected components. At a fixed low threshold the
  components of neighbouring atoms merge when the signal is strong, and the
  merged centroid lands between two atoms.
- No Otsu on the residual. It has no noise reference and can fall to the
  noise floor.
- The gate. Without it, `k_B` = 4 (its value on the emulated geometry,
  used for every number below) over the shipped sphere admits 2 to 4 noise
  maxima per particle (4.2 for white and band-limited noise, 2.3 for
  strongly correlated; 1.5 within the particle envelope; section 8.8,
  current variant of the levels; the first version of this note gave 6 to
  13 from an earlier variant), and in one single-species run of that earlier
  variant 13 of them were returned as a second species holding 2.4% of the
  atoms. That is how a species gets invented. With the gate, `k_B` = 4 gave
  0.2 to 0.4 false atoms per particle and raised the recall of weak, broad
  atoms from the 25% of stage A alone to 61% on a weak surface shell and 58%
  for weak adatoms, and from 59% to 78% on a weak-only domain.

**What the gate assumes and what it costs.** It is a local lattice prior: it
needs three found neighbours in a shell of 0.85 to 1.15 `d_NN`, so it
assumes a locally close-packed neighbourhood whose spacing is within 15% of
the particle's median. It uses only `d_NN`, so it does not assume one
lattice across the particle, but an amorphous or heavily strained surface
would fail it, and nothing in the emulation tests that. It cannot find a
weak atom with fewer than three accepted neighbours and a z below `k_A`; an
isolated adatom is found by stage A or not at all. A region made only of the
weak species is reached from its boundary inwards, one shell per level, and
in the check of section 8.8 the interior of such a domain was found to 99%.
Measured at the same `k`, the gate costs recall against an ungated search:
about six points on the surface shell and nine in the coordination 5 to 7
bin of the weak-only domain. Measured at the same false-positive rate, which
is the comparison that counts, the gated stage B at `k_B` = 4 beat an
ungated search at `k` = 4.5 by 7 to 15 points in every case and every
unsaturated coordination bin (the interior bin is 99% for both; 61 against
46% on the surface shell, 78
against 71% on the weak-only domain, 58 against 43% for weak adatoms of
coordination 2 to 4). The gate is kept, as a switch, on that evidence.

This is the one step that matches against a shape, and section 7 treats what
that means at the surface.

### 3.4 Density fit, stage 1: free amplitude

Three sweeps over all atoms. For atom `i` the data are the map minus the
background minus the fitted density of every other atom, inside a sphere of
radius `d_NN / 2`. The amplitude is linear; `B_i` is found by a
one-dimensional search that minimises

    SSE(B) / sigma_n^2 + ((ln B - m_i) / tau_i)^2

`SSE` is the sum of squared residuals over the sphere. `m_i` and `tau_i` are the median and the robust spread (floor 0.1) of `ln B`
over the strong atoms (amplitude above `10 s_A`) of the same radial shell,
within `d_NN` in radius; the global values are used when fewer than five
qualify. After each sweep the background is the residual smoothed with a
Gaussian of standard deviation `d_NN`.

A strong atom keeps its own width, since its data term dominates. A weak atom
cannot: with a free amplitude its width is not determined by its own density
(in a comparison of single-atom estimators made while drafting, not
reported here, the fitted sigma of an atom at 0.175 of
the strong amplitude ranged from 0.26 to 0.67 A around a true 0.35), so at
this stage it borrows the width of the strong atoms at its radius. Stage 2
removes the need for that.

Subtracting the neighbours is what allows the sphere to reach `d_NN / 2`
without a per-atom background term. An earlier variant, with the sphere of
`calc_isotropic_disp` (three quarters of that radius) and a free background
per atom, lost part of the integral of broad atoms to the background.

Positions were held fixed in the emulation. In the implementation each sweep
should also move the centre to the centroid of the atom's own residual.

### 3.5 Aperture intensity

    I_i = smpd^3 sum over |r - r_i| < d_NN/2 of [ map - b - sum_{j /= i} rho_j ] / F(d_NN / 2 sigma_i)

`F(t) = erf(t / sqrt 2) - sqrt(2 / pi) t exp(-t^2 / 2)` is the fraction of a
Gaussian inside `t` standard deviations, and `sigma_i^2 = B_i / 8 pi^2`. The
correction is 0.2% for the narrowest atoms of the emulation and 6% at twice
their variance. The factor `smpd^3` puts the voxel sum in the units of the
continuous integral `A_i (B_i / 4 pi)^(3/2)` of section 1, which the
stage-2 penalty of 3.7 compares it with; the first version of this note
omitted it (the mixture of 3.6 is scale-free, so the emulation did not see
the error, but the tie of 3.7 would have been off by `smpd^-3`, about 22 at
0.358 A). Every intensity in this note, in the products and in the tests is
in those units: the voxel sum of a pseudo-atom of intensity `q` is `q /
smpd^3` (section 10, phase 1).

This is the classification statistic. It is a plain sum, so it involves no
matching; it measures the quantity that carries the species; and it depends
on `B_i` only through the small correction. Its cost is noise: about 2.3
times that of the matched amplitude for white noise, and more when the noise
is concentrated at low frequency (case 8b of section 8.1).

### 3.6 Species discovery

Input: the `I_i` and `s_I`. A one-dimensional Gaussian mixture is fitted for
`K` = 1, 2, 3 as an extreme deconvolution (`simple_xd_gmm` of the clustering
library, after Bovy, Hogg and Roweis 2011): every intensity is a class value
plus measurement noise of spread `s_I`, so a class variance is the class's
intrinsic spread and its observed spread is that plus `s_I^2`. This is the
"variance floored at `s_I^2`" of the first version made exact, in the same
likelihood family; the fit is expectation-maximisation from two
deterministic starts (sorted quantiles, as in `sortmeans`, and the `K-1`
largest gaps of the sorted sample), keeping the likelier. Each fit gets the
Bayesian information criterion `BIC = -2 ln L + (3K - 1) ln N`, with `L` the
likelihood and `N` the number of atoms.

A `K` is admissible when every class holds at least `max(8, 0.02 N)` atoms and
adjacent classes are separated by `D >= 3`, with
`D = |mu_a - mu_b| / sqrt((var_a + var_b) / 2)`. `K` = 1 is always admissible.
The admissible `K` with the lowest BIC is taken. Each atom gets the posterior
probability of every class, which is stored, and is labelled by the largest.
Classes are numbered by decreasing intensity.

The labels use no spatial information. Whether the species are mixed,
ordered or segregated is the result, so no neighbour may influence a label.

`nspecies` fixes `K` when the user knows it, and a species list does the same
(section 5).

### 3.7 Density fit, stage 2: the B factor inside the mixture

With the classes known, every atom is refitted with its integral tied to its
class. For atom `i` of class `k` the one-dimensional search over `B` now
minimises

    SSE(A, B) / sigma_n^2 + ((A (B / 4 pi)^(3/2) - I_k) / tau_k)^2

with the amplitude given in closed form for each `B`. `tau_k` is the
intrinsic spread of the class from the deconvolution fit (the observed
`var_k` less `s_I^2`, under the root) with a floor of
`0.05 I_k`. Two rounds: fit, recompute the aperture intensities and the
classes, fit again.

The tie makes the width of a weak atom measurable. Its peak height is now
informative about its width, because a broader atom of the same species must
be lower by the cube of the width ratio. No width is borrowed from the strong
class any more, so a light species that is broader than the heavy one at the
same radius comes out broader (case R6 in section 8.2).

The cost is an assumption: all atoms of a class have the same integral to
within `tau_k`. An atom that is present part of the time will be fitted as a
broader atom of its class. The stage-1 widths are kept in the per-atom table
so that the two can be compared.

### 3.8 Radial profiles

The per-atom table (section 6) holds every fit parameter with the atom's
radius and coordination number, so any dependence can be plotted directly.
`detect_atoms` also writes the profiles: per class, in equal-count radial
shells of at least 20 atoms, the number of atoms, the mean radius, the width
(root mean square sigma and its standard error, and the same as `B`), the
mean intensity with its standard error, the mean coordination, and the
predicted detection signal-to-noise of section 7.

### 3.9 Pruning and validation

The order is: prune level 1, then fit it provisionally, then recover, then
prune the recovered atoms on the atom table, then fit and classify. No
per-atom fit survives a `discard_atoms` pass, so fitting starts after
level-1 pruning; the recovered atoms are pruned without `discard_atoms`,
which keeps their fits.

Level 1 is pruned as today: `discard_small_ccs`, `split_atoms`,
`validate_atoms` and `discard_atoms` run on the level-1 set, unchanged in
logic, with the length scales of section 2 in place of the element's. This
is the production path and its products (`_ATMS.pdb`, `_BIN`, `_CC`, `_MSK`,
`_SIM`) are written from it exactly as now. The first version of this note
placed pruning after the fits; the constraint of section 0 rules that out,
since `discard_atoms` rebuilds `atominfo` on every pass.

Residual-level atoms have no component at the level-1 threshold, so the
present pruning cannot apply to them and must not be imitated with
model-generated components (the size test would then be circular). They
are pruned on the atom table alone: the contact-score rule of
`discard_atoms` (fewer contacts than the threshold among the 15% of atoms
farthest from the centre), computed on the merged set, plus the exclusion
distances of 3.3, in a routine that touches neither `_BIN` nor `_CC`.

The stage-1 and stage-2 fits (3.4, 3.7), the aperture intensities (3.5) and
the species call (3.6) then run on the merged, pruned set. `valid_corr` of
the level-1 atoms is the one measured today against the equal-atom
simulation; for the recovered atoms the same correlation against a
simulation in which every atom has its class intensity and its own width is
written to the species table as a diagnostic. Once recovered atoms are
promoted into `_ATMS.pdb` (section 10.7), `valid_corr` moves to that
simulation for every atom, and two penalties of the present equal-atom,
equal-width simulation disappear: the one on a weak atom next to strong
neighbours, and the one on a broad atom at the surface.

## 4. Simulation

`atoms%convolve` gets a pseudo-atom branch: for a sentinel `Z` the kernel is
the Gaussian of section 1 with `B` from `beta` and amplitude
`q (4 pi / B)^(3/2)`, `q` from `occupancy`. An atom with `q` = 1 has unit
integrated intensity whatever its width.

At promotion (10.7; in the prototype `_SIM.mrc` is the present equal-atom
simulation) `_SIM.mrc` is rendered with the class intensity for `q` and the
fitted `B_i` of each atom. Against the noise-free generating density of the multi-species
cases this correlated at 0.99 or better, the same as rendering every atom
with its own fitted amplitude, so the class constraint costs nothing. One
width per class gave 0.94 to 0.99, and the equal-atom simulation used today
0.86 to 0.96 (section 8.1).

Because `occupancy` and `beta` are PDB columns, the PDB defines the
simulation. `simulate_nanoparticle pdbfile=` then is the species-free map
simulator, with per-atom amplitude and variance, and also the fixture
generator for the tests of section 10.

**B factors in the element kernels** (phase 4, ruling of 2026-10-08). In
the five-Gaussian parametrisation the scattering factor is
`f(s) = sum_j a_j exp(-b_j s^2)` and the Debye-Waller factor multiplies it by
`exp(-B s^2)`, so the product is the same sum with `b_j + B` in every term:
in real space each Gaussian of the atom has its width parameter raised by
`B`. `convolve` already does this for the resolution blur,
`b = b + bfac` with `bfac = (4 lp)^2` and `lp` at least `2 smpd`
(`simple_atoms.f90:943-945`, `:1255`, equation B.6 of Rullgard et al.), so
every element atom carries at least `(8 smpd)^2` = 8.2 A^2 at 0.358 A. The
change adds the atom's own `beta(i)` to the same sum, behind an optional
logical `bfac_pdb` of `convolve`, default false; a negative `beta` with the
flag is an error, and the pseudo-atom branch ignores the flag (its `beta`
is its width already). The integral of an atom is independent of `B`, so
the voxel sum is conserved and the peak falls as
`sum_j a_j / (b_j + bfac + B)^1.5` over `sum_j a_j / (b_j + bfac)^1.5`,
0.243 for Pt at `B` = 20 A^2 and the default blur. It is opt-in because
`simulate_atoms` and `write_centers` write coordination numbers and
`valid_corr` into the B column before rendering
(`simple_nanoparticle.f90`, the `which` cases of `write_centers_2` and the
`set_beta` calls of `simulate_atoms`); every existing call of `convolve` is
unchanged. `simulate_nanoparticle` exposes it as the key `pdb_bfac`
(`yes|no`, default `no`), accepted only with `pdbfile=`; a fixture is then
fully specified by a PDB with elements, coordinates and B factors.

## 5. When the species are given

**The maximum over per-species correlation maps.** `phase_corr` removes the
scale of the template, and at this sampling the shapes of the element
templates are nearly the same: the Pt and Ni templates correlate at 0.998. In
the emulation the correlation map made with the Pt template and the one made
with the Ni template correlated at 0.9955, and the maximum of the two found
the same atoms as either: all of the strong class and, at an intensity ratio
of 6:1, none of the weak class (section 8.5). The weak atoms are missed
because the correlation map is linear in atom amplitude and the first
threshold sits above them. A template of the weak species has the same
shape, so it does not raise them. Residual levels do (3.3).

The same holds for widths: a template of the wrong width loses little. For an
atom with twice the variance of the template the matched signal falls by 8%.

**What a species list does buy.** Three things, none of them in the
segmentation: `K` is fixed, so model selection cannot go wrong; the classes
get element names; and the species PDB (and `_ATMS.pdb` at promotion)
carries real element symbols. The class intensities and widths still come
from the data. They cannot come from the scattering-factor table, which
describes the electrostatic potential (Pt over Ni is 2.0 at the peak of the
templates and 1.65 in the integral) and not the HAADF signal.

**The species list** (phase 4, ruling of 2026-10-08: the species of a
synthesised particle are known and are given). `element=Pt,Ni`, a
comma-separated list of element symbols, brightest first. Without a comma
the key keeps today's meaning exactly, one element or one of the compound
selectors (`CdSeW`, `CdSeZ`, `CdSeR`; section 9 says how these are both a
lattice and a species statement), and every present path executes the
same branches.

- The order given is the order of decreasing expected intensity, and the
  program does not reorder it: classes found by decreasing intensity take
  the symbols in the order given. Atomic number is not a safe proxy even
  in the simulator (in the five-Gaussian model Pt integrates above Au),
  and in HAADF data the user knows which species is brighter. The symbols
  must be physical elements (the pseudo-atom symbols `X1` to `X3` are
  rejected) and distinct. The list has any length: `nspecies` is set from
  it, and everything sized by the number of classes is sized by `nspecies`
  (the mixture of 3.6, the posterior columns of `_species.csv`, the chains
  of `_species.pdb`), not by a constant. `MAX_NSPECIES` (3) remains only
  the ceiling of the automatic search of the blind mode.
- Parsing reuses the comma-list machinery of `simple_string_utils`:
  `list_of_ints2arr` is what `clustinds` goes through, and the species list
  goes through its string sibling, added beside it (blanks around an entry
  and empty entries ignored, as there). `parameters%element` becomes
  `character(len=STDLEN)` like `clustinds`, so any length fits. The element
  checks sit in a non-fatal helper, `parse_element_list(str, symbols,
  errmsg)` in `simple_atoms`, where `element_exists` lives, returning a
  logical and a message so that its rejections can be unit-tested in
  process; the fatal stop happens in the caller. `simple_parameters_phases`,
  where the element is validated today, calls it when the value holds a
  comma, stores the symbols in a new `species(:)` (`character(len=2)`,
  allocatable), sets `nspecies` to their count, a new logical
  `l_species_list` and `l_discover_species` true, and sets `element` itself
  to the first symbol, so that every present use of `params%element` sees
  one element. A list is accepted by `detect_atoms` only (`prg` is in
  `parameters`); every other program stops with a message, so that nothing
  runs silently as the list's first element.
- `exec_detect_atoms` normalises the raw command line before anything is
  checked: when `element` holds a comma it sets `discover_species=yes` on
  the command line (and rejects `discover_species=no`), so that the checks
  that follow (the half maps, `min_nbrs` and `nspecies` need discovery) see
  the implied value. `nspecies` is derived from the list, so giving it with
  a list is an error rather than something to reconcile. The combinations:

| `element` | `discover_species` | `nspecies` | `min_nbrs` | half maps | result |
| --- | --- | --- | --- | --- | --- |
| one symbol, or absent | absent or `no` | given | any | any | error: needs `discover_species=yes` |
| one symbol, or absent | absent or `no` | absent | given or half maps given | | error: needs `discover_species=yes` |
| one symbol, or absent | absent or `no` | absent | absent | absent | present path, no discovery |
| one symbol, or absent | `yes` | 0 or 1 to 3 | any | both or neither | blind discovery; `K` automatic or fixed |
| list | absent or `yes` | absent | any | both or neither | discovery with `K` = list length, classes named |
| list | `no` | | | | error |
| list | | given | | | error: `nspecies` is the list's length |
| compound selector | absent or `no` | absent | absent | absent | present path on the selector's lattice, no discovery |
| compound selector | `yes` | absent | any | both or neither | discovery with `K` = 2, classes named by the selector's symbols |
| compound selector | `yes` | given | | | error: `nspecies` is the selector's |
| list, with any program but `detect_atoms` | | | | | error |

- Level 1 then runs as `element=<first symbol>` runs today: that element's
  template and lattice table, the present pruning, the present products,
  byte-identical to the one-element run (10.5 checks it). The maximum over
  per-species maps above is why this costs nothing.
- The species call uses the list: `K` is its length, the admissibility
  test is still computed and reported (with a fixed `K` it is the
  diagnostic that says whether the data support the separation the list
  asserts; an inadmissible fit is reported in `_species.txt`, not an
  error), and `_species.pdb` carries the symbols in its element column and
  the chains `A`, `B`, `C` by class. The nanoparticle keeps `species(:)`
  for naming only; nothing in detection or in the fits reads it.
- The blind mode, `discover_species=yes` without a list, is unchanged, and
  `nspecies` applies there only. Without `element` and without the key the
  present flow runs with the Gaussian template and the `d_NN` length
  scales, and nothing more.
- The shared help text of `element` (`simple_ui_params_common.f90`) stays
  as it is, since the list form is accepted by one program; the program
  description of `detect_atoms` in `single_ui_atom.f90` documents it.

**The mixture model with B factors** is sections 3.6 and 3.7: class
intensities, per-atom class probabilities, and a width per atom estimated
inside the class model. A fuller form would compute the class probabilities
from the density likelihood of each atom under each class with `B` profiled
out, replacing the aperture sum. It uses the same model and the same code
path as stage 2; it was not emulated, and the aperture sum has the advantage
that it depends on the atom's own fit only through the small correction of
3.5.

## 6. Products

In the prototype the present products are written by the present path from
the level-1 atoms and are not changed: `_ATMS.pdb` (B column `valid_corr`,
as today), `_BIN.mrc`, `_CC.mrc`, `_MSK.mrc` and `_SIM.mrc` (equal atoms of
the reference width). Everything the prototype finds goes into four new
files, so that a run with `discover_species=yes` differs from one without it
only by their presence:

- `_species.csv`, one row per atom at full precision, level-1 and recovered
  atoms alike: position, radius, coordination number, stage and level that
  accepted it, detection z, amplitude, `B` of stage 1 and of stage 2,
  aperture intensity and its signal-to-noise, class, class probabilities,
  `valid_corr`, and whether the gate was needed to accept it.
- `_species.pdb`, the same atoms for a viewer: element column and chain by
  class (`X1` and chain `A` the strongest class, `X2` and `B` the next;
  with a species list, the symbols given, in their order),
  occupancy the class intensity over the strongest class's, B column the
  stage-2 `B_i`. It is the layout promotion gives `_ATMS.pdb`, written to a
  file of its own so that `_ATMS.pdb` stays as it is.
- `_species_radial.csv`: the profiles of 3.8.
- `_species.txt`: `K`; class intensities, ratios and fractions; `D`; the BIC
  table with admissibility; the expected misclassification rate; `d_NN`,
  `B_ref`, the noise numbers and their source (half maps or region); the
  calibrated thresholds and the expected false count at each; atoms added
  per stage and level; the diagnostics of section 7.

Promotion, after the validation of section 10.7, adds the recovered atoms
and the class information to the present products:

- `_ATMS.pdb`. Element column: the class symbol (`X1`, `X2`, `X3`, registered
  in `get_element_Z_and_radius` with sentinel `Z`, precedent `CDSE`, 999; or
  the element names when a species list is given). Chain: `A`, `B`, `C` by
  class, so a class is one selection in a viewer. Occupancy: `I_k / I_1`.
  B column: `B_i`. Today the B column holds `valid_corr`
  (`simple_nanoparticle.f90:1821`); within SIMPLE only the unit tester reads
  it back, and `valid_corr_in_bfac_field.pdb` keeps it for `atoms_stats`.
- `_SIM.mrc` as in section 4.
- `_CC.mrc`, `_BIN.mrc`, `_MSK.mrc`. Atoms found at residual levels have no
  connected component of the level-1 threshold. Their components are derived
  from the data, never from the model: the connected region of the residual
  z map above `k_B / 2` around the accepted maximum, split between
  neighbouring maxima by a watershed on the z map (the merging of components
  at a low threshold is what moved detection to local maxima, so the split
  is required, not optional). The first version of this note generated the
  component from the fitted Gaussian; that would make the size test of
  `discard_atoms` circular and change what `atoms_stats` measures, and is
  withdrawn. `atoms_stats` reads `_CC.mrc` with `_ATMS.pdb`, so promoted
  atoms need a component. Not emulated.

## 7. Diagnostics and limits

- **Core and surface.** For an atom of integrated intensity `I` and width
  sigma the best possible detection signal-to-noise scales as
  `I sigma^(-3/2)`. An atom with twice the variance is detected at 0.59 of the
  signal-to-noise of the same atom in the core. That loss belongs to the data
  and no method avoids it; fitting the atom needs it to be found first. Using
  a core-width template on it costs a further 8%. A bank of two template
  widths was tried and dropped: it changed the number of atoms found by a few
  in either direction and admitted more noise peaks. So detection is less
  complete at the surface, for the light species first.
  The design does three things about it. The species call of an atom that
  was found does not depend on detection sensitivity, and neither does its
  fitted width. The predicted detection signal-to-noise is reported per class
  and radial shell, computed from the class intensity and the fitted width of
  that shell. And the shells where it falls below about 6.5 are flagged:
  there the emulation found 56-81% of the light atoms against 88-100% above
  (section 8.2), and the atoms that were found were the narrower ones, so the
  shell average of the light class is biased low in width.
- **Detection limit overall.** In the emulation a class whose found atoms had
  a median fitted amplitude of 6.6 or more times `s_A` was found to 88-100%;
  at 6.4 it was 84%, at 5.0 it was 65% and at 4.1 it was 45%. The same number
  can be read off an existing map before any code is written: the peak of the
  correlation map at the weakest visible atoms over its standard deviation
  away from the particle.
- **Missing species.** When a class is below the limit the program reports
  one class and cannot know better from intensities. The coordination of
  interior atoms can: with a quarter of the sites undetected the inner 30% of
  the atoms had a mean coordination of 8.9 where 12 is expected (case 8b).
  Reported as the interior coordination deficit.
- **What a class is.** An intensity class. A site that is occupied part of
  the time, or an atom that moved during acquisition, integrates to less than
  its species and can form a class of its own. The per-class coordination and
  radial profile are reported so that a class confined to the surface is
  visible. It is a flag and not a verdict, since surface segregation of the
  light species gives the same picture.
- **Marginal atoms.** Stage B occasionally admits a noise peak at a plausible
  surface site (section 8.6 gives the count). Each atom carries its detection
  z, its stage and level, and whether the gate was needed to accept it.
- **Expected false atoms.** The calibration of 3.2 gives, for each stage,
  the expected number of noise maxima above its threshold in its search
  volume. Both numbers are reported with the thresholds and the noise source
  (half-map difference or region), together with the ratio of the region's
  robust to plain spread, which flags a contaminated region.
- **Half-map agreement.** With `vol_even` and `vol_odd` given, the aperture
  intensities are measured at the final positions in each half map, and the
  label agreement and the intensity correlation are reported. The halves
  share their alignment, so this measures the effect of noise on the labels
  and not model bias. Not emulated.
- **Model bias.** Labels are found from scratch at every call; nothing is
  carried between iterations. Since `_SIM.mrc` only aligns (section 0), a
  wrong label can reach the next map only through the orientations. A direct
  test is cheap and worth having once `_SIM.mrc` carries the labels (after
  promotion, 10.7): render a random 10% of the atoms at the pooled intensity
  and check that they separate as well as the rest in the next iteration.
  Not emulated.

## 8. Emulation

Particle: FCC, a = 3.84 A, diameter 24 A, 528 atoms, box 128 at 0.358 A.
Gaussian atoms, sigma 0.35 A in the core with 5% scatter, integrated
intensity 1 for the strongest class. Two models of broadening, both at
conserved integral: "surface x2" doubles the variance of the 195 atoms with
coordination 9 or less; "radial x2" (or x3) lets the variance grow with the
square of the radius to twice (three times) its central value at the surface.
Noise is given as its standard deviation per voxel over the peak of a core
atom of the strongest class. One particle per case unless stated.

### 8.1 Detection, discovery and simulation

"Level 1" is the present flow. Counts are atoms found per class, strongest
first. "SNR" is the median fitted amplitude over `s_A` per class. The last
column is the correlation with the noise-free generating density of the
equal-atom simulation used today and of the proposed one.

| # | Case | Level 1 | All levels | False | K | Class intensities | SNR | Simulation: today, proposed |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3:1 atoms, intensity 3:1, noise 0.03 | 396, 0 | 396/396, 132/132 | 1 | 2 | 1 : 0.330 | 79, 27 | 0.937, 0.998 |
| 2 | 3:1 atoms, intensity 6:1, noise 0.03 | 396, 0 | 396/396, 132/132 | 0 | 2 | 1 : 0.167 | 76, 12 | 0.903, 0.999 |
| 3 | 6:1 atoms, intensity 3:1, noise 0.05 | 453, 0 | 453/453, 75/75 | 0 | 2 | 1 : 0.336 | 46, 15 | 0.960, 0.997 |
| 4 | 6:1 atoms, intensity 6:1, noise 0.05 | 453, 0 | 453/453, 73/75 | 0 | 2 | 1 : 0.161 | 46, 7.8 | 0.943, 0.998 |
| 5 | as 2, noise 0.05, surface x2 | 396, 0 | 396/396, 116/132 | 0 | 2 | 1 : 0.169 | 42, 7.0 | 0.871, 0.997 |
| 6 | as 5, light species at the surface | 396, 0 | 396/396, 60/132 | 0 | 2 | 1 : 0.155 | 47, 4.1 | 0.918, 0.997 |
| 7 | as 5, light species in the core | 396, 0 | 396/396, 130/132 | 0 | 2 | 1 : 0.159 | 37, 8.6 | 0.862, 0.997 |
| 8 | as 2, noise 0.05 band-limited at 1 A | 396, 0 | 396/396, 132/132 | 0 | 2 | 1 : 0.169 | 49, 7.9 | 0.901, 0.998 |
| 8b | as 2, strongly correlated noise 0.05 | 395, 0 | 396/396, 1/132 | 0 | 1 | - | 12, 3.1 | 0.981, 0.984 |
| 9 | single species, noise 0.05 | 528 | 528/528 | 0 | 1 | - | 41 | 0.994, 0.998 |
| 10 | single species, noise 0.05, surface x2 | 528 | 528/528 | 0 | 1 | - | 45 | 0.952, 0.998 |
| 11 | as 10, noise band-limited at 1 A | 528 | 528/528 | 0 | 1 | - | 41 | 0.948, 0.998 |
| 12 | three species 2:1:1, intensities 1 : 0.5 : 0.167, noise 0.03 | 264, 132, 0 | all 528 | 0 | 3 | 1 : 0.495 : 0.166 | 79, 40, 13 | 0.877, 0.999 |
| R1 | as 2, noise 0.05, radial x2 | 396, 0 | 396/396, 111/132 | 0 | 2 | 1 : 0.173 | 33, 6.4 | 0.896, 0.996 |
| R2 | as R1, light species at the surface | 396, 0 | 396/396, 86/132 | 0 | 2 | 1 : 0.166 | 34, 5.0 | 0.909, 0.996 |
| R3 | single species, noise 0.05, radial x2 | 528 | 528/528 | 0 | 1 | - | 32 | 0.976, 0.997 |
| R4 | 3:1 atoms, intensity 3:1, noise 0.05, radial x3 | 394, 0 | 396/396, 128/132 | 1 | 2 | 1 : 0.314 | 25, 8.5 | 0.895, 0.995 |
| R5 | 6:1 atoms, intensity 6:1, noise 0.03, radial x2 | 453, 0 | 453/453, 74/75 | 0 | 2 | 1 : 0.170 | 52, 8.4 | 0.928, 0.999 |
| R6 | as R1 at noise 0.03, light species 1.3x the variance of the heavy one at every radius | 396, 0 | 396/396, 129/132 | 0 | 2 | 1 : 0.163 | 52, 6.6 | 0.881, 0.998 |
| R7 | single species, noise 0.05, radial x3 | 528 | 528/528 | 0 | 1 | - | 26 | 0.954, 0.996 |

Every atom that was found and matched to the generating model received the
right label, in every case. `D` was 10.1 or more wherever more than one class
was found. Class fractions equalled the fractions among the found atoms;
where weak atoms were missed (cases 5, 6, R1, R2) the light fraction is
underestimated accordingly (0.227, 0.132, 0.219 and 0.178 against 0.25).

In the multi-species cases the simulation with one width per class
correlated at 0.941-0.994 and the one with per-atom amplitudes at
0.995-0.999, against 0.995-0.999 for the proposed one.

### 8.2 Radial profiles

Case R6, where the light species is broader than the heavy one at every
radius. Widths are root-mean-square sigma per shell in Angstroms: fitted with
its standard error, then the generating value. "Stage 1" is the width before
the class tie (3.4).

| Class | Atoms | Mean radius (A) | Stage 1 | Stage 2 | Generating | Intensity / I_1 | Predicted SNR | Found in shell |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 80 | 5.3 | 0.383 | 0.384 +- 0.002 | 0.384 | 0.996 | 66 | 100% |
| 1 | 79 | 8.0 | 0.417 | 0.417 +- 0.003 | 0.417 | 1.002 | 59 | 100% |
| 1 | 79 | 9.6 | 0.448 | 0.448 +- 0.002 | 0.448 | 1.003 | 53 | 100% |
| 1 | 79 | 10.8 | 0.473 | 0.475 +- 0.003 | 0.477 | 1.001 | 48 | 100% |
| 1 | 79 | 11.7 | 0.489 | 0.492 +- 0.003 | 0.492 | 0.998 | 46 | 100% |
| 2 | 26 | 5.7 | 0.416 | 0.434 +- 0.007 | 0.438 | 0.174 | 9.0 | 100% |
| 2 | 26 | 8.4 | 0.447 | 0.475 +- 0.010 | 0.477 | 0.162 | 7.8 | 100% |
| 2 | 26 | 9.7 | 0.466 | 0.498 +- 0.007 | 0.516 | 0.161 | 7.3 | 100% |
| 2 | 26 | 10.7 | 0.466 | 0.504 +- 0.008 | 0.532 | 0.157 | 7.1 | 96% |
| 2 | 25 | 11.7 | 0.479 | 0.528 +- 0.008 | 0.548 | 0.161 | 6.6 | 92% |

What the profiles of all radial cases show:

- **Strongest class.** The fitted shell width was within 0.006 A (1.1%) of
  the generating value in every shell of every case, including the single
  species cases R3 and R7 and the factor-three gradient of R4 and R7.
- **Light class.** With stage 2 its own profile is recovered: in R6 it is
  broader than the heavy class in every shell, by 6% to 14% in sigma at
  similar radii where 14% was generated. Stage 1 alone pulls it to the heavy
  profile (first width column). Over all radial cases the fitted width of the
  light class was within 3% of the generating value in the inner shells, up
  to 6% low in outer shells where detection was complete or nearly so, and up
  to 10% low where it was not (R2, 56-70% found): there the atoms that are
  found are the narrower ones.
- **Intensity.** The class intensity is flat in radius for the strongest
  class (0.987 to 1.008). For the light class it is flat within its errors in
  the inner shells and falls by up to 10% in the outermost ones (0.157-0.161
  against 0.174 here).
- **Detection.** Over the light-class shells of R1, R2, R4, R5 and R6, a
  predicted signal-to-noise of 6.6 or more went with 88-100% of the atoms
  found, and 6.1 or less with 56-81%.

### 8.3 What the cases show

- Cases 1-4: both readings of "3-6:1", and both together, are within reach at
  these noise levels.
- Cases 5-7 and R1-R6: widths that grow towards the surface do not confuse
  the labels, wherever the light species sits, and the widths are recovered
  per class. They do limit detection: a light atom that is also broad is the
  hardest object (cases 6 and R2).
- Case 8: the thresholds hold when the noise is not white.
- Case 8b: when the light species is below the detection limit the answer is
  one class, with a coordination deficit. The aperture statistic is also at
  its worst here, because noise at low frequency adds coherently over the
  sphere.
- Cases 9-11, R3 and R7: no species is invented from broadening or from
  noise.
- Case 12: `K` is not tied to two.

### 8.4 Choice of the classification statistic

Three statistics were run through the same class discovery.

| Cases | Amplitude at fixed width | Integral from a single-atom fit with a width prior | Aperture intensity (3.5) |
| --- | --- | --- | --- |
| 1-4, 8, 12, R2 | right K, all labels right | right K, all labels right | right K, all labels right |
| R1, R5, R6 | right K, 98.9-99.2% | right K, 99.8-100% | right K, 100% |
| 5 | K = 3 | K = 2, 99.8% | K = 2, 100% |
| 6 | K = 3 | K = 2, 85.7% | K = 2, 100% |
| 7 | K = 3 | K = 1 | K = 2, 100% |
| R4 | K = 1: the light species is not seen | K = 1 | K = 2, 100% |
| 10, single species | K = 2: 37% of the atoms called a second species | K = 2 | K = 1 |
| 11, single species | K = 2 | K = 1 | K = 1 |
| 9, R3, R7, single species | K = 1 | K = 1 | K = 1 |
| 8b | K = 1: the light species was not detected | K = 1 | K = 1 |

The amplitude at a fixed width has the lowest noise and is the wrong
statistic: it reads a broad heavy atom as a lighter species. An earlier
review of the first draft (2026-10-07, before the one recorded in section
15) recommended it on the evidence of particles with uniform widths; this
table withdraws that recommendation.

### 8.5 Per-species correlation maps

The same synthetic maps, correlated with the Pt template, with the Ni
template, and the voxel-wise maximum of the two. "Gate" is the Otsu threshold
of the present flow.

| Case | Correlation of the two maps | Weak peak / strong peak | Gate / strong peak | Weak atoms above the gate: Pt map, Ni map, maximum |
| --- | --- | --- | --- | --- |
| 1 (intensity 3:1) | 0.9955 | 0.34 | 0.25 | 132, 132, 132 of 132 |
| 2 (intensity 6:1) | 0.9957 | 0.17 | 0.25 | 0, 0, 0 of 132 |
| 4 (intensity 6:1) | 0.9955 | 0.17 | 0.25 | 0, 0, 0 of 75 |

The three maps behave as one. At 3:1 the weak atoms clear the gate in all of
them, and the threshold criterion of the present flow then drops them: level 1
of case 1 returned none.

### 8.6 Repeats

Five seeds each for cases 2, 4, 6, 9, 10, R1, R2 and R3 (40 runs): the right
`K` in 40 of 40 and every label right. In the 15 single-species runs `K` was 1
throughout and one noise peak was accepted in one run. Recall of the light
species ranged over 73-75 of 75 (case 4), 60-74 of 132 (case 6), 94-107 of 132
(R1) and 75-89 of 132 (R2), and its intensity estimate over 0.162-0.187 where
0.167 was generated.

False atoms over all 60 runs of sections 8.1 and 8.6: five, one each in five
runs. The one that was examined (case 1) was accepted at stage B, at a surface
site with four accepted neighbours.

### 8.7 Not covered

`split_atoms`, `discard_atoms`, `validate_atoms`, the soft mask and the
iterated-conditional-modes denoising that a map may have been through
(`icm_imgfile`) were not emulated. Atoms are exact Gaussians, the background is zero,
positions are not refined, the particle is a sphere, and the noise is
stationary. The noise levels are guesses; what transfers is the
signal-to-noise scale of section 7. `D` = 3, the class-size floor and the 5%
floor on `tau_k` are starting values; the detection thresholds are
calibrated per map (3.2), and the target of 0.2 false atoms per particle is
a starting value.

### 8.8 Checks of 2026-10-08

A second implementation of the generator of section 8 and of the detector of
3.3 (`codex_checks.py`, numpy only, written without reference to
`emul_species.py`) was run to test the review's points about the thresholds
and the gate. Same box, sampling, template, lattice and diameter; its own cut
of the sphere gave 531 atoms, 210 of them with coordination 9 or less (528
and 195 in the emulation of 8.1); the level-1 set was taken as the strong
class at its true positions jittered by 0.05 A and fitted on the map, as the
emulation found it in every case; stage-1 fits with positions fixed and no
background; six noise maps per colour, five seeds per case.

Noise maxima above `k` in a pure-noise z map, standardised with the region's
mean and standard deviation, per map (search over the shipped sphere of
radius 21 A, 845,000 voxels; over the particle envelope of radius 14.7 A,
291,000 voxels; and inside the noise region, 233,000 voxels):

| Noise | k = 3.5: sphere, envelope, region | k = 4 | k = 4.5 | k = 5 | Gaussian tail x sphere voxels at k = 4 |
| --- | --- | --- | --- | --- | --- |
| white | 29.0, 10.3, 7.2 | 4.2, 1.5, 1.2 | 0.33, 0, 0 | 0, 0, 0 | 26.8 |
| band-limited at 1 A | 27.3, 10.2, 6.7 | 4.2, 1.5, 1.3 | 0.67, 0.17, 0 | 0, 0, 0 | 26.8 |
| strongly correlated (2.5 A) | 11.3, 4.2, 2.3 | 2.3, 0.67, 0.83 | 0.33, 0, 0.33 | 0, 0, 0 | 26.8 |

The region's count scaled by volume predicted the sphere's count at `k` = 4
to within a factor 0.98, 0.86 and 0.77. The calibration of 3.2 gave `k` =
5.01, 5.00 and 4.73 for a target of 0.2 false atoms per particle, with 0,
0 and 0.17 observed at those thresholds. For the tester of phase 1: the
white-noise region counts above 2.5, 3.0 and 3.5 were 154.5, 38.0 and 7.2
per map, the search-to-region volume ratio 3.628, the fitted `C` 2105
(the Euler-characteristic value for white noise smoothed by the template,
`V / FWHM^3 (4 ln 2)^(3/2) / (2 pi)^2`, is 1656), and the resulting `k`
4.99. `detect_peak_thres_fdr` of
`simple_segmentation.f90` (the false-discovery-rate routine: median and
robust spread of the peak values as the null, Benjamini-Hochberg selection
bounding the expected fraction of false discoveries at q = 0.05), applied to
all local maxima of the sphere, selected 0.2 to 0.3 noise peaks per map: the distribution of maxima
(mean z 1.84, spread 0.73) has a heavier right tail than its median and
spread describe.

Recall of the weak atoms and false atoms per particle at noise 0.05 with
surface broadening, three cases: the surface shell of case 6 (132 weak atoms
at intensity 1/6, coordination 5 to 9); a weak-only domain (118 atoms with
x > 4 A at intensity 1/6, 57 of them interior); 40 weak adatoms at
intensity 1/6 on sites of coordination 2 to 4 outside a single-species
particle. "Calibrated" is the threshold of 3.2 (5.0 within 0.02 on these
maps).

| Variant | Shell: recall, false | Domain: recall, false | Adatoms: recall, false | Domain interior (coordination 10-12) | Domain surface (5-7) |
| --- | --- | --- | --- | --- | --- |
| stage A only, k = 5 | 25.5%, 0.2 | 58.6%, 0 | 24.5%, 0 | 93% | 26% |
| A at k = 5, B gated at k = 4 | 60.6%, 0.4 | 77.6%, 0.4 | 57.5%, 0.2 | 99% | 55% |
| A at k = 5, B ungated at k = 4 | 66.2%, 3.8 | 82.4%, 3.2 | 68.0%, 3.8 | 100% | 64% |
| A at k = 5, B ungated at k = 4.5 | 45.9%, 0.6 | 71.0%, 0.6 | 43.0%, 0.4 | 99% | 43% |
| calibrated A, gated B at k = 4 | 60.6%, 0.4 | 77.6%, 0.4 | 57.0%, 0.2 | 99% | 55% |
| `detect_peak_thres_fdr`, q = 0.05, ungated | 55.5%, 1.6 | 76.8%, 1.4 | 52.5%, 0.4 | 99% | 54% |
| A at k = 5, B gated at k = 3.5 | 72.3%, 0.8 | 84.9%, 0.8 | 70.0%, 0.2 | 99% | 70% |

The ungated search at `k` = 4.5 is the one with a false count close to the
gated stage B's; the gate wins against it everywhere. The FDR routine
controls what it promises (its false discoveries are 2% of its discoveries
here) but admits three to four times the gate's false atoms at lower recall.
The gated stage B at `k` = 3.5 shows the gate can be pushed; its false count
doubles and it was not adopted.

Contamination of the region by an undetected weak-only domain, with the
level-1 atoms excluded by 1.5 `d_NN`: standard deviation measured over true,
1.054 at a mask radius of 21 A (337,000 voxels, 29 undetected atoms inside
the region), 1.356 at 16 A (30,500 voxels), 1.398 at 14 A (12,100 voxels);
robust spread 1.019, 1.221, 1.298. With the surface shell undetected
instead: 1.000 to 1.002, no weak atom inside the region.

What these checks do not change: they use the same idealised generator as
section 8, so they say nothing about reconstructions, backgrounds, the
five-Gaussian kernels or position refinement. That is what section 10.7 is
for.

## 9. Code map

The prototype touches `detect_atoms` and the modules under it, and nothing
else that runs today. Every new behaviour sits behind a condition that is
false for a run as it is given today (`element` present, `discover_species`
absent), so the present path is the same code executing the same branches.

New: `src/main/nano/simple_nano_species.f90`, array numerics with no image
dependence: the mixture fit, BIC, admissibility and the choice of `K`
(`fit_species_mixture(x, s_meas, nspecies_in, K, labels, post, mu, var)`; the
routine is not named after the key; the fit itself is the extreme
deconvolution `xd_gmm` of `src/utils/clustering/`, driven in one dimension
from the two starts of 3.6, so the module owns no EM of its own), the
enclosed-fraction function `F` of 3.5, and the threshold calibration of 3.2
from an array of region maxima counts. Tester
`simple_nano_species_tester.f90`, sub-suite `species` of `unit_single`
(table `suites_single` in `simple_commanders_test_class.f90`, list in
`simple_test_ui_class.f90`, checked by `scripts/check_test_registry.py`).

`src/defs/simple_defs_atoms.f90`: the pseudo-atom symbols `X1`, `X2`, `X3`
in `get_element_Z_and_radius`, each with its own sentinel `Z` in a reserved
range held in named constants (precedent `CDSE`, 999; the range must not
collide with it), radius 1.0 until the data give `0.4 d_NN`.

`src/main/nano/simple_atoms.f90`: `convolve` gets the pseudo-atom branch of
section 4 for `Z` in the reserved range, reading `beta` as `B` and
`occupancy` as `q`; `element_exists`, `set_element`, `guess_an_element` and
`Z_and_radius_from_name` accept the symbols. The per-element table, `epot`
and every existing branch stay as they are. `simple_atoms_tester.f90`
gains the closed-form tests of phase 1.

`src/main/nano/simple_nanoparticle_utils.f90`: `phasecorr_one_atom` gets an
optional kernel width (`B_ref`) and, when it is given, correlates with the
pseudo-atom instead of the element; `est_nn_dist(centers)`; `find_rMax` and
`calc_contact_scores` get an optional `d_NN` that replaces the element
lookup. Callers that do not pass the new arguments behave as before.

`src/main/nano/simple_nanoparticle.f90`:

- a logical `l_species_free`, set in `new` when `params%element` is empty.
  `new` then does no `Z` lookup; `theoretical_radius`, the neighbour cutoff
  and the contact-score ceiling are set from `d_NN` once the first
  binarisation has produced centres (section 2), in `identify_atomic_pos`,
  before `discard_small_ccs`. With an element, `new` is unchanged.
- `atom_stats`: `amp`, `bfac`, `bfac_stage1`, `aper_int`, `aper_snr`,
  `det_z`, `det_stage`, `det_level`, `gate_used`, `species`, `species_post`.
  `cendist`, `u_iso` and `isocorr` exist.
- `identify_atomic_pos`: the present sequence, with the pseudo-atom template
  and the `d_NN` length scales substituted at the element entry points when
  `l_species_free` (`phasecorr_one_atom` at `:587`, `t2c` through
  `write_centers` and `atoms%convolve`, `calc_contact_scores` at `:1045`,
  `:1121`, `:1156`, `get_lattice_params` at `:1032-1033`, `simulate_atoms`
  at `:1775-1776`). After `discard_atoms`, and only when
  `params%l_discover_species`: `fit_atoms_joint` on the level-1 atoms
  (stage 1, the provisional fit that the subtraction of 3.3 needs),
  `est_noise_region`, `calibrate_threshold`, `detect_residual_levels`,
  `prune_recovered`, `fit_atoms_joint` (both stages, merged set),
  `calc_aperture_int`, `assign_species`, a diagnostic render of the
  class-intensity, per-width simulation into a scratch image through the
  pseudo-atom `convolve` (for the recovered atoms' `valid_corr` column; not
  written to disk), `write_species_report`, `write_radial_profiles`; then
  the present output block, unchanged, from the level-1 atoms. The
  recovered atoms live in a second atom table (`recovered(:)` of the same
  type) so that nothing that writes the present products can see them until
  promotion. `l_species_free` belongs to the nanoparticle (set in `new`);
  `l_discover_species` and `min_nbrs` belong to `parameters`, derived from
  or validated with the keys in `simple_parameters_phases.f90`.
- `simulate_atoms`, `write_centers`, `kill`: unchanged in the prototype,
  except that `kill` frees the new table.

`calc_isotropic_disp` stays for `atoms_stats`. Its log-linear estimator over
positive voxels is biased in noise (at noise 0.03 it returned amplitude 0.91
and sigma 0.371 A for true 1 and 0.350); `atoms_stats` can take `u_iso` from
the stage-2 fit once that exists.

Parameters and UI, prototype: `element` becomes optional
(`required_override=.false.`) in `detect_atoms` only
(`single_ui_atom.f90:226-227`); three new keys on `detect_atoms`:
`discover_species` (`yes|no`, default `no`: residual recovery and the species
call, diagnostics only), `nspecies` (integer, default 0 = automatic `K`),
`min_nbrs` (integer, default 3, 0 = no neighbour requirement); `vol_even`
and `vol_odd` (existing
keys) accepted as the optional half maps of 3.2. The thresholds other than
the calibrated ones stay module constants. `exec_detect_atoms`
(`simple_commanders_atoms.f90:405-423`) validates the combination before
`params%new`: `nspecies` and `min_nbrs` need `discover_species=yes`; both
half maps or neither. `atoms_stats`, `autorefine3D_nano`,
`conv_atom_denoise` and `analysis2D_nano` keep requiring `element` until
promotion.

Later, with promotion (section 10.7): `element` optional in the other four
programs; `simulate_atoms` sets `occupancy` and `beta` per atom and
`write_centers` writes symbol, chain, occupancy and `B`; data-derived
components for the recovered atoms (section 6); `autorefine3D_nano` passes
the half maps to `detect_atoms` and copies the species files with the
other per-iteration products; `analysis2D_nano` without an element cannot
build the ideal-lattice start (Q4).

Phase 4 adds, in the same areas: `parse_element_list` in `simple_atoms.f90`
with its tests in the atoms tester; the `bfac_pdb` argument of `convolve`
(section 4); `pdb_bfac` through the whole parameter lifecycle, declaration
and default in `simple_parameters.f90` with its logical, registration in
`simple_parameters_parse.f90`, `yes|no` validation and the `pdbfile`
dependence in `simple_parameters_phases.f90`, the key on
`simulate_nanoparticle` in `single_ui_atom.f90`, and the flag passed to
`convolve` in `exec_simulate_nanoparticle` (`simple_commanders_sim.f90`);
`species(:)`, `l_species_list` and the list parsing in the parameters
(section 5); the list rules in `exec_detect_atoms`; `species(:)` copied
into the nanoparticle in `new` for `write_species_pdb`.

**The compound selectors `CdSeR`, `CdSeW`, `CdSeZ`** (phase 4; the
maintainer asked on 2026-10-08 that the present species mixtures be handled
correctly). Today `nanoparticle%element` is `character(len=2)`, so the
selector is cut to `Cd` before every lattice lookup, which then falls back
to fcc with `a` = 3.76 (`simple_nanoparticle.f90:216`, `:300`, `:1115`,
`:1804`, `:2459`; the table is `get_lattice_params` in
`simple_defs_atoms.f90:223-271`, which knows `CDSER`, `CDSEW` and `CDSEZ`),
and the utilities `fit_lattice`, `run_cn_analysis`, `strain_analysis`,
`calc_contact_scores` and `find_rMax` take the element as
`character(len=2)` and cut it the same way. Widening those arguments alone
is not a repair: `find_rMax` uses the one value for the lattice table,
where `CDSEW` is valid, and for `get_element_Z_and_radius`, where it is not
and stops (`simple_nanoparticle_utils.f90:271` and following), and
`strain_analysis` uses it for the lattice and for PDB names. The repair
keeps two values apart everywhere: the selector as given (`element_key`,
five characters, for every lattice lookup) and the two-letter symbol of
the first species (for radii, atom names and `set_element`). The utilities
take the selector as `character(len=*)` and derive the symbol as its first
two characters themselves, as `simulate_nanoparticle` does (`el1`,
`simple_commanders_sim.f90:439`); the nanoparticle holds both. For a
one-element selector the two values coincide and nothing changes, so the
Pt regression stays byte-identical; for the compound selectors the lattice
becomes the one the table defines, a change on a path that was wrong.
Regressions: `find_rMax('CdSeW')` and `get_lattice_params('CDSEW')` in a
unit test; and `single_atoms_stats` run a second time by CTest with
`element=CdSeW` (the program already takes `element`), with the Pt floors.
If its assertions prove specific to the fcc coordination, they are
parameterised by the crystal system, not skipped.

A compound selector is also a species statement: `CdSeW` says the particle
holds cadmium and selenium. The parser of section 5 therefore returns, for
a selector, the symbols it names (`Cd`, `Se`, in the order written, which
for the three selectors is the brighter first) together with the selector
itself; without `discover_species=yes` nothing more happens than today
(now on the right lattice), and with it the species call runs with `K` = 2
and those names, exactly as `element=Cd,Se` would, but with the wurtzite,
zincblende or rocksalt lattice of the selector instead of the close-packed
cutoff of the first symbol. A comma list and a compound selector together
(`CdSeW,Pt`) are rejected.

## 10. Phases and the tests that gate them

Sections 10.1 to 10.4 are phases 0 to 3, the first run (done 2026-10-08);
10.5 and 10.6 are phases 4 and 5, the second run; 10.7 is what comes after
them.

### 10.0 Rules of coexistence

The present detection is working production code with a regression test,
and it stays that way throughout. Concretely:

- No existing routine changes its behaviour for a call as it is made today.
  New behaviour is reached only through `l_species_free` (no `element`) or
  `l_discover_species` (`discover_species=yes`), both false for every
  existing caller and every existing test. Optional arguments added to
  utilities default to the present behaviour.
- `single_atoms_stats` is run on every phase's build and must pass with its
  present floors; its atom count, recall, precision, root-mean-square error
  and correlation are recorded in Progress each time and compared with
  phase 0. Any drift is a stop condition, not something to explain.
- The present products are written by the present code from the level-1
  atoms; the prototype writes only its own files (section 6).
- No new CMake target, no per-file compile option, no new global flag
  (`doc/policies/compile_time_policy.md`); a new tester goes into the
  existing `unit_single` executable through the suite table; a new CTest
  entry raises `SIMPLE_CTEST_BUDGET` in `production/CMakeLists.txt` by one.
- Release 4 rules: nothing is aliased or kept for compatibility; the only
  on-disk format that must keep reading is the `.simple` project.
- The fast gate (`scripts/run_fast_gate.sh`: the test-registry, description
  and `scripts/check_flex_dag.py` checks, then the `fast` and `provisional`
  CTest labels) passes at the end of every phase, in Debug.

### 10.1 Phase 0: baseline

Build the copy as it was handed over (Debug with `BUILD_TESTS=ON`, the
machine's compiler and environment recorded), run the fast gate and
`single_atoms_stats`, and record its numbers. Then run the same simulation
and detection by hand (`simulate_nanoparticle element=Pt moldiam=20
box=160 smpd=0.358`, `detect_atoms element=Pt`) and keep the input volume,
`simatms.pdb` and `_ATMS.pdb` under the run's `scratch/keep/`: they are the
reference of phase 2's equivalence test. Record the atom count of
`simatms.pdb` and of `_ATMS.pdb`, the positions' root-mean-square
difference and the `_SIM.mrc` correlation. Exit: fast gate passes;
`single_atoms_stats` passes; the reference is kept and its numbers are in
Progress.

### 10.2 Phase 1: foundations that change nothing

The pseudo-atom kernel in `convolve` with its symbols and the PDB round
trip through `occupancy` and `beta`, so that `simulate_nanoparticle
pdbfile=` renders pseudo-atoms; and `simple_nano_species` with its tester.
Exit: unit tests on one and two pseudo-atoms against the closed forms (peak
`q (4 pi / B)^(3/2)`, voxel sum `q / smpd^3` within 1% at a cutoff of
`4 sigma`, two atoms adding linearly), and a round trip PDB to `atoms` to PDB
preserving symbol, occupancy and `B`; `discover_species` returns the right
`K`, labels and BIC ordering on fixed samples (one class; two classes at
`D` of 2.5 and 4, the first inadmissible; a class below the size floor;
three classes; `fit_species_mixture`) and the calibration returns 5.0
within 0.1 from the region counts of section 8.8 (white noise) scaled to
its geometry; the `species`
sub-suite is registered and `scripts/check_test_registry.py` is clean; fast
gate; `single_atoms_stats` unchanged against phase 0.

### 10.3 Phase 2: level 1 without an element

`d_NN`, the derived length scales, the Gaussian template, `element`
optional in `detect_atoms`. Exit: `detect_atoms` on the phase 0 input
volume without `element` gives the same atom count as the reference within
1%, positions matched within 0.5 A with a root-mean-square difference of
at most 0.1 A, and the same `_SIM.mrc` correlation within 0.01; the same
command with `element=Pt` gives a `_ATMS.pdb` identical to the reference
(same count, positions within 1e-3 A: the element path did not move); the
measured `d_NN` is within 2% of 2.766 A (Pt, `a` = 3.912); fast gate;
`single_atoms_stats` unchanged. If the species-free count differs by more
than 1% the phase stops and reports which step (template, threshold
ladder, splitting, pruning) lost or gained the atoms, with the counts after
each.

### 10.4 Phase 3: residual recovery and species call, diagnostics only

Noise region and half-map reference, calibrated thresholds, the two stages
with the gate switch, pruning of the recovered atoms on the atom table,
both fit stages, aperture intensities, species assignment, radial profiles,
the files of section 6. The fixtures are rendered by
`simulate_nanoparticle pdbfile=` from pseudo-atom PDBs written by the test
itself (the generating models of cases 2, 9, R3 and R6 of section 8 on the
Pt lattice of phase 0's particle, so that level 1 is known to work on them),
with Gaussian noise added in process at a fixed seed and at the section 8
noise levels relative to the peak of a core atom of the strongest class.
Exit, as a new high-level test `species_discovery` (registered like
`single_atoms_stats`, budget raised by one): for the two-species fixtures,
recall of the weak class at least 0.90 where its predicted detection
signal-to-noise (section 7) is 6.5 or more, no more than one false atom per
particle, `K` = 2, every label of a found atom right, the intensity ratio
within 10%, the per-shell width of the strongest class within 3% and of the
light class within 10%; for the single-species fixtures, `K` = 1 and no
recovered atom; for every fixture, `_ATMS.pdb`, `_BIN.mrc`, `_CC.mrc`,
`_MSK.mrc` and `_SIM.mrc` identical to the run of the same command without
`discover_species` (the products did not move); fast gate;
`single_atoms_stats` unchanged. The floors are those of the emulation and
are set before the first run; a floor the code does not reach stops the
phase with the measured value, it is not lowered.

### 10.5 Phase 4: B factors in the element kernels and the species list

The second run starts from the maintainer's checkout after the first run
and its follow-ups (section 16 and the report list them). Two rulings of
2026-10-08 define it: B factors go into the element kernels now, as an
opt-in (section 4), and the species of a synthesised particle are given on
the command line as a list (section 5). The reason is the review's main
objection, model-matched evidence: the test of phase 3 renders pseudo-atoms
that the detector fits exactly, so a pass shows that the code does what the
emulation did; the five-Gaussian element kernels are the independent model
SIMPLE has had all along, and with B factors they can render everything
the design is about.

Work: the `bfac_pdb` argument of `convolve` and the `pdb_bfac` key of
`simulate_nanoparticle` (section 4); the list parser `parse_element_list`
in `simple_atoms`, the parameters and the commander rules of section 5, the
species names in `_species.pdb` (section 6); the compound-selector repair
and the selectors as species statements (end of section 9).

Exit, besides the checks every phase makes (10.0):

- `simple_atoms_tester`: an element atom rendered with `bfac_pdb` at `B` =
  0 and at `B` = 20 A^2 conserves its voxel sum within 1%, and its peak
  falls by `sum_j a_j / (b_j + bfac + B)^1.5` over
  `sum_j a_j / (b_j + bfac)^1.5` within 1%, with `bfac = (4 lp)^2` the blur
  `convolve` always adds (`(8 smpd)^2` = 8.2 A^2 at 0.358 A without `lp`;
  for Pt and `B` = 20 the ratio is 0.243, not the 0.007 one gets by leaving
  the blur out); the default leaves the B column unread (same voxels as
  before the change); a pseudo-atom is unaffected by the flag.
  `parse_element_list`: `Pt,Ni` gives `Pt`, `Ni` in that order; `Ni,Pt`
  gives `Ni`, `Pt`; `Pt` alone is one symbol and not a list; `Pt,Pt`,
  `Pt,`, `,Pt`, `Pt,Xx` and `Pt,X1` are rejected with a message, without
  stopping, since a fatal error cannot be asserted in process; a list of
  five symbols is accepted and sizes `nspecies` to five.
- `simulate_nanoparticle pdbfile= pdb_bfac=yes` on a PDB of three Pt atoms
  with B factors 0, 10 and 20 A^2 gives three peaks in that order with
  equal voxel sums (recorded); without the key the three atoms are
  identical.
- The phase 0 reference regenerated with the commands of its Progress entry
  (the first run's directory is not relied on): `element=Pt` reproduces its
  numbers; `element=Pt,Ni` on the same volume gives the five present
  products byte-identical to `element=Pt` and a `_species.pdb` with `Pt`
  and `Ni` in the element column.
- The parameter truth table of section 5 holds, checked by hand for every
  row with a fatal outcome and recorded.
- `find_rMax('CdSeW')` returns the wurtzite cutoff and
  `get_lattice_params('CDSEW')` the wurtzite lattice in a unit test;
  `parse_element_list('CdSeW')` returns `Cd`, `Se` and the selector;
  `single_atoms_stats` with `element=CdSeW` passes as a second CTest entry
  (`SIMPLE_CTEST_BUDGET` raised by one), its numbers recorded next to the
  Pt ones; `element=CdSeW discover_species=yes` on that fixture writes a
  `_species.pdb` with `Cd` and `Se` in the element column.

### 10.6 Phase 5: the species test on element ground truth

`species_discovery` is rebuilt. No pseudo-atom fixture remains in it; the
pseudo-atoms keep their unit tests and their role as the detection
template.

Fixtures. The Pt lattice that `simulate_nanoparticle element=Pt moldiam=20
box=160 smpd=0.358` writes (285 atoms); a chosen set of atoms relabelled
with a second element; every atom given `B_i = B_SURF (r_i / r_max)^2` with
`B_SURF` = 10 A^2, an added displacement variance of `B / 8 pi^2` = 0.127
A^2 at the surface (what doubled the variance of a 0.35 A Gaussian in the
emulation; on an element kernel with the default blur the effective
broadening is smaller and element-dependent, and the fitted widths are
recorded rather than compared with a generating sigma); rendered with
`simulate_nanoparticle pdbfile=model.pdb pdb_bfac=yes` at box 160 and 0.358
A (no `lp`); Gaussian noise added in process at a fixed seed, with a
standard deviation of 0.05 of the peak of a core Pt atom measured on the
clean render. Four cases:

| Case | Composition | Runs |
| --- | --- | --- |
| pure | Pt only | `element=Pt`; `element=Pt discover_species=yes`; `element=Pt,Ni`; no element; no element with `discover_species=yes` |
| alloy | Pt, a random quarter Ni | `element=Pt`; `element=Pt,Ni`; `element=Pt,Ni` with half maps; no element; no element with `discover_species=yes` |
| shell | Pt core, the outermost quarter Ni | `element=Pt`; `element=Pt,Ni` |
| light | Pt, a random quarter Al | `element=Pt`; `element=Pt,Al`; no element; no element with `discover_species=yes` |

The runs without an element are the blind mode on element kernels, which
the first run covered only on pseudo-atoms. Half maps are two renders with
independent noise at sqrt(2) times the standard deviation, whose average is
the map.

Ground truth, measured: one isolated atom of each element at `B` = 0,
rendered the same way in a box of 64, integrated (voxel sum times
`smpd^3`), gives the generating intensities and their ratios. From the
coefficients in `convolve` (the integral is proportional to the sum of the
`a_j`), Pt over Ni is 1.65 in the integral and 2.0 in the peak at the
default blur, closer than HAADF contrast will give, so the harder case for
the species call; at that contrast the Ni atoms are expected to be found at
level 1, so the alloy and shell cases test the species call, the fits and
the widths on element kernels, not the residual recovery. The light case is
there for the recovery and for the design's central claim: aluminium (Z =
13) is a diffuse atom whose integral is only 1.84 times below platinum's
while its peak is 3.5 times below. The peak is what the first threshold
sees, so this is the regime of case 1 of section 8, where the first
threshold dropped the light atoms and the residual levels found them; the
integral is what the species call sees, so the classes separate on the
integral ratio, not the peak ratio.

Eligibility for the recall floors is defined without the fit, so that it
cannot be circular: a generating atom is eligible when the clean render,
filtered with the detection template (a Gaussian at `B_ref`), divided at
the atom by the standard deviation of the added noise filtered the same
way, is 6.5 or more. The test computes both from the maps it made.

Floors, written into the test before its first run, none lowered:

- Within a case, the five present products are identical across every run
  with the same `element` value (`files_identical`): the list, the key and
  the half maps do not move level 1. This holds for `Pt` against `Pt,Ni`
  and `Pt,Al`, and for no element against no element with the key.
- Pure case: `K` = 1 and no recovered atom with `element=Pt
  discover_species=yes` and with no element and the key; with
  `element=Pt,Ni`, no recovered atom and the two-class fit reported
  inadmissible in `_species.txt`.
- Alloy, shell and light, with the list: every found atom within `0.3
  d_NN` of a generating atom carries that atom's element; at most one
  false atom (a found atom with no generating atom within `0.3 d_NN`, or a
  second atom on one site); the two-class fit admissible; recall of the
  light element of at least 0.90 over the eligible atoms; the fitted
  intensity ratio within 10% of the measured single-atom ratio;
  `_species.pdb` with one atom per table row and the element column equal
  to the symbol of the table's class.
- Light case: the `element=Pt` run finds fewer than 90% of the Al atoms,
  otherwise the fixture does not test the recovery and the phase stops to
  have it redesigned (a lighter element, or a lower noise); the
  `element=Pt,Al` run finds at least 90% of the eligible ones, so the
  difference is the residual levels' work, recorded by stage and level.
- Alloy and light without an element, with the key: `K` = 2 and every
  label of a found atom right, class 1 being the heavier element.
- Widths: the stage-2 `B` of the Pt class rises from the inner to the
  outer radial shell of `_species_radial.csv`, and the rise lies between
  0.5 and 1.5 of the generated rise (`B_SURF` times the difference of the
  shells' mean `(r / r_max)^2` over the generating atoms in them). A wide
  consistency floor, deliberately; the measured value goes to Progress.
- With half maps: `noise_source` is `half_maps` and the label agreement
  between the halves is at least 0.95.
- The test removes its fixture tree when it ends, pass or fail.

`detect_atoms` runs with the test's thread count. The CTest entry keeps its
name, label and budget; if it exceeds its timeout in Debug the timeout is
raised, not the cases cut. If a floor is not met, the phase stops with the
measured values and, for a false atom, its position, its distance to the
nearest Pt atom and the residual z map's value there: this is the first
contact between the method and non-Gaussian atoms, and the maintainer
decides.

### 10.7 After these phases (not part of a run yet)

Validation beyond the element model: particles passed through the
reconstruction's low-pass filtering and soft mask, with background and
shell-dependent noise; cross-half reproducibility of the labels on real
half maps; an A/B of `detect_atoms` with and without discovery on a real
Pt/Ni map at matched false-positive rate, judged by the maintainer. Then
promotion: recovered atoms into `_ATMS.pdb` with data-derived components in
`_CC.mrc`, class intensities and per-atom widths into `_SIM.mrc`, `element`
optional in the other four programs, `autorefine3D_nano` end to end on a
simulated two-species trajectory, `atoms_stats` without the lattice table.


## 11. Relation to `single_pcg_integration.md`

With pseudo-atoms the columns of `G` are the Gaussians of section 1, `Z`
drops out, and the least-squares amplitudes `c = (G^T G)^-1 G^T x` are the
per-atom amplitudes with neighbour overlap removed exactly. The fit of 3.4 is
a Gauss-Seidel approximation of that solve. Two points in that note need
changing for Z contrast. Its section 5, step 5 makes atoms with
`c_i < 0.2 median(c)` removal candidates, which would flag a whole light
species; the rule should be referenced to the noise. And its statement that
trailing reconstruction blends the simulated map away does not match the
accumulator-domain contract cited in section 0.

## 12. Deliberately absent

- **Element names from the data.** With two classes the intensity ratio
  constrains `(Z_2 / Z_1)^n` for an exponent between about 1.5 and 2, and the
  fitted lattice constant with the class fractions constrains the pair
  through Vegard's law. Both depend on calibrations (the exponent, the pixel
  size) that the map does not carry. Names come from the user (section 5) or
  not at all.
- **Spatial priors on the labels** (3.6).
- **A lattice model.** The gate and the length scales use `d_NN` only. Testing
  the empty sites of a fitted lattice at a lower threshold would reach weaker
  atoms, at the price of assuming a single crystal.
- **Anisotropic widths** in the fit. `calc_anisotropic_disp` stays in
  `atoms_stats`.
- **Per-atom amplitudes in `_SIM.mrc`** (section 4).
- **A bank of template widths or of element templates** in detection
  (sections 5 and 7).
- **Model-generated connected components** (section 6, withdrawn).
- **Fixed detection thresholds** (3.2): they are calibrated from the map.
- **Any change to the present products or the present pruning** before the
  validation of 10.7.

## 13. Open questions

- **Q1.** Is 3-6:1 the atom ratio or the intensity ratio? And what is the
  signal-to-noise of the weakest atoms you can see in a current map
  (section 7)? Together they say whether the light species is in the regime of
  cases 1-4 or of cases 6 and R2.
- **Q2.** Class symbols in the PDB element column (`X1`, `X2`, ...) with the
  class also in the chain: acceptable, or one symbol with the class in the
  chain only?
- **Q3.** Radius from the centre as the radial coordinate, or depth below the
  surface?
- **Q4.** Starting volume without an element: `solve3D_nano`, or the
  ideal-lattice start with a user-supplied lattice constant?
- **Q5.** Is the close-packed neighbour cutoff an acceptable default until the
  pair-distance form is in?
- **Q6.** Stage 2 assumes equal integrals within a class. Is that acceptable
  for surface atoms, or should the stage-1 width be the reported one where
  the two disagree?
- **Q7.** The target of 0.2 expected false atoms per particle behind the
  calibrated thresholds, and the gate on by default: both are prototype
  defaults to be judged on the Pt/Ni A/B of 10.7.

## 14. File table

The files a run may change, and in which phase. A phase edits only its
files; a file the work turns out to need is added here, with its phase and
a one-line reason, before it is edited. Testers of a listed Fortran file
(`*_tester.f90`) are covered by its row. Phase 0 edits nothing but this
note. Phases 0 to 3 are the first run (done); 4 and 5 the second.

| File | What changes | Phases |
| --- | --- | --- |
| `src/main/nano/simple_atoms.f90` | pseudo-atom branch in `convolve`; the symbols accepted by `element_exists`, `set_element`, `guess_an_element`, `Z_and_radius_from_name`; closed-form tests in its tester (phase 1). `bfac_pdb` in `convolve`; `parse_element_list`; their tests (phase 4) | 1, 4 |
| `src/defs/simple_defs_atoms.f90` | pseudo-atom symbols and their sentinel range in `get_element_Z_and_radius` | 1 |
| `src/main/nano/simple_nano_species.f90` | new: mixture fit, BIC, admissibility, `discover_species`, enclosed fraction, threshold calibration; its tester | 1, 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | `species` sub-suite in `suites_single` | 1 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | `unit_single` suite list | 1 |
| `doc/policies/test_environment_policy.md` | the `unit_single` row of its table of fast sub-suites names `species`, so the policy lists what the gate runs (phase 1); the CTest budget and the table of high-level entries name `species_discovery` (phase 3); the `CdSeW` entry (phase 4); the test's description if it changes (phase 5) | 1, 3, 4, 5 |
| `src/main/nano/simple_nanoparticle_utils.f90` | optional kernel width in `phasecorr_one_atom`; `est_nn_dist`; optional `d_NN` in `find_rMax` and `calc_contact_scores` (phase 2); the selector as `character(len=*)` with the symbol derived inside, and the unit test of the compound lattice (phase 4) | 2, 4 |
| `src/main/nano/simple_nanoparticle.f90` | `l_species_free` and the `d_NN` length scales (phase 2); the residual recovery, fits, species call and the files behind `l_discover_species` (phase 3); `species(:)` for the names in `_species.pdb`, `element_key` for the lattice lookups (phase 4) | 2, 3, 4 |
| `src/main/ui/single/single_ui_atom.f90` | `element` optional in `detect_atoms` (phase 2); `discover_species`, `nspecies`, `min_nbrs`, `vol_even`, `vol_odd` on `detect_atoms` (phase 3); `pdb_bfac` on `simulate_nanoparticle` and the list in the description of `detect_atoms` (phase 4) | 2, 3, 4 |
| `src/main/ui/simple_ui_params_common.f90` | only if a new key has to be shared rather than declared on a program | 3, 4 |
| `src/main/params/simple_parameters.f90` | `discover_species`, `nspecies`, `min_nbrs`, `l_discover_species` (phase 3); `element` to `STDLEN`, `species(:)`, `l_species_list`, `pdb_bfac` and `l_pdb_bfac` (phase 4) | 3, 4 |
| `src/main/params/simple_parameters_parse.f90` | registration of the new keys | 3, 4 |
| `src/main/params/simple_parameters_phases.f90` | phase 2 only if the element validation at lines 937-941 needs it for an absent key; phase 3 the derived logicals and the validation of the new keys; phase 4 the list parsing and the validation of `pdb_bfac` | 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_atoms.f90` | `exec_detect_atoms`: command-line validation, half maps (phases 2, 3); the list's normalisation and rules (phase 4) | 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_sim.f90` | `exec_simulate_nanoparticle` passes `bfac_pdb` to `convolve` | 4 |
| `src/main/commanders/test/simple_commanders_test_single.f90` | the high-level test `species_discovery` (phase 3); `single_atoms_stats` parameterised by the crystal system only if its assertions prove fcc-specific (phase 4); `species_discovery` rebuilt on element ground truth (phase 5) | 3, 4, 5 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | its program entry (phase 3); its description if it changes (phase 5) | 3, 5 |
| `src/main/exec/simple_test_exec_single.f90` | its router case | 3 |
| `production/CMakeLists.txt` | its CTest entry; `SIMPLE_CTEST_BUDGET` raised by one (phase 3); the `single_atoms_stats` entry with `element=CdSeW`, budget raised by one (phase 4); the `species_discovery` timeout only if Debug needs it (phase 5) | 3, 4, 5 |
| `doc/implementation_notes/planned/species_discovery.md` | this note: Progress, file-table rows, rulings | each |

## 15. Review of 2026-10-07 and what it changed

The first version of this note was reviewed on 2026-10-07 by a second
reviewer and by the maintainer, who declined to approve it as written and
approved a small, opt-in prototype instead. The findings, what was checked,
and what the note now says:

The pipeline order broke the pruning contract. Fitting and classifying
before `discard_atoms` would have lost every per-atom quantity, because
pruning rebuilds the atom table (section 0, last constraint but one), and
components generated from the model would have made the size test of the
pruning circular. Confirmed from the source. Pruning now comes first and
the recovered atoms are pruned on the atom table (3.9); model-generated
components are withdrawn in favour of data-derived ones at promotion
(section 6).

The aperture intensity lacked a factor `smpd^3` against the continuous
integral it is compared with. Confirmed; added (3.5).

A fixed threshold on local maxima of a z map is not statistically
controlled, the voxel standard deviation does not calibrate the maxima, and
undetected weak atoms can contaminate the noise region. Checked by
simulation (8.8): the voxel tail indeed overpredicts the number of noise
maxima six times, but the region's own maxima calibrate it to within 0.8 to
1.0, and the calibrated threshold equals the fixed one for the emulated
geometry (5.0) and differs for strongly correlated noise (4.7). The
threshold is now calibrated per map (3.2). Contamination is real for a
weak-only domain thicker than 1.5 `d_NN` inside a tight mask (36 to 40% on
the standard deviation) and absent for a surface shell; the half-map
difference is now the noise reference when the halves are given, and the
robust-to-plain spread ratio is reported otherwise. The existing
false-discovery-rate routine of `simple_segmentation.f90` was benchmarked on
the same maps: it admits three to four times the gate's false atoms at
lower recall, so it is not adopted.

The three-neighbour gate is a lattice prior that will disadvantage
surfaces, defects, interfaces and weak-only domains. Checked by simulation
(8.8): at the same threshold the gate costs about six points of recall at
low coordination (six to nine points); at the same false-positive rate it gains 7 to 15 points
everywhere, including the low-coordination bins and a weak-only domain. Its
assumptions and its blind spot are now stated (3.3), and it is a switch with
its use recorded per atom.

The evidence is model-matched: the generator produces the Gaussian atoms the
detector assumes, several production steps were not emulated, and some
threshold numbers came from an earlier variant of the algorithm. Agreed. The
numbers from the earlier variant are replaced (3.3), the independent check
is described for what it is (the preamble, 8.8), and validation on an
independent forward model with a real Pt/Ni comparison gates promotion
(10.7).

What the review asked for and the note now does: the existing path is
unchanged and remains the production path with its regression test; the
recovery runs behind an experimental key and produces diagnostics without
feeding `_SIM.mrc`; components stay data-derived; pruning precedes fitting
and classification; the forward-model validation and the Pt/Ni comparison
at matched false-positive rate come before promotion.

### The plan review of 2026-10-08 (phases 4 and 5)

The plan for the second run was first written as a separate note and
reviewed. The review found seven blocking points, all accepted, and the
plan was folded into this note (the repository keeps one living note per
development) with these changes:

- A list could name five species where the mixture, the posterior columns
  and the chains were sized for three. The review offered a cap or a
  generalisation; the maintainer ruled for the latter: the list has any
  length, `nspecies` is its length, and everything sized by the number of
  classes is sized by `nspecies` (section 5, section 9).
- The peak-ratio check for the B factors left out the resolution blur that
  `convolve` always adds; the stated formula would have failed its own 1%
  test (0.007 against the rendered 0.243 for Pt at `B` = 20). Corrected
  (section 4, 10.5).
- Eligibility for the recall floor used the predicted signal-to-noise of
  section 7, which depends on the fitted width, so a floor on it would have
  been circular on kernels that have no generating width. Eligibility is
  now computed by the test from the clean render and the added noise,
  filtered with the detection template (10.6).
- A list implied discovery only inside `parameters`, after the commander
  had already rejected half maps given without `discover_species=yes`. The
  commander now normalises the raw command line first, and section 5
  carries the truth table of every combination; `nspecies` with a list is
  an error in every case, where the first draft contradicted itself.
- Widening the utilities' element arguments would not have repaired the
  compound-selector truncation: `find_rMax` uses one value for the lattice
  table, where `CDSEW` is valid, and for the atomic radius, where it is not
  and stops. The repair keeps two values apart, the selector for the
  lattice and the symbol for radii and names, with its own regressions;
  first deferred, then brought into phase 4 when the maintainer asked that
  the present species mixtures be handled correctly (end of section 9).
- Removing every pseudo-atom fixture would have removed the only
  high-level coverage of the blind mode and of automatic `K`. The blind
  mode now runs on the element-kernel maps of the alloy and light cases
  without an element (10.6).
- Assigning symbols to classes by atomic number can mislabel: in the
  five-Gaussian model Pt integrates above Au, and the sentinel numbers of
  the pseudo-atoms would sort first. The order given is the order of
  decreasing expected intensity, the program does not reorder, and the
  pseudo-atom symbols are rejected (section 5).

Also taken from the review: the purposes of the light case and of the half
maps are floors, not Progress entries (the `element=Pt` run must miss Al
atoms that the list run finds; admissibility for the mixed cases; a 0.95
half-map agreement floor; product identity for every `element` value, not
only `Pt,Ni`); `B_SURF` is described as an added displacement variance
rather than a doubling that holds only for a standalone Gaussian; the
parser is a non-fatal helper so that its rejections are unit-testable;
the shared help text of `element` stays and the list is documented on
`detect_atoms`; `pdb_bfac` goes through the whole parameter lifecycle; the
test removes its fixture tree.

## 16. Progress

Each phase records here: the date, the files changed, the evidence with log
paths under the run directory, and each exit item of section 10 with how it
was met. The `single_atoms_stats` numbers are recorded for every phase next
to phase 0's.

### Phase 0, baseline (2026-10-08)

Files changed: this note only (this entry). No source file was touched; the
build is of the base revision `92d684c6f` as handed over, with no other
difference in the working tree.

Build and environment. The machine is dell (Oracle Linux 8.10, 24 cores, no
GPU). Compiler gfortran 15.2.1 from gcc-toolset-15
(`/opt/rh/gcc-toolset-15/root/usr/bin/gfortran`, enabled in `~/.bashrc`),
CMake 3.26.5; `python3` on the path is 3.6.8, under which the three gate
scripts run clean. The build is `./compile_debug.sh`: a Debug build
(`-O0 -g -fbounds-check -fcheck=all`) with `BUILD_TESTS=ON`, followed by the
fast gate and the install into `build/`. Every shell that ran SIMPLE set
`SIMPLE_PATH` to that build and took the maintainer's checkout off `PATH`.
Log: `scratch/build_debug.log` under the run directory
(`/home/elmlundho/agent_runs/single_species_proto`); the build, gate and
install took 41 s.

Fast gate (`scripts/run_fast_gate.sh`: test-registry, description and
dependency checks, then the `fast` and `provisional` CTest labels): passed,
15 of 15 entries, 5.3 s of tests. Log: `repo/build/test_runs/ctest_fast.log`.

`single_atoms_stats` (`ctest -R '^single_atoms_stats$'`, with
`SIMPLE_SEED=20260923` set by CTest), Debug build, machine otherwise idle
(load average 3.2 at the start, the tail of the build): passed in 42.6 s of
wall time, far inside its 3600 s timeout, so no Release build is needed.
Log: `scratch/single_atoms_stats_p0.log`. The numbers every later phase is
compared with:

| Quantity | Phase 0 |
| --- | --- |
| atoms in the simulated model (`simatms.pdb`) | 285 |
| atoms detected (`outvol_ATMS.pdb`) | 285 |
| recall within 1 A | 1.0000 |
| precision within 1 A | 1.0000 |
| root-mean-square position error | 0.0148 A |
| correlation of `outvol_SIM.mrc` with the simulated volume | 0.9997 |
| `atoms_stats` rows (atoms) / anisotropic atoms / diameter | 285 / 285 / 19.527 A |

Reference for the equivalence test of phase 2, kept in
`scratch/keep/phase0_ref/` (script `run_ref.sh`, log `run_ref.log`, MD5 sums
of every file in `md5sums.txt`). Both programs are run by `single_exec`
(`simple_exec` refuses them), with `SIMPLE_SEED=20260923`:

    single_exec prg=simulate_nanoparticle element=Pt moldiam=20 box=160 smpd=0.358 outvol=ptvol.mrc pdbout=ptatoms.pdb nthr=8
    single_exec prg=detect_atoms vol1=ptvol.mrc smpd=0.358 element=Pt nthr=8

Kept: the input volume `ptvol.mrc`, the simulated model `ptatoms.pdb` (the
file section 10.1 calls `simatms.pdb`; it is byte-identical to the
`simatms.pdb` of the `single_atoms_stats` run, and `ptvol.mrc` and
`ptvol_ATMS.pdb` are byte-identical to that run's `outvol.mrc` and
`outvol_ATMS.pdb`), and every product of the detection (`ptvol_ATMS.pdb`,
`ptvol_BIN.mrc`, `ptvol_CC.mrc`, `ptvol_MSK.mrc`, `ptvol_SIM.mrc`, and the
intermediate `split_ccs.mrc`). Numbers, from `scratch/keep/compare_pdb.py`
and `scratch/keep/corr_mrc.py` (the same nearest-neighbour matching within
1 A and whole-box correlation as the test):

| Quantity | Reference |
| --- | --- |
| atoms in `ptatoms.pdb` | 285 |
| atoms in `ptvol_ATMS.pdb` | 285 |
| matched within 1 A, each way | 285 and 285 |
| root-mean-square position difference (maximum) | 0.0148 A (0.0340 A) |
| correlation of `ptvol_SIM.mrc` with `ptvol.mrc` | 0.99969 |
| median nearest-neighbour distance, model / detected | 2.7662 A / 2.7486 A |

The detection is deterministic: run again in another directory without
`SIMPLE_SEED`, it wrote `_ATMS.pdb`, `_SIM.mrc`, `_CC.mrc` and `_BIN.mrc`
byte-identical to the kept ones. The median nearest-neighbour distance of the
detected atoms, the quantity phase 2 measures as `d_NN`, is 0.6% below the
2.766 A of the plan's exit for phase 2.

Exit items of 10.1: the fast gate passes (above); `single_atoms_stats`
passes with its present floors (above); the reference is kept under
`scratch/keep/phase0_ref/` and its numbers are recorded here.

### Phase 1, foundations that change nothing (2026-10-08)

Files changed: `src/defs/simple_defs_atoms.f90`,
`src/main/nano/simple_atoms.f90`, `src/main/nano/simple_atoms_tester.f90`,
`src/main/commanders/test/simple_commanders_test_class.f90`,
`src/main/ui/simple_test/simple_test_ui_class.f90`,
`doc/policies/test_environment_policy.md` (added to the file table first:
its list of the fast sub-suites of `unit_single` now names `species`), this
note; new: `src/main/nano/simple_nano_species.f90` and
`src/main/nano/simple_nano_species_tester.f90`.

What was done.

- Pseudo-atom symbols. `X1`, `X2` and `X3` are registered in
  `get_element_Z_and_radius` with the sentinel atomic numbers 901, 902 and
  903, held in the named constants `Z_PSEUDO_FIRST` and `Z_PSEUDO_LAST`
  (clear of the 999 of `CDSE`), radius 1.0 A. `element_exists`,
  `set_element`, `set_name`, `guess_an_element` and
  `Z_and_radius_from_name` all go through that table, so they accept the
  symbols with no change of their own.
- Pseudo-atom kernel. `atoms%convolve` renders an atom whose atomic number
  is in the reserved range as `q (4 pi / B)^(3/2) exp(-4 pi^2 r^2 / B)`, with
  `q` from `occupancy` and `B` from `beta`, inside the same window and cutoff
  logic as the elements. Decision recorded here: the resolution blur that
  `convolve` adds to every element (its `lp` argument, by default a B of
  `(8 smpd)^2`) is not added to a pseudo-atom, whose width is its `B`; the
  closed forms of the phase 1 exit need exactly that, and section 4 defines
  the kernel by `beta` alone. A pseudo-atom with `B` of zero or less stops
  the program. The element branches compute what they did before: the
  constant of eq. B.6 now reaches the sum through a per-atom scale that holds
  the same value, and the same Pt volume and detection are byte-identical to
  phase 0 (below). The unused single-Gaussian helper `egau` and its two
  variables were deleted as dead code. `atoms` gained the getter
  `get_occupancy`.
- PDB round trip. `readPDB` and `writePDB` already carry the element column,
  `occupancy` and `beta`; both numeric columns have two decimals, so a
  pseudo-atom PDB keeps `q` and `B` to 0.005. A fixture whose class intensity
  is 1/6 is therefore rendered at 0.17; phase 3 judges against the values the
  PDB holds.
- `simulate_nanoparticle pdbfile=` renders pseudo-atoms with no change to
  that commander (it reads the PDB with `atoms%new` and calls `convolve`), and
  element PDBs take the unchanged branches.
- `simple_nano_species`: free procedures on arrays, no image or parameters.
  `fit_species_mixture(x, s_meas, nspecies_in, K, labels, post, mu, var
  [, bic, admissible])` fits the one-dimensional Gaussian mixture of 3.6 by
  expectation-maximisation for `K` = 1 to 3, each component variance floored
  at `s_meas^2`, from the sorted-quantile and the largest-gap starts, keeping
  the start of higher likelihood; `BIC = -2 ln L + (3K - 1) ln N`; a `K`
  above 1 is admissible when every class holds at least `max(8, 0.02 N)`
  atoms (counted by largest posterior) and adjacent classes are separated by
  `D >= 3`; the admissible `K` of lowest BIC is taken, or `nspecies_in` when
  it is positive (then even an inadmissible `K`); classes are numbered by
  decreasing mean. Also `class_separation` (`D`), `enclosed_fraction` (`F`
  of 3.5), `maxima_shape` (`(u^2 - 1) exp(-u^2 / 2)`), `calibrate_threshold`
  and `expected_false_count`. The calibration takes `C` as the summed region
  counts over the summed model shape, scaled by the search-to-region volume
  ratio, which is the estimator that reproduces the `C` of 2105 in 8.8, and
  finds `k` by bisection above `sqrt 3`, where the shape falls
  monotonically.

Tests. `simple_test_exec test=unit_single suite=species`: 36 assertions,
all pass (log `scratch/p1_unit/species.log` under the run directory
`/home/elmlundho/agent_runs/single_species_proto`). The samples are drawn at
fixed seeds from the generating parameters each test states, and the fits
are judged against them:

- one class of 400 at spread 0.05: `K` = 1, every label 1, the mean within
  four standard errors, `BIC(1)` equal to the closed form `-2 ln L + 2 ln N`
  computed in the test, `BIC(1)` the lowest; `nspecies` = 2 forces `K` = 2;
- two classes of 200: at `D` = 2.5 the two-class fit is inadmissible and
  `K` = 1; at `D` = 4 it is admissible, `K` = 2, `BIC(2) < BIC(1)`, at least
  95% of labels right (the Bayes error is 2.3%), means within 0.02;
- size floor at `N` = 400: a distinct class of 5 atoms is inadmissible and
  `K` = 1, a class of 8 (the floor) is admitted with every label right;
- three classes (case 12: 264, 132, 132 atoms at 1, 0.5, 0.167, drawn in a
  shuffled order): `K` = 3, `BIC(3) < BIC(2) < BIC(1)`, every label right,
  the class intensities within 0.02;
- `enclosed_fraction` equal to a Simpson integral of the radial density to
  1e-5;
- calibration from the white-noise region counts of 8.8 (154.5, 38.0, 7.2
  above 2.5, 3.0, 3.5; volume ratio 3.628): `C` = 2105.3, `k` = 4.99,
  within 0.1 of 5.0 as the exit asks; the expected false count at `k` equals
  the target; counts that follow the model return its `C`; a doubled search
  volume raises `k`.

`suite=atoms`: 48 assertions, all pass (log `scratch/p1_unit/atoms.log`),
of which the new ones: the three symbols exist, have distinct Z in the
reserved range and radius 1; one pseudo-atom (`q` = 2, `B` = 13.9 A^2,
0.358 A voxels) peaks at `q (4 pi / B)^(3/2)` to 1e-4, reaches 1.43 A and
stops before 2.03 A at a cutoff of `4 sigma` = 1.68 A; off the grid its voxel
sum times `smpd^3` is `q` within 1%; `lp` does not change a pseudo-atom; two
overlapping pseudo-atoms (`B` 13.9 and 20, 2.4 A apart) render as the sum of
the two rendered alone and integrate to `q1 + q2` within 1%; a PDB of `X1`,
`X2`, `X3` written by `writepdb` and read by `atoms%new` keeps symbol,
atomic number, occupancy and `B` (to 0.005) and coordinates (to 1e-3 A).

End to end (`scratch/p1_e2e/`): `simulate_nanoparticle pdbfile=` on a
three-pseudo-atom PDB (box 64, 0.358 A) gave a map whose integral is the sum
of the `q` (1.67000 against 1.67, relative difference -1.1e-7) and whose value
at `X1` is the closed form including the neighbours' tails (0.859593 against
0.859592); its `pdbout` keeps the symbols, occupancies and `B`. The phase 0
commands (`scratch/keep/phase0_ref/run_ref.sh`) run with this build gave
`ptatoms.pdb`, `ptvol.mrc` and all five detection products byte-identical
to phase 0.

A second reader reviewed the diff against this note and found no defect;
its minor points were taken: the log-likelihood returned when the EM
iteration cap is reached now belongs to the returned parameters, a forced
`nspecies` that leaves a class empty is documented at the routine, the
calibration comment states the estimator, and the cutoff check samples a
voxel inside the atom's window.

Noted for phase 3, not changed: `atoms%extract_atom` copies `beta` but not
`occupancy`, so a pseudo-atom extracted with it renders with `q` = 0.

Build and gates. `./compile_debug.sh` (clean Debug build with tests, fast
gate, install): passed, 15 of 15 fast entries, `unit_single` 0.6 s (log
`scratch/p1_build_debug2.log`; gate log `repo/build/test_runs/ctest_fast.log`).
`scripts/check_test_registry.py` and `scripts/check_descr.py` clean; no
compiler warning from the new or touched files; no mode change in
`git diff --summary`.

`single_atoms_stats` on the final build (log `scratch/single_atoms_stats_p1.log`,
load average 4.0 at the start, the tail of the build), against phase 0:

| Quantity | Phase 0 | Phase 1 |
| --- | --- | --- |
| atoms simulated / detected | 285 / 285 | 285 / 285 |
| recall within 1 A | 1.0000 | 1.0000 |
| precision within 1 A | 1.0000 | 1.0000 |
| root-mean-square position error | 0.0148 A | 0.0148 A |
| `_SIM.mrc` correlation | 0.9997 | 0.9997 |
| `atoms_stats` atoms / diameter | 285 / 19.527 A | 285 / 19.527 A |
| wall time (Debug) | 42.6 s | 32.8 s |

Exit items of 10.2: unit tests on one and two pseudo-atoms against the
closed forms (peak, voxel sum within 1% at `4 sigma`, linear addition) pass;
the PDB round trip preserves symbol, occupancy and `B`; the mixture returns
the right `K`, labels and BIC ordering on the fixed samples (one class; two
classes at `D` 2.5, inadmissible, and 4; a class below the size floor; three
classes); the calibration returns 4.99 from the counts of 8.8; the `species`
sub-suite is registered and `scripts/check_test_registry.py` is clean; the
fast gate passes; `single_atoms_stats` is unchanged against phase 0.

### Phase 2, level 1 without an element (2026-10-08)

Files changed: `src/main/nano/simple_nanoparticle.f90`,
`src/main/nano/simple_nanoparticle_utils.f90`,
`src/main/ui/single/single_ui_atom.f90`, this note. Not needed in this
phase, so not touched: `simple_parameters_phases.f90` (its element check
runs only when the key is given) and `simple_commanders_atoms.f90` (with
`element` absent there is nothing to normalise before `params%new`; the
validation of the new keys comes with them in phase 3).

What was done.

- `element` is optional in `detect_atoms` only (`required_override=.false.`).
  The other programs that take it still require it.
- `nanoparticle%new` sets `l_species_free` when `element` is empty. It then
  looks up no atomic number, names the atoms `X1` (the pseudo-atom symbol of
  phase 1) and leaves the atom radius unset until the nearest-neighbour
  distance is measured. With an element, `new` runs exactly as before.
- Template. `phasecorr_one_atom` takes an optional B factor and, when it is
  given, correlates with a unit pseudo-atom of that width instead of the
  element; the species-free path passes `B_ref` = 13.9 A^2 (sigma 0.42 A).
  In the threshold search of the binarisation (`t2c`), which writes the
  trial centres to a PDB and renders them, the species-free path sets the
  B column of the re-read atoms to `B_ref`, because `write_centers` puts the
  per-atom correlation there; `simulate_atoms` does the same, so `_SIM.mrc`
  of a species-free run is the equal-atom simulation at the reference width,
  as section 6 states.
- Length scales. After the first binarisation and before
  `discard_small_ccs`, the species-free path measures `d_NN`, the median
  over atoms of the distance to the nearest other centre (`est_nn_dist`, new
  in the utilities), and sets the atom radius to `0.4 d_NN`. The neighbour
  cutoff of the contact scores is `1.267 d_NN` (an optional `d_nn` argument
  of `find_rMax` and `calc_contact_scores`, which replaces the lattice
  lookup when present), and the contact-score ceiling of `discard_atoms` is
  12, the close-packed value. The exclusion radius of the splitting step is
  `0.21 d_NN` in Angstroms, as decided for this run (section 2: today's
  effective value; the element path's test multiplies squared voxel
  distances by `smpd` once where `smpd^2` is meant, so its written `0.9 r` acts as
  0.59 A for Pt at 0.358 A). Both paths now take the squared radius from one
  variable set at the top of `split_atoms`; for an element it is
  `(0.9 r)^2`, the same expression as before, so the comparisons are the
  same floating-point operations.
- `kill` resets the two new fields. `conv_denoise`, `identify_lattice_params`
  and `fillin_atominfo`, which only the programs that still require
  `element` reach, stop with a message if they are ever called without one,
  instead of falling back silently to the default lattice (fcc,
  `a` = 3.76 A) for the symbol `X1`.

Equivalence test (section 10.3), run by hand on the phase 0 input volume
`scratch/keep/phase0_ref/ptvol.mrc` with the final build, in
`scratch/keep/phase2/` (script `run_phase2.sh`, logs `with_pt/detect.log`
and `no_element/detect.log`, the two `_ATMS.pdb` kept, comparisons in the
`cmp_*.txt` and `corr_*.txt` files), with `SIMPLE_SEED=20260923`:

    single_exec prg=detect_atoms vol1=ptvol.mrc smpd=0.358 element=Pt nthr=8
    single_exec prg=detect_atoms vol1=ptvol.mrc smpd=0.358 nthr=8

| Exit item | Required | Measured |
| --- | --- | --- |
| `element=Pt`: `_ATMS.pdb` against phase 0 | same count, positions within 1e-3 A | byte-identical file (285 atoms); `_BIN`, `_CC`, `_MSK` and `_SIM.mrc` byte-identical too |
| no `element`: atom count against the reference (285) | within 1% | 285 |
| no `element`: positions against the reference `_ATMS.pdb` | matched within 0.5 A, root mean square at most 0.1 A | 285 of 285 matched each way, root mean square 0.0139 A, largest 0.036 A |
| no `element`: `_SIM.mrc` correlation with the input against the reference's 0.99969 | within 0.01 | 0.99488 (difference 0.0048) |
| measured `d_NN` against 2.766 A (Pt, `a` = 3.912) | within 2% | 2.7498 A (0.6% below) |

Against the generating model `ptatoms.pdb` the species-free detection
matches all 285 atoms with a root-mean-square error of 0.0121 A (largest
0.026 A); the element path's is 0.0148 A. The step counts agree on both
paths: no atom split, 285 after splitting, none discarded by contact score
or size, 285 final. The binarisation chose a slightly different threshold
(0.01916 against 0.01929 in map units, a different template), and the
correlation of the trial simulation with the map was 0.867 against 0.861.
The species-free `_SIM.mrc` correlates a little less with the input
(0.9949 against 0.9997) because the input was rendered with the Pt
five-Gaussian kernel and the species-free simulation uses one Gaussian.
The neighbour cutoff came out at 3.484 A against the element's 3.504 A.

A second reader reviewed the diff against this note: the element path
executes as before, nothing on the species-free path reads a lattice or an
atomic number for `X1`, and the radius is set before its first use. Its one
point, that a program still requiring `element` would silently use the
default lattice if it ever got none, was closed by the three stops above;
the equivalence runs were repeated after that change with the same files.

Build and gates. `./compile_debug.sh` on the final source: passed, 15 of 15
fast entries (log `scratch/p2_build_debug2.log`);
`scripts/check_test_registry.py` and `scripts/check_descr.py` clean; no
compiler warning from the touched files; no mode change.

`single_atoms_stats` on the final build (log `scratch/single_atoms_stats_p2.log`,
load average 4.3 at the start, the tail of the build):

| Quantity | Phase 0 | Phase 2 |
| --- | --- | --- |
| atoms simulated / detected | 285 / 285 | 285 / 285 |
| recall within 1 A | 1.0000 | 1.0000 |
| precision within 1 A | 1.0000 | 1.0000 |
| root-mean-square position error | 0.0148 A | 0.0148 A |
| `_SIM.mrc` correlation | 0.9997 | 0.9997 |
| `atoms_stats` atoms / diameter | 285 / 19.527 A | 285 / 19.527 A |
| wall time (Debug) | 42.6 s | 41.7 s |

Exit items of 10.3: the species-free count equals the reference's (285,
within 1%); positions matched within 0.5 A with a root-mean-square
difference of 0.014 A (at most 0.1 A); the `_SIM.mrc` correlation is
within 0.0048 of the reference's (0.01 allowed); `element=Pt` reproduces
the reference `_ATMS.pdb` byte for byte; `d_NN` is 2.7498 A, 0.6% from
2.766 A; the fast gate passes; `single_atoms_stats` is unchanged.

### Phase 3, residual recovery and species call, diagnostics only (2026-10-08)

Files changed: `src/main/nano/simple_nanoparticle.f90` (the discovery),
`src/main/nano/simple_nano_species.f90` and its tester (the width fit, the
separable Gaussian filter, local maxima, robust spread, with tests),
`src/main/params/simple_parameters.f90`, `simple_parameters_parse.f90`,
`simple_parameters_phases.f90` (the keys `discover_species`, `nspecies`,
`nn_gate` and the logicals `l_discover_species`, `l_nn_gate`),
`src/main/ui/single/single_ui_atom.f90` (the three keys and `vol_even`,
`vol_odd` on `detect_atoms`), `src/main/commanders/simple/simple_commanders_atoms.f90`
(`exec_detect_atoms`), the test files `simple_commanders_test_single.f90`,
`simple_test_ui_highlevel.f90`, `simple_test_exec_single.f90`,
`production/CMakeLists.txt` (CTest entry, `SIMPLE_CTEST_BUDGET` 34 to 35),
`doc/policies/test_environment_policy.md` (its file-table row extended to
phase 3 first: the budget count and the table of high-level entries), and
this note. `simple_ui_params_common.f90` was not needed: the new keys are
declared on `detect_atoms` itself.

What was built. With `discover_species=yes`, `detect_atoms` runs the present
detection and writes its five products exactly as without the key; only then
does `discover_species` (a type-bound procedure of the nanoparticle whose
steps are internal procedures named as in section 9: `fit_atoms_joint`,
`build_region`/`noise_stats` for the noise region, `calibrate`,
`detect_residual_level`, `prune_recovered`, `calc_aperture_int`,
`assign_species`, `write_species_report`, `write_radial_profiles`) work on
a copy of the level-1 atom table and the map. The recovered atoms are a
second table (`recovered(:)` of `atom_stats`, freed by `kill`); the level-1
atoms get the new discovery fields in `atominfo` after the products are
written. Nothing reaches `_ATMS.pdb`, `_BIN`, `_CC`, `_MSK` or `_SIM.mrc`.
Its three files are `<map>_species.csv` (one row per atom, level-1 and
recovered: position, radius, coordination, stage and level, detection z,
amplitude, `B` of stage 1 and stage 2, aperture intensity and its
signal-to-noise, class and class posteriors, `valid_corr`, whether the gate
accepted it), `<map>_species_radial.csv` (the profiles of 3.8 per class and
equal-count shell, with the predicted detection signal-to-noise and a flag
below 6.5) and `<map>_species.txt` (the report of section 6 as `key = value`
lines: `K`, class intensities, ratios, fractions and spreads, `D`, the BIC
table with admissibility, the expected misclassification, `d_NN`, `B_ref`
and the stage-1 median `B` of the strong atoms with its ratio to `B_ref`,
the noise source, `sigma_n`, `s_A`, `s_I`, the robust-to-plain spread ratio,
the region size, `k_A`, `k_B`, the search volumes and the expected false
counts, atoms added per stage and level, pruned atoms, the interior
coordination and its deficit, flagged shells, and with half maps the label
agreement and intensity correlation of the two halves).

Decisions taken where the plan leaves the detail open, each recorded here:

- Order. Discovery runs after the present output block rather than before
  it (section 9 lists it before): the products are on disk before discovery
  starts, so their identity with and without the key holds by construction.
- Half maps. `parameters` refuses `vol1` together with `vol_even` and
  `vol_odd`, so `exec_detect_atoms` checks the combination (the three
  discovery keys and the half maps need `discover_species=yes`; both half
  maps or neither), takes the half maps off the command line before
  `params%new` and hands them to the nanoparticle (`set_half_maps`). They
  are masked as the map is.
- `d_NN` on the element path. With an element and `discover_species=yes`,
  `d_NN` is measured at the same point as on the species-free path (after
  the first binarisation) without touching the atom radius.
- Positions are held fixed, as in the emulation: level-1 atoms at the
  present centres, recovered atoms at the centroid of the positive residual
  found at detection (section 3.4 suggests moving them each sweep; not done).
- Width prior. Before any fit there are no strong atoms, so the first
  stage-1 sweep uses a prior centred on `ln B_ref` with spread 0.5; later
  sweeps use the median and robust spread of `ln B` of the strong atoms
  (amplitude above `10 s_A`) within `d_NN` in radius, all strong atoms when
  fewer than five qualify.
- Noise. The region is the voxels of full mask weight (the same
  `mask3D_soft` call as `nanoparticle%new`, on a unit volume) farther than
  `1.5 d_NN` from every atom, rebuilt before each level; the program stops
  below 20 000 voxels. With half maps `sigma_n`, the standardisation of the z
  map and the calibration counts come from the half-map difference over the
  region, as 3.2 says; `s_A` and `s_I` always come from 1 000 phantom sites of
  the map's region, taken at a fixed stride rather than at random so that
  the run stays deterministic. A map without noise stops with a message.
- Calibration once per map, after the level-1 fit, from the plain local
  maxima of the noise field in the region above 2.5, 3.0 and 3.5 (the
  candidates themselves need 3 of their 27 voxels above the threshold, so
  the calibrated count is conservative); `k_B = k_A - 1`.
- Gate. Stage B with the gate counts the atoms accepted before the level,
  so a level adds one shell; without the gate stage B uses the `0.7 d_NN`
  exclusion of stage A. `gate_used` marks the atoms accepted through the
  gate at stage B.
- Pruning: the contact-score rule of `discard_atoms` (its neighbour cutoff
  and ceiling, the outer 15% by distance from the centre of the atom
  positions, fewer than threshold minus one contacts), one pass, on the
  recovered atoms only.
- Predicted detection signal-to-noise (section 7): the amplitude of an atom
  of the class intensity and the shell width over `s_A`, scaled from the
  template width to the atom's as `(sigma / sigma_ref)^(3/2)`, so that it goes
  as `I sigma^(-3/2)`; this reproduces the values of the table in 8.2. The
  first version of the code used the amplitude over `s_A` alone, which goes
  as `I sigma^(-3)`; that was found on reading the first test run (it held
  only 25 of the 71 light atoms of R6 to the recall floor) and corrected in
  the code and the test before any further run. The floors were not
  touched.

The test `species_discovery` (label `highlevel`, its own CTest entry
as section 10.4 asks). It simulates the Pt lattice of phase 0 (285 atoms,
`d_NN` 2.766 A) and builds four pseudo-atom models on it: case 2 (a random
quarter of the atoms light at intensity 1/6, sigma 0.35 A with 5% scatter,
noise 0.03 of the peak of a core strong atom), case 9 (single species, noise
0.05), R3 (single species, variance growing with the square of the radius to
twice its central value at the surface, noise 0.05) and R6 (case 2 with that
radial growth and the light class at 1.3 times the variance of the heavy one
at every radius, noise 0.03). Each model is written as a PDB (`X1`/`X2`,
occupancy the intensity, B column `8 pi^2 sigma^2`; the PDB keeps two
decimals, so the light intensity is 0.17), rendered by
`simulate_nanoparticle pdbfile=`, and given Gaussian noise in process at a
fixed seed: two half maps with independent noise of `sqrt 2` times the
target, whose average is the map. `detect_atoms` (no element) runs three
times per case in directories of its own: plain, with `discover_species=yes`,
and with the half maps as well. The generating model is the PDB as written;
a found atom matches the nearest generating atom within `0.3 d_NN`, an
unmatched found atom or a second found atom on one site is false. The
floors of 10.4 are named constants written before the first run.

A defect of the present detection met on the way, not fixed. In the second
run of the test two cases (case 9 with half maps, R3 with discovery) wrote
an `_ATMS.pdb` and `_SIM.mrc` that differed from the plain run's: one atom
centre moved by 0.29 A in case 9 and by 0.02 A in R3, while `_BIN` and `_CC`
were identical. The centres come from `nanoparticle%find_centers`, whose
OpenMP loop adds every voxel into the shared sums `m(:,i)` and `sum_mass(i)`
of its component without a reduction (the same loop at the base revision);
threads that share a component lose updates. The code that runs before the
products are written is the same with and without `discover_species`, so
this is the present detection disagreeing with itself at 24 threads.
Six further plain runs at 24 threads agreed with each other (the race is
intermittent), and runs at 8 and at 24 threads differ in atom order and by
up to 0.034 A in position (`scratch/p3_race/`). Fixing the loop would change
an existing routine for its existing callers, which this run does not do;
the test runs `detect_atoms` on one thread, where the comparison isolates
the effect of `discover_species`. Recorded for the maintainer. (Fixed on the
Mac after the run, 2026-10-08: the loop carries `reduction(+:m,sum_mass)`
and the test runs `detect_atoms` with its thread count again; the report
has the details.)

Results of the final CTest run (`ctest -R '^species_discovery$'`,
Debug, 505 s of wall time, load average 4.3 at the start, the tail of the
build; log `scratch/species_discovery_ctest.log`; the species files
of every run and the generating PDBs in `scratch/p3_final_evidence/`):

| Fixture, run | Level 1 / recovered (stage A, B) | False | K | Weak recall where predicted SNR >= 6.5 | Ratio (generating 0.170) | Strong shells, largest width error | Light shells, largest width error |
| --- | --- | --- | --- | --- | --- | --- | --- |
| case 2, map | 212 / 73 (73, 0) | 0 | 2 | 71 of 71 = 1.000 | 0.1715 | 0.25% | 1.8% |
| case 2, half maps | 212 / 73 (73, 0) | 0 | 2 | 71 of 71 = 1.000 | 0.1715 | 0.25% | 1.8% |
| case 9, map | 285 / 0 | 0 | 1 | - | - | - | - |
| case 9, half maps | 285 / 0 | 0 | 1 | - | - | - | - |
| R3, map | 285 / 0 | 0 | 1 | - | - | - | - |
| R3, half maps | 285 / 0 | 0 | 1 | - | - | - | - |
| R6, map | 207 / 76 (66, 10) | 0 | 2 | 59 of 60 = 0.983 | 0.1699 | 0.47% | 2.5% |
| R6, half maps | 207 / 74 (59, 15) | 0 | 2 | 63 of 66 = 0.955 | 0.1687 | 0.41% | 3.4% |

Every found atom of the two-species fixtures carries the right label. The
level-1 detection found the strong class only (212 of 214 in case 2, 207 of
214 in R6, the rest of the strong class recovered at stage A) and none of
the light class, as section 3.1 says. Over all light atoms, whatever their
predicted signal-to-noise, 71 of 71 were found in case 2 and 69 and 67 of 71
in R6 (map, half maps); in R6 stage B (the gate) added 10 and 15 of them. The light-class widths come out up to
3.4% low in the outermost shell of R6, where the broadest light atoms are
the ones missed (section 7). The calibrated `k_A` was 5.30 to 5.32 on every
map (search volume 2.14 million voxels against a region of about 1.9
million), `k_B` 4.30 to 4.32; the expected number of noise maxima above
`k_B` in the gated stage-B volume before the gate is applied was 1.1 to 1.2.
`sigma_n` matched the noise put in (0.0444 and 0.0740 against 0.04443 and
0.07405), the robust-to-plain spread ratio was 1.000 to 1.002, the stage-1
median `B` of the strong atoms over `B_ref` was 0.70 for the narrow fixtures
(generating `B` 9.67 A^2 against 13.9) and 1.10 for the radial ones, the
interior coordination 12.0 (11.98 in R6), and with half maps the labels of
the two halves agreed for every atom (intensity correlation 0.98 in the
two-species cases). The five present products were identical with and
without `discover_species` in every fixture and both runs.

A further run with `element=Pt discover_species=yes nn_gate=no nspecies=2`
(the gate switch was renamed `min_nbrs` on 2026-10-08, `nn_gate=no` being
`min_nbrs=0`)
on a noisy case 2 map (`scratch/p3_smoke2/`) found the same 212 level-1
atoms, added 73 atoms at stage A and 5 at stage B without the gate, pruned
those 5 by contact score and returned `K` = 2 with ratio 0.173;
`nspecies=2` without `discover_species=yes` stops with
`nspecies needs discover_species=yes`.

A second reader reviewed the code against this note and found no defect
that changes a result; its points were taken before the final runs: the
gate counted atoms accepted earlier in the same level (now only those of
earlier levels), the pruning used the close-packed ceiling with an element
(now that of `discard_atoms`), non-maximum suppression acted on centroids
(now on the maxima), the voxel count of a candidate excluded the maximum
(now 3 of the 27 voxels), and the test could skip the recall floor when no
light atom reached the predicted signal-to-noise and did not check which
noise source a run used (both now fail the test).

Build and gates. `./compile_debug.sh` on the final source: passed, 15 of 15
fast entries, `unit_single` 0.6 s, the `species` sub-suite 48 assertions
(log `scratch/p3_build_debug.log`); 35 CTest entries, matching the budget;
`scripts/check_test_registry.py` and `scripts/check_descr.py` clean; no
compiler warning from the touched files; no mode change.

`single_atoms_stats` on the final build (log `scratch/single_atoms_stats_p3.log`,
load average 0.9):

| Quantity | Phase 0 | Phase 3 |
| --- | --- | --- |
| atoms simulated / detected | 285 / 285 | 285 / 285 |
| recall within 1 A | 1.0000 | 1.0000 |
| precision within 1 A | 1.0000 | 1.0000 |
| root-mean-square position error | 0.0148 A | 0.0148 A |
| `_SIM.mrc` correlation | 0.9997 | 0.9997 |
| `atoms_stats` atoms / diameter | 285 / 19.527 A | 285 / 19.527 A |
| wall time (Debug) | 42.6 s | 33.0 s |

Exit items of 10.4: for the two-species fixtures the recall of the weak
class is at least 0.90 where its predicted signal-to-noise is 6.5 or more
(1.000, 1.000, 0.983, 0.955), no more than one false atom per particle (0
everywhere), `K` = 2, every label of a found atom right, the intensity ratio
within 10% (largest error 0.9%), the per-shell width of the strongest class
within 3% (largest 0.47%) and of the light class within 10% (largest 3.4%);
for the single-species fixtures `K` = 1 and no recovered atom; for every
fixture the five present products identical to the run without
`discover_species`; the fast gate passes; `single_atoms_stats` is unchanged.

## Appendix. Emulation scripts

`emul_detect.py`: `atoms%convolve`, `phase_corr`, `otsu`, `sortmeans`,
26-connected components, raw-weighted centroids and the bisection of
`binarize_and_find_centers`, written from the source. `emul_species.py`: the
generating model, sections 3.1 to 3.8 and 4 of this note, the three statistics
of 8.4, and the case list (`python3 emul_species.py`, or with case labels as
arguments; about 25 s per case). Seeds are fixed. The template-width bank of
section 7 and the "13 noise peaks returned as a species" observation of 3.3
come from an earlier variant of the residual levels (connected components,
exclusion at `0.5 d_NN`); the gate numbers of 3.3 are from the current
variant (8.8).

`codex_checks.py` (numpy only, kept with the other two): the independent
generator and detector of 8.8, with `python3 codex_checks.py null | contam |
gate | all` (about two minutes for `all`). It ports `detect_peak_thres_fdr`
and implements the threshold calibration of 3.2 as `calibrate_k`.
