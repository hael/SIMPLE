# Atomic B factors in SINGLE model building and refinement

Implementation note, 2026-10-09. Design for review; nothing here is implemented. It builds on
the species-discovery work (`species_discovery.md`), whose second run adds the pieces this note
uses, and it is meant to start after that run is merged.

**Status, after the checks of section 1a: the fitting method is open.** Per-atom B factors
fitted to a reconstructed map are not trustworthy yet, in any form tried. The rendering,
product and refinement parts of this note stand; how B is measured has to be settled first
(section 1a ends with the options).

## 1. Why

Atoms in a nanoparticle become more mobile towards the surface. `detect_atoms` knows where the
atoms are but draws every atom identically when it builds the simulated map `_SIM.mrc`, the
reference that `autorefine3D_nano` aligns the particle images against in the next round. The
mobility gradient that every reconstruction shows is therefore absent from the model and from
the refinement, and the atomic model that SINGLE publishes (`_ATMS.pdb`) carries no B factor at
all: its B column holds a per-atom validation correlation.

The aim is to fit one isotropic B factor per atom, with an occupancy, render the reference with
them, publish them in the model, and carry them through `autorefine3D_nano`. Anisotropic
displacement is not part of this note (section 9).

What a real map shows. On the Pt map used during the species work (`recvol_state01_iter005`,
476 atoms, 0.358 A per voxel), a Gaussian fitted to each atom with its neighbours subtracted
gives a width that rises from B = 19 A^2 in the inner fifth of the atoms to 22.5 A^2 in the outer
fifth (both include the map's resolution), and the outer atoms integrate to 0.91 of the inner
ones. Against that map, a reference drawn with these per-atom widths and intensities correlates
better than one drawn with equal atoms in the bands 5-3.3 A (0.111 against 0.065), 2.5-2.0 A
(0.820 against 0.809) and 2.0-1.6 A (0.426 against 0.403); with equal atoms, the shape of the
atom (the Pt kernel or a Gaussian) makes no difference to these numbers. Section 1a shows that
the gradient read this way is not a reliable measurement.

## 1a. Checks of 2026-10-09: what a per-atom fit measures on a filtered map

Reconstructions are filtered: on the real map the low frequencies are suppressed (its amplitude
against the simulated atoms is 0.37 beyond 20 A and 0.50 at 10-20 A) and the signal ends near
1.1-1.3 A. Every atom therefore carries a negative halo, and in the interior the halos of twelve
neighbours overlap. Numpy checks (scripts in ~/agent_runs/runs/species_validation/numpy_checks, outside the repository):

| Map | Fit | Fitted rise of B, inner to outer fifth | True rise |
| --- | --- | --- | --- |
| Pt particle, same B on every atom, no filtering | Gaussian, neighbours subtracted | -0.7 A^2 | 0 |
| the same, with the real map's filtering | the same | +4.4 A^2 (coordination 12: 19.8, coordination 6: 23.8) | 0 |
| B rising with radius, real map's filtering | the same | +8.1 A^2 | +4.6 A^2 |
| same B on every atom, real map's filtering | Pt kernel through a transfer estimated per shell against the complete model | +0.5 A^2 | 0 |
| B rising with radius, real map's filtering | the same | +5.1 A^2 | +4.6 A^2 |
| the real map, its 476 atoms | Pt kernel, no transfer | +3.0 A^2 | unknown |
| the same | through a transfer estimated against the 476 atoms | +15.8 A^2, surface amplitude 1.8 times the core's | unknown |
| the same | through a transfer estimated against the 476 atoms plus 179 partly occupied surface sites | +7.8 A^2, surface amplitude 1.16 times | unknown |

What this says:

- A per-atom fit on a filtered map, with neighbours modelled without the filtering, reads the
  overlap of negative halos as a mobility gradient: on a particle with no gradient it reports
  one of 4.4 A^2. The gradient read from the real map in section 1 is of that size.
- Fitting through the map's transfer function removes that artefact when the model is complete
  (simulated particles, every atom known).
- On real data the model is not complete (the partly occupied surface sites that detection
  prunes on purpose, and whatever else the map holds), and the transfer estimated against it
  absorbs the missing density at the lattice frequencies: the real map's gradient comes out as
  3.0, 7.8 or 15.8 A^2 depending only on what the model contains. This is the failure of the
  withdrawn phase 6 of the species work, on a single-species particle.

So no per-atom fit on the reconstructed map gives a trustworthy B profile yet. Options, to be
decided before section 4 is built:

1. **Measure B in the refinement, not on the map.** The profile of section 3 has three
   parameters per particle. Instead of fitting atoms to a filtered map, choose the profile by how
   well references drawn with it match the particle images in `refine3D_nano` (its scores against
   the images, before any reconstruction filtering). Three parameters judged against the raw data
   avoid both the halo artefact and the incomplete-model problem. It is the direction this note
   would take.
2. **Complete the model and fit through the transfer**, with the partly occupied sites as atoms
   of their own occupancy. The table shows this moves the answer a long way, but there is no
   check on when the model is complete enough.
3. **Profile only, from the map, with an artefact correction** computed by fitting a simulated
   particle of the same lattice, size and filtering with uniform B, and subtracting its apparent
   gradient. Cheap, but it assumes the filtering is known and the same everywhere.

- `atoms%convolve` can add each atom's B factor to the element kernel (the argument `bfac_pdb`):
  in the five-Gaussian parametrisation the Debye-Waller factor adds B to every width term, the
  same operation `convolve` already uses for the resolution blur. It does not yet scale an
  element atom by its occupancy.
- `simulate_nanoparticle pdbfile= pdb_bfac=yes` renders a PDB with per-atom B factors, which
  makes simulated particles with a mobility gradient possible.
- The species-discovery branch of `detect_atoms` fits each atom's width and amplitude with its
  neighbours subtracted, with a prior that ties a weak atom's width to the strong atoms at the
  same radius. It runs only under `discover_species=yes`, and its results reach only the
  diagnostic `_species` files. (A fit through an estimated transfer function of the map was
  tried in that work and withdrawn: a transfer estimated against an incomplete model absorbs
  unseen atoms' density, and the extra parameters outweighed the gain. This note does not use
  one.)
- `atoms_stats` estimates an isotropic displacement per atom (`calc_isotropic_disp`) by a
  log-linear fit to the positive voxels within three quarters of half the nearest-neighbour
  distance. That estimator is biased in noise (at a noise of 0.03 of the peak it returned a
  width 6% too large in the emulation of the species note). Its result is written as `u_iso`,
  a variance in A^2, under the CSV column `BFAC`.
- On the present path `simulate_atoms` writes the generalised coordination number into the B
  column of the atoms it renders and `convolve` ignores it; `_ATMS.pdb` carries `valid_corr` in
  its B column, which `valid_corr_in_bfac_field.pdb` also keeps. Nothing in SIMPLE reads the B
  column of `_ATMS.pdb` back (`atoms_stats` reads only the coordinates).

## 3. The model

Each atom i has a position (from detection, unchanged), an isotropic B factor `B_i` and an
occupancy `q_i`, at most 1. With an element the atom is the element's five-Gaussian kernel
with `B_i` added to each width; without one it is the Gaussian of the species note at width
`B_i`; in both cases scaled by `q_i`. No transfer function of the map enters the model: the fit
compares this kernel with the map directly, as the hand fit of section 1 did. `B_i` therefore
includes the map's resolution as well as the atom's mobility; what the fit reports reliably is
the gradient from core to surface, and the absolute value with an element is an upper bound on
the Debye-Waller factor (section 7).

Occupancy and B are separable to the extent that width and integral are: B sets the width, the
occupancy the integral relative to a fully occupied atom of the species. Fitting B alone would
read a partly occupied surface site, which real maps have (the outer shell of the species work),
as a more mobile one, and overstate the gradient. Both are fitted; the reference of a fully
occupied atom is the median integral of the inner 30% of the atoms (for one species) or the
class intensity (with several).

**The radial profile.** A single weak surface atom's width is poorly determined: in the
emulation of the species note the fitted width of a weak atom ranged over a factor of 2.5. The
physics says the mobility rises towards the surface, so each particle gets a smooth profile,
`ln B(d) = c0 + c1 d + c2 d^2` in the depth `d` below the surface (the distance from the atom to
the convex hull of the level-1 atoms, which follows a non-spherical particle better than the
distance from the centre), fitted to the strong atoms by weighted least squares, and each atom's
`ln B_i` is drawn towards it with a Gaussian prior whose spread is the robust spread of the
strong atoms around the profile (floor 0.1). A strong atom's own density dominates its prior; a
weak one's does not. The profile is reported: it is the per-particle summary of the mobility
gradient.

## 4. The fit

After level-1 pruning, on the present path, opt-in (section 6):

1. Gauss-Seidel sweeps over the atoms (three, as in the species note): for each atom, the map
   minus the background minus every other atom's current model, inside a sphere of half the
   nearest-neighbour distance; the amplitude in closed form and `ln B_i` by a one-dimensional
   search with the profile prior; `q_i` from the amplitude, capped at 1. After each sweep the
   background is re-estimated by smoothing the residual, and the profile is refitted.
2. Positions stay as detection gives them. Refining them jointly with B is possible but not part
   of this note (section 9).

The species-discovery branch uses the same routines; when both are on, the species fit (with its
class tie) supersedes this one for the atoms it covers, and the reported B factors are the
species fit's.

## 5. The products

- `_SIM.mrc`: every atom drawn with its own `B_i` and `q_i`. With an element through `convolve` with `bfac_pdb` and the occupancy scaling of
  section 2 added to the element branch; without one through the pseudo-atom branch, which already
  reads both.
- `_ATMS.pdb`: B column `B_i` in A^2, occupancy column `q_i`. `valid_corr` stays in
  `valid_corr_in_bfac_field.pdb`. Release 4 keeps no compatibility, and no reader in SIMPLE is
  affected.
- `_bfactors.txt`: the profile coefficients, the profile evaluated at the core and at the
  surface, and the spread of the strong atoms around the profile.
- `atoms_stats`: takes `u_iso = B_i / 8 pi^2` from `_ATMS.pdb` when its B column holds fitted
  factors, instead of `calc_isotropic_disp`, which is deleted. The anisotropic estimate
  (`calc_anisotropic_disp`) is unchanged.

## 6. Keys and the refinement loop

- `detect_atoms fit_bfactors=yes|no`, default `no` until section 8 has run. With `no` every
  product is byte-identical to today.
- `autorefine3D_nano fit_bfactors=yes|no` passes the key to each round's `detect_atoms`; each
  round's reference then carries that round's B profile. Every round's `_bfactors.txt` is copied
  with the other per-round products, so the profile's trajectory over the rounds is visible.
- Each round fits from scratch; nothing is carried between rounds. The risk is the usual one of
  model-based refinement: a reference with a strong gradient pulls alignment towards maps with
  that gradient. Section 8 measures it against the truth instead of assuming it away.

## 7. Limits

- `B_i` includes the map's resolution: the fit compares the kernel with the map as it is. The
  gradient from core to surface is what the fit measures reliably; the absolute value with an
  element is an upper bound on the Debye-Waller factor, and without an element it is a width
  relative to the Gaussian of the detection template.
- Occupancy and B separate only as well as width and integral do at the map's resolution. Near
  the surface, where both change at once, the two are correlated; the per-atom covariance from the
  fit is reported with the values.
- The surface layer is where the map is worst. The profile prior is what keeps surface B factors
  sensible; atoms the present pruning removes get no B factor.

## 8. Validation

1. **Simulated refinements with a known gradient.** The `single_workflow` high-level test already
   simulates a Pt particle, projects it, adds a trajectory's noise, runs `analysis2D_nano` and
   `autorefine3D_nano`, and measures the pair pose error against the true orientations. A new
   case renders its particle with a B profile rising towards the surface (with `pdb_bfac=yes`) and
   runs `autorefine3D_nano` twice, with `fit_bfactors=no` and `yes`. Measured against the truth:
   the pair pose error of both runs (small noisy tests are judged by orientation errors against
   the truth, not by FSC), and the fitted B profile and occupancies against the generating ones.
   Floors, set before the first run: the fitted profile within 15% of the generating one at the
   core and at the surface; the pose error with B factors no worse than without.
2. **The 66 real Pt maps** (`/data/NanoX/zenodo`, through the species validation run): per map
   the profile and the occupancies; per particle, whether the profile is stable across its time
   points; the reference correlation by resolution band with and without B factors.
3. **Present path:** with `fit_bfactors=no`, products byte-identical to today; the fast gate;
   `single_atoms_stats` unchanged.

## 9. Deliberately absent

- Anisotropic displacement parameters in the fit or the reference.
- Joint refinement of positions with B factors.
- A B-factor-dependent pruning rule: pruning stays as it is.
- Carrying B factors between rounds of `autorefine3D_nano`.

## 10. Phases

1. Occupancy scaling in the element branch of `convolve` (with `bfac_pdb`), with unit tests; the
   profile fit and the per-atom fit as routines on arrays in `simple_nano_species`, with unit
   tests on samples drawn from a known profile.
2. `fit_bfactors` in `detect_atoms`: the fit on the present path, the products of section 5,
   `atoms_stats` on the fitted factors. Exit: with `no` byte-identical; with `yes`, on a simulated
   particle with a known profile, the profile within 15% at the core and the surface and the
   occupancies within 0.1.
3. `autorefine3D_nano`: the key passed through, the per-round products, and the `single_workflow`
   case of section 8. Exit: the floors of section 8.
4. The real maps, inside the species validation run.

## 11. Open questions

- Is depth below the convex hull the right coordinate for the profile, or radius from the centre
  for the near-spherical particles SINGLE sees most?
- Should the reference also carry occupancy, or only B (occupancy shown in the model only)?
- When the validation is done, should `fit_bfactors=yes` become the default for
  `autorefine3D_nano`?
