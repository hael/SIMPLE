# The atomic density model: species-agnostic atoms that render like the map

Design note, 2026-10-25. Nothing here is implemented. It is the first step of the
species-agnostic rebuild of `detect_atoms` that the maintainer set as the top priority
on 2026-10-25: get the atom detection and the simulated atomic density right on the map
alone, before any species call, refinement loop or reconstruction prior is touched.

Terms. SINGLE is the nanoparticle workflow of SIMPLE; `detect_atoms` finds atoms in a
3D map and writes `_ATMS.pdb` (the atoms), `_SIM.mrc` (the map rendered from them) and
other products. The B factor `B` of an atom sets its width (a Gaussian of variance
`B / 8 pi^2` per axis). The occupancy of a site is the fraction of the time an atom is
there; in a map it scales the atom's integral. The coordination number `c` of an atom
is its number of neighbours within the bond cutoff. `d_NN` is the median
nearest-neighbour distance measured from the first detection. The halo is the negative
ring every atom carries in a reconstruction whose low frequencies are suppressed. The
noise region is the part of the map far from every atom.

## 1. What is asked, and the rulings that shape it

The maintainer's rulings of 2026-10-25:

- `detect_atoms` works on a map alone. Half maps are typically not available in this
  field, and the code must take reconstructions made outside SINGLE (GENFIRE, other
  tomography codes). Half maps may enter the later species step, as a non-independent
  split aligned to one reference, never as a requirement.
- Per-atom intensities are not measured on the images. In liquid-cell EM they are too
  noisy; for tomographic maps the reconstruction is often all we get. The map is the
  measurement domain. This overrides the "measure B in the refinement" option of
  `atom_bfactor_model.md`, which is withdrawn by this note.
- Order of work: species-agnostic detection and the simulated atomic density first;
  the species call (`detect_atom_species`), `autorefine3D_nano` and the atomicity prior
  only once the rendered model fits the experimental maps.
- The model is the minimal one that carries every parameter species detection needs,
  and a map rendered from it must show the radial fade of the atomic densities that
  every experimental map shows.

Rulings on the first draft, same day:

- The coordination number is the covariate of the profiles (not the depth below the
  surface).
- The gap between the species-free detection and the Pt kernel on the real maps (47
  against 56 of 66 maps matched) is tested two ways in the validation: a second
  detection pass with the fitted kernel, and the mask factor (1.2 times the measured
  diameter).
- On a map without an empty region the noise reference is the residual between the
  atoms.
- The residual recovery of the species branch joins the present path in this step:
  this workflow must identify and model every atom. The atoms it adds are labelled
  abstractly, as weak atoms, with no species meaning.

Two measured facts from the species work set the bar (`species_discovery_status.md`,
`atom_bfactor_model.md`): on the 66 real Pt maps an atom's intensity rises by about 4 to
6 percent per neighbour while its fitted width stays flat, and a per-atom Gaussian fit
on a filtered map invents a core-to-surface B rise of 4.4 A^2 on a particle generated
with one B, because the halos of twelve neighbours overlap in the interior and of five
at the surface. A model that ignores the halo cannot fit the map, and one that fits the
halo freely per shell absorbs missing atoms (phase 6 of the species work, withdrawn).

## 2. The model

The map as read, `m(x)`, is modelled as

    m(x) = b(x) + sum_i I_i [ K(x - r_i; B_i) - a K(x - r_i; B_i + B_h) ] + n(x)

    K(x; B) = (4 pi / B)^(3/2) exp(-4 pi^2 |x|^2 / B)        (unit integral)

Each atom `i` has a position `r_i`, an integrated intensity `I_i` and a width `B_i`.
The bracket is one atom as the reconstruction shows it: its Gaussian less a wider,
shallower Gaussian. In Fourier space the atom's transform is
`I_i exp(-B_i s^2 / 4) H(s)` with the transfer

    H(s) = 1 - a exp(-B_h s^2 / 4)

which is 1 at the lattice frequencies and `1 - a` at zero frequency. Two numbers, `a`
and `B_h`, describe the reconstruction's suppression of the low frequencies and the
halo it leaves between atoms. On the Pt map used in the species work the amplitude
against simulated atoms was 0.37 beyond 20 A and 0.50 at 10 to 20 A, which this form
reproduces with `a` about 0.6 and `B_h` about 200 A^2 (a halo of standard deviation 1.6
A, 0.6 `d_NN`); the interstitial wells at minus two noise standard deviations are the
sum of the neighbours' halos. `a` and `B_h` are global to the map, and they are
estimated on the interior atoms only, where every neighbour is present, so that
missing surface atoms cannot leak into them. That is the difference from the withdrawn
per-shell transfer.

`b(x)` is a smooth background, the residual smoothed at the scale of `d_NN`. `n(x)`
is noise of standard deviation `sigma_n`.

The per-atom parameters have structure across the particle, and that structure is what
species detection needs. With one species of intensity `S`:

    ln I_i = ln S + beta (c_i - 12) + e_i,      e_i ~ N(0, tau_I(c_i)^2)
    ln B_i = ln B_p + gamma (12 - c_i) + f_i,   f_i ~ N(0, tau_B^2)

`beta` is the environment slope: the relative loss of intensity per missing neighbour,
the same for every species (on the FePt tomogram Fe and Pt lose the same 5.2 percent per
neighbour; the multiplicative form is therefore the right one, and it replaces the
additive term now in `fit_species_mixture`). `exp(beta (c_i - 12) + e_i)` is the atom's
occupancy relative to a fully coordinated atom of its species. `tau_I(c)` is the spread
of occupancy at coordination `c`, larger at the surface than in the interior; it is a
measured quantity, and it is the error bar every surface species call will carry.
`gamma` is the width slope towards the surface and `B_p` the width of a fully
coordinated atom; `B_i` includes the map's resolution, so only its variation across the
particle is physical. If the data say `gamma` is zero the model loses nothing.

With `K` species the only change is `S -> S_k` with class fractions `pi_k`; `beta`,
`gamma`, `tau_I` and `tau_B` stay shared. The species posterior of an atom is then

    P(k | I_i, c_i) ~ pi_k N( ln I_i ; ln S_k + beta (c_i - 12), tau_I(c_i)^2 + v_i )

with `v_i` the measurement variance of `ln I_i`. That belongs to `detect_atom_species`
and is written here only to show that nothing is missing from the table `detect_atoms`
must produce: `r_i`, `I_i`, `B_i`, `c_i`, `v_i`, and the particle's `beta`, `gamma`,
`B_p`, `tau_I(c)`, `tau_B`, `S`.

Parameter count: five per atom, five per particle (`S`, `beta`, `gamma`, `B_p`, and the
spreads), two for the reconstruction (`a`, `B_h`), the background field and `sigma_n`.
Every atom, strong or weak, is one instance of the same five parameters; the label
records how it was found and enters no equation.

## 3. What is species-agnostic in it

Nothing in section 2 names an element. The kernel is a Gaussian; the length scales
come from `d_NN` (atom radius 0.4 `d_NN`, bond cutoff 1.267 `d_NN`, the close-packed
contact ceiling of 12), as the species-free path of `detect_atoms` already does; the
intensity scale `S` is the median `I` of the interior atoms. The five-Gaussian element
kernels of `atoms%convolve` stay for `simulate_nanoparticle` and for `atoms_stats`,
which still needs an element for the lattice fit; they are not used in detection or in
the model. What an element would add, the absolute Debye-Waller factor and a Z-scaled
intensity, is exactly what the map cannot calibrate without one.

## 4. The fit

The fit runs on the present path after the level-1 pruning (`discard_atoms`), with no
element, and it now includes the residual recovery, so that every atom in the map is
identified and modelled. The order is: level 1 and its pruning as today; a provisional
fit of the model to the level-1 atoms; the residual recovery; the pruning of the
recovered atoms; the full fit; the products. The atoms of level 1 are the strong
atoms; the atoms the recovery adds are the weak atoms. The labels say how an atom was
found and nothing about its species.

The recovery is the species branch's, moved as it is: the residual (map minus
background minus the model) is matched-filtered with the kernel and standardised to a
z map; candidates are its local maxima above a threshold calibrated from the count of
noise maxima in the noise region, so that the expected number of false atoms is 0.2 per
particle; a first stage searches the whole mask, a second stage, one unit lower, only
where a candidate has at least `min_nbrs` found atoms within 1.15 `d_NN`; each new atom
is fitted on the residual and the levels repeat until nothing is added. One difference
from the species branch is expected to help: with the halo in the kernel, the residual
no longer holds the halos of the strong atoms, which were part of the noise the
matched filter saw. Weak atoms are pruned by the present policy, by the ruling of
2026-10-09 of the species work: a weak atom's contact score counts strong atoms only,
and the outer zone and the threshold are those `discard_atoms` derived for the
level-1 set. Weak atoms have no connected component; `_BIN.mrc` and `_CC.mrc` are the
level-1 products.

The fit reuses the routines the species branch of `nanoparticle` has
(`fit_gauss_width`, the Gauss-Seidel sweeps of `fit_atoms_joint`, `update_background`,
the noise region, `noise_stats`, `compute_zmap`, `calibrate`, `detect_residual_level`,
`prune_recovered`), with these changes:

1. **Kernel.** The atom rendered and fitted is the difference of Gaussians of section 2,
   not the plain Gaussian. `render`, `gather` and `fit_gauss_width` take `a` and `B_h`.
2. **Globals first.** `a` and `B_h` are fitted by a two-dimensional search of the
   squared residual over the union of spheres of radius `d_NN` around the interior atoms
   (`c_i = 12`), the atoms held at their current parameters; once after the first sweep
   of the provisional fit and once in the full fit. Start: `a = 0.5`,
   `B_h = 8 pi^2 (0.6 d_NN)^2` (about 215 A^2 for Pt).
3. **Per-atom sweeps.** Three sweeps, as now, over strong and weak atoms alike: for
   atom `i` the data are the map minus the background minus every other atom's model,
   in a sphere of radius `d_NN / 2`; the amplitude is linear, `ln B_i` is a
   one-dimensional search, and the position moves to the centroid of the atom's own
   residual (bounded by 0.35 `d_NN`). The width prior is the profile of section 2,
   `ln B_p + gamma (12 - c_i)` with spread `tau_B`, in place of the shell median of the
   species branch. A weak atom's own density rarely fixes its width; the profile does,
   which is what the species branch borrowed from the strong atoms of the same shell.
4. **Profiles after each sweep.** `beta`, `S` and `tau_I(c)` by weighted regression of
   `ln I_i` on `c_i` over the atoms with amplitude above ten times its noise, the
   spread in bins of `c` (at least 20 atoms per bin, the interior bin as the floor);
   `gamma`, `B_p`, `tau_B` the same way from `ln B_i`. Coordination counts are the
   contact scores on the current positions of all atoms, strong and weak, so that a
   strong surface atom next to a recovered weak one is counted as the neighbour it has.
5. **Noise.** `sigma_n` and the measurement noise of `I_i` and `B_i` come from the map
   itself: the standard deviation of the noise region when the map has one, and the
   same estimators applied at phantom sites in it, as now. When a map has no empty
   region (a support-masked tomogram that fills its box, the FePt case), the noise
   reference is the residual between the atoms: the robust spread of the fit residual
   inside the envelope farther than `d_NN / 2` from every atom, and the phantom sites
   are drawn there. It is reported as such; misfit adds to it, so it is an upper
   bound, and the recovery thresholds calibrated from it are on the conservative side.

Cost: a sphere of `d_NN / 2` holds about 250 voxels at 0.358 A, so three sweeps over
500 atoms take about a second; the FePt particle (23,000 atoms) stays well under a
minute. The global search evaluates the model on the interior union, of the order of
a hundred thousand voxels, a few dozen times.

## 5. Products

Release 4 keeps no compatibility, so the products change in place:

- `_ATMS.pdb`: every atom, strong and weak; occupancy column
  `exp(beta (c_i - 12) + e_i) = I_i / S`, B column `B_i`; chain `A` for the strong
  atoms and `B` for the weak ones, so that a viewer selects either set by chain
  (the label is also a column of the table). `valid_corr` leaves the B column; it
  stays in `valid_corr_in_bfac_field.pdb` and in the per-atom table.
- `_SIM.mrc`: `sum_i I_i K(x - r_i; B_i)`, the atoms without the transfer and without
  the background. This is the map that shows the fade: fainter and broader atoms
  towards the surface because `I_i` and `B_i` say so. It is the reference
  `autorefine3D_nano` will align against later (refinement filters itself).
- `_FIT.mrc`: the full model of section 2, `b + sum_i I_i [K - a K]`, to be compared
  with the input; `_RES.mrc = m - _FIT`, the residual, where missing atoms show as peaks
  and false ones as wells.
- `_atoms.csv`: one row per atom with `r_i`, `I_i`, `B_i`, `c_i`, the measurement
  spreads of `ln I_i` and `ln B_i`, the occupancy, `valid_corr`, the label (strong or
  weak), and for a weak atom the stage, level and detection z that accepted it.
  `_model.txt`: `S`, `beta`, `gamma`, `B_p`, `tau_I(c)` by bin, `tau_B`, `a`, `B_h`,
  `sigma_n` and where it was measured, the residual-to-noise ratio inside the
  envelope, `d_NN`, the counts of strong and weak atoms, the recovery thresholds and
  the expected false count behind them, the pruning zone and threshold.

`_BIN.mrc`, `_CC.mrc` and `_MSK.mrc` are the level-1 products, unchanged. The
`_species` files stay with `discover_species` until `detect_atom_species` takes them
over; the recovery no longer belongs to that key, so `discover_species=yes` reduces to
the class call on the table this step writes.

## 6. Validation

Written before the run, as the acceptance checks.

1. **Synthetic recovery.** Particles rendered from the model with known `beta`,
   `gamma`, `a`, `B_h`, with white and with low-frequency noise at the real maps'
   level, including the uniform-B case (`gamma = 0`). The fit returns `beta` within
   0.005 per neighbour, `gamma` within the equivalent of 1 A^2 from core to surface
   (the per-atom fit without the transfer gives 4.4), `a` within 0.05 and `B_h` within
   20 percent. Unit test on arrays in the new tester; a high-level case on a rendered
   particle through `detect_atoms`.
2. **The 66 real Pt maps** (`/data/NanoX/zenodo` on dell, maps and stored coordinates,
   no half maps). Per map: the residual inside the envelope over `sigma_n` (target below
   1.3; the equal-atom render of today is the baseline to beat); the shell-mean
   intensity profile of `_FIT.mrc` against the map's within noise (the fade
   reproduced); the depth of the interstitial wells reproduced. The strong atoms
   against the stored coordinates: unchanged from today's species-free path or every
   difference named. The weak atoms reported on their own: how many per map, where
   they sit (distance to the nearest stored atom, coordination, radius), their
   occupancy and width; the species validation found about 20 per map one bond length
   from a stored atom, mostly in the outer half, and this run says what they are. Per
   particle, across its 2 to 8 time points: `beta`, `gamma`, `a`, `B_h` stable, and the
   per-atom `I_i` of matched atoms correlated beyond what noise allows.
3. **The present detection gap.** Today the species-free level 1 matches the stored
   coordinates on 47 of 66 maps, the Pt kernel on 56. Two causes are tested in the
   same run, each on all 66 maps and on the ten short maps in particular: a second
   detection pass with the fitted kernel (the difference of Gaussians at `B_p`, `a`,
   `B_h`, replacing the Gaussian at `B_REF_PSEUDO`); and the mask factor, 1.2 times
   the measured diameter today, run at 1.0, 1.1, 1.2 and 1.3 (on ESC1_time28s the 1.2
   mask matched 85 percent of the stored atoms and a smaller one 99.8 percent). The
   two are crossed, a table of eight, and the matched fraction per map decides.
4. **FePt** (`~/Downloads/TomoFromPeter/for_simple`): the fit on a map from outside
   SINGLE with no empty region, the noise from the residual between the atoms; `beta`
   against the 5.2 percent per neighbour measured on their model; the weak atoms and
   the residual map against their traced-but-removed peaks.

## 7. Code map

- New module `src/main/nano/simple_nano_density.f90`: the kernel (`K` and its
  difference form), the per-atom fit (from `fit_gauss_width`), the profile regressions,
  the global transfer search, the rendering of `_SIM` and `_FIT`; with
  `simple_nano_density_tester.f90`. The array routines of `simple_nano_species.f90`
  that it takes over move; `simple_nano_species` keeps the mixture and the threshold
  calibration for the later workflow.
- `simple_nanoparticle.f90`: a procedure `fit_density_model` on the present path after
  `discard_atoms`, before the products, holding the provisional fit, the recovery
  (the contained routines of `discover_species` that do it, moved), the pruning of
  the weak atoms and the full fit; the atom table (`atominfo`) grows by the weak
  atoms, each with its label, so that every writer sees one set. `simulate_atoms`
  renders from the per-atom table through the pseudo-atom branch of
  `atoms%convolve`, which already reads an occupancy and a B per atom (the halo term
  is a second Gaussian in that branch). `discover_species` loses its stage 1 and its
  recovery and keeps the class call on the fitted `I_i` (the aperture sum stays in
  the table as a cross-check during validation) until `detect_atom_species` takes it.
- `write_centers` and the csv writers: the new columns and files.
- Tests: `unit_single` gets the density sub-suite; `single_atoms_stats` entries and
  `species_discovery` are updated to the new products (CTest budget unchanged).

## 8. Deliberately absent

- Element kernels in detection or in the model.
- Intensities from the particle images or the projections.
- Half maps as a requirement.
- A free per-shell transfer (withdrawn phase 6); the transfer has two global numbers.
- Anisotropic widths, a lattice model, spatial priors on labels.
- The species call, `autorefine3D_nano` and the atomicity prior: they come after this
  model fits the 66 maps.

## 9. Open questions

- The weak-atom label in `_ATMS.pdb`: the chain letter (proposed above) or a different
  carrier; the per-atom table carries it either way.
- Whether the neighbour gate of the second recovery stage stays on by default
  (`min_nbrs` = 3) on the present path. It was the species prototype's default; the
  false-atom count of the 66-map run is where it gets judged.

## 10. Phases of the unattended run

One run on dell, phases 0 to 4, each gated by the tests named with it. Floors are
written here before the run; none is relaxed by the run. The real maps are
`/data/NanoX/zenodo` (read only), run with `single_exec prg=detect_atoms vol1=<map>
smpd=0.358` and nothing else unless a phase says so (the mask is the program's own).
"Matched" means the fraction of a map's stored atoms with a detected atom within 1 A,
as the species validation counted it.

**Phase 0, baseline.** No source change. Build the base (`./compile_clean.sh`, the fast
gate, the install), run `species_discovery` and the four `single_atoms_stats` entries
and record their numbers and times. Keep the base build for the later comparisons.
Then the baseline on the 66 maps, with the base build: `detect_atoms` without an
element and with `element=Pt`, matched fraction and atom count per map against the
stored coordinates, wall time per map; and the mask-factor scan on the ten maps the
species validation found short (matched below 0.98 without an element): the measured
diameter from the log, then `mskdiam` set by hand to 1.0, 1.1 and 1.3 times it (1.2 is
the program's own), without an element. Exit: the fast gate and the five entries pass;
the two baseline tables and the scan table are in the run's keep directory and
summarised in Progress, with the 47 and 56 of the species validation reproduced or
every difference named (the size rule of `discard_atoms` was fixed since, so removals
may differ).

**Phase 1, the model as a module.** `simple_nano_density.f90` with its tester: the
kernel and its difference form, rendering of a set of atoms into an array, the
per-atom fit with the difference kernel (from `fit_gauss_width`), the position update,
the profile regressions with their spreads, the global transfer search, the residual
noise estimate. The halo term in the pseudo-atom branch of `atoms%convolve`, with its
test. The `density` sub-suite of `unit_single`, registered in the policy's table.
Unit tests on arrays, expected values from closed forms and from particles rendered by
the module itself with known parameters and added noise: the uniform-B particle
(`gamma = 0`) returns `gamma` within the equivalent of 1 A^2 from twelve to five
neighbours, `beta` within 0.005 per neighbour, `a` within 0.05, `B_h` within 20
percent; a particle with `gamma > 0` and `beta < 0` returns both within the same
tolerances; a fit without the halo term on the halo'd particle reproduces the known
artefact (the 4.4 A^2 of the species work within a factor of two), which documents why
the term exists. Exit: the sub-suite passes, the fast gate passes, `species_discovery`
and the `single_atoms_stats` entries are unchanged from phase 0.

**Phase 2, the fit on the present path.** `fit_density_model` in `nanoparticle` after
`discard_atoms`, on the level-1 atoms: the noise region as the species branch builds
it, with the residual between the atoms when the map has no empty region; the
provisional fit, the globals, the sweeps, the profiles; the products of section 5 with
the strong atoms only (`_ATMS.pdb` columns, `_SIM.mrc` from the per-atom table,
`_FIT.mrc`, `_RES.mrc`, `_atoms.csv`, `_model.txt`); `valid_corr` out of the B column.
`single_atoms_stats` and `species_discovery` read the new products; the high-level
check of the synthetic recovery through `detect_atoms` joins the `species_discovery`
test program as a case (no new CTest entry; the program keeps its name until the
species workflow takes it). Exit: the synthetic floors of phase 1 hold through
`detect_atoms`; the strong atoms on the 66 maps are those of phase 0's species-free
run (count and positions, or every difference named); on every map the residual
inside the envelope over `sigma_n` is below the equal-atom render's, and the
shell-mean intensity profile of `_FIT.mrc` follows the map's within the noise; the
`single_atoms_stats` floors hold.

**Phase 3, the weak atoms.** The residual recovery and its pruning on the present path,
after the provisional fit; the full fit over strong and weak atoms; the labels in
every product; `discover_species` reduced to the class call on the fitted
intensities of the table (the aperture sum kept in the table as a cross-check), its
`_species` files unchanged in meaning. `species_discovery` keeps every floor it has
today (weak recall, false atoms, labels) on the new path; the `species` sub-suite
adapted to the routines that moved. Exit: those floors, the fast gate, the
`single_atoms_stats` entries; on the pure synthetic case no weak atom beyond the
expected false count.

**Phase 4, the real maps and the report.** The validation of section 6 on the 66
maps with the phase-3 build, against phase 0: the residual over noise, the fade
profile and the wells, the strong atoms, the weak atoms described, the time-point
stability of the particle parameters; the detection-gap table, the second pass with
the fitted kernel crossed with the mask factor at 1.0, 1.1, 1.2 and 1.3 on all 66 maps;
FePt if its volume is on the machine (the noise from the residual between atoms,
`beta` against 5.2 percent per neighbour, the weak atoms against their traced peaks),
otherwise reported as not run. The second pass is a code change if the table says it
helps: then it is made in this phase, with the file-table row, and the table is
remeasured with it. The report, `atomic_density_model_report.md` beside this note,
reads on its own: what the model reproduces of the real maps, what it does not, the
gap table and what decided it, and what is left open. The status note of the
species work gets a pointer to the moved recovery.

## 11. File table

A phase edits only its files; a file the work turns out to need is added here, with
its phase and a one-line reason, before it is edited. Testers (`*_tester.f90`) of a
listed Fortran file are covered by its row.

| File | What changes | Phases |
| --- | --- | --- |
| `doc/implementation_notes/planned/atomic_density_model.md` | this note: Progress, file-table rows, rulings | each |
| `doc/implementation_notes/planned/atomic_density_model_report.md` | the review report | 4 |
| `doc/implementation_notes/planned/species_discovery_status.md` | pointer to the recovery's new home and the reduced `discover_species` | 4 |
| `doc/implementation_notes/planned/species_discovery_unfinished.md` | the same pointer | 4 |
| `src/main/nano/simple_nano_density.f90` | new: the model, its fit and rendering, with its tester; the fitted kernel as a detection template if the gap table calls for it (phase 4) | 1, 2, 3, 4 |
| `src/main/nano/simple_nano_species.f90` | the array routines the density module takes over leave; the class call on the fitted intensity | 1, 3 |
| `src/main/nano/simple_atoms.f90` | the halo term in the pseudo-atom branch of `convolve` | 1 |
| `src/defs/simple_defs_atoms.f90` | only if the pseudo-atom definitions need it | 1 |
| `src/main/nano/simple_nanoparticle.f90` | `fit_density_model`, the products (phase 2); the recovery and the weak atoms on the present path (phase 3); the second detection pass with the fitted kernel, kept or removed by the gap table (phase 4) | 2, 3, 4 |
| `src/main/nano/simple_nanoparticle_utils.f90` | only if a length-scale or contact routine needs the weak atoms (phases 2, 3); `phasecorr_one_atom` with the fitted kernel for the second pass (phase 4) | 2, 3, 4 |
| `src/main/commanders/simple/simple_commanders_atoms.f90` | `exec_detect_atoms`: the key rules once the recovery no longer needs `discover_species` | 2, 3 |
| `src/main/params/simple_parameters.f90` | the meaning of `discover_species` and `min_nbrs` on the present path | 3 |
| `src/main/params/simple_parameters_phases.f90` | their validation | 3 |
| `src/main/ui/single/single_ui_atom.f90` | the descriptions of `detect_atoms` and its keys | 2, 3 |
| `src/main/commanders/test/simple_commanders_test_single.f90` | `single_atoms_stats` and `species_discovery` on the new products; the synthetic case | 2, 3 |
| `src/main/commanders/test/simple_commanders_test_class.f90` | the `density` sub-suite of `unit_single` | 1, 3 |
| `src/main/ui/simple_test/simple_test_ui_class.f90` | the `unit_single` suite list | 1 |
| `src/main/ui/simple_test/simple_test_ui_highlevel.f90` | the description of `species_discovery` if it changes | 2, 3 |
| `production/CMakeLists.txt` | the comments and timeouts of the affected entries; the budget is unchanged | 2, 3 |
| `doc/policies/test_environment_policy.md` | the `unit_single` row names `density`; the descriptions of the affected high-level entries | 1, 2, 3 |

## 12. Progress

Each phase records here: the date, the files changed, the evidence with log paths
under the run directory, and each exit item of section 10 with how it was met.
