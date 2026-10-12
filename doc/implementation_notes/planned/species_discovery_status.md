# Species discovery: where it stands and where to go next

Status note, 2026-10-12, at master `618d709ec`. It is the entry point for
picking the work up again after a break and does not repeat the design.
Three other documents hold the detail:

- `species_discovery.md` (the note): the design, every ruling, and the
  record of each phase.
- `species_discovery_report.md` (the report): the reports of the
  unattended runs.
- `atom_bfactor_model.md`: the follow-on work on per-atom B factors and
  occupancy.

## What exists

`detect_atoms` (SINGLE) has an opt-in species discovery. A run without the
new keys behaves as before.

- `element` may be left out. Atoms are then detected with a Gaussian
  template instead of an element kernel.
- `discover_species=yes` adds two steps after the present detection.
  - First, weaker atoms are looked for in the residual map (the map minus
    the model of the atoms already found). They are kept only if they pass
    the same contact-score pruning as the present detection. These are the
    recovered atoms.
  - Second, every atom's aperture intensity is measured: the density summed
    in a sphere of half the nearest-neighbour distance. The atoms are then
    sorted into one to three intensity classes, and each class is a
    candidate species. `nspecies` fixes the number of classes; 0 lets the
    data choose.
- `element=Pt,Ni` gives the species as a list, brightest first. The number
  of classes is then fixed, and the names go into the products.
- A two-class result is admissible when both classes are large enough and
  their means are at least three pooled spreads apart (`D >= 3`).
- Since `618d709ec` the class fit has one more term. Intensity rises
  linearly with the atom's coordination number (its number of neighbours),
  with one slope shared by all classes. The slope is fitted with the classes
  and written to `_species.txt` as `coord_slope`.
- Discovery writes four extra files: `_species.csv` (one row per atom),
  `_species_radial.csv`, `_species.txt` (the summary) and `_species.pdb`
  (class in the element column). They are diagnostics only: the five
  existing products, `_ATMS.pdb`, `_CC.mrc` and `_SIM.mrc` among them, are
  identical with and without discovery.

The same work also produced these changes:

- B factors in the simulator's element kernels (`pdb_bfac=yes` on
  `simulate_nanoparticle`).
- Compound-crystal lattices for `CdSeW`, `CdSeZ` and `CdSeR`.
- A fix of the automatic mask of `detect_atoms`, which covered the whole box
  on unmasked maps.
- The test seed fixed inside `simple_test_exec`.

The commits:

- `d3480ae79`: phases 0 to 3, the first pass.
- `29808bf3e`: phases 4 and 5, species lists, element kernels, compound
  lattices.
- `0029ed77e`: aperture intensity at one common width; the test seed.
- `c28bfc101`: the mask fix.
- `618d709ec`: the coordination term; the core case's label check back;
  `extract_atom` copies every field.

Code:

- `src/main/nano/simple_nano_species.f90`: the class fit,
  `fit_species_mixture`, on the clustering library's extreme-deconvolution
  mixture `simple_xd_gmm`.
- `src/main/nano/simple_nanoparticle.f90`: `discover_species`, which runs
  the recovery, the intensities and the products.
- `src/main/nano/simple_atoms.f90`: pseudo-atoms, species lists, lattices.

Tests:

- the `species` sub-suite of `unit_single`;
- the high-level `species_discovery`;
- four `single_atoms_stats` entries, one for Pt and one for each CdSe
  selector.

## What the tests and the real maps show

**Synthetic particles** (`species_discovery`). Each case is a 285-atom
particle, rendered with element kernels, filtered as a reconstruction is,
with noise shaped like a real map's. The four cases are pure Pt, Pt/Ni mixed
at random, a Ni core under a Pt skin, and Pt with a light species (Al).

On 2026-10-12 every check passed:

- no false atoms;
- recall of 0.986 to 1.000 for both species;
- every atom labelled with its element in all three mixed cases;
- a two-class fit forced on pure Pt not admissible (`D` 0.55).

The fitted intensity ratio of the classes is biased low: by 23% for Pt/Ni
mixed, 6% for the Ni core and 43% for Pt/Al. The classes still come out in
the right order. This bias comes from measuring on filtered maps and is a
documented limitation. Phase 6 tried to remove it by fitting with a
map-filtered kernel; that made things worse and was withdrawn on 2026-10-09.
Don't try that form again.

**Real maps.** The archive holds 66 Pt maps: 15 particles, each with 2 to 8
time points, all one species. The second validation run (2026-10-11, before
the coordination term) found:

- **False species are rare.** One class was found on 95.5% of the maps with
  `element=Pt discover_species=yes`, and on 97.0% with discovery but no
  element.
  `element=Pt,Ni` was rejected as a two-class fit on 62 of 66 maps.
- **The present detection** matched the stored coordinates on 56 of 66 maps
  with `element=Pt` and on 47 without an element. A map matches when 98% of
  its stored atoms are found within 1 Å. The ten maps that fall short match
  85% to 98%.
- **Sensitivity is poor.** A weaker class was planted into 20% of the atoms
  at 0.6 of the Pt density, about Ni. It was called a second class on 4
  maps; at 0.4 it was called on 21. The planted atoms were detected; they did
  not separate as a class.
- **The reason is the spread of real intensities.** Real Pt atoms differ in
  aperture intensity by about twice the measurement noise. Most of that
  follows the atom's environment:
  - +3.8% per neighbour at the median map, 26% more at 12 neighbours than
    at 5;
  - the fitted width flat.

  Missing intensity without broadening is what partial occupancy of the
  outer sites gives. The synthetic particles show the opposite trend, a halo
  effect of the filtering, so they do not reproduce this.

**The coordination term** was chosen on 2026-10-12 to absorb that trend.
The numpy check used one real map (`Pt_example`, 476 atoms), with classes
planted as above and eight seeds per level:

- the 0.4 class was found on 8 of 8 seeds (recall 0.96, about 3 wrong
  labels);
- the 0.6 class was found on 5 of 8 (recall 0.72);
- no second class was called on the unplanted map or the noise control.

On the synthetic particles it changed no label. In SIMPLE it has so far
passed only the synthetic tests. **It has not run on the real maps.**

## Next steps, in order

1. **Phase 9 on dell: the coordination term on the 66 maps.** This is a
   measurement only, with no code change. The plan is in the note, section
   10.11.
   - The run spec is `~/agent_runs/runs/species_validation_coord/`
     (PHASES "9"). It reuses the second validation run's scripts, with
     `LOCAL_REPO` = `~/src/SIMPLE` (now at or after `618d709ec`).
   - It repeats the five `detect_atoms` runs per map: `Pt`, `Pt` with
     discovery, `Pt,Ni`, no element, and no element with discovery. It also
     repeats the planting at 0.6, 0.4 and 0.25 and the noise control, with
     the same seeds.
   - The expectations were stated before the run:
     - one class, and `Pt,Ni` rejected, on at least 95% of the unplanted
       maps;
     - the 0.4 class found with recall of at least 0.9 on at least 90% of
       the maps, with false labels on at most 2% of the unplanted atoms;
     - the noise control never called a class.
   - The products stay in the run's `scratch/keep/validation/` on dell.
     Rsync them back and make the figures on the Mac; dell has no numpy.
   - The run also compares the `species_discovery` metrics with phase 7's
     copy. That copy is on dell and used the same test seeding.
2. **Decide on phase 9's result.**
   - **If the expectations hold,** the validation before promotion comes
     next (note, section 10.8). It has three parts:
     - particles passed through a real reconstruction's filtering and
       masking;
     - label reproducibility across real half maps;
     - an A/B of `detect_atoms` with and without discovery on a real Pt/Ni
       map, at matched false-positive rate.

     **No real two-species map is at hand**, because the archive is all Pt.
     Getting one is what blocks promotion. Promotion itself comes after
     that: recovered atoms into `_ATMS.pdb`, classes into `_SIM.mrc`,
     `element` optional in the other SINGLE programs, and
     `autorefine3D_nano` on a two-species trajectory.
   - **If the expectations don't hold,** make occupancy explicit instead.
     The simple version calls species only on fully coordinated atoms (12
     neighbours), where full occupancy can be assumed. It reports every
     other atom's intensity as an occupancy of its likeliest class. In the
     numpy check it was as sensitive as the coordination term on the atoms
     it labels, about half of them. The full version belongs to the
     B-factor and occupancy work.
3. **B factors and occupancy** (`atom_bfactor_model.md`), after species.
   The fitting method is open. Per-atom Gaussian fits on filtered maps
   invent a core-to-surface B rise: +4.4 Å² on a particle generated with
   one B. The preferred route measures a smooth B profile with depth in
   `refine3D_nano`. There references are matched to the particle images, so
   no map filtering is involved. An occupancy profile from that work could
   take over the role of the coordination slope.

## Smaller open items

- **Recovered atoms off the stored model.** About 20 per real map sit one
  bond length from a stored atom, mostly in the outer half of the particle.
  They may be partly occupied sites or noise. Real half maps would tell,
  and the archive has none. On these boxes the expected number of noise
  maxima above the second-stage threshold is 2 to 4 per map. The
  first-stage threshold meets the target of 0.2.
- **The present detection is short on ten maps.** ESC1_time28s depends on
  the mask. The automatic 40.1 Å (1.2 times the measured particle diameter)
  matches 85% of its stored atoms; 37.8 Å matched 99.8%. The factor 1.2 may
  deserve a look on non-spherical particles. On two other maps the measured
  diameter is 19% and 27% above the coordinate extent.
- **`Pt,Ni` was accepted on four real maps**: EFC4_time03s, ESC4_time01s,
  ESIP1_time08s and ESIP2_time05s. That is one more than the 5% allowed.
  Phase 9 shows whether the coordination term removes them.
- **A light species confined to the surface layer is not a target**
  (ruling of 2026-10-09). It cannot be told from weak, partly occupied
  surface sites.
- **The three CdSe `single_atoms_stats` entries** were not run by hand after
  the coordination term. The next high-level ctest covers them.
- **`single_workflow_wurtzite` crashed in `refine3D`** on 2026-10-09
  (`simple_strategy3D_srch.f90`, `subspace_inds` past its bound). It crashed
  at the base revision too. It is not species work and has not been
  rechecked since.
- **The random seed.** Production `seed_rnd` still honours the `SIMPLE_SEED`
  environment variable, against policy. The plan for an `rnd_seed` key is
  `doc/refactoring_notes/planned/rnd_seed_parameter.md`. It awaits review
  and needs no dell run.
- **The open questions of the note's section 13** are still open:
  - the atom ratio against the intensity ratio of the target alloys, and
    the signal-to-noise of the weakest visible atoms;
  - the class symbols in the PDB;
  - radius against depth as the radial coordinate;
  - the starting volume without an element;
  - the close-packed neighbour cutoff;
  - equal integrals within a class at the surface;
  - the false-atom target and the default of the neighbour gate
    (`min_nbrs`).

  The prototype took the recommended defaults. The Pt/Ni A/B is where they
  get judged.

## Lessons worth keeping

- **Compare the synthetic input with a real map before blaming the code on
  a synthetic failure.** Phase 5 stopped when the present detection found 3
  of 285 atoms on a fixture, although it works on real maps. The fixture was
  the defect: its B factors, surface-to-core peak ratio and noise spectrum
  were off, and it lacked the negative wells between atoms.
- **The contact-score pruning removes real densities on purpose.** These
  are weak, poorly coordinated, high-B surface atoms. Discovery must not put
  them back, so recovered atoms are pruned against the present detection's
  atoms only.
- **Don't add fitted machinery to fix a secondary output.** Phase 6 did, and
  made the intensity ratios worse.
- **Count the mixture fits before wrapping a search around them.** Each
  tight fit can run hundreds of iterations, and the first version of the
  coordination term timed out `unit_single`.

## Where things are

- **The maps:** `/data/NanoX/zenodo` on dell, read-only. It holds 66 pairs,
  `mrcs/<name>.mrc` and `pdbs/<name>.pdb`, at 0.358 Å, in boxes of 120, 160
  and 200. The `._*` files are macOS debris.
- **One real map with its coordinates on the Mac:**
  `~/agent_runs/runs/species_element_model/Pt_example/` (476 atoms). Its
  `_even`/`_odd` maps are not half reconstructions.
- **The numpy emulation and checks:**
  `~/agent_runs/runs/species_validation/numpy_checks/`:
  - the coordination-term checks, `real_fit.py`, `classfit.py` and
    `fixture_check.py`;
  - the mixture-fit cost model, `em_cost*.py`;
  - the pipeline emulation, `emul_species_cases.py`.
- **The run specs** are under `~/agent_runs/runs/`:
  - `single_species_proto`: phases 0 to 3;
  - `species_element_model`: phases 4 to 6, with `rulings.md`;
  - `species_validation`: the first validation, which found the mask bug;
  - `species_validation_r2`: phases 7 and 8;
  - `species_validation_coord`: phase 9, not yet run.
