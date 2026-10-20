# Species discovery: open work

Status (2026-10-20): the opt-in prototype is done and in the tree; what is open is one
measurement run (phase 9), the validation that must follow it, promotion of species results
into the standard SINGLE products, and a list of open questions. The living entry point for
picking the work up is [`species_discovery_status.md`](species_discovery_status.md) (it holds the
current numbers, the FePt findings, the data locations and the lessons). The design, method,
rulings and the record of phases 0 to 8 are in
[`../completed/species_discovery_completed.md`](../completed/species_discovery_completed.md); the
run reports are in [`../completed/species_discovery_report.md`](../completed/species_discovery_report.md).

Terms. SINGLE is the cryo-EM-style workflow for nanoparticles. `detect_atoms` finds atoms in a
3D map. A "species" here is an intensity class of atom (HAADF-STEM contrast rises with atomic
number; the data can say there are two kinds of atom, not which elements they are). The
"coordination term" is the linear rise of an atom's intensity with its number of neighbours,
fitted together with the classes. A "recovered atom" is a weaker atom found in the residual map
(the map minus the model of the atoms already found). "Promotion" means making species results
part of the standard products instead of diagnostics beside them.

The learning so far. Species discovery does not invent species on single-species real maps
(one class on about 95 percent of 66 Pt maps), but it is insensitive: a class planted at 0.6 of
the Pt density was called on 4 of 66 maps, at 0.4 on 21. The limit is how much real atoms of one
element differ in intensity (about twice the measurement noise, mostly by environment). The
coordination term was added to absorb that trend. It has passed the synthetic tests and a numpy
check on one real map; it has not run on the real maps.

## 1. Phase 9: the coordination term on the 66 real Pt maps (not run)

A measurement only. No source change. It decides whether promotion goes ahead.

Preconditions, all met: the term is in the checkout (`fit_species_mixture` in
`src/main/nano/simple_nano_species.f90`, `discover_species` in `simple_nanoparticle.f90`,
`coord_slope` in `_species.txt`, unit test `test_coordination_term`); the maintainer built it and
the fast gate, `species_discovery` and the `single_atoms_stats` Pt entry passed (2026-10-12).
Not yet built after the later fixes of `extract_atom`, `atoms%copy` and the restored core label
floor (see section 4).

What to run, on the Dell, maps in `/data/NanoX/zenodo` (66 pairs, 0.358 Å):

- The five `detect_atoms` runs per map of the second validation: `element=Pt`; `Pt` with
  discovery; `Pt,Ni`; no element; no element with discovery.
- Runs 1 and 4 do not touch the class fit. Since 2026-10-15 the size rule of `discard_atoms`
  works after the contact-score pass (before, the pass renumbered atoms and the rule read the old
  numbering), so their numbers against the stored coordinates can differ from the second
  validation's, but only by atoms the size rule now removes. Name every difference. An atom moved
  or added stops the run.
- Runs 2, 3 and 5: number of classes `K`, admissibility, class separation `D`, class sizes and
  ratios, `coord_slope` per map, and per-coordination intensities after the term is removed. The
  mixture fits also stop earlier than in phases 7 and 8, so the comparison measures both changes
  together.
- The planting of phase 8 at 0.6, 0.4, 0.25 and the noise control, same seeds; recall, precision
  and false labels per map against phase 8's.
- Compare the `species_discovery` metrics with phase 7's copy (on the Dell, same test seeding).
- The run spec is `~/agent_runs/runs/species_validation_coord/` (PHASES "9"), with
  `LOCAL_REPO` = `~/src/SIMPLE`. Products stay in its `scratch/keep/validation/`. Rsync them back
  and make figures on the Mac (the Dell has no numpy).

Expectations, stated before the run (they are the acceptance checks):

1. One class in runs 2 and 5, and `Pt,Ni` inadmissible, on at least 95 percent of the unplanted
   maps.
2. The class planted at 0.4 called (the forced two-class fit admissible) with recall at least
   0.9 on at least 90 percent of the maps, with false labels on at most 2 percent of the
   unplanted atoms (the stricter "at most one" of phase 8 reported beside it).
3. The noise control never calls a second class.
4. 0.6 and 0.25 are reported, not gated.

Write the result in the run's section of the report. Four real maps accepted `Pt,Ni` in phase 7
(EFC4_time03s, ESC4_time01s, ESIP1_time08s, ESIP2_time05s), one more than the 5 percent allowed;
phase 9 shows whether the coordination term removes them.

## 2. Decision after phase 9

- **If the expectations hold:** do the validation of section 3, then promotion (section 5).
- **If they do not hold:** make occupancy explicit. The simple version calls species only on
  fully coordinated atoms (12 neighbours), where full occupancy can be assumed, and reports every
  other atom's intensity as an occupancy of its likeliest class. In the numpy check it was as
  sensitive as the coordination term on the atoms it labels (about half of them). The full
  version belongs to the B-factor and occupancy work (`atom_bfactor_model.md`).

## 3. Validation before promotion (not part of a run yet)

Three parts, judged by the maintainer:

1. Particles passed through a real reconstruction's low-pass filtering and soft mask, with
   background and shell-dependent noise (the synthetic fixtures do not reproduce the real
   environment trend; they show the opposite, a halo effect).
2. Label reproducibility across real half maps. The archive has none; the half-map route exists
   in `detect_atoms` (`vol_even`, `vol_odd`).
3. An A/B of `detect_atoms` with and without discovery on a real Pt/Ni map, at matched
   false-positive rate. The archive is all Pt. Real two-species data now come from atomic
   electron tomography (AET): the FePt reconstruction of Yang et al. (Nature 542, 75, 2017) in
   `~/Downloads/TomoFromPeter` (`for_simple/` has the volume as MRC and their models as PDB), and
   possibly the PtNi tomograms of Jeong et al. (Nat. Commun. 16, 8026, 2025), to be asked for.
   First step, needing no Dell and possible before phase 9: run `detect_atoms` with
   `element=Pt,Fe` on the FePt volume, with and without the coordination term, and compare with
   their model by coordination number. Findings so far (status note): Pt over Fe intensity 1.7,
   class separation `D` 2.6 (our gate `D >= 3` would reject the split), and the surface
   composition depends on how the environment is modelled (86 percent Fe with a plain fit, 20
   percent with the coordination term, 65 to 68 percent in their model). Species calls at low
   coordination are model-dependent and must be reported as such.

## 4. Small items still open

- Build and test the three same-day fixes not yet built: `atoms%extract_atom` copies every
  per-atom field (it dropped the occupancy); `atoms%copy` copies the charge column; the core
  case's label floor restored in `species_discovery`.
- The three CdSe `single_atoms_stats` entries (`CdSeW`, `CdSeZ`, `CdSeR`) were not run by hand
  after the coordination term; the next high-level ctest covers them.
- Recovered atoms off the stored model: about 20 per real map sit one bond length from a stored
  atom, mostly in the outer half. They may be partly occupied sites or noise. Real half maps
  would tell. Expected noise maxima above the second-stage threshold are 2 to 4 per map.
- The present detection is short on ten of 66 maps. ESC1_time28s depends on the mask: the
  automatic 40.1 Å (1.2 times the measured diameter) matches 85 percent of its stored atoms,
  37.8 Å matches 99.8 percent. The factor 1.2 may need a look for non-spherical particles; on two
  other maps the measured diameter is 19 and 27 percent above the coordinate extent.
- `single_workflow_wurtzite` crashed in `refine3D` on 2026-10-09 (`simple_strategy3D_srch.f90`,
  `subspace_inds` past its bound), also at the base revision. Not species work; not rechecked.
- Production `seed_rnd` still honours the `SIMPLE_SEED` environment variable, against policy. The
  plan for an `rnd_seed` key is `../../refactoring_notes/planned/rnd_seed_parameter.md`.
- A light species confined to the surface layer is not a target (ruling of 2026-10-09): it cannot
  be told from weak, partly occupied surface sites.
- B factors and occupancy follow after species: `atom_bfactor_model.md`. Per-atom Gaussian fits
  on filtered maps invent a core-to-surface B rise (+4.4 Å² on a particle generated with one B).
  The preferred route measures a smooth B profile with depth in `refine3D_nano`, where references
  are matched to the particle images. An occupancy profile from that work could replace the
  coordination slope.

## 5. Promotion (waits for sections 1 to 3)

What promotion contains, from the design:

- Recovered atoms into `_ATMS.pdb`, with data-derived components in `_CC.mrc`.
- Class intensities and per-atom widths into `_SIM.mrc`.
- `element` optional in the other four SINGLE programs.
- `autorefine3D_nano` end to end on a simulated two-species trajectory, then on a real one.
- `atoms_stats` without the lattice table.

Until then the five existing products (`_ATMS.pdb`, `_CC.mrc`, `_SIM.mrc` among them) must stay
identical with and without discovery. The validation of section 3 comes before the B-factor
work. Two references in the design need amending when this lands: the single_pcg_integration
note's rule that atoms with amplitude below 0.2 of the median are removal candidates (it would
flag a whole light species; reference it to the noise instead), and its statement that trailing
reconstruction blends the simulated map away (it does not match the accumulator-domain contract).

## 6. Open questions (all still open; the prototype took the recommended defaults)

- Is 3-6:1 (Pt:Ni) the atom ratio or the intensity ratio? What is the signal-to-noise of the
  weakest atoms you can see in a current map? Together they say whether the light species is in
  the regime of the easy cases or of the hard ones.
- Class symbols in the PDB element column (`X1`, `X2`, ...) with the class also in the chain, or
  one symbol with the class in the chain only?
- Radius from the centre, or depth below the surface, as the radial coordinate of the profiles?
- Starting volume without an element: `solve3D_nano`, or the ideal-lattice start with a
  user-supplied lattice constant?
- Is the close-packed neighbour cutoff acceptable as a default until the pair-distance form is
  in?
- Stage 2 of the fit assumes equal integrals within a class. Acceptable for surface atoms, or
  report the stage-1 width where the two disagree?
- The target of 0.2 expected false atoms per particle behind the calibrated thresholds, and the
  neighbour gate (`min_nbrs`) on by default: both are prototype defaults, to be judged on the
  Pt/Ni A/B of section 3.

## 7. Bindings

- Bound rulings: phase 6 (fitting with a map-filtered kernel) is withdrawn; do not try it again
  (it worsened the intensity ratios). The contact-score pruning removes weak, poorly coordinated
  surface atoms on purpose, so recovered atoms are pruned against the present detection's atoms
  only. Do not add fitted machinery to fix a secondary output. Count mixture fits before wrapping
  a search around them (the first coordination term timed out `unit_single`).
- Exclusions of the design (section 12 of the completed note): no element names from the data,
  no spatial priors on labels, no lattice model, no anisotropic widths, no per-atom amplitudes in
  `_SIM.mrc`, no change to the present products or pruning before the validation.
- Compare synthetic fixtures with a real map before blaming code for a synthetic failure.
