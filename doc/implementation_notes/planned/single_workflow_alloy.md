# Random-site Pt/Ni SINGLE workflow

## Contract

Add an independently registered `single_workflow_alloy` high-level gate. Its truth is the established
20 A FCC Pt cluster with the 71 Ni sites from the first validated fixed-seed draw frozen as a fixture;
the other 214 sites are Pt. The alloy follows the hand-staged atomic-prior workflow used by the wurtzite suite;
Pt remains the FCC lattice selector used to initialize `analysis2D_nano`. Reconstruction does not
receive `element=Pt,Ni`: atom detection and the PCG atomic prior remain species-free and use fitted
positions, widths and amplitudes.

The gate checks that the written truth PDB preserves all 285 planted site elements and their 71/214 counts, that the atomic-prior
products exist, and that the existing SINGLE map and pose contracts hold in cubic symmetry.

## Implementation

- Generate the normal FCC Pt lattice and PDB with `simulate_nanoparticle`.
- Apply the frozen 71-site Ni mask to the canonical 285-site FCC PDB and rerender it with the
  element-specific atom kernels. This keeps the original scattered alloy layout independent of the
  workflow RNG stream.
- Add the `alloy` suite to the test UI, registry and CTest registration, with its own working directory.
- Raise the fixed process budget for the explicitly requested independent workflow entry and update
  the policy inventory and card.
- Add the suite to `SIMPLE_data_testing/submit_data_testing.sh` and give it a dedicated report label;
  that existing harness supplies the independent campaign run, metric collection and 5 A report gate.

## Validation

- Static registry, CLI, descriptor, FLEX-DAG and diff checks passed.
- The incremental build and install passed. The fast gate passed 14 of 15 entries; its only failure
  was the pre-existing four-thread HHCNN bit-for-bit reduction check in `unit_numerics`, reproduced
  in an isolated rerun (1,925 of 1,926 assertions passed).
- `simple_test_exec test=single_workflow suite=list` listed `fcc`, `alloy` and `wurtzite`.
- The complete `single_workflow_alloy` run, launched in the environment used by
  `SIMPLE_data_testing/submit_data_testing.sh`, passed all of its internal gates in 966.7 s at one
  thread. The planted PDB contained 285 atoms: exactly 71 Ni and 214 Pt sites. Its direct and
  mirrored docking correlations were 0.998925 and 0.998740, the selected soft-mask map correlation
  was 0.978985, and its truth resolutions were 6.36 A at FSC=0.5 and 3.81867 A at FSC=0.143.
- `SIMPLE_data_testing/make_report.sh` completed on the run after two BSD Bash/awk portability
  corrections in that repository. Its formal report criterion is FSC=0.5 at no worse than 5.0 A,
  so the report card is **FAIL** at 6.36 A even though the workflow's internal FSC=0.143 ceiling
  passes at 3.82 A. A 20-seed profile remains required before changing an inherited workflow floor.
- After the alloy mask was frozen, the complete workflow passed again in 925.2 s. The generated
  285-site PDB was byte-identical to the first validated run (SHA-256
  `a2c5485617971a9376485fe3c032d367e0085f49f8fd6e2176862c28cb1acb38`), all three fixture checks
  passed, and the map, FSC, docking and pose metrics were unchanged.
