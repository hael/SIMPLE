Refine3D Iteration
    ↓
Produces partition volumes
    ↓
Calls volassemble
    ├─ Assemble volumes from partitions
    ├─ Generate even/odd pairs
    ├─ Execute automasking (if automsk != 'no' AND frequency met)
    ├─ Apply non-uniform filtering (if enabled)
    └─ Return volumes + mask file
         ↓
      Mask used in:
      • FSC calculation
      • Centering logic
      • Reference filtering

