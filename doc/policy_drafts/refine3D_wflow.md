# Refine3D Workflow Draft

Superseded by:

- `doc/policies/refine3D_policy.md`
- `doc/policies/automasking_policy.md`
- `doc/policies/nonuniform_filtering_policy.md`

Current high-level flow:

1. `refine3D` produces partial reconstructions.
2. `volassemble` assembles state volumes and even/odd pairs.
3. `volassemble` runs the shared volume-domain postprocessing policy:
   - automasking when enabled
   - nonuniform filtering when enabled
4. The reconstructor consumes compatible state masks for FSC when present; otherwise it falls back to the spherical mask.

The important architectural change is that shared-memory volume work is centralized in `volassemble` for both shared-memory and distributed refine3D paths.
