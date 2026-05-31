# Abinitio3D Cavgs Policy

This document records the current policy for `abinitio3D_cavgs`, the
class-average route for generating an initial 3D model. It is separate from
the particle-based [abinitio3D_policy.md](abinitio3D_policy.md) policy and
from the base [refine3D_policy.md](refine3D_policy.md) policy.

## 1. Scope

`abinitio3D_cavgs` builds an initial 3D model from selected 2D or 3D class
averages. It creates a temporary project containing even and odd class-average
entries, runs staged `refine3D` over that temporary project, maps the selected
class orientations back to the original project, and optionally writes
validation artifacts.

Only refinement and reconstruction children are distributed. The
`abinitio3D_cavgs` master itself rejects direct worker execution with `part`.

## 2. Defaults

The route sets:

- `sigma_est=global`
- initial `oritype=out`, then `ptcl3D` for the temporary project
- `bfac=0`
- `filt_mode=none`
- `automsk=no`
- `nu_refine=no`

When unset, it supplies:

- `mkdir=yes`
- `objfun=euclid`
- `overlap=0.95`
- `prob_athres=90`
- `cenlp` from the ab initio controller default
- `imgkind=cavg`
- `noise_norm=no`
- `lpstart=20`
- `lpstop=8`
- `gauref=yes`

NU filtering and automasking are intentionally disabled for this route.

## 3. Temporary Project

The workflow writes a temporary project named
`abinitio3D_cavgs_tmpproj.simple`.

The temporary project:

- preserves project info, compute environment, and job-process metadata
- registers the even and odd class-average stacks as two stacks
- expands each class average into an even entry and an odd entry
- copies class/state/orientation metadata into both entries
- sets even/odd flags and stack indices for the temporary `ptcl3D` segment

The temporary project is deleted at the end of the workflow.

## 4. Inputs and Low-Pass Schedule

`imgkind=cavg` reads state labels from `cls2D`. `imgkind=cavg3D` is accepted
as an input mode, but the current implementation still retrieves the selected
class-average stack and then uses class state information for the temporary
entries.

The class-average stack plus `_even` and `_odd` companion stacks must exist.
If all class-average states are zero, the workflow stops.

The number of particles in the temporary project is twice the number of class
averages. If distributed execution requested more partitions than temporary
entries, `nparts` is reduced to the number of even/odd class-average entries.

Stage low-pass limits are derived from class FRCs by default. Explicit
`lpstart_ini3D` and `lpstop_ini3D` override that schedule together; supplying
only one is an error.

## 5. Staged Refinement

The number of ini3D stages is capped by `abinitio_nstages_ini3D_max()`. A user
`nstages` value can shorten the route up to that cap.

Before staged refinement, `rndstart` randomizes orientations, zeros shifts,
randomizes states uniformly for multi-state runs, and reconstructs starting
volumes. Starting volumes and half maps are renamed to the standard
`refine3D` start-volume names, including `_unfil` copies of the half maps.

Each stage is configured through the shared ab initio stage controller with
`l_cavgs=.true.`. In this mode:

- `envfsc=no`
- `filt_mode=none`
- `automsk=no`
- `update_frac` is deleted from emitted refine3D child commands
- `snr_noise_reg` is set from the stage controls
- early Gaussian reference filtering remains available through stage policy

At the symmetry-search stage, the workflow runs the shared symmetry handling
used by ab initio workflows.

## 6. Mapping Back

After staged refinement, the temporary `ptcl3D` and `out` segments are read
back. Temporary stack-index and even/odd fields are removed.

For each original class, the even or odd temporary entry with the better
correlation is selected. Its correlation, projection index, Euler angles,
2D shift, and state are copied into the original project's `cls3D` segment.

The resulting `cls3D` orientations are mapped back to particles with
`map2ptcls`.

## 7. Validation and Outputs

When the route runs the maximum ini3D stage count, it produces validation
artifacts:

- final original-sampling reconstructions without automatic postprocess
- final raw and low-pass diagnostic volumes
- `vol_cavg` entries in `os_out` for final class-average volumes
- `final_oris.txt`
- reprojections of the final low-pass volumes
- an alternating `cavgs_reprojs.mrc` stack
- a shifted class-average stack registered as `cavg_shifted`
- optional class-average ranking by cavg-vs-reprojection correlation

Even/odd convergence diagnostics report average angular distance. Multi-state
runs also report even/odd state overlap.
