# TDO list
**Date:** February 22, 2026  
**Status:** Incomplete

## Preparation for block-tree 2D and 3D refinement applications
- Remove global pftc instances, ensure pft instantiatable, and consider pftc as part of builder to reduce dummy argument passing

### Global Dependency Issues
- Remove eulc_sigma2_glob pointer and replace with appropriate alternative

### Encapsulation issues
- Remove mixes of public and private derived type declarations and do proper encapsulation
- Specifically: euclid_sigma2, qsys_env

### Quality issues
- Replace all polarizations and reprojections with oversampled ones, remove non-oversampled routines

### Execution issues
- Look into dedicated workflows for shared-memory excution that overcomes the need for disk interaction, see simple_commanders_refine3D_inmem for example code
- The code paths for shared-memory and distributed execution are too divergent, we neeed to find ways of reducing duplicated logic unlless duplication provides some real value, like avoiding disk interactions
