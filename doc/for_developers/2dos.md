# TDO list
**Date:** February 22, 2026  
**Status:** Incomplete

### Encapsulation issues
- Remove mixes of public and private derived type declarations and do proper encapsulation
- Specifically: qsys_env

### Quality issues
- Replace all polarizations and reprojections with oversampled ones, remove non-oversampled routines

### Execution issues
- Look into dedicated workflows for shared-memory excution that overcomes the need for disk interaction, see simple_commanders_refine3D_inmem for example code
- The code paths for shared-memory and distributed execution are too divergent, we neeed to find ways of reducing duplicated logic unless duplication provides some real value, like avoiding disk interactions
