# UI Placeholder Review

This document inventories placeholder values by actual parameter key and program context.
It is generated review evidence from the registered Fortran UI, not a source of runtime metadata.
After review, apply approved wording to the owning Fortran descriptor and regenerate this inventory.

Registry source: `/Users/elmlundho/src/SIMPLE/build/production/simple_private_exec prg=print_ui_json`.

Each heading is one parameter key. Each row is a separate `ui_program_input` instance, so a key may have context-dependent wording.
Edit **Placeholder** directly. Leave choice, binary, and hidden inputs empty unless the Fortran descriptor rule changes.

Inventory: 343 unique parameter keys and 1382 parameter instances.

## Parameter `acf`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| stackops | parameter input/output | binary | Autocorrelation, A * conjg(A) | *(empty)* |

## Parameter `algorithm`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | search controls | multi | Algorithm for motion correction | *(empty)* |
| preproc | parameter input/output | multi | Algorithm for motion correction | *(empty)* |
| preprocess | search controls | multi | Algorithm for motion correction | *(empty)* |
| tseries_motion_correct | search controls | multi | Algorithm for motion correction | *(empty)* |

## Parameter `amsklp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| auto_spher_mask | filter controls | num | Low-pass limit for envelope mask generation | e.g. 10 |
| automask | filter controls | num | Low-pass limit for envelope mask generation | e.g. 10 |
| automask2D | filter controls | num | Low-pass limit for envelope mask generation | e.g. 10 |
| estimate_diam | filter controls | num | Automask low-pass limit | e.g. 10 |
| gen_pickrefs | filter controls | num | Automask low-pass limit | e.g. 10 |
| refine3D | filter controls | num | Low-pass limit for envelope mask generation | e.g. 10 |
| refine3D_auto | filter controls | num | Low-pass limit for envelope mask generation | e.g. 10 |

## Parameter `angastunit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | parameter input/output | multi | Angle of astigmatism unit | *(empty)* |

## Parameter `angerr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_oris | parameter input/output | num | Rotation angle error half-width | e.g. 10 |
| oriops | parameter input/output | num | Rotation angle error half-width | e.g. 10 |

## Parameter `angstep`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| stackops | parameter input/output | num | Angular stepsize | e.g. 10 |

## Parameter `append`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| selection | parameter input/output | binary | Append selection to existing | *(empty)* |

## Parameter `astigerr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_particles | parameter input/output | num | Astigmatism error | e.g. 10 |

## Parameter `astigthreshold`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick_extract | filter controls | num | Astigmatism rejection threshold | e.g. 10 |

## Parameter `astigtol`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctf_estimate | search controls | num | Expected astigmatism | e.g. 10 |
| preprocess | search controls | num | Expected astigmatism | e.g. 10 |

## Parameter `athres`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fsc_area_score | parameter input/output | num | Cone half-angle | e.g. 10 |
| trajectory_make_projavgs | parameter input/output | num | Angular threshold (degrees) | e.g. 10 |

## Parameter `automsk`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | mask controls | multi | Perform envelope masking | *(empty)* |
| automask | mask controls | multi | Perform envelope masking | *(empty)* |
| estimate_diam | mask controls | multi | Perform envelope masking | *(empty)* |
| fsc_area_score | mask controls | multi | Perform envelope masking | *(empty)* |
| nu_filt3D | mask controls | multi | Perform envelope masking | *(empty)* |
| refine3D | mask controls | multi | Perform envelope masking | *(empty)* |
| refine3D_auto | mask controls | multi | Perform envelope masking | *(empty)* |
| refine3D_multi | mask controls | multi | Perform envelope masking | *(empty)* |

## Parameter `autoscale`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | binary | Automatic down-scaling | *(empty)* |
| refine3D_auto | search controls | binary | Automatic down-scaling | *(empty)* |
| refine3D_multi | search controls | binary | Automatic down-scaling | *(empty)* |

## Parameter `avg`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| stackops | parameter input/output | binary | Average stack | *(empty)* |

## Parameter `backgr_subtr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract | parameter input/output | binary | Perform micrograph background subtraction(new picker only) | *(empty)* |
| pick | search controls | binary | Perform micrograph background subtraction(new picker only) | *(empty)* |
| pick_extract | parameter input/output | binary | Perform micrograph background subtraction(new picker only) | *(empty)* |
| reextract | parameter input/output | binary | Perform micrograph background subtraction(new picker only) | *(empty)* |

## Parameter `balance`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| selection | parameter input/output | binary | Balanced selection of particles across classes | *(empty)* |

## Parameter `bandwidth_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | multi | Diffusion-map bandwidth mode | *(empty)* |
| denoise_project | filter controls | multi | Diffusion-map bandwidth mode | *(empty)* |
| flex_analysis | filter controls | multi | Diffusion-map bandwidth mode | *(empty)* |
| ppca_denoise | filter controls | multi | Diffusion-map bandwidth mode | *(empty)* |

## Parameter `bandwidth_tune`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | num | Ferguson bandwidth multiplier (default 1) | e.g. 10 |
| denoise_project | filter controls | num | Ferguson bandwidth multiplier (default 1) | e.g. 10 |
| flex_analysis | filter controls | num | Ferguson bandwidth multiplier (default 1) | e.g. 10 |
| ppca_denoise | filter controls | num | Ferguson bandwidth multiplier (default 1) | e.g. 10 |

## Parameter `beamtilt`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| assign_optics_groups | parameter input/output | binary | Use beamtilts in optics group assignment | *(empty)* |
| preproc | search controls | binary | Use beamtilts in optics group assignment | *(empty)* |

## Parameter `bfac`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctfops | filter controls | num | CTF B-factor | e.g. 10 |
| filter | filter controls | num | B-factor of Gaussian low-/high-pass filter | e.g. 10 |
| motion_correct | search controls | num | B-factor applied to frames | e.g. 10 |
| postprocess | filter controls | num | B-factor for sharpening | e.g. 10 |
| preproc | filter controls | num | B-factor for sharpening | e.g. 10 |
| preprocess | search controls | num | B-factor applied to frames | e.g. 10 |
| simulate_movie | filter controls | num | CTF B-factor | e.g. 10 |
| simulate_particles | filter controls | num | CTF B-factor | e.g. 10 |
| tseries_make_pickavg | search controls | num | B-factor applied to frames | e.g. 10 |
| tseries_motion_correct | search controls | num | B-factor applied to frames | e.g. 10 |
| volops | filter controls | num | B-factor for sharpening | e.g. 10 |

## Parameter `bfacerr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_particles | filter controls | num | B-factor error | e.g. 10 |

## Parameter `binwidth`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask | mask controls | num | Envelope binary layers width | e.g. 10 |

## Parameter `box`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| crys_score | parameter input/output | num | Particle box size | e.g. 10 |
| extract | parameter input/output | num | Particle box size | e.g. 10 |
| noisevol | parameter input/output | num | Particle box size | e.g. 10 |
| print_dose_weights | parameter input/output | num | Particle box size | e.g. 10 |
| print_fsc | parameter input/output | num | Particle box size | e.g. 10 |
| print_magic_boxes | parameter input/output | num | Particle box size | e.g. 10 |
| reextract | parameter input/output | num | Particle box size | e.g. 10 |
| simulate_nanoparticle | parameter input/output | num | Particle box size | e.g. 10 |
| simulate_noise | parameter input/output | num | Particle box size | e.g. 10 |

## Parameter `box_coarse`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Box size for coarse sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Coarse-pass box size | e.g. 10 |

## Parameter `box_extract`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | parameter input/output | int | Force box size (px, optional) | e.g. 10 |
| pick_extract | parameter input/output | num | Extracted particle image size | e.g. 10 |
| tseries_extractor | parameter input/output | num | Extracted particle image size | e.g. 10 |

## Parameter `box_fine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Box size for fine sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Fine-pass box size | e.g. 10 |

## Parameter `boxes`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| print_project_field | parameter input/output | binary | output coordinates in JSON format | *(empty)* |

## Parameter `boxfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| track_particles | file input/output | file | List of particle coordinates | e.g. coords.box |
| tseries_motion_correct | file input/output | file | List of particle coordinates | e.g. coords.box |

## Parameter `boxtab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_boxes | file input/output | file | List of box files | e.g. boxes.txt |
| import_movies | file input/output | file | List of box files | e.g. boxes.txt |

## Parameter `cavg_ini`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | binary | 3D initialization on class averages | *(empty)* |

## Parameter `cavg_ini_ext`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | binary | External class-average 3D initialization | *(empty)* |

## Parameter `cenlp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | num | Centering low-pass limit | e.g. 10 |
| abinitio2D_chunks | filter controls | num | Centering low-pass limit | e.g. 10 |
| abinitio3D | filter controls | num | Centering low-pass limit | e.g. 10 |
| abinitio3D_cavgs | filter controls | num | Centering low-pass limit | e.g. 10 |
| abinitio3D_nano | filter controls | num | Centering low-pass limit | e.g. 10 |
| autorefine3D_nano | filter controls | num | Centering low-pass limit | e.g. 10 |
| center | filter controls | num | Centering low-pass limit | e.g. 10 |
| cluster2D_nano | filter controls | num | Centering low-pass limit | e.g. 10 |
| refine3D | filter controls | num | Centering low-pass limit | e.g. 10 |
| refine3D_nano | filter controls | num | Centering low-pass limit | e.g. 10 |
| symaxis_search | filter controls | num | Centering low-pass limit | e.g. 10 |
| symmetrize_map | filter controls | num | Centering low-pass limit | e.g. 10 |
| symmetry_test | filter controls | num | Centering low-pass limit | e.g. 10 |
| track_particles | filter controls | num | Centering low-pass limit | e.g. 10 |

## Parameter `center`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | binary | Center class averages | *(empty)* |
| abinitio2D_chunks | search controls | binary | Center class averages | *(empty)* |
| abinitio2D_stream | search controls | binary | Center class averages | *(empty)* |
| abinitio3D | search controls | binary | Center reference volume(s) | *(empty)* |
| abinitio3D_cavgs | search controls | binary | Center reference volume(s) | *(empty)* |
| autorefine3D_nano | search controls | binary | Center reference volume(s) | *(empty)* |
| center2D_nano | search controls | binary | Center class averages | *(empty)* |
| cluster2D_nano | search controls | binary | Center class averages | *(empty)* |
| mask | search controls | binary | Center input volume | *(empty)* |
| refine3D | search controls | binary | Center reference volume(s) | *(empty)* |
| refine3D_auto | search controls | binary | Center reference volume(s) | *(empty)* |
| refine3D_multi | search controls | binary | Center reference volume(s) | *(empty)* |
| refine3D_nano | search controls | binary | Center reference volume(s) | *(empty)* |
| symaxis_search | search controls | binary | Center input volume | *(empty)* |
| symmetrize_map | search controls | binary | Center input volume | *(empty)* |
| symmetry_test | search controls | binary | Center input volume | *(empty)* |

## Parameter `center_pdb`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pdb2mrc | parameter input/output | binary | Whether to move the PDB atomic center to the center of the box | *(empty)* |

## Parameter `chunk_count_penalty`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Chunk-count penalty | e.g. 10 |

## Parameter `chunk_max_len`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Maximum latent chunk length | e.g. 10 |

## Parameter `chunk_max_shift`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Maximum boundary shift | e.g. 10 |

## Parameter `chunk_min_len`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Minimum latent chunk length | e.g. 10 |

## Parameter `chunk_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | multi | Trajectory chunking mode | *(empty)* |

## Parameter `ciffile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cif2mrc | file input/output | file | PDBx/mmCIF input coordinates file | e.g. molecule.cif |
| cif2pdb | file input/output | file | PDBx/mmCIF input coordinates file | e.g. molecule.cif |

## Parameter `class`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | parameter input/output | num | Optional class index to split | e.g. 10 |
| extract_subproj | parameter input/output | num | 2D class index | e.g. 10 |
| stackops | parameter input/output | num | Class index | e.g. 10 |

## Parameter `classtats`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oristats | parameter input/output | binary | Class statistics | *(empty)* |

## Parameter `clip`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| scale | parameter input/output | num | Clipped box size | e.g. 10 |
| stack | parameter input/output | num | Clipped box size | e.g. 10 |

## Parameter `cls_init`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | multi | Scheme for initial class generation | *(empty)* |

## Parameter `clust_crit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cluster_cavgs | search controls | multi | Clustering criterion | *(empty)* |
| cluster_cavgs_selection | search controls | multi | Clustering criterion | *(empty)* |
| cluster_stack | search controls | multi | Clustering criterion | *(empty)* |
| match_cavgs | search controls | multi | Clustering criterion | *(empty)* |
| match_stacks | search controls | multi | Clustering criterion | *(empty)* |

## Parameter `clustind`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract_subproj | parameter input/output | num | Cluster index | e.g. 10 |
| select_clusters | parameter input/output | num | Cluster index | e.g. 10 |

## Parameter `clustinds`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| select_clusters | parameter input/output | str | Comma separated cluster indices | e.g. value |

## Parameter `cn_stop`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| symmetry_test | search controls | num | Rotational symmetry order stop index | e.g. 10 |

## Parameter `combine_eo`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_rec3D | filter controls | binary | Whether e/o references are combined for final alignment(yes\|no){no} | *(empty)* |
| refine3D | filter controls | binary | Whether e/o references are combined for final alignment(yes\|no){no} | *(empty)* |
| refine3D_auto | filter controls | binary | Whether e/o references are combined for final alignment(yes\|no){no} | *(empty)* |

## Parameter `conical_fsc`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | filter controls | binary | Conical FSC regularization | *(empty)* |
| abinitio3D_cavgs | filter controls | binary | Conical FSC regularization | *(empty)* |
| refine3D | filter controls | binary | Conical FSC regularization | *(empty)* |

## Parameter `continue`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D | search controls | binary | Continue previous refinement | *(empty)* |
| refine3D_auto | search controls | binary | Continue previous refinement | *(empty)* |
| refine3D_multi | search controls | binary | Continue previous refinement | *(empty)* |
| refine3D_nano | search controls | binary | Continue previous refinement | *(empty)* |

## Parameter `cs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | num | Spherical aberration | e.g. 10 |
| import_movies | parameter input/output | num | Spherical aberration | e.g. 10 |
| import_particles | parameter input/output | num | Spherical aberration | e.g. 10 |
| master | parameter input/output | float | Spherical aberration (mm) | e.g. 10 |
| mini_stream | parameter input/output | num | Spherical aberration | e.g. 10 |
| preproc | parameter input/output | num | Spherical aberration | e.g. 10 |
| simulate_movie | parameter input/output | num | Spherical aberration | e.g. 10 |
| simulate_particles | parameter input/output | num | Spherical aberration | e.g. 10 |
| tseries_import | parameter input/output | num | Spherical aberration | e.g. 10 |

## Parameter `ctf`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctfops | filter controls | multi | CTF status | *(empty)* |
| import_movies | parameter input/output | multi | CTF status | *(empty)* |
| import_particles | parameter input/output | multi | CTF status | *(empty)* |
| oriops | parameter input/output | multi | CTF status | *(empty)* |
| reimport_particles | parameter input/output | multi | CTF status | *(empty)* |
| simulate_movie | parameter input/output | multi | CTF status | *(empty)* |
| simulate_particles | parameter input/output | multi | CTF status | *(empty)* |

## Parameter `ctf_correct_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctf_correct | parameter input/output | multi | CTF correction mode | *(empty)* |

## Parameter `ctfpatch`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | binary | Patch CTF estimation | *(empty)* |
| ctf_estimate | parameter input/output | binary | Patch CTF estimation | *(empty)* |
| mini_stream | parameter input/output | binary | Patch CTF estimation | *(empty)* |
| preproc | parameter input/output | binary | Patch CTF estimation | *(empty)* |
| preprocess | parameter input/output | binary | Patch CTF estimation | *(empty)* |

## Parameter `ctfresthreshold`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick_extract | filter controls | num | CTF Resolution rejection threshold | e.g. 10 |
| preproc | filter controls | num | CTF Resolution rejection threshold | e.g. 10 |
| selection | parameter input/output | num | CTF Resolution rejection threshold | e.g. 10 |

## Parameter `ctfstats`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oristats | parameter input/output | binary | CTF statistics | *(empty)* |

## Parameter `defocus`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_movie | parameter input/output | num | Underfocus | e.g. 10 |
| simulate_particles | parameter input/output | num | Underfocus | e.g. 10 |

## Parameter `deftab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctfops | file input/output | file | CTF parameter file | e.g. deftab.txt |
| import_movies | file input/output | file | Pre-determined per-micrograph CTF parameters | e.g. deftab.txt |
| import_particles | file input/output | file | CTF parameter file | e.g. deftab.txt |
| import_trajectory | file input/output | file | CTF parameter file | e.g. deftab.txt |
| simulate_particles | file input/output | file | CTF parameter file | e.g. deftab.txt |

## Parameter `deselfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| selection | file input/output | file | File with deselection indices | e.g. deselected.txt |

## Parameter `dferr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | num | Underfocus error half-width | e.g. 10 |
| simulate_particles | parameter input/output | num | Underfocus error half-width | e.g. 10 |

## Parameter `dfmax`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | Expected maximum defocus | e.g. 10 |
| ctf_estimate | search controls | num | Expected maximum defocus | e.g. 10 |
| mini_stream | search controls | num | Expected maximum defocus | e.g. 10 |
| preproc | search controls | num | Expected maximum defocus | e.g. 10 |
| preprocess | search controls | num | Expected maximum defocus | e.g. 10 |

## Parameter `dfmin`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | Expected minimum defocus | e.g. 10 |
| ctf_estimate | search controls | num | Expected minimum defocus | e.g. 10 |
| mini_stream | search controls | num | Expected minimum defocus | e.g. 10 |
| preproc | search controls | num | Expected minimum defocus | e.g. 10 |
| preprocess | search controls | num | Expected minimum defocus | e.g. 10 |
| selection | parameter input/output | num | Expected minimum defocus | e.g. 10 |

## Parameter `dfunit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | parameter input/output | multi | Underfocus unit | *(empty)* |

## Parameter `dir`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | file input/output | dir | Project directory | e.g. project/ |
| pick | file input/output | dir | Output directory | e.g. pick/ |

## Parameter `dir_box`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract | file input/output | dir | Box files directory | e.g. boxes/ |

## Parameter `dir_exec`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D_stream | file input/output | file | Previous run directory | e.g. 3_abinitio2D_stream |
| pick_extract | file input/output | file | Previous run directory | e.g. 2_pick_extract |
| sieve_cavgs | file input/output | file | Previous run directory | e.g. 3_sieve_cavgs |

## Parameter `dir_meta`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | file input/output | dir | Input metadata directory | e.g. /path/to/folder |
| preproc | file input/output | dir | Directory containing per-movie metadata in XML format | e.g. /dataset/metadata |

## Parameter `dir_movies`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_movies | file input/output | dir | Input movies directory | e.g. /cryodata/ |
| master | file input/output | dir | Input movies directory | e.g. /path/to/folder |
| preproc | file input/output | dir | Input movies directory | e.g. /cryodata/ |

## Parameter `dir_preprocess`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | file input/output | hidden_dir | Pre-existing preprocessing directory | *(empty)* |

## Parameter `dir_prev`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| preproc | file input/output | file | Previous run directory | e.g. 2_preproc |

## Parameter `dir_target`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D_stream | file input/output | file | Target directory | e.g. 2_pick_extract |
| assign_optics | file input/output | file | Target directory | e.g. 1_preproc |
| gen_pickrefs | file input/output | file | Target directory | e.g. 1_preproc |
| pick_extract | file input/output | file | Target directory | e.g. 1_preproc |
| sieve_cavgs | file input/output | file | Target directory | e.g. 2_pick_extract |

## Parameter `dm_alpha`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | num | Diffusion-map density normalization (default 0) | e.g. 10 |
| denoise_project | filter controls | num | Diffusion-map density normalization (default 0) | e.g. 10 |
| ppca_denoise | filter controls | num | Diffusion-map density normalization (default 0) | e.g. 10 |

## Parameter `e1`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | num | Rotation along Phi | e.g. 10 |
| volops | parameter input/output | num | Rotation along Phi | e.g. 10 |

## Parameter `e2`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | num | Rotation along Theta | e.g. 10 |
| volops | parameter input/output | num | Rotation along Theta | e.g. 10 |

## Parameter `e3`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | num | Rotation along Psi | e.g. 10 |
| volops | parameter input/output | num | Rotation along Psi | e.g. 10 |

## Parameter `edge`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask | mask controls | num | Envelope mask soft edge | e.g. 10 |
| automask2D | mask controls | num | Envelope mask soft edge | e.g. 10 |
| gen_pickrefs | mask controls | num | Automask soft edge | e.g. 10 |
| mask | mask controls | num | Envelope mask soft edge | e.g. 10 |

## Parameter `eer_fraction`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | parameter input/output | num | # of EER frames to fraction together | e.g. 10 |
| preproc | parameter input/output | num | # of EER frames to fraction together | e.g. 10 |
| preprocess | parameter input/output | num | # of EER frames to fraction together | e.g. 10 |

## Parameter `element`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| analysis2D_nano | parameter input/output | str | Atom element name: Au, Pt etc. | e.g. value |
| atoms_rmsd | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| atoms_stats | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| autorefine3D_nano | parameter input/output | str | Atom element name: Au, Pt etc. | e.g. value |
| conv_atom_denoise | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| core_atoms_analysis | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| crys_score | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| detect_atoms | filter controls | str | Atom element name: Au, Pt etc. | e.g. value |
| simulate_nanoparticle | parameter input/output | str | Atom element name: Au, Pt etc. | e.g. value |

## Parameter `envfsc`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| reconstruct3D | filter controls | binary | Envelope mask e/o maps for FSC | *(empty)* |
| refine3D | filter controls | binary | Envelope mask e/o maps for FSC | *(empty)* |

## Parameter `eo_stage`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | binary | Use even/odd stage | *(empty)* |

## Parameter `even`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_oris | parameter input/output | binary | Generate even projections | *(empty)* |
| simulate_particles | parameter input/output | binary | Generate even projections | *(empty)* |

## Parameter `extr_lim`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | num | Extreme-particle iteration limit | e.g. 10 |

## Parameter `fbody`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fsc_area_score | parameter input/output | string | Output file body | e.g. value |
| motion_correct | parameter input/output | string | Template output micrograph name | e.g. value |
| preprocess | parameter input/output | string | Template output micrograph name | e.g. value |
| track_particles | parameter input/output | string | Template output tracked series | e.g. value |
| tseries_extractor | parameter input/output | string | Output file body | e.g. value |

## Parameter `filetab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | file input/output | file | List of files | e.g. input_files.txt |
| import_movies | file input/output | file | List of movie files | e.g. input_files.txt |
| mini_stream | file input/output | file | List of files | e.g. input_files.txt |
| model_cavgs_rejection | file input/output | file | Analysis file table | e.g. input_files.txt |
| scale | file input/output | file | Stacks list | e.g. input_files.txt |
| stack | file input/output | file | Stacks list | e.g. input_files.txt |
| tsegmaps_core_finder | file input/output | file | Volumes list | e.g. input_files.txt |
| tseries_import | file input/output | file | List of individual movie frame files | e.g. input_files.txt |
| volanalyze | file input/output | file | Volumes list | e.g. input_files.txt |
| volcluster | file input/output | file | Volumes list | e.g. input_files.txt |

## Parameter `fill_holes`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| binarize | parameter input/output | binary | Fill holes | *(empty)* |

## Parameter `filt_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | filter controls | multi | Filtering mode | *(empty)* |
| refine3D_auto | filter controls | multi | Filtering mode | *(empty)* |
| refine3D_multi | filter controls | multi | Filtering mode | *(empty)* |

## Parameter `filter`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | multi | Filter type(bs\|nlmean\|no){no} | *(empty)* |
| track_particles | filter controls | multi | Alternative filter for particle tracking | *(empty)* |

## Parameter `fit_phshift`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | binary | Fit CTF phase shift | *(empty)* |
| ctf_estimate | parameter input/output | binary | Fit CTF phase shift | *(empty)* |
| master | parameter input/output | binary | Phase plate | *(empty)* |
| mini_stream | parameter input/output | binary | Fit CTF phase shift | *(empty)* |
| preproc | parameter input/output | binary | Fit CTF phase shift | *(empty)* |
| preprocess | parameter input/output | binary | Fit CTF phase shift | *(empty)* |

## Parameter `flipgain`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fractionate_movies | parameter input/output | multi | Flip the gain reference | *(empty)* |
| master | image input/output | multi | Gain processing | *(empty)* |
| motion_correct | parameter input/output | multi | Flip the gain reference | *(empty)* |
| preproc | parameter input/output | multi | Flip the gain reference | *(empty)* |
| preprocess | parameter input/output | multi | Flip the gain reference | *(empty)* |

## Parameter `fname`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_register | file input/output | file | PDB file list | e.g. pdb_files.txt |
| crys_score | file input/output | file | PDB file list | e.g. np_pdbs.txt |
| info_image | file input/output | file | Name of image file | e.g. image.mrc |
| model_cavgs_rejection | file input/output | file | Quality model output | e.g. cavgs_quality_output.txt |
| write_mic_filetab | file input/output | file | Filename micrograph list | e.g. mics.txt |

## Parameter `force_lp_range`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | filter controls | binary | Force low-pass range | *(empty)* |

## Parameter `frac_best`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_cavgs | parameter input/output | num | Anchor fraction | e.g. 10 |
| sample_classes | parameter input/output | num | Fraction of best particles to sample from | e.g. 10 |

## Parameter `frac_diam`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_rmsd | parameter input/output | num | Fraction of atomic diameter | e.g. 10 |
| core_atoms_analysis | parameter input/output | num | Fraction of atomic diameter | e.g. 10 |

## Parameter `frac_worst`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| sample_classes | parameter input/output | num | Fraction of worst particles to sample from | e.g. 10 |

## Parameter `fraca`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| import_movies | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| import_particles | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| master | parameter input/output | float | Amplitude contrast fraction | e.g. 10 |
| mini_stream | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| preproc | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| simulate_movie | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| simulate_particles | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |
| tseries_import | parameter input/output | num | Amplitude contrast fraction | e.g. 10 |

## Parameter `fraction_dose_target`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | parameter input/output | num | EER fraction dose target (e/Ang^2) | e.g. 10 |
| preproc | parameter input/output | num | EER fraction dose target (e/Ang^2) | e.g. 10 |
| preprocess | parameter input/output | num | EER fraction dose target (e/Ang^2) | e.g. 10 |

## Parameter `frcs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | str | Projection FRCs file | e.g. value |
| print_fsc | parameter input/output | str | Projection FRCs file | e.g. value |

## Parameter `fromf`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fractionate_movies | parameter input/output | num | Starting frame | e.g. 10 |
| track_particles | parameter input/output | num | Frame to start tracking from | e.g. 10 |
| tseries_make_pickavg | parameter input/output | num | Frame to start averaging from | e.g. 10 |

## Parameter `fromp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract_subproj | parameter input/output | num | From index | e.g. 10 |
| extract_substk | parameter input/output | num | From index | e.g. 10 |
| scale | parameter input/output | num | First stack index | e.g. 10 |
| stackops | parameter input/output | num | From particle index | e.g. 10 |
| trajectory_reconstruct3D | parameter input/output | num | From particle index | e.g. 10 |

## Parameter `fsc`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctf_correct | filter controls | file | FSC file | e.g. fsc.bin |
| filter | filter controls | file | FSC file | e.g. fsc.bin |
| postprocess | filter controls | file | FSC file | e.g. fsc.bin |
| print_fsc | file input/output | file | FSC file | e.g. fsc.bin |

## Parameter `gainref`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | file input/output | file | Gain reference | e.g. gainref.mrc |
| motion_correct | file input/output | file | Gain reference | e.g. gainref.mrc |
| preproc | file input/output | file | Gain reference | e.g. gainref.mrc |
| preprocess | file input/output | file | Gain reference | e.g. gainref.mrc |

## Parameter `graph`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | multi | Class split graph | *(empty)* |
| denoise_project | filter controls | multi | Diffusion graph | *(empty)* |

## Parameter `greedy_sampling`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| sample_classes | parameter input/output | binary | Greedy balanced selection | *(empty)* |
| selection | parameter input/output | binary | Greedy balanced selection | *(empty)* |

## Parameter `guinier`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| volops | parameter input/output | binary | Guinier plot | *(empty)* |

## Parameter `hp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | num | High-pass limit | e.g. 10 |
| abinitio2D_chunks | filter controls | num | High-pass limit | e.g. 10 |
| abinitio3D | filter controls | num | High-pass limit | e.g. 10 |
| abinitio3D_cavgs | filter controls | num | High-pass limit | e.g. 10 |
| abinitio3D_nano | filter controls | num | High-pass limit | e.g. 10 |
| autorefine3D_nano | filter controls | num | High-pass limit | e.g. 10 |
| check_refpick | filter controls | num | High-pass limit | e.g. 10 |
| cluster2D_nano | filter controls | num | High-pass limit | e.g. 10 |
| cluster_cavgs | filter controls | num | High-pass limit | e.g. 10 |
| cluster_cavgs_selection | filter controls | num | High-pass limit | e.g. 10 |
| cluster_stack | filter controls | num | High-pass limit | e.g. 10 |
| ctf_estimate | filter controls | num | High-pass limit | e.g. 10 |
| dock_volpair | filter controls | num | High-pass limit | e.g. 10 |
| filter | filter controls | num | High-pass limit | e.g. 10 |
| fsc | filter controls | num | High-pass limit | e.g. 10 |
| match_cavgs | filter controls | num | High-pass limit | e.g. 10 |
| match_stacks | filter controls | num | High-pass limit | e.g. 10 |
| mini_stream | filter controls | num | High-pass limit | e.g. 10 |
| ppca_denoise | filter controls | num | High-pass limit | e.g. 10 |
| refine3D | filter controls | num | High-pass limit | e.g. 10 |
| refine3D_nano | filter controls | num | High-pass limit | e.g. 10 |
| symaxis_search | filter controls | num | High-pass limit | e.g. 10 |
| symmetrize_map | filter controls | num | High-pass limit | e.g. 10 |
| symmetry_test | filter controls | num | High-pass limit | e.g. 10 |
| track_particles | filter controls | num | High-pass limit | e.g. 10 |
| trajectory_denoise | filter controls | num | High-pass limit | e.g. 10 |
| volanalyze | filter controls | num | High-pass limit | e.g. 10 |
| volcluster | filter controls | num | High-pass limit | e.g. 10 |
| volops | filter controls | num | High-pass limit | e.g. 10 |

## Parameter `hp_ctf_estimate`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| preproc | filter controls | num | CTF estimation high-pass limit | e.g. 10 |
| preprocess | filter controls | num | High-pass limit for CTF parameter estimation | e.g. 10 |

## Parameter `icefracthreshold`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick_extract | filter controls | num | Ice Fraction rejection threshold | e.g. 10 |
| selection | parameter input/output | num | Ice Fraction rejection threshold | e.g. 10 |

## Parameter `icm`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | filter controls | binary | Whether to perform ICM filtering of reference(s) | *(empty)* |
| flex_analysis | filter controls | binary | Automatic diffusion-mode selection | *(empty)* |
| refine3D_nano | filter controls | binary | Whether to perform ICM filtering of reference(s) | *(empty)* |
| uniform_filter3D | filter controls | binary | Whether to perform ICM filtering of reference(s) | *(empty)* |

## Parameter `imgkind`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| postprocess | parameter input/output | str | Volume kind | e.g. value |

## Parameter `import_dir`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_starproject | file input/output | dir | Import directory | e.g. MotionCorr/job001 |

## Parameter `import_type`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_starproject | parameter input/output | multi | Import type | *(empty)* |

## Parameter `infile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| model_cavgs_rejection | file input/output | file | Quality model input | e.g. cavgs_quality_model.txt |
| selection | file input/output | file | File with selection state (0/1) flags | e.g. selection.txt |
| tseries_extractor | file input/output | file | Selected particle trajectory | e.g. trajectory.txt |

## Parameter `job_memory_per_task`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | Memory per computing node | e.g. value |
| update_project | computer controls | str | Memory per computing node | e.g. value |

## Parameter `json`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| print_project_field | parameter input/output | binary | output in JSON format | *(empty)* |

## Parameter `k_nn`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | parameter input/output | num | Diffusion graph neighbors (default 10; try 5-30) | e.g. 10 |
| denoise_project | filter controls | num | Diffusion graph neighbors (default 10; try 5-30) | e.g. 10 |
| flex_analysis | filter controls | num | Nearest neighbors (default 100) | e.g. 10 |
| ppca_denoise | filter controls | num | Diffusion graph neighbors (default 5; try 5-30) | e.g. 10 |
| trajectory_denoise | filter controls | num | Diffusion graph neighbors (default 5; try 5-30) | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Nearest neighbors | e.g. 10 |

## Parameter `kpca_backend`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | multi | Kernel PCA backend | *(empty)* |
| ppca_denoise_classes | filter controls | multi | Kernel PCA backend | *(empty)* |
| trajectory_denoise | filter controls | multi | Kernel PCA backend | *(empty)* |

## Parameter `kpca_ker`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | multi | Kernel PCA kernel | *(empty)* |
| ppca_denoise_classes | filter controls | multi | Kernel PCA kernel | *(empty)* |
| trajectory_denoise | filter controls | multi | Kernel PCA kernel | *(empty)* |

## Parameter `kpca_nystrom_local_nbrs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | num | Nyström max local support neighbors (default 96; try 96, 128) | e.g. 10 |
| ppca_denoise_classes | filter controls | num | Nyström max local support neighbors (default 96; try 96, 128) | e.g. 10 |
| trajectory_denoise | filter controls | num | Nyström max local support neighbors (default 96; try 96, 128) | e.g. 10 |

## Parameter `kpca_nystrom_npts`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | num | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | e.g. 10 |
| ppca_denoise_classes | filter controls | num | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | e.g. 10 |
| trajectory_denoise | filter controls | num | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | e.g. 10 |

## Parameter `kpca_rbf_gamma`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | num | RBF gamma (0 => auto) | e.g. 10 |
| ppca_denoise_classes | filter controls | num | RBF gamma (0 => auto) | e.g. 10 |
| trajectory_denoise | filter controls | num | RBF gamma (0 => auto) | e.g. 10 |

## Parameter `kv`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | num | Acceleration voltage | e.g. 10 |
| import_movies | parameter input/output | num | Acceleration voltage | e.g. 10 |
| import_particles | parameter input/output | num | Acceleration voltage | e.g. 10 |
| master | parameter input/output | int | Acceleration voltage (kV) | e.g. 10 |
| mini_stream | parameter input/output | num | Acceleration voltage | e.g. 10 |
| preproc | parameter input/output | num | Acceleration voltage | e.g. 10 |
| print_dose_weights | parameter input/output | num | Acceleration voltage | e.g. 10 |
| simulate_movie | parameter input/output | num | Acceleration voltage | e.g. 10 |
| simulate_particles | parameter input/output | num | Acceleration voltage | e.g. 10 |
| tseries_import | parameter input/output | num | Acceleration voltage | e.g. 10 |

## Parameter `lambda`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | num | BS smoother lambda | e.g. 10 |
| icm2D | filter controls | num | ICM lambda regularization parameter | e.g. 10 |
| icm3D | filter controls | num | ICM lambda regularization parameter | e.g. 10 |
| refine3D_nano | filter controls | num | ICM lambda regularization parameter | e.g. 10 |
| uniform_filter3D | filter controls | num | ICM lambda regularization parameter | e.g. 10 |

## Parameter `lp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | num | Low-pass limit | e.g. 10 |
| abinitio2D_chunks | filter controls | num | Low-pass limit | e.g. 10 |
| abinitio3D | filter controls | num | Low-pass limit | e.g. 10 |
| autorefine3D_nano | filter controls | num | Initial low-pass limit | e.g. 10 |
| check_refpick | filter controls | num | Low-pass limit | e.g. 10 |
| cluster2D_nano | filter controls | num | Static low-pass limit | e.g. 10 |
| cluster_cavgs | filter controls | num | Low-pass limit | e.g. 10 |
| cluster_cavgs_selection | filter controls | num | Low-pass limit | e.g. 10 |
| cluster_stack | filter controls | num | Low-pass limit | e.g. 10 |
| ctf_estimate | filter controls | num | Low-pass limit | e.g. 10 |
| dock_volpair | filter controls | num | Low-pass limit | e.g. 10 |
| estimate_diam | filter controls | num | low-pass limit in Angstroms{7.} | e.g. 10 |
| filter | filter controls | num | Low-pass limit | e.g. 10 |
| flex_analysis | filter controls | num | Low-pass limit | e.g. 10 |
| fsc | filter controls | num | Low-pass limit | e.g. 10 |
| match_cavgs | filter controls | num | Low-pass limit | e.g. 10 |
| match_stacks | filter controls | num | Low-pass limit | e.g. 10 |
| mini_stream | filter controls | num | Low-pass limit | e.g. 10 |
| pick | filter controls | num | Low-pass limit | e.g. 10 |
| postprocess | filter controls | num | Low-pass limit for map filtering | e.g. 10 |
| ppca_denoise | filter controls | num | Low-pass limit | e.g. 10 |
| refine3D | filter controls | num | Static low-pass limit | e.g. 10 |
| refine3D_nano | filter controls | num | Static low-pass limit | e.g. 10 |
| symaxis_search | filter controls | num | Low-pass limit | e.g. 10 |
| symmetrize_map | filter controls | num | Low-pass limit | e.g. 10 |
| symmetry_test | filter controls | num | Low-pass limit | e.g. 10 |
| track_particles | filter controls | num | Low-pass limit in Angs | e.g. 10 |
| trajectory_denoise | filter controls | num | Low-pass limit | e.g. 10 |
| trajectory_reconstruct3D | filter controls | num | Low-pass limit | e.g. 10 |
| volanalyze | filter controls | num | Low-pass limit | e.g. 10 |
| volcluster | filter controls | num | Low-pass limit | e.g. 10 |
| volops | filter controls | num | Low-pass limit | e.g. 10 |

## Parameter `lp_backgr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| mask | filter controls | num | Background low-pass resolution | e.g. 10 |
| refine3D | filter controls | num | Background low-pass resolution | e.g. 10 |

## Parameter `lp_ctf_estimate`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| preproc | filter controls | num | CTF estimation low-pass limit | e.g. 10 |
| preprocess | filter controls | num | Low-pass limit for CTF parameter estimation | e.g. 10 |

## Parameter `lp_pick`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick_extract | filter controls | num | Low-pass limit for picking | e.g. 10 |

## Parameter `lplim_crit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fsc_area_score | filter controls | num | Low-pass limit FSC criterion | e.g. 10 |
| refine3D | filter controls | num | Low-pass limit FSC criterion | e.g. 10 |

## Parameter `lpstart`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | num | Initial low-pass limit | e.g. 10 |
| abinitio3D | filter controls | num | Starting low-pass limit | e.g. 10 |
| abinitio3D_cavgs | filter controls | num | Starting low-pass limit | e.g. 10 |
| abinitio3D_nano | filter controls | num | Starting low-pass limit | e.g. 10 |
| cluster2D_nano | filter controls | num | Initial low-pass limit | e.g. 10 |
| motion_correct | filter controls | num | Initial low-pass limit | e.g. 10 |
| particle_sieving | filter controls | num | Initial low-pass limit | e.g. 10 |
| preproc | filter controls | num | Motion-correction low-pass start | e.g. 10 |
| preprocess | filter controls | num | Initial low-pass limit for movie alignment | e.g. 10 |
| sieve_cavgs | filter controls | num | Initial sieving low-pass limit | e.g. 10 |
| tseries_make_pickavg | filter controls | num | Initial low-pass limit | e.g. 10 |
| tseries_motion_correct | filter controls | num | Initial low-pass limit | e.g. 10 |
| uniform_filter2D | filter controls | num | Starting resolution limit | e.g. 10 |
| uniform_filter3D | filter controls | num | Starting resolution limit | e.g. 10 |

## Parameter `lpstart_ini3D`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | filter controls | num | Starting low-pass limit ini3D | e.g. 10 |

## Parameter `lpstop`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | num | Final low-pass limit | e.g. 10 |
| abinitio2D_chunks | filter controls | num | Final low-pass limit | e.g. 10 |
| abinitio3D | filter controls | num | Final low-pass limit | e.g. 10 |
| abinitio3D_cavgs | filter controls | num | Final low-pass limit | e.g. 10 |
| abinitio3D_nano | filter controls | num | Final low-pass limit | e.g. 10 |
| cluster2D_nano | filter controls | num | Final low-pass limit | e.g. 10 |
| motion_correct | filter controls | num | Final low-pass limit | e.g. 10 |
| preproc | filter controls | num | Motion-correction low-pass stop | e.g. 10 |
| preprocess | filter controls | num | Final low-pass limit for movie alignment | e.g. 10 |
| refine3D | filter controls | num | Low-pass limit for frequency limited refinement | e.g. 10 |
| refine3D_multi | filter controls | num | Low-pass limit for frequency limited refinement | e.g. 10 |
| refine3D_nano | filter controls | num | Low-pass limit for frequency limited refinement | e.g. 10 |
| tseries_make_pickavg | filter controls | num | Final low-pass limit | e.g. 10 |
| tseries_motion_correct | filter controls | num | Final low-pass limit | e.g. 10 |
| uniform_filter2D | filter controls | num | Stopping resolution limit | e.g. 10 |
| uniform_filter3D | filter controls | num | Stopping resolution limit | e.g. 10 |

## Parameter `lpstop_coarse`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | filter controls | num | Final low-pass limit for coarse sieving | e.g. 10 |
| sieve_cavgs | filter controls | num | Coarse-pass low-pass limit | e.g. 10 |

## Parameter `lpstop_fine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | filter controls | num | Final low-pass limit for fine sieving | e.g. 10 |
| sieve_cavgs | filter controls | num | Fine-pass low-pass limit | e.g. 10 |

## Parameter `lpstop_ini3D`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | filter controls | num | Final low-pass limit ini3D | e.g. 10 |

## Parameter `makemovie`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| stackops | parameter input/output | binary | Whether to make a movie | *(empty)* |

## Parameter `max_dose`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | parameter input/output | num | Maximum dose threshold(e/A2) | e.g. 10 |
| preproc | parameter input/output | num | Maximum dose threshold(e/A2) | e.g. 10 |
| preprocess | parameter input/output | num | Maximum dose threshold(e/A2) | e.g. 10 |

## Parameter `maxits`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | search controls | num | Max iterations | e.g. 10 |
| cluster2D_nano | search controls | num | Max iterations | e.g. 10 |
| reconstruct3D_pcg | filter controls | num | Max iterations | e.g. 10 |
| refine3D | search controls | num | Max iterations | e.g. 10 |
| refine3D_auto | search controls | num | Max iterations | e.g. 10 |
| refine3D_multi | search controls | num | Max iterations | e.g. 10 |
| refine3D_nano | search controls | num | Max iterations | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Max iterations | e.g. 10 |

## Parameter `maxnchunks`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Max # of chunks to process | e.g. 10 |
| sieve_cavgs | search controls | num | Number of subsets after which 2D analysis ends | e.g. 10 |

## Parameter `maxpop`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| assign_optics_groups | parameter input/output | num | Maximum number of movies/micrographs/stacks in each optics group | e.g. 10 |

## Parameter `mcconvention`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fractionate_movies | parameter input/output | str | Movie alignment convention | e.g. value |
| motion_correct | search controls | str | Frame of reference during movie alignment | e.g. value |
| preproc | parameter input/output | str | Frame of reference during movie alignment | e.g. value |
| preprocess | search controls | str | Frame of reference during movie alignment | e.g. value |

## Parameter `mcpatch`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | search controls | binary | Patch-based motion correction | *(empty)* |
| preproc | parameter input/output | binary | Patch-based motion correction | *(empty)* |
| preprocess | search controls | binary | Patch-based motion correction | *(empty)* |
| tseries_make_pickavg | search controls | binary | Patch-based motion correction | *(empty)* |
| tseries_motion_correct | search controls | binary | Patch-based motion correction | *(empty)* |

## Parameter `mcpatch_thres`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | search controls | binary | Use motion correction patch threshold | *(empty)* |
| preproc | parameter input/output | binary | Use motion correction patch threshold | *(empty)* |
| preprocess | search controls | binary | Use motion correction patch threshold | *(empty)* |

## Parameter `minits`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D_auto | search controls | num | Minimum automatic iterations | e.g. 10 |

## Parameter `mirr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | multi | Mirror orientations | *(empty)* |
| postprocess | filter controls | multi | Perform mirroring | *(empty)* |
| stackops | parameter input/output | multi | Perform mirroring | *(empty)* |
| volops | parameter input/output | multi | Perform mirroring | *(empty)* |

## Parameter `ml_reg`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | filter controls | binary | ML regularization | *(empty)* |
| refine3D | filter controls | binary | ML regularization | *(empty)* |
| refine3D_multi | filter controls | binary | ML regularization | *(empty)* |

## Parameter `moldiam`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cluster2D_nano | parameter input/output | num | Molecular diameter | e.g. 10 |
| pick | parameter input/output | num | Molecular diameter | e.g. 10 |
| pick_extract | parameter input/output | num | Molecular diameter | e.g. 10 |
| print_magic_boxes | parameter input/output | num | Molecular diameter | e.g. 10 |
| simulate_nanoparticle | parameter input/output | num | Molecular diameter | e.g. 10 |

## Parameter `moldiam_max`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| mini_stream | parameter input/output | num | Max molecular diameter | e.g. 10 |
| pick | parameter input/output | num | Max molecular diameter | e.g. 10 |
| pick_extract | parameter input/output | num | Max molecular diameter | e.g. 10 |

## Parameter `mskdiam`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | mask controls | num | Mask diameter | e.g. 10 |
| abinitio2D_chunks | mask controls | num | Mask diameter | e.g. 10 |
| abinitio2D_stream | mask controls | num | Mask diameter | e.g. 10 |
| abinitio3D | mask controls | num | Mask diameter | e.g. 10 |
| abinitio3D_cavgs | mask controls | num | Mask diameter | e.g. 10 |
| abinitio3D_nano | mask controls | num | Mask diameter | e.g. 10 |
| automask2D | mask controls | num | Mask diameter | e.g. 10 |
| autorefine3D_nano | mask controls | num | Mask diameter | e.g. 10 |
| bootstrap_rec3D | mask controls | num | Mask diameter | e.g. 10 |
| cavgseoproc_nano | mask controls | num | Mask diameter | e.g. 10 |
| cavgsproc_nano | mask controls | num | Mask diameter | e.g. 10 |
| cls_split | mask controls | num | Mask diameter | e.g. 10 |
| cluster2D_nano | mask controls | num | Mask diameter | e.g. 10 |
| cluster_cavgs | mask controls | num | Mask diameter | e.g. 10 |
| cluster_cavgs_selection | mask controls | num | Mask diameter | e.g. 10 |
| cluster_stack | mask controls | num | Mask diameter | e.g. 10 |
| conv_atom_denoise | mask controls | num | Mask diameter | e.g. 10 |
| denoise_project | mask controls | num | Mask diameter | e.g. 10 |
| detect_atoms | mask controls | num | Mask diameter | e.g. 10 |
| dock_volpair | mask controls | num | Mask diameter | e.g. 10 |
| estimate_diam | mask controls | num | Mask diameter | e.g. 10 |
| flex_analysis | mask controls | num | Mask diameter | e.g. 10 |
| fsc | mask controls | num | Mask diameter | e.g. 10 |
| fsc_area_score | mask controls | num | Mask diameter | e.g. 10 |
| mask | mask controls | num | Mask diameter | e.g. 10 |
| match_cavgs | mask controls | num | Mask diameter | e.g. 10 |
| match_stacks | mask controls | num | Mask diameter | e.g. 10 |
| model_cavgs_rejection | mask controls | num | Mask diameter | e.g. 10 |
| normalize | mask controls | num | Mask diameter | e.g. 10 |
| nu_filt3D | mask controls | num | Mask diameter | e.g. 10 |
| postprocess | mask controls | num | Mask diameter | e.g. 10 |
| ppca_denoise | mask controls | num | Mask diameter | e.g. 10 |
| ppca_volvar | mask controls | num | Mask diameter | e.g. 10 |
| ptclsproc_nano | mask controls | num | Mask diameter | e.g. 10 |
| reconstruct3D | mask controls | num | Mask diameter | e.g. 10 |
| reconstruct3D_pcg | mask controls | num | Mask diameter | e.g. 10 |
| refine3D | mask controls | num | Mask diameter | e.g. 10 |
| refine3D_auto | mask controls | num | Mask diameter | e.g. 10 |
| refine3D_multi | mask controls | num | Mask diameter | e.g. 10 |
| refine3D_nano | mask controls | num | Mask diameter | e.g. 10 |
| reproject | mask controls | num | Mask diameter | e.g. 10 |
| sieve_cavgs | mask controls | num | Mask diameter | e.g. 10 |
| simulate_particles | mask controls | num | Mask diameter | e.g. 10 |
| symaxis_search | mask controls | num | Mask diameter | e.g. 10 |
| symmetrize_map | mask controls | num | Mask diameter | e.g. 10 |
| symmetry_test | mask controls | num | Mask diameter | e.g. 10 |
| trajectory_denoise | mask controls | num | Mask diameter | e.g. 10 |
| trajectory_make_projavgs | mask controls | num | Mask diameter | e.g. 10 |
| trajectory_reconstruct3D | mask controls | num | Mask diameter | e.g. 10 |
| uniform_filter2D | mask controls | num | Mask diameter | e.g. 10 |
| uniform_filter3D | mask controls | num | Mask diameter | e.g. 10 |
| validate_cavgs_vs_model | mask controls | num | Mask diameter | e.g. 10 |
| volanalyze | mask controls | num | Mask diameter | e.g. 10 |
| volcluster | mask controls | num | Mask diameter | e.g. 10 |
| volops | mask controls | num | Mask diameter | e.g. 10 |

## Parameter `mskdiam_detect`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | mask controls | num | Detect-atoms mask diameter | e.g. 10 |

## Parameter `mul`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_cavgs | parameter input/output | num | Shift multiplication factor | e.g. 10 |
| oriops | parameter input/output | num | Shift multiplication factor | e.g. 10 |
| volops | parameter input/output | num | Multiplication factor | e.g. 10 |

## Parameter `multi_moldiams`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick | parameter input/output | str | Comma-separated molecular diameters with which to execute multiple gaussian pick | e.g. value |

## Parameter `multivol_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | multi | Multi-volume ab initio mode | *(empty)* |
| abinitio3D_cavgs | search controls | multi | Multi-volume class-average ab initio mode | *(empty)* |
| refine3D_multi | search controls | multi | Multi-volume refinement mode | *(empty)* |

## Parameter `nang_nbrs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| flex_analysis | filter controls | num | Angular candidate cap (default 1000) | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Angular candidate cap | e.g. 10 |

## Parameter `nboxes_max`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | Max # boxes per micrograph | e.g. 10 |
| pick | search controls | num | Max # boxes per micrograph | e.g. 10 |

## Parameter `nchunks`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D_chunks | parameter input/output | num | Number of chunks | e.g. 10 |
| particle_sieving | computer controls | num | Legacy concurrent-chunk alias | e.g. 10 |
| sieve_cavgs | computer controls | num | Number of subsets to classify simultaneously | e.g. 10 |
| trajectory_reconstruct3D | parameter input/output | num | Number of temporal chunks | e.g. 10 |

## Parameter `nchunks_max`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Maximum automatic chunk count | e.g. 10 |

## Parameter `nchunks_min`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Minimum automatic chunk count | e.g. 10 |

## Parameter `nchunksperset`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| sieve_cavgs | search controls | num | Number of subsets to group | e.g. 10 |

## Parameter `ncls`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | num | Number of 2D clusters | e.g. 10 |
| abinitio2D_stream | search controls | num | Maximum number of 2D clusters | e.g. 10 |
| center2D_nano | search controls | num | Number of 2D clusters | e.g. 10 |
| cls_split | parameter input/output | num | Fixed number of subclasses (0 => auto) | e.g. 10 |
| cluster_cavgs | parameter input/output | num | Number of clusters | e.g. 10 |
| cluster_stack | parameter input/output | num | Number of clusters | e.g. 10 |
| make_cavgs | parameter input/output | num | Number of 2D clusters | e.g. 10 |
| make_oris | parameter input/output | num | Number of random class labels | e.g. 10 |
| pick | search controls | num | Cluster input pickrefs before template generation | e.g. 10 |
| sieve_cavgs | search controls | num | Number of 2D clusters | e.g. 10 |
| volcluster | parameter input/output | num | Number of 2D clusters | e.g. 10 |

## Parameter `ncls_coarse`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Number of clusters for coarse sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Coarse-pass class count | e.g. 10 |

## Parameter `ncls_fine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Number of clusters for fine sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Fine-pass class count | e.g. 10 |

## Parameter `ncunits`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| track_particles | computer controls | num | Concurrent particle trackers | e.g. 10 |

## Parameter `ndev`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| binarize | parameter input/output | num | Binarization threshold | e.g. 10 |
| pick | search controls | num | # of sigmas for outlier detection | e.g. 10 |

## Parameter `ndiscrete`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_oris | parameter input/output | num | Number of discrete projection directions | e.g. 10 |
| oriops | parameter input/output | num | Number of discrete projection directions | e.g. 10 |

## Parameter `neg`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ctfops | parameter input/output | binary | Invert contrast | *(empty)* |
| reproject | parameter input/output | binary | Invert contrast | *(empty)* |
| stackops | parameter input/output | binary | Invert contrast | *(empty)* |
| track_particles | parameter input/output | binary | Invert contrast | *(empty)* |
| tseries_extractor | parameter input/output | binary | Invert contrast | *(empty)* |
| volops | parameter input/output | binary | Invert contrast | *(empty)* |

## Parameter `neigs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | num | Number of eigencomponents (0 => auto scan; default 200) | e.g. 10 |
| denoise_project | filter controls | num | Number of eigencomponents (0 => auto scan; default 200) | e.g. 10 |
| flex_analysis | filter controls | num | Maximum number of diffusion modes (default 15) | e.g. 10 |
| ppca_denoise | filter controls | num | Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160) | e.g. 10 |
| ppca_denoise_classes | filter controls | num | # eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160) | e.g. 10 |
| ppca_volvar | filter controls | num | # eigenvecs | e.g. 10 |
| trajectory_denoise | filter controls | num | Number of diffusion-map components (0 => auto; default 0) | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Flex latent dimensions | e.g. 10 |

## Parameter `newbox`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| scale | parameter input/output | num | Scaled box size | e.g. 10 |

## Parameter `nframes`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| print_dose_weights | parameter input/output | num | Number of frames | e.g. 10 |
| simulate_movie | parameter input/output | num | Number of frames | e.g. 10 |
| simulate_particles | parameter input/output | num | # of particle frames | e.g. 10 |

## Parameter `nframesgrp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| stackops | parameter input/output | num | Number of stack entries to group & average | e.g. 10 |
| track_particles | search controls | num | Number of contigous frames to average | e.g. 10 |
| tseries_make_pickavg | parameter input/output | num | # contigous frames to average | e.g. 10 |
| tseries_motion_correct | search controls | num | # frames in time moving time window | e.g. 10 |

## Parameter `ngrow`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask2D | mask controls | num | # layers to grow | e.g. 10 |
| gen_pickrefs | mask controls | num | Automask growth layers | e.g. 10 |

## Parameter `nicedispid`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | parameter input/output | hidden_int | Optics group offset delta multiplier | *(empty)* |

## Parameter `ninipick`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| preproc | parameter input/output | num | Number of micrographs to perform initial picking preprocessing on | e.g. 10 |

## Parameter `nits_per_stage`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | num | Iterations per stage | e.g. 10 |

## Parameter `nmics`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| gen_pickrefs | parameter input/output | num | Number of micrographs to import | e.g. 10 |
| particle_sieving | search controls | num | Max # of micrographs per chunk | e.g. 10 |
| sieve_cavgs | parameter input/output | num | Micrographs per sieve cycle | e.g. 10 |

## Parameter `nmoldiams`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick | parameter input/output | num | Number of molecular diameters to investigate | e.g. 10 |
| pick_extract | parameter input/output | num | Number of molecular diameters to investigate | e.g. 10 |

## Parameter `noise_norm`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| normalize | parameter input/output | binary | Noise normalize | *(empty)* |

## Parameter `norm`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| normalize | parameter input/output | binary | Normalize | *(empty)* |

## Parameter `nparts`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | computer controls | num | Number of computing nodes | e.g. 10 |
| abinitio2D_chunks | computer controls | num | Number of computing nodes | e.g. 10 |
| abinitio2D_stream | computer controls | num | Number of computing nodes | e.g. 10 |
| abinitio3D | computer controls | num | Number of computing nodes | e.g. 10 |
| abinitio3D_cavgs | computer controls | num | Number of computing nodes | e.g. 10 |
| abinitio3D_nano | computer controls | num | Number of computing nodes | e.g. 10 |
| bootstrap_cavgs | computer controls | num | Number of computing nodes | e.g. 10 |
| bootstrap_rec3D | computer controls | num | Number of computing nodes | e.g. 10 |
| cls_split | computer controls | num | Number of computing nodes | e.g. 10 |
| cluster2D_nano | computer controls | num | Number of computing nodes | e.g. 10 |
| ctf_estimate | computer controls | num | Number of computing nodes | e.g. 10 |
| denoise_project | computer controls | num | Number of computing nodes | e.g. 10 |
| extract | computer controls | num | Number of computing nodes | e.g. 10 |
| flex_analysis | computer controls | num | Number of computing nodes | e.g. 10 |
| fractionate_movies | computer controls | num | Number of computing nodes | e.g. 10 |
| gen_pspecs_and_thumbs | computer controls | num | Number of computing nodes | e.g. 10 |
| make_cavgs | computer controls | num | Number of computing nodes | e.g. 10 |
| motion_correct | computer controls | num | Number of computing nodes | e.g. 10 |
| particle_sieving | computer controls | num | Number of chunks classified simultaneously | e.g. 10 |
| pick | computer controls | num | Number of computing nodes | e.g. 10 |
| pick_extract | computer controls | num | Number of computing nodes | e.g. 10 |
| preproc | computer controls | num | Number of computing nodes | e.g. 10 |
| preprocess | computer controls | num | Number of computing nodes | e.g. 10 |
| prune_project | computer controls | num | Number of computing nodes | e.g. 10 |
| reconstruct3D | computer controls | num | Number of computing nodes | e.g. 10 |
| reextract | computer controls | num | Number of computing nodes | e.g. 10 |
| refine3D | computer controls | num | Number of computing nodes | e.g. 10 |
| refine3D_auto | computer controls | num | Number of computing nodes | e.g. 10 |
| refine3D_multi | computer controls | num | Number of computing nodes | e.g. 10 |
| refine3D_nano | computer controls | num | Number of computing nodes | e.g. 10 |
| sample_classes | parameter input/output | num | Number of partitions in balancing | e.g. 10 |
| selection | parameter input/output | num | Number of partitions in balancing | e.g. 10 |
| sieve_cavgs | computer controls | num | Number of computing nodes | e.g. 10 |
| split | computer controls | num | Number of computing nodes | e.g. 10 |
| split_stack | parameter input/output | num | Number of parts balanced splitting of the stack | e.g. 10 |
| trajectory_reconstruct3D | computer controls | num | Number of computing nodes | e.g. 10 |
| tseries_motion_correct | computer controls | num | Number of computing nodes | e.g. 10 |

## Parameter `npreimages`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| flex_analysis | filter controls | num | Representative state volumes (default 8) | e.g. 10 |

## Parameter `nptcls`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_oris | parameter input/output | num | Number of per-particle orientations | e.g. 10 |
| simulate_noise | parameter input/output | num | Number of particles | e.g. 10 |
| simulate_particles | parameter input/output | num | Number of particles | e.g. 10 |

## Parameter `nptcls_coarse`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Target # of particles per coarse chunk | e.g. 10 |
| sieve_cavgs | search controls | num | Target coarse-pass particle count | e.g. 10 |

## Parameter `nptcls_fine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Target # of particles per fine chunk | e.g. 10 |
| sieve_cavgs | search controls | num | Target fine-pass particle count | e.g. 10 |

## Parameter `nptcls_per_cls`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D_chunks | search controls | num | Number of particles per cluster | e.g. 10 |
| analysis2D_nano | search controls | num | Number of particles per cluster | e.g. 10 |
| check_refpick | search controls | num | Number of particles per class | e.g. 10 |
| cluster2D_nano | search controls | num | Number of particles per cluster | e.g. 10 |
| gen_pickrefs | search controls | num | Number of particles per cluster | e.g. 10 |
| mini_stream | search controls | num | Number of particles per class | e.g. 10 |
| sieve_cavgs | search controls | num | Number of particles per cluster | e.g. 10 |

## Parameter `nptcls_per_part`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| sample_classes | parameter input/output | num | Number of ptcls per part to select when balancing | e.g. 10 |
| selection | parameter input/output | num | Number of ptcls per part to select when balancing | e.g. 10 |

## Parameter `nran`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| selection | parameter input/output | num | Number of random samples | e.g. 10 |
| stackops | parameter input/output | num | Number of random samples | e.g. 10 |

## Parameter `nrestarts`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | search controls | num | Number of restarts | e.g. 10 |

## Parameter `nsample`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | num | Number of particles to sample | e.g. 10 |
| abinitio3D | search controls | num | Number of particles to sample | e.g. 10 |
| abinitio3D_nano | search controls | num | Number of particles to sample | e.g. 10 |
| refine3D_auto | search controls | num | Projection samples | e.g. 10 |
| refine3D_multi | search controls | num | Particle sample target | e.g. 10 |
| sample_classes | parameter input/output | num | Number of particles to sample | e.g. 10 |

## Parameter `nsample_coarse`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Number of particles to sample in coarse sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Coarse-pass sample count | e.g. 10 |

## Parameter `nsample_fine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | num | Number of particles to sample in fine sieving | e.g. 10 |
| sieve_cavgs | search controls | num | Fine-pass sample count | e.g. 10 |

## Parameter `nspace`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | search controls | num | Number of projection directions | e.g. 10 |
| denoise_project | search controls | num | Number of projection directions | e.g. 10 |
| flex_analysis | search controls | num | Number of projection directions | e.g. 10 |
| fsc_area_score | parameter input/output | num | Number of cone axes | e.g. 10 |
| make_cavgs | parameter input/output | num | Number of projection directions | e.g. 10 |
| oristats | parameter input/output | num | Number of projection directions | e.g. 10 |
| refine3D | search controls | num | Number of projection directions | e.g. 10 |
| refine3D_nano | search controls | num | Number of projection directions | e.g. 10 |
| reproject | search controls | num | Number of projection directions | e.g. 10 |
| trajectory_make_projavgs | parameter input/output | num | Number of projection directions | e.g. 10 |
| vizoris | parameter input/output | num | Number of projection directions | e.g. 10 |

## Parameter `nspace_sub`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| denoise_project | search controls | num | SO3 mixture subspace size | e.g. 10 |

## Parameter `nstages`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Last ab initio stage to run | e.g. 10 |
| estimate_lpstages | parameter input/output | num | Number of low-pass limit stages | e.g. 10 |

## Parameter `nstates`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Number of states | e.g. 10 |
| abinitio3D_cavgs | search controls | num | Number of states | e.g. 10 |
| bootstrap_rec3D | search controls | num | Number of states | e.g. 10 |
| make_oris | parameter input/output | num | Number of random state labels | e.g. 10 |
| noisevol | parameter input/output | num | Number states | e.g. 10 |
| oriops | parameter input/output | num | Number of random state labels | e.g. 10 |
| ptcl3D_state_consensus | search controls | num | Number of states | e.g. 10 |
| refine3D | search controls | num | Number of states | e.g. 10 |
| refine3D_multi | search controls | num | Number of states | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Number of states | e.g. 10 |

## Parameter `nsubcls_max`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | parameter input/output | num | Maximum subclass trial count in auto mode (default 10) | e.g. 10 |

## Parameter `nsubcls_min`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | parameter input/output | num | Minimum subclass trial count in auto mode (default 3) | e.g. 10 |

## Parameter `nthr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| abinitio2D_chunks | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| abinitio2D_stream | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| abinitio3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| abinitio3D_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| abinitio3D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| analysis2D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| assign_optics | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| atoms_register | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| atoms_stats | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| auto_spher_mask | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| automask | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| automask2D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| autorefine3D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| binarize | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| bootstrap_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| bootstrap_rec3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cavgseoproc_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cavgsproc_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| center | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| center2D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| check_refpick | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cls_split | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cluster2D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cluster_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cluster_cavgs_selection | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| cluster_stack | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| conv_atom_denoise | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| crys_score | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ctf_correct | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ctf_estimate | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ctfops | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| denoise_project | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| detect_atoms | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| dock_volpair | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| estimate_diam | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| extract | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| filter | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| flex_analysis | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| fractionate_movies | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| fsc | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| fsc_area_score | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| gen_pickrefs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| gen_pspecs_and_thumbs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| graphene_subtr | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| icm2D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| icm3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| make_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| make_oris | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| mask | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| match_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| match_stacks | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| mini_stream | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| model_cavgs_rejection | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| motion_correct | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| normalize | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| nu_filt3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| oristats | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| particle_sieving | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| pick | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| pick_extract | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| postprocess | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ppca_denoise | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ppca_denoise_classes | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ppca_volvar | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| preproc | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| preprocess | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| ptclsproc_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| reconstruct3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| reconstruct3D_pcg | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| reextract | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| refine3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| refine3D_auto | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| refine3D_multi | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| refine3D_nano | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| reproject | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| sample_classes | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| scale | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| sieve_cavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| simulate_movie | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| simulate_nanoparticle | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| simulate_particles | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| stackops | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| symaxis_search | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| symmetrize_map | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| symmetry_test | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| track_particles | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| trajectory_denoise | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| trajectory_make_projavgs | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| trajectory_reconstruct3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| tseries_extractor | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| tseries_make_pickavg | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| tseries_motion_correct | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| tseries_prep4tracking | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| uniform_filter2D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| uniform_filter3D | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| validate_cavgs_vs_model | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| volanalyze | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| volcluster | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |
| volops | computer controls | num | Number of threads per computing node, give 0 if unsure | e.g. 10 |

## Parameter `nthr_ini3D`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | computer controls | num | Number of threads for ini3D phase, give 0 if unsure | e.g. 10 |

## Parameter `nu_refine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D | filter controls | binary | NU resolution expansion refinement | *(empty)* |
| refine3D_auto | filter controls | binary | NU resolution expansion refinement | *(empty)* |

## Parameter `numlen`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| preprocess | parameter input/output | num | Length of number string | e.g. 10 |

## Parameter `nxpatch`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | search controls | num | # of patches along x-axis | e.g. 10 |
| preprocess | search controls | num | # of patches along x-axis | e.g. 10 |
| tseries_make_pickavg | search controls | num | # of patches along x-axis | e.g. 10 |
| tseries_motion_correct | search controls | num | # of patches along x-axis | e.g. 10 |

## Parameter `nypatch`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | search controls | num | # of patches along y-axis | e.g. 10 |
| preprocess | search controls | num | # of patches along y-axis | e.g. 10 |
| tseries_make_pickavg | search controls | num | # of patches along y-axis | e.g. 10 |
| tseries_motion_correct | search controls | num | # of patches along y-axis | e.g. 10 |

## Parameter `objfun`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| reconstruct3D_pcg | parameter input/output | multi | Objective function | *(empty)* |
| refine3D | search controls | multi | Objective function | *(empty)* |

## Parameter `objfun_den`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | binary | Denoised objective | *(empty)* |
| refine3D | search controls | binary | Denoised objective | *(empty)* |

## Parameter `objfun_den_w`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Denoised objective weight | e.g. 10 |
| refine3D | search controls | num | Denoised objective weight | e.g. 10 |

## Parameter `offset`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| track_particles | search controls | num | Shift half-width search bound | e.g. 10 |

## Parameter `optics_dir`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| gen_pickrefs | file input/output | dir | Target directory for optics import | e.g. optics_assignment |

## Parameter `optics_offset`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| assign_optics_groups | parameter input/output | num | Numbering offset to apply to optics groups | e.g. 10 |
| export_relion | parameter input/output | num | Offset to apply to optics group numbering | e.g. 10 |

## Parameter `oritab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| center | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| ctfops | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| import_particles | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| mask | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| oriops | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| oristats | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| reproject | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| simulate_particles | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| stackops | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |
| vizoris | file input/output | file | Orientation and CTF parameter file | e.g. oritab.txt |

## Parameter `oritab2`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oristats | file input/output | file | 2nd orientation and CTF parameter file | e.g. oritab2.txt |

## Parameter `oritype`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | parameter input/output | multi | Particle type to split | *(empty)* |
| flex_analysis | parameter input/output | str | Particle orientation segment | e.g. value |
| make_oris | parameter input/output | multi | Oritype segment in project | *(empty)* |
| oriops | parameter input/output | multi | Oritype segment in project | *(empty)* |
| oristats | parameter input/output | multi | Oritype segment in project | *(empty)* |
| print_project_field | parameter input/output | multi | Oritype segment in project | *(empty)* |
| reconstruct3D_pcg | parameter input/output | multi | Oritype segment in project | *(empty)* |
| reextract | parameter input/output | multi | Oritype segment in project | *(empty)* |
| replace_project_field | parameter input/output | multi | Oritype segment in project | *(empty)* |
| selection | parameter input/output | multi | Oritype segment in project | *(empty)* |
| vizoris | parameter input/output | multi | Oritype segment in project | *(empty)* |

## Parameter `osmpl_fac`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_cavgs | parameter input/output | num | Oversampling factor | e.g. 10 |

## Parameter `outfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_rec3D | file input/output | file | Resolution output prefix | e.g. resolution |
| center | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| dock_volpair | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| flex_analysis | file input/output | file | Discrete-state project | e.g. flex_cluster_states.simple |
| make_oris | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| mask | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| oriops | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| simulate_particles | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| stackops | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| volcluster | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |
| volops | file input/output | file | Output orientation and CTF parameter file | e.g. orientations.txt |

## Parameter `outside`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract | parameter input/output | binary | Extract outside stage boundaries | *(empty)* |
| reextract | parameter input/output | binary | Extract outside stage boundaries | *(empty)* |

## Parameter `outstk`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| convert | image input/output | file | Output stack name | e.g. output.mrcs |
| ctf_correct | image input/output | file | Output stack name | e.g. output.mrcs |
| ctfops | image input/output | file | Output stack name | e.g. output.mrcs |
| filter | image input/output | file | Output stack name | e.g. output.mrcs |
| graphene_subtr | image input/output | file | Output stack name | e.g. output.mrcs |
| ppca_denoise | image input/output | file | Output stack name | e.g. output.mrcs |
| ppca_volvar | image input/output | file | Output stack name | e.g. output.mrcs |
| reproject | image input/output | file | Output stack name | e.g. output.mrcs |
| scale | image input/output | file | Output stack name | e.g. output.mrcs |
| simulate_particles | image input/output | file | Output stack name | e.g. output.mrcs |
| stack | image input/output | file | Output stack name | e.g. output.mrcs |
| stackops | image input/output | file | Output stack name | e.g. output.mrcs |
| trajectory_denoise | image input/output | file | Output stack name | e.g. output.mrcs |

## Parameter `outvol`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| conv_atom_denoise | image input/output | file | Output volume name | e.g. volume.mrc |
| convert | image input/output | file | Output volume name | e.g. volume.mrc |
| dock_volpair | image input/output | file | Output volume name | e.g. volume.mrc |
| filter | image input/output | file | Output volume name | e.g. volume.mrc |
| flex_analysis | image input/output | file | Output volume name | e.g. volume.mrc |
| nu_filt3D | image input/output | file | Output volume name | e.g. volume.mrc |
| pdb2mrc | image input/output | file | Output volume name | e.g. volume.mrc |
| scale | image input/output | file | Output volume name | e.g. volume.mrc |
| simulate_nanoparticle | image input/output | file | Output volume name | e.g. volume.mrc |
| symmetrize_map | image input/output | file | Output volume name | e.g. volume.mrc |
| trajectory_reconstruct3D | image input/output | file | Output volume name | e.g. volume.mrc |
| volops | image input/output | file | Output volume name | e.g. volume.mrc |

## Parameter `overlap`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Convergence overlap target | e.g. 10 |
| abinitio3D_cavgs | search controls | num | Convergence overlap target | e.g. 10 |
| autorefine3D_nano | search controls | num | Convergence overlap target | e.g. 10 |
| refine3D_multi | search controls | num | State-overlap convergence target | e.g. 10 |

## Parameter `particle_density`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | multi | Particle density in micrographs | *(empty)* |
| pick | search controls | multi | Particle density in micrographs | *(empty)* |

## Parameter `pca_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cls_split | filter controls | multi | Class split embedding method | *(empty)* |
| ppca_denoise | filter controls | multi | PCA methods: PPCA, PPCA plus residual kPCA, standard SVD PCA, kernel PCA, or diffusion maps | *(empty)* |
| ppca_denoise_classes | filter controls | multi | PCA methods: PPCA, standard SVD PCA or kernel PCA | *(empty)* |
| trajectory_denoise | filter controls | multi | PCA methods: diffusion maps, PPCA, PPCA plus residual kPCA, standard SVD PCA, or kernel PCA | *(empty)* |

## Parameter `pcgop`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| reconstruct3D_pcg | search controls | multi | PCG normal operator | *(empty)* |

## Parameter `pcontrast`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | multi | Input particle contrast | *(empty)* |
| extract | parameter input/output | multi | Input particle contrast | *(empty)* |
| mini_stream | parameter input/output | multi | Input particle contrast | *(empty)* |
| pick | parameter input/output | multi | Input particle contrast | *(empty)* |
| pick_extract | parameter input/output | multi | Input particle contrast | *(empty)* |
| reextract | parameter input/output | multi | Input particle contrast | *(empty)* |

## Parameter `pdbfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_stats | file input/output | file | PDB | e.g. molecule.pdb |
| mask | mask controls | file | PDB for 3D envelope masking | e.g. molecule.pdb |
| model_validate | file input/output | file | PDB input coordinates file | e.g. molecule.pdb |
| pdb2mrc | file input/output | file | PDB input coordinates file | e.g. molecule.pdb |
| simulate_nanoparticle | file input/output | file | PDB | e.g. molecule.pdb |

## Parameter `pdbfile2`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_stats | file input/output | file | PDB | e.g. molecule.pdb |

## Parameter `pdbfiles`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_rmsd | file input/output | file | txt | e.g. pdb_files.txt |
| core_atoms_analysis | file input/output | file | txt | e.g. pdb_files.txt |

## Parameter `pdbout`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pdb2mrc | file input/output | file | Output PDB volume-centered molecule | e.g. output.pdb |
| simulate_nanoparticle | file input/output | file | Output PDB | e.g. output.pdb |

## Parameter `period`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| tseries_make_pickavg | parameter input/output | num | Period for repeated averaging windows | e.g. 10 |

## Parameter `pgrp`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | str | Point-group symmetry | e.g. value |
| abinitio3D_cavgs | search controls | str | Point-group symmetry | e.g. value |
| autorefine3D_nano | search controls | str | Point-group symmetry | e.g. value |
| bootstrap_rec3D | search controls | str | Point-group symmetry | e.g. value |
| cavgseoproc_nano | search controls | str | Point-group symmetry | e.g. value |
| cavgsproc_nano | search controls | str | Point-group symmetry | e.g. value |
| make_oris | parameter input/output | str | Point-group symmetry | e.g. value |
| oriops | parameter input/output | str | Point-group symmetry | e.g. value |
| oristats | parameter input/output | str | Point-group symmetry | e.g. value |
| pick_extract | search controls | str | Point-group symmetry | e.g. value |
| ptclsproc_nano | search controls | str | Point-group symmetry | e.g. value |
| reconstruct3D | search controls | str | Point-group symmetry | e.g. value |
| refine3D | search controls | str | Point-group symmetry | e.g. value |
| refine3D_auto | search controls | str | Point-group symmetry | e.g. value |
| refine3D_multi | search controls | str | Point-group symmetry | e.g. value |
| refine3D_nano | search controls | str | Point-group symmetry | e.g. value |
| reproject | search controls | str | Point-group symmetry | e.g. value |
| simulate_particles | search controls | str | Point-group symmetry | e.g. value |
| symaxis_search | search controls | str | Point-group symmetry | e.g. value |
| symmetrize_map | search controls | str | Point-group symmetry | e.g. value |
| trajectory_reconstruct3D | search controls | str | Point-group symmetry | e.g. value |
| validate_cavgs_vs_model | search controls | str | Point-group symmetry | e.g. value |
| vizoris | parameter input/output | str | Point-group symmetry | e.g. value |

## Parameter `pgrp_start`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | str | Initital point-group symmetry | e.g. value |
| abinitio3D_cavgs | search controls | str | Initital point-group symmetry | e.g. value |

## Parameter `phrand`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | binary | Phase randomization | *(empty)* |

## Parameter `phshift_max`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | Maximum CTF phase shift | e.g. 10 |
| ctf_estimate | search controls | num | Maximum CTF phase shift | e.g. 10 |
| mini_stream | search controls | num | Maximum CTF phase shift | e.g. 10 |
| preproc | search controls | num | Maximum CTF phase shift | e.g. 10 |
| preprocess | search controls | num | Maximum CTF phase shift | e.g. 10 |

## Parameter `phshift_min`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | Minimum CTF phase shift | e.g. 10 |
| ctf_estimate | search controls | num | Minimum CTF phase shift | e.g. 10 |
| mini_stream | search controls | num | Minimum CTF phase shift | e.g. 10 |
| preproc | search controls | num | Minimum CTF phase shift | e.g. 10 |
| preprocess | search controls | num | Minimum CTF phase shift | e.g. 10 |

## Parameter `phshift_step`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | num | CTF phase-shift step | e.g. 10 |
| ctf_estimate | search controls | num | CTF phase-shift step | e.g. 10 |
| mini_stream | search controls | num | CTF phase-shift step | e.g. 10 |
| preproc | search controls | num | CTF phase-shift step | e.g. 10 |
| preprocess | search controls | num | CTF phase-shift step | e.g. 10 |

## Parameter `phshiftunit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | parameter input/output | multi | Phase-shift unit | *(empty)* |

## Parameter `pick_roi`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | search controls | binary | Artefactual regions exclusion(new picker only) | *(empty)* |
| gen_pickrefs | search controls | binary | Artefactual regions exclusion(new picker only) | *(empty)* |
| mini_stream | search controls | binary | Artefactual regions exclusion(new picker only) | *(empty)* |
| pick | search controls | binary | Artefactual regions exclusion(new picker only) | *(empty)* |
| pick_extract | search controls | binary | Artefactual regions exclusion(new picker only) | *(empty)* |

## Parameter `picker`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pick | parameter input/output | multi | Which picker to use | *(empty)* |

## Parameter `pickrefs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | image input/output | file | Stack of class-averages/reprojections for picking | e.g. pickrefs.mrc |
| master | image input/output | file | 2D averages for use as picking references (optional) | e.g. pickrefs.mrc |
| pick | image input/output | file | Stack of class-averages/reprojections for picking | e.g. pickrefs.mrc |
| pick_extract | image input/output | file | Stack of class-averages/reprojections for picking | e.g. pickrefs.mrc |

## Parameter `plaintexttab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | file input/output | file | Plain text file of input parameters | e.g. params.txt |

## Parameter `platonic`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| symmetry_test | search controls | binary | Search for Platonic symmetries | *(empty)* |

## Parameter `plot_key`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| print_project_field | parameter input/output | string | plot plot_key on , sort on x | e.g. value |

## Parameter `positive`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask2D | mask controls | binary | Consider only positive pixels | *(empty)* |

## Parameter `postprocess`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_rec3D | filter controls | binary | Postprocess final map | *(empty)* |
| reconstruct3D | filter controls | binary | Postprocess final map | *(empty)* |

## Parameter `ppca_kpca_resid_alpha`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise | filter controls | num | Residual hybrid damping (0 => PPCA only; default 0.5) | e.g. 10 |
| trajectory_denoise | filter controls | num | Residual hybrid damping (0 => PPCA only; default 0.5) | e.g. 10 |

## Parameter `pre_norm`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise_classes | parameter input/output | binary | Pre-normalize images | *(empty)* |

## Parameter `preimage_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| flex_analysis | filter controls | multi | Diffusion-map output mode | *(empty)* |

## Parameter `preimage_ndim`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| flex_analysis | filter controls | num | Local-linear pre-image design dimension (default 2) | e.g. 10 |

## Parameter `prob_neigh_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D | search controls | multi | Prob-neigh neighborhood mode | *(empty)* |
| refine3D_multi | search controls | multi | Prob-neigh neighborhood mode | *(empty)* |

## Parameter `projfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| estimate_lpstages | file input/output | file | Project file | e.g. input.simple |
| export_manifoldem_starproject | file input/output | file | Project file | e.g. input.simple |
| extract_subproj | file input/output | file | Project file | e.g. input.simple |
| extract_substk | file input/output | file | Project file | e.g. input.simple |
| map_params_from_den | file input/output | file | Project file | e.g. input.simple |
| match_cavgs | file input/output | file | Project file | e.g. input.simple |
| new_project | file input/output | file | Project file | e.g. input.simple |
| ptcl3D_state_consensus | file input/output | file | Project file | e.g. input.simple |
| replace_project_field | file input/output | file | Project file | e.g. input.simple |
| scale | file input/output | file | Project file | e.g. input.simple |
| selection | file input/output | file | Project file | e.g. input.simple |
| update_project | file input/output | file | Project file | e.g. input.simple |
| validate_projfile | file input/output | file | Project file | e.g. input.simple |
| write_mic_filetab | file input/output | file | Project file | e.g. input.simple |

## Parameter `projfile_den`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| map_params_from_den | file input/output | file | Denoised child project file | e.g. input.simple |

## Parameter `projfile_merged`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| merge_projects | file input/output | file | Merged output project file | e.g. input.simple |

## Parameter `projfile_orig`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| unbootstrap_cavgs | file input/output | file | Original project file | e.g. input.simple |

## Parameter `projfile_raw`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| map_params_from_den | file input/output | file | Raw project file | e.g. input.simple |

## Parameter `projfile_ref`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| match_cavgs | file input/output | file | Reference project file | e.g. input.simple |

## Parameter `projfile_target`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| replace_project_field | file input/output | file | Another project file | e.g. input.simple |

## Parameter `projname`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | parameter input/output | str | Project name | e.g. value |

## Parameter `projrec`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | binary | Projection-direction reconstruction | *(empty)* |
| reconstruct3D | search controls | binary | Projection-direction reconstruction | *(empty)* |
| refine3D | search controls | binary | Projection-direction reconstruction | *(empty)* |

## Parameter `projstats`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oristats | parameter input/output | binary | Projection statistics | *(empty)* |

## Parameter `projtab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| merge_projects | file input/output | file | Project file table | e.g. projtab.txt |
| ptcl3D_state_consensus | file input/output | file | Project file table | e.g. projtab.txt |

## Parameter `prune`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cluster_cavgs | parameter input/output | binary | Automated particles pruning | *(empty)* |
| map_cavgs_selection | parameter input/output | binary | Automated particles pruning | *(empty)* |
| match_cavgs | parameter input/output | binary | Automated particles pruning | *(empty)* |
| model_cavgs_rejection | parameter input/output | binary | Automated particles pruning | *(empty)* |
| ptcl3D_state_consensus | parameter input/output | binary | Automated particles pruning | *(empty)* |
| select_clusters | parameter input/output | binary | Automated particles pruning | *(empty)* |
| selection | parameter input/output | binary | Automated particles pruning | *(empty)* |

## Parameter `pspecsz`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| check_refpick | parameter input/output | num | Size of power spectrum | e.g. 10 |
| ctf_estimate | parameter input/output | num | Size of power spectrum | e.g. 10 |
| gen_pspecs_and_thumbs | parameter input/output | num | Size of power spectrum | e.g. 10 |
| mini_stream | parameter input/output | num | Size of power spectrum | e.g. 10 |
| motion_correct | parameter input/output | num | Size of power spectrum | e.g. 10 |
| preproc | parameter input/output | num | Size of power spectrum | e.g. 10 |
| preprocess | parameter input/output | num | Size of power spectrum | e.g. 10 |

## Parameter `ptcl_src`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | multi | Particle source | *(empty)* |
| reconstruct3D | search controls | multi | Particle source | *(empty)* |
| refine3D | search controls | multi | Particle source | *(empty)* |
| refine3D_auto | search controls | multi | Particle source | *(empty)* |
| refine3D_multi | search controls | multi | Particle source | *(empty)* |

## Parameter `qsys_name`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | multi | Queue system kind | *(empty)* |
| update_project | computer controls | multi | Queue system kind | *(empty)* |

## Parameter `qsys_partition`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | Name of SLURM/PBS/LSF partition | e.g. value |
| update_project | computer controls | str | Name of SLURM/PBS/LSF partition | e.g. value |

## Parameter `qsys_qos`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | Schedule priority | e.g. value |
| update_project | computer controls | str | Schedule priority | e.g. value |

## Parameter `qsys_reservation`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | Name of reserved partition | e.g. value |
| update_project | computer controls | str | Name of reserved partition | e.g. value |

## Parameter `quality_context`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| model_cavgs_rejection | parameter input/output | multi | Class-average quality context | *(empty)* |

## Parameter `quality_mode`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| model_cavgs_rejection | parameter input/output | multi | Class-average quality mode | *(empty)* |

## Parameter `quality_model`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| model_cavgs_rejection | parameter input/output | multi | Class-average quality model | *(empty)* |

## Parameter `real_filter`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | multi | Real-space filter | *(empty)* |

## Parameter `ref_ind`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| volanalyze | parameter input/output | num | Reference volume index | e.g. 10 |

## Parameter `refine`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | multi | Refinement mode | *(empty)* |
| abinitio2D_chunks | search controls | multi | Refinement mode | *(empty)* |
| refine3D | search controls | multi | Refinement mode | *(empty)* |
| refine3D_nano | search controls | multi | Refinement mode | *(empty)* |

## Parameter `refs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_cavgs | image input/output | file | Bootstrap reference class averages | e.g. references.mrc |
| make_cavgs | image input/output | file | Output 2D references | e.g. references.mrc |
| particle_sieving | search controls | file | Reference class averages to initialise size compatibility model | e.g. references.mrc |

## Parameter `reliongroups`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| export_relion | parameter input/output | num | Number of Relion groups based on defocus | e.g. 10 |

## Parameter `remap_cls`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_cavgs | parameter input/output | binary | Whether to remap 2D clusters | *(empty)* |

## Parameter `res_target`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D_auto | filter controls | num | Resolution target (in A) | e.g. 10 |

## Parameter `res_threshold`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| selection | parameter input/output | num | Class resolution threshold | e.g. 10 |

## Parameter `rmsd_file`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_stats | file input/output | file | bin | e.g. atom_rmsd.bin |

## Parameter `roavg`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| estimate_diam | parameter input/output | binary | Rotationally average | *(empty)* |
| stackops | parameter input/output | binary | Rotationally average | *(empty)* |

## Parameter `rtol`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| reconstruct3D_pcg | filter controls | num | PCG relative residual tolerance | e.g. 10 |

## Parameter `scale`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| scale | parameter input/output | num | Scaling ratio | e.g. 10 |

## Parameter `script`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| analysis2D_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| autorefine3D_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| cavgseoproc_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| cavgsproc_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| center2D_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| cluster2D_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| ptclsproc_nano | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |
| validate_cavgs_vs_model | computer controls | binary | Generate script for shared-mem exec on cluster | *(empty)* |

## Parameter `select_flag`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| select_clusters | parameter input/output | multi | flag to use for selection | *(empty)* |

## Parameter `sherr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| make_oris | parameter input/output | num | Shift error half-width | e.g. 10 |
| oriops | parameter input/output | num | Shift error half-width | e.g. 10 |
| simulate_particles | parameter input/output | num | Shift error half-width | e.g. 10 |

## Parameter `sigma`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | num | sigma, for Gaussian generation | e.g. 10 |

## Parameter `sigma_est`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | multi | Sigma estimation method | *(empty)* |
| refine3D | search controls | multi | Sigma estimation method | *(empty)* |
| refine3D_auto | search controls | multi | Sigma estimation method | *(empty)* |
| refine3D_multi | search controls | multi | Sigma estimation method | *(empty)* |
| refine3D_nano | search controls | multi | Sigma estimation method | *(empty)* |

## Parameter `single_pass`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | str | Coarse pass only | e.g. value |
| sieve_cavgs | search controls | binary | Run coarse sieving pass only | *(empty)* |

## Parameter `skip_rejection`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| cluster_cavgs | search controls | binary | Skip class-average rejection | *(empty)* |

## Parameter `smpd`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_rmsd | parameter input/output | num | Sampling distance | e.g. 10 |
| atoms_stats | parameter input/output | num | Sampling distance | e.g. 10 |
| auto_spher_mask | parameter input/output | num | Sampling distance | e.g. 10 |
| automask | parameter input/output | num | Sampling distance | e.g. 10 |
| automask2D | parameter input/output | num | Sampling distance | e.g. 10 |
| autorefine3D_nano | parameter input/output | num | Sampling distance | e.g. 10 |
| cavgseoproc_nano | parameter input/output | num | Sampling distance | e.g. 10 |
| cavgsproc_nano | parameter input/output | num | Sampling distance | e.g. 10 |
| center | parameter input/output | num | Sampling distance | e.g. 10 |
| check_refpick | parameter input/output | num | Sampling distance | e.g. 10 |
| cif2mrc | parameter input/output | num | Sampling distance | e.g. 10 |
| conv_atom_denoise | parameter input/output | num | Sampling distance | e.g. 10 |
| core_atoms_analysis | parameter input/output | num | Sampling distance | e.g. 10 |
| crys_score | parameter input/output | num | Sampling distance | e.g. 10 |
| ctf_correct | parameter input/output | num | Sampling distance | e.g. 10 |
| ctfops | parameter input/output | num | Sampling distance | e.g. 10 |
| detect_atoms | parameter input/output | num | Sampling distance | e.g. 10 |
| dock_volpair | parameter input/output | num | Sampling distance | e.g. 10 |
| estimate_diam | parameter input/output | num | Sampling distance | e.g. 10 |
| filter | parameter input/output | num | Sampling distance | e.g. 10 |
| fsc | parameter input/output | num | Sampling distance | e.g. 10 |
| fsc_area_score | parameter input/output | num | Sampling distance | e.g. 10 |
| graphene_subtr | parameter input/output | num | Sampling distance | e.g. 10 |
| icm2D | parameter input/output | num | Sampling distance | e.g. 10 |
| icm3D | parameter input/output | num | Sampling distance | e.g. 10 |
| import_cavgs | parameter input/output | num | Sampling distance | e.g. 10 |
| import_movies | parameter input/output | num | Sampling distance | e.g. 10 |
| import_particles | parameter input/output | num | Sampling distance | e.g. 10 |
| import_trajectory | parameter input/output | num | Sampling distance | e.g. 10 |
| mask | parameter input/output | num | Sampling distance | e.g. 10 |
| master | parameter input/output | float | Pixel size (A) | e.g. 10 |
| mini_stream | parameter input/output | num | Sampling distance | e.g. 10 |
| model_validate | parameter input/output | num | Sampling distance | e.g. 10 |
| noisevol | parameter input/output | num | Sampling distance | e.g. 10 |
| normalize | parameter input/output | num | Sampling distance | e.g. 10 |
| nu_filt3D | parameter input/output | num | Sampling distance | e.g. 10 |
| pdb2mrc | parameter input/output | num | Sampling distance | e.g. 10 |
| ppca_denoise | parameter input/output | num | Sampling distance | e.g. 10 |
| ppca_volvar | parameter input/output | num | Sampling distance | e.g. 10 |
| preproc | parameter input/output | num | Sampling distance | e.g. 10 |
| print_dose_weights | parameter input/output | num | Sampling distance | e.g. 10 |
| print_fsc | parameter input/output | num | Sampling distance | e.g. 10 |
| print_magic_boxes | parameter input/output | num | Sampling distance | e.g. 10 |
| ptclsproc_nano | parameter input/output | num | Sampling distance | e.g. 10 |
| reproject | parameter input/output | num | Sampling distance | e.g. 10 |
| scale | parameter input/output | num | Sampling distance | e.g. 10 |
| simulate_movie | parameter input/output | num | Sampling distance | e.g. 10 |
| simulate_nanoparticle | parameter input/output | num | Sampling distance | e.g. 10 |
| simulate_particles | parameter input/output | num | Sampling distance | e.g. 10 |
| split | parameter input/output | num | Sampling distance | e.g. 10 |
| stack | parameter input/output | num | Sampling distance | e.g. 10 |
| stackops | parameter input/output | num | Sampling distance | e.g. 10 |
| symaxis_search | parameter input/output | num | Sampling distance | e.g. 10 |
| symmetrize_map | parameter input/output | num | Sampling distance | e.g. 10 |
| symmetry_test | parameter input/output | num | Sampling distance | e.g. 10 |
| trajectory_denoise | parameter input/output | num | Sampling distance | e.g. 10 |
| tsegmaps_core_finder | parameter input/output | num | Sampling distance | e.g. 10 |
| tseries_import | parameter input/output | num | Sampling distance | e.g. 10 |
| uniform_filter2D | parameter input/output | num | Sampling distance | e.g. 10 |
| uniform_filter3D | parameter input/output | num | Sampling distance | e.g. 10 |
| validate_cavgs_vs_model | parameter input/output | num | Sampling distance | e.g. 10 |
| volanalyze | parameter input/output | num | Sampling distance | e.g. 10 |
| volcluster | parameter input/output | num | Sampling distance | e.g. 10 |
| volops | parameter input/output | num | Sampling distance | e.g. 10 |

## Parameter `smpd_downscale`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | parameter input/output | hidden_float | Downscaled pixel size (A) | *(empty)* |
| motion_correct | parameter input/output | num | Sampling distance after downscale | e.g. 10 |
| preproc | parameter input/output | num | Sampling distance after downscale | e.g. 10 |
| preprocess | parameter input/output | num | Sampling distance after downscale | e.g. 10 |

## Parameter `smpd_target`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| model_validate | parameter input/output | num | Target sampling distance | e.g. 10 |

## Parameter `snr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_movie | parameter input/output | num | SNR | e.g. 10 |
| simulate_particles | parameter input/output | num | SNR | e.g. 10 |
| stackops | parameter input/output | num | Apply noise to give SNR | e.g. 10 |
| volops | parameter input/output | num | SNR | e.g. 10 |

## Parameter `sort`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| print_project_field | parameter input/output | string | sort oris on key | e.g. value |

## Parameter `split_stage`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Docked-state split stage | e.g. 10 |
| abinitio3D_cavgs | search controls | num | Docked-state split stage | e.g. 10 |

## Parameter `starfile`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| export_manifoldem_starproject | file input/output | file | STAR-format file name | e.g. input.star |
| import_particles | file input/output | file | Particles Metadata starfile | e.g. input.star |
| import_starproject | file input/output | file | Metadata starfile | e.g. input.star |

## Parameter `state`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | search controls | num | Continuation state label | e.g. 10 |
| extract_substk | parameter input/output | num | State index | e.g. 10 |
| oriops | parameter input/output | num | State to modify | e.g. 10 |
| postprocess | parameter input/output | num | State to postprocess | e.g. 10 |
| prune_project | parameter input/output | num | State index | e.g. 10 |
| reconstruct3D_pcg | parameter input/output | num | State to reconstruct | e.g. 10 |
| selection | parameter input/output | num | State number | e.g. 10 |
| stackops | parameter input/output | num | State index | e.g. 10 |

## Parameter `stats`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D | search controls | binary | Write class-average statistics | *(empty)* |
| info_image | parameter input/output | binary | Output statistics | *(empty)* |
| stackops | parameter input/output | binary | Provide statistics | *(empty)* |

## Parameter `stepsz`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| trajectory_reconstruct3D | parameter input/output | num | Time window size (# frames){500} | e.g. 10 |

## Parameter `stk`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask2D | image input/output | file | Particle image stack | e.g. particles.mrcs |
| binarize | image input/output | file | Stack | e.g. particles.mrcs |
| center | image input/output | file | Particle image stack | e.g. particles.mrcs |
| cluster_cavgs | image input/output | file | Particle image stack | e.g. particles.mrcs |
| cluster_stack | image input/output | file | Particle image stack | e.g. particles.mrcs |
| convert | image input/output | file | Stack | e.g. particles.mrcs |
| ctfops | image input/output | file | Particle image stack | e.g. particles.mrcs |
| estimate_diam | image input/output | file | Particle image stack | e.g. particles.mrcs |
| filter | image input/output | file | Stack to filter | e.g. particles.mrcs |
| icm2D | image input/output | file | Odd stack | e.g. particles.mrcs |
| import_cavgs | image input/output | file | Stack of class averages | e.g. particles.mrcs |
| import_particles | image input/output | file | Stack of particles | e.g. particles.mrcs |
| import_trajectory | image input/output | file | Particle image stack | e.g. particles.mrcs |
| map_cavgs_selection | image input/output | file | Stack of cavgs to select from | e.g. particles.mrcs |
| mask | image input/output | file | Particle image stack | e.g. particles.mrcs |
| match_stacks | image input/output | file | Particle image stack | e.g. particles.mrcs |
| normalize | image input/output | file | Stack to normalize | e.g. particles.mrcs |
| ppca_denoise | image input/output | file | Stack to denoise | e.g. particles.mrcs |
| reimport_particles | image input/output | file | Denoised particle stack | e.g. particles.mrcs |
| scale | image input/output | file | Particle image stack | e.g. particles.mrcs |
| simulate_movie | image input/output | file | Particle image stack | e.g. particles.mrcs |
| split | image input/output | file | Particle image stack | e.g. particles.mrcs |
| stackops | image input/output | file | Particle image stack | e.g. particles.mrcs |
| trajectory_denoise | image input/output | file | Stack to denoise | e.g. particles.mrcs |
| trajectory_swap_stack | image input/output | file | Particle image stack | e.g. particles.mrcs |
| uniform_filter2D | image input/output | file | Odd stack | e.g. particles.mrcs |

## Parameter `stk2`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| icm2D | image input/output | file | Even stack | e.g. particles2.mrcs |
| map_cavgs_selection | image input/output | file | Stack of selected cavgs | e.g. particles2.mrcs |
| match_stacks | image input/output | file | Second Particle image stack | e.g. particles2.mrcs |
| uniform_filter2D | image input/output | file | Even stack | e.g. particles2.mrcs |

## Parameter `stk_backgr`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| graphene_subtr | image input/output | file | background power spectra stack, eg NP_X_background_pspec.mrc | e.g. background_pspec.mrcs |

## Parameter `stk_den`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | image input/output | file | Denoised particle image stack | e.g. denoised.mrcs |

## Parameter `stk_traj`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| graphene_subtr | image input/output | file | Tracked NP image stack | e.g. trajectory.mrcs |

## Parameter `stktab`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | file input/output | file | List of per-micrograph particle stacks | e.g. stktab.txt |
| info_stktab | file input/output | file | List of per-micrograph particle stacks | e.g. stktab.txt |

## Parameter `stktab_den`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| import_particles | file input/output | file | List of denoised particle stacks | e.g. stktab_den.txt |

## Parameter `subprojname`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract_subproj | parameter input/output | str | Subproject name | e.g. value |

## Parameter `symrnd`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | binary | Randomize over subgroubs of point-group | *(empty)* |

## Parameter `taper_edges`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| mask | mask controls | binary | Taper edges | *(empty)* |

## Parameter `tester`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| track_particles | parameter input/output | multi | Write periodic tester-mode outputs | *(empty)* |

## Parameter `thres`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask | mask controls | num | Volume threshold | e.g. 10 |
| master | parameter input/output | hidden_float | Distance threshold for peak picking(A) | *(empty)* |
| pick_extract | search controls | num | Peak-picking distance threshold | e.g. 10 |

## Parameter `tilt_thres`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| assign_optics_groups | parameter input/output | num | Threshold for hierarchical clustering of beamtilts | e.g. 10 |
| export_relion | parameter input/output | num | Distance threshold | e.g. 10 |
| preproc | search controls | num | Threshold for hierarchical clustering of beamtilts | e.g. 10 |

## Parameter `tiltgroupmax`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| export_relion | parameter input/output | num | Max movies in a tilt/optics group | e.g. 10 |

## Parameter `time_per_image`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | num | Time per image | e.g. 10 |
| update_project | computer controls | num | Time per image | e.g. 10 |

## Parameter `tof`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| fractionate_movies | parameter input/output | num | Final frame | e.g. 10 |

## Parameter `top`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| extract_subproj | parameter input/output | num | To index | e.g. 10 |
| extract_substk | parameter input/output | num | To index | e.g. 10 |
| scale | parameter input/output | num | Last stack index | e.g. 10 |
| stackops | parameter input/output | num | To particle index | e.g. 10 |
| trajectory_reconstruct3D | parameter input/output | num | To particle index | e.g. 10 |

## Parameter `total_dose`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| master | parameter input/output | float | Total exposure dose (e/A2) | e.g. 10 |
| motion_correct | parameter input/output | num | Total exposure dose (e/Ang^2) | e.g. 10 |
| preproc | parameter input/output | num | Total exposure dose (e/Ang^2) | e.g. 10 |
| preprocess | parameter input/output | num | Total exposure dose (e/Ang^2) | e.g. 10 |
| print_dose_weights | parameter input/output | num | Total exposure dose (e/Ang^2) | e.g. 10 |

## Parameter `transp_pca`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| ppca_denoise_classes | filter controls | binary | transpose for pixel-wise learning | *(empty)* |

## Parameter `trs`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| autorefine3D_nano | search controls | num | Maximum translational shift | e.g. 10 |
| center2D_nano | search controls | num | Maximum translational shift | e.g. 10 |
| cluster2D_nano | search controls | num | Maximum translational shift | e.g. 10 |
| cluster_cavgs | search controls | num | Maximum translational shift | e.g. 10 |
| cluster_cavgs_selection | search controls | num | Maximum translational shift | e.g. 10 |
| cluster_stack | search controls | num | Maximum translational shift | e.g. 10 |
| dock_volpair | search controls | num | Maximum translational shift | e.g. 10 |
| match_cavgs | search controls | num | Maximum translational shift | e.g. 10 |
| match_stacks | search controls | num | Maximum translational shift | e.g. 10 |
| motion_correct | search controls | num | Max shift per iter in pixels{10.} | e.g. 10 |
| preproc | search controls | num | Max shift per iter in pixels{10.} | e.g. 10 |
| preprocess | search controls | num | Max shift per iter in pixels{10.} | e.g. 10 |
| reconstruct3D | search controls | num | Maximum translational shift | e.g. 10 |
| refine3D | search controls | num | Maximum translational shift | e.g. 10 |
| refine3D_nano | search controls | num | Maximum translational shift | e.g. 10 |
| simulate_movie | parameter input/output | num | Maximum translational shift | e.g. 10 |
| trajectory_reconstruct3D | search controls | num | Maximum translational shift | e.g. 10 |
| tseries_make_pickavg | search controls | num | Max shift per iter in pixels{10.} | e.g. 10 |
| tseries_motion_correct | search controls | num | Max shift per iter in pixels{10.} | e.g. 10 |

## Parameter `trsstats`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oristats | parameter input/output | binary | Shift statistics | *(empty)* |

## Parameter `tseries`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| vizoris | parameter input/output | binary | Time series analysis | *(empty)* |

## Parameter `update_frac`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D | search controls | num | Fractional update per iteration | e.g. 10 |
| refine3D_nano | search controls | num | Fractional update per iteration | e.g. 10 |

## Parameter `use_model`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| particle_sieving | search controls | str | Use model for class-average rejection in sieving | e.g. value |
| sieve_cavgs | search controls | binary | Use class-average rejection model | *(empty)* |

## Parameter `user_account`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | User account name in SLURM/PBS/LSF | e.g. value |
| update_project | computer controls | str | User account name in SLURM/PBS/LSF | e.g. value |

## Parameter `user_email`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | parameter input/output | str | Your e-mail address | e.g. value |
| update_project | parameter input/output | str | Your e-mail address | e.g. value |

## Parameter `user_project`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| new_project | computer controls | str | User project name in SLURM/PBS/LSF | e.g. value |
| update_project | computer controls | str | User project name in SLURM/PBS/LSF | e.g. value |

## Parameter `view_balance`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| flex_analysis | filter controls | binary | Correct uneven view distribution | *(empty)* |

## Parameter `vis`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| info_image | parameter input/output | binary | Visualize image | *(empty)* |
| stackops | parameter input/output | binary | Visualize stack images | *(empty)* |

## Parameter `vol1`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio3D | image input/output | file | Starting template volume | e.g. volume.mrc |
| atoms_stats | image input/output | file | Raw volume | e.g. volume.mrc |
| auto_spher_mask | image input/output | file | Odd volume | e.g. volume.mrc |
| automask | image input/output | file | Odd volume | e.g. volume.mrc |
| autorefine3D_nano | image input/output | file | FCC reference volume | e.g. volume.mrc |
| binarize | image input/output | file | Volume | e.g. volume.mrc |
| cavgseoproc_nano | image input/output | file | Volume | e.g. volume.mrc |
| cavgsproc_nano | image input/output | file | Volume | e.g. volume.mrc |
| center | image input/output | file | Volume | e.g. volume.mrc |
| conv_atom_denoise | image input/output | file | Volume | e.g. volume.mrc |
| convert | image input/output | file | Volume | e.g. volume.mrc |
| detect_atoms | image input/output | file | Volume | e.g. volume.mrc |
| dock_volpair | image input/output | file | Volume | e.g. volume.mrc |
| filter | image input/output | file | Volume to filter | e.g. volume.mrc |
| flex_analysis | image input/output | file | Mean volume | e.g. volume.mrc |
| fsc | image input/output | file | Odd volume | e.g. volume.mrc |
| fsc_area_score | image input/output | file | Odd volume | e.g. volume.mrc |
| icm3D | image input/output | file | Odd volume | e.g. volume.mrc |
| mask | image input/output | file | Volume | e.g. volume.mrc |
| model_validate | image input/output | file | Experimental volume | e.g. volume.mrc |
| normalize | image input/output | file | Volume to normalize | e.g. volume.mrc |
| nu_filt3D | image input/output | file | Odd volume | e.g. volume.mrc |
| postprocess | image input/output | file | Volume override | e.g. volume.mrc |
| ppca_volvar | image input/output | file | Volume | e.g. volume.mrc |
| refine3D | image input/output | file | Reference volume | e.g. volume.mrc |
| refine3D_auto | image input/output | file | Starting template volume | e.g. volume.mrc |
| refine3D_nano | image input/output | file | FCC reference volume | e.g. volume.mrc |
| reproject | image input/output | file | Volume | e.g. volume.mrc |
| scale | image input/output | file | Input volume | e.g. volume.mrc |
| simulate_particles | image input/output | file | Volume | e.g. volume.mrc |
| symaxis_search | image input/output | file | C1 Volume to identify symmetry axis of | e.g. volume.mrc |
| symmetrize_map | image input/output | file | Volume to symmetrize | e.g. volume.mrc |
| symmetry_test | image input/output | file | C1 Volume to identify symmetry of | e.g. volume.mrc |
| trajectory_reconstruct3D | image input/output | file | Mean volume for latent chunking | e.g. volume.mrc |
| uniform_filter3D | image input/output | file | Odd volume | e.g. volume.mrc |
| validate_cavgs_vs_model | image input/output | file | Volume | e.g. volume.mrc |
| volops | image input/output | file | Volume | e.g. volume.mrc |

## Parameter `vol2`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_stats | image input/output | file | Connected components volume | e.g. volume2.mrc |
| automask | image input/output | file | Even volume | e.g. volume2.mrc |
| dock_volpair | image input/output | file | Volume | e.g. volume2.mrc |
| fsc | image input/output | file | Even volume | e.g. volume2.mrc |
| fsc_area_score | image input/output | file | Even volume | e.g. volume2.mrc |
| icm3D | image input/output | file | Even volume | e.g. volume2.mrc |
| nu_filt3D | image input/output | file | Even volume | e.g. volume2.mrc |
| uniform_filter3D | image input/output | file | Even volume | e.g. volume2.mrc |

## Parameter `vol3`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| atoms_stats | image input/output | file | Volume | e.g. volume3.mrc |

## Parameter `vol_dim`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| pdb2mrc | parameter input/output | num | Simulated volume dimensions | e.g. 10 |

## Parameter `vol_even`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D_nano | image input/output | file | Even volume | e.g. even.mrc |

## Parameter `vol_odd`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| refine3D_nano | image input/output | file | Odd volume | e.g. odd.mrc |

## Parameter `walltime`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| abinitio2D_chunks | computer controls | num | Walltime | e.g. 10 |
| abinitio2D_stream | computer controls | num | Walltime | e.g. 10 |
| new_project | computer controls | num | Walltime | e.g. 10 |
| particle_sieving | computer controls | num | Walltime | e.g. 10 |
| pick_extract | computer controls | num | Walltime | e.g. 10 |
| preproc | computer controls | num | Walltime | e.g. 10 |
| sieve_cavgs | computer controls | num | Walltime | e.g. 10 |
| update_project | computer controls | num | Walltime | e.g. 10 |

## Parameter `wcrit`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| motion_correct | filter controls | multi | Correlation to weights conversion scheme | *(empty)* |
| tseries_make_pickavg | filter controls | multi | Correlation to weights conversion scheme | *(empty)* |
| tseries_motion_correct | filter controls | multi | Correlation to weights conversion scheme | *(empty)* |

## Parameter `which_iter`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| bootstrap_rec3D | parameter input/output | num | Sigma iteration index | e.g. 10 |

## Parameter `width`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| filter | filter controls | num | Cosine low-pass filter falloff | e.g. 10 |
| mask | mask controls | num | Falloff of inner mask | e.g. 10 |

## Parameter `winsz`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| automask2D | filter controls | num | Window size for median filter | e.g. 10 |
| binarize | parameter input/output | num | Half-window size | e.g. 10 |
| cluster2D_nano | search controls | num | Half-window size | e.g. 10 |
| filter | filter controls | num | Half-window size | e.g. 10 |
| gen_pickrefs | mask controls | num | Automask window size | e.g. 10 |
| pick | search controls | num | Window size for sauvola | e.g. 10 |
| stackops | parameter input/output | num | Window size for local sdev estimation | e.g. 10 |

## Parameter `xdim`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_movie | parameter input/output | num | x-dimension | e.g. 10 |

## Parameter `xmldir`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| assign_optics_groups | file input/output | dir | Directory containing per movie EPU XML files | e.g. /data/datasetid/xml |

## Parameter `xmlloc`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| export_relion | file input/output | file | Pathname of EPU XML files | e.g. /data/xml |

## Parameter `xsh`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| volops | parameter input/output | num | Translation along x-axis | e.g. 10 |

## Parameter `ydim`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| simulate_movie | parameter input/output | num | y-dimension | e.g. 10 |

## Parameter `ysh`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| volops | parameter input/output | num | Translation along y-axis | e.g. 10 |

## Parameter `zero`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| oriops | parameter input/output | binary | Zero shifts | *(empty)* |

## Parameter `zsh`

| Used by program | Section | Type | Parameter label | Placeholder |
| --- | --- | --- | --- | --- |
| volops | parameter input/output | num | Translation along z-axis | e.g. 10 |
