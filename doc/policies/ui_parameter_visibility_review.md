# Parameter Instance Visibility Review

This document inventories parameter visibility, not program visibility. Program visibility is intentionally excluded.

Each heading is one unique parameter key. Every row beneath it is a separate `ui_program_input` instance of that parameter in the named program. The value in **Parameter visibility** belongs only to that parameter instance and may differ between program contexts.

Edit **Parameter visibility** to `standard`, `advanced`, or `developer`. Required parameter instances must remain `standard`; optional instances may use any of the three values.

The initial parameter-instance assignments were inferred from the code:

1. Required parameter instances are `standard`.
2. Previously explicit parameter-instance visibility choices were preserved.
3. Otherwise-unclassified optional parameter instances in Developer programs were provisionally assigned `developer`.
4. Other otherwise-unclassified optional parameter instances were provisionally assigned `advanced`.

Inventory: 337 unique parameter keys and 1284 parameter instances. Parameter visibility totals: 455 Standard, 683 Advanced, and 146 Developer. 92 parameter keys currently have context-dependent visibility.

## Parameter `acf`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| stackops | Autocorrelation, A * conjg(A) | no | advanced |

## Parameter `algorithm`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | Algorithm for motion correction | no | advanced |
| preprocess | Algorithm for motion correction | no | advanced |
| tseries_motion_correct | Algorithm for motion correction | no | advanced |

## Parameter `amsklp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| auto_spher_mask | Low-pass limit for envelope mask generation | yes | standard |
| automask | Low-pass limit for envelope mask generation | no | advanced |
| automask2D | Low-pass limit for envelope mask generation | no | advanced |
| refine3D | Low-pass limit for envelope mask generation | no | advanced |
| refine3D_auto | Low-pass limit for envelope mask generation | no | advanced |

## Parameter `angastunit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | Angle of astigmatism unit | no | advanced |

## Parameter `angerr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_oris | Rotation angle error half-width | no | developer |
| oriops | Rotation angle error half-width | no | developer |

## Parameter `angstep`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| stackops | Angular stepsize | no | advanced |

## Parameter `append`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| selection | Append selection to existing | no | advanced |

## Parameter `astigerr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_particles | Astigmatism error | no | advanced |

## Parameter `astigthreshold`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick_extract | Astigmatism rejection threshold | no | advanced |

## Parameter `astigtol`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_estimate | Expected astigmatism | no | advanced |
| preprocess | Expected astigmatism | no | advanced |

## Parameter `athres`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fsc_area_score | Cone half-angle | no | advanced |
| trajectory_make_projavgs | Angular threshold (degrees) | no | developer |

## Parameter `automsk`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Perform envelope masking | no | standard |
| automask | Perform envelope masking | no | advanced |
| fsc_area_score | Perform envelope masking | no | advanced |
| nu_filt3D | Perform envelope masking | no | advanced |
| refine3D | Perform envelope masking | no | advanced |
| refine3D_multi | Perform envelope masking | no | advanced |

## Parameter `autoscale`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Automatic down-scaling | no | advanced |
| refine3D_auto | Automatic down-scaling | no | advanced |
| refine3D_multi | Automatic down-scaling | no | advanced |

## Parameter `avg`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| stackops | Average stack | no | advanced |

## Parameter `backgr_subtr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract | Perform micrograph background subtraction(new picker only) | no | advanced |
| pick | Perform micrograph background subtraction(new picker only) | no | advanced |
| reextract | Perform micrograph background subtraction(new picker only) | no | advanced |

## Parameter `balance`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| selection | Balanced selection of particles across classes | no | advanced |

## Parameter `bandwidth_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Diffusion-map bandwidth mode | no | advanced |
| denoise_project | Diffusion-map bandwidth mode | no | advanced |
| flex_analysis | Diffusion-map bandwidth mode | no | advanced |
| ppca_denoise | Diffusion-map bandwidth mode | no | advanced |

## Parameter `bandwidth_tune`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Ferguson bandwidth multiplier (default 1) | no | advanced |
| denoise_project | Ferguson bandwidth multiplier (default 1) | no | advanced |
| flex_analysis | Ferguson bandwidth multiplier (default 1) | no | advanced |
| ppca_denoise | Ferguson bandwidth multiplier (default 1) | no | advanced |

## Parameter `beamtilt`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| assign_optics_groups | Use beamtilts in optics group assignment | no | developer |
| preproc | Use beamtilts in optics group assignment | no | developer |

## Parameter `bfac`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctfops | CTF B-factor | no | advanced |
| filter | B-factor of Gaussian low-/high-pass filter | no | advanced |
| motion_correct | B-factor applied to frames | no | advanced |
| postprocess | B-factor for sharpening | no | advanced |
| preprocess | B-factor applied to frames | no | advanced |
| simulate_movie | CTF B-factor | no | advanced |
| simulate_particles | CTF B-factor | no | advanced |
| tseries_make_pickavg | B-factor applied to frames | no | advanced |
| tseries_motion_correct | B-factor applied to frames | no | advanced |
| volops | B-factor for sharpening | no | advanced |

## Parameter `bfacerr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_particles | B-factor error | no | advanced |

## Parameter `binwidth`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask | Envelope binary layers width | no | advanced |

## Parameter `box`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| crys_score | Particle box size | yes | standard |
| extract | Particle box size | no | standard |
| noisevol | Particle box size | yes | standard |
| print_dose_weights | Particle box size | yes | standard |
| print_fsc | Particle box size | no | advanced |
| print_magic_boxes | Particle box size | yes | standard |
| reextract | Particle box size | no | advanced |
| simulate_nanoparticle | Particle box size | yes | standard |
| simulate_noise | Particle box size | yes | standard |

## Parameter `box_coarse`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Box size for coarse sieving | no | advanced |

## Parameter `box_extract`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Force box size (px, optional) | no | developer |
| pick_extract | Extracted particle image size | no | advanced |
| tseries_extractor | Extracted particle image size | yes | standard |

## Parameter `box_fine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Box size for fine sieving | no | advanced |

## Parameter `boxes`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| print_project_field | output coordinates in JSON format | no | developer |

## Parameter `boxfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| track_particles | List of particle coordinates | yes | standard |
| tseries_motion_correct | List of particle coordinates | no | advanced |

## Parameter `boxtab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_boxes | List of box files | yes | standard |
| import_movies | List of box files | no | advanced |

## Parameter `cavg_ini`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | 3D initialization on class averages | no | advanced |

## Parameter `cavg_ini_ext`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | External class-average 3D initialization | no | advanced |

## Parameter `cenlp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Centering low-pass limit | no | advanced |
| abinitio2D_chunks | Centering low-pass limit | no | advanced |
| abinitio3D | Centering low-pass limit | no | advanced |
| abinitio3D_cavgs | Centering low-pass limit | no | advanced |
| abinitio3D_nano | Centering low-pass limit | no | advanced |
| autorefine3D_nano | Centering low-pass limit | no | advanced |
| center | Centering low-pass limit | no | developer |
| cluster2D_nano | Centering low-pass limit | no | advanced |
| refine3D | Centering low-pass limit | no | advanced |
| refine3D_nano | Centering low-pass limit | no | advanced |
| symaxis_search | Centering low-pass limit | no | advanced |
| symmetrize_map | Centering low-pass limit | no | advanced |
| symmetry_test | Centering low-pass limit | no | advanced |
| track_particles | Centering low-pass limit | no | advanced |

## Parameter `center`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Center class averages | no | advanced |
| abinitio2D_chunks | Center class averages | no | advanced |
| abinitio3D | Center reference volume(s) | no | advanced |
| abinitio3D_cavgs | Center reference volume(s) | no | advanced |
| autorefine3D_nano | Center reference volume(s) | no | advanced |
| cluster2D_nano | Center class averages | no | advanced |
| mask | Center input volume | no | advanced |
| refine3D | Center reference volume(s) | no | advanced |
| refine3D_nano | Center reference volume(s) | no | advanced |
| symaxis_search | Center input volume | no | advanced |
| symmetrize_map | Center input volume | no | advanced |
| symmetry_test | Center input volume | no | advanced |

## Parameter `center_pdb`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pdb2mrc | Whether to move the PDB atomic center to the center of the box | no | advanced |

## Parameter `chunk_count_penalty`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Chunk-count penalty | no | developer |

## Parameter `chunk_max_len`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Maximum latent chunk length | no | developer |

## Parameter `chunk_max_shift`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Maximum boundary shift | no | developer |

## Parameter `chunk_min_len`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Minimum latent chunk length | no | developer |

## Parameter `chunk_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Trajectory chunking mode | no | developer |

## Parameter `ciffile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cif2mrc | PDBx/mmCIF input coordinates file | yes | standard |
| cif2pdb | PDBx/mmCIF input coordinates file | yes | standard |

## Parameter `class`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Optional class index to split | no | advanced |
| extract_subproj | 2D class index | no | advanced |
| stackops | Class index | no | advanced |

## Parameter `classtats`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oristats | Class statistics | no | developer |

## Parameter `clip`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| scale | Clipped box size | no | advanced |
| stack | Clipped box size | no | advanced |

## Parameter `cls_init`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Scheme for initial class generation | no | advanced |

## Parameter `clust_crit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cluster_cavgs | Clustering criterion | no | advanced |
| cluster_cavgs_selection | Clustering criterion | no | advanced |
| cluster_stack | Clustering criterion | no | advanced |
| match_cavgs | Clustering criterion | no | developer |
| match_stacks | Clustering criterion | no | developer |

## Parameter `clustind`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract_subproj | Cluster index | no | advanced |
| select_clusters | Cluster index | no | developer |

## Parameter `clustinds`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| select_clusters | Comma separated cluster indices | no | developer |

## Parameter `cn_stop`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| symmetry_test | Rotational symmetry order stop index | no | advanced |

## Parameter `combine_eo`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | Whether e/o references are combined for final alignment(yes\|no){no} | no | advanced |
| refine3D_auto | Whether e/o references are combined for final alignment(yes\|no){no} | no | advanced |

## Parameter `conical_fsc`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Conical FSC regularization | no | advanced |
| abinitio3D_cavgs | Conical FSC regularization | no | advanced |
| refine3D | Conical FSC regularization | no | advanced |

## Parameter `continue`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | Continue previous refinement | no | advanced |
| refine3D_auto | Continue previous refinement | no | advanced |
| refine3D_multi | Continue previous refinement | no | advanced |
| refine3D_nano | Continue previous refinement | no | advanced |

## Parameter `cs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Spherical aberration | yes | standard |
| import_movies | Spherical aberration | yes | standard |
| import_particles | Spherical aberration | yes | standard |
| master | Spherical aberration (mm) | yes | standard |
| mini_stream | Spherical aberration | yes | standard |
| preproc | Spherical aberration | yes | standard |
| simulate_movie | Spherical aberration | no | advanced |
| simulate_particles | Spherical aberration | no | advanced |
| tseries_import | Spherical aberration | yes | standard |

## Parameter `ctf`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctfops | CTF status | yes | standard |
| import_movies | CTF status | no | advanced |
| import_particles | CTF status | no | advanced |
| oriops | CTF status | no | developer |
| reimport_particles | CTF status | no | developer |
| simulate_particles | CTF status | yes | standard |

## Parameter `ctf_correct_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_correct | CTF correction mode | no | advanced |

## Parameter `ctfpatch`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_estimate | Patch CTF estimation | no | advanced |
| preprocess | Patch CTF estimation | no | advanced |

## Parameter `ctfresthreshold`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick_extract | CTF Resolution rejection threshold | no | advanced |
| selection | CTF Resolution rejection threshold | no | advanced |

## Parameter `ctfstats`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oristats | CTF statistics | no | developer |

## Parameter `defocus`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_movie | Underfocus | no | advanced |
| simulate_particles | Underfocus | no | advanced |

## Parameter `deftab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctfops | CTF parameter file | no | advanced |
| import_movies | Pre-determined per-micrograph CTF parameters | no | advanced |
| import_particles | CTF parameter file | no | advanced |
| import_trajectory | CTF parameter file | no | advanced |
| simulate_particles | CTF parameter file | no | advanced |

## Parameter `deselfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| selection | File with deselection indices | no | advanced |

## Parameter `dferr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Underfocus error half-width | no | developer |
| simulate_particles | Underfocus error half-width | no | advanced |

## Parameter `dfmax`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_estimate | Expected maximum defocus | no | advanced |
| preproc | Expected maximum defocus | no | developer |
| preprocess | Expected maximum defocus | no | advanced |

## Parameter `dfmin`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_estimate | Expected minimum defocus | no | advanced |
| preproc | Expected minimum defocus | no | developer |
| preprocess | Expected minimum defocus | no | advanced |
| selection | Expected minimum defocus | no | advanced |

## Parameter `dfunit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | Underfocus unit | no | advanced |

## Parameter `dir`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Project directory | no | advanced |
| pick | Output directory | no | advanced |

## Parameter `dir_box`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract | Box files directory | no | advanced |

## Parameter `dir_exec`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D_stream | Previous run directory | no | developer |
| pick_extract | Previous run directory | no | advanced |
| sieve_cavgs | Previous run directory | no | advanced |

## Parameter `dir_meta`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Input metadata directory | no | developer |
| preproc | Directory containing per-movie metadata in XML format | no | standard |

## Parameter `dir_movies`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_movies | Input movies directory | no | advanced |
| master | Input movies directory | yes | standard |
| preproc | Input movies directory | yes | standard |

## Parameter `dir_preprocess`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Pre-existing preprocessing directory | no | developer |

## Parameter `dir_prev`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| preproc | Previous run directory | no | developer |

## Parameter `dir_target`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D_stream | Target directory | yes | standard |
| assign_optics | Target directory | yes | standard |
| gen_pickrefs | Target directory | yes | standard |
| pick_extract | Target directory | yes | standard |
| sieve_cavgs | Target directory | yes | standard |

## Parameter `dm_alpha`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Diffusion-map density normalization (default 0) | no | advanced |
| denoise_project | Diffusion-map density normalization (default 0) | no | advanced |
| ppca_denoise | Diffusion-map density normalization (default 0) | no | advanced |

## Parameter `e1`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Rotation along Phi | no | developer |
| volops | Rotation along Phi | no | advanced |

## Parameter `e2`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Rotation along Theta | no | developer |
| volops | Rotation along Theta | no | advanced |

## Parameter `e3`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Rotation along Psi | no | developer |
| volops | Rotation along Psi | no | advanced |

## Parameter `edge`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask | Envelope mask soft edge | no | advanced |
| automask2D | Envelope mask soft edge | no | advanced |
| mask | Envelope mask soft edge | no | advanced |

## Parameter `eer_fraction`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | # of EER frames to fraction together | no | advanced |
| preproc | # of EER frames to fraction together | no | developer |
| preprocess | # of EER frames to fraction together | no | advanced |

## Parameter `element`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| analysis2D_nano | Atom element name: Au, Pt etc. | yes | standard |
| atoms_rmsd | Atom element name: Au, Pt etc. | yes | standard |
| atoms_stats | Atom element name: Au, Pt etc. | yes | standard |
| autorefine3D_nano | Atom element name: Au, Pt etc. | yes | standard |
| conv_atom_denoise | Atom element name: Au, Pt etc. | yes | standard |
| core_atoms_analysis | Atom element name: Au, Pt etc. | yes | standard |
| crys_score | Atom element name: Au, Pt etc. | yes | standard |
| detect_atoms | Atom element name: Au, Pt etc. | yes | standard |
| simulate_nanoparticle | Atom element name: Au, Pt etc. | no | advanced |

## Parameter `envfsc`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| reconstruct3D | Envelope mask e/o maps for FSC | no | advanced |
| refine3D | Envelope mask e/o maps for FSC | no | advanced |

## Parameter `even`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_oris | Generate even projections | no | developer |
| simulate_particles | Generate even projections | no | advanced |

## Parameter `fbody`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fsc_area_score | Output file body | no | advanced |
| motion_correct | Template output micrograph name | no | advanced |
| preprocess | Template output micrograph name | no | advanced |
| track_particles | Template output tracked series | yes | standard |

## Parameter `filetab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | List of files | yes | standard |
| import_movies | List of movie files | no | advanced |
| mini_stream | List of files | yes | standard |
| model_cavgs_rejection | Analysis file table | no | advanced |
| scale | Stacks list | no | advanced |
| stack | Stacks list | yes | standard |
| tsegmaps_core_finder | Volumes list | yes | standard |
| tseries_import | List of individual movie frame files | yes | standard |
| volanalyze | Volumes list | yes | standard |
| volcluster | Volumes list | yes | standard |

## Parameter `fill_holes`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| binarize | Fill holes | no | developer |

## Parameter `filt_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Filtering mode | no | advanced |
| refine3D_auto | Filtering mode | no | advanced |
| refine3D_multi | Filtering mode | no | advanced |

## Parameter `filter`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | Filter type(bs\|nlmean\|no){no} | no | advanced |
| track_particles | Alternative filter for particle tracking | no | advanced |

## Parameter `fit_phshift`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Fit CTF phase shift | no | developer |
| ctf_estimate | Fit CTF phase shift | no | advanced |
| master | Phase plate | no | developer |
| mini_stream | Fit CTF phase shift | no | advanced |
| preproc | Fit CTF phase shift | no | standard |
| preprocess | Fit CTF phase shift | no | advanced |

## Parameter `flipgain`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fractionate_movies | Flip the gain reference | no | advanced |
| master | Gain processing | no | developer |
| motion_correct | Flip the gain reference | no | advanced |
| preproc | Flip the gain reference | no | developer |
| preprocess | Flip the gain reference | no | advanced |

## Parameter `fname`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_register | PDB file list | yes | standard |
| crys_score | PDB file list | yes | standard |
| info_image | Name of image file | yes | standard |
| model_cavgs_rejection | Quality model output | no | advanced |
| write_mic_filetab | Filename micrograph list | yes | standard |

## Parameter `force_lp_range`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Force low-pass range | no | advanced |

## Parameter `frac_best`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| bootstrap_cavgs | Anchor fraction | no | developer |
| sample_classes | Fraction of best particles to sample from | no | advanced |

## Parameter `frac_diam`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_rmsd | Fraction of atomic diameter | no | advanced |
| core_atoms_analysis | Fraction of atomic diameter | no | developer |

## Parameter `frac_worst`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| sample_classes | Fraction of worst particles to sample from | no | advanced |

## Parameter `fraca`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Amplitude contrast fraction | no | developer |
| import_movies | Amplitude contrast fraction | yes | standard |
| import_particles | Amplitude contrast fraction | yes | standard |
| master | Amplitude contrast fraction | yes | standard |
| mini_stream | Amplitude contrast fraction | no | advanced |
| preproc | Amplitude contrast fraction | yes | standard |
| simulate_movie | Amplitude contrast fraction | no | advanced |
| simulate_particles | Amplitude contrast fraction | no | advanced |
| tseries_import | Amplitude contrast fraction | yes | standard |

## Parameter `fraction_dose_target`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | EER fraction dose target (e/Ang^2) | no | advanced |
| preproc | EER fraction dose target (e/Ang^2) | no | standard |
| preprocess | EER fraction dose target (e/Ang^2) | no | advanced |

## Parameter `frcs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | Projection FRCs file | no | advanced |
| print_fsc | Projection FRCs file | no | advanced |

## Parameter `fromf`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fractionate_movies | Starting frame | no | advanced |
| track_particles | Frame to start tracking from | no | advanced |
| tseries_make_pickavg | Frame to start averaging from | no | advanced |

## Parameter `fromp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract_subproj | From index | no | advanced |
| extract_substk | From index | no | advanced |
| scale | First stack index | no | advanced |
| stackops | From particle index | no | advanced |
| trajectory_reconstruct3D | From particle index | no | developer |

## Parameter `fsc`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_correct | FSC file | no | advanced |
| filter | FSC file | no | advanced |
| postprocess | FSC file | no | advanced |
| print_fsc | FSC file | no | advanced |

## Parameter `gainref`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Gain reference | no | developer |
| motion_correct | Gain reference | no | standard |
| preproc | Gain reference | no | standard |
| preprocess | Gain reference | no | advanced |

## Parameter `graph`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Class split graph | no | advanced |
| denoise_project | Diffusion graph | no | advanced |

## Parameter `greedy_sampling`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| sample_classes | Greedy balanced selection | no | advanced |
| selection | Greedy balanced selection | no | advanced |

## Parameter `guinier`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| volops | Guinier plot | no | advanced |

## Parameter `hp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | High-pass limit | no | advanced |
| abinitio2D_chunks | High-pass limit | no | advanced |
| abinitio3D | High-pass limit | no | advanced |
| abinitio3D_cavgs | High-pass limit | no | advanced |
| abinitio3D_nano | High-pass limit | no | advanced |
| autorefine3D_nano | High-pass limit | no | advanced |
| cluster2D_nano | High-pass limit | no | advanced |
| cluster_cavgs | High-pass limit | no | advanced |
| cluster_cavgs_selection | High-pass limit | no | advanced |
| cluster_stack | High-pass limit | no | advanced |
| ctf_estimate | High-pass limit | no | advanced |
| dock_volpair | High-pass limit | no | advanced |
| filter | High-pass limit | no | advanced |
| fsc | High-pass limit | no | advanced |
| match_cavgs | High-pass limit | no | developer |
| match_stacks | High-pass limit | yes | standard |
| ppca_denoise | High-pass limit | no | advanced |
| refine3D | High-pass limit | no | advanced |
| refine3D_nano | High-pass limit | no | advanced |
| symaxis_search | High-pass limit | no | advanced |
| symmetrize_map | High-pass limit | no | advanced |
| symmetry_test | High-pass limit | no | advanced |
| track_particles | High-pass limit | no | advanced |
| trajectory_denoise | High-pass limit | no | advanced |
| volanalyze | High-pass limit | yes | standard |
| volcluster | High-pass limit | no | advanced |
| volops | High-pass limit | no | advanced |

## Parameter `hp_ctf_estimate`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| preprocess | High-pass limit for CTF parameter estimation | no | advanced |

## Parameter `icefracthreshold`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick_extract | Ice Fraction rejection threshold | no | advanced |
| selection | Ice Fraction rejection threshold | no | advanced |

## Parameter `icm`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Whether to perform ICM filtering of reference(s) | no | advanced |
| flex_analysis | Automatic diffusion-mode selection | no | advanced |
| refine3D_nano | Whether to perform ICM filtering of reference(s) | no | advanced |
| uniform_filter3D | Whether to perform ICM filtering of reference(s) | no | advanced |

## Parameter `imgkind`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| postprocess | Volume kind | no | advanced |

## Parameter `import_dir`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_starproject | Import directory | yes | standard |

## Parameter `import_type`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_starproject | Import type | no | standard |

## Parameter `infile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| model_cavgs_rejection | Quality model input | no | advanced |
| selection | File with selection state (0/1) flags | no | advanced |
| tseries_extractor | Selected particle trajectory | yes | standard |

## Parameter `job_memory_per_task`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Memory per computing node | no | advanced |
| update_project | Memory per computing node | no | advanced |

## Parameter `json`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| print_project_field | output in JSON format | no | developer |

## Parameter `k_nn`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Diffusion graph neighbors (default 10; try 5-30) | no | advanced |
| denoise_project | Diffusion graph neighbors (default 10; try 5-30) | no | advanced |
| flex_analysis | Nearest neighbors (default 100) | no | advanced |
| ppca_denoise | Diffusion graph neighbors (default 5; try 5-30) | no | advanced |
| trajectory_denoise | Diffusion graph neighbors (default 5; try 5-30) | no | advanced |

## Parameter `kpca_backend`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | Kernel PCA backend | no | advanced |
| ppca_denoise_classes | Kernel PCA backend | no | advanced |
| trajectory_denoise | Kernel PCA backend | no | advanced |

## Parameter `kpca_ker`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | Kernel PCA kernel | no | advanced |
| ppca_denoise_classes | Kernel PCA kernel | no | advanced |
| trajectory_denoise | Kernel PCA kernel | no | advanced |

## Parameter `kpca_nystrom_local_nbrs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | Nyström max local support neighbors (default 96; try 96, 128) | no | advanced |
| ppca_denoise_classes | Nyström max local support neighbors (default 96; try 96, 128) | no | advanced |
| trajectory_denoise | Nyström max local support neighbors (default 96; try 96, 128) | no | advanced |

## Parameter `kpca_nystrom_npts`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | no | advanced |
| ppca_denoise_classes | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | no | advanced |
| trajectory_denoise | Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512) | no | advanced |

## Parameter `kpca_rbf_gamma`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | RBF gamma (0 => auto) | no | advanced |
| ppca_denoise_classes | RBF gamma (0 => auto) | no | advanced |
| trajectory_denoise | RBF gamma (0 => auto) | no | advanced |

## Parameter `kv`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Acceleration voltage | yes | standard |
| import_movies | Acceleration voltage | yes | standard |
| import_particles | Acceleration voltage | yes | standard |
| master | Acceleration voltage (kV) | yes | standard |
| mini_stream | Acceleration voltage | yes | standard |
| preproc | Acceleration voltage | yes | standard |
| print_dose_weights | Acceleration voltage | no | advanced |
| simulate_movie | Acceleration voltage | no | advanced |
| simulate_particles | Acceleration voltage | no | advanced |
| tseries_import | Acceleration voltage | yes | standard |

## Parameter `lambda`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | BS smoother lambda | no | advanced |
| icm2D | ICM lambda regularization parameter | no | advanced |
| icm3D | ICM lambda regularization parameter | no | advanced |
| uniform_filter3D | ICM lambda regularization parameter | no | advanced |

## Parameter `lp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Low-pass limit | no | advanced |
| abinitio2D_chunks | Low-pass limit | no | advanced |
| abinitio3D | Low-pass limit | no | advanced |
| autorefine3D_nano | Initial low-pass limit | yes | standard |
| cluster2D_nano | Static low-pass limit | no | advanced |
| cluster_cavgs | Low-pass limit | no | advanced |
| cluster_cavgs_selection | Low-pass limit | no | advanced |
| cluster_stack | Low-pass limit | yes | standard |
| ctf_estimate | Low-pass limit | no | advanced |
| dock_volpair | Low-pass limit | no | advanced |
| estimate_diam | low-pass limit in Angstroms{7.} | no | advanced |
| filter | Low-pass limit | no | advanced |
| flex_analysis | Low-pass limit | no | standard |
| fsc | Low-pass limit | no | advanced |
| match_cavgs | Low-pass limit | no | developer |
| match_stacks | Low-pass limit | yes | standard |
| pick | Low-pass limit | no | advanced |
| postprocess | Low-pass limit for map filtering | no | advanced |
| ppca_denoise | Low-pass limit | no | advanced |
| refine3D | Static low-pass limit | no | advanced |
| refine3D_nano | Static low-pass limit | no | advanced |
| symaxis_search | Low-pass limit | no | advanced |
| symmetrize_map | Low-pass limit | no | advanced |
| symmetry_test | Low-pass limit | no | advanced |
| trajectory_denoise | Low-pass limit | no | advanced |
| trajectory_reconstruct3D | Low-pass limit | no | developer |
| volanalyze | Low-pass limit | yes | standard |
| volcluster | Low-pass limit | yes | standard |
| volops | Low-pass limit | no | advanced |

## Parameter `lp_backgr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| mask | Background low-pass resolution | no | advanced |
| refine3D | Background low-pass resolution | no | advanced |

## Parameter `lp_ctf_estimate`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| preprocess | Low-pass limit for CTF parameter estimation | no | advanced |

## Parameter `lp_pick`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick_extract | Low-pass limit for picking | no | advanced |
| track_particles | Low-pass limit in Angs | no | advanced |

## Parameter `lplim_crit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fsc_area_score | Low-pass limit FSC criterion | no | advanced |
| refine3D | Low-pass limit FSC criterion | no | advanced |

## Parameter `lpstart`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Initial low-pass limit | no | advanced |
| abinitio3D | Starting low-pass limit | no | advanced |
| abinitio3D_cavgs | Starting low-pass limit | no | advanced |
| abinitio3D_nano | Starting low-pass limit | no | advanced |
| cluster2D_nano | Initial low-pass limit | no | advanced |
| motion_correct | Initial low-pass limit | no | advanced |
| particle_sieving | Initial low-pass limit | no | advanced |
| preprocess | Initial low-pass limit for movie alignment | no | advanced |
| tseries_make_pickavg | Initial low-pass limit | no | advanced |
| tseries_motion_correct | Initial low-pass limit | no | advanced |
| uniform_filter2D | Starting resolution limit | yes | standard |
| uniform_filter3D | Starting resolution limit | yes | standard |

## Parameter `lpstart_ini3D`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Starting low-pass limit ini3D | no | advanced |

## Parameter `lpstop`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Final low-pass limit | no | advanced |
| abinitio2D_chunks | Final low-pass limit | no | advanced |
| abinitio3D | Final low-pass limit | no | advanced |
| abinitio3D_cavgs | Final low-pass limit | no | advanced |
| abinitio3D_nano | Final low-pass limit | no | advanced |
| cluster2D_nano | Final low-pass limit | no | advanced |
| motion_correct | Final low-pass limit | no | advanced |
| preprocess | Final low-pass limit for movie alignment | no | advanced |
| refine3D | Low-pass limit for frequency limited refinement | no | advanced |
| refine3D_multi | Low-pass limit for frequency limited refinement | no | advanced |
| tseries_make_pickavg | Final low-pass limit | no | advanced |
| tseries_motion_correct | Final low-pass limit | no | advanced |
| uniform_filter2D | Stopping resolution limit | yes | standard |
| uniform_filter3D | Stopping resolution limit | yes | standard |

## Parameter `lpstop_coarse`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Final low-pass limit for coarse sieving | no | advanced |

## Parameter `lpstop_fine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Final low-pass limit for fine sieving | no | advanced |

## Parameter `lpstop_ini3D`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Final low-pass limit ini3D | no | advanced |

## Parameter `makemovie`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| stackops | Whether to make a movie | no | advanced |

## Parameter `max_dose`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | Maximum dose threshold(e/A2) | no | advanced |
| preproc | Maximum dose threshold(e/A2) | no | developer |
| preprocess | Maximum dose threshold(e/A2) | no | advanced |

## Parameter `maxits`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Max iterations | no | advanced |
| cluster2D_nano | Max iterations | no | advanced |
| reconstruct3D_pcg | Max iterations | no | developer |
| refine3D | Max iterations | no | advanced |
| refine3D_auto | Max iterations | no | advanced |
| refine3D_multi | Max iterations | no | advanced |
| refine3D_nano | Max iterations | no | advanced |
| trajectory_reconstruct3D | Max iterations | no | developer |

## Parameter `maxnchunks`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Max # of chunks to process | no | advanced |

## Parameter `maxpop`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| assign_optics_groups | Maximum number of movies/micrographs/stacks in each optics group | no | developer |

## Parameter `mcconvention`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fractionate_movies | Movie alignment convention | no | advanced |
| motion_correct | Frame of reference during movie alignment | no | advanced |
| preprocess | Frame of reference during movie alignment | no | advanced |

## Parameter `mcpatch`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | Patch-based motion correction | no | advanced |
| preprocess | Patch-based motion correction | no | advanced |
| tseries_make_pickavg | Patch-based motion correction | no | advanced |
| tseries_motion_correct | Patch-based motion correction | no | advanced |

## Parameter `mcpatch_thres`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | Use motion correction patch threshold | no | advanced |
| preprocess | Use motion correction patch threshold | no | advanced |

## Parameter `mirr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Mirror orientations | no | developer |
| postprocess | Perform mirroring | no | advanced |
| stackops | Perform mirroring | no | advanced |
| volops | Perform mirroring | no | advanced |

## Parameter `ml_reg`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | ML regularization | no | advanced |
| refine3D_multi | ML regularization | no | advanced |

## Parameter `moldiam`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cluster2D_nano | Molecular diameter | no | advanced |
| pick | Molecular diameter | no | advanced |
| pick_extract | Molecular diameter | no | advanced |
| print_magic_boxes | Molecular diameter | no | developer |
| simulate_nanoparticle | Molecular diameter | no | advanced |

## Parameter `moldiam_max`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| mini_stream | Max molecular diameter | no | advanced |
| pick | Max molecular diameter | no | advanced |
| pick_extract | Max molecular diameter | no | advanced |

## Parameter `mskdiam`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Mask diameter | yes | standard |
| abinitio2D_chunks | Mask diameter | yes | standard |
| abinitio2D_stream | Mask diameter | no | standard |
| abinitio3D | Mask diameter | yes | standard |
| abinitio3D_cavgs | Mask diameter | yes | standard |
| abinitio3D_nano | Mask diameter | yes | standard |
| automask2D | Mask diameter | yes | standard |
| autorefine3D_nano | Mask diameter | yes | standard |
| bootstrap_rec3D | Mask diameter | yes | standard |
| cavgseoproc_nano | Mask diameter | yes | standard |
| cavgsproc_nano | Mask diameter | yes | standard |
| cls_split | Mask diameter | no | standard |
| cluster2D_nano | Mask diameter | yes | standard |
| cluster_cavgs | Mask diameter | yes | standard |
| cluster_cavgs_selection | Mask diameter | yes | standard |
| cluster_stack | Mask diameter | yes | standard |
| conv_atom_denoise | Mask diameter | yes | standard |
| denoise_project | Mask diameter | no | standard |
| detect_atoms | Mask diameter | no | advanced |
| dock_volpair | Mask diameter | yes | standard |
| estimate_diam | Mask diameter | yes | standard |
| flex_analysis | Mask diameter | no | standard |
| fsc | Mask diameter | yes | standard |
| fsc_area_score | Mask diameter | yes | standard |
| mask | Mask diameter | no | advanced |
| match_cavgs | Mask diameter | yes | standard |
| match_stacks | Mask diameter | yes | standard |
| model_cavgs_rejection | Mask diameter | yes | standard |
| normalize | Mask diameter | yes | standard |
| nu_filt3D | Mask diameter | yes | standard |
| postprocess | Mask diameter | yes | standard |
| ppca_denoise | Mask diameter | no | advanced |
| ppca_volvar | Mask diameter | no | developer |
| ptclsproc_nano | Mask diameter | yes | standard |
| reconstruct3D | Mask diameter | yes | standard |
| reconstruct3D_pcg | Mask diameter | no | developer |
| refine3D | Mask diameter | yes | standard |
| refine3D_auto | Mask diameter | yes | standard |
| refine3D_multi | Mask diameter | yes | standard |
| refine3D_nano | Mask diameter | yes | standard |
| reproject | Mask diameter | no | advanced |
| sieve_cavgs | Mask diameter | no | standard |
| simulate_particles | Mask diameter | yes | standard |
| symaxis_search | Mask diameter | yes | standard |
| symmetrize_map | Mask diameter | yes | standard |
| symmetry_test | Mask diameter | yes | standard |
| trajectory_denoise | Mask diameter | no | advanced |
| trajectory_make_projavgs | Mask diameter | yes | standard |
| trajectory_reconstruct3D | Mask diameter | yes | standard |
| uniform_filter2D | Mask diameter | yes | standard |
| uniform_filter3D | Mask diameter | yes | standard |
| validate_cavgs_vs_model | Mask diameter | yes | standard |
| volanalyze | Mask diameter | yes | standard |
| volcluster | Mask diameter | yes | standard |
| volops | Mask diameter | no | advanced |

## Parameter `mskdiam_detect`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Detect-atoms mask diameter | no | advanced |

## Parameter `mul`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_cavgs | Shift multiplication factor | no | advanced |
| oriops | Shift multiplication factor | no | developer |
| volops | Multiplication factor | no | advanced |

## Parameter `multi_moldiams`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick | Comma-separated molecular diameters with which to execute multiple gaussian pick | no | advanced |

## Parameter `multivol_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Multi-volume ab initio mode | no | advanced |
| abinitio3D_cavgs | Multi-volume class-average ab initio mode | no | advanced |
| refine3D_multi | Multi-volume refinement mode | no | advanced |

## Parameter `nang_nbrs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| flex_analysis | Angular candidate cap (default 1000) | no | advanced |

## Parameter `nboxes_max`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Max # boxes per micrograph | no | developer |
| pick | Max # boxes per micrograph | no | advanced |

## Parameter `nchunks`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D_chunks | Number of chunks | no | advanced |
| sieve_cavgs | Number of subsets to classify simultaneously | yes | standard |
| trajectory_reconstruct3D | Number of temporal chunks | no | developer |

## Parameter `nchunks_max`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Maximum automatic chunk count | no | developer |

## Parameter `nchunks_min`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Minimum automatic chunk count | no | developer |

## Parameter `nchunksperset`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| sieve_cavgs | Number of subsets to group | no | standard |

## Parameter `ncls`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Number of 2D clusters | yes | standard |
| abinitio2D_stream | Maximum number of 2D clusters | yes | standard |
| center2D_nano | Number of 2D clusters | yes | standard |
| cls_split | Fixed number of subclasses (0 => auto) | no | advanced |
| cluster_cavgs | Number of clusters | no | advanced |
| cluster_stack | Number of clusters | no | advanced |
| make_cavgs | Number of 2D clusters | no | advanced |
| make_oris | Number of random class labels | no | developer |
| pick | Cluster input pickrefs before template generation | no | advanced |
| sieve_cavgs | Number of 2D clusters | yes | standard |
| volcluster | Number of 2D clusters | no | advanced |

## Parameter `ncls_coarse`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Number of clusters for coarse sieving | no | advanced |

## Parameter `ncls_fine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Number of clusters for fine sieving | no | advanced |

## Parameter `ncunits`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| track_particles | Concurrent particle trackers | yes | standard |

## Parameter `ndev`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| binarize | Binarization threshold | no | developer |
| pick | # of sigmas for outlier detection | no | standard |

## Parameter `ndiscrete`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_oris | Number of discrete projection directions | no | developer |
| oriops | Number of discrete projection directions | no | developer |

## Parameter `neg`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctfops | Invert contrast | no | advanced |
| reproject | Invert contrast | no | advanced |
| stackops | Invert contrast | no | advanced |
| track_particles | Invert contrast | no | advanced |
| volops | Invert contrast | no | advanced |

## Parameter `neigs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Number of eigencomponents (0 => auto scan; default 200) | no | advanced |
| denoise_project | Number of eigencomponents (0 => auto scan; default 200) | no | advanced |
| flex_analysis | Maximum number of diffusion modes (default 15) | no | advanced |
| ppca_denoise | Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160) | yes | standard |
| ppca_denoise_classes | # eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160) | no | advanced |
| ppca_volvar | # eigenvecs | yes | standard |
| trajectory_denoise | Number of diffusion-map components (0 => auto; default 0) | no | advanced |
| trajectory_reconstruct3D | Flex latent dimensions | no | developer |

## Parameter `newbox`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| scale | Scaled box size | no | advanced |

## Parameter `nframes`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| print_dose_weights | Number of frames | yes | standard |
| simulate_movie | Number of frames | no | advanced |
| simulate_particles | # of particle frames | no | advanced |

## Parameter `nframesgrp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| stackops | Number of stack entries to group & average | no | advanced |
| track_particles | Number of contigous frames to average | no | advanced |
| tseries_make_pickavg | # contigous frames to average | no | advanced |
| tseries_motion_correct | # frames in time moving time window | no | advanced |

## Parameter `ngrow`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask2D | # layers to grow | no | advanced |

## Parameter `nicedispid`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Optics group offset delta multiplier | no | developer |

## Parameter `ninipick`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| preproc | Number of micrographs to perform initial picking preprocessing on | no | developer |

## Parameter `nmics`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| gen_pickrefs | Number of micrographs to import | no | advanced |
| particle_sieving | Max # of micrographs per chunk | no | advanced |

## Parameter `nmoldiams`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick | Number of molecular diameters to investigate | no | advanced |
| pick_extract | Number of molecular diameters to investigate | no | advanced |

## Parameter `noise_norm`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| normalize | Noise normalize | no | advanced |

## Parameter `norm`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| normalize | Normalize | no | advanced |

## Parameter `nparts`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Number of computing nodes | no | standard |
| abinitio2D_chunks | Number of computing nodes | no | standard |
| abinitio2D_stream | Number of computing nodes | yes | standard |
| abinitio3D | Number of computing nodes | no | standard |
| abinitio3D_cavgs | Number of computing nodes | no | standard |
| abinitio3D_nano | Number of computing nodes | no | advanced |
| bootstrap_cavgs | Number of computing nodes | no | developer |
| bootstrap_rec3D | Number of computing nodes | no | developer |
| cls_split | Number of computing nodes | no | standard |
| cluster2D_nano | Number of computing nodes | no | advanced |
| ctf_estimate | Number of computing nodes | yes | standard |
| denoise_project | Number of computing nodes | no | standard |
| extract | Number of computing nodes | yes | standard |
| flex_analysis | Number of computing nodes | no | standard |
| fractionate_movies | Number of computing nodes | yes | standard |
| gen_pspecs_and_thumbs | Number of computing nodes | yes | standard |
| make_cavgs | Number of computing nodes | yes | standard |
| motion_correct | Number of computing nodes | yes | standard |
| particle_sieving | Number of chunks classified simultaneously | no | standard |
| pick | Number of computing nodes | yes | standard |
| pick_extract | Number of computing nodes | yes | standard |
| preproc | Number of computing nodes | yes | standard |
| preprocess | Number of computing nodes | yes | standard |
| prune_project | Number of computing nodes | yes | standard |
| reconstruct3D | Number of computing nodes | no | advanced |
| reextract | Number of computing nodes | yes | standard |
| refine3D | Number of computing nodes | no | standard |
| refine3D_auto | Number of computing nodes | yes | standard |
| refine3D_multi | Number of computing nodes | yes | standard |
| refine3D_nano | Number of computing nodes | no | advanced |
| sample_classes | Number of partitions in balancing | no | advanced |
| selection | Number of partitions in balancing | no | advanced |
| sieve_cavgs | Number of computing nodes | yes | standard |
| split | Number of computing nodes | yes | standard |
| split_stack | Number of parts balanced splitting of the stack | yes | standard |
| trajectory_reconstruct3D | Number of computing nodes | no | developer |
| tseries_motion_correct | Number of computing nodes | yes | standard |

## Parameter `npreimages`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| flex_analysis | Representative state volumes (default 8) | no | advanced |

## Parameter `nptcls`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_oris | Number of per-particle orientations | yes | standard |
| simulate_noise | Number of particles | yes | standard |
| simulate_particles | Number of particles | yes | standard |

## Parameter `nptcls_coarse`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Target # of particles per coarse chunk | no | advanced |

## Parameter `nptcls_fine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Target # of particles per fine chunk | no | advanced |

## Parameter `nptcls_per_cls`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D_chunks | Number of particles per cluster | no | standard |
| analysis2D_nano | Number of particles per cluster | no | advanced |
| check_refpick | Number of particles per class | no | developer |
| cluster2D_nano | Number of particles per cluster | no | advanced |
| mini_stream | Number of particles per class | no | advanced |
| sieve_cavgs | Number of particles per cluster | yes | standard |

## Parameter `nptcls_per_part`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| sample_classes | Number of ptcls per part to select when balancing | no | advanced |
| selection | Number of ptcls per part to select when balancing | no | advanced |

## Parameter `nran`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| selection | Number of random samples | no | advanced |
| stackops | Number of random samples | no | advanced |

## Parameter `nrestarts`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Number of restarts | no | advanced |

## Parameter `nsample`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Number of particles to sample | no | advanced |
| abinitio3D | Number of particles to sample | no | standard |
| abinitio3D_nano | Number of particles to sample | no | advanced |
| sample_classes | Number of particles to sample | no | advanced |

## Parameter `nsample_coarse`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Number of particles to sample in coarse sieving | no | advanced |

## Parameter `nsample_fine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Number of particles to sample in fine sieving | no | advanced |

## Parameter `nspace`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Number of projection directions | no | advanced |
| denoise_project | Number of projection directions | no | advanced |
| flex_analysis | Number of projection directions | yes | standard |
| fsc_area_score | Number of cone axes | no | advanced |
| make_cavgs | Number of projection directions | no | advanced |
| oristats | Number of projection directions | no | developer |
| refine3D | Number of projection directions | no | advanced |
| refine3D_nano | Number of projection directions | no | advanced |
| reproject | Number of projection directions | yes | standard |
| trajectory_make_projavgs | Number of projection directions | no | developer |
| vizoris | Number of projection directions | no | advanced |

## Parameter `nspace_sub`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| denoise_project | SO3 mixture subspace size | no | advanced |

## Parameter `nstages`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Last ab initio stage to run | no | advanced |
| estimate_lpstages | Number of low-pass limit stages | yes | standard |

## Parameter `nstates`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Number of states | no | standard |
| abinitio3D_cavgs | Number of states | no | standard |
| bootstrap_rec3D | Number of states | no | developer |
| make_oris | Number of random state labels | no | developer |
| noisevol | Number states | no | advanced |
| oriops | Number of random state labels | no | developer |
| ptcl3D_state_consensus | Number of states | no | advanced |
| refine3D | Number of states | no | advanced |
| refine3D_multi | Number of states | no | advanced |

## Parameter `nsubcls_max`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Maximum subclass trial count in auto mode (default 10) | no | advanced |

## Parameter `nsubcls_min`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Minimum subclass trial count in auto mode (default 3) | no | advanced |

## Parameter `nthr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Number of threads per computing node, give 0 if unsure | yes | standard |
| abinitio2D_chunks | Number of threads per computing node, give 0 if unsure | yes | standard |
| abinitio2D_stream | Number of threads per computing node, give 0 if unsure | yes | standard |
| abinitio3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| abinitio3D_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| abinitio3D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| analysis2D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| assign_optics | Number of threads per computing node, give 0 if unsure | yes | standard |
| atoms_register | Number of threads per computing node, give 0 if unsure | yes | standard |
| atoms_stats | Number of threads per computing node, give 0 if unsure | yes | standard |
| auto_spher_mask | Number of threads per computing node, give 0 if unsure | yes | standard |
| automask | Number of threads per computing node, give 0 if unsure | yes | standard |
| automask2D | Number of threads per computing node, give 0 if unsure | yes | standard |
| autorefine3D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| binarize | Number of threads per computing node, give 0 if unsure | yes | standard |
| bootstrap_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| bootstrap_rec3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| cavgseoproc_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| cavgsproc_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| center | Number of threads per computing node, give 0 if unsure | yes | standard |
| center2D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| check_refpick | Number of threads per computing node, give 0 if unsure | yes | standard |
| cls_split | Number of threads per computing node, give 0 if unsure | yes | standard |
| cluster2D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| cluster_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| cluster_cavgs_selection | Number of threads per computing node, give 0 if unsure | yes | standard |
| cluster_stack | Number of threads per computing node, give 0 if unsure | yes | standard |
| conv_atom_denoise | Number of threads per computing node, give 0 if unsure | yes | standard |
| crys_score | Number of threads per computing node, give 0 if unsure | yes | standard |
| ctf_correct | Number of threads per computing node, give 0 if unsure | yes | standard |
| ctf_estimate | Number of threads per computing node, give 0 if unsure | yes | standard |
| ctfops | Number of threads per computing node, give 0 if unsure | yes | standard |
| denoise_project | Number of threads per computing node, give 0 if unsure | yes | standard |
| detect_atoms | Number of threads per computing node, give 0 if unsure | yes | standard |
| dock_volpair | Number of threads per computing node, give 0 if unsure | yes | standard |
| estimate_diam | Number of threads per computing node, give 0 if unsure | yes | standard |
| extract | Number of threads per computing node, give 0 if unsure | yes | standard |
| filter | Number of threads per computing node, give 0 if unsure | yes | standard |
| flex_analysis | Number of threads per computing node, give 0 if unsure | yes | standard |
| fractionate_movies | Number of threads per computing node, give 0 if unsure | yes | standard |
| fsc | Number of threads per computing node, give 0 if unsure | yes | standard |
| fsc_area_score | Number of threads per computing node, give 0 if unsure | yes | standard |
| gen_pickrefs | Number of threads per computing node, give 0 if unsure | yes | standard |
| gen_pspecs_and_thumbs | Number of threads per computing node, give 0 if unsure | yes | standard |
| graphene_subtr | Number of threads per computing node, give 0 if unsure | yes | standard |
| icm2D | Number of threads per computing node, give 0 if unsure | yes | standard |
| icm3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| make_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| make_oris | Number of threads per computing node, give 0 if unsure | yes | standard |
| mask | Number of threads per computing node, give 0 if unsure | yes | standard |
| match_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| match_stacks | Number of threads per computing node, give 0 if unsure | yes | standard |
| mini_stream | Number of threads per computing node, give 0 if unsure | yes | standard |
| model_cavgs_rejection | Number of threads per computing node, give 0 if unsure | yes | standard |
| motion_correct | Number of threads per computing node, give 0 if unsure | yes | standard |
| normalize | Number of threads per computing node, give 0 if unsure | yes | standard |
| nu_filt3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| oristats | Number of threads per computing node, give 0 if unsure | yes | standard |
| particle_sieving | Number of threads per computing node, give 0 if unsure | yes | standard |
| pick | Number of threads per computing node, give 0 if unsure | yes | standard |
| pick_extract | Number of threads per computing node, give 0 if unsure | yes | standard |
| postprocess | Number of threads per computing node, give 0 if unsure | yes | standard |
| ppca_denoise | Number of threads per computing node, give 0 if unsure | yes | standard |
| ppca_denoise_classes | Number of threads per computing node, give 0 if unsure | yes | standard |
| ppca_volvar | Number of threads per computing node, give 0 if unsure | yes | standard |
| preproc | Number of threads per computing node, give 0 if unsure | yes | standard |
| preprocess | Number of threads per computing node, give 0 if unsure | yes | standard |
| ptclsproc_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| reconstruct3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| reconstruct3D_pcg | Number of threads per computing node, give 0 if unsure | yes | standard |
| reextract | Number of threads per computing node, give 0 if unsure | yes | standard |
| refine3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| refine3D_auto | Number of threads per computing node, give 0 if unsure | yes | standard |
| refine3D_multi | Number of threads per computing node, give 0 if unsure | yes | standard |
| refine3D_nano | Number of threads per computing node, give 0 if unsure | yes | standard |
| reproject | Number of threads per computing node, give 0 if unsure | yes | standard |
| sample_classes | Number of threads per computing node, give 0 if unsure | yes | standard |
| scale | Number of threads per computing node, give 0 if unsure | yes | standard |
| sieve_cavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| simulate_movie | Number of threads per computing node, give 0 if unsure | yes | standard |
| simulate_nanoparticle | Number of threads per computing node, give 0 if unsure | yes | standard |
| simulate_particles | Number of threads per computing node, give 0 if unsure | yes | standard |
| stackops | Number of threads per computing node, give 0 if unsure | yes | standard |
| symaxis_search | Number of threads per computing node, give 0 if unsure | yes | standard |
| symmetrize_map | Number of threads per computing node, give 0 if unsure | yes | standard |
| symmetry_test | Number of threads per computing node, give 0 if unsure | yes | standard |
| track_particles | Number of threads per computing node, give 0 if unsure | yes | standard |
| trajectory_denoise | Number of threads per computing node, give 0 if unsure | yes | standard |
| trajectory_make_projavgs | Number of threads per computing node, give 0 if unsure | yes | standard |
| trajectory_reconstruct3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| tseries_extractor | Number of threads per computing node, give 0 if unsure | yes | standard |
| tseries_make_pickavg | Number of threads per computing node, give 0 if unsure | yes | standard |
| tseries_motion_correct | Number of threads per computing node, give 0 if unsure | yes | standard |
| tseries_prep4tracking | Number of threads per computing node, give 0 if unsure | yes | standard |
| uniform_filter2D | Number of threads per computing node, give 0 if unsure | yes | standard |
| uniform_filter3D | Number of threads per computing node, give 0 if unsure | yes | standard |
| validate_cavgs_vs_model | Number of threads per computing node, give 0 if unsure | yes | standard |
| volanalyze | Number of threads per computing node, give 0 if unsure | yes | standard |
| volcluster | Number of threads per computing node, give 0 if unsure | no | advanced |
| volops | Number of threads per computing node, give 0 if unsure | yes | standard |

## Parameter `nthr_ini3D`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Number of threads for ini3D phase, give 0 if unsure | no | standard |

## Parameter `nu_refine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | NU resolution expansion refinement | no | advanced |
| refine3D_auto | NU resolution expansion refinement | no | advanced |

## Parameter `numlen`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| preprocess | Length of number string | no | advanced |

## Parameter `nxpatch`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | # of patches along x-axis | no | advanced |
| preprocess | # of patches along x-axis | no | advanced |
| tseries_make_pickavg | # of patches along x-axis | no | advanced |
| tseries_motion_correct | # of patches along x-axis | no | advanced |

## Parameter `nypatch`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | # of patches along y-axis | no | advanced |
| preprocess | # of patches along y-axis | no | advanced |
| tseries_make_pickavg | # of patches along y-axis | no | advanced |
| tseries_motion_correct | # of patches along y-axis | no | advanced |

## Parameter `objfun`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| reconstruct3D_pcg | Objective function | no | developer |
| refine3D | Objective function | no | advanced |

## Parameter `objfun_den`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Denoised objective | no | advanced |
| refine3D | Denoised objective | no | advanced |

## Parameter `objfun_den_w`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Denoised objective weight | no | advanced |
| refine3D | Denoised objective weight | no | advanced |

## Parameter `offset`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| track_particles | Shift half-width search bound | no | advanced |

## Parameter `optics_dir`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| gen_pickrefs | Target directory for optics import | no | advanced |

## Parameter `optics_offset`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| assign_optics_groups | Numbering offset to apply to optics groups | no | developer |
| export_relion | Offset to apply to optics group numbering | no | developer |

## Parameter `oritab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| center | Orientation and CTF parameter file | no | developer |
| ctfops | Orientation and CTF parameter file | no | advanced |
| import_particles | Orientation and CTF parameter file | no | advanced |
| mask | Orientation and CTF parameter file | no | advanced |
| oriops | Orientation and CTF parameter file | yes | standard |
| oristats | Orientation and CTF parameter file | yes | standard |
| reproject | Orientation and CTF parameter file | no | advanced |
| simulate_particles | Orientation and CTF parameter file | no | advanced |
| stackops | Orientation and CTF parameter file | no | advanced |
| vizoris | Orientation and CTF parameter file | yes | standard |

## Parameter `oritab2`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oristats | 2nd orientation and CTF parameter file | no | developer |

## Parameter `oritype`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Particle type to split | no | advanced |
| flex_analysis | Particle orientation segment | no | advanced |
| make_oris | Oritype segment in project | no | developer |
| oriops | Oritype segment in project | no | developer |
| oristats | Oritype segment in project | no | developer |
| print_project_field | Oritype segment in project | yes | standard |
| reconstruct3D_pcg | Oritype segment in project | no | developer |
| reextract | Oritype segment in project | no | advanced |
| replace_project_field | Oritype segment in project | yes | standard |
| selection | Oritype segment in project | no | advanced |
| vizoris | Oritype segment in project | no | advanced |

## Parameter `osmpl_fac`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| bootstrap_cavgs | Oversampling factor | no | developer |

## Parameter `outfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| bootstrap_rec3D | Resolution output prefix | no | developer |
| center | Output orientation and CTF parameter file | no | developer |
| dock_volpair | Output orientation and CTF parameter file | no | advanced |
| flex_analysis | Discrete-state project | no | advanced |
| make_oris | Output orientation and CTF parameter file | no | developer |
| mask | Output orientation and CTF parameter file | no | advanced |
| oriops | Output orientation and CTF parameter file | no | developer |
| simulate_particles | Output orientation and CTF parameter file | no | advanced |
| stackops | Output orientation and CTF parameter file | no | advanced |
| volcluster | Output orientation and CTF parameter file | no | advanced |
| volops | Output orientation and CTF parameter file | no | advanced |

## Parameter `outside`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract | Extract outside stage boundaries | no | advanced |
| reextract | Extract outside stage boundaries | no | advanced |

## Parameter `outstk`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| convert | Output stack name | no | advanced |
| ctfops | Output stack name | no | advanced |
| filter | Output stack name | no | advanced |
| graphene_subtr | Output stack name | no | advanced |
| ppca_denoise | Output stack name | no | advanced |
| reproject | Output stack name | no | advanced |
| scale | Output stack name | no | advanced |
| simulate_particles | Output stack name | no | advanced |
| stack | Output stack name | no | advanced |
| stackops | Output stack name | no | advanced |
| trajectory_denoise | Output stack name | no | advanced |

## Parameter `outvol`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| convert | Output volume name | no | advanced |
| dock_volpair | Output volume name | no | advanced |
| filter | Output volume name | no | advanced |
| flex_analysis | Output volume name | no | advanced |
| nu_filt3D | Output volume name | no | advanced |
| pdb2mrc | Output volume name | no | advanced |
| scale | Output volume name | no | advanced |
| simulate_nanoparticle | Output volume name | no | advanced |
| symmetrize_map | Output volume name | no | advanced |
| volops | Output volume name | no | advanced |

## Parameter `particle_density`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Particle density in micrographs | no | developer |
| pick | Particle density in micrographs | no | advanced |

## Parameter `pca_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cls_split | Class split embedding method | no | advanced |
| ppca_denoise | PCA methods: PPCA, PPCA plus residual kPCA, standard SVD PCA, kernel PCA, or diffusion maps | no | advanced |
| ppca_denoise_classes | PCA methods: PPCA, standard SVD PCA or kernel PCA | no | advanced |
| trajectory_denoise | PCA methods: diffusion maps, PPCA, PPCA plus residual kPCA, standard SVD PCA, or kernel PCA | no | advanced |

## Parameter `pcgop`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| reconstruct3D_pcg | PCG normal operator | no | developer |

## Parameter `pcontrast`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Input particle contrast | no | developer |
| extract | Input particle contrast | no | standard |
| mini_stream | Input particle contrast | no | advanced |
| pick | Input particle contrast | no | advanced |
| pick_extract | Input particle contrast | no | advanced |
| reextract | Input particle contrast | no | advanced |

## Parameter `pdbfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_stats | PDB | yes | standard |
| mask | PDB for 3D envelope masking | no | advanced |
| model_validate | PDB input coordinates file | yes | standard |
| pdb2mrc | PDB input coordinates file | yes | standard |
| simulate_nanoparticle | PDB | no | advanced |

## Parameter `pdbfile2`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_stats | PDB | no | advanced |

## Parameter `pdbfiles`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_rmsd | txt | yes | standard |
| core_atoms_analysis | txt | yes | standard |

## Parameter `pdbout`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pdb2mrc | Output PDB volume-centered molecule | no | advanced |

## Parameter `period`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| tseries_make_pickavg | Period for repeated averaging windows | no | advanced |

## Parameter `pgrp`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Point-group symmetry | yes | standard |
| abinitio3D_cavgs | Point-group symmetry | yes | standard |
| autorefine3D_nano | Point-group symmetry | yes | standard |
| bootstrap_rec3D | Point-group symmetry | yes | standard |
| cavgseoproc_nano | Point-group symmetry | yes | standard |
| cavgsproc_nano | Point-group symmetry | yes | standard |
| make_oris | Point-group symmetry | no | developer |
| oriops | Point-group symmetry | no | developer |
| oristats | Point-group symmetry | no | developer |
| pick_extract | Point-group symmetry | no | standard |
| ptclsproc_nano | Point-group symmetry | yes | standard |
| reconstruct3D | Point-group symmetry | yes | standard |
| refine3D | Point-group symmetry | yes | standard |
| refine3D_auto | Point-group symmetry | yes | standard |
| refine3D_multi | Point-group symmetry | yes | standard |
| refine3D_nano | Point-group symmetry | yes | standard |
| reproject | Point-group symmetry | yes | standard |
| simulate_particles | Point-group symmetry | yes | standard |
| symaxis_search | Point-group symmetry | yes | standard |
| symmetrize_map | Point-group symmetry | yes | standard |
| trajectory_reconstruct3D | Point-group symmetry | yes | standard |
| validate_cavgs_vs_model | Point-group symmetry | yes | standard |
| vizoris | Point-group symmetry | yes | standard |

## Parameter `pgrp_start`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Initital point-group symmetry | no | advanced |
| abinitio3D_cavgs | Initital point-group symmetry | no | advanced |

## Parameter `phrand`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | Phase randomization | no | advanced |

## Parameter `phshift_max`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Maximum CTF phase shift | no | developer |
| ctf_estimate | Maximum CTF phase shift | no | advanced |
| mini_stream | Maximum CTF phase shift | no | advanced |
| preproc | Maximum CTF phase shift | no | developer |
| preprocess | Maximum CTF phase shift | no | advanced |

## Parameter `phshift_min`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Minimum CTF phase shift | no | developer |
| ctf_estimate | Minimum CTF phase shift | no | advanced |
| mini_stream | Minimum CTF phase shift | no | advanced |
| preproc | Minimum CTF phase shift | no | developer |
| preprocess | Minimum CTF phase shift | no | advanced |

## Parameter `phshift_step`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | CTF phase-shift step | no | developer |
| ctf_estimate | CTF phase-shift step | no | advanced |
| mini_stream | CTF phase-shift step | no | advanced |
| preproc | CTF phase-shift step | no | developer |
| preprocess | CTF phase-shift step | no | advanced |

## Parameter `phshiftunit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | Phase-shift unit | no | advanced |

## Parameter `pick_roi`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Artefactual regions exclusion(new picker only) | no | developer |
| pick | Artefactual regions exclusion(new picker only) | no | advanced |

## Parameter `picker`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pick | Which picker to use | no | advanced |

## Parameter `pickrefs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| check_refpick | Stack of class-averages/reprojections for picking | yes | standard |
| master | 2D averages for use as picking references (optional) | no | developer |
| pick | Stack of class-averages/reprojections for picking | no | standard |
| pick_extract | Stack of class-averages/reprojections for picking | no | standard |

## Parameter `plaintexttab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | Plain text file of input parameters | no | advanced |

## Parameter `platonic`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| symmetry_test | Search for Platonic symmetries | no | advanced |

## Parameter `plot_key`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| print_project_field | plot plot_key on , sort on x | no | developer |

## Parameter `positive`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask2D | Consider only positive pixels | no | advanced |

## Parameter `postprocess`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| bootstrap_rec3D | Postprocess final map | no | developer |
| reconstruct3D | Postprocess final map | no | advanced |

## Parameter `ppca_kpca_resid_alpha`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise | Residual hybrid damping (0 => PPCA only; default 0.5) | no | advanced |
| trajectory_denoise | Residual hybrid damping (0 => PPCA only; default 0.5) | no | advanced |

## Parameter `pre_norm`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise_classes | Pre-normalize images | no | advanced |

## Parameter `preimage_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| flex_analysis | Diffusion-map output mode | no | advanced |

## Parameter `preimage_ndim`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| flex_analysis | Local-linear pre-image design dimension (default 2) | no | advanced |

## Parameter `prob_neigh_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | Prob-neigh neighborhood mode | no | advanced |
| refine3D_multi | Prob-neigh neighborhood mode | no | advanced |

## Parameter `projfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| estimate_lpstages | Project file | yes | standard |
| export_manifoldem_starproject | Project file | yes | standard |
| extract_subproj | Project file | yes | standard |
| extract_substk | Project file | yes | standard |
| map_params_from_den | Project file | no | advanced |
| match_cavgs | Project file | yes | standard |
| new_project | Project file | no | advanced |
| ptcl3D_state_consensus | Project file | yes | standard |
| replace_project_field | Project file | yes | standard |
| scale | Project file | no | advanced |
| selection | Project file | yes | standard |
| update_project | Project file | yes | standard |
| validate_projfile | Project file | yes | standard |
| write_mic_filetab | Project file | yes | standard |

## Parameter `projfile_den`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| map_params_from_den | Denoised child project file | yes | standard |

## Parameter `projfile_merged`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| merge_projects | Merged output project file | yes | standard |

## Parameter `projfile_orig`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| unbootstrap_cavgs | Original project file | yes | standard |

## Parameter `projfile_raw`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| map_params_from_den | Raw project file | yes | standard |

## Parameter `projfile_ref`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| match_cavgs | Reference project file | yes | standard |

## Parameter `projfile_target`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| replace_project_field | Another project file | yes | standard |

## Parameter `projname`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Project name | no | advanced |

## Parameter `projrec`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Projection-direction reconstruction | no | advanced |
| reconstruct3D | Projection-direction reconstruction | no | advanced |
| refine3D | Projection-direction reconstruction | no | advanced |

## Parameter `projstats`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oristats | Projection statistics | no | developer |

## Parameter `projtab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| merge_projects | Project file table | yes | standard |
| ptcl3D_state_consensus | Project file table | yes | standard |

## Parameter `prune`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cluster_cavgs | Automated particles pruning | no | advanced |
| map_cavgs_selection | Automated particles pruning | no | advanced |
| match_cavgs | Automated particles pruning | no | developer |
| model_cavgs_rejection | Automated particles pruning | no | advanced |
| ptcl3D_state_consensus | Automated particles pruning | no | advanced |
| select_clusters | Automated particles pruning | no | developer |
| selection | Automated particles pruning | no | advanced |

## Parameter `pspecsz`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ctf_estimate | Size of power spectrum | no | advanced |
| gen_pspecs_and_thumbs | Size of power spectrum | no | developer |
| motion_correct | Size of power spectrum | no | advanced |
| preprocess | Size of power spectrum | no | advanced |

## Parameter `ptcl_src`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Particle source | no | advanced |
| reconstruct3D | Particle source | no | advanced |
| refine3D | Particle source | no | advanced |
| refine3D_auto | Particle source | no | advanced |
| refine3D_multi | Particle source | no | advanced |

## Parameter `qsys_name`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Queue system kind | no | advanced |
| update_project | Queue system kind | no | advanced |

## Parameter `qsys_partition`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Name of SLURM/PBS/LSF partition | no | advanced |
| update_project | Name of SLURM/PBS/LSF partition | no | advanced |

## Parameter `qsys_qos`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Schedule priority | no | advanced |
| update_project | Schedule priority | no | advanced |

## Parameter `qsys_reservation`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Name of reserved partition | no | advanced |
| update_project | Name of reserved partition | no | advanced |

## Parameter `quality_context`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| model_cavgs_rejection | Class-average quality context | no | advanced |

## Parameter `quality_mode`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| model_cavgs_rejection | Class-average quality mode | no | advanced |

## Parameter `quality_model`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| model_cavgs_rejection | Class-average quality model | no | advanced |

## Parameter `real_filter`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | Real-space filter | no | advanced |

## Parameter `ref_ind`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| volanalyze | Reference volume index | no | advanced |

## Parameter `refine`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Refinement mode | no | advanced |
| abinitio2D_chunks | Refinement mode | no | advanced |
| refine3D | Refinement mode | no | advanced |

## Parameter `refs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_cavgs | Output 2D references | no | advanced |
| particle_sieving | Reference class averages to initialise size compatibility model | no | advanced |

## Parameter `reliongroups`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| export_relion | Number of Relion groups based on defocus | no | developer |

## Parameter `remap_cls`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_cavgs | Whether to remap 2D clusters | no | advanced |

## Parameter `res_target`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D_auto | Resolution target (in A) | no | advanced |

## Parameter `res_threshold`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| selection | Class resolution threshold | no | advanced |

## Parameter `rmsd_file`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_stats | bin | no | advanced |

## Parameter `roavg`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| estimate_diam | Rotationally average | no | advanced |
| stackops | Rotationally average | no | advanced |

## Parameter `rtol`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| reconstruct3D_pcg | PCG relative residual tolerance | no | developer |

## Parameter `scale`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| scale | Scaling ratio | no | advanced |

## Parameter `script`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| analysis2D_nano | Generate script for shared-mem exec on cluster | no | advanced |
| autorefine3D_nano | Generate script for shared-mem exec on cluster | no | advanced |
| cavgseoproc_nano | Generate script for shared-mem exec on cluster | no | developer |
| cavgsproc_nano | Generate script for shared-mem exec on cluster | no | developer |
| center2D_nano | Generate script for shared-mem exec on cluster | no | advanced |
| cluster2D_nano | Generate script for shared-mem exec on cluster | no | advanced |
| ptclsproc_nano | Generate script for shared-mem exec on cluster | no | developer |
| validate_cavgs_vs_model | Generate script for shared-mem exec on cluster | no | developer |

## Parameter `select_flag`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| select_clusters | flag to use for selection | no | developer |

## Parameter `sherr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| make_oris | Shift error half-width | no | developer |
| oriops | Shift error half-width | no | developer |
| simulate_particles | Shift error half-width | no | advanced |

## Parameter `sigma`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | sigma, for Gaussian generation | no | advanced |

## Parameter `sigma_est`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D | Sigma estimation method | no | advanced |
| refine3D | Sigma estimation method | no | advanced |

## Parameter `single_pass`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Coarse pass only | no | advanced |

## Parameter `skip_rejection`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| cluster_cavgs | Skip class-average rejection | no | advanced |

## Parameter `smpd`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_rmsd | Sampling distance | yes | standard |
| atoms_stats | Sampling distance | yes | standard |
| auto_spher_mask | Sampling distance | yes | standard |
| automask | Sampling distance | yes | standard |
| automask2D | Sampling distance | yes | standard |
| autorefine3D_nano | Sampling distance | yes | standard |
| cavgseoproc_nano | Sampling distance | yes | standard |
| cavgsproc_nano | Sampling distance | yes | standard |
| center | Sampling distance | yes | standard |
| check_refpick | Sampling distance | yes | standard |
| cif2mrc | Sampling distance | no | advanced |
| conv_atom_denoise | Sampling distance | yes | standard |
| core_atoms_analysis | Sampling distance | yes | standard |
| crys_score | Sampling distance | yes | standard |
| ctf_correct | Sampling distance | yes | standard |
| ctfops | Sampling distance | yes | standard |
| detect_atoms | Sampling distance | yes | standard |
| dock_volpair | Sampling distance | yes | standard |
| estimate_diam | Sampling distance | yes | standard |
| filter | Sampling distance | yes | standard |
| fsc | Sampling distance | yes | standard |
| fsc_area_score | Sampling distance | yes | standard |
| graphene_subtr | Sampling distance | yes | standard |
| icm2D | Sampling distance | yes | standard |
| icm3D | Sampling distance | yes | standard |
| import_cavgs | Sampling distance | yes | standard |
| import_movies | Sampling distance | yes | standard |
| import_particles | Sampling distance | yes | standard |
| import_trajectory | Sampling distance | yes | standard |
| mask | Sampling distance | yes | standard |
| master | Pixel size (A) | yes | standard |
| mini_stream | Sampling distance | yes | standard |
| model_validate | Sampling distance | yes | standard |
| noisevol | Sampling distance | yes | standard |
| normalize | Sampling distance | yes | standard |
| nu_filt3D | Sampling distance | yes | standard |
| pdb2mrc | Sampling distance | no | advanced |
| ppca_denoise | Sampling distance | yes | standard |
| ppca_volvar | Sampling distance | yes | standard |
| preproc | Sampling distance | yes | standard |
| print_dose_weights | Sampling distance | yes | standard |
| print_fsc | Sampling distance | no | advanced |
| print_magic_boxes | Sampling distance | yes | standard |
| ptclsproc_nano | Sampling distance | yes | standard |
| reproject | Sampling distance | yes | standard |
| scale | Sampling distance | yes | standard |
| simulate_movie | Sampling distance | yes | standard |
| simulate_nanoparticle | Sampling distance | yes | standard |
| simulate_particles | Sampling distance | yes | standard |
| split | Sampling distance | yes | standard |
| stack | Sampling distance | yes | standard |
| stackops | Sampling distance | yes | standard |
| symaxis_search | Sampling distance | yes | standard |
| symmetrize_map | Sampling distance | yes | standard |
| symmetry_test | Sampling distance | yes | standard |
| trajectory_denoise | Sampling distance | yes | standard |
| tsegmaps_core_finder | Sampling distance | yes | standard |
| tseries_import | Sampling distance | yes | standard |
| uniform_filter2D | Sampling distance | yes | standard |
| uniform_filter3D | Sampling distance | yes | standard |
| validate_cavgs_vs_model | Sampling distance | yes | standard |
| volanalyze | Sampling distance | yes | standard |
| volcluster | Sampling distance | yes | standard |
| volops | Sampling distance | no | advanced |

## Parameter `smpd_downscale`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Downscaled pixel size (A) | no | developer |
| motion_correct | Sampling distance after downscale | no | advanced |
| preproc | Sampling distance after downscale | no | standard |
| preprocess | Sampling distance after downscale | no | advanced |

## Parameter `smpd_target`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| model_validate | Target sampling distance | yes | standard |

## Parameter `snr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_movie | SNR | no | advanced |
| simulate_particles | SNR | yes | standard |
| stackops | Apply noise to give SNR | no | advanced |
| volops | SNR | no | advanced |

## Parameter `sort`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| print_project_field | sort oris on key | no | developer |

## Parameter `starfile`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| export_manifoldem_starproject | STAR-format file name | no | developer |
| import_particles | Particles Metadata starfile | no | standard |
| import_starproject | Metadata starfile | no | standard |

## Parameter `state`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Continuation state label | no | advanced |
| extract_substk | State index | no | advanced |
| oriops | State to modify | no | developer |
| postprocess | State to postprocess | no | advanced |
| prune_project | State index | no | developer |
| reconstruct3D_pcg | State to reconstruct | no | developer |
| selection | State number | no | advanced |
| stackops | State index | no | advanced |

## Parameter `stats`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| info_image | Output statistics | no | advanced |
| stackops | Provide statistics | no | advanced |

## Parameter `stepsz`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| trajectory_reconstruct3D | Time window size (# frames){500} | no | developer |

## Parameter `stk`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask2D | Particle image stack | yes | standard |
| binarize | Stack | no | developer |
| center | Particle image stack | no | developer |
| cluster_cavgs | Particle image stack | no | advanced |
| cluster_stack | Particle image stack | yes | standard |
| convert | Stack | no | advanced |
| ctfops | Particle image stack | no | advanced |
| estimate_diam | Particle image stack | yes | standard |
| filter | Stack to filter | no | advanced |
| icm2D | Odd stack | yes | standard |
| import_cavgs | Stack of class averages | yes | standard |
| import_particles | Stack of particles | no | advanced |
| import_trajectory | Particle image stack | yes | standard |
| map_cavgs_selection | Stack of cavgs to select from | no | advanced |
| mask | Particle image stack | no | advanced |
| match_stacks | Particle image stack | yes | standard |
| normalize | Stack to normalize | no | advanced |
| ppca_denoise | Stack to denoise | yes | standard |
| reimport_particles | Denoised particle stack | yes | standard |
| scale | Particle image stack | no | advanced |
| simulate_movie | Particle image stack | no | advanced |
| split | Particle image stack | yes | standard |
| stackops | Particle image stack | yes | standard |
| trajectory_denoise | Stack to denoise | yes | standard |
| trajectory_swap_stack | Particle image stack | yes | standard |
| uniform_filter2D | Odd stack | yes | standard |

## Parameter `stk2`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| icm2D | Even stack | yes | standard |
| map_cavgs_selection | Stack of selected cavgs | yes | standard |
| match_stacks | Second Particle image stack | yes | standard |
| uniform_filter2D | Even stack | yes | standard |

## Parameter `stk_backgr`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| graphene_subtr | background power spectra stack, eg NP_X_background_pspec.mrc | yes | standard |

## Parameter `stk_den`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | Denoised particle image stack | no | advanced |

## Parameter `stk_traj`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| graphene_subtr | Tracked NP image stack | yes | standard |

## Parameter `stktab`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | List of per-micrograph particle stacks | no | advanced |
| info_stktab | List of per-micrograph particle stacks | yes | standard |

## Parameter `stktab_den`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| import_particles | List of denoised particle stacks | no | advanced |

## Parameter `subprojname`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract_subproj | Subproject name | yes | standard |

## Parameter `symrnd`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Randomize over subgroubs of point-group | no | developer |

## Parameter `taper_edges`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| mask | Taper edges | no | advanced |

## Parameter `tester`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| track_particles | Write periodic tester-mode outputs | no | advanced |

## Parameter `thres`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask | Volume threshold | no | advanced |
| master | Distance threshold for peak picking(A) | no | developer |

## Parameter `tilt_thres`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| assign_optics_groups | Threshold for hierarchical clustering of beamtilts | no | developer |
| export_relion | Distance threshold | no | developer |
| preproc | Threshold for hierarchical clustering of beamtilts | no | developer |

## Parameter `tiltgroupmax`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| export_relion | Max movies in a tilt/optics group | no | developer |

## Parameter `time_per_image`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Time per image | no | advanced |
| update_project | Time per image | no | advanced |

## Parameter `tof`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| fractionate_movies | Final frame | no | advanced |

## Parameter `top`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| extract_subproj | To index | no | advanced |
| extract_substk | To index | no | advanced |
| scale | Last stack index | no | advanced |
| stackops | To particle index | no | advanced |
| trajectory_reconstruct3D | To particle index | no | developer |

## Parameter `total_dose`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| master | Total exposure dose (e/A2) | yes | standard |
| motion_correct | Total exposure dose (e/Ang^2) | no | standard |
| preproc | Total exposure dose (e/Ang^2) | no | standard |
| preprocess | Total exposure dose (e/Ang^2) | no | advanced |
| print_dose_weights | Total exposure dose (e/Ang^2) | yes | standard |

## Parameter `transp_pca`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| ppca_denoise_classes | transpose for pixel-wise learning | no | advanced |

## Parameter `trs`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| autorefine3D_nano | Maximum translational shift | no | advanced |
| center2D_nano | Maximum translational shift | no | advanced |
| cluster2D_nano | Maximum translational shift | no | advanced |
| dock_volpair | Maximum translational shift | no | advanced |
| motion_correct | Max shift per iter in pixels{10.} | no | advanced |
| preprocess | Max shift per iter in pixels{10.} | no | advanced |
| refine3D | Maximum translational shift | no | advanced |
| refine3D_nano | Maximum translational shift | no | advanced |
| simulate_movie | Maximum translational shift | no | advanced |
| tseries_make_pickavg | Max shift per iter in pixels{10.} | no | advanced |
| tseries_motion_correct | Max shift per iter in pixels{10.} | no | advanced |

## Parameter `trsstats`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oristats | Shift statistics | no | developer |

## Parameter `tseries`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| vizoris | Time series analysis | no | advanced |

## Parameter `update_frac`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D | Fractional update per iteration | no | advanced |
| refine3D_nano | Fractional update per iteration | no | advanced |

## Parameter `use_model`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| particle_sieving | Use model for class-average rejection in sieving | no | advanced |

## Parameter `user_account`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | User account name in SLURM/PBS/LSF | no | advanced |
| update_project | User account name in SLURM/PBS/LSF | no | advanced |

## Parameter `user_email`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | Your e-mail address | no | advanced |
| update_project | Your e-mail address | no | advanced |

## Parameter `user_project`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| new_project | User project name in SLURM/PBS/LSF | no | advanced |
| update_project | User project name in SLURM/PBS/LSF | no | advanced |

## Parameter `view_balance`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| flex_analysis | Correct uneven view distribution | no | advanced |

## Parameter `vis`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| info_image | Visualize image | no | advanced |
| stackops | Visualize stack images | no | advanced |

## Parameter `vol1`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio3D | Starting template volume | no | advanced |
| atoms_stats | Raw volume | yes | standard |
| auto_spher_mask | Odd volume | yes | standard |
| automask | Odd volume | yes | standard |
| autorefine3D_nano | FCC reference volume | yes | standard |
| binarize | Volume | no | developer |
| cavgseoproc_nano | Volume | yes | standard |
| cavgsproc_nano | Volume | yes | standard |
| center | Volume | no | developer |
| conv_atom_denoise | Volume | yes | standard |
| convert | Volume | no | advanced |
| detect_atoms | Volume | yes | standard |
| dock_volpair | Volume | yes | standard |
| filter | Volume to filter | no | advanced |
| flex_analysis | Mean volume | yes | standard |
| fsc | Odd volume | yes | standard |
| fsc_area_score | Odd volume | yes | standard |
| icm3D | Odd volume | yes | standard |
| mask | Volume | no | advanced |
| model_validate | Experimental volume | yes | standard |
| normalize | Volume to normalize | no | advanced |
| nu_filt3D | Odd volume | yes | standard |
| postprocess | Volume override | no | advanced |
| ppca_volvar | Volume | yes | standard |
| refine3D | Reference volume | no | advanced |
| refine3D_auto | Starting template volume | no | advanced |
| refine3D_nano | FCC reference volume | no | advanced |
| reproject | Volume | yes | standard |
| scale | Input volume | no | advanced |
| simulate_particles | Volume | yes | standard |
| symaxis_search | C1 Volume to identify symmetry axis of | yes | standard |
| symmetrize_map | Volume to symmetrize | yes | standard |
| symmetry_test | C1 Volume to identify symmetry of | yes | standard |
| trajectory_reconstruct3D | Mean volume for latent chunking | no | developer |
| uniform_filter3D | Odd volume | yes | standard |
| validate_cavgs_vs_model | Volume | yes | standard |
| volops | Volume | yes | standard |

## Parameter `vol2`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_stats | Connected components volume | yes | standard |
| automask | Even volume | yes | standard |
| dock_volpair | Volume | yes | standard |
| fsc | Even volume | yes | standard |
| fsc_area_score | Even volume | yes | standard |
| icm3D | Even volume | yes | standard |
| nu_filt3D | Even volume | yes | standard |
| uniform_filter3D | Even volume | yes | standard |

## Parameter `vol3`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| atoms_stats | Volume | no | advanced |

## Parameter `vol_dim`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| pdb2mrc | Simulated volume dimensions | no | advanced |

## Parameter `vol_even`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D_nano | Even volume | no | advanced |

## Parameter `vol_odd`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| refine3D_nano | Odd volume | no | advanced |

## Parameter `walltime`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| abinitio2D_chunks | Walltime | no | advanced |
| abinitio2D_stream | Walltime | no | developer |
| new_project | Walltime | no | advanced |
| particle_sieving | Walltime | no | advanced |
| pick_extract | Walltime | no | advanced |
| preproc | Walltime | no | developer |
| sieve_cavgs | Walltime | no | advanced |
| update_project | Walltime | no | advanced |

## Parameter `wcrit`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| motion_correct | Correlation to weights conversion scheme | no | advanced |
| tseries_make_pickavg | Correlation to weights conversion scheme | no | advanced |
| tseries_motion_correct | Correlation to weights conversion scheme | no | advanced |

## Parameter `which_iter`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| bootstrap_rec3D | Sigma iteration index | no | developer |

## Parameter `width`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| filter | Cosine low-pass filter falloff | no | advanced |
| mask | Falloff of inner mask | no | advanced |

## Parameter `winsz`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| automask2D | Window size for median filter | no | advanced |
| binarize | Half-window size | no | developer |
| cluster2D_nano | Half-window size | no | advanced |
| filter | Half-window size | no | advanced |
| pick | Window size for sauvola | no | advanced |
| stackops | Window size for local sdev estimation | no | advanced |

## Parameter `xdim`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_movie | x-dimension | no | advanced |

## Parameter `xmldir`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| assign_optics_groups | Directory containing per movie EPU XML files | no | developer |

## Parameter `xmlloc`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| export_relion | Pathname of EPU XML files | no | developer |

## Parameter `xsh`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| volops | Translation along x-axis | no | advanced |

## Parameter `ydim`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| simulate_movie | y-dimension | no | advanced |

## Parameter `ysh`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| volops | Translation along y-axis | no | advanced |

## Parameter `zero`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| oriops | Zero shifts | no | developer |

## Parameter `zsh`

| Used by program | Parameter label | Required in this program | Parameter visibility |
| --- | --- | --- | --- |
| volops | Translation along z-axis | no | advanced |
