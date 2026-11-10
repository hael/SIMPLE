!@descr: public parameters type and interfaces for parameter parsing and derivation phases
! New parameter: declare here (scripts/simple_args_generator.pl derives the accepted-argument list from these
! declarations), register in simple_parameters_parse.f90, derive/validate in simple_parameters_phases.f90
! (phase order: new()), expose in src/main/ui. Dynamic type(string) defaults: init_dynamic_defaults (core).
module simple_parameters
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_ui_program, only: ui_program
use simple_ui,         only: get_prg_ptr, get_test_prg_ptr
use simple_atoms,      only: atoms
use simple_decay_funs
use simple_parameters_registry, only: param_registry
implicit none

public :: parameters
private
#include "simple_local_flags.inc"

type :: parameters
    ! pointer 2 program UI
    type(ui_program), pointer :: ptr2prg => null()
    ! yes/no decision variables in ascending alphabetical order
    character(len=3)          :: acf='no'             !< calculate autocorrelation function(yes|no){no}
    character(len=3)          :: append='no'          !< append selection (yes|no){no}
    character(len=3)          :: assign_all='no'      !< 2D: converge only once every active particle has a class (fractional updates)(yes|no){no}
    character(len=3)          :: async='no'           !< asynchronous (yes|no){no}
    character(len=3)          :: atom_thres='yes'     !< do atomic thresholding or not (yes|no){yes}
    character(len=3)          :: autoscale='no'       !< automatic down-scaling(yes|no){yes}
    character(len=3)          :: avg='no'             !< calculate average (yes|no){no}
    character(len=3)          :: backgr_subtr='no'    !< Whether to perform micrograph background subtraction(yes|no){no}
    character(len=3)          :: beamtilt='no'        !< use beamtilt values when generating optics groups(yes|no){no}
    character(len=4)          :: bfac_model='line'    !< dose_filter3D: shrink each window's B-factor towards the dose curve or keep its own(line|free){line}
    character(len=3)          :: bin='no'             !< binarize image(yes|no){no}
    character(len=3)          :: boxes='no'           !< add box coordinates to JSON output(yes|no){no}
    character(len=3)          :: umap='yes'           !< flex_pca: UMAP plot coordinates and figures of the final embedding(yes|no){yes}
    character(len=3)          :: cavg_ini='no'        !< use class averages for initialization(yes|no){no}
    character(len=3)          :: cavg_ini_ext='no'    !< use class averages for (external) initialization(yes|no){no}
    character(len=3)          :: center='yes'         !< center image(s)/class average(s)/volume(s)(yes|no){no}
    character(len=3)          :: center_pdb='no'      !< move PDB atomic center to the center of the box(yes|no){no}
    character(len=3)          :: chunk='no'           !< indicates whether we are within a chunk(yes|no){no}
    character(len=3)          :: classtats='no'       !< calculate class population statistics(yes|no){no}
    character(len=3)          :: clear='no'           !< clear exising processing upon start (stream)(yes|no){no}
    character(len=3)          :: cluster_only='no'    !< perform clustering only (id junk but no rejection)in cluster_cavgs(yes|no){no}
    character(len=3)          :: combine_eo='no'      !< Whether combined e/o volumes have been used for alignment(yes|no){no}
    character(len=3)          :: continue='no'        !< continue previous refinement(yes|no){no}
    character(len=3)          :: ctfstats='no'        !< calculate ctf statistics(yes|no){no}
    character(len=3)          :: ctfpatch='yes'       !< whether to perform patched CTF estimation(yes|no){yes}
    character(len=3)          :: doprint='no' !< (yes|no){no}
    character(len=8)          :: dose_weighting='gg'  !< reextract: weight the movie fractions by Grant & Grigorieff or by the measured dose curve(gg|measured){gg}
    character(len=3)          :: downscale='yes'      !< wheter to downscale or not in motion correction(yes|no){yes}
    character(len=3)          :: dw='yes'             !< Whether dose weighted micrographs will be generated, for use outside of the motion correction path(yes|no){yes}
    character(len=3)          :: dynreslim='no'       !< Whether the alignement resolution limit should be dynamic in streaming(yes|no){no}
    character(len=3)          :: exit_collapse='no'   !< early exit on state collapse(yes|no){no}
    character(len=3)          :: envfsc='no'          !< envelope solvent correction for FSC calculation(yes|no){no}
    character(len=3)          :: pcg_order='no'       !< PCG order prior on the base solve: real-space ridge on what lacks NU evidence of order (solvent, detergent), then the same cold solve again(yes|no){no}
    character(len=3)          :: atm_prior='no'       !< PCG atomicity prior: penalize the map outside the span of the atom kernels of pdbfile, after a prior-free base solve(yes|no){no}
    character(len=3)          :: eo_stage='yes'       !< Whether the last stage of solve2D uses a resolution limit determined with e/o pairs(yes|no){yes}
    character(len=3)          :: eo_rebalance='yes'   !< reassign selected particle halves by direction and score(yes|no){yes}
    character(len=3)          :: eo_pool='yes'        !< trajectory_denoise flex: pool only frames of the same row parity(yes|no){yes}
    character(len=3)          :: even='no'            !< even orientation distribution(yes|no){no}
    character(len=3)          :: ewald='neg'          !< Ewald-sphere curvature of the wave models(neg|pos|no){neg}
    character(len=3)          :: extract='yes'        !< whether to extract particles after picking (streaming only)(yes|no){yes}
    character(len=3)          :: extractfrommov='no'  !< whether to extract particles from the movie(yes|no){no}
    character(len=3)          :: fill_holes='no'      !< fill the holes post binarisation(yes|no){no}
    character(len=3)          :: fillin='no'          !< fillin particle sampling (3D: every fifth iteration the particles never updated)(yes|no){no}
    character(len=3)          :: update_missing='no'  !< update only active particles with updatecnt==0(yes|no){no}
    character(len=3)          :: force_lp_range='no'  !< force solve3D low-pass stages to use lpstart/lpstop directly(yes|no){no}
    character(len=3)          :: gauref='no'          !< Whether to apply a gaussian filter to the polar reference(yes|no){no}
    character(len=3)          :: gen_model='no'       !< fit generative model and write generated/denoised particles(yes|no){no}
    character(len=3)          :: guinier='no'         !< calculate Guinier plot(yes|no){no}
    character(len=3)          :: graphene_filt='no'   !< filter out graphene bands in correlation search(yes|no){no}
    character(len=3)          :: greedy_sampling='yes' !< greedy class sampling or not (referring to objective function)(yes|no){yes}
    character(len=3)          :: hist='no'            !< whether to print histogram(yes|no){no}
    character(len=3)          :: preimage_auto='no'   !< flex_pca automatic state count: raise the npreimages ceiling and merge down(yes|no){no}
    character(len=3)          :: rec_states='yes'     !< flex_pca: reconstruct the state volumes, or deliver the embedding only(yes|no){yes}
    character(len=3)          :: incrreslim='no'      !< Whether to add ten shells to the FSC resolution limit(yes|no){no}
    character(len=3)          :: interactive='no'     !< Whether job is interactive(yes|no){no}
    character(len=3)          :: iterstats='no'       !< Whether to keep track alignment stats throughout iterations(yes|no){no}
    character(len=3)          :: json='no'            !< Print in json format (mainly for nice)(yes|no){no}
    character(len=3)          :: keepvol='no'         !< dev flag for preserving iterative volumes in refine3d(yes|no){no}
    character(len=3)          :: lam_anneal='no'      !< anneal lambda parameter(yes|no){no}
    character(len=3)          :: loc_sdev='no'        !< Whether to calculate local standard deviations(yes|no){no}
    character(len=STDLEN)     :: filt_mode='none'     !< filtering mode(none|uniform|fsc|nonuniform|nonuniform_lpset){none}
    character(len=3)          :: fit_phshift='no'     !< fit additive CTF phase shift(yes|no){no}
    character(len=3)          :: makemovie='no' !< (yes|no){no}
    character(len=3)          :: mcpatch='yes'        !< whether to perform patch-based alignment during motion correction(yes|no){yes}
    character(len=3)          :: mcpatch_thres='yes'  !< whether to use the threshold for motion correction patch solution(yes|no){yes}
    character(len=3)          :: memreport='no'       !< enable per-process runtime memory telemetry(yes|no){no}
    character(len=3)          :: mirr='no'            !< mirror(no|x|y){no}
    character(len=3)          :: mirr_proj='no'       !< generate paired 3D references by mirroring(no|yes){no}
    character(len=3)          :: mkdir='no'           !< make auto-named execution directory(yes|no){no}
    character(len=3)          :: unfil_unmasked='no' !< internal: export unmasked unfiltered half maps(yes|no){no}
    character(len=3)          :: ml_reg='yes'         !< apply ML regularization to class averages or volume(yes|no){yes}
    character(len=3)          :: ml_reg_chunk='no'    !< apply ML regularization to class averages or volume in chunks(yes|no){no}
    character(len=3)          :: ml_reg_pool='no'     !< apply ML regularization to class averages or volume in pool(yes|no){no}
    character(len=3)          :: neg='no'             !< invert contrast of images(yes|no){no}
    character(len=3)          :: noise_norm ='yes'    !< image normalization based on background/foreground standardization(yes|no){yes}
    character(len=3)          :: norm='no'            !< do statistical normalisation avg(yes|no){no}
    character(len=3)          :: nu_envmsk='no'       !< derive an envelope mask from the nonuniform filter evidence margin(yes|no){no}
    character(len=3)          :: omit_neg='no'        !< omit negative pixels(yes|no){no}
    character(len=3)          :: regpass='yes'        !< refine3D_auto: one global registration pass at the FSC=regpass_fsc band before the neighbourhood iterations(yes|no){yes}
    character(len=3)          :: outside='no'         !< extract boxes outside the micrograph boundaries(yes|no){no}
    character(len=3)          :: pad='no' !< (yes|no){no}
    character(len=3)          :: pdb_bfac='no'        !< simulate_nanoparticle pdbfile=: the B column broadens the element kernels(yes|no){no}
    character(len=3)          :: partition='no'       !< selection: states are partitions of classes, written to cls2D cluster(yes|no){no}
    character(len=3)          :: phrand='no'          !< phase randomize(yes|no){no}
    character(len=3)          :: pick_roi='yes' !< (yes|no){yes}
    character(len=3)          :: platonic='yes'       !< platonic symmetry or not(yes|no){yes}
    character(len=3)          :: postprocess='yes'    !< postprocess reconstruction output(yes|no){yes}
    character(len=3)          :: print_states='no'     !< exporting states during the refinement(yes|no){no}
    character(len=3)          :: proj_is_class='no'   !< intepret projection directions as classes(yes|no){no}
    character(len=3)          :: addon_diag='no'      !< solve3D_addon cohort-only diagnostic reconstruction(yes|no){no}
    character(len=3)          :: projstats='no' !< (yes|no){no}
    character(len=3)          :: prune='no' !< (yes|no){no}
    character(len=3)          :: prob_inpl='no'       !< probabilistic in-plane search in refine=neigh mode(yes|no){no}
    character(len=3)          :: rank_cavgs='yes'     !< Whether to rank class averages(yes|no)
    character(len=3)          :: ranked_parts='yes'   !< generate ranked rather than balanced partitions in class sampling(yes|no){yes}
    character(len=3)          :: recthres='no'        !< reconstruction angular threshold (yes|no){no}
    character(len=3)          :: reject_mics='no'     !< whether to reject micrographs based on ctfres/icefrac(yes|no){no}
    character(len=3)          :: reseed_states='no'   !< solve3D_cavgs: reseed a state that empties from the largest(yes|no){no}
    character(len=3)          :: remap_cls='no' !< (yes|no){no}
    character(len=3)          :: remove_chunks='yes'  !< whether to remove chunks after completion (yes|no){yes}
    character(len=3)          :: reset_boxfiles='no'  !< whether to remove existing boxfiles and set boxfile in current dir (yes|no){no}
    character(len=3)          :: restore_cavgs='yes'  !< Whether to restore images to class averages after orientation search (yes|no){yes}
    character(len=3)          :: ring='no'            !< whether to use ring in interactive stream pick(yes|no){no}
    character(len=3)          :: roavg='no'           !< rotationally average images in stack(yes|no){no}
    character(len=3)          :: trust_header='no'    !< whether to trust the header information in the input files(yes|no){no}
    character(len=3)          :: script='no'          !< do not execute but generate a script for submission to the queue(yes|no){no}
    character(len=3)          :: score_states='no'    !< report quality-model scores grouped by cls3D state(yes|no){no}
    character(len=3)          :: single_pass='no'     !< only run coarse pass of sieving(yes|no){no}
    character(len=3)          :: skip_rejection='no'  !< skip class-average rejection/update path(yes|no){no}
    character(len=3)          :: sort_asc='yes'       !< sort oris ascending(yes|no){yes}
    character(len=3)          :: stepwise='no'        !< stream pool 2D: import only the sets that reach the particle threshold(yes|no){no}
    character(len=3)          :: cc_emit_sigma='no' !< internal: emit Euclidean residual sigmas after CC assignment(yes|no){no}
    character(len=3)          :: cohort_sampling='no' !< internal: draw at a stage's first iteration, rescore that cohort after(yes|no){no}
    character(len=3)          :: scan_exit='no'       !< internal: multi-state convergence also on the mean % scanned above fracsrch alone(yes|no){no}
    character(len=3)          :: prob_post_diag='no'  !< internal: the dense probability table workers write per-particle pose-posterior statistics (prob_posterior_*.txt)(yes|no){no}
    character(len=3)          :: shift_prior='yes'    !< the probability tables' shift searches carry a Gaussian prior on the stored shift's displacement from the seed (the penalty follows the stored shift at every rotation), width from the previous iteration's shift increments; the window search's, off under the polar shift map (yes|no){yes}
    character(len=5)          :: balance='none'      !< fractional-update sampling units(none|class|cavg){none}
    character(len=16)         :: cont_prior='none' !< prior for continuous local refinement (none|shift|pose){none}
    real                      :: cont_prior_shift_sigma=0. !< fixed shift-prior standard deviation in native pixels; zero uses all-active increments
    real                      :: cont_prior_rot_sigma=0. !< fixed per-axis rotation-prior standard deviation in degrees; required for cont_prior=pose
    character(len=16)         :: cont_route='joint' !< internal: LM route of a continuous Cartesian pose solve(joint|shift_then_joint){joint}
    character(len=3)          :: sigma_commit_deferred='no' !< internal: caller commits the final canonical sigma update(yes|no){no}
    character(len=4)          :: ref_pose_init='none' !< external-reference pose initialization(cc|none){none}
    character(len=3)          :: stream='no'          !< stream (real time) execution mode(yes|no){no}
    character(len=3)          :: stream2d='no'        !< indicates streaming 2D clustering(yes|no){no}
    character(len=3)          :: symrnd='no'          !< randomize over symmetry operations(yes|no){no}
    character(len=3)          :: taper_edges='no'     !< self-explanatory(yes|no){no}
    character(len=3)          :: tester='no'          !< write periodic tester-mode tracking outputs(yes|no){no}
    character(len=3)          :: tophat='no'          !< tophat filter(yes|no){no}
    character(len=3)          :: trail_rec='no'       !< trailing (weighted average) reconstruction when update_frac=yes(yes|no){no}
    character(len=3)          :: trsstats='no'        !< provide origin shift statistics(yes|no){no}
    character(len=3)          :: tseries='no'         !< images represent a time-series(yes|no){no}
    character(len=3)          :: updated='no'         !< whether parameters has been updated; not a command-line key
    character(len=3)          :: use_model='yes'      !< use model for class-average rejection in sieve(yes|no){yes}
    character(len=3)          :: vis='no'             !< visualise(yes|no)
    character(len=3)          :: verbose_exit='yes'   !< Whether to write a indicator file when task completes(yes|no){no}
    character(len=3)          :: volrec='yes'         !< volume reconstruction in 3D(yes|no){yes}
    character(len=3)          :: wfloat16='no'        !< write MRC stacks in mode 12 float16 format(yes|no){no}
    character(len=3)          :: write_imgarr='no'    !< write out cavgs(yes|no){no}
    character(len=3)          :: zero='no'            !< zeroing(yes|no){no}
    character(len=3)          :: zerodose='no'        !< dose_filter3D: also write the zero-dose extrapolation(yes|no){no}
    ! files & directories strings in ascending alphabetical order
    ! default initialization is done in the init_strings method below
    type(string)              :: cache_dir            !< where the flex PCA plane cache and part files live (a fast local disk)
    type(string)              :: angtab               !< per-projection angles of a tilt series(.txt)
    type(string)              :: boxfile              !< file with EMAN particle coordinates(.txt)
    type(string)              :: boxtab               !< table (text file) of files with EMAN particle coordinates(.txt)
    type(string)              :: cavg_new_root        !< new root for relocated class-average data
    type(string)              :: cavg_old_root        !< old root stored for class-average data
    type(string)              :: ciffile              !< xPDB/mmCIF file
    type(string)              :: class_assignment     !< text file listing class ids assigned to a worker
    type(string)              :: classdoc             !< doc with per-class stats(.txt)
    type(string)              :: cwd !< not a command-line key
    type(string)              :: deftab               !< file with CTF info(.txt|.simple)
    type(string)              :: deselfile            !< file with indices to be deselected(.txt)
    type(string)              :: dir                  !< directory
    type(string)              :: dir_box
    type(string)              :: dir_exec             !< name of execution directory
    type(string)              :: dir_meta             !< grab xml files from here
    type(string)              :: dir_movies           !< grab mrc mrcs files from here
    type(string)              :: dir_preprocess       !< grab preprocessed files from here
    type(string)              :: dir_prev             !< grab previous projects for streaming
    type(string)              :: dir_refine           !< refinement directory
    type(string)              :: dir_reject           !< move rejected files to here{rejected}
    type(string)              :: dir_select           !< move selected files to here{selected}
    type(string)              :: dir_target           !< put output here
    type(string)              :: dose_curve           !< dose curve written by dose_filter3D(.txt)
    type(string)              :: exec_dir             !< auto-named execution directory; not a command-line key
    type(string)              :: executable           !< name of executable; not a command-line key
    type(string)              :: ext                  !< file extension{.mrc}
    type(string)              :: fbody                !< file body
    type(string)              :: filetab              !< list of files(.txt)
    type(string)              :: fname                !< file name
    type(string)              :: frcs                 !< binary file with per-class/proj Fourier Ring Correlations(.bin)
    type(string)              :: fsc                  !< binary file with FSC info{fsc_state01.bin}
    type(string)              :: pcg_mskfile          !< real-space [0,1] mask volume constraining the PCG solve (hard support projection P H P; experimental focused/support mode, pcg_priors_history.md dev item 5); accepted on both backends; spherical mskdiam support when absent
    type(string)              :: gainref              !< gain reference for movie alignment
    type(string)              :: import_dir           !< dir to import .star files from for import_starproject
    type(string)              :: infile               !< file with inputs(.txt)
    type(string)              :: infile2              !< file with inputs(.txt)
    type(string)              :: pindfile             !< distributed worker: particle rows written by the master (flex_pca, reconstruct3D)
    type(string)              :: keys
    type(string)              :: last_prev_dir        !< last previous execution directory; not a command-line key
    type(string)              :: mic_new_root         !< new root for relocated micrograph data
    type(string)              :: mic_old_root         !< old root stored for micrograph data
    type(string)              :: msklist              !< table (text file) of mask volume files(.txt)
    type(string)              :: mskvols(MAXS)
    type(string)              :: new_root             !< global new root for relocated project data
    type(string)              :: niceserver           !< address and port of nice server for comms
    type(string)              :: old_root             !< global old root stored in project data paths
    type(string)              :: optics_dir           !< directory containing stream optics
    type(string)              :: oritab               !< table  of orientations(.txt|.simple)
    type(string)              :: oritab2              !< 2nd table of orientations(.txt|.simple)
    type(string)              :: outdir               !< manually set output directory name
    type(string)              :: outfile              !< output document
    type(string)              :: outstk               !< output image stack
    type(string)              :: outvol               !< output volume{outvol.ext}
    type(string)              :: pdbfile              !< PDB file
    type(string)              :: pdbfile2             !< PDB file, another one
    type(string)              :: transfer_tab         !< detect_atoms: the transfer table (_transfer.tsv) to use instead of measuring it
    type(string)              :: pdbfiles             !< list of PDB files
    type(string)              :: pdbout               !< PDB output file
    type(string)              :: pickrefs             !< picking references
    type(string)              :: plaintexttab         !< plain text file of input parameters
    type(string)              :: prg                  !< SIMPLE program being executed
    type(string)              :: test                 !< SIMPLE TEST program being executed
    type(string)              :: projfile             !< SIMPLE *.simple project file
    type(string)              :: projfile_orig        !< original SIMPLE *.simple project file (unbootstrap source)
    type(string)              :: projfile_merged      !< merged SIMPLE *.simple project file output
    type(string)              :: projfile_optics      !< SIMPLE *.simple project file containing optics group definitions
    type(string)              :: projfile_out         !< output SIMPLE project file for non-destructive updates
    type(string)              :: projfile_ref         !< SIMPLE project file containing reference assignments
    type(string)              :: projfile_target      !< another SIMPLE *.simple project file
    type(string)              :: projfile_frozen      !< solve3D_addon: the frozen solution's SIMPLE *.simple project file
    type(string)              :: projname             !< SIMPLE  project name
    type(string)              :: projtab              !< table of SIMPLE *.simple project files
    type(string)              :: ptcl_new_root        !< new root for relocated particle data
    type(string)              :: ptcl_old_root        !< old root stored for particle data
    type(string)              :: refs                 !< initial2Dreferences.ext
    type(string)              :: refs_even !< not a command-line key
    type(string)              :: refs_odd !< not a command-line key
    type(string)              :: snapshot             !< path to write snapshot project file to
    type(string)              :: star_datadir         !< STAR-generated data directory
    type(string)              :: star_mic             !< STAR-formatted EM file (micrographs.star)
    type(string)              :: star_model           !< STAR-formatted EM file (model.star)
    type(string)              :: star_ptcl            !< STAR-formatted EM file (data.star)
    type(string)              :: starfile             !< STAR-formatted EM file (proj.star)
    type(string)              :: stk                  !< particle stack with all images(ptcls.ext)
    type(string)              :: stk2                 !< 2nd stack(in selection map: selected(cavgs).ext)
    type(string)              :: stk3                 !< 3d stack (in selection map (cavgs)2selectfrom.ext)
    type(string)              :: stk_backgr           !< stack with image for background subtraction
    type(string)              :: stktab               !< list of per-micrograph stacks
    type(string)              :: subprojname          !< SIMPLE  subproject name
    type(string)              :: verbose_exit_fname   !< File name of indicator file when task completes(TASK_FINISHED)
    type(string)              :: vol !< not a command-line key
    type(string)              :: vol_even             !< even reference volume
    type(string)              :: vol_new_root         !< new root for relocated volume data
    type(string)              :: vol_odd              !< odd  reference volume
    type(string)              :: vol_old_root         !< old root stored for volume data
    type(string)              :: vols(MAXS)
    type(string)              :: vols_even(MAXS)
    type(string)              :: vols_odd(MAXS)
    type(string)              :: worker_priority      !< priority to submit jobs with
    type(string)              :: worker_server        !< address and port of worker server for job submission
    type(string)              :: xmldir
    type(string)              :: xmlloc
    ! other character variables in ascending alphabetical order
    character(len=STDLEN)     :: algorithm=''         !< algorithm to be used
    character(len=STDLEN)     :: angastunit='degrees' !< angle of astigmatism unit (radians|degrees){degrees}
    character(len=STDLEN)     :: angconv='simple'     !< convention of the angtab angles (simple|tilt|genfire){simple}
    character(len=5)          :: automsk='no'         !< automatic envelope masking (yes|nu|tight|no){no}
    character(len=STDLEN)     :: center_type='mass'   !< Centering scheme used(mass|seg|params)
    character(len=STDLEN)     :: cls_init='ptcl'      !< Scheme to generate initial references for 2D analysis(ptcl|randcls|rand|prev){ptcl}; prev: solve2D seeded restart from the previous 2D clustering
    character(len=STDLEN)     :: clustinds=''         !< comma-separated cluster indices
    character(len=STDLEN)     :: clust_crit='hybrid'  !< clustering criterion (fm|pow|hist|hybrid){hybrid}
    character(len=STDLEN)     :: ctf='no'             !< ctf flag(yes|no|flip)
    character(len=STDLEN)     :: ctf_correct_mode='phaseflip' !< CTF correction operation(phaseflip|wiener)
    character(len=STDLEN)     :: dfunit='microns'     !< defocus unit (A|microns){microns}
    character(len=STDLEN)     :: element=''           !< atom kind
    character(len=2), allocatable :: species(:)       ! detect_atomic_species: species named by a list or compound selector, brightest first
    character(len=STDLEN)     :: filter='no'          !< filter type{no}
    character(len=STDLEN)     :: flag='dummy'         !< convenience flag for testing purpose
    character(len=STDLEN)     :: flipgain='no'        !< gain reference flipping (no|x|y|xy|yx)
    character(len=STDLEN)     :: graph='euc'          !< class split graph metric(euc|ori){euc}
    character(len=STDLEN)     :: multivol_mode='single' !< multivolume mode
    character(len=STDLEN)     :: imgkind='ptcl'       !< type of image(ptcl|cavg|mic|movie){ptcl}
    character(len=STDLEN)     :: import_type='auto'   !< type of import(auto|mic|ptcl2D|ptcl3D){auto}
    character(len=STDLEN)     :: mcconvention='simple'!< which frame of reference convention to use for motion correction(simple|unblur|relion){simple}
    character(len=STDLEN)     :: multi_moldiams=''    !< list of molecular diameters to be used for multiple gaussian pick
    character(len=7)          :: objfun='euclid'      !< objective function(euclid|cc){euclid}
    character(len=3)          :: inpl_cont='yes'      !< joint continuous in-plane refinement; internal, not a command-line key (no only under refine=cont)
    character(len=4)          :: pose_cont='no'       !< continuous Cartesian pose polishing(no|yes){no}
    character(len=STDLEN)     :: oritype='ptcl3D'     !< SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D)
    character(len=STDLEN)     :: rec_backend='gridding' !< 3D reconstruction backend(gridding|pcg){gridding}
    character(len=STDLEN)     :: m_estimator='no'     !< reconstruction M-estimator: flex weighs particles by the project's state weight set(no|flex){no}
    character(len=STDLEN)     :: rec_states_backend='gridding' !< flex_pca final state-map backend; independent of rec_backend, which governs the M-step(gridding|pcg){gridding}
    character(len=STDLEN)     :: state_placement='kcenter' !< flex_pca state targets at state_axis=0: diffusion k-center or equal-occupancy path(kcenter|equal_occ){kcenter}
    character(len=STDLEN)     :: pcontrast='black'    !< particle contrast(black|white){black}
    character(len=STDLEN)     :: pickkind='gau'       !< Picking quasi-template(gau|ring|disc){gau}
    character(len=STDLEN)     :: pgrp='c1'            !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN)     :: pgrp_start='c1'      !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN)     :: phshiftunit='radians'!< additional phase-shift unit (radians|degrees){radians}
    character(len=STDLEN)     :: particle_density='optimal' !< particle density level (low|optimal|high){optimal}
    character(len=STDLEN)     :: picker='new'         !< which picker to use (new|segdiam){new}
    character(len=STDLEN)     :: plot_key=''          !< plot using plot_key on y axis, sort on x
    character(len=STDLEN)     :: subset_mode='quality_coverage' !< select_micrographs objective(quality_coverage|particles)
    character(len=STDLEN)     :: subset_score='corr'   !< particle subset score(corr|corr_cart|none); none is micrograph-only
    character(len=STDLEN)     :: pose_policy='global' !< state-refinement pose policy(local|global){global}
    character(len=STDLEN)     :: prob_neigh_mode='state' !< prob_neigh neighborhood mode(state|geom|shc|snhc){state}
    character(len=STDLEN)     :: qsys_name='local'    !< name of queue system (local|coarray|slurm|pbs|lsf|sge)
    character(len=STDLEN)     :: qsys_partition2D=''  !< partition name for streaming 2D analysis
    character(len=STDLEN)     :: quality_mode='apply' !< class-average quality mode(apply|analyze|learn|evaluate|promote){apply}
    ! class-average quality context: sieve=small hard-gate-only chunks, chunk=larger learned chunks, pool=pre-3D pooled classes
    character(len=STDLEN)     :: quality_context='chunk' !< analysis hard-gate context(chunk|pool|sieve)
    ! class-average quality model preset(chunk100mics|sieve|pool){chunk100mics}
    character(len=STDLEN)     :: quality_model='chunk100mics'
    character(len=STDLEN)     :: real_filter=''
    character(len=STDLEN)     :: refine='shc'         !< refinement mode; 3D (shc|neigh|greedy|prob|prob_neigh|cont); 2D (snhc|snhc_smpl|snhc_smpl_many|prob|prob_neigh){shc}
    character(len=STDLEN)     :: select_flag='cluster' !< which flag to use for cluster selection (cluster|class){cluster}
    character(len=STDLEN)     :: sort=''              !< key to sort oris on
    character(len=STDLEN)     :: speckind='sqrt'      !< power spectrum kind(real|power|sqrt|log|phase){sqrt}
    character(len=STDLEN)     :: split_mode='even'
    character(len=STDLEN)     :: startype=''          !< export type for STAR format (micrograph|select|extract|class2d|initmodel|refine3d|post){all}
    character(len=STDLEN)     :: states=''             !< comma-separated list of states
    character(len=STDLEN)     :: state_pops=''         !< simulate_particles: comma-separated population fractions of the states{equal}
    character(len=STDLEN)     :: stats='no'           !< provide statistics(yes|no|print){no}
    character(len=STDLEN)     :: suite=''             !< one sub-suite of a unit-test area to run alone (simple_test_exec test=unit_<area> suite=<name>)
    character(len=STDLEN)     :: tag=''               !< just a tag
    character(len=STDLEN)     :: wave_model='rytov'   !< wave model of the image formation(projection|first_born|rytov|multislice){rytov}
    character(len=STDLEN)     :: wcrit = 'no'         !< correlation weighting scheme (softmax|zscore|sum|cen|exp|uniformno){sum}
    ! special integer kinds
    integer(kind(ENUM_ORISEG))     :: spproj_iseg = PTCL3D_SEG    !< sp-project segments that b%a points to; not a command-line key
    integer(kind(ENUM_OBJFUN))     :: cc_objfun   = OBJFUN_EUCLID !< objective function(OBJFUN_CC = 0, OBJFUN_EUCLID = 1); not a command-line key
    integer(kind=kind(ENUM_WCRIT)) :: wcrit_enum  = CORRW_CRIT    !< criterium for correlation-based weights; not a command-line key
    ! integer variables in ascending alphabetical order
    integer :: angstep=5
    integer :: binwidth=1          !< binary layers grown for molecular envelope(in pixels){1}
    integer :: box=0               !< square image size(in pixels)
    integer :: box_crop=0          !< square image size(in pixels), relates to Fourier cropped references
    integer :: box_croppd=0        !< square image size(in pixels), relates to Fourier cropped references and padded; not a command-line key
    integer :: box_extract
    integer :: box_coarse=0
    integer :: box_fine=0
    integer :: boxpd=0 !< not a command-line key
    integer :: class=1             !< 2D class identity
    integer :: clip=0              !< clipped image box size(in pixels)
    integer :: clustind=0          !< cluster index
    integer :: cn=8                !< fixed std coord number for atoms in nanos
    integer :: cn_max=12           !< max std coord number for atoms in nanos
    integer :: cn_min=4            !< min std coord number for atoms in nanos
    integer :: cn_stop=10          !< rotational symmetry order stop index{10}
    integer :: device=-1           !< Device id for OpenMP offloading
    integer :: edge=6              !< edge size for softening molecular envelope(in pixels)
    integer :: eer_fraction=20     !< # of EER raw frames to sum into a movie fraction
    integer :: eer_upsampling=1    !< Controls the movie desired output sampling: 1 = 4K x 4K pixels; 2 = 8K x 8K pixels
    integer :: extr_iter=1
    integer :: extr_lim=MAX_EXTRLIM2D
    integer :: find=1              !< Fourier index
    integer :: nframesgrp=0        !< # frames to group before motion_correct(Falcon 3){0}
    integer :: fromp=1             !< start ptcl index
    integer :: fromf=1             !< frame start index
    integer :: grow=0              !< # binary layers to grow(in pixels)
    integer :: halfwin=0           !< trajectory_denoise flex: window half-width in frames (0: the whole trajectory)
    integer :: hpind_fsc           !< high-pass Fourier index for FSC; not a command-line key
    integer :: iptcl=1
    integer :: istart=0
    integer :: job_memory_per_task2D=JOB_MEMORY_PER_TASK_DEFAULT
    integer :: kfromto(2)
    integer :: ldim(3)=0
    integer :: maxits=100          !< maximum # iterations
    integer :: maxits_pcg=2        !< maximum # PCG reconstruction iterations{2}
    integer :: maxits_atm=30       !< cap of the PCG iterations of the atomicity-prior solve{30}
    integer :: maxits_ml=0         !< coupled PCG iterations of the regularized system from the closed-form start; 0 = closed form only{0}
    integer :: maxits_glob=100     !< maximum # iterations, global
    integer :: maxits_between=30   !< maximum # iterations in between model building steps
    integer :: maxits_sh=60        !< maximum # iterations of shifting lbfgsb
    integer :: maxnchunks=0
    integer :: maxpop=0            !< max population of an optics group
    integer :: memreport_interval=30 !< minimum runtime memory report interval in seconds
    integer :: minits=0            !< minimum # iterations
    integer :: minpop=25           !< min # particles per subclass of cls_expansion; a class needs 2x to be split{25}
    integer :: nboxes_max=0
    integer :: nchunks=0
    integer :: nchunks_in_parallel=1 !< max # of chunks run in parallel (0 = all at once)
    integer :: nchunksperset=0
    integer :: ncunits=0           !< # computing units, can be < nparts{nparts}
    integer :: ncls=500            !< # clusters
    integer :: ncls_coarse=0       !< # coarse clusters
    integer :: ncls_fine=0         !< # fine clusters
    integer :: nclust=20           !< # class-average groups of balance=cavg sampling and cls_expansion{20}
    integer :: nptcls_coarse=0     !< # coarse particles
    integer :: nptcls_fine=0       !< # fine particles
    integer :: nsubcls_min=3       !< minimum subclasses per parent class for class splitting
    integer :: nsubcls_max=10      !< maximum subclasses per parent class for class splitting
    integer :: ncls_start=10       !< minimum # clusters for 2D streaming
    integer :: ndiscrete=0         !< # discrete orientations
    integer :: neigs=0             !< # of eigenvectors {0=>auto for Nyström kPCA}
    integer :: newbox=0            !< new box for scaling (by Fourier padding/clipping)
    integer :: nframes=0           !< # frames{30}
    integer :: ngrow=0             !< # of white pixel layers to grow in binary image
    integer :: niceprocid=0        !< # id of process in nice database
    integer :: nicedispid=0        !< # display id of process in nice database
    integer :: ninipick=0          !< # of micrographs to run inipick preprocessing on in preprocess
    integer :: nits_per_stage=5    !< # of iterations per stage
    integer :: nmics=0             !< # micrographs
    integer :: subset_max_nodes=100000 !< select_micrographs branch-and-bound node budget
    integer :: nmoldiams=1         !< # moldiams
    integer :: noris=0
    integer :: nparts=1            !< # partitions in distributed execution
    integer :: nparts_chunk=1      !< # partitions in chunks distributed execution
    integer :: nparts_pool =1      !< # partitions for pool distributed execution
    integer :: nparts3D=0          !< # partitions of the stream's 3D jobs (0: the stage's default)
    integer :: npeaks=NPEAKS_DEFAULT !< # of greedy subspace peaks to construct multi-neighborhood search spaces from
    integer :: npeaks_inpl=NPEAKS_INPL_DEFAULT !< # peaks refined with L-BFGS: the SHC strategies' in-plane search, the neighbourhood table per state, the 2D dense table (with a tenth of the classes)
    integer :: npix=0              !< # pixles/voxels in binary representation
    integer :: nptcls=1            !< # images in stk/# orientations in oritab
    integer :: nptcls_base=0       !< # particles in the initial solve3D snapshot
    integer :: nptcls_per_cls=500  !< # images in stk/# orientations in oritab
    integer :: nptcls_per_part=0   !< # particles per part in balanced selection
    integer :: nptcls_refs=0       !< # selected particles solve3D_pyramid makes its reference volumes from (0: its default)
    integer :: npreimages=0        !< # state volumes; 0 = DISCOVER it (over-provision, then two-gate merge)
    integer :: niter=5             !< # alternating low-rank covariance fitting iterations
    integer :: min_neff=2000       !< minimum effective particle count for covariance state kernels
    integer :: state_axis=0        !< covariance latent coordinate for state targets; 0=diffusion k-center over all comps
    integer :: nkern=0             !< # leading latent components used for state placement/kernel; 0=all (neigs)
    integer :: nbins=1             !< # kernel bandwidth bins for cross-validated state selection (1=off)
    integer :: column_separation=2 !< minimum grid separation between selected covariance columns
    integer :: n_probe_iters=0     !< EM/probe subspace-iteration refinements of the flex_pca column basis (0=off; commander sets 5)
    integer :: nquanta=0           !< # quanta in quantization
    integer :: nran=0              !< # random images to select
    integer :: nrefs=100           !< # references used for picking{100}
    integer :: nrestarts=1         !< # of restarts for eligible applications
    integer :: nrestarts_collapse=1   !< # of restarts for solve3D_cavgs when states collapse
    integer :: nrots=0             !< number of in-plane rotations in greedy Cartesian search
    integer :: nsample=0           !< # particles to sample in refinement with fractional update
    integer :: nsample_max=0       !< maximum # particles to sample in refinement with fractional update
    integer :: nsample_start=0     !< # particles to sample in refinement with fractional update, lower bound
    integer :: nsample_stop=0      !< # particles to sample in refinement with fractional update, upper bound
    integer :: nsample_coarse=0    !< # particles to sample in refinement with fractional update, coarse search
    integer :: nsample_fine=0      !< # particles to sample in refinement with fractional update, fine search
    integer :: nsnapshots=0        !< # cumulative project snapshots
    integer :: nspace=2500         !< # projection directions
    integer :: nspecies=0          !< detect_atomic_species: number of species, 0 = from the data{0}
    integer :: min_nbrs=3          !< detect_atoms: already found neighbours a weak atom of the gated residual stage needs, 0 = no requirement{3}
    integer :: nspace_sub=500      !< # projection directions in subspace
    integer :: nspace_pickrefs=0   !< # reprojections of the stream's picking references (0: the initial analysis' default)
    integer :: nstages=8           !< # low-pass limit stages
    integer :: nu_stage=0          !< solve3D stage NU filtering switches on at(0: the built-in default){0}
    integer :: gold_std_stage=0    !< first gold-standard stage in solve2D/solve3D (0: workflow default, 2/6){0}
    integer :: nstages_pickrefs=0  !< # solve3D_cavgs stages of the stream's picking references (0: the initial analysis' default)
    integer :: nstates=1           !< # states to reconstruct
    integer :: nstates_pickrefs=0  !< # solve3D_cavgs states of the stream's picking references (0: the initial analysis' default)
    integer :: nsym=1 !< not a command-line key
    integer :: nthr=1              !< # OpenMP threads{1}
    integer :: nthr2D=1            !< # OpenMP threads{1}
    integer :: nthr_ini3D=1        !< # OpenMP threads{1}
    integer :: nthr3D=0            !< # OpenMP threads of the stream's 3D jobs (0: the stage's default)
    integer :: nthr3D_pickrefs=0   !< # OpenMP threads of the solve3D_cavgs of the stream's picking references (0: the initial analysis' default)
    integer :: numlen=0            !< length of number string
    integer :: nxpatch=MC_NPATCH   !< # of patches along x for motion correction{5}
    integer :: nypatch=MC_NPATCH   !< # of patches along y for motion correction{5}
    integer :: offset=20           !< pixels offset{20}
    integer :: optics_offset=0
    integer :: osmpl_fac=2         !< class-average oversampling factor for bootstrap_cavgs
    integer :: part=1
    integer :: pcafit=0            !< flex_pca worker-side fit selector, set by the master (0=all|1=halfset A|2=halfset B)
    integer :: nfits=1             !< flex_pca worker: fits served by this round (1 single, 2 paired halves)
    integer :: period=0           !< periodic window step in frames (0 means disabled)
    integer :: pid=0               !< process ID; not a command-line key
    integer :: pftsz=0             !< Desired size of polarft_calc object (half the # of rotations)
    integer :: pspecsz=512         !< size of power spectrum(in pixels)
    integer :: ptcl=1 !< not a command-line key
    integer :: ref_ind=0           !> reference index
    integer :: reliongroups=0
    integer :: startit=1           !< start iterating from here
    integer :: stage=0             !< current stage carried to staged workflow workers
    integer :: state=1             !< state to extract
    integer :: stepf=3             !< Incremental frame step size{3}
    integer :: stepsz=1            !< size of step{1}
    integer :: tofny=0 !< not a command-line key
    integer :: top=1
    integer :: topk=3              !< refine3D/solve3D prob modes: poses a particle enters the iteration's reconstruction at, the committed draw and up to topk-1 extras, the refined entries of its row of highest posterior at the particle's own temperature within the insertion cone of the committed pose (rotation and shift moving the rim by less than the map's FSC=0.5 resolution), each with its share of that posterior; the search and the final reconstruction see the committed pose alone; 1 = the committed pose alone, the single-draw run bit for bit{3}
    integer :: tof=0               !< end index
    integer :: vol_dim=0           !< input simulated pdb2mrc volume dimensions
    integer :: walltime=WALLTIME_DEFAULT  !< Walltime in seconds for workload management
    integer :: which_iter=0        !< iteration nr
    integer :: workers=0           !< # workers for distributed execution
    integer :: worker_nthr=0       !< # threads per worker for distributed execution
    integer :: worker_server_nthr=0 !< # threads per worker of the worker_server (0: not given)
    integer :: xcoord=0            !< x coordinate{0}
    integer :: ycoord=0            !< y coordinate{0}
    integer :: xdim=0              !< x dimension(in pixles)
    integer :: xdimpd=0 !< not a command-line key
    integer :: ydim=0              !< y dimension(in pixles)
    ! real variables in ascending alphabetical order
    real    :: amsklp=8.           !< low-pass limit for NU/standalone envelope generation(in A)
    real    :: angerr=0.           !< angular error(in degrees){0}
    real    :: angthres_mi_proj = ANGTHRES_MI_PROJ_DEFAULT !< 4 convergence checking
    real    :: astigerr=0.         !< astigmatism error(in microns)
    real    :: astigthreshold=ASTIG_THRESHOLD !< ice fraction threshold{1.0}
    real    :: astigtol=0.05       !< expected (tolerated) astigmatism(in microns){0.05}
    real    :: athres=10.          !< angular threshold(in degrees)
    real    :: athres_cont=10.     !< local search radius: refine=cont's total rotation from the seed, prob_neigh_mode=geom's cone of directions about the previous one(in degrees){10}
    real    :: beta=0.0            !< Convenience parameter{0.0}
    real    :: bfac=200            !< bfactor for sharpening/low-pass filtering(in A**2){200.}
    real    :: bfacerr=50.         !< bfactor error in simulated images(in A**2){0}
    real    :: bw_ratio=0.3        !< ratio between foreground-background pixel desired in edge detection
    real    :: cenlp=20.           !< low-pass limit for binarisation in centering(in A){30 A}
    real    :: cont_conv_frac=CONT_CONV_FRAC    !< refine=cont convergence: stable fraction per state and streak coverage{0.9}
    real    :: cont_conv_rot=CONT_CONV_ROT_DEG  !< refine=cont convergence: largest rotation of a stable particle(in degrees){0.5}
    real    :: cont_conv_shift=CONT_CONV_SHIFT_A !< refine=cont convergence: largest shift of a stable particle(in A){1.0}
    real    :: cs=2.7              !< spherical aberration constant(in mm){2.7}
    real    :: ctfresthreshold=CTFRES_THRESHOLD !< ctf resolution threshold{30A}
    real    :: defocus=2.          !< defocus(in microns){2.}
    real    :: dferr=1.            !< defocus error(in microns){1.0}
    real    :: dfmax=DFMAX_DEFAULT !< maximum expected defocus(in microns)
    real    :: dfmin=DFMIN_DEFAULT !< minimum expected defocus(in microns)
    real    :: dstep=0. !< not a command-line key
    real    :: dsteppd=0. !< not a command-line key
    real    :: e1=0.               !< 1st Euler(in degrees){0}
    real    :: e2=0.               !< 2nd Euler(in degrees){0}
    real    :: e3=0.               !< 3d Euler(in degrees){0}
    real    :: envmsklp=ENVMSKLP_DEFAULT !< low-pass limit for envfsc density-envelope generation(in A)
    real    :: eps=0.5             !< learning rate
    real    :: eps_bounds(2) = [0.5,1.0]
    real    :: eullims(3,2)=0.
    real    :: extr_init=EXTRINITHRES !< initial extremal ratio (0-1)
    real    :: fny=0. !< not a command-line key
    real    :: fraca=0.1           !< fraction of amplitude contrast used for fitting CTF{0.1}
    real    :: fracdeadhot=0.05    !< fraction of dead or hot pixels{0.01}
    real    :: frac_best=1.0       !< fraction of best particles to sample from per class with sampling units (balance=class|cavg)
    real    :: frac_worst=1.0      !< fraction of worst particles to sample from per class in balanced selection
    real    :: frac_diam=0.5       !< fraction of atomic diameter
    real    :: fracsrch=0.9        !< fraction of serach space scanned for convergence
    real    :: fraction_dose_target=FRACTION_DOSE_TARGET_DEFAULT !< dose (in e/A2)
    real    :: frac_outliers=0.
    real    :: ftol=1e-6
    real    :: gaufreq=-1.0         ! Full width at half maximum frequency for the gaussian filter
    real    :: hp=100.             !< high-pass limit(in A)
    real    :: hp_ctf_estimate=HP_CTF_ESTIMATE !< high-pass limit 4 ctf_estimate(in A)
    real    :: icefracthreshold=ICEFRAC_THRESHOLD !< ice fraction threshold{1.0}
    real    :: kv=300.             !< acceleration voltage(in kV){300.}
    real    :: lambda=1.0
    real    :: lam_bounds(2) = [0.05,1.0]
    real    :: lp=20.              !< low-pass limit(in A)
    real    :: lp2D=20.            !< low-pass limit(in A)
    real    :: lp_backgr=20.       !< low-pass for solvent blurring (in A)
    real    :: lp_ctf_estimate=LP_CTF_ESTIMATE !< low-pass limit 4 ctf_estimate(in A)
    real    :: lp_pick=PICK_LP_DEFAULT !< low-pass limit 4 picker(in A)
    real    :: lplim_crit=0.143    !< FSC criterion low-pass limit assignment(0.143-0.5){0.143}
    real    :: lplims2D(3)
    real    :: lpstart=0.          !< start low-pass limit(in A){15}
    real    :: lpstop=8.0          !< stop low-pass limit(in A){8}
    real    :: lpstart_ini3D=0.    !< start low-pass limit(in A){15}
    real    :: lpstop_ini3D=8.0    !< stop low-pass limit(in A){8}
    real    :: lpstop_pickrefs=0.  !< solve3D_cavgs low-pass stop of the stream's picking references (in A; 0: the initial analysis' default)
    real    :: lpstop2D=8.0        !< stop low-pass limit(in A){8}
    real    :: lpstop_coarse=0.    !< stop low-pass limit(in A){0}
    real    :: lpexit=0.           !< internal: end the stage once FSC=0.5 resolves finer than this (in A; 0: off)
    real    :: lpstop_fine=0.      !< stop low-pass limit(in A){0}
    real    :: lpthres=STREAM_RES_THRESHOLD
    real    :: local_ang_bound=-1. !< refine3D_states local projection-angle override in degrees; negative is automatic
    real    :: local_shift_bound=-1. !< refine3D_states local shift override in pixels; negative is automatic
    real    :: min_state_frac=0.   !< flex_pca population floor: fraction of embedded particles every delivered state must hold; 0=off
    real    :: max_dose=0.         !< maximum dose threshold (e/A2)
    real    :: max_rad=0.          !< particle longest  dim (in pixels); not a command-line key
    real    :: min_rad=100.        !< particle shortest dim (in pixels); not a command-line key
    real    :: moldiam=140.        !< molecular diameter(in A)
    real    :: moldiam_max=500.    !< upper bound molecular diameter(in A)
    real    :: moldiam_refine=0.   !< upper bound molecular diameter(in A)
    real    :: moldiam_ring=0.     !< ring picker diameter(in A)
    real    :: model_gate_nsd=0.   !< trajectory_solve3D model gate: cut in robust spreads below the median (0 = K_MODEL)
    real    :: moment=0. !< not a command-line key
    real    :: msk=0.              !< mask radius(in pixels)
    real    :: msk_crop=0.         !< mask radius(in pixels)
    real    :: msk_box=0.          !< the box's mask radius (in pixels; the polar shift map's disk particles); not a command-line key
    real    :: msk_box_crop=0.     !< the box's mask radius at the cropped box (in pixels); not a command-line key
    real    :: mskdiam=0.          !< mask diameter(in Angstroms)
    real    :: mskdiam_detect=0.   !< detect_atoms mask diameter(in Angstroms)
    real    :: mul=1.              !< origin shift multiplication factor{1}
    real    :: noise_sdev=0.       !< additive Gaussian detector-noise standard deviation{0}
    real    :: ndev=2.5            !< # deviations in one-cluster clustering
    real    :: ndev2D=CLS_REJECT_STD    !< # deviations for 2D class selection/rejection
    real    :: nsig=2.5            !< # sigmas
    real    :: nu_msk_sig=3.0      !< NU evidence envelope mask threshold, in MADs above the solvent null{3.0}
    real    :: osmpd=0.            !< target output pixel size
    real    :: overlap=0.9         !< required parameters overlap for convergence
    real    :: phshift_max=180.    !< maximum fitted phase shift(in degrees){180}
    real    :: phshift_min=0.      !< minimum fitted phase shift(in degrees){0}
    real    :: phshift_step=10.    !< phase-shift grid step(in degrees){10}
    real    :: prob_athres=10.     !< angular bound of the SHC strategies' in-plane samplers (refine=neigh, the 2D *_smpl modes) and of the dense table's refined set (projs_ns with the mean pose change); the tables' draws take no threshold
    real    :: prob_beta=PROB_BETA !< set by the workflows, not a UI input: the iteration's tempering factor on the particle temperature of every posterior draw of the probability tables (in-plane, direction, class, neighbourhood); PROB_BETA flat, ramped to prob_beta_end by the annealing tail{0.3}
    real    :: prob_beta_end=0.    !< set by the workflows, not a UI input: the factor the run anneals to once its resolution stops gaining: a geometric ramp from prob_beta over the last PROB_BETA_TAIL iterations, then the run ends (1 = the physical posterior); 0 = prob_beta, no annealing; solve2D and solve3D set 1 on their last stage{0.0}
    real    :: rec_athres=10.      !< angle threshold for reconstruction
    real    :: regpass_fsc=0.143   !< refine3D_auto: FSC value of the startup pair defining the registration-pass band{0.143}
    real    :: res_target = 3.     !< resolution target in A
    real    :: res_threshold=-1.   !< resolution threshold in A (-1 means no threshold)
    real    :: rtol=0.             !< PCG relative residual tolerance; <=0 runs exactly maxits_pcg{0}
    real    :: lambda_atm=0.       !< atomicity prior strength relative to the data precision per voxel; 0 = the moment estimate{0}
    real    :: bfac_data=0.        !< atomicity prior: B (A^2) added to the narrowest kernel; fitted on the base map when not given
    real    :: rtol_atm=0.         !< atomicity-prior solve: stop at this true L2 relative residual; <=0 runs exactly maxits_atm{0}
    real    :: pcg_order_lambda=1.0 !< PCG order prior ridge coefficient relative to the data precision per voxel; chosen each iteration when not given (pcg_order=yes){auto}
    real    :: scale=1.            !< image scale factor{1}
    real    :: scale_movies=1.     !< movie scale factor; not a command-line key
    real    :: sherr=0.            !< shift error(in pixels){2}
    real    :: sigma=1.0           !< for gaussian function generation {1.}
    real    :: shift_sigma=10.     !< the polar shift map: sigma (in A) of the Gaussian prior on a particle's shift from its current origin, the picking error's width; a particle not yet centred, or whose prior width exceeds a third of the window, searches a disk of 3 sigma with a coarse map, the others the window, and the map initialises every polish; on in the probability tables and the classic polar modes with a shift search, off elsewhere and under objfun=cc, which refuses an explicit key(0: the window search){10}
    real    :: smpd=1.3            !< sampling distance; same as EMANs apix(in A)
    real    :: smpd_downscale      !< sampling distance for movie downscaling; native smpd when omitted
    real    :: smpd_target=0.5     !< target sampling distance; same as EMANs apix(in A) refers to paddep cavg/volume
    real    :: smpd_crop=2.        !< sampling distance; same as EMANs apix(in A) refers to cropped cavg/volume
    real    :: smpd_targets2D(2)
    real    :: snr=0.              !< signal-to-noise ratio
    real    :: snr_noise_reg=0.    !< signal to noise ratio of noise regularization
    real    :: tau=TAU_DEFAULT     !< fudge factor controlling ML regularization strength
    real    :: tilt_thres=0.05
    real    :: thres=0.            !< threshold (binarisation: 0-1; distance filer: in pixels)
    real    :: tiltgroupmax=0
    real    :: total_dose
    real    :: trs=0.              !< maximum halfwidth shift(in pixels)
    real    :: update_frac = 1.
    real    :: ufrac_trec  = 1.    !< explicit single-state trailing rec update-fraction override
    real    :: width=10.           !< falloff of mask(in pixels){10}
    real    :: winsz=KBWINSZ
    real    :: xsh=0.              !< x shift(in pixels){0}
    real    :: ysh=0.              !< y shift(in pixels){0}
    real    :: zsh=0.              !< z shift(in pixels){0}
    ! logical variables in (roughly) ascending alphabetical order
    logical :: l_autoscale       = .false. !< not a command-line key
    logical :: l_bfac            = .false. !< not a command-line key
    logical :: l_umap            = .false. !< not a command-line key
    logical :: l_distr_worker    = .false. !< not a command-line key
    logical :: l_dose_weight     = .false. !< not a command-line key
    logical :: l_doshift         = .false. !< not a command-line key
    logical :: l_eer_fraction    = .false. !< not a command-line key
    logical :: l_envfsc          = .false. !< not a command-line key
    logical :: l_pcg_order     = .false. !< not a command-line key
    logical :: l_pcg_order_lambda_auto = .true. !< pcg_order_lambda not given: chosen per state and iteration; not a command-line key
    logical :: l_atm_prior       = .false. !< not a command-line key
    logical :: l_bfac_data_fit   = .true.  !< bfac_data not given: fitted on the base map; not a command-line key
    logical :: l_fillin          = .false. !< not a command-line key
    logical :: l_assign_all      = .false. !< assign_all=yes under fractional updates: refine2D converges only when every active particle has a class; not a command-line key
    logical :: l_greedy_smpl     = .true. !< not a command-line key
    logical :: l_frac_best       = .false. !< not a command-line key
    logical :: l_update_frac     = .false. !< not a command-line key
    logical :: l_gauref          = .false. !< not a command-line key
    logical :: l_preimage_auto   = .false. !< not a command-line key
    logical :: l_incrreslim      = .false. !< not a command-line key
    logical :: l_lam_anneal      = .false. !< not a command-line key
    logical :: l_lpauto          = .false. !< not a command-line key
    logical :: l_lpset           = .false. !< not a command-line key
    logical :: l_ml_reg          = .true. !< not a command-line key
    logical :: l_gold_standard   = .false. !< solve2D worker stage >= its resolved gold_std_stage; not a command-line key
    logical :: l_gold_standard_resolved = .false. !< staged 2D reference policy supplied; not a command-line key
    logical :: l_noise_reg       = .false. !< not a command-line key
    logical :: l_neigh           = .false. !< not a command-line key
    logical :: l_pdb_bfac        = .false. !< not a command-line key
    logical :: l_species_list    = .false. !< element= is a comma-separated species list (detect_atomic_species); not a command-line key
    logical :: l_nonuniform      = .false. !< not a command-line key
    logical :: l_nonuniform_lpset = .false. !< not a command-line key
    logical :: l_regpass         = .true. !< not a command-line key
    logical :: l_prob_inpl       = .false. !< not a command-line key
    logical :: l_prob_align_mode = .false. !< not a command-line key
    logical :: l_cart_refine     = .false. !< a Cartesian (continuous-pose) pass (refine=cont); not a command-line key
    logical :: l_cont_polish     = .false. !< the Cartesian pass the polish schedules after a discrete one; not a command-line key
    logical :: l_cont_shift_first = .false. !< cont_route=shift_then_joint: the shift-only LM stage first; not a command-line key
    logical :: l_cohort_sampling = .false. !< not a command-line key
    logical :: l_scan_exit       = .false. !< scan_exit=yes: solve3D's multi-state shc stage 2 ends on the scan alone; not a command-line key
    logical :: l_prob_post_diag  = .false. !< prob_post_diag=yes: the dense table workers write pose-posterior statistics; not a command-line key
    logical :: l_shift_prior     = .false. !< shift_prior=yes; not a command-line key
    logical :: l_shift_map       = .false. !< shift_sigma > 0 under objfun=euclid: the polar shift map's key is on (the workflows pass it to their stages); not a command-line key
    logical :: l_centring        = .false. !< l_shift_map in a run that centres: a probability table or a classic polar mode with a shift search (the searches, the preparation and the consumers read this one); not a command-line key
    logical :: l_trail_rec       = .false. !< not a command-line key
    logical :: l_ufrac_trec_defined = .false. !< explicit ufrac_trec override was provided; not a command-line key
    logical :: l_m_estimator_flex   = .false. !< m_estimator=flex: reconstruction weighs particles by the state weight set; not a command-line key
    logical :: l_state_defined      = .false. !< state= was given: reconstruction restricted to that state; not a command-line key
    logical :: l_update_missing  = .false. !< not a command-line key
    logical :: l_remap_cls       = .false. !< not a command-line key
    logical :: sp_required       = .false. !< not a command-line key
contains
    procedure, private :: init_dynamic_defaults
    procedure, private :: parse_inputs
    procedure, private :: bind_input_carg
    procedure, private :: bind_input_files
    procedure, private :: bind_input_dirs
    procedure, private :: bind_input_iarg
    procedure, private :: bind_input_rarg
    procedure, private :: setup_execution_context
    procedure, private :: resolve_parameter_sources
    procedure, private :: derive_sampling_settings
    procedure, private :: derive_parallel_settings
    procedure, private :: derive_image_settings
    procedure, private :: validate_parameter_consistency
    procedure, private :: set_img_format
    procedure          :: new
    procedure          :: is_final_planned_iter
end type parameters

interface
    module subroutine init_dynamic_defaults(self)
        class(parameters), intent(inout) :: self
    end subroutine init_dynamic_defaults

    module subroutine parse_inputs(self, cline, cntfile, checkupfile)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        integer,           intent(inout) :: cntfile
        character(len=1),  intent(inout) :: checkupfile(:)
    end subroutine parse_inputs

    module subroutine bind_input_carg(self, reg)
        class(parameters),   target, intent(inout) :: self
        class(param_registry),       intent(inout) :: reg
    end subroutine bind_input_carg

    module subroutine bind_input_files(self, reg)
        class(parameters),   target, intent(inout) :: self
        class(param_registry),       intent(inout) :: reg
    end subroutine bind_input_files

    module subroutine bind_input_dirs(self, reg)
        class(parameters),   target, intent(inout) :: self
        class(param_registry),       intent(inout) :: reg
    end subroutine bind_input_dirs

    module subroutine bind_input_iarg(self, reg)
        class(parameters),   target, intent(inout) :: self
        class(param_registry),       intent(inout) :: reg
    end subroutine bind_input_iarg

    module subroutine bind_input_rarg(self, reg)
        class(parameters),   target, intent(inout) :: self
        class(param_registry),       intent(inout) :: reg
    end subroutine bind_input_rarg

    module subroutine new(self, cline, silent)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: silent
    end subroutine new

    module subroutine setup_execution_context(self, cline, silent)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical,           intent(in)    :: silent
    end subroutine setup_execution_context

    module subroutine resolve_parameter_sources(self, cline, cntfile, checkupfile)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        integer,           intent(in)    :: cntfile
        character(len=1),  intent(in)    :: checkupfile(:)
    end subroutine resolve_parameter_sources

    module subroutine derive_sampling_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
    end subroutine derive_sampling_settings

    module subroutine derive_parallel_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
    end subroutine derive_parallel_settings

    module subroutine derive_image_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
    end subroutine derive_image_settings

    module subroutine validate_parameter_consistency(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
    end subroutine validate_parameter_consistency

    module subroutine set_img_format(self, ext)
        class(parameters), intent(inout) :: self
        character(len=*),  intent(in)    :: ext
    end subroutine set_img_format

    module function is_final_planned_iter(self) result(final_planned)
        class(parameters), intent(in) :: self
        logical :: final_planned
    end function is_final_planned_iter

end interface

end module simple_parameters
