!@descr: public parameters type and interfaces for parameter parsing and derivation phases
! for dummies:
! 1. To add a new command-line parameter, first declare it in type(parameters) below.
!    Put it in the same type-ordered/alphabetical section as similar parameters and
!    give it a declaration-time default when Fortran allows that.
! 2. To make the parameter parseable from the command line, register it in
!    src/main/params/simple_parameters_parse.f90 using the matching registry call:
!    add_char for fixed-length character or type(string), add_int for integer,
!    add_real for real, add_file for file arguments with format checks, and
!    add_dir for directory arguments.
! 3. If the parameter affects other parameters, derive those values in
!    src/main/params/simple_parameters_phases.f90. Keep the logic in the phase
!    that best matches the meaning of the parameter:
!    setup_execution_context, resolve_parameter_sources, derive_sampling_settings,
!    derive_parallel_settings, derive_image_settings, or
!    validate_parameter_consistency.
! 4. If the parameter needs sanity checks or mode restrictions, add them to
!    validate_parameter_consistency in simple_parameters_phases.f90.
! 5. If the parameter is a dynamic type(string) field that needs a non-empty or
!    special runtime default, set that in init_dynamic_defaults in
!    src/main/params/simple_parameters_core.f90.
! 6. The new constructor flow is:
!    reset defaults -> init_dynamic_defaults -> parse_inputs ->
!    setup_execution_context -> resolve_parameter_sources ->
!    derive_sampling_settings -> derive_parallel_settings ->
!    derive_image_settings ->
!    validate_parameter_consistency.
! 7. Rule of thumb:
!    declare in this file, parse in simple_parameters_parse.f90, derive or validate
!    in simple_parameters_phases.f90, and only touch simple_parameters_core.f90 for
!    dynamic-string defaults or shared utilities.
module simple_parameters
use simple_core_module_api
use simple_cmdline,    only: cmdline
use simple_ui_program, only: ui_program
use simple_ui,         only: get_prg_ptr
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
    character(len=3)          :: async='no'           !< asynchronous (yes|no){no}
    character(len=3)          :: atom_thres='yes'     !< do atomic thresholding or not (yes|no){yes}
    character(len=3)          :: autoscale='no'       !< automatic down-scaling(yes|no){yes}
    character(len=3)          :: avg='no'             !< calculate average (yes|no){no}
    character(len=3)          :: backgr_subtr='no'    !< Whether to perform micrograph background subtraction
    character(len=3)          :: balance='no'         !< Balance class populations to smallest selected
    character(len=3)          :: beamtilt='no'        !< use beamtilt values when generating optics groups
    character(len=3)          :: bin='no'             !< binarize image(yes|no){no}
    character(len=3)          :: boxes='no'           !< add box coordinates to JSON output(yes|no){no}
    character(len=3)          :: cavg_ini='no'        !< use class averages for initialization(yes|no){no}
    character(len=3)          :: cavg_ini_ext='no'    !< use class averages for (external) initialization(yes|no){no}
    character(len=3)          :: cavgw='no'           !< use class averages weights during 3D ab initio(yes|no){no}
    character(len=3)          :: center='yes'         !< center image(s)/class average(s)/volume(s)(yes|no){no}
    character(len=3)          :: center_pdb='no'      !< move PDB atomic center to the center of the box(yes|no){no}
    character(len=3)          :: chunk='no'           !< indicates whether we are within a chunk(yes|no){no}
    character(len=3)          :: classtats='no'       !< calculate class population statistics(yes|no){no}
    character(len=3)          :: clear='no'           !< clear exising processing upon start (stream)
    character(len=3)          :: combine_eo='no'      !< Whether combined e/o volumes have been used for alignment(yes|no){no}
    character(len=3)          :: continue='no'        !< continue previous refinement(yes|no){no}
    character(len=3)          :: ctfstats='no'        !< calculate ctf statistics(yes|no){no}
    character(len=3)          :: ctfpatch='yes'       !< whether to perform patched CTF estimation(yes|no){yes}
    character(len=3)          :: doprint='no'
    character(len=3)          :: downscale='yes'      !< wheter to downscale or not in motion correction
    character(len=3)          :: dw='yes'             !< Whether dose weighted micrographs will be generated, for use outside of the motion correction path(yes|no){yes}
    character(len=3)          :: dynreslim='no'       !< Whether the alignement resolution limit should be dynamic in streaming(yes|no){no}
    character(len=3)          :: empty3Dcavgs='yes'   !< whether 3D empty cavgs are okay(yes|no){yes}
    character(len=3)          :: envfsc='no'          !< envelope mask even/odd pairs for FSC calculation(yes|no){no}
    character(len=3)          :: eo_stage='yes'       !< Whether the last stage of abinitio2D uses a resolution limit determined with e/o pairs(yes|no){yes}
    character(len=3)          :: even='no'            !< even orientation distribution(yes|no){no}
    character(len=3)          :: extract='yes'        !< whether to extract particles after picking (streaming only)
    character(len=3)          :: extractfrommov='no'  !< whether to extract particles from the movie(yes|no){no}
    character(len=3)          :: fill_holes='no'      !< fill the holes post binarisation(yes|no){no}
    character(len=3)          :: fillin='no'          !< fillin particle sampling
    character(len=3)          :: force_lp_range='no'  !< force abinitio3D low-pass stages to use lpstart/lpstop directly(yes|no){no}
    character(len=3)          :: gauref='no'          !< Whether to apply a gaussian filter to the polar reference(yes|no){no}
    character(len=3)          :: guinier='no'         !< calculate Guinier plot(yes|no){no}
    character(len=3)          :: graphene_filt='no'   !< filter out graphene bands in correlation search
    character(len=3)          :: greedy_sampling='yes' !< greedy class sampling or not (referring to objective function)
    character(len=3)          :: hist='no'            !< whether to print histogram
    character(len=3)          :: icm='no'             !< whether to apply ICM filter to reference
    character(len=3)          :: incrreslim='no'      !< Whether to add ten shells to the FSC resolution limit
    character(len=3)          :: interactive='no'     !< Whether job is interactive
    character(len=3)          :: iterstats='no'       !< Whether to keep track alignment stats throughout iterations
    character(len=3)          :: json='no'            !< Print in json format (mainly for nice)
    character(len=3)          :: keepvol='no'         !< dev flag for preserving iterative volumes in refine3d
    character(len=3)          :: lam_anneal='no'      !< anneal lambda parameter
    character(len=3)          :: linethres='no'       !< whether to consider angular threshold in common lines (yes|no){no}
    character(len=3)          :: loc_sdev='no'        !< Whether to calculate local standard deviations(yes|no){no}
    character(len=10)         :: filt_mode='none'     !< filtering mode(none|uniform|fsc|nonuniform){none}
    character(len=3)          :: makemovie='no'
    character(len=3)          :: mcpatch='yes'        !< whether to perform patch-based alignment during motion correction
    character(len=3)          :: mcpatch_thres='yes'  !< whether to use the threshold for motion correction patch solution(yes|no){yes}
    character(len=3)          :: mirr='no'            !< mirror(no|x|y){no}
    character(len=3)          :: mirr_proj='no'       !< mirror projections(no|x|y){no}
    character(len=3)          :: mkdir='no'           !< make auto-named execution directory(yes|no){no}
    character(len=3)          :: ml_reg='yes'         !< apply ML regularization to class averages or volume
    character(len=3)          :: ml_reg_chunk='no'    !< apply ML regularization to class averages or volume in chunks
    character(len=3)          :: ml_reg_pool='no'     !< apply ML regularization to class averages or volume in pool
    character(len=3)          :: neg='no'             !< invert contrast of images(yes|no){no}
    character(len=3)          :: neigs_per='no'       !< using neigs as percentage of the total dimension(yes|no){no}
    character(len=3)          :: noise_norm ='yes'    !< image normalization based on background/foreground standardization(yes|no){yes}
    character(len=3)          :: norm='no'            !< do statistical normalisation avg
    character(len=3)          :: omit_neg='no'        !< omit negative pixels(yes|no){no}
    character(len=3)          :: outside='no'         !< extract boxes outside the micrograph boundaries(yes|no){no}
    character(len=3)          :: pad='no'
    character(len=3)          :: partition='no'
    character(len=3)          :: pca_img_ori='no'     !< original (no rotation/shifting within classes) ptcl stack to pca(yes|no){no}
    character(len=3)          :: pca_ori_stk='no'     !< output denoised particle stack in the original order and shifted/rotated back(yes|no){no}
    character(len=3)          :: phaseplate='no'      !< images obtained with Volta phaseplate(yes|no){no}
    character(len=3)          :: phrand='no'          !< phase randomize(yes|no){no}
    character(len=3)          :: pick_roi='yes'
    character(len=3)          :: platonic='yes'       !< platonic symmetry or not(yes|no){yes}
    character(len=3)          :: potts_prior='no'     !< ordered-label/Potts prior for nonuniform filter(yes|no){no}
    character(len=3)          :: postprocess='yes'    !< postprocess reconstruction output(yes|no){yes}
    character(len=3)          :: pre_norm='no'        !< pre-normalize images for PCA analysis
    character(len=3)          :: print_states='no'     !< exporting states during the refinement(yes|no){no}
    character(len=3)          :: proj_is_class='no'   !< intepret projection directions as classes
    character(len=3)          :: projstats='no'
    character(len=3)          :: prune='no'
    character(len=3)          :: prob_inpl='no'       !< probabilistic in-plane search in refine=neigh mode(yes|no){no}
    character(len=3)          :: randomise='no'       !< whether to randomise particle order
    character(len=3)          :: rank_cavgs='yes'     !< Whether to rank class averages(yes|no)
    character(len=3)          :: ranked_parts='yes'   !< generate ranked rather than balanced partitions in class sampling
    character(len=3)          :: recthres='no'        !< reconstruction angular threshold (yes|no){no}
    character(len=3)          :: reject_mics='no'     !< whether to reject micrographs based on ctfres/icefrac
    character(len=3)          :: remap_cls='no'
    character(len=3)          :: remove_chunks='yes'  !< whether to remove chunks after completion (yes|no){yes}
    character(len=3)          :: reset_boxfiles='no'  !< whether to remove existing boxfiles and set boxfile in current dir (yes|no){no}
    character(len=3)          :: restore_cavgs='yes'  !< Whether to restore images to class averages after orientation search (yes|no){yes}
    character(len=3)          :: ring='no'            !< whether to use ring in interactive stream pick
    character(len=3)          :: roavg='no'           !< rotationally average images in stack
    character(len=3)          :: transp_pca='no'
    character(len=3)          :: script='no'          !< do not execute but generate a script for submission to the queue
    character(len=3)          :: shbarrier='yes'      !< use shift search barrier constraint(yes|no){yes}
    character(len=3)          :: sort_asc='yes'       !< sort oris ascending
    character(len=3)          :: srch_oris='yes'      !< whether to search orientations in multivolume assignment(yes|no){yes} 
    character(len=3)          :: stream='no'          !< stream (real time) execution mode(yes|no){no}
    character(len=3)          :: stream2d='no'        !< indicates streaming 2D clustering(yes|no){no}
    character(len=3)          :: symrnd='no'          !< randomize over symmetry operations(yes|no){no}
    character(len=3)          :: taper_edges='no'     !< self-explanatory
    character(len=3)          :: tester='no'          !< write periodic tester-mode tracking outputs(yes|no){no}
    character(len=3)          :: tophat='no'          !< tophat filter(yes|no){no}
    character(len=3)          :: trail_rec='no'       !< trailing (weighted average) reconstruction when update_frac=yes 
    character(len=3)          :: trsstats='no'        !< provide origin shift statistics(yes|no){no}
    character(len=3)          :: tseries='no'         !< images represent a time-series(yes|no){no}
    character(len=3)          :: updated='no'         !< whether parameters has been updated
    character(len=3)          :: use_thres='yes'      !< Use contact-based thresholding(yes|no){yes}
    character(len=3)          :: vis='no'             !< visualise(yes|no)
    character(len=3)          :: verbose_exit='yes'   !< Whether to write a indicator file when task completes(yes|no){no}
    character(len=3)          :: volrec='yes'         !< volume reconstruction in 3D(yes|no){yes}
    character(len=3)          :: write_imgarr='no'    !< write out cavgs
    character(len=3)          :: zero='no'            !< zeroing(yes|no){no}
    ! files & directories strings in ascending alphabetical order
    ! default initialization is done in the init_strings method below
    type(string)              :: boxfile              !< file with EMAN particle coordinates(.txt)
    type(string)              :: boxtab               !< table (text file) of files with EMAN particle coordinates(.txt)
    type(string)              :: ciffile              !< xPDB/mmCIF file          
    type(string)              :: class_assignment     !< text file listing class ids assigned to a worker
    type(string)              :: classdoc             !< doc with per-class stats(.txt)
    type(string)              :: cwd
    type(string)              :: deftab               !< file with CTF info(.txt|.simple)
    type(string)              :: deselfile            !< file with indices to be deselected(.txt)
    type(string)              :: dir                  !< directory
    type(string)              :: dir_box
    type(string)              :: dir_exec             !< name of execution directory
    type(string)              :: dir_meta             !< grab xml files from here
    type(string)              :: dir_movies           !< grab mrc mrcs files from here
    type(string)              :: dir_prev             !< grab previous projects for streaming
    type(string)              :: dir_refine           !< refinement directory
    type(string)              :: dir_reject           !< move rejected files to here{rejected}
    type(string)              :: dir_select           !< move selected files to here{selected}
    type(string)              :: dir_target           !< put output here
    type(string)              :: exec_dir             !< auto-named execution directory
    type(string)              :: executable           !< name of executable
    type(string)              :: ext                  !< file extension{.mrc}
    type(string)              :: fbody                !< file body
    type(string)              :: filetab              !< list of files(.txt)
    type(string)              :: fname                !< file name
    type(string)              :: frcs                 !< binary file with per-class/proj Fourier Ring Correlations(.bin)
    type(string)              :: fsc                  !< binary file with FSC info{fsc_state01.bin}
    type(string)              :: gainref              !< gain reference for movie alignment
    type(string)              :: import_dir           !< dir to import .star files from for import_starproject
    type(string)              :: infile               !< file with inputs(.txt)
    type(string)              :: infile2              !< file with inputs(.txt)
    type(string)              :: last_prev_dir        !< last previous execution directory
    type(string)              :: msklist              !< table (text file) of mask volume files(.txt)
    type(string)              :: mskvols(MAXS)
    type(string)              :: niceserver           !< address and port of nice server for comms
    type(string)              :: optics_dir           !< directory containing stream optics
    type(string)              :: oritab               !< table  of orientations(.txt|.simple)
    type(string)              :: oritab2              !< 2nd table of orientations(.txt|.simple)
    type(string)              :: outdir               !< manually set output directory name
    type(string)              :: outfile              !< output document
    type(string)              :: outstk               !< output image stack
    type(string)              :: outvol               !< output volume{outvol.ext}
    type(string)              :: pdbfile              !< PDB file
    type(string)              :: pdbfile2             !< PDB file, another one
    type(string)              :: pdbfiles             !< list of PDB files
    type(string)              :: pdbout               !< PDB output file
    type(string)              :: pdfile
    type(string)              :: pickrefs             !< picking references
    type(string)              :: plaintexttab         !< plain text file of input parameters
    type(string)              :: prg                  !< SIMPLE program being executed
    type(string)              :: test                 !< SIMPLE TEST program being executed
    type(string)              :: projfile             !< SIMPLE *.simple project file
    type(string)              :: projfile_merged      !< merged SIMPLE *.simple project file output
    type(string)              :: projfile_optics      !< SIMPLE *.simple project file containing optics group definitions
    type(string)              :: projfile_ref         !< SIMPLE project file containing reference assignments
    type(string)              :: projfile_target      !< another SIMPLE *.simple project file
    type(string)              :: projname             !< SIMPLE  project name
    type(string)              :: refs                 !< initial2Dreferences.ext
    type(string)              :: refs_even
    type(string)              :: refs_odd
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
    type(string)              :: vol
    type(string)              :: vol_even             !< even reference volume
    type(string)              :: vol_odd              !< odd  reference volume
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
    character(len=4)          :: automatic='no'       !< automatic thres for edge detect (yes|no){no}
    character(len=5)          :: automsk='no'         !< automatic envelope masking (yes|tight|no){no}
    character(len=STDLEN)     :: center_type='mass'   !< Centering scheme used(mass|seg|params)
    character(len=STDLEN)     :: cls_init='ptcl'      !< Scheme to generate initial references for 2D analysis(ptcl|randcls|rand)
    character(len=STDLEN)     :: clustinds=''         !< comma-separated cluster indices
    character(len=STDLEN)     :: clust_crit='hybrid'  !< clustering criterion (fm|pow|hist|hybrid){hybrid}
    character(len=STDLEN)     :: cn_type='cn_std'     !< generalised coordination number (cn_gen) or stardard (cn_std)
    character(len=STDLEN)     :: ctf='no'             !< ctf flag(yes|no|flip)
    character(len=STDLEN)     :: detector='bin'       !< detector for edge detection (sobel|bin|otsu)
    character(len=STDLEN)     :: dfunit='microns'     !< defocus unit (A|microns){microns}
    character(len=5)          :: element ='     '     !< atom kind
    character(len=STDLEN)     :: filter='no'          !< filter type{no}
    character(len=STDLEN)     :: flag='dummy'         !< convenience flag for testing purpose
    character(len=STDLEN)     :: flipgain='no'        !< gain reference flipping (no|x|y|xy|yx)
    character(len=STDLEN)     :: inivol='sphere'      !< Different schemes for random volume generation(rand|rand_scaled|sphere){sphere}
    character(len=STDLEN)     :: multivol_mode='single' !< multivolume abinitio3D mode(single|independent|docked|input_oris_start|input_oris_fixed){single}
    character(len=STDLEN)     :: imgkind='ptcl'       !< type of image(ptcl|cavg|mic|movie){ptcl}
    character(len=STDLEN)     :: import_type='auto'   !< type of import(auto|mic|ptcl2D|ptcl3D){auto}
    character(len=STDLEN)     :: mcconvention='simple'!< which frame of reference convention to use for motion correction(simple|unblur|relion){simple}
    character(len=STDLEN)     :: multi_moldiams=''    !< list of molecular diameters to be used for multiple gaussian pick
    character(len=7)          :: objfun='euclid'      !< objective function(euclid|cc){euclid}
    character(len=STDLEN)     :: opt='bfgs'           !< optimiser (bfgs|simplex){bfgs}
    character(len=STDLEN)     :: oritype='ptcl3D'     !< SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D)
    character(len=STDLEN)     :: pca_mode='ppca' !< PCA mode(ppca|ppca_kpca_resid|pca_svd|kpca|diffusion_maps|steerable_diff_map){ppca}
    character(len=STDLEN)     :: steerable_denoise_mode='coeffproj' !< Steerable denoise mode(coeffproj|transport){coeffproj}
    character(len=STDLEN)     :: kpca_backend='nystrom' !< kPCA backend(exact|nystrom){nystrom}
    character(len=STDLEN)     :: kpca_ker='rbf'       !< kPCA kernel(rbf|cosine){rbf}
    character(len=STDLEN)     :: pcontrast='black'    !< particle contrast(black|white){black}
    character(len=STDLEN)     :: pickkind='gau'       !< Picking quasi-template(gau|ring|disc){gau}
    character(len=STDLEN)     :: pgrp='c1'            !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN)     :: pgrp_start='c1'      !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN)     :: phshiftunit='radians'!< additional phase-shift unit (radians|degrees){radians}
    character(len=STDLEN)     :: particle_density='optimal' !< particle density level (low|optimal|high){optimal}
    character(len=STDLEN)     :: picker='new'         !< which picker to use (old|new|segdiam){new}
    character(len=STDLEN)     :: plot_key=''          !< plot using plot_key on y axis, sort on x
    character(len=STDLEN)     :: protocol=''          !< generic option
    character(len=STDLEN)     :: ptclw='no'           !< use particle weights(yes|no){no}
    character(len=STDLEN)     :: qsys_name='local'    !< name of queue system (local|slurm|pbs|lsf)
    character(len=STDLEN)     :: qsys_partition2D=''  !< partition name for streaming 2D analysis
    character(len=STDLEN)     :: quality_mode='apply' !< class-average quality mode(apply|analyze|learn|promote){apply}
    character(len=STDLEN)     :: quality_model='chunk_default_v2' !< class-average quality model preset(chunk_default_v2|pool_default_v1){chunk_default_v2}
    character(len=STDLEN)     :: real_filter=''
    character(len=STDLEN)     :: refine='shc'         !< refinement mode(snhc|shc|neigh|shc_neigh|prob|prob_state|prob_neigh){shc}
    character(len=STDLEN)     :: refine_type='3D'     !< refinement mode(3D|2D|hybrid){3D}
    character(len=STDLEN)     :: select_flag='cluster' !< which flag to use for cluster selection (cluster|class){cluster}
    character(len=STDLEN)     :: sigma_est='group'    !< sigma estimation kind (group|global){group}
    character(len=STDLEN)     :: sort=''              !< key to sort oris on
    character(len=STDLEN)     :: speckind='sqrt'      !< power spectrum kind(real|power|sqrt|log|phase){sqrt}
    character(len=STDLEN)     :: split_mode='even'
    character(len=STDLEN)     :: startype=''          !< export type for STAR format (micrograph|select|extract|class2d|initmodel|refine3d|post){all}
    character(len=STDLEN)     :: stats='no'           !< provide statistics(yes|no|print){no}
    character(len=STDLEN)     :: tag=''               !< just a tag
    character(len=STDLEN)     :: wcrit = 'no'         !< correlation weighting scheme (softmax|zscore|sum|cen|exp|uniformno){sum}
    ! special integer kinds
    integer(kind(ENUM_ORISEG))     :: spproj_iseg = PTCL3D_SEG    !< sp-project segments that b%a points to
    integer(kind(ENUM_OBJFUN))     :: cc_objfun   = OBJFUN_EUCLID !< objective function(OBJFUN_CC = 0, OBJFUN_EUCLID = 1)
    integer(kind=kind(ENUM_WCRIT)) :: wcrit_enum  = CORRW_CRIT    !< criterium for correlation-based weights
    ! integer variables in ascending alphabetical order
    integer :: angstep=5
    integer :: binwidth=1          !< binary layers grown for molecular envelope(in pixels){1}
    integer :: box=0               !< square image size(in pixels)
    integer :: box_crop=0          !< square image size(in pixels), relates to Fourier cropped references
    integer :: box_croppd=0        !< square image size(in pixels), relates to Fourier cropped references and padded
    integer :: box_extract
    integer :: boxpd=0
    integer :: class=1             !< 2D class identity
    integer :: clip=0              !< clipped image box size(in pixels)
    integer :: clustind=0          !< cluster index
    integer :: cn=8                !< fixed std coord number for atoms in nanos
    integer :: cn_max=12           !< max std coord number for atoms in nanos
    integer :: cn_min=4            !< min std coord number for atoms in nanos
    integer :: cn_stop=10          !< rotational symmetry order stop index{10}
    integer :: cs_thres=2          !< contact score threshold for discarding atoms during autorefine3D_nano
    integer :: device=-1           !< Device id for OpenMP offloading
    integer :: edge=6              !< edge size for softening molecular envelope(in pixels)
    integer :: eer_fraction=20     !< # of EER raw frames to sum into a movie fraction
    integer :: eer_upsampling=1    !< Controls the movie desired output sampling: 1 = 4K x 4K pixels; 2 = 8K x 8K pixels
    integer :: extr_iter=1
    integer :: extr_lim=MAX_EXTRLIM2D
    integer :: find=1              !< Fourier index
    integer :: gauref_last_stage=0 !< When to switch off gaussian filtering{0}
    integer :: nframesgrp=0        !< # frames to group before motion_correct(Falcon 3){0}
    integer :: fromp=1             !< start ptcl index
    integer :: fromf=1             !< frame start index
    integer :: grow=0              !< # binary layers to grow(in pixels)
    integer :: hpind_fsc           !< high-pass Fourier index for FSC
    integer :: icm_stage=0
    integer :: iptcl=1
    integer :: istart=0
    integer :: job_memory_per_task2D=JOB_MEMORY_PER_TASK_DEFAULT
    integer :: kfromto(2)
    integer :: ldim(3)=0
    integer :: maxits=100          !< maximum # iterations
    integer :: maxits_glob=100     !< maximum # iterations, global
    integer :: maxits_between=30   !< maximum # iterations in between model building steps
    integer :: maxits_sh=60        !< maximum # iterations of shifting lbfgsb
    integer :: maxnchunks=0
    integer :: maxnptcls=0
    integer :: maxpop=0            !< max population of an optics group
    integer :: minits=0            !< minimum # iterations
    integer :: nboxes_max=0
    integer :: nchunks=0
    integer :: nchunksperset=0
    integer :: ncunits=0           !< # computing units, can be < nparts{nparts}
    integer :: ncls=500            !< # clusters
    integer :: ncls_sub=10         !< # sub-clusters
    integer :: nsubcls_min=3       !< minimum subclasses per parent class for class splitting
    integer :: nsubcls_max=10      !< maximum subclasses per parent class for class splitting
    integer :: ncls_start=10       !< minimum # clusters for 2D streaming
    integer :: ndiscrete=0         !< # discrete orientations
    integer :: neigs=0             !< # of eigenvectors {0=>auto for Nyström kPCA}
    integer :: kpca_nystrom_npts=512 !< # of Nyström landmarks
    integer :: kpca_nystrom_local_nbrs=96 !< max extra local support neighbors for Nyström reconstruction
    integer :: k_nn=10              !< local nearest-neighbor count for graph-based diffusion splitting
    integer :: steerable_nmodes=4   !< angular Fourier modes for cls_split steerable diffusion maps
    integer :: newbox=0            !< new box for scaling (by Fourier padding/clipping)
    integer :: nframes=0           !< # frames{30}
    integer :: ngrow=0             !< # of white pixel layers to grow in binary image
    integer :: niceprocid=0        !< # id of process in nice database
    integer :: ninipick=0          !< # of micrographs to run inipick preprocessing on in preprocess
    integer :: ninit=3             !< # of micrographs to use during diameter estimation global search
    integer :: nits_per_stage=5    !< # of iterations per stage
    integer :: nmics=0             !< # micrographs
    integer :: nmoldiams=1         !< # moldiams
    integer :: noris=0
    integer :: nparts=1            !< # partitions in distributed execution
    integer :: nparts_per_part=1   !< # partitions in distributed execution of balanced parts
    integer :: nparts_chunk=1      !< # partitions in chunks distributed execution
    integer :: nparts_pool =1      !< # partitions for pool distributed execution
    integer :: npeaks=NPEAKS_DEFAULT !< # of greedy subspace peaks to construct multi-neighborhood search spaces from
    integer :: npeaks_inpl=NPEAKS_INPL_DEFAULT !< # multi-neighborhood search peaks to refine with L-BFGS
    integer :: npix=0              !< # pixles/voxels in binary representation
    integer :: nptcls=1            !< # images in stk/# orientations in oritab
    integer :: nptcls_per_cls=500  !< # images in stk/# orientations in oritab
    integer :: nptcls_per_subcls=300 !< target particle count per subclass for class splitting
    integer :: nptcls_per_part=0   !< # particles per part in balanced selection
    integer :: nquanta=0           !< # quanta in quantization
    integer :: nran=0              !< # random images to select
    integer :: nrefs=100           !< # references used for picking{100}
    integer :: nrestarts=1
    integer :: nrots=0             !< number of in-plane rotations in greedy Cartesian search
    integer :: nsample=0           !< # particles to sample in refinement with fractional update
    integer :: nsample_max=0       !< maximum # particles to sample in refinement with fractional update
    integer :: nsample_start=0     !< # particles to sample in refinement with fractional update, lower bound
    integer :: nsample_stop=0      !< # particles to sample in refinement with fractional update, upper bound
    integer :: nsearch=40          !< # search grid points{40}
    integer :: nspace=2500         !< # projection directions
    integer :: nspace_sub=500      !< # projection directions in subspace
    integer :: nspace_max=1500     !< Maximum # of projection directions
    integer :: nstages=8           !< # low-pass limit stages
    integer :: nstates=1           !< # states to reconstruct
    integer :: nsym=1
    integer :: nthr=1              !< # OpenMP threads{1}
    integer :: nthr2D=1            !< # OpenMP threads{1}
    integer :: nthr_ini3D=1        !< # OpenMP threads{1}
    integer :: numlen=0            !< length of number string
    integer :: nxpatch=MC_NPATCH   !< # of patches along x for motion correction{5}
    integer :: nypatch=MC_NPATCH   !< # of patches along y for motion correction{5}
    integer :: offset=20           !< pixels offset{20}
    integer :: optics_offset=0
    integer :: part=1
    integer :: period=0           !< periodic window step in frames (0 means disabled)
    integer :: pid=0               !< process ID
    integer :: pftsz=0             !< Desired size of polarft_calc object (half the # of rotations)
    integer :: pspecsz=512         !< size of power spectrum(in pixels)
    integer :: ptcl=1
    integer :: ref_ind=0           !> reference index
    integer :: reliongroups=0
    integer :: shift_stage=0
    integer :: split_stage=7       !< splitting stage when multivol_mode==docked
    integer :: startit=1           !< start iterating from here
    integer :: stage=0
    integer :: state=1             !< state to extract
    integer :: stepsz=1            !< size of step{1}
    integer :: tofny=0
    integer :: top=1
    integer :: tof=0               !< end index
    integer :: vol_dim=0           !< input simulated pdb2mrc volume dimensions
    integer :: walltime=WALLTIME_DEFAULT  !< Walltime in seconds for workload management
    integer :: which_iter=0        !< iteration nr
    integer :: workers=0           !< # workers for distributed execution
    integer :: worker_nthr=0       !< # threads per worker for distributed execution
    integer :: xcoord=0            !< x coordinate{0}
    integer :: ycoord=0            !< y coordinate{0}
    integer :: xdim=0              !< x dimension(in pixles)
    integer :: xdimpd=0
    integer :: ydim=0              !< y dimension(in pixles)
    ! real variables in ascending alphabetical order
    real    :: amsklp=8.           !< low-pass limit for envelope mask generation(in A)
    real    :: angerr=0.           !< angular error(in degrees){0}
    real    :: angthres_mi_proj = ANGTHRES_MI_PROJ_DEFAULT !< 4 convergence checking
    real    :: ares=7.
    real    :: astigerr=0.         !< astigmatism error(in microns)
    real    :: astigthreshold=ASTIG_THRESHOLD !< ice fraction threshold{1.0}
    real    :: astigtol=0.05       !< expected (tolerated) astigmatism(in microns){0.05}
    real    :: athres=10.          !< angular threshold(in degrees)
    real    :: beta=0.0            !< Convenience parameter{0.0}
    real    :: bfac=200            !< bfactor for sharpening/low-pass filtering(in A**2){200.}
    real    :: bfacerr=50.         !< bfactor error in simulated images(in A**2){0}
    real    :: bw_ratio=0.3        !< ratio between foreground-background pixel desired in edge detection
    real    :: cenlp=20.           !< low-pass limit for binarisation in centering(in A){30 A}
    real    :: cs=2.7              !< spherical aberration constant(in mm){2.7}
    real    :: corr_thres=0.5      !< per-atom validation correlation threshold for discarding atoms
    real    :: ctfresthreshold=CTFRES_THRESHOLD !< ctf resolution threshold{30A}
    real    :: defocus=2.          !< defocus(in microns){2.}
    real    :: dferr=1.            !< defocus error(in microns){1.0}
    real    :: dfmax=DFMAX_DEFAULT !< maximum expected defocus(in microns)
    real    :: dfmin=DFMIN_DEFAULT !< minimum expected defocus(in microns)
    real    :: dfsdev=0.1
    real    :: dstep=0.
    real    :: dsteppd=0.
    real    :: e1=0.               !< 1st Euler(in degrees){0}
    real    :: e2=0.               !< 2nd Euler(in degrees){0}
    real    :: e3=0.               !< 3d Euler(in degrees){0}
    real    :: eps=0.5             !< learning rate
    real    :: eps_bounds(2) = [0.5,1.0]
    real    :: eullims(3,2)=0.
    real    :: extr_init=EXTRINITHRES !< initial extremal ratio (0-1)
    real    :: fny=0.
    real    :: frac=1.             !< fraction of ptcls(0-1){1}
    real    :: fraca=0.1           !< fraction of amplitude contrast used for fitting CTF{0.1}
    real    :: fracdeadhot=0.05    !< fraction of dead or hot pixels{0.01}
    real    :: frac_best=1.0       !< fraction of best particles to sample from per class when balance=yes
    real    :: frac_worst=1.0      !< fraction of worst particles to sample from per class when balance=yes
    real    :: frac_diam=0.5       !< fraction of atomic diameter
    real    :: fracsrch=0.9        !< fraction of serach space scanned for convergence
    real    :: fraction_dose_target=FRACTION_DOSE_TARGET_DEFAULT !< dose (in e/A2)
    real    :: frac_outliers=0.
    real    :: fraczero=0.
    real    :: ftol=1e-6
    real    :: gaufreq=-1.0         ! Full width at half maximum frequency for the gaussian filter
    real    :: hp=100.             !< high-pass limit(in A)
    real    :: hp_ctf_estimate=HP_CTF_ESTIMATE !< high-pass limit 4 ctf_estimate(in A)
    real    :: icefracthreshold=ICEFRAC_THRESHOLD !< ice fraction threshold{1.0}
    real    :: kv=300.             !< acceleration voltage(in kV){300.}
    real    :: kpca_cosine_weight_power=1.5 !< cosine local-weight sharpening power
    real    :: kpca_rbf_gamma=0.   !< RBF gamma (0=>auto)
    real    :: ppca_kpca_resid_alpha=0.5 !< damping for residual kPCA correction on top of PPCA
    real    :: lambda=1.0
    real    :: lam_bounds(2) = [0.05,1.0]
    real    :: lp=20.              !< low-pass limit(in A)
    real    :: lp2D=20.            !< low-pass limit(in A)
    real    :: lp_backgr=20.       !< low-pass for solvent blurring (in A)
    real    :: lp_discrete=20.     !< low-pass for discrete search used for peak detection (in A)
    real    :: lp_ctf_estimate=LP_CTF_ESTIMATE !< low-pass limit 4 ctf_estimate(in A)
    real    :: lpstart_nonuni= 30. !< optimization(search)-based low-pass limit lower bound
    real    :: lp_pick=PICK_LP_DEFAULT !< low-pass limit 4 picker(in A)
    real    :: lplim_crit=0.143    !< FSC criterion low-pass limit assignment(0.143-0.5){0.143}
    real    :: lplims2D(3)
    real    :: lpstart=0.          !< start low-pass limit(in A){15}
    real    :: lpstop=8.0          !< stop low-pass limit(in A){8}
    real    :: lpstart_ini3D=0.    !< start low-pass limit(in A){15}
    real    :: lpstop_ini3D=8.0    !< stop low-pass limit(in A){8}
    real    :: lpstop2D=8.0        !< stop low-pass limit(in A){8}
    real    :: lpthres=STREAM_RES_THRESHOLD
    real    :: max_dose=0.         !< maximum dose threshold (e/A2)
    real    :: max_rad=0.          !< particle longest  dim (in pixels)
    real    :: min_rad=100.        !< particle shortest dim (in pixels)
    real    :: moldiam=140.        !< molecular diameter(in A)
    real    :: moldiam_max=300.    !< upper bound molecular diameter(in A)
    real    :: moldiam_refine=0.   !< upper bound molecular diameter(in A)
    real    :: moldiam_ring=0.     !< ring picker diameter(in A)
    real    :: moment=0.
    real    :: msk=0.              !< mask radius(in pixels)
    real    :: msk_crop=0.         !< mask radius(in pixels)
    real    :: mskdiam=0.          !< mask diameter(in Angstroms)
    real    :: mskdiam_detect=0.   !< detect_atoms mask diameter(in Angstroms)
    real    :: mul=1.              !< origin shift multiplication factor{1}
    real    :: ndev=2.5            !< # deviations in one-cluster clustering
    real    :: ndev2D=CLS_REJECT_STD    !< # deviations for 2D class selection/rejection
    real    :: nsig=2.5            !< # sigmas
    real    :: osmpd=0.            !< target output pixel size
    real    :: overlap=0.9         !< required parameters overlap for convergence
    real    :: phranlp=35.         !< low-pass phase randomize(yes|no){no}
    real    :: prob_athres=10.     !< angle threshold for prob distribution samplings
    real    :: rec_athres=10.      !< angle threshold for reconstruction
    real    :: res_target = 3.     !< resolution target in A
    real    :: scale=1.            !< image scale factor{1}
    real    :: scale_movies=1.     !< movie scale factor
    real    :: sherr=0.            !< shift error(in pixels){2}
    real    :: sigma=1.0           !< for gaussian function generation {1.}
    real    :: smpd=1.3            !< sampling distance; same as EMANs apix(in A)
    real    :: smpd_downscale=SMPD4DOWNSCALE !< sampling distance for movie downscaling
    real    :: smpd_pickrefs       !< sampling distance of pickrefs
    real    :: smpd_target=0.5     !< target sampling distance; same as EMANs apix(in A) refers to paddep cavg/volume
    real    :: smpd_crop=2.        !< sampling distance; same as EMANs apix(in A) refers to cropped cavg/volume
    real    :: smpd_targets2D(2)
    real    :: snr=0.              !< signal-to-noise ratio
    real    :: snr_noise_reg=0.    !< signal to noise ratio of noise regularization
    real    :: tau=TAU_DEFAULT     !< for empirical scaling of cc-based particle weights
    real    :: tilt_thres=0.05
    real    :: thres=0.            !< threshold (binarisation: 0-1; distance filer: in pixels)
    real    :: thres_low=0.        !< lower threshold for canny edge detection
    real    :: thres_up=1.         !< upper threshold for canny edge detection
    real    :: tiltgroupmax=0
    real    :: total_dose
    real    :: trs=0.              !< maximum halfwidth shift(in pixels)
    real    :: update_frac = 1.
    real    :: ufrac_trec  = 1.    !< update frac trailing rec
    real    :: width=10.           !< falloff of mask(in pixels){10}
    real    :: winsz=KBWINSZ
    real    :: xsh=0.              !< x shift(in pixels){0}
    real    :: ysh=0.              !< y shift(in pixels){0}
    real    :: zsh=0.              !< z shift(in pixels){0}
    ! logical variables in (roughly) ascending alphabetical order
    logical :: l_autoscale       = .false.
    logical :: l_bfac            = .false.
    logical :: l_corrw           = .false.
    logical :: l_distr_worker    = .false.
    logical :: l_dose_weight     = .false.
    logical :: l_doshift         = .false.
    logical :: l_eer_fraction    = .false.
    logical :: l_envfsc          = .false.
    logical :: l_fillin          = .false.
    logical :: l_greedy_smpl     = .true.
    logical :: l_frac_best       = .false.
    logical :: l_frac_worst      = .false.
    logical :: l_update_frac     = .false.
    logical :: l_gauref          = .false.
    logical :: l_graphene        = .false.
    logical :: l_icm             = .false.
    logical :: l_incrreslim      = .false.
    logical :: l_lam_anneal      = .false.
    logical :: l_lpauto          = .false.
    logical :: l_lpset           = .false.
    logical :: l_ml_reg          = .true.
    logical :: l_noise_reg       = .false.
    logical :: l_neigh           = .false.
    logical :: l_phaseplate      = .false.
    logical :: l_potts_prior     = .false.
    logical :: l_prob_inpl       = .false.
    logical :: l_prob_align_mode = .false.
    logical :: l_sigma_glob      = .false.
    logical :: l_trail_rec       = .false.
    logical :: l_remap_cls       = .false.
    logical :: sp_required       = .false.
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
