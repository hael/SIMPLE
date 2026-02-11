module simple_ui_program
use simple_core_module_api
use simple_ansi_ctrls
implicit none
! public :: ui_program, set_param, set_ui_params
! private
#include "simple_local_flags.inc"

! common input parameter type
type ui_param
    type(string) :: key
    type(string) :: keytype ! (binary|multi|num|str|file|dir)
    type(string) :: descr_short
    type(string) :: descr_long
    type(string) :: descr_placeholder
    type(string) :: gui_submenu
    type(string) :: active_flags
    type(string) :: exclusive_group
    type(string) :: cval_default
    real    :: rval_default = 0.
    logical :: required = .true.
    logical :: advanced = .true.
    logical :: online   = .false.
end type ui_param

! production-level program interface for simple_exec, single_exec & simple_stream executables
type :: ui_program
    type(string) :: name
    type(string) :: descr_short
    type(string) :: descr_long
    type(string) :: executable
    type(string) :: gui_submenu_list
    logical      :: advanced = .true.
    ! image input/output
    type(ui_param), allocatable :: img_ios(:)
    ! parameter input/output
    type(ui_param), allocatable :: parm_ios(:)
    ! alternative inputs
    type(ui_param), allocatable :: alt_ios(:)
    ! search controls
    type(ui_param), allocatable :: srch_ctrls(:)
    ! filter controls
    type(ui_param), allocatable :: filt_ctrls(:)
    ! mask controls
    type(ui_param), allocatable :: mask_ctrls(:)
    ! computer controls
    type(ui_param), allocatable :: comp_ctrls(:)
    ! sp_project required flag
    logical :: sp_required = .true.
    ! existence flag
    logical :: exists = .false.
  contains
    procedure          :: new
    procedure, private :: set_input_1
    procedure, private :: set_input_2
    procedure, private :: set_input_3
    generic            :: set_input => set_input_1, set_input_2, set_input_3
    procedure          :: print_ui
    procedure          :: print_cmdline
    procedure          :: print_prg_descr_long
    procedure          :: write2json
    procedure          :: get_name
    procedure          :: get_executable
    procedure          :: get_nrequired_keys
    procedure          :: get_required_keys
    procedure          :: requires_sp_project
    procedure, private :: kill
end type ui_program

! declare common params here, with name same as flag
type(ui_param), public :: algorithm
type(ui_param), public :: angerr
type(ui_param), public :: astigthreshold
type(ui_param), public :: astigtol
type(ui_param), public :: automsk
type(ui_param), public :: autosample
type(ui_param), public :: backgr_subtr
type(ui_param), public :: bfac
type(ui_param), public :: box
type(ui_param), public :: box_extract
type(ui_param), public :: cc_iters
type(ui_param), public :: center_pdb
type(ui_param), public :: clip
type(ui_param), public :: cls_init
type(ui_param), public :: clust_crit
type(ui_param), public :: cn
type(ui_param), public :: cn_max
type(ui_param), public :: cn_min
type(ui_param), public :: combine_eo
type(ui_param), public :: cs
type(ui_param), public :: ctf
type(ui_param), public :: ctf_yes
type(ui_param), public :: ctfpatch
type(ui_param), public :: ctfresthreshold
type(ui_param), public :: deftab
type(ui_param), public :: dferr
type(ui_param), public :: dfmax
type(ui_param), public :: dfmin
type(ui_param), public :: dir_movies
type(ui_param), public :: e1, e2, e3
type(ui_param), public :: eer_fraction
type(ui_param), public :: element
type(ui_param), public :: envfsc
type(ui_param), public :: eo
type(ui_param), public :: flipgain
type(ui_param), public :: frac
type(ui_param), public :: fraca
type(ui_param), public :: fraction_dose_target
type(ui_param), public :: frcs
type(ui_param), public :: gainref
type(ui_param), public :: graphene_filt
type(ui_param), public :: hp
type(ui_param), public :: icefracthreshold
type(ui_param), public :: icm
type(ui_param), public :: job_memory_per_task
type(ui_param), public :: kv
type(ui_param), public :: lp
type(ui_param), public :: lp_backgr
type(ui_param), public :: lp_pick
type(ui_param), public :: lplim_crit
type(ui_param), public :: lpstart_nonuni
type(ui_param), public :: lpthres
type(ui_param), public :: max_dose
type(ui_param), public :: max_rad
type(ui_param), public :: maxits
type(ui_param), public :: maxnchunks
type(ui_param), public :: mcconvention
type(ui_param), public :: mcpatch
type(ui_param), public :: mcpatch_thres
type(ui_param), public :: min_rad
type(ui_param), public :: mirr
type(ui_param), public :: ml_reg
type(ui_param), public :: ml_reg_chunk
type(ui_param), public :: ml_reg_pool
type(ui_param), public :: moldiam
type(ui_param), public :: moldiam_max
type(ui_param), public :: mskdiam
type(ui_param), public :: mskfile
type(ui_param), public :: mul
type(ui_param), public :: nboxes_max
type(ui_param), public :: nchunks
type(ui_param), public :: nchunksperset
type(ui_param), public :: ncls
type(ui_param), public :: ncls_start
type(ui_param), public :: neg
type(ui_param), public :: niceprocid
type(ui_param), public :: niceserver
type(ui_param), public :: nparts
type(ui_param), public :: nparts_chunk
type(ui_param), public :: nparts_pool
type(ui_param), public :: nptcls
type(ui_param), public :: nptcls_per_cls
type(ui_param), public :: nran
type(ui_param), public :: nrestarts
type(ui_param), public :: nsample
type(ui_param), public :: nsearch
type(ui_param), public :: nsig
type(ui_param), public :: nspace
type(ui_param), public :: nstates
type(ui_param), public :: nthr
type(ui_param), public :: numlen
type(ui_param), public :: nxpatch
type(ui_param), public :: nypatch
type(ui_param), public :: objfun
type(ui_param), public :: oritab
type(ui_param), public :: oritab2
type(ui_param), public :: oritype
type(ui_param), public :: outdir
type(ui_param), public :: outfile
type(ui_param), public :: outside
type(ui_param), public :: outstk
type(ui_param), public :: outvol
type(ui_param), public :: particle_density
type(ui_param), public :: pcontrast
type(ui_param), public :: pdbout
type(ui_param), public :: pgrp
type(ui_param), public :: pgrp_start
type(ui_param), public :: phaseplate
type(ui_param), public :: pick_roi
type(ui_param), public :: picker
type(ui_param), public :: pickrefs
type(ui_param), public :: projfile
type(ui_param), public :: projfile_merged
type(ui_param), public :: projfile_target
type(ui_param), public :: projname
type(ui_param), public :: prune
type(ui_param), public :: pspecsz
type(ui_param), public :: qsys_name
type(ui_param), public :: qsys_partition
type(ui_param), public :: qsys_qos
type(ui_param), public :: qsys_reservation
type(ui_param), public :: reject_cls
type(ui_param), public :: remap_cls
type(ui_param), public :: remove_chunks
type(ui_param), public :: script
type(ui_param), public :: sherr
type(ui_param), public :: sigma
type(ui_param), public :: sigma_est
type(ui_param), public :: smpd
type(ui_param), public :: smpd_downscale
type(ui_param), public :: smpd_target
type(ui_param), public :: star_datadir
type(ui_param), public :: star_mic
type(ui_param), public :: star_model
type(ui_param), public :: star_ptcl
type(ui_param), public :: starfile
type(ui_param), public :: startit
type(ui_param), public :: startype
type(ui_param), public :: stepsz
type(ui_param), public :: stk
type(ui_param), public :: stk2
type(ui_param), public :: stktab
type(ui_param), public :: time_per_image
type(ui_param), public :: total_dose
type(ui_param), public :: trs
type(ui_param), public :: tseries
type(ui_param), public :: update_frac
type(ui_param), public :: user_account
type(ui_param), public :: user_email
type(ui_param), public :: user_project
type(ui_param), public :: vol_dim
type(ui_param), public :: walltime
type(ui_param), public :: wcrit
type(ui_param), public :: width
type(ui_param), public :: wiener
type(ui_param), public :: winsz

interface set_param
   module procedure set_param_1
   module procedure set_param_2
end interface set_param

contains

    subroutine set_param_1( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                intent(in)    :: required
        real,                   intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%descr_short       = trim(descr_short)
        self%descr_long        = trim(descr_long)
        self%descr_placeholder = trim(descr_placeholder)
        self%required = required
        if( .not. self%required ) self%rval_default = default_value
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                intent(in)    :: required
        character(len=*),       intent(in)    :: default_value
        self%key               = trim(key)
        self%keytype           = trim(keytype)
        self%descr_short       = trim(descr_short)
        self%descr_long        = trim(descr_long)
        self%descr_placeholder = trim(descr_placeholder)
        self%required = required
        if( .not. self%required ) self%cval_default = trim(default_value)
    end subroutine set_param_2

    subroutine set_ui_params
        call set_param(algorithm,      'algorithm',       'multi',  'Algorithm for motion correction','Algorithm for motion correction(iso|patch|patch_refine){patch}','(iso|patch|patch_refine){patch}', .false.,'patch')
        call set_param(angerr,         'angerr',          'num',    'Rotation angle error half-width', 'Uniform rotation angle shift error half-width(in degrees)', 'rotation error in degrees', .false., 0.)
        call set_param(astigthreshold, 'astigthreshold',  'num',    'Astigmatism rejection threshold', 'Micrographs with astigmatism (%) above the threshold will be ignored from further processing{10.0}', 'Astigmatism threshold{10.0}', .false., 10.0)
        call set_param(astigtol,       'astigtol',        'num',    'Expected astigmatism', 'expected (tolerated) astigmatism(in microns){0.05}', 'in microns{0.05}',  .false., 0.05)
        call set_param(automsk,        'automsk',         'multi',  'Perform envelope masking', 'Whether to generate/apply an envelope mask(yes|tight|no){no}', '(yes|tight|no){no}', .false., 'no')
        call set_param(autosample,     'autosample',      'binary', 'Automated particles sampling scheme', 'Use automated sampling scheme to select particles subsets(yes|no){no}' , '(yes|no){no}', .false., 'no')
        call set_param(backgr_subtr,   'backgr_subtr',    'binary', 'Perform micrograph background subtraction(new picker only)', 'Perform micrograph background subtraction before picking/extraction(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(bfac,           'bfac',            'num',    'B-factor for sharpening','B-factor for sharpening in Angstroms^2', 'B-factor in Angstroms^2', .false., 200.)
        call set_param(box,            'box',             'num',    'Particle box size', 'Particle box size(in pixels)', '# pixels of box', .true., 0.)
        call set_param(box_extract,    'box_extract',     'num',    'Extracted particle image size', 'Extracted particle image size(in pixels)', 'Extracted particle image size', .false., 0.)
        call set_param(cc_iters,       'cc_iters',        'num',    'Number of correlation iterations before switching to ML', 'Number of correlation iterations before switching to ML{10}', '# of iterations{10}', .false., 10.)
        call set_param(center_pdb,     'center_pdb',      'binary', 'Whether to move the PDB atomic center to the center of the box', 'Whether to move the PDB atomic center to the center of the box (yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(clip,           'clip',            'num',    'Clipped box size', 'Target box size for clipping in pixels', 'in pixels', .false., 0.)
        call set_param(cls_init,       'cls_init',        'multi',  'Scheme for initial class generation', 'Initiate 2D analysis from raw images|random classes|noise images(ptcl|randcls|rand){ptcl}', '(ptcl|randcls|rand){ptcl}', .false., 'ptcl')
        call set_param(clust_crit,     'clust_crit',      'multi',   'Clustering criterion', 'Clustering criterion(sig|sig_clust|cc|res|hybrid){hybrid}', '(sig|sig_clust|cc|res|hybrid){hybrid}', .false., 'hybrid')
        call set_param(cn,             'cn',              'num',    'Fixed std coordination number', 'Minimum std cn to consider for dipole calc ', '8',  .false., 8.)
        call set_param(cn_max,         'cn_max',          'num',    'Maximum std coordination number', 'Maximum std cn to consider ', '12', .false., 12.)
        call set_param(cn_min,         'cn_min',          'num',    'Minimum std coordination number', 'Minimum std cn to consider ', '4',  .false., 4.)
        call set_param(combine_eo,     'combine_eo',      'binary', 'Whether e/o references are combined for final alignment(yes|no){no}', 'whether e/o references are combined for final alignment(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(cs,             'cs',              'num',    'Spherical aberration', 'Spherical aberration constant(in mm){2.7}', 'in mm{2.7}', .false., 2.7)
        call set_param(ctf,            'ctf',             'multi',  'CTF status', 'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){no}', '(yes|no|flip){no}', .true., 'no')
        call set_param(ctf_yes,        'ctf',             'multi',  'CTF status', 'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){yes}', '(yes|no|flip){yes}', .false., 'yes')
        call set_param(ctfpatch,       'ctfpatch',        'binary', 'Patch CTF estimation', 'Whether to perform patch CTF estimation(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(ctfresthreshold,'ctfresthreshold', 'num',    'CTF Resolution rejection threshold', 'Micrographs with a CTF resolution above the threshold (in Angs) will be ignored from further processing{50}', 'CTF resolution threshold(in Angstroms){50}', .false., 50.0)
        call set_param(deftab,         'deftab',          'file',   'CTF parameter file', 'CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format with dfx, dfy and angast values', '.simple|.txt parameter file', .false., 'deftab'//trim(METADATA_EXT))
        call set_param(dferr,          'dferr',           'num',    'Underfocus error half-width',  'Uniform underfoucs error half-width(in microns)',  'defocus error in microns', .false., 1.)
        call set_param(dfmax,          'dfmax',           'num',    'Expected maximum defocus', 'Expected maximum defocus in microns{5.0}', 'in microns{5.0}', .false., DFMAX_DEFAULT)
        call set_param(dfmin,          'dfmin',           'num',    'Expected minimum defocus', 'Expected minimum defocus in microns{0.2}', 'in microns{0.2}', .false., DFMIN_DEFAULT)
        call set_param(dir_movies,     'dir_movies',      'dir',    'Input movies directory', 'Where the movies to process are located or will squentially appear', 'e.g. /cryodata/', .true., 'preprocess/')
        call set_param(e1,             'e1',              'num',    'Rotation along Phi',  'Phi Euler angle',   'in degrees', .false., 0.)
        call set_param(e2,             'e2',              'num',    'Rotation along Theta','Theat Euler angle', 'in degrees', .false., 0.)
        call set_param(e3,             'e3',              'num',    'Rotation along Psi',  'Psi Euler angle',   'in degrees', .false., 0.)
        call set_param(eer_fraction,   'eer_fraction',    'num',    '# of EER frames to fraction together', 'Number of raw EER frames in a movie fraction', '# EER frames{20}', .false., 20.)
        call set_param(element,        'element',         'str',    'Atom element name: Au, Pt etc.', 'Atom element name: Au, Pt etc.', 'atom composition e.g. Pt', .true., '  ')
        call set_param(envfsc,         'envfsc',          'binary', 'Envelope mask e/o maps for FSC', 'Envelope mask even/odd pairs prior to FSC calculation(yes|no){yes}',  '(yes|no){yes}',  .false., 'yes')
        call set_param(eo,             'eo',              'binary', 'Gold-standard FSC for filtering and resolution estimation', 'Gold-standard FSC for filtering and resolution estimation(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(flipgain,       'flipgain',        'multi',  'Flip the gain reference', 'Flip the gain reference along the provided axis(no|x|y|xy|yx){no}', '(no|x|y|xy|yx){no}', .false., 'no')
        call set_param(frac,           'frac',            'num',    'Fraction of particles to include', 'Fraction of particles to include based on spectral score (median of FRC between reference and particle)', 'fraction of particles(0.1-0.9){1.0}', .false., 1.0)
        call set_param(fraca,          'fraca',           'num',    'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.1}', 'fraction{0.1}', .false., 0.1)
        call set_param(fraction_dose_target,'fraction_dose_target','num','EER fraction dose target (e/Ang^2)', 'EER fraction dose target, used to determine how many EER frames are included in each movie fraction(e/Ang^2)', 'in e/Ang^2', .false., 1.)
        call set_param(frcs,           'frcs',            'str',    'Projection FRCs file', 'Projection FRCs file', 'e.g. frcs.bin', .false., '')
        call set_param(gainref,        'gainref',         'file',   'Gain reference', 'Gain reference image', 'input image e.g. gainref.mrc', .false., '')
        call set_param(graphene_filt,  'graphene_filt',   'binary', 'Omit graphene bands from corr calc', 'Omit graphene bands from corr calc(yes|no){no}',  '(yes|no){no}',  .false., 'no')
        call set_param(hp,             'hp',              'num',    'High-pass limit', 'High-pass resolution limit', 'high-pass limit in Angstroms', .false., 100.)
        call set_param(icefracthreshold,'icefracthreshold','num',   'Ice Fraction rejection threshold', 'Micrographs with an ice ring/1st pspec maxima fraction above the threshold will be ignored from further processing{1.0}', 'Ice fraction threshold{1.0}', .false., 1.0)
        call set_param(icm,            'icm',             'binary', 'Whether to perform ICM filtering of reference(s)', 'Whether to perform ICM filtering of reference(s)(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(job_memory_per_task, 'job_memory_per_task','str', 'Memory per computing node', 'Memory in MB per part/computing node in distributed execution{16000}', 'MB per part{16000}', .false., 16000.)
        call set_param(kv,             'kv',              'num',    'Acceleration voltage', 'Acceleration voltage in kV{300}', 'in kV{300}', .false., 300.)
        call set_param(lp,             'lp',              'num',    'Low-pass limit', 'Low-pass resolution limit', 'low-pass limit in Angstroms', .false., 20.)
        call set_param(lp_backgr,      'lp_backgr',       'num',    'Background low-pass resolution', 'Low-pass resolution for solvent blurring', 'low-pass limit in Angstroms', .false., 20.)
        call set_param(lp_pick,        'lp_pick',         'num',    'Low-pass limit for picking', 'Low-pass limit for picking in Angstroms{20}', 'in Angstroms{20}', .false., 20.)
        call set_param(lplim_crit,     'lplim_crit',      'num',    'Low-pass limit FSC criterion', 'FSC criterion for determining the low-pass limit(0.143-0.5){0.143}', 'low-pass FSC criterion(0.143-0.5){0.143}', .false., 0.143)
        call set_param(lpthres,        'lpthres',         'num',    'Resolution rejection threshold', 'Classes with lower resolution are iteratively rejected in Angstroms{30}', 'give rejection threshold in angstroms{30}', .false., 30.)
        call set_param(max_dose,       'max_dose',        'num',    'Maximum dose threshold(e/A2)', 'Threshold for maximum dose and number of frames used during movie alignment(e/A2), if <=0 all frames are used{0.0}','{0.0}',.false., 0.0)
        call set_param(max_rad,        'max_rad',         'num',    'Maximum radius in A', 'Maximum radius in A {300.}', '{300.}', .true., 300.)
        call set_param(maxits,         'maxits',          'num',    'Max iterations', 'Maximum number of iterations', 'Max # iterations', .false., 100.)
        call set_param(maxnchunks,     'maxnchunks',      'num',    'Number of subsets after which 2D analysis ends', 'After this number of subsets has been classified all processing will stop(0=no end){0}','{0}',.false., 0.0)
        call set_param(mcconvention,   'mcconvention',    'str',    'Frame of reference during movie alignment', 'Frame of reference during movie alignment; simple/unblur:central; relion/motioncorr:first(simple|unblur|relion|motioncorr){simple}', '(simple|unblur|relion|motioncorr){simple}', .false., 'simple')
        call set_param(mcpatch,        'mcpatch',         'binary', 'Patch-based motion correction', 'Whether to perform Patch-based motion correction(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(mcpatch_thres,  'mcpatch_thres',   'binary', 'Use motion correction patch threshold', 'Whether to use the threshold for motion correction patch solution(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(min_rad,        'min_rad',         'num',    'Minimum radius in A', 'Minimum radius in A {50.} ', '{50.}',  .true., 50.)
        call set_param(mirr,           'mirr',            'multi',  'Perform mirroring', 'Whether to mirror and along which axis(no|x|y){no}', '(no|x|y){no}', .false., 'no')
        call set_param(ml_reg,         'ml_reg',          'binary', 'ML regularization', 'Regularization (ML-style) based on the signal power(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(ml_reg_chunk,   'ml_reg_chunk',    'binary', 'Subset ML regularization', 'Subset Regularization (ML-style) based on the signal power(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(ml_reg_pool,    'ml_reg_pool',     'binary', 'Pool ML regularization', 'Pool Regularization (ML-style) based on the signal power(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(moldiam,        'moldiam',         'num',    'Molecular diameter', 'Molecular diameter(in Angstroms)','In Angstroms',.false., 0.)
        call set_param(moldiam_max,    'moldiam_max',     'num',    'Max molecular diameter', 'Max molecular diameter(in Angstroms)','In Angstroms',.false., 300.)
        call set_param(mskdiam,        'mskdiam',         'num',    'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to remove background noise', 'mask diameter in A', .true., 0.)
        call set_param(mskfile,        'mskfile',         'file',   'Input mask file', 'Input mask file to apply to reference volume(s) before projection', 'e.g. automask.mrc from postprocess', .false., 'mskfile.mrc')
        call set_param(mul,            'mul',             'num',    'Multiplication factor', 'Multiplication factor{1.}','{1.}',.false., 1.)
        call set_param(nboxes_max,     'nboxes_max',      'num',    'Max # boxes per micrograph', 'Max # boxes per micrograph', 'Max # boxes', .false., 500.)
        call set_param(nchunks,        'nchunks',         'num',    'Number of subsets to classify simultaneously', 'Maximum number of particles subsets (chunks) to classify simultaneously', '# of subsets', .true., 2.)
        call set_param(nchunksperset,  'nchunksperset',   'num',    'Number of subsets to group', 'Number of particles subsets (chunks) to group into a independent set', '# of subsets per set', .false., 1.)
        call set_param(ncls,           'ncls',            'num',    'Number of 2D clusters', 'Number of groups to sort the particles into prior to averaging to create 2D class averages with improved SNR', '# 2D clusters', .true., 200.)
        call set_param(ncls_start,     'ncls_start',      'num',    'Number of 2D clusters per subset of particles', 'Number of class averages used in the independent 2D analysis of each subset of particles', '# 2D clusters / subset', .true., 50.)
        call set_param(neg,            'neg',             'binary', 'Invert contrast','Invert contrast(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(nparts,         'nparts',          'num',    'Number of computing nodes', 'Number of partitions for distributed memory execution. One part typically corresponds to one CPU socket in the distributed system. On a single-socket machine there may be speed benefits to dividing the jobs into a few (2-4) partitions, depending on memory capacity', 'divide job into # parts', .true., 1.0)
        call set_param(nparts_chunk,   'nparts_chunk',    'num',    'Number of computing nodes per subset', 'Number of computing nodes allocated to 2D analysis of each particles subset (chunk){1}', '# of nodes per subset{1}', .false., 1.0)
        call set_param(nparts_pool,    'nparts_pool',     'num',    'Number of computing nodes for the pooled subsets', 'Number of computing nodes allocated to 2D analysis of the pooled particles subsets', '# of nodes for the pooled subsets', .false., 2.0)
        call set_param(nptcls,         'nptcls',          'num',    'Number of particles', 'Number of particle images', '# particles', .true., 0.)
        call set_param(nptcls_per_cls, 'nptcls_per_cls',  'num',    'Number of particles per cluster', 'Initial number of particles per cluster{35}', '# initial particles per cluster{35}', .false., 500.)
        call set_param(nran,           'nran',            'num',    'Number of random samples', 'Number of entries to randomly sample', '# random samples', .false., 0.)        
        call set_param(nrestarts,      'nrestarts',       'num',    'Number of restarts', 'Number of program restarts to execute{1}', '# restarts{1}', .false., 1.0)
        call set_param(nsample,        'nsample',         'num',    'Number of particles to sample', 'Number of particles to sample each iteration', '# particles to sample', .false., 0.)
        call set_param(nsig,           'nsig',            'num',    'Number of sigmas for outlier removal', 'Number of standard deviations threshold for pixel outlier removal{6}', '# standard deviations{6}', .false., 6.)
        call set_param(nspace,         'nspace',          'num',    'Number of projection directions', 'Number of projection directions used', '# projections', .false., 2500.)
        call set_param(nstates,        'nstates',         'num',    'Number of states', 'Number of conformational/compositional states to reconstruct', '# states to reconstruct', .false., 1.0)
        call set_param(nthr,           'nthr',            'num',    'Number of threads per computing node, give 0 if unsure', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of logical threads in a socket.', '# shared-memory CPU threads', .true., 0.)
        call set_param(numlen,         'numlen',          'num',    'Length of number string', 'Length of number string', '# characters', .false., 5.0)
        call set_param(nxpatch,        'nxpatch',         'num',    '# of patches along x-axis', 'Motion correction # of patches along x-axis', '# x-patches{5}', .false., 5.)
        call set_param(nypatch,        'nypatch',         'num',    '# of patches along y-axis', 'Motion correction # of patches along y-axis', '# y-patches{5}', .false., 5.)
        call set_param(objfun,         'objfun',          'multi',  'Objective function', 'Objective function(euclid|cc|prob){euclid}', '(euclid|cc|prob){euclid}', .false., 'euclid')
        call set_param(oritab,         'oritab',          'file',   'Orientation and CTF parameter file', 'Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', '.simple|.txt parameter file', .false., 'oritab'//trim(METADATA_EXT))
        call set_param(oritab2,        'oritab2',         'file',   '2nd orientation and CTF parameter file', '2nd orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', '.simple|.txt parameter file', .false., 'oritab2'//trim(METADATA_EXT))
        call set_param(oritype,        'oritype',         'multi',  'Oritype segment in project',  'Oritype segment in project(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}', '(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}', .false., 'ptcl3D')
        call set_param(outfile,        'outfile',         'file',   'Output orientation and CTF parameter file', 'Output Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', '.simple|.txt parameter file', .false., 'outfile'//trim(METADATA_EXT))
        call set_param(outside,        'outside',         'binary', 'Extract outside stage boundaries', 'Extract boxes outside the micrograph boundaries(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(outstk,         'outstk',          'file',   'Output stack name', 'Output images stack name', 'e.g. outstk.mrc', .false., '')
        call set_param(outvol,         'outvol',          'file',   'Output volume name', 'Output volume name', 'e.g. outvol.mrc', .false., '')
        call set_param(particle_density,'particle_density','multi', 'Particle density in micrographs', 'Particle density in micrographs(low|optimal|high){optimal}', '(low|optimal|high){optimal}', .false., 'optimal')
        call set_param(pcontrast,      'pcontrast',       'multi',  'Input particle contrast', 'Input particle contrast(black|white){black}', '(black|white){black}', .false., 'black')
        call set_param(pdbout,         'pdbout',          'file',   'Output PDB volume-centered molecule', 'Output coordinates file in PDB format for the volume-centered molecule', 'e.g. output.pdb', .false., 'pdbout.pdb')
        call set_param(pgrp,           'pgrp',            'str',    'Point-group symmetry', 'Point-group symmetry of particle(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .true., 'c1')
        call set_param(pgrp_start,     'pgrp_start',      'str',    'Initital point-group symmetry', 'Initial point-group symmetry(cn|dn|t|o|i){c1}', 'point-group(cn|dn|t|o|i){c1}', .false., 'c1')
        call set_param(phaseplate,     'phaseplate',      'binary', 'Phase-plate images', 'Images obtained with Volta phase-plate(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(pick_roi,       'pick_roi',        'binary', 'Artefactual regions exclusion(new picker only)', 'Whether to exclude regions of disinterest(carbon, thick ice, new picker only)(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(picker,         'picker',          'multi',  'Which picker to use', 'Which picker to use(old|new|seg){new}', '(old|new|seg){new}', .false., 'new')
        call set_param(pickrefs,       'pickrefs',        'file',   'Stack of class-averages/reprojections for picking', 'Stack of class-averages/reprojections for picking', 'e.g. pickrefs.mrc', .false., '')
        call set_param(projfile,       'projfile',        'file',   'Project file', 'SIMPLE projectfile', 'e.g. myproject.simple', .true., '')
        call set_param(projfile_target,'projfile_target', 'file',   'Another project file', 'SIMPLE projectfile', 'e.g. myproject2.simple', .true., '')
        call set_param(projfile_merged,'projfile_merged', 'file',   'Merged output project file', 'SIMPLE projectfile', 'e.g. merged.simple', .true., '')
        call set_param(projname,       'projname',        'str',    'Project name', 'Name of project to create ./myproject/myproject.simple file for', 'e.g. to create ./myproject/myproject.simple', .true., '')
        call set_param(prune,          'prune',           'binary', 'Automated particles pruning', 'Whether to prune deselected particles(yes|no){no}', 'Automated particles pruning(yes|no){no}', .false., 'no')
        call set_param(pspecsz,        'pspecsz',         'num',    'Size of power spectrum', 'Size of power spectrum in pixels{512}', 'give # pixels{512}', .false., 512.)
        call set_param(qsys_name,      'qsys_name',       'multi',  'Queue system kind', 'Queue system kind(local|slurm|pbs|lsf)', '(local|slurm|pbs|lsf)', .false., 'local')
        call set_param(qsys_partition, 'qsys_partition',  'str',    'Name of SLURM/PBS/LSF partition', 'Name of target partition of distributed computer system (SLURM/PBS/LSF)', 'give partition name', .false., '')
        call set_param(qsys_qos,       'qsys_qos',        'str',    'Schedule priority', 'Job scheduling priority (SLURM/PBS/LSF)', 'give priority', .false., '')
        call set_param(qsys_reservation, 'qsys_reservation', 'str', 'Name of reserved partition', 'Name of reserved target partition of distributed computer system (SLURM/PBS/LSF)', 'give your part', .false., '')
        call set_param(reject_cls,     'reject_cls',      'binary', 'Whether to reject class averages', 'Whether to automatically reject 2D clusters and associated particles(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(remap_cls,      'remap_cls',       'binary', 'Whether to remap 2D clusters', 'Whether to remap the number of 2D clusters(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(remove_chunks,  'remove_chunks',   'binary', 'Whether to remove subsets', 'Whether to remove subsets after completion(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call set_param(script,         'script',          'binary', 'Generate script for shared-mem exec on cluster', 'Generate script for shared-mem exec on cluster(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(sherr,          'sherr',           'num',    'Shift error half-width', 'Uniform rotational origin shift error half-width(in pixels)', 'shift error in pixels', .false., 0.)
        call set_param(sigma_est,      'sigma_est',       'multi',  'Sigma estimation method', 'Sigma estimation method(group|global){group}', '(group|global){group}', .false., 'group')
        call set_param(smpd,           'smpd',            'num',    'Sampling distance', 'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms', .true., 1.0)
        call set_param(smpd_downscale, 'smpd_downscale',  'num',    'Sampling distance after downscale', 'Distance between neighbouring pixels in Angstroms after downscale', 'pixel size in Angstroms{1.3}', .false., 1.3)
        call set_param(smpd_target,    'smpd_target',     'num',    'Target sampling distance', 'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms', .true., 1.0)
        call set_param(star_datadir,   'star_datadir',    'file',   'STAR project data directory', 'Pathname of STAR image/data files', 'e.g. Micrographs', .false., '')
        call set_param(star_mic,       'star_mic',        'file',   'Micrographs STAR file name', 'Micrographs STAR-formatted filename', 'e.g. micrographs.star', .true., '')
        call set_param(star_model,     'star_model',      'file',   'Model STAR file name', 'Model STAR-formatted filename', 'e.g. model.star', .true., '')
        call set_param(star_ptcl,      'star_ptcl',       'file',   'Particles STAR file name', 'Particles STAR-formatted filename', 'e.g. particles.star', .true., '')
        call set_param(starfile,       'starfile',        'file',   'STAR-format file name', 'STAR-formatted filename', 'e.g. proj.star', .false., '')
        call set_param(startit,        'startit',         'num',    'First iteration', 'Index of first iteration when starting from a previous run', 'start iterations from here', .false., 1.0)
        call set_param(startype,       'startype',        'str',    'STAR-format export type', 'STAR experiment type used to define variables in export file', 'e.g. micrographs or class2d or refine3d', .false., '')
        call set_param(stepsz,         'stepsz',          'num',    'Steps size in A', 'Step size in A {10.} ', '{10.}',  .true., 10.)
        call set_param(stk,            'stk',             'file',   'Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., '')
        call set_param(stk2,           'stk2',            'file',   'Second Particle image stack', 'Particle image stack', 'xxx.mrc file with particles', .false., 'stk2.mrc')
        call set_param(stktab,         'stktab',          'file',   'List of per-micrograph particle stacks', 'List of per-micrograph particle stacks', 'stktab.txt file containing file names', .false., 'stktab.txt')
        call set_param(time_per_image, 'time_per_image',  'num',    'Time per image', 'Estimated time per image in seconds for forecasting total execution time{100}', 'in seconds{100}', .false., 100.)
        call set_param(total_dose,     'total_dose',      'num',    'Total exposure dose (e/Ang^2)', 'Total exposure dose (e/Ang^2)', 'in e/Ang^2', .false., 50.)
        call set_param(trs,            'trs',             'num',    'Maximum translational shift', 'Maximum half-width for bund-constrained search of rotational origin shifts', 'max shift per iteration in pixels{5}', .false., 5.0)
        call set_param(tseries,        'tseries',         'binary', 'Stack is time-series', 'Stack is time-series(yes|no){no}', '(yes|no){no}', .false., 'no')
        call set_param(update_frac,    'update_frac',     'num',    'Fractional update per iteration', 'Fraction of particles to update per iteration in incremental learning scheme for accelerated convergence rate(0.1-0.5){1.}', 'update this fraction per iter(0.1-0.5){1.0}', .false., 1.0)
        call set_param(user_account,   'user_account',    'str',    'User account name in SLURM/PBS/LSF', 'User account name in SLURM/PBS/LSF system', 'e.g. Account084', .false., '')
        call set_param(user_email,     'user_email',      'str',    'Your e-mail address', 'Your e-mail address', 'e.g. myname@uni.edu', .false., '')
        call set_param(user_project,   'user_project',    'str',    'User project name in SLURM/PBS/LSF', 'User project name in SLURM/PBS/LSF system', 'e.g. Project001', .false., '')
        call set_param(vol_dim,        'vol_dim',         'num',    'Simulated volume dimensions', 'Dimensions of the simulated volume in voxels', '# dimensions of the simulated volume', .false., 0.)
        call set_param(walltime,       'walltime',        'num',    'Walltime', 'Maximum execution time for job scheduling and management(23h59mins){86340}', 'in seconds(23h59mins){86340}', .false., 86340.)
        call set_param(wcrit,          'wcrit',           'multi',  'Correlation to weights conversion scheme', 'Correlation to weights conversion scheme(softmax|zscore|sum|cen|exp|inv|no){softmax}',  '(softmax|zscore|sum|cen|exp|inv|no){softmax}',  .false., 'softmax')
        call set_param(width,          'width',           'num',    'Falloff of inner mask', 'Number of cosine edge pixels of inner mask in pixels', '# pixels cosine edge{10}', .false., 10.)
        call set_param(wiener,         'wiener',          'multi',  'Wiener restoration', 'Wiener restoration, full or partial (full|partial){full}','(full|partial){full}', .false., 'full')
    end subroutine set_ui_params

    ! type-bound procedures for ui_program

    subroutine new( self, name, descr_short, descr_long, executable, n_img_ios, n_parm_ios,&
        &n_alt_ios, n_srch_ctrls, n_filt_ctrls, n_mask_ctrls, n_comp_ctrls, sp_required, gui_advanced, gui_submenu_list)
        class(ui_program),      intent(inout) :: self
        character(len=*),           intent(in)    :: name, descr_short, descr_long, executable
        integer,                    intent(in)    :: n_img_ios, n_parm_ios, n_alt_ios, n_srch_ctrls
        integer,                    intent(in)    :: n_filt_ctrls, n_mask_ctrls, n_comp_ctrls
        logical,                    intent(in)    :: sp_required
        logical,          optional, intent(in)    :: gui_advanced
        character(len=*), optional, intent(in)    :: gui_submenu_list
        call self%kill
        self%name        = trim(name)
        self%descr_short = trim(descr_short)
        self%descr_long  = trim(descr_long)
        self%executable  = trim(executable)
        if( n_img_ios    > 0 ) allocate(self%img_ios(n_img_ios)      )
        if( n_parm_ios   > 0 ) allocate(self%parm_ios(n_parm_ios)    )
        if( n_alt_ios    > 0 ) allocate(self%alt_ios(n_alt_ios)      )
        if( n_srch_ctrls > 0 ) allocate(self%srch_ctrls(n_srch_ctrls))
        if( n_filt_ctrls > 0 ) allocate(self%filt_ctrls(n_filt_ctrls))
        if( n_mask_ctrls > 0 ) allocate(self%mask_ctrls(n_mask_ctrls))
        if( n_comp_ctrls > 0 ) allocate(self%comp_ctrls(n_comp_ctrls))
        self%sp_required = sp_required
        self%exists      = .true.
        if(present(gui_advanced)) self%advanced = gui_advanced
        if(present(gui_submenu_list)) self%gui_submenu_list = gui_submenu_list
    end subroutine new

    subroutine set_input_1( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder,&
        &required, default_value, gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online)
        class(ui_program), target,   intent(inout) :: self
        character(len=*),                intent(in)    :: which
        integer,                         intent(in)    :: i
        character(len=*),                intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                         intent(in)    :: required
        real,                            intent(in)    :: default_value
        character(len=*),      optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,               optional, intent(in)    :: gui_advanced, gui_online
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; ui_program :: set_input_1')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)    :: i
                type(ui_param), intent(inout) :: arr(:)
                arr(i)%key               = trim(key)
                arr(i)%keytype           = trim(keytype)
                arr(i)%descr_short       = trim(descr_short)
                arr(i)%descr_long        = trim(descr_long)
                arr(i)%descr_placeholder = trim(descr_placeholder)
                arr(i)%required = required
                if( .not. arr(i)%required ) arr(i)%rval_default = default_value
                ! GUI options
                if( present(gui_submenu) ) then
                    arr(i)%gui_submenu = trim(gui_submenu)
                endif
                if( present (gui_exclusive_group) ) then
                    arr(i)%exclusive_group =trim(gui_exclusive_group)
                endif
                if( present(gui_active_flags) ) then
                    arr(i)%active_flags =trim(gui_active_flags)
                endif
                if( present(gui_online) ) then
                    arr(i)%online = gui_online
                endif
                if( present(gui_advanced) ) then
                    arr(i)%advanced = gui_advanced
                endif
            end subroutine set

    end subroutine set_input_1

    subroutine set_input_2( self, which, i, key, keytype, descr_short, descr_long, descr_placeholder, required,&
        &default_value, gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online)
        class(ui_program), target,   intent(inout) :: self
        character(len=*),                intent(in)    :: which
        integer,                         intent(in)    :: i
        character(len=*),                intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                         intent(in)    :: required
        character(len=*),                intent(in)    :: default_value
        character(len=*),      optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,               optional, intent(in)    :: gui_advanced, gui_online
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; ui_program :: set_input_2')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)    :: i
                type(ui_param), intent(inout) :: arr(:)
                arr(i)%key               = trim(key)
                arr(i)%keytype           = trim(keytype)
                arr(i)%descr_short       = trim(descr_short)
                arr(i)%descr_long        = trim(descr_long)
                arr(i)%descr_placeholder = trim(descr_placeholder)
                arr(i)%required = required
                if( .not. arr(i)%required ) arr(i)%cval_default = trim(default_value)
                ! GUI options
                if( present(gui_submenu) ) then
                    arr(i)%gui_submenu     = trim(gui_submenu)
                endif
                if( present (gui_exclusive_group) ) then
                    arr(i)%exclusive_group = trim(gui_exclusive_group)
                endif
                if( present(gui_active_flags) ) then
                    arr(i)%active_flags    = trim(gui_active_flags)
                endif
                if( present(gui_online) ) then
                    arr(i)%online = gui_online
                endif
                if( present(gui_advanced) ) then
                    arr(i)%advanced = gui_advanced
                endif
            end subroutine set

    end subroutine set_input_2

    subroutine set_input_3( self, which, i, param, gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online )
        class(ui_program),    target,   intent(inout) :: self
        character(len=*),                   intent(in)    :: which
        integer,                            intent(in)    :: i
        type(ui_param),           intent(in)    :: param
        character(len=*),         optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,                  optional, intent(in)    :: gui_advanced, gui_online
        select case(trim(which))
            case('img_ios')
                call set(self%img_ios, i)
            case('parm_ios')
                call set(self%parm_ios, i)
            case('alt_ios')
                call set(self%alt_ios, i)
            case('srch_ctrls')
                call set(self%srch_ctrls, i)
            case('filt_ctrls')
                call set(self%filt_ctrls, i)
            case('mask_ctrls')
                call set(self%mask_ctrls, i)
            case('comp_ctrls')
                call set(self%comp_ctrls, i)
            case DEFAULT
                THROW_HARD('which field selector: '//trim(which)//' is unsupported; ui_program :: set_input_3')
        end select

        contains

            subroutine set( arr, i )
                integer,                  intent(in)    :: i
                type(ui_param), intent(inout) :: arr(:)
                arr(i)%key               = param%key
                arr(i)%keytype           = param%keytype
                arr(i)%descr_short       = param%descr_short
                arr(i)%descr_long        = param%descr_long
                arr(i)%descr_placeholder = param%descr_placeholder
                arr(i)%required = param%required
                if( .not. arr(i)%required ) then
                    if(arr(i)%keytype%to_char() == "num") then
                        arr(i)%rval_default = param%rval_default
                    else if( param%cval_default%is_allocated() ) then
                        arr(i)%cval_default = param%cval_default
                    endif
                endif
                ! GUI options
                if( present(gui_submenu) ) then
                    arr(i)%gui_submenu = trim(gui_submenu)
                endif
                if( present (gui_exclusive_group) ) then
                    arr(i)%exclusive_group = trim(gui_exclusive_group)
                endif
                if( present(gui_active_flags) ) then
                    arr(i)%active_flags = trim(gui_active_flags)
                endif
                if( present(gui_online) ) then
                    arr(i)%online = gui_online
                endif
                if( present(gui_advanced) ) then
                    arr(i)%advanced = gui_advanced
                endif
            end subroutine set

    end subroutine set_input_3

    subroutine print_ui( self )
        class(ui_program), intent(in) :: self
        type(chash) :: ch
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name%to_char())
        call ch%push('descr_short', self%descr_short%to_char())
        call ch%push('descr_long',  self%descr_long%to_char())
        call ch%push('executable',  self%executable%to_char())
        call ch%print_key_val_pairs(logfhandle)
        call ch%kill
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash_local(self%img_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash_local(self%parm_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash_local(self%alt_ios)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash_local(self%srch_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash_local(self%filt_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash_local(self%mask_ctrls)
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash_local(self%comp_ctrls)

        contains

            subroutine print_param_hash_local( arr )
                type(ui_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        write(logfhandle,'(a,1x,i3)') '>>> PARAMETER #', i
                        call ch%new(6)
                        call ch%push('key',               arr(i)%key%to_char())
                        call ch%push('keytype',           arr(i)%keytype%to_char())
                        call ch%push('descr_short',       arr(i)%descr_short%to_char())
                        call ch%push('descr_long',        arr(i)%descr_long%to_char())
                        call ch%push('descr_placeholder', arr(i)%descr_placeholder%to_char())
                        if( arr(i)%required )then
                            call ch%push('required', 'T')
                        else
                            call ch%push('required', 'F')
                        endif
                        call ch%print_key_val_pairs(logfhandle)
                        call ch%kill
                    end do
                endif
            end subroutine print_param_hash_local

    end subroutine print_ui

    ! supporting the print_cmdline routines (above)
    subroutine print_param_hash( arr )
        type(ui_param), allocatable, intent(in) :: arr(:)
        character(len=KEYLEN),    allocatable :: sorted_keys(:), rearranged_keys(:)
        logical,                  allocatable :: required(:)
        integer,                  allocatable :: inds(:)
        type(chash) :: ch
        integer     :: i, nparams, nreq, iopt
        if( allocated(arr) )then
            nparams = size(arr)
            call ch%new(nparams)
            allocate(sorted_keys(nparams), rearranged_keys(nparams), required(nparams))
            do i=1,nparams
                call ch%push(arr(i)%key%to_char(), arr(i)%descr_short%to_char()//'; '//arr(i)%descr_placeholder%to_char())
                sorted_keys(i) = arr(i)%key%to_char()
                required(i)    = arr(i)%required
            end do
            call lex_sort(sorted_keys, inds=inds)
            required = required(inds)
            if( any(required) )then
                ! fish out the required ones
                nreq = 0
                do i=1,nparams
                    if( required(i) )then
                        nreq = nreq + 1
                        rearranged_keys(nreq) = sorted_keys(i)
                    endif
                enddo
                ! fish out the optional ones
                iopt = nreq
                do i=1,nparams
                    if( .not. required(i) )then
                        iopt = iopt + 1
                        rearranged_keys(iopt) = sorted_keys(i)
                    endif
                end do
                ! replace string array
                sorted_keys = rearranged_keys
                ! modify logical mask
                required(:nreq)     = .true.
                required(nreq + 1:) = .false.
            endif
            call ch%print_key_val_pairs(logfhandle, sorted_keys, mask=required)
            call ch%kill
            deallocate(sorted_keys, required)
        endif
    end subroutine print_param_hash

    subroutine print_cmdline( self )
        class(ui_program), intent(in) :: self
        write(logfhandle,'(a)') format_str('USAGE', C_UNDERLINED)
        write(logfhandle,'(a)') format_str('bash-3.2$ simple_exec prg=' //self%name%to_char() // ' key1=val1 key2=val2 ...', C_ITALIC)
        write(logfhandle,'(a)') 'Required input parameters in ' // format_str('bold', C_BOLD) // ' (ensure terminal support)'
        if( allocated(self%img_ios) )    write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
        call print_param_hash(self%img_ios)
        if( allocated(self%parm_ios) )   write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
        call print_param_hash(self%parm_ios)
        if( allocated(self%alt_ios) )    write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
        call print_param_hash(self%alt_ios)
        if( allocated(self%srch_ctrls) ) write(logfhandle,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%srch_ctrls)
        if( allocated(self%filt_ctrls) ) write(logfhandle,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
        call print_param_hash(self%filt_ctrls)
        if( allocated(self%mask_ctrls) ) write(logfhandle,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
        call print_param_hash(self%mask_ctrls)
        if( allocated(self%comp_ctrls) ) write(logfhandle,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
        call print_param_hash(self%comp_ctrls)
    end subroutine print_cmdline

    subroutine print_prg_descr_long( self )
        class(ui_program), intent(in) :: self
        write(logfhandle,'(a)') self%descr_long%to_char()
    end subroutine print_prg_descr_long

    subroutine write2json( self )
        use json_module
        class(ui_program), intent(in) :: self
        type(json_core)           :: json
        type(json_value), pointer :: program_entry, program
        ! JSON init
        call json%initialize()
        call json%create_object(program_entry,'')
        call json%create_object(program, self%name%to_char())
        call json%add(program_entry, program)
        ! program section
        call json%add(program, 'name',        self%name%to_char())
        call json%add(program, 'descr_short', self%descr_short%to_char())
        call json%add(program, 'descr_long',  self%descr_long%to_char())
        call json%add(program, 'executable',  self%executable%to_char())
        ! all sections
        call create_section( 'image input/output',     self%img_ios )
        call create_section( 'parameter input/output', self%parm_ios )
        call create_section( 'alternative inputs',     self%alt_ios )
        call create_section( 'search controls',        self%srch_ctrls )
        call create_section( 'filter controls',        self%filt_ctrls )
        call create_section( 'mask controls',          self%mask_ctrls )
        call create_section( 'computer controls',      self%comp_ctrls )
        ! write & clean
        call json%print(program_entry, self%name%to_char()//'.json')
        if( json%failed() )then
            THROW_HARD('json input/output error for program: '//self%name%to_char())
        endif
        call json%destroy(program_entry)

        contains

            subroutine create_section( name, arr )
                character(len=*),          intent(in) :: name
                type(ui_param), allocatable, intent(in) :: arr(:)
                type(json_value), pointer :: entry, section
                character(len=STDLEN)     :: options_str, before
                character(len=KEYLEN)     :: args(8)
                integer                   :: i, j, sz, nargs
                logical :: found, param_is_multi, param_is_binary, exception
                call json%create_array(section, trim(name))
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        call json%create_object(entry,            arr(i)%key%to_char())
                        call json%add(entry, 'key',               arr(i)%key%to_char())
                        call json%add(entry, 'keytype',           arr(i)%keytype%to_char())
                        call json%add(entry, 'descr_short',       arr(i)%descr_short%to_char())
                        call json%add(entry, 'descr_long',        arr(i)%descr_long%to_char())
                        call json%add(entry, 'descr_placeholder', arr(i)%descr_placeholder%to_char())
                        call json%add(entry, 'required',          arr(i)%required)
                        param_is_multi  = arr(i)%keytype%to_char().eq.'multi'
                        param_is_binary = arr(i)%keytype%to_char().eq.'binary'
                        if( param_is_multi .or. param_is_binary )then
                            options_str = arr(i)%descr_placeholder%to_char()
                            call split( options_str, '(', before )
                            call split( options_str, ')', before )
                            call parsestr(before, '|', args, nargs)
                            exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 2)
                            if( exception )then
                                write(logfhandle,*)'Poorly formatted options string for entry ', arr(i)%key%to_char()
                                write(logfhandle,*) arr(i)%descr_placeholder%to_char()
                                stop
                            endif
                            call json%add(entry, 'options', args(1:nargs))
                            do j = 1, nargs
                                call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                            enddo
                        endif
                        call json%add(section, entry)
                    enddo
                endif
                call json%add(program_entry, section)
            end subroutine create_section

    end subroutine write2json

     function get_name( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%name
    end function get_name

    function get_executable( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%executable
    end function get_executable

    integer function get_nrequired_keys( self )
        class(ui_program), intent(in) :: self
        get_nrequired_keys = nreq_counter(self%img_ios) + nreq_counter(self%parm_ios) +&
        &nreq_counter(self%srch_ctrls) + nreq_counter(self%filt_ctrls) +&
        &nreq_counter(self%mask_ctrls) + nreq_counter(self%comp_ctrls)
        if( get_nrequired_keys == 0 .and. allocated(self%alt_ios) ) get_nrequired_keys = 1

        contains

            function nreq_counter( arr ) result( nreq )
                type(ui_param), allocatable, intent(in) :: arr(:)
                integer :: nreq, i
                nreq = 0
                if( allocated(arr) )then
                    do i=1,size(arr)
                        if( arr(i)%required ) nreq = nreq + 1
                    end do
                endif
            end function nreq_counter

    end function get_nrequired_keys

    function get_required_keys( self ) result( keys )
        class(ui_program), intent(in) :: self
        type(string), allocatable :: keys(:)
        integer :: nreq, ireq
        ! count # required
        nreq = self%get_nrequired_keys()
        ! extract keys
        if( nreq > 0 )then
            allocate(keys(nreq))
            ireq = 0
            call key_extractor(self%img_ios)
            call key_extractor(self%parm_ios)
            call key_extractor(self%alt_ios)
            call key_extractor(self%srch_ctrls)
            call key_extractor(self%filt_ctrls)
            call key_extractor(self%mask_ctrls)
            call key_extractor(self%comp_ctrls)
        endif

        contains

            subroutine key_extractor( arr )
                type(ui_param), allocatable, intent(in) :: arr(:)
                integer :: i
                if( allocated(arr) )then
                    do i=1,size(arr)
                        if( arr(i)%required )then
                            ireq = ireq + 1
                            keys(ireq) = arr(i)%key
                        endif
                    end do
                endif
            end subroutine key_extractor

    end function get_required_keys

    logical function requires_sp_project( self )
        class(ui_program), intent(in) :: self
        requires_sp_project = self%sp_required
    end function requires_sp_project

    subroutine kill( self )
        class(ui_program), intent(inout) :: self
        integer :: i, sz
        if( self%exists )then
            call self%name%kill
            call self%descr_short%kill
            call self%descr_long%kill
            call self%executable%kill
            call dealloc_field(self%img_ios)
            call dealloc_field(self%parm_ios)
            call dealloc_field(self%alt_ios)
            call dealloc_field(self%srch_ctrls)
            call dealloc_field(self%filt_ctrls)
            call dealloc_field(self%mask_ctrls)
            call dealloc_field(self%comp_ctrls)
            self%exists = .false.
        endif

        contains

            subroutine dealloc_field( arr )
                type(ui_param), allocatable, intent(inout) :: arr(:)
                if( allocated(arr) )then
                    sz = size(arr)
                    do i=1,sz
                        call arr(i)%key%kill
                        call arr(i)%keytype%kill
                        call arr(i)%descr_short%kill
                        call arr(i)%descr_long%kill
                        call arr(i)%descr_placeholder%kill
                        call arr(i)%gui_submenu%kill
                        call arr(i)%active_flags%kill
                        call arr(i)%exclusive_group%kill
                        call arr(i)%cval_default%kill
                    end do
                    deallocate(arr)
                endif
            end subroutine dealloc_field

    end subroutine kill

end module simple_ui_program
