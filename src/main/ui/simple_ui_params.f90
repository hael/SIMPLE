module simple_ui_params
use simple_core_module_api
implicit none
public

! common input parameter type
type simple_ui_param
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
end type simple_ui_param

interface set_param
   module procedure set_param_1
   module procedure set_param_2
end interface set_param

! declare common params here, with name same as flag
type(simple_ui_param) :: algorithm
type(simple_ui_param) :: angerr
type(simple_ui_param) :: astigthreshold
type(simple_ui_param) :: astigtol
type(simple_ui_param) :: automsk
type(simple_ui_param) :: autosample
type(simple_ui_param) :: backgr_subtr
type(simple_ui_param) :: bfac
type(simple_ui_param) :: box
type(simple_ui_param) :: box_extract
type(simple_ui_param) :: cc_iters
type(simple_ui_param) :: center_pdb
type(simple_ui_param) :: clip
type(simple_ui_param) :: cls_init
type(simple_ui_param) :: clust_crit
type(simple_ui_param) :: cn
type(simple_ui_param) :: cn_max
type(simple_ui_param) :: cn_min
type(simple_ui_param) :: combine_eo
type(simple_ui_param) :: cs
type(simple_ui_param) :: ctf
type(simple_ui_param) :: ctf_yes
type(simple_ui_param) :: ctfpatch
type(simple_ui_param) :: ctfresthreshold
type(simple_ui_param) :: deftab
type(simple_ui_param) :: dferr
type(simple_ui_param) :: dfmax
type(simple_ui_param) :: dfmin
type(simple_ui_param) :: dir_movies
type(simple_ui_param) :: e1, e2, e3
type(simple_ui_param) :: eer_fraction
type(simple_ui_param) :: element
type(simple_ui_param) :: envfsc
type(simple_ui_param) :: eo
type(simple_ui_param) :: flipgain
type(simple_ui_param) :: frac
type(simple_ui_param) :: fraca
type(simple_ui_param) :: fraction_dose_target
type(simple_ui_param) :: frcs
type(simple_ui_param) :: gainref
type(simple_ui_param) :: graphene_filt
type(simple_ui_param) :: hp
type(simple_ui_param) :: icefracthreshold
type(simple_ui_param) :: icm
type(simple_ui_param) :: job_memory_per_task
type(simple_ui_param) :: kv
type(simple_ui_param) :: lp
type(simple_ui_param) :: lp_backgr
type(simple_ui_param) :: lp_pick
type(simple_ui_param) :: lplim_crit
type(simple_ui_param) :: lpstart_nonuni
type(simple_ui_param) :: lpthres
type(simple_ui_param) :: max_dose
type(simple_ui_param) :: max_rad
type(simple_ui_param) :: maxits
type(simple_ui_param) :: maxnchunks
type(simple_ui_param) :: mcconvention
type(simple_ui_param) :: mcpatch
type(simple_ui_param) :: mcpatch_thres
type(simple_ui_param) :: min_rad
type(simple_ui_param) :: mirr
type(simple_ui_param) :: ml_reg
type(simple_ui_param) :: ml_reg_chunk
type(simple_ui_param) :: ml_reg_pool
type(simple_ui_param) :: moldiam
type(simple_ui_param) :: moldiam_max
type(simple_ui_param) :: mskdiam
type(simple_ui_param) :: mskfile
type(simple_ui_param) :: mul
type(simple_ui_param) :: nboxes_max
type(simple_ui_param) :: nchunks
type(simple_ui_param) :: nchunksperset
type(simple_ui_param) :: ncls
type(simple_ui_param) :: ncls_start
type(simple_ui_param) :: neg
type(simple_ui_param) :: niceprocid
type(simple_ui_param) :: niceserver
type(simple_ui_param) :: nparts
type(simple_ui_param) :: nparts_chunk
type(simple_ui_param) :: nparts_pool
type(simple_ui_param) :: nptcls
type(simple_ui_param) :: nptcls_per_cls
type(simple_ui_param) :: nran
type(simple_ui_param) :: nrestarts
type(simple_ui_param) :: nsample
type(simple_ui_param) :: nsearch
type(simple_ui_param) :: nsig
type(simple_ui_param) :: nspace
type(simple_ui_param) :: nstates
type(simple_ui_param) :: nthr
type(simple_ui_param) :: numlen
type(simple_ui_param) :: nxpatch
type(simple_ui_param) :: nypatch
type(simple_ui_param) :: objfun
type(simple_ui_param) :: oritab
type(simple_ui_param) :: oritab2
type(simple_ui_param) :: oritype
type(simple_ui_param) :: outdir
type(simple_ui_param) :: outfile
type(simple_ui_param) :: outside
type(simple_ui_param) :: outstk
type(simple_ui_param) :: outvol
type(simple_ui_param) :: particle_density
type(simple_ui_param) :: pcontrast
type(simple_ui_param) :: pdbout
type(simple_ui_param) :: pgrp
type(simple_ui_param) :: pgrp_start
type(simple_ui_param) :: phaseplate
type(simple_ui_param) :: pick_roi
type(simple_ui_param) :: picker
type(simple_ui_param) :: pickrefs
type(simple_ui_param) :: projfile
type(simple_ui_param) :: projfile_merged
type(simple_ui_param) :: projfile_target
type(simple_ui_param) :: projname
type(simple_ui_param) :: prune
type(simple_ui_param) :: pspecsz
type(simple_ui_param) :: qsys_name
type(simple_ui_param) :: qsys_partition
type(simple_ui_param) :: qsys_qos
type(simple_ui_param) :: qsys_reservation
type(simple_ui_param) :: reject_cls
type(simple_ui_param) :: remap_cls
type(simple_ui_param) :: remove_chunks
type(simple_ui_param) :: script
type(simple_ui_param) :: sherr
type(simple_ui_param) :: sigma
type(simple_ui_param) :: sigma_est
type(simple_ui_param) :: smpd
type(simple_ui_param) :: smpd_downscale
type(simple_ui_param) :: smpd_target
type(simple_ui_param) :: star_datadir
type(simple_ui_param) :: star_mic
type(simple_ui_param) :: star_model
type(simple_ui_param) :: star_ptcl
type(simple_ui_param) :: starfile
type(simple_ui_param) :: startit
type(simple_ui_param) :: startype
type(simple_ui_param) :: stepsz
type(simple_ui_param) :: stk
type(simple_ui_param) :: stk2
type(simple_ui_param) :: stktab
type(simple_ui_param) :: time_per_image
type(simple_ui_param) :: total_dose
type(simple_ui_param) :: trs
type(simple_ui_param) :: tseries
type(simple_ui_param) :: update_frac
type(simple_ui_param) :: user_account
type(simple_ui_param) :: user_email
type(simple_ui_param) :: user_project
type(simple_ui_param) :: vol_dim
type(simple_ui_param) :: walltime
type(simple_ui_param) :: wcrit
type(simple_ui_param) :: width
type(simple_ui_param) :: wiener
type(simple_ui_param) :: winsz

contains

    subroutine set_param_1( self, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value )
        class(simple_ui_param), intent(inout) :: self
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
        class(simple_ui_param), intent(inout) :: self
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

end module simple_ui_params
