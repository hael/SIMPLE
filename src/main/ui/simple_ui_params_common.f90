!@descr: module defining the common parameters for all simple_ui_program interfaces
module simple_ui_params_common
use simple_core_module_api
use simple_ui_param, only: ui_param
implicit none
! declare common params here, with name same as flag
type(ui_param) :: algorithm
type(ui_param) :: angerr
type(ui_param) :: astigthreshold
type(ui_param) :: astigtol
type(ui_param) :: automsk
type(ui_param) :: autosample
type(ui_param) :: backgr_subtr
type(ui_param) :: bfac
type(ui_param) :: box
type(ui_param) :: box_extract
type(ui_param) :: center_pdb
type(ui_param) :: clip
type(ui_param) :: cls_init
type(ui_param) :: clust_crit
type(ui_param) :: cn
type(ui_param) :: cn_max
type(ui_param) :: cn_min
type(ui_param) :: combine_eo
type(ui_param) :: cs
type(ui_param) :: ctf
type(ui_param) :: ctf_yes
type(ui_param) :: ctfpatch
type(ui_param) :: ctfresthreshold
type(ui_param) :: deftab
type(ui_param) :: dferr
type(ui_param) :: dfmax
type(ui_param) :: dfmin
type(ui_param) :: dir_movies
type(ui_param) :: e1, e2, e3
type(ui_param) :: eer_fraction
type(ui_param) :: element
type(ui_param) :: envfsc
type(ui_param) :: eo
type(ui_param) :: flipgain
type(ui_param) :: frac
type(ui_param) :: fraca
type(ui_param) :: fraction_dose_target
type(ui_param) :: frcs
type(ui_param) :: gainref
type(ui_param) :: graphene_filt
type(ui_param) :: hp
type(ui_param) :: icefracthreshold
type(ui_param) :: icm
type(ui_param) :: job_memory_per_task
type(ui_param) :: kv
type(ui_param) :: lp
type(ui_param) :: lp_backgr
type(ui_param) :: lp_pick
type(ui_param) :: lp_track
type(ui_param) :: lplim_crit
type(ui_param) :: lpstart_nonuni
type(ui_param) :: lpthres
type(ui_param) :: max_dose
type(ui_param) :: max_rad
type(ui_param) :: maxits
type(ui_param) :: maxnchunks
type(ui_param) :: mcconvention
type(ui_param) :: mcpatch
type(ui_param) :: mcpatch_thres
type(ui_param) :: min_rad
type(ui_param) :: mirr
type(ui_param) :: ml_reg
type(ui_param) :: ml_reg_chunk
type(ui_param) :: ml_reg_pool
type(ui_param) :: moldiam
type(ui_param) :: moldiam_max
type(ui_param) :: mskdiam
type(ui_param) :: mskfile
type(ui_param) :: mul
type(ui_param) :: nboxes_max
type(ui_param) :: nchunks
type(ui_param) :: nchunksperset
type(ui_param) :: ncls
type(ui_param) :: ncls_start
type(ui_param) :: neg
type(ui_param) :: niceprocid
type(ui_param) :: niceserver
type(ui_param) :: nparts
type(ui_param) :: nparts_chunk
type(ui_param) :: nparts_pool
type(ui_param) :: nptcls
type(ui_param) :: nptcls_per_cls
type(ui_param) :: nptcls_per_cls_cleanup2D
type(ui_param) :: nran
type(ui_param) :: nrestarts
type(ui_param) :: nsample
type(ui_param) :: nsearch
type(ui_param) :: nsig
type(ui_param) :: nspace
type(ui_param) :: nstates
type(ui_param) :: nthr
type(ui_param) :: numlen
type(ui_param) :: nxpatch
type(ui_param) :: nypatch
type(ui_param) :: objfun
type(ui_param) :: oritab
type(ui_param) :: oritab2
type(ui_param) :: oritype
type(ui_param) :: outdir
type(ui_param) :: outfile
type(ui_param) :: outside
type(ui_param) :: outstk
type(ui_param) :: outvol
type(ui_param) :: particle_density
type(ui_param) :: pcontrast
type(ui_param) :: pdbout
type(ui_param) :: pgrp
type(ui_param) :: pgrp_start
type(ui_param) :: phaseplate
type(ui_param) :: pick_roi
type(ui_param) :: picker
type(ui_param) :: pickrefs
type(ui_param) :: projfile
type(ui_param) :: projfile_merged
type(ui_param) :: projfile_ref
type(ui_param) :: projfile_target
type(ui_param) :: projname
type(ui_param) :: prune
type(ui_param) :: pspecsz
type(ui_param) :: qsys_name
type(ui_param) :: qsys_partition
type(ui_param) :: qsys_qos
type(ui_param) :: qsys_reservation
type(ui_param) :: remap_cls
type(ui_param) :: remove_chunks
type(ui_param) :: script
type(ui_param) :: sherr
type(ui_param) :: sigma
type(ui_param) :: sigma_est
type(ui_param) :: smpd
type(ui_param) :: smpd_downscale
type(ui_param) :: smpd_target
type(ui_param) :: star_datadir
type(ui_param) :: star_mic
type(ui_param) :: star_model
type(ui_param) :: star_ptcl
type(ui_param) :: starfile
type(ui_param) :: startit
type(ui_param) :: startype
type(ui_param) :: stepsz
type(ui_param) :: stk
type(ui_param) :: stk_backgr
type(ui_param) :: stk_traj
type(ui_param) :: stk2
type(ui_param) :: stktab
type(ui_param) :: time_per_image
type(ui_param) :: total_dose
type(ui_param) :: trs
type(ui_param) :: trs_mc
type(ui_param) :: tseries
type(ui_param) :: update_frac
type(ui_param) :: user_account
type(ui_param) :: user_email
type(ui_param) :: user_project
type(ui_param) :: vol_dim
type(ui_param) :: walltime
type(ui_param) :: wcrit
type(ui_param) :: width
type(ui_param) :: wiener
type(ui_param) :: winsz

contains

subroutine set_ui_params
    call algorithm%set_param(      'algorithm',       'multi',  'Algorithm for motion correction', &
                                   'Algorithm for motion correction(iso|patch|patch_refine){patch}', &
                                   '(iso|patch|patch_refine){patch}', .false., 'patch')

    call angerr%set_param(         'angerr',          'num',    'Rotation angle error half-width', &
                                   'Uniform rotation angle shift error half-width(in degrees)', &
                                   'rotation error in degrees', .false., 0.)

    call astigthreshold%set_param( 'astigthreshold',  'num',    'Astigmatism rejection threshold', &
                                   'Micrographs with astigmatism (%) above the threshold will be ignored from further processing{10.0}', &
                                   'Astigmatism threshold{10.0}', .false., 10.0)

    call astigtol%set_param(       'astigtol',        'num',    'Expected astigmatism', &
                                   'expected (tolerated) astigmatism(in microns){0.05}', &
                                   'in microns{0.05}', .false., 0.05)

    call automsk%set_param(        'automsk',         'multi',  'Perform envelope masking', &
                                   'Whether to generate/apply an envelope mask(yes|tight|no){no}', &
                                   '(yes|tight|no){no}', .false., 'no')

    call autosample%set_param(     'autosample',      'binary', 'Automated particles sampling scheme', &
                                   'Use automated sampling scheme to select particles subsets(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call backgr_subtr%set_param(   'backgr_subtr',    'binary', 'Perform micrograph background subtraction(new picker only)', &
                                   'Perform micrograph background subtraction before picking/extraction(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call bfac%set_param(           'bfac',            'num',    'B-factor for sharpening', &
                                   'B-factor for sharpening in Angstroms^2', &
                                   'B-factor in Angstroms^2', .false., 200.)

    call box%set_param(            'box',             'num',    'Particle box size', &
                                   'Particle box size(in pixels)', &
                                   '# pixels of box', .true., 0.)

    call box_extract%set_param(    'box_extract',     'num',    'Extracted particle image size', &
                                   'Extracted particle image size(in pixels)', &
                                   'Extracted particle image size', .false., 0.)

    call center_pdb%set_param(     'center_pdb',      'binary', 'Whether to move the PDB atomic center to the center of the box', &
                                   'Whether to move the PDB atomic center to the center of the box (yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call clip%set_param(           'clip',            'num',    'Clipped box size', &
                                   'Target box size for clipping in pixels', &
                                   'in pixels', .false., 0.)

    call cls_init%set_param(       'cls_init',        'multi',  'Scheme for initial class generation', &
                                   'Initiate 2D analysis from raw images|random classes|noise images(ptcl|randcls|rand){ptcl}', &
                                   '(ptcl|randcls|rand){ptcl}', .false., 'rand')

    call clust_crit%set_param(     'clust_crit',      'multi',  'Clustering criterion', &
                                   'Clustering criterion(sig|sig_clust|cc|res|hybrid){hybrid}', &
                                   '(sig|sig_clust|cc|res|hybrid){hybrid}', .false., 'hybrid')

    call cn%set_param(             'cn',              'num',    'Fixed std coordination number', &
                                   'Minimum std cn to consider for dipole calc ', &
                                   '8', .false., 8.)

    call cn_max%set_param(         'cn_max',          'num',    'Maximum std coordination number', &
                                   'Maximum std cn to consider ', &
                                   '12', .false., 12.)

    call cn_min%set_param(         'cn_min',          'num',    'Minimum std coordination number', &
                                   'Minimum std cn to consider ', &
                                   '4', .false., 4.)

    call combine_eo%set_param(     'combine_eo',      'binary', 'Whether e/o references are combined for final alignment(yes|no){no}', &
                                   'whether e/o references are combined for final alignment(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call cs%set_param(             'cs',              'num',    'Spherical aberration', &
                                   'Spherical aberration constant(in mm){2.7}', &
                                   'in mm{2.7}', .false., 2.7)

    call ctf%set_param(            'ctf',             'multi',  'CTF status', &
                                   'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){no}', &
                                   '(yes|no|flip){no}', .true., 'no')

    call ctf_yes%set_param(        'ctf',             'multi',  'CTF status', &
                                   'Contrast Transfer Function status; flip indicates that images have been phase-flipped prior(yes|no|flip){yes}', &
                                   '(yes|no|flip){yes}', .false., 'yes')

    call ctfpatch%set_param(       'ctfpatch',        'binary', 'Patch CTF estimation', &
                                   'Whether to perform patch CTF estimation(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call ctfresthreshold%set_param('ctfresthreshold', 'num',    'CTF Resolution rejection threshold', &
                                   'Micrographs with a CTF resolution above the threshold (in Angs) will be ignored from further processing{6.}', &
                                   'CTF resolution threshold(in Angstroms){6.}', .false., 6.0)

    call deftab%set_param(         'deftab',          'file',   'CTF parameter file', &
                                   'CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format with dfx, dfy and angast values', &
                                   '.simple|.txt parameter file', .false., 'deftab'//trim(METADATA_EXT))

    call dferr%set_param(          'dferr',           'num',    'Underfocus error half-width', &
                                   'Uniform underfoucs error half-width(in microns)', &
                                   'defocus error in microns', .false., 1.)

    call dfmax%set_param(          'dfmax',           'num',    'Expected maximum defocus', &
                                   'Expected maximum defocus in microns{5.0}', &
                                   'in microns{5.0}', .false., DFMAX_DEFAULT)

    call dfmin%set_param(          'dfmin',           'num',    'Expected minimum defocus', &
                                   'Expected minimum defocus in microns{0.2}', &
                                   'in microns{0.2}', .false., DFMIN_DEFAULT)

    call dir_movies%set_param(     'dir_movies',      'dir',    'Input movies directory', &
                                   'Where the movies to process are located or will squentially appear', &
                                   'e.g. /cryodata/', .true., 'preprocess/')

    call e1%set_param(             'e1',              'num',    'Rotation along Phi', &
                                   'Phi Euler angle', &
                                   'in degrees', .false., 0.)

    call e2%set_param(             'e2',              'num',    'Rotation along Theta', &
                                   'Theat Euler angle', &
                                   'in degrees', .false., 0.)

    call e3%set_param(             'e3',              'num',    'Rotation along Psi', &
                                   'Psi Euler angle', &
                                   'in degrees', .false., 0.)

    call eer_fraction%set_param(   'eer_fraction',    'num',    '# of EER frames to fraction together', &
                                   'Number of raw EER frames in a movie fraction', &
                                   '# EER frames{20}', .false., 20.)

    call element%set_param(        'element',         'str',    'Atom element name: Au, Pt etc.', &
                                   'Atom element name: Au, Pt etc.', &
                                   'atom composition e.g. Pt', .true., '  ')

    call envfsc%set_param(         'envfsc',          'binary', 'Envelope mask e/o maps for FSC', &
                                   'Envelope mask even/odd pairs prior to FSC calculation(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call eo%set_param(             'eo',              'binary', 'Gold-standard FSC for filtering and resolution estimation', &
                                   'Gold-standard FSC for filtering and resolution estimation(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call flipgain%set_param(       'flipgain',        'multi',  'Flip the gain reference', &
                                   'Flip the gain reference along the provided axis(no|x|y|xy|yx){no}', &
                                   '(no|x|y|xy|yx){no}', .false., 'no')

    call frac%set_param(           'frac',            'num',    'Fraction of particles to include', &
                                   'Fraction of particles to include based on spectral score (median of FRC between reference and particle)', &
                                   'fraction of particles(0.1-0.9){1.0}', .false., 1.0)

    call fraca%set_param(          'fraca',           'num',    'Amplitude contrast fraction', &
                                   'Fraction of amplitude contrast used for fitting CTF{0.1}', &
                                   'fraction{0.1}', .false., 0.1)

    call fraction_dose_target%set_param('fraction_dose_target','num','EER fraction dose target (e/Ang^2)', &
                                        'EER fraction dose target, used to determine how many EER frames are included in each movie fraction(e/Ang^2)', &
                                        'in e/Ang^2', .false., 1.)

    call frcs%set_param(           'frcs',            'str',    'Projection FRCs file', &
                                   'Projection FRCs file', &
                                   'e.g. frcs.bin', .false., '')

    call gainref%set_param(        'gainref',         'file',   'Gain reference', &
                                   'Gain reference image', &
                                   'input image e.g. gainref.mrc', .false., '')

    call graphene_filt%set_param(  'graphene_filt',   'binary', 'Omit graphene bands from corr calc', &
                                   'Omit graphene bands from corr calc(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call hp%set_param(             'hp',              'num',    'High-pass limit', &
                                   'High-pass resolution limit', &
                                   'high-pass limit in Angstroms', .false., 100.)

    call icefracthreshold%set_param('icefracthreshold','num',   'Ice Fraction rejection threshold', &
                                    'Micrographs with an ice ring/1st pspec maxima fraction above the threshold will be ignored from further processing{1.0}', &
                                    'Ice fraction threshold{1.0}', .false., 1.0)

    call icm%set_param(            'icm',             'binary', 'Whether to perform ICM filtering of reference(s)', &
                                   'Whether to perform ICM filtering of reference(s)(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call job_memory_per_task%set_param('job_memory_per_task','str', 'Memory per computing node', &
                                       'Memory in MB per part/computing node in distributed execution{16000}', &
                                       'MB per part{16000}', .false., '16000')

    call kv%set_param(             'kv',              'num',    'Acceleration voltage', &
                                   'Acceleration voltage in kV{300}', &
                                   'in kV{300}', .false., 300.)

    call lp%set_param(             'lp',              'num',    'Low-pass limit', &
                                   'Low-pass resolution limit', &
                                   'low-pass limit in Angstroms', .false., 20.)

    call lp_backgr%set_param(      'lp_backgr',       'num',    'Background low-pass resolution', &
                                   'Low-pass resolution for solvent blurring', &
                                   'low-pass limit in Angstroms', .false., 20.)

    call lp_pick%set_param(        'lp_pick',         'num',    'Low-pass limit for picking', &
                                   'Low-pass limit for picking in Angstroms{20}', &
                                   'in Angstroms{20}', .false., 20.)

    call lp_track%set_param(       'lp_pick',         'num',    'Low-pass limit in Angs', &
                                   'Low-pass limit in Angs{2.3}', &
                                   'lp in Angs{2.3}', .false., 2.3)


    call lplim_crit%set_param(     'lplim_crit',      'num',    'Low-pass limit FSC criterion', &
                                   'FSC criterion for determining the low-pass limit(0.143-0.5){0.143}', &
                                   'low-pass FSC criterion(0.143-0.5){0.143}', .false., 0.143)

    call lpthres%set_param(        'lpthres',         'num',    'Resolution rejection threshold', &
                                   'Classes with lower resolution are iteratively rejected in Angstroms{30}', &
                                   'give rejection threshold in angstroms{30}', .false., 30.)

    call max_dose%set_param(       'max_dose',        'num',    'Maximum dose threshold(e/A2)', &
                                   'Threshold for maximum dose and number of frames used during movie alignment(e/A2), if <=0 all frames are used{0.0}', &
                                   '{0.0}', .false., 0.0)

    call max_rad%set_param(        'max_rad',         'num',    'Maximum radius in A', &
                                   'Maximum radius in A {300.}', &
                                   '{300.}', .true., 300.)

    call maxits%set_param(         'maxits',          'num',    'Max iterations', &
                                   'Maximum number of iterations', &
                                   'Max # iterations', .false., 100.)

    call maxnchunks%set_param(     'maxnchunks',      'num',    'Number of subsets after which 2D analysis ends', &
                                   'After this number of subsets has been classified all processing will stop(0=no end){0}', &
                                   '{0}', .false., 0.0)

    call mcconvention%set_param(   'mcconvention',    'str',    'Frame of reference during movie alignment', &
                                   'Frame of reference during movie alignment; simple/unblur:central; relion/motioncorr:first(simple|unblur|relion|motioncorr){simple}', &
                                   '(simple|unblur|relion|motioncorr){simple}', .false., 'simple')

    call mcpatch%set_param(        'mcpatch',         'binary', 'Patch-based motion correction', &
                                   'Whether to perform Patch-based motion correction(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call mcpatch_thres%set_param(  'mcpatch_thres',   'binary', 'Use motion correction patch threshold', &
                                   'Whether to use the threshold for motion correction patch solution(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call min_rad%set_param(        'min_rad',         'num',    'Minimum radius in A', &
                                   'Minimum radius in A {50.} ', &
                                   '{50.}', .true., 50.)

    call mirr%set_param(           'mirr',            'multi',  'Perform mirroring', &
                                   'Whether to mirror and along which axis(no|x|y){no}', &
                                   '(no|x|y){no}', .false., 'no')

    call ml_reg%set_param(         'ml_reg',          'binary', 'ML regularization', &
                                   'Regularization (ML-style) based on the signal power(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call ml_reg_chunk%set_param(   'ml_reg_chunk',    'binary', 'Subset ML regularization', &
                                   'Subset Regularization (ML-style) based on the signal power(yes|no){no}', &
                                   '(yes|no){no}', .false., 'yes')

    call ml_reg_pool%set_param(    'ml_reg_pool',     'binary', 'Pool ML regularization', &
                                   'Pool Regularization (ML-style) based on the signal power(yes|no){no}', &
                                   '(yes|no){no}', .false., 'yes')

    call moldiam%set_param(        'moldiam',         'num',    'Molecular diameter', &
                                   'Molecular diameter(in Angstroms)', &
                                   'In Angstroms', .false., 0.)

    call moldiam_max%set_param(    'moldiam_max',     'num',    'Max molecular diameter', &
                                   'Max molecular diameter(in Angstroms)', &
                                   'In Angstroms', .false., 300.)

    call mskdiam%set_param(        'mskdiam',         'num',    'Mask diameter', &
                                   'Mask diameter (in A) for application of a soft-edged circular mask to remove background noise', &
                                   'mask diameter in A', .true., 0.)

    call mskfile%set_param(        'mskfile',         'file',   'Input mask file', &
                                   'Input mask file to apply to reference volume(s) before projection', &
                                   'e.g. automask.mrc from postprocess', .false., 'mskfile.mrc')

    call mul%set_param(            'mul',             'num',    'Multiplication factor', &
                                   'Multiplication factor{1.}', &
                                   '{1.}', .false., 1.)

    call nboxes_max%set_param(     'nboxes_max',      'num',    'Max # boxes per micrograph', &
                                   'Max # boxes per micrograph', &
                                   'Max # boxes', .false., 500.)

    call nchunks%set_param(        'nchunks',         'num',    'Number of subsets to classify simultaneously', &
                                   'Maximum number of particles subsets (chunks) to classify simultaneously', &
                                   '# of subsets', .true., 2.)

    call nchunksperset%set_param(  'nchunksperset',   'num',    'Number of subsets to group', &
                                   'Number of particles subsets (chunks) to group into a independent set', &
                                   '# of subsets per set', .false., 1.)

    call ncls%set_param(           'ncls',            'num',    'Number of 2D clusters', &
                                   'Number of groups to sort the particles into prior to averaging to create 2D class averages with improved SNR', &
                                   '# 2D clusters', .true., 200.)

    call ncls_start%set_param(     'ncls_start',      'num',    'Number of 2D clusters per subset of particles', &
                                   'Number of class averages used in the independent 2D analysis of each subset of particles', &
                                   '# 2D clusters / subset', .true., 50.)

    call neg%set_param(            'neg',             'binary', 'Invert contrast', &
                                   'Invert contrast(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call nparts%set_param(         'nparts',          'num',    'Number of computing nodes', &
                                   'Number of partitions for distributed memory execution. One part typically corresponds to one CPU socket in the distributed system. On a single-socket machine there may be speed benefits to dividing the jobs into a few (2-4) partitions, depending on memory capacity', &
                                   'divide job into # parts', .true., 1.0)

    call nparts_chunk%set_param(   'nparts_chunk',    'num',    'Number of computing nodes per subset', &
                                   'Number of computing nodes allocated to 2D analysis of each particles subset (chunk){1}', &
                                   '# of nodes per subset{1}', .false., 1.0)

    call nparts_pool%set_param(    'nparts_pool',     'num',    'Number of computing nodes for the pooled subsets', &
                                   'Number of computing nodes allocated to 2D analysis of the pooled particles subsets', &
                                   '# of nodes for the pooled subsets', .false., 2.0)

    call nptcls%set_param(         'nptcls',          'num',    'Number of particles', &
                                   'Number of particle images', &
                                   '# particles', .true., 0.)

    call nptcls_per_cls%set_param( 'nptcls_per_cls',  'num',    'Number of particles per cluster', &
                                   'Initial number of particles per cluster{35}', &
                                   '# initial particles per cluster{35}', .false., 35.)

    call nptcls_per_cls_cleanup2D%set_param( 'nptcls_per_cls',  'num',    'Number of particles per cluster', &
                                   'Number of particles per cluster{500}', &
                                   '# particles per cluster{500}', .false., 500.)

    call nran%set_param(           'nran',            'num',    'Number of random samples', &
                                   'Number of entries to randomly sample', &
                                   '# random samples', .false., 0.)

    call nrestarts%set_param(      'nrestarts',       'num',    'Number of restarts', &
                                   'Number of program restarts to execute{1}', &
                                   '# restarts{1}', .false., 1.0)

    call nsample%set_param(        'nsample',         'num',    'Number of particles to sample', &
                                   'Number of particles to sample each iteration', &
                                   '# particles to sample', .false., 0.)

    call nsig%set_param(           'nsig',            'num',    'Number of sigmas for outlier removal', &
                                   'Number of standard deviations threshold for pixel outlier removal{6}', &
                                   '# standard deviations{6}', .false., 6.)

    call nspace%set_param(         'nspace',          'num',    'Number of projection directions', &
                                   'Number of projection directions used', &
                                   '# projections', .false., 2500.)

    call nstates%set_param(        'nstates',         'num',    'Number of states', &
                                   'Number of conformational/compositional states to reconstruct', &
                                   '# states to reconstruct', .false., 1.0)

    call nthr%set_param(           'nthr',            'num',    'Number of threads per computing node, give 0 if unsure', &
                                   'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of logical threads in a socket.', &
                                   '# shared-memory CPU threads', .true., 0.)

    call numlen%set_param(         'numlen',          'num',    'Length of number string', &
                                   'Length of number string', &
                                   '# characters', .false., 5.0)

    call nxpatch%set_param(        'nxpatch',         'num',    '# of patches along x-axis', &
                                   '# of patches along x-axis(3)', &
                                   '# x-patches{3}', .false., 3.)

    call nypatch%set_param(        'nypatch',         'num',    '# of patches along y-axis', &
                                   '# of patches along y-axis(3)', &
                                   '# y-patches{3}', .false., 3.)

    call objfun%set_param(         'objfun',          'multi',  'Objective function', &
                                   'Objective function(euclid|cc|prob){euclid}', &
                                   '(euclid|cc|prob){euclid}', .false., 'euclid')

    call oritab%set_param(         'oritab',          'file',   'Orientation and CTF parameter file', &
                                   'Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', &
                                   '.simple|.txt parameter file', .false., 'oritab'//trim(METADATA_EXT))

    call oritab2%set_param(        'oritab2',         'file',   '2nd orientation and CTF parameter file', &
                                   '2nd orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', &
                                   '.simple|.txt parameter file', .false., 'oritab2'//trim(METADATA_EXT))

    call oritype%set_param(        'oritype',         'multi',  'Oritype segment in project', &
                                   'Oritype segment in project(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}', &
                                   '(mic|stk|ptcl2D|cls2D|cls3D|ptcl3D|out|projinfo|jobproc|compenv){ptcl3D}', .false., 'ptcl3D')

    call outfile%set_param(        'outfile',         'file',   'Output orientation and CTF parameter file', &
                                   'Output Orientation and CTF parameter file in plain text (.txt) or SIMPLE project (*.simple) format', &
                                   '.simple|.txt parameter file', .false., 'outfile'//trim(METADATA_EXT))

    call outside%set_param(        'outside',         'binary', 'Extract outside stage boundaries', &
                                   'Extract boxes outside the micrograph boundaries(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call outstk%set_param(         'outstk',          'file',   'Output stack name', &
                                   'Output images stack name', &
                                   'e.g. outstk.mrc', .false., '')

    call outvol%set_param(         'outvol',          'file',   'Output volume name', &
                                   'Output volume name', &
                                   'e.g. outvol.mrc', .false., '')

    call particle_density%set_param('particle_density','multi', 'Particle density in micrographs', &
                                    'Particle density in micrographs(low|optimal|high){optimal}', &
                                    '(low|optimal|high){optimal}', .false., 'optimal')

    call pcontrast%set_param(      'pcontrast',       'multi',  'Input particle contrast', &
                                   'Input particle contrast(black|white){black}', &
                                   '(black|white){black}', .false., 'black')

    call pdbout%set_param(         'pdbout',          'file',   'Output PDB volume-centered molecule', &
                                   'Output coordinates file in PDB format for the volume-centered molecule', &
                                   'e.g. output.pdb', .false., 'pdbout.pdb')

    call pgrp%set_param(           'pgrp',            'str',    'Point-group symmetry', &
                                   'Point-group symmetry of particle(cn|dn|t|o|i){c1}', &
                                   'point-group(cn|dn|t|o|i){c1}', .true., 'c1')

    call pgrp_start%set_param(     'pgrp_start',      'str',    'Initital point-group symmetry', &
                                   'Initial point-group symmetry(cn|dn|t|o|i){c1}', &
                                   'point-group(cn|dn|t|o|i){c1}', .false., 'c1')

    call phaseplate%set_param(     'phaseplate',      'binary', 'Phase-plate images', &
                                   'Images obtained with Volta phase-plate(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call pick_roi%set_param(       'pick_roi',        'binary', 'Artefactual regions exclusion(new picker only)', &
                                   'Whether to exclude regions of disinterest(carbon, thick ice, new picker only)(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call picker%set_param(         'picker',          'multi',  'Which picker to use', &
                                   'Which picker to use(old|new|segdiam){new}', &
                                   '(old|new|segdiam){new}', .false., 'new')

    call pickrefs%set_param(       'pickrefs',        'file',   'Stack of class-averages/reprojections for picking', &
                                   'Stack of class-averages/reprojections for picking', &
                                   'e.g. pickrefs.mrc', .false., '')

    call projfile%set_param(       'projfile',        'file',   'Project file', &
                                   'SIMPLE projectfile', &
                                   'e.g. myproject.simple', .true., '')

    call projfile_target%set_param('projfile_target', 'file',   'Another project file', &
                                   'SIMPLE projectfile', &
                                   'e.g. myproject2.simple', .true., '')

    call projfile_merged%set_param('projfile_merged', 'file',   'Merged output project file', &
                                   'SIMPLE projectfile', &
                                   'e.g. merged.simple', .true., '')

    call projfile_ref%set_param(   'projfile_ref',    'file',   'Reference project file', &
                                   'SIMPLE projectfile', &
                                   'e.g. myproject2.simple', .true., '')

    call projname%set_param(       'projname',        'str',    'Project name', &
                                   'Name of project to create ./myproject/myproject.simple file for', &
                                   'e.g. to create ./myproject/myproject.simple', .true., '')

    call prune%set_param(          'prune',           'binary', 'Automated particles pruning', &
                                   'Whether to prune deselected particles(yes|no){no}', &
                                   'Automated particles pruning(yes|no){no}', .false., 'no')

    call pspecsz%set_param(        'pspecsz',         'num',    'Size of power spectrum', &
                                   'Size of power spectrum in pixels{512}', &
                                   'give # pixels{512}', .false., 512.)

    call qsys_name%set_param(      'qsys_name',       'multi',  'Queue system kind', &
                                   'Queue system kind(local|slurm|pbs|lsf)', &
                                   '(local|slurm|pbs|lsf)', .false., 'local')

    call qsys_partition%set_param( 'qsys_partition',  'str',    'Name of SLURM/PBS/LSF partition', &
                                   'Name of target partition of distributed computer system (SLURM/PBS/LSF)', &
                                   'give partition name', .false., '')

    call qsys_qos%set_param(       'qsys_qos',        'str',    'Schedule priority', &
                                   'Job scheduling priority (SLURM/PBS/LSF)', &
                                   'give priority', .false., '')

    call qsys_reservation%set_param('qsys_reservation','str',   'Name of reserved partition', &
                                    'Name of reserved target partition of distributed computer system (SLURM/PBS/LSF)', &
                                    'give your part', .false., '')

    call remap_cls%set_param(      'remap_cls',       'binary', 'Whether to remap 2D clusters', &
                                   'Whether to remap the number of 2D clusters(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call remove_chunks%set_param(  'remove_chunks',   'binary', 'Whether to remove subsets', &
                                   'Whether to remove subsets after completion(yes|no){yes}', &
                                   '(yes|no){yes}', .false., 'yes')

    call script%set_param(         'script',          'binary', 'Generate script for shared-mem exec on cluster', &
                                   'Generate script for shared-mem exec on cluster(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call sherr%set_param(          'sherr',           'num',    'Shift error half-width', &
                                   'Uniform rotational origin shift error half-width(in pixels)', &
                                   'shift error in pixels', .false., 0.)

    call sigma_est%set_param(      'sigma_est',       'multi',  'Sigma estimation method', &
                                   'Sigma estimation method(group|global){group}', &
                                   '(group|global){group}', .false., 'group')

    call smpd%set_param(           'smpd',            'num',    'Sampling distance', &
                                   'Distance between neighbouring pixels in Angstroms', &
                                   'pixel size in Angstroms', .true., 1.0)

    call smpd_downscale%set_param( 'smpd_downscale',  'num',    'Sampling distance after downscale', &
                                   'Distance between neighbouring pixels in Angstroms after downscale', &
                                   'pixel size in Angstroms{1.3}', .false., 1.3)

    call smpd_target%set_param(    'smpd_target',     'num',    'Target sampling distance', &
                                   'Distance between neighbouring pixels in Angstroms', &
                                   'pixel size in Angstroms', .true., 1.0)

    call star_datadir%set_param(   'star_datadir',    'file',   'STAR project data directory', &
                                   'Pathname of STAR image/data files', &
                                   'e.g. Micrographs', .false., '')

    call star_mic%set_param(       'star_mic',        'file',   'Micrographs STAR file name', &
                                   'Micrographs STAR-formatted filename', &
                                   'e.g. micrographs.star', .true., '')

    call star_model%set_param(     'star_model',      'file',   'Model STAR file name', &
                                   'Model STAR-formatted filename', &
                                   'e.g. model.star', .true., '')

    call star_ptcl%set_param(      'star_ptcl',       'file',   'Particles STAR file name', &
                                   'Particles STAR-formatted filename', &
                                   'e.g. particles.star', .true., '')

    call starfile%set_param(       'starfile',        'file',   'STAR-format file name', &
                                   'STAR-formatted filename', &
                                   'e.g. proj.star', .false., '')

    call startit%set_param(        'startit',         'num',    'First iteration', &
                                   'Index of first iteration when starting from a previous run', &
                                   'start iterations from here', .false., 1.0)

    call startype%set_param(       'startype',        'str',    'STAR-format export type', &
                                   'STAR experiment type used to define variables in export file', &
                                   'e.g. micrographs or class2d or refine3d', .false., '')

    call stepsz%set_param(         'stepsz',          'num',    'Steps size in A', &
                                   'Step size in A {10.} ', &
                                   '{10.}', .true., 10.)

    call stk%set_param(            'stk',             'file',   'Particle image stack', &
                                   'Particle image stack', &
                                   'xxx.mrc file with particles', .false., '')


    call stk_backgr%set_param(      'stk_backgr',      'file',   'background power spectra stack, eg NP_X_background_pspec.mrc', &
                                   'background power spectra stack', &
                                   'xxx.mrc file with bg spectra', .true., '')


    call stk_traj%set_param(       'stk_traj',         'file',   'Tracked NP image stack', &
                                   'xxx.mrc file with NP time-trajectory', &
                                   'trajectory.mrcs', .true., '')

    call stk2%set_param(           'stk2',            'file',   'Second Particle image stack', &
                                   'Particle image stack', &
                                   'xxx.mrc file with particles', .false., 'stk2.mrc')

    call stktab%set_param(         'stktab',          'file',   'List of per-micrograph particle stacks', &
                                   'List of per-micrograph particle stacks', &
                                   'stktab.txt file containing file names', .false., 'stktab.txt')

    call time_per_image%set_param( 'time_per_image',  'num',    'Time per image', &
                                   'Estimated time per image in seconds for forecasting total execution time{100}', &
                                   'in seconds{100}', .false., 100.)

    call total_dose%set_param(     'total_dose',      'num',    'Total exposure dose (e/Ang^2)', &
                                   'Total exposure dose (e/Ang^2)', &
                                   'in e/Ang^2', .false., 50.)

    call trs%set_param(            'trs',             'num',    'Maximum translational shift', &
                                   'Maximum half-width for bund-constrained search of rotational origin shifts', &
                                   'max shift per iteration in pixels{5}', .false., 5.0)

    call trs_mc%set_param(         'trs_mc',           'num',    'Max shift per iter in pixels{10.}', &
                                   'Maximum half-width for bund-constrained search of movie frame shifts', &
                                   'max shift per iter in pixels{10}', .false., 5.0)

    call tseries%set_param(        'tseries',         'binary', 'Stack is time-series', &
                                   'Stack is time-series(yes|no){no}', &
                                   '(yes|no){no}', .false., 'no')

    call update_frac%set_param(    'update_frac',     'num',    'Fractional update per iteration', &
                                   'Fraction of particles to update per iteration in incremental learning scheme for accelerated convergence rate(0.1-0.5){1.}', &
                                   'update this fraction per iter(0.1-0.5){1.0}', .false., 1.0)

    call user_account%set_param(   'user_account',    'str',    'User account name in SLURM/PBS/LSF', &
                                   'User account name in SLURM/PBS/LSF system', &
                                   'e.g. Account084', .false., '')

    call user_email%set_param(     'user_email',      'str',    'Your e-mail address', &
                                   'Your e-mail address', &
                                   'e.g. myname@uni.edu', .false., '')

    call user_project%set_param(   'user_project',    'str',    'User project name in SLURM/PBS/LSF', &
                                   'User project name in SLURM/PBS/LSF system', &
                                   'e.g. Project001', .false., '')

    call vol_dim%set_param(        'vol_dim',         'num',    'Simulated volume dimensions', &
                                   'Dimensions of the simulated volume in voxels', &
                                   '# dimensions of the simulated volume', .false., 0.)

    call walltime%set_param(       'walltime',        'num',    'Walltime', &
                                   'Maximum execution time for job scheduling and management(23h59mins){86340}', &
                                   'in seconds(23h59mins){86340}', .false., 86340.)

    call wcrit%set_param(          'wcrit',           'multi',  'Correlation to weights conversion scheme', &
                                   'Correlation to weights conversion scheme(softmax|zscore|sum|cen|exp|inv|no){softmax}', &
                                   '(softmax|zscore|sum|cen|exp|inv|no){softmax}', .false., 'softmax')

    call width%set_param(          'width',           'num',    'Falloff of inner mask', &
                                   'Number of cosine edge pixels of inner mask in pixels', &
                                   '# pixels cosine edge{10}', .false., 10.)

    call wiener%set_param(         'wiener',          'multi',  'Wiener restoration', &
                                   'Wiener restoration, full or partial (full|partial){full}', &
                                   '(full|partial){full}', .false., 'full')
end subroutine set_ui_params

end module simple_ui_params_common
