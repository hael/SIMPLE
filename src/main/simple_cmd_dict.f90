! command line dictonary
module simple_cmd_dict
include 'simple_lib.f08'
implicit none

public :: print_cmdline, test_cmd_dict, print_cmd_key_descr
private
#include "simple_local_flags.inc"
integer, parameter :: NMAX=300
type(chash)        :: chdict
logical            :: initialised=.false.

contains

    subroutine init_cmd_dict
        call chdict%new(NMAX)
        call chdict%push('acf',           'calculate autocorrelation function(yes|no){no}')
        call chdict%push('alpha',         'oversampling factor{sqrt(2)}')
        call chdict%push('amsklp',        'low-pass limit for envelope mask generation(in A)')
        call chdict%push('angastunit',    'angle of astigmatism unit (radians|degrees){degrees}')
        call chdict%push('angerr',        'angular error(in degrees){0}')
        call chdict%push('append',        'append in context of files(yes|no){no}')
        call chdict%push('astigerr',      'astigmatism error(in microns)')
        call chdict%push('astigtol',      'expected (tolerated) astigmatism(in microns){0.05}')
        call chdict%push('async',         'asynchronous mode of operation(yes|no){no}')
        call chdict%push('automatic',     'automatic threshold selection for edge detection(yes|no)')
        call chdict%push('automsk',       'envelope masking(yes|no|cavg){no}')
        call chdict%push('autoscale',     'automatic down-scaling(yes|no)')
        call chdict%push('avg',           'calculate average(yes|no)')
        call chdict%push('bfac',          'bfactor for sharpening/low-pass filtering(in A**2){200.}')
        call chdict%push('bfacerr',       'bfactor error in simulated images(in A**2){0}')
        call chdict%push('balance',       'max pop for balancing restraint{0}')
        call chdict%push('bin',           'binarise image(yes|no){no}')
        call chdict%push('binwidth',      'binary layers grown for molecular envelope(in pixels){1}')
        call chdict%push('box',           'square image size(in pixels)')
        call chdict%push('box_extract',   'square image box size for extraction(in pixels)')
        call chdict%push('boxconvsz',     'size of box used for box-convolution(in pixels)')
        call chdict%push('boxfile',       'file with EMAN particle coordinates(.txt)')
        call chdict%push('boxtab',        'table (text file) of files with EMAN particle coordinates(.txt)')
        call chdict%push('bw_ratio',      'ratio between foreground-background pixels to obtain in edge detection')
        call chdict%push('cenlp',         'low-pass limit for binarisation in centering(in A){30 A}')
        call chdict%push('center',        'center image(s)/class average(s)/volume(s)(yes|no){no}')
        call chdict%push('chunksz',       '# images/orientations in chunk')
        call chdict%push('class',         'cluster identity')
        call chdict%push('classdoc',      'doc with per-class stats(.txt)')
        call chdict%push('classtats',     'calculate class population statistics(yes|no){no}')
        call chdict%push('clip',          'clipped image box size(in pixels)')
        call chdict%push('clustvalid',    'validate clustering(yes|homo|no){no}')
        call chdict%push('compare',       'do comparison(yes|no){no}')
        call chdict%push('corner',        'corner size(in pixels){0}')
        call chdict%push('countvox',      'count # voxels(yes|no){no}')
        call chdict%push('cs',            'spherical aberration constant(in mm){2.7}')
        call chdict%push('ctf',           'ctf flag(yes|no|flip)')
        call chdict%push('ctfreslim',     'ctf resolution (A) for setting state flag {8 A}')
        call chdict%push('ctfsq',         'apply ctf**2 to the images(yes|no){no}')
        call chdict%push('ctfstats',      'calculate ctf statistics(yes|no){no}')
        call chdict%push('projstats',     'calculate projection direction population statistics(yes|no){no}')
        call chdict%push('cube',          'side size(in pixels){0}')
        call chdict%push('dcrit_rel',     'critical distance relative to box(0-1){0.5}')
        call chdict%push('defocus',       'defocus(in microns){3.}')
        call chdict%push('deftab',        'file with CTF info(.txt|.simple)')
        call chdict%push('detector',      'edge detection algorithm (sobel|canny)')
        call chdict%push('dev',           'development flag for experimental code(yes|no){no}')
        call chdict%push('dfclose',       'close to focus limit for state flag(in microns){1}')
        call chdict%push('dffar',         'far from focus limit for state flag(in microns){4}')
        call chdict%push('dferr',         'defocus error(in microns){1.0}')
        call chdict%push('dfmax',         'maximum expected defocus(in microns){5.0}')
        call chdict%push('dfmin',         'minimum expected defocus(in microns){0.5}')
        call chdict%push('dfunit',        'defocus unit (A|microns){microns}')
        call chdict%push('dfstep',        'step size for defocus grid search(in microns){0.05}')
        call chdict%push('dir',           'directory')
        call chdict%push('dir_movies',    'grab .mrc/.mrcs files from here')
        call chdict%push('dir_ptcls',     'grab pre-micrograph stacks and docs from here')
        call chdict%push('dir_reject',    'move rejected files to here{rejected}')
        call chdict%push('dir_select',    'move selected files to here{selected}')
        call chdict%push('dir_target',    'put output here')
        call chdict%push('discrete',      'discrete(yes|no){no}')
        call chdict%push('diverse',       'diverse or not flag (yes|no){no}')
        call chdict%push('doclist',       'list of oritabs for different states')
        call chdict%push('dockmode',      'mode of docking(rot|shift|rotshift|refine){rotshift}')
        call chdict%push('doprint',       'whether to print(yes|no){no}')
        call chdict%push('dose_rate',     'dose rate(in e/A2/s)')
        call chdict%push('dynlp',         'automatic resolution limit update(yes|no){yes}')
        call chdict%push('e1',            'Euler 1 (in degrees){0}')
        call chdict%push('e2',            'Euler 2 (in degrees){0}')
        call chdict%push('e3',            'Euler 3 (in degrees){0}')
        call chdict%push('edge',          'edge size for softening molecular envelope(in pixels){6}')
        call chdict%push('endian',        'endiannesss of files(big|little|native){native}')
        call chdict%push('eps',           'learning rate{0.003}')
        call chdict%push('eo',            'use FSC for filtering and low-pass limit update(yes|aniso|no)')
        call chdict%push('errify',        'introduce error(yes|no){no}')
        call chdict%push('even',          'even orientation distribution(yes|no){no}')
        call chdict%push('exp_doc',       'specifying exp_time and dose_rate per tomogram')
        call chdict%push('startype',      'specifying startype for STAR formated projects (movies|micrographs|mcmicrographs|ctf_estimation|select|extract|class2d|init3dmodel|refine3d|post|all)')
        call chdict%push('exp_time',      'exposure time(in s)')
        call chdict%push('extr_init',     'initial extremal ratio(0-1){0.5}')
        call chdict%push('ext',           'file extension{.mrc}')
        call chdict%push('fbody',         'file body')
        call chdict%push('filetab',       'list of files(.txt)')
        call chdict%push('find',          'Fourier index')
        call chdict%push('fname',         'file name')
        call chdict%push('for3D',         'for 3D analysis(yes|no){yes}')
        call chdict%push('focusmsk',      'focused mask radius(in pixels)')
        call chdict%push('frac',          'fraction of ptcls(0-1){1}')
        call chdict%push('frac_outliers', 'fraction of outliers(0-1){0.0}')
        call chdict%push('fraca',         'fraction of amplitude contrast used for fitting CTF{0.1}')
        call chdict%push('fracdeadhot',   'fraction of dead or hot pixels{0.01}')
        call chdict%push('fromf',         'start frame index')
        call chdict%push('fromp',         'start ptcl index')
        call chdict%push('frcs',          'binary file with per-class/proj Fourier Ring Correlations(.bin)')
        call chdict%push('fsc',           'binary file with FSC info{fsc_state01.bin}')
        call chdict%push('ft2img',        'convert Fourier transform to real image of power(yes|no){no}')
        call chdict%push('gainref',       'gain reference')
        call chdict%push('grow',          '# binary layers to grow(in pixels)')
        call chdict%push('guinier',       'calculate Guinier plot(yes|no){no}')
        call chdict%push('hfun',          'function used for normalization(sigm|tanh|lin){sigm}')
        call chdict%push('hist',          'give variable for histogram plot')
        call chdict%push('hp',            'high-pass limit(in A)')
        call chdict%push('hp_fsc',        'FSC high-pass limit(in A)')
        call chdict%push('hp_ctf_estimate', 'high-pass limit 4 ctf_estimate(in A)')
        call chdict%push('iares',         'integer angular resolution{10}')
        call chdict%push('imgkind',       'type of image(ptcl|cavg|mic|movie){ptcl}')
        call chdict%push('infile',        'file with inputs(.txt)')
        call chdict%push('infile2',       'file with inputs(.txt)')
        call chdict%push('inner',         'inner mask radius(in pixels)')
        call chdict%push('job_memory_per_task', 'memory in MB per task in distributed exec (typically memory per socket)')
        call chdict%push('jumpsz',        'size of contigous segment')
        call chdict%push('keys',          'keys of values to print')
        call chdict%push('kv',            'acceleration voltage(in kV){300.}')
        call chdict%push('label',         'discrete label(class|state){class}')
        call chdict%push('lp',            'low-pass limit(in A)')
        call chdict%push('lp_ctf_estimate', 'low-pass limit 4 ctf_estimate(in A)')
        call chdict%push('lp_pick',       'low-pass limit 4 picker(in A)')
        call chdict%push('lplim_crit',    'corr criterion low-pass limit assignment(0.143-0.5){0.3}')
        call chdict%push('lpstart',       'start low-pass limit(in A){15}')
        call chdict%push('lpstop',        'stop low-pass limit(in A){8}')
        call chdict%push('masscen',       'center using binarisation and mass centering(yes|no){no}')
        call chdict%push('match_filt',    'matched filter on (yes|no){yes}')
        call chdict%push('maxits',        'maximum # iterations')
        call chdict%push('minp',          'minimum cluster population')
        call chdict%push('mirr',          'mirror(no|x|y){no}')
        call chdict%push('moldiam',       'molecular diameter(in A)')
        call chdict%push('msk',           'mask radius(in pixels)')
        call chdict%push('mskfile',       'maskfile.ext')
        call chdict%push('msktype',       'type of mask(hard|soft){soft}')
        call chdict%push('mul',           'origin shift multiplication factor{1}')
        call chdict%push('mw',            'molecular weight(in kD)')
        call chdict%push('msklist',       'table (text file) of mask volume files(.txt)')
        call chdict%push('ncls',          '# clusters')
        call chdict%push('ncls_start',    '# clusters required to start prime2D streaming')
        call chdict%push('ncunits',       '# computing units, can be < nparts{nparts}')
        call chdict%push('ndiscrete',     '# discrete orientations')
        call chdict%push('ndev',          '# deviations in one-cluster clustering{2.0}')
        call chdict%push('ndocs',         '# documents')
        call chdict%push('neg',           'invert contrast of images(yes|no)')
        call chdict%push('neigh',         'nearest neighbor refinement(yes|no)')
        call chdict%push('newbox',        'new box for scaling (by Fourier padding/clipping')
        call chdict%push('nframes',       '# frames{30}')
        call chdict%push('nframesgrp',    '# frames to group before motion_correct(Falcon 3){0}')
        call chdict%push('nnn',           '# nearest neighbors{200}')
        call chdict%push('noise',         'noise initialisation(yes|no){no}')
        call chdict%push('noise_norm',    'normalise based on sdev of background(yes|no){no}')
        call chdict%push('norec',         'do not reconstruct volume(s)(yes|no){no}')
        call chdict%push('norm',          'do statistical normalisation avg')
        call chdict%push('nparts',        '# partitions in distributed exection')
        call chdict%push('npix',          '# pixles/voxels in binary representation')
        call chdict%push('nptcls',        '# images in stk/# orientations in oritab')
        call chdict%push('nptcls_per_cls','# images per class for 2D streaming')
        call chdict%push('nran',          '# random images to select')
        call chdict%push('nrefs',         '# references used for picking{100}')
        call chdict%push('nrepeats',      '# times to restart workflow{1}')
        call chdict%push('nsig',          '# sigmas')
        call chdict%push('nspace',        '# projection directions')
        call chdict%push('nstates',       '# states to reconstruct')
        call chdict%push('nsub',          '# proj dirs in coarse grid search{300}')
        call chdict%push('nthr',          '# OpenMP threads{1}')
        call chdict%push('numlen',        'length of number string')
        call chdict%push('numlen_tomo',   'length of number string tomo series index{3}')
        call chdict%push('nvox',          '# voxels{0}')
        call chdict%push('objfun',        'objective function(cc|ccres|euclid){cc}')
        call chdict%push('offset',        'pixels offset{10}')
        call chdict%push('opt',           'optimiser (bfgs|simplex){bfgs}')
        call chdict%push('order',         'order ptcls according to correlation(yes|no){no}')
        call chdict%push('oritab',        'table of orientations(.txt|.simple)')
        call chdict%push('oritab2',       'table of orientations 2(.txt|.simple)')
        call chdict%push('oritab3D',      'table of 3D orientations(.txt|.simple)')
        call chdict%push('oritype',       'SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D|projinfo|jobproc|compenv)')
        call chdict%push('outer',         'outer mask radius(in pixels)')
        call chdict%push('outfile',       'output document')
        call chdict%push('outside',       'extract boxes outside the micrograph boundaries(yes|no){no}')
        call chdict%push('outstk',        'output image stack')
        call chdict%push('outstk2',       'output image stack 2')
        call chdict%push('outvol',        'output volume{outvol.ext}')
        call chdict%push('pcontrast',     'particle contrast(black|white){black}')
        call chdict%push('pdbfile',       'input PDB formatted file')
        call chdict%push('pgrp',          'point-group symmetry(cn|dn|t|o|i)')
        call chdict%push('phaseplate',    'images obtained with Volta phaseplate(yes|no){no}')
        call chdict%push('phrand',        'phase randomize(yes|no){no}')
        call chdict%push('phshiftunit',   'additional phase-shift unit (radians|degrees){radians}')
        call chdict%push('plaintexttab',  'plain text file of input parameters')
        call chdict%push('plot',          'make plot(yes|no){no}')
        call chdict%push('prg',           'SIMPLE program to execute')
        call chdict%push('projfile',      'SIMPLE *.simple project file')
        call chdict%push('projfile_target', 'another SIMPLE *.simple project file')
        call chdict%push('projname',      'Project name (for creation of projname.simple)')
        call chdict%push('pspecsz',       'size of power spectrum(in pixels)')
        call chdict%push('clsfrcs',       'Whether the projection frcs are to be read from class(yes|no){no}')
        call chdict%push('qsys_partition', 'Name of target partition of distributed computer system (SLURM/PBS)')
        call chdict%push('qsys_qos',      'job scheduling priority (SLURM/PBS)')
        call chdict%push('qsys_reservation', 'Name of reserved target partition of distributed computer system (SLURM/PBS)')
        call chdict%push('real_filter',   'real-space filter kind(median|average|bman)')
        call chdict%push('recvol_sigma',  'noise(sigma)-weighted volume reconstruction strategy(yes|no){no}')
        call chdict%push('refine',        'refinement mode(snhc|single|multi|greedy_single|greedy_multi|cluster|clustersym|cont_single){no}')
        call chdict%push('refs',          'initial2Dreferences.ext')
        call chdict%push('remap_cls',    'remove empty and renumber and/or expand # clusters(yes|no){no}')
        call chdict%push('rm_outliers',   'remove outliers{yes}')
        call chdict%push('rnd',           'random(yes|no){no}')
        call chdict%push('rrate',         'randomization rate{0.8}')
        call chdict%push('scale',         'image scale factor{1}')
        call chdict%push('scale2',        'image scale factor 2{1}')
        call chdict%push('shalgn',        'do 2D shift alignment(yes|no){no}')
        call chdict%push('shbarrier',     'use shift search barrier constraint(yes|no){yes}')
        call chdict%push('shell_norm',    'normalise based on power spectrum (yes|no){no}')
        call chdict%push('shellw',        'shell-weight reconstruction (yes|no){no}')
        call chdict%push('sherr',         'shift error(in pixels){2}')
        call chdict%push('sigma2_fudge',  'fudge factor for sigma2_noise {1}')
        call chdict%push('single',        'simulate a single image(yes|no){no}')
        call chdict%push('smpd',          'sampling distance, same as EMANs apix(in A)')
        call chdict%push('snr',           'signal-to-noise ratio')
        call chdict%push('soften',        'soften envelope with cosine edge(yes|no){no}')
        call chdict%push('speckind',      'power spectrum kind(real|power|sqrt|log|phase){sqrt}')
        call chdict%push('startit',       'start iterating from here')
        call chdict%push('starfile',      'STAR-formatted project file')
        call chdict%push('state',         'state to extract')
        call chdict%push('state2split',   'state group to split')
        call chdict%push('stats',         'provide statistics(yes|no){yes}')
        call chdict%push('stepsz',        'size of step{0}')
        call chdict%push('stk',           'particle stack with all images(ptcls.ext)')
        call chdict%push('stktab',        'list of per-micrograph stacks')
        call chdict%push('stk2',          'stack 2 (in selection map: selected(cavgs).ext)')
        call chdict%push('stk3',          'stack 3 (in selection map: (cavgs)2selectfrom.ext)')
        call chdict%push('stk_backgr',    'stack with image for background subtraction')
        call chdict%push('stream',        'sream (real time) execution mode(yes|no){no}')
        call chdict%push('symrnd',        'randomize over symmetry operations(yes|no){no}')
        call chdict%push('szsn',          'size of stochastic neighborhood{5}')
        call chdict%push('taper_edges',   'to taper edges(yes|no){no}')
        call chdict%push('time_inactive', 'end time limit(in seconds){7200}')
        call chdict%push('time_per_image', 'Estimated time per image in seconds for forecasting total execution time')
        call chdict%push('thres',         'threshold (binarisation: 0-1; distance filer: # pixels; post-proc: pix val)')
        call chdict%push('thres_low',     'lower threshold for canny edge detection')
        call chdict%push('thres_up',      'upper threshold for canny edge detection')
        call chdict%push('tof',           'stop frame index')
        call chdict%push('tomo',          'tomography mode(yes|no){no}')
        call chdict%push('tomoseries',    'filetable of filetables of tomograms')
        call chdict%push('top',           'stop particle index')
        call chdict%push('trs',           'maximum halfwidth shift(in pixels)')
        call chdict%push('trsstats',      'provide origin shift statistics(yes|no){no}')
        call chdict%push('tseries',       'images represent a time-series(yes|no){no}')
        call chdict%push('unidoc',        'unified resources and orientations doc')
        call chdict%push('update_frac',   'fraction of particles to update(0.-1.){1.}')
        call chdict%push('user_account',  'Name of user account in distributed computer system (SLURM/PBS)')
        call chdict%push('user_email',    'Your e-mail address')
        call chdict%push('user_project',  'Name of project in distributed computer system (SLURM/PBS)')
        call chdict%push('verbose',       'verbosity flag (yes|no){no}')
        call chdict%push('vis',           'visualise(yes|no)')
        call chdict%push('vol1',          'input volume no1(invol1.ext)')
        call chdict%push('vol2',          'input volume no2(invol2.ext)')
        call chdict%push('vol_filt',      'input filter volume(vol_filt.ext)')
        call chdict%push('vollist',       'table (text file) of volume files(.txt)')
        call chdict%push('voltab',        'table (text file) of volume files(.txt)')
        call chdict%push('voltab2',       'table 2(text file) of volume files(.txt)')
        call chdict%push('which_iter',    'iteration nr')
        call chdict%push('width',         'falloff of inner mask or filter(in pixels){10}')
        call chdict%push('wiener',        'Wiener restoration mode(full|highres){highres}')
        call chdict%push('winsz',         'half-width of window for real-space filter(in pixels)')
        call chdict%push('wtype',         'type of orientation weights (factorial|flat){factorial}')
        call chdict%push('xcoord',        'x coordinate{0}')
        call chdict%push('xdim',          'x dimension(in pixles)')
        call chdict%push('xsh',           'x shift(in pixels){0}')
        call chdict%push('ycoord',        'y coordinate{0}')
        call chdict%push('ydim',          'y dimension(in pixles)')
        call chdict%push('ysh',           'y shift(in pixels){0}')
        call chdict%push('zero',          'zeroing(yes|no){no}')
        call chdict%push('zsh',           'z shift(in pixels){0}')
        initialised = .true.
    end subroutine init_cmd_dict

    subroutine print_cmd_key_descr( fname )
        character(len=*), intent(in) :: fname
        ! initialise if needed
        if( .not. initialised ) call init_cmd_dict
        ! sort hash
        call chdict%sort
        ! write hash
        call chdict%write(fname)
    end subroutine print_cmd_key_descr

    subroutine print_cmdline( keys_required, keys_optional, fhandle, distr )
        character(len=KEYLEN), optional, intent(in) :: keys_required(:), keys_optional(:)
        integer,               optional, intent(in) :: fhandle
        logical,               optional, intent(in) :: distr
        character(len=KEYLEN), allocatable :: sorted_keys(:)
        integer :: nreq, nopt
        logical :: ddistr
        ddistr = .false.
        if( present(distr) ) ddistr = distr
        ! initialise if needed
        if( .not. initialised ) call init_cmd_dict
        write(logfhandle,'(a)') 'USAGE:'
        if( ddistr )then
            write(logfhandle,'(a)') 'bash-3.2$ simple_distr_exec prg=simple_program key1=val1 key2=val2 ...'
        else
            write(logfhandle,'(a)') 'bash-3.2$ simple_exec prg=simple_program key1=val1 key2=val2 ...'
        endif
        ! print required
        if( present(keys_required) )then
            nreq =  size(keys_required)
            DebugPrint  '# required keys provided: ', nreq
            if( nreq > 0 )then
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'REQUIRED'
                allocate(sorted_keys(nreq), source=keys_required, stat=alloc_stat)
                if(alloc_stat /= 0)call allocchk("simple_cmd_dict::print_cmdline  sorted_keys 1",alloc_stat)
                call lexSort(sorted_keys)
                call chdict%print_key_val_pairs(sorted_keys, fhandle)
                deallocate(sorted_keys)
            endif
        endif
        ! print optionals
        if( present(keys_optional) )then
            nopt = size(keys_optional)
            DebugPrint  '# optional keys provided: ', nopt
            if( nopt > 0 )then
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'OPTIONAL'
                allocate(sorted_keys(nopt), source=keys_optional, stat=alloc_stat)
                if(alloc_stat /= 0)call allocchk("simple_cmd_dict::print_cmdline  sorted_keys 2",alloc_stat)
                call lexSort(sorted_keys)
                call chdict%print_key_val_pairs(sorted_keys, fhandle)
                deallocate(sorted_keys)
            endif
        endif
        write(logfhandle,'(a)') ''
    end subroutine print_cmdline

    subroutine test_cmd_dict
        character(len=32) :: keys_required(5), keys_optional(5)
        keys_required(1) = 'vol1'
        keys_required(2) = 'stk'
        keys_required(3) = 'smpd'
        keys_required(4) = 'msk'
        keys_required(5) = 'oritab'
        keys_optional(1) = 'refine'
        keys_optional(2) = 'automsk'
        keys_optional(3) = 'trs'
        keys_optional(4) = 'ctf'
        keys_optional(5) = 'pgrp'
        call print_cmdline
        call print_cmdline( keys_required )
        call print_cmdline( keys_optional=keys_optional )
        call print_cmdline( keys_required, keys_optional )
        call print_cmd_key_descr('cmd_key_descr_from_test_cmd_dict.txt')
    end subroutine test_cmd_dict

end module simple_cmd_dict
