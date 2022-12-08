module simple_private_prgs
include 'simple_lib.f08'
implicit none

public :: make_private_user_interface, print_private_cmdline, get_private_keys_required, get_n_private_keys_required, print_cmdline_oldschool
private
#include "simple_local_flags.inc"

! private program instance
type simple_private_prg
    private
    character(len=:), allocatable :: name
    character(len=KEYLEN)         :: keys_required(MAXNKEYS)
    character(len=KEYLEN)         :: keys_optional(MAXNKEYS)
    integer                       :: nreq=0, nopt=0
  contains
    procedure :: set_name
    procedure :: push_req_key
    procedure :: push_opt_key
end type simple_private_prg

! array of simple_private_exec program specifications
integer, parameter       :: NMAX_PRIVATE_PRGS = 100
integer                  :: n_private_prgs    = 0
type(simple_private_prg) :: private_prgs(NMAX_PRIVATE_PRGS)

! old-school cmd_dict (command line dictionary), used by simple_private_exec
integer, parameter :: NMAX_CMD_DICT = 300
type(chash)        :: cmd_dict
logical            :: cmd_dict_initialised = .false.

contains

    ! instance methods

    subroutine set_name( self, name )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: name
        allocate(self%name, source=trim(name))
    end subroutine set_name

    subroutine push_req_key( self, key )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: key
        self%nreq = self%nreq + 1
        self%keys_required(self%nreq) = trim(key)
    end subroutine push_req_key

    subroutine push_opt_key( self, key )
        class(simple_private_prg), intent(inout) :: self
        character(len=*),          intent(in)    :: key
        self%nopt = self%nopt + 1
        self%keys_optional(self%nopt) = trim(key)
    end subroutine push_opt_key

    ! class methods

    subroutine make_private_user_interface
        call init_cmd_dict
        call new_private_prgs
    end subroutine make_private_user_interface

    function get_n_private_keys_required( prg ) result( nreq )
        character(len=*), intent(in)  :: prg
        integer :: nreq, iprg, i
        nreq = 0
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nreq = private_prgs(iprg)%nreq
    end function get_n_private_keys_required

    function get_private_keys_required( prg ) result( keys_required )
        character(len=*), intent(in)  :: prg
        type(str4arr), allocatable    :: keys_required(:)
        integer :: iprg, i, nreq
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nreq = private_prgs(iprg)%nreq
        if( nreq == 0 ) return
        allocate(keys_required(nreq))
        do i=1,nreq
            allocate(keys_required(i)%str, source=trim(private_prgs(iprg)%keys_required(i)))
        end do
    end function get_private_keys_required

    function get_n_private_keys_optional( prg ) result( nopt )
        character(len=*), intent(in)  :: prg
        integer :: nopt, iprg, i
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nopt = private_prgs(iprg)%nopt
    end function get_n_private_keys_optional

    function get_private_keys_optional( prg, nopt ) result( keys_optional )
        character(len=*), intent(in)  :: prg
        integer,          intent(out) :: nopt
        type(str4arr), allocatable    :: keys_optional(:)
        integer :: iprg, i
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 ) return
        nopt = private_prgs(iprg)%nopt
        if( nopt == 0 ) return
        allocate(keys_optional(nopt))
        do i=1,nopt
            allocate(keys_optional(i)%str, source=trim(private_prgs(iprg)%keys_optional(i)))
        end do
    end function get_private_keys_optional

    subroutine print_private_cmdline( prg )
        character(len=*), intent(in) :: prg
        character(len=KEYLEN), allocatable :: sorted_keys(:)
        integer :: iprg, i, nreq, nopt
        iprg = 0
        do i=1,n_private_prgs
            if( trim(prg) .eq. private_prgs(i)%name )then
                iprg = i
                exit
            endif
        end do
        if( iprg == 0 )then
            THROW_WARN(trim(prg)//' lacks description in the private_prgs class')
            return
        endif
        write(logfhandle,'(a)') 'USAGE:'
        write(logfhandle,'(a)') 'bash-3.2$ simple_private_exec prg=simple_program key1=val1 key2=val2 ...'
        ! print required
        nreq = private_prgs(iprg)%nreq
        if( nreq > 0 )then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'REQUIRED'
            allocate(sorted_keys(nreq), source=private_prgs(iprg)%keys_required(:nreq))
            call lexSort(sorted_keys)
            call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
            deallocate(sorted_keys)
        endif
        ! print optionals
        nopt = private_prgs(iprg)%nopt
        if( nopt > 0 )then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') 'OPTIONAL'
            allocate(sorted_keys(nopt), source=private_prgs(iprg)%keys_optional(:nopt))
            call lexSort(sorted_keys)
            call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
            deallocate(sorted_keys)
        endif
        write(logfhandle,'(a)') ''
    end subroutine print_private_cmdline

    subroutine print_cmdline_oldschool( keys_required, keys_optional, distr )
        character(len=KEYLEN), optional, intent(in) :: keys_required(:), keys_optional(:)
        logical,               optional, intent(in) :: distr
        character(len=KEYLEN), allocatable :: sorted_keys(:)
        integer :: nreq, nopt
        logical :: ddistr
        ddistr = .false.
        if( present(distr) ) ddistr = distr
        ! initialise if needed
        if( .not. cmd_dict_initialised ) call init_cmd_dict
        write(logfhandle,'(a)') 'USAGE:'
        if( ddistr )then
            write(logfhandle,'(a)') 'bash-3.2$ simple_distr_exec prg=simple_program key1=val1 key2=val2 ...'
        else
            write(logfhandle,'(a)') 'bash-3.2$ simple_exec prg=simple_program key1=val1 key2=val2 ...'
        endif
        ! print required
        if( present(keys_required) )then
            nreq =  size(keys_required)
            if( nreq > 0 )then
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'REQUIRED'
                allocate(sorted_keys(nreq), source=keys_required)
                call lexSort(sorted_keys)
                call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
                deallocate(sorted_keys)
            endif
        endif
        ! print optionals
        if( present(keys_optional) )then
            nopt = size(keys_optional)
            if( nopt > 0 )then
                write(logfhandle,'(a)') ''
                write(logfhandle,'(a)') 'OPTIONAL'
                allocate(sorted_keys(nopt), source=keys_optional)
                call lexSort(sorted_keys)
                call cmd_dict%print_key_val_pairs(logfhandle, sorted_keys)
                deallocate(sorted_keys)
            endif
        endif
        write(logfhandle,'(a)') ''
    end subroutine print_cmdline_oldschool

    subroutine init_cmd_dict
        call cmd_dict%new(NMAX_CMD_DICT)
        call cmd_dict%push('acf',           'calculate autocorrelation function(yes|no){no}')
        call cmd_dict%push('algorithm',     'Algorithm to be used')
        call cmd_dict%push('alpha',         'oversampling factor{sqrt(2)}')
        call cmd_dict%push('amsklp',        'low-pass limit for envelope mask generation(in A)')
        call cmd_dict%push('angastunit',    'angle of astigmatism unit (radians|degrees){degrees}')
        call cmd_dict%push('angerr',        'angular error(in degrees){0}')
        call cmd_dict%push('append',        'append in context of files(yes|no){no}')
        call cmd_dict%push('astigerr',      'astigmatism error(in microns)')
        call cmd_dict%push('astigtol',      'expected (tolerated) astigmatism(in microns){0.05}')
        call cmd_dict%push('async',         'asynchronous mode of operation(yes|no){no}')
        call cmd_dict%push('automatic',     'automatic threshold selection for edge detection(yes|no)')
        call cmd_dict%push('automsk',       'envelope masking(yes|tight|no){no}')
        call cmd_dict%push('autoscale',     'automatic down-scaling(yes|no)')
        call cmd_dict%push('avg',           'calculate average(yes|no)')
        call cmd_dict%push('bfac',          'bfactor for sharpening/low-pass filtering(in A**2){200.}')
        call cmd_dict%push('bfacerr',       'bfactor error in simulated images(in A**2){0}')
        call cmd_dict%push('balance',       'max pop for balancing restraint{0}')
        call cmd_dict%push('bin',           'binarize image(yes|no){no}')
        call cmd_dict%push('binwidth',      'binary layers grown for molecular envelope(in pixels){1}')
        call cmd_dict%push('box',           'square image size(in pixels)')
        call cmd_dict%push('box_extract',   'square image box size for extraction(in pixels)')
        call cmd_dict%push('boxfile',       'file with EMAN particle coordinates(.txt)')
        call cmd_dict%push('boxtab',        'table (text file) of files with EMAN particle coordinates(.txt)')
        call cmd_dict%push('bw_ratio',      'ratio between foreground-background pixels to obtain in edge detection')
        call cmd_dict%push('cenlp',         'low-pass limit for binarisation in centering(in A){30 A}')
        call cmd_dict%push('center',        'center image(s)/class average(s)/volume(s)(yes|no){no}')
        call cmd_dict%push('chunksz',       '# images/orientations in chunk')
        call cmd_dict%push('class',         'cluster identity')
        call cmd_dict%push('classdoc',      'doc with per-class stats(.txt)')
        call cmd_dict%push('classtats',     'calculate class population statistics(yes|no){no}')
        call cmd_dict%push('clip',          'clipped image box size(in pixels)')
        call cmd_dict%push('clustvalid',    'validate clustering(yes|homo|no){no}')
        call cmd_dict%push('combine_eo',    'whether combined e/o references have been combined for alignment(yes|no){no}')
        call cmd_dict%push('compare',       'do comparison(yes|no){no}')
        call cmd_dict%push('corner',        'corner size(in pixels){0}')
        call cmd_dict%push('countvox',      'count # voxels(yes|no){no}')
        call cmd_dict%push('cs',            'spherical aberration constant(in mm){2.7}')
        call cmd_dict%push('ctf',           'ctf flag(yes|no|flip)')
        call cmd_dict%push('ctfpatch',      'whether to perform patch-based ctf estimation(yes|no){yes}')
        call cmd_dict%push('ctfresthreshold','ctf resolution (A) rejection theshold{30 A}')
        call cmd_dict%push('ctfsq',         'apply ctf**2 to the images(yes|no){no}')
        call cmd_dict%push('ctfstats',      'calculate ctf statistics(yes|no){no}')
        call cmd_dict%push('projstats',     'calculate projection direction population statistics(yes|no){no}')
        call cmd_dict%push('cube',          'side size(in pixels){0}')
        call cmd_dict%push('dcrit_rel',     'critical distance relative to box(0-1){0.5}')
        call cmd_dict%push('defocus',       'defocus(in microns){3.}')
        call cmd_dict%push('deftab',        'file with CTF info(.txt|.simple)')
        call cmd_dict%push('detector',      'edge detection algorithm (sobel|canny)')
        call cmd_dict%push('dev',           'development flag for experimental code(yes|no){no}')
        call cmd_dict%push('dfclose',       'close to focus limit for state flag(in microns){1}')
        call cmd_dict%push('dffar',         'far from focus limit for state flag(in microns){4}')
        call cmd_dict%push('dferr',         'defocus error(in microns){1.0}')
        call cmd_dict%push('dfmax',         'maximum expected defocus(in microns){5.0}')
        call cmd_dict%push('dfmin',         'minimum expected defocus(in microns){0.2}')
        call cmd_dict%push('dfunit',        'defocus unit (A|microns){microns}')
        call cmd_dict%push('dir',           'directory')
        call cmd_dict%push('dir_movies',    'grab .mrc/.mrcs files from here')
        call cmd_dict%push('dir_prev',      'grab previous project info & files')
        call cmd_dict%push('dir_ptcls',     'grab pre-micrograph stacks and docs from here')
        call cmd_dict%push('dir_refine',    'directory with opeaks_part*.bin; e.g. 1_refine3D')
        call cmd_dict%push('dir_reject',    'move rejected files to here{rejected}')
        call cmd_dict%push('dir_select',    'move selected files to here{selected}')
        call cmd_dict%push('dir_target',    'put output here')
        call cmd_dict%push('discrete',      'discrete(yes|no){no}')
        call cmd_dict%push('diverse',       'diverse or not flag (yes|no){no}')
        call cmd_dict%push('doclist',       'list of oritabs for different states')
        call cmd_dict%push('doprint',       'whether to print(yes|no){no}')
        call cmd_dict%push('draw_color',    'color of the cross that identify the picked particle {white|black}')
        call cmd_dict%push('dynlp',         'automatic resolution limit update(yes|no){yes}')
        call cmd_dict%push('e1',            'Euler 1 (in degrees){0}')
        call cmd_dict%push('e2',            'Euler 2 (in degrees){0}')
        call cmd_dict%push('e3',            'Euler 3 (in degrees){0}')
        call cmd_dict%push('edge',          'edge size for softening molecular envelope(in pixels){6}')
        call cmd_dict%push('eer_fraction',  'Number of raw frames to fraction together{30}')
        call cmd_dict%push('eer_upsampling','Desired output sampling: 1=4K; 2=8K{1}')
        call cmd_dict%push('element',       'atom element')
        call cmd_dict%push('endian',        'endiannesss of files(big|little|native){native}')
        call cmd_dict%push('eps',           'learning rate{0.003}')
        call cmd_dict%push('eo',            'use FSC for filtering and low-pass limit update(yes|aniso|no)')
        call cmd_dict%push('even',          'even orientation distribution(yes|no){no}')
        call cmd_dict%push('startype',      'specifying startype for STAR formated projects (movies|micrographs|mcmicrographs|ctf_estimation|select|extract|class2d|init3dmodel|refine3d|post|all)')
        call cmd_dict%push('extr_init',     'initial extremal ratio(0-1){0.5}')
        call cmd_dict%push('ext',           'file extension{.mrc}')
        call cmd_dict%push('fbody',         'file body')
        call cmd_dict%push('fill_holes',    'fill holes(yes|no){no}')
        call cmd_dict%push('filetab',       'list of files(.txt)')
        call cmd_dict%push('filter',        'filter type{no}')
        call cmd_dict%push('find',          'Fourier index')
        call cmd_dict%push('fname',         'file name')
        call cmd_dict%push('for3D',         'for 3D analysis(yes|no){yes}')
        call cmd_dict%push('focusmskdiam',  'focused mask diameter(in A)')
        call cmd_dict%push('frac',          'fraction of ptcls(0-1){1}')
        call cmd_dict%push('frac_outliers', 'fraction of outliers(0-1){0.0}')
        call cmd_dict%push('fraca',         'fraction of amplitude contrast used for fitting CTF{0.1}')
        call cmd_dict%push('fracdeadhot',   'fraction of dead or hot pixels{0.01}')
        call cmd_dict%push('fraction_dose_target','EER fraction dose target(in e/A^2)')
        call cmd_dict%push('fromp',         'start ptcl index')
        call cmd_dict%push('frcs',          'binary file with per-class/proj Fourier Ring Correlations(.bin)')
        call cmd_dict%push('fsc',           'binary file with FSC info{fsc_state01.bin}')
        call cmd_dict%push('ft2img',        'convert Fourier transform to real image of power(yes|no){no}')
        call cmd_dict%push('gainref',       'gain reference')
        call cmd_dict%push('groupframes',   'whether to perform patch weighted frame averaging during patch-based motion correction(yes)|no){no}')
        call cmd_dict%push('grow',          '# binary layers to grow(in pixels)')
        call cmd_dict%push('guinier',       'calculate Guinier plot(yes|no){no}')
        call cmd_dict%push('hfun',          'function used for normalization(sigm|tanh|lin){sigm}')
        call cmd_dict%push('hist',          'give variable for histogram plot')
        call cmd_dict%push('hp',            'high-pass limit(in A)')
        call cmd_dict%push('hp_fsc',        'FSC high-pass limit(in A)')
        call cmd_dict%push('hp_ctf_estimate', 'high-pass limit 4 ctf_estimate(in A)')
        call cmd_dict%push('iares',         'integer angular resolution{10}')
        call cmd_dict%push('imgkind',       'type of image(ptcl|cavg|mic|movie){ptcl}')
        call cmd_dict%push('infile',        'file with inputs(.txt)')
        call cmd_dict%push('infile2',       'file with inputs(.txt)')
        call cmd_dict%push('interpfun',     'interpolation function(kb|linear){kb}')
        call cmd_dict%push('job_memory_per_task', 'memory in MB per task in distributed exec (typically memory per socket)')
        call cmd_dict%push('job_memory_per_task2D', 'memory in MB per 2D task in distributed exec (typically memory per socket)')
        call cmd_dict%push('jumpsz',        'size of contigous segment')
        call cmd_dict%push('keepvol',       'whether to keep iterative volumes')
        call cmd_dict%push('keys',          'keys of values to print')
        call cmd_dict%push('kv',            'acceleration voltage(in kV){300.}')
        call cmd_dict%push('label',         'discrete label(class|state){class}')
        call cmd_dict%push('lp',            'low-pass limit(in A)')
        call cmd_dict%push('lp2D',          'low-pass limit(in A)')
        call cmd_dict%push('lp_backgr',     'low-pass limit 4 solvent blurring(in A)')
        call cmd_dict%push('lp_ctf_estimate', 'low-pass limit 4 ctf_estimate(in A)')
        call cmd_dict%push('lp_iters',      '# of iterations for low-pass limited refinement')
        call cmd_dict%push('lp_pick',       'low-pass limit 4 picker(in A)')
        call cmd_dict%push('lplim_crit',    'corr criterion low-pass limit assignment(0.143-0.5){0.5}')
        call cmd_dict%push('lpstart',       'start low-pass limit(in A){15}')
        call cmd_dict%push('lpstop',        'stop low-pass limit(in A){8}')
        call cmd_dict%push('lpstop2D',      'stop low-pass limit(in A){8}')
        call cmd_dict%push('lpthres',      'resolution rejection limit(in A){30}')
        call cmd_dict%push('makemovie',     'produces files to generate a movie with ffmpeg(yes|no){no}')
        call cmd_dict%push('masscen',       'center using binarisation and mass centering(yes|no){no}')
        call cmd_dict%push('maxits',        'maximum # iterations')
        call cmd_dict%push('max_dose',      'maximum dose threshold (in e-/A2)')
        call cmd_dict%push('max_rad',       'particle longest  dim (in pixels)')
        call cmd_dict%push('mcpatch',       'Whether to perform Patch-based motion correction(yes|no){no}')
        call cmd_dict%push('mcpatch_thres','Whether to use the threshold for motion correction patch solution(yes|no){yes}')
        call cmd_dict%push('min_rad',       'particle shortest dim (in pixels)')
        call cmd_dict%push('minp',          'minimum cluster population')
        call cmd_dict%push('mirr',          'mirror(no|x|y){no}')
        call cmd_dict%push('moldiam',       'molecular diameter(in A)')
        call cmd_dict%push('mskdiam',       'mask diameter(in A)')
        call cmd_dict%push('mskfile',       'maskfile.ext')
        call cmd_dict%push('msktype',       'type of mask(hard|soft){soft}')
        call cmd_dict%push('mul',           'origin shift multiplication factor{1}')
        call cmd_dict%push('mw',            'molecular weight(in kD)')
        call cmd_dict%push('msklist',       'table (text file) of mask volume files(.txt)')
        call cmd_dict%push('ncls',          '# clusters')
        call cmd_dict%push('ncls_start',    '# clusters required to start prime2D streaming')
        call cmd_dict%push('ncunits',       '# computing units, can be < nparts{nparts}')
        call cmd_dict%push('nchunks',       '# of chunks')
        call cmd_dict%push('ndiscrete',     '# discrete orientations')
        call cmd_dict%push('ndev',          '# deviations in one-cluster clustering{2.0}')
        call cmd_dict%push('ndev2D',        '# deviations for 2D classes selection/rejection{1.5}')
        call cmd_dict%push('ndocs',         '# documents')
        call cmd_dict%push('needs_sigma',   'whether to calculate sigma2 during 3Drefinement')
        call cmd_dict%push('neg',           'invert contrast of images(yes|no)')
        call cmd_dict%push('newbox',        'new box for scaling (by Fourier padding/clipping')
        call cmd_dict%push('nframes',       '# frames{30}')
        call cmd_dict%push('nframesgrp',    '# frames to group before motion_correct(Falcon 3){0}')
        call cmd_dict%push('noise',         'noise initialisation(yes|no){no}')
        call cmd_dict%push('noise_norm',    'normalise based on sdev of background(yes|no){no}')
        call cmd_dict%push('norm',          'do statistical normalisation avg')
        call cmd_dict%push('nparts',        '# partitions in distributed exection')
        call cmd_dict%push('nparts_chunk',  '# partitions for chunks distributed exection')
        call cmd_dict%push('nparts_pool',   '# partitions for pool distributed exection')
        call cmd_dict%push('npix',          '# pixles/voxels in binary representation')
        call cmd_dict%push('nptcls',        '# images in stk/# orientations in oritab')
        call cmd_dict%push('nptcls_per_cls','# images per class for 2D streaming')
        call cmd_dict%push('nran',          '# random images to select')
        call cmd_dict%push('nrefs',         '# references used for picking{100}')
        call cmd_dict%push('nrepeats',      '# times to restart workflow{1}')
        call cmd_dict%push('nsig',          '# sigmas')
        call cmd_dict%push('nspace',        '# projection directions')
        call cmd_dict%push('nstates',       '# states to reconstruct')
        call cmd_dict%push('nsub',          '# proj dirs in coarse grid search{300}')
        call cmd_dict%push('nthr',          '# OpenMP threads{1}')
        call cmd_dict%push('nthr2D',        '# OpenMP threads for 2D classification{1}')
        call cmd_dict%push('numlen',        'length of number string')
        call cmd_dict%push('nvox',          '# voxels{0}')
        call cmd_dict%push('nxpatch',       'Motion correction # of patches along x-axis{5}')
        call cmd_dict%push('nypatch',       'Motion correction # of patches along y-axis{5}')
        call cmd_dict%push('objfun',        'objective function(cc|euclid){cc}')
        call cmd_dict%push('offset',        'pixels offset{10}')
        call cmd_dict%push('opt',           'optimiser (bfgs|simplex){bfgs}')
        call cmd_dict%push('order',         'order ptcls according to correlation(yes|no){no}')
        call cmd_dict%push('oritab',        'table of orientations(.txt|.simple)')
        call cmd_dict%push('oritab2',       'table of orientations 2(.txt|.simple)')
        call cmd_dict%push('oritab3D',      'table of 3D orientations(.txt|.simple)')
        call cmd_dict%push('oritype',       'SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D|projinfo|jobproc|compenv)')
        call cmd_dict%push('outer',         'outer mask radius(in pixels)')
        call cmd_dict%push('outfile',       'output document')
        call cmd_dict%push('outside',       'extract boxes outside the micrograph boundaries(yes|no){no}')
        call cmd_dict%push('outstk',        'output image stack')
        call cmd_dict%push('outvol',        'output volume{outvol.ext}')
        call cmd_dict%push('part',          'partition in distributed execution{1}')
        call cmd_dict%push('pcontrast',     'particle contrast(black|white){black}')
        call cmd_dict%push('pdbfile',       'input PDB formatted file')
        call cmd_dict%push('pgrp',          'point-group symmetry(cn|dn|t|o|i)')
        call cmd_dict%push('phaseplate',    'images obtained with Volta phaseplate(yes|no){no}')
        call cmd_dict%push('phrand',        'phase randomize(yes|no){no}')
        call cmd_dict%push('phshiftunit',   'additional phase-shift unit (radians|degrees){radians}')
        call cmd_dict%push('plaintexttab',  'plain text file of input parameters')
        call cmd_dict%push('plot',          'make plot(yes|no){no}')
        call cmd_dict%push('prg',           'SIMPLE program to execute')
        call cmd_dict%push('projfile',      'SIMPLE *.simple project file')
        call cmd_dict%push('projfile_target', 'another SIMPLE *.simple project file')
        call cmd_dict%push('projname',      'Project name (for creation of projname.simple)')
        call cmd_dict%push('prune',         'Whether to perform particles pruning')
        call cmd_dict%push('pspecsz',       'size of power spectrum(in pixels)')
        call cmd_dict%push('qsys_partition', 'Name of target partition of distributed computer system (SLURM/PBS)')
        call cmd_dict%push('qsys_partition2D', 'Name of target partition of distributed computer system dedictaed to 2D classification(SLURM/PBS)')
        call cmd_dict%push('qsys_qos',      'job scheduling priority (SLURM/PBS)')
        call cmd_dict%push('qsys_reservation', 'Name of reserved target partition of distributed computer system (SLURM/PBS)')
        call cmd_dict%push('real_filter',   'real-space filter kind(median|average|bman)')
        call cmd_dict%push('recvol_sigma',  'noise(sigma)-weighted volume reconstruction strategy(yes|no){no}')
        call cmd_dict%push('refine',        'refinement mode(shc|neigh|cluster|clustersym){shc}')
        call cmd_dict%push('refs',          'initial2Dreferences.ext')
        call cmd_dict%push('remap_cls',     'remove empty and renumber and/or expand # clusters(yes|no){no}')
        call cmd_dict%push('rnd',           'random(yes|no){no}')
        call cmd_dict%push('rrate',         'randomization rate{0.8}')
        call cmd_dict%push('scale',         'image scale factor{1}')
        call cmd_dict%push('shalgn',        'do 2D shift alignment(yes|no){no}')
        call cmd_dict%push('shbarrier',     'use shift search barrier constraint(yes|no){yes}')
        call cmd_dict%push('sherr',         'shift error(in pixels){2}')
        call cmd_dict%push('single',        'simulate a single image(yes|no){no}')
        call cmd_dict%push('smpd',          'sampling distance, same as EMANs apix(in A)')
        call cmd_dict%push('snr',           'signal-to-noise ratio')
        call cmd_dict%push('ptclw',         'soft particle weights(yes|no){yes}')
        call cmd_dict%push('speckind',      'power spectrum kind(real|power|sqrt|log|phase){sqrt}')
        call cmd_dict%push('startit',       'start iterating from here')
        call cmd_dict%push('starfile',      'STAR-formatted project file')
        call cmd_dict%push('state',         'state to extract')
        call cmd_dict%push('state2split',   'state group to split')
        call cmd_dict%push('stats',         'provide statistics(yes|no){yes}')
        call cmd_dict%push('stepsz',        'size of step{0}')
        call cmd_dict%push('stk',           'particle stack with all images(ptcls.ext)')
        call cmd_dict%push('stktab',        'list of per-micrograph stacks')
        call cmd_dict%push('stk2',          'stack 2 (in selection map: selected(cavgs).ext)')
        call cmd_dict%push('stk3',          'stack 3 (in selection map: (cavgs)2selectfrom.ext)')
        call cmd_dict%push('stk_backgr',    'stack with image for background subtraction')
        call cmd_dict%push('stream',        'sream (real time) execution mode(yes|no){no}')
        call cmd_dict%push('symrnd',        'randomize over symmetry operations(yes|no){no}')
        call cmd_dict%push('szsn',          'size of stochastic neighborhood{5}')
        call cmd_dict%push('taper_edges',   'to taper edges(yes|no){no}')
        call cmd_dict%push('thres',         'threshold (binarisation: 0-1; distance filer: # pixels; post-proc: pix val)')
        call cmd_dict%push('thres_low',     'lower threshold for canny edge detection')
        call cmd_dict%push('thres_up',      'upper threshold for canny edge detection')
        call cmd_dict%push('top',           'stop particle index')
        call cmd_dict%push('total_dose',    'Total movie dose (in e/A2)')
        call cmd_dict%push('trs',           'maximum halfwidth shift(in pixels)')
        call cmd_dict%push('trsstats',      'provide origin shift statistics(yes|no){no}')
        call cmd_dict%push('tseries',       'images represent a time-series(yes|no){no}')
        call cmd_dict%push('unidoc',        'unified resources and orientations doc')
        call cmd_dict%push('update_frac',   'fraction of particles to update(0.-1.){1.}')
        call cmd_dict%push('user_account',  'Name of user account in distributed computer system (SLURM/PBS)')
        call cmd_dict%push('user_email',    'Your e-mail address')
        call cmd_dict%push('user_project',  'Name of project in distributed computer system (SLURM/PBS)')
        call cmd_dict%push('verbose',       'verbosity flag (yes|no){no}')
        call cmd_dict%push('vis',           'visualise(yes|no)')
        call cmd_dict%push('vol1',          'input volume no1(invol1.ext)')
        call cmd_dict%push('vol2',          'input volume no2(invol2.ext)')
        call cmd_dict%push('vol_filt',      'input filter volume(vol_filt.ext)')
        call cmd_dict%push('vollist',       'table (text file) of volume files(.txt)')
        call cmd_dict%push('voltab',        'table (text file) of volume files(.txt)')
        call cmd_dict%push('voltab2',       'table 2(text file) of volume files(.txt)')
        call cmd_dict%push('walltime',      'Walltime in seconds for workload management systems{86340}')
        call cmd_dict%push('wcrit',         'Correlation to weights conversion scheme(softmax|zscore|sum|cen|exp|inv|no){softmax}')
        call cmd_dict%push('which_iter',    'iteration nr')
        call cmd_dict%push('width',         'falloff of inner mask or filter(in pixels){10}')
        call cmd_dict%push('wiener',        'Wiener restoration mode(yes|no){no}')
        call cmd_dict%push('winsz',         'half-width of window for real-space filter(in pixels)')
        call cmd_dict%push('wtype',         'type of orientation weights (factorial|flat){factorial}')
        call cmd_dict%push('xcoord',        'x coordinate{0}')
        call cmd_dict%push('xdim',          'x dimension(in pixles)')
        call cmd_dict%push('xsh',           'x shift(in pixels){0}')
        call cmd_dict%push('ycoord',        'y coordinate{0}')
        call cmd_dict%push('ydim',          'y dimension(in pixles)')
        call cmd_dict%push('ysh',           'y shift(in pixels){0}')
        call cmd_dict%push('zero',          'zeroing(yes|no){no}')
        call cmd_dict%push('zsh',           'z shift(in pixels){0}')
        cmd_dict_initialised = .true.
    end subroutine init_cmd_dict

    subroutine new_private_prgs
        private_prgs(:)%nreq = 0
        private_prgs(:)%nopt = 0

        ! AUTOMASK, for volumetric envelope masking
        call private_prgs(1)%set_name('automask')
        ! required keys
        call private_prgs(1)%push_req_key('mskdiam')
        call private_prgs(1)%push_req_key('amsklp')
        call private_prgs(1)%push_req_key('mw')
        call private_prgs(1)%push_req_key('thres')
        call private_prgs(1)%push_req_key('vol1')
        call private_prgs(1)%push_req_key('smpd')
        ! optional keys
        call private_prgs(1)%push_opt_key('edge')
        call private_prgs(1)%push_opt_key('binwidth')
        call private_prgs(1)%push_opt_key('nthr')

        ! binarize, for binarisation of stacks and volumes
        call private_prgs(2)%set_name('binarize')
        ! optional keys
        call private_prgs(2)%push_opt_key('nthr')
        call private_prgs(2)%push_opt_key('stk')
        call private_prgs(2)%push_opt_key('vol1')
        call private_prgs(2)%push_opt_key('thres')
        call private_prgs(2)%push_opt_key('npix')
        call private_prgs(2)%push_opt_key('grow')
        call private_prgs(2)%push_opt_key('edge')
        call private_prgs(2)%push_opt_key('neg')
        call private_prgs(2)%push_opt_key('outvol')
        call private_prgs(2)%push_opt_key('outstk')
        call private_prgs(2)%push_opt_key('fill_holes')

        ! CAVGASSEMBLE, for assembling class averages
        call private_prgs(3)%set_name('cavgassemble')
        ! required keys
        call private_prgs(3)%push_req_key('projfile')
        call private_prgs(3)%push_req_key('nparts')
        call private_prgs(3)%push_req_key('ncls')
        ! optional keys
        call private_prgs(3)%push_opt_key('nthr')
        call private_prgs(3)%push_opt_key('refs')

        ! CHECK_2DCONV, for convergence checking and run-time stats printing (2D)
        call private_prgs(4)%set_name('check_2Dconv')
        ! required keys
        call private_prgs(4)%push_req_key('projfile')

        ! CHECK_3DCONV, for convergence checking and run-time stats printing (3D)
        call private_prgs(5)%set_name('check_3Dconv')
        ! required keys
        call private_prgs(5)%push_req_key('projfile')
        call private_prgs(5)%push_req_key('pgrp')
        ! optional keys
        call private_prgs(5)%push_opt_key('lp')
        call private_prgs(5)%push_opt_key('nstates')
        call private_prgs(5)%push_opt_key('nspace')
        call private_prgs(5)%push_opt_key('refine')

        ! CHECK_BOX
        call private_prgs(6)%set_name('check_box')
        ! optional keys
        call private_prgs(6)%push_opt_key('stk')
        call private_prgs(6)%push_opt_key('vol1')

        ! CHECK_NPTCLS
        call private_prgs(7)%set_name('check_nptcls')
        ! required keys
        call private_prgs(7)%push_req_key('stk')

        ! EDGE_DETECT
        call private_prgs(8)%set_name('edge_detect')
        ! required keys
        call private_prgs(8)%push_req_key('detector')
        call private_prgs(8)%push_req_key('stk')
        call private_prgs(8)%push_req_key('automatic')
        call private_prgs(8)%push_req_key('smpd')
        ! optional keys
        call private_prgs(8)%push_opt_key('outstk')
        call private_prgs(8)%push_opt_key('thres')
        call private_prgs(8)%push_opt_key('bw_ratio')
        call private_prgs(8)%push_opt_key('thres_low')
        call private_prgs(8)%push_opt_key('thres_up')
        call private_prgs(8)%push_opt_key('lp')

        ! EXPORT_CAVGS
        call private_prgs(9)%set_name('export_cavgs')
        ! required keys
        call private_prgs(9)%push_req_key('projfile')
        ! optional keys
        call private_prgs(9)%push_opt_key('outstk')

        ! KSTEST, Kolmogorov-Smirnov test to deduce equivalence or
        ! non-equivalence between two distributions in a non-parametric manner
        call private_prgs(10)%set_name('kstest')
        ! required keys
        call private_prgs(10)%push_req_key('infile')
        call private_prgs(10)%push_req_key('infile2')

        ! MAKE_PICKREFS, for preparing templates for particle picking
        call private_prgs(11)%set_name('make_pickrefs')
        ! optional keys
        call private_prgs(11)%push_opt_key('nthr')
        call private_prgs(11)%push_opt_key('refs')
        call private_prgs(11)%push_opt_key('vol1')
        call private_prgs(11)%push_opt_key('pcontrast')
        call private_prgs(11)%push_opt_key('pgrp')

        ! MAP_CAVGS_SELECTION, for mapping class average selection to project
        call private_prgs(12)%set_name('map_cavgs_selection')
        ! required keys
        call private_prgs(12)%push_req_key('stk')
        call private_prgs(12)%push_req_key('stk2')
        call private_prgs(12)%push_req_key('projfile')
        call private_prgs(12)%push_req_key('prune')

        ! MASSCEN, for centering images acccording to their centre of mass
        call private_prgs(13)%set_name('masscen')
        ! required keys
        call private_prgs(13)%push_req_key('stk')
        call private_prgs(13)%push_req_key('smpd')
        call private_prgs(13)%push_req_key('lp')
        ! optional keys
        call private_prgs(13)%push_opt_key('mskdiam')
        call private_prgs(13)%push_opt_key('neg')
        call private_prgs(13)%push_opt_key('outstk')

        ! PICK_EXTRACT, for template-based particle picking & extraction
        call private_prgs(14)%set_name('pick_extract')
        ! required keys
        call private_prgs(14)%push_req_key('projfile')
        call private_prgs(14)%push_req_key('refs')
        ! optional keys
        call private_prgs(14)%push_opt_key('nthr')
        call private_prgs(14)%push_opt_key('lp')
        call private_prgs(14)%push_opt_key('thres')
        call private_prgs(14)%push_opt_key('ndev')
        call private_prgs(14)%push_opt_key('box_extract')
        call private_prgs(14)%push_opt_key('pcontrast')
        call private_prgs(14)%push_opt_key('outside')

        ! PRINT_PROJECT_VALS, for printing specific values in a project file field
        call private_prgs(15)%set_name('print_project_vals')
        ! required keys
        call private_prgs(15)%push_req_key('projfile')
        call private_prgs(15)%push_req_key('keys')
        call private_prgs(15)%push_req_key('oritype')

        ! RANK_CAVGS, for ranking class averages
        call private_prgs(16)%set_name('rank_cavgs')
        ! required keys
        call private_prgs(16)%push_req_key('projfile')
        call private_prgs(16)%push_req_key('stk')
        ! set optional keys
        call private_prgs(16)%push_opt_key('outstk')

        ! ROTMATS2ORIS, for converting a text file (9 records per line) describing rotation matrices into a SIMPLE oritab
        call private_prgs(17)%set_name('rotmats2oris')
        ! required keys
        call private_prgs(17)%push_req_key('infile')
        ! optional keys
        call private_prgs(17)%push_opt_key('outfile')
        call private_prgs(17)%push_opt_key('oritype')

        ! SPLIT, for splitting of image stacks into partitions for parallel execution
        call private_prgs(18)%set_name('split')
        ! required keys
        call private_prgs(18)%push_req_key('stk')
        call private_prgs(18)%push_req_key('smpd')
        call private_prgs(18)%push_req_key('nparts')

        ! VOLASSEMBLE, for asssembling subvolumes generated in distributed execution
        call private_prgs(19)%set_name('volassemble')
        ! required keys
        call private_prgs(19)%push_req_key('nparts')
        call private_prgs(19)%push_req_key('projfile')
        call private_prgs(19)%push_req_key('mskdiam')
        ! optional keys
        call private_prgs(19)%push_opt_key('nthr')
        call private_prgs(19)%push_opt_key('state')
        call private_prgs(19)%push_opt_key('nstates')
        call private_prgs(19)%push_opt_key('mskfile')

        ! CALC_PSPEC, for asssembling power spectra for refine3D
        call private_prgs(20)%set_name('calc_pspec')
        ! required keys
        call private_prgs(20)%push_req_key('nparts')
        call private_prgs(20)%push_req_key('projfile')
        call private_prgs(20)%push_req_key('nthr')

        ! TSERIES_MOTION_CORRECT
        call private_prgs(21)%set_name('tseries_motion_correct')
        ! required keys
        call private_prgs(21)%push_req_key('nthr')
        call private_prgs(21)%push_req_key('projfile')
        call private_prgs(21)%push_req_key('fromp')
        call private_prgs(21)%push_req_key('top')
        call private_prgs(21)%push_req_key('part')
        call private_prgs(21)%push_req_key('nparts')
        ! optional keys
        call private_prgs(21)%push_opt_key('nframesgrp')
        call private_prgs(21)%push_opt_key('mcpatch')
        call private_prgs(21)%push_opt_key('nxpatch')
        call private_prgs(21)%push_opt_key('nypatch')
        call private_prgs(21)%push_opt_key('trs')
        call private_prgs(21)%push_opt_key('lpstart')
        call private_prgs(21)%push_opt_key('lpstop')
        call private_prgs(21)%push_opt_key('bfac')
        call private_prgs(21)%push_opt_key('groupframes')
        call private_prgs(21)%push_opt_key('wcrit')

        ! tseries_track_particles
        call private_prgs(22)%set_name('tseries_track_particles')
        ! required keys
        call private_prgs(22)%push_req_key('nthr')
        call private_prgs(22)%push_req_key('projfile')
        call private_prgs(22)%push_req_key('fromp')
        call private_prgs(22)%push_req_key('top')
        call private_prgs(22)%push_req_key('part')
        call private_prgs(22)%push_req_key('nparts')
        call private_prgs(22)%push_req_key('box')
        call private_prgs(22)%push_req_key('xcoord')
        call private_prgs(22)%push_req_key('ycoord')
        call private_prgs(22)%push_req_key('ind')
        call private_prgs(22)%push_req_key('numlen')
        call private_prgs(22)%push_req_key('nframesgrp')
        ! optional keys
        call private_prgs(22)%push_opt_key('offset')
        call private_prgs(22)%push_opt_key('hp')
        call private_prgs(22)%push_opt_key('lp')
        call private_prgs(22)%push_opt_key('cenlp')
        call private_prgs(22)%push_opt_key('neg')
        call private_prgs(22)%push_opt_key('filter')

        ! PSPEC_INT_RANK
        call private_prgs(23)%set_name('pspec_int_rank')
        ! required keys
        call private_prgs(23)%push_req_key('nthr')
        call private_prgs(23)%push_req_key('stk')
        call private_prgs(23)%push_req_key('smpd')
        call private_prgs(23)%push_req_key('moldiam')
        ! optional keys
        call private_prgs(23)%push_opt_key('lp_backgr')

        ! CALC_GROUP_SIGMAS, for asssembling sigmas for refine3D
        call private_prgs(24)%set_name('calc_group_sigmas')
        ! required keys
        call private_prgs(24)%push_req_key('nparts')
        call private_prgs(24)%push_req_key('projfile')
        call private_prgs(24)%push_req_key('nthr')
        call private_prgs(24)%push_req_key('which_iter')

        n_private_prgs = 24
    end subroutine new_private_prgs

end module simple_private_prgs
