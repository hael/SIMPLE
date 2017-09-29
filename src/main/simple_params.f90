! provides global distribution of constants and derived constants
module simple_params
#include "simple_lib.f08"
use simple_ori,            only: ori
use simple_cmdline,        only: cmdline
use simple_magic_boxes,    only: find_magic_box
use simple_imghead,        only: find_ldim_nptcls
use simple_binoris,        only: binoris
use simple_stktab_handler, only: stktab_handler
implicit none

public :: params
private
#include "simple_local_flags.inc"

!> global parameters
type :: params
    ! global objects
    type(ori)             :: ori_glob
    type(ctfplan)         :: tfplan
    type(stktab_handler)  :: stkhandle
    ! yes/no decision variables in ascending alphabetical order
    character(len=3)      :: acf='no'             !< calculate autocorrelation function(yes|no){no}
    character(len=3)      :: append='no'          !< append in context of files(yes|no){no}
    character(len=3)      :: async='no'           !< asynchronous (yes|no){no}
    character(len=3)      :: autoscale='no'       !< automatic down-scaling(yes|no){yes}
    character(len=3)      :: avg='no'             !< calc average automatic (yes|no){no}
    character(len=3)      :: bin='no'             !< binarise image(yes|no){no}
    character(len=3)      :: center='yes'         !< center image(s)/class average(s)/volume(s)(yes|no){no}
    character(len=3)      :: classtats='no'       !< calculate class population statistics(yes|no){no}
    character(len=3)      :: clustvalid='no'      !< validate clustering(yes|homo|no){no}
    character(len=3)      :: compare='no'         !< do comparison(yes|no){no}
    character(len=3)      :: countvox='no'        !< count # voxels(yes|no){no}
    character(len=3)      :: ctfstats='no'        !< calculate ctf statistics(yes|no){no}
    character(len=3)      :: cure='no'
    character(len=3)      :: dev='no'             !< development flag for experimental code(yes|no){no}
    character(len=3)      :: discrete='no'        !< be discrete(yes|no){no}
    character(len=3)      :: diverse='no'         !< diverse or not flag (yes|no){no}
    character(len=3)      :: doalign='yes'
    character(len=3)      :: dopca='yes'
    character(len=3)      :: dopick='yes'
    character(len=3)      :: doprint='no'
    character(len=3)      :: dynlp='yes'          !< automatic resolution limit update(yes|no){yes}
    character(len=3)      :: errify='no'          !< introduce error(yes|no){no}
    character(len=3)      :: even='no'            !< even orientation distribution(yes|no){no}
    character(len=3)      :: ft2img='no'          !< convert Fourier transform to real image of power(yes|no){no}
    character(len=3)      :: guinier='no'         !< calculate Guinier plot(yes|no){no}
    character(len=3)      :: kmeans='yes'
    character(len=3)      :: local='no'
    character(len=3)      :: masscen='no'         !< center using binarisation and mass centering(yes|no){no}
    character(len=3)      :: match_filt='yes'     !< matched filter on (yes|no){yes}
    character(len=3)      :: merge='no'
    character(len=3)      :: mirr='no'            !< mirror(no|x|y){no}
    character(len=3)      :: neg='no'             !< invert contrast of images(yes|no)
    character(len=3)      :: noise_norm ='no'
    character(len=3)      :: noise='no'           !< noise initialisation(yes|no){no}
    character(len=3)      :: norec='no'           !< do not reconstruct volume(s)(yes|no){no}
    character(len=3)      :: norm='no'            !< do statistical normalisation avg
    character(len=3)      :: order='no'           !< order ptcls according to correlation(yes|no){no}
    character(len=3)      :: outside='no'         !< extract boxes outside the micrograph boundaries(yes|no){no}
    character(len=3)      :: pad='no'
    character(len=3)      :: pgrp_known='no'      !< point-group known a priori(yes|no){no}
    character(len=3)      :: phaseplate='no'      !< images obtained with phaseplate(yes|no){no}
    character(len=3)      :: phrand='no'          !< phase randomize(yes|no){no}
    character(len=3)      :: plot='no'            !< make plot(yes|no){no}
    character(len=3)      :: pproc='yes'          !< whether to perform volume post-processing(yes|no){yes}
    character(len=3)      :: projstats='no'
    character(len=3)      :: readwrite='no'
    character(len=3)      :: remap_classes='no'
    character(len=3)      :: restart='no'
    character(len=3)      :: rnd='no'             !< random(yes|no){no}
    character(len=3)      :: rm_outliers='yes'    !< remove outliers{yes}
    character(len=3)      :: roalgn='no'
    character(len=3)      :: round='no'
    character(len=3)      :: shalgn='no'          !< do 2D shift alignment(yes|no){no}
    character(len=3)      :: shellnorm='no'
    character(len=3)      :: shbarrier='yes'      !< use shift search barrier constraint(yes|no){yes}
    character(len=3)      :: single='no'          !< simulate a single image(yes|no){no}
    character(len=3)      :: soften='no'          !< soften envelope with cosine edge(yes|no){no}
    character(len=3)      :: stats='no'           !< provide statistics(yes|no){yes}
    character(len=3)      :: stream='no'          !< sream (real time) execution mode(yes|no){no}
    character(len=3)      :: swap='no'
    character(len=3)      :: taper_edges='no'     !< self-explanatory
    character(len=3)      :: test='no'
    character(len=3)      :: tomo='no'            !< tomography mode(yes|no){no}
    character(len=3)      :: time='no'
    character(len=3)      :: trsstats='no'        !< provide origin shift statistics(yes|no){no}
    character(len=3)      :: tseries='no'         !< images represent a time-series(yes|no){no}
    character(len=3)      :: vis='no'             !< visualise(yes|no)
    character(len=3)      :: weights2D='no'
    character(len=3)      :: zero='no'            !< zeroing(yes|no){no}
    ! other fixed length character variables in ascending alphabetical order
    character(len=STDLEN) :: angastunit='degrees' !< angle of astigmatism unit (radians|degrees){degrees}
    character(len=4)      :: automsk='no'
    character(len=STDLEN) :: boxfile=''           !< file with EMAN particle coordinates(.txt)
    character(len=STDLEN) :: boxtab=''            !< table (text file) of files with EMAN particle coordinates(.txt)
    character(len=STDLEN) :: boxtype='eman'
    character(len=STDLEN) :: chunktag=''
    character(len=STDLEN) :: classdoc=''          !< doc with per-class stats(.txt)
    character(len=STDLEN) :: comlindoc=''         !< shc_clustering_nclsX.txt
    character(len=STDLEN) :: ctf='no'             !< ctf flag(yes|no|flip)
    character(len=STDLEN) :: cwd=''
    character(len=STDLEN) :: deftab=''            !< file with CTF info(.txt|.bin)
    character(len=STDLEN) :: dfunit='microns'     !< defocus unit (A|microns){microns}
    character(len=STDLEN) :: dir=''               !< directory
    character(len=STDLEN) :: dir_movies=''        !< grab mrc mrcs files from here
    character(len=STDLEN) :: dir_reject='rejected'!< move rejected files to here{rejected}
    character(len=STDLEN) :: dir_select='selected'!< move selected files to here{selected}
    character(len=STDLEN) :: dir_target=''        !< put output here
    character(len=STDLEN) :: dir_ptcls=''
    character(len=STDLEN) :: dir_mics=''          !< grab micrographs from here
    character(len=STDLEN) :: dockmode='eul'       !< volume docking mode(eul|shift|eulshift|all){eul}
    character(len=STDLEN) :: doclist=''           !< list of oritabs for different states
    character(len=STDLEN) :: eo='yes'             !< use FSC for filtering and low-pass limit update(yes|aniso|no){no}
    character(len=STDLEN) :: exec_abspath=''
    character(len=STDLEN) :: exp_doc=''           !< specifying exp_time and dose_rate per tomogram
    character(len=4)      :: ext='.mrc'           !< file extension{.mrc}
    character(len=4)      :: ext_meta=''          !< meta data file extension(.txt|.bin)
    character(len=STDLEN) :: extrmode='all'
    character(len=STDLEN) :: fbody=''             !< file body
    character(len=STDLEN) :: featstk='expecstk.bin'
    character(len=STDLEN) :: filetab=''           !< list of files(.txt)
    character(len=STDLEN) :: fname=''             !< file name
    character(len=STDLEN) :: frcs=''              !< binary file with per-class/proj Fourier Ring Correlations(.bin) 
    character(len=STDLEN) :: fsc='fsc_state01.bin'!< binary file with FSC info{fsc_state01.bin}
    character(len=STDLEN) :: hfun='sigm'          !< function used for normalization(sigm|tanh|lin){sigm}
    character(len=STDLEN) :: hist='corr'          !< give variable for histogram plot
    character(len=STDLEN) :: infile=''            !< file with inputs(.txt|.bin)
    character(len=STDLEN) :: label='class'        !< discrete label(class|state){class}
    character(len=STDLEN) :: mskfile=''           !< maskfile.ext
    character(len=STDLEN) :: msktype='soft'       !< type of mask(hard|soft){soft}
    character(len=STDLEN) :: opt='simplex'        !< optimiser (powell|simplex|oasis|bforce|pso|de){simplex}
    character(len=STDLEN) :: oritab=''            !< table  of orientations(.txt|.bin)
    character(len=STDLEN) :: oritab2=''           !< 2nd table of orientations(.txt|.bin)
    character(len=STDLEN) :: oritab3D=''          !< table of 3D orientations(.txt|.bin)
    character(len=STDLEN) :: outfile=''           !< output document
    character(len=STDLEN) :: outstk=''            !< output image stack
    character(len=STDLEN) :: outstk2=''           !< output image stack 2nd
    character(len=STDLEN) :: outvol=''            !< output volume{outvol.ext}
    character(len=STDLEN) :: ctffind_doc=''       !< per-micrograph CTF parameters to transfer
    character(len=STDLEN) :: pcastk='pcavecinstk.bin'
    character(len=STDLEN) :: pcontrast='black'    !< particle contrast(black|white){black}
    character(len=STDLEN) :: pdbfile=''           !< PDB file
    character(len=STDLEN) :: pdfile='pdfile.bin'
    character(len=STDLEN) :: pgrp='c1'            !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN) :: plaintexttab=''      !< plain text file of input parameters
    character(len=STDLEN) :: prg=''               !< SIMPLE program to execute
    character(len=STDLEN) :: real_filter=''
    character(len=STDLEN) :: refine='no'
    character(len=STDLEN) :: refs_msk=''
    character(len=STDLEN) :: refs=''              !< initial2Dreferences.ext
    character(len=STDLEN) :: speckind='sqrt'      !< power spectrum kind(amp|square|phase|real|log|sqrt){sqrt}
    character(len=STDLEN) :: split_mode='even'
    character(len=STDLEN) :: stk_part=''
    character(len=STDLEN) :: stk=''               !< particle stack with all images(ptcls.ext)
    character(len=STDLEN) :: stktab=''            !< list of per-micrograph stacks
    character(len=STDLEN) :: stk2=''              !< 2nd stack(in map2ptcls/select: selected(cavgs).ext)
    character(len=STDLEN) :: stk3=''              !< 3d stack (in map2ptcls/select: (cavgs)2selectfrom.ext)
    character(len=STDLEN) :: stk_backgr=''        !< stack with image for background subtraction
    character(len=STDLEN) :: tomoseries=''        !< filetable of filetables of tomograms
    character(len=STDLEN) :: unidoc=''            !< unified resources and orientations doc
    character(len=STDLEN) :: vol=''
    character(len=STDLEN) :: vol_filt=''          !< input filter volume(vol_filt.ext)
    character(len=STDLEN) :: vollist=''           !< table (text file) of volume files(.txt)
    character(len=STDLEN) :: vols(MAXS)=''
    character(len=STDLEN) :: voltab=''            !< table (text file) of volume files(.txt)
    character(len=STDLEN) :: voltab2=''           !< 2nd table (text file) of volume files(.txt)
    character(len=STDLEN) :: wfun='kb'
    ! integer variables in ascending alphabetical order
    integer :: astep=1
    integer :: avgsz=0
    integer :: batchsz=0
    integer :: balance=0           !< max pop for balancing restraint{0}
    integer :: binwidth=1          !< binary layers grown for molecular envelope(in pixels){1}
    integer :: box=0               !< square image size(in pixels)
    integer :: boxconvsz=256       !< size of box used for box-convolution(in pixels)
    integer :: boxmatch=0
    integer :: box_original
    integer :: boxpd=0
    integer :: chunk=0
    integer :: chunksz=0           !< # images/orientations in chunk
    integer :: class=1             !< cluster identity
    integer :: clip=0              !< clipped image box size(in pixels)
    integer :: corner=0            !< corner size(in pixels){0}
    integer :: cube=0              !< side size(in pixels){0}
    integer :: edge=6              !< edge size for softening molecular envelope(in pixels)
    integer :: extr_iter=1
    integer :: find=1              !< Fourier index
    integer :: nframesgrp=0        !< # frames to group before unblur(Falcon 3){0}
    integer :: fromf=1             !< start frame index
    integer :: fromp=1             !< start ptcl index
    integer :: fromm=1             !< start movie index
    integer :: fstep=1
    integer :: grow=0              !< # binary layers to grow(in pixels)
    integer :: iares=10            !< integer angular resolution{10}
    integer :: ind=0
    integer :: iptcl=1
    integer :: jptcl=1
    integer :: jumpsz=0            !< size of contigous segment
    integer :: kfromto(2)
    integer :: kstop_grid=0
    integer :: ldim(3)=0
    integer :: maxits=500          !< maximum # iterations
    integer :: maxp=0
    integer :: minp=10             !< minimum cluster population
    integer :: mrcmode=2
    integer :: navgs=1
    integer :: ncunits=0           !< # computing units, can be < nparts{nparts}
    integer :: nbest=10
    integer :: nboot=0
    integer :: ncls=500            !< # clusters
    integer :: ncls_start=10       !< minimum # clusters for 2D streaming
    integer :: ncomps=0
    integer :: ndiscrete=0         !< # discrete orientations
    integer :: ndocs=0             !< # documents
    integer :: newbox=0            !< new box for scaling (by Fourier padding/clipping
    integer :: newbox2=0
    integer :: nframes=0           !< # frames{30}
    integer :: nmembers=0
    integer :: nnn=50              !< # nearest neighbors{500}
    integer :: nmics=0             !< # micographs
    integer :: noris=0
    integer :: nparts=1            !< # partitions in distributed exection
    integer :: npeaks=1            !< # nonzero orientation weights{1}
    integer :: npix=0              !< # pixles/voxels in binary representation
    integer :: nptcls=1            !< # images in stk/# orientations in oritab
    integer :: nptcls_per_cls=400  !< # images in stk/# orientations in oritab
    integer :: nran=0              !< # random images to select
    integer :: nrefs=100           !< # references used for picking{100}
    integer :: nrestarts=1
    integer :: nrepeats=HETNREPEATS !< # repeats in het_ensemble worklow
    integer :: nrots=0
    integer :: nspace=1000         !< # projection directions
    integer :: nstates=1           !< # states to reconstruct
    integer :: nsub=300
    integer :: nsym=1
    integer :: nthr=1              !< # OpenMP threads{1}
    integer :: numlen=0            !< length of number string
    integer :: numlen_tomo=3       !< length of number string tomo series index{3}
    integer :: nvalid=0
    integer :: nvars=30
    integer :: nvox=0              !< # voxels{0}
    integer :: offset=10           !< pixels offset{10}
    integer :: part=1
    integer :: pcasz=0
    integer :: ppca=0
    integer :: pspecsz=512         !< size of power spectrum(in pixels)
    integer :: pspecsz_unblur=512  !< size of power spectrum 4 unblur(in pixels)
    integer :: pspecsz_ctffind=1024
    integer :: ptcl=1
    integer :: ring1=2
    integer :: ring2=0
    integer :: spec=0
    integer :: startit=1           !< start iterating from here
    integer :: state=1             !< state to extract
    integer :: state2split=0       !< state group to split
    integer :: stepsz=1            !< size of step{0}
    integer :: szsn=SZSN_INIT      !< size of stochastic neighborhood{5}
    integer :: tofny=0
    integer :: tof=1               !< stop frame index
    integer :: top=1
    integer :: tos=1
    integer :: trsstep=1
    integer :: update=1000
    integer :: which_iter=0        !< iteration nr
    integer :: xcoord=0            !< x coordinate{0}
    integer :: ycoord=0            !< y coordinate{0}
    integer :: xdim=0              !< x dimension(in pixles)
    integer :: xdimpd=0
    integer :: ydim=0              !< y dimension(in pixles)
    ! real variables in ascending alphabetical order
    real    :: alpha=2.
    real    :: amsklp=20.          !< low-pass limit for envelope mask generation(in A)
    real    :: angerr=0.           !< angular error(in degrees){0}
    real    :: ares=7.
    real    :: astigerr=0.         !< astigmatism error(in microns)
    real    :: astigtol=0.05       !< expected (tolerated) astigmatism(in microns){0.1}
    real    :: athres=0.           !< angular threshold(in degrees)
    real    :: batchfrac=1.0
    real    :: bfac=200            !< bfactor for sharpening/low-pass filtering(in A**2){200.}
    real    :: bfacerr=50.         !< bfactor error in simulated images(in A**2){0}
    real    :: cenlp=30.           !< low-pass limit for binarisation in centering(in A){30 A}
    real    :: cs=2.7              !< spherical aberration constant(in mm){2.7}
    real    :: ctfreslim=8.
    real    :: dcrit_rel=0.5       !< critical distance relative to box(0-1){0.5}
    real    :: deflim=4.
    real    :: defocus=3.          !< defocus(in microns){3.}
    real    :: dens=0.
    real    :: dfclose=1.
    real    :: dffar=4.
    real    :: dferr=1.            !< defocus error(in microns){1.0}
    real    :: dfmax=5.0           !< maximum expected defocus(in microns)
    real    :: dfmin=0.5           !< minimum expected defocus(in microns)
    real    :: dfsdev=0.1
    real    :: dfstep=0.05         !< defocus step size for grid search
    real    :: dose_rate=30.0      !< dose rate(in e/A2/s)
    real    :: dstep=0.
    real    :: dsteppd=0.
    real    :: e1=0.               !< 1st Euler(in degrees){0}
    real    :: e2=0.               !< 2nd Euler(in degrees){0}
    real    :: e3=0.               !< 3d Euler(in degrees){0}
    real    :: eps=0.003           !< learning rate{0.003}
    real    :: eullims(3,2)=0.
    real    :: exp_time=2.0        !< exposure time(in s)
    real    :: filwidth=0.         !< width of filament (in A)
    real    :: fny=0.
    real    :: frac=1.             !< fraction of ptcls(0-1){1}
    real    :: fraca=0.07          !< fraction of amplitude contrast used for fitting CTF{0.07}
    real    :: fracdeadhot=0.05    !< fraction of dead or hot pixels{0.01}
    real    :: frac_outliers=0.
    real    :: fraczero=0.
    real    :: ftol=1e-6
    real    :: gw=0.5
    real    :: hp=100.             !< high-pass limit(in A)
    real    :: hp_ctffind=30.      !< high-pass limit 4 ctffind(in A)
    real    :: inner=0.            !< inner mask radius(in pixels)
    real    :: kv=300.             !< acceleration voltage(in kV){300.}
    real    :: lam=0.5
    real    :: lp_dyn=20.
    real    :: lp_grid=20.
    real    :: lp=20.              !< low-pass limit(in A)
    real    :: lp_ctffind=5.0      !< low-pass limit 4 ctffind(in A)
    real    :: lp_pick=20.         !< low-pass limit 4 picker(in A)
    real    :: lplims2D(3)
    real    :: lpmed=20.
    real    :: lpstart=0.          !< start low-pass limit(in A){15}
    real    :: lpstop=7.0          !< stop low-pass limit(in A){8}
    real    :: lpvalid=20.
    real    :: moldiam=140.        !< molecular diameter(in A)
    real    :: moment=0.
    real    :: msk=0.              !< mask radius(in pixels)
    real    :: mul=1.              !< origin shift multiplication factor{1}
    real    :: mw=0.               !< molecular weight(in kD)
    real    :: ndev=2.0            !< # deviations in one-cluster clustering
    real    :: neigh=0.2
    real    :: nsig=2.5            !< # sigmas
    real    :: optlims(7,2)=0.
    real    :: outer=0.            !< outer mask radius(in pixels)
    real    :: phranlp=35.         !< low-pass phase randomize(yes|no){no}
    real    :: power=2.
    real    :: scale=1.            !< image scale factor{1}
    real    :: scale2=1.           !< image scale factor 2nd{1}
    real    :: sherr=0.            !< shift error(in pixels){2}
    real    :: smpd=2.             !< sampling distance, same as EMANs apix(in A)
    real    :: smpd_targets2D(2)
    real    :: snr=0.              !< signal-to-noise ratio
    real    :: thres=0.            !< threshold (binarisation: 0-1; distance filer: in pixels)
    real    :: time_per_image=200.
    real    :: time_per_frame=0.
    real    :: trs=0.              !< maximum halfwidth shift(in pixels)
    real    :: update_frac = 1.
    real    :: width=10.           !< falloff of inner mask(in pixels){10}
    real    :: winsz=1.
    real    :: xsh=0.              !< x shift(in pixels){0}
    real    :: ysh=0.              !< y shift(in pixels){0}
    real    :: zsh=0.              !< z shift(in pixels){0}
    ! logical variables in ascending alphabetical order
    logical :: cyclic(7)       = .false.
    logical :: l_distr_exec    = .false.
    logical :: l_chunk_distr   = .false.
    logical :: l_match_filt    = .true.
    logical :: doshift         = .false.
    logical :: l_envmsk        = .false.
    logical :: l_autoscale     = .false.
    logical :: l_dose_weight   = .false. 
    logical :: l_frac_update   = .false.
    logical :: l_innermsk      = .false. 
    logical :: l_pick          = .false.
    logical :: l_remap_classes = .false.
    logical :: l_stktab_input  = .false.
  contains
    procedure :: new
end type params

interface params
    module procedure constructor
end interface

contains

    !> \brief  is a constructor
    function constructor( cline, checkdistr, allow_mix ) result( self )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: checkdistr, allow_mix
        type(params) :: self
        call self%new( cline, checkdistr, allow_mix )
    end function constructor

    !> \brief  is a constructor
    subroutine new( self, cline, checkdistr, allow_mix )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: round2even
        use simple_map_reduce  ! use all in there
        class(params),     intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: checkdistr, allow_mix
        type(binoris)                    :: bos
        character(len=STDLEN)            :: cwd_local, debug_local, verbose_local
        character(len=STDLEN)            :: stk_part_fname_sc, stk_part_fname
        character(len=1)                 :: checkupfile(50)
        character(len=:), allocatable    :: conv
        integer                          :: i, ncls, ifoo, lfoo(3), cntfile
        logical                          :: nparts_set, vol_defined(MAXS), ccheckdistr, aamix
        nparts_set        = .false.
        vol_defined(MAXS) = .false.
        debug_local       = 'no'
        verbose_local     = 'no'
        ! take care of optionals
        ccheckdistr = .true.
        if( present(checkdistr) ) ccheckdistr = checkdistr
        aamix = .false.
        if( present(allow_mix) ) aamix = allow_mix
        ! file counter
        cntfile = 0
        ! make global ori
        call self%ori_glob%new
        ! get cwd
        call simple_getcwd(self%cwd)
        cwd_local = self%cwd
        ! get absolute path of executable
        call getarg(0,self%exec_abspath)
        ! take care of debug/verbose flags
        call check_carg('debug', debug_local)
        if( debug_local == 'yes' )then
            global_debug = .true. ! from simple_params
            debug        = .true. ! from simple_local_flags.inc
        end if
        call check_carg('verbose', verbose_local)
        if(verbose_local == 'yes')then
            global_verbose = .true.
            verbose        = .true.
        end if
        ! default initialisations that depend on meta-data file format
        self%outfile = 'outfile'//METADATEXT
        ! checkers in ascending alphabetical order
        call check_carg('acf',            self%acf)
        call check_carg('angastunit',     self%angastunit)
        call check_carg('append',         self%append)
        call check_carg('async',          self%async)
        call check_carg('automsk',        self%automsk)
        call check_carg('autoscale',      self%autoscale)
        call check_carg('avg',            self%avg)
        call check_carg('bin',            self%bin)
        call check_carg('boxtype',        self%boxtype)
        call check_carg('center',         self%center)
        call check_carg('chunktag',       self%chunktag)
        call check_carg('classtats',      self%classtats)
        call check_carg('clustvalid',     self%clustvalid)
        call check_carg('compare',        self%compare)
        call check_carg('countvox',       self%countvox)
        call check_carg('ctf',            self%ctf)
        call check_carg('ctfstats',       self%ctfstats)
        call check_carg('cure',           self%cure)
        call check_carg('dfunit',         self%dfunit)
        call check_carg('dir',            self%dir)
        call check_carg('dir_movies',     self%dir_movies)
        call check_carg('dir_ptcls',      self%dir_ptcls)
        call check_carg('dir_reject',     self%dir_reject)
        call check_carg('dir_select',     self%dir_select)
        call check_carg('dir_target',     self%dir_target)
        call check_carg('dir_mics',       self%dir_mics)
        call check_carg('discrete',       self%discrete)
        call check_carg('diverse',        self%diverse)
        call check_carg('doalign',        self%doalign)
        call check_carg('dockmode',       self%dockmode)
        call check_carg('dev',            self%dev)
        call check_carg('dopca',          self%dopca)
        call check_carg('dopick',         self%dopick)
        call check_carg('doprint',        self%doprint)
        call check_carg('dynlp',          self%dynlp)
        call check_carg('eo',             self%eo)
        call check_carg('errify',         self%errify)
        call check_carg('even',           self%even)
        call check_carg('exp_doc',        self%exp_doc)
        call check_carg('extrmode',       self%extrmode)
        call check_carg('fbody',          self%fbody)
        call check_carg('ft2img',         self%ft2img)
        call check_carg('guinier',        self%guinier)
        call check_carg('hfun',           self%hfun)
        call check_carg('hist',           self%hist)
        call check_carg('kmeans',         self%kmeans)
        call check_carg('label',          self%label)
        call check_carg('local',          self%local)
        call check_carg('masscen',        self%masscen)
        call check_carg('match_filt',     self%match_filt)
        call check_carg('merge',          self%merge)
        call check_carg('mirr',           self%mirr)
        call check_carg('msktype',        self%msktype)
        call check_carg('neg',            self%neg)
        call check_carg('noise_norm',     self%noise_norm)
        call check_carg('noise',          self%noise)
        call check_carg('norec',          self%norec)
        call check_carg('norm',           self%norm)
        call check_carg('opt',            self%opt)
        call check_carg('order',          self%order)
        call check_carg('outside',        self%outside)
        call check_carg('pad',            self%pad)
        call check_carg('pcontrast',      self%pcontrast)
        call check_carg('pgrp',           self%pgrp)
        call check_carg('pgrp_known',     self%pgrp_known)
        call check_carg('phaseplate',     self%phaseplate)
        call check_carg('phrand',         self%phrand)
        call check_carg('prg',            self%prg)
        call check_carg('projstats',      self%projstats)
        call check_carg('readwrite',      self%readwrite)
        call check_carg('real_filter',    self%real_filter)
        call check_carg('refine',         self%refine)
        call check_carg('refs',           self%refs)
        call check_carg('remap_classes',  self%remap_classes)
        call check_carg('restart',        self%restart)
        call check_carg('rm_outliers',    self%rm_outliers)
        call check_carg('rnd',            self%rnd)
        call check_carg('roalgn',         self%roalgn)
        call check_carg('round',          self%round)
        call check_carg('shalgn',         self%shalgn)
        call check_carg('shbarrier',      self%shbarrier)
        call check_carg('shellnorm',      self%shellnorm)
        call check_carg('single',         self%single)
        call check_carg('soften',         self%soften)
        call check_carg('speckind',       self%speckind)
        call check_carg('stats',          self%stats)
        call check_carg('stream',         self%stream)
        call check_carg('swap',           self%swap)
        call check_carg('taper_edges',    self%taper_edges)
        call check_carg('test',           self%test)
        call check_carg('time',           self%time)
        call check_carg('tomo',           self%tomo)
        call check_carg('tomoseries',     self%tomoseries)
        call check_carg('trsstats',       self%trsstats)
        call check_carg('tseries',        self%tseries)
        call check_carg('vis',            self%vis)
        call check_carg('vol',            self%vol)
        call check_carg('wfun',           self%wfun)
        call check_carg('weights2D',      self%weights2D)
        call check_carg('zero',           self%zero)
        ! File args
        call check_file('boxfile',        self%boxfile,      'T')
        call check_file('boxtab',         self%boxtab,       'T')
        call check_file('classdoc',       self%classdoc,     'T')
        call check_file('ctffind_doc',    self%ctffind_doc,  'T', 'B')
        call check_file('comlindoc',      self%comlindoc,    'T')
        call check_file('deftab',         self%deftab,       'T', 'B')
        call check_file('doclist',        self%doclist,      'T')
        call check_file('ext',            self%ext,          notAllowed='T')
        call check_file('ext_meta',       self%ext_meta,     'T', 'B')
        call check_file('filetab',        self%filetab,      'T')
        call check_file('fname',          self%fname)
        call check_file('frcs',           self%frcs,         'B')
        call check_file('fsc',            self%fsc,          'B')
        call check_file('infile',         self%infile)
        call check_file('mskfile',        self%mskfile,      notAllowed='T')
        call check_file('oritab',         self%oritab,       'T', 'B')
        call check_file('oritab2',        self%oritab2,      'T', 'B')
        call check_file('oritab3D',       self%oritab3D,     'T', 'B')
        call check_file('outfile',        self%outfile,      'T', 'B')
        call check_file('outstk',         self%outstk,       notAllowed='T')
        call check_file('outstk2',        self%outstk2,      notAllowed='T')
        call check_file('outvol',         self%outvol,       notAllowed='T')
        call check_file('pdbfile',        self%pdbfile)
        call check_file('plaintexttab',   self%plaintexttab, 'T')
        call check_file('stk',            self%stk,          notAllowed='T')
        call check_file('stktab',         self%stktab,       'T')
        call check_file('stk2',           self%stk2,         notAllowed='T')
        call check_file('stk3',           self%stk3,         notAllowed='T')
        call check_file('stk_backgr',     self%stk_backgr,   notAllowed='T')
        call check_file('unidoc',         self%unidoc,       'T')
        call check_file('vol_filt',       self%vol_filt,     notAllowed='T')
        call check_file('vollist',        self%vollist,      'T')
        call check_file('voltab',         self%voltab,       'T')
        call check_file('voltab2',        self%voltab2,      'T')
        ! Integer args
        call check_iarg('astep',          self%astep)
        call check_iarg('avgsz',          self%avgsz)
        call check_iarg('balance',        self%balance)
        call check_iarg('binwidth',       self%binwidth)
        call check_iarg('box',            self%box)
        call check_iarg('boxconvsz',      self%boxconvsz)
        call check_iarg('chunk',          self%chunk)
        call check_iarg('chunksz',        self%chunksz)
        call check_iarg('clip',           self%clip)
        call check_iarg('corner',         self%corner)
        call check_iarg('cube',           self%cube)
        call check_iarg('edge',           self%edge)
        call check_iarg('extr_iter',      self%extr_iter)
        call check_iarg('find',           self%find)
        call check_iarg('nframesgrp',     self%nframesgrp)
        call check_iarg('fromf',          self%fromf)
        call check_iarg('fromp',          self%fromp)
        call check_iarg('fromm',          self%fromm)
        call check_iarg('fstep',          self%fstep)
        call check_iarg('grow',           self%grow)
        call check_iarg('iares',          self%iares)
        call check_iarg('ind',            self%ind)
        call check_iarg('jumpsz',         self%jumpsz)
        call check_iarg('maxits',         self%maxits)
        call check_iarg('maxp',           self%maxp)
        call check_iarg('minp',           self%minp)
        call check_iarg('mrcmode',        self%mrcmode)
        call check_iarg('navgs',          self%navgs)
        call check_iarg('nbest',          self%nbest)
        call check_iarg('nboot',          self%nboot)
        call check_iarg('ncls',           self%ncls)
        call check_iarg('ncls_start',     self%ncls_start)
        call check_iarg('ncunits',        self%ncunits)
        call check_iarg('ndiscrete',      self%ndiscrete)
        call check_iarg('ndocs',          self%ndocs)
        call check_iarg('newbox',         self%newbox)
        call check_iarg('nframes',        self%nframes)
        call check_iarg('nmembers',       self%nmembers)
        call check_iarg('nnn',            self%nnn)
        call check_iarg('noris',          self%noris)
        call check_iarg('nran',           self%nran)
        call check_iarg('nrefs',          self%nrefs)
        call check_iarg('nrepeats',       self%nrepeats)
        call check_iarg('nrestarts',      self%nrestarts)
        call check_iarg('nspace',         self%nspace)
        call check_iarg('nsub',           self%nsub)
        call check_iarg('nstates',        self%nstates)
        call check_iarg('class',          self%class)
        call check_iarg('nparts',         self%nparts)
        call check_iarg('npeaks',         self%npeaks)
        call check_iarg('npix',           self%npix)
        call check_iarg('nptcls',         self%nptcls)
        call check_iarg('nptcls_per_cls', self%nptcls_per_cls)
        call check_iarg('nthr',           self%nthr)
        call check_iarg('numlen',         self%numlen)
        call check_iarg('numlen_tomo',    self%numlen_tomo)
        call check_iarg('nvars',          self%nvars)
        call check_iarg('nvox',           self%nvox)
        call check_iarg('offset',         self%offset)
        call check_iarg('part',           self%part)
        call check_iarg('ppca',           self%ppca)
        call check_iarg('pspecsz',        self%pspecsz)
        call check_iarg('pspecsz_unblur', self%pspecsz_unblur)
        call check_iarg('pspecsz_ctffind', self%pspecsz_ctffind)
        call check_iarg('ring1',          self%ring1)
        call check_iarg('ring2',          self%ring2)
        call check_iarg('startit',        self%startit)
        call check_iarg('state',          self%state)
        call check_iarg('state2split',    self%state2split)
        call check_iarg('stepsz',         self%stepsz)
        call check_iarg('szsn',           self%szsn)
        call check_iarg('tof',            self%tof)
        call check_iarg('top',            self%top)
        call check_iarg('tos',            self%tos)
        call check_iarg('trsstep',        self%trsstep)
        call check_iarg('update',         self%update)
        call check_iarg('which_iter',     self%which_iter)
        call check_iarg('xdim',           self%xdim)
        call check_iarg('xcoord',         self%xcoord)
        call check_iarg('ycoord',         self%ycoord)
        call check_iarg('ydim',           self%ydim)
        ! Float args
        call check_rarg('alpha',          self%alpha)
        call check_rarg('amsklp',         self%amsklp)
        call check_rarg('angerr',         self%angerr)
        call check_rarg('ares',           self%ares)
        call check_rarg('astigerr',       self%astigerr)
        call check_rarg('astigtol',       self%astigtol)
        call check_rarg('athres',         self%athres)
        call check_rarg('batchfrac',      self%batchfrac)
        call check_rarg('bfac',           self%bfac)
        call check_rarg('bfacerr',        self%bfacerr)
        call check_rarg('cenlp',          self%cenlp)
        call check_rarg('cs',             self%cs)
        call check_rarg('ctfreslim',      self%ctfreslim)
        call check_rarg('dcrit_rel',      self%dcrit_rel)
        call check_rarg('deflim',         self%deflim)
        call check_rarg('defocus',        self%defocus)
        call check_rarg('dfclose',        self%dfclose)
        call check_rarg('dffar',          self%dffar)
        call check_rarg('dens',           self%dens)
        call check_rarg('dferr',          self%dferr)
        call check_rarg('dfmax',          self%dfmax)
        call check_rarg('dfmin',          self%dfmin)
        call check_rarg('dfsdev',         self%dfsdev)
        call check_rarg('dfstep',         self%dfstep)
        call check_rarg('dose_rate',      self%dose_rate)
        call check_rarg('e1',             self%e1)
        call check_rarg('e2',             self%e2)
        call check_rarg('e3',             self%e3)
        call check_rarg('eps',            self%eps)
        call check_rarg('exp_time',       self%exp_time)
        call check_rarg('filwidth',       self%filwidth)
        call check_rarg('frac',           self%frac)
        call check_rarg('fraca',          self%fraca)
        call check_rarg('fracdeadhot',    self%fracdeadhot)
        call check_rarg('frac_outliers',  self%frac_outliers)
        call check_rarg('fraczero',       self%fraczero)
        call check_rarg('ftol',           self%ftol)
        call check_rarg('gw',             self%gw)
        call check_rarg('hp',             self%hp)
        call check_rarg('hp_ctffind',     self%hp_ctffind)
        call check_rarg('inner',          self%inner)
        call check_rarg('kv',             self%kv)
        call check_rarg('lam',            self%lam)
        call check_rarg('lp',             self%lp)
        call check_rarg('lp_ctffind',     self%lp_ctffind)
        call check_rarg('lp_grid',        self%lp_grid)
        call check_rarg('lp_pick',        self%lp_pick)
        call check_rarg('lpstart',        self%lpstart)
        call check_rarg('lpstop',         self%lpstop)
        call check_rarg('lpvalid',        self%lpvalid)
        call check_rarg('moldiam',        self%moldiam)
        call check_rarg('msk',            self%msk)
        call check_rarg('mul',            self%mul)
        call check_rarg('mw',             self%mw)
        call check_rarg('ndev',           self%ndev)
        call check_rarg('neigh',          self%neigh)
        call check_rarg('nsig',           self%nsig)
        call check_rarg('outer',          self%outer)
        call check_rarg('phranlp',        self%phranlp)
        call check_rarg('power',          self%power)
        call check_rarg('scale',          self%scale)
        call check_rarg('scale2',         self%scale2)
        call check_rarg('sherr',          self%sherr)
        call check_rarg('smpd',           self%smpd)
        call check_rarg('snr',            self%snr)
        call check_rarg('thres',          self%thres)
        call check_rarg('time_per_image', self%time_per_image)
        call check_rarg('trs',            self%trs)
        call check_rarg('update_frac',    self%update_frac)
        call check_rarg('width',          self%width)
        call check_rarg('winsz',          self%winsz)
        call check_rarg('xsh',            self%xsh)
        call check_rarg('ysh',            self%ysh)
        call check_rarg('zsh',            self%zsh)

!>>> START, SANITY CHECKING AND PARAMETER EXTRACTION FROM ORITAB(S)/VOL(S)/STACK(S)
        ! put unidoc (if defined as oritab)
        if( cline%defined('unidoc') )then
            if( .not. cline%defined('oritab') )then
                call cline%set('oritab', self%unidoc)
                self%oritab = self%unidoc
            else
                write(*,*) 'WARNING! Could not set unidoc to oritab because oritab is defined'
            endif
        endif
        ! put ctffind_doc (if defined) as oritab
        if( cline%defined('ctffind_doc') )then
            if( .not. cline%defined('oritab') )then
                call cline%set('oritab', self%ctffind_doc)
                self%oritab = self%ctffind_doc
            else
                write(*,*) 'WARNING! Could not set ctffind_doc to oritab because oritab is defined'
            endif
        endif
        ! make all programs have the simple_prefix
        if(cline%defined('prg'))then
            if(.not. str_has_substr(self%prg, 'simple_')) self%prg = 'simple_'//trim(self%prg)
        endif
        ! check nr of states
        if( cline%defined('nstates') )then
            self%nstates = nint(cline%get_rarg('nstates'))
            if( self%nstates > MAXS )then
                write(*,'(a)') 'ERROR, the number of states is limited to 20!'
                write(*,'(a)') 'In: constructor, module: simple_params.f90'
                stop
            endif
        endif
        ! determines whether at least one volume is on the cmdline
        do i=1,self%nstates
            if( cline%defined( trim('vol')//int2str(i) ))vol_defined(i) = .true.
        enddo
        ! check inputted vols
        if( cline%defined('vollist') )then
            if( nlines(self%vollist)< MAXS )then
                call read_vols
            endif
        else
            if( cline%defined('vol') )then
                self%vols(1) = self%vol
            endif
            if( cline%defined('vol') .or. any(vol_defined) )then
                do i=1,MAXS
                    call check_vol( i )
                end do
            endif
        endif
        ! no stack given, get ldim from volume if present
        if( self%stk .eq. '' .and. self%vols(1) .ne. '' )then
            call find_ldim_nptcls(self%vols(1), self%ldim, ifoo)
            self%box  = self%ldim(1)
            DebugPrint 'found logical dimension of volume: ', self%ldim
        endif
        ! take care of nptcls
        if( self%stk .ne. '' )then
            if( file_exists(self%stk) )then
                if( .not. cline%defined('nptcls') )then
                    ! get number of particles from stack
                     call find_ldim_nptcls(self%stk, lfoo, self%nptcls)
                     DebugPrint 'found logical dimension of stack: ', lfoo
                     DebugPrint 'found nr of ptcls from stack: ', self%nptcls
                endif
            else
                write(*,'(a,1x,a)') 'Inputted stack file does not exist!', trim(self%stk)
                stop
            endif
        else if( self%oritab .ne. '' )then
            ! get nr of particles from oritab
            if( .not. cline%defined('nptcls') )then
                if( str_has_substr(self%oritab,'.txt') )then
                    self%nptcls = nlines(self%oritab)
                else
                    ! needed because use of binoris_io causes circular dependency
                    ! since params is used by prime3D_srch
                    call bos%open(self%oritab)
                    self%nptcls = bos%get_n_records()
                    call bos%close
                endif            
            endif
        else if( self%refs .ne. '' )then
            if( file_exists(self%refs) )then
                if( .not. cline%defined('box') )then
                    call find_ldim_nptcls(self%refs, self%ldim, ifoo)
                    self%ldim(3) = 1
                    DebugPrint 'found logical dimension of refs: ', self%ldim
                    self%box = self%ldim(1)
                endif
            else
                write(*,'(a,1x,a)') 'Inputted stack file does not exist!', trim(self%refs)
                stop
            endif
        endif
        ! check file formats
        call check_file_formats(aamix)
        call double_check_file_formats
        ! make file names
        call mkfnames
        ! check box
        if( self%box > 0 .and. self%box < 26 ) stop 'box size need to be larger than 26; simple_params'
        ! fractional search and volume update
        if( self%update_frac <= .99)then
            if( self%update_frac < 0.01 )stop 'UPDATE_FRAC is too small 1; simple_params :: constructor'
            if( nint(self%update_frac*real(self%nptcls)) < 1 )&
                &stop 'UPDATE_FRAC is too small 2; simple_params :: constructor'
            self%l_frac_update = .true.
        endif
!<<< END, SANITY CHECKING AND PARAMETER EXTRACTION FROM VOL(S)/STACK(S)

!>>> START, PARALLELISATION-RELATED
        ! set split mode (even)
        self%split_mode = 'even'
        nparts_set      = .false.
        if( cline%defined('nparts') ) nparts_set  = .true.
        ! set execution mode (shmem or distr)
        if( cline%defined('part') .and. cline%defined('nparts') )then
            self%l_distr_exec = .true.
        else
            self%l_distr_exec = .false.
        endif
        l_distr_exec_glob = self%l_distr_exec
        ! set lenght of number string for zero padding
        if( .not. cline%defined('numlen') )then
            if( nparts_set ) self%numlen = len(int2str(self%nparts))
        endif
        self%l_stktab_input = .false.
        if( cline%defined('stktab') )then
            if( cline%defined('stk') )&
                &stop 'stk and stktab cannot simultaneously be defined; params :: new'
            ! prepare stktab handler and set everything accordingly in this instance
            call self%stkhandle%new(trim(self%stktab))
            self%nmics   = self%stkhandle%get_nmics()
            self%nptcls  = self%stkhandle%get_nptcls()
            self%ldim    = self%stkhandle%get_ldim()
            self%ldim(3) = 1
            self%box     = self%ldim(1)
            self%l_stktab_input = .true.
        else
            ! set name of partial files in parallel execution
            stk_part_fname_sc = trim(STKPARTFBODY_SC)//int2str_pad(self%part,self%numlen)//self%ext
            stk_part_fname    = trim(STKPARTFBODY)//int2str_pad(self%part,self%numlen)//self%ext
            self%stk_part     = stk_part_fname
            if( self%autoscale .eq. 'yes' )then
                if( file_exists(stk_part_fname_sc) )then
                    self%stk_part = stk_part_fname_sc
                endif
            endif
            call set_ldim_box_from_stk( self%stk )
            self%box_original = self%box
            if( file_exists(self%stk_part) .and. cline%defined('nparts') )then
                call set_ldim_box_from_stk( self%stk_part )
                if( cline%defined('stk') .and. self%autoscale .eq. 'no' )then
                    if( self%box /= self%box_original )then
                        write(*,*) 'original box:                ', self%box_original
                        write(*,*) 'box read from partial stack: ', self%box
                        stop 'dim mismatch; simple_params :: new'
                    endif
                endif
            endif
            ! Check for the existance of this file if part is defined on the command line
            if( cline%defined('part') )then
                if( ccheckdistr )then
                    if(trim(self%prg).eq.'simple_symsrch')then
                        ! no need for split stack with prg=symsrch
                    else if( .not. file_exists(self%stk_part) )then
                        write(*,*) 'Need partial stacks to be generated for parallel execution'
                        write(*,*) 'Use simple_exec prg=split'
                        stop
                    endif
                endif
            endif
        endif
        ! if we are doing chunk-based parallelisation...
        self%l_chunk_distr = .false.
        if( cline%defined('chunksz') )then
            self%l_chunk_distr = .true.
        endif
        if( .not. cline%defined('ncunits') )then
            ! we assume that the number of computing units is equal to the number of partitions
            self%ncunits = self%nparts
        endif
        ! OpenMP threads
        if( cline%defined('nthr') )then
!$          call omp_set_num_threads(self%nthr)
        else
!$          self%nthr = omp_get_max_threads()
!$          call omp_set_num_threads(self%nthr)
        endif
        nthr_glob = self%nthr

!<<< END, PARALLELISATION-RELATED

!>>> START, IMAGE-PROCESSING-RELATED
        if( .not. cline%defined('xdim') ) self%xdim = self%box/2
        self%xdimpd = round2even(self%alpha*real(self%box/2))
        self%boxpd  = 2*self%xdimpd
        ! set derived Fourier related variables
        self%dstep   = real(self%box-1)*self%smpd                  ! first wavelength of FT
        self%dsteppd = real(self%boxpd-1)*self%smpd                ! first wavelength of padded FT
        if( .not. cline%defined('hp') ) self%hp = 0.7*self%dstep   ! high-pass limit
        self%fny = 2.*self%smpd                                    ! Nyqvist limit
        if( .not. cline%defined('lpstop') )then                    ! default highest resolution lp
            self%lpstop = self%fny                                 ! deafult lpstop
        endif
        if( self%fny > 0. ) self%tofny = nint(self%dstep/self%fny) ! Nyqvist Fourier index
        if( cline%defined('lp') ) self%dynlp = 'no'                ! override dynlp=yes and lpstop
        ! set 2D low-pass limits and smpd_targets 4 scaling
        self%lplims2D(1)       = self%lpstart
        self%lplims2D(2)       = self%lplims2D(1) - (self%lpstart - self%lpstop)/2.
        self%lplims2D(3)       = self%lpstop
        self%smpd_targets2D(1) = self%lplims2D(2)*LP2SMPDFAC
        self%smpd_targets2D(2) = self%lplims2D(3)*LP2SMPDFAC
        ! set default ring2 value
        if( .not. cline%defined('ring2') )then
            if( cline%defined('msk') )then
                self%ring2 = round2even(self%msk)
            else
                self%ring2 = round2even(real(self%box/2-2))
            endif
        else
            self%ring2 = round2even(real(min(self%box/2-2,self%ring2)))
        endif
        ! set default msk value
        if( .not. cline%defined('msk') )then
            if( cline%defined('ring2') )then
                self%msk = self%ring2
            else
                self%msk = self%box/2
            endif
        endif
        ! set mode of masking
        if( cline%defined('inner') )then
            self%l_innermsk = .true.
        else
            self%l_innermsk = .false.
        endif
        ! set nr of rotations
        self%nrots = round2even(twopi*real(self%ring2))
        if( self%box > 2*(nint(self%msk)+10) )then
            self%boxmatch = find_magic_box(2*(nint(self%msk)+10))
            if( self%boxmatch > self%box ) self%boxmatch = self%box
        else
            self%boxmatch = self%box
        endif
        ! set default outer mask value
        if( .not. cline%defined('outer') ) self%outer = self%msk
        ! matched filter flag
        select case(self%match_filt)
            case('no')
                self%l_match_filt = .false.
            case DEFAULT
                self%l_match_filt = .true.
        end select
        ! checks automask related values
        self%l_envmsk = .false.
        if( self%automsk .ne. 'no' ) self%l_envmsk = .true.
        if( cline%defined('mskfile') )then
            if( .not. file_exists(self%mskfile) )then
                write(*,*) 'file: ', trim(self%mskfile)
                stop 'input mask file not in cwd'
            endif
            self%l_envmsk = .true.
        endif
        ! scaling stuff
        self%l_autoscale = .false.
        if( self%autoscale .eq. 'yes' ) self%l_autoscale = .true.
        if( .not. cline%defined('newbox') )then
            ! set newbox if scale is defined
            if( cline%defined('scale') )then
                self%newbox = find_magic_box(nint(self%scale*real(self%box)))
            endif
        endif
        ! set newbox if scale is defined
        if( cline%defined('scale2') )then
            self%newbox2 = find_magic_box(nint(self%scale2*real(self%box)))
        endif
        self%kfromto(1) = max(2,int(self%dstep/self%hp)) ! high-pass Fourier index set according to hp
        self%kfromto(2) = int(self%dstep/self%lp)        ! low-pass Fourier index set according to lp
        self%lp         = max(self%fny,self%lp)          ! lowpass limit
        self%lp_dyn     = self%lp                        ! dynamic lowpass limit
        self%lpmed      = self%lp                        ! median lp
        if( .not. cline%defined('ydim') ) self%ydim = self%xdim
        ! set ldim
        if( cline%defined('xdim') ) self%ldim = [self%xdim,self%ydim,1]
        ! if CTF refinement
        if( self%ctf .eq. 'refine' ) stop 'CTF refinement not yet implemented'
        ! set eullims
        self%eullims(:,1) = 0.
        self%eullims(:,2) = 359.99
        self%eullims(2,2) = 180.
        ! set default size of random sample
        if( .not. cline%defined('nran') )then
            self%nran = self%nptcls
        endif
        ! define clipped box if not given
        if( .not. cline%defined('clip') )then
            self%clip = self%box
        endif
        ! error check ncls
        if( file_exists(self%refs) )then
            ! get number of particles from stack
            call find_ldim_nptcls(self%refs, lfoo, ncls)
            DebugPrint 'found ncls from refs: ', ncls
            if( cline%defined('ncls') )then
                if( ncls /= self%ncls ) stop 'inputtend number of clusters (ncls) not&
                &consistent with the number of references in stack (p%refs)'
            else
                self%ncls = ncls
            endif
        endif
        ! set remap_classes flag
        self%l_remap_classes = .false.
        if( self%remap_classes .eq. 'yes' ) self%l_remap_classes = .true.
        ! set nbest to 20% of ncls (if present) or 20% of NSPACE_BALANCE
        if( cline%defined('ncls') )then
            self%nbest = max(1,nint(real(self%ncls)*0.2))
        else
            self%nbest = max(1,nint(real(NSPACE_BALANCE)*0.2))
        endif
        ! set to particle index if not defined in cmdlin
        if( .not. cline%defined('top') ) self%top = self%nptcls
        ! set the number of input orientations
        if( .not. cline%defined('noris') ) self%noris = self%nptcls
        ! fix translation param
        self%trs = abs(self%trs)
        if( .not. cline%defined('trs') )then
            if( self%refine.eq.'no' .or. self%refine.eq.'snhc' )then
                self%trs = 0.
            else
                self%trs = 1.
            endif
        endif
        self%doshift = .true.
        if( self%trs < 0.5 )then
            self%trs = 0.00001
            self%doshift = .false.
        endif
        ! set optlims
        self%optlims(:3,:)  = self%eullims
        self%optlims(4:5,1) = -self%trs
        self%optlims(4:5,2) = self%trs
        self%optlims(6:7,1) = -self%dfsdev
        self%optlims(6:7,2) = self%dfsdev
        ! Set first three in cyclic true
        self%cyclic(1) = .true.
        self%cyclic(2) = .true.
        self%cyclic(3) = .true.
        ! Set molecular diameter
        if( .not. cline%defined('moldiam') )then
            self%moldiam = 2. * self%msk * self%smpd
        endif
        ! check if we are dose-weighting or not
        self%l_dose_weight = .false.
        if( cline%defined('exp_time') .or. cline%defined('dose_rate') )then
            if( cline%defined('exp_time') .and. .not. cline%defined('dose_rate') )then
                stop 'need dose_rate to be part of the command line! simple_params :: new'
            endif
            if( .not. cline%defined('exp_time') .and. cline%defined('dose_rate') )then
                stop 'need exp_time to be part of the command line! simple_params :: new'
            endif
            self%l_dose_weight = .true.
        endif
        ! prepare CTF plan
        self%tfplan%flag = self%ctf
        ! set logical pick flag
        self%l_pick = .false.
        if( self%dopick .eq. 'yes' ) self%l_pick = .true.
!>>> END, IMAGE-PROCESSING-RELATED

        write(*,'(A)') '>>> DONE PROCESSING PARAMETERS'

        contains

            subroutine check_vol( i )
                integer, intent(in) :: i
                character(len=STDLEN) :: nam
                nam = 'vol'//int2str(i)
                if( cline%defined(nam) )then
                    call check_file(nam, self%vols(i), notAllowed='T')
                    if( .not. file_exists(self%vols(i)) )then
                        write(*,*) 'Input volume:', self%vols(i), 'does not exist!'
                        stop
                    endif
                    DebugPrint nam, '=', self%vols(i)
                endif
            end subroutine check_vol

            subroutine read_vols
                character(len=STDLEN) :: filenam, nam
                integer :: nl, fnr, i, io_stat
                filenam = cline%get_carg('vollist')
                nl      = nlines(filenam)
                call fopen(fnr, file=filenam, iostat=io_stat)
                if(io_stat /= 0) call fileio_errmsg("params ; read_vols error opening "//trim(filenam), io_stat)
                do i=1,nl
                    read(fnr,*, iostat=io_stat) nam
                    if(io_stat /= 0) call fileio_errmsg("params ; read_vols error reading "//trim(filenam), io_stat)
                    if( nam .ne. '' )then
                        self%vols(i) = nam
                    endif
                end do
                call fclose(fnr,errmsg="params ; read_vols error closing "//trim(filenam))
            end subroutine read_vols

            subroutine check_file( file, var, allowed1, allowed2, notAllowed )
                character(len=*),           intent(in)  :: file
                character(len=*),           intent(out) :: var
                character(len=1), optional, intent(in)  :: allowed1, allowed2, notAllowed
                character(len=1) :: file_descr
                logical          :: raise_exception
                if( cline%defined(file) )then
                    var = cline%get_carg(file)
                    DebugPrint 'var = ', var
                    file_descr = fname2format(var)
                    DebugPrint 'file_descr = ', file_descr
                    raise_exception = .false.
                    if( present(allowed1) )then
                        if( present(allowed2) )then
                            if( allowed1 == file_descr .or. allowed2 == file_descr )then
                                ! all good
                            else
                                raise_exception = .true.
                            endif
                        else
                            if( allowed1 /= file_descr ) raise_exception = .true.
                        endif
                    endif
                    if( present(notAllowed) )then
                        if( notAllowed == file_descr ) raise_exception = .true.
                    endif
                    if( raise_exception )then
                        write(*,*) 'This format: ', file_descr, ' is not allowed for this file: ', var
                        stop
                    endif
                    select case(file_descr)
                        case ('I')
                            stop 'Support for IMAGIC files is not yet implemented!'
                        case ('M')
                            ! MRC files are supported
                            cntfile = cntfile+1
                            checkupfile(cntfile) = 'M'
                        case ('S')
                            ! SPIDER files are supported
                            cntfile = cntfile+1
                            checkupfile(cntfile) = 'S'
                        case ('N')
                            write(*,*) 'file: ', trim(file)
                            stop 'This file format is not supported by SIMPLE; simple_params::check_file'
                        case ('T')
                            ! text files are supported
                        case ('B')
                            ! binary files are supported
                        case ('P')
                            ! PDB files
                        case DEFAULT
                            write(*,*) 'file: ', trim(file)
                            stop 'This file format is not supported by SIMPLE; simple_params::check_file'
                    end select
                    DebugPrint file, '=', var
                endif
            end subroutine check_file


            subroutine check_file_formats( allow_mix )
                logical, intent(in) :: allow_mix
                integer :: i, j
                if( cntfile > 0 )then
                    do i=1,cntfile
                        do j=1,cntfile
                            if( i == j ) cycle
                            if( checkupfile(i) == checkupfile(j) )then
                                ! all ok
                            else
                                if( .not. allow_mix )then
                                    stop 'The inputted file names have nonconforming format (mixed formats not yet allowed)!'
                                endif
                            endif
                        end do
                    end do
                    select case(checkupfile(1))
                        case('M')
                            self%ext = '.mrc'
                        case('S')
                            self%ext = '.spi'
                        case('D')
                            self%ext = '.mrc'
                        case('B')
                            self%ext = '.mrc'
                        case DEFAULT
                            stop 'This file format is not supported by SIMPLE; check_file_formats; simple_params'
                    end select
                endif
            end subroutine check_file_formats

            subroutine double_check_file_formats
                character(len=STDLEN) :: fname
                character(len=1)      :: form
                integer :: funit, io_stat
                if( cntfile == 0 )then
                    if( cline%defined('filetab') )then
                        call fopen(funit, status='old', file=self%filetab, iostat=io_stat)
                        call fileio_errmsg("In params:: double_check_file_formats fopen failed "//trim(self%filetab) , io_stat)
                        read(funit,'(a256)') fname
                        form = fname2format(fname)
                        call fclose(funit, &
                            errmsg="In params:: double_check_file_formats fclose failed "//trim(self%filetab) )
                        select case(form)
                            case('M')
                                self%ext = '.mrc'
                            case('S')
                                self%ext = '.spi'
                            case('D')
                                self%ext = '.mrc'
                            case('B')
                                self%ext = '.mrc'
                            case DEFAULT
                                write(*,*) 'format string is ', form
                                stop 'This file format is not supported by SIMPLE; double_check_file_formats; simple_params'
                        end select
                    endif
                endif
            end subroutine double_check_file_formats

            subroutine mkfnames
                if( .not. cline%defined('outstk')  ) self%outstk  = 'outstk'//self%ext
                if( .not. cline%defined('outstk2') ) self%outstk2 = 'outstk2'//self%ext
                if( .not. cline%defined('outvol')  ) self%outvol  = 'outvol'//self%ext
            end subroutine mkfnames

            subroutine check_carg( carg, var )
                character(len=*), intent(in)  :: carg
                character(len=*), intent(out) :: var
                if( cline%defined(carg) )then
                    var = cline%get_carg(carg)
                    DebugPrint carg, '=', var
                endif
            end subroutine check_carg

            subroutine check_iarg( iarg, var )
                character(len=*), intent(in)  :: iarg
                integer, intent(out) :: var
                if( cline%defined(iarg) )then
                    var = nint(cline%get_rarg(iarg))
                    DebugPrint iarg, '=', var
                endif
            end subroutine check_iarg

            subroutine check_larg( larg, var )
                character(len=*), intent(in)  :: larg
                logical, intent(out) :: var
                integer :: tmp
                if( cline%defined(larg) )then
                    tmp =  NINT( cline%get_rarg(larg) )
                    var = tmp /= 0
                    DebugPrint larg, '=', var
                endif
            end subroutine check_larg

            subroutine check_rarg( rarg, var )
                character(len=*), intent(in)  :: rarg
                real, intent(out) :: var
                if( cline%defined(rarg) )then
                    var = cline%get_rarg(rarg)
                    DebugPrint rarg, '=', var
                endif
            end subroutine check_rarg

            subroutine set_ldim_box_from_stk( stkfname )
                character(len=*), intent(in) :: stkfname
                if( stkfname .ne. '' )then
                    if( file_exists(stkfname) )then
                        if( cline%defined('box') )then
                        else
                            call find_ldim_nptcls(stkfname, self%ldim, ifoo)
                            self%ldim(3) = 1
                            DebugPrint 'found logical dimension of stack: ', self%ldim
                            self%box     = self%ldim(1)
                        endif
                    else
                        write(*,'(a)')      'simple_params :: set_ldim_box_from_stk'
                        write(*,'(a,1x,a)') 'Stack file does not exist!', trim(stkfname)
                        stop
                    endif
                endif
            end subroutine set_ldim_box_from_stk

    end subroutine new

end module simple_params
