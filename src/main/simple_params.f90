! provides global distribution of constants and derived constants
module simple_params
include 'simple_lib.f08'
use simple_ori,     only: ori
use simple_cmdline, only: cmdline
use simple_binoris, only: binoris
use simple_user_interface
!$ use omp_lib
!$ use omp_lib_kinds
implicit none

public :: params
private
#include "simple_local_flags.inc"

!> global parameters
type :: params
    ! global objects
    type(ori)                      :: ori_glob
    type(ctfplan)                  :: tfplan
    class(simple_program), pointer :: ptr2prg
    ! yes/no decision variables in ascending alphabetical order
    character(len=:), allocatable :: acf             !< calculate autocorrelation function(yes|no){no}
    character(len=:), allocatable :: append          !< append in context of files(yes|no){no}
    character(len=:), allocatable :: async           !< asynchronous (yes|no){no}
    character(len=:), allocatable :: autoscale       !< automatic down-scaling(yes|no){yes}
    character(len=:), allocatable :: avg             !< calculate average (yes|no){no}
    character(len=:), allocatable :: bin             !< binarise image(yes|no){no}
    character(len=:), allocatable :: center          !< center image(s)/class average(s)/volume(s)(yes|no){no}
    character(len=:), allocatable :: classtats       !< calculate class population statistics(yes|no){no}
    character(len=:), allocatable :: clustvalid      !< validate clustering(yes|homo|no){no}
    character(len=:), allocatable :: compare         !< do comparison(yes|no){no}
    character(len=:), allocatable :: countvox        !< count # voxels(yes|no){no}
    character(len=:), allocatable :: ctfstats        !< calculate ctf statistics(yes|no){no}
    character(len=:), allocatable :: cure
    character(len=:), allocatable :: dev             !< development flag for experimental code(yes|no){no}
    character(len=:), allocatable :: discrete        !< be discrete(yes|no){no}
    character(len=:), allocatable :: diverse         !< diverse or not flag (yes|no){no}
    character(len=:), allocatable :: doalign
    character(len=:), allocatable :: dopca
    character(len=:), allocatable :: doprint
    character(len=:), allocatable :: dryrun          !< no 3D search and reconstruction, for profiling only
    character(len=:), allocatable :: dynlp           !< automatic resolution limit update(yes|no){yes}
    character(len=:), allocatable :: dyncls          !< dynamic class update(yes|no){yes}
    character(len=:), allocatable :: errify          !< introduce error(yes|no){no}
    character(len=:), allocatable :: even            !< even orientation distribution(yes|no){no}
    character(len=:), allocatable :: ft2img          !< convert Fourier transform to real image of power(yes|no){no}
    character(len=:), allocatable :: for3D           !< for 3D analysis(yes|no){yes}
    character(len=:), allocatable :: guinier         !< calculate Guinier plot(yes|no){no}
    character(len=:), allocatable :: kmeans
    character(len=:), allocatable :: local
    character(len=:), allocatable :: masscen         !< center using binarisation and mass centering(yes|no){no}
    character(len=:), allocatable :: match_filt      !< matched filter on (yes|no){yes}
    character(len=:), allocatable :: merge
    character(len=:), allocatable :: mirr            !< mirror(no|x|y){no}
    character(len=:), allocatable :: mkdir           !< make auto-named directory for execution and gathering output(yes|no){no}
    character(len=:), allocatable :: neg             !< invert contrast of images(yes|no)
    character(len=:), allocatable :: neigh           !< neighbourhood refinement(yes|no){no}
    character(len=:), allocatable :: noise_norm
    character(len=:), allocatable :: norec           !< do not reconstruct volume(s)(yes|no){no}
    character(len=:), allocatable :: norm            !< do statistical normalisation avg
    character(len=:), allocatable :: order           !< order ptcls according to correlation(yes|no){no}
    character(len=:), allocatable :: outside         !< extract boxes outside the micrograph boundaries(yes|no){no}
    character(len=:), allocatable :: pad
    character(len=:), allocatable :: pgrp_known      !< point-group known a priori(yes|no){no}
    character(len=:), allocatable :: phaseplate      !< images obtained with Volta phaseplate(yes|no){no}
    character(len=:), allocatable :: phrand          !< phase randomize(yes|no){no}
    character(len=:), allocatable :: plot            !< make plot(yes|no){no}
    character(len=:), allocatable :: projstats
    character(len=:), allocatable :: readwrite
    character(len=:), allocatable :: remap_cls
    character(len=:), allocatable :: restart
    character(len=:), allocatable :: rnd             !< random(yes|no){no}
    character(len=:), allocatable :: rm_outliers     !< remove outliers{yes}
    character(len=:), allocatable :: roalgn
    character(len=:), allocatable :: round
    character(len=:), allocatable :: shalgn          !< do 2D shift alignment(yes|no){no}
    character(len=:), allocatable :: shellnorm
    character(len=:), allocatable :: shbarrier       !< use shift search barrier constraint(yes|no){yes}
    character(len=:), allocatable :: soften          !< soften envelope with cosine edge(yes|no){no}
    character(len=:), allocatable :: stats           !< provide statistics(yes|no){yes}
    character(len=:), allocatable :: stream          !< sream (real time) execution mode(yes|no){no}
    character(len=:), allocatable :: symrnd          !< randomize over symmetry operations(yes|no){no}
    character(len=:), allocatable :: swap
    character(len=:), allocatable :: taper_edges     !< self-explanatory
    character(len=:), allocatable :: test
    character(len=:), allocatable :: tomo            !< tomography mode(yes|no){no}
    character(len=:), allocatable :: time
    character(len=:), allocatable :: trsstats        !< provide origin shift statistics(yes|no){no}
    character(len=:), allocatable :: tseries         !< images represent a time-series(yes|no){no}
    character(len=:), allocatable :: vis             !< visualise(yes|no)
    character(len=:), allocatable :: weights2D       !< to use 2d spectral weights(yes|no){no}
    character(len=:), allocatable :: weights3D       !< to use 3d spectral weights(yes|no){yes}
    character(len=:), allocatable :: zero            !< zeroing(yes|no){no}
    ! other character variables in ascending alphabetical order
    character(len=:), allocatable :: angastunit      !< angle of astigmatism unit (radians|degrees){degrees}
    character(len=:), allocatable :: automsk
    character(len=:), allocatable :: boxfile         !< file with EMAN particle coordinates(.txt)
    character(len=:), allocatable :: boxtab          !< table (text file) of files with EMAN particle coordinates(.txt)
    character(len=:), allocatable :: boxtype
    character(len=:), allocatable :: classdoc        !< doc with per-class stats(.txt)
    character(len=:), allocatable :: comlindoc       !< clustering_shc_nclsX.txt
    character(len=:), allocatable :: ctf             !< ctf flag(yes|no|flip)
    character(len=:), allocatable :: cwd
    character(len=:), allocatable :: deftab          !< file with CTF info(.txt|.simple)
    character(len=:), allocatable :: dfunit          !< defocus unit (A|microns){microns}
    character(len=:), allocatable :: dir             !< directory
    character(len=:), allocatable :: dir_movies      !< grab mrc mrcs files from here
    character(len=:), allocatable :: dir_reject      !< move rejected files to here{rejected}
    character(len=:), allocatable :: dir_select      !< move selected files to here{selected}
    character(len=:), allocatable :: dir_target      !< put output here
    character(len=:), allocatable :: dir_ptcls
    character(len=:), allocatable :: dir_mics        !< grab micrographs from here
    character(len=:), allocatable :: dockmode        !< volume docking mode(eul|shift|eulshift|all){eul}
    character(len=:), allocatable :: doclist         !< list of oritabs for different states
    character(len=:), allocatable :: eo              !< use FSC for filtering and low-pass limit update(yes|aniso|no){no}
    character(len=:), allocatable :: exec_dir
    character(len=:), allocatable :: executable
    character(len=:), allocatable :: exp_doc         !< specifying exp_time and dose_rate per tomogram
    character(len=:), allocatable :: ext             !< file extension{.mrc}
    character(len=:), allocatable :: extrmode
    character(len=:), allocatable :: fbody           !< file body
    character(len=:), allocatable :: featstk
    character(len=:), allocatable :: filetab         !< list of files(.txt)
    character(len=:), allocatable :: fname           !< file name
    character(len=:), allocatable :: frcs            !< binary file with per-class/proj Fourier Ring Correlations(.bin)
    character(len=:), allocatable :: fsc             !< binary file with FSC info{fsc_state01.bin}
    character(len=:), allocatable :: hfun            !< function used for normalization(sigm|tanh|lin){sigm}
    character(len=:), allocatable :: hist            !< give variable for histogram plot
    character(len=:), allocatable :: imgkind         !< type of image(ptcl|cavg|mic|movie){ptcl}
    character(len=:), allocatable :: infile          !< file with inputs(.txt|.simple)
    character(len=:), allocatable :: label           !< discrete label(class|state){class}
    character(len=:), allocatable :: mskfile         !< maskfile.ext
    character(len=:), allocatable :: msklist         !< table (text file) of mask volume files(.txt)
    character(len=:), allocatable :: msktype         !< type of mask(hard|soft){soft}
    character(len=:), allocatable :: objfun          !< objective function(cc|ccres){cc}
    character(len=:), allocatable :: opt             !< optimiser (bfgs|simplex){bfgs}
    character(len=:), allocatable :: oritab          !< table  of orientations(.txt|.simple)
    character(len=:), allocatable :: oritab2         !< 2nd table of orientations(.txt|.simple)
    character(len=:), allocatable :: oritab3D        !< table of 3D orientations(.txt|.simple)
    character(len=:), allocatable :: oritype         !< SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D)
    character(len=:), allocatable :: outfile         !< output document
    character(len=:), allocatable :: outstk          !< output image stack
    character(len=:), allocatable :: outstk2         !< output image stack 2nd
    character(len=:), allocatable :: outvol          !< output volume{outvol.ext}
    character(len=:), allocatable :: ctf_estimate_doc!< per-micrograph CTF parameters to transfer
    character(len=:), allocatable :: pcastk
    character(len=:), allocatable :: pcontrast       !< particle contrast(black|white){black}
    character(len=:), allocatable :: pdbfile         !< PDB file
    character(len=:), allocatable :: pdfile
    character(len=:), allocatable :: pgrp            !< point-group symmetry(cn|dn|t|o|i)
    character(len=:), allocatable :: phshiftunit     !< additional phase-shift unit (radians|degrees){radians}
    character(len=:), allocatable :: plaintexttab    !< plain text file of input parameters
    character(len=:), allocatable :: prg             !< SIMPLE program being executed
    character(len=:), allocatable :: projfile        !< SIMPLE *.simple project file
    character(len=:), allocatable :: projfile_target !< another SIMPLE *.simple project file
    character(len=:), allocatable :: projname        !< SIMPLE  project name
    character(len=:), allocatable :: real_filter
    character(len=:), allocatable :: refine          !< refinement mode(snhc|single|multi|greedy_single|greedy_multi|cluster|clustersym){no}
    character(len=:), allocatable :: refs            !< initial2Dreferences.ext
    character(len=:), allocatable :: refs_even
    character(len=:), allocatable :: refs_odd
    character(len=:), allocatable :: speckind        !< power spectrum kind(real|power|sqrt|log|phase){sqrt}
    character(len=:), allocatable :: split_mode
    character(len=:), allocatable :: stk_part
    character(len=:), allocatable :: stk             !< particle stack with all images(ptcls.ext)
    character(len=:), allocatable :: stktab          !< list of per-micrograph stacks
    character(len=:), allocatable :: stk2            !< 2nd stack(in map2ptcls/select: selected(cavgs).ext)
    character(len=:), allocatable :: stk3            !< 3d stack (in map2ptcls/select: (cavgs)2selectfrom.ext)
    character(len=:), allocatable :: stk_backgr      !< stack with image for background subtraction
    character(len=:), allocatable :: tomoseries      !< filetable of filetables of tomograms
    character(len=:), allocatable :: unidoc          !< unified resources and orientations doc
    character(len=:), allocatable :: vol_filt        !< input filter volume(vol_filt.ext)
    character(len=:), allocatable :: vollist         !< table (text file) of volume files(.txt)
    character(len=:), allocatable :: voltab          !< table (text file) of volume files(.txt)
    character(len=:), allocatable :: voltab2         !< 2nd table (text file) of volume files(.txt)
    character(len=:), allocatable :: wfun
    type(str4arr),    allocatable :: vols(:)
    type(str4arr),    allocatable :: vols_even(:)
    type(str4arr),    allocatable :: vols_odd(:)
    type(str4arr),    allocatable :: mskvols(:)
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
    integer :: box_extract
    integer :: boxpd=0
    integer :: chunksz=0           !< # images/orientations in chunk
    integer :: class=1             !< cluster identity
    integer :: clip=0              !< clipped image box size(in pixels)
    integer :: corner=0            !< corner size(in pixels){0}
    integer :: cube=0              !< side size(in pixels){0}
    integer :: edge=6              !< edge size for softening molecular envelope(in pixels)
    integer :: extr_iter=1
    integer :: find=1              !< Fourier index
    integer :: nframesgrp=0        !< # frames to group before motion_correct(Falcon 3){0}
    integer :: fromf=1             !< start frame index
    integer :: fromp=1             !< start ptcl index
    integer :: fstep=1
    integer :: grow=0              !< # binary layers to grow(in pixels)
    integer :: hpind_fsc           !< high-pass Fourier index for FSC
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
    integer :: newbox=0            !< new box for scaling (by Fourier padding/clipping)
    integer :: newbox2=0
    integer :: nframes=0           !< # frames{30}
    integer :: nmembers=0
    integer :: nnn=200             !< # nearest neighbors{200}
    integer :: nmics=0             !< # micographs
    integer :: noris=0
    integer :: nparts=1            !< # partitions in distributed exection
    integer :: npix=0              !< # pixles/voxels in binary representation
    integer :: nptcls=1            !< # images in stk/# orientations in oritab
    integer :: nptcls_per_cls=400  !< # images in stk/# orientations in oritab
    integer :: nran=0              !< # random images to select
    integer :: nrefs=100           !< # references used for picking{100}
    integer :: nrestarts=1
    integer :: nrepeats=HETNREPEATS !< # repeats in het_ensemble worklow
    integer :: nrots=0
    integer :: nspace=2500         !< # projection directions
    integer :: nstates=1           !< # states to reconstruct
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
    integer :: pid=0               !< process ID
    integer :: ppca=0
    integer :: pspecsz=512         !< size of power spectrum(in pixels)
    integer :: ptcl=1
    integer :: recl_cgrid=-1
    integer :: ring1=2
    integer :: ring2=0
    integer :: spec=0
    integer :: spproj_a_seg=0      !< sp-project segments that b%a points to
    integer :: startit=1           !< start iterating from here
    integer :: state=1             !< state to extract
    integer :: state2split=0       !< state group to split
    integer :: stepsz=1            !< size of step{0}
    integer :: szsn=SZSN_INIT      !< size of stochastic neighborhood{5}
    integer :: time_inactive=7200  !< end time limit(secs)
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
    real    :: alpha=KBALPHA
    real    :: amsklp=15.          !< low-pass limit for envelope mask generation(in A)
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
    real    :: defocus=2.          !< defocus(in microns){2.}
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
    real    :: fny=0.
    real    :: focusmsk=0.         !< spherical msk for use with focused refinement
    real    :: frac=1.             !< fraction of ptcls(0-1){1}
    real    :: fraca=0.1           !< fraction of amplitude contrast used for fitting CTF{0.1}
    real    :: fracdeadhot=0.05    !< fraction of dead or hot pixels{0.01}
    real    :: frac_outliers=0.
    real    :: fraczero=0.
    real    :: ftol=1e-6
    real    :: motion_correctftol = 1e-2   !< tolerance (gradient) for motion_correctrer
    real    :: motion_correctgtol = 1e-2   !< tolerance (function value) for motion_correctrer
    real    :: gw=0.5
    real    :: hp=100.             !< high-pass limit(in A)
    real    :: hp_fsc=0.           !< FSC high-pass limit(in A)
    real    :: hp_ctf_estimate=30. !< high-pass limit 4 ctf_estimate(in A)
    real    :: inner=0.            !< inner mask radius(in pixels)
    real    :: kv=300.             !< acceleration voltage(in kV){300.}
    real    :: lam=0.5
    real    :: lp_dyn=20.
    real    :: lp=20.              !< low-pass limit(in A)
    real    :: lp_ctf_estimate=5.0 !< low-pass limit 4 ctf_estimate(in A)
    real    :: lp_pick=20.         !< low-pass limit 4 picker(in A)
    real    :: lplim_crit=0.3      !< corr criterion low-pass limit assignment(0.143-0.5){0.3}
    real    :: lplims2D(3)
    real    :: lpmed=20.
    real    :: lpstart=0.          !< start low-pass limit(in A){15}
    real    :: lpstop=8.0          !< stop low-pass limit(in A){8}
    real    :: lpvalid=20.
    real    :: moldiam=140.        !< molecular diameter(in A)
    real    :: moment=0.
    real    :: msk=0.              !< mask radius(in pixels)
    real    :: mul=1.              !< origin shift multiplication factor{1}
    real    :: mw=0.               !< molecular weight(in kD)
    real    :: ndev=2.0            !< # deviations in one-cluster clustering
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
    real    :: tau=TAU_DEFAULT     !< tau fudge factor, controls the sharpness of the orientation weight distribution,
                                   !! smaller number means sharper distribution
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
    logical :: cyclic(7)      = .false.
    logical :: l_distr_exec   = .false.
    logical :: l_match_filt   = .true.
    logical :: l_doshift      = .false.
    logical :: l_focusmsk     = .false.
    logical :: l_autoscale    = .false.
    logical :: l_dose_weight  = .false.
    logical :: l_frac_update  = .false.
    logical :: l_innermsk     = .false.
    logical :: l_remap_cls    = .false.
    logical :: l_cc_objfun    = .true.
    logical :: l_cc_bfac      = .true.
    logical :: l_dev          = .false.
  contains
    procedure, private :: set_char_defaults
    procedure          :: new
end type params

interface params
    module procedure constructor
end interface

contains

    subroutine set_char_defaults( self )
        class(params), intent(inout) :: self
        allocate(self%acf,           source='no')
        allocate(self%append,        source='no')
        allocate(self%async,         source='no')
        allocate(self%autoscale,     source='no')
        allocate(self%avg,           source='no')
        allocate(self%bin,           source='no')
        allocate(self%center,        source='yes')
        allocate(self%classtats,     source='no')
        allocate(self%clustvalid,    source='no')
        allocate(self%compare,       source='no')
        allocate(self%countvox,      source='no')
        allocate(self%ctfstats,      source='no')
        allocate(self%cure,          source='no')
        allocate(self%dev,           source='no')
        allocate(self%discrete,      source='no')
        allocate(self%diverse,       source='no')
        allocate(self%doalign,       source='yes')
        allocate(self%dopca,         source='yes')
        allocate(self%doprint,       source='no')
        allocate(self%dryrun,        source='no')
        allocate(self%dynlp,         source='yes')
        allocate(self%dyncls,        source='yes')
        allocate(self%errify,        source='no')
        allocate(self%even,          source='no')
        allocate(self%ft2img,        source='no')
        allocate(self%for3D,         source='yes')
        allocate(self%guinier,       source='no')
        allocate(self%kmeans,        source='yes')
        allocate(self%local,         source='no')
        allocate(self%masscen,       source='no')
        allocate(self%match_filt,    source='yes')
        allocate(self%merge,         source='no')
        allocate(self%mirr,          source='no')
        allocate(self%mkdir,         source='no')
        allocate(self%neg,           source='no')
        allocate(self%neigh,         source='no')
        allocate(self%noise_norm,    source='no')
        allocate(self%norec,         source='no')
        allocate(self%norm,          source='no')
        allocate(self%order,         source='no')
        allocate(self%outside,       source='no')
        allocate(self%pad,           source='no')
        allocate(self%pgrp_known,    source='no')
        allocate(self%phaseplate,    source='no')
        allocate(self%phrand,        source='no')
        allocate(self%plot,          source='no')
        allocate(self%projstats,     source='no')
        allocate(self%readwrite,     source='no')
        allocate(self%remap_cls,     source='no')
        allocate(self%restart,       source='no')
        allocate(self%rnd,           source='no')
        allocate(self%rm_outliers,   source='yes')
        allocate(self%roalgn,        source='no')
        allocate(self%round,         source='no')
        allocate(self%shalgn,        source='no')
        allocate(self%shellnorm,     source='no')
        allocate(self%shbarrier,     source='yes')
        allocate(self%soften,        source='no')
        allocate(self%stats,         source='no')
        allocate(self%stream,        source='no')
        allocate(self%symrnd,        source='no')
        allocate(self%swap,          source='no')
        allocate(self%taper_edges,   source='no')
        allocate(self%test,          source='no')
        allocate(self%tomo,          source='no')
        allocate(self%time,          source='no')
        allocate(self%trsstats,      source='no')
        allocate(self%tseries,       source='no')
        allocate(self%vis,           source='no')
        allocate(self%weights2D,     source='no')
        allocate(self%weights3D,     source='no')
        allocate(self%zero,          source='no')
        allocate(self%angastunit,    source='degrees')
        allocate(self%automsk,       source='no')
        allocate(self%boxtype,       source='eman')
        allocate(self%ctf,           source='no')
        allocate(self%dfunit,        source='microns')
        allocate(self%dir_reject,    source='rejected')
        allocate(self%dir_select,    source='selected')
        allocate(self%dockmode,      source='eul')
        allocate(self%eo,            source='yes')
        allocate(self%extrmode,      source='all')
        allocate(self%hfun,          source='sigm')
        allocate(self%hist,          source='corr')
        allocate(self%imgkind,       source='ptcl')
        allocate(self%label,         source='class')
        allocate(self%msktype,       source='soft')
        allocate(self%objfun,        source='cc')
        allocate(self%opt,           source='bfgs')
        allocate(self%pcontrast,     source='black')
        allocate(self%pgrp,          source='c1')
        allocate(self%phshiftunit,   source='radians')
        allocate(self%refine,        source='single')
        allocate(self%speckind,      source='sqrt')
        allocate(self%split_mode,    source='even')
        allocate(self%wfun,          source='kb')
    end subroutine set_char_defaults

    !> \brief  is a constructor
    function constructor( cline, allow_mix, del_scaled, spproj_a_seg ) result( self )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: allow_mix, del_scaled
        integer, optional, intent(in)    :: spproj_a_seg
        type(params) :: self
        call self%new( cline, allow_mix, del_scaled, spproj_a_seg)
    end function constructor

    !> \brief  is a constructor
    subroutine new( self, cline, allow_mix, del_scaled, spproj_a_seg )
        class(params),     intent(inout)   :: self
        class(cmdline),    intent(inout)   :: cline
        logical, optional, intent(in)      :: allow_mix, del_scaled
        integer, optional, intent(in)      :: spproj_a_seg
        character(len=STDLEN), allocatable :: sp_files(:)
        character(len=:),      allocatable :: stk_part_fname_sc, phaseplate, ctfflag
        logical,               allocatable :: vol_defined(:)
        character(len=:), allocatable :: debug_local, verbose_local
        character(len=STDLEN) :: stk_part_fname, cwd_tmp
        character(len=1)      :: checkupfile(50)
        type(binoris)         :: bos
        type(ori)             :: o
        integer               :: i, ncls, ifoo, lfoo(3), cntfile, istate
        integer               :: spproj_a_seg_inputted, idir, nsp_files
        logical               :: nparts_set, aamix, ddel_scaled
        ! set character defaults
        call self%set_char_defaults
        nparts_set        = .false.
        debug_local       = 'no'
        verbose_local     = 'no'
        aamix = .false.
        if( present(allow_mix) ) aamix = allow_mix
        ddel_scaled = .false.
        if( present(del_scaled) ) ddel_scaled = del_scaled
        ! sp_project field corresponding to b%a
        self%spproj_a_seg = GENERIC_SEG
        if( present(spproj_a_seg) ) self%spproj_a_seg = spproj_a_seg
        ! file counter
        cntfile = 0
        ! make global ori
        call self%ori_glob%new
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
        self%outfile = 'outfile'//trim(METADATA_EXT)
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
        call check_carg('doprint',        self%doprint)
        call check_carg('dryrun',         self%dryrun)
        call check_carg('dynlp',          self%dynlp)
        call check_carg('dyncls',         self%dyncls)
        call check_carg('eo',             self%eo)
        call check_carg('errify',         self%errify)
        call check_carg('even',           self%even)
        call check_carg('exp_doc',        self%exp_doc)
        call check_carg('extrmode',       self%extrmode)
        call check_carg('fbody',          self%fbody)
        call check_carg('for3D',          self%for3D)
        call check_carg('ft2img',         self%ft2img)
        call check_carg('guinier',        self%guinier)
        call check_carg('hfun',           self%hfun)
        call check_carg('hist',           self%hist)
        call check_carg('imgkind',        self%imgkind)
        call check_carg('kmeans',         self%kmeans)
        call check_carg('label',          self%label)
        call check_carg('local',          self%local)
        call check_carg('masscen',        self%masscen)
        call check_carg('match_filt',     self%match_filt)
        call check_carg('merge',          self%merge)
        call check_carg('mirr',           self%mirr)
        call check_carg('mkdir',          self%mkdir)
        call check_carg('msktype',        self%msktype)
        call check_carg('neg',            self%neg)
        call check_carg('neigh',          self%neigh)
        call check_carg('noise_norm',     self%noise_norm)
        ! call check_carg('noise',          self%noise)
        call check_carg('norec',          self%norec)
        call check_carg('norm',           self%norm)
        call check_carg('objfun',         self%objfun)
        call check_carg('opt',            self%opt)
        call check_carg('order',          self%order)
        call check_carg('oritype',        self%oritype)
        call check_carg('outside',        self%outside)
        call check_carg('pad',            self%pad)
        call check_carg('pcontrast',      self%pcontrast)
        call check_carg('pgrp',           self%pgrp)
        call check_carg('pgrp_known',     self%pgrp_known)
        call check_carg('phaseplate',     self%phaseplate)
        call check_carg('phrand',         self%phrand)
        call check_carg('phshiftunit',    self%phshiftunit)
        call check_carg('prg',            self%prg)
        call check_carg('projname',       self%projname)
        call check_carg('projstats',      self%projstats)
        call check_carg('readwrite',      self%readwrite)
        call check_carg('real_filter',    self%real_filter)
        call check_carg('refine',         self%refine)
        call check_carg('refs',           self%refs)
        call check_carg('remap_cls',      self%remap_cls)
        call check_carg('restart',        self%restart)
        call check_carg('rm_outliers',    self%rm_outliers)
        call check_carg('rnd',            self%rnd)
        call check_carg('roalgn',         self%roalgn)
        call check_carg('round',          self%round)
        call check_carg('shalgn',         self%shalgn)
        call check_carg('shbarrier',      self%shbarrier)
        call check_carg('shellnorm',      self%shellnorm)
        call check_carg('soften',         self%soften)
        call check_carg('speckind',       self%speckind)
        call check_carg('stats',          self%stats)
        call check_carg('stream',         self%stream)
        call check_carg('symrnd',         self%symrnd)
        call check_carg('swap',           self%swap)
        call check_carg('taper_edges',    self%taper_edges)
        call check_carg('test',           self%test)
        call check_carg('time',           self%time)
        call check_carg('tomo',           self%tomo)
        call check_carg('tomoseries',     self%tomoseries)
        call check_carg('trsstats',       self%trsstats)
        call check_carg('tseries',        self%tseries)
        call check_carg('vis',            self%vis)
        call check_carg('wfun',           self%wfun)
        call check_carg('weights2D',      self%weights2D)
        call check_carg('weights3D',      self%weights3D)
        call check_carg('zero',           self%zero)
        ! File args
        call check_file('boxfile',        self%boxfile,      'T')
        call check_file('boxtab',         self%boxtab,       'T')
        call check_file('classdoc',       self%classdoc,     'T')
        call check_file('ctf_estimate_doc', self%ctf_estimate_doc,  'T', 'O')
        call check_file('comlindoc',      self%comlindoc,    'T')
        call check_file('deftab',         self%deftab,       'T', 'O')
        call check_file('doclist',        self%doclist,      'T')
        call check_file('ext',            self%ext,          notAllowed='T')
        call check_file('filetab',        self%filetab,      'T')
        call check_file('fname',          self%fname)
        call check_file('frcs',           self%frcs,         'B')
        call check_file('fsc',            self%fsc,          'B')
        call check_file('infile',         self%infile)
        call check_file('mskfile',        self%mskfile,      notAllowed='T')
        call check_file('oritab',         self%oritab,       'T', 'O')
        call check_file('oritab2',        self%oritab2,      'T', 'O')
        call check_file('oritab3D',       self%oritab3D,     'T', 'O')
        call check_file('outfile',        self%outfile,      'T', 'O')
        call check_file('outstk',         self%outstk,       notAllowed='T')
        call check_file('outstk2',        self%outstk2,      notAllowed='T')
        call check_file('outvol',         self%outvol,       notAllowed='T')
        call check_file('pdbfile',        self%pdbfile)
        call check_file('plaintexttab',   self%plaintexttab, 'T')
        call check_file('projfile',       self%projfile,     'O')
        call check_file('projfile_target',self%projfile_target,'O')
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
        call check_iarg('box_extract',    self%box_extract)
        call check_iarg('boxconvsz',      self%boxconvsz)
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
        call check_iarg('nstates',        self%nstates)
        call check_iarg('class',          self%class)
        call check_iarg('nparts',         self%nparts)
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
        call check_iarg('ring1',          self%ring1)
        call check_iarg('ring2',          self%ring2)
        call check_iarg('startit',        self%startit)
        call check_iarg('state',          self%state)
        call check_iarg('state2split',    self%state2split)
        call check_iarg('stepsz',         self%stepsz)
        call check_iarg('szsn',           self%szsn)
        call check_iarg('time_inactive',  self%time_inactive)
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
        call check_rarg('focusmsk',       self%focusmsk)
        call check_rarg('frac',           self%frac)
        call check_rarg('fraca',          self%fraca)
        call check_rarg('fracdeadhot',    self%fracdeadhot)
        call check_rarg('frac_outliers',  self%frac_outliers)
        call check_rarg('fraczero',       self%fraczero)
        call check_rarg('ftol',           self%ftol)
        call check_rarg('gw',             self%gw)
        call check_rarg('hp',             self%hp)
        call check_rarg('hp_ctf_estimate',self%hp_ctf_estimate)
        call check_rarg('hp_fsc',         self%hp_fsc)
        call check_rarg('inner',          self%inner)
        call check_rarg('kv',             self%kv)
        call check_rarg('lam',            self%lam)
        call check_rarg('lp',             self%lp)
        call check_rarg('lp_ctf_estimate',self%lp_ctf_estimate)
        call check_rarg('lp_pick',        self%lp_pick)
        call check_rarg('lplim_crit',     self%lplim_crit)
        call check_rarg('lpstart',        self%lpstart)
        call check_rarg('lpstop',         self%lpstop)
        call check_rarg('lpvalid',        self%lpvalid)
        call check_rarg('moldiam',        self%moldiam)
        call check_rarg('msk',            self%msk)
        call check_rarg('mul',            self%mul)
        call check_rarg('mw',             self%mw)
        call check_rarg('ndev',           self%ndev)
        call check_rarg('nsig',           self%nsig)
        call check_rarg('outer',          self%outer)
        call check_rarg('phranlp',        self%phranlp)
        call check_rarg('power',          self%power)
        call check_rarg('scale',          self%scale)
        call check_rarg('scale2',         self%scale2)
        call check_rarg('sherr',          self%sherr)
        call check_rarg('smpd',           self%smpd)
        call check_rarg('snr',            self%snr)
        call check_rarg('tau',            self%tau)
        call check_rarg('thres',          self%thres)
        call check_rarg('time_per_image', self%time_per_image)
        call check_rarg('trs',            self%trs)
        call check_rarg('motion_correctftol', self%motion_correctftol)
        call check_rarg('motion_correctgtol', self%motion_correctgtol)
        call check_rarg('update_frac',    self%update_frac)
        call check_rarg('width',          self%width)
        call check_rarg('winsz',          self%winsz)
        call check_rarg('xsh',            self%xsh)
        call check_rarg('ysh',            self%ysh)
        call check_rarg('zsh',            self%zsh)

!>>> START, EXECUTION RELATED
        ! get cwd
        call simple_getcwd(cwd_tmp)
        self%cwd = trim(cwd_tmp)
        ! get process ID
        self%pid = get_process_id()
        ! get name of of executable
        call getarg(0,self%executable)
        ! get pointer to program user interface
        call get_prg_ptr(self%prg, self%ptr2prg)
        if( self%ptr2prg%requires_sp_project() )then
            ! sp_project file required, go look for it
            sp_files  = simple_list_files(glob='*.simple')
            if( allocated(sp_files) )then
                nsp_files = size(sp_files)
            else
                nsp_files = 0
            endif
            select case(nsp_files)
                case(0)
                    write(*,*) 'program: ', trim(self%prg), ' requires a project file!'
                    write(*,*) 'cwd:     ', trim(self%cwd)
                    stop 'ERROR! no *.simple project file identified; simple_params :: new'
                case(1)
                    ! good, we found a single monolithic project file
                    ! set projfile and projname fields
                    self%projfile = trim(sp_files(1))
                    self%projname = basename(self%projfile)
                case DEFAULT
                    write(*,*) 'Multiple *simple project files detected in ', trim(self%cwd)
                    do i=1,nsp_files
                        write(*,*) trim(sp_files(i))
                    end do
                    stop 'ERROR! a unique *.simple project could NOT be identified; simple_params :: new'
            end select
        endif
        L_MKDIR_EXEC = .false.
        if( self%mkdir .eq. 'yes' )then
            if( .not. str_has_substr(self%executable,'private') )then
                ! get next directory number
                idir          = find_next_int_dir_prefix(self%cwd)
                self%exec_dir = int2str(idir)//'_'//trim(self%prg)
                ! make execution directory
                call simple_mkdir('./'//trim(self%exec_dir))
                ! change to execution directory directory
                call simple_chdir('./'//trim(self%projname))
                if( self%ptr2prg%requires_sp_project() )then
                    ! copy the project file from upstairs
                    call syslib_copy_file('../'//trim(self%projfile), './'//trim(self%projfile))
                    ! cwd of SP-project will be update in the builder
                endif
                ! update cwd and set CWD_ORIGINAL
                CWD_ORIGINAL = self%cwd
                ! get new cwd
                call simple_getcwd(cwd_tmp)
                self%cwd = trim(cwd_tmp)
                ! flag what we are up to
                L_MKDIR_EXEC = .true.
            endif
        endif
!>>> END, EXECUTION RELATED

!>>> START, SANITY CHECKING AND PARAMETER EXTRACTION FROM ORITAB(S)/VOL(S)/STACK(S)
        ! put ctf_estimate_doc (if defined) as oritab
        if( cline%defined('ctf_estimate_doc') )then
            if( .not. cline%defined('oritab') )then
                call cline%set('oritab', self%ctf_estimate_doc)
                self%oritab = self%ctf_estimate_doc
            else
                write(*,*) 'WARNING! Could not set ctf_estimate_doc to oritab because oritab is defined'
            endif
        endif
        ! array allocation
        allocate(vol_defined(self%nstates), source=.false.)
        allocate(self%vols(self%nstates), self%vols_even(self%nstates), self%vols_odd(self%nstates), self%mskvols(self%nstates))
        ! determines whether at least one volume is on the cmdline
        do i=1,self%nstates
            if( cline%defined( trim('vol')//int2str(i) )) vol_defined(i) = .true.
        enddo
        ! check inputted vols
        if( cline%defined('vollist') )then
            if( nlines(self%vollist)< MAXS )then
                call read_vols
            endif
        else
            if( any(vol_defined) )then
                do i=1,self%nstates
                    call check_vol( i )
                end do
            endif
        endif
        ! check inputted mask vols
        if( cline%defined('msklist') )then
            if( nlines(self%msklist)< MAXS ) call read_masks
        endif
        ! no stack given, get ldim from volume if present
        if( self%stk .eq. '' .and. allocated(self%vols(1)%str) )then
            call find_ldim_nptcls(self%vols(1)%str, self%ldim, ifoo)
            self%box  = self%ldim(1)
            DebugPrint 'found logical dimension of volume: ', self%ldim
        endif
        ! no stack given,not vol given, get ldim from mskfile if present
        if( self%stk .eq. '' .and. .not. allocated(self%vols(1)%str) .and. self%mskfile .ne. '' )then
            call find_ldim_nptcls(self%mskfile, self%ldim, ifoo)
            self%box  = self%ldim(1)
            DebugPrint 'found logical dimension of volume: ', self%ldim
        endif
        if( cline%defined('oritype') )then ! oritype overrides any other spproj_a_seg setter
            ! this determines the spproj_a_seg
            if( present(spproj_a_seg) ) spproj_a_seg_inputted = self%spproj_a_seg
            select case(trim(self%oritype))
                case('mic')
                    self%spproj_a_seg = MIC_SEG
                case('stk')
                    self%spproj_a_seg = STK_SEG
                case('ptcl2D')
                    self%spproj_a_seg = PTCL2D_SEG
                case('cls2D')
                    self%spproj_a_seg = CLS2D_SEG
                case('cls3D')
                    self%spproj_a_seg = CLS3D_SEG
                case('ptcl3D')
                    self%spproj_a_seg = PTCL3D_SEG
                case('out')
                    self%spproj_a_seg = OUT_SEG
                case('projinfo')
                    self%spproj_a_seg = PROJINFO_SEG
                case('jobproc')
                    self%spproj_a_seg = JOBPROC_SEG
                case('compenv')
                    self%spproj_a_seg = COMPENV_SEG
                case DEFAULT
                    write(*,*) 'oritype: ', trim(self%oritype)
                    stop 'unsupported oritype; simple_params :: new'
            end select
            if( present(spproj_a_seg) )then
                if( spproj_a_seg_inputted .ne. self%spproj_a_seg )then
                    write(*,*) 'oritype: ', trim(self%oritype)
                    write(*,*) 'inputted spproj_a_seg: ', spproj_a_seg_inputted
                    stop 'ERROR! oritype and inputted spproj_a_seg inconsistent'
                endif
            endif
        else
            select case(self%spproj_a_seg)
                case(MIC_SEG)
                    self%oritype = 'mic'
                case(STK_SEG)
                    self%oritype = 'stk'
                case(PTCL2D_SEG)
                    self%oritype = 'ptcl2D'
                case(CLS2D_SEG)
                    self%oritype = 'cls2D'
                case(CLS3D_SEG)
                    self%oritype = 'cls3D'
                case(PTCL3D_SEG)
                    self%oritype = 'ptcl3D'
                case(PROJINFO_SEG)
                    self%oritype = 'projinfo'
                case(JOBPROC_SEG)
                    self%oritype = 'jobproc'
                case(COMPENV_SEG)
                    self%oritype = 'compenv'
                case DEFAULT
                    stop 'unknown spproj_a_seg index; simple_params :: new'
            end select
        endif
        ! take care of nptcls etc.
        ! project is in charge
        if( cline%defined('projname') )then
            self%projfile = trim(self%projname)//'.simple'
            call cline%set('projfile', trim(self%projfile))
        endif
        if( file_exists(trim(self%projfile)) )then
            ! get nptcls/box/smpd from project file
            call bos%open(trim(self%projfile))
            if( self%stream.eq.'no' )then
                if( .not. cline%defined('nptcls') ) self%nptcls = bos%get_n_records(self%spproj_a_seg)
                call o%new_ori_clean
                select case(self%spproj_a_seg)
                    case(MIC_SEG)
                        call bos%read_first_segment_record(MIC_SEG, o)
                    case(STK_SEG, PTCL2D_SEG, PTCL3D_SEG, CLS2D_SEG)
                        call bos%read_first_segment_record(STK_SEG, o)
                    case(OUT_SEG)
                        call bos%read_first_segment_record(OUT_SEG, o)
                end select
                if( o%isthere('smpd') .and. .not. cline%defined('smpd') ) self%smpd = o%get('smpd')
                call bos%read_first_segment_record(STK_SEG, o)
                if( o%isthere('box') .and. .not. cline%defined('box')   ) self%box = nint(o%get('box'))
                call o%kill
            else
                ! nothing to do for streaming, values set at runtime
            endif
            ! CTF plan
            select case(trim(self%oritype))
                case('ptcl2D', 'ptcl3D')
                    call bos%read_first_segment_record(STK_SEG, o)
            end select
            if( .not. cline%defined('ctf') )then
                self%tfplan%flag = 'no'
                if( o%exists() )then
                    if( o%isthere('ctf') )then
                        call o%getter('ctf', ctfflag)
                        self%tfplan%flag = trim(ctfflag)
                    endif
                endif
                self%ctf = trim(self%tfplan%flag)
            else
                self%tfplan%flag = self%ctf
            endif
            if( .not. cline%defined('phaseplate') )then
                if( o%isthere('phaseplate') )then
                    call o%getter('phaseplate', phaseplate)
                    self%tfplan%l_phaseplate = trim(phaseplate).eq.'yes'
                else
                    self%tfplan%l_phaseplate = .false.
                endif
            else
                self%tfplan%l_phaseplate = trim(self%phaseplate).eq.'yes'
            endif
            call bos%close
        else if( self%stk .ne. '' )then
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
            if( .not. cline%defined('nptcls') .and. str_has_substr(self%oritab, '.txt') )then
                ! get # particles from oritab
                self%nptcls = nlines(self%oritab)
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
        ! set refs_even and refs_odd
        if( cline%defined('refs') )then
            self%refs_even = add2fbody(self%refs, self%ext, '_even')
            self%refs_odd  = add2fbody(self%refs, self%ext, '_odd')
        endif
        ! set vols_even and vols_odd
        if( cline%defined('vol1') )then
            do istate=1,self%nstates
                self%vols_even(istate)%str = add2fbody(self%vols(istate)%str, self%ext, '_even')
                self%vols_odd(istate)%str  = add2fbody(self%vols(istate)%str, self%ext, '_odd' )
            end do
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
        if( ddel_scaled )then
            ! delete possibly pre-existing scaled stack parts
            call del_files(trim(STKPARTFBODY), self%nparts, ext=self%ext, suffix='_sc')
        endif
        ! set name of partial files in parallel execution
        stk_part_fname    = trim(STKPARTFBODY)//int2str_pad(self%part,self%numlen)//self%ext
        stk_part_fname_sc = add2fbody(stk_part_fname, self%ext, '_sc')
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
                    call simple_stop('dim mismatch; simple_params :: new')
                endif
            endif
        endif
        ! fractional search and volume update
        if( self%update_frac <= .99)then
            if( self%update_frac < 0.01 )stop 'UPDATE_FRAC is too small 1; simple_params :: constructor'
            if( nint(self%update_frac*real(self%nptcls)) < 1 )&
                &stop 'UPDATE_FRAC is too small 2; simple_params :: constructor'
            self%l_frac_update = .true.
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
        self%hpind_fsc = 0                                         ! high-pass Fouirer index FSC
        if( cline%defined('hp_fsc') ) self%hpind_fsc = nint(self%dstep/self%hp_fsc)
        if( cline%defined('lp') ) self%dynlp = 'no'                ! override dynlp=yes and lpstop
        ! set 2D low-pass limits and smpd_targets 4 scaling
        self%lplims2D(1)       = max(self%fny, self%lpstart)
        self%lplims2D(2)       = max(self%fny, self%lplims2D(1) - (self%lpstart - self%lpstop)/2.)
        self%lplims2D(3)       = max(self%fny, self%lpstop)
        self%smpd_targets2D(1) = self%lplims2D(2)*LP2SMPDFAC2D
        self%smpd_targets2D(2) = self%lplims2D(3)*LP2SMPDFAC2D
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
        ! boxmatch
        self%boxmatch = find_boxmatch(self%box, self%msk)
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
        if( cline%defined('mskfile') )then
            if( .not. file_exists(self%mskfile) )then
                write(*,*) 'file: ', trim(self%mskfile)
                call simple_stop('input mask file not in cwd')
            endif
        endif
        ! focused masking
        self%l_focusmsk = .false.
        if( cline%defined('focusmsk') )then
            if( .not.cline%defined('mskfile') )stop 'mskfile must be provided together with focusmsk'
            if( .not.cline%defined('msk') )then
                stop 'msk must be provided together with focusmsk'
            else
                if(self%focusmsk >= self%msk)stop 'focusmsk should be smaller than msk'
            endif
            self%l_focusmsk = .true.
        else
            self%focusmsk = self%msk
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
        ! box on the command line overrides all other ldim setters
        if( cline%defined('box')  ) self%ldim = [self%box,self%box,1]
        ! if CTF refinement
        if( self%ctf .eq. 'refine' ) call simple_stop('CTF refinement not yet implemented')
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
                if( ncls /= self%ncls ) call simple_stop('inputtend number of clusters (ncls) not&
                &consistent with the number of references in stack (p%refs)')
            else
                self%ncls = ncls
            endif
        endif
        ! set remap_clusters flag
        self%l_remap_cls = .false.
        if( self%remap_cls .eq. 'yes' ) self%l_remap_cls = .true.
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
            select case(trim(self%refine))
            case('single', 'multi', 'snhc')
                    self%trs = 0.
                case DEFAULT
                    self%trs = 1.
            end select
        endif
        self%l_doshift = .true.
        if( self%trs < 0.5 )then
            self%trs = 0.00001
            self%l_doshift = .false.
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
                call simple_stop('need dose_rate to be part of the command line! simple_params :: new')
            endif
            if( .not. cline%defined('exp_time') .and. cline%defined('dose_rate') )then
                call simple_stop('need exp_time to be part of the command line! simple_params :: new')
            endif
            self%l_dose_weight = .true.
        endif
        ! prepare CTF plan
        if( .not.cline%defined('projfile') .and. .not.cline%defined('projname') )then
            self%tfplan%flag         = self%ctf
            self%tfplan%l_phaseplate = .false.
            if( self%phaseplate .eq. 'yes' ) self%tfplan%l_phaseplate = .true.
        endif
        ! set default nnn value to 10% of the search space
        if( .not. cline%defined('nnn') )then
            self%nnn = nint(0.1 * real(self%nspace))
        endif
        ! objective function used in prime2D/3D
        select case(trim(self%objfun))
            case('cc')
                self%l_cc_objfun = .true.
            case('ccres')
                self%l_cc_objfun = .false.
            case DEFAULT
                write(*,*) 'objfun flag: ', trim(self%objfun)
                stop 'unsupported objective function; params :: new'
        end select
        ! B-factor weighted corr or not
        self%l_cc_bfac = .false.
        if( cline%defined('bfac') ) self%l_cc_bfac = .true.
        ! global dev (development) flag
        self%l_dev = .false.
        if( trim(self%dev) .eq. 'yes' ) self%l_dev = .true.
        ! sanity check imgkind
        select case(trim(self%imgkind))
            case('movie', 'mic','ptcl', 'cavg')
                ! alles gut!
            case DEFAULT
                write(*,*) 'imgkind: ', trim(self%imgkind)
                stop 'unsupported imgkind; params :: new'
        end select
!>>> END, IMAGE-PROCESSING-RELATED

        write(*,'(A)') '>>> DONE PROCESSING PARAMETERS'

        contains

            subroutine check_vol( i )
                integer, intent(in)   :: i
                character(len=STDLEN) :: nam
                nam = 'vol'//int2str(i)
                if( cline%defined(nam) )then
                    call check_file(nam, self%vols(i)%str, notAllowed='T')
                    if( .not. file_exists(self%vols(i)%str) )then
                        write(*,*) 'Input volume:', self%vols(i)%str, 'does not exist!'
                        stop
                    endif
                    DebugPrint nam, '=', self%vols(i)%str
                endif
            end subroutine check_vol

            subroutine read_vols
                character(len=STDLEN) :: filenam, nam
                integer :: nl, fnr, i, io_stat
                filenam = cline%get_carg('vollist')
                nl      = nlines(filenam)
                call fopen(fnr, file=filenam, iostat=io_stat)
                if(io_stat /= 0) call fileiochk("params ; read_vols error opening "//trim(filenam), io_stat)
                do i=1,nl
                    read(fnr,*, iostat=io_stat) nam
                    if(io_stat /= 0) call fileiochk("params ; read_vols error reading "//trim(filenam), io_stat)
                    if( nam .ne. '' ) self%vols(i)%str = trim(nam)
                end do
                call fclose(fnr,errmsg="params ; read_vols error closing "//trim(filenam))
            end subroutine read_vols

            subroutine read_masks
                character(len=STDLEN) :: filenam, nam
                integer :: nl, fnr, i, io_stat
                filenam = cline%get_carg('msklist')
                nl      = nlines(filenam)
                call fopen(fnr, file=filenam, iostat=io_stat)
                if(io_stat /= 0) call fileiochk("params ; read_masks error opening "//trim(filenam), io_stat)
                do i=1,nl
                    read(fnr,*, iostat=io_stat) nam
                    if(io_stat /= 0) call fileiochk("params ; read_masks error reading "//trim(filenam), io_stat)
                    if( nam .ne. '' ) self%mskvols(i)%str = trim(nam)
                end do
                call fclose(fnr,errmsg="params ; read_masks error closing "//trim(filenam))
            end subroutine read_masks

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
                        write(*,*) 'File:', trim(file)
                        stop
                    endif
                    select case(file_descr)
                        case ('I')
                            call simple_stop('Support for IMAGIC files is not yet implemented!')
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
                            call simple_stop('This file format is not supported by SIMPLE; simple_params::check_file')
                        case ('T','B','P','O')
                            ! text files are supported
                            ! binary files are supported
                            ! PDB files are supported
                            ! *.simple project files are supported
                        case DEFAULT
                            write(*,*) 'file: ', trim(file)
                            call simple_stop('This file format is not supported by SIMPLE; simple_params::check_file')
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
                                    call simple_stop('The input file names have nonconforming format (mixed formats not allowed)')
                                endif
                            endif
                        end do
                    end do
                    select case(checkupfile(1))
                        case('M','D','B')
                            self%ext = '.mrc'
                        case('S')
                            self%ext = '.spi'
                        case DEFAULT
                            call simple_stop('This file format is not supported by SIMPLE; check_file_formats; simple_params')
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
                        call fileiochk("In params:: double_check_file_formats fopen failed "//trim(self%filetab) , io_stat)
                        read(funit,'(a256)') fname
                        form = fname2format(fname)
                        call fclose(funit, &
                            errmsg="In params:: double_check_file_formats fclose failed "//trim(self%filetab) )
                        select case(form)
                        case('M','D','B')
                                self%ext = '.mrc'
                            case('S')
                                self%ext = '.spi'
                            case DEFAULT
                                write(*,*) 'format string is ', form
                                call simple_stop('File format is not supported by SIMPLE; double_check_file_formats; simple_params')
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
                character(len=*),              intent(in)    :: carg
                character(len=:), allocatable, intent(inout) :: var
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
                        call simple_stop(" In simple_params set_ldim_box_from_stk")
                    endif
                endif
            end subroutine set_ldim_box_from_stk

    end subroutine new

end module simple_params
