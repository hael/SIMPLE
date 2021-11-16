! provides global distribution of constants and derived constants
module simple_parameters
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_user_interface, only: simple_program, get_prg_ptr
use simple_atoms,          only: atoms
!$ use omp_lib
!$ use omp_lib_kinds
implicit none

public :: parameters, params_glob
private
#include "simple_local_flags.inc"

type :: parameters
    ! pointer 2 program UI
    type(simple_program), pointer :: ptr2prg => null()
    ! yes/no decision variables in ascending alphabetical order
    character(len=3)      :: acf='no'             !< calculate autocorrelation function(yes|no){no}
    character(len=3)      :: anneal='yes'         !< use annealing or not
    character(len=3)      :: append='no'          !< append in context of files(yes|no){no}
    character(len=3)      :: async='no'           !< asynchronous (yes|no){no}
    character(len=3)      :: autoscale='no'       !< automatic down-scaling(yes|no){yes}
    character(len=3)      :: avg='no'             !< calculate average (yes|no){no}
    character(len=3)      :: bin='no'             !< binarise image(yes|no){no}
    character(len=3)      :: center='yes'         !< center image(s)/class average(s)/volume(s)(yes|no){no}
    character(len=3)      :: classtats='no'       !< calculate class population statistics(yes|no){no}
    character(len=3)      :: clustvalid='no'      !< validate clustering(yes|homo|no){no}
    character(len=3)      :: compare='no'         !< do comparison(yes|no){no}
    character(len=3)      :: continue='no'        !< continue previous refinement(yes|no){no}
    character(len=3)      :: countvox='no'        !< count # voxels(yes|no){no}
    character(len=3)      :: ctfstats='no'        !< calculate ctf statistics(yes|no){no}
    character(len=3)      :: ctfpatch='yes'       !< whether to perform patched CTF estimation(yes|no){yes}
    character(len=3)      :: cure='no'
    character(len=3)      :: dev='no'             !< development flag for experimental code(yes|no){no}
    character(len=3)      :: dihedral='no'        !< dihedral symmetry or not(yes|no){no}
    character(len=3)      :: discrete='no'        !< be discrete(yes|no){no}
    character(len=3)      :: diverse='no'         !< diverse or not flag (yes|no){no}
    character(len=3)      :: doalign='yes'
    character(len=3)      :: dodock='no'          !< perform atomic position registration (yes|no){no}
    character(len=3)      :: dopca='yes'
    character(len=3)      :: doprint='no'
    character(len=3)      :: dorec='yes'
    character(len=3)      :: elongated='no'       !< elongated particles(yes|no){no}
    character(len=3)      :: even='no'            !< even orientation distribution(yes|no){no}
    character(len=3)      :: fill_holes='no'      !< fill the holes post binarisation(yes|no){no}
    character(len=3)      :: ft2img='no'          !< convert Fourier transform to real image of power(yes|no){no}
    character(len=3)      :: for3D='yes'          !< for 3D analysis(yes|no){yes}
    character(len=3)      :: guinier='no'         !< calculate Guinier plot(yes|no){no}
    character(len=3)      :: graphene_filt='no'   !< filter out graphene bands in correcation search
    character(len=3)      :: gridding='no'         !< to test gridding correction
    character(len=3)      :: groupframes='no'     !< Whether to perform weighted frames averaging during motion correction(yes|no){no}
    character(len=3)      :: keepvol='no'         !< dev flag for preserving iterative volumes in refine3d
    character(len=3)      :: kmeans='yes'
    character(len=3)      :: local='no'
    character(len=3)      :: locres='no'          !< filter based on local resolution or not(yes|no){no}
    character(len=3)      :: makemovie='no'
    character(len=3)      :: masscen='no'         !< center to center of gravity(yes|no){no}
    character(len=3)      :: match_filt='yes'     !< matched filter on (yes|no){yes}
    character(len=3)      :: mcpatch='yes'        !< whether to perform patch-based alignment during motion correction
    character(len=3)      :: merge='no'
    character(len=3)      :: mirr='no'            !< mirror(no|x|y){no}
    character(len=3)      :: mkdir='no'           !< make auto-named execution directory(yes|no){no}
    character(len=3)      :: envfsc='yes'         !< envelope mask even/odd pairs for FSC calculation(yes|no){yes}
    character(len=3)      :: needs_sigma='no'     !< invert contrast of images(yes|no)
    character(len=3)      :: neg='no'             !< invert contrast of images(yes|no)
    character(len=3)      :: noise_norm ='no'
    character(len=3)      :: norm='no'            !< do statistical normalisation avg
    character(len=3)      :: order='no'           !< order ptcls according to correlation(yes|no){no}
    character(len=3)      :: outside='no'         !< extract boxes outside the micrograph boundaries(yes|no){no}
    character(len=3)      :: pad='no'
    character(len=3)      :: phaseplate='no'      !< images obtained with Volta phaseplate(yes|no){no}
    character(len=3)      :: phrand='no'          !< phase randomize(yes|no){no}
    character(len=3)      :: platonic='yes'       !< platonic symmetry or not(yes|no){yes}
    character(len=3)      :: plot='no'            !< make plot(yes|no){no}
    character(len=3)      :: proj_is_class='no'   !< intepret projection directions as classes
    character(len=3)      :: projstats='no'
    character(len=3)      :: roavg='no'           !< rotationally average images in stack
    character(len=3)      :: clsfrcs='no'
    character(len=3)      :: readwrite='no'
    character(len=3)      :: remap_cls='no'
    character(len=3)      :: restart='no'
    character(len=3)      :: rnd='no'             !< random(yes|no){no}
    character(len=3)      :: roalgn='no'
    character(len=3)      :: round='no'
    character(len=3)      :: shalgn='no'          !< do 2D shift alignment(yes|no){no}
    character(len=3)      :: shellnorm='no'
    character(len=3)      :: shbarrier='yes'      !< use shift search barrier constraint(yes|no){yes}
    character(len=3)      :: specstats='no'
    character(len=3)      :: stream='no'          !< sream (real time) execution mode(yes|no){no}
    character(len=3)      :: subtr_backgr='no'
    character(len=3)      :: symrnd='no'          !< randomize over symmetry operations(yes|no){no}
    character(len=3)      :: swap='no'
    character(len=3)      :: taper_edges='no'     !< self-explanatory
    character(len=3)      :: test='no'
    character(len=3)      :: tomo='no'            !< tomography mode(yes|no){no}
    character(len=3)      :: tophat='no'          !< tophat filter(yes|no){no}
    character(len=3)      :: time='no'
    character(len=3)      :: trsstats='no'        !< provide origin shift statistics(yes|no){no}
    character(len=3)      :: tseries='no'         !< images represent a time-series(yes|no){no}
    character(len=3)      :: vis='no'             !< visualise(yes|no)
    character(len=3)      :: zero='no'            !< zeroing(yes|no){no}
    ! files & directories strings in ascending alphabetical order
    character(len=LONGSTRLEN) :: boxfile=''           !< file with EMAN particle coordinates(.txt)
    character(len=LONGSTRLEN) :: boxtab=''            !< table (text file) of files with EMAN particle coordinates(.txt)
    character(len=LONGSTRLEN) :: classdoc=''          !< doc with per-class stats(.txt)
    character(len=LONGSTRLEN) :: cwd=''
    character(len=LONGSTRLEN) :: deftab=''            !< file with CTF info(.txt|.simple)
    character(len=LONGSTRLEN) :: dir=''               !< directory
    character(len=LONGSTRLEN) :: dir_movies=''        !< grab mrc mrcs files from here
    character(len=LONGSTRLEN) :: dir_prev=''          !< grab previous projects for streaming
    character(len=LONGSTRLEN) :: dir_refine=''        !< refinement directory
    character(len=LONGSTRLEN) :: dir_reject='rejected'!< move rejected files to here{rejected}
    character(len=LONGSTRLEN) :: dir_select='selected'!< move selected files to here{selected}
    character(len=LONGSTRLEN) :: dir_target=''        !< put output here
    character(len=LONGSTRLEN) :: dir_ptcls=''
    character(len=LONGSTRLEN) :: dir_box=''
    character(len=LONGSTRLEN) :: doclist=''           !< list of oritabs for different states
    character(len=LONGSTRLEN) :: exec_dir='./'        !< auto-named execution directory
    character(len=LONGSTRLEN) :: filetab=''           !< list of files(.txt)
    character(len=LONGSTRLEN) :: fname=''             !< file name
    character(len=LONGSTRLEN) :: frcs=trim(FRCS_FILE) !< binary file with per-class/proj Fourier Ring Correlations(.bin)
    character(len=LONGSTRLEN) :: fsc='fsc_state01.bin'!< binary file with FSC info{fsc_state01.bin}
    character(len=LONGSTRLEN) :: gainref=''           !< gain reference for movie alignment
    character(len=LONGSTRLEN) :: infile=''            !< file with inputs(.txt)
    character(len=LONGSTRLEN) :: infile2=''           !< file with inputs(.txt)
    character(len=LONGSTRLEN) :: mskfile=''           !< maskfile.ext
    character(len=LONGSTRLEN) :: msklist=''           !< table (text file) of mask volume files(.txt)
    character(len=LONGSTRLEN) :: mskvols(MAXS)=''
    character(len=LONGSTRLEN) :: oritab=''            !< table  of orientations(.txt|.simple)
    character(len=LONGSTRLEN) :: oritab2=''           !< 2nd table of orientations(.txt|.simple)
    character(len=LONGSTRLEN) :: oritab3D=''          !< table of 3D orientations(.txt|.simple)
    character(len=LONGSTRLEN) :: outfile=''           !< output document
    character(len=LONGSTRLEN) :: outstk=''            !< output image stack
    character(len=LONGSTRLEN) :: outvol=''            !< output volume{outvol.ext}
    character(len=LONGSTRLEN) :: pdbfile=''           !< PDB file
    character(len=LONGSTRLEN) :: pdbfile2=''          !< PDB file, another one
    character(len=LONGSTRLEN) :: pdbfiles=''          !< list of PDB files
    character(len=LONGSTRLEN) :: pdfile='pdfile.bin'
    character(len=LONGSTRLEN) :: plaintexttab=''      !< plain text file of input parameters
    character(len=LONGSTRLEN) :: projfile=''          !< SIMPLE *.simple project file
    character(len=LONGSTRLEN) :: projfile_target=''   !< another SIMPLE *.simple project file
    character(len=LONGSTRLEN) :: refs=''              !< initial2Dreferences.ext
    character(len=LONGSTRLEN) :: refs_even=''
    character(len=LONGSTRLEN) :: refs_odd=''
    character(len=LONGSTRLEN) :: star_datadir=''      !< STAR-generated data directory
    character(len=LONGSTRLEN) :: starfile=''          !< STAR-formatted EM file (proj.star)
    character(len=LONGSTRLEN) :: star_mic=''          !< STAR-formatted EM file (micrographs.star)
    character(len=LONGSTRLEN) :: star_model=''        !< STAR-formatted EM file (model.star)
    character(len=LONGSTRLEN) :: star_ptcl=''         !< STAR-formatted EM file (data.star)
    character(len=LONGSTRLEN) :: stk=''               !< particle stack with all images(ptcls.ext)
    character(len=LONGSTRLEN) :: stktab=''            !< list of per-micrograph stacks
    character(len=LONGSTRLEN) :: stk2=''              !< 2nd stack(in selection map: selected(cavgs).ext)
    character(len=LONGSTRLEN) :: stk3=''              !< 3d stack (in selection map (cavgs)2selectfrom.ext)
    character(len=LONGSTRLEN) :: stk_backgr=''        !< stack with image for background subtraction
    character(len=LONGSTRLEN) :: vol=''
    character(len=LONGSTRLEN) :: vol_filt=''          !< input filter volume(vol_filt.ext)
    character(len=LONGSTRLEN) :: vollist=''           !< table (text file) of volume files(.txt)
    character(len=LONGSTRLEN) :: vols(MAXS)=''
    character(len=LONGSTRLEN) :: vols_even(MAXS)=''
    character(len=LONGSTRLEN) :: vols_odd(MAXS)=''
    character(len=LONGSTRLEN) :: voltab=''            !< table (text file) of volume files(.txt)
    character(len=LONGSTRLEN) :: voltab2=''           !< 2nd table (text file) of volume files(.txt)
    character(len=LONGSTRLEN) :: xmlloc=''
    ! other character variables in ascending alphabetical order
    character(len=STDLEN) :: algorithm=''         !< algorithm to be used
    character(len=STDLEN) :: cn_type='cn_std'     !< generalised coordination number (cn_gen) or stardard (cn_std)
    character(len=STDLEN) :: angastunit='degrees' !< angle of astigmatism unit (radians|degrees){degrees}
    character(len=4)      :: automatic='no'       !< automatic thres for edge detect (yes|no) {no}
    character(len=4)      :: automsk='no'
    character(len=STDLEN) :: boxtype='eman'
    character(len=STDLEN) :: wcrit = 'no'         !< correlation weighting scheme (softmax|zscore|sum|cen|exp|no){sum}
    character(len=STDLEN) :: clustermode = 'ar'   !< feature used for clustering (ar|dist|ang|maxint){ar}
    character(len=STDLEN) :: ctf='no'             !< ctf flag(yes|no|flip)
    character(len=STDLEN) :: detector='bin'       !< detector for edge detection (sobel|bin|otsu)
    character(len=STDLEN) :: dfunit='microns'     !< defocus unit (A|microns){microns}
    character(len=STDLEN) :: dir_exec=''          !< name of execution directory
    character(len=STDLEN) :: dockmode='rotshift'  !< mode for docking (rot|shift|rotshift)
    character(len=STDLEN) :: draw_color='white'   !< color in which to identify the picked particle
    character(len=STDLEN) :: executable=''        !< name of executable
    character(len=STDLEN) :: exp_doc=''           !< specifying exp_time and dose_rate per tomogram
    character(len=STDLEN) :: startype=''          !< export type for STAR format (micrograph|select|extract|class2d|initmodel|refine3d|post){all}
    character(len=4)      :: element ='    '      !< atom kind
    character(len=4)      :: ext='.mrc'           !< file extension{.mrc}
    character(len=STDLEN) :: fbody=''             !< file body
    character(len=STDLEN) :: filter='no'          !< filter type{no}
    character(len=STDLEN) :: hfun='sigm'          !< function used for normalization(sigm|tanh|lin){sigm}
    character(len=STDLEN) :: hist='corr'          !< give variable for histogram plot
    character(len=STDLEN) :: imgkind='ptcl'       !< type of image(ptcl|cavg|mic|movie){ptcl}
    character(len=STDLEN) :: keys=''
    character(len=STDLEN) :: label='class'        !< discrete label(class|state){class}
    character(len=STDLEN) :: mcconvention='simple'!< which frame of reference convention to use for motion correction(simple|unblur|relion){simple}
    character(len=STDLEN) :: msktype='soft'       !< type of mask(hard|soft){soft}
    character(len=7)      :: objfun='cc'          !< objective function(cc|euclid){cc}
    character(len=STDLEN) :: opt='bfgs'           !< optimiser (bfgs|simplex){bfgs}
    character(len=STDLEN) :: oritype='ptcl3D'     !< SIMPLE project orientation type(stk|ptcl2D|cls2D|cls3D|ptcl3D)
    character(len=STDLEN) :: pcontrast='black'    !< particle contrast(black|white){black}
    character(len=STDLEN) :: pgrp='c1'            !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN) :: pgrp_start='c1'      !< point-group symmetry(cn|dn|t|o|i)
    character(len=STDLEN) :: phshiftunit='radians'!< additional phase-shift unit (radians|degrees){radians}
    character(len=STDLEN) :: picker='phasecorr'   !< picker approach (seg|phasecorr|old_school){phasecorr}
    character(len=STDLEN) :: prg=''               !< SIMPLE program being executed
    character(len=STDLEN) :: projname=''          !< SIMPLE  project name
    character(len=STDLEN) :: ptclw='yes'          !< use particle weights(yes|no){yes}
    character(len=STDLEN) :: qsys_name='local'    !< name of queue system (local|slurm|pbs)
    character(len=STDLEN) :: real_filter=''
    character(len=STDLEN) :: refine='shc'         !< refinement mode(snhc|shc|greedy|neigh|cont|cluster|clustersym){shc}
    character(len=STDLEN) :: speckind='sqrt'      !< power spectrum kind(real|power|sqrt|log|phase){sqrt}
    character(len=STDLEN) :: split_mode='even'
    character(len=STDLEN) :: stats='no'           !< provide statistics(yes|no|print){no}
    character(len=STDLEN) :: stk_part=''
    character(len=STDLEN) :: tomoseries=''        !< filetable of filetables of tomograms
    character(len=STDLEN) :: wfun='kb'
    character(len=:), allocatable :: last_prev_dir !< last previous execution directory
    ! special integer kinds
    integer(kind(ENUM_ORISEG))     :: spproj_iseg  = PTCL3D_SEG    !< sp-project segments that b%a points to
    integer(kind(ENUM_OBJFUN))     :: cc_objfun    = OBJFUN_CC     !< objective function(OBJFUN_CC = 0, OBJFUN_EUCLID = 1)
    integer(kind=kind(ENUM_WCRIT)) :: wcrit_enum   = CORRW_CRIT    !< criterium for correlation-based weights
    ! integer variables in ascending alphabetical order
    integer :: angstep=5
    integer :: avgsz=0
    integer :: batchsz=0
    integer :: balance=0           !< max pop for balancing restraint{0}
    integer :: binwidth=1          !< binary layers grown for molecular envelope(in pixels){1}
    integer :: box=0               !< square image size(in pixels)
    integer :: box_original
    integer :: box_extract
    integer :: boxpd=0
    integer :: chunksz=0           !< # images/orientations in chunk
    integer :: class=1             !< cluster identity
    integer :: clip=0              !< clipped image box size(in pixels)
    integer :: cn=8                !< fixed std coord number for atoms in nanos
    integer :: cn_max=12           !< max std coord number for atoms in nanos
    integer :: cn_min=4            !< min std coord number for atoms in nanos
    integer :: cn_stop=10          !< rotational symmetry order stop index{10}
    integer :: corner=0            !< corner size(in pixels){0}
    integer :: cube=0              !< side size(in pixels){0}
    integer :: edge=6              !< edge size for softening molecular envelope(in pixels)
    integer :: eer_fraction=20     !< # of eer raw frames to fraction together
    integer :: eer_upsampling=1    !< eer up-sampling
    integer :: extr_iter=1
    integer :: find=1              !< Fourier index
    integer :: nframesgrp=0        !< # frames to group before motion_correct(Falcon 3){0}
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
    integer :: kstop=0
    integer :: ldim(3)=0
    integer :: lp_iters=1          !< # iters low-pass limited refinement
    integer :: maxits=500          !< maximum # iterations
    integer :: maxp=0
    integer :: minp=10             !< minimum cluster population
    integer :: mrcmode=2
    integer :: navgs=1
    integer :: nchunks=0          !< # computing units, can be < nparts{nparts}
    integer :: ncunits=0           !< # computing units, can be < nparts{nparts}
    integer :: nboot=0
    integer :: ncls=500            !< # clusters
    integer :: ncls_start=10       !< minimum # clusters for 2D streaming
    integer :: ncomps=0
    integer :: ndiscrete=0         !< # discrete orientations
    integer :: ndocs=0             !< # documents
    integer :: newbox=0            !< new box for scaling (by Fourier padding/clipping)
    integer :: nframes=0           !< # frames{30}
    integer :: ngrow=0             !< # of white pixel layers to grow in binary image
    integer :: nmics=0             !< # micographs
    integer :: nmovies_trial=0     !< # of movies after which preprocess_stream will stop
    integer :: noris=0
    integer :: nparts=1            !< # partitions in distributed execution
    integer :: nparts_chunk=1      !< # partitions in chunks distributed execution
    integer :: npix=0              !< # pixles/voxels in binary representation
    integer :: nptcls=1            !< # images in stk/# orientations in oritab
    integer :: nptcls_per_cls=400  !< # images in stk/# orientations in oritab
    integer :: nran=0              !< # random images to select
    integer :: nrefs=100           !< # references used for picking{100}
    integer :: nrestarts=1
    integer :: nspace=2500         !< # projection directions
    integer :: nstates=1           !< # states to reconstruct
    integer :: nsym=1
    integer :: nthr=1              !< # OpenMP threads{1}
    integer :: nptcls_trial=0      !< # of particles after which preprocess_stream will stop
    integer :: numlen=0            !< length of number string
    integer :: numlen_tomo=3       !< length of number string tomo series index{3}
    integer :: nvalid=0
    integer :: nvars=30
    integer :: nvox=0              !< # voxels{0}
    integer :: nxpatch=MC_NPATCH   !< # of patches along x for motion correctionn{5}
    integer :: nypatch=MC_NPATCH   !< # of patches along y for motion correctionn{5}
    integer :: offset=10           !< pixels offset{10}
    integer :: optics_offset=0
    integer :: part=1
    integer :: pcasz=0
    integer :: pid=0               !< process ID
    integer :: pspecsz=512         !< size of power spectrum(in pixels)
    integer :: ptcl=1
    integer :: recl_cgrid=-1
    integer :: reliongroups=0
    integer :: spec=0
    integer :: startit=1           !< start iterating from here
    integer :: state=1             !< state to extract
    integer :: state2split=0       !< state group to split
    integer :: stepsz=1            !< size of step{0}
    integer :: szsn=SZSN_INIT      !< size of stochastic neighborhood{5}
    integer :: time_inactive=120   !< end time limit(mins)
    integer :: tofny=0
    integer :: top=1
    integer :: tos=1
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
    real    :: astigtol=0.05       !< expected (tolerated) astigmatism(in microns){0.05}
    real    :: athres=7.           !< angular threshold(in degrees)
    real    :: batchfrac=1.0
    real    :: bfac=200            !< bfactor for sharpening/low-pass filtering(in A**2){200.}
    real    :: bfacerr=50.         !< bfactor error in simulated images(in A**2){0}
    real    :: bw_ratio=0.3        !< ratio between foreground-background pixel desired in edge detection
    real    :: cenlp=20.           !< low-pass limit for binarisation in centering(in A){30 A}
    real    :: cs=2.7              !< spherical aberration constant(in mm){2.7}
    real    :: cn_thres=5.         !< threshold (outliers removal based on coordination number)
    real    :: ctfreslim=8.
    real    :: dcrit_rel=0.5       !< critical distance relative to box(0-1){0.5}
    real    :: deflim=4.
    real    :: defocus=2.          !< defocus(in microns){2.}
    real    :: dens=0.
    real    :: dfclose=1.
    real    :: dffar=4.
    real    :: dferr=1.            !< defocus error(in microns){1.0}
    real    :: dfmax=5.0           !< maximum expected defocus(in microns)
    real    :: dfmin=0.3           !< minimum expected defocus(in microns)
    real    :: dfsdev=0.1
    real    :: dose_rate=30.0      !< dose rate(in e/A2/s)
    real    :: dstep=0.
    real    :: dsteppd=0.
    real    :: e1=0.               !< 1st Euler(in degrees){0}
    real    :: e2=0.               !< 2nd Euler(in degrees){0}
    real    :: e3=0.               !< 3d Euler(in degrees){0}
    real    :: eps=0.003           !< learning rate{0.003}
    real    :: eullims(3,2)=0.
    real    :: exp_time=2.0        !< exposure time(in s)
    real    :: extr_init=EXTRINITHRESH !< initial extremal ratio (0-1)
    real    :: fny=0.
    real    :: focusmsk=0.         !< spherical msk for use with focused refinement (radius in pixels)
    real    :: focusmskdiam=0.     !< spherical msk for use with focused refinement (diameter in Angstroms)
    real    :: frac=1.             !< fraction of ptcls(0-1){1}
    real    :: fraca=0.1           !< fraction of amplitude contrast used for fitting CTF{0.1}
    real    :: fracdeadhot=0.05    !< fraction of dead or hot pixels{0.01}
    real    :: frac_outliers=0.
    real    :: fraczero=0.
    real    :: ftol=1e-6
    real    :: motion_correctftol = 1e-6   !< tolerance (gradient) for motion_correct
    real    :: motion_correctgtol = 1e-6   !< tolerance (function value) for motion_correct
    real    :: hp=100.             !< high-pass limit(in A)
    real    :: hp_fsc=0.           !< FSC high-pass limit(in A)
    real    :: hp_ctf_estimate=30. !< high-pass limit 4 ctf_estimate(in A)
    real    :: inner=0.            !< inner mask radius(in pixels)
    real    :: innerdiam=0.        !< inner mask diameter(in A)
    real    :: kv=300.             !< acceleration voltage(in kV){300.}
    real    :: lambda=0.5
    real    :: lp=20.              !< low-pass limit(in A)
    real    :: lp_backgr=20.       !< low-pass for solvent blurring (in A)
    real    :: lp_ctf_estimate=5.0 !< low-pass limit 4 ctf_estimate(in A)
    real    :: lp_pick=20.         !< low-pass limit 4 picker(in A)
    real    :: lplim_crit=0.5      !< corr criterion low-pass limit assignment(0.143-0.5){0.5}
    real    :: lplims2D(3)
    real    :: lpmed=20.
    real    :: lpstart=0.          !< start low-pass limit(in A){15}
    real    :: lpstop=8.0          !< stop low-pass limit(in A){8}
    real    :: lpthresh=30.
    real    :: max_rad=0.          !< particle longest  dim (in pixels)
    real    :: min_rad=100.        !< particle shortest dim (in pixels)
    real    :: moldiam=140.        !< molecular diameter(in A)
    real    :: moment=0.
    real    :: msk=0.              !< mask radius(in pixels)
    real    :: mskdiam=0.          !< mask diameter(in Angstroms)
    real    :: mul=1.              !< origin shift multiplication factor{1}
    real    :: mw=0.               !< molecular weight(in kD)
    real    :: ndev=2.0            !< # deviations in one-cluster clustering
    real    :: nsig=2.5            !< # sigmas
    real    :: phranlp=35.         !< low-pass phase randomize(yes|no){no}
    real    :: power=2.
    real    :: scale=1.            !< image scale factor{1}
    real    :: sherr=0.            !< shift error(in pixels){2}
    real    :: sigma=1.0           !< for gaussian function generation {1.}
    real    :: smpd=2.             !< sampling distance, same as EMANs apix(in A)
    real    :: smpd_targets2D(2)
    real    :: snr=0.              !< signal-to-noise ratio
    real    :: tilt_thres=0.05
    real    :: thres=0.            !< threshold (binarisation: 0-1; distance filer: in pixels)
    real    :: thres_low=0.        !< lower threshold for canny edge detection
    real    :: thres_up=1.         !< upper threshold for canny edge detection
    real    :: tiltgroupmax=0
    real    :: time_per_image=200.
    real    :: time_per_frame=0.
    real    :: trs=0.              !< maximum halfwidth shift(in pixels)
    real    :: update_frac = 1.
    real    :: width=10.           !< falloff of inner mask(in pixels){10}
    real    :: winsz=RECWINSZ
    real    :: xsh=0.              !< x shift(in pixels){0}
    real    :: ysh=0.              !< y shift(in pixels){0}
    real    :: zsh=0.              !< z shift(in pixels){0}
    ! logical variables in (roughly) ascending alphabetical order
    logical :: l_autoscale      = .false.
    logical :: l_corrw          = .false.
    logical :: l_distr_exec     = .false.
    logical :: l_dev            = .false.
    logical :: l_dose_weight    = .false.
    logical :: l_doshift        = .false.
    logical :: l_envfsc         = .false.
    logical :: l_focusmsk       = .false.
    logical :: l_frac_update    = .false.
    logical :: l_graphene       = .false.
    logical :: l_innermsk       = .false.
    logical :: l_lpset          = .false.
    logical :: l_locres         = .false.
    logical :: l_match_filt     = .true.
    logical :: l_needs_sigma    = .false.
    logical :: l_phaseplate     = .false.
    logical :: l_remap_cls      = .false.
    logical :: l_wglob          = .true.
    logical :: sp_required      = .false.
  contains
    procedure          :: new
    procedure, private :: set_img_format
end type parameters

class(parameters), pointer :: params_glob => null()

contains

    subroutine new( self, cline, silent )
        use simple_ori,        only: ori
        use simple_sp_project, only: sp_project
        use simple_binoris,    only: binoris
        class(parameters), target, intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        logical,         optional, intent(in)    :: silent
        character(len=LONGSTRLEN), allocatable   :: sp_files(:)
        character(len=:),          allocatable   :: phaseplate, ctfflag
        logical                       :: vol_defined(MAXS)
        character(len=1)              :: checkupfile(50)
        character(len=:), allocatable :: absname
        type(binoris)    :: bos
        type(sp_project) :: spproj
        type(ori)        :: o
        type(atoms)      :: atoms_obj
        real             :: smpd, mskdiam_default, msk_default
        integer          :: i, ncls, ifoo, lfoo(3), cntfile, istate
        integer          :: idir, nsp_files, box, nptcls, nthr
        logical          :: nparts_set, ssilent
        ssilent = .false.
        if( present(silent) ) ssilent = silent
        ! seed random number generator
        call seed_rnd
        ! constants
        nparts_set = .false.
        ! file counter
        cntfile = 0
        ! default initialisations that depend on meta-data file format
        self%outfile = 'outfile'//trim(METADATA_EXT)
        ! checkers in ascending alphabetical order
        call check_carg('acf',            self%acf)
        call check_carg('algorithm',      self%algorithm)
        call check_carg('anneal',         self%anneal)
        call check_carg('angastunit',     self%angastunit)
        call check_carg('append',         self%append)
        call check_carg('async',          self%async)
        call check_carg('automsk',        self%automsk)
        call check_carg('automatic',      self%automatic)
        call check_carg('autoscale',      self%autoscale)
        call check_carg('avg',            self%avg)
        call check_carg('bin',            self%bin)
        call check_carg('boxtype',        self%boxtype)
        call check_carg('center',         self%center)
        call check_carg('classtats',      self%classtats)
        call check_carg('clustermode',    self%clustermode)
        call check_carg('clustvalid',     self%clustvalid)
        call check_carg('cn_type',        self%cn_type)
        call check_carg('compare',        self%compare)
        call check_carg('continue',       self%continue)
        call check_carg('countvox',       self%countvox)
        call check_carg('ctf',            self%ctf)
        call check_carg('ctfpatch',       self%ctfpatch)
        call check_carg('ctfstats',       self%ctfstats)
        call check_carg('cure',           self%cure)
        call check_carg('detector',       self%detector)
        call check_carg('dfunit',         self%dfunit)
        call check_carg('dir_exec',       self%dir_exec)
        call check_carg('discrete',       self%discrete)
        call check_carg('diverse',        self%diverse)
        call check_carg('doalign',        self%doalign)
        call check_carg('dockmode',       self%dockmode)
        call check_carg('dodock',         self%dodock)
        call check_carg('dorec',          self%dorec)
        call check_carg('dev',            self%dev)
        call check_carg('dihedral',       self%dihedral)
        call check_carg('dopca',          self%dopca)
        call check_carg('doprint',        self%doprint)
        call check_carg('draw_color',     self%draw_color)
        call check_carg('element',        self%element)
        call check_carg('elongated',      self%elongated)
        call check_carg('even',           self%even)
        call check_carg('exp_doc',        self%exp_doc)
        call check_carg('startype',       self%startype)
        call check_carg('fbody',          self%fbody)
        call check_carg('fill_holes',     self%fill_holes)
        call check_carg('filter',         self%filter)
        call check_carg('for3D',          self%for3D)
        call check_carg('groupframes',    self%groupframes)
        call check_carg('ft2img',         self%ft2img)
        call check_carg('guinier',        self%guinier)
        call check_carg('graphene_filt',  self%graphene_filt)
        call check_carg('gridding',        self%gridding)
        call check_carg('hfun',           self%hfun)
        call check_carg('hist',           self%hist)
        call check_carg('imgkind',        self%imgkind)
        call check_carg('keepvol',        self%keepvol)
        call check_carg('keys',           self%keys)
        call check_carg('kmeans',         self%kmeans)
        call check_carg('label',          self%label)
        call check_carg('local',          self%local)
        call check_carg('locres',         self%locres)
        call check_carg('makemovie',      self%makemovie)
        call check_carg('masscen',        self%masscen)
        call check_carg('match_filt',     self%match_filt)
        call check_carg('mcpatch',        self%mcpatch)
        call check_carg('merge',          self%merge)
        call check_carg('mirr',           self%mirr)
        call check_carg('mkdir',          self%mkdir)
        call check_carg('envfsc',         self%envfsc)
        call check_carg('msktype',        self%msktype)
        call check_carg('mcconvention',   self%mcconvention)
        call check_carg('needs_sigma',    self%needs_sigma)
        call check_carg('neg',            self%neg)
        call check_carg('noise_norm',     self%noise_norm)
        call check_carg('norm',           self%norm)
        call check_carg('objfun',         self%objfun)
        call check_carg('opt',            self%opt)
        call check_carg('order',          self%order)
        call check_carg('oritype',        self%oritype)
        call check_carg('outside',        self%outside)
        call check_carg('pad',            self%pad)
        call check_carg('pcontrast',      self%pcontrast)
        call check_carg('pgrp',           self%pgrp)
        call check_carg('pgrp_start',     self%pgrp_start)
        call check_carg('phaseplate',     self%phaseplate)
        call check_carg('phrand',         self%phrand)
        call check_carg('phshiftunit',    self%phshiftunit)
        call check_carg('picker',         self%picker)
        call check_carg('platonic',       self%platonic)
        call check_carg('prg',            self%prg)
        call check_carg('projname',       self%projname)
        call check_carg('proj_is_class',  self%proj_is_class)
        call check_carg('projstats',      self%projstats)
        call check_carg('ptclw',          self%ptclw)
        call check_carg('clsfrcs',        self%clsfrcs)
        call check_carg('qsys_name',      self%qsys_name)
        call check_carg('readwrite',      self%readwrite)
        call check_carg('real_filter',    self%real_filter)
        call check_carg('refine',         self%refine)
        call check_carg('remap_cls',      self%remap_cls)
        call check_carg('restart',        self%restart)
        call check_carg('rnd',            self%rnd)
        call check_carg('roalgn',         self%roalgn)
        call check_carg('roavg',          self%roavg)
        call check_carg('round',          self%round)
        call check_carg('shalgn',         self%shalgn)
        call check_carg('shbarrier',      self%shbarrier)
        call check_carg('shellnorm',      self%shellnorm)
        call check_carg('speckind',       self%speckind)
        call check_carg('specstats',      self%specstats)
        call check_carg('stats',          self%stats)
        call check_carg('stream',         self%stream)
        call check_carg('subtr_backgr',   self%subtr_backgr)
        call check_carg('symrnd',         self%symrnd)
        call check_carg('swap',           self%swap)
        call check_carg('taper_edges',    self%taper_edges)
        call check_carg('test',           self%test)
        call check_carg('time',           self%time)
        call check_carg('tomo',           self%tomo)
        call check_carg('tomoseries',     self%tomoseries)
        call check_carg('tophat',         self%tophat)
        call check_carg('trsstats',       self%trsstats)
        call check_carg('tseries',        self%tseries)
        call check_carg('vis',            self%vis)
        call check_carg('wcrit',          self%wcrit)
        call check_carg('wfun',           self%wfun)
        call check_carg('zero',           self%zero)
        ! File args
        call check_file('boxfile',        self%boxfile,      'T')
        call check_file('boxtab',         self%boxtab,       'T')
        call check_file('classdoc',       self%classdoc,     'T')
        call check_file('deftab',         self%deftab,       'T', 'O')
        call check_file('doclist',        self%doclist,      'T')
        call check_file('ext',            self%ext,          notAllowed='T')
        call check_file('filetab',        self%filetab,      'T')
        call check_file('fname',          self%fname)
        call check_file('frcs',           self%frcs,         'B')
        call check_file('fsc',            self%fsc,          'B')
        call check_file('gainref',        self%gainref)
        call check_file('infile',         self%infile,       'T')
        call check_file('infile2',        self%infile2,      'T')
        call check_file('mskfile',        self%mskfile,      notAllowed='T')
        call check_file('oritab',         self%oritab,       'T', 'O')
        call check_file('oritab2',        self%oritab2,      'T', 'O')
        call check_file('oritab3D',       self%oritab3D,     'T', 'O')
        call check_file('outfile',        self%outfile,      'T', 'O')
        call check_file('outstk',         self%outstk,       notAllowed='T')
        call check_file('outvol',         self%outvol,       notAllowed='T')
        call check_file('pdbfile',        self%pdbfile)
        call check_file('pdbfile2',       self%pdbfile2)
        call check_file('pdbfiles',       self%pdbfiles,     'T')
        call check_file('plaintexttab',   self%plaintexttab, 'T')
        call check_file('projfile',       self%projfile,     'O')
        call check_file('projfile_target',self%projfile_target,'O')
        call check_file('refs',           self%refs,         notAllowed='T')
        call check_file('starfile',       self%starfile,     'R')  ! R for relion, S taken by SPIDER
        call check_file('star_mic',       self%star_mic,     'R')  ! R for relion, S taken by SPIDER
        call check_file('star_model',     self%star_model,   'R')  ! R for relion, S taken by SPIDER
        call check_file('star_ptcl',      self%star_ptcl,    'R')  ! R for relion, S taken by SPIDER
        call check_file('stk',            self%stk,          notAllowed='T')
        call check_file('stktab',         self%stktab,       'T')
        call check_file('stk2',           self%stk2,         notAllowed='T')
        call check_file('stk3',           self%stk3,         notAllowed='T')
        call check_file('stk_backgr',     self%stk_backgr,   notAllowed='T')
        call check_file('vol_filt',       self%vol_filt,     notAllowed='T')
        call check_file('vollist',        self%vollist,      'T')
        call check_file('voltab',         self%voltab,       'T')
        call check_file('voltab2',        self%voltab2,      'T')
        ! Dir args
        call check_dir('dir',            self%dir)
        call check_dir('dir_box',        self%dir_box)
        call check_dir('dir_movies',     self%dir_movies)
        call check_dir('dir_prev',       self%dir_prev)
        call check_dir('dir_ptcls',      self%dir_ptcls)
        call check_dir('dir_refine',     self%dir_refine)
        call check_dir('dir_reject',     self%dir_reject)
        call check_dir('dir_select',     self%dir_select)
        call check_dir('dir_target',     self%dir_target)
        call check_dir('star_datadir',   self%star_datadir)
        ! Integer args
        call check_iarg('angstep',        self%angstep)
        call check_iarg('avgsz',          self%avgsz)
        call check_iarg('balance',        self%balance)
        call check_iarg('binwidth',       self%binwidth)
        call check_iarg('box',            self%box)
        call check_iarg('box_extract',    self%box_extract)
        call check_iarg('chunksz',        self%chunksz)
        call check_iarg('clip',           self%clip)
        call check_iarg('cn',             self%cn)
        call check_iarg('cn_max',         self%cn_max)
        call check_iarg('cn_min',         self%cn_min)
        call check_iarg('cn_stop',        self%cn_stop)
        call check_iarg('corner',         self%corner)
        call check_iarg('cube',           self%cube)
        call check_iarg('edge',           self%edge)
        call check_iarg('eer_fraction',   self%eer_fraction)
        call check_iarg('eer_upsampling', self%eer_upsampling)
        call check_iarg('extr_iter',      self%extr_iter)
        call check_iarg('find',           self%find)
        call check_iarg('nframesgrp',     self%nframesgrp)
        call check_iarg('fromp',          self%fromp)
        call check_iarg('fstep',          self%fstep)
        call check_iarg('grow',           self%grow)
        call check_iarg('iares',          self%iares)
        call check_iarg('ind',            self%ind)
        call check_iarg('jumpsz',         self%jumpsz)
        call check_iarg('lp_iters',       self%lp_iters)
        call check_iarg('maxits',         self%maxits)
        call check_iarg('maxp',           self%maxp)
        call check_iarg('minp',           self%minp)
        call check_iarg('mrcmode',        self%mrcmode)
        call check_iarg('navgs',          self%navgs)
        call check_iarg('nboot',          self%nboot)
        call check_iarg('nchunks',        self%nchunks)
        call check_iarg('ncls',           self%ncls)
        call check_iarg('ncls_start',     self%ncls_start)
        call check_iarg('ncunits',        self%ncunits)
        call check_iarg('ndiscrete',      self%ndiscrete)
        call check_iarg('ndocs',          self%ndocs)
        call check_iarg('newbox',         self%newbox)
        call check_iarg('nframes',        self%nframes)
        call check_iarg('ngrow',          self%ngrow)
        call check_iarg('nmovies_trial',  self%nmovies_trial)
        call check_iarg('noris',          self%noris)
        call check_iarg('nran',           self%nran)
        call check_iarg('nrefs',          self%nrefs)
        call check_iarg('nrestarts',      self%nrestarts)
        call check_iarg('nspace',         self%nspace)
        call check_iarg('nstates',        self%nstates)
        call check_iarg('class',          self%class)
        call check_iarg('nparts',         self%nparts)
        call check_iarg('nparts_chunk',   self%nparts_chunk)
        call check_iarg('npix',           self%npix)
        call check_iarg('nptcls',         self%nptcls)
        call check_iarg('nptcls_per_cls', self%nptcls_per_cls)
        call check_iarg('nptcls_trial',   self%nptcls_trial)
        call check_iarg('nthr',           self%nthr)
        call check_iarg('numlen',         self%numlen)
        call check_iarg('numlen_tomo',    self%numlen_tomo)
        call check_iarg('nvars',          self%nvars)
        call check_iarg('nvox',           self%nvox)
        call check_iarg('nxpatch',        self%nxpatch)
        call check_iarg('nypatch',        self%nypatch)
        call check_iarg('offset',         self%offset)
        call check_iarg('optics_offset',  self%optics_offset)
        call check_iarg('part',           self%part)
        call check_iarg('pspecsz',        self%pspecsz)
        call check_iarg('startit',        self%startit)
        call check_iarg('state',          self%state)
        call check_iarg('state2split',    self%state2split)
        call check_iarg('stepsz',         self%stepsz)
        call check_iarg('szsn',           self%szsn)
        call check_iarg('time_inactive',  self%time_inactive)
        call check_iarg('top',            self%top)
        call check_iarg('tos',            self%tos)
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
        call check_rarg('bw_ratio',       self%bw_ratio)
        call check_rarg('cenlp',          self%cenlp)
        call check_rarg('cs',             self%cs)
        call check_rarg('cn_thres',       self%cn_thres)
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
        call check_rarg('dose_rate',      self%dose_rate)
        call check_rarg('e1',             self%e1)
        call check_rarg('e2',             self%e2)
        call check_rarg('e3',             self%e3)
        call check_rarg('eps',            self%eps)
        call check_rarg('exp_time',       self%exp_time)
        call check_rarg('extr_init',      self%extr_init)
        call check_rarg('focusmsk',       self%focusmsk)
        call check_rarg('focusmskdiam',   self%focusmskdiam)
        call check_rarg('frac',           self%frac)
        call check_rarg('fraca',          self%fraca)
        call check_rarg('fracdeadhot',    self%fracdeadhot)
        call check_rarg('frac_outliers',  self%frac_outliers)
        call check_rarg('fraczero',       self%fraczero)
        call check_rarg('ftol',           self%ftol)
        call check_rarg('hp',             self%hp)
        call check_rarg('hp_ctf_estimate',self%hp_ctf_estimate)
        call check_rarg('hp_fsc',         self%hp_fsc)
        call check_rarg('inner',          self%inner)
        call check_rarg('innerdiam',      self%innerdiam)
        call check_rarg('kv',             self%kv)
        call check_rarg('lambda',         self%lambda)
        call check_rarg('lp',             self%lp)
        call check_rarg('lp_backgr',      self%lp_backgr)
        call check_rarg('lp_ctf_estimate',self%lp_ctf_estimate)
        call check_rarg('lp_pick',        self%lp_pick)
        call check_rarg('lplim_crit',     self%lplim_crit)
        call check_rarg('lpstart',        self%lpstart)
        call check_rarg('lpstop',         self%lpstop)
        call check_rarg('lpthresh',       self%lpthresh)
        call check_rarg('max_rad',        self%max_rad)
        call check_rarg('min_rad',        self%min_rad)
        call check_rarg('moldiam',        self%moldiam)
        call check_rarg('msk',            self%msk)
        call check_rarg('mskdiam',        self%mskdiam)
        call check_rarg('mul',            self%mul)
        call check_rarg('mw',             self%mw)
        call check_rarg('ndev',           self%ndev)
        call check_rarg('nsig',           self%nsig)
        call check_rarg('phranlp',        self%phranlp)
        call check_rarg('power',          self%power)
        call check_rarg('scale',          self%scale)
        call check_rarg('sherr',          self%sherr)
        call check_rarg('smpd',           self%smpd)
        call check_rarg('sigma',          self%sigma)
        call check_rarg('snr',            self%snr)
        call check_rarg('tilt_thres',     self%tilt_thres)
        call check_rarg('thres',          self%thres)
        call check_rarg('thres_low',      self%thres_low)
        call check_rarg('thres_up',       self%thres_up)
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
        call simple_getcwd(self%cwd)
        ! update CWD globals in defs
        if( .not. allocated(cwd_glob_orig) ) allocate(cwd_glob_orig, source=trim(self%cwd))
        cwd_glob = trim(self%cwd)
        ! get process ID
        self%pid = get_process_id()
        ! get name of executable
        call get_command_argument(0,self%executable)
        if(len_trim(self%executable) == 0) THROW_HARD('get_command_argument failed; new')
        ! set execution mode (shmem or distr)
        if( cline%defined('part') .and. cline%defined('nparts') )then
            self%l_distr_exec = .true.
            part_glob = self%part
        else
            self%l_distr_exec = .false.
            part_glob = 1
        endif
        l_distr_exec_glob = self%l_distr_exec
        ! get pointer to program user interface
        call get_prg_ptr(self%prg, self%ptr2prg)
        ! look for the last previous execution directory and get next directory number
        if( allocated(self%last_prev_dir) ) deallocate(self%last_prev_dir)
        idir = find_next_int_dir_prefix(self%cwd, self%last_prev_dir)
        ! look for a project file
        if( .not.cline%defined('projfile') )then
            if( allocated(self%last_prev_dir) )then
                call simple_list_files(self%last_prev_dir//'/*.simple', sp_files)
            else
                call simple_list_files('*.simple', sp_files)
            endif
            if( allocated(sp_files) )then
                nsp_files = size(sp_files)
            else
                nsp_files = 0
            endif
        else
            nsp_files = 0
        endif
        if( associated(self%ptr2prg) .and. .not. str_has_substr(self%executable,'private') )then
            ! the associated(self%ptr2prg) condition required for commanders executed within
            ! distributed workflows without invoking simple_private_exec
            self%sp_required = self%ptr2prg%requires_sp_project()
            if( .not.cline%defined('projfile') .and. self%sp_required )then
                ! so that distributed execution can deal with multiple project files (eg scaling)
                ! and simple_exec can also not require a project file
                if( nsp_files > 1 )then
                    write(logfhandle,*) 'Multiple *simple project files detected'
                    do i=1,nsp_files
                        write(logfhandle,*) trim(sp_files(i))
                    end do
                    THROW_HARD('a unique *.simple project could NOT be identified; new')
                endif
            endif
        else
            self%sp_required = .false.
        endif
        if( nsp_files == 0 .and. self%sp_required .and. .not. cline%defined('projfile') )then
            write(logfhandle,*) 'program: ', trim(self%prg), ' requires a project file!'
            write(logfhandle,*) 'cwd:     ', trim(self%cwd)
            THROW_HARD('no *.simple project file identified; new')
        endif
        if( nsp_files == 1 .and. self%sp_required )then
            ! good, we found a single monolithic project file
            ! set projfile and projname fields unless given on command line
            if( .not. cline%defined('projfile') ) self%projfile = trim(sp_files(1))
            self%projname = get_fbody(basename(self%projfile), 'simple')
            ! update command line so that prototype copy in distr commanders transfers the info
            if( .not. cline%defined('projfile') ) call cline%set('projfile', trim(self%projname)//'.simple')
            if( .not. cline%defined('projname') ) call cline%set('projname', trim(self%projname))
        endif
        ! put a full path on projfile
        if( self%projfile .ne. '' )then
            if( file_exists(self%projfile) )then
                absname       = simple_abspath(self%projfile,errmsg='simple_parameters::new 1')
                self%projfile = trim(absname)
                self%projname = get_fbody(basename(self%projfile), 'simple')
            endif
        endif
        ! execution directory
        if( self%mkdir .eq. 'yes' )then
            if( associated(self%ptr2prg) .and. .not. str_has_substr(self%executable,'private') )then
                if( trim(self%prg) .eq. 'mkdir' )then
                    self%exec_dir = int2str(idir)//'_'//trim(self%dir)
                else if( cline%defined('dir_exec') )then
                    self%exec_dir = int2str(idir)//'_'//trim(self%dir_exec)
                else
                    self%exec_dir = int2str(idir)//'_'//trim(self%prg)
                endif
                ! make execution directory
                call simple_mkdir( filepath(PATH_HERE, trim(self%exec_dir)), errmsg="parameters:: new 2")
                if( trim(self%prg) .eq. 'mkdir' ) return
                write(logfhandle,'(a)') '>>> EXECUTION DIRECTORY: '//trim(self%exec_dir)
                ! change to execution directory directory
                call simple_chdir( filepath(PATH_HERE, trim(self%exec_dir)), errmsg="parameters:: new 3")
                if( self%sp_required )then
                    ! copy the project file
                    call simple_copy_file(trim(self%projfile), filepath(PATH_HERE, basename(self%projfile)))
                    ! update the projfile/projname
                    self%projfile = filepath(PATH_HERE, basename(self%projfile))
                    absname       = simple_abspath(self%projfile,errmsg='simple_parameters::new 4')
                    self%projfile = trim(absname)
                    self%projname = get_fbody(basename(self%projfile), 'simple')
                    ! so the projfile in cline and parameters are consistent
                    call cline%set('projfile', self%projfile)
                    ! cwd of SP-project will be updated in the builder
                endif
                ! get new cwd
                call simple_getcwd(self%cwd)
                cwd_glob = trim(self%cwd)
            endif
        endif
        ! open log file for write
        if( .not. str_has_substr(self%executable,'private') .and. STDOUT2LOG )then
            ! this is NOT a simple_private_exec execution but a simple_exec or simple_distr_exec one
            ! open the log file
            call fopen(logfhandle, status='replace', file=LOGFNAME, action='write')
        endif
        !>>> END, EXECUTION RELATED

        !>>> START, SANITY CHECKING AND PARAMETER EXTRACTION FROM ORITAB(S)/VOL(S)/STACK(S)
        ! determines whether at least one volume is on the cmdline
        vol_defined = .false.
        do i=1,size(vol_defined)
            if( cline%defined( trim('vol')//int2str(i) ))then
                vol_defined(i) = .true.
            endif
        enddo
        ! check inputted vols
        if( cline%defined('vollist') )then
            if( nlines(self%vollist)< MAXS )then
                call read_vols
            endif
        else
            if( any(vol_defined) )then
                do i=1,size(vol_defined)
                    if(vol_defined(i))call check_vol( i )
                end do
            endif
        endif
        ! check inputted mask vols
        if( cline%defined('msklist') )then
            if( nlines(self%msklist)< MAXS ) call read_masks
        endif
        ! no stack given, get ldim from volume if present
        if( self%stk .eq. '' .and. vol_defined(1) )then
            call find_ldim_nptcls(self%vols(1), self%ldim, ifoo)
            self%box  = self%ldim(1)
        endif
        ! no stack given,not vol given, get ldim from mskfile if present
        if( self%stk .eq. '' .and. .not. vol_defined(1) .and. self%mskfile .ne. '' )then
            call find_ldim_nptcls(self%mskfile, self%ldim, ifoo)
            self%box  = self%ldim(1)
        endif
        ! directories
        if( self%mkdir.eq.'yes' )then
            if( self%dir_movies(1:1).ne.PATH_SEPARATOR )self%dir_movies = PATH_PARENT//trim(self%dir_movies)
            if( self%dir_target(1:1).ne.PATH_SEPARATOR )self%dir_target = PATH_PARENT//trim(self%dir_target)
        endif
        ! project file segment
        if( cline%defined('oritype') )then
            ! this determines the spproj_iseg
            select case(trim(self%oritype))
                case('mic')
                    self%spproj_iseg = MIC_SEG
                case('stk')
                    self%spproj_iseg = STK_SEG
                case('ptcl2D')
                    self%spproj_iseg = PTCL2D_SEG
                case('cls2D')
                    self%spproj_iseg = CLS2D_SEG
                case('cls3D')
                    self%spproj_iseg = CLS3D_SEG
                case('ptcl3D')
                    self%spproj_iseg = PTCL3D_SEG
                case('out')
                    self%spproj_iseg = OUT_SEG
                case('projinfo')
                    self%spproj_iseg = PROJINFO_SEG
                case('jobproc')
                    self%spproj_iseg = JOBPROC_SEG
                case('compenv')
                    self%spproj_iseg = COMPENV_SEG
                case DEFAULT
                    write(logfhandle,*) 'oritype: ', trim(self%oritype)
                    THROW_HARD('unsupported oritype; new')
            end select
        else
            self%oritype     = 'ptcl3D'
            self%spproj_iseg = PTCL3D_SEG
        endif
        ! take care of nptcls etc.
        if( file_exists(trim(self%projfile)) )then ! existence should be the only requirement here (not sp_required)
            ! or the private_exec programs don't get what they need
            ! get nptcls/box/smpd from project file
            if( self%stream.eq.'no' )then
                if( self%spproj_iseg==OUT_SEG .or. self%spproj_iseg==CLS3D_SEG)then
                    call spproj%read_segment('out', self%projfile)
                    call spproj%get_imginfo_from_osout(smpd, box, nptcls)
                    call spproj%kill
                    if( .not.cline%defined('smpd'))   self%smpd   = smpd
                    if( .not.cline%defined('box'))    self%box    = box
                    if( .not.cline%defined('nptcls')) self%nptcls = nptcls
                else
                    call bos%open(trim(self%projfile)) ! projfile opened here
                    ! nptcls
                    if( .not. cline%defined('nptcls') ) self%nptcls = bos%get_n_records(self%spproj_iseg)
                    ! smpd/box
                    call o%new(is_ptcl=.false.)
                    select case(self%spproj_iseg)
                        case(MIC_SEG)
                            call bos%read_first_segment_record(MIC_SEG, o)
                        case DEFAULT
                            call bos%read_first_segment_record(STK_SEG, o)
                    end select
                    if( o%isthere('smpd') .and. .not. cline%defined('smpd') ) self%smpd = o%get('smpd')
                    if( o%isthere('box')  .and. .not. cline%defined('box')  ) self%box  = nint(o%get('box'))
                    call o%kill
                endif
            else
                ! nothing to do for streaming, values set at runtime
            endif
            if( .not.bos%is_opened() ) call bos%open(trim(self%projfile)) ! projfile opened here
            ! CTF plan
            select case(trim(self%oritype))
                case('ptcl2D', 'ptcl3D')
                    call bos%read_first_segment_record(STK_SEG, o)
            end select
            if( .not. cline%defined('ctf') )then
                if( o%exists() )then
                    if( o%isthere('ctf') )then
                        call o%getter('ctf', ctfflag)
                        self%ctf = trim(ctfflag)
                    endif
                endif
            endif
            if( .not. cline%defined('phaseplate') )then
                if( o%isthere('phaseplate') )then
                    call o%getter('phaseplate', phaseplate)
                    self%phaseplate   = trim(phaseplate)
                    self%l_phaseplate = trim(self%phaseplate).eq.'yes'
                else
                    self%phaseplate   = 'no'
                    self%l_phaseplate = .false.
                endif
            else
                self%l_phaseplate = trim(self%phaseplate).eq.'yes'
            endif
            call bos%close
        else if( self%stk .ne. '' )then
            if( file_exists(self%stk) )then
                if( .not. cline%defined('nptcls') )then
                    ! get number of particles from stack
                    call find_ldim_nptcls(self%stk, lfoo, self%nptcls)
                endif
            else
                write(logfhandle,'(a,1x,a)') 'Inputted stack (stk) file does not exist!', trim(self%stk)
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
                    self%box = self%ldim(1)
                endif
            else
                ! we don't check for existence of refs as they can be output as well as input (cavgassemble)
            endif
        endif
        ! check file formats
        call check_file_formats
        call double_check_file_formats
        ! make file names
        call mkfnames
        ! check box
        if( self%box > 0 .and. self%box < 26 ) THROW_HARD('box size need to be larger than 26')
        ! set refs_even and refs_odd
        if( cline%defined('refs') )then
            self%refs_even = add2fbody(self%refs, self%ext, '_even')
            self%refs_odd  = add2fbody(self%refs, self%ext, '_odd')
        endif
        ! set vols_even and vols_odd
        if( cline%defined('vol1') )then
            do istate=1,self%nstates
                self%vols_even(istate) = add2fbody(self%vols(istate), self%ext, '_even')
                self%vols_odd(istate)  = add2fbody(self%vols(istate), self%ext, '_odd' )
            end do
        endif
        !<<< END, SANITY CHECKING AND PARAMETER EXTRACTION FROM VOL(S)/STACK(S)

        !>>> START, PARALLELISATION-RELATED
        ! set split mode (even)
        self%split_mode = 'even'
        nparts_set      = .false.
        if( cline%defined('nparts') ) nparts_set  = .true.
        ! set length of number string for zero padding
        if( .not. cline%defined('numlen') )then
            if( nparts_set ) self%numlen = len(int2str(self%nparts))
        endif
        ! ldim & box from stack
        call set_ldim_box_from_stk
        ! fractional search and volume update
        if( self%update_frac <= .99) self%l_frac_update = .true.
        if( .not. cline%defined('ncunits') )then
            ! we assume that the number of computing units is equal to the number of partitions
            self%ncunits = self%nparts
        endif
        ! OpenMP threads
        if( cline%defined('nthr') )then
            nthr = nint(cline%get_rarg('nthr'))
            if( nthr == 0 )then
                !$ self%nthr = omp_get_max_threads()
                !$ call omp_set_num_threads(self%nthr)
            else
                !$ call omp_set_num_threads(self%nthr)
            endif
        else
            !$ self%nthr = omp_get_max_threads()
            !$ call omp_set_num_threads(self%nthr)
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
        if( .not. cline%defined('hp') ) self%hp = 0.5*self%dstep   ! high-pass limit
        self%fny = 2.*self%smpd                                    ! Nyqvist limit
        if( .not. cline%defined('lpstop') )then                    ! default highest resolution lp
            self%lpstop = self%fny                                 ! deafult lpstop
        endif
        if( self%fny > 0. ) self%tofny = nint(self%dstep/self%fny) ! Nyqvist Fourier index
        self%hpind_fsc = 0                                         ! high-pass Fouirer index FSC
        if( cline%defined('hp_fsc') ) self%hpind_fsc = nint(self%dstep/self%hp_fsc)
        ! set 2D low-pass limits and smpd_targets 4 scaling
        self%lplims2D(1)       = max(self%fny, self%lpstart)
        self%lplims2D(2)       = max(self%fny, self%lplims2D(1) - (self%lpstart - self%lpstop)/2.)
        self%lplims2D(3)       = max(self%fny, self%lpstop)
        self%smpd_targets2D(1) = self%lplims2D(2)*LP2SMPDFAC2D
        self%smpd_targets2D(2) = self%lplims2D(3)*LP2SMPDFAC2D
        ! check scale factor sanity
        if( self%scale < 0.00001 ) THROW_HARD('scale out if range, should be > 0; new')
        ! set default msk value
        if( cline%defined('msk') )then
            THROW_HARD('msk (mask radius in pixels) is depreciated! Use mskdiam (mask diameter in A)')
        endif
        mskdiam_default = (real(self%box) - COSMSKHALFWIDTH) * self%smpd
        msk_default     = round2even((real(self%box) - COSMSKHALFWIDTH) / 2.)
        if( cline%defined('mskdiam') )then
            self%msk = round2even((self%mskdiam / self%smpd) / 2.)

            print *, 'did set msk from mskdiam: ', self%msk

            if( self%msk > msk_default )then
                THROW_WARN('Mask diameter too large, falling back on default value')
                self%mskdiam = mskdiam_default
                self%msk     = msk_default

                print *, 'did set msk from msk default: ', self%msk

            endif
        else
            self%mskdiam = mskdiam_default
            self%msk     = msk_default
        endif

        print *, '*************************************'
        print *, 'smpd                 : ', self%smpd
        print *, 'box size             : ', self%box
        print *, 'mskdiam_default      : ', mskdiam_default
        print *, 'msk_default          : ', msk_default
        print *, 'mask diameter in A   : ', self%mskdiam
        print *, 'mask radius in pixels: ', self%msk
        print *, '*************************************'

        ! set mode of masking
        if( cline%defined('inner') )then
            THROW_HARD('inner (inner mask radius in pixels) is depreciated! Use innerdiam (inner mask diameter in A)')
        endif
        if( cline%defined('innerdiam') )then
            self%inner = round2even((self%innerdiam / self%smpd) / 2.)
            self%l_innermsk = .true.
        else
            self%l_innermsk = .false.
        endif
        ! set lpset flag
        self%l_lpset  = cline%defined('lp')
        ! set envfsc flag
        self%l_envfsc = self%envfsc .ne. 'no'
        ! set correlation weighting scheme
        self%l_corrw = self%wcrit .ne. 'no'
        if( self%l_corrw )then
            select case(trim(self%wcrit))
                case('softmax')
                    self%wcrit_enum = CORRW_CRIT
                case('zscore')
                    self%wcrit_enum = CORRW_ZSCORE_CRIT
                case('sum')
                    self%wcrit_enum = RANK_SUM_CRIT
                case('cen')
                    self%wcrit_enum = RANK_CEN_CRIT
                case('exp')
                    self%wcrit_enum = RANK_EXP_CRIT
                case('inv')
                    self%wcrit_enum = RANK_INV_CRIT
                case DEFAULT
                    THROW_HARD('unsupported correlation weighting method')
            end select
        endif
        ! set graphene flag
        self%l_graphene = self%graphene_filt .ne. 'no'
        ! checks automask related values
        if( cline%defined('mskfile') )then
            if( .not. file_exists(self%mskfile) )then
                write(logfhandle,*) 'file: ', trim(self%mskfile)
                THROW_HARD('input mask file not in cwd')
            endif
        endif
        ! focused masking
        self%l_focusmsk = .false.
        if( cline%defined('focusmsk') )then
            THROW_HARD('focusmsk (focused mask radius in pixels) is depreciated! Use focusmskdiam (focused mask diameter in A)')
        endif
        if( cline%defined('focusmskdiam') )then
            if( .not.cline%defined('mskfile') )THROW_HARD('mskfile must be provided together with focusmskdiam')
            if( .not.cline%defined('mskdiam') )then
                THROW_HARD('mskdiam must be provided together with focusmskdiam')
            endif
            self%focusmsk = round2even((self%focusmskdiam / self%smpd) / 2.)
            if(self%focusmsk >= self%msk) THROW_HARD('focusmskdiam should be smaller than mskdiam')
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
        self%kfromto(1) = max(2,int(self%dstep/self%hp)) ! high-pass Fourier index set according to hp
        self%kfromto(2) = int(self%dstep/self%lp)        ! low-pass Fourier index set according to lp
        self%kstop      = self%kfromto(2)                ! -"-
        self%lp         = max(self%fny,self%lp)          ! lowpass limit
        self%lpmed      = self%lp                        ! median lp
        if( .not. cline%defined('ydim') ) self%ydim = self%xdim
        ! set ldim
        if( cline%defined('xdim') ) self%ldim = [self%xdim,self%ydim,1]
        ! box on the command line overrides all other ldim setters
        if( cline%defined('box')  ) self%ldim = [self%box,self%box,1]
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
            if( cline%defined('ncls') )then
                if( ncls /= self%ncls )then
                    write(logfhandle,*)'ncls in ',trim(self%refs),' : ',ncls
                    write(logfhandle,*)'self%ncls : ',self%ncls
                    THROW_HARD('input number of clusters (ncls) not consistent with the number of references in stack (p%refs)')
                endif
            else
                self%ncls = ncls
            endif
        endif
        ! set remap_clusters flag
        self%l_remap_cls = .false.
        if( self%remap_cls .eq. 'yes' ) self%l_remap_cls = .true.
        ! set to particle index if not defined in cmdlin
        if( .not. cline%defined('top') ) self%top = self%nptcls
        ! set the number of input orientations
        if( .not. cline%defined('noris') ) self%noris = self%nptcls
        ! fix translation param
        self%trs = abs(self%trs)
        if( .not. cline%defined('trs') )then
            select case(trim(self%refine))
                case('snhc')
                    self%trs = 0.
                case DEFAULT
                    self%trs = MINSHIFT
            end select
        endif
        self%l_doshift = .true.
        if( self%trs < 0.5 )then
            self%trs = 0.00001
            self%l_doshift = .false.
        endif
        ! Set molecular diameter
        if( .not. cline%defined('moldiam') )then
            self%moldiam = 2. * self%msk * self%smpd
        endif
        ! check if we are dose-weighting or not
        self%l_dose_weight = .false.
        if( cline%defined('exp_time') .or. cline%defined('dose_rate') )then
            if( cline%defined('exp_time') .and. .not. cline%defined('dose_rate') )then
                THROW_HARD('need dose_rate to be part of the command line!')
            endif
            if( .not. cline%defined('exp_time') .and. cline%defined('dose_rate') )then
                THROW_HARD('need exp_time to be part of the command line!')
            endif
            self%l_dose_weight = .true.
        endif
        ! objective function used in prime2D/3D
        select case(trim(self%objfun))
            case('cc')
                self%cc_objfun = OBJFUN_CC
            case('euclid')
                self%cc_objfun = OBJFUN_EUCLID
            case DEFAULT
                write(logfhandle,*) 'objfun flag: ', trim(self%objfun)
                THROW_HARD('unsupported objective function; new')
        end select
        ! eer sanity checks
        select case(self%eer_upsampling)
            case(1,2)
            case DEFAULT
                THROW_HARD('eer_upsampling not supported: '//int2str(self%eer_upsampling))
        end select
        if( self%eer_fraction < 1 )THROW_HARD('invalid eer_fraction: '//int2str(self%eer_fraction))
        ! FILTERS
        ! matched filter and sigma needs flags
        select case(self%cc_objfun)
            case(OBJFUN_EUCLID)
                self%l_match_filt  = .false.
                self%l_needs_sigma = .true.
            case(OBJFUN_CC)
                self%l_match_filt  = (trim(self%match_filt) .eq.'yes') .and. (.not.self%l_lpset)
                self%l_needs_sigma = (trim(self%needs_sigma).eq.'yes')
                if( self%l_needs_sigma )then
                    self%l_match_filt = .false.
                endif
        end select
        ! atoms
        if( cline%defined('element') )then
            if( .not. atoms_obj%element_exists(self%element) )then
                THROW_HARD('Element: '//trim(self%element)//' unsupported for now')
            endif
        endif
        ! local resolution for filtering or  not
        self%l_locres = .false.
        if( trim(self%locres) .eq. 'yes' ) self%l_locres = .true.
        ! global dev (development) flag
        self%l_dev = .false.
        if( trim(self%dev) .eq. 'yes' ) self%l_dev = .true.
        ! sanity check imgkind
        select case(trim(self%imgkind))
            case('movie','mic','ptcl','cavg','vol','vol_cavg')
                ! alles gut!
            case DEFAULT
                write(logfhandle,*) 'imgkind: ', trim(self%imgkind)
                THROW_HARD('unsupported imgkind; new')
        end select
        ! refinement modes
        select case(trim(self%refine))
            case('snhc')
            case('shc')
            case('greedy')
            case('neigh')
            case('greedy_neigh')
            case('cont')
            case('cluster','clustersym')
            case('eval')
            case DEFAULT
                THROW_HARD('refinement mode: '//trim(self%refine)//' unsupported')
        end select
        if( str_has_substr(self%refine, 'neigh') )then
            if( .not. cline%defined('nspace')    ) self%nspace = 10000
            if( .not. cline%defined('athres')    ) self%athres = 10.
        endif
        ! motion correction
        if( self%tomo .eq. 'yes' ) self%mcpatch = 'no'
        if( self%mcpatch.eq.'yes' .and. self%nxpatch*self%nypatch<=1 ) self%mcpatch = 'no'
        if( self%mcpatch.eq.'no' )then
            self%nxpatch = 0
            self%nypatch = 0
        endif
        select case(trim(self%mcconvention))
        case('simple','unblur','motioncorr','relion','first','central')
        case DEFAULT
            THROW_HARD('Invalid entry for MCCONVENTION='//trim(self%mcconvention))
        end select
        !>>> END, IMAGE-PROCESSING-RELATED
        ! set global pointer to instance
        ! first touch policy here
        if( .not. associated(params_glob) ) params_glob => self
        if( .not. ssilent ) write(logfhandle,'(A)') '>>> DONE PROCESSING PARAMETERS'

    contains

        subroutine check_vol( i )
            integer,           intent(in) :: i
            character(len=STDLEN)         :: key
            character(len=LONGSTRLEN)     :: vol
            character(len=:), allocatable :: abs_fname
            key = 'vol'//int2str(i)
            if( cline%defined(key) )then
                vol = trim(cline%get_carg(key))
                if( vol(1:1).eq.PATH_SEPARATOR )then
                    ! already in absolute path format
                    call check_file(key, self%vols(i), notAllowed='T')
                    if( .not. file_exists(self%vols(i)) )then
                        write(logfhandle,*) 'Input volume:', trim(self%vols(i)), ' does not exist! 1'
                        stop
                    endif
                else
                    if( self%mkdir .eq. 'yes' )then
                        ! with respect to parent folder
                        ! needs to be done here because not part of the check_file list
                        vol = PATH_PARENT//trim(vol)
                        call cline%set(key, vol)
                    endif
                    call check_file(key, self%vols(i), notAllowed='T')
                    if( .not. file_exists(self%vols(i)) )then
                        write(logfhandle,*) 'Input volume:', trim(self%vols(i)), ' does not exist! 2'
                        stop
                    else
                        abs_fname = simple_abspath(self%vols(i),'parameters :: check_vol', check_exists=.false.)
                        if( len_trim( abs_fname) > LONGSTRLEN )then
                            THROW_HARD('argument too long: '//trim( abs_fname)//' new :: check_vol')
                        endif
                        self%vols(i) = trim(abs_fname)
                        call cline%set(key, trim(self%vols(i)))
                        deallocate(abs_fname)
                    endif
                endif
            endif
        end subroutine check_vol

        subroutine read_vols
            character(len=LONGSTRLEN)     :: filename, name
            character(len=:), allocatable :: abs_name
            integer                       :: nl, fnr, i, io_stat
            filename = cline%get_carg('vollist')
            if( filename(1:1).ne.PATH_SEPARATOR )then
                if( self%mkdir.eq.'yes' ) filename = PATH_PARENT//trim(filename)
            endif
            nl = nlines(filename)
            call fopen(fnr, file=filename, iostat=io_stat)
            if(io_stat /= 0) call fileiochk("parameters ; read_vols error opening "//trim(filename), io_stat)
            do i=1,nl
                read(fnr,*, iostat=io_stat) name
                if(io_stat /= 0) call fileiochk("parameters ; read_vols error reading "//trim(filename), io_stat)
                if( name .ne. '' )then
                    abs_name = simple_abspath(name,'parameters :: read_vols', check_exists=.false.)
                    self%vols(i) = trim(abs_name)
                    deallocate(abs_name)
                endif
            end do
            call fclose(fnr,errmsg="parameters ; read_vols error closing "//trim(filename))
        end subroutine read_vols

        subroutine read_masks
            character(len=LONGSTRLEN)     :: filename, name
            character(len=:), allocatable :: abs_name
            integer                       :: nl, fnr, i, io_stat
            filename = cline%get_carg('msklist')
            if( filename(1:1).ne.PATH_SEPARATOR )then
                if( self%mkdir.eq.'yes' ) filename = PATH_PARENT//trim(filename)
            endif
            nl = nlines(filename)
            call fopen(fnr, file=filename, iostat=io_stat)
            if(io_stat /= 0) call fileiochk("parameters ; read_masks error opening "//trim(filename), io_stat)
            do i=1,nl
                read(fnr,*, iostat=io_stat) name
                if(io_stat /= 0) call fileiochk("parameters ; read_masks error reading "//trim(filename), io_stat)
                if( name .ne. '' )then
                    abs_name = simple_abspath(name,errmsg='parameters :: read_masks', check_exists=.false.)
                    self%mskvols(i) = trim(abs_name)
                    deallocate(abs_name)
                endif
            end do
            call fclose(fnr,errmsg="parameters ; read_masks error closing "//trim(filename))
        end subroutine read_masks

        subroutine check_file( file, var, allowed1, allowed2, notAllowed )
            character(len=*),           intent(in)    :: file
            character(len=*),           intent(inout) :: var
            character(len=1), optional, intent(in)    :: allowed1, allowed2, notAllowed
            character(len=:), allocatable :: abspath_file
            character(len=1)              :: file_descr
            logical                       :: raise_exception
            if( cline%defined(file) )then
                var             = trim(cline%get_carg(file))
                file_descr      = fname2format(var)
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
                    write(logfhandle,*) 'This format: ', file_descr, ' is not allowed for this file: ', trim(var)
                    write(logfhandle,*) 'flag:', trim(file)
                    stop
                endif
                select case(file_descr)
                    case ('I')
                        THROW_HARD('Support for IMAGIC files is not implemented!')
                    case ('M')
                        ! MRC files are supported
                        cntfile = cntfile+1
                        checkupfile(cntfile) = 'M'
                    case ('S')
                        ! SPIDER files are supported
                        cntfile = cntfile+1
                        checkupfile(cntfile) = 'S'
                    case ('N')
                        write(logfhandle,*) 'file: ', trim(var)
                        THROW_HARD('This file format is not supported by SIMPLE')
                    case ('T','B','P','O', 'R')
                        ! text files are supported
                        ! binary files are supported
                        ! PDB files are supported
                        ! *.simple project files are supported
                        ! R=*.star format -- in testing
#ifdef USING_TIFF
                    case('J')
                        ! TIFF
                        cntfile = cntfile+1
                        checkupfile(cntfile) = 'J'
                    case('L')
                        ! .gain, which is a single tiff image
                        cntfile = cntfile+1
                        checkupfile(cntfile) = 'L'
#endif
                    case DEFAULT
                        write(logfhandle,*) 'file: ', trim(var)
                        THROW_HARD('This file format is not supported by SIMPLE')
                end select
                if( file_exists(var) )then
                    ! updates name to include absolute path
                    abspath_file = simple_abspath(var,errmsg='parameters :: check_file', check_exists=.false.)
                    if( len_trim(abspath_file) > LONGSTRLEN )then
                        THROW_HARD('argument too long: '//trim(abspath_file)//' new :: checkfile')
                    endif
                    var = trim(abspath_file)
                    call cline%set(file,trim(var))
                    deallocate(abspath_file)
                endif
            endif
        end subroutine check_file

        subroutine check_dir( dir, var )
            character(len=*), intent(in)    :: dir
            character(len=*), intent(inout) :: var
            character(len=:), allocatable   :: abspath_dir
            if( cline%defined(dir) )then
                var = trim(cline%get_carg(dir))
                if( file_exists(var) )then
                    ! updates name to include absolute path
                    abspath_dir = simple_abspath(var,errmsg='parameters :: check_dir', check_exists=.false.)
                    if( len_trim(abspath_dir) > LONGSTRLEN )then
                        THROW_HARD('argument too long: '//trim(abspath_dir)//' new :: checkdir')
                    endif
                    var = trim(abspath_dir)
                    call cline%set(dir,trim(var))
                    deallocate(abspath_dir)
                endif
            endif
        end subroutine check_dir

        subroutine check_file_formats
            integer :: i, j
            if( cntfile > 0 )then
                do i=1,cntfile
                    do j=1,cntfile
                        if( i == j ) cycle
                        if( checkupfile(i) == checkupfile(j) ) cycle ! all ok
                    end do
                end do
                call self%set_img_format(checkupfile(1))
            endif
        end subroutine check_file_formats

        subroutine double_check_file_formats
            character(len=STDLEN) :: fname
            character(len=1)      :: form
            integer :: funit, io_stat
            if( cntfile == 0 )then
                if( cline%defined('filetab') )then
                    call fopen(funit, status='old', file=self%filetab, iostat=io_stat)
                    call fileiochk("In parameters:: double_check_file_formats fopen failed "//trim(self%filetab) , io_stat)
                    read(funit,'(a256)') fname
                    form = fname2format(fname)
                    call fclose(funit, &
                        errmsg="In parameters:: double_check_file_formats fclose failed "//trim(self%filetab) )
                    call self%set_img_format(form)
                endif
            endif
        end subroutine double_check_file_formats

        subroutine mkfnames
            if( .not. cline%defined('outstk')  ) self%outstk  = 'outstk'//self%ext
            if( .not. cline%defined('outvol')  ) self%outvol  = 'outvol'//self%ext
        end subroutine mkfnames

        subroutine check_carg( carg, var )
            character(len=*), intent(in)    :: carg
            character(len=*), intent(inout) :: var
            if( cline%defined(carg) )then
                var = cline%get_carg(carg)
            endif
        end subroutine check_carg

        subroutine check_iarg( iarg, var )
            character(len=*), intent(in)  :: iarg
            integer,          intent(out) :: var
            if( cline%defined(iarg) )then
                var = nint(cline%get_rarg(iarg))
            endif
        end subroutine check_iarg

        subroutine check_rarg( rarg, var )
            character(len=*), intent(in)  :: rarg
            real, intent(out) :: var
            if( cline%defined(rarg) )then
                var = cline%get_rarg(rarg)
            endif
        end subroutine check_rarg

        subroutine set_ldim_box_from_stk
            if( cline%defined('stk') )then
                if( file_exists(self%stk) )then
                    if( cline%defined('box') )then
                    else
                        call find_ldim_nptcls(self%stk, self%ldim, ifoo)
                        self%ldim(3) = 1
                        self%box     = self%ldim(1)
                    endif
                else
                    write(logfhandle,'(a)')      'simple_parameters :: set_ldim_box_from_stk'
                    write(logfhandle,'(a,1x,a)') 'Stack file does not exist!', trim(self%stk)
                    THROW_HARD("set_ldim_box_from_stk")
                endif
            endif
        end subroutine set_ldim_box_from_stk

    end subroutine new

    subroutine set_img_format( self, ext )
        class(parameters), intent(inout) :: self
        character(len=*),  intent(in)    :: ext
        select case(trim(ext))
            case('M','D','B')
                self%ext = '.mrc'
            case('S')
                self%ext = '.spi'
#ifdef USING_TIFF
            case('J','K','L')
                ! for tiff/eer/gain we set .mrc as preferred output format
                self%ext = '.mrc'
#endif
            case DEFAULT
                write(logfhandle,*)'format: ', trim(ext)
                THROW_HARD('This file format is not supported by SIMPLE')
        end select
    end subroutine set_img_format

end module simple_parameters
