module simple_defs_fname
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! command line
integer,          parameter :: MAXNKEYS=100, KEYLEN=32
! GLOBAL STRINGS CONSTANTS
integer,          parameter :: SHORTSTRLEN          = 32   !< shorter string length
integer,          parameter :: STDLEN               = 512  !< standard string length
integer,          parameter :: LONGSTRLEN           = 1024 !< longer string length
integer,          parameter :: XLONGSTRLEN          = 8192 !< extra longer string length
! GLOBAL FILE CONSTANTS
character(len=*), parameter :: SIMPLE_SUBPROC_OUT   = 'SIMPLE_SUBPROC_OUTPUT'
character(len=*), parameter :: JOB_FINISHED_FBODY   = 'JOB_FINISHED_'
character(len=*), parameter :: TASK_FINISHED        = 'TASK_FINISHED'
! EXTENSIONS
character(len=*), parameter :: TXT_EXT              = '.txt'
character(len=*), parameter :: MAP_EXT              = '.map'
character(len=*), parameter :: BIN_EXT              = '.bin'
character(len=*), parameter :: BOX_EXT              = '.box'
character(len=*), parameter :: METADATA_EXT         = '.simple'
character(len=*), parameter :: JPG_EXT              = '.jpg'
character(len=*), parameter :: STK_EXT              = '.mrcs'
character(len=*), parameter :: STAR_EXT             = '.star'
! SUFFIXES
character(len=*), parameter :: SCALE_SUFFIX         = '_sc'
character(len=*), parameter :: THUMBNAIL_SUFFIX     = '_thumb'
character(len=*), parameter :: INTGMOV_SUFFIX       = '_intg'
character(len=*), parameter :: FORCTF_SUFFIX        = '_forctf'
character(len=*), parameter :: POWSPEC_SUFFIX       = '_pspec'
character(len=*), parameter :: LP_SUFFIX            = '_lp'
character(len=*), parameter :: PPROC_SUFFIX         = '_pproc'
character(len=*), parameter :: MIRR_SUFFIX          = '_mirr'
character(len=*), parameter :: DEN_SUFFIX           = '_den'
character(len=*), parameter :: TOPO_SUFFIX          = '_topo'
character(len=*), parameter :: BIN_SUFFIX           = '_bin'
character(len=*), parameter :: DIAMS_SUFFIX         = '_diams'
! STACK PART RELATED AND FILE FORMAT CONSTANTS
character(len=*), parameter :: STKPARTSDIR          = 'stack_parts'
character(len=*), parameter :: STKPARTFBODY         = trim(STKPARTSDIR)//'/stack_part'
character(len=*), parameter :: STKDENPARTFBODY      = trim(STKPARTSDIR)//'/stack_den_part'
character(len=*), parameter :: DEFAULT_FILE_FORMAT  = 'M'
! CLUSTER2D
character(len=*), parameter :: CLUSTER2D_ITER_FBODY = 'cluster2Ddoc_'
character(len=*), parameter :: CAVGS_ITER_FBODY     = 'cavgs_iter'
character(len=*), parameter :: MAKECAVGS_FINISHED   = 'MAKECAVGS_FINISHED'
character(len=*), parameter :: CLUSTER2D_FINISHED   = 'CLUSTER2D_FINISHED'
character(len=*), parameter :: ABINITIO2D_FINISHED  = 'ABINITIO2D_FINISHED'
character(len=*), parameter :: WFILT_SUFFIX         = '_wfilt'
character(len=*), parameter :: CLUSTER2D_ITER_THUMB = 'cls2D_thumbnail.jpeg'
character(len=*), parameter :: CLS2D_STARFBODY      = 'clusters2D'
! AUTOMASK2D
character(len=*), parameter :: BIN_OTSU             = 'binarized_otsu.mrc'
character(len=*), parameter :: BIN_OTSU_GROWN       = 'binarized_otsu_grown.mrc'
character(len=*), parameter :: BIN_OTSU_MED         = 'binarized_otsu_median.mrc'
character(len=*), parameter :: BIN_OTSU_HOLES_FILL  = 'binarized_otsu_holes_fill.mrc'
character(len=*), parameter :: MSK_OTSU             = 'masks_otsu.mrc'
character(len=*), parameter :: AMSK_OTSU            = 'automasked_otsu.mrc'
! AUTOMASK3D
character(len=*), parameter :: MSKVOL_FILE          = 'automask3D.mrc'
! REFINE3D
character(len=*), parameter :: REFINE3D_ITER_FBODY  = 'refine3Ddoc_'
character(len=*), parameter :: STARTVOL_FBODY       = 'startvol_state'
character(len=*), parameter :: VOL_FBODY            = 'recvol_state'
! 2D/3D COMMON
character(len=*), parameter :: FSC_FBODY            = 'fsc_state'
character(len=*), parameter :: FRCS_FILE            = 'frcs'//BIN_EXT
character(len=*), parameter :: ALGN_FBODY           = 'algndoc_'
character(len=*), parameter :: ARRAY_SCRIPT         = 'simple_script_array'
character(len=*), parameter :: POLARIZED_PTCLS      = 'polar_ptcls'
character(len=*), parameter :: POLARIZED_CTFS       = 'polar_ctfs'
character(len=*), parameter :: POLAR_REFS_FBODY     = 'polar_refs'
! STATS
character(len=*), parameter :: STATS_FILE           = 'simple_stats'//trim(TXT_EXT)
character(len=*), parameter :: ITERSTATS_FILE       = 'simple_iter_stats'//trim(TXT_EXT)
! SAMPLING
character(len=*), parameter :: CLASS_SAMPLING_FILE  = 'clssmp.bin'
character(len=*), parameter :: BALPROJPARTFBODY     = 'balanced_xvalid_group'
character(len=*), parameter :: RANKPROJPARTFBODY    = 'rank_group'
! PREPROCESSING
character(len=*), parameter :: PICKREFS_FBODY       = 'pickrefs'
character(len=*), parameter :: EXTRACT_STK_FBODY    = 'ptcls_from_'
character(len=*), parameter :: EXTRACT_PARAMS_FBODY = 'extract_params_'
character(len=*), parameter :: SHAPE_RANKED_CAVGS_MRCNAME = 'shaped_ranked_cavgs.mrcs'
character(len=*), parameter :: SHAPE_RANKED_CAVGS_JPGNAME = 'shaped_ranked_cavgs.jpg'
! ML
character(len=*), parameter :: SIGMA2_FBODY         = 'sigma2_noise_part'
character(len=*), parameter :: SIGMA2_GROUP_FBODY   = 'sigma2_it_'
! OLD DIRECTORIES
character(len=*), parameter :: STDERROUT_DIR        = 'stderrout/'
! NEW DIRECTORIES
character(len=*), parameter :: DIR_CTF_ESTIMATE     = 'ctf_estimate/'
character(len=*), parameter :: DIR_MOTION_CORRECT   = 'motion_correct/'
character(len=*), parameter :: DIR_INIPICK_PREPROC  = 'pick_preprocessing/'
character(len=*), parameter :: DIR_EXTRACT          = 'extract/'
character(len=*), parameter :: DIR_PICKER           = 'picker/'
character(len=*), parameter :: DIR_PREPROC          = './'
! REG CORR/ASSIGNMENT
character(len=*), parameter :: DIST_FBODY           = 'dist_part'
character(len=*), parameter :: ASSIGNMENT_FBODY     = 'assignment_part'
! STREAMING
character(len=*), parameter :: PREPROCESS_PREFIX    = 'preprocess_'
character(len=*), parameter :: STREAM_SPPROJFILES   = './stream_spprojfiles.txt'
character(len=*), parameter :: TERM_STREAM          = './SIMPLE_TERM_STREAM'
character(len=*), parameter :: PAUSE_STREAM         = './SIMPLE_PAUSE_STREAM'
character(len=*), parameter :: STREAM_REJECT_CLS    = './SIMPLE_REJECT_CLS'
character(len=*), parameter :: STREAM_SELECTED_REFS = './selected_references'
character(len=*), parameter :: STREAM_MOLDIAM       = 'moldiam.txt'
character(len=*), parameter :: STREAM_NMICS         = 'nmics.txt'
character(len=*), parameter :: STREAM_CHUNKSIZE     = 'chunksize.txt'
character(len=*), parameter :: CALCPSPEC_FINISHED   = 'CALCPSPEC_FINISHED'
character(len=*), parameter :: ABINITIO3D_FINISHED  = 'ABINITIO3D_FINISHED'
character(len=*), parameter :: DIR_SNAPSHOT         = './snapshots/'
character(len=*), parameter :: SNAPSHOT_REQUEST     = 'SNAPSHOT'
character(len=*), parameter :: GUISTATS_FILE        = '.guistats'
character(len=*), parameter :: POOLSTATS_FILE       = '.poolstats'
character(len=*), parameter :: USER_PARAMS          = 'stream_user_params.txt'
character(len=*), parameter :: SIGMAS_DIR           = './sigma2/'
character(len=*), parameter :: DIR_CHUNK            = 'chunk_'
character(len=*), parameter :: DIR_SET              = 'set_'
character(len=*), parameter :: USER_PARAMS2D        = 'stream2D_user_params.txt'
character(len=*), parameter :: OPTICS_MAP_PREFIX    = 'optics_map_'
! STARFILES
character(len=*), parameter :: MICS_STAR_BODY       = 'micrographs'
character(len=*), parameter :: PTCL2D_STAR_BODY     = 'particles2D'
! MISCELLANEOUS
character(len=3), parameter :: NIL                  = 'nil'
character(len=*), parameter :: STDERR2STDOUT        = '2>&1'
character(len=*), parameter :: CLUST_MEDIODS_FNAME  = 'clust_medoids.bin'
! character constants
character(len=*), parameter :: NEWLINE              = new_line('a')
character(len=*), parameter :: SUPPRESS_MSG         = '2>/dev/null'
character(len=*), parameter :: CSV_DELIM            = ', '
end module simple_defs_fname
