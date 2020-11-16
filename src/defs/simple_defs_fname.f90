module simple_defs_fname
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! command line
integer, parameter :: MAXNKEYS=100, KEYLEN=32
! GLOBAL STRINGS CONSTANTS
integer, parameter :: SHORTSTRLEN  = 12   !< shorter string length
integer, parameter :: STDLEN       = 256  !< standard string length
integer, parameter :: LONGSTRLEN   = 1024 !< longer string length
integer, parameter :: XLONGSTRLEN  = 4096 !< extra longer string length
! GLOBAL FILE CONSTANTS
character(len=*), parameter :: SIMPLE_SUBPROC_OUT   = 'SIMPLE_SUBPROC_OUTPUT'
character(len=*), parameter :: O_PEAKS_FBODY        = 'opeaks_part'
character(len=*), parameter :: PROJ_WEIGHTS_FBODY   = 'projection_weights_part'
! EXTENSIONS
character(len=*), parameter :: TXT_EXT              = '.txt'
character(len=*), parameter :: BIN_EXT              = '.bin'
character(len=*), parameter :: METADATA_EXT         = '.simple'
character(len=*), parameter :: JPG_EXT              = '.jpg'
! SUFFIXES
character(len=*), parameter :: SCALE_SUFFIX         = '_sc'
character(len=*), parameter :: THUMBNAIL_SUFFIX     = '_thumb'
character(len=*), parameter :: INTGMOV_SUFFIX       = '_intg'
character(len=*), parameter :: FORCTF_SUFFIX        = '_forctf'
character(len=*), parameter :: POWSPEC_SUFFIX       = '_pspec'
character(len=*), parameter :: PPROC_SUFFIX         = '_pproc'
character(len=*), parameter :: MIRR_SUFFIX          = '_mirr'
! STACK PART RELATED AND FILE FORMAT CONSTANTS
character(len=*), parameter :: STKPARTSDIR          = 'stack_parts'
character(len=*), parameter :: STKPARTFBODY         = trim(STKPARTSDIR)//'/stack_part'
character(len=*), parameter :: DEFAULT_FILE_FORMAT  = 'M'
! CLUSTER2D
character(len=*), parameter :: CLUSTER2D_ITER_FBODY = 'cluster2Ddoc_'
character(len=*), parameter :: CAVGS_ITER_FBODY     = 'cavgs_iter'
! AUTOMASK2D
character(len=*), parameter :: BIN_OTSU             = 'binarised_otsu.mrc'
character(len=*), parameter :: BIN_OTSU_GROW        = 'binarised_otsu_grown.mrc'
character(len=*), parameter :: BIN_OTSU_GROW_MED    = 'binarised_otsu_grown_median.mrc'
character(len=*), parameter :: MSK_OTSU             = 'masks_otsu.mrc'
character(len=*), parameter :: AMSK_OTSU            = 'automasked_otsu.mrc'
! REFINE3D
character(len=*), parameter :: REFINE3D_ITER_FBODY  = 'refine3Ddoc_'
character(len=*), parameter :: STARTVOL_FBODY       = 'startvol_state'
character(len=*), parameter :: VOL_FBODY            = 'recvol_state'
character(len=*), parameter :: ANISOLP_FBODY        = 'aniso_optlp_state'
character(len=*), parameter :: SNHCVOL              = 'snhc_recvol_state'
! 2D/3D COMMON
character(len=*), parameter :: FSC_FBODY            = 'fsc_state'
character(len=*), parameter :: PSSNR_FBODY          = 'pssnr_state'
character(len=*), parameter :: FRCS_FILE            = 'frcs'//BIN_EXT
character(len=*), parameter :: ALGN_FBODY           = 'algndoc_'
! STATS
character(len=*), parameter :: STATS_FILE           = 'simple_stats'//trim(TXT_EXT)
! LOCAL RESOLUTION
character(len=*), parameter :: LOCRESMAP3D_FILE     = 'locresmap3D_finds.bin'
! PREPROCESSING
character(len=*), parameter :: PICKREFS             = 'pickrefs'
character(len=*), parameter :: EXTRACT_STK_FBODY    = 'ptcls_from_'
character(len=*), parameter :: EXTRACT_PARAMS_FBODY = 'extract_params_'
! CLUSTER3D
character(len=*), parameter :: CLUSTER3D_FSC        = 'mixed_fsc'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_FRCS       = 'mixed_frcs'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_ANISOLP    = 'mixed_aniso_optlp'
character(len=*), parameter :: CLUSTER3D_VOL        = 'mixed_recvol'
! ML
character(len=*), parameter :: SIGMA2_FBODY         = 'sigma2_noise_part'
! OLD DIRECTORIES
character(len=*), parameter :: STDERROUT_DIR        = 'stderrout/'
! NEW DIRECTORIES
character(len=*), parameter :: DIR_CTF_ESTIMATE     = 'ctf_estimate/'
character(len=*), parameter :: DIR_MOTION_CORRECT   = 'motion_correct/'
character(len=*), parameter :: DIR_EXTRACT          = 'extract/'
character(len=*), parameter :: DIR_PICKER           = 'picker/'
character(len=*), parameter :: DIR_PREPROC          = './'
! STREAMING
character(len=*), parameter :: PREPROCESS_PREFIX    = 'preprocess_'
character(len=*), parameter :: STREAM_SPPROJFILES   = './stream_spprojfiles.txt'
character(len=*), parameter :: STREAM_SELMICS       = './stream_selected_mics.txt'
character(len=*), parameter :: TERM_STREAM          = './SIMPLE_TERM_STREAM'
character(len=*), parameter :: PAUSE_STREAM         = './SIMPLE_PAUSE_STREAM'
! MISCELLANEOUS
character(len=3), parameter :: NIL                  = 'nil'
character(len=*), parameter :: STDERR2STDOUT        = '2>&1'
! character constants
character(len=*), parameter :: NEWLINE              = new_line('a')
character(len=*), parameter :: SUPPRESS_MSG         = '2>/dev/null'
end module simple_defs_fname
