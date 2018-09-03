module simple_defs_fname
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! command line
integer, parameter :: MAXNKEYS=100, KEYLEN=32
! GLOBAL STRINGS CONSTANTS
integer, parameter :: SHORTSTRLEN = 12   !< shorter string length
integer, parameter :: STDLEN      = 256  !< standard string length
integer, parameter :: LONGSTRLEN  = 1024 !< longer string length
! GLOBAL FILE CONSTANTS
character(len=*), parameter :: SIMPLE_SUBPROC_OUT   = 'SIMPLE_SUBPROC_OUTPUT'
character(len=*), parameter :: O_PEAKS_FBODY        = 'oridistributions_part'
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
! stack part related and file format constants
character(len=*), parameter :: STKPARTSDIR          = 'stack_parts'
character(len=*), parameter :: STKPARTFBODY         = trim(STKPARTSDIR)//'/stack_part'
character(len=*), parameter :: DEFAULT_FILE_FORMAT  = 'M'
! CLUSTER2D
character(len=*), parameter :: CLUSTER2D_ITER_FBODY = 'cluster2Ddoc_'
character(len=*), parameter :: CAVGS_ITER_FBODY     = 'cavgs_iter'
! REFINE3D
character(len=*), parameter :: REFINE3D_ITER_FBODY  = 'refine3Ddoc_'
character(len=*), parameter :: STARTVOL_FBODY       = 'startvol_state'
character(len=*), parameter :: VOL_FBODY            = 'recvol_state'
character(len=*), parameter :: ANISOLP_FBODY        = 'aniso_optlp_state'
character(len=*), parameter :: SNHCDOC              = 'snhc_oris'//trim(METADATA_EXT)
character(len=*), parameter :: SNHCVOL              = 'snhc_recvol_state'
! 2D/3D COMMON
character(len=*), parameter :: FSC_FBODY            = 'fsc_state'
character(len=*), parameter :: FRCS_FBODY           = 'frcs_state'
character(len=*), parameter :: FRCS_FILE            = 'frcs'//BIN_EXT
character(len=*), parameter :: ALGN_FBODY           = 'algndoc_'
! EXTRACT
character(len=*), parameter :: EXTRACT_STK_FBODY    = 'ptcls_from_'
character(len=*), parameter :: EXTRACT_PARAMS_FBODY = 'extract_params_'
! CLUSTER3D
character(len=*), parameter :: CLUSTER3D_FSC        = 'mixed_fsc'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_FRCS       = 'mixed_frcs'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_ANISOLP    = 'mixed_aniso_optlp'
character(len=*), parameter :: CLUSTER3D_VOL        = 'mixed_recvol'
! OLD DIRECTORIES
character(len=*), parameter :: STDERROUT_DIR        = 'stderrout/'
! NEW DIRECTORIES
character(len=*), parameter :: DIR_CTF_ESTIMATE     = 'ctf_estimate/'
character(len=*), parameter :: DIR_MOTION_CORRECT   = 'motion_correct/'
character(len=*), parameter :: DIR_EXTRACT          = 'extract/'
character(len=*), parameter :: DIR_PICKER           = 'picker/'
character(len=*), parameter :: DIR_PREPROC          = './'
! STREAMING
character(len=*), parameter :: STREAM_SPPROJFILES   = './stream_spprojfiles.txt'
character(len=*), parameter :: TERM_STREAM          = './SIMPLE_TERM_STREAM'
character(len=*), parameter :: PAUSE_STREAM         = './SIMPLE_PAUSE_STREAM'
! MISCELLANEOUS
character(len=*), parameter :: NIL                  = 'nil'
character(len=*), parameter :: STDERR2STDOUT        = '2>&1'
! oritype enumeration
enum, bind(c)
    enumerator :: ENUM_ORISEG=0, MIC_SEG=1, STK_SEG=2, PTCL2D_SEG=3, CLS2D_SEG=4,&
         &CLS3D_SEG=5, PTCL3D_SEG=6, OUT_SEG=7, PROJINFO_SEG=11, JOBPROC_SEG=12,&
         &COMPENV_SEG=13, GENERIC_SEG=6
end enum
! export (to STAR) type enumeration
enum, bind(c)
    enumerator :: ENUM_STARTYPE=0, MOV_STAR=1, MIC_STAR=2, STK_STAR=3, PTCL_STAR=4,&
         &CLSAVG_STAR=5, OUT_STAR=6, PROJINFO_STAR=11, GENERIC_STAR=4
end enum

end module simple_defs_fname
