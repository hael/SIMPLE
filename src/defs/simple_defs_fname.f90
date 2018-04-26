module simple_defs_fname
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! command line
integer, parameter :: MAXNKEYS=100, KEYLEN=32
! GLOBAL STRINGS CONSTANTS
integer, parameter :: STDLEN     = 256  !< standard string length
integer, parameter :: LONGSTRLEN = 2048 !< longer string length
! EXTENSIONS
character(len=*), parameter :: TXT_EXT              = '.txt'
character(len=*), parameter :: BIN_EXT              = '.bin'
character(len=*), parameter :: METADATA_EXT         = '.simple'
character(len=*), parameter :: JPG_EXT              = '.jpg'
! SUFFIXES
character(len=*), parameter :: SCALE_SUFFIX         = '_sc'
character(len=*), parameter :: THUMBNAIL_SUFFIX     = '_thumb'
character(len=*), parameter :: INTGMOV_SUFFIX       = '_intg'
character(len=*), parameter :: POWSPEC_SUFFIX       = '_pspec'
! stack part related and file format constants
character(len=*), parameter :: STKPARTSDIR          = 'stack_parts'
character(len=*), parameter :: STKPARTFBODY         = trim(STKPARTSDIR)//'/stack_part'
character(len=*), parameter :: DEFAULT_FILE_FORMAT  = 'M'
! CLUSTER2D
character(len=*), parameter :: CLUSTER2D_ITER_FBODY = 'cluster2Ddoc_'
character(len=*), parameter :: CAVGS_ITER_FBODY     = 'cavgs_iter'
! CLUSTER2D STREAM
character(len=*), parameter :: SCSTK_DIR            = 'stacks_sc/'
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
character(len=*), parameter :: FRCS_ITER_FBODY      = 'frcs_iter'
character(len=*), parameter :: ALGN_FBODY           = 'algndoc_'
! EXTRACT
character(len=*), parameter :: EXTRACT_STK_FBODY    = 'ptcls_from_'
character(len=*), parameter :: EXTRACT_PARAMS_FBODY = 'extract_params_'
! UNIDOC
character(len=*), parameter :: UNIDOC_FBODY         = 'unidoc_'
character(len=*), parameter :: UNIDOC_OUTPUT        = 'unidoc_output_'
character(len=*), parameter :: SIMPLE_UNIDOC        = 'simple_unidoc'//trim(METADATA_EXT)
! CLUSTER3D
character(len=*), parameter :: CLUSTER3D_FSC        = 'cluster3D_fsc'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_FRCS       = 'cluster3D_frcs'//BIN_EXT
character(len=*), parameter :: CLUSTER3D_ANISOLP    = 'cluster3D_aniso_optlp'
character(len=*), parameter :: CLUSTER3D_VOL        = 'cluster3D_mixed_recvol'
! OLD DIRECTORIES
character(len=*), parameter :: STDERROUT_DIR        = 'stderrout/'
! NEW DIRECTORIES
character(len=*), parameter :: DIR_CTF_ESTIMATE     = 'ctf_estimate/'
character(len=*), parameter :: DIR_MOTION_CORRECT   = 'motion_correct/'
character(len=*), parameter :: DIR_EXTRACT          = 'extract/'
character(len=*), parameter :: DIR_UNIDOC           = 'unidocs/'
character(len=*), parameter :: DIR_PICKER           = 'picker/'
character(len=*), parameter :: DIR_PREPROC          = './'
character(len=*), parameter :: DIR_PREPROC_STREAM   = './'
! STREAMING
character(len=*), parameter :: TERM_STREAM          = './SIMPLE_TERM_STREAM'
! oritype enumeration
enum, bind(c)
    enumerator :: MIC_SEG=1, STK_SEG=2, PTCL2D_SEG=3, CLS2D_SEG=4, CLS3D_SEG=5,&
    &PTCL3D_SEG=6, OUT_SEG=7, FRCS_SEG=9, FSCS_SEG=10, PROJINFO_SEG=11, JOBPROC_SEG=12,&
    &COMPENV_SEG=13, GENERIC_SEG=6
end enum
end module simple_defs_fname
