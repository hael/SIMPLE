module simple_defs_fname
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! command line
integer, parameter :: MAXNKEYS=100, KEYLEN=32
! GLOBAL STRINGS CONSTANTS
integer,  parameter :: STDLEN     = 256  !< standard string length
integer,  parameter :: LONGSTRLEN = 2048 !< longer string length
! EXTENSIONS
character(len=4),      parameter :: BIN_EXT                   = '.bin'
character(len=7),      parameter :: METADATA_EXT              = '.txt'
! SUFFIXES
character(len=3),      parameter :: SCALE_SUFFIX              = '_sc'
character(len=6),      parameter :: THUMBNAIL_SUFFIX          = '_thumb'
character(len=5),      parameter :: INTGMOV_SUFFIX            = '_intg'
character(len=6),      parameter :: POWSPEC_SUFFIX            = '_pspec'
! stack part related and file format constants
character(len=KEYLEN), parameter :: STKPARTSDIR               = 'stack_parts'
character(len=KEYLEN), parameter :: STKPARTFBODY              = trim(STKPARTSDIR)//'/stack_part'
character(len=1),      parameter :: DEFAULT_FILE_FORMAT       = 'M'
! PRIME2D
character(len=KEYLEN), parameter :: CLUSTER2D_ITER_FBODY      = 'cluster2Ddoc_'
character(len=KEYLEN), parameter :: CAVGS_ITER_FBODY          = 'cavgs_iter'
! PRIME2D STREAM
character(len=KEYLEN), parameter :: SCSTK_DIR                 = 'stacks_sc/'
! PRIME3D
character(len=KEYLEN), parameter :: REFINE3D_ITER_FBODY       = 'refine3Ddoc_'
character(len=KEYLEN), parameter :: STARTVOL_FBODY            = 'startvol_state'
character(len=KEYLEN), parameter :: VOL_FBODY                 = 'recvol_state'
character(len=KEYLEN), parameter :: ANISOLP_FBODY             = 'aniso_optlp_state'
character(len=KEYLEN), parameter :: SNHCDOC                   = 'snhc_oris'//trim(METADATA_EXT)
character(len=KEYLEN), parameter :: SNHCVOL                   = 'snhc_recvol_state'
! PRIME3D COMMON
character(len=KEYLEN), parameter :: FSC_FBODY                 = 'fsc_state'
character(len=KEYLEN), parameter :: FRCS_FBODY                = 'frcs_state'
character(len=KEYLEN), parameter :: FRCS_ITER_FBODY           = 'frcs_iter'
character(len=KEYLEN), parameter :: ALGN_FBODY                = 'algndoc_'
! EXTRACT
character(len=KEYLEN), parameter :: EXTRACT_STK_FBODY         = 'ptcls_from_'
character(len=KEYLEN), parameter :: EXTRACT_PARAMS_FBODY      = 'extract_params_'
! UNIDOC
character(len=KEYLEN), parameter :: UNIDOC_OUTPUT             = 'unidoc_output_'
! CLUSTER3D
character(len=KEYLEN), parameter :: CLUSTER3D_FSC             = 'cluster3D_fsc'//BIN_EXT
character(len=KEYLEN), parameter :: CLUSTER3D_FRCS            = 'cluster3D_frcs'//BIN_EXT
character(len=KEYLEN), parameter :: CLUSTER3D_ANISOLP         = 'cluster3D_aniso_optlp'
character(len=KEYLEN), parameter :: CLUSTER3D_VOL             = 'cluster3D_mixed_recvol'
! DIRECTORIES
character(len=KEYLEN), parameter :: STDERROUT_DIR             = 'stderrout/'
character(len=KEYLEN), parameter :: PREPROCESS_STREAM_DIR     = 'pipeline/'
character(len=KEYLEN), parameter :: MOTION_CORRECT_STREAM_DIR = trim(PREPROCESS_STREAM_DIR)//'micrographs/'
character(len=KEYLEN), parameter :: CTF_STREAM_DIR            = trim(PREPROCESS_STREAM_DIR)//'ctf/'
character(len=KEYLEN), parameter :: PICK_STREAM_DIR           = trim(PREPROCESS_STREAM_DIR)//'boxes/'
character(len=KEYLEN), parameter :: EXTRACT_STREAM_DIR        = trim(PREPROCESS_STREAM_DIR)//'particles/'
character(len=KEYLEN), parameter :: UNIDOC_STREAM_DIR         = trim(PREPROCESS_STREAM_DIR)//'unidocs/'
! oritype enumeration
enum, bind(c)
    enumerator :: STK_SEG = 1, PTCL2D_SEG = 2, CLS2D_SEG = 3,&
    &CLS3D_SEG = 4, PTCL3D_SEG = 5, PROJINFO_SEG=11, JOBPROC_SEG = 12
end enum
end module simple_defs_fname
