module simple_defs_fname
! EXTENSIONS
character(len=4),  parameter :: BIN_EXT            = '.bin'
character(len=4),  parameter :: METADATA_EXT       = '.txt'
! PRIME2D
character(len=11), parameter :: PRIME2D_ITER_FBODY = 'prime2Ddoc_'
character(len=10), parameter :: CAVGS_ITER_FBODY   = 'cavgs_iter'
! PRIME2D STREAM
character(len=10), parameter :: SCSTK_DIR          = 'stacks_sc/'
! PRIME3D
character(len=11), parameter :: PRIME3D_ITER_FBODY = 'prime3Ddoc_'
character(len=14), parameter :: STARTVOL_FBODY     = 'startvol_state'
character(len=12), parameter :: VOL_FBODY          = 'recvol_state'
character(len=17), parameter :: ANISOLP_FBODY      = 'aniso_optlp_state'
! PRIME COMMON
character(len=9),  parameter :: FSC_FBODY          = 'fsc_state'
character(len=10), parameter :: FRCS_FBODY         = 'frcs_state'
character(len=9),  parameter :: FRCS_ITER_FBODY    = 'frcs_iter'
character(len=8),  parameter :: ALGN_FBODY         = 'algndoc_'
! EXTRACT
character(len=11), parameter :: EXTRACT_STK_FBODY    = 'ptcls_from_'
character(len=20), parameter :: EXTRACT_PARAMS_FBODY = 'extract_params_'
! HET
character(len=11), parameter :: HET_FSC     = 'het_fsc'//BIN_EXT
character(len=12), parameter :: HET_FRCS    = 'het_frcs'//BIN_EXT
character(len=15), parameter :: HET_ANISOLP = 'het_aniso_optlp'
character(len=16), parameter :: HET_VOL     = 'het_mixed_recvol'
! SUFFIXES
character(len=3),  parameter :: SCALE_SUFFIX       = '_sc'
character(len=6),  parameter :: THUMBNAIL_SUFFIX   = '_thumb'
character(len=5),  parameter :: INTGMOV_SUFFIX     = '_intg'
character(len=6),  parameter :: POWSPEC_SUFFIX     = '_pspec'
! DIRECTORIES
character(len=10), parameter :: STDERROUT_DIR      = 'stderrout/'
character(len=9),  parameter :: PREPROC_STREAM_DIR = 'pipeline/'
character(len=21), parameter :: UNBLUR_STREAM_DIR  = PREPROC_STREAM_DIR//'micrographs/'
character(len=13), parameter :: CTF_STREAM_DIR     = PREPROC_STREAM_DIR//'ctf/'
character(len=15), parameter :: PICK_STREAM_DIR    = PREPROC_STREAM_DIR//'boxes/'
character(len=19), parameter :: EXTRACT_STREAM_DIR = PREPROC_STREAM_DIR//'particles/'
character(len=17), parameter :: UNIDOC_STREAM_DIR  = PREPROC_STREAM_DIR//'unidocs/'
! UNIDOC
character(len=14), parameter :: UNIDOC_OUTPUT      = 'unidoc_output_'
end module simple_defs_fname
