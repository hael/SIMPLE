module simple_defs_environment
use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char
! SLURM ENVIRONMENT VARIABLES
character(len=*), parameter :: SIMPLE_DEFAULT_PARTITION_TIME   = 'SIMPLE_DEFAULT_PARTITION_TIME'
! STREAM ENVIRONMENT VARIABLES
character(len=*), parameter :: SIMPLE_STREAM_PREPROC_NTHR      = 'SIMPLE_STREAM_PREPROC_NTHR'
character(len=*), parameter :: SIMPLE_STREAM_PREPROC_PARTITION = 'SIMPLE_STREAM_PREPROC_PARTITION'
character(len=*), parameter :: SIMPLE_STREAM_PREPROC_NPARTS    = 'SIMPLE_STREAM_PREPROC_NPARTS'
character(len=*), parameter :: SIMPLE_STREAM_PICK_PARTITION    = 'SIMPLE_STREAM_PICK_PARTITION'
character(len=*), parameter :: SIMPLE_STREAM_PICK_NTHR         = 'SIMPLE_STREAM_PICK_NTHR'
character(len=*), parameter :: SIMPLE_STREAM_CHUNK_PARTITION   = 'SIMPLE_STREAM_CHUNK_PARTITION'
character(len=*), parameter :: SIMPLE_STREAM_CHUNK_NTHR        = 'SIMPLE_STREAM_CHUNK_NTHR'
character(len=*), parameter :: SIMPLE_STREAM_POOL_PARTITION    = 'SIMPLE_STREAM_POOL_PARTITION'
character(len=*), parameter :: SIMPLE_STREAM_POOL_NTHR         = 'SIMPLE_STREAM_POOL_NTHR'
character(len=*), parameter :: SIMPLE_STREAM_REFGEN_PARTITION  = 'SIMPLE_STREAM_REFGEN_PARTITION'
character(len=*), parameter :: SIMPLE_STREAM_REFGEN_NTHR       = 'SIMPLE_STREAM_REFGEN_NTHR'
end module simple_defs_environment
