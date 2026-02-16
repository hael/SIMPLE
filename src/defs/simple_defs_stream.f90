!@descr: constants used in streaming
module simple_defs_stream
character(len=*), parameter :: CHUNK_CLS_REJECTED         = 'cls_rejected_chunks.mrc'
character(len=*), parameter :: CHUNK_PROJNAME             = 'chunk'
character(len=*), parameter :: CLASS2D_JOB_NAME           = 'classification_2D'       ! name of 2D classification job. also used for folder name
character(len=*), parameter :: DIR_STREAM                 = './spprojs/'              ! location for projects to be processed
character(len=*), parameter :: DIR_STREAM_COMPLETED       = './spprojs_completed/'    ! location for projects processed
character(len=*), parameter :: MICSPPROJ_FNAME            = './streamdata.simple'
character(len=*), parameter :: OPENING2D_JOB_NAME         = 'opening_2D'              ! name of opening 2D job. also used for folder name
character(len=*), parameter :: OPTICS_JOB_NAME            = 'optics_assignment'       ! name of optics assignment job. also used for folder name
character(len=*), parameter :: POOL_DIR                   = ''                        ! should be './pool/' for tidyness but difficult with gui
character(len=*), parameter :: POOL_DISTR_EXEC_FNAME      = './distr_cluster2D_pool'
character(len=*), parameter :: POOL_LOGFILE               = 'simple_log_cluster2D_pool'
character(len=*), parameter :: POOL_PROJFILE              = 'cluster2D.simple'
character(len=*), parameter :: PREPROC_JOB_NAME           = 'preprocessing'           ! name of preproc job. also used for folder name
character(len=*), parameter :: REFPICK_JOB_NAME           = 'reference_based_picking' ! name of reference based picking job. also used for folder name
character(len=*), parameter :: REJECTED_CLS_STACK         = './rejected_cls.mrc'
character(len=*), parameter :: SIEVING_JOB_NAME           = 'particle_sieving'        ! name of particle sieving job. also used for folder name
character(len=*), parameter :: SIEVING_REFS_FNAME         = 'sieving_references'
character(len=*), parameter :: STREAM_DEFAULT_CS          = '2.7'
character(len=*), parameter :: STREAM_DEFAULT_FRACA       = '0.1'
integer(kind=8),  parameter :: FLUSH_TIMELIMIT            = 900                       ! time (secs) after which leftover particles join the pool IF the 2D analysis is paused
integer,          parameter :: CHUNK_CC_ITERS             = 8                         ! maximum number of correlatiion-based iterations for chunks
integer,          parameter :: CHUNK_EXTR_ITER            = 3                         ! starting extremal iteration for chunks
integer,          parameter :: CHUNK_MINITS               = 13                        ! minimum number of iterations for chunks
integer,          parameter :: CHUNK_MAXITS               = CHUNK_MINITS + 2          ! maximum number of iterations for chunks
integer,          parameter :: CHUNK_MINBOXSZ             = 128                       ! minimum boxsize for scaling
integer,          parameter :: CLASS2D_NCLS               = 200
integer,          parameter :: CLASS2D_NPARTS             = 10
integer,          parameter :: CLASS2D_NTHR               = 8
integer,          parameter :: DEFAULT_NTHR_MASTER        = 4                         ! number of threads requested from queue system for master processes
integer,          parameter :: INACTIVE_TIME              = 900                       ! inactive time trigger for writing project file
integer,          parameter :: LONGTIME                   = 60                        ! time lag after which a movie/project is processed
integer,          parameter :: NMICS_DELTA                = 100                       ! number of micrographs to increment nmics by when user requests more particles to be used in reference generation
integer,          parameter :: OPENING2D_NTHR             = 32                        ! number of threads requested from queue system for master process. overrides default
integer,          parameter :: PAUSE_NITERS               = 5                         ! # of iterations after which 2D analysis is paused
integer,          parameter :: PAUSE_TIMELIMIT            = 600                       ! time (secs) after which 2D analysis is paused
integer,          parameter :: POOL_FREQ_REJECTION        = 5                         ! pool class rejection performed every POOL_FREQ_REJECTION iteration
integer,          parameter :: POOL_NPREV_RES             = 5                         ! # of previous resolution resolutions to store for resolution update (>=2)
integer,          parameter :: PREPROC_NINIPICK           = 500                       ! number of micrographs to perform pick preprocessing on
integer,          parameter :: PREPROC_NPARTS             = 10
integer,          parameter :: PREPROC_NTHR               = 4
integer,          parameter :: REFPICK_NPARTS             = 10
integer,          parameter :: REFPICK_NTHR               = 4
integer,          parameter :: SHORTWAIT                  = 2                         ! movie folder watched every SHORTTIME seconds in shmem
integer,          parameter :: SIEVING_MATCH_CAVGS_MAX    = 20
integer,          parameter :: SIEVING_REF_CAVGS_MAX      = 100
integer,          parameter :: SIEVING_NCHUNKS            = 2
integer,          parameter :: SIEVING_NCLS               = 100
integer,          parameter :: SIEVING_NPARTS             = 8
integer,          parameter :: SIEVING_NPTCLS_PER_CLASS   = 0
integer,          parameter :: SIEVING_NTHR               = 8
integer,          parameter :: STREAM_DEFAULT_KV          = 300
integer,          parameter :: STREAM_NMOVS_SET           = 5                         !< number of movies processed at once (>1)
integer,          parameter :: STREAM_NMOVS_SET_TIFF      = 3                         !< number of TIFF movies processed at once (>1)
integer,          parameter :: STREAM_NPTCLS_MAX          = 500000                    !< 2D analysis: cap for adjusting update_frac in 2D streaming
integer,          parameter :: STREAM_SRCHLIM             = 5                         !< 2D analysis: maximum # of systematic iterations for streaming 2D pool
integer,          parameter :: WAITTIME                   = 10                        ! movie folder watched every WAITTIME seconds
real,             parameter :: ABS_VAR_THRESHOLD          = 1.5                       !< class rejection: image variance threshold (absolute value)
real,             parameter :: CLS_REJECT_STD             = 2.5                       !< class rejection: # deviations for 2D class selection/rejection
real,             parameter :: FRAC_SKIP_REJECTION        = 0.7                       !< 2D analysis: When the number of classes to reject is too high rejection is skipped
real,             parameter :: LOWRES_REJECT_THRESHOLD    = 199.                      !< class rejection: Deactivates resolution-based rejection when lpthres > LOWRES_REJECT_THRESHOLD
real,             parameter :: MEAN_THRESHOLD             = -8.0                      !< class rejection: image mean     threshold (Mahalabonis distance)
real,             parameter :: MINMAX_THRESHOLD           = 2.0                       !< class rejection: image min & max threshold (absolute value)
real,             parameter :: POOL_SMPD_HARD_LIMIT       = 1.5                       ! Pixel size hard limit -> max resolution=3Angs
real,             parameter :: REL_VAR_THRESHOLD          = 6.0                       !< class rejection: image variance threshold (Mahalabonis distance)
real,             parameter :: STREAM_ASTIG_THRESHOLD     = 10.0                      !< preprocessing: Stream astigmatism rejection threshold
real,             parameter :: STREAM_CTFRES_THRESHOLD    = 10.0                      !< preprocessing: Stream ctfres rejection threshold (Angstroms)
real,             parameter :: STREAM_ICEFRAC_THRESHOLD   = 1.0                       !< preprocessing: Stream icefrac rejection threshold
real,             parameter :: STREAM_RES_THRESHOLD       = 35.0                      !< class rejection: Default streaming best resolution rejection threshold
real,             parameter :: TVD_THRESHOLD              = 0.55                      !< class rejection: Total Variation Distance of image distributions

real, parameter, dimension(21)  :: ASTIG_BINS    = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(19)  :: CTFRES_BINS   = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(21)  :: ICESCORE_BINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

type scaled_dims
    real    :: smpd=0., msk=0.
    integer :: box=0, boxpd=0
end type scaled_dims

end module simple_defs_stream