module simple_defs_stream

    integer,            parameter :: NMICS_DELTA                = 100   ! number of micrographs to increment nmics by when user requests more particles to be used in reference generation
    integer,            parameter :: DEFAULT_NTHR_MASTER        = 4     ! number of threads requested from queue system for master processes
    integer,            parameter :: PREPROC_NPARTS             = 10
    integer,            parameter :: PREPROC_NTHR               = 4
    integer,            parameter :: PREPROC_NINIPICK           = 500   ! number of micrographs to perform pick preprocessing on
    integer,            parameter :: OPENING2D_NTHR             = 32    ! number of threads requested from queue system for master process. overrides default
    integer,            parameter :: REFPICK_NPARTS             = 10
    integer,            parameter :: REFPICK_NTHR               = 4
    integer,            parameter :: SIEVING_NCLS               = 100
    integer,            parameter :: SIEVING_NPTCLS_PER_CLASS   = 0
    integer,            parameter :: SIEVING_NCHUNKS            = 2
    integer,            parameter :: SIEVING_NPARTS             = 8
    integer,            parameter :: SIEVING_NTHR               = 8
    integer,            parameter :: CLASS2D_NCLS               = 200
    integer,            parameter :: CLASS2D_NPARTS             = 10
    integer,            parameter :: CLASS2D_NTHR               = 8
    integer,            parameter :: STREAM_DEFAULT_KV          = 300
    character(len=*),   parameter :: STREAM_DEFAULT_CS          = '2.7'
    character(len=*),   parameter :: STREAM_DEFAULT_FRACA       = '0.1'
    character(len=*),   parameter :: PREPROC_JOB_NAME           = 'preprocessing'           ! name of preproc job. also used for folder name
    character(len=*),   parameter :: OPTICS_JOB_NAME            = 'optics_assignment'       ! name of optics assignment job. also used for folder name
    character(len=*),   parameter :: OPENING2D_JOB_NAME         = 'opening_2D'              ! name of opening 2D job. also used for folder name
    character(len=*),   parameter :: REFPICK_JOB_NAME           = 'reference_based_picking' ! name of reference based picking job. also used for folder name
    character(len=*),   parameter :: SIEVING_JOB_NAME           = 'particle_sieving'        ! name of particle sieving job. also used for folder name
    character(len=*),   parameter :: CLASS2D_JOB_NAME           = 'classification_2D'       ! name of 2D classification job. also used for folder name

end module simple_defs_stream