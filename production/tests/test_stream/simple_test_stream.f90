program simple_test_stream
! include 'simple_lib.f08'
! use simple_cmdline,              only: cmdline
! use simple_parameters,           only: parameters
! use simple_commanders_project,   only: commander_new_project
! use simple_qsys_env,             only: qsys_env
! implicit none
! integer, parameter :: NMOLDIAMS = 5
! type(parameters)            :: params
! type(cmdline)               :: cline
! type(commander_new_project) :: xnew_project
! type(cmdline)               :: cline_new_project, cline_preproc
! type(cmdline)               :: cline_estimate_moldiam, cline_gaupick
! type(qsys_env)              :: qenv
! ! Parsing
! if( command_argument_count() < 13 )then
!     write(logfhandle,'(a)') 'ERROR! See source code for usage and required parameters'
!     call exit(-1)
! endif
! call cline%parse_oldschool
! ! New project
! call cline%checkvar('projname',        1)
! ! MC & CTF
! call cline%checkvar('smpd',            2)
! call cline%checkvar('cs',              3)
! call cline%checkvar('kv',              4)
! call cline%checkvar('fraca',           5)
! call cline%checkvar('total_dose',      6)
! call cline%checkvar('dir_movies',      7)
! call cline%checkvar('gainref',         8)
! ! Compute
! call cline%checkvar('nparts',          9)
! call cline%checkvar('nthr',           10)
! ! Diameter estimation
! call cline%checkvar('moldiam',        11)
! call cline%checkvar('moldiam_max',    12)
! ! Gaussian picking
! call cline%checkvar('moldiam_refine', 13)
! call cline%check
! call params%new(cline)
! ! Initial project in currect directory, must be empty
! call cline_new_project%set('projname', params%projname)
! call cline_new_project%set('dir',      PATH_HERE)
! if( cline%defined('qsys_partition') )then
!     call cline_new_project%set('qsys_partition', cline%get_carg('qsys_partition'))
! endif
! params%projfile = params%projname//METADATA_EXT
! ! MC & CTF estimation
! call cline_preproc%set('prg',       'preproc')
! call cline_preproc%set('projfile',   params%projfile)
! call cline_preproc%set('cs',         params%cs)
! call cline_preproc%set('kv',         params%kv)
! call cline_preproc%set('fraca',      params%fraca)
! call cline_preproc%set('smpd',       params%smpd)
! call cline_preproc%set('nparts',     params%nparts)
! call cline_preproc%set('nthr',       params%nthr)
! call cline_preproc%set('total_dose', params%total_dose)
! call cline_preproc%set('gainref',    params%gainref)
! call cline_preproc%set('dir_movies', params%dir_movies)
! if( cline%defined('smpd_downscale') ) call cline_preproc%set('smpd_downscale', params%smpd_downscale)
! ! Diameter estimation
! call cline_estimate_moldiam%set('prg',         'pick_extract')
! call cline_estimate_moldiam%set('projfile',    params%projfile)
! call cline_estimate_moldiam%set('dir_target',  '1_preproc')
! call cline_estimate_moldiam%set('moldiam',     params%moldiam)
! call cline_estimate_moldiam%set('moldiam_max', params%moldiam_max)
! call cline_estimate_moldiam%set('nmoldiams',   NMOLDIAMS)
! call cline_estimate_moldiam%set('nparts',      params%nparts)
! call cline_estimate_moldiam%set('nthr',        params%nthr)
! if( .not.cline%defined('ctfresthreshold')  ) call cline_estimate_moldiam%set('ctfresthreshold',  8.0)
! if( .not.cline%defined('icefracthreshold') ) call cline_estimate_moldiam%set('icefracthreshold', 1.0)
! ! Gaussian picking
! call cline_gaupick%set('prg',        'pick_extract')
! call cline_gaupick%set('projfile',   params%projfile)
! call cline_gaupick%set('dir_target', '1_preproc')
! call cline_gaupick%set('moldiam',    params%moldiam_refine)
! call cline_gaupick%set('nparts',     params%nparts)
! call cline_gaupick%set('nthr',       params%nthr)
! if( .not.cline%defined('ctfresthreshold')  ) call cline_gaupick%set('ctfresthreshold',  8.0)
! if( .not.cline%defined('icefracthreshold') ) call cline_gaupick%set('icefracthreshold', 1.0)
! ! Executions
! call qenv%new(1)
! call xnew_project%execute_safe(cline_new_project)
! call qenv%exec_simple_prg_in_queue_async(cline_preproc,          'preproc_script',          'PREPROC_OUTPUT',    exec_bin='simple_stream')
! call qenv%exec_simple_prg_in_queue_async(cline_estimate_moldiam, 'estimate_moldiam_script', 'ESTMOLDIAM_OUTPUT', exec_bin='simple_stream')
! call sleep(900)
! call simple_touch('2_pick_extract/SIMPLE_TERM_STREAM')
! call qenv%exec_simple_prg_in_queue_async(cline_gaupick,          'gaupick_script',          'GAUPICK_OUTPUT',    exec_bin='simple_stream')
! ! End after 30mins
! call sleep(900)
! call simple_touch('1_preproc/SIMPLE_TERM_STREAM')
! call simple_touch('3_pick_extract/SIMPLE_TERM_STREAM')
end program simple_test_stream