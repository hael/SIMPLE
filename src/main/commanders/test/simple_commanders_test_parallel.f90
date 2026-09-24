!@descr: the coarray test through the job queue (qsys=coarray), run by the coarray CI job
module simple_commanders_test_parallel
use simple_commanders_api
use simple_commanders_sim, only: commander_simulate_noise
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_coarrays
  contains
    procedure :: execute      => exec_test_coarrays
end type commander_test_coarrays

contains

subroutine exec_test_coarrays( self, cline )
    class(commander_test_coarrays), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
#ifdef USE_COARRAYS
    type(parameters) :: params
    type(qsys_env)   :: qenv
    type(chash)      :: job_descr
    type(cmdline)    :: cline_sim
    type(string)     :: simulated_stk
    type(commander_simulate_noise) :: xsim_noise
    integer          :: requested_nparts, requested_ncunits, requested_nptcls
    if( .not. cline%defined('nparts') ) call cline%set('nparts', 2)
    requested_nparts = cline%get_iarg('nparts')
    if( requested_nparts < 1 ) THROW_HARD('coarrays test requires nparts >= 1')
    if( .not. cline%defined('ncunits') ) call cline%set('ncunits', requested_nparts)
    requested_ncunits = cline%get_iarg('ncunits')
    if( requested_ncunits < 1 ) THROW_HARD('coarrays test requires ncunits >= 1')
    if( .not. cline%defined('nptcls') ) call cline%set('nptcls', requested_nparts)
    requested_nptcls = max(cline%get_iarg('nptcls'), requested_nparts)
    call cline%set('nptcls', requested_nptcls)
    call cline%set('qsys_name', 'coarray')
    call params%new(cline)
    simulated_stk = 'coarray_simulated_noise.mrc'
    call del_file(simulated_stk)
    call cline_sim%set('prg',      'simulate_noise')
    call cline_sim%set('mkdir',                'no')
    call cline_sim%set('box',                    32)
    call cline_sim%set('smpd',                  1.0)
    call cline_sim%set('nptcls',      params%nptcls)
    call cline_sim%set('outstk',      simulated_stk)
    call xsim_noise%execute(cline_sim)
    call cline_sim%kill
    call job_descr%new(10)
    call job_descr%set('prg',     'check_nptcls')
    call job_descr%set('stk',     simulated_stk)
    call job_descr%set('nparts',  int2str(params%nparts))
    call job_descr%set('ncunits', int2str(params%ncunits))
    call job_descr%set('nptcls',  int2str(params%nptcls))
    write(logfhandle,'(A,I0,A,I0,A)') '>>> running simulated coarray test with nparts=', &
        params%nparts, ' ncunits=', params%ncunits, ' via qsys=coarray'
    call qenv%new(params, params%nparts)
    call qenv%gen_scripts_and_schedule_jobs(job_descr, extra_params=params)
    call qenv%kill
    call job_descr%kill
    write(logfhandle,'(A)') '>>> coarray check_nptcls jobs completed on all partitions'
    call simple_end('**** SIMPLE_TEST_COARRAYS_WORKFLOW NORMAL STOP ****')
#else
    THROW_HARD('exec_test_coarrays requires a USE_COARRAYS build')
#endif
end subroutine exec_test_coarrays

end module simple_commanders_test_parallel
