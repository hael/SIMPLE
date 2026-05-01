!@descr: Queue-system execution environment with optional persistent worker dispatch
module simple_qsys_env
use simple_core_module_api
use simple_defs_environment
use simple_cmdline,           only: cmdline
use simple_qsys_funs,         only: qsys_watcher, qsys_cleanup
use simple_qsys_factory,      only: qsys_factory
use simple_qsys_base,         only: qsys_base
use simple_qsys_local,        only: qsys_local
use simple_qsys_slurm,        only: qsys_slurm
use simple_qsys_persistent_worker, only: qsys_persistent_worker
use simple_qsys_lsf,          only: qsys_lsf
use simple_qsys_pbs,          only: qsys_pbs
use simple_qsys_sge,          only: qsys_sge
use simple_qsys_ctrl,         only: qsys_ctrl
use simple_parameters,        only: parameters
use simple_persistent_worker_server, only: persistent_worker
implicit none

public :: qsys_env
private
#include "simple_local_flags.inc"

type :: qsys_env
    ! --- public state (read directly by callers) ---
    integer, allocatable,        public  :: parts(:,:)             !< particle ranges [ipart, start/end]
    type(qsys_ctrl),             public  :: qscripts               !< script controller for the dispatch path
    type(chash),                 public  :: qdescr                 !< key/value queue-job description
    ! --- private implementation details ---
    type(qsys_ctrl),             private :: base_qscripts          !< script controller for the base/worker-launch path
    type(string),                private :: simple_exec_bin        !< fully resolved path to submission executable
    type(qsys_factory),          private :: qsys_fac_base          !< factory that owns base_qsys
    type(qsys_factory),          private :: qsys_fac_dispatch      !< factory that owns dispatch_qsys
    class(qsys_base),   pointer, private :: base_qsys     => null() !< base backend (local/slurm/lsf/…)
    class(qsys_base),   pointer, private :: dispatch_qsys => null() !< dispatch backend (persistent_worker)
    integer,                     private :: nparts                 !< number of job partitions
    logical,                     private :: existence   = .false.  !< .true. after new(); .false. after kill()
    contains
    procedure :: new
    procedure :: exists
    procedure :: gen_script
    procedure :: gen_scripts_and_schedule_jobs
    procedure :: gen_subproject_scripts_and_schedule
    procedure :: schedule_subproject_jobs
    procedure :: exec_simple_prg_in_queue
    procedure :: exec_simple_prg_in_queue_async
    procedure :: exec_simple_prgs_in_queue_async
    procedure :: start_persistent_workers
    procedure :: get_exec_bin
    procedure :: get_qsys
    procedure :: get_navail_computing_units
    procedure :: kill
end type qsys_env

contains

    !> Build queue description from runtime parameters and environment variables.
    !! Falls back to shell environment variables (SIMPLE_PATH, SIMPLE_QSYS, …)
    !! when the corresponding parameter fields are empty.
    subroutine init_qdescr_from_runtime( qdescr, params, qsys_name, qsys_partition )
        type(chash),               intent(inout) :: qdescr           !< populated with scheduler metadata on return
        class(parameters),         intent(in)    :: params           !< run-time parameter set
        class(string), optional,   intent(in)    :: qsys_name        !< explicit scheduler name; overrides params and environment
        class(string), optional,   intent(in)    :: qsys_partition   !< explicit partition; overrides environment
        type(string)          :: env_var, qsys_name_here, job_name
        integer               :: iostat
        ! Resolve SIMPLE installation path – mandatory.
        call qdescr%set('simple_path', simple_getenv('SIMPLE_PATH', iostat))
        if( iostat /= 0 ) THROW_HARD('SIMPLE_PATH is not defined in your environment; qsys_env::init_qdescr_from_runtime')
        ! Resolve scheduler name: argument > params%qsys_name > SIMPLE_QSYS env var.
        if( present(qsys_name) ) then
            qsys_name_here = qsys_name
        else
            qsys_name_here = string(trim(params%qsys_name))
            if( qsys_name_here%strlen_trim() == 0 ) then
                env_var = simple_getenv('SIMPLE_QSYS', iostat)
                if( iostat /= 0 ) THROW_HARD('SIMPLE_QSYS is not defined in your environment; qsys_env::init_qdescr_from_runtime')
                qsys_name_here = env_var
            end if
        end if
        call qdescr%set('qsys_name', qsys_name_here)
        ! User e-mail falls back to a placeholder when not configured.
        env_var = simple_getenv('SIMPLE_EMAIL', iostat)
        if( iostat /= 0 ) env_var = 'my.name@uni.edu'
        call qdescr%set('user_email', env_var)
        call qdescr%set('time_per_image', string(TIME_PER_IMAGE_DEFAULT))
        call qdescr%set('walltime', string(params%walltime))
        ! Partition: argument > SIMPLE_QSYS_PARTITION env var; omitted when neither is set.
        if( present(qsys_partition) ) then
            call qdescr%set('qsys_partition', qsys_partition)
        else
            env_var = simple_getenv('SIMPLE_QSYS_PARTITION', iostat)
            if( iostat == 0 ) call qdescr%set('qsys_partition', env_var)
        end if
        job_name = 'simple_'//params%prg%to_char()
        call qdescr%set('job_name', job_name)
        call qdescr%set('job_memory_per_task', string(JOB_MEMORY_PER_TASK_DEFAULT))
        call job_name%kill
        call qsys_name_here%kill
        call env_var%kill
    end subroutine init_qdescr_from_runtime

    !> Initialize qsys environment state, script controllers, and optional worker backend.
    !! A qsys_name ending in '_worker' (e.g. 'local_worker', 'slurm_worker') activates the
    !! TCP persistent-worker dispatch path; all other names use the standard backend directly.
    subroutine new( self, params, nparts, stream, numlen, nptcls, exec_bin, qsys_name, qsys_nthr, qsys_partition )
        use simple_sp_project, only: sp_project
        class(qsys_env),         intent(inout) :: self
        class(parameters),       intent(in)    :: params             !< run-time parameter set
        integer,                 intent(in)    :: nparts             !< number of parallel job partitions
        logical,       optional, intent(in)    :: stream             !< .true. for streaming / continuous mode
        integer,       optional, intent(in)    :: numlen             !< zero-padding width for script numbering
        integer,       optional, intent(in)    :: nptcls             !< override params%nptcls particle count
        integer,       optional, intent(in)    :: qsys_nthr          !< override per-job thread count
        class(string), optional, intent(in)    :: exec_bin           !< override default submission executable
        class(string), optional, intent(in)    :: qsys_name          !< override scheduler name from params/env
        class(string), optional, intent(in)    :: qsys_partition     !< override scheduler partition from params/env
        type(ori)             :: compenv_o
        type(sp_project)      :: spproj
        type(string)          :: qsnam, tpi, hrs_str, mins_str, secs_str
        integer               :: nthr_workers, n_workers
        character(len=STDLEN) :: default_time_env
        integer               :: partsz, hrs, mins, secs, nptcls_here, envlen
        real                  :: rtpi, tot_time_sec
        logical               :: sstream
        integer, parameter    :: MAXENVKEYS = 30
        call self%kill
        sstream = .false.
        if( present(stream) ) sstream = stream
        self%nparts = nparts
        nptcls_here = params%nptcls
        if( present(nptcls) ) nptcls_here = nptcls
        select case(trim(params%split_mode))
            case('even')
                self%parts = split_nobjs_even(nptcls_here, self%nparts)
                partsz     = self%parts(1,2) - self%parts(1,1) + 1
            case('singles')
                allocate(self%parts(nptcls_here,2))
                self%parts(:,:) = 1
                partsz          = 1
            case('stream')
                self%nparts = params%ncunits
                allocate(self%parts(nptcls_here,2)) ! unused
                self%parts(:,:) = 1                 ! unused
                partsz          = 1                 ! unused
            case DEFAULT
                THROW_HARD('split_mode: '//trim(params%split_mode)//' is unsupported; new')
        end select
        ! Retrieve queue metadata from project compenv when available, otherwise
        ! synthesize it directly from runtime parameters and shell environment.
        call self%qdescr%new(MAXENVKEYS)
        if( params%projfile /= '' ) then
            call spproj%read_segment('compenv', params%projfile)
            call spproj%compenv%get_ori(1, compenv_o)
            self%qdescr = compenv_o%ori2chash()
        else
            call init_qdescr_from_runtime(self%qdescr, params, qsys_name, qsys_partition)
        end if
        ! Derive per-job runtime budget.
        call get_environment_variable(SIMPLE_DEFAULT_PARTITION_TIME, default_time_env, envlen)
        if( envlen > 0 .and. trim(default_time_env) .eq. 'true' ) then
            if( self%qdescr%isthere('job_time') ) then
                call self%qdescr%delete('job_time')
            end if
        else if( self%qdescr%isthere('time_per_image') ) then
            tpi          = self%qdescr%get('time_per_image')
            rtpi         = str2real(tpi)
            tot_time_sec = rtpi*real(partsz)
            if( self%qdescr%isthere('walltime') ) then
                tot_time_sec = min(tot_time_sec, str2real(self%qdescr%get('walltime')))
            else
                tot_time_sec = min(tot_time_sec, real(WALLTIME_DEFAULT))
            end if
            tot_time_sec = min(tot_time_sec, real(params%walltime)) ! command line override
            hrs          = int(tot_time_sec/3600.)
            hrs_str      = int2str(hrs)
            mins         = int((tot_time_sec - 3600.*real(hrs))/60.)
            mins_str     = int2str(mins)
            secs         = int(tot_time_sec - 3600.*real(hrs) - 60.*real(mins))
            secs_str     = int2str(secs)
            call self%qdescr%set('job_time', '0-'//hrs_str%to_char()//':'//mins_str%to_char()//':'//secs_str%to_char())
        end if

        if( present(exec_bin) ) then
            self%simple_exec_bin = filepath(self%qdescr%get('simple_path'),string('bin'),exec_bin)
        else
            self%simple_exec_bin = filepath(self%qdescr%get('simple_path'),string('bin'),string('simple_private_exec'))
        end if
        if( present(qsys_nthr) ) then
            call self%qdescr%set('job_cpus_per_task', int2str(qsys_nthr))                  ! overrides env file and params
        else
            call self%qdescr%set('job_cpus_per_task', int2str(params%nthr))                ! overrides env file
        end if
        if( present(qsys_partition) ) call self%qdescr%set('qsys_partition', qsys_partition) ! overrides env file and params
        call self%qdescr%set('job_nparts', int2str(params%nparts))                         ! overrides env file
        if( present(qsys_name) ) call self%qdescr%set('qsys_name', qsys_name)
        qsnam = self%qdescr%get('qsys_name')
        ! The '_worker' suffix on a qsys name (e.g. 'local_worker', 'slurm_worker')
        ! opts into the TCP persistent-worker backend.  Strip the suffix so the
        ! underlying qsys factory still resolves the base scheduler (local, slurm, …).
        if( qsnam%has_substr(string('_worker')) ) then
            qsnam = qsnam%substr_remove(string('_worker'))
            ! Resolve worker concurrency: explicit params override ncunits; qsys_nthr overrides params%nthr.
            nthr_workers = params%nthr
            n_workers    = max(1, params%ncunits)
            if( params%workers > 0     ) n_workers    = params%workers
            if( params%worker_nthr > 0 ) nthr_workers = params%worker_nthr
            if( present(qsys_nthr)     ) nthr_workers = qsys_nthr
            call self%qdescr%set('job_cpus_per_task', int2str(nthr_workers)) ! workers use the main thread count for their own per-task CPU allocation; nthr_workers is only for the worker server
            ! Build two separate backends: base backend launches workers via the underlying scheduler;
            ! dispatch backend routes tasks to the running TCP worker pool.
            call self%qsys_fac_base%new(qsnam, self%base_qsys)
            call self%qsys_fac_dispatch%new(string('persistent_worker'), self%dispatch_qsys)
            ! qscripts dispatches through the TCP worker path.
            call self%qscripts%new(self%simple_exec_bin, self%dispatch_qsys, self%parts,&
                &[1, self%nparts], params%ncunits, sstream, numlen, nthr_worker=nthr_workers)
            ! base_qscripts submits directly via the base scheduler (used to launch worker processes).
            call self%base_qscripts%new(self%simple_exec_bin, self%base_qsys, self%parts,&
                &[1, self%nparts], params%ncunits, sstream, numlen)
            ! Reuse a live server when configuration is compatible; otherwise start a fresh one.
            if( associated(persistent_worker%server) ) then
                if( .not. persistent_worker%server%is_running()      ) THROW_HARD('cannot reuse existing worker server that is not running;')
                if( persistent_worker%launch_backend /= qsnam        ) THROW_HARD('cannot reuse existing worker server with different backend; kill the server or use a persistent qsys name')
                if( nthr_workers > persistent_worker%nthr_per_worker ) THROW_HARD('cannot reuse existing worker server with lower nthr_per_worker than requested;')
                if( n_workers > persistent_worker%n_workers          ) THROW_HARD('cannot reuse existing worker server with lower n_workers than requested;')
            else
                persistent_worker%launch_backend  = qsnam
                persistent_worker%nthr_per_worker = nthr_workers
                persistent_worker%n_workers       = n_workers
                allocate(persistent_worker%server)
                call persistent_worker%server%new(persistent_worker%nthr_per_worker)
                call self%start_persistent_workers()
            end if
        else
            ! Standard path: a single backend handles both script generation and dispatch.
            call self%qsys_fac_base%new(qsnam, self%base_qsys)
            call self%qsys_fac_dispatch%new(qsnam, self%dispatch_qsys)
            call self%qscripts%new(self%simple_exec_bin, self%dispatch_qsys, self%parts,&
            &[1, self%nparts], params%ncunits, sstream, numlen)
        end if
        ! Release temporary strings and project objects.
        call qsnam%kill
        if( params%projfile /= '' ) call compenv_o%kill
        call spproj%kill
        self%existence = .true.
    end subroutine new

    !> Return .true. when the environment has been initialized.
    function exists( self ) result( is )
        class(qsys_env) :: self
        logical         :: is
        is = self%existence
    end function exists

    !> Generate a single submission script.
    subroutine gen_script( self, cline, script_name, prg_output )
        use simple_cmdline, only: cmdline
        class(qsys_env)              :: self
        class(cmdline)               :: cline
        class(string),    intent(in) :: script_name, prg_output
        call self%qscripts%generate_script(cline, self%qdescr, script_name, prg_output)
    end subroutine gen_script

    !> Generate scripts and schedule standard or array jobs.
    subroutine gen_scripts_and_schedule_jobs( self, job_descr, part_params, algnfbody, array, extra_params )
        class(qsys_env)                        :: self
        class(chash)                           :: job_descr
        class(chash),     optional             :: part_params(self%nparts)
        type(parameters), optional, intent(in) :: extra_params
        class(string),    optional             :: algnfbody
        logical,          optional             :: array
        logical :: aarray
        aarray = .false.
        if( present(array) ) aarray = array
        ! we only support array execution by SLURM and distributed local shared-memory
        select type(pbase_qsys => self%base_qsys)
            class is(qsys_local)
                aarray = .false.
            class is(qsys_slurm)
                ! keep aarray value
            class is(qsys_lsf)
                aarray = .false.
            class is(qsys_sge)
                aarray = .false.
            class is(qsys_pbs)
                aarray = .false.
            class is(qsys_persistent_worker)
                aarray = .false.
        end select
        if( present(extra_params) ) then
            call qsys_cleanup(extra_params)
        end if
        if( aarray ) then
            call self%qscripts%generate_array_script(job_descr, string(MRC_EXT), self%qdescr, outfile_body=algnfbody, part_params=part_params)
            call self%qscripts%schedule_array_jobs
        else
            call self%qscripts%generate_scripts(job_descr, string(MRC_EXT), self%qdescr, outfile_body=algnfbody, part_params=part_params, extra_params=extra_params)
            call self%qscripts%schedule_jobs
        end if
    end subroutine gen_scripts_and_schedule_jobs

    !>  \brief  Generate scripts for subprojects and schedule them for parallel execution.
    !  Each element of jobs_descr(:) describes one independent subproject.
    !  The qsys_ctrl machinery generates one script per subproject and
    !  schedule_jobs() dispatches them in parallel respecting ncunits concurrency.
    !  Fully qsys-agnostic: works with local, SLURM, LSF, PBS, SGE.
    subroutine gen_subproject_scripts_and_schedule( self, jobs_descr, exec_bin, subproj_dirs )
        class(qsys_env)                     :: self
        type(chash),          intent(inout) :: jobs_descr(:)   !< one job description per subproject
        class(string), optional, intent(in) :: exec_bin        !< override executable binary
        class(string), optional, intent(in) :: subproj_dirs(:) !< per-subproject working directories
        call self%qscripts%generate_scripts_subprojects(jobs_descr, self%qdescr, exec_bin, subproj_dirs)
        call self%qscripts%schedule_subproject_jobs
    end subroutine gen_subproject_scripts_and_schedule

    !> Schedule already-generated subproject jobs.
    subroutine schedule_subproject_jobs( self )
        class(qsys_env), intent(inout) :: self
        call self%qscripts%schedule_subproject_jobs
    end subroutine schedule_subproject_jobs

    !> Generate and submit a single script, then block until its finish indicator file appears.
    subroutine exec_simple_prg_in_queue( self, cline, finish_indicator )
        use simple_cmdline, only: cmdline
        class(qsys_env), intent(inout) :: self
        class(cmdline),  intent(inout) :: cline
        class(string),   intent(in)    :: finish_indicator !< path of the completion sentinel file
        character(len=*), parameter    :: SCRIPT_NAME = 'simple_script_single'
        type(chash)                    :: job_descr
        call del_file(finish_indicator)
        call cline%gen_job_descr(job_descr)
        call self%qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, string(SCRIPT_NAME))
        call wait_for_closure(string(SCRIPT_NAME))
        call self%qscripts%submit_script(string(SCRIPT_NAME))
        call qsys_watcher(finish_indicator)
        call del_file(finish_indicator)
        call job_descr%kill
    end subroutine exec_simple_prg_in_queue

    !> Generate and submit one script without blocking; caller is responsible for monitoring.
    !! When base=.true. the script is routed through base_qscripts (the underlying
    !! scheduler) rather than through the dispatch (TCP worker) path.
    subroutine exec_simple_prg_in_queue_async( self, cline, script_name, outfile, exec_bin, base )
        use simple_cmdline, only: cmdline
        class(qsys_env),         intent(inout) :: self
        class(cmdline),          intent(in)    :: cline        !< command-line parameters for the job
        class(string),           intent(in)    :: script_name  !< path of the script to generate and submit
        class(string),           intent(in)    :: outfile      !< stdout/stderr log path for the job
        class(string), optional, intent(in)    :: exec_bin     !< override submission executable
        logical,       optional, intent(in)    :: base         !< .true. to route via base (not dispatch) controller
        type(chash) :: job_descr
        logical     :: l_base
        l_base = .false.
        if( present(base) ) l_base = base
        call cline%gen_job_descr(job_descr)
        if( present(exec_bin) ) then
            if( l_base ) then
                call self%base_qscripts%generate_script(job_descr, self%qdescr, exec_bin, script_name, outfile=outfile)
            else
                call self%qscripts%generate_script(job_descr, self%qdescr, exec_bin, script_name, outfile=outfile)
            end if
        else
            if( l_base ) then
                call self%base_qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, script_name, outfile=outfile)
            else
                call self%qscripts%generate_script(job_descr, self%qdescr, self%simple_exec_bin, script_name, outfile=outfile)
            end if
        end if
        call wait_for_closure(script_name)
        if( l_base ) then
            call self%base_qscripts%submit_script(script_name)
        else
            call self%qscripts%submit_script(script_name)
        end if
        call job_descr%kill
    end subroutine exec_simple_prg_in_queue_async

    !> Generate and submit multiple scripts as a single batch without blocking.
    !! When base=.true. all scripts are routed through base_qscripts.
    subroutine exec_simple_prgs_in_queue_async( self, clines, script_name, outfile, exec_bins, base )
        use simple_cmdline, only: cmdline
        class(qsys_env),            intent(inout) :: self
        type(cmdline), allocatable, intent(in)    :: clines(:)      !< one command line per job
        class(string),              intent(in)    :: script_name    !< shared script path prefix
        class(string),              intent(in)    :: outfile        !< shared stdout/stderr log path
        class(string),    optional, intent(in)    :: exec_bins(:)   !< per-job executable overrides
        logical,          optional, intent(in)    :: base           !< .true. to route via base controller
        type(chash), allocatable :: jobs_descr(:)
        integer :: i, njobs
        logical :: l_base
        njobs = size(clines)
        l_base = .false.
        if( present(base) ) l_base = base
        allocate(jobs_descr(njobs))
        do i = 1, njobs
            call clines(i)%gen_job_descr(jobs_descr(i))
        end do
        if( l_base ) then
            call self%base_qscripts%generate_script(jobs_descr, self%qdescr, self%simple_exec_bin, script_name, outfile, exec_bins=exec_bins)
        else
            call self%qscripts%generate_script(jobs_descr, self%qdescr, self%simple_exec_bin, script_name, outfile, exec_bins=exec_bins)
        end if
        call wait_for_closure(script_name)
        if( l_base ) then
            call self%base_qscripts%submit_script(script_name)
        else
            call self%qscripts%submit_script(script_name)
        end if
        do i = 1, njobs
            call jobs_descr(i)%kill
        end do
        deallocate(jobs_descr)
    end subroutine exec_simple_prgs_in_queue_async

    !> Submit one simple_persistent_worker process per slot through the base scheduler.
    !! Each worker is given the server host/port at launch so it can connect back.
    !! Scripts and logs are placed under WORKER_DIR for easy monitoring.
    subroutine start_persistent_workers( self )
        class(qsys_env), intent(inout) :: self
        type(cmdline)                  :: cline
        type(string)                   :: script_name, outfile, executable, host_ips
        integer                        :: iworker, ncunits, nthr, port
        if( .not. associated(persistent_worker%server) ) THROW_HARD('cannot start persistent workers without an allocated server;')
        ! Cache server properties so each worker script can connect back.
        ncunits  = persistent_worker%n_workers
        nthr     = persistent_worker%nthr_per_worker
        port     = persistent_worker%server%get_port()
        host_ips = persistent_worker%server%get_host_ips()
        write(*,*) 'STARTING # WORKERS: ', ncunits
        ! Ensure per-worker output directories exist before submitting scripts.
        call simple_mkdir(WORKER_DIR)
        call simple_mkdir(WORKER_DIR//STDERROUT_DIR)
        executable = 'simple_persistent_worker'
        ! Common connection parameters shared by every worker script.
        call cline%set('nthr',   int2str(nthr))
        call cline%set('port',   int2str(port))
        call cline%set('server', host_ips%to_char())
        ! Submit one script per worker slot via the base (non-dispatch) controller.
        do iworker = 1, ncunits
            call cline%set('worker_id', int2str(iworker))
            script_name = WORKER_DIR//WORKER_SCRIPT_PREFIX//int2str_pad(iworker, 4)//SCRIPT_EXT
            outfile     = WORKER_DIR//STDERROUT_DIR//int2str_pad(iworker, 4)//LOG_EXT
            call self%exec_simple_prg_in_queue_async(cline, script_name, outfile, executable, base=.true.)
        end do
        call cline%kill
    end subroutine start_persistent_workers

    !> Get fully resolved executable path used for generated scripts.
    function get_exec_bin( self ) result( exec_bin )
        class(qsys_env), intent(in) :: self
        type(string)                :: exec_bin
        exec_bin = self%simple_exec_bin
    end function get_exec_bin

    !> Return canonical base qsys backend name.
    function get_qsys( self ) result( qsys )
        class(qsys_env), intent(in)   :: self
        character(len=:), allocatable :: qsys
        select type(pbase_qsys => self%base_qsys)
            class is(qsys_local)
                qsys = 'local'
            class is(qsys_slurm)
                qsys = 'slurm'
            class is(qsys_lsf)
                qsys = 'lsf'
            class is(qsys_sge)
                qsys = 'sge'
            class is(qsys_pbs)
                qsys = 'pbs'
            class is(qsys_persistent_worker)
                qsys = 'persistent_worker'
        end select
    end function get_qsys

    !> Report currently available computing units from scheduler controller.
    integer function get_navail_computing_units( self )
        class(qsys_env), intent(in) :: self
        get_navail_computing_units = self%qscripts%get_ncomputing_units_avail()
    end function get_navail_computing_units

    !> Release scripts/controllers and reset internal state.
    subroutine kill( self )
        class(qsys_env) :: self
        if( self%existence ) then
            if( allocated(self%parts) ) deallocate(self%parts)
            call self%qscripts%kill
            call self%base_qscripts%kill
            call self%qdescr%kill
            call self%qsys_fac_base%kill
            call self%qsys_fac_dispatch%kill
            self%base_qsys => null()
            self%dispatch_qsys => null()
            self%existence = .false.
        end if
    end subroutine kill

end module simple_qsys_env
