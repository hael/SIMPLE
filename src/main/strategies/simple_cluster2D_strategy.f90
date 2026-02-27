module simple_cluster2D_strategy
use simple_core_module_api
use simple_builder,     only: builder
use simple_parameters,  only: parameters
use simple_cmdline,     only: cmdline
use simple_qsys_env,    only: qsys_env
use simple_convergence, only: convergence
use simple_gui_utils,   only: mrc2jpeg_tiled
use simple_progress,    only: progressfile_update
implicit none

public :: cluster2D_strategy, cluster2D_inmem_strategy, cluster2D_distr_strategy, create_cluster2D_strategy
private
#include "simple_local_flags.inc"

!> Minimal strategy interface - only divergent operations
type, abstract :: cluster2D_strategy
contains
    procedure(init_interface),          deferred :: initialize
    procedure(exec_iter_interface),     deferred :: execute_iteration
    procedure(finalize_iter_interface), deferred :: finalize_iteration
    procedure(finalize_run_interface),  deferred :: finalize_run
    procedure(cleanup_interface),       deferred :: cleanup
end type cluster2D_strategy

type, extends(cluster2D_strategy) :: cluster2D_inmem_strategy
    type(convergence) :: conv
contains
    procedure :: initialize         => inmem_initialize
    procedure :: execute_iteration  => inmem_execute_iteration
    procedure :: finalize_iteration => inmem_finalize_iteration
    procedure :: finalize_run       => inmem_finalize_run
    procedure :: cleanup            => inmem_cleanup
end type cluster2D_inmem_strategy

type, extends(cluster2D_strategy) :: cluster2D_distr_strategy
    type(qsys_env)     :: qenv
    type(chash)        :: job_descr
    type(convergence)  :: conv
    integer            :: nthr_master
contains
    procedure :: initialize         => distr_initialize
    procedure :: execute_iteration  => distr_execute_iteration
    procedure :: finalize_iteration => distr_finalize_iteration
    procedure :: finalize_run       => distr_finalize_run
    procedure :: cleanup            => distr_cleanup
end type cluster2D_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: cluster2D_strategy, parameters, builder, cmdline
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_iter_interface(self, params, build, cline, which_iter, converged)
        import :: cluster2D_strategy, parameters, builder, cmdline
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        logical,                   intent(out)   :: converged
    end subroutine exec_iter_interface

    subroutine finalize_iter_interface(self, params, build, which_iter)
        import :: cluster2D_strategy, parameters, builder
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        integer,                   intent(in)    :: which_iter
    end subroutine finalize_iter_interface

    subroutine finalize_run_interface(self, params, build, cline, last_iter)
        import :: cluster2D_strategy, parameters, builder, cmdline
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
        integer,                   intent(in)    :: last_iter
    end subroutine finalize_run_interface

    subroutine cleanup_interface(self, params)
        import :: cluster2D_strategy, parameters
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    function create_cluster2D_strategy(params, cline) result(strategy)
        type(parameters), intent(in) :: params
        class(cmdline),   intent(in) :: cline
        class(cluster2D_strategy), allocatable :: strategy
        if( (params%nparts > 1) .and. (.not.cline%defined('part')) )then
            allocate(cluster2D_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(cluster2D_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_cluster2D_strategy

    ! ========================================================================
    ! SHARED-MEMORY IMPLEMENTATION
    ! ========================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer :: startit
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
    end subroutine inmem_initialize

    subroutine inmem_execute_iteration(self, params, build, cline, which_iter, converged)
        use simple_strategy2D_matcher, only: cluster2D_exec
        use simple_starproject,        only: starproject
        use simple_commanders_euclid,  only: commander_calc_group_sigmas
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: which_iter
        logical,                         intent(out)   :: converged
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(starproject) :: starproj
        if( .not. params%l_worker_distr )then
            call cline%set('startit',    which_iter)
            call cline%set('which_iter', which_iter)
            call cline%set('extr_iter',  params%extr_iter)
            call cline%set('outfile', ALGN_FBODY//int2str_pad(params%part,params%numlen)//METADATA_EXT)
        endif
        ! Execute alignment (cluster2D_exec handles everything: refs prep, alignment, cavgs)
        call cluster2D_exec(params, build, cline, which_iter, converged)
        if( params%l_worker_distr )then
            converged = .true. ! Worker subprocesses execute one iteration and exit
            return
        endif
        ! Euclid sigma2 consolidation for next iteration
        if( params%l_needs_sigma )then
            call cline%set('which_iter', which_iter + 1)
            call xcalc_group_sigmas%execute(cline)
            call cline%set('which_iter', which_iter)
        endif
        ! Write starfile
        call starproj%export_cls2D(build%spproj, which_iter)
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration(self, params, build, which_iter)
        class(cluster2D_inmem_strategy), intent(inout)  :: self
        type(parameters),                intent(in)     :: params
        type(builder),                   intent(inout)  :: build
        integer,                         intent(in)     :: which_iter
        if( .not. params%l_worker_distr ) call gen_jpeg( which_iter ) 
    end subroutine inmem_finalize_iteration

    subroutine inmem_cleanup(self, params)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        ! No cleanup needed for in-memory strategy, but could add here if needed in the future
    end subroutine inmem_cleanup

    subroutine inmem_finalize_run(self, params, build, cline, last_iter)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: last_iter
        if( .not. params%l_worker_distr ) call build%spproj%write_segment_inside(params%oritype, params%projfile)
    end subroutine inmem_finalize_run

    ! ========================================================================
    ! DISTRIBUTED IMPLEMENTATION
    ! ========================================================================

    subroutine distr_initialize(self, params, build, cline)
        use simple_exec_helpers, only: set_master_num_threads
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        call set_master_num_threads(self%nthr_master, string('CLUSTER2D'))
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
        call build%spproj%split_stk(params%nparts)
    end subroutine distr_initialize

    subroutine distr_execute_iteration(self, params, build, cline, which_iter, converged)
        use simple_stream_utils,         only: terminate_stream
        use simple_commanders_mkcavgs,   only: commander_cavgassemble
        use simple_commanders_euclid,    only: commander_calc_group_sigmas
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: which_iter
        logical,                         intent(out)   :: converged
        type(commander_cavgassemble)      :: xcavgassemble
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline)                     :: cline_cavgassemble, cline_calc_sigma
        type(string)                      :: str_iter, refs, refs_even, refs_odd
        real                              :: frac_srch_space
        logical                           :: l_stream
        ! Update job description
        call cline%set('nparts',     params%nparts)
        call cline%set('startit',    which_iter)
        call cline%set('which_iter', which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        call self%job_descr%set('refs',       params%refs)
        call self%job_descr%set('nparts',     int2str(params%nparts))
        call self%job_descr%set('startit',    int2str(which_iter))
        call self%job_descr%set('which_iter', int2str(which_iter))
        call self%job_descr%set('extr_iter',  int2str(params%extr_iter))
        call self%job_descr%set('frcs',       FRCS_FILE)
        ! Schedule distributed jobs
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, &
                                                     algnfbody=string(ALGN_FBODY), &
                                                     array=L_USE_SLURM_ARR, &
                                                     extra_params=params)
        call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 1')
        ! Merge alignment docs
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
        ! Assemble class averages
        if( trim(params%restore_cavgs) .eq. 'yes' )then
            str_iter   = int2str_pad(which_iter, 3)
            refs       = CAVGS_ITER_FBODY // str_iter%to_char()            // MRC_EXT
            refs_even  = CAVGS_ITER_FBODY // str_iter%to_char() // '_even' // MRC_EXT
            refs_odd   = CAVGS_ITER_FBODY // str_iter%to_char() // '_odd'  // MRC_EXT
            cline_cavgassemble = cline
            call cline_cavgassemble%set('prg',  'cavgassemble')
            call cline_cavgassemble%delete('which_iter')
            call cline_cavgassemble%set('refs', refs)
            call cline_cavgassemble%set('nthr', self%nthr_master)
            call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 2')
            call xcavgassemble%execute(cline_cavgassemble)
            params%refs      = refs
            params%refs_even = refs_even
            params%refs_odd  = refs_odd
            call cline%set('refs', refs)
        endif
        ! Sigma2 consolidation
        if( params%l_needs_sigma )then
            cline_calc_sigma = cline
            call cline_calc_sigma%set('prg',        'calc_group_sigmas')
            call cline_calc_sigma%set('which_iter', which_iter + 1)
            call cline_calc_sigma%set('nthr',       self%nthr_master)
            call xcalc_group_sigmas%execute(cline_calc_sigma)
        endif
        ! Check convergence
        converged = self%conv%check_conv2D(params, cline, build%spproj_field, &
                                           build%spproj_field%get_n('class'), params%msk)
        l_stream = trim(params%stream) .eq. 'yes'
        if( .not. l_stream ) call progressfile_update(self%conv%get('progress'))
        frac_srch_space = 0.
        if( which_iter > 1 ) frac_srch_space = self%conv%get('frac_srch')
        ! Activate shift search if needed
        if( which_iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. params%l_doshift) )then
            if( .not. self%job_descr%isthere('trs') )then
                call self%job_descr%set('trs', real2str(params%trs))
            endif
        endif
        converged = (which_iter >= params%minits) .and. converged
        converged = converged .or. (which_iter >= params%maxits)
    end subroutine distr_execute_iteration

    subroutine distr_finalize_iteration(self, params, build, which_iter)
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        integer,                         intent(in)    :: which_iter
        type(string) :: str_iter, fbody
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call gen_jpeg(which_iter)
    end subroutine distr_finalize_iteration

    subroutine distr_cleanup(self, params)
        use simple_qsys_funs, only: qsys_cleanup
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        call qsys_cleanup(params)
        call self%job_descr%kill
    end subroutine distr_cleanup

    subroutine distr_finalize_run(self, params, build, cline, last_iter)
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: last_iter
    end subroutine distr_finalize_run

    ! private helpers

    subroutine gen_jpeg( which_iter ) 
        integer, intent(in) :: which_iter
        type(string) :: str_iter, fbody, mrc, jpg
        ! Generate JPEG
        str_iter = int2str_pad(which_iter, 3)
        fbody = string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()
        mrc   = fbody//MRC_EXT
        jpg   = fbody//JPG_EXT
        if( file_exists(mrc) )then
            call mrc2jpeg_tiled(mrc, jpg)
            write(logfhandle,'(A,A)') '>>> JPEG ', jpg%to_char()
        endif
        call str_iter%kill
        call fbody%kill
        call mrc%kill
        call jpg%kill
    end subroutine gen_jpeg

end module simple_cluster2D_strategy