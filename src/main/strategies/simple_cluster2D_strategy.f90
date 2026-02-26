module simple_cluster2D_strategy
use simple_core_module_api
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_commanders_api, only: commander_base
use simple_qsys_env,       only: qsys_env
use simple_convergence,    only: convergence
use simple_gui_utils,      only: mrc2jpeg_tiled
use simple_progress,       only: progressfile_update
implicit none

public :: cluster2D_strategy, cluster2D_inmem_strategy, cluster2D_distr_strategy, create_strategy
private
#include "simple_local_flags.inc"

!> Minimal strategy interface - only divergent operations
type, abstract :: cluster2D_strategy
contains
    procedure(init_interface),          deferred :: initialize
    procedure(probtab_interface),       deferred :: build_probability_table
    procedure(exec_iter_interface),     deferred :: execute_iteration
    procedure(finalize_iter_interface), deferred :: finalize_iteration
    procedure(finalize_run_interface),  deferred :: finalize_run
    procedure(cleanup_interface),       deferred :: cleanup
end type cluster2D_strategy

type, extends(cluster2D_strategy) :: cluster2D_inmem_strategy
    type(convergence) :: conv
contains
    procedure :: initialize              => inmem_initialize
    procedure :: build_probability_table => inmem_build_probability_table
    procedure :: execute_iteration       => inmem_execute_iteration
    procedure :: finalize_iteration      => inmem_finalize_iteration
    procedure :: finalize_run            => inmem_finalize_run
    procedure :: cleanup                 => inmem_cleanup
end type cluster2D_inmem_strategy

type, extends(cluster2D_strategy) :: cluster2D_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    integer        :: nthr_master
contains
    procedure :: initialize              => distr_initialize
    procedure :: build_probability_table => distr_build_probability_table
    procedure :: execute_iteration       => distr_execute_iteration
    procedure :: finalize_iteration      => distr_finalize_iteration
    procedure :: finalize_run            => distr_finalize_run
    procedure :: cleanup                 => distr_cleanup
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

    subroutine probtab_interface(self, params, cline, which_iter, xprob_tab2D, xprob_tab2D_distr)
        import :: cluster2D_strategy, parameters, cmdline, commander_base
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(cmdline),             intent(inout) :: cline
        integer,                   intent(in)    :: which_iter
        class(commander_base),     intent(inout) :: xprob_tab2D
        class(commander_base),     intent(inout) :: xprob_tab2D_distr
    end subroutine probtab_interface

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

    function create_strategy(params, cline) result(strategy)
        type(parameters), intent(in) :: params
        class(cmdline),   intent(in) :: cline
        class(cluster2D_strategy), allocatable :: strategy
        if( cline%defined('nparts') )then
            allocate(cluster2D_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(cluster2D_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_strategy

    ! ========================================================================
    ! SHARED-MEMORY IMPLEMENTATION
    ! ========================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        if( params%startit == 1 .and. (.not.str_has_substr(params%refine,'prob')) )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
    end subroutine inmem_initialize

    subroutine inmem_build_probability_table(self, params, cline, which_iter, xprob_tab2D, xprob_tab2D_distr)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: which_iter
        class(commander_base),           intent(inout) :: xprob_tab2D
        class(commander_base),           intent(inout) :: xprob_tab2D_distr
        type(cmdline) :: cline_prob_tab2D
        cline_prob_tab2D = cline
        call cline_prob_tab2D%set('prg',       'prob_tab2D')
        call cline_prob_tab2D%set('which_iter', which_iter)
        call cline_prob_tab2D%set('refs',       params%refs)
        call cline_prob_tab2D%set('frcs',       FRCS_FILE)
        call cline_prob_tab2D%set('startit',    which_iter)
        call cline_prob_tab2D%set('extr_iter',  params%extr_iter)
        call xprob_tab2D%execute(cline_prob_tab2D)
    end subroutine inmem_build_probability_table

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
        logical :: l_stream
        ! Execute alignment (cluster2D_exec handles everything: refs prep, alignment, cavgs)
        call cluster2D_exec(params, build, cline, which_iter, converged)
        ! Euclid sigma2 consolidation for next iteration (legacy sequencing)
        if( params%l_needs_sigma )then
            call cline%set('which_iter', which_iter + 1)
            call xcalc_group_sigmas%execute(cline)
            call cline%set('which_iter', which_iter)
        endif
        ! Check convergence
        converged = self%conv%check_conv2D(params, cline, build%spproj_field, &
                                           build%spproj_field%get_n('class'), params%msk)
        converged = converged .and. (which_iter >= params%minits)
        converged = converged .or.  (which_iter >= params%maxits)
        ! Update progress
        l_stream = trim(params%stream) .eq. 'yes'
        if( .not. l_stream ) call progressfile_update(self%conv%get('progress'))
        ! Write starfile
        call starproj%export_cls2D(build%spproj, which_iter)
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration(self, params, build, which_iter)
        class(cluster2D_inmem_strategy), intent(inout)  :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        integer,                         intent(in)    :: which_iter
        type(string) :: str_iter, cavg_mrc, cavg_jpg
        real, allocatable :: states(:)
        ! Update project segments
        call build%spproj%os_cls3D%new(params%ncls, is_ptcl=.false.)
        states = build%spproj%os_cls2D%get_all('state')
        call build%spproj%os_cls3D%set_all('state', states)
        call build%spproj%write_segment_inside('cls2D', params%projfile)
        call build%spproj%write_segment_inside('cls3D', params%projfile)
        deallocate(states)
        ! Generate JPEG
        str_iter = int2str_pad(which_iter, 3)
        cavg_mrc = string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//MRC_EXT
        cavg_jpg = string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//JPG_EXT
        if( file_exists(cavg_mrc) )then
            call mrc2jpeg_tiled(cavg_mrc, cavg_jpg)
            write(logfhandle,'(A,A)') '>>> JPEG ', cavg_jpg%to_char()
        endif
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
        type(string) :: finalcavgs
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            call build%spproj%write_segment_inside('out', params%projfile)
            if( .not. params%l_polar )then
                finalcavgs = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//params%ext%to_char()
                call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
            endif
        endif
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

    subroutine distr_build_probability_table(self, params, cline, which_iter, xprob_tab2D, xprob_tab2D_distr)
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(cmdline),                   intent(inout) :: cline
        integer,                         intent(in)    :: which_iter
        class(commander_base),           intent(inout) :: xprob_tab2D
        class(commander_base),           intent(inout) :: xprob_tab2D_distr
        type(cmdline) :: cline_prob_tab2D
        type(string)  :: refs_prob
        cline_prob_tab2D = cline
        call cline_prob_tab2D%set('prg',       'prob_tab2D')
        call cline_prob_tab2D%set('which_iter', which_iter)
        call cline_prob_tab2D%set('refs',       params%refs)
        call cline_prob_tab2D%set('frcs',       FRCS_FILE)
        call cline_prob_tab2D%set('startit',    which_iter)
        call cline_prob_tab2D%set('extr_iter',  params%extr_iter)
        call xprob_tab2D_distr%execute(cline_prob_tab2D)
    end subroutine distr_build_probability_table

    subroutine distr_execute_iteration(self, params, build, cline, which_iter, converged)
        use simple_stream_utils,         only: terminate_stream
        use simple_commanders_mkcavgs,   only: commander_cavgassemble
        use simple_commanders_euclid,    only: commander_calc_group_sigmas
        use simple_strategy2D_matcher,   only: cluster2D_exec
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                      intent(inout) :: params
        type(builder),                         intent(inout) :: build
        type(cmdline),                         intent(inout) :: cline
        integer,                               intent(in)    :: which_iter
        logical,                               intent(out)   :: converged
        type(commander_cavgassemble)      :: xcavgassemble
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline)                     :: cline_cavgassemble, cline_calc_sigma
        type(string)                      :: str_iter, refs, refs_even, refs_odd
        type(convergence)                 :: conv
        real                              :: frac_srch_space
        logical                           :: l_stream
        ! Update job description
        call self%job_descr%set('refs',      params%refs)
        call self%job_descr%set('startit',   int2str(which_iter))
        call self%job_descr%set('which_iter',int2str(which_iter))
        call self%job_descr%set('extr_iter', int2str(params%extr_iter))
        call self%job_descr%set('frcs',      FRCS_FILE)
        ! Schedule distributed jobs
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, &
                                                     algnfbody=string(ALGN_FBODY), &
                                                     array=L_USE_SLURM_ARR, &
                                                     extra_params=params)
        call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 1')
        ! Merge alignment docs
        call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
        ! Assemble class averages for this iteration
        if( trim(params%restore_cavgs) .eq. 'yes' )then
            str_iter   = int2str_pad(which_iter, 3)
            refs       = CAVGS_ITER_FBODY // str_iter%to_char()            // params%ext%to_char()
            refs_even  = CAVGS_ITER_FBODY // str_iter%to_char() // '_even' // params%ext%to_char()
            refs_odd   = CAVGS_ITER_FBODY // str_iter%to_char() // '_odd'  // params%ext%to_char()
            cline_cavgassemble = cline
            call cline_cavgassemble%set('prg',  'cavgassemble')
            if( cline_cavgassemble%defined('which_iter') ) call cline_cavgassemble%delete('which_iter')
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
        converged = conv%check_conv2D(params, cline, build%spproj_field, &
                                      build%spproj_field%get_n('class'), params%msk)
        l_stream = trim(params%stream) .eq. 'yes'
        if( .not. l_stream ) call progressfile_update(conv%get('progress'))
        frac_srch_space = 0.
        if( which_iter > 1 ) frac_srch_space = conv%get('frac_srch')
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
        type(string) :: str_iter, cavg_mrc, cavg_jpg
        ! Print correlations if requested
        if( trim(params%print_corrs) .eq. 'yes' )then
            call build%spproj_field%write(string('ptcl2D_'//int2str_pad(which_iter,2)//'.txt'))
        endif
        ! Generate JPEG
        if( trim(params%restore_cavgs) .eq. 'yes' )then
            str_iter = int2str_pad(which_iter, 3)
            cavg_mrc = string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//params%ext%to_char()
            cavg_jpg = string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//JPG_EXT
            if( file_exists(cavg_mrc) )then
                call mrc2jpeg_tiled(cavg_mrc, cavg_jpg)
                write(logfhandle,'(A,A)') '>>> JPEG ', cavg_jpg%to_char()
            endif
        endif
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
        type(string) :: finalcavgs
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            call build%spproj%write_segment_inside('out', params%projfile)
            if( .not. params%l_polar )then
                finalcavgs = CAVGS_ITER_FBODY//int2str_pad(last_iter,3)//params%ext%to_char()
                call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
            endif
        endif
    end subroutine distr_finalize_run

end module simple_cluster2D_strategy
