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
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_iter_interface(self, params, build, cline, converged)
        import :: cluster2D_strategy, parameters, builder, cmdline
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
        logical,                   intent(out)   :: converged
    end subroutine exec_iter_interface

    subroutine finalize_iter_interface(self, params, build)
        import :: cluster2D_strategy, parameters, builder
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
    end subroutine finalize_iter_interface

    subroutine finalize_run_interface(self, params, build, cline)
        import :: cluster2D_strategy, parameters, builder, cmdline
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        type(cmdline),             intent(inout) :: cline
    end subroutine finalize_run_interface

    subroutine cleanup_interface(self, params)
        import :: cluster2D_strategy, parameters
        class(cluster2D_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
    end subroutine cleanup_interface
end interface

contains

    !> Strategy selection based on command-line shape.
    function create_cluster2D_strategy(cline) result(strategy)
        class(cmdline),   intent(in) :: cline
        class(cluster2D_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not.cline%defined('part')) )then
            allocate(cluster2D_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(cluster2D_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_cluster2D_strategy

    ! ======================================================================
    ! SHARED-MEMORY STRATEGY METHODS
    ! ======================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        integer :: startit
        call cline%set('outfile', ALGN_FBODY//int2str_pad(params%part,params%numlen)//METADATA_EXT)
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call cline%set('mkdir', 'no')
        params%which_iter = max(1, params%startit)
        call init_cluster2D_refs(cline, params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( startit == 1 )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
    end subroutine inmem_initialize

    subroutine inmem_execute_iteration( self, params, build, cline, converged)
        use simple_strategy2D_matcher, only: cluster2D_exec
        use simple_starproject,        only: starproject
        use simple_commanders_euclid,  only: commander_calc_group_sigmas
        use simple_commanders_prob,    only: commander_prob_align2D
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        logical,                         intent(out)   :: converged
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_prob_align2D)      :: xprob_align2D
        type(starproject) :: starproj
        type(cmdline)     :: cline_prob_align
        call self%conv%print_iteration(params%which_iter)
        call cline%set('startit',    params%startit)
        call cline%set('which_iter', params%which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            call cline_prob_align%set('prg', 'prob_align2D')
            call cline_prob_align%set('which_iter', params%which_iter)
            call cline_prob_align%set('startit',    params%startit)
            call build%spproj%write_segment_inside(params%oritype)
            call xprob_align2D%execute(cline_prob_align)
            call build%spproj%read_segment(params%oritype, params%projfile)
        endif
        ! main clustering/alignment step
        call cluster2D_exec(params, build, cline, params%which_iter, converged)
        ! Euclid sigma2 consolidation for next iteration
        if( params%cc_objfun==OBJFUN_EUCLID )then
            call cline%set('which_iter', params%which_iter + 1)
            call xcalc_group_sigmas%execute(cline)
            call cline%set('which_iter', params%which_iter)
        endif
        ! Write starfile
        call starproj%export_cls2D(build%spproj, params%which_iter)
    end subroutine inmem_execute_iteration

    subroutine inmem_finalize_iteration( self, params, build)
        class(cluster2D_inmem_strategy), intent(inout)  :: self
        type(parameters),                intent(in)     :: params
        type(builder),                   intent(inout)  :: build
        call gen_jpeg(params%which_iter) 
    end subroutine inmem_finalize_iteration

    subroutine inmem_cleanup( self, params)
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        ! No cleanup needed for in-memory strategy, but could add here if needed in the future
    end subroutine inmem_cleanup

    subroutine inmem_finalize_run( self, params, build, cline )
        class(cluster2D_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        type(string) :: finalcavgs
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            if( .not. params%l_polar )then
                finalcavgs = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//MRC_EXT
                call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
                call finalcavgs%kill
            endif
            call build%spproj%write_segment_inside('out', params%projfile)
        endif
        call cline%set('endit', params%which_iter)
    end subroutine inmem_finalize_run

    ! ======================================================================
    ! DISTRIBUTED STRATEGY METHODS
    ! ======================================================================

    subroutine distr_initialize( self, params, build, cline )
        use simple_exec_helpers, only: set_master_num_threads
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call cline%set('mkdir', 'no')
        params%which_iter = max(1, params%startit)
        call init_cluster2D_refs(cline, params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        if( params%startit == 1 )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        call set_master_num_threads(self%nthr_master, string('CLUSTER2D'))
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
        call build%spproj%split_stk(params%nparts)
    end subroutine distr_initialize

    subroutine distr_execute_iteration( self, params, build, cline, converged )
        use simple_stream_utils,         only: terminate_stream
        use simple_commanders_mkcavgs,   only: commander_cavgassemble
        use simple_commanders_euclid,    only: commander_calc_group_sigmas
        use simple_commanders_prob,      only: commander_prob_align2D
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        logical,                         intent(out)   :: converged
        type(commander_cavgassemble)      :: xcavgassemble
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_prob_align2D)      :: xprob_align2D
        type(cmdline)                     :: cline_cavgassemble, cline_calc_sigma, cline_prob_align
        type(string)                      :: str_iter
        real                              :: frac_srch_space
        call self%conv%print_iteration(params%which_iter)
        ! Update job description
        call cline%set('nparts',     params%nparts)
        call cline%set('startit',    params%startit)
        call cline%set('which_iter', params%which_iter)
        call cline%set('extr_iter',  params%extr_iter)
        call self%job_descr%set('refs',       params%refs)
        call self%job_descr%set('nparts',     int2str(params%nparts))
        call self%job_descr%set('startit',    int2str(params%startit))
        call self%job_descr%set('which_iter', int2str(params%which_iter))
        call self%job_descr%set('extr_iter',  int2str(params%extr_iter))
        call self%job_descr%set('frcs',       FRCS_FILE)
        if( params%l_prob_align_mode )then
            cline_prob_align = cline
            call cline_prob_align%set('prg', 'prob_align2D')
            call cline_prob_align%set('which_iter', params%which_iter)
            call cline_prob_align%set('startit',    params%startit)
            call build%spproj%write_segment_inside(params%oritype)
            call xprob_align2D%execute(cline_prob_align)
            call build%spproj%read_segment(params%oritype, params%projfile)
        endif
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
            str_iter   = int2str_pad(params%which_iter, 3)
            params%refs       = CAVGS_ITER_FBODY // str_iter%to_char()            // MRC_EXT
            params%refs_even  = CAVGS_ITER_FBODY // str_iter%to_char() // '_even' // MRC_EXT
            params%refs_odd   = CAVGS_ITER_FBODY // str_iter%to_char() // '_odd'  // MRC_EXT
            call cline%set('refs', params%refs)
            cline_cavgassemble = cline
            call cline_cavgassemble%set('prg',  'cavgassemble')
            call cline_cavgassemble%delete('which_iter')
            call cline_cavgassemble%set('refs', params%refs)
            call cline_cavgassemble%set('nthr', self%nthr_master)
            call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 2')
            call xcavgassemble%execute(cline_cavgassemble)
        endif
        ! Sigma2 consolidation
        if( params%cc_objfun==OBJFUN_EUCLID )then
            cline_calc_sigma = cline
            call cline_calc_sigma%set('prg',        'calc_group_sigmas')
            call cline_calc_sigma%set('which_iter', params%which_iter + 1)
            call cline_calc_sigma%set('nthr',       self%nthr_master)
            call xcalc_group_sigmas%execute(cline_calc_sigma)
        endif
        ! Check convergence
        converged = self%conv%check_conv2D(params, cline, build%spproj_field, &
                                           build%spproj_field%get_n('class'), params%msk)
        if( trim(params%stream2d).eq.'no' ) call progressfile_update(self%conv%get('progress'))
        frac_srch_space = 0.
        if( params%which_iter > 1 ) frac_srch_space = self%conv%get('frac_srch')
        ! Activate shift search if needed
        if( params%which_iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. params%l_doshift) )then
            if( .not. self%job_descr%isthere('trs') )then
                call self%job_descr%set('trs', real2str(params%trs))
            endif
        endif
        converged = (params%which_iter >= params%minits) .and. converged
        converged = converged .or. (params%which_iter >= params%maxits)
    end subroutine distr_execute_iteration

    subroutine distr_finalize_iteration( self, params, build )
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        call gen_jpeg(params%which_iter)
    end subroutine distr_finalize_iteration

    subroutine distr_cleanup( self, params )
        use simple_qsys_funs, only: qsys_cleanup
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        call qsys_cleanup(params)
        call self%job_descr%kill
    end subroutine distr_cleanup

    subroutine distr_finalize_run( self, params, build, cline )
        class(cluster2D_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        type(cmdline),                   intent(inout) :: cline
        type(string) :: finalcavgs
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) )then
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            endif
            if( .not. params%l_polar )then
                finalcavgs = CAVGS_ITER_FBODY//int2str_pad(params%which_iter,3)//MRC_EXT
                call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
                call finalcavgs%kill
            endif
            call build%spproj%write_segment_inside('out', params%projfile)
        endif
        call cline%set('endit', params%which_iter)
    end subroutine distr_finalize_run

    ! private helpers

    !> Initialize references for cluster2D (used by inmem and distr modes).
    subroutine init_cluster2D_refs( cline, params, build )
        use simple_procimgstk,         only: copy_imgfile
        use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr
        use simple_commanders_imgops,  only: commander_scale
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(commander_scale) :: xscale
        type(cmdline)         :: cline_make_cavgs, cline_scalerefs
        type(string)          :: refs_sc
        logical               :: l_scale_inirefs
        if( cline%defined('refs') ) return
        cline_make_cavgs = cline
        params%refs      = 'start2Drefs'     //MRC_EXT
        params%refs_even = 'start2Drefs_even'//MRC_EXT
        params%refs_odd  = 'start2Drefs_odd' //MRC_EXT
        l_scale_inirefs  = .false.
        if( build%spproj%is_virgin_field('ptcl2D') .or. params%which_iter <= 1 )then
            if( params%tseries .eq. 'yes' )then
                call init_tseries_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            else
                call init_standard_refs(cline, params, build, cline_make_cavgs, l_scale_inirefs)
            endif
        else
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        endif
        if( l_scale_inirefs )then
            refs_sc = 'refs'//SCALE_SUFFIX//MRC_EXT
            call cline_scalerefs%set('stk',    params%refs)
            call cline_scalerefs%set('outstk', refs_sc)
            call cline_scalerefs%set('smpd',   params%smpd)
            call cline_scalerefs%set('newbox', params%box_crop)
            call xscale%execute(cline_scalerefs)
            call simple_rename(refs_sc, params%refs)
        endif
        call copy_imgfile(params%refs, params%refs_even, params%smpd_crop, [1,params%ncls])
        call copy_imgfile(params%refs, params%refs_odd,  params%smpd_crop, [1,params%ncls])
        call cline%set('refs', params%refs)
        call cline_make_cavgs%kill
        call cline_scalerefs%kill
        call refs_sc%kill
    end subroutine init_cluster2D_refs

    subroutine init_tseries_refs( cline, params, build, cline_make_cavgs, l_scale_inirefs )
        use simple_procimgstk, only: selection_from_tseries_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: cnt, iptcl, ptclind
        if( cline%defined('nptcls_per_cls') )then
            if( build%spproj%os_ptcl2D%any_state_zero() )then
                THROW_HARD('cluster2D_nano does not allow state=0 particles, prune project before execution')
            endif
            cnt = 0
            do iptcl=1,params%nptcls,params%nptcls_per_cls
                cnt = cnt + 1
                params%ncls = cnt
                do ptclind=iptcl,min(params%nptcls, iptcl + params%nptcls_per_cls - 1)
                    call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                end do
            end do
            call cline%set('ncls', params%ncls)
            call cline_make_cavgs%set('ncls', params%ncls)
            call cline_make_cavgs%set('refs', params%refs)
            call execute_make_cavgs(cline_make_cavgs, cline, params)
            l_scale_inirefs = .false.
        else
            if( trim(params%refine).eq.'inpl' )then
                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                call cline%set('ncls', params%ncls)
                call cline_make_cavgs%set('ncls', params%ncls)
                call cline_make_cavgs%delete('tseries')
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            else
                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            endif
        endif
    end subroutine init_tseries_refs

    subroutine init_standard_refs( cline, params, build, cline_make_cavgs, l_scale_inirefs )
        use simple_procimgstk, only: random_selection_from_imgfile, noise_imgfile
        class(cmdline),   intent(inout) :: cline
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(cmdline),    intent(inout) :: cline_make_cavgs
        logical,          intent(out)   :: l_scale_inirefs
        integer :: iptcl
        select case(trim(params%cls_init))
            case('ptcl')
                call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                l_scale_inirefs = .true.
            case('rand')
                call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                l_scale_inirefs = .false.
            case('randcls')
                if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS')
                do iptcl=1,params%nptcls
                    if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                    call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                    call build%spproj_field%set(iptcl, 'w',     1.0)
                    call build%spproj_field%e3set(iptcl,ran3()*360.0)
                end do
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                call cline_make_cavgs%set('refs', params%refs)
                call execute_make_cavgs(cline_make_cavgs, cline, params)
                l_scale_inirefs = .false.
            case DEFAULT
                THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
        end select
    end subroutine init_standard_refs

    subroutine execute_make_cavgs( cline_make_cavgs, cline, params)
        use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr
        type(cmdline),    intent(inout) :: cline_make_cavgs
        class(cmdline),   intent(in)    :: cline
        type(parameters), intent(in)    :: params
        type(commander_make_cavgs)       :: xmake_cavgs
        type(commander_make_cavgs_distr) :: xmake_cavgs_distr
        call cline_make_cavgs%set('prg', 'make_cavgs')
        if( (params%nparts > 1) .and. (.not.cline%defined('part')) )then
            call xmake_cavgs_distr%execute(cline_make_cavgs)
        else
            call xmake_cavgs%execute(cline_make_cavgs)
        endif
    end subroutine execute_make_cavgs

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
