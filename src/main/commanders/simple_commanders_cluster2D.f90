!@descr: simultaneous 2D alignment and clustering of single-particle images
module simple_commanders_cluster2D
use simple_commanders_api
use simple_pftc_srch_api
use simple_classaverager
use simple_commanders_cavgs,   only: commander_rank_cavgs
use simple_commanders_mkcavgs, only: commander_make_cavgs, commander_make_cavgs_distr, commander_cavgassemble
use simple_commanders_euclid,  only: commander_calc_pspec_distr, commander_calc_group_sigmas
use simple_gui_utils,          only: mrc2jpeg_tiled
use simple_procimgstk,         only: selection_from_tseries_imgfile, random_selection_from_imgfile, copy_imgfile, noise_imgfile
use simple_progress,           only: progressfile_init, progressfile_update
use simple_commanders_imgops,  only: commander_scale
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cluster2D_distr
  contains
    procedure :: execute      => exec_cluster2D_distr
end type commander_cluster2D_distr

type, extends(commander_base) :: commander_cluster2D
  contains
    procedure :: execute      => exec_cluster2D
end type commander_cluster2D

type, extends(commander_base) :: commander_prob_tab2D
  contains
    procedure :: execute      => exec_prob_tab2D
end type commander_prob_tab2D

type, extends(commander_base) :: commander_prob_tab2D_distr
  contains
    procedure :: execute      => exec_prob_tab2D_distr
end type commander_prob_tab2D_distr

type, extends(commander_base) :: commander_ppca_denoise_classes
  contains
    procedure :: execute      => exec_ppca_denoise_classes
end type commander_ppca_denoise_classes

contains

    subroutine exec_cluster2D_unified(self, cline)
        use simple_cluster2D_strategy
        use simple_cluster2D_common
        class(commander_base),   intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)                       :: params
        type(builder)                          :: build
        class(cluster2D_strategy), allocatable :: strategy
        type(simple_nice_communicator)         :: nice_communicator
        type(string)                           :: finalcavgs
        logical                                :: converged, l_stream, l_worker_distr
        integer                                :: iter
        ! Initialize
        call cline%set('prg','cluster2D')

        !!!!!!!!!!!!! Replace this when this becomes the only workflow for cluster2D
        call set_cluster2D_defaults(cline)

        
        l_stream       = .false.
        if( cline%defined('stream') ) l_stream = cline%get_carg('stream')=='yes'
        if( l_stream ) call cline%set('stream','no')
        call params%new(cline)
        ! flag subprocess executed through simple_private_exec 
        params%l_worker_distr = cline%defined('part')
        if( str_has_substr(params%refine, 'prob') )then
            THROW_HARD('REFINE=prob* is unsupported in cluster2D')
        endif
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        call build%build_strategy2D_tbox(params)
        if( l_stream ) call cline%set('stream','yes')
        if( build%spproj%get_nptcls() == 0 ) THROW_HARD('no particles found!')
        call handle_objfun(params, cline)
        call cline%set('mkdir', 'no')
        if( .not. params%l_worker_distr ) then
            ! Nice communicator
            call nice_communicator%init(params%niceprocid, params%niceserver)
            nice_communicator%stat_root%stage = "initialising"
            call nice_communicator%cycle()
            if(cline%defined("niceserver")) call cline%delete('niceserver')
            if(cline%defined("niceprocid")) call cline%delete('niceprocid')
        endif
        ! This needs to be before init_cluster2D_refs since it uses which_iter
        params%which_iter = max(1, params%startit)
        iter = params%which_iter - 1
        call init_cluster2D_refs(cline, params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! Create strategy
        strategy = create_cluster2D_strategy(params, cline)
        call strategy%initialize(params, build, cline)
        if(.not. l_stream) call progressfile_init()
        ! Main loop
        do
            iter = iter + 1
            params%which_iter = iter
            params%extr_iter  = params%which_iter
            call cline%set('which_iter', int2str(params%which_iter))
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
            write(logfhandle,'(A)')   '>>>'
            if( .not. params%l_worker_distr ) then
                nice_communicator%stat_root%stage = "iteration " // int2str(params%which_iter)
                call nice_communicator%cycle()
            endif
            ! Strategy handles everything: alignment + cavgs + convergence
            call strategy%execute_iteration(params, build, cline, iter, converged)
            call strategy%finalize_iteration(params, build, iter)
            if( converged .or. iter >= params%maxits ) exit
        end do
        if( .not. params%l_worker_distr )then
            ! Cleanup
            nice_communicator%stat_root%stage = "terminating"
            call nice_communicator%cycle()
            call strategy%finalize_run(params, build, cline, iter)
            ! At the end:
            if( trim(params%restore_cavgs).eq.'yes' )then
                if( file_exists(FRCS_FILE) )then
                    call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
                endif
                call build%spproj%write_segment_inside('out', params%projfile)
                if( .not. params%l_polar )then
                    finalcavgs = CAVGS_ITER_FBODY//int2str_pad(iter,3)//MRC_EXT
                    call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
                endif
            endif
            call strategy%cleanup(params)
            call cline%set('endit', iter)
            call nice_communicator%terminate()
            call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        endif
        if (allocated(strategy)) deallocate(strategy)
        call build%kill_general_tbox()
        call build%kill_strategy2D_tbox()
        call simple_touch(CLUSTER2D_FINISHED)
    end subroutine exec_cluster2D_unified

    ! Replace exec_cluster2D
    ! subroutine exec_cluster2D(self, cline)
    !     class(commander_cluster2D), intent(inout) :: self
    !     class(cmdline),             intent(inout) :: cline
    !     call exec_cluster2D_unified(self, cline)
    ! end subroutine exec_cluster2D

    ! Replace exec_cluster2D_distr
    ! subroutine exec_cluster2D_distr(self, cline)
    !     class(commander_cluster2D_distr), intent(inout) :: self
    !     class(cmdline),                   intent(inout) :: cline
    !     call exec_cluster2D_unified(self, cline)
    ! end subroutine exec_cluster2D_distr

    subroutine exec_cluster2D_distr( self, cline )
        class(commander_cluster2D_distr), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(commander_make_cavgs_distr)  :: xmake_cavgs
        type(commander_scale)             :: xscale
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_cavgassemble)      :: xcavgassemble
        type(commander_prob_tab2D_distr)  :: xprob_tab2D_distr
         type(simple_nice_communicator)   :: nice_communicator
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        type(cmdline) :: cline_calc_sigma
        type(cmdline) :: cline_scalerefs
        type(cmdline) :: cline_prob_tab2D_distr
        integer(timer_int_kind)   :: t_init,   t_scheduled,  t_merge_algndocs,  t_cavgassemble,  t_tot
        real(timer_int_kind)      :: rt_init, rt_scheduled, rt_merge_algndocs, rt_cavgassemble, rt_tot
        type(string)              :: benchfname
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        type(chash)               :: job_descr
        type(string)              :: refs, refs_even, refs_odd, str, str_iter, finalcavgs, refs_sc
        real                      :: frac_srch_space
        integer                   :: nthr_here, iter, cnt, iptcl, ptclind, fnr
        logical                   :: l_stream, l_converged, l_scale_inirefs
        call cline%set('prg','cluster2D')
        call set_cluster2D_defaults( cline )
        ! streaming
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = cline%get_carg('stream')=='yes'
        endif
        call cline%set('stream','no') ! for parameters determination
        ! deal with # threads for the master process
        call set_master_num_threads(nthr_here, string('CLUSTER2D'))
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        nice_communicator%stat_root%stage = "initialising"
        call nice_communicator%cycle()
        if(cline%defined("niceserver")) call cline%delete('niceserver')
        if(cline%defined("niceprocid")) call cline%delete('niceprocid')
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        if( l_stream ) call cline%set('stream','yes')
        ! objective functions
        select case( params%cc_objfun )
        case( OBJFUN_CC )
            ! making sure euclid options are turned off
            params%l_needs_sigma = .false.
            params%needs_sigma   = 'no'
            call cline%set('ml_reg','no')
            params%ml_reg        = 'no'
            params%l_ml_reg      = .false.
        case( OBJFUN_EUCLID )
            ! all set
        end select
        ! sanity check
        if( build%spproj%get_nptcls() == 0 )then
            THROW_HARD('no particles found! exec_cluster2D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! splitting
        call build%spproj%split_stk(params%nparts)
        ! prepare command lines from prototype master
        cline_check_2Dconv      = cline
        cline_cavgassemble      = cline
        cline_make_cavgs        = cline ! ncls is transferred here
        cline_calc_sigma        = cline
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_make_cavgs%set(  'prg', 'make_cavgs')
        call cline_calc_sigma%set(  'prg', 'calc_group_sigmas')
        ! Initial references
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs'//params%ext%to_char()
            params%refs      = refs
            params%refs_even = 'start2Drefs_even'//params%ext%to_char()
            params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
            l_scale_inirefs  = .false.
            if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                if( params%tseries .eq. 'yes' )then
                    if( cline%defined('nptcls_per_cls') )then
                        if( build%spproj%os_ptcl2D%any_state_zero() )then
                            THROW_HARD('cluster2D_nano does not allow state=0 particles, prune project before execution; exec_cluster2D_distr')
                        endif
                        cnt = 0
                        do iptcl=1,params%nptcls,params%nptcls_per_cls
                            cnt = cnt + 1
                            params%ncls = cnt
                            do ptclind=iptcl,min(params%nptcls, iptcl + params%nptcls_per_cls - 1)
                                call build%spproj%os_ptcl2D%set(ptclind, 'class', cnt)
                            end do
                        end do
                        call job_descr%set('ncls',int2str(params%ncls))
                        call cline%set('ncls', params%ncls)
                        call cline_make_cavgs%set('ncls', params%ncls)
                        call cline_cavgassemble%set('ncls', params%ncls)
                        call cline_make_cavgs%set('refs', params%refs)
                        call xmake_cavgs%execute(cline_make_cavgs)
                        l_scale_inirefs = .false.
                    else
                        if( trim(params%refine).eq.'inpl' )then
                            params%ncls = build%spproj%os_ptcl2D%get_n('class')
                            call job_descr%set('ncls',int2str(params%ncls))
                            call cline%set('ncls', params%ncls)
                            call cline_make_cavgs%set('ncls', params%ncls)
                            call cline_make_cavgs%delete('tseries')
                            call cline_cavgassemble%set('ncls', params%ncls)
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs = .false.
                        else
                            call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                            l_scale_inirefs = .true.
                        endif
                    endif
                else
                    select case(trim(params%cls_init))
                        case('ptcl')
                            ! initialization from raw images
                            call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                            l_scale_inirefs = .true.
                        case('rand')
                            ! from noise
                            call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                            l_scale_inirefs = .false.
                        case('randcls')
                            if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS')
                            ! initialization from random classes
                            do iptcl=1,params%nptcls
                                if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                                call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                                call build%spproj_field%set(iptcl, 'w',     1.0)
                                call build%spproj_field%e3set(iptcl,ran3()*360.0)
                            end do
                            call build%spproj%write_segment_inside(params%oritype, params%projfile)
                            call cline_make_cavgs%set('refs', params%refs)
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs = .false.
                        case DEFAULT
                            THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
                    end select
                endif
            else
                call cline_make_cavgs%set('refs', params%refs)
                call xmake_cavgs%execute(cline_make_cavgs)
                l_scale_inirefs = .false.
            endif
            ! scale references to box_crop
            if( l_scale_inirefs )then
                refs_sc = 'refs'//SCALE_SUFFIX//params%ext%to_char()
                call cline_scalerefs%set('stk',    params%refs)
                call cline_scalerefs%set('outstk', refs_sc)
                call cline_scalerefs%set('smpd',   params%smpd)
                call cline_scalerefs%set('newbox', params%box_crop)
                call cline_scalerefs%set('nthr',   nthr_here)
                call xscale%execute(cline_scalerefs)
                call simple_rename(refs_sc, params%refs)
            endif
            call copy_imgfile(params%refs, params%refs_even, params%smpd_crop, [1,params%ncls])
            call copy_imgfile(params%refs, params%refs_odd,  params%smpd_crop, [1,params%ncls])
        else
            refs = params%refs
        endif
        ! variable neighbourhood size
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! deal with eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! initialise progress monitor
        if(.not. l_stream) call progressfile_init()
        ! main loop
        iter = params%startit - 1
        do
            if( L_BENCH_GLOB )then
                t_init = tic()
                t_tot  = t_init
            endif
            iter = iter + 1
            params%which_iter = iter
            call cline%set(     'which_iter', int2str(params%which_iter))
            call job_descr%set( 'which_iter', int2str(params%which_iter))
            str_iter = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
            write(logfhandle,'(A)')   '>>>'
            ! nice
            nice_communicator%stat_root%stage = "iteration " // int2str(params%which_iter)
            call nice_communicator%cycle()
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', int2str(params%extr_iter))
            call cline%set('extr_iter', params%extr_iter)
            ! build probability table
            if( str_has_substr(params%refine, 'prob') )then
                cline_prob_tab2D_distr = cline
                call cline_prob_tab2D_distr%set('refs',      refs)
                call cline_prob_tab2D_distr%set('frcs',      FRCS_FILE)
                call cline_prob_tab2D_distr%set('startit',   iter)
                call cline_prob_tab2D_distr%set('extr_iter', params%extr_iter)
                call xprob_tab2D_distr%execute(cline_prob_tab2D_distr)
            endif
            ! updates
            call job_descr%set('refs', refs)
            call job_descr%set('startit', int2str(iter))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', FRCS_FILE)
            ! schedule
            if( L_BENCH_GLOB )then
                rt_init = toc(t_init)
                t_scheduled = tic()
            endif
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR, extra_params=params)
            call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 1')
            ! assemble alignment docs
            if( L_BENCH_GLOB )then
                rt_scheduled = toc(t_scheduled)
                t_merge_algndocs = tic()
            endif
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            if( L_BENCH_GLOB )then
                rt_merge_algndocs = toc(t_merge_algndocs)
                t_cavgassemble = tic()
            endif
            ! assemble class averages
            if( trim(params%restore_cavgs).eq.'yes' )then
                refs      = CAVGS_ITER_FBODY // str_iter%to_char()            // params%ext%to_char()
                refs_even = CAVGS_ITER_FBODY // str_iter%to_char() // '_even' // params%ext%to_char()
                refs_odd  = CAVGS_ITER_FBODY // str_iter%to_char() // '_odd'  // params%ext%to_char()
                call cline_cavgassemble%set('refs', refs)
                call cline_cavgassemble%set('nthr', nthr_here)
                call terminate_stream(params, 'SIMPLE_DISTR_CLUSTER2D HARD STOP 2')
                call xcavgassemble%execute(cline_cavgassemble)
                if( L_BENCH_GLOB ) rt_cavgassemble = toc(t_cavgassemble)
            endif
            ! objfun=euclid, part 4: sigma2 consolidation
            if( params%l_needs_sigma )then
                call cline_calc_sigma%set('which_iter', params%which_iter+1)
                call cline_calc_sigma%set('nthr',       nthr_here)
                call xcalc_group_sigmas%execute(cline_calc_sigma)
            endif
            ! print out particle parameters per iteration
            if( trim(params%print_corrs).eq.'yes' )then
                call build%spproj_field%write(string('ptcl2D_'//int2str_pad(params%which_iter,2)//'.txt'))
            endif
            ! check convergence
            call check_2Dconv(cline_check_2Dconv, build%spproj_field)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check_2Dconv%get_rarg('frac_srch')
            ! the below activates shifting
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check_2Dconv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check_2Dconv%get_rarg('trs'))
                    call job_descr%set('trs', str)
                endif
            endif
            l_converged = (iter >= params%minits) .and. (cline_check_2Dconv%get_carg('converged').eq.'yes')
            if( l_converged .or. iter==params%maxits ) then
                exit
            else
                call cline_check_2Dconv%delete('converged')
            endif
            if( L_BENCH_GLOB )then
                rt_tot  = toc(t_init)
                benchfname = 'CLUSTER2D_DISTR_BENCH_ITER'//int2str_pad(iter,3)//'.txt'
                call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', rt_scheduled
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', rt_merge_algndocs
                write(fnr,'(a,1x,f9.2)') 'cavgassemble    : ', rt_cavgassemble
                write(fnr,'(a,1x,f9.2)') 'total time      : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', (rt_init/rt_tot)           * 100.
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', (rt_scheduled/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', (rt_merge_algndocs/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'cavgassemble    : ', (rt_cavgassemble/rt_tot)   * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for : ',&
                    &((rt_init+rt_scheduled+rt_merge_algndocs+rt_cavgassemble)/rt_tot) * 100.
                call fclose(fnr)
            endif
            ! write jpeg image - used in gui to display image in logfile
            call mrc2jpeg_tiled(string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//params%ext%to_char(), string(CWD_GLOB)//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//JPG_EXT)
            write(logfhandle,'(A,A)')'>>> JPEG ', CWD_GLOB//'/'//CAVGS_ITER_FBODY//str_iter%to_char()//JPG_EXT
        end do
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        ! updates os_out
        if( trim(params%restore_cavgs).eq.'yes' )then
            if( file_exists(FRCS_FILE) ) call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
            call build%spproj%write_segment_inside('out', params%projfile)
            if( .not. params%l_polar )then
                finalcavgs = CAVGS_ITER_FBODY//int2str_pad(iter,3)//params%ext%to_char()
                call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
            endif
        endif
        call qsys_cleanup(params)
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', iter)
        ! end gracefully
        call build%spproj_field%kill
        call nice_communicator%terminate()
        call simple_touch(CLUSTER2D_FINISHED)
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_distr

    subroutine exec_cluster2D( self, cline )
        use simple_strategy2D_matcher, only: cluster2D_exec
        class(commander_cluster2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(commander_make_cavgs)        :: xmake_cavgs
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(commander_scale)             :: xscale
        type(commander_prob_tab2D_distr)  :: xprob_tab2D_distr
        type(simple_nice_communicator)    :: nice_communicator
        type(cmdline)             :: cline_make_cavgs, cline_scalerefs, cline_prob_tab2D
        type(parameters)          :: params
        type(builder),     target :: build
        type(starproject)         :: starproj
        type(string)              :: finalcavgs, refs_sc, fname, fname_sigma
        real,         allocatable :: corrs(:), corrs_all(:)
        integer,      allocatable :: order(:), class_cnt(:), class_all(:)
        integer :: startit, i, cnt, iptcl, ptclind
        integer :: iter_switch2euclid, j, io_stat, funit, class_ind, class_max
        logical :: converged, l_stream, l_switch2euclid, l_ml_reg, l_scale_inirefs
        call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('maxits') ) call cline%set('maxits', 30)
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        ! Streaming flag
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = cline%get_carg('stream')=='yes'
        endif
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        nice_communicator%stat_root%stage = "initialising"
        call nice_communicator%cycle()
        startit = 1
        if( cline%defined('startit') )startit = params%startit
        if( (startit == 1) .and. (.not.str_has_substr(params%refine,'prob')) )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        if( params%l_distr_exec )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
            call cluster2D_exec( params, build, cline, startit, converged )
            ! end gracefully
            call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
            call qsys_job_finished(params, string('simple_commanders_cluster2D :: exec_cluster2D'))
        else
            ! Initial references
            if( .not. cline%defined('refs') )then
                cline_make_cavgs = cline ! ncls is transferred here
                params%refs      = 'start2Drefs'//params%ext%to_char()
                params%refs_even = 'start2Drefs_even'//params%ext%to_char()
                params%refs_odd  = 'start2Drefs_odd'//params%ext%to_char()
                l_scale_inirefs  = .false.
                if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                    if( params%tseries .eq. 'yes' )then
                        if( cline%defined('nptcls_per_cls') )then
                            if( build%spproj%os_ptcl2D%any_state_zero() )then
                                THROW_HARD('cluster2D_nano does not allow state=0 particles, prune project before execution; exec_cluster2D_distr')
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
                            call xmake_cavgs%execute(cline_make_cavgs)
                            l_scale_inirefs  = .false.
                        else
                            if( trim(params%refine).eq.'inpl' )then
                                params%ncls = build%spproj%os_ptcl2D%get_n('class')
                                call cline%set('ncls', params%ncls)
                                call cline_make_cavgs%set('ncls', params%ncls)
                                call cline_make_cavgs%delete('tseries')
                                call cline_make_cavgs%set('refs', params%refs)
                                call xmake_cavgs%execute(cline_make_cavgs)
                                l_scale_inirefs  = .false.
                            else
                                call selection_from_tseries_imgfile(build%spproj, params%refs, params%box, params%ncls)
                                l_scale_inirefs  = .true.
                            endif
                        endif
                    else
                        select case(trim(params%cls_init))
                            case('ptcl')
                                ! initialization from raw images
                                call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                                l_scale_inirefs  = .true.
                            case('rand')
                                ! from noise
                                call noise_imgfile(params%refs, params%ncls, params%box_crop, params%smpd_crop)
                                l_scale_inirefs  = .false.
                            case('randcls')
                                if(.not.cline%defined('ncls')) THROW_HARD('NCLS must be provide with CLS_INIT=RANDCLS ')
                                ! initialization from random classes
                                do iptcl=1,params%nptcls
                                    if( build%spproj_field%get_state(iptcl) == 0 ) cycle
                                    call build%spproj_field%set(iptcl, 'class', irnd_uni(params%ncls))
                                    call build%spproj_field%set(iptcl, 'w',     1.0)
                                    call build%spproj_field%e3set(iptcl,ran3()*360.0)
                                end do
                                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                                call cline_make_cavgs%set('refs', params%refs)
                                call xmake_cavgs%execute(cline_make_cavgs)
                                l_scale_inirefs  = .false.
                            case DEFAULT
                                THROW_HARD('Unsupported mode of initial class generation CLS_INIT='//trim(params%cls_init))
                        end select
                    endif
                    ! scale references to box_crop
                    if( l_scale_inirefs )then
                        refs_sc = 'refs'//SCALE_SUFFIX//params%ext%to_char()
                        call cline_scalerefs%set('stk',    params%refs)
                        call cline_scalerefs%set('outstk', refs_sc)
                        call cline_scalerefs%set('smpd',   params%smpd)
                        call cline_scalerefs%set('newbox', params%box_crop)
                        call xscale%execute(cline_scalerefs)
                        call simple_rename(refs_sc, params%refs)
                    endif
                    call copy_imgfile(params%refs, params%refs_even, params%smpd_crop, [1,params%ncls])
                    call copy_imgfile(params%refs, params%refs_odd,  params%smpd_crop, [1,params%ncls])
                else
                    call cline_make_cavgs%set('refs', params%refs)
                    call xmake_cavgs%execute(cline_make_cavgs)
                endif
                call cline%set('refs', params%refs)
            endif
            params%startit = startit
            params%outfile = 'algndoc'//METADATA_EXT
            ! variable neighbourhood size
            if( cline%defined('extr_iter') )then
                ! all good
            else
                params%extr_iter = params%startit
            endif
            ! objective functions
            select case( params%cc_objfun )
            case( OBJFUN_CC )
                params%l_needs_sigma = .false.
                params%needs_sigma   = 'no'
                call cline%set('ml_reg','no')
                params%ml_reg        = 'no'
                params%l_ml_reg      = .false.
            case( OBJFUN_EUCLID )
                fname_sigma = sigma2_star_from_iter(params%startit)
                if( .not.file_exists(fname_sigma) )then
                    THROW_HARD('Sigma2 file does not exists: '//fname_sigma%to_char())
                endif
            end select
            ! initialise progress monitor
            if(.not. l_stream) call progressfile_init()
            ! Main loop
            do i = 1, params%maxits
                params%which_iter = params%startit
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
                write(logfhandle,'(A)')   '>>>'
                nice_communicator%stat_root%stage = "iteration " // int2str(params%which_iter)
                call nice_communicator%cycle()
                call cline%set('which_iter', params%which_iter)
                ! refine=prob
                if( str_has_substr(params%refine, 'prob') )then
                    cline_prob_tab2D = cline
                    call cline_prob_tab2D%set('prg',        'prob_tab2D')
                    call cline_prob_tab2D%set('which_iter', params%which_iter)
                    call cline_prob_tab2D%set('startit',    params%which_iter)
                    call cline_prob_tab2D%set('refs',       params%refs)
                    call cline_prob_tab2D%set('frcs',       FRCS_FILE)
                    call cline_prob_tab2D%set('extr_iter',  params%extr_iter)
                    call xprob_tab2D_distr%execute( cline_prob_tab2D )
                    ! read back sampling
                    call build%spproj%read_segment(params%oritype, params%projfile)
                endif
                ! stochastic search
                call cluster2D_exec( params, build, cline, params%startit, converged )
                ! objective functions
                if( params%l_needs_sigma )then
                    params%which_iter = params%which_iter + 1
                    call cline%set('which_iter', params%which_iter)
                    call xcalc_group_sigmas%execute(cline)
                    params%which_iter = params%which_iter - 1
                    call cline%set('which_iter', params%which_iter)
                endif
                ! cooling of the randomization rate
                params%extr_iter = params%extr_iter + 1
                ! update project with the new orientations (this needs to be here rather than after convergence for the current implementation to work)
                call build%spproj%write_segment_inside(params%oritype)
                ! write cavgs starfile for iteration
                call starproj%export_cls2D(build%spproj, params%which_iter)
                ! exit condition
                if( converged )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', params%startit)
                    call del_file(params%outfile)
                    ! update os_out
                    if( trim(params%restore_cavgs).eq.'yes' )then
                        if( file_exists(FRCS_FILE) ) call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
                        if( .not. params%l_polar )then ! no Cartesian class averages in the polar version
                            finalcavgs = CAVGS_ITER_FBODY//int2str_pad(params%startit,3)//params%ext%to_char()
                            call build%spproj%add_cavgs2os_out(finalcavgs, build%spproj%get_smpd(), imgkind='cavg')
                            if( trim(params%chunk).eq.'yes' )then
                                call cavger_write_eo(params%refs_even, params%refs_odd)
                                call cavger_readwrite_partial_sums('write')
                                call cavger_kill
                            else
                               if( trim(params%tseries).eq.'yes' )then
                                    call cavger_write_eo(params%refs_even, params%refs_odd)
                                    call cavger_kill
                                else
                                    call cavger_kill(dealloccavgs=.false.)
                                endif
                            endif
                        endif
                        call build%spproj%write_segment_inside('out', params%projfile)
                    endif
                    exit
                endif
                ! update iteration counter
                params%startit = params%startit + 1
            end do
            nice_communicator%stat_root%stage = "terminating"
            call nice_communicator%cycle()
            ! print CSV file of correlation vs particle number
            if( trim(params%print_corrs).eq.'yes' )then
                class_all = build%spproj%os_ptcl2D%get_all_asint('class')
                class_max = maxval(class_all)
                ! counting the number of particles in each class
                allocate(class_cnt(class_max))
                do class_ind = 1, class_max
                    class_cnt(class_ind) = 0
                    do j = 1, size(class_all)
                        if( class_all(j) == class_ind ) class_cnt(class_ind) = class_cnt(class_ind) + 1
                    enddo
                enddo
                ! print all sorted corrs
                corrs_all = build%spproj%os_ptcl2D%get_all('corr')
                order     = (/(j,j=1,size(corrs_all))/)
                call hpsort(corrs_all, order)
                fname = 'ptcls_vs_cavgs_corrs_iter'// int2str(params%which_iter) //'.csv'
                call fopen(funit, fname, 'replace', 'unknown', iostat=io_stat, form='formatted')
                call fileiochk('cluster2D fopen failed: '//fname%to_char(), io_stat)
                write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
                do j = 1,size(corrs_all)
                    write(funit,*) int2str(order(j))//CSV_DELIM//real2str(corrs_all(j))
                end do
                call fclose(funit)
                ! print sorted corrs for each class
                do class_ind = 1, class_max
                    if( allocated(corrs) ) deallocate(corrs)
                    allocate(corrs(class_cnt(class_ind)), source=0.)
                    cnt = 0
                    do j = 1, size(corrs_all)
                        if( class_all(j) == class_ind )then
                            cnt = cnt + 1
                            corrs(cnt) = corrs_all(j)
                        endif
                    enddo
                    if( allocated(order) ) deallocate(order)
                    order = (/(j,j=1,class_cnt(class_ind))/)
                    call hpsort(corrs, order)
                    fname = 'ptcls_vs_cavgs_corrs_cls_'// int2str(class_ind) //'.csv'
                    call fopen(funit, fname, 'replace', 'unknown', iostat=io_stat, form='formatted')
                    call fileiochk('cluster2D fopen failed: '//fname%to_char(), io_stat)
                    write(funit,*) 'PTCL_INDEX'//CSV_DELIM//'CORR'
                    do j = 1,size(corrs)
                        write(funit,*) int2str(order(j))//CSV_DELIM//real2str(corrs(j))
                    end do
                    call fclose(funit)
                enddo
            endif
            ! end gracefully
            call nice_communicator%terminate()
            call simple_touch(CLUSTER2D_FINISHED)
            call simple_end('**** SIMPLE_CLUSTER2D NORMAL STOP ****')
        endif
    end subroutine exec_cluster2D

    subroutine exec_prob_tab2D_distr( self, cline )
        use simple_eul_prob_tab2D,     only: eul_prob_tab2D
        class(commander_prob_tab2D_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(builder)              :: build
        type(parameters)           :: params
        type(commander_prob_tab2D) :: xprob_tab2D
        type(eul_prob_tab2D)       :: eulprob
        type(cmdline)              :: cline_prob_tab2D
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        integer :: nptcls
        logical :: l_maxpop, l_stream
        l_stream = .false.
        if(cline%defined('stream')) l_stream = cline%get_carg('stream').eq.'yes'
        call cline%set('mkdir',   'no')
        call cline%set('stream',  'no')
        call cline%set('oritype', 'ptcl2D')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( l_stream ) call cline%set('stream', 'yes')
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! Whether to weight based-on the top maxpop particles
        l_maxpop = cline%defined('maxpop') .and. (params%maxpop > 0)
        ! sample particles
        if( params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([1,params%nptcls], params%update_frac, nptcls, pinds, .true.)
        else
            call build%spproj_field%sample4update_all([1,params%nptcls], nptcls, pinds, .true.)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        ! more prep
        call eulprob%new(params, build, pinds)
        ! generating all scores
        cline_prob_tab2D = cline
        call cline_prob_tab2D%set('prg', 'prob_tab2D' )
        ! execution
        if( .not.cline_prob_tab2D%defined('nparts') )then
            ! shared memory
            call xprob_tab2D%execute(cline_prob_tab2D)
        else
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab2D%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! reading scores from all parts
        call eulprob%read_table_parts_to_glob
        ! perform assignment
        if( l_stream )then
            select case(trim(params%refine))
                case('prob_smpl')
                    call eulprob%assign_smpl(build%spproj_field, l_maxpop)
                case DEFAULT
                    THROW_HARD('Unsupported REFINE flag: '//trim(params%refine))
            end select
        else
            if( params%which_iter == 1 )then
                ! Always greedy assignement with first iteration
                call eulprob%assign_greedy(l_maxpop)
            else
                select case(trim(params%refine))
                case('prob')
                    if( params%extr_iter <= params%extr_lim )then
                        call eulprob%normalize_table
                        call eulprob%assign_prob(build%spproj_field, l_maxpop)
                    else
                        call eulprob%assign_smpl(build%spproj_field, l_maxpop)
                    endif
                case('prob_smpl')
                    call eulprob%assign_smpl(build%spproj_field, l_maxpop)
                case('prob_greedy')
                    call eulprob%assign_greedy(l_maxpop)
                case DEFAULT
                    THROW_HARD('Unsupported REFINE flag: '//trim(params%refine))
                end select
            endif
        endif
        ! write
        call eulprob%write_assignment(string(ASSIGNMENT_FBODY)//'.dat')
        ! cleanup
        call eulprob%kill
        call cline_prob_tab2D%kill
        call qenv%kill
        call job_descr%kill
        call qsys_job_finished(params, string('simple_commanders_cluster2D :: exec_prob_tab2D_distr'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_TAB2D_DISTR NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D_distr

    subroutine exec_prob_tab2D( self, cline )
        use simple_strategy2D_matcher
        use simple_strategy2D3D_common, only: set_bp_range2D
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        type(string)                  :: fname
        type(polarft_calc)            :: pftc
        type(builder)                 :: build
        type(parameters)              :: params
        type(eul_prob_tab2D)          :: eulprob
        real    :: frac_srch_space
        integer :: nptcls
        logical :: l_ctf, l_stream, l_alloc_read_cavgs
        l_stream = .false.
        if(cline%defined('stream')) l_stream = cline%get_carg('stream') .eq.'yes'
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
        call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
        ! Nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_tab2D_distr, above)
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab2D requires prior particle sampling')
        endif
        ! Resolution range
        frac_srch_space = build%spproj_field%get_avg('frac')
        if( file_exists(params%frcs) ) call build%clsfrcs%read(params%frcs)
        call set_bp_range2D( params, build, cline, params%which_iter, frac_srch_space )
        ! Read references
        l_alloc_read_cavgs = .true.
        if( .not.l_distr_exec_glob )then
            l_alloc_read_cavgs = params%which_iter==1
        endif
        if( l_alloc_read_cavgs )then
            if( .not. cline%defined('refs') )then
                THROW_HARD('need refs to be part of command line for cluster2D execution')
            endif
            call cavger_new(params, build, pinds, alloccavgs=.true.)
            call cavger_read_all
        else
            call cavger_new(params, build, pinds, alloccavgs=.false.)
        endif
        ! init scorer & prep references
        call preppftc4align2D(nptcls, params%which_iter, l_stream)
        ! minor cleanup
        call cavger_kill(dealloccavgs=l_distr_exec_glob)
        ! prep particles
        l_ctf = build%spproj%get_ctfflag('ptcl2D',iptcl=params%fromp).ne.'no'
        call prep_batch_particles2D(nptcls)
        call build_batch_particles2D(nptcls, pinds)
        ! init prob table
        call eulprob%new(params, build, pinds)
        fname = DIST_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
        ! Fill probability table
        if( l_stream )then
            select case(trim(params%refine))
                case('prob_smpl')
                    call eulprob%fill_table_smpl_stream(build%spproj_field)
                case DEFAULT
                    THROW_HARD('Unsupported REFINE flag: '//trim(params%refine))
            end select
        else
            if( params%which_iter == 1 )then
                ! always greedy in-plane in first iteration
                call eulprob%fill_table_greedy
            else
                select case(trim(params%refine))
                    case('prob','prob_smpl')
                        call eulprob%fill_table_smpl
                    case('prob_greedy')
                        call eulprob%fill_table_greedy
                    case DEFAULT
                        THROW_HARD('Unsupported REFINE flag: '//trim(params%refine))
                end select
            endif
        endif
        call pftc%kill
        ! call build%esig%kill
        call clean_batch_particles2D
        ! write
        call eulprob%write_table(fname)
        ! clean & end
        call eulprob%kill
        call build%kill_general_tbox
        call build%kill_strategy2D_tbox
        call qsys_job_finished(params, string('simple_commanders_cluster2D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_ppca_denoise_classes( self, cline )
        use simple_imgproc,       only: make_pcavecs
        use simple_pca,           only: pca
        use simple_pca_svd,       only: pca_svd
        use simple_kpca_svd,      only: kpca_svd
        use simple_ppca_inmem,    only: ppca_inmem
        class(commander_ppca_denoise_classes), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        integer,          parameter   :: MAXPCAITS = 15
        class(pca),       pointer     :: pca_ptr  => null()
        type(parameters)              :: params
        type(builder)                 :: build
        type(image),      allocatable :: imgs(:), imgs_ori(:)
        type(image)                   :: cavg, img, timg
        type(oris)                    :: os
        type(sp_project), target      :: spproj
        type(string)                  :: label, fname, fname_denoised, fname_cavgs, fname_cavgs_denoised
        type(string)                  :: fname_oris, fname_denoised_ori, fname_ori, fname_class_ptcls_den
        integer,          allocatable :: cls_inds(:), pinds(:), cls_pops(:), ori_map(:)
        real,             allocatable :: avg(:), avg_pix(:), pcavecs(:,:), tmpvec(:)
        real    :: shift(2), loc(2), dist(2), e3, kw, mat(2,2), mat_inv(2,2)
        complex :: fcompl, fcompll
        integer :: npix, i, j, ncls, nptcls, cnt1, cnt2, neigs, h, k, win_corner(2),&
                  &l, ll, m, mm, phys(2), logi_lims(3,2), cyc_lims(3,2), cyc_limsR(2,2), errflg
        logical :: l_phflip, l_transp_pca, l_pre_norm ! pixel-wise learning
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        if( .not. cline%defined('neigs')   ) call cline%set('neigs',    10)
        call build%init_params_and_build_general_tbox(cline, params, do3d=(trim(params%oritype) .eq. 'ptcl3D'))
        call spproj%read(params%projfile)
        select case(trim(params%oritype))
            case('ptcl2D')
                label = 'class'
            case('ptcl3D')
                label = 'proj'
                call build%spproj_field%proj2class
            case DEFAULT
                THROW_HARD('ORITYPE not supported!')
        end select
        l_transp_pca = (trim(params%transp_pca) .eq. 'yes')
        l_pre_norm   = (trim(params%pre_norm)   .eq. 'yes')
        l_phflip     = .false.
        select case( spproj%get_ctfflag_type(params%oritype) )
            case(CTFFLAG_NO)
                THROW_WARN('No CTF information could be found, phase flipping is deactivated')
            case(CTFFLAG_FLIP)
                THROW_WARN('Images have already been phase-flipped, phase flipping is deactivated')
            case(CTFFLAG_YES)
                l_phflip = .true.
            case DEFAULT
                THROW_HARD('UNSUPPORTED CTF FLAG')
        end select
        cls_inds = build%spproj_field%get_label_inds(label%to_char())
        if( cline%defined('class') .and. cline%defined('ncls') )then
            THROW_HARD('EITHER class OR ncls CAN BE DEFINED')
        endif
        if( cline%defined('class') )then
            cls_inds = pack(cls_inds, mask=(cls_inds == params%class))
        endif
        ncls = size(cls_inds)
        if( cline%defined('ncls') )then
            ncls     = params%ncls
            cls_inds = cls_inds(1:ncls)
        endif
        allocate(cls_pops(ncls), source=0)
        do i = 1, ncls
            call build%spproj_field%get_pinds(cls_inds(i), label%to_char(), pinds)
            if( allocated(pinds) )then
                cls_pops(i) = size(pinds)
                nptcls = nptcls + cls_pops(i)
                deallocate(pinds)
            endif
        end do
        cls_inds = pack(cls_inds, mask=cls_pops > 2)
        nptcls   = sum(cls_pops,  mask=cls_pops > 2)
        ncls     = size(cls_inds)
        if( trim(params%pca_ori_stk) .eq. 'yes' ) allocate(ori_map(nptcls))
        call os%new(nptcls, is_ptcl=.true.)
        fname                = 'ptcls.mrcs'
        fname_denoised       = 'ptcls_denoised.mrcs'
        fname_cavgs          = 'cavgs.mrcs'
        fname_cavgs_denoised = 'cavgs_denoised.mrcs'
        fname_oris           = 'oris_denoised.txt'
        fname_ori            = 'ptcls_ori_order.mrcs'
        fname_denoised_ori   = 'ptcls_denoised_ori_order.mrcs'
        cnt1 = 0
        cnt2 = 0
        ! pca allocation
        select case(trim(params%pca_mode))
            case('ppca')
                allocate(ppca_inmem :: pca_ptr)
            case('pca_svd')
                allocate(pca_svd    :: pca_ptr)
            case('kpca')
                allocate(kpca_svd   :: pca_ptr)
        end select
        do i = 1, ncls
            call progress_gfortran(i,ncls)
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
                do j = 1, size(imgs)
                    call imgs(j)%copy_fast(imgs_ori(j))
                enddo
            else
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg)
            endif
            nptcls = size(imgs)
            if( trim(params%neigs_per).eq.'yes' )then
                if( params%neigs >= 99 )then
                    THROW_WARN('neigs is greater than 99% the number of particles within this class. All eigens are used now!')
                    neigs = nptcls - 1
                else
                    neigs = max(2, nint(real(params%neigs * nptcls) / 100.))
                endif
            else
                neigs = params%neigs
                if( neigs >= nptcls )then
                    THROW_WARN('neigs is greater than the number of particles within this class. All eigens are used now!')
                    neigs = nptcls - 1
                endif
            endif
            if( l_pre_norm )then
                do j = 1, nptcls
                    call imgs(j)%norm
                end do
            endif
            do j = 1, nptcls
                cnt1 = cnt1 + 1
                call imgs(j)%write(fname, cnt1)
            end do
            call cavg%write(fname_cavgs, i)
            ! performs ppca
            if( trim(params%projstats).eq.'yes' )then
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca, avg_pix=avg_pix)
            else
                call make_pcavecs(imgs, npix, avg, pcavecs, transp=l_transp_pca)
            endif
            if( allocated(tmpvec) ) deallocate(tmpvec)
            if( l_transp_pca )then
                call pca_ptr%new(npix, nptcls, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(nptcls))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, npix
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
                pcavecs = transpose(pcavecs)
            else
                call pca_ptr%new(nptcls, npix, neigs)
                call pca_ptr%master(pcavecs, MAXPCAITS)
                allocate(tmpvec(npix))
                !$omp parallel do private(j,tmpvec) default(shared) proc_bind(close) schedule(static)
                do j = 1, nptcls
                    call pca_ptr%generate(j, avg, tmpvec)
                    pcavecs(:,j) = tmpvec
                end do
                !$omp end parallel do
            endif
            if( trim(params%projstats).eq.'yes' )then
                call cavg%unserialize(avg_pix)
                call cavg%write(string('cavgs_unserialized.mrcs'), i)
            endif
            ! output
            call cavg%zero_and_unflag_ft
            if( trim(params%pca_img_ori) .eq. 'yes' )then
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs_ori(j)%unserialize(pcavecs(:,j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs_ori(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                end do
                call transform_ptcls(params, build, spproj, params%oritype, cls_inds(i), imgs, pinds, phflip=l_phflip, cavg=cavg, imgs_ori=imgs_ori)
            else
                fname_class_ptcls_den = 'class'//int2str_pad(i,4)//'ptcls.mrcs'
                do j = 1, nptcls
                    cnt2 = cnt2 + 1
                    call imgs(j)%unserialize(pcavecs(:,j))
                    call cavg%add(imgs(j))
                    call os%transfer_ori(cnt2, build%spproj_field, pinds(j))
                    call imgs(j)%write(fname_class_ptcls_den, j)
                    ! call imgs(j)%write(fname_denoised, cnt2)
                    if( trim(params%pca_ori_stk) .eq. 'yes' ) ori_map(pinds(j)) = cnt2
                    call imgs(j)%kill
                end do
                call cavg%div(real(nptcls))
            endif
            call cavg%write(fname_cavgs_denoised, i)
        end do
        if( trim(params%pca_ori_stk) .eq. 'yes' )then
            call  img%new([params%boxpd,params%boxpd,1],params%smpd, wthreads=.false.)
            call timg%new([params%boxpd,params%boxpd,1],params%smpd, wthreads=.false.)
            logi_lims      = img%loop_lims(2)
            cyc_lims       = img%loop_lims(3)
            cyc_limsR(:,1) = cyc_lims(1,:)
            cyc_limsR(:,2) = cyc_lims(2,:)
            do i = 1, cnt2
                shift = build%spproj_field%get_2Dshift(i)
                e3    = build%spproj_field%e3get(i)
                do j = 1, 2
                    call  img%zero_and_flag_ft
                    call timg%zero_and_flag_ft
                    call cavg%zero_and_unflag_ft
                    if( j == 1 )then
                        call cavg%read(fname_denoised, ori_map(i))
                    else
                        call cavg%read(fname,          ori_map(i))
                    endif
                    call cavg%pad_fft(img)
                    ! particle rotation
                    call rotmat2d(-e3, mat)
                    call matinv(mat, mat_inv, 2, errflg)
                    !$omp parallel do collapse(2) private(h,k,loc,win_corner,dist,l,ll,m,mm,phys,kw,fcompl,fcompll) default(shared) proc_bind(close) schedule(static)
                    do h = logi_lims(1,1),logi_lims(1,2)
                        do k = logi_lims(2,1),logi_lims(2,2)
                            ! Rotation
                            loc        = matmul(real([h,k]),mat_inv)
                            win_corner = floor(loc) ! bottom left corner
                            dist       = loc - real(win_corner)
                            ! Bi-linear interpolation
                            l      = cyci_1d(cyc_limsR(:,1), win_corner(1))
                            ll     = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                            m      = cyci_1d(cyc_limsR(:,2), win_corner(2))
                            mm     = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                            ! l, bottom left corner
                            phys   = img%comp_addr_phys(l,m)
                            kw     = (1.-dist(1))*(1.-dist(2))   ! interpolation kernel weight
                            fcompl = kw * img%get_cmat_at(phys(1), phys(2),1)
                            ! l, bottom right corner
                            phys   = img%comp_addr_phys(l,mm)
                            kw     = (1.-dist(1))*dist(2)
                            fcompl = fcompl + kw * img%get_cmat_at(phys(1), phys(2),1)
                            if( l < 0 ) fcompl = conjg(fcompl) ! conjugation when required!
                            ! ll, upper left corner
                            phys    = img%comp_addr_phys(ll,m)
                            kw      = dist(1)*(1.-dist(2))
                            fcompll = kw * img%get_cmat_at(phys(1), phys(2),1)
                            ! ll, upper right corner
                            phys    = img%comp_addr_phys(ll,mm)
                            kw      = dist(1)*dist(2)
                            fcompll = fcompll + kw * img%get_cmat_at(phys(1), phys(2),1)
                            if( ll < 0 ) fcompll = conjg(fcompll) ! conjugation when required!
                            ! update with interpolated values
                            phys = img%comp_addr_phys(h,k)
                            call timg%set_cmat_at(phys(1),phys(2),1, fcompl + fcompll)
                        end do
                    end do
                    !$omp end parallel do
                    ! shift
                    call timg%shift2Dserial(shift)
                    call timg%ifft
                    call timg%clip(cavg)
                    if( j == 1 )then
                        call cavg%write(fname_denoised_ori, i)
                    else
                        call cavg%write(fname_ori, i)
                    endif
                enddo
            enddo
        endif
        call os%zero_inpl
        call os%write(fname_oris)
        if( trim(params%projstats).eq.'yes' ) call build%spproj_field%write(string('ptcl_field.txt'))
        ! cleanup
        deallocate(imgs)
        call build%kill_general_tbox
        call os%kill
        ! end gracefully
        call simple_end('**** SIMPLE_PPCA_DENOISE_CLASSES NORMAL STOP ****')
    end subroutine exec_ppca_denoise_classes

    ! UTILITIES

    subroutine check_2Dconv( cline, os )
        use simple_convergence, only: convergence
        class(cmdline), intent(inout) :: cline
        class(oris),    intent(inout) :: os
        type(parameters)  :: params
        type(convergence) :: conv
        logical :: converged, l_stream
        call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        l_stream = .false.
        if( cline%defined('stream') )then
            l_stream = cline%get_carg('stream') == 'yes'
        endif
        ! convergence check
        converged = conv%check_conv2D(params, cline, os, os%get_n('class'), params%msk)
        ! Update progress file
        if(.not. l_stream) call progressfile_update(conv%get('progress'))
        call cline%set('frac_srch', conv%get('frac_srch'))
        ! activates shift search
        if( params%l_doshift ) call cline%set('trs', params%trs)
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_CHECK_2DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine check_2Dconv

end module simple_commanders_cluster2D
