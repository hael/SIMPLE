! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_builder,          only: builder
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_parameters,       only: parameters
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_qsys_env,         only: qsys_env
use simple_cluster_seed,     only: gen_labelling
use simple_commander_volops, only: postprocess_commander
use simple_commander_euclid, only: calc_pspec_commander_distr
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_commander_distr
public :: refine3D_commander
public :: check_3Dconv_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: nspace_commander
 contains
   procedure :: execute      => exec_nspace
end type nspace_commander

type, extends(commander_base) :: refine3D_commander_distr
  contains
    procedure :: execute      => exec_refine3D_distr
end type refine3D_commander_distr

type, extends(commander_base) :: refine3D_commander
  contains
    procedure :: execute      => exec_refine3D
end type refine3D_commander

type, extends(commander_base) :: check_3Dconv_commander
  contains
    procedure :: execute      => exec_check_3Dconv
end type check_3Dconv_commander

contains

    subroutine exec_nspace(self,cline)
        class(nspace_commander), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: o
        real             :: ares
        integer          :: i
        call params%new(cline)
        do i=500,5000,500
            o = oris(i, is_ptcl=.false.)
            call o%spiral
            ares = o%find_angres()
            write(logfhandle,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, params%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_refine3D_distr( self, cline )
        use simple_commander_rec, only: reconstruct3D_commander_distr
        use simple_fsc,           only: plot_fsc
        class(refine3D_commander_distr), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        type(check_3Dconv_commander)        :: xcheck_3Dconv
        type(postprocess_commander)         :: xpostprocess
        type(refine3D_commander)            :: xrefine3D_shmem
        ! command lines
        type(cmdline)    :: cline_reconstruct3D_distr
        type(cmdline)    :: cline_calc_pspec_distr
        type(cmdline)    :: cline_calc_group_sigmas
        type(cmdline)    :: cline_check_3Dconv
        type(cmdline)    :: cline_volassemble
        type(cmdline)    :: cline_postprocess
        integer(timer_int_kind) :: t_init,   t_scheduled,  t_merge_algndocs,  t_volassemble,  t_tot
        real(timer_int_kind)    :: rt_init, rt_scheduled, rt_merge_algndocs, rt_volassemble, rt_tot
        character(len=STDLEN)   :: benchfname
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv, qenv_cls3D
        type(chash)      :: job_descr, job_descr_cavgs
        character(len=:),          allocatable :: vol_fname, prev_refine_path, target_name
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        integer,                   allocatable :: state_pops(:), tmp_iarr(:)
        real,                      allocatable :: res(:), tmp_rarr(:), fsc(:)
        character(len=STDLEN)     :: vol, vol_iter, str, str_iter, fsc_templ
        character(len=STDLEN)     :: vol_even, vol_odd, str_state, fsc_file, volpproc, vollp
        character(len=LONGSTRLEN) :: volassemble_output
        logical :: err, vol_defined, have_oris, do_abinitio, converged, fall_over
        logical :: l_projmatch, l_switch2eo, l_switch2euclid, l_continue, l_multistates
        logical :: l_ptclw, l_combine_eo, l_lpset, l_griddingset
        real    :: corr, corr_prev, smpd, lplim
        integer :: ldim(3), i, state, iter, box, nfiles, niters, iter_switch2euclid, ifoo
        integer :: ncls, icls, ind, fnr
        if( .not. cline%defined('nparts') )then
            call xrefine3D_shmem%execute(cline)
            return
        endif
        l_multistates = cline%defined('nstates')
        l_lpset       = cline%defined('lp')
        l_griddingset = cline%defined('gridding')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('ptclw')   ) call cline%set('ptclw',       'no')
        if( .not. cline%defined('lp_iters')) call cline%set('lp_iters',      1.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! objfun=euclid logics, part 1
        l_switch2euclid  = .false.
        if( cline%defined('objfun') )then
            l_continue = .false.
            if( cline%defined('continue') ) l_continue = trim(cline%get_carg('continue')).eq.'yes'
            if( (trim(cline%get_carg('objfun')).eq.'euclid') .and. .not.l_continue )then
                l_switch2euclid = .true.
                call cline%set('objfun','cc')
                l_ptclw = trim(cline%get_carg('ptclw')).eq.'yes'
                call cline%set('ptclw', 'no')
            endif
        endif
        ! init
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
            case DEFAULT
                write(logfhandle,*)'Unsupported ORITYPE; simple_commander_refine3D :: exec_refine3D_distr'
        end select
        if( fall_over ) THROW_HARD('no particles found! exec_refine3D_distr')
        if( .not. l_multistates .and. params%nstates >  1 ) THROW_HARD('nstates > 1 but refine mode is single')
        ! take care of automask flag
        if( cline%defined('automsk') ) call cline%delete('automsk')
        if( params%l_automsk .and. l_multistates ) THROW_HARD('automsk.ne.no not currenty supported for multi-state refinement')
        ! switch from low-pass to e/o refinement
        l_switch2eo = params%lp_iters >= 1
        if( .not.l_lpset ) l_switch2eo = .false. ! is already e/o
        ! final iteration with combined e/o
        l_combine_eo = .false.
        if( cline%defined('combine_eo') )then
            l_combine_eo = .true.
            call cline%set('combine_eo','no')
            params%combine_eo = 'no'
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' ) call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr = cline
        cline_calc_pspec_distr    = cline
        cline_check_3Dconv        = cline
        cline_volassemble         = cline
        cline_postprocess         = cline
        cline_calc_group_sigmas   = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' )    ! required for distributed call
        call cline_calc_pspec_distr%set(    'prg', 'calc_pspec' )       ! required for distributed call
        call cline_postprocess%set(         'prg', 'postprocess' )      ! required for local call
        call cline_calc_group_sigmas%set(   'prg', 'calc_group_sigmas' )! required for local call
        if( trim(params%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set('pgrp', 'c1')
        call cline_postprocess%set('mirr',    'no')
        call cline_postprocess%set('mkdir',   'no')
        call cline_postprocess%set('imgkind', 'vol')
        if( trim(params%oritype).eq.'cls3D' ) call cline_postprocess%set('imgkind', 'vol_cavg')
        call cline_postprocess%delete('l_nonuniform') ! to save time in the iterations
        ! for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates))
        ! removes unnecessary volume keys and generates volassemble finished names
        do state = 1,params%nstates
            vol = 'vol'//int2str( state )
            call cline_check_3Dconv%delete( vol )
            call cline_volassemble%delete( vol )
            call cline_postprocess%delete( vol )
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        if( l_switch2euclid ) call cline_volassemble%set('needs_sigma','yes')
        ! E/O PARTITIONING
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! STATE LABEL INIT
        if( l_multistates )then
            if( build%spproj_field%get_n('state') /= params%nstates )then
                call gen_labelling(build%spproj_field, params%nstates, 'squared_uniform')
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! GENERATE INITIAL NOISE POWER ESTIMATES
        if( l_switch2euclid .and. params%continue.ne.'yes' )then
            call build%spproj_field%set_all2single('w', 1.0)
            call build%spproj%write_segment_inside(params%oritype)
            call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
        endif
        ! GENERATE STARTING MODELS & ORIENTATIONS
        if( params%continue .eq. 'yes' )then
            ! we are continuing from a previous refinement round,
            ! i.e. projfile is fetched from a X_refine3D dir
            ! set starting volume(s), iteration number & previous refinement path
            do state=1,params%nstates
                ! volume(s)
                vol = 'vol' // int2str(state)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%get_vol('vol_cavg', state, vol_fname, smpd, box)
                else
                    call build%spproj%get_vol('vol', state, vol_fname, smpd, box)
                endif
                call cline%set(trim(vol), vol_fname)
                params%vols(state) = vol_fname
            end do
            prev_refine_path = get_fpath(vol_fname)
            ! carry over FSCs
            ! one FSC file per state
            do state=1,params%nstates
                str_state = int2str_pad(state,2)
                fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                call simple_copy_file(trim(prev_refine_path)//trim(fsc_file), fsc_file)
            end do
            ! one FRC file for all states
            if( file_exists(trim(prev_refine_path)//trim(FRCS_FILE)) )then
                call simple_copy_file(trim(prev_refine_path)//trim(FRCS_FILE), trim(FRCS_FILE))
            endif
            ! if we are doing fractional volume update, partial reconstructions need to be carried over
            if( params%l_frac_update )then
                call simple_list_files(prev_refine_path//'*recvol_state*part*', list)
                nfiles = size(list)
                err = params%nparts * 4 /= nfiles
                if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
                deallocate(list)
            endif
            ! if we are doing objfun=euclid the sigm estimates need to be carried over
            if( trim(params%objfun) .eq. 'euclid' )then
                call cline%set('needs_sigma','yes')
                call cline_reconstruct3D_distr%set('needs_sigma','yes')
                call cline_volassemble%set('needs_sigma','yes')
                if( .not.l_griddingset ) call cline%set('gridding','yes')
                call simple_list_files(prev_refine_path//trim(SIGMA2_FBODY)//'*', list)
                nfiles = size(list)
                if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
                deallocate(list)
            endif
        endif
        vol_defined = .false.
        do state = 1,params%nstates
            vol = 'vol' // int2str(state)
            if( cline%defined(trim(vol)) )then
                vol_defined = .true.
                call find_ldim_nptcls(trim(params%vols(state)),ldim,ifoo)
                if( ldim(1) /= params%box )then
                    THROW_HARD('Incompatible dimensions between input volume and images: '//params%vols(state))
                endif
            endif
        enddo
        have_oris   = .not. build%spproj%is_virgin_field(params%oritype)
        do_abinitio = .not. have_oris .and. .not. vol_defined
        if( do_abinitio )then
            call build%spproj_field%rnd_oris
            call build%spproj_field%zero_shifts
            have_oris = .true.
            call build%spproj%write_segment_inside(params%oritype)
        endif
        l_projmatch = .false.
        if( have_oris .and. .not. vol_defined )then
            ! reconstructions needed
            call xreconstruct3D_distr%execute( cline_reconstruct3D_distr )
            do state = 1,params%nstates
                ! rename volumes and update cline
                str_state = int2str_pad(state,2)
                vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                str       = trim(STARTVOL_FBODY)//trim(str_state)//params%ext
                call      simple_rename( trim(vol), trim(str) )
                vol       = 'vol'//trim(int2str(state))
                call      cline%set( trim(vol), trim(str) )
                vol_even  = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                str       = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params%ext
                call      simple_rename( trim(vol_even), trim(str) )
                vol_odd   = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                str       = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params%ext
                call      simple_rename( trim(vol_odd), trim(str) )
            enddo
        else if( vol_defined .and. params%continue .ne. 'yes' )then
            ! projection matching
            l_projmatch = .true.
            if( .not. have_oris )then
                if( str_has_substr(params%refine, 'neigh')) then
                    THROW_HARD('neigh refinement mode requires input orientations')
                endif
            endif
            if( .not.l_lpset )then
                THROW_HARD('LP needs be defined for the first step of projection matching!')
                call cline%delete('update_frac')
            endif
            if( (.not.str_has_substr(params%refine, 'neigh')) .and. (trim(params%refine).ne.'cont') )then
                ! this forces the first round of alignment on the starting model(s)
                ! to be greedy and the subseqent ones to be whatever the refinement flag is set to
                call build%spproj%os_ptcl3D%delete_3Dalignment(keepshifts=.true.)
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! EXTREMAL DYNAMICS
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! objfun=euclid logics, part 2
        iter_switch2euclid = -1
        if( l_switch2euclid )then
            iter_switch2euclid = 1
            if( cline%defined('update_frac') ) iter_switch2euclid = ceiling(1./(params%update_frac+0.001))
            if( l_projmatch .and. l_switch2eo ) iter_switch2euclid = params%lp_iters
            call cline%set('needs_sigma','yes')
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        niters       = 0
        iter         = params%startit - 1
        corr         = -1.
        do
            if( L_BENCH_GLOB )then
                t_init = tic()
                t_tot  = t_init
            endif
            niters            = niters + 1
            iter              = iter + 1
            params%which_iter = iter
            str_iter          = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
            if( l_switch2euclid .or. trim(params%objfun).eq.'euclid' )then
                call cline_calc_group_sigmas%set('which_iter',real(iter))
                call qenv%exec_simple_prg_in_queue(cline_calc_group_sigmas, 'CALC_GROUP_SIGMAS_FINISHED')
            endif
            if( have_oris .or. iter > params%startit )then
                call build%spproj%read(params%projfile)
                if( params%refine .eq. 'snhc' )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = build%spproj_field%get_avg('corr')
                    if( iter > 1 .and. corr <= corr_prev )then
                        params%szsn = min(SZSN_MAX,params%szsn + SZSN_STEP)
                    endif
                    call job_descr%set('szsn', int2str(params%szsn))
                    call cline%set('szsn', real(params%szsn))
                endif
            endif
            ! exponential cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set( 'extr_iter',  trim(int2str(params%extr_iter)))
            call cline%set(     'extr_iter',  real(params%extr_iter))
            call job_descr%set( 'which_iter', trim(int2str(params%which_iter)))
            call cline%set(     'which_iter', real(params%which_iter))
            call job_descr%set( 'startit',    trim(int2str(iter)))
            call cline%set(     'startit',    real(iter))
            ! FRCs
            if( cline%defined('frcs') )then
                ! all good
            else
                if( l_multistates )call job_descr%set('frcs', trim(FRCS_FILE))
            endif
            ! schedule
            if( L_BENCH_GLOB )then
                rt_init = toc(t_init)
                t_scheduled = tic()
            endif
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
            ! assemble alignment docs
            if( L_BENCH_GLOB )then
                rt_scheduled = toc(t_scheduled)
                t_merge_algndocs = tic()
            endif
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY)
            if( L_BENCH_GLOB )then
                rt_merge_algndocs = toc(t_merge_algndocs)
                t_volassemble = tic()
            endif
            ! ASSEMBLE VOLUMES
            select case(trim(params%refine))
                case('eval')
                    ! nothing to do
                case DEFAULT
                    call cline_volassemble%set( 'prg', 'volassemble' ) ! required for cmdline exec
                    call cline_volassemble%set( 'which_iter', int2str(params%which_iter) )
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
                        if( str_has_substr(params%refine,'snhc') )then
                            volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
                        else
                            volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
                        endif
                        call cline_volassemble%set( 'state', real(state) )
                        if( params%nstates>1 )call cline_volassemble%set('part', real(state))
                        call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                        &'simple_script_state'//trim(str_state), volassemble_output)
                    end do
                    call qsys_watcher(state_assemble_finished)
                    ! rename & add volumes to project & update job_descr
                    call build%spproj_field%get_pops(state_pops, 'state')
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
                        if( state_pops(state) == 0 )then
                            ! cleanup for empty state
                            vol = 'vol'//trim(int2str(state))
                            call cline%delete( vol )
                            call job_descr%delete( vol )
                        else
                            ! rename state volume
                            vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                            if( params%refine .eq. 'snhc' )then
                                vol_iter = trim(SNHCVOL)//trim(str_state)//params%ext
                                call simple_rename( vol, vol_iter )
                            else
                                vol_iter = trim(vol)
                            endif
                            fsc_file = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                            call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                            ! generate FSC pdf
                            res = get_resarr(params%box, params%smpd)
                            fsc = file2rarr(fsc_file)
                            if( str_has_substr(params%refine,'snhc') )then
                                fsc_templ = 'fsc_state'//trim(str_state)
                            else
                                fsc_templ = 'fsc_state'//trim(str_state)//'_iter'//trim(str_iter)
                            endif
                            call plot_fsc(size(fsc), fsc, res, params%smpd, fsc_templ)
                            ! add state volume to os_out
                            if( trim(params%oritype).eq.'cls3D' )then
                                call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol_cavg')
                            else
                                call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol')
                            endif
                            ! updates cmdlines & job description
                            vol = 'vol'//trim(int2str(state))
                            call job_descr%set( vol, vol_iter )
                            call cline%set(vol, vol_iter)
                            if( params%keepvol.ne.'no' )then
                                call simple_copy_file(vol_iter,trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//params%ext)
                                ! keeping odd/even vols
                                call simple_copy_file(vol_odd, trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//'_odd' //params%ext)
                                call simple_copy_file(vol_even,trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//'_even'//params%ext)
                            endif
                        endif
                    enddo
                    ! volume mask, one for all states
                    if( cline%defined('mskfile') )then
                        if( file_exists(trim(params%mskfile)) )then
                            call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                        endif
                    endif
                    ! writes os_out
                    call build%spproj%write_segment_inside('out')
                    ! automasking in postprocess
                    if( params%l_automsk )then
                        if( mod(iter,AUTOMSK_FREQ) == 0 .or. iter == params%startit )then
                            call cline_postprocess%delete('mskfile')
                            call cline_postprocess%set('automsk', trim(params%automsk))
                        endif
                    endif
                    ! per state post-process
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
                        if( state_pops(state) == 0 ) cycle
                        call cline_postprocess%set('state', real(state))
                        if( cline%defined('lp') ) call cline_postprocess%set('lp', params%lp)
                        call xpostprocess%execute(cline_postprocess)
                        ! for gui visualization
                        if( params%refine .ne. 'snhc' )then
                            volpproc = trim(VOL_FBODY)//trim(str_state)//PPROC_SUFFIX//params%ext
                            vollp    = trim(VOL_FBODY)//trim(str_state)//LP_SUFFIX//params%ext
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//PPROC_SUFFIX//params%ext
                            call simple_copy_file(volpproc, vol_iter)
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//LP_SUFFIX//params%ext
                            call simple_copy_file(vollp, vol_iter)
                            if( iter > 1 )then
                                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter-1,3)//PPROC_SUFFIX//params%ext
                                call del_file(vol_iter)
                                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter-1,3)//LP_SUFFIX//params%ext
                                call del_file(vol_iter)
                            endif
                        endif
                    enddo
                    ! update command-lines to use the mskfile for the next AUTOMSK_FREQ - 1 iterations
                    if( params%l_automsk )then
                        if( mod(iter,AUTOMSK_FREQ) == 0 .or. iter == params%startit )then
                            params%mskfile = 'automask'//params%ext
                            call cline_postprocess%set('mskfile', trim(params%mskfile))
                            call cline_postprocess%delete('automsk')
                            call cline%set('mskfile', trim(params%mskfile))
                            call job_descr%set('mskfile', trim(params%mskfile))
                        endif
                    endif
            end select
            if( L_BENCH_GLOB ) rt_volassemble = toc(t_volassemble)
            ! CONVERGENCE
            converged = .false.
            select case(trim(params%refine))
                case('eval')
                    ! nothing to do
                case DEFAULT
                    if( str_has_substr(params%refine,'cluster')) call cline_check_3Dconv%delete('update_res')
                    call xcheck_3Dconv%execute(cline_check_3Dconv)
                    if( iter >= params%startit + 2 )then
                        ! after a minimum of 2 iterations
                        if( cline_check_3Dconv%get_carg('converged') .eq. 'yes' ) converged = .true.
                    endif
            end select
            if( iter >= params%maxits ) converged = .true.
            if ( l_combine_eo .and. converged )then
                converged         = .false.
                l_combine_eo      = .false.
                params%combine_eo = 'yes'
                params%maxits     = iter + 1
                params%lplim_crit = min(0.143,params%lplim_crit)
                call cline%set('lplim_crit',params%lplim_crit)
                call job_descr%set('lplim_crit',real2str(params%lplim_crit))
                call cline_volassemble%set('combine_eo', 'yes')
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> PERFORMING FINAL ITERATION WITH COMBINED EVEN/ODD VOLUMES & RESOLUTION LIMIT BEYOND FSC=0.143'
                call simple_copy_file(trim(VOL_FBODY)//trim(str_state)//params%ext, trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext)
                call simple_copy_file(trim(VOL_FBODY)//trim(str_state)//params%ext, trim(VOL_FBODY)//trim(str_state)//'_odd'//params%ext)
            endif
            if( converged )then
                if(trim(params%oritype).eq.'cls3D') call build%spproj%map2ptcls
                ! safest to write the whole thing here as multiple fields updated
                call build%spproj%write
                exit ! main loop
            endif
            ! ITERATION DEPENDENT UPDATES
            if( cline_check_3Dconv%defined('trs') .and. .not.job_descr%isthere('trs') )then
                ! activates shift search if frac_srch >= 90
                str = real2str(cline_check_3Dconv%get_rarg('trs'))
                call job_descr%set( 'trs', trim(str) )
                call cline%set( 'trs', cline_check_3Dconv%get_rarg('trs') )
            endif
            if( l_switch2eo .and. (niters == params%lp_iters ) )then
                ! e/o projection matching
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO EVEN/ODD RESOLUTION LIMIT'
                l_projmatch = .false.
                if( cline%defined('match_filt') )then
                    if( cline%get_carg('match_filt').eq.'no' )then
                        ! flags are kept so match_filt is not used
                        call job_descr%set('match_filt','no')
                    else
                        call cline%delete('lp')
                        call job_descr%delete('lp')
                        call cline_postprocess%delete('lp')
                        call cline%delete('lp_iters')
                        call job_descr%delete('lp_iters')
                    endif
                else
                    call cline%delete('lp')
                    call job_descr%delete('lp')
                    call cline_postprocess%delete('lp')
                    call cline%delete('lp_iters')
                    call job_descr%delete('lp_iters')
                endif
            endif
            if( l_projmatch )then
                if( params%l_frac_update )then
                    call job_descr%set('update_frac', real2str(params%update_frac))
                    call cline%set('update_frac', params%update_frac)
                    call cline_check_3Dconv%set('update_frac', params%update_frac)
                    call cline_volassemble%set('update_frac', params%update_frac)
                endif
            endif
            ! objfun=euclid, part 3: actual switch
            if( l_switch2euclid .and. niters.eq.iter_switch2euclid )then
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO OBJFUN=EUCLID'
                call cline%set('objfun',    'euclid')
                call cline%set('match_filt','no')
                if(.not.l_griddingset )then
                    call cline%set('gridding',     'yes')
                    call job_descr%set('gridding', 'yes')
                endif
                call cline%delete('lp')
                call job_descr%set('objfun',    'euclid')
                call job_descr%set('match_filt','no')
                call job_descr%delete('lp')
                call cline_volassemble%set('objfun','euclid')
                call cline_postprocess%delete('lp')
                if( l_ptclw )then
                    call cline%set('ptclw',    'yes')
                    call job_descr%set('ptclw','yes')
                endif
                params%objfun    = 'euclid'
                params%cc_objfun = OBJFUN_EUCLID
                l_switch2euclid  = .false.
            endif
            if( L_BENCH_GLOB )then
                rt_tot  = toc(t_init)
                benchfname = 'REFINE3D_DISTR_BENCH_ITER'//int2str_pad(iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', rt_scheduled
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', rt_merge_algndocs
                write(fnr,'(a,1x,f9.2)') 'volassemble     : ', rt_volassemble
                write(fnr,'(a,1x,f9.2)') 'total time      : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation  : ', (rt_init/rt_tot)           * 100.
                write(fnr,'(a,1x,f9.2)') 'scheduled jobs  : ', (rt_scheduled/rt_tot)      * 100.
                write(fnr,'(a,1x,f9.2)') 'merge_algndocs  : ', (rt_merge_algndocs/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'volassemble     : ', (rt_volassemble/rt_tot)    * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for : ',&
                    &((rt_init+rt_scheduled+rt_merge_algndocs+rt_volassemble)/rt_tot) * 100.
                call fclose(fnr)
            endif
        end do
        ! put back automsk flag if needed
        if( params%l_automsk ) call cline%set('automsk', trim(params%automsk))
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_DISTR_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D_distr

    subroutine exec_refine3D( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(refine3D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(parameters)         :: params
        type(builder)            :: build
        integer                  :: startit, i, state
        character(len=STDLEN)    :: str_state, fsc_file, vol, vol_iter
        logical                  :: converged, l_automsk
        real                     :: corr, corr_prev
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        l_automsk = params%l_automsk
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        if( startit == 1 ) call build%spproj_field%clean_updatecnt
        if( params%l_distr_exec )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
            call refine3D_exec(cline, startit, converged)
        else
            if( trim(params%objfun) == 'euclid'   ) THROW_HARD('shared-memory implementation of refine3D does not support objfun=euclid')
            if( trim(params%continue) == 'yes'    ) THROW_HARD('shared-memory implementation of refine3D does not support continue=yes')
            if( .not. file_exists(params%vols(1)) ) THROW_HARD('shared-memory implementation of refine3D requires starting volume(s) input')
            ! take care of automask flag
            if( cline%defined('automsk') ) call cline%delete('automsk')
            if( params%l_automsk .and. params%nstates > 1 ) THROW_HARD('automsk.ne.no not currenty supported for multi-state refinement')
            params%startit     = startit
            params%outfile     = 'algndoc'//METADATA_EXT
            params%extr_iter   = params%startit - 1
            corr               = -1.
            do i = 1, params%maxits
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%startit
                write(logfhandle,'(A)')   '>>>'
                if( params%refine .eq. 'snhc' .and. params%startit > 1 )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = build%spproj_field%get_avg('corr')
                    if( corr <= corr_prev ) params%szsn = min(SZSN_MAX,params%szsn + SZSN_STEP)
                endif
                ! exponential cooling of the randomization rate
                params%extr_iter = params%extr_iter + 1
                ! to control the masking behaviour in simple_strategy2D3D_common :: norm_struct_facts
                if( l_automsk )then
                    if( mod(params%startit,AUTOMSK_FREQ) == 0 .or. i == 1 )then
                        call cline%set('automsk', trim(params%automsk))
                        params%l_automsk = .true.
                    else
                        call cline%delete('automsk')
                        params%l_automsk = .false.
                    endif
                endif
                ! in strategy3D_matcher:
                call refine3D_exec(cline, params%startit, converged)
                if( converged .or. i == params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', real(params%startit))
                    ! update project with the new orientations
                    call build%spproj%write_segment_inside(params%oritype)
                    call del_file(params%outfile)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        fsc_file = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                        call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                        ! add state volume to os_out
                        vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                        if( params%refine .eq. 'snhc' )then
                            vol_iter = trim(SNHCVOL)//trim(str_state)//params%ext
                        else
                            vol_iter = trim(vol)
                        endif
                        if( trim(params%oritype).eq.'cls3D' )then
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol_cavg')
                        else
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol')
                        endif! volume mask, one for all states
                    end do
                    if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                    call build%spproj%write_segment_inside('out')
                    exit
                endif
                ! update iteration counter
                params%startit = params%startit + 1
            end do
            ! put back automsk flag if needed
            if( l_automsk ) call cline%set('automsk', trim(params%automsk))
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    subroutine exec_check_3Dconv( self, cline )
        use simple_convergence, only: convergence
        use simple_parameters,  only: params_glob
        class(check_3Dconv_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)  :: params
        type(builder)     :: build
        type(convergence) :: conv
        real, allocatable :: maplp(:)
        integer           :: istate, loc(1)
        logical           :: converged, update_res
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        update_res = .false.
        allocate( maplp(params%nstates))
        maplp = 0.
        do istate=1,params%nstates
            if( build%spproj_field%get_pop( istate, 'state' ) == 0 )cycle ! empty state
            params%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
            if( file_exists(params%fsc) )then
                build%fsc(istate,:) = file2rarr(params%fsc)
                maplp(istate)   = max(build%img%get_lp(get_lplim_at_corr(build%fsc(istate,:),params%lplim_crit)),2.*params%smpd)
            else
                THROW_HARD('tried to check the fsc file: '//trim(params%fsc)//' but it does not exist!')
            endif
        enddo
        loc     = maxloc( maplp )
        params%state = loc(1)                ! state with worst low-pass
        params%lp    = maplp( params%state ) ! worst lp
        params%fsc   = 'fsc_state'//int2str_pad(params%state,2)//'.bin'
        deallocate(maplp)
        ! check convergence
        if( cline%defined('update_res') )then
            update_res = .false.
            if( cline%get_carg('update_res').eq.'yes' )update_res = .true.
            if( cline%get_carg('update_res').eq.'no' .and. str_has_substr(params%refine,'cluster') )then
                converged = conv%check_conv_cluster(cline)
            else
                if( params%l_cartesian )then
                    converged = conv%check_conv3Dc(cline, params%msk)
                else
                    converged = conv%check_conv3D(cline, params%msk)
                endif
            endif
        else
            select case(params%refine)
                case('cluster','clustersym')
                    converged = conv%check_conv_cluster(cline)
                case('shcc','neighc','greedyc')
                    converged = conv%check_conv3Dc(cline, params%msk)
                case DEFAULT
                    converged = conv%check_conv3D(cline, params%msk)
            end select
        endif
        ! reports convergence, shift activation, resolution update and
        ! fraction of search space scanned to the distr commander
        if( params_glob%l_doshift )then
            call cline%set('trs', params_glob%trs) ! activates shift search
        endif
        if( converged )then
            call cline%set('converged', 'yes')
        else
            call cline%set('converged', 'no')
        endif
        if( update_res )then
            call cline%set('update_res', 'yes') ! fourier index to be updated in distr commander
        else
            call cline%set('update_res', 'no')
        endif
        call cline%set('frac_srch', conv%get('frac_srch'))
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CHECK_3DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_3Dconv

end module simple_commander_refine3D
