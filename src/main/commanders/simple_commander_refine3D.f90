! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_builder,          only: builder, build_glob
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_parameters,       only: parameters
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_qsys_env,         only: qsys_env
use simple_cluster_seed,     only: gen_labelling
use simple_commander_volops, only: postprocess_commander
use simple_starproject,      only: starproject
use simple_commander_euclid
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_commander_distr
public :: refine3D_commander
public :: check_3Dconv_commander
public :: check_align_commander
public :: check_align_inpl_commander
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

type, extends(commander_base) :: check_align_commander
  contains
    procedure :: execute      => exec_check_align
end type check_align_commander

type, extends(commander_base) :: check_align_inpl_commander
  contains
    procedure :: execute      => exec_check_align_inpl
end type check_align_inpl_commander

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
        type(reconstruct3D_commander_distr)   :: xreconstruct3D_distr
        type(calc_pspec_commander_distr)      :: xcalc_pspec_distr
        type(check_3Dconv_commander)          :: xcheck_3Dconv
        type(postprocess_commander)           :: xpostprocess
        type(refine3D_commander)              :: xrefine3D_shmem
        type(estimate_first_sigmas_commander) :: xfirst_sigmas
        ! command lines
        type(cmdline)    :: cline_reconstruct3D_distr
        type(cmdline)    :: cline_calc_pspec_distr
        type(cmdline)    :: cline_calc_sigma
        type(cmdline)    :: cline_check_3Dconv
        type(cmdline)    :: cline_volassemble
        type(cmdline)    :: cline_postprocess
        integer(timer_int_kind) :: t_init,   t_scheduled,  t_merge_algndocs,  t_volassemble,  t_tot
        real(timer_int_kind)    :: rt_init, rt_scheduled, rt_merge_algndocs, rt_volassemble, rt_tot
        character(len=STDLEN)   :: benchfname
        ! other variables
        type(parameters)    :: params
        type(builder)       :: build
        type(qsys_env)      :: qenv
        type(chash)         :: job_descr
        type(starproject)   :: starproj
        character(len=:),          allocatable :: vol_fname, prev_refine_path, target_name
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        integer,                   allocatable :: state_pops(:)
        real,                      allocatable :: res(:), fsc(:)
        character(len=STDLEN)     :: vol, vol_iter, str, str_iter, fsc_templ, orig_objfun
        character(len=STDLEN)     :: vol_even, vol_odd, str_state, fsc_file, volpproc, vollp
        character(len=LONGSTRLEN) :: volassemble_output
        logical :: err, vol_defined, have_oris, do_abinitio, converged, fall_over
        logical :: l_projmatch, l_switch2eo, l_switch2euclid, l_continue, l_multistates
        logical :: l_combine_eo, l_lpset, l_griddingset ! l_ptclw,
        real    :: corr, corr_prev, smpd
        integer :: ldim(3), i, state, iter, box, nfiles, niters, iter_switch2euclid, ifoo
        integer :: fnr
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
            if( trim(cline%get_carg('objfun')).eq.'euclid' .and. .not.l_continue )then
                orig_objfun     = trim(cline%get_carg('objfun'))
                l_switch2euclid = .true.
                call cline%set('objfun','cc')
                ! l_ptclw = trim(cline%get_carg('ptclw')).eq.'yes'
                ! call cline%set('ptclw', 'no')
            endif
        endif
        ! init
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
        ! randomized oris and zero shifts when reg_ref is on
        if( params%l_reg_init )then
            call build%spproj_field%rnd_oris
            call build%spproj_field%zero_shifts
            write(logfhandle,'(A)')   '>>> APPLYING RANDOMIZED ORIS AND ZERO SHIFTS'
        endif
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
        cline_calc_sigma          = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' )     ! required for distributed call
        call cline_calc_pspec_distr%set(    'prg', 'calc_pspec' )        ! required for distributed call
        call cline_postprocess%set(         'prg', 'postprocess' )       ! required for local call
        call cline_calc_sigma%set(          'prg', 'calc_group_sigmas' ) ! required for local call
        if( trim(params%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set('pgrp', 'c1')
        call cline_postprocess%set('mirr',    'no')
        call cline_postprocess%set('mkdir',   'no')
        call cline_postprocess%set('imgkind', 'vol')
        if( trim(params%oritype).eq.'cls3D' ) call cline_postprocess%set('imgkind', 'vol_cavg')
        call cline_postprocess%delete('nonuniform') ! to save time in the iterations
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
            if( trim(params%objfun).eq.'euclid' )then
                call cline%set('needs_sigma','yes')
                call cline_reconstruct3D_distr%set('needs_sigma','yes')
                call cline_volassemble%set('needs_sigma','yes')
                if( .not.l_griddingset .and. .not.params%l_cartesian ) call cline%set('gridding','yes')
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
                if( (ldim(1) /= params%box_crop) .and.  (ldim(1) /= params%box) )then
                    THROW_HARD('Incompatible dimensions between input volume and images: '//params%vols(state))
                endif
            endif
        enddo
        have_oris   = .not. build%spproj%is_virgin_field(params%oritype)
        do_abinitio = .not. have_oris .and. .not. vol_defined
        if( do_abinitio )then
            call build%spproj_field%rnd_oris
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
            if( l_switch2euclid )then
                ! first, estimate group sigmas
                call cline_calc_sigma%set('which_iter', real(params%startit))
                call qenv%exec_simple_prg_in_queue(cline_calc_sigma, 'CALC_GROUP_SIGMAS_FINISHED')
                ! then, estimate first sigmas given reconstructed starting volumes(s) and previous orientations
                if( .not.cline%defined('nspace') ) call cline%set('nspace', real(params%nspace))
                if( .not.cline%defined('athres') ) call cline%set('athres', real(params%athres))
                call xfirst_sigmas%execute(cline)
                ! update command lines
                call cline%set('needs_sigma','yes')
                call cline_volassemble%set('needs_sigma','yes')
                call cline_reconstruct3D_distr%set('needs_sigma','yes')
                call cline%set('objfun', orig_objfun)
                params%objfun   = trim(orig_objfun)
                l_switch2euclid = .false.
            endif
        else if( vol_defined .and. params%continue .ne. 'yes' )then
            ! projection matching
            l_projmatch = .true.
            if( .not.l_lpset )then
                THROW_HARD('LP needs be defined for the first step of projection matching!')
                call cline%delete('update_frac')
            endif
            if( (.not.str_has_substr(params%refine, 'neigh')) )then
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
                call cline_calc_sigma%set('which_iter',real(iter))
                call qenv%exec_simple_prg_in_queue(cline_calc_sigma, 'CALC_GROUP_SIGMAS_FINISHED')
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
                            res = get_resarr(params%box_crop, params%smpd_crop)
                            fsc = file2rarr(fsc_file)
                            if( str_has_substr(params%refine,'snhc') )then
                                fsc_templ = 'fsc_state'//trim(str_state)
                            else
                                fsc_templ = 'fsc_state'//trim(str_state)//'_iter'//trim(str_iter)
                            endif
                            call plot_fsc(size(fsc), fsc, res, params%smpd_crop, fsc_templ)
                            ! add state volume to os_out
                            if( trim(params%oritype).eq.'cls3D' )then
                                call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                            else
                                call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                            endif
                            ! updates cmdlines & job description
                            vol = 'vol'//trim(int2str(state))
                            call job_descr%set( vol, vol_iter )
                            call cline%set(vol, vol_iter)
                            if( params%keepvol.ne.'no' )then
                                call simple_copy_file(vol_iter,trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//params%ext)
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
                        if( mod(niters,AUTOMSK_FREQ) == 0 .or. iter == params%startit )then
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
                        if( mod(niters,AUTOMSK_FREQ) == 0 .or. iter == params%startit )then
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
            if( niters == params%maxits ) converged = .true.
            if ( l_combine_eo .and. converged )then
                converged            = .false.
                l_combine_eo         = .false.
                params%combine_eo    = 'yes'
                params%l_frac_update = .false.
                params%update_frac   = 1.0
                params%maxits        = niters + 1
                params%lplim_crit    = min(0.143,params%lplim_crit)
                call cline%set('lplim_crit',params%lplim_crit)
                call cline%set('update_frac',1.0)
                call job_descr%set('lplim_crit',real2str(params%lplim_crit))
                call job_descr%set('update_frac',real2str(1.0))
                call cline_volassemble%set('combine_eo', 'yes')
                call cline_volassemble%set('update_frac', 1.0)
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
                call cline%delete('lp')
                call job_descr%delete('lp')
                call cline_postprocess%delete('lp')
                call cline%delete('lp_iters')
                call job_descr%delete('lp_iters')
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
                call cline%set('objfun', orig_objfun)
                if(.not.l_griddingset .and. .not.params%l_cartesian )then
                    call cline%set('gridding',     'yes')
                    call job_descr%set('gridding', 'yes')
                endif
                call job_descr%set('objfun', orig_objfun)
                call cline_volassemble%set('objfun', orig_objfun)
                if( l_switch2eo )then
                    ! delete resolution limit
                    call cline%delete('lp')
                    call job_descr%delete('lp')
                    call cline_postprocess%delete('lp')
                endif
                ! if( l_ptclw )then
                !     call cline%set('ptclw',    'yes')
                !     call job_descr%set('ptclw','yes')
                ! endif
                params%objfun = orig_objfun
                select case(trim(params%objfun))
                    case('euclid')
                        params%cc_objfun = OBJFUN_EUCLID
                    case('prob')
                        params%cc_objfun = OBJFUN_PROB
                end select
                l_switch2euclid = .false.
            endif
            ! write per iteration star file
            call starproj%export_iter3D(build%spproj, params%nstates,  params%which_iter)
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
        use simple_regularizer_inpl,   only: reg_params
        class(refine3D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(reg_params),          allocatable   :: cur_tab(:,:,:)
        type(calc_group_sigmas_commander) :: xcalc_group_sigmas
        type(calc_pspec_commander_distr)  :: xcalc_pspec_distr
        type(parameters)                  :: params
        type(builder)                     :: build
        type(cmdline)                     :: cline_calc_sigma, cline_calc_pspec_distr
        character(len=STDLEN)             :: str_state, fsc_file, vol, vol_iter, orig_objfun
        integer                           :: startit, i, state
        real                              :: corr, corr_prev
        logical                           :: converged, l_automsk, l_sigma, l_switch2euclid
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        l_automsk = params%l_automsk
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        if( startit == 1 ) call build%spproj_field%clean_updatecnt
        if( params%l_distr_exec )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
            call refine3D_exec(cline, startit, converged)
        else
            if( trim(params%continue) == 'yes'    ) THROW_HARD('shared-memory implementation of refine3D does not support continue=yes')
            if( .not. file_exists(params%vols(1)) ) THROW_HARD('shared-memory implementation of refine3D requires starting volume(s) input')
            ! objfun=euclid|prob
            orig_objfun     = trim(params%objfun)
            l_sigma         = .false.
            l_switch2euclid = .false.
            select case(trim(orig_objfun))
            case('euclid')
                l_sigma = .true.
                call cline%set('needs_sigma','yes')
                params%l_needs_sigma = .true.
                cline_calc_sigma     = cline
                if( file_exists(trim(SIGMA2_GROUP_FBODY)//trim(int2str(params%which_iter))//'.star') )then
                    ! it is assumed that we already have precalculted sigmas2 and all corresponding flags have been set
                else
                    ! sigma2 not provided & are calculated
                    if( build%spproj_field%get_nevenodd() == 0 )then
                        ! make sure we have e/o partitioning prior to calc_pspec_distr
                        call build%spproj_field%partition_eo
                        call build%spproj%write_segment_inside(params%oritype)
                    endif
                    params%objfun    = 'cc'
                    params%cc_objfun = OBJFUN_CC
                    cline_calc_pspec_distr = cline
                    call cline_calc_pspec_distr%set('prg', 'calc_pspec' )
                    call xcalc_pspec_distr%execute( cline_calc_pspec_distr )
                    l_switch2euclid = .true.
                endif
            case('cc')
                ! nothing to do
            end select
            ! take care of automask flag
            if( cline%defined('automsk') ) call cline%delete('automsk')
            if( params%l_automsk .and. params%nstates > 1 ) THROW_HARD('automsk.ne.no not currenty supported for multi-state refinement')
            params%startit     = startit
            params%which_iter  = params%startit
            params%outfile     = 'algndoc'//METADATA_EXT
            params%extr_iter   = params%startit - 1
            corr               = -1.
            allocate(cur_tab(params%fromp:params%top,params%nspace,params%reg_nrots))
            cur_tab(:,:,:)%prob = 0.
            do i = 1, params%maxits
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
                write(logfhandle,'(A)')   '>>>'
                if( params%refine .eq. 'snhc' .and. params%which_iter > 1 )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = build%spproj_field%get_avg('corr')
                    if( corr <= corr_prev ) params%szsn = min(SZSN_MAX,params%szsn + SZSN_STEP)
                endif
                ! exponential cooling of the randomization rate
                params%extr_iter = params%extr_iter + 1
                ! to control the masking behaviour in simple_strategy2D3D_common :: norm_struct_facts
                if( l_automsk )then
                    if( mod(params%which_iter,AUTOMSK_FREQ) == 0 .or. i == 1 )then
                        call cline%set('automsk', trim(params%automsk))
                        params%l_automsk = .true.
                    else
                        call cline%delete('automsk')
                        params%l_automsk = .false.
                    endif
                endif
                if( l_sigma )then
                    call cline_calc_sigma%set('which_iter',real(i))
                    call xcalc_group_sigmas%execute(cline_calc_sigma)
                endif
                ! in strategy3D_matcher:
                call refine3D_exec(cline, params%which_iter, converged, cur_tab)
                if( converged .or. i == params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', real(params%which_iter))
                    ! update project with the new orientations
                    call build%spproj%write_segment_inside(params%oritype)
                    call del_file(params%outfile)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        fsc_file = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                        call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                        ! add state volume to os_out
                        vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                        if( params%refine .eq. 'snhc' )then
                            vol_iter = trim(SNHCVOL)//trim(str_state)//params%ext
                        else
                            vol_iter = trim(vol)
                        endif
                        if( trim(params%oritype).eq.'cls3D' )then
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                        else
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                        endif! volume mask, one for all states
                    end do
                    if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                    call build%spproj%write_segment_inside('out')
                    exit
                endif
                ! update iteration counter
                params%which_iter = params%which_iter + 1
                ! whether to switch objective function
                if( l_switch2euclid )then
                    params%objfun = trim(orig_objfun)
                    select case(trim(params%objfun))
                        case('euclid')
                            params%cc_objfun = OBJFUN_EUCLID
                        case('prob')
                            params%cc_objfun = OBJFUN_PROB
                    end select
                    if( .not.cline%defined('gridding') )then
                        call cline%set('gridding', 'yes')
                        params%gridding = 'yes'
                    endif
                    l_switch2euclid = .false.
                endif
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
        integer           :: istate
        logical           :: converged, update_res
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        update_res = .false.
        allocate( maplp(params%nstates), source=0.)
        do istate=1,params%nstates
            if( build%spproj_field%get_pop( istate, 'state' ) == 0 )cycle ! empty state
            params%fsc = 'fsc_state'//int2str_pad(istate,2)//'.bin'
            if( file_exists(params%fsc) )then
                build%fsc(istate,:) = file2rarr(params%fsc)
                maplp(istate) = calc_lowpass_lim(get_lplim_at_corr(build%fsc(istate,:),params%lplim_crit), params_glob%box_crop, params_glob%smpd_crop)
                maplp(istate) = max(maplp(istate), 2.*params%smpd_crop)
            else
                THROW_HARD('tried to check the fsc file: '//trim(params%fsc)//' but it does not exist!')
            endif
        enddo
        params%state = maxloc(maplp, dim=1)  ! state with worst low-pass
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
                case('shcc','neighc','greedyc','hybrid')
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

    subroutine exec_check_align( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, prepimg4align, calcrefvolshift_and_mapshifts2ptcls,&
                    &read_and_filter_refvols, preprefvol, preprecvols, norm_struct_facts, killrecvols, grid_ptcl
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_fplane,              only: fplane
        use simple_regularizer,         only: regularizer
        use simple_image
        class(check_align_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        integer,          parameter   :: MAXITS = 60
        integer,          allocatable :: pinds(:), best_ip(:), best_ir(:)
        logical,          allocatable :: ptcl_mask(:)
        complex,          allocatable :: cmat(:,:)
        real,             allocatable :: sigma2_noise(:,:)
        type(image),      allocatable :: tmp_imgs(:)
        type(fplane),     allocatable :: fpls(:)
        type(ctfparams),  allocatable :: ctfparms(:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(ori)                     :: o_tmp
        type(image)                   :: img
        type(regularizer)             :: reg_obj
        type(ori)                     :: orientation
        integer  :: nptcls, iptcl, j, s, iref, box, loc, pind_here, ithr
        logical  :: l_ctf, do_center
        real     :: xyz(3), euls(3), shvec(2), sdev
        call cline%set('mkdir',    'yes')
        call cline%set('oritype',  'ptcl3D')
        call cline%set('center',   'no')
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        call build%spproj%update_projinfo(cline)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        allocate(ptcl_mask(params%fromp:params%top))
        call build_glob%spproj_field%sample4update_and_incrcnt([params%fromp,params%top],&
            &1.0, nptcls, pinds, ptcl_mask)
        print *, 'nptcls = ', nptcls, '; fromp = ', params%fromp, '; top = ', params%top
        call pftcc%new(params%nspace, [1,nptcls], params%kfromto)
        call pftcc%reallocate_ptcls(nptcls, pinds)
        call reg_obj%new(pftcc)
        print *, 'Preparing the references ...'
        ! e/o partioning
        if( build%spproj%os_ptcl3D%get_nevenodd() == 0 )then
            call build%spproj%os_ptcl3D%partition_eo
            call build%spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        do s=1,params%nstates
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params%vols(s), do_center, xyz)
            call read_and_filter_refvols( cline, params%vols(s), params%vols(s))
            ! PREPARE E/O VOLUMES
            call preprefvol(cline, s, do_center, xyz, .false.)
            call preprefvol(cline, s, do_center, xyz, .true.)
            ! PREPARE REFERENCES
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1, params%nspace
                call build_glob%eulspace%get_ori(iref, o_tmp)
                call build_glob%vol_odd%fproject_polar((s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.false., mask=build_glob%l_resmsk)
                call build_glob%vol%fproject_polar(    (s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.true.,  mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
        end do
        call pftcc%memoize_refs
        ! PREPARATION OF PARTICLES
        print *, 'Preparing the particles ...'
        call prepimgbatch(params%top-params%fromp+1)
        call read_imgbatch([params%fromp,params%top], ptcl_mask)
        call build%img_crop_polarizer%init_polarizer(pftcc, params%alpha)
        allocate(tmp_imgs(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call tmp_imgs(ithr)%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(iptcl,ithr) schedule(static) proc_bind(close)
        do iptcl = params%fromp,params%top
            if( .not.ptcl_mask(iptcl) ) cycle
            ithr = omp_get_thread_num()+1
            ! prep
            call prepimg4align(iptcl, build%imgbatch(iptcl), tmp_imgs(ithr))
            call build%imgbatch(iptcl)%ifft ! for reconstruction
            ! transfer to polar coordinates
            call build%img_crop_polarizer%polarize(pftcc, tmp_imgs(ithr), iptcl, .true., .true., mask=build%l_resmsk)
            ! e/o flags
            ! call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
            call pftcc%set_eo(iptcl, .true.)
        enddo
        !$omp end parallel do
        ! getting the ctfs
        l_ctf = build%spproj%get_ctfflag('ptcl2D',iptcl=pinds(1)).ne.'no'
        ! make CTFs
        if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, 'ptcl2D')

        ! ALIGNMENT OF PARTICLES
        print *, 'Aligning the particles ...'
        ! using gencorrs (cc-based to estimate the sigma)
        if( params%l_needs_sigma )then
            allocate( sigma2_noise(pftcc%kfromto(1):pftcc%kfromto(2), params%fromp:params%top), source=1. )
            ! do j = pftcc%kfromto(1),pftcc%kfromto(2)
            !     sigma2_noise(j,:) = real(j)
            ! enddo
            call pftcc%assign_sigma2_noise(sigma2_noise)
            ! call pftcc%memoize_ptcls
            ! params%cc_objfun = OBJFUN_CC
            ! !$omp parallel do default(shared) private(j,iref,ithr,iptcl,inpl_corrs,cxy,max_corr,max_iref,max_sh,max_loc,loc,corr,sh) proc_bind(close) schedule(static)
            ! do j = 1, nptcls
            !     max_corr = 0.
            !     do iref = 1, params%nspace
            !         ithr  = omp_get_thread_num() + 1
            !         iptcl = pinds(j)
            !         ! find best irot/shift for this pair of iref, iptcl
            !         call pftcc%gencorrs( iref, iptcl, inpl_corrs )
            !         loc = maxloc(inpl_corrs, dim=1)
            !         call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
            !         cxy = grad_shsrch_obj(ithr)%minimize(irot=loc)
            !         if( loc > 0 )then
            !             corr = cxy(1)
            !             sh   = cxy(2:3)
            !         else
            !             loc  = maxloc(inpl_corrs, dim=1)
            !             corr = inpl_corrs(loc)
            !             sh   = 0.
            !         endif
            !         if( corr > max_corr )then
            !             max_corr = corr
            !             max_loc  = loc
            !             max_sh   = sh
            !             max_iref = iref
            !         endif
            !     enddo
            !     call pftcc%update_sigma( max_iref, iptcl, max_sh, max_loc )
            ! enddo
            ! !$omp end parallel do
            ! params%cc_objfun = OBJFUN_EUCLID
        endif
        ! actual alignment using the defined cost function
        ! scaling by the ctf
        if( params%l_reg_scale )then
            call pftcc%reg_scale
            !$omp parallel do default(shared) private(j) proc_bind(close) schedule(static)
            do j = 1, nptcls
                call pftcc%memoize_sqsum_ptcl(pinds(j))
            enddo
            !$omp end parallel do
        endif
        ! Memoize particles FFT parameters
        call pftcc%memoize_ptcls
        call reg_obj%init_tab
        call reg_obj%fill_tab(pinds)
        print *, 'Assembling the class averages with'
        select case(trim(params%reg_mode))
            case('tab')
                print *, 'soft-sorting the tab...'
                call reg_obj%sort_tab_no_norm
                call reg_obj%ref_reg_cc_tab
            case('normtab')
                print *, 'normalizing and soft-sorting...'
                call reg_obj%sort_tab
                call reg_obj%ref_reg_cc_tab
            case('hard')
                print *, 'cluster-hard-sorting the tab...'
                allocate(best_ir(params%fromp:params%top), best_ip(params%fromp:params%top))
                call reg_obj%cluster_sort_tab(best_ip, best_ir)
                call reg_obj%uniform_cavgs(best_ip, best_ir)
            case('unihard')
                print *, 'uniformly-hard-sorting the tab...'
                allocate(best_ir(params%fromp:params%top), best_ip(params%fromp:params%top))
                call reg_obj%uniform_sort_tab(best_ip, best_ir)
                call reg_obj%uniform_cavgs(best_ip, best_ir)
        end select
        ! descaling
        if( params%l_reg_scale ) call pftcc%reg_descale
        call reg_obj%regularize_refs
        select case(trim(params%reg_mode))
            case('tab')
            case('normtab')
                print *, 'Reconstructing the 3D volume ...'
                ! init volumes
                call preprecvols
                ! prep img, fpls, ctfparms
                allocate(fpls(params%fromp:params%top),ctfparms(params%fromp:params%top))
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sdev)
                do iptcl = params%fromp,params%top
                    if( .not.ptcl_mask(iptcl) ) cycle
                    call build%imgbatch(iptcl)%norm_noise(build%lmsk, sdev)
                    call build%imgbatch(iptcl)%fft
                    call fpls(iptcl)%new(build%imgbatch(iptcl))
                    ctfparms(iptcl) = build_glob%spproj%get_ctfparams(params%oritype, iptcl)
                    call fpls(iptcl)%gen_planes(build%imgbatch(iptcl), ctfparms(iptcl), iptcl=iptcl)
                enddo
                !$omp end parallel do
                do iref = 1, params%nspace
                    euls = build_glob%eulspace%get_euler(iref)
                    do j = 1, params%reg_num
                        pind_here = params%fromp + j - 1
                        if( reg_obj%ref_ptcl_tab(pind_here, iref)%prob < TINY ) cycle
                        iptcl = reg_obj%ref_ptcl_tab(pind_here, iref)%iptcl
                        call build_glob%spproj_field%get_ori(iptcl, orientation)
                        if( orientation%isstatezero() ) cycle
                        ! getting the particle orientation
                        shvec = orientation%get_2Dshift() + reg_obj%ref_ptcl_tab(pind_here,iref)%sh
                        call orientation%set_shift(shvec)
                        loc     = reg_obj%ref_ptcl_tab(pind_here, iref)%loc
                        euls(3) = 360. - pftcc%get_rot(loc)
                        call orientation%set_euler(euls)
                        call orientation%set('w', reg_obj%ref_ptcl_tab(pind_here, iref)%prob)
                        ! insert
                        call grid_ptcl(fpls(iptcl), build_glob%pgrpsyms, orientation)
                    enddo
                enddo
                ! normalise structure factors
                call norm_struct_facts( cline )
                ! destruct
                call killrecvols()
                call orientation%kill
            case('hard','unihard')
                print *, 'Reconstructing the 3D volume (unihard-alignment) ...'
                ! init volumes
                call preprecvols
                ! prep img, fpls, ctfparms
                allocate(fpls(params%fromp:params%top),ctfparms(params%fromp:params%top))
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sdev)
                do iptcl = params%fromp,params%top
                    if( .not.ptcl_mask(iptcl) ) cycle
                    call build%imgbatch(iptcl)%norm_noise(build%lmsk, sdev)
                    call build%imgbatch(iptcl)%fft
                    call fpls(iptcl)%new(build%imgbatch(iptcl))
                    ctfparms(iptcl) = build_glob%spproj%get_ctfparams(params%oritype, iptcl)
                    call fpls(iptcl)%gen_planes(build%imgbatch(iptcl), ctfparms(iptcl), iptcl=iptcl)
                enddo
                !$omp end parallel do
                do j = params%fromp,params%top
                    iref  = best_ir(j)
                    iptcl = best_ip(j)
                    if( reg_obj%ref_ptcl_tab(iptcl, iref)%prob < TINY ) cycle
                    euls = build_glob%eulspace%get_euler(iref)
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( orientation%isstatezero() ) cycle
                    ! getting the particle orientation
                    shvec = orientation%get_2Dshift() + reg_obj%ref_ptcl_tab(iptcl,iref)%sh
                    call orientation%set_shift(shvec)
                    loc     = reg_obj%ref_ptcl_tab(iptcl, iref)%loc
                    euls(3) = 360. - pftcc%get_rot(loc)
                    call orientation%set_euler(euls)
                    call orientation%set('w', reg_obj%ref_ptcl_tab(iptcl, iref)%prob)
                    ! insert
                    call grid_ptcl(fpls(iptcl), build_glob%pgrpsyms, orientation)
                enddo
                ! normalise structure factors
                call norm_struct_facts( cline )
                ! destruct
                call killrecvols()
                call orientation%kill
        end select
        call simple_end('**** SIMPLE_CHECK_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_align

    subroutine exec_check_align_inpl( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_strategy2D3D_common, only: read_imgbatch, prepimgbatch, prepimg4align, calcrefvolshift_and_mapshifts2ptcls,&
                    &read_and_filter_refvols, preprefvol, preprecvols, norm_struct_facts, killrecvols, grid_ptcl
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_fplane,              only: fplane
        use simple_regularizer_inpl,    only: regularizer_inpl
        use simple_image
        class(check_align_inpl_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,          parameter   :: MAXITS = 60
        integer,          allocatable :: pinds(:), best_ip(:), best_ir(:), best_irot(:)
        logical,          allocatable :: ptcl_mask(:)
        complex,          allocatable :: cmat(:,:)
        real,             allocatable :: sigma2_noise(:,:)
        type(image),      allocatable :: tmp_imgs(:)
        type(fplane),     allocatable :: fpls(:)
        type(ctfparams),  allocatable :: ctfparms(:)
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(ori)                     :: o_tmp
        type(image)                   :: img
        type(regularizer_inpl)        :: reg_inpl
        type(ori)                     :: orientation
        integer  :: nptcls, iptcl, j, s, iref, box, loc, pind_here, ithr, irot
        logical  :: l_ctf, do_center
        real     :: xyz(3), euls(3), shvec(2), sdev
        call cline%set('mkdir',    'yes')
        call cline%set('oritype',  'ptcl3D')
        call cline%set('center',   'no')
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        call build%spproj%update_projinfo(cline)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        allocate(ptcl_mask(params%fromp:params%top))
        call build_glob%spproj_field%sample4update_and_incrcnt([params%fromp,params%top],&
            &1.0, nptcls, pinds, ptcl_mask)
        print *, 'nptcls = ', nptcls, '; fromp = ', params%fromp, '; top = ', params%top
        call pftcc%new(params%nspace, [1,nptcls], params%kfromto)
        call pftcc%reallocate_ptcls(nptcls, pinds)
        call reg_inpl%new(pftcc)
        print *, 'Preparing the references ...'
        ! e/o partioning
        if( build%spproj%os_ptcl3D%get_nevenodd() == 0 )then
            call build%spproj%os_ptcl3D%partition_eo
            call build%spproj%write_segment_inside(params%oritype,params%projfile)
        endif
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        do s=1,params%nstates
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params%vols(s), do_center, xyz)
            call read_and_filter_refvols( cline, params%vols(s), params%vols(s))
            ! PREPARE E/O VOLUMES
            call preprefvol(cline, s, do_center, xyz, .false.)
            call preprefvol(cline, s, do_center, xyz, .true.)
            ! PREPARE REFERENCES
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1, params%nspace
                call build_glob%eulspace%get_ori(iref, o_tmp)
                call build_glob%vol_odd%fproject_polar((s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.false., mask=build_glob%l_resmsk)
                call build_glob%vol%fproject_polar(    (s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.true.,  mask=build_glob%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
        end do
        call pftcc%memoize_refs
        ! PREPARATION OF PARTICLES
        print *, 'Preparing the particles ...'
        call prepimgbatch(params%top-params%fromp+1)
        call read_imgbatch([params%fromp,params%top], ptcl_mask)
        call build%img_crop_polarizer%init_polarizer(pftcc, params%alpha)
        allocate(tmp_imgs(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call tmp_imgs(ithr)%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(iptcl,ithr) schedule(static) proc_bind(close)
        do iptcl = params%fromp,params%top
            if( .not.ptcl_mask(iptcl) ) cycle
            ithr = omp_get_thread_num()+1
            ! prep
            call prepimg4align(iptcl, build%imgbatch(iptcl), tmp_imgs(ithr))
            call build%imgbatch(iptcl)%ifft ! for reconstruction
            ! transfer to polar coordinates
            call build%img_crop_polarizer%polarize(pftcc, tmp_imgs(ithr), iptcl, .true., .true., mask=build%l_resmsk)
            ! e/o flags
            ! call pftcc%set_eo(iptcl, nint(build_glob%spproj_field%get(iptcl,'eo'))<=0 )
            call pftcc%set_eo(iptcl, .true.)
        enddo
        !$omp end parallel do
        ! getting the ctfs
        l_ctf = build%spproj%get_ctfflag('ptcl2D',iptcl=pinds(1)).ne.'no'
        ! make CTFs
        if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, 'ptcl2D')

        ! ALIGNMENT OF PARTICLES
        print *, 'Aligning the particles ...'
        ! using gencorrs (cc-based to estimate the sigma)
        if( params%l_needs_sigma )then
            allocate( sigma2_noise(pftcc%kfromto(1):pftcc%kfromto(2), params%fromp:params%top), source=1. )
            ! do j = pftcc%kfromto(1),pftcc%kfromto(2)
            !     sigma2_noise(j,:) = real(j)
            ! enddo
            call pftcc%assign_sigma2_noise(sigma2_noise)
            ! call pftcc%memoize_ptcls
            ! params%cc_objfun = OBJFUN_CC
            ! !$omp parallel do default(shared) private(j,iref,ithr,iptcl,inpl_corrs,cxy,max_corr,max_iref,max_sh,max_loc,loc,corr,sh) proc_bind(close) schedule(static)
            ! do j = 1, nptcls
            !     max_corr = 0.
            !     do iref = 1, params%nspace
            !         ithr  = omp_get_thread_num() + 1
            !         iptcl = pinds(j)
            !         ! find best irot/shift for this pair of iref, iptcl
            !         call pftcc%gencorrs( iref, iptcl, inpl_corrs )
            !         loc = maxloc(inpl_corrs, dim=1)
            !         call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
            !         cxy = grad_shsrch_obj(ithr)%minimize(irot=loc)
            !         if( loc > 0 )then
            !             corr = cxy(1)
            !             sh   = cxy(2:3)
            !         else
            !             loc  = maxloc(inpl_corrs, dim=1)
            !             corr = inpl_corrs(loc)
            !             sh   = 0.
            !         endif
            !         if( corr > max_corr )then
            !             max_corr = corr
            !             max_loc  = loc
            !             max_sh   = sh
            !             max_iref = iref
            !         endif
            !     enddo
            !     call pftcc%update_sigma( max_iref, iptcl, max_sh, max_loc )
            ! enddo
            ! !$omp end parallel do
            ! params%cc_objfun = OBJFUN_EUCLID
        endif
        ! actual alignment using the defined cost function
        ! scaling by the ctf
        if( params%l_reg_scale )then
            call pftcc%reg_scale
            !$omp parallel do default(shared) private(j) proc_bind(close) schedule(static)
            do j = 1, nptcls
                call pftcc%memoize_sqsum_ptcl(pinds(j))
            enddo
            !$omp end parallel do
        endif
        ! Memoize particles FFT parameters
        call pftcc%memoize_ptcls
        call reg_inpl%init_tab
        call reg_inpl%fill_tab(pinds)
        print *, 'Assembling the class averages with'
        select case(trim(params%reg_mode))
            case('tab')
                print *, 'soft-sorting the tab...'
                call reg_inpl%sort_tab_no_norm
                call reg_inpl%ref_reg_cc_tab
            case('normtab')
                print *, 'normalizing and soft-sorting...'
                call reg_inpl%sort_tab
                call reg_inpl%ref_reg_cc_tab
            case('hard')
                print *, 'cluster-hard-sorting the tab...'
                allocate(best_ir(params%fromp:params%top),&
                        &best_ip(params%fromp:params%top),&
                        &best_irot(params%fromp:params%top))
                call reg_inpl%cluster_sort_tab(best_ip, best_ir, best_irot)
                call reg_inpl%uniform_cavgs(best_ip, best_ir, best_irot)
            case('unihard')
                print *, 'uniformly-hard-sorting the tab...'
                allocate(best_ir(params%fromp:params%top),&
                        &best_ip(params%fromp:params%top),&
                        &best_irot(params%fromp:params%top))
                call reg_inpl%uniform_sort_tab(best_ip, best_ir, best_irot)
                call reg_inpl%uniform_cavgs(best_ip, best_ir, best_irot)
        end select
        ! descaling
        if( params%l_reg_scale ) call pftcc%reg_descale
        call reg_inpl%regularize_refs
        select case(trim(params%reg_mode))
            case('tab')
            case('normtab')
                print *, 'Reconstructing the 3D volume ...'
                ! init volumes
                call preprecvols
                ! prep img, fpls, ctfparms
                allocate(fpls(params%fromp:params%top),ctfparms(params%fromp:params%top))
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sdev)
                do iptcl = params%fromp,params%top
                    if( .not.ptcl_mask(iptcl) ) cycle
                    call build%imgbatch(iptcl)%norm_noise(build%lmsk, sdev)
                    call build%imgbatch(iptcl)%fft
                    call fpls(iptcl)%new(build%imgbatch(iptcl))
                    ctfparms(iptcl) = build_glob%spproj%get_ctfparams(params%oritype, iptcl)
                    call fpls(iptcl)%gen_planes(build%imgbatch(iptcl), ctfparms(iptcl), iptcl=iptcl)
                enddo
                !$omp end parallel do
                do irot = 1, reg_inpl%reg_nrots
                    do iref = 1, params%nspace
                        euls = build_glob%eulspace%get_euler(iref)
                        do j = 1, params%reg_num
                            pind_here = params%fromp + j - 1
                            if( reg_inpl%ref_ptcl_tab(pind_here, iref, irot)%prob < TINY ) cycle
                            iptcl = reg_inpl%ref_ptcl_tab(pind_here, iref, irot)%iptcl
                            call build_glob%spproj_field%get_ori(iptcl, orientation)
                            if( orientation%isstatezero() ) cycle
                            ! getting the particle orientation
                            shvec = orientation%get_2Dshift() + reg_inpl%ref_ptcl_tab(pind_here,iref,irot)%sh
                            call orientation%set_shift(shvec)
                            loc     = reg_inpl%ref_ptcl_tab(pind_here, iref, irot)%loc
                            euls(3) = 360. - pftcc%get_rot(loc)
                            call orientation%set_euler(euls)
                            call orientation%set('w', reg_inpl%ref_ptcl_tab(pind_here, iref, irot)%prob)
                            ! insert
                            call grid_ptcl(fpls(iptcl), build_glob%pgrpsyms, orientation)
                        enddo
                    enddo
                enddo
                ! normalise structure factors
                call norm_struct_facts( cline )
                ! destruct
                call killrecvols()
                call orientation%kill
            case('hard', 'unihard')
                print *, 'Reconstructing the 3D volume (unihard-alignment) ...'
                ! init volumes
                call preprecvols
                ! prep img, fpls, ctfparms
                allocate(fpls(params%fromp:params%top),ctfparms(params%fromp:params%top))
                !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sdev)
                do iptcl = params%fromp,params%top
                    if( .not.ptcl_mask(iptcl) ) cycle
                    call build%imgbatch(iptcl)%norm_noise(build%lmsk, sdev)
                    call build%imgbatch(iptcl)%fft
                    call fpls(iptcl)%new(build%imgbatch(iptcl))
                    ctfparms(iptcl) = build_glob%spproj%get_ctfparams(params%oritype, iptcl)
                    call fpls(iptcl)%gen_planes(build%imgbatch(iptcl), ctfparms(iptcl), iptcl=iptcl)
                enddo
                !$omp end parallel do
                do j = params%fromp,params%top
                    iref  = best_ir(j)
                    iptcl = best_ip(j)
                    irot  = best_irot(j)
                    if( reg_inpl%ref_ptcl_tab(iptcl, iref, irot)%prob < TINY ) cycle
                    euls = build_glob%eulspace%get_euler(iref)
                    call build_glob%spproj_field%get_ori(iptcl, orientation)
                    if( orientation%isstatezero() ) cycle
                    ! getting the particle orientation
                    shvec = orientation%get_2Dshift() + reg_inpl%ref_ptcl_tab(iptcl,iref,irot)%sh
                    call orientation%set_shift(shvec)
                    loc     = reg_inpl%ref_ptcl_tab(iptcl, iref, irot)%loc
                    euls(3) = 360. - pftcc%get_rot(loc)
                    call orientation%set_euler(euls)
                    call orientation%set('w', reg_inpl%ref_ptcl_tab(iptcl, iref, irot)%prob)
                    ! insert
                    call grid_ptcl(fpls(iptcl), build_glob%pgrpsyms, orientation)
                enddo
                ! normalise structure factors
                call norm_struct_facts( cline )
                ! destruct
                call killrecvols()
                call orientation%kill
        end select
        call simple_end('**** SIMPLE_CHECK_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_align_inpl

end module simple_commander_refine3D
