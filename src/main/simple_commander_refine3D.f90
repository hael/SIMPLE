! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_sigma2_binfile, only: sigma2_binfile
use simple_qsys_env,       only: qsys_env
use simple_euclid_sigma2,  only: write_groups_starfile
use simple_cluster_seed,   only: gen_labelling
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_commander_distr
public :: refine3D_commander
public :: check_3Dconv_commander
public :: calc_pspec_commander_distr
public :: calc_pspec_commander
public :: calc_pspec_assemble_commander
public :: calc_group_sigmas_commander
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
type, extends(commander_base) :: calc_pspec_commander_distr
  contains
    procedure :: execute      => exec_calc_pspec_distr
end type calc_pspec_commander_distr
type, extends(commander_base) :: calc_pspec_commander
  contains
    procedure :: execute      => exec_calc_pspec
end type calc_pspec_commander
type, extends(commander_base) :: calc_pspec_assemble_commander
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type calc_pspec_assemble_commander
type, extends(commander_base) :: calc_group_sigmas_commander
  contains
    procedure :: execute      => exec_calc_group_sigmas
end type calc_group_sigmas_commander

type :: sigma_array
    character(len=:), allocatable :: fname
    real,             allocatable :: sigma2(:,:)
end type sigma_array

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
        use simple_commander_volops, only: postprocess_commander
        use simple_commander_rec,    only: reconstruct3D_commander_distr
        use simple_estimate_ssnr,    only: plot_fsc
        class(refine3D_commander_distr), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(calc_pspec_commander_distr)    :: xcalc_pspec_distr
        type(check_3Dconv_commander)        :: xcheck_3Dconv
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline)    :: cline_reconstruct3D_distr
        type(cmdline)    :: cline_calc_pspec_distr
        type(cmdline)    :: cline_cls3D_distr
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
        logical :: l_projmatch, l_lp_iters, l_switch2euclid, l_continue, l_multistates
        real    :: corr, corr_prev, smpd, lplim
        integer :: ldim(3), i, state, iter, box, nfiles, niters, iter_switch2euclid, ifoo
        integer :: ncls, icls, ind, fnr
        l_multistates = .false.
        if( cline%defined('nstates') ) l_multistates = .true.
        l_lp_iters = cline%defined('lp_iters')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('ptclw')   ) call cline%set('ptclw',       'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! objfun=euclid logics, part 1
        l_switch2euclid  = .false.
        if( cline%defined('objfun') )then
            l_continue = .false.
            if( cline%defined('continue') ) l_continue = trim(cline%get_carg('continue')).eq.'yes'
            if( (trim(cline%get_carg('objfun')).eq.'euclid') .and. .not.l_continue )then
                l_switch2euclid = .true.
                call cline%set('objfun','cc')
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
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' ) call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr = cline
        cline_cls3D_distr         = cline
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
                call cline_cls3D_distr%set(trim(vol), vol_fname)
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
                call find_ldim_nptcls(params%vols(state),ldim,ifoo)
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
            if( .not.cline%defined('lp') )then
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
            if( l_projmatch .and. l_lp_iters ) iter_switch2euclid = params%lp_iters
            call cline%set('needs_sigma','yes')
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        niters = 0
        iter   = params%startit - 1
        corr   = -1.
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
                ! use of l_needs_sigma?
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
                call job_descr%set('frcs', trim(FRCS_FILE))
            endif
            ! schedule
            if( L_BENCH_GLOB )then
                rt_init = toc(t_init)
                t_scheduled = tic()
            endif
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
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
                            call cline_cls3D_distr%set(vol, vol_iter)
                            if( params%keepvol.ne.'no' )then
                                call simple_copy_file(vol_iter,trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//params%ext)
                            endif
                        endif
                    enddo
                    ! volume mask, one for all states
                    if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                    ! writes os_out
                    call build%spproj%write_segment_inside('out')
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
            if( l_lp_iters .and. (niters == params%lp_iters ) )then
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
                call cline%set('objfun','euclid')
                call cline%set('match_filt','no')
                call cline%delete('lp')
                call job_descr%set('objfun','euclid')
                call job_descr%set('match_filt','no')
                call job_descr%delete('lp')
                call cline_volassemble%set('objfun','euclid')
                call cline_postprocess%delete('lp')
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
        type(parameters) :: params
        type(builder)    :: build
        integer :: startit
        logical :: converged
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        if( startit == 1 ) call build%spproj_field%clean_updatecnt
        if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
        call refine3D_exec( cline, startit, converged) ! partition or not, depending on 'part'
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
                converged = conv%check_conv3D(cline, params%msk)
            endif
        else
            select case(params%refine)
            case('cluster','clustersym')
                    converged = conv%check_conv_cluster(cline)
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

    subroutine exec_calc_pspec_distr( self, cline )
        class(calc_pspec_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=STDLEN), parameter :: PSPEC_FBODY = 'pspec_'
        ! command lines
        type(cmdline)    :: cline_calc_pspec
        type(cmdline)    :: cline_calc_pspec_assemble
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        logical          :: fall_over
        if( .not. cline%defined('mkdir')    ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec_distr')
        endif
        ! init
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl2D','ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case DEFAULT
                write(logfhandle,*)'Unsupported ORITYPE; simple_commander_refine3D :: exec_refine3D_distr'
        end select
        if( fall_over )then
            THROW_HARD('no particles found! :exec_refine3D_distr')
        endif
        if( build%spproj_field%get_nevenodd() == 0 )then
            THROW_HARD('no eve/odd flag found! :calc_pspec__distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! prepare command lines from prototype master
        cline_calc_pspec          = cline
        cline_calc_pspec_assemble = cline
        ! initialise static command line parameters and static job description parameter
        call cline_calc_pspec%set('prg', 'calc_pspec' )                   ! required for distributed call
        call cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble' ) ! required for local call
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble
        call qenv%exec_simple_prg_in_queue(cline_calc_pspec_assemble, 'CALC_PSPEC_FINISHED')
        ! end gracefully
        call qsys_cleanup
        call build%spproj%kill
        call simple_end('**** SIMPLE_DISTR_CALC_PSPEC NORMAL STOP ****')
    end subroutine exec_calc_pspec_distr

    subroutine exec_calc_pspec( self, cline )
        use simple_parameters,          only: params_glob
        use simple_strategy2D3D_common, only: prepimgbatch, read_imgbatch
        use simple_image,               only: image
        class(calc_pspec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: sum_img
        type(builder)        :: build
        complex, pointer     :: cmat(:,:,:), cmat_sum(:,:,:)
        real,    allocatable :: pspecs(:,:), pspec(:)
        logical, allocatable :: mask(:)
        real                 :: sdev_noise
        integer              :: batchlims(2),iptcl,iptcl_batch,imatch,nyq,nptcls_part,batchsz_max
        real, allocatable    :: sigma2(:,:)
        type(sigma2_binfile) :: binfile
        integer              :: kfromto(2)
        character(len=:), allocatable :: binfname
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! init
        nptcls_part = params%top-params%fromp+1
        nyq         = build%img%get_nyq()
        batchsz_max = 10 * nthr_glob
        allocate(mask(batchsz_max),source=.false.)
        allocate(pspecs(nyq,nptcls_part),source=0.)
        call prepimgbatch(batchsz_max)
        call sum_img%new([params%box,params%box,1],params%smpd)
        call sum_img%zero_and_flag_ft
        call sum_img%get_cmat_ptr(cmat_sum)
        do iptcl_batch=params_glob%fromp,params_glob%top,batchsz_max
            batchlims = [iptcl_batch, min(params_glob%top,iptcl_batch + batchsz_max - 1)]
            ! mask
            do iptcl = batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                mask(imatch) = .not. (build%spproj_field%get_state(iptcl) == 0)
            enddo
            ! read
            call read_imgbatch( batchlims )
            ! preprocess
            !$omp parallel do default(shared) private(iptcl,imatch,pspec)&
            !$omp schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                if( .not. mask(imatch) ) cycle
                ! normalize
                call build%imgbatch(imatch)%noise_norm(build%lmsk, sdev_noise)
                !  mask
                if( params%l_innermsk )then
                    call build%imgbatch(imatch)%mask(params%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
                else
                    if( params%l_focusmsk )then
                        call build%imgbatch(imatch)%mask(params%focusmsk, 'soft')
                    else
                        call build%imgbatch(imatch)%mask(params%msk, 'soft')
                    endif
                endif
                call build%imgbatch(imatch)%fft
                ! power spectrum
                call build%imgbatch(imatch)%spectrum('power',pspec,norm=.true.)
                pspecs(:,iptcl-params_glob%fromp+1) = pspec
            end do
            !$omp end parallel do
            ! global average
            do iptcl=batchlims(1),batchlims(2)
                imatch = iptcl - batchlims(1) + 1
                if( .not. mask(imatch) ) cycle
                call build%imgbatch(imatch)%get_cmat_ptr(cmat)
                !$omp workshare
                cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmat(:,:,:)
                !$omp end workshare
            enddo
        end do
        call sum_img%write('sum_img_part'//int2str_pad(params%part,params%numlen)//params%ext)
        ! write to disk
        kfromto(1) = 1
        kfromto(2) = nyq
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        allocate(sigma2(nyq,params%fromp:params%top))
        do iptcl = params%fromp, params%top
            sigma2(:,iptcl) = pspecs(:,iptcl-params%fromp+1)
        end do
        call binfile%new(binfname,params%fromp,params%top,kfromto)
        call binfile%write(sigma2)
        ! end gracefully
        call qsys_job_finished('simple_commander_refine3D :: exec_calc_pspec')
        call simple_end('**** SIMPLE_CALC_PSPEC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec

    subroutine exec_calc_pspec_assemble( self, cline )
        use simple_image,               only: image
        class(calc_pspec_assemble_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=STDLEN), parameter :: PSPEC_FBODY = 'pspec_'
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_binfile)             :: binfile
        type(sigma_array), allocatable   :: sigma2_arrays(:)
        character(len=:),  allocatable   :: part_fname,starfile_fname,outbin_fname
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nyq,pspec_l,pspec_u
        real,              allocatable   :: group_pspecs(:,:,:),pspec_ave(:),pspecs(:,:),sigma2_output(:,:)
        integer,           allocatable   :: group_weights(:,:)
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! set Fourier index range
        params%kfromto(1) = max(2, calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) =        calc_fourier_index(2.*params%smpd, params%box, params%smpd)
        ! generate average power spectrum
        nptcls     = build%spproj_field%get_noris(consider_state=.false.)
        nptcls_sel = build%spproj_field%get_noris(consider_state=.true.)
        call avg_img%new([params%box,params%box,1], params%smpd)
        call avg_img%zero_and_flag_ft
        do ipart = 1,params%nparts
            call build%img%zero_and_flag_ft
            part_fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext
            call build%img%read(part_fname)
            call avg_img%add(build%img)
            call del_file(part_fname)
        enddo
        call avg_img%div(real(nptcls_sel))
        ! calculate power spectrum
        call avg_img%spectrum('power',pspec_ave,norm=.true.)
        nyq = avg_img%get_nyq()
        ! read power spectra of particles
        allocate(pspecs(nyq,params%nptcls),sigma2_arrays(params%nparts))
        do ipart = 1,params%nparts
            sigma2_arrays(ipart)%fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new_from_file(sigma2_arrays(ipart)%fname)
            call binfile%read(sigma2_arrays(ipart)%sigma2)
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                THROW_HARD('commander_refine3d; exec_calc_pspec_assemble; file ' // sigma2_arrays(ipart)%fname // ' has ptcl range ' // int2str(pspec_l) // '-' // int2str(pspec_u))
            end if
            pspecs(:,pspec_l:pspec_u) = sigma2_arrays(ipart)%sigma2(:,:)
        end do
        ! generate group averages & write
        ngroups = 0
        !$omp parallel do default(shared) private(iptcl,igroup)&
        !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
            ngroups = max(igroup,ngroups)
        enddo
        !$omp end parallel do
        allocate(group_pspecs(2,ngroups,nyq),source=0.)
        allocate(group_weights(2,ngroups),source=0)
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
            igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            group_pspecs(eo+1,igroup,:) = group_pspecs(eo+1,igroup,:) + pspecs(:, iptcl)
            group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)  + 1
        enddo
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < 1 ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / real(group_weights(eo,igroup))
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) - pspec_ave(:)
                call remove_negative_sigmas(eo, igroup)
            end do
        end do
        ! write group sigmas to starfile
        starfile_fname = 'sigma2_it_1.star'
        call write_groups_starfile(starfile_fname, group_pspecs, ngroups)
        ! update sigmas in binfiles to match averages
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
            igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            pspecs(:,iptcl) = group_pspecs(eo+1,igroup,:)
        enddo
        ! write updated sigmas to disc
        do ipart = 1,params%nparts
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( allocated(sigma2_output) ) deallocate(sigma2_output)
            allocate(sigma2_output(params%kfromto(1):params%kfromto(2),pspec_l:pspec_u))
            do iptcl = pspec_l, pspec_u
                sigma2_output(params%kfromto(1):params%kfromto(2),iptcl) = pspecs(params%kfromto(1):params%kfromto(2),iptcl)
            end do
            outbin_fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call binfile%new(outbin_fname, fromp=pspec_l, top=pspec_u, kfromto=(/params%kfromto(1), params%kfromto(2)/))
            call binfile%write(sigma2_output)
        end do
        ! end gracefully
        do ipart = 1,params%nparts
            deallocate(sigma2_arrays(ipart)%fname)
            deallocate(sigma2_arrays(ipart)%sigma2)
        end do
        call simple_touch('CALC_PSPEC_FINISHED',errmsg='In: commander_refine3D::calc_pspec_assemble')
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    contains

    subroutine remove_negative_sigmas(eo, igroup)
        integer, intent(in) :: eo, igroup
        logical :: is_positive
        logical :: fixed_from_prev
        integer :: nn, idx
        ! remove any negative sigma2 noise values: replace by positive neighboring value
        do idx = 1, size(group_pspecs, 3)
            if( group_pspecs(eo,igroup,idx) < 0. )then
                ! first try the previous value
                fixed_from_prev = .false.
                if( idx - 1 >= 1 )then
                    if( group_pspecs(eo,igroup,idx-1) > 0. )then
                        group_pspecs(eo,igroup,idx) = group_pspecs(eo,igroup,idx-1)
                        fixed_from_prev = .true.
                    end if
                end if
                if( .not. fixed_from_prev )then
                    is_positive = .false.
                    nn          = idx
                    do while (.not. is_positive)
                        nn = nn + 1
                        if( nn > size(group_pspecs,3) )then
                            THROW_HARD('BUG! Cannot find positive values in sigma2 noise spectrum; eo=' // trim(int2str(eo)) // ', igroup=' // trim(int2str(igroup)))
                        end if
                        if( group_pspecs(eo,igroup,nn) > 0. )then
                            is_positive = .true.
                            group_pspecs(eo,igroup,idx) = group_pspecs(eo,igroup,nn)
                        end if
                    end do
                end if
            end if
        end do
    end subroutine remove_negative_sigmas

    end subroutine exec_calc_pspec_assemble

    subroutine exec_calc_group_sigmas( self, cline )
        use simple_image,               only: image
        class(calc_group_sigmas_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=STDLEN), parameter :: PSPEC_FBODY = 'pspec_'
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_binfile)             :: binfile
        type(sigma_array), allocatable   :: sigma2_arrays(:)
        character(len=:),  allocatable   :: part_fname,starfile_fname,outbin_fname
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup,nyq,pspec_l,pspec_u
        real,              allocatable   :: group_pspecs(:,:,:),pspec_ave(:),pspecs(:,:),sigma2_output(:,:)
        integer,           allocatable   :: group_weights(:,:)
        call cline%set('mkdir', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! set Fourier index range
        params%kfromto(1) = max(2, calc_fourier_index(params%hp, params%box, params%smpd))
        params%kfromto(2) =        calc_fourier_index(2.*params%smpd, params%box, params%smpd)
        ! read sigmas from binfiles
        allocate(pspecs(params%kfromto(1):params%kfromto(2),params%nptcls),sigma2_arrays(params%nparts))
        do ipart = 1,params%nparts
            sigma2_arrays(ipart)%fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call binfile%new_from_file(sigma2_arrays(ipart)%fname)
            call binfile%read(sigma2_arrays(ipart)%sigma2)
            pspec_l = lbound(sigma2_arrays(ipart)%sigma2,2)
            pspec_u = ubound(sigma2_arrays(ipart)%sigma2,2)
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                THROW_HARD('commander_refine3d; exec_calc_group_sigmas; file ' // sigma2_arrays(ipart)%fname // ' has ptcl range ' // int2str(pspec_l) // '-' // int2str(pspec_u))
            end if
            pspecs(:,pspec_l:pspec_u) = sigma2_arrays(ipart)%sigma2(:,:)
        end do
        ngroups = 0
        !$omp parallel do default(shared) private(iptcl,igroup)&
        !$omp schedule(static) proc_bind(close) reduction(max:ngroups)
        do iptcl = 1,params%nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            igroup  = nint(build%spproj_field%get(iptcl,'stkind'))
            ngroups = max(igroup,ngroups)
        enddo
        !$omp end parallel do
        allocate(group_pspecs(2,ngroups,params%kfromto(1):params%kfromto(2)),source=0.)
        allocate(group_weights(2,ngroups),source=0)
        do iptcl = 1,params%nptcls
            if( build%spproj_field%get_state(iptcl) == 0 ) cycle
            eo     = nint(build%spproj_field%get(iptcl,'eo'    )) ! 0/1
            igroup = nint(build%spproj_field%get(iptcl,'stkind'))
            group_pspecs(eo+1,igroup,:) = group_pspecs (eo+1,igroup,:) + pspecs(:, iptcl)
            group_weights(eo+1,igroup)  = group_weights(eo+1,igroup)   + 1
        enddo
        do eo = 1,2
            do igroup = 1,ngroups
                if( group_weights(eo,igroup) < 1 ) cycle
                group_pspecs(eo,igroup,:) = group_pspecs(eo,igroup,:) / real(group_weights(eo,igroup))
            end do
        end do
        ! write group sigmas to starfile
        starfile_fname = 'sigma2_it_' // trim(int2str(params%which_iter)) // '.star'
        call write_groups_starfile(starfile_fname, group_pspecs, ngroups)
        call simple_touch('CALC_GROUP_SIGMAS_FINISHED',errmsg='In: commander_refine3D::calc_group_sigmas')
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CALC_GROUP_SIGMAS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_group_sigmas

end module simple_commander_refine3D
