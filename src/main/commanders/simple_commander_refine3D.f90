! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_builder,          only: builder, build_glob
use simple_cmdline,          only: cmdline
use simple_commander_base,   only: commander_base
use simple_parameters,       only: parameters, params_glob
use simple_sigma2_binfile,   only: sigma2_binfile
use simple_qsys_env,         only: qsys_env
use simple_cluster_seed,     only: gen_labelling
use simple_commander_volops, only: postprocess_commander
use simple_commander_mask,   only: automask_commander
use simple_decay_funs,       only: inv_cos_decay
use simple_commander_euclid
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_commander_distr
public :: refine3D_commander
public :: estimate_first_sigmas_commander
public :: check_3Dconv_commander
public :: prob_tab_commander
public :: prob_align_commander
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

type, extends(commander_base) :: estimate_first_sigmas_commander
  contains
    procedure :: execute      => exec_estimate_first_sigmas
end type estimate_first_sigmas_commander

type, extends(commander_base) :: check_3Dconv_commander
  contains
    procedure :: execute      => exec_check_3Dconv
end type check_3Dconv_commander

type, extends(commander_base) :: prob_tab_commander
  contains
    procedure :: execute      => exec_prob_tab
end type prob_tab_commander

type, extends(commander_base) :: prob_align_commander
  contains
    procedure :: execute      => exec_prob_align
end type prob_align_commander

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
        type(prob_align_commander)            :: xprob_align
        type(automask_commander)              :: xautomask
        ! command lines
        type(cmdline) :: cline_reconstruct3D_distr
        type(cmdline) :: cline_calc_pspec_distr
        type(cmdline) :: cline_calc_sigma
        type(cmdline) :: cline_check_3Dconv
        type(cmdline) :: cline_volassemble
        type(cmdline) :: cline_postprocess
        type(cmdline) :: cline_prob_align
        type(cmdline) :: cline_tmp
        type(cmdline) :: cline_automask
        integer(timer_int_kind) :: t_init,   t_scheduled,  t_merge_algndocs,  t_volassemble,  t_tot
        real(timer_int_kind)    :: rt_init, rt_scheduled, rt_merge_algndocs, rt_volassemble, rt_tot
        character(len=STDLEN)   :: benchfname
        ! other variables
        type(parameters)    :: params
        type(builder)       :: build
        type(qsys_env)      :: qenv
        type(chash)         :: job_descr
        character(len=:),          allocatable :: vol_fname, prev_refine_path, target_name, fname_automasked
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        integer,                   allocatable :: state_pops(:)
        real,                      allocatable :: res(:), fsc(:)
        character(len=LONGSTRLEN) :: vol, vol_iter, str, str_iter, fsc_templ
        character(len=STDLEN)     :: vol_even, vol_odd, str_state, fsc_file, volpproc, vollp
        character(len=LONGSTRLEN) :: volassemble_output
        logical :: err, vol_defined, have_oris, converged, fall_over
        logical :: l_continue, l_multistates, l_automsk
        logical :: l_combine_eo, l_griddingset, do_automsk
        real    :: corr, corr_prev, smpd
        integer :: ldim(3), i, state, iter, box, nfiles, niters, ifoo
        integer :: fnr
        601 format(A,1X,F12.3)
        if( .not. cline%defined('nparts') )then
            call xrefine3D_shmem%execute(cline)
            return
        endif
        l_multistates = cline%defined('nstates')
        l_griddingset = cline%defined('gridding')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
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
        ! final iteration with combined e/o
        l_combine_eo = .false.
        if( cline%defined('combine_eo') )then
            l_combine_eo = .true.
            call cline%set('combine_eo','no')
            params%combine_eo = 'no'
        endif
        ! automasking
        l_automsk = trim(params%automsk).ne.'no'
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
        cline_prob_align          = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' )     ! required for distributed call
        call cline_calc_pspec_distr%set(    'prg', 'calc_pspec' )        ! required for distributed call
        call cline_prob_align%set(          'prg', 'prob_align' )        ! required for distributed call
        call cline_postprocess%set(         'prg', 'postprocess' )       ! required for local call
        call cline_calc_sigma%set(          'prg', 'calc_group_sigmas' ) ! required for local call
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
        if( trim(params%objfun).eq.'euclid' ) call cline_volassemble%set('needs_sigma','yes')
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
        if( params%continue .eq. 'yes' )then
            ! we are continuing from a previous refinement round,
            ! i.e. projfile is fetched from a X_refine3D dir
            ! set starting volume(s), iteration number & previous refinement path...
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
            if( trim(simple_abspath(prev_refine_path,check_exists=.false.)) .eq. trim(cwd_glob) )then
                ! ...unless we operate in the same folder
                do state=1,params%nstates
                    str_state = int2str_pad(state,2)
                    fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                    if( .not.file_exists(fsc_file)) THROW_HARD('Missing file: '//trim(fsc_file))
                end do
                if( params%l_frac_update )then
                    call simple_list_files(prev_refine_path//'*recvol_state*part*', list)
                    nfiles = size(list)
                    err = params%nparts * 4 /= nfiles
                    if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
                    deallocate(list)
                endif
                if( trim(params%objfun).eq.'euclid' )then
                    call cline%set('needs_sigma','yes')
                    call cline_reconstruct3D_distr%set('needs_sigma','yes')
                    call cline_volassemble%set('needs_sigma','yes')
                    if( .not.l_griddingset ) call cline%set('gridding','yes')
                    call simple_list_files(prev_refine_path//trim(SIGMA2_FBODY)//'*', list)
                    nfiles = size(list)
                    if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                    deallocate(list)
                endif
            else
                ! carry over FSCs
                ! one FSC file per state
                do state=1,params%nstates
                    str_state = int2str_pad(state,2)
                    fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                    call simple_copy_file(trim(prev_refine_path)//trim(fsc_file), fsc_file)
                end do
                ! if we are doing moving average volume update, partial reconstructions need to be carried over
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
        else
            ! generate initial noise power estimates
            call build%spproj_field%set_all2single('w', 1.0)
            call build%spproj%write_segment_inside(params%oritype)
            call xcalc_pspec_distr%execute(cline_calc_pspec_distr)
            ! check if we have input volume(s) and/or 3D orientations
            vol_defined = .false.
            do state = 1,params%nstates
                vol = 'vol' // int2str(state)
                if( cline%defined(trim(vol)) )then
                    vol_defined = .true.
                    !!!!!!!!!!!!!!!!!!!! cropping done in strategy2D3D_common :: read_and_filter_refvvols
                    ! call find_ldim_nptcls(trim(params%vols(state)),ldim,ifoo)
                    ! if( (ldim(1) /= params%box_crop) .and.  (ldim(1) /= params%box) )then
                    !     THROW_HARD('Incompatible dimensions between input volume and images: '//params%vols(state))
                    ! endif
                    !!!!!!!!!!!!!!!!!!!! cropping done in strategy2D3D_common :: read_and_filter_refvvols
                endif
            enddo
            have_oris = .not. build%spproj%is_virgin_field(params%oritype)
            if( .not. have_oris )then
                call build%spproj_field%rnd_oris
                have_oris = .true.
                call build%spproj%write_segment_inside(params%oritype)
            endif
            if( .not. vol_defined )then
                ! reconstructions needed
                cline_tmp = cline_reconstruct3D_distr
                call cline_tmp%delete('objfun')
                call cline_tmp%delete('needs_sigma')
                call cline_tmp%delete('sigma_est')
                call cline_tmp%set('objfun', 'cc') ! ugly, but this is how it works in parameters 
                call xreconstruct3D_distr%execute( cline_tmp )
                do state = 1,params%nstates
                    ! rename volumes and update cline
                    str_state = int2str_pad(state,2)
                    vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//params%ext
                    call      simple_rename( trim(vol), trim(str) )
                    params%vols(state) = trim(str)
                    vol       = 'vol'//trim(int2str(state))
                    call      cline%set( trim(vol), trim(str) )
                    vol_even  = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//'_even_unfil'//params%ext
                    call      simple_copy_file( trim(vol_even), trim(str) )
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params%ext
                    call      simple_rename( trim(vol_even), trim(str) )
                    vol_odd   = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//'_odd_unfil'//params%ext
                    call      simple_copy_file( trim(vol_odd), trim(str) )
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params%ext
                    call      simple_rename( trim(vol_odd), trim(str) )
                enddo
                vol_defined = .true.
            endif
            ! at this stage, we have both volume and 3D orientations (either random or previously estimated)
            if( trim(params%objfun).eq.'euclid' )then
                ! first, estimate group sigmas
                call cline_calc_sigma%set('which_iter', params%startit)
                call qenv%exec_simple_prg_in_queue(cline_calc_sigma, 'CALC_GROUP_SIGMAS_FINISHED')
                ! then, estimate first sigmas given reconstructed starting volumes(s) and previous orientations
                if( .not.cline%defined('nspace') ) call cline%set('nspace', real(params%nspace))
                if( .not.cline%defined('athres') ) call cline%set('athres', real(params%athres))
                call xfirst_sigmas%execute(cline)
                ! update command lines
                call cline%set('needs_sigma','yes')
                call cline_volassemble%set('needs_sigma','yes')
                call cline_reconstruct3D_distr%set('needs_sigma','yes')
            endif
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
            if( params%l_noise_reg )then
                ! set annealing parameter
                params%eps = inv_cos_decay(iter, params%maxits_glob, params%eps_bounds)
                write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
            endif
            if( trim(params%objfun).eq.'euclid' )then
                call cline_calc_sigma%set('which_iter', iter)
                call qenv%exec_simple_prg_in_queue(cline_calc_sigma, 'CALC_GROUP_SIGMAS_FINISHED')
            endif
            if( have_oris .or. iter > params%startit )then
                call build%spproj%read(params%projfile)
            endif
            if( str_has_substr(params%refine, 'prob') )then
                ! generate all corrs
                call cline_prob_align%set('which_iter', params%which_iter)
                call cline_prob_align%set('startit',    iter)
                do state = 1,params%nstates
                    vol = 'vol'//trim(int2str(state))
                    call cline_prob_align%set(vol, cline%get_carg(vol))
                enddo
                if( cline%defined('lp') ) call cline_prob_align%set('lp',params%lp)
                ! reading corrs from all parts into one table
                call xprob_align%execute_shmem( cline_prob_align )
            endif
            call job_descr%set( 'which_iter', trim(int2str(params%which_iter)))
            call cline%set(     'which_iter', params%which_iter)
            call job_descr%set( 'startit',    trim(int2str(iter)))
            call cline%set(     'startit',    iter)
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
                    call cline_volassemble%set( 'which_iter', params%which_iter)
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
                        volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
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
                            vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                            vol_iter  = trim(vol)
                            fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                            call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                            ! generate FSC pdf
                            res       = get_resarr(params%box_crop, params%smpd_crop)
                            fsc       = file2rarr(fsc_file)
                            fsc_templ = 'fsc_state'//trim(str_state)//'_iter'//trim(str_iter)
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
                    ! per state post-process
                    do state = 1,params%nstates
                        str_state = int2str_pad(state,2)
                        if( state_pops(state) == 0 ) cycle
                        call cline_postprocess%set('state', real(state))
                        if( cline%defined('lp') ) call cline_postprocess%set('lp', params%lp)
                        call xpostprocess%execute(cline_postprocess)         
                        volpproc = trim(VOL_FBODY)//trim(str_state)//PPROC_SUFFIX//params%ext
                        vollp    = trim(VOL_FBODY)//trim(str_state)//LP_SUFFIX//params%ext
                        if( l_automsk )then
                            if( niters == 1 .and. .not.params%l_filemsk )then
                                do_automsk = .true.
                            else if( mod(iter,AMSK_FREQ)==0 )then
                                do_automsk = .true.
                            else 
                                do_automsk = .false. 
                            endif
                            if( do_automsk )then
                                cline_automask = cline
                                call cline_automask%set('vol1', trim(vollp))
                                call cline%set('smpd', params%smpd_crop)
                                call xautomask%execute(cline_automask)
                                params%mskfile   = 'automask'//params%ext
                                params%l_filemsk = .true.
                                call cline%set('mskfile', trim(params%mskfile))
                                call job_descr%set('mskfile', trim(params%mskfile))
                                fname_automasked = basename(add2fbody(trim(vollp), params%ext, '_automsk'))
                                call del_file(fname_automasked)
                            endif
                        endif
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
        use simple_image,              only: image
        class(refine3D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(estimate_first_sigmas_commander) :: xfirst_sigmas
        type(calc_group_sigmas_commander)     :: xcalc_group_sigmas
        type(calc_pspec_assemble_commander)   :: xcalc_pspec_assemble
        type(calc_pspec_commander)            :: xcalc_pspec
        type(prob_align_commander)            :: xprob_align
        type(parameters)                      :: params
        type(builder)                         :: build
        type(cmdline)                         :: cline_calc_sigma, cline_prob_align
        type(cmdline)                         :: cline_calc_pspec, cline_first_sigmas
        type(image)                           :: noisevol
        character(len=STDLEN)                 :: str_state, fsc_file, vol, vol_iter
        integer                               :: startit, i, state, s
        real                                  :: corr, corr_prev
        logical                               :: converged, l_sigma
        601 format(A,1X,F12.3)
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        select case(trim(params%refine))
            case('prob')
                ! random sampling and updatecnt dealt with in prob_align
            case DEFAULT
                if( startit == 1 ) call build%spproj_field%clean_updatecnt_sampled
        end select
        if( params%l_distr_exec )then
            if( .not. cline%defined('outfile') ) THROW_HARD('need unique output file for parallel jobs')
            call refine3D_exec(cline, startit, converged)
        else
            if( trim(params%continue) == 'yes'    ) THROW_HARD('shared-memory implementation of refine3D does not support continue=yes')
            if( .not. file_exists(params%vols(1)) ) then
                THROW_HARD('shared-memory implementation of refine3D needs starting volume input')
            endif
            if( build%spproj%is_virgin_field(params%oritype) )then  ! we don't have orientations, so randomize
                call build%spproj_field%rnd_oris
            endif
            ! objfun=euclid
            l_sigma = .false.
            if( trim(params%objfun) == 'euclid' )then
                l_sigma = .true.
                call cline%set('needs_sigma','yes')
                params%l_needs_sigma = .true.
                cline_calc_sigma     = cline
                if( file_exists(trim(SIGMA2_GROUP_FBODY)//trim(int2str(params%which_iter))//'.star') )then
                    ! it is assumed that we already have precalculted sigmas2 and all corresponding flags have been set
                else
                    ! sigma2 not provided & are calculated
                    if( build%spproj_field%get_nevenodd() == 0 )then
                        ! make sure we have e/o partitioning prior to calc_pspec
                        call build%spproj_field%partition_eo
                        call build%spproj%write_segment_inside(params%oritype)
                    endif
                    cline_calc_pspec   = cline
                    cline_first_sigmas = cline
                    call xcalc_pspec%execute_shmem( cline_calc_pspec )
                    call cline_calc_sigma%set('which_iter', startit)
                    call xcalc_pspec_assemble%execute_shmem(cline_calc_sigma)
                    if( .not.cline_first_sigmas%defined('nspace') ) call cline_first_sigmas%set('nspace', real(params%nspace))
                    if( .not.cline_first_sigmas%defined('athres') ) call cline_first_sigmas%set('athres', real(params%athres))
                    call xfirst_sigmas%execute_shmem(cline)
                endif
            endif
            params%startit    = startit
            params%which_iter = params%startit
            params%outfile    = 'algndoc'//METADATA_EXT
            corr              = -1.
            do i = 1, params%maxits
                write(logfhandle,'(A)')   '>>>'
                write(logfhandle,'(A,I6)')'>>> ITERATION ', params%which_iter
                write(logfhandle,'(A)')   '>>>'
                if( params%l_noise_reg )then
                    ! set annealing parameter
                    params%eps = inv_cos_decay(params%which_iter, params%maxits_glob, params%eps_bounds)
                    write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
                endif
                if( l_sigma )then
                    call cline_calc_sigma%set('which_iter', params%which_iter)
                    call xcalc_group_sigmas%execute(cline_calc_sigma)
                endif
                if( str_has_substr(params%refine, 'prob') )then
                    cline_prob_align = cline
                    call cline_prob_align%set('prg',       'prob_align')
                    call cline_prob_align%set('which_iter', params%which_iter)
                    call cline_prob_align%set('vol1',       params%vols(1))
                    call xprob_align%execute( cline_prob_align )
                endif
                ! in strategy3D_matcher:
                call refine3D_exec(cline, params%which_iter, converged)
                if( converged .or. i == params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', real(params%which_iter))
                    ! update project with the new orientations
                    call build%spproj%write_segment_inside(params%oritype)
                    call del_file(params%outfile)
                    do state = 1, params%nstates
                        str_state = int2str_pad(state,2)
                        fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                        call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                        ! add state volume to os_out
                        vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                        vol_iter  = trim(vol)
                        if( trim(params%oritype).eq.'cls3D' )then
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                        else
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                        endif! volume mask, one for all states
                    end do
                    if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                    call build%spproj%write_segment_inside('out')
                    if( l_sigma )then
                        ! so final sigma2 can be used for a subsequent refine3D run
                        call cline_calc_sigma%set('which_iter',params%which_iter+1)
                        call xcalc_group_sigmas%execute_shmem(cline_calc_sigma)
                    endif
                    exit
                endif
                ! update iteration counter
                params%which_iter = params%which_iter + 1
            end do
        endif
        ! end gracefully
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        if( .not.params%l_distr_exec ) call simple_touch(JOB_FINISHED_FBODY)
        call simple_end('**** SIMPLE_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D

    subroutine exec_estimate_first_sigmas( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(estimate_first_sigmas_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        ! command lines
        type(cmdline)    :: cline_first_sigmas
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        logical          :: l_shmem, converged
        if( .not. cline%defined('vol1')     ) THROW_HARD('starting volume is needed for first sigma estimation')
        if( .not. cline%defined('pgrp')     ) THROW_HARD('point-group symmetry (pgrp) is needed for first sigma estimation')
        if( .not. cline%defined('mskdiam')  ) THROW_HARD('mask diameter (mskdiam) is needed for first sigma estimation')
        if( .not. cline%defined('nthr')     ) THROW_HARD('number of threads (nthr) is needed for first sigma estimation')
        if( .not. cline%defined('projfile') ) THROW_HARD('missing project file entry; exec_estimate_first_sigmas')
        if( .not. cline%defined('oritype')  ) call cline%set('oritype', 'ptcl3D')
        l_shmem = .not.cline%defined('nparts')
        cline_first_sigmas = cline
        call cline_first_sigmas%set('prg', 'refine3D')
        call cline_first_sigmas%set('center',    'no')
        call cline_first_sigmas%set('continue',  'no')
        call cline_first_sigmas%set('maxits',     1.0)
        call cline_first_sigmas%set('which_iter', 1.0)
        call cline_first_sigmas%set('objfun','euclid')
        call cline_first_sigmas%set('refine', 'sigma')
        call cline_first_sigmas%delete('update_frac') ! all particles neeed to contribute
        call cline_first_sigmas%delete('hp')
        call cline_first_sigmas%set('mkdir', 'no')    ! generate the sigma files in the root refine3D dir
        ! init
        if( l_shmem )then
            call build%init_params_and_build_strategy3D_tbox(cline_first_sigmas, params )
            call refine3D_exec(cline_first_sigmas, params%which_iter, converged)
            call build%kill_strategy3D_tbox
        else
            call build%init_params_and_build_spproj(cline_first_sigmas, params)
            call build%spproj%update_projinfo(cline_first_sigmas)
            call build%spproj%write_segment_inside('projinfo')
            ! setup the environment for distributed execution
            call qenv%new(params%nparts)
            ! prepare job description
            call cline_first_sigmas%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
            ! end gracefully
            call qsys_cleanup
        endif
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_ESTIMATE_FIRST_SIGMAS NORMAL STOP ****')
    end subroutine exec_estimate_first_sigmas

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
                maplp(istate) = calc_lowpass_lim(get_find_at_corr(build%fsc(istate,:),params%lplim_crit), params_glob%box_crop, params_glob%smpd_crop)
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
            converged = conv%check_conv3D(cline, params%msk)
        else
            converged = conv%check_conv3D(cline, params%msk)
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

    subroutine exec_prob_tab( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_strategy2D3D_common, only: prepimgbatch, prepimg4align, calcrefvolshift_and_mapshifts2ptcls,killimgbatch,&
                                             &read_and_filter_refvols, preprefvol, discrete_read_imgbatch, set_bp_range
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_eul_prob_tab,        only: eul_prob_tab
        use simple_euclid_sigma2,       only: euclid_sigma2
        use simple_image
        class(prob_tab_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        logical,          allocatable :: ptcl_mask(:)
        type(image),      allocatable :: tmp_imgs(:)
        character(len=:), allocatable :: fname
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(ori)                     :: o_tmp
        type(eul_prob_tab)            :: eulprob_obj_part
        type(euclid_sigma2)           :: eucl_sigma
        integer  :: nptcls, iptcl, s, ithr, iref, i
        logical  :: l_ctf, do_center
        real     :: xyz(3)
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        allocate(ptcl_mask(params%fromp:params%top))        
        if( params%l_frac_update )then
            if( build%spproj_field%has_been_sampled() )then
                call build%spproj_field%sample4update_reprod([params%fromp,params%top],&
                &nptcls, pinds, ptcl_mask)
            else
                call build%spproj_field%sample4update_rnd([params%fromp,params%top],&
                &params%update_frac, nptcls, pinds, ptcl_mask, .false.) ! no increment of sampled
            endif
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top],&
            &nptcls, pinds, ptcl_mask, .false.) ! no increement of sampled
        endif
        ! more prep
        call set_bp_range( cline )
        call pftcc%new(params%nspace * params%nstates, [1,nptcls], params%kfromto)
        call eulprob_obj_part%new(pinds)
        call prepimgbatch(nptcls)
        call discrete_read_imgbatch( nptcls, pinds, [1,nptcls] )
        call pftcc%reallocate_ptcls(nptcls, pinds)
        call build%img_crop_polarizer%init_polarizer(pftcc, params%alpha)
        allocate(tmp_imgs(nthr_glob))
        !$omp parallel do default(shared) private(ithr) schedule(static) proc_bind(close)
        do ithr = 1,nthr_glob
            call tmp_imgs(ithr)%new([params%box_crop,params%box_crop,1], params%smpd_crop, wthreads=.false.)
        enddo
        !$omp end parallel do
        ! PREPARATION OF REFERENCES IN PFTCC
        if( params%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call eucl_sigma%new(fname, params%box)
            call eucl_sigma%read_part(  build%spproj_field, ptcl_mask)
            call eucl_sigma%read_groups(build%spproj_field, ptcl_mask)
        end if
        ! read reference volumes and create polar projections
        do s=1,params%nstates
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params%vols(s), do_center, xyz, map_shift=.true.)
            call read_and_filter_refvols(s)
            ! PREPARE E/O VOLUMES
            call preprefvol(cline, s, do_center, xyz, .false.)
            call preprefvol(cline, s, do_center, xyz, .true.)
            ! PREPARE REFERENCES
            !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
            do iref=1, params%nspace
                call build%eulspace%get_ori(iref, o_tmp)
                call build%vol_odd%fproject_polar((s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.false., mask=build%l_resmsk)
                call build%vol%fproject_polar(    (s - 1) * params%nspace + iref,&
                    &o_tmp, pftcc, iseven=.true.,  mask=build%l_resmsk)
                call o_tmp%kill
            end do
            !$omp end parallel do
        end do
        call pftcc%memoize_refs
        ! PREPARATION OF PARTICLES
        !$omp parallel do default(shared) private(i,iptcl,ithr) schedule(static) proc_bind(close)
        do i = 1,nptcls
            ithr  = omp_get_thread_num()+1
            iptcl = pinds(i)
            if( .not.ptcl_mask(iptcl) ) cycle
            ! prep
            call prepimg4align(iptcl, build%imgbatch(i), tmp_imgs(ithr))
            ! transfer to polar coordinates
            call build%img_crop_polarizer%polarize(pftcc, tmp_imgs(ithr), iptcl, .true., .true., mask=build%l_resmsk)
            ! e/o flags
            call pftcc%set_eo(iptcl, mod(iptcl, 2)==0)
        enddo
        !$omp end parallel do
        ! getting the ctfs
        l_ctf = build%spproj%get_ctfflag(params%oritype,iptcl=pinds(1)).ne.'no'
        ! make CTFs
        if( l_ctf ) call pftcc%create_polar_absctfmats(build%spproj, params%oritype)
        call pftcc%memoize_ptcls
        call eulprob_obj_part%fill_tab(pftcc)
        fname = trim(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        call killimgbatch
        call pftcc%kill
        call build%kill_general_tbox
        call qsys_job_finished('simple_commander_refine3D :: exec_prob_tab')
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_align( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_eul_prob_tab, only: eul_prob_tab
        use simple_image
        class(prob_align_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        logical,          allocatable :: ptcl_mask(:)
        character(len=:), allocatable :: fname
        type(builder)                 :: build
        type(parameters)              :: params
        type(prob_tab_commander)      :: xprob_tab
        type(eul_prob_tab)            :: eulprob_obj_glob
        type(cmdline)                 :: cline_prob_tab
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        integer :: nptcls, ipart
        if( associated(build_glob) )then
            if( .not.associated(params_glob) )then
                THROW_HARD('Builder & parameters must be associated for shared memory execution!')
            endif
        else
            call cline%set('mkdir',  'no')
            call cline%set('stream', 'no')
            call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        endif
        allocate(ptcl_mask(1:params_glob%nptcls))
        if( params_glob%startit == 1 ) call build_glob%spproj_field%clean_updatecnt_sampled
        if( params_glob%l_frac_update )then
            if( params_glob%l_stoch_update )then
                call build_glob%spproj_field%sample4update_rnd([1,params_glob%nptcls],&
                &params_glob%update_frac, nptcls, pinds, ptcl_mask, .true.) ! sampled incremented
            else
                call build_glob%spproj_field%sample4update_rnd2([1,params_glob%nptcls],&
                &params_glob%update_frac, nptcls, pinds, ptcl_mask, .true.) ! sampled incremented
            endif
        else                                                    ! we sample all state > 0
            call build_glob%spproj_field%sample4update_all([1,params_glob%nptcls],&
            &nptcls, pinds, ptcl_mask, .true.) ! sampled incremented
        endif
        ! increment update counter
        call build_glob%spproj_field%incr_updatecnt([1,params_glob%nptcls], ptcl_mask)
        ! communicate to project file
        call build_glob%spproj%write_segment_inside(params_glob%oritype)        
        ! more prep
        call eulprob_obj_glob%new(pinds)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute_shmem(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params_glob%nparts, nptcls=params_glob%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! reading corrs from all parts
        do ipart = 1, params_glob%nparts
            fname = trim(DIST_FBODY)//int2str_pad(ipart,params_glob%numlen)//'.dat'
            call eulprob_obj_glob%read_tab_to_glob(fname)
        enddo
        call eulprob_obj_glob%prob_assign
        ! write the iptcl->(iref,istate) assignment
        fname = trim(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call qsys_job_finished('simple_commander_refine3D :: exec_prob_align')
        call qsys_cleanup
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

end module simple_commander_refine3D
