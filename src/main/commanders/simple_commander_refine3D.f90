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
use simple_decay_funs,       only: inv_cos_decay, cos_decay
use simple_image,            only: image
use simple_masker,           only: masker
use simple_exec_helpers,     only: set_master_num_threads
use simple_commander_euclid
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_auto_commander
public :: refine3D_distr_commander
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

type, extends(commander_base) :: refine3D_auto_commander
  contains
    procedure :: execute      => exec_refine3D_auto
end type refine3D_auto_commander

type, extends(commander_base) :: refine3D_distr_commander
  contains
    procedure :: execute      => exec_refine3D_distr
end type refine3D_distr_commander

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

    subroutine exec_refine3D_auto( self, cline )
        use simple_commander_rec, only: reconstruct3D_commander_distr
        use simple_sp_project,    only: sp_project
        class(refine3D_auto_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(cmdline)               :: cline_reconstruct3D_distr
        type(postprocess_commander) :: xpostprocess
        type(parameters)            :: params
        real,    parameter :: LP2SMPD_TARGET   = 1./3.
        real,    parameter :: SMPD_TARGET_MIN  = 1.3
        logical, parameter :: DEBUG            = .true.
        integer, parameter :: MINBOX           = 256
        integer, parameter :: NPDIRS4BAL       = 300
        character(len=:),   allocatable :: str_state
        real             :: smpd_target, smpd_crop, scale, trslim
        integer          :: box_crop, maxits_phase1, maxits_phase2, iter
        logical          :: l_autoscale
        ! commanders
        type(reconstruct3D_commander_distr) :: xreconstruct3D_distr
        type(refine3D_distr_commander)      :: xrefine3D_distr
        ! hard defaults
        call cline%set('balance',         'no') ! balanced particle sampling based on available 3D solution
        call cline%set('greedy_sampling', 'no') ! stochastic within-class selection without consideration to objective function value
        call cline%set('trail_rec',      'yes') ! trailing average 3D reconstruction
        call cline%set('refine',       'neigh') ! greedy multi-neighborhood 3D refinement 
        call cline%set('icm',            'yes') ! ICM regularization to maximize map connectivity
        call cline%set('automsk',        'yes') ! envelope masking for background flattening
        call cline%set('sh_first',       'yes') ! estimate shifts before rotational search
        call cline%set('overlap',         0.99) ! convergence if overlap > 99%
        call cline%set('nstates',            1) ! only single-state refinement is supported
        call cline%set('objfun',      'euclid') ! the objective function is noise-normalized Euclidean distance
        call cline%set('envfsc',         'yes') ! we use the envelope mask when calculating an FSC plot
        call cline%set('lplim_crit',     0.143) ! we use the 0.143 criterion for low-pass limitation
        call cline%set('lam_anneal',     'yes') ! we conduct deterministic annealing of the lambda regularization parameter
        call cline%set('incrreslim',      'no') ! if anything 'yes' makes it slightly worse, but no real difference right now
        ! overridable defaults
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',        'yes')
        if( .not. cline%defined('center')      ) call cline%set('center',        'no') ! 4 now, probably fine
        if( .not. cline%defined('sigma_est')   ) call cline%set('sigma_est', 'global') ! 4 now, probably fine
        if( .not. cline%defined('combine_eo')  ) call cline%set('combine_eo',    'no') ! 4 now, to allow more rapid testing
        if( .not. cline%defined('prob_inpl')   ) call cline%set('prob_inpl',    'yes') ! no difference at this stage, so prefer 'yes'
        if( .not. cline%defined('update_frac') ) call cline%set('update_frac',    0.1) ! 4 now, needs testing/different logic (nsample?)
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',       'yes') ! better map with ml_reg='yes'
        if( .not. cline%defined('lp_auto')     ) call cline%set('lp_auto',       'no') ! works, should be considered if the defaults are not satisfactory
        if( .not. cline%defined('maxits')      ) call cline%set('maxits',          20) ! ~2 passes over particles, sufficient ?
        if( .not. cline%defined('keepvol')     ) call cline%set('keepvol',       'no') ! we do not keep volumes for each iteration by deafult
        call params%new(cline)
        call cline%set('maxits_glob', params%maxits) ! needed for correct lambda annealing
        call cline%set('mkdir', 'no') ! to avoid nested directory structure
        if( mod(params%maxits,2) /= 0) THROW_HARD('Maximum number of iterations (MAXITS) need to be divisible with 2')
        maxits_phase1 = params%maxits / 2
        maxits_phase2 = maxits_phase1
        if( params%box <= MINBOX )then
            smpd_target = params%smpd
            smpd_crop   = params%smpd
            box_crop    = params%box
            scale       = 1.0
            l_autoscale = .false.
        else
            smpd_target = max(SMPD_TARGET_MIN, params%res_target * LP2SMPD_TARGET)
            call autoscale(params%box, params%smpd, smpd_target, box_crop, smpd_crop, scale, minbox=MINBOX)
            l_autoscale = box_crop < params%box
        endif
        trslim = min(8.,max(2.0, AHELIX_WIDTH / smpd_crop))
        if( DEBUG )then
            print *, 'smpd_target: ', smpd_target
            print *, 'box:         ', params%box
            print *, 'box_crop:    ', box_crop
            print *, 'smpd:        ', params%smpd
            print *, 'smpd_crop:   ', smpd_crop
            print *, 'scale:       ', scale
            print *, 'trslim:      ', trslim
            print *, 'l_autoscale: ', l_autoscale
        endif
        call cline%set('trs', trslim)
        if( l_autoscale )then
            call cline%set('box_crop',  box_crop)
            call cline%set('smpd_crop', smpd_crop)
        endif
        ! generate an initial 3D reconstruction
        cline_reconstruct3D_distr = cline
        call cline_reconstruct3D_distr%set('prg', 'reconstruct3D') ! required for distributed call
        call cline_reconstruct3D_distr%delete('trail_rec')
        call cline_reconstruct3D_distr%delete('objfun')
        call cline_reconstruct3D_distr%delete('needs_sigma')
        call cline_reconstruct3D_distr%delete('sigma_est')
        call cline_reconstruct3D_distr%delete('update_frac')
        call cline_reconstruct3D_distr%set('objfun', 'cc') ! ugly, but this is how it works in parameters
        call xreconstruct3D_distr%execute_safe(cline_reconstruct3D_distr)
        ! 3D refinement, phase1
        str_state = int2str_pad(1,2)
        call cline%set('vol1', VOL_FBODY//str_state//params_glob%ext)
        params%mskfile = MSKVOL_FILE
        call cline%set('mskfile',           MSKVOL_FILE)
        call cline%set('prg',                'refine3D')
        call cline%set('ufrac_trec', params%update_frac)
        call cline%set('maxits',          maxits_phase1)
        call cline%set('lp_auto',                 'yes')
        call xrefine3D_distr%execute_safe(cline)
        ! iteration number bookkeeping
        iter = 0
        if( cline%defined('endit') )then
            iter = cline%get_iarg('endit')
            call cline%delete('endit')
        endif
        iter = iter + 1
        ! re-reconstruct from all particle images
        call cline_reconstruct3D_distr%set('mskfile', MSKVOL_FILE)
        call xreconstruct3D_distr%execute_safe(cline_reconstruct3D_distr)
        ! 3D refinement, phase2
        call cline%set('vol1', VOL_FBODY//str_state//params_glob%ext)
        params%mskfile = MSKVOL_FILE
        call cline%set('mskfile',           MSKVOL_FILE)
        call cline%set('maxits',          maxits_phase1)
        call cline%set('lp_auto',  trim(params%lp_auto))
        call cline%set('startit',                  iter)
        call cline%set('which_iter',               iter)
        call xrefine3D_distr%execute_safe(cline)
        ! re-reconstruct from all particle images
        call xreconstruct3D_distr%execute_safe(cline_reconstruct3D_distr)
        ! postprocess
        call cline%set('prg', 'postprocess')
        call cline%set('mkdir', 'yes')
        call xpostprocess%execute_safe(cline)

        !**************TESTS2DO
        ! check so that we have a starting 3D alignment
        ! check so that all states are 0 or 1 and fall over otherwise
        ! should test 40k projection directions and see how that performs (params class needs modification)
        
    end subroutine exec_refine3D_auto

    subroutine exec_refine3D_distr( self, cline )
        use simple_commander_rec,    only: reconstruct3D_commander_distr, volassemble_commander
        use simple_fsc,              only: plot_fsc
        use simple_commander_euclid, only: calc_group_sigmas_commander
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_polarops
        class(refine3D_distr_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(reconstruct3D_commander_distr)   :: xreconstruct3D_distr
        type(calc_pspec_commander_distr)      :: xcalc_pspec_distr
        type(check_3Dconv_commander)          :: xcheck_3Dconv
        type(postprocess_commander)           :: xpostprocess
        type(refine3D_commander)              :: xrefine3D_shmem
        type(calc_group_sigmas_commander)     :: xcalc_group_sigmas
        type(estimate_first_sigmas_commander) :: xfirst_sigmas_distr
        type(prob_align_commander)            :: xprob_align_distr
        type(volassemble_commander)           :: xvolassemble
        type(polarft_corrcalc)                :: pftcc
        ! command lines
        type(cmdline) :: cline_reconstruct3D_distr
        type(cmdline) :: cline_calc_pspec_distr
        type(cmdline) :: cline_prob_align_distr
        type(cmdline) :: cline_calc_group_sigmas
        type(cmdline) :: cline_check_3Dconv
        type(cmdline) :: cline_volassemble
        type(cmdline) :: cline_postprocess
        type(cmdline) :: cline_tmp
        integer(timer_int_kind) :: t_init,   t_scheduled,  t_merge_algndocs,  t_volassemble,  t_tot
        real(timer_int_kind)    :: rt_init, rt_scheduled, rt_merge_algndocs, rt_volassemble
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        type(masker)     :: mskvol
        type(image)      :: vol_e, vol_o    
        character(len=:),          allocatable :: prev_refine_path, target_name, fname_vol, fname_even, fname_odd
        character(len=LONGSTRLEN), allocatable :: list(:)
        integer,                   allocatable :: state_pops(:)
        real,                      allocatable :: res(:), fsc(:)
        character(len=LONGSTRLEN) :: vol, vol_iter, str, str_iter, fsc_templ
        character(len=STDLEN)     :: vol_even, vol_odd, str_state, fsc_file, volpproc, vollp
        logical :: err, vol_defined, have_oris, converged, fall_over, l_multistates, l_automsk
        logical :: l_combine_eo, l_griddingset, do_automsk, l_polar
        real    :: corr, smpd
        integer :: i, state, iter, box, nfiles, niters, nthr_here
        601 format(A,1X,F12.3)
        if( .not. cline%defined('nparts') )then
            call xrefine3D_shmem%execute(cline)
            return
        endif
        ! deal with # threads for the master process
        call set_master_num_threads( nthr_here, 'REFINE3D' )
        ! local options & flags
        l_multistates = cline%defined('nstates')
        l_griddingset = cline%defined('gridding')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',      'yes')
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp',        30.)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! init
        call build%init_params_and_build_spproj(cline, params)
        l_polar = trim(params%polar).eq.'yes'
        if( l_polar )then
            call build%build_general_tbox(params, cline, do3d=.true.)
            call pftcc%new(1, [1,1], params%kfromto)
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
        ! final iteration with combined e/o
        l_combine_eo = .false.
        if( trim(params%combine_eo).eq.'yes' )then
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
        cline_prob_align_distr    = cline
        cline_check_3Dconv        = cline
        cline_postprocess         = cline
        cline_calc_group_sigmas   = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' )     ! required for distributed call
        call cline_calc_pspec_distr%set(    'prg', 'calc_pspec' )        ! required for distributed call
        call cline_prob_align_distr%set(    'prg', 'prob_align' )        ! required for distributed call
        call cline_postprocess%set(         'prg', 'postprocess' )       ! required for local call
        call cline_calc_group_sigmas%set(   'prg', 'calc_group_sigmas' ) ! required for local call
        call cline_postprocess%set('mirr',    'no')
        call cline_postprocess%set('mkdir',   'no')
        call cline_postprocess%set('imgkind', 'vol')
        if( trim(params%oritype).eq.'cls3D' ) call cline_postprocess%set('imgkind', 'vol_cavg')
        ! removes unnecessary volume keys and generates volassemble finished names
        do state = 1,params%nstates
            vol = 'vol'//int2str( state )
            call cline_check_3Dconv%delete( vol )
            call cline_postprocess%delete( vol )
        enddo
        if( trim(params%objfun).eq.'euclid' ) call cline%set('needs_sigma','yes')
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
                    call build%spproj%get_vol('vol_cavg', state, fname_vol, smpd, box)
                else
                    call build%spproj%get_vol('vol', state, fname_vol, smpd, box)
                endif
                call cline%set(trim(vol), fname_vol)
                params%vols(state) = fname_vol
            end do
            prev_refine_path = get_fpath(fname_vol)
            if( trim(simple_abspath(prev_refine_path,check_exists=.false.)) .eq. trim(cwd_glob) )then
                ! ...unless we operate in the same folder
                do state=1,params%nstates
                    str_state = int2str_pad(state,2)
                    fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                    if( .not.file_exists(fsc_file)) THROW_HARD('Missing file: '//trim(fsc_file))
                end do
                if( params%l_update_frac )then
                    call simple_list_files(prev_refine_path//'*recvol_state*part*', list)
                    nfiles = size(list)
                    err = params%nparts * 4 /= nfiles
                    if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
                    deallocate(list)
                endif
                if( trim(params%objfun).eq.'euclid' )then
                    call cline%set('needs_sigma','yes')
                    call cline_reconstruct3D_distr%set('needs_sigma','yes')
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
                ! if we are doing objfun=euclid the sigma estimates need to be carried over
                if( trim(params%objfun).eq.'euclid' )then
                    call cline%set('needs_sigma','yes')
                    call cline_reconstruct3D_distr%set('needs_sigma','yes')
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
            call xcalc_pspec_distr%execute_safe(cline_calc_pspec_distr)
            ! check if we have input volume(s) and/or 3D orientations
            vol_defined = .false.
            do state = 1,params%nstates
                vol_defined = cline%defined('vol'//int2str(state))
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
                call cline_tmp%delete('trail_rec')
                call cline_tmp%delete('objfun')
                call cline_tmp%delete('needs_sigma')
                call cline_tmp%delete('sigma_est')
                call cline_tmp%set('objfun', 'cc') ! ugly, but this is how it works in parameters 
                call xreconstruct3D_distr%execute_safe( cline_tmp )
                do state = 1,params%nstates
                    ! rename volumes and update cline
                    str_state = int2str_pad(state,2)
                    vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                    str       = trim(STARTVOL_FBODY)//trim(str_state)//params%ext
                    call      simple_rename( trim(vol), trim(str) )
                    ! update command line
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
                call cline_calc_group_sigmas%set('which_iter', params%startit)
                call cline_calc_group_sigmas%set('nthr', nthr_here)
                call xcalc_group_sigmas%execute_safe(cline_calc_group_sigmas)
                ! then, estimate first sigmas given reconstructed starting volumes(s) and previous orientations
                if( .not.cline%defined('nspace') ) call cline%set('nspace', real(params%nspace))
                if( .not.cline%defined('athres') ) call cline%set('athres', real(params%athres))
                call xfirst_sigmas_distr%execute_safe(cline)
                ! update command lines
                call cline%set('needs_sigma','yes')
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
                params%eps = inv_cos_decay(iter, params%maxits_glob, params%eps_bounds)
                write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
            endif
            if( params%l_lam_anneal )then
                params%lambda = cos_decay(iter, params%maxits_glob, params%lam_bounds)
                write(logfhandle,601) '>>> LAMBDA, MAP CONNECTIVITY ANNEALING        ', params%lambda
            endif
            if( trim(params%objfun).eq.'euclid' )then
                call cline_calc_group_sigmas%set('which_iter', iter)
                call xcalc_group_sigmas%execute_safe(cline_calc_group_sigmas)
            endif
            if( have_oris .or. iter > params%startit )then
                call build%spproj%read(params%projfile)
            endif
            if( str_has_substr(params%refine, 'prob') )then
                cline_prob_align_distr = cline
                call cline_prob_align_distr%set('which_iter', params%which_iter)
                call cline_prob_align_distr%set('startit',    iter)
                call xprob_align_distr%execute_safe( cline_prob_align_distr )
            endif
            call job_descr%set( 'which_iter', trim(int2str(params%which_iter)))
            call cline%set(     'which_iter', params%which_iter)
            call job_descr%set( 'startit',    trim(int2str(iter)))
            call cline%set(     'startit',    iter)
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
            if( L_BENCH_GLOB ) rt_merge_algndocs = toc(t_merge_algndocs)
            if( l_polar )then
                params%refs = trim(CAVGS_ITER_FBODY)//int2str_pad(params%which_iter,3)//params%ext
                call polar_cavger_new(pftcc, .true., nrefs=params%nspace)
                call polar_cavger_calc_pops(build%spproj)
                call polar_cavger_assemble_sums_from_parts(reforis=build_glob%eulspace)
                call build%clsfrcs%new(params%nspace, params%box_crop, params%smpd_crop, params%nstates)
                call polar_cavger_calc_and_write_frcs_and_eoavg(FRCS_FILE)
                call polar_cavger_writeall(get_fbody(params%refs,params%ext,separator=.false.))
                call polar_cavger_write_cartrefs(pftcc, get_fbody(params%refs,params%ext,separator=.false.), 'merged')
                call polar_cavger_kill
            endif
            ! CONVERGENCE
            converged = .false.
            select case(trim(params%refine))
                case('eval')
                    ! nothing to do
                case DEFAULT
                    call cline_check_3Dconv%set('nthr', nthr_here)
                    call xcheck_3Dconv%execute_safe(cline_check_3Dconv)
                    if( iter >= params%startit + 2 )then
                        ! after a minimum of 2 iterations
                        if( cline_check_3Dconv%get_carg('converged') .eq. 'yes' ) converged = .true.
                    endif
            end select
            if( niters == params%maxits ) converged = .true.
            if( L_BENCH_GLOB ) t_volassemble = tic()
            if( (trim(params%volrec).eq.'yes') )then
                ! ASSEMBLE VOLUMES
                select case(trim(params%refine))
                    case('eval')
                        ! nothing to do
                    case DEFAULT
                        cline_volassemble = cline ! transfer run-time command-line
                        call cline_volassemble%set('which_iter', params%which_iter)
                        call cline_volassemble%set('nthr',       nthr_here)
                        call xvolassemble%execute_safe(cline_volassemble)
                        ! rename & add volumes to project & update job_descr
                        call build%spproj_field%get_pops(state_pops, 'state')
                        do state = 1,params%nstates
                            str_state = int2str_pad(state,2)
                            if( state_pops(state) == 0 )then
                                ! cleanup for empty state
                                vol = 'vol'//trim(int2str(state))
                                call cline%delete( vol )
                                call job_descr%delete( vol )
                                if( trim(params%oritype).eq.'cls3D' )then
                                    call build%spproj%remove_entry_from_osout('vol_cavg', state)
                                else
                                    call build%spproj%remove_entry_from_osout('vol', state)
                                endif
                                call build%spproj%remove_entry_from_osout('fsc', state)
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
                            call cline_postprocess%set('state',    state)
                            call cline_postprocess%set('nthr', nthr_here)
                            if( cline%defined('lp') ) call cline_postprocess%set('lp', params%lp)
                            call xpostprocess%execute_safe(cline_postprocess)
                            volpproc = trim(VOL_FBODY)//trim(str_state)//PPROC_SUFFIX//params%ext
                            vollp    = trim(VOL_FBODY)//trim(str_state)//LP_SUFFIX//params%ext
                            if( l_automsk )then
                                do_automsk = .false.
                                if( niters == 1 .and. .not.params%l_filemsk )then
                                    do_automsk = .true.
                                else if( mod(iter,AMSK_FREQ)==0 )then
                                    do_automsk = .true.
                                endif
                                if( do_automsk )then
                                    call build%spproj%get_vol('vol', state, fname_vol, smpd, box)
                                    fname_even = add2fbody(trim(fname_vol), params%ext, '_even')
                                    fname_odd  = add2fbody(trim(fname_vol), params%ext, '_odd' )
                                    call vol_e%new([box,box,box], smpd)
                                    call vol_e%read(fname_even)
                                    call vol_o%new([box,box,box], smpd)
                                    call vol_o%read(fname_odd)
                                    if( cline%defined('thres') )then
                                        call mskvol%automask3D(vol_e, vol_o, trim(params%automsk).eq.'tight', params%thres)
                                    else
                                        call mskvol%automask3D(vol_e, vol_o, trim(params%automsk).eq.'tight')
                                    endif
                                    call mskvol%write(MSKVOL_FILE)
                                    params%mskfile   = MSKVOL_FILE
                                    params%l_filemsk = .true.
                                    call cline%set('mskfile', MSKVOL_FILE)
                                    call job_descr%set('mskfile', MSKVOL_FILE)
                                    call mskvol%kill_bimg
                                    call vol_e%kill
                                    call vol_o%kill
                                endif
                            endif
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//PPROC_SUFFIX//params%ext
                            call simple_copy_file(volpproc, vol_iter)
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//LP_SUFFIX//params%ext
                            call simple_copy_file(vollp, vol_iter)
                            if( iter > 1 .and. params%keepvol.eq.'no' )then
                                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter-1,3)//PPROC_SUFFIX//params%ext
                                call del_file(vol_iter)
                                vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter-1,3)//LP_SUFFIX//params%ext
                                call del_file(vol_iter)
                            endif
                        enddo
                end select
            endif
            if( L_BENCH_GLOB ) rt_volassemble = toc(t_volassemble)
            if ( l_combine_eo .and. converged )then
                converged            = .false.
                l_combine_eo         = .false.
                params%combine_eo    = 'yes'
                params%l_update_frac = .false.
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
                write(logfhandle,'(A)')'>>> PERFORMING FINAL ITERATION WITH COMBINED EVEN/ODD VOLUMES'
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
        call pftcc%kill
        call simple_end('**** SIMPLE_DISTR_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D_distr

    ! this routine should be kept minimal and clean of all automasking, postprocessing, etc. for performance
    subroutine exec_refine3D( self, cline )
        use simple_strategy3D_matcher, only: refine3D_exec
        class(refine3D_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        type(estimate_first_sigmas_commander) :: xfirst_sigmas
        type(calc_group_sigmas_commander)     :: xcalc_group_sigmas
        type(calc_pspec_assemble_commander)   :: xcalc_pspec_assemble
        type(calc_pspec_commander)            :: xcalc_pspec
        type(prob_align_commander)            :: xprob_align
        type(parameters)                      :: params
        type(builder)                         :: build
        type(cmdline)                         :: cline_calc_group_sigmas, cline_prob_align
        type(cmdline)                         :: cline_calc_pspec, cline_first_sigmas
        character(len=STDLEN)                 :: str_state, fsc_file, vol, vol_iter
        integer                               :: startit, i, state
        real                                  :: corr
        logical                               :: converged, l_sigma
        601 format(A,1X,F12.3)
        call build%init_params_and_build_strategy3D_tbox(cline,params)
        startit = 1
        if( cline%defined('startit') ) startit = params%startit
        select case(trim(params%refine))
            case('prob', 'prob_state')
                ! random sampling and updatecnt dealt with in prob_align
            case DEFAULT
                if( startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
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
                params%l_needs_sigma    = .true.
                cline_calc_group_sigmas = cline
                if( file_exists(trim(SIGMA2_GROUP_FBODY)//trim(int2str(params%which_iter))//trim(STAR_EXT)) )then
                    ! it is assumed that we already have precalculated sigmas2 and all corresponding flags have been set
                else
                    ! sigma2 not provided & are calculated
                    if( build%spproj_field%get_nevenodd() == 0 )then
                        ! make sure we have e/o partitioning prior to calc_pspec
                        call build%spproj_field%partition_eo
                        call build%spproj%write_segment_inside(params%oritype)
                    endif
                    if( startit == 1 )then
                        ! make sure we have weights for first_sigmas
                        call build%spproj_field%set_all2single('w', 1.0)
                        call build%spproj%write_segment_inside(params%oritype)
                    endif
                    cline_calc_pspec   = cline
                    cline_first_sigmas = cline
                    call xcalc_pspec%execute_safe( cline_calc_pspec )
                    call cline_calc_group_sigmas%set('which_iter', startit)
                    call xcalc_pspec_assemble%execute_safe(cline_calc_group_sigmas)
                    if( .not.cline_first_sigmas%defined('nspace') ) call cline_first_sigmas%set('nspace', params%nspace)
                    if( .not.cline_first_sigmas%defined('athres') ) call cline_first_sigmas%set('athres', params%athres)
                    call xfirst_sigmas%execute_safe(cline)
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
                    params%eps = inv_cos_decay(params%which_iter, params%maxits_glob, params%eps_bounds)
                    write(logfhandle,601) '>>> SNR, WHITE NOISE REGULARIZATION           ', params%eps
                endif
                if( params%l_lam_anneal )then
                    params%lambda = cos_decay(params%which_iter, params%maxits_glob, params%lam_bounds)
                    write(logfhandle,601) '>>> LAMBDA, MAP CONNECTIVITY ANNEALING        ', params%lambda
                endif
                if( l_sigma )then
                    call cline_calc_group_sigmas%set('which_iter', params%which_iter)
                    call xcalc_group_sigmas%execute(cline_calc_group_sigmas)
                endif
                if( str_has_substr(params%refine, 'prob') )then
                    cline_prob_align = cline
                    call cline_prob_align%set('prg',       'prob_align')
                    call cline_prob_align%set('which_iter', params%which_iter)
                    do state = 1, params%nstates
                        call cline%set('vol'//int2str(state), params%vols(state))
                    enddo
                    call xprob_align%execute( cline_prob_align )
                endif
                ! in strategy3D_matcher:
                call refine3D_exec(cline, params%which_iter, converged)
                ! convergence
                if( converged .or. i == params%maxits )then
                    ! report the last iteration on exit
                    call cline%delete( 'startit' )
                    call cline%set('endit', real(params%which_iter))
                    ! update project with the new orientations
                    call build%spproj%write_segment_inside(params%oritype)
                    call del_file(params%outfile)
                    if( trim(params%volrec) .eq. 'yes' )then
                        do state = 1, params%nstates
                            if( build%spproj_field%get_pop(state, 'state') == 0 )then
                                ! cleanup empty state
                                if( trim(params%oritype).eq.'cls3D' )then
                                    call build%spproj%remove_entry_from_osout('vol_cavg', state)
                                else
                                    call build%spproj%remove_entry_from_osout('vol', state)
                                endif
                                call build%spproj%remove_entry_from_osout('fsc', state)
                            else
                                ! add state volume, fsc to os_out
                                str_state = int2str_pad(state,2)
                                fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                                call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                                vol       = trim(VOL_FBODY)//trim(str_state)//params%ext
                                vol_iter  = trim(vol)
                                if( trim(params%oritype).eq.'cls3D' )then
                                    call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol_cavg')
                                else
                                    call build%spproj%add_vol2os_out(vol_iter, params%smpd_crop, state, 'vol')
                                endif
                            endif
                        end do
                        ! volume mask, one for all states
                        if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                        call build%spproj%write_segment_inside('out')
                    endif
                    if( l_sigma )then
                        ! so final sigma2 can be used for a subsequent refine3D run
                        call cline_calc_group_sigmas%set('which_iter',params%which_iter+1)
                        call xcalc_group_sigmas%execute_safe(cline_calc_group_sigmas)
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
        call cline_first_sigmas%set('maxits',       1)
        call cline_first_sigmas%set('which_iter',   1)
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
        logical           :: converged
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! check convergence
        converged = conv%check_conv3D(cline, params%msk)
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
        call cline%set('frac_srch', conv%get('frac_srch'))
        ! end gracefully
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CHECK_3DCONV NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_check_3Dconv

    subroutine exec_prob_tab( self, cline )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_strategy2D3D_common
        use simple_polarft_corrcalc,    only: polarft_corrcalc
        use simple_eul_prob_tab,        only: eul_prob_tab
        use simple_euclid_sigma2,       only: euclid_sigma2
        class(prob_tab_commander), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        type(image),      allocatable :: tmp_imgs(:)
        character(len=:), allocatable :: fname
        type(polarft_corrcalc)        :: pftcc
        type(builder)                 :: build
        type(parameters)              :: params
        type(eul_prob_tab)            :: eulprob_obj_part
        type(euclid_sigma2)           :: eucl_sigma
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( cline )
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, POLARIZER, PTCLS
        call prepare_refs_sigmas_ptcls( pftcc, cline, eucl_sigma, tmp_imgs, nptcls, params%which_iter,&
                                        do_polar=(trim(params%polar).eq.'yes' .and. params%which_iter>1) )
        ! Build polar particle images
        call build_batch_particles(pftcc, nptcls, pinds, tmp_imgs)
        ! Filling prob table in eul_prob_tab
        call eulprob_obj_part%new(pinds)
        fname = trim(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( str_has_substr(params%refine, 'prob_state') )then
            call eulprob_obj_part%fill_tab_state_only(pftcc)
            call eulprob_obj_part%write_state_tab(fname)
        else
            call eulprob_obj_part%fill_tab(pftcc)
            call eulprob_obj_part%write_tab(fname)
        endif
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
        use simple_eul_prob_tab,        only: eul_prob_tab
        use simple_strategy2D3D_common, only: sample_ptcls4fillin, sample_ptcls4update
        class(prob_align_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,            allocatable :: pinds(:)
        character(len=:),   allocatable :: fname
        type(builder)                   :: build
        type(parameters)                :: params
        type(prob_tab_commander)        :: xprob_tab
        type(eul_prob_tab)              :: eulprob_obj_glob
        type(cmdline)                   :: cline_prob_tab
        type(qsys_env)                  :: qenv
        type(chash)                     :: job_descr
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
        if( params_glob%startit == 1 ) call build_glob%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params_glob%l_fillin .and. mod(params_glob%startit,5) == 0 )then
            call sample_ptcls4fillin([1,params_glob%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update([1,params_glob%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build_glob%spproj%write_segment_inside(params_glob%oritype)
        ! more prep
        call eulprob_obj_glob%new(pinds)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute_safe(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params_glob%nparts, nptcls=params_glob%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! reading corrs from all parts
        if( str_has_substr(params%refine, 'prob_state') )then
            do ipart = 1, params_glob%nparts
                fname = trim(DIST_FBODY)//int2str_pad(ipart,params_glob%numlen)//'.dat'
                call eulprob_obj_glob%read_state_tab(fname)
            enddo
            call eulprob_obj_glob%state_assign
        else
            do ipart = 1, params_glob%nparts
                fname = trim(DIST_FBODY)//int2str_pad(ipart,params_glob%numlen)//'.dat'
                call eulprob_obj_glob%read_tab_to_glob(fname)
            enddo
            call eulprob_obj_glob%ref_assign
        endif
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
