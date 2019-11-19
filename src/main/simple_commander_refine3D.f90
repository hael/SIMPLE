! concrete commander: refine3D for ab initio 3D reconstruction and 3D refinement
module simple_commander_refine3D
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_parameters,     only: parameters
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
implicit none

public :: nspace_commander
public :: refine3D_commander_distr
public :: refine3D_commander
public :: check_3Dconv_commander
public :: calc_pspec_commander_distr
public :: calc_pspec_commander
public :: calc_pspec_assemble_commander
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
            o = oris(i)
            call o%spiral
            ares = o%find_angres()
            write(logfhandle,'(A,1X,I7,1X,A,1X,F5.2)') 'NR OF PROJDIRS:', i, 'RESOLUTION:', resang(ares, params%moldiam)
        end do
        call simple_end('**** SIMPLE_NSPACE NORMAL STOP ****')
    end subroutine exec_nspace

    subroutine exec_refine3D_distr( self, cline )
        use simple_commander_volops, only: postprocess_commander
        use simple_commander_rec,    only: reconstruct3D_commander_distr
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
        type(cmdline)    :: cline_check_3Dconv
        type(cmdline)    :: cline_volassemble
        type(cmdline)    :: cline_postprocess
        ! other variables
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        character(len=:),          allocatable :: vol_fname, prev_refine_path, target_name
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        integer,                   allocatable :: state_pops(:)
        character(len=STDLEN)     :: vol, vol_iter, str, str_iter, optlp_file
        character(len=STDLEN)     :: vol_even, vol_odd, str_state, fsc_file, volpproc
        character(len=LONGSTRLEN) :: volassemble_output
        real    :: corr, corr_prev, smpd
        integer :: ldim(3), i, state, iter, iostat, box, nfiles, niters, iter_switch2euclid, ifoo
        logical :: err, vol_defined, have_oris, do_abinitio, converged, fall_over
        logical :: l_projection_matching, l_switch2euclid, l_continue
        if( .not. cline%defined('refine') )then
            call cline%set('refine',  'single')
        else
            if( cline%get_carg('refine').eq.'multi' .and. .not. cline%defined('nstates') )then
                THROW_HARD('refine=MULTI requires specification of NSTATES')
            endif
        endif
        if( .not. cline%defined('cenlp')   ) call cline%set('cenlp', 30.)
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
        if( fall_over )then
            THROW_HARD('no particles found! :exec_refine3D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' )then
            call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        endif
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr = cline
        cline_calc_pspec_distr    = cline
        cline_check_3Dconv        = cline
        cline_volassemble         = cline
        cline_postprocess         = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' ) ! required for distributed call
        call cline_calc_pspec_distr%set( 'prg', 'calc_pspec' )       ! required for distributed call
        call cline_postprocess%set('prg', 'postprocess' )            ! required for local call
        if( trim(params%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set( 'pgrp', 'c1' )
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
        ! E/O PARTITIONING
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
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
                params%vols(state) = vol_fname
            end do
            prev_refine_path = get_fpath(vol_fname)
            ! carry over FRCs/FSCs
            ! one FSC file per state
            do state=1,params%nstates
                str_state = int2str_pad(state,2)
                fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                call simple_copy_file(trim(prev_refine_path)//trim(fsc_file), fsc_file)
            end do
            ! one FRC file for all states
            call simple_copy_file(trim(prev_refine_path)//trim(FRCS_FILE), trim(FRCS_FILE))
            ! carry over the oridistributions_part* files
            call simple_list_files(prev_refine_path//'oridistributions_part*', list)
            nfiles = size(list)
            err    = params%nparts /= nfiles
            if( err ) THROW_HARD('# partitions not consistent with previous refinement round')
            do i=1,nfiles
                target_name = PATH_HERE//basename(trim(list(i)))
                call simple_copy_file(trim(list(i)), target_name)
            end do
            deallocate(list)
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
        l_projection_matching = .false.
        if( have_oris .and. .not. vol_defined )then
            ! reconstructions needed
            call xreconstruct3D_distr%execute( cline_reconstruct3D_distr )
            do state = 1,params%nstates
                ! rename volumes and update cline
                str_state = int2str_pad(state,2)
                vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//params%ext
                iostat = simple_rename( trim(vol), trim(str) )
                vol = 'vol'//trim(int2str(state))
                call cline%set( trim(vol), trim(str) )
                vol_even = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params%ext
                iostat= simple_rename( trim(vol_even), trim(str) )
                vol_odd  = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params%ext
                iostat =  simple_rename( trim(vol_odd), trim(str) )
            enddo
        else if( vol_defined .and. params%continue .ne. 'yes' )then
            ! projection matching
            l_projection_matching = .true.
            if( .not. have_oris )then
                select case( params%neigh )
                    case( 'yes' )
                        THROW_HARD('refinement method requires input orientations')
                    case DEFAULT
                        ! all good
                end select
            endif
            if( .not.cline%defined('lp') ) THROW_HARD('LP needs be defined for the first step of projection matching!')
            call cline%delete('update_frac')
            if( params%neigh .ne. 'yes' )then
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
            if( l_projection_matching .and. cline%defined('lp_iters') ) iter_switch2euclid = params%lp_iters
            call cline%set('needs_sigma','yes')
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        niters = 0
        iter   = params%startit - 1
        corr   = -1.
        do
            niters            = niters + 1
            iter              = iter + 1
            params%which_iter = iter
            str_iter          = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
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
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            call job_descr%set('which_iter', trim(int2str(params%which_iter)))
            call cline%set('which_iter', real(params%which_iter))
            call job_descr%set( 'startit', trim(int2str(iter)))
            call cline%set('startit', real(iter))
            ! switch to refine=greedy_* when frac >= 99 and iter >= 5
            if( cline_check_3Dconv%defined('frac_srch') )then
                if( iter >= MIN_ITERS_SHC )then
                    if( cline_check_3Dconv%get_rarg('frac_srch') >= FRAC_GREEDY_LIM )then
                        select case(trim(params%refine))
                            case('single')
                                params%refine = 'greedy_single'
                            case('multi')
                                params%refine = 'greedy_multi'
                        end select
                        call job_descr%set( 'refine', params%refine )
                        call cline%set('refine', params%refine)
                        call cline_check_3Dconv%set('refine',params%refine)
                    endif
                endif
            endif
            ! FRCs
            if( cline%defined('frcs') )then
                ! all good
            else
                call job_descr%set('frcs', trim(FRCS_FILE))
            endif
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
            ! ASSEMBLE ALIGNMENT DOCS
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, params%oritype, ALGN_FBODY)
            ! ASSEMBLE VOLUMES
            select case(trim(params%refine))
            case('eval')
                ! nothing to do
            case DEFAULT
                call cline_volassemble%set( 'prg', 'volassemble' ) ! required for cmdline exec
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
                            iostat   = simple_rename( vol, vol_iter )
                        else
                            vol_iter = trim(vol)
                        endif
                        fsc_file   = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                        optlp_file = ANISOLP_FBODY//trim(str_state)//params%ext
                        ! add filters to os_out
                        call build%spproj%add_fsc2os_out(fsc_file, state, params%box)
                        call build%spproj%add_vol2os_out(optlp_file, params%smpd, state, 'vol_filt', box=params%box)
                        ! add state volume to os_out
                        if( trim(params%oritype).eq.'cls3D' )then
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol_cavg')
                        else
                            call build%spproj%add_vol2os_out(vol_iter, params%smpd, state, 'vol')
                        endif
                        ! updates cmdlines & job description
                        vol = 'vol'//trim(int2str(state))
                        call job_descr%set( vol, vol_iter )
                        call cline%set(vol, vol_iter )
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
                    if( state_pops(state) == 0 )cycle
                    call cline_postprocess%set('state', real(state))
                    if( cline%defined('lp') ) call cline_postprocess%set('lp', params%lp)
                    call xpostprocess%execute(cline_postprocess)
                    ! for gui visualization
                    if( params%refine .ne. 'snhc' )then
                        volpproc = trim(VOL_FBODY)//trim(str_state)//PPROC_SUFFIX//params%ext
                        vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter,3)//PPROC_SUFFIX//params%ext
                        call simple_copy_file(volpproc, vol_iter) ! for GUI visualization
                        if( iter > 1 )then
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//int2str_pad(iter-1,3)//PPROC_SUFFIX//params%ext
                            call del_file(vol_iter)
                        endif
                    endif
                enddo
            end select
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
            if( l_projection_matching .and. cline%defined('lp_iters') .and. (niters == params%lp_iters ) )then
                ! e/o projection matching
                write(logfhandle,'(A)')'>>>'
                write(logfhandle,'(A)')'>>> SWITCHING TO EVEN/ODD RESOLUTION LIMIT'
                write(logfhandle,'(A)')'>>>'
                l_projection_matching = .false.
                if( cline%defined('match_filt') )then
                    if( cline%get_carg('match_filt').eq.'no' )then
                        ! flags are kept so match_filt is not used
                        call job_descr%set('match_filt','no')
                    else
                        call cline%delete('lp')
                        call job_descr%delete('lp')
                        call cline_postprocess%delete('lp')
                    endif
                else
                    call cline%delete('lp')
                    call job_descr%delete('lp')
                    call cline_postprocess%delete('lp')
                endif
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
        if( cline%defined('startit') )startit = params%startit
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
        allocate( maplp(params%nstates), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In simple_commander_refine3D:: exec_check3D_conv", alloc_stat)
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
            case('cluster','clustersym','clustersoft')
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
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec_distr')
        endif
        ! init
        call build%init_params_and_build_spproj(cline, params)
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
        complex, allocatable :: cmat(:,:,:), cmat_sum(:,:,:)
        ! real,    allocatable :: pspecs(:,:), pspec(:)
        logical, allocatable :: mask(:)
        real                 :: sdev_noise
        integer              :: batchlims(2),iptcl,iptcl_batch,imatch,nyq,nptcls_part,batchsz_max
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! init
        nptcls_part = params%top-params%fromp+1
        nyq         = build%img%get_nyq()
        batchsz_max = 10 * nthr_glob
        allocate(mask(batchsz_max),source=.false.)
        ! allocate(pspecs(nyq,nptcls_part),source=0.)
        call prepimgbatch(batchsz_max)
        call sum_img%new([params%boxmatch,params%boxmatch,1],params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_sum = sum_img%get_cmat()
        cmat     = sum_img%get_cmat()
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
            cmat = cmplx(0.,0.)
            !$omp parallel do default(shared) private(iptcl,imatch,cmat)&
            !$omp schedule(static) proc_bind(close) reduction(+:cmat_sum)
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
                ! call build%imgbatch(imatch)%spectrum('power',pspec,norm=.true.)
                ! pspecs(:,iptcl) = pspec
                ! global average
                call build%imgbatch(imatch)%get_cmat_sub(cmat)
                cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmat(:,:,:)
            end do
            !$omp end parallel do
        end do
        call sum_img%set_cmat(cmat_sum)
        call sum_img%write('sum_img_part'//int2str_pad(params%part,params%numlen)//params%ext)
        !! debug
        ! call sum_img%ifft
        ! call sum_img%write('realspace_sum_img_part'//int2str(params%part)//params%ext)
        !! debug
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
        character(len=:), allocatable    :: part_fname
        integer :: iptcl,ipart,nptcls,nptcls_sel,eo,ngroups,igroup
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
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
        !! debug
        ! call avg_img%ifft
        ! call avg_img%write('avg_img'//params%ext)
        !! debug
        ! calculate power spectrum

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
        ! do iptcl = 1,nptcls
        !     if( build%spproj_field%get_state(iptcl) == 0 ) cycle
        !     eo     = nint(build%spproj_field%get(iptcl,'eo')) ! 0/1
        !     igroup = nint(build%spproj_field%get(iptcl,'stkind'))
        !
        ! enddo
        ! end gracefully
        call simple_touch('CALC_PSPEC_FINISHED',errmsg='In: commander_refine3D::calc_pspec_assemble')
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_calc_pspec_assemble

end module simple_commander_refine3D
