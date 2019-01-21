! concrete commander: distributed workflows
module simple_commander_distr_wflows
include 'simple_lib.f08'
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup, qsys_watcher
use simple_commander_base, only: commander_base
use simple_sp_project,     only: sp_project
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_builder,        only: builder
implicit none

public :: preprocess_distr_commander
public :: motion_correct_distr_commander
public :: gen_pspecs_and_thumbs_distr_commander
public :: motion_correct_tomo_distr_commander
public :: ctf_estimate_distr_commander
public :: pick_distr_commander
public :: make_cavgs_distr_commander
public :: cluster2D_distr_commander
public :: refine3D_init_distr_commander
public :: refine3D_distr_commander
public :: reconstruct3D_distr_commander
public :: tseries_track_distr_commander
public :: scale_project_distr_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_distr_commander
  contains
    procedure :: execute      => exec_preprocess_distr
end type preprocess_distr_commander
type, extends(commander_base) :: motion_correct_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_distr
end type motion_correct_distr_commander
type, extends(commander_base) :: gen_pspecs_and_thumbs_distr_commander
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs_distr
end type gen_pspecs_and_thumbs_distr_commander
type, extends(commander_base) :: motion_correct_tomo_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_tomo_distr
end type motion_correct_tomo_distr_commander
type, extends(commander_base) :: ctf_estimate_distr_commander
  contains
    procedure :: execute      => exec_ctf_estimate_distr
end type ctf_estimate_distr_commander
type, extends(commander_base) :: pick_distr_commander
  contains
    procedure :: execute      => exec_pick_distr
end type pick_distr_commander
type, extends(commander_base) :: make_cavgs_distr_commander
  contains
    procedure :: execute      => exec_make_cavgs_distr
end type make_cavgs_distr_commander
type, extends(commander_base) :: cluster2D_distr_commander
  contains
    procedure :: execute      => exec_cluster2D_distr
end type cluster2D_distr_commander
type, extends(commander_base) :: refine3D_init_distr_commander
  contains
    procedure :: execute      => exec_refine3D_init_distr
end type refine3D_init_distr_commander
type, extends(commander_base) :: refine3D_distr_commander
  contains
    procedure :: execute      => exec_refine3D_distr
end type refine3D_distr_commander
type, extends(commander_base) :: reconstruct3D_distr_commander
  contains
    procedure :: execute      => exec_reconstruct3D_distr
end type reconstruct3D_distr_commander
type, extends(commander_base) :: tseries_track_distr_commander
  contains
    procedure :: execute      => exec_tseries_track_distr
end type tseries_track_distr_commander
type, extends(commander_base) :: scale_project_distr_commander
  contains
    procedure :: execute      => exec_scale_project_distr
end type scale_project_distr_commander

contains

    subroutine exec_preprocess_distr( self, cline )
        use simple_commander_preprocess, only: preprocess_commander
        class(preprocess_distr_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(sp_project)              :: spproj
        logical                       :: l_pick
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! picking
        if( cline%defined('refs') )then
            l_pick = .true.
        else
            l_pick = .false.
        endif
        ! read in movies
        call spproj%read(params%projfile)
        ! DISTRIBUTED EXECUTION
        params%nptcls = spproj%get_nmovies()
        if( params%nptcls == 0 )then
            THROW_HARD('no movie to process! exec_preprocess_distr')
        endif
        if( params%nparts > params%nptcls ) THROW_HARD('# partitions (nparts) must be < number of entries in filetable')
        ! deal with numlen so that length matches JOB_FINISHED indicator files
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess_distr

    subroutine exec_motion_correct_distr( self, cline )
        class(motion_correct_distr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        call cline%set('oritype', 'mic')
        call params%new(cline)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nmovies() ==0 )then
            THROW_HARD('no movie to process! exec_motion_correct_distr')
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct_distr

    subroutine exec_gen_pspecs_and_thumbs_distr( self, cline )
        class(gen_pspecs_and_thumbs_distr_commander), intent(inout) :: self
        class(cmdline),                               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: nintgs
        call cline%set('oritype', 'mic')
        call params%new(cline)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nintgs = spproj%get_nintgs()
        if( nintgs ==0 )then
            THROW_HARD('no integrated movies to process! exec_gen_pspecs_and_thumbs_distr')
        endif
        if( params%nparts > nintgs )then
            call cline%set('nparts', real(nintgs))
            params%nparts = nintgs
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs_distr

    subroutine exec_motion_correct_tomo_distr( self, cline )
        use simple_oris, only: oris
        class(motion_correct_tomo_distr_commander), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: tomonames(:)
        type(parameters)         :: params
        type(oris)               :: exp_doc
        integer                  :: nseries, ipart
        type(qsys_env)           :: qenv
        character(len=KEYLEN)    :: str
        type(chash)              :: job_descr
        type(chash), allocatable :: part_params(:)
        call cline%set('prg', 'motion_correct')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('tomoseries') )then
            call read_filetable(params%tomoseries, tomonames)
        else
            THROW_HARD('need tomoseries (filetable of filetables) to be part of the command line when tomo=yes')
        endif
        nseries = size(tomonames)
        call exp_doc%new(nseries)
        if( cline%defined('exp_doc') )then
            if( file_exists(params%exp_doc) )then
                call exp_doc%read(params%exp_doc)
            else
                THROW_HARD('the required parameter file (flag exp_doc): '//trim(params%exp_doc)//' not in cwd')
            endif
        else
            THROW_HARD('need exp_doc (line: exp_time=X dose_rate=Y) to be part of the command line when tomo=yes')
        endif
        params%nparts = nseries
        params%nptcls = nseries
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts), stat=alloc_stat) ! -1. is default excluded value
        if(alloc_stat.ne.0)call allocchk("simple_commander_distr_wflows::motion_correct_tomo_moview_distr ", alloc_stat)
        do ipart=1,params%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('filetab', trim(tomonames(ipart)))
            call part_params(ipart)%set('fbody', 'tomo'//int2str_pad(ipart,params%numlen_tomo))
            str = real2str(exp_doc%get(ipart,'exp_time'))
            call part_params(ipart)%set('exp_time', trim(str))
            str = real2str(exp_doc%get(ipart,'dose_rate'))
            call part_params(ipart)%set('dose_rate', trim(str))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT_TOMO NORMAL STOP ****')
    end subroutine exec_motion_correct_tomo_distr

    subroutine exec_ctf_estimate_distr( self, cline )
        class(ctf_estimate_distr_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(chash)                   :: job_descr
        type(qsys_env)                :: qenv
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('no micrograph to process! exec_ctf_estimate_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call qsys_cleanup
        ! graceful ending
        call simple_end('**** SIMPLE_DISTR_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate_distr

    subroutine exec_pick_distr( self, cline )
        class(pick_distr_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('No micrograph to process! exec_pick_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call qsys_cleanup
        ! graceful exit
        call simple_end('**** SIMPLE_DISTR_PICK NORMAL STOP ****')
    end subroutine exec_pick_distr

    subroutine exec_make_cavgs_distr( self, cline )
        class(make_cavgs_distr_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(cmdline)    :: cline_cavgassemble
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure the use of all resources in assembly
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble class averages
        call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_cluster2D_distr( self, cline )
        use simple_procimgfile,         only: random_selection_from_imgfile, copy_imgfile
        use simple_commander_cluster2D, only: check_2Dconv_commander
        class(cluster2D_distr_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        ! commanders
        type(check_2Dconv_commander)     :: xcheck_2Dconv
        type(make_cavgs_distr_commander) :: xmake_cavgs
        ! command lines
        type(cmdline) :: cline_check_2Dconv
        type(cmdline) :: cline_cavgassemble
        type(cmdline) :: cline_make_cavgs
        ! other variables
        type(parameters)          :: params
        type(builder)             :: build
        type(qsys_env)            :: qenv
        character(len=LONGSTRLEN) :: refs, refs_even, refs_odd, str, str_iter
        integer                   :: iter
        type(chash)               :: job_descr
        real                      :: frac_srch_space
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! streaming gets its own logics because it is an exception to rules in parameters
        call cline%set('stream','no')
        ! builder & params
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        if( build%spproj%get_nptcls() == 0 )then
            THROW_HARD('no particles found! exec_cluster2D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! splitting
        call build%spproj%split_stk(params%nparts)
        ! prepare command lines from prototype master
        cline_check_2Dconv = cline
        cline_cavgassemble = cline
        cline_make_cavgs   = cline ! ncls is transferred here
        ! initialise static command line parameters and static job description parameters
        call cline_cavgassemble%set('prg', 'cavgassemble')
        call cline_cavgassemble%set('nthr', 0.) ! to ensure use of all resources in assembly
        call cline_make_cavgs%set('prg',   'make_cavgs')
        ! execute initialiser
        if( .not. cline%defined('refs') )then
            refs             = 'start2Drefs' // params%ext
            params%refs      = trim(refs)
            params%refs_even = 'start2Drefs_even'//params%ext
            params%refs_odd  = 'start2Drefs_odd'//params%ext
            if( build%spproj%is_virgin_field('ptcl2D') .or. params%startit == 1 )then
                call random_selection_from_imgfile(build%spproj, params%refs, params%box, params%ncls)
                call copy_imgfile(trim(params%refs), trim(params%refs_even), params%smpd, [1,params%ncls])
                call copy_imgfile(trim(params%refs), trim(params%refs_odd),  params%smpd, [1,params%ncls])
            else
                call cline_make_cavgs%set('refs', params%refs)
                call xmake_cavgs%execute(cline_make_cavgs)
            endif
        else
            refs = trim(params%refs)
        endif
        ! variable neighbourhood size
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! deal with eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 )then
            if( params%tseries .eq. 'yes' )then
                call build%spproj_field%partition_eo(tseries=.true.)
            else
                call build%spproj_field%partition_eo
            endif
            call build%spproj%write_segment_inside(params%oritype)
        endif
        ! main loop
        iter = params%startit - 1
        do
            iter = iter + 1
            str_iter = int2str_pad(iter,3)
            write(logfhandle,'(A)')   '>>>'
            write(logfhandle,'(A,I6)')'>>> ITERATION ', iter
            write(logfhandle,'(A)')   '>>>'
            ! cooling of the randomization rate
            params%extr_iter = params%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(params%extr_iter)))
            call cline%set('extr_iter', real(params%extr_iter))
            ! updates
            call job_descr%set('refs', trim(refs))
            call job_descr%set('startit', int2str(iter))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', trim(FRCS_FILE))
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
            ! merge orientation documents
            call build%spproj%merge_algndocs(params%nptcls, params%nparts, 'ptcl2D', ALGN_FBODY)
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // params%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // params%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // params%ext
            call cline_cavgassemble%set('refs', trim(refs))
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE_FINISHED')
            ! check convergence
            call xcheck_2Dconv%execute(cline_check_2Dconv)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check_2Dconv%get_rarg('frac')
            ! the below activates shifting & automasking
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check_2Dconv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check_2Dconv%get_rarg('trs'))
                    call job_descr%set('trs', trim(str) )
                endif
            endif
            if( cline_check_2Dconv%get_carg('converged').eq.'yes' .or. iter==params%maxits )then
                if( cline_check_2Dconv%get_carg('converged').eq.'yes' )call cline%set('converged','yes')
                exit
            endif
        end do
        call qsys_cleanup
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_cluster2D_distr

    subroutine exec_refine3D_init_distr( self, cline )
        class(refine3D_init_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)      :: params
        type(builder)         :: build
        type(cmdline)         :: cline_volassemble
        type(qsys_env)        :: qenv
        character(len=STDLEN) :: vol
        type(chash)           :: job_descr
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_spproj(cline, params)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! init
        if( cline%defined('vol1') )then
            vol = trim(params%vols(1))
        else
            vol = 'startvol_state01'//params%ext
        endif
        ! splitting
        call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! prepare command lines from prototype master
        cline_volassemble = cline
        call cline_volassemble%set( 'outvol',  vol)
        if( params%l_eo )then
            call cline_volassemble%set( 'prg', 'volassemble_eo')
        else
            call cline_volassemble%set( 'eo', 'no')
            call cline_volassemble%set( 'prg', 'volassemble')
        endif
        call cline_volassemble%set('nthr', 0.) ! to ensure use of all resources in assembly
        call qenv%gen_scripts_and_schedule_jobs( job_descr)
        call qenv%exec_simple_prg_in_queue(cline_volassemble, 'VOLASSEMBLE_FINISHED')
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_REFINE3D_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_refine3D_init_distr

    subroutine exec_refine3D_distr( self, cline )
        use simple_commander_refine3D, only: check_3Dconv_commander
        use simple_commander_volops,   only: postprocess_commander
        class(refine3D_distr_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        ! commanders
        type(reconstruct3D_distr_commander) :: xreconstruct3D_distr
        type(check_3Dconv_commander)        :: xcheck_3Dconv
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline)    :: cline_reconstruct3D_distr
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
        character(len=STDLEN)     :: vol, vol_even, vol_odd, vol_iter, vol_iter_even
        character(len=STDLEN)     :: vol_iter_odd, str, str_iter, optlp_file
        character(len=STDLEN)     :: str_state, fsc_file
        character(len=LONGSTRLEN) :: volassemble_output
        real    :: corr, corr_prev, smpd
        integer :: i, state, iter, iostat, box, nfiles
        logical :: err, vol_defined, have_oris, do_abinitio, converged, fall_over
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            write(logfhandle,*)'Unsupported ORITYPE; simple_commander_distr_wflows::exec_refine3D_distr'
        end select
        if( fall_over )then
            THROW_HARD('no particles found! :exec_refine3D_distr')
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! splitting
        call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr = cline
        cline_check_3Dconv        = cline
        cline_volassemble         = cline
        cline_postprocess         = cline
        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' ) ! required for distributed call
        call cline_postprocess%set('prg', 'postprocess' )            ! required for local call
        if( trim(params%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set( 'pgrp', 'c1' )
        call cline_postprocess%set('mirr',    'no')
        call cline_postprocess%set('mkdir',   'no')
        call cline_postprocess%set('imgkind','vol')
        call cline_volassemble%set('nthr', 0.)  ! to ensure use of all resources in assembly
        ! for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates) , stat=alloc_stat)
        if(alloc_stat /= 0)call allocchk("simple_commander_distr_wflows::exec_refine3D_distr state_assemble ",alloc_stat)
        ! removes unnecessary volume keys and generates volassemble finished names
        do state = 1,params%nstates
            vol = 'vol'//int2str( state )
            call cline_check_3Dconv%delete( trim(vol) )
            call cline_volassemble%delete( trim(vol) )
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        DebugPrint ' In exec_refine3D_distr; begin starting models'
        ! GENERATE STARTING MODELS & ORIENTATIONS
        if( params%continue .eq. 'yes' )then
            ! we are continuing from a previous refinement round,
            ! i.e. projfile is fetched from a X_refine3D dir
            do state=1,params%nstates
                vol = 'vol' // int2str(state)
                call build%spproj%get_vol('vol', state, vol_fname, smpd, box)
                call cline%set(trim(vol), vol_fname)
                params%vols(state) = vol_fname
                if( state == 1 )then
                    ! get the iteration number
                    iter = fname2iter(basename(vol_fname))
                    ! startit becomes the next
                    params%startit = iter + 1
                    call cline%set('startit', real(params%startit))
                endif
            end do
            prev_refine_path = get_fpath(vol_fname)
            ! if we are doing fractional volume update, partial reconstructions need to be carried over
            if( params%l_frac_update )then
                call simple_list_files(prev_refine_path//'*recvol_state*part*', list)
                nfiles = size(list)
                err    = .false.
                select case(trim(params%eo))
                    case('no')
                        if( params%nparts * 2 /= nfiles ) err = .true.
                    case DEFAULT
                        if( params%nparts * 4 /= nfiles ) err = .true.
                end select
                if( err )then
                    THROW_HARD('# partitions not consistent with previous refinement round')
                endif
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
            endif
            ! if we are doing objfun=euclid the sigm estimates need to be carried over
            if( trim(params%objfun) .eq. 'euclid' )then
                call simple_list_files(prev_refine_path//'sigma2_noise_part*', list)
                nfiles = size(list)
                if( nfiles /= params%nparts ) THROW_HARD('# partitions not consistent with previous refinement round')
                do i=1,nfiles
                    target_name = PATH_HERE//basename(trim(list(i)))
                    call simple_copy_file(trim(list(i)), target_name)
                end do
            endif
        endif
        vol_defined = .false.
        do state = 1,params%nstates
            vol = 'vol' // int2str(state)
            if( cline%defined(trim(vol)) ) vol_defined = .true.
        enddo
        have_oris   = .not. build%spproj%is_virgin_field(params%oritype)
        do_abinitio = .not. have_oris .and. .not. vol_defined
        if( do_abinitio )then
            call build%spproj_field%rnd_oris
            call build%spproj_field%zero_shifts
            ! take care of E/O partitioning
            if( params%l_eo )then
                if( build%spproj_field%get_nevenodd() == 0 )then
                    call build%spproj_field%partition_eo
                endif
            endif
            have_oris = .true.
            call build%spproj%write_segment_inside(params%oritype)
        endif
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
                if( params%eo .ne. 'no' )then
                    vol_even = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                    str = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//params%ext
                    iostat= simple_rename( trim(vol_even), trim(str) )
                    vol_odd  = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                    str = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//params%ext
                    iostat =  simple_rename( trim(vol_odd), trim(str) )
                endif
            enddo
        else if( .not. have_oris .and. vol_defined )then
            ! projection matching
            select case( params%neigh )
                case( 'yes' )
                    THROW_HARD('refinement method requires input orientations')
                case DEFAULT
                    ! all good
            end select
        endif
        ! EXTREMAL DYNAMICS
        if( cline%defined('extr_iter') )then
            params%extr_iter = params%extr_iter - 1
        else
            params%extr_iter = params%startit - 1
        endif
        ! EO PARTITIONING
        DebugPrint ' In exec_refine3D_distr; begin partition_eo'
        if( params%eo .ne. 'no' )then
            if( build%spproj_field%get_nevenodd() == 0 )then
                if( params%tseries .eq. 'yes' )then
                    call build%spproj_field%partition_eo(tseries=.true.)
                else
                    call build%spproj_field%partition_eo
                endif
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        iter = params%startit - 1
        corr = -1.
        do
            iter = iter + 1
            params%which_iter = iter
            str_iter = int2str_pad(iter,3)
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
            ! FRCs
            if( cline%defined('frcs') )then
                ! all good
            else
                call job_descr%set('frcs', trim(FRCS_FBODY)//'01'//BIN_EXT)
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
                if( params%eo.ne.'no' )then
                    call cline_volassemble%set( 'prg', 'volassemble_eo' ) ! required for cmdline exec
                else
                    call cline_volassemble%set( 'prg', 'volassemble' )    ! required for cmdline exec
                endif
                do state = 1,params%nstates
                    str_state = int2str_pad(state,2)
                    if( params%eo .ne. 'no' ) volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
                    call cline_volassemble%set( 'state', real(state) )
                    if( params%nstates>1 )call cline_volassemble%set('part', real(state))
                    if( params%eo .ne. 'no' )then
                        call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                        &'simple_script_state'//trim(str_state), volassemble_output)
                    else
                        call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                        &'simple_script_state'//trim(str_state))
                    endif
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
                        call job_descr%delete( trim(vol) )
                    else
                        ! rename state volume
                        vol = trim(VOL_FBODY)//trim(str_state)//params%ext
                        if( params%refine .eq. 'snhc' )then
                            vol_iter = trim(SNHCVOL)//trim(str_state)//params%ext
                        else
                            vol_iter = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//params%ext
                        endif
                        iostat = simple_rename( trim(vol), trim(vol_iter) )
                        if( params%eo .ne. 'no' )then
                            vol_even      = trim(VOL_FBODY)//trim(str_state)//'_even'//params%ext
                            vol_odd       = trim(VOL_FBODY)//trim(str_state)//'_odd' //params%ext
                            vol_iter_even = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_even'//params%ext
                            vol_iter_odd  = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_odd' //params%ext
                            iostat        = simple_rename( trim(vol_even), trim(vol_iter_even) )
                            iostat        = simple_rename( trim(vol_odd),  trim(vol_iter_odd)  )
                            fsc_file      = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                            optlp_file    = ANISOLP_FBODY//trim(str_state)//params%ext
                            ! add filters to os_out
                            call build%spproj%add_fsc2os_out(trim(fsc_file), state, params%box)
                            call build%spproj%add_vol2os_out(trim(optlp_file), params%smpd, state, 'vol_filt', box=params%box)
                        endif
                        ! add state volume to os_out
                        call build%spproj%add_vol2os_out(trim(vol_iter), params%smpd, state, 'vol')
                        ! updates cmdlines & job description
                        vol = 'vol'//trim(int2str(state))
                        call job_descr%set( trim(vol), trim(vol_iter) )
                        call cline%set( trim(vol), trim(vol_iter) )
                    endif
                enddo
                ! volume mask, one for all states
                if( cline%defined('mskfile') )call build%spproj%add_vol2os_out(trim(params%mskfile), params%smpd, 1, 'vol_msk')
                ! writes os_out
                call build%spproj%write_segment_inside('out')
                ! per state post-process
                do state = 1,params%nstates
                    if( state_pops(state) == 0 )cycle
                    call cline_postprocess%set('state', real(state))
                    call cline_postprocess%set('lp', params%lp)
                    if( params%eo .ne. 'no' )call cline_postprocess%delete('lp')
                    call xpostprocess%execute(cline_postprocess)
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
                ! activates shift search if frac >= 90
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
        call simple_end('**** SIMPLE_DISTR_REFINE3D NORMAL STOP ****')
    end subroutine exec_refine3D_distr

    subroutine exec_reconstruct3D_distr( self, cline )
        class(reconstruct3D_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=:),          allocatable :: target_name
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        character(len=LONGSTRLEN) :: refine_path
        character(len=STDLEN)     :: volassemble_output, str_state, fsc_file, optlp_file
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        integer          :: state, ipart
        logical          :: fall_over
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over )then
            THROW_HARD('No images found!')
        endif
        ! soft reconstruction from o_peaks in dir_refine?
        if( params%l_rec_soft )then
            call make_relativepath(CWD_GLOB,params%dir_refine,refine_path)
            call simple_list_files(refine_path//'oridistributions_part*.bin', list)
            if( size(list) /= params%nparts )then
                THROW_HARD('# partitions not consistent with that in '//trim(params%dir_refine))
            endif
            ! copy the orientation peak distributions
            do ipart=1,params%nparts
                target_name = PATH_HERE//basename(trim(list(ipart)))
                call simple_copy_file(trim(list(ipart)), target_name)
            end do
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! splitting
        call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! eo partitioning
        if( params%eo .ne. 'no' )then
            if( build%spproj_field%get_nevenodd() == 0 )then
                if( params%tseries .eq. 'yes' )then
                    call build%spproj_field%partition_eo(tseries=.true.)
                else
                    call build%spproj_field%partition_eo
                endif
                call build%spproj%write_segment_inside(params%oritype)
            endif
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr)
        ! assemble volumes
        ! this is for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates) )
        do state = 1, params%nstates
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        cline_volassemble = cline
        if( params%eo .ne. 'no' )then
            call cline_volassemble%set('prg', 'volassemble_eo')
        else
            call cline_volassemble%set('prg', 'volassemble')
        endif
        call cline_volassemble%set('nthr', 0.) ! to ensure the use of all resources in assembly
        ! parallel assembly
        do state = 1,params%nstates
            str_state = int2str_pad(state,2)
            if( params%eo .ne. 'no' ) volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
            call cline_volassemble%set( 'state', real(state) )
            if( params%nstates>1 )call cline_volassemble%set('part', real(state))
            if( params%eo .ne. 'no' )then
                call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                'simple_script_state'//trim(str_state), trim(volassemble_output))
            else
                call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
                'simple_script_state'//trim(str_state))
            endif
        end do
        call qsys_watcher(state_assemble_finished)
        ! updates project file only if called from another workflow
        if( params%mkdir.eq.'yes' )then
            do state = 1,params%nstates
                if( params%eo .ne. 'no' )then
                    fsc_file      = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                    optlp_file    = ANISOLP_FBODY//trim(str_state)//params%ext
                    call build%spproj%add_fsc2os_out(trim(fsc_file), state, params%box)
                    call build%spproj%add_vol2os_out(trim(optlp_file), params%smpd, state, 'vol_filt', box=params%box)
                endif
                call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd, state, 'vol')
            enddo
            call build%spproj%write_segment_inside('out',params%projfile)
        endif
        ! termination
        call qsys_cleanup
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D_distr

    subroutine exec_tseries_track_distr( self, cline )
        use simple_nrtxtfile,         only: nrtxtfile
        class(tseries_track_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, alloc_stat, j, orig_box, ipart
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( .not. file_exists(params%boxfile) ) THROW_HARD('inputted boxfile does not exist in cwd')
        if( nlines(params%boxfile) > 0 )then
            call boxfile%new(params%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk('In: simple_commander_tseries :: exec_tseries_track', alloc_stat)
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    THROW_HARD('Only square windows allowed!')
                endif
            end do
        else
            THROW_HARD('inputted boxfile is empty; exec_tseries_track')
        endif
        call boxfile%kill
        call cline%delete('boxfile')
        params%nptcls = ndatlines
        params%nparts = params%nptcls
        if( params%ncunits > params%nparts )&
        &THROW_HARD('# computational units (ncunits) mjust be <= number of entries in boxfiles')
        ! box and numlen need to be part of command line
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(3)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('ind',    int2str(ipart))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track_distr

    recursive subroutine exec_scale_project_distr( self, cline )
        class(scale_project_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(scale_project_distr_commander) :: xscale_distr
        type(qsys_env)                     :: qenv
        type(chash)                        :: job_descr
        type(cmdline)                      :: cline_scale
        type(chash),               allocatable :: part_params(:)
        character(len=LONGSTRLEN), allocatable :: part_stks(:)
        type(parameters)              :: params
        type(builder)                 :: build
        character(len=:), allocatable :: projfile_sc
        character(len=STDLEN) :: filetab
        integer, allocatable  :: parts(:,:)
        real                  :: smpd, smpd_target
        integer               :: istk, ipart, nparts, nstks, cnt, partsz, box, newbox
        logical               :: gen_sc_project
        ! mkdir=yes: a new *_sc project + stacks are generated
        ! mkdir=no : only stacks are scaled
        gen_sc_project = cline%get_carg('mkdir').eq.'yes'
        ! make parameters and project
        call build%init_params_and_build_spproj(cline, params)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! copy command line
        cline_scale = cline
        ! prepare part-dependent parameters
        nstks = build%spproj%os_stk%get_noris()
        if( nstks == 0 ) THROW_HARD('os_stk field of spproj empty; exec_scale_distr')
        if( cline%defined('nparts') )then
            nparts = min(params%nparts, nstks)
            call cline_scale%set('nparts', real(nparts))
        else
            nparts = 1
        endif
        smpd = build%spproj%get_smpd()
        box  = build%spproj%get_box()
        call cline_scale%set('smpd', smpd)
        call cline_scale%set('box',  real(box))
        if( gen_sc_project )then
            ! make new project & scales
            smpd_target = max(smpd, smpd * real(box)/real(params%newbox))
            call simple_mkdir(filepath(PATH_PARENT,'stack_parts_sc'), errmsg="commander_distr_wflows::exec_scale_project_distr ")
            call build%spproj%scale_projfile(smpd_target, projfile_sc, cline, cline_scale,&
                dir=filepath(PATH_PARENT,'stack_parts_sc'))
            newbox = nint(cline_scale%get_rarg('newbox'))
            if( newbox == box )then
                write(logfhandle,*)'Inconsistent input dimensions: from ',box,' to ',newbox
                THROW_HARD('inconsistent input dimensions; exec_scale_project_distr')
            endif
            call cline_scale%set('newbox', real(newbox))
            call xscale_distr%execute( cline_scale )
            ! delete copy in working directory
            call del_file(params%projfile)
        else
            ! scaling only
            params%nparts = nparts
            parts = split_nobjs_even(nstks, nparts)
            allocate(part_params(nparts))
            cnt = 0
            do ipart=1,nparts
                call part_params(ipart)%new(1)
                partsz = parts(ipart,2) - parts(ipart,1) + 1
                allocate(part_stks(partsz))
                ! creates part filetab
                filetab = 'scale_stktab_part'//int2str(ipart)//trim(TXT_EXT)
                do istk=1,partsz
                    cnt = cnt + 1
                    part_stks(istk) = build%spproj%get_stkname(cnt)
                enddo
                ! write part filetab & update part parameters
                call write_filetable( filetab, part_stks )
                call part_params(ipart)%set('filetab', filetab)
                deallocate(part_stks)
            end do
            deallocate(parts)
            ! setup the environment for distributed execution
            call qenv%new(nparts)
            ! prepare job description
            call cline_scale%gen_job_descr(job_descr)
            call job_descr%set('prg', 'scale')
            call job_descr%set('autoscale', 'no')
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
            ! clean
            call qsys_cleanup
            ! removes temporary split stktab lists
            do ipart=1,nparts
                filetab = 'scale_stktab_part'//int2str(ipart)//trim(TXT_EXT)
                call del_file( filetab )
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_SCALE NORMAL STOP ****')
    end subroutine exec_scale_project_distr

end module simple_commander_distr_wflows
