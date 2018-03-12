! concrete commander: distributed workflows
module simple_commander_distr_wflows
#include "simple_lib.f08"
use simple_cmdline,             only: cmdline
use simple_chash,               only: chash
use simple_qsys_env,            only: qsys_env
use simple_build,               only: build
use simple_params,              only: params
use simple_commander_base,      only: commander_base
use simple_commander_preprocess ! use all in there
use simple_commander_cluster2D  ! use all in there
use simple_commander_distr      ! use all in there
use simple_commander_mask       ! use all in there
use simple_commander_distr      ! use all in there
use simple_qsys_funs            ! use all in there
use simple_binoris_io           ! use all in there
implicit none

public :: motion_correct_ctf_estimate_distr_commander
public :: motion_correct_distr_commander
public :: motion_correct_tomo_distr_commander
public :: powerspecs_distr_commander
public :: ctf_estimate_distr_commander
public :: pick_distr_commander
public :: make_cavgs_distr_commander
public :: cluster2D_distr_commander
public :: prime3D_init_distr_commander
public :: prime3D_distr_commander
public :: reconstruct3D_distr_commander
public :: tseries_track_distr_commander
public :: symsrch_distr_commander
public :: scale_stk_parts_commander
private

type, extends(commander_base) :: motion_correct_ctf_estimate_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_ctf_estimate_distr
end type motion_correct_ctf_estimate_distr_commander
type, extends(commander_base) :: motion_correct_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_distr
end type motion_correct_distr_commander
type, extends(commander_base) :: motion_correct_tomo_distr_commander
  contains
    procedure :: execute      => exec_motion_correct_tomo_distr
end type motion_correct_tomo_distr_commander
type, extends(commander_base) :: powerspecs_distr_commander
  contains
    procedure :: execute      => exec_powerspecs_distr
end type powerspecs_distr_commander
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
type, extends(commander_base) :: prime3D_init_distr_commander
  contains
    procedure :: execute      => exec_prime3D_init_distr
end type prime3D_init_distr_commander
type, extends(commander_base) :: prime3D_distr_commander
  contains
    procedure :: execute      => exec_prime3D_distr
end type prime3D_distr_commander
type, extends(commander_base) :: reconstruct3D_distr_commander
  contains
    procedure :: execute      => exec_reconstruct3D_distr
end type reconstruct3D_distr_commander
type, extends(commander_base) :: tseries_track_distr_commander
  contains
    procedure :: execute      => exec_tseries_track_distr
end type tseries_track_distr_commander
type, extends(commander_base) :: symsrch_distr_commander
  contains
    procedure :: execute      => exec_symsrch_distr
end type symsrch_distr_commander
type, extends(commander_base) :: scale_stk_parts_commander
  contains
    procedure :: execute      => exec_scale_stk_parts
end type scale_stk_parts_commander

contains

    subroutine exec_motion_correct_ctf_estimate_distr( self, cline )
        use simple_commander_preprocess
        class(motion_correct_ctf_estimate_distr_commander), intent(inout) :: self
        class(cmdline),                                     intent(inout) :: cline
        type(merge_algndocs_commander) :: xmerge_algndocs
        type(cmdline)                  :: cline_merge_algndocs
        type(qsys_env)                 :: qenv
        type(params)                   :: p_master
        type(chash)                    :: job_descr
        character(len=:), allocatable  :: output_dir, fbody
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master = params(cline)
        p_master%nptcls = nlines(p_master%filetab)
        if( p_master%nparts > p_master%nptcls ) stop 'nr of partitions (nparts) mjust be < number of entries in filetable'
        ! deal with numlen so that length matches JOB_FINISHED indicator files
        p_master%numlen = len(int2str(p_master%nptcls))
        call cline%set('numlen', real(p_master%numlen))
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p_master%dir)
        else
            output_dir = trim(DIR_MOTION_CORRECT)
            call cline%set('dir', trim(output_dir))
        endif
        call mkdir(output_dir)
        allocate(fbody, source=trim(output_dir)//trim(UNIDOC_FBODY))
        ! prepare merge_algndocs command line
        cline_merge_algndocs = cline
        call cline_merge_algndocs%set( 'nthr',     1.                   )
        call cline_merge_algndocs%set( 'fbody',    trim(fbody)          )
        call cline_merge_algndocs%set( 'nptcls',   real(p_master%nptcls))
        call cline_merge_algndocs%set( 'ndocs',    real(p_master%nparts))
        call cline_merge_algndocs%set( 'outfile',  trim(SIMPLE_UNIDOC)  )
        call cline_merge_algndocs%set( 'numlen',   real(p_master%numlen))
        ! setup the environment for distributed execution
        call qenv%new(p_master, numlen=p_master%numlen)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, algnfbody=trim(UNIDOC_FBODY))
        ! merge docs
        call xmerge_algndocs%execute( cline_merge_algndocs )
        ! clean
        call qsys_cleanup(p_master)
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_motion_correct_ctf_estimate_distr

    subroutine exec_motion_correct_distr( self, cline )
        class(motion_correct_distr_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(merge_algndocs_commander) :: xmerge_algndocs
        type(cmdline)                  :: cline_merge_algndocs
        type(qsys_env)                 :: qenv
        type(params)                   :: p_master
        type(chash)                    :: job_descr
        character(len=:), allocatable  :: output_dir, fbody
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master        = params(cline)
        p_master%nptcls = nlines(p_master%filetab)
        if( p_master%nparts > p_master%nptcls ) stop 'nr of partitions (nparts) mjust be < number of entries in filetable'
        ! deal with numlen so that length matches JOB_FINISHED indicator files
        p_master%numlen = len(int2str(p_master%nptcls))
        call cline%set('numlen', real(p_master%numlen))
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p_master%dir)
        else
            output_dir = trim(DIR_MOTION_CORRECT)
            call cline%set('dir', trim(output_dir))
        endif
        call mkdir(output_dir)
        allocate(fbody, source=trim(output_dir)//trim(UNIDOC_FBODY))
        ! prepare merge_algndocs command line
        cline_merge_algndocs = cline
        call cline_merge_algndocs%set( 'nthr',     1.)
        call cline_merge_algndocs%set( 'fbody',    trim(fbody))
        call cline_merge_algndocs%set( 'nptcls',   real(p_master%nptcls))
        call cline_merge_algndocs%set( 'ndocs',    real(p_master%nparts))
        call cline_merge_algndocs%set( 'outfile',  trim(output_dir)//trim(SIMPLE_UNIDOC))
        call cline_merge_algndocs%set( 'numlen',   real(p_master%numlen))
        ! setup the environment for distributed execution
        call qenv%new(p_master, numlen=p_master%numlen)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, algnfbody=trim(UNIDOC_FBODY))
        ! merge docs
        call xmerge_algndocs%execute( cline_merge_algndocs )
        ! clean
        call qsys_cleanup(p_master)
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct_distr

    subroutine exec_motion_correct_tomo_distr( self, cline )
        use simple_oris,           only: oris
        class(motion_correct_tomo_distr_commander), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        character(len=STDLEN), allocatable :: tomonames(:)
        type(oris)               :: exp_doc
        integer                  :: nseries, ipart
        type(qsys_env)           :: qenv
        type(params)             :: p_master
        character(len=KEYLEN)    :: str
        type(chash)              :: job_descr
        type(chash), allocatable :: part_params(:)
        call cline%set('prg', 'motion_correct')
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master = params(cline)
        if( cline%defined('tomoseries') )then
            call read_filetable(p_master%tomoseries, tomonames)
        else
            stop 'need tomoseries (filetable of filetables) to be part of the command line when tomo=yes'
        endif
        nseries = size(tomonames)
        call exp_doc%new(nseries)
        if( cline%defined('exp_doc') )then
            if( file_exists(p_master%exp_doc) )then
                call exp_doc%read(p_master%exp_doc)
            else
                write(*,*) 'the required parameter file (flag exp_doc): ', trim(p_master%exp_doc)
                stop 'not in cwd'
            endif
        else
            stop 'need exp_doc (line: exp_time=X dose_rate=Y) to be part of the command line when tomo=yes'
        endif
        p_master%nparts = nseries
        p_master%nptcls = nseries
        ! prepare part-dependent parameters
        allocate(part_params(p_master%nparts), stat=alloc_stat) ! -1. is default excluded value
        allocchk("simple_commander_distr_wflows::motion_correct_tomo_moview_distr ")
        do ipart=1,p_master%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('filetab', trim(tomonames(ipart)))
            call part_params(ipart)%set('fbody', 'tomo'//int2str_pad(ipart,p_master%numlen_tomo))
            str = real2str(exp_doc%get(ipart,'exp_time'))
            call part_params(ipart)%set('exp_time', trim(str))
            str = real2str(exp_doc%get(ipart,'dose_rate'))
            call part_params(ipart)%set('dose_rate', trim(str))
        end do
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, part_params=part_params)
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT_TOMO NORMAL STOP ****')
    end subroutine exec_motion_correct_tomo_distr

    subroutine exec_powerspecs_distr( self, cline )
        class(powerspecs_distr_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(qsys_env) :: qenv
        type(params)   :: p_master
        type(chash)    :: job_descr
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master = params(cline)
        p_master%nptcls = nlines(p_master%filetab)
        if( p_master%nparts > p_master%nptcls ) stop 'nr of partitions (nparts) mjust be < number of entries in filetable'
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! clean
        call qsys_cleanup(p_master)
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_POWERSPECS NORMAL STOP ****')
    end subroutine exec_powerspecs_distr

    subroutine exec_ctf_estimate_distr( self, cline )
        class(ctf_estimate_distr_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(merge_algndocs_commander) :: xmerge_algndocs
        type(cmdline)                  :: cline_merge_algndocs
        type(params)                   :: p_master
        type(chash)                    :: job_descr
        type(qsys_env)                 :: qenv
        character(len=:), allocatable  :: output_dir, fbody
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master = params(cline)
        p_master%nptcls = nlines(p_master%filetab)
        if( p_master%nparts > p_master%nptcls ) stop 'nr of partitions (nparts) mjust be < number of entries in filetable'
        ! output directory
        if( cline%defined('dir') )then
            output_dir = trim(p_master%dir)
        else
            output_dir = trim(DIR_CTF_ESTIMATE)
            call cline%set('dir', trim(output_dir))
        endif
        call mkdir(output_dir)
        allocate(fbody, source=trim(output_dir)//trim('ctf_estimate_output_part'))
        ! prepare merge_algndocs command line
        cline_merge_algndocs = cline
        call cline_merge_algndocs%set( 'nthr',     1.)
        call cline_merge_algndocs%set( 'fbody',    trim(fbody) )
        call cline_merge_algndocs%set( 'nptcls',   real(p_master%nptcls) )
        call cline_merge_algndocs%set( 'ndocs',    real(p_master%nparts) )
        call cline_merge_algndocs%set( 'outfile',  trim(output_dir)//'ctf_estimate_output_merged'//trim(METADATA_EXT))
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! merge docs
        call xmerge_algndocs%execute( cline_merge_algndocs )
        ! clean
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_DISTR_ctf_estimate NORMAL STOP ****')
    end subroutine exec_ctf_estimate_distr

    subroutine exec_pick_distr( self, cline )
        class(pick_distr_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(qsys_env) :: qenv
        type(params)   :: p_master
        type(chash)    :: job_descr
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'stk')
        ! make master parameters
        p_master = params(cline)
        p_master%nptcls = nlines(p_master%filetab)
        if( p_master%nparts > p_master%nptcls ) stop 'nr of partitions (nparts) mjust be < number of entries in filetable'
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_DISTR_PICK NORMAL STOP ****')
    end subroutine exec_pick_distr

    subroutine exec_make_cavgs_distr( self, cline )
        class(make_cavgs_distr_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(split_commander) :: xsplit
        type(cmdline)         :: cline_cavgassemble
        type(qsys_env)        :: qenv
        type(params)          :: p_master
        type(chash)           :: job_descr
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set os segment in sp_project
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! make master parameters
        p_master = params(cline)
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! prepare command lines from prototype master
        cline_cavgassemble = cline
        call cline_cavgassemble%set('nthr',1.)
        call cline_cavgassemble%set('prg', 'cavgassemble')
        if( cline%defined('outfile') )then
            ! because outfile is output from distributed exec of make_cavgs
            call cline_cavgassemble%set('oritab', p_master%outfile)
        else
            ! because cluster2D_startdoc.METADATA_EXT is default output in the absence of outfile
            call cline_cavgassemble%set('oritab', 'cluster2D_startdoc'//trim(METADATA_EXT))
        endif
        if( .not. cline%defined('stktab') )then
            ! split stack
            call xsplit%execute(cline)
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! assemble class averages
        call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE', 'CAVGASSEMBLE_FINISHED')
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_make_cavgs_distr

    subroutine exec_cluster2D_distr( self, cline )
        use simple_procimgfile, only: random_selection_from_imgfile, copy_imgfile
        class(cluster2D_distr_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! commanders
        type(check2D_conv_commander)         :: xcheck2D_conv
        type(merge_algndocs_commander)       :: xmerge_algndocs
        type(split_commander)                :: xsplit
        type(make_cavgs_distr_commander)      :: xmake_cavgs
        ! command lines
        type(cmdline)         :: cline_check2D_conv
        type(cmdline)         :: cline_cavgassemble
        type(cmdline)         :: cline_merge_algndocs
        type(cmdline)         :: cline_make_cavgs
        ! other variables
        type(qsys_env)        :: qenv
        type(params)          :: p_master
        type(build)           :: b
        character(len=STDLEN) :: refs, refs_even, refs_odd, oritab, str, str_iter
        integer               :: iter
        type(chash)           :: job_descr
        real                  :: frac_srch_space
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        ! make master parameters
        p_master = params(cline)
        ! make builder
        call b%build_general_tbox(p_master, cline, do3d=.false.)
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! initialise starting references, orientations
        if( cline%defined('oritab') )then
            oritab=trim(p_master%oritab)
        else
            oritab=''
        endif
        if( cline%defined('refs') )then
            refs = trim(p_master%refs)
        else
            refs = 'start2Drefs' // p_master%ext
        endif

        ! prepare command lines from prototype master
        cline_check2D_conv     = cline
        cline_cavgassemble     = cline
        cline_merge_algndocs   = cline
        cline_make_cavgs       = cline

        ! initialise static command line parameters and static job description parameters
        call cline_merge_algndocs%set('fbody',  trim(ALGN_FBODY)     )
        call cline_merge_algndocs%set('nptcls', real(p_master%nptcls))
        call cline_merge_algndocs%set('ndocs',  real(p_master%nparts))
        call cline_check2D_conv%set('box',      real(p_master%box)   )
        call cline_check2D_conv%set('nptcls',   real(p_master%nptcls))
        call cline_cavgassemble%set('prg',      'cavgassemble'       )
        call cline_make_cavgs%set('prg',        'make_cavgs'         )
        if( job_descr%isthere('automsk') ) call job_descr%delete('automsk')

        if( .not. cline%defined('stktab') )then
            ! split stack
            call xsplit%execute(cline)
        endif

        ! execute initialiser
        if( .not. cline%defined('refs') )then
            p_master%refs      = 'start2Drefs'//p_master%ext
            p_master%refs_even = 'start2Drefs_even'//p_master%ext
            p_master%refs_odd  = 'start2Drefs_odd'//p_master%ext
            if( cline%defined('oritab') )then
                call cline_make_cavgs%set('refs', p_master%refs)
                call xmake_cavgs%execute(cline_make_cavgs)
            else
                if( cline%defined('stktab') )then
                    call random_selection_from_imgfile(p_master%stktab, p_master%refs,&
                        &p_master%ncls, p_master%box, p_master%smpd)
                else
                    call random_selection_from_imgfile(p_master%stk, p_master%refs,&
                        &p_master%ncls, p_master%box, p_master%smpd)
                endif
                call copy_imgfile(trim(p_master%refs), trim(p_master%refs_even), p_master%smpd, [1,p_master%ncls])
                call copy_imgfile(trim(p_master%refs), trim(p_master%refs_odd),  p_master%smpd, [1,p_master%ncls])
            endif
        endif

        ! extremal dynamics
        if( cline%defined('extr_iter') )then
            p_master%extr_iter = p_master%extr_iter - 1
        else
            p_master%extr_iter = p_master%startit - 1
        endif

        ! deal with eo partitioning
        if( b%a%get_nevenodd() == 0 )then
            if( p_master%tseries .eq. 'yes' )then
                call b%a%partition_eo(tseries=.true.)
            else
                call b%a%partition_eo
            endif
            if( cline%defined('oritab') )then
                call binwrite_oritab(p_master%oritab, b%spproj, b%a, [1,p_master%nptcls])
            else if( cline%defined('deftab') )then
                call binwrite_oritab(p_master%deftab, b%spproj, b%a, [1,p_master%nptcls])
            else
                p_master%deftab = 'deftab_from_distr_wflow'//trim(METADATA_EXT)
                call binwrite_oritab(p_master%deftab, b%spproj, b%a, [1,p_master%nptcls])
                call job_descr%set('deftab', trim(p_master%deftab))
                call cline%set('deftab', trim(p_master%deftab))
            endif
        endif

        ! main loop
        iter = p_master%startit - 1
        do
            iter = iter + 1
            str_iter = int2str_pad(iter,3)
            write(*,'(A)')   '>>>'
            write(*,'(A,I6)')'>>> ITERATION ', iter
            write(*,'(A)')   '>>>'
            ! cooling of the randomization rate
            p_master%extr_iter = p_master%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(p_master%extr_iter)))
            call cline%set('extr_iter', real(p_master%extr_iter))
            ! updates
            if( oritab .ne. '' ) call job_descr%set('oritab',  trim(oritab))
            call job_descr%set('refs', trim(refs))
            call job_descr%set('startit', int2str(iter))
            ! the only FRC we have is from the previous iteration, hence the iter - 1
            call job_descr%set('frcs', trim(FRCS_ITER_FBODY)//int2str_pad(iter - 1,3)//BIN_EXT)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, algnfbody=trim(ALGN_FBODY))
            ! merge orientation documents
            oritab = trim(CLUSTER2D_ITER_FBODY)//trim(str_iter)//trim(METADATA_EXT)
            call cline_merge_algndocs%set('outfile', trim(oritab))
            call xmerge_algndocs%execute(cline_merge_algndocs)
            ! assemble class averages
            refs      = trim(CAVGS_ITER_FBODY) // trim(str_iter)            // p_master%ext
            refs_even = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_even' // p_master%ext
            refs_odd  = trim(CAVGS_ITER_FBODY) // trim(str_iter) // '_odd'  // p_master%ext
            call cline_cavgassemble%set('oritab', trim(oritab))
            call cline_cavgassemble%set('which_iter', real(iter))
            call qenv%exec_simple_prg_in_queue(cline_cavgassemble, 'CAVGASSEMBLE', 'CAVGASSEMBLE_FINISHED')
            ! remapping of empty classes
            call remap_empty_cavgs
            ! check convergence
            call cline_check2D_conv%set('oritab', trim(oritab))
            call xcheck2D_conv%execute(cline_check2D_conv)
            frac_srch_space = 0.
            if( iter > 1 ) frac_srch_space = cline_check2D_conv%get_rarg('frac')
            ! the below activates shifting & automasking
            if( iter > 3 .and. (frac_srch_space >= FRAC_SH_LIM .or. cline_check2D_conv%defined('trs')) )then
                if( .not.job_descr%isthere('trs') )then
                    ! activates shift search
                    str = real2str(cline_check2D_conv%get_rarg('trs'))
                    call job_descr%set('trs', trim(str) )
                endif
                if( cline%defined('automsk') )then
                    ! activates masking
                    if( cline%get_carg('automsk') .ne. 'no' ) call job_descr%set('automsk','cavg')
                endif
            endif
            if( cline_check2D_conv%get_carg('converged').eq.'yes' .or. iter==p_master%maxits ) exit
        end do
        call qsys_cleanup(p_master)
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D NORMAL STOP ****')

        contains

            subroutine remap_empty_cavgs
                use simple_image,           only: image
                use simple_projection_frcs, only: projection_frcs
                type(image)           :: img_cavg
                type(projection_frcs) :: frcs
                integer, allocatable  :: fromtocls(:,:)
                character(len=STDLEN) :: frcs_iter
                integer               :: icls, state
                if( p_master%dyncls.eq.'yes' )then
                    call binread_oritab(oritab, b%spproj, b%a, [1,p_master%nptcls])
                    call b%a%fill_empty_classes(p_master%ncls, fromtocls)
                    if( allocated(fromtocls) )then
                        ! updates document
                        call binwrite_oritab(oritab, b%spproj, b%a, [1,p_master%nptcls])
                        ! updates refs
                        call img_cavg%new([p_master%box,p_master%box,1], p_master%smpd)
                        do icls = 1, size(fromtocls, dim=1)
                            call img_cavg%read(trim(refs), fromtocls(icls, 1))
                            call img_cavg%write(trim(refs), fromtocls(icls, 2))
                        enddo
                        call img_cavg%read(trim(refs), p_master%ncls)
                        call img_cavg%write(trim(refs), p_master%ncls)     ! to preserve size
                        do icls = 1, size(fromtocls, dim=1)
                            call img_cavg%read(trim(refs_even), fromtocls(icls, 1))
                            call img_cavg%write(trim(refs_even), fromtocls(icls, 2))
                        enddo
                        call img_cavg%read(trim(refs_even), p_master%ncls)
                        call img_cavg%write(trim(refs_even), p_master%ncls) ! to preserve size
                        do icls = 1, size(fromtocls, dim=1)
                            call img_cavg%read(trim(refs_odd), fromtocls(icls, 1))
                            call img_cavg%write(trim(refs_odd), fromtocls(icls, 2))
                        enddo
                        call img_cavg%read(trim(refs_odd), p_master%ncls)
                        call img_cavg%write(trim(refs_odd), p_master%ncls)  ! to preserve size
                        if( p_master%match_filt.eq.'yes')then
                            ! updates FRCs
                            state     = 1
                            frcs_iter = trim(FRCS_ITER_FBODY)//int2str_pad(iter,3)//'.bin'
                            call frcs%new(p_master%ncls, p_master%box, p_master%smpd, state)
                            call frcs%read(frcs_iter)
                            do icls = 1, size(fromtocls, dim=1)
                                call frcs%set_frc( fromtocls(icls,2),&
                                &frcs%get_frc(fromtocls(icls,1), p_master%box, state), state)
                            enddo
                            call frcs%write(frcs_iter)
                        endif
                    endif
                endif
            end subroutine remap_empty_cavgs

    end subroutine exec_cluster2D_distr

    subroutine exec_prime3D_init_distr( self, cline )
        use simple_commander_prime3D
        use simple_commander_rec
        class(prime3D_init_distr_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(split_commander) :: xsplit
        type(cmdline)         :: cline_volassemble
        type(qsys_env)        :: qenv
        type(params)          :: p_master
        character(len=STDLEN) :: vol
        type(chash)           :: job_descr
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! make master parameters
        p_master = params(cline)
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! init
        if( cline%defined('vol1') )then
            vol = trim(p_master%vols(1))
        else
            vol = 'startvol_state01'//p_master%ext
        endif
        if( .not. cline%defined('stktab') )then
            ! split stack
            call xsplit%execute(cline)
        endif
        ! prepare command lines from prototype master
        cline_volassemble = cline
        call cline_volassemble%set( 'outvol',  vol)
        call cline_volassemble%set( 'eo',     'no')
        call cline_volassemble%set( 'prg',    'volassemble')
        call cline_volassemble%set( 'oritab', 'prime3D_startdoc'//trim(METADATA_EXT))
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        call qenv%exec_simple_prg_in_queue(cline_volassemble, 'VOLASSEMBLE', 'VOLASSEMBLE_FINISHED')
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_DISTR_PRIME3D_INIT NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prime3D_init_distr

    subroutine exec_prime3D_distr( self, cline )
        use simple_commander_prime3D
        use simple_commander_mask
        use simple_commander_rec
        use simple_commander_volops
        use simple_oris,    only: oris
        use simple_math,    only: calc_fourier_index, calc_lowpass_lim
        use simple_strings, only: real2str
        class(prime3D_distr_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        ! commanders
        type(prime3D_init_distr_commander)   :: xprime3D_init_distr
        type(reconstruct3D_distr_commander)  :: xreconstruct3D_distr
        type(merge_algndocs_commander)       :: xmerge_algndocs
        type(check3D_conv_commander)         :: xcheck3D_conv
        type(split_commander)                :: xsplit
        type(postprocess_commander)         :: xpostprocess
        ! command lines
        type(cmdline)         :: cline_reconstruct3D_distr
        type(cmdline)         :: cline_refine3D_init
        type(cmdline)         :: cline_check3D_conv
        type(cmdline)         :: cline_merge_algndocs
        type(cmdline)         :: cline_volassemble
        type(cmdline)         :: cline_postprocess
        ! other variables
        type(qsys_env)        :: qenv
        type(params)          :: p_master
        type(build)           :: b
        type(chash)           :: job_descr
        character(len=STDLEN), allocatable :: state_assemble_finished(:)
        character(len=STDLEN) :: vol, vol_even, vol_odd, vol_iter, vol_iter_even
        character(len=STDLEN) :: vol_iter_odd, oritab, str, str_iter, optlp_file
        character(len=STDLEN) :: str_state, fsc_file, volassemble_output
        real                  :: frac_srch_space, corr, corr_prev
        integer               :: s, state, iter
        logical               :: vol_defined
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! make master parameters
        p_master = params(cline)
        ! make general builder to the the orientations in
        call b%build_general_tbox(p_master, cline, do3d=.false.)
        ! setup the environment for distributed execution
        call qenv%new(p_master)

        ! initialise
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
        call cline%set( 'box', real(p_master%box) )
        ! prepare command lines from prototype master
        cline_reconstruct3D_distr   = cline
        cline_refine3D_init   = cline
        cline_check3D_conv   = cline
        cline_merge_algndocs = cline
        cline_volassemble    = cline
        cline_postprocess   = cline

        ! initialise static command line parameters and static job description parameter
        call cline_reconstruct3D_distr%set( 'prg', 'reconstruct3D' )       ! required for distributed call
        call cline_refine3D_init%set( 'prg', 'prime3D_init' ) ! required for distributed call
        if( trim(p_master%refine).eq.'clustersym' ) call cline_reconstruct3D_distr%set( 'pgrp', 'c1' )
        call cline_merge_algndocs%set('nthr',   1.)
        call cline_merge_algndocs%set('fbody',  trim(ALGN_FBODY))
        call cline_merge_algndocs%set('nptcls', real(p_master%nptcls))
        call cline_merge_algndocs%set('ndocs',  real(p_master%nparts))
        call cline_check3D_conv%set('box',      real(p_master%box)   )
        call cline_check3D_conv%set('nptcls',   real(p_master%nptcls))
        call cline_postprocess%set('nstates',   1.)
        call cline_postprocess%set('mirr',      'no')

        ! for parallel volassemble over states
        allocate(state_assemble_finished(p_master%nstates) )

        ! removes unnecessary volume keys and generates volassemble finished names
        do state = 1,p_master%nstates
            vol = 'vol'//int2str( state )
            call cline_check3D_conv%delete( trim(vol) )
            call cline_merge_algndocs%delete( trim(vol) )
            call cline_volassemble%delete( trim(vol) )
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        ! split stack
        if( .not. cline%defined('stktab') )then
            call xsplit%execute(cline)
        endif

        ! GENERATE STARTING MODELS & ORIENTATIONS
        ! Orientations
        if( cline%defined('oritab') )then
            oritab=trim(p_master%oritab)
        else
            oritab='prime3D_startdoc'//trim(METADATA_EXT)
        endif
        ! Models
        vol_defined = .false.
        do state = 1,p_master%nstates
            vol = 'vol' // int2str(state)
            if( cline%defined(trim(vol)) )vol_defined = .true.
        enddo
        if( .not.cline%defined('oritab') .and. .not.vol_defined )then
            ! ab-initio
            call xprime3D_init_distr%execute( cline_refine3D_init )
            call cline%set( 'vol1', 'startvol_state01'//p_master%ext )
            call cline%set( 'oritab', oritab )
        else if( cline%defined('oritab') .and. .not.vol_defined )then
            ! reconstructions needed
            call xreconstruct3D_distr%execute( cline_reconstruct3D_distr )
            do state = 1,p_master%nstates
                ! rename volumes and update cline
                str_state = int2str_pad(state,2)
                vol = trim(VOL_FBODY)//trim(str_state)//p_master%ext
                str = trim(STARTVOL_FBODY)//trim(str_state)//p_master%ext
                call simple_rename( trim(vol), trim(str) )
                vol = 'vol'//trim(int2str(state))
                call cline%set( trim(vol), trim(str) )
                if( p_master%eo .ne. 'no' )then
                    vol_even = trim(VOL_FBODY)//trim(str_state)//'_even'//p_master%ext
                    str = trim(STARTVOL_FBODY)//trim(str_state)//'_even'//p_master%ext
                    call simple_rename( trim(vol_even), trim(str) )
                    vol_odd  = trim(VOL_FBODY)//trim(str_state)//'_odd' //p_master%ext
                    str = trim(STARTVOL_FBODY)//trim(str_state)//'_odd'//p_master%ext
                    call simple_rename( trim(vol_odd), trim(str) )
                endif
            enddo
        else if( .not.cline%defined('oritab') .and. vol_defined )then
            ! projection matching
            select case( p_master%neigh )
                case( 'yes' )
                    stop 'refinement method requires input orientation document'
                case DEFAULT
                    ! all good
            end select
        else
            ! all good
        endif
        ! EXTREMAL DYNAMICS
        if( cline%defined('extr_iter') )then
            p_master%extr_iter = p_master%extr_iter - 1
        else
            p_master%extr_iter = p_master%startit - 1
        endif
        ! EO PARTITIONING
        if( p_master%eo .ne. 'no' )then
            if( b%a%get_nevenodd() == 0 )then
                if( p_master%tseries .eq. 'yes' )then
                    call b%a%partition_eo(tseries=.true.)
                else
                    call b%a%partition_eo
                endif
                if( cline%defined('oritab') )then
                    call binwrite_oritab(p_master%oritab, b%spproj, b%a, [1,p_master%nptcls])
                else if( cline%defined('deftab') )then
                    call binwrite_oritab(p_master%deftab, b%spproj, b%a, [1,p_master%nptcls])
                else
                    p_master%deftab = 'deftab_from_distr_wflow'//trim(METADATA_EXT)
                    call binwrite_oritab(p_master%deftab, b%spproj, b%a, [1,p_master%nptcls])
                    call job_descr%set('deftab', trim(p_master%deftab))
                    call cline%set('deftab', trim(p_master%deftab))
                endif
            endif
        endif
        ! prepare Prime3D job description
        call cline%gen_job_descr(job_descr)
        ! MAIN LOOP
        iter = p_master%startit-1
        corr = -1.
        do
            iter = iter + 1
            str_iter = int2str_pad(iter,3)
            write(*,'(A)')   '>>>'
            write(*,'(A,I6)')'>>> ITERATION ', iter
            write(*,'(A)')   '>>>'
            if( cline%defined('oritab') )then
                call binread_oritab(trim(cline%get_carg('oritab')), b%spproj, b%a, [1,p_master%nptcls])
                frac_srch_space = b%a%get_avg('frac')
                call job_descr%set( 'oritab', trim(oritab) )
                if( p_master%refine .eq. 'snhc' )then
                    ! update stochastic neighborhood size if corr is not improving
                    corr_prev = corr
                    corr      = b%a%get_avg('corr')
                    if( iter > 1 .and. corr <= corr_prev )then
                        p_master%szsn = min(SZSN_MAX,p_master%szsn + SZSN_STEP)
                    endif
                    call job_descr%set('szsn', int2str(p_master%szsn))
                    call cline%set('szsn', real(p_master%szsn))
                endif
            endif
            ! exponential cooling of the randomization rate
            p_master%extr_iter = p_master%extr_iter + 1
            call job_descr%set('extr_iter', trim(int2str(p_master%extr_iter)))
            call cline%set('extr_iter', real(p_master%extr_iter))
            call job_descr%set( 'startit', trim(int2str(iter)) )
            call cline%set( 'startit', real(iter) )
            ! FRCs
            if( cline%defined('frcs') )then
                ! all good
            else
                call job_descr%set('frcs', trim(FRCS_FBODY)//'01'//BIN_EXT)
            endif
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, algnfbody=trim(ALGN_FBODY))
            ! ASSEMBLE ALIGNMENT DOCS
            if( p_master%refine .eq. 'snhc' )then
                oritab = trim(SNHCDOC)
            else
                oritab = trim(REFINE3D_ITER_FBODY)//trim(str_iter)//trim(METADATA_EXT)
            endif
            call cline%set( 'oritab', oritab )
            call cline_merge_algndocs%set( 'outfile', trim(oritab) )
            call xmerge_algndocs%execute( cline_merge_algndocs )
            ! ASSEMBLE VOLUMES
            call cline_volassemble%set( 'oritab', trim(oritab) )
            if( p_master%eo.ne.'no' )then
                call cline_volassemble%set( 'prg', 'volassemble_eo' ) ! required for cmdline exec
            else
                call cline_volassemble%set( 'prg', 'volassemble' )    ! required for cmdline exec
            endif
            do state = 1,p_master%nstates
                str_state = int2str_pad(state,2)
                if( p_master%eo .ne. 'no' )then
                    volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
                else
                    volassemble_output = ''
                endif
                call cline_volassemble%set( 'state', real(state) )
                call qenv%exec_simple_prg_in_queue(cline_volassemble, trim(volassemble_output),&
                    &script_name='simple_script_state'//trim(str_state))
            end do
            call qsys_watcher(state_assemble_finished)
            ! rename volumes, postprocess & update job_descr
            call binread_oritab(trim(oritab), b%spproj, b%a, [1,p_master%nptcls])
            do state = 1,p_master%nstates
                str_state = int2str_pad(state,2)
                if( b%a%get_pop( state, 'state' ) == 0 )then
                    ! cleanup for empty state
                    vol = 'vol'//trim(int2str(state))
                    call cline%delete( vol )
                    call job_descr%delete( trim(vol) )
                else
                    if( p_master%nstates>1 )then
                        ! cleanup postprocessing cmdline as it only takes one volume at a time
                        do s = 1,p_master%nstates
                            vol = 'vol'//int2str(s)
                            call cline_postprocess%delete( trim(vol) )
                        enddo
                    endif
                    ! rename state volume
                    vol      = trim(VOL_FBODY)//trim(str_state)//p_master%ext
                    vol_even = trim(VOL_FBODY)//trim(str_state)//'_even'//p_master%ext
                    vol_odd  = trim(VOL_FBODY)//trim(str_state)//'_odd' //p_master%ext
                    if( p_master%refine .eq. 'snhc' )then
                        vol_iter  = trim(SNHCVOL)//trim(str_state)//p_master%ext
                    else
                        vol_iter      = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//p_master%ext
                        vol_iter_even = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_even'//p_master%ext
                        vol_iter_odd  = trim(VOL_FBODY)//trim(str_state)//'_iter'//trim(str_iter)//'_odd' //p_master%ext
                    endif
                    call simple_rename( trim(vol), trim(vol_iter) )
                    if( p_master%eo .ne. 'no' )then
                        call simple_rename( trim(vol_even), trim(vol_iter_even) )
                        call simple_rename( trim(vol_odd),  trim(vol_iter_odd)  )
                    endif
                    ! post-process
                    vol = 'vol'//trim(int2str(state))
                    call cline_postprocess%set( 'vol1', trim(vol_iter))
                    fsc_file   = 'fsc_state'//trim(str_state)//'.bin'
                    optlp_file = 'aniso_optlp_state'//trim(str_state)//p_master%ext
                    if( file_exists(optlp_file) .and. p_master%eo .ne. 'no' )then
                        call cline_postprocess%delete('lp')
                        call cline_postprocess%set('fsc', trim(fsc_file))
                        call cline_postprocess%set('vol_filt', trim(optlp_file))
                    else if( file_exists(fsc_file) .and. p_master%eo .ne. 'no' )then
                        call cline_postprocess%delete('lp')
                        call cline_postprocess%set('fsc', trim(fsc_file))
                    else
                        call cline_postprocess%delete('fsc')
                        call cline_postprocess%set('lp', p_master%lp)
                    endif
                    call xpostprocess%execute(cline_postprocess)
                    ! updates cmdlines & job description
                    vol = 'vol'//trim(int2str(state))
                    call job_descr%set( trim(vol), trim(vol_iter) )
                    call cline%set( trim(vol), trim(vol_iter) )
                endif
            enddo
            ! CONVERGENCE
            call cline_check3D_conv%set( 'oritab', trim(oritab) )
            if(p_master%refine.eq.'cluster') call cline_check3D_conv%delete('update_res')
            call xcheck3D_conv%execute( cline_check3D_conv )
            if( iter >= p_master%startit+2 )then
                ! after a minimum of 2 iterations
                if( cline_check3D_conv%get_carg('converged') .eq. 'yes' ) exit
            endif
            if( iter >= p_master%maxits ) exit
            ! ITERATION DEPENDENT UPDATES
            if( cline_check3D_conv%defined('trs') .and. .not.job_descr%isthere('trs') )then
                ! activates shift search if frac >= 90
                str = real2str(cline_check3D_conv%get_rarg('trs'))
                call job_descr%set( 'trs', trim(str) )
                call cline%set( 'trs', cline_check3D_conv%get_rarg('trs') )
            endif
        end do
        call qsys_cleanup(p_master)
        ! report the last iteration on exit
        call cline%delete( 'startit' )
        call cline%set('endit', real(iter))
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_PRIME3D NORMAL STOP ****')
    end subroutine exec_prime3D_distr

    subroutine exec_reconstruct3D_distr( self, cline )
        use simple_commander_rec
        use simple_oris, only: oris
        class(reconstruct3D_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(split_commander)              :: xsplit
        type(qsys_env)                     :: qenv
        type(params)                       :: p_master
        type(cmdline)                      :: cline_volassemble
        character(len=STDLEN)              :: volassemble_output, str_state
        character(len=STDLEN), allocatable :: state_assemble_finished(:)
        type(build)                        :: b
        type(chash)                        :: job_descr
        integer                            :: state
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! make master parameters
        call cline%delete('refine')
        p_master = params(cline)
        ! make general builder to get oris in
        call b%build_general_tbox(p_master, cline, do3d=.false.)
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        call cline%gen_job_descr(job_descr)
        if( .not. cline%defined('stktab') )then
            ! split stack
            call xsplit%execute(cline)
        endif
        ! eo partitioning
        if( p_master%eo .ne. 'no' )then
            call binread_oritab(p_master%oritab, b%spproj, b%a, [1,p_master%nptcls])
            if( b%a%get_nevenodd() == 0 )then
                if( p_master%tseries .eq. 'yes' )then
                    call b%a%partition_eo(tseries=.true.)
                else
                    call b%a%partition_eo
                endif
                call binwrite_oritab(p_master%oritab, b%spproj, b%a, [1,p_master%nptcls])
            endif
        endif
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! assemble volumes
        ! this is for parallel volassemble over states
        allocate(state_assemble_finished(p_master%nstates) )
        do state = 1, p_master%nstates
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        cline_volassemble = cline
        if( p_master%eo .ne. 'no' )then
            call cline_volassemble%set('prg', 'volassemble_eo')
        else
            call cline_volassemble%set('prg', 'volassemble')
        endif
        ! parallel assembly
        do state = 1, p_master%nstates
            str_state = int2str_pad(state,2)
            if( p_master%eo .ne. 'no' )then
                volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
            else
                volassemble_output = ''
            endif
            call cline_volassemble%set( 'state', real(state) )
            call qenv%exec_simple_prg_in_queue(cline_volassemble, trim(volassemble_output),&
                &script_name='simple_script_state'//trim(str_state))
        end do
        call qsys_watcher(state_assemble_finished)
        ! termination
        call qsys_cleanup(p_master)
        call simple_end('**** SIMPLE_reconstruct3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D_distr

    subroutine exec_tseries_track_distr( self, cline )
        use simple_nrtxtfile,         only: nrtxtfile
        use simple_strings,           only: real2str
        class(tseries_track_distr_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(qsys_env)                :: qenv
        type(params)                  :: p_master
        type(chash)                   :: job_descr
        type(nrtxtfile)               :: boxfile
        real,        allocatable      :: boxdata(:,:)
        type(chash), allocatable      :: part_params(:)
        integer :: ndatlines, numlen, alloc_stat, j, orig_box, ipart
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! make master parameters
        p_master = params(cline)
        if( .not. file_exists(p_master%boxfile)  ) stop 'inputted boxfile does not exist in cwd'
        if( nlines(p_master%boxfile) > 0 )then
            call boxfile%new(p_master%boxfile, 1)
            ndatlines = boxfile%get_ndatalines()
            numlen    = len(int2str(ndatlines))
            allocate( boxdata(ndatlines,boxfile%get_nrecs_per_line()), stat=alloc_stat)
            allocchk('In: simple_commander_tseries :: exec_tseries_track')
            do j=1,ndatlines
                call boxfile%readNextDataLine(boxdata(j,:))
                orig_box = nint(boxdata(j,3))
                if( nint(boxdata(j,3)) /= nint(boxdata(j,4)) )then
                    stop 'Only square windows are currently allowed!'
                endif
            end do
        else
            stop 'inputted boxfile is empty; simple_commander_tseries :: exec_tseries_track'
        endif
        call boxfile%kill
        call cline%delete('boxfile')
        p_master%nptcls = ndatlines
        p_master%nparts = p_master%nptcls
        if( p_master%ncunits > p_master%nparts )&
        &stop 'nr of computational units (ncunits) mjust be <= number of entries in boxfiles'
        ! box and numlen need to be part of command line
        call cline%set('box',    real(orig_box))
        call cline%set('numlen', real(numlen)  )
        ! prepare part-dependent parameters
        allocate(part_params(p_master%nparts))
        do ipart=1,p_master%nparts
            call part_params(ipart)%new(3)
            call part_params(ipart)%set('xcoord', real2str(boxdata(ipart,1)))
            call part_params(ipart)%set('ycoord', real2str(boxdata(ipart,2)))
            call part_params(ipart)%set('ind',    int2str(ipart))
        end do
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! schedule & clean
        call cline%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, part_params=part_params)
        call qsys_cleanup(p_master)
        ! end gracefully
        call simple_end('**** SIMPLE_TSERIES_TRACK NORMAL STOP ****')
    end subroutine exec_tseries_track_distr

    subroutine exec_symsrch_distr( self, cline )
        use simple_comlin_srch,    only: comlin_srch_get_nproj
        use simple_commander_misc, only: sym_aggregate_commander
        use simple_sym,     only: sym
        use simple_ori,     only: ori
        use simple_oris,    only: oris
        use simple_strings, only: int2str_pad, int2str,real2str
        class(symsrch_distr_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(merge_algndocs_commander) :: xmerge_algndocs
        type(cmdline)                  :: cline_merge_algndocs
        type(cmdline)                  :: cline_gridsrch
        type(cmdline)                  :: cline_srch
        type(cmdline)                  :: cline_aggregate
        type(qsys_env)                 :: qenv
        type(params)                   :: p_master
        type(build)                    :: b
        type(chash)                    :: job_descr
        type(oris)                     :: tmp_os, sympeaks, o_shift, grid_symaxes
        type(ori)                      :: symaxis
        type(sym)                      :: syme
        integer,    allocatable        :: order(:)
        real,       allocatable        :: corrs(:)
        real                           :: shvec(3)
        integer                        :: i, comlin_srch_nproj, nl,  nbest_here
        integer                        :: bestloc(1), cnt, numlen
        character(len=STDLEN)          :: part_tab
        character(len=*), parameter :: GRIDSYMFBODY = 'grid_symaxes_part'           !<
        character(len=*), parameter :: GRIDSYMTAB   = 'grid_symaxes'//trim(TXT_EXT) !<
        character(len=*), parameter :: SYMFBODY     = 'symaxes_part'                !< symmetry axes doc (distributed mode)
        character(len=*), parameter :: SYMTAB       = 'symaxes'//trim(TXT_EXT)      !<
        character(len=*), parameter :: SYMPEAKSTAB  = 'sympeaks'//trim(TXT_EXT)     !< symmetry peaks to refine
        character(len=*), parameter :: SYMSHTAB     = 'sym_3dshift'//trim(TXT_EXT)  !< volume 3D shift
        character(len=*), parameter :: SYMPROJSTK   = 'sym_projs.mrc'               !< volume reference projections
        character(len=*), parameter :: SYMPROJTAB   = 'sym_projs'//trim(TXT_EXT)    !< volume reference projections doc
        integer,          parameter :: NBEST = 30
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls3D')
        ! make master parameters
        p_master          = params(cline)
        comlin_srch_nproj = comlin_srch_get_nproj( pgrp=trim(p_master%pgrp) )
        p_master%nptcls   = comlin_srch_nproj
        if( p_master%nparts > p_master%nptcls )then
            stop 'number of partitions (npart) > nr of jobs, adjust!'
        endif

        ! 1. GRID SEARCH
        cline_gridsrch = cline
        call cline_gridsrch%set('prg', 'symsrch')
        call cline_gridsrch%set('refine', 'no') !!
        call cline_gridsrch%set('fbody', trim(GRIDSYMFBODY))
        call qenv%new(p_master)
        call cline_gridsrch%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! consolidate grid search
        cline_merge_algndocs = cline
        call cline_merge_algndocs%set( 'nthr',     1.                      )
        call cline_merge_algndocs%set( 'fbody',    trim(GRIDSYMFBODY)      )
        call cline_merge_algndocs%set( 'nptcls',   real(comlin_srch_nproj) )
        call cline_merge_algndocs%set( 'ndocs',    real(p_master%nparts)   )
        call cline_merge_algndocs%set( 'outfile',  trim(GRIDSYMTAB)        )
        call xmerge_algndocs%execute( cline_merge_algndocs )

        ! 2. SELECTION OF SYMMETRY PEAKS TO REFINE
        nl = nlines(trim(GRIDSYMTAB))
        nbest_here = min(NBEST, nl)
        call grid_symaxes%new(nl)
        call grid_symaxes%read(trim(GRIDSYMTAB), [1,nl])
        call del_file(trim(GRIDSYMTAB))
        order = grid_symaxes%order_corr()
        call tmp_os%new(nbest_here)
        cnt = 0
        do i = nl, nl-nbest_here+1, -1
            cnt = cnt + 1
            call tmp_os%set_ori(cnt, grid_symaxes%get_ori(order(i)))
        enddo
        grid_symaxes = tmp_os
        call grid_symaxes%write(trim(GRIDSYMTAB), [1,nbest_here])
        deallocate(order)
        call tmp_os%kill

        ! 3. REFINEMENT
        cline_srch = cline
        call qsys_cleanup(p_master)
        call cline_srch%set('prg', 'symsrch')
        call cline_srch%set('refine', 'yes') !!
        call cline_srch%set('nptcls', real(comlin_srch_nproj))
        call cline_srch%set('oritab', trim(GRIDSYMTAB))
        call cline_srch%set('fbody',  trim(SYMFBODY))
        ! switch to collection of single threaded jobs
        p_master%ncunits = p_master%nparts * p_master%nthr
        p_master%nparts  = nbest_here
        p_master%nthr    = 1
        call qenv%new(p_master)
        call cline_srch%gen_job_descr(job_descr)
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr)
        ! merge_algndocs command line
        cline_merge_algndocs = cline
        call cline_merge_algndocs%set( 'nthr',     1.                    )
        call cline_merge_algndocs%set( 'fbody',    trim(SYMFBODY)        )
        call cline_merge_algndocs%set( 'nptcls',   real(nbest_here)      )
        call cline_merge_algndocs%set( 'ndocs',    real(p_master%nparts) )
        call cline_merge_algndocs%set( 'outfile',  trim(SYMTAB)          )
        call xmerge_algndocs%execute( cline_merge_algndocs )

        ! 4. REAL-SPACE EVALUATION
        cline_aggregate = cline
        call cline_aggregate%set( 'prg' ,   'sym_aggregate' )
        call cline_aggregate%set( 'nptcls',  real(p_master%nspace) )
        call cline_aggregate%set( 'oritab' , trim(SYMPROJTAB) )
        call cline_aggregate%set( 'oritab2', trim(SYMTAB) )
        call cline_aggregate%set( 'stk' ,    trim(SYMPROJSTK) )
        call cline_aggregate%set( 'outfile', trim(SYMPEAKSTAB) )
        call cline_aggregate%set( 'eo',      'no' )
        call qenv%exec_simple_prg_in_queue(cline_aggregate,&
        &'SYM_AGGREGATE', 'SYM_AGGREGATE_FINISHED')
        ! read and pick best
        nl = nlines(trim(SYMPEAKSTAB))
        call sympeaks%new(nl)
        call sympeaks%read(trim(SYMPEAKSTAB), [1,nl])
        corrs   = sympeaks%get_all('corr')
        bestloc = maxloc(corrs)
        symaxis = sympeaks%get_ori(bestloc(1))
        write(*,'(A)') '>>> FOUND SYMMETRY AXIS ORIENTATION:'
        call symaxis%print_ori()
        ! output
        if( cline%defined('oritab') )then
            ! transfer shift and symmetry to input orientations
            call syme%new(p_master%pgrp)
            call o_shift%new(1)
            ! retrieve shift
            call o_shift%read(trim(SYMSHTAB), [1,1])
            shvec(1) = o_shift%get(1,'x')
            shvec(2) = o_shift%get(1,'y')
            shvec(3) = o_shift%get(1,'z')
            shvec    = -1. * shvec ! the sign is right
            ! rotate the orientations & transfer the 3d shifts to 2d
            ! begin with making a general builder to support *.simple oritab input
            call b%build_general_tbox(p_master, cline, do3d=.false.)
            ! oritab is now read in, whatever format it had (.txt or .simple)
            nl = binread_nlines(p_master, p_master%oritab)
            call binread_oritab(p_master%oritab, b%spproj, b%a, [1,nl])
            if( cline%defined('state') )then
                call syme%apply_sym_with_shift(b%a, symaxis, shvec, p_master%state )
            else
                call syme%apply_sym_with_shift(b%a, symaxis, shvec )
            endif
            call binwrite_oritab(p_master%outfile, b%spproj, b%a, [1,nl])
            call b%kill_general_tbox
        endif

        ! THE END
        call qsys_cleanup(p_master)
        call del_file(trim(SYMSHTAB))
        numlen =  len(int2str(nbest_here))
        do i = 1, nbest_here
            part_tab = trim(SYMFBODY)//int2str_pad(i, numlen)//trim(TXT_EXT)
            call del_file(trim(part_tab))
        enddo
        p_master%nparts = nint(cline%get_rarg('nparts'))
        numlen =  len(int2str(p_master%nparts))
        do i = 1, p_master%nparts
            part_tab = trim(GRIDSYMFBODY)//int2str_pad(i, numlen)//trim(TXT_EXT)
            call del_file(trim(part_tab))
        enddo
        call del_file('SYM_AGGREGATE')
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_SYMSRCH NORMAL STOP ****')
    end subroutine exec_symsrch_distr

    subroutine exec_scale_stk_parts( self, cline )
        use simple_map_reduce, only: split_nobjs_even
        class(scale_stk_parts_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(qsys_env)                     :: qenv
        type(params)                       :: p_master
        type(chash)                        :: job_descr
        type(cmdline)                      :: cline_scale
        type(chash),           allocatable :: part_params(:)
        character(len=:),      allocatable :: stkname, outstkname
        character(len=STDLEN), allocatable :: part_stks(:)
        character(len=STDLEN) :: filetab
        integer, allocatable  :: parts(:,:)
        integer               :: ipart, nparts, nstks, cnt, istk, partsz
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! make master parameters
        p_master = params(cline)
        ! copy command line
        cline_scale = cline
        ! prepare part-dependent parameters
        if( cline%defined('stktab') )then
            call cline_scale%delete('stktab')
            nstks = p_master%stkhandle%get_nmics()
            if(cline%defined('nparts'))then
                nparts = min(p_master%nparts, nstks)
                call cline_scale%set('nparts',real(nparts))
                p_master%nparts = nparts
            else
                nparts = 1
            endif
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
                    part_stks(istk) = p_master%stkhandle%get_stkname(cnt)
                enddo
                ! write part filetab & update part parameters
                call write_filetable( filetab, part_stks )
                call part_params(ipart)%set('filetab', filetab)
                deallocate(part_stks)
            end do
            deallocate(parts)
        else
            ! prepare part-dependent parameters
            allocate(part_params(p_master%nparts))
            do ipart=1,p_master%nparts
                call part_params(ipart)%new(2)
                allocate(stkname, source=trim(STKPARTFBODY)//int2str_pad(ipart,p_master%numlen)//p_master%ext)
                outstkname = add2fbody(stkname, p_master%ext, '_sc')
                call part_params(ipart)%set('stk',    stkname)
                call part_params(ipart)%set('outstk', outstkname)
                deallocate(stkname, outstkname)
            end do
        endif
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline_scale%gen_job_descr(job_descr)
        call job_descr%set('prg', 'scale')
        call job_descr%set('autoscale', 'no')
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(p_master, job_descr, part_params=part_params)
        ! clean
        call qsys_cleanup(p_master)
        if( cline%defined('stktab') )then
            ! removes temporary split stktab lists
            do ipart=1,nparts
                filetab = 'scale_stktab_part'//int2str(ipart)//trim(TXT_EXT)
                call del_file( filetab )
            end do
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_SCALE_STK_PARTS NORMAL STOP ****')
    end subroutine exec_scale_stk_parts

end module simple_commander_distr_wflows
