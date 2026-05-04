! object-oriented strategy pattern for motion_correct
!
! Goals
!   - One unified workflow selected via runtime polymorphism
!   - Worker/shared-memory execution lives in inmem strategy
!   - Distributed-master scheduling + merge lives in distr strategy
!   - Hook-less: strategies implement the work directly
!   - No separate common module: shared helpers are private procedures below
!
module simple_motion_correct_strategy
use simple_commanders_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none

public :: motion_correct_strategy
public :: motion_correct_inmem_strategy
public :: motion_correct_distr_strategy
public :: create_motion_correct_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: motion_correct_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
    procedure(endmsg_interface),   deferred :: end_message
end type motion_correct_strategy


! Worker/shared-memory implementation
type, extends(motion_correct_strategy) :: motion_correct_inmem_strategy
contains
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type motion_correct_inmem_strategy


! Distributed-master implementation
type, extends(motion_correct_strategy) :: motion_correct_distr_strategy
    type(sp_project) :: spproj
    type(qsys_env)   :: qenv
    type(chash)      :: job_descr
contains
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type motion_correct_distr_strategy


abstract interface
    subroutine init_interface(self, params, cline)
        import :: motion_correct_strategy, parameters, cmdline
        class(motion_correct_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: motion_correct_strategy, parameters, cmdline
        class(motion_correct_strategy), intent(inout) :: self
        type(parameters),              intent(inout) :: params
        class(cmdline),                intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: motion_correct_strategy, parameters, cmdline
        class(motion_correct_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: motion_correct_strategy, parameters, cmdline
        class(motion_correct_strategy), intent(inout) :: self
        type(parameters),              intent(in)    :: params
        class(cmdline),                intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: motion_correct_strategy
        class(motion_correct_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------

    function create_motion_correct_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(motion_correct_strategy), allocatable :: strategy
        logical :: is_master
        ! Master heuristic: nparts defined, but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(motion_correct_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED MOTION_CORRECT (MASTER)'
        else
            allocate(motion_correct_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> MOTION_CORRECT (SHARED-MEMORY / WORKER)'
        endif
    end function create_motion_correct_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults (private; no separate common module)
    ! --------------------------------------------------------------------

    subroutine set_motion_correct_defaults(cline)
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('trs')            ) call cline%set('trs',               20.)
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',            8.)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',             5.)
        if( .not. cline%defined('bfac')           ) call cline%set('bfac',              50.)
        if( .not. cline%defined('mcconvention')   ) call cline%set('mcconvention', 'simple')
        if( .not. cline%defined('wcrit')          ) call cline%set('wcrit',       'softmax')
        if( .not. cline%defined('eer_upsampling') ) call cline%set('eer_upsampling',      1)
        if( .not. cline%defined('mcpatch')        ) call cline%set('mcpatch',         'yes')
        if( .not. cline%defined('mcpatch_thres')  ) call cline%set('mcpatch_thres',   'yes')
        if( .not. cline%defined('algorithm')      ) call cline%set('algorithm',     'patch')
        ! Always motion_correct micrographs/movies segment
        call cline%set('oritype', 'mic')
    end subroutine set_motion_correct_defaults

    ! ====================================================================
    ! MOTION_CORRECT (SHARED-MEMORY / WORKER)
    ! ====================================================================

    subroutine inmem_initialize(self, params, cline)
        class(motion_correct_inmem_strategy), intent(inout) :: self
        type(parameters),                     intent(inout) :: params
        class(cmdline),                       intent(inout) :: cline
        call set_motion_correct_defaults(cline)
        call params%new(cline)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_binoris_io,          only: binwrite_oritab
        class(motion_correct_inmem_strategy), intent(inout) :: self
        type(parameters),                     intent(inout) :: params
        class(cmdline),                       intent(inout) :: cline
        type(motion_correct_iter) :: mciter
        type(ctfparams)           :: ctfvars
        type(sp_project)          :: spproj
        type(ori)                 :: o
        type(string)              :: output_dir, moviename, fbody
        integer                   :: nmovies, fromto(2), imovie, ntot, frame_counter, cnt
        call spproj%read(params%projfile)
        ! sanity check
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('No movie to process!')
        endif
        if( params%scale_movies > 1.01 )then
            THROW_HARD('scale_movies cannot be > 1; exec_motion_correct')
        endif
        if( cline%defined('gainref') )then
            if(.not.file_exists(params%gainref) )then
                THROW_HARD('gain reference: '//params%gainref%to_char()//' not found; motion_correct')
            endif
        endif
        ! output directory & names
        output_dir = PATH_HERE
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = ''
        endif
        ! determine loop range
        if( cline%defined('fromp') .and. cline%defined('top') )then
            fromto = [params%fromp, params%top]
        else
            THROW_HARD('fromp & top args need to be defined in parallel execution; motion_correct')
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! align
        frame_counter = 0
        cnt = 0
        do imovie = fromto(1), fromto(2)
            call spproj%os_mic%get_ori(imovie, o)
            if( o%isthere('imgkind') )then
                if( o%isthere('movie') .or. o%isthere('mic') )then
                    cnt = cnt + 1
                    call o%getter('movie', moviename)
                    ctfvars = spproj%get_micparams(imovie)
                    if( cline%defined('gainref') )then
                        call mciter%iterate(params, cline, ctfvars, o, fbody, frame_counter, moviename, &
                            &output_dir, gainref_fname=params%gainref)
                    else
                        call mciter%iterate(params, cline, ctfvars, o, fbody, frame_counter, moviename, output_dir)
                    endif
                    call spproj%os_mic%set_ori(imovie, o)
                    write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
                endif
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! cleanup locals
        call o%kill
        call spproj%kill
        call output_dir%kill
        call moviename%kill
        call fbody%kill
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(motion_correct_inmem_strategy), intent(inout) :: self
        type(parameters),                     intent(in)    :: params
        class(cmdline),                       intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_preprocess :: exec_motion_correct'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(motion_correct_inmem_strategy), intent(inout) :: self
        type(parameters),                     intent(in)    :: params
        class(cmdline),                       intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(motion_correct_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_MOTION_CORRECT NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED MOTION_CORRECT (MASTER)
    ! ====================================================================

    subroutine distr_initialize(self, params, cline)
        use simple_motion_correct_utils, only: flip_gain
        class(motion_correct_distr_strategy), intent(inout) :: self
        type(parameters),                     intent(inout) :: params
        class(cmdline),                       intent(inout) :: cline
        integer :: nmovies
        call set_motion_correct_defaults(cline)
        call params%new(cline)
        ! Preserve original behavior: set numlen from params as parsed/derived
        call cline%set('numlen', params%numlen)
        ! sanity check: movie segment exists and has entries
        call self%spproj%read_segment(params%oritype, params%projfile)
        nmovies = self%spproj%get_nmovies()
        if( nmovies == 0 ) THROW_HARD('no movies to process! exec_motion_correct_distr')
        ! (Helpful consistency fix) ensure merge uses a valid count
        params%nptcls = nmovies
        if( params%nparts > params%nptcls )then
            THROW_HARD('# partitions (nparts) must be < number of entries in filetable')
        endif
        call self%spproj%kill
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! setup distributed environment
        call self%qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(motion_correct_distr_strategy), intent(inout) :: self
        type(parameters),                     intent(inout) :: params
        class(cmdline),                       intent(inout) :: cline
        ! schedule & clean
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        ! merge docs
        call self%spproj%read(params%projfile)
        call self%spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call self%spproj%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(motion_correct_distr_strategy), intent(inout) :: self
        type(parameters),                     intent(in)    :: params
        class(cmdline),                       intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(motion_correct_distr_strategy), intent(inout) :: self
        type(parameters),                     intent(in)    :: params
        class(cmdline),                       intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(motion_correct_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****'
    end function distr_end_message

end module simple_motion_correct_strategy