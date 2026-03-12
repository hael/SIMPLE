! object-oriented strategy pattern for gen_pspecs_and_thumbs
!
! Goals
!   - One unified workflow selected via runtime polymorphism
!   - Worker/shared-memory execution lives in inmem strategy
!   - Distributed-master scheduling + merge lives in distr strategy
!   - Hook-less: strategies implement the work directly
!   - No separate common module: shared helpers are private procedures below
!
module simple_gen_pspecs_and_thumbs_strategy
use simple_commanders_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none

public :: gen_pspecs_and_thumbs_strategy
public :: gen_pspecs_and_thumbs_inmem_strategy
public :: gen_pspecs_and_thumbs_distr_strategy
public :: create_gen_pspecs_and_thumbs_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: gen_pspecs_and_thumbs_strategy
contains
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
end type gen_pspecs_and_thumbs_strategy


! Worker/shared-memory implementation
type, extends(gen_pspecs_and_thumbs_strategy) :: gen_pspecs_and_thumbs_inmem_strategy
contains
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type gen_pspecs_and_thumbs_inmem_strategy


! Distributed-master implementation
type, extends(gen_pspecs_and_thumbs_strategy) :: gen_pspecs_and_thumbs_distr_strategy
    type(sp_project) :: spproj
    type(qsys_env)   :: qenv
    type(chash)      :: job_descr
contains
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type gen_pspecs_and_thumbs_distr_strategy

abstract interface

    subroutine init_interface(self, params, cline)
        import :: gen_pspecs_and_thumbs_strategy, parameters, cmdline
        class(gen_pspecs_and_thumbs_strategy), intent(inout) :: self
        type(parameters),                      intent(inout) :: params
        class(cmdline),                        intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: gen_pspecs_and_thumbs_strategy, parameters, cmdline
        class(gen_pspecs_and_thumbs_strategy), intent(inout) :: self
        type(parameters),                      intent(inout) :: params
        class(cmdline),                        intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: gen_pspecs_and_thumbs_strategy, parameters, cmdline
        class(gen_pspecs_and_thumbs_strategy), intent(inout) :: self
        type(parameters),                      intent(in)    :: params
        class(cmdline),                        intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: gen_pspecs_and_thumbs_strategy, parameters, cmdline
        class(gen_pspecs_and_thumbs_strategy), intent(inout) :: self
        type(parameters),                      intent(in)    :: params
        class(cmdline),                        intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: gen_pspecs_and_thumbs_strategy
        class(gen_pspecs_and_thumbs_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------

    function create_gen_pspecs_and_thumbs_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(gen_pspecs_and_thumbs_strategy), allocatable :: strategy
        logical :: is_master
        ! Master heuristic: nparts defined, but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(gen_pspecs_and_thumbs_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED GEN_PSPECS_AND_THUMBS (MASTER)'
        else
            allocate(gen_pspecs_and_thumbs_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> GEN_PSPECS_AND_THUMBS (SHARED-MEMORY / WORKER)'
        endif
    end function create_gen_pspecs_and_thumbs_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults (private; no separate common module)
    ! --------------------------------------------------------------------

    subroutine set_gen_pspecs_and_thumbs_defaults(cline)
        class(cmdline), intent(inout) :: cline
        call cline%set('oritype', 'mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
    end subroutine set_gen_pspecs_and_thumbs_defaults

    ! ====================================================================
    ! GEN_PSPECS_AND_THUMBS (SHARED-MEMORY / WORKER)
    ! ====================================================================

    subroutine inmem_initialize(self, params, cline)
        class(gen_pspecs_and_thumbs_inmem_strategy), intent(inout) :: self
        type(parameters),                            intent(inout) :: params
        class(cmdline),                              intent(inout) :: cline
        call set_gen_pspecs_and_thumbs_defaults(cline)
        call params%new(cline)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_pspec_thumb_iter, only: pspec_thumb_iter
        use simple_binoris_io,       only: binwrite_oritab
        class(gen_pspecs_and_thumbs_inmem_strategy), intent(inout) :: self
        type(parameters),                            intent(inout) :: params
        class(cmdline),                              intent(inout) :: cline
        type(pspec_thumb_iter) :: ptiter
        type(sp_project)       :: spproj
        type(ori)              :: o
        type(string)           :: output_dir, moviename_intg, imgkind
        integer                :: nintgs, fromto(2), iintg, ntot, cnt
        call spproj%read(params%projfile)
        ! sanity check
        nintgs = spproj%get_nintgs()
        if( nintgs == 0 )then
            THROW_HARD('No integrated movies to process!')
        endif
        ! output directory
        output_dir = PATH_HERE
        ! determine loop range
        if( params%l_distr_exec )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; gen_pspecs_and_thumbs')
            endif
        else
            fromto = [1, nintgs]
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! iterate
        cnt = 0
        do iintg = fromto(1), fromto(2)
            call spproj%os_mic%get_ori(iintg, o)
            if( o%isthere('imgkind') .and. o%isthere('intg') )then
                cnt = cnt + 1
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', moviename_intg)
                call ptiter%iterate(params, o, moviename_intg, output_dir)
                call spproj%os_mic%set_ori(iintg, o)
                write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the integrated movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        call spproj%kill
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(gen_pspecs_and_thumbs_inmem_strategy), intent(inout) :: self
        type(parameters),                            intent(in)    :: params
        class(cmdline),                              intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_preprocess :: exec_gen_pspecs_and_thumbs'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(gen_pspecs_and_thumbs_inmem_strategy), intent(inout) :: self
        type(parameters),                            intent(in)    :: params
        class(cmdline),                              intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(gen_pspecs_and_thumbs_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_GEN_PSPECS_AND_THUMBS NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED GEN_PSPECS_AND_THUMBS (MASTER)
    ! ====================================================================

    subroutine distr_initialize(self, params, cline)
        class(gen_pspecs_and_thumbs_distr_strategy), intent(inout) :: self
        type(parameters),                             intent(inout) :: params
        class(cmdline),                               intent(inout) :: cline
        integer :: nintgs
        call set_gen_pspecs_and_thumbs_defaults(cline)
        call params%new(cline)
        ! numlen must match number of partitions (JOB_FINISHED naming)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! sanity check: integrated movies exist
        call self%spproj%read_segment(params%oritype, params%projfile)
        nintgs = self%spproj%get_nintgs()
        if( nintgs == 0 )then
            THROW_HARD('no integrated movies to process! exec_gen_pspecs_and_thumbs_distr')
        endif
        ! clamp nparts to nintgs (same behavior as original, but also refresh numlen)
        if( params%nparts > nintgs )then
            call cline%set('nparts', nintgs)
            params%nparts = nintgs
            params%numlen = len(int2str(params%nparts))
            call cline%set('numlen', params%numlen)
        endif
        ! ensure merge uses correct entry count (original code used params%nptcls but never set it)
        params%nptcls = nintgs
        call self%spproj%kill
        ! setup the environment for distributed execution
        call self%qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(gen_pspecs_and_thumbs_distr_strategy), intent(inout) :: self
        type(parameters),                             intent(inout) :: params
        class(cmdline),                               intent(inout) :: cline
        ! schedule & clean
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        ! merge docs
        call self%spproj%read(params%projfile)
        call self%spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call self%spproj%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(gen_pspecs_and_thumbs_distr_strategy), intent(inout) :: self
        type(parameters),                             intent(in)    :: params
        class(cmdline),                               intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(gen_pspecs_and_thumbs_distr_strategy), intent(inout) :: self
        type(parameters),                             intent(in)    :: params
        class(cmdline),                               intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(gen_pspecs_and_thumbs_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_GEN_PSPECS_AND_THUMBS NORMAL STOP ****'
    end function distr_end_message

end module simple_gen_pspecs_and_thumbs_strategy
