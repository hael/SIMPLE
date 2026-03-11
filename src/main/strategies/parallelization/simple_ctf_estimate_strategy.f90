! object-oriented strategy pattern for ctf_estimate
!
! Goals
!   - One unified workflow selected via runtime polymorphism
!   - Worker/shared-memory execution lives in inmem strategy
!   - Distributed-master scheduling + merge lives in distr strategy
!   - Hook-less: strategies implement the work directly
!   - No separate common module: shared helpers are private procedures below
!
module simple_ctf_estimate_strategy
use simple_commanders_api
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_qsys_env,       only: qsys_env
use simple_sp_project,     only: sp_project
implicit none
public :: ctf_estimate_strategy
public :: ctf_estimate_inmem_strategy
public :: ctf_estimate_distr_strategy
public :: create_ctf_estimate_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: ctf_estimate_strategy
contains
    procedure(apply_defaults_interface), deferred :: apply_defaults
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
end type ctf_estimate_strategy

! Worker/shared-memory implementation
type, extends(ctf_estimate_strategy) :: ctf_estimate_inmem_strategy
contains
    procedure :: apply_defaults => inmem_apply_defaults
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type ctf_estimate_inmem_strategy


! Distributed-master implementation
type, extends(ctf_estimate_strategy) :: ctf_estimate_distr_strategy
    type(sp_project) :: spproj
    type(chash)      :: job_descr
    type(qsys_env)   :: qenv
    integer          :: nintgs = 0
contains
    procedure :: apply_defaults => distr_apply_defaults
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type ctf_estimate_distr_strategy


abstract interface
    subroutine apply_defaults_interface(self, cline)
        import :: ctf_estimate_strategy, cmdline
        class(ctf_estimate_strategy), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
    end subroutine apply_defaults_interface

    subroutine init_interface(self, params, cline)
        import :: ctf_estimate_strategy, parameters, cmdline
        class(ctf_estimate_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        class(cmdline),              intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: ctf_estimate_strategy, parameters, cmdline
        class(ctf_estimate_strategy), intent(inout) :: self
        type(parameters),            intent(inout) :: params
        class(cmdline),              intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: ctf_estimate_strategy, parameters, cmdline
        class(ctf_estimate_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        class(cmdline),              intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: ctf_estimate_strategy, parameters, cmdline
        class(ctf_estimate_strategy), intent(inout) :: self
        type(parameters),            intent(in)    :: params
        class(cmdline),              intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: ctf_estimate_strategy
        class(ctf_estimate_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------

    function create_ctf_estimate_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(ctf_estimate_strategy), allocatable :: strategy
        logical :: is_master
        ! Master heuristic: nparts defined, but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(ctf_estimate_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED CTF_ESTIMATE (MASTER)'
        else
            allocate(ctf_estimate_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> CTF_ESTIMATE (SHARED-MEMORY / WORKER)'
        endif
    end function create_ctf_estimate_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults (private; no separate common module)
    ! --------------------------------------------------------------------

    subroutine set_ctf_estimate_defaults(cline)
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz',   512)
        if( .not. cline%defined('hp')      ) call cline%set('hp',        30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         5.)
        if( .not. cline%defined('dfmin')   ) call cline%set('dfmin', DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')   ) call cline%set('dfmax', DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')) call cline%set('ctfpatch', 'yes')
        ! Always work on the micrograph/movies segment
        call cline%set('oritype', 'mic')
    end subroutine set_ctf_estimate_defaults

    ! ====================================================================
    ! CTF_ESTIMATE (SHARED-MEMORY / WORKER)
    ! ====================================================================

    subroutine inmem_apply_defaults(self, cline)
        class(ctf_estimate_inmem_strategy), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        call set_ctf_estimate_defaults(cline)
    end subroutine inmem_apply_defaults

    subroutine inmem_initialize(self, params, cline)
        class(ctf_estimate_inmem_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        call params%new(cline)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_ctf_estimate_iter, only: ctf_estimate_iter
        use simple_binoris_io,        only: binwrite_oritab
        class(ctf_estimate_inmem_strategy), intent(inout) :: self
        type(parameters),                   intent(inout) :: params
        class(cmdline),                     intent(inout) :: cline
        type(sp_project)        :: spproj
        type(ctf_estimate_iter) :: ctfiter
        type(ctfparams)         :: ctfvars
        type(ori)               :: o
        type(string)            :: intg_forctf, output_dir, imgkind
        integer                 :: fromto(2), imic, ntot, cnt, state
        logical                 :: l_gen_thumb, l_del_forctf
        call spproj%read(params%projfile)
        ! sanity check
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        ! output directory
        output_dir = PATH_HERE
        ! parameters & loop range
        if( params%stream .eq. 'yes' )then
            ! streaming mode processes one at a time (as in original)
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_ctf_estimate')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! loop over micrographs
        cnt = 0
        do imic = fromto(1), fromto(2)
            cnt = cnt + 1
            call spproj%os_mic%get_ori(imic, o)
            state = 1
            if( o%isthere('state') ) state = o%get_state()
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                l_del_forctf = .false.
                if( o%isthere('forctf') )then
                    call o%getter('forctf', intg_forctf)
                    if( file_exists(intg_forctf) )then
                        l_del_forctf = .true.
                    else
                        if( o%isthere('intg') )then
                            call o%getter('intg', intg_forctf)
                        endif
                    endif
                else if( o%isthere('intg') )then
                    call o%getter('intg', intg_forctf)
                else
                    THROW_HARD('no image available (forctf|intg) for CTF fittings :: exec_ctf_estimate')
                endif
                l_gen_thumb = .not. o%isthere('thumb')
                ctfvars     = o%get_ctfvars()
                call ctfiter%iterate(params, ctfvars, intg_forctf, o, output_dir, l_gen_thumb)
                ! delete file after estimation
                if( l_del_forctf )then
                    call o%delete_entry('forctf')
                    call del_file(intg_forctf)
                endif
                ! update project
                call spproj%os_mic%set_ori(imic, o)
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        call spproj%kill
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(ctf_estimate_inmem_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_preprocess :: exec_ctf_estimate'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(ctf_estimate_inmem_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(ctf_estimate_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_CTF_ESTIMATE NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED CTF_ESTIMATE (MASTER)
    ! ====================================================================

    subroutine distr_apply_defaults(self, cline)
        class(ctf_estimate_distr_strategy), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        call set_ctf_estimate_defaults(cline)
    end subroutine distr_apply_defaults

    subroutine distr_initialize(self, params, cline)
        class(ctf_estimate_distr_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        integer :: nintgs_local
        call params%new(cline)
        ! sanity check
        call self%spproj%read_segment(params%oritype, params%projfile)
        nintgs_local = self%spproj%get_nintgs()
        if( nintgs_local == 0 )then
            THROW_HARD('no micrograph to process! exec_ctf_estimate_distr')
        endif
        call self%spproj%kill
        self%nintgs = nintgs_local
        ! (small robustness fix) ensure merge count is defined
        params%nptcls = self%nintgs
        ! (optional robustness) clamp nparts if too large
        if( params%nparts > self%nintgs )then
            call cline%set('nparts', self%nintgs)
            params%nparts = self%nintgs
        endif
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! setup environment
        call self%qenv%new(params, params%nparts)
        ! prepare job description
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(ctf_estimate_distr_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        ! schedule
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        ! merge docs
        call self%spproj%read(params%projfile)
        call self%spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call self%spproj%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(ctf_estimate_distr_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(ctf_estimate_distr_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(ctf_estimate_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_CTF_ESTIMATE NORMAL STOP ****'
    end function distr_end_message

end module simple_ctf_estimate_strategy
