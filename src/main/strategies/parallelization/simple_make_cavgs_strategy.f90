module simple_make_cavgs_strategy
use simple_core_module_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_builder,    only: builder
use simple_classaverager, only: cavger_new, cavger_transf_oridat, cavger_read_euclid_sigma2, cavger_assemble_sums, &
                              cavger_restore_cavgs, cavger_gen2Dclassdoc, cavger_write_all, cavger_kill, &
                              cavger_readwrite_partial_sums
use simple_qsys_env,   only: qsys_env
use simple_qsys_funs,  only: qsys_job_finished, qsys_cleanup
use simple_exec_helpers, only: set_shmem_flag, set_master_num_threads
implicit none

public :: make_cavgs_hooks
public :: make_cavgs_strategy
public :: make_cavgs_shmem_strategy
public :: make_cavgs_worker_strategy
public :: make_cavgs_master_strategy
public :: create_make_cavgs_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Commander-provided hook(s)
! --------------------------------------------------------------------
abstract interface
    subroutine cavgassemble_hook(cline, nthr)
        import :: cmdline
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: nthr
    end subroutine cavgassemble_hook
end interface

type :: make_cavgs_hooks
    procedure(cavgassemble_hook), pointer, nopass :: run_cavgassemble => null()
end type make_cavgs_hooks

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------
type, abstract :: make_cavgs_strategy
    type(make_cavgs_hooks) :: hooks
    logical :: l_from_distr_cmd = .false.  ! entrypoint was exec_make_cavgs_distr
    logical :: l_touch_async    = .false.  ! touch MAKECAVGS_FINISHED when async=yes
contains
    procedure(apply_defaults_interface), deferred :: apply_defaults
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
    procedure :: after_end => base_after_end
end type make_cavgs_strategy

type, extends(make_cavgs_strategy) :: make_cavgs_shmem_strategy
    type(builder) :: build
contains
    procedure :: apply_defaults => shmem_apply_defaults
    procedure :: initialize     => shmem_initialize
    procedure :: execute        => shmem_execute
    procedure :: finalize_run   => shmem_finalize_run
    procedure :: cleanup        => shmem_cleanup
    procedure :: end_message    => shmem_end_message
end type make_cavgs_shmem_strategy

type, extends(make_cavgs_strategy) :: make_cavgs_worker_strategy
    type(builder) :: build
contains
    procedure :: apply_defaults => worker_apply_defaults
    procedure :: initialize     => worker_initialize
    procedure :: execute        => worker_execute
    procedure :: finalize_run   => worker_finalize_run
    procedure :: cleanup        => worker_cleanup
    procedure :: end_message    => worker_end_message
end type make_cavgs_worker_strategy

type, extends(make_cavgs_strategy) :: make_cavgs_master_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    integer        :: nthr_master = 1
contains
    procedure :: apply_defaults => master_apply_defaults
    procedure :: initialize     => master_initialize
    procedure :: execute        => master_execute
    procedure :: finalize_run   => master_finalize_run
    procedure :: cleanup        => master_cleanup
    procedure :: end_message    => master_end_message
end type make_cavgs_master_strategy

abstract interface
    subroutine apply_defaults_interface(self, cline)
        import :: make_cavgs_strategy, cmdline
        class(make_cavgs_strategy), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
    end subroutine apply_defaults_interface

    subroutine init_interface(self, params, cline)
        import :: make_cavgs_strategy, parameters, cmdline
        class(make_cavgs_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: make_cavgs_strategy, parameters, cmdline
        class(make_cavgs_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: make_cavgs_strategy, parameters, cmdline
        class(make_cavgs_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: make_cavgs_strategy, parameters, cmdline
        class(make_cavgs_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: make_cavgs_strategy
        class(make_cavgs_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! ====================================================================
    ! Base "after_end" action (preserve exec_make_cavgs_distr async touch)
    ! ====================================================================
    subroutine base_after_end(self, params, cline)
        class(make_cavgs_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
        if( self%l_touch_async )then
            if( trim(params%async) .eq. 'yes' ) call simple_touch(MAKECAVGS_FINISHED)
        endif
    end subroutine base_after_end

    ! ====================================================================
    ! Factory
    ! ====================================================================
    function create_make_cavgs_strategy(cline, hooks, from_distr_cmd) result(strategy)
        class(cmdline),         intent(in) :: cline
        type(make_cavgs_hooks), intent(in) :: hooks
        logical, optional,      intent(in) :: from_distr_cmd
        class(make_cavgs_strategy), allocatable :: strategy
        integer :: nparts
        logical :: is_worker, is_master, l_from
        l_from = .false.
        if( present(from_distr_cmd) ) l_from = from_distr_cmd
        if( cline%defined('nparts') )then
            nparts = cline%get_iarg('nparts')
        else
            nparts = 1
        endif
        is_worker = cline%defined('part') .or. (cline%defined('fromp') .and. cline%defined('top') .and. (nparts > 1))
        is_master = (nparts > 1) .and. (.not. is_worker)
        if( is_master )then
            allocate(make_cavgs_master_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED MAKE_CAVGS (MASTER)'
        else if( is_worker )then
            allocate(make_cavgs_worker_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> MAKE_CAVGS (DISTRIBUTED WORKER)'
        else
            allocate(make_cavgs_shmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> MAKE_CAVGS (SHARED-MEMORY)'
        endif
        strategy%hooks            = hooks
        strategy%l_from_distr_cmd = l_from
        strategy%l_touch_async    = l_from
    end function create_make_cavgs_strategy

    ! ====================================================================
    ! Shared helpers (private) — THIS IS THE OVERLAP
    ! ====================================================================

    subroutine apply_distr_entry_defaults(cline)
        class(cmdline), intent(inout) :: cline
        call cline%set('wiener', 'full')
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('ml_reg') ) call cline%set('ml_reg','no')
    end subroutine apply_distr_entry_defaults

    subroutine normalize_ncls_nspace_oritype(cline, set_ncls_from_nspace)
        class(cmdline), intent(inout) :: cline
        logical,        intent(in)    :: set_ncls_from_nspace
        ! if( (cline%defined('ncls')) .and. cline%defined('nspace') )then
        !     THROW_HARD('NCLS and NSPACE cannot be both defined!')
        ! endif
        ! if( cline%defined('nspace') )then
        !     if( cline%defined('oritype') )then
        !         if( cline%get_carg('oritype') .eq. 'ptcl2D' )then
        !             THROW_HARD('NSPACE & PTCL2D are incompatible!')
        !         endif
        !     endif
        !     call cline%set('oritype', 'ptcl3D')
        !     if( set_ncls_from_nspace ) call cline%set('ncls', cline%get_iarg('nspace'))
        ! else
            call cline%set('oritype', 'ptcl2D')
        ! endif
    end subroutine normalize_ncls_nspace_oritype

    subroutine generate_cluster_centers(build, params, cline)
        type(builder),     intent(inout) :: build
        type(parameters),  intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        integer :: ncls_here
        if( L_VERBOSE_GLOB ) write(logfhandle,'(a)') '>>> GENERATING CLUSTER CENTERS'
        if( trim(params%oritype) .eq. 'ptcl3D' )then
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
            call build%spproj%os_ptcl3D%set_projs(build%eulspace)
            call build%spproj%os_ptcl3D%proj2class
        else
            ncls_here = build%spproj_field%get_n('class')
            if( .not. cline%defined('ncls') ) params%ncls = build%spproj_field%get_n('class')

            if( params%l_remap_cls )then
                call build%spproj_field%remap_cls()
                if( cline%defined('ncls') )then
                    if( params%ncls < ncls_here ) THROW_HARD('inputted ncls < ncls_in_oritab not allowed!')
                    if( params%ncls > ncls_here )then
                        call build%spproj_field%expand_classes(params%ncls)
                    endif
                endif
            else if( params%tseries .eq. 'yes' )then
                if( .not. cline%defined('ncls') )then
                    THROW_HARD('# class averages (ncls) need to be part of command line when tseries=yes')
                endif
                call build%spproj_field%ini_tseries(params%ncls, 'class')
                call build%spproj_field%partition_eo
            else if( params%proj_is_class .eq. 'yes' )then
                call build%spproj_field%proj2class
            endif
        endif
    end subroutine generate_cluster_centers

    subroutine setup_weights_and_evenodd(build, params)
        type(builder),    intent(inout) :: build
        type(parameters), intent(inout) :: params
        call build%spproj_field%calc_hard_weights2D(params%frac, params%ncls)
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
    end subroutine setup_weights_and_evenodd

    subroutine write_oritype_segment_for_compute(build, params, l_shmem)
        type(builder),    intent(inout) :: build
        type(parameters), intent(inout) :: params
        logical,          intent(in)    :: l_shmem
        if( l_shmem )then
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        else
            if( params%part .eq. 1 ) call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
    end subroutine write_oritype_segment_for_compute

    subroutine cavger_prepare_and_assemble(params, build)
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        call cavger_new(params, build)
        call cavger_read_euclid_sigma2
        call cavger_assemble_sums(.false.)
    end subroutine cavger_prepare_and_assemble

    subroutine shmem_bookkeeping(build, params, cline)
        type(builder),    intent(inout) :: build
        type(parameters), intent(inout) :: params
        class(cmdline),   intent(inout) :: cline
        integer :: icls
        select case(trim(params%oritype))
            case('ptcl2D')
                call build%spproj%write_segment_inside('cls2D', params%projfile)
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc2D')
                call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg')
                call build%spproj%write_segment_inside('out', params%projfile)
            case('ptcl3D')
                do icls = 1, params%nspace
                    call build%spproj%os_cls3D%set_euler(icls, build%eulspace%get_euler(icls))
                enddo
                if( cline%defined('outfile') )then
                    call build%spproj%os_cls3D%write(params%outfile)
                else
                    call build%spproj%os_cls3D%write(string('cls3D_oris.txt'))
                endif
                call build%spproj%write_segment_inside('cls3D', params%projfile)
                call build%spproj%add_frcs2os_out(string(FRCS_FILE), 'frc3D')
                call build%spproj%add_cavgs2os_out(params%refs, build%spproj%get_smpd(), imgkind='cavg3D')
                call build%spproj%write_segment_inside('out', params%projfile)
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE: '//trim(params%oritype))
        end select
    end subroutine shmem_bookkeeping

    ! ====================================================================
    ! SHMEM STRATEGY — uses shared helpers + shmem tail
    ! ====================================================================

    subroutine shmem_apply_defaults(self, cline)
        class(make_cavgs_shmem_strategy), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        if( self%l_from_distr_cmd ) call apply_distr_entry_defaults(cline)
        call normalize_ncls_nspace_oritype(cline, set_ncls_from_nspace=.true.)
    end subroutine shmem_apply_defaults

    subroutine shmem_initialize(self, params, cline)
        class(make_cavgs_shmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        logical :: l_shmem
        l_shmem = set_shmem_flag(cline)
        if( l_shmem .and. (.not. cline%defined('refs')) )then
            THROW_HARD('need input refs (filename) for shared-memory execution')
        endif
        call self%build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
    end subroutine shmem_initialize

    subroutine shmem_execute(self, params, cline)
        class(make_cavgs_shmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        ! ---- SHARED BODY (overlap with worker) ----
        call generate_cluster_centers(self%build, params, cline)
        call setup_weights_and_evenodd(self%build, params)
        call write_oritype_segment_for_compute(self%build, params, l_shmem=.true.)
        call cavger_prepare_and_assemble(params, self%build)
        ! ------------------------------------------
        ! ---- SHMEM TAIL (differs from worker) ----
        call cavger_restore_cavgs(params%frcs)
        call cavger_gen2Dclassdoc
        call cavger_write_all(params%refs, params%refs_even, params%refs_odd)
        call cavger_kill
        call shmem_bookkeeping(self%build, params, cline)
        ! ------------------------------------------
    end subroutine shmem_execute

    subroutine shmem_finalize_run(self, params, cline)
        class(make_cavgs_shmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        ! no-op
    end subroutine shmem_finalize_run

    subroutine shmem_cleanup(self, params, cline)
        class(make_cavgs_shmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        call self%build%kill_strategy2D_tbox
        call self%build%kill_general_tbox
    end subroutine shmem_cleanup

    function shmem_end_message(self) result(msg)
        class(make_cavgs_shmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_MAKE_CAVGS NORMAL STOP ****'
    end function shmem_end_message

    ! ====================================================================
    ! WORKER STRATEGY — uses shared helpers + worker tail
    ! ====================================================================

    subroutine worker_apply_defaults(self, cline)
        class(make_cavgs_worker_strategy), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        if( self%l_from_distr_cmd ) call apply_distr_entry_defaults(cline)
        call normalize_ncls_nspace_oritype(cline, set_ncls_from_nspace=.true.)
    end subroutine worker_apply_defaults

    subroutine worker_initialize(self, params, cline)
        class(make_cavgs_worker_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        call self%build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.true.)
    end subroutine worker_initialize

    subroutine worker_execute(self, params, cline)
        class(make_cavgs_worker_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        ! ---- SHARED BODY (overlap with shmem) ----
        call generate_cluster_centers(self%build, params, cline)
        call setup_weights_and_evenodd(self%build, params)
        call write_oritype_segment_for_compute(self%build, params, l_shmem=.false.)
        call cavger_prepare_and_assemble(params, self%build)
        ! ------------------------------------------
        ! ---- WORKER TAIL (differs from shmem) ----
        call cavger_readwrite_partial_sums('write')
        call cavger_kill
        ! ------------------------------------------
    end subroutine worker_execute

    subroutine worker_finalize_run(self, params, cline)
        class(make_cavgs_worker_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_cluster2D :: exec_make_cavgs'))
    end subroutine worker_finalize_run

    subroutine worker_cleanup(self, params, cline)
        class(make_cavgs_worker_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        call self%build%kill_strategy2D_tbox
        call self%build%kill_general_tbox
    end subroutine worker_cleanup

    function worker_end_message(self) result(msg)
        class(make_cavgs_worker_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_MAKE_CAVGS NORMAL STOP ****'
    end function worker_end_message

    ! ====================================================================
    ! MASTER STRATEGY — unchanged conceptually
    ! ====================================================================

    subroutine master_apply_defaults(self, cline)
        class(make_cavgs_master_strategy), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        call apply_distr_entry_defaults(cline)
        call normalize_ncls_nspace_oritype(cline, set_ncls_from_nspace=.false.)
    end subroutine master_apply_defaults

    subroutine master_initialize(self, params, cline)
        class(make_cavgs_master_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        type(builder) :: build_tmp
        integer       :: ncls_here
        logical       :: l_shmem
        ! Master threads (original intent)
        call set_master_num_threads(self%nthr_master, string('CLUSTER2D'))
        ! Parse parameters & project-field to set ncls default for ptcl2D
        call build_tmp%init_params_and_build_spproj(cline, params)
        if( .not. cline%defined('nspace') )then
            ncls_here = build_tmp%spproj_field%get_n('class')
            if( .not. cline%defined('ncls') )then
                call cline%set('ncls', ncls_here)
                params%ncls = ncls_here
            endif
        endif
        call cline%set('mkdir', 'no')  ! avoid nested dirs in workers
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
    end subroutine master_initialize

    subroutine master_execute(self, params, cline)
        class(make_cavgs_master_strategy), intent(inout) :: self
        type(parameters),                  intent(inout) :: params
        class(cmdline),                    intent(inout) :: cline
        type(cmdline) :: cline_cavgassemble
        if( .not. associated(self%hooks%run_cavgassemble) )then
            THROW_HARD('make_cavgs_master_strategy: run_cavgassemble hook not set')
        endif
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        cline_cavgassemble = cline
        if( trim(params%oritype) .eq. 'ptcl3D' )then
            call cline_cavgassemble%set('ncls', params%nspace)
        endif
        call self%hooks%run_cavgassemble(cline_cavgassemble, self%nthr_master)
        call cline_cavgassemble%kill
    end subroutine master_execute

    subroutine master_finalize_run(self, params, cline)
        class(make_cavgs_master_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        ! no-op
    end subroutine master_finalize_run

    subroutine master_cleanup(self, params, cline)
        class(make_cavgs_master_strategy), intent(inout) :: self
        type(parameters),                  intent(in)    :: params
        class(cmdline),                    intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine master_cleanup

    function master_end_message(self) result(msg)
        class(make_cavgs_master_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_MAKE_CAVGS NORMAL STOP ****'
    end function master_end_message

end module simple_make_cavgs_strategy
