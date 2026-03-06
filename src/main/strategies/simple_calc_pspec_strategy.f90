module simple_calc_pspec_strategy
use simple_core_module_api
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none

public :: calc_pspec_strategy, calc_pspec_inmem_strategy, calc_pspec_distr_strategy, create_calc_pspec_strategy
private
#include "simple_local_flags.inc"

!> Abstract strategy interface
type, abstract :: calc_pspec_strategy
contains
    procedure(execute_interface), deferred :: execute
end type calc_pspec_strategy

!> In-memory (shared memory) implementation of the strategy
type, extends(calc_pspec_strategy) :: calc_pspec_inmem_strategy
contains
    procedure :: execute => inmem_execute
end type calc_pspec_inmem_strategy

!> Distributed memory implementation of the strategy
type, extends(calc_pspec_strategy) :: calc_pspec_distr_strategy
contains
    procedure :: execute => distr_execute
end type calc_pspec_distr_strategy

abstract interface
    subroutine execute_interface(self, params, spproj, cline)
        import :: calc_pspec_strategy, parameters, sp_project, cmdline
        class(calc_pspec_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        type(sp_project),           intent(inout) :: spproj
        class(cmdline),             intent(inout) :: cline
    end subroutine execute_interface
end interface

contains

    !> Factory function to create the appropriate strategy based on parameters.
    function create_calc_pspec_strategy(params, cline) result(strategy)
        type(parameters), intent(in) :: params
        class(cmdline),   intent(in) :: cline
        class(calc_pspec_strategy), allocatable :: strategy
        if( (params%nparts > 1) .and. (.not.cline%defined('part')) )then
            allocate(calc_pspec_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(calc_pspec_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY EXECUTION'
        endif
    end function create_calc_pspec_strategy

    ! ========================================================================
    ! SHARED-MEMORY IMPLEMENTATION
    ! ========================================================================
    subroutine inmem_execute(self, params, spproj, cline)
        use simple_calc_pspec_common, only: calc_pspec_exec
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(sp_project),                 intent(inout) :: spproj
        class(cmdline),                   intent(inout) :: cline
        type(builder) :: build
        ! build environment and execute
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        ! call build%init_from_spproj(spproj)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        call calc_pspec_exec(params, build)
        ! cleanup
        call build%kill_general_tbox
    end subroutine inmem_execute

    ! ========================================================================
    ! DISTRIBUTED IMPLEMENTATION
    ! ========================================================================
    subroutine distr_execute(self, params, spproj, cline)
        use simple_qsys_funs,               only: qsys_cleanup
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble, commander_calc_pspec_distr_worker
        class(calc_pspec_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(sp_project),                 intent(inout) :: spproj ! not used, but required by interface
        class(cmdline),                   intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        type(cmdline)        :: cline_calc_pspec_worker
        type(cmdline)        :: cline_calc_pspec_assemble
        type(qsys_env)       :: qenv
        type(chash)          :: job_descr
        ! Prepare command lines from prototype master
        cline_calc_pspec_worker   = cline
        cline_calc_pspec_assemble = cline
        ! The worker program is the specific worker commander
        call cline_calc_pspec_worker%set('prg', 'calc_pspec_distr_worker' )
        call cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble' )
        ! Setup the environment for distributed execution
        call qenv%new(params, params%nparts)
        call cline_calc_pspec_worker%gen_job_descr(job_descr)
        ! Schedule jobs
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! Assemble results from workers
        call xcalc_pspec_assemble%execute(cline_calc_pspec_assemble)
        ! End gracefully
        call cline_calc_pspec_worker%kill
        call cline_calc_pspec_assemble%kill
        call qenv%kill
        call job_descr%kill
        call qsys_cleanup(params)
    end subroutine distr_execute

end module simple_calc_pspec_strategy
