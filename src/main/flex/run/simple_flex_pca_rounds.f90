!@descr: flex_pca distribution contract: the role/round object handed down by the strategy
!!
!! Domain modules (model, em, rec3D) receive a `flex_pca_rounds` and ask it only what a
!! distributable phase needs to know: whether this process is the distributed master or a
!! worker, how many parts a round spans, and how to run one qsys round for a typed
!! `flex_stage_request` (simple_flex_pca_stages). The strategies in strategies/parallelization
!! extend this type; the master implements the rounds on top of its qsys context, shared memory
!! and workers refuse them. Part naming lives in simple_flex_pca_artifacts.
module simple_flex_pca_rounds
use simple_core_module_api
use simple_parameters,      only: parameters
use simple_flex_pca_stages, only: flex_stage_request, FLEX_FIT_ALL
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_rounds, flex_pca_rounds_shmem

type, abstract :: flex_pca_rounds
    logical :: l_master   = .false.
    logical :: l_worker   = .false.
    integer :: nparts_run = 1
    integer :: fit_sel    = FLEX_FIT_ALL   !< the fit every scheduled round serves (two-fit harnesses)
contains
    procedure(plan_iface), deferred :: plan_partitions
    procedure(run_iface),  deferred :: run_stage
    procedure :: is_master   => rounds_is_master
    procedure :: is_worker   => rounds_is_worker
    procedure :: nparts      => rounds_nparts
    procedure :: distributed => rounds_distributed
end type flex_pca_rounds

!> Shared memory and workers: no partitions, no rounds. Workers set l_worker.
type, extends(flex_pca_rounds) :: flex_pca_rounds_shmem
contains
    procedure :: plan_partitions => shmem_plan_partitions
    procedure :: run_stage       => shmem_run_stage
end type flex_pca_rounds_shmem

abstract interface
    !> Partition the master's particle selection into one list per part (no-op unless distributed)
    subroutine plan_iface( self, params, pinds )
        import :: flex_pca_rounds, parameters
        class(flex_pca_rounds), intent(inout) :: self
        type(parameters),       intent(inout) :: params
        integer,                intent(in)    :: pinds(:)
    end subroutine plan_iface
    !> One qsys round: every worker runs the requested stage over its own particle list and writes
    !! its part file; returns when all parts are on disk.
    subroutine run_iface( self, params, req )
        import :: flex_pca_rounds, parameters, flex_stage_request
        class(flex_pca_rounds),   intent(inout) :: self
        type(parameters),         intent(in)    :: params
        type(flex_stage_request), intent(in)    :: req
    end subroutine run_iface
end interface

contains

    logical function rounds_is_master( self )
        class(flex_pca_rounds), intent(in) :: self
        rounds_is_master = self%l_master
    end function rounds_is_master

    logical function rounds_is_worker( self )
        class(flex_pca_rounds), intent(in) :: self
        rounds_is_worker = self%l_worker
    end function rounds_is_worker

    integer function rounds_nparts( self )
        class(flex_pca_rounds), intent(in) :: self
        rounds_nparts = self%nparts_run
    end function rounds_nparts

    !> A distributable phase fans out iff this is the master of a multi-part run
    logical function rounds_distributed( self )
        class(flex_pca_rounds), intent(in) :: self
        rounds_distributed = self%l_master .and. self%nparts_run > 1
    end function rounds_distributed

    subroutine shmem_plan_partitions( self, params, pinds )
        class(flex_pca_rounds_shmem), intent(inout) :: self
        type(parameters),             intent(inout) :: params
        integer,                      intent(in)    :: pinds(:)
    end subroutine shmem_plan_partitions

    subroutine shmem_run_stage( self, params, req )
        class(flex_pca_rounds_shmem), intent(inout) :: self
        type(parameters),             intent(in)    :: params
        type(flex_stage_request),     intent(in)    :: req
        THROW_HARD('flex_pca round requested outside the distributed master: '//trim(req%label))
    end subroutine shmem_run_stage

end module simple_flex_pca_rounds
