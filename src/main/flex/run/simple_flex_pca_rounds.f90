!@descr: flex_pca distribution contract: the role/round object handed down by the strategy
!! Domain modules (model, em, rec3D) ask a flex_pca_rounds only for the master/worker role, the part
!! count and how to run one qsys round for a typed flex_stage_request (simple_flex_pca_stages). The
!! strategies extend it: the master implements rounds on its qsys context, shared memory and workers
!! refuse them. Part naming lives in simple_flex_pca_artifacts.
module simple_flex_pca_rounds
use simple_core_module_api,    only: simple_exception, string
use simple_parameters,         only: parameters
use simple_flex_pca_stages,    only: flex_stage_request, FLEX_FIT_ALL
use simple_flex_pca_artifacts, only: flex_pca_artifact_catalog
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_rounds, flex_pca_rounds_shmem

type, abstract :: flex_pca_rounds
    logical :: l_master   = .false.
    logical :: l_worker   = .false.
    integer :: nparts_run = 1
    integer :: fit_sel    = FLEX_FIT_ALL   !< the fit every scheduled round serves (two-fit harnesses)
    type(flex_pca_artifact_catalog), private :: artifacts
  contains
    procedure(plan_iface), deferred :: plan_partitions
    procedure(run_iface),  deferred :: run_stage
    procedure :: is_master   => rounds_is_master
    procedure :: is_worker   => rounds_is_worker
    procedure :: nparts      => rounds_nparts
    procedure :: distributed => rounds_distributed
    procedure :: part_fname  => rounds_part_fname
    procedure :: part_path   => rounds_part_path
    procedure :: set_part_dir => rounds_set_part_dir
    procedure :: kill_artifacts => rounds_kill_artifacts
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

    function rounds_part_fname( self, prefix, part, numlen ) result( fname )
        class(flex_pca_rounds), intent(in) :: self
        character(len=*),       intent(in) :: prefix
        integer,                intent(in) :: part, numlen
        type(string) :: fname
        fname = self%artifacts%part_fname(prefix, part, numlen)
    end function rounds_part_fname

    function rounds_part_path( self, name ) result( path )
        class(flex_pca_rounds), intent(in) :: self
        character(len=*),       intent(in) :: name
        type(string) :: path
        path = self%artifacts%part_path(name)
    end function rounds_part_path

    subroutine rounds_set_part_dir( self, dir )
        class(flex_pca_rounds), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        call self%artifacts%new(dir)
    end subroutine rounds_set_part_dir

    subroutine rounds_kill_artifacts( self )
        class(flex_pca_rounds), intent(inout) :: self
        call self%artifacts%kill
    end subroutine rounds_kill_artifacts

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
