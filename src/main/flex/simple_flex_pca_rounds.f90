!@descr: flex_pca distribution contract: the role/round object handed down by the strategy
!!
!! Domain modules (model, em, rec3D) receive a `flex_pca_rounds` and ask it only what a
!! distributable phase needs to know: whether this process is the distributed master or a
!! worker, how many parts a round spans, and how to run one qsys round. The strategies in
!! strategies/parallelization extend this type; the master implements the rounds on top of its
!! qsys context, shared memory and workers refuse them. This leaf module also holds the stage
!! and fit identifiers that travel to workers through job_descr, and the ONE mod-4 half rule
!! the paired engine and the two-job half harness share, so both partition identically.
module simple_flex_pca_rounds
use simple_core_module_api
use simple_parameters, only: parameters
implicit none
private
#include "simple_local_flags.inc"

!> part-directory prefix ('<dir>/'), unallocated = the run directory (see flex_pca_part_path)
character(len=:), allocatable :: part_dir_prefix

public :: flex_pca_rounds, flex_pca_rounds_shmem
public :: PCA_STAGE_EMBED, PCA_STAGE_STATES, PCA_STAGE_PROBE, PCA_STAGE_POLISH
public :: FLEX_FIT_ALL, FLEX_FIT_A, FLEX_FIT_B
public :: flex_pca_half_of, flex_pca_part_fname, flex_pca_part_path
public :: flex_pca_set_part_dir, flex_pca_local_part_dir
public :: FLEX_PCA_PART_MAGIC, PROBE_PART_VERSION, PROBE_PART_VERSION5, EMBED_STATS_VERSION

! Stage selector carried to the worker in params%stage. 1..3 were the moment estimator's SNR /
! column / reduced-solve rounds; that path is gone and the numbering is kept so an old part file
! cannot be mistaken for a current one.
integer, parameter :: PCA_STAGE_EMBED  = 4
integer, parameter :: PCA_STAGE_STATES = 5
! One qsys round per probe EM iteration: the basis changes every iteration, so workers are
! re-launched against the master's refreshed flex_pca_pc*.mrc rather than looping locally.
integer, parameter :: PCA_STAGE_PROBE  = 6
! The joint (polish) fit after the paired merge: the same E-step pass over ALL particles against
! the merged basis under the polished namespace. A stage, not a string stamp: the stage says what
! to compute and which basis namespace to load.
integer, parameter :: PCA_STAGE_POLISH = 7

! Which fit a round serves, carried to the worker in params%pcafit beside the stage: the stage
! says what to compute, this says over which particles.
integer, parameter :: FLEX_FIT_ALL = 0
integer, parameter :: FLEX_FIT_A   = 1
integer, parameter :: FLEX_FIT_B   = 2

type, abstract :: flex_pca_rounds
    logical :: l_master   = .false.
    logical :: l_worker   = .false.
    integer :: nparts_run = 1
    integer :: fit_sel    = FLEX_FIT_ALL
contains
    procedure(plan_iface), deferred :: plan_partitions
    procedure(run_iface),  deferred :: run_stage
    procedure :: is_master   => rounds_is_master
    procedure :: is_worker   => rounds_is_worker
    procedure :: nparts      => rounds_nparts
    procedure :: distributed => rounds_distributed
    procedure :: set_fit     => rounds_set_fit
    procedure :: fit         => rounds_fit
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
    !> One qsys round: every worker runs `stage` over its own particle list and writes its part
    !! file; returns when all parts are on disk. which_iter/maxits/nfits are the round's control
    !! state, carried to the workers in job_descr under registered keys.
    subroutine run_iface( self, params, stage_id, label, which_iter, maxits, nfits )
        import :: flex_pca_rounds, parameters
        class(flex_pca_rounds), intent(inout) :: self
        type(parameters),       intent(in)    :: params
        integer,                intent(in)    :: stage_id
        character(len=*),       intent(in)    :: label
        integer, optional,      intent(in)    :: which_iter, maxits, nfits
    end subroutine run_iface
end interface

! ---- part-file contract: every part file is magic + version + shape header + payload, written
! to a .tmp and renamed, so a master that finds the final name is guaranteed a complete file ----
integer, parameter :: EMBED_STATS_VERSION = 1
integer, parameter :: PROBE_PART_VERSION  = 12  ! rho rows are always the full packed triangle; trailing PCG kernel + rhs blocks on the shared band list (slot for slot, no per-part index list)
integer, parameter :: PROBE_PART_VERSION5 = 10  ! v5 layout, payloads band-boxed (nonzero bounding box per lattice) + trailing PCG kernel + rhs blocks per fit
integer, parameter :: FLEX_PCA_PART_MAGIC = 1180053590

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

    !> Master-side: stamp every subsequent round as serving this fit (two-fit harnesses)
    subroutine rounds_set_fit( self, ifit )
        class(flex_pca_rounds), intent(inout) :: self
        integer,                intent(in)    :: ifit
        self%fit_sel = ifit
    end subroutine rounds_set_fit

    integer function rounds_fit( self )
        class(flex_pca_rounds), intent(in) :: self
        rounds_fit = self%fit_sel
    end function rounds_fit

    subroutine shmem_plan_partitions( self, params, pinds )
        class(flex_pca_rounds_shmem), intent(inout) :: self
        type(parameters),             intent(inout) :: params
        integer,                      intent(in)    :: pinds(:)
    end subroutine shmem_plan_partitions

    subroutine shmem_run_stage( self, params, stage_id, label, which_iter, maxits, nfits )
        class(flex_pca_rounds_shmem), intent(inout) :: self
        type(parameters),             intent(in)    :: params
        integer,                      intent(in)    :: stage_id
        character(len=*),             intent(in)    :: label
        integer, optional,            intent(in)    :: which_iter, maxits, nfits
        THROW_HARD('flex_pca round requested outside the distributed master: '//trim(label))
    end subroutine shmem_run_stage

    !> The ONE mod-4 split rule, shared by the two-job pcafit harness (validate_covariance_inputs)
    !! and the paired engine's driver -- so the two instruments partition the selection identically
    !! by construction. Pairing 1 (default) puts row residues {0,1} in half A; pairing 3 puts
    !! {0,3}. Both pair one even-row residue with one odd-row residue, so each half keeps both
    !! internal e/o classes under the row-alternating project eo split. Pairing 2 ({0,2}|{1,3})
    !! is eo-degenerate by construction and is REFUSED where SIMPLE_COV_MOD4_PAIRING is read --
    !! never silently mapped here.
    pure integer function flex_pca_half_of( pind, vpair ) result( ifit )
        integer, intent(in) :: pind, vpair
        integer :: r4, a2
        a2 = 1
        if( vpair == 3 ) a2 = 3
        r4 = mod(pind, 4)
        if( r4 == 0 .or. r4 == a2 )then
            ifit = FLEX_FIT_A
        else
            ifit = FLEX_FIT_B
        endif
    end function flex_pca_half_of

    function flex_pca_part_fname( prefix, part, numlen ) result( fname )
        character(len=*), intent(in) :: prefix
        integer,          intent(in) :: part, numlen
        type(string) :: fname
        fname = flex_pca_part_path('flex_pca_'//prefix//'_part'//int2str_pad(part, numlen)//'.bin')
    end function flex_pca_part_fname

    !> Every part file (probe/embed parts, state part volumes) lives under one directory: the
    !! run directory by default, or a node-local scratch directory when the master decided so
    !! (flex_pca_local_part_dir). Producers and consumers name parts only through this path.
    function flex_pca_part_path( name ) result( path )
        character(len=*), intent(in) :: name
        type(string) :: path
        if( allocated(part_dir_prefix) )then
            path = string(part_dir_prefix)//name
        else
            path = string(name)
        endif
    end function flex_pca_part_path

    subroutine flex_pca_set_part_dir( dir )
        character(len=*), intent(in) :: dir
        if( len_trim(dir) == 0 )then
            if( allocated(part_dir_prefix) ) deallocate(part_dir_prefix)
        else
            part_dir_prefix = trim(dir)//'/'
        endif
    end subroutine flex_pca_set_part_dir

    !> Node-local part directory for the local queue system: parts are written and reduced once
    !! per round (5-54 s per iteration and ~90 s in the states stage over the network on
    !! 2026-09-08), so they go to the disk the user already declared local through cache_dir.
    !! Master and workers derive the SAME name from the run directory (no new key travels), which
    !! is safe because local workers run on the master's node in the master's directory. Empty
    !! (= run directory) for any other queue system, when cache_dir is not given, or in shared
    !! memory.
    function flex_pca_local_part_dir( params, nparts ) result( dir )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nparts
        character(len=:), allocatable :: dir
        type(string) :: cwd
        character(len=:), allocatable :: cwdc
        integer(kind=8) :: h
        integer :: i
        dir = ''
        if( nparts < 2 ) return
        if( trim(params%qsys_name) /= 'local' ) return
        if( params%cache_dir%is_blank() ) return
        call simple_getcwd(cwd)
        cwdc = cwd%to_char()
        h = 1469598103934665603_8            ! FNV-1a over the run directory path
        do i = 1, len_trim(cwdc)
            h = ieor(h, int(ichar(cwdc(i:i)),8))
            h = h * 1099511628211_8
        end do
        h = iand(h, 9223372036854775807_8)
        dir = params%cache_dir%to_char()//'/flex_pca_parts_'//trim(int2str(int(mod(h, 1000000007_8))))
        call cwd%kill
    end function flex_pca_local_part_dir

end module simple_flex_pca_rounds
