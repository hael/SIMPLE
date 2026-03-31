!@descr: strategy for distributed reprojection of 3D reference volumes to polar FTs
!        Worker: reads, masks and filters each state volume; calls vol_pad2ref_pfts_write_range
!                to project and write the assigned [fromp,top] range of orientations for
!                every state into per-state/per-part binary files.
!        Master: assembles the per-part files back into a full pftc ready for alignment.
module simple_reproject_strategy
use simple_core_module_api
use simple_builder,              only: builder
use simple_parameters,           only: parameters
use simple_cmdline,              only: cmdline
use simple_qsys_env,             only: qsys_env
use simple_exec_helpers,         only: set_master_num_threads
use simple_polarft_calc,         only: vol_pad2ref_pfts_write_range
use simple_strategy2D3D_common,  only: read_mask_filter_refvols, calcrefvolshift_and_mapshifts2ptcls
implicit none

public  :: reproject_strategy, reproject_inmem_strategy, reproject_distr_strategy, create_reproject_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: reproject_strategy
contains
    procedure(init_interface),      deferred :: initialize
    procedure(exec_interface),      deferred :: execute
    procedure(finalize_interface),  deferred :: finalize_run
    procedure(cleanup_interface),   deferred :: cleanup
end type reproject_strategy

! Worker / shared-memory (runs one part)
type, extends(reproject_strategy) :: reproject_inmem_strategy
contains
    procedure :: initialize   => inmem_initialize
    procedure :: execute      => inmem_execute
    procedure :: finalize_run => inmem_finalize_run
    procedure :: cleanup      => inmem_cleanup
end type reproject_inmem_strategy

! Distributed master (schedules workers, then assembles)
type, extends(reproject_strategy) :: reproject_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    integer        :: nthr_master = 1
contains
    procedure :: initialize   => distr_initialize
    procedure :: execute      => distr_execute
    procedure :: finalize_run => distr_finalize_run
    procedure :: cleanup      => distr_cleanup
end type reproject_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: reproject_strategy, parameters, builder, cmdline
        class(reproject_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: reproject_strategy, parameters, builder, cmdline
        class(reproject_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: reproject_strategy, parameters, builder, cmdline
        class(reproject_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, build, cline)
        import :: reproject_strategy, parameters, builder, cmdline
        class(reproject_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Strategy selection
    ! --------------------------------------------------------------------

    function create_reproject_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(reproject_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            allocate(reproject_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED REPROJECT EXECUTION'
        else
            allocate(reproject_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY REPROJECT EXECUTION'
        endif
    end function create_reproject_strategy

    ! =====================================================================
    ! WORKER (inmem) IMPLEMENTATION
    ! =====================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(reproject_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        ! Build the general toolbox: volumes, oris, pftc geometry etc.
        call build%init_params_and_build_general_tbox(cline, params)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, build, cline)
        use simple_qsys_funs, only: qsys_job_finished
        class(reproject_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        type(string) :: part_tmpl
        real         :: xyz(3)
        integer      :: s, nrefs
        logical      :: do_center
        ! Build a full pftc: nrefs covers all states × all orientations so that
        ! vol_pad2ref_pfts_write_range can address the correct global ref positions.
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1, 1], params%kfromto)
        do s = 1, params%nstates
            ! Optional centering and particle-shift mapping.  vol_pad is not yet allocated
            ! here so calcrefvolshift_and_mapshifts2ptcls reads the volume internally.
            call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                & do_center, xyz, map_shift=.false.)
            ! Read, mask, and filter the average volume into build%vol (Fourier space on return).
            ! For the reproject case we use a single filtered volume for both E/O halves,
            ! equivalent to the l_lpset path in read_mask_filter_reproject_refvols.
            call read_mask_filter_refvols(params, build, s)
            ! Prepare even padded volume
            call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol%fft()
                call build%vol%shift(xyz)
            endif
            call build%vol%ifft()
            call build%vol%pad_fft(build%vol_pad)
            call build%vol_pad%expand_cmat(params%box)
            ! Project the assigned orientation range and write part files for this state.
            ! The tmpl produces: polar_refs_s{ss}_part{NNN}_even.bin and _odd.bin
            part_tmpl = string(POLAR_REFS_FBODY)//'_s'//int2str_pad(s, 2) &
                & //'_part'//int2str_pad(params%part, params%numlen)
            call vol_pad2ref_pfts_write_range(build%pftc, build%vol_pad, build%eulspace, &
                & s, params%fromp, params%top, build%l_resmsk, part_tmpl)
            call build%vol_pad%kill
            call build%vol_pad%kill_expanded
            call part_tmpl%kill
        end do
        call qsys_job_finished(params, string('simple_reproject_strategy :: reproject_exec'))
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, build, cline)
        class(reproject_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        ! No-op
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, build, cline)
        class(reproject_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        ! Nothing to clean up beyond what the builder owns.
    end subroutine inmem_cleanup

    ! =====================================================================
    ! DISTRIBUTED MASTER IMPLEMENTATION
    ! =====================================================================

    subroutine distr_initialize(self, params, build, cline)
        class(reproject_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call set_master_num_threads(self%nthr_master, string('reproject'))
        ! Parse parameters and read the project (no particle stacking needed).
        call build%init_params_and_build_spproj(cline, params)
        ! Avoid nested directory structure for job scripts.
        call cline%set('mkdir', 'no')
        ! Partition nspace orientations across nparts workers: fromp/top will be
        ! projection-direction indices (1..nspace), not particle indices.
        call self%qenv%new(params, params%nparts, nptcls=params%nspace)
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, build, cline)
        class(reproject_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        integer :: nrefs
        ! -- WORKER PHASE --
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! -- ASSEMBLY PHASE --
        ! Build a minimal pftc (no particles) sized for all states × all orientations.
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1, 1], params%kfromto)
        ! Read each worker's per-state binary slices and place them at the correct
        ! global reference positions in pfts_refs_even/odd.
        call build%pftc%assemble_projected_refs_from_parts(params%nparts, params%numlen)
        ! Memoize the assembled references for immediate use in alignment.
        call build%pftc%memoize_refs
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, build, cline)
        class(reproject_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, build, cline)
        use simple_qsys_funs, only: qsys_cleanup
        class(reproject_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        type(builder),                   intent(inout) :: build
        class(cmdline),                  intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

end module simple_reproject_strategy
