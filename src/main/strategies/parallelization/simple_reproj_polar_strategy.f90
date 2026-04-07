!@descr: strategy for distributed reprojection of 3D reference volumes to polar FTs
!        Worker: reads, masks and filters each state volume; calls vol_pad2ref_pfts_write_range
!                to project and write the assigned [fromp,top] range of orientations for
!                every state into per-state/per-part binary files.
!        Master: assembles the per-part files back into a full pftc ready for alignment.
module simple_reproj_polar_strategy
use simple_core_module_api
use simple_builder,              only: builder
use simple_parameters,           only: parameters
use simple_cmdline,              only: cmdline
use simple_qsys_env,             only: qsys_env
use simple_exec_helpers,         only: set_master_num_threads
use simple_polarft_calc,         only: vol_pad2ref_pfts, vol_pad2ref_pfts_write_range
use simple_reproj_refvol_utils,  only: read_mask_filter_refvols, calcrefvolshift_and_mapshifts2ptcls
implicit none

public  :: reproj_polar_strategy, reproject_inmem_strategy, reproject_distr_strategy, create_reproj_polar_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: reproj_polar_strategy
contains
    procedure(init_interface),    deferred :: initialize
    procedure(exec_interface),    deferred :: execute
    procedure(cleanup_interface), deferred :: cleanup
end type reproj_polar_strategy

! Worker / shared-memory (runs one part)
type, extends(reproj_polar_strategy) :: reproject_inmem_strategy
contains
    procedure :: initialize => inmem_initialize
    procedure :: execute    => inmem_execute
    procedure :: cleanup    => inmem_cleanup
end type reproject_inmem_strategy

! Distributed master (schedules workers, then assembles)
type, extends(reproj_polar_strategy) :: reproject_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    integer        :: nthr_master = 1
contains
    procedure :: initialize => distr_initialize
    procedure :: execute    => distr_execute
    procedure :: cleanup    => distr_cleanup
end type reproject_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: reproj_polar_strategy, parameters, builder, cmdline
        class(reproj_polar_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: reproj_polar_strategy, parameters, builder, cmdline
        class(reproj_polar_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine exec_interface

    subroutine cleanup_interface(self, params, build, cline)
        import :: reproj_polar_strategy, parameters, builder, cmdline
        class(reproj_polar_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        type(builder),             intent(inout) :: build
        class(cmdline),            intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Strategy selection
    ! --------------------------------------------------------------------

    function create_reproj_polar_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(reproj_polar_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            allocate(reproject_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED REPROJECT EXECUTION'
        else
            allocate(reproject_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY REPROJECT EXECUTION'
        endif
    end function create_reproj_polar_strategy

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
        if( cline%defined('part') )then
            call inmem_execute_worker(params, build)
            call qsys_job_finished(params, string('simple_reproj_polar_strategy :: reproject_exec'))
        else
            call inmem_execute_full(params, build)
        endif
    end subroutine inmem_execute

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
        ! Parse parameters and build the in-memory project representation.
        ! Unlike refinement/reconstruction workflows, reproj_polar can operate
        ! directly from volumes and does not inherently require an on-disk
        ! project file.
        call params%new(cline)
        call build%build_spproj(params, cline)
        call distr_prepare_runtime(self, params, cline)
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
    end subroutine distr_execute

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

    subroutine inmem_execute_worker( params, build )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        type(string) :: part_tmpl
        real         :: xyz(3)
        integer      :: s, nrefs
        logical      :: do_center
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1, 1], params%kfromto)
        do s = 1, params%nstates
            call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                & do_center, xyz, map_shift=.false.)
            call read_mask_filter_refvols(params, build, s)
            call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol%fft()
                call build%vol%shift(xyz)
            endif
            call build%vol%ifft()
            call build%vol%pad_fft(build%vol_pad)
            call build%vol_pad%expand_cmat(params%box)
            part_tmpl = string(POLAR_REFS_FBODY)//'_s'//int2str_pad(s, 2) &
                & //'_part'//int2str_pad(params%part, params%numlen)
            call vol_pad2ref_pfts_write_range(build%pftc, build%vol_pad, build%eulspace, &
                & s, params%fromp, params%top, build%l_resmsk, part_tmpl)
            call build%vol_pad%kill
            call build%vol_pad%kill_expanded
            call part_tmpl%kill
        end do
        call build%pftc%kill
    end subroutine inmem_execute_worker

    subroutine inmem_execute_full( params, build )
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        real    :: xyz(3)
        integer :: s, nrefs
        logical :: do_center, l_prob_align_mode
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1, 1], params%kfromto)
        select case(trim(params%refine))
            case('prob','prob_state','prob_neigh')
                l_prob_align_mode = .true.
            case DEFAULT
                l_prob_align_mode = .false.
        end select
        do s = 1, params%nstates
            if( l_prob_align_mode )then
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                    & do_center, xyz, map_shift=l_distr_worker_glob)
            else
                call calcrefvolshift_and_mapshifts2ptcls(params, build, s, params%vols(s), &
                    & do_center, xyz, map_shift=.true.)
            endif
            call read_mask_filter_refvols(params, build, s)
            call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol%fft()
                call build%vol%shift(xyz)
            endif
            call build%vol%ifft()
            call build%vol%pad_fft(build%vol_pad)
            call build%vol_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_pad, build%eulspace, s, .true., build%l_resmsk)
            call build%vol_pad%kill
            call build%vol_pad%kill_expanded
            call build%vol_odd_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                & params%smpd_crop, wthreads=.true.)
            if( do_center )then
                call build%vol_odd%fft()
                call build%vol_odd%shift(xyz)
            endif
            call build%vol_odd%ifft()
            call build%vol_odd%pad_fft(build%vol_odd_pad)
            call build%vol_odd_pad%expand_cmat(params%box)
            call vol_pad2ref_pfts(build%pftc, build%vol_odd_pad, build%eulspace, s, .false., build%l_resmsk)
            call build%vol_odd_pad%kill
            call build%vol_odd_pad%kill_expanded
        end do
    end subroutine inmem_execute_full

    subroutine distr_prepare_runtime( self, params, cline )
        class(reproject_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
        call set_master_num_threads(self%nthr_master, string('reproject'))
        call cline%set('nspace', params%nspace)
        if( .not. cline%defined('top') ) call cline%set('top', params%nspace)
        call cline%set('mkdir', 'no')
        call self%qenv%new(params, params%nparts, nptcls=params%nspace)
        call cline%gen_job_descr(self%job_descr, string('reproj_polar'))
        call self%job_descr%delete('test')
        call self%job_descr%move_key_to_front('prg')
    end subroutine distr_prepare_runtime

end module simple_reproj_polar_strategy
