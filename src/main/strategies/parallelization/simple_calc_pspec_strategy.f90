module simple_calc_pspec_strategy
use simple_core_module_api
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup, qsys_job_finished
use simple_image,          only: image
use simple_sigma2_binfile, only: sigma2_binfile
implicit none

public :: calc_pspec_strategy, calc_pspec_inmem_strategy, calc_pspec_distr_strategy, create_calc_pspec_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: calc_pspec_strategy
    character(len=:), allocatable :: end_msg
    logical :: end_print_simple = .true.
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
end type calc_pspec_strategy

type, extends(calc_pspec_strategy) :: calc_pspec_inmem_strategy
    type(cmdline) :: cline_calc_pspec_assemble
contains
    procedure :: initialize   => inmem_initialize
    procedure :: execute      => inmem_execute
    procedure :: finalize_run => inmem_finalize_run
    procedure :: cleanup      => inmem_cleanup
end type calc_pspec_inmem_strategy

type, extends(calc_pspec_strategy) :: calc_pspec_distr_strategy
    type(qsys_env) :: qenv
    type(chash)    :: job_descr
    type(cmdline)  :: cline_calc_pspec
    type(cmdline)  :: cline_calc_pspec_assemble
contains
    procedure :: initialize   => distr_initialize
    procedure :: execute      => distr_execute
    procedure :: finalize_run => distr_finalize_run
    procedure :: cleanup      => distr_cleanup
end type calc_pspec_distr_strategy

abstract interface
    subroutine init_interface(self, params, build, cline)
        import :: calc_pspec_strategy, parameters, builder, cmdline
        class(calc_pspec_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        type(builder),              intent(inout) :: build
        class(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, build, cline)
        import :: calc_pspec_strategy, parameters, builder, cmdline
        class(calc_pspec_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        type(builder),              intent(inout) :: build
        class(cmdline),             intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, build, cline)
        import :: calc_pspec_strategy, parameters, builder, cmdline
        class(calc_pspec_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        type(builder),              intent(inout) :: build
        class(cmdline),             intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, build, cline)
        import :: calc_pspec_strategy, parameters, builder, cmdline
        class(calc_pspec_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        type(builder),              intent(inout) :: build
        class(cmdline),             intent(inout) :: cline
    end subroutine cleanup_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------

    function create_calc_pspec_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(calc_pspec_strategy), allocatable :: strategy
        if( cline%defined('nparts') .and. (.not. cline%defined('part')) )then
            allocate(calc_pspec_distr_strategy :: strategy)
            strategy%end_msg          = '**** SIMPLE_DISTR_CALC_PSPEC NORMAL STOP ****'
            strategy%end_print_simple = .true.
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED EXECUTION'
        else
            allocate(calc_pspec_inmem_strategy :: strategy)
            strategy%end_msg          = '**** SIMPLE_CALC_PSPEC NORMAL STOP ****'
            strategy%end_print_simple = .false.
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY/WORKER EXECUTION'
        endif
    end function create_calc_pspec_strategy

    ! =====================================================================
    ! SHARED-MEMORY / DISTRIBUTED WORKER IMPLEMENTATION
    ! =====================================================================

    subroutine inmem_initialize(self, params, build, cline)
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        ! Worker/shared-memory execution must not create nested directories
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        ! Parse params + build general toolbox (2D)
        call params%new(cline)
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        self%cline_calc_pspec_assemble = cline
        call self%cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble')
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, build, cline)
        use simple_strategy2D3D_common,     only: prepimgbatch, discrete_read_imgbatch
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        type(image)                         :: sum_img
        type(sigma2_binfile)                :: binfile
        type(string)                        :: binfname
        complex(dp), allocatable :: cmat_thr_sum(:,:,:)
        complex,     allocatable :: cmat_sum(:,:,:)
        integer,     allocatable :: pinds(:)
        real,        allocatable :: sigma2(:,:), sigma2_batch(:,:)
        integer :: batchlims(2), kfromto(2)
        integer :: i, iptcl, imatch, nyq, nptcls_part_sel, batchsz_max, nbatch
        logical :: l_scale_update_frac
        real    :: sig2_mul
        ! Sampling
        ! Because this is always run prior to reconstruction/search, sampling is not always informed
        ! or may change with workflows. Instead of setting a sampling for the following operations when
        ! l_update_frac, we sample uniformly AND do not write the corresponding field
        l_scale_update_frac = .false.
        if( params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], params%update_frac, nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
            sig2_mul = 1.0 / (2.0 * params%update_frac)
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls_part_sel, pinds, .false.)
            sig2_mul = 0.5
        endif
        ! Init
        nyq = build%img%get_nyq()
        allocate(sigma2(nyq,params%fromp:params%top), source=0.)
        batchsz_max = min(nptcls_part_sel, 50 * nthr_glob)
        call prepimgbatch(params, build, batchsz_max)
        call sum_img%new([params%box,params%box,1], params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_sum = sum_img%allocate_cmat()
        allocate(cmat_thr_sum(size(cmat_sum,dim=1), size(cmat_sum,dim=2), 1))
        allocate(sigma2_batch(nyq,batchsz_max), source=0.)
        ! mask memoization
        call build%imgbatch(1)%memoize_mask_coords
        do i = 1, nptcls_part_sel, batchsz_max
            batchlims = [i, min(i+batchsz_max-1, nptcls_part_sel)]
            nbatch    = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nbatch, pinds(batchlims(1):batchlims(2)), [1,nbatch])
            ! allocate contigous local sigma2 array for optimal caching in parallell loop
            cmat_thr_sum = dcmplx(0.d0, 0.d0)
            !$omp parallel do default(shared) private(iptcl,imatch)&
            !$omp schedule(static) proc_bind(close) reduction(+:cmat_thr_sum)
             do imatch = 1, nbatch
                iptcl = pinds(batchlims(1) + imatch - 1)
                call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, sigma2_batch(:,imatch))
                sigma2_batch(:,imatch) = sigma2_batch(:,imatch) * sig2_mul
                ! thread average
                call build%imgbatch(imatch)%add_dble_cmat2mat(cmat_thr_sum(:,:,:))
            end do
            !$omp end parallel do
            ! global average
            cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmplx(cmat_thr_sum(:,:,:), kind=sp)
            ! update non-contiguous sigma2 array to provide the correct geometry on disk
            do imatch = 1, nbatch
                iptcl = pinds(batchlims(1) + imatch - 1)
                sigma2(:,iptcl) = sigma2_batch(:,imatch)
            end do
        end do
        deallocate(sigma2_batch)
        call sum_img%set_cmat(cmat_sum)
        call sum_img%write(string('sum_img_part')//int2str_pad(params%part,params%numlen)//params%ext%to_char())
        ! write sigma2 to disk
        kfromto  = [1, nyq]
        binfname = 'init_pspec_part'//trim(int2str(params%part))//'.dat'
        call binfile%new(binfname, params%fromp, params%top, kfromto)
        call binfile%write(sigma2)
        ! destruct local objects
        call binfile%kill
        call sum_img%kill
        ! Group averaging in non-distributed execution
        if( .not.cline%defined('nparts') .and. params%part==1 )then
            call xcalc_pspec_assemble%execute(self%cline_calc_pspec_assemble)
        endif
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, build, cline)
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_euclid :: exec_calc_pspec'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, build, cline)
        use simple_strategy2D3D_common, only: killimgbatch
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%cline_calc_pspec_assemble%kill
        call build%kill_general_tbox
        call killimgbatch(build)
    end subroutine inmem_cleanup

    ! =====================================================================
    ! DISTRIBUTED MASTER IMPLEMENTATION
    ! =====================================================================

    subroutine distr_initialize(self, params, build, cline)
        class(calc_pspec_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call cline%set('stream', 'no')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec_distr')
        endif
        ! Parse parameters (mkdir default kept as provided/assigned by commander-level defaults)
        call params%new(cline)
        ! Build project ONLY (no general toolbox)
        call build%build_spproj(params, cline)
        ! Sanity check + ensure EO partition exists and is persisted
        call sanity_check_calc_pspec_input(params, build)
        call ensure_calc_pspec_eo_partition(params, build)
        ! Set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! Prepare command lines from prototype master
        self%cline_calc_pspec          = cline
        self%cline_calc_pspec_assemble = cline
        ! Required program keys
        call self%cline_calc_pspec%set('prg', 'calc_pspec')
        call self%cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble')
        ! Setup the environment for distributed execution
        call self%qenv%new(params, params%nparts)
        call self%cline_calc_pspec%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, build, cline)
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(calc_pspec_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        ! schedule
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! assemble (local)
        call xcalc_pspec_assemble%execute(self%cline_calc_pspec_assemble)
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, build, cline)
        class(calc_pspec_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, build, cline)
        class(calc_pspec_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%cline_calc_pspec%kill
        call self%cline_calc_pspec_assemble%kill
        call self%qenv%kill
        call self%job_descr%kill
        call qsys_cleanup(params)
        ! Best-effort cleanup of project-related resources
        call build%spproj_field%kill
        call build%spproj%kill
        call simple_touch(CALCPSPEC_FINISHED)
    end subroutine distr_cleanup

    ! =====================================================================
    ! Internal helpers
    ! =====================================================================

    subroutine sanity_check_calc_pspec_input(params, build)
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        logical :: fall_over
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl2D','ptcl3D','cls3D')
                fall_over = build%spproj%get_nptcls() == 0
            case DEFAULT
                THROW_HARD('Unsupported ORITYPE; simple_commanders_euclid :: exec_calc_pspec_distr')
        end select
        if( fall_over )then
            THROW_HARD('no particles found! :exec_calc_pspec_distr')
        endif
    end subroutine sanity_check_calc_pspec_input

    subroutine ensure_calc_pspec_eo_partition(params, build)
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
    end subroutine ensure_calc_pspec_eo_partition

end module simple_calc_pspec_strategy