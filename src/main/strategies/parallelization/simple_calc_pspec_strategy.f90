module simple_calc_pspec_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup, qsys_job_finished
use simple_image,          only: image, unmemoize_powspec_coords
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
        use simple_matcher_ptcl_io,         only: prepimgbatch, discrete_read_imgbatch, discrete_read_imgbatch_source
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        type(image)                         :: sum_img
        type(sigma2_binfile)                :: binfile
        type(string)                        :: binfname
        complex(dp), allocatable :: cmat_thr_sum(:,:,:,:)
        complex,     allocatable :: cmat_sum(:,:,:)
        integer,     allocatable :: pinds(:)
        real,        allocatable :: sigma2(:,:), sigma2_batch(:,:)
        integer :: batchlims(2), kfromto(2)
        integer :: i, iptcl, imatch, ithr, nyq, nptcls_part_sel, nptcls_active_tot, batchsz_max, nbatch
        integer :: ninvalid_sigma2
        logical :: l_scale_update_frac
        logical :: l_add_to_sum
        real    :: sig2_mul, update_frac_eff
        ! Sampling
        ! Because this is always run prior to reconstruction/search, sampling is not always informed
        ! or may change with workflows. Instead of setting a sampling for the following operations when
        ! l_update_frac, we sample uniformly AND do not write the corresponding field. Group sigmas need
        ! direct spectra for every active particle in each stack group; only the global estimate uses
        ! the sampled/bootstrap path.
        l_scale_update_frac = .false.
        if( params%l_sigma_glob .and. params%l_update_frac )then
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], params%update_frac, &
                nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
            sig2_mul = 1.0 / (2.0 * params%update_frac)
        else if( params%l_sigma_glob .and. params%nsample > 0 )then
            nptcls_active_tot = build%spproj%count_state_gt_zero()
            update_frac_eff = min(1.0, real(params%nsample) / real(max(1, nptcls_active_tot)))
            call build%spproj_field%sample4update_rnd([params%fromp,params%top], update_frac_eff, nptcls_part_sel, pinds, .false. )
            l_scale_update_frac = .true.
            sig2_mul = 1.0 / (2.0 * max(update_frac_eff, TINY))
        else
            call build%spproj_field%sample4update_all([params%fromp,params%top], nptcls_part_sel, pinds, .false.)
            sig2_mul = 0.5
        endif
        nyq = build%img%get_nyq()
        allocate(sigma2(nyq,params%fromp:params%top), source=0.)
        call sum_img%new([params%box,params%box,1], params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_sum = sum_img%allocate_cmat()
        allocate(cmat_thr_sum(size(cmat_sum,dim=1), size(cmat_sum,dim=2), 1, nthr_glob))
        call compute_pspec_channel('init_pspec_part', 'sum_img_part', 'calc_pspec')
        ! destruct local objects
        call binfile%kill
        call sum_img%kill
        ! Group averaging in non-distributed execution
        if( .not.cline%defined('nparts') .and. params%part==1 )then
            call xcalc_pspec_assemble%execute(self%cline_calc_pspec_assemble)
        endif

    contains

        subroutine compute_pspec_channel(init_fbody, sum_fbody, warning_label)
            character(len=*), intent(in) :: init_fbody, sum_fbody, warning_label
            sigma2           = 0.
            ninvalid_sigma2  = 0
            call sum_img%zero_and_flag_ft
            cmat_sum = cmplx(0., 0.)
            if( nptcls_part_sel > 0 )then
                batchsz_max = min(nptcls_part_sel, 50 * nthr_glob)
                call prepimgbatch(params, build, batchsz_max)
                if( allocated(sigma2_batch) ) deallocate(sigma2_batch)
                allocate(sigma2_batch(nyq,batchsz_max), source=0.)
                ! mask and radial-shell memoization for the calc_pspec fused kernel
                call build%imgbatch(1)%memoize_powspec_coords(params%msk)
                do i = 1, nptcls_part_sel, batchsz_max
                    batchlims = [i, min(i+batchsz_max-1, nptcls_part_sel)]
                    nbatch    = batchlims(2) - batchlims(1) + 1
                    if( params%l_ptcl_src_den )then
                        call discrete_read_imgbatch_source(params, build, 'den', nbatch, pinds(batchlims(1):batchlims(2)), &
                            [1,nbatch], build%imgbatch(:nbatch))
                    else
                        call discrete_read_imgbatch(params, build, nbatch, pinds(batchlims(1):batchlims(2)), [1,nbatch])
                    endif
                    ! allocate contigous local sigma2 array for optimal caching in parallell loop
                    cmat_thr_sum = dcmplx(0.d0, 0.d0)
                    !$omp parallel do default(shared) private(iptcl,imatch,ithr,l_add_to_sum)&
                    !$omp schedule(static) proc_bind(close) reduction(+:ninvalid_sigma2)
                     do imatch = 1, nbatch
                        ithr  = omp_get_thread_num() + 1
                        iptcl = pinds(batchlims(1) + imatch - 1)
                        call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, sigma2_batch(:,imatch))
                        sigma2_batch(:,imatch) = sigma2_batch(:,imatch) * sig2_mul
                        l_add_to_sum = all(ieee_is_finite(sigma2_batch(:,imatch))) .and. any(sigma2_batch(:,imatch) > real(DTINY))
                        if( sanitize_computed_sigma2(sigma2_batch(:,imatch)) ) ninvalid_sigma2 = ninvalid_sigma2 + 1
                        ! thread average
                        if( l_add_to_sum ) call build%imgbatch(imatch)%add_dble_cmat2mat(cmat_thr_sum(:,:,:,ithr))
                    end do
                    !$omp end parallel do
                    ! global average
                    do ithr = 1, nthr_glob
                        cmat_sum(:,:,:) = cmat_sum(:,:,:) + cmplx(cmat_thr_sum(:,:,:,ithr), kind=sp)
                    end do
                    ! update non-contiguous sigma2 array to provide the correct geometry on disk
                    do imatch = 1, nbatch
                        iptcl = pinds(batchlims(1) + imatch - 1)
                        sigma2(:,iptcl) = sigma2_batch(:,imatch)
                    end do
                end do
                deallocate(sigma2_batch)
            endif
            if( ninvalid_sigma2 > 0 )then
                write(logfhandle,*) '>>> WARNING: '//trim(warning_label)//&
                    ' floored invalid computed sigma spectra to 1.0; part/fromp/top/count/selected: ', &
                    params%part, params%fromp, params%top, ninvalid_sigma2, nptcls_part_sel
            endif
            call sum_img%set_cmat(cmat_sum)
            call sum_img%write(string(sum_fbody)//int2str_pad(params%part,params%numlen)//params%ext%to_char())
            kfromto  = [1, nyq]
            binfname = trim(init_fbody)//trim(int2str(params%part))//'.dat'
            call binfile%new(binfname, params%fromp, params%top, kfromto)
            call binfile%write(sigma2)
            call binfile%kill
        end subroutine compute_pspec_channel

        logical function sanitize_computed_sigma2(sigma2_curve)
            real, intent(inout) :: sigma2_curve(:)
            logical :: invalid
            invalid = .false.
            if( .not. all(ieee_is_finite(sigma2_curve)) )then
                where( .not. ieee_is_finite(sigma2_curve) ) sigma2_curve = 1.0
                invalid = .true.
            endif
            if( .not. any(sigma2_curve > real(DTINY)) )then
                sigma2_curve = 1.0
                invalid = .true.
            endif
            sanitize_computed_sigma2 = invalid
        end function sanitize_computed_sigma2

    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, build, cline)
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_euclid :: exec_calc_pspec'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, build, cline)
        use simple_matcher_ptcl_io, only: killimgbatch
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%cline_calc_pspec_assemble%kill
        call build%kill_general_tbox
        call killimgbatch(build)
        call unmemoize_powspec_coords
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
        call cleanup_calc_pspec_outputs(params)
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

    subroutine cleanup_calc_pspec_outputs(params)
        type(parameters), intent(in) :: params
        type(string) :: fname
        integer :: ipart
        do ipart = 1,params%nparts
            fname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call del_file(fname)
            fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext%to_char()
            call del_file(fname)
            fname = SIGMA2_FBODY//int2str_pad(ipart,params%numlen)//'.dat'
            call del_file(fname)
        enddo
        call del_file('CALC_PSPEC_FINISHED')
        call del_file(CALCPSPEC_FINISHED)
        call fname%kill
    end subroutine cleanup_calc_pspec_outputs

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
