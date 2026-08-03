module simple_calc_pspec_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,        only: builder
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_qsys_funs,      only: qsys_job_finished
use simple_image,          only: image, unmemoize_powspec_coords
use simple_sigma2_binfile, only: sigma2_binfile
use simple_ran_tabu,       only: ran_tabu
implicit none

public :: calc_pspec_strategy, calc_pspec_inmem_strategy, calc_pspec_partitioned_strategy, create_calc_pspec_strategy
private
#include "simple_local_flags.inc"

integer,          parameter :: SIGMA2_BOOTSTRAP_MAX_PARTICLES = 25000
character(len=*), parameter :: SIGMA2_BOOTSTRAP_SELECTION_FNAME = 'sigma2_bootstrap_selection.dat'

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

! Computes every partition of the sigma2 bootstrap in a single process. The
! artifacts written are the ones nparts distributed workers would have produced,
! so downstream consumers of init_pspec_part*/sum_img_part*/sigma2_noise_part*
! see an unchanged on-disk contract.
type, extends(calc_pspec_strategy) :: calc_pspec_partitioned_strategy
    type(cmdline)        :: cline_calc_pspec_assemble
    integer, allocatable :: parts(:,:)
contains
    procedure :: initialize   => partitioned_initialize
    procedure :: execute      => partitioned_execute
    procedure :: finalize_run => partitioned_finalize_run
    procedure :: cleanup      => partitioned_cleanup
end type calc_pspec_partitioned_strategy

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
            allocate(calc_pspec_partitioned_strategy :: strategy)
            strategy%end_msg          = '**** SIMPLE_CALC_PSPEC NORMAL STOP ****'
            strategy%end_print_simple = .true.
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> SHARED-MEMORY PARTITIONED EXECUTION'
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
        if( cline%defined('part') )then
            if( .not.file_exists(SIGMA2_BOOTSTRAP_SELECTION_FNAME) ) &
                &THROW_HARD('Missing global sigma2 bootstrap selection for calc_pspec worker')
        else
            call write_sigma2_bootstrap_selection(params, build)
        endif
        self%cline_calc_pspec_assemble = cline
        call self%cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble')
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, build, cline)
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        integer, allocatable :: sel_pinds(:)
        real :: sig2_mul(2)
        ! Sigma2 is bootstrapped once from a capped, random global E/O sample.
        ! The resulting two curves are propagated to every particle record;
        ! matcher residual updates establish any requested group structure later.
        call read_sigma2_bootstrap_selection(params, build, sel_pinds, sig2_mul)
        call compute_pspec_partitions(params, build, [params%part],&
            &reshape([params%fromp, params%top], [1,2]), sel_pinds, sig2_mul)
        deallocate(sel_pinds)
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
        use simple_matcher_ptcl_io, only: killimgbatch
        class(calc_pspec_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        type(builder),                    intent(inout) :: build
        class(cmdline),                   intent(inout) :: cline
        call self%cline_calc_pspec_assemble%kill
        call build%kill_general_tbox
        call killimgbatch(build)
        call unmemoize_powspec_coords
        if( .not.cline%defined('part') ) call del_file(SIGMA2_BOOTSTRAP_SELECTION_FNAME)
    end subroutine inmem_cleanup

    ! =====================================================================
    ! SHARED-MEMORY PARTITIONED IMPLEMENTATION
    ! =====================================================================

    subroutine partitioned_initialize(self, params, build, cline)
        class(calc_pspec_partitioned_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        call cline%set('stream', 'no')
        if( .not. cline%defined('projfile') )then
            THROW_HARD('Missing project file entry; exec_calc_pspec')
        endif
        ! Parse parameters (mkdir default kept as provided/assigned by commander-level defaults)
        call params%new(cline)
        call cleanup_calc_pspec_outputs(params)
        ! Project + 2D toolbox, this process does the particle work itself
        call build%build_spproj(params, cline, wthreads=.true.)
        call build%build_general_tbox(params, cline, do3d=.false.)
        ! Sanity check + ensure EO partition exists and is persisted
        call sanity_check_calc_pspec_input(params, build)
        call ensure_calc_pspec_eo_partition(params, build)
        call write_sigma2_bootstrap_selection(params, build)
        ! Set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! Partition ranges must be the ones the downstream distributed consumers
        ! derive, i.e. qsys_env's split_mode='even' partitioning of nptcls.
        self%parts = split_nobjs_even(params%nptcls, params%nparts)
        self%cline_calc_pspec_assemble = cline
        call self%cline_calc_pspec_assemble%set('prg', 'calc_pspec_assemble')
    end subroutine partitioned_initialize

    subroutine partitioned_execute(self, params, build, cline)
        use simple_matcher_ptcl_io,         only: killimgbatch
        use simple_commanders_euclid_distr, only: commander_calc_pspec_assemble
        class(calc_pspec_partitioned_strategy), intent(inout) :: self
        type(parameters),                       intent(inout) :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        type(commander_calc_pspec_assemble) :: xcalc_pspec_assemble
        integer, allocatable :: sel_pinds(:), part_ids(:)
        real    :: sig2_mul(2)
        integer :: ipart
        call read_sigma2_bootstrap_selection(params, build, sel_pinds, sig2_mul)
        part_ids = [(ipart, ipart=1,params%nparts)]
        call compute_pspec_partitions(params, build, part_ids, self%parts, sel_pinds, sig2_mul)
        deallocate(sel_pinds, part_ids)
        ! Release the compute toolbox before the assemble commander builds its own
        call killimgbatch(build)
        call unmemoize_powspec_coords
        call build%kill_general_tbox
        call xcalc_pspec_assemble%execute(self%cline_calc_pspec_assemble)
    end subroutine partitioned_execute

    subroutine partitioned_finalize_run(self, params, build, cline)
        class(calc_pspec_partitioned_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        ! No-op
    end subroutine partitioned_finalize_run

    subroutine partitioned_cleanup(self, params, build, cline)
        class(calc_pspec_partitioned_strategy), intent(inout) :: self
        type(parameters),                       intent(in)    :: params
        type(builder),                          intent(inout) :: build
        class(cmdline),                         intent(inout) :: cline
        call self%cline_calc_pspec_assemble%kill
        if( allocated(self%parts) ) deallocate(self%parts)
        ! kill_general_tbox is idempotent; execute already released it
        call build%kill_general_tbox
        call del_file(SIGMA2_BOOTSTRAP_SELECTION_FNAME)
        call simple_touch(CALCPSPEC_FINISHED)
    end subroutine partitioned_cleanup

    !>  Computes the bootstrap sigma2 spectra and Fourier-space sums for a set of
    !!  partitions and writes the two artifacts per partition that
    !!  calc_pspec_assemble consumes. part_ids supplies the file index of each
    !!  partition, part_ranges its [fromp,top].
    !!
    !!  Particle images are read once, in batches spanning the whole selection,
    !!  independently of where the partition boundaries fall. Accumulation is
    !!  partition-local: a batch is split into the contiguous runs it shares with
    !!  each partition (sel_pinds is sorted and the ranges are disjoint), and each
    !!  run reduces into that partition's sum, so no two partitions ever contend.
    !!  Particles outside the bootstrap sample keep a zero spectrum; the assemble
    !!  step compensates via sig2_mul.
    subroutine compute_pspec_partitions( params, build, part_ids, part_ranges, sel_pinds, sig2_mul )
        use simple_matcher_ptcl_io, only: prepimgbatch, discrete_read_imgbatch, discrete_read_imgbatch_source
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        integer,          intent(in)    :: part_ids(:)
        integer,          intent(in)    :: part_ranges(:,:)
        integer,          intent(in)    :: sel_pinds(:)
        real,             intent(in)    :: sig2_mul(2)
        type(image)          :: sum_img
        type(sigma2_binfile) :: binfile
        type(string)         :: binfname
        complex(dp), allocatable :: cmat_thr_sum(:,:,:,:)
        complex,     allocatable :: cmat_part_sum(:,:,:,:), cmat_dims(:,:,:)
        integer,     allocatable :: work_pinds(:), work_slot(:), slot_from(:), slot_to(:)
        real,        allocatable :: sigma2(:,:), sigma2_sel(:,:), sigma2_batch(:,:)
        integer :: batchlims(2)
        integer :: i, nyq, nselected, n_work, batchsz_max, nbatch
        integer :: islot, nslots, ipart, fromp, top, run_lo, run_hi
        integer :: ninvalid_sigma2
        nslots          = size(part_ids)
        nselected       = size(sel_pinds)
        ninvalid_sigma2 = 0
        nyq             = build%img%get_nyq()
        ! Map the shared selection onto the requested partitions, dropping any
        ! particle none of them owns (the single-partition worker case).
        allocate(work_pinds(nselected), work_slot(nselected), source=0)
        n_work = 0
        do i = 1,nselected
            islot = slot_of(sel_pinds(i))
            if( islot == 0 ) cycle
            n_work             = n_work + 1
            work_pinds(n_work) = sel_pinds(i)
            work_slot(n_work)  = islot
        end do
        allocate(slot_from(nslots), source=0)
        allocate(slot_to(nslots),   source=-1)
        do i = 1,n_work
            islot = work_slot(i)
            if( slot_from(islot) == 0 ) slot_from(islot) = i
            slot_to(islot) = i
        end do
        call sum_img%new([params%box,params%box,1], params%smpd)
        call sum_img%zero_and_flag_ft
        cmat_dims = sum_img%allocate_cmat()
        allocate(cmat_part_sum(size(cmat_dims,dim=1), size(cmat_dims,dim=2), 1, nslots), source=cmplx(0.,0.))
        allocate(cmat_thr_sum( size(cmat_dims,dim=1), size(cmat_dims,dim=2), 1, nthr_glob))
        allocate(sigma2_sel(nyq,max(1,n_work)), source=0.)
        if( n_work > 0 )then
            batchsz_max = min(n_work, 50 * nthr_glob)
            call prepimgbatch(params, build, batchsz_max)
            allocate(sigma2_batch(nyq,batchsz_max), source=0.)
            ! mask and radial-shell memoization for the calc_pspec fused kernel
            call build%imgbatch(1)%memoize_powspec_coords(params%msk)
            do i = 1, n_work, batchsz_max
                batchlims = [i, min(i+batchsz_max-1, n_work)]
                nbatch    = batchlims(2) - batchlims(1) + 1
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', nbatch, work_pinds(batchlims(1):batchlims(2)), &
                        [1,nbatch], build%imgbatch(:nbatch))
                else
                    call discrete_read_imgbatch(params, build, nbatch, work_pinds(batchlims(1):batchlims(2)), [1,nbatch])
                endif
                ! reduce the batch one partition-run at a time
                do islot = 1,nslots
                    if( slot_from(islot) == 0 ) cycle
                    run_lo = max(slot_from(islot), batchlims(1))
                    run_hi = min(slot_to(islot),   batchlims(2))
                    if( run_hi < run_lo ) cycle
                    call accumulate_run(islot, run_lo, run_hi, batchlims(1))
                end do
            end do
            deallocate(sigma2_batch)
        endif
        if( ninvalid_sigma2 > 0 )then
            write(logfhandle,*) '>>> WARNING: calc_pspec'//&
                ' floored invalid computed sigma spectra to 1.0; count/selected: ', &
                ninvalid_sigma2, n_work
        endif
        ! scatter the compact spectra into partition-local geometry and write
        do islot = 1,nslots
            ipart = part_ids(islot)
            fromp = part_ranges(islot,1)
            top   = part_ranges(islot,2)
            allocate(sigma2(nyq,fromp:top), source=0.)
            if( slot_from(islot) > 0 )then
                do i = slot_from(islot), slot_to(islot)
                    sigma2(:,work_pinds(i)) = sigma2_sel(:,i)
                end do
            endif
            call sum_img%zero_and_flag_ft
            call sum_img%set_cmat(cmat_part_sum(:,:,:,islot))
            call sum_img%write(string('sum_img_part')//int2str_pad(ipart,params%numlen)//params%ext%to_char())
            binfname = 'init_pspec_part'//trim(int2str(ipart))//'.dat'
            call binfile%new(binfname, fromp, top, [1,nyq])
            call binfile%write(sigma2)
            call binfile%kill
            deallocate(sigma2)
        end do
        call sum_img%kill
        call binfname%kill
        deallocate(work_pinds, work_slot, slot_from, slot_to)
        deallocate(sigma2_sel, cmat_part_sum, cmat_thr_sum, cmat_dims)

    contains

        !>  Index of the partition owning iptcl, 0 when none does
        integer function slot_of( iptcl_here )
            integer, intent(in) :: iptcl_here
            integer :: jslot
            slot_of = 0
            do jslot = 1,nslots
                if( iptcl_here >= part_ranges(jslot,1) .and. iptcl_here <= part_ranges(jslot,2) )then
                    slot_of = jslot
                    return
                endif
            end do
        end function slot_of

        !>  Reduces work items [run_lo,run_hi], all owned by partition islot, into
        !!  that partition's Fourier sum. batch_base maps work index to batch slot.
        subroutine accumulate_run( islot_here, run_lo_here, run_hi_here, batch_base )
            integer, intent(in) :: islot_here, run_lo_here, run_hi_here, batch_base
            integer :: j, iwork, imatch, nrun, eo_here, ithr_here, iptcl_here
            logical :: l_add_to_sum
            nrun         = run_hi_here - run_lo_here + 1
            cmat_thr_sum = dcmplx(0.d0, 0.d0)
            !$omp parallel do default(shared) private(j,iwork,imatch,iptcl_here,ithr_here,eo_here,l_add_to_sum)&
            !$omp schedule(static) proc_bind(close) reduction(+:ninvalid_sigma2)
            do j = 1, nrun
                iwork      = run_lo_here + j - 1
                imatch     = iwork - batch_base + 1
                ithr_here  = omp_get_thread_num() + 1
                iptcl_here = work_pinds(iwork)
                call build%imgbatch(imatch)%norm_noise_mask_fft_powspec(build%lmsk, params%msk, sigma2_batch(:,imatch))
                eo_here = build%spproj_field%get_eo(iptcl_here)
                sigma2_batch(:,imatch) = sigma2_batch(:,imatch) * sig2_mul(eo_here+1)
                l_add_to_sum = all(ieee_is_finite(sigma2_batch(:,imatch))) .and. any(sigma2_batch(:,imatch) > real(DTINY))
                if( sanitize_computed_sigma2(sigma2_batch(:,imatch)) ) ninvalid_sigma2 = ninvalid_sigma2 + 1
                ! thread average
                if( l_add_to_sum ) call build%imgbatch(imatch)%add_dble_cmat2mat(cmat_thr_sum(:,:,:,ithr_here))
            end do
            !$omp end parallel do
            ! partition-local average
            do ithr_here = 1, nthr_glob
                cmat_part_sum(:,:,:,islot_here) = cmat_part_sum(:,:,:,islot_here)&
                    &+ cmplx(cmat_thr_sum(:,:,:,ithr_here), kind=sp)
            end do
            ! keep the compact spectra, scattered to disk geometry after the read pass
            do j = 1, nrun
                iwork  = run_lo_here + j - 1
                imatch = iwork - batch_base + 1
                sigma2_sel(:,iwork) = sigma2_batch(:,imatch)
            end do
        end subroutine accumulate_run

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

    end subroutine compute_pspec_partitions

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
        call del_file(SIGMA2_BOOTSTRAP_SELECTION_FNAME)
        call fname%kill
    end subroutine cleanup_calc_pspec_outputs

    subroutine write_sigma2_bootstrap_selection(params, build)
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(ran_tabu) :: rt
        integer, allocatable :: candidates(:,:), pinds(:)
        integer :: nactive(2), nsample(2), nfilled(2), nselected, iptcl, eo, i, ios, fd, offset
        logical :: is_open
        nactive = 0
        do iptcl = 1,params%nptcls
            if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
            eo = build%spproj_field%get_eo(iptcl)
            if( eo < 0 .or. eo > 1 ) THROW_HARD('Invalid even/odd label in sigma2 bootstrap selection')
            nactive(eo+1) = nactive(eo+1) + 1
        enddo
        if( sum(nactive) < 1 ) THROW_HARD('No active particles for sigma2 bootstrap selection')
        nsample(1) = min(nactive(1), SIGMA2_BOOTSTRAP_MAX_PARTICLES / 2)
        nsample(2) = min(nactive(2), SIGMA2_BOOTSTRAP_MAX_PARTICLES - nsample(1))
        if( sum(nsample) < SIGMA2_BOOTSTRAP_MAX_PARTICLES )then
            nsample(1) = nsample(1) + min(SIGMA2_BOOTSTRAP_MAX_PARTICLES - sum(nsample), nactive(1) - nsample(1))
            nsample(2) = nsample(2) + min(SIGMA2_BOOTSTRAP_MAX_PARTICLES - sum(nsample), nactive(2) - nsample(2))
        endif
        nselected = sum(nsample)
        allocate(candidates(maxval(nactive),2), source=0)
        nfilled = 0
        do iptcl = 1,params%nptcls
            if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
            eo = build%spproj_field%get_eo(iptcl) + 1
            nfilled(eo) = nfilled(eo) + 1
            candidates(nfilled(eo),eo) = iptcl
        enddo
        allocate(pinds(nselected), source=0)
        offset = 0
        do eo = 1,2
            if( nsample(eo) == 0 ) cycle
            if( nsample(eo) < nactive(eo) )then
                rt = ran_tabu(nactive(eo))
                call rt%shuffle(candidates(:nactive(eo),eo))
                call rt%kill
            endif
            pinds(offset+1:offset+nsample(eo)) = candidates(:nsample(eo),eo)
            offset = offset + nsample(eo)
        enddo
        call hpsort(pinds)
        is_open = .false.
        open(newunit=fd, file=SIGMA2_BOOTSTRAP_SELECTION_FNAME, status='replace', action='write', iostat=ios)
        if( ios /= 0 ) THROW_HARD('Cannot open sigma2 bootstrap selection for writing')
        is_open = .true.
        write(fd,'(4(I0,1X))',iostat=ios) nactive(1), nactive(2), nsample(1), nsample(2)
        if( ios /= 0 ) THROW_HARD('Cannot write sigma2 bootstrap selection header')
        do i = 1,nselected
            iptcl = pinds(i)
            write(fd,'(I0,1X,I0)',iostat=ios) iptcl, build%spproj_field%get_eo(iptcl)
            if( ios /= 0 ) THROW_HARD('Cannot write sigma2 bootstrap selection record')
        enddo
        if( is_open ) close(fd)
        write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8)') &
            '>>> SIGMA2 GLOBAL BOOTSTRAP ACTIVE E/O: ', nactive(1), '/', nactive(2), &
            '; SAMPLED E/O: ', nsample(1), '/', nsample(2), '; CAP: ', SIGMA2_BOOTSTRAP_MAX_PARTICLES
        deallocate(candidates, pinds)
    end subroutine write_sigma2_bootstrap_selection

    !>  Reads the global bootstrap selection shared by every partition. The
    !!  returned indices are sorted and span the whole particle range; callers
    !!  slice them to their own [fromp,top].
    subroutine read_sigma2_bootstrap_selection(params, build, sel_pinds, sig2_mul)
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        integer, allocatable, intent(out) :: sel_pinds(:)
        real,             intent(out)   :: sig2_mul(2)
        integer :: nactive(2), nsample(2), nselected, eo, i, ios, fd
        logical :: is_open
        if( .not.file_exists(SIGMA2_BOOTSTRAP_SELECTION_FNAME) ) &
            &THROW_HARD('Missing sigma2 bootstrap selection for calc_pspec worker')
        is_open = .false.
        open(newunit=fd, file=SIGMA2_BOOTSTRAP_SELECTION_FNAME, status='old', action='read', iostat=ios)
        if( ios /= 0 ) THROW_HARD('Cannot open sigma2 bootstrap selection for reading')
        is_open = .true.
        read(fd,*,iostat=ios) nactive(1), nactive(2), nsample(1), nsample(2)
        if( ios /= 0 ) THROW_HARD('Cannot read sigma2 bootstrap selection header')
        nselected = sum(nsample)
        if( nselected < 1 .or. nselected > SIGMA2_BOOTSTRAP_MAX_PARTICLES ) &
            &THROW_HARD('Invalid sigma2 bootstrap selection size')
        if( any(nsample < 0) .or. any(nsample > nactive) ) &
            &THROW_HARD('Invalid sigma2 bootstrap selection counts')
        allocate(sel_pinds(nselected), source=0)
        do i = 1,nselected
            read(fd,*,iostat=ios) sel_pinds(i), eo
            if( ios /= 0 ) THROW_HARD('Cannot read sigma2 bootstrap selection record')
            if( sel_pinds(i) < 1 .or. sel_pinds(i) > params%nptcls .or. eo < 0 .or. eo > 1 ) &
                &THROW_HARD('Invalid sigma2 bootstrap selection record')
            if( build%spproj_field%get_eo(sel_pinds(i)) /= eo ) &
                &THROW_HARD('Sigma2 bootstrap selection even/odd mismatch')
        enddo
        if( is_open ) close(fd)
        sig2_mul = 0.
        do eo = 1,2
            if( nsample(eo) > 0 ) sig2_mul(eo) = real(nactive(eo)) / (2. * real(nsample(eo)))
        enddo
    end subroutine read_sigma2_bootstrap_selection

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
