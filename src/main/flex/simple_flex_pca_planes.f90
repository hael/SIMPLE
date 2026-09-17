!@descr: flex_pca resident planes: prepped particle Fourier planes kept in memory across E-step passes
!!
!! Every E-step pass (mean, paired half-set iterations, joint fit, embedding) reads its particles
!! from disk and runs the same prep (normalise, pad, FFT, CTF and whitening into an fplane_type).
!! The prep depends only on the particle, its pose record, the sigma2 table and the run's band and
!! mask, none of which change between passes of one flex_pca run. In a shared-memory run the
!! process lives across all passes, so the planes can be prepared ONCE and served afterwards.
!!
!! The store is a plain array of fplane_type indexed by project row. A row is held when its
!! cmplx_plane is allocated. planes_batch_load is the single entry point every pass uses: a batch
!! whose rows are all held is served by copy (the E-step subtracts the mean projection in place,
!! so the caller always works on its own copy); any other batch is read and prepped exactly as
!! before and then stored, until the memory budget is reached. Distributed workers never enable
!! the store (each round is a fresh process), so their path is the read+prep branch unchanged.
!!
!! Budget: a quarter of MemAvailable from /proc/meminfo at enable time, SIMPLE_COV_RESIDENT_GB overrides
!! it, SIMPLE_COV_RESIDENT=0 disables the store (A/B switch only).
module simple_flex_pca_planes
use simple_core_module_api
use simple_builder,                       only: builder
use simple_parameters,                    only: parameters
use simple_matcher_ptcl_io,               only: discrete_read_imgbatch
use simple_flex_pca_plane_cache,          only: plane_cache_read_batch
use simple_flex_reconstructor_latent_ops, only: prep_imgs4projected_model
implicit none
private
#include "simple_local_flags.inc"

public :: planes_enable, planes_kill, planes_batch_load, planes_enabled

type(fplane_type), allocatable :: held(:)        !< by project row; held when cmplx_plane is allocated
logical  :: l_on         = .false.
logical  :: l_full       = .false.               !< budget reached: no more rows are stored
logical  :: l_key_set    = .false.
logical  :: l_cached_key = .false.               !< read path the held planes came from
integer  :: nheld        = 0
integer  :: nserved      = 0                     !< batches served from the store
integer  :: nprepped     = 0                     !< batches read and prepped
real(dp) :: gb_held      = 0.d0
real(dp) :: gb_budget    = 0.d0

contains

    !> Turn the store on for a project of nrows particle rows (shared-memory runs only). Only the
    !! cropped grid is worth holding: a plane prepared from the full-box particle lives on the
    !! full padded lattice (5 MB at box 360) while the cache serves it on the box_crop lattice
    !! (~40x smaller), so the store is tied to the particle cache being in use.
    subroutine planes_enable( nrows, cropped )
        integer, intent(in) :: nrows
        logical, intent(in) :: cropped
        character(len=32) :: envval
        integer  :: ln, stat, ival, ios
        real(dp) :: gb_avail
        if( l_on ) return
        if( .not. cropped )then
            write(logfhandle,'(A)') '>>> FLEX_PCA RESIDENT PLANES OFF: planes are held only on the cropped grid &
                &(run with cache=yes to enable)'
            return
        endif
        gb_avail  = mem_available_gb()
        gb_budget = 0.25d0 * gb_avail
        if( gb_budget <= 0.d0 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA RESIDENT PLANES OFF (no memory budget: MemAvailable &
                &unreadable and SIMPLE_COV_RESIDENT_GB unset)'
            return
        endif
        allocate(held(nrows))
        l_on = .true.
        write(logfhandle,'(A,I0,A,F8.1,A,F8.1,A)') '>>> FLEX_PCA RESIDENT PLANES ON: ', nrows, &
            &' rows addressable, budget ', gb_budget, ' GB (MemAvailable ', gb_avail, ' GB)'
        call flush(logfhandle)
    end subroutine planes_enable

    logical function planes_enabled()
        planes_enabled = l_on
    end function planes_enabled

    !> Serve batchlims of the pinds list from the store when every row is held; otherwise read
    !! (plane cache or stacks) and prep exactly as the passes always did, then store what fits.
    subroutine planes_batch_load( params, build, n, pinds, batchlims, fpls, mskrad, cached, sec_read, sec_prep )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(fplane_type), intent(inout) :: fpls(:)
        real,              intent(in)    :: mskrad
        logical,           intent(in)    :: cached
        real(timer_int_kind), optional, intent(inout) :: sec_read, sec_prep
        integer(timer_int_kind) :: t
        integer :: i, batchsz
        batchsz = batchlims(2) - batchlims(1) + 1
        if( l_on )then
            if( .not. l_key_set )then
                l_cached_key = cached
                l_key_set    = .true.
            elseif( cached .neqv. l_cached_key )then
                THROW_HARD('resident planes were prepared from a different read path than requested')
            endif
            if( all_held(batchsz, pinds(batchlims(1):batchlims(2))) )then
                t = tic()
                !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
                do i = 1, batchsz
                    fpls(i) = held(pinds(batchlims(1)+i-1))
                end do
                !$omp end parallel do
                if( present(sec_read) ) sec_read = sec_read + toc(t)
                nserved = nserved + 1
                return
            endif
        endif
        t = tic()
        if( cached )then
            call plane_cache_read_batch(params, n, pinds, batchlims)
        else
            call discrete_read_imgbatch(params, build, n, pinds, batchlims)
        endif
        if( present(sec_read) ) sec_read = sec_read + toc(t)
        t = tic()
        call prep_imgs4projected_model(params, build, batchsz, build%imgbatch(:batchsz), &
            &pinds(batchlims(1):batchlims(2)), fpls(:batchsz), mskrad=mskrad, cached=cached)
        if( present(sec_prep) ) sec_prep = sec_prep + toc(t)
        nprepped = nprepped + 1
        if( l_on ) call store(batchsz, pinds(batchlims(1):batchlims(2)), fpls)
    end subroutine planes_batch_load

    logical function all_held( batchsz, rows )
        integer, intent(in) :: batchsz, rows(batchsz)
        integer :: i
        all_held = .false.
        do i = 1, batchsz
            if( .not. allocated(held(rows(i))%cmplx_plane) ) return
        end do
        all_held = .true.
    end function all_held

    !> Copy the freshly prepped planes of a batch into the store, stopping at the budget.
    subroutine store( batchsz, rows, fpls )
        integer,           intent(in) :: batchsz, rows(batchsz)
        type(fplane_type), intent(in) :: fpls(:)
        real(dp) :: gb
        integer  :: i
        if( l_full ) return
        do i = 1, batchsz
            if( allocated(held(rows(i))%cmplx_plane) ) cycle
            gb = plane_gb(fpls(i))
            if( gb_held + gb > gb_budget )then
                l_full = .true.
                write(logfhandle,'(A,I0,A,F8.1,A)') '>>> FLEX_PCA RESIDENT PLANES budget reached: ', &
                    &nheld, ' rows held, ', gb_held, ' GB; the remaining rows are read every pass'
                call flush(logfhandle)
                return
            endif
            held(rows(i)) = fpls(i)
            nheld   = nheld + 1
            gb_held = gb_held + gb
        end do
    end subroutine store

    real(dp) function plane_gb( fpl )
        type(fplane_type), intent(in) :: fpl
        real(dp) :: bytes
        bytes = 0.d0
        if( allocated(fpl%cmplx_plane) )    bytes = bytes + 8.d0 * real(size(fpl%cmplx_plane),dp)
        if( allocated(fpl%ctfsq_plane) )    bytes = bytes + 4.d0 * real(size(fpl%ctfsq_plane),dp)
        if( allocated(fpl%transfer_plane) ) bytes = bytes + 8.d0 * real(size(fpl%transfer_plane),dp)
        plane_gb = bytes / 1.d9
    end function plane_gb

    !> Release the store and report what it did.
    subroutine planes_kill()
        integer :: i
        if( .not. l_on ) return
        write(logfhandle,'(A,I0,A,F8.1,A,I0,A,I0,A)') '>>> FLEX_PCA RESIDENT PLANES: held ', nheld, &
            &' rows (', gb_held, ' GB), batches served=', nserved, ' read+prepped=', nprepped
        call flush(logfhandle)
        do i = 1, size(held)
            if( allocated(held(i)%cmplx_plane) )    deallocate(held(i)%cmplx_plane)
            if( allocated(held(i)%ctfsq_plane) )    deallocate(held(i)%ctfsq_plane)
            if( allocated(held(i)%transfer_plane) ) deallocate(held(i)%transfer_plane)
        end do
        deallocate(held)
        l_on = .false.; l_full = .false.; l_key_set = .false.
        nheld = 0; nserved = 0; nprepped = 0; gb_held = 0.d0; gb_budget = 0.d0
    end subroutine planes_kill

    !> MemAvailable from /proc/meminfo in GB; zero when unreadable (non-Linux).
    real(dp) function mem_available_gb()
        character(len=256) :: line
        integer  :: funit, ios, kb
        mem_available_gb = 0.d0
        open(newunit=funit, file='/proc/meminfo', status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(funit, '(A)', iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, 'MemAvailable:') == 1 )then
                read(line(14:), *, iostat=ios) kb
                if( ios == 0 ) mem_available_gb = real(kb,dp) / 1.d6
                exit
            endif
        end do
        close(funit)
    end function mem_available_gb

end module simple_flex_pca_planes
