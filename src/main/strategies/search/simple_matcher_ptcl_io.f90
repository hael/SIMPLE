!@descr: particle image batch I/O routines shared by matcher workflows
module simple_matcher_ptcl_io
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_stack_read_session, only: stack_read_session, stack_read_counts
use simple_syslib,            only: get_peak_rss_bytes, get_node_mem_in_use_bytes, get_node_mem_total_bytes
use simple_image_slab,        only: image_slab, SLAB_OK, SLAB_OVER_BUDGET, SLAB_NO_MEMORY
use simple_image,             only: image
use, intrinsic :: iso_fortran_env, only: int64
implicit none
#include "simple_local_flags.inc"

public :: prepimgbatch, killimgbatch, read_imgbatch, ptcl_read_session
public :: prep_rec_observation
public :: ptcl_read_counts, ptcl_read_counts_now, ptcl_read_counts_since, add_ptcl_cache_read_time
public :: write_ptcl_read_bench_context, write_mem_bench_context
public :: ptcl_buffer_budget, fill_pass_buffer, release_pass_buffer, pass_buffer_status_name, pass_buffer_bytes
public :: BUFFER_FILLED, BUFFER_OFF, BUFFER_OVER_BUDGET, BUFFER_NO_MEMORY
private

! outcomes of filling the particle buffer
integer, parameter :: BUFFER_FILLED = 0, BUFFER_OFF = 1, BUFFER_OVER_BUDGET = 2, BUFFER_NO_MEMORY = 3

interface read_imgbatch
    module procedure read_imgbatch_1
    module procedure read_imgbatch_2
    module procedure read_imgbatch_3
end interface read_imgbatch

type(stack_io) :: stkio_r

!> The particle reads of every read session, cumulative over the process and updated by the
!! master thread; a pass's counts are the difference of two snapshots (ptcl_read_counts_since).
!! The downscaled 2D cache has its own reader: only its read time is kept.
type :: ptcl_read_counts
    integer(int64)       :: nopens     = 0_int64 !< stacks opened
    integer(int64)       :: nrequests  = 0_int64 !< read requests
    integer(int64)       :: bytes_req  = 0_int64 !< file bytes of the images delivered
    integer(int64)       :: bytes_read = 0_int64 !< file bytes read, the gaps spans bridge and the window included
    integer(int64)       :: bytes_prefetched    = 0_int64 !< file bytes read ahead into the prefetch window
    integer(int64)       :: bytes_prefetch_used = 0_int64 !< window bytes images were decoded from
    integer(int64)       :: bytes_advised       = 0_int64 !< file bytes hinted to the kernel beyond the window
    integer(int64)       :: nbuffered  = 0_int64 !< images copied from the particle buffer, no file touched
    real(timer_int_kind) :: rt_read    = 0.      !< wall time reading stacks (L_BENCH_GLOB)
    real(timer_int_kind) :: rt_cache   = 0.      !< wall time reading the downscaled 2D cache (L_BENCH_GLOB)
end type ptcl_read_counts

type(ptcl_read_counts) :: read_totals

! The particle buffer of one 3D search pass: its share's images, read once into an immutable slab and copied
! from it by the search batches and the reconstruction of the same process (section 3.1 of the plan).
type(image_slab)     :: pass_slab
integer, allocatable :: pass_slot(:)    !< slab slot of particle p at pass_slot(p - pass_lo + 1), 0 when absent
integer              :: pass_lo  = 0
integer              :: pass_box = 0
logical              :: l_pass_ready = .false.

!> The particle images of one pass (a search, a reconstruction, a probability table), read through a
!! stack read session. One per pass: new before its first batch, kill after its last. After each batch
!! the session reads the next one ahead: the batch after batchlims when read is given the pass's whole
!! index list, or the next batch of the list given to expect.
type :: ptcl_read_session
    private
    type(stack_read_session) :: stacks
    integer, allocatable     :: pass(:)   !< the pass's particles in read order (expect)
    integer                  :: cursor = 1
contains
    procedure :: new    => new_ptcl_read_session
    procedure :: expect => expect_ptcl_pass
    procedure :: read   => read_ptcl_batch
    procedure :: kill   => kill_ptcl_read_session
end type ptcl_read_session

contains

    type(ptcl_read_counts) function ptcl_read_counts_now()
        ptcl_read_counts_now = read_totals
    end function ptcl_read_counts_now

    !> The reads since the snapshot t0
    type(ptcl_read_counts) function ptcl_read_counts_since( t0 ) result( d )
        type(ptcl_read_counts), intent(in) :: t0
        d%nopens     = read_totals%nopens     - t0%nopens
        d%nrequests  = read_totals%nrequests  - t0%nrequests
        d%bytes_req  = read_totals%bytes_req  - t0%bytes_req
        d%bytes_read = read_totals%bytes_read - t0%bytes_read
        d%bytes_prefetched    = read_totals%bytes_prefetched    - t0%bytes_prefetched
        d%bytes_prefetch_used = read_totals%bytes_prefetch_used - t0%bytes_prefetch_used
        d%bytes_advised       = read_totals%bytes_advised       - t0%bytes_advised
        d%nbuffered           = read_totals%nbuffered           - t0%nbuffered
        d%rt_read    = read_totals%rt_read    - t0%rt_read
        d%rt_cache   = read_totals%rt_cache   - t0%rt_cache
    end function ptcl_read_counts_since

    !> The 2D cache has its own reader (simple_ptcl_cache); its callers time it
    subroutine add_ptcl_cache_read_time( rt )
        real(timer_int_kind), intent(in) :: rt
        read_totals%rt_cache = read_totals%rt_cache + rt
    end subroutine add_ptcl_cache_read_time

    !> Context lines of a pass's particle reads for a bench file
    subroutine write_ptcl_read_bench_context( fnr, prefix, counts )
        integer,                intent(in) :: fnr
        character(len=*),       intent(in) :: prefix
        type(ptcl_read_counts), intent(in) :: counts
        write(fnr,'(a,i0)') prefix//' particle stacks opened      : ', counts%nopens
        write(fnr,'(a,i0)') prefix//' particle read requests      : ', counts%nrequests
        write(fnr,'(a,i0)') prefix//' particle bytes requested    : ', counts%bytes_req
        write(fnr,'(a,i0)') prefix//' particle bytes read         : ', counts%bytes_read
        write(fnr,'(a,i0)') prefix//' particle bytes prefetched   : ', counts%bytes_prefetched
        write(fnr,'(a,i0)') prefix//' particle bytes prefetch used: ', counts%bytes_prefetch_used
        write(fnr,'(a,i0)') prefix//' particle bytes advised      : ', counts%bytes_advised
        write(fnr,'(a,i0)') prefix//' particle images from buffer : ', counts%nbuffered
    end subroutine write_ptcl_read_bench_context

    !> Context lines of the process's peak RSS and the node's memory in use, sampled together
    subroutine write_mem_bench_context( fnr, prefix )
        integer,          intent(in) :: fnr
        character(len=*), intent(in) :: prefix
        write(fnr,'(a,i0)') prefix//' peak RSS (bytes)            : ', get_peak_rss_bytes()
        write(fnr,'(a,i0)') prefix//' node memory in use (bytes)  : ', get_node_mem_in_use_bytes()
    end subroutine write_mem_bench_context

    subroutine prepimgbatch( params, build, batchsz, box, smpd )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: batchsz
        integer, optional, intent(in)    :: box
        real,    optional, intent(in)    :: smpd
        integer :: currsz, ibatch, box_here, ldim_curr(3)
        real    :: smpd_here
        logical :: doprep
        box_here  = params%box
        smpd_here = params%smpd
        if( present(box)  ) box_here  = box
        if( present(smpd) ) smpd_here = smpd
        if( .not. allocated(build%imgbatch) )then
            doprep = .true.
        else
            currsz    = size(build%imgbatch)
            ldim_curr = build%imgbatch(1)%get_ldim()
            ! a batch read from the downscaled cache needs box_crop buffers at
            ! smpd_crop, so both have to be honoured on reuse, not just the batch size
            if( batchsz > currsz .or. ldim_curr(1) /= box_here .or. &
               &abs(build%imgbatch(1)%get_smpd() - smpd_here) > 1.e-6 )then
                call killimgbatch(build)
                doprep = .true.
            else
                doprep = .false.
            endif
        endif
        if( doprep )then
            allocate(build%imgbatch(batchsz))
            !$omp parallel do default(shared) private(ibatch) schedule(static) proc_bind(close)
            do ibatch = 1,batchsz
                call build%imgbatch(ibatch)%new([box_here, box_here, 1], smpd_here, wthreads=.false.)
            end do
            !$omp end parallel do
        endif
    end subroutine prepimgbatch

    subroutine killimgbatch( build )
        type(builder), intent(inout) :: build
        integer :: ibatch
        if( allocated(build%imgbatch) )then
            do ibatch=1,size(build%imgbatch)
                call build%imgbatch(ibatch)%kill
            end do
            deallocate(build%imgbatch)
        endif
    end subroutine killimgbatch

    subroutine read_imgbatch_1( params, build, fromptop )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: fromptop(2)
        type(string) :: stkname
        integer :: iptcl, ind_in_batch, ind_in_stk
        do iptcl=fromptop(1),fromptop(2)
            ind_in_batch = iptcl - fromptop(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ind_in_batch))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_1

    subroutine read_imgbatch_2( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(string) :: stkname
        integer :: ind_in_stk, i, ii
        do i=batchlims(1),batchlims(2)
            ii = i - batchlims(1) + 1
            call build%spproj%get_stkname_and_ind(params%oritype, pinds(i), stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%imgbatch(ii))
        end do
        call stkio_r%close
    end subroutine read_imgbatch_2

    subroutine read_imgbatch_3( params, build, iptcl, img )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: iptcl
        class(image),      intent(inout) :: img
        type(string) :: stkname
        integer :: ind_in_stk
        call build%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
        if( .not. stkio_r%stk_is_open() )then
            call stkio_r%open(stkname, params%smpd, 'read')
        else if( .not. stkio_r%same_stk(stkname, [params%box,params%box,1]) )then
            call stkio_r%close
            call stkio_r%open(stkname, params%smpd, 'read')
        endif
        call stkio_r%read(ind_in_stk, img)
        call stkio_r%close
    end subroutine read_imgbatch_3

    subroutine new_ptcl_read_session( self, params )
        class(ptcl_read_session), intent(inout) :: self
        class(parameters),        intent(in)    :: params
        call self%stacks%new(params%box, params%smpd, report=params%part == 1)
        if( allocated(self%pass) ) deallocate(self%pass)
        self%cursor = 1
    end subroutine new_ptcl_read_session

    !> The pass's particles in the order its batches will ask for them, for reading ahead
    subroutine expect_ptcl_pass( self, pinds )
        class(ptcl_read_session), intent(inout) :: self
        integer,                  intent(in)    :: pinds(:)
        self%pass   = pinds
        self%cursor = 1
    end subroutine expect_ptcl_pass

    !> build%imgbatch(1:batchlims(2)-batchlims(1)+1) receives the particles pinds(batchlims(1):batchlims(2))
    subroutine read_ptcl_batch( self, params, build, n, pinds, batchlims )
        class(ptcl_read_session), intent(inout) :: self
        class(parameters),        intent(in)    :: params
        class(builder),           intent(inout) :: build
        integer,                  intent(in)    :: n, pinds(n), batchlims(2)
        type(stack_read_counts)   :: c0, c1
        integer(timer_int_kind)   :: t_read
        integer :: nbatch, nexp, i
        if( batchlims(1) < 1 .or. batchlims(2) > n .or. batchlims(1) > batchlims(2) )then
            write(logfhandle,*) 'batchlims: ', batchlims
            write(logfhandle,*) 'n        : ', n
            THROW_HARD('invalid batchlims; read_ptcl_batch')
        endif
        if( L_BENCH_GLOB ) t_read = tic()
        nbatch = batchlims(2) - batchlims(1) + 1
        if( buffer_covers(params, pinds(batchlims(1):batchlims(2))) )then
            do i = batchlims(1), batchlims(2)
                call pass_slab%get(pass_slot(pinds(i) - pass_lo + 1), build%imgbatch(i - batchlims(1) + 1))
            end do
            read_totals%nbuffered = read_totals%nbuffered + int(nbatch,int64)
            if( L_BENCH_GLOB ) read_totals%rt_read = read_totals%rt_read + toc(t_read)
            return
        endif
        c0 = self%stacks%get_counts()
        call read_particles(pinds(batchlims(1):batchlims(2)), .false.)
        ! the next batch: the rest of a whole-pass list, or the next stretch of the expected pass
        if( allocated(self%pass) )then
            nexp = size(self%pass)
            ! a batch off the expected order finds its place in the pass again
            if( self%cursor + nbatch - 1 > nexp .or. self%cursor < 1 )then
                self%cursor = findloc(self%pass, pinds(batchlims(1)), dim=1)
            else if( self%pass(self%cursor) /= pinds(batchlims(1)) )then
                self%cursor = findloc(self%pass, pinds(batchlims(1)), dim=1)
            endif
            if( self%cursor >= 1 .and. self%cursor + nbatch - 1 <= nexp )then
                if( all(self%pass(self%cursor:self%cursor+nbatch-1) == pinds(batchlims(1):batchlims(2))) )then
                    self%cursor = self%cursor + nbatch
                    if( self%cursor <= nexp ) call read_particles(self%pass(self%cursor:min(nexp,self%cursor+nbatch-1)), .true.)
                endif
            endif
        else if( batchlims(2) < n )then
            call read_particles(pinds(batchlims(2)+1:min(n,batchlims(2)+nbatch)), .true.)
        endif
        c1 = self%stacks%get_counts()
        read_totals%nopens     = read_totals%nopens     + c1%nopens     - c0%nopens
        read_totals%nrequests  = read_totals%nrequests  + c1%nrequests  - c0%nrequests
        read_totals%bytes_req  = read_totals%bytes_req  + c1%bytes_req  - c0%bytes_req
        read_totals%bytes_read = read_totals%bytes_read + c1%bytes_read - c0%bytes_read
        read_totals%bytes_prefetched    = read_totals%bytes_prefetched    + c1%bytes_prefetched    - c0%bytes_prefetched
        read_totals%bytes_prefetch_used = read_totals%bytes_prefetch_used + c1%bytes_prefetch_used - c0%bytes_prefetch_used
        read_totals%bytes_advised       = read_totals%bytes_advised       + c1%bytes_advised       - c0%bytes_advised
        if( L_BENCH_GLOB ) read_totals%rt_read = read_totals%rt_read + toc(t_read)

    contains

        !> reads the particles into build%imgbatch, or (ahead) opens their stacks and reads them into the window
        subroutine read_particles( pinds_here, l_ahead )
            integer, intent(in) :: pinds_here(:)
            logical, intent(in) :: l_ahead
            type(string), allocatable :: stknames(:)
            integer,      allocatable :: inds_in_stk(:)
            integer :: i
            allocate(stknames(size(pinds_here)), inds_in_stk(size(pinds_here)))
            do i = 1, size(pinds_here)
                call build%spproj%get_stkname_and_ind(params%oritype, pinds_here(i), stknames(i), inds_in_stk(i))
            end do
            if( l_ahead )then
                call self%stacks%prefetch(stknames, inds_in_stk)
            else
                call self%stacks%read(stknames, inds_in_stk, build%imgbatch(1:size(pinds_here)))
            endif
            call stknames(:)%kill
        end subroutine read_particles

    end subroutine read_ptcl_batch

    !> The particle buffer's budget in bytes: half the task memory a queue system grants (job_memory_per_task);
    !! otherwise the node's physical memory over twice the processes it runs (the parts and their master, or one
    !! process in shared memory)
    integer(int64) function ptcl_buffer_budget( params, build ) result( budget )
        class(parameters), intent(in) :: params
        class(builder),    intent(in) :: build
        type(string)      :: qsys
        integer(int64)    :: total
        integer           :: nprocs
        if( build%spproj%compenv%get_noris() > 0 )then
            if( build%spproj%compenv%isthere(1, 'qsys_name') .and. build%spproj%compenv%isthere(1, 'job_memory_per_task') )then
                qsys = build%spproj%compenv%get_str(1, 'qsys_name')
                if( .not. (qsys .eq. 'local' .or. qsys .eq. 'coarray') )then
                    budget = int(build%spproj%compenv%get_int(1, 'job_memory_per_task'),int64) * 1024_int64**2 / 2_int64
                    return
                endif
            endif
        endif
        nprocs = 1
        if( params%nparts > 1 ) nprocs = params%nparts + 1
        total  = get_node_mem_total_bytes()
        budget = max(0_int64, total) / int(2 * nprocs,int64)
    end function ptcl_buffer_budget

    !> Reads the pass's share pinds(1:n) once, in order through reads, into the particle buffer; then every
    !! batch of these particles is copied from it. Returns BUFFER_FILLED, or BUFFER_OFF (no budget or no particle),
    !! BUFFER_OVER_BUDGET or BUFFER_NO_MEMORY with the buffer left empty and the batches read from the stacks.
    integer function fill_pass_buffer( params, build, reads, n, pinds, budget ) result( status )
        class(parameters),        intent(in)    :: params
        class(builder),           intent(inout) :: build
        class(ptcl_read_session), intent(inout) :: reads
        integer,                  intent(in)    :: n, pinds(n)
        integer(int64),           intent(in)    :: budget
        integer :: slab_status, i, k, nb, ldim(3)
        call release_pass_buffer
        status = BUFFER_OFF
        if( n < 1 .or. budget <= 0_int64 ) return
        status = BUFFER_FILLED
        if( .not. allocated(build%imgbatch) ) THROW_HARD('the image batch is not allocated; fill_pass_buffer')
        ldim = build%imgbatch(1)%get_ldim()
        if( ldim(1) /= params%box ) THROW_HARD('the image batch is not at the particle box; fill_pass_buffer')
        call pass_slab%new(params%box, n, budget, slab_status)
        select case(slab_status)
            case(SLAB_OVER_BUDGET)
                status = BUFFER_OVER_BUDGET
                return
            case(SLAB_NO_MEMORY)
                status = BUFFER_NO_MEMORY
                return
        end select
        pass_lo = minval(pinds)
        allocate(pass_slot(maxval(pinds) - pass_lo + 1), source=0)
        nb = size(build%imgbatch)
        do i = 1, n, nb
            call reads%read(params, build, n, pinds, [i, min(n, i + nb - 1)])
            do k = i, min(n, i + nb - 1)
                call pass_slab%put(k, build%imgbatch(k - i + 1))
                pass_slot(pinds(k) - pass_lo + 1) = k
            end do
        end do
        call pass_slab%seal
        pass_box     = params%box
        l_pass_ready = .true.
    end function fill_pass_buffer

    subroutine release_pass_buffer
        call pass_slab%kill
        if( allocated(pass_slot) ) deallocate(pass_slot)
        pass_lo      = 0
        pass_box     = 0
        l_pass_ready = .false.
    end subroutine release_pass_buffer

    integer(int64) function pass_buffer_bytes()
        pass_buffer_bytes = pass_slab%bytes()
    end function pass_buffer_bytes

    function pass_buffer_status_name( status ) result( name )
        integer, intent(in) :: status
        character(len=:), allocatable :: name
        select case(status)
            case(BUFFER_FILLED);      name = 'filled'
            case(BUFFER_OFF);         name = 'off, batch reads'
            case(BUFFER_OVER_BUDGET); name = 'over budget, batch reads'
            case default;             name = 'allocation failed, batch reads'
        end select
    end function pass_buffer_status_name

    !> Does the particle buffer hold every one of these particles at this box
    logical function buffer_covers( params, pinds )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: pinds(:)
        integer :: i, j
        buffer_covers = .false.
        if( .not. l_pass_ready ) return
        if( params%box /= pass_box ) return
        do i = 1, size(pinds)
            j = pinds(i) - pass_lo + 1
            if( j < 1 .or. j > size(pass_slot) ) return
            if( pass_slot(j) < 1 ) return
        end do
        buffer_covers = .true.
    end function buffer_covers

    subroutine kill_ptcl_read_session( self )
        class(ptcl_read_session), intent(inout) :: self
        call self%stacks%kill
        if( allocated(self%pass) ) deallocate(self%pass)
        self%cursor = 1
    end subroutine kill_ptcl_read_session

    !> The reconstruction observation of both backends. Cropped: noise-normalized at the native box against lmsk,
    !! Fourier-cropped to obs, back in real space, then tapered at the cropped box if l_taper; uncropped: tapered,
    !! then normalized, as in the fused gridding routine. Cropping and tapering do not commute, so all backends use it.
    subroutine prep_rec_observation( ptcl_img, lmsk, obs, l_taper )
        class(image), intent(inout) :: ptcl_img
        logical,      intent(in)    :: lmsk(:,:,:)
        class(image), intent(inout) :: obs
        logical,      intent(in)    :: l_taper
        integer :: ldim(3), ldim_obs(3)
        real    :: sdev_noise, edge_mean
        ldim     = ptcl_img%get_ldim()
        ldim_obs = obs%get_ldim()
        if( ldim_obs(1) > ldim(1) ) THROW_HARD('observation box exceeds the particle box; prep_rec_observation')
        if( ldim_obs(1) < ldim(1) )then
            call ptcl_img%norm_noise_fft_clip_shift(lmsk, obs, [0.,0.])
            call obs%ifft
            if( l_taper ) call obs%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
        else
            call obs%copy(ptcl_img)
            if( l_taper ) call obs%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
            call obs%norm_noise(lmsk, sdev_noise)
        endif
    end subroutine prep_rec_observation

end module simple_matcher_ptcl_io
