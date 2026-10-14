!@descr: pass-scoped reading of images from many MRC stacks: plain descriptors, coalesced spans, a prefetch window
! A session lives for one pass (a search, a reconstruction, a table). Each read stream opens, reads and closes
! its own stacks inside the parallel region (plain read-only descriptors, no Fortran units); the stack at a
! batch boundary stays open for the next batch; the images needed from a stack are read in spans, one
! positioned request each, and decoded into the images. Told the next batch, the session opens its stacks and
! a producer thread reads them into a bounded window while the caller computes. Policy: section 3.1 of
! particle_io_few_particle_stacks.md.
module simple_stack_read_session
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_image,             only: image
use simple_imghead,           only: MrcImgHead, MRCHEADSZ, MRC_MODE_FLOAT16
use simple_imgfile,           only: simple_float16_to_float32
use simple_discrete_stack_io, only: dstack_io
use simple_syslib,            only: io_read_nstreams, ro_open, ro_close, pread_full, fadvise_willneed, &
    &prefetch_new, prefetch_window, prefetch_submit, prefetch_wait, prefetch_drain, prefetch_free
implicit none

public :: stack_read_session, stack_read_counts, plan_spans
public :: SPAN_GAP_BUDGET, SPAN_SCRATCH_CAP, PREFETCH_WINDOW_BYTES, PREFETCH_MAX_STACKS
private
#include "simple_local_flags.inc"

integer(int64), parameter :: SPAN_GAP_BUDGET       = 2_int64   * 1024_int64**2 !< file bytes a span may bridge between needed images
integer(int64), parameter :: SPAN_SCRATCH_CAP      = 64_int64  * 1024_int64**2 !< raw span scratch of a session, all streams
integer(int64), parameter :: PREFETCH_WINDOW_BYTES = 128_int64 * 1024_int64**2 !< bytes read ahead into the window
integer,        parameter :: PREFETCH_MAX_STACKS   = 64                         !< stacks held open for the next batch
! outcomes of opening and reading a stack, reported after the parallel region
integer, parameter :: ST_OK = 0, ST_OPEN = 1, ST_HEADER = 2, ST_DIMS = 3, ST_INDEX = 4, ST_SHORT = 5, ST_READ = 6

type :: stack_read_counts
    integer(int64) :: nopens              = 0_int64 !< stacks opened
    integer(int64) :: nrequests           = 0_int64 !< read requests, the producer's included
    integer(int64) :: bytes_req           = 0_int64 !< file bytes of the images delivered
    integer(int64) :: bytes_read          = 0_int64 !< file bytes read, the gaps spans bridge and the window included
    integer(int64) :: bytes_prefetched    = 0_int64 !< file bytes the producer read into the window
    integer(int64) :: bytes_prefetch_used = 0_int64 !< window bytes images were decoded from
    integer(int64) :: bytes_advised       = 0_int64 !< file bytes hinted to the kernel beyond the window
end type stack_read_counts

!> an MRC stack open on a plain descriptor, with what its header says
type :: open_stack
    type(string)   :: name
    integer        :: fd = -1, nimgs = 0, mode = -1
    integer(int64) :: data_off = 0_int64, ibytes = 0_int64 !< file offset of image 1, bytes per image
    logical        :: l_signed = .true.
end type open_stack

!> images first..last of held stack istk, read by the producer to window byte boff+1
type :: window_span
    integer        :: istk = 0, first = 0, last = 0, job = -1
    integer(int64) :: boff = 0_int64
end type window_span

type :: raw_scratch
    integer(kind=1), allocatable :: b(:)
end type raw_scratch

type :: stack_read_session
    private
    type(open_stack),  allocatable :: held(:)      !< stacks kept open: the boundary stack and the window's
    type(window_span), allocatable :: wspans(:)    !< the spans in the window
    type(raw_scratch), allocatable :: bufs(:)      !< raw span scratch, one per stream
    integer(kind=1),   pointer     :: wbuf(:) => null() !< the window, owned by the producer
    type(stack_read_counts) :: counts
    type(c_ptr) :: producer = c_null_ptr
    real    :: smpd = 0.
    integer :: box = 0, nstreams = 0, nheld = 0, nwsp = 0
    logical :: l_report = .false., exists = .false.
contains
    procedure :: new
    procedure :: read
    procedure :: prefetch
    procedure :: get_counts
    procedure :: kill
    procedure, private :: setup_streams
    procedure, private :: span_cap
    procedure, private :: find_held
    procedure, private :: release_held
end type stack_read_session

contains

    !> report: the device notice of io_read_nstreams goes to the log (one process of a distributed run)
    subroutine new( self, box, smpd, report )
        class(stack_read_session), intent(inout) :: self
        integer,                   intent(in)    :: box
        real,                      intent(in)    :: smpd
        logical,                   intent(in)    :: report
        call self%kill
        self%box      = box
        self%smpd     = smpd
        self%l_report = report
        allocate(self%held(PREFETCH_MAX_STACKS + 1), self%wspans(0))
        self%exists   = .true.
    end subroutine new

    !> The stream count, from the device of the first stack read, fixes the scratch and the span cap
    subroutine setup_streams( self, stkname )
        class(stack_read_session), intent(inout) :: self
        class(string),             intent(in)    :: stkname
        self%nstreams = io_read_nstreams(stkname, huge(1), self%l_report)
        allocate(self%bufs(self%nstreams))
    end subroutine setup_streams

    !> Images per span within a stream's share of SPAN_SCRATCH_CAP
    integer function span_cap( self, ibytes )
        class(stack_read_session), intent(in) :: self
        integer(int64),            intent(in) :: ibytes
        span_cap = int(max(1_int64, SPAN_SCRATCH_CAP / (int(self%nstreams,int64) * ibytes)))
    end function span_cap

    integer function find_held( self, name )
        class(stack_read_session), intent(in) :: self
        class(string),             intent(in) :: name
        integer :: i
        find_held = 0
        do i = 1, self%nheld
            if( self%held(i)%name .eq. name )then
                find_held = i
                return
            endif
        end do
    end function find_held

    !> Closes the held stacks, all but the one named keep (when given)
    subroutine release_held( self, keep )
        class(stack_read_session), intent(inout) :: self
        type(string), optional,    intent(in)    :: keep
        type(open_stack) :: kept
        integer :: i, rc
        logical :: l_kept
        l_kept = .false.
        do i = 1, self%nheld
            if( present(keep) .and. .not. l_kept )then
                if( self%held(i)%name .eq. keep )then
                    kept   = self%held(i)
                    l_kept = .true.
                    cycle
                endif
            endif
            rc = ro_close(self%held(i)%fd)
        end do
        self%nheld = 0
        if( l_kept )then
            self%nheld   = 1
            self%held(1) = kept
        endif
        self%nwsp = 0
    end subroutine release_held

    !> imgs(i) receives image inds(i) of stack stknames(i), byte for byte as a single read of it
    subroutine read( self, stknames, inds, imgs )
        class(stack_read_session), target, intent(inout) :: self
        type(string),                      intent(in)    :: stknames(:)
        integer,                           intent(in)    :: inds(:)
        class(image),                      intent(inout) :: imgs(:)
        type(string),      allocatable :: gname(:)
        type(open_stack),  allocatable :: gst(:)
        integer,           allocatable :: grp(:), gptr(:), members(:), gheld(:), gstat(:), gerr(:)
        integer,           allocatable :: sp_first(:), sp_last(:), sp_mfrom(:), sp_mto(:)
        logical,           allocatable :: l_fallback(:), l_wok(:)
        type(open_stack) :: st
        type(dstack_io)  :: fallback_reader
        type(c_ptr)      :: src
        character(len=:), allocatable :: what
        integer(int64)   :: nopens, nreq, breq, bread, bused, r, nb, woff
        integer :: n, ng, i, g, k, glast, ldim(3), nthr, j, nsp, isp, iw, m, first, last, rc
        n = size(inds)
        if( n < 1 ) return
        if( .not. self%exists ) THROW_HARD('session not created; read')
        if( size(stknames) /= n .or. size(imgs) < n ) THROW_HARD('stack names, indices and images differ in number; read')
        do i = 1, n
            ldim = imgs(i)%get_ldim()
            if( any(ldim /= [self%box,self%box,1]) ) THROW_HARD('images do not match the session box; read')
        end do
        if( .not. allocated(self%bufs) ) call self%setup_streams(stknames(1))
        ! the window's reads end before its buffer is decoded or its descriptors closed
        allocate(l_wok(self%nwsp), source=.false.)
        do iw = 1, self%nwsp
            nb = int(self%wspans(iw)%last - self%wspans(iw)%first + 1,int64) * self%held(self%wspans(iw)%istk)%ibytes
            r  = prefetch_wait(self%producer, self%wspans(iw)%job)
            l_wok(iw) = r == nb
            if( .not. l_wok(iw) )then
                self%counts%bytes_read       = self%counts%bytes_read       - nb + max(r, 0_int64)
                self%counts%bytes_prefetched = self%counts%bytes_prefetched - nb + max(r, 0_int64)
            endif
        end do
        if( c_associated(self%producer) ) call prefetch_drain(self%producer)
        call group_requests(stknames, inds, ng, gname, grp, gptr, members)
        glast = grp(n)
        allocate(gst(ng), gheld(ng), gstat(ng), gerr(ng), l_fallback(ng))
        gstat = ST_OK
        gerr  = 0
        l_fallback = .false.
        do g = 1, ng
            gheld(g) = self%find_held(gname(g))
        end do
        nopens = 0_int64; nreq = 0_int64; breq = 0_int64; bread = 0_int64; bused = 0_int64
        nthr = min(self%nstreams, ng)
        !$omp parallel do default(shared) private(g,j,st,m,nsp,isp,first,last,nb,iw,woff,r,src,k,i,rc,&
        !$omp& sp_first,sp_last,sp_mfrom,sp_mto) schedule(dynamic,1) proc_bind(close) num_threads(nthr) if(nthr > 1) &
        !$omp& reduction(+:nopens,nreq,breq,bread,bused)
        do g = 1, ng
            j = omp_get_thread_num() + 1
            if( gheld(g) > 0 )then
                st = self%held(gheld(g))
            else
                call open_stack_file(gname(g), self%box, st, gstat(g), gerr(g))
                if( st%fd >= 0 ) nopens = nopens + 1_int64
                if( gstat(g) /= ST_OK ) cycle
            endif
            gst(g) = st
            m = gptr(g+1) - gptr(g)
            if( .not. reads_by_descriptor(st%mode) )then
                l_fallback(g) = .true.
                cycle
            endif
            if( inds(members(gptr(g))) < 1 )then
                gstat(g) = ST_INDEX
                gerr(g)  = inds(members(gptr(g)))
                cycle
            else if( inds(members(gptr(g+1)-1)) > st%nimgs )then
                gstat(g) = ST_INDEX
                gerr(g)  = inds(members(gptr(g+1)-1))
                cycle
            endif
            breq = breq + int(m,int64) * st%ibytes
            allocate(sp_first(m), sp_last(m), sp_mfrom(m), sp_mto(m))
            call plan_spans(inds(members(gptr(g):gptr(g+1)-1)), st%ibytes, self%span_cap(st%ibytes), nsp, &
                &sp_first, sp_last, sp_mfrom, sp_mto)
            do isp = 1, nsp
                first = sp_first(isp)
                last  = sp_last(isp)
                nb    = int(last - first + 1,int64) * st%ibytes
                ! a window span covering this one serves it from the window
                woff = -1_int64
                if( gheld(g) > 0 )then
                    do iw = 1, self%nwsp
                        if( self%wspans(iw)%istk /= gheld(g) .or. .not. l_wok(iw) ) cycle
                        if( self%wspans(iw)%first <= first .and. self%wspans(iw)%last >= last )then
                            woff = self%wspans(iw)%boff + int(first - self%wspans(iw)%first,int64) * st%ibytes
                            exit
                        endif
                    end do
                endif
                if( woff >= 0_int64 )then
                    bused = bused + nb
                else
                    if( allocated(self%bufs(j)%b) )then
                        if( size(self%bufs(j)%b,kind=int64) < nb ) deallocate(self%bufs(j)%b)
                    endif
                    if( .not. allocated(self%bufs(j)%b) ) allocate(self%bufs(j)%b(nb))
                    r = pread_full(st%fd, c_loc(self%bufs(j)%b), nb, st%data_off + int(first - 1,int64) * st%ibytes)
                    nreq  = nreq + 1_int64
                    bread = bread + max(r, 0_int64)
                    if( r < 0_int64 )then
                        gstat(g) = ST_READ
                        gerr(g)  = int(-r)
                        exit
                    else if( r < nb )then
                        gstat(g) = ST_SHORT
                        gerr(g)  = first + int(r / st%ibytes)
                        exit
                    endif
                endif
                do k = gptr(g) - 1 + sp_mfrom(isp), gptr(g) - 1 + sp_mto(isp)
                    i = members(k)
                    if( woff >= 0_int64 )then
                        src = c_loc(self%wbuf(woff + int(inds(i) - first,int64) * st%ibytes + 1_int64))
                    else
                        src = c_loc(self%bufs(j)%b(int(inds(i) - first,int64) * st%ibytes + 1_int64))
                    endif
                    call decode_image(st, src, self%box, imgs(i))
                end do
            end do
            deallocate(sp_first, sp_last, sp_mfrom, sp_mto)
            ! a transient descriptor closes at once, the boundary stack's stays open
            if( gheld(g) == 0 .and. g /= glast )then
                rc = ro_close(st%fd)
                gst(g)%fd = -1
            endif
        end do
        !$omp end parallel do
        self%counts%nopens              = self%counts%nopens              + nopens
        self%counts%nrequests           = self%counts%nrequests           + nreq
        self%counts%bytes_req           = self%counts%bytes_req           + breq
        self%counts%bytes_read          = self%counts%bytes_read          + bread
        self%counts%bytes_prefetch_used = self%counts%bytes_prefetch_used + bused
        do g = 1, ng
            if( gstat(g) == ST_OK ) cycle
            select case(gstat(g))
                case(ST_OPEN);   what = 'cannot open, errno '//int2str(gerr(g))
                case(ST_HEADER); what = 'no complete MRC header, errno '//int2str(gerr(g))
                case(ST_DIMS);   what = 'images of '//int2str(gerr(g))//' pixels, the box is '//int2str(self%box)
                case(ST_INDEX);  what = 'image index '//int2str(gerr(g))//' out of range'
                case(ST_SHORT);  what = 'end of file within image '//int2str(gerr(g))//' (truncated stack)'
                case default;    what = 'read failed, errno '//int2str(gerr(g))
            end select
            write(logfhandle,'(A)') '>>> PARTICLE STACK '//gname(g)%to_char()//': '//what
            THROW_HARD('a particle stack could not be read; read')
        end do
        ! modes the descriptor path does not decode keep the image file reader
        do g = 1, ng
            if( .not. l_fallback(g) ) cycle
            if( gst(g)%fd >= 0 ) rc = ro_close(gst(g)%fd)
            call fallback_reader%new(self%smpd, self%box)
            do k = gptr(g), gptr(g+1) - 1
                i = members(k)
                call fallback_reader%read(gname(g), inds(i), imgs(i))
            end do
            call fallback_reader%kill
            self%counts%nopens     = self%counts%nopens     + 1_int64
            self%counts%nrequests  = self%counts%nrequests  + int(gptr(g+1) - gptr(g),int64)
            self%counts%bytes_req  = self%counts%bytes_req  + int(gptr(g+1) - gptr(g),int64) * gst(g)%ibytes
            self%counts%bytes_read = self%counts%bytes_read + int(gptr(g+1) - gptr(g),int64) * gst(g)%ibytes
        end do
        ! the boundary stack stays open for the next batch, the other held stacks close
        if( gheld(glast) > 0 .and. .not. l_fallback(glast) )then
            call self%release_held(keep=gname(glast))
        else
            call self%release_held()
            if( .not. l_fallback(glast) .and. gst(glast)%fd >= 0 )then
                self%nheld   = 1
                self%held(1) = gst(glast)
            endif
        endif
        call gname(:)%kill
    end subroutine read

    !> The images the next batch will ask for: their stacks are opened, and a producer thread reads their
    !! spans into the window (WILLNEED hints beyond it) while the caller computes. A wrong guess costs reads,
    !! never bytes: read only serves spans the window covers.
    subroutine prefetch( self, stknames, inds )
        class(stack_read_session), target, intent(inout) :: self
        type(string),                      intent(in)    :: stknames(:)
        integer,                           intent(in)    :: inds(:)
        type(string),      allocatable :: gname(:)
        type(open_stack),  allocatable :: gst(:)
        type(window_span), allocatable :: wsp(:)
        integer,           allocatable :: grp(:), gptr(:), members(:), gheld(:), gstat(:), gerr(:), gopen(:)
        integer,           allocatable :: sp_first(:), sp_last(:), sp_mfrom(:), sp_mto(:)
        type(c_ptr)    :: wptr
        integer(int64) :: off, nb, need, boff
        integer :: n, ng, g, nopen, io, nsp, isp, m, rc, job, nw
        n = size(inds)
        if( n < 1 .or. .not. allocated(self%bufs) ) return
        if( size(stknames) /= n ) THROW_HARD('stack names and indices differ in number; prefetch')
        if( .not. c_associated(self%producer) ) self%producer = prefetch_new()
        if( .not. c_associated(self%producer) ) return
        call prefetch_drain(self%producer)
        self%nwsp = 0
        call group_requests(stknames, inds, ng, gname, grp, gptr, members)
        allocate(gst(ng), gheld(ng), gstat(ng), gerr(ng), gopen(ng))
        gstat = ST_OK
        nopen = 0
        do g = 1, ng
            gheld(g) = self%find_held(gname(g))
            if( gheld(g) == 0 .and. self%nheld + nopen < size(self%held) )then
                nopen        = nopen + 1
                gopen(nopen) = g
            endif
        end do
        ! the window's stacks open concurrently; one that fails is left to the read, which reports it
        !$omp parallel do default(shared) private(io,g) schedule(dynamic,1) proc_bind(close) &
        !$omp& num_threads(max(1,min(self%nstreams,nopen))) if(nopen > 1)
        do io = 1, nopen
            g = gopen(io)
            call open_stack_file(gname(g), self%box, gst(g), gstat(g), gerr(g))
        end do
        !$omp end parallel do
        do io = 1, nopen
            g = gopen(io)
            if( gst(g)%fd >= 0 ) self%counts%nopens = self%counts%nopens + 1_int64
            if( gstat(g) == ST_OK .and. reads_by_descriptor(gst(g)%mode) )then
                self%nheld = self%nheld + 1
                self%held(self%nheld) = gst(g)
                gheld(g) = self%nheld
            else if( gst(g)%fd >= 0 )then
                rc = ro_close(gst(g)%fd)
            endif
        end do
        ! the window's spans in batch order, the window bytes allocated before any read is queued
        need = 0_int64
        nw   = 0
        do g = 1, ng
            if( gheld(g) == 0 ) cycle
            if( .not. reads_by_descriptor(self%held(gheld(g))%mode) ) cycle
            m = gptr(g+1) - gptr(g)
            nw = nw + m
            need = need + int(m,int64) * 2_int64 * self%held(gheld(g))%ibytes
        end do
        if( nw < 1 ) return
        need = min(need, PREFETCH_WINDOW_BYTES)
        wptr = prefetch_window(self%producer, need)
        if( .not. c_associated(wptr) ) return
        call c_f_pointer(wptr, self%wbuf, [need])
        allocate(wsp(nw))
        boff = 0_int64
        do g = 1, ng
            if( gheld(g) == 0 ) cycle
            associate( st => self%held(gheld(g)) )
                if( .not. reads_by_descriptor(st%mode) ) cycle
                m = gptr(g+1) - gptr(g)
                if( minval(inds(members(gptr(g):gptr(g+1)-1))) < 1 .or. &
                    &maxval(inds(members(gptr(g):gptr(g+1)-1))) > st%nimgs ) cycle
                allocate(sp_first(m), sp_last(m), sp_mfrom(m), sp_mto(m))
                call plan_spans(inds(members(gptr(g):gptr(g+1)-1)), st%ibytes, self%span_cap(st%ibytes), nsp, &
                    &sp_first, sp_last, sp_mfrom, sp_mto)
                do isp = 1, nsp
                    nb  = int(sp_last(isp) - sp_first(isp) + 1,int64) * st%ibytes
                    off = st%data_off + int(sp_first(isp) - 1,int64) * st%ibytes
                    job = -1
                    if( boff + nb <= size(self%wbuf,kind=int64) ) job = prefetch_submit(self%producer, st%fd, off, nb, &
                        &c_loc(self%wbuf(boff + 1_int64)))
                    if( job >= 0 )then
                        self%nwsp = self%nwsp + 1
                        wsp(self%nwsp) = window_span(gheld(g), sp_first(isp), sp_last(isp), job, boff)
                        boff = boff + nb
                        self%counts%nrequests        = self%counts%nrequests        + 1_int64
                        self%counts%bytes_read       = self%counts%bytes_read       + nb
                        self%counts%bytes_prefetched = self%counts%bytes_prefetched + nb
                    else
                        if( fadvise_willneed(st%fd, off, nb) == 0 ) self%counts%bytes_advised = self%counts%bytes_advised + nb
                    endif
                end do
                deallocate(sp_first, sp_last, sp_mfrom, sp_mto)
            end associate
        end do
        call move_alloc(wsp, self%wspans)
        call gname(:)%kill
    end subroutine prefetch

    !> Groups a request by stack file in order of first appearance (stack rows sharing a file share a group);
    !! members(gptr(g):gptr(g+1)-1) are group g's requests by index in the stack; grp(i) is request i's group
    subroutine group_requests( stknames, inds, ng, gname, grp, gptr, members )
        type(string),              intent(in)  :: stknames(:)
        integer,                   intent(in)  :: inds(:)
        integer,                   intent(out) :: ng
        type(string), allocatable, intent(out) :: gname(:)
        integer,      allocatable, intent(out) :: grp(:), gptr(:), members(:)
        integer, allocatable :: gfill(:)
        integer :: n, i, g, k, k2, mtmp
        n = size(inds)
        allocate(grp(n), gname(n))
        ng = 0
        do i = 1, n
            g = 0
            if( ng > 0 )then
                if( stknames(i) .eq. gname(ng) )then
                    g = ng
                else
                    do k = 1, ng - 1
                        if( stknames(i) .eq. gname(k) )then
                            g = k
                            exit
                        endif
                    end do
                endif
            endif
            if( g == 0 )then
                ng        = ng + 1
                gname(ng) = stknames(i)
                g         = ng
            endif
            grp(i) = g
        end do
        ! insertion sort by index in the stack: runs usually arrive sorted
        allocate(gptr(ng+1), gfill(ng), source=0)
        do i = 1, n
            gfill(grp(i)) = gfill(grp(i)) + 1
        end do
        gptr(1) = 1
        do g = 1, ng
            gptr(g+1) = gptr(g) + gfill(g)
        end do
        gfill = gptr(1:ng)
        allocate(members(n))
        do i = 1, n
            members(gfill(grp(i))) = i
            gfill(grp(i)) = gfill(grp(i)) + 1
        end do
        do g = 1, ng
            do k = gptr(g) + 1, gptr(g+1) - 1
                mtmp = members(k)
                k2   = k - 1
                do while( k2 >= gptr(g) )
                    if( inds(members(k2)) <= inds(mtmp) ) exit
                    members(k2+1) = members(k2)
                    k2 = k2 - 1
                end do
                members(k2+1) = mtmp
            end do
        end do
    end subroutine group_requests

    !> Real-space pixel modes the descriptor path decodes: bytes, 16-bit integers, float32 and float16
    pure logical function reads_by_descriptor( mode )
        integer, intent(in) :: mode
        reads_by_descriptor = any(mode == [0, 1, 2, 6, MRC_MODE_FLOAT16])
    end function reads_by_descriptor

    !> Opens a stack read-only and decodes its header once; st%fd stays -1 when the open failed. Thread safe.
    subroutine open_stack_file( name, box, st, stat, err )
        type(string),     intent(in)  :: name
        integer,          intent(in)  :: box
        type(open_stack), intent(out) :: st
        integer,          intent(out) :: stat, err
        integer(kind=1), target :: hbytes(MRCHEADSZ)
        type(MrcImgHead) :: hed
        integer(int64)   :: r
        integer          :: dims(3)
        stat    = ST_OK
        err     = 0
        st%name = name
        st%fd   = ro_open(name)
        if( st%fd < 0 )then
            stat  = ST_OPEN
            err   = -st%fd
            st%fd = -1
            return
        endif
        r = pread_full(st%fd, c_loc(hbytes), int(MRCHEADSZ,int64), 0_int64)
        if( r /= int(MRCHEADSZ,int64) )then
            stat = ST_HEADER
            err  = int(max(0_int64, -r))
            return
        endif
        call hed%new
        call hed%decode(hbytes)
        dims        = hed%getDims()
        st%mode     = hed%getMode()
        st%l_signed = hed%pixIsSigned()
        st%nimgs    = dims(3)
        st%data_off = int(hed%firstDataByte(),int64) - 1_int64
        if( dims(1) /= box .or. dims(2) /= box )then
            stat = ST_DIMS
            err  = dims(1)
            call hed%kill
            return
        endif
        st%ibytes = int(box,int64)**2 * int(hed%bytesPerPix(),int64)
        call hed%kill
    end subroutine open_stack_file

    !> One image from its bytes in the file into img, padding zeroed: the values and bits rSlices gives
    subroutine decode_image( st, src, box, img )
        type(open_stack), intent(in)    :: st
        type(c_ptr),      intent(in)    :: src
        integer,          intent(in)    :: box
        class(image),     intent(inout) :: img
        real(c_float),      pointer :: rmat(:,:,:), f32(:,:)
        integer(c_int16_t), pointer :: i16(:,:)
        integer(c_int8_t),  pointer :: i8(:,:)
        integer :: j
        call img%get_rmat_ptr_padded(rmat)
        select case(st%mode)
            case(2)
                call c_f_pointer(src, f32, [box,box])
                rmat(1:box,1:box,1) = f32
            case(MRC_MODE_FLOAT16)
                call c_f_pointer(src, i16, [box,box])
                do j = 1, box
                    call simple_float16_to_float32(c_loc(i16(1,j)), c_loc(rmat(1,j,1)), int(box,c_size_t))
                end do
            case(1, 6)
                call c_f_pointer(src, i16, [box,box])
                if( st%l_signed )then
                    rmat(1:box,1:box,1) = real(i16)
                else
                    rmat(1:box,1:box,1) = real(iand(int(i16,kind=4), 65535))
                endif
            case(0)
                call c_f_pointer(src, i8, [box,box])
                if( st%l_signed )then
                    rmat(1:box,1:box,1) = real(i8)
                else
                    rmat(1:box,1:box,1) = real(iand(int(i8,kind=4), 255))
                endif
        end select
        rmat(box+1:,:,1) = 0.
        call img%set_ft(.false.)
    end subroutine decode_image

    !> The spans of one stack: inds ascending (repeats allowed), image_bytes per image in the file,
    !! nspan_max images at most per span. A gap is bridged within SPAN_GAP_BUDGET when the span's
    !! over-read stays at most the bytes it needs. Span s reads first(s)..last(s) for inds(mfrom(s):mto(s)).
    subroutine plan_spans( inds, image_bytes, nspan_max, nspans, first, last, mfrom, mto )
        integer,        intent(in)    :: inds(:)
        integer(int64), intent(in)    :: image_bytes
        integer,        intent(in)    :: nspan_max
        integer,        intent(out)   :: nspans
        integer,        intent(inout) :: first(:), last(:), mfrom(:), mto(:)
        integer :: k, v, len, need
        nspans = 0
        if( size(inds) < 1 ) return
        call open_span(1)
        need = 1
        do k = 2, size(inds)
            v = inds(k)
            if( v < last(nspans) ) THROW_HARD('indices not ascending; plan_spans')
            if( v == last(nspans) ) cycle
            len = v - first(nspans) + 1
            if( int(v - last(nspans) - 1,int64) * image_bytes <= SPAN_GAP_BUDGET .and. len <= nspan_max .and. &
                &len - (need + 1) <= need + 1 )then
                last(nspans) = v
                need = need + 1
            else
                mto(nspans) = k - 1
                call open_span(k)
                need = 1
            endif
        end do
        mto(nspans) = size(inds)

    contains

        subroutine open_span( k )
            integer, intent(in) :: k
            nspans        = nspans + 1
            first(nspans) = inds(k)
            last(nspans)  = inds(k)
            mfrom(nspans) = k
        end subroutine open_span

    end subroutine plan_spans

    type(stack_read_counts) function get_counts( self )
        class(stack_read_session), intent(in) :: self
        get_counts = self%counts
    end function get_counts

    subroutine kill( self )
        class(stack_read_session), intent(inout) :: self
        if( c_associated(self%producer) ) call prefetch_free(self%producer)
        self%producer = c_null_ptr
        if( allocated(self%held) )then
            call self%release_held()
            deallocate(self%held)
        endif
        if( allocated(self%wspans) ) deallocate(self%wspans)
        if( allocated(self%bufs)   ) deallocate(self%bufs)
        nullify(self%wbuf)
        self%counts   = stack_read_counts()
        self%nstreams = 0
        self%nheld    = 0
        self%nwsp     = 0
        self%box      = 0
        self%smpd     = 0.
        self%exists   = .false.
    end subroutine kill

end module simple_stack_read_session
