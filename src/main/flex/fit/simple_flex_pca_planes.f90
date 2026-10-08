!@descr: flex_pca plane store: application-owned resident planes and disk-cache session
!! Shared-memory runs retain prepped Fourier planes when the disk cache is in use. Fetches return
!! a held batch by copy because passes modify planes in place; rows fill 0.25*MemAvailable at most.
module simple_flex_pca_planes
use simple_core_module_api, only: dp, fplane_type, logfhandle, simple_exception, tic, timer_int_kind, toc
use simple_builder,              only: builder
use simple_image,                only: image
use simple_parameters,           only: parameters
use simple_flex_pca_plane_cache, only: flex_pca_plane_cache
implicit none
private
#include "simple_local_flags.inc"

public :: flex_plane_store

type :: flex_plane_store
    private
    type(fplane_type), allocatable :: held(:)        !< by project row; held when cmplx_plane is allocated
    type(flex_pca_plane_cache), allocatable :: cache
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
    procedure, public :: ensure_cache      => plane_store_ensure_cache
    procedure, public :: adopt_cache       => plane_store_adopt_cache
    procedure, public :: cache_in_use      => plane_store_cache_in_use
    procedure, public :: read_cache_batch  => plane_store_read_cache_batch
    procedure, public :: fill_cached_image => plane_store_fill_cached_image
    procedure, public :: enable            => plane_store_enable
    procedure, public :: fetch             => plane_store_fetch
    procedure, public :: store             => plane_store_store
    procedure, public :: kill              => plane_store_kill
    procedure, public :: release           => plane_store_release
end type flex_plane_store

contains

    !> Turn the store on for a project of nrows particle rows (shared-memory runs only). Only planes on
    !! the box_crop lattice are worth holding, so the store is tied to the plane cache being in use.
    subroutine plane_store_ensure_cache( self, params, build, pinds, nptcls )
        class(flex_plane_store), intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        integer,                 intent(in)    :: pinds(:), nptcls
        if( .not. allocated(self%cache) ) allocate(self%cache)
        call self%cache%ensure(params, build, pinds, nptcls)
    end subroutine plane_store_ensure_cache

    subroutine plane_store_adopt_cache( self, params, build, pinds, nptcls )
        class(flex_plane_store), intent(inout) :: self
        class(parameters),       intent(in)    :: params
        class(builder),          intent(inout) :: build
        integer,                 intent(in)    :: pinds(:), nptcls
        if( .not. allocated(self%cache) ) allocate(self%cache)
        call self%cache%adopt(params, build, pinds, nptcls)
    end subroutine plane_store_adopt_cache

    logical function plane_store_cache_in_use( self )
        class(flex_plane_store), intent(in) :: self
        plane_store_cache_in_use = .false.
        if( allocated(self%cache) ) plane_store_cache_in_use = self%cache%available()
    end function plane_store_cache_in_use

    subroutine plane_store_read_cache_batch( self, params, n, pinds, batchlims )
        class(flex_plane_store), intent(inout) :: self
        class(parameters),       intent(in)    :: params
        integer,                 intent(in)    :: n, pinds(n), batchlims(2)
        if( .not. allocated(self%cache) ) THROW_HARD('plane cache read requested before cache initialization')
        call self%cache%read_batch(params, n, pinds, batchlims)
    end subroutine plane_store_read_cache_batch

    subroutine plane_store_fill_cached_image( self, i, img )
        class(flex_plane_store), intent(in)    :: self
        integer,                 intent(in)    :: i
        class(image),            intent(inout) :: img
        if( .not. allocated(self%cache) ) THROW_HARD('plane cache fill requested before cache initialization')
        call self%cache%fill(i, img)
    end subroutine plane_store_fill_cached_image

    subroutine plane_store_enable( self, nrows )
        class(flex_plane_store), intent(inout) :: self
        integer, intent(in) :: nrows
        real(dp) :: gb_avail
        if( self%l_on ) return
        if( .not. self%cache_in_use() )then
            write(logfhandle,'(A)') '>>> FLEX_PCA RESIDENT PLANES OFF: planes are held only on the cropped grid &
                &(run with cache=yes to enable)'
            return
        endif
        gb_avail  = mem_available_gb()
        self%gb_budget = 0.25d0 * gb_avail
        if( self%gb_budget <= 0.d0 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA RESIDENT PLANES OFF (no memory budget: MemAvailable &
                &unreadable)'
            return
        endif
        allocate(self%held(nrows))
        self%l_on = .true.
        write(logfhandle,'(A,I0,A,F8.1,A,F8.1,A)') '>>> FLEX_PCA RESIDENT PLANES ON: ', nrows, &
            &' rows addressable, budget ', self%gb_budget, ' GB (MemAvailable ', gb_avail, ' GB)'
        call flush(logfhandle)
    end subroutine plane_store_enable

    !> Fetch rows when the whole batch is resident. The read-path key prevents reuse across
    !! incompatible cached/native preparations.
    subroutine plane_store_fetch( self, rows, fpls, cached, found, sec_read )
        class(flex_plane_store),        intent(inout) :: self
        integer,                        intent(in)    :: rows(:)
        type(fplane_type),              intent(inout) :: fpls(:)
        logical,                        intent(in)    :: cached
        logical,                        intent(out)   :: found
        real(timer_int_kind), optional, intent(inout) :: sec_read
        integer(timer_int_kind) :: t
        integer :: i
        found = .false.
        if( .not. self%l_on ) return
        if( .not. self%l_key_set )then
            self%l_cached_key = cached
            self%l_key_set    = .true.
        elseif( cached .neqv. self%l_cached_key )then
            THROW_HARD('resident planes were prepared from a different read path than requested')
        endif
        if( .not. all_held(self, size(rows), rows) ) return
        t = tic()
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, size(rows)
            fpls(i) = self%held(rows(i))
        end do
        !$omp end parallel do
        if( present(sec_read) ) sec_read = sec_read + toc(t)
        self%nserved = self%nserved + 1
        found = .true.
    end subroutine plane_store_fetch

    logical function all_held( self, batchsz, rows )
        class(flex_plane_store), intent(in) :: self
        integer, intent(in) :: batchsz, rows(batchsz)
        integer :: i
        all_held = .false.
        do i = 1, batchsz
            if( .not. allocated(self%held(rows(i))%cmplx_plane) ) return
        end do
        all_held = .true.
    end function all_held

    !> Copy the freshly prepped planes of a batch into the store, stopping at the budget.
    subroutine plane_store_store( self, rows, fpls )
        class(flex_plane_store), intent(inout) :: self
        integer,                 intent(in)    :: rows(:)
        type(fplane_type),       intent(in)    :: fpls(:)
        real(dp) :: gb
        integer  :: i
        if( .not. self%l_on ) return
        self%nprepped = self%nprepped + 1
        if( self%l_full ) return
        do i = 1, size(rows)
            if( allocated(self%held(rows(i))%cmplx_plane) ) cycle
            gb = plane_gb(fpls(i))
            if( self%gb_held + gb > self%gb_budget )then
                self%l_full = .true.
                write(logfhandle,'(A,I0,A,F8.1,A)') '>>> FLEX_PCA RESIDENT PLANES budget reached: ', &
                    &self%nheld, ' rows held, ', self%gb_held, ' GB; the remaining rows are read every pass'
                call flush(logfhandle)
                return
            endif
            self%held(rows(i)) = fpls(i)
            self%nheld   = self%nheld + 1
            self%gb_held = self%gb_held + gb
        end do
    end subroutine plane_store_store

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
    subroutine plane_store_kill( self )
        class(flex_plane_store), intent(inout) :: self
        integer :: i
        if( self%l_on )then
            write(logfhandle,'(A,I0,A,F8.1,A,I0,A,I0,A)') '>>> FLEX_PCA RESIDENT PLANES: held ', self%nheld, &
                &' rows (', self%gb_held, ' GB), batches served=', self%nserved, ' read+prepped=', self%nprepped
            call flush(logfhandle)
            do i = 1, size(self%held)
                if( allocated(self%held(i)%cmplx_plane) )    deallocate(self%held(i)%cmplx_plane)
                if( allocated(self%held(i)%ctfsq_plane) )    deallocate(self%held(i)%ctfsq_plane)
                if( allocated(self%held(i)%transfer_plane) ) deallocate(self%held(i)%transfer_plane)
            end do
            deallocate(self%held)
        endif
        if( allocated(self%cache) )then
            call self%cache%kill
            deallocate(self%cache)
        endif
        self%l_on = .false.; self%l_full = .false.; self%l_key_set = .false.; self%l_cached_key = .false.
        self%nheld = 0; self%nserved = 0; self%nprepped = 0; self%gb_held = 0.d0; self%gb_budget = 0.d0
    end subroutine plane_store_kill

    !> Release the store once the embedding exists: the resident planes are freed and the disk cache this
    !! run built is deleted (the caller owns the cache: the master or the shared-memory run)
    subroutine plane_store_release( self )
        class(flex_plane_store), intent(inout) :: self
        if( allocated(self%cache) ) call self%cache%delete
        call self%kill
    end subroutine plane_store_release

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
