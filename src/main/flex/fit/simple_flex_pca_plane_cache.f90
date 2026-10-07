!@descr: flex_pca disk plane cache with project, geometry, selection and worker-coverage validation
!! cache=yes and box_crop<box only. Record 1 is the validated contract; record p+1 is project row p.
!! The master requires the exact selection. A worker requires matching provenance/geometry and rows covered
!! by the completed master cache. The payload is the box_croppd Fourier block used before gen_fplane4rec.
module simple_flex_pca_plane_cache
!$ use omp_lib, only: omp_get_thread_num
use simple_core_module_api, only: dp, fclose, file_exists, fopen, fplane_type, int2str, logfhandle, maximgbatchsz, &
    &simple_abspath, simple_exception, simple_file_stat, simple_mkdir, stdlen, string, tic, timer_int_kind, toc
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_image,           only: image
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
implicit none
private
#include "simple_local_flags.inc"

public :: flex_pca_plane_cache
public :: plane_cache_contract_header, plane_cache_master_matches, plane_cache_worker_matches

integer(8),       parameter :: PLANE_CACHE_MAGIC   = 7134151962_8
integer(8),       parameter :: PLANE_CACHE_VERSION = 2_8
integer,          parameter :: NHEAD               = 11
integer,          parameter :: IH_NROWS = 5, IH_SMPD = 8, IH_NSEL = 9, IH_HIGHROW = 11
character(len=*), parameter :: CACHE_DIR_ENV       = 'SIMPLE_PTCL_CACHE_DIR'

type :: flex_pca_plane_cache
    private
    complex, allocatable :: blk(:,:,:,:)  !< (nph, npd, 1, MAXIMGBATCHSZ)
    integer      :: funit    = 0
    integer      :: nph      = 0
    integer      :: npd      = 0
    integer      :: reclen   = 0
    logical      :: l_probed = .false.
    logical      :: l_avail  = .false.
    type(string) :: fname
  contains
    procedure, public :: new      => cache_new
    procedure, public :: kill     => cache_kill
    procedure, public :: pristine => cache_pristine
    procedure, public :: available => cache_available
    procedure, public :: ensure    => plane_cache_ensure
    procedure, public :: adopt     => plane_cache_adopt
    procedure, public :: read_batch => plane_cache_read_batch
    procedure, public :: fill      => plane_cache_fill
end type flex_pca_plane_cache

contains

    subroutine cache_new( self, fname )
        class(flex_pca_plane_cache), intent(inout) :: self
        type(string),                intent(in)    :: fname
        call self%kill
        self%fname    = fname
        self%l_probed = .true.
    end subroutine cache_new

    subroutine cache_kill( self )
        class(flex_pca_plane_cache), intent(inout) :: self
        if( self%funit /= 0 ) call fclose(self%funit)
        self%funit = 0
        if( allocated(self%blk) ) deallocate(self%blk)
        call self%fname%kill
        self%nph = 0; self%npd = 0; self%reclen = 0
        self%l_probed = .false.; self%l_avail = .false.
    end subroutine cache_kill

    logical function cache_pristine( self )
        class(flex_pca_plane_cache), intent(in) :: self
        cache_pristine = self%funit == 0 .and. self%nph == 0 .and. self%npd == 0 .and. &
            &self%reclen == 0 .and. .not. self%l_probed .and. .not. self%l_avail .and. &
            &.not. allocated(self%blk) .and. self%fname%is_blank()
    end function cache_pristine

    logical function cache_available( self )
        class(flex_pca_plane_cache), intent(in) :: self
        cache_available = self%l_avail
    end function cache_available

    ! ---------------------------------------------------------------- naming and key

    function plane_cache_dir( params ) result( dir )
        class(parameters), intent(in) :: params
        type(string)          :: dir
        character(len=STDLEN) :: envdir
        integer               :: envlen, envstat
        if( .not. params%cache_dir%is_blank() )then
            dir = params%cache_dir
        else
            call get_environment_variable(CACHE_DIR_ENV, envdir, envlen, envstat)
            if( envstat == 0 .and. envlen > 0 )then
                dir = string(envdir(:envlen))
            else
                dir = string('.')
            endif
        endif
    end function plane_cache_dir

    function plane_cache_fname( params ) result( fname )
        class(parameters), intent(in) :: params
        type(string) :: fname, dir, project_path
        integer(8)   :: h
        dir = plane_cache_dir(params)
        project_path = simple_abspath(params%projfile)
        h = 1469598103934665603_8
        call fold_str(h, project_path%to_char())
        fname = string(dir%to_char()//'/flex_pca_planes_b'//int2str(params%box_crop)//'_'//&
            &int2str(int(abs(mod(h, 1000000007_8))))//'.bin')
        call project_path%kill
    end function plane_cache_fname

    !> FNV-1a style fold of a string into a 64-bit hash
    pure subroutine fold_str( h, s )
        integer(8),       intent(inout) :: h
        character(len=*), intent(in)    :: s
        integer :: i
        do i = 1, len_trim(s)
            h = ieor(h, int(ichar(s(i:i)), 8))
            h = h * 1099511628211_8
        end do
    end subroutine fold_str

    pure integer(8) function selection_hash( pinds ) result( h )
        integer, intent(in) :: pinds(:)
        integer :: i
        h = 1469598103934665603_8
        do i = 1, size(pinds)
            h = ieor(h, int(pinds(i), 8))
            h = h * 1099511628211_8
        end do
    end function selection_hash

    !> Encode the persistent header contract from values independent of the cache implementation.
    function plane_cache_contract_header( project_path, project_mtime, project_rows, box, box_crop, smpd, pinds ) &
        &result( key )
        character(len=*), intent(in) :: project_path
        integer(8),       intent(in) :: project_mtime
        integer,          intent(in) :: project_rows, box, box_crop
        real,             intent(in) :: smpd
        integer,          intent(in) :: pinds(:)
        integer(8), allocatable :: key(:)
        integer(8) :: path_hash
        integer :: highrow
        path_hash = 1469598103934665603_8
        call fold_str(path_hash, project_path)
        highrow = 0
        if( size(pinds) > 0 ) highrow = maxval(pinds)
        allocate(key(NHEAD))
        key = [PLANE_CACHE_MAGIC, PLANE_CACHE_VERSION, path_hash, project_mtime, int(project_rows,8), &
            &int(box,8), int(box_crop,8), int(nint(smpd * 1.e6),8), int(size(pinds),8), &
            &selection_hash(pinds), int(highrow,8)]
    end function plane_cache_contract_header

    pure function plane_cache_master_matches( stored, expected ) result( matches )
        integer(8), intent(in) :: stored(:), expected(:)
        logical :: matches
        matches = size(stored) == NHEAD .and. size(expected) == NHEAD
        if( matches ) matches = all(stored == expected)
    end function plane_cache_master_matches

    pure function plane_cache_worker_matches( stored, expected, pinds ) result( matches )
        integer(8), intent(in) :: stored(:), expected(:)
        integer,    intent(in) :: pinds(:)
        logical :: matches
        matches = .false.
        if( size(stored) /= NHEAD .or. size(expected) /= NHEAD .or. size(pinds) < 1 ) return
        if( .not. all(stored(1:IH_SMPD) == expected(1:IH_SMPD)) ) return
        if( stored(IH_NSEL) < 1_8 .or. stored(IH_NSEL) > stored(IH_NROWS) ) return
        if( stored(IH_HIGHROW) < 1_8 .or. stored(IH_HIGHROW) > stored(IH_NROWS) ) return
        if( int(size(pinds),8) > stored(IH_NSEL) ) return
        if( any(pinds < 1) .or. maxval(pinds) > stored(IH_HIGHROW) ) return
        matches = .true.
    end function plane_cache_worker_matches

    subroutine project_provenance( params, path, mtime )
        class(parameters), intent(in)  :: params
        type(string),      intent(out) :: path
        integer(8),        intent(out) :: mtime
        integer, allocatable :: statbuf(:)
        integer :: stat_status
        path = simple_abspath(params%projfile)
        stat_status = 0
        call simple_file_stat(path, stat_status, statbuf)
        if( stat_status /= 0 ) THROW_HARD('plane cache could not stat the project file')
        mtime = int(statbuf(10),8)
        deallocate(statbuf)
    end subroutine project_provenance

    !> The master key includes every provenance, geometry and selection field.
    function header_key( params, build, pinds, nptcls ) result( key )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls
        integer(8) :: key(NHEAD), mtime
        integer(8), allocatable :: encoded(:)
        type(string) :: project_path
        call project_provenance(params, project_path, mtime)
        encoded = plane_cache_contract_header(project_path%to_char(), mtime, &
            &build%spproj_field%get_noris(), params%box, params%box_crop, params%smpd, pinds(1:nptcls))
        key = encoded
        deallocate(encoded)
        call project_path%kill
    end function header_key

    logical function plane_cache_applicable( params )
        class(parameters), intent(in) :: params
        plane_cache_applicable = params%l_cache .and. params%box_crop < params%box
    end function plane_cache_applicable

    ! ---------------------------------------------------------------- probe / open

    !> Probe once with the master selection or a worker partition; later calls query that result.
    function plane_cache_probe( self, params, build, pinds, nptcls, worker ) result( in_use )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer, optional, intent(in)    :: pinds(:), nptcls
        logical, optional, intent(in)    :: worker
        integer(8) :: stored(NHEAD), expected(NHEAD), fsize
        integer    :: io, u
        logical    :: in_use, l_worker
        if( .not. plane_cache_applicable(params) )then
            in_use = .false.
            return
        endif
        if( .not. self%l_probed )then
            if( .not. (present(pinds) .and. present(nptcls)) )then
                in_use = .false.
                return
            endif
            if( nptcls < 1 .or. nptcls > size(pinds) )then
                in_use = .false.
                return
            endif
            call self%new(plane_cache_fname(params))
            l_worker = .false.
            if( present(worker) ) l_worker = worker
            if( file_exists(self%fname) )then
                call set_geometry(self, params)
                call fopen(u, self%fname, 'old', 'read', iostat=io, access='direct', &
                    &form='unformatted', recl=self%reclen)
                if( io == 0 )then
                    read(u, rec=1, iostat=io) stored
                    if( io == 0 )then
                        expected = header_key(params, build, pinds, nptcls)
                        if( l_worker )then
                            self%l_avail = plane_cache_worker_matches(stored, expected, pinds(1:nptcls))
                        else
                            self%l_avail = plane_cache_master_matches(stored, expected)
                        endif
                        inquire(unit=u, size=fsize)
                        if( self%l_avail ) self%l_avail = &
                            &fsize >= (stored(IH_HIGHROW) + 1_8) * int(self%reclen,8)
                    endif
                    call fclose(u)
                endif
            endif
        endif
        in_use = self%l_avail
    end function plane_cache_probe

    !> A worker may use the master cache only when its whole partition is covered.
    subroutine plane_cache_adopt( self, params, build, pinds, nptcls )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls
        logical :: adopted
        adopted = plane_cache_probe(self, params, build, pinds, nptcls, worker=.true.)
        if( adopted )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PLANE CACHE adopted by worker: '//self%fname%to_char()
            call flush(logfhandle)
        else if( plane_cache_applicable(params) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PLANE CACHE refused by worker; using native particle reads'
            call flush(logfhandle)
        endif
    end subroutine plane_cache_adopt

    subroutine set_geometry( self, params )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(in) :: params
        complex :: probe(1)
        integer :: one_len
        self%npd = params%box_croppd
        self%nph = self%npd/2 + 1
        if( allocated(self%blk) ) deallocate(self%blk)
        allocate(self%blk(self%nph, self%npd, 1, MAXIMGBATCHSZ), source=cmplx(0.,0.))
        inquire(iolength=one_len) probe
        self%reclen = one_len * self%nph * self%npd
        if( self%reclen < one_len * 2 * NHEAD ) THROW_HARD('plane cache record too small for its header')
    end subroutine set_geometry

    subroutine open_for_read( self, params )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(in) :: params
        integer :: io
        if( self%funit /= 0 ) return
        if( self%npd == 0 ) call set_geometry(self, params)
        call fopen(self%funit, self%fname, 'old', 'read', iostat=io, access='direct', &
            &form='unformatted', recl=self%reclen)
        if( io /= 0 ) THROW_HARD('could not open the plane cache for reading: '//self%fname%to_char())
    end subroutine open_for_read

    ! ---------------------------------------------------------------- build

    !> Build the cache for the run's selection unless a valid one is already there.
    subroutine plane_cache_ensure( self, params, build, pinds, nptcls )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls
        type(fplane_type), allocatable :: fpls(:)
        type(image) :: probe_img, tmp_img
        complex, allocatable :: blk1(:,:,:)
        integer(8)  :: key(NHEAD)
        integer(timer_int_kind) :: t0
        integer :: ibatch, batchlims(2), batchsz, i, ithr, h, k, kp, io, u
        real    :: maxdiff
        if( .not. plane_cache_applicable(params) ) return
        if( plane_cache_probe(self, params, build, pinds, nptcls) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PLANE CACHE adopted: '//self%fname%to_char()
            call flush(logfhandle)
            return
        endif
        t0 = tic()
        call set_geometry(self, params)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA BUILDING PLANE CACHE: ', nptcls, &
            &' particles, box ', params%box, ' -> padded transform block on the box_crop grid (', self%npd, ')'
        call flush(logfhandle)
        ! full-box prep buffers: padded heap at boxpd, raw image batch at box
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call simple_mkdir(plane_cache_dir(params))
        call fopen(u, self%fname, 'replace', 'readwrite', iostat=io, access='direct', &
            &form='unformatted', recl=self%reclen)
        if( io /= 0 ) THROW_HARD('could not create the plane cache: '//self%fname%to_char())
        key = 0_8
        write(u, rec=1) key          ! placeholder: the real header is written when every record is in
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch + MAXIMGBATCHSZ - 1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            if( ibatch == 1 )then
                ! layout self-check on a COPY of the first particle (the prep tapers its input in place):
                ! a box_croppd image loaded from the extracted block must return the same components the
                ! padded full-box image returns at every frequency of the block
                call tmp_img%copy(build%imgbatch(1))
                call tmp_img%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(1))
                allocate(blk1(self%nph, self%npd, 1))
                do k = -self%npd/2, self%npd/2 - 1
                    kp = merge(k + 1, k + self%npd + 1, k >= 0)
                    do h = 0, self%nph - 1
                        blk1(h+1, kp, 1) = build%img_pad_heap(1)%get_fcomp2D(h, k)
                    end do
                end do
                call probe_img%new([self%npd, self%npd, 1], params%smpd_crop)
                call probe_img%set_cmat(blk1)
                call probe_img%set_ft(.true.)
                maxdiff = 0.
                do k = -self%npd/2, self%npd/2 - 1
                    do h = 0, self%nph - 1
                        maxdiff = max(maxdiff, abs(probe_img%get_fcomp2D(h,k) - build%img_pad_heap(1)%get_fcomp2D(h,k)))
                    end do
                end do
                call probe_img%kill
                call tmp_img%kill
                deallocate(blk1)
                if( maxdiff > 0. ) THROW_HARD('plane cache layout self-check failed: the re-loaded block &
                    &does not reproduce the padded transform')
            endif
            !$omp parallel do default(shared) private(i,ithr,h,k,kp) schedule(static) proc_bind(close)
            do i = 1, batchsz
                ithr = omp_get_thread_num() + 1
                ! the full-box prep, exactly as prep_imgs4projected_model runs it
                call build%imgbatch(i)%norm_noise_taper_edge_pad_fft(build%lmsk, build%img_pad_heap(ithr))
                ! the box_croppd block, laid out as that grid's own cmat (h >= 0 first index,
                ! negative k wrapped to the top)
                do k = -self%npd/2, self%npd/2 - 1
                    kp = merge(k + 1, k + self%npd + 1, k >= 0)
                    do h = 0, self%nph - 1
                        self%blk(h+1, kp, 1, i) = build%img_pad_heap(ithr)%get_fcomp2D(h, k)
                    end do
                end do
            end do
            !$omp end parallel do
            do i = 1, batchsz
                write(u, rec=pinds(batchlims(1)+i-1) + 1) self%blk(:,:,:,i)
            end do
            if( mod(batchlims(2), 20*MAXIMGBATCHSZ) == 0 .or. batchlims(2) == nptcls )then
                write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA PLANE CACHE PARTICLES: ', batchlims(2), ' / ', nptcls
                call flush(logfhandle)
            endif
        end do
        key = header_key(params, build, pinds, nptcls)
        write(u, rec=1) key
        call fclose(u)
        call cleanup_rec_buffers(build, fpls)
        self%l_probed = .true.
        self%l_avail  = .true.
        write(logfhandle,'(A,F8.1,A,F6.2,A)') '>>> FLEX_PCA PLANE CACHE WRITTEN: '//self%fname%to_char()// &
            &'  seconds=', toc(t0), '  (', real(nptcls,dp)*real(self%reclen,dp)/1.d9, ' GB)'
        call flush(logfhandle)
    end subroutine plane_cache_ensure

    ! ---------------------------------------------------------------- read

    !> Read the blocks of batchlims of the pinds list into the batch buffer (slot i = batch position).
    subroutine plane_cache_read_batch( self, params, n, pinds, batchlims )
        class(flex_pca_plane_cache), intent(inout) :: self
        class(parameters), intent(in) :: params
        integer,           intent(in) :: n, pinds(n), batchlims(2)
        integer :: i, io
        call open_for_read(self, params)
        do i = batchlims(1), batchlims(2)
            read(self%funit, rec=pinds(i) + 1, iostat=io) &
                &self%blk(:,:,:,i-batchlims(1)+1)
            if( io /= 0 ) THROW_HARD('plane cache read failed at project row '//int2str(pinds(i)))
        end do
    end subroutine plane_cache_read_batch

    !> Load batch slot i into a box_croppd image as its Fourier transform.
    subroutine plane_cache_fill( self, i, img )
        class(flex_pca_plane_cache), intent(in) :: self
        integer,      intent(in)    :: i
        class(image), intent(inout) :: img
        call img%set_cmat(self%blk(:,:,:,i))
        call img%set_ft(.true.)
    end subroutine plane_cache_fill

end module simple_flex_pca_plane_cache
