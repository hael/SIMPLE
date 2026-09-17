!@descr: flex_pca plane cache: the full-box prep's padded transform, restricted to the box_crop grid, kept on disk per particle
!!
!! The full-box prep of a particle (noise-normalise, taper, gridding-pad to boxpd, FFT scaled by
!! 1/boxpd^2) is followed by gen_fplane4rec, which only ever reads the padded transform at the
!! frequencies of the box_crop grid. Both padded grids (boxpd and box_croppd) span the same
!! physical field of view, so the Fourier index h means the same spatial frequency on either and
!! the 1/n^2 scaling makes the coefficients identical. The cache therefore stores, per particle,
!! the complex block of the boxpd transform at |h| <= box_croppd/2 -- exactly the cmat a
!! box_croppd image would carry -- and a cached run injects that block into the padded heap image
!! and continues through gen_fplane4rec unchanged. No second normalisation, taper or crop happens
!! on the cached path, so its planes equal the full-box planes to floating-point rounding.
!! (The shared particle cache stores a real-space Fourier-cropped image instead; re-tapering and
!! padding that 64-pixel image is a different prep: rank cut 7 -> 13-15 on the 20k harness.)
!!
!! Layout: one direct-access unformatted file under cache_dir (or SIMPLE_PTCL_CACHE_DIR, else the
!! run directory), record 1 = header key, record p+1 = the block of project row p. Only the rows
!! of the run's selection are written; a header mismatch (rows, box, box_crop, smpd, selection
!! hash) rebuilds. Built once by the process that owns the run (shared memory or the distributed
!! master) before any pass; workers open it read-only.
module simple_flex_pca_plane_cache
use simple_core_module_api
use simple_builder,         only: builder
use simple_parameters,      only: parameters
use simple_image,           only: image
use simple_matcher_3Drec,   only: init_rec, cleanup_rec_buffers
use simple_matcher_ptcl_io, only: discrete_read_imgbatch, prepimgbatch
implicit none
private
#include "simple_local_flags.inc"

public :: plane_cache_in_use, plane_cache_ensure, plane_cache_read_batch, plane_cache_fill, plane_cache_close

integer(8),       parameter :: PLANE_CACHE_MAGIC   = 7134151962_8
integer(8),       parameter :: PLANE_CACHE_VERSION = 1_8
integer,          parameter :: NHEAD               = 9
character(len=*), parameter :: CACHE_DIR_ENV       = 'SIMPLE_PTCL_CACHE_DIR'

complex, allocatable :: blk(:,:,:,:)      !< (nph, npd, 1, MAXIMGBATCHSZ): one box_croppd cmat per batch slot
integer      :: funit    = 0              !< open direct-access unit, 0 = closed
integer      :: nph      = 0              !< fdim(box_croppd)
integer      :: npd      = 0              !< box_croppd
integer      :: reclen   = 0
logical      :: l_probed = .false.
logical      :: l_avail  = .false.
type(string) :: fname_glob

contains

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
        type(string) :: fname, dir
        integer(8)   :: h
        dir = plane_cache_dir(params)
        h = 1469598103934665603_8
        if( allocated(CWD_GLOB) ) call fold_str(h, CWD_GLOB)
        call fold_str(h, params%projfile%to_char())
        fname = string(dir%to_char()//'/flex_pca_planes_b'//int2str(params%box_crop)//'_'//&
            &int2str(int(abs(mod(h, 1000000007_8))))//'.bin')
    end function plane_cache_fname

    !> FNV-1a style fold of a string into a 64-bit hash
    subroutine fold_str( h, s )
        integer(8),       intent(inout) :: h
        character(len=*), intent(in)    :: s
        integer :: i
        do i = 1, len_trim(s)
            h = ieor(h, int(ichar(s(i:i)), 8))
            h = h * 1099511628211_8
        end do
    end subroutine fold_str

    !> the header record: magic, version, project rows, box, box_crop, smpd (x1e6), selection count, selection hash,
    !! highest project row written (so a truncated file can be recognised from its size)
    function header_key( params, build, pinds, nptcls ) result( key )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nptcls
        integer(8) :: key(NHEAD), h
        integer    :: i
        h = 1469598103934665603_8
        do i = 1, nptcls
            h = ieor(h, int(pinds(i), 8))
            h = h * 1099511628211_8
        end do
        key = [PLANE_CACHE_MAGIC, PLANE_CACHE_VERSION, int(build%spproj_field%get_noris(), 8), &
            &int(params%box, 8), int(params%box_crop, 8), int(nint(params%smpd * 1.e6), 8), &
            &int(nptcls, 8), h, int(maxval(pinds(1:nptcls)), 8)]
    end function header_key

    logical function plane_cache_applicable( params )
        class(parameters), intent(in) :: params
        plane_cache_applicable = params%l_cache .and. params%box_crop < params%box
    end function plane_cache_applicable

    ! ---------------------------------------------------------------- probe / open

    !> True when cache=yes, box_crop < box, and a file with a matching header exists (probed once).
    !! Does not need the selection: a worker adopts whatever the master built.
    logical function plane_cache_in_use( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer(8) :: key(NHEAD), fsize
        integer    :: io, u
        if( .not. plane_cache_applicable(params) )then
            plane_cache_in_use = .false.
            return
        endif
        if( .not. l_probed )then
            l_probed   = .true.
            l_avail    = .false.
            fname_glob = plane_cache_fname(params)
            if( file_exists(fname_glob) )then
                call set_geometry(params)
                call fopen(u, fname_glob, 'old', 'read', iostat=io, access='direct', form='unformatted', recl=reclen)
                if( io == 0 )then
                    read(u, rec=1, iostat=io) key
                    if( io == 0 )then
                        l_avail = key(1) == PLANE_CACHE_MAGIC .and. key(2) == PLANE_CACHE_VERSION .and. &
                            &key(3) == int(build%spproj_field%get_noris(), 8) .and. key(4) == int(params%box, 8) .and. &
                            &key(5) == int(params%box_crop, 8) .and. key(6) == int(nint(params%smpd * 1.e6), 8)
                        ! a file whose header was written but whose records were not is not a cache
                        inquire(unit=u, size=fsize)
                        if( l_avail ) l_avail = fsize >= (key(9) + 1_8) * int(reclen, 8)
                    endif
                    call fclose(u)
                endif
            endif
        endif
        plane_cache_in_use = l_avail
    end function plane_cache_in_use

    subroutine set_geometry( params )
        class(parameters), intent(in) :: params
        complex :: probe(1)
        integer :: one_len
        npd = params%box_croppd
        nph = npd/2 + 1
        if( allocated(blk) ) deallocate(blk)
        allocate(blk(nph, npd, 1, MAXIMGBATCHSZ), source=cmplx(0.,0.))
        inquire(iolength=one_len) probe
        ! header and particle records share one length; the header is NHEAD integer(8) = 16 complex
        reclen = one_len * nph * npd
        if( reclen < one_len * 2 * NHEAD ) THROW_HARD('plane cache record too small for its header')
    end subroutine set_geometry

    subroutine open_for_read( params )
        class(parameters), intent(in) :: params
        integer :: io
        if( funit /= 0 ) return
        if( npd == 0 ) call set_geometry(params)
        call fopen(funit, fname_glob, 'old', 'read', iostat=io, access='direct', form='unformatted', recl=reclen)
        if( io /= 0 ) THROW_HARD('could not open the plane cache for reading: '//fname_glob%to_char())
    end subroutine open_for_read

    subroutine plane_cache_close()
        if( funit /= 0 ) call fclose(funit)
        funit = 0
        if( allocated(blk) ) deallocate(blk)
    end subroutine plane_cache_close

    ! ---------------------------------------------------------------- build

    !> Build the cache for the run's selection unless a valid one is already there.
    subroutine plane_cache_ensure( params, build, pinds, nptcls )
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
        if( plane_cache_in_use(params, build) )then
            write(logfhandle,'(A)') '>>> FLEX_PCA PLANE CACHE adopted: '//fname_glob%to_char()
            call flush(logfhandle)
            return
        endif
        t0 = tic()
        call set_geometry(params)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA BUILDING PLANE CACHE: ', nptcls, &
            &' particles, box ', params%box, ' -> padded transform block on the box_crop grid (', npd, ')'
        call flush(logfhandle)
        ! full-box prep buffers: padded heap at boxpd, raw image batch at box
        call init_rec(params, build, MAXIMGBATCHSZ, fpls)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        call simple_mkdir(plane_cache_dir(params))
        call fopen(u, fname_glob, 'replace', 'readwrite', iostat=io, access='direct', form='unformatted', recl=reclen)
        if( io /= 0 ) THROW_HARD('could not create the plane cache: '//fname_glob%to_char())
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
                allocate(blk1(nph, npd, 1))
                do k = -npd/2, npd/2 - 1
                    kp = merge(k + 1, k + npd + 1, k >= 0)
                    do h = 0, nph - 1
                        blk1(h+1, kp, 1) = build%img_pad_heap(1)%get_fcomp2D(h, k)
                    end do
                end do
                call probe_img%new([npd, npd, 1], params%smpd_crop)
                call probe_img%set_cmat(blk1)
                call probe_img%set_ft(.true.)
                maxdiff = 0.
                do k = -npd/2, npd/2 - 1
                    do h = 0, nph - 1
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
                do k = -npd/2, npd/2 - 1
                    kp = merge(k + 1, k + npd + 1, k >= 0)
                    do h = 0, nph - 1
                        blk(h+1, kp, 1, i) = build%img_pad_heap(ithr)%get_fcomp2D(h, k)
                    end do
                end do
            end do
            !$omp end parallel do
            do i = 1, batchsz
                write(u, rec=pinds(batchlims(1)+i-1) + 1) blk(:,:,:,i)
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
        l_probed = .true.
        l_avail  = .true.
        write(logfhandle,'(A,F8.1,A,F6.2,A)') '>>> FLEX_PCA PLANE CACHE WRITTEN: '//fname_glob%to_char()// &
            &'  seconds=', toc(t0), '  (', real(nptcls,dp)*real(reclen,dp)/1.d9, ' GB)'
        call flush(logfhandle)
    end subroutine plane_cache_ensure

    ! ---------------------------------------------------------------- read

    !> Read the blocks of batchlims of the pinds list into the batch buffer (slot i = batch position).
    subroutine plane_cache_read_batch( params, n, pinds, batchlims )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: n, pinds(n), batchlims(2)
        integer :: i, io
        call open_for_read(params)
        do i = batchlims(1), batchlims(2)
            read(funit, rec=pinds(i) + 1, iostat=io) blk(:,:,:,i-batchlims(1)+1)
            if( io /= 0 ) THROW_HARD('plane cache read failed at project row '//int2str(pinds(i)))
        end do
    end subroutine plane_cache_read_batch

    !> Load batch slot i into a box_croppd image as its Fourier transform.
    subroutine plane_cache_fill( i, img )
        integer,      intent(in)    :: i
        class(image), intent(inout) :: img
        call img%set_cmat(blk(:,:,:,i))
        call img%set_ft(.true.)
    end subroutine plane_cache_fill

end module simple_flex_pca_plane_cache
