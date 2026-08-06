!@descr: downscaled particle cache shared by the 2D matcher workflows
!
! Every iteration the samplers draw a fresh ~nsample subset of the particles and
! each selected particle is read at full box and Fourier-cropped to box_crop in
! prepimg4align. Nothing is reused: sample4update_cnt exhausts the lowest
! updatecnt tier first, so consecutive iterations draw near-disjoint subsets and
! sweep the whole dataset before revisiting anything. A partial cache would
! therefore never hit; only a cache covering all particles is worth anything.
!
! What is cacheable is the iteration-independent prefix of prepimg4align: noise
! normalization against the full-box lmsk, forward FFT, and the clip to box_crop.
! The shift is not cacheable (it is re-read from the project every iteration) and
! neither is anything after it, since the mask must follow the shift. The CTF
! phase flip is iteration-independent and commutes with the shift, but is left
! out so the cache does not depend on CTF parameters.
!
! An entry is thus the noise-normalized, Fourier-cropped particle, stored as a
! real-space box_crop image: Fourier crop -> inverse FFT -> forward FFT is an
! exact round trip, so reading an entry back reproduces the same coefficients
! prepimg4align would have computed. At box_crop=88 that is ~30 KB against
! ~379 KB at box=308.
!
! abinitio2D holds box_crop fixed across all its stages, so one cache serves a
! whole run and there is nothing to invalidate between stages.
!
! Both consumers read from it: prob_tab2D, which needs only the alignment product,
! and cluster2D_exec, which additionally restores class averages from the raw
! images. The latter is not a pure I/O substitution -- the edge taper and the
! gridding source grid then live at box_crop rather than at box -- and it requires
! cavger_init_online(cropped_ptcls=.true.), which switches the class-averager work
! images to ldim_croppd, the mask to lmsk_crop, and scales the shift and the CTF
! pixel size to the cropped grid. Note in particular that a cached particle is
! already noise-normalized, so restoration must not normalize it a second time:
! cropping discards most of the noise power, and re-normalizing would rescale the
! particle by the crop factor while leaving ctfsq_plane and sigma2 untouched.
module simple_ptcl_cache
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_discrete_stack_io, only: dstack_io
use simple_stack_io,          only: stack_io
use simple_imgarr_utils,      only: alloc_imgarr, dealloc_imgarr
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, killimgbatch
use simple_syslib,            only: io_read_nstreams, simple_mkdir, dir_exists, simple_file_stat, simple_rename
implicit none
#include "simple_local_flags.inc"

public :: ptcl_cache_in_use, ptcl_cache_assert_ready, ptcl_cache_ensure
public :: ptcl_cache_read_batch, ptcl_cache_reset
private

character(len=*), parameter :: CACHE_FBODY   = 'ptcl_cache_'
character(len=*), parameter :: CACHE_DIR_ENV = 'SIMPLE_PTCL_CACHE_DIR'
character(len=*), parameter :: TMP_TAG       = '_part'

! resolved once per process by ptcl_cache_in_use
logical :: l_probed    = .false.
logical :: l_available = .false.
! cache record index per global particle index, 0 for particles that were not cached.
! Deselected particles are never sampled, so caching them would waste both a full-size
! read at build time and box_crop^2*4 bytes on disk; the map buys that back while
! keeping lookups O(1). It is monotone over cached particles, so a sorted batch of
! pinds still maps to a forward scan of the file.
integer, allocatable :: cache_ind(:)

contains

    !>  Where the cache lives. It is large (nsel * box_crop^2 * 4 bytes), so when the
    !!  project sits on a slow disk it usually wants to be somewhere faster: cache_dir
    !!  on the command line, else the SIMPLE_PTCL_CACHE_DIR environment variable, else
    !!  the execution directory. An empty result means the execution directory.
    function cache_dir( params ) result( dirname )
        class(parameters), intent(in) :: params
        type(string)          :: dirname
        character(len=STDLEN) :: envdir
        integer               :: envlen, envstat
        if( .not. params%cache_dir%is_blank() )then
            dirname = params%cache_dir
            return
        endif
        call get_environment_variable(CACHE_DIR_ENV, envdir, envlen, envstat)
        if( envstat == 0 .and. envlen > 0 )then
            dirname = string(envdir(:envlen))
        else
            dirname = string('')
        endif
    end function cache_dir

    !>  Cache files for this run. The project name is part of the basename because
    !!  cache_dir is explicitly meant to be shared -- pointing several runs at the same
    !!  fast disk must not have them delete or overwrite each other. box_crop is in
    !!  there too so a run with a different crop cannot pick up a stale cache.
    !!
    !!  tmp=.true. returns the name used while the cache is being written. The tag goes
    !!  before the extension, never after: image format is inferred from the extension
    !!  (fname2format), so a name ending in the tag would be rejected as an unsupported
    !!  format the moment anything tried to read its header.
    function cache_fname( params, ext, tmp ) result( fname )
        class(parameters), intent(in) :: params
        character(len=*),  intent(in) :: ext
        logical, optional, intent(in) :: tmp
        type(string) :: fname, dirname, basename_here
        character(len=:), allocatable :: tag
        allocate(tag, source='')
        if( present(tmp) )then
            if( tmp ) deallocate(tag)
            if( tmp ) allocate(tag, source=TMP_TAG)
        endif
        if( params%projname%is_blank() )then
            basename_here = string(CACHE_FBODY//int2str(params%box_crop)//tag//ext)
        else
            basename_here = string(CACHE_FBODY//params%projname%to_char()//&
                &'_'//int2str(params%box_crop)//tag//ext)
        endif
        dirname = cache_dir(params)
        if( dirname%is_blank() )then
            fname = basename_here
        else
            fname = filepath(dirname, basename_here)
        endif
        call dirname%kill
        call basename_here%kill
    end function cache_fname

    function cache_stkname( params ) result( fname )
        class(parameters), intent(in) :: params
        type(string) :: fname
        fname = cache_fname(params, STK_EXT)
    end function cache_stkname

    function cache_keyname( params ) result( fname )
        class(parameters), intent(in) :: params
        type(string) :: fname
        fname = cache_fname(params, TXT_EXT)
    end function cache_keyname

    function cache_idxname( params ) result( fname )
        class(parameters), intent(in) :: params
        type(string) :: fname
        fname = cache_fname(params, '_idx'//BIN_EXT)
    end function cache_idxname

    !>  Persist the particle -> record map. It has to be stored rather than recomputed
    !!  from the current states: if a particle is deselected after the cache is built,
    !!  recomputing would silently shift every later particle onto the wrong record.
    subroutine write_cache_index( idxfile, nptcls )
        class(string), intent(in) :: idxfile
        integer,       intent(in) :: nptcls
        integer :: funit, io_stat
        call fopen(funit, idxfile, 'replace', 'write', io_stat, access='stream', form='unformatted')
        call fileiochk('write_cache_index; simple_ptcl_cache', io_stat)
        write(funit) nptcls
        write(funit) cache_ind(1:nptcls)
        call fclose(funit)
    end subroutine write_cache_index

    !>  Load the map and check it still covers every active particle. A particle that
    !!  became active after the cache was written has no record, and rather than read
    !!  a wrong one we declare the cache stale and fall back to the originals.
    logical function read_cache_index( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: idxfile
        integer      :: funit, io_stat, nptcls_file, nptcls, iptcl
        read_cache_index = .false.
        nptcls  = cache_nptcls(build)
        idxfile = cache_idxname(params)
        if( .not. file_exists(idxfile) )then
            call idxfile%kill
            return
        endif
        call fopen(funit, idxfile, 'old', 'read', io_stat, access='stream', form='unformatted')
        if( io_stat /= 0 )then
            call idxfile%kill
            return
        endif
        read(funit, iostat=io_stat) nptcls_file
        if( io_stat /= 0 .or. nptcls_file /= nptcls )then
            call fclose(funit)
            call idxfile%kill
            return
        endif
        if( allocated(cache_ind) ) deallocate(cache_ind)
        allocate(cache_ind(nptcls), source=0)
        read(funit, iostat=io_stat) cache_ind(1:nptcls)
        call fclose(funit)
        call idxfile%kill
        if( io_stat /= 0 )then
            deallocate(cache_ind)
            return
        endif
        do iptcl = 1, nptcls
            if( build%spproj_field%get_state(iptcl) > 0 .and. cache_ind(iptcl) < 1 )then
                write(logfhandle,'(A)') '>>> PARTICLE CACHE: active particles are missing from it, rebuilding'
                deallocate(cache_ind)
                return
            endif
        end do
        read_cache_index = .true.
    end function read_cache_index

    !>  Number of records in the cache, which is the global particle count: entries are
    !!  indexed by global iptcl so that a distributed worker holding only [fromp,top]
    !!  addresses the same file as the master that wrote it. params%nptcls is not used
    !!  here because it is not guaranteed to carry the same meaning in both roles.
    integer function cache_nptcls( build )
        class(builder), intent(inout) :: build
        cache_nptcls = build%spproj_field%get_noris()
    end function cache_nptcls

    !>  Fold a string into two independent modular rolling hashes. Two are used rather
    !!  than one because a collision here means silently serving stale pixels; together
    !!  they give ~2^62 of key space at a handful of integer ops per character.
    pure subroutine fold_str( h1, h2, str )
        integer(kind=8),  intent(inout) :: h1, h2
        character(len=*), intent(in)    :: str
        integer(kind=8), parameter :: M1 = 2147483647_8, A1 = 16807_8  ! 2^31-1
        integer(kind=8), parameter :: M2 = 2147483629_8, A2 = 48271_8  ! largest prime < 2^31-1
        integer :: i
        do i = 1, len_trim(str)
            h1 = mod(h1 * A1 + int(iachar(str(i:i)), 8), M1)
            h2 = mod(h2 * A2 + int(iachar(str(i:i)), 8), M2)
        end do
        ! separator, so that concatenating fields cannot alias a different split
        h1 = mod(h1 * A1 + 255_8, M1)
        h2 = mod(h2 * A2 + 255_8, M2)
    end subroutine fold_str

    !>  Geometry the cached pixels depend on, all of it already in memory. Cheap enough
    !!  for every rank to recompute on every probe.
    function cache_key_geom( params, build ) result( key )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: key
        key = string('box='//int2str(params%box)//&
            &' box_crop='//int2str(params%box_crop)//&
            &' smpd='//trim(real2str(params%smpd))//&
            &' smpd_crop='//trim(real2str(params%smpd_crop))//&
            &' msk='//trim(real2str(params%msk))//&
            &' oritype='//trim(params%oritype)//&
            &' nptcls='//int2str(cache_nptcls(build))//&
            &' nstks='//int2str(build%spproj%os_stk%get_noris()))
    end function cache_key_geom

    !>  Fingerprint of the particle sources themselves.
    !!
    !!  Every stack is folded in project order together with its particle range and its
    !!  physical size and modification time: a name-only or endpoints-only key would miss
    !!  a stack replaced in place, a changed stack in the middle, or an altered
    !!  particle-to-stack mapping, any of which would serve stale pixels. Size and mtime
    !!  make the check conservative -- a touched file forces a rebuild -- which is the
    !!  right way to be wrong.
    !!
    !!  This costs one stat per stack, so only the rank that decides whether to rebuild
    !!  pays it. See cache_key_matches.
    function cache_key_stkfp( build ) result( key )
        class(builder), intent(inout) :: build
        type(string) :: key, stkname
        integer, allocatable :: statbuf(:)
        integer(kind=8) :: h1, h2
        integer :: nstks, istk, stat_status
        nstks = build%spproj%os_stk%get_noris()
        h1 = 1_8
        h2 = 1_8
        do istk = 1, nstks
            stkname = build%spproj%os_stk%get_str(istk, 'stk')
            call fold_str(h1, h2, stkname%to_char())
            call fold_str(h1, h2, int2str(build%spproj%os_stk%get_fromp(istk)))
            call fold_str(h1, h2, int2str(build%spproj%os_stk%get_top(istk)))
            if( build%spproj%os_stk%isthere(istk, 'nptcls_stk') )then
                call fold_str(h1, h2, int2str(build%spproj%os_stk%get_int(istk, 'nptcls_stk')))
            endif
            call simple_file_stat(stkname, stat_status, statbuf)
            if( stat_status == 0 )then
                call fold_str(h1, h2, int2str(statbuf(8)))   ! size
                call fold_str(h1, h2, int2str(statbuf(10)))  ! mtime
            else
                call fold_str(h1, h2, 'nostat')
            endif
            call stkname%kill
        end do
        if( allocated(statbuf) ) deallocate(statbuf)
        key = string('stkfp='//int2str(int(h1))//'-'//int2str(int(h2)))
    end function cache_key_stkfp

    !>  Compare the key file against this run. The key has two lines: geometry, then the
    !!  source fingerprint.
    !!
    !!  full=.true. recomputes both and is used only by ptcl_cache_ensure, which is the
    !!  rank that decides whether to rebuild. Consumers pass full=.false. and take the
    !!  fingerprint line on trust, because recomputing it would mean every worker
    !!  process stat-ing every stack -- a metadata burst on a parallel filesystem, and
    !!  redundant besides, since the builder validated the same sources moments earlier
    !!  and rewrites the key whenever they change.
    logical function cache_key_matches( params, build, full )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        logical,           intent(in)    :: full
        type(string) :: keyfile, key, stored_geom, stored_fp
        integer      :: funit, io_stat
        cache_key_matches = .false.
        keyfile = cache_keyname(params)
        if( .not. file_exists(keyfile) )then
            call keyfile%kill
            return
        endif
        call fopen(funit, keyfile, 'old', 'read', io_stat)
        if( io_stat /= 0 )then
            call keyfile%kill
            return
        endif
        call stored_geom%readline(funit, io_stat)
        if( io_stat == 0 ) call stored_fp%readline(funit, io_stat)
        call fclose(funit)
        if( io_stat == 0 )then
            key = cache_key_geom(params, build)
            cache_key_matches = stored_geom .eq. key
            call key%kill
            if( cache_key_matches .and. full )then
                key = cache_key_stkfp(build)
                cache_key_matches = stored_fp .eq. key
                if( .not. cache_key_matches ) write(logfhandle,'(A)') &
                    &'>>> PARTICLE CACHE: particle stacks have changed since it was built'
                call key%kill
            endif
        endif
        call stored_geom%kill
        call stored_fp%kill
        call keyfile%kill
    end function cache_key_matches

    subroutine write_cache_key( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: keyfile, key_geom, key_fp
        integer      :: funit, io_stat
        keyfile  = cache_keyname(params)
        key_geom = cache_key_geom(params, build)
        key_fp   = cache_key_stkfp(build)
        call fopen(funit, keyfile, 'replace', 'write', io_stat)
        call fileiochk('write_cache_key; simple_ptcl_cache', io_stat)
        write(funit,'(A)') key_geom%to_char()
        write(funit,'(A)') key_fp%to_char()
        call fclose(funit)
        call keyfile%kill
        call key_geom%kill
        call key_fp%kill
    end subroutine write_cache_key

    !>  Is a usable cache present for this run? Probed once per process. A missing or
    !!  stale cache is not an error here; ptcl_cache_assert_ready is what makes it one
    !!  when the user asked for the cache.
    logical function ptcl_cache_in_use( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: stkname
        if( .not. params%l_cache )then
            ptcl_cache_in_use = .false.
            return
        endif
        if( .not. l_probed )then
            l_probed = .true.
            stkname  = cache_stkname(params)
            l_available = file_exists(stkname) .and. cache_key_matches(params, build, full=.false.)
            if( l_available ) l_available = read_cache_index(params, build)
            call stkname%kill
        endif
        ptcl_cache_in_use = l_available
    end function ptcl_cache_in_use

    !>  Refuse to run uncached when the user asked for the cache.
    !!
    !!  Restoring class averages from cropped particles is deliberately not the same
    !!  preprocessing as restoring from full-size ones -- the edge taper and the gridding
    !!  source grid live at box_crop. If some ranks found the cache and others silently
    !!  fell back, their partial sums would be assembled into the same class averages
    !!  from two different preprocessing paths. That is very easy to hit with a
    !!  node-local cache_dir, so availability is mandatory and uniform, not best-effort.
    subroutine ptcl_cache_assert_ready( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: stkname
        if( .not. params%l_cache ) return
        if( params%box_crop >= params%box ) return  ! nothing to cache, see ptcl_cache_ensure
        if( ptcl_cache_in_use(params, build) ) return
        stkname = cache_stkname(params)
        write(logfhandle,'(A)') '>>> PARTICLE CACHE: expected but not usable: '//stkname%to_char()
        write(logfhandle,'(A)') '>>> If cache_dir is node-local it must be visible to every rank,'
        write(logfhandle,'(A)') '>>> otherwise ranks would mix cropped and full-size class-average restoration.'
        call stkname%kill
        THROW_HARD('cache=yes but no usable particle cache; ptcl_cache_assert_ready')
    end subroutine ptcl_cache_assert_ready

    !>  Forget the probe result, e.g. after building the cache.
    subroutine ptcl_cache_reset
        l_probed    = .false.
        l_available = .false.
        if( allocated(cache_ind) ) deallocate(cache_ind)
    end subroutine ptcl_cache_reset

    !>  Build the cache if it is missing or stale. Single sequential pass over the
    !!  particles in index order, which is also the only wholly sequential read of
    !!  the originals the pipeline ever performs. Master-side, call before workers.
    subroutine ptcl_cache_ensure( params, build )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        type(image), allocatable :: cache_imgs(:)
        type(stack_io) :: stkio_w
        type(string)   :: stkname, dirname, keyfile, tmpstk, tmpidx
        integer, allocatable :: pinds(:)
        integer :: batchsz, nbatches, ibatch, batch_start, batch_end, nbatch, i, iptcl, nptcls, nsel
        integer :: ldim_check(3), nptcls_check
        if( .not. params%l_cache ) return
        if( params%box_crop >= params%box )then
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: box_crop == box, nothing to gain, cache disabled'
            return
        endif
        call ptcl_cache_reset
        ! This is the rank that decides whether to rebuild, so it is the one that pays
        ! for the full source fingerprint. Consumers check only the geometry line and
        ! inherit this verdict, which keeps the per-stack stat sweep to once per
        ! cluster2D invocation instead of once per worker process per iteration.
        stkname = cache_stkname(params)
        if( file_exists(stkname) )then
            if( cache_key_matches(params, build, full=.true.) )then
                if( read_cache_index(params, build) )then
                    write(logfhandle,'(A)') '>>> PARTICLE CACHE: up to date'
                    call ptcl_cache_reset
                    call stkname%kill
                    return
                endif
            endif
        endif
        call stkname%kill
        ! the cache directory may well be on another device, so it need not exist yet
        dirname = cache_dir(params)
        if( .not. dirname%is_blank() )then
            if( .not. dir_exists(dirname) ) call simple_mkdir(dirname)
        endif
        call dirname%kill
        stkname = cache_stkname(params)
        keyfile = cache_keyname(params)
        tmpstk  = cache_fname(params, STK_EXT,          tmp=.true.)
        tmpidx  = cache_fname(params, '_idx'//BIN_EXT,  tmp=.true.)
        ! Invalidate first: the key is the commit record, so from here until it is
        ! rewritten no reader can accept whatever state the cache files are in.
        call del_file(keyfile)
        call del_file(tmpstk)
        call del_file(tmpidx)
        nptcls = cache_nptcls(build)
        ! Only active particles are cached: a deselected one is never sampled, so
        ! caching it would cost a full-size read now and box_crop^2*4 bytes forever.
        ! cache_ind carries the resulting particle -> record mapping.
        if( allocated(cache_ind) ) deallocate(cache_ind)
        allocate(cache_ind(nptcls), source=0)
        nsel = 0
        do iptcl = 1, nptcls
            if( build%spproj_field%get_state(iptcl) > 0 )then
                nsel             = nsel + 1
                cache_ind(iptcl) = nsel
            endif
        end do
        if( nsel < 1 )then
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: no active particles, cache disabled'
            deallocate(cache_ind)
            call stkname%kill
            call keyfile%kill
            call tmpstk%kill
            call tmpidx%kill
            return
        endif
        write(logfhandle,'(A,I8,A,I8,A,I4,A,I4)') '>>> BUILDING PARTICLE CACHE: cached=', nsel, &
            &'/', nptcls, ' box=', params%box, ' -> box_crop=', params%box_crop
        batchsz = min(nsel, max(1, params%nthr) * BATCHTHRSZ)
        allocate(pinds(nsel))
        do iptcl = 1, nptcls
            if( cache_ind(iptcl) > 0 ) pinds(cache_ind(iptcl)) = iptcl
        end do
        nbatches = ceiling(real(nsel) / real(batchsz))
        call prepimgbatch(params, build, batchsz)
        call alloc_imgarr(batchsz, [params%box_crop, params%box_crop, 1], params%smpd_crop, cache_imgs)
        call stkio_w%open(tmpstk, params%smpd_crop, 'write', box=params%box_crop, is_ft=.false.)
        do ibatch = 1, nbatches
            batch_start = (ibatch - 1) * batchsz + 1
            batch_end   = min(batch_start + batchsz - 1, nsel)
            nbatch      = batch_end - batch_start + 1
            call discrete_read_imgbatch(params, build, nsel, pinds, [batch_start, batch_end])
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1, nbatch
                ! iteration-independent prefix of prepimg4align, with no shift, then
                ! back to real space; the read path re-applies the forward transform
                call build%imgbatch(i)%norm_noise_fft_clip_shift(build%lmsk, cache_imgs(i), [0.,0.])
                call cache_imgs(i)%ifft
            end do
            !$omp end parallel do
            do i = 1, nbatch
                call stkio_w%write(batch_start + i - 1, cache_imgs(i))
            end do
            call progress(ibatch, nbatches)
        end do
        call stkio_w%close
        call dealloc_imgarr(cache_imgs)
        call killimgbatch(build)
        deallocate(pinds)
        call write_cache_index(tmpidx, nptcls)
        ! Validate what actually landed on disk before publishing it. A short write, a
        ! full filesystem or an interrupted run must not end up wearing a valid key.
        call find_ldim_nptcls(tmpstk, ldim_check, nptcls_check)
        if( ldim_check(1) /= params%box_crop .or. ldim_check(2) /= params%box_crop .or. &
           &nptcls_check /= nsel )then
            write(logfhandle,*) 'ldim/nptcls written : ', ldim_check(1:2), nptcls_check
            write(logfhandle,*) 'ldim/nptcls expected: ', params%box_crop, params%box_crop, nsel
            call del_file(tmpstk)
            call del_file(tmpidx)
            THROW_HARD('particle cache failed validation; ptcl_cache_ensure')
        endif
        ! Publish: data first, key last, so a reader either sees no key or sees a key
        ! backed by a complete stack and index.
        call simple_rename(tmpstk, stkname)
        call simple_rename(tmpidx, cache_idxname(params))
        call write_cache_key(params, build)
        call ptcl_cache_reset
        write(logfhandle,'(A,A)') '>>> PARTICLE CACHE WRITTEN: ', stkname%to_char()
        call stkname%kill
        call keyfile%kill
        call tmpstk%kill
        call tmpidx%kill
    end subroutine ptcl_cache_ensure

    !>  Fill build%imgbatch (sized box_crop) from the cache. pinds arrive sorted from
    !!  the samplers, so each reader walks its slice of the file forward.
    subroutine ptcl_cache_read_batch( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(dstack_io), allocatable :: dstkios(:)
        type(string) :: stkname
        integer :: nbatch, nstreams, istream, chunk, from_i, to_i, i, ii, irec
        if( batchlims(1) < 1 .or. batchlims(2) > n .or. batchlims(1) > batchlims(2) )then
            write(logfhandle,*) 'batchlims: ', batchlims
            write(logfhandle,*) 'n        : ', n
            THROW_HARD('invalid batchlims; ptcl_cache_read_batch')
        endif
        if( .not. allocated(cache_ind) ) THROW_HARD('cache index not loaded; ptcl_cache_read_batch')
        ! ptcl_cache_in_use has already established that every active particle has a
        ! record, so a zero here means the caller sampled a deselected particle
        do i = batchlims(1), batchlims(2)
            if( cache_ind(pinds(i)) < 1 )then
                write(logfhandle,*) 'iptcl: ', pinds(i)
                THROW_HARD('particle absent from the cache; ptcl_cache_read_batch')
            endif
        end do
        nbatch   = batchlims(2) - batchlims(1) + 1
        stkname  = cache_stkname(params)
        nstreams = io_read_nstreams(stkname, nbatch)
        chunk    = ceiling(real(nbatch) / real(nstreams))
        nstreams = ceiling(real(nbatch) / real(chunk))
        allocate(dstkios(nstreams))
        do istream = 1, nstreams
            call dstkios(istream)%new(params%smpd_crop, params%box_crop)
            call dstkios(istream)%cache_stack_info(stkname, &
                &[params%box_crop, params%box_crop, 1], maxval(cache_ind))
            call dstkios(istream)%open(stkname)
        end do
        !$omp parallel do default(shared) private(istream,from_i,to_i,i,ii,irec) schedule(static) &
        !$omp& proc_bind(close) num_threads(nstreams) if(nstreams > 1)
        do istream = 1, nstreams
            from_i = (istream - 1) * chunk + 1
            to_i   = min(from_i + chunk - 1, nbatch)
            do ii = from_i, to_i
                i    = batchlims(1) + ii - 1
                irec = cache_ind(pinds(i))
                call dstkios(istream)%read(stkname, irec, build%imgbatch(ii))
            end do
        end do
        !$omp end parallel do
        do istream = 1, nstreams
            call dstkios(istream)%kill
        end do
        deallocate(dstkios)
        call stkname%kill
    end subroutine ptcl_cache_read_batch

end module simple_ptcl_cache
