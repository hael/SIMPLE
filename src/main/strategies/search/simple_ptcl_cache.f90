!@descr: downscaled particle cache shared by the 2D and 3D matcher workflows
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
! whole 2D run and there is nothing to invalidate between stages. abinitio3D
! changes box_crop per stage, so its refine3D children rebuild once per stage
! boundary and the ownership handoff in ptcl_cache_own deletes the previous
! stage's files. The cache lives exactly as long as the run that owns it: the
! process that builds (or adopts) it removes the files on normal exit and on
! hard exception, via the cache_cleanup_glob hook in simple_defs. Only
! SIGKILL-class deaths can leave files behind, and a rerun in the same execution
! directory finds and reclaims those (see cache_run_token).
!
! 2D consumers: prob_tab2D, which needs only the alignment product, and
! cluster2D_exec, which additionally restores class averages from the raw
! images. The latter is not a pure I/O substitution -- the edge taper and the
! gridding source grid then live at box_crop rather than at box -- and it requires
! cavger_init_online(cropped_ptcls=.true.), which switches the class-averager work
! images to ldim_croppd, the mask to lmsk_crop, and scales the shift and the CTF
! pixel size to the cropped grid. Note in particular that a cached particle is
! already noise-normalized, so restoration must not normalize it a second time:
! cropping discards most of the noise power, and re-normalizing would rescale the
! particle by the crop factor while leaving ctfsq_plane and sigma2 untouched.
!
! 3D consumers (refine3D search, prob_tab, prob_tab_neigh) take alignment reads
! only: the cached entry is exactly the iteration-independent prefix of their
! prepimg4align as well. The 3D reconstruction phase is deliberately excluded --
! its preprocessing (norm_noise_taper_edge_pad_fft in simple_matcher_3Drec)
! tapers and pads at full box, so it keeps re-reading the originals. A denoised
! primary source (ptcl_src=den) disqualifies the cache outright; a denoised
! objective (objfun_den=yes) does not, since those images are read separately.
module simple_ptcl_cache
use simple_pftc_srch_api
use simple_builder,           only: builder
use simple_cmdline,           only: cmdline
use simple_discrete_stack_io, only: dstack_io
use simple_stack_io,          only: stack_io
use simple_imgarr_utils,      only: alloc_imgarr, dealloc_imgarr
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, killimgbatch
use simple_syslib,            only: simple_mkdir, dir_exists, simple_file_stat, simple_rename, fs_avail_bytes
implicit none
#include "simple_local_flags.inc"

public :: ptcl_cache_in_use, ptcl_cache_assert_ready, ptcl_cache_ensure
public :: ptcl_cache_read_batch, ptcl_cache_reset, ptcl_cache_cleanup
private

character(len=*), parameter :: CACHE_FBODY   = 'ptcl_cache_'
character(len=*), parameter :: CACHE_DIR_ENV = 'SIMPLE_PTCL_CACHE_DIR'
character(len=*), parameter :: TMP_TAG       = '_part'
! the cache may not claim more than this fraction of the free space at cache_dir
integer(kind=8),  parameter :: FREE_SPACE_DENOM = 4_8
! rolling-hash constants shared by fold_str/fold_int
integer(kind=8),  parameter :: HASH_M1 = 2147483647_8, HASH_A1 = 16807_8  ! 2^31-1
integer(kind=8),  parameter :: HASH_M2 = 2147483629_8, HASH_A2 = 48271_8  ! largest prime < 2^31-1

! resolved once per process by ptcl_cache_in_use
logical :: l_probed    = .false.
logical :: l_available = .false.
! Exit-time cleanup: the process that builds -- or adopts, see ptcl_cache_ensure --
! the cache owns its files and removes them on normal exit or hard exception, via
! the cache_cleanup_glob hook that simple_exception and the exec programs call.
! Workers never take ownership, so a dying worker cannot delete the cache out from
! under the other ranks or a resubmitted part. Paths are resolved at ownership
! time so that cleanup does not depend on params still being alive.
logical      :: l_owner = .false.
type(string) :: owned_files(5)
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

    !>  Per-run discriminator folded into the cache basename: a hash of the absolute
    !!  execution directory. cache_dir is explicitly meant to be shared, so several
    !!  concurrent runs pointing at the same fast disk must not delete or overwrite
    !!  each other, and projname alone cannot guarantee that. The execution directory
    !!  can: no two live runs share one, every rank can recompute it locally (the
    !!  distributed workers cd into the master's directory before they start, see
    !!  simple_qsys_ctrl), and unlike a PID it survives a restart, so a rerun in the
    !!  same directory finds -- and can adopt or rebuild -- what a killed predecessor
    !!  left behind instead of orphaning it.
    function cache_run_token( ) result( tok )
        type(string) :: tok
        integer(kind=8) :: h1, h2
        h1 = 1_8
        h2 = 1_8
        if( allocated(CWD_GLOB) ) call fold_str(h1, h2, CWD_GLOB)
        tok = string(int2str(int(h1))//'-'//int2str(int(h2)))
    end function cache_run_token

    !>  Cache files for this run. The run token keeps concurrent runs writing to a
    !!  shared cache_dir apart; the project name is kept in the basename for the
    !!  human reading a directory listing, and box_crop so a run with a different
    !!  crop cannot pick up a stale cache. A leftover with the same name -- same
    !!  directory recreated after deletion, or a killed run -- is not an error:
    !!  ptcl_cache_ensure validates it against the key (which records the execution
    !!  directory verbatim, so even a token collision cannot validate) and either
    !!  adopts or rebuilds in place.
    !!
    !!  tmp=.true. returns the name used while the cache is being written. The tag goes
    !!  before the extension, never after: image format is inferred from the extension
    !!  (fname2format), so a name ending in the tag would be rejected as an unsupported
    !!  format the moment anything tried to read its header.
    function cache_fname( params, ext, tmp ) result( fname )
        class(parameters), intent(in) :: params
        character(len=*),  intent(in) :: ext
        logical, optional, intent(in) :: tmp
        type(string) :: fname, dirname, basename_here, token
        character(len=:), allocatable :: tag
        allocate(tag, source='')
        if( present(tmp) )then
            if( tmp ) deallocate(tag)
            if( tmp ) allocate(tag, source=TMP_TAG)
        endif
        token = cache_run_token()
        if( params%projname%is_blank() )then
            basename_here = string(CACHE_FBODY//int2str(params%box_crop)//&
                &'_'//token%to_char()//tag//ext)
        else
            basename_here = string(CACHE_FBODY//params%projname%to_char()//&
                &'_'//int2str(params%box_crop)//'_'//token%to_char()//tag//ext)
        endif
        call token%kill
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
        integer :: i
        do i = 1, len_trim(str)
            h1 = mod(h1 * HASH_A1 + int(iachar(str(i:i)), 8), HASH_M1)
            h2 = mod(h2 * HASH_A2 + int(iachar(str(i:i)), 8), HASH_M2)
        end do
        ! separator, so that concatenating fields cannot alias a different split
        h1 = mod(h1 * HASH_A1 + 255_8, HASH_M1)
        h2 = mod(h2 * HASH_A2 + 255_8, HASH_M2)
    end subroutine fold_str

    !>  Integer variant of fold_str, without the cost of an intermediate string; used
    !!  where millions of values are folded (the particle -> stack mapping).
    pure subroutine fold_int( h1, h2, val )
        integer(kind=8), intent(inout) :: h1, h2
        integer,         intent(in)    :: val
        h1 = mod(h1 * HASH_A1 + int(val, 8), HASH_M1)
        h2 = mod(h2 * HASH_A2 + int(val, 8), HASH_M2)
    end subroutine fold_int

    !>  Geometry the cached pixels depend on, all of it already in memory. Cheap enough
    !!  for every rank to recompute on every probe. The execution directory is recorded
    !!  verbatim so that a cache left by a different run whose directory happens to hash
    !!  to the same token can never validate; all ranks share the directory (workers cd
    !!  into it), so the line compares equal within a run.
    function cache_key_geom( params, build ) result( key )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: key
        character(len=:), allocatable :: dir_here
        if( allocated(CWD_GLOB) )then
            allocate(dir_here, source=CWD_GLOB)
        else
            allocate(dir_here, source='')
        endif
        key = string('dir='//dir_here//&
            &' box='//int2str(params%box)//&
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
    !!  a stack replaced in place or a changed stack in the middle, either of which would
    !!  serve stale pixels. Size and mtime make the check conservative -- a touched file
    !!  forces a rebuild -- which is the right way to be wrong.
    !!
    !!  The particle -> image mapping is folded as well: the image a particle reads is
    !!  determined by its stkind and its explicit indstk when present (see
    !!  map_ptcl_ind2stk_ind), with the range-derived fallback covered by the per-stack
    !!  fromp/top already in the hash. Without this, a project reordered or remapped
    !!  over unchanged stacks -- e.g. adopting a SIGKILL orphan after such an edit --
    !!  would silently associate particles with the wrong images.
    !!
    !!  This costs one stat per stack plus O(nptcls) integer folds, so only the rank
    !!  that decides whether to rebuild pays it, once per build decision. See
    !!  cache_key_matches and the ownership fast path in ptcl_cache_ensure.
    function cache_key_stkfp( build ) result( key )
        class(builder), intent(inout) :: build
        type(string) :: key, stkname
        integer, allocatable :: statbuf(:)
        integer(kind=8) :: h1, h2
        integer :: nstks, istk, stat_status, nptcls, iptcl, stkind, indstk
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
        nptcls = cache_nptcls(build)
        do iptcl = 1, nptcls
            stkind = 0
            indstk = 0
            if( build%spproj_field%isthere(iptcl, 'stkind') ) stkind = build%spproj_field%get_int(iptcl, 'stkind')
            if( build%spproj_field%isthere(iptcl, 'indstk') ) indstk = build%spproj_field%get_int(iptcl, 'indstk')
            call fold_int(h1, h2, stkind)
            call fold_int(h1, h2, indstk)
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

    !>  Entries are derived from the raw particle stacks, so a workflow whose primary
    !!  particle source is the denoised stack (ptcl_src=den) must not be served from
    !!  the cache: the pixels would simply be the wrong ones. A denoised *objective*
    !!  (objfun_den=yes) is different: its images are read separately after the raw
    !!  primary batch (see build_batch_particles3D), so the raw alignment read may
    !!  still come from the cache. Checked wherever cache state is decided, so even
    !!  a mis-plumbed worker command line cannot read wrong data.
    logical function primary_src_is_den( params )
        class(parameters), intent(in) :: params
        primary_src_is_den = params%l_ptcl_src_den
    end function primary_src_is_den

    !>  The stack fingerprint walks os_stk, so only workflows whose images come from
    !!  the imported particle stacks can be cached. cls3D "particles" are class
    !!  averages living in os_out; a replaced cavg stack would evade staleness
    !!  detection, so those workflows are refused rather than served wrong pixels.
    logical function oritype_cacheable( params )
        class(parameters), intent(in) :: params
        oritype_cacheable = trim(params%oritype) == 'ptcl2D' .or. trim(params%oritype) == 'ptcl3D'
    end function oritype_cacheable

    !>  Is a usable cache present for this run? Probed once per process. A missing or
    !!  stale cache is not an error here; ptcl_cache_assert_ready is what makes it one
    !!  when the user asked for the cache.
    logical function ptcl_cache_in_use( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: stkname
        if( .not. params%l_cache .or. params%box_crop >= params%box .or. &
            &primary_src_is_den(params) .or. .not. oritype_cacheable(params) )then
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
        if( params%box_crop >= params%box )     return  ! nothing to cache, see ptcl_cache_ensure
        if( primary_src_is_den(params) )        return  ! refused uniformly, see ptcl_cache_ensure
        if( .not. oritype_cacheable(params) )   return  ! refused uniformly, see ptcl_cache_ensure
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

    !>  Take ownership of the on-disk cache files: resolve their names now and arm the
    !!  exit-time cleanup hook. Called by ptcl_cache_ensure only, i.e. by the process
    !!  that builds or adopts the cache, before the first byte is written, so that an
    !!  exception mid-build also sweeps the temporaries.
    !!
    !!  A staged workflow (abinitio3D) re-owns once per stage with a stage-specific
    !!  box_crop, so the names change under one process. The previous stage's cache
    !!  is dead weight the moment a new one is owned -- no later stage can use a
    !!  different crop -- so it is deleted at the handoff rather than left to
    !!  accumulate until exit.
    subroutine ptcl_cache_own( params )
        class(parameters), intent(in) :: params
        type(string) :: newkey
        newkey = cache_keyname(params)
        if( l_owner )then
            ! files only: ensure has live probe/index state at this point, so the
            ! full cleanup (which resets it) must not run here
            if( .not. (owned_files(1) .eq. newkey) ) call delete_owned_files
        endif
        call newkey%kill
        owned_files(1) = cache_keyname(params)
        owned_files(2) = cache_idxname(params)
        owned_files(3) = cache_stkname(params)
        owned_files(4) = cache_fname(params, '_idx'//BIN_EXT, tmp=.true.)
        owned_files(5) = cache_fname(params, STK_EXT,         tmp=.true.)
        l_owner = .true.
        cache_cleanup_glob => ptcl_cache_cleanup
    end subroutine ptcl_cache_own

    !>  Delete the owned files and drop ownership, leaving probe/index state alone.
    !!  Key file first, so a partially-completed deletion can never leave a cache
    !!  that still validates.
    subroutine delete_owned_files
        integer :: i
        if( .not. l_owner ) return
        l_owner = .false.
        do i = 1, size(owned_files)
            call del_file(owned_files(i))
            call owned_files(i)%kill
        end do
    end subroutine delete_owned_files

    !>  Delete the owned cache files; a no-op in every process that never took
    !!  ownership. Runs on normal exit (called from the exec programs) and on hard
    !!  exception (via cache_cleanup_glob in simple_exception). The hook is disarmed
    !!  first so a throw during deletion cannot re-enter, and the key file goes first
    !!  so a partially-completed cleanup can never leave a cache that still validates.
    !!
    !!  The memoized probe/index state goes with the files: in a staged in-memory run
    !!  the same process consumes the next stage's cache, and a stale l_available or
    !!  cache_ind from this one would have it read a deleted file.
    subroutine ptcl_cache_cleanup
        nullify(cache_cleanup_glob)
        call ptcl_cache_reset
        call delete_owned_files
    end subroutine ptcl_cache_cleanup

    !>  Uniform fallback to uncached execution, decided master-side before any worker
    !!  is scheduled. Flipping cache=no on the command line matters as much as on
    !!  params: worker command lines are generated from the cline, and cache
    !!  availability must be uniform across ranks -- ptcl_cache_assert_ready explains
    !!  why mixed modes are forbidden.
    subroutine disable_cache( params, cline )
        class(parameters), intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        params%l_cache = .false.
        params%cache   = 'no'
        call cline%set('cache', 'no')
    end subroutine disable_cache

    !>  Build the cache if it is missing or stale. Single sequential pass over the
    !!  particles in index order, which is also the only wholly sequential read of
    !!  the originals the pipeline ever performs. Master-side, call before workers;
    !!  cline is the command line the worker jobs will be generated from, so that a
    !!  fallback decision here reaches every rank.
    !!
    !!  This is also where cache ownership is taken: whichever process builds the
    !!  cache -- or adopts a valid one left in place, e.g. by a killed predecessor in
    !!  the same execution directory -- arms the exit-time cleanup that removes the
    !!  files on normal termination or hard exception.
    subroutine ptcl_cache_ensure( params, build, cline )
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        type(image), allocatable :: cache_imgs(:)
        type(image)     :: mskimg
        type(stack_io)  :: stkio_w
        type(string)    :: stkname, dirname, keyfile, idxname, tmpstk, tmpidx
        integer(kind=8) :: want_bytes, avail_bytes
        logical, allocatable :: lmsk(:,:,:)
        integer, allocatable :: pinds(:)
        integer :: batchsz, nbatches, ibatch, batch_start, batch_end, nbatch, i, iptcl, nptcls, nsel
        integer :: ldim_check(3), nptcls_check
        if( .not. params%l_cache ) return
        if( params%box_crop >= params%box )then
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: box_crop == box, nothing to gain, running without cache'
            call disable_cache(params, cline)
            ! in a staged workflow a previous stage's cache may still be owned; it is
            ! unusable from here on (crops only grow), so release it now
            call ptcl_cache_cleanup
            return
        endif
        if( primary_src_is_den(params) )then
            ! ptcl_src is staged, so a denoised stage may follow a cached one;
            ! release the previous stage's cache -- crops only grow, it is dead
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: denoised particle source in use, running without cache'
            call disable_cache(params, cline)
            call ptcl_cache_cleanup
            return
        endif
        if( .not. oritype_cacheable(params) )then
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: only particle stacks can be cached (oritype='//&
                &trim(params%oritype)//'), running without cache'
            call disable_cache(params, cline)
            call ptcl_cache_cleanup
            return
        endif
        ! Fast path for the per-iteration prob_align calls: this process already owns
        ! exactly this cache -- it built or adopted it earlier in the run, and per-run
        ! naming means no other process may touch it -- so the full revalidation
        ! (source stat sweep, index reload, active-particle scan) is skipped and the
        ! memoized probe state stays valid. A particle activated mid-run past this
        ! point fails loudly at read time (ptcl_cache_read_batch) rather than being
        ! served wrong data; one stat guards against out-of-band file deletion.
        if( l_owner )then
            keyfile = cache_keyname(params)
            if( (owned_files(1) .eq. keyfile) .and. file_exists(keyfile) )then
                call keyfile%kill
                return
            endif
            call keyfile%kill
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
                    ! adopt it: this run is now responsible for removing it on exit
                    call ptcl_cache_own(params)
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
        idxname = cache_idxname(params)
        tmpstk  = cache_fname(params, STK_EXT,          tmp=.true.)
        tmpidx  = cache_fname(params, '_idx'//BIN_EXT,  tmp=.true.)
        ! Invalidate first: the key is the commit record, so from here until it is
        ! rewritten no reader can accept whatever state the cache files are in. A
        ! stale stack and index go now rather than being overwritten by the publish
        ! renames: without the key they are unreachable anyway, and freeing them
        ! first lets the free-space budget below measure what the rebuild can
        ! actually use.
        call del_file(keyfile)
        call del_file(idxname)
        call del_file(stkname)
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
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: no active particles, running without cache'
            call disable_cache(params, cline)
            ! releases any previous stage's cache and deallocates cache_ind via reset
            call ptcl_cache_cleanup
            call stkname%kill
            call keyfile%kill
            call idxname%kill
            call tmpstk%kill
            call tmpidx%kill
            return
        endif
        ! The cache may claim at most a quarter of the free space at its destination,
        ! so a large dataset cannot starve whatever else lives there. The predicted
        ! size is exact: nsel records of box_crop^2 reals, the stack header, and the
        ! index file. Over budget, the whole run falls back to uncached execution --
        ! uniformly, via disable_cache, before any worker command line is generated.
        ! An unknown free space (fs_avail_bytes < 0, unsupported filesystem) is
        ! treated as no verdict rather than as zero.
        want_bytes = int(nsel,8) * int(params%box_crop,8)**2 * 4_8 &
            &+ 1024_8 + 4_8 * int(nptcls + 1, 8)
        dirname = cache_dir(params)
        if( dirname%is_blank() ) dirname = string('.')
        avail_bytes = fs_avail_bytes(dirname)
        call dirname%kill
        if( avail_bytes >= 0_8 .and. want_bytes > avail_bytes / FREE_SPACE_DENOM )then
            write(logfhandle,'(A,F12.2,A,F12.2,A)') '>>> PARTICLE CACHE: needs ', &
                &real(want_bytes)/real(1024**3), ' GB, more than 25% of the ', &
                &real(avail_bytes)/real(1024**3), ' GB free at its destination'
            write(logfhandle,'(A)') '>>> PARTICLE CACHE: running without cache; point cache_dir at more storage to enable it'
            call disable_cache(params, cline)
            ! releases any previous stage's cache and deallocates cache_ind via reset
            call ptcl_cache_cleanup
            call stkname%kill
            call keyfile%kill
            call idxname%kill
            call tmpstk%kill
            call tmpidx%kill
            return
        endif
        ! own the files from here on: an exception at any point during the build must
        ! sweep the temporaries as well as the published names
        call ptcl_cache_own(params)
        write(logfhandle,'(A,I8,A,I8,A,I4,A,I4)') '>>> BUILDING PARTICLE CACHE: cached=', nsel, &
            &'/', nptcls, ' box=', params%box, ' -> box_crop=', params%box_crop
        batchsz = min(nsel, max(1, params%nthr) * BATCHTHRSZ)
        allocate(pinds(nsel))
        do iptcl = 1, nptcls
            if( cache_ind(iptcl) > 0 ) pinds(cache_ind(iptcl)) = iptcl
        end do
        nbatches = ceiling(real(nsel) / real(batchsz))
        ! The full-box noise-normalization mask. The refine3D distributed master
        ! builds only the project, not the general toolbox, so build%lmsk need not
        ! exist in the process that builds the cache; construct the identical disc
        ! mask locally in that case (same recipe as build_general_tbox).
        if( allocated(build%lmsk) )then
            lmsk = build%lmsk
        else
            call mskimg%disc([params%box, params%box, 1], params%smpd, params%msk, lmsk)
            call mskimg%kill
        endif
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
                call build%imgbatch(i)%norm_noise_fft_clip_shift(lmsk, cache_imgs(i), [0.,0.])
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
        deallocate(pinds, lmsk)
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
        call simple_rename(tmpidx, idxname)
        call write_cache_key(params, build)
        call ptcl_cache_reset
        write(logfhandle,'(A,A)') '>>> PARTICLE CACHE WRITTEN: ', stkname%to_char()
        call stkname%kill
        call keyfile%kill
        call idxname%kill
        call tmpstk%kill
        call tmpidx%kill
    end subroutine ptcl_cache_ensure

    !>  Fill build%imgbatch (sized box_crop) from the cache. pinds arrive sorted from
    !!  the samplers, so each reader walks its slice of the file forward.
    subroutine ptcl_cache_read_batch( params, build, n, pinds, batchlims )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: n, pinds(n), batchlims(2)
        type(dstack_io) :: dstkio
        type(string) :: stkname
        integer :: i, ii, irec
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
        stkname = cache_stkname(params)
        ! One handle, read serially.
        !
        ! The whole cache is a single file, and Fortran allows a file to be connected to
        ! at most one unit at a time, so the several-readers-over-one-file arrangement
        ! used for particle stacks is not available here -- a second open of the same
        ! path fails. Nor would it buy anything: libgfortran serializes access per unit,
        ! so threads sharing one connection queue up regardless.
        !
        ! This is much less costly than it sounds. Records are box_crop-sized and pinds
        ! arrive sorted, so this is a forward scan over a file an order of magnitude
        ! smaller than the originals it replaces.
        call dstkio%new(params%smpd_crop, params%box_crop)
        call dstkio%cache_stack_info(stkname, &
            &[params%box_crop, params%box_crop, 1], maxval(cache_ind))
        call dstkio%open(stkname)
        do i = batchlims(1), batchlims(2)
            ii   = i - batchlims(1) + 1
            irec = cache_ind(pinds(i))
            call dstkio%read(stkname, irec, build%imgbatch(ii))
        end do
        call dstkio%kill
        call stkname%kill
    end subroutine ptcl_cache_read_batch

end module simple_ptcl_cache
