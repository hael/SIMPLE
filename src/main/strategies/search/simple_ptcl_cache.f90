!@descr: on-disk cache of noise-normalized, Fourier-cropped particles for the 2D matcher workflows
! An entry is the iteration-independent prefix of prepimg4align (noise norm vs lmsk, FFT, clip to box_crop),
! stored in real space (exact round trip). Covers every active particle. Read by prob_tab2D and refine2D_exec;
! restoration then runs at box_crop with no second normalization (cavger_init_online(cropped_ptcls=.true.)).
! refine3D and solve3D reject cache=yes. Contract: doc/policies/2D/particle_cache_policy.md
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
! Owner = the process that built or adopted the cache; workers never own. Paths are resolved at
! ownership time so exit cleanup (cache_cleanup_glob) does not depend on params.
logical      :: l_owner = .false.
type(string) :: owned_files(5)
! cache record per global particle index (0 = not cached: only active particles are); monotone, so
! sorted pinds read the file forward. cache_nrecs = maxval(cache_ind), memoized at index load.
integer, allocatable :: cache_ind(:)
integer              :: cache_nrecs = 0

contains

    !>  Where the cache lives: cache_dir, else $SIMPLE_PTCL_CACHE_DIR, else the execution directory (empty result).
    !!  It is large (nsel * box_crop^2 * 4 bytes), so a project on a slow disk wants it somewhere faster.
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

    !>  Per-run token in the cache basename: a hash of the execution directory, so concurrent runs sharing cache_dir
    !!  stay apart. Every rank recomputes it (workers cd into the master's directory, simple_qsys_ctrl), and unlike a
    !!  PID it survives a restart: a rerun in the same directory adopts or rebuilds what a killed run left.
    function cache_run_token( ) result( tok )
        type(string) :: tok
        integer(kind=8) :: h1, h2
        h1 = 1_8
        h2 = 1_8
        if( allocated(CWD_GLOB) ) call fold_str(h1, h2, CWD_GLOB)
        tok = string(int2str(int(h1))//'-'//int2str(int(h2)))
    end function cache_run_token

    !>  Cache file name from projname (for listings), box_crop (no stale crop) and the run token; a same-name leftover
    !!  is validated against the key by ptcl_cache_ensure, then adopted or rebuilt. tmp=.true. gives the name used while
    !!  writing, tagged before the extension because fname2format infers the image format from the extension.
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
        cache_nrecs      = maxval(cache_ind)
        read_cache_index = .true.
    end function read_cache_index

    !>  Global particle count, the length of the particle -> record map: global iptcl indexing lets a worker holding
    !!  [fromp,top] address the master's file. params%nptcls is not guaranteed to mean the same in both roles.
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

    !>  Geometry line of the key, from memory, cheap for every rank on every probe. The execution directory is
    !!  recorded verbatim, so a run whose directory hashes to the same token never validates; all ranks run in
    !!  that directory, so the line compares equal within a run.
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

    !>  Source fingerprint: each stack (name, range, nptcls_stk, size, mtime) plus every particle's stkind/indstk,
    !!  so a stack replaced in place or a remapped project cannot validate a stale cache. One stat per
    !!  stack + O(nptcls) folds; only the rank that decides whether to rebuild computes it.
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

    !>  Compares the key file (geometry line, then source fingerprint) with this run. full=.true. (ptcl_cache_ensure,
    !!  the deciding rank) recomputes both; consumers trust the fingerprint, since recomputing it has every worker stat
    !!  every stack, and ptcl_cache_ensure validated the sources moments earlier and rebuilds when they change.
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

    !>  The stack fingerprint walks os_stk, so only particle oritypes are cacheable: cls3D "particles"
    !!  are class averages in os_out, where a replaced cavg stack would evade staleness detection.
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
            &.not. oritype_cacheable(params) )then
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

    !>  Refuses to run uncached when cache=yes. Cropped restoration preprocesses differently (taper and gridding grid
    !!  at box_crop), so ranks mixing cache and fallback would sum two paths into one class average; a node-local
    !!  cache_dir makes that easy, so availability is mandatory and uniform.
    subroutine ptcl_cache_assert_ready( params, build )
        class(parameters), intent(in)    :: params
        class(builder),    intent(inout) :: build
        type(string) :: stkname
        if( .not. params%l_cache ) return
        if( params%box_crop >= params%box )     return  ! nothing to cache, see ptcl_cache_ensure
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
        cache_nrecs = 0
        if( allocated(cache_ind) ) deallocate(cache_ind)
    end subroutine ptcl_cache_reset

    !>  Takes ownership of the cache files (ptcl_cache_ensure only, before the first byte is written): resolves their
    !!  names and arms the exit cleanup hook, so an exception mid-build sweeps the temporaries too. When a staged
    !!  workflow changes box_crop, the previous stage's owned files are deleted at the handoff.
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

    !>  Deletes the owned files (none unless owner) and resets the probe/index state, which the next stage of a staged
    !!  run would otherwise read stale. Runs via cache_cleanup_glob on exit (simple_exec, single_exec) and on hard
    !!  exceptions (simple_error); the hook is disarmed first so a throw during deletion cannot re-enter.
    subroutine ptcl_cache_cleanup
        nullify(cache_cleanup_glob)
        call ptcl_cache_reset
        call delete_owned_files
    end subroutine ptcl_cache_cleanup

    !>  Uniform fallback to uncached execution, decided master-side before workers are scheduled. cache=no goes on the
    !!  cline too, as worker command lines are generated from it (ptcl_cache_assert_ready: why modes must not mix).
    subroutine disable_cache( params, cline )
        class(parameters), intent(inout) :: params
        class(cmdline),    intent(inout) :: cline
        params%l_cache = .false.
        params%cache   = 'no'
        call cline%set('cache', 'no')
    end subroutine disable_cache

    !>  Master-side, before workers: build if missing/stale or adopt a valid leftover (the builder/adopter
    !!  owns the files); any fallback sets cache=no on cline so all ranks agree.
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
            ! release any cache still owned from an earlier stage of a staged workflow
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
        ! Fast path for per-iteration prob_align2D calls: this process already owns this exact cache, so skip
        ! full revalidation (late activations fail in ptcl_cache_read_batch); the key stat catches deletion.
        if( l_owner )then
            keyfile = cache_keyname(params)
            if( (owned_files(1) .eq. keyfile) .and. file_exists(keyfile) )then
                call keyfile%kill
                return
            endif
            call keyfile%kill
        endif
        call ptcl_cache_reset
        ! the deciding rank pays for the full source fingerprint; consumers check only the geometry line
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
        ! Invalidate first: the key is the commit record. Freeing the stale stack/index now also lets the
        ! free-space budget below measure what the rebuild can use.
        call del_file(keyfile)
        call del_file(idxname)
        call del_file(stkname)
        call del_file(tmpstk)
        call del_file(tmpidx)
        nptcls = cache_nptcls(build)
        ! only active particles are cached (a deselected one is never sampled)
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
        ! At most 1/FREE_SPACE_DENOM of the destination's free space; over budget the whole run falls back
        ! uniformly. An unknown free space (fs_avail_bytes < 0) is no verdict, not zero.
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
        ! Full-box noise mask. Every current caller builds the general toolbox (build%lmsk); the
        ! disc fallback (same recipe as build_general_tbox) is defensive only.
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
        if( cache_nrecs < 1 )            THROW_HARD('cache record count not set; ptcl_cache_read_batch')
        ! ptcl_cache_in_use has already established that every active particle has a
        ! record, so a zero here means the caller sampled a deselected particle
        do i = batchlims(1), batchlims(2)
            if( cache_ind(pinds(i)) < 1 )then
                write(logfhandle,*) 'iptcl: ', pinds(i)
                THROW_HARD('particle absent from the cache; ptcl_cache_read_batch')
            endif
        end do
        stkname = cache_stkname(params)
        ! One handle, read serially: Fortran connects a file to at most one unit, and with small records
        ! and sorted pinds this is a forward scan of a file far smaller than the originals.
        call dstkio%new(params%smpd_crop, params%box_crop)
        call dstkio%cache_stack_info(stkname, &
            &[params%box_crop, params%box_crop, 1], cache_nrecs)
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
