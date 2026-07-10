!@descr: helper routines for sieve sub-project partitioning and chunk-map generation
module simple_ptcl_sieve_utils
use simple_stream_api
use simple_string,          only: string
use simple_cmdline,         only: cmdline
use simple_sp_project,      only: sp_project
use simple_parameters,      only: parameters
use simple_string_utils,    only: int2str_pad
use simple_core_module_api, only: simple_mkdir, simple_abspath

public :: generate_sieve_projects
private
#include "simple_local_flags.inc"

contains

    ! Build sieve sub-projects from a global project in three passes:
    !   pass 1 - count selected particles per active stack and determine chunk count
    !   pass 2 - record [first_stk, last_stk] stack boundaries for each chunk
    !   pass 3 - write one self-contained sub-project per chunk and register
    !            each active micrograph in project_list for ptcl_sieve chunk tracking
    ! Only stacks with state > 0 and nptcls > 0 are considered; inactive stacks
    ! are skipped in all three passes.
    subroutine generate_sieve_projects( project_list, params, spproj_glob, ntot_chunks )
        type(rec_list),       intent(inout) :: project_list
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(in)    :: spproj_glob
        integer,              intent(out)   :: ntot_chunks
        character(len=*),     parameter     :: DIR_PROJS = trim(PATH_HERE)//'spprojs/'
        type(sp_project)     :: spproj
        type(project_rec)    :: prec
        integer, allocatable :: stk_nptcls(:), stk_all_nptcls(:), chunks_map(:,:)
        type(string)         :: fname, absfname, path, projname, projfile
        integer              :: cnt, cnt_stk, ichunk, istk, iptcl, jptcl, kptcl
        integer              :: fromp, top, cfromp, ctop, n, nstks_chunk, nptcls_chunk, nptcls
        integer              :: tail_start, tail_end, nptcls_tot, nstks, n_chunks_available, nptcls_chunks_total
        logical              :: has_mic, l_truncate
        call simple_mkdir(DIR_PROJS)
        nstks  = spproj_glob%os_stk%get_noris()
        nptcls = spproj_glob%os_ptcl2D%get_noris()
        if( nstks == 0  ) THROW_HARD('No stacks found in project file: '//params%projfile%to_char())
        if( nptcls == 0 ) THROW_HARD('No particles found in project file: '//params%projfile%to_char())
        ! --- pass 1: count selected particles per active stack; tally chunks ---
        allocate(stk_all_nptcls(nstks), stk_nptcls(nstks), source=0)
        ntot_chunks = 0
        cnt         = 0
        cnt_stk     = 0
        do istk = 1, nstks
            if( spproj_glob%os_stk%get_state(istk) == 0       ) cycle
            if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0) cycle
            fromp = spproj_glob%os_stk%get_fromp(istk)
            top   = spproj_glob%os_stk%get_top(istk)
            stk_all_nptcls(istk) = top - fromp + 1          ! all particles (any state)
            do iptcl = fromp, top                           ! selected particles (state > 0)
                if( spproj_glob%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                stk_nptcls(istk) = stk_nptcls(istk) + 1
            enddo
            cnt     = cnt     + stk_nptcls(istk)
            cnt_stk = cnt_stk + 1
            if( (cnt_stk >= params%nmics) .or. (cnt >= params%nptcls_coarse) )then
                ntot_chunks = ntot_chunks + 1
                cnt         = 0
                cnt_stk     = 0
            endif
        enddo
        ! Tail policy: if leftover selected particles are large enough,
        ! create a new chunk; otherwise merge tail into the last chunk.
        if( cnt > 0 ) then
            if( ntot_chunks == 0 ) then
                ntot_chunks = 1
            else if( real(cnt) > 0.5 * real(params%nptcls_coarse) ) then
                ntot_chunks = ntot_chunks + 1
            endif
        endif
        nptcls_tot = sum(stk_nptcls)
        write(logfhandle,'(A,I8)') '>>> # OF STACKS               : ', nstks
        write(logfhandle,'(A,I8)') '>>> # OF SELECTED PARTICLES   : ', nptcls_tot
        write(logfhandle,'(A,I8)') '>>> # OF AVAILABLE CHUNKS     : ', ntot_chunks
        n_chunks_available = ntot_chunks
        l_truncate = .false.
        if( params%maxnchunks > 0 ) ntot_chunks = min(params%maxnchunks, ntot_chunks)
        if( params%maxnchunks > 0 ) l_truncate = ntot_chunks < n_chunks_available
        if( ntot_chunks == 0 ) return
        if( l_truncate ) &
            write(logfhandle,'(A,I8)') '>>> # OF REQUESTED CHUNKS     : ', ntot_chunks

        ! --- pass 2: record [first_stk, last_stk] boundaries for each chunk ---
        ! Tail handling follows pass-1 policy: large tail becomes a new chunk,
        ! otherwise it is merged into the last complete chunk.
        allocate(chunks_map(ntot_chunks,2), source=0)
        cnt     = 0
        cnt_stk = 0
        ichunk  = 0
        tail_start = 0
        tail_end   = 0
        do istk = 1, nstks
            if( spproj_glob%os_stk%get_state(istk) == 0       ) cycle
            if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0) cycle
            if( cnt == 0 ) tail_start = istk
            tail_end = istk
            cnt     = cnt     + stk_nptcls(istk)
            cnt_stk = cnt_stk + 1
            if( (cnt_stk >= params%nmics) .or. (cnt >= params%nptcls_coarse) )then
                if( ichunk >= ntot_chunks ) then
                    if( l_truncate ) then
                        exit
                    else
                        THROW_HARD('generate_sieve_projects: chunk mapping overflow (internal chunk count mismatch)')
                    end if
                end if
                ichunk = ichunk + 1
                chunks_map(ichunk,1) = tail_start            ! first active stack in chunk
                chunks_map(ichunk,2) = istk                  ! last active stack in chunk
                cnt     = 0
                cnt_stk = 0
                tail_start = 0
                if( l_truncate .and. ichunk == ntot_chunks ) exit
            endif
        enddo
        if( cnt > 0 .and. .not. l_truncate ) then
            if( ichunk == 0 ) then
                chunks_map(1,1) = tail_start
                chunks_map(1,2) = tail_end
            else if( real(cnt) > 0.5 * real(params%nptcls_coarse) .and. ichunk < ntot_chunks ) then
                chunks_map(ichunk+1,1) = tail_start
                chunks_map(ichunk+1,2) = tail_end
            else
                chunks_map(ichunk,2) = tail_end
            endif
        endif

        ! --- pass 3: build and write one sub-project per chunk ---
        has_mic        = spproj_glob%os_mic%get_noris() > 0
        spproj%compenv = spproj_glob%compenv
        spproj%jobproc = spproj_glob%jobproc
        call spproj%projinfo%new(1, is_ptcl=.false.)
        path = CWD_GLOB//'/'//DIR_PROJS
        call spproj%projinfo%set(1,'cwd', path)
        do ichunk = 1, ntot_chunks
            projname = int2str_pad(ichunk,6)
            projfile = path//projname//METADATA_EXT
            call spproj%projinfo%set(1,'projname', projname)
            call spproj%projinfo%set(1,'projfile', projfile)
            ! Size orientation containers for this chunk's active stacks only.
            nstks_chunk = 0
            do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                if( spproj_glob%os_stk%get_state(istk) == 0        ) cycle
                if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0 ) cycle
                nstks_chunk = nstks_chunk + 1
            end do
            nptcls = sum(stk_all_nptcls(chunks_map(ichunk,1):chunks_map(ichunk,2)))
            call spproj%os_stk%new(nstks_chunk, is_ptcl=.false.)
            call spproj%os_mic%new(nstks_chunk, is_ptcl=.false.)
            call spproj%os_ptcl2D%new(nptcls,   is_ptcl=.true.)
            cnt   = 0
            ctop  = 0
            jptcl = 0
            do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                if( spproj_glob%os_stk%get_state(istk) == 0        ) cycle
                if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0 ) cycle
                cnt = cnt + 1
                n   = stk_all_nptcls(istk)
                ! Transfer micrograph metadata (or synthesize from stack if absent).
                if( has_mic )then
                    call spproj%os_mic%transfer_ori(cnt, spproj_glob%os_mic, istk)
                else
                    call spproj%os_mic%set(cnt,'nptcls', n)
                    call spproj%os_mic%set_state(cnt, spproj_glob%os_stk%get_state(istk))
                endif
                ! Transfer stack metadata; convert stack path to absolute.
                call spproj%os_stk%transfer_ori(cnt, spproj_glob%os_stk, istk)
                call spproj%os_stk%getter(cnt,'stk', fname)
                absfname = simple_abspath(fname)
                call spproj%os_stk%set(cnt,'stk', absfname)
                ! Transfer particle orientations, re-indexing into the local project.
                fromp = spproj_glob%os_stk%get_fromp(istk)
                top   = spproj_glob%os_stk%get_top(istk)
                !$omp parallel do private(iptcl,kptcl) default(shared)
                do iptcl = fromp, top
                    kptcl = jptcl + iptcl - fromp + 1
                    call spproj%os_ptcl2D%transfer_ori(kptcl, spproj_glob%os_ptcl2D, iptcl)
                    call spproj%os_ptcl2D%set_stkind(kptcl, cnt)
                enddo
                !$omp end parallel do
                ! Update local fromp/top in the sub-project stack record.
                jptcl  = jptcl  + n
                cfromp = ctop   + 1
                ctop   = cfromp + n - 1
                call spproj%os_stk%set(cnt,'fromp', cfromp)
                call spproj%os_stk%set(cnt,'top',   ctop)
                ! Register micrograph in project_list for ptcl_sieve chunk tracking.
                prec%id         = project_list%size() + 1
                prec%projname   = projfile
                prec%micind     = cnt
                prec%nptcls     = n
                prec%nptcls_sel = stk_nptcls(istk)
                prec%included   = .false.
                call project_list%push_back(prec)
            enddo
            ! Strip any prior 2D results before writing the sub-project.
            call spproj%os_ptcl2D%delete_2Dclustering(keepshifts=.false., keepcls=.false.)
            call spproj%write(projfile)
            ! Report selected-particle count for this chunk.
            nptcls_chunk = 0
            do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                nptcls_chunk = nptcls_chunk + stk_nptcls(istk)
            end do
        enddo
        nptcls_chunks_total = 0
        do ichunk = 1, ntot_chunks
            do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                nptcls_chunks_total = nptcls_chunks_total + stk_nptcls(istk)
            end do
        end do
        write(logfhandle,'(A,I8,A,I8,A,I8,A)') '>>> ', ntot_chunks,' CHUNKS CREATED CONTAINING ', nptcls_chunks_total, ' / ', nptcls_tot, ' PARTICLES'
        call spproj%kill
        deallocate(stk_all_nptcls, stk_nptcls, chunks_map)
    end subroutine generate_sieve_projects

end module simple_ptcl_sieve_utils