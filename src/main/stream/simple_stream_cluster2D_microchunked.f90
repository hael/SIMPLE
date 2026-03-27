!@descr: offline incremental 2D classification via microchunked processing
!
! This module drives post-collection 2D classification using the microchunked2D
! framework. The dataset is partitioned into fixed-size chunks, each independently
! classified, then progressively merged in pipeline order. The result is a
! fully-connected chain of pass-1 → pass-2 → reference → match classifications
! with no global reclustering step.
!
! Workflow
! --------
! (1) The global project is scanned; only active stacks (state > 0, nptcls > 0)
!     are considered. Chunks are bounded by nmics (micrographs) and maxnptcls
!     (selected particles). Each chunk is written as a self-contained sub-project
!     under spprojs/ and registered in project_list.
! (2) microchunked2D%cycle() drives the multi-tier processing loop. The loop
!     exits when get_finished() returns true (all tiers complete).
! (3) Completed match-chunk projects are combined into a single output project.
!
! Chunk tiers (handled by microchunked2D)
! ----------------------------------------
!   pass-1  : per-chunk ab-initio 2D    (up to MICROCHUNK_P1_THRESHOLD particles)
!   pass-2  : merged pass-1 results     (up to MICROCHUNK_P2_THRESHOLD particles)
!   refchunk: merged pass-2 results     — establishes the reference class averages
!   match   : remaining chunks matched  against the refchunk class averages
!
module simple_stream_cluster2D_microchunked
use simple_stream_api
use simple_microchunked2D, only: microchunked2D
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_cluster2D_microchunked
  contains
    procedure :: execute => exec_stream_cluster2D_microchunked
end type stream_cluster2D_microchunked

contains

    subroutine exec_stream_cluster2D_microchunked( self, cline )
        class(stream_cluster2D_microchunked), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=*), parameter :: DIR_PROJS = trim(PATH_HERE)//'spprojs/'
        integer,          parameter :: WAITTIME  = 5
        type(rec_list)       :: project_list
        type(parameters)     :: params
        type(sp_project)     :: spproj_glob
        type(microchunked2D) :: chunked_2D
        integer              :: nstks, nptcls, nptcls_tot, ntot_chunks
        ! --- default command-line parameters ---
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',    'yes')
        if( .not. cline%defined('nmics')     ) call cline%set('nmics',      5)
        if( .not. cline%defined('maxnptcls') ) call cline%set('maxnptcls', 1000)
        ! --- initialise parameters and read project ---
        call params%new(cline)
        call spproj_glob%read_non_data_segments(params%projfile)
        call spproj_glob%read_segment('mic',    params%projfile)
        call spproj_glob%read_segment('stk',    params%projfile)
        call spproj_glob%read_segment('ptcl2D', params%projfile)
        ! --- sanity checks ---
        nstks  = spproj_glob%os_stk%get_noris()
        nptcls = spproj_glob%get_nptcls()
        if( spproj_glob%os_mic%get_noris() > 0 )then
            if( spproj_glob%os_mic%get_noris() /= nstks )&
                THROW_HARD('Inconsistent # of micrographs and stacks, use prune_project.')
        endif
        if( nptcls == 0 )&
            THROW_HARD('No particles found in project file: '//params%projfile%to_char()//'; exec_cluster2d_subsets')
        ! --- partition dataset into per-chunk sub-projects ---
        call generate_chunk_projects()
        if( ntot_chunks == 0 )then
            call simple_end('**** NO CHUNKS GENERATED — INSUFFICIENT PARTICLES OR MICROGRAPHS ****')
        endif
        ! --- drive multi-tier classification loop until all tiers are complete ---
        call simple_mkdir(PATH_HERE // DIR_STREAM_COMPLETED)
        call chunked_2D%new(params, string(PATH_HERE // DIR_STREAM_COMPLETED))
        do
            call chunked_2D%cycle(project_list)
            if( chunked_2D%get_finished() ) exit
            call sleep(WAITTIME)
        end do
        ! --- combine match results and finalise ---
            ! Build combined projfile path and skip if already written

        call chunked_2D%combine_completed_match_chunks(string('microchunks_match_combined' // METADATA_EXT))
        call chunked_2D%kill()
        call spproj_glob%kill()
        call simple_rmdir(STDERROUT_DIR)
        call simple_rmdir(DIR_PROJS)
        call simple_rmdir(DIR_SNAPSHOT)
        call del_file(POOL_DIR//POOL_PROJFILE)
        call simple_rmdir(SIGMAS_DIR)
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_CLUSTER2D_SUBSETS NORMAL STOP ****')

    contains

        ! Scans the global project in three passes:
        !   pass 1 — counts active particles per stack; determines ntot_chunks
        !   pass 2 — records [first_stk, last_stk] boundaries for each chunk
        !   pass 3 — writes one self-contained sub-project per chunk and registers
        !             each micrograph in project_list for microchunked2D tracking
        ! Only stacks with state > 0 and nptcls > 0 are considered; zero-state
        ! stacks are skipped in all three passes.
        subroutine generate_chunk_projects
            type(sp_project)     :: spproj
            type(project_rec)    :: prec
            integer, allocatable :: stk_nptcls(:), stk_all_nptcls(:), chunks_map(:,:)
            type(rec_list)       :: project_list_slice
            type(string)         :: fname, absfname, path, projname, projfile
            integer              :: cnt, cnt_stk, ichunk, istk, iptcl, jptcl, kptcl
            integer              :: fromp, top, cfromp, ctop, n, nstks_chunk, nptcls_chunk
            logical              :: has_mic
            call simple_mkdir(DIR_PROJS)

            ! --- pass 1: count active particles per stack; tally chunks ---
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
                do iptcl = fromp, top                            ! selected particles (state > 0)
                    if( spproj_glob%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    stk_nptcls(istk) = stk_nptcls(istk) + 1
                enddo
                cnt     = cnt     + stk_nptcls(istk)
                cnt_stk = cnt_stk + 1
                if( (cnt_stk >= params%nmics) .or. (cnt >= params%maxnptcls) )then
                    ntot_chunks = ntot_chunks + 1
                    cnt         = 0
                    cnt_stk     = 0
                endif
            enddo
            nptcls_tot = sum(stk_nptcls)
            write(logfhandle,'(A,I8)') '>>> # OF STACKS          : ', nstks
            write(logfhandle,'(A,I8)') '>>> # OF PARTICLES       : ', nptcls_tot
            write(logfhandle,'(A,I8)') '>>> # OF AVAILABLE CHUNKS: ', ntot_chunks
            if( cline%defined('maxnchunks') ) ntot_chunks = min(params%maxnchunks, ntot_chunks)
            if( ntot_chunks == 0 ) return

            ! --- pass 2: record [first_stk, last_stk] boundaries for each chunk ---
            ! The tail (stacks beyond the last complete chunk) is intentionally abandoned.
            allocate(chunks_map(ntot_chunks,2), source=0)
            cnt     = 0
            cnt_stk = 0
            ichunk  = 1
            do istk = 1, nstks
                if( ichunk > ntot_chunks ) exit
                if( spproj_glob%os_stk%get_state(istk) == 0       ) cycle
                if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0) cycle
                if( cnt == 0 ) chunks_map(ichunk,1) = istk      ! first active stack in chunk
                cnt     = cnt     + stk_nptcls(istk)
                cnt_stk = cnt_stk + 1
                if( (cnt_stk >= params%nmics) .or. (cnt >= params%maxnptcls) )then
                    chunks_map(ichunk,2) = istk                  ! last active stack in chunk
                    ichunk  = ichunk + 1
                    cnt     = 0
                    cnt_stk = 0
                endif
            enddo

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
                ! size the oris objects for this chunk's stack range
                nstks_chunk = chunks_map(ichunk,2) - chunks_map(ichunk,1) + 1
                nptcls      = sum(stk_all_nptcls(chunks_map(ichunk,1):chunks_map(ichunk,2)))
                call spproj%os_stk%new(nstks_chunk, is_ptcl=.false.)
                call spproj%os_mic%new(nstks_chunk, is_ptcl=.false.)
                call spproj%os_ptcl2D%new(nptcls,   is_ptcl=.true.)
                cnt   = 0
                ctop  = 0
                jptcl = 0
                do istk = chunks_map(ichunk,1), chunks_map(ichunk,2)
                    cnt = cnt + 1
                    n   = stk_all_nptcls(istk)
                    ! transfer micrograph metadata (or synthesise from stack if absent)
                    if( has_mic )then
                        call spproj%os_mic%transfer_ori(cnt, spproj_glob%os_mic, istk)
                    else
                        call spproj%os_mic%set(cnt,'nptcls', n)
                        call spproj%os_mic%set_state(cnt, spproj_glob%os_stk%get_state(istk))
                    endif
                    ! transfer stack metadata; convert stack path to absolute
                    call spproj%os_stk%transfer_ori(cnt, spproj_glob%os_stk, istk)
                    call spproj%os_stk%getter(cnt,'stk', fname)
                    absfname = simple_abspath(fname)
                    call spproj%os_stk%set(cnt,'stk', absfname)
                    ! transfer particle orientations, re-indexing into the local project
                    fromp = spproj_glob%os_stk%get_fromp(istk)
                    top   = spproj_glob%os_stk%get_top(istk)
                    !$omp parallel do private(iptcl,kptcl) default(shared)
                    do iptcl = fromp, top
                        kptcl = jptcl + iptcl - fromp + 1
                        call spproj%os_ptcl2D%transfer_ori(kptcl, spproj_glob%os_ptcl2D, iptcl)
                        call spproj%os_ptcl2D%set_stkind(kptcl, cnt)
                    enddo
                    !$omp end parallel do
                    ! update local fromp/top in the sub-project stack record
                    jptcl  = jptcl  + n
                    cfromp = ctop   + 1
                    ctop   = cfromp + n - 1
                    call spproj%os_stk%set(cnt,'fromp', cfromp)
                    call spproj%os_stk%set(cnt,'top',   ctop)
                    ! register micrograph in project_list for microchunked2D tracking
                    prec%id         = project_list%size() + 1
                    prec%projname   = projfile
                    prec%micind     = cnt
                    prec%nptcls     = n
                    prec%nptcls_sel = stk_nptcls(istk)
                    prec%included   = .false.
                    call project_list%push_back(prec)
                enddo
                ! strip any prior 2D results before writing the sub-project
                call spproj%os_ptcl2D%delete_2Dclustering(keepshifts=.false., keepcls=.false.)
                call spproj%write(projfile)
                ! report selected-particle count for this chunk
                call project_list%slice(chunks_map(ichunk,1), chunks_map(ichunk,2), project_list_slice)
                nptcls_chunk = project_list_slice%get_nptcls_sel_tot()
                write(logfhandle,'(A,I6,A,I8,A,I8)') &
                    '>>> CHUNK ', ichunk, ' : ', nptcls_chunk, ' / ', nptcls_tot, ' SELECTED PARTICLES'
                call project_list_slice%kill
            enddo
            call spproj%kill
            deallocate(stk_all_nptcls, stk_nptcls, chunks_map)
        end subroutine generate_chunk_projects

    end subroutine exec_stream_cluster2D_microchunked

end module simple_stream_cluster2D_microchunked
