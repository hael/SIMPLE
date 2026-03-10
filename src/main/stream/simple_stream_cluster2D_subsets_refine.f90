!@descr: high-level batch testing of 2D refinement within clusters of 2D averages from subsets
!
! This module implements an incremental 2D clustering policy with the following guarantees:
!
! (1) Particles are processed in bounded subsets (“chunks”), rather than as a monolithic dataset.
! (2) Chunks are constructed deterministically from stacks/micrographs, subject to:
!         * a maximum number of micrographs (nmics)
!         * a maximum number of particles (maxnptcls)
! (3) Each chunk is independently 2D-clustered, using identical clustering parameters.
! (4) Results are merged sequentially until either:
!         * the total number of non-junk class averages surpasses a limit set by NCLS_SET,
!         * all chunks have been merged
! (5) Only one merge frontier exists at any time (no parallel merges)
! (6) The merged solution is clustered using cluster_cavgs
! (7) Particles within each cluster are independently 2D-clustered, using identical clustering parameters except for:
!         * the number of classes is set to the maximum of NCLS_MULT * number classes in the cluster and MIN_NCLS
! (8) The resulting clustering solutions are merged — the final state contains the class assignments from the per cluster 2D
!
! In short:
!
! “Split → independently cluster → merge until threshold met -> cluster cavgs -> independently cluster particles within each cluster -> report final cluster membership in project file”
!
! This is a testing policy, not a global optimization.
!
module simple_stream_cluster2D_subsets_refine
use simple_stream_api
use simple_commanders_cavgs,      only: commander_cluster_cavgs, commander_match_cavgs
use simple_commanders_abinitio2D, only: commander_abinitio2D
use simple_imgarr_utils,          only: read_cavgs_into_imgarr, dealloc_imgarr
use simple_strategy2D_utils,      only: flag_non_junk_cavgs
use simple_projfile_utils,        only: merge_chunk_projfiles
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_cluster2D_subsets_refine
  contains
    procedure :: execute => exec_stream_cluster2D_subsets_refine
end type stream_cluster2D_subsets_refine

contains

    ! cluster2D_subsets: splits into chunks & analyze2D them
    subroutine exec_stream_cluster2D_subsets_refine( self, cline )
        class(stream_cluster2D_subsets_refine), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        character(len=*),                           parameter :: DIR_PROJS = trim(PATH_HERE)//'spprojs/'
        integer,                                    parameter :: WAITTIME    = 5
        real,                                       parameter :: NCLS_MULT   = 1.2
        integer,                                    parameter :: MIN_NCLS    = 5
        integer,                                    parameter :: NCLS_SET    = 200
        integer,                                  allocatable :: global_ptcl_state_map(:)
        integer,                                  allocatable :: nptcls_per_chunk_vec(:), cluster_assignments(:)
        type(rec_list)   :: project_list
        type(string)     :: fname
        type(parameters) :: params
        type(sp_project) :: spproj_glob
        type(rec_list)   :: chunkslist
        type(oris)       :: cls2d_merged
        type(class_frcs) :: frcs
        integer          :: ichunk, nstks, nptcls, nptcls_tot, ntot_chunks, n_non_junk, ic, id, nc, icls, box4frc, iter, subproc_fhandle, ios
        logical          :: all_chunks_submitted
        call cline%set('oritype',      'ptcl2D')
        call cline%set('wiener',       'full')
        call cline%set('autoscale',    'yes')
        call cline%set('remove_chunks','no')
        call cline%set('objfun',       'euclid')
        call cline%set('numlen',       5)
        call cline%set('sigma_est',    'global')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('nchunks',      1)
        call cline%set('nthr2D',       cline%get_iarg('nthr'))
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('center')         ) call cline%set('center',        'yes')
        if( .not. cline%defined('center_type')    ) call cline%set('center_type',   'seg')
        if( .not. cline%defined('walltime')       ) call cline%set('walltime',      29*60) ! 29 minutes
        if( .not. cline%defined('rank_cavgs')     ) call cline%set('rank_cavgs',    'no')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',         100)
        if( .not. cline%defined('maxnptcls')      ) call cline%set('maxnptcls',     100000)
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls',200)
        if( .not. cline%defined('nparts')         ) call cline%set('nparts',        1)
        ! parse
        call params%new(cline)
        ! read strictly required fields
        call spproj_glob%read_non_data_segments(params%projfile)
        call spproj_glob%read_segment('mic',   params%projfile)
        call spproj_glob%read_segment('stk',   params%projfile)
        call spproj_glob%read_segment('ptcl2D',params%projfile)
        ! clear previous clustering
        call clear_previous()
        ! sanity checks
        nstks  = spproj_glob%os_stk%get_noris()
        nptcls = spproj_glob%get_nptcls()
        if( spproj_glob%os_mic%get_noris() > 0 )then
            if( spproj_glob%os_mic%get_noris() /= nstks )then
                THROW_HARD('Inconsistent # of micrographs and stacks, use prune_project.')
            endif
        endif
        if( nptcls == 0 )then
            THROW_HARD('No particles found in project file: '//params%projfile%to_char()//'; exec_cluster2d_subsets')
        endif
        ! projects packaging
        call generate_chunk_projects
        ! Update to global parameters prior to 2D inititalization
        nptcls_per_chunk = nint(real(sum(nptcls_per_chunk_vec)) / real(ntot_chunks))    ! average value
        params%ncls      = floor(real(nptcls_per_chunk) / real(params%nptcls_per_cls))
        ! General streaming initialization
        call init_chunk_clustering( params, cline, spproj_glob )
        ! Updates folllowing streaming init
        numlen = params%numlen
        call del_file(POOL_DIR//CLUSTER2D_FINISHED)
        call cline_cluster2D_chunk%set('center', params%center)
        if( cline%defined('center_type') ) call cline_cluster2D_chunk%set('center_type', params%center_type)
        call cline_cluster2D_chunk%delete('minits')
        call cline_cluster2D_chunk%delete('maxits')
        call cline_cluster2D_chunk%delete('extr_iter')
        call cline_cluster2D_chunk%delete('extr_lim')
        call cline_cluster2D_chunk%set('rank_cavgs', params%rank_cavgs)
        ! re-init with updated command-lines
        do ichunk = 1,params%nchunks
            call chunks(ichunk)%kill
            call chunks(ichunk)%init_chunk(params, cline_cluster2D_chunk, ichunk, spproj_glob)
        enddo
        params%nthr2D = params%nthr ! ?? cf. Joe
        ! Main loop
        ichunk     = 0  ! # of chunks that have been submitted
        n_non_junk = 0  ! # cumulative number of non junk classes
        iter       = 0
        all_chunks_submitted = .false.
        do
            ! sequential chunk prep & submission
            if( .not.all_chunks_submitted )then
                if( chunks(1)%is_available() )then
                    ichunk = ichunk + 1
                    nptcls_per_chunk = nptcls_per_chunk_vec(ichunk) ! is a variable
                    call analyze2D_new_chunks(params, project_list, .false.)
                    all_chunks_submitted = ichunk == ntot_chunks
                endif
            endif
            ! convergence
            call check_completed_chunks
            if( allocated(converged_chunks) )then
                ! # of converged chunks
                nc = size(converged_chunks)
                ! update global list
                do ic = 1,nc
                    fname = converged_chunks(ic)%get_projfile_fname()
                    id    = converged_chunks(ic)%get_id()
                    ! cleanup
                    call converged_chunks(ic)%kill
                    ! reject junk averages from chunk
                    call reject_junk(id, fname)
                    ! update chunkslist and set processed to true
                    call chunkslist%push2chunk_list(fname, id, .true.)
                enddo
                deallocate(converged_chunks)
            endif
            ! Completion
            if( n_non_junk > NCLS_SET ) then
                write(logfhandle,'(A,I8)')'>>> TARGET # NON-JUNK CLASSES REACHED : ', n_non_junk
                exit
            else if( chunkslist%size() == ntot_chunks )then
                if( all(chunkslist%get_processed_flags()) ) exit
            endif
            call sleep(WAITTIME)
        end do
        ! merge processed chunks
        call merge_processed_chunks()
        ! allocate global_ptcl_state_map
        global_ptcl_state_map = spproj_glob%os_ptcl2D%get_all_asint('state')
        ! cluster merged cavgs
        call cluster_merged_cavgs()
        call print_stats()
        ! open subproc_fhandle
        open(UNIT=subproc_fhandle, FILE=SIMPLE_SUBPROC_OUT, IOSTAT=ios, ACTION='WRITE', STATUS='NEW', POSITION='APPEND')
        do while(iter<2)
            iter = iter + 1
            ! run abinitio2D on particles within each cluster separately
            icls = 1
            cluster_assignments = spproj_glob%os_cls2D%get_all_asint('cluster')
            ! remove existing class assignments
            call spproj_glob%os_ptcl2D%set_all2single('class', 0)
            call spproj_glob%os_ptcl3D%set_all2single('class', 0)
            call spproj_glob%write()
            !------------------------------------------------------------------
            ! will combine, optimise and move these to external subroutine once methodology converges
            do ic=1, maxval(cluster_assignments)
                write(logfhandle, '(A,I8)')'>>> RUNNING ABINITIO2D ON PARTICLES IN CLUSTER ', ic
                logfhandle = subproc_fhandle
                call abinitio2D_cluster(ic)
                logfhandle = OUTPUT_UNIT
            enddo
            do ic=1, maxval(cluster_assignments)
                call merge_cluster_cavgs(ic)
            enddo
            do ic=1, maxval(cluster_assignments)
                call merge_cluster_frcs(ic)
            enddo
            !------------------------------------------------------------------
            ! write project
            call spproj_glob%write()
            ! recluster merged cavgs
            call cluster_merged_cavgs()
            ! write iteration project
            call spproj_glob%write()
            ! print some stats
            call print_stats()
            ! cleanup
            call cls2d_merged%kill()
        end do
        ! cleanup
        call flush(subproc_fhandle)
        close(subproc_fhandle)
        call spproj_glob%kill
        call chunkslist%kill
        call simple_rmdir(STDERROUT_DIR)
        call simple_rmdir(DIR_PROJS)
        call simple_rmdir(DIR_SNAPSHOT)
        call del_file(POOL_DIR//POOL_PROJFILE)
        call simple_rmdir(SIGMAS_DIR)
        call qsys_cleanup(params)
        ! graceful end
        call simple_end('**** SIMPLE_CLUSTER2D_SUBSETS NORMAL STOP ****')

    contains

        subroutine check_completed_chunks
            type(stream_chunk), allocatable :: tmpchunks(:)
            integer :: ichunk, jchunk, nthr2D, n
            logical :: chunk_complete
            if( .not. l_stream2D_active ) return
            do ichunk = 1,params%nchunks
                if( chunks(ichunk)%is_available() ) cycle
                chunk_complete = .false.
                if( chunks(ichunk)%to_analyze2D() )then
                    chunk_complete = chunks(ichunk)%has_converged()
                else
                    THROW_HARD('Fatal Error 1')
                endif
                if( chunk_complete )then
                    ! read & print out info
                    call chunks(ichunk)%display_iter
                    ! book-keeping
                    if( allocated(converged_chunks) )then
                        n = size(converged_chunks)
                        allocate(tmpchunks(n+1),source=[converged_chunks(:), chunks(ichunk)])
                        do jchunk = 1,n
                            call converged_chunks(jchunk)%kill
                        enddo
                        deallocate(converged_chunks)
                        allocate(converged_chunks(n+1),source=tmpchunks)
                        do jchunk = 1,n+1
                            call tmpchunks(jchunk)%kill
                        enddo
                        deallocate(tmpchunks)
                    else
                        ! first item
                        allocate(converged_chunks(1),source=[chunks(ichunk)])
                    endif
                    ! reinit and deal with nthr2D != nthr
                    glob_chunk_id = glob_chunk_id + 1
                    ! deal with nthr2d .ne. nthr
                    nthr2D = params%nthr2D
                    params%nthr2D = cline_cluster2D_chunk%get_iarg('nthr')
                    call chunks(ichunk)%init_chunk(params, cline_cluster2D_chunk, glob_chunk_id, pool_proj)
                    params%nthr2D = nthr2D
                endif
            enddo
        end subroutine check_completed_chunks

        subroutine generate_chunk_projects
            type(sp_project)     :: spproj
            type(project_rec)    :: prec
            integer, allocatable :: stk_nptcls(:), stk_all_nptcls(:), chunks_map(:,:)
            type(rec_list)       :: project_list_slice
            type(string)         :: fname,absfname,path,projname,projfile
            integer :: cnt, ichunk, istk, iptcl,jptcl,kptcl,fromp,top,cfromp,ctop,n,cnt_stk
            call simple_mkdir(DIR_PROJS)
            ! Mapping particles/stacks, number of chunks
            allocate(stk_all_nptcls(nstks),stk_nptcls(nstks),source=0)
            ntot_chunks = 0
            cnt         = 0
            cnt_stk     = 0
            do istk = 1,nstks
                if( (spproj_glob%os_stk%get_state(istk)==0) .or. (spproj_glob%os_stk%get_int(istk,'nptcls')==0) )cycle
                fromp = spproj_glob%os_stk%get_fromp(istk)
                top   = spproj_glob%os_stk%get_top(istk)
                do iptcl = fromp,top
                    stk_all_nptcls(istk) = stk_all_nptcls(istk) + 1 ! including state=0
                    if( spproj_glob%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    stk_nptcls(istk) = stk_nptcls(istk) + 1 ! excluding state=0
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
            write(logfhandle,'(A,I8)')'>>> # OF STACKS          : ', nstks
            write(logfhandle,'(A,I8)')'>>> # OF PARTICLES       : ', nptcls_tot
            write(logfhandle,'(A,I8)')'>>> # OF AVAILABLE CHUNKS: ', ntot_chunks
            if( cline%defined('maxnchunks') ) ntot_chunks = min(params%maxnchunks, ntot_chunks)
            ! chunks map, leftovers are abandoned
            allocate(chunks_map(ntot_chunks,2),nptcls_per_chunk_vec(ntot_chunks),source=0)
            cnt     = 0
            cnt_stk = 0
            ichunk  = 1
            do istk = 1,nstks
                if( ichunk > ntot_chunks ) exit
                if( cnt==0 ) chunks_map(ichunk,1) = istk
                cnt     = cnt + stk_nptcls(istk)
                cnt_stk = cnt_stk + 1
                if( (cnt_stk >= params%nmics) .or. (cnt >= params%maxnptcls) )then
                    chunks_map(ichunk,2)         = istk
                    nptcls_per_chunk_vec(ichunk) = cnt
                    ichunk  = ichunk + 1
                    cnt     = 0
                    cnt_stk = 0
                endif
            enddo
            write(logfhandle,'(A)')'>>> CHUNKS MAP: '
            spproj%compenv = spproj_glob%compenv
            spproj%jobproc = spproj_glob%jobproc
            call spproj%projinfo%new(1, is_ptcl=.false.)
            path = CWD_GLOB//'/'//DIR_PROJS
            call spproj%projinfo%set(1,'cwd', path)
            ! stacks/ptcls transfer
            do ichunk = 1,ntot_chunks
                projname = int2str_pad(ichunk,6)
                projfile = path//projname//METADATA_EXT
                call spproj%projinfo%set(1,'projname', projname)
                call spproj%projinfo%set(1,'projfile', projfile)
                nstks  = chunks_map(ichunk,2) - chunks_map(ichunk,1) + 1
                nptcls = sum(stk_all_nptcls(chunks_map(ichunk,1):chunks_map(ichunk,2)))
                call spproj%os_stk%new(nstks, is_ptcl=.false.)
                call spproj%os_mic%new(nstks, is_ptcl=.false.)
                call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                cnt   = 0
                ctop  = 0
                jptcl = 0 ! particle index in local project
                do istk = chunks_map(ichunk,1),chunks_map(ichunk,2)
                    cnt = cnt + 1
                    n   = stk_all_nptcls(istk)
                    ! micrograph
                    if( spproj_glob%os_mic%get_noris() > 0 )then
                        call spproj%os_mic%transfer_ori(cnt, spproj_glob%os_mic, istk)
                    else
                        call spproj%os_mic%set(cnt,'nptcls', n)
                        call spproj%os_mic%set_state(cnt,spproj_glob%os_stk%get_state(istk))
                    endif
                    ! stack
                    call spproj%os_stk%transfer_ori(cnt, spproj_glob%os_stk, istk)
                    call spproj%os_stk%getter(cnt,'stk',fname)
                    absfname = simple_abspath(fname)
                    call spproj%os_stk%set(cnt,'stk',absfname)
                    ! particle
                    fromp = spproj_glob%os_stk%get_fromp(istk)
                    top   = spproj_glob%os_stk%get_top(istk)
                    !$omp parallel do private(iptcl,kptcl) default(shared)
                    do iptcl = fromp,top
                        kptcl = jptcl + iptcl-fromp+1
                        call spproj%os_ptcl2D%transfer_ori(kptcl, spproj_glob%os_ptcl2D, iptcl)
                        call spproj%os_ptcl2D%set_stkind(kptcl, cnt)
                    enddo
                    !$omp end parallel do
                    jptcl  = jptcl + n
                    cfromp = ctop + 1
                    ctop   = cfromp + n - 1
                    call spproj%os_stk%set(cnt,'fromp',cfromp)
                    call spproj%os_stk%set(cnt,'top',  ctop)
                    prec%projname   = projfile
                    prec%micind     = cnt
                    prec%nptcls     = n
                    prec%nptcls_sel = stk_nptcls(istk)
                    prec%included   = .false.
                    call project_list%push_back(prec)
                enddo
                ! remove previous parameters
                call spproj%os_ptcl2D%delete_2Dclustering(keepshifts=.false., keepcls=.false.)
                call spproj%write(projfile)
                ! generate a slice of project_list corresponding to the chunk
                call project_list%slice(chunks_map(ichunk,1),chunks_map(ichunk,2),project_list_slice)
                ! get total # particles in chunk
                nptcls = project_list_slice%get_nptcls_sel_tot()
                write(logfhandle,'(A,I8,A,I8,A,I8)')'>>> CHUNK ID; # OF PARTICLES  : ',  ichunk, ' ; ',nptcls,' / ',sum(stk_nptcls)
                call project_list_slice%kill
            enddo
            call spproj%kill
        end subroutine generate_chunk_projects

        ! subjects first chunk of to cluster_cavgs
        subroutine cluster_merged_cavgs
            type(commander_cluster_cavgs) :: xcluster_cavgs
            type(cmdline)                 :: cline_cluster_cavgs
            type(string)                  :: dir, cwd, path
            integer,          allocatable :: clusters(:), states(:)
            ! updates directory structure
            if( iter == 0 ) then
                dir  = "cluster_cavgs"
            else
                dir = "cluster_cavgs_" // int2str(iter)
            endif
            path = string('../') // basename(params%projfile)
            call simple_mkdir(dir)
            ! execute
            cline_cluster_cavgs = cline
            call cline_cluster_cavgs%set('mkdir',   'no')
            call cline_cluster_cavgs%set('prg',     'cluster_cavgs')
            call cline_cluster_cavgs%set('projfile', path)
            call cline_cluster_cavgs%delete('nparts')
            call simple_chdir(dir)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            call xcluster_cavgs%execute(cline_cluster_cavgs)
            call simple_chdir('..')
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            ! reread master project
            call spproj_glob%kill()
            call spproj_glob%read(basename(params%projfile))
            ! update global_ptcl_state_map
            clusters = spproj_glob%os_ptcl2D%get_all_asint('cluster')
            where( clusters == 0 ) global_ptcl_state_map = 0
            ! reapply global_ptcl_state_map
            call spproj_glob%os_ptcl2D%set_all('state', global_ptcl_state_map)
            call spproj_glob%os_ptcl3D%set_all('state', global_ptcl_state_map)
            deallocate(clusters)
            clusters = spproj_glob%os_cls2D%get_all_asint('cluster')
            allocate(states(size(clusters)))
            states = 1
            where( clusters == 0 ) states = 0
            call spproj_glob%os_cls2D%set_all('state', states)
            ! cleanup
            deallocate(clusters, states)
            call cline_cluster_cavgs%kill
        end subroutine cluster_merged_cavgs

        ! Sets junk classes to state=0 and updates global n_non_junk
        subroutine reject_junk( id, fname )
            integer,          intent(in) :: id
            type(string),     intent(in) :: fname
            type(sp_project)             :: spproj_rj
            type(image),     allocatable :: cavg_imgs(:)
            logical,         allocatable :: l_non_junk(:)
            integer,         allocatable :: states(:)
            real,              parameter :: LP_BIN = 20.
            real                         :: smpd, mskrad
            integer                      :: ldim(3), box, n_junk
            ! read chunk project file
            call spproj_rj%read(fname)
            ! load cavgs and associated parameters
            cavg_imgs  = read_cavgs_into_imgarr(spproj_rj)
            smpd       = cavg_imgs(1)%get_smpd()
            ldim       = cavg_imgs(1)%get_ldim()
            box        = ldim(1)
            mskrad     = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam/smpd)
            ! flag non-junk in cavgs
            n_junk = 0
            call flag_non_junk_cavgs( cavg_imgs, LP_BIN, mskrad, l_non_junk )
            if( allocated(l_non_junk) ) then
                n_junk     = size(l_non_junk) - count(l_non_junk)
                n_non_junk = n_non_junk + count(l_non_junk)
                ! retrieve states
                states = spproj_rj%os_cls2D%get_all_asint('state')
                ! set state=0 where junk
                where( .not. l_non_junk ) states = 0
                ! map states to particles
                call spproj_rj%map_cavgs_selection(states)
                ! write project
                call spproj_rj%write()
                ! cleanup
                deallocate(l_non_junk)
                deallocate(states)
            endif
            write(logfhandle,'(A,I8,A,I8)')'>>> ', n_junk, ' JUNK CLASSES REJECTED FROM CHUNK ', id
            ! cleanup
            call spproj_rj%kill()
            call dealloc_imgarr(cavg_imgs)
        end subroutine reject_junk

        subroutine merge_processed_chunks()
            type(string),   allocatable :: projfiles(:)
            logical, allocatable        :: mask(:)
            mask      = chunkslist%get_processed_flags()
            projfiles = chunkslist%get_projfiles([1, chunkslist%size()])
            projfiles = pack(projfiles, mask)
            call merge_chunk_projfiles(projfiles, string('.'), spproj_glob)
            call spproj_glob%write(basename(params%projfile))
        end subroutine merge_processed_chunks

        subroutine abinitio2D_cluster(ic)
            integer,                   intent(in) :: ic
            type(commander_abinitio2D)            :: xabinitio2D
            type(sp_project)                      :: spproj_cluster
            type(cmdline)                         :: cline_abinitio2D
            type(string)                          :: projfile, dir, cwd, path
            integer,                  allocatable :: states(:), clusters(:)
            write(logfhandle,'(A,I8)')'>>> RUNNING ABINITIO_2D ON PARTICLES FROM CLUSTER ', ic
            ! set dir and projfile
            if(iter == 0) then
                dir = 'cluster_'//int2str(ic)
            else
                dir = 'cluster_'//int2str(ic)//'_iter'//int2str(iter)
            endif
            projfile = 'cluster_'//int2str(ic)//'.simple'
            path     = string('../') // projfile
            ! copy global project file
            call spproj_cluster%copy(spproj_glob)
            ! retrieve clusters field for particles
            clusters = spproj_cluster%os_ptcl2D%get_all_asint('cluster')
            ! allocate states
            allocate(states(size(clusters)))
            states = 0
            ! update states to match ic
            where( clusters == ic ) states = 1
            ! apply global states
            where( global_ptcl_state_map == 0 ) states = 0
            ! update ptcl2D and 3D
            call spproj_cluster%os_ptcl2D%set_all('state',  states)
            call spproj_cluster%os_ptcl3D%set_all('state',  states)
            ! cleanup
            deallocate(clusters, states)
            ! retrieve clusters and states fields for classes
            clusters = spproj_cluster%os_cls2D%get_all_asint('cluster')
            allocate(states(size(clusters)))
            states = 1
            ! update states to match 
            where( clusters /= ic ) states = 0
            ! destroy existing cls2D
            call spproj_cluster%os_cls2D%kill()
            call spproj_cluster%os_cls3D%kill()
            call spproj_cluster%os_ptcl2D%set_all2single('class',  0)
            call spproj_cluster%os_ptcl3D%set_all2single('class',  0)
            call spproj_cluster%os_ptcl2D%set_all2single('e1', 0)
            call spproj_cluster%os_ptcl2D%set_all2single('e2', 0)
            call spproj_cluster%os_ptcl2D%set_all2single('e3', 0)
            call spproj_cluster%os_ptcl2D%set_all2single('x',  0)
            call spproj_cluster%os_ptcl2D%set_all2single('y',  0)
            call spproj_cluster%os_ptcl3D%set_all2single('e1', 0)
            call spproj_cluster%os_ptcl3D%set_all2single('e2', 0)
            call spproj_cluster%os_ptcl3D%set_all2single('e3', 0)
            call spproj_cluster%os_ptcl3D%set_all2single('x',  0)
            call spproj_cluster%os_ptcl3D%set_all2single('y',  0)
            ! write cluster projfile
            call spproj_cluster%write(projfile)
            ! set abinitio2D parameters
            call cline_abinitio2D%set('prg',                                       'abinitio2D')
            call cline_abinitio2D%set('mkdir',                                             'no')
            call cline_abinitio2D%set('projfile',                                          path)
            call cline_abinitio2D%set('mskdiam',                                 params%mskdiam)
            call cline_abinitio2D%set('nthr',                                       params%nthr)
            call cline_abinitio2D%set('nparts',                                   params%nparts)
            call cline_abinitio2D%set('ncls',     max(ceiling(count(states == 1) * NCLS_MULT), MIN_NCLS))
            call cline_abinitio2D%printline()
            ! run abinitio2D
            call simple_mkdir(dir)
            ! execute
            call simple_chdir(dir)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            call xabinitio2D%execute(cline_abinitio2D)
            call simple_chdir('..')
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            ! cleanup
            deallocate(clusters, states)
            call spproj_cluster%kill()
            call cline_abinitio2D%kill()
        end subroutine abinitio2D_cluster

        subroutine merge_cluster_cavgs(ic)
            integer,                   intent(in) :: ic
            type(image)                           :: img
            type(sp_project)                      :: spproj_cluster
            type(oris)                            :: cls2d_backup
            type(string)                          :: projfile
            type(string)                          :: stkname, evenname, oddname, cavgs
            integer,                  allocatable :: classes_glob(:), classes(:)
            real,                     allocatable :: e3s_glob(:), e3s(:), xs_glob(:), xs(:), ys_glob(:), ys(:)
            integer                               :: ldim(3), ncls, i
            real                                  :: smpd
            write(logfhandle,'(A,I8)')'>>> MERGING CLS2D FROM CLUSTER ', ic
            ! projfile
            projfile = 'cluster_'//int2str(ic)//'.simple'
            if(iter == 0) then
                cavgs = 'merged_cluster_cavgs'//MRC_EXT
            else
                cavgs = 'merged_cluster_cavgs_iter'//int2str(iter)//MRC_EXT
            endif
            call spproj_cluster%read(projfile)
            call spproj_cluster%get_cavgs_stk(stkname, ncls, smpd, imgkind='cavg')
            call find_ldim_nptcls(stkname, ldim, ncls)
            ldim(3) = 1
            ncls    = spproj_cluster%os_cls2D%get_noris()
            call img%new(ldim, smpd)
            evenname = add2fbody(stkname, MRC_EXT, '_even')
            oddname  = add2fbody(stkname, MRC_EXT, '_odd')
            call cls2d_backup%copy(cls2d_merged)
            call cls2d_merged%kill()
            call cls2d_merged%new(cls2d_backup%get_noris() + ncls, is_ptcl=.false.)
            do i=1, cls2d_backup%get_noris()
                call cls2d_merged%transfer_ori(i, cls2d_backup, i)
            enddo
            classes      = spproj_cluster%os_ptcl2D%get_all_asint('class')
            classes_glob = spproj_glob%os_ptcl2D%get_all_asint('class')
            e3s          = spproj_cluster%os_ptcl2D%get_all('e3')
            e3s_glob     = spproj_glob%os_ptcl2D%get_all('e3')
            xs           = spproj_cluster%os_ptcl2D%get_all('x')
            xs_glob      = spproj_glob%os_ptcl2D%get_all('x')
            ys           = spproj_cluster%os_ptcl2D%get_all('y')
            ys_glob      = spproj_glob%os_ptcl2D%get_all('y')
            do i=1, ncls
                call img%read(stkname,i)
                call img%write(cavgs,icls)
                call img%read(evenname,i)
                call img%write(get_fbody(basename(cavgs), fname2ext(cavgs))//'_even'//MRC_EXT, icls)
                call img%read(oddname,i)
                call img%write(get_fbody(basename(cavgs), fname2ext(cavgs))//'_odd'//MRC_EXT, icls)
                where(classes == i) 
                    classes_glob = icls
                    e3s_glob     = e3s
                    xs_glob      = xs
                    ys_glob      = ys
                end where
                call cls2d_merged%transfer_ori(icls, spproj_cluster%os_cls2D, i)
                call cls2d_merged%set_class(icls, icls)
                icls = icls+1   
            enddo
            call spproj_glob%os_ptcl2D%set_all('class',          classes_glob)
            call spproj_glob%os_ptcl3D%set_all('class',          classes_glob)
            call spproj_glob%os_ptcl2D%set_all('e3',                 e3s_glob)
            call spproj_glob%os_ptcl3D%set_all('e3',                 e3s_glob)
            call spproj_glob%os_ptcl2D%set_all('x',                   xs_glob)
            call spproj_glob%os_ptcl3D%set_all('x',                   xs_glob)
            call spproj_glob%os_ptcl2D%set_all('y',                   ys_glob)
            call spproj_glob%os_ptcl3D%set_all('y',                   ys_glob)
            call spproj_glob%os_ptcl2D%set_all('state', global_ptcl_state_map)
            call spproj_glob%os_ptcl3D%set_all('state', global_ptcl_state_map)
            call spproj_glob%add_cavgs2os_out(cavgs, smpd, imgkind='cavg')
            call spproj_glob%os_cls2D%copy(cls2d_merged)
            ! cleanup
            call cls2d_backup%kill()
        end subroutine merge_cluster_cavgs

        subroutine merge_cluster_frcs(ic)
            integer,                   intent(in) :: ic
            type(sp_project)                      :: spproj_cluster
            type(class_frcs)                      :: frcs_cluster
            type(string)                          :: projfile, frc_fname
            integer                               :: ncls, i
            real                                  :: smpd
            write(logfhandle,'(A,I8)')'>>> MERGING FRCS FROM CLUSTER ', ic
            projfile = 'cluster_'//int2str(ic)//'.simple'
            call spproj_cluster%read(projfile)
            call spproj_cluster%get_frcs(frc_fname, 'frc2D')
            call frcs_cluster%read(frc_fname)
            smpd = spproj_cluster%get_smpd()
            ncls = spproj_cluster%os_cls2D%get_noris()
            if( ic == 1 )then
                icls = 1
                box4frc = frcs_cluster%get_box()
                call frcs%new(spproj_glob%os_cls2D%get_noris(), box4frc, smpd)
            endif
            do i=1, ncls
                if( spproj_cluster%os_cls2D%get_state(i) > 0 ) then
                    call frcs%set_frc(icls, frcs_cluster%get_frc(i,  box4frc))
                endif
                icls = icls+1
            enddo
            call frcs%write(string(trim(FRCS_FILE)))
            call spproj_glob%add_frcs2os_out(string(trim(FRCS_FILE)), 'frc2D')
            call frcs_cluster%kill()
            call spproj_cluster%kill()
        end subroutine merge_cluster_frcs

        subroutine print_stats()
            integer, allocatable :: classes_glob(:), clusters_glob(:)
            integer              :: icls2d
            clusters_glob = spproj_glob%os_cls2D%get_all_asint('cluster')
            classes_glob  = spproj_glob%os_ptcl2D%get_all_asint('class')
            write(logfhandle, '(A,I6,A,I6)')'NPTCL:', spproj_glob%os_ptcl2D%get_noris(), ' NPTCL REJECTED:', count(global_ptcl_state_map == 0)
            do icls2d=1, spproj_glob%os_cls2D%get_noris()
                write(logfhandle, '(A,I6,A,I6,A,I6)')'CLS:', icls2d, ' POP:', count(classes_glob == icls2d), ' CLUSTER:', clusters_glob(icls2d)
            end do
        end subroutine print_stats

        subroutine clear_previous()
            ! nuke previous alignments
            call spproj_glob%os_ptcl2D%set_all2single('e1', 0)
            call spproj_glob%os_ptcl2D%set_all2single('e2', 0)
            call spproj_glob%os_ptcl2D%set_all2single('e3', 0)
            call spproj_glob%os_ptcl2D%set_all2single('x',  0)
            call spproj_glob%os_ptcl2D%set_all2single('y',  0)
            call spproj_glob%os_ptcl3D%set_all2single('e1', 0)
            call spproj_glob%os_ptcl3D%set_all2single('e2', 0)
            call spproj_glob%os_ptcl3D%set_all2single('e3', 0)
            call spproj_glob%os_ptcl3D%set_all2single('x',  0)
            call spproj_glob%os_ptcl3D%set_all2single('y',  0)
        end subroutine clear_previous

    end subroutine exec_stream_cluster2D_subsets_refine

end module simple_stream_cluster2D_subsets_refine