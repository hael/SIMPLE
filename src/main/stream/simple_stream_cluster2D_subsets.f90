module simple_stream_cluster2D_subsets
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: params_glob
use simple_projfile_utils, only: merge_chunk_projfiles
use simple_sp_project,     only: sp_project
use simple_stream_chunk,   only: stream_chunk
use simple_stream_chunk2D_utils
use simple_stream2D_state
use simple_commanders_cluster2D
use simple_qsys_funs
use simple_rec_list
use simple_stream_utils
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_cluster2D_subsets
  contains
    procedure :: execute      => exec_stream_cluster2D_subsets
end type stream_cluster2D_subsets

contains

    ! cluster2D_subsets: splits into chunks & analyze2D them
    subroutine exec_stream_cluster2D_subsets( self, cline )
        class(stream_cluster2D_subsets), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        character(len=*), parameter :: DIR_PROJS   = trim(PATH_HERE)//'spprojs/'
        integer,          parameter :: WAITTIME    = 5
        integer,        allocatable :: nptcls_per_chunk_vec(:)
        type(rec_list)      :: project_list
        type(string)        :: fname
        type(parameters)    :: params
        type(sp_project)    :: spproj_glob
        type(rec_list)      :: chunkslist
        integer             :: ichunk, nstks, nptcls, nptcls_tot, ntot_chunks, ic, id, nc
        logical             :: all_chunks_submitted
        call cline%set('oritype',      'ptcl2D')
        call cline%set('wiener',       'full')
        call cline%set('autoscale',    'yes')
        call cline%set('remove_chunks','no')
        call cline%set('reject_cls',   'no')
        call cline%set('objfun',       'euclid')
        call cline%set('numlen',       5)
        call cline%set('sigma_est',    'global')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('algorithm',    'abinitio2D')
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
        params_glob%ncls = floor(real(nptcls_per_chunk) / real(params_glob%nptcls_per_cls))
        ! General streaming initialization
        call init_chunk_clustering( cline, spproj_glob )
        ! Updates folllowing streaming init
        numlen = params%numlen
        call del_file(POOL_DIR//CLUSTER2D_FINISHED)
        call cline_cluster2D_chunk%set('center', params%center)
        if( cline%defined('center_type') ) call cline_cluster2D_chunk%set('center_type', params%center_type)
        call cline_cluster2D_chunk%delete('minits')
        call cline_cluster2D_chunk%delete('maxits')
        call cline_cluster2D_chunk%delete('extr_iter')
        call cline_cluster2D_chunk%delete('extr_lim')
        call cline_cluster2D_chunk%delete('cc_iters')
        call cline_cluster2D_chunk%set('rank_cavgs', params%rank_cavgs)
        ! re-init with updated command-lines
        do ichunk = 1,params_glob%nchunks
            call chunks(ichunk)%kill
            call chunks(ichunk)%init_chunk(ichunk, cline_cluster2D_chunk, spproj_glob)
        enddo
        params_glob%nthr2D = params_glob%nthr ! ?? cf. Joe
        ! Main loop
        ichunk = 0  ! # of chunks that have been submitted
        all_chunks_submitted = .false.
        do
            ! sequential chunk prep & submission
            if( .not.all_chunks_submitted )then
                if( chunks(1)%is_available() )then
                    ichunk = ichunk + 1
                    nptcls_per_chunk = nptcls_per_chunk_vec(ichunk) ! is a variable
                    call analyze2D_new_chunks(project_list, .false.)
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
                    call chunkslist%push2chunk_list(fname, id, .false.)
                    ! cleanup
                    call converged_chunks(ic)%kill
                    ! run cluster_cavgs on first item of chunkslist
                    call cluster_chunk_cavgs
                    ! run match_cavgs on last item of chunkslist
                    call match_sets
                enddo
                deallocate(converged_chunks)
            endif
            ! Completion
            if( chunkslist%size() == ntot_chunks )then
                if( all(chunkslist%get_processed_flags()) ) exit
            endif
            call sleep(WAITTIME)
        end do
        ! no final project
        ! cleanup
        call spproj_glob%kill
        call chunkslist%kill
        call simple_rmdir(STDERROUT_DIR)
        call simple_rmdir(DIR_PROJS)
        call simple_rmdir(DIR_SNAPSHOT)
        call del_file(POOL_DIR//POOL_PROJFILE)
        call simple_rmdir(SIGMAS_DIR)
        call qsys_cleanup
        ! graceful end
        call simple_end('**** SIMPLE_CLUSTER2D_SUBSETS NORMAL STOP ****')
    contains

        subroutine check_completed_chunks
            type(stream_chunk), allocatable :: tmpchunks(:)
            integer :: ichunk, jchunk, nthr2D, n
            logical :: chunk_complete
            if( .not. l_stream2D_active ) return
            do ichunk = 1,params_glob%nchunks
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
                    nthr2D = params_glob%nthr2D
                    params_glob%nthr2D = cline_cluster2D_chunk%get_iarg('nthr')
                    call chunks(ichunk)%init_chunk(glob_chunk_id, cline_cluster2D_chunk, pool_proj)
                    params_glob%nthr2D = nthr2D
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
        subroutine cluster_chunk_cavgs
            type(commander_cluster_cavgs) :: xcluster_cavgs
            type(cmdline)                 :: cline_cluster_cavgs
            type(sp_project)              :: spproj
            type(string)                  :: path, tmpl, cwd
            logical, allocatable          :: l_processed(:)
            type(chunk_rec)               :: crec
            if( chunkslist%size() /= 1 ) return
            l_processed = chunkslist%get_processed_flags()
            if( l_processed(1)         )return
            ! updates directory structure
            tmpl = trim(DIR_SET)//int2str(1)
            call simple_mkdir(tmpl)
            call chunkslist%at(1, crec)
            call merge_chunk_projfiles([crec%projfile], tmpl, spproj, projname_out=tmpl, write_proj=.true.)
            ! update path
            crec%projfile = simple_abspath(tmpl//'/'//tmpl//METADATA_EXT)
            ! execute
            cline_cluster_cavgs = cline
            call cline_cluster_cavgs%set('mkdir',   'no')
            call cline_cluster_cavgs%set('prg',     'cluster_cavgs')
            call cline_cluster_cavgs%set('projfile', basename(crec%projfile))
            call cline_cluster_cavgs%delete('nparts')
            path = stemname(crec%projfile)
            call simple_chdir(path)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            call xcluster_cavgs%execute_safe(cline_cluster_cavgs)
            call simple_chdir('..')
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            crec%processed = .true.
            call chunkslist%replace_at(1, crec)
            ! cleanup
            call spproj%kill
            call cline_cluster_cavgs%kill
            if( allocated(l_processed) ) deallocate(l_processed)
        end subroutine cluster_chunk_cavgs

        ! subjects last chunk(s) chunk of to match_cavgs
        subroutine match_sets
            type(commander_match_cavgs) :: xmatch_cavgs
            type(cmdline)               :: cline_match_cavgs
            logical, allocatable        :: l_processed(:)
            type(chunk_rec) :: crec_id, crec_id_min_1
            type(string)  :: tmpl, cwd
            integer       :: id
            if( chunkslist%size() < 2 )return
            id = chunkslist%size()
            l_processed = chunkslist%get_processed_flags()
            if(      l_processed(id)    ) return
            if( .not. l_processed(id-1) ) return
            ! updates directory structure
            tmpl = trim(DIR_SET)//int2str(id)
            call simple_mkdir(tmpl)
            ! execute
            call chunkslist%at(id,   crec_id)
            call chunkslist%at(id-1, crec_id_min_1)
            cline_match_cavgs = cline
            call cline_match_cavgs%set('mkdir',           'no')
            call cline_match_cavgs%set('prg',             'match_cavgs')
            call cline_match_cavgs%set('projfile',        simple_abspath(crec_id_min_1%projfile))
            call cline_match_cavgs%set('projfile_target', simple_abspath(crec_id%projfile))
            call cline_match_cavgs%set('projfile_merged', 'set_'//int2str(id)//METADATA_EXT)
            call cline_match_cavgs%delete('nparts')
            call simple_chdir(tmpl)
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            call xmatch_cavgs%execute_safe(cline_match_cavgs)
            call simple_chdir('..')
            call simple_getcwd(cwd)
            CWD_GLOB = cwd%to_char()
            crec_id%processed = .true.
            crec_id%projfile  = simple_abspath(tmpl//'/'//'set_'//int2str(id)//METADATA_EXT)
            call chunkslist%replace_at(id, crec_id)
            ! cleanup
            call cline_match_cavgs%kill
        end subroutine match_sets

    end subroutine exec_stream_cluster2D_subsets

end module simple_stream_cluster2D_subsets