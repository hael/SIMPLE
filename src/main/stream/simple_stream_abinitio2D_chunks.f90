!@descr: split a project into independent abinitio2D subset jobs
!
! This module implements a bounded manual-selection workflow:
!
! (1) The input project is split into NCHUNKS stack-bound subsets, balanced by
!     active particle count as closely as stack boundaries allow.
! (2) Each subset is written as a temporary project and submitted as an
!     independent abinitio2D job.
! (3) No automatic class-average selection, matching, merging, or global
!     reclustering is performed.
! (4) The final artifacts are the solved chunk projects:
!         chunk_1/chunk.simple, chunk_2/chunk.simple, ...
!
! In short:
!
! Split project -> run independent abinitio2D jobs -> leave solved chunk projects.
!
module simple_stream_abinitio2D_chunks
use simple_stream_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_abinitio2D_chunks
  contains
    procedure :: execute => exec_stream_abinitio2D_chunks
end type stream_abinitio2D_chunks

contains

    subroutine exec_stream_abinitio2D_chunks( self, cline )
        class(stream_abinitio2D_chunks), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        character(len=*), parameter :: DIR_PROJS = trim(PATH_HERE)//'spprojs/'
        integer,          parameter :: WAITTIME  = 5
        integer, allocatable :: nptcls_per_chunk_vec(:), chunk_rec_fromto(:,:)
        type(rec_list)   :: project_list
        type(string)     :: fname
        type(parameters) :: params
        type(sp_project) :: spproj_glob
        integer          :: ichunk, nstks, nptcls, nptcls_tot, ntot_chunks, ncls_chunk
        call cline%set('oritype',      'ptcl2D')
        call cline%set('autoscale',    'yes')
        call cline%set('remove_chunks','no')
        call cline%set('objfun',       'euclid')
        call cline%set('numlen',       5)
        call cline%set('sigma_est',    'global')
        if( .not. cline%defined('refine') ) call cline%set('refine', 'snhc_smpl')
        call cline%set('nthr2D',       cline%get_iarg('nthr'))
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('center')         ) call cline%set('center',        'yes')
        if( .not. cline%defined('center_type')    ) call cline%set('center_type',   'seg')
        if( .not. cline%defined('walltime')       ) call cline%set('walltime',      29*60) ! 29 minutes
        if( .not. cline%defined('rank_cavgs')     ) call cline%set('rank_cavgs',    'no')
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls',300)
        if( .not. cline%defined('nparts')         ) call cline%set('nparts',        1)
        ! parse
        call params%new(cline)
        if( params%nchunks < 1 ) THROW_HARD('nchunks must be >= 1; exec_abinitio2D_chunks')
        if( params%nptcls_per_cls < 1 ) THROW_HARD('nptcls_per_cls must be >= 1; exec_abinitio2D_chunks')
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
            THROW_HARD('No particles found in project file: '//params%projfile%to_char()//'; exec_abinitio2D_chunks')
        endif
        ! subset project packaging
        call generate_chunk_projects
        ! initialize common abinitio2D chunk command line
        nptcls_per_chunk = nint(real(sum(nptcls_per_chunk_vec)) / real(ntot_chunks))
        params%ncls      = max(1, floor(real(nptcls_per_chunk) / real(params%nptcls_per_cls)))
        call init_chunk_clustering( params, cline, spproj_glob )
        numlen = params%numlen
        call del_file(POOL_DIR//CLUSTER2D_FINISHED)
        call cline_cluster2D_chunk%set('center', params%center)
        if( cline%defined('center_type') ) call cline_cluster2D_chunk%set('center_type', params%center_type)
        call cline_cluster2D_chunk%delete('minits')
        call cline_cluster2D_chunk%delete('maxits')
        call cline_cluster2D_chunk%delete('extr_iter')
        call cline_cluster2D_chunk%delete('extr_lim')
        call cline_cluster2D_chunk%set('rank_cavgs', params%rank_cavgs)
        ! Reinitialize each chunk with a class count derived from its own active-particle count.
        do ichunk = 1,ntot_chunks
            ncls_chunk = max(1, floor(real(nptcls_per_chunk_vec(ichunk)) / real(params%nptcls_per_cls)))
            call init_one_chunk(ichunk, ncls_chunk)
        enddo
        params%nthr2D = params%nthr ! ?? cf. Joe
        ! Submit independent abinitio2D subset jobs sequentially.
        do ichunk = 1,ntot_chunks
            call submit_one_chunk(ichunk)
            call wait_one_chunk(ichunk)
        enddo
        write(logfhandle,'(A)')'>>> INDEPENDENT ABINITIO2D SUBSET PROJECTS:'
        do ichunk = 1,ntot_chunks
            fname = chunks(ichunk)%get_projfile_fname()
            write(logfhandle,'(A,I6,A,A)')'>>> SUBSET ', ichunk, ' PROJECT: ', fname%to_char()
        enddo
        ! cleanup; chunk_N folders are the requested outputs and are intentionally retained.
        do ichunk = 1,ntot_chunks
            call chunks(ichunk)%kill
        enddo
        call spproj_glob%kill
        call project_list%kill
        call simple_rmdir(STDERROUT_DIR)
        call simple_rmdir(DIR_PROJS)
        call simple_rmdir(DIR_SNAPSHOT)
        call del_file(POOL_DIR//POOL_PROJFILE)
        call simple_rmdir(SIGMAS_DIR)
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_ABINITIO2D_CHUNKS NORMAL STOP ****')

    contains

        subroutine init_one_chunk( ichunk, ncls_chunk )
            integer, intent(in) :: ichunk, ncls_chunk
            type(cmdline) :: cline_chunk
            cline_chunk = cline_cluster2D_chunk
            call cline_chunk%set('ncls', ncls_chunk)
            call chunks(ichunk)%kill
            call chunks(ichunk)%init_chunk(params, cline_chunk, ichunk, spproj_glob)
            call cline_chunk%kill
        end subroutine init_one_chunk

        subroutine submit_one_chunk( ichunk )
            integer, intent(in) :: ichunk
            type(rec_list) :: project_list_slice
            call project_list%slice(chunk_rec_fromto(ichunk,1), chunk_rec_fromto(ichunk,2), project_list_slice)
            call chunks(ichunk)%generate(project_list_slice)
            call chunks(ichunk)%analyze2D(makecavgs=.false.)
            call project_list_slice%kill
        end subroutine submit_one_chunk

        subroutine wait_one_chunk( ichunk )
            integer, intent(in) :: ichunk
            do
                call chunks(ichunk)%display_iter
                if( chunks(ichunk)%has_converged() ) exit
                call sleep(WAITTIME)
            enddo
        end subroutine wait_one_chunk

        subroutine generate_chunk_projects
            type(sp_project)     :: spproj
            type(project_rec)    :: prec
            integer, allocatable :: stk_nptcls(:), stk_all_nptcls(:), active_stks(:), active_sel(:), chunks_map(:,:)
            type(string)         :: fname, absfname, path, projname, projfile
            integer :: active_count, assigned, cnt, end_active, end_limit, iactive, ichunk
            integer :: istk, iptcl, jptcl, kptcl, fromp, top, cfromp, ctop, n, nstks_chunk
            integer :: remaining_chunks, start_active, target, next_sel
            logical :: has_mic
            call simple_mkdir(DIR_PROJS)
            allocate(stk_all_nptcls(nstks), stk_nptcls(nstks), source=0)
            active_count = 0
            do istk = 1,nstks
                if( spproj_glob%os_stk%get_state(istk) == 0 ) cycle
                if( spproj_glob%os_stk%get_int(istk,'nptcls') == 0 ) cycle
                fromp = spproj_glob%os_stk%get_fromp(istk)
                top   = spproj_glob%os_stk%get_top(istk)
                stk_all_nptcls(istk) = top - fromp + 1
                do iptcl = fromp,top
                    if( spproj_glob%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    stk_nptcls(istk) = stk_nptcls(istk) + 1
                enddo
                if( stk_nptcls(istk) > 0 ) active_count = active_count + 1
            enddo
            nptcls_tot = sum(stk_nptcls)
            write(logfhandle,'(A,I8)')'>>> # OF STACKS            : ', nstks
            write(logfhandle,'(A,I8)')'>>> # OF ACTIVE STACKS     : ', active_count
            write(logfhandle,'(A,I8)')'>>> # OF ACTIVE PARTICLES  : ', nptcls_tot
            if( nptcls_tot == 0 ) THROW_HARD('No active particles found in project; exec_abinitio2D_chunks')
            if( params%nchunks > active_count )then
                write(logfhandle,*) 'requested nchunks / active stacks: ', params%nchunks, active_count
                THROW_HARD('nchunks cannot exceed the number of active stacks for stack-bound subset splitting')
            endif
            ntot_chunks = params%nchunks
            allocate(active_stks(active_count), active_sel(active_count))
            cnt = 0
            do istk = 1,nstks
                if( stk_nptcls(istk) == 0 ) cycle
                cnt = cnt + 1
                active_stks(cnt) = istk
                active_sel(cnt)  = stk_nptcls(istk)
            enddo
            allocate(chunks_map(ntot_chunks,2), nptcls_per_chunk_vec(ntot_chunks), chunk_rec_fromto(ntot_chunks,2), source=0)
            iactive = 1
            assigned = 0
            do ichunk = 1,ntot_chunks
                start_active     = iactive
                remaining_chunks = ntot_chunks - ichunk + 1
                if( ichunk == ntot_chunks )then
                    end_active = active_count
                else
                    target     = nint(real(nptcls_tot - assigned) / real(remaining_chunks))
                    cnt        = 0
                    end_active = start_active - 1
                    end_limit  = active_count - (remaining_chunks - 1)
                    do while( end_active < end_limit )
                        next_sel = active_sel(end_active + 1)
                        if( cnt > 0 )then
                            if( abs(cnt - target) <= abs(cnt + next_sel - target) ) exit
                        endif
                        cnt        = cnt + next_sel
                        end_active = end_active + 1
                        if( cnt >= target ) exit
                    enddo
                    if( end_active < start_active ) end_active = start_active
                endif
                chunks_map(ichunk,1)         = active_stks(start_active)
                chunks_map(ichunk,2)         = active_stks(end_active)
                nptcls_per_chunk_vec(ichunk) = sum(active_sel(start_active:end_active))
                assigned = assigned + nptcls_per_chunk_vec(ichunk)
                iactive  = end_active + 1
            enddo
            write(logfhandle,'(A,I8)')'>>> # OF REQUESTED CHUNKS  : ', ntot_chunks
            spproj%compenv = spproj_glob%compenv
            spproj%jobproc = spproj_glob%jobproc
            call spproj%projinfo%new(1, is_ptcl=.false.)
            path = CWD_GLOB//'/'//DIR_PROJS
            call spproj%projinfo%set(1,'cwd', path)
            has_mic = spproj_glob%os_mic%get_noris() > 0
            do ichunk = 1,ntot_chunks
                projname = 'subset_'//int2str_pad(ichunk,6)
                projfile = path//projname//METADATA_EXT
                call spproj%projinfo%set(1,'projname', projname)
                call spproj%projinfo%set(1,'projfile', projfile)
                nstks_chunk = 0
                nptcls      = 0
                do istk = chunks_map(ichunk,1),chunks_map(ichunk,2)
                    if( stk_nptcls(istk) == 0 ) cycle
                    nstks_chunk = nstks_chunk + 1
                    nptcls      = nptcls + stk_all_nptcls(istk)
                enddo
                call spproj%os_stk%new(nstks_chunk, is_ptcl=.false.)
                call spproj%os_mic%new(nstks_chunk, is_ptcl=.false.)
                call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                cnt   = 0
                ctop  = 0
                jptcl = 0
                chunk_rec_fromto(ichunk,1) = project_list%size() + 1
                do istk = chunks_map(ichunk,1),chunks_map(ichunk,2)
                    if( stk_nptcls(istk) == 0 ) cycle
                    cnt = cnt + 1
                    n   = stk_all_nptcls(istk)
                    if( has_mic )then
                        call spproj%os_mic%transfer_ori(cnt, spproj_glob%os_mic, istk)
                    else
                        call spproj%os_mic%set(cnt,'nptcls', n)
                        call spproj%os_mic%set_state(cnt, spproj_glob%os_stk%get_state(istk))
                    endif
                    call spproj%os_stk%transfer_ori(cnt, spproj_glob%os_stk, istk)
                    call spproj%os_stk%getter(cnt,'stk',fname)
                    absfname = simple_abspath(fname)
                    call spproj%os_stk%set(cnt,'stk',absfname)
                    fromp = spproj_glob%os_stk%get_fromp(istk)
                    top   = spproj_glob%os_stk%get_top(istk)
                    !$omp parallel do private(iptcl,kptcl) default(shared)
                    do iptcl = fromp,top
                        kptcl = jptcl + iptcl - fromp + 1
                        call spproj%os_ptcl2D%transfer_ori(kptcl, spproj_glob%os_ptcl2D, iptcl)
                        call spproj%os_ptcl2D%set_stkind(kptcl, cnt)
                    enddo
                    !$omp end parallel do
                    jptcl  = jptcl + n
                    cfromp = ctop + 1
                    ctop   = cfromp + n - 1
                    call spproj%os_stk%set(cnt,'fromp',cfromp)
                    call spproj%os_stk%set(cnt,'top',  ctop)
                    prec%id         = project_list%size() + 1
                    prec%projname   = projfile
                    prec%micind     = cnt
                    prec%nptcls     = n
                    prec%nptcls_sel = stk_nptcls(istk)
                    prec%included   = .false.
                    call project_list%push_back(prec)
                enddo
                chunk_rec_fromto(ichunk,2) = project_list%size()
                call spproj%os_ptcl2D%delete_2Dclustering(keepshifts=.false., keepcls=.false.)
                call spproj%write(projfile)
                write(logfhandle,'(A,I6,A,I8,A,I8,A,I6,A,I6)')'>>> SUBSET ', ichunk, ' : ', &
                    &nptcls_per_chunk_vec(ichunk), ' / ', nptcls_tot, ' ACTIVE PARTICLES; STACKS ', &
                    &chunks_map(ichunk,1), ' - ', chunks_map(ichunk,2)
            enddo
            call spproj%kill
            deallocate(stk_all_nptcls, stk_nptcls, active_stks, active_sel, chunks_map)
        end subroutine generate_chunk_projects

    end subroutine exec_stream_abinitio2D_chunks

end module simple_stream_abinitio2D_chunks
