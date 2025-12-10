module simple_commanders_cavgs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_builder,        only: builder, build_glob
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_stack_io,       only: stack_io
use simple_strategy2D_utils
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_rank_cavgs
    contains
        procedure :: execute      => exec_rank_cavgs
end type commander_rank_cavgs

type, extends(commander_base) :: commander_cluster_cavgs
    contains
        procedure :: execute      => exec_cluster_cavgs
end type commander_cluster_cavgs

type, extends(commander_base) :: commander_select_clusters
    contains
        procedure :: execute      => exec_select_clusters
end type commander_select_clusters

type, extends(commander_base) :: commander_match_cavgs
    contains
        procedure :: execute      => exec_match_cavgs
end type commander_match_cavgs

contains

    subroutine exec_rank_cavgs( self, cline )
        class(commander_rank_cavgs), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(sp_project), target  :: spproj
        class(oris),      pointer :: os_ptr => null()
        type(parameters)     :: params
        type(oris)           :: clsdoc_ranked
        type(stack_io)       :: stkio_r, stkio_w
        type(image)          :: img
        integer, allocatable :: order(:)
        real,    allocatable :: vals(:), rstates(:)
        logical, allocatable :: mask(:)
        integer              :: ldim(3), ncls, icls, i, rank
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'cls2D')
        if( .not. cline%defined('flag') )    call cline%set('flag',    'res')
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        call find_ldim_nptcls(params%stk, ldim, ncls)
        if( trim(params%oritype) .eq. 'cls2D' )then
            os_ptr => spproj%os_cls2D
        else
            os_ptr => spproj%os_cls3D
        endif
        params%ncls = ncls
        if( os_ptr%get_noris() == params%ncls )then
            ! all we need to do is fetch from classdoc in projfile &
            ! order according to resolution
            call img%new([params%box,params%box,1], params%smpd)
            call clsdoc_ranked%new(params%ncls, is_ptcl=.false.)
            select case(trim(params%flag))
                case('res','corr','pop')
                    ! supported
                case DEFAULT
                    THROW_HARD('Unsupported cavg flag: '//trim(params%flag))
            end select
            vals    = os_ptr%get_all(params%flag)
            order   = (/(icls,icls=1,params%ncls)/)
            rstates = os_ptr%get_all('state')
            mask    = rstates > 0.5
            call hpsort(vals, order)
            select case(trim(params%flag))
                case('corr')
                    where( rstates < 0.5 ) vals = -1.0
                    call reverse(order)
                case DEFAULT
                    ! done
            end select
            call stkio_r%open(params%stk, params%smpd, 'read', bufsz=params%ncls)
            call stkio_r%read_whole ! because need asynchronous access
            call stkio_w%open(params%outstk, params%smpd, 'write', box=ldim(1), bufsz=params%ncls)
            rank = 0
            do icls=1,params%ncls
                i = order(icls)
                if( mask(i) )then
                    rank = rank + 1
                    call clsdoc_ranked%set(rank, 'class',    i)
                    call clsdoc_ranked%set(rank, 'rank',     rank)
                    select case(trim(params%flag))
                    case('corr')
                        call clsdoc_ranked%set(rank, 'corr', os_ptr%get(i, 'corr'))
                        write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,1x,f6.2)') 'CLASS:', i,&
                        &'RANK:', rank ,'CORR:', os_ptr%get(i, 'corr')
                    case DEFAULT
                        call clsdoc_ranked%set(rank, 'pop',  os_ptr%get(i,  'pop'))
                        call clsdoc_ranked%set(rank, 'res',  os_ptr%get(i,  'res'))
                        call clsdoc_ranked%set(rank, 'corr', os_ptr%get(i, 'corr'))
                        call clsdoc_ranked%set(rank, 'w',    os_ptr%get(i,    'w'))
                        write(logfhandle,'(a,1x,i5,1x,a,1x,i5,1x,a,i5,1x,a,1x,f6.2)') 'CLASS:', i,&
                        &'RANK:', rank ,'POP:', os_ptr%get_int(i, 'pop'),&
                        &'RES:', os_ptr%get(i, 'res')
                    end select
                    call flush(logfhandle)
                endif
                call stkio_r%get_image(i, img)
                call stkio_w%write(icls, img)
            end do
            call stkio_r%close
            call stkio_w%close
            call clsdoc_ranked%write(string('classdoc_ranked.txt'))
        else
            ! nothing to do
        endif
        ! end gracefully
        call clsdoc_ranked%kill
        call img%kill
        call spproj%kill
        nullify(os_ptr)
        call simple_end('**** SIMPLE_RANK_CAVGS NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rank_cavgs

    subroutine exec_cluster_cavgs( self, cline )
        use simple_clustering_utils, only: cluster_dmat
        class(commander_cluster_cavgs), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        logical,          parameter   :: DEBUG = .true.
        real,             parameter   :: SCORE_THRES_INCL = 75.
        integer,          parameter   :: NCLUST_MAX = 65
        type(image),      allocatable :: cavg_imgs(:)
        type(image)                   :: img_msk
        real,             allocatable :: mm(:,:), dmat(:,:), resvals_tmp(:), resvals(:)
        logical,          allocatable :: l_non_junk(:), l_msk(:,:,:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), inds(:)
        integer,          allocatable :: clspops(:), states(:), labels4write(:), inds_glob(:)
        type(clust_info), allocatable :: clust_info_arr(:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        integer                       :: ncls, ncls_sel, icls, cnt, nptcls
        integer                       :: i, nclust, iclust, nptcls_good
        integer                       :: ldim(3), pop, box
        real                          :: frac_good, mskrad
        real                          :: oa_min, oa_max, smpd
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',       10.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune')   ) call cline%set('prune',    'no')
        ! master parameters
        call params%new(cline)
        ! read project file
        call spproj%read(params%projfile)
        ncls = spproj%os_cls2D%get_noris()
        ! prep class average stack
        if( trim(params%have_clustering).eq.'yes' )then
            labels    = spproj%os_cls2D%get_all_asint('cluster')
            allocate(l_non_junk(size(labels)), source=labels > 0)
            labels    = pack(labels,mask=l_non_junk)
            ncls_sel  = count(l_non_junk)
            cavg_imgs = read_cavgs_into_imgarr(spproj, l_non_junk)
            smpd      = cavg_imgs(1)%get_smpd()
            ldim      = cavg_imgs(1)%get_ldim()
            box       = ldim(1)
            mskrad    = min(real(box/2) - COSMSKHALFWIDTH - 1., 0.5 * params%mskdiam/smpd)
            clspops   = spproj%os_cls2D%get_all_asint('pop')
            clspops   = pack(clspops, mask=l_non_junk)
            clsinds   = pack((/(i,i=1,ncls)/), mask=l_non_junk)
            ! create the stuff needed in the loop
            allocate(mm(ncls_sel,2), source=0.)
            ! prep mask
            call img_msk%new([box,box,1], smpd)
            img_msk = 1.
            call img_msk%mask(mskrad, 'hard')
            l_msk = img_msk%bin2logical()
            call img_msk%kill
            !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
            do i = 1, ncls_sel
                ! normalization
                call cavg_imgs(i)%norm_within(l_msk)
                ! mask
                call cavg_imgs(i)%mask(mskrad, 'soft', backgr=0.)
                ! stash minmax
                mm(i,:) = cavg_imgs(i)%minmax(mskrad)
            end do
            !$omp end parallel do
            deallocate(l_msk)
        else
            call prep_cavgs4clustering(spproj, cavg_imgs, params%mskdiam, clspops, clsinds, l_non_junk, mm )
        endif
        ncls_sel = size(cavg_imgs)
        smpd     = cavg_imgs(1)%get_smpd()
        ldim     = cavg_imgs(1)%get_ldim()
        allocate(resvals(ncls_sel), source=0.)
        do i = 1, ncls_sel
            resvals(i) = spproj%os_cls2D%get(clsinds(i), 'res')
        enddo
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! calculate overall minmax
        oa_min      = minval(mm(:,1))
        oa_max      = maxval(mm(:,2))
        ! calculate distance matrix
        dmat = calc_cluster_cavgs_dmat(params, cavg_imgs, [oa_min,oa_max], params%clust_crit)
        ! cluster
        if( trim(params%have_clustering).eq.'yes' )then
            nclust = maxval(labels)
            call cluster_dmat(dmat, 'refine', nclust, i_medoids, labels)
        else if( cline%defined('ncls') )then
            nclust = params%ncls
            call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
        else
            call cluster_dmat( dmat, 'aprop', nclust, i_medoids, labels, nclust_max=NCLUST_MAX)
            if( nclust < 3 )then
                nclust = 3
                call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
            endif
        endif
        clust_info_arr = align_and_score_cavg_clusters( params, dmat, cavg_imgs, clspops, i_medoids, labels )
        ! re-create cavg_imgs
        call dealloc_imgarr(cavg_imgs)
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        ! write aligned clusters
        call write_aligned_cavgs(labels, cavg_imgs, clust_info_arr, 'cluster_aligned', params%ext%to_char())
        ! write un-aligned clusters
        call write_imgarr(ncls_sel, cavg_imgs, labels, 'cluster', params%ext%to_char() )
        ! update project
        call spproj%os_ptcl2D%transfer_class_assignment(spproj%os_ptcl3D)
        call spproj%os_cls2D%set_all2single('cluster',  0)
        call spproj%os_cls3D%set_all2single('cluster',  0)
        call spproj%os_cls2D%set_all2single('accept',   0)
        call spproj%os_cls3D%set_all2single('accept',   0)
        call spproj%os_ptcl2D%set_all2single('cluster', 0)
        call spproj%os_ptcl3D%set_all2single('cluster', 0)
        do iclust = 1, nclust
            do icls = 1, ncls_sel 
                if( labels(icls) == iclust )then
                    call spproj%os_cls2D%set(clsinds(icls),'cluster',iclust)                          ! 2D class field
                    call spproj%os_cls3D%set(clsinds(icls),'cluster',iclust)                          ! 3D class field
                    call spproj%os_cls2D%set(clsinds(icls),'accept', clust_info_arr(iclust)%good_bad) ! 2D class accepted field
                    call spproj%os_cls3D%set(clsinds(icls),'accept', clust_info_arr(iclust)%good_bad) ! 3D class accepted field
                    call spproj%os_ptcl2D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 2D particle field
                    call spproj%os_ptcl3D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 3D particle field
                endif
            enddo
        enddo
        ! generate ranked class averages by ordering according to class resolution within clusters
        allocate(inds_glob(ncls_sel), source=0)
        cnt = 0
        do iclust = 1, nclust
            pop         = count(labels == iclust)
            inds        = mask2inds(labels == iclust) 
            resvals_tmp = pack(resvals, mask=labels == iclust)
            call hpsort(resvals_tmp, inds)
            do i = 1, pop
                cnt = cnt + 1
                inds_glob(cnt) = inds(i)
            enddo
            deallocate(inds, resvals_tmp)
        enddo
        ! write ranked_cavgs
        call write_imgarr(cavg_imgs, string('ranked_cavgs')//params%ext%to_char(), inds_glob)
        deallocate(inds_glob)
        ! report cluster info
        do iclust = 1, nclust
            write(logfhandle,'(A,A,f5.1,A,f5.1,A,f5.1,A,f5.1,A,f5.1,A,I3)') 'cluster_ranked'//int2str_pad(iclust,2)//'.mrc',&
            &' resolution(A) ',   clust_info_arr(iclust)%res,& 
            &' resscore(%) ',     clust_info_arr(iclust)%resscore,& 
            &' homogeneity(%) ',  clust_info_arr(iclust)%homogeneity,&
            &' clustscore(%) ',   clust_info_arr(iclust)%clustscore,&
            &' jointscore(%) ',   clust_info_arr(iclust)%jointscore,&
            &' good_bad_assign ', clust_info_arr(iclust)%good_bad
        end do
        ! check number of particles selected
        nptcls      = sum(clust_info_arr(:)%nptcls)
        nptcls_good = sum(clust_info_arr(:)%nptcls, mask=clust_info_arr(:)%good_bad == 1)
        frac_good   = real(nptcls_good)  / real(nptcls)
        write(logfhandle,'(a,1x,f8.2)') '% PARTICLES CLASSIFIED AS 1ST RATE: ', frac_good  * 100.
        ! translate to state array
        allocate(states(ncls), source=0)
        do icls = 1, ncls_sel
            if( clust_info_arr(labels(icls))%good_bad == 1 ) states(clsinds(icls)) = 1
        end do
        ! write selection
        allocate(labels4write(ncls_sel), source=0)
        do icls = 1, ncls_sel
            labels4write(icls) = clust_info_arr(labels(icls))%good_bad
        end do
        ! write selection
        call write_selected_cavgs(ncls_sel, cavg_imgs, labels4write, params%ext%to_char())
        ! map selection to project
        call spproj%map_cavgs_selection(states)
        ! optional pruning
        if( trim(params%prune).eq.'yes') call spproj%prune_particles
        ! this needs to be a full write as many segments are updated
        call spproj%write(params%projfile)
        ! destruct
        call spproj%kill
        call dealloc_imgarr(cavg_imgs)
        ! deallocate anything not specifically allocated above
        deallocate(clust_info_arr, l_non_junk, labels, clsinds, i_medoids)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')
    end subroutine exec_cluster_cavgs

    subroutine exec_select_clusters( self, cline )
        class(commander_select_clusters), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer, allocatable :: clustinds(:)
        type(string) :: selflag
        integer :: iclust, nclust_sel, nclust_max
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('prune')       ) call cline%set('prune',           'yes')
        if( .not. cline%defined('select_flag') ) call cline%set('select_flag', 'cluster')
        ! master parameters
        call params%new(cline)
        if( cline%defined('clustind') )then
            nclust_sel = 1
            allocate(clustinds(nclust_sel), source=params%clustind)
        else
            clustinds = list_of_ints2arr(params%clustinds)
            nclust_sel = size(clustinds)
        endif
        do iclust = 1, nclust_sel
            print *, 'selected cluster: ', clustinds(iclust)
        end do
        ! check flag
        selflag = trim(params%select_flag)
        select case(selflag%to_char())
            case('cluster','class')
                ! all good
            case DEFAULT
                THROW_HARD('Unsupported flag for cluster selection (SELECT_FLAG)')
        end select
        ! read project file
        call spproj%read(params%projfile)
        nclust_max = spproj%os_ptcl2D%get_n(selflag%to_char())
        if( any(clustinds > nclust_max) ) THROW_HARD('Maximum cluster index value: '//int2str(nclust_max)//' exceeded!')
        do iclust = 1, nclust_max
            if( any(clustinds == iclust) )then
                ! set state=1 to flag inclusion
                call spproj%os_cls2D%set_field2single(selflag%to_char(),  iclust, 'state', 1) ! 2D class field
                call spproj%os_cls3D%set_field2single(selflag%to_char(),  iclust, 'state', 1) ! 3D class field
                call spproj%os_ptcl2D%set_field2single(selflag%to_char(), iclust, 'state', 1) ! 2D particle field
                call spproj%os_ptcl3D%set_field2single(selflag%to_char(), iclust, 'state', 1) ! 3D particle field
            else
                ! set state=0 to flag exclusion
                call spproj%os_cls2D%set_field2single(selflag%to_char(),  iclust, 'state', 0) ! 2D class field
                call spproj%os_cls3D%set_field2single(selflag%to_char(),  iclust, 'state', 0) ! 3D class field
                call spproj%os_ptcl2D%set_field2single(selflag%to_char(), iclust, 'state', 0) ! 2D particle field
                call spproj%os_ptcl3D%set_field2single(selflag%to_char(), iclust, 'state', 0) ! 3D particle field
            endif
        end do
        ! prune
        if( trim(params%prune).eq.'yes') call spproj%prune_particles
        ! write project
        call spproj%write(params%projfile)
        ! destruct
        call spproj%kill
        ! end gracefully
        call simple_end('**** SIMPLE_SELECT_CLUSTERS_CAVGS NORMAL STOP ****')
    end subroutine exec_select_clusters

    subroutine exec_match_cavgs( self, cline )
        use simple_projfile_utils, only: merge_chunk_projfiles
        class(commander_match_cavgs), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj_ref, spproj_match, spproj_merged
        type(cmdline)    :: cline_cluster_cavgs
        type(commander_cluster_cavgs)  :: xcluster_cavgs
        character(len=*),  parameter   :: TMPPROJFILE = 'tmp_projfile_match_cavgs.simple'
        type(string)                   :: chunk_fnames(2), folder
        type(image),       allocatable :: cavg_imgs_ref(:), cavg_imgs_match(:)
        integer,           allocatable :: clspops_ref(:), clsinds_ref(:), clspops_match(:), clsinds_match(:), labels(:)
        integer,           allocatable :: states(:), labels_match(:), states_map(:)
        logical,           allocatable :: l_non_junk_ref(:), l_non_junk_match(:)
        real,              allocatable :: mm_ref(:, :), mm_match(:, :), dmat_clust(:, :), dmat(:, :)
        integer :: nmatch, nrefs, ldim(3), i, ncls_match, nclust, icls, iclust, imatch
        real    :: smpd, oa_minmax(2)
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',       10.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune')   ) call cline%set('prune',    'no')
        ! master parameters
        call params%new(cline)
        ! read base project file
        call spproj_ref%read(params%projfile)
        if( .not. spproj_ref%os_cls2D%isthere('cluster') ) THROW_HARD('Reference project lacks clustering information in cls2D field')
        ! read match project file
        call spproj_match%read(params%projfile_target)
        ! prepare class averages
        call prep_cavgs4clustering(spproj_ref,   cavg_imgs_ref,   params%mskdiam, clspops_ref,   clsinds_ref,   l_non_junk_ref,   mm_ref)
        call prep_cavgs4clustering(spproj_match, cavg_imgs_match, params%mskdiam, clspops_match, clsinds_match, l_non_junk_match, mm_match)
        nrefs        = size(cavg_imgs_ref)
        nmatch       = size(cavg_imgs_match)
        smpd         = cavg_imgs_ref(1)%get_smpd()
        ldim         = cavg_imgs_ref(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd  = smpd
        params%box   = ldim(1)
        params%msk   = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! extract nonzero cluster labels
        labels       = spproj_ref%os_cls2D%get_all_asint('cluster')
        labels       = pack(labels, mask=l_non_junk_ref)
        nclust       = maxval(labels)
        states       = spproj_ref%os_cls2D%get_all_asint('state')
        states       = pack(states, mask=l_non_junk_ref)
        ! calculate overall minmax
        oa_minmax(1) = minval(mm_ref(:,1))
        oa_minmax(2) = maxval(mm_ref(:,2))
        ! generate matching distance matrix
        dmat         = calc_match_cavgs_dmat(params, cavg_imgs_ref, cavg_imgs_match, oa_minmax, params%clust_crit)
        ! genrate cluster distance matrix
        allocate(dmat_clust(nclust,nmatch), source=0.)
        do iclust = 1, nclust
            do imatch = 1, nmatch
                dmat_clust(iclust,imatch) = sum(dmat(:,imatch), mask=labels == iclust) / real(count(labels == iclust))
            end do
        end do
        allocate(labels_match(nmatch), source=0)
        labels_match = minloc(dmat_clust, dim=1)
        call  write_imgarr(nmatch, cavg_imgs_match, labels_match, 'cluster_match', params%ext%to_char())
        ! update project
        call spproj_match%os_ptcl2D%transfer_class_assignment(spproj_match%os_ptcl3D)
        call spproj_match%os_cls2D%set_all2single('cluster',  0)
        call spproj_match%os_cls3D%set_all2single('cluster',  0)
        call spproj_match%os_cls2D%set_all2single('accept',   0)
        call spproj_match%os_cls3D%set_all2single('accept',   0)
        call spproj_match%os_ptcl2D%set_all2single('cluster', 0)
        call spproj_match%os_ptcl3D%set_all2single('cluster', 0)
        do iclust = 1, nclust
            do icls = 1, nmatch
                if( labels_match(icls) == iclust )then
                    call spproj_match%os_cls2D%set(clsinds_match(icls),'cluster',iclust)                               ! 2D class field
                    call spproj_match%os_cls3D%set(clsinds_match(icls),'cluster',iclust)                               ! 3D class field
                    call spproj_match%os_cls2D%set(clsinds_match(icls),'accept', find_label_state(labels_match(icls))) ! 2D class accepted field
                    call spproj_match%os_cls3D%set(clsinds_match(icls),'accept', find_label_state(labels_match(icls))) ! 3D class accepted field
                    call spproj_match%os_ptcl2D%set_field2single('class', clsinds_match(icls), 'cluster', iclust)      ! 2D particle field
                    call spproj_match%os_ptcl3D%set_field2single('class', clsinds_match(icls), 'cluster', iclust)      ! 3D particle field
                endif
            enddo
        enddo
        if( cline%defined('projfile_merged') )then
            ! merge spproj_ref & spproj_match projects
            call spproj_match%write(string(TMPPROJFILE))
            chunk_fnames(1) = params%projfile
            chunk_fnames(2) = simple_abspath(TMPPROJFILE)
            folder = PATH_HERE
            call merge_chunk_projfiles(chunk_fnames, folder, spproj_merged)
            call del_file(chunk_fnames(2))
            call spproj_merged%write(params%projfile_merged)
            ! refine joint solution
            call cline_cluster_cavgs%set('mkdir',                      'no')
            call cline_cluster_cavgs%set('have_clustering',           'yes')
            call cline_cluster_cavgs%set('prune',              params%prune)
            call cline_cluster_cavgs%set('clust_crit',    params%clust_crit)
            call cline_cluster_cavgs%set('hp',                    params%hp)
            call cline_cluster_cavgs%set('lp',                    params%lp)
            call cline_cluster_cavgs%set('mskdiam',          params%mskdiam)
            call cline_cluster_cavgs%set('nthr',                params%nthr)
            call cline_cluster_cavgs%set('projfile', params%projfile_merged)
            call xcluster_cavgs%execute_safe(cline_cluster_cavgs)
            call cline_cluster_cavgs%kill
        else
            ! translate to state array
            ncls_match = spproj_match%os_cls2D%get_noris()
            allocate(states_map(ncls_match), source=0)
            do icls = 1, nmatch
                states_map(clsinds_match(icls)) = find_label_state(labels_match(icls))
            end do
            ! map selection to project
            call spproj_match%map_cavgs_selection(states_map)
            ! optional pruning
            if( trim(params%prune).eq.'yes') call spproj_match%prune_particles
            ! this needs to be a full write as many segments are updated
            call spproj_match%write(params%projfile_target)
        endif
        ! cleanup
        call spproj_match%kill
        call spproj_ref%kill
        call spproj_merged%kill
        call dealloc_imgarr(cavg_imgs_ref)
        call dealloc_imgarr(cavg_imgs_match)
        ! end gracefully
        call simple_end('**** SIMPLE_MATCH_CAVGS NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')

    contains
        
        function find_label_state( label ) result( state )
            integer, intent(in) :: label
            integer :: state
            state = 1
            do i = 1, nrefs
                if( labels(i) == label )then
                    state = states(i)
                    return
                endif
            end do
        end function find_label_state

    end subroutine exec_match_cavgs

end module simple_commanders_cavgs
