!@descr: commanders for analysis of class averages
module simple_commanders_cavgs
use simple_commander_module_api
use simple_strategy2D_utils
use simple_imgarr_utils, only: read_cavgs_into_imgarr, dealloc_imgarr, write_imgarr, extract_imgarr, write_selected_cavgs, join_imgarrs
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

type, extends(commander_base) :: commander_cluster_cavgs_selection
    contains
        procedure :: execute      => exec_cluster_cavgs_selection
end type commander_cluster_cavgs_selection

type, extends(commander_base) :: commander_select_clusters
    contains
        procedure :: execute      => exec_select_clusters
end type commander_select_clusters

type, extends(commander_base) :: commander_match_cavgs
    contains
        procedure :: execute      => exec_match_cavgs
end type commander_match_cavgs

type, extends(commander_base) :: commander_map_cavgs_selection
  contains
    procedure :: execute      => exec_map_cavgs_selection
end type commander_map_cavgs_selection

type, extends(commander_base) :: commander_map_cavgs_states
  contains
    procedure :: execute      => exec_map_cavgs_states
end type commander_map_cavgs_states

type, extends(commander_base) :: commander_shape_rank_cavgs
  contains
    procedure :: execute      => exec_shape_rank_cavgs
end type commander_shape_rank_cavgs

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
        real,             allocatable :: mm(:,:), dmat(:,:), resvals_tmp(:), resvals(:)
        logical,          allocatable :: l_non_junk(:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), inds(:)
        integer,          allocatable :: clspops(:), states(:), labels4write(:), inds_glob(:)
        type(clust_info), allocatable :: clust_info_arr(:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        integer                       :: ncls, ncls_sel, icls, cnt, nptcls
        integer                       :: i, nclust, iclust, nptcls_good
        integer                       :: ldim(3), pop
        real                          :: frac_good
        real                          :: oa_min, oa_max, smpd
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
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
        call id_junk_and_prep_cavgs4clust(spproj, cavg_imgs, params%mskdiam, clspops, clsinds, l_non_junk, mm )
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
        if( cline%defined('ncls') )then
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
        ! communicate medoid indices to cls2D field of project (this have to be after scoring & ranking)
        call spproj%os_cls2D%set_all2single('medoid_ind',  0)
        do iclust = 1, nclust
            if( i_medoids(iclust) > 0 )then
                call spproj%os_cls2D%set(clsinds(i_medoids(iclust)), 'medoid_ind', iclust)
            endif
        end do
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

    ! to create medoid representatives of good/bad selection in cls2D field for fast matching later
    subroutine exec_cluster_cavgs_selection( self, cline )
        use simple_clustering_utils, only: cluster_dmat
        class(commander_cluster_cavgs_selection), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        logical,          parameter   :: DEBUG = .true.
        type(image),      allocatable :: cavg_imgs_state1(:), cavg_imgs_state0(:)
        real,             allocatable :: mm_state1(:,:), mm_state0(:,:), dmat_state1(:,:), dmat_state0(:,:)
        integer,          allocatable :: labels_state1(:), clsinds_state1(:), i_medoids_state1(:), clspops_state1(:)
        integer,          allocatable :: labels_state0(:), clsinds_state0(:), i_medoids_state0(:), clspops_state0(:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), clspops(:), labels_state0_copy(:)
        integer,          allocatable :: states(:)
        type(parameters) :: params
        type(sp_project) :: spproj
        integer          :: ncls, ncls_state1, ncls_state0, icls
        integer          :: nclust, nclust_state1, nclust_state0, iclust
        integer          :: ldim(3)
        real             :: oa_minmax_state1(2), oa_minmax_state0(2), smpd
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')   ) call cline%set('trs',       10.)
        if( .not. cline%defined('lp')    ) call cline%set('lp',         6.)
        ! master parameters
        call params%new(cline)
        ! read project file
        call spproj%read(params%projfile)
        ncls        = spproj%os_cls2D%get_noris()
        states      = spproj%os_cls2D%get_all_asint('state')
        call prep_cavgs4clust(spproj, cavg_imgs_state1, params%mskdiam, clspops_state1, clsinds_state1, states  > 0, mm_state1)
        call prep_cavgs4clust(spproj, cavg_imgs_state0, params%mskdiam, clspops_state0, clsinds_state0, states == 0, mm_state0)
        ncls_state1 = size(cavg_imgs_state1)
        ncls_state0 = size(cavg_imgs_state0)
        smpd        = cavg_imgs_state1(1)%get_smpd()
        ldim        = cavg_imgs_state1(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! calculate overall minmax
        oa_minmax_state1(1) = minval(mm_state1(:,1))
        oa_minmax_state1(2) = maxval(mm_state1(:,2))
        oa_minmax_state0(1) = minval(mm_state0(:,1))
        oa_minmax_state0(2) = maxval(mm_state0(:,2))
        ! calculate distance matrices
        dmat_state1 = calc_cluster_cavgs_dmat(params, cavg_imgs_state1, oa_minmax_state1, params%clust_crit)
        dmat_state0 = calc_cluster_cavgs_dmat(params, cavg_imgs_state0, oa_minmax_state0, params%clust_crit)
        ! cluster
        call cluster_dmat(dmat_state1, 'aprop', nclust_state1, i_medoids_state1, labels_state1)
        call cluster_dmat(dmat_state0, 'aprop', nclust_state0, i_medoids_state0, labels_state0)
        ! merge solution
        nclust = nclust_state1 + nclust_state0
        labels_state0_copy = labels_state0
        do iclust = 1, nclust_state0
            where(labels_state0 == iclust) labels_state0_copy = nclust_state1 + iclust
        end do
        allocate(labels(ncls),      source=[labels_state1,labels_state0_copy])
        allocate(clsinds(ncls),     source=[clsinds_state1,clsinds_state0])
        allocate(i_medoids(nclust), source=[i_medoids_state1,i_medoids_state0])
        allocate(clspops(ncls),     source=[clspops_state1,clspops_state0])
        deallocate(labels_state0_copy)
        ! communicate medoid indices to cls2D field of project
        call spproj%os_cls2D%set_all2single('medoid_ind',  0)
        do iclust = 1, nclust
            if( i_medoids(iclust) > 0 )then
                call spproj%os_cls2D%set(clsinds(i_medoids(iclust)), 'medoid_ind', iclust)
            endif
        end do
        ! update project
        call spproj%os_ptcl2D%transfer_class_assignment(spproj%os_ptcl3D)
        call spproj%os_cls2D%set_all2single('cluster',  0)
        call spproj%os_cls3D%set_all2single('cluster',  0)
        call spproj%os_cls2D%set_all2single('accept',   0)
        call spproj%os_cls3D%set_all2single('accept',   0)
        call spproj%os_ptcl2D%set_all2single('cluster', 0)
        call spproj%os_ptcl3D%set_all2single('cluster', 0)
        do iclust = 1, nclust
            do icls = 1, ncls
                if( labels(icls) == iclust )then
                    call spproj%os_cls2D%set(clsinds(icls),'cluster',iclust)                          ! 2D class field
                    call spproj%os_cls3D%set(clsinds(icls),'cluster',iclust)                          ! 3D class field
                    call spproj%os_cls2D%set(clsinds(icls),'accept', find_label_state(labels(icls)))  ! 2D class accepted field
                    call spproj%os_cls3D%set(clsinds(icls),'accept', find_label_state(labels(icls)))  ! 3D class accepted field
                    call spproj%os_ptcl2D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 2D particle field
                    call spproj%os_ptcl3D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 3D particle field
                endif
            enddo
        enddo
        ! this needs to be a full write as many segments are updated
        call spproj%write(params%projfile)
        ! destruct
        call spproj%kill
        call dealloc_imgarr(cavg_imgs_state0)
        call dealloc_imgarr(cavg_imgs_state1)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS_SELECTION NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')

        contains

        function find_label_state( label ) result( state )
            integer, intent(in) :: label
            integer :: i,state
            state = 1
            do i = 1, ncls
                if( labels(i) == label )then
                    state = states(i)
                    return
                endif
            end do
        end function find_label_state

    end subroutine exec_cluster_cavgs_selection

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
        use simple_projfile_utils,   only: merge_chunk_projfiles
        use simple_clustering_utils, only: cluster_dmat
        use simple_srch_sort_loc,    only: unique
        class(commander_match_cavgs), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        integer,           parameter   :: NCLUST_MAX = 65
        type(parameters)               :: params
        type(sp_project)               :: spproj_ref, spproj_match, spproj_merged
        type(cmdline)                  :: cline_cluster_cavgs
        type(commander_cluster_cavgs)  :: xcluster_cavgs
        type(string)                   :: chunk_fnames(2), folder, cavgs_merged
        type(image),       allocatable :: cavg_imgs_ref(:), cavg_imgs_match(:), medoid_imgs_ref(:), medoid_imgs_match(:)
        integer,           allocatable :: clspops_ref(:), clsinds_ref(:), clspops_match(:), clsinds_match(:), labels(:)
        integer,           allocatable :: i_medoids_match(:), states(:), labels_match(:), labels_med(:), labels_map(:)
        integer,           allocatable :: inds(:), medoid_inds(:), medoid_map(:), medoid_inds_unique(:)
        logical,           allocatable :: l_non_junk_ref(:), l_non_junk_match(:)
        real,              allocatable :: mm_ref(:,:), mm_match(:,:), dmat(:,:), dmat_med(:,:)
        integer :: nmatch, nrefs, ldim(3), i, nclust_match, nclust, icls, iclust, iclust_match, nmerged
        real    :: smpd, oa_minmax(2), oa_minmax_match(2)
        logical :: l_recluster
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')   ) call cline%set('trs',       10.)
        if( .not. cline%defined('lp')    ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune') ) call cline%set('prune',    'no')
        ! master parameters
        call params%new(cline)
        ! read base project file
        call spproj_ref%read(params%projfile)
        if( .not. spproj_ref%os_cls2D%isthere('cluster')    ) THROW_HARD('Reference project lacks clustering information in cls2D field')
        if( .not. spproj_ref%os_cls2D%isthere('medoid_ind') ) THROW_HARD('Reference project lacks medoid information in cls2D field')
        ! read match project file
        call spproj_match%read(params%projfile_target)
        ! prepare class averages
        call id_junk_and_prep_cavgs4clust(spproj_ref,   cavg_imgs_ref,   params%mskdiam, clspops_ref,   clsinds_ref,   l_non_junk_ref,   mm_ref)
        call id_junk_and_prep_cavgs4clust(spproj_match, cavg_imgs_match, params%mskdiam, clspops_match, clsinds_match, l_non_junk_match, mm_match)
        nrefs              = size(cavg_imgs_ref)
        nmatch             = size(cavg_imgs_match)
        nmerged            = nrefs + nmatch
        if( nmerged <= SIEVING_MATCH_CAVGS_MAX )then
            ! merge & recluster with affinity propagation 
            l_recluster = .true. 
        else
            ! cluster target (spproj_match) with affinity propagation and map onto reference solution
            l_recluster = .false.         
        endif
        smpd               = cavg_imgs_ref(1)%get_smpd()
        ldim               = cavg_imgs_ref(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd        = smpd
        params%box         = ldim(1)
        params%msk         = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! extract nonzero cluster labels
        labels             = spproj_ref%os_cls2D%get_all_asint('cluster')
        labels             = pack(labels, mask=l_non_junk_ref)
        nclust             = maxval(labels)
        states             = spproj_ref%os_cls2D%get_all_asint('state')
        states             = pack(states, mask=l_non_junk_ref)
        ! calculate overall minmax
        oa_minmax(1)       = min(minval(mm_ref(:,1)),minval(mm_match(:,1)))
        oa_minmax(2)       = max(maxval(mm_ref(:,2)),maxval(mm_match(:,2)))
        oa_minmax_match(1) = minval(mm_match(:,1))
        oa_minmax_match(2) = maxval(mm_match(:,2))
        ! set folder & file names
        folder             = PATH_HERE
        chunk_fnames(1)    = params%projfile
        chunk_fnames(2)    = params%projfile_target
        cavgs_merged       = string('cavgs_merged')//params%ext
        ! set cluster_cavgs command line
        call cline_cluster_cavgs%set('mkdir',                   'no')
        call cline_cluster_cavgs%set('prune',           params%prune)
        call cline_cluster_cavgs%set('clust_crit', params%clust_crit)
        call cline_cluster_cavgs%set('hp',                 params%hp)
        call cline_cluster_cavgs%set('lp',                 params%lp)
        call cline_cluster_cavgs%set('mskdiam',       params%mskdiam)
        call cline_cluster_cavgs%set('nthr',             params%nthr)
        if( l_recluster )then
            ! merge spproj_ref & spproj_match projects
            call merge_chunk_projfiles(chunk_fnames, folder, spproj_merged, cavgs_out=cavgs_merged)
            call spproj_merged%write(params%projfile_merged)
            ! replace solution by reclustered solution
            call cline_cluster_cavgs%set('projfile', params%projfile_merged)
            call xcluster_cavgs%execute_safe(cline_cluster_cavgs)
        else
            ! identify medoids in reference set
            inds         = (/(i,i=1,nrefs)/)
            medoid_inds  = spproj_ref%os_cls2D%get_all_asint('medoid_ind') ! all medoid_inds
            medoid_inds  = pack(medoid_inds, mask=l_non_junk_ref)          ! all included medoid_inds after junk rejection
            medoid_map   = pack(inds,        mask=medoid_inds > 0)         ! mapping the medoids to physical class indices
            medoid_inds  = pack(medoid_inds, mask=medoid_inds > 0)         ! making medoid_inds congruent with medoid_map
            call unique(medoid_inds, medoid_inds_unique)
            if( minval(medoid_inds) < 1            ) THROW_HARD('medoid_inds does not meet expectation (<1)')
            if( maxval(medoid_inds) > nclust       ) THROW_HARD('medoid_inds does not meet expectation (>nclust)')
            if( size(medoid_inds_unique) /= nclust ) THROW_HARD('medoid_inds does not meeet expectation, but contains duplicates')
            if( allocated(medoid_inds_unique) ) deallocate(medoid_inds_unique)
            medoid_imgs_ref = extract_imgarr(cavg_imgs_ref, medoid_map)
            if( size(medoid_imgs_ref) /= nclust ) THROW_HARD('size(medoid_imgs_ref) /= nclust')
            ! cluster without any reference to the previous solution
            ! calculate distance matrix
            dmat = calc_cluster_cavgs_dmat(params, cavg_imgs_match, oa_minmax_match, params%clust_crit)
            ! cluster
            call cluster_dmat( dmat, 'aprop', nclust_match, i_medoids_match, labels_match, nclust_max=NCLUST_MAX)
            if( nclust_match < 3 )then
                nclust_match = 3
                call cluster_dmat(dmat, 'kmed', nclust_match, i_medoids_match, labels_match)
            endif
            ! extract medoid images for matching
            medoid_imgs_match = extract_imgarr(cavg_imgs_match, i_medoids_match)
            ! calculate distance matrix between medoids
            dmat_med = calc_match_cavgs_dmat(params, medoid_imgs_ref, medoid_imgs_match, oa_minmax, params%clust_crit)
            ! find closest matching reference medoids
            labels_med = minloc(dmat_med, dim=1)
            ! map new solution to reference solution
            allocate(labels_map(nmatch), source=0)
            do iclust_match = 1, nclust_match
                where( labels_match == iclust_match ) labels_map = medoid_inds(labels_med(iclust_match))
            end do
            labels_match = labels_map
            deallocate(labels_map)
            ! write class averages
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
            ! write project file
            call spproj_match%write(params%projfile_target)
            ! merge
            call merge_chunk_projfiles(chunk_fnames, folder, spproj_merged, cavgs_out=cavgs_merged)
            call spproj_merged%write(params%projfile_merged)
        endif
        ! cleanup
        call cline_cluster_cavgs%kill
        call spproj_match%kill
        call spproj_ref%kill
        call spproj_merged%kill
        call dealloc_imgarr(cavg_imgs_ref)
        call dealloc_imgarr(cavg_imgs_match)
        call dealloc_imgarr(medoid_imgs_ref)
        call dealloc_imgarr(medoid_imgs_match)
        ! end gracefully
        call simple_end('**** SIMPLE_MATCH_CAVGS NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')

    contains
        
        function find_label_state( label ) result( state )
            integer, intent(in) :: label
            integer :: i,state
            state = 1
            do i = 1, size(labels)
                if( labels(i) == label )then
                    state = states(i)
                    return
                endif
            end do
        end function find_label_state

    end subroutine exec_match_cavgs

    subroutine exec_map_cavgs_selection( self, cline )
        use simple_corrmat, only: calc_cartesian_corrmat
        class(commander_map_cavgs_selection), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(image),      allocatable :: imgs_sel(:), imgs_all(:)
        integer,          allocatable :: states(:)
        real,             allocatable :: correlations(:,:)
        type(string) :: cavgstk
        integer      :: iimg, isel, nall, nsel, loc(1), lfoo(3)
        real         :: smpd
        call cline%set('dir_exec', 'selection')
        call cline%set('mkdir',    'yes')
        if( .not.cline%defined('prune')   ) call cline%set('prune',   'no')
        if( .not.cline%defined('imgkind') ) call cline%set('imgkind', 'cavg')
        call build%init_params_and_build_spproj(cline,params)
        ! find number of selected cavgs
        call find_ldim_nptcls(params%stk2, lfoo, nsel)
        if( cline%defined('ares') ) nsel = int(params%ares)
        ! find number of original cavgs
        if( .not. cline%defined('stk' ) )then
            call build%spproj%get_cavgs_stk(cavgstk, nall, smpd, imgkind=params%imgkind)
            params%stk = cavgstk
        else
            call find_ldim_nptcls(params%stk, lfoo, nall)
        endif
        ! read images
        allocate(imgs_sel(nsel), imgs_all(nall))
        do isel=1,nsel
            call imgs_sel(isel)%new([params%box,params%box,1], params%smpd)
            call imgs_sel(isel)%read(params%stk2, isel)
        end do
        do iimg=1,nall
            call imgs_all(iimg)%new([params%box,params%box,1], params%smpd)
            call imgs_all(iimg)%read(params%stk, iimg)
        end do
        write(logfhandle,'(a)') '>>> CALCULATING CORRELATIONS'
        call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
        ! create the states array for mapping the selection
        allocate(states(nall), source=0)
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            states(loc(1)) = 1
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! optional pruning
        if( trim(params%prune).eq.'yes') call build%spproj%prune_particles
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_CAVGS_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_selection

    subroutine exec_map_cavgs_states( self, cline )
        use simple_corrmat, only: calc_cartesian_corrmat
        class(commander_map_cavgs_states), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline !< command line input
        type(parameters)          :: params
        type(builder)             :: build
        type(image),  allocatable :: imgs_sel(:), imgs_all(:)
        integer,      allocatable :: states(:)
        real,         allocatable :: correlations(:,:)
        type(string), allocatable :: stkfnames(:)
        type(string)              :: cavgstk, fname
        integer :: iimg, isel, nall, nsel, loc(1), lfoo(3), s
        real    :: smpd
        call cline%set('dir_exec', 'state_mapping')
        call cline%set('mkdir',    'yes')
        call build%init_params_and_build_spproj(cline,params)
        call read_filetable(params%stktab, stkfnames)
        ! find number of original cavgs
        if( .not. cline%defined('stk' ) )then
            call build%spproj%get_cavgs_stk(cavgstk, nall, smpd)
            params%stk = cavgstk
        else
            call find_ldim_nptcls(params%stk, lfoo, nall)
        endif
        ! read images
        allocate(imgs_all(nall))
        do iimg=1,nall
            call imgs_all(iimg)%new([params%box,params%box,1], params%smpd)
            call imgs_all(iimg)%read(params%stk, iimg)
        end do
        ! create the states array for mapping the selection
        allocate(states(nall), source=0)
        do s = 1,size(stkfnames)
            ! find number of selected cavgs
            fname = '../'//stkfnames(s)%to_char()
            call find_ldim_nptcls(fname, lfoo, nsel)
            ! read images
            allocate(imgs_sel(nsel))
            do isel=1,nsel
                call imgs_sel(isel)%new([params%box,params%box,1], params%smpd)
                call imgs_sel(isel)%read(fname, isel)
            end do
            call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
            do isel=1,nsel
                loc = maxloc(correlations(isel,:))
                states(loc(1)) = s
            end do
            ! destruct
            do isel=1,nsel
                call imgs_sel(isel)%kill
            end do
            deallocate(imgs_sel)
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_CAVGS_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_states

    subroutine exec_shape_rank_cavgs( self, cline )
        class(commander_shape_rank_cavgs), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        real,        parameter   :: LP_BIN = 20.
        type(parameters)         :: params
        type(sp_project)         :: spproj
        type(image), allocatable :: cavg_imgs(:), mask_imgs(:), masked_imgs(:)
        logical,     allocatable :: l_non_junk(:)
        real,        allocatable :: diams(:), shifts(:,:), ints(:)
        integer,     allocatable :: pops(:), order(:), clsinds(:), cavg_inds(:)
        integer :: ncls_sel, icls, ldim(3), i, irank, xtiles, ytiles
        real    :: mskrad
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        ! set defaults
        call set_automask2D_defaults(cline)
        ! parse parameters
        call params%new(cline)
        ! read project
        call spproj%read(params%projfile)
        ! extract cavgs from project
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        ! flag non-junk cavgs
        ldim   = cavg_imgs(1)%get_ldim()
        mskrad = real(ldim(1)/2) - COSMSKHALFWIDTH - 1.
        call flag_non_junk_cavgs(cavg_imgs, LP_BIN, mskrad, l_non_junk)
        clsinds = (/(icls,icls=1,size(cavg_imgs))/)
        clsinds = pack(clsinds, mask=l_non_junk)
        ! re-read non-junk cavg_imgs
        call dealloc_imgarr(cavg_imgs)
        cavg_imgs   = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        mask_imgs   = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        masked_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        ncls_sel  = size(cavg_imgs)
        ! extract class populations
        pops      = spproj%os_cls2D%get_all_asint('pop')
        pops      = pack(pops, mask=l_non_junk)
        ! Automasking
        call automask2D(mask_imgs, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        ! calc integrated intesities and shift
        allocate(ints(ncls_sel), source=0.)
        do icls = 1, ncls_sel
            call masked_imgs(icls)%mul(mask_imgs(icls))
            ints(icls) = masked_imgs(icls)%get_sum_int()
            call cavg_imgs(icls)%shift([shifts(icls,1),shifts(icls,2),0.])
        end do
        ! order
        order = (/(i,i=1,ncls_sel)/)
        call hpsort(order, p1_lt_p2 )
        call reverse(order) ! largest first
        ! communicate ranks to project file
        call spproj%os_cls2D%set_all2single('shape_rank', 0)
        do irank = 1, ncls_sel
            icls = clsinds(order(irank))
            call spproj%os_cls2D%set(icls, 'shape_rank', irank)
        end do
        ! write ranks to project file
        call spproj%write_segment_inside('cls2D')
        ! write class averages
        call write_imgarr(cavg_imgs, string(SHAPE_RANKED_CAVGS_MRCNAME), order)
        call spproj%shape_ranked_cavgs2jpg(cavg_inds, string(SHAPE_RANKED_CAVGS_JPGNAME), xtiles, ytiles)
        ! kill
        call spproj%kill
        call dealloc_imgarr(cavg_imgs)
        call dealloc_imgarr(mask_imgs)
        if( allocated(cavg_inds)  ) deallocate(cavg_inds)
        if( allocated(clsinds)    ) deallocate(clsinds)
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        if( allocated(diams)      ) deallocate(diams)
        if( allocated(shifts)     ) deallocate(shifts)
        if( allocated(pops)       ) deallocate(pops)
        call simple_end('**** SIMPLE_SHAPE_RANK_CAVGS NORMAL STOP ****')

        contains

            function p1_lt_p2( p1, p2 ) result( val )
                integer, intent(in) :: p1, p2
                logical :: val
                val = .false.
                if( ints(p1) < ints(p2) )then
                    val = .true.
                else if( diams(p1) < diams(p2) )then
                    val = .true.
                endif
            end function p1_lt_p2

    end subroutine exec_shape_rank_cavgs

end module simple_commanders_cavgs
