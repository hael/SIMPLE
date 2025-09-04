module simple_commander_cavgs
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_builder,           only: builder, build_glob
use simple_cmdline,           only: cmdline
use simple_commander_base,    only: commander_base
use simple_parameters,        only: parameters, params_glob
use simple_sp_project,        only: sp_project
use simple_image,             only: image
use simple_stack_io,          only: stack_io
use simple_strategy2D_utils
implicit none

public :: rank_cavgs_commander
public :: cluster_cavgs_commander
public :: select_clusters_commander
public :: match_cavgs_commander
public :: match_cavgs2afm_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: rank_cavgs_commander
    contains
        procedure :: execute      => exec_rank_cavgs
end type rank_cavgs_commander

type, extends(commander_base) :: cluster_cavgs_commander
    contains
        procedure :: execute      => exec_cluster_cavgs
end type cluster_cavgs_commander

type, extends(commander_base) :: select_clusters_commander
    contains
        procedure :: execute      => exec_select_clusters
end type select_clusters_commander

type, extends(commander_base) :: match_cavgs_commander
    contains
        procedure :: execute      => exec_match_cavgs
end type match_cavgs_commander

type, extends(commander_base) :: match_cavgs2afm_commander
    contains
        procedure :: execute      => exec_match_cavgs2afm
end type match_cavgs2afm_commander

contains

    subroutine exec_rank_cavgs( self, cline )
        class(rank_cavgs_commander), intent(inout) :: self
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
            call clsdoc_ranked%write('classdoc_ranked.txt')
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
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        logical,          parameter   :: DEBUG = .true.
        real,             parameter   :: SCORE_THRES_INCL = 75.
        integer,          parameter   :: NCLUST_MAX = 65
        type(image),      allocatable :: cavg_imgs(:), cluster_imgs(:)
        real,             allocatable :: frc(:), mm(:,:), jointscores(:), dmat(:,:), dmat_sel(:,:)
        real,             allocatable :: resvals(:), res_bad(:), res_good(:), res_maybe(:), clustscores(:)
        logical,          allocatable :: l_non_junk(:), good_mask(:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), inds(:), cluster_inds(:)
        integer,          allocatable :: clspops(:), clspops_sel(:), states(:), labels4write(:)
        type(clust_info), allocatable :: clust_info_arr(:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(stats_struct)            :: res_stats
        integer                       :: ncls, ncls_sel, icls, cnt, rank, nptcls, nptcls_good, loc(1)
        integer                       :: i, j, ii, jj, nclust, iclust, pop_good, pop_bad, nclust_sel
        integer                       :: ngood, minv_labels, ind, ldim(3), cnt_clust, pop
        real                          :: fsc_res, rfoo, frac_good, best_res, worst_res, res_max
        real                          :: oa_min, oa_max, dist_rank, dist_rank_best, smpd, simsum
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',       10.)
        if( .not. cline%defined('kweight')    ) call cline%set('kweight', 'all')
        if( .not. cline%defined('lp')         ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune')      ) call cline%set('prune',    'no')
        ! master parameters
        call params%new(cline)
        ! read project file
        call spproj%read(params%projfile)
        ncls        = spproj%os_cls2D%get_noris()
        ! prep class average stack
        call prep_cavgs4clustering(spproj, cavg_imgs, params%mskdiam, clspops, clsinds, l_non_junk, mm )
        ncls_sel    = size(cavg_imgs)
        smpd        = cavg_imgs(1)%get_smpd()
        ldim        = cavg_imgs(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! calculate overall minmax
        oa_min      = minval(mm(:,1))
        oa_max      = maxval(mm(:,2))
        if( trim(params%have_clustering).eq.'yes' )then
            labels  = spproj%os_cls2D%get_all_asint('cluster')
            labels  = pack(labels, mask=l_non_junk)
            if( all(labels == 0) ) THROW_HARD('Invalid clustering solution in cls2D field')
            nclust  = maxval(labels)
            allocate(i_medoids(nclust),   source=0)
            allocate(clustscores(nclust), source=0.)
            allocate(inds(ncls_sel), source=(/(i,i=1,ncls_sel)/))
            do iclust = 1, nclust
                pop = count(labels == iclust)
                if( pop > 0 )then
                    cluster_inds = pack(inds, mask=labels == iclust)
                    cluster_imgs = pack_imgarr(cavg_imgs, mask=labels == iclust)
                    dmat         = calc_cluster_cavgs_dmat(params, cluster_imgs, [oa_min,oa_max])
                    call medoid_from_dmat(dmat, ind)
                    clustscores(iclust) = sum(dmat(ind,:)) / real(pop)
                    i_medoids(iclust) = cluster_inds(ind)
                    deallocate(cluster_inds, dmat)
                    call dealloc_imgarr(cluster_imgs)
                endif
            end do
            clust_info_arr = align_and_score_cavg_clusters( params, dmat, cavg_imgs, clspops, i_medoids, labels, clustscores )
        else
            ! calculate distance matrix
            dmat = calc_cluster_cavgs_dmat(params, cavg_imgs, [oa_min,oa_max])
            ! cluster
            call cluster_dmat( dmat, 'aprop', nclust, i_medoids, labels, nclust_max=NCLUST_MAX)
            if( cline%defined('ncls') )then
                if( nclust > params%ncls )then
                    nclust = params%ncls
                    deallocate(i_medoids, labels)
                    call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
                endif
            else
                if( nclust > 5 .and. nclust < 20 )then
                    nclust = 20
                    deallocate(i_medoids, labels)
                    call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
                endif
            endif
            clust_info_arr = align_and_score_cavg_clusters( params, dmat, cavg_imgs, clspops, i_medoids, labels )
        endif
        ! re-create cavg_imgs
        call dealloc_imgarr(cavg_imgs)
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
        ! write ranked clusters
        call write_aligned_cavgs(labels, cavg_imgs, clust_info_arr, 'cluster_ranked', trim(params%ext))
        ! update project
        call spproj%os_ptcl2D%transfer_class_assignment(spproj%os_ptcl3D)
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
        ! calculate resolution statistics for good/bad classes
        res_good    = pack(clust_info_arr(:)%res, mask=clust_info_arr(:)%good_bad == 1)
        res_bad     = pack(clust_info_arr(:)%res, mask=clust_info_arr(:)%good_bad == 0)
        pop_good    = count(clust_info_arr(:)%good_bad == 1)
        pop_bad     = count(clust_info_arr(:)%good_bad == 0)
        if( pop_good > 1 )then
            write(logfhandle,'(A)') 'RESOLUTION STATS FOR GOOD PARTITION'
            call calc_stats(res_good, res_stats)
            write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_stats%minv
            write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_stats%maxv
            write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_stats%avg
            write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_stats%med
            write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', res_stats%sdev
        else if( pop_good == 1 )then
            write(logfhandle,'(A)') 'RESOLUTION STATS FOR GOOD PARTITION'
            write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_good(1)
            write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_good(1)
            write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_good(1)
            write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_good(1)
            write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', 0.
        endif
        if( pop_bad > 1 )then
            write(logfhandle,'(A)') 'RESOLUTION STATS FOR BAD  PARTITION'
            call calc_stats(res_bad, res_stats)
            write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_stats%minv
            write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_stats%maxv
            write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_stats%avg
            write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_stats%med
            write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', res_stats%sdev
        else if( pop_bad == 1 )then
            write(logfhandle,'(A)') 'RESOLUTION STATS FOR BAD  PARTITION'
            write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_bad(1)
            write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_bad(1)
            write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_bad(1)
            write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_bad(1)
            write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', 0.
        endif
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
        call write_selected_cavgs(ncls_sel, cavg_imgs, labels4write, params%ext)
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
        class(select_clusters_commander), intent(inout) :: self
        class(cmdline),                   intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        integer, allocatable :: clustinds(:)
        integer :: iclust, nclust_sel, nclust_max
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('prune') ) call cline%set('prune', 'yes')
        ! master parameters
        call params%new(cline)
        if( cline%defined('clustind') )then
            nclust_sel = 1
            allocate(clustinds(nclust_sel), source=params%clustind)
        else
            clustinds = listofints2arr(params%clustinds)
            nclust_sel = size(clustinds)
        endif
        do iclust = 1, nclust_sel
            print *, 'selected cluster: ', clustinds(iclust)
        end do
        ! read project file
        call spproj%read(params%projfile)
        nclust_max = spproj%os_ptcl2D%get_n('cluster')
        if( any(clustinds > nclust_max) ) THROW_HARD('Maximum cluster index value: '//int2str(nclust_max)//' exceeded!')
        do iclust = 1, nclust_max
            if( any(clustinds == iclust) )then
                ! set state=1 to flag inclusion
                call spproj%os_cls2D%set_field2single('cluster',  iclust, 'state', 1) ! 2D class field
                call spproj%os_cls3D%set_field2single('cluster',  iclust, 'state', 1) ! 3D class field
                call spproj%os_ptcl2D%set_field2single('cluster', iclust, 'state', 1) ! 2D particle field
                call spproj%os_ptcl3D%set_field2single('cluster', iclust, 'state', 1) ! 3D particle field
            else
                ! set state=0 to flag exclusion
                call spproj%os_cls2D%set_field2single('cluster',  iclust, 'state', 0) ! 2D class field
                call spproj%os_cls3D%set_field2single('cluster',  iclust, 'state', 0) ! 3D class field
                call spproj%os_ptcl2D%set_field2single('cluster', iclust, 'state', 0) ! 2D particle field
                call spproj%os_ptcl3D%set_field2single('cluster', iclust, 'state', 0) ! 3D particle field
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
        class(match_cavgs_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj_ref, spproj_match
        type(image),       allocatable :: cavg_imgs_ref(:), cavg_imgs_match(:)
        integer,           allocatable :: clspops_ref(:), clsinds_ref(:), clspops_match(:), clsinds_match(:), labels(:)
        integer,           allocatable :: i_medoids_ref(:), states(:), labels_match(:), states_map(:)
        logical,           allocatable :: l_non_junk_ref(:), l_non_junk_match(:), l_med_msk(:)
        real,              allocatable :: mm_ref(:,:), mm_match(:,:), corrmat(:,:), dmat_clust(:,:), dmat(:,:)
        integer :: nmatch, nrefs, ldim(3), i, j, ncls_match, nclust, icls, iclust, imatch
        real    :: smpd, oa_minmax(2)
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',       10.)
        if( .not. cline%defined('kweight')    ) call cline%set('kweight', 'all')
        if( .not. cline%defined('lp')         ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune')      ) call cline%set('prune',    'no')
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
        nrefs       = size(cavg_imgs_ref)
        nmatch      = size(cavg_imgs_match)
        smpd        = cavg_imgs_ref(1)%get_smpd()
        ldim        = cavg_imgs_ref(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! extract nonzero cluster labels
        labels = spproj_ref%os_cls2D%get_all_asint('cluster')
        labels = pack(labels, mask=l_non_junk_ref)
        nclust = maxval(labels)
        states = spproj_ref%os_cls2D%get_all_asint('state')
        states = pack(states, mask=l_non_junk_ref)
        ! calculate overall minmax
        oa_minmax(1) = minval(mm_ref(:,1))
        oa_minmax(2) = maxval(mm_ref(:,2))
        ! generate matching distance matrix
        dmat = calc_match_cavgs_dmat(params, cavg_imgs_ref, cavg_imgs_match, oa_minmax)
        ! genrate cluster distance matrix
        allocate(dmat_clust(nclust,nmatch), source=0.)
        do iclust = 1, nclust
            do imatch = 1, nmatch
                dmat_clust(iclust,imatch) = sum(dmat(:,imatch), mask=labels == iclust) / real(count(labels == iclust))
            end do
        end do
        allocate(labels_match(nmatch), source=0)
        labels_match = minloc(dmat_clust, dim=1)
        call  write_cavgs(nmatch, cavg_imgs_match, labels_match, 'cluster_match', params%ext)
        ! update project
        call spproj_match%os_ptcl2D%transfer_class_assignment(spproj_match%os_ptcl3D)
        do iclust = 1, nclust
            do icls = 1, nmatch
                if( labels_match(icls) == iclust )then
                    call spproj_match%os_cls2D%set(clsinds_match(icls),'cluster',iclust)                          ! 2D class field
                    call spproj_match%os_cls3D%set(clsinds_match(icls),'cluster',iclust)                          ! 3D class field
                    call spproj_match%os_ptcl2D%set_field2single('class', clsinds_match(icls), 'cluster', iclust) ! 2D particle field
                    call spproj_match%os_ptcl3D%set_field2single('class', clsinds_match(icls), 'cluster', iclust) ! 3D particle field
                endif
            enddo
        enddo
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
        ! cleanup
        call spproj_match%kill
        call spproj_ref%kill
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

    subroutine exec_match_cavgs2afm( self, cline )
        class(match_cavgs2afm_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj_ref
        type(stack_io)   :: stkio_r
        type(image)      :: img_msk   
        type(image),       allocatable :: cavg_imgs_ref(:), cavg_imgs_rescale(:), afm_imgs(:)
        logical,           allocatable :: l_msk(:,:,:)
        integer,           allocatable :: clspops_ref(:), clsinds_ref(:), clspops_match(:), clsinds_match(:), cavg_labels(:), afm_labels(:)
        integer,           allocatable :: i_medoids_ref(:), states(:), labels_match(:), states_map(:)
        real,              allocatable :: mm_ref(:,:), mm_afm(:,:), corrmat(:,:), dmat_clust(:,:), dmat(:,:), rank_cavgs(:), rank_afm(:)
        integer :: nmatch, nrefs, ldim(3), ldim_afm(3), i, j, ncls_match, nclust, nafm, icls, iclust, imatch, box_new, dmat_shape(2)
        real    :: smpd, smpd_target, oa_minmax(2), mskrad
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',       10.)
        if( .not. cline%defined('kweight')    ) call cline%set('kweight', 'all')
        if( .not. cline%defined('lp')         ) call cline%set('lp',         6.)
        if( .not. cline%defined('prune')      ) call cline%set('prune',    'no')
        ! master parameters
        call params%new(cline)
        ! read base project file, write to image array
        call spproj_ref%read(params%projfile)
        if( .not. spproj_ref%os_cls2D%isthere('cluster') ) THROW_HARD('Reference project lacks clustering information in cls2D field')
        cavg_imgs_ref = read_cavgs_into_imgarr(spproj_ref)
        ! read afm stack
        call find_ldim_nptcls(params%stk, ldim_afm, nafm)
        call stkio_r%open(params%stk, params%smpd_target, 'read', box=ldim_afm(1), bufsz=nafm)
        call stkio_r%read_whole 
        allocate(afm_imgs(nafm))
        do i = 1, nafm
            call afm_imgs(i)%new([ldim_afm(1), ldim_afm(2), 1],params%smpd_target,.false.)
            call stkio_r%read(i, afm_imgs(i))
        end do 
        nrefs       = size(cavg_imgs_ref)
        smpd        = cavg_imgs_ref(1)%get_smpd()
        ldim        = cavg_imgs_ref(1)%get_ldim()
        allocate(cavg_imgs_rescale(nrefs))
        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        box_new        = ldim_afm(1)
        params%mskdiam = nint(0.9 * real(box_new))
        params%msk     = params%mskdiam / 2.
        call img_msk%new([box_new,box_new,1], smpd_target)
        img_msk = 1.
        call img_msk%mask(params%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        allocate(mm_afm(nafm,2))
        allocate(mm_ref(nrefs,2))
        ! downscale and prepare cavgs 
        do i = 1, nrefs 
            call cavg_imgs_rescale(i)%new([box_new, box_new, 1], params%smpd_target)
            call cavg_imgs_ref(i)%fft
            call cavg_imgs_rescale(i)%fft
            call cavg_imgs_ref(i)%lp(nint(params%lp))
            call cavg_imgs_ref(i)%clip(cavg_imgs_rescale(i))
            call cavg_imgs_rescale(i)%set_smpd(params%smpd_target)
            call cavg_imgs_rescale(i)%ifft
            call cavg_imgs_rescale(i)%norm_within(l_msk)
            call cavg_imgs_rescale(i)%mask(params%msk, 'soft', backgr=0.)
            mm_ref(i,:) = cavg_imgs_ref(i)%minmax(params%msk)
        end do 
        ! prepare afm 
        do i = 1, nafm 
            call afm_imgs(i)%fft
            call afm_imgs(i)%lp(nint(params%lp))
            call afm_imgs(i)%ifft
            call afm_imgs(i)%norm_within(l_msk)
            call afm_imgs(i)%mask(params%msk, 'soft', backgr=0.)
            mm_afm(i,:) = afm_imgs(i)%minmax(params%msk)
        end do 
        oa_minmax(1) = (minval(mm_ref(:,1)) + minval(mm_afm(:,1))) / 2.
        oa_minmax(2) = (maxval(mm_ref(:,2)) + minval(mm_afm(:,2))) / 2.
        ! generate matching distance matrix
        dmat = calc_match_cavgs_dmat(params, afm_imgs, cavg_imgs_rescale, oa_minmax)
        allocate(rank_cavgs(nrefs))
        allocate(cavg_labels(nrefs))
        do i = 1, nrefs 
            rank_cavgs(i)  = sum(dmat(:,i)) / real(nafm)
            cavg_labels(i) = i
        end do 
        ! rank afm 
        allocate(rank_afm(nafm))
        allocate(afm_labels(nafm))
        do i = 1, nafm 
            rank_afm(i)     = sum(dmat(i,:)) / real(nrefs)
            afm_labels(i)   = i
        end do 
        ! sort + write 
        call hpsort(rank_cavgs, cavg_labels)
        call hpsort(rank_afm, afm_labels)
        do i = 1, nrefs 
            call cavg_imgs_rescale(cavg_labels(i))%write('ranked_cavg.mrc', i)
        end do 
        do i = 1, nafm 
            call afm_imgs(afm_labels(i))%write('ranked_afm.mrc', i)
        end do 
        call dealloc_imgarr(cavg_imgs_ref)
        call dealloc_imgarr(afm_imgs)
        call dealloc_imgarr(cavg_imgs_rescale)
        ! end gracefully
        call simple_end('**** SIMPLE_match_cavgs2afmS NORMAL STOP ****', verbose_exit=trim(params%verbose_exit).eq.'yes')
    end subroutine exec_match_cavgs2afm

end module simple_commander_cavgs
