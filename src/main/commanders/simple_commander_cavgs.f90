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
        type(image),      allocatable :: cavg_imgs(:), cavg_imgs_sel(:)
        real,             allocatable :: frc(:), mm(:,:), jointscores(:), dmat(:,:), dmat_sel(:,:)
        real,             allocatable :: resvals(:), res_bad(:), res_good(:), res_maybe(:)
        logical,          allocatable :: l_non_junk(:), good_mask(:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), inds(:), inds_sel(:)
        integer,          allocatable :: labels_sel(:), clsinds_sel(:), i_medoids_sel(:)
        integer,          allocatable :: clspops(:), clspops_sel(:), states(:), labels4write(:)
        type(clust_info), allocatable :: clust_info_arr(:), clust_info_arr_sel(:), clust_info_arr_merged(:)
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(stats_struct)            :: res_stats
        integer                       :: ncls, ncls_sel, icls, cnt, rank, nptcls, nptcls_good, loc(1), ldim(3), cnt_clust
        integer                       :: i, j, ii, jj, nclust, iclust, pop_good, pop_bad, nclust_sel, ngood, minv_labels
        real                          :: fsc_res, rfoo, frac_good, best_res, worst_res, res_max
        real                          :: oa_min, oa_max, dist_rank, dist_rank_best, smpd, simsum
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        call cline%set('sh_inv',    'yes') ! shift invariant search
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',         'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',             10.)
        if( .not. cline%defined('kweight')    ) call cline%set('kweight',       'all')
        if( .not. cline%defined('lp')         ) call cline%set('lp',               6.)
        if( .not. cline%defined('prune')      ) call cline%set('prune',          'no')
        if( .not. cline%defined('clust_crit') ) call cline%set('clust_crit', 'hybrid')
        ! master parameters
        call params%new(cline)
        ! read project file
        call spproj%read(params%projfile)
        ncls           = spproj%os_cls2D%get_noris()
        ! prep class average stack
        call prep_cavgs4clustering(spproj, cavg_imgs, params%mskdiam, clspops, clsinds, l_non_junk, mm )
        ncls_sel       = size(cavg_imgs)
        smpd           = cavg_imgs(1)%get_smpd()
        ldim           = cavg_imgs(1)%get_ldim()
        ! ensure correct smpd/box in params class
        params%smpd    = smpd
        params%box     = ldim(1)
        params%msk     = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! calculate overall minmax
        oa_min         = minval(mm(:,1))
        oa_max         = maxval(mm(:,2))
        ! calculate distance matrrix
        dmat           = calc_cluster_cavgs_dmat(params, cavg_imgs, [oa_min,oa_max])
        ! cluster
        call cluster_cavg_imgs(dmat, nclust, i_medoids, labels)
        if( nclust <= 5 )then
            clust_info_arr_merged = align_and_score_cavg_clusters(params, dmat, cavg_imgs, clspops, i_medoids, labels, l_prelim=.false.)
        else
            clust_info_arr = align_and_score_cavg_clusters(params, dmat, cavg_imgs, clspops, i_medoids, labels, l_prelim=.true.)
            jointscores    = pack(clust_info_arr(:)%jointscore, mask=clust_info_arr(:)%good_bad == 1)
            if( all(jointscores >= SCORE_THRES_INCL) )then
                ! if all the selected class averages score that high, we are done
                clust_info_arr_merged = clust_info_arr
            else
                ! create good mask
                allocate(good_mask(ncls_sel), source=.false.)
                do icls = 1, ncls_sel
                    if( clust_info_arr(labels(icls))%good_bad == 1 )then
                        good_mask(icls) = .true.
                    endif
                end do
                ngood = count(good_mask)
                ! re-cluster
                allocate(inds(ncls_sel),source=(/(i,i=1,ncls_sel)/))
                inds_sel           = pack(inds,             mask=good_mask)
                cavg_imgs_sel      = pack_imgarr(cavg_imgs, mask=good_mask)
                clspops_sel        = pack(clspops,          mask=good_mask)
                dmat_sel           = calc_cluster_cavgs_dmat(params, cavg_imgs_sel, [oa_min,oa_max])
                call cluster_cavg_imgs(dmat_sel, nclust_sel, i_medoids_sel, labels_sel)

                print *, 'size(cavg_imgs_sel) ', size(cavg_imgs_sel)
                print *, 'size(labels_sel)    ', size(labels_sel)

                clust_info_arr_sel = align_and_score_cavg_clusters(params, dmat_sel, cavg_imgs_sel, clspops_sel, i_medoids_sel, labels_sel, l_prelim=.false.)
                ! index adjustment of deselected
                minv_labels = minval(labels, mask=.not.good_mask)
                where(     good_mask ) labels = 0
                where(.not.good_mask ) labels = labels + nclust_sel - minv_labels + 1
                ! merge labels
                do i = 1, ngood
                    if( labels(inds_sel(i)) == 0 )then
                        labels(inds_sel(i)) = labels_sel(i)
                    else
                        THROW_HARD('Index logic error')
                    endif
                end do
                allocate(clust_info_arr_merged(maxval(labels)))
                cnt_clust = 0
                do i = 1, nclust_sel
                    cnt_clust = cnt_clust + 1
                    clust_info_arr_merged(cnt_clust) = clust_info_arr_sel(i)
                    print *, 'taking index: '//int2str(i)//' from selected solution'
                end do
                do i = count(clust_info_arr(:)%good_bad == 1) + 1, nclust
                    cnt_clust = cnt_clust + 1
                    clust_info_arr_merged(cnt_clust) = clust_info_arr(i)
                    print *, 'taking index: '//int2str(i)//' from initial solution'
                end do
                print *, 'cnt_clust: ', cnt_clust, 'maxval(labels): ', maxval(labels)
            endif
        endif

        print *, labels

        nclust = maxval(labels)

        print *, 'nclust ', nclust

        select case(trim(params%clust_crit))
            case('powfm','fm','histfm','hybrid')
                ! re-create cavg_imgs
                call dealloc_imgarr(cavg_imgs)
                cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
                ! write ranked clusters
                call write_aligned_cavgs(labels, cavg_imgs, clust_info_arr_merged, 'cluster_ranked', trim(params%ext))
                ! update project
                call spproj%os_ptcl2D%transfer_class_assignment(spproj%os_ptcl3D)
                do iclust = 1, nclust
                    do icls = 1, ncls_sel 
                        if( labels(icls) == iclust )then
                            call spproj%os_cls2D%set(clsinds(icls),'cluster',iclust)                          ! 2D class field
                            call spproj%os_cls3D%set(clsinds(icls),'cluster',iclust)                          ! 3D class field
                            call spproj%os_ptcl2D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 2D particle field
                            call spproj%os_ptcl3D%set_field2single('class', clsinds(icls), 'cluster', iclust) ! 3D particle field
                        endif
                    enddo
                enddo
                ! report cluster info
                do iclust = 1, nclust
                    write(logfhandle,'(A,A,f5.1,A,f5.1,A,f5.1,A,f5.1,A,f5.1,A,I3)') 'cluster_ranked'//int2str_pad(iclust,2)//'.mrc',&
                    &' resolution(A) ',   clust_info_arr_merged(iclust)%res,& 
                    &' resscore(%) ',     clust_info_arr_merged(iclust)%resscore,& 
                    &' homogeneity(%) ',  clust_info_arr_merged(iclust)%homogeneity,&
                    &' clustscore(%) ',   clust_info_arr_merged(iclust)%clustscore,&
                    &' jointscore(%) ',   clust_info_arr_merged(iclust)%jointscore,&
                    &' good_bad_assign ', clust_info_arr_merged(iclust)%good_bad
                end do
                ! check number of particles selected
                nptcls      = sum(clust_info_arr_merged(:)%nptcls)
                nptcls_good = sum(clust_info_arr_merged(:)%nptcls, mask=clust_info_arr_merged(:)%good_bad == 1)
                frac_good   = real(nptcls_good)  / real(nptcls)
                write(logfhandle,'(a,1x,f8.2)') '% PARTICLES CLASSIFIED AS 1ST RATE: ', frac_good  * 100.
                ! calculate resolution statistics for good/bad classes
                res_good    = pack(clust_info_arr_merged(:)%res, mask=clust_info_arr_merged(:)%good_bad == 1)
                res_bad     = pack(clust_info_arr_merged(:)%res, mask=clust_info_arr_merged(:)%good_bad == 0)
                pop_good    = count(clust_info_arr_merged(:)%good_bad == 1)
                pop_bad     = count(clust_info_arr_merged(:)%good_bad == 0)
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
                    if( clust_info_arr_merged(labels(icls))%good_bad == 1 ) states(clsinds(icls)) = 1
                end do
                ! write selection
                allocate(labels4write(ncls_sel), source=0)
                do icls = 1, ncls_sel
                    labels4write(icls) = clust_info_arr_merged(labels(icls))%good_bad
                end do
                ! write selection
                call write_selected_cavgs(ncls_sel, cavg_imgs, labels4write, params%ext)
                ! map selection to project
                call spproj%map_cavgs_selection(states)
                ! optional pruning
                if( trim(params%prune).eq.'yes') call spproj%prune_particles
                ! this needs to be a full write as many segments are updated
                call spproj%write(params%projfile)
            case DEFAULT
                ! re-create cavg_imgs
                call dealloc_imgarr(cavg_imgs)
                cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
                call  write_cavgs(ncls_sel, cavg_imgs, labels, 'cluster', params%ext)
        end select
        ! destruct
        call spproj%kill
        ! call pows%kill
        do icls=1,ncls_sel
            call cavg_imgs(icls)%kill
        end do
        ! deallocate anything not specifically allocated above
        deallocate(cavg_imgs, l_non_junk, labels, clsinds, i_medoids)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****')


    contains

        subroutine cluster_cavg_imgs( dmat, nclust, i_medoids, labels )
            real,                 intent(in)    :: dmat(:,:)
            integer,              intent(inout) :: nclust
            integer, allocatable, intent(inout) :: i_medoids(:), labels(:)
            call cluster_dmat( dmat, 'aprop', nclust, i_medoids, labels, nclust_max=NCLUST_MAX)
            if( nclust > 5 .and. nclust < 20 )then
                nclust = 20
                deallocate(i_medoids, labels)
                call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
            endif
        end subroutine cluster_cavg_imgs

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
        integer,           allocatable :: clspops_ref(:), clsinds_ref(:), clspops_match(:), clsinds_match(:)
        logical,           allocatable :: l_non_junk_ref(:), l_non_junk_match(:)
        real,              allocatable :: mm_ref(:,:), mm_match(:,:)
        type(inpl_struct), allocatable :: algninfo(:,:)
        integer :: nmatch, nrefs, ldim(3), i
        real    :: smpd
        ! defaults
        call cline%set('oritype', 'cls2D')
        call cline%set('ctf',        'no')
        call cline%set('objfun',     'cc')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',       10.)
        if( .not. cline%defined('kweight') ) call cline%set('kweight', 'all')
        if( .not. cline%defined('lp')      ) call cline%set('lp',         6.)
        ! master parameters
        call params%new(cline)
        ! read base project file
        call spproj_ref%read(params%projfile)
        if( .not. spproj_ref%os_cls2D%isthere('cluster') ) THROW_HARD('Base project lacks clustering information in cls2D field')
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
        ! do the matching
        algninfo = match_imgs2refs(nrefs, nmatch, params%hp, params%lp, params%trs, cavg_imgs_ref, cavg_imgs_match)

    end subroutine exec_match_cavgs

end module simple_commander_cavgs
