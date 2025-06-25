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
use simple_pspecs,            only: pspecs
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
        use simple_corrmat,          only: calc_inpl_invariant_fm
        use simple_histogram,        only: histogram
        use simple_clustering_utils, only: cluster_dmat
        class(cluster_cavgs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        real,             parameter   :: HP_SPEC = 20., LP_SPEC = 6., RES_THRES=6., SCORE_THRES=65., SCORE_THRES_REJECT=50., SCORE_THRES_INCL=75.
        integer,          parameter   :: NHISTBINS = 128, NQUANTA=32
        logical,          parameter   :: DEBUG = .true.
        type(image),      allocatable :: cavg_imgs(:)
        type(histogram),  allocatable :: hists(:)
        real,             allocatable :: frc(:), mm(:,:)
        real,             allocatable :: corrmat(:,:), dmat_pow(:,:), smat_pow(:,:), dmat_tvd(:,:), smat_tvd(:,:), dmat_joint(:,:)
        real,             allocatable :: smat_joint(:,:), dmat(:,:), res_bad(:), res_good(:), res_maybe(:), dmat_jsd(:,:), smat_jsd(:,:)
        real,             allocatable :: dmat_hd(:,:), dmat_hist(:,:), dmat_fm(:,:), smat(:,:)
        real,             allocatable :: resvals(:)
        logical,          allocatable :: l_msk(:,:,:), l_non_junk(:)
        integer,          allocatable :: labels(:), clsinds(:), i_medoids(:), labels_copy(:), i_medoids_copy(:)
        integer,          allocatable :: clspops(:), states(:), labels4write(:)
        type(clust_info), allocatable :: clust_info_arr(:), clust_info_arr_copy(:)
        type(parameters)   :: params
        type(sp_project)   :: spproj
        type(image)        :: img_msk
        type(pspecs)       :: pows
        type(stats_struct) :: res_stats
        integer            :: ncls, ncls_sel, icls, cnt, rank, nptcls, nptcls_good, loc(1), ldim(3)
        integer            :: i, j, ii, jj, nclust, iclust, nptcls_maybe, pop_good, pop_bad, pop_maybe
        real               :: fsc_res, rfoo, frac_good, best_res, worst_res, frac_maybe
        real               :: oa_min, oa_max, dist_rank, dist_rank_best, smpd, simsum
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
        ! generate histograms
        allocate(hists(ncls_sel))
        do i = 1, ncls_sel
            call hists(i)%new(cavg_imgs(i), NHISTBINS, minmax=[oa_min,oa_max], radius=params%msk)
        end do
        ! generate power spectra and associated distance/similarity matrix
        ! create pspecs object
        call pows%new(cavg_imgs, spproj%os_cls2D, params%msk, HP_SPEC, LP_SPEC, params%ncls_spec, l_exclude_junk=.false.)
        ! create a joint similarity matrix for clustering based on spectral profile and in-plane invariant correlation
        call pows%calc_distmat(is_l1=trim(params%dist_type).eq.'l1')
        dmat_pow = pows%get_distmat()
        call normalize_minmax(dmat_pow)
        smat_pow = dmat2smat(dmat_pow)
        ! calculate inpl_invariant_fm corrmat
        ! pairwise correlation through Fourier-Mellin + shift search
        write(logfhandle,'(A)') '>>> PAIRWISE CORRELATIONS THROUGH FOURIER-MELLIN & SHIFT SEARCH'
        call calc_inpl_invariant_fm(cavg_imgs, params%hp, params%lp, params%trs, corrmat)
        dmat_fm = smat2dmat(corrmat)
        ! calculate histogram-based distance matrices
        allocate(dmat_tvd(ncls_sel,ncls_sel), dmat_jsd(ncls_sel,ncls_sel), dmat_hd(ncls_sel,ncls_sel), source=0.)
        !$omp parallel do default(shared) private(i,j)&
        !$omp schedule(dynamic) proc_bind(close)
        do i = 1, ncls_sel - 1
            do j = i + 1, ncls_sel
                dmat_tvd(i,j) = hists(i)%TVD(hists(j))
                dmat_tvd(j,i) = dmat_tvd(i,j)
                dmat_jsd(i,j) = hists(i)%JSD(hists(j))
                dmat_jsd(j,i) = dmat_jsd(i,j)
                dmat_hd(i,j)  = hists(i)%HD(hists(j))
                dmat_hd(j,i)  = dmat_hd(i,j)
            end do
        end do
        !$omp end parallel do
        call hists(:)%kill
        dmat_hist = merge_dmats(dmat_tvd, dmat_jsd, dmat_hd) ! the different histogram distances are given equal weight
        ! set appropriate distance matrix for the clustering criterion given
        select case(trim(params%clust_crit))
            case('pow')
                dmat = dmat_pow 
            case('powfm')
                dmat = 0.5 * dmat_pow + 0.5 * dmat_fm
            case('fm')
                dmat = dmat_fm
            case('hist')
                dmat = dmat_hist
            case('histfm')
                dmat = 0.5 * dmat_hist + 0.5 * dmat_fm
            case('hybrid')
                dmat = 0.2 * dmat_hist + 0.4 * dmat_pow + 0.4 * dmat_fm       
            case DEFAULT
                THROW_HARD('Unsupported clustering criterion: '//trim(params%clust_crit))
        end select
        ! cluster
        call cluster_dmat( dmat, 'aprop', nclust, i_medoids, labels)
        if( nclust > 5 .and. nclust < 20 )then
            nclust = 20
            deallocate(i_medoids, labels)
            call cluster_dmat(dmat, 'kmed', nclust, i_medoids, labels)
        endif
        ! prep mask
        call img_msk%new([params%box,params%box,1], params%smpd)
        img_msk = 1.
        call img_msk%mask(params%msk, 'hard')
        l_msk = img_msk%bin2logical()
        call img_msk%kill
        select case(trim(params%clust_crit))
            case('powfm','fm','histfm','hybrid')
                ! align clusters to medoids and gather information
                clust_info_arr = align_clusters2medoids( labels, i_medoids, cavg_imgs, params%hp, params%lp, params%trs, l_msk )
                ! set particle populations
                do iclust = 1, nclust
                    clust_info_arr(iclust)%nptcls = sum(clspops, mask=labels == iclust)
                end do
                ! calculate scores
                call calc_scores
                ! rank clusters
                call rank_clusters
                ! calculate discretized joint score based on a second pass of AP clustering of score vecs
                call calc_jointscore
                ! re-rank clusters
                call rank_clusters
                ! good/bad cluster identification through tresholding
                call identify_good_bad_clusters
                write(logfhandle,'(A)') '>>> ROTATING & SHIFTING UNMASKED, UNFILTERED CLASS AVERAGES'
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
                nptcls       = sum(clust_info_arr(:)%nptcls)
                nptcls_good  = sum(clust_info_arr(:)%nptcls, mask=clust_info_arr(:)%good_bad == 1)
                nptcls_maybe = sum(clust_info_arr(:)%nptcls, mask=clust_info_arr(:)%good_bad == 2)
                frac_good    = real(nptcls_good)  / real(nptcls)
                frac_maybe   = real(nptcls_maybe) / real(nptcls)
                write(logfhandle,'(a,1x,f8.2)') '% PARTICLES CLASSIFIED AS 1ST RATE: ', frac_good  * 100.
                write(logfhandle,'(a,1x,f8.2)') '% PARTICLES CLASSIFIED AS 2ND RATE: ', frac_maybe * 100.
                ! calculate resolution statistics for good/maybe/bad classes
                res_good     = pack(clust_info_arr(:)%res, mask=clust_info_arr(:)%good_bad == 1)
                res_maybe    = pack(clust_info_arr(:)%res, mask=clust_info_arr(:)%good_bad == 2)
                res_bad      = pack(clust_info_arr(:)%res, mask=clust_info_arr(:)%good_bad == 0)
                pop_good     = count(clust_info_arr(:)%good_bad == 1)
                pop_maybe    = count(clust_info_arr(:)%good_bad == 2)
                pop_bad      = count(clust_info_arr(:)%good_bad == 0)
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
                if( pop_maybe > 1 )then
                    write(logfhandle,'(A)') 'RESOLUTION STATS FOR MAYBE PARTITION'
                    call calc_stats(res_maybe, res_stats)
                    write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_stats%minv
                    write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_stats%maxv
                    write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_stats%avg
                    write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_stats%med
                    write(logfhandle,'(a,1x,f8.2)') 'SDEV    RES: ', res_stats%sdev
                else if( pop_maybe == 1 )then
                    write(logfhandle,'(A)') 'RESOLUTION STATS FOR MAYBE PARTITION'
                    write(logfhandle,'(a,1x,f8.2)') 'MINIMUM RES: ', res_maybe(1)
                    write(logfhandle,'(a,1x,f8.2)') 'MAXIMUM RES: ', res_maybe(1)
                    write(logfhandle,'(a,1x,f8.2)') 'AVERAGE RES: ', res_maybe(1)
                    write(logfhandle,'(a,1x,f8.2)') 'MEDIAN  RES: ', res_maybe(1)
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
                    if( clust_info_arr(labels(icls))%good_bad == 2 ) states(clsinds(icls)) = 2
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
            case DEFAULT
                ! re-create cavg_imgs
                call dealloc_imgarr(cavg_imgs)
                cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_non_junk)
                call  write_cavgs(ncls_sel, cavg_imgs, labels, 'cluster', params%ext)
        end select
        ! destruct
        call spproj%kill
        call pows%kill
        do icls=1,ncls_sel
            call cavg_imgs(icls)%kill
        end do
        ! deallocate anything not specifically allocated above
        deallocate(cavg_imgs, dmat, l_msk, l_non_junk, labels, clsinds, i_medoids)
        if( allocated(corrmat)    ) deallocate(corrmat)
        if( allocated(dmat_pow)   ) deallocate(dmat_pow)
        if( allocated(smat_pow)   ) deallocate(smat_pow)
        if( allocated(smat_joint) ) deallocate(smat_joint)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER_CAVGS NORMAL STOP ****')

    contains

        subroutine calc_scores
            integer :: iclust, icls, cnt
            real    :: euclid_max, res_max, corrfmscore_min, clustscore_min
            ! HOMOGENEITY SCORE
            clust_info_arr(:)%homogeneity = clust_info_arr(:)%euclid
            euclid_max = maxval(clust_info_arr(:)%euclid)
            where( clust_info_arr(:)%pop < 2 ) clust_info_arr(:)%homogeneity = euclid_max
            call dists2scores_percen(clust_info_arr(:)%homogeneity)
            ! RESOLUTION SCORE
            clust_info_arr(:)%resscore = clust_info_arr(:)%res
            res_max = maxval(clust_info_arr(:)%res)
            where( clust_info_arr(:)%pop < 2 ) clust_info_arr(:)%resscore = res_max
            call dists2scores_percen(clust_info_arr(:)%resscore)
            ! FM CORR, CLUSTSCORE
            do iclust = 1, nclust
                clust_info_arr(iclust)%corrfmscore = 0.
                clust_info_arr(iclust)%clustscore  = 0.
                cnt = 0
                do icls = 1, ncls_sel 
                    if( labels(icls) == iclust )then
                        clust_info_arr(iclust)%corrfmscore = clust_info_arr(iclust)%corrfmscore +   corrmat(icls,i_medoids(iclust))
                        clust_info_arr(iclust)%clustscore  = clust_info_arr(iclust)%clustscore  +      dmat(icls,i_medoids(iclust))
                        cnt = cnt + 1
                    endif
                enddo
                clust_info_arr(iclust)%corrfmscore = clust_info_arr(iclust)%corrfmscore / real(cnt)
                clust_info_arr(iclust)%clustscore  = clust_info_arr(iclust)%clustscore  / real(cnt)
            end do
            corrfmscore_min = minval(clust_info_arr(:)%corrfmscore)
            clustscore_min  = minval(clust_info_arr(:)%clustscore)
            where( clust_info_arr(:)%pop < 2 )
                clust_info_arr(:)%corrfmscore = corrfmscore_min
                clust_info_arr(:)%clustscore  = clustscore_min
            endwhere
            call scores2scores_percen(clust_info_arr(:)%corrfmscore)
            call dists2scores_percen(clust_info_arr(:)%clustscore)
            ! calculate joint score
            clust_info_arr(:)%jointscore = 0.35 * clust_info_arr(:)%homogeneity + 0.5 * clust_info_arr(:)%resscore + 0.15 * clust_info_arr(:)%clustscore
            call scores2scores_percen(clust_info_arr(:)%jointscore)
        end subroutine calc_scores

        subroutine copy_clustering
            clust_info_arr_copy = clust_info_arr
            labels_copy         = labels
            i_medoids_copy      = i_medoids
        end subroutine copy_clustering

        subroutine put_back_clustering_copy
            clust_info_arr = clust_info_arr_copy 
            labels         = labels_copy
            i_medoids      = i_medoids_copy 
        end subroutine put_back_clustering_copy

        subroutine rank_clusters
            integer              :: i_medoids_ranked(nclust), rank_assign(ncls_sel)
            type(clust_info)     :: clust_info_arr_ranked(nclust)
            integer              :: rank, icls
            integer, allocatable :: clust_order(:)
            clust_order = scores2order(clust_info_arr(:)%jointscore)
            do rank = 1, nclust
                i_medoids_ranked(rank) = i_medoids(clust_order(rank))
                clust_info_arr_ranked(rank) = clust_info_arr(clust_order(rank))
                do icls = 1, ncls_sel
                    if( labels(icls) == clust_order(rank) ) rank_assign(icls) = rank
                end do
            end do
            i_medoids      = i_medoids_ranked
            clust_info_arr = clust_info_arr_ranked
            labels         = rank_assign
            deallocate(clust_order)
        end subroutine rank_clusters

        subroutine identify_good_bad_clusters
            real, allocatable :: jointscores_2(:), resvals_2(:), homogeneity_2(:), clustscores_2(:)
            integer           :: scoreclust_1, scoreclust_2
            real              :: avgscore_2
            logical           :: l_incl_rank2
            clust_info_arr(:)%good_bad = 0
            if( nclust <= 3 )then
                clust_info_arr(:)%good_bad      = 1
                clust_info_arr(nclust)%good_bad = 0
            else
                scoreclust_1  = clust_info_arr(1)%scoreclust
                where( clust_info_arr(:)%scoreclust == scoreclust_1 ) clust_info_arr(:)%good_bad = 1
                scoreclust_2  = clust_info_arr(count(clust_info_arr(:)%good_bad == 1) + 1)%scoreclust
                jointscores_2 = pack(clust_info_arr(:)%jointscore,  mask=clust_info_arr(:)%scoreclust == scoreclust_2)
                homogeneity_2 = pack(clust_info_arr(:)%homogeneity, mask=clust_info_arr(:)%scoreclust == scoreclust_2)
                clustscores_2 = pack(clust_info_arr(:)%clustscore,  mask=clust_info_arr(:)%scoreclust == scoreclust_2)
                resvals_2     = pack(clust_info_arr(:)%res,         mask=clust_info_arr(:)%scoreclust == scoreclust_2)
                avgscore_2    = sum(jointscores_2) / real(count(clust_info_arr(:)%scoreclust == scoreclust_2))
                l_incl_rank2  = .false. 
                if( any(resvals_2 <= RES_THRES) )            l_incl_rank2 = .true.
                if( avgscore_2 >= SCORE_THRES   )            l_incl_rank2 = .true.
                if( any(homogeneity_2 >= SCORE_THRES_INCL) .and.&
                    any(clustscores_2 >= SCORE_THRES_INCL) ) l_incl_rank2 = .true.
                if( l_incl_rank2 )then
                    where( clust_info_arr(:)%scoreclust == scoreclust_2 ) clust_info_arr(:)%good_bad = 1
                else if( avgscore_2 < SCORE_THRES_REJECT )then
                    where( clust_info_arr(:)%scoreclust == scoreclust_2 ) clust_info_arr(:)%good_bad = 0
                else
                    where( clust_info_arr(:)%scoreclust == scoreclust_2 ) clust_info_arr(:)%good_bad = 2
                endif
            endif
            if( DEBUG )then
                print *, 'found '//int2str(count(clust_info_arr(:)%good_bad == 1))//' 1st rate cluster(s) of class averages'
                print *, 'found '//int2str(count(clust_info_arr(:)%good_bad == 2))//' 2nd rate cluster(s) of class averages'
            endif
        end subroutine identify_good_bad_clusters

        subroutine calc_jointscore
            real,    allocatable :: score_vecs(:,:), dmat_scores(:,:), jointscores(:)
            integer, allocatable :: i_medoids_score(:), labels_score(:)
            integer :: iclust, jclust, nclust_score
            if( nclust <= 3 ) return
            ! extract feature vectors    
            allocate(score_vecs(nclust,4))
            do iclust = 1, nclust
                score_vecs(iclust,1) = clust_info_arr(iclust)%homogeneity
                score_vecs(iclust,2) = clust_info_arr(iclust)%resscore
                score_vecs(iclust,3) = clust_info_arr(iclust)%clustscore
                score_vecs(iclust,4) = clust_info_arr(iclust)%jointscore
            end do
            ! create distance matrix
            allocate(dmat_scores(nclust,nclust), source=0.)
            do iclust = 1, nclust - 1
                do jclust = iclust + 1, nclust
                    dmat_scores(iclust,jclust) = euclid(score_vecs(iclust,:3),score_vecs(jclust,:3))
                    dmat_scores(jclust,iclust) = dmat_scores(iclust,jclust)
                end do
            end do
            call normalize_minmax(dmat_scores)
            ! cluster
            call cluster_dmat(dmat_scores, 'aprop', nclust_score, i_medoids_score, labels_score)
            if( nclust_score < 3 )then
                nclust_score = 3
                call cluster_dmat(dmat_scores, 'kmed', nclust_score, i_medoids_score, labels_score)
            endif
            ! calculate discretized joint scores
            allocate(jointscores(nclust_score), source=0.)
            do iclust = 1, nclust_score
                jointscores(iclust) = sum(score_vecs(:,4), mask=labels_score == iclust) / real(count(labels_score == iclust))
            end do
            ! set discretized score
            do iclust = 1, nclust
                clust_info_arr(iclust)%jointscore = jointscores(labels_score(iclust))
                clust_info_arr(iclust)%scoreclust = labels_score(iclust)
            end do
        end subroutine calc_jointscore
 
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
        integer :: nmatch, nrefs, ldim(3)
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

        print *, 'nrefs   ', nrefs
        print *, 'nmatch  ', nmatch
        print *, 'smpd    ', smpd
        print *, 'ldim(1) ', ldim(1)

        ! ensure correct smpd/box in params class
        params%smpd = smpd
        params%box  = ldim(1)
        params%msk  = min(real(params%box/2)-COSMSKHALFWIDTH-1., 0.5*params%mskdiam /params%smpd)
        ! do the matching
        algninfo = match_imgs2refs(nrefs, nmatch, params%hp, params%lp, params%trs, cavg_imgs_ref, cavg_imgs_match)
        
        print *, algninfo

    end subroutine exec_match_cavgs

end module simple_commander_cavgs
