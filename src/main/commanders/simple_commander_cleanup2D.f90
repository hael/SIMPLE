!@descr: cleanup for molecule with large heterogeneity 
module simple_commanders_cleanup2D
use simple_commanders_api
use simple_commanders_project_core, only: commander_extract_subproj
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cavgs,        only: commander_cluster_cavgs
use simple_segmentation,            only: otsu_img
use simple_projfile_utils,          only: merge_chunk_projfiles
use simple_imgarr_utils,            only: read_cavgs_into_imgarr
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cleanup2D
    contains
        procedure :: execute      => exec_cleanup2D
end type commander_cleanup2D

contains 
    subroutine exec_cleanup2D( self, cline )
        class(commander_cleanup2D), intent(inout)  :: self
        class(cmdline),              intent(inout) :: cline
        type(cmdline), allocatable      :: cline_clust_js(:), cline_clust_pop(:)
        type(commander_abinitio2D)      :: xabinitio2D
        type(commander_cluster_cavgs)   :: xcluster_cavgs
        type(commander_extract_subproj) :: xextract_subproj
        type(image),  allocatable       :: stk(:)
        type(string), allocatable       :: chunk_fnames_js(:), chunk_fnames_pop(:)
        integer, allocatable            :: labels(:), pinds(:)
        logical, allocatable            :: inds_msk_js(:), inds_msk_pop(:)
        real,    allocatable            :: joint_scores(:), diams_arr(:), shifts(:,:), fin_mskdiams(:)
        type(sp_project)     :: spproj, spproj_merged_1, spproj_merged_fin
        type(parameters)     :: params
        type(stats_struct)   :: diam_stats
        type(string)         :: folder
        integer :: icls, ncls, pop, icls_ccavgs, ncls_ccavgs, jcls, ncls_js 
        integer :: icls2, ncls2, jcls2, n_cavgs_subproj, n_ptcls_subproj
        real    :: mskdiam
        call cline%set('polar', 'yes')
        call cline%set('ncls', 10)
        call cline%set('autoscale', 'yes')
        call params%new(cline)
        call spproj%read_segment('mic',   params%projfile)
        call spproj%read_segment('stk',   params%projfile)
        call spproj%read_segment('ptcl2D',params%projfile)
        call xabinitio2D%execute_safe(cline)
        call cline%set('prune', 'no')
        ! if ncls set, defaults to kmedoids clustering
        call cline%delete('ncls')
        call xcluster_cavgs%execute_safe(cline)
        call params%new(cline)
        call spproj%read(params%projfile)
        labels = spproj%os_cls2D%get_all_asint('cluster')
        ncls = maxval(labels)
        allocate(joint_scores(ncls))
        do icls = 1, ncls
            joint_scores(icls) = spproj%os_cls2D%get(icls, 'jointscore')
            joint_scores(icls) = 50. 
        end do
        allocate(inds_msk_js(ncls), source = .false.)
        where(joint_scores >= 40.) inds_msk_js = .true.
        deallocate(joint_scores)
        allocate(cline_clust_js(count(inds_msk_js)), source=cline)
        allocate(chunk_fnames_js(count(inds_msk_js)))
        allocate(fin_mskdiams(count(inds_msk_js)))
        print *,'# of clusters of cavgs', ncls
        jcls = 0
        do icls = 1, ncls
            if(.not. inds_msk_js(icls)) cycle
            print *, 'cluster of cavgs idx:', icls
            jcls = jcls + 1
            call cline_clust_js(jcls)%set('mkdir', 'yes')
            call cline_clust_js(jcls)%set('prune', 'yes')
            call cline_clust_js(jcls)%set('clustind', icls)
            call cline_clust_js(jcls)%set('subprojname', 'subproj' // int2str(icls))
            call xextract_subproj%execute_safe(cline_clust_js(jcls))
            call params%new(cline_clust_js(jcls))
            ! original project file, updated params
            call spproj%os_cls2D%get_pinds(params%clustind, 'cluster', pinds)
            n_cavgs_subproj = size(pinds)
            deallocate(pinds)
            call spproj%os_ptcl2D%get_pinds(params%clustind, 'cluster', pinds)
            n_ptcls_subproj = size(pinds)
            deallocate(pinds)
            ! ! read cavgs and calculate mask
            print *, 'n_cavgs / subproj', n_cavgs_subproj, 'n_ptcls / subproj', n_ptcls_subproj

            allocate(stk(n_cavgs_subproj))
            stk = read_cavgs_into_imgarr(spproj, labels == jcls)
            call automask2D(stk, params%ngrow, nint(params%winsz), params%edge, diams_arr, shifts)
            deallocate(stk)
            call calc_stats(diams_arr, diam_stats)
            mskdiam = diam_stats%med
            fin_mskdiams(jcls) = mskdiam
            print *, 'mskdiam before', params%mskdiam 
            call cline_clust_js(jcls)%set('mskdiam', mskdiam)
            print *, 'mskdiam after', mskdiam

            ! compute abinitio2D and cluster cavgs with more suitable mask
            call xabinitio2D%execute_safe(cline_clust_js(jcls))
            call xcluster_cavgs%execute_safe(cline_clust_js(jcls))
            call params%new(cline_clust_js(jcls))
            call spproj%read(params%projfile)

            ! ! logic to filter subsubcls
            ncls_ccavgs = spproj%os_cls2D%get_noris()
            allocate(labels(ncls_ccavgs), inds_msk_pop(ncls_ccavgs))
            labels = spproj%os_cls2D%get_all_asint('class')
            inds_msk_pop = .true.
            ! ! first 300 clusters, only keep clusters with max 300 ptcls
            do icls_ccavgs = 1, min(300, ncls_ccavgs)
                ! pop = spproj%os_cls2D%get(icls_ccavgs, 'pop')
                if (pop > 300) inds_msk_pop(icls_ccavgs) = .false.
            end do
            allocate(cline_clust_pop(count(inds_msk_pop)), source = cline_clust_js(jcls))
            allocate(chunk_fnames_pop(count(inds_msk_pop)))
            jcls2 = 0
            ncls2 = size(spproj%os_cls2D%get_all_asint('class'))
            do icls2 = 1, ncls2
                if(.not. inds_msk_pop(icls2)) cycle
                jcls2 = jcls2 + 1
                call cline_clust_pop(jcls2)%set('clustind', icls2)
                call cline_clust_pop(jcls2)%set('subprojname', 'subsubproj' // int2str(icls2))
                call xextract_subproj%execute_safe(cline_clust_pop(jcls2))
                call cline_clust_pop(jcls2)%set('prune', 'yes')
                call params%new(cline_clust_pop(jcls2))
                chunk_fnames_pop(jcls2) = params%projfile
            end do
            deallocate(cline_clust_pop)
            folder         = "."

            ! first merge
            call merge_chunk_projfiles(chunk_fnames_pop, folder, spproj_merged_1)
            call spproj_merged_1%write(params%projfile_merged)
            deallocate(chunk_fnames_pop)
        end do
        ! set mskdiam to maxval of subprojects and prune subprojects
        jcls = 0
        do icls = 1, ncls
            if(.not. inds_msk_js(icls)) cycle
            jcls = jcls + 1
            call cline_clust_js(jcls)%set('prune', 'yes')
            call cline_clust_js(jcls)%set('mskdiam', maxval(fin_mskdiams))
            call params%new(cline_clust_js(jcls))
            chunk_fnames_js(jcls) = params%projfile
        end do
        folder             = "."
        ! merge
        call merge_chunk_projfiles(chunk_fnames_js, folder, spproj_merged_fin)
        call spproj_merged_fin%write(params%projfile_merged)
        print *, spproj_merged_fin%os_ptcl2D%get_noris()
        deallocate(inds_msk_js)
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****', print_simple = .false.)
    end subroutine exec_cleanup2D    
end module simple_commanders_cleanup2D
    
