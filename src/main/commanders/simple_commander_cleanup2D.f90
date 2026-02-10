!@ descr cleanup for molecule with large heterogeneity 
module simple_commanders_cleanup2D
use simple_commanders_api
use simple_commanders_project_core, only: commander_extract_subproj
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cavgs,        only: commander_cluster_cavgs
use simple_segmentation,            only: otsu_img
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cleanup2D
    contains
        procedure :: execute      => exec_cleanup2D
end type commander_cleanup2D

contains 
    subroutine exec_cleanup2D( self, cline )
        real,              parameter     :: SMPD_SHRINK1  = 4.0,  SIGMA_CRIT = 2., SIGMA_CRIT_MSK = 2.5
        integer,           parameter     :: BOXFAC = 3, NQ_DIAMS = 10
        class(commander_cleanup2D), intent(inout)  :: self
        class(cmdline),              intent(inout) :: cline
        type(cmdline), allocatable      :: cline_clust_js(:), cline_clust_pop(:)
        type(cmdline)                   :: cline_cluster_cavgs
        type(commander_abinitio2D)      :: xabinitio2D
        type(commander_cluster_cavgs)   :: xcluster_cavgs
        type(commander_extract_subproj) :: xextract_subproj
        type(image),          allocatable :: stk(:)
        type(image_bin), allocatable      :: binptcl_stk(:)
        integer,     allocatable :: labels(:), inds(:), pinds(:), diam_labels(:)
        logical,     allocatable :: inds_msk_js(:), inds_msk_pop(:)
        real,        allocatable :: joint_scores(:), diams_arr(:), diam_means(:), abs_z_scores(:), tmp(:), diams_arr_ts(:)
        type(sp_project), target :: spproj, spproj_fin
        type(parameters)     :: params
        type(string)         :: projname
        type(stats_struct)   :: diam_stats
        type(string)         :: selflag, ptcl_stk
        integer :: icls, ncls, iptcl, nptcls, pop, i, box_in_pix, icls_ccavgs, ncls_ccavgs
        real    :: mad, mskdiam
        ! need to set inputs in user interface
        call cline%set('polar', 'yes')
        call cline%set('ncls', 10)
        call cline%set('autoscale', 'yes')
        call cline%set('mskdiam', 100)
        call params%new(cline)
        call spproj%read_segment('mic',   params%projfile)
        call spproj%read_segment('stk',   params%projfile)
        call spproj%read_segment('ptcl2D',params%projfile)
        call xabinitio2D%execute_safe(cline)
        call cline%set('prune', 'no')
        call cline%delete('ncls')
        call xcluster_cavgs%execute_safe(cline)
        call cline%set('select_flag', 'class')
        call params%new(cline)
        call spproj%read(params%projfile)
        labels = spproj%os_cls2D%get_all_asint('class')
        ncls = size(labels)
        allocate(joint_scores(ncls))
        do icls = 1, ncls
            joint_scores(icls) = spproj%os_cls2D%get(icls, 'jointscore')
        end do
        allocate(inds_msk_js(ncls), source = .false.)
        where(joint_scores >= 40.) inds_msk_js = .true.
        deallocate(joint_scores, labels)
        allocate(cline_clust_js(count(inds_msk_js)), source = cline)
        do icls = 1, ncls 
            if(.not. inds_msk_js(icls)) cycle
            call cline_clust_js(icls)%set('clustind', icls)
            call cline_clust_js(icls)%set('subprojname', 'subproj' // int2str(icls))
            call xextract_subproj%execute_safe(cline_clust_js(icls))
            call spproj%read(params%projfile)
            nptcls = spproj%os_ptcl2D%get_noris()
            allocate(stk(nptcls))
            if( .not. cline%defined('stk' ) )then
                call spproj%get_cavgs_stk(ptcl_stk, nptcls, params%smpd, imgkind='ptcl')
                params%stk = ptcl_stk
            end if 
            allocate(binptcl_stk(nptcls), diams_arr(nptcls))
            do iptcl = 1, nptcls
                call stk(iptcl)%new([params%box,params%box,1], params%smpd)
                call stk(iptcl)%read(params%stk, iptcl)
                call binptcl_stk(i)%copy(stk(iptcl))
                call otsu_img(binptcl_stk(iptcl))
                call binptcl_stk(iptcl)%erode()
                call binptcl_stk(iptcl)%erode()
                call binptcl_stk(iptcl)%diameter_bin(diams_arr(iptcl))
            end do 
            deallocate(stk, binptcl_stk)
            ! calc mskdiams in same way as stream
            diams_arr = diams_arr + 2. * SMPD_SHRINK1
            call calc_stats(diams_arr, diam_stats)
            allocate( diam_means(NQ_DIAMS), diam_labels(NQ_DIAMS) )
            call sortmeans(diams_arr, NQ_DIAMS, diam_means, diam_labels)
            mad = mad_gau(diams_arr, diam_stats%med)
            allocate(abs_z_scores(NQ_DIAMS), source=abs((diam_means - diam_stats%med) / mad))
            do i = 2, NQ_DIAMS - 1
                if( abs_z_scores(i) < SIGMA_CRIT )then
                    pop = pop + count(diam_labels == i)
                    tmp = pack(diams_arr, mask=diam_labels == i)
                    if( allocated(diams_arr_ts) )then
                        diams_arr_ts = [diams_arr_ts(:), tmp(:)]
                        deallocate(tmp)
                    else
                        allocate(diams_arr_ts(size(tmp)), source=tmp)
                        deallocate(tmp)
                    endif
                endif
            end do
            call calc_stats(diams_arr_ts, diam_stats)
            ! what to do with box... 
            box_in_pix = find_magic_box(BOXFAC * nint(diam_stats%med/params%smpd))
            mskdiam = min((real(box_in_pix) - COSMSKHALFWIDTH) * params%smpd, diam_stats%maxv)
            mskdiam = min(diam_stats%avg + SIGMA_CRIT_MSK * diam_stats%sdev, mskdiam)
            call cline_clust_js(icls)%set('mskdiam', mskdiam)
            deallocate(stk, diams_arr)
        end do 
        ! do in parallel
        do icls = 1, ncls 
            if(.not. inds_msk_js(icls)) cycle
            call xabinitio2D%execute_safe(cline_clust_js(icls))
            call xcluster_cavgs%execute_safe(cline_clust_js(icls))
            call params%new(cline_clust_js(icls))
            call spproj%read(params%projfile)
            ncls_ccavgs = spproj%os_ptcl2D%get_noris()
            allocate(labels(ncls_ccavgs), inds_msk_pop(ncls_ccavgs))
            labels = spproj%os_cls2D%get_all_asint('class')
            inds_msk_pop = .true. 
            ! 300 ptcls/ clusters and 300 clusters
            do while(icls_ccavgs < 300) 
                pop = spproj%os_cls2D%get(icls_ccavgs, 'pop')
                if(pop > 300) then
                    inds_msk_pop(icls_ccavgs) = .false.  
                    cycle 
                end if
                if(icls_ccavgs == ncls_ccavgs) exit 
                icls_ccavgs = icls_ccavgs + 1
            end do
            allocate(cline_clust_pop(count(inds_msk_pop)))
            do icls_ccavgs = 1, ncls_ccavgs
                if(.not. inds_msk_pop(icls_ccavgs)) cycle
                
            end do
            ! extract subprojects based on inds_msk_pop 
            ! merge the ones that satify the the pop condition
            ! merge again based on jointscore condition 
        
            call params%new(cline_clust_js(icls))
            call spproj%read(params%projfile)
        end do 
        ! prune from rest of ptcls
        deallocate(inds_msk_js)
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****', print_simple = .false.)
    end subroutine exec_cleanup2D

    end module simple_commanders_cleanup2D
    
