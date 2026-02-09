!@ descr cleanup for molecule with large heterogeneity 
module simple_commanders_cleanup2D
use simple_commanders_api
use simple_commanders_project_core, only: commander_extract_subproj
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cavgs,        only: commander_cluster_cavgs
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
        type(cmdline)  :: cline_cluster_cavgs
        type(commander_abinitio2D)      :: xabinitio2D
        type(commander_cluster_cavgs)   :: xcluster_cavgs
        type(commander_extract_subproj) :: xextract_subproj
        type(sp_project), target  :: spproj
        type(parameters)     :: params
        type(string)         :: projname
        type(stats_struct)   :: diam_stats
        type(string)         :: selflag
        type(image), allocatable :: stk(:)
        integer,     allocatable :: labels(:), inds(:), pinds(:)
        logical,     allocatable :: inds_msk(:)
        real,        allocatable :: joint_scores(:), diams_arr(:)
        integer :: icls, ncls, iptcls, nptcls
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
        allocate(inds_msk(ncls), source = .false.)
        where(joint_scores >= 40.) inds_msk = .true.
        deallocate(joint_scores, labels)
        do icls = 1, ncls 
            if(.not. inds_msk(icls)) cycle
            call cline%set('clustind', icls)
            call cline%set('subprojname', 'subproj' // int2str(icls))
            call xextract_subproj%execute_safe(cline)
            call spproj%read(params%projfile)
            nptcls = spproj%os_ptcl2D%get_noris()
            allocate(stk(nptcls))
            allocate(diams_arr(nptcls))
           
            call calc_stats(diams_arr, diam_stats)
            deallocate(stk, diams_arr)
        end do 

        deallocate(inds_msk)
        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****', print_simple = .false.)
    end subroutine exec_cleanup2D

    end module simple_commanders_cleanup2D
    
