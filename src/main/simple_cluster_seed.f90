! seeding methods for state labelling
module simple_cluster_seed
include 'simple_lib.f08'
use simple_oris, only: oris
implicit none

public :: gen_labelling
private

integer              :: nptcls      !< total # of particles
integer              :: nincl_ptcls !< # of particles with nonzero state
integer, allocatable :: states(:)   !< states for all particles

contains

    subroutine gen_labelling( os, nlabels, method )
        type(oris),       intent(inout) :: os
        integer,          intent(in)    :: nlabels
        character(len=*), intent(in)    :: method
        character(len=STDLEN) :: method_here
        ! init
        nptcls = os%get_noris()
        if(allocated(states))deallocate(states)
        states      = nint(os%get_all('state'))
        nincl_ptcls = count(states>0)
        ! switch
        if(.not.os%isthere('corr'))then
            ! defaults to random uniform
            method_here = 'uniform'
        else
            method_here = trim(method)
        endif
        select case(trim(method))
            case('uniform')
                call draw_uniform(os, nlabels)
            case('uniform_corr')
                call draw_uniform_corr(os, nlabels)
            case('dsquared')
                stop 'Unsupported method; simple_cluster_seed::gen_labelling'
            case DEFAULT
                stop 'Unsupported method; simple_cluster_seed::gen_labelling'
        end select
        ! updates labelling
        call os%set_all('state',real(states))
        ! cleanup
        deallocate(states)
    end subroutine gen_labelling

    !> random uniform labelling
    subroutine draw_uniform(os, nlabels)
        type(oris),           intent(inout) :: os
        integer,              intent(in)    :: nlabels
        integer, allocatable :: config(:), tmp(:)
        type(ran_tabu)       :: rt
        integer              :: iptcl, alloc_stat, cnt
        write(*,'(A)') '>>> GENERATING DIVERSE LABELING'
        allocate(config(nptcls), tmp(nincl_ptcls), source=0, stat=alloc_stat )
        if(alloc_stat /= 0) call allocchk('In: commander_hlev_wflows::diverse_labeling ', alloc_stat)
        rt = ran_tabu(nincl_ptcls)
        call rt%balanced(nlabels, tmp)
        cnt = 0
        do iptcl=1,nptcls
            if(states(iptcl)>0 .and. states(iptcl)<=nlabels)then
                cnt = cnt + 1
                config(iptcl) = tmp(cnt)
            else
                config(iptcl) = states(iptcl)
            endif
        enddo
        deallocate(tmp,states,config)
        call rt%kill
    end subroutine draw_uniform

    !>  partitions have a uniform distribution of correlations
    subroutine draw_uniform_corr(os, nlabels)
        type(oris),           intent(inout) :: os
        integer,              intent(in)    :: nlabels
        integer, allocatable :: tmp(:), order(:), config(:)
        real,    allocatable :: corrs(:)
        type(ran_tabu)       :: rt
        integer              :: iptcl, nptcls,nonzero_nptcls, alloc_stat, cnt, s, ind
        allocate(config(nptcls), order(nptcls), tmp(nlabels), source=0, stat=alloc_stat )
        corrs = os%get_all('corr')
        tmp   = (/(s,s=1,nlabels)/)
        where( states<=0 ) corrs = -1.
        order = (/(iptcl,iptcl=1,nptcls)/)
        call hpsort( corrs, order )
        call reverse(order)
        call reverse(corrs)
        rt = ran_tabu(nlabels)
        do iptcl = 1, nonzero_nptcls+nlabels, nlabels
            call rt%reset
            call rt%shuffle(tmp)
            do s = 1, nlabels
                ind = iptcl + s - 1
                if(ind > nptcls)exit
                config(order(ind)) = tmp(s)
            enddo
        enddo
        where( states<=0 ) config = 0
        deallocate(corrs,order,tmp,config)
        call rt%kill
    end subroutine draw_uniform_corr

    subroutine draw_dsquared(os, nlabels)
        type(oris),           intent(inout) :: os
        integer,              intent(in)    :: nlabels
        integer, allocatable :: tmp(:), config(:)
        real,    allocatable :: corrs(:)
        type(ran_tabu)       :: rt
        integer              :: iptcl, alloc_stat, cnt, s
        ! first partition: random uniform
        allocate(config(nptcls), tmp(nincl_ptcls), source=0, stat=alloc_stat )
        rt = ran_tabu(nincl_ptcls)
        call rt%balanced(nlabels, tmp)
        cnt = 0
        do iptcl=1,nptcls
            if(states(iptcl)>0 .and. states(iptcl)<=nlabels)then
                cnt = cnt + 1
                if(tmp(cnt)==1) config(iptcl) = 1
            endif
        enddo
        deallocate(tmp,config)
        call rt%kill
        !! TODO
    end subroutine draw_dsquared

end module simple_cluster_seed
