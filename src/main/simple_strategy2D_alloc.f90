! array allocation for concrete strategy2D extensions to improve caching and reduce alloc overheads
module simple_strategy2D_alloc
include 'simple_lib.f08'
use simple_singletons
implicit none

public :: s2D, clean_strategy2D, prep_strategy2D
private
type strategy2D_alloc
integer, allocatable :: cls_pops(:)
integer, allocatable :: srch_order(:,:)
contains
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D
contains

    subroutine clean_strategy2D
        if( allocated(s2D%cls_pops)   ) deallocate(s2D%cls_pops)
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
    end subroutine clean_strategy2D

    subroutine prep_strategy2D( ptcl_mask, which_iter )
        logical,        intent(in)    :: ptcl_mask(p%fromp:p%top)
        integer,        intent(in)    :: which_iter
        type(ran_tabu) :: rt
        integer        :: iptcl, prev_class, cnt, nptcls
        ! gather class populations
        if( b%a%isthere('class') )then
            if( p%weights2D .eq. 'yes' .and. which_iter > 3 )then
                call b%a%get_pops(s2D%cls_pops, 'class', consider_w=.true., maxn=p%ncls)
            else
                call b%a%get_pops(s2D%cls_pops, 'class', consider_w=.false., maxn=p%ncls)
            endif
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(p%ncls), source=MINCLSPOPLIM+1, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("simple_strategy2D_matcher; prime2D_exec cls_pops")
        endif
        ! make random reference direction order
        nptcls = count(ptcl_mask)
        allocate(s2D%srch_order(1:nptcls, p%ncls), source=0)
        rt  = ran_tabu(p%ncls)
        cnt = 0
        do iptcl = p%fromp, p%top
            if(ptcl_mask(iptcl))then
                cnt = cnt + 1
                call rt%reset
                call rt%ne_ran_iarr( s2D%srch_order(cnt,:) )
                prev_class = nint(b%a%get(iptcl,'class'))
                call put_last(prev_class, s2D%srch_order(cnt,:))
            endif
        end do
        if( any(s2D%srch_order == 0) ) stop 'Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch'
        call rt%kill
    end subroutine prep_strategy2D

end module simple_strategy2D_alloc
