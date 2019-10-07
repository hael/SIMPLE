! array allocation for concrete strategy2D extensions to improve caching and reduce alloc overheads
module simple_strategy2D_alloc
include 'simple_lib.f08'
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
implicit none

public :: s2D, clean_strategy2D, prep_strategy2D
private
#include "simple_local_flags.inc"

type strategy2D_alloc
    ! per class
    integer, allocatable :: cls_pops(:)
    integer, allocatable :: cls_chunk(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
    ! time series
    integer, allocatable :: ptcls2neigh(:,:)
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    subroutine clean_strategy2D
        ! per class
        if( allocated(s2D%cls_pops)  ) deallocate(s2D%cls_pops)
        if( allocated(s2D%cls_chunk) ) deallocate(s2D%cls_chunk)
        ! per particle
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)) deallocate(s2D%do_inplsrch)
        ! time series
        if( allocated(s2D%ptcls2neigh)) deallocate(s2D%ptcls2neigh)
    end subroutine clean_strategy2D

    subroutine prep_strategy2D( ptcl_mask, which_iter )
        logical, intent(in)  :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer, intent(in)  :: which_iter
        type(ran_tabu)       :: rt
        integer              :: iptcl,prev_class,cnt,nptcls
        nptcls = count(ptcl_mask)
        ! gather class populations
        if( build_glob%spproj_field%isthere('class') )then
            call build_glob%spproj_field%get_pops(s2D%cls_pops, 'class', consider_w=params_glob%l_ptclw, maxn=params_glob%ncls)
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("simple_strategy2D_srch :: prep4strategy2D_srch")
        endif
        ! chunks, 1 by default
        allocate(s2D%cls_chunk(params_glob%ncls), source=1)
        ! in plane search
        allocate(s2D%do_inplsrch(1:nptcls), source=(params_glob%l_doshift .and. which_iter>2))
        ! stochastic search order
        allocate(s2D%srch_order(1:nptcls, params_glob%ncls), source=0)
        rt  = ran_tabu(params_glob%ncls)
        cnt = 0
        do iptcl = params_glob%fromp, params_glob%top
            if(ptcl_mask(iptcl))then
                cnt = cnt + 1
                call rt%ne_ran_iarr(s2D%srch_order(cnt,:))
                prev_class = nint(build_glob%spproj_field%get(iptcl,'class'))
                call put_last(prev_class, s2D%srch_order(cnt,:))
            endif
        end do
        call rt%kill
        if( any(s2D%srch_order == 0) ) THROW_HARD('Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch')
        if( trim(params_glob%tseries) .eq. 'yes')&
        &call build_glob%spproj_field%get_tseries_neighs(nint(params_glob%winsz),s2D%ptcls2neigh)
    end subroutine prep_strategy2D

end module simple_strategy2D_alloc
