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
    logical, allocatable :: cls_withdraw(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    subroutine clean_strategy2D
        ! per class
        if( allocated(s2D%cls_pops)    ) deallocate(s2D%cls_pops)
        if( allocated(s2D%cls_chunk)   ) deallocate(s2D%cls_chunk)
        if( allocated(s2D%cls_withdraw)) deallocate(s2D%cls_withdraw)
        ! per particle
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)) deallocate(s2D%do_inplsrch)
    end subroutine clean_strategy2D

    subroutine prep_strategy2D( ptcl_mask, which_iter, l_stream )
        logical, intent(in)  :: ptcl_mask(params_glob%fromp:params_glob%top), l_stream
        integer, intent(in)  :: which_iter
        type(ran_tabu)       :: rt
        integer, allocatable :: finds(:)
        integer              :: icls, iptcl, prev_class, cnt, nptcls
        real                 :: res0143, res05, ave, sdev, var
        logical              :: err
        nptcls = count(ptcl_mask)
        ! gather class populations
        allocate(s2D%cls_withdraw(params_glob%ncls), source=.true., stat=alloc_stat)
        if( build_glob%spproj_field%isthere('class') )then
            call build_glob%spproj_field%get_pops(s2D%cls_pops, 'class', consider_w=.false., maxn=params_glob%ncls)
            allocate(finds(params_glob%ncls))
            do icls=1,params_glob%ncls
                call build_glob%projfrcs%estimate_res(icls, res05, res0143)
                finds(icls) = calc_fourier_index(res0143, params_glob%boxmatch, params_glob%smpd)
            enddo
            call moment( real(finds), ave, sdev, var, err )
            do icls=1,params_glob%ncls
                s2D%cls_withdraw(icls) = finds(icls) > max(6,floor(ave))
            enddo
            deallocate(finds)
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("simple_strategy2D_srch :: prep4strategy2D_srch")
        endif
        ! chunks, 1 by default
        allocate(s2D%cls_chunk(params_glob%ncls), source=1)
        if( l_stream )then
            do icls=1,params_glob%ncls
                s2D%cls_chunk(icls) = ceiling(real(icls)/real(params_glob%ncls_start))
            enddo
        endif
        ! in plane search
        allocate(s2D%do_inplsrch(1:nptcls), source=params_glob%l_doshift)
        if( l_stream )then
            cnt = 0
            do iptcl = params_glob%fromp, params_glob%top
                if(ptcl_mask(iptcl))then
                    cnt = cnt + 1
                    s2D%do_inplsrch(cnt) = (build_glob%spproj_field%get(iptcl,'updatecnt') > 10.)
                endif
            end do
        endif
        ! make random reference direction order
        allocate(s2D%srch_order(1:nptcls, params_glob%ncls), source=0)
        rt  = ran_tabu(params_glob%ncls)
        cnt = 0
        do iptcl = params_glob%fromp, params_glob%top
            if(ptcl_mask(iptcl))then
                cnt = cnt + 1
                call rt%ne_ran_iarr( s2D%srch_order(cnt,:) )
                prev_class = nint(build_glob%spproj_field%get(iptcl,'class'))
                call put_last(prev_class, s2D%srch_order(cnt,:))
            endif
        end do
        if( any(s2D%srch_order == 0) ) THROW_HARD('Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch')
        call rt%kill
    end subroutine prep_strategy2D

end module simple_strategy2D_alloc
