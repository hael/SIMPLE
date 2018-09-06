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
    logical, allocatable :: cls_topseg(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    subroutine clean_strategy2D
        ! per class
        if( allocated(s2D%cls_pops)  ) deallocate(s2D%cls_pops)
        if( allocated(s2D%cls_chunk) ) deallocate(s2D%cls_chunk)
        if( allocated(s2D%cls_topseg)) deallocate(s2D%cls_topseg)
        ! per particle
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)) deallocate(s2D%do_inplsrch)
    end subroutine clean_strategy2D

    subroutine prep_strategy2D( ptcl_mask, which_iter, l_stream )
        use simple_image, only: image
        logical, intent(in)  :: ptcl_mask(params_glob%fromp:params_glob%top), l_stream
        integer, intent(in)  :: which_iter
        type(ran_tabu)       :: rt, rt_best, rt_worst
        type(image)          :: img
        real,    allocatable :: corrs(:)
        integer, allocatable :: inds(:)
        real                 :: ave, sdev, threshold
        integer              :: n, icls, iptcl, prev_class, cnt, nptcls, n_worst, updatecnt
        nptcls = count(ptcl_mask)
        ! gather class populations
        if( build_glob%spproj_field%isthere('class') )then
            call build_glob%spproj_field%get_pops(s2D%cls_pops, 'class', consider_w=.false., maxn=params_glob%ncls)
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1, stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("simple_strategy2D_srch :: prep4strategy2D_srch")
        endif
        ! chunks, 1 by default
        allocate(s2D%cls_chunk(params_glob%ncls), source=1)
        if( l_stream )then
            do icls=1,params_glob%ncls
                s2D%cls_chunk(icls) = nint(build_glob%spproj%os_cls2D%get(icls,'chunk'))
            enddo
        endif
        ! in plane search
        allocate(s2D%do_inplsrch(1:nptcls), source=params_glob%l_doshift)
        if( l_stream )then
            cnt = 0
            do iptcl = params_glob%fromp, params_glob%top
                if(ptcl_mask(iptcl))then
                    cnt = cnt + 1
                    s2D%do_inplsrch(cnt) = (build_glob%spproj_field%get(iptcl,'updatecnt') > 5.)
                endif
            end do
        endif
        ! make random reference direction order
        allocate(s2D%srch_order(1:nptcls, params_glob%ncls), source=0)
        allocate(s2D%cls_topseg(params_glob%ncls), source=.true.)
        rt  = ran_tabu(params_glob%ncls)
        if( l_stream .and. build_glob%spproj%os_cls2D%isthere('corr') )then
            ! determine threshold
            corrs     = build_glob%spproj%os_cls2D%get_all('corr')
            n         = count(corrs>0.0001)
            ave       = sum(corrs,mask=(corrs>0.0001))/real(n)
            sdev      = sqrt(sum((corrs-ave)**2.,mask=(corrs>0.0001))/real(n))
            threshold = ave-sdev
            n_worst = count(corrs<threshold .and. corrs>0.0001)
            ! set top ranking mask
            cnt = 0
            call img%new([params_glob%box,params_glob%box,1],params_glob%smpd)
            do icls=1,params_glob%ncls
                if( corrs(icls) > 0.0001 )s2D%cls_topseg(icls) = (corrs(icls) > threshold)
                if( params_glob%part==1 )then
                    if(.not.s2D%cls_topseg(icls))then
                        cnt = cnt +1
                        call img%read(params_glob%refs,icls)
                        call img%write('rejects_'//int2str(which_iter)//'.mrc',cnt)
                    endif
                endif
            enddo
            call img%kill
            ! generates search orders
            inds = (/(icls,icls=1,params_glob%ncls)/)
            where( s2D%cls_topseg )
                corrs = 1.
            else where
                corrs = 0.
            end where
            call hpsort(corrs,inds)
            rt_best  = ran_tabu(params_glob%ncls-n_worst)
            rt_worst = ran_tabu(n_worst)
            cnt = 0
            do iptcl=params_glob%fromp,params_glob%top
                if(.not.ptcl_mask(iptcl)) cycle
                cnt = cnt + 1
                if(build_glob%spproj_field%get(iptcl,'state')<0.5)then
                    call rt%ne_ran_iarr(s2D%srch_order(cnt,:))
                    cycle
                endif
                updatecnt  = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                prev_class = nint(build_glob%spproj_field%get(iptcl,'class'))
                if(updatecnt>=2 .and. updatecnt<=10)then
                    s2D%srch_order(cnt,:) = inds
                    call rt_worst%shuffle(s2D%srch_order(cnt,1:n_worst))
                    call rt_best%shuffle(s2D%srch_order(cnt,n_worst+1:params_glob%ncls))
                    if(s2D%cls_topseg(prev_class)) call reverse(s2D%srch_order(cnt,:))
                else
                    call rt%ne_ran_iarr(s2D%srch_order(cnt,:))
                endif
                call put_last(prev_class, s2D%srch_order(cnt,:))
            enddo
            deallocate(corrs,inds)
            call rt_best%kill
            call rt_worst%kill
        else
            cnt = 0
            do iptcl = params_glob%fromp, params_glob%top
                if(ptcl_mask(iptcl))then
                    cnt = cnt + 1
                    call rt%ne_ran_iarr(s2D%srch_order(cnt,:))
                    prev_class = nint(build_glob%spproj_field%get(iptcl,'class'))
                    call put_last(prev_class, s2D%srch_order(cnt,:))
                endif
            end do
        endif
        call rt%kill
        if( any(s2D%srch_order == 0) ) THROW_HARD('Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch')
    end subroutine prep_strategy2D

end module simple_strategy2D_alloc
