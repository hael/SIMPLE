! array allocation for concrete strategy2D extensions to improve caching and reduce alloc overheads
module simple_strategy2D_alloc
include 'simple_lib.f08'
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
implicit none

public :: s2D, clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob
private
#include "simple_local_flags.inc"

type strategy2D_alloc
    ! per class
    integer, allocatable :: cls_pops(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
    ! per thread
    real,    allocatable :: cls_corrs(:,:)
    logical, allocatable :: cls_searched(:,:)
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    subroutine clean_strategy2D
        ! per class
        if( allocated(s2D%cls_pops)   ) deallocate(s2D%cls_pops)
        ! per particle
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)) deallocate(s2D%do_inplsrch)
        ! per thread
        if( allocated(s2D%cls_corrs))    deallocate(s2D%cls_corrs)
        if( allocated(s2D%cls_searched)) deallocate(s2D%cls_searched)
    end subroutine clean_strategy2D

    !>  prep class & global parameters
    subroutine prep_strategy2D_glob
        logical :: zero_oris, ncls_diff
        zero_oris = build_glob%spproj%os_cls2D%get_noris() == 0
        ncls_diff = build_glob%spproj%os_cls2D%get_noris() /= params_glob%ncls
        ! gather class populations
        if( build_glob%spproj_field%isthere('class') )then
            if( zero_oris .or. ncls_diff  )then
                ! the ==0    is to overcome bug in shared-memory version
                ! the ==ncls is to be able to restart after having run cleanup with fewer classes
                if( zero_oris )then
                    call build_glob%spproj%os_ptcl2D%get_pops(s2D%cls_pops, 'class', consider_w=.false., maxn=params_glob%ncls)
                    return
                endif
                if( ncls_diff )then
                    allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1)
                    return
                endif
            endif
            s2D%cls_pops = nint(build_glob%spproj%os_cls2D%get_all('pop'))
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1)
        endif
        if( all(s2D%cls_pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        ! per-thread allocations
        allocate(s2D%cls_corrs(params_glob%ncls,nthr_glob),source=0.0)
        allocate(s2D%cls_searched(params_glob%ncls,nthr_glob),source=.false.)
    end subroutine prep_strategy2D_glob

    !>  prep batch related parameters (particles level)
    subroutine prep_strategy2D_batch( pftcc, which_iter, nptcls, pinds )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        type(polarft_corrcalc), intent(in) :: pftcc
        integer,                intent(in) :: which_iter
        integer,                intent(in) :: nptcls        ! # of particles in batch
        integer,                intent(in) :: pinds(nptcls)
        type(ran_tabu) :: rt
        integer        :: i,iptcl,prev_class
        logical        :: l_alloc
        l_alloc = .true.
        if( allocated(s2D%do_inplsrch) )then
            l_alloc = size(s2D%do_inplsrch) /= nptcls
            if( l_alloc ) deallocate(s2D%do_inplsrch,s2D%srch_order)
        endif
        if( l_alloc ) allocate(s2D%do_inplsrch(1:nptcls), s2D%srch_order(1:nptcls, params_glob%ncls))
        ! in plane search
        s2D%do_inplsrch = params_glob%l_doshift .and. which_iter>2
        ! stochastic search order
        s2D%srch_order = 0
        rt = ran_tabu(params_glob%ncls)
        do i = 1,nptcls
            iptcl = pinds(i)
            call rt%ne_ran_iarr(s2D%srch_order(i,:))
            prev_class = nint(build_glob%spproj_field%get(iptcl,'class'))
            call put_last(prev_class, s2D%srch_order(i,:))
        enddo
        call rt%kill
        if( any(s2D%srch_order == 0) ) THROW_HARD('Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch')
    end subroutine prep_strategy2D_batch

end module simple_strategy2D_alloc
