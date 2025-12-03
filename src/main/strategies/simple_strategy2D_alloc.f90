! array allocation for concrete strategy2D extensions to improve caching and reduce alloc overheads
module simple_strategy2D_alloc
include 'simple_lib.f08'
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
use simple_polarft_calc, only: pftc_glob
use simple_eul_prob_tab2D,   only: eul_prob_tab2D, neighfrac2nsmpl
implicit none

public :: s2D, clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob
private
#include "simple_local_flags.inc"

type strategy2D_alloc
    ! global parameters
    real                 :: power            = EXTR_POWER ! power for the sampling function
    integer              :: snhc_nrefs_bound = 0        ! refine=snhc
    integer              :: snhc_smpl_ncls   = 0        !       =snhc_smpl
    integer              :: snhc_smpl_ninpl  = 0        !       =snhc_smpl
    integer              :: smpl_ninpl       = 0        !       =greedy_smpl
    integer              :: smpl_ncls        = 0        !       =greedy_smpl
    class(eul_prob_tab2D), pointer :: probtab => null() ! pointer to probabibility table (refine=prob*)
    ! per class
    integer, allocatable :: cls_pops(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    subroutine clean_strategy2D
        if( associated(s2D%probtab)   ) nullify(s2D%probtab)
        ! per class
        if( allocated(s2D%cls_pops)   ) deallocate(s2D%cls_pops)
        ! per particle
        if( allocated(s2D%srch_order) ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)) deallocate(s2D%do_inplsrch)
    end subroutine clean_strategy2D

    !>  prep class & global parameters
    subroutine prep_strategy2D_glob( neigh_frac )
        use simple_eul_prob_tab, only: calc_athres
        real,           intent(in) :: neigh_frac
        real    :: overlap, avg_dist_inpl
        logical :: zero_oris, ncls_diff
        ! gather class populations
        zero_oris = build_glob%spproj%os_cls2D%get_noris() == 0
        ncls_diff = build_glob%spproj%os_cls2D%get_noris() /= params_glob%ncls
        if( build_glob%spproj_field%isthere('class') )then
            if( zero_oris .or. ncls_diff  )then
                ! the ==0    is to overcome bug in shared-memory version
                ! the ==ncls is to be able to restart after having run cleanup with fewer classes
                if( zero_oris )then
                    call build_glob%spproj%os_ptcl2D%get_pops(s2D%cls_pops, 'class', maxn=params_glob%ncls)
                    return
                endif
                if( ncls_diff )then
                    allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1)
                    return
                endif
            endif
            s2D%cls_pops = nint(build_glob%spproj%os_cls2D%get_all('pop'))
            where(s2D%cls_pops < 2 )
                ! ignoring classes with one particle
                s2D%cls_pops = 0
            end where
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params_glob%ncls), source=MINCLSPOPLIM+1)
        endif
        if( all(s2D%cls_pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        ! snhc
        s2D%snhc_nrefs_bound = min(params_glob%ncls, nint(real(params_glob%ncls)*(1.-neigh_frac)))
        s2D%snhc_nrefs_bound = max(2, s2D%snhc_nrefs_bound)
        ! snhc_smpl
        s2D%snhc_smpl_ncls  = neighfrac2nsmpl(neigh_frac, params_glob%ncls)
        s2D%snhc_smpl_ninpl = neighfrac2nsmpl(neigh_frac, pftc_glob%get_nrots())
        select case(trim(params_glob%refine))
        case('greedy_smpl','inpl_smpl')
            overlap        = build_glob%spproj_field%get_avg('mi_class', state=1)
            avg_dist_inpl  = calc_athres('dist_inpl', state=1)
            avg_dist_inpl  = avg_dist_inpl * (1.-overlap)
            s2D%smpl_ninpl = max(2,nint(avg_dist_inpl*real(pftc_glob%get_nrots())/180.))
            s2D%smpl_ncls  = nint(real(params_glob%ncls) * (1.-overlap)**2)
            s2D%smpl_ncls  = max(1,min(s2D%smpl_ncls, ceiling(params_glob%prob_athres/180.*real(params_glob%ncls))))
        end select
        ! Sampling function power parameter
        if( (trim(params_glob%stream)=='yes') )then
            s2D%power = EXTR_POWER
        else
            if( (trim(params_glob%refine)=='snhc_smpl') .and. (params_glob%extr_iter>params_glob%extr_lim) )then
                s2D%power = POST_EXTR_POWER
            endif
        endif
    end subroutine prep_strategy2D_glob

    !>  prep batch related parameters (particles level)
    subroutine prep_strategy2D_batch( pftc, which_iter, nptcls, pinds )
        use simple_polarft_calc, only: polarft_calc
        type(polarft_calc), intent(in) :: pftc
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
