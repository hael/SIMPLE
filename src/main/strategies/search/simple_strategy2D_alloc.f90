!@descr: array allocation for concrete strategy2D extensions to improve caching and reduce alloc overheads
module simple_strategy2D_alloc
use simple_pftc_srch_api
implicit none

public :: s2D, clean_strategy2D, prep_strategy2D_batch, prep_strategy2D_glob
public :: prep_strategy2D_thread, set_strategy2D_stoch_bound, is_fresh_2D_start
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
    ! per class
    integer, allocatable :: cls_pops(:)
    ! per particle
    integer, allocatable :: srch_order(:,:)
    logical, allocatable :: do_inplsrch(:)
    ! per thread
    real,    allocatable :: class_space_corrs(:,:)     !< class vs. particle correlations
    real,    allocatable :: class_space_e3s(:,:)       !< in-plane Euler angles (e3)
    integer, allocatable :: class_space_inplinds(:,:)  !< in-plane rotation indices
end type strategy2D_alloc

type(strategy2D_alloc) :: s2D

contains

    logical function is_fresh_2D_start( params, which_iter )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: which_iter
        is_fresh_2D_start = params%startit <= 1 .and. which_iter <= params%startit &
            &.and. trim(params%continue) /= 'yes' .and. .not. params%l_fillin
    end function is_fresh_2D_start

    subroutine clean_strategy2D
        ! per class
        if( allocated(s2D%cls_pops)            ) deallocate(s2D%cls_pops)
        ! per particle
        if( allocated(s2D%srch_order)          ) deallocate(s2D%srch_order)
        if( allocated(s2D%do_inplsrch)         ) deallocate(s2D%do_inplsrch)
        ! per thread
        if( allocated(s2D%class_space_corrs)   ) deallocate(s2D%class_space_corrs)
        if( allocated(s2D%class_space_e3s)     ) deallocate(s2D%class_space_e3s)
        if( allocated(s2D%class_space_inplinds)) deallocate(s2D%class_space_inplinds)
    end subroutine clean_strategy2D

    ! Stochastic bound used by snhc* strategies
    subroutine set_strategy2D_stoch_bound( ncls, neigh_frac )
        integer, intent(in)  :: ncls
        real,    intent(in)  :: neigh_frac
        s2D%snhc_nrefs_bound = min(ncls, nint(real(ncls)*(1.-neigh_frac)))
        s2D%snhc_nrefs_bound = max(2, s2D%snhc_nrefs_bound)
    end subroutine set_strategy2D_stoch_bound

    !>  prep class & global parameters
    subroutine prep_strategy2D_glob( params, spproj, nrots, neigh_frac )
        use simple_eul_prob_tab_utils, only: calc_athres
        class(parameters), intent(in)    :: params
        type(sp_project),  intent(inout) :: spproj
        integer,           intent(in)    :: nrots
        real,              intent(in)    :: neigh_frac
        real    :: overlap, avg_dist_inpl
        logical :: l_fresh_start
        ! gather class populations
        l_fresh_start = is_fresh_2D_start(params, params%which_iter)
        if( l_fresh_start )then
            ! A fresh run with freshly generated references must not inherit
            ! class populations from a previous/streamed 2D assignment.
            allocate(s2D%cls_pops(params%ncls), source=MINCLSPOPLIM+1)
            ! Previous shift are also zeroed prior to and consistently with search
            call spproj%os_ptcl2D%zero_shifts
        else if( spproj%os_ptcl2D%isthere('class') )then
            ! The particle assignment table is the authoritative state for the
            ! next search. cls2D%pop can lag behind distributed assembly and
            ! would make valid previous classes look empty for overlap reporting.
            call spproj%os_ptcl2D%get_pops(s2D%cls_pops, 'class', maxn=params%ncls)
            where(s2D%cls_pops < 2 )
                ! ignoring classes with one particle
                s2D%cls_pops = 0
            end where
        else
            ! first iteration, no class assignment: all classes are up for grab
            allocate(s2D%cls_pops(params%ncls), source=MINCLSPOPLIM+1)
        endif
        if( all(s2D%cls_pops == 0) ) THROW_HARD('All class pops cannot be zero!')
        ! snhc/snhc_smpl
        call set_strategy2D_stoch_bound(params%ncls, neigh_frac)
        ! snhc_smpl
        s2D%snhc_smpl_ncls  = neighfrac2nsmpl(neigh_frac, params%ncls)
        s2D%snhc_smpl_ninpl = neighfrac2nsmpl(neigh_frac, nrots)
        select case(trim(params%refine))
        case('greedy_smpl','inpl_smpl')
            overlap        = spproj%os_ptcl2D%get_avg('mi_class', state=1)
            avg_dist_inpl  = calc_athres(os=spproj%os_ptcl2D, field_str='dist_inpl', prob_athres=params%prob_athres, state=1)
            avg_dist_inpl  = avg_dist_inpl * (1.-overlap)
            s2D%smpl_ninpl = max(2,nint(avg_dist_inpl*real(nrots)/180.))
            s2D%smpl_ncls  = nint(real(params%ncls) * (1.-overlap)**2)
            s2D%smpl_ncls  = max(1,min(s2D%smpl_ncls, ceiling(params%prob_athres/180.*real(params%ncls))))
        end select
        ! Sampling function power parameter
        if( (trim(params%stream2d)=='yes') )then
            s2D%power = EXTR_POWER
        else
            if( (params%extr_iter>params%extr_lim) )then
                if( (trim(params%refine)=='snhc_smpl') .or.&
                    &(trim(params%refine)=='snhc_smpl_many') )then
                    s2D%power = POST_EXTR_POWER
                endif
            endif
        endif
        ! per-thread class-space arrays
        if( allocated(s2D%class_space_corrs)    ) deallocate(s2D%class_space_corrs)
        if( allocated(s2D%class_space_e3s)      ) deallocate(s2D%class_space_e3s)
        if( allocated(s2D%class_space_inplinds) ) deallocate(s2D%class_space_inplinds)
        allocate(s2D%class_space_corrs(   params%ncls, nthr_glob),&
                &s2D%class_space_e3s(     params%ncls, nthr_glob),&
                &s2D%class_space_inplinds(params%ncls, nthr_glob))
        s2D%class_space_corrs    = -huge(1.0)
        s2D%class_space_e3s      = 0.
        s2D%class_space_inplinds = 0
    end subroutine prep_strategy2D_glob

    !>  prep batch related parameters (particles level)
    subroutine prep_strategy2D_batch( params, spproj, which_iter, nptcls, pinds )
        class(parameters),  intent(in) :: params  
        type(sp_project),   intent(in) :: spproj
        integer,            intent(in) :: which_iter
        integer,            intent(in) :: nptcls        ! # of particles in batch
        integer,            intent(in) :: pinds(nptcls)
        type(ran_tabu) :: rt
        integer        :: i,iptcl,prev_class
        logical        :: l_alloc, l_fresh_start
        l_alloc = .true.
        l_fresh_start = is_fresh_2D_start(params, which_iter)
        if( allocated(s2D%do_inplsrch) )then
            l_alloc = size(s2D%do_inplsrch) /= nptcls
            if( l_alloc ) deallocate(s2D%do_inplsrch,s2D%srch_order)
        endif
        if( l_alloc ) allocate(s2D%do_inplsrch(1:nptcls), s2D%srch_order(params%ncls, 1:nptcls))
        ! in plane search
        s2D%do_inplsrch = params%l_doshift .and. which_iter>2
        ! stochastic search order
        s2D%srch_order = 0
        rt = ran_tabu(params%ncls)
        do i = 1,nptcls
            iptcl = pinds(i)
            call rt%ne_ran_iarr(s2D%srch_order(:,i))
            if( l_fresh_start )then
                prev_class = 0
            else
                prev_class = nint(spproj%os_ptcl2D%get(iptcl,'class'))
            endif
            if( prev_class > 0 .and. prev_class <= params%ncls ) call put_last(prev_class, s2D%srch_order(:,i))
        enddo
        call rt%kill
        if( any(s2D%srch_order == 0) ) THROW_HARD('Invalid index in srch_order; simple_strategy2D_srch :: prep4strategy2D_srch')
    end subroutine prep_strategy2D_batch

    !> init thread-specific class-space search arrays (call per-particle before searching)
    subroutine prep_strategy2D_thread( ithr )
        integer, intent(in) :: ithr
        if( .not.allocated(s2D%class_space_corrs) .or. .not.allocated(s2D%class_space_e3s) .or. .not.allocated(s2D%class_space_inplinds) )then
            THROW_HARD('prep_strategy2D_thread called before class-space arrays were allocated')
        endif
        if( ithr < 1 .or. ithr > size(s2D%class_space_corrs,2) )then
            THROW_HARD('prep_strategy2D_thread received invalid thread index')
        endif
        s2D%class_space_corrs(   :,ithr) = -huge(1.0)
        s2D%class_space_e3s(     :,ithr) = 0.
        s2D%class_space_inplinds(:,ithr) = 0
    end subroutine prep_strategy2D_thread

end module simple_strategy2D_alloc
