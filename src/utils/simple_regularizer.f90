! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_image
implicit none

public :: regularizer
private
#include "simple_local_flags.inc"

type reg_params
    integer :: iptcl        !< iptcl index
    integer :: iref         !< iref index
    integer :: loc          !< inpl index
    real    :: prob, sh(2)  !< probability, shift
end type reg_params

type :: regularizer
    integer              :: nrots
    integer              :: nrefs
    integer              :: pftsz
    integer              :: kfromto(2)
    integer              :: inpl_ns, refs_ns
    
    real,    allocatable :: ref_ptcl_cor(:,:)           !< 2D corr table
    integer, allocatable :: ptcl_ref_map(:)             !< hard-alignment tab
    real,    allocatable :: inpl_corr(:,:), refs_corr(:,:)
    integer, allocatable :: inpl_inds(:,:), refs_inds(:,:)
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(reg_params),        allocatable :: ref_ptcl_tab(:,:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: fill_tab_inpl_smpl
    procedure          :: tab_normalize
    procedure          :: shift_search
    procedure          :: nonuni_tab_align
    procedure, private :: ref_multinomal, inpl_multinomal
    ! DESTRUCTOR
    procedure          :: kill
end type regularizer

contains
    ! CONSTRUCTORS

    subroutine new( self, pftcc )
        class(regularizer),      target, intent(inout) :: self
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        integer :: iptcl, iref
        self%nrots   = pftcc%nrots
        self%nrefs   = pftcc%nrefs
        self%pftsz   = pftcc%pftsz
        self%kfromto = pftcc%kfromto
        self%inpl_ns = params_glob%reg_athres * self%nrots / 180.
        self%refs_ns = params_glob%reg_athres * self%nrefs / 180.
        self%pftcc => pftcc
        allocate(self%ref_ptcl_cor(self%nrefs,params_glob%fromp:params_glob%top),&
                &self%refs_corr(self%nrefs,params_glob%nthr), self%inpl_corr(self%nrots,params_glob%nthr), source=0.)
        allocate(self%inpl_inds(self%nrots,params_glob%nthr), self%refs_inds(self%nrefs,params_glob%nthr), source=0)
        allocate(self%ref_ptcl_tab(self%nrefs,params_glob%fromp:params_glob%top))
        allocate(self%ptcl_ref_map(params_glob%fromp:params_glob%top))
        do iptcl = params_glob%fromp,params_glob%top
            do iref = 1,self%nrefs
                self%ref_ptcl_tab(iref,iptcl)%iptcl = iptcl
                self%ref_ptcl_tab(iref,iptcl)%iref  = iref
                self%ref_ptcl_tab(iref,iptcl)%loc   = 0
                self%ref_ptcl_tab(iref,iptcl)%prob  = 0.
                self%ref_ptcl_tab(iref,iptcl)%sh    = 0.
            enddo
        enddo
    end subroutine new

    subroutine fill_tab_inpl_smpl( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        integer :: i, iref, iptcl, irnd
        real    :: inpl_corrs(self%nrots)
        call seed_rnd
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,inpl_corrs,irnd) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                iptcl = glob_pinds(i)
                ! sampling the inpl rotation
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                irnd = self%inpl_multinomal(inpl_corrs)
                self%ref_ptcl_tab(iref,iptcl)%sh  = 0.
                self%ref_ptcl_tab(iref,iptcl)%loc = irnd
                self%ref_ptcl_cor(iref,iptcl)     = inpl_corrs(irnd)
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_tab_inpl_smpl

    subroutine tab_normalize( self )
        use simple_builder,           only: build_glob
        class(regularizer), intent(inout) :: self
        type(ori) :: o
        integer   :: iref, iptcl, iref2
        real      :: sum_corr(self%nrefs)
        logical   :: lnns(self%nrefs), ref_neigh_tab(self%nrefs,self%nrefs)
        ref_neigh_tab = .false.
        do iref = 1, self%nrefs
            lnns = .false.
            call build_glob%eulspace%get_ori(iref, o)
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, params_glob%reg_athres, lnns)
            do iref2 = 1, self%nrefs
                if( iref2 /= iref .and. lnns(iref2) )then
                    ref_neigh_tab(iref,  iref2) = .true.
                    ref_neigh_tab(iref2, iref ) = .true.
                endif
            enddo
            ref_neigh_tab(iref, iref) = .true.
        enddo
        ! normalize so prob of each ptcl is between [0,1] for all refs
        if( params_glob%l_reg_norm )then
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_corr,iref)
            do iptcl = params_glob%fromp, params_glob%top
                do iref = 1, self%nrefs
                    sum_corr(iref) = sum(self%ref_ptcl_cor(:,iptcl), mask=ref_neigh_tab(:,iref))
                enddo
                do iref = 1, self%nrefs
                    if( sum_corr(iref) < TINY )then
                        self%ref_ptcl_cor(iref,iptcl) = 0.
                    else
                        self%ref_ptcl_cor(iref,iptcl) = self%ref_ptcl_cor(iref,iptcl) / sum_corr(iref)
                    endif
                enddo
            enddo
            !$omp end parallel do
        endif
        if( maxval(self%ref_ptcl_cor) < TINY )then
            self%ref_ptcl_cor = 0.
        else
            self%ref_ptcl_cor = self%ref_ptcl_cor / maxval(self%ref_ptcl_cor)
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(iref,iptcl)
        do iptcl = params_glob%fromp,params_glob%top
            do iref = 1, self%nrefs
                self%ref_ptcl_tab(iref,iptcl)%prob = self%ref_ptcl_cor(iref,iptcl)
            enddo
        enddo
        !$omp end parallel do
    end subroutine tab_normalize

    subroutine shift_search( self )
        use simple_pftcc_shsrch_reg, only: pftcc_shsrch_reg
        class(regularizer), intent(inout) :: self
        type(pftcc_shsrch_reg) :: grad_shsrch_obj(params_glob%nthr)
        integer :: iref, iptcl, ithr, irot
        real    :: lims(2,2), cxy(3)
        lims(1,1) = -params_glob%trs
        lims(1,2) =  params_glob%trs
        lims(2,1) = -params_glob%trs
        lims(2,2) =  params_glob%trs
        do ithr = 1, params_glob%nthr
            call grad_shsrch_obj(ithr)%new(lims, opt_angle=params_glob%l_reg_opt_ang)
        enddo
        !$omp parallel do default(shared) private(iref,iptcl,irot,ithr,cxy) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            iptcl = self%ptcl_ref_map(iref)
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                ithr = omp_get_thread_num() + 1
                call grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                irot = self%ref_ptcl_tab(iref,iptcl)%loc
                cxy  = grad_shsrch_obj(ithr)%minimize(irot)
                if( irot > 0 )then
                    self%ref_ptcl_tab(iref,iptcl)%sh  = cxy(2:3)
                    self%ref_ptcl_tab(iref,iptcl)%loc = irot
                endif
            endif
        enddo
        !$omp end parallel do
    end subroutine shift_search

    subroutine nonuni_tab_align( self )
        class(regularizer), intent(inout) :: self
        integer :: ir, min_ind_ir, min_ind_ip, min_ip(self%nrefs)
        real    :: min_ir(self%nrefs)
        logical :: mask_ip(params_glob%fromp:params_glob%top)
        self%ptcl_ref_map = 1   
        mask_ip           = .true.
        call seed_rnd
        do while( any(mask_ip) )
            min_ir = 0.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ir)
            do ir = 1, self%nrefs
                min_ip(ir) = params_glob%fromp + minloc(self%ref_ptcl_cor(ir,:), dim=1, mask=mask_ip) - 1
                min_ir(ir) = self%ref_ptcl_cor(ir,min_ip(ir))
            enddo
            !$omp end parallel do
            min_ind_ir = self%ref_multinomal(min_ir)
            min_ind_ip = min_ip(min_ind_ir)
            self%ptcl_ref_map(min_ind_ip) = min_ind_ir
            mask_ip(min_ind_ip) = .false.
        enddo
    end subroutine nonuni_tab_align

    !>  \brief  generates a multinomal 1-of-K random number according to the
    !!          distribution in pvec
    function ref_multinomal( self, pvec ) result( which )
        class(regularizer), intent(inout) :: self
        real,               intent(in)    :: pvec(:) !< probabilities
        integer :: i, which, ithr
        real    :: rnd, bound
        ithr = omp_get_thread_num() + 1
        self%refs_corr(:,ithr) = pvec
        self%refs_inds(:,ithr) = (/(i,i=1,self%nrefs)/)
        call hpsort(self%refs_corr(:,ithr), self%refs_inds(:,ithr) )
        rnd = ran3()
        if( sum(self%refs_corr(1:self%refs_ns,ithr)) < TINY )then
            ! uniform sampling
            which = 1 + floor(real(self%refs_ns) * rnd)
        else
            ! normalizing within the hard-limit
            self%refs_corr(1:self%refs_ns,ithr) = self%refs_corr(1:self%refs_ns,ithr) / sum(self%refs_corr(1:self%refs_ns,ithr))
            do which=1,self%refs_ns
                bound = sum(self%refs_corr(1:which, ithr))
                if( rnd >= bound )exit
            enddo
            which = self%refs_inds(min(which,self%refs_ns),ithr)
        endif
    end function ref_multinomal

    ! inpl multinomal based on unnormalized pvec
    function inpl_multinomal( self, pvec ) result( which )
        class(regularizer), intent(inout) :: self
        real,               intent(in)    :: pvec(:) !< probabilities
        integer :: i, which, ithr
        real    :: rnd, bound
        ithr = omp_get_thread_num() + 1
        self%inpl_corr(:,ithr) = pvec
        self%inpl_inds(:,ithr) = (/(i,i=1,self%nrots)/)
        call hpsort(self%inpl_corr(:,ithr), self%inpl_inds(:,ithr) )
        rnd = ran3()
        if( sum(self%inpl_corr(1:self%inpl_ns,ithr)) < TINY )then
            ! uniform sampling
            which = 1 + floor(real(self%inpl_ns) * rnd)
        else
            ! normalizing within the hard-limit
            self%inpl_corr(1:self%inpl_ns,ithr) = self%inpl_corr(1:self%inpl_ns,ithr) / sum(self%inpl_corr(1:self%inpl_ns,ithr))
            do which=1,self%inpl_ns
                bound = sum(self%inpl_corr(1:which, ithr))
                if( rnd >= bound )exit
            enddo
            which = self%inpl_inds(min(which,self%inpl_ns),ithr)
        endif
    end function inpl_multinomal

    ! DESTRUCTOR

    subroutine kill( self )
        class(regularizer), intent(inout) :: self
        deallocate(self%ref_ptcl_cor,self%ref_ptcl_tab,self%ptcl_ref_map,self%inpl_corr,self%refs_corr,self%inpl_inds,self%refs_inds)
    end subroutine kill
end module simple_regularizer
