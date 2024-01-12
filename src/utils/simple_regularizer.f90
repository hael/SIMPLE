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
    real,    allocatable :: ref_ptcl_cor(:,:)           !< 2D corr table
    integer, allocatable :: ptcl_ref_map(:)             !< hard-alignment tab
    real,    allocatable :: irot_corr(:,:), refs_corr(:,:)
    integer, allocatable :: irot_inds(:,:), refs_inds(:,:)
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
    procedure, private :: calc_raw_frc, calc_pspec, ref_multinomal, inpl_multinomal
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
        self%pftcc => pftcc
        allocate(self%ref_ptcl_cor(self%nrefs,params_glob%fromp:params_glob%top),&
                &self%refs_corr(self%nrefs,params_glob%nthr), self%irot_corr(self%nrots,params_glob%nthr), source=0.)
        allocate(self%irot_inds(self%nrots,params_glob%nthr), self%refs_inds(self%nrefs,params_glob%nthr), source=0)
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
        integer :: i, j, iref, iptcl, irnd, ithr
        real    :: inpl_corrs(self%nrots), rnd_num, inpl_corrs_norm(self%nrots)
        call seed_rnd
        !$omp parallel do collapse(2) default(shared) private(i,j,iref,iptcl,inpl_corrs,inpl_corrs_norm,irnd,rnd_num,ithr) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                ithr  = omp_get_thread_num() + 1
                iptcl = glob_pinds(i)
                ! find best irot/shift for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                if( maxval(inpl_corrs) < TINY )then
                    self%irot_inds(:,ithr) = (/(j,j=1,self%nrots)/)
                    call hpsort(inpl_corrs, self%irot_inds(:,ithr))
                    call random_number(rnd_num)
                    irnd = 1 + floor(real(self%nrots) * rnd_num)
                    self%ref_ptcl_tab(iref,iptcl)%sh  = 0.
                    self%ref_ptcl_tab(iref,iptcl)%loc = self%irot_inds(irnd, ithr)
                    self%ref_ptcl_cor(iref,iptcl)     =     inpl_corrs(irnd)
                else
                    inpl_corrs_norm = 1. - inpl_corrs/maxval(inpl_corrs)
                    inpl_corrs_norm = inpl_corrs_norm/sum(inpl_corrs_norm)
                    irnd            = self%inpl_multinomal(inpl_corrs_norm)
                    self%ref_ptcl_tab(iref,iptcl)%sh  = 0.
                    self%ref_ptcl_tab(iref,iptcl)%loc = irnd
                    self%ref_ptcl_cor(iref,iptcl)     = inpl_corrs(irnd)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_tab_inpl_smpl

    subroutine tab_normalize( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, iptcl
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        if( params_glob%l_reg_norm )then
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,sum_corr)
            do iptcl = params_glob%fromp, params_glob%top
                sum_corr = sum(self%ref_ptcl_cor(:,iptcl))
                if( sum_corr < TINY )then
                    self%ref_ptcl_cor(:,iptcl) = 0.
                else
                    self%ref_ptcl_cor(:,iptcl) = self%ref_ptcl_cor(:,iptcl) / sum_corr
                endif
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
            min_ir     = 1. - min_ir / maxval(min_ir)
            min_ir     = min_ir / sum(min_ir)
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
        do which=self%nrefs,1,-1
            bound = sum(self%refs_corr(which:self%nrefs, ithr))
            if( rnd <= bound )exit
        enddo
        which = self%refs_inds(max(which,1),ithr)
    end function ref_multinomal

    !>  \brief  generates a multinomal 1-of-K random number according to the
    !!          distribution in pvec
    function inpl_multinomal( self, pvec ) result( which )
        class(regularizer), intent(inout) :: self
        real,               intent(in)    :: pvec(:) !< probabilities
        integer :: i, which, ithr
        real    :: rnd, bound
        ithr = omp_get_thread_num() + 1
        self%irot_corr(:,ithr) = pvec
        self%irot_inds(:,ithr) = (/(i,i=1,self%nrots)/)
        call hpsort(self%irot_corr(:,ithr), self%irot_inds(:,ithr) )
        rnd = ran3()
        do which=self%nrots,1,-1
            bound = sum(self%irot_corr(which:self%nrots, ithr))
            if( rnd <= bound )exit
        enddo
        which = self%irot_inds(max(which,1),ithr)
    end function inpl_multinomal

    ! Calculates frc between two PFTs, rotation, shift & ctf are not factored in
    subroutine calc_raw_frc( self, pft1, pft2, frc )
        class(regularizer), intent(inout) :: self
        complex(sp),        intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(sp),        intent(in)    :: pft2(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,               intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        real(dp) :: num, denom
        integer  :: k
        do k = self%kfromto(1),self%kfromto(2)
            num   = real(sum(pft1(:,k)*conjg(pft2(:,k))),dp)
            denom = real(sum(pft1(:,k)*conjg(pft1(:,k))),dp) * real(sum(pft2(:,k)*conjg(pft2(:,k))),dp)
            if( denom > DTINY )then
                frc(k) = real(num / dsqrt(denom))
            else
                frc(k) = 0.0
            endif
        end do
    end subroutine calc_raw_frc

    ! Calculates normalized PFT power spectrum
    subroutine calc_pspec( self, pft, pspec )
        class(regularizer), intent(inout) :: self
        complex(dp),        intent(in)    :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,               intent(out)   :: pspec(self%kfromto(1):self%kfromto(2))
        integer :: k
        do k = self%kfromto(1),self%kfromto(2)
            pspec(k) = real( real(sum(pft(:,k)*conjg(pft(:,k))),dp) / real(self%pftsz,dp) )
        end do
    end subroutine calc_pspec

    ! DESTRUCTOR

    subroutine kill( self )
        class(regularizer), intent(inout) :: self
        deallocate(self%ref_ptcl_cor,self%ref_ptcl_tab,self%ptcl_ref_map,self%irot_corr,self%refs_corr,self%irot_inds,self%refs_inds)
    end subroutine kill
end module simple_regularizer
