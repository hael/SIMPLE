! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_ori,               only: geodesic_frobdev
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
implicit none

public :: regularizer
private
#include "simple_local_flags.inc"

type :: regularizer
    integer               :: nrots
    integer               :: nrefs
    integer               :: pftsz
    integer               :: kfromto(2)
    real(dp), allocatable :: regs_even(:,:,:)        !< -"-, reg terms, even
    real(dp), allocatable :: regs_odd(:,:,:)         !< -"-, reg terms, odd
    real(dp), allocatable :: regs(:,:,:)             !< -"-, reg terms
    real(dp), allocatable :: regs_neigh(:,:,:)       !< -"-, neighborhood reg terms
    real(dp), allocatable :: regs_denom_even(:,:,:)  !< -"-, reg denom, even
    real(dp), allocatable :: regs_denom_odd(:,:,:)   !< -"-, reg denom, odd
    real(dp), allocatable :: regs_denom(:,:,:)       !< -"-, reg denom
    real(dp), allocatable :: regs_denom_neigh(:,:,:) !< -"-, neighborhood reg denom
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_obj(:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: ref_reg_cc
    procedure          :: regularize_refs
    procedure          :: reset_regs
    procedure, private :: calc_raw_frc, calc_pspec
    procedure, private :: rotate_polar_real, rotate_polar_complex, rotate_polar_test
    generic            :: rotate_polar => rotate_polar_real, rotate_polar_complex, rotate_polar_test
    ! DESTRUCTOR
    procedure          :: kill
end type regularizer

contains
    ! CONSTRUCTORS

    subroutine new( self, pftcc )
        class(regularizer),      target, intent(inout) :: self
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        integer, parameter :: MAXITS = 60
        real    :: lims(2,2), lims_init(2,2)
        integer :: ithr
        self%nrots   = pftcc%nrots
        self%nrefs   = pftcc%nrefs
        self%pftsz   = pftcc%pftsz
        self%kfromto = pftcc%kfromto
        ! allocation
        allocate(self%regs_denom_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom_neigh(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_neigh(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%grad_shsrch_obj(params_glob%nthr),&
                &self%regs(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs))
        self%regs_even        = 0.d0
        self%regs_odd         = 0.d0
        self%regs             = 0.d0
        self%regs_neigh       = 0.d0
        self%regs_denom_even  = 0.d0
        self%regs_denom_odd   = 0.d0
        self%regs_denom_neigh = 0.d0
        self%regs_denom       = 0.d0
        self%pftcc            => pftcc
        lims(:,1)             = -params_glob%reg_minshift
        lims(:,2)             =  params_glob%reg_minshift
        lims_init(:,1)        = -SHC_INPL_TRSHWDTH
        lims_init(:,2)        =  SHC_INPL_TRSHWDTH
        do ithr = 1, params_glob%nthr
            call self%grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init,&
            &shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=params_glob%l_reg_opt_ang)
        enddo
    end subroutine new

    ! accumulating reference reg terms for each batch of particles, with cc-based global objfunc
    subroutine ref_reg_cc( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(regularizer), intent(inout) :: self
        type(oris),         intent(in)    :: eulspace
        type(oris),         intent(in)    :: ptcl_eulspace
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        complex(sp),        pointer       :: shmat(:,:)
        integer  :: i, iref, iptcl, loc, ithr
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls), cur_corr
        real     :: euls(3), euls_ref(3), theta, cxy(3)
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), init_xy(2)
        ptcl_ctf = real(self%pftcc%pfts_ptcls * self%pftcc%ctfmats)
        !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot,theta,cur_corr,init_xy,ithr,shmat,cxy) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                ithr     = omp_get_thread_num() + 1
                iptcl    = glob_pinds(i)
                euls_ref = eulspace%get_euler(iref)
                euls     = ptcl_eulspace%get_euler(iptcl)
                ! projection direction distance, euler_dist could be used instead
                euls_ref(3)   = 0.
                euls(3)       = 0.
                ptcl_ref_dist = geodesic_frobdev(euls_ref,euls)
                ! find best irot/shift for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                loc = maxloc(inpl_corrs, dim=1)
                call self%grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                cxy = self%grad_shsrch_obj(ithr)%minimize(irot=loc)
                if( loc > 0 )then
                    cur_corr = cxy(1)
                    init_xy  = cxy(2:3)
                else
                    loc      = maxloc(inpl_corrs, dim=1)
                    cur_corr = inpl_corrs(loc)
                    init_xy  = 0.
                endif
                if( cur_corr < TINY ) cycle
                ! distance & correlation weighing
                ptcl_ref_dist = cur_corr / ( 1. + exp(ptcl_ref_dist) )
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                shmat => self%pftcc%heap_vars(ithr)%shmat
                call self%pftcc%gen_shmat(ithr, -real(init_xy), shmat)
                ptcl_ctf_rot = ptcl_ctf_rot * shmat
                self%regs(:,:,iref)       = self%regs(:,:,iref)       + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + ctf_rot**2
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc

    subroutine regularize_refs( self, ref_freq_in )
        class(regularizer), intent(inout) :: self
        real, optional,     intent(in)    :: ref_freq_in
        integer :: iref, k
        real    :: ref_freq
        ref_freq = 0.
        if( present(ref_freq_in) ) ref_freq = ref_freq_in
        !$omp parallel default(shared) private(k,iref) proc_bind(close)
        !$omp do schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            where( abs(self%regs_denom(:,k,:)) < TINY )
                self%regs(:,k,:) = 0._dp
            elsewhere
                self%regs(:,k,:) = real(k, dp) * self%regs(:,k,:) / self%regs_denom(:,k,:)
            endwhere
        enddo
        !$omp end do
        !$omp do schedule(static)
        do iref = 1, self%nrefs
            if( ran3() < ref_freq )then
                ! keep the refs
            else
                ! using the reg terms as refs
                self%pftcc%pfts_refs_even(:,:,iref) = real(self%regs(:,:,iref))
                self%pftcc%pfts_refs_odd( :,:,iref) = real(self%regs(:,:,iref))
            endif
        enddo
        !$omp end do
        !$omp end parallel
        call self%pftcc%memoize_refs
    end subroutine regularize_refs
    
    subroutine reset_regs( self )
        class(regularizer), intent(inout) :: self
        self%regs_even        = 0._dp
        self%regs_odd         = 0._dp
        self%regs             = 0._dp
        self%regs_neigh       = 0._dp
        self%regs_denom_even  = 0._dp
        self%regs_denom_odd   = 0._dp
        self%regs_denom_neigh = 0._dp
        self%regs_denom       = 0._dp
    end subroutine reset_regs

    subroutine rotate_polar_real( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(regularizer), intent(inout) :: self
        real(sp),           intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),           intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,            intent(in)    :: irot
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        ! just need the realpart
        if( irot == 1 .or. irot == self%pftsz + 1 )then
            ptcl_ctf_rot = ptcl_ctf
        else
            ptcl_ctf_rot(  1:rot-1    , :) = ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:)
            ptcl_ctf_rot(rot:self%pftsz,:) = ptcl_ctf(               1:self%pftsz-rot+1,:)
        end if
    end subroutine rotate_polar_real

    subroutine rotate_polar_complex( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(regularizer), intent(inout) :: self
        complex(dp),        intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),        intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,            intent(in)    :: irot
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        if( irot == 1 )then
            ptcl_ctf_rot = ptcl_ctf
        else if( irot <= self%pftsz )then
            ptcl_ctf_rot(rot:self%pftsz,:) =       ptcl_ctf(               1:self%pftsz-rot+1,:)
            ptcl_ctf_rot(  1:rot-1     ,:) = conjg(ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:))
        else if( irot == self%pftsz + 1 )then
            ptcl_ctf_rot = conjg(ptcl_ctf)
        else
            ptcl_ctf_rot(rot:self%pftsz,:) = conjg(ptcl_ctf(               1:self%pftsz-rot+1,:))
            ptcl_ctf_rot(  1:rot-1     ,:) =       ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:)
        end if
    end subroutine rotate_polar_complex

    subroutine rotate_polar_test( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(regularizer), intent(inout) :: self
        complex(sp),        intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),           intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,            intent(in)    :: irot
        integer :: rot
        if( irot >= self%pftsz + 1 )then
            rot = irot - self%pftsz
        else
            rot = irot
        end if
        ! just need the realpart
        if( irot == 1 .or. irot == self%pftsz + 1 )then
            ptcl_ctf_rot = real(ptcl_ctf, dp)
        else
            ptcl_ctf_rot(  1:rot-1    , :) = real(ptcl_ctf(self%pftsz-rot+2:self%pftsz      ,:), dp)
            ptcl_ctf_rot(rot:self%pftsz,:) = real(ptcl_ctf(               1:self%pftsz-rot+1,:), dp)
        end if
    end subroutine rotate_polar_test

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
        deallocate(self%regs_even,self%regs_odd,self%regs_denom_even,self%regs_denom_odd,&
                  &self%regs,self%regs_denom,self%regs_neigh,self%regs_denom_neigh)
    end subroutine kill
end module simple_regularizer