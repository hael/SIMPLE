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
    integer                  :: nrots
    integer                  :: nrefs
    integer                  :: pftsz
    integer                  :: kfromto(2)
    complex(dp), allocatable :: regs(:,:,:)             !< -"-, reg terms
    real(dp),    allocatable :: regs_denom(:,:,:)       !< -"-, reg denom
    real,        allocatable :: ref_corr(:)             !< total ref corr sum
    real,        allocatable :: ref_ptcl_prob(:,:)      !< 2D corr/prob table
    integer,     allocatable :: ref_ptcl_loc(:,:)       !< 2D in-plane table
    real,        allocatable :: ref_ptcl_sh(:,:,:)      !< 2D sh table
    integer,     allocatable :: ref_ptcl_ind(:,:)       !< 2D index table
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_obj(:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: init_tab
    procedure          :: fill_tab
    procedure          :: sort_tab
    procedure          :: ref_reg_cc_tab
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
        allocate(self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%grad_shsrch_obj(params_glob%nthr),self%ref_corr(self%nrefs),&
                &self%regs(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs))
        self%regs       = 0.d0
        self%regs_denom = 0.d0
        self%ref_corr   = 0.
        self%pftcc      => pftcc
        lims(:,1)       = -params_glob%reg_minshift
        lims(:,2)       =  params_glob%reg_minshift
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        do ithr = 1, params_glob%nthr
            call self%grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init,&
            &shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=params_glob%l_reg_opt_ang)
        enddo
    end subroutine new

    subroutine init_tab( self, pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: pinds(params_glob%fromp,params_glob%top)
        allocate(self%ref_ptcl_prob(params_glob%fromp:params_glob%top,self%nrefs),self%ref_ptcl_sh(2,params_glob%fromp:params_glob%top,self%nrefs), source=0.)
        allocate(self%ref_ptcl_loc( params_glob%fromp:params_glob%top,self%nrefs), self%ref_ptcl_ind(params_glob%fromp:params_glob%top,self%nrefs), source=0)
    end subroutine init_tab

    ! filling prob/corr 2D table
    subroutine fill_tab( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        integer  :: i, iref, iptcl, ithr
        real     :: inpl_corrs(self%nrots), cxy(3)
        !$omp parallel do collapse(2) default(shared) private(i,iref,ithr,iptcl,inpl_corrs,cxy) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                ithr  = omp_get_thread_num() + 1
                iptcl = glob_pinds(i)
                ! find best irot/shift for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                self%ref_ptcl_loc(iptcl,iref) = maxloc(inpl_corrs, dim=1)
                call self%grad_shsrch_obj(ithr)%set_indices(iref, iptcl)
                cxy = self%grad_shsrch_obj(ithr)%minimize(irot=self%ref_ptcl_loc(iptcl,iref))
                if( self%ref_ptcl_loc(iptcl,iref) > 0 )then
                    self%ref_ptcl_prob(iptcl,iref) = cxy(1)
                    self%ref_ptcl_sh(:,iptcl,iref) = cxy(2:3)
                else
                    self%ref_ptcl_loc(iptcl,iref)  = maxloc(inpl_corrs, dim=1)
                    self%ref_ptcl_prob(iptcl,iref) = inpl_corrs(self%ref_ptcl_loc(iptcl,iref))
                    self%ref_ptcl_sh(:,iptcl,iref) = 0.
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_tab

    subroutine sort_tab( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, j, iptcl
        real    :: sum_prob
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_prob)
        do iptcl = params_glob%fromp, params_glob%top
            sum_prob = sum(self%ref_ptcl_prob(iptcl,:))
            if( sum_prob < TINY )then
                self%ref_ptcl_prob(iptcl,:) = 0.
            else
                self%ref_ptcl_prob(iptcl,:) = self%ref_ptcl_prob(iptcl,:) / sum_prob
            endif
        enddo
        !$omp end parallel do
        ! sorting the normalized prob for each iref, to sample only the best #ptcls/#refs ptcls for each iref
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,j)
        do iref = 1, self%nrefs
            self%ref_ptcl_ind(:,iref) = (/(j,j=params_glob%fromp,params_glob%top)/)
            call hpsort(self%ref_ptcl_prob(:,iref), self%ref_ptcl_ind(:,iref))
            call reverse(self%ref_ptcl_ind( :,iref))
            call reverse(self%ref_ptcl_prob(:,iref))
        enddo
        !$omp end parallel do
    end subroutine sort_tab

    ! accumulating reference reg terms for each batch of particles, with cc-based global objfunc
    subroutine ref_reg_cc_tab( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        complex(sp),        pointer       :: shmat(:,:)
        integer     :: i, iptcl, iref, ithr, ninds, loc, pind_here
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        real        :: weight
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        ninds    = size(self%ref_ptcl_prob, 1)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2)&
        !$omp private(iref,ithr,i,iptcl,loc,ptcl_ctf_rot,ctf_rot,shmat,pind_here,weight)
        do iref = 1, self%nrefs
            ! taking top sorted corrs/probs
            do i = params_glob%fromp,(params_glob%fromp + int(ninds / self%nrefs))
                if( self%ref_ptcl_prob(i, iref) < TINY ) cycle
                ithr  = omp_get_thread_num() + 1
                iptcl = self%ref_ptcl_ind(i, iref)
                if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                    pind_here = self%pftcc%pinds(iptcl)
                    ! computing the reg terms as the gradients w.r.t 2D references of the probability
                    loc = self%ref_ptcl_loc(iptcl, iref)
                    loc = (self%nrots+1)-(loc-1)
                    if( loc > self%nrots ) loc = loc - self%nrots
                    shmat => self%pftcc%heap_vars(ithr)%shmat
                    call self%pftcc%gen_shmat(ithr, -real(self%ref_ptcl_sh(:,iptcl,iref)), shmat)
                    call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp), ptcl_ctf_rot, loc)
                    call self%rotate_polar(self%pftcc%ctfmats(:,:,pind_here),                    ctf_rot, loc)
                    weight = self%ref_ptcl_prob(i, iref)
                    self%regs(:,:,iref)       = self%regs(:,:,iref)       + weight * ptcl_ctf_rot
                    self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + weight * ctf_rot**2
                    self%ref_corr(iref)       = self%ref_corr(iref)       + self%ref_ptcl_prob(i, iref)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc_tab

    subroutine regularize_refs( self, ref_freq_in )
        use simple_image
        class(regularizer), intent(inout) :: self
        real, optional,     intent(in)    :: ref_freq_in
        real,               parameter     :: REF_FRAC = 1
        integer,            allocatable   :: ref_ind(:)
        complex,            allocatable   :: cmat(:,:)
        type(image) :: calc_cavg
        integer :: iref, k, box
        real    :: ref_freq
        ref_freq = 0.
        if( present(ref_freq_in) ) ref_freq = ref_freq_in
        !$omp parallel default(shared) private(k) proc_bind(close)
        !$omp do schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            where( abs(self%regs_denom(:,k,:)) < TINY )
                self%regs(:,k,:) = 0._dp
            elsewhere
                self%regs(:,k,:) = self%regs(:,k,:) / self%regs_denom(:,k,:)
            endwhere
        enddo
        !$omp end do
        !$omp end parallel
        ! sort ref_corr to only change refs to regs for high-score cavgs
        ref_ind = (/(iref,iref=1,self%nrefs)/)
        call hpsort(self%ref_corr, ref_ind)
        call reverse(ref_ind)
        ! output images for debugging
        if( params_glob%l_reg_debug )then
            do k = 1, int(self%nrefs * REF_FRAC)
                iref = ref_ind(k)
                call self%pftcc%polar2cartesian(cmplx(self%regs(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%new([box,box,1], params_glob%smpd * real(params_glob%box)/real(box))
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('polar_cavg_reg_'//int2str(params_glob%which_iter)//'.mrc', k)
                call self%pftcc%polar2cartesian(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('polar_cavg_'//int2str(params_glob%which_iter)//'.mrc', k)
            enddo
        endif
        !$omp parallel default(shared) private(k,iref) proc_bind(close)
        !$omp do schedule(static)
        do k = 1, int(self%nrefs * REF_FRAC)
            iref = ref_ind(k)
            if( ran3() < ref_freq )then
                ! keep the refs
            else
                ! using the reg terms as refs
                self%pftcc%pfts_refs_even(:,:,iref) = self%regs(:,:,iref)
                self%pftcc%pfts_refs_odd( :,:,iref) = self%regs(:,:,iref)
            endif
        enddo
        !$omp end do
        !$omp end parallel
        call self%pftcc%memoize_refs
        call calc_cavg%kill
    end subroutine regularize_refs
    
    subroutine reset_regs( self )
        class(regularizer), intent(inout) :: self
        self%regs       = 0._dp
        self%regs_denom = 0._dp
        self%ref_corr   = 0.
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
        real(dp),           intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
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
        deallocate(self%regs,self%regs_denom,self%grad_shsrch_obj,self%ref_corr)
        if(allocated(self%ref_ptcl_prob)) deallocate(self%ref_ptcl_prob,self%ref_ptcl_loc,self%ref_ptcl_sh,self%ref_ptcl_ind)
    end subroutine kill
end module simple_regularizer
