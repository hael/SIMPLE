! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,       only: params_glob
use simple_ori,              only: geodesic_frobdev
use simple_polarft_corrcalc, only: polarft_corrcalc
implicit none

public :: regularizer
private
#include "simple_local_flags.inc"

type :: regularizer
    integer               :: nrots
    integer               :: nrefs
    integer               :: nptcls
    integer               :: pftsz
    integer               :: kfromto(2)
    real(dp), allocatable :: refs_reg_even(:,:,:)        !< -"-, reference reg terms, even
    real(dp), allocatable :: refs_reg_odd(:,:,:)         !< -"-, reference reg terms, odd
    real(dp), allocatable :: refs_reg(:,:,:)             !< -"-, reference reg terms
    real(dp), allocatable :: regs_denom_even(:,:,:)      !< -"-, even
    real(dp), allocatable :: regs_denom_odd(:,:,:)       !< -"-, odd
    real(dp), allocatable :: regs_denom(:,:,:)           !< -"-
    class(polarft_corrcalc), pointer :: pftcc => null()
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: ref_reg_cc_2D
    procedure          :: ref_reg_cc
    procedure          :: ref_reg_cc_neigh
    procedure          :: ref_reg_cc_dev
    procedure          :: regularize_refs
    procedure          :: regularize_refs_2D
    procedure          :: reset_regs
    procedure, private :: calc_raw_frc, calc_pspec
    procedure, private :: rotate_polar_real, rotate_polar_complex
    generic            :: rotate_polar => rotate_polar_real, rotate_polar_complex
    ! DESTRUCTOR
    procedure          :: kill
end type regularizer

contains
    ! CONSTRUCTORS

    subroutine new( self, pftcc )
        class(regularizer),      target, intent(inout) :: self
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        self%nrots   = pftcc%nrots
        self%nrefs   = pftcc%nrefs
        self%nptcls  = pftcc%nptcls
        self%pftsz   = pftcc%pftsz
        self%kfromto = pftcc%kfromto
        ! allocation
        allocate(self%regs_denom_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%refs_reg_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%refs_reg_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%refs_reg(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs))
        self%refs_reg_even   = 0.d0
        self%refs_reg_odd    = 0.d0
        self%refs_reg        = 0.d0
        self%regs_denom_even = 0.d0
        self%regs_denom_odd  = 0.d0
        self%pftcc           => pftcc
    end subroutine new

    ! 2D accumulating reference reg terms for each batch of particles
    subroutine ref_reg_cc_2D( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc
        real     :: inpl_corrs(self%nrots), ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pftcc%pfts_ptcls * self%pftcc%ctfmats)
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%nptcls
                iptcl = glob_pinds(i)
                ! find best irot for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                loc = maxloc(inpl_corrs, dim=1)
                if( inpl_corrs(loc) < TINY ) cycle
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                if( self%pftcc%iseven(i) )then
                    self%refs_reg_even(  :,:,iref) = self%refs_reg_even(  :,:,iref) + ptcl_ctf_rot * real(inpl_corrs(loc), dp)
                    self%regs_denom_even(:,:,iref) = self%regs_denom_even(:,:,iref) + ctf_rot**2
                else
                    self%refs_reg_odd(  :,:,iref) = self%refs_reg_odd(  :,:,iref) + ptcl_ctf_rot * real(inpl_corrs(loc), dp)
                    self%regs_denom_odd(:,:,iref) = self%regs_denom_odd(:,:,iref) + ctf_rot**2
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc_2D

    ! accumulating reference reg terms for each batch of particles, with cc-based global objfunc
    subroutine ref_reg_cc( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(regularizer), intent(inout) :: self
        type(oris),         intent(in)    :: eulspace
        type(oris),         intent(in)    :: ptcl_eulspace
        integer,            intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real     :: euls(3), euls_ref(3)
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pftcc%pfts_ptcls * self%pftcc%ctfmats)
        ! even/odd only when lpset is .false.
        if( params_glob%l_lpset )then
            !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot, ctf_rot) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    iptcl    = glob_pinds(i)
                    euls_ref = eulspace%get_euler(iref)
                    euls     = ptcl_eulspace%get_euler(iptcl)
                    ! projection direction distance, euler_dist could be used instead
                    euls_ref(3)   = 0.
                    euls(3)       = 0.
                    ptcl_ref_dist = geodesic_frobdev(euls_ref,euls)
                    ! find best irot for this pair of iref, iptcl
                    call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                    loc = maxloc(inpl_corrs, dim=1)
                    if( inpl_corrs(loc) < TINY ) cycle
                    ! distance & correlation weighing
                    ptcl_ref_dist = inpl_corrs(loc) / ( 1. + ptcl_ref_dist )
                    ! computing the reg terms as the gradients w.r.t 2D references of the probability
                    loc = (self%nrots+1)-(loc-1)
                    if( loc > self%nrots ) loc = loc - self%nrots
                    call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                    call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                    self%refs_reg(  :,:,iref) = self%refs_reg(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                    self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + ctf_rot**2
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot, ctf_rot) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    iptcl    = glob_pinds(i)
                    euls_ref = eulspace%get_euler(iref)
                    euls     = ptcl_eulspace%get_euler(iptcl)
                    ! projection direction distance, euler_dist could be used instead
                    euls_ref(3)   = 0.
                    euls(3)       = 0.
                    ptcl_ref_dist = geodesic_frobdev(euls_ref,euls)
                    ! find best irot for this pair of iref, iptcl
                    call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                    loc = maxloc(inpl_corrs, dim=1)
                    if( inpl_corrs(loc) < TINY ) cycle
                    ! distance & correlation weighing
                    ptcl_ref_dist = inpl_corrs(loc) / ( 1. + ptcl_ref_dist )
                    ! computing the reg terms as the gradients w.r.t 2D references of the probability
                    loc = (self%nrots+1)-(loc-1)
                    if( loc > self%nrots ) loc = loc - self%nrots
                    call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                    call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                    if( self%pftcc%iseven(i) )then
                        self%refs_reg_even(  :,:,iref) = self%refs_reg_even(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                        self%regs_denom_even(:,:,iref) = self%regs_denom_even(:,:,iref) + ctf_rot**2
                    else
                        self%refs_reg_odd(  :,:,iref) = self%refs_reg_odd(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                        self%regs_denom_odd(:,:,iref) = self%regs_denom_odd(:,:,iref) + ctf_rot**2
                    endif
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine ref_reg_cc

    subroutine ref_reg_cc_neigh( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(regularizer), intent(inout) :: self
        type(oris),         intent(in)    :: eulspace
        type(oris),         intent(in)    :: ptcl_eulspace
        integer,            intent(in)    :: glob_pinds(self%nptcls)
        integer  :: i, iref, iptcl, loc
        real     :: inpl_corrs(self%nrots), ptcl_ref_dist, ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%nptcls)
        real     :: euls(3), euls_ref(3), theta
        real(dp) :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2)), ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = real(self%pftcc%pfts_ptcls * self%pftcc%ctfmats)
        ! even/odd only when lpset is .false.
        if( params_glob%l_lpset )then
            !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot,theta) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    iptcl    = glob_pinds(i)
                    euls_ref =      eulspace%get_euler(iref)  * pi / 180.
                    euls     = ptcl_eulspace%get_euler(iptcl) * pi / 180.
                    theta    = acos(cos(euls_ref(2))*cos(euls(2)) + sin(euls_ref(2))*sin(euls(2))*cos(euls_ref(1) - euls(1)))
                    if( theta <= params_glob%arc_thres*pi/180.  .and. theta >= 0. )then
                        ! find best irot for this pair of iref, iptcl
                        call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                        loc = maxloc(inpl_corrs, dim=1)
                        if( inpl_corrs(loc) < TINY ) cycle
                        ! distance & correlation weighing
                        ptcl_ref_dist = inpl_corrs(loc) / ( 1. + theta )
                        ! computing the reg terms as the gradients w.r.t 2D references of the probability
                        loc = (self%nrots+1)-(loc-1)
                        if( loc > self%nrots ) loc = loc - self%nrots
                        call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                        call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                        self%refs_reg(  :,:,iref) = self%refs_reg(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                        self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + ctf_rot**2
                    endif
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) default(shared) private(i,iref,euls_ref,euls,ptcl_ref_dist,iptcl,inpl_corrs,loc,ptcl_ctf_rot,ctf_rot,theta) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    iptcl    = glob_pinds(i)
                    euls_ref =      eulspace%get_euler(iref)  * pi / 180.
                    euls     = ptcl_eulspace%get_euler(iptcl) * pi / 180.
                    theta    = acos(cos(euls_ref(2))*cos(euls(2)) + sin(euls_ref(2))*sin(euls(2))*cos(euls_ref(1) - euls(1)))
                    if( theta <= params_glob%arc_thres*pi/180.  .and. theta >= 0. )then
                        ! find best irot for this pair of iref, iptcl
                        call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                        loc = maxloc(inpl_corrs, dim=1)
                        if( inpl_corrs(loc) < TINY ) cycle
                        ! distance & correlation weighing
                        ptcl_ref_dist = inpl_corrs(loc) / ( 1. + theta )
                        ! computing the reg terms as the gradients w.r.t 2D references of the probability
                        loc = (self%nrots+1)-(loc-1)
                        if( loc > self%nrots ) loc = loc - self%nrots
                        call self%rotate_polar(          ptcl_ctf(:,:,i), ptcl_ctf_rot, loc)
                        call self%rotate_polar(self%pftcc%ctfmats(:,:,i),      ctf_rot, loc)
                        if( self%pftcc%iseven(i) )then
                            self%refs_reg_even(  :,:,iref) = self%refs_reg_even(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                            self%regs_denom_even(:,:,iref) = self%regs_denom_even(:,:,iref) + ctf_rot**2
                        else
                            self%refs_reg_odd(  :,:,iref) = self%refs_reg_odd(  :,:,iref) + ptcl_ctf_rot * real(ptcl_ref_dist, dp)
                            self%regs_denom_odd(:,:,iref) = self%regs_denom_odd(:,:,iref) + ctf_rot**2
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine ref_reg_cc_neigh

    ! for testing new ideas
    subroutine ref_reg_cc_dev( self, eulspace, ptcl_eulspace, glob_pinds )
        use simple_oris
        class(regularizer), intent(inout) :: self
        type(oris),         intent(in)    :: eulspace
        type(oris),         intent(in)    :: ptcl_eulspace
        integer,            intent(in)    :: glob_pinds(self%nptcls)
    end subroutine ref_reg_cc_dev

    subroutine regularize_refs_2D( self )
        class(regularizer), intent(inout) :: self
        integer :: iref, k
        !$omp parallel default(shared) private(k,iref) proc_bind(close)
        !$omp do schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            where( abs(self%regs_denom_even(:,k,:)) < TINY )
                self%refs_reg_even(:,k,:) = real(k, dp) * self%refs_reg_even(:,k,:)
            elsewhere
                self%refs_reg_even(:,k,:) = real(k, dp) * self%refs_reg_even(:,k,:) / self%regs_denom_even(:,k,:)
            endwhere
            where( abs(self%regs_denom_odd(:,k,:)) < TINY )
                self%refs_reg_odd(:,k,:) = real(k, dp) * self%refs_reg_odd(:,k,:)
            elsewhere
                self%refs_reg_odd(:,k,:) = real(k, dp) * self%refs_reg_odd(:,k,:) / self%regs_denom_odd(:,k,:)
            endwhere
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do iref = 1, self%nrefs
            self%pftcc%pfts_refs_even(:,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_even(:,:,iref) + params_glob%eps * real(self%refs_reg_even(:,:,iref))
            self%pftcc%pfts_refs_odd( :,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_odd( :,:,iref) + params_glob%eps * real(self%refs_reg_odd(:,:,iref))
        enddo
        !$omp end do nowait
        !$omp end parallel
        call self%pftcc%memoize_refs
    end subroutine regularize_refs_2D

    subroutine regularize_refs( self )
        use simple_opt_filter, only: butterworth_filter
        class(regularizer), intent(inout) :: self
        integer :: iref, k, find
        real    :: filt(self%kfromto(1):self%kfromto(2))
        ! even/odd only when lpset is .false.
        if( params_glob%l_lpset )then
            ! generating butterworth filter at cut-off = lp
            find = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
            call butterworth_filter(find, self%kfromto, filt)
            !$omp parallel default(shared) private(k,iref) proc_bind(close)
            !$omp do schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                where( abs(self%regs_denom(:,k,:)) < TINY )
                    self%refs_reg(:,k,:) = real(k, dp) * self%refs_reg(:,k,:)
                elsewhere
                    self%refs_reg(:,k,:) = real(k, dp) * self%refs_reg(:,k,:) / self%regs_denom(:,k,:)
                endwhere
            enddo
            !$omp end do nowait
            !$omp do schedule(static)
            do iref = 1, self%nrefs
                self%pftcc%pfts_refs_even(:,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_even(:,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref))
                self%pftcc%pfts_refs_odd( :,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_odd( :,:,iref) + params_glob%eps * real(self%refs_reg(:,:,iref))
            enddo
            !$omp end do nowait
            !$omp do schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                self%pftcc%pfts_refs_even(:,k,:) = filt(k) * self%pftcc%pfts_refs_even(:,k,:)
                self%pftcc%pfts_refs_odd( :,k,:) = filt(k) * self%pftcc%pfts_refs_odd( :,k,:)
            enddo
            !$omp end do nowait
            !$omp end parallel
        else
            !$omp parallel default(shared) private(k,iref,filt) proc_bind(close)
            !$omp do schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                where( abs(self%regs_denom_even(:,k,:)) < TINY )
                    self%refs_reg_even(:,k,:) = real(k, dp) * self%refs_reg_even(:,k,:)
                elsewhere
                    self%refs_reg_even(:,k,:) = real(k, dp) * self%refs_reg_even(:,k,:) / self%regs_denom_even(:,k,:)
                endwhere
                where( abs(self%regs_denom_odd(:,k,:)) < TINY )
                    self%refs_reg_odd(:,k,:) = real(k, dp) * self%refs_reg_odd(:,k,:)
                elsewhere
                    self%refs_reg_odd(:,k,:) = real(k, dp) * self%refs_reg_odd(:,k,:) / self%regs_denom_odd(:,k,:)
                endwhere
            enddo
            !$omp end do nowait
            !$omp do schedule(static)
            do iref = 1, self%nrefs
                self%pftcc%pfts_refs_even(:,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_even(:,:,iref) + params_glob%eps * real(self%refs_reg_even(:,:,iref))
                self%pftcc%pfts_refs_odd( :,:,iref) = (1. - params_glob%eps) * self%pftcc%pfts_refs_odd( :,:,iref) + params_glob%eps * real(self%refs_reg_odd(:,:,iref))
            enddo
            !$omp end do nowait
            ! applying frc filter
            !$omp do schedule(static)
            do iref = 1, self%nrefs
                call self%calc_raw_frc(self%pftcc%pfts_refs_even(:,:,iref), self%pftcc%pfts_refs_odd(:,:,iref), filt)
                do k = self%kfromto(1),self%kfromto(2)
                    self%pftcc%pfts_refs_even(:,k,iref) = filt(k) * self%pftcc%pfts_refs_even(:,k,iref)
                    self%pftcc%pfts_refs_odd( :,k,iref) = filt(k) * self%pftcc%pfts_refs_odd( :,k,iref)
                enddo
            enddo
            !$omp end do nowait
            !$omp end parallel
        endif
        call self%pftcc%memoize_refs
    end subroutine regularize_refs

    subroutine reset_regs( self )
        class(regularizer), intent(inout) :: self
        self%refs_reg_even   = 0._dp
        self%refs_reg_odd    = 0._dp
        self%refs_reg        = 0._dp
        self%regs_denom_even = 0._dp
        self%regs_denom_odd  = 0._dp
        self%regs_denom      = 0._dp
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
        deallocate(self%refs_reg_even,self%refs_reg_odd,self%regs_denom_even,self%regs_denom_odd,self%refs_reg,self%regs_denom)
    end subroutine kill
end module simple_regularizer