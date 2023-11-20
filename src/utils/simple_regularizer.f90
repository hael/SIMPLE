! regularizer of the cluster2D and refine3D
module simple_regularizer
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_builder,           only: build_glob
use simple_ori,               only: geodesic_frobdev
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
implicit none

public :: regularizer
private
#include "simple_local_flags.inc"

type reg_params
    integer :: iptcl            !< iptcl index
    integer :: iref             !< iref index
    integer :: loc              !< inpl index
    real    :: prob, sh(2), w   !< probability, shift, and weight
    real    :: sum
end type reg_params

type :: regularizer
    integer                  :: nrots
    integer                  :: nrefs
    integer                  :: nneighs
    integer                  :: pftsz
    integer                  :: kfromto(2)
    complex(dp), allocatable :: regs_odd(:,:,:)             !< -"-, reg terms
    complex(dp), allocatable :: regs_even(:,:,:)            !< -"-, reg terms
    complex(dp), allocatable :: regs_grad_odd(:,:,:)        !< -"-, reg terms
    complex(dp), allocatable :: regs_grad_even(:,:,:)       !< -"-, reg terms
    real(dp),    allocatable :: regs_denom_odd(:,:,:)       !< -"-, reg denom
    real(dp),    allocatable :: regs_denom_even(:,:,:)      !< -"-, reg denom
    real,        allocatable :: ref_ptcl_corr(:,:)      !< 2D corr table
    integer,     allocatable :: ptcl_ref_map(:)         !< hard-alignment tab
    integer,     allocatable :: ref_neigh_map(:)        !< mapping ref to neighborhood
    integer,     allocatable :: prev_ptcl_ref(:)        !< prev ptcl-ref map
    logical,     allocatable :: ref_neigh_tab(:,:)      !< athres-neighborhood map
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(reg_params),        allocatable :: ref_ptcl_tab(:,:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: init_tab
    procedure          :: fill_tab_noshift
    procedure          :: partition_refs
    procedure          :: make_neigh_tab
    procedure          :: map_ptcl_ref
    procedure          :: uniform_cluster_sort
    procedure          :: find_closest_iref
    procedure          :: uniform_cluster_sort_neigh
    procedure          :: uniform_cluster_sort_dyn
    procedure          :: reg_uniform_cluster
    procedure          :: reg_lap
    procedure          :: prev_cavgs
    procedure          :: form_cavgs
    procedure          :: compute_grad_ptcl
    procedure          :: compute_grad_cavg
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
        self%nrots   = pftcc%nrots
        self%nrefs   = pftcc%nrefs
        self%pftsz   = pftcc%pftsz
        self%kfromto = pftcc%kfromto
        ! allocation
        allocate(self%regs_denom_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_denom_odd( self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_grad_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_grad_odd( self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_odd( self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ref_neigh_tab(self%nrefs,self%nrefs),self%prev_ptcl_ref(params_glob%fromp:params_glob%top))
        self%regs_odd        = 0.d0
        self%regs_even       = 0.d0
        self%regs_grad_odd   = 0.d0
        self%regs_grad_even  = 0.d0
        self%regs_denom_odd  = 0.d0
        self%regs_denom_even = 0.d0
        self%pftcc      => pftcc
        if( params_glob%l_reg_neigh )then
            call self%partition_refs
            call self%make_neigh_tab
        endif
    end subroutine new

    ! setting up ref neigh tab
    subroutine make_neigh_tab( self, athres_in )
        class(regularizer), target, intent(inout) :: self
        real,             optional, intent(in)    :: athres_in
        type(ori) :: o
        real      :: athres
        integer   :: iref, iref2
        logical   :: lnns(self%nrefs)
        self%ref_neigh_tab = .false.
        athres             = params_glob%athres
        if( present(athres_in) ) athres = athres_in
        do iref = 1, self%nrefs
            lnns = .false.
            call build_glob%eulspace%get_ori(iref, o)
            call build_glob%pgrpsyms%nearest_proj_neighbors(build_glob%eulspace, o, athres, lnns)
            do iref2 = 1, self%nrefs
                if( iref2 /= iref .and. lnns(iref2) )then
                    self%ref_neigh_tab(iref,  iref2) = .true.
                    self%ref_neigh_tab(iref2, iref ) = .true.
                endif
            enddo
            self%ref_neigh_tab(iref, iref) = .true.
        enddo
    end subroutine make_neigh_tab

    subroutine init_tab( self )
        class(regularizer), intent(inout) :: self
        integer :: iptcl, iref
        if( .not.(allocated(self%ref_ptcl_corr)) )then
            allocate(self%ref_ptcl_corr(params_glob%fromp:params_glob%top,self%nrefs), source=0.)
            allocate(self%ref_ptcl_tab( self%nrefs,params_glob%fromp:params_glob%top))
            allocate(self%ptcl_ref_map( params_glob%fromp:params_glob%top))
        endif
        do iref = 1,self%nrefs
            do iptcl = params_glob%fromp,params_glob%top
                self%ref_ptcl_tab(iref,iptcl)%iptcl = iptcl
                self%ref_ptcl_tab(iref,iptcl)%iref  = iref
                self%ref_ptcl_tab(iref,iptcl)%loc   = 0
                self%ref_ptcl_tab(iref,iptcl)%prob  = 0.
                self%ref_ptcl_tab(iref,iptcl)%sh    = 0.
                self%ref_ptcl_tab(iref,iptcl)%w     = 0.
                self%ref_ptcl_tab(iref,iptcl)%sum   = 0.
            enddo
        enddo
    end subroutine init_tab

    subroutine fill_tab_noshift( self, glob_pinds )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: glob_pinds(self%pftcc%nptcls)
        integer   :: i, iref, iptcl
        real      :: inpl_corrs(self%nrots)
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,inpl_corrs) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                iptcl = glob_pinds(i)
                ! find best irot/shift for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, inpl_corrs )
                self%ref_ptcl_tab(iref,iptcl)%sh  = 0.
                self%ref_ptcl_tab(iref,iptcl)%loc = maxloc(inpl_corrs, dim=1)
                self%ref_ptcl_corr(iptcl,iref)    = max(0.,inpl_corrs(self%ref_ptcl_tab(iref,iptcl)%loc))
            enddo
        enddo
        !$omp end parallel do
    end subroutine fill_tab_noshift

    subroutine partition_refs( self )
        class(regularizer), intent(inout) :: self
        integer :: i1, i2, iref, sqn, cnt
        real    :: eul2, eul1
        self%nneighs = params_glob%reg_nneighs
        if( self%nneighs > self%nrefs ) THROW_HARD('reg partition_refs: nneighs > nrefs')
        if( .not.(allocated(self%ref_neigh_map)) ) allocate(self%ref_neigh_map(self%nrefs))
        sqn = int(sqrt(real(self%nneighs)))
        cnt = 1
        do i1 = 1, sqn
            do i2 = 1, sqn
                do iref = 1, self%nrefs
                    eul1 = build_glob%eulspace%e1get(iref) * PI / 180.
                    eul2 = build_glob%eulspace%e2get(iref) * PI / 180.
                    if( ( cos(eul2) >= (-1. + 2.*(i1-1.)/sqn) ) .and. &
                        ( cos(eul2) <= (-1. + 2.*(i1-0.)/sqn) ) .and. &
                        ( eul1      >= (2. * (i2-1.) * PI/sqn)) .and. &
                        ( eul1      <= (2. * (i2-0.) * PI/sqn)) )then
                        self%ref_neigh_map(iref) = cnt
                    endif
                enddo
                cnt = cnt + 1
            enddo
        enddo
    end subroutine partition_refs

    subroutine form_cavgs( self, best_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer     :: iptcl, iref, loc, pind_here
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real        :: grad_w
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iptcl = params_glob%fromp, params_glob%top
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                iref = best_ir(iptcl)
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = self%ref_ptcl_tab(iref, iptcl)%loc
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here), kind=dp), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%pftcc%ctfmats(:,:,pind_here),            ctf_rot, loc)
                grad_w = 1./self%ref_ptcl_tab(iref, iptcl)%sum - self%ref_ptcl_tab(iref, iptcl)%w/self%ref_ptcl_tab(iref, iptcl)%sum**2
                if( self%pftcc%ptcl_iseven(iptcl) )then
                    self%regs_even(:,:,iref)       = self%regs_even(:,:,iref)       + ptcl_ctf_rot
                    self%regs_denom_even(:,:,iref) = self%regs_denom_even(:,:,iref) + ctf_rot**2
                    self%regs_grad_even(:,:,iref)  = self%regs_grad_even(:,:,iref)  + ptcl_ctf_rot * grad_w
                else
                    self%regs_odd(:,:,iref)       = self%regs_odd(:,:,iref)       + ptcl_ctf_rot
                    self%regs_denom_odd(:,:,iref) = self%regs_denom_odd(:,:,iref) + ctf_rot**2
                    self%regs_grad_odd(:,:,iref)  = self%regs_grad_odd(:,:,iref)  + ptcl_ctf_rot * grad_w
                endif
            endif
        enddo
    end subroutine form_cavgs

    subroutine compute_grad_ptcl( self, best_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer     :: iptcl, iref, loc, pind_here, iref2
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iref = 1, self%nrefs
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl,loc,pind_here,iref2,ptcl_ctf_rot)
            do iptcl = params_glob%fromp, params_glob%top
                if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                    iref2 = best_ir(iptcl)
                    if( iref2 /= iref )then
                        pind_here = self%pftcc%pinds(iptcl)
                        loc = self%ref_ptcl_tab(iref2, iptcl)%loc
                        loc = (self%nrots+1)-(loc-1)
                        if( loc > self%nrots ) loc = loc - self%nrots
                        call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here), kind=dp), ptcl_ctf_rot, loc)
                        if( self%pftcc%ptcl_iseven(iptcl) )then
                            self%regs_grad_even(:,:,iref) = self%regs_grad_even(:,:,iref) - ptcl_ctf_rot * &
                                &self%ref_ptcl_tab(iref2, iptcl)%w/self%ref_ptcl_tab(iref2, iptcl)%sum**2
                        else
                            self%regs_grad_odd(:,:,iref) = self%regs_grad_odd(:,:,iref) - ptcl_ctf_rot * &
                                &self%ref_ptcl_tab(iref2, iptcl)%w/self%ref_ptcl_tab(iref2, iptcl)%sum**2
                        endif
                    endif
                endif
            enddo
            !$omp end parallel do
        enddo
    end subroutine compute_grad_ptcl

    subroutine compute_grad_cavg( self )
        class(regularizer), intent(inout) :: self
        complex(dp) :: cavgs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                      &cavgs_odd( self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                      &sum_cavgs_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                      &sum_cavgs_odd( self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs)
        integer     :: iref, iref2
        real(dp)    :: cc_odd(self%nrefs), cc_even(self%nrefs), sum_cc_odd(self%nrefs), sum_cc_even(self%nrefs)
        where( abs(self%regs_denom_odd) < TINY )
            cavgs_odd = 0._dp
        elsewhere
            cavgs_odd = self%regs_odd / self%regs_denom_odd
        endwhere
        where( abs(self%regs_denom_even) < TINY )
            cavgs_even = 0._dp
        elsewhere
            cavgs_even = self%regs_even / self%regs_denom_even
        endwhere
        ! computing cc between avgs and refs
        cc_odd  = 0.d0
        cc_even = 0.d0
        !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            cc_odd(iref)  = cc_odd(iref)  + sum(real(self%pftcc%pfts_refs_odd( :,:,iref) * conjg(cavgs_odd( :,:,iref)),dp))
            cc_even(iref) = cc_even(iref) + sum(real(self%pftcc%pfts_refs_even(:,:,iref) * conjg(cavgs_even(:,:,iref)),dp))
        enddo
        !$omp end parallel do
        sum_cc_odd  = 0.d0
        sum_cc_even = 0.d0
        do iref = 1, self%nrefs
            do iref2 = 1, self%nrefs
                sum_cc_odd( iref) = sum_cc_odd( iref) + sum(real(self%pftcc%pfts_refs_odd( :,:,iref2) * conjg(cavgs_odd( :,:,iref)),dp))
                sum_cc_even(iref) = sum_cc_even(iref) + sum(real(self%pftcc%pfts_refs_even(:,:,iref2) * conjg(cavgs_even(:,:,iref)),dp))
            enddo
        enddo
        ! computing the gradient
        sum_cavgs_even = 0.
        sum_cavgs_odd  = 0.
        do iref = 1, self%nrefs
            do iref2 = 1, self%nrefs
                if( iref2 /= iref )then
                    if( sum_cc_even(iref2) > DTINY )then
                        sum_cavgs_even(:,:,iref) = sum_cavgs_even(:,:,iref) + cc_even(iref2) * cavgs_even(:,:,iref2) / sum_cc_even(iref2)**2
                    endif
                    if( sum_cc_odd(iref2) > DTINY )then
                        sum_cavgs_odd( :,:,iref) = sum_cavgs_odd( :,:,iref) + cc_odd( iref2) * cavgs_odd( :,:,iref2) / sum_cc_odd( iref2)**2
                    endif
                endif
            enddo
        enddo
        do iref = 1, self%nrefs
            if( sum_cc_even(iref) > DTINY )then
                self%regs_grad_odd( :,:,iref) = (cc_odd( iref) - sum_cc_odd( iref)) * cavgs_odd( :,:,iref)/sum_cc_odd( iref)**2 - sum_cavgs_odd( :,:,iref)
            endif
            if( sum_cc_odd(iref) > DTINY )then
                self%regs_grad_even(:,:,iref) = (cc_even(iref) - sum_cc_even(iref)) * cavgs_even(:,:,iref)/sum_cc_even(iref)**2 - sum_cavgs_even(:,:,iref)
            endif
        enddo
    end subroutine compute_grad_cavg

    subroutine prev_cavgs( self )
        class(regularizer), intent(inout) :: self
        type(ori)   :: o_prev
        integer     :: iptcl, iref, loc, pind_here
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iptcl = params_glob%fromp, params_glob%top
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                call build_glob%spproj_field%get_ori(iptcl, o_prev)     ! previous ori
                iref  = build_glob%eulspace%find_closest_proj(o_prev)   ! previous projection direction
                self%prev_ptcl_ref(iptcl) = iref
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = self%pftcc%get_roind(360.-o_prev%e3get())
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here), kind=dp), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%pftcc%ctfmats(:,:,pind_here), ctf_rot, loc)
                if( self%pftcc%ptcl_iseven(iptcl) )then
                    self%regs_even(:,:,iref)       = self%regs_even(:,:,iref)       + ptcl_ctf_rot
                    self%regs_denom_even(:,:,iref) = self%regs_denom_even(:,:,iref) + ctf_rot**2
                else
                    self%regs_odd(:,:,iref)       = self%regs_odd(:,:,iref)       + ptcl_ctf_rot
                    self%regs_denom_odd(:,:,iref) = self%regs_denom_odd(:,:,iref) + ctf_rot**2
                endif
            endif
        enddo
    end subroutine prev_cavgs

    subroutine map_ptcl_ref( self, best_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer :: iptcl
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)
        do iptcl = params_glob%fromp, params_glob%top
            self%ptcl_ref_map(iptcl) = best_ir(iptcl)
        enddo
        !$omp end parallel do
    end subroutine map_ptcl_ref

    subroutine reg_uniform_cluster( self, out_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(inout) :: out_ir(params_glob%fromp:params_glob%top)
        integer :: iref, iptcl
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            sum_corr = sum(self%ref_ptcl_corr(iptcl,:))
            if( sum_corr < TINY )then
                self%ref_ptcl_tab(:,iptcl)%sum = 1.
                self%ref_ptcl_tab(:,iptcl)%w   = 0.
                self%ref_ptcl_corr(iptcl,:)    = 0.
            else
                self%ref_ptcl_tab(:,iptcl)%sum = sum_corr
                self%ref_ptcl_tab(:,iptcl)%w   = self%ref_ptcl_corr(iptcl,:)
                self%ref_ptcl_corr(iptcl,:)    = self%ref_ptcl_corr(iptcl,:) / sum_corr
            endif
        enddo
        !$omp end parallel do
        self%ref_ptcl_corr = self%ref_ptcl_corr / maxval(self%ref_ptcl_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(iref,iptcl)
        do iref = 1, self%nrefs
            do iptcl = params_glob%fromp,params_glob%top
                self%ref_ptcl_tab(iref,iptcl)%prob = self%ref_ptcl_corr(iptcl,iref)
            enddo
        enddo
        !$omp end parallel do
        ! sorted clustering
        out_ir = 1
        if( params_glob%l_reg_neigh )then
            call self%uniform_cluster_sort_dyn(self%nrefs, out_ir)
        else
            call self%uniform_cluster_sort(self%nrefs, out_ir)
        endif
    end subroutine reg_uniform_cluster

    subroutine reg_lap( self, ncols, cur_ir )
        use simple_lap, only: lap
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: ncols
        integer,            intent(inout) :: cur_ir(params_glob%fromp:params_glob%top)
        type(lap) :: lap_obj
        integer   :: sol(params_glob%fromp:params_glob%top), iptcl, iref, n_dup, idup
        real      :: mat(self%nrefs,params_glob%fromp:params_glob%top)
        real      :: cost_mat(params_glob%fromp:params_glob%top,params_glob%fromp:params_glob%top)
        n_dup = ((params_glob%top - params_glob%fromp) + 1)/self%nrefs
        do iref = 1, self%nrefs
            do iptcl = params_glob%fromp, params_glob%top
                mat(iref, iptcl) = self%ref_ptcl_tab(iref,iptcl)%prob
            enddo
        enddo
        do idup = 1, n_dup
            cost_mat((idup-1)*self%nrefs+1:idup*self%nrefs,:) = mat
        enddo
        call lap_obj%new(cost_mat)
        call lap_obj%solve_lap(sol)
        do iptcl = params_glob%fromp, params_glob%top
            cur_ir(iptcl) = mod(sol(iptcl)-1, self%nrefs) + 1
        enddo
    end subroutine reg_lap

    subroutine uniform_cluster_sort( self, ncols, cur_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: ncols
        integer,            intent(inout) :: cur_ir(params_glob%fromp:params_glob%top)
        integer   :: ir, ip, max_ind_ir, max_ind_ip, max_ip(ncols), next_ir, back_ir
        real      :: max_ir(ncols)
        logical   :: mask_ir(ncols), mask_ip(params_glob%fromp:params_glob%top)
        mask_ir = .true.
        mask_ip = .true.
        max_ir  = -1.
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ir,ip)
        do ir = 1, ncols
            if( mask_ir(ir) )then
                do ip = params_glob%fromp, params_glob%top
                    if( mask_ip(ip) .and. self%ref_ptcl_tab(ir, ip)%prob > max_ir(ir) )then
                        max_ir(ir) = self%ref_ptcl_tab(ir, ip)%prob
                        max_ip(ir) = ip
                    endif
                enddo
            endif
        enddo
        !$omp end parallel do
        max_ind_ir = maxloc(max_ir, dim=1, mask=mask_ir)
        max_ind_ip = max_ip(max_ind_ir)
        cur_ir( max_ind_ip) = max_ind_ir
        mask_ip(max_ind_ip) = .false.
        mask_ir(max_ind_ir) = .false.
        back_ir = max_ind_ir
        next_ir = self%find_closest_iref(max_ind_ir, mask_ir)
        do
            if( .not.(any(mask_ip)) ) return
            if( .not.(any(mask_ir)) )then
                mask_ir = .true.
                next_ir = back_ir
            endif
            max_ir(next_ir) = -1
            do ip = params_glob%fromp, params_glob%top
                if( mask_ip(ip) .and. self%ref_ptcl_tab(next_ir, ip)%prob > max_ir(next_ir) )then
                    max_ir(next_ir) = self%ref_ptcl_tab(next_ir, ip)%prob
                    max_ind_ip = ip
                endif
            enddo
            cur_ir( max_ind_ip) = next_ir
            mask_ip(max_ind_ip) = .false.
            mask_ir(next_ir)    = .false.
            next_ir = self%find_closest_iref(next_ir, mask_ir)
        enddo
    end subroutine uniform_cluster_sort

    function find_closest_iref( self, iref, mask_ir ) result( closest )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: iref
        logical,            intent(in)    :: mask_ir(self%nrefs)
        real    :: dists(self%nrefs)
        integer :: closest, i
        type(ori) :: oi, oiref
        call build_glob%eulspace%get_ori(iref, oiref)
        dists = huge(dists(1))
        do i = 1, self%nrefs
            if( i /= iref .and. mask_ir(i) )then
                call build_glob%eulspace%get_ori(i, oi)
                dists(i) = oi.euldist.oiref
            endif
        end do
        closest = minloc( dists, dim=1 )
    end function find_closest_iref

    subroutine uniform_cluster_sort_neigh( self, ncols, cur_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: ncols
        integer,            intent(inout) :: cur_ir(params_glob%fromp:params_glob%top)
        logical,            allocatable   :: mask_neigh(:)
        integer :: ir, ip, max_ind_ir, max_ind_ip, max_ip(ncols)
        real    :: max_ir(ncols)
        logical :: mask_ir(ncols), mask_ip(params_glob%fromp:params_glob%top)
        allocate(mask_neigh(self%nneighs), source=.false.)
        mask_ip = .true.
        mask_ir = .false.
        do
            if( .not.(any(mask_ip)) )    return
            if( .not.(any(mask_neigh)) ) mask_neigh = .true.
            if( .not.(any(mask_ir)) )    mask_ir    = .true.
            max_ir = -1.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ir,ip)
            do ir = 1, ncols
                if( mask_neigh(self%ref_neigh_map(ir)) .and. mask_ir(ir) )then
                    do ip = params_glob%fromp, params_glob%top
                        if( mask_ip(ip) .and. self%ref_ptcl_tab(ir, ip)%prob > max_ir(ir) )then
                            max_ir(ir) = self%ref_ptcl_tab(ir, ip)%prob
                            max_ip(ir) = ip
                        endif
                    enddo
                endif
            enddo
            !$omp end parallel do
            max_ind_ir = maxloc(max_ir, dim=1, mask=mask_ir)
            max_ind_ip = max_ip(max_ind_ir)
            cur_ir( max_ind_ip) = max_ind_ir
            mask_ip(max_ind_ip) = .false.
            mask_ir(max_ind_ir) = .false.
            mask_neigh(self%ref_neigh_map(max_ind_ir)) = .false.
        enddo
    end subroutine uniform_cluster_sort_neigh

    subroutine uniform_cluster_sort_dyn( self, ncols, cur_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: ncols
        integer,            intent(inout) :: cur_ir(params_glob%fromp:params_glob%top)
        integer :: ir, ip, max_ind_ir, max_ind_ip, max_ip(ncols), iref
        real    :: max_ir(ncols)
        logical :: mask_ir(ncols), mask_neigh(ncols), mask_ip(params_glob%fromp:params_glob%top)
        mask_ip    = .true.
        mask_ir    = .false.
        mask_neigh = .false.
        do
            if( .not.(any(mask_ip)) )    return
            if( .not.(any(mask_neigh)) ) mask_neigh = .true.
            if( .not.(any(mask_ir)) )    mask_ir    = .true.
            max_ir = -1.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ir,ip)
            do ir = 1, ncols
                if( mask_neigh(ir) .and. mask_ir(ir) )then
                    do ip = params_glob%fromp, params_glob%top
                        if( mask_ip(ip) .and. self%ref_ptcl_tab(ir, ip)%prob > max_ir(ir) )then
                            max_ir(ir) = self%ref_ptcl_tab(ir, ip)%prob
                            max_ip(ir) = ip
                        endif
                    enddo
                endif
            enddo
            !$omp end parallel do
            max_ind_ir = maxloc(max_ir, dim=1, mask=mask_ir)
            max_ind_ip = max_ip(max_ind_ir)
            cur_ir( max_ind_ip) = max_ind_ir
            mask_ip(max_ind_ip) = .false.
            mask_ir(max_ind_ir) = .false.
            mask_neigh(max_ind_ir) = .false.
            ! flag all the neighbors of max_ind_ir
            do iref = 1, self%nrefs
                if( self%ref_neigh_tab(max_ind_ir, iref) ) mask_neigh(iref) = .false.
            enddo
        enddo
    end subroutine uniform_cluster_sort_dyn


    subroutine regularize_refs( self )
        use simple_image
        use simple_opt_filter, only: butterworth_filter
        class(regularizer), intent(inout) :: self
        complex,            allocatable   :: cmat(:,:)
        type(image) :: calc_cavg
        integer :: iref, k, box, find
        real    :: eps, filt(self%kfromto(1):self%kfromto(2))
        ! form the cavgs
        where( abs(self%regs_denom_odd) < TINY )
            self%regs_odd = 0._dp
        elsewhere
            self%regs_odd = self%regs_odd / self%regs_denom_odd
        endwhere
        where( abs(self%regs_denom_even) < TINY )
            self%regs_even = 0._dp
        elsewhere
            self%regs_even = self%regs_even / self%regs_denom_even
        endwhere
        ! output images for debugging
        if( params_glob%l_reg_debug )then
            do iref = 1, self%nrefs
                ! odd
                call self%pftcc%polar2cartesian(cmplx(self%regs_odd(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%new([box,box,1], params_glob%smpd * real(params_glob%box)/real(box))
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('odd_polar_cavg_reg_'//int2str(params_glob%which_iter)//'.mrc', iref)
                call self%pftcc%polar2cartesian(cmplx(self%pftcc%pfts_refs_odd(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('odd_polar_cavg_'//int2str(params_glob%which_iter)//'.mrc', iref)
                !even
                call self%pftcc%polar2cartesian(cmplx(self%regs_even(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%new([box,box,1], params_glob%smpd * real(params_glob%box)/real(box))
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('even_polar_cavg_reg_'//int2str(params_glob%which_iter)//'.mrc', iref)
                call self%pftcc%polar2cartesian(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=sp), cmat, box)
                call calc_cavg%zero_and_flag_ft
                call calc_cavg%set_cmat(cmat)
                call calc_cavg%shift_phorig()
                call calc_cavg%ifft
                call calc_cavg%write('even_polar_cavg_'//int2str(params_glob%which_iter)//'.mrc', iref)
            enddo
        endif
        ! k-weight
        !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            self%regs_odd( :,k,:)      = real(k) * self%regs_odd( :,k,:)
            self%regs_even(:,k,:)      = real(k) * self%regs_even(:,k,:)
            self%regs_grad_odd( :,k,:) = real(k) * self%regs_grad_odd( :,k,:)
            self%regs_grad_even(:,k,:) = real(k) * self%regs_grad_even(:,k,:)
        enddo
        !$omp end parallel do

        ! applying butterworth filter at cut-off = lp
        find = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
        call butterworth_filter(find, self%kfromto, filt)
        !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            self%regs_odd( :,k,:)      = filt(k) * self%regs_odd( :,k,:)
            self%regs_even(:,k,:)      = filt(k) * self%regs_even(:,k,:)
            self%regs_grad_odd( :,k,:) = filt(k) * self%regs_grad_odd( :,k,:)
            self%regs_grad_even(:,k,:) = filt(k) * self%regs_grad_even(:,k,:)
        enddo
        !$omp end parallel do

        ! taking the real part only (since the global cost function takes only real part)
        self%regs_odd       = real(self%regs_odd,  dp)
        self%regs_even      = real(self%regs_even, dp)
        self%regs_grad_odd  = real(self%regs_grad_odd,  dp)
        self%regs_grad_even = real(self%regs_grad_even, dp)

        ! annealing and different grad styles
        if( params_glob%l_reg_anneal ) eps = min(1., real(params_glob%which_iter) / real(params_glob%reg_iters))
        if( params_glob%l_reg_anneal )then
            if( params_glob%l_reg_grad )then
                !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
                do iref = 1, self%nrefs
                    self%pftcc%pfts_refs_even(:,:,iref) = eps * self%pftcc%pfts_refs_even(:,:,iref) + (1. - eps) * self%regs_grad_even(:,:,iref)
                    self%pftcc%pfts_refs_odd( :,:,iref) = eps * self%pftcc%pfts_refs_odd( :,:,iref) + (1. - eps) * self%regs_grad_odd( :,:,iref)
                enddo
                !$omp end parallel do
            else
                !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
                do iref = 1, self%nrefs
                    self%pftcc%pfts_refs_even(:,:,iref) = eps * self%pftcc%pfts_refs_even(:,:,iref) + (1. - eps) * self%regs_even(:,:,iref)
                    self%pftcc%pfts_refs_odd( :,:,iref) = eps * self%pftcc%pfts_refs_odd( :,:,iref) + (1. - eps) * self%regs_odd( :,:,iref)
                enddo
                !$omp end parallel do
            endif
        else
            if( params_glob%l_reg_grad )then
                !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
                do iref = 1, self%nrefs
                    self%pftcc%pfts_refs_even(:,:,iref) = self%pftcc%pfts_refs_even(:,:,iref) + self%regs_grad_even(:,:,iref)
                    self%pftcc%pfts_refs_odd( :,:,iref) = self%pftcc%pfts_refs_odd( :,:,iref) + self%regs_grad_odd( :,:,iref)
                enddo
                !$omp end parallel do
            else
                !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
                do iref = 1, self%nrefs
                    self%pftcc%pfts_refs_even(:,:,iref) = self%pftcc%pfts_refs_even(:,:,iref) + self%regs_even(:,:,iref)
                    self%pftcc%pfts_refs_odd( :,:,iref) = self%pftcc%pfts_refs_odd( :,:,iref) + self%regs_odd( :,:,iref)
                enddo
                !$omp end parallel do
            endif
        endif
        call self%pftcc%memoize_refs
        call calc_cavg%kill
    end subroutine regularize_refs
    
    subroutine reset_regs( self )
        class(regularizer), intent(inout) :: self
        self%regs_odd        = 0._dp
        self%regs_even       = 0._dp
        self%regs_grad_odd   = 0._dp
        self%regs_grad_even  = 0._dp
        self%regs_denom_odd  = 0._dp
        self%regs_denom_even = 0._dp
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
        deallocate(self%regs_odd, self%regs_denom_odd, self%regs_grad_odd,self%ref_neigh_tab,&
                  &self%regs_even,self%regs_denom_even,self%regs_grad_even,self%prev_ptcl_ref)
        if(allocated(self%ref_neigh_map)) deallocate(self%ref_neigh_map)
        if(allocated(self%ref_ptcl_corr)) deallocate(self%ref_ptcl_corr,self%ref_ptcl_tab,self%ptcl_ref_map)
    end subroutine kill
end module simple_regularizer
