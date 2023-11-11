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
    complex(dp), allocatable :: regs(:,:,:)             !< -"-, reg terms
    complex(dp), allocatable :: regs_grad(:,:,:)        !< -"-, reg terms
    real(dp),    allocatable :: regs_denom(:,:,:)       !< -"-, reg denom
    real,        allocatable :: ref_ptcl_corr(:,:)      !< 2D corr table
    integer,     allocatable :: ptcl_ref_map(:)         !< hard-alignment tab
    integer,     allocatable :: ref_neigh_map(:)        !< mapping ref to neighborhood
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
    procedure          :: uniform_cluster_sort_neigh
    procedure          :: uniform_cluster_sort_dyn
    procedure          :: reg_uniform_cluster
    procedure          :: prev_cavgs
    procedure          :: form_cavgs
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
        allocate(self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs_grad(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%regs(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ref_neigh_tab(self%nrefs,self%nrefs))
        self%regs       = 0.d0
        self%regs_grad  = 0.d0
        self%regs_denom = 0.d0
        self%pftcc      => pftcc
        if( params_glob%l_reg_neigh )then
            call self%partition_refs
            call self%make_neigh_tab
        else
            self%ref_neigh_tab = .true.
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
        type(ori) :: o_prev
        !$omp parallel do collapse(2) default(shared) private(i,iref,iptcl,inpl_corrs, o_prev) proc_bind(close) schedule(static)
        do iref = 1, self%nrefs
            do i = 1, self%pftcc%nptcls
                iptcl = glob_pinds(i)
                call build_glob%spproj_field%get_ori(iptcl, o_prev) ! previous ori
                ! find best irot/shift for this pair of iref, iptcl
                call self%pftcc%gencorrs( iref, iptcl, -o_prev%get_2Dshift(), inpl_corrs )
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
        complex(sp),        pointer       :: shmat(:,:)
        integer     :: iptcl, iref, ithr, loc, pind_here, iref2
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        type(ori)   :: o_prev
        real        :: grad_w
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iptcl = params_glob%fromp, params_glob%top
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                iref = best_ir(iptcl)
                ithr = omp_get_thread_num() + 1
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = self%ref_ptcl_tab(iref, iptcl)%loc
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                shmat => self%pftcc%heap_vars(ithr)%shmat
                call build_glob%spproj_field%get_ori(iptcl, o_prev)           ! previous ori
                call self%pftcc%gen_shmat(ithr, o_prev%get_2Dshift(), shmat)
                call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp), ptcl_ctf_rot, loc)
                call self%rotate_polar(self%pftcc%ctfmats(:,:,pind_here),                    ctf_rot, loc)
                grad_w = 1./self%ref_ptcl_tab(iref, iptcl)%sum - self%ref_ptcl_tab(iref, iptcl)%w/self%ref_ptcl_tab(iref, iptcl)%sum**2
                self%regs(:,:,iref)       = self%regs(:,:,iref)       + ptcl_ctf_rot
                self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + ctf_rot**2
                self%regs_grad(:,:,iref)  = self%regs_grad(:,:,iref)  + ptcl_ctf_rot * grad_w
                do iref2 = 1, self%nrefs
                    if( iref2 /= iref )then
                        loc = self%ref_ptcl_tab(iref2, iptcl)%loc
                        loc = (self%nrots+1)-(loc-1)
                        if( loc > self%nrots ) loc = loc - self%nrots
                        call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp), ptcl_ctf_rot, loc)
                        self%regs_grad(:,:,iref)  = self%regs_grad(:,:,iref) - ptcl_ctf_rot * self%ref_ptcl_tab(iref2, iptcl)%w/self%ref_ptcl_tab(iref2, iptcl)%sum**2
                    endif
                enddo
            endif
        enddo
    end subroutine form_cavgs

    subroutine prev_cavgs( self )
        class(regularizer), intent(inout) :: self
        complex(sp),        pointer       :: shmat(:,:)
        type(ori)   :: o_prev
        integer     :: iptcl, iref, ithr, loc, pind_here
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iptcl = params_glob%fromp, params_glob%top
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                call build_glob%spproj_field%get_ori(iptcl, o_prev)     ! previous ori
                iref  = build_glob%eulspace%find_closest_proj(o_prev)   ! previous projection direction
                ithr  = omp_get_thread_num() + 1
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                loc = self%pftcc%get_roind(360.-o_prev%e3get())
                loc = (self%nrots+1)-(loc-1)
                if( loc > self%nrots ) loc = loc - self%nrots
                shmat => self%pftcc%heap_vars(ithr)%shmat
                call self%pftcc%gen_shmat(ithr, o_prev%get_2Dshift(), shmat)
                call self%rotate_polar(cmplx(ptcl_ctf(:,:,pind_here), kind=dp), ptcl_ctf_rot, loc)
                ptcl_ctf_rot = ptcl_ctf_rot * shmat
                call self%rotate_polar(self%pftcc%ctfmats(:,:,pind_here),            ctf_rot, loc)
                self%regs(:,:,iref)       = self%regs(:,:,iref)       + ptcl_ctf_rot
                self%regs_denom(:,:,iref) = self%regs_denom(:,:,iref) + ctf_rot**2
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
        integer :: iref, iptcl, jref
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(iptcl,iref,jref,sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            do iref = 1, self%nrefs
                sum_corr = 0.
                do jref = 1, self%nrefs
                    if( self%ref_neigh_tab(iref, jref) )then
                        sum_corr = sum_corr + self%ref_ptcl_corr(iptcl, jref)
                    endif
                enddo
                if( sum_corr < TINY )then
                    self%ref_ptcl_tab(iref,iptcl)%sum = 1.
                    self%ref_ptcl_tab(iref,iptcl)%w   = 0.
                    self%ref_ptcl_corr(iptcl,iref)    = 0.
                else
                    self%ref_ptcl_tab(iref,iptcl)%sum = sum_corr
                    self%ref_ptcl_tab(iref,iptcl)%w   = self%ref_ptcl_corr(iptcl,iref)
                    self%ref_ptcl_corr(iptcl,iref)    = self%ref_ptcl_corr(iptcl,iref) / sum_corr
                endif
            enddo
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

    subroutine uniform_cluster_sort( self, ncols, cur_ir )
        class(regularizer), intent(inout) :: self
        integer,            intent(in)    :: ncols
        integer,            intent(inout) :: cur_ir(params_glob%fromp:params_glob%top)
        integer :: ir, ip, max_ind_ir, max_ind_ip, max_ip(ncols)
        real    :: max_ir(ncols)
        logical :: mask_ir(ncols), mask_ip(params_glob%fromp:params_glob%top)
        mask_ir = .false.
        mask_ip = .true.
        do
            if( .not.(any(mask_ip)) ) return
            if( .not.(any(mask_ir)) ) mask_ir = .true.
            max_ir = -1.
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
        enddo
    end subroutine uniform_cluster_sort

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
        class(regularizer), intent(inout) :: self
        complex,            allocatable   :: cmat(:,:)
        type(image) :: calc_cavg
        integer :: iref, k, box
        real    :: eps
        !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            where( abs(self%regs_denom(:,k,:)) < TINY )
                self%regs(:,k,:) = 0._dp
            elsewhere
                self%regs(:,k,:) = self%regs(:,k,:) / self%regs_denom(:,k,:)
            endwhere
        enddo
        !$omp end parallel do
        ! output images for debugging
        if( params_glob%l_reg_debug )then
            do iref = 1, self%nrefs
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
        if( params_glob%l_reg_anneal )then
            eps = real(params_glob%which_iter) / real(params_glob%reg_iters)
            eps = min(1., eps)
            if( params_glob%l_reg_grad ) self%regs = eps * self%regs + (1. - eps) * self%regs_grad
            !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                self%pftcc%pfts_refs_even(:,:,iref) = eps * self%pftcc%pfts_refs_even(:,:,iref) + (1. - eps) * self%regs(:,:,iref)
                self%pftcc%pfts_refs_odd( :,:,iref) = eps * self%pftcc%pfts_refs_odd( :,:,iref) + (1. - eps) * self%regs(:,:,iref)
            enddo
            !$omp end parallel do
        else
            if( params_glob%l_reg_grad ) self%regs = self%regs + self%regs_grad
            !$omp parallel do default(shared) private(iref) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                self%pftcc%pfts_refs_even(:,:,iref) = self%regs(:,:,iref)
                self%pftcc%pfts_refs_odd( :,:,iref) = self%regs(:,:,iref)
            enddo
            !$omp end parallel do
        endif
        call self%pftcc%memoize_refs
        call calc_cavg%kill
    end subroutine regularize_refs
    
    subroutine reset_regs( self )
        class(regularizer), intent(inout) :: self
        self%regs       = 0._dp
        self%regs_grad  = 0._dp
        self%regs_denom = 0._dp
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
        deallocate(self%regs,self%regs_denom,self%regs_grad,self%ref_neigh_tab)
        if(allocated(self%ref_neigh_map)) deallocate(self%ref_neigh_map)
        if(allocated(self%ref_ptcl_corr)) deallocate(self%ref_ptcl_corr,self%ref_ptcl_tab,self%ptcl_ref_map)
    end subroutine kill
end module simple_regularizer
