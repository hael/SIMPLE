! regularizer (with inpl discretization) of the cluster2D and refine3D
module simple_regularizer_inpl
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,        only: params_glob
use simple_ori,               only: geodesic_frobdev
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_pftcc_shsrch_grad, only: pftcc_shsrch_grad  ! gradient-based in-plane angle and shift search
implicit none

public :: regularizer_inpl, reg_params
private
#include "simple_local_flags.inc"

type reg_params
    integer :: iptcl            !< iptcl index
    integer :: iref             !< iref index
    integer :: loc              !< inpl index
    real    :: prob, sh(2), w   !< probability, shift, and weight
    real    :: sum
end type reg_params

type :: regularizer_inpl
    integer                  :: nrots
    integer                  :: nrefs
    integer                  :: pftsz
    integer                  :: reg_nrots
    integer                  :: kfromto(2)
    integer,     allocatable :: rot_inds(:)
    complex(dp), allocatable :: regs(:,:,:,:)           !< -"-, reg terms
    real(dp),    allocatable :: regs_denom(:,:,:,:)     !< -"-, reg denom
    real,        allocatable :: ref_corr(:)             !< total ref corr sum
    real,        allocatable :: ref_ptcl_corr(:,:,:)    !< 2D corr table
    real,        allocatable :: ptcl_ref_map(:)         !< hard-alignment tab
    real,        allocatable :: ptcl_loc_map(:)         !< hard-alignment tab
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(pftcc_shsrch_grad), allocatable :: grad_shsrch_obj(:)
    type(reg_params),        allocatable :: ref_ptcl_tab(:,:,:), ref_ptcl_ori(:,:,:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: init_tab
    procedure          :: fill_tab
    procedure          :: sort_tab
    procedure          :: sort_tab_ptcl
    procedure          :: map_ptcl_ref
    procedure          :: sort_tab_no_norm
    procedure          :: reg_uniform_sort
    procedure          :: uniform_sort_tab
    procedure          :: cluster_sort_tab
    procedure          :: reg_cluster_sort
    procedure          :: uniform_cavgs
    procedure          :: ref_reg_cc_tab
    procedure          :: regularize_refs
    procedure          :: reset_regs
    procedure, private :: calc_raw_frc, calc_pspec
    procedure, private :: rotate_polar_real, rotate_polar_complex, rotate_polar_test
    generic            :: rotate_polar => rotate_polar_real, rotate_polar_complex, rotate_polar_test
    ! DESTRUCTOR
    procedure          :: kill
end type regularizer_inpl

contains
    ! CONSTRUCTORS

    subroutine new( self, pftcc )
        class(regularizer_inpl), target, intent(inout) :: self
        class(polarft_corrcalc), target, intent(inout) :: pftcc
        integer, parameter :: MAXITS = 60
        real    :: lims(2,2), lims_init(2,2)
        integer :: ithr, irot
        self%pftcc   => pftcc
        self%nrots   = pftcc%nrots
        self%nrefs   = pftcc%nrefs
        self%pftsz   = pftcc%pftsz
        self%kfromto = pftcc%kfromto
        self%reg_nrots = params_glob%reg_nrots
        if( self%reg_nrots > self%pftcc%nrots ) self%reg_nrots = self%pftcc%nrots
        allocate(self%rot_inds(self%reg_nrots))
        do irot = 1, self%reg_nrots
            self%rot_inds(irot) = 1 + (irot-1) * int(pftcc%nrots / self%reg_nrots)
        enddo
        ! allocation
        allocate(self%regs_denom(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs,self%reg_nrots),&
                &self%grad_shsrch_obj(params_glob%nthr),self%ref_corr(self%nrefs),&
                &self%regs(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs,self%reg_nrots))
        self%regs       = 0.d0
        self%regs_denom = 0.d0
        self%ref_corr   = 0.
        lims(:,1)       = -params_glob%reg_minshift
        lims(:,2)       =  params_glob%reg_minshift
        lims_init(:,1)  = -SHC_INPL_TRSHWDTH
        lims_init(:,2)  =  SHC_INPL_TRSHWDTH
        do ithr = 1, params_glob%nthr
            call self%grad_shsrch_obj(ithr)%new(lims, lims_init=lims_init,&
            &shbarrier=params_glob%shbarrier, maxits=MAXITS, opt_angle=params_glob%l_reg_opt_ang)
        enddo
    end subroutine new

    subroutine init_tab( self )
        class(regularizer_inpl), intent(inout) :: self
        integer :: iptcl, iref, irot
        if( .not.(allocated(self%ref_ptcl_corr)) )then
            allocate(self%ref_ptcl_corr(params_glob%fromp:params_glob%top,self%nrefs,self%reg_nrots), source=0.)
            allocate(self%ref_ptcl_tab( params_glob%fromp:params_glob%top,self%nrefs,self%reg_nrots))
            allocate(self%ref_ptcl_ori( params_glob%fromp:params_glob%top,self%nrefs,self%reg_nrots))
            allocate(self%ptcl_ref_map( params_glob%fromp:params_glob%top),self%ptcl_loc_map( params_glob%fromp:params_glob%top))
        endif
        !$omp parallel do collapse(3) default(shared) private(iptcl,irot,iref) proc_bind(close) schedule(static)
        do irot = 1,self%reg_nrots
            do iref = 1,self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%iptcl = iptcl
                    self%ref_ptcl_tab(iptcl,iref,irot)%iref  = iref
                    self%ref_ptcl_tab(iptcl,iref,irot)%loc   = self%rot_inds(irot)
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob  = 0.
                    self%ref_ptcl_tab(iptcl,iref,irot)%sh    = 0.
                    self%ref_ptcl_tab(iptcl,iref,irot)%w     = 0.
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine init_tab

    ! filling prob/corr 2D table
    subroutine fill_tab( self, glob_pinds, use_reg )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: glob_pinds(self%pftcc%nptcls)
        logical,       optional, intent(in)    :: use_reg
        complex(dp) :: ref_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer     :: i, iref, iptcl, irot, loc
        real        :: eps
        eps = min(1., real(params_glob%which_iter) / real(params_glob%reg_iters))
        if( present(use_reg) .and. use_reg )then
            !$omp parallel do collapse(3) default(shared) private(i,irot,iref,iptcl,loc,ref_rot) proc_bind(close) schedule(static)
            do irot = 1, self%reg_nrots
                do iref = 1, self%nrefs
                    do i = 1, self%pftcc%nptcls
                        iptcl = glob_pinds(i)
                        loc   = self%rot_inds(irot)
                        self%ref_ptcl_tab(iptcl,iref,irot)%loc = loc
                        self%ref_ptcl_tab(iptcl,iref,irot)%sh  = 0.
                        if( params_glob%l_reg_grad )then
                            if( self%pftcc%iseven(i) )then
                                call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=dp), loc, ref_rot)
                            else
                                call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_odd( :,:,iref), kind=dp), loc, ref_rot)
                            endif
                            self%ref_ptcl_corr(iptcl,iref,irot) = max(0.,&
                                &real(self%pftcc%gencorr_for_rot_8(iref, iptcl, [0._dp,0._dp], loc,&
                                     &eps*ref_rot + (1. - eps)*self%regs(:,:,iref,irot))))
                        else
                            self%ref_ptcl_corr(iptcl,iref,irot) = max(0., &
                                &real(self%pftcc%gencorr_for_rot_8(iref, iptcl, [0._dp,0._dp], loc,&
                                     &self%regs(:,:,iref,irot))))
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) default(shared) private(i,irot,iref,iptcl,loc) proc_bind(close) schedule(static)
            do irot = 1, self%reg_nrots
                do iref = 1, self%nrefs
                    do i = 1, self%pftcc%nptcls
                        iptcl = glob_pinds(i)
                        loc   = self%rot_inds(irot)
                        self%ref_ptcl_tab(iptcl,iref,irot)%loc = loc
                        self%ref_ptcl_tab(iptcl,iref,irot)%sh  = 0.
                        self%ref_ptcl_corr(iptcl,iref,irot)    = max(0., real(self%pftcc%gencorr_for_rot_8(iref, iptcl, [0._dp,0._dp], loc)))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_tab

    subroutine sort_tab( self )
        class(regularizer_inpl), intent(inout) :: self
        integer :: iref, iptcl, irot
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            sum_corr = sum(self%ref_ptcl_corr(iptcl,:,:))
            if( sum_corr < TINY )then
                self%ref_ptcl_corr(iptcl,:,:) = 0.
            else
                self%ref_ptcl_corr(iptcl,:,:) = self%ref_ptcl_corr(iptcl,:,:) / sum_corr
            endif
        enddo
        !$omp end parallel do
        self%ref_ptcl_corr = self%ref_ptcl_corr / maxval(self%ref_ptcl_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%w    = self%ref_ptcl_corr(iptcl,iref,irot)
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%ref_ptcl_ori = self%ref_ptcl_tab
        ! sorting the normalized prob for each iref, to sample only the best #ptcls/#refs ptcls for each iref
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(irot,iref)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                call reg_hpsort(self%ref_ptcl_tab(:,iref,irot))
            enddo
        enddo
        !$omp end parallel do
    end subroutine sort_tab

    subroutine uniform_cavgs( self, best_ip, best_ir, best_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,            intent(in)    :: best_ip(params_glob%fromp:params_glob%top)
        integer,            intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer,            intent(in)    :: best_irot(params_glob%fromp:params_glob%top)
        complex(sp),        pointer       :: shmat(:,:)
        integer     :: i, iptcl, iref, ithr, pind_here, irot
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        real        :: weight
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do i = params_glob%fromp, params_glob%top
            iref  = best_ir(i)
            irot  = best_irot(i)
            iptcl = best_ip(i)
            if( self%ref_ptcl_tab(iptcl, iref, irot)%prob < TINY ) cycle
            ithr  = omp_get_thread_num() + 1
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                shmat => self%pftcc%heap_vars(ithr)%shmat
                call self%pftcc%gen_shmat(ithr, -real(self%ref_ptcl_tab(iptcl, iref, irot)%sh), shmat)
                ptcl_ctf_rot = cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp)
                ctf_rot      = self%pftcc%ctfmats(:,:,pind_here)
                if( params_glob%l_reg_grad )then
                    weight = 1./self%ref_ptcl_tab(iptcl, iref, irot)%w - 1./self%ref_ptcl_tab(iptcl, iref, irot)%sum
                else
                    weight = self%ref_ptcl_tab(iptcl, iref, irot)%prob
                endif
                self%regs(:,:,iref,irot)       = self%regs(:,:,iref,irot)       + weight * ptcl_ctf_rot
                self%regs_denom(:,:,iref,irot) = self%regs_denom(:,:,iref,irot) + weight * ctf_rot**2
                self%ref_corr(iref)            = self%ref_corr(iref)       + self%ref_ptcl_tab(iptcl, iref, irot)%prob
            endif
        enddo
    end subroutine uniform_cavgs

    subroutine map_ptcl_ref( self, best_ip, best_ir, best_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: best_ip(params_glob%fromp:params_glob%top)
        integer,                 intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer,                 intent(in)    :: best_irot(params_glob%fromp:params_glob%top)
        integer :: i
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i)
        do i = params_glob%fromp, params_glob%top
            self%ptcl_ref_map(best_ip(i)) = best_ir(i)
            self%ptcl_loc_map(best_ip(i)) = best_irot(i)
        enddo
        !$omp end parallel do
    end subroutine map_ptcl_ref

    subroutine uniform_sort_tab( self, out_ip, out_ir, out_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(inout) :: out_ip(params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: out_ir(params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: out_irot(params_glob%fromp:params_glob%top)
        integer, allocatable :: best_ip(:), best_ir(:), best_irot(:)
        logical, allocatable :: mask_ir(:,:)
        integer :: iref, iptcl, np, irot
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            sum_corr = sum(self%ref_ptcl_corr(iptcl,:,:))
            if( sum_corr < TINY )then
                self%ref_ptcl_corr(iptcl,:,:) = 0.
            else
                self%ref_ptcl_corr(iptcl,:,:) = self%ref_ptcl_corr(iptcl,:,:) / sum_corr
            endif
        enddo
        !$omp end parallel do
        self%ref_ptcl_corr = self%ref_ptcl_corr / maxval(self%ref_ptcl_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%ref_ptcl_ori = self%ref_ptcl_tab
        ! sorted clustering
        np       = params_glob%top-params_glob%fromp+1
        best_ip  = (/(iptcl, iptcl=params_glob%fromp,params_glob%top)/)
        allocate(best_ir(np), best_irot(np), mask_ir(self%nrefs,self%reg_nrots))
        mask_ir   = .false.
        best_ir   = 1
        best_irot = 1
        call self%reg_uniform_sort(np, self%nrefs, self%reg_nrots, best_ip, best_ir, best_irot, mask_ir)
        ! rearranging the tab
        out_ip   = best_ip
        out_ir   = best_ir
        out_irot = best_irot
    end subroutine uniform_sort_tab

    subroutine cluster_sort_tab( self, out_ip, out_ir, out_irot, cur_tab )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(inout) :: out_ip(params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: out_ir(params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: out_irot(params_glob%fromp:params_glob%top)
        type(reg_params), optional, intent(inout) :: cur_tab(:,:,:)
        integer, allocatable :: best_ip(:), best_ir(:), best_irot(:)
        logical, allocatable :: mask_ir(:,:)
        integer :: iref, iptcl, np, from_ind, to_ind, ind, irot
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            sum_corr = sum(self%ref_ptcl_corr(iptcl,:,:))
            self%ref_ptcl_tab(iptcl,:,:)%sum = sum_corr
            self%ref_ptcl_tab(iptcl,:,:)%w   = self%ref_ptcl_corr(iptcl,:,:)
            if( sum_corr < TINY )then
                self%ref_ptcl_corr(iptcl,:,:) = 0.
            else
                self%ref_ptcl_corr(iptcl,:,:) = self%ref_ptcl_corr(iptcl,:,:) / sum_corr
            endif
        enddo
        !$omp end parallel do
        self%ref_ptcl_corr = self%ref_ptcl_corr / maxval(self%ref_ptcl_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
        if( present(cur_tab) )then
            !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
            do irot = 1, self%reg_nrots
                do iref = 1, self%nrefs
                    do iptcl = params_glob%fromp,params_glob%top
                        if( cur_tab(iptcl,iref,irot)%prob > self%ref_ptcl_tab(iptcl,iref,irot)%prob )then
                            self%ref_ptcl_tab(iptcl,iref,irot)%prob = cur_tab(iptcl,iref,irot)%prob
                        else
                            cur_tab(iptcl,iref,irot)%prob = self%ref_ptcl_tab(iptcl,iref,irot)%prob
                        endif
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
        self%ref_ptcl_ori = self%ref_ptcl_tab
        ! sorted clustering
        np       = params_glob%top-params_glob%fromp+1
        best_ip  = (/(iptcl, iptcl=params_glob%fromp,params_glob%top)/)
        allocate(best_ir(np), best_irot(np), mask_ir(self%nrefs,self%reg_nrots))
        mask_ir  = .false.
        best_ir  = 1
        best_irot = 1
        call self%reg_cluster_sort(np, self%nrefs, self%reg_nrots, best_ip, best_ir, best_irot, mask_ir)
        out_ip   = best_ip
        out_ir   = best_ir
        out_irot = best_irot
    end subroutine cluster_sort_tab

    ! recursively sort the columns of a 2D table, w.r.t the sum of the best nrows/ncols
    ! entries of each column
    ! for example:
    !    original table:
    !                   r1          r2          r3
    !           ------------------------------------
    !             c1  |  2           0           6
    !             c2  |  4           0           7
    !             c3  |  3           1           0
    !             c4  |  3           0           2
    !             c5  |  0           4           7
    !             c6  |  6           4           2
    !   sorted table (smallest to largest):
    !                   r2          r1          r3
    !           ------------------------------------
    !             c1  |  0           2           6
    !             c4  |  0           3           2
    !             c3  |  1           3           0
    !             c6  |  4           6           2
    !             c5  |  4           0           7
    !             c2  |  0           4           7
    !
    ! based on this uniformly sorted table, one can cluster:
    !            (c2, c5) -> r3, (c6, c3) -> r1, (c1, c4) -> r2
    subroutine reg_cluster_sort( self, nrows, ncols, nz, cur_id, cur_ir, cur_irot, mask_ir )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: nrows, ncols, nz
        integer,                 intent(inout) :: cur_id(nrows)
        integer,                 intent(inout) :: cur_ir(nrows)
        integer,                 intent(inout) :: cur_irot(nrows)
        logical,                 intent(inout) :: mask_ir(ncols, nz)
        integer :: ir, tmp_id(nrows), max_ind_ir(2), to_ii, num, irot
        real    :: max_ir(ncols,nz)
        num     = params_glob%reg_num
        to_ii   = nrows
        mask_ir = .true.
        cur_ir  = 1
        do
            if( to_ii < 1 ) return
            if( .not.(any(mask_ir)) ) return
            !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(ir,irot,tmp_id)
            do irot = 1, nz
                do ir = 1, ncols
                    if( mask_ir(ir, irot) )then
                        ! sum of 'num' best ptcls
                        tmp_id(1:to_ii) = cur_id(1:to_ii)
                        call reg_hpsort_ind(tmp_id(1:to_ii), self%ref_ptcl_tab(:,ir,irot))
                        max_ir(ir,irot) = sum(self%ref_ptcl_tab(tmp_id(to_ii-num+1:to_ii),ir,irot)%prob)
                    endif
                enddo
            enddo
            !$omp end parallel do
            max_ind_ir      = maxloc(max_ir, mask=mask_ir)
            tmp_id(1:to_ii) = cur_id(1:to_ii)
            call reg_hpsort_ind(tmp_id(1:to_ii), self%ref_ptcl_tab(:,max_ind_ir(1),max_ind_ir(2)))
            mask_ir(max_ind_ir(1),max_ind_ir(2)) = .false.
            cur_id(1:to_ii)                      = tmp_id(1:to_ii)
            cur_ir(to_ii-num+1:to_ii)            = max_ind_ir(1)
            cur_irot(to_ii-num+1:to_ii)          = max_ind_ir(2)
            to_ii = to_ii - num
        enddo
    end subroutine reg_cluster_sort

    subroutine reg_uniform_sort( self, nrows, ncols, nz, cur_id, cur_ir, cur_irot, mask_irr )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: nrows, ncols, nz
        integer,                 intent(inout) :: cur_id(nrows)
        integer,                 intent(inout) :: cur_ir(nrows)
        integer,                 intent(inout) :: cur_irot(nrows)
        logical,                 intent(inout) :: mask_irr(ncols, nz)
        integer :: irot, ir, ip, tmp_i, max_ind_ir(2), max_ind_ip, to_ii
        real    :: max_ir(ncols,nz), max_ip(ncols,nz)
        to_ii = nrows
        do
            if( to_ii < 1 ) return
            if( .not.(any(mask_irr)) ) mask_irr = .true.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(irot,ir,ip)
            do ir = 1, ncols
                do irot = 1, nz
                    if( mask_irr(ir,irot) )then
                        max_ir(ir,irot) = 0.
                        do ip = 1, to_ii
                            if( self%ref_ptcl_tab(cur_id(ip), ir, irot)%prob > max_ir(ir,irot) )then
                                max_ir(ir,irot) = self%ref_ptcl_tab(cur_id(ip), ir, irot)%prob
                                max_ip(ir,irot) = ip
                            endif
                        enddo
                    endif
                enddo
            enddo
            !$omp end parallel do
            max_ind_ir = maxloc(max_ir, mask=mask_irr)
            max_ind_ip = max_ip(max_ind_ir(1), max_ind_ir(2))
            ! swapping max_ind_ip
            tmp_i                = cur_id(to_ii)
            cur_id(to_ii)        = cur_id(max_ind_ip)
            cur_id(max_ind_ip)   = tmp_i
            cur_ir(to_ii)        = max_ind_ir(1)
            cur_irot(to_ii)      = max_ind_ir(2)
            mask_irr(max_ind_ir(1), max_ind_ir(2)) = .false.
            to_ii = to_ii - 1
        enddo
    end subroutine reg_uniform_sort

    ! sorting each ref column without normalization
    subroutine sort_tab_no_norm( self )
        class(regularizer_inpl), intent(inout) :: self
        integer :: iref, iptcl, irot
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%w    = self%ref_ptcl_corr(iptcl,iref,irot)
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%ref_ptcl_ori = self%ref_ptcl_tab
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(irot,iref)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                call reg_hpsort(self%ref_ptcl_tab(:,iref,irot))
            enddo
        enddo
        !$omp end parallel do
    end subroutine sort_tab_no_norm

    subroutine sort_tab_ptcl( self )
        class(regularizer_inpl), intent(inout) :: self
        integer :: iref, iptcl, irot
        real    :: sum_corr
        ! normalize so prob of each weight is between [0,1] for all ptcls
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref, sum_corr)
        do iref = 1, self%nrefs
            sum_corr = sum(self%ref_ptcl_corr(:,iref,:))
            if( sum_corr < TINY )then
                self%ref_ptcl_corr(:,iref,:) = 0.
            else
                self%ref_ptcl_corr(:,iref,:) = self%ref_ptcl_corr(:,iref,:) / sum_corr
            endif
        enddo
        !$omp end parallel do
        self%ref_ptcl_corr = self%ref_ptcl_corr / maxval(self%ref_ptcl_corr)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(3) private(irot,iref,iptcl)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iptcl,iref,irot)%w    = self%ref_ptcl_corr(iptcl,iref,irot)
                    self%ref_ptcl_tab(iptcl,iref,irot)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine sort_tab_ptcl

    ! sorting rarr, but only keep the sorted indeces
    subroutine reg_hpsort_ind( iarr, rarr )
        integer,          intent(inout) :: iarr(:)
        type(reg_params), intent(in)    :: rarr(:)
        type(reg_params) :: ra
        integer          :: i, ir, j, l, ia, n
        n = size(iarr)
        if( n < 2 ) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ia = iarr(l)
                ra = rarr(ia)
            else
                ia = iarr(ir)
                ra = rarr(ia)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ia
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(iarr(j))%prob < rarr(iarr(j+1))%prob) j = j+1
                endif
                if(ra%prob < rarr(iarr(j))%prob)then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ia
            end do
        end do
    end subroutine reg_hpsort_ind

    ! reg_params heapsort from hpsort_4 (smallest last)
    subroutine reg_hpsort( rarr )
        type(reg_params), intent(inout) :: rarr(:)
        integer          :: i, ir, j, l, n
        type(reg_params) :: ra
        n = size(rarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
            else
                ra       = rarr(ir)
                rarr(ir) = rarr(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j)%prob > rarr(j+1)%prob) j = j+1
                endif
                if(ra%prob > rarr(j)%prob)then
                    rarr(i) = rarr(j)
                    i       = j
                    j       = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
            end do
        end do
    end subroutine reg_hpsort

    ! accumulating reference reg terms for each batch of particles, with cc-based global objfunc
    subroutine ref_reg_cc_tab( self, np )
        class(regularizer_inpl), intent(inout) :: self
        integer,       optional, intent(in)    :: np
        complex(sp),        pointer       :: shmat(:,:)
        integer     :: i, iptcl, iref, ithr, ninds, pind_here, irot
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        real        :: weight
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        if( present(np) )then
            ninds = np
        else
            ninds = size(self%ref_ptcl_corr, 1)
        endif
        !$omp parallel do default(shared) proc_bind(close) schedule(static)&
        !$omp private(irot,iref,ithr,i,iptcl,ptcl_ctf_rot,ctf_rot,shmat,pind_here,weight)
        do irot = 1, self%reg_nrots
            do iref = 1, self%nrefs
                ! taking top sorted corrs/probs
                do i = params_glob%fromp,(params_glob%fromp + params_glob%reg_num)
                    if( self%ref_ptcl_tab(i, iref, irot)%prob < TINY ) cycle
                    ithr  = omp_get_thread_num() + 1
                    iptcl = self%ref_ptcl_tab(i, iref, irot)%iptcl
                    if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                        pind_here = self%pftcc%pinds(iptcl)
                        shmat => self%pftcc%heap_vars(ithr)%shmat
                        call self%pftcc%gen_shmat(ithr, -real(self%ref_ptcl_tab(i, iref, irot)%sh), shmat)
                        ptcl_ctf_rot = cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp)
                        ctf_rot      = self%pftcc%ctfmats(:,:,pind_here)
                        weight       = self%ref_ptcl_tab(i, iref, irot)%prob
                        self%regs(:,:,iref,irot)       = self%regs(:,:,iref,irot)       + weight * ptcl_ctf_rot
                        self%regs_denom(:,:,iref,irot) = self%regs_denom(:,:,iref,irot) + weight * ctf_rot**2
                        self%ref_corr(iref)            = self%ref_corr(iref)            + self%ref_ptcl_tab(i, iref,irot)%prob
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine ref_reg_cc_tab

    subroutine regularize_refs( self, ref_freq_in )
        use simple_image
        use simple_opt_filter, only: butterworth_filter
        class(regularizer_inpl), intent(inout) :: self
        real,      optional,     intent(in)    :: ref_freq_in
        real,        parameter   :: REF_FRAC = 1
        integer,     allocatable :: ref_ind(:)
        complex,     allocatable :: cmat(:,:)
        complex(dp), allocatable :: regs_tmp(:,:)
        type(image) :: calc_cavg
        integer     :: iref, k, box, irot, cnt, find
        real        :: ref_freq, filt(self%kfromto(1):self%kfromto(2))
        ref_freq = 0.
        if( present(ref_freq_in) ) ref_freq = ref_freq_in
        if( params_glob%l_reg_grad )then
            ! keep regs
        else
            !$omp parallel default(shared) private(k) proc_bind(close)
            !$omp do schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                where( abs(self%regs_denom(:,k,:,:)) < TINY )
                    self%regs(:,k,:,:) = 0._dp
                elsewhere
                    self%regs(:,k,:,:) = self%regs(:,k,:,:) / self%regs_denom(:,k,:,:)
                endwhere
            enddo
            !$omp end do
            !$omp end parallel
        endif
        ! applying butterworth filter at cut-off = lp
        find = calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
        call butterworth_filter(find, self%kfromto, filt)
        !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
        do k = self%kfromto(1),self%kfromto(2)
            self%regs(:,k,:,:) = filt(k) * self%regs(:,k,:,:)
        enddo
        !$omp end parallel do
        ! sort ref_corr to only change refs to regs for high-score cavgs
        ref_ind = (/(iref,iref=1,self%nrefs)/)
        ! call hpsort(self%ref_corr, ref_ind)
        ! call reverse(ref_ind)
        ! output images for debugging
        if( params_glob%l_reg_debug )then
            allocate(regs_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)))
            cnt = 1
            do k = 1, int(self%nrefs * REF_FRAC)
                iref = ref_ind(k)
                do irot = 1, self%reg_nrots
                    call self%pftcc%polar2cartesian(cmplx(self%regs(:,:,iref,irot), kind=sp), cmat, box)
                    call calc_cavg%new([box,box,1], params_glob%smpd * real(params_glob%box)/real(box))
                    call calc_cavg%zero_and_flag_ft
                    call calc_cavg%set_cmat(cmat)
                    call calc_cavg%shift_phorig()
                    call calc_cavg%ifft
                    call calc_cavg%write('polar_cavg_reg_'//int2str(params_glob%which_iter)//'.mrc', cnt)
                    call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=dp), self%rot_inds(irot), regs_tmp)
                    call self%pftcc%polar2cartesian(cmplx(regs_tmp, kind=sp), cmat, box)
                    call calc_cavg%zero_and_flag_ft
                    call calc_cavg%set_cmat(cmat)
                    call calc_cavg%shift_phorig()
                    call calc_cavg%ifft
                    call calc_cavg%write('polar_cavg_'//int2str(params_glob%which_iter)//'.mrc', cnt)
                    cnt = cnt + 1
                enddo
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
                self%pftcc%pfts_refs_even(:,:,iref) = self%regs(:,:,iref,1) ! CONSTANT FIX
                self%pftcc%pfts_refs_odd( :,:,iref) = self%regs(:,:,iref,1) ! CONSTANT FIX
            endif
        enddo
        !$omp end do
        !$omp end parallel
        call self%pftcc%memoize_refs
        call calc_cavg%kill
    end subroutine regularize_refs
    
    subroutine reset_regs( self )
        class(regularizer_inpl), intent(inout) :: self
        self%regs       = 0._dp
        self%regs_denom = 0._dp
        self%ref_corr   = 0.
    end subroutine reset_regs

    subroutine rotate_polar_real( self, ptcl_ctf, ptcl_ctf_rot, irot )
        class(regularizer_inpl), intent(inout) :: self
        real(sp),                intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),                intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
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
        class(regularizer_inpl), intent(inout) :: self
        complex(dp),             intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),             intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
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
        class(regularizer_inpl), intent(inout) :: self
        real(dp),                intent(in)    :: ptcl_ctf(    self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),                intent(inout) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,                 intent(in)    :: irot
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
        class(regularizer_inpl), intent(inout) :: self
        complex(sp),             intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(sp),             intent(in)    :: pft2(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,                    intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
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
        class(regularizer_inpl), intent(inout) :: self
        complex(dp),             intent(in)    :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        real,                    intent(out)   :: pspec(self%kfromto(1):self%kfromto(2))
        integer :: k
        do k = self%kfromto(1),self%kfromto(2)
            pspec(k) = real( real(sum(pft(:,k)*conjg(pft(:,k))),dp) / real(self%pftsz,dp) )
        end do
    end subroutine calc_pspec

    ! DESTRUCTOR

    subroutine kill( self )
        class(regularizer_inpl), intent(inout) :: self
        deallocate(self%regs,self%regs_denom,self%grad_shsrch_obj,self%ref_corr,self%rot_inds)
        if(allocated(self%ref_ptcl_corr)) deallocate(self%ref_ptcl_corr,self%ref_ptcl_tab,self%ref_ptcl_ori,self%ptcl_ref_map,self%ptcl_loc_map)
    end subroutine kill
end module simple_regularizer_inpl
