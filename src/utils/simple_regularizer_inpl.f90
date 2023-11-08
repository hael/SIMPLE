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
    real,        allocatable :: ref_ptcl_corr(:,:,:)    !< 2D corr table
    integer,     allocatable :: ptcl_ref_map(:)         !< hard-alignment tab
    integer,     allocatable :: ptcl_loc_map(:)         !< hard-alignment tab
    class(polarft_corrcalc), pointer     :: pftcc => null()
    type(reg_params),        allocatable :: ref_ptcl_tab(:,:,:)
    contains
    ! CONSTRUCTOR
    procedure          :: new
    ! PROCEDURES
    procedure          :: init_tab
    procedure          :: fill_tab
    procedure          :: map_ptcl_ref
    procedure          :: reg_uniform_cluster
    procedure          :: uniform_cluster_sort
    procedure          :: form_cavgs
    procedure          :: compute_regs
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
        integer :: irot
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
                &self%regs(self%pftsz,self%kfromto(1):self%kfromto(2),self%nrefs,self%reg_nrots))
        self%regs       = 0.d0
        self%regs_denom = 0.d0
    end subroutine new

    subroutine init_tab( self )
        class(regularizer_inpl), intent(inout) :: self
        integer :: iptcl, iref, irot
        if( .not.(allocated(self%ref_ptcl_corr)) )then
            allocate(self%ref_ptcl_corr(params_glob%fromp:params_glob%top,self%nrefs,self%reg_nrots), source=0.)
            allocate(self%ref_ptcl_tab( self%nrefs,self%reg_nrots,params_glob%fromp:params_glob%top))
            allocate(self%ptcl_ref_map( params_glob%fromp:params_glob%top),self%ptcl_loc_map( params_glob%fromp:params_glob%top))
        endif
        !$omp parallel do collapse(3) default(shared) private(iptcl,irot,iref) proc_bind(close) schedule(static)
        do irot = 1,self%reg_nrots
            do iref = 1,self%nrefs
                do iptcl = params_glob%fromp,params_glob%top
                    self%ref_ptcl_tab(iref,irot,iptcl)%iptcl = iptcl
                    self%ref_ptcl_tab(iref,irot,iptcl)%iref  = iref
                    self%ref_ptcl_tab(iref,irot,iptcl)%loc   = self%rot_inds(irot)
                    self%ref_ptcl_tab(iref,irot,iptcl)%prob  = 0.
                    self%ref_ptcl_tab(iref,irot,iptcl)%sh    = 0.
                    self%ref_ptcl_tab(iref,irot,iptcl)%w     = 0.
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
                        self%ref_ptcl_tab(iref,irot,iptcl)%loc = loc
                        self%ref_ptcl_tab(iref,irot,iptcl)%sh  = 0.
                        if( params_glob%l_reg_grad )then
                            if( self%pftcc%iseven(i) )then
                                call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=dp), loc, ref_rot)
                            else
                                call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_odd( :,:,iref), kind=dp), loc, ref_rot)
                            endif
                            self%ref_ptcl_corr(iptcl,iref,irot) = max(0.,&
                                &real(self%pftcc%gencorr_for_rot_8(iref, iptcl, [0._dp,0._dp], loc,&
                                     &ref_rot + self%regs(:,:,iref,irot))))
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
            do i = 1, self%pftcc%nptcls
                do iref = 1, self%nrefs
                    do irot = 1, self%reg_nrots
                        iptcl = glob_pinds(i)
                        loc   = self%rot_inds(irot)
                        self%ref_ptcl_tab(iref,irot,iptcl)%loc = loc
                        self%ref_ptcl_tab(iref,irot,iptcl)%sh  = 0.
                        self%ref_ptcl_corr(iptcl,iref,irot)    = max(0., real(self%pftcc%gencorr_for_rot_8(iref, iptcl, [0._dp,0._dp], loc)))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine fill_tab

    subroutine form_cavgs( self, best_ir, best_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,            intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer,            intent(in)    :: best_irot(params_glob%fromp:params_glob%top)
        complex(sp),        pointer       :: shmat(:,:)
        integer     :: iptcl, iref, ithr, pind_here, irot
        complex     :: ptcl_ctf(self%pftsz,self%kfromto(1):self%kfromto(2),self%pftcc%nptcls)
        real        :: weight
        complex(dp) :: ptcl_ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf_rot(self%pftsz,self%kfromto(1):self%kfromto(2))
        ptcl_ctf = self%pftcc%pfts_ptcls * self%pftcc%ctfmats
        do iptcl = params_glob%fromp, params_glob%top
            iref  = best_ir(iptcl)
            irot  = best_irot(iptcl)
            if( self%ref_ptcl_tab(iref, irot, iptcl)%prob < TINY ) cycle
            ithr  = omp_get_thread_num() + 1
            if( iptcl >= self%pftcc%pfromto(1) .and. iptcl <= self%pftcc%pfromto(2))then
                pind_here = self%pftcc%pinds(iptcl)
                ! computing the reg terms as the gradients w.r.t 2D references of the probability
                shmat => self%pftcc%heap_vars(ithr)%shmat
                call self%pftcc%gen_shmat(ithr, -real(self%ref_ptcl_tab(iref, irot, iptcl)%sh), shmat)
                ptcl_ctf_rot = cmplx(ptcl_ctf(:,:,pind_here) * shmat, kind=dp)
                ctf_rot      = self%pftcc%ctfmats(:,:,pind_here)
                if( params_glob%l_reg_grad )then
                    weight = 1./self%ref_ptcl_tab(iref, irot, iptcl)%sum - self%ref_ptcl_tab(iref, irot, iptcl)%w/self%ref_ptcl_tab(iref, irot, iptcl)%sum**2
                else
                    weight = self%ref_ptcl_tab(iref, irot, iptcl)%prob
                endif
                self%regs(:,:,iref,irot)       = self%regs(:,:,iref,irot)       + weight * ptcl_ctf_rot
                self%regs_denom(:,:,iref,irot) = self%regs_denom(:,:,iref,irot) + weight * ctf_rot**2
            endif
        enddo
    end subroutine form_cavgs

    subroutine map_ptcl_ref( self, best_ir, best_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: best_ir(params_glob%fromp:params_glob%top)
        integer,                 intent(in)    :: best_irot(params_glob%fromp:params_glob%top)
        integer :: iptcl
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl)
        do iptcl = params_glob%fromp, params_glob%top
            self%ptcl_ref_map(iptcl) = best_ir(iptcl)
            self%ptcl_loc_map(iptcl) = best_irot(iptcl)
        enddo
        !$omp end parallel do
    end subroutine map_ptcl_ref

    subroutine reg_uniform_cluster( self, out_ir, out_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(inout) :: out_ir(params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: out_irot(params_glob%fromp:params_glob%top)
        integer :: iref, iptcl, irot
        real    :: sum_corr
        ! normalize so prob of each ptcl is between [0,1] for all refs
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iptcl, sum_corr)
        do iptcl = params_glob%fromp, params_glob%top
            sum_corr = sum(self%ref_ptcl_corr(iptcl,:,:))
            self%ref_ptcl_tab(:,:,iptcl)%sum = sum_corr
            self%ref_ptcl_tab(:,:,iptcl)%w   = self%ref_ptcl_corr(iptcl,:,:)
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
                    self%ref_ptcl_tab(iref,irot,iptcl)%prob = self%ref_ptcl_corr(iptcl,iref,irot)
                enddo
            enddo
        enddo
        !$omp end parallel do
        ! sorted clustering
        out_ir   = 1
        out_irot = 1
        call self%uniform_cluster_sort(self%nrefs, self%reg_nrots, out_ir, out_irot)
    end subroutine reg_uniform_cluster

    subroutine uniform_cluster_sort( self, ncols, nz, cur_ir, cur_irot )
        class(regularizer_inpl), intent(inout) :: self
        integer,                 intent(in)    :: ncols, nz
        integer,                 intent(inout) :: cur_ir(  params_glob%fromp:params_glob%top)
        integer,                 intent(inout) :: cur_irot(params_glob%fromp:params_glob%top)
        integer :: irot, ir, ip, max_ind_ir(2), max_ind_ip, max_ip(ncols,nz)
        real    :: max_ir(ncols,nz)
        logical :: mask_ir(ncols,nz), mask_ip(params_glob%fromp:params_glob%top)
        mask_ir = .false.
        mask_ip = .true.
        do
            if( .not.(any(mask_ip)) ) return
            if( .not.(any(mask_ir)) ) mask_ir = .true.
            max_ir = -1.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) collapse(2) private(irot,ir,ip)
            do ir = 1, ncols
                do irot = 1, nz
                    if( mask_ir(ir,irot) )then
                        do ip = params_glob%fromp, params_glob%top
                            if( mask_ip(ip) .and. self%ref_ptcl_tab(ir, irot, ip)%prob > max_ir(ir,irot) )then
                                max_ir(ir,irot) = self%ref_ptcl_tab(ir, irot, ip)%prob
                                max_ip(ir,irot) = ip
                            endif
                        enddo
                    endif
                enddo
            enddo
            !$omp end parallel do
            max_ind_ir = maxloc(max_ir, mask=mask_ir)
            max_ind_ip = max_ip(max_ind_ir(1), max_ind_ir(2))
            ! swapping max_ind_ip
            cur_ir(  max_ind_ip) = max_ind_ir(1)
            cur_irot(max_ind_ip) = max_ind_ir(2)
            mask_ip( max_ind_ip) = .false.
            mask_ir( max_ind_ir(1), max_ind_ir(2)) = .false.
        enddo
    end subroutine uniform_cluster_sort

    subroutine compute_regs( self, ref_freq_in )
        use simple_image
        use simple_opt_filter, only: butterworth_filter
        class(regularizer_inpl), intent(inout) :: self
        real,      optional,     intent(in)    :: ref_freq_in
        integer,     allocatable :: ref_ind(:)
        complex,     allocatable :: cmat(:,:)
        complex(dp), allocatable :: regs_tmp(:,:)
        type(image) :: calc_cavg
        integer     :: iref, k, box, irot, cnt, find
        real        :: ref_freq, filt(self%kfromto(1):self%kfromto(2))
        ref_freq = 0.
        if( present(ref_freq_in) ) ref_freq = ref_freq_in
        if( params_glob%l_reg_grad )then
            ! keep regs as the gradient
            !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                self%regs(:,k,:,:) = real(self%regs(:,k,:,:), kind=dp)
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(k) proc_bind(close) schedule(static)
            do k = self%kfromto(1),self%kfromto(2)
                where( abs(self%regs_denom(:,k,:,:)) < TINY )
                    self%regs(:,k,:,:) = 0._dp
                elsewhere
                    self%regs(:,k,:,:) = self%regs(:,k,:,:) / self%regs_denom(:,k,:,:)
                endwhere
            enddo
            !$omp end parallel do
        endif
        ! output images for debugging
        if( params_glob%l_reg_debug )then
            allocate(regs_tmp(self%pftsz,self%kfromto(1):self%kfromto(2)))
            cnt = 1
            do iref = 1, self%nrefs
                do irot = 1, self%reg_nrots
                    call self%pftcc%polar2cartesian(cmplx(self%regs(:,:,iref,irot), kind=sp), cmat, box)
                    call calc_cavg%new([box,box,1], params_glob%smpd * real(params_glob%box)/real(box))
                    call calc_cavg%zero_and_flag_ft
                    call calc_cavg%set_cmat(cmat)
                    call calc_cavg%shift_phorig()
                    call calc_cavg%ifft
                    call calc_cavg%write('inpl_polar_cavg_reg_'//int2str(params_glob%which_iter)//'.mrc', cnt)
                    call self%pftcc%rotate_ref(cmplx(self%pftcc%pfts_refs_even(:,:,iref), kind=dp), self%rot_inds(irot), regs_tmp)
                    call self%pftcc%polar2cartesian(cmplx(regs_tmp, kind=sp), cmat, box)
                    call calc_cavg%zero_and_flag_ft
                    call calc_cavg%set_cmat(cmat)
                    call calc_cavg%shift_phorig()
                    call calc_cavg%ifft
                    call calc_cavg%write('inpl_polar_cavg_'//int2str(params_glob%which_iter)//'.mrc', cnt)
                    cnt = cnt + 1
                enddo
            enddo
        endif
        call calc_cavg%kill
    end subroutine compute_regs
    
    subroutine reset_regs( self )
        class(regularizer_inpl), intent(inout) :: self
        self%regs       = 0._dp
        self%regs_denom = 0._dp
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
        deallocate(self%regs,self%regs_denom,self%rot_inds)
        if(allocated(self%ref_ptcl_corr)) deallocate(self%ref_ptcl_corr,self%ref_ptcl_tab,self%ptcl_ref_map,self%ptcl_loc_map)
    end subroutine kill
end module simple_regularizer_inpl
