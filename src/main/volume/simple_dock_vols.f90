module simple_dock_vols
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,            only: image
use simple_projector,        only: projector
use simple_simple_volinterp, only: rotvol
use simple_volpft_corrcalc,  only: volpft_corrcalc
implicit none

public :: dock_vols
private
#include "simple_local_flags.inc"

type dock_vols
    private
    integer               :: box = 0, box_clip = 0, ldim(3) = [0,0,0], ldim_clip(3) = [0,0,0]
    real                  :: smpd, hp, lp, msk, smpd_clip, scale
    real                  :: eul(3)=0., shift(3)=0., cc
    type(projector)       :: vol_ref, vol
    type(oris)            :: eulspace, eulspace_sub
    type(volpft_corrcalc) :: vpcc
    type(sym)             :: pgrpsym
    logical, allocatable  :: lnns(:)
    logical               :: mag
contains
    procedure             :: new
    procedure             :: get_dock_info
    procedure             :: set_dock_info
    procedure, private    :: set_ref
    procedure, private    :: set_target
    procedure, private    :: setup_srch_spaces
    procedure, private    :: srch_rots
    procedure, private    :: rotpeak_interp
    procedure, private    :: srch_shift
    procedure             :: srch
    procedure             :: rotate_target
    procedure             :: kill
end type dock_vols

character(len=*), parameter :: PGRP           = 'c1'
integer,          parameter :: NSPACE         = 20000
integer,          parameter :: NSPACE_SUB     = 500
integer,          parameter :: NBEST          = 3
real,             parameter :: LP2SMPD_TARGET = 1./3.

contains

    subroutine new( self, vol_ref_fname, vol_fname, smpd, hp, lp, mskdiam, mag )
        class(dock_vols),  intent(inout) :: self
        class(string),     intent(in)    :: vol_ref_fname, vol_fname
        real,              intent(in)    :: smpd, hp, lp, mskdiam
        logical, optional, intent(in)    :: mag
        self%mag = .true.
        if( present(mag) ) self%mag = mag
        call self%set_ref(vol_ref_fname, smpd, hp, lp, mskdiam)
        call self%set_target(vol_fname)
        call self%setup_srch_spaces
    end subroutine new

    subroutine get_dock_info( self, eul, shift, cc )
        class(dock_vols), intent(in)  :: self
        real,             intent(out) :: eul(3), shift(3), cc
        eul   = self%eul
        shift = self%shift
        cc    = self%cc
    end subroutine get_dock_info

    subroutine set_dock_info( self, eul, shift )
        class(dock_vols), intent(inout) :: self
        real,             intent(in)    :: eul(3), shift(3)
        self%eul   = eul
        self%shift = shift
    end subroutine set_dock_info

    ! (1)
    subroutine set_ref( self, vol_ref_fname, smpd, hp, lp, mskdiam )
        class(dock_vols), intent(inout) :: self
        class(string),    intent(in)    :: vol_ref_fname
        real,             intent(in)    :: smpd, hp, lp, mskdiam
        real    :: smpd_target, smpd_here
        integer :: ifoo
        self%smpd    = smpd
        self%hp      = hp
        self%lp      = lp
        self%msk     = mskdiam/self%smpd/2.
        call find_ldim_nptcls(vol_ref_fname, self%ldim, ifoo, smpd=smpd_here)
        ! HE, I would not trust the smpd from the header
        if( self%ldim(3) /= self%ldim(1) ) THROW_HARD('Only for volumes')
        self%box = self%ldim(1)
        call self%vol_ref%new(self%ldim, self%smpd)
        call self%vol_ref%read(vol_ref_fname)
        call self%vol_ref%mask(self%msk, 'soft')
        smpd_target = max(self%smpd, (self%lp * LP2SMPD_TARGET))
        call autoscale(self%box, self%smpd, smpd_target, self%box_clip, self%smpd_clip, self%scale, minbox=64)
        self%ldim_clip = [self%box_clip,self%box_clip,self%box_clip]
        ! clip
        call self%vol_ref%fft
        call self%vol_ref%clip_inplace(self%ldim_clip)
        call self%vol_ref%ifft
        call self%vol_ref%set_smpd(self%smpd_clip) ! safety
    end subroutine set_ref

    ! (2)
    subroutine set_target( self, vol_fname )
        class(dock_vols), intent(inout) :: self
        class(string),    intent(in)    :: vol_fname
        real    :: smpd_here
        integer :: ldim_here(3), ifoo
        call find_ldim_nptcls(vol_fname, ldim_here, ifoo, smpd=smpd_here)
        if( any(ldim_here /= self%ldim ) ) THROW_HARD('Nonconforming volume dimensions')
        call self%vol%new(self%ldim, self%smpd)
        call self%vol%read(vol_fname)
        call self%vol%mask(self%msk, 'soft')
        ! clip
        call self%vol%fft
        call self%vol%clip_inplace(self%ldim_clip)
        call self%vol%ifft
        call self%vol%set_smpd(self%smpd_clip) ! safety
        ! create volpft_corrcalc object
        call self%vpcc%new(self%vol_ref, self%hp, self%lp, KBALPHA, self%vol)
    end subroutine set_target

    ! (3)
    subroutine setup_srch_spaces( self )
        class(dock_vols), intent(inout) :: self
        call self%pgrpsym%new(PGRP)
        call self%eulspace    %new(NSPACE,     is_ptcl=.false.)
        call self%eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
        call self%pgrpsym%build_refspiral(self%eulspace)
        call self%pgrpsym%build_refspiral(self%eulspace_sub)
        allocate(self%lnns(NSPACE), source=.false.)
    end subroutine setup_srch_spaces

    subroutine srch( self )
        class(dock_vols), intent(inout) :: self
        call self%srch_rots
        if( self%mag ) call self%srch_shift
    end subroutine srch

    subroutine srch_rots( self )
        class(dock_vols), intent(inout) :: self
        type(oris)        :: eulspace_refine
        type(ori)         :: e
        real              :: rmats(3,3,36*NSPACE_SUB), eul(3), e3, e3_new, ccs(36*NSPACE_SUB)
        real, allocatable :: rmats_refine(:,:,:), ccs_refine(:)
        integer           :: inpl, cnt, iproj, i, nloc(NBEST), nspace_refine
        ! fill-in the rotation matrices
        cnt = 0
        do iproj = 1, NSPACE_SUB
            do inpl = 1, 36
                cnt    = cnt + 1
                e3     = self%eulspace_sub%e3get(iproj)
                e3_new = real(inpl-1) * 10.
                call self%eulspace_sub%e3set(iproj, e3_new)
                rmats(:,:,cnt) = self%eulspace_sub%get_mat(iproj)
                call self%eulspace_sub%e3set(iproj, e3)
            end do
        end do
        ! search
        ccs(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 36*NSPACE_SUB
            if( self%mag )then
                ccs(i) = self%vpcc%corr_mag(rmats(:,:,i))
            else
                ccs(i) = self%vpcc%corr(rmats(:,:,i))
            endif
        end do
        !$omp end parallel do
        nloc = maxnloc(ccs, NBEST)
        ! construct multi-neighborhood search space from subspace peaks
        self%lnns = .false.
        call e%new_ori(is_ptcl=.false.)
        do i = 1, NBEST 
            call e%set_euler(m2euler(rmats(:,:,nloc(i))))
            call self%pgrpsym%nearest_proj_neighbors(self%eulspace, e, 10., self%lnns)
        end do
        nspace_refine = count(self%lnns)
        call self%eulspace%extract_subspace(self%lnns, eulspace_refine)
        ! fill-in the rotation matrices
        allocate(rmats_refine(3,3,72*nspace_refine), ccs_refine(72*nspace_refine), source=0.)
        cnt = 0
        do iproj = 1, nspace_refine 
            do inpl = 1, 72
                cnt    = cnt + 1
                e3     = eulspace_refine%e3get(iproj)
                e3_new = real(inpl-1) * 5.
                call eulspace_refine%e3set(iproj, e3_new)
                rmats_refine(:,:,cnt) = eulspace_refine%get_mat(iproj)
                call eulspace_refine%e3set(iproj, e3)
            end do
        end do
        ! search
        ccs_refine(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 72*nspace_refine
            if( self%mag )then
                ccs_refine(i) = self%vpcc%corr_mag(rmats_refine(:,:,i))
            else
                ccs_refine(i) = self%vpcc%corr(rmats_refine(:,:,i))
            endif
        end do
        !$omp end parallel do
        nloc = maxnloc(ccs_refine, NBEST)
        ! final refinement step
        ! top orientation from previous step
        eul = m2euler(rmats_refine(:,:,nloc(1)))
        call e%set_euler(eul)
        ! fill-in the rotation matrices
        deallocate(rmats_refine, ccs_refine)
        allocate(rmats_refine(3,3,360), ccs_refine(360), source=0.)
        do inpl = 1, 360
            e3     = e%e3get()
            e3_new = real(inpl-1)
            call e%e3set(e3_new)
            rmats_refine(:,:,inpl) = e%get_mat()
            call e%e3set(e3)
        end do
        ! search with in-plane angular step of 1 degree
        ccs_refine(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 360
            if( self%mag )then
                ccs_refine(i) = self%vpcc%corr_mag(rmats_refine(:,:,i))
            else
                ccs_refine(i) = self%vpcc%corr(rmats_refine(:,:,i))
            endif
        end do
        !$omp end parallel do
        nloc     = maxnloc(ccs_refine, NBEST)
        self%eul = m2euler(rmats_refine(:,:,nloc(1)))
        self%cc  = ccs_refine(nloc(1))
        ! in-plane peak interpolation using previous angular step of 1 degree
        call self%rotpeak_interp(1., ccs_refine, e3_new, self%cc)
        call e%e3set(e3_new)
        self%eul = e%get_euler()
        ! destruct
        call e%kill
        call eulspace_refine%kill
    end subroutine srch_rots

    subroutine rotpeak_interp( self, angstep, corrs, peak_ang, peak_cc )
        class(dock_vols), intent(in) :: self
        real, intent(in)             :: angstep, corrs(:)
        real, intent(out)            :: peak_ang, peak_cc
        real    :: denom, alpha,beta,gamma
        integer :: nangs, maxpos
        nangs  = size(corrs)
        maxpos = maxloc(corrs,dim=1)
        beta   = corrs(maxpos)
        if( maxpos == 1 )then
            alpha = corrs(nangs)
            gamma = corrs(2)
        else if( maxpos == nangs )then
            alpha = corrs(nangs-1)
            gamma = corrs(1)
        else
            alpha = corrs(maxpos-1)
            gamma = corrs(maxpos+1)
        endif
        peak_ang = 0.
        if( alpha<beta .and. gamma<beta )then
            denom = alpha + gamma - 2.*beta
            if( abs(denom) > TINY ) peak_ang = 0.5 * (alpha-gamma) / denom
        endif
        peak_cc  = max(-1.0,min(1.0,beta - 0.25 * peak_ang * (alpha-gamma)))
        peak_ang = angstep * (real(maxpos-1) + peak_ang)
    end subroutine rotpeak_interp

    subroutine srch_shift( self )
        class(dock_vols), intent(inout) :: self
        type(image) :: volrot
        type(ori)   :: e1, e2
        real        :: offset1(3), offset2(3), cc1, cc2, trs
        ! shift limit
        trs = real(self%ldim_clip(1)) / 6.
        call self%vol_ref%fft
        ! first orientation
        e1  = ori(is_ptcl=.true.)
        call e1%set_euler(self%eul)
        call self%vol%ifft
        volrot = rotvol(self%vol, e1)
        call volrot%fft
        call volrot%fcorr_shift3D( self%vol_ref, trs, offset1, peak_interp=.true.)
        offset1 = -matmul(offset1,e1%get_mat())
        cc1     = self%vpcc%corr(e1%get_mat(), offset1)
        ! second orientation
        e2 = e1
        call e2%transp
        volrot = rotvol(self%vol, e2)
        call volrot%fft
        call volrot%fcorr_shift3D( self%vol_ref, trs, offset2, peak_interp=.true.)
        offset2 = -matmul(offset2,e2%get_mat())
        cc2     = self%vpcc%corr(e2%get_mat(), offset2)
        ! solution
        if( cc1 > cc2 )then
            self%shift = (self%smpd_clip / self%smpd) * offset1
            self%cc    = cc1
        else
            self%eul   = e2%get_euler()
            self%shift = (self%smpd_clip / self%smpd) * offset2
            self%cc    = cc2
        endif
        call volrot%kill
        call e1%kill
        call e2%kill
    end subroutine srch_shift

    subroutine rotate_target( self, vol_fname, vol_rot_fname )
        class(dock_vols), intent(inout) :: self
        class(string),    intent(in)    :: vol_fname, vol_rot_fname
        type(image) :: vol_rot
        type(ori)   :: e
        real        :: smpd_here
        integer     :: ldim_here(3), ifoo
        call find_ldim_nptcls(vol_fname, ldim_here, ifoo, smpd=smpd_here)
        if( any(ldim_here /= self%ldim ) ) THROW_HARD('Nonconforming volume dimensions')
        call self%vol%new(self%ldim, self%smpd)
        call self%vol%read(vol_fname)
        call e%new_ori(is_ptcl=.false.)
        call e%set_euler(self%eul)
        if( self%mag )then
            vol_rot = rotvol(self%vol, e, shvec=self%shift)
        else
            vol_rot = rotvol(self%vol, e)
        endif
        call vol_rot%write(vol_rot_fname)
        call vol_rot%kill
        call e%kill
    end subroutine rotate_target

    subroutine kill( self )
        class(dock_vols), intent(inout) :: self
        call self%vol_ref%kill
        call self%vol%kill
        call self%eulspace%kill
        call self%eulspace_sub%kill
        call self%pgrpsym%kill
        if( allocated(self%lnns) ) deallocate(self%lnns)
    end subroutine kill

end module simple_dock_vols


