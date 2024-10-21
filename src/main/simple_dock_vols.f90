module simple_dock_vols
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,           only: image
use simple_projector,       only: projector
use simple_projector_hlev,  only: rotvol
use simple_volpft_corrcalc, only: volpft_corrcalc
implicit none

public :: dock_vols
private
#include "simple_local_flags.inc"

type dock_vols
    private
    integer               :: box = 0, box_clip = 0, ldim(3) = [0,0,0], ldim_clip(3) = [0,0,0]
    real                  :: smpd, hp, lp, mskdiam, smpd_clip, scale, eul(3), cc
    type(projector)       :: vol_ref, vol
    type(oris)            :: eulspace, eulspace_sub
    type(volpft_corrcalc) :: vpcc
    type(sym)             :: pgrpsym
    logical, allocatable  :: lnns(:)
contains
    procedure             :: new
    procedure, private    :: set_ref
    procedure, private    :: set_target
    procedure, private    :: setup_srch_spaces
    procedure             :: srch_rots
    procedure             :: rotate_target
end type dock_vols

character(len=*), parameter :: PGRP           = 'c2'
integer,          parameter :: NSPACE         = 20000
integer,          parameter :: NSPACE_SUB     = 500
integer,          parameter :: NBEST          = 3
real,             parameter :: LP2SMPD_TARGET = 1./3.

contains

    subroutine new( self, vol_ref_fname, vol_fname, smpd, hp, lp, mskdiam )
        class(dock_vols), intent(inout) :: self
        character(len=*), intent(in)    :: vol_ref_fname, vol_fname
        real,             intent(in)    :: smpd, hp, lp, mskdiam
        call self%set_ref(vol_ref_fname, smpd, hp, lp, mskdiam)
        call self%set_target(vol_fname)
        call self%setup_srch_spaces
    end subroutine new

    ! (1)
    subroutine set_ref( self, vol_ref_fname, smpd, hp, lp, mskdiam )
        class(dock_vols), intent(inout) :: self
        character(len=*), intent(in)    :: vol_ref_fname
        real,             intent(in)    :: smpd, hp, lp, mskdiam
        real    :: smpd_target, smpd_here
        integer :: ifoo
        self%smpd    = smpd
        self%hp      = hp
        self%lp      = lp
        self%mskdiam = mskdiam
        call find_ldim_nptcls(vol_ref_fname, self%ldim, ifoo, smpd=smpd_here)
        ! HE, I would not trust the smpd from the header
        if( self%ldim(3) /= self%ldim(1) ) THROW_HARD('Only for volumes')
        self%box = self%ldim(1)
        call self%vol_ref%new(self%ldim, self%smpd)
        call self%vol_ref%read(vol_ref_fname)
        call self%vol_ref%mask(self%mskdiam, 'soft')
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
        character(len=*), intent(in)    :: vol_fname
        real    :: smpd_here
        integer :: ldim_here(3), ifoo
        call find_ldim_nptcls(vol_fname, ldim_here, ifoo, smpd=smpd_here)
        if( any(ldim_here /= self%ldim ) ) THROW_HARD('Nonconforming volume dimensions')
        call self%vol%new(self%ldim, self%smpd)
        call self%vol%read(vol_fname)
        call self%vol%mask(self%mskdiam, 'soft')
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

    subroutine srch_rots( self )
        class(dock_vols), intent(inout) :: self
        type(oris)        :: eulspace_refine
        type(ori)         :: e
        real              :: rmats(36*NSPACE_SUB,3,3), eul(3), e3, e3_new, ccs(36*NSPACE_SUB), cc
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
                eul = self%eulspace_sub%get_euler(iproj)
                rmats(cnt,:,:) = euler2m(eul)
                call self%eulspace_sub%e3set(iproj, e3)
            end do
        end do
        ! search
        ccs(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 36*NSPACE_SUB
            ccs(i) = self%vpcc%corr(rmats(i,:,:))
        end do
        !$omp end parallel do
        nloc = maxnloc(ccs, NBEST)
        ! do i = 1, NBEST
        !     eul = m2euler(rmats(nloc(i),:,:))
        !     cc  = ccs(nloc(i))
        !     print *, eul, cc
        ! end do
        ! construct multi-neighborhood search space from subspace peaks
        self%lnns = .false.
        call e%new_ori(is_ptcl=.false.)
        do i = 1, NBEST 
            call e%set_euler(m2euler(rmats(nloc(i),:,:)))
            call self%pgrpsym%nearest_proj_neighbors(self%eulspace, e, 10., self%lnns)
        end do
        nspace_refine = count(self%lnns)
        call self%eulspace%extract_subspace(self%lnns, eulspace_refine)
        ! fill-in the rotation matrices
        allocate(rmats_refine(72*nspace_refine,3,3), ccs_refine(72*nspace_refine), source=0.)
        cnt = 0
        do iproj = 1, nspace_refine 
            do inpl = 1, 72
                cnt    = cnt + 1
                e3     = eulspace_refine%e3get(iproj)
                e3_new = real(inpl-1) * 5.
                call eulspace_refine%e3set(iproj, e3_new)
                eul = eulspace_refine%get_euler(iproj)
                rmats_refine(cnt,:,:) = euler2m(eul)
                call eulspace_refine%e3set(iproj, e3)
            end do
        end do
        ! search
        ccs_refine(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 72*nspace_refine
            ccs_refine(i) = self%vpcc%corr(rmats_refine(i,:,:))
        end do
        !$omp end parallel do
        nloc = maxnloc(ccs_refine, NBEST)
        ! do i = 1, NBEST
        !     eul = m2euler(rmats_refine(nloc(i),:,:))
        !     cc  = ccs_refine(nloc(i))
        !     print *, eul, cc
        ! end do
        ! final refinement step
        ! top orientation from previous step
        eul = m2euler(rmats_refine(nloc(1),:,:))
        call e%set_euler(eul)
        ! fill-in the rotation matrices
        deallocate(rmats_refine, ccs_refine)
        allocate(rmats_refine(360,3,3), ccs_refine(360), source=0.)
        do inpl = 1, 360
            e3     = e%e3get()
            e3_new = real(inpl-1)
            call e%e3set(e3_new)
            rmats_refine(inpl,:,:) = euler2m(e%get_euler())
            call e%e3set(e3)
        end do
        ! search
        ccs_refine(:) = -1.
        !$omp parallel do schedule(static) default(shared) private(i) proc_bind(close)
        do i = 1, 360
            ccs_refine(i) = self%vpcc%corr(rmats_refine(i,:,:))
        end do
        !$omp end parallel do
        nloc     = maxnloc(ccs_refine, NBEST)
        self%eul = m2euler(rmats_refine(nloc(1),:,:))
        self%cc  = ccs_refine(nloc(1))
        ! do i = 1, NBEST
        !     eul = m2euler(rmats_refine(nloc(i),:,:))
        !     cc  = ccs_refine(nloc(i))
        !     print *, eul, cc
        ! end do
        ! destruct
        call e%kill
        call eulspace_refine%kill
    end subroutine srch_rots

    subroutine rotate_target( self, vol_fname, vol_rot_fname )
        class(dock_vols), intent(inout) :: self
        character(len=*), intent(in)    :: vol_fname, vol_rot_fname
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
        vol_rot = rotvol(self%vol, e)
        call vol_rot%write(vol_rot_fname)
        call vol_rot%kill
        call e%kill
    end subroutine rotate_target

end module simple_dock_vols


