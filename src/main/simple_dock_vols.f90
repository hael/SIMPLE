module simple_dock_vols
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,           only: image
use simple_projector,       only: projector
use simple_projector_hlev,  only: rotvol_slim
use simple_volpft_corrcalc, only: volpft_corrcalc
implicit none

public :: dock_vols
private
#include "simple_local_flags.inc"

type dock_vols
    private
    integer               :: box = 0, box_clip = 0, ldim(3) = [0,0,0], ldim_clip(3) = [0,0,0]
    real                  :: smpd, hp, lp, mskdiam, smpd_clip, scale
    type(projector)       :: vol_ref, vol
    type(oris)            :: eulspace, eulspace_sub
    type(volpft_corrcalc) :: vpcc
contains
    procedure             :: new
    procedure, private    :: set_ref
    procedure, private    :: set_target
    procedure, private    :: setup_srch_spaces
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
        type(sym) :: pgrpsym 
        call pgrpsym%new(PGRP)
        call self%eulspace    %new(NSPACE,     is_ptcl=.false.)
        call self%eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
        call pgrpsym%build_refspiral(self%eulspace)
        call pgrpsym%build_refspiral(self%eulspace_sub)
        call pgrpsym%kill
    end subroutine setup_srch_spaces

end module simple_dock_vols


