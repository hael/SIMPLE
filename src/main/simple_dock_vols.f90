module simple_dock_vols
include 'simple_lib.f08'
!$ use omp_lib
!$ use omp_lib_kinds
use simple_image,          only: image
use simple_projector,      only: projector
use simple_projector_hlev, only: rotvol_slim
implicit none

public :: set_ref
private
#include "simple_local_flags.inc"

character(len=*), parameter :: PGRP           = 'c2'
integer,          parameter :: NSPACE         = 20000
integer,          parameter :: NSPACE_SUB     = 500
integer,          parameter :: NBEST          = 3
real,             parameter :: LP2SMPD_TARGET = 1./3.


integer         :: nvols = 0, box = 0, box_clip = 0, box_pad = 0
integer         :: ldim(3) = [0,0,0], ldim_clip(3) = [0,0,0], ldim_pad(3) = [0,0,0]
real            :: lp, scale, smpd, smpd_clip
type(projector) :: vol_pad
type(image)     :: vol_rot_pad, vol_rot, vol_ref, vol


type(oris)         :: eulspace, eulspace_sub

contains

    ! (1)
    subroutine set_ref( vol_ref_fname, smpd_in, lp_in )
        character(len=*), intent(in) :: vol_ref_fname
        real,             intent(in) :: smpd_in, lp_in
        real    :: smpd_target, smpd_here
        integer :: ifoo
        smpd        = smpd_in
        lp          = lp_in
        call find_ldim_nptcls(vol_ref_fname, ldim, ifoo, smpd=smpd_here)
        ! HE, I would not trust the smpd from the header
        if( ldim(3) /= ldim(1) ) THROW_HARD('Only for volumes')
        box = ldim(1)
        call vol_ref%new(ldim, smpd)
        call vol_ref%read(vol_ref_fname)
        smpd_target = max(smpd, (lp * LP2SMPD_TARGET))
        call autoscale(box, smpd, smpd_target, box_clip, smpd_clip, scale, minbox=64)
        ldim_clip = [box_clip,box_clip,box_clip]
        ! clip
        call vol_ref%fft
        call vol_ref%clip_inplace(ldim_clip)
        call vol_ref%ifft
        call vol_ref%set_smpd(smpd_clip) ! safety
    end subroutine set_ref

    ! (2)
    subroutine set_target( vol_fname )
        character(len=*), intent(in) :: vol_fname
        real    :: smpd_here
        integer :: ldim_here(3), ifoo
        call find_ldim_nptcls(vol_fname, ldim_here, ifoo, smpd=smpd_here)
        if( any(ldim_here /= ldim ) ) THROW_HARD('Nonconforming volume dimensions')
        call vol%new(ldim, smpd)
        call vol%read(vol_fname)
        ! clip
        call vol%fft
        call vol%clip_inplace(ldim_clip)
        call vol%ifft
        call vol%set_smpd(smpd_clip) ! safety
        ! create rotated version
        call vol_rot%new(ldim_clip, smpd_clip)
        ! create real-space padded versions
        box_pad  = 2 * round2even(KBALPHA * real(box/2))
        ldim_pad = [box_pad, box_pad,box_pad]
        call vol_pad%new(ldim_pad, smpd_clip)
        call vol_rot_pad%new(ldim_pad, smpd_clip)
        call vol%pad(vol_pad)
        ! FFT and expand for interpolation
        call vol_pad%fft
        call vol_pad%expand_cmat(KBALPHA)
    end subroutine set_target

    ! (3)
    subroutine setup_srch_spaces
        type(sym) :: pgrpsym 
        call pgrpsym%new(PGRP)
        call eulspace    %new(NSPACE,     is_ptcl=.false.)
        call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
        call pgrpsym%build_refspiral(eulspace)
        call pgrpsym%build_refspiral(eulspace_sub)
        call pgrpsym%kill
    end subroutine setup_srch_spaces

end module simple_dock_vols


