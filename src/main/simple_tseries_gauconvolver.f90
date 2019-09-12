module simple_tseries_gauconvolver
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tseries_gauconvolver, tseries_gauconvolve, kill_tseries_gauconvolver
private
#include "simple_local_flags.inc"

real, parameter :: GAU_CUTOFF = 5.0

type(image), allocatable :: ptcl_imgs(:)        ! all particles in the time-series
type(image)              :: ptcl_twin           ! particles extracted over time window and put as real-space sections in volume
type(image)              :: ptcl_twin_corr      ! convolved particle images (correlation)
type(image)              :: ptcl_cen_sect       ! central section of Gaussian convolved volume
type(image)              :: gau_img             ! Gaussian for convolution
integer                  :: nz = 0              ! size of time window
logical                  :: existence = .false. ! to flag existence

contains

    subroutine init_tseries_gauconvolver
        integer :: fromto(2), i
        ! first, kill pre-existing
        call kill_tseries_gauconvolver
        ! create image objects
        fromto(1) = 1 - params_glob%nframesgrp/2
        fromto(2) = 1 + params_glob%nframesgrp/2 - 1
        nz        = fromto(2) - fromto(1) + 1
        if( .not. is_even(nz) ) THROW_HARD('Z-dim: '//int2str(nz)//' of time window volume must be even, please change nframesgrp; init_tseries_gauconvolver')
        call ptcl_twin%new(     [params_glob%box,params_glob%box,nz], params_glob%smpd)
        call ptcl_twin_corr%new([params_glob%box,params_glob%box,nz], params_glob%smpd)
        call ptcl_twin_corr%set_ft(.true.)
        call ptcl_cen_sect%new( [params_glob%box,params_glob%box,1],  params_glob%smpd)
        call gau_img%new(       [params_glob%box,params_glob%box,nz], params_glob%smpd)
        call gau_img%gauimg3D(params_glob%sigma, params_glob%sigma, params_glob%sigma, cutoff=GAU_CUTOFF)
        call gau_img%fft
        allocate(ptcl_imgs(params_glob%nptcls))
        do i=1,params_glob%nptcls
             call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
             call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        ! flag existence
        existence = .true.
    end subroutine init_tseries_gauconvolver

    subroutine tseries_gauconvolve
        integer :: fromto(2), slice_last, i, iframe, slice
        do iframe=1,params_glob%nptcls
            call progress(iframe, params_glob%nptcls)
            ! set time window
            fromto(1) = iframe - params_glob%nframesgrp/2
            fromto(2) = iframe + params_glob%nframesgrp/2 - 1
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > params_glob%nptcls)
                fromto = fromto - 1
            end do
            ! fill up the time window volume
            call ptcl_twin%zero_and_unflag_ft
            do i=fromto(1),fromto(2)
                slice = i - fromto(1) + 1
                call ptcl_twin%set_slice(slice, ptcl_imgs(i))
            end do
            ! convolve
            call ptcl_twin%fft
            call ptcl_twin%phase_corr(gau_img,ptcl_twin_corr,params_glob%smpd)
            call ptcl_twin_corr%get_slice(nz/2, ptcl_cen_sect)
            call ptcl_twin_corr%zero_and_flag_ft
            ! write
            call ptcl_cen_sect%write(params_glob%outstk, iframe)
        end do
    end subroutine tseries_gauconvolve

    subroutine kill_tseries_gauconvolver
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            deallocate(ptcl_imgs)
            call ptcl_twin%kill
            call ptcl_twin_corr%kill
            call ptcl_cen_sect%kill
            call gau_img%kill
        endif
    end subroutine kill_tseries_gauconvolver

end module simple_tseries_gauconvolver
