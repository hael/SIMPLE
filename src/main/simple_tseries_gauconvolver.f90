module simple_tseries_gauconvolver
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tseries_gauconvolver, tseries_gauconvolve, kill_tseries_gauconvolver
private
#include "simple_local_flags.inc"

real,    parameter :: GAU_CUTOFF = 5.0
logical, parameter :: DEBUG_HERE = .false.

type(image), allocatable :: ptcl_imgs(:)        ! all particles in the time-series
type(image), allocatable :: ptcls_twin_corr(:)  ! individual de-noised particle images
type(image)              :: ptcl_twin           ! particles extracted over time window and put as real-space sections in volume
type(image)              :: ptcl_slices_conv    ! convolved particle images (correlation)
type(image)              :: ptcl_sect           ! section of Gaussian convolved volume
type(image)              :: gau_img             ! Gaussian for convolution
integer                  :: nz = 0              ! size of time window
logical                  :: existence = .false. ! to flag existence

contains

    subroutine init_tseries_gauconvolver
        integer     :: fromto(2), i
        ! first, kill pre-existing
        call kill_tseries_gauconvolver
        ! create image objects
        fromto(1) = 1 - params_glob%nframesgrp/2
        fromto(2) = 1 + params_glob%nframesgrp/2 - 1
        nz        = fromto(2) - fromto(1) + 1
        if( .not. is_even(nz) ) THROW_HARD('Z-dim: '//int2str(nz)//' of time window volume must be even, please change nframesgrp; init_tseries_gauconvolver')
        call ptcl_twin%new(     [params_glob%box,params_glob%box,nz], params_glob%smpd)
        call ptcl_slices_conv%new([params_glob%box,params_glob%box,nz], params_glob%smpd)
        call ptcl_slices_conv%set_ft(.true.)
        call ptcl_sect%new( [params_glob%box,params_glob%box,1],  params_glob%smpd)
        call gau_img%new(       [params_glob%box,params_glob%box,nz], params_glob%smpd)
        call gau_img%gauimg3D(params_glob%sigma, params_glob%sigma, params_glob%sigma, cutoff=GAU_CUTOFF)
        call gau_img%fft
        allocate(ptcl_imgs(params_glob%nptcls), ptcls_twin_corr(nz))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
            call ptcl_imgs(i)%read(params_glob%stk, i)
            if( i <= nz ) call ptcls_twin_corr(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
        end do
        ! flag existence
        existence = .true.
    end subroutine init_tseries_gauconvolver

    subroutine tseries_gauconvolve
        real, allocatable :: weights(:)
        integer :: fromto(2), fromtowavg(2), i, iframe, slice, sec_ind
        real    :: sxx, sumw
        logical :: renorm
        do iframe=1,params_glob%nptcls
            call progress(iframe, params_glob%nptcls)
            ! set time window
            fromto(1) = iframe - params_glob%nframesgrp/2
            fromto(2) = iframe + params_glob%nframesgrp/2 - 1
            if(fromto(1) < 1)then
                ! move upward in Z-direction when fromto(1) negative (lower index)
                sec_ind = nz/2 + fromto(1)
            else if(fromto(2) > params_glob%nptcls )then
                ! move downward in Z-direction when fromto(2) > nptcls (higher index)
                sec_ind = nz/2 + (params_glob%nptcls - fromto(2))
            else
                sec_ind = nz / 2
            endif
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > params_glob%nptcls)
                fromto = fromto - 1
            end do
            if( DEBUG_HERE ) print *, 'time window: ', fromto(1), fromto(2), ' sec_ind: ', sec_ind
            ! fill up the time window volume
            call ptcl_twin%zero_and_unflag_ft
            !$omp parallel do default(shared) private(i,slice) schedule(static) proc_bind(close)
            do i=fromto(1),fromto(2)
                slice = i - fromto(1) + 1
                call ptcl_twin%set_slice(slice, ptcl_imgs(i))
            end do
            !$omp end parallel do
            if( DEBUG_HERE ) call ptcl_twin%write('ptcls_twin'//int2str_pad(iframe,6)//'.mrcs')
            ! convolve
            call ptcl_twin%fft
            call ptcl_twin%phase_corr(gau_img,ptcl_slices_conv,2.*params_glob%smpd)
            if( DEBUG_HERE )then
                do i=fromto(1),fromto(2)
                    slice = i - fromto(1) + 1
                    call ptcl_slices_conv%get_slice(slice, ptcl_sect)
                    call ptcl_sect%write('ptcls_denoised'//int2str_pad(iframe,6)//'.mrcs', slice)
                end do
            endif

            ! calculate correlations to central section of time window
            ! call ptcl_sect%prenorm4real_corr(sxx, l_msk)
            ! !$omp parallel do default(shared) private(i,slice) schedule(static) proc_bind(close)
            ! do i=fromto(1),fromto(2)
            !     slice = i - fromto(1) + 1
            !     call ptcl_slices_conv%get_slice(slice, ptcls_twin_corr(slice))
            !     corrs(slice) = ptcl_sect%real_corr_prenorm(ptcls_twin_corr(slice), sxx, l_msk)
            ! end do
            ! !$omp end parallel do
            ! ! derive weights from correlations
            ! weights = corrs2weights(corrs, CORRW_ZSCORE_CRIT, RANK_INV_CRIT)
            ! ! check weights backward in time
            ! renorm     = .false.
            ! fromtowavg = fromto
            ! do i=fromto(1) + nz/2,fromto(1),-1
            !     slice = i - fromto(1) + 1
            !     if( weights(slice) <= TINY )then
            !         weights(:slice) = 0.
            !         fromtowavg(1) = slice
            !         renorm = .true.
            !         exit
            !     endif
            ! end do
            ! ! check weights forward in time
            ! do i=fromto(1) + nz/2,fromto(2)
            !     slice = i - fromto(1) + 1
            !     if( weights(slice) <= TINY )then
            !         weights(slice:) = 0.
            !         fromtowavg(2) = slice
            !         renorm = .true.
            !         exit
            !     endif
            ! end do
            ! sumw = sum(weights)
            ! if( renorm ) weights = weights / sumw
            ! ! create weighted average
            ! ! if( DEBUG_HERE ) print *, 'range for weighted avg calc: ', fromtowavg(1), fromtowavg(2)
            ! call ptcl_avg%zero_and_unflag_ft
            ! ! would  have to do a proper reductio with pointers here to parallelise
            ! do i=fromtowavg(1),fromtowavg(2)
            !     slice = i - fromto(1) + 1
            !     call ptcl_avg%add(ptcls_twin_corr(slice), weights(slice))
            ! end do
            ! call ptcl_avg%div(sumw)
            ! if( DEBUG_HERE )then
            !     do i=fromto(1),fromto(2)
            !         slice = i - fromto(1) + 1
            !         ! print *, 'slice/corr/weight: ', slice, corrs(slice), weights(slice)
            !     end do
            ! endif

            ! write
            call ptcl_slices_conv%get_slice(sec_ind, ptcl_sect)
            call ptcl_sect%write(params_glob%outstk, iframe)
            ! re-initialize for the next iteration
            call ptcl_slices_conv%zero_and_flag_ft
        end do
    end subroutine tseries_gauconvolve

    subroutine kill_tseries_gauconvolver
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            do i=1,size(ptcls_twin_corr)
                call ptcls_twin_corr(i)%kill
            end do
            deallocate(ptcl_imgs,ptcls_twin_corr)
            call ptcl_twin%kill
            call ptcl_slices_conv%kill
            call ptcl_sect%kill
            call gau_img%kill
        endif
    end subroutine kill_tseries_gauconvolver

end module simple_tseries_gauconvolver
