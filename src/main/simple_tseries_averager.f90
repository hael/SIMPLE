module simple_tseries_averager
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tseries_averager, tseries_average, kill_tseries_averager
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE = .false.
integer, parameter :: MAXITS = 5

type(image), allocatable :: ptcl_imgs(:)        ! all particles in the time-series
type(image)              :: ptcl_avg            ! average over time window
logical,     allocatable :: corr_mask(:,:,:)    ! logical mask for corr calc
real,        allocatable :: corrs(:)            ! correlations to weighted average over time window
integer                  :: nz = 0              ! size of time window
logical                  :: existence = .false. ! to flag existence

contains

    subroutine init_tseries_averager
        integer     :: i, fromto(2)
        type(image) :: img_tmp
        ! first, kill pre-existing
        call kill_tseries_averager
        ! create image objects & arrays
        fromto(1) = 1 - params_glob%nframesgrp/2
        fromto(2) = 1 + params_glob%nframesgrp/2 - 1
        nz        = fromto(2) - fromto(1) + 1
        if( .not. is_even(nz) ) THROW_HARD('Z-dim: '//int2str(nz)//' of time window volume must be even, please change nframesgrp; init_tseries_averager')
        allocate(ptcl_imgs(params_glob%nptcls), corrs(nz))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        call ptcl_avg%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        ! make logical mask for real-space corr calc
        call img_tmp%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        img_tmp   = 1.0
        call img_tmp%mask(params_glob%msk, 'hard')
        corr_mask = img_tmp%bin2logical()
        call img_tmp%kill
        ! flag existence
        existence = .true.
        if( DEBUG_HERE ) print *, 'init_tseries_averager, done'
    end subroutine init_tseries_averager

    subroutine tseries_average
        real, allocatable :: weights(:)
        integer :: fromto(2), i, iframe, ind, ref_ind
        real    :: w, sumw
        do iframe=1,params_glob%nptcls
            call progress(iframe, params_glob%nptcls)
            ! set time window
            fromto(1) = iframe - params_glob%nframesgrp/2
            fromto(2) = iframe + params_glob%nframesgrp/2 - 1
            ! shift the window if it's outside the time-series
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > params_glob%nptcls)
                fromto = fromto - 1
            end do
            if( DEBUG_HERE ) print *, 'time window: ', fromto(1), fromto(2)
            ! set average to the particle in the current frame to initialize the process
            call ptcl_avg%copy(ptcl_imgs(iframe))
            ! de-noise through weighted averaging in time window
            do i=1,MAXITS
                ! correlate to average
                call calc_corrs
                ! calculate weights
                call calc_weights
                ! calculate weighted average
                call calc_wavg
            end do
            call ptcl_avg%write(params_glob%outstk, iframe)
        end do


        contains

            subroutine calc_corrs
                integer :: i, ind
                real    :: sxx
                call ptcl_avg%prenorm4real_corr(sxx, corr_mask)
                !$omp parallel do default(shared) private(i,ind) schedule(static) proc_bind(close)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    corrs(ind) = ptcl_avg%real_corr_prenorm(ptcl_imgs(i), sxx, corr_mask)
                end do
                !$omp end parallel do
            end subroutine calc_corrs

            subroutine calc_weights
                logical :: renorm
                integer :: i, ind
                if( DEBUG_HERE )then
                    print *, 'corrw_crit: ',  params_glob%ccw_crit
                    print *, 'rankw_crit: ',  params_glob%rankw_crit
                endif
                ! calculate weights
                if( params_glob%l_rankw )then
                    weights = corrs2weights(corrs, params_glob%ccw_crit, params_glob%rankw_crit, norm_sigm=.false.)
                else
                    weights = corrs2weights(corrs, params_glob%ccw_crit, norm_sigm=.false.)
                endif
                ! check weights backward in time
                renorm = .false.
                do i=fromto(1) + nz/2,fromto(1),-1
                    ind = i - fromto(1) + 1
                    if( weights(ind) <= TINY )then
                        weights(:ind) = 0.
                        renorm = .true.
                        exit
                    endif
                end do
                ! check weights forward in time
                do i=fromto(1) + nz/2,fromto(2)
                    ind = i - fromto(1) + 1
                    if( weights(ind) <= TINY )then
                        weights(ind:) = 0.
                        renorm = .true.
                        exit
                    endif
                end do
                sumw = sum(weights)
                if( renorm ) weights = weights / sumw
                if( DEBUG_HERE )then
                    print *, '*********************'
                    do i=fromto(1),fromto(2)
                        ind = i - fromto(1) + 1
                        print *, 'i/corr/weight: ', i, corrs(ind), weights(ind)
                    end do
                endif
            end subroutine calc_weights

            subroutine calc_wavg
                integer :: i, ind
                call ptcl_avg%zero_and_unflag_ft
                ! HAVE TO DO A PROPER REDUCTION WITH PTRS TO PARALLELIZE THIS ONE
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    call ptcl_avg%add(ptcl_imgs(i), weights(ind))
                end do
                call ptcl_avg%div(sumw)
                ! if( DEBUG_HERE ) call ptcl_avg%write('ptcl_avg.mrc')
            end subroutine calc_wavg

    end subroutine tseries_average

    subroutine kill_tseries_averager
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            deallocate(ptcl_imgs, corr_mask, corrs)
            call ptcl_avg%kill
        endif
    end subroutine kill_tseries_averager

end module simple_tseries_averager
