module simple_tseries_averager
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tseries_averager, tseries_average, kill_tseries_averager
private
#include "simple_local_flags.inc"

type(image), allocatable :: ptcl_imgs(:)        ! particle extracted from individual frame images over time window
type(image)              :: ptcl_avg            ! average over time window
logical,     allocatable :: corr_mask(:,:,:)    ! logical mask for corr calc
real,        allocatable :: corrs(:)            ! correlations to weightied average over time window
logical                  :: existence = .false. ! to flag existence

integer, parameter :: MAXITS=3

contains

    subroutine init_tseries_averager
        integer     :: iframe, n
        type(image) :: img_tmp
        ! first, kill pre-existing
        call kill_tseries_averager
        ! create image objects & arrays
        n = params_glob%nframesgrp + 2
        allocate(ptcl_imgs(n), corrs(n))
        do iframe=1,params_glob%nframesgrp + 2
            call ptcl_imgs(iframe)%new([params_glob%box,params_glob%box,1], params_glob%smpd, wthreads=.false.)
        end do
        call ptcl_avg%new([params_glob%box,params_glob%box,1], params_glob%smpd, wthreads=.false.)
        ! make logical mask for real-space corr calc
        call img_tmp%new([params_glob%box,params_glob%box,1], params_glob%smpd, wthreads=.false.)
        img_tmp   = 1.0
        call img_tmp%mask(params_glob%msk, 'hard')
        corr_mask = img_tmp%bin2logical()
        call img_tmp%kill
        ! flag existence
        existence = .true.
    end subroutine init_tseries_averager

    subroutine tseries_average
        real, allocatable :: weights(:)
        integer :: fromto(2), ind_last, i, iframe, ind
        real    :: w, sumw
        do iframe=1,params_glob%nptcls
            ! set time window
            fromto(1) = iframe - params_glob%nframesgrp/2
            fromto(2) = iframe + params_glob%nframesgrp/2 - 1
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > params_glob%nptcls)
                fromto = fromto - 1
            end do
            ! create first average
            call ptcl_avg%zero_and_unflag_ft
            w        = 1. / real(fromto(2) - fromto(1) + 1)
            sumw     = 0.
            ind_last = fromto(2) - fromto(1) + 1
            do i=fromto(1),fromto(2)
                ind = i - fromto(1) + 1
                call ptcl_imgs(ind)%read(params_glob%stk, i)
                call ptcl_avg%add(ptcl_imgs(ind), w)
                sumw = sumw + w
            end do
            call ptcl_avg%div(sumw)
            ! de-noise through weighted averaging in time window
            do i=1,MAXITS
                ! correlate to avergae
                call calc_corrs
                ! calculate weights
                weights = corrs2weights(corrs(:ind_last), params_glob%ccw_crit, params_glob%rankw_crit)
                sumw    = sum(weights)
                ! calculate weighted average
                call calc_wavg
            end do
            call ptcl_avg%write(params_glob%outstk)
        end do


        contains

            subroutine calc_corrs
                integer :: i, ind
                real    :: sxx
                call ptcl_avg%prenorm4real_corr(sxx, corr_mask)
                !$omp parallel do default(shared) private(i,ind) schedule(static) proc_bind(close)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    corrs(ind) = ptcl_avg%real_corr_prenorm(ptcl_imgs(ind), sxx, corr_mask)
                end do
                !$omp end parallel do
            end subroutine calc_corrs

            subroutine calc_wavg
                integer :: i, ind
                call ptcl_avg%zero_and_unflag_ft
                !$omp parallel do default(shared) private(i,ind) schedule(static) proc_bind(close)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    call ptcl_avg%add(ptcl_imgs(ind), weights(ind))
                end do
                !$omp end parallel do
                call ptcl_avg%div(sumw)
            end subroutine calc_wavg

    end subroutine tseries_average

    subroutine kill_tseries_averager
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            deallocate(ptcl_imgs)
            call ptcl_avg%kill
        endif
    end subroutine kill_tseries_averager

end module simple_tseries_averager
