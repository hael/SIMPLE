module simple_tseries_averager
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: init_tseries_averager, tseries_average, kill_tseries_averager
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE = .true.
integer, parameter :: MAXITS     = 5
real,    parameter :: NSIGMAS    = 6.           !< # standard deviations for outliers detection

type(image), allocatable :: imgs(:)             ! all particles/frames in the time-series
type(image)              :: avg                 ! average over time window
type(stats_struct)       :: cstats              ! correlation statistics
type(stats_struct)       :: wstats              ! weight statistics
logical,     allocatable :: corr_mask(:,:,:)    ! logical mask for corr calc
real,        allocatable :: corrs(:)            ! correlations to weighted average over time window
real,        allocatable :: rmat_sum(:,:,:)     ! for OpenMP reduction
integer                  :: ldim(3)             ! logical dimension of 2D image
integer                  :: nz = 0              ! size of time window
integer                  :: fromto(2)           ! time window range, global
logical                  :: existence = .false. ! to flag existence

contains

    subroutine init_tseries_averager( fnames )
        character(len=LONGSTRLEN), optional, intent(in) :: fnames(:)
        integer     :: i, fromto(2), n, nframes, cnt
        type(image) :: img_tmp
        ! first, kill pre-existing
        call kill_tseries_averager
        if( params_glob%top > 1 )then
            if( .not. present(fnames) ) THROW_HARD('fnames dummy argument required when top (last frame index) is present')
            ! create image objects & arrays
            fromto(1) = params_glob%fromp
            fromto(2) = params_glob%top
        else
            ! create image objects & arrays
            fromto(1) = 1 - params_glob%nframesgrp/2
            fromto(2) = 1 + params_glob%nframesgrp/2 - 1
        endif
        nz = fromto(2) - fromto(1) + 1
        if( .not. is_even(nz) ) THROW_HARD('Z-dim: '//int2str(nz)//' of time window volume must be even, please change nframesgrp; init_tseries_averager')
        if( params_glob%top > 1 )then
            allocate(imgs(nz), corrs(nz))
            nframes = size(fnames)
            call find_ldim_nptcls(fnames(1),ldim,n)
            if( n == 1 .and. ldim(3) == 1 )then
                ! all ok
            else
                write(logfhandle,*) 'ldim(3): ', ldim(3)
                write(logfhandle,*) 'nframes: ', n
                THROW_HARD('init_tseries_averager; assumes one frame per file')
            endif
            do i=1,nz
                call imgs(i)%new(ldim, params_glob%smpd)
                call imgs(i)%read(fnames(i), 1)
                call imgs(i)%norm
            end do
            call avg%new(ldim, params_glob%smpd)
            ! make logical mask for real-space corr calc
            call img_tmp%new(ldim, params_glob%smpd)
        else
            allocate(imgs(params_glob%nptcls), corrs(nz))
            do i=1,params_glob%nptcls
                call imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
                call imgs(i)%read(params_glob%stk, i)
            end do
            call avg%new([params_glob%box,params_glob%box,1], params_glob%smpd)
            ! make logical mask for real-space corr calc
            call img_tmp%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        endif
        img_tmp   = 1.0
        call img_tmp%mask(params_glob%msk, 'hard')
        corr_mask = img_tmp%bin2logical()
        call img_tmp%kill
        ! allocate real matrix for OpenMP reduction
        ldim = avg%get_ldim()
        allocate(rmat_sum(ldim(1),ldim(2),ldim(3)), source=0.)
        ! flag existence
        existence = .true.
    end subroutine init_tseries_averager

    subroutine tseries_average
        real,    allocatable :: weights(:)
        logical, allocatable :: outliers(:,:)
        integer :: fromto_loc(2), i, iframe, ind, ref_ind, n_nonzero, ncured, deadhot(2)
        real    :: w, sumw
        604 format(A,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.3)
        if( params_glob%top > 1 )then
            allocate(weights(nz), source=1.0/real(nz))
            fromto_loc = fromto
            call calc_wavg
            call avg%cure_outliers(ncured, NSIGMAS, deadhot, outliers)
            write(logfhandle,'(a,1x,i7)') '>>> # DEAD PIXELS:', deadhot(1)
            write(logfhandle,'(a,1x,i7)') '>>> # HOT  PIXELS:', deadhot(2)
            call avg%write(params_glob%outstk, 1)
        else
            do iframe=1,params_glob%nptcls
                ! set time window
                fromto_loc(1) = iframe - params_glob%nframesgrp/2
                fromto_loc(2) = iframe + params_glob%nframesgrp/2 - 1
                ! shift the window if it's outside the time-series
                do while(fromto_loc(1) < 1)
                    fromto_loc = fromto_loc + 1
                end do
                do while(fromto_loc(2) > params_glob%nptcls)
                    fromto_loc = fromto_loc - 1
                end do
                ! set average to the particle in the current frame to initialize the process
                call avg%copy(imgs(iframe))
                ! de-noise through weighted averaging in time window
                do i=1,MAXITS
                    ! correlate to average
                    call calc_corrs
                    ! calculate weights
                    call calc_weights
                    ! calculate weighted average
                    call calc_wavg
                end do
                n_nonzero = count(weights > TINY)
                call calc_stats(corrs,   cstats)
                call calc_stats(weights, wstats)
                write(logfhandle,'(A,1X,I7)') '>>> FRAME', iframe
                write(logfhandle,604)         '>>> CORR    AVG/SDEV/MIN/MAX:', cstats%avg, cstats%sdev, cstats%minv, cstats%maxv
                write(logfhandle,604)         '>>> WEIGHT  AVG/SDEV/MIN/MAX:', wstats%avg, wstats%sdev, wstats%minv, wstats%maxv
                write(logfhandle,'(A,1X,I5)') '>>> # NONZERO WEIGHTS:       ', n_nonzero
                call avg%write(params_glob%outstk, iframe)
            end do
        endif

        contains

            subroutine calc_corrs
                integer :: i, ind
                real    :: sxx
                call avg%prenorm4real_corr(sxx, corr_mask)
                !$omp parallel do default(shared) private(i,ind) schedule(static) proc_bind(close)
                do i=fromto_loc(1),fromto_loc(2)
                    ind = i - fromto_loc(1) + 1
                    corrs(ind) = avg%real_corr_prenorm(imgs(i), sxx, corr_mask)
                end do
                !$omp end parallel do
            end subroutine calc_corrs

            subroutine calc_weights
                logical :: renorm
                integer :: i, ind
                ! calculate weights
                weights = corrs2weights(corrs, params_glob%wcrit_enum, norm_sigm=.false.)
                ! check weights backward in time
                renorm = .false.
                do i=fromto_loc(1) + nz/2,fromto_loc(1),-1
                    ind = i - fromto_loc(1) + 1
                    if( weights(ind) <= TINY )then
                        weights(:ind) = 0.
                        renorm = .true.
                        exit
                    endif
                end do
                ! check weights forward in time
                do i=fromto_loc(1) + nz/2,fromto_loc(2)
                    ind = i - fromto_loc(1) + 1
                    if( weights(ind) <= TINY )then
                        weights(ind:) = 0.
                        renorm = .true.
                        exit
                    endif
                end do
                sumw = sum(weights)
                if( renorm ) weights = weights / sumw
            end subroutine calc_weights

            subroutine calc_wavg
                integer :: i, ind
                real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
                rmat_sum = 0.
                if( params_glob%top > 1 )then
                    !$omp parallel do default(shared) private(i,ind,rmat_ptr) proc_bind(close) schedule(static) reduction(+:rmat_sum)
                    do i=fromto_loc(1),fromto_loc(2)
                        ind = i - fromto_loc(1) + 1
                        call imgs(ind)%get_rmat_ptr(rmat_ptr)
                        rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3)) * weights(ind)
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) private(i,ind,rmat_ptr) proc_bind(close) schedule(static) reduction(+:rmat_sum)
                    do i=fromto_loc(1),fromto_loc(2)
                        ind = i - fromto_loc(1) + 1
                        call imgs(i)%get_rmat_ptr(rmat_ptr)
                        rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3)) * weights(ind)
                    end do
                    !$omp end parallel do
                endif
                call avg%set_rmat(rmat_sum)
            end subroutine calc_wavg

    end subroutine tseries_average

    subroutine kill_tseries_averager
        integer :: i
        if( existence )then
            do i=1,size(imgs)
                call imgs(i)%kill
            end do
            deallocate(imgs, corr_mask, corrs, rmat_sum)
            call avg%kill
        endif
    end subroutine kill_tseries_averager

end module simple_tseries_averager
