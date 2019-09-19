module simple_tseries_preproc
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_segmentation, only: otsu_img
use simple_procimgfile,  only: corrfilt_tseries_imgfile
implicit none

public :: init_tseries_preproc, kill_tseries_preproc, tseries_center_and_mask
private
#include "simple_local_flags.inc"

logical,          parameter :: DEBUG_HERE = .true.
integer,          parameter :: CHUNKSZ    = 500
integer,          parameter :: STEPSZ     = 20
real,             parameter :: EXTRA_EDGE = 6.0
character(len=*), parameter :: DENOISED_TSERIES = 'denoised_tseries.mrcs'


type(image), allocatable :: ptcl_imgs(:)        ! all particles in the time-series
type(image), allocatable :: ptcl_imgs_den(:)    ! denoised particles
type(image)              :: ptcl_avg            ! average over time window
type(image)              :: cc_img              ! connected components image
real,        allocatable :: shifts(:,:)         ! shifts obtained throough center of mass centering
real,        allocatable :: rmat_sum(:,:,:)     ! for OpenMP reduction
integer                  :: ldim(3)             ! logical dimension of 2D image
logical                  :: existence = .false. ! to flag existence

contains

    subroutine init_tseries_preproc
        integer     :: i, fromto(2)
        type(image) :: img_tmp
        ! first, kill pre-existing
        call kill_tseries_preproc
        ! denoise
        write(logfhandle,'(A)') '>>> DENOISING TIME-SERIES USING GLOBAL SPACE-TIME GAUSSIAN CONVOLUTION'
        call corrfilt_tseries_imgfile(params_glob%stk, params_glob%sigma, DENOISED_TSERIES, params_glob%smpd, params_glob%lp)
        ! create image objects & arrays
        allocate(ptcl_imgs(params_glob%nptcls), ptcl_imgs_den(params_glob%nptcls), shifts(params_glob%nptcls,3))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_imgs_den(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_imgs(i)%read(params_glob%stk, i)
            call ptcl_imgs_den(i)%read(DENOISED_TSERIES, i)
        end do
        call ptcl_avg%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        ! allocate real matrix for OpenMP reduction
        ldim = ptcl_avg%get_ldim()
        allocate(rmat_sum(ldim(1),ldim(2),ldim(3)), source=0.)
        ! prepare thread safe images
        call ptcl_avg%construct_thread_safe_tmp_imgs(nthr_glob)
        ! flag existence
        existence = .true.
    end subroutine init_tseries_preproc

    subroutine tseries_center_and_mask
        integer, allocatable :: ccsizes(:)
        real,    allocatable :: diams(:)
        logical, allocatable :: include_mask(:)
        real    :: diam_ave, diam_sdev, diam_var, mask_radius, nndiams(2)
        real    :: one_sigma_thresh, boundary_avg, mskrad, mskrad_max
        integer :: iframe, i, j, fromto(2), cnt, loc(1), lb, rb
        logical :: err, l_nonzero(2)
        ! allocate diameters array
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
        end do
        allocate(diams(cnt), source=0.)
        write(logfhandle,'(A)') '>>> ESTIMATING PARTICLE DIAMETERS THROUGH BINARY IMAGE PROCESSING'
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            call progress(iframe,params_glob%nptcls)
            cnt = cnt + 1
            ! set time window
            fromto(1) = iframe - CHUNKSZ/2
            fromto(2) = iframe + CHUNKSZ/2 - 1
            ! shift the window if it's outside the time-series
            do while(fromto(1) < 1)
                fromto = fromto + 1
            end do
            do while(fromto(2) > params_glob%nptcls)
                fromto = fromto - 1
            end do
            call calc_avg
            call ptcl_avg%bp(0., params_glob%cenlp)
            call ptcl_avg%mask(real(params_glob%box/2) - EXTRA_EDGE, 'soft')
            if( DEBUG_HERE ) call ptcl_avg%write('chunk_avgs.mrcs', cnt)
            call otsu_img(ptcl_avg)
            if( DEBUG_HERE ) call ptcl_avg%write('otsu.mrcs', cnt)
            call ptcl_avg%find_connected_comps(cc_img)
            ccsizes = cc_img%size_connected_comps()
            loc = maxloc(ccsizes)
            call cc_img%diameter_cc(loc(1), diams(cnt))
        end do
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> CORRECTING FOR OUTLIERS WITH 1-SIGMA THRESHOLDING'
        ! 1-sigma-thresh
        call moment(diams, diam_ave, diam_sdev, diam_var, err)
        one_sigma_thresh = diam_ave + diam_sdev
        include_mask = .not. diams > one_sigma_thresh
        ! correct for the outliers
        do i=2,cnt
            if( .not. include_mask(i) .and. include_mask(i-1) )then ! boundary
                lb = i-1 ! left bound
                rb = lb  ! right bound
                do j=i,cnt
                    if( .not. include_mask(j) )then
                        cycle
                    else
                        rb = j
                        exit
                    endif
                end do
                boundary_avg = (diams(lb) + diams(rb)) / 2.
                do j=lb+1,rb-1
                    diams(j) = boundary_avg
                end do
            endif
        enddo
        diams(1) = diams(2)
        write(logfhandle,'(A)') '>>> ESTIMATING MASS CENTERS'
        mskrad_max = maxval(diams) /2. + EXTRA_EDGE
        do iframe=1,params_glob%nptcls,STEPSZ
            fromto(1) = iframe
            fromto(2) = min(params_glob%nptcls, iframe + STEPSZ - 1)
            !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
            do i=fromto(1),fromto(2)
                shifts(i,:) = ptcl_imgs(i)%calc_shiftcen_serial(params_glob%cenlp, mskrad_max)
            end do
            !$omp end parallel do
        end do
        write(logfhandle,'(A)') '>>> CENTERING AND MASKING THE ORIGINAL PARTICLE IMAGES'
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
            fromto(1) = iframe
            fromto(2) = min(params_glob%nptcls, iframe + STEPSZ - 1)
            !$omp parallel do default(shared) private(i,mskrad) proc_bind(close) schedule(static)
            do i=fromto(1),fromto(2)
                call ptcl_imgs(i)%fft()
                call ptcl_imgs(i)%shift2Dserial(shifts(i,1:2))
                call ptcl_imgs(i)%ifft()
                mskrad = diams(cnt) / 2. + EXTRA_EDGE
                call ptcl_imgs(i)%mask(mskrad, 'soft')
            end do
            !$omp end parallel do
        end do
        write(logfhandle,'(A)') '>>> WRITING MASKED IMAGES TO DISK'
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%write(params_glob%outstk, i)
        end do

        contains

            subroutine calc_avg
                integer :: i, ind
                real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
                rmat_sum = 0.
                !$omp parallel do default(shared) private(i,ind,rmat_ptr) proc_bind(close) schedule(static) reduction(+:rmat_sum)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    call ptcl_imgs_den(i)%get_rmat_ptr(rmat_ptr)
                    rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3))
                end do
                !$omp end parallel do
                call ptcl_avg%set_rmat(rmat_sum)
                call ptcl_avg%div(real(fromto(2) - fromto(1) + 1))
            end subroutine calc_avg

    end subroutine tseries_center_and_mask

    subroutine kill_tseries_preproc
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
                call ptcl_imgs_den(i)%kill
            end do
            deallocate(ptcl_imgs,ptcl_imgs_den,shifts,rmat_sum)
            call ptcl_avg%kill
            call cc_img%kill
        endif
    end subroutine kill_tseries_preproc

end module simple_tseries_preproc
