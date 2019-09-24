module simple_tseries_preproc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_segmentation, only: otsu_img, otsu_img_robust
use simple_procimgfile,  only: corrfilt_tseries_imgfile
implicit none

public :: init_tseries_preproc, kill_tseries_preproc, tseries_center_and_mask
private
#include "simple_local_flags.inc"

! constants
logical,          parameter   :: DEBUG_HERE         = .true.
integer,          parameter   :: CHUNKSZ            = 500
integer,          parameter   :: STEPSZ             = 1000
! integer,          parameter   :: STEPSZ             = 20
integer,          parameter   :: BIN_GROW           = 2
real,             parameter   :: EXTRA_EDGE         = 6.0
real,             parameter   :: TRS                = 20.0
integer,          parameter   :: COS_EDGE_FALLOFF   = 12
character(len=*), parameter   :: TWIN_AVGS_RAW      = 'time_window_avgs_raw.mrc'
character(len=*), parameter   :: TWIN_AVGS          = 'time_window_avgs.mrc'
character(len=*), parameter   :: TWIN_MASKS         = 'time_window_masks.mrc'
character(len=*), parameter   :: SHIFTED_PTCLS      = 'tseries_shifted.mrc'
character(len=*), parameter   :: BACKGR_SUBTR_PTCLS = 'tseries_backgr_subtr.mrc'
character(len=*), parameter   :: MASKED_PTCLS       = 'tseries_masked.mrc'

! module varables
type(image),      allocatable :: ptcl_imgs(:)         ! all particles in the time-series
type(image),      allocatable :: ptcl_avgs(:)         ! averages over time window
type(image),      allocatable :: ptcl_masks(:)        ! masks for the particle images
type(image),      allocatable :: corr_imgs(:)         ! images for phase corr-based shift search
type(image)                   :: cc_img               ! connected components image
integer,          allocatable :: avg_inds(:)          ! indices that map particles to averages
real,             allocatable :: shifts4avgs(:,:)     ! shifts obtained through center of mass centering
real,             allocatable :: shifts4ptcls(:,:)    ! shifts obtained through phase corr-based shift search
real,             allocatable :: rmat_sum(:,:,:)      ! for OpenMP reduction
integer                       :: ldim(3)              ! logical dimension of 2D image
logical                       :: existence = .false.  ! to flag existence

contains

    subroutine init_tseries_preproc
        integer :: i
        ! first, kill pre-existing
        call kill_tseries_preproc
        allocate(ptcl_imgs(params_glob%nptcls), corr_imgs(nthr_glob), avg_inds(params_glob%nptcls))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        do i=1,nthr_glob
            call corr_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
        end do
        ! allocate real matrices for OpenMP reduction
        ldim = ptcl_imgs(1)%get_ldim()
        allocate(rmat_sum(ldim(1),ldim(2),ldim(3)), shifts4ptcls(params_glob%nptcls,2), source=0.)
        ! prepare thread safe images
        call ptcl_imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        ! flag existence
        existence = .true.
    end subroutine init_tseries_preproc

    ! (1) create time window averages and binarise them using 2D Otsu, identify the largest connected component,
    !     estimate its diameter and turn it into a binary image (mask) for further processing
    ! (2) remove diameter outliers through 1-sigma thresholding
    ! (3) create soft-edged binary masks and apply them to the low-pass filtered time-window averages
    ! (4) use the time window masks to remove the background of the particle images
    ! (5) obtain shifts through mass centering of the masked time window averages
    ! (6) shift the background subtracted images according to the mass center of their corresponding average
    ! (7) shift the raw particle images in the same way
    ! (8) create shifted and spherically masked particle images based on the estimated diameters
    ! OUTPUTS: time_window_avgs.mrc, time_window_masks.mrc, tseries_shifted.mrc, tseries_backgr_subtr.mrc, tseries_masked.mrc
    subroutine tseries_center_and_mask
        integer, allocatable :: ccsizes(:)
        real,    allocatable :: diams(:)
        type(image) :: img_tmp
        real        :: diam_ave, diam_sdev, diam_var, mask_radius
        real        :: one_sigma_thresh, mskrad, mskrad_max
        integer     :: iframe, i, j, fromto(2), cnt, loc(1), ithr
        logical     :: err
        ! allocate diameters, ptcl_avgs & ptcl_masks arrays
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
        end do
        allocate(diams(cnt), ptcl_avgs(cnt), ptcl_masks(cnt), shifts4avgs(cnt,3))
        diams       = 0.
        shifts4avgs = 0.
        ! construct particle averages and masks
        do i=1,cnt
            call ptcl_avgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd)
            call ptcl_masks(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        end do
        write(logfhandle,'(A)') '>>> ESTIMATING PARTICLE DIAMETERS THROUGH BINARY IMAGE PROCESSING'
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            call progress(iframe,params_glob%nptcls)
            cnt = cnt + 1
            ! set indices that map particles to averages
            avg_inds(iframe:) = cnt
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
            ! calculate and write time window average
            call calc_avg
            call ptcl_avgs(cnt)%write(TWIN_AVGS_RAW, cnt)
            ! low-pass filter
            call ptcl_avgs(cnt)%bp(0., params_glob%cenlp)
            ! make a copy
            call img_tmp%copy(ptcl_avgs(cnt))
            ! binarise
            call otsu_img_robust(ptcl_avgs(cnt))
            ! identify connected components
            call ptcl_avgs(cnt)%find_connected_comps(cc_img)
            ! find the largest connected component
            ccsizes = cc_img%size_connected_comps()
            loc     = maxloc(ccsizes)
            ! estimate its diameter
            call cc_img%diameter_cc(loc(1), diams(cnt))
            ! turn it into a binary image for mask creation
            call cc_img%cc2bin(loc(1))
            call ptcl_masks(cnt)%copy(cc_img)
            ! put back the original low-pass filtered average
            call ptcl_avgs(cnt)%copy(img_tmp)
        end do
        call img_tmp%kill
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> REMOVING DIAMETER OUTLIERS WITH 1-SIGMA THRESHOLDING'
        ! 1-sigma-thresh
        call moment(diams, diam_ave, diam_sdev, diam_var, err)
        one_sigma_thresh = diam_ave + diam_sdev
        write(logfhandle,'(A,F6.1)') '>>> % DIAMETERS ABOVE THRESHOLD: ', (real(count(diams > one_sigma_thresh))/real(size(diams)))*100.
        where( diams > one_sigma_thresh ) diams = one_sigma_thresh
        ! some reporting
        mskrad_max = maxval(diams) / 2. + EXTRA_EDGE
        write(logfhandle,'(A,F6.1)') '>>> AVERAGE DIAMETER      (IN PIXELS): ', sum(diams) / real(size(diams))
        write(logfhandle,'(A,F6.1)') '>>> MAX MASK RADIUS (MSK) (IN PIXELS): ', mskrad_max
        write(logfhandle,'(A)') '>>> CREATING SOFT-EDGED BINARY MASKS'
        do i=1,size(ptcl_masks)
            ! mask the binary masks
            mskrad = (diams(i)/2.) + real(BIN_GROW)
            call ptcl_masks(i)%mask(mskrad, 'hard')
            ! grow them a BIN_GROW layers
            call ptcl_masks(i)%grow_bins(BIN_GROW)
            ! apply a cosine edge for softening to avoid Fourier artefacts
            call ptcl_masks(i)%cos_edge(COS_EDGE_FALLOFF)
            ! mask and write the time window averages
            call ptcl_avgs(i)%norm
            call ptcl_avgs(i)%mul(ptcl_masks(i))
            call ptcl_avgs(i)%write(TWIN_AVGS, i)
            call ptcl_masks(i)%write(TWIN_MASKS, i)
        end do
        write(logfhandle,'(A)') '>>> CREATING BACKGROUND SUBTRACTED IMAGES'
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%read(params_glob%stk, i)
            call ptcl_imgs(i)%norm
            call ptcl_imgs(i)%mul(ptcl_masks(avg_inds(i)))
            call ptcl_imgs(i)%write(BACKGR_SUBTR_PTCLS, i)
        end do
        ! turn off the threading of the ptcl_avgs and read the time window avgs
        do i=1,size(ptcl_avgs)
            call ptcl_avgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_avgs(i)%read(TWIN_AVGS, i)
        end do





        write(logfhandle,'(A)') '>>> CENTERING PARTICLE IMAGES'
        !$omp parallel default(shared) private(i,ithr) proc_bind(close)
        !$omp do schedule(static)
        do i=1,size(ptcl_avgs)
            call ptcl_avgs(i)%norm_bin
            call ptcl_avgs(i)%masscen(shifts4avgs(i,:))
            call ptcl_avgs(i)%fft
            call ptcl_avgs(i)%shift2Dserial(shifts4avgs(i,1:2))
            call ptcl_avgs(i)%ifft
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do i=1,params_glob%nptcls
            ! get thread index
            ithr = omp_get_thread_num() + 1
            call corr_imgs(ithr)%copy(ptcl_avgs(avg_inds(i)))
            call corr_imgs(ithr)%fft
            call ptcl_imgs(i)%fft
            call corr_imgs(ithr)%fcorr_shift(ptcl_imgs(i), TRS, shifts4ptcls(i,:))
            call ptcl_imgs(i)%shift2Dserial(shifts4ptcls(i,:))
            call ptcl_imgs(i)%ifft()
        end do
        !$omp end do
        !$omp end parallel

        ! REPLACED THE AVG SHIFT MAPPING WITH PHASECORR BASED SHIFT SEARCH
        ! do i=1,params_glob%nptcls
        !     call ptcl_imgs(i)%fft
        !     call ptcl_imgs(i)%shift2Dserial(shifts4avgs(avg_inds(i),1:2))
        !     call ptcl_imgs(i)%ifft()
        ! end do

        ! read back in the unmasked images
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        ! apply shifts
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%fft()
            ! REPLACED THE AVG SHIFT MAPPING WITH PHASECORR BASED SHIFT SEARCH
            ! call ptcl_imgs(i)%shift2Dserial(shifts4avgs(avg_inds(i),1:2))
            call ptcl_imgs(i)%shift2Dserial(shifts4ptcls(i,:))
            call ptcl_imgs(i)%ifft()
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '>>> WRITING SHIFTED IMAGES TO DISK'
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%write(SHIFTED_PTCLS, i)
        end do
        !$omp parallel do default(shared) private(i,mskrad) proc_bind(close) schedule(static)
        do i=1,params_glob%nptcls
            mskrad = (diams(avg_inds(i))/2.) + EXTRA_EDGE
            call ptcl_imgs(i)%norm
            call ptcl_imgs(i)%mask(mskrad, 'soft')
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '>>> WRITING SPHERICALLY MASKED IMAGES TO DISK'
        do i=1,size(ptcl_avgs)
            call ptcl_avgs(i)%write(TWIN_AVGS, i)
        end do
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%write(MASKED_PTCLS, i)
        end do

        contains

            subroutine calc_avg
                integer :: i, ind
                real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
                rmat_sum = 0.
                !$omp parallel do default(shared) private(i,ind,rmat_ptr) proc_bind(close) schedule(static) reduction(+:rmat_sum)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    call ptcl_imgs(i)%get_rmat_ptr(rmat_ptr)
                    rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3))
                end do
                !$omp end parallel do
                call ptcl_avgs(cnt)%set_rmat(rmat_sum)
                call ptcl_avgs(cnt)%div(real(fromto(2) - fromto(1) + 1))
            end subroutine calc_avg

    end subroutine tseries_center_and_mask

    subroutine kill_tseries_preproc
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            do i=1,size(ptcl_avgs)
                call ptcl_avgs(i)%kill
                call ptcl_masks(i)%kill
            end do
            do i=1,size(corr_imgs)
                call corr_imgs(i)%kill
            end do

            call cc_img%kill
            deallocate(ptcl_imgs,ptcl_avgs,ptcl_masks,avg_inds,shifts4avgs,shifts4ptcls,rmat_sum)
        endif
    end subroutine kill_tseries_preproc

end module simple_tseries_preproc
