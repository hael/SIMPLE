module simple_tseries_preproc
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_segmentation, only: otsu_img, otsu_robust_fast
use simple_procimgfile,  only: corrfilt_tseries_imgfile
implicit none

public :: init_tseries_preproc, kill_tseries_preproc, tseries_estimate_diam
private
#include "simple_local_flags.inc"

! constants
logical,          parameter   :: DEBUG_HERE         = .true.
integer,          parameter   :: CHUNKSZ            = 500
integer,          parameter   :: STEPSZ             = 20
real,             parameter   :: EXTRA_EDGE         = 3.0
character(len=*), parameter   :: TWIN_AVGS_RAW      = 'time_window_avgs_raw.mrc'
character(len=*), parameter   :: TWIN_AVGS          = 'time_window_avgs.mrc'

! module varables
type(image),      allocatable :: ptcl_imgs(:)         ! all particles in the time-series
type(image),      allocatable :: ptcl_avgs(:)         ! averages over time window
type(image)                   :: cc_img               ! connected components image
real,             allocatable :: rmat_sum(:,:,:)      ! for OpenMP reduction
integer                       :: ldim(3)              ! logical dimension of 2D image
logical                       :: existence = .false.  ! to flag existence

contains

    subroutine init_tseries_preproc
        integer :: i
        ! first, kill pre-existing
        call kill_tseries_preproc
        allocate(ptcl_imgs(params_glob%nptcls))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        ! allocate real matrices for OpenMP reduction
        ldim = ptcl_imgs(1)%get_ldim()
        allocate(rmat_sum(ldim(1),ldim(2),ldim(3)), source=0.)
        ! prepare thread safe images
        call ptcl_imgs(1)%construct_thread_safe_tmp_imgs(nthr_glob)
        ! flag existence
        existence = .true.
    end subroutine init_tseries_preproc

    subroutine tseries_estimate_diam( avg_diam, med_diam, max_diam, mskrad )
        real,    intent(out) :: avg_diam, med_diam, max_diam, mskrad
        integer, allocatable :: ccsizes(:)
        real,    allocatable :: diams(:)
        type(image) :: img_tmp
        real        :: diam_ave, diam_sdev, diam_var, mask_radius
        real        :: one_sigma_thresh, thresh(3)
        integer     :: iframe, i, j, fromto(2), cnt, loc(1), ithr
        logical     :: err
        ! allocate diameters & ptcl_avgs arrays
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
        end do
        allocate(diams(cnt), ptcl_avgs(cnt))
        diams = 0.
        ! construct particle averages and masks
        do i=1,cnt
            call ptcl_avgs(i)%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        end do
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
            ! calculate and write time window average
            call calc_avg
            call ptcl_avgs(cnt)%write(TWIN_AVGS_RAW, cnt)
            ! low-pass filter
            call ptcl_avgs(cnt)%bp(0., params_glob%cenlp)
            ! make a copy
            call img_tmp%copy(ptcl_avgs(cnt))
            ! binarise
            call otsu_robust_fast(ptcl_avgs(cnt), is2D=.false., noneg=.false., thresh=thresh)
            ! identify connected components
            call ptcl_avgs(cnt)%find_connected_comps(cc_img)
            ! find the largest connected component
            ccsizes = cc_img%size_connected_comps()
            loc     = maxloc(ccsizes)
            ! estimate its diameter
            call cc_img%diameter_cc(loc(1), diams(cnt))
            ! turn it into a binary image for mask creation
            call cc_img%cc2bin(loc(1))
            ! put back the original low-pass filtered average
            call ptcl_avgs(cnt)%copy(img_tmp)
        end do
        call img_tmp%kill
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> REMOVING DIAMETER OUTLIERS WITH 1-SIGMA THRESHOLDING'
        ! 1-sigma-thresh
        call moment(diams, diam_ave, diam_sdev, diam_var, err)
        one_sigma_thresh = diam_ave + diam_sdev
        write(logfhandle,'(A,F6.1)') '>>> % DIAMETERS ABOVE THRESHOLD: ',&
        &(real(count(diams > one_sigma_thresh))/real(size(diams)))*100.
        where( diams > one_sigma_thresh ) diams = one_sigma_thresh
        ! output
        avg_diam = sum(diams) / real(size(diams))
        med_diam = median(diams)
        max_diam = one_sigma_thresh
        mskrad   = max_diam/2. + EXTRA_EDGE
        write(logfhandle,'(A,F6.1)') '>>> AVERAGE DIAMETER      (IN PIXELS): ', avg_diam
        write(logfhandle,'(A,F6.1)') '>>> MEDIAN  DIAMETER      (IN PIXELS): ', med_diam
        write(logfhandle,'(A,F6.1)') '>>> MAXIMUM DIAMETER      (IN PIXELS): ', max_diam
        write(logfhandle,'(A,F6.1)') '>>> MAX MASK RADIUS (MSK) (IN PIXELS): ', mskrad

        contains

            subroutine calc_avg
                integer :: i, ind
                real(kind=c_float), pointer :: rmat_ptr(:,:,:) => null()
                rmat_sum = 0.
                !$omp parallel do default(shared) private(i,ind,rmat_ptr)&
                !$omp proc_bind(close) schedule(static) reduction(+:rmat_sum)
                do i=fromto(1),fromto(2)
                    ind = i - fromto(1) + 1
                    call ptcl_imgs(i)%get_rmat_ptr(rmat_ptr)
                    rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3))
                end do
                !$omp end parallel do
                call ptcl_avgs(cnt)%set_rmat(rmat_sum)
                call ptcl_avgs(cnt)%div(real(fromto(2) - fromto(1) + 1))
            end subroutine calc_avg

    end subroutine tseries_estimate_diam

    subroutine kill_tseries_preproc
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            do i=1,size(ptcl_avgs)
                call ptcl_avgs(i)%kill
            end do
            call cc_img%kill
            deallocate(ptcl_imgs,ptcl_avgs,rmat_sum)
        endif
    end subroutine kill_tseries_preproc

end module simple_tseries_preproc
