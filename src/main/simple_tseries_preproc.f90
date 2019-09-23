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


! GROW A COUPLE OF LAYERS

logical,          parameter :: DEBUG_HERE       = .true.
integer,          parameter :: CHUNKSZ          = 500
integer,          parameter :: STEPSZ           = 20
integer,          parameter :: WINSZ_MAXFILT    = 5
real,             parameter :: EXTRA_EDGE       = 6.0
integer,          parameter :: COS_EDGE_FALLOFF = 12
character(len=*), parameter :: SHIFTED_PTCLS    = 'shifted_tseries.mrc'
character(len=*), parameter :: MASKED_PATCLS    = 'backgr_subtr_tseries.mrc'
character(len=*), parameter :: TWIN_AVGS        = 'time_window_avgs.mrc'


type(image), allocatable :: ptcl_imgs(:)         ! all particles in the time-series
type(image), allocatable :: ptcl_avgs(:)         ! averages over time window
type(image), allocatable :: ptcl_masks(:)        ! masks for the particle images
type(image)              :: cc_img               ! connected components image
integer,     allocatable :: avg_inds(:)          ! indices that map particles to averages
real,        allocatable :: shifts4avgs(:,:)     ! shifts obtained through center of mass centering
real,        allocatable :: rmat_sum(:,:,:)      ! for OpenMP reduction
integer                  :: ldim(3)              ! logical dimension of 2D image
logical                  :: existence = .false.  ! to flag existence

contains

    subroutine init_tseries_preproc
        integer :: i
        ! first, kill pre-existing
        call kill_tseries_preproc
        allocate(ptcl_imgs(params_glob%nptcls), avg_inds(params_glob%nptcls))
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

    subroutine tseries_center_and_mask
        integer, allocatable :: ccsizes(:)
        real,    allocatable :: diams(:)
        logical, allocatable :: include_mask(:)
        type(image) :: img_tmp
        real        :: diam_ave, diam_sdev, diam_var, mask_radius, nndiams(2)
        real        :: one_sigma_thresh, boundary_avg, mskrad, mskrad_max
        integer     :: iframe, i, j, fromto(2), cnt, loc(1), lb, rb, n_incl
        logical     :: err, l_nonzero(2)
        ! allocate diameters, ptcl_avgs & ptcl_masks arrays
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
        end do
        allocate(diams(cnt), ptcl_avgs(cnt), ptcl_masks(cnt), shifts4avgs(cnt,3))
        diams = 0.
        ! construct particle averages (used later for centering the individual particles)
        do i=1,cnt
            call ptcl_avgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
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
            call calc_avg
            call ptcl_avgs(cnt)%bp(0., params_glob%cenlp)
            call ptcl_avgs(cnt)%mask(real(params_glob%box/2) - EXTRA_EDGE, 'soft')
            call img_tmp%copy(ptcl_avgs(cnt))
            call otsu_img_robust(ptcl_avgs(cnt))
            call ptcl_avgs(cnt)%find_connected_comps(cc_img)
            ccsizes = cc_img%size_connected_comps()
            loc     = maxloc(ccsizes)
            call cc_img%diameter_cc(loc(1), diams(cnt))
            call cc_img%cc2bin(loc(1))
            call cc_img%mask(real(params_glob%box/2) - EXTRA_EDGE, 'hard')
            call cc_img%cos_edge(COS_EDGE_FALLOFF)
            call ptcl_masks(cnt)%copy(cc_img)
            call img_tmp%mul(ptcl_masks(cnt))
            call img_tmp%write(TWIN_AVGS, cnt)
        end do
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> MASKING THE PARICLE IMAGES'
        do i=1,params_glob%nptcls
            call progress(i,params_glob%nptcls)
            call ptcl_imgs(i)%read(params_glob%stk, i)
            call ptcl_imgs(i)%mul(ptcl_masks(avg_inds(i)))
        end do
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> REMOVING DIAMETER OUTLIERS WITH 1-SIGMA THRESHOLDING'
        ! 1-sigma-thresh
        call moment(diams, diam_ave, diam_sdev, diam_var, err)
        one_sigma_thresh = diam_ave + diam_sdev
        include_mask = .not. diams > one_sigma_thresh
        n_incl = count(include_mask)
        write(logfhandle,'(A,F6.1)') '>>> % DIAMETER OUTLIERS: ', (real(n_incl) / real(size(diams))) * 100.
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
        ! max filter the diams array to avoid too tight masking
        do i=1,cnt,WINSZ_MAXFILT
            fromto(1) = i
            fromto(2) = i + WINSZ_MAXFILT - 1
            diams(fromto(1):fromto(2)) = maxval(diams(fromto(1):fromto(2)))
        end do
        ! some reporting
        mskrad_max = maxval(diams) / 2. + EXTRA_EDGE
        write(logfhandle,'(A,F6.1)') '>>> AVERAGE DIAMETER      (IN PIXELS): ', sum(diams) / real(size(diams))
        write(logfhandle,'(A,F6.1)') '>>> MAX MASK RADIUS (MSK) (IN PIXELS): ', mskrad_max
        ! read the time window averages back in and turn off the threading
        do i=1,cnt
            call ptcl_avgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd, wthreads=.false.)
            call ptcl_avgs(i)%read(TWIN_AVGS, i)
        end do
        write(logfhandle,'(A)') '>>> MASS CENTERING TIME WINDOW AVERAGES'
        !$omp parallel default(shared) private(i) proc_bind(close)
        !$omp do schedule(static)
        do i=1,cnt
            call ptcl_avgs(i)%norm_bin
            call ptcl_avgs(i)%masscen(shifts4avgs(i,:))
            call ptcl_avgs(i)%fft
            call ptcl_avgs(i)%shift2Dserial(shifts4avgs(i,1:2))
            call ptcl_avgs(i)%ifft
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%fft
            call ptcl_imgs(i)%shift2Dserial(shifts4avgs(avg_inds(i),1:2))
            call ptcl_imgs(i)%ifft()
        end do
        !$omp end do
        !$omp end parallel
        ! write the background subtracted and shifted images and read back in the unmasked images
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%write(MASKED_PATCLS, i)
            call progress(i,params_glob%nptcls)
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        ! apply shifts
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%fft()
            call ptcl_imgs(i)%shift2Dserial(shifts4avgs(avg_inds(i),1:2))
            call ptcl_imgs(i)%ifft()
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '>>> WRITING SHIFTED IMAGES TO DISK'
        do i=1,size(ptcl_avgs)
            call ptcl_avgs(i)%write(SHIFTED_PTCLS, i)
        end do
        !$omp parallel do default(shared) private(i,mskrad) proc_bind(close) schedule(static)
        do i=1,params_glob%nptcls
            mskrad = diams(avg_inds(i)) / 2. + EXTRA_EDGE
            call ptcl_imgs(i)%norm
            call ptcl_imgs(i)%mask(mskrad, 'soft')
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '>>> WRITING SPHERICALLY MASKED IMAGES TO DISK'
        do i=1,size(ptcl_avgs)
            call ptcl_avgs(i)%write(TWIN_AVGS, i)
        end do
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
                    call ptcl_imgs(i)%get_rmat_ptr(rmat_ptr)
                    rmat_sum = rmat_sum + rmat_ptr(:ldim(1),:ldim(2),:ldim(3))
                end do
                !$omp end parallel do
                call ptcl_avgs(cnt)%set_rmat(rmat_sum)
                call ptcl_avgs(cnt)%div(real(fromto(2) - fromto(1) + 1))
            end subroutine calc_avg

    end subroutine tseries_center_and_mask

    subroutine kill_tseries_preproc
        integer :: i, cnt
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            cnt = 0
            do i=1,params_glob%nptcls,STEPSZ
                cnt = cnt + 1
                call ptcl_avgs(cnt)%kill
            end do
            deallocate(ptcl_imgs,ptcl_avgs,avg_inds,shifts4avgs,rmat_sum)
            call cc_img%kill
        endif
    end subroutine kill_tseries_preproc

end module simple_tseries_preproc
