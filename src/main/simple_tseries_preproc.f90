module simple_tseries_preproc
include 'simple_lib.f08'
use simple_parameters,   only: params_glob
use simple_image,        only: image
use simple_segmentation, only: otsu_img
implicit none

public :: init_tseries_preproc, kill_tseries_preproc, tseries_approx_mask_radius
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG_HERE = .true.
integer, parameter :: CHUNKSZ    = 1000
integer, parameter :: STEPSZ     = 50
integer, parameter :: WINSZ      = 5
real,    parameter :: EXTRA_EDGE = 3.0


type(image), allocatable :: ptcl_imgs(:)        ! all particles in the time-series
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
        ! create image objects & arrays
        allocate(ptcl_imgs(params_glob%nptcls), shifts(params_glob%nptcls,3))
        do i=1,params_glob%nptcls
            call ptcl_imgs(i)%new([params_glob%box,params_glob%box,1],  params_glob%smpd)
            call ptcl_imgs(i)%read(params_glob%stk, i)
        end do
        call ptcl_avg%new([params_glob%box,params_glob%box,1], params_glob%smpd)
        ! allocate real matrix for OpenMP reduction
        ldim = ptcl_avg%get_ldim()
        allocate(rmat_sum(ldim(1),ldim(2),ldim(3)), source=0.)
        ! flag existence
        existence = .true.
    end subroutine init_tseries_preproc

    subroutine tseries_approx_mask_radius
        integer, allocatable :: ccsizes(:)
        real,    allocatable :: diams(:)
        logical, allocatable :: include_mask(:)
        real    :: diam_ave, diam_sdev, diam_var, mask_radius, nndiams(2)
        real    :: one_sigma_thresh, boundary_avg, mskrad
        integer :: iframe, i, j, fromto(2), cnt, loc(1), lb, rb
        logical :: err, l_nonzero(2)
        ! allocate diameters array
        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
        end do
        allocate(diams(cnt), source=0.)
        ! determine diameters
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
            call otsu_img(ptcl_avg)
            call ptcl_avg%find_connected_comps(cc_img)
            ccsizes = cc_img%size_connected_comps()
            loc = maxloc(ccsizes)
            call cc_img%diameter_cc(loc(1), diams(cnt))
            if( DEBUG_HERE ) call ptcl_avg%write('chunk_avgs.mrcs', cnt)
        end do
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

        cnt = 0
        do iframe=1,params_glob%nptcls,STEPSZ
            cnt = cnt + 1
            fromto(1) = iframe
            fromto(2) = min(params_glob%nptcls, iframe + STEPSZ - 1)
            do i=fromto(1),fromto(2)
                mskrad = diams(cnt) / 2. + EXTRA_EDGE

                ! center and mask here

            end do
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
                call ptcl_avg%set_rmat(rmat_sum)
                call ptcl_avg%div(real(fromto(2) - fromto(1) + 1))
            end subroutine calc_avg

    end subroutine tseries_approx_mask_radius

    subroutine kill_tseries_preproc
        integer :: i
        if( existence )then
            do i=1,size(ptcl_imgs)
                call ptcl_imgs(i)%kill
            end do
            deallocate(shifts,ptcl_imgs)
        endif
    end subroutine kill_tseries_preproc

end module simple_tseries_preproc
