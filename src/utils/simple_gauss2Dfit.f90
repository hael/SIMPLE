! Module for fitting the intensity distributions of 2D particle image arrays
! with 2D multivariate Gaussians.  There are two implementations.
! The second implementation includes sorting by correlation
! and an output array of 2D images displaying the Gaussian fits.
module simple_gauss2Dfit

include 'simple_lib.f08'
use simple_image,        only: image
implicit none

public  :: batch_gauss2Dfit, gauss2Dfit
private

interface batch_gauss2Dfit
    module procedure batch_gauss2Dfit_1, batch_gauss2Dfit_2
end interface batch_gauss2Dfit

contains

! NO SORTING OR VISUALIZATION OF FITS.
! Takes in an array of particle images, 'refs.' Fits each particle with a 2D
! multivariate Gaussian distribution.  Outputs same-sized arrays of the fit
! centers, fit covariance matrices, and real-space correlations of fits and refs.
! Centers and covs are reported as pixels or pixels^2.
! The standard deviations of the fits are the square roots of the eigenvalues
! of covs, and the principal axes are the corresponding eigenvectors. 
subroutine batch_gauss2Dfit_1(refs, centers, covs, corrs)
    type(image), allocatable,   intent(inout)   :: refs(:)
    real, allocatable,          intent(out)     :: centers(:,:), covs(:,:,:)
    real, allocatable,          intent(out)     :: corrs(:)
    type(image)         :: fit
    real                :: smpd
    integer             :: nrefs, ldim(3), iref

    ! Setup
    if (.not. allocated(refs)) then
        write(logfhandle,'(a)') "gauss2Dfit: 'refs' must be allocated"
        return
    end if
    nrefs = size(refs)
    ldim = refs(1)%get_ldim()
    if (ldim(3) > 1) then
        print *, "gauss2Dfit: argument 'refs' must be an array of 2D images!"
        return
    end if
    smpd = refs(1)%get_smpd()
    allocate(centers(2,nrefs), covs(2,2,nrefs), corrs(nrefs), source=0.)

    ! Main Loop
    do iref = 1,nrefs
        call gauss2Dfit(refs(iref), centers(:2,iref), covs(:2,:2,iref), &
            &corrs(iref), fit)
    end do
end subroutine batch_gauss2Dfit_1

! INCLUDES VISUALIZATIONS OF
! Takes in an array of particle images, 'refs.' Fits each particle with a 2D
! multivariate Gaussian distribution.  Outputs same-sized arrays of the fit
! centers, fit covariance matrices, real-space correlations of fits and refs,
! original indices of refs, and fit images. All inout/out arrays are sorted by
! corrs (largest last). Centers and covs are reported as pixels or pixels^2.
! The standard deviations of the fits are the square roots of the eigenvalues
! of covs, and the principal axes are the corresponding eigenvectors. 
subroutine batch_gauss2Dfit_2(refs, centers, covs, corrs, indices, fits)
    type(image), allocatable,   intent(inout)   :: refs(:)
    type(image), allocatable,   intent(out)     :: fits(:)  
    real, allocatable,          intent(out)     :: centers(:,:), covs(:,:,:)
    real, allocatable,          intent(out)     :: corrs(:)
    integer, allocatable,       intent(out)     :: indices(:)
    real, allocatable   :: centers_copy(:,:), covs_copy(:,:,:)
    real                :: smpd
    integer             :: nrefs, ldim(3), iref

    ! Setup
    if (.not. allocated(refs)) then
        write(logfhandle,'(a)') "gauss2Dfit: 'refs' must be allocated"
        return
    end if
    nrefs = size(refs)
    ldim = refs(1)%get_ldim()
    if (ldim(3) > 1) then
        print *, "gauss2Dfit: argument 'refs' must be an array of 2D images!"
        return
    end if
    smpd = refs(1)%get_smpd()
    allocate(centers_copy(2,nrefs), covs_copy(2,2,nrefs), corrs(nrefs), source=0.)
    allocate(indices(nrefs), source=0)
    allocate(fits(nrefs))

    ! Main Loop
    do iref = 1,nrefs
        indices(iref) = iref
        call gauss2Dfit(refs(iref), centers_copy(:2,iref), covs_copy(:2,:2,iref), &
            &corrs(iref), fits(iref))
    end do

    ! Sort output
    call hpsort(corrs, indices)
    allocate(centers(2,nrefs), covs(2,2,nrefs))
    do iref=1, nrefs
        centers(:2,iref) = centers_copy(:2,indices(iref))
        covs(:2,:2,iref) = covs_copy(:2,:2,indices(iref))
    end do
end subroutine batch_gauss2Dfit_2

! Fits a ref image with a 2D multivariate Gaussian by estimating
! the center of mass and the covariance matrix.  Outputs the 
! estimated center and cov matrix, the correlation corr between
! the ref image and the fit image, and the fit image.
subroutine gauss2Dfit(ref, center, cov, corr, fit)
    type(image), intent(inout)   :: ref
    type(image), intent(out)     :: fit  
    real,        intent(out)     :: center(2), cov(2,2), corr
    real, pointer       :: rmat(:,:,:)
    real, allocatable   :: rmat_fit(:,:,:)
    real    :: cov_inv(2,2), displ(2,1), displ_T(1,2), lambda(1,1)
    real    :: smpd, prob, minint, totint, totprob
    integer :: ldim(3), i, j, errflg

    ldim = ref%get_ldim()
    smpd = ref%get_smpd()
    call ref%get_rmat_ptr(rmat)
    allocate(rmat_fit(ldim(1), ldim(2), 1), source = 0.)

    minint = minval(rmat(:ldim(1),:ldim(2),1)) ! So all probs > 0
    totint = sum(rmat(:ldim(1),:ldim(2),1)) ! For normalization
    totint = totint + ldim(1) * ldim(2) * abs(minint)

    ! Find center
    center = 0.
    do j=1, ldim(2)
        do i=1, ldim(1)
            prob = (rmat(i,j,1) + abs(minint)) / totint
            center(1) = center(1) + i * prob
            center(2) = center(2) + j * prob
        end do
    end do

    ! Find covariance matrix
    do j=1, ldim(2)
        do i=1, ldim(1)
            prob = (rmat(i,j,1) + abs(minint)) / totint
            cov(1,1) = cov(1,1) + prob * (i - center(1)) ** 2
            cov(2,2) = cov(2,2) + prob * (j - center(2)) ** 2
            cov(1,2) = cov(1,2) + prob * (i - center(1))*(j - center(2))
        end do
    end do
    cov(2,1) = cov(1,2)
    call matinv(cov, cov_inv, 2, errflg)

    ! Generate fit image
    if (fit%exists()) call fit%kill()
    call fit%new(ldim, smpd)
    totprob = 0.
    do j=1, ldim(2)
        do i=1, ldim(1)
            displ(1,1) = i - center(1)
            displ(2,1) = j - center(2)
            displ_T = transpose(displ)
            lambda = -0.5 * matmul(matmul(displ_T, cov_inv), displ)
            prob = exp(lambda(1, 1))
            totprob = totprob + prob
            rmat_fit(i, j, 1) = prob
        end do
    end do
    ! Normalize and convert probs -> intensities
    rmat_fit = rmat_fit / totprob * totint + minint
    call fit%set_rmat(rmat_fit, .false.)
    corr = fit%real_corr(ref)
end subroutine gauss2Dfit

end module simple_gauss2Dfit