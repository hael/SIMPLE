! Module for fitting the intensity distributions of 2D particle image arrays
! with 2D multivariate Gaussians.  There are two implementations.
! The second implementation is slower but includes sorting by correlation
! and an output array of 2D images displaying the Gaussian fits.
module simple_gauss2Dfit

include 'simple_lib.f08'
use simple_picker_utils, only: picker_utils
use simple_image,        only: image
implicit none

public  :: gauss2Dfit
private

interface gauss2Dfit
    module procedure gauss2Dfit_1, gauss2Dfit_2
end interface gauss2Dfit

contains

! FASTER IMPLEMENTATION BUT NO VISUALIZATION OF FITS OR SORTING
! Takes in an array of particle images, 'refs.' Fits each particle with a 2D
! multivariate Gaussian distribution.  Outputs same-sized arrays of the fit
! centers, fit covariance matrices, and real-space correlations of fits and
! refs. Centers and covs are reported as pixels or pixels^2. The standard
! deviations of the fits are the square roots of the eigenvalues
! of covs, and the principal axes are the corresponding eigenvectors. 
subroutine gauss2Dfit_1(refs, centers, covs, corrs)
    type(image), allocatable,   intent(inout)   :: refs(:)
    real, allocatable,          intent(out)     :: centers(:,:), covs(:,:,:)
    real, allocatable,          intent(out)     :: corrs(:)
    type(image)        :: fit
    real, pointer      :: rmat(:,:,:)
    real, allocatable  :: rmat_fit(:,:,:)
    real               :: minint, totint, prob, totprob, beta(1,1), mu(2), smpd
    real               :: sigma(2,2), sigma_inv(2,2), displ(2,1), displ_T(1,2)
    integer            :: nrefs, ldim(3), iref, i, j, errflg

    ! Setup
    if (.not. allocated(refs)) then
        write(logfhandle,'(a)') "gauss2Dfit: 'refs' must be allocated"
        return
    end if
    nrefs = size(refs)
    ldim = refs(1)%get_ldim()
    if (ldim(3) > 1) then
        write(logfhandle,'(a)') "gauss2Dfit: 'refs' must be array of 2D images"
        return
    end if
    allocate(centers(2,nrefs), covs(2,2,nrefs), corrs(nrefs), source=0.)
    allocate(rmat_fit(ldim(1), ldim(2), 1), source = 0.)

    ! Main Loop
    do iref = 1,nrefs
        call refs(iref)%get_rmat_ptr(rmat)
        minint = minval(rmat(:ldim(1),:ldim(2),1)) ! So all probs > 0
        if (minint > 0.) minint = 0.
        totint = sum(rmat(:ldim(1),:ldim(2),1)) ! For normalization
        totint = totint + ldim(1) * ldim(2) * abs(minint)
    
        ! Find center
        mu(:2) = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                prob = (rmat(i,j,1) + abs(minint)) / totint
                mu(1) = mu(1) + i * prob
                mu(2) = mu(2) + j * prob
            end do
        end do
        centers(:2,iref) = mu(:2)
    
        ! Find 2x2 covariance matrix
        sigma(:2,:2) = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                prob = (rmat(i,j,1) + abs(minint)) / totint
                sigma(1,1) = sigma(1,1) + prob * (i - mu(1)) ** 2
                sigma(2,2) = sigma(2,2) + prob * (j - mu(2)) ** 2
                sigma(1,2) = sigma(1,2) + prob * (i - mu(1))*(j - mu(2))
            end do
        end do
        sigma(2,1) = sigma(1,2) ! The matrix is symmetric
        covs(:2,:2,iref) = sigma(:2,:2)

        ! Generate fit image
        call fit%new(ldim, smpd)
        call matinv(sigma, sigma_inv, 2, errflg)
        totprob = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                displ(1,1) = i - mu(1)
                displ(2,1) = j - mu(2)
                displ_T = transpose(displ)
                beta = -0.5 * matmul(matmul(displ_T, sigma_inv), displ)
                prob = exp(beta(1, 1))
                totprob = totprob + prob
                rmat_fit(i, j, 1) = prob
            end do
        end do
        ! Normalize and convert probs -> intensities
        rmat_fit = rmat_fit / totprob * totint + minint
        call fit%set_rmat(rmat_fit, .false.)
        corrs(iref) = fit%real_corr(refs(iref))
        call fit%kill()
    end do
end subroutine gauss2Dfit_1

! SLOWER IMPLEMENTATION THAT ALLOWS VISUALIZATION OF FITS AND SORTING BY CORR.
! Takes in an array of particle images, 'refs.' Fits each particle with a 2D
! multivariate Gaussian distribution.  Outputs same-sized arrays of the fit
! centers, fit covariance matrices, real-space correlations of fits and refs,
! original indices of refs, and fit images. All inout/out arrays are sorted by
! corrs (largest last). Centers and covs are reported as pixels or pixels^2.
! The standard deviations of the fits are the square roots of the eigenvalues
! of covs, and the principal axes are the corresponding eigenvectors. 
subroutine gauss2Dfit_2(refs, centers, covs, corrs, indices, fits)
    type(image), allocatable,   intent(inout)   :: refs(:)
    type(image), allocatable,   intent(out)     :: fits(:)  
    real, allocatable,          intent(out)     :: centers(:,:), covs(:,:,:)
    real, allocatable,          intent(out)     :: corrs(:)
    integer, allocatable,       intent(out)     :: indices(:)
    real, pointer       :: rmat(:,:,:)
    real, allocatable   :: rmat_fit(:,:,:), centers_copy(:,:), covs_copy(:,:,:)
    real                :: smpd, minint, totint, prob, totprob, beta(1,1), mu(2)
    real                :: sigma(2,2), sigma_inv(2,2), displ(2,1), displ_T(1,2)
    integer             :: nrefs, ldim(3), iref, i, j, errflg

    ! Setup
    if (.not. allocated(refs)) then
        print *, "gauss2Dfit: argument 'refs' must be allocated!"
        stop
    end if
    nrefs = size(refs)
    ldim = refs(1)%get_ldim()
    if (ldim(3) > 1) then
        print *, "gauss2Dfit: argument 'refs' must be an array of 2D images!"
        stop
    end if
    smpd = refs(1)%get_smpd()
    allocate(centers_copy(2,nrefs), covs_copy(2,2,nrefs), corrs(nrefs), source=0.)
    allocate(indices(nrefs), source=0)
    allocate(fits(nrefs))
    allocate(rmat_fit(ldim(1), ldim(2), 1), source = 0.)

    ! Main Loop
    do iref = 1,nrefs
        indices(iref) = iref
        call refs(iref)%get_rmat_ptr(rmat)
        minint = minval(rmat(:ldim(1),:ldim(2),1)) ! So all probs > 0
        if (minint > 0.) minint = 0.
        totint = sum(rmat(:ldim(1),:ldim(2),1)) ! For normalization
        totint = totint + ldim(1) * ldim(2) * abs(minint)
    
        ! Find center
        mu(:2) = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                prob = (rmat(i,j,1) + abs(minint)) / totint
                mu(1) = mu(1) + i * prob
                mu(2) = mu(2) + j * prob
            end do
        end do
        centers_copy(:2,iref) = mu(:2)
    
        ! Find 2x2 covariance matrix
        sigma(:2,:2) = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                prob = (rmat(i,j,1) + abs(minint)) / totint
                sigma(1,1) = sigma(1,1) + prob * (i - mu(1)) ** 2
                sigma(2,2) = sigma(2,2) + prob * (j - mu(2)) ** 2
                sigma(1,2) = sigma(1,2) + prob * (i - mu(1))*(j - mu(2))
            end do
        end do
        sigma(2,1) = sigma(1,2) ! The matrix is symmetric
        covs_copy(:2,:2,iref) = sigma(:2,:2)
    
        ! Generate fit image
        call fits(iref)%new(ldim, smpd)
        call matinv(sigma, sigma_inv, 2, errflg)
        totprob = 0.
        do j=1, ldim(2)
            do i=1, ldim(1)
                displ(1,1) = i - mu(1)
                displ(2,1) = j - mu(2)
                displ_T = transpose(displ)
                beta = -0.5 * matmul(matmul(displ_T, sigma_inv), displ)
                prob = exp(beta(1, 1))
                totprob = totprob + prob
                rmat_fit(i, j, 1) = prob
            end do
        end do
        ! Normalize and convert probs -> intensities
        rmat_fit = rmat_fit / totprob * totint + minint
        call fits(iref)%set_rmat(rmat_fit, .false.)
        corrs(iref) = fits(iref)%real_corr(refs(iref))
    end do

    ! Sort output
    call hpsort(corrs, indices)
    allocate(centers(2,nrefs), covs(2,2,nrefs))
    do iref=1, nrefs
        centers(:2,iref) = centers_copy(:2,indices(iref))
        covs(:2,:2,iref) = covs_copy(:2,:2,indices(iref))
    end do
end subroutine gauss2Dfit_2

end module simple_gauss2Dfit