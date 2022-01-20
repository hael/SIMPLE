module simple_cluster_cavgs
include 'simple_lib.f08'
use simple_polarizer,        only: polarizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_parameters,       only: params_glob
use simple_aff_prop,         only: aff_prop
implicit none
#include "simple_local_flags.inc"

contains
    
    subroutine cluster_cavgs( cavg_imgs, centers, labels )
        class(polarizer), target, intent(inout) :: cavg_imgs(:)
        integer, allocatable,     intent(inout) :: centers(:), labels(:)
        type(polarft_corrcalc) :: pftcc
        type(aff_prop)         :: aprop
        real, allocatable      :: corrs(:), corrmat(:,:)
        integer :: icls, i, j, ncls
        real    :: simsum
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES FOR MATCHING'
        ncls = size(cavg_imgs)
        ! create the polarft_corrcalc object
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
        params_glob%kfromto(2) =        calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
        params_glob%kstop      = params_glob%kfromto(2)
        call pftcc%new(ncls, [1,ncls])
        ! initialize polarizer for the first image, then copy it to the rest
        call cavg_imgs(1)%init_polarizer(pftcc, params_glob%alpha)
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            if( icls /= 1 ) call cavg_imgs(icls)%copy_polarizer(cavg_imgs(1))
            call cavg_imgs(icls)%fft()
            call cavg_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true.) ! 2 polar coords
            call pftcc%cp_even2odd_ref(icls)
            call pftcc%cp_even_ref2ptcl(icls, icls)
        end do
        !$omp end parallel do
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        write(logfhandle,'(A)') '>>> CALCULATING CORRELATION MATRIX'
        allocate(corrs(pftcc%get_nrots()), corrmat(ncls,ncls), source=-1.)
        !$omp parallel do default(shared) private(i,j,corrs) schedule(dynamic) proc_bind(close)
        do i = 1, ncls - 1
            do j = i + 1, ncls
                call pftcc%gencorrs(i,j,corrs)
                corrmat(i,j) = maxval(corrs)
                corrmat(j,i) = corrmat(i,j)
            enddo
            corrmat(i,i) = 1.
        enddo
        !$omp end parallel do
        call pftcc%kill
        write(logfhandle,'(A)') '>>> PERFORMING CLUSTERING WITH AFFINITY PROPAGATION'
        call aprop%new(ncls, corrmat)
        call aprop%propagate(centers, labels, simsum)
        call aprop%kill
        deallocate(corrmat)
    end subroutine cluster_cavgs

    subroutine cluster_cavgs_comlin(  cavg_imgs, centers, labels )
        class(polarizer), target, intent(inout) :: cavg_imgs(:)
        integer, allocatable,     intent(inout) :: centers(:), labels(:)
        type(polarft_corrcalc) :: pftcc
        type(aff_prop)         :: aprop
        real, allocatable      :: corrs(:), corrmat(:,:)
        integer :: icls, i, j, ncls
        real    :: simsum
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES FOR COMMON-LINE CALCULATION'
        ncls = size(cavg_imgs)
        ! create the polarft_corrcalc object
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
        params_glob%kfromto(2) =        calc_fourier_index(params_glob%lp, params_glob%box, params_glob%smpd)
        params_glob%kstop      = params_glob%kfromto(2)
        call pftcc%new(ncls, [1,1])
        ! initialize polarizer for the first image, then copy it to the rest
        call cavg_imgs(1)%init_polarizer(pftcc, params_glob%alpha)
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            if( icls /= 1 ) call cavg_imgs(icls)%copy_polarizer(cavg_imgs(1))
            ! move to Fourier space
            call cavg_imgs(icls)%fft()
            call cavg_imgs(icls)%polarize(pftcc, icls, isptcl=.false., iseven=.true.) ! 2 polar coords
        end do
        !$omp end parallel do
        write(logfhandle,'(A)') '>>> CALCULATING COMMON-LINE CORRELATION MATRIX'
        allocate(corrs(pftcc%get_nrots()), corrmat(ncls,ncls), source=-1.)
        !$omp parallel do default(shared) private(i,j) schedule(dynamic) proc_bind(close)
        do i = 1, ncls - 1
            do j = i + 1, ncls
                corrmat(i,j) = pftcc%genmaxcorr_comlin(i,j)
                corrmat(j,i) = corrmat(i,j)
            enddo
            corrmat(i,i) = 1.
        enddo
        !$omp end parallel do
        call pftcc%kill
        write(logfhandle,'(A)') '>>> PERFORMING CLUSTERING WITH AFFINITY PROPAGATION'
        call aprop%new(ncls, corrmat)
        call aprop%propagate(centers, labels, simsum)
        call aprop%kill
        deallocate(corrmat)
    end subroutine cluster_cavgs_comlin

end module simple_cluster_cavgs
