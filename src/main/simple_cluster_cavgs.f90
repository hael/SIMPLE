module simple_cluster_cavgs
include 'simple_lib.f08'
use simple_polarizer,        only: polarizer
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_parameters,       only: params_glob
implicit none
#include "simple_local_flags.inc"

contains

    subroutine cluster_cavgs(  cavg_imgs, centers, labels )
        use simple_aff_prop, only: aff_prop
        class(polarizer), target, intent(inout) :: cavg_imgs(:)
        integer, allocatable,     intent(inout) :: centers(:), labels(:)
        type(polarft_corrcalc) :: pftcc
        type(aff_prop)         :: aprop
        real, allocatable      :: corrs(:), corrmat(:,:), tmparr(:)
        integer :: icls, i, j, ncls
        real    :: corr_med, simsum
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES FOR MATCHING'
        ncls = size(cavg_imgs)
        do icls = 1, ncls
            if( cavg_imgs(icls)%is_ft() )                 THROW_HARD('cavgs assumed not be FTed; cluster_cavgs')
            if( .not. params_glob%l_prenormpremsk )       call cavg_imgs(icls)%norm
            if( .not. params_glob%l_prenormpremsk )       call cavg_imgs(icls)%mask(params_glob%msk, 'soft')
            if( params_glob%boxmatch /= params_glob%box ) call cavg_imgs(icls)%clip_inplace([params_glob%boxmatch,params_glob%boxmatch,1])
        end do
        write(logfhandle,'(A)') '>>> PREPARING REFERENCES IN POLAR REPRESENTATION'
        ! create the polarft_corrcalc object
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%boxmatch, params_glob%smpd))
        params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%boxmatch, params_glob%smpd)
        params_glob%kstop      = params_glob%kfromto(2)
        call pftcc%new(ncls, [1,ncls])
        ! initialize polarizer for the first image, then copy it to the rest
        call cavg_imgs(1)%init_polarizer(pftcc, params_glob%alpha)
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            if( icls /= 1 ) call cavg_imgs(icls)%copy_polarizer(cavg_imgs(1))
            ! gridding prep
            if( params_glob%griddev.eq.'yes' )then
                call cavg_imgs(icls)%div_by_instrfun
            endif
            ! move to Fourier space
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
        ! find median correlation
        tmparr   = pack(corrmat, mask=.true.)
        corr_med = median_nocopy(tmparr)
        deallocate(tmparr)
        write(logfhandle,'(A)') '>>> PERFORMING CLUSTERING WITH AFFINITY PROPAGATION'
        call aprop%new(ncls, corrmat, pref=1.5*corr_med)
        call aprop%propagate(centers, labels, simsum)
        call aprop%kill
        deallocate(corrmat)
    end subroutine cluster_cavgs

end module simple_cluster_cavgs
