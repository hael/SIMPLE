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
        use gnufor2, only: hist
        class(polarizer), target, intent(inout) :: cavg_imgs(:)
        integer, allocatable,     intent(inout) :: centers(:), labels(:)
        type(polarft_corrcalc) :: pftcc
        type(aff_prop)         :: aprop
        real,    allocatable   :: corrs(:), corrmat_comlin(:,:), corrs_median(:)
        real,    allocatable   :: corrmat(:,:), corrs_rot(:), specmat(:,:)
        logical, allocatable   :: mask(:,:), mask_median(:), mask_otsu(:)
        integer, allocatable   :: order(:), nloc(:)
        integer :: icls, i, j, ncls, pop1, pop2, loc(1), nsel
        real    :: simsum
        write(logfhandle,'(A)') '>>> PREPARING CLASS AVERAGES'
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
            call cavg_imgs(icls)%ifft()
        end do
        !$omp end parallel do
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        ! some allocations needed
        nsel = ceiling(0.1 * real(ncls))
        allocate(nloc(nsel), source=0)
        allocate(corrmat(ncls,ncls), specmat(ncls,ncls), corrmat_comlin(ncls,ncls),&
        &corrs_median(ncls), corrs_rot(pftcc%get_nrots()), source=-1.)
        allocate(mask(ncls,ncls), mask_median(ncls), source=.true.)
        allocate(order(ncls), source=0)
        write(logfhandle,'(A)') '>>> CALCULATING CORRELATION MATRICES'
        !$omp parallel do default(shared) private(i,j,corrs_rot,loc) schedule(dynamic) proc_bind(close)
        do i = 1, ncls - 1
            do j = i + 1, ncls
                ! 2D correlations
                call pftcc%gencorrs(i,j,corrs_rot)
                loc = maxloc(corrs_rot)
                corrmat(i,j) = corrs_rot(loc(1))
                corrmat(j,i) = corrmat(i,j)
                ! 2D specscores
                specmat(i,j) = pftcc%specscore(i,j,loc(1))
                specmat(j,i) = specmat(i,j)
                ! Common line orrelations
                ! corrmat_comlin(i,j) = pftcc%genmaxcorr_comlin(i,j)
                corrmat_comlin(i,j) = pftcc%genmaxspecscore_comlin(i,j)
                corrmat_comlin(j,i) = corrmat_comlin(i,j)
                mask(i,j)    = .true.
            enddo
            corrmat(i,i)        = 1.
            specmat(i,i)        = 1.
            corrmat_comlin(i,i) = 1.
        enddo
        !$omp end parallel do
        ! take care of the last diagonal element
        corrmat(ncls,ncls)        = 1.
        specmat(ncls,ncls)        = 1.
        corrmat_comlin(ncls,ncls) = 1.
        ! histogram, remove first lime of command file and execute 'cat command_file.txt | gnuplot --persist'
        ! corrs = pack(corrmat_comlin, mask)
        ! call hist(corrs, ncls/5)


        do i = 1, ncls
            mask_median    = .true.
            mask_median(i) = .false. ! to remove the diagonal element
            corrs = pack(specmat(i,:), mask_median)


            nloc = maxnloc(corrs, nsel)
            corrs_median(i) = 0.
            do j = 1, nsel
                corrs_median(i) = corrs_median(i) + corrs(nloc(j))
            end do
            corrs_median(i) = corrs_median(i) / real(nsel)


            ! corrs_median(i) = median_nocopy(corrs)
            ! corrs_median(i) = maxval(corrs)

        end do

        !>>>>>>>>>>> OTSU DOESNT CUT IT, NEED MORE BIAS TOWARDS THE GOOD ONES
        call otsu(corrs_median, mask_otsu)
        pop1 = count(      mask_otsu)
        pop2 = count(.not. mask_otsu)
        write(*,*) 'average corr cluster 1: ', sum(corrs_median, mask=      mask_otsu) / real(pop1), ' pop ', pop1
        write(*,*) 'average corr cluster 2: ', sum(corrs_median, mask=.not. mask_otsu) / real(pop2), ' pop ', pop2
        pop1 = 0
        pop2 = 0
        do i = 1, ncls
            if( mask_otsu(i) )then
                pop1 = pop1 + 1
                call cavg_imgs(i)%write('good.mrc', pop1)
            else
                pop2 = pop2 + 1
                call cavg_imgs(i)%write('bad.mrc',  pop2)
            endif
        end do
        !<<<<<<<<<<< OTSU DOESNT CUT IT, NEED MORE BIAS TOWARDS THE GOOD ONES

        ! rank
        order = (/(i,i=1,ncls)/)
        call hpsort(corrs_median, order)
        do i = 1, ncls
            print *, corrs_median(i)
            call cavg_imgs(order(i))%write('ranked.mrc', i)
        end do



        write(logfhandle,'(A)') '>>> CLUSTERING WITH AFFINITY PROPAGATION'
        call aprop%new(ncls, corrmat_comlin)
        call aprop%propagate(centers, labels, simsum)



        ! deallocate

        call pftcc%kill
        call aprop%kill


    end subroutine cluster_cavgs

end module simple_cluster_cavgs
