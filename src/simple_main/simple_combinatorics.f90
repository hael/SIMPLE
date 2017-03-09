module simple_combinatorics
use simple_ran_tabu, only: ran_tabu
use simple_jiffys,   only: progress
implicit none

contains

    function diverse_labeling( nptcls, nlabels, ndiverse ) result( configs_diverse )
        !$ use omp_lib
        !$ use omp_lib_kinds
        integer, intent(in)  :: nptcls, nlabels, ndiverse
        integer, allocatable :: configs_diverse(:,:), configs_trial(:,:), tmp(:)
        real,    allocatable :: dists(:,:), rmat(:,:)
        type(ran_tabu)       :: rt
        integer :: isample, idiv, curr_pop, loc(2), nsamples, backcnt
        write(*,'(a)') '>>> GENERATING DIVERSE LABELING'
        nsamples = 10*ndiverse*ndiverse
        allocate(configs_diverse(ndiverse,nptcls), configs_trial(nsamples,nptcls),&
        rmat(nsamples,nptcls), dists(nsamples,ndiverse), tmp(nptcls))
        curr_pop = 1
        ! generate random configurations
        call random_number(rmat)
        !$omp parallel workshare
        configs_trial = ceiling(rmat*real(nlabels))
        where(configs_trial < 1       ) configs_trial = 1
        where(configs_trial > nlabels ) configs_trial = nlabels
        !$omp end parallel workshare
        ! set the first solution
        curr_pop = 1
        backcnt  = nsamples
        configs_diverse(curr_pop,:) = configs_trial(curr_pop,:)
        do while(curr_pop < ndiverse)
            ! calculate Hamming distances
            dists = 0.
            !$omp parallel do schedule(auto) default(shared) private(isample,idiv)
            do isample=1,10*ndiverse
                do idiv=1,curr_pop
                    dists(isample,idiv) = real(count(configs_diverse(idiv,:) /= configs_trial(isample,:)))
                end do
            end do
            !$omp end parallel do
            ! find the most diverse labeling solution
            loc = maxloc(dists)
            curr_pop = curr_pop + 1
            configs_diverse(curr_pop,:) = configs_trial(loc(1),:)
            ! swap with last in the trial set (using the backward counter)
            tmp = configs_trial(backcnt,:)
            configs_trial(backcnt,:) = configs_trial(loc(1),:)
            configs_trial(loc(1),:) = tmp
            backcnt = backcnt - 1
        end do
        deallocate(configs_trial, rmat, dists, tmp)
    end function diverse_labeling

end module simple_combinatorics