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

    subroutine aggregate( nptcls, nrepeats, labels, oris_out )
        use simple_shc_cluster
        use simple_oris, only: oris
        integer,              intent(in)    :: nptcls, nrepeats, labels(nrepeats,nptcls)
        class(oris),          intent(inout) :: oris_out
        type(shc_cluster) :: cluster
        type(oris) :: myoris
        real,    allocatable :: simmat(:,:)
        real    :: sim
        integer :: inds(nptcls)
        integer :: iptcl, jptcl, irep, icnt, jcnt, nexisting, nlabels, label
        logical :: zeroes(nptcls)
        ! sanity check
        if( oris_out%get_noris() /= nptcls ) stop 'nonconforming number of&
        &oris (oris_out); simple_combinatorics :: aggregate'
        zeroes  = labels(1,:) == 0
        nlabels = maxval(labels)
        do iptcl=1,nptcls
            icnt = count(labels(:,iptcl) == 0)
            if( icnt == 0 .or. icnt == nrepeats )then
                ! all good
            else
                stop 'the pattern of zeroes (exclusions) not consistent&
                &between solutions; simple_combinatorics :: aggregate'
            endif
        end do
        ! similarity matrix init
        inds      = 0
        nexisting = 0
        do iptcl=1,nptcls
            if( labels(1,iptcl) > 0 )then
                nexisting = nexisting + 1
                inds(iptcl) = nexisting
            endif
        enddo
        allocate( simmat(nexisting,nexisting) )
        simmat = 0.
        icnt = 0
        do iptcl=1,nptcls-1
            if(inds(iptcl)==0)cycle
            icnt = icnt + 1
            jcnt = 0
            do jptcl=iptcl+1,nptcls
                if(inds(jptcl)==0)cycle
                jcnt = jcnt + 1
                do irep = 1,nrepeats
                    if( labels(irep,iptcl) == labels(irep,jptcl) )simmat(icnt,jcnt) =simmat(icnt,jcnt)+1.
                enddo
            enddo
        enddo
        simmat = simmat / real(nrepeats)
        forall(iptcl=1:nexisting)simmat(iptcl,iptcl)=1.
        do iptcl=1,nexisting-1
            do jptcl=iptcl+1,nexisting
                simmat(jptcl,iptcl) =simmat(iptcl,jptcl)
            enddo
        enddo        
        ! clustering
        call myoris%new(nexisting)
        call cluster%new( nexisting, nlabels, simmat, myoris )
        call cluster%shc( .true.,'state',sim )
        call cluster%kill
        ! reporting
        icnt = 0
        do iptcl=1,nptcls
            if(inds(iptcl) > 0)then
                icnt = icnt + 1
                label = nint(myoris%get(icnt,'state'))
                call oris_out%set(iptcl,'state',real(label))
            endif
        enddo
        deallocate(simmat)
    end subroutine aggregate

end module simple_combinatorics