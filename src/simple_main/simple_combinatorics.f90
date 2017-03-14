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

    ! subroutine aggregate( nptcls, nrepeats, labels, oris_out )
    !     use simple_shc_cluster
    !     use simple_oris, only: oris
    !     integer,              intent(in)    :: nptcls, nrepeats, labels(nrepeats,nptcls)
    !     class(oris),          intent(inout) :: oris_out
    !     type(shc_cluster) :: cluster
    !     type(oris) :: myoris
    !     real,    allocatable :: simmat(:,:)
    !     real    :: sim
    !     integer :: inds(nptcls)
    !     integer :: iptcl, jptcl, irep, icnt, jcnt, nexisting, nlabels, label
    !     logical :: zeroes(nptcls)
    !     ! sanity check
    !     if( oris_out%get_noris() /= nptcls ) stop 'nonconforming number of&
    !     &oris (oris_out); simple_combinatorics :: aggregate'
    !     zeroes  = labels(1,:) == 0
    !     nlabels = maxval(labels)
    !     do iptcl=1,nptcls
    !         icnt = count(labels(:,iptcl) == 0)
    !         if( icnt == 0 .or. icnt == nrepeats )then
    !             ! all good
    !         else
    !             stop 'the pattern of zeroes (exclusions) not consistent&
    !             &between solutions; simple_combinatorics :: aggregate'
    !         endif
    !     end do
    !     ! similarity matrix init
    !     inds      = 0
    !     nexisting = 0
    !     do iptcl=1,nptcls
    !         if( labels(1,iptcl) > 0 )then
    !             nexisting = nexisting + 1
    !             inds(iptcl) = nexisting
    !         endif
    !     enddo
    !     allocate( simmat(nexisting,nexisting) )
    !     simmat = 0.
    !     icnt = 0
    !     do iptcl=1,nptcls-1
    !         if(inds(iptcl)==0)cycle
    !         icnt = icnt + 1
    !         jcnt = 0
    !         do jptcl=iptcl+1,nptcls
    !             if(inds(jptcl)==0)cycle
    !             jcnt = jcnt + 1
    !             do irep = 1,nrepeats
    !                 if( labels(irep,iptcl) == labels(irep,jptcl) )simmat(icnt,jcnt) =simmat(icnt,jcnt)+1.
    !             enddo
    !         enddo
    !     enddo
    !     simmat = simmat / real(nrepeats)
    !     forall(iptcl=1:nexisting)simmat(iptcl,iptcl)=1.
    !     do iptcl=1,nexisting-1
    !         do jptcl=iptcl+1,nexisting
    !             simmat(jptcl,iptcl) =simmat(iptcl,jptcl)
    !         enddo
    !     enddo        
    !     ! clustering
    !     call myoris%new(nexisting)
    !     call cluster%new( nexisting, nlabels, simmat, myoris )
    !     call cluster%shc( .true.,'state',sim )
    !     call cluster%kill
    !     ! reporting
    !     icnt = 0
    !     do iptcl=1,nptcls
    !         if(inds(iptcl) > 0)then
    !             icnt = icnt + 1
    !             label = nint(myoris%get(icnt,'state'))
    !             call oris_out%set(iptcl,'state',real(label))
    !         endif
    !     enddo
    !     deallocate(simmat)
    ! end subroutine aggregate

    subroutine shc_aggregation( nrepeats, nptcls, labels, consensus, score_out )
        use simple_rnd, only: irnd_uni, irnd_uni_pair
        integer, intent(in)    :: nrepeats, nptcls
        integer, intent(inout) :: labels(nrepeats,nptcls), consensus(nptcls)
        real,    intent(out)   :: score_out
        integer, parameter     :: MAXITS=1000
        real,    parameter     :: TINYTINY=1e-40
        integer, allocatable   :: counts(:)
        integer :: nlabels, loc(1), rp(2), it, irnd, irep, ilab, iptcl
        real    :: scores(nrepeats), s, naccepted, norm
        nlabels = maxval(labels)
        norm    = real((nrepeats-1)*nptcls)
        ! obtain an initial solution using a greedy appraoch
        call greedy_init
        ! initialize scores
        do irep=1,nrepeats
            scores = score( irep )
        end do
        ! stochastic hill-climbing
        naccepted = 0.
        do it=1,MAXITS
            ! pick a random solution
            irnd = irnd_uni( nrepeats )
            ! swap a random pair of labels
            call swap_labels( irnd )
            ! evaluate the score
            s = score(irnd)
            if( s >= scores(irnd) )then
                ! solution accepted, update scores
                scores(irnd) = s
                do irep=1,nrepeats
                    if( irep == irnd ) cycle
                    scores = score( irep )
                end do
                naccepted = 0.95*naccepted + 0.05
            else
                ! swap back
                call swap_labels( irnd, rp )
                naccepted = 0.95*naccepted
            endif
            if( naccepted <= TINYTINY ) exit
        end do
        score_out = sum(scores)/norm
        write(*,'(a,1x,f7.2)') '>>> SHC AGGREGATION, FINAL SCORE:', score_out
        ! report the consensus solution
        allocate(counts(nlabels))
        do iptcl=1,nptcls
            do ilab=1,nlabels
                counts(ilab) = count(labels(:,iptcl) == ilab)
            end do
            loc = maxloc(counts)
            consensus(iptcl) = loc(1)
        end do

        contains

            subroutine greedy_init
                integer :: irep, iswap, jswap, swap_best(2)
                real    :: smax, s
                ! (1) we fix the first solution in its original configuration
                ! (2) we loop over the remainding solutions
                do irep=2,nrepeats
                    swap_best = 0
                    smax      = 0.
                    ! (3) we evaluate all possible swaps for the current labeling solution (irep)
                    do iswap=1,nlabels-1
                        do jswap=iswap+1,nlabels
                            call swap_labels( irep, [iswap,jswap] )
                            ! the score is now over all pairs of solutions included in the set
                            ! this is the greedy aspect of the initialisation
                            s = score_pairs( irep )
                            if( s > smax )then
                                ! search update
                                swap_best = [iswap,jswap]
                                smax      = s
                            endif
                            ! put back labels
                            call swap_labels( irep, [iswap,jswap] )
                        end do
                    end do
                    ! (4) store the best greedy solution obtained
                    call swap_labels( irep, swap_best )
                end do
            end subroutine greedy_init

            subroutine swap_labels( irep, rp_in )
                integer,           intent(in) :: irep
                integer, optional, intent(in) :: rp_in(2)
                integer :: iptcl
                if( present(rp_in) )then
                    rp = rp_in
                else
                    do 
                        rp = irnd_uni_pair(nlabels)
                        if( count(labels(irep,:) == rp(1)) > 0 .and.&
                            count(labels(irep,:) == rp(2)) > 0 ) exit
                    end do
                endif
                !$omp parallel default(shared) private(iptcl)
                !$omp do schedule(auto)
                do iptcl=1,nptcls
                    if( labels(irep,iptcl) == rp(2) ) labels(irep,iptcl) = nlabels + 1
                    if( labels(irep,iptcl) == rp(1) ) labels(irep,iptcl) = rp(2)
                end do
                !$omp end do nowait
                !$omp do
                do iptcl=1,nptcls
                    if( labels(irep,iptcl) == nlabels + 1 ) labels(irep,iptcl) = rp(1)
                end do
                !$omp end do nowait
                !$omp end parallel
            end subroutine swap_labels

            real function score( irep )
                integer, intent(in) :: irep
                integer :: jrep
                score = 0.
                !$omp parallel do default(shared) private(jrep) reduction(+:score) schedule(auto)
                do jrep=1,nrepeats
                    if( jrep /= irep ) score = score + count(labels(irep,:) == labels(jrep,:))
                end do
                !$omp end parallel do
            end function score

            real function score_pairs( n )
                integer, intent(in) :: n
                integer :: irep, jrep, npairs
                npairs = n*(n-1)/2
                score_pairs = 0.
                !$omp parallel do default(shared) private(irep,jrep) reduction(+:score_pairs) schedule(auto)
                do irep=1,n-1
                    do jrep=irep+1,n
                        score_pairs = score_pairs + count(labels(irep,:) == labels(jrep,:))
                    end do
                end do
                !$omp end parallel do
                score_pairs = score_pairs/real(npairs)
            end function score_pairs

    end subroutine shc_aggregation

end module simple_combinatorics