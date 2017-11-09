! combinatorics module
module simple_combinatorics
use simple_math
use simple_rnd,      only: seed_rnd, irnd_uni_pair, irnd_uni
use simple_ran_tabu, only: ran_tabu
use simple_syslib,   only: alloc_errchk
implicit none

contains

    function merge_into_disjoint_set( nprojs, nnn, nnmat, which_nns ) result( disjoint )
        integer, intent(in)  :: nprojs, nnn, nnmat(nprojs,nnn), which_nns(:)
        integer, allocatable :: disjoint(:)
        logical :: isthere(nprojs)
        integer :: cnt, i, inn, ndisjoint
        isthere   = .false.
        ndisjoint = size(which_nns)
        ! count the number of set members
        cnt = 0
        do i=1,ndisjoint
            do inn=1,nnn
                if( isthere(nnmat(which_nns(i),inn)) )then
                    ! element already in disjoint set
                else
                    cnt = cnt + 1
                    isthere(nnmat(which_nns(i),inn)) = .true.
                endif
            end do
        end do
        ! set them
        allocate(disjoint(cnt), stat=alloc_stat )
        if(alloc_stat /= 0) call alloc_errchk('In: combinatorics:: merge_into_disjoint_set ', alloc_stat)
        isthere = .false.
        cnt     = 0
        do i=1,ndisjoint
            do inn=1,nnn
                if( isthere(nnmat(which_nns(i),inn)) )then
                    ! element already in disjoint set
                else
                    cnt = cnt + 1
                    disjoint(cnt) = nnmat(which_nns(i),inn)
                    isthere(disjoint(cnt)) = .true.
                endif
            end do
        end do
    end function merge_into_disjoint_set

    subroutine shc_aggregation( nrepeats, nptcls, labels, consensus, best_ind )
        integer, intent(in)    :: nrepeats, nptcls
        integer, intent(inout) :: labels(nrepeats,nptcls), consensus(nptcls)
        integer, intent(in)    :: best_ind
        integer, parameter     :: MAXITS   = 1000
        real,    parameter     :: TINYTINY = 1e-20               !! Intel warn too small for kind=4 real
        logical, parameter     :: DOPRINT  = .true.
        integer, allocatable   :: counts(:), labels_consensus(:,:), labels_backup(:,:)
        integer :: nlabels, loc(1), rp(2), it, irnd, inds(nrepeats)
        integer :: irep, ilab, iptcl, irestart, restart_winner
        real    :: scores(nrepeats), s, naccepted, norm, score_best
        real    :: score_curr, dists(nrepeats)
        call seed_rnd
        nlabels    = maxval(labels)
        norm       = real((nrepeats-1)*nptcls)
        score_best = 0.0
        allocate(labels_consensus(nrepeats,nptcls),counts(nlabels), stat=alloc_stat )
        if(alloc_stat /= 0) call alloc_errchk('In: combinatorics::shc_aggregation ', alloc_stat)
        labels_backup = labels
        do irestart=1,nrepeats*nlabels
            if( DOPRINT ) write(*,'(a,1x,I5)') '>>> SHC AGGREGATION, RESTART ROUND:', irestart
            ! change the first solution in the row (will affect the greedy initial solution)
            labels = labels_backup
            call change_first( best_ind )
            ! obtain an initial solution using a greedy appraoch
            call greedy_init
            ! initialize scores
            do irep=1,nrepeats
                scores(irep) = score( irep )
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
            score_curr = sum(scores)/norm
            if( score_curr > score_best )then
                score_best       = score_curr
                labels_consensus = labels                      !! reallocation warning
                restart_winner   = irestart
            endif
            if( DOPRINT ) write(*,'(a,1x,f7.2)') '>>> SHC AGGREGATION, SCORE:', score_curr
        end do
        if( DOPRINT ) write(*,'(a,1x,f7.2)') '>>> SHC AGGREGATION, FINAL SCORE:', score_best
        ! report the consensus solution
        do iptcl=1,nptcls
            do ilab=1,nlabels
                counts(ilab) = count(labels_consensus(:,iptcl) == ilab)
            end do
            loc = maxloc(counts)
            consensus(iptcl) = loc(1)
        end do

        if(allocated(labels_consensus)) then
            deallocate(labels_consensus, stat=alloc_stat )
            if(alloc_stat /= 0) call alloc_errchk('In: combinatorics::shc_aggregation deallocating labels_consensus ', alloc_stat)
        end if
        if(allocated(counts)) then
            deallocate(counts, stat=alloc_stat )
            if(alloc_stat /= 0) call alloc_errchk('In: combinatorics::shc_aggregation deallocating counts', alloc_stat)
        end if
        if(allocated(labels_backup)) then
            deallocate(labels_backup, stat=alloc_stat )
            if(alloc_stat /= 0) call alloc_errchk('In: combinatorics::shc_aggregation deallocating labels_backup', alloc_stat)
        end if
        contains

            subroutine change_first( irep )
                integer, intent(in) :: irep
                integer :: tmp(nptcls)
                if( irep == 1 ) return
                tmp            = labels(1,:)
                labels(1,:)    = labels(irep,:)
                labels(irep,:) = tmp
            end subroutine change_first

            subroutine greedy_init
                integer :: irep, iswap, jswap, swap_best(2)
                real    :: smax, s
                ! (1) we fix the first solution in its original configuration
                ! (2) we loop over the remaining solutions
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

            real function score_pairs( n )
                integer, intent(in) :: n
                integer :: irep, jrep, npairs
                npairs = n*(n-1)/2
                score_pairs = 0.
                do irep=1,n-1
                    do jrep=irep+1,n
                        score_pairs = score_pairs + count(labels(irep,:) == labels(jrep,:))
                    end do
                end do
                score_pairs = score_pairs/real(npairs)
            end function score_pairs

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
                do iptcl=1,nptcls
                    if( labels(irep,iptcl) == rp(2) ) labels(irep,iptcl) = nlabels + 1
                    if( labels(irep,iptcl) == rp(1) ) labels(irep,iptcl) = rp(2)
                end do
                do iptcl=1,nptcls
                    if( labels(irep,iptcl) == nlabels + 1 ) labels(irep,iptcl) = rp(1)
                end do
            end subroutine swap_labels

            real function score( irep )
                integer, intent(in) :: irep
                integer :: jrep
                score = 0.
                do jrep=1,nrepeats
                    if( jrep /= irep ) score = score + count(labels(irep,:) == labels(jrep,:))
                end do
            end function score

    end subroutine shc_aggregation

end module simple_combinatorics
