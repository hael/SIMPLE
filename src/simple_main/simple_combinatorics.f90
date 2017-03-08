module simple_combinatorics
use simple_stat,        only: rand_index
use simple_rnd,         only: irnd_uni
use simple_ran_tabu,    only: ran_tabu
use simple_jiffys,      only: alloc_err
use simple_shc_cluster, only: shc_cluster
use simple_oris,        only: oris
use simple_jiffys,      only: progress
implicit none


contains

    function balanced_diverse_labeling( nptcls, nlabels, ndiverse ) result( configs )
        integer, intent(in)  :: nptcls, nlabels, ndiverse
        integer, allocatable :: configs(:,:)
        integer :: step, iconf, labcnt, i, istart, istop
        allocate(configs(ndiverse,nptcls))
        step = 1
        do iconf=1,ndiverse
            labcnt = 1
            do i=1,nptcls,step
                istart = i
                istop  = min(nptcls, i + step - 1)
                configs(iconf,istart:istop) = labcnt
                labcnt = labcnt + 1
                if( labcnt > nlabels ) labcnt = 1
            end do
            if( istop < nptcls )then
                do i=istop,nptcls,1
                    configs(iconf,i) = labcnt
                    labcnt = labcnt + 1
                    if( labcnt > nlabels ) labcnt = 1
                end do
            endif
            step = step + 1
        end do
    end function balanced_diverse_labeling


    ! subroutine balanced_diverse_labeling( nptcls, nlabels, ndiverse, configs_diverse, nsamples )
    !     integer,              intent(in) :: nptcls, nlabels, ndiverse
    !     integer, allocatable, intent(out) :: configs_diverse(:,:)
    !     integer, optional,    intent(in) :: nsamples
    !     real, allocatable    :: smat(:,:)
    !     type(ran_tabu)       :: rt
    !     type(shc_cluster)    :: shcc
    !     type(oris)           :: os
    !     logical, parameter   :: DOPRINT=.true., DEBUG=.true.
    !     integer, allocatable :: ptcls_in_which(:), configs(:,:)
    !     integer              :: nnsamples, alloc_stat, isample, i, j, icls, cand, npairs
    !     real                 :: sim, rix, avgrand, minrand, maxrand
    !     nnsamples = 10*ndiverse
    !     if( present(nsamples) ) nnsamples = nsamples
    !     ! build
    !     allocate( configs(nnsamples,nptcls), smat(nnsamples,nnsamples),&
    !               configs_diverse(ndiverse,nptcls), stat=alloc_stat )
    !     call alloc_err( "In: simple_combinatorics :: balanced_diverse_labeling", alloc_stat )
    !     call os%new(nnsamples)
    !     call shcc%new(nnsamples, ndiverse, smat, os, minsim=0.)
    !     if( DEBUG ) print *, 'generating random balanced configurations'
    !     ! generate random balanced configurations
    !     do isample=1,nnsamples
    !         rt = ran_tabu(nptcls)
    !         call rt%balanced(nlabels, configs(isample,:)) 
    !         call rt%kill
    !     end do
    !     if( DEBUG ) print *, 'generating similarity matrix'
    !     ! generate a similarity matrix based on the rand index similarity metric
    !     avgrand = 0.
    !     maxrand = 0.
    !     minrand = 1.
    !     !$omp parallel do schedule(auto) default(shared) private(i,j)
    !     do i=1,nnsamples - 1
    !         do j=i + 1,nnsamples
    !             smat(i,j) = rand_index(configs(i,:), configs(j,:))
    !             smat(j,i) = smat(i,j)
    !         end do
    !     end do
    !     !$omp end parallel do
    !     if( DEBUG ) print *, 'calculating similarity matrix stats'
    !     ! calculate stats
    !     do i=1,nnsamples - 1
    !         do j=i + 1,nnsamples
    !             avgrand   = avgrand + smat(i,j)
    !             if( smat(i,j) > maxrand ) maxrand = smat(i,j)
    !             if( smat(i,j) < minrand ) minrand = smat(i,j)
    !         end do
    !     end do
    !     npairs  = (nnsamples*(nnsamples - 1)) / 2.
    !     avgrand = avgrand/real(npairs)
    !     ! set the diagonal elements to one
    !     forall( i=1:nnsamples ) smat(i,i) = 1.0
    !     if( DEBUG ) print *, 'clustering configurations'
    !     ! cluster the configurations
    !     call shcc%shc(DOPRINT, 'class', sim)
    !     if( DEBUG ) print *, 'selecting candidates'
    !     ! select a candidate at random from each cluster of configurations
    !     do icls=1,ndiverse
    !         ptcls_in_which = os%get_cls_pinds(icls)
    !         configs_diverse(icls,:) = configs(irnd_uni(size(ptcls_in_which)),:)
    !         deallocate(ptcls_in_which)
    !     end do
    !     if( DEBUG ) print *, 'validating'
    !     ! validate the obtained set of solutions
    !     write(*,'(a)') '>>> RANDOM SET STATS'
    !     write(*,'(a,F8.2)') 'AVG:', avgrand
    !     write(*,'(a,F8.2)') 'MAX:', maxrand
    !     write(*,'(a,F8.2)') 'MAX:', minrand
    !     avgrand = 0.
    !     maxrand = 0.
    !     minrand = 1.
    !     do i=1,ndiverse - 1
    !         do j=i + 1,ndiverse
    !             rix = rand_index(configs_diverse(i,:), configs_diverse(j,:))
    !             avgrand = avgrand + rix 
    !             if( rix > maxrand ) maxrand = rix 
    !             if( rix < minrand ) minrand = rix
    !         end do
    !     end do
    !     npairs  = (ndiverse*(ndiverse - 1)) / 2.
    !     avgrand = avgrand/real(npairs)
    !     write(*,'(a)') '>>> DIVERSE SET STATS'
    !     write(*,'(a,F8.2)') 'AVG:', avgrand
    !     write(*,'(a,F8.2)') 'MAX:', maxrand
    !     write(*,'(a,F8.2)') 'MIN:', minrand
    !     ! destruct
    !     call os%kill
    !     call shcc%kill
    !     deallocate( configs, smat )
    ! end subroutine balanced_diverse_labeling



end module simple_combinatorics