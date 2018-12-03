! seeding methods for state labelling
module simple_cluster_seed
include 'simple_lib.f08'
use simple_oris, only: oris
implicit none

public :: gen_labelling
private
#include "simple_local_flags.inc"

integer              :: nptcls      = 0 !< total # of particles
integer              :: nincl_ptcls = 0 !< # of particles with nonzero state
integer              :: nlabels     = 0 !< # of partitions
integer, allocatable :: states(:)       !< states for all particles

contains

    subroutine gen_labelling( os, nclasses, method, power)
        type(oris),        intent(inout) :: os
        integer,           intent(in)    :: nclasses
        character(len=*),  intent(in)    :: method
        real, optional,    intent(in)    :: power
        character(len=STDLEN) :: method_here
        real                  :: power_here
        if(nclasses<=1) THROW_HARD('Inconstsistent number of classes; gen_labelling')
        ! init
        call seed_rnd
        nlabels = nclasses
        nptcls  = os%get_noris()
        if(allocated(states))deallocate(states)
        states      = nint(os%get_all('state'))
        nincl_ptcls = count(states>0)
        ! switch
        if(.not.os%isthere('corr'))then
            ! defaults to random uniform
            method_here = 'uniform'
        else
            method_here = trim(method)
        endif
        power_here = 2.
        if( present(power) )power_here = power
        select case(trim(method))
            case('uniform')
                call draw_uniform
            case('uniform_corr')
                call draw_uniform_corr(os)
            case('squared')
                call draw_squared(os, 2.)
            case('ranked_uniform')
                call draw_ranked_uniform(os)
            case('ranked')
                call draw_ranked(os)
            case('squared_uniform')
                call draw_squared_uniform(os, power_here)
            case DEFAULT
                THROW_HARD('Unsupported method; gen_labelling')
        end select
        ! updates labelling
        call os%set_all('state',real(states))
        ! cleanup
        deallocate(states)
    end subroutine gen_labelling

    !> random uniform labelling
    subroutine draw_uniform
        type(ran_tabu) :: rt
        integer        :: tmp(nincl_ptcls)
        integer        :: iptcl, cnt
        write(logfhandle,'(A)') '>>> GENERATING DIVERSE LABELLING'
        rt = ran_tabu(nincl_ptcls)
        call rt%balanced(nlabels, tmp)
        cnt = 0
        do iptcl=1,nptcls
            if(states(iptcl)>0 .and. states(iptcl)<=nlabels)then
                cnt = cnt + 1
                states(iptcl) = tmp(cnt)
            endif
        enddo
        call rt%kill
    end subroutine draw_uniform

    !>  partitions have a uniform distribution of correlations
    subroutine draw_uniform_corr(os)
        type(oris), intent(inout) :: os
        type(ran_tabu) :: rt
        integer        :: tmp(nlabels), order(nptcls), config(nptcls)
        real           :: corrs(nptcls)
        integer        :: iptcl, s, ind
        write(logfhandle,'(A)') '>>> GENERATING DIVERSE LABELLING WITH RESPECT TO OBJECTIVE FUNCTION'
        config = 0
        corrs = os%get_all('corr')
        tmp   = (/(s,s=1,nlabels)/)
        where( states<=0 ) corrs = -1.
        order = (/(iptcl,iptcl=1,nptcls)/)
        call hpsort( corrs, order )
        call reverse(order)
        call reverse(corrs)
        rt = ran_tabu(nlabels)
        do iptcl = 1, nincl_ptcls+nlabels, nlabels
            call rt%reset
            call rt%shuffle(tmp)
            do s = 1, nlabels
                ind = iptcl + s - 1
                if(ind > nptcls)exit
                config(order(ind)) = tmp(s)
            enddo
        enddo
        ! update states
        where((states > 0) .and. (states <= nlabels)) states = config
        ! cleanup
        call rt%kill
    end subroutine draw_uniform_corr

    !>  partitions have a uniform distribution of correlations
    subroutine draw_ranked(os)
        type(oris), intent(inout) :: os
        type(ran_tabu) :: rt
        integer        :: tmp(nlabels), order(nptcls), config(nptcls), pops(nlabels)
        real           :: corrs(nptcls)
        integer        :: iptcl, s, ind
        write(logfhandle,'(A)') '>>> GENERATING ranked'
        ! even partitions
        pops = floor(real(nincl_ptcls)/real(nlabels))
        do s=1,nlabels
            if(sum(pops)==nincl_ptcls)exit
            pops(s) = pops(s)+1
        enddo
        config = 0
        corrs = os%get_all('corr')
        where( states<=0 ) corrs = -1.
        order = (/(iptcl,iptcl=1,nptcls)/)
        call hpsort( corrs, order )
        call reverse(order)
        iptcl = 0
        do s=1,nlabels
            do ind=1,pops(s)
                iptcl = iptcl+1
                config(order(iptcl)) = s
            enddo
        enddo
        ! update states
        where((states > 0) .and. (states <= nlabels)) states = config
        ! cleanup
    end subroutine draw_ranked

    !>  first partition squared, all others uniform, simultaneous sampling
    subroutine draw_squared_uniform(os, power)
        type(oris), intent(inout) :: os
        real,       intent(in)    :: power
        integer        :: tmp(nlabels), order(nptcls), config(nptcls)
        real           :: corrs(nptcls), rnincl
        integer        :: iptcl, s, ind, i, n95
        logical        :: mask(nincl_ptcls)
        write(logfhandle,'(A)') '>>> MIXED SQUARED & UNIFORM SAMPLING'
        config = 0
        rnincl = real(nincl_ptcls)
        mask   = .true.
        corrs  = os%get_all('corr')
        order = (/(iptcl,iptcl=1,nptcls)/)
        where( states<=0 ) corrs = -1.
        call hpsort( corrs, order )
        call reverse(order)
        ! 95% of the data is sampled following: first parttion follos squared distribution sampling,
        ! all others are uniformly sampled
        n95   = nint(0.95*real(nincl_ptcls)) + nlabels
        iptcl = 0
        do i=1,n95,nlabels
            if(i>nincl_ptcls)exit
            ind = ceiling((ran3()**power)*rnincl)
            do while(.not.mask(ind))
                ind = ceiling((ran3()**power)*rnincl)
            enddo
            config(order(ind)) = 1
            mask(ind)          = .false.
            iptcl              = i
            do s=2,nlabels
                iptcl = iptcl+1
                if(iptcl>nincl_ptcls)exit
                ind = ceiling(ran3()*rnincl)
                do while(.not.mask(ind))
                    ind = ceiling(ran3()*rnincl)
                enddo
                config(order(ind)) = s
                mask(ind)          = .false.
            enddo
        enddo
        ! the remaining 5% are drawn uniformly to avoid overhead in big datasets
        s = 0
        do iptcl=1,nptcls
            if( states(iptcl) == 0 )cycle
            if( config(iptcl) /= 0 )cycle
            s = s + 1
            if( s > nlabels ) s = 1
            config(iptcl) = s
        enddo
        ! update states
        where((states > 0) .and. (states <= nlabels)) states = config
        ! cleanup
    end subroutine draw_squared_uniform

    ! Loose adaptation of k++ seeding procedure
    subroutine draw_squared(os, power)
        real,       intent(in)    :: power
        type(oris), intent(inout) :: os
        real     :: dists(nptcls), dists_part(nptcls)
        integer  :: inds(nptcls), pops(nlabels), config(nptcls)
        logical  :: mask(nptcls)
        real     :: areal
        integer  :: iptcl, s, n_drawn, ind, n_avail
        write(logfhandle,'(A)') '>>> GENERATING SKEWED LABELLING WITH RESPECT TO OBJECTIVE FUNCTION'
        mask   = (states > 0) .and. (states <= nlabels)
        dists  = 1. - os%get_all('corr')
        config = 0
        ! even partitions
        pops = floor(real(nincl_ptcls)/real(nlabels))
        do s=1,nlabels
            if(sum(pops)==nincl_ptcls)exit
            pops(s) = pops(s)+1
        enddo
        ! draw following squared distribution
        do s=1,nlabels-1
            n_avail    = count(mask)
            inds       = (/(iptcl, iptcl=1,nptcls)/)
            dists_part = huge(areal)
            if( s==1 )then
                ! first partition: sample lowest distances
                where(mask) dists_part = dists
            else
                ! other partitions: sample highest distances
                where(mask) dists_part = 1. - dists
            endif
            call hpsort(dists_part,inds)
            n_drawn = 0
            do while(n_drawn < pops(s))
                ind   = ceiling((ran3()**power)*real(n_avail))
                iptcl = inds(ind)
                if(mask(iptcl))then
                    mask(iptcl)   = .false.
                    n_drawn       = n_drawn + 1
                    config(iptcl) = s
                endif
            enddo
        enddo
        ! last partition: leftovers
        where(mask)config = nlabels
        ! updates states
        where((states > 0) .and. (states <= nlabels)) states = config
    end subroutine draw_squared

    ! Bottom half randomization followed by sorted assignment for the top half
    subroutine draw_ranked_uniform(os)
        type(oris), intent(inout) :: os
        real     :: dists(nptcls), dists_part(nptcls)
        integer  :: inds(nptcls), pops(nlabels), config(nptcls)
        logical  :: mask(nptcls)
        real     :: areal
        integer  :: i,iptcl, cnt, s, n_avail, half_pop
        write(logfhandle,'(A)') '>>> GENERATING MIXED UNIFORM & RANKED LABELLING'
        mask   = (states > 0) .and. (states <= nlabels)
        dists  = 1. - os%get_all('corr')
        config = 0
        ! assign uniformly the bottom half worst scoring particles
        ! this is only for consistency with the first iteration of extremal optimization
        half_pop   = ceiling(real(nincl_ptcls)/2.)
        inds       = (/(iptcl, iptcl=1,nptcls)/)
        dists_part = -huge(areal)
        where(mask) dists_part = dists
        call hpsort(dists_part,inds) ! highest distances last
        s = irnd_uni(nlabels)
        do i=nptcls,half_pop,-1
            iptcl         = inds(i)
            config(iptcl) = s
            mask(iptcl)   = .false.
            s = s + 1
            if(s > nlabels) s = 1
        enddo
        ! even partitions for the top half
        n_avail = count(mask)
        pops    = floor(real(n_avail)/real(nlabels))
        do s=1,nlabels
            if(sum(pops)==n_avail)exit
            pops(s) = pops(s)+1
        enddo
        ! deterministic ranked assignement for the top half
        n_avail    = count(mask)
        inds       = (/(iptcl, iptcl=1,nptcls)/)
        dists_part = huge(areal)
        where(mask) dists_part = dists
        call hpsort(dists_part,inds)
        s   = 1
        cnt = 0
        do i=1,n_avail
            cnt = cnt + 1
            if(cnt > pops(s))then
                s   = s + 1
                cnt = 0
            endif
            iptcl         = inds(i)
            config(iptcl) = s
        enddo
        ! updates states
        where((states > 0) .and. (states <= nlabels)) states = config
    end subroutine draw_ranked_uniform

end module simple_cluster_seed
