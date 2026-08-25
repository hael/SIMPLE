!@descr: seeding methods for state labelling
module simple_cluster_seed
use simple_core_module_api

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
        if(nclasses<=1) THROW_HARD('Inconsistent number of classes; gen_labelling')
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
        select case(trim(method_here))
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
            case('squared_uniform_proj')
                if( .not.os%isthere('proj') )then
                    write(logfhandle,'(A)') '>>> WARNING: no projection direction information; falling back to uniform'
                    call draw_uniform
                else
                    call draw_squared_uniform_projdir(os, power_here)
                endif
            case DEFAULT
                THROW_HARD('Unsupported method; gen_labelling')
        end select
        ! updates labelling
        call os%set_all('state', states)
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
        integer        :: order(nptcls), config(nptcls), pops(nlabels)
        real           :: corrs(nptcls)
        integer        :: iptcl, s, ind
        write(logfhandle,'(A)') '>>> GENERATING RANKED'
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
        integer        :: order(nptcls), config(nptcls), pops(nlabels)
        real           :: corrs(nptcls), rnincl_1
        integer        :: iptcl, s, ind, i, n99, pop
        logical        :: mask(nincl_ptcls)
        write(logfhandle,'(A)') '>>> MIXED SQUARED & UNIFORM SAMPLING'
        config   = 0
        rnincl_1 = real(nincl_ptcls-1)
        mask     = .true.
        corrs    = os%get_all('corr')
        order    = (/(iptcl,iptcl=1,nptcls)/)
        where( states<=0 ) corrs = -1.
        call hpsort( corrs, order )
        call reverse(order)
        ! even populations
        pops = floor(real(nincl_ptcls)/real(nlabels))
        do s=1,nlabels
            if(sum(pops)==nincl_ptcls)exit
            pops(s) = pops(s)+1
        enddo
        ! 99% of the data is sampled following: first partition follows squared distribution sampling,
        ! all others are uniformly sampled
        n99   = nint(0.99*real(nincl_ptcls)) + nlabels
        iptcl = 0
        do i=1,n99,nlabels
            if(i>nincl_ptcls)exit
            ind = ceiling(ran3()**power * rnincl_1) + 1
            do while(.not.mask(ind))
                ind = ceiling(ran3()**power * rnincl_1) + 1
            enddo
            config(order(ind)) = 1
            mask(ind)          = .false.
            iptcl              = i
            do s=2,nlabels
                iptcl = iptcl+1
                if(iptcl>nincl_ptcls)exit
                ind = ceiling(ran3() * rnincl_1) + 1
                do while(.not.mask(ind))
                    ind = ceiling(ran3() * rnincl_1) + 1
                enddo
                config(order(ind)) = s
                mask(ind)          = .false.
            enddo
        enddo
        ! remaining 1%:
        ! greedy assignment for first partition
        pop = count(config==1)
        do ind=1,nincl_ptcls
            if( mask(ind) )then
                pop = pop+1
                if( pop > pops(1) )exit
                iptcl = order(ind)
                config(order(ind)) = 1
                mask(ind)          = .false.
            endif
        enddo
        ! uniform assignement for all others
        s = 1
        do iptcl=1,nptcls
            if( states(iptcl) == 0 )cycle
            if( config(iptcl) /= 0 )cycle
            s = s + 1
            if( s > nlabels ) s = 2
            config(iptcl) = s
        enddo
        ! update states
        where((states > 0) .and. (states <= nlabels)) states = config
        ! cleanup
    end subroutine draw_squared_uniform

    !>  first partition squared, all others uniform, simultaneous sampling by projection direction
    subroutine draw_squared_uniform_projdir(os, power)
        type(oris), intent(inout) :: os
        real,       intent(in)    :: power
        integer, allocatable      :: order(:), config(:), projs(:), fen(:)
        integer, allocatable      :: bucket(:), bcount(:), boffset(:), cursor(:)
        real,    allocatable      :: corrs(:), projcorrs(:)
        integer        :: iptcl, s, i, n95, ninproj, iproj, nprojs, ofs, navail
        write(logfhandle,'(A)') '>>> MIXED SQUARED & UNIFORM PROJECTION DIRECTION SAMPLING'
        allocate(order(nptcls), config(nptcls), corrs(nptcls), projs(nptcls), projcorrs(nptcls), fen(nptcls))
        config = 0
        corrs  = os%get_all('corr')
        projs  = int(os%get_all('proj'))
        ! particles not eligible for relabelling are excluded from every proj bucket and sort
        where( .not.((states > 0) .and. (states <= nlabels)) )
            corrs = -1.
            projs = 0
        end where
        nprojs = maxval(projs)
        ! bucket particles by proj id once (counting sort), replacing the per-projection
        ! pack() that made the original loop O(nptcls * nprojs)
        allocate(bucket(nptcls), bcount(nprojs), boffset(0:nprojs), cursor(nprojs))
        bcount = 0
        do iptcl = 1,nptcls
            if( projs(iptcl) > 0 ) bcount(projs(iptcl)) = bcount(projs(iptcl)) + 1
        enddo
        boffset(0) = 0
        do iproj = 1,nprojs
            boffset(iproj) = boffset(iproj-1) + bcount(iproj)
        enddo
        cursor = boffset(0:nprojs-1)
        do iptcl = 1,nptcls
            iproj = projs(iptcl)
            if( iproj == 0 ) cycle
            cursor(iproj) = cursor(iproj) + 1
            bucket(cursor(iproj)) = iptcl
        enddo
        do iproj = 1,nprojs
            ninproj = bcount(iproj)
            if( ninproj == 0 ) cycle
            ofs = boffset(iproj-1)
            if( ninproj <= 2*nlabels )then
                ! random uniform
                s = irnd_uni(nlabels)
                do i = 1, ninproj
                    iptcl = bucket(ofs+i)
                    s = s + 1
                    if( s > nlabels ) s = 1
                    config(iptcl) = s
                enddo
            else
                ! order/projcorrs are preallocated once above and reused per iproj via slicing,
                ! avoiding the per-iteration allocate/deallocate the original incurred
                order(1:ninproj)     = bucket(ofs+1:ofs+ninproj)
                projcorrs(1:ninproj) = corrs(order(1:ninproj))
                call hpsort(projcorrs(1:ninproj), order(1:ninproj))
                call reverse(order(1:ninproj))
                ! Fenwick (BIT) tree over ranks 1..ninproj tracks which ranks are still
                ! available, so draw() can pick-and-remove the k-th available rank in
                ! O(log ninproj) instead of the original rejection loop, whose cost blows
                ! up (up to O(ninproj) per draw) once most low ranks have been consumed
                fen(1:ninproj) = 1
                do i = 1,ninproj
                    if( i + iand(i,-i) <= ninproj ) fen(i+iand(i,-i)) = fen(i+iand(i,-i)) + fen(i)
                enddo
                navail = ninproj
                ! First 95%, stochastic rank-powered for first partition
                n95    = min(ninproj-1, nint(0.95*real(ninproj)) + nlabels)
                iptcl  = 0
                do i = 1, n95, nlabels
                    if( i > ninproj )exit
                    call draw(power, 1)
                    iptcl = i
                    do s = 2, nlabels
                        iptcl = iptcl+1
                        if( iptcl > ninproj )exit
                        call draw(1.0, s)
                    enddo
                enddo
                ! Leftovers: random uniform
                s = irnd_uni(nlabels)
                do i = 1, ninproj
                    iptcl = bucket(ofs+i)
                    if( config(iptcl) == 0 )then
                        s = s + 1
                        if( s > nlabels ) s = 1
                        config(iptcl) = s
                    endif
                enddo
            endif
        enddo
        where((states > 0) .and. (states <= nlabels)) states = config
        ! cleanup
        deallocate(config,corrs,projs,order,projcorrs,fen,bucket,bcount,boffset,cursor)
    contains

        subroutine draw( p, s )
            real,    intent(in) :: p
            integer, intent(in) :: s
            integer :: k, rank
            k    = ceiling(ran3()**p * real(navail-1)) + 1
            rank = bit_select(k)
            config(order(rank)) = s
            call bit_remove(rank)
            navail = navail - 1
        end subroutine draw

        ! finds the position of the k-th rank still flagged available in fen(1:ninproj)
        integer function bit_select(k)
            integer, intent(in) :: k
            integer :: remaining, pw, idx
            remaining = k
            idx       = 0
            pw        = 1
            do while( pw*2 <= ninproj )
                pw = pw*2
            enddo
            do while( pw > 0 )
                if( idx+pw <= ninproj )then
                    if( fen(idx+pw) < remaining )then
                        idx       = idx + pw
                        remaining = remaining - fen(idx)
                    endif
                endif
                pw = pw/2
            enddo
            bit_select = idx + 1
        end function bit_select

        ! marks rank as no longer available in the Fenwick tree
        subroutine bit_remove(rank)
            integer, intent(in) :: rank
            integer :: pos
            pos = rank
            do while( pos <= ninproj )
                fen(pos) = fen(pos) - 1
                pos = pos + iand(pos,-pos)
            enddo
        end subroutine bit_remove

    end subroutine draw_squared_uniform_projdir

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
                ind   = ceiling((ran3()**power)*real(n_avail - 1)) + 1
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
