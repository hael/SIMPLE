!@descr: sampling and updatecnt related routines for oris
submodule (simple_oris) simple_oris_sampling
use simple_ori_api
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine select_particles_set( self, nptcls, inds )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nptcls
        integer, allocatable, intent(inout) :: inds(:)
        integer        :: i,n
        if( allocated(inds) ) deallocate(inds)
        inds = (/(i,i=1,self%n)/)
        inds = pack(inds, mask=self%get_all('state') > 0.5)
        n    = size(inds)
        if( n < nptcls )then
            ! if less than desired particles select all
        else
            call partial_shuffle(inds, nptcls)
            inds = inds(1:nptcls)
            call hpsort(inds)
        endif
    end subroutine select_particles_set

    ! Active rows of fromto to reconstruct: those with updatecnt > 0 when any
    ! active row of the whole project has been updated, else every active row.
    ! The coverage decision is global, so that all partitions of a distributed
    ! reconstruction select from the same population.
    module subroutine sample4rec( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        integer, allocatable :: states(:), updatecnts(:)
        integer :: i, cnt, nptcls
        logical :: l_any_updated
        l_any_updated = any_active_updated(self)
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), updatecnts(nptcls), inds(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_int('updatecnt')
            inds(cnt)       = i
        end do
        if( l_any_updated )then
            nsamples = count(states > 0 .and. updatecnts > 0)
            inds     = pack(inds, mask=states > 0 .and. updatecnts > 0)
        else
            nsamples = count(states > 0)
            inds     = pack(inds, mask=states > 0)
        endif
    end subroutine sample4rec

    ! Per-state populations of the rows sample4rec reconstructs over the whole
    ! project: the population a full reconstruction (a trailing-chain seed) represents
    module subroutine get_state_rec_pops( self, nstates, pops )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: nstates
        integer, allocatable, intent(inout) :: pops(:)
        integer :: i, s
        logical :: l_any_updated
        if( allocated(pops) ) deallocate(pops)
        allocate(pops(nstates), source=0)
        l_any_updated = any_active_updated(self)
        do i = 1, self%n
            s = self%o(i)%get_state()
            if( s < 1 .or. s > nstates ) cycle
            if( l_any_updated .and. self%o(i)%get_int('updatecnt') <= 0 ) cycle
            pops(s) = pops(s) + 1
        end do
    end subroutine get_state_rec_pops

    !> Population-rule weights of one group: new = s*current + w*previous, f = n/N, u = ufrac (in [0,1]) or f,
    !! s = u/f, w = (1-u)*N/M (0 if M = 0), mnew = s*n + w*M (= N when M > 0). n = 0 keeps the previous sums
    !! (s = 0, w = N/M); N = 0 gives zeros. nrep = N active updated rows, nsmp = n sampled, mrep = M stored mass.
    elemental module subroutine population_blend_weights( nrep, nsmp, mrep, s, w, mnew, ufrac )
        integer,        intent(in)  :: nrep, nsmp
        real,           intent(in)  :: mrep
        real,           intent(out) :: s, w, mnew
        real, optional, intent(in)  :: ufrac
        real :: f, u
        s    = 0.
        w    = 0.
        mnew = 0.
        if( nrep <= 0 ) return
        if( nsmp <= 0 )then
            ! nothing sampled in the group: keep the previous sums at mass N
            if( mrep > 0. ) w = real(nrep) / mrep
        else
            f = real(min(nsmp, nrep)) / real(nrep)
            u = f
            if( present(ufrac) ) u = max(0., min(1., ufrac))
            s = u / f
            if( mrep > 0. ) w = (1. - u) * real(nrep) / mrep
        endif
        mnew = s * real(min(nsmp, nrep)) + w * mrep
    end subroutine population_blend_weights

    ! whether any active row of the project has been updated (sample4rec's coverage decision)
    logical function any_active_updated( self )
        class(oris), intent(in) :: self
        integer :: i
        any_active_updated = .false.
        do i = 1, self%n
            if( self%o(i)%get_state() > 0 .and. self%o(i)%get_int('updatecnt') > 0 )then
                any_active_updated = .true.
                return
            endif
        end do
    end function any_active_updated

    module subroutine sample4update_all( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: states(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt      = 0
        nsamples = 0
        do i = fromto(1), fromto(2)
            cnt         = cnt + 1
            states(cnt) = self%o(i)%get_state()
            inds(cnt)   = i
            if( states(cnt) > 0 ) nsamples = nsamples + 1
        end do
        inds = pack(inds, mask=states > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_all

    module subroutine sample4update_rnd( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: states(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt    = 0
        nptcls = 0
        do i = fromto(1), fromto(2)
            cnt         = cnt + 1
            states(cnt) = self%o(i)%get_state()
            inds(cnt)   = i
            if( states(cnt) > 0 ) nptcls = nptcls + 1
        end do
        if( nptcls == 0 ) THROW_HARD('no active particles to sample')
        inds     = pack(inds, mask=states > 0)
        nsamples = min(nptcls, max(1, nint(update_frac * real(nptcls))))
        call partial_shuffle(inds, nsamples)
        inds = inds(1:nsamples)
        call hpsort(inds)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_rnd

    module subroutine sample4update_cnt( self, fromto, update_frac, nsamples, inds, incr_sampled, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        logical, optional,    intent(in)    :: allow_empty
        integer, allocatable :: states(:), updatecnts(:), candidates(:), selected(:)
        integer :: i, cnt, nptcls, ucnt
        integer :: ncandidates, nfill, nselected
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), updatecnts(nptcls), inds(nptcls), source=0)
        cnt    = 0
        nptcls = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_int('updatecnt')
            inds(cnt)       = i
            if( states(cnt) > 0 ) nptcls = nptcls + 1
        end do
        if( nptcls == 0 )then
            ! a distributed partition may hold no active particle
            if( .not. empty_allowed(allow_empty) ) THROW_HARD('no active particles to sample')
            nsamples = 0
            deallocate(inds)
            allocate(inds(0))
            return
        endif
        inds       = pack(inds,       mask=states > 0)
        updatecnts = pack(updatecnts, mask=states > 0)
        deallocate(states)
        nsamples   = min(nptcls, max(1, nint(update_frac * real(nptcls))))
        allocate(selected(nsamples), source=0)
        nselected = 0
        ! Exhaust lower update-count tiers first. At the cutoff tier, draw
        ! uniformly without replacement.
        ucnt = minval(updatecnts)
        do
            candidates  = pack(inds, mask=updatecnts == ucnt)
            ncandidates = size(candidates)
            nfill = min(nsamples - nselected, ncandidates)
            if( nfill < ncandidates ) call partial_shuffle(candidates, nfill)
            selected(nselected + 1:nselected + nfill) = candidates(:nfill)
            nselected = nselected + nfill
            if( nselected == nsamples ) exit
            ucnt = minval(updatecnts, mask=updatecnts > ucnt)
        end do
        if( nselected /= nsamples ) THROW_HARD('insufficient update-count sampling candidates')
        call move_alloc(selected, inds)
        call hpsort(inds)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_cnt

    !> Nested equal quota: nint(update_frac * active) is shared equally over the groups (capped at their
    !! populations), each group's share equally over its units (the remainder to the least-updated units),
    !! and a unit draws its quota lowest updatecnt first. Units with group 0 are their own group.
    module subroutine sample4update_class( self, clssmp, fromto, update_frac, nsamples, inds, incr_sampled, l_greedy, &
        &frac_best, sampled_only, allow_empty )
        class(oris),          intent(inout) :: self
        type(class_sample),   intent(inout) :: clssmp(:)
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled, l_greedy
        real,    optional,    intent(in)    :: frac_best
        logical, optional,    intent(in)    :: sampled_only, allow_empty
        integer, allocatable :: states(:), eligible(:), updatecnts(:), sampleinds(:), inds_pool(:), inds_fill(:)
        real,    allocatable :: rstates(:)
        integer :: i, j, cnt, nptcls, nsamples_class, states_bal(self%n)
        integer :: nbest, nfill, ucnt, ucnt_min, ucnt_max
        logical :: l_sampled_only
        l_sampled_only = .false.
        if( present(sampled_only) ) l_sampled_only = sampled_only
        rstates        = self%get_all('state')
        nsamples_class = nint(update_frac * real(count(rstates > 0.5)))
        deallocate(rstates)
        call alloc_unit_quotas(self, clssmp, nsamples_class)
        if( l_sampled_only )then
            ! The equal group allocation may overshoot by up to one particle
            ! per group. Cohort sampling is also the exact-K split
            ! reconstruction contract, so trim that overshoot deterministically.
            do i = size(clssmp), 1, -1
                if( sum(clssmp(:)%nsample) == nsamples_class ) exit
                if( clssmp(i)%nsample > 0 ) clssmp(i)%nsample = clssmp(i)%nsample - 1
            enddo
        endif
        states_bal = 0
        do i = 1, size(clssmp)
            if( clssmp(i)%nsample < 1 ) cycle
            if( present(frac_best) )then
                nbest = max(clssmp(i)%nsample, nint(frac_best * real(clssmp(i)%pop)))
                nbest = min(nbest, clssmp(i)%pop)
                eligible = clssmp(i)%pinds(:nbest)
            else if( l_greedy .and. .not. l_sampled_only )then
                do j = 1, clssmp(i)%nsample
                    states_bal(clssmp(i)%pinds(j)) = 1
                end do
                cycle
            else
                eligible = clssmp(i)%pinds
            endif
            if( l_sampled_only )then
                allocate(sampleinds(size(eligible)), source=0)
                do j = 1, size(eligible)
                    sampleinds(j) = self%o(eligible(j))%get_sampled()
                enddo
                eligible = pack(eligible, mask=sampleinds > 0)
                deallocate(sampleinds)
                if( size(eligible) < clssmp(i)%nsample )then
                    THROW_HARD('insufficient previously sampled class-balanced update candidates')
                endif
            endif
            allocate(updatecnts(size(eligible)), source=0)
            do j = 1, size(eligible)
                updatecnts(j) = self%o(eligible(j))%get_int('updatecnt')
            end do
            ucnt_min = minval(updatecnts)
            ucnt_max = maxval(updatecnts)
            allocate(inds_pool(0), source=0)
            do ucnt = ucnt_min, ucnt_max
                if( size(inds_pool) >= clssmp(i)%nsample ) exit
                inds_fill = pack(eligible, mask=updatecnts == ucnt)
                if( size(inds_fill) == 0 ) cycle
                call shuffle(inds_fill)
                nfill = min(clssmp(i)%nsample - size(inds_pool), size(inds_fill))
                inds_pool = [inds_pool, inds_fill(1:nfill)]
            end do
            if( size(inds_pool) < clssmp(i)%nsample ) THROW_HARD('insufficient class-balanced update candidates')
            do j = 1, clssmp(i)%nsample
                states_bal(inds_pool(j)) = 1
            end do
            if( allocated(eligible)  ) deallocate(eligible)
            if( allocated(updatecnts)) deallocate(updatecnts)
            if( allocated(inds_pool) ) deallocate(inds_pool)
            if( allocated(inds_fill) ) deallocate(inds_fill)
        end do
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(states(nptcls), inds(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt          = cnt + 1
            states(cnt)  = states_bal(i)
            inds(cnt)    = i
        end do
        nsamples = count(states > 0)
        ! the class-balanced sample is drawn over the whole project, so a distributed
        ! partition may receive none of it
        if( nsamples == 0 .and. .not. empty_allowed(allow_empty) ) THROW_HARD('no active particles to sample')
        inds     = pack(inds, mask=states > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_class

    !> Integer unit quotas (clssmp%nsample) of the nested equal quota. Inside a group the remainder of
    !! the equal split goes to the units with the lowest mean updatecnt (ties: larger population, then
    !! lower index), so it is the same in every partition and rotates over the iterations.
    subroutine alloc_unit_quotas( self, clssmp, ntarget )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: ntarget
        integer, allocatable :: gquota(:), members(:), open_units(:)
        real,    allocatable :: ucnt_mean(:)
        logical, allocatable :: l_taken(:)
        integer :: gind(size(clssmp)), ngroups, nunits, g, i, k, remaining, nopen, ibest
        nunits = size(clssmp)
        clssmp(:)%nsample = 0
        if( nunits < 1 ) return
        call unit_groups(clssmp, gind, ngroups)
        call group_quotas(clssmp, gind, ngroups, ntarget, gquota)
        allocate(ucnt_mean(nunits), source=0.)
        do i = 1, nunits
            if( clssmp(i)%pop < 1 ) cycle
            do k = 1, clssmp(i)%pop
                ucnt_mean(i) = ucnt_mean(i) + real(self%o(clssmp(i)%pinds(k))%get_int('updatecnt'))
            end do
            ucnt_mean(i) = ucnt_mean(i) / real(clssmp(i)%pop)
        end do
        do g = 1, ngroups
            members   = pack([(i, i=1,nunits)], mask=gind == g)
            remaining = gquota(g)
            do
                open_units = pack(members, mask=clssmp(members)%nsample < clssmp(members)%pop)
                nopen      = size(open_units)
                if( remaining <= 0 .or. nopen == 0 ) exit
                if( remaining >= nopen )then
                    clssmp(open_units)%nsample = clssmp(open_units)%nsample + 1
                    remaining = remaining - nopen
                else
                    allocate(l_taken(nopen), source=.false.)
                    do k = 1, remaining
                        ibest = 0
                        do i = 1, nopen
                            if( l_taken(i) ) cycle
                            if( ibest == 0 )then
                                ibest = i
                            else if( takes_remainder(open_units(i), open_units(ibest)) )then
                                ibest = i
                            endif
                        end do
                        l_taken(ibest) = .true.
                        clssmp(open_units(ibest))%nsample = clssmp(open_units(ibest))%nsample + 1
                    end do
                    deallocate(l_taken)
                    remaining = 0
                endif
            end do
        end do

    contains

        logical function takes_remainder( a, b )
            integer, intent(in) :: a, b
            if( ucnt_mean(a) /= ucnt_mean(b) )then
                takes_remainder = ucnt_mean(a) < ucnt_mean(b)
            else if( clssmp(a)%pop /= clssmp(b)%pop )then
                takes_remainder = clssmp(a)%pop > clssmp(b)%pop
            else
                takes_remainder = a < b
            endif
        end function takes_remainder

    end subroutine alloc_unit_quotas

    !> group index of every unit; a unit with group 0 is a group of its own
    subroutine unit_groups( clssmp, gind, ngroups )
        type(class_sample), intent(in)  :: clssmp(:)
        integer,            intent(out) :: gind(size(clssmp)), ngroups
        integer, allocatable :: labels(:)
        integer :: i, k
        allocate(labels(0))
        do i = 1, size(clssmp)
            k = 0
            if( clssmp(i)%group > 0 ) k = findloc(labels, clssmp(i)%group, dim=1)
            if( k == 0 )then
                if( clssmp(i)%group > 0 )then
                    labels = [labels, clssmp(i)%group]
                else
                    labels = [labels, -i]
                endif
                k = size(labels)
            endif
            gind(i) = k
        end do
        ngroups = size(labels)
    end subroutine unit_groups

    !> equal group quotas capped at the group populations; every open group gets one more particle
    !! per round, so the total may exceed ntarget by fewer than the number of groups
    subroutine group_quotas( clssmp, gind, ngroups, ntarget, gquota )
        type(class_sample),   intent(in)    :: clssmp(:)
        integer,              intent(in)    :: gind(:), ngroups, ntarget
        integer, allocatable, intent(inout) :: gquota(:)
        integer :: gpop(ngroups), g
        do g = 1, ngroups
            gpop(g) = sum(clssmp(:)%pop, mask=gind == g)
        end do
        if( allocated(gquota) ) deallocate(gquota)
        allocate(gquota(ngroups), source=0)
        do while( sum(gquota) < ntarget .and. any(gquota < gpop) )
            where( gquota < gpop ) gquota = gquota + 1
        end do
    end subroutine group_quotas

    module subroutine class_sample_quotas( clssmp, ntarget, quotas )
        type(class_sample), intent(in)  :: clssmp(:)
        integer,            intent(in)  :: ntarget
        real,               intent(out) :: quotas(size(clssmp))
        integer, allocatable :: gquota(:), members(:), open_units(:), capped(:)
        integer :: gind(size(clssmp)), ngroups, nunits, g, i
        real    :: remaining, share
        nunits = size(clssmp)
        quotas = 0.
        if( nunits < 1 ) return
        call unit_groups(clssmp, gind, ngroups)
        call group_quotas(clssmp, gind, ngroups, ntarget, gquota)
        do g = 1, ngroups
            members    = pack([(i, i=1,nunits)], mask=gind == g .and. clssmp(:)%pop > 0)
            open_units = members
            remaining  = real(gquota(g))
            do
                if( size(open_units) == 0 .or. remaining <= 0. ) exit
                share  = remaining / real(size(open_units))
                capped = pack(open_units, mask=real(clssmp(open_units)%pop) <= share)
                if( size(capped) == 0 )then
                    quotas(open_units) = share
                    exit
                endif
                quotas(capped) = real(clssmp(capped)%pop)
                remaining      = remaining - real(sum(clssmp(capped)%pop))
                open_units     = pack(open_units, mask=real(clssmp(open_units)%pop) > share)
            end do
        end do
    end subroutine class_sample_quotas

    module function class_sample_sweep( clssmp, ntarget ) result( sweep )
        type(class_sample), intent(in) :: clssmp(:)
        integer,            intent(in) :: ntarget
        integer :: sweep
        real    :: quotas(size(clssmp))
        integer :: i
        call class_sample_quotas(clssmp, ntarget, quotas)
        sweep = 1
        do i = 1, size(clssmp)
            if( clssmp(i)%pop < 1 ) cycle
            if( quotas(i) <= 0. )then
                sweep = huge(sweep)
                return
            endif
            sweep = max(sweep, ceiling(real(clssmp(i)%pop) / quotas(i) - 1.e-4))
        end do
    end function class_sample_sweep

    !> The particles of fromto the previous sampling selected. allow_empty
    !! accepts a range without any (a distributed partition of state-0 rows
    !! only, e.g. the frozen rows of solve3D_addon): nsamples is 0 and inds
    !! empty, for a caller that emits empty partition outputs
    module subroutine sample4update_reprod( self, fromto, nsamples, inds, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical, optional,    intent(in)    :: allow_empty
        integer, allocatable :: sampled(:)
        integer :: i, cnt, nptcls, sample_ind
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), sampled(nptcls), source=0)
        cnt        = 0
        sample_ind = self%get_sample_ind(.false.)
        if( sample_ind == 0 ) THROW_HARD('requires previous sampling')
        do i = fromto(1), fromto(2)
            cnt          = cnt + 1
            inds(cnt)    = i
            sampled(cnt) = self%o(i)%get_sampled()
        end do
        nsamples = count(sampled == sample_ind)
        if( nsamples == 0 .and. .not. empty_allowed(allow_empty) ) THROW_HARD('no particles sampled in previous sampling')
        inds     = pack(inds, mask=sampled == sample_ind)
    end subroutine sample4update_reprod

    !> Cohort rescoring: the rows of fromto in the latest sampling round, as sample4update_reprod
    !! returns them, with sampled and updatecnt advanced once, so the next call finds the same set
    module subroutine sample4update_rescore( self, fromto, nsamples, inds, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical, optional,    intent(in)    :: allow_empty
        call self%sample4update_reprod(fromto, nsamples, inds, allow_empty)
        if( nsamples > 0 ) call self%incr_sampled_updatecnt(inds, .true.)
    end subroutine sample4update_rescore

    module subroutine sample4update_updated( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: updatecnts(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), updatecnts(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            inds(cnt)       = i
            updatecnts(cnt) = self%o(i)%get_updatecnt()
        end do
        if( .not. any(updatecnts > 0) ) THROW_HARD('requires previous update')
        nsamples = count(updatecnts > 0)
        if( nsamples == 0 ) THROW_HARD('no particles updated in previous update')
        inds     = pack(inds, mask=updatecnts > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_updated

    module subroutine sample4update_fillin( self, fromto, update_frac, nsamples, inds, incr_sampled, allow_empty )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        logical, optional,    intent(in)    :: allow_empty
        integer, allocatable :: updatecnts(:), states(:), updatecnts_active(:), inds_pool(:), inds_fill(:)
        integer :: i, cnt, nptcls, nptcls_active, maxucnt, nfill, ucnt
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), states(nptcls), updatecnts(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            inds(cnt)       = i
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_updatecnt()
        end do
        nptcls_active = count(states > 0)
        if( nptcls_active == 0 )then
            ! a distributed partition may hold no active particle
            if( .not. empty_allowed(allow_empty) ) THROW_HARD('no active particles to sample for fill-in')
            nsamples = 0
            deallocate(inds)
            allocate(inds(0))
            return
        endif
        updatecnts_active = pack(updatecnts, mask=states > 0)
        maxucnt           = maxval(updatecnts_active)
        nsamples          = min(nptcls_active, max(1, nint(update_frac * real(nptcls_active))))
        allocate(inds_pool(0), source=0)
        ! Prefer never-updated particles first (updatecnt==0), then 1, 2, ...
        do ucnt = 0, maxucnt
            if( size(inds_pool) >= nsamples ) exit
            inds_fill = pack(inds, mask=states > 0 .and. updatecnts == ucnt)
            if( size(inds_fill) == 0 ) cycle
            call shuffle(inds_fill)
            nfill = min(nsamples - size(inds_pool), size(inds_fill))
            inds_pool = [inds_pool, inds_fill(1:nfill)]
        enddo
        if( size(inds_pool) < nsamples )then
            THROW_HARD('insufficient active fill-in candidates')
        endif
        inds = inds_pool(1:nsamples)
        call hpsort(inds)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_fillin

    module subroutine sample4update_missing( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: states(:), updatecnts(:)
        integer :: i, cnt, nptcls
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), states(nptcls), updatecnts(nptcls), source=0)
        cnt = 0
        do i = fromto(1), fromto(2)
            cnt             = cnt + 1
            inds(cnt)       = i
            states(cnt)     = self%o(i)%get_state()
            updatecnts(cnt) = self%o(i)%get_updatecnt()
        end do
        nsamples = count(states > 0 .and. updatecnts == 0)
        inds     = pack(inds, mask=states > 0 .and. updatecnts == 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_missing

    module subroutine sample_balanced_1( self, clssmp, nptcls, l_greedy, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        logical,            intent(in)    :: l_greedy
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds_left(:)
        integer        :: i, j
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0
        if( l_greedy )then
            do i = 1, size(clssmp)
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = 1
                end do
            end do
        else
            do i = 1, size(clssmp)
                allocate(pinds_left(clssmp(i)%pop), source=clssmp(i)%pinds)
                call shuffle(pinds_left)
                do j = 1, clssmp(i)%nsample
                    states(pinds_left(j)) = 1
                end do
                deallocate(pinds_left)
            end do
        endif
    end subroutine sample_balanced_1

    module subroutine sample_balanced_2( self, clssmp, nptcls, frac_best, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        real,               intent(in)    :: frac_best
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds2sample(:)
        integer        :: i, j, nbest
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0    
        do i = 1, size(clssmp)
            nbest = max(clssmp(i)%nsample, nint(frac_best * real(clssmp(i)%pop)))
            if( nbest == clssmp(i)%nsample )then
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = 1
                end do
            else
                allocate(pinds2sample(nbest), source=clssmp(i)%pinds(:nbest))
                call shuffle(pinds2sample)
                do j = 1, clssmp(i)%nsample
                    states(pinds2sample(j)) = 1
                end do
                deallocate(pinds2sample)
            endif
        end do
    end subroutine sample_balanced_2

    module subroutine sample_balanced_inv( self, clssmp, nptcls, frac_worst, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        real,               intent(in)    :: frac_worst
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds_rev(:), pinds2sample(:)
        integer        :: i, j, nworst
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0    
        do i = 1, size(clssmp)
            nworst = max(clssmp(i)%nsample, nint(frac_worst * real(clssmp(i)%pop)))
            if( nworst == clssmp(i)%nsample )then
                do j = clssmp(i)%nsample,1,-1
                    states(clssmp(i)%pinds(j)) = 1
                end do
            else
                allocate(pinds_rev(clssmp(i)%pop), source=clssmp(i)%pinds)
                call reverse(pinds_rev)
                allocate(pinds2sample(nworst), source=pinds_rev(:nworst))
                call shuffle(pinds2sample)
                do j = 1, clssmp(i)%nsample
                    states(pinds2sample(j)) = 1
                end do
                deallocate(pinds_rev, pinds2sample)
            endif
        end do
    end subroutine sample_balanced_inv

    module subroutine sample_balanced_parts( self, clssmp, nparts, states, nptcls_per_part )
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
        integer :: nptcls, i, j, k, nptcls_eff
        nptcls_eff = self%count_state_gt_zero() 
        if( present(nptcls_per_part) )then
            nptcls = min(nptcls_eff, nparts * nptcls_per_part)
        else
            nptcls = nptcls_eff
        endif
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0
        do i = 1, size(clssmp)
            j = 1
            do while( j < clssmp(i)%nsample )
                if( j + nparts - 1 > clssmp(i)%pop ) exit
                do k = 1, nparts
                    states(clssmp(i)%pinds(j)) = k
                    j = j + 1
                end do
            end do
        end do
    end subroutine sample_balanced_parts

    module subroutine sample_ranked_parts( self, clssmp, nparts, states, nptcls_per_part )
        use simple_map_reduce, only: split_nobjs_even
        class(oris),        intent(inout) :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nparts
        integer,            intent(inout) :: states(self%n)
        integer, optional,  intent(in)    :: nptcls_per_part
        integer, allocatable :: parts(:,:)
        integer :: nptcls, i, j, nptcls_eff, ipart
        nptcls_eff = self%count_state_gt_zero() 
        if( present(nptcls_per_part) )then
            nptcls = min(nptcls_eff, nparts * nptcls_per_part)
        else
            nptcls = nptcls_eff
        endif
        clssmp(:)%nsample = 0
        do while( sum(clssmp(:)%nsample) < nptcls )
            where( clssmp(:)%nsample < clssmp(:)%pop ) clssmp(:)%nsample = clssmp(:)%nsample + 1
        end do
        states = 0
        do i = 1, size(clssmp)
            if( clssmp(i)%nsample <= nparts )then
                do j = 1, clssmp(i)%nsample
                    states(clssmp(i)%pinds(j)) = j
                end do
            else
                parts = split_nobjs_even(clssmp(i)%nsample, nparts)
                do ipart = 1, nparts
                    do j = parts(ipart,1), parts(ipart,2)
                        states(clssmp(i)%pinds(j)) = ipart
                    end do
                end do
                deallocate(parts)
            endif
        enddo
    end subroutine sample_ranked_parts

    module function get_sample_ind( self, incr_sampled ) result( sample_ind )
        class(oris), intent(in) :: self
        logical,     intent(in) :: incr_sampled
        integer :: i, sample_ind
        sample_ind = 0
        do i = 1, self%n
            if( self%o(i)%get_state() > 0 ) sample_ind = max(sample_ind, self%o(i)%get_sampled())
        end do
        if( incr_sampled ) sample_ind = sample_ind + 1
    end function get_sample_ind

    module subroutine incr_sampled_updatecnt( self, inds, incr_sampled )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: inds(:)
        logical,     intent(in)    :: incr_sampled
        integer :: i, iptcl, sample_ind
        real    :: val
        sample_ind = self%get_sample_ind(incr_sampled)
        do i = 1, size(inds)
            iptcl = inds(i)
            val   = self%o(iptcl)%get('updatecnt')
            call self%o(iptcl)%set('updatecnt', val + 1.0)
            call self%o(iptcl)%set('sampled',   sample_ind)
        end do
    end subroutine incr_sampled_updatecnt

    module subroutine set_nonzero_updatecnt( self, updatecnt  )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: updatecnt
        integer :: i
        do i = 1,self%n
            if( self%o(i)%get('updatecnt') > 0 )then
                call self%o(i)%set('updatecnt', updatecnt)
            endif
        enddo
    end subroutine set_nonzero_updatecnt

    module subroutine set_updatecnt( self, updatecnt, pinds )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: updatecnt, pinds(:)
        integer :: i, n
        do i = 1, self%n
            call self%o(i)%set('updatecnt', 0)
        enddo
        n = size(pinds)
        do i = 1, n
            call self%o(pinds(i))%set('updatecnt', updatecnt)
        enddo
    end subroutine set_updatecnt

    module subroutine clean_entry( self, varflag1, varflag2 )
        class(oris),                 intent(inout) :: self
        character (len=*),           intent(in)    :: varflag1
        character (len=*), optional, intent(in)    :: varflag2
        logical :: varflag2_present
        integer :: i
        varflag2_present = present(varflag2)
        do i = 1,self%n
            call self%o(i)%delete_entry(varflag1)
            if( varflag2_present ) call self%o(i)%delete_entry(varflag2)
        enddo
    end subroutine clean_entry

    !> an absent allow_empty keeps the hard stop on an empty sample; a distributed
    !! partition passes .true., since its range may hold nothing to update
    pure logical function empty_allowed( allow_empty )
        logical, optional, intent(in) :: allow_empty
        empty_allowed = .false.
        if( present(allow_empty) ) empty_allowed = allow_empty
    end function empty_allowed

end submodule simple_oris_sampling
