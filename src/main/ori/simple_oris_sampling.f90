!@descr: sampling and updatecnt related routines for oris
submodule (simple_oris) simple_oris_sampling
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine select_particles_set( self, nptcls, inds )
        class(oris),          intent(in)    :: self
        integer,              intent(in)    :: nptcls
        integer, allocatable, intent(inout) :: inds(:)
        type(ran_tabu) :: rt
        integer        :: i,n
        if( allocated(inds) ) deallocate(inds)
        inds = (/(i,i=1,self%n)/)
        inds = pack(inds, mask=self%get_all('state') > 0.5)
        n    = size(inds)
        if( n < nptcls )then
            ! if less than desired particles select all
        else
            rt   = ran_tabu(n)
            call rt%shuffle(inds)
            call rt%kill
            inds = inds(1:nptcls)
            call hpsort(inds)
        endif
    end subroutine select_particles_set

    module subroutine sample4rec( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        integer, allocatable :: states(:), updatecnts(:)
        integer :: i, cnt, nptcls
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
        if( any(updatecnts > 0) )then
            nsamples = count(states > 0 .and. updatecnts > 0)
            inds     = pack(inds, mask=states > 0 .and. updatecnts > 0)
        else
            nsamples = count(states > 0)
            inds     = pack(inds, mask=states > 0)
        endif
    end subroutine sample4rec

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
        type(ran_tabu) :: rt
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
        inds     = pack(inds,       mask=states > 0)
        nsamples = min(nptcls, nint(update_frac * real(nptcls)))
        rt       = ran_tabu(nptcls)
        call rt%shuffle(inds)
        call rt%kill
        inds = inds(1:nsamples)
        call hpsort(inds)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_rnd

    module subroutine sample4update_cnt( self, fromto, update_frac, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        type(ran_tabu)       :: rt
        integer, allocatable :: states(:), updatecnts(:)
        integer              :: i, cnt, nptcls, ucnt, ucnt_min, ucnt_max, ucnt_lim
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
        inds        = pack(inds,       mask=states > 0)
        updatecnts  = pack(updatecnts, mask=states > 0)
        ucnt_min    = minval(updatecnts)
        ucnt_max    = maxval(updatecnts)
        nsamples    = min(nptcls, nint(update_frac * real(nptcls)))
        ucnt_lim    = 0
        if( ucnt_max > ucnt_min )then
            do ucnt = ucnt_min,ucnt_max
                if( count(updatecnts < ucnt) >= nsamples )then
                    ucnt_lim = ucnt
                    exit
                endif
            end do
        endif
        if( ucnt_lim > 0 )then
            inds   = pack(inds, mask=updatecnts < ucnt)
            nptcls = size(inds)
        endif
        rt = ran_tabu(nptcls)
        call rt%shuffle(inds)
        call rt%kill
        inds = inds(1:nsamples)
        call hpsort(inds)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_cnt

    module subroutine sample4update_class( self, clssmp, fromto, update_frac, nsamples, inds, incr_sampled, l_greedy, frac_best )
        class(oris),          intent(inout) :: self
        type(class_sample),   intent(inout) :: clssmp(:)
        integer,              intent(in)    :: fromto(2)
        real,                 intent(in)    :: update_frac
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled, l_greedy
        real,    optional,    intent(in)    :: frac_best
        integer, allocatable :: states(:)
        real,    allocatable :: rstates(:)
        integer :: i, cnt, nptcls, nsamples_class, states_bal(self%n)
        rstates        = self%get_all('state')
        nsamples_class = nint(update_frac * real(count(rstates > 0.5)))
        deallocate(rstates)
        if( present(frac_best) )then
            call self%sample_balanced(clssmp, nsamples_class, frac_best,  states_bal)
        else
            call self%sample_balanced(clssmp, nsamples_class, l_greedy, states_bal)
        endif
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
        inds     = pack(inds, mask=states > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_class

    module subroutine sample4update_reprod( self, fromto, nsamples, inds )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        integer, allocatable :: sampled(:)
        integer :: i, cnt, nptcls, sample_ind
        nptcls = fromto(2) - fromto(1) + 1
        if( allocated(inds) ) deallocate(inds)
        allocate(inds(nptcls), sampled(nptcls), source=0)
        cnt        = 0
        sample_ind = self%get_sample_ind(.false.)
        do i = fromto(1), fromto(2)
            cnt          = cnt + 1
            inds(cnt)    = i
            sampled(cnt) = self%o(i)%get_sampled()
        end do
        if( sample_ind == 0 ) THROW_HARD('requires previous sampling')
        nsamples = count(sampled == sample_ind)
        inds     = pack(inds, mask=sampled == sample_ind)
    end subroutine sample4update_reprod

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
        inds     = pack(inds, mask=updatecnts > 0)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_updated

    module subroutine sample4update_fillin( self, fromto, nsamples, inds, incr_sampled )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: fromto(2)
        integer,              intent(inout) :: nsamples
        integer, allocatable, intent(inout) :: inds(:)
        logical,              intent(in)    :: incr_sampled
        integer, allocatable :: updatecnts(:), states(:)
        integer :: i, cnt, nptcls, minucnt
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
        updatecnts = pack(updatecnts, mask=states > 0)
        inds       = pack(inds,       mask=states > 0)
        minucnt    = minval(updatecnts)
        nsamples   = count(updatecnts == minucnt)
        inds       = pack(inds, mask=updatecnts == minucnt)
        call self%incr_sampled_updatecnt(inds, incr_sampled)
    end subroutine sample4update_fillin

    module subroutine sample_balanced_1( self, clssmp, nptcls, l_greedy, states )
        class(oris),        intent(in)    :: self
        type(class_sample), intent(inout) :: clssmp(:)
        integer,            intent(in)    :: nptcls
        logical,            intent(in)    :: l_greedy
        integer,            intent(inout) :: states(self%n)
        integer,            allocatable   :: pinds_left(:)
        type(ran_tabu) :: rt
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
                rt = ran_tabu(clssmp(i)%pop)
                call rt%shuffle(pinds_left)
                call rt%kill
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
        type(ran_tabu) :: rt
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
                rt = ran_tabu(nbest)
                call rt%shuffle(pinds2sample)
                call rt%kill
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
        type(ran_tabu) :: rt
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
                rt = ran_tabu(nworst)
                call rt%shuffle(pinds2sample)
                call rt%kill
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

    module subroutine balance_ptcls_within_cls( self, nptcls, pinds, maxpop, nparts )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: nptcls, maxpop, nparts
        integer,     intent(in)    :: pinds(1:nptcls)
        real,    allocatable :: cls_scores(:)
        integer, allocatable :: cls_inds(:)
        integer :: classes(nptcls), i,icls,pop,n2reject,ncls,mpop
        if( maxpop < 1 )return
        classes = self%o(pinds(:))%get_class()
        ncls    = maxval(classes)
        mpop    = floor(real(maxpop)/real(nparts))
        do icls = 1,ncls
            cls_inds = pack((/(i,i=1,nptcls)/), mask=classes==icls)
            if( .not.allocated(cls_inds) ) cycle
            pop = size(cls_inds)
            if( pop <= mpop )cycle
            cls_scores = self%o(pinds(cls_inds(:)))%get('corr')
            call hpsort(cls_scores, cls_inds)
            n2reject = pop - mpop
            do i = 1,n2reject
                call self%o(pinds(cls_inds(i)))%set('w',0.)
            enddo
            deallocate(cls_inds,cls_scores)
        enddo
    end subroutine balance_ptcls_within_cls

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

end submodule simple_oris_sampling
