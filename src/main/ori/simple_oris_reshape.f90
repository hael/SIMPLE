!@descr: reshaping/partitioning/remapping routines for oris object
submodule (simple_oris) simple_oris_reshape
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine compress( self, mask )
        class(oris), intent(inout) :: self
        logical,     intent(in)    :: mask(:)
        type(oris) :: os_tmp
        integer    :: i, cnt
        logical    :: is_ptcl
        if( size(mask) /= self%n )then
            write(logfhandle,*) 'self%n:     ', self%n
            write(logfhandle,*) 'size(mask): ', size(mask)
            THROW_HARD('nonconforming mask size; compress')
        endif
        is_ptcl = self%o(1)%is_particle()
        call os_tmp%new(count(mask), is_ptcl)
        cnt = 0
        do i=1,self%n
            if( mask(i) )then
                cnt = cnt + 1
                os_tmp%o(cnt) = self%o(i)
            endif
        end do
        call self%copy(os_tmp, is_ptcl)
        call os_tmp%kill
    end subroutine compress

    module subroutine split_state( self, which )
        use simple_ran_tabu, only: ran_tabu
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: states(:)
        type(ran_tabu)       :: rt
        integer              ::  n, nstates, iptcl
        nstates = self%get_n('state')
        if( which < 1 .or. which > nstates )then
            THROW_HARD('which (state) is out of range; split_state')
        endif
        call self%get_pinds(which, 'state', ptcls_in_which)
        n = size(ptcls_in_which)
        allocate(states(n))
        rt = ran_tabu(n)
        call rt%balanced(2, states)
        do iptcl=1,n
            if( states(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set_state(nstates+1)
            endif
        end do
        call rt%kill
        deallocate(ptcls_in_which,states)
    end subroutine split_state

    module subroutine split_class( self, which )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: which
        integer, allocatable :: ptcls_in_which(:)
        integer, allocatable :: members(:)
        type(ran_tabu)       :: rt
        integer              ::  n, nmembers, iptcl
        nmembers = self%get_n('class')
        if( which < 1 .or. which > nmembers )then
            THROW_HARD('which member is out of range; split_class')
        endif
        call self%get_pinds(which, 'class', ptcls_in_which)
        n = size(ptcls_in_which)
        allocate(members(n))
        rt = ran_tabu(n)
        call rt%balanced(2, members)
        do iptcl=1,n
            if( members(iptcl) == 1 )then
                ! do nothing, leave this state as is
            else
                call self%o(ptcls_in_which(iptcl))%set_class(nmembers+1)
            endif
        end do
        call rt%kill
        deallocate(ptcls_in_which,members)
    end subroutine split_class

    module subroutine expand_classes( self, ncls_target )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: ncls_target
        integer, allocatable       :: pops(:)
        integer :: ncls, loc(1), myncls, icls
        ncls = self%get_n('class')
        if( ncls_target <= ncls ) THROW_HARD('nr of target classes cannot be <= original number')
        allocate(pops(ncls_target))
        pops = 0
        do icls=1,ncls
            pops(icls) = self%get_pop(icls, 'class')
        end do
        myncls = ncls
        do while( myncls < ncls_target )
            loc = maxloc(pops)
            call self%split_class(loc(1))
            myncls = myncls+1
            pops(loc(1)) = self%get_pop(loc(1), 'class')
            pops(myncls) = self%get_pop(myncls, 'class')
        end do
    end subroutine expand_classes

    module subroutine remap_cls( self )
        class(oris), intent(inout) :: self
        integer :: ncls, clsind_remap, pop, icls, iptcl, old_cls
        integer , allocatable :: clspops(:)
        ncls = self%get_n('class')
        allocate(clspops(ncls))
        do icls=1,ncls
            clspops(icls) = self%get_pop(icls, 'class')
        end do
        if( any(clspops == 0) )then
            clsind_remap = ncls
            do icls=1,ncls
                pop = clspops(icls)
                if( pop > 1 )then
                    clsind_remap = clsind_remap + 1
                    do iptcl=1,self%n
                        old_cls = self%o(iptcl)%get_class()
                        if( old_cls == icls ) call self%o(iptcl)%set_class(clsind_remap)
                    end do
                else
                    do iptcl=1,self%n
                        old_cls = self%o(iptcl)%get_class()
                        if( old_cls == icls )then
                            call self%o(iptcl)%set_class(0)
                            call self%o(iptcl)%set_state(0)
                        endif
                    end do
                endif
            end do
            do iptcl=1,self%n
                old_cls = self%o(iptcl)%get_class()
                if( old_cls /= 0 ) call self%o(iptcl)%set_class(old_cls-ncls)
            end do
        endif
        deallocate(clspops)
    end subroutine remap_cls

    !> Re-seed a K-class partition from a previous clustering, metadata only.
    !! clsinds are the accepted parent classes. Each parent receives a number
    !! of seed classes proportional to its active population (largest-remainder
    !! allocation; at least one when ncls_target >= size(clsinds)), so every
    !! seed class has ~nptcls/ncls_target particles (balanced) and the seed
    !! set represents the previous view distribution (representative). A
    !! parent with more than one seed class is split by rank interleaving on
    !! corr (best-to-worst, dealt round-robin), so the children are equal in
    !! size and in objective-value distribution. When ncls_target < size(clsinds)
    !! the least populous parents receive no seed class. Every active particle
    !! outside the seeded parents (dropped, rejected, or unlabelled) gets
    !! class=0; e3/shift are left untouched. State=0 particles are ignored.
    !! parent_of_seed(k) is the parent class of seed class k, seed_pops(k) its
    !! population, ndropped the number of accepted parents left without a seed.
    module subroutine reseed_classes( self, clsinds, ncls_target, parent_of_seed, seed_pops, ndropped )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: clsinds(:)
        integer,              intent(in)    :: ncls_target
        integer, allocatable, intent(inout) :: parent_of_seed(:), seed_pops(:)
        integer,              intent(out)   :: ndropped
        type(class_sample), allocatable :: clssmp(:)
        integer, allocatable :: nseeds(:)
        real,    allocatable :: quota(:), remainder(:)
        integer :: nparents, nptcls_tot, i, j, k, iseed, offset, nalloc, loc(1)
        nparents = size(clsinds)
        if( nparents   < 1 ) THROW_HARD('reseed_classes requires at least one accepted parent class')
        if( ncls_target < 1 ) THROW_HARD('reseed_classes requires ncls_target >= 1')
        if( .not. self%isthere('class') ) THROW_HARD('reseed_classes requires a previous class assignment')
        ! parents: active particles sorted best-to-worst by corr
        call self%get_class_sample_stats(clsinds, clssmp)
        nptcls_tot = 0
        do i = 1, nparents
            if( .not. allocated(clssmp(i)%pinds) ) clssmp(i)%pop = 0
            nptcls_tot = nptcls_tot + clssmp(i)%pop
        end do
        if( nptcls_tot < 1 ) THROW_HARD('reseed_classes: accepted parent classes hold no active particles')
        ! largest-remainder allocation of seed classes to parents
        allocate(nseeds(nparents), quota(nparents), remainder(nparents))
        do i = 1, nparents
            quota(i)  = real(ncls_target) * real(clssmp(i)%pop) / real(nptcls_tot)
            nseeds(i) = floor(quota(i))
            if( ncls_target >= nparents .and. clssmp(i)%pop > 0 ) nseeds(i) = max(1, nseeds(i))
        end do
        nalloc = sum(nseeds)
        do while( nalloc < ncls_target )
            ! largest remainder first, ties by larger population
            remainder = quota - real(nseeds)
            where( clssmp(:)%pop <= nseeds(:) ) remainder = -huge(1.0)  ! cannot seed more classes than particles
            remainder = remainder + 1.e-6 * real(clssmp(:)%pop) / real(nptcls_tot)
            loc = maxloc(remainder)
            if( remainder(loc(1)) <= -huge(1.0)/2. ) exit
            nseeds(loc(1)) = nseeds(loc(1)) + 1
            nalloc = nalloc + 1
        end do
        do while( nalloc > ncls_target )
            ! the floor-to-one guarantee overshot: take from the parent that least deserves its last seed
            remainder = quota - real(nseeds)
            where( nseeds <= 1 ) remainder = huge(1.0)
            remainder = remainder + 1.e-6 * real(clssmp(:)%pop) / real(nptcls_tot) ! ties: take from the smaller parent
            loc = minloc(remainder)
            if( remainder(loc(1)) >= huge(1.0)/2. ) THROW_HARD('reseed_classes: cannot reduce the seed allocation')
            nseeds(loc(1)) = nseeds(loc(1)) - 1
            nalloc = nalloc - 1
        end do
        if( nalloc /= ncls_target ) THROW_HARD('reseed_classes: seed allocation does not match ncls_target')
        ndropped = count(nseeds == 0)
        ! every particle starts unassigned; seeded parents are relabelled below
        call self%set_all2single('class', 0)
        if( allocated(parent_of_seed) ) deallocate(parent_of_seed)
        if( allocated(seed_pops)      ) deallocate(seed_pops)
        allocate(parent_of_seed(ncls_target), seed_pops(ncls_target), source=0)
        offset = 0
        do i = 1, nparents
            if( nseeds(i) == 0 ) cycle
            do k = 1, nseeds(i)
                parent_of_seed(offset + k) = clsinds(i)
            end do
            ! rank interleaving: rank j (best first) -> child mod(j-1,nseeds)+1
            do j = 1, clssmp(i)%pop
                iseed = offset + mod(j - 1, nseeds(i)) + 1
                call self%o(clssmp(i)%pinds(j))%set_class(iseed)
                seed_pops(iseed) = seed_pops(iseed) + 1
            end do
            offset = offset + nseeds(i)
        end do
        do i = 1, nparents
            if( allocated(clssmp(i)%pinds) ) deallocate(clssmp(i)%pinds)
            if( allocated(clssmp(i)%ccs)   ) deallocate(clssmp(i)%ccs)
        end do
        deallocate(clssmp, nseeds, quota, remainder)
    end subroutine reseed_classes

    module subroutine merge_classes( self, class_merged, class )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class_merged, class
        integer                    :: i, clsnr
        do i=1,self%n
            clsnr = self%get_class(i)
            if(clsnr == class) call self%set_class(i, class_merged)
        end do
    end subroutine merge_classes

    module subroutine discretize( self, n )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: n
        type(oris) :: d
        integer    :: closest, i
        if( n < self%n )then
            call d%new(n, self%o(1)%is_particle())
            call d%spiral
            do i=1,self%n
                closest = d%find_closest_proj(self%o(i))
                call self%o(i)%e1set(d%e1get(closest))
                call self%o(i)%e2set(d%e2get(closest))
                call self%o(i)%set_class(closest)
            end do
        else
            THROW_HARD('the number of discrete oris is too large; discretize')
        endif
    end subroutine discretize

    module subroutine extract_subspace( self, lnns, subself )
        class(oris), intent(in)    :: self
        logical,     intent(in)    :: lnns(self%n)
        class(oris), intent(inout) :: subself
        integer :: n, cnt, i
        n = count(lnns)
        if( n < 1 ) THROW_HARD('logical array for subspace generation empty')
        call subself%new(n, is_ptcl=self%o(1)%is_particle())
        cnt = 0
        do i = 1, self%n
            if( lnns(i) )then
                cnt            = cnt + 1
                subself%o(cnt) = self%o(i)
            endif
        enddo
    end subroutine extract_subspace

end submodule simple_oris_reshape
