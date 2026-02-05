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
