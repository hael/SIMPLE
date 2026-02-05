! @descr: getters/queries for oris object
submodule (simple_oris) simple_oris_getters
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module pure function exists( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        exists = self%o(i)%exists()
    end function exists

    module pure function e1get( self, i ) result( e1 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e1
        e1 = self%o(i)%e1get()
    end function e1get

    module pure function e2get( self, i ) result( e2 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e2
        e2 = self%o(i)%e2get()
    end function e2get

    module pure function e3get( self, i ) result( e3 )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: e3
        e3 = self%o(i)%e3get()
    end function e3get

    module pure function get_euler( self, i ) result( euls )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: euls(3)
        euls = self%o(i)%get_euler()
    end function get_euler

    module pure function get_noris( self, consider_state ) result( n )
        class(oris),       intent(in) :: self
        logical, optional, intent(in) :: consider_state
        integer :: i, n
        logical :: consider_state_here
        consider_state_here = .false.
        if(present(consider_state)) consider_state_here = consider_state
        if( consider_state_here )then
            n = 0
            do i=1,self%n
                if( self%o(i)%isthere('state') )then
                    if( self%o(i)%get('state') > 0.5 ) n = n+1
                else
                    n = n+1 ! included by default
                endif
            enddo
        else
            n = self%n
        endif
    end function get_noris

    module subroutine get_ori( self, i, o )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        type(ori),   intent(inout) :: o
        if( self%n == 0 ) THROW_HARD('oris object does not exist; get_ori')
        if( i > self%n .or. i < 1 )then
            write(logfhandle,*) 'trying to get ori: ', i, ' among: ', self%n, ' oris'
            THROW_HARD('i out of range; get_ori')
        endif
        o = self%o(i)
    end subroutine get_ori

    module pure function get( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        get = self%o(i)%get(key)
    end function get

    module pure function get_int( self, i, key )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        call self%o(i)%getter(key, get_int)
    end function get_int

    module function get_str( self, i, key ) result( val )
        class(oris),      intent(in) :: self
        integer,          intent(in) :: i
        character(len=*), intent(in) :: key
        type(string) :: val
        call self%getter_1(i, key, val)
    end function get_str

    module pure subroutine get_static( self, i, key, val )
        class(oris),      intent(in)  :: self
        integer,          intent(in)  :: i
        character(len=*), intent(in)  :: key
        character(len=*), intent(out) :: val
        call self%o(i)%get_static(key, val)
    end subroutine get_static

    module subroutine getter_1( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        type(string),     intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_1

    module subroutine getter_2( self, i, key, val )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        real,             intent(inout) :: val
        call self%o(i)%getter(key, val)
    end subroutine getter_2

    module subroutine getter_3( self, i, key, ival )
        class(oris),      intent(in)    :: self
        integer,          intent(in)    :: i
        character(len=*), intent(in)    :: key
        integer,          intent(inout) :: ival
        call self%o(i)%getter(key, ival)
    end subroutine getter_3

    module function get_all( self, key, fromto, nonzero ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        logical, optional, intent(in) :: nonzero
        real,    allocatable :: arr(:)
        integer :: ffromto(2)
        logical :: nnonzero
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        nnonzero = .false.
        if( present(nonzero) ) nnonzero = nonzero
        allocate( arr(ffromto(1):ffromto(2)), source=self%o(ffromto(1):ffromto(2))%get(key) )
        if( nnonzero )then
            arr = pack(arr, mask=self%o(ffromto(1):ffromto(2))%get_state()>0)
        endif
    end function get_all

    module function get_all_asint( self, key, fromto ) result( iarr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: fromto(2)
        integer, allocatable :: iarr(:)
        integer :: ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        allocate(iarr(ffromto(1):ffromto(2)), source=nint(self%o(ffromto(1):ffromto(2))%get(key)))
    end function get_all_asint

    module function gen_ptcl_mask( self, key, ival, frac_best ) result( mask )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer,           intent(in) :: ival
        real,    optional, intent(in) :: frac_best
        logical, allocatable :: mask(:)
        integer, allocatable :: ivals(:)
        real,    allocatable :: corrs(:), tmp(:)
        real    :: corr_t
        integer :: ncorrs, nsample
        allocate(mask(self%n))
        allocate(ivals(self%n), source=self%o(1:self%n)%get_int(key))
        if( present(frac_best) )then
            allocate(corrs(self%n), source=self%o(self%n)%get('corr'))
            tmp     = pack(corrs, mask=ivals == ival)
            ncorrs  = size(tmp)
            nsample = ceiling(frac_best * real(ncorrs))
            call hpsort(tmp)
            corr_t  = tmp(ncorrs - nsample)
            deallocate(tmp)
            where(ivals == ival .and. corrs >= corr_t )
                mask = .true.
            elsewhere
                mask = .false.
            endwhere
        else
            where(ivals == ival )
                mask = .true.
            elsewhere
                mask = .false.
            endwhere
        endif
    end function gen_ptcl_mask

    module function get_all_sampled( self, key, state, lowerbound ) result( arr )
        class(oris),       intent(in) :: self
        character(len=*),  intent(in) :: key
        integer, optional, intent(in) :: state
        real,    optional, intent(in) :: lowerbound
        real,    allocatable :: arr(:), sampled(:)
        integer, allocatable :: states(:)
        integer :: i
        real    :: lb
        allocate(arr(self%n), sampled(self%n), source=0.)
        if( present(state) ) states = nint(self%get_all('state'))
        do i=1,self%n
            arr(i)     = self%o(i)%get(key)
            sampled(i) = self%o(i)%get('sampled')
        enddo
        if( present(lowerbound) )then
            lb = lowerbound
        else
            lb = maxval(sampled) - 0.5
        endif
        if( present(state) )then
            arr = pack(arr, mask=(sampled > lb .and. states == state))
        else
            arr = pack(arr, mask=sampled > lb)
        endif
    end function get_all_sampled

    module pure function get_mat( self, i ) result( mat )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: mat(3,3)
        mat = self%o(i)%get_mat()
    end function get_mat

    module pure function get_normal( self, i ) result( normal )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: normal(3)
        normal = self%o(i)%get_normal()
    end function get_normal

    module pure function get_2Dshift( self, i )  result( shvec )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        real :: shvec(2)
        shvec = self%o(i)%get_2Dshift()
    end function get_2Dshift

    module pure function get_dfx( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_dfx = self%o(i)%get_dfx()
    end function get_dfx

    module pure function get_dfy( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_dfy = self%o(i)%get_dfy()
    end function get_dfy

    module pure function get_state( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_state = self%o(i)%get_state()
    end function get_state

    module pure function get_class( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_class = self%o(i)%get_class()
    end function get_class

    module pure function get_proj( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_proj = self%o(i)%get_proj()
    end function get_proj

    module function get_label_inds( self, label ) result( inds )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: i, icls, icls_max
        integer, allocatable :: inds(:), inds_all(:)
        logical, allocatable :: isthere(:)
        icls_max = self%get_n(label)
        inds_all = (/(i,i=1,icls_max)/)
        allocate(isthere(icls_max), source=.false.)
        do i = 1, self%n
            icls = self%get_class(i)
            isthere(icls) = .true.
        end do
        inds = pack(inds_all, mask=isthere)
    end function get_label_inds

    module pure function get_eo( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_eo = self%o(i)%get_eo()
    end function get_eo

    module pure function get_fromp( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_fromp = self%o(i)%get_fromp()
    end function get_fromp

    module pure function get_top( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_top = self%o(i)%get_top()
    end function get_top

    module pure function get_sampled( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_sampled = self%o(i)%get_sampled()
    end function get_sampled

    module pure function get_updatecnt( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        get_updatecnt = self%o(i)%get_updatecnt()
    end function get_updatecnt

    module subroutine get_tseries_neighs( self, nsz, ptcls2neigh )
       class(oris),          intent(in)    :: self
       integer,              intent(in)    :: nsz
       integer, allocatable, intent(inout) :: ptcls2neigh(:,:)
       integer :: i, j, cls1, cls2
       if( allocated(ptcls2neigh) ) deallocate(ptcls2neigh)
       allocate(ptcls2neigh(self%n,2), source=0)
       do i=1,self%n-nsz
           cls1 = self%o(i)%get_class()
           cls2 = self%o(i+1)%get_class()
           if( cls2 == cls1 )then
               cycle
           else
               do j=max(1,i-nsz+1),min(self%n,i+nsz)
                   ptcls2neigh(j,1) = cls1
                   ptcls2neigh(j,2) = cls2
               end do
           endif
       end do
    end subroutine get_tseries_neighs

    module pure function isthere_1( self, key ) result( is )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: key
        logical :: is
        integer :: i
        is = .false.
        do i=1,self%n
            is = self%o(i)%isthere(key)
            if( is ) exit
        end do
    end function isthere_1

    module pure function isthere_2( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
        is = self%o(i)%isthere(key)
    end function isthere_2

    module function ischar( self, i, key ) result( is )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: i
        character(len=*),  intent(in) :: key
        logical :: is
        is = self%o(i)%ischar(key)
    end function ischar

    module pure function is_particle( self ) result( t )
        class(oris), intent(in) :: self
        logical :: t
        t = self%o(1)%is_particle()
    end function is_particle

    module function max_ori_strlen_trim( self )
        class(oris), intent(in) :: self
        integer :: i
        max_ori_strlen_trim = 0
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)&
        !$omp reduction(max:max_ori_strlen_trim)
        do i=1,self%n
            max_ori_strlen_trim = max(max_ori_strlen_trim,self%o(i)%ori_strlen_trim())
        end do
        !$omp end parallel do
    end function max_ori_strlen_trim

    module function get_n( self, label ) result( n )
        class(oris),      intent(in) :: self
        character(len=*), intent(in) :: label
        integer :: i, n, ival
        n = 1
        do i=1,self%n
            ival = self%o(i)%get_int(label)
            if( ival > n ) n = ival
        end do
    end function get_n

    module function get_pop_1( self, ind, label, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: ind
        character(len=*),  intent(in) :: label
        integer, optional, intent(in) :: eo
        integer :: pop, i
        logical :: consider_eo
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        pop = 0
        if( consider_eo )then
            !$omp parallel do private(i) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero()  ) cycle
                if( self%o(i)%get_eo() /= eo ) cycle
                if( self%o(i)%get_int(label) == ind )  pop = pop + 1
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(i) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero() ) cycle
                if( self%o(i)%get_int(label) == ind )  pop = pop + 1
            end do
            !$omp end parallel do
        endif
    end function get_pop_1

    module function get_pop_2( self, inds, labels, eo ) result( pop )
        class(oris),       intent(in) :: self
        integer,           intent(in) :: inds(:)
        character(len=*),  intent(in) :: labels(:)
        integer, optional, intent(in) :: eo
        integer :: pop, i, j, n
        logical :: consider_eo, l_okay
        consider_eo = .false.
        if( present(eo) ) consider_eo = .true.
        n   = size(inds)
        pop = 0
        if( consider_eo )then
            !$omp parallel do private(i,j,l_okay) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero()  ) cycle
                if( self%o(i)%get_eo() /= eo ) cycle
                l_okay = .true.
                do j=1,n
                    if( self%o(i)%get_int(labels(j)) /= inds(j) )then
                        l_okay = .false.
                        exit
                    endif
                enddo
                if( l_okay )  pop = pop + 1
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(i,j,l_okay) default(shared) proc_bind(close) reduction(+:pop)
            do i=1,self%n
                if( self%o(i)%isstatezero() ) cycle
                l_okay = .true.
                do j=1,n
                    if( self%o(i)%get_int(labels(j)) /= inds(j) )then
                        l_okay = .false.
                        exit
                    endif
                enddo
                if( l_okay )  pop = pop + 1
            end do
            !$omp end parallel do
        endif
    end function get_pop_2

    module subroutine get_pops( self, pops, label, maxn, weight )
        class(oris),          intent(in)    :: self
        integer, allocatable, intent(out)   :: pops(:)
        character(len=*),     intent(in)    :: label
        integer, optional,    intent(in)    :: maxn
        logical, optional,    intent(in)    :: weight
        integer :: i, val, n
        logical :: consider_w
        consider_w = .false.
        if( present(weight) ) consider_w = weight
        n = self%get_n(label)
        if( present(maxn) ) n = max(n, maxn)
        if(allocated(pops)) deallocate(pops)
        allocate(pops(n),source=0)
        !$omp parallel do private(i,val) default(shared) proc_bind(close) reduction(+:pops)
        do i = 1,self%n
            if( self%o(i)%isstatezero() ) cycle
            if( consider_w )then
                if( self%o(i)%get('w') < 1.e-6 ) cycle
            endif
            val = self%o(i)%get_int(label)
            if( val > 0 ) pops(val) = pops(val) + 1
        end do
        !$omp end parallel do
    end subroutine get_pops

    module subroutine get_pinds( self, ind, label, indices, l_shuffle, l_require_updated )
        class(oris),          intent(in)  :: self
        character(len=*),     intent(in)  :: label
        integer,              intent(in)  :: ind
        integer, allocatable, intent(out) :: indices(:)
        logical, optional,    intent(in)  :: l_shuffle, l_require_updated
        type(ran_tabu)       :: rt
        logical, allocatable :: mask(:)
        integer :: pop, i
        logical :: ll_shuffle, ll_require_updated
        ll_shuffle = .false.
        if( present(l_shuffle) ) ll_shuffle = l_shuffle
        ll_require_updated = .false.
        if( present(l_require_updated) ) ll_require_updated = l_require_updated
        if( allocated(indices) )deallocate(indices)
        allocate(indices(self%n),mask(self%n))
        !$omp parallel do private(i) default(shared) proc_bind(close)
        do i = 1,self%n
            if( self%o(i)%isstatezero() )then
                mask(i) = .false.
            else
                if( ll_require_updated )then
                    if( self%o(i)%get_int('updatecnt') == 0 )then
                        mask(i) = .false.
                    else
                        mask(i) = self%o(i)%get_int(label) == ind
                    endif
                else
                    mask(i) = self%o(i)%get_int(label) == ind
                endif
            endif
            if( mask(i) ) indices(i) = i
        end do
        !$omp end parallel do
        pop = count(mask)
        if( pop > 0 )then
            indices = pack(indices, mask=mask)
            deallocate(mask)
            if( ll_shuffle )then
                rt  = ran_tabu(pop)
                call rt%shuffle(indices)
                call rt%kill
            endif
        else
            deallocate(indices,mask)
        endif
    end subroutine get_pinds

    module subroutine gen_mask( self, state, ind, label, l_mask, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state, ind
        character(len=*),     intent(in)    :: label
        logical, allocatable, intent(out)   :: l_mask(:)
        integer, optional,    intent(in)    :: fromto(2)
        integer :: i, mystate, myval, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if( allocated(l_mask) ) deallocate(l_mask)
        allocate(l_mask(ffromto(1):ffromto(2)))
        l_mask = .false.
        do i=ffromto(1),ffromto(2)
            mystate = self%o(i)%get_state()
            myval   = self%o(i)%get_int(label)
            if( mystate == state .and. myval == ind ) l_mask(i) = .true.
        end do
    end subroutine gen_mask

    module subroutine mask_from_state( self, state, l_mask, pinds, fromto )
        class(oris),          intent(inout) :: self
        integer,              intent(in)    :: state
        logical, allocatable, intent(inout) :: l_mask(:)
        integer, allocatable, intent(inout) :: pinds(:)
        integer, optional,    intent(in)    :: fromto(2)
        integer :: i, cnt, ffromto(2)
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        if( allocated(l_mask) ) deallocate(l_mask)
        if( allocated(pinds)  ) deallocate(pinds)
        allocate(l_mask(ffromto(1):ffromto(2)))
        l_mask = .false.
        do i=ffromto(1),ffromto(2)
            if( self%o(i)%get_state() == state ) l_mask(i)=.true.
        end do
        allocate(pinds(1:count(l_mask)))
        cnt = 0
        do i=ffromto(1),ffromto(2)
            if( .not.l_mask(i) ) cycle
            cnt = cnt + 1
            pinds(cnt) = i
        end do
    end subroutine mask_from_state

    module function get_all_normals( self ) result( normals )
        class(oris), intent(inout) :: self
        real, allocatable :: normals(:,:)
        integer :: i
        allocate(normals(self%n,3))
        do i=1,self%n
            normals(i,:) = self%o(i)%get_normal()
        end do
    end function get_all_normals

    module function states_exist( self, nstates, thres ) result( state_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer, optional, intent(in)    :: thres
        integer :: i, min_pop
        logical :: state_exists(nstates)
        min_pop = 0
        if( present(thres) ) min_pop = thres
        do i=1,nstates
            state_exists(i) = (self%get_pop(i, 'state') > min_pop)
        end do
    end function states_exist

    module function projs_exist( self, nstates, nprojs, thres ) result( proj_exists )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: nstates
        integer,           intent(in)    :: nprojs
        integer, optional, intent(in)    :: thres
        integer :: i, j, min_pop
        logical :: proj_exists(nprojs,nstates)
        min_pop = 0
        if( present(thres) ) min_pop = thres
        do i=1,nstates
            do j=1,nprojs
                proj_exists(j,i) = (self%get_pop([i,j], ['state', 'proj ']) > min_pop)
            enddo
        end do
    end function projs_exist

    module function get_arr( self, which, class, state ) result( vals )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        real, allocatable :: vals(:)
        integer :: pop, cnt, clsnr, i, mystate
        real    :: val
        logical :: class_present, state_present
        class_present = present(class)
        state_present = present(state)
        if( class_present )then
            pop = self%get_pop(class, 'class')
        else if( state_present )then
            pop = self%get_pop(state, 'state')
        else
            pop = self%n
        endif
        if( pop > 0 )then
            allocate( vals(pop) )
            cnt = 0
            do i=1,self%n
                val = self%get(i, which)
                if( class_present )then
                    mystate = self%o(i)%get_state()
                    if( mystate > 0 )then
                        clsnr   = self%o(i)%get_class()
                        if( clsnr == class )then
                            cnt = cnt+1
                            vals(cnt) = val
                        endif
                    endif
                else if( state_present )then
                    mystate = self%o(i)%get_state()
                    if( mystate == state )then
                        cnt = cnt+1
                        vals(cnt) = val
                    endif
                else
                    vals(i) = val
                endif
            end do
        endif
    end function get_arr

    module subroutine calc_sum( self, which, sum, cnt, class, state, fromto, mask )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: sum
        integer,           intent(out)   :: cnt
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: clsnr, i, mystate, istart, istop
        real    :: val
        logical :: proceed, class_present, state_present, mask_present
        class_present = present(class)
        state_present = present(state)
        mask_present  = present(mask)
        cnt = 0
        sum = 0.
        if( mask_present )then
            if( count(mask)==0 )then
                write(logfhandle,*)'Empty mask; simple_oris :: calc_sum'
                return
            endif
        endif
        if( present(fromto) )then
            istart = fromto(1)
            istop  = fromto(2)
        else
            istart = 1
            istop  = self%n
        endif
        do i=istart,istop
            mystate = self%o(i)%get_state()
            if( mystate == 0 ) cycle
            if( mask_present )then
                proceed = mask(i)
            else
                proceed = .true.
            endif
            if( proceed )then
                val = self%get(i, which)
                if( .not. is_a_number(val) ) val=0.
                if( class_present )then
                    clsnr = self%o(i)%get_class()
                    if( clsnr == class )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else if( state_present )then
                    if( mystate == state )then
                        cnt = cnt+1
                        sum = sum+val
                    endif
                else
                    cnt = cnt+1
                    sum = sum+val
                endif
            endif
        end do
    end subroutine calc_sum

    module function get_sum( self, which, class, state, fromto, mask) result( sum )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: cnt
        real    :: sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
    end function get_sum

    module function get_avg( self, which, class, state, fromto, mask ) result( avg )
        class(oris),       intent(inout) :: self
        character(len=*),  intent(in)    :: which
        integer, optional, intent(in)    :: class
        integer, optional, intent(in)    :: state
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: mask(self%n)
        integer :: cnt
        real    :: avg, sum
        call self%calc_sum(which, sum, cnt, class, state, fromto, mask)
        avg = sum/real(cnt)
    end function get_avg

    module function included( self ) result( incl )
        class(oris),       intent(inout) :: self
        logical, allocatable :: incl(:)
        integer :: i
        if(.not.allocated(incl))allocate(incl(self%n))
        incl = .false.
        do i=1,self%n
            incl(i) = self%o(i)%get_state() > 0
        end do
    end function included

    module function get_neven( self )
        class(oris), intent(inout) :: self
        integer :: i
        get_neven = 0
        do i=1,self%n
            if(self%o(i)%isthere('eo'))then
                if( self%o(i)%get_eo()==0) get_neven = get_neven + 1
            endif
        enddo
    end function get_neven

    module function get_nodd( self )
        class(oris), intent(inout) :: self
        integer, allocatable :: eopart(:)
        eopart = nint(self%get_all('eo'))
        get_nodd = count(eopart == 1)
    end function get_nodd

    module function get_nevenodd( self )
        class(oris),       intent(inout) :: self
        get_nevenodd = self%get_neven() + self%get_nodd()
    end function get_nevenodd

    module function get_all_rmats( self ) result( mat )
        class(oris), intent(in) :: self
        real, allocatable       :: mat(:,:,:)
        integer :: i,n
        n = self%n
        if(allocated(mat))deallocate(mat)
        allocate(mat(n,3,3),source=0.)
        do i=1,self%n
            mat(i,:,:) = self%o(i)%get_mat()
        end do
    end function get_all_rmats

    module function count_state_gt_zero( self ) result( cnt )
        class(oris), intent(in) :: self
        integer :: i, cnt
        cnt = 0
        do i = 1, self%n
            if( self%o(i)%get_state() > 0 ) cnt = cnt + 1
        end do
    end function count_state_gt_zero

    module function ori2str( self, i ) result( str )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(string) :: str
        str = self%o(i)%ori2str()
    end function ori2str

    module subroutine ori2json( self, i, json_ori, boxes)
        class(oris),                  intent(in)    :: self
        integer,                      intent(in)    :: i
        logical,          optional,   intent(in)    :: boxes
        type(json_value), pointer,    intent(inout) :: json_ori
        type(json_core)                             :: json
        logical :: l_boxes = .false.
        if(present(boxes)) l_boxes = boxes
        call self%o(i)%ori2json(json_ori, boxes=l_boxes)
        call json%add(json_ori, 'n', i) 
    end subroutine ori2json

    module subroutine ori2prec( self, i, prec )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: i
        real,        intent(inout) :: prec(N_PTCL_ORIPARAMS)
        call self%o(i)%ori2prec(prec)
    end subroutine ori2prec

    module subroutine prec2ori( self, i, prec )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        real,        intent(in)    :: prec(N_PTCL_ORIPARAMS)
        call self%o(i)%prec2ori(prec)
    end subroutine prec2ori

    module function get_ctfvars( self, i ) result( ctfvars )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        type(ctfparams) :: ctfvars
        ctfvars = self%o(i)%get_ctfvars()
    end function get_ctfvars

    module function has_been_sampled( self )
        class(oris), intent(inout) :: self
        integer :: i
        has_been_sampled =.false.
        do i = 1,self%n
            if( self%o(i)%get_state() > 0 )then
                if( self%o(i)%get_sampled() > 0 )then
                    has_been_sampled = .true.
                    exit
                endif
            endif
        end do
    end function has_been_sampled

    module function has_been_searched( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        has_been_searched = self%o(i)%has_been_searched()
    end function has_been_searched

    module function any_state_zero( self )
        class(oris), intent(in) :: self
        integer :: i
        any_state_zero = .false.
        do i=1,self%n
            if( self%o(i)%get_state() == 0 )then
                any_state_zero = .true.
                return
            endif
        end do
    end function any_state_zero

    module function is_first_update( self, iter, iptcl )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: iter, iptcl
        is_first_update = (self%o(iptcl)%get_int('updatecnt') == 1) .and. (iter > 1)
    end function is_first_update

    module function get_update_frac( self ) result( update_frac )
        class(oris), intent(inout) :: self
        integer :: updatecnts(self%n), sampled(self%n), states(self%n), updatecnt_max, sampled_max, i
        real    :: update_frac
        sampled_max   = 0
        updatecnt_max = 0
        do i = 1, self%n
            updatecnts(i) = self%o(i)%get_updatecnt()
            sampled(i)    = self%o(i)%get_sampled()
            states(i)     = self%o(i)%get_state()
            sampled_max   = max(sampled_max,sampled(i))
            updatecnt_max = max(updatecnt_max,updatecnts(i))
        end do
        if( sampled_max   == 0 ) THROW_HARD('requires previous sampling')
        if( updatecnt_max == 0 ) THROW_HARD('requires previous update')
        update_frac = real(count(sampled == sampled_max .and. states > 0)) / real(count(updatecnts > 0 .and. states > 0))
    end function get_update_frac

    module subroutine get_class_sample_stats( self, clsinds, clssmp, label )
        class(oris),                     intent(inout) :: self
        integer,                         intent(in)    :: clsinds(:)
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        character(len=*),      optional, intent(in)    :: label
        character(len=:), allocatable :: flag
        integer :: n, i, j, nc
        if( present(label) )then
            flag = trim(label)
        else
            flag = 'class'
        endif
        n = size(clsinds)
        if( allocated(clssmp) )then
            nc = size(clssmp)
            do i = 1, nc
                if( allocated(clssmp(i)%ccs)   ) deallocate(clssmp(i)%ccs)
                if( allocated(clssmp(i)%pinds) ) deallocate(clssmp(i)%pinds)
            end do
            deallocate(clssmp)
        endif
        allocate(clssmp(n))
        do i = 1, n
            call self%get_pinds(clsinds(i), flag, clssmp(i)%pinds)
            if( allocated(clssmp(i)%pinds) )then
                clssmp(i)%clsind = clsinds(i)
                clssmp(i)%pop    = size(clssmp(i)%pinds)
                allocate(clssmp(i)%ccs(clssmp(i)%pop), source=0.)
                do j = 1, clssmp(i)%pop
                    clssmp(i)%ccs(j) = self%o(clssmp(i)%pinds(j))%get('corr')
                end do
                call hpsort(clssmp(i)%ccs, clssmp(i)%pinds)
                call reverse(clssmp(i)%ccs)
                call reverse(clssmp(i)%pinds)
            endif
        end do
    end subroutine get_class_sample_stats

    module subroutine get_proj_sample_stats( self, eulspace, clssmp )
        class(oris),                     intent(inout) :: self
        class(oris),                     intent(in)    :: eulspace
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        integer, allocatable :: tmpinds(:), clsinds(:), clspops(:)
        integer              :: ncls, icls
        type(oris)           :: self_copy
        ncls = eulspace%get_noris()
        call self_copy%copy(self)
        call self_copy%set_projs(eulspace)
        call self_copy%proj2class
        allocate(clspops(ncls))
        do icls=1,ncls
            clspops(icls) = self_copy%get_pop(icls, 'class')
        end do
        tmpinds = (/(icls,icls=1,ncls)/)
        clsinds = pack(tmpinds, mask= clspops > 0)
        call self_copy%get_class_sample_stats(clsinds, clssmp)
        call self_copy%kill
    end subroutine get_proj_sample_stats

end submodule simple_oris_getters
