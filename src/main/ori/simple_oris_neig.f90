!@descr: routines for identifying neighbors of oris and calculating correlations/overlaps/distances between oris
submodule (simple_oris) simple_oris_neigh
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module function find_closest_proj( self, o_in ) result( closest )
        class(oris), intent(in) :: self
        class(ori),  intent(in) :: o_in
        real    :: dists(self%n)
        integer :: closest, i
        do i=1,self%n
            dists(i)=self%o(i).euldist.o_in
        end do
        closest = minloc( dists, dim=1 )
    end function find_closest_proj

    module subroutine nearest_proj_neighbors_1( self, k, nnmat )
        class(oris), intent(in)    :: self
        integer,     intent(in)    :: k
        integer,     intent(inout) :: nnmat(k,self%n)
        real      :: dists(self%n)
        integer   :: inds(self%n), i, j
        if( k >= self%n ) THROW_HARD('need to identify fewer nearest_proj_neighbors')
        !$omp parallel do default(shared) proc_bind(close) private(i,j,inds,dists)
        do i=1,self%n
            do j=1,self%n
                inds(j)  = j
                dists(j) = self%o(j).euldist.self%o(i)
            end do
            call hpsort(dists, inds)
            do j=1,k
                nnmat(j,i) = inds(j)
            end do
        end do
        !$omp end parallel do
    end subroutine nearest_proj_neighbors_1

    module subroutine nearest_proj_neighbors_2( self, o, euldist_thres, lnns )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o
        real,        intent(in)    :: euldist_thres
        logical,     intent(inout) :: lnns(self%n)
        real    :: dists(self%n), euldist_thres_rad
        integer :: j
        euldist_thres_rad = deg2rad(euldist_thres)
        do j=1,self%n
            dists(j) = self%o(j).euldist.o
        end do
        where( dists <= euldist_thres_rad ) lnns = .true.
    end subroutine nearest_proj_neighbors_2

    module subroutine nearest_proj_neighbors_3( self, o, k, lnns )
        class(oris), intent(in)    :: self
        class(ori),  intent(in)    :: o
        integer,     intent(in)    :: k
        logical,     intent(inout) :: lnns(self%n)
        real    :: dists(self%n)
        integer :: inds(self%n), i, j
        if( k >= self%n ) THROW_HARD('need to identify fewer nearest_proj_neighbors')
        do i=1,self%n
            do j=1,self%n
                inds(j)  = j
                dists(j) = self%o(j).euldist.o
            end do
            call hpsort(dists, inds)
            do j=1,k
                lnns(inds(j)) = .true.
            end do
        end do
    end subroutine nearest_proj_neighbors_3

    module subroutine replace_with_closest( self, other )
        class(oris), intent(inout) :: self
        class(oris), intent(in)    :: other
        integer :: i, ind_closest
        if( other%n <= self%n ) THROW_HARD('other must have more oris than self; replace_with_closest')
        do i=1,self%n
            ind_closest = other%find_closest_proj(self%o(i))
            self%o(i) = other%o(ind_closest)
        end do
    end subroutine replace_with_closest

    module function corr_oris( self1, self2 ) result( corr )
        class(oris), intent(inout) :: self1, self2
        real :: arr1(5), arr2(5), corr
        integer :: i
        corr = 0.
        do i=1,self1%n
            arr1(1:3) = self1%get_euler(i)
            arr1(4)   = self1%get(i,'x')
            arr1(5)   = self1%get(i,'y')
            arr2(1:3) = self2%get_euler(i)
            arr2(4)   = self2%get(i,'x')
            arr2(5)   = self2%get(i,'y')
            corr = corr+pearsn(arr1,arr2)
        end do
        corr = corr/real(self1%n)
    end function corr_oris

    module subroutine diststat_1( self, sumd, avgd, sdevd, mind, maxd )
        class(oris), intent(in)  :: self
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
        integer :: i, j, cnt
        real    :: dists((self%n*(self%n-1))/2), vard
        logical :: err
        cnt  = 0
        do i=1,self%n-1
           do j=i+1,self%n
              cnt = cnt+1
              dists(cnt) = self%o(i).euldist.self%o(j)
           end do
        end do
        mind = minval(dists)
        maxd = maxval(dists)
        sumd = sum(dists)
        call moment(dists, avgd, sdevd, vard, err )
    end subroutine diststat_1

    module subroutine diststat_2( self1, self2, sumd, avgd, sdevd, mind, maxd )
        use simple_linalg, only: vector_angle_norm
        class(oris), intent(inout)  :: self1, self2
        real,        intent(out) :: mind, maxd, avgd, sdevd, sumd
        real, allocatable :: onormals1(:,:),onormals2(:,:)
        real    :: dists(self1%n), vard, x
        integer :: i
        logical :: err
        if( self1%n /= self2%n )then
            THROW_HARD('cannot calculate distance between sets of different size; euldist_2')
        endif
        mind = huge(x)
        maxd = -mind
        sumd = 0.
        onormals1 = self1%get_all_normals()
        onormals2 = self2%get_all_normals()
        !$omp parallel do default(shared) private(i) proc_bind(close) schedule(static)&
        !$omp reduction(+:sumd) reduction(min:mind) reduction(max:maxd)
         do i=1,self1%n
            dists(i) = vector_angle_norm(onormals1(i,:),onormals2(i,:))
            if( dists(i) < mind ) mind = dists(i)
            if( dists(i) > maxd ) maxd = dists(i)
            sumd = sumd+dists(i)
        end do
        !$omp end parallel do
        call moment(dists, avgd, sdevd, vard, err )
        deallocate(onormals1,onormals2)
    end subroutine diststat_2

    module real function overlap( self1, self2, which, state )
        class(oris),      intent(inout) :: self1, self2
        character(len=*), intent(in)    :: which
        integer,          intent(in)    :: state
        real,    allocatable :: arr1(:), arr2(:), tmp(:), ows(:), states(:)
        logical, allocatable :: mask1(:), mask2(:)
        integer :: n1, n2, sz, n_min, i
        overlap = 0.
        if( self1%n == 0 .or. self2%n == 0 ) return
        if( .not. self1%isthere(trim(which)) )then
            THROW_WARN('key: '//trim(which)//' not present in self1; overlap')
            return
        endif
        if( .not. self2%isthere(trim(which)) )then
            THROW_WARN('key: '//trim(which)//' not present in self2; overlap')
            return
        endif
        if( .not. self1%isthere('state') )     THROW_HARD('key: state not present in self1; overlap')
        if( .not. self2%isthere('state') )     THROW_HARD('key: state not present in self2; overlap')
        if( .not. self1%isthere('ow') )        THROW_HARD('key: ow not present in self1; overlap')
        if( .not. self2%isthere('ow') )        THROW_HARD('key: ow not present in self2; overlap')
        ows    = self1%get_all('ow')
        states = self1%get_all('state')
        where(nint(states) .ne. state) ows = 0.
        if( .not. any(ows > TINY) ) return
        tmp    = self1%get_all(trim(which))
        arr1   = pack(tmp, mask=ows > TINY)
        n1     = size(arr1)
        ows    = self2%get_all('ow')
        states = self2%get_all('state')
        where(nint(states) .ne. state) ows = 0.
        if( .not. any(ows > TINY) ) return
        tmp    = self2%get_all(trim(which))
        arr2   = pack(tmp, mask=ows > TINY)
        n2     = size(arr2)
        sz = nint(max(maxval(arr1), maxval(arr2)))
        allocate(mask1(sz), mask2(sz), source=.false.)
        forall(i=1:n1) mask1(nint(arr1(i))) = .true.
        forall(i=1:n2) mask2(nint(arr2(i))) = .true.
        n_min   = min(count(mask1),count(mask2))
        overlap = real(count(mask1 .and. mask2)) / real(n_min)
    end function overlap

end submodule simple_oris_neigh
