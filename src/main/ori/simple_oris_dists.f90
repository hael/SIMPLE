! @descr: distance calculators and related functions for oris
submodule (simple_oris) simple_oris_dists
implicit none
#include "simple_local_flags.inc"

contains

    module pure real function euldist_1( self, i, j )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i, j
        euldist_1 = self%o(i).euldist.self%o(j)
    end function euldist_1

    module pure real function euldist_2( self, i, o )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        class(ori),  intent(in) :: o
        euldist_2 = self%o(i).euldist.o
    end function euldist_2

    module subroutine min_euldist( self, o_in, mindist )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: o_in
        real,        intent(inout) :: mindist
        real      :: dists(self%n), x
        integer   :: inds(self%n), i, loc(1)
        type(ori) :: o
        dists = huge(x)
        do i=1,self%n
            inds(i) = i
            call self%get_ori(i, o)
            dists(i) = o.euldist.o_in
        end do
        loc = minloc(dists)
        mindist = rad2deg(dists(loc(1)))
    end subroutine min_euldist

    module function find_angres( self ) result( res )
        class(oris), intent(in) :: self
        real    :: dists(self%n), dists_max(self%n), x, nearest3(3), res
        integer :: i, j
        !$omp parallel do default(shared) proc_bind(close) private(j,i,dists,nearest3)
        do j=1,self%n
            do i=1,self%n
                if( i == j )then
                    dists(i) = huge(x)
                else
                    dists(i) = self%o(i).euldist.self%o(j)
                endif
            end do
            nearest3     = min3(dists)
            dists_max(j) = maxval(nearest3)
        end do
        !$omp end parallel do
        res = rad2deg(maxval(dists_max))
    end function find_angres

end submodule simple_oris_dists