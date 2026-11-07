!@descr: distance calculators and related functions for oris
submodule (simple_oris) simple_oris_dists
implicit none
#include "simple_local_flags.inc"

contains

    module pure real function euldist_inds( self, i, j )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i, j
        euldist_inds = self%o(i).euldist.self%o(j)
    end function euldist_inds

    module pure real function euldist_ori( self, i, o )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        class(ori),  intent(in) :: o
        euldist_ori = self%o(i).euldist.o
    end function euldist_ori

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

end submodule simple_oris_dists