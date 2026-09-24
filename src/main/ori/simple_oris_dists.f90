!@descr: distance calculators and related functions for oris
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

    !> the largest angle (degrees) from a direction of the set to its third-nearest neighbour;
    !! simple_exec prg=measure_projspace_angres measures it for the 3D search space (with its
    !! point group). Values for the full-sphere spiral (oris%spiral) of n directions, from the
    !! former angres test (500 to 20000 in steps of 500; 3D searches use at most 20000):
    !!     n  angres     n  angres     n  angres     n  angres
    !!   500  9.9531  5500  3.0090 10500  2.1878 15500  1.8039
    !!  1000  7.0143  6000  2.8865 11000  2.1318 16000  1.7749
    !!  1500  5.7937  6500  2.7653 11500  2.0926 16500  1.7466
    !!  2000  4.9782  7000  2.6863 12000  2.0507 17000  1.7185
    !!  2500  4.4927  7500  2.5882 12500  2.0077 17500  1.6967
    !!  3000  4.0799  8000  2.5080 13000  1.9652 18000  1.6737
    !!  3500  3.7839  8500  2.4351 13500  1.9348 18500  1.6469
    !!  4000  3.5436  9000  2.3635 14000  1.8960 19000  1.6281
    !!  4500  3.3389  9500  2.3016 14500  1.8642 19500  1.6091
    !!  5000  3.1713 10000  2.2460 15000  1.8334 20000  1.5846
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