! utility for interpolating and integrating a real space grid
module simple_intg_atompeak
use simple_defs        ! use all in there
use simple_image,      only: image
implicit none

public :: set_intgvol, intg_nn, intg_shell

private
#include "simple_local_flags.inc"

real,     allocatable :: rmat(:,:,:)
real                  :: smpd        = 0.
integer               :: rmat_dim(3) = 0
logical               :: vol_set     = .false.

contains

    !>  \brief  
    subroutine set_intgvol( vol )
        class(image), intent(inout) :: vol
        call reset_vol
        if(.not.vol%is_3d())stop 'volume only! simple_intg_atom_peak%set_vol'
        if(vol%is_ft())stop 'Real space volume only! simple_intg_atom_peak%set_vol'
        rmat_dim = vol%get_ldim()
        smpd     = vol%get_smpd()
        rmat     = vol%get_rmat()
        vol_set  = .true.
    end subroutine set_intgvol

    !>  \brief  nearest neighbour interpolation
    real function intg_nn( xyz_in )
        real, intent(in) :: xyz_in(3)
        real    :: xyz(3), dists_sq(8), vec(3), vals(8)
        integer :: cnt, i,j,k, loc(1), fxyz(3)
        if(.not.vol_set)stop 'no volume to interpolate; simple_intg_atom_peak%set_vol'
        xyz  = xyz_in / smpd
        fxyz = floor(xyz) + 1
        cnt  = 0
        do i = fxyz(1), fxyz(1)+1
            do j = fxyz(2), fxyz(2)+1
                do k = fxyz(3), fxyz(3)+1
                    cnt = cnt + 1
                    vec = [real([i,j,k]) - xyz]
                    dists_sq(cnt) = dot_product(vec,vec)
                    vals(cnt) = rmat(i,j,k)
                enddo
            enddo
        enddo
        loc     = minloc(dists_sq)
        intg_nn = vals(loc(1))
    end function intg_nn

    !>  \brief  sum within shell
    real function intg_shell( xyz_in, rad_in )
        real, intent(in) :: xyz_in(3)
        real, intent(in) :: rad_in
        real    :: xyz(3), radius, d_sq, vec(3)
        integer :: cnt, i,j,k, fxyz(3), left(3), right(3), cradius
        xyz     = xyz_in / smpd
        fxyz    = floor(xyz) + 1
        radius  = rad_in / smpd
        cradius = ceiling(radius)
        left    = fxyz - ceiling(radius)
        where(left < 1) left = 1
        right = fxyz + cradius + 1
        right(1) = min(right(1), rmat_dim(1))
        right(2) = min(right(2), rmat_dim(2))
        right(3) = min(right(3), rmat_dim(3))
        cnt = 0
        intg_shell = 0.
        do i = left(1), right(1)
            do j = left(2), right(2)
                do k = left(3), right(3)
                    vec  = real([i,j,k]) - xyz
                    d_sq = dot_product(vec, vec)
                    if( d_sq > radius )cycle
                    cnt = cnt + 1
                    intg_shell = intg_shell + rmat(i,j,k)
                enddo
            enddo
        enddo
    end function intg_shell

    ! DESTRUCTORS

    !>  \brief  
    subroutine reset_vol
        if(allocated(rmat))deallocate(rmat)
        rmat_dim = 0
        smpd = 0.
        vol_set  = .false.
    end subroutine reset_vol

end module simple_intg_atompeak
