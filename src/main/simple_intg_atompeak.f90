! utility for interpolating and integrating a real space grid
module simple_intg_atompeak
include 'simple_lib.f08'
use simple_atoms
use simple_image,      only: image
implicit none

public :: set_intgvol, intg_nn, intg_shell, find_peaks
public :: test_intg_atompeak
private
#include "simple_local_flags.inc"

real, allocatable :: rmat(:,:,:)
real              :: smpd        = 0.
real              :: msk         = 0.
integer           :: rmat_dim(3) = 0
logical           :: vol_set     = .false.

contains

    subroutine set_intgvol( vol, rmsk )
        class(image),   intent(inout) :: vol
        real, optional, intent(in)    :: rmsk
        type(image) :: tmpvol
        call reset_vol
        if(.not.vol%is_3d()) THROW_HARD('volume only! set_intgvol')
        if(vol%is_ft())      THROW_HARD('Real space volume only! set_intgvol')
        rmat_dim = vol%get_ldim()
        smpd     = vol%get_smpd()
        if( present(rmsk) )then
            tmpvol = vol
            call tmpvol%mask(rmsk, 'soft')
            rmat = tmpvol%get_rmat()
        else
            rmat = vol%get_rmat()
        endif
        vol_set  = .true.
    end subroutine set_intgvol

    !>  \brief  nearest neighbour interpolation
    real function intg_nn( xyz_in )
        real, intent(in) :: xyz_in(3)
        real    :: xyz(3), dists_sq(8), vec(3), vals(8)
        integer :: cnt, i,j,k, loc(1), fxyz(3)
        if(.not.vol_set) THROW_HARD('no volume to interpolate; intg_nn')
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
        real    :: xyz(3), radius, radius_sq, d_sq, vec(3)
        integer :: cnt, i,j,k, fxyz(3), left(3), right(3), cradius
        xyz       = xyz_in / smpd
        fxyz      = floor(xyz) + 1
        radius    = rad_in / smpd
        radius_sq = radius*radius
        cradius   = ceiling(radius)
        left      = fxyz - ceiling(radius)
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
                    if( d_sq > radius_sq )cycle
                    cnt = cnt + 1
                    intg_shell = intg_shell + rmat(i,j,k)
                enddo
            enddo
        enddo
    end function intg_shell

    subroutine find_peaks( npeaks, radius, fname )
        integer,          intent(in) :: npeaks
        real,             intent(in) :: radius
        character(len=*), intent(in) :: fname
        type(atoms)           :: mol
        logical, allocatable  :: mask(:,:,:)
        integer, allocatable  :: coords(:,:)
        character(len=STDLEN) :: fbody
        real                  :: grid_radius, grid_radius_sq, d_sq, vec(3)
        integer               :: i,j,k,loc(3), left(3), right(3), ipeak, peak_cnt
        grid_radius    = 2.*radius / smpd
        grid_radius_sq = grid_radius**2.
        allocate(mask(rmat_dim(1),rmat_dim(2),rmat_dim(3)), source=.true.)
        allocate(coords(npeaks,3), source=-1)
        where(rmat==0.)mask = .false.
        do ipeak = 1, npeaks
            call progress(ipeak, npeaks)
            ! current highest
            loc = maxloc(rmat, mask=mask)
            if( rmat(loc(1),loc(2),loc(3)) < TINY )exit
            ! updates coodinates
            coords(ipeak,:) = loc
            peak_cnt = ipeak
            ! updates mask
            left  = floor(real(loc) - grid_radius)
            where( left < 1 )left = 1
            right = ceiling(real(loc) + grid_radius)
            right(1) = min(right(1), rmat_dim(1))
            right(2) = min(right(2), rmat_dim(2))
            right(3) = min(right(3), rmat_dim(3))
            do i = left(1), right(1)
                do j = left(2), right(2)
                    do k = left(3), right(3)
                        vec = real([i,j,k] - loc)
                        d_sq = dot_product(vec, vec)
                        if( d_sq > grid_radius_sq )cycle
                        mask(i,j,k) = .false.
                    enddo
                enddo
            enddo
            if(count(mask) < nint(grid_radius**3.)) exit
        enddo
        ! output
        fbody = trim(get_fbody(trim(fname), trim('pdb')))
        call mol%new(peak_cnt)
        do i = 1, peak_cnt
            call mol%set_chain(i, 'A')
            call mol%set_num(i, i)
            call mol%set_resnum(i, i)
            call mol%set_coord(i, real(coords(i,:)-1)*smpd)
        enddo
        call mol%writepdb(fbody)
    end subroutine find_peaks

    subroutine test_intg_atompeak
        type(image) :: vol
        type(atoms) :: mol
        integer     :: i, natoms
        natoms = 128
        call vol%new([128,128,128], 2.)
        call vol%square(32)
        call mol%new(natoms)
        do i = 1, natoms
            call mol%set_coord(i,2.*real([i-1,64,64]))
        enddo
        call set_intgvol(vol)
        call vol%write('cube.mrc')
        do i = 1, natoms
            if(intg_nn(mol%get_coord(i)).ne.rmat(i,64,64))write(logfhandle,*) 'interpolation error for atom:', i
        enddo
        call reset_vol
        write(logfhandle,'(A)')'>>> FINISHED INTGPEAK TEST'
    end subroutine test_intg_atompeak

    ! DESTRUCTORS

    !>  \brief
    subroutine reset_vol
        if(allocated(rmat))deallocate(rmat)
        rmat_dim = 0
        smpd     = 0.
        vol_set  = .false.
    end subroutine reset_vol

end module simple_intg_atompeak
