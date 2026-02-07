!@descr: utilities for convolution interpolation (gridding)
module simple_gridding
use simple_core_module_api
use simple_image,     only: image
use simple_projector, only: projector
implicit none

public :: prep2D_inv_instrfun4mul, prep3D_inv_instrfun4mul
private
#include "simple_local_flags.inc"

contains

    ! THIS SEEMS TO BE THE CORRECT WAY OF PREPARING THE IMAGE FOR INSTRUMENT FUNCTION DIVISION
    ! INSTRUMENT FUNCTION DIVISION OUGHT TO NOT MATTER AS LONG AS WE STAY IN THE FOURIER DOMAIN
    ! HOWEVER WHEN THE IMAGES THAT HAVE BEEN INTERPOLATED WITH THE KB APODIZATION FUNCTION AND NEED TO
    ! BE BACK-TRANSFORMED TO CREATE REAL-SPACE OUTPUT, DIVISION WITH THIS IMAGE IS NECESSARY POST IFFT
    ! OF AVERAGED INTERPOLATED OUTPUT
    ! >  \brief  corrects for KB Fourier interpolation
    ! subroutine prep_instrfun4div( img )
    !     class(image), intent(inout) :: img
    !     type(kbinterpol) :: kbwin
    !     real    :: center(3),dist(2),pid,sinc,pad_sc,kbzero,vj,v
    !     integer :: i,j
    !     call img%new(ldim_crop,smpd_crop)
    !     center = real(ldim_crop/2 + 1)
    !     pad_sc = 1. / real(ldim_croppd(1))
    !     kbwin  = kbinterpol(KBWINSZ, KBALPHA2D)
    !     kbzero = kbwin%instr(0.)
    !     do j = 1, ldim_crop(2)
    !         dist(2) = pad_sc * (real(j) - center(2))
    !         vj      = kbwin%instr(dist(2)) / kbzero**2
    !         do i = 1, ldim_crop(1)
    !             dist(1) = pad_sc * (real(i) - center(1))
    !             call img%set([i,j,1], kbwin%instr(dist(1)) * vj)
    !         enddo
    !     enddo
    ! end subroutine prep_instrfun4div

    !=============================================================
    ! 2D: SIMPLE convention => ldim_crop(3)=1, square (n x n)
    ! Returns inverse instrument function for multiplication.
    !=============================================================
    function prep2D_inv_instrfun4mul( ldim_crop, ldim_croppd, smpd_crop ) result( img )
        integer,      intent(in) :: ldim_crop(3), ldim_croppd(3)
        real,         intent(in) :: smpd_crop
        type(image)              :: img
        real(c_float), parameter :: EPS_DIV = 1.0e-8_c_float
        type(kbinterpol)         :: kbwin
        real(c_float)            :: center, dist, pad_sc, kbzero, kv
        integer                  :: i, j, n
        real(c_float)            :: inv_1d(ldim_crop(1))
        ! Determine single dimension (assumes square)
        n = ldim_crop(1)
        if (ldim_crop(2) /= n .or. ldim_crop(3) /= 1) then
            THROW_HARD('prep2D_inv_instrfun4mul: expected ldim_crop = [n,n,1]')
        end if
        ! Allocate target image to expected size
        call img%new(ldim_crop, smpd_crop)
        ! centre coordinate (Fortran 1-based)
        center = real(n/2 + 1, c_float)
        ! pad scaling: mirror original intent (use provided ldim_croppd(1))
        pad_sc = 1.0_c_float / real(ldim_croppd(1), c_float)
        ! create KB window object (assumes KBWINSZ, KBALPHA2D in scope)
        kbwin = kbinterpol(KBWINSZ, KBALPHA2D)
        ! normalization constant at zero (guard small values)
        kbzero = kbwin%instr(0.0_c_float)
        if ( abs(kbzero) < EPS_DIV ) kbzero = 1.0_c_float
        ! Precompute 1-D inverse vector
        do i = 1, n
            dist = pad_sc * ( real(i, c_float) - center )
            kv   = kbwin%instr(dist)
            if ( abs(kv) < EPS_DIV ) then
                inv_1d(i) = 0.0_c_float
            else
                inv_1d(i) = kbzero / kv
            end if
        end do
        ! Fill image with outer product inv_1d(i) * inv_1d(j)
        do j = 1, n
            do i = 1, n
                call img%set([i, j, 1], inv_1d(i) * inv_1d(j))
            end do
        end do
    end function prep2D_inv_instrfun4mul

    !=============================================================
    ! 3D: cube (n x n x n)
    ! Returns inverse instrument function for multiplication.
    !=============================================================
    function prep3D_inv_instrfun4mul( ldim_crop, ldim_croppd, smpd_crop ) result( img )
        integer,      intent(in) :: ldim_crop(3), ldim_croppd(3)
        real,         intent(in) :: smpd_crop
        type(image)              :: img
        real(c_float), parameter :: EPS_DIV = 1.0e-8_c_float
        type(kbinterpol)         :: kbwin
        real(c_float)            :: center, dist, pad_sc, kbzero, kv, inv_1d(ldim_crop(1))
        integer                  :: i, j, k, n
        ! Determine single dimension (assumes cube)
        n = ldim_crop(1)
        if (ldim_crop(2) /= n .or. ldim_crop(3) /= n) then
            THROW_HARD('prep3D_inv_instrfun4mul: ldim_crop must be cube (n x n x n)')
        end if
        ! Allocate target image to expected size
        call img%new(ldim_crop, smpd_crop)
        ! centre coordinate (Fortran 1-based)
        center = real(n/2 + 1, c_float)
        ! pad scaling: mirror original intent (use provided ldim_croppd(1))
        pad_sc = 1.0_c_float / real(ldim_croppd(1), c_float)
        ! create KB window object (assumes KBWINSZ, KBALPHA3D in scope)
        kbwin = kbinterpol(KBWINSZ, KBALPHA3D)
        ! normalization constant at zero (guard small values)
        kbzero = kbwin%instr(0.0_c_float)
        if ( abs(kbzero) < EPS_DIV ) kbzero = 1.0_c_float
        ! Precompute 1-D inverse vector
        do i = 1, n
            dist = pad_sc * ( real(i, c_float) - center )
            kv   = kbwin%instr(dist)
            if ( abs(kv) < EPS_DIV ) then
                inv_1d(i) = 0.0_c_float
            else
                inv_1d(i) = kbzero / kv
            end if
        end do
        ! Fill image with triple product inv_1d(i) * inv_1d(j) * inv_1d(k)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    call img%set([i, j, k], inv_1d(i) * inv_1d(j) * inv_1d(k))
                end do
            end do
        end do
    end function prep3D_inv_instrfun4mul

end module simple_gridding
