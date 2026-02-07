!@descr: utilities for convolution interpolation (gridding)
module simple_gridding
use simple_core_module_api
use simple_image,     only: image
use simple_projector, only: projector
implicit none

public :: mul_w_instr
private
#include "simple_local_flags.inc"

contains

    ! THIS SEEMS TO BE THE CORRECT WAY OF PREPARING THE IMAGE FOR INSTRUMENT FUNCTION DIVISION
    ! INSTRUMENT FUNCTION DIVISION OUGHT TO NOT MATTER AS LONG AS WE STAY IN THE FOURIER DOMAIN
    ! HOWEVER WHEN THE IMAGES THAT HAVE BEEN INTERPOLATED WITH THE KB APODIZATION FUNCTION AND NEED TO
    ! BE BACK-TRANSFORMED TO CREATE REAL-SPACE OUTPUT, DIVISION WITH THIS IMAGE IS NECESSARY POST IFFT
    ! OF AVERAGED INTERPOLATED OUTPUT
    !>  \brief  corrects for Fourier domain bilinear interpolation
    ! subroutine cavger_prep_gridding_correction( img )
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
    ! end subroutine cavger_prep_gridding_correction

    !>  \brief generates instrument function image for division of real-space images
    subroutine gen_instrfun_img( instrfun_img, kbwin )
        class(image),      intent(inout) :: instrfun_img
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: ldim(3), center(3), i,j,k, iarg, dim
        ldim = instrfun_img%get_ldim()
        if( any(ldim==0) .or. instrfun_img%is_ft() .or. .not.instrfun_img%square_dims() )then
            THROW_HARD('Erroneous image in gen_instrfun_img')
        endif
        center = ldim/2+1
        dim    = ldim(1)
        ! Kaiser-bessel window
        allocate(w(ldim(1)),source=1.)
        do i = 1,ldim(1)
            arg  = real(i-center(1))/real(dim)
            w(i) = kbwin%instr(arg)
        end do
        w = w / kbwin%instr(0.)
        if( instrfun_img%is_2d() )then
            !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close) schedule(static)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    call instrfun_img%set([i,j,1], w(i)*w(j))
                enddo
            enddo
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) private(i,j,k) default(shared) proc_bind(close) schedule(static)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    do k = 1,ldim(3)
                        call instrfun_img%set([i,j,k], w(i)*w(j)*w(k))
                    enddo
                enddo
            enddo
            !$omp end parallel do
        endif
    end subroutine gen_instrfun_img

    !> \brief  for dividing a real image with the instrument function
    subroutine mul_w_instr( vol, kbwin )
        class(image),      intent(inout) :: vol
        class(kbinterpol), intent(in)    :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: ldim(3), center(3), i,j,k, iarg
        if( vol%is_ft() ) THROW_HARD('Volume must not be FTed')
        ! get the limits
        ldim = vol%get_ldim()
        center = ldim/2+1
        allocate(w(ldim(1)),source=1.)
        ! kaiser-besel window
        do i = 1,ldim(1)
            arg  = real(i-center(1))/real(ldim(1))
            w(i) = kbwin%instr(arg)
        end do
        w = w / kbwin%instr(0.)
        !$omp parallel do collapse(3) private(i,j,k) default(shared) proc_bind(close) schedule(static)
        do i = 1,ldim(1)
            do j = 1,ldim(2)
                do k = 1,ldim(3)
                    call vol%mul([i,j,k], w(i)*w(j)*w(k))
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine mul_w_instr

end module simple_gridding
