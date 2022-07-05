! utilities for convolution interpolation (gridding)
module simple_gridding
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image,      only: image
use simple_projector,  only: projector
use simple_kbinterpol, only: kbinterpol
implicit none

public :: mul_w_instr, gen_instrfun_img
private
#include "simple_local_flags.inc"

contains

    !>  \brief generates instrument function image for division of real-space images
    subroutine gen_instrfun_img( instrfun_img, interpfun, kbwin, padded_dim )
        class(image),                intent(inout) :: instrfun_img
        character(len=*),            intent(in)    :: interpfun
        class(kbinterpol), optional, intent(in)    :: kbwin
        integer,           optional, intent(in)    :: padded_dim
        real, allocatable :: w(:)
        real    :: g, arg
        integer :: ldim(3), center(3), i,j,k, iarg, dim
        ldim = instrfun_img%get_ldim()
        if( any(ldim==0) .or. instrfun_img%is_ft() .or. .not.instrfun_img%square_dims() )then
            THROW_HARD('Erroneous image in gen_instrfun_img')
        endif
        center = ldim/2+1
        dim    = ldim(1)
        if( present(padded_dim) ) dim = padded_dim
        select case(trim(interpfun))
        case('kb')
            ! Kaiser-bessel window
            if(.not.present(kbwin)) THROW_HARD('KB interpolator must be gigen with interpfun=kb')
            allocate(w(ldim(1)),source=1.)
            do i = 1,ldim(1)
                arg  = real(i-center(1))/real(dim)
                w(i) = kbwin%instr(arg)
            end do
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
        case('linear')
            if( instrfun_img%is_2d() )then
                !$omp parallel do collapse(2) private(i,j,iarg,arg) default(shared) proc_bind(close) schedule(static)
                do i = 1,ldim(1)
                    do j = 1,ldim(2)
                        iarg = sum(([i,j]-center(1:2))**2)
                        if( iarg == 0 )then
                            arg = 1.0
                        else
                            arg = PI*sqrt(real(iarg)) / dim
                            arg = sin(arg)/arg
                            arg = arg*arg ! normalised sinc^2
                        endif
                        call instrfun_img%set([i,j,1], arg)
                    enddo
                enddo
                !$omp end parallel do
            else
                !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
                do i = 1,ldim(1)
                    do j = 1,ldim(2)
                        do k = 1,ldim(3)
                            iarg = sum(([i,j,k]-center)**2)
                            if( iarg == 0 )then
                                arg = 1.0
                            else
                                arg = PI * sqrt(real(iarg)) / dim
                                arg = sin(arg)/arg
                                arg = arg*arg ! normalised sinc^2
                            endif
                            call instrfun_img%set([i,j,k], arg)
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            endif
        case DEFAULT
            THROW_HARD('Unsupported interpolation method: '//trim(interpfun))
        end select
    end subroutine gen_instrfun_img


    !> \brief  for dividing a real image with the instrument function
    subroutine mul_w_instr( vol, interpfun, kbwin )
        class(image),                intent(inout) :: vol
        character(len=*),            intent(in)    :: interpfun
        class(kbinterpol), optional, intent(in)    :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: ldim(3), center(3), i,j,k, iarg
        if( vol%is_ft() ) THROW_HARD('Volume must not be FTed')
        ! get the limits
        ldim = vol%get_ldim()
        center = ldim/2+1
        select case(trim(interpfun))
        case('kb')
            if(.not.present(kbwin)) THROW_HARD('KB interpolator must be gigen with interpfun=kb')
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
        case('linear')
            ! tri-linear interpolation
            !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
            do i = 1,ldim(1)
                do j = 1,ldim(2)
                    do k = 1,ldim(3)
                        iarg = sum(([i,j,k]-center)**2)
                        if( iarg == 0 )cycle
                        arg = PI * (sqrt(real(iarg)) / ldim(1))
                        arg = sin(arg)/arg
                        call vol%mul([i,j,k], arg*arg)
                    enddo
                enddo
            enddo
            !$omp end parallel do
        case DEFAULT
            THROW_HARD('Unsupported interpolation method: '//trim(interpfun))
        end select
    end subroutine mul_w_instr

end module simple_gridding
